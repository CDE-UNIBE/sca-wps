import os
import sys
import time
import numpy as np
from ConfigParser import ConfigParser
import rpy2.rinterface as rinterface
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as npri
import osgeo.osr as osr
try:
    import simplejson as json
except ImportError:
    import json

# Prevent R outputs to stdout
f = open(os.devnull, 'w')
sys.stdout = f

# Zoo service return values
SERVICE_FAILED = 4
SERVICE_SUCCEEDED = 3

def BufferStatistics(conf, inputs, outputs):

    mimeType = outputs['bufferstatistics']['mimeType']
    
    lon = float(inputs['lon']['value'])
    lat = float(inputs['lat']['value'])
    
    try:
        epsg = int(inputs['epsg']['value'])
    except ValueError:
        conf["lenv"]["message"] = "Parameter \"epsg\" is not valid or not supported."
        return SERVICE_FAILED
    if epsg != 4326:
        conf["lenv"]["message"] = "Request CRS is not supported."
        return SERVICE_FAILED
    
    coords = reprojectCoordinates((lon,lat), epsg)
    
    config = ConfigParser()
    config.read('Layers.ini')
    alllayers = config.sections()
    
    inputlayer = str(inputs['layer']['value'])     
    if inputlayer == 'all' or inputlayer == 'NULL':
        layers = alllayers
    elif inputlayer in alllayers:
        layers = [inputlayer]
    else:
        conf["lenv"]["message"] = "Request layer (" + inputlayer + ") is not supported"
        return SERVICE_FAILED
    
    try:
        buffer=int(inputs['buffer']['value'])
        if buffer not in [100,500,5000,50000,500000]:
            conf["lenv"]["message"] = "Request buffer (" + str(buffer) + ") is not supported"
            return SERVICE_FAILED
    except:
        if str(inputs['buffer']['value']) == 'NULL':
            buffer = 5000
        else:
            conf["lenv"]["message"] = "Request buffer (" + str(inputs['buffer']['value']) + ") is not supported"
            return SERVICE_FAILED
        
    layerstatistics = []
    for layer in layers:    
        values = getRasterValues(coords[0],coords[1],config.get(layer, 'filepath'), buffer)
        PixelCount = len(values)
        NoDataValues = int(np.isnan(values).sum())
        DataValues = PixelCount - NoDataValues
        valuesnanremove = values[~np.isnan(values)]
        if config.get(layer, 'type') == "Values":
            layerstatistics.append({'layername': config.get(layer, 'name'),
                                    'description': config.get(layer, 'description').replace('\n', ' '),
                                    'datatype': config.get(layer, 'type'), 
                                    'unit': config.get(layer, 'unit'),
                                    'pixel': PixelCount,
                                    'datapixel': DataValues,
                                    'nodatapixel': NoDataValues,
                                    'bufferradius': buffer,
                                    'statistics': calculateStatistics(valuesnanremove)})
        elif config.get(layer, 'type') == "Classes":
            classstring = config.get(layer, 'classes')
            classes = {}
            for cl in classstring.split('\n'):
                c = cl.split(',',1)
                classes[c[0]] = c[1]
            layerstatistics.append({'layername': config.get(layer, 'name'), 
                                    'description': config.get(layer, 'description').replace('\n', ' '),
                                    'datatype': config.get(layer, 'type'), 
                                    'unit': config.get(layer, 'unit'),
                                    'pixel': PixelCount,
                                    'datapixel': DataValues,
                                    'nodatapixel': NoDataValues,
                                    'bufferradius': buffer,
                                    'classes': calculateAreaShare(values,classes)})

                                    
    outputs['bufferstatistics']['value'] = json.dumps({'layers':layerstatistics})
    
#    return SERVICE_FAILED
    return SERVICE_SUCCEEDED

    
def getRasterValues(lon, lat, layer, buffervalue):
    startTime = time.time()
    rinterface.initr()
    r = robjects.r
    r.require('raster')
    ras = r.raster(layer)
    rasvalues = r.extract(ras, r.cbind(lon,lat), buffer=buffervalue, small=True)
    values = npri.ri2numpy(rasvalues[0])
    endTime = time.time()
    print(str(endTime - startTime))
    return values
    
def calculateStatistics(values):
    statistics = [{'name': "Minimum", 'value': values.min()}, 
                  {'name': "Maximum", 'value': values.max()}, 
                  {'name': "Mean", 'value': values.mean()}, 
                  {'name': "Sum", 'value': values.sum()},
                  {'name': "Standard deviation", 'value': values.std()},
                  {'name': "Variance", 'value': values.var()}]
    return statistics

def calculateAreaShare(values, classes):
    PixelCount = len(values)
    NoDataValues = int(np.isnan(values).sum())
    values = values[~np.isnan(values)]
    if NoDataValues > 0:
        areashares = [{'value' : 'No data', 
                       'name' : 'No data', 
                       'frequency': NoDataValues, 
                       'areashare': float(NoDataValues)/PixelCount*100.0}]
    else:
        areashares = []
    for kvp in countUnique(values):
        value = kvp['value']
        frequency = kvp['frequency']
        try:
            classname = classes[str(value)]
        except:
            classname = "Classname unknown"
        areashares.append({'value' : value, 
                           'name' : classname, 
                           'frequency': frequency, 
                           'areashare': float(frequency)/PixelCount*100.0})
    return areashares

def countUnique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    values = uniq_keys.tolist() 
    frequency = np.bincount(bins).tolist()
    result = []
    for i in range(len(values)):
        result.append({'value':int(values[i]),'frequency':frequency[i]})
    return result
    
def reprojectCoordinates(coords, epsg_code):
    """
    Reproject the requested coordinates to Mollweide projection. EPSG
    code of input CRS must be known to GDAL.
    """
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg_code))
    mollweideSrs = osr.SpatialReference()
    # From spatialreference.org: http://spatialreference.org/ref/sr-org/7/
    mollweideSrs.ImportFromProj4("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
    ct = osr.CoordinateTransformation(srs, mollweideSrs)
    (x, y, z) = ct.TransformPoint(float(coords[0]), float(coords[1]))
    return (x, y)