#!/usr/bin/env python

'''
pull MODIS MOD09/MYD09 data from Google Earth Engine
including filtering and calculation of angular kernels

'''

try:
  import ee
except:
  print 'you need to install the google earth engine packages'

import wget
import zipfile
import gdal
from linearBRDFBase import linearBRDFBase

import os
import numpy as np
import sys

unload = lambda self,name,kwargs: ((name in kwargs and kwargs[name]) or \
                              (hasattr(self,name)) and getattr(self,name))

class getModisEE(linearBRDFBase):
  '''
  Class to pull MODIS data from Google Earth Engine

  '''         
  def __init__(self,*args,**kwargs):
    '''
    Initialise class
    '''
    linearBRDFBase.__init__(self,*args,**kwargs)
    try:
      ee.Initialize()
    except:
      self.err('ee.Initialize() Failed: check earth engine login',fatal=True)

    self.parser(*args,**kwargs)
    


  def parser(self,*args,**kwargs):
    '''
    Call to interpret args and keyword args
    '''
    try:
      self.args += args
    except:
      self.args = args
    try:
      self.kwargs.update(kwargs)
    except:
      self.kwargs = kwargs

    # set default area in degrees
    self.centre = unload(self,'centre',kwargs) or [-0.030663,51.547376]
    self.extent = unload(self,'extent',kwargs) or [-0.1,0.1]

    self.recover = unload(self,'recover',kwargs) or False
    self.scale = unload(self,'scale',kwargs) or 500
    self.dumpFreq = unload(self,'dumpFreq',kwargs) or 50
    self.maxn  = unload(self,'maxn',kwargs) or 10000   # max number of datasets to pull
    try:
      self.localParser(*args,**kwargs)
    except:
      pass


  def localParser(self,*args,**kwargs):
    '''
    Application specific parsing
    '''
    self.sortField = unload(self,'sortField',kwargs) or 'system:time_start'
    self.oname = unload(self,'oname',kwargs) or 'MODIS'
    self.maps = unload(self,'maps',kwargs) or [self.maskEmptyPixels,self.maskClouds,\
                                          self.makeVariables,self.addTime,\
                                          self.subtractZero]
    self.modis = unload(self,'modis',kwargs) or ['MOD09GA','MYD09GA']
    self.dates = unload(self,'dates',kwargs) or ['2000-01-01', '2020-03-01']


  def setAOI(self,centre,extent):
    '''
    set area of interest in lat/lon degrees
    from definition of centre and extent (degrees)

    or just set self.topLeft and self.bottomRight manually

    '''
    self.centre = np.array(centre or self.centre)
    self.extent = np.array(extent or self.extent)
    
    self.topLeft     = self.centre - self.extent * 0.5
    self.bottomRight = self.centre + self.extent * 0.5 

  def getModisCollections(self,*args,**kwargs):
    '''
    Specify datasets and time period
    '''
    self.parser(*args,**kwargs)
    collection = ee.ImageCollection(self.modis[0]).filterDate(self.dates[0],self.dates[1])
    for i in self.modis[1:]:
      collectionb   = ee.ImageCollection(i).filterDate(self.dates[0],self.dates[1])
      collection    = ee.ImageCollection(collection.merge(collectionb))
    self.collection = collection

    for m in self.maps:
      self.collection = self.collection.map(m)

    try:
      self.collection = self.collection.sort(self.sortField)
    except:
      self.err('Failed to apply sort field %s'%self.sortField)

    return self.collection
 
  def pullData(self,image,count=0,clean=True,centre=None,extent=None):
  
    count = count or (hasattr(self,'counter') and self.counter) or 0

    # force calculation of topLeft and bottomRight
    # 
    if (centre) or (extent) or (not hasattr(self, 'topLeft')) or (not hasattr(self, 'bottomRight')):
      self.setAOI(centre,extent)

    # define AOI rectangle
    # and associated boundingBox and region
    # then clip to region

    definedAoi  = ee.Geometry.Rectangle(self.topLeft[0],self.topLeft[1],\
                                    self.bottomRight[0],self.bottomRight[1])

    boundingBox = definedAoi.bounds(1)
    region = ee.Geometry(boundingBox.getInfo())\
                        .toGeoJSONString()

    image = image.clip(definedAoi)

    url = image.getDownloadURL({\
      'name':'%s_%06d'%(self.oname,count),\
      'crs': 'EPSG:4326',\
      'scale': self.scale,\
      'region':region\
    })

    if self.verbose: 
      print(url)
      print self.centre,self.extent
      print region

    # op
    filename = wget.download(url, bar=wget.bar_thermometer)

    zf = zipfile.ZipFile(filename, 'r')
    data = {}
    # ensure odir exists
    self.mkdir(self.odir)

    # look for bad data
    badData = False
    for name in zf.namelist():
      if self.verbose: print name
      if name.split('.')[-1] == 'tif':
        # its a tif file
        f = open(self.odir + os.sep + name,'w+b')
        f.write(zf.read(name))
        f.close()
        data[name] = gdal.Open(self.odir + os.sep + name).ReadAsArray()
        if 'sur_refl' in name and data[name].sum() == 0:
          badData = True
          print 'no data'
          break
      else:
        f = open(self.odir+os.sep+name,'w+')
        f.write(zf.read(name))
        f.close()
        data[name] = np.loadtxt(self.odir+os.sep+name)

    # clean up
    if clean: os.remove(filename)
    if badData: return None
    if not hasattr(self,'data'):
      self.data = {}
      for k in data.keys():
        if data[k].ndim == 1:
          self.data[k] = np.atleast_2d(data[k])
        else:
          self.data[k] = np.atleast_3d(data[k].T).T

    else:
      # append data from each key
      # using numpy hstack
      for k in data.keys():
        # make sure all lists
        try:
          if data[k].ndim == 2:
            self.data[k] = np.vstack([self.data[k],np.atleast_3d(data[k].T).T])
          elif data[k].ndim == 1:
            self.data[k] = np.vstack([self.data[k],np.atleast_2d(data[k])])
        except:
          self.data[k] = data[k] 
    return data

  def get(self,*args,**kwargs):
    '''
    download and load into self.data
    '''
    self.getModisCollections(*args,**kwargs)
    dumper = 'dump_%s.tmp'%str(os.getpid())
    if self.verbose: print 'dumper file',dumper
    # how many items?
    print ee.ImageCollection.size

    try:
      for i in xrange(maxn):
        if self.verbose: print i
        self.pullData(ee.ImageCollection([self.collection.toList(1,i).get(-1)]).min(),count=0)
        if i%self.dumpFreq == 0:
          if self.verbose:
            print 'temp dump at count %d to'%i,dumper
          self.save(dumper)
    except:
      pass

  def load(self,*args,**kwargs):
    '''
    load dataset as pickle for now
    '''
    import pickle
    ofile = (len(args) and type(args[0] == 'str') and args[0]) or \
            (('pfile' in kwargs) and kwargs['pfile']) or \
            (self.odir + os.sep + self.oname + '.pkl')
    output = open(ofile, 'rb') 
    self.data = pickle.load(output)
    output.close()

  def save(self,*args,**kwargs):
    '''
    Save dataset as pickle for now
    '''
    import pickle
    if not hasattr(self,'data'):
      return

    ofile = (len(args) and type(args[0] == 'str') and args[0]) or \
            (('pfile' in kwargs) and kwargs['pfile']) or \
            (self.odir + os.sep + self.oname + '.pkl')
    output = open(ofile, 'wb')
    pickle.dump(self.data, output)
    output.close()

  def getSingle(self,*args,**kwargs):
    '''
    download and load into self.data
    '''
    self.counter = unload('counter',kwargs) or (hasattr(self,'counter')) \
                     or (len(args) and int(args[0])) or 0
    if not hasattr(self,'collection'):
      self.getModisCollections()
    self.pullData(ee.ImageCollection([self.collection.toList(1,count).get(-1)]).min(),count=0)
    self.counter += 1
    return self.data

  def addTime(self,image):
    return image.addBands(image.metadata('system:time_start').float().divide(1000 * 60 * 60 * 24));

  def addSensor(self,image):
    return image.addBands();

  def subtractZero(self,image):
  
    # in degrees
    vza_0 = 0
    sza_0 = 0
    vaa_0 = 0
    saa_0 = 0
  
    # in degrees * 100
    zero_angles = ee.Image([ee.Image(sza_0/0.01).rename(['SolarZenith']),\
                                        ee.Image(saa_0/0.01).rename(['SolarAzimuth']),\
                                        ee.Image(vza_0/0.01).rename(['SensorZenith']),\
                                        ee.Image(vaa_0/0.01).rename(['SensorAzimuth']),\
                                        ])
    # get the kernels
    kernels_0 = self.makeVariables(zero_angles)

    image.addBands(image.select('ross').subtract(kernels_0.select('ross')),['ross'],True)
    image.addBands(image.select('li').subtract(kernels_0.select('li')),['li'],True)

    return image

  def makeVariables(self,image):
    # get metedata list
    properties = image.propertyNames()
    print properties

    # linear kernel models code:
    # after: https://github.com/profLewis/modisPriors/blob/master/python/kernels.py
    BR = 1.0;
    HB = 2.0;
    d2r = np.pi / 180.0;
    zthresh = 0.00001
  
    # interpret view and illumination angles
    sza = image.select('SolarZenith').float().multiply(ee.Number(0.01*d2r));
    vza = image.select('SensorZenith').float().multiply(ee.Number(0.01*d2r));
    vaa = image.select('SensorAzimuth').float().multiply(ee.Number(0.01*d2r));
    saa = image.select('SolarAzimuth').float().multiply(ee.Number(0.01*d2r));
    raa = vaa.subtract(saa)
    raa_plus = raa.add(ee.Number(np.pi))
 
    # correct bounds
    w = vza.lt(0);
    vza = ee.Image(vza.where(w,vza.multiply(ee.Number(-1)))).rename(['vza']);
    raa = ee.Image(raa.where(w,raa_plus));
    w = sza.lt(0);
    sza = ee.Image(sza.where(w,sza.multiply(ee.Number(-1)))).rename(['sza']);
    raa = ee.Image(raa.where(w,raa_plus)); 
    raa = ee.Image(0).expression('raa % (2*pi)',{'raa': raa,'pi':np.pi}).rename(['raa']);     

    # trig functions
    cos_vza = vza.cos().rename(['cos_vza'])
    sin_vza = vza.sin().rename(['sin_vza'])
    cos_sza = sza.cos().rename(['cos_sza'])
    sin_sza = sza.sin().rename(['sin_sza'])
    cos_raa = raa.cos().rename(['cos_raa'])
    sin_raa = raa.sin().rename(['sin_raa'])
    tanti   = sza.tan().rename(['tanti'])
    tantv   = vza.tan().rename(['tantv'])
  
    # trig GO corrected angles: illumination
    tan1    = ee.Image(ee.Image(0).expression('BR*tan1',{'tan1': tanti,'BR':BR}));
    angp1 = tan1.atan();
    sin1 = angp1.sin();
    cos1 = angp1.cos();
    w = cos1.lte(zthresh);
    cos1 = cos1.where(w,zthresh); 
  
    # trig GO corrected angles: view
    tan2 = ee.Image(ee.Image(0).expression('BR*tan1',{'tan1': tantv,'BR':BR}));
    angp2 = tan2.atan();
    sin2 = angp2.sin();
    cos2 = angp2.cos();
  
    # avoid cos == 0 by setting threshold zthresh
    w = cos2.lte(zthresh);
    cos2 = cos2.where(w,zthresh);   
  
  
    # phase angle
    cdict = {'cos1': cos_vza,'sin1': sin_vza,'cos2': cos_sza,'sin2': sin_sza,'cos3':cos_raa};
    cosphaang = ee.Image(0).expression('cos1*cos2 + sin1*sin2*cos3',cdict);   
    # make sure limited -1 to 1
    w = cosphaang.lte(-1);
    cosphaang = ee.Image(cosphaang.where(w,-1));
    w = cosphaang.gte(1);
    cosphaang = ee.Image(cosphaang.where(w,1)).rename(['cos_phaang']);
    phaang = cosphaang.acos().rename(['phaang']);
    sinphaang = phaang.sin().rename(['sin_phaang']);
  
  
    # ross kernel
    cdict = {'cosphaang': cosphaang,'sinphaang': sinphaang,'pi': np.pi, 'phaang':phaang,
        'cos1': cos_vza, 'cos2': cos_sza};
    ross = ee.Image(0).expression('((pi/2. - phaang)*cosphaang+sinphaang)/(cos1+cos2)',cdict).rename(['ross']);
  
    # Li kernel
    cdict = {'tan1': tan1,'tan2': tan2,'cos3':cos_raa};
    temp = ee.Image(0).expression('tan1*tan1 + tan2*tan2*cos3',cdict);
    w = temp.lte(0);
    temp = temp.where(w,0);
    distance = temp.sqrt();

    cdict = {'cos1': cos1,'sin1': sin1,'cos2': cos2,'sin2': sin2,'cos3':cos_raa};
    temp = ee.Image(0).expression('1./cos1 + 1./cos2',cdict);
  
    cdict = {'tan1': tan1,'tan2': tan2,'cos3':cos_raa,'HB':HB,'distance':distance,'sin3':sin_raa,'temp':temp};
    cost = ee.Image(0).expression('HB * sqrt(distance * distance + tan1 * tan1 * tan2 * tan2 * sin3 * sin3) / temp',cdict);
    w = cost.lte(-1);
    cost = cost.where(w,-1);
    w = cost.gte(1);
    cost = cost.where(w,1);
    tvar = cost.acos();
    sint = tvar.sin();
  
    cdict = {'tvar': tvar,'sint': sint,'cost':cost,'pi':np.pi, 'temp':temp};
    overlap = ee.Image(0).expression('(1/pi) * (tvar - sint * cost) * temp',cdict);
    w = overlap.lte(0);
    overlap = overlap.where(w,0).rename(['overlap']);
  
    cdict = {'overlap': overlap,'cosphaang': cosphaang,'cos1':cos1,'cos2':cos2, 'temp':temp};
    li = ee.Image(0).expression('overlap - temp + 0.5 * (1. + cosphaang) / cos1 / cos2',cdict).rename(['li'])
    isotropic = ee.Image(1.0).rename(['isotropic'])
  
    return image.select().addBands(isotropic).addBands(ross).addBands(li)\
            .addBands(image.select('sur_refl_b01').float().multiply(ee.Number(1./10000.)))\
            .addBands(image.select('sur_refl_b02').float().multiply(ee.Number(1./10000.)))\
            .addBands(image.select('sur_refl_b03').float().multiply(ee.Number(1./10000.)))\
            .addBands(image.select('sur_refl_b04').float().multiply(ee.Number(1./10000.)))\
            .addBands(image.select('sur_refl_b05').float().multiply(ee.Number(1./10000.)))\
            .addBands(image.select('sur_refl_b06').float().multiply(ee.Number(1./10000.)))\
            .addBands(image.select('sur_refl_b07').float().multiply(ee.Number(1./10000.))).toFloat()
                
          #.addBands(image).toFloat();

  def getQABits(self,image, start, end, newName):
    pattern = 0;
    #for (var i = start; i <= end; i++) {
    for i in xrange(start,end+1):
       pattern += 2**i 
    #}
    return image.select([0], [newName])\
                  .bitwiseAnd(pattern)\
                  .rightShift(start);

  def maskEmptyPixels(self,image):
    withObs = image.select('num_observations_1km').gt(0);
    return image.updateMask(withObs);

  def maskClouds(self,image):
    getQABits = self.getQABits
    QA = image.select('state_1km');
    land = getQABits(QA, 3, 5, 'land/sea');
    internalCloud = getQABits(QA, 10, 10, 'internal_cloud_algorithm_flag');
    MOD35SnowIce = getQABits(QA, 12, 12, 'MOD35_snow_ice_flag');
    internalSnow = getQABits(QA, 15, 15, 'internal_snow_mask');
    cirrus_detected = getQABits(QA, 8, 9, 'cirrus_detected');
    aerosol_quality = getQABits(QA, 6, 7, 'aerosol_quality');
    valid = image.select("sur_refl_b01").gt(0)\
          .And(image.select("sur_refl_b02").gt(0))\
         .And(image.select("sur_refl_b03").gt(0))\
         .And(image.select("sur_refl_b04").gt(0))\
         .And(image.select("sur_refl_b05").gt(0))\
         .And(image.select("sur_refl_b06").gt(0))\
         .And(image.select("sur_refl_b07").gt(0))\
         .And(image.select("sur_refl_b01").lte(1.1/0.0001))\
         .And(image.select("sur_refl_b02").lte(1.1/0.0001))\
         .And(image.select("sur_refl_b03").lte(1.1/0.0001))\
         .And(image.select("sur_refl_b04").lte(1.1/0.0001))\
         .And(image.select("sur_refl_b05").lte(1.1/0.0001))\
         .And(image.select("sur_refl_b06").lte(1.1/0.0001))\
         .And(image.select("sur_refl_b07").lte(1.1/0.0001));

    nosnow = ((MOD35SnowIce.eq(0)).And(internalSnow.eq(0))).rename(['nosnow']);

    masker = valid.And(nosnow).And(aerosol_quality.gte(1))\
                  .And(cirrus_detected.lte(2))\
                  .And(land.neq(7))\
                  .And(internalCloud.eq(0));
    snow = ((MOD35SnowIce.eq(1)).Or(internalSnow.eq(1)))
    snow = snow.rename(['snow'])
    image = image.select()\
      .addBands(image)\
      .addBands(land.eq(1).rename(['land']))\
      .toFloat();
    image = image.updateMask(masker)
    #image = ee.Algorithms.If(ee.Number(image.sum().gt(0),image,ee.Image(0))
    
    return image

import sys
import getopt

def main():
  # parse command line options
  try:
      opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
  except getopt.error, msg:
      print msg
      print "for help use --help"
      sys.exit(2)
  # process options
  for o, a in opts:
      if o in ("-h", "--help"):
          print __doc__
          sys.exit(0)

  # use eg http://www.latlong.net
  centre = [-17.52,15.42]
  extent = [0.02,0.02]
  options = {'verbose':True,'centre':centre,'extent':extent,\
             'oname':'angola','scale':500,'maxn':100000}
  self = getModisEE(**options)

  try:
    self.load()
  except:
    self.maxn = 10; 
    self.get()
    self.save()

if __name__ == "__main__":
  # execute only if run as a script
  main()
