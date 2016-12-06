#!/usr/bin/env python

'''
pull MODIS MOD09/MYD09 data from Google Earth Engine
including filtering and calculation of angular kernels

'''

try:
  import ee
except:
  print 'you need to install the google earth engine packages'

from linearBRDFBase import linearBRDFBase
from mapsModisEE import mapsModisEE 
from reduceModisEE import reduceModisEE

import wget
import zipfile
import gdal

import os
import numpy as np
import sys

unload = lambda self,name,kwargs: ((name in kwargs and kwargs[name]) or \
                              (hasattr(self,name)) and getattr(self,name))

class getModisEE(linearBRDFBase,mapsModisEE,reduceModisEE):
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
                                          self.makeBRDFKernels,self.addTime]
    self.sensors = unload(self,'sensors',kwargs) or ['MOD09GA','MYD09GA']
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

  def addSensor(self,image):
    return image.addBands();

  def getCollections(self,*args,**kwargs):
    '''
    Specify datasets and time period
    '''
    self.parser(*args,**kwargs)
    collection = ee.ImageCollection(self.sensors[0]).filterDate(self.dates[0],self.dates[1])
    for i in self.sensors[1:]:
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

  def pullData(self,image,count,clean=True,centre=None,extent=None,show=True,pull=False):
 
    if self.verbose: print '##'*8,'\n',count,'\n','##'*8 
    # force calculation of topLeft and bottomRight
    # 
    if (centre) or (extent) or \
       (not hasattr(self, 'topLeft')) or (not hasattr(self, 'bottomRight')):
      self.setAOI(centre,extent)
    definedAoi  = ee.Geometry.Rectangle(self.topLeft[0],self.topLeft[1],\
                                    self.bottomRight[0],self.bottomRight[1])

    boundingBox = definedAoi.bounds(1)
    region = ee.Geometry(boundingBox.getInfo())\
                        .toGeoJSONString()

    image = image.clip(definedAoi)
    if show:
      print 1
    if pull:
      self.pullData2(image,count,clean=clean,centre=centre,extent=extent)

 
  def pullData2(self,image,count,clean=True,centre=None,extent=None):
 
    if self.verbose: print '##'*8,'\n',count,'\n','##'*8 
    # force calculation of topLeft and bottomRight
    # 
    if (centre) or (extent) or \
       (not hasattr(self, 'topLeft')) or (not hasattr(self, 'bottomRight')):
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
      'name':'%s'%(self.oname),\
      'crs': 'EPSG:4326',\
      'scale': self.scale,\
      'region':region\
    })

    if self.verbose: 
      print(url)
      print self.centre,self.extent
      #print region

    # op
    filename = wget.download(url, bar=wget.bar_thermometer)

    zf = zipfile.ZipFile(filename, 'r')
    data = {}
    # ensure odir exists
    self.mkdir(self.odir)

    # look for bad data
    badData = False
    for name in zf.namelist():
      if self.verbose: print name,
      if name.split('.')[-1] == 'tif':
        # its a tif file
        f = open(self.odir + os.sep + name,'w+b')
        f.write(zf.read(name))
        f.close()
        data[name] = gdal.Open(self.odir + os.sep + name).ReadAsArray()
        if (('sur_refl' in name) or ('isotropic' in name)) \
					and data[name].sum() == 0:
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
      self.data['count'] = count
      for k in data.keys():
        if data[k].ndim == 1:
          self.data[k] = np.atleast_2d(data[k])
        else:
          self.data[k] = np.atleast_3d(data[k].T).T

    else:
      self.data['count'] = count
      # append data from each key
      # using numpy hstack
      for k in data.keys():
        # make sure all lists
        try:
          if data[k].ndim == 2:
            # try to stack on array, if not load
            try:
              self.data[k] = np.vstack([self.data[k],np.atleast_3d(data[k].T).T])
            except:
              self.data[k] = np.atleast_3d(data[k].T).T

          elif data[k].ndim == 1:
            # try to stack on array, if not load
            try:
              self.data[k] = np.vstack([self.data[k],np.atleast_2d(data[k])])
            except:
              self.data[k] = np.atleast_2d(data[k])

        except:
          self.data[k] = data[k] 
    return data

  def get(self,*args,**kwargs):
    '''
    download and load into self.data
    '''
    self.getCollections(*args,**kwargs)
    dumper = 'dump_%s.tmp'%str(os.getpid())
    if self.verbose: print 'dumper file',dumper

    # try to recover first?
    start = 0
    if self.recover:
      try:
        self.load()
        start = self.data['count'] + 1
      except:
        start = 0

    if self.verbose: print 'starting at',start

    try:
      for i in xrange(start,self.maxn+start):
        if self.verbose: print i
        self.pullData(ee.ImageCollection([self.collection.toList(1,i).get(-1)]).min(),i)
        try:
          if i%self.dumpFreq == 0:
            if self.verbose:
              print 'temp dump at count %d to'%i,dumper
            self.save(dumper)
        except:
          pass
    except:
      if self.verbose: print 'stopping at',i
      try:
        self.data['count'] = i
      except:
        pass
      pass
    self.data['count'] = i

  def load(self,*args,**kwargs):
    '''
    load dataset as pickle for now
    '''
    '''
    TODO:
      check we can mmap, else use other object
    '''
    import pickle
    ofile = (len(args) and type(args[0] == 'str') and args[0]) or \
            (('pfile' in kwargs) and kwargs['pfile']) or \
            (self.odir + os.sep + self.oname + '.pkl')
    output = open(ofile, 'rb') 
    self.data = pickle.load(output)
    if not 'count' in self.data:
      self.data['count'] = 0
    output.close()

  def save(self,*args,**kwargs):
    '''
    Save dataset as pickle for now
    '''
    import pickle
    if not hasattr(self,'data'):
      self.data = {}
      self.data['count'] = 0

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
      self.getCollections()
    self.pullData(ee.ImageCollection([self.collection.toList(1,count).get(-1)]).min(),0)
    self.counter += 1
    return self.data


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
