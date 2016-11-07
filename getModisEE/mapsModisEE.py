#!/usr/bin/env python

'''
pull MODIS MOD09/MYD09 data from Google Earth Engine
including filtering and calculation of angular kernels


Here, we define a library of maps

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

class mapsModisEE():
  '''
  maps
  '''         
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
    #properties = image.propertyNames()
    #print properties

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

