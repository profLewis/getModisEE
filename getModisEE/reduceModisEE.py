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

class reduceModisEE():
  '''
  reduce
  '''         
  def conOne(self):
    pass


