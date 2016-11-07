#!/usr/bin/env python

from getModisEE import getModisEE
from getModisEE import mapsModisEE 

import sys

# [long,lat]
try:
  year = sys.argv[1]
except:
  year = 2005

print year

centre = [0.675659,52.438432]
extent = [2.0,2.0]
options = {'verbose':True,'centre':centre,'extent':extent,\
           'oname':'coef_norfolk_%4d'%int(year),'scale':500,'maxn':10,\
           'sensors':['MOD09GA','MYD09GA'],'recover':True,\
           'dates':['%4d-01-01'%int(year), '%4d-12-31'%int(year)]}


self = getModisEE(**options)
self.maps = [self.maskEmptyPixels,\
             self.maskClouds,\
             self.makeBRDFKernels,\
             self.addTime,
             self.makeCoefficients]

self.get()
self.save()

