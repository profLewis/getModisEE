#!/usr/bin/env python

from getModisEE import getModisEE
import sys

# [long,lat]
year = sys.argv[1]

print year

centre = [0.675659,52.438432]
extent = [2.0,2.0]
options = {'verbose':True,'centre':centre,'extent':extent,\
           'oname':'norfolk_%4d'%int(year),'scale':500,'maxn':5,\
           'modis':['MOD09GA','MYD09GA'],'recover':True,\
           'dates':['%4d-01-01'%int(year), '%4d-12-31'%int(year)]}



self = getModisEE(**options)
self.maps = [self.maskEmptyPixels,\
             self.maskClouds,\
             self.makeVariables,\
             self.addTime,\
             self.subtractZero]

self.get()
self.save()

