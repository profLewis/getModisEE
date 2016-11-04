# getModisEE
A class for pulling MODIS MOD09 and MYD09 from Google Earth Engine to the local computer
for a set geographic region (in lat/lon) and time period (default all)

The class uses EE to pre-process the dataset by masking clouds etc
and calculated angular BRDF kernels from the satellite angle information.

It loads the dataset to an array (self.data) and also can store it in a pickle file (in odir by default)

If the pickle files exists, you can load that with self.load()

It has only limited documentation at present and the functionality could be extended if required.

Also, should add a grid of geographic locations for pixels.

Also, should keep track of which sample we are on to be tolerant to failure.

# install

           python setup.py install

# Example use #1


           from getModisEE import getModisEE
           centre = [-17.52,15.42]
           extent = [0.02,0.02]
           options = {'verbose':True,'centre':centre,'extent':extent,\
                      'oname':'angola','scale':500,'maxn':100000}
           self = getModisEE(**options)
           self.get()
           self.save()

# Example use #2

           options = {'verbose':True,'centre':centre,'extent':extent,\
                      'oname':'angola','scale':500,'maxn':100000}
           self = getModisEE(**options)
           self.load()
