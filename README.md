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

First, get Google Earth Engine sorted:

           see: https://developers.google.com/earth-engine/python_install

Then

           pip install -e git+https://github.com/profLewis/getModisEE.git#egg=getModisEE
           
or

           git clone https://github.com/profLewis/getModisEE.git
           cd getModisEE
           python setup.py install

# Example use #1


           from getModisEE import getModisEE
           
           # [long,lat]
           centre = [0.675659,52.438432]
           extent = [1.0,1.0]
           options = {'verbose':True,'centre':centre,'extent':extent,\
                      'oname':'norfolk','scale':500,'maxn':100000}
           self = getModisEE(**options)
           
           loadOptions = {'modis':['MOD09GA','MYD09GA'],\
                          'dates':['2001-01-01', '2020-01-01'],\
                          'maps':[self.maskEmptyPixels,self.maskClouds,\
                                          self.makeVariables,self.addTime,\
                                          self.subtractZero]}
           self.get(**loadOptions)
           self.save()
           
You can check what the current dataset looks like, e.g. with:

           import gdal
           f = gdal.Open('odir/norfolk_000000.sur_refl_b02.tif').ReadAsArray()
           plt.clf()
           plt.imshow(f,interpolation='nearest')
           plt.colorbar()
           
Note that if `sur_refl_b01` is all zero values (or rather, if sum is zero) then we ignore the dataset.

# Example use #2

           options = {'verbose':True,'centre':centre,'extent':extent,\
                      'oname':'angola','scale':500,'maxn':100000}
           self = getModisEE(**options)
           self.load()
