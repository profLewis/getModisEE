# getModisEE
A class for pulling MODIS MOD09 and MYD09 from Google Earth Engine to the local computer
for a set geographic region (in lat/lon) and time period (default all)

The class uses EE to pre-process the dataset by masking clouds etc
and calculated angular BRDF kernels from the satellite angle information.

It loads the dataset to a dictionary of arrays (self.data) and also can store it in a pickle file (in odir by default)

If the pickle files exists, you can load that with self.load()

It has only limited documentation at present and the functionality could be extended if required.

Also, should add a grid of geographic locations for pixels.

Also, should keep track of which sample we are on to be tolerant to failure.

This is not a super-speedy way of downloading MODIS data, but it is convenient for grabbing subsets and applying a consistent masking. It is also convenient for further processing to BRDF information.

# install

First, get Google Earth Engine, see: [https://developers.google.com/earth-engine/python_install](https://developers.google.com/earth-engine/python_install)

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
           extent = [2.0,2.0]
           
           year = 2016
           
           options = {'verbose':True,'centre':centre,'extent':extent,\
                      'oname':'norfolk_%4d'%int(year),'scale':500,'maxn':2000,\
                      'sensors':['MOD09GA','MYD09GA'],'recover':True,\
                      'dates':['%4d-01-01'%int(year), '%4d-12-31'%int(year)]}



           self = getModisEE(**options)
           self.maps = [self.maskEmptyPixels,\
                        self.maskClouds,\
                        self.makeBRDFKernels,\
                        self.addTime]

           self.get()
           self.save()
           
Note that we have the 'recover' flag set True, which means that the process will attempt to recover from a previous attempt at downloading the data.
           
You can check what the current dataset looks like, e.g. with:

           import gdal
           f = gdal.Open('odir/norfolk_000000.sur_refl_b02.tif').ReadAsArray()
           plt.clf()
           plt.imshow(f,interpolation='nearest')
           plt.colorbar()
           
Note that if `sur_refl_b01` is all zero values (or rather, if sum is zero) then we ignore the dataset.

Note also that, at present, we do not store the coordinate information. This is however accessible from the tif files.

You will find a temporary dump file in the running directory called e.g. `download.3l3c3f.tmp`. This contains a data dump every e.g. 50 samples. The frequency of dump is controlled by `dumpFreq=50` in `options`. In case of a bad recovery, you can replace the output data file by the dump file.

You can also use the script [getData.py](getData.py) to access data for a particular year. If you want a more complex parser, you could build one.

In practice, it seems you can run around 5 processes simultaneously, i.e. set 5 lots of [getData.py](getData.py) going with different years, at the same time. For a 2 degree by 2 degree area, this takes around one hour per year to process and download. e.g.:

           python getData.py 2001
           
or just:

           ./getData.py 2001

 The output files would be of the order of 8GB. If they end up significantly less than this, check the dates of the downloaded data. Ocassionaly, EE can just fail on you, so some error checking would be of value.  

# Example use #2

           options = {'verbose':True,'centre':centre,'extent':extent,\
                      'oname':'angola','scale':500,'maxn':100000}
           self = getModisEE(**options)
           self.load()
