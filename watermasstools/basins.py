import numpy as np
from scipy.interpolate import RegularGridInterpolator, NearestNDInterpolator
import netCDF4

class BasinIdentifier(object):
    
    def __init__(self):
        self.nc = netCDF4.Dataset(
          'http://data.nodc.noaa.gov/' +
          'thredds/dodsC/woa/WOA05nc/other/basin.nc')
          
        lon = self.nc.variables['lon'][:]
        lat = self.nc.variables['lat'][:]
        basin = self.nc.variables['basin'][0,0]

        # wrap in longitude
        lon = np.hstack([2*lon[0]-lon[1], lon, 2*lon[-1]-lon[-2] ])
        basin = np.ma.hstack(
            [basin[:,-1][:,np.newaxis], basin, basin[:,0][:,np.newaxis]])
            
        xx,yy = np.meshgrid(lon,lat)
        xx = np.ma.masked_array(xx, basin.mask)
        yy = np.ma.masked_array(yy, basin.mask)
        points = np.vstack([xx.compressed(), yy.compressed()]).T
        values = basin.compressed()
        print points.shape
        print values.shape
        
        #self.interpolator = RegularGridInterpolator(
        #                      (lat,lon), basin, bounds_error=False,
        #                      fill_value=None, method='nearest')
        self.interpolator = NearestNDInterpolator(
            points, values)

    def __call__(self, lon, lat):
        assert lon.shape == lat.shape
        shape = lon.shape
        #res = self.interpolator((lat.ravel(), lon.ravel()), method='nearest')
        res = self.interpolator((lon.ravel(), lat.ravel()))
        res.shape = shape
        return res


def get_basin_id(name):
    return BASINS.index(name)+1
                              
BASINS = [
    'Atlantic Ocean',
    'Pacific Ocean',
    'Indian Ocean',
    'Mediterranean Sea',
    'Baltic Sea',
    'Black Sea',
    'Red Sea',
    'Persian Gulf',
    'Hudson Bay',
    'Southern Ocean',
    'Arctic Ocean',
    'Sea of Japan',
    'Kara Sea',
    'Sulu Sea',
    'Baffin Bay',
    'East Mediterranean',
    'West Mediterranean',
    'Sea of Okhotsk',
    'Banda Sea',
    'Caribbean Sea',
    'Andaman Basin',
    'North Caribbean',
    'Gulf of Mexico',
    'Beaufort Sea',
    'South China Sea',
    'Barents Sea',
    'Celebes Sea',
    'Aleutian Basin',
    'Fiji Basin',
    'North American Basin',
    'West European Basin',
    'Southeast Indian Basin',
    'Coral Sea',
    'East Indian Basin',
    'Central Indian Basin',
    'Southwest Atlantic Basin',
    'Southeast Atlantic Basin',
    'Southeast Pacific Basin',
    'Guatemala Basin',
    'East Caroline Basin',
    'Marianas Basin',
    'Philippine Sea',
    'Arabian Sea',
    'Chile Basin',
    'Somali Basin',
    'Mascarene Basin',
    'Crozet Basin',
    'Guinea Basin',
    'Brazil Basin',
    'Argentine Basin',
    'Tasman Sea',
    'Atlantic Indian Basin',
    'Caspian Sea',
    'Sulu Sea II',
    'Venezuela Basin',
    'Bay of Bengal',
    'Java Sea',
    'East Indian Atlantic Basin' ]
