import numpy as np
import pop_model
from scipy import ndimage

POP_REGIONS = {
    'Southern Ocean': 1,
    'Pacific Ocean': 2,
    'Indian Ocean': 3,
    'Red Sea': 4,
    'Atlantic Ocean': 6,
    'Mediterranean Sea': 7,
    'Labrador Sea': 8,
    'GIN Seas': 9,
    'Arctic Ocean': 10
}

class WaterMassRegion(object):
    """For defining water mass transformation regions"""
    
    def __init__(self, name, rholevs=None, basin_names=None, area=1.,
                 lonmin=None, lonmax=None, latmin=None, latmax=None):
        self.name = name
        self.rholevs = rholevs
        self.basin_names = basin_names
        self.area = 1.
        self.lonmin, self.lonmax = lonmin, lonmax
        self.latmin, self.latmax = latmin, latmax
        self.mask = None
                
    def initialize_mask(self, pmodel, extra_mask=None):
        """Create a store a mask for the region based on the shape."""
        assert isinstance(pmodel, pop_model.POPFile)
        
        lon = pmodel.nc.variables['TLONG'][:]
        lat = pmodel.nc.variables['TLAT'][:]  
        
        mask = pmodel.mask
        
        if self.basin_names is not None:
            basin_mask = np.zeros_like(mask, dtype='bool')
            if isinstance(self.basin_names, str):
                self.basin_names = [self.basin_names,]
            for bname in self.basin_names:
                bid = POP_REGIONS[bname]
                basin_mask += (pmodel.nc.variables['REGION_MASK'][:]==bid)
            mask += (~basin_mask)
        
        if self.lonmin is not None:
            mask += (lon <= self.lonmin)
        if self.latmin is not None:
            mask += (lat <= self.latmin)
        if self.lonmax is not None:
            mask += (lon > self.lonmax)
        if self.latmax is not None:
            mask += (lat > self.latmax)
        if extra_mask is not None:
            mask += extra_mask
            
        self.mask = mask
        self.area = self.mask_field(pmodel.nc.variables['TAREA'][:]/1e4)
        
    def mask_field(self, field):
        if self.mask is None:
            raise Exception('Mask is not initialized.')
        return np.ma.masked_array(field, self.mask)
        
    def calculate_rholevs(self, rho, nlevs=100, linear=False,
                            rhomin=None, rhomax=None):
        rho = self.mask_field(rho)
        if rhomin is None:
            rhomin = rho.min()
        if rhomax is None:
            rhomax = rho.max()
        rholin = np.linspace(rhomax, rhomin, nlevs)
        if linear:
            self.rholevs = rholin
        else:
            acum = np.cumsum(
              sum_inside_contours(rho, rholin, self.area))
            alevs = np.linspace(acum.min(), acum.max(), nlevs)
            self.rholevs = np.interp(alevs, acum, rholin)
            
    def calculate_transformation_rate(self, rho, *args):
        if self.mask is None:
            raise Exception('Mask is not initialized.')
        if self.mask is None:
            raise Exception('Mask is not initialized.')
        return np.array(sum_inside_contours(
            self.mask_field(rho), self.rholevs,
            args, area=self.area))[:,1:] / np.diff(self.rholevs)

def sum_inside_contours(index_field, index_levels, fields, area=1.):
    """Sum the arrays in args within countours
    of index_field. Bins are specifiied by index levels.
    index_levels must be monotonically increasing or decreasing."""
    di = np.diff(index_levels)
    sign = np.sign(di[0])
    assert np.all(np.sign(di)==sign)
    shape = index_field.shape
    labels = np.digitize(index_field.ravel(), index_levels,
                         right=(sign==1))
    labels.shape = shape
    
    res = []
    
    # wrap it in a list if we only got one argument
    if isinstance(fields, np.ndarray):
        fields = [fields,]
    for field in fields:
        assert field.shape == shape
        field = np.ma.masked_array(field, index_field.mask)
        res.append( ndimage.sum(
            (area*field).filled(0.),
            labels=labels, index=np.arange(len(index_levels))
        ))
    if len(res)==1:
        return res[0]
    else:
        return res
        