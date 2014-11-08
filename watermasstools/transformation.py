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
    
    def __init__(self, rholevs=None, basin_names=None, area=1.,
                 lonmin=None, lonmax=None, latmin=None, latmax=None):

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
        
        mask = pmodel.mask.copy()   
        
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
        
    def calculate_rholevs(self, rho=None, nlevs=100, linear=False,
                            rhomin=None, rhomax=None):
        if rho is None:
            assert linear
            assert rhomin is not None
            assert rhomax is not None
        else:
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
            
    def set_rholevs(self, rholevs):
        self.rholevs = rholevs
            
    def calculate_transformation_rate(self, rho, *args):
        if self.mask is None:
            raise Exception('Mask is not initialized.')
        return np.array(sum_inside_contours(
            self.mask_field(rho), self.rholevs,
            args, area=self.area))[:,1:] / np.diff(self.rholevs)

class TSWaterMassRegion(WaterMassRegion):
    
    def __init__(self, Tlevs, Slevs, **kwargs):
        """Calculate transformation rates in regions *bounded* by
        Tlevs and Slevs. The output of the calculations will have
        dimensions of len(Tlevs)-1, len(Slevs)-1. The center-point
        coordinates can be accessed via the .Tc and .Sc variables.
        """
        self.Tlevs = Tlevs
        self.Slevs = Slevs
        self.DT = np.diff(self.Tlevs)
        self.DS = np.diff(self.Slevs)
        self.DTDS = (self.DT[np.newaxis,:,np.newaxis] * 
                     self.DS[np.newaxis,np.newaxis,:])
        
        # center point coordinates
        self.Tc = self.Tlevs[:-1] + 0.5*self.DT
        self.Sc = self.Slevs[:-1] + 0.5*self.DS
        
        # cell edge coordinates (for pcolor)
        self.S, self.T = np.meshgrid(Slevs,Tlevs)
        WaterMassRegion.__init__(self, **kwargs)
        
    def calculate_transformation_rate(self, T, S, Tforcing, Sforcing):
        assert isinstance(Tforcing, list)
        assert isinstance(Sforcing, list)
        assert len(Tforcing) == len(Sforcing)
        
        Nforc = len(Tforcing)
        if self.mask is None:
            raise Exception('Mask is not initialized.')

        # it should be faster if we do T and S together
        F = Tforcing + Sforcing
        A = sum_inside_contours_2D(
            self.mask_field(S), self.mask_field(T),
            self.Slevs, self.Tlevs, F, area=self.area)
            
        # discarding the first point means discarding all points with
        # values lower than Tlevs[0] and Slevs[0]. This seems consistent
        # with the API.
        AT = np.array(A[:Nforc])[:,1:,1:] / self.DTDS
        AS = np.array(A[Nforc:])[:,1:,1:] / self.DTDS
        
        return AT, AS
   
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
        
def sum_inside_contours_2D(idx1_field, idx2_field, idx1_levels, idx2_levels,
                           fields, area=1.):
    shape = idx1_field.shape
    assert idx2_field.shape == shape
    N1 = len(idx1_levels)
    N2 = len(idx2_levels)
    labels1 = np.digitize(idx1_field.ravel(), idx1_levels)
    labels2 = np.digitize(idx2_field.ravel(), idx2_levels)
    labels = labels1 + N1*labels2
    labels.shape = shape
    
    res = []
    # wrap it in a list if we only got one argument
    if isinstance(fields, np.ndarray):
        fields = [fields,]
    for field in fields:
        assert field.shape == shape
        field = np.ma.masked_array(field, idx1_field.mask)
        ndsum = ndimage.sum(
            (area*field).filled(0.),
            labels=labels, index=np.arange(N1*N2))
        ndsum.shape = N2, N1
        res.append(ndsum)
    if len(res)==1:
        return res[0]
    else:
        return res
    
    
    
########### Grid for T/S Transformation Calculation ########
#
#
#  T_1+1/2 x-------v_2,0-------x-------v_2,1-------x
#          |                   |                   |
#          |                   |                   |
#          |                   |                   |
#  T_1    u_1,0      X        u_1,1      X       u_1,2
#          |                   |                   |
#          |                   |                   |         
#          |                   |                   |
#  T_1/2   x-------v_1,0-------x-------v_1,1---------x
#          |                   |                   |
#          |                   |                   |
#          |->                 |                   |
#  T_0    u_0,0      X        u_0,1      X       u_0,2
#          |->                 |                   |
#          |                   |                   |
#          |        ^^^        |                   |
#  T_-1/2  x-------v_0,0-------x-------v_0,1---------x
#
#        S_-1/2     S_0      S_1/2      S_1      S_1+1/2 
#
#
###########################################################
#
#  We actually want the T and S transformation rates at different
#  points on this grid? Or not...
#
#
#



    
        
