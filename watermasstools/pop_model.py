import numpy as np
import netCDF4
from seawater import eos80

class POPFile(object):
    
    def __init__(self, fname):
        """Wrapper for POP model netCDF files"""
        self.nc = netCDF4.Dataset(fname)
        self.Ny, self.Nx = self.nc.variables['TAREA'].shape
        
        # mask
        self.mask = self.nc.variables['KMT'][:] <= 1
        
        self.alpha = Alpha(self)
        self.beta = Beta(self)
        self.cp = Cp(self)
        self.rho = Rho(self)
        self.dens_flux = DensFlux(self)
        
    def initialize_gradient_operator(self):
        # raw grid geometry
        work1 = (self.nc.variables['HTN'][:] /
                 self.nc.variables['HUW'][:])
        tarea = self.nc.variables['TAREA'][:]
        tarea_r = np.ma.masked_invalid(tarea**-1).filled(0.)
        dtn = work1*tarea_r
        dts = np.roll(work1,-1,axis=0)*tarea_r
        
        work1 = (self.nc.variables['HTE'][:] /
                 self.nc.variables['HUS'][:])
        dte = work1*tarea_r
        dtw = np.roll(work1,-1,axis=1)*tarea_r
        
        # boundary conditions
        kmt = self.nc.variables['KMT'][:] > 1
        kmtn = np.roll(kmt,-1,axis=0)
        kmts = np.roll(kmt,1,axis=0)
        kmte = np.roll(kmt,-1,axis=1)
        kmtw = np.roll(kmt,1,axis=1)
        self._cn = np.where( kmt & kmtn, dtn, 0.)
        self._cs = np.where( kmt & kmts, dts, 0.)
        self._ce = np.where( kmt & kmte, dte, 0.)
        self._cw = np.where( kmt & kmte, dtw, 0.)
        self._cc = -(self._cn + self._cs + self._ce + self._cw)
        
        # mixing coefficients
        self._ah = -0.2e20*(1280.0/self.Nx)
        j_eq = np.argmin(self.nc.variables['ULAT'][:,0]**2)
        self._ahf = (tarea / self.nc.variables['UAREA'][j_eq,0])**1.5
    
    def laplacian(self, T):
        return (
            self._cc*T +
            self._cn*np.roll(T,-1,axis=0) +
            self._cs*np.roll(T,1,axis=0) +
            self._ce*np.roll(T,-1,axis=1) +
            self._cw*np.roll(T,1,axis=1)          
        )
    
    def biharmonic_tendency(self, T):
        """Caclulate tendency due to biharmonic diffusion of T."""
        d2tk = self._ahf * self.laplacian(T)
        return self._ah * self.laplacian(d2tk)
    
class EOSCalculator(object):
    def __init__(self, parent, p=0.):
        self.p = p
        self.nc = parent.nc
        self.parent = parent


# from forcing.F90
# !-----------------------------------------------------------------------
# !
# !  convert fresh water flux (kg/m^2/s) to virtual salt flux (msu*cm/s):
# !  --------------------------------------------------------------------
# !    ocean reference salinity in (o/oo=psu)
# !    density of freshwater rho_fw = 1.0 (g/cm^3)
# !    h2o_flux in (kg/m^2/s) = 0.1 (g/cm^2/s)
# !
# !    salt_flux            = - h2o_flux * ocn_ref_salinity / rho_fw
# !    salt_flux (msu*cm/s) = - h2o_flux (kg/m^2/s)
# !                           * ocn_ref_salinity (psu)
# !                           * 1.e-3 (msu/psu)
# !                           * 0.1 (g/cm^2/s)/(kg/m^2/s)
# !                           / 1.0 (g/cm^3)
# !                         = - h2o_flux (kg/m^2/s)
# !                           * ocn_ref_salinity (psu)
# !                           * fwflux_factor (cm/s)(msu/psu)/(kg/m^2/s)
# !
# !    ==>  fwflux_factor = 1.e-4
# !
# !    salt_flux(msu*cm/s) = h2oflux(kg/m^2/s) * salinity_factor
# !
# !    ==> salinity_factor = - ocn_ref_salinity(psu) * fwflux_factor
# !
# !-----------------------------------------------------------------------
#
ocn_ref_salinity = 34.7
#fwflux_factor   = 1e-4
fwflux_factor = 1.
salinity_factor = -ocn_ref_salinity*fwflux_factor

# !-----------------------------------------------------------------------
# !
# !  convert heat, solar flux (W/m^2) to temperature flux (C*cm/s):
# !  --------------------------------------------------------------
# !    heat_flux in (W/m^2) = (J/s/m^2) = 1000(g/s^3)
# !    density of seawater rho_sw in (g/cm^3)
# !    specific heat of seawater cp_sw in (erg/g/C) = (cm^2/s^2/C)
# !
# !    temp_flux          = heat_flux / (rho_sw*cp_sw)
# !    temp_flux (C*cm/s) = heat_flux (W/m^2)
# !                         * 1000 (g/s^3)/(W/m^2)
# !                         / [(rho_sw*cp_sw) (g/cm/s^2/C)]
# !
# !                       = heat_flux (W/m^2)
# !                         * hflux_factor (C*cm/s)/(W/m^2)
# !
# !    ==>  hflux_factor = 1000/(rho_sw*cp_sw)
# !
# !-----------------------------------------------------------------------

cp_sw = 3.996e7
rho_sw = 4.1/3.996
hflux_factor = 1000.0/(rho_sw*cp_sw)

class DensFlux(EOSCalculator):
    def __getitem__(self, i):
        S0 = self.nc.variables['SSS'].__getitem__(i)
        T0 = self.nc.variables['SST'].__getitem__(i)
        Ffw = self.nc.variables['SFWF_2'].__getitem__(i)
        Qhf = self.nc.variables['SHF_2'].__getitem__(i)

        alpha = eos80.alpha(S0, T0, self.p, pt=True)
        beta = eos80.alpha(S0, T0, self.p, pt=True) 
        cp = eos80.cp(S0, T0, self.p)
        rho0 = eos80.dens0(S0, T0)
        
        #return -alpha*Qhf/cp, rho0*beta*Fsalt*S0/(1 - S0)
        return -alpha*Qhf/cp, rho0*beta*Ffw/salinity_factor

class Alpha(EOSCalculator):
    def __getitem__(self, i):
        """Calculate alpha from SST and SSS"""
        return eos80.alpha(
            self.nc.variables['SSS'].__getitem__(i),
            self.nc.variables['SST'].__getitem__(i),
            self.p, pt=True)
        
class Beta(EOSCalculator):
    def __getitem__(self, i):
        """Calculate beta from SST and SSS"""
        return eos80.beta(
            self.nc.variables['SSS'].__getitem__(i),
            self.nc.variables['SST'].__getitem__(i),
            self.p, pt=True)

class Cp(EOSCalculator):
    def __getitem__(self, i):
        """Calculate Cp from SST and SSS"""
        return eos80.cp(
            self.nc.variables['SSS'].__getitem__(i),
            self.nc.variables['SST'].__getitem__(i),
            self.p)

class Rho(EOSCalculator):
    def __getitem__(self, i):
        """Calculate Cp from SST and SSS"""
        return eos80.dens0(
            self.nc.variables['SSS'].__getitem__(i),
            self.nc.variables['SST'].__getitem__(i))


        