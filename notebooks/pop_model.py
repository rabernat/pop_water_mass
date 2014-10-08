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

class DensFlux(EOSCalculator):
    def __getitem__(self, i):
        S0 = self.nc.variables['SSS'].__getitem__(i)
        T0 = self.nc.variables['SST'].__getitem__(i)
        Fsalt = self.nc.variables['SFWF_2'].__getitem__(i)
        Qheat = self.nc.variables['SHF_2'].__getitem__(i)

        alpha = eos80.alpha(S0, T0, self.p, pt=True)
        beta = eos80.alpha(S0, T0, self.p, pt=True) 
        cp = eos80.cp(S0, T0, self.p)
        rho0 = eos80.dens0(S0, T0)
        
        return -alpha*Qheat/cp, rho0*beta*Fsalt*S0/(1 - S0)

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


        