import numpy as np
import netCDF4
#from seawater import eos80
import jmd95

class POPFile(object):
    
    def __init__(self, fname, hmax=None, hconst=None, pref=0., ah=-3e17):
        """Wrapper for POP model netCDF files"""
        self.nc = netCDF4.Dataset(fname)
        self.Ny, self.Nx = self.nc.variables['TAREA'].shape
        self.pref = pref        

        # mask
        self.mask = self.nc.variables['KMT'][:] <= 1
        
        self.alpha = Alpha(self)
        self.beta = Beta(self)
        self.cp = Cp(self)
        self.rho = Rho(self)
        self.dens_forcing = DensForcing(self, hmax=hmax, hconst=hconst, p=pref)
        self.ts_forcing = TSForcing(self, hmax=hmax, hconst=hconst)
        
        self._ah = ah
        
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
        self._cw = np.where( kmt & kmtw, dtw, 0.)
        self._cc = -(self._cn + self._cs + self._ce + self._cw)
        
        # mixing coefficients
        #self._ah = -0.2e20*(1280.0/self.Nx)
        j_eq = np.argmin(self.nc.variables['ULAT'][:,0]**2)
        self._ahf = (tarea / self.nc.variables['UAREA'][j_eq,0])**1.5
        self._ahf[self.mask] = 0.   
        
        # stuff for gradient
        # reciprocal of dx and dy (in meters)
        self._dxtr = 100.*self.nc.variables['DXT'][:]**-1
        self._dytr = 100.*self.nc.variables['DYT'][:]**-1
        self._kmaske = np.where(kmt & kmte, 1., 0.)
        self._kmaskn = np.where(kmt & kmtn, 1., 0.)
                
    def laplacian(self, T):
        return (
            self._cc*T +
            self._cn*np.roll(T,-1,axis=0) +
            self._cs*np.roll(T,1,axis=0) +
            self._ce*np.roll(T,-1,axis=1) +
            self._cw*np.roll(T,1,axis=1)          
        )
        
    def gradient_modulus(self, T):
        dTx = self._kmaske * (np.roll(T,-1,axis=0) - T)
        dTy = self._kmaskn * (np.roll(T,-1,axis=1) - T)
        
        return np.sqrt( 0.5 *
                    (dTx**2 + np.roll(dTx,1,axis=0)**2) * self._dxtr**2
                   +(dTy**2 + np.roll(dTy,1,axis=1)**2) * self._dytr**2
        )        
        
    def biharmonic_tendency(self, T):
        """Caclulate tendency due to biharmonic diffusion of T."""
        d2tk = self._ahf * self.laplacian(T)
        return self._ah * self.laplacian(d2tk)
    
class EOSCalculator(object):
    def __init__(self, parent, p=0., hmax=None, hconst=None):
        self.p = p
        self.hmax = hmax
        self.hconst = hconst
        self.nc = parent.nc
        self.parent = parent
        if self.nc.variables.has_key('SHF'):
            self.hfname = 'SHF'
            self.fwfname = 'SFWF'
            self.mlname = 'HMXL'
        else:
            self.hfname = 'SHF_2'
            self.fwfname = 'SFWF_2'
            self.mlname = 'HMXL_2'




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
# using PSU, kg, m as units
fwflux_factor = 1e-3
rho_fw = 1e3
#fwflux_factor   = 1e-4
#fwflux_factor = 1.  
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
hflux_factor = 1000.0/(rho_sw*cp_sw) / 100.

def get_surface_ts(nc, i):
    try:
        S0 = nc.variables['SSS'].__getitem__(i)
        T0 = nc.variables['SST'].__getitem__(i)
    except KeyError:
        S0 = nc.variables['SALT'][:,0,:,:].__getitem__(i)
        T0 = nc.variables['TEMP'][:,0,:,:].__getitem__(i)    
    return T0, S0

class TSForcing(EOSCalculator):
    def __getitem__(self, i):
        T0, S0 = get_surface_ts(self.nc, i)
        Ffw = self.nc.variables[self.fwfname].__getitem__(i)
        Qhf = self.nc.variables[self.hfname].__getitem__(i)

        if self.hconst is not None:
            H_ml = self.hconst
        else:
            H_ml = self.nc.variables[self.mlname].__getitem__(i)/100.
            if self.hmax is not None:
                H_ml = np.ma.masked_greater(H_ml, self.hmax).filled(self.hmax)
        FT_forc = hflux_factor * Qhf
        FS_forc = salinity_factor * Ffw
        FT_mix = H_ml*self.parent.biharmonic_tendency(T0)
        FS_mix = H_ml*self.parent.biharmonic_tendency(S0)
        
        return [ np.ma.masked_array(F, self.parent.mask) 
                 for F in [T0, S0, FT_forc, FT_mix, FS_forc, FS_mix] ]        

class DensForcing(EOSCalculator):
    def __getitem__(self, i):
        T0, S0 = get_surface_ts(self.nc, i)
        Ffw = self.nc.variables[self.fwfname].__getitem__(i)
        Qhf = self.nc.variables[self.hfname].__getitem__(i)
        if self.hconst is not None:
            H_ml = self.hconst
        else:
            H_ml = self.nc.variables[self.mlname].__getitem__(i)/100.
            if self.hmax is not None:
                H_ml = np.ma.masked_greater(H_ml, self.hmax).filled(self.hmax)
        
        # average the variables if we got multiple time elements
        if isinstance(i, slice):
            T0, S0, Ffw, Qhf = (T0.mean(axis=0), S0.mean(axis=0),
                                Ffw.mean(axis=0), Qhf.mean(axis=0))
            if isinstance(H_ml, np.ndarray):
                H_ml = H_ml.mean(axis=0)

        if self.p == 0.:
            rho, drhodT, drhodS = jmd95.eos.state_surface(T0, S0)
        else:
            rho, drhodT, drhodS = jmd95.eos.state(self.p, T0, S0)
        Fdens_heat = drhodT * hflux_factor * Qhf
        Fdens_salt = drhodS * salinity_factor * Ffw

        # mixing and cabbelling
        mixing_rho = H_ml * self.parent.biharmonic_tendency(rho)
        mixing_T = H_ml * self.parent.biharmonic_tendency(T0)
        mixing_S = H_ml * self.parent.biharmonic_tendency(S0)
        cab = drhodT*mixing_T + drhodS*mixing_S - mixing_rho
        
        return [ np.ma.masked_array(fld, self.parent.mask) for fld in
                 [rho, Fdens_heat, Fdens_salt, mixing_rho, cab ] ]

        #return (np.ma.masked_array(rho, self.parent.mask),
        #        np.ma.masked_array(Fdens_heat, self.parent.mask),
        #        np.ma.masked_array(Fdens_salt, self.parent.mask),
        #        np.ma.masked_array(Fdens_mix, self.parent.mask))
        
        #return -alpha*Qhf/cp, rho0*beta*Fsalt*S0/(1 - S0)
        #return -alpha*Qhf/cp, rho0*beta*Ffw/salinity_factor
        #return -alpha*Qhf/cp, rho0*beta*salinity_factor*Ffw

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
        T0, S0 = get_surface_ts(self.nc, i)
        
        # average the variables if we got multiple time elements
        if isinstance(i, slice):
            T0, S0, = T0.mean(axis=0), S0.mean(axis=0)
        if self.p == 0.:
            rho, drhodT, drhodS = jmd95.eos.state_surface(T0, S0)
        else:
            rho, drhodT, drhodS = jmd95.eos.state(self.p, T0, S0)
        return rho

        
