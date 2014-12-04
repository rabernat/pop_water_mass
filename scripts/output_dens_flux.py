from IPython.parallel import Client, interactive
import time

#####################
### This Analysis ###
#####################
hconst = 50. # assume a surface layer of const depth 50 m

##############################
## Set Up Parallel Engines ###
##############################

# give engines time to load
time.sleep(20)

c = Client()
dview = c.direct_view()
lview = c.load_balanced_view()

with dview.sync_imports():
    import numpy
    from watermasstools import pop_model, transformation
    import netCDF4
    import os

# where to find the data
ddir = '/glade/p/ncgd0001/hybrid_v5_rel04_BC5_ne120_t12_pop62/ocn-hist'
fprefix = 'hybrid_v5_rel04_BC5_ne120_t12_pop62.pop.h.nday1'
fnames = []
for year in xrange(47,86):
    for month in xrange(1,13):
        fnames.append(
         '%s/%s.%04d-%02d-01.nc' % (ddir, fprefix, year, month)
        )

# where to output the data
output_dir = '/glade/scratch/rpa/hybrid_v5_rel04_BC5_ne120_t12_pop62'
dview.push(dict(output_dir=output_dir))
              
#############################
## Engine Worker Function ###
#############################
@interactive
def output_dens_flux(pop_fname):
    p = pop_model.POPFile(pop_fname, hconst=50.)
    p.initialize_gradient_operator()
    Nt = len(p.nc.variables['time'])
    
    # create new netCDF file
    output_fname = os.path.join(output_dir, os.path.basename(pop_fname))
    rootgrp = netCDF4.Dataset(output_fname, 'w', format='NETCDF4')
    # create dimensions
    dim_time = rootgrp.createDimension('time', None)
    dim_nlat = rootgrp.createDimension('nlat', p.Ny)
    dim_nlon = rootgrp.createDimension('nlon', p.Nx)

    var_time = rootgrp.createVariable('time','f8',('time',))
    var_time.units = "days since 0000-01-01 00:00:00"
    var_time[:] = p.nc.variables['time'][:]

    mv = p.nc.variables['TLONG'].missing_value
    var_tlong = rootgrp.createVariable('TLONG', 'f8', ('nlat','nlon'), fill_value=mv)
    var_tlong.long_name = "array of t-grid longitudes"
    var_tlong.units = "degrees_east"
    #var_tlong.missing_value = mv
    var_tlong[:] = p.nc.variables['TLONG'][:]

    var_tlat = rootgrp.createVariable('TLAT', 'f8', ('nlat','nlon'))
    var_tlat.long_name = "array of t-grid latitudes"
    var_tlat.units = "degrees_north"
    #var_tlat.missing_value = mv
    var_tlat[:] = p.nc.variables['TLAT'][:]
    
    # for missing values
    mval = numpy.finfo(numpy.float32).max

    # create variables
    varnames = ['rho', 'dens_flux_heat', 'dens_flux_freshwater',
                'dens_flux_mixing', 'dens_flux_cabbeling']
    ncvars = dict()
    for v in varnames:
        myv = rootgrp.createVariable(v, 'f4', 
                        ('time', 'nlat', 'nlon'), fill_value = mval)
        myv.coordinates = "TLONG TLAT time"
        myv.grid_loc = "2110"
        #myv.missing_value = mval
        ncvars[v] = myv

    for n in range(Nt):
        print n
        
        # calculate density fluxes
        rho, Fheat, Fsalt, Fmix, Fcab = p.dens_forcing[n]
        ncvars['rho'][n] = rho.filled(mval)
        ncvars['dens_flux_heat'][n] = Fheat.filled(mval)
        ncvars['dens_flux_freshwater'][n] = Fsalt.filled(mval)
        ncvars['dens_flux_mixing'][n] = Fmix.filled(mval)
        ncvars['dens_flux_cabbeling'][n] = Fcab.filled(mval)

    rootgrp.close()
    return output_fname

#######################
## Apply on Engines ###
#######################

res = lview.map(output_dens_flux, fnames)

while not res.ready():
    print 'progress', res.progress / float(len(res))
    time.sleep(60)

assert res.successful()

