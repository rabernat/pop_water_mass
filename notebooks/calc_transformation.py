from IPython.parallel import Client

c = Client()
dview = c.direct_view()
lview = c.load_balanced_view()

with dview.sync_imports():
    import numpy
    from watermasstools import pop_model, transformation

dview.push('data_dir', data_dir)
    
# define basins
natl = transformation.WaterMassRegion('natl',
                    basin_names=['Atlantic Ocean', 'GIN Seas', 'Labrador Sea'], latmin=15)
natl.initialize_mask(p)
natl.set_rholevs(np.load('../data/rholevs_natl.npy'))
# push to engines
dview.push('region_dict', {'natl': natl})
    
def calc_transformation_rates(pop_fname):
    p = pop_model.POPFile(ddir + fname)
    p.initialize_gradient_operator()
    Nt = len(p.nc.variables['time'])
    
    result = dict()
    for rname in region_dict:
        result[rname] = 0.
    
    for n in range(Nt):
        F_temp, F_salt = p.dens_flux[n]
        
        H_ml = p.nc.variables['HMXL_2'][0]/100.
        rho0 = p.rho[n]
        bt = np.ma.masked_invalid(
            p.biharmonic_tendency(p.nc.variables['SST'][0]))
        bt.mask += mask
        
        bs = np.ma.masked_invalid(
            p.biharmonic_tendency(p.nc.variables['SSS'][0]))
        bs.mask += mask
        D_salt = np.ma.masked_array(p.beta[n] * rho0 * H_ml * bs, mask)
        D_temp = np.ma.masked_array(-p.alpha[n] * rho0 * H_ml * bt, mask)
        
        for rname, reg in region_dict.iteritems():
            A = reg.calculate_transformation_rate(rho0,
                              F_temp, F_salt, D_temp, D_salt)
            result[rname] += A

    for rname in region_dict:
        result[rname] ./ Nt
        
    return result
    
