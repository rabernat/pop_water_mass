from IPython.parallel import Client

c = Client()
dview = c.direct_view()
lview = c.load_balanced_view()

with dview.sync_imports():
    import numpy
    from poputils import pop_model, transformation
    
# define basins and push
dview.push('data_dir', data_dir)
dview.push('region_dict', region_dict)
    
def calc_transformation_rates(pop_fname):
    p = pop_model.POPFile(ddir + fname)
    p.initialize_gradient_operator()
    Nt = len(p.nc.variables['time'])
    
    result = dict()
    
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
        
        for r in region_dict:
            region_mask = region_dict[r]['mask']
            rholevs = region_dict[r]['rholevs']
            rho_reg = np.ma.masked_array(rho0, region_mask)
            result[r] += np.array(
                transformation.sum_inside_contours(
                rho_reg, rholevs,
                [Fheat, Fsalt, F_Tmix, F_Smix],
                area=area)
            )
    
