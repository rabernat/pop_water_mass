from IPython.parallel import Client
import time
import numpy
from watermasstools import pop_model, transformation
from IPython.parallel import interactive

#############################
## Engine Worker Function ###
#############################

@interactive
def calc_transformation_rates(pop_fname):
    """Worker function for transformation analysis.
    Takes a single file name as argument.
    The following variables MUST be set in globals:
        region_dict - a dictionary of water mass regions
        hconst - the assumed surface layer depth (in m)
        pref - the reference pressure for EOS
    """

    assert 'region_dict' in globals()
    assert 'hconst' in globals()
    assert 'pref' in globals()
    assert 'monthly_mean' in globals()

    p = pop_model.POPFile(pop_fname, hconst=hconst, pref=pref)
    p.initialize_gradient_operator()
    Nt = len(p.nc.variables['time'])
    
    result = dict()
    for rname in region_dict:
        result[rname] = 0.
    
    if monthly_mean:
    # take the monthly mean of T, S, fluxes first
        rho, Fheat, Fsalt, Fmix, Fcab = p.dens_forcing[:Nt]
        for rname, reg in region_dict.iteritems():
            A = reg.calculate_transformation_rate(
                                rho,Fheat, Fsalt, Fmix, Fcab)
            result[rname] += A
    else:
    # calculate transformation for each daily snapshot
        for n in range(Nt):
            print n
            rho, Fheat, Fsalt, Fmix, Fcab = p.dens_forcing[n]
            
            for rname, reg in region_dict.iteritems():
                A = reg.calculate_transformation_rate(
                                    rho,Fheat, Fsalt, Fmix, Fcab)
                result[rname] += A

        for rname in region_dict:
            result[rname] /= Nt

    return result

# This is too repetetive
# I guess we need a more generic framework
@interactive
def calc_Fd(pop_fname):
    """Worker function for Nakamura analysis.
    Takes a single file name as argument.
    The following variables MUST be set in globals:
        region_dict - a dictionary of water mass regions
        hconst - the assumed surface layer depth (in m)
        pref - the reference pressure for EOS
    """

    assert 'region_dict' in globals()
    assert 'hconst' in globals()
    assert 'pref' in globals()
    assert 'monthly_mean' in globals()

    p = pop_model.POPFile(pop_fname, hconst=hconst, pref=pref)
    p.initialize_gradient_operator()
    Nt = len(p.nc.variables['time'])
    
    result = dict()
    for rname in region_dict:
        result[rname] = 0.
    
    if monthly_mean:
        rho_mm = p.rho[:Nt]
        Nt = 1
        
    # calculate transformation for each daily snapshot
    for n in range(Nt):
        print n
        if monthly_mean:
            rho = rho_mm
        else:
            rho = p.rho[n]
        diss = p._ah * p._ahf * p.laplacian(rho)**2
        for rname, reg in region_dict.iteritems():
            dA, diss_sum = reg.sum_in_rholevs(rho, numpy.ones_like(diss), diss)
            dDiss_dA = diss_sum[1:] / dA[1:]
            dq_dA = numpy.ma.masked_invalid(numpy.diff(reg.rholevs) / dA[1:])
            A = numpy.cumsum(dA)[1:]

            Fd = numpy.ma.masked_invalid(dDiss_dA / dq_dA).filled(0.)
            result[rname] += numpy.array([A, Fd, dq_dA])

    for rname in region_dict:
        result[rname] /= Nt

    return result


def wmt_rho(aname, ddir, fprefix, years, pref=0, hconst=50.,
            task='calc_transformation_rates',
            monthly_mean=False, fsuffix=''):
    """Perform water mass analysis on a specific POP model run.
    aname - (string) the nickname of this specific analysis,
                used when saving data
        ddir - the directory where the run lives
        fprefix - the string that begins the file names
            (e.g. hybrid_v5_rel04_BC5_ne120_t12_pop62.pop.h.nday1)
        years - a list of years to analyze
        fsuffix - a trailing suffix (before .nc)
        pref - reference pressure for analysis
        hconst - depth of assumed surface layer
        monthly_mean - whether to calculate the
            transformation rate on the monthly mean instead of
            of the daily snapshots
    """

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

    #####################################
    ## Define Regions for Calculation ###
    #####################################

    fnames = []
    for year in years:
        for month in xrange(1,13):
            fname = '%s/%s.%04d-%02d%s.nc' % (ddir, fprefix, year, month, fsuffix)
            fnames.append(fname)
    for f in fnames:
        print f 

    # load a test file
    p = pop_model.POPFile(fnames[0], pref=pref)

    # define basins
    if pref==0:
        natl = transformation.WaterMassRegion(
                            basin_names=['Atlantic Ocean', 'GIN Seas', 'Labrador Sea'], latmin=15)
        natl.initialize_mask(p)
        natl.calculate_rholevs(rhomin=1022, rhomax=1028.5, nlevs=120, linear=True)

        npac = transformation.WaterMassRegion(
                            basin_names=['Pacific Ocean'], latmin=15)
        npac.initialize_mask(p)
        npac.calculate_rholevs(rhomin=1020, rhomax=1027, nlevs=120, linear=True)

        so = transformation.WaterMassRegion(
                            basin_names=['Southern Ocean'], latmax=30)
        so.initialize_mask(p)
        so.calculate_rholevs(rhomin=1022, rhomax=1030, nlevs=120, linear=True)

        globe = transformation.WaterMassRegion()
        globe.initialize_mask(p)
        globe.calculate_rholevs(rhomin=1018, rhomax=1030, nlevs=120, linear=True)
    elif pref==2000:  
        natl = transformation.WaterMassRegion(
                            basin_names=['Atlantic Ocean', 'GIN Seas', 'Labrador Sea'], latmin=15)
        natl.initialize_mask(p)
        natl.calculate_rholevs(rhomin=1028, rhomax=1037.8, nlevs=120, linear=True)
         
        npac = transformation.WaterMassRegion(
                            basin_names=['Pacific Ocean'], latmin=15)
        npac.initialize_mask(p)
        npac.calculate_rholevs(rhomin=1027, rhomax=1036.5, nlevs=120, linear=True)
         
        so = transformation.WaterMassRegion(
                            basin_names=['Southern Ocean'], latmax=30)
        so.initialize_mask(p)
        so.calculate_rholevs(rhomin=1030, rhomax=1037.8, nlevs=120, linear=True)

        globe = transformation.WaterMassRegion()
        globe.initialize_mask(p)
        globe.calculate_rholevs(rhomin=1026, rhomax=1037.8, nlevs=120, linear=True)
    else:
        raise ValueError('Invalid pressure level %g specified' % pref)


    region_dict = {'natl': natl, 'npac': npac, 'so': so, 'globe': globe}

    # push to engines
    dview.push(dict(hconst=hconst, pref=pref, monthly_mean=monthly_mean))
    dview.push(region_dict)
    dview.execute("region_dict = {'natl': natl, 'npac': npac, 'so': so, 'globe': globe}")
    # check
    dview.execute('a = region_dict.keys()[0]')
    a = dview.gather('a')
    for r in a.get():
        assert r=='npac'

    #######################
    ## Apply on Engines ###
    #######################

    if task=='calc_transformation_rates':
        mapfunc = calc_transformation_rates
        prefix = 'WMT'
    elif task=='calc_Fd':
        mapfunc = calc_Fd
        prefix = 'FD'

    res = lview.map(calc_transformation_rates, fnames)

    while not res.ready():
        print 'progress %3.2f%%' % 100*res.progress/float(len(res))
        time.sleep(60)

    assert res.successful()

    ###################
    ## Save Results ###
    ###################

    all_res = dict()
    for k in region_dict:
        all_res[k] = []
    for r in res :
        for k in all_res:
            all_res[k].append(r[k])
    for k in all_res:
        all_res[k] = numpy.array(all_res[k])
        numpy.savez('../data/%s_%s_sigma%1d_hconst%03d_%s.npz' % 
                        (prefix, aname, pref/1000, hconst, k),
                A=all_res[k], rholevs=region_dict[k].rholevs)
