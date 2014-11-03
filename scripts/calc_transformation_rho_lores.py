from IPython.parallel import Client
import time
import os

#####################
### This Analysis ###
#####################
hconst = 50. # assume a surface layer of const depth 50 m
anal_name = 'lores'

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

# where to find the data
ddir = '/glade/scratch/enewsom/HR_ANALYSIS/DATA/LRfiles/'
fprefix = 'LRC01.pop.h1'
fnames = []
years = range(148,168)
years.remove(159) # for some reason this year is screwed up
for year in years:
    for month in xrange(1,13):
        fname = '%s/%s.%04d-%02d.nc' % (ddir, fprefix, year, month)
        fnames.append(fname)
        if not os.path.exists(fname):
            print 'does not exist: ', fname
print fnames 

# referene pressure
pref = 2000.              


# load a test file
p = pop_model.POPFile(fnames[0], pref=pref)
# this doesn't actually get used
rho = p.mask

# define basins
natl = transformation.WaterMassRegion(
                    basin_names=['Atlantic Ocean', 'GIN Seas', 'Labrador Sea'], latmin=15)
natl.initialize_mask(p)
natl.calculate_rholevs(rho, rhomin=1028, rhomax=1037.8, nlevs=120, linear=True)
 
npac = transformation.WaterMassRegion(
                    basin_names=['Pacific Ocean'], latmin=15)
npac.initialize_mask(p)
npac.calculate_rholevs(rho, rhomin=1027, rhomax=1036.5, nlevs=120, linear=True)
 
so = transformation.WaterMassRegion(
                    basin_names=['Southern Ocean'], latmax=30)
so.initialize_mask(p)
so.calculate_rholevs(rho, rhomin=1030, rhomax=1037.8, nlevs=120, linear=True)

globe = transformation.WaterMassRegion()
globe.initialize_mask(p)
globe.calculate_rholevs(rho, rhomin=1026, rhomax=1037.8, nlevs=120, linear=True)

region_dict = {'natl': natl, 'npac': npac, 'so': so, 'globe': globe}

# push to engines
dview.push(region_dict)
dview.execute("region_dict = {'natl': natl, 'npac': npac, 'so': so, 'globe': globe}")
# check
dview.execute('a = region_dict.keys()[0]')
a = dview.gather('a')
for r in a.get():
    assert r=='npac'

#############################
## Engine Worker Function ###
#############################

def calc_transformation_rates(pop_fname):
    p = pop_model.POPFile(pop_fname, hconst=50., pref=2000.)
    p.initialize_gradient_operator()
    Nt = len(p.nc.variables['time'])
    
    result = dict()
    for rname in region_dict:
        result[rname] = 0.
    
    for n in range(Nt):
        print n
        
        # calculate density fluxes
        rho, Fheat, Fsalt, Fmix = p.dens_forcing[n]
        Fmix = numpy.ma.masked_greater(numpy.ma.masked_invalid(Fmix),1.)
        
        for rname, reg in region_dict.iteritems():
            A = reg.calculate_transformation_rate(
                                rho,Fheat, Fsalt, Fmix)
            result[rname] += A

    for rname in region_dict:
        result[rname] /= Nt
    
    return result

#######################
## Apply on Engines ###
#######################

res = lview.map(calc_transformation_rates, fnames)

while not res.ready():
    print 'progress', res.progress / float(len(res))
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
    numpy.savez('../data/A_rho_%s_%s.npz' % (anal_name, k), A=all_res[k], rholevs=region_dict[k].rholevs)
