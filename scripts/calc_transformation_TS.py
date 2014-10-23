from IPython.parallel import Client
import time

#####################
### This Analysis ###
#####################
hconst = 50. # assume a surface layer of const depth 50 m
anal_name = 'TS_hconst%03d' % hconst

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
ddir = '/glade/p/ncgd0001/hybrid_v5_rel04_BC5_ne120_t12_pop62/ocn-hist'
fprefix = 'hybrid_v5_rel04_BC5_ne120_t12_pop62.pop.h.nday1'
fnames = []
for year in xrange(47,86):
    for month in xrange(1,13):
        fnames.append(
         '%s/%s.%04d-%02d-01.nc' % (ddir, fprefix, year, month)
        )
              
# load a test file
p = pop_model.POPFile(fnames[0])

# define TS space
NT = 100
NS = 98
Tlevs = numpy.linspace(-4,34,NT)
Slevs = numpy.linspace(28,41,NS)

# define basins
natl = transformation.TSWaterMassRegion(Tlevs=Tlevs, Slevs=Slevs,
                    basin_names=['Atlantic Ocean', 'GIN Seas', 'Labrador Sea'], latmin=15)
natl.initialize_mask(p)                    
 
npac = transformation.TSWaterMassRegion(Tlevs=Tlevs, Slevs=Slevs,
                    basin_names=['Pacific Ocean'], latmin=15)
npac.initialize_mask(p)
 
so = transformation.TSWaterMassRegion(Tlevs=Tlevs, Slevs=Slevs,
                    basin_names=['Southern Ocean'], latmax=30)
so.initialize_mask(p)

globe = transformation.TSWaterMassRegion(Tlevs=Tlevs, Slevs=Slevs)
globe.initialize_mask(p)

region_dict = {'natl': natl, 'npac': npac, 'so': so, 'globe': globe}

# push to engines
dview.push(region_dict)
dview.execute(
    "region_dict = {'natl': natl, 'npac': npac, 'so': so, 'globe': globe}"
)
# check
dview.execute('a = region_dict.keys()[0]')
a = dview.gather('a')
for r in a.get():
    assert r=='npac'

#############################
## Engine Worker Function ###
#############################

def calc_transformation_rates(pop_fname):
    p = pop_model.POPFile(pop_fname, hconst=50.)
    p.initialize_gradient_operator()
    Nt = len(p.nc.variables['time'])
    
    result = dict()
    for rname in region_dict:
        result[rname] = 0.
    
    for n in range(Nt):
        print n
        
        # calculate density fluxes
        T, S, FT_forc, FT_mix, FS_forc, FS_mix = p.ts_forcing[n]
        FS_mix = numpy.ma.masked_greater(numpy.ma.masked_invalid(FS_mix),1.)
        FT_mix = numpy.ma.masked_greater(numpy.ma.masked_invalid(FT_mix),1.)
        
        for rname, reg in region_dict.iteritems():
            AT, AS = reg.calculate_transformation_rate(
                                T, S, [FT_forc, FT_mix], [FS_forc, FS_mix])
            result[rname] += numpy.array([AT, AS])

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
    numpy.savez('../data/A_%s_%s.npz' % (anal_name, k), A=all_res[k])
