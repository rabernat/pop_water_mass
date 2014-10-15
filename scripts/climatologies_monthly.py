import os
from subprocess import call, Popen
import time

odir = '/glade/p/work/rpa/pop/hybrid_v5_rel04_BC5_ne120_t12_pop62_DERIVED'
ddir = '/glade/p/ncgd0001/hybrid_v5_rel04_BC5_ne120_t12_pop62/ocn-hist/'
prefix = 'hybrid_v5_rel04_BC5_ne120_t12_pop62.pop.h.nday1'
years = range(46,87)
seasons = {'DJF': (12,1,2), 'MAM': (3,4,5), 'JJA': (6,7,8), 'SON': (9,10,11)}

np = 1 # multithreading doesn't work, must not be compiled into nco
	
procs = []
for m in range(1,13):
    filelist = []
    for y in years:
        fname = '%s/%s.%04d-%02d-01.nc' % (ddir, prefix, y, m)
        assert os.path.exists(fname)
        filelist.append(fname)
    output_fname = '%s/climatology.%02d.nc' % (odir,m)
    print 'Launching', output_fname
    procs.append(Popen(['ncra', '-O', '-t', str(np), '-o', output_fname] + filelist))


finished = False
while not finished:
    stat = [p.poll() for p in procs]
    print 'STAT: ', stat
    working = [s is None for s in stat]
    if not any(working):
        finished = True
    time.sleep(60)
