from parallel_analysis import wmt_rho

hconst = 50. # assume a surface layer of const depth 50 m
pref = 2000
anal_name = 'hconst%03d_new' % hconst
# where to find the data
ddir = '/glade/p/ncgd0001/hybrid_v5_rel04_BC5_ne120_t12_pop62/ocn-hist'
fprefix = 'hybrid_v5_rel04_BC5_ne120_t12_pop62.pop.h.nday1'
fsuffix = '-01'
years = range(47,86)

wmt_rho(anal_name, ddir, fprefix, years, fsuffix=fsuffix,
        pref=pref, hconst=hconst)
