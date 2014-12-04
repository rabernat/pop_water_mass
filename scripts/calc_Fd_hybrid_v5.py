from parallel_analysis import wmt_rho

hconst = 50. # assume a surface layer of const depth 50 m
pref = 0

name = 'hybrid_v5'
ddir = '/glade/p/ncgd0001/hybrid_v5_rel04_BC5_ne120_t12_pop62/ocn-hist'
fprefix = 'hybrid_v5_rel04_BC5_ne120_t12_pop62.pop.h.nday1'
fsuffix = '-01'
years = range(47,86)

wmt_rho(name, ddir, fprefix, years, fsuffix=fsuffix,
        task='calc_Fd',
        pref=pref, hconst=hconst)
