from parallel_analysis import wmt_rho

# general params
hconst = 50
pref = 2000

# hires version
name = 'hires_CO2'
ddir = '/glade/scratch/enewsom/HR_ANALYSIS/DATA/HRfiles/'
fprefix = 'HRC03.branch.CO2ramp.pop.h1'
fsuffix = ''
years = range(147,168)

wmt_rho(name, ddir, fprefix, years, fsuffix=fsuffix,
        pref=pref, hconst=hconst)

# lores version
name = 'lores_CO2'
ddir = '/glade/scratch/enewsom/HR_ANALYSIS/DATA/LRfiles/'
fprefix = 'LRC01.CO2ramp.pop.h1'
years = range(147,168)

wmt_rho(name, ddir, fprefix, years, fsuffix=fsuffix,
        pref=pref, hconst=hconst)
