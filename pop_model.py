import numpy as np
import netCDF4

class POPFile(netCDF4.Dataset):
    
    def __init__(self, fname):
        """Wrapper for POP model netCDF files"""
        netCDF4.Dataset.__init__(self, fname)
    