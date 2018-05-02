from glob import glob
from os.path import isfile, split, join
import numpy as np
from ._sate_param import _data_

try:
    import h5py
    READABLE = False
except ImportError:
    print('Note: Unable to find h5py required for SNPP VIIRS.')
    READABLE = False

class SNPP_HDF:

    def __init__(fid):
        self.fid = fid

    @staticmethod
    def search(filedir):
        filelist = glob(join(filedir, 'G*-SV*_npp_d*_t*_e*_b*_c*_noaa_ops.h5'))
        options = []
        minor_idlist = []
        for entire_fname in filelist:
            fname = split(entire_fname)[1]
            if fname[16:34] not in minor_idlist:
                fid = {}
                fid['time'] = fname[17:25] + fname[27:31]
                fid['area'] = ''
                fid['name'] = entire_fname
                fid['sate'] = 'S-NPP'
                fid['type'] = 'snpp_hdf'
                fid['fdir'] = filedir
                minor_idlist.append(fname[16:34])
                options.append(fid)
        return options
