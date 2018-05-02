from glob import glob
from os.path import isfile, split, join
import numpy as np
from datetime import datetime
from ._sate_param import _data_

try:
    from netCDF4 import Dataset
    READABLE = True
except ImportError:
    print('Note: Unable to find netCDF4 required for GOES.')
    READABLE = False

class GOES_NC:

    remap = True
    composite = {}
    channels = {'vis':('vis',), 'nir':None, 'swir':('swir',), 'ir2':('ir2',), 'wv':('wv',), 'ir':('ir',)}
    resolution = {'vis':1., 'ir':4., 'ir2':4., 'wv':4., 'swir':4.}
    crossdict = {'vis':1, 'ir':4, 'swir':2, 'wv':3, 'ir2':5}

    def __init__(self, fid):
        self.fid = fid
        self.ncfiles = dict()

    @staticmethod
    def search(filedir):
        filelist = glob(join(filedir, 'goes*.*.*.*.BAND_0*.nc'))
        options = []
        minor_idlist = []
        for entire_fname in filelist:
            fname = split(entire_fname)[1]
            time = datetime.strptime(fname[7:22], '%Y.%j.%H%M%S').strftime('%Y%m%d%H%M')
            sate = 'GOES-' + fname[4:6]
            if sate + time not in minor_idlist:
                fid = dict(time=time, area='', name=entire_fname, sate=sate, type='goes_nc', fdir=filedir)
                minor_idlist.append(sate + time)
                options.append(fid)
        return options

    @staticmethod
    def get_channel_filename(fname, channel):
        return fname[:28] + '%02d.nc' % (GOES_NC.crossdict[channel])

    def get_channel_info(self):
        filedir, fname = split(self.fid['name'])
        channel_flag = {}
        for channel in self.crossdict.keys():
            channel_flag[channel] = isfile(join(filedir, GOES_NC.get_channel_filename(fname, channel)))
        if not channel_flag['ir2']:
            self.crossdict['ir2'] = 6
            channel_flag['ir2'] = isfile(join(filedir, GOES_NC.get_channel_filename(fname, 'ir2')))
        return channel_flag

    def extract(self, channel, georange):
        filedir, fname = split(self.fid['name'])
        latmin, latmax, lonmin, lonmax = georange
        if channel in self.ncfiles.keys():
            ncfile = self.ncfiles[channel]
        else:
            ncfile = Dataset(join(filedir, GOES_NC.get_channel_filename(fname, channel)))
            self.ncfiles[channel] = ncfile
        data = ncfile.variables['data'][:][0]
        if channel in ('vis', 'swir', 'ir2'):
            return data / 32
        else:
            FK1, FK2, BC1, BC2 = _data_[self.fid['sate']]['planck_const'][channel]
            offset, factor = _data_[self.fid['sate']]['linear_const'][channel]
            return (FK2 / np.log1p(FK1 / ((data // 32 - factor) / offset)) - BC1) / BC2 - 273.15

    def geocoord(self, channel):
        ncfile = self.ncfiles[channel]
        lon, lat = np.ma.masked_greater(ncfile.variables['lon'][:], 360), ncfile.variables['lat'][:]
        lon[lon < 0] = 360 + lon[lon < 0]
        return lon, lat
