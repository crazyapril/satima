from glob import glob
from os.path import isfile, split, join
import numpy as np
from ._sate_param import _data_

READABLE = True

class MTSAT_GEOSS:

    remap = False
    crossdict = {'vis':'VIS', 'ir':'IR1', 'ir2':'IR2', 'wv':'IR3', 'swir':'IR4'}
    composite = {'RGB': ('vis', 'vis', 'ir'), 'IR-RGB': ('ir2', 'ir2', 'ir')}
    channels = {'vis':('vis',), 'nir':None, 'swir':('swir',), 'ir2':('ir2',), 'wv':('wv',), 'ir':('ir',)}
    resolution = {'vis':1., 'ir':4., 'ir2':4., 'wv':4., 'swir':4.}

    def __init__(self, fid):
        self.fid = fid
        
    @staticmethod
    def search(filedir):
        filelist = glob(join(filedir, 'IMG_DK0*.geoss'))
        options = []
        minor_idlist = []
        for entire_fname in filelist:
            fname = split(entire_fname)[1]
            time = fname[12:24]
            if time not in minor_idlist:
                fid = {}
                fid['time'] = time
                fid['area'] = {'1':'Full', '2':'NH', '3':'SH'}[fname[7]]
                fid['name'] = entire_fname
                if int(time) < 201007010300: fid['sate'] = 'MTSAT-1R'
                else: fid['sate'] = 'MTSAT-2'
                fid['type'] = 'mtsat_geoss'
                fid['fdir'] = filedir
                minor_idlist.append(time)
                options.append(fid)
        return options

    @staticmethod
    def get_channel_filename(fname, channel):
        return fname[:8] + MTSAT_GEOSS.crossdict[channel] + fname[11:]

    def get_channel_info(self):
        filedir, fname = split(self.fid['name'])
        channel_flag = {}
        for channel in self.crossdict.keys():
            channel_flag[channel] = isfile(join(filedir, MTSAT_GEOSS.get_channel_filename(fname, channel)))
        return channel_flag

    def extract(self, channel, georange):
        latmin, latmax, lonmin, lonmax = georange
        visflag = channel == 'vis'
        filedir, fname = split(self.fid['name'])
        channel_filename = join(filedir, MTSAT_GEOSS.get_channel_filename(fname, channel))
        pixels = 12000 if visflag else 3000
        pixel_step = 0.01 if visflag else 0.04
        rawdata = np.fromfile(channel_filename, dtype=np.dtype('>H')).reshape(pixels, pixels)
        satlatmax, satlatmin, satlonmax, satlonmin = _data_[self.fid['sate']]['georange']
        pixel_lonmin = (lonmin - satlonmin) // pixel_step
        pixel_lonmax = (lonmax - satlonmin) // pixel_step
        pixel_latmin = (satlatmax - latmin) // pixel_step
        pixel_latmax = (satlatmax - latmax) // pixel_step
        rawdata = rawdata[pixel_latmax:pixel_latmin, pixel_lonmin:pixel_lonmax]
        if channel in ['vis', 'ir2', 'swir']:
            rawdata = rawdata / 1024
        else:
            slope, intercept = _data_[self.fid['sate']]['linear_const'][channel]
            FK1, FK2, BC1, BC2 = _data_[self.fid['sate']]['planck_const'][channel]
            rawdata = (FK2 / np.log1p(FK1 / (slope * rawdata + intercept) + 1) - BC1) / BC2 - 273.15
        return rawdata
