from glob import glob
from os.path import isfile, split, join
import numpy as np
from datetime import datetime
from ._sate_param import _data_
from scipy.interpolate import CubicSpline
from numpy import loadtxt,savetxt,zeros

try:
    from netCDF4 import Dataset
    READABLE = True
except ImportError:
    print('Note: Unable to find netCDF4 required for GOES.')
    READABLE = False

class MFG:

    remap = True
    composite = {}
    channels = {'vis':('vis',), 'nir':None, 'swir':None, 'ir2':None, 'wv':('wv',), 'ir':('ir',)}
    resolution = {'vis':5., 'ir':5.,'wv':5.}
    crossdict = {'vis':5, 'ir':5, 'wv':5}

    def __init__(self, fid):
        self.fid = fid
        self.ncfiles = dict()

    @staticmethod
    def search(filedir):
        filelist = glob(join(filedir, 'W_XX-EUMETSAT-Darmstadt,VIS+IR+IMAGERY,MET*+MVIRI_C_EUMS_*.nc'))
        options = []
        minor_idlist = []
        for entire_fname in filelist:
            fname = split(entire_fname)[1]
            print(fname[57:71])
            time = datetime.strptime(fname[57:71], '%Y%m%d%H%M%S').strftime('%Y%m%d%H%M')
            sate = 'METEOSAT-' + fname[42:43]
            if sate + time not in minor_idlist:
                fid = dict(time=time, area='', name=entire_fname, sate=sate, type='MFG', fdir=filedir)
                minor_idlist.append(sate + time)
                options.append(fid)
        return options

    @staticmethod
    def get_channel_filename(fname, channel):
        return fname[:]

    def get_channel_info(self):
        filedir, fname = split(self.fid['name'])
        channel_flag = {}
        for channel in self.crossdict.keys():
            channel_flag[channel] = isfile(join(filedir, MFG.get_channel_filename(fname, channel)))
        return channel_flag

    def extract(self, channel, georange):
        filedir, fname = split(self.fid['name'])
        latmin, latmax, lonmin, lonmax = georange
        print(channel,self.fid['sate'])
        if channel in self.ncfiles.keys():
            ncfile = self.ncfiles[channel]
        else:
            ncfile = Dataset(join(filedir, MFG.get_channel_filename(fname, channel)))
            self.ncfiles[channel] = ncfile

        if channel in ('ir'):
            data = ncfile.variables['ch2']
            data.set_auto_scale(False)
            count = data[:]
            scale_factor=ncfile.variables['ch2'].scale_factor
            add_offset=ncfile.variables['ch2'].add_offset
            radiance=scale_factor*(count - add_offset)
            if self.fid['sate'] == 'METEOSAT-5' :
             TB=loadtxt('./satellite/METEO-5IR1.txt')[:,0]
             R=loadtxt('./satellite/METEO-5IR1.txt')[:,1]
            elif self.fid['sate'] == 'METEOSAT-6' :
             TB=loadtxt('./satellite/METEO-6IR1.txt')[:,0]
             R=loadtxt('./satellite/METEO-6IR1.txt')[:,1]             
            elif self.fid['sate'] == 'METEOSAT-7' :
             TB=loadtxt('./satellite/METEO-7IR2.txt')[:,0]
             R=loadtxt('./satellite/METEO-7IR2.txt')[:,1]  
            cs = CubicSpline(R, TB)
            return cs(radiance)-273.15
        elif channel in ('wv'):
            data = ncfile.variables['ch3']
            data.set_auto_scale(False)
            count = data[:]
            scale_factor=ncfile.variables['ch3'].scale_factor
            add_offset=ncfile.variables['ch3'].add_offset
            radiance=scale_factor*(count - add_offset)
            if self.fid['sate'] == 'METEOSAT-5' :
             TB=loadtxt('./satellite/METEO-5WV2.txt')[:,0]
             R=loadtxt('./satellite/METEO-5WV2.txt')[:,1]
            elif self.fid['sate'] == 'METEOSAT-6' :
             TB=loadtxt('./satellite/METEO-6WV2.txt')[:,0]
             R=loadtxt('./satellite/METEO-6WV2.txt')[:,1]             
            elif self.fid['sate'] == 'METEOSAT-7' :
             TB=loadtxt('./satellite/METEO-7WV2.txt')[:,0]
             R=loadtxt('./satellite/METEO-7WV2.txt')[:,1]  
            cs = CubicSpline(R, TB)
            return cs(radiance)-273.15         
        elif channel in ('vis'):
            data = ncfile.variables['ch1']
            data.set_auto_scale(False)
            count = data[:]
            scale_factor=ncfile.variables['ch1'].scale_factor
            add_offset=ncfile.variables['ch1'].add_offset
            radiance=scale_factor*(count - add_offset)
            return radiance / 255	#16384


    def geocoord(self, channel):
        ncfile = self.ncfiles[channel]
        if channel == 'vis' :
         lon, lat = np.ma.masked_greater(ncfile.variables['lon_vis'][:], 360), ncfile.variables['lat_vis'][:]
        else :
         lon, lat = np.ma.masked_greater(ncfile.variables['lon'][:], 360), ncfile.variables['lat'][:]
        lon[lon < 0] = 360 + lon[lon < 0]
        return lon, lat
