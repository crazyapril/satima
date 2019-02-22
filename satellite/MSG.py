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

class MSG:
    remap = True
    composite = {}
    channels = {'vis':('band01','band02',), 'nir':('band03',), 'swir':('band04','band07','band08',), 'ir2':('band11',), 'wv':('band05','band06',), 'ir':('band09','band10',)}
    resolution = {'band01':3., 'band02':3, 'band03':3, 'band04':3., 'band05':3., 'band06':3.,
                  'band07':3., 'band08':3., 'band09':3., 'band10':3., 'band11':3.}
    crossdict = {'band01':1, 'band02':2, 'band03':3, 'band04':4, 'band05':5, 'band06':6, 'band07':7, 'band08':8, 'band09':9, 'band10':10, 'band11':11}
    def __init__(self, fid):
        self.fid = fid
        self.ncfiles = dict()


    @staticmethod
    def search(filedir):
        filelist = glob(join(filedir, 'W_XX-EUMETSAT-Darmstadt,VIS+IR+IMAGERY,MSG*+SEVIRI_C_EUMG_*.nc'))
        options = []
        minor_idlist = []
        for entire_fname in filelist:                  
            fname = split(entire_fname)[1]
            print(fname[58:72])
            time = datetime.strptime(fname[58:72], '%Y%m%d%H%M%S').strftime('%Y%m%d%H%M')
            sate = 'METEOSAT-' + str((int(fname[42:43])+7))
            if sate + time not in minor_idlist:
                fid = dict(time=time, area='', name=entire_fname, sate=sate, type='MSG', fdir=filedir)
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
            channel_flag[channel] = isfile(join(filedir, MSG.get_channel_filename(fname, channel)))
        return channel_flag

    def extract(self, channel, georange):
        filedir, fname = split(self.fid['name'])
        latmin, latmax, lonmin, lonmax = georange
        print(channel,self.fid['sate'])
        if channel in self.ncfiles.keys():
            ncfile = self.ncfiles[channel]
        else:
            ncfile = Dataset(join(filedir, MSG.get_channel_filename(fname, channel)))
            self.ncfiles[channel] = ncfile            
            
        channelindex = int(channel[4:6])

        c=299792458        
        h=6.62606957*(10**-34.0)
        k=1.3806488*(10**-23.0)
        '''
        C1=1.191044E-5
        C2=1.438769'''
        C1=2*h*(c**2)*10**11
        C2=h*c/k*100
        if channel == 'band01' or channel == 'band02' or channel == 'band03' :
         data = ncfile.variables['ch%s'%(channelindex)]
         data.set_auto_scale(False)
         count = data[:]
         return count / 900
        else : 
         data = ncfile.variables['ch%s'%(channelindex)]
         data.set_auto_scale(False)
         count = data[:]
         scale_factor=ncfile.variables['ch%s'%(channelindex)].scale_factor
         add_offset=ncfile.variables['ch%s'%(channelindex)].add_offset
         radiance=(scale_factor*count) + add_offset
         alpha = _data_[self.fid['sate']]['alpha'][channel]
         beta = _data_[self.fid['sate']]['beta'][channel]
         wnc =  _data_[self.fid['sate']]['wnc'][channel]
         print(alpha,beta,wnc,C1,C2)
         TBB=((C2*wnc)/np.log(1.+(C1*wnc**3)/radiance)-beta)/alpha-273.15  
  
         return TBB

    def geocoord(self, channel):
        ncfile = self.ncfiles[channel]
        lon, lat = np.ma.masked_greater(ncfile.variables['lon'][:], 360), ncfile.variables['lat'][:]        
        return lon, lat
         
         
    
    
    
    
    