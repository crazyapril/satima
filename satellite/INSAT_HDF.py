from glob import glob
from os.path import isfile, split, join
import numpy as np
from ._sate_param import _data_
import struct
from datetime import datetime, timedelta
from numpy import loadtxt,savetxt,zeros
import h5py
READABLE = True

class INSAT_HDF:

    remap = True
    crossdict = {'vis':'VIS', 'ir':'IR1', 'ir2':'IR2', 'wv':'IR3', 'swir':'IR4', 'mir':'MIR'}
    composite = {'RGB': ('vis', 'vis', 'ir'), 'IR-RGB': ('ir2', 'ir2', 'ir')}
    channels = {'vis':('vis',), 'nir':None, 'swir':('swir',), 'wv':('wv',), 'ir':('ir','ir2'), 'mir':('mir',)}
    resolution = {'vis':1., 'ir':4., 'ir2':4., 'wv':8., 'swir':4., 'mir':4.}

    def __init__(self, fid):
        self.fid = fid
        self.hrit_dict = {}
        
    @staticmethod      
    def search(filedir):
        filelist = glob(join(filedir, '3*IMG_*_L1B_STD.h5'))
        options = []
        minor_idlist = []
        for entire_fname in filelist:
         fname = split(entire_fname)[1]
         month = {'JAN':'01', 'FEB':'02','MAR':'03','APR':'04','MAY':'05','JUN':'06','JUL':'07','AUG':'08','SEP':'09','OCT':'10','NOV':'11','DEC':'12'}[fname[8:11]]
         time = fname[11:15] + month + fname[6:8] + fname[16:20]
         if time not in minor_idlist:
            fid = {}
            fid['time'] = time
            fid['name'] = entire_fname
            if fname[1:2]=='D' : 
             fid['sate'] = 'INSAT-3D'
            elif fname[1:2]=='R' :  
             fid['sate'] = 'INSAT-3DR'
            fid['area'] = ''
            fid['type'] = 'INSAT_HDF'
            fid['fdir'] = filedir         
            minor_idlist.append(time)
            options.append(fid)
        return options        
        
    def get_channel_filename(fname, channel):  
        return fname[:]

    def get_channel_info(self):
        filedir, fname = split(self.fid['name'])
        channel_flag = {}
        for channel in self.crossdict.keys():
            channel_flag[channel] = isfile(join(filedir, INSAT_HDF.get_channel_filename(fname, channel)))
        return channel_flag
  
    def extract(self, channel, georange):
        filedir, fname = split(self.fid['name'])
        latmin, latmax, lonmin, lonmax = georange 
        filename = join(filedir, INSAT_HDF.get_channel_filename(fname, channel))
        count = _HDFReader(filename,channel)
        print('COUNT.SHAPE',count.shape,channel)
        TBB=INSAT_Calibration(count,channel,fname)
        print(TBB.shape,np.amax(TBB))
        if channel  == 'vis' :  
           return TBB/100
        else :     
           return TBB-273.15
        
    def geocoord(self, channel):
        filedir, fname = split(self.fid['name'])
        filename = join(filedir, INSAT_HDF.get_channel_filename(fname, channel))    
        lon,lat=_GetLonLat(filename,channel)
        print(np.amax(lon),np.amin(lon),np.amax(lat),np.amin(lat))
        lon[lon < 0] = 360 + lon[lon < 0]
        return lon, lat   
        
def _GetLonLat(fname,channel):
    hdf5_obj=h5py.File(fname, 'r')
    if channel == 'vis' :
     lats = hdf5_obj['Latitude_VIS'][:]/1000
     lons = hdf5_obj['Longitude_VIS'][:]/1000
    elif channel == 'wv' :
     lats = hdf5_obj['Latitude_WV'][:]/100
     lons = hdf5_obj['Longitude_WV'][:]/100
    else :
     lats = hdf5_obj['Latitude'][:]/100
     lons = hdf5_obj['Longitude'][:]/100
    return lons, lats

        

def read_cds_time(buf):
    days = read_uint2(buf[:2])
    msecs = read_uint4(buf[2:6])
    return datetime(1958, 1, 1) + timedelta(days=days, milliseconds=msecs)    

def _HDFReader(fname,channel):
    hdf5_obj=h5py.File(fname, 'r')
    if channel == 'ir' :
     raw = hdf5_obj['IMG_TIR1'][0,:,:]
    elif channel == 'ir2' :
     raw = hdf5_obj['IMG_TIR2'][0,:,:]
    elif channel == 'wv' :
     raw = hdf5_obj['IMG_WV'][0,:,:]
    elif channel == 'swir' :  
     raw = hdf5_obj['IMG_SWIR'][0,:,:]
    elif channel == 'mir' :  
     raw = hdf5_obj['MIR'][0,:,:]   
    elif channel == 'vis' :  
     raw = hdf5_obj['IMG_VIS'][0,:,:]
    
    return raw
        
def INSAT_Calibration(count,channel,fname):
    hdf5_obj=h5py.File(fname, 'r')
    DN=np.arange(1024)
    if channel == 'ir' :
     CTB = hdf5_obj['IMG_TIR1_TEMP'][:]
    elif channel == 'ir2' :
     CTB = hdf5_obj['IMG_TIR2_TEMP'][:]
    elif channel == 'wv' :
     CTB = hdf5_obj['IMG_WV_TEMP'][:]
    elif channel == 'swir' :  
     CTB = hdf5_obj['IMG_SWIR_RADIANCE'][:]
    elif channel == 'mir' :  
     CTB = hdf5_obj['MIR_TEMP'][:]   
    elif channel == 'vis' :  
     CTB = hdf5_obj['IMG_VIS_ALBEDO'][:]   
    dictionary = dict(zip(DN,CTB))
    TBB=np.vectorize(dictionary.__getitem__)(count)   
    return TBB 
         