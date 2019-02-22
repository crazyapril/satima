from glob import glob
from os.path import isfile, split, join
import numpy as np
from ._sate_param import _data_
from datetime import datetime, timedelta
from numpy import loadtxt,savetxt,zeros
import h5py
READABLE = True

class COMS_HDF:

    remap = True
    crossdict = {'vis':'VIS', 'ir':'IR1', 'ir2':'IR2', 'wv':'IR3', 'swir':'IR4'}
    composite = {'RGB': ('vis', 'vis', 'ir'), 'IR-RGB': ('ir2', 'ir2', 'ir')}
    channels = {'vis':('vis',), 'nir':None, 'swir':('swir',), 'ir2':('ir2',), 'wv':('wv',), 'ir':('ir',)}
    resolution = {'vis':1., 'ir':4., 'ir2':4., 'wv':4., 'swir':4.}

    def __init__(self, fid):
        self.fid = fid
        self.hrit_dict = {}
        
    @staticmethod      
    def search(filedir):
        filelist = glob(join(filedir, 'coms_mi_le1b_zzz_*.he5'))
        options = []
        minor_idlist = []
        for entire_fname in filelist:
         fname = split(entire_fname)[1]
         time = fname[20:32]
         if time not in minor_idlist:
            fid = {}
            fid['time'] = time
            fid['name'] = entire_fname
            fid['sate'] = 'COMS-1'
            fid['area'] = {'cf':'FD', 'cn':'ENH'}[fname[17:19]]
            fid['type'] = 'coms_hdf'
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
            channel_flag[channel] = isfile(join(filedir, COMS_HDF.get_channel_filename(fname, channel)))
        return channel_flag
  
    def extract(self, channel, georange):
        filedir, fname = split(self.fid['name'])
        latmin, latmax, lonmin, lonmax = georange 
        filename = join(filedir, COMS_HDF.get_channel_filename(fname, channel))
        count = _HDFReader(filename,channel)

        TBB=COMS1_Calibration(count,channel)
        print(TBB.shape)
        if channel  == 'vis' :  
           return TBB/100
        else :     
           return TBB-273.15
        
    def geocoord(self, channel):
        filedir, fname = split(self.fid['name'])
        filename = join(filedir, COMS_HDF.get_channel_filename(fname, channel))    
        lon,lat=_GetLonLat(filename,channel)
        print(np.amax(lon),np.amin(lon),np.amax(lat),np.amin(lat))
        lon[lon < 0] = 360 + lon[lon < 0]
        return lon, lat   
           

        
def _GetLonLat(fname,channel):
    hdf5_obj=h5py.File(fname, 'r')
    if channel == 'ir' :
     NL,NC = hdf5_obj['HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/IR1 band Image Pixel Values'][:].shape
    elif channel == 'ir2' :
     NL,NC = hdf5_obj['HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/IR2 band Image Pixel Values'][:].shape
    elif channel == 'wv' :
     NL,NC = hdf5_obj['HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/WV band Image Pixel Values'][:].shape
    elif channel == 'swir' :  
     NL,NC = hdf5_obj['HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/SWIR band Image Pixel Values'][:].shape
    elif channel == 'vis' :  
     NL,NC = hdf5_obj['HDFEOS/GRIDS/Visible Image Pixel Value/Data Fields/Visible Image Pixel Values'][:].shape
    if channel == 'vis' :
     proj=list(hdf5_obj['HDFEOS/POINTS/Navigation Result/Data/NAVI MI/'].attrs['Visible'])
     CFAC = proj[1]
     LFAC = proj[3]*(-1)
     COFF = proj[0]
     LOFF = proj[2]
    else :
     proj=list(hdf5_obj['HDFEOS/POINTS/Navigation Result/Data/NAVI MI/'].attrs['IR'])
     CFAC = proj[1]
     LFAC = proj[3]*(-1)
     COFF = proj[0]
     LOFF = proj[2]     

    SubLon = list(hdf5_obj['HDFEOS/POINTS/Satellite Status/Data/Satellite Definition/'].attrs['Longitude'])[0]


    EarthPolarRadius=list(hdf5_obj['HDFEOS/POINTS/Geometric Processing/Data/Earth Model/'].attrs['Polar Radius of Earth Ellipsoid'])[0]/1000
    EarthEquatorRadius=list(hdf5_obj['HDFEOS/POINTS/Geometric Processing/Data/Earth Model/'].attrs['Equitorial Radius of Earth Ellipsoid'])[0]/1000
    EarthConst1=(EarthEquatorRadius**2-EarthPolarRadius**2)/EarthEquatorRadius**2
    EarthConst2=EarthPolarRadius**2/EarthEquatorRadius**2
    EarthConst3=EarthEquatorRadius**2/EarthPolarRadius**2
    Distance=list(hdf5_obj['HDFEOS/POINTS/Satellite Status/Data/Satellite_Orbit/'].attrs['Altitude at Scene Center Time'])[0]
    EarthConstStd=Distance**2-EarthEquatorRadius**2
    
    DEGTORAD = np.pi / 180.
    RADTODEG = 180. / np.pi
    SCLUNIT = np.power(2., -16)
    DIS = Distance
    CON = EarthConst3
    #Calculation
    lines = np.arange(1,NL+1)
    columns = np.arange(1,NC+1)
    print(lines,columns)
    xx, yy = np.meshgrid(columns, lines)
    del lines, columns
    x = DEGTORAD * (xx - COFF) / (SCLUNIT * CFAC)
    y = DEGTORAD * (yy - LOFF) / (SCLUNIT * LFAC)
    del xx, yy
    Sd = np.sqrt(np.square(DIS * np.cos(x) * np.cos(y)) - (np.square(np.cos(y)) + CON * np.square(np.sin(y))) * EarthConstStd)
    Sn = (DIS * np.cos(x) * np.cos(y) - Sd) / (np.square(np.cos(y)) + CON * np.square(np.sin(y)))
    S1 = DIS - Sn * np.cos(x) * np.cos(y)
    S2 = Sn * np.sin(x) * np.cos(y)
    S3 = -Sn * np.sin(y)
    Sxy = np.sqrt(np.square(S1) + np.square(S2))
    del x, y, Sd, Sn
    lons = RADTODEG * np.arctan2(S2, S1) + SubLon
    lats = np.ma.masked_outside(RADTODEG * np.arctan(CON * S3 / Sxy), -90., 90.)
    return np.ma.masked_invalid(lons), np.ma.masked_invalid(lats)  
        

def read_cds_time(buf):
    days = read_uint2(buf[:2])
    msecs = read_uint4(buf[2:6])
    return datetime(1958, 1, 1) + timedelta(days=days, milliseconds=msecs)    

def _HDFReader(fname,channel):
    hdf5_obj=h5py.File(fname, 'r')
    if channel == 'ir' :
     raw = hdf5_obj['HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/IR1 band Image Pixel Values'][:]
    elif channel == 'ir2' :
     raw = hdf5_obj['HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/IR2 band Image Pixel Values'][:]
    elif channel == 'wv' :
     raw = hdf5_obj['HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/WV band Image Pixel Values'][:]
    elif channel == 'swir' :  
     raw = hdf5_obj['HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/SWIR band Image Pixel Values'][:]
    elif channel == 'vis' :  
     raw = hdf5_obj['HDFEOS/GRIDS/Visible Image Pixel Value/Data Fields/Visible Image Pixel Values'][:]  
    
    return raw
        
def COMS1_Calibration(count,channel):
    DN=loadtxt('./satellite/COMS-1.txt')[:,0]
    if channel == 'ir' :
     CTB=loadtxt('./satellite/COMS-1.txt')[:,1]
    elif channel == 'ir2' :
     CTB=loadtxt('./satellite/COMS-1.txt')[:,2]
    elif channel == 'wv' :
     CTB=loadtxt('./satellite/COMS-1.txt')[:,3]
    elif channel == 'swir' :  
     CTB=loadtxt('./satellite/COMS-1.txt')[:,4]
    elif channel == 'vis' :  
     CTB=loadtxt('./satellite/COMS-1.txt')[:,5]     
    dictionary = dict(zip(DN,CTB))
    TBB=np.vectorize(dictionary.__getitem__)(count)   
    return TBB
    
