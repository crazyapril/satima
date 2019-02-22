from glob import glob
from os.path import isfile, split, join
import numpy as np
from datetime import datetime
import h5py
from ._sate_param import _FY4_data_for_seg_calc as _data_

READABLE = True



class FY4_HDF:
    '''channels = {'vis':('band01', 'band02', 'band03'), 'nir':('band04', 'band05', 'band06'),
                'swir':('band07', 'band11', 'band12'), 'wv':('band08', 'band09', 'band10'),
                'ir':('band13', 'band14', 'band15'), 'ir2':('band16',)}'''
    remap = True
    
    composite = {}
    '''
    channels = {'vis':('vis',), 'nir':None, 'swir':('swir',), 'ir2':('ir2',), 'wv':('wv',), 'ir':('ir',)}
    resolution = {'vis':1., 'ir':4., 'ir2':4., 'wv':4., 'swir':4.}
    crossdict = {'vis':1, 'ir':13, 'swir':2, 'wv':3, 'ir2':5}'''
    channels = {'vis':('band01','band02','band03',), 'nir':('band04','band05','band06',), 'swir':('band07','band8','band11',), 'ir2':('band14',), 'wv':('band09','band10',), 'ir':('band12','band13',)}
    resolution = {'band01':1., 'band02':0.5, 'band03':1, 'band04':2., 'band05':2., 'band06':2.,
                  'band07':2., 'band08':4., 'band09':4., 'band10':4., 'band11':4., 'band12':4.,
                  'band13':4., 'band14':4.}
    crossdict = {'band01':1, 'band02':2, 'band03':3, 'band04':4, 'band05':5, 'band06':6, 'band07':7, 'band08':8, 'band09':9, 'band10':10, 'band11':11, 'band12':12, 'band13':13, 'band14':14}
    def __init__(self, fid):
        self.fid = fid
        self.ncfiles = dict()
        self.georange = ''
        self.lineN = 0
        self.lineS = 0
        self.lineW = 0
        self.lineE = 0

#FY4A-_AGRI--_N_DISK_1047E_L1-_FDI-_MULT_NOM_20180321050000_20180321051459_0500M_V0001.HDF
#FY4A-_AGRI--_N_REGC_1047E_L1-_FDI-_MULT_NOM_20180321053000_20180321053416_0500M_V0001.HDF
        
    @staticmethod
    def search(filedir):
        filelist = glob(join(filedir, 'FY4*-_AGRI--_N_*_*_L1-_FDI-_MULT_NOM_*_*_*_V0001.HDF'))
        options = []
        minor_idlist = []
        for entire_fname in filelist:
            fname = split(entire_fname)[1]            
            if fname[15:19] == 'DISK' :
             time = datetime.strptime(fname[44:58], '%Y%m%d%H%M%S').strftime('%Y%m%d%H%M')             
             sate = 'FY-' + fname[2:4]
             area = 'FULL DISK '+ fname[74:78] + ' M'
             fid = dict(time=time, area=area, name=entire_fname, sate=sate, type='FY4_HDF', fdir=filedir)
             minor_idlist.append(sate + time)
             options.append(fid)
            elif fname[15:19] == 'REGC' :
             time = datetime.strptime(fname[44:58], '%Y%m%d%H%M%S').strftime('%Y%m%d%H%M')             
             sate = 'FY-' + fname[2:4]
             area = 'CHINA REGION '+ fname[74:78] + ' M'
             fid = dict(time=time, area=area, name=entire_fname, sate=sate, type='FY4_HDF', fdir=filedir)
             minor_idlist.append(sate + time)
             options.append(fid)         
        
             '''
            if sate + time not in minor_idlist:
                fid = dict(time=time, area=area, name=entire_fname, sate=sate, type='FY4_HDF', fdir=filedir)
                minor_idlist.append(sate + time)
                options.append(fid)'''
        return options
        
    @staticmethod
    def get_channel_filename(fname, channel):
        C0500 = {'band02':2}
        C1000 = {'band01':1, 'band02':2, 'band03':3}
        C2000 = {'band01':1, 'band02':2, 'band03':3, 'band04':4, 'band05':5, 'band06':6, 'band07':7}
        C4000 = {'band01':1, 'band02':2, 'band03':3, 'band04':4, 'band05':5, 'band06':6, 'band07':7, 'band08':8, 'band09':9, 'band10':10, 'band11':11, 'band12':12, 'band13':13, 'band14':14}
        if fname[74:78] == '0500' :
         if channel in C0500.keys() :
          return fname[:] 
         else :
          return 'X'
        elif fname[74:78] == '1000' :
         if channel in C1000.keys() :
          return fname[:] 
         else :
          return 'X'
        elif fname[74:78] == '2000' :
         if channel in C2000.keys() :
          return fname[:] 
         else :
          return 'X' 
        elif fname[74:78] == '4000' :
         if channel in C4000.keys() :
          return fname[:] 
         else :
          return 'X'           
        #return fname[:]   

    def get_channel_info(self):
        filedir, fname = split(self.fid['name'])
        channel_flag = {}
        for channel in self.crossdict.keys():
            channel_flag[channel] = isfile(join(filedir, FY4_HDF.get_channel_filename(fname, channel)))
        return channel_flag

    def extract(self, channel, georange):
        filedir, fname = split(self.fid['name'])
        latmin, latmax, lonmin, lonmax = georange
        self.georange = georange
        filename = join(filedir, FY4_HDF.get_channel_filename(fname, channel))
        data = _HDFReader(fname,channel,georange,self)
        return data
        

    def geocoord(self, channel):
        filedir, fname = split(self.fid['name'])
        filename = join(filedir, FY4_HDF.get_channel_filename(fname, channel))    
        lon,lat=_GetLonLat(channel,fname,self)
        print(np.amax(lon),np.amin(lon),np.amax(lat),np.amin(lat))
        #lon[lon < 0] = 360 + lon[lon < 0]
        return lon, lat   


         
def _HDFReader(fname,channel,georange,self):
    hdf5_obj=h5py.File(fname, 'r')
    #print('ddddddddddd',channel[-2:])
    SubLon = hdf5_obj.attrs['NOMCenterLon'][0]
    BL = hdf5_obj.attrs['Begin Line Number'][0]
    EL = hdf5_obj.attrs['End Line Number'][0]
    BC = hdf5_obj.attrs['Begin Pixel Number'][0]
    EC = hdf5_obj.attrs['End Pixel Number'][0]
    self.lineN,self.lineS,self.lineW,self.lineE=_CalcSegNeeded(georange,SubLon,channel,fname)
    
    y1=self.lineN-BL
    y2=self.lineS-BL
    x1=self.lineW-BC
    x2=self.lineE-BC
    if y1<BL:
     y1=BL
    if y2>EL:
     y2=EL
    if x1<BC:
     x1=BC
    if x2>EC:
     x2=EC   
    #print(y1,y2,x1,x2)    
    raw = hdf5_obj['NOMChannel%s'%(channel[-2:])][y1:y2,x1:x2]
    #print(raw.shape,fname[15:19])
    raw[raw > 65530] = 0
    table = hdf5_obj['CALChannel%s'%(channel[-2:])][:]
    count=np.arange(0,4096,1)
    dictionary = dict(zip(count,table))  
    
    if channel == 'band01' or channel == 'band02' or channel == 'band03' :
     data=np.vectorize(dictionary.__getitem__)(raw)
    else:
     data=np.vectorize(dictionary.__getitem__)(raw) - 273.15
    del raw     
    return data

def _GetLonLat(channel,fname,self):
    hdf5_obj=h5py.File(fname, 'r')
    BL = hdf5_obj.attrs['Begin Line Number'][0]
    EL = hdf5_obj.attrs['End Line Number'][0]
    BC = hdf5_obj.attrs['Begin Pixel Number'][0]
    EC = hdf5_obj.attrs['End Pixel Number'][0]
    print(fname[74:78],fname,'aads')
    if fname[74:78] == '0500' :
     CFAC = 81865099
     LFAC = 81865099
     COFF = 10991.5
     LOFF = 10991.5
    elif fname[74:78] == '1000' :
     CFAC = 40932549
     LFAC = 40932549
     COFF = 5495.5
     LOFF = 5495.5
    elif fname[74:78] == '2000' :
     CFAC = 20466274
     LFAC = 20466274
     COFF = 2747.5
     LOFF = 2747.5   
    elif fname[74:78] == '4000' :
     LOFF=1373.5
     LFAC=10233137
     COFF=1373.5
     CFAC=10233137

    SubLon = hdf5_obj.attrs['NOMCenterLon'][0]
    print(LOFF,LFAC,COFF,CFAC,SubLon)
    

    EarthPolarRadius=6356.7523
    EarthEquatorRadius=6378.137
    EarthConst1=(EarthEquatorRadius**2-EarthPolarRadius**2)/EarthEquatorRadius**2
    EarthConst2=EarthPolarRadius**2/EarthEquatorRadius**2
    EarthConst3=EarthEquatorRadius**2/EarthPolarRadius**2
    #Distance=(hdf5_obj.attrs['NOMSatHeight'][0]/1000 + hdf5_obj.attrs['dEA'][0])
    Distance = 42164
    EarthConstStd=Distance**2-EarthEquatorRadius**2
    
    DEGTORAD = np.pi / 180.
    RADTODEG = 180. / np.pi
    SCLUNIT = np.power(2., -16)
    DIS = Distance
    CON = EarthConst3
    #Calculation
    y1=self.lineN
    y2=self.lineS
    x1=self.lineW
    x2=self.lineE
    #print(y1,y2,x1,x2) 
    if y1<BL:
     y1=BL
    if y2>EL:
     y2=EL
    if x1<BC:
     x1=BC
    if x2>EC:
     x2=EC 
    lines = np.arange(y1,y2)
    columns = np.arange(x1,x2)
    
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
    print(lats.shape,lons.shape)
    return np.ma.masked_invalid(lons), np.ma.masked_invalid(lats)  

def _CalcSegNeeded(georange,SubLon,channel,fname):
    latmin, latmax, lonmin, lonmax = georange
    DEGTORAD = np.pi / 180.
    RADTODEG = 180. / np.pi
    SCLUNIT = np.power(2., -16)
    lons, lats = np.meshgrid(np.array([lonmin, lonmax],dtype='float'), np.array([latmin, latmax],dtype='float'))
    lons *= DEGTORAD
    lats *= DEGTORAD
    phi = np.arctan(_data_['EarthConst2'] * np.tan(lats))
    Re = _data_['EarthPolarRadius'] / np.sqrt(1 - _data_['EarthConst1'] * np.square(np.cos(phi)))
    r1 = _data_['Distance'] - Re * np.cos(phi) * np.cos(lons - SubLon * DEGTORAD)
    r2 = - Re * np.cos(phi) * np.sin(lons - SubLon * DEGTORAD)
    r3 = Re * np.sin(phi)
    #seeable = r1 * (r1 - _data_['Distance']) + np.square(r2) + np.square(r3)
    del lats, lons, phi, Re
    rn = np.sqrt(np.square(r1) + np.square(r2) + np.square(r3))
    x = np.arctan2(-r2, r1) * RADTODEG
    y = np.arcsin(-r3 / rn) * RADTODEG
    del r1, r2, r3, rn
    columns = _data_['COFF'] + x * SCLUNIT * _data_['CFAC']
    lines = _data_['LOFF'] + y * SCLUNIT * _data_['LFAC']
    lineN, lineS = np.amin(lines), np.amax(lines)
    columnW, columnE = np.amin(columns), np.amax(columns)
    if fname[74:78] == '0500' :
     scale = 4
     edge = 100     
    elif fname[74:78] == '1000' :
     scale = 2
     edge = 100
    elif fname[74:78] == '2000' :
     scale = 1
     edge = 100
    elif fname[74:78] == '4000' :
     scale = 0.5
     edge = 100
    print(int(lineN*scale-edge), int(lineS*scale+edge), int(columnW*scale-edge), int(columnE*scale+edge))
    return int(lineN*scale-edge), int(lineS*scale+edge), int(columnW*scale-edge), int(columnE*scale+edge)
    
'''        
    def VISGEOCALC(self, channel, georange):
        latmin, latmax, lonmin, lonmax = georange
        ncfile = self.ncfiles[channel]
        latnorth = ncfile.variables['geospatial_lat_lon_extent'].geospatial_northbound_latitude
        latsouth = ncfile.variables['geospatial_lat_lon_extent'].geospatial_southbound_latitude
        latconstant = 21696/(latnorth-latsouth)
        lonwest = ncfile.variables['geospatial_lat_lon_extent'].geospatial_westbound_longitude
        if lonwest < 0 :
         lonwest = 360 + lonwest
        loneast = ncfile.variables['geospatial_lat_lon_extent'].geospatial_eastbound_longitude
        if loneast < 0 :
         loneast = 360 + loneast
        lonconstant = 21696/(loneast-lonwest)
        if channel == 'band02' :
         (latmax-latmin)*latconstant
         lata = (latnorth-latmax)*latconstant - 3000
         latb = (latnorth-latmin)*latconstant + 3000
         lona = (lonmin-lonwest)*latconstant - 3000
         lonb = (lonmax-lonwest)*latconstant + 3000
         print('loooooooooon',int(lona),int(lonb),int(lata),int(latb))
         return int(lona),int(lonb),int(lata),int(latb)'''    
    
    
    
    