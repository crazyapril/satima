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

class GOES16_NC:
    '''channels = {'vis':('band01', 'band02', 'band03'), 'nir':('band04', 'band05', 'band06'),
                'swir':('band07', 'band11', 'band12'), 'wv':('band08', 'band09', 'band10'),
                'ir':('band13', 'band14', 'band15'), 'ir2':('band16',)}'''
    remap = True
    
    composite = {}
    '''
    channels = {'vis':('vis',), 'nir':None, 'swir':('swir',), 'ir2':('ir2',), 'wv':('wv',), 'ir':('ir',)}
    resolution = {'vis':1., 'ir':4., 'ir2':4., 'wv':4., 'swir':4.}
    crossdict = {'vis':1, 'ir':13, 'swir':2, 'wv':3, 'ir2':5}'''
    channels = {'vis':('band01','band02',), 'nir':('band03','band04','band05','band06',), 'swir':('band07','band11','band12',), 'ir2':('band16',), 'wv':('band08','band09','band10',), 'ir':('band13','band14','band15',)}
    resolution = {'band01':1., 'band02':0.5, 'band03':1, 'band04':2., 'band05':1., 'band06':2.,
                  'band07':2., 'band08':2., 'band09':2., 'band10':2., 'band11':2., 'band12':2.,
                  'band13':2., 'band14':2., 'band15':2., 'band16':2.}
    crossdict = {'band01':1, 'band02':2, 'band03':3, 'band04':4, 'band05':5, 'band06':6, 'band07':7, 'band08':8, 'band09':9, 'band10':10, 'band11':11, 'band12':12, 'band13':13, 'band14':14, 'band15':15, 'band16':16}
    def __init__(self, fid):
        self.fid = fid
        self.ncfiles = dict()
        self.lats = ''
        self.lons = ''


    @staticmethod
    def search(filedir):
        filelist = glob(join(filedir, 'OR_ABI-L1b-Rad*-M*C*_G*_s*_e*_c*.nc'))
        options = []
        minor_idlist = []
        for entire_fname in filelist:
            fname = split(entire_fname)[1]            
            if fname[14:16] == 'C-' :
             time = datetime.strptime(fname[27:40], '%Y%j%H%M%S').strftime('%Y%m%d%H%M')             
             sate = 'GOES-' + fname[23:25]
             area = '(CONUS)'+ ' BAND'+ fname[19:21]
             flag = []
             fid = dict(time=time, area=area, name=entire_fname, sate=sate, type='GOES16_NC', fdir=filedir,flag=flag,channel='band'+ fname[19:21])
             minor_idlist.append(sate + time)
             options.append(fid)
            elif fname[14:16] == 'F-' :
             time = datetime.strptime(fname[27:40], '%Y%j%H%M%S').strftime('%Y%m%d%H%M')
             sate = 'GOES-' + fname[23:25]
             area = 'BAND'+ fname[19:21]
             flag = []
             fid = dict(time=time, area=area, name=entire_fname, sate=sate, type='GOES16_NC', fdir=filedir,flag=flag,channel='band'+ fname[19:21])
             minor_idlist.append(sate + time)
             options.append(fid)             
            elif fname[14:16] == 'M1' :
             time = datetime.strptime(fname[28:41], '%Y%j%H%M%S').strftime('%Y%m%d%H%M')
             sate = 'GOES-' + fname[24:26]
             area = '(MESO1)'+ ' BAND'+ fname[20:22]
             flag = ['AUTOGEO']
             fid = dict(time=time, area=area, name=entire_fname, sate=sate, type='GOES16_NC', fdir=filedir,flag=flag,channel='band'+ fname[20:22])
             minor_idlist.append(sate + time)
             options.append(fid)             
            elif fname[14:16] == 'M2' :
             time = datetime.strptime(fname[28:41], '%Y%j%H%M%S').strftime('%Y%m%d%H%M')
             sate = 'GOES-' + fname[24:26]
             area = '(MESO2)'+ ' BAND'+ fname[20:22]
             flag = ['AUTOGEO']
             fid = dict(time=time, area=area, name=entire_fname, sate=sate, type='GOES16_NC', fdir=filedir,flag=flag,channel='band'+ fname[20:22])
             minor_idlist.append(sate + time)
             options.append(fid)             
             '''
            if sate + time not in minor_idlist:
                fid = dict(time=time, area=area, name=entire_fname, sate=sate, type='GOES16_NC', fdir=filedir)
                minor_idlist.append(sate + time)
                options.append(fid)'''
        return options
        
    @staticmethod
    def get_channel_filename(fname, channel):
        #print(channel[-1:])
        #print(channel[-2:],'ssssssssssss')
        if fname[14:16] == 'C-' :
           return fname[:19] + channel[-2:] + fname[21:]
        elif fname[14:16] == 'F-' :
           return fname[:19] + channel[-2:] + fname[21:]
        elif fname[14:16] == 'M1' :
           return fname[:20] + channel[-2:] + fname[22:]
        elif fname[14:16] == 'M2' :
           return fname[:20] + channel[-2:] + fname[22:]       

    def get_channel_info(self):
        filedir, fname = split(self.fid['name'])
        channel_flag = {}
        for channel in self.crossdict.keys():
            channel_flag[channel] = isfile(join(filedir, GOES16_NC.get_channel_filename(fname, channel)))
        return channel_flag

    def extract(self, channel, georange):
        filedir, fname = split(self.fid['name'])
        latmin, latmax, lonmin, lonmax = georange
        print(georange)
        if channel in self.ncfiles.keys():
            ncfile = self.ncfiles[channel]
        else:
            ncfile = Dataset(join(filedir, GOES16_NC.get_channel_filename(fname, channel)))
            self.ncfiles[channel] = ncfile
        print(channel)
        if channel == 'band02' and fname[14:16] == 'F-':
            #lona,lonb,lata,latb = self.VISGEOCALC(channel, georange)            
            lona,lonb,lata,latb = 1500,10500,1500,10500
            data=ncfile.variables['Rad'][lona:lonb,lata:latb]
            #data=ncfile.variables['Rad']
            kappa=ncfile.variables['kappa0'][:]
            raw=data*kappa
            print(raw.shape)
            del data
            return raw
        elif channel == 'band02' and fname[14:16] != 'F-':
            #lona,lonb,lata,latb = self.VISGEOCALC(channel, georange)            
            #lona,lonb,lata,latb = 15000,22000,0,7000
            data=ncfile.variables['Rad'][:]
            #data=ncfile.variables['Rad']
            kappa=ncfile.variables['kappa0'][:]
            raw=data*kappa
            print(raw.shape)
            del data
            return raw            
        elif channel == 'band01' or channel == 'band03'  or channel == 'band04'  or channel == 'band05'  or channel == 'band06' :
            #lona,lonb,lata,latb = self.VISGEOCALC(channel, georange)            
            data=ncfile.variables['Rad']
            kappa=ncfile.variables['kappa0'][:]
            raw=data*kappa
            del data
            print(raw.shape)
            return raw            
        else :
            data=ncfile.variables['Rad'][:]
            fk1=ncfile.variables['planck_fk1'][:]
            fk2=ncfile.variables['planck_fk2'][:]
            bc1=ncfile.variables['planck_bc1'][:]
            bc2=ncfile.variables['planck_bc2'][:]
            raw=(fk2/np.log((fk1/data)+1)-bc1)/bc2-273.15
            print(raw.shape)
            return raw

    def geocoord(self, channel,AUTO=False):
        if AUTO:
         filedir, fname = split(self.fid['name'])
         ncfile = Dataset(join(filedir, GOES16_NC.get_channel_filename(fname, channel)))
        else:
         ncfile = self.ncfiles[channel]     
         filedir, fname = split(self.fid['name'])
        if channel == 'band02' and fname[14:16] == 'F-':
          #lona,lonb,lata,latb = self.VISGEOCALC(channel, georange)          
          lona,lonb,lata,latb = 1500,10500,1500,10500
          x=ncfile.variables['x'][lona:lonb]
          y=ncfile.variables['y'][lata:latb] 
          #x=ncfile.variables['x'][:]
          #y=ncfile.variables['y'][:]             
          x,y=np.meshgrid(x,y)
       
        else :
          x=ncfile.variables['x'][:]
          y=ncfile.variables['y'][:]
          x,y=np.meshgrid(x,y)
        #Variables
        req = ncfile.variables['goes_imager_projection'].semi_major_axis
        flatten = ncfile.variables['goes_imager_projection'].inverse_flattening
        rpol = ncfile.variables['goes_imager_projection'].semi_minor_axis
        e = 0.0818191910435
        H = ncfile.variables['goes_imager_projection'].perspective_point_height + req
        a = (np.sin(x))**2 + ((np.cos(x))**2)*((np.cos(y))**2 + ((req/rpol)**2)*(np.sin(y))**2)
        b = -2*H*np.cos(x)*np.cos(y)
        c = H**2 - req**2
        rs = (-b-np.sqrt(b**2-4*a*c))/2*a
        del a,b,c
        Sx = rs*np.cos(x)*np.cos(y)
        Sy = -rs*np.sin(x)
        Sz = rs*np.cos(x)*np.sin(y)
        constant1 = req**2/rpol**2
        constant2 = H-Sx
        del x,y
        lamda0 = np.radians(ncfile.variables['goes_imager_projection'].longitude_of_projection_origin)
        lats2 = np.degrees(np.arctan((constant1*Sz)/np.sqrt(constant2**2+Sy**2)))
        del Sx,Sz
        lons2 = np.degrees(lamda0 - np.arctan(Sy/constant2))
        lons2[lons2 < 0]= 360 + lons2[lons2 < 0 ]

        print(np.nanmax(lons2),np.nanmin(lons2),np.nanmax(lats2),np.nanmin(lats2))
        print(lons2.shape,lats2.shape)
        return lons2, lats2
        
    def get_auto_georange(self):
        channel = self.fid['channel']
        lons, lats = self.geocoord(channel,AUTO=True)
        return lats.min(), lats.max(), lons.min(), lons.max()        
        
        
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
         
         
         
    
    
    
    
    