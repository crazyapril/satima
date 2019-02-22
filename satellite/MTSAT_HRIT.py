from glob import glob
from os.path import isfile, split, join
import numpy as np
from ._sate_param import _data_
import struct
from datetime import datetime, timedelta
from numpy import loadtxt,savetxt,zeros
READABLE = True

class MTSAT_HRIT:

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
        filelist = glob(join(filedir, 'HRIT_MTSAT*'))
        options = []
        minor_idlist = []
        for entire_fname in filelist:
         fname = split(entire_fname)[1]
         time = fname[12:20] + fname[21:25]
         if time not in minor_idlist:
            fid = {}
            fid['time'] = time
            fid['name'] = entire_fname
            if int(fname[10:11])==1 : 
             fid['sate'] = 'MTSAT-1R'
            else : 
             fid['sate'] = 'MTSAT-2'
            fid['area'] = {'1':'Full', '2':'NH', '3':'SH'}[fname[29:30]]
            fid['type'] = 'mtsat_hrit'
            fid['fdir'] = filedir         
            minor_idlist.append(time)
            options.append(fid)
        return options

        
        
    def get_channel_filename(fname, channel):  
        return fname[:30] + MTSAT_HRIT.crossdict[channel]    


    def get_channel_info(self):
        filedir, fname = split(self.fid['name'])
        channel_flag = {}
        for channel in self.crossdict.keys():
            channel_flag[channel] = isfile(join(filedir, MTSAT_HRIT.get_channel_filename(fname, channel)))
        return channel_flag
  
    def extract(self, channel, georange):
        filedir, fname = split(self.fid['name'])
        latmin, latmax, lonmin, lonmax = georange 
        filename = join(filedir, MTSAT_HRIT.get_channel_filename(fname, channel))
        hrit, count = _HRITReader(filename)
        self.hrit_dict[channel] = hrit
        if self.fid['sate'] == 'MTSAT-2' and int(self.fid['time']) < 201003150000 :
           TBB=MTSAT2_Calibration_1(count,channel)
        elif self.fid['sate'] == 'MTSAT-2' and int(self.fid['time']) >= 201003150000 :
           TBB=MTSAT2_Calibration_2(count,channel)
        elif self.fid['sate'] == 'MTSAT-1R':
           TBB=MTSAT1R_Calibration(count,channel)
        print(TBB.shape)
        if channel  == 'vis' :  
           return TBB
        else :     
           return TBB-273.15
        
    def geocoord(self, channel):
        hrit = self.hrit_dict[channel]
        lon,lat=_GetLonLat(hrit)
        print(np.amax(lon),np.amin(lon),np.amax(lat),np.amin(lat))
        lon[lon < 0] = 360 + lon[lon < 0]
        return lon, lat   

#Header Type 0 - Primary Header
_Primary_Header = np.dtype([('Header_Type', '>u1'),
                            ('Header_Record_Length', '>u2'),
                            ('File_Type_Code', '>u1'),
                            ('Total_Header_Length', '>u4'),
                            ('Data_Field_Length', '>u8')])
#Header Type 1 - Image Structure
_Image_Structure = np.dtype([('Header_Type', '>u1'),
                             ('Header_Record_Length', '>u2'),
                             ('NB', '>u1'),
                             ('NC', '>u2'),
                             ('NL', '>u2'),
                             ('Compression_Flag', '>u1')])
#Header Type 2 - Image Navigation                             
_Image_Navigation = np.dtype([('Header_Type', '>u1'),
                              ('Header_Record_Length', '>u2'),
                              ('Projection_Name', '>a32'),
                              ('CFAC', '>i4'),
                              ('LFAC', '>i4'),
                              ('COFF', '>i4'),
                              ('LOFF', '>i4')])
#Header Type 3 - Image Data Function                     
_Image_Data_Function = np.dtype([('Header_Type', '>u1'),
                                 ('Header_Record_Length', '>u2')])
            
#Header Type 4 - Annotation
_Annotation = np.dtype([('Header_Type', '>u1'),
                        ('Header_Record_Length', '>u2')])                     
#Header Type 5 - Time Stamp                     
_Time_Stamp = np.dtype([('Header_Type', '>u1'),
                        ('Header_Record_Length', '>u2'),
                        ('CDS_P_Field', '>u1'),
                        ('CDS_T_Field_1', '>u2'),
                        ('CDS_T_Field_2', '>u4'),])

#Header Type 128 - Image Segment Identification                     
_Image_Segment_Identification = np.dtype([('Header_Type', '>u1'),
                                          ('Header_Record_Length', '>u2'),
                                          ('Image_Segm_Seq_No', '>u1'),
                                          ('Total_No_Image_Segm', '>u1'),
                                          ('Line_No_Image_Segm', '>u2')])
#Header Type 130 – Image Compensation Information Header
_Image_Compensation_Information_Header = np.dtype([('Header_Type', '>u1'),
                                                  ('Header_Record_Length', '>u2')])                      
#Header Type 131 – Image Observation Time Header
_Image_Observation_Time_Header = np.dtype([('Header_Type', '>u1'),
                                          ('Header_Record_Length', '>u2')]) 
#Header Type 132 – Image Quality Information Header
_Image_Quality_Information_Header = np.dtype([('Header_Type', '>u1'),
                                             ('Header_Record_Length', '>u2')])                   

        
def _GetLonLat(HRIT):

    NC = HRIT['Image_Structure']['NC']
    NL = HRIT['Image_Structure']['NL']
    CFAC = HRIT['Image_Navigation']['CFAC']
    LFAC = HRIT['Image_Navigation']['LFAC']
    COFF = HRIT['Image_Navigation']['COFF']
    LOFF = HRIT['Image_Navigation']['LOFF']
    proj_name = HRIT['Image_Navigation']['Projection_Name'][0].decode('utf-8')
    
    i1 = proj_name.find('(')
    i2 = proj_name.find(')')
    if i1 != -1 and i2 != -1:
       ssp = float(proj_name[i1+1:i2])
    else:
       ssp = None
    SubLon = ssp
    print(SubLon)

    EarthPolarRadius=6356.5838
    EarthEquatorRadius=6378.1690
    EarthConst1=(EarthEquatorRadius**2-EarthPolarRadius**2)/EarthEquatorRadius**2
    EarthConst2=EarthPolarRadius**2/EarthEquatorRadius**2
    EarthConst3=EarthEquatorRadius**2/EarthPolarRadius**2
    Distance=42164.0
    EarthConstStd=Distance**2-EarthEquatorRadius**2
    
    DEGTORAD = np.pi / 180.
    RADTODEG = 180. / np.pi
    SCLUNIT = np.power(2., -16)
    DIS = 42164.0
    CON = EarthConst3
    #Calculation
    lines = np.arange(1,NL+1)
    columns = np.arange(1,NC+1)
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

def _HRITReader(fname):
    HRIT = {}
    f = open(fname, 'rb')
    #0
    HRIT['Primary_Header'] = np.fromstring(f.read(16), dtype=_Primary_Header)
    #1
    HRIT['Image_Structure'] = np.fromstring(f.read(9), dtype=_Image_Structure)
    #2
    HRIT['Image_Navigation'] = np.fromstring(f.read(51), dtype=_Image_Navigation)
    #3
    HRIT['Image_Data_Function'] = np.fromstring(f.read(3), dtype=_Image_Data_Function)
    HRIT['Data_Definition_Block'] = f.read(HRIT['Image_Data_Function']['Header_Record_Length'][0]-3)
    #4
    HRIT['Annotation'] = np.fromstring(f.read(3), dtype=_Annotation)
    HRIT['Annotation_Text'] = f.read(HRIT['Annotation']['Header_Record_Length'][0]-3)
    print(HRIT['Annotation_Text'].decode('utf-8'))
    #5
    HRIT['Time_Stamp'] = np.fromstring(f.read(10), dtype=_Time_Stamp)
    #128
    HRIT['Image_Segment_Identification'] = np.fromstring(f.read(7), dtype=_Image_Segment_Identification)
    #130
    HRIT['Image_Compensation_Information_Header'] = np.fromstring(f.read(3), dtype=_Image_Compensation_Information_Header)
    HRIT['Image_Compensation_Information'] = f.read(HRIT['Image_Compensation_Information_Header']['Header_Record_Length'][0]-3)
    #131
    HRIT['Image_Observation_Time_Header'] = np.fromstring(f.read(3), dtype=_Image_Observation_Time_Header)
    HRIT['Image_Observation_Time'] = f.read(HRIT['Image_Observation_Time_Header']['Header_Record_Length'][0]-3)
    #132
    HRIT['Image_Quality_Information_Header'] = np.fromstring(f.read(3), dtype=_Image_Quality_Information_Header)
    HRIT['Image_Quality_Information'] = f.read(HRIT['Image_Quality_Information_Header']['Header_Record_Length'][0]-3)
    NC = HRIT['Image_Structure']['NC'].item()
    NL = HRIT['Image_Structure']['NL'].item()
    raw = np.fromstring(f.read(), dtype='>u2').reshape(NL, NC)
    f.close()
    return HRIT, raw
        
def MTSAT1R_Calibration(count,channel):
    DN=loadtxt('./satellite/MTSAT-1R.txt')[:,0]
    if channel == 'ir' :
     CTB=loadtxt('./satellite/MTSAT-1R.txt')[:,1]
    elif channel == 'ir2' :
     CTB=loadtxt('./satellite/MTSAT-1R.txt')[:,2]
    elif channel == 'wv' :
     CTB=loadtxt('./satellite/MTSAT-1R.txt')[:,3]
    elif channel == 'swir' :  
     CTB=loadtxt('./satellite/MTSAT-1R.txt')[:,4]
    elif channel == 'vis' :  
     CTB=loadtxt('./satellite/MTSAT-1R.txt')[:,5]     
    dictionary = dict(zip(DN,CTB))
    TBB=np.vectorize(dictionary.__getitem__)(count)   
    return TBB
    
def MTSAT2_Calibration_1(count,channel):
    DN=loadtxt('./satellite/MTSAT-2-1.txt')[:,0]
    if channel == 'ir' :
     CTB=loadtxt('./satellite/MTSAT-2-1.txt')[:,1]
    elif channel == 'ir2' :
     CTB=loadtxt('./satellite/MTSAT-2-1.txt')[:,2]
    elif channel == 'wv' :
     CTB=loadtxt('./satellite/MTSAT-2-1.txt')[:,3]
    elif channel == 'swir' :  
     CTB=loadtxt('./satellite/MTSAT-2-1.txt')[:,4]
    elif channel == 'vis' :  
     CTB=loadtxt('./satellite/MTSAT-2-1.txt')[:,5]     
    dictionary = dict(zip(DN,CTB))
    TBB=np.vectorize(dictionary.__getitem__)(count)   
    return TBB
    
def MTSAT2_Calibration_2(count,channel):
    DN=loadtxt('./satellite/MTSAT-2-2.txt')[:,0]
    if channel == 'ir' :
     CTB=loadtxt('./satellite/MTSAT-2-2.txt')[:,1]
    elif channel == 'ir2' :
     CTB=loadtxt('./satellite/MTSAT-2-2.txt')[:,2]
    elif channel == 'wv' :
     CTB=loadtxt('./satellite/MTSAT-2-2.txt')[:,3]
    elif channel == 'swir' :  
     CTB=loadtxt('./satellite/MTSAT-2-2.txt')[:,4]
    elif channel == 'vis' :  
     CTB=loadtxt('./satellite/MTSAT-2-2.txt')[:,5]     
    dictionary = dict(zip(DN,CTB))
    TBB=np.vectorize(dictionary.__getitem__)(count)   
    return TBB  
         