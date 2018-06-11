import bz2
import datetime
from glob import glob
from os.path import isfile, join, split

import numpy as np

from ._sate_param import _hmw_data_for_seg_calc as _data_
from .utils import reglob

__all__ = ['READABLE', 'HIMAWARI_HSD']

READABLE = True

class HIMAWARI_HSD:

    remap = True
    channels = {'vis':('band01', 'band02', 'band03'), 'nir':('band04', 'band05', 'band06'),
                'swir':('band07', 'band11', 'band12'), 'wv':('band08', 'band09', 'band10'),
                'ir':('band13', 'band14', 'band15'), 'ir2':('band16',)}
    resolution = {'band01':1., 'band02':1., 'band03':0.5, 'band04':1., 'band05':2., 'band06':2.,
                  'band07':2., 'band08':2., 'band09':2., 'band10':2., 'band11':2., 'band12':2.,
                  'band13':2., 'band14':2., 'band15':2., 'band16':2.}
    composite = {'True Color':('band03', 'band02', 'band01')}

    def __init__(self, fid):
        self.fid = fid
        self.hsd_dict = {}
        np.seterr(invalid='ignore')
		
    @classmethod
    def search(cls, filedir):
        #HS_H0[0:8/9]_[1:YYYYMMDD]_[2:HHMM]_B[3:Band]_[4:Area]_R[5:Res]_S[6:Seg].DAT[.bz2]
        namefmt = 'HS_H0*_*_*_B*_*_R*_S*.DAT'
        filelist = reglob(filedir, namefmt) + reglob(filedir, namefmt+'.bz2')
        idlist = dict()
        for fname, info in filelist:
            if info[4] == 'FLDK':
                area = ''
            elif info[4][:2] == 'R3':
                area = 'TARGET'
            elif info[4][:2] == 'JP':
                area = 'JAPAN'
            time = datetime.datetime.strptime(info[1] + info[2], '%Y%m%d%H%M')
            if area:
                time = time + datetime.timedelta(minutes=2.5*int(info[4][2:])-2.5)
                flag = ['AUTOGEO']
            else:
                flag = []
            band = 'band' + info[3]
            minor_id = time, area
            if minor_id not in idlist:
                fid = dict(time=time.strftime('%Y%m%d%H%M'), fdir=filedir, type='himawari_hsd',
                           flag=flag, fcls=cls, sate='HIMAWARI-'+info[0][-1], channel={band:fname},
                           area=area, special_channel=dict(), name=join(filedir, fname))
                idlist[minor_id] = fid
            else:
                idlist[minor_id]['channel'][band] = fname
        return list(idlist.values())

    def get_channel_info(self):
        filedir, fname = split(self.fid['name'])
        channel_flag = {}
        for channel in self.resolution.keys():
            for seg in range(1,11):
                cflag = isfile(join(filedir, _get_segment_filename(_get_channel_filename(fname, channel), seg)))
                channel_flag[channel] = cflag
                if cflag:
                    break
        return channel_flag

    def extract(self, channel, georange):
        if 'AUTOGEO' in self.fid.get('flag', []):
            hsd, raw = _HSDReader(join(self.fid['fdir'], self.fid['channel'][channel]))
        else:
            seglist, corner = _CalcSegNeeded(georange)
            print(channel, corner)
            hsd, raw = _CombineHSD(self.fid, seglist, corner, channel)
        if channel not in self.hsd_dict:
            self.hsd_dict[channel] = hsd
        return _Calibration(self.hsd_dict[channel], raw)

    def geocoord(self, channel):
        return _GetLonLat(self.hsd_dict[channel])

    def get_auto_georange(self):
        channel = next(iter(self.fid['channel']))
        hsd, raw = _HSDReader(join(self.fid['fdir'], self.fid['channel'][channel]))
        self.hsd_dict[channel] = hsd
        lons, lats = _GetLonLat(hsd)
        return lats.min(), lats.max(), lons.min(), lons.max()

_BLOCK_01 = np.dtype([('HeaderBlockNumber', 'u1'),
                     ('BlockLength', 'u2'),
                     ('TotalNumberOfHeaderBlocks', 'u2'),
                     ('ByteOrder', 'u1'),
                     ('SatelliteName', 'a16'),
                     ('ProcessingCenterName', 'a16'),
                     ('ObservationArea', 'a4'),
                     ('OtherObservationInformation', 'a2'),
                     ('ObservationTimeline', 'u2'),
                     ('ObservationStartTime', 'f8'),
                     ('ObservationEndTime', 'f8'),
                     ('FileCreationTime', 'f8'),
                     ('TotalHeaderLength', 'u4'),
                     ('TotalDataLength', 'u4'),
                     ('QualityFlag1', 'u1'),
                     ('QualityFlag2', 'u1'),
                     ('QualityFlag3', 'u1'),
                     ('QualityFlag4', 'u1'),
                     ('FileFormatVersion', 'a32'),
                     ('FileName', 'a128'),
                     ('Spare', 'a40')])

_BLOCK_02 = np.dtype([('HeaderBlockNumber', 'u1'),
                     ('BlockLength', 'u2'),
                     ('NumberOfBitsPerPixel', 'u2'),
                     ('NumberOfColumns', 'u2'),
                     ('NumberOfLines', 'u2'),
                     ('CompressionFlag', 'u1'),
                     ('Spare', 'a40')])

_BLOCK_03 = np.dtype([('HeaderBlockNumber', 'u1'),
                     ('BlockLength', 'u2'),
                     ('SubLon', 'f8'),
                     ('CFAC', 'u4'),
                     ('LFAC', 'u4'),
                     ('COFF', 'f4'),
                     ('LOFF', 'f4'),
                     ('Distance', 'f8'),
                     ('EarthEquatorialRadius', 'f8'),
                     ('EarthPolarRadius', 'f8'),
                     ('EarthConst1', 'f8'),
                     ('EarthConst2', 'f8'),
                     ('EarthConst3', 'f8'),
                     ('EarthConstStd', 'f8'),
                     ('ResamplingTypes', 'u2'),
                     ('ResamplingSize', 'u2'),
                     ('Spare', 'a40')])

_BLOCK_05 = np.dtype([('HeaderBlockNumber', 'u1'),
                     ('BlockLength', 'u2'),
                     ('BandNumber', 'u2'),
                     ('CentralWaveLength', 'f8'),
                     ('ValidNumberOfBitsPerPixel', 'u2'),
                     ('CountValueOfErrorPixels', 'u2'),
                     ('CountValueOfPixelsOutsideScanArea', 'u2'),
                     ('Gain', 'f8'),
                     ('Constant', 'f8')])

_InfraredBand = np.dtype([('c0', 'f8'),
                         ('c1', 'f8'),
                         ('c2', 'f8'),
                         ('C0', 'f8'),
                         ('C1', 'f8'),
                         ('C2', 'f8'),
                         ('c', 'f8'),
                         ('h', 'f8'),
                         ('k', 'f8'),
                         ('Spare', 'a40')])

_VisibleBand = np.dtype([('c*', 'f8'),
                        ('Spare', 'a104')])

_BLOCK_07 = np.dtype([('HeaderBlockNumber', 'u1'),
                     ('BlockLength', 'u2'),
                     ('TotalNumberOfSegments', 'u1'),
                     ('SegmentSequenceNumber', 'u1'),
                     ('FirstLineNumber', 'u2'),
                     ('Spare', 'a40')])

_Header = np.dtype([('HeaderBlockNumber', 'u1'),
                ('BlockLength', 'u2')])

def _get_channel_filename(fname, channel):
        return fname[:22] + channel[-2:] + fname[24:31] + '%02d' % (HIMAWARI_HSD.resolution[channel]*10) + fname[33:]

def _get_segment_filename(fname, seg_no):
        return fname[:35] + '%02d' % (seg_no) + fname[37:]

def _HSDReader(fname, corner=None):
    HSD = {}
    if not isfile(fname):
        raise IOError('Require ' + fname)
    if fname.endswith('.bz2'):
        f = bz2.open(fname, mode='rb')
    else:
        f = open(fname, 'rb')
    HSD['BLOCK_01'] = np.fromstring(f.read(282), dtype=_BLOCK_01)
    HSD['BLOCK_02'] = np.fromstring(f.read(50), dtype=_BLOCK_02)
    HSD['BLOCK_03'] = np.fromstring(f.read(127), dtype=_BLOCK_03)
    _LeapBlock(f, 1)
    HSD['BLOCK_05'] = np.fromstring(f.read(35), dtype=_BLOCK_05)
    if HSD['BLOCK_05']['BandNumber'] <= 6:
        HSD['VisibleBand'] = np.fromstring(f.read(112), dtype=_VisibleBand)
    else:
        HSD['InfraredBand'] = np.fromstring(f.read(112), dtype=_InfraredBand)
    _LeapBlock(f, 1)
    HSD['BLOCK_07'] = np.fromstring(f.read(47), dtype=_BLOCK_07)
    _LeapBlock(f, 4)
    lines = HSD['BLOCK_02']['NumberOfLines'].item()
    columns = HSD['BLOCK_02']['NumberOfColumns'].item()
    raw = np.ma.masked_greater(np.fromstring(f.read(), dtype='uint16').reshape((lines, columns)), 65530)
    if corner:
        column_west = int(raw.shape[1] * corner[0])
        column_east = int(raw.shape[1] * corner[1])
        print(column_west, column_east)
    else:
        column_west = 0
        column_east = raw.shape[1] - 1
    HSD['ColumnBoundary'] = (column_west, column_east)
    f.close()
    return HSD, raw[:,column_west:column_east]

def _LeapBlock(f, n):
    for i in range(n):
        tmparr = np.fromstring(f.read(3), dtype=_Header)
        f.seek(tmparr['BlockLength'].item()-3, 1)
    del tmparr

def _Calibration(HSD, raw):
    if HSD['BLOCK_05']['BandNumber'] <= 6:
        return _VISCalibration(HSD, raw)
    else:
        return _IRCalibration(HSD, raw)

def _IRCalibration(HSD, raw):
    lam = HSD['BLOCK_05']['CentralWaveLength'] * 1e-6
    gain = HSD['BLOCK_05']['Gain']
    const = HSD['BLOCK_05']['Constant']
    c = HSD['InfraredBand']['c']
    k = HSD['InfraredBand']['k']
    h = HSD['InfraredBand']['h']
    c0 = HSD['InfraredBand']['c0']
    c1 = HSD['InfraredBand']['c1']
    c2 = HSD['InfraredBand']['c2']
    const1 = h * c / (k * lam)
    const2 = 2 * h * np.power(c, 2) * np.power(lam, -5)
    I = (gain * raw + const) * 1e6
    EBT = const1 / np.log1p(const2 / I)
    return c0 + c1 * EBT + c2 * np.power(EBT, 2) - 273.15

def _VISCalibration(HSD, raw):
    gain = HSD['BLOCK_05']['Gain']
    const = HSD['BLOCK_05']['Constant']
    c = HSD['VisibleBand']['c*']
    return c * gain * raw + c * const

def _GetLonLat(HSD):
    #参考JMA提供代码编写
    #Constants
    DEGTORAD = np.pi / 180.
    RADTODEG = 180. / np.pi
    SCLUNIT = np.power(2., -16)
    DIS = HSD['BLOCK_03']['Distance']
    CON = HSD['BLOCK_03']['EarthConst3']
    #Calculation
    lines = np.arange(HSD['BLOCK_07']['FirstLineNumber'], HSD['BLOCK_07']['FirstLineNumber'] + HSD['BLOCK_02']['NumberOfLines'])
    columns = np.arange(*HSD['ColumnBoundary'])
    xx, yy = np.meshgrid(columns, lines)
    del lines, columns
    x = DEGTORAD * (xx - HSD['BLOCK_03']['COFF']) / (SCLUNIT * HSD['BLOCK_03']['CFAC'])
    y = DEGTORAD * (yy - HSD['BLOCK_03']['LOFF']) / (SCLUNIT * HSD['BLOCK_03']['LFAC'])
    del xx, yy
    Sd = np.sqrt(np.square(DIS * np.cos(x) * np.cos(y)) - (np.square(np.cos(y)) + CON * np.square(np.sin(y))) * HSD['BLOCK_03']['EarthConstStd'])
    Sn = (DIS * np.cos(x) * np.cos(y) - Sd) / (np.square(np.cos(y)) + CON * np.square(np.sin(y)))
    S1 = DIS - Sn * np.cos(x) * np.cos(y)
    S2 = Sn * np.sin(x) * np.cos(y)
    S3 = -Sn * np.sin(y)
    Sxy = np.sqrt(np.square(S1) + np.square(S2))
    del x, y, Sd, Sn
    lons = RADTODEG * np.arctan2(S2, S1) + HSD['BLOCK_03']['SubLon']
    lats = np.ma.masked_outside(RADTODEG * np.arctan(CON * S3 / Sxy), -90., 90.)
    return np.ma.masked_invalid(lons), np.ma.masked_invalid(lats),

def _CalcSegNeeded(georange):
    latmin, latmax, lonmin, lonmax = georange
    DEGTORAD = np.pi / 180.
    RADTODEG = 180. / np.pi
    SCLUNIT = np.power(2., -16)
    lons, lats = np.meshgrid(np.array([lonmin, lonmax],dtype='float'), np.array([latmin, latmax],dtype='float'))
    lons *= DEGTORAD
    lats *= DEGTORAD
    phi = np.arctan(_data_['EarthConst2'] * np.tan(lats))
    Re = _data_['EarthPolarRadius'] / np.sqrt(1 - _data_['EarthConst1'] * np.square(np.cos(phi)))
    r1 = _data_['Distance'] - Re * np.cos(phi) * np.cos(lons - _data_['SubLon'] * DEGTORAD)
    r2 = - Re * np.cos(phi) * np.sin(lons - _data_['SubLon'] * DEGTORAD)
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
    return list(range(int(lineN-1)//550+1, int(lineS-1)//550+2)), [columnW/5500-0.001, columnE/5500+0.001]

def _CombineHSD(fid, seg_range, corner, channel):
    filedir, fname = split(fid['name'])
    datalist = [_HSDReader(join(filedir, _get_channel_filename(_get_segment_filename(fname, seg), channel)), corner) for seg in seg_range]
    NewHSD = datalist[0][0]
    NewHSD['BLOCK_02']['NumberOfLines'] = np.sum([data[0]['BLOCK_02']['NumberOfLines'] for data in datalist])
    #NewHSD['BLOCK_07']['FirstLineNumber'] = np.amin([HSD['BLOCK_07']['FirstLineNumber'] for HSD in HSDList])
    #HSDSequence = []
    #for HSD in HSDList:
    #    HSDSequence.append([HSD['BLOCK_07']['FirstLineNumber'], HSD['Data']])
    #HSDSequence.sort()
    raw = np.concatenate(tuple((data[1] for data in datalist)))
    return NewHSD, raw
