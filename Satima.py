from os import system
from os.path import join
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from datetime import datetime
from _FileSearch import FileSearch
from _Settings import Settings
from satellite import  _sate_param
from colormap import Colormap
from _confs import (convert_time_format, get_multi_choices_template, is_ascii,
                    safe_str2value, safe_input, unique_elements)

DEBUG = False

class Satima:

    def __init__(self, FileSearch):
        self.fid = FileSearch.choice
        self.settings = Settings().settings
        self.fileclass = FileSearch.fileclass[self.fid['type'].upper()](self.fid)
        self.channel_flag = self.fileclass.get_channel_info()
        self.common_imoptions = {'Common':[], 'Composite':[]}
        self.spec_imoptions = {'WV强化':{'channels':[],'colormaps':[]}, 'IR强化':{'channels':[],'colormaps':[]}}
        self.raw = {}
        self.run()

    def run(self):
        self.get_options()
        self.get_choices()
        self.get_georange()
        self.get_size()
        self.get_title()
        self.set_map()
        self.processor()

    def get_options(self):
        wv_colormaps, ir_colormaps = Colormap.get_colormap_list()
        for imtype in self.channel_flag.keys():
            if self.channel_flag[imtype]:
                self.common_imoptions['Common'].append(imtype)
        for imtype in self.fileclass.composite.keys():
            composite_flag = [self.channel_flag[channel_require] for channel_require in self.fileclass.composite[imtype]]
            if all(composite_flag):
                self.common_imoptions['Composite'].append(imtype)
        for wv in self.fileclass.channels['wv']:
            if self.channel_flag[wv]:
                self.spec_imoptions['WV强化']['channels'].append(wv)
                self.spec_imoptions['WV强化']['colormaps'] = wv_colormaps
        for ir in self.fileclass.channels['ir']:
            if self.channel_flag[ir]:
                self.spec_imoptions['IR强化']['channels'].append(ir)
                self.spec_imoptions['IR强化']['colormaps'] = ir_colormaps

    def get_choices(self):
        count = 0
        choicedict = {}
        print()
        for imcategory, imtypes in self.common_imoptions.items():
            if len(imtypes) == 0:
                continue
            print('[{}]'.format(imcategory), end='')
            for imtype in imtypes:
                count += 1
                print(' {0:d}.{1:s}'.format(count, imtype.upper()), end='')
                choicedict[count] = (imcategory, imtype)
            print()
        for imcategory, iminfo in self.spec_imoptions.items():
            num_of_channels = len(iminfo['channels'])
            if num_of_channels == 0:
                continue
            print('[{}]'.format(imcategory), end='')
            for imchannel in iminfo['channels']:
                if num_of_channels > 1:
                    print('[{}]'.format(imchannel), end='')
                for imtype in iminfo['colormaps']:
                    count += 1
                    print(' {0:d}.{1:s}'.format(count, imtype.upper()), end='')
                    choicedict[count] = (imchannel, imtype)
                print()
        raw_choices = get_multi_choices_template(choicedict.keys(), '输入要制作的图像种类：', sort=True)
        self.choices = [choicedict[choice] for choice in raw_choices]
        print()

    def get_georange(self):
        while True:
            latmin = safe_input('请输入最南纬度: ')
            latmax = safe_input('请输入最北纬度: ')
            lonmin = safe_input('请输入最西经度: ')
            lonmax = safe_input('请输入最东经度: ')
            if lonmin < 0:
                lonmin += 360
            if lonmax < 0:
                lonmax += 360
            if latmin >= latmax:
                print('纬度无效。')
                continue
            if lonmin >= lonmax:
                print('经度无效。')
                continue
            if not self.fileclass.remap:
                saterange = _sate_param._data_[self.fid['sate']]['georange']
                if latmax > saterange[0] or latmin < saterange[1] or lonmax > saterange[2] or lonmin < saterange[3]:
                    print('输入范围无效。')
                    print('参考范围：纬度 {}~{}，经度 {}~{}'.format(saterange[1], saterange[0], saterange[3], saterange[2]))
                    continue
            self.georange = (latmin, latmax, lonmin, lonmax)
            print()
            return

    def get_size(self):
        latmin, latmax, lonmin, lonmax = self.georange
        ratio = (latmax - latmin) / (lonmax - lonmin)
        print('默认图片大小为2000x{:d}。你可以自定义图片宽度，其高度将会按比例计算。'.format(int(2000*ratio)))
        width = safe_input('输入自定义图片宽度(500~4000)，回车使用默认值：', func=int, vrange=(500, 4000), default=2000)
        self.dpi = width / Settings.widthsize
        self.figsize = (Settings.widthsize, Settings.widthsize*ratio)
        print()

    def get_title(self):
        if self.settings['switch']['title']:
            self.title = input('输入图片标题：')
            ascii_flag = is_ascii(self.title)
            self.title_font = self.settings['font']['title'] if ascii_flag else self.settings['font']['title_nonascii']
            self.title_weight = self.settings['font']['title_weight'] if ascii_flag else self.settings['font']['title_nonascii_weight']
        self.basic_info = convert_time_format(self.fid['time']) + '    {} '.format(self.fid['sate'])
        print()

    def get_colormap(self, choice):
        imcategory, imtype = choice
        if imcategory == 'Common':
            if imtype in self.fileclass.channels['ir'] or imtype in self.fileclass.channels['wv']:
                return 'gray_r'
            return 'gray'
        cmap = Colormap.get_colormap(imtype)
        if cmap == 0:
            print(imtype.upper()+'色阶存在问题。将使用灰度色阶。')
            return 'gray_r'
        return cmap

    def get_imtype(self, choice):
        imcategory, imtype = choice
        if imcategory in ('Common', 'Composite'):
            return imtype
        return imcategory + '-' + imtype

    def set_map(self):
        if self.settings['switch']['coastline'] or self.settings['switch']['borderline'] or self.settings['switch']['latlon']:
            latmin, latmax, lonmin, lonmax = self.georange
            resolution = self.settings['image']['resolution']
            if resolution not in ('c', 'l', 'i', 'h', 'f'):
                resolution = 'i'
            self.map = Basemap(projection='cyl', llcrnrlat=latmin, urcrnrlat=latmax, llcrnrlon=lonmin, urcrnrlon=lonmax, resolution=resolution)

    def get_imager_kwargs(self, choice):
        imcategory, imtype = choice
        if imcategory == 'Composite':
            return {}
        cmap = self.get_colormap(choice)
        if imcategory == 'Common':
            if imtype in self.fileclass.channels['vis']:
                return {'vmin': 0., 'vmax': 1., 'cmap':cmap}
            elif imtype in self.fileclass.channels['ir']:
                return {'vmin': -80., 'vmax': 50., 'cmap':cmap}
        if imcategory not in ('Common', 'Composite'):
            return {'vmin': -100., 'vmax': 50., 'cmap':cmap}
        return {}

    def find_corresponding_channel(self, choice):
        imcategory, imtype = choice
        if imcategory == 'Common':
            return (imtype,)
        if imcategory == 'Composite':
            return self.fileclass.composite[imtype]
        return (imcategory,)

    def get_required_channels(self):
        required_channels = []
        for choice in self.choices:
            channels = self.find_corresponding_channel(choice)
            required_channels.extend(list(channels))
        return unique_elements(required_channels)

    def clear_unnecessary_channel(self):
        required_channels = self.get_required_channels()
        to_be_del = []
        for channel in self.raw.keys():
            if channel not in required_channels:
                to_be_del.append(channel)
        for channel in to_be_del:
            del self.raw[channel]

    def load_channel(self, channel):
        if channel not in self.raw.keys():
            self.raw[channel] = self.fileclass.extract(channel, self.georange)
        return self.raw[channel]

    def load_channel_8bit(self, channel):
        data = self.load_channel(channel)
        if channel in self.fileclass.channels['vis']:
            return np.uint8(np.clip(data, 0, 1)*255)
        if channel in self.fileclass.channels['ir']:
            return np.uint8((np.clip(data, -100, 50) + 100)  * (255/150))
        data_min, data_max = np.amin(data), np.amax(data)
        return np.uint8((data - data_min) * (255/(data_max - data_min)))

    def load_raw(self, choice):
        channel = self.find_corresponding_channel(choice)
        if len(channel) == 1:
            return self.load_channel(channel[0])
        elif self.fileclass.remap:
            return self.remap_rgb(channel)
        else:
            return self.merger(channel)

    def load_geocoord(self, choice):
        channel = self.find_corresponding_channel(choice)
        if len(channel) == 1:
            ch = channel[0]
        else:
            ch = channel[np.argmin([self.fileclass.resolution[c] for c in channel])]
        return self.fileclass.geocoord(ch)

    def load_masked_raw_geocoord(self, choice, _8bit=False):
        EDGE = 0.05
        if _8bit:
            raw = self.load_channel_8bit(choice[1])
        else:
            raw = self.load_raw(choice)
        lon, lat = self.load_geocoord(choice)
        print(np.mean(lon), np.amin(lon), np.amax(lon))
        latmin, latmax, lonmin, lonmax = self.georange
        bool_arr_ind = np.where((lat > latmin - EDGE) & (lat < latmax + EDGE) & (lon > lonmin - EDGE) & (lon < lonmax + EDGE))
        ymin, ymax = np.amin(bool_arr_ind[0]), np.amax(bool_arr_ind[0])
        xmin, xmax = np.amin(bool_arr_ind[1]), np.amax(bool_arr_ind[1])
        raw = raw[ymin:ymax, xmin:xmax]
        lat = lat[ymin:ymax, xmin:xmax]
        lon = lon[ymin:ymax, xmin:xmax]
        return raw, lon, lat

    def merger(self, channel):
        finest_resolution = min([self.fileclass.resolution[c] for c in channel])
        temp_raw = {}
        for c in channel:
            if c not in temp_raw.keys():
                if self.fileclass.resolution[c] != finest_resolution:
                    multifier = self.fileclass.resolution[c] // finest_resolution
                    temp_raw[c] = np.kron(self.load_channel_8bit(c), np.ones((multifier, multifier)))
                else:
                    temp_raw[c] = self.load_channel_8bit(c)
        return np.dstack((temp_raw[channel[0]], temp_raw[channel[1]], temp_raw[channel[2]]))

    def remap_rgb(self, channel):
        for i in range(3):
            fig = plt.figure(figsize=self.figsize)
            choice = 'Common', channel[i]
            mat, lon, lat = self.load_masked_raw_geocoord(choice)
            self.map.pcolormesh(lon, lat, mat, latlon=True, cmap=self.get_colormap(choice))
            plt.axis('off')
            plt.tight_layout(pad=0.)
            plt.savefig(join(self.fid['fdir'], 'TMP%d.png' % (i)), dpi=self.dpi)
            plt.clf()
        return np.dstack(tuple([plt.imread(join(self.fid['fdir'], 'TMP%d.png' % (i)))[:,:,0] for i in range(3)]))

    def processor(self):
        for i in range(len(self.choices)):
            ch = self.choices.pop()
            print('正在制作' + self.get_imtype(ch).upper() + '...')
            self.imager(ch)
            self.clear_unnecessary_channel()
            
    def imager(self, choice):
        fig = plt.figure(figsize=self.figsize)
        ax = plt.gca()
        if not self.fileclass.remap:
            mat = self.load_raw(choice)
            plt.imshow(mat, interpolation='nearest', extent=[self.georange[2], self.georange[3], self.georange[0], self.georange[1]],
                       **self.get_imager_kwargs(choice))
        elif choice[0] != 'Composite':
            mat, lon, lat = self.load_masked_raw_geocoord(choice)
            self.map.pcolormesh(lon, lat, mat, latlon=True, **self.get_imager_kwargs(choice))
            del lon, lat
        else:
            mat = self.load_raw(choice)
            plt.imshow(mat, interpolation='nearest', extent=[self.georange[2], self.georange[3], self.georange[0], self.georange[1]])
        if choice[0] not in ('Common', 'Composite') and self.settings['switch']['maxtbb']:
            print('(仅供参考)最高云顶温度: %.2fdegC' % (np.amax(mat)))
        del mat
        if self.settings['switch']['title']:
            plt.annotate(self.title, xy=(0.04,0.04), va='bottom', xycoords='axes fraction',
                         bbox=dict(facecolor='w', edgecolor='none', alpha=0.5), fontsize=12, family=self.title_font, weight=self.title_weight)
        plt.annotate((self.basic_info+self.get_imtype(choice)).upper(), xy=(0.5,0), va='bottom', ha='center', xycoords='axes fraction',
                     bbox=dict(facecolor='w', edgecolor='none', alpha=0.7), fontsize=6, family=self.settings['font']['info'], weight=self.settings['font']['info_weight'])
        if self.settings['switch']['coastline']:
            self.map.drawcoastlines(**self.settings['coastline'])
        if self.settings['switch']['borderline']:
            self.map.drawcountries(**self.settings['borderline'])
        if self.settings['switch']['latlon']:
            self.map.drawparallels(range(-90,90,self.settings['latlon']['step']), linewidth=self.settings['latlon']['width'],
                                   dashes=(None, None), color=self.settings['latlon']['color'])
            self.map.drawmeridians(range(0,360,self.settings['latlon']['step']), linewidth=self.settings['latlon']['width'],
                                   dashes=(None, None), color=self.settings['latlon']['color'])
        plt.axis('off')
        plt.tight_layout(pad=0.)
        plt.savefig(join(self.fid['fdir'], self.get_imtype(choice).upper() + self.fid['time'] + '.png'), dpi=self.dpi, facecolor=self.settings['image']['background'])
        plt.clf()

def launcher(path='.'):
    FS = FileSearch(path)
    Satima(FS)

launcher()
