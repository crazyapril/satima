from configparser import ConfigParser

class Settings:

    dpi = 400
    widthsize = 5.
    
    def __init__(self, fname='Satima.ini'):
        self.settings = {'switch':{}, 'image':{}, 'font':{}}
        self.config = ConfigParser()
        self.config.read(fname)
        self.parse_switch_settings()
        self.parse_other_settings()

    def parse_switch_settings(self):
        for attr in self.config.options('switch'):
            self.settings['switch'][attr] = self.config.getboolean('switch', attr)
            if self.settings['switch'][attr] and attr not in ('title', 'maxtbb'):
                self.settings[attr] = {}
                for attribute in self.config.options(attr):
                    if attribute == 'color': self.settings[attr][attribute] = self.config.get(attr, attribute)
                    else: self.settings[attr][attribute] = eval(self.config.get(attr, attribute))

    def parse_other_settings(self):
        for attr in self.config.options('image'):
            self.settings['image'][attr] = self.config.get('image', attr)
        for attr in self.config.options('font'):
            self.settings['font'][attr] = self.config.get('font', attr)
