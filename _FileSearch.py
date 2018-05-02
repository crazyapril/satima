from glob import glob
from os.path import isfile, split, join
import pkgutil
from tkinter import Tk
from tkinter.filedialog import askdirectory
from _confs import convert_time_format, get_single_choice_template

class FileSearch:

    def __init__(self, filedir='.'):
        self.options = []
        self.num_of_options = 0
        self.choice = None
        self.filedir = filedir
        self.find_available_fileclass()
        self.run()

    def run(self):
        self.searcher()
        self.display_options()
        self.get_choice()

    def searcher(self):
        for fileclass in self.fileclass.values():
            self.options.extend(fileclass.search(self.filedir))
        self.num_of_options = len(self.options)

    def get_choice(self):
        results = get_single_choice_template(list(range(self.num_of_options + 2)), '输入要读取的数据编号：')
        if results == self.num_of_options + 1:
            is_filedir_changed = self.change_filedir()
            if is_filedir_changed: self.run()
            else: self.get_choice()
        elif results == 0:
            exit()
        else:
            self.choice = self.options[results - 1]
        
    def display_options(self):
        if self.filedir != '.':
            print('\n目前目录：{}'.format(self.filedir))
        for i, fid in enumerate(self.options, 1):
            print('{:d}. [{}] {} {}'.format(i, fid['sate'], convert_time_format(fid['time']), fid['area']))
        print('{:d}. 选择文件夹...'.format(self.num_of_options + 1))

    def change_filedir(self):
        root = Tk()
        root.withdraw()
        tmpdir = askdirectory(initialdir=self.filedir, title='选择要打开的文件夹')
        if tmpdir != '':
            self.filedir = tmpdir
            return True
        else:
            return False

    def find_available_fileclass(self):
        import satellite
        self.fileclass = {}
        prefix = satellite.__name__ + '.'
        for importer, modname, ispkg in pkgutil.iter_modules(satellite.__path__, prefix):
            module = __import__(modname, fromlist='dummy')
            mname = modname.split('.')[-1]
            if module.READABLE:
                self.fileclass[mname] = getattr(module, mname)
