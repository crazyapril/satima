from configparser import ConfigParser
from matplotlib.colors import LinearSegmentedColormap as LSCMAP

def parse_colormap_line(line):
    lineinfo = {'temp':255, 'r':[], 'b':[], 'g':[], 'status':0}
    if line[-1] == '\n':
        line = line[:-1]
    data = line.split(' ')
    datalen = len(data)
    if datalen not in (4, 5, 7):
        print('Colormap Error Code: 1\nError Line:'+line)
        return lineinfo
    if datalen == 5:
        if data[-1] == '~':
            data = data[:-1]
            data.extend(data[1:])
        else:
            print('Colormap Error Code: 2\nError Line:'+line)
            return lineinfo
    try:
        intdata = [float(d) for d in data]
    except (SyntaxError, ValueError, NameError):
        print('Colormap Error Code: 3\nError Line:'+line)
        return lineinfo
    lineinfo['temp'] = intdata[0]
    for i, num in enumerate(intdata[1:]):
        if num < 0 or num > 255:
            print('Colormap Error Code: 4\nError Line:')
            return lineinfo
        if i % 3 == 0:
            lineinfo['r'].append(num)
        elif i % 3 == 1:
            lineinfo['g'].append(num)
        else:
            lineinfo['b'].append(num)
    if datalen == 4:
        lineinfo['status'] = 1
    else:
        lineinfo['status'] = 2
    return lineinfo
    
def parse_colormap_data(file):
    data = {'temp':[], 'r':[], 'g':[], 'b':[]}
    with open(file, 'r') as f:
        line = f.readline()
        lineinfo = parse_colormap_line(line)
        if lineinfo['status'] != 1:
            return 0
        pasttemp = lineinfo['temp']
        data['temp'].append(lineinfo['temp'])
        for color in 'rgb':
            data[color].append(lineinfo[color])
        while True:
            line = f.readline()
            lineinfo = parse_colormap_line(line)
            if lineinfo['status'] == 0:
                return 0
            nowtemp = lineinfo['temp']
            if nowtemp <= pasttemp:
                print('Colormap Error Code: 4\nError Temp:'+nowtemp+' <= '+pasttemp)
                return 0
            pasttemp = nowtemp
            data['temp'].append(lineinfo['temp'])
            for color in 'rgb':
                data[color].append(lineinfo[color])
            if lineinfo['status'] == 1:
                if f.read(1) != '':
                    return 0
                break
    return data

def parse_colormap(data):
    colormap = {'red':[], 'green':[], 'blue':[]}
    template = (0, 0, 0)
    vmin, vmax = data['temp'][0], data['temp'][-1]
    span = vmax - vmin
    colormap['red'].append((0., 0., data['r'][0][0]/255))
    colormap['green'].append((0., 0., data['g'][0][0]/255))
    colormap['blue'].append((0., 0., data['b'][0][0]/255))
    for t, r, g, b in zip(data['temp'][1:-1], data['r'][1:-1], data['g'][1:-1], data['b'][1:-1]):
        rel_t = (t - vmin) / span
        colormap['red'].append((rel_t, r[0]/255, r[1]/255))
        colormap['green'].append((rel_t, g[0]/255, g[1]/255))
        colormap['blue'].append((rel_t, b[0]/255, b[1]/255))
    colormap['red'].append((1., data['r'][-1][0]/255, 0.))
    colormap['green'].append((1., data['g'][-1][0]/255, 0.))
    colormap['blue'].append((1., data['b'][-1][0]/255, 0.))
    return vmin, vmax, colormap

def get_colormap_list():
    __colormaps__ = ConfigParser()
    __colormaps__.read('colormap/colormaps')
    return __colormaps__.get('wv', 'wv').upper().replace(' ', '').split(','), __colormaps__.get('ir', 'ir').upper().replace(' ', '').split(',')

def get_colormap(name):
    data = parse_colormap_data('colormap/'+name.lower() + '.txt')
    if not data:
        return 0
    vmin, vmax, colormap = parse_colormap(data)
    return LSCMAP(name, colormap)
