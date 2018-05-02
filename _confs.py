#Convert time format from YYYYMMDDhhmm to YYYY/MM/DD hhmmZ
def convert_time_format(stdtime):
    return '{}/{}/{} {}Z'.format(stdtime[:4], stdtime[4:6], stdtime[6:8], stdtime[8:])

def get_single_choice_template(options, prompt):
    while True:
        choice = int(input(prompt))
        if choice in options:
            return choice
    
def get_multi_choices_template(options, prompt, sort=False):
    choices = []
    while len(choices) == 0:
        raw_choices = input(prompt).split()
        for raw_ch in raw_choices:
            ch = int(raw_ch)
            if ch in options:
                if ch not in choices:
                    choices.append(ch)
            else:
                print('{}不存在，未录入。'.format(raw_ch))
        if len(choices) == 0:
            print('未录入任何信息，请重新输入。\n')
    if sort:
        choices.sort()
    return choices

def safe_str2value(instr, func=eval):
    try:
        return func(instr)
    except (NameError, SyntaxError, ValueError):
        print('{:s}不是有效的数字。'.format(instr))
        return None

def safe_input(prompt, func=eval, vrange=None, default=None):
    while True:
        strval = input(prompt)
        if len(strval) == 0 and default is not None:
            return default
        value = safe_str2value(strval, func)
        if value is not None:
            if vrange is None:
                return value
            vmin, vmax = vrange
            if value >= vmin and value <= vmax:
                return value
            else:
                print('数字不在有效范围内。')

def unique_elements(iterable):
    unique_list = []
    iterable = list(iterable)
    for i in iterable:
        if i not in unique_list:
            unique_list.append(i)
    return unique_list

def is_ascii(s):
    return all(ord(c) < 128 for c in s)
