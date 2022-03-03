# readfile
import pandas as pd
import numpy as np


def readpars(filename):
    with open(filename) as f:
        lines = f.readlines()
    params = []
    values = []

    for i in range(len(lines) - 1):
        param = ''
        value = ''
        j = 0
        val = 0
        while True:
            try:
                char = lines[i][j]
            except IndexError:
                break

            if (char != '=') & val == 0:
                param += char

            if val == 1:
                value += char

            if char == '=':
                val = 1

            if char == '#':
                break

            j += 1

        param = param.rstrip('=')
        value = value.rstrip('\n')
        param = param.strip()
        value.replace(" ", "")
        param.replace(" ", "")
        if param != '#':
            params.append(param)
            values.append(value)

    return params, values


def readpars_class(filename):
    with open(filename) as f:
        lines = f.readlines()
    params = []
    values = []

    for i in range(len(lines) - 1):
        param = ''
        value = ''
        j = 0
        val = 0
        while True:
            try:
                char = lines[i][j]
            except IndexError:
                break

            if (char != '=') & val == 0:
                param += char

            if val == 1:
                value += char

            if char == '=':
                val = 1

            if char == '#':
                break

            j += 1

        param = param.rstrip('=')
        value = value.rstrip('\n')
        param = param.strip()
        value.replace(" ", "")
        value.strip(' ')
        param.replace(" ", "")
        if param != '#':
            params.append('self.' + param)
            values.append(value)

    return params, values


def readdecfile(filename, sign):
    '''
            Reads T2jump.out files from deceleration simulation
            required for creating decelerator burst sequences!
    '''

    rowsskip = (0, 1, 2, 3)
    data = pd.read_csv(filename, delim_whitespace=True, skiprows=rowsskip, error_bad_lines=False, header=None,
                       skipfooter=3, engine='python')
    data = np.array(data)
    # data=np.insert(data,0,[1010,'0x0010','!',0],axis=0)
    #    print(data)
    #    print(data.T[1][data.T[1] == '0x0010'].shape)
    #    print(data.T[1][data.T[1] == '0x0020'].shape)

    '''
        since this line (first of decelerator sequence) is always the same and it does not have
        the same shape as the rest, it is skipped when reading the file and then added.
    '''

    duration = []
    # print(data)
    # data=pd.DataFrame(data)
    ns = 1e-9
    unit = ns

    if sign == 'H+':
        time = np.array([data[i][0] - data[0][0] for i in range(len(data)) if data[i][1][1] == '1'])
        duration = np.array([(data[i + 1][0] - data[i][0]) for i in range(len(data) - 1) if data[i][1][1] == '1'])
        # time=np.array(data[data[1]=='0x0010'][0]-data[0][0])
        # duration=np.array([(data[0][2*i+1]-data[0][2*i]) for i in range(time.shape[0])])
        # print(time)
        # print(duration)

        time = np.array(time) * unit
        duration = np.array(duration) * unit

    if sign == 'H-':
        time = np.array([data[i][0] - data[0][0] for i in range(len(data)) if data[i][1][2] == '1'])
        duration = np.array([(data[i + 1][0] - data[i][0]) for i in range(len(data) - 1) if data[i][1][2] == '1'])
        # print(time)
        # print(duration)

        time = np.array(time) * unit
        duration = np.array(duration) * unit

    if sign == 'V+':
        time = np.array([data[i][0] - data[0][0] for i in range(len(data)) if data[i][1][3] == '1'])
        duration = np.array([(data[i + 1][0] - data[i][0]) for i in range(len(data) - 1) if data[i][1][3] == '1'])
        # print(time)
        # print(duration)

        time = np.array(time) * unit
        duration = np.array(duration) * unit

    if sign == 'V-':
        time = np.array([data[i][0] - data[0][0] for i in range(len(data)) if data[i][1][4] == '1'])
        duration = np.array([(data[i + 1][0] - data[i][0]) for i in range(len(data) - 1) if data[i][1][4] == '1'])
        # print(time)
        # print(duration)

        time = np.array(time) * unit
        duration = np.array(duration) * unit
    # if sign=='V+':
    #     time=(data[data[1]=='0x0020'][0]-data[0][0])
    #     duration=np.array([data[0][2*i+2]-data[0][2*i+1] for i in range(time.shape[0])])
    #     time=np.array(time)*unit
    #     duration=np.array(duration)*unit

    return time, duration


# ~ for i in range (len(params)):
# ~ try:
# ~ exec(params[i]+"=%s"%(values[i]))
# ~ except (NameError, SyntaxError):
# ~ exec(params[i]+"='%s'"%(values[i]))

# ~ print(scan_step)

if __name__ == '__main__':
    #    t2jumpfile = './input/sequences_FM/T2jump_fm.out'
    t2jumpfile = r'.\input\sequencesFM_vel_12p5kV_124switches\dec_450_30\T2jump.out'
    dectime, decdur = readdecfile(t2jumpfile, 'V+')
#    print('dectime:', dectime)
#    print('decdur:', decdur)
    # filename = '/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/53p69deg/T2jump.dat'
    # params, values = readpars(filename)
    # time_H, duration_H = readdecfile('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/53p69deg/T2jump.out','H')
    # time_V, duration_V = readdecfile('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/53p69deg/T2jump.out','V')

    # save = False
    # if save:
    #     np.save('53p69deg_time_H',time_H)
    #     np.save('53p69deg_duration_H',duration_H)
    #     np.save('53p69deg_time_V',time_V)
    #     np.save('53p69deg_duration_V',duration_V)
