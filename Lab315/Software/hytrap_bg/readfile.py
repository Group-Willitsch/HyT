#readfile
import pandas as pd
import numpy as np

def readpars(filename):
    with open(filename) as f:
        lines = f.readlines()
    params=[]
    values=[]

    for i in range (len(lines)-1):
        param=''
        value=''
        j=0
        val=0
        while True:
            try:
                char=lines[i][j]
            except IndexError:
                break

            if (char is not '=')&val==0:
                param+=char

            if val==1:
                value+=char

            if char=='=':
                val=1

            if char=='#':
                break

            j+=1

        param = param.rstrip('=')
        value = value.rstrip('\n')
        param = param.strip()
        value.replace(" ","")
        param.replace(" ","")
        if param is not '#':
            params.append(param)
            values.append(value)

    return params, values

def readpars_class(filename):
    with open(filename) as f:
        lines = f.readlines()
    params=[]
    values=[]

    for i in range (len(lines)-1):
        param=''
        value=''
        j=0
        val=0
        while True:
            try:
                char=lines[i][j]
            except IndexError:
                break

            if (char is not '=')&val==0:
                param+=char

            if val==1:
                value+=char

            if char=='=':
                val=1

            if char=='#':
                break

            j+=1

        param = param.rstrip('=')
        value = value.rstrip('\n')
        param = param.strip()
        value.replace(" ","")
        value.strip(' ')
        param.replace(" ","")
        if param is not '#':
            params.append('self.'+param)
            values.append(value)

    return params, values






def readdecfile(filename,sign):
    '''
            Reads T2jump.out files from deceleration simulation
            required for creating decelerator burst sequences!
    '''

    rowsskip= (0,1,2,3,4,128,129,130,131) #skip unnecessary rows, s.t. files can be read in easily
    data=pd.read_csv(filename,delim_whitespace=True,skiprows = rowsskip,error_bad_lines=False, header = None)
    data=np.array(data)
    data=np.insert(data,0,[1010,'0x0010','!',0],axis=0)
#    print(data)
#    print(data.T[1][data.T[1] == '0x0010'].shape)
#    print(data.T[1][data.T[1] == '0x0020'].shape)

    '''
        since this line (first of decelerator sequence) is always the same and it does not have
        the same shape as the rest, it is skipped when reading the file and then added.
    '''

    duration=[]
    data=pd.DataFrame(data)
    ns=1e-9
    unit=ns
    if sign=='V':
        time=np.array(data[data[1]=='0x0010'][0]-data[0][0])
        duration=np.array([(data[0][2*i+1]-data[0][2*i]) for i in range(time.shape[0])])
        time=np.array(time)*unit
        duration=np.array(duration)*unit

    if sign=='H':
        time=(data[data[1]=='0x0020'][0]-data[0][0])
        duration=np.array([data[0][2*i+2]-data[0][2*i+1] for i in range(time.shape[0])])
        time=np.array(time)*unit
        duration=np.array(duration)*unit

    return time, duration

#~ for i in range (len(params)):
    #~ try:
        #~ exec(params[i]+"=%s"%(values[i]))
    #~ except (NameError, SyntaxError):
        #~ exec(params[i]+"='%s'"%(values[i]))

#~ print(scan_step)

if __name__ == '__main__':
    filename = '/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/53p69deg/T2jump.dat'
    params, values = readpars(filename)
    time_H, duration_H = readdecfile('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/53p69deg/T2jump.out','H')
    time_V, duration_V = readdecfile('/Users/thomas/ownCloud/Lab/Lab315/Simulations/Stark_dec/switching/t2jump_example/53p69deg/T2jump.out','V')

    save = True
    if save:
        np.save('53p69deg_time_H',time_H)
        np.save('53p69deg_duration_H',duration_H)
        np.save('53p69deg_time_V',time_V)
        np.save('53p69deg_duration_V',duration_V)
