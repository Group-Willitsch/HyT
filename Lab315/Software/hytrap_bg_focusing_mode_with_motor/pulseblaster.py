import spinapi as sp
import time
import signal
import sys
import pandas as pd
import numpy as np
import readfile as rf


def decfile(directory,filename):
    '''
            Reads T2jump.out files from deceleration simulation
            required for creating decelerator burst sequences!
    '''
    
#    rowsskip= (0,1,2,3,4,129,130,131) #skip unnecessary rows, s.t. files can be read in easily
    rowsskip= (0,1,2,3,4,130,131) #skip unnecessary rows, s.t. files can be read in easily

    data=pd.read_csv(directory+filename,skiprows=rowsskip,delim_whitespace=True,header=None,error_bad_lines=False)
    data=np.array(data)
    data=np.insert(data,0,[1010,'0x0010','!',0],axis=0) 
    ''' 
        since this line (first of decelerator sequence) is alwasy the same and it does not have
        the same shape as the rest, it is skipped when reading the file and added later.
    '''
    return data.T
    
def channelflag(numbers):
    '''
        takes a list of channelnumbers as input. Creates a string of 24 zeroes 
        (1 for each available channel on pulseblaster) and then sets the string of each channel
        number from the input to 1. Required for Spinapi instruction functions.
        Example:
            use for channels [1,3]
            output will be
            '101000000000000000000000'
    '''
    maxno=24
    flags=list('0'*maxno)
    output=''
    for num in numbers:
        flags[num-1]='1'
    
    for flag in flags:
        output+=flag
        
    if numbers==[0]:
        output='0'*maxno
    return output

def pulseprog():
    ''''testprogram for pulseblaster'''
    directory = './input/31102017-470-470-0deg-431us-10kv/'
    filename='T2jump.out'
    dec_data=decfile(directory,filename)
    #~ print(dec_data)
    
    clock=100 #MHz
    dec_chn = {'0x0010':6, '0x0020':7, '0x0000':0} 
    
    '''decelerator channels  in t2jump.out 
                0x0010 for 6, dec(+), 
                0x0020 for 7, dec(-)'''
        
    if (sp.pb_init() != 0):
        print('Error initializing board: %s\n'%sp.pb_get_error())
    print('spinpts version: %s\n'%sp.spinpts_get_version())
    print('firmware id: %s\n'%sp.pb_get_firmware_id())
    print('board count: %s\n'%sp.pb_count_boards())
    #~ sp.pb_set_debug(1)
    
    sp.pb_core_clock(clock)
    sp.pb_select_board(0)
    sp.pb_start_programming(sp.PULSE_PROGRAM)
    
    #~ decelerator
    pulse=0
    reprate=10 #Hz
    trep=1/reprate *sp.s
    tsequence=0
    print('reprate: %f, time rep in ms: %f'%(reprate,trep*sp.ms))

    #first pulse
    timeunit=sp.ns
    duration=(dec_data[0][pulse+1]-dec_data[0][pulse])*timeunit
    channels=channelflag([dec_chn[dec_data[1][pulse]]])
    start = sp.pb_inst_pbonly(channels,sp.CONTINUE,0,duration)
    tsequence+=duration
    
    for pulse in range(1,len(dec_data[0])-1):
        channels=channelflag([dec_chn[dec_data[1][pulse]]])
        duration=(dec_data[0][pulse+1]-dec_data[0][pulse])*timeunit
        sp.pb_inst_pbonly(channels,sp.CONTINUE,0,duration)
        tsequence+=duration
        
    #~ last pulse
    pulse=len(dec_data[0])-1

    duration=trep-tsequence
    print('waitingtime for next pulse to start in ms: %f'%(duration/sp.ms))
    channels=channelflag([dec_chn[dec_data[1][pulse]]])    
    stop = sp.pb_inst_pbonly(0x0,sp.BRANCH,start,duration)
    
    sp.pb_stop_programming()
    print(sp.pb_read_status())
    
def pulseblaster_program(inputs):
    timeunit=sp.s
    
    try:
        if (sp.pb_init() != 0):
            output = 'Error initializing board: %s\n'%sp.pb_get_error()
            print(output)
            return output
    except:
        print(str(sys.exc_info())+'\n sp.pb_init() failed')
#    print('spinpts version: %s\n'%sp.spinpts_get_version())
#    print('firmware id: %s\n'%sp.pb_get_firmware_id())
#    print('board count: %s\n'%sp.pb_count_boards())
    #~ sp.pb_set_debug(1)
    clock = 100.0
    sp.pb_core_clock(clock)
    sp.pb_select_board(0)
    sp.pb_start_programming(sp.PULSE_PROGRAM)
    sumz = 0
#    print(inputs)
    start = sp.pb_inst_pbonly(inputs[0][0],sp.CONTINUE,0,float(inputs[1][0])*timeunit)
    sumz += float(inputs[1][0])*timeunit
#    print(inputs[1][-1],inputs[0][-1])
    for i in range(1,len(inputs[0])-1):
        try:
            sp.pb_inst_pbonly(inputs[0][i],sp.CONTINUE,0,round(float(inputs[1][i])*timeunit,8))
            sumz += float(inputs[1][i])*timeunit
        except:
            print(str(sys.exc_info()))
            print('problem at: ',inputs[0][i],inputs[1][i],'\n')
            print(round(float(inputs[1][i])*timeunit,8))
    sumz += float(inputs[1][-1])*timeunit
#    print(sumz)
    stop = sp.pb_inst_pbonly(inputs[0][-1],sp.BRANCH,start,float(inputs[1][-1])*timeunit)
    sp.pb_stop_programming()

    return sp.pb_read_status()
    
def pulseblaster_start():
    # try:
    sp.pb_start()
    # except Exception as err:
    #     raise Exception(err.args[0])
    return

def pulseblaster_stop():
    # try:
    sp.pb_stop()
    sp.pb_close()
    # except Exception as err:
    #     raise Exception(err.args[0])
    return

def pulseblaster_status():
    return sp.pb_read_status()
    
if __name__ == "__main__":
    if (sp.pb_init() != 0):
        print('Error initializing board: %s\n'%sp.pb_get_error())
    print('spinpts version: %s\n'%sp.spinpts_get_version())
    print('firmware id: %s\n'%sp.pb_get_firmware_id())
    print('board count: %s\n'%sp.pb_count_boards())
#    sp.pb_core_clock(100.0)
#
#    sp.pb_select_board(0)
#    sp.pb_start_programming(sp.PULSE_PROGRAM)
#    trigtime = 100*sp.us
#    start = sp.pb_inst_pbonly('10000',sp.CONTINUE,0,trigtime)
#    stop = sp.pb_inst_pbonly('00000',sp.BRANCH,start,(100*sp.ms-trigtime))
#    sp.pb_stop_programming()
#
##    pulseprog()
#    sp.pb_start()
#    print('press enter to stop')
#    input()
#    sp.pb_stop()
    sp.pb_close()


