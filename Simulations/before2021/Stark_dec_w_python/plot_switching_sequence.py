#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:45:12 2041

@author: thomas
"""
import numpy as np
import matplotlib.pyplot as plt

#This program is used to check the calculated switching sequences from
#the python code and compare it to Dominiks sequences
####-------Notes: Vertical is switched first

plt.close('all')

velo  = 425
phase = 55.47

### d_h_t = Dominik Horizontal Time, d_h_d = Dominik Horizontal Duration
# d_v_t = np.load('switching/t2jump_example/0deg_'+str(velo)+'/v_time_'+str(velo)+'.npy',allow_pickle=True)
# d_v_d = np.load('switching/t2jump_example/0deg_'+str(velo)+'/v_duration_'+str(velo)+'.npy')
# d_h_t = np.load('switching/t2jump_example/0deg_'+str(velo)+'/h_time_'+str(velo)+'.npy',allow_pickle=True)
# d_h_d = np.load('switching/t2jump_example/0deg_'+str(velo)+'/h_duration_'+str(velo)+'.npy')
#
# #
# d_v_t = np.load('switching/t2jump_example/0deg_425/0deg_time_V.npy',allow_pickle=True)
# d_v_d = np.load('switching/t2jump_example/0deg_425/0deg_duration_V.npy')
# d_h_t = np.load('switching/t2jump_example/0deg_425/0deg_time_H.npy',allow_pickle=True)
# d_h_d = np.load('switching/t2jump_example/0deg_425/0deg_duration_H.npy')
# #
#d_v_t = np.load('switching/t2jump_example/53p69deg/53p69deg_time_V.npy',allow_pickle=True)
#d_v_d = np.load('switching/t2jump_example/53p69deg/53p69deg_duration_V.npy')
#d_h_t = np.load('switching/t2jump_example/53p69deg/53p69deg_time_H.npy',allow_pickle=True)
#d_h_d = np.load('switching/t2jump_example/53p69deg/53p69deg_duration_H.npy')

d_v_t = np.load('switching/t2jump_example/55p486deg/55p486deg_time_V.npy',allow_pickle=True)
d_v_d = np.load('switching/t2jump_example/55p486deg/55p486deg_duration_V.npy')
d_h_t = np.load('switching/t2jump_example/55p486deg/55p486deg_time_H.npy',allow_pickle=True)
d_h_d = np.load('switching/t2jump_example/55p486deg/55p486deg_duration_H.npy')


####p_v_t = python vertical time
p_v_t = np.load('/Users/thomas/ownCloud/Lab/Lab315/Lab315/Software/hytrap_bg/input/input_test_py_code/v_time_'+str(velo)+'_'+str(phase)+'.npy',allow_pickle=True)
p_v_d = np.load('/Users/thomas/ownCloud/Lab/Lab315/Lab315/Software/hytrap_bg/input/input_test_py_code/v_duration_'+str(velo)+'_'+str(phase)+'.npy')
p_h_t = np.load('/Users/thomas/ownCloud/Lab/Lab315/Lab315/Software/hytrap_bg/input/input_test_py_code/h_time_'+str(velo)+'_'+str(phase)+'.npy',allow_pickle=True)
p_h_d = np.load('/Users/thomas/ownCloud/Lab/Lab315/Lab315/Software/hytrap_bg/input/input_test_py_code/h_duration_'+str(velo)+'_'+str(phase)+'.npy')

d_horizontal = []
d_vertical   = []
p_horizontal = []
p_vertical   = []
for i in range(len(d_v_d)):
    d_v_t[i] = np.round(d_v_t[i],9)
    d_v_d[i] = np.round(d_v_d[i],9)
    d_vertical = np.append(d_vertical,d_v_t[i])
    d_vertical = np.append(d_vertical,d_v_t[i]+d_v_d[i])

    p_v_t[i] = np.round(p_v_t[i],9)
    p_v_d[i] = np.round(p_v_d[i],9)
    p_vertical = np.append(p_vertical,p_v_t[i])
    p_vertical = np.append(p_vertical,p_v_t[i]+p_v_d[i])

for i in range(len(d_h_d)):
    d_h_t[i] = np.round(d_h_t[i],9)
    d_h_d[i] = np.round(d_h_d[i],9)
    d_horizontal = np.append(d_horizontal,d_h_t[i])
    d_horizontal = np.append(d_horizontal,d_h_t[i]+d_h_d[i])

    p_h_t[i] = np.round(p_h_t[i],9)
    p_h_d[i] = np.round(p_h_d[i],9)
    p_horizontal = np.append(p_horizontal,p_h_t[i])
    p_horizontal = np.append(p_horizontal,p_h_t[i]+p_h_d[i])

dominik_collected_time = []
dominik_collected_amp = []

python_collected_time = []
python_collected_amp = []

####---- Getting the times for the plot... super easy...
####---- Getting the times for the plot... super easy...
def get_times(time,amp,v,h):
    for i in range(61):
        time = np.append(time,v[i*2])
        amp = np.append(amp,-1)
        time = np.append(time,v[(i*2)+1])
        amp = np.append(amp,-1)

        time = np.append(time,h[i*2])
        amp = np.append(amp,1)
        time = np.append(time,h[(i*2)+1])
        amp = np.append(amp,1)

    time = np.append(time,v[61*2])
    amp = np.append(amp,-1)
    time = np.append(time,v[(61*2)+1])
    amp = np.append(amp,-1)
    return(time,amp)

dominik_collected_time,dominik_collected_amp = get_times(dominik_collected_time,dominik_collected_amp,d_vertical,d_horizontal)
python_collected_time,python_collected_amp = get_times(python_collected_time,python_collected_amp,p_vertical,p_horizontal)




plt.figure()
plt.title(str(velo)+' m/s')
plt.plot(dominik_collected_time,dominik_collected_amp,label='Fortran')
plt.plot(python_collected_time,python_collected_amp,label='Python',alpha=0.5)
plt.legend()
plt.show()
