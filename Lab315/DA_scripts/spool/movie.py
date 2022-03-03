# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:35:41 2020
@author: T Kierspel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
import sys
import h5py

from scipy.ndimage import gaussian_filter1d

plt.close('all')    

location = "C:/Users/Gruppe Willitsch/Desktop/data/2020-09-16/"
file     = ['2020-09-16_001.h5']

data = location + file[0]


h5_data = h5py.File(data,'r')

frames = int(len(list(h5_data.keys())))


data_mean = np.zeros(np.shape(h5_data['data0'][:]))

data_norm_up = np.zeros(frames)
data_norm_mid = np.zeros(frames)
data_norm_low = np.zeros(frames)

for i in range(frames):
    data_new = h5_data['data%i'%i][:]
    data_mean += data_new

    norm_fac  = np.mean(data_new[800:,800])
    data_new  = data_new/norm_fac

    data_norm_up[i]  = np.mean(data_new[:,750:])
    data_norm_mid[i] = np.mean(data_new[:,250:750])
    data_norm_low[i] = np.mean(data_new[:,0:250])

data_mean = data_mean/frames




plt.figure(0)
plt.imshow(data_mean.T,origin='lower',interpolation='none',cmap='gnuplot',vmax=1000)

plt.plot(np.mean(data_mean,axis=0))

plt.figure(1)
plt.plot(data_norm_up)
plt.plot(data_norm_mid)
plt.plot(data_norm_low)

plt.figure(2)
plt.plot(data_norm_up/np.max(data_norm_up))
plt.plot(data_norm_mid/np.max(data_norm_mid))
plt.plot(data_norm_low/np.max(data_norm_low))

plt.figure(3)
plt.plot(gaussian_filter1d(data_norm_up,2))
plt.plot(gaussian_filter1d(data_norm_mid,2))
plt.plot(gaussian_filter1d(data_norm_low,2))


plt.show()
