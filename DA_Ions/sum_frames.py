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
import matplotlib.dates as mdates
import time
import datetime
from lmfit import Model

from scipy.ndimage import gaussian_filter1d
plt.close('all')

#/Volumes/sw_group/pc-swplab-pc08/2020-08-06/2020-08-06_001.h5
###File to analyze
loc     =["/Volumes/sw_group/pc-swplab-pc08/"]
date     = "2020-05-11"
file     = ['007']




for i in range(len(loc)):
    try:
        link = loc[i]+date+'/'+date+'_'+file[0]+'.h5'
        h5_data = h5py.File(link,'r')
        print('Analyzing file: '+link)
    except:
        print('file reading failed')
        pass

nbr_frames = len(list(h5_data.keys()))
print('Number of frames:\t\t', nbr_frames)
frames = []
for i in range(0, nbr_frames):
    frames = np.append(frames,'data'+str((i)))

sum_frames = np.zeros(np.shape(h5_data[frames[0]]))

for i in range(len(frames)):
    sum_frames = sum_frames + h5_data[frames[0]][:]

print(np.sum(sum_frames))
sum_frames = sum_frames/len(frames)
print(np.sum(sum_frames))

#for presentation purposes
sum_frames[226:228, 166] = 500


plt.figure()
plt.imshow(sum_frames.T, origin='lower',cmap='inferno', vmin=200, vmax=1800)
plt.colorbar()


#print nice image
resize_im = sum_frames[0:499, 100:400].T #change to ur taste
resize_im[ resize_im > 1800] = 500
fig = plt.figure(frameon=False)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)

ax.imshow(resize_im, cmap='inferno', aspect='auto')
plt.savefig('print_ion_temp.png', format='png', cmap='inferno')



plt.show()
