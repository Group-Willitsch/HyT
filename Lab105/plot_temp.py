# -*- coding: utf-8 -*-
"""
@author: T to da Kierspel !!!!!
"""

##############
## Script listens to serial port and writes contents into a file
##############
import numpy as np
from matplotlib import pyplot as plt
import sys
import os
from pathlib import Path
from scipy.ndimage.filters import gaussian_filter1d

####---- General preparation of the script
plt.close('all')
computer = 0       ###Makes it easier later with the different computer


####---- Get the right folder for the stored data on the working computer
try:
    if os.environ['HOME']=='/Users/thomas':
        data_folder = '/Users/thomas/ownCloud/Lab/Data/'
        computer = 0
except:
    home = str(Path.home())
    if home == 'C:\\Users\\Gruppe Willitsch':
        data_folder = 'C:\\Users\\Gruppe Willitsch\\Desktop\\data\\temp_log\\'
        computer = 1
    else:
        print('Unknown directory')
        sys.exit()

####---- Which days do you want to analyse
days = ['2020-02-10','2020-02-11']

####---- Collect the data
for i in range(len(days)):
    if computer == 0:
        link = data_folder+days[i]+'/'+'temp_'+days[i]+'.txt'
    if computer == 1:
        link = data_folder+'temp_'+days[i]+'.txt'
        
    if i == 0:
        data = np.genfromtxt(link,skip_header=1)
    else:
        data = np.append(data,np.genfromtxt(link,skip_header=1),axis=0)
    


####---- Plot the data
fgc = 0             ###Figure Counter
plt.figure(fgc)
#plt.plot(np.convolve(data[:,4],1000))
plt.plot(gaussian_filter1d(data[:,4],2))



plt.show()
