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

###File to analyze
loc     =["C:/Users/Gruppe Willitsch/Desktop/data/","/Users/thomas/ownCloud/Lab/Data/"]
date     = "2020-07-03"
file     = ['004']


for i in range(len(loc)):
    try:
        link = loc[i]+date+'/'+date+'_'+file[0]+'.h5'
        h5_data = h5py.File(link,'r')
        print('Analyzing file: '+link)
    except:
        pass

nbr_frames = len(list(h5_data.keys()))
frames = []
for i in range(0,nbr_frames):
    frames = np.append(frames,'data'+str((i)))


###Relevant variables
pixels = np.shape(h5_data[frames[0]][:])
roi = np.zeros(pixels)
frame_intensities = []
calibration = True
timestamp = []

####--Find the conncted frames
higher = True
intensity_threshold = 520
intervals = 0

###Relevant variables after calibration
region_oi = [110,380,195,295]
frames_to_anal = [0,7,82,125,201,249,328,409,489,546,652,693]


###Loop over all files --- Only requried to set the correct areas for DA and calirate stuff
if calibration:
    for i in range(0,nbr_frames):
        current_frame = h5_data[frames[i]][:]
        roi = roi+current_frame
        frame_intensities = np.append(frame_intensities,np.mean(current_frame[110:380,195:295]))
        if frame_intensities[i]<intensity_threshold and higher:
            intervals = np.append(intervals,i)
            higher = False

        if frame_intensities[i]>intensity_threshold and higher == False:
            intervals = np.append(intervals,i)
            higher = True

    ###Find interesting frames
    if len(intervals)%2 == 1:
        print('yellow sub')
        intervals = np.append(intervals,nbr_frames-1)

    print(intervals)
    frames_to_anal = intervals

    ###Normalize
    roi = roi/nbr_frames



    plt.figure(0)
    plt.imshow(roi[110:380,195:295].T,origin='lower')

    plt.figure(1)
    plt.plot(frame_intensities)
    plt.grid()
    for i in range(len(intervals)):
        if i%2 == 0:
            color = 'red'
        if i%2 == 1:
            color = 'blue'
        plt.plot((intervals[i],intervals[i]),(500,600),color=color)
    #plt.show()
    #sys.exit()

###--- Get data from ROI
data_array = h5_data[frames[0]][region_oi[0]:region_oi[1],region_oi[2]:region_oi[3]]

for i in range(1,nbr_frames):
    current_frame = h5_data[frames[i]][region_oi[0]:region_oi[1],region_oi[2]:region_oi[3]]
    data_array = np.dstack((data_array,current_frame))

    ###Get corresponding attibutes
    ###NOTE--- Since loop starts with i = 1 I have to put some corrections in...:()
    if i == 1:
            attributes = h5_data[frames[i-1]].attrs
            time_att = attributes['timestamp']
            time_att = time_att.astype(str) ###---Convert to string
            time_in_sec = time.mktime(datetime.datetime.strptime(time_att,"%Y-%m-%d %H:%M:%S").timetuple())  ###---Convert string to seconcs
            timestamp = np.append(timestamp,time_in_sec)

    attributes = h5_data[frames[i]].attrs
    ###Get time in seconds
    time_att = attributes['timestamp']
    time_att = time_att.astype(str) ###---Convert to string
    time_in_sec = time.mktime(datetime.datetime.strptime(time_att,"%Y-%m-%d %H:%M:%S").timetuple())  ###---Convert string to seconcs
    timestamp = np.append(timestamp,time_in_sec)



# ###--- Extract individual crystalls
mean_to_analyse = np.zeros(np.shape(current_frame))
timestamp_falling_edge = [0]
timestamp_rising_edge = [0]
timestamp = timestamp-timestamp[0]
projection_horizontal = np.zeros(np.shape(mean_to_analyse[:,0]))
projection_vertical   = np.zeros(np.shape(mean_to_analyse[0]))


for i in range(int(len(frames_to_anal)/2)):
    to_analyse = np.mean(data_array[:,:,frames_to_anal[i*2]:frames_to_anal[(i*2+1)]],axis=2)
    mean_to_analyse = np.dstack((mean_to_analyse,to_analyse))

    ###Get the time
    timestamp_falling_edge = np.append(timestamp_falling_edge,timestamp[frames_to_anal[i*2+1]])
    timestamp_rising_edge = np.append(timestamp_rising_edge,timestamp[frames_to_anal[i*2]])

    projection_vertical = np.vstack((projection_vertical,np.mean(to_analyse,axis=0)))
    projection_horizontal = np.vstack((projection_horizontal,np.mean(to_analyse,axis=1)))

###---Get the dark time
t_dark_delta = 0
t_dark_tot   = 0
for i in range(len(timestamp_falling_edge)-2):
    t_dark_delta = np.append(t_dark_delta,timestamp_rising_edge[i+2]-timestamp_falling_edge[i+1])
    t_dark_tot = np.append(t_dark_tot,np.sum(t_dark_delta))


#####Remove dummy image
mean_to_analyse=mean_to_analyse[:,:,1:]
projection_vertical = projection_vertical[1:,:]
projection_horizontal = projection_horizontal[1:,:]
#####Remove the last entry
timestamp_falling_edge = timestamp_falling_edge[:-1]


###---- Offset correction
projection_vertical = projection_vertical-np.mean(projection_vertical[:,90:])
projection_horizontal = projection_horizontal-np.mean(projection_horizontal[:,255:])
###---- Normalization
for i in range(int(len(frames_to_anal)/2)):
    projection_vertical[i,:]   = projection_vertical[i,:]/np.mean( projection_vertical[i,:])
    projection_horizontal[i,:] = projection_horizontal[i,:]/np.mean(projection_horizontal[i,:])

colors = ['black','red','green','blue','orange','yellow','gray']
colors = np.append(colors,colors)
colors = np.append(colors,colors)


###--- Show images images
plt.figure()
plt.title('Integrated Signal')
plt.plot(timestamp,np.mean(np.mean(data_array,axis=0),axis=0))
plt.xlabel('')

for i in range(int(len(frames_to_anal)/2)):
    plt.figure()
    plt.title(str(i))
    plt.imshow(mean_to_analyse[:,:,i].T,origin='lower',vmax=1000,vmin=500)

for i in range(int(len(frames_to_anal)/2)):
    plt.figure(96)
    plt.title('Vertical')
    plt.plot(projection_vertical[i,:],'o',color=colors[i])

for i in range(int(len(frames_to_anal)/2)):
    plt.figure(97)
    plt.title('Horizontal')
    plt.plot(projection_horizontal[i,:],'o',color=colors[i])


#####Make Fits!!
#####Make Fits!!
sigma_vertical = []
sigma_horizontal = []
x_axis = timestamp_falling_edge*1
def gaussian(x, amp, cen, wid):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))
gmodel = Model(gaussian)


for i in range(int(len(frames_to_anal)/2)):
    tf = projection_vertical[i,:]
    x = np.linspace(0,len(tf)-1,len(tf))
    result = gmodel.fit(tf, x=x, amp=np.max(tf), cen=np.mean(x), wid=10)
    sigma_vertical = np.append(sigma_vertical,result.best_values['wid'])
    plt.figure(96)
    plt.title('Vertical')
    plt.plot(x, result.best_fit, '-', label='best fit',color=colors[i])

for i in range(int(len(frames_to_anal)/2)):
    tf = projection_horizontal[i,:]
    x = np.linspace(0,len(tf)-1,len(tf))
    result = gmodel.fit(tf, x=x, amp=np.max(tf), cen=np.mean(x), wid=10)
    sigma_horizontal = np.append(sigma_horizontal,result.best_values['wid'])
    plt.figure(97)
    plt.title('Horizontal')
    plt.plot(x, result.best_fit, '-', label='best fit',color=colors[i])

# ###normalize width
# sigma_vertical = sigma_vertical/np.max(sigma_vertical)
# sigma_horizontal = sigma_horizontal/np.max(sigma_horizontal)
#
def linear(x,m,n):
    """f(x) = mx + n"""
    return (m*x+n)
gmodel = Model(linear)


plt.figure(98)
result = gmodel.fit(sigma_vertical, x=t_dark_tot, m=-5, n=1)
n = result.best_values['n']
m_norm =  result.best_values['m']/n

plt.title('Vertical projection')
plt.plot(t_dark_tot,sigma_vertical/n,'o',color='blue')
plt.plot(t_dark_tot,linear(t_dark_tot,m_norm,1),'-', label='best fit',color='blue')
plt.xlabel('time (s)')
plt.ylabel('Norm. radial change')
print('Fit slope vertical:    '+str(m_norm))

plt.figure(99)
result = gmodel.fit(sigma_horizontal, x=t_dark_tot, m=-5, n=1)
n = result.best_values['n']
m_norm =  result.best_values['m']/n

plt.title('Horizontal projection')
plt.plot(t_dark_tot,sigma_horizontal/n,'o',color='blue')
plt.plot(t_dark_tot,linear(t_dark_tot,m_norm,1),'-', label='best fit',color='blue')
plt.xlabel('time (s)')
plt.ylabel('Norm. radial change')
print('Fit slope horizontal:  '+str(m_norm))





plt.show()
