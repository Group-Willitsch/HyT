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
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.ndimage import gaussian_filter1d


def bug():
    plt.show()
    sys.exit()

def gaussian(x, amp, cen, wid):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))
gmodel = Model(gaussian)
plt.close('all')
colors = ['black','red','green','blue','orange','yellow','gray']
colors = np.append(colors,colors)
colors = np.append(colors,colors)


####This is a code to make a movie of the ion cloud and do some magig.

###File to analyze
loc     =["/Volumes/sw_group/pc-swplab-pc08/"]
date     = "2020-07-01"
file     = ['002']

###Try to open the file you would like to analyze
for i in range(len(loc)):
    try:
        link = loc[i]+date+'/'+date+'_'+file[0]+'.h5'
        h5_data = h5py.File(link,'r')
        print('Analyzing file: '+link)
    except:
        print('Failed in getting the file')
        pass

fgc = 0                                         ####Local variable for the figures
nbr_frames = len(list(h5_data.keys()))          ####Get the number of frames
frames = []
for i in range(0,nbr_frames):                   ####This loop is used to correct order of frames. list(h5_data.keys())
    frames = np.append(frames,'data'+str((i)))  ####produces an order of, for example, frame0,frame10,frame100,frame101 ect


########Background subtraction?----
########The background subtraction is done by making an average of the last frames--assuming that the signal is
########Gone in this
background = False
bg = np.zeros(np.shape(h5_data[str(frames[0])][:]))
si = np.zeros(np.shape(h5_data[str(frames[0])][:]))
bg_frames = 5          ###----Nbr of BG frames
if background:
    bg_frame_list = frames[-bg_frames:]                ####Have a look at the last x-frames for background
    for i in range(len(bg_frame_list)):
        bg = bg + h5_data[str(bg_frame_list[i])][:]
    bg = bg/bg_frames

    bg_frame_list = frames[:bg_frames]                 ####Have a look at the first x-frames for signal
    for i in range(len(bg_frame_list)):
        si = si + h5_data[str(bg_frame_list[i])][:]
    si = si/bg_frames

    fgc+=1
    plt.figure(fgc)
    plt.title('Background')
    plt.imshow(bg.T,origin='lower')
    plt.colorbar()

    fgc+=1
    plt.figure(fgc)
    plt.title('Raw Signal')
    plt.imshow(si.T,origin='lower')
    plt.colorbar()

    print('Please note, you subtract this backgroud')




########Loop over the frames----
########Loop over the frames----
sigma_horizontal = []   ###Variable to store the value of the fit for the width
sigma_vertical   = []   ###Variable to store the value of the fit for the width
mean = []
time = []
t_start_measurment = h5_data['data0'].attrs['timestamp'][:19]
frames_for_movie = []

print(np.datetime64(t_start_measurment))

for i in range(0,len(frames),4):
    ####---Get the image of the frame
    data = h5_data[str(frames[i])][:]-bg
    ####---Get the time of the frame with respect to t_start_measurment
    frames_for_movie.append(data[100:400,170:320])

    mean.append(np.mean(data[100:400,170:320]))

    if i > 0:
        time = np.append(time,(np.datetime64(h5_data[str(frames[i])].attrs['timestamp'][:19])-np.datetime64(t_start_measurment))/np.timedelta64(1, 's'))
    else:
        time = np.append(time,0)

    projection_horizontal = np.mean(data,axis=1)
    projection_horizontal = gaussian_filter1d(projection_horizontal,1)-np.mean(projection_horizontal)
    projection_vertical   = np.mean(data,axis=0)
    projection_vertical   = gaussian_filter1d(projection_vertical,1)-np.mean(projection_vertical)
    #
    #
    #
    # tf = projection_horizontal
    # x = np.linspace(0,len(tf)-1,len(tf))
    # result = gmodel.fit(tf, x=x, amp=np.max(tf), cen=np.argmax(tf), wid=2)
    # sigma_horizontal = np.append(sigma_horizontal,result.best_values['wid'])
    #
    # tf = projection_vertical
    # x = np.linspace(0,len(tf)-1,len(tf))
    # result = gmodel.fit(tf, x=x, amp=np.max(tf), cen=np.argmax(tf), wid=2)
    # sigma_vertical = np.append(sigma_vertical,result.best_values['wid'])


fgc+=1
plt.figure(fgc)
plt.plot(time,mean)


fgc+=1
gs = plt.GridSpec(1,1)
fig = plt.figure(fgc)#,figsize=(18,5))

pos_y, pos_x = np.meshgrid(np.linspace(0,len(frames_for_movie[0][0])-1,len(frames_for_movie[0][0])),np.linspace(0,len(frames_for_movie[0])-1,len(frames_for_movie[0])))

def update_graph(num):
	movie_frame_raw = frames_for_movie[num][:][:]
	movie_frame,x,y = np.histogram2d(pos_x.ravel(),pos_y.ravel(),weights=movie_frame_raw.ravel(),bins=25)
	movie_frame = movie_frame/72
    #print(shape(movie_frame))

	plt.title('timestamp: %.3f'%(time[num]))
	im1 = plt.imshow(movie_frame.T,origin='lower',vmax=900,vmin=550,cmap='Blues')
    #plt.imshow(data[100:400,170:320],origin='lower',vmax=800,vmin=600,cmap='Blues')


#anim = FuncAnimation(fig, update_graph, frames=np.arange(0, moviedat.shape[0]),interval=5, blit=False)
anim = FuncAnimation(fig, update_graph,frames=np.arange(0,len(time)),interval=1,blit=False)



anim.save('2020-07-01.gif', writer='imagemagick', fps=20)
#plt.show()
# fgc+=1
# plt.figure(fgc)
# if background:
#     plt.title('Signal-Back')
# else:
#     plt.title('Raw data')
# plt.imshow(data.T,origin='lower')
# plt.colorbar()
#
# fgc+=1
# plt.figure(fgc)
# plt.plot(time[4:],test[4:])
# plt.grid()
# plt.xlabel('Time (s)')
# plt.ylabel('Intensity (arb. u.)')
bug()
