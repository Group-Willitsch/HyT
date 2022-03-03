# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 14:51:27 2020

@author: Gruppe Willitsch
"""
import numpy as np
from matplotlib import pyplot as plt
import h5py

plt.close('all')
day = '2020-08-06'
location = "C:/Users/Gruppe Willitsch/Desktop/data/"+day+"/"
file     = [day+'_007.h5',day+'_008.h5',day+'_001.h5',day+'_001.h5']
#### ca, oven, bg, ls   ---- Historic labelling of the order 


laser = 0
bg    = 0
oven  = 0
ca = 0

all_frames = True
nbr_frms_to_look = 50


for i in range(len(file)):
    filename = location+file[i]
    data = h5py.File(filename, "r")
    frames = len(list(data.keys()))
    
    if i == 0:
        frames_ca = frames
        #if all_frames:
        for j in range(int(frames)):
        #for j in range(0,nbr_frms_to_look):
            if j == 0:
                ca = data['data%s'%j][:]
            else:
                ca = ca + data['data%s'%j][:]
    if i == 1:
        frames_ov = frames
        for j in range(int(frames)):
        #for j in range(0,nbr_frms_to_look):
            if j == 0:
                oven = data['data%s'%j][:]
            else:
                oven = oven + data['data%s'%j][:]
    
    if i == 2:
        frames_bg = frames
        for j in range(int(frames)):
        #for j in range(0,nbr_frms_to_look):
            if j == 0:
                bg = data['data%s'%j][:]
            else:
                bg = bg + data['data%s'%j][:]
                
    if i == 3:
        frames_la = frames
        for j in range(int(frames)):
        #for j in range(0,nbr_frms_to_look):
            if j == 0:
                laser = data['data%s'%j][:]
            else:
                laser = laser + data['data%s'%j][:]

laser = laser/frames_ca
bg = bg/frames_ov
oven  = oven/frames_bg
ca = ca/frames_la

ca_diff = ca # + oven/2 #- laser + bg
#ca_diff = ca-bg
print(np.mean(ca_diff))

#### ca, oven, bg, ls   ---- Historic labelling of the order 
plt.figure(1)
plt.imshow((ca).T,origin='lower',vmax=700)#,vmin=1000)
plt.colorbar()

plt.figure(2)
plt.imshow((oven).T,origin='lower',vmax=700)#,vmin=1000)
plt.colorbar()

plt.figure(3)
plt.imshow((bg).T,origin='lower',vmax=900)#,vmin=1000)
plt.colorbar()

plt.figure(5)
plt.plot(np.mean(ca,axis=0)/np.max(np.sum(ca)))
plt.plot(np.mean(oven,axis=0)/np.sum(oven))


plt.figure(6)
plt.plot(np.mean(ca,axis=1)/np.max(np.sum(ca)))
plt.plot(np.mean(oven,axis=1)/np.sum(oven))


#plt.figure(1)
#plt.imshow((np.flipud(ca)).T,origin='lower',vmax=700,vmin=500)
#plt.colorbar()
#
#plt.figure(2)
#plt.imshow((oven).T,origin='lower',vmax=700,vmin=500)
#plt.colorbar()
#plt.figure(3)
#plt.imshow((np.fliplr(bg)).T,origin='lower',vmax=700,vmin=500)
#plt.colorbar()
#
#plt.figure(4)
#plt.imshow((laser).T,origin='lower',vmax=700,vmin=500)
#plt.colorbar()

#
#
#plt.figure(2)
#plt.hist(ca_diff.ravel(),bins=200)
#plt.yscale('log')
#
#
x,y = np.meshgrid(np.linspace(0,999,len(ca)),np.linspace(0,999,len(ca)))
scale = 10
plt.figure(4)
plt.hist2d(y.ravel(),x.ravel(),weights=(ca-oven).ravel(),bins=len(ca)/scale,vmax=200,vmin=-200)#,cmap='RdBu_r',vmax=20,vmin=-20)
plt.colorbar()
#
#
#plt.figure(4)

#plt.plot(np.mean(ca[200:340,200:300],axis=1),label='001')#,cmap='RdBu_r',vmax=800,vmin=600)
#plt.plot(np.mean(oven[200:340,200:300],axis=1),label='002')#,cmap='RdBu_r',vmax=800,vmin=600)
#plt.plot(np.mean(bg[200:340,200:300],axis=1),label='003')#,cmap='RdBu_r',vmax=800,vmin=600)
#plt.plot(np.mean(laser[200:340,200:300],axis=1),label='004')#,cmap='RdBu_r',vmax=800,vmin=600)
#plt.legend()
#
#
#print('Signal:  '+str(np.round(np.mean(ca_diff[0:300,400:600]),2)))
#print('Back:  '+str(np.round(np.mean(oven[0:300,400:600]),2)))
#print("SNR:  "+(str(np.round(np.mean(ca_diff[0:300,400:600])/np.mean(oven[0:300,400:600]),6))))
plt.show()
  