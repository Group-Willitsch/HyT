<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 15:35:41 2019
@author: Gruppe Willitsch
"""

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
import sys

plt.close('all')


def shop_bin(data, bindeg):
#    bins = 1000/bindeg
    camera_res = 1000
    bins = camera_res/bindeg
    x_mesh,y_mesh = np.meshgrid(np.linspace(0,camera_res-1,camera_res),np.linspace(0,camera_res-1,camera_res))
    x_mesh = x_mesh.ravel()
    y_mesh = y_mesh.ravel()
    data = data[:camera_res,:camera_res]
    data = data*(data>0).astype(int)
    data = data*(data<100000).astype(int)
    data_hist,x,y = np.histogram2d(x_mesh,y_mesh,weights=data.ravel(),bins=(bins,bins))
    data_hist = data_hist.T
    return data_hist/(bindeg*bindeg)


def breakfast_at_tiffanys(filename):
    im = Image.open(filename)
    im.seek(0)

    imag   = np.array(im)

    bins = 4

    imag   = shop_bin(imag,bins)
    print(imag.shape)

    images = [imag]
    for i in range(1,int(1e6)):
        #if i%100 != 0:
        #    continue# for dem fasterz
        try:
            im.seek(i)

#           images = np.vstack(((images), np.array(im)))
            cut_im = shop_bin(np.array(im),bins)

            images.append( np.array(cut_im))
            #print(i)

        except:
            print('%i images in this tiff'%i)
            break
    return np.array(images),cut_im

def runningaverage(avgno, array):
    runnings = []
    for i in range(array.shape[0]-avgno):
        runnings.append(np.mean(array[i:i+avgno], axis = 0))
    return np.array(runnings)


#   filename = file1






directory = 'C:/Users/Gruppe Willitsch/Desktop/data/2020-01-17/'
file1 = directory + '000.tif'
images1,cut_im = breakfast_at_tiffanys(file1)
#np.save('C:/Users/Gruppe Willitsch/Desktop/data/2020-01-14/4.npy',images1)
#np.save('C:/Users/Gruppe Willitsch/Desktop/data/2020-01-14/4_cut_im.npy',cut_im)
#
#
#file_npy = directory + '4.npy'
#cuty_npy = directory + '4_cut_im.npy'
#images1 = np.load(file_npy)
#cut_im  = np.load(cuty_npy)
#
#
#import h5py
#hf = h5py.File('C:/Users/Gruppe Willitsch/Desktop/data/2020-01-14/4.h5', 'w')
#hf.create_dataset('data', data=images1,compression = 'gzip')
#hf.create_dataset('cut_image', data=cut_im,compression = 'gzip')
#hf.close()

#hf = h5py.File(file1, 'r')
#
#images1 = hf['data'][:]
#cut_im = hf['cut_image'][:]


sumim1 = np.mean(images1, axis = 0)
runningwater = runningaverage(20, images1)
plt.figure(5)
gs = plt.GridSpec(1,3)
fig = plt.figure(5)
ax1 = fig.add_subplot(gs[0,0])
plt.title('Up down')
plt.imshow(np.mean(runningwater,axis=1),aspect='auto',interpolation='none')
plt.colorbar()
ax1 = fig.add_subplot(gs[0,1])
plt.title('Left Right')
plt.imshow(np.mean(runningwater,axis=2),aspect='auto',interpolation='none')
plt.colorbar()

#normalization = np.mean(np.mean(runningwater[:,:15,-15:],axis=2),axis=1)
#normpic = runningwater[:,:15,-15:]
#m, nh, nw = normpic.shape
#normpic = normpic.reshape(m, nh * nw)
#normalization = np.mean(normpic, axis = 1)

normalization = np.mean(runningwater[:,-15:,:15],axis=(2,1))
runningwater_mean_norm = np.mean(runningwater,axis=2)
for i in range(len(normalization)):
    runningwater_mean_norm[i,:] = runningwater_mean_norm[i,:]/normalization[i]


ax1 = fig.add_subplot(gs[0,2])
plt.title('Left Right Norm')
plt.imshow(runningwater_mean_norm,aspect='auto',interpolation='none')
plt.colorbar()

plt.figure(6)
#plt.plot(np.mean(runningwater_mean_norm[:,100:150],axis=1))
#plt.plot(np.mean(runningwater_mean_norm[:,150:200],axis=1))
#plt.plot(np.mean(runningwater_mean_norm[:,0:40],axis=1))

plt.plot(np.mean(runningwater_mean_norm[:,25:125],axis=1))#/np.max(np.mean(runningwater_mean_norm[:,0:125],axis=1)))
plt.plot(np.mean(runningwater_mean_norm[:,125:225],axis=1))#/np.max(np.mean(runningwater_mean_norm[:,125:],axis=1)))



#fig = plt.figure(2)
#ax = plt.axes()
#im = ax.imshow(runningwater[0])#, vmax=800,vmin=0,cmap='cool')
#
#def init():
#    im.set_data(np.zeros(runningwater[0].shape))
#    return
#
#def update(frame):
#    print(frame)
#    imagedat = runningwater[frame,:,:]
##    ax.imshow(imagedat, vmax=800,vmin=300,cmap='hot')
#    im.set_data(imagedat)#, vmax=800,vmin=300,cmap='hot')
##    ax.colorbar()
#    return
#
#ani = animation.FuncAnimation(fig, update, init_func=init, frames = runningwater.shape[0], interval = 200, blit = False)
#fig.show()
#
#print('hola')

#plt.hist(images1[5].reshape([int(images1.shape[1]*images1.shape[2]),1]), bins = 200)
#plt.yscale('Log')
#plt.imshow(runningwater[3,:,:],norm=LogNorm())
#plt.colorbar()



'''
    positive minus negative images

'''





a  = (700,900)
b  = (500,700)

plt.axvspan(a[0],a[1], alpha=0.2, color='red')
plt.axvspan(b[0],b[1], alpha=0.2, color='blue')


#165-180 positives
positives = np.mean(images1[a[0]:a[1],:,:], axis = 0)
#245 - 260 negatives
negatives = np.mean(images1[b[0]:b[1],:,:], axis = 0)

y_cord,x_cord = np.meshgrid(np.linspace(0,249,250),np.linspace(0,249,250),indexing='ij')

posnegs = positives - negatives
plt.figure(7)
#plt.imshow(posnegs,cmap='bwr',vmin=-150,vmax=150)
fuCL = 150
plt.hist2d(x_cord.ravel(),y_cord.ravel(),weights=posnegs.ravel(),cmap='bwr',bins=50,vmax=fuCL*1,vmin=-(fuCL*1))
plt.colorbar()



=======
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 15:35:41 2019
@author: Gruppe Willitsch
"""

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
import sys

plt.close('all')


def shop_bin(data, bindeg):
#    bins = 1000/bindeg
    camera_res = 1000
    bins = camera_res/bindeg
    x_mesh,y_mesh = np.meshgrid(np.linspace(0,camera_res-1,camera_res),np.linspace(0,camera_res-1,camera_res))
    x_mesh = x_mesh.ravel()
    y_mesh = y_mesh.ravel()
    data = data[:camera_res,:camera_res]
    data = data*(data>0).astype(int)
    data = data*(data<100000).astype(int)
    data_hist,x,y = np.histogram2d(x_mesh,y_mesh,weights=data.ravel(),bins=(bins,bins))
    data_hist = data_hist.T
    return data_hist/(bindeg*bindeg)


def breakfast_at_tiffanys(filename):
    im = Image.open(filename)
    im.seek(0)

    imag   = np.array(im)

    bins = 4

    imag   = shop_bin(imag,bins)
    print(imag.shape)

    images = [imag]
    for i in range(1,int(1e6)):
        #if i%100 != 0:
        #    continue# for dem fasterz
        try:
            im.seek(i)

#           images = np.vstack(((images), np.array(im)))
            cut_im = shop_bin(np.array(im),bins)

            images.append( np.array(cut_im))
            #print(i)

        except:
            print('%i images in this tiff'%i)
            break
    return np.array(images),cut_im

def runningaverage(avgno, array):
    runnings = []
    for i in range(array.shape[0]-avgno):
        runnings.append(np.mean(array[i:i+avgno], axis = 0))
    return np.array(runnings)


#   filename = file1




directory = 'C:/Users/Gruppe Willitsch/Desktop/data/2020-01-13/'
file1 = directory + '10.tif'
images1,cut_im = breakfast_at_tiffanys(file1)


sumim1 = np.mean(images1, axis = 0)
runningwater = runningaverage(20, images1)
plt.figure(5)
gs = plt.GridSpec(1,3)
fig = plt.figure(5)
ax1 = fig.add_subplot(gs[0,0])
plt.title('Up down')
plt.imshow(np.mean(runningwater,axis=1),aspect='auto',interpolation='none')
plt.colorbar()
ax1 = fig.add_subplot(gs[0,1])
plt.title('Left Right')
plt.imshow(np.mean(runningwater,axis=2),aspect='auto',interpolation='none')
plt.colorbar()

#normalization = np.mean(np.mean(runningwater[:,:15,-15:],axis=2),axis=1)
#normpic = runningwater[:,:15,-15:]
#m, nh, nw = normpic.shape
#normpic = normpic.reshape(m, nh * nw)
#normalization = np.mean(normpic, axis = 1)

normalization = np.mean(runningwater[:,:15,-15:],axis=(2,1))
runningwater_mean_norm = np.mean(runningwater,axis=2)
#for i in range(len(normalization)):
#    runningwater_mean_norm[i,:] = runningwater_mean_norm[i,:]/normalization[i]


ax1 = fig.add_subplot(gs[0,2])
plt.title('Left Right Norm')
plt.imshow(runningwater_mean_norm,aspect='auto',interpolation='none')
plt.colorbar()

plt.figure(6)
#plt.plot(np.mean(runningwater_mean_norm[:,100:150],axis=1))
#plt.plot(np.mean(runningwater_mean_norm[:,150:200],axis=1))
#plt.plot(np.mean(runningwater_mean_norm[:,0:40],axis=1))

plt.plot(np.mean(runningwater_mean_norm[:,0:125],axis=1))#/np.max(np.mean(runningwater_mean_norm[:,0:125],axis=1)))
plt.plot(np.mean(runningwater_mean_norm[:,125:],axis=1))#/np.max(np.mean(runningwater_mean_norm[:,125:],axis=1)))

plt.show()


#fig = plt.figure(2)
#ax = plt.axes()
#im = ax.imshow(runningwater[0])#, vmax=800,vmin=0,cmap='cool')
#
#def init():
#    im.set_data(np.zeros(runningwater[0].shape))
#    return
#
#def update(frame):
#    print(frame)
#    imagedat = runningwater[frame,:,:]
##    ax.imshow(imagedat, vmax=800,vmin=300,cmap='hot')
#    im.set_data(imagedat)#, vmax=800,vmin=300,cmap='hot')
##    ax.colorbar()
#    return
#
#ani = animation.FuncAnimation(fig, update, init_func=init, frames = runningwater.shape[0], interval = 200, blit = False)
#fig.show()
#
#print('hola')

#plt.hist(images1[5].reshape([int(images1.shape[1]*images1.shape[2]),1]), bins = 200)
#plt.yscale('Log')
#plt.imshow(runningwater[3,:,:],norm=LogNorm())
#plt.colorbar()



'''
    positive minus negative images

'''


#165-180 positives
positives = np.mean(images1[700:800,:,:], axis = 0)
#245 - 260 negatives
negatives = np.mean(images1[500:600,:,:], axis = 0)

y_cord,x_cord = np.meshgrid(np.linspace(0,249,250),np.linspace(0,249,250),indexing='ij')

posnegs = positives - negatives
plt.figure(7)
#plt.imshow(posnegs,cmap='bwr',vmin=-150,vmax=150)
fuCL = 20
plt.hist2d(x_cord.ravel(),y_cord.ravel(),weights=posnegs.ravel(),cmap='bwr',bins=50,vmax=fuCL*150,vmin=-(fuCL*150))
plt.colorbar()

plt.show()
