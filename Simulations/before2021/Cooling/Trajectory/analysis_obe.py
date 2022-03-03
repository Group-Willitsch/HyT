#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 13:29:13 2019

@author: vonplanta
"""
import numpy as np
import math
import matplotlib as mpl
import scipy.constants as con

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.axes3d as p3
import sys

def bug():
    plt.show()
    sys.exit()
    
 
print('analysis...')
plt.close('all')
mpl.rcParams['agg.path.chunksize'] = 10000

####---- "Global variables"
save = False
fgc  = 0         ##fgc = (manual) figure counter




#~ std3d = np.load('./xyz_std.npy')
#~ plot3ddata = np.load('./xyz.npy')
rffreq = 8e6     
wavelength_resonance = 396.950150e-9
freq_resonance = con.c / wavelength_resonance
detuning = - 20e6
freq_abs = freq_resonance + detuning
wavelength_abs = con.c / freq_abs
kz = 2 * np.pi / wavelength_abs
linewidth =  132e6#20.7e6 * 2*self.localPi
v = np.sqrt(3 * con.k)
plots = True
movie = False
to_um = 1e6
#~ Fit an ellipsoid to the crystal
#~ plot3ddata = plot3ddata * to_um


if plots:    
    
    gs = plt.GridSpec(2,3)
    fig = plt.figure(fgc,figsize=(13,7))
    #plt.title('Positions')
    
    posxlist = np.load('./xpos.npy')
    ax1 = fig.add_subplot(gs[0,0])
    plt.ylabel('xpos')
    plt.title('xpos')
    plt.plot(posxlist*to_um)
    plt.grid()

    posylist = np.load('./ypos.npy')
    ax1 = fig.add_subplot(gs[0,1])
    plt.ylabel('ypos')
    plt.plot(posylist*to_um)
    plt.title('y-pos')
    plt.grid()
    
    poszlist = np.load('./zpos.npy')
    ax1 = fig.add_subplot(gs[0,2])
    plt.ylabel('zpos')
    plt.plot(poszlist*to_um)
    plt.title('z-pos')
    plt.grid()
    
    velxlist = np.load('./vx.npy')
    ax1 = fig.add_subplot(gs[1,0])
    plt.ylabel('xvel (m/s)')
    plt.plot(velxlist)
    plt.title('x-vel')
    plt.grid()
    
    velylist = np.load('./vy.npy')
    ax1 = fig.add_subplot(gs[1,1])
    plt.ylabel('yvel (m/s)')
    plt.plot(velylist)
    plt.title('y-vel')
    plt.grid()
 
    
    velzlist = np.load('./vz.npy')
    ax1 = fig.add_subplot(gs[1,2])
    plt.ylabel('zvel (m/s)')
    plt.plot(velzlist)
    plt.title('z-vel')
    plt.grid()

    if save:
        plt.tight_layout()
        fig.savefig('pos_vel.png')
        

    fgc+=1
    velxz = np.sqrt(velxlist**2 + velzlist**2)
    k = 1/np.sqrt(2) * np.array([1,0,1])
    velxz = np.dot(k, np.array([velxlist,np.zeros_like(velxlist),velzlist]))
    fig = plt.figure(fgc)
    plt.ylabel('xzvel (m/s)')
    plt.plot(velxz)
    plt.grid()
    if save:
        plt.tight_layout()
        fig.savefig('xzvel.png')


    fgc+=1
    pstateslist = np.load('./pstates.npy')
    fig = plt.figure(fgc)
    plt.ylabel('pstates')
    plt.plot(pstateslist)
    plt.grid()
    if save:
        plt.tight_layout()
        fig.savefig('./pstates.png', dpi = 300, bbox_inches='tight')
        



    fgc+=1
    numphot = np.load('./numphot.npy')
    fig = plt.figure(fgc)
    plt.ylabel('numphot')
    plt.plot(numphot)
    if save:
        plt.tight_layout()
        fig.savefig('./numphot.png', dpi = 300, bbox_inches='tight')


    fgc+=1
    fig = plt.figure()
    d1s = np.load('./d1s.npy')
    plt.ylabel('detunings 1')
    plt.plot(d1s)
    if save:
        plt.tight_layout()
        fig.savefig('./d1s.png', dpi = 300, bbox_inches='tight')

    fgc+=1
    fig = plt.figure()    
    d2s = np.load('./d2s.npy')
    plt.ylabel('detunings 2')
    plt.plot(d2s)
    if save:
        plt.tight_layout()
        fig.savefig('./d2s.png', dpi = 300, bbox_inches='tight')


    fgc+=1
    fig = plt.figure()    
    ekin = np.load('./ekin.npy')
    plt.ylabel('ekin')
    plt.plot(ekin)
    if save:
        plt.tight_layout()
        fig.savefig('./ekin.png', dpi = 300, bbox_inches='tight')



    fgc+=1
    fig = plt.figure()    
    etot = np.load('./etot.npy')
    plt.ylabel('etot')
    plt.plot(etot)
    if save:
        plt.tight_layout()
        fig.savefig('./etot.png', dpi = 300, bbox_inches='tight')

    
    dt = 5e-9
    rffreq = 8e6
    period = 1/rffreq
    cyclesteps = int(np.floor(period / dt))
    length = len(etot)
    numcycles = int(np.floor(length / cyclesteps))
    cycleenergy = [np.mean(etot[i*cyclesteps:(i+1)*cyclesteps]) for i in range(numcycles)]
    cyclekin = [np.mean(ekin[i*cyclesteps:(i+1)*cyclesteps]) for i in range(numcycles)]
    
    fgc+=1
    fig = plt.figure()    
    plt.ylabel('etot cycle')
    plt.plot(cycleenergy)
    if save: 
        plt.tight_layout()
        fig.savefig('./etot_cycle.png', dpi = 300, bbox_inches='tight')

    fgc+=1
    fig = plt.figure()    
    plt.ylabel('ekin cycle')
    plt.plot(cyclekin)
    if save:
        plt.tight_layout()
        fig.savefig('./ekin_cycle.png', dpi = 300, bbox_inches='tight')

    
if movie:
	from matplotlib.animation import FuncAnimation

	moviedat = np.load('./positions.npy')
	moviedat = moviedat * 1e-3
	print('shape of moviedat', moviedat.shape)
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	
	# Setting the axes properties
	#~ ax.set_xlim3d([0.0, 1.0])
	#~ ax.set_xlabel('X')

	#~ ax.set_ylim3d([0.0, 1.0])
	#~ ax.set_ylabel('Y')

	#~ ax.set_zlim3d([0.0, 1.0])
	#~ ax.set_zlabel('Z')

	#~ ax.set_title('3D Test')

	ploet = ax.scatter(moviedat[0,:,0], moviedat[0,:,1], moviedat[0,:,2], marker = 'o', alpha = 0.06)
	plt.title('timestamp: %.3f'%(50 * 0) + 'us')

	fig.set_tight_layout(True)
	#~ lines = [ax.plot(dat[:,0], dat[:,1], dat[:,2])[0] for dat in moviedat]

	#~ plt.show()
	def update_graph(num):
		plt.title('timestamp: %.3f'%(50 * num) + 'us')
		ploet._offsets3d = (moviedat[num,:,0], moviedat[num,:,1], moviedat[num,:,2])

    
	anim = FuncAnimation(fig, update_graph, frames=np.arange(0, moviedat.shape[0]),
                                   interval=200, blit=False)
	plt.show()
	
