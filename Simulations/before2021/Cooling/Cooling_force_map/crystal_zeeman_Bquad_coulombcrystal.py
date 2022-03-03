#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 15:38:09 2020

@author: claudiovonplanta
"""

import OBE
from numpy.linalg import norm
import numpy as np
from cazeeman_import import Bmodel
from cazeeman_import import Bvecs
from matplotlib import pyplot as plt
import scipy.constants as con
from numba import jit, njit
import time
#import h5py
import subprocess

# frequencies in MHz
# positions in mm
# seed for RNG
np.random.seed(0)

# some constants
muB = 9.27400994e-24    # J/T
amtoJ = 4.359744e-18    # conversion au to J
hbar = con.hbar
to_mm = 1e3
kb = con.k
wavelength_397 = 397e-9
k397 = 2*np.pi / wavelength_397
wavelength_866 = 866e-9
k866 = 2*np.pi / wavelength_866


laser397_pol = np.array([1/np.sqrt(2),0,1/np.sqrt(2)]) 
laser866_pol = np.array([1/np.sqrt(2),0,1/np.sqrt(2)]) 

plt.close('all')

'''

    model B field form cazeeman

'''

def sigv(T, m):
    return np.sqrt(kb * T/m)

def angle_3d(v1, v2):
    return np.arccos(np.dot(v1,v2) / (norm(v1) * norm(v2)))

def seady_state_OBE(u, alpha, beta, D1, D2, s1, s2, l1, l2):
    Ln = OBE.L_Matrix_oberst(u, alpha, beta, D1, D2, s1, s2, l1, l2)
    pops = OBE.steadystate_populations(Ln)
    return pops

def u_from_B(B):
    return B * muB/hbar / (2*np.pi)

def doppler(kz, vz):
    return -np.dot(kz,vz) / (2*np.pi)



positions = np.load('./positions_last.npy') * 1e-9 #convert to m
velocities = np.load('./velocities_last.npy') * 1e3 #convert to m/s
velocities = np.zeros_like(velocities)

xg = positions[:,0]
yg = positions[:,1]
zg = positions[:,2]

vx = velocities[:,0]
vy = velocities[:,1]
vz = velocities[:,2]

velocities = velocities.T

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.axes3d as p3

fig = plt.figure()
    
mins = np.min(positions* 1e6)
maxs = np.max(positions* 1e6)

ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(mins,maxs)
ax.set_ylim(mins,maxs)
ax.set_zlim(mins,maxs)

ax.set_xlabel('x /um')
ax.set_ylabel('y /um')
ax.set_zlabel('z /um')
ploet = ax.scatter(positions[:,0] * 1e6, positions[:,1]* 1e6, positions[:,2]* 1e6, marker = 'o', alpha = 0.7)
#ploet.set_sizes(1*std3d)
plt.tight_layout()
fig.savefig('xyz.png')
plt.show()

m = 40e-3 / con.N_A
# Bar magnets
#Bvec = Bvecs(xg,yg,zg).T
#Bmag = Bmodel(xg,yg,zg)

Bvec = Bvecs(xg,zg,yg).T
Bmag = Bmodel(xg,zg,yg)

us = u_from_B(Bmag) * 1e-6 # in  MHz    ###Zeeman shift
d1 = +40
d2 = 0
l1 = 0.3
l2 = 2
s1 = 1
s2 = 5
alphas = angle_3d(Bvec, laser397_pol)
betas = angle_3d(Bvec, laser866_pol)
Gamma = 2*np.pi*20.7e6      ###Life time p state
dt = 1e-9 #1e-9


dE_list = []

def steadystates(us, alphas, betas, d1s, d2s, s1, s2, l1, l2):
    pops = np.zeros((us.shape[0],8))

    for i in range(us.shape[0]):
        popslist = seady_state_OBE(us[i], alphas[i], betas[i], d1s[i], d2s[i], s1, s2, l1, l2)
        pops[i,:] = popslist
    return pops

i=0
#Dlist = np.linspace(-1000,1000,1)

kz = 2 * np.pi / 397e-9
kzvecn = np.array([1, 0, 1]) * 1/np.sqrt(2) 
kzvec = kz * kzvecn

kz2 = 2 * np.pi / 866e-9
kzvecn2 = np.array([1, 0, 1]) * 1/np.sqrt(2) 
kzvec2 = kz2 * kzvecn

Dlist = [0]
for d1 in Dlist:   
    # effective detunings with Doppler shifts
    d1s = d1 + doppler(kzvec, velocities)*1e-6 * 1
    d2s = d2 + doppler(kzvec2, velocities)*1e-6 * 1
    
    print('Nr %i in %i'%(i, len(Dlist)))
    i+=1
    
    time0 = time.time()
    pops = steadystates(us, alphas, betas, d1s, d2s, s1, s2, l1, l2)
    print('computation time:', time.time()-time0)
    if np.sum(pops<0)>1:
        print('numerical error')
    pstates = pops[:,2] + pops[:,3]
    

    negative_vz = vz<=0
    positive_vz = np.logical_not(negative_vz)
    
    pos_mean = np.mean(pstates[positive_vz])
    neg_mean = np.mean(pstates[negative_vz])
    
#    print('steady state population p-states with positive velocity', pos_mean)
#    print('steady state population p-states with negative velocity', neg_mean)
#    print('ratio neg/pos', neg_mean/pos_mean)
    
    vnew = velocities + pstates * Gamma * dt *  kzvec.reshape([3,1]) * hbar / m
    
    vold = velocities
    
    Ekin = 1/2 * m * np.linalg.norm(vnew, axis = 0)**2
    Ekin_old = 1/2* m * np.linalg.norm(vold, axis = 0)**2
    
    Ediff = Ekin - Ekin_old
    
    vmag = np.linalg.norm(vold, axis = 0)
    vmag_new = np.linalg.norm(vnew, axis = 0)
    
#    print('old average velocity magnitude', np.mean(vmag))
#    print('new average velocity magnitude', np.mean(vmag_new))
#    print('difference', np.mean(vmag_new)-np.mean(vmag))
    
#    print('mean vz * pstates', np.mean(test))
    Enew = 1/2 * m * vmag_new**2
    Eold = 1/2 * m * vmag**2
    
    Enew_mean = np.mean(Enew)
    Eold_mean = np.mean(Eold)
    dE = Enew_mean-Eold_mean
    dE_percent = dE/Eold_mean * 100
    
    print('Average kinetic energy before:', Eold_mean)
    print('Average kinetic energy after:', Enew_mean)
    print('new - old', dE)
    print('new - old percent', np.abs(dE)/Eold_mean * 100)
    print('\n')
    
    dE_list.append(dE_percent)


import matplotlib.cm as cmx
import matplotlib

crit = positions[:,1] > 0

cs = Ediff / np.std(Ediff)
cs = cs[crit]

colorsMap='jet'
cm = plt.get_cmap(colorsMap)
cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(mins,maxs)
ax.set_ylim(mins,maxs)
ax.set_zlim(mins,maxs)



ax.set_xlabel('x /um')
ax.set_ylabel('y /um')
ax.set_zlabel('z /um')
ploet = ax.scatter(positions[:,0][crit] * 1e6, positions[:,1][crit]* 1e6, positions[:,2][crit]* 1e6, marker = 'o', alpha = 0.7, c=scalarMap.to_rgba(cs))

scalarMap.set_array(cs)
fig.colorbar(scalarMap)
plt.show()



out = subprocess.Popen(['rm', '-R', './out2'], 
               stdout=subprocess.PIPE, 
               stderr=subprocess.STDOUT)


out = subprocess.Popen(['cp', '-R', './out', './out2'], 
               stdout=subprocess.PIPE, 
               stderr=subprocess.STDOUT)

out = subprocess.Popen(['rm', '-R', './out'], 
               stdout=subprocess.PIPE, 
               stderr=subprocess.STDOUT)

out = subprocess.Popen(['mkdir', './out'], 
               stdout=subprocess.PIPE, 
               stderr=subprocess.STDOUT)

fig, ax = plt.subplots()
plt.plot(Dlist, dE_list, 'ro')
plt.xlabel('detuning 1/MHz')
plt.ylabel('efficiency')
#plt.xscale('log')
plt.tight_layout()
fig.savefig('./out/eff_d1.png')
plt.show()


np.save('./out/Dlist',Dlist)
np.save('./out/eff', dE_list)
np.save('./out/ulist', us)
np.save('./out/Bmag', Bmag)
np.save('./out/Bvec', Bvec)
np.save('./out/alphas', alphas)
np.save('./out/betas', betas)


out = subprocess.Popen(['cp', 'crystal_zeeman_Bhom.py', './out/'], 
               stdout=subprocess.PIPE, 
               stderr=subprocess.STDOUT)
stdout,stderr = out.communicate()