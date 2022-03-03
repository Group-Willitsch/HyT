#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:31:42 2019

@author: claudiovonplanta
"""
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import numpy as np
from numba import jit
import time

numpoints = 250

muB = 9.274009994e-24
hbar = 1.054571800e-34
massCa = 40e-3 / 6.022e23

L2 = 1
S2 = 1/2
J2 = 1/2
mJ2 = 1/2

L1 = 0
S1 = 1/2
mJ1 = 1/2
J1 = 1/2

rad = 0.018


wr = 396.950150e-9
cw = 299792458
fr = cw/wr
kr = 2 * np.pi / wr
kb = 1.38064852e-23
mCa = 40e-3 / 6.022e23
Temp = 50e-3
detuning = -20e6
fa = fr + detuning
wa = cw / fa
ka = 2 * np.pi / wa


'''

    Read all the parameters from the mathematica fit to the magnetic field function in a cube of
    +/- 0.5 mm x,y,z with steps of 0.05 mm
    x,y,z polynominal 2nd degree
    
    Note that y and z axes have to be exchanged to get the same coordinate system (as seen in plot)

'''

end = 2

Bxpars = [-4.166797649282299e-17,4.499791541407458e-17,-2.0095955963780472e-16,-1.0339992589091577e-16,2.265386737793366e-16,
     5.430314042260869e-16,-7.515927877902493e-17,-1.3032753701426085e-15,1.9891947039050155e-15,-0.172394683436446,
     1.2270844829714065e-16,-0.05933537890945301,-1.258548187662981e-16,-1.6627391415388717e-16,-8.1991940486209e-16,
     0.007019312353927059,-1.1273891816853738e-15,0.0006236809383367347,-2.9593966019241066e-16,-3.4909161700248443e-16,
     2.0657021925167467e-15,1.3187905531204967e-16,-3.459034989261942e-15,-2.274263263503015e-15,4.590449316703882e-16,
     6.696441831425545e-15,-1.931421090480607e-14]

index = 0
for i in range(end+1):
        for j in range(end+1):
            for k in range(end+1):
                exec('ax%i%i%i = Bxpars[index]'%(i,j,k))
                index += 1
                
Bypars = [-2.6431822114139612e-17,0.2505758825134559,3.417333416015495e-16,-7.814337511744594e-17,0.,0.,2.466163834936756e-16,
     -0.0010885493672717998,-1.071104840564239e-15,9.145065278463794e-17,-2.0136771002607696e-16,-5.585465872039751e-16,
     8.809837313640867e-17,3.3254782830777434e-16,-3.1003202496347777e-15,-1.2722450041868321e-15,1.6910837725280605e-15,
     4.2958306088390296e-15,1.691083772528061e-16,-0.058545155723230824,-3.978389407810031e-15,-1.3963664680099378e-16,
     7.686744420582094e-16,3.285046936171022e-15,0.,0.0006069539174226405,0.]

index = 0
for i in range(end+1):
        for j in range(end+1):
            for k in range(end+1):
                exec('ay%i%i%i = Bypars[index]'%(i,j,k))
                index += 1
  
Bzpars = [3.0308293500748444e-17,5.667041606062843e-17,-1.4953570172530723e-16,-0.07236407296250044,-2.01367710026077e-16,
     -0.0013274599507121045,-3.4761166435299033e-16,-6.555164808157763e-16,1.8313980086433195e-15,-3.810443866026581e-18,
     -5.191511274109796e-17,-1.3575785105652173e-17,-3.2092978785406016e-16,0.,2.0241760307532845e-15,2.947884765798758e-16,
     2.562248140194031e-16,-5.053918363340035e-16,-7.985673370271399e-17,-4.111523489140372e-16,0.,0.007085118636322758,
     5.1244962803880624e-17,0.0006534903811673893,1.3771347950111646e-15,4.1694826497555284e-15,-6.230390614453571e-15]
              
index = 0
for i in range(end+1):
        for j in range(end+1):
            for k in range(end+1):
                exec('az%i%i%i = Bzpars[index]'%(i,j,k))
                index += 1
#@jit             
def Bxmodel(x,y,z):
    ret = (ax000 + z*ax001 + z**2*ax002 + y*ax010 + y*z*ax011 + 
        y*z**2*ax012 + y**2*ax020 + y**2*z*ax021 + y**2*z**2*ax022 + 
        x*ax100 + x*z*ax101 + x*z**2*ax102 + x*y*ax110 + 
        x*y*z*ax111 + x*y*z**2*ax112 + x*y**2*ax120 + x*y**2*z*ax121 + 
        x*y**2*z**2*ax122 + x**2*ax200 + x**2*z*ax201 + x**2*z**2*ax202 + 
        x**2*y*ax210 + x**2*y*z*ax211 + x**2*y*z**2*ax212 + x**2*y**2*ax220 + 
        x**2*y**2*z*ax221 + x**2*y**2*z**2*ax222)     
    return ret

#@jit             
def Bymodel(x,y,z):
    ret = (ay000 + z*ay001 + z**2*ay002 + y*ay010 + y*z*ay011 + 
        y*z**2*ay012 + y**2*ay020 + y**2*z*ay021 + y**2*z**2*ay022 + 
        x*ay100 + x*z*ay101 + x*z**2*ay102 + x*y*ay110 + 
        x*y*z*ay111 + x*y*z**2*ay112 + x*y**2*ay120 + x*y**2*z*ay121 + 
        x*y**2*z**2*ay122 + x**2*ay200 + x**2*z*ay201 + x**2*z**2*ay202 + 
        x**2*y*ay210 + x**2*y*z*ay211 + x**2*y*z**2*ay212 + x**2*y**2*ay220 + 
        x**2*y**2*z*ay221 + x**2*y**2*z**2*ay222)     
    return ret

#@jit             
def Bzmodel(x,y,z):
    ret = (az000 + z*az001 + z**2*az002 + y*az010 + y*z*az011 + 
        y*z**2*az012 + y**2*az020 + y**2*z*az021 + y**2*z**2*az022 + 
        x*az100 + x*z*az101 + x*z**2*az102 + x*y*az110 + 
        x*y*z*az111 + x*y*z**2*az112 + x*y**2*az120 + x*y**2*z*az121 + 
        x*y**2*z**2*az122 + x**2*az200 + x**2*z*az201 + x**2*z**2*az202 + 
        x**2*y*az210 + x**2*y*z*az211 + x**2*y*z**2*az212 + x**2*y**2*az220 + 
        x**2*y**2*z*az221 + x**2*y**2*z**2*az222)     
    return ret


#@jit
def Bvecs(x,y,z):
    return np.array([Bxmodel(x,y,z), Bymodel(x,y,z), Bzmodel(x,y,z)])
#@jit
def Bmodel(x,y,z):
    return np.sqrt(Bxmodel(x,y,z)**2 + Bymodel(x,y,z)**2 + Bzmodel(x,y,z)**2)

#@jit
def gJ(L,S,J):
    return 1 + (J * (J+1) + S * (S+1) - L * (L+1))/(2 * J * (J+1))

#@jit
def Zeemanfactor(L,S,J,mJ):
    return muB * gJ(L,S,J) * mJ

#@jit
def Zeemanfactor2lv(L1, S1, J1, mJ1, L2, S2, J2, mJ2):
    return Zeemanfactor(L2, S2, J2, mJ2) - Zeemanfactor(L1, S1, J1, mJ1)





###note b




#print(gJ(L1,S1,J1))
#print(gJ(L2,S2,J2))
#
#
#print(Zeemanfactor(L1,S1,J1,mJ1)/hbar * 1e-4)
#print(Zeemanfactor(L2,S2,J2,mJ2)/hbar * 1e-4)
#
#print(Zeemanfactor2lv(L1, S1, J1, mJ1, L2, S2, J2, mJ2)/hbar * 1e-4)
#print(Zeemanfactor2lv(L1, S1, J1, -mJ1, L2, S2, J2, mJ2)/hbar * 1e-4)
#print(Zeemanfactor2lv(L1, S1, J1, mJ1, L2, S2, J2, -mJ2)/hbar * 1e-4)
#print(Zeemanfactor2lv(L1, S1, J1, -mJ1, L2, S2, J2, -mJ2)/hbar * 1e-4)



#@jit
def sigv(T, m):
    return np.sqrt(kb * T/m)
#@jit
def doppler(k, v):
    return -k * v

def find_surface(numpoints, rad):
    xpoints = np.linspace(0,rad,numpoints)    
    ypoints = np.linspace(0,rad,numpoints)    
    zpoints = np.linspace(0,rad,numpoints)      
    xp,yp,zp = np.meshgrid(xpoints, ypoints, zpoints, indexing = 'ij')
    
    Zeemanenergies = Bmodel(xp, yp, zp) *  Zeemanfactor2lv(L1, S1, J1, -mJ1, L2, S2, J2, mJ2)/hbar 
    
    '''
    
        Dopplershifts
    
    '''    
    vel = 2 * sigv(Temp,mCa)
#    dr = doppler(kr, vel)
    da = doppler(ka, vel)
    
    gate = Zeemanenergies < np.abs(da)
    
    xgate, ygate, zgate = xp[gate], yp[gate], zp[gate]
    
    #rspher = np.sqrt(xgate**2 + ygate**2 + zgate**2)
    #phi = np.arctan(ygate / xgate)
    #theta = np.arccos(zgate / rspher)
    
    #print('xmin', np.min(np.abs(xgate)))
    print('xmax /um', np.max(np.abs(xgate))*1e3)
    
    #print('ymin', np.min(np.abs(ygate)))
    print('ymax /um', np.max(np.abs(ygate))*1e3)
    
    #print('zmin', np.min(np.abs(zgate)))
    print('zmax /um', np.max(np.abs(zgate))*1e3)
    
    
    #surface
    xy = np.unique(np.vstack((xgate,ygate)).T, axis = 0)
    volume = np.vstack((xgate,ygate,zgate))
    
    z = np.zeros(xy.shape[0])
    for i in range(len(xy)):
        z[i] = np.max(volume[2][np.logical_and(volume[0]==xy[i,0], volume[1] == xy[i,1])])
        
    return volume, xy, z, np.max(np.abs(xgate)), np.max(np.abs(ygate)), np.max(np.abs(zgate))

'''

    Plot volume where Zeemanshift < Dopplershift

'''
if __name__ == "__main__":

    #plt.plot(xmaxes)
    #plt.plot(ymaxes)
    #plt.plot(zmaxes)
    #plt.show()
    #
    #
    #Plot all points#
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #
    #mins = np.min([np.min(xgate), np.min(ygate), np.min(zgate)])
    #maxs = np.max([np.max(xgate), np.max(ygate), np.max(zgate)])
    #
    #ax.set_xlim(mins * 1e3 , maxs * 1e3)
    #ax.set_ylim(mins * 1e3 , maxs * 1e3)
    #ax.set_zlim(mins * 1e3 , maxs * 1e3)
    #
    #ax.set_xlabel('x /um')
    #ax.set_ylabel('z /um')
    #ax.set_zlabel('y /um')
    #
    #ploet = ax.scatter(xgate * 1e3, ygate * 1e3, zgate * 1e3, marker = 'o', alpha = 0.7)
    ##ploet.set_sizes(1e-4)
    #plt.tight_layout()
    #plt.show()
    
    #find surface and plot
    time0 = time.time()
    volume, xy, z, xmaxes, ymaxes, zmaxes = find_surface(numpoints, rad)
    print('time passed:', time.time()-time0)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    maxs = np.max([xmaxes, ymaxes, zmaxes])
    
    ax.set_xlim(0 , maxs * 1e3)
    ax.set_ylim(0 , maxs * 1e3)
    ax.set_zlim(0 , maxs * 1e3)
    
    ax.set_xlabel('x /um')
    ax.set_ylabel('z /um')
    ax.set_zlabel('y /um')
    
    ploet = ax.scatter(xy.T[0]*1e3, xy.T[1]*1e3, z*1e3, marker = 'o', alpha = 0.7)
    #ploet.set_sizes(1e-4)
    plt.tight_layout()
    plt.show()