# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 13:56:48 2019

@author: Gruppe Willitsch
"""
#    adapted from Gregor Hegi: Ion_Trap_Stability.vi


import numpy as np
from matplotlib import pyplot as plt
# MHz, mm, V, u
Vrf = 86.5*110/65
Uend = 2.
UDC = 0
r0 = np.sqrt(2)
z0 = (25-(2*3))/2
kappa = 0.485165
Omega = 2*np.pi*11.
const = 9.77602310141 * 4 * np.pi**2# 4 e/(amu 4 pi^2 MHz^2 mm^{-2})
mass = 40
charge = 1
phase = 1.0 # single: 0.5, two: 1.0
NA = 6.022e23

# for old equation
e_c =  1.602176634e-19
mass_o = mass*1e-3/NA;
Omega_o = Omega*1e6
r0_o = r0*1e-3
z0_o = z0*1e-3

def mathieuparams_greg(mass,UDC,Vrf,Omega,Uend):
#    adapted from Gregor Hegi: Ion_Trap_Stability.vi

#    qx = phase * const * charge * Vrf / ( mass * Omega**2 * r0**2 )
#    ax = -const * charge / ( mass * Omega**2 ) * ( 2 * UDC / r0**2 + kappa * Uend / z0**2)
#    az = const * 2 * charge * kappa * Uend / ( mass * Omega**2 * z0**2)
    
    qx = phase * const * charge * Vrf / ( mass * Omega**2 * r0**2 )
    ax = -const * charge / ( mass * Omega**2 ) * ( 2 * UDC / r0**2 + kappa * Uend / z0**2)
    az = const * 2 * charge * kappa * Uend / ( mass * Omega**2 * z0**2)
    
    vr = Omega / 2 * 1000 * np.sqrt(ax + qx * qx / 2)
    vz = Omega / 2 * 1000 * np.sqrt(az)
    return ax, qx

def mathieuparams(mass_o,Vrf,Omega_o,Uend):
    ax =  -4*kappa * Uend * e_c * charge / (mass_o * Omega_o**2 * z0_o**2)
    qx = 4 * Vrf * e_c * charge / (mass_o * Omega_o**2 * r0_o**2)   
    return ax, qx

plt.close('all')

# plot
fig, axs = plt.subplots()

plt.ylim([-0.3, 0.3])
plt.xlabel('qx')
plt.ylabel('ax')



# stability diagram from Gregor
# Assumes a standard quadrupole with infinite hyperboloid RF-electrodes (ideal trap) and non interfering end-caps!
ranges = np.linspace(0,0.90891,100)
stab1 = np.array([-x**2/2 + 7*x**4/128 + 29*x**6/2304 for x in ranges])
stab2 = -stab1

stab3 = np.array([1-x-x**2/8+ x**3/64 + x**4/1536 for x in ranges])
stab4 = -stab3



#fill area
higher = np.maximum(stab1, stab4)
lower = np.minimum(stab2, stab3)

plt.plot(ranges, lower, c = 'b')
plt.plot(ranges, higher, c = 'b')
plt.plot([ranges[0], ranges[-1]], [0,0], c = 'b')

plt.fill_between(ranges, lower, higher, alpha = 0.4)

# calculate ax, qx according to Arezoos thesis
axold , qxold = mathieuparams(mass_o,Vrf,Omega_o,Uend)
print(axold, qxold)
plt.scatter([qxold], [axold], c = 'r')
plt.grid()
plt.show()

# qx, qx from Gregors program
#ax, qx = mathieuparams_greg(mass,UDC,Vrf,Omega,Uend)
#print(ax, qx)
#plt.scatter([qx], [ax], c = 'b')

#rfdat = np.array([mathieuparams(mass,UDC,x,Omega,Uend) for x in np.linspace(100,1000,100)]).T
#plt.plot(rfdat[1],rfdat[0], c = 'red', label = 'rf')
#
#Omegadat =  np.array([mathieuparams(mass,UDC,Vrf,x,Uend) for x in np.linspace(2*np.pi*0.8,2*np.pi*5,100)]).T
#plt.plot(Omegadat[1],Omegadat[0], c = 'blue', label = '$\Omega$')

#plt.legend()