#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:31:42 2019

@author: T. Kierpsel copied from claudiovonplanta file
"""
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import numpy as np
#from numba import jit

plt.close('all')

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
                print(index)

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

def Bxmodel(x,y,z):
    ret = (ax000 + z*ax001 + z**2*ax002 + y*ax010 + y*z*ax011 +
        y*z**2*ax012 + y**2*ax020 + y**2*z*ax021 + y**2*z**2*ax022 +
        x*ax100 + x*z*ax101 + x*z**2*ax102 + x*y*ax110 +
        x*y*z*ax111 + x*y*z**2*ax112 + x*y**2*ax120 + x*y**2*z*ax121 +
        x*y**2*z**2*ax122 + x**2*ax200 + x**2*z*ax201 + x**2*z**2*ax202 +
        x**2*y*ax210 + x**2*y*z*ax211 + x**2*y*z**2*ax212 + x**2*y**2*ax220 +
        x**2*y**2*z*ax221 + x**2*y**2*z**2*ax222)
    return ret

def Bymodel(x,y,z):
    ret = (ay000 + z*ay001 + z**2*ay002 + y*ay010 + y*z*ay011 +
        y*z**2*ay012 + y**2*ay020 + y**2*z*ay021 + y**2*z**2*ay022 +
        x*ay100 + x*z*ay101 + x*z**2*ay102 + x*y*ay110 +
        x*y*z*ay111 + x*y*z**2*ay112 + x*y**2*ay120 + x*y**2*z*ay121 +
        x*y**2*z**2*ay122 + x**2*ay200 + x**2*z*ay201 + x**2*z**2*ay202 +
        x**2*y*ay210 + x**2*y*z*ay211 + x**2*y*z**2*ay212 + x**2*y**2*ay220 +
        x**2*y**2*z*ay221 + x**2*y**2*z**2*ay222)
    return ret

def Bzmodel(x,y,z):
    ret = (az000 + z*az001 + z**2*az002 + y*az010 + y*z*az011 +
        y*z**2*az012 + y**2*az020 + y**2*z*az021 + y**2*z**2*az022 +
        x*az100 + x*z*az101 + x*z**2*az102 + x*y*az110 +
        x*y*z*az111 + x*y*z**2*az112 + x*y**2*az120 + x*y**2*z*az121 +
        x*y**2*z**2*az122 + x**2*az200 + x**2*z*az201 + x**2*z**2*az202 +
        x**2*y*az210 + x**2*y*z*az211 + x**2*y*z**2*az212 + x**2*y**2*az220 +
        x**2*y**2*z*az221 + x**2*y**2*z**2*az222)
    return ret

def Bmodel(x,y,z):
    return np.sqrt(Bxmodel(x,y,z)**2 + Bymodel(x,y,z)**2 + Bzmodel(x,y,z)**2)


#if __name__ == "__main__":



plt.tight_layout()
plt.show()
