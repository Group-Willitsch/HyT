import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

plt.close('all')

radius  = 0.003         #### Use radius in cm
energy  = 0.5e-3         #### Use SI units
pulse_d = 10e-9              #### Use SI units
power  = energy/pulse_d

sigma = radius

ft    = True   ###Flat Top
gauss = True

####--- Flat Top---###
if ft:
    area = radius**2*np.pi   ###100*100 damit nicht m^2 aber cm^2
    intensity = power/area
    print('Peak-Intensity for flat-Top is: '+str(np.round(intensity,4))+' W/cm^2')

####--- Gauss ---###
if gauss:
    ##https://en.wikipedia.org/wiki/Gaussian_function
        #E = Integral I dr
        #E = I0 * ( (2pi)^3/2) * sigma x sigma y sigma t)
        #Io = E / ( (2pi)^3/2) * sigma x sigma y sigma t)

    I0 = power/(2*np.pi*sigma**2)
    print('Peak-Intensity for for Gauss is: '+str(np.round(I0,4))+' W/cm^2')


###--- Plot something ^^
x = np.linspace(-0.5,0.5,100001)

y_square = ((x>=-radius)*(x<=radius)).astype(int)
y_sq_n   = y_square/np.sum(y_square)


y_gauss  = 1*np.exp(-x**2/(2*sigma**2))
y_gauss_n = y_gauss/np.sum(y_gauss)

y_sq_n   = y_sq_n/np.max(y_gauss_n)
y_gauss_n = y_gauss_n/np.max(y_gauss_n)




plt.plot(x,y_sq_n)
plt.plot(x,y_gauss_n)
plt.grid()
plt.show()

print((intensity/I0))
print((np.max(y_sq_n)))
