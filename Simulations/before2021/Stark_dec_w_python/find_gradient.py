#Code written by Thomas Kierspel
#Problems please contact Claudio von Planta
#coding: utf-8

import numpy as np
from numba import njit
from numba import jit
import matplotlib.pyplot as plt
import sys
from matplotlib.colors import LogNorm
import random
import time
from scipy.interpolate import RegularGridInterpolator  ###http://lagrange.univ-lyon1.fr/docs/scipy/0.17.1/generated/scipy.interpolate.RegularGridInterpolator.html
from find_gradient_ast import accelerate

def bug():
    plt.show()
    sys.exit()

def get_fields():
    read_raw_acc = False              ######-----True: Create 3D fields and save it, False: Read the calculated fields
    if read_raw_acc:
        ax = np.genfromtxt('10kVNormal/outax.dat',skip_header=1)
        ay = np.genfromtxt('10kVNormal/outay.dat',skip_header=1)
        az = np.genfromtxt('10kVNormal/outaz.dat',skip_header=1)

        ax_3d = np.zeros((151,41,41),dtype=np.float64)
        ay_3d = np.zeros((151,41,41),dtype=np.float64)
        az_3d = np.zeros((151,41,41),dtype=np.float64)

        t1 = time.time()

        @jit(nopython=True)
        def make_3d(oneDX,oneDY,oneDZ,threeDX,threeDY,threeDZ):
            cnt = 0
            for i in range(0,151):
                for j in range(0,41):
                    for k in range(0,41):
                        threeDX[i,j,k] = oneDX[cnt,3]
                        threeDY[i,j,k] = oneDY[cnt,3]
                        threeDZ[i,j,k] = -oneDZ[cnt,3]   ####---No clue why...:(
                        cnt+=1
            return(threeDX,threeDY,threeDZ)

        ax_3d,ay_3d,az_3d = make_3d(ax,ay,az,ax_3d,ay_3d,az_3d)
        np.save('10kVNormal/ax.npy',ax_3d)
        np.save('10kVNormal/ay.npy',ay_3d)
        np.save('10kVNormal/az.npy',az_3d)

        print('Set bla to zero')
        bug()

    else:
        ax_3d = np.load('10kVNormal/ax.npy')    ###ax = acceleration x, Arrays with a shape of [151, 41, 41]
        ay_3d = np.load('10kVNormal/ay.npy')    ###ay = acceleration y, Arrays with a shape of [151, 41, 41]
        az_3d = np.load('10kVNormal/az.npy')    ###az = acceleration z, Arrays with a shape of [151, 41, 41]
    return(ax_3d,ay_3d,az_3d)



def get_fields_interpolated():
    ax_3d,ay_3d,az_3d = get_fields()
    ####--- Local variables
    x_lab     = np.linspace(-1e-3,6.5e-3,151)
    x_lab_r   = np.linspace(-6.5e-3,1e-3,151)   ####--- NOTE I HAVE TO INVERSE THE SIGN OF ACCELERATION
    y_lab     = np.linspace(-2e-3,2e-3,41)
    z_lab     = np.linspace(-2e-3,2e-3,41)


    ####--- Make the interpolation fuctions
    x_interp_std = RegularGridInterpolator((x_lab, y_lab, z_lab), ax_3d)
    y_interp_std = RegularGridInterpolator((x_lab, y_lab, z_lab), ay_3d)
    z_interp_std = RegularGridInterpolator((x_lab, y_lab, z_lab), az_3d)

    ####--- Make the interpolation fuctions
    x_interp_rev = RegularGridInterpolator((x_lab_r, y_lab, z_lab), -ax_3d[::-1,:,:])
    y_interp_rev = RegularGridInterpolator((x_lab_r, y_lab, z_lab), ay_3d[::-1,:,:])
    z_interp_rev = RegularGridInterpolator((x_lab_r, y_lab, z_lab), az_3d[::-1,:,:])

    return(x_interp_std,y_interp_std,z_interp_std,x_interp_rev,y_interp_rev,z_interp_rev)



###--- Check which switching cycle
def calculate_acc(d,dec_input,phase,fish,cycle,interpolation_fields,sequence_calculation,switching_sequence):
###Note not working for a temporally spread molecular beam
    acc_inter = np.zeros(np.shape(d))
    if sequence_calculation:
        switching_cycle = (d[0,:]-dec_input-phase)//(fish/2)
        ppos_x = d[0,:]-dec_input-(switching_cycle+1)*fish/2

    else:
        switching_cycle = (d[0,:]-dec_input-phase)//(fish/2)
        correction = switching_cycle-switching_sequence+2       ######Not tested but hopefully corrects the calcuation if a molecular beam spreads over more than one fish
        to_correct = abs(correction%2)==1                       ###Identify the ones which are in the other switching cycle
        ppos_x = d[0,:]-dec_input-(switching_cycle+1)*fish/2#-correction*fish/2
        switching_cycle = switching_cycle + to_correct.astype(int)          ####Shift the switching cycle

        ppos_x[to_correct] = -ppos_x[to_correct]                ####Correct the x-position for the ones which are off
        y_old = d[1,:]*1                                        ####Rotate the y-z frame in the next three lines
        #d[1,to_correct] = d[2,to_correct]
        #d[2,to_correct] = y_old[to_correct]

    for i in range(len(d[0,:])):
        #print(i)
        if switching_cycle[i]%2 == 0 and switching_cycle[i]<=122 and switching_cycle[i]>=0:     #######This is defined now as horizontal
            cycle = 1                                                                 #######Thus cycle = 1 means voltage on the horizontal electrodes
            if ppos_x[i]>=0:
                acc_inter[0,i]=  interpolation_fields[0](np.array([ppos_x[i],d[1,i],d[2,i]]))
                acc_inter[1,i]=  interpolation_fields[1](np.array([ppos_x[i],d[1,i],d[2,i]]))
                acc_inter[2,i]=  interpolation_fields[2](np.array([ppos_x[i],d[1,i],d[2,i]]))

            if ppos_x[i]<0:
                acc_inter[0,i]=  interpolation_fields[3](np.array([ppos_x[i],d[1,i],d[2,i]]))
                acc_inter[1,i]=  interpolation_fields[4](np.array([ppos_x[i],d[1,i],d[2,i]]))
                acc_inter[2,i]=  interpolation_fields[5](np.array([ppos_x[i],d[1,i],d[2,i]]))


        elif switching_cycle[i]%2 == 1 and switching_cycle[i]<=121 and switching_cycle[i]>=-1: #######This is defined now as vertical
            cycle = -1
            if ppos_x[i]>=0:
                acc_inter[0,i]=  interpolation_fields[0](np.array([ppos_x[i],d[1,i],d[2,i]]))
                acc_inter[1,i]=  interpolation_fields[1](np.array([ppos_x[i],d[1,i],d[2,i]]))
                acc_inter[2,i]=  interpolation_fields[2](np.array([ppos_x[i],d[1,i],d[2,i]]))

            if ppos_x[i]<0:
                acc_inter[0,i]=  interpolation_fields[3](np.array([ppos_x[i],d[1,i],d[2,i]]))
                acc_inter[1,i]=  interpolation_fields[4](np.array([ppos_x[i],d[1,i],d[2,i]]))
                acc_inter[2,i]=  interpolation_fields[5](np.array([ppos_x[i],d[1,i],d[2,i]]))
        #else:    #######Used where my code is bad, aka the 1/2 switching stuff at the beginning and at the end of the decelarator
            #if d[0,:] < dec_input+fish:      ####We need the entrance electrodes on ground -> switching the vertical electrodes to 1
        #    print(str(d[0,:])+'\t'+str(ppos_x)+'\t'+str(switching_cycle[i])+'\t'+str(switching_cycle[i]%2 == 1))
            #x = 0

    # ppos_x[to_correct] = -ppos_x[to_correct]                ####Correct the x-position for the ones which are off
    # y_old = d[1,:]*1                                        ####Rotate the y-z frame in the next three lines
    # d[1,to_correct] = d[2,to_correct]
    # d[2,to_correct] = y_old[to_correct]


    return(acc_inter,cycle)




if __name__ == '__main__':
    plt.close('all')
    ###--- Geometrie
    dec_input = 237.3e-3
    fish      = 11.e-3
    phase     = fish/4    ####phi 0

    ####Notes:
    #### At x = 2.75 -> phi = 0
    #### At x = 5.5  -> phi = 90 deg

    ####--- Create dummy particles
    d      = np.empty((3,20))
    for i in range(len(d[0,])):
        d[0,i] = dec_input+np.random.normal(loc=5.5e-3*1+1e-3, scale=1e-3, size=1)
        d[1,i] = np.random.random_sample(1)*1e-3*random.choice((-1, 1))
        d[2,i] = np.random.random_sample(1)*1e-3*random.choice((-1, 1))

    interpolation_fields = get_fields_interpolated()
    acc_inter,cycle = calculate_acc(d,dec_input,phase,fish,True,interpolation_fields)

    plt.show()
