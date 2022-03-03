#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:45:12 2019

@author: thomas
"""
import random
from numba import jit
import numpy as np
import matplotlib.pyplot as plt


@jit(nopython=True)
def get_index_and_delta(pos,lab):
    diff_old = 1.
    ind_low  = len(lab)
    ind_high = 0.
    flag_end = False

    for i in range(len(lab)):
        if pos-lab[i]<0 and not flag_end:
            ind_high = i-0
            ind_low  = i-1
            flag_end = True

    if not flag_end:
        ind_high = int(len(lab)-1)
        ind_low  = int(ind_high-1)

    delta_pos = (pos-lab[int(ind_low)])/(lab[int(ind_high)]-lab[int(ind_low)])
    return(ind_low,ind_high,delta_pos)

@jit(nopython=True)
def final_tri(acc,scale_x,scale_y,scale_z):
    c00 = acc[0]*(1-scale_x)+acc[4]*scale_x
    c01 = acc[1]*(1-scale_x)+acc[5]*scale_x
    c10 = acc[2]*(1-scale_x)+acc[6]*scale_x
    c11 = acc[3]*(1-scale_x)+acc[7]*scale_x
    c0  = c00*(1-scale_y)+c10*scale_y
    c1  = c01*(1-scale_y)+c11*scale_y
    c   = c0*(1-scale_z)+c1*scale_z
    return(c)


@jit(nopython=True)
def accelerate(ax_3d,ay_3d,az_3d,ax_3d_n,ay_3d_n,az_3d_n,x_lab,x_lab_r,y_lab,z_lab,d,dec_input,phase,fish,cycle,sw_cy):
    acc_inter = np.zeros(np.shape(d))

    switching_auto = False
    if switching_auto:
        #####---Here I find the switching cycle manually.
        switching_cycle = (np.mean(d[0,:])-dec_input-phase)//(fish/2)
        ppos_x = d[0,:]-dec_input-(switching_cycle+1)*fish/2
        switching_cycle_particles = switching_cycle*np.ones(np.shape(d[0,:]))

    else:
        switching_cycle = sw_cy
        if switching_cycle%2 == 0 and sw_cy>=0:
            switching_cycle_particles = (d[0,:]-dec_input)//(fish)
            ppos_x = d[0,:]-dec_input-switching_cycle_particles*fish-fish/2

        if switching_cycle%2 == 1 and sw_cy>=0:
            switching_cycle_particles = (d[0,:]-dec_input-fish/2)//(fish)
            ppos_x = d[0,:]-dec_input-switching_cycle_particles*fish-fish/2-fish/2


    ###--- Bug finder---###
    #min_pppos_x = np.min(ppos_x)
    #max_pppos_x = np.max(ppos_x)
    #diff = max_pppos_x-min_pppos_x
    #print(max_pppos_x)

    acc_x = np.zeros((8))
    acc_y = np.zeros((8))
    acc_z = np.zeros((8))

    if switching_cycle%2 == 0 and switching_cycle<121 and switching_cycle>=0:
        cycle = True
        ###---- Loop over all positions
        for i in range(0,len(ppos_x)):
            #if switching_auto:
            if switching_cycle_particles[i]<121 and switching_cycle_particles[i]>=0:
            #if 1 == 1:
                #####---------------------    TRUE AND POS ----------------------
                if ppos_x[i]>=0:
                    x_ind_l,x_ind_h,scale_x = get_index_and_delta(ppos_x[i],x_lab)
                    y_ind_l,y_ind_h,scale_y = get_index_and_delta(d[1,i],y_lab)
                    z_ind_l,z_ind_h,scale_z = get_index_and_delta(d[2,i],z_lab)

                    cnt=0
                    for x in range(int(x_ind_l),int(x_ind_h+1)):
                        for y in range(int(y_ind_l),int(y_ind_h+1)):
                            for z in range(int(z_ind_l),int(z_ind_h+1)):
                                acc_x[cnt] = ax_3d[x,y,z]
                                acc_y[cnt] = ay_3d[x,y,z]
                                acc_z[cnt] = az_3d[x,y,z]
                                cnt+=1

                    acc_inter[0,i] = final_tri(acc_x,scale_x,scale_y,scale_z)
                    acc_inter[1,i] = final_tri(acc_y,scale_x,scale_y,scale_z)
                    acc_inter[2,i] = final_tri(acc_z,scale_x,scale_y,scale_z)



                #####---------------------    TRUE AND Neg ----------------------
                if ppos_x[i]<0:
                    x_ind_l,x_ind_h,scale_x = get_index_and_delta(ppos_x[i],x_lab_r)
                    y_ind_l,y_ind_h,scale_y = get_index_and_delta(d[1,i],y_lab)
                    z_ind_l,z_ind_h,scale_z = get_index_and_delta(d[2,i],z_lab)

                    cnt=0
                    for x in range(int(x_ind_l),int(x_ind_h+1)):
                        for y in range(int(y_ind_l),int(y_ind_h+1)):
                            for z in range(int(z_ind_l),int(z_ind_h+1)):
                                acc_x[cnt] = ax_3d_n[x,y,z]
                                acc_y[cnt] = ay_3d_n[x,y,z]
                                acc_z[cnt] = az_3d_n[x,y,z]
                                cnt+=1

                    acc_inter[0,i] = final_tri(acc_x,scale_x,scale_y,scale_z)
                    acc_inter[1,i] = final_tri(acc_y,scale_x,scale_y,scale_z)
                    acc_inter[2,i] = final_tri(acc_z,scale_x,scale_y,scale_z)
            ###--- In case particles are not yet in the decc
            ###--- Or they have already left it
            else:
                acc_inter[0,i] = acc_inter[0,i]
                acc_inter[1,i] = acc_inter[1,i]
                acc_inter[2,i] = acc_inter[2,i]


    if switching_cycle%2 == 1 and switching_cycle<121 and switching_cycle>=0:
        cycle = False
        ####------ reverse the y and z positions
        y_old  = d[1,:]*1
        d[1,:] = d[2,:]*1
        d[2,:] = y_old

        for i in range(0,len(ppos_x)):
            #if 1 == 1:
            if switching_cycle_particles[i]<121 and switching_cycle_particles[i]>=0:

                #####---------------------    False AND POS ----------------------
                if ppos_x[i]>=0:
                    x_ind_l,x_ind_h,scale_x = get_index_and_delta(ppos_x[i],x_lab)
                    y_ind_l,y_ind_h,scale_y = get_index_and_delta(d[1,i],y_lab)
                    z_ind_l,z_ind_h,scale_z = get_index_and_delta(d[2,i],z_lab)
                    cnt=0

                    for x in range(int(x_ind_l),int(x_ind_h+1)):
                        for y in range(int(y_ind_l),int(y_ind_h+1)):
                            for z in range(int(z_ind_l),int(z_ind_h+1)):
                                acc_x[cnt] = ax_3d[x,y,z]
                                acc_y[cnt] = ay_3d[x,y,z]
                                acc_z[cnt] = az_3d[x,y,z]
                                cnt+=1

                    acc_inter[0,i] = final_tri(acc_x,scale_x,scale_z,scale_y)
                    acc_inter[1,i] = final_tri(acc_y,scale_x,scale_z,scale_y)
                    acc_inter[2,i] = final_tri(acc_z,scale_x,scale_z,scale_y)


                if ppos_x[i]<0:
                    #####---------------------    False AND Neg ----------------------
                    x_ind_l,x_ind_h,scale_x = get_index_and_delta(ppos_x[i],x_lab_r)
                    y_ind_l,y_ind_h,scale_z = get_index_and_delta(d[1,i],y_lab)
                    z_ind_l,z_ind_h,scale_y = get_index_and_delta(d[2,i],z_lab)

                    cnt=0
                    for x in range(int(x_ind_l),int(x_ind_h+1)):
                        for y in range(int(y_ind_l),int(y_ind_h+1)):
                            for z in range(int(z_ind_l),int(z_ind_h+1)):
                                acc_x[cnt] = ax_3d_n[x,y,z]
                                acc_y[cnt] = ay_3d_n[x,y,z]
                                acc_z[cnt] = az_3d_n[x,y,z]
                                cnt+=1

                    acc_inter[0,i] = final_tri(acc_x,scale_x,scale_z,scale_y)
                    acc_inter[1,i] = final_tri(acc_y,scale_x,scale_z,scale_y)
                    acc_inter[2,i] = final_tri(acc_z,scale_x,scale_z,scale_y)
            ###--- In case particles are not yet in the decc
            ###--- Or they have already left it
            else:
                acc_inter[0,i] = 0
                acc_inter[1,i] = 0
                acc_inter[2,i] = 0


        ####------ reverse the y and z positions
        y_old  = d[1,:]*1
        d[1,:] = d[2,:]*1
        d[2,:] = y_old*1

        old_acc = acc_inter[2,:]*1
        acc_inter[2,:] = acc_inter[1,:]*1
        acc_inter[1,:] = old_acc*1

    return(acc_inter,cycle,switching_cycle)

if __name__ == '__main__':
    plt.close('all')
    ####--- Stuff which should only be once calculated...
    ax_3d = np.load('10kVNormal/ax.npy')
    ay_3d = np.load('10kVNormal/ay.npy')
    az_3d = np.load('10kVNormal/az.npy')


    x_lab     = np.linspace(-1e-3,6.5e-3,151)
    x_lab_r   = np.linspace(-6.5e-3,1e-3,151)   ####--- NOTE I HAVE TO INVERSE THE SIGN OF ACCELERATION
    y_lab     = np.linspace(-2e-3,2e-3,41)
    z_lab     = np.linspace(-2e-3,2e-3,41)


    ###--- Geometrie
    dec_input = 237.3e-3
    fish      = 11.e-3
    phase     = fish/4    ####phi 0

    ####--- Create dummy particles
    d      = np.empty((3,10))
    for i in range(len(d[0,])):
        d[0,i] = dec_input+np.random.normal(loc=5.5e-3*2+0.2e-3, scale=1e-3, size=1)
        d[1,i] = np.random.random_sample(1)*1e-3*random.choice((-1, 1))
        d[2,i] = np.random.random_sample(1)*1e-3*random.choice((-1, 1))

    acc_inter,cycle = accelerate(ax_3d,ay_3d,az_3d,x_lab,x_lab_r,y_lab,z_lab,d,dec_input,phase,fish,True)
