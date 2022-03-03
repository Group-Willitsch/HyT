#Code written by Dr. Thomas Kierspel
#Problems? Please contact Dr. Claudio von Planta
#coding: utf-8

import numpy as np
import time
import matplotlib.pyplot as plt
import sys
from numba import jit
from numba import jitclass
import scipy.constants as sci_con
from find_gradient import get_fields
from find_gradient import calculate_acc
from find_gradient import get_fields_interpolated
from find_gradient_ast_2 import accelerate


#####-------- TO DOs:
#####-------- Fix everything which is not working...
#####-------- Get a propper TOF (at position x)
#####-------- Check the fields in transverse direction (Hopefully working now for ast)
#####-------- Replace Euler/Verlet by Runge-Kutta-Fehlberg
#####-------- Check that only within the decelerator acceleration is calculated
#####-------- Solve the problem with extended beam in longditudinal direction
#####-------- Remove shit from memory aka not-used particles
#####-------- Improve calculation speed by moving to jit-calculation, variable step size etc



def bug():
    plt.show()
    sys.exit()

def plot_particles(fgc,particles):
    fgc+=1
    gs = plt.GridSpec(3,3)
    fig = plt.figure(fgc)
    ax1 = fig.add_subplot(gs[0,0],aspect='auto')
    plt.hist2d(particles[0,:].ravel(),particles[1,:].ravel(),bins=40)#,range=[[-pos_range,pos_range],[-pos_range,pos_range]])

    ax1 = fig.add_subplot(gs[0,1],aspect='auto')
    plt.hist2d(particles[0,:].ravel(),particles[1,:].ravel(),bins=40)#,range=[[-pos_range,pos_range],[-pos_range,pos_range]])
    ax1 = fig.add_subplot(gs[0,2],aspect='auto')
    plt.hist2d(particles[1,:].ravel(),particles[2,:].ravel(),bins=40)#,range=[[-pos_range,pos_range],[-pos_range,pos_range]])

    ax1 = fig.add_subplot(gs[1,0],aspect='auto')
    plt.hist2d(particles[3,:].ravel(),particles[4,:].ravel(),bins=40)#,range=[[-3e-3,3e-3],[-3e-3,3e-3]])
    ax1 = fig.add_subplot(gs[1,1],aspect='auto')
    plt.hist2d(particles[3,:].ravel(),particles[5,:].ravel(),bins=40)#,range=[[-3e-3,3e-3],[-3e-3,3e-3]])
    ax1 = fig.add_subplot(gs[1,2],aspect='auto')
    plt.hist2d(particles[4,:].ravel(),particles[5,:].ravel(),bins=40)#,range=[[-3e-3,3e-3],[-3e-3,3e-3]])

    ax1 = fig.add_subplot(gs[2,0],aspect='auto')
    plt.hist2d(particles[0,:].ravel(),particles[3,:].ravel(),bins=40)#,range=[[-3e-3,3e-3],[-3e-3,3e-3]])
    ax1 = fig.add_subplot(gs[2,1],aspect='auto')
    plt.hist2d(particles[1,:].ravel(),particles[4,:].ravel(),bins=40)#,range=[[-3e-3,3e-3],[-3e-3,3e-3]])
    ax1 = fig.add_subplot(gs[2,2],aspect='auto')
    plt.hist2d(particles[2,:].ravel(),particles[5,:].ravel(),bins=40)#,range=[[-3e-3,3e-3],[-3e-3,3e-3]])
    return(fgc)


class source:
    def __init__(self):
        self.N      = int(3e2)    ###--- Number of paritcles, 2e7 enough for imac, 1e8 particles ≈ 10...15 GB ram...b4 they fly... 5e8 resulted in > 70 GB

        self.mean_x = 425.        ###---  m/s
        self.mean_y = 0.          ###---  m/s
        self.mean_z = 0.          ###---  m/s
        self.dv_x   = 10. #self.mean_x/2.3548*0.13            ###--- FWHM Velocity spread in m/s
        self.dv_y   = 0. #self.mean_x/2.3548/2*0.13/2        ###--- Velocity spread in m/s
        self.dv_z   = 0. #self.mean_x/2.3548/2*0.13/2        ###--- Velocity spread in m/s

        ### --- Virtual source ---
        self.virt_x =  1.e-3                 ####---mean position X
        self.virt_y =  0. #0.5e-3         ####---mean position Y
        self.virt_z =  0. #0.5e-3         ####---mean position Z
        self.dvirt_x = 0. #11.5e-3/2.3548/2      ####--- "sigma" X if gauss is chosen
        self.dvirt_y = 0.#1e-3 #0.15e-3/2.3548/2      ####--- "sigma" Y if gauss is chosen
        self.dvirt_z = 0.#1e-3 #0.15e-3/2.3548/2      ####--- "sigma" Z if gauss is chosen

        ###--- Position vectors for paritcles ---###
        self.X    = 0.  ### Position X lab-frame
        self.Y    = 0.  ### Position Y lab-frame
        self.Z    = 0.  ### Position Z lab-frame
        self.v_x  = 0.  ### Velocity X lab-frame
        self.v_y  = 0.  ### Velocity Y lab-frame
        self.v_z  = 0.  ### Velocity Z lab-frame

    def gauss(self):
        self.X = np.random.normal(self.virt_x,self.dvirt_x,self.N)
        self.Y = np.random.normal(self.virt_y,self.dvirt_y,self.N)
        self.Z = np.random.normal(self.virt_z,self.dvirt_z,self.N)
        self.v_x = np.random.normal(self.mean_x,self.dv_x,self.N)
        self.v_y = np.random.normal(self.mean_y,self.dv_y,self.N)
        self.v_z = np.random.normal(self.mean_z,self.dv_z,self.N)
        return(np.array((self.X,self.Y,self.Z,self.v_x,self.v_y,self.v_z)))


##### jitclass: https://numba.pydata.org/numba-doc/dev/user/jitclass.html#numba.jitclass
#@jitclass()
class integrator:
    def __init__(self,particles,apertures_geom):
        self.time = 0.0
        self.dt   = 5e-7             #5e-9  1e-6

        self.particles = particles
        self.pos  = self.particles[:3,:]           ### Position X lab-frame
        self.vel  = self.particles[3:,:]
        self.acc  = np.zeros(np.shape(self.vel))
        self.apertures = apertures_geom

        ####----- Flags, and local varialbes ----
        self.before_dec = True
        self.after_dec  = False
        self.no_fore    = 0                        ### Used as multplier for free-flight
        self.free_flight = True
        self.mesh       = np.linspace(-3,3,101)
        self.interpolation_fields = get_fields_interpolated()
        self.ax_3d,self.ay_3d,self.az_3d = get_fields()#,self.ax_3d_n,self.ay_3d_n,self.az_3d_n = get_fields()
        self.ax_3d_n = -self.ax_3d[::-1,:,:]
        self.ay_3d_n = self.ay_3d[::-1,:,:]
        self.az_3d_n = self.az_3d[::-1,:,:]
        self.cycle      = 0
        self.switching_cycle = -999
        self.list_switching_cycle = -999
        self.switching_sequence = 0
        self.particle_position_x  = 0           ####Used for debugging
        self.particle_position_y  = 0           ####Used for debugging
        self.particle_position_z  = 0           ####Used for debugging
        self.particle_velocity_x  = 0           ####Used for debugging
        self.particle_velocity_y  = 0           ####Used for debugging
        self.particle_velocity_z  = 0           ####Used for debugging


        self.particle_acc       = 0           ####Used for debugging
        self.x_lab     = np.linspace(-1e-3,6.5e-3,151)
        self.x_lab_r   = np.linspace(-6.5e-3,1e-3,151)   ####--- NOTE I HAVE TO INVERSE THE SIGN OF ACCELERATION
        self.y_lab     = np.linspace(-2e-3,2e-3,41)
        self.z_lab     = np.linspace(-2e-3,2e-3,41)


        ####---- Calculating switching sequence OR trajectories with given sequence
        self.sequence_calculation = False                #####If true calculate switching sequence
        self.save_sequence_for_experiment = False       #####If true save the calculated sequence
        self.read_switching_cycle_list    = np.load('switching/switching_cycle_sim_10kV_5e-07_0.0deg.npy')  ####Read list indipendent of self.sequence_calculation

        self.fish      = 11.e-3
        self.phase_deg = 0.0 #55.47      #55.468   #55.468#50.468#55.468
        self.phase     = self.fish/4+self.phase_deg*self.fish/(4*90)    ####phi


    def eulerint_force(self):
        ###--- Check if particles are still b4 the decelerator
        if self.before_dec:
            self.no_force = self.pos[0,]>self.apertures.dec_input     ####Check if any particel is already at the decelerator input
            if np.sum(self.no_force)==len(self.pos[0,]):              ####If non
                self.before_dec = False

            ###--- Used for debugging
            self.particle_position_x  = np.append(self.particle_position_x,self.pos[0,0])
            self.particle_position_y  = np.append(self.particle_position_y,self.pos[1,0])
            self.particle_position_z  = np.append(self.particle_position_z,self.pos[2,0])
            self.switching_sequence   = np.append(self.switching_sequence,self.cycle)
            self.list_switching_cycle = np.append(self.list_switching_cycle,-999)

            self.particle_velocity_x  = np.append(self.particle_velocity_x,self.vel[0,0])
            self.particle_velocity_y  = np.append(self.particle_velocity_y,self.vel[1,0])
            self.particle_velocity_z  = np.append(self.particle_velocity_z,self.vel[2,0])
            self.particle_acc       = np.append(self.particle_acc,0)

            self.pos += self.vel * self.dt
            self.time+=1

        ###--- Check if particles already left the decelerator
        if self.after_dec:
            ###--- Used for debugging
            self.particle_position_x  = np.append(self.particle_position_x,self.pos[0,0])
            self.particle_position_y  = np.append(self.particle_position_y,self.pos[1,0])
            self.particle_position_z  = np.append(self.particle_position_z,self.pos[2,0])
            self.switching_sequence   = np.append(self.switching_sequence,0)
            self.list_switching_cycle = np.append(self.list_switching_cycle,999)
            self.particle_velocity_x  = np.append(self.particle_velocity_x,self.vel[0,0])
            self.particle_velocity_y  = np.append(self.particle_velocity_y,self.vel[1,0])
            self.particle_velocity_z  = np.append(self.particle_velocity_z,self.vel[2,0])
            self.particle_acc       = np.append(self.particle_acc,0)

            self.pos += self.vel * self.dt
            self.time+=1

        ###--- Check if particles already left the decelerator
        if self.before_dec == False:
            if self.after_dec == False:
                ###--- Remove particels which are not within the Z-Y boarders within the decelerator
                af_dec = self.particles[0,:]<self.apertures.dec_output
                self.particles = self.particles[:,(((abs(self.particles[1,:])>self.apertures.dec_width).astype(int)+af_dec.astype(int))<2)]
                b4_dec = self.particles[0,:]>self.apertures.dec_input
                af_dec = self.particles[0,:]<self.apertures.dec_output
                self.particles = self.particles[:,(((abs(self.particles[2,:])>self.apertures.dec_width).astype(int)+af_dec.astype(int))<2)]
                self.pos = self.particles[:3,:]
                self.vel = self.particles[3:,:]

                verlet = True
                if verlet:
                    ###From https://en.wikipedia.org/wiki/Verlet_integration "Eliminating the half-step velocity, this algorithm may be shortened to:"
                    acc_old = self.acc
                    self.pos += self.vel*self.dt + 0.5*acc_old*self.dt*self.dt
                    self.acc,self.cycle = calculate_acc(self.pos,self.apertures.dec_input,self.phase,self.fish ,self.cycle,self.interpolation_fields,self.sequence_calculation,self.read_switching_cycle_list[1][int(self.time)])
                    #self.acc = self.acc*(1-electrode*0.000335)    ####--- Correction factor for the possibility that the electrodes have different width/distances
                    #print(electrode+1)
                    self.vel += 0.5*(acc_old+self.acc)*self.dt


                else:   ###Euler
                    print('Old stuff, should not use it anymore')
                    bug()
                self.time+=1


                ###--- Used for debugging
                self.switching_sequence = np.append(self.switching_sequence,self.cycle)
                self.list_switching_cycle = np.append(self.list_switching_cycle,self.switching_cycle)
                self.particle_position_x  = np.append(self.particle_position_x,self.pos[0,0])
                self.particle_position_y  = np.append(self.particle_position_y,self.pos[1,0])
                self.particle_position_z  = np.append(self.particle_position_z,self.pos[2,0])
                self.particle_velocity_x  = np.append(self.particle_velocity_x,self.vel[0,0])
                self.particle_velocity_y  = np.append(self.particle_velocity_y,self.vel[1,0])
                self.particle_velocity_z  = np.append(self.particle_velocity_z,self.vel[2,0])
                self.particle_acc       = np.append(self.particle_acc,self.acc[0,0])


                self.no_force = self.pos[0,:]<self.apertures.dec_output
                if np.sum(self.no_force)==0:
                    self.after_dec = True

        return(np.append(self.pos,self.vel,axis=0))


class apertures:
    def __init__(self,particles):
        self.pos  = particles[:3,:]  ### Position X lab-frame
        self.vel  = particles[3:,:]  ### Position Y lab-frame

        self.skimmer_1_dia = 3e-3    ### circle
        self.skimmer_1_pos = 100e-3

        self.skimmer_2_dia = 2e-3    ### circle
        self.skimmer_2_pos = 210e-3

        self.dec_input     = 237.3e-3   ### square
        self.dec_width     = 2e-3       ### 2x2 mm
        self.dec_output    = self.dec_input+123*11e-3/2


    def skimmer_source_1_2_dec(self):
        ###---Calculate the time when particles arrive at the apertures
        t_skimmer_1 = self.pos[0,:]+self.skimmer_1_pos/self.vel[0,:]
        t_skimmer_2 = self.pos[0,:]+self.skimmer_2_pos/self.vel[0,:]
        t_dec       = self.pos[0,:]+self.dec_input/self.vel[0,:]

        ###---Calculate the transversal position of the particles at the apertures
        rad_skimmer_1 = np.sqrt((self.pos[1,:]+self.vel[1,:]*t_skimmer_1)**2+(self.pos[2,:]+self.vel[2,:]*t_skimmer_1)**2)
        rad_skimmer_2 = np.sqrt((self.pos[1,:]+self.vel[1,:]*t_skimmer_2)**2+(self.pos[2,:]+self.vel[2,:]*t_skimmer_2)**2)
        dec_y         = self.pos[1,:]+self.vel[1,:]*t_dec
        dec_z         = self.pos[2,:]+self.vel[2,:]*t_dec

        gate_skimmer_1 = (rad_skimmer_1<self.skimmer_1_dia/2.).astype(int)
        gate_skimmer_2 = (rad_skimmer_2<self.skimmer_2_dia/2.).astype(int)
        gate_dec       = (abs(dec_y)<self.dec_width/2.).astype(int)*(abs(dec_z)<self.dec_width/2.).astype(int)
        return(gate_skimmer_1,gate_skimmer_2,gate_dec,gate_skimmer_1*gate_skimmer_2*gate_dec)





if __name__ == '__main__':
    plt.close('all')
    t1 = time.time()                                   ####--- Save the time to check how long things take  ####--- Start a timer
    fgc=0

    source = source()                                  ####---- Initiallize the source
    particles = source.gauss()                         ####---- Get the particles inital values (X,Y,Z,v_x,v_y,v_z)
    apertures_1   = apertures(particles)               ####---- Initiallize the apertures
    gates = apertures_1.skimmer_source_1_2_dec()       ####---- Get the gates for the apertures

                                                       ####---- U have to decide here!? :)

    ###----- Remove all particles which don't make it to the Decelerator
    particles = particles[:,gates[3]==1]
    gates = 0

    ###--- Particles That make it into the decelarator)
    source.X = particles[0,:]
    source.Y = particles[1,:]
    source.Z = particles[2,:]
    source.v_x = particles[3,:]
    source.v_y = particles[4,:]
    source.v_z = particles[5,:]

    apertures_1   = apertures(particles)
    trajectories  = integrator(particles,apertures_1)    ####---- Initiallize the integrator


    ###--- Prepare the propagation plot
    bins = np.linspace(-0.1,2,100000)
    plt.figure(fgc)
    i = 0

    ###---- Calculate the trajectories ---
    old_mean = -5                                        ###---- Local variable to start the loop
    import operator
    t1 = time.time()
    pos_long = 0
    position_h = 0
    position_v = 0

    ###---- Calculate the switching sequence ---
    ###---- Calculate the switching sequence ---
    ###---- Calculate the switching sequence ---
    if trajectories.sequence_calculation:
        steps = 0
        print('Calculating your switching cycles! U r whalecome!')
        while(particles[0,0]<apertures_1.dec_output+11.5e-3):
            i+=1
            old_mean = np.mean(particles[0,:])
            particles = trajectories.eulerint_force()


            ####---- Stuff for testing
            if i%5 == 0:
                #plt.hist(particles[0,:],bins=bins,histtype='step',color='black')
                pos_long = np.append(pos_long,particles[0,0])
                position_h = np.append(position_h,particles[1,0])
                position_v = np.append(position_v,particles[2,0])


        time_of_flight = np.linspace(0,trajectories.dt*len(trajectories.switching_sequence),len(trajectories.switching_sequence))
        cnt_rising = 0
        cnt_falling = 0
        time_rise = []
        time_fall = []
        for i in range(len(time_of_flight)-1):
            if trajectories.switching_sequence[i] < trajectories.switching_sequence[i+1]:    ###rising edge
                cnt_rising +=1
                time_rise = np.append(time_rise,time_of_flight[i])
            if trajectories.switching_sequence[i] > trajectories.switching_sequence[i+1]:    ###falling edge
                time_fall = np.append(time_fall,time_of_flight[i])
                cnt_falling +=1


        duration_rising  = abs(time_rise-time_fall[1:])   ####--We start the sequence with a falling edge
        duration_falling = abs(time_fall[:-1]-time_rise)

        if trajectories.save_sequence_for_experiment:
            t0 = time_fall[0]
            time_rise = time_rise-t0
            time_fall = time_fall-t0
            if trajectories.dt > 5e-9:
                print("!!!!!WARNING WARNING YOUR TIME STEPS FOR INTEGRATRION ARE TOOOOOOO LONG")

            print('Save switching sequence for experiment')
            np.save('/Users/thomas/ownCloud/Lab/Lab315/Lab315/Software/hytrap_bg/input/input_test_py_code/v_time_'+str(int(source.mean_x))+'_'+str(trajectories.phase_deg),time_fall)
            np.save('/Users/thomas/ownCloud/Lab/Lab315/Lab315/Software/hytrap_bg/input/input_test_py_code/v_duration_'+str(int(source.mean_x))+'_'+str(trajectories.phase_deg),duration_falling)
            np.save('/Users/thomas/ownCloud/Lab/Lab315/Lab315/Software/hytrap_bg/input/input_test_py_code/h_time_'+str(int(source.mean_x))+'_'+str(trajectories.phase_deg),time_rise)
            np.save('/Users/thomas/ownCloud/Lab/Lab315/Lab315/Software/hytrap_bg/input/input_test_py_code/h_duration_'+str(int(source.mean_x))+'_'+str(trajectories.phase_deg),duration_rising)


        ###count switching-cycles:
        get_cycles = np.zeros(np.shape(trajectories.switching_sequence))
        cnt = 0
        for i in range(len(trajectories.switching_sequence)-1):
            if trajectories.switching_sequence[i] != trajectories.switching_sequence[i+1]:
                cnt+=1
                get_cycles[i] = cnt
            else:
                get_cycles[i] = get_cycles[i-1]
        to_save_switching = trajectories.switching_sequence,get_cycles
        np.save('switching/switching_cycle_sim_10kV_'+str(trajectories.dt)+'_'+str(trajectories.phase_deg)+'deg',to_save_switching)

        #####----Plots and prints----
        plt.figure(1)
        plt.plot(time_of_flight,trajectories.switching_sequence)

        plt.figure(2)
        plt.plot(trajectories.particle_position_y)
        plt.plot(trajectories.particle_position_z)

        print('Velo: '+str(np.mean(particles[3,:])))
        t2 = time.time()
        print(str(t2-t1)+' s')
        bug()


    ###---- Calculate the molecules TOF for a given switching sequence ---
    ###---- Calculate the molecules TOF for a given switching sequence ---
    ###---- Calculate the molecules TOF for a given switching sequence ---
    else:
        print('Read from switching cycle list')
        print('Check that the time steps are correct!!!!!!!!!!!')
        while(trajectories.time<1e6):

            ###--- Save the old position of the particles
            old_mean = np.mean(particles[0,:])

            if trajectories.time == len(trajectories.read_switching_cycle_list[0])-1:
                print('this is the end')
                break

            ###--- Calculate the next step
            particles = trajectories.eulerint_force()

            ####---- Stuff for testing
            if i%5 == 0:
                #plt.hist(particles[0,:],bins=bins,histtype='step',color='black')
                pos_long = np.append(pos_long,particles[0,0])
                position_h = np.append(position_h,particles[1,0])
                position_v = np.append(position_v,particles[2,0])
            i+=1



    #print(np.shape(particles))
    t2 = time.time()


    #print('Velo: '+str(np.mean(particles[3,:])))
    fgc = plot_particles(fgc,particles)


    ###count switching-cycles:
    cnt_sw_cy = np.zeros(np.shape(trajectories.switching_sequence))
    for i in range(len(trajectories.switching_sequence)-1):
        if trajectories.switching_sequence[i] != trajectories.switching_sequence[i+1]:
            cnt_sw_cy[i] = True


    fgc+=1
    plt.figure(fgc)
    plt.hist2d(particles[1,:].ravel(),particles[2,:].ravel(),bins=40,range=[[-3e-3,3e-3],[-3e-3,3e-3]])

    fgc+=1
    plt.figure(fgc)
    #plt.plot(trajectories.particle_position,trajectories.switching_sequence)
    plt.plot(trajectories.particle_position_x,trajectories.particle_velocity_x)
    fgc+=1
    plt.figure(fgc)
    plt.plot(trajectories.particle_position_x,trajectories.particle_position_y)
    plt.plot(trajectories.particle_position_x,trajectories.particle_position_z)

    fgc+=1
    plt.figure(fgc)
    plt.plot(pos_long,position_h)
    plt.plot(pos_long,position_v)

    fgc+=1
    plt.figure(fgc)
    plt.hist2d(particles[1,:],particles[2,:])
    #plt.hist(particles[0,:])

    fgc+=1
    plt.figure(fgc)
    plt.hist(particles[0,:])
    # print('Particle starting parameter')
    # print('#'+str(trajectories.particle_position_x[1]))           ####Used for debugging
    # print('#'+str(trajectories.particle_position_y[1]))              ####Used for debugging
    # print('#'+str(trajectories.particle_position_z[1]))              ####Used for debugging
    # print('#'+str(trajectories.particle_velocity_x[1]))              ####Used for debugging
    # print('#'+str(trajectories.particle_velocity_y[1]))               ####Used for debugging
    # print('#'+str(trajectories.particle_velocity_z[1]))              ####Used for debugging


    plt.show()
    t2 = time.time()
    print(t2-t1)
