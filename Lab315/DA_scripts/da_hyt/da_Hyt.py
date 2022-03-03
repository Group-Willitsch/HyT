#### Developed for and by Claudantus von Plantus & T2TheS Kierspel
import matplotlib
import matplotlib.cm as cm #colormaps
import numpy as np
import h5py
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import time
import funky_fcts as ff # analysis:ignore
import glob
import os
from scipy.signal import argrelextrema
from scipy.signal import argrelmax
from numpy import shape
from matplotlib.pyplot import *
from scipy.optimize import curve_fit
#matplotlib.use('TkAgg')
'''
To do's:
- Check how the gating on the extreme values influences the error bars (artifically smaller/larger)
- Check alternative ways for photon counting
'''


#####General stuff
#####General stuff
#plt.close('all')               ###Close all open plots
fgc=0                          ###Start a figure counter
t1 = time.time()               ###Start a timer

###############-----Check config file--------#############
try:
    from lab_config_da import *
except:
    try: 
        from mac_config_da import *
    except:
        print('U need a config file...here')
        ff.bug()



###############-----Check the files_mode and the following file structure--------#############
###############-----Check the files_mode and the following file structure--------#############
###############-----Check the files_mode and the following file structure--------#############
###############-----Check the files_mode and the following file structure--------#############
if files_mode != 0 and files_mode != 1 and files_mode != 2 and files_mode != 3 and files_mode != 666:
    print('Screw you')
    ff.bug()
if da_method != 0 and da_method != 1 and da_method != 2 and da_method != 666:
    print('Screw me')
    ff.bug()


###############-----Make list of files to analyse--------#############
###############-----Check/Correct for input errors--------#############
if files_mode == 0:
    print('Files Mode: Certain files')
    file_list = []
    for i in range(len(files)):
        file_list = np.append(file_list,str(location[0])+str(dates[0]+str(files[i])+str(ending[0])))
if files_mode == 1 and start == end:
    files_mode = 0
if two_d_plots and files_mode == 1:
    for i in range(len(shots_from_2d)):
        if shots_from_2d[i]<start or shots_from_2d[i]>end:
            print('Your input does not make sense. Check \'shots_from_2d\' and your file list')
            print('I update ur list')
            shots_from_2d[i] = start
elif files_mode == 1:
    print('Files Mode: Range')
    file_list = []
    ###Check the name of all files in the dirname_pressure_dataset and take just the normal structure -- super stupid solution ^^
    h5_files = glob.glob(str(location[0])+str(dates[0])+'*.h5')
    file_name = h5_files[0][-17:-6]
    range_files   = list(range(start,end,1))
    for i in range(len(range_files)):
        file_list.append((str(location[0])+str(dates[0])+str(file_name)+str("%03d"%range_files[i])+str(ending[0])))
elif files_mode == 2:
    print("You are selecting from a created list")
    file_list = files
else:
    print('No no no not implemented. If you need it implement it!')



###---- Initilize the class to extract the data----
ev_data = ff.get_data(file_list)
ev_data.lif_gate_sig = lif_gate_sig
ev_data.lif_gate_bac = lif_gate_bac
ev_data.laser = laser


###----Extract data----
###----Extract data----
###----Extract data----
#### 0 = classic, 1 = variable gates, 2 = photon counting, 666 = test mode aka directly plot the data
if da_method == 0:               
    y_arr,signal,back,laser_arr_sig,laser_arr_bac,scan_var,scan1_t_att,scan2_t_att,check_t1,check_t2 = ev_data.classic()
if da_method == 1:
    y_arr,signal,back,laser_arr_sig,laser_arr_bac,scan_var,scan1_t_att,scan2_t_att,check_t1,check_t2 = ev_data.int_lif()
if da_method == 2:
    y_arr,signal,back,laser_arr_sig,laser_arr_bac,scan_var,scan1_t_att,scan2_t_att,check_t1,check_t2 = ev_data.pho_cnt()
if da_method == 666:
    y_arr,signal,back = ev_data.classic_only_data()
    if len(file_list) == 1:
        print('1D File')
        plt.figure(fgc)
        plt.plot(y_arr)
        ff.bug()
    else:                            ####Try if 2D file, else just plot the 1D
        print('2D File')
        plt.figure(fgc)
        gs = plt.GridSpec(1,3)
        fig = plt.figure(fgc)
        ax1 = fig.add_subplot(gs[0,0])
        plt.title('Difference')
        plt.imshow(y_arr,aspect='auto')
        plt.xlabel('laser valve delay')
        plt.colorbar()

        ax1 = fig.add_subplot(gs[0,1])
        plt.title('Signal')
        plt.imshow(signal,aspect='auto')
        plt.xlabel('test')
        plt.colorbar()
        
        ax1 = fig.add_subplot(gs[0,2])
        plt.title('Background')
        plt.imshow(back,aspect='auto')
        plt.colorbar()
        ff.bug()    ### bug is just a function used to stop the programm and show the plots




##########-------You did a fake scan? Scan, for instance,x times the same delay with 10 averages before moving to the next delay
if fake_scan:
    x_fake = np.linspace(start,stop,points)
    y = np.zeros(len(x_fake))
    std = np.zeros(len(x_fake))

    for i in range(points):
        y[i] = np.mean(y_arr[i*len(y_arr)//points:(i+1)*len(y_arr)//points])
        std[i] = np.std(y_arr[i*len(y_arr)//points:(i+1)*len(y_arr)//points])

    fgc+=1
    plt.figure(fgc)
    plt.errorbar(x_fake,y,yerr=std)


####--- Sort the laser shots ---####
if laser:
    ev_laser = ff.sort_laser(laser_arr_sig,laser_arr_bac)
    laser_t  = ev_laser.sort()
else:
    laser_t  = 0


###--- Gating---####
if gating:
    two_d_gate = ff.gating_two_d(y_arr,signal,back,laser_arr_sig,laser_arr_bac,laser_t,sigma_g)
    #signal[signal>6] = np.nan
    #back[back>4] = np.nan

    if gate_data:
        y_arr,signal,back = two_d_gate.gate_data()

    if gate_laser:
        gate_laser_sig,gate_laser_bac,gate_laser_t = two_d_gate.gate_laser()
        ###--
        laser_arr_sig[np.logical_not(gate_laser_sig)] = np.nan
        laser_arr_bac[np.logical_not(gate_laser_bac)] = np.nan

        laser_t_gated = laser_t*1
        laser_t_gated[np.logical_not(gate_laser_t).ravel()] = np.nan

        y_arr[np.logical_not(gate_laser_bac*gate_laser_sig)] = np.nan



if binning:
    print('Binning mode')
    #y_arr,signal,back,laser_arr_sig,laser_arr_bac,scan_var,scan1_t_att,scan2_t_att,check_t1,check_t2
    if len(scan_var)//bin_val*bin_val != len(scan_var):
        cut = (len(scan_var)-len(scan_var)//bin_val*bin_val)
        print('The last '+str(cut)+' value(s) will be removed')
        print(type(scan_var))

        if two_d_plots == True:            ########Depending on if we have a single or mulitple files...
            print("Yummy")
            print('cow')
            y_arr = y_arr[:,:-cut]
            signal = signal[:,:-cut]

            back = back[:,:-cut]
            laser_arr_sig = laser_arr_sig[:,:-cut]
            laser_arr_bac = laser_arr_bac[:,:-cut]
            scan_var = scan_var[:-cut]

            x_mesh,y_mesh = np.meshgrid(scan_var,np.linspace(0,len(scan1_t_att)-1,len(scan1_t_att)))
            y_arr,bla,ueberfluessig = np.histogram2d(x_mesh.ravel(),y_mesh.ravel(),weights=y_arr.ravel(),bins=(len(scan_var)//bin_val,len(scan1_t_att)))
            back,bla,ueberfluessig = np.histogram2d(x_mesh.ravel(),y_mesh.ravel(),weights=back.ravel(),bins=(len(scan_var)//bin_val,len(scan1_t_att)))
            signal,bla,ueberfluessig = np.histogram2d(x_mesh.ravel(),y_mesh.ravel(),weights=signal.ravel(),bins=(len(scan_var)//bin_val,len(scan1_t_att)))
            laser_arr_sig,bla,ueberfluessig = np.histogram2d(x_mesh.ravel(),y_mesh.ravel(),weights=laser_arr_sig.ravel(),bins=(len(scan_var)//bin_val,len(scan1_t_att)))
            laser_arr_bac,bla,ueberfluessig = np.histogram2d(x_mesh.ravel(),y_mesh.ravel(),weights=laser_arr_bac.ravel(),bins=(len(scan_var)//bin_val,len(scan1_t_att)))
            scan_var = np.mean(np.array(scan_var).reshape(-1,bin_val), axis=1)

            y_arr = y_arr.T
            back  = back.T
            signal = signal.T
            laser_arr_sig = laser_arr_sig.T
            laser_arr_bac = laser_arr_bac.T

            print(shape(y_arr))
            print(shape(signal))


        elif two_d_plots == False:         ########Depending on if we have a single or mulitple files...
            print('Yeamensxn')
            y_arr = y_arr[:-cut]
            signal = signal[:-cut]
            back = back[:-cut]
            laser_arr_sig = laser_arr_sig[:-cut]
            laser_arr_bac = laser_arr_bac[:-cut]
            scan_var = scan_var[:-cut]


            y_arr    = np.mean(y_arr.reshape(-1,bin_val), axis=1)
            signal    = np.mean(signal.reshape(-1,bin_val), axis=1)
            back    = np.mean(back.reshape(-1,bin_val), axis=1)
            laser_arr_sig    = np.mean(laser_arr_sig.reshape(-1,bin_val), axis=1)
            laser_arr_bac    = np.mean(laser_arr_bac.reshape(-1,bin_val), axis=1)
            scan_var = np.mean(np.array(scan_var).reshape(-1,bin_val), axis=1)

        else:
            print("Yummy")
            print('cowwy')

        scan_var = scan_var.tolist()



if smooth_2d:
    y_arr_smoothed = ff.smooth_2d_row(y_arr,gating,smooth_val,norm)

scan_variables_numpy =np.array(scan_var)

if plot_pietro == True:
    print('\nPietro loves these plots\n')
    fgc+=1
    plt.figure(fgc)
     #get a proper data vector, not a stupid list
    
    plt.subplot(1, 2, 1)
    plt.plot(scan_variables_numpy*1000000, y_arr*1000)
#    for i in range(0, y_arr.shape[0]):
 #       plt.plot(scan_variables_numpy*1000000, y_arr[i, :]*1000, label=str(files[0][-3:]))

    plt.xlabel('Laser - discharge delay (us)')
    plt.ylabel('Intensity*1000 (arb unit / volts)')
#    labels = ['2 bar', '2.5 bar', '3 bar', '3.5 bar', '1.5 bar']
    #labels = ['V=5.70, T=5.60', 'V8.00, T4.00', 'V8.40, T4.00', 'V6.80, T4.00', 'V6.80, T4.00, 2 real bar', 'V6.80, T4.00, real 3.5 bar, \n aka sanity check', '7', '8', '9', '10', '11', '12']
    labels = ['3.5', '3.0', '2.5', '2.0', '1.5', '1.0', '0.5']
    plt.legend(labels)
    
    # repeat the smoothed plot
    plt.subplot(1, 2, 2)
    data = ff.smooth_1d_fct(y_arr, smooth_val)
    plt.plot(scan_variables_numpy*1000000, data*1000)
#    plt.plot(scan_variables_numpy*1000, np.abs( data - np.nanmax( data ) / 2.) )
#    for i in range(0, y_arr.shape[0]):
#        smooth_val=1.8
#        data = ff.smooth_1d_fct(y_arr,smooth_val)
#        plt.plot(scan_variables_numpy*1000000, data[i, :]*1000, label=str(files[0][-3:]), linewidth=3, alpha=0.8)
 #       plt.plot(scan_variables_numpy*1000, np.abs( data[i, :] - np.nanmax( data[i, :] ) / 2.) )

    plt.xlabel('Laser - discharge delay (us)')
    plt.ylabel('Intensity*1000 (arb unit / volts)')
    plt.legend(labels)
    
    
    
###---- Plot a single file
if plot_format == 0 and files_mode == 0 and plot_pietro == False:
    fgc+=1
    plt.figure(fgc)
    plt.plot(scan_var,y_arr,label=str(files[0][-3:]))
    plt.xlabel("ee")
    #plt.legend()
    #plt.hist(scan_var[0:161],weights=y_arr[0:161],bins=80,histtype='step',color='black')
    plt.grid()
    if smooth_1d == True:
        print('smoothing 1D')
        data = ff.smooth_1d_fct(y_arr,smooth_val)
        np.amax(data)
        plt.plot(scan_var, data)
        plt.plot(scan_var, np.abs(data-np.nanmax(data)/2.))
        

    #fgc= fgc+1
    #plt.figure(fgc)
    #plt.plot(scan_var,y_arr,label=str(files[0][-3:]))
    #plt.legend()
    #plt.hist(y_arr,bins=20,histtype='step',color='black')
    #plt.grid()
    



####--- Plot all data as a 1D Plot
if plot_ente:

    data_ente = ff.analyse(y_arr)
    data_ente.ente =  entee

    nte_arr,nte_arr_err = data_ente.get_ente()
    #nte_arr = nte_arr - np.nanmean(back)

#    if two_d_plots == False:
#        fgc+=1
#        plt.figure(fgc)
#        plt.xlabel('parameter set')
#        plt.ylabel('photons')
#        plt.errorbar(np.linspace(1,entee,len(nte_arr)),nte_arr,yerr=nte_arr_err)
#
#        fgc+=1
#        fig, ax = plt.subplots()
#        plt.plot(signal,label='sig')
#        plt.plot(back,label='bac')
#        xticks = [str(i+1) for i in range(entee)] * int(np.ceil((len(signal)/entee)))
#        ax.set_xticks(range(len(signal)))
#        ax.set_xticklabels(xticks)
#        plt.xlabel('parameter set')
#        plt.ylabel('photons')
#        plt.legend()
#        plt.grid()

    if smooth_2d:
        nte_arr = ff.smooth_2d_row(nte_arr,0,5,0)


    fgc+=1
    plt.figure(fgc)
    gs = plt.GridSpec(1,2)
    fig = plt.figure(fgc)
    ax1 = fig.add_subplot(gs[0,0])
    plt.imshow(nte_arr,interpolation='none',origin='lower',cmap='hot',extent=(np.min(scan_var),np.max(scan_var),np.min(scan2_t_att),np.max(scan2_t_att)),aspect='auto')
    plt.colorbar()

    offset = 0.0005
    ax1 = fig.add_subplot(gs[0,1])
    for i in range(len(nte_arr)):
        plt.errorbar(scan_var,nte_arr[i,:]+offset*i,yerr=nte_arr_err[i,:],label=str(i),fmt='o',capsize=5,color=colors[i])
        plt.plot(scan_var,np.ones(len(scan_var))*offset*i,colors[i])
    plt.grid()
    plt.legend()


if plot_format == 1:
    data = np.nanmean(y_arr,axis=0)
    data_error = np.nanstd(y_arr,axis=0)/np.sqrt(len(y_arr))

    if smooth_1d:
        data = ff.smooth_1d_fct(data,smooth_val)
        data_error = ff.smooth_1d_fct(data_error,smooth_val)

    if norm:
        data = data/np.max(data)
        data_error = data_error/np.max(data)


    fgc+=1
    gs = plt.GridSpec(1,3)
    fig = plt.figure(fgc,figsize=(17,5))
    ax1 = fig.add_subplot(gs[0,0])
    plt.title('Signal')
    plt.imshow(signal,interpolation='none',origin='lower',aspect='auto',extent=[np.min(scan_var),np.max(scan_var),start,end])
    plt.colorbar()
    ax1 = fig.add_subplot(gs[0,1])
    plt.title('Back')
    plt.imshow(back,interpolation='none',origin='lower',aspect='auto',extent=[np.min(scan_var),np.max(scan_var),start,end])
    plt.colorbar()

    n_s = np.sqrt(signal.shape[0] - np.isnan(signal).sum(axis = 0))
    std_s =  np.nanstd(signal, axis = 0)
    ste_s = std_s / n_s

    n_b = np.sqrt(back.shape[0] - np.isnan(back).sum(axis = 0))
    std_b = np.nanstd(back, axis = 0)
    ste_b = std_b / n_b
    error    = np.sqrt(ste_s**2+ste_b**2)

    sum_back = np.nanmean(back,axis=0)
    sum_sig  = np.nanmean(signal,axis=0)
    ax1 = fig.add_subplot(gs[0,2])

    plt.errorbar(scan_var,sum_sig-sum_back,yerr=error,fmt='o',capsize=5)
    plt.grid()
    print('sum of errors:', np.sum(error))
    print('sum of std sig:', np.sum(std_s))
    ###Added fit
    def expfit(x,b,c):
        return b*np.exp(-x/c)

    from lmfit import Model

    gmodel = Model(expfit)
    result = gmodel.fit(sum_sig-sum_back, x=scan_var,b = 2,c = 0.5)#, weights = error)

    b = result.best_values['b']
    c = result.best_values['c']
    print("1/e lifetime is "+str(c))

#    try:
#        popt0,pcov0=curve_fit(expfit,scan_var,sum_sig-sum_back, sigma = error)
#        print(popt0)
#
#        plt.plot(np.linspace(0,2,100),expfit(np.linspace(0,2,100),b,c))
#        print(result.fit_report())
#        plt.grid()
#    except:
#        print("No fit possible")


    fgc+=1
    gs = plt.GridSpec(1,2)
    fig = plt.figure(fgc)
    ax1 = fig.add_subplot(gs[0,0])
    plt.imshow(y_arr,interpolation='none',origin='lower',aspect='auto',extent=[np.min(scan_var),np.max(scan_var),start,end])
    plt.colorbar()
    ax1 = fig.add_subplot(gs[0,1])
    plt.title('Mean of all')
    

    if smooth_1d:
        data = ff.smooth_1d_fct(data,smooth_val)
        data_error = ff.smooth_1d_fct(data_error,smooth_val)

    plt.errorbar(scan_var,data,yerr=data_error)
    plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
    plt.grid()


####--- Get scan paraemters from a 2D scan
if plot_format == 0 and files_mode == 1:
    threshold = np.max(y_arr)*0.2
    gate = y_arr>threshold
    gate_1d = (np.sum(gate,axis=1)>0)

    fgc+=1
    gs = plt.GridSpec(1,2)
    fig = plt.figure(fgc)
    ax1 = fig.add_subplot(gs[0,0])
    plt.imshow(y_arr,interpolation='none',origin='lower',aspect='auto',extent=[np.min(scan_var),np.max(scan_var),start,end])
    plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
    ax1 = fig.add_subplot(gs[0,1])
    plt.imshow(gate,interpolation='none',origin='lower',aspect='auto',extent=[np.min(scan_var),np.max(scan_var),start,end])
    plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))

    scan1_t_att_2d = scan1_t_att[gate_1d]
    scan2_t_att_2d = scan2_t_att[gate_1d]

    for i in range(len(scan1_t_att_2d)):    
        print(str(i+2-start)+str(scan1_t_att_2d[i])+str(scan2_t_att_2d[i]))
    print('\n')



####--- Plot data  ---
if two_d_plots and files_mode < 3:
    y_arr[ np.isnan(y_arr) ] = 0. #replace nan with 0s. Very very dangerous. Where are nans coming from??
    
    fgc+=1
    plt.figure(fgc)
    plt.imshow(y_arr,interpolation='none',origin='lower',aspect='auto',cmap='gnuplot2',extent=[np.min(scan_var),np.max(scan_var),start,end],vmin=0)
    plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
    plt.colorbar()
    
    fgc+=1
    plt.figure(fgc)
    plt.subplot(2,2,1)
    #check for nan
    tt = np.mean(y_arr, axis=0)
    np.isnan(np.sum(tt))
    tt[np.isnan(tt)] = 0. #very dangerous
    
    plt.imshow(y_arr,interpolation='none',origin='lower',aspect='auto',cmap='gnuplot2',extent=[np.min(scan_var),np.max(scan_var),start,end],vmin=0)
    plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
    plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
    plt.subplot(2,3,3)
    plt.plot(np.mean(y_arr, axis=0))
    plt.subplot(3,2,5)
    plt.plot(np.mean(y_arr, axis=1))
    
    #slice dataset in three sets
    if False:
        tt_1 = np.mean(y_arr[0:12, :], axis=0)
        tt_1[np.isnan(tt_1)] = 0.
        tt_2 = np.mean(y_arr[13:24, :], axis=0)
        tt_2[np.isnan(tt_2)] = 0.
        tt_3 = np.mean(y_arr[25:end, :], axis=0)
        tt_3[np.isnan(tt_3)] = 0.
        fgc+=1
        plt.figure(fgc)
        plt.plot(scan_variables_numpy*1e6, ff.smooth_1d_fct(tt_1, smooth_val), linewidth=5)
        plt.plot(scan_variables_numpy*1e6, ff.smooth_1d_fct(tt_2, smooth_val), linewidth=5)
        plt.plot(scan_variables_numpy*1e6, ff.smooth_1d_fct(tt_3, smooth_val), linewidth=5)
        plt.plot(scan_variables_numpy*1e6, ff.smooth_1d_fct(tt, smooth_val), 'k', linewidth=5)
        labels = ['inc times from 350 to 400 us', 'inc times from 400 to 450 us', 'inc times from 450 to 500 us', 'sum']
        plt.legend(labels)


    fgc+=1
    plt.figure(fgc)
    plt.plot(scan_variables_numpy*1000000, tt, '--k', linewidth=2,)
    plt.plot(scan_variables_numpy*1e6, ff.smooth_1d_fct(tt, smooth_val), '-r', linewidth=5)
    plt.xlabel('laser delay (us)')
    
    fgc+=1
    plt.figure(fgc)
    for i in range(0,len(y_arr)):
        
        #plt.plot(scan_var, y_arr[i,:], label=str(i))
        plt.subplot(1,2,1)
        plt.plot(scan_variables_numpy*1000, y_arr[i, :], label=str(i), linewidth=2, alpha=0.9)
        plt.subplot(1,2,2)
        plt.plot(scan_variables_numpy*1000, ff.smooth_1d_fct(y_arr[i, :],smooth_val), label=str(i), linewidth=5, alpha=0.6)
        
    plt.subplot(1,2,1)
    plt.grid()
    plt.subplot(1,2,2)
    plt.grid()
    plt.legend()        
    plt.grid()

    if len(shots_from_2d)>0:
        fgc+=1
        print('shots from 2d:\t\t', shots_from_2d)
        plt.figure(fgc)
        offset = 0.002
        for i in range(len(shots_from_2d)):
            print((int(shots_from_2d[i])))
            plt.plot(scan_var,y_arr[(int(shots_from_2d[i])),:]+i*offset,label=(str(int(shots_from_2d[i]))+' '+str(scan1_t_att[(int(shots_from_2d[i]))])+str(scan2_t_att[(int(shots_from_2d[i]))])),color=colors[i])
            plt.plot(scan_var,np.zeros((len(scan_var)))+i*offset,color=colors[i])
            gate_int = (np.array(scan_var)>2.075e-3)*(np.array(scan_var)<2.17e-3)
            print(np.sum(y_arr[(int(shots_from_2d[i])),gate_int]))
        plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
        plt.legend()



    scan_var_2d = np.array(scan_var*len(y_arr)).reshape(len(y_arr),len(scan_var))

    gates_scan_var = np.array([sumgate[0],sumgate[1]])*1e-6

    super_gate = (scan_var_2d>=gates_scan_var[0])*(scan_var_2d<=gates_scan_var[1])

    y_arr_sun = np.nansum(y_arr*super_gate,axis=1)


    fgc+=1
    plt.figure(fgc)
    plt.hist2d(scan1_t_att.ravel()/1e-6,scan2_t_att.ravel()/1e-6,weights=y_arr_sun,bins=(len(np.unique(scan1_t_att)),len(np.unique(scan2_t_att))),cmap='hot')#,norm=LogNorm())#,vmin=0,vmax=0.08)
    plt.xlabel('Discharge Delay')
    plt.ylabel('Discharge Duration')
    plt.colorbar()


try:
    print('Signal gain:  '+str(np.mean(y_arr[0,:])/(np.mean(y_arr[0,:]))))
    print('Background gain:  '+str(np.mean(back[0,:])/(np.mean(back[0,:]))))
except:
    print('')



####--- Plot Laser  ---
if laser:
    fgc+=1
    gs = plt.GridSpec(1,3)
    fig = plt.figure(fgc,figsize=(17,5))
    ax1 = fig.add_subplot(gs[0,0])
    plt.title('Laser Signal')
    plt.imshow(laser_arr_sig,origin='lower',extent=[np.min(scan_var),np.max(scan_var),start,end],aspect='auto')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
    ax1 = fig.add_subplot(gs[0,1])
    plt.title('Laser Background')
    plt.imshow(laser_arr_bac,origin='lower',extent=[np.min(scan_var),np.max(scan_var),start,end],aspect='auto')
    plt.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))

    ax1 = fig.add_subplot(gs[0,2])
    plt.title('Laser stability')
    plt.scatter(np.linspace(0,len(laser_t)-1,len(laser_t)),laser_t,alpha=0.04,label='All')
    plt.scatter(np.linspace(0,len(laser_t)-1,len(laser_t)),laser_t_gated,alpha=0.04,color='black',label='Gated')
    plt.legend()


if gif:
    import imageio

    block = 20

    data_ani = np.zeros((len(y_arr/block),len(y_arr[0,])))
    err_ani = np.zeros((len(y_arr/block),len(y_arr[0,])))

    for i in range(len(y_arr)//block):
        n_s = np.sqrt(signal[i*block:(i+1)*block].shape[0] - np.isnan(signal[i*block:(i+1)*block]).sum(axis = 0))
        ste_s =  np.nanstd(signal[i*block:(i+1)*block], axis = 0) / n_s

        n_b = np.sqrt(back[i*block:(i+1)*block].shape[0] - np.isnan(back[i*block:(i+1)*block]).sum(axis = 0))
        ste_b = np.nanstd(back[i*block:(i+1)*block], axis = 0) / n_b
        error    = np.sqrt(ste_s**2+ste_b**2)

        sum_back = np.nanmean(back[i*block:(i+1)*block],axis=0)
        sum_sig  = np.nanmean(signal[i*block:(i+1)*block],axis=0)

        data_ani[i,:] = sum_sig-sum_back
        err_ani[i,:] = error

    def plot_for_offset(i):
        fig, ax = plt.subplots(figsize=(10,5))
        ax.errorbar(scan_var,data_ani[i,:],yerr=err_ani[i,:],fmt='o',capsize=5,label=str(i))
        plt.legend()
        plt.ylim(-1,1)
        plt.grid()
        fig.canvas.draw()       # draw the canvas, cache the renderer
        image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        #plt.close(fig)
        return image


    kwargs_write = {'fps':2.0, 'quantizer':'nq'}
    imageio.mimsave('./powers.gif', [plot_for_offset(i) for i in range(len(y_arr)//block)], fps=2)




t2 = time.time()
print(str(t2-t1)+' s')
ff.bug()
