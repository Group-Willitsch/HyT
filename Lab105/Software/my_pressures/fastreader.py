# -*- coding: utf-8 -*-
"""
Created on Mon May 27 16:01:27 2019

@author: Gruppe Willitsch
"""

from PfeifferVacuum import MaxiGauge, MaxiGaugeError
import serial
import datetime
import numpy as np
import pandas as pd
import time
from matplotlib import pyplot as plt

plt.close("all")


class FastReader():
    def __init__(self):
        self.com = 'COM19'
        self.channel = 3
        self.channel_dex = 2
        self.wait = 0
        self.file = './fastread.txt'
        self.points = 50000

        
    def connectPfeiffer(self):
        try:
            self.mg = MaxiGauge(self.com)
            self.connected = True
            print('Connected to Pfeiffer')
        except:
            msg = 'No connection to port %s'%self.com
            print(msg)
            return msg
        
    def pressure(self, channel):
        msg = self.mg.my_pressure(channel)
        pressure = msg[0]
    #        status = msg[1]
        return pressure
    
    def start_reading(self):
        time0 =  datetime.datetime.now()

        self.pressures_dex = []
        self.pressures = []
        self.times = []
        self.times_dex = []
        for i in range(self.points):
            now = (datetime.datetime.now() - time0).total_seconds()
            pressure = self.pressure(self.channel)
            self.times.append(now)
            self.pressures.append(pressure)

            
            now = (datetime.datetime.now() - time0).total_seconds()
            pressure_dex = self.pressure(self.channel_dex)
            self.times_dex.append(now)

            self.pressures_dex.append(pressure_dex)
            time.sleep(self.wait)
            if i%100 == 0:
                print('point %i of %i'%(i, self.points))
        return
    
    
    def disconnect(self):
        self.mg.disconnect()
        return
    

def linear(x,a,b):
    return a * x + b

def optimize(datax, datay, rep, drep, nrep):
    
    costs = []
    ymaxes = []
    xmaxes = []
    
    repreal = np.linspace(rep - drep, rep + drep, nrep)
    
    for i in range(len(repreal)):
        times = []
        for time in datax:
            no = time - np.floor(time / repreal[i]) * repreal[i]
            times.append(no)
            
        datay = datay + np.abs(np.min(datay))
        crit = datay > 0.90 * np.max(datay)
        
#        print(datay)
        times = np.array(times)
        xmax = times[crit]
        ymax = datay[crit]
        
#        print(xmax)
#        print(crit)
#        print(ymax)
#        print(times)
#        print(repreal)
        cost = (np.max(xmax) - np.min(xmax))/repreal[i]
        
        costs.append(cost)
#        ymaxes.append(ymax)
#        xmaxes.append(xmax)
        
    return costs, xmaxes, ymaxes,repreal
    
    

if __name__ == '__main__':
    
    repstr = '_0p1Hz'
    
    directory = 'C:/Users/Gruppe Willitsch/Desktop/data/2019-06-07/'

###################################################################
#    reader = FastReader()
#    reader.connectPfeiffer()
#    reader.start_reading()
#    reader.mg.disconnect()
##    
#    datay = reader.pressures
#    datay_pd = pd.DataFrame(datay)
#    
#    datax = reader.times
#    datax_pd = pd.DataFrame(datax)
#    
#    datay_dex = reader.pressures_dex
#    datay_dex_pd = pd.DataFrame(datay_dex)
#    
#    datax_dex = reader.times_dex
#    datax_dex_pd = pd.DataFrame(datax_dex)
##
#    datay_pd.to_csv('pressures' +repstr +'.txt')
#    datay_dex_pd.to_csv('pressure_dex' +repstr +'.txt')
#    datax_pd.to_csv('times' +repstr +'.txt')
#    datax_dex_pd.to_csv('times_dex' +repstr +'.txt')

###################################################################

    points = 1000
    bino = 1000
    repstart = 5
    drep = 0.1
    
    cut = 1000
    
    datay_pd = pd.read_csv(directory + 'pressures' +repstr +'.txt')['0'][cut:]
    datay_dex_pd = pd.read_csv(directory + 'pressure_dex' +repstr +'.txt')['0'][cut:]
    datax_pd = pd.read_csv(directory + 'times' +repstr +'.txt')['0'][cut:]
    datax_dex_pd = pd.read_csv(directory + 'times_dex' +repstr +'.txt')['0'][cut:]
    
    pressures = []
    timess = []
    
    pressures_tot = []
    times_tot = []
    
    pressures_dex = []
    timess_dex = []
    
    pressures_tot_dex = []
    times_tot_dex = []
    
    datax = np.array(datax_pd).reshape(datax_pd.shape[0])
    datay = np.array(datay_pd).reshape(datax_pd.shape[0])
    datay_dex = np.array(datay_dex_pd).reshape(datay_dex_pd.shape[0])
    datax_dex = np.array(datax_dex_pd).reshape(datax_dex_pd.shape[0])

    
    from scipy.optimize import curve_fit

    popt,pcov=curve_fit(linear,datax,datay)
    popt_dex, pcov_dex = curve_fit(linear,datax_dex,datay_dex)
    
    corr = linear(datax, *popt)
    datay_corr = datay-corr
    
    corr_dex = linear(datax_dex, *popt_dex)
    datay_dex_corr = datay_dex-corr
    
#    y_fft = np.fft.fft(datay_corr)
#    y_fft_power = np.abs(y_fft)
    
#    
#
#    costs, xmaxes, ymaxes,repreal = optimize(datax, datay_corr, repstart, drep, points)
#    plt.plot(repreal,costs)
#    
#    repreal_opt = repreal[np.where(costs == np.min(costs))]
##    repreal_opt = 10
#    for times in datax:
#        no = times - np.floor(times / repreal_opt) * repreal_opt
#        timess.append(no)
#    
#    timess = np.array(timess).reshape(len(timess))
##    plt.hist(timess, weights = datay_corr, bins = bino, histtype='step')
#    bla,bla1 = np.histogram(timess, weights = datay_corr, bins = bino)
#    cor,cor1 = np.histogram(timess, weights = np.ones(len(datay_corr)), bins = bino)
#    
#    dx = cor1[1] - cor1[0]
#    xnew = np.linspace(np.min(cor1) - dx, np.max(cor1) + dx, bino)
#    yvals = bla/cor
#    yvals = yvals - np.min(yvals)
#    yvals = yvals / np.max(yvals)
#    
#    plt.plot(xnew, yvals)
#    plt.xlabel('time /s')
#    plt.ylabel('pressure /mbar')
#    
#    
#    print(' now optimizing dex part')
##    
#    repstart_dex = 5
#    drep_dex = 0.1
#    costs_dex, xmaxes_dex, ymaxes_dex, repreal_dex = optimize(datax_dex, datay_dex_corr, repstart_dex, drep_dex, points)
#    plt.plot(repreal_dex,costs_dex)
#    
#    repreal_opt_dex = repreal_dex[np.where(costs_dex == np.min(costs_dex))]
##    repreal_opt = 0.09
#    for times in datax_dex:
#        no = times - np.floor(times / repreal_opt_dex) * repreal_opt_dex
#        timess_dex.append(no)
#    
#    timess_dex = np.array(timess_dex).reshape(len(timess_dex))
##    plt.hist(timess_dex, weights = datay_dex_corr, bins = bino, histtype='step')
#    plt.xlabel('time /s')
#    plt.ylabel('pressure /mbar')
#    
#    bla,bla1 = np.histogram(timess_dex, weights = datay_dex_corr, bins = bino)
#    cor,cor1 = np.histogram(timess_dex, weights = np.ones(len(datay_dex_corr)), bins = bino)
#    
#    yvals_dex = bla/cor
#    yvals_dex = yvals_dex - np.min(yvals_dex)
#    yvals_dex = yvals_dex / np.max(yvals_dex)
#    
#    dx = cor1[1] - cor1[0]
#    xnew_dex = np.linspace(np.min(cor1) - dx, np.max(cor1) + dx, bino)
#    plt.plot(xnew_dex, yvals_dex)
#    
#    xoffset = np.abs(1.72924 - 1.77869)
#    xnew_dex_offs = xnew_dex - xoffset
#    
#    cutting = (xnew_dex_offs < 4.5) * (xnew_dex_offs > 1.0)
#    
#    plt.plot(xnew_dex_offs[cutting], yvals_dex[cutting], label = 'decelerator')
#    plt.plot(xnew[cutting], yvals[cutting], label = 'trap')
#    
#    plt.plot(xnew_dex_offs[cutting], yvals[cutting] - yvals_dex[cutting], label = 'diff')
#    
#    plt.legend()
#    
#    plt.xlabel('time /s')
#    plt.ylabel('pressure /arb. u.')
#    
    
    
    plt.plot()