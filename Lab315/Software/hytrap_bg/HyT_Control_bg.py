#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QThread, pyqtSignal
import sys
import readfile as rf
import gui
import pulseblaster as pb
#~ import lcwr_scope as scope
import plotsequence as plotseq
import pulseseq as ps
import numpy as np
import time
from matplotlib import pyplot as plt
import LeCroy_visa as lc
import pyqtgraph as pg
import datetime
import os
import h5py
import pandas as pd
import gc
import pulsare_comm as pulsare
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d # smoothing
from scipy.signal._peak_finding import find_peaks  # peak finding algorithm



units={'ns':1e-9, 'us':1e-6, 'ms':1e-3, 's':1}

'''

    13.05.2019
    
    restructuring such that class variables are inherited by subthreads
    starting with the new FakeScan

'''

#start by reading the config file
version = r'0.1 α' 
devmode = 0     #if 1: developer mode which means no communication to other devices
filename = './input/HyT-input.dat'    #config file
params, values=rf.readpars(filename)
for i in range (len(params)):
    try:
        exec(params[i]+"=%s"%(values[i]))
    except (NameError, SyntaxError):
        if values[i] == '':
            print('shid')
            continue
        print(params[i], values[i])
        exec(params[i]+"='%s'"%(values[i].strip()))



class TrigPlotThread(QThread):
    TrigLog = pyqtSignal(object)
    ScanLog = pyqtSignal(object)

    PlotSeq = pyqtSignal(object,object,object)

    def __init__(self,output,checkboxvals):
        QThread.__init__(self)
        self.output = output
        self.checkboxvals = checkboxvals
        self.plotsequence = plotseq.plotsequence

    def __del__(self):
        self.wait()
        
    def run(self):
        output = self.output
        checkboxvals = self.checkboxvals
        xvals,yvals,newlabels=self.plotsequence(output, checkboxvals)
        self.TrigLog.emit('Done')   
        self.PlotSeq.emit(xvals,yvals,newlabels)   
        #~ self.terminate()          
        return
  
class DelayScanThread(QThread):
 
    def __init__(self, units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, t2jumpfile,
                    scanchannel,relchannel,bgchannel,delay,numberofsegments,scope,channel):
        QThread.__init__(self)

        self.units1 = units1
        self.chno1 = chno1
        self.time1 = time1
        
        self.units2 = units2
        self.chno2 = chno2
        self.time2 = time2
        
        self.freqs = freqs
        self.pulsenature = pulsenature
        self.checkboxvals = checkboxvals
        
        self.posnegvals = posnegvals
        self.t2jumpfile = t2jumpfile
        self.scanchannel = scanchannel
        self.relchannel = relchannel
        self.bgchannel = bgchannel
        self.delay = delay
        self.numberofsegments = numberofsegments
        self.scope = scope
        self.channel = channel

        
    def __del__(self):
        self.wait()


    def tof_scan(self): 
        
        
        try:
            pb.pulseblaster_stop()
        except:
            print('pulseblaster not yet running')
            
        pulseseq = ps.sequence(self.units1, self.chno1, self.time1, self.units2,
                               self.chno2, self.time2, self.freqs, self.pulsenature,
                               self.checkboxvals, self.posnegvals, self.scanchannel,
                               self.relchannel, self.bgchannel, self.delay, self.t2jumpfile)
        
        if 'b' in self.scanchannel:
            index = int(self.scanchannel[:-1])-1
            pulseseq.chno1[index] = self.relchannel
            pulseseq.time1[index] = self.delay
            
        elif 'e' in self.scanchannel:
            index = int(self.scanchannel[:-1])-1
            pulseseq.chno2[index] = self.relchannel
            pulseseq.time2[index] = self.delay
            
        pulseseq.seq()
        output = pulseseq.output
                
        pb.pulseblaster_program(output)
        time.sleep(1/np.min(self.freqs) * 1.1)

        pb.pulseblaster_start()
    
        count = 0
        
        try:
            while count < self.numberofsegments:
                count = self.scope.segment_counter(self.channel)[0]
                time.sleep(5e-2)
        except:
            print('tofscan, numberofsegments')
            
        return

    def run(self):
        self.tof_scan()
        return
    
class ListenerThread(QThread):
    Plot = pyqtSignal(object,object,object)
    Log =  pyqtSignal(object)
    Bar = pyqtSignal(object)
    SetWaveLength = pyqtSignal(object)
    ScanTime = pyqtSignal(object)
    FailSafe = pyqtSignal(object)
    Freerun_Points = pyqtSignal(object)
    
    def __init__(self, parent):
        QThread.__init__(self)
        self.guithread = parent
        self.numberofsegments = self.guithread.numberofsegments
        
        now = datetime.datetime.now()
        self.file_pre = '%04d-%02d-%02d-'%(now.year, now.month, now.day)
        self.file_end = '.h5'
        
        self.directory = 'C:/Users/Gruppe Willitsch/Desktop/data/%04d-%02d-%02d/'%(now.year, now.month, now.day)
        if not os.path.exists(self.directory):
            self.Log.emit('creating directory \n %s'%self.directory)
            os.makedirs(self.directory)
            
            
        self.count = 0
            
        for file in os.listdir(self.directory):
            if file.endswith('.h5'):
                number = int(file[-6:-3])
                if number>=self.count:
                    self.count = number
                    
        self.scanmode = self.guithread.ScanModeBox.currentText()
        self.delaymode = self.guithread.DelayModeBox.currentText()
        self.bgmode = self.guithread.BGModeBox.currentText()
        self.repetitions = self.guithread.ScanRepsBox.value()
        self.numberofsegments = self.guithread.AverageSBox.value()
        
        self.scanchannel = 'No'
        self.relchannel = 'No'
        
        self.bgchannel = int(self.guithread.ScanBGChannelBox.currentText())
        
        
        self.channel = 1  #signal channel scope
        self.pdchannel = 2
        self.scopewait = 2e-2
        
        
        self.x1 = 0.195e-6 #int time start
        self.x2 = 4.0e-6 #integration time end
        self.xb1 = 0.195e-6 #background integration
        self.xb2 = 4.0e-6 #
        
        self.xpd1 = 2.5e-7 # probably integration times of the detector
        self.xpd2 = 3.5e-7
        self.xpdb1 = 3.5e-7
        self.xpdb2 = 6e-7
        
        self.bgchannel = int(self.guithread.ScanBGChannelBox.currentText())
        self.freerun_points = 10
        self.number = 0
        
        self.unit = units[self.guithread.ScanUnitBox.currentText()]


    def __del__(self):
        self.wait()
        
    def filename_update(self):
        self.count += 1
        self.filename = self.file_pre + '%03d'%(self.count) + self.file_end
        return
    
    def file_head3d(self):
        self.f = h5py.File(self.directory + self.filename, 'w')
                            
        grp_s = self.f.create_group("signal")
        grp_bg = self.f.create_group("background")
        self.f.attrs['scanchannel'] = self.scanchannel
        self.f.attrs['relchannel'] = self.relchannel
        self.f.attrs['bgchannel'] = self.bgchannel
        self.f.attrs['averages'] = self.numberofsegments
        
        
#                    
        inputs = ['self.guithread.units1', 'self.guithread.chno1', 'self.guithread.time1', 'self.guithread.units2',
                      'self.guithread.chno2', 'self.guithread.time2', 'self.guithread.freqs', 'self.guithread.pulsenature', 
                      'self.guithread.checkboxvals', 'self.guithread.posnegvals', 'self.guithread.t2jumpfile',
                      'self.scanchannel', 'self.relchannel', 'self.bgchannel',
                      'self.numberofsegments'] 
        
        
        for inputz in inputs:
            lelist = eval(inputz)
            word1 = 'self.'
            word2 = 'guithread.'
            words = [word1, word2]
            
            length = 0
            for word in words:
                if word in inputz:
                    length += len(word)
                    
            name = inputz[length:]
            
            
            try:
                asciiList = [n.encode("ascii", "ignore") for n in lelist]
                self.f.create_dataset('inputs/%s'%name, (len(asciiList),1),'S10', asciiList)#, compression = 'gzip')

            except:
                asciiList = eval(inputz)
                self.f.create_dataset('inputs/%s'%name, data = asciiList)#, compression = 'gzip')

#                        
        
#        print(self.scanchannels,self.relchannels, self.scantimes)
        for i in range(2):
            try:
                self.f.attrs['scan%i'%i] = self.scanchannels[i]
                self.f.attrs['scanrel%i'%i] = self.relchannels[i]

            except:
#                print('no scan%i'%i)
                a = 1
                
        self.f.close()
    
    def Scan3d(self):
        self.Log.emit('Starting grid Scan.')
        
        
        inputfile = self.guithread.Scan3dText.text()
        self.data = pd.read_csv(self.guithread.Scan3dText.text(), delim_whitespace = True, header = None,comment='#')
                           
        self.Log.emit('Scan starts with inputfile %s'%(inputfile))
        
        self.values0 = np.arange(self.data[2][0], self.data[3][0]+0.1*self.data[4][0], self.data[4][0])
        self.values1 = np.arange(self.data[2][1], self.data[3][1]+0.1*self.data[4][1], self.data[4][1])
        
        if self.delaymode == 'list':
            inputfile = self.guithread.Scan3dListText.text()
            self.Log.emit('Using list file for delays: %s'%inputfile)
            self.values2 = np.array(pd.read_csv(inputfile, delim_whitespace = True, header = None,comment='#')[0])
                                           

        elif self.delaymode == 'normal':
            self.values2 = np.arange(self.data[2][2], self.data[3][2] + 0.1 * self.data[4][2], self.data[4][2])
            self.Log.emit('Using start, stop and step values from %s'%inputfile)

        self.vals1,self.vals2,self.vals3 = np.array(np.meshgrid(self.values1,self.values2,self.values0))
                
        self.scanchannels = np.array(self.data[0])
        self.relchannels = np.array(self.data[1])
        
        self.scanchannel = self.scanchannels[-1]
        self.relchannel = self.relchannels[-1]
        
    def Scan1d(self):
        self.Log.emit('Starting 1d Scan. %s'%self.scanmode)
        
        self.scanchannel = self.guithread.ScanChannelBox.currentText()
        self.relchannel = self.guithread.ScanRelChannelBox.currentText()
        
        if self.delaymode == 'list':
            inputfile = self.guithread.Scan3dListText.text()
            self.Log.emit('Using list file for delays: %s'%inputfile)

            self.values2 = np.array(pd.read_csv(inputfile, delim_whitespace = True, header = None,comment='#')[0])
                                           

        elif self.delaymode == 'normal':
            self.values2 = np.arange(self.guithread.ScanStartBox.value() * self.unit , self.guithread.ScanStopBox.value() * self.unit  + 0.1 * self.guithread.ScanStepBox.value() * self.unit , self.guithread.ScanStepBox.value() * self.unit )
            self.Log.emit('Using start, stop and step values from GUI')
            
        if self.scanmode == 'wavelength':
            self.scanchannel = 'No'
            self.relchannel = 'No'
            self.values2 = np.arange(self.guithread.ScanStartBoxLaser.value(), self.guithread.ScanStopBoxLaser.value(), self.guithread.ScanStepBoxLaser.value())
            self.Log.emit('Wavelength Scan')

        self.values0 = [1e-6] # some dummy values, not used in 1d scan!
        self.values1 = [1e-6]
        self.vals1,self.vals2,self.vals3 = np.array(np.meshgrid(self.values1,self.values2,self.values0))
                
        self.scanchannels = np.array(['No','No',self.scanchannel])
        self.relchannels = np.array(['No','No',self.relchannel])
        
    def FreeRun(self):
        self.Log.emit('Start continuous mode')
        datax = []
        datay = []
        
        delay = 0
        
        datax_bg = []
        datay_bg = []
        
        dataypd = []
        dataypd_bg = []
        count = 0
        while True:
            self.Freerun_Points.emit(self)
            
            self.scope.clear_sweeps()
            

            while count < self.numberofsegments:
                count = self.scope.segment_counter(self.channel)[0]
                time.sleep(5e-2)
            count = 0
            self.data = self.scope.get_waveform(self.channel)
            self.photodiode = self.scope.get_waveform(self.pdchannel) #photodiode signal
            
            
            ypdb = np.mean(self.photodiode[1][np.logical_and(self.photodiode[0]>self.xpdb1, self.photodiode[0]<self.xpdb2)])
            self.photodiode[1] = self.photodiode[1]-ypdb
        
            y = np.mean(self.data[1][np.logical_and(self.data[0]>self.x1, self.data[0]<self.x2)])
            yb = np.mean(self.data[1][np.logical_and(self.data[0]>self.xb1, self.data[0]<self.xb2)])
#                    datay.append(y-yb)

            if len(datax) >= self.freerun_points:
                datax.pop(0)
                datay.pop(0)
                dataypd.pop(0)
            datax.append(delay)
            
            datay.append(y) # here might be what u are looking for

            
            #photodiode
            ypd = np.sum(self.photodiode[1][np.logical_and((self.photodiode[0])>self.xpd1, (self.photodiode[0])<self.xpd2)])
            dataypd.append(ypd)
            
            
            #plot
            self.Plot.emit(self.data[0],self.data[1],1)
            self.Plot.emit(datax,datay,2) 
            self.Plot.emit(self.photodiode[0],self.photodiode[1],5)
            delay += 1

        
    def ScanFake(self):
        datax = []
        datay = []
        
        delay = 0
        
        datax_bg = []
        datay_bg = []
        
        dataypd = []
        dataypd_bg = []
        
        try:
            pb.pulseblaster_stop()
        except:
            print('pulseblaster not yet running')
        
        pointreps = self.guithread.ScanFakePointRepBox.value()
        totalreps = self.guithread.ScanFakeTotalRepBox.value()
        
        count = 0
        
        
        self.scanchannel = self.guithread.ScanChannelBox.currentText()
        self.relchannel = self.guithread.ScanRelChannelBox.currentText()
        
        self.scanchannels = np.array(['No','No',self.scanchannel])
        self.relchannels = np.array(['No','No',self.relchannel])
        
        self.scantimes = ['No','No']

        self.filename_update()
        self.file_head3d()
        
        for total in range(totalreps):
            self.Fail = 'None'
            
            for point in range(pointreps):
                self.scope.clear_sweeps()
                
                self.Log.emit('Scanning Fake: point %i'%point + ' in position %i'%total)
                self.scan = DelayScanThread(self.guithread.units1, self.guithread.chno1, self.guithread.time1, self.guithread.units2,
                                            self.guithread.chno2, self.guithread.time2, self.guithread.freqs, self.guithread.pulsenature, 
                                            self.guithread.checkboxvals, self.guithread.posnegvals, self.guithread.t2jumpfile,
                                            'No', 'No', 'No',
                                            delay, self.numberofsegments,
                                            self.scope, self.channel)
                self.scan.start()
                        
                while self.scan.isRunning():
                    time.sleep(self.scopewait)
                  

                datax.append(count)

                self.data = self.scope.get_waveform(self.channel)
                self.photodiode = self.scope.get_waveform(self.pdchannel) #photodiode signal
                
                
                ypdb = np.mean(self.photodiode[1][np.logical_and(self.photodiode[0]>self.xpdb1, self.photodiode[0]<self.xpdb2)])
                self.photodiode[1] = self.photodiode[1]-ypdb
            
                y = np.mean(self.data[1][np.logical_and(self.data[0]>self.x1, self.data[0]<self.x2)])
                yb = np.mean(self.data[1][np.logical_and(self.data[0]>self.xb1, self.data[0]<self.xb2)])
#                    datay.append(y-yb)

                # y will contain at the end the number of photons per laser shot detected.
                # data[0] contains the timebase
                # data[1] contains the signal from the oscilloscope
                # x1, x2 are the gating times for the signal, currently 1.95e-7, 5e-6
                # xb1, xb2 are the gating times from the background, currently 1.95e-7, 5e-6
                photon_counting_in_fake_mode = False
                if photon_counting_in_fake_mode:  # implement photon counting in during fake scans
                    gate_for_signal = np.logical_and(self.data[0] >= self.x1,
                                                     self.data[0] <= self.x2)
                    gate_for_background = np.logical_and(self.data[0] >= self.xb1,
                                                         self.data[0] <= self.xb2)

                    # get average-corrected and time-gated oscilloscope traces
                    signal_trace = self.data[1][gate_for_signal] * self.numberofsegments  # very imp: trace is
                    # multiplied by the number of laser shot, s.t. it works with all number of averages
                    background_trace = self.data[1][gate_for_background] * self.numberofsegments

                    # gaussian smoothing
                    signal_trace = gaussian_filter1d(signal_trace, 2.5)  # smoothing
                    background_trace = gaussian_filter1d(background_trace, 2.5)

                    my_prominence = 80e-3  # in Volts! Typ. a peak is 50 mV to 200 mV high
                    my_width = 10  # in vector index units, typ. FWHM around 10 points
                    [signal_peaks, _] = find_peaks(signal_trace, prominence=my_prominence,
                                                   width=my_width)  # peak finding function call
                    [background_peaks, _] = find_peaks(background_trace, prominence=my_prominence, width=my_width)
                    y = len(signal_peaks)
                    yb = len(background_peaks)
                    print('number of photons y :', y)
                    print('number of photons yb:', yb)

                datay.append(y) # here might be what u are looking for

                
                #photodiode
                ypd = np.sum(self.photodiode[1][np.logical_and((self.photodiode[0])>self.xpd1, (self.photodiode[0])<self.xpd2)])
                dataypd.append(ypd)
                
                #plot
                self.Plot.emit(self.data[0],self.data[1],1)
                self.Plot.emit(datax,datay,2) 
                self.Plot.emit(self.photodiode[0],self.photodiode[1],5)
#                    self.Plot.emit(datax,dataypd,6)
                
                self.f = h5py.File(self.directory+self.filename, 'a')
                dset_s = self.f.create_dataset('signal/'+'%03d-'%total+'%03d'%point, data = np.array(self.data[1]), compression = 'gzip')
                dset_s.attrs['x1'] = self.x1
                dset_s.attrs['x2'] = self.x2
                dset_s.attrs['y'] = y
                dset_s.attrs['xb1'] = self.xb1
                dset_s.attrs['xb2'] = self.xb2
                dset_s.attrs['y_oc'] = y-yb
                
                
                pd_bd = np.histogram(self.photodiode[0],weights=np.array(self.photodiode[1]),bins=200)
                
                # photodiode
                dset_s = self.f.create_dataset('signal_pd/'+'%03d-'%total+'%03d'%point, data = np.array(pd_bd[0]), compression = 'gzip')
                dset_s.attrs['xpd1'] = self.xpd1
                dset_s.attrs['xpd2'] = self.xpd2
                dset_s.attrs['ypd'] = ypd
                dset_s.attrs['xpdb1'] = self.xpdb1
                dset_s.attrs['xpdb2'] = self.xpdb2
                                   
#                                            
#                    begin bg scan
  
                self.scan = DelayScanThread(self.guithread.units1, self.guithread.chno1, self.guithread.time1, self.guithread.units2,
                                            self.guithread.chno2, self.guithread.time2, self.guithread.freqs, self.guithread.pulsenature, 
                                            self.guithread.checkboxvals, self.guithread.posnegvals, self.guithread.t2jumpfile,
                                            'No', 'No', self.bgchannel,
                                            delay, self.numberofsegments,
                                            self.scope, self.channel)
                
                self.scope.clear_sweeps()
                self.scan.start()
                
                                
                while self.scan.isRunning():
                    time.sleep(self.scopewait)
                    
                self.data = self.scope.get_waveform(self.channel)
                self.photodiode = self.scope.get_waveform(self.pdchannel) #photodiode signal
                ypdb = np.mean(self.photodiode[1][np.logical_and(self.photodiode[0]>self.xpdb1, self.photodiode[0]<self.xpdb2)])
                self.photodiode[1] = self.photodiode[1]-ypdb

                
                datax_bg.append(count)
                y = np.mean(self.data[1][np.logical_and(self.data[0]>self.x1, self.data[0]<self.x2)])
                yb = np.mean(self.data[1][np.logical_and(self.data[0]>self.xb1, self.data[0]<self.xb2)])
#                    datay_bg.append(y-yb)
                datay_bg.append(y) # super XXX


                ypd = np.sum(self.photodiode[1][np.logical_and((self.photodiode[0])>self.xpd1, (self.photodiode[0])<self.xpd2)])
                dataypd_bg.append(ypd)
                
                self.Plot.emit(self.data[0],self.data[1],1)
                self.Plot.emit(datax,datay_bg,3)  
                self.Plot.emit(np.array(datax),np.array(datay)-np.array(datay_bg),4)
                
                self.Plot.emit(self.photodiode[0],self.photodiode[1],5)
                self.Plot.emit(datax,dataypd_bg,6)
                
                dset_bg = self.f.create_dataset('background/'+'%03d-'%total+'%03d'%point, data = np.array(self.data[1]), compression = 'gzip')
                dset_bg.attrs['x1'] = self.x1
                dset_bg.attrs['x2'] = self.x2
                dset_bg.attrs['y'] = y
                dset_bg.attrs['xb1'] = self.xb1
                dset_bg.attrs['xb2'] = self.xb2
                dset_bg.attrs['y_oc'] = y-yb

                # photodiode
                pd_bd = np.histogram(self.photodiode[0],weights=np.array(self.photodiode[1]),bins=200)

                dset_s = self.f.create_dataset('background_pd/'+'%03d-'%total+'%03d'%point, data = np.array(pd_bd[0]), compression = 'gzip')
                dset_s.attrs['xpd1'] = self.xpd1
                dset_s.attrs['xpd2'] = self.xpd2
                dset_s.attrs['ypd'] = ypd
                dset_s.attrs['xpdb1'] = self.xpdb1
                dset_s.attrs['xpdb2'] = self.xpdb2

                self.f.close()
                try:   
                    self.f.close()
                except:
                    print('already closed')
                    
                count += 1
                    
            self.FailSafe.emit(self)
            
            while self.Fail == 'None':
                time.sleep(1e-2)
                
            print(self.Fail)
                
            if not self.Fail:
                self.Log.emit('Scanning Fake aborted')
                break
            else:
                self.Log.emit('Next position')
                
            
        self.Log.emit('Scanning Fake finished')
        
    def StartScan(self):
        
        self.number = 0
        self.totallength = self.vals3.T.shape[0] * self.vals3.T.shape[1] * self.repetitions
        letime = 2*np.shape(self.vals2.T[0][0])[0]/np.min(self.guithread.freqs)*self.numberofsegments*(self.totallength-self.number)
        day = letime // (24*3600)
        letime = letime % (24*3600)
        hour = letime // 3600
        letime %= 3600
        minutes = letime // 60
        letime %= 60
        seconds = letime
        letime = "dd:hh:mm:ss-> %02d:%02d:%02d:%02d" % (day, hour, minutes, seconds)
                
        self.ScanTime.emit('Approx. time left: %s'%(letime))
        
        for repetition in range(self.repetitions):
        
            for index0 in range (self.vals3.T.shape[0]):
                for index1 in range (self.vals3.T.shape[1]):
                    scanval1 = self.vals3.T[index0][index1][0]
                    scanval2 = self.vals1.T[index0][index1][0]
                    self.scantimes = [scanval1,scanval2]
                    
                    delays = self.vals2.T[index0][index1]
    
                    if self.guithread.ScanModeBox.currentText() == '3d-scan':
                        for index in range(2):
                            chno1r = int(self.scanchannels[index][:-1])
                            
                            if self.scanchannels[index][-1] == 'b':
                                self.guithread.chno1[chno1r-1] = self.relchannels[index]
                                self.guithread.time1[chno1r-1] = self.scantimes[index]
                                
                            elif self.scanchannels[index][-1] == 'e':
                                self.guithread.chno2[chno1r-1] = self.relchannels[index]
                                self.guithread.time2[chno1r-1] = self.scantimes[index]

                    for file in os.listdir(self.directory):
                        
                        if file.endswith('.h5'):
                            number = int(file[-6:-3])
                            if number>=self.count:
                                self.count = number
    #                        print('next')  
    
                    self.filename_update()
                    self.file_head3d()
                    self.f = h5py.File(self.directory + self.filename, 'a')
                    
                    for i in range(2):
                        self.f.attrs['scan%i_time'%i] = self.scantimes[i]
                    self.f.close()
                    

                    
            
                    self.Log.emit('Scanning %s rel to %s with bg %s'%(self.scanchannel, self.relchannel,self.bgchannel))
                    self.Log.emit('%s'%self.scanchannels)
                    self.Log.emit('%s'%self.relchannels)
                    self.Log.emit('%s'%self.scantimes)
                    datax = []
                    datay = []
                    
                    datax_bg = []
                    datay_bg = []
                    
                    dataypd = []
                    dataypd_bg = []
                    
                    self.number +=1
                    
    #                print(delays)
                    
                    for delay in delays:
                        now = datetime.datetime.now()
                        tstamp = '%04d-%02d-%02d-%02d:%02d%02d'%(now.year, now.month, now.day, now.hour, now.minute, now.second)
                                       
    #                            signal scan
                        time0 = time.clock()
                        self.scope.clear_sweeps()
    #                    scope.arm()
                        
                        if self.scanmode == 'wavelength':
                            self.SetWaveLength.emit(delay)
                            time.sleep(0.5)
                            
                            
                        self.scan = DelayScanThread(self.guithread.units1, self.guithread.chno1, self.guithread.time1, self.guithread.units2,
                                                    self.guithread.chno2, self.guithread.time2, self.guithread.freqs, self.guithread.pulsenature, 
                                                    self.guithread.checkboxvals, self.guithread.posnegvals, self.guithread.t2jumpfile,
                                                    self.scanchannel, self.relchannel, 'None',
                                                    delay, self.numberofsegments,
                                                    self.scope, self.channel)
                        self.scan.start()
                        
                        while self.scan.isRunning():
                            time.sleep(self.scopewait)
                          
                        datax.append(delay)
#                        print(delay)
    
                        self.data = self.scope.get_waveform(self.channel)
                        self.photodiode = self.scope.get_waveform(self.pdchannel) #photodiode signal
                        
                        
                        ypdb = np.mean(self.photodiode[1][np.logical_and(self.photodiode[0]>self.xpdb1, self.photodiode[0]<self.xpdb2)])
                        self.photodiode[1] = self.photodiode[1]-ypdb
                    
                        y = np.mean(self.data[1][np.logical_and(self.data[0]>self.x1, self.data[0]<self.x2)])
                        yb = np.mean(self.data[1][np.logical_and(self.data[0]>self.xb1, self.data[0]<self.xb2)])
    #                    datay.append(y-yb)
                        datay.append(y) # here might be what u are looking for
    
                        
                        #photodiode
                        ypd = np.sum(self.photodiode[1][np.logical_and((self.photodiode[0])>self.xpd1, (self.photodiode[0])<self.xpd2)])
                        dataypd.append(ypd)
                        
                        #plot
                        self.Plot.emit(self.data[0],self.data[1],1)
                        self.Plot.emit(datax,datay,2) 
                        self.Plot.emit(self.photodiode[0],self.photodiode[1],5)
    #                    self.Plot.emit(datax,dataypd,6)
                        
                        self.f = h5py.File(self.directory + self.filename, 'a')
                        dset_s = self.f.create_dataset('signal/'+'%s'%round(delay,8), data = np.array(self.data[1]), compression = 'gzip')
                        dset_s.attrs['x1'] = self.x1
                        dset_s.attrs['x2'] = self.x2
                        dset_s.attrs['y'] = y
                        dset_s.attrs['xb1'] = self.xb1
                        dset_s.attrs['xb2'] = self.xb2
                        dset_s.attrs['y_oc'] = y-yb
                        dset_s.attrs['tstamp'] = tstamp

                        
                        if delay==delays[0]:
                            self.f.attrs['osci_tstart'] = self.data[0][0]
                            self.f.attrs['osci_tstop'] = self.data[0][-1]
                            self.f.attrs['osci_tbins'] = len(self.data[0])
    
    
    
                        
                        pd_bd = np.histogram(self.photodiode[0],weights=np.array(self.photodiode[1]),bins=200)
                        
                        # photodiode
                        dset_s = self.f.create_dataset('signal_pd/'+'%s'%round(delay,8), data = np.array(pd_bd[0]), compression = 'gzip')
                        dset_s.attrs['xpd1'] = self.xpd1
                        dset_s.attrs['xpd2'] = self.xpd2
                        dset_s.attrs['ypd'] = ypd
                        dset_s.attrs['xpdb1'] = self.xpdb1
                        dset_s.attrs['xpdb2'] = self.xpdb2
                                           
    #                                            
    #                    begin bg scan
                        now = datetime.datetime.now()
                        tstamp = '%04d-%02d-%02d-%02d:%02d%02d'%(now.year, now.month, now.day, now.hour, now.minute, now.second)
            
                        self.scan = DelayScanThread(self.guithread.units1, self.guithread.chno1, self.guithread.time1, self.guithread.units2,
                                                    self.guithread.chno2, self.guithread.time2, self.guithread.freqs, self.guithread.pulsenature, 
                                                    self.guithread.checkboxvals, self.guithread.posnegvals, self.guithread.t2jumpfile,
                                                    self.scanchannel, self.relchannel, self.bgchannel,
                                                    delay, self.numberofsegments,
                                                    self.scope, self.channel)
                        
                        self.scan.start()
                        
                        self.scope.clear_sweeps()
    #                    scope.arm()
                        
                        while self.scan.isRunning():
                            time.sleep(self.scopewait)
                            
                        self.data = self.scope.get_waveform(self.channel)
                        self.photodiode = self.scope.get_waveform(self.pdchannel) #photodiode signal
                        ypdb = np.mean(self.photodiode[1][np.logical_and(self.photodiode[0]>self.xpdb1, self.photodiode[0]<self.xpdb2)])
                        self.photodiode[1] = self.photodiode[1]-ypdb
    
                        
                        datax_bg.append(delay)
                        y = np.mean(self.data[1][np.logical_and(self.data[0]>self.x1, self.data[0]<self.x2)])
                        yb = np.mean(self.data[1][np.logical_and(self.data[0]>self.xb1, self.data[0]<self.xb2)])
    #                    datay_bg.append(y-yb)
                        datay_bg.append(y) # super XXX
    
    
                        ypd = np.sum(self.photodiode[1][np.logical_and((self.photodiode[0])>self.xpd1, (self.photodiode[0])<self.xpd2)])
                        dataypd_bg.append(ypd)
                        
                        self.Plot.emit(self.data[0],self.data[1],1)
                        self.Plot.emit(datax_bg,datay_bg,3)  
                        self.Plot.emit(np.array(datax),np.array(datay)-np.array(datay_bg),4)
                        
                        self.Plot.emit(self.photodiode[0],self.photodiode[1],5)
                        self.Plot.emit(datax,dataypd_bg,6)
                        
                        dset_bg = self.f.create_dataset('background/'+'%s'%round(delay,8), data = np.array(self.data[1]), compression = 'gzip')
                        dset_bg.attrs['x1'] = self.x1
                        dset_bg.attrs['x2'] = self.x2
                        dset_bg.attrs['y'] = y
                        dset_bg.attrs['xb1'] = self.xb1
                        dset_bg.attrs['xb2'] = self.xb2
                        dset_bg.attrs['y_oc'] = y-yb
                        dset_s.attrs['tstamp'] = tstamp

    
                        # photodiode
                        pd_bd = np.histogram(self.photodiode[0],weights=np.array(self.photodiode[1]),bins=200)
    
                        dset_s = self.f.create_dataset('background_pd/'+'%s'%round(delay,8), data = np.array(pd_bd[0]), compression = 'gzip')
                        dset_s.attrs['xpd1'] = self.xpd1
                        dset_s.attrs['xpd2'] = self.xpd2
                        dset_s.attrs['ypd'] = ypd
                        dset_s.attrs['xpdb1'] = self.xpdb1
                        dset_s.attrs['xpdb2'] = self.xpdb2
    
                        self.f.close()
                        try:   
                            self.f.close()
                        except:
                            print('already closed')
                            
                    time1 = time.clock()
    #                letime = (time1-time0)*(self.totallength-self.number)
                    letime = 2*np.shape(delays)[0]/np.min(self.guithread.freqs)*self.numberofsegments*(self.totallength-self.number)
                    day = letime // (24*3600)
                    letime = letime % (24*3600)
                    hour = letime // 3600
                    letime %= 3600
                    minutes = letime // 60
                    letime %= 60
                    seconds = letime
                    
                    letime = "dd:hh:mm:ss-> %02d:%02d:%02d:%02d" % (day, hour, minutes, seconds)
                    
                    
                    
                    self.ScanTime.emit('Approx. time left: %s'%(letime))
                    self.Bar.emit(round(self.number / self.totallength * 100, 8))
#
#            
    
    def run(self):
        self.scope = lc.WaveRunner64MXsB()
        self.scope.numberofsegments = self.numberofsegments
        
        self.scope.connect()
        msg = self.scope.initialize([self.channel,self.pdchannel])
        self.Log.emit('%s'%msg)
        
        if self.scanmode == '3d-scan':
            self.Scan3d()
            self.StartScan()
            
        if self.scanmode == 'normal':
            self.Scan1d()
            self.StartScan()
            
        if self.scanmode == 'wavelength':
            self.Scan1d()
            self.StartScan()
            self.SetWaveLength.emit(float(self.guithread.WavelengthSB.value()))

            
        if self.scanmode == 'fake':
            self.ScanFake()
            
        if self.scanmode == 'free':
            self.scanchannels = []
            self.FreeRun()
            print('freee')
            
            
            




    

class GuiThread(QtWidgets.QMainWindow, gui.Ui_MainWindow):
    plotsequence = pyqtSignal(object,object)
    tof_scan = pyqtSignal(object,object,object,object, object,object, object,object,object,object, 
                         object,object,object,object)
    testfunc = pyqtSignal()

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        
        
        
        #new : read in parameters as class attributes
        
        filename = './input/HyT-input.dat'    #config file
        params, values=rf.readpars_class(filename)
        for i in range (len(params)):
            try:
                exec(params[i]+"=%s"%(values[i]))
            except (NameError, SyntaxError):
                if values[i] == '':
                    print('shid')
                    continue
                print(params[i], values[i])
                exec(params[i]+"='%s'"%(values[i].strip()))

        
        self.RunningStatusBox.setStyleSheet("QTextBrowser {background-color: rgb(255,0,0);}")

        # connect functions to signals (GUI)
        self.LoadBoardButton.clicked.connect(self.OnLoadBoard)
        self.SaveConfigButton.clicked.connect(self.OnSaveConfig)
        self.PlotSequenceButton.clicked.connect(self.OnPlotSeq)
        self.T2JumpButton.clicked.connect(self.OnT2Jump)
        self.ConfigButton.clicked.connect(self.OnConfig)
        self.StartTriggerButton.clicked.connect(self.OnStartTrigger)
        self.StopTriggerButton.clicked.connect(self.OnStopTrigger)
        self.ScanStopButton.clicked.connect(self.OnStopScan)
        self.Scan3dOpenButton.clicked.connect(self.OnScan3dFile)
        self.Scan3dOpenListButton.clicked.connect(self.OnScan3dListFile)
        self.ScopeUpdateButton.clicked.connect(self.OnScopeUpdate)
        self.StartScanButton.clicked.connect(self.OnScan)
#        self.ScanListStartButton.clicked.connect(self.OnListScan)


        # initialize the trigger control tabs as defined in config file 
        channels=12
        
        index=self.ch_combo_chno1.findText(ch1b)
        self.ch_combo_chno1.setCurrentIndex(index)
        
        index=self.ch_combo_chno2.findText(ch1e)
        self.ch_combo_chno2.setCurrentIndex(index)
        
        self.ch_time1_sbox.setValue(ch1t1)
        self.ch_time2_sbox.setValue(ch1t2)
        
        index=self.ch_time1_unitbox.findText(ch1u1)
        self.ch_time1_unitbox.setCurrentIndex(index)
        
        index=self.ch_time2_unitbox.findText(ch1u2)
        self.ch_time2_unitbox.setCurrentIndex(index)

        index=self.ch_pulsebox.findText(ch1pulse)
        self.ch_pulsebox.setCurrentIndex(index)
        
        self.ch_freqbox.setValue(ch1freq)
        
        index=self.ch_posneg.findText(ch1trig)
        self.ch_posneg.setCurrentIndex(index)
        
        self.ch_checkBox.setCheckState(ch1cb)
        
        self.T2jumpFile = T2jumpfile

        for chno in range(2,channels+1):

            try:
                index=eval("self.ch_combo_chno1_%i.findText(ch%ib)"%(chno,chno))
                exec('self.ch_combo_chno1_%i.setCurrentIndex(%i)'%(chno,index))
                
                index=eval("self.ch_combo_chno2_%i.findText(ch%ie)"%(chno,chno))
                exec('self.ch_combo_chno2_%i.setCurrentIndex(%i)'%(chno,index))
                
                exec("self.ch_time1_sbox_%i.setValue(ch%it1)"%(chno,chno))
                exec("self.ch_time2_sbox_%i.setValue(ch%it2)"%(chno,chno))
                
                #~ print(ch1u1)
                index=eval("self.ch_time1_unitbox_%i.findText(ch%iu1)"%(chno,chno))
                #~ print(index)
                exec("self.ch_time1_unitbox_%i.setCurrentIndex(%i)"%(chno,index))
                index=eval("self.ch_time2_unitbox_%i.findText(ch%iu2)"%(chno,chno))
                exec("self.ch_time2_unitbox_%i.setCurrentIndex(%i)"%(chno,index))
                
                index=eval("self.ch_pulsebox_%i.findText(ch%ipulse)"%(chno,chno))
                exec("self.ch_pulsebox_%i.setCurrentIndex(%i)"%(chno,index))
                
                exec("self.ch_freqbox_%i.setValue(ch%ifreq)"%(chno,chno))

                index=eval("self.ch_posneg_%i.findText(ch%itrig)"%(chno,chno))
                exec("self.ch_posneg_%i.setCurrentIndex(%i)"%(chno,index))
                
                exec("self.ch_checkBox_%i.setCheckState(ch%icb)"%(chno,chno))

            except NameError:
                #~ print('NameError')
                msg = 0
                
           
        self.T2JumpLine.setText(self.T2jumpFile)
        
        self.configfile = './input/HyT-input.dat'
        self.ConfigLine.setText(self.configfile)
        
        self.units={'ns':1e-9, 'us':1e-6, 'ms':1e-3, 's':1}


        index = self.ScanChannelBox.findText(scan_ch)
        self.ScanChannelBox.setCurrentIndex(index)
        
        index = self.ScanRelChannelBox.findText(rel_ch)
        self.ScanRelChannelBox.setCurrentIndex(index)
        
        index = self.ScanBGChannelBox.findText('%s'%bg_ch)
        self.ScanBGChannelBox.setCurrentIndex(index)


        self.ScanStartBox.setValue(scan_start)
        self.ScanStopBox.setValue(scan_stop)
        self.ScanStepBox.setValue(scan_step)
        self.AverageSBox.setValue(numberofsegments)
        
        index = self.ScanUnitBox.findText(scan_start_unit)
        self.ScanUnitBox.setCurrentIndex(index)
        
        self.Scan3dText.setText(scan3dfile)
        self.Scan3dListText.setText(listscanfile)
        
#        XXX todo scope settings
        index = self.ScopeTrigSource.findText(str(trig_source))
        self.ScopeTrigSource.setCurrentIndex(index)
        
        index = self.ScopeTrigEdge.findText(str(trig_slope))
        self.ScopeTrigEdge.setCurrentIndex(index)
        
        self.ScopeTrigLevel.setValue(trig_level)
        
        self.WavelengthSB.setValue(wavelength)

        channels = 4
        for ch in range (1,channels+1):
            exec('self.ScopeCB%i.setCheckState(ch%i_cb)'%(ch,ch))
            exec('index = self.ScopeCP%i.findText(ch%i_verticalcoupling)'%(ch,ch))
            exec('self.ScopeCP%i.setCurrentIndex(index)'%(ch))
            exec('self.ScopeCP%i_act.setText(str(ch%i_verticalcoupling))'%(ch,ch))
            exec('self.ScopeVScale%i.setValue(ch%i_scale)'%(ch,ch))
            exec('self.ScopeVScale%i_act.setText(str(ch%i_scale))'%(ch,ch))
            exec('self.ScopeVOffs%i.setValue(ch%i_verticaloffset)'%(ch,ch))
            exec('self.ScopeVOffs%i_act.setText(str(ch%i_verticaloffset))'%(ch,ch))
            exec('self.ScopeInv%i.setCheckState(ch%i_invert)'%(ch,ch))
            exec('self.ScopeInv%i_act.setText(str(ch%i_invert))'%(ch,ch))



            
        
        
        
        for frames in range(1,7):
            exec("self.ScanPlot%i = pg.PlotWidget()"%frames)
            exec("self.ScanPlot%i.setObjectName('ScanPlot%i')"%(frames,frames))
            exec("self.ScanPlotFrame%i.addWidget(self.ScanPlot%i)"%(frames,frames))
            exec("self.ScanPlt%i = self.ScanPlot%i.plot()"%(frames,frames))
            exec("self.ScanPlot%i.setLabel('bottom', text='time', units='s')"%(frames))

        
        
        self.start_working()
    
    def FailSafe(self, thread):
        choice = QtGui.QMessageBox.question(self, 'Failsafe',
                                    'Continue ?',
                                    QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        
        thread.Fail = choice==QtGui.QMessageBox.Yes
        return choice==QtGui.QMessageBox.Yes
     
    def TrigLog(self,msg):
        self.TriggerLog.append(msg)
        return
    
    def fScanLog(self,msg):
        self.ScanLog.append(msg)
        return
    
    def ScanBar(self,msg):
        self.scan_progressBar.setValue(msg)
        return
    
    def ScanTime(self,msg):
        self.scan_timeleftText.setText(msg)
        return
    
    def Freerun_Points(self, thread):
        thread.freerun_points = self.FreeRunPoints.value()
        return
    
    def OnScan(self):
        self.readboard_class()
        self.scan = ListenerThread(self)
#                
        self.scan.Plot.connect(self.PlotScan)
        self.scan.Log.connect(self.fScanLog)
        self.scan.Bar.connect(self.ScanBar)
        self.scan.ScanTime.connect(self.ScanTime)
        self.scan.FailSafe.connect(self.FailSafe)
        self.scan.SetWaveLength.connect(self.SetWaveLength)
        self.scan.Freerun_Points.connect(self.Freerun_Points)
        self.scan.start()
        return

        
    def OnStartTrigger(self):
        self.TriggerLog.append('Starting trigger.')
        self.RunningStatusBox.setStyleSheet("QTextBrowser {background-color: rgb(0,255,0);}")
        units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, t2jumpfile = self.readboard()
        
        self.TriggerLog.append('Starting Pulse Blaster!')
        scantime = 'None'
        scanchannel = 'None'
        relchannel = 'None'
        bgchannel = 'None'
        
        
        pulseseq = ps.sequence(units1, chno1, time1, units2, chno2, time2, freqs, 
                               pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, bgchannel,
                               scantime, t2jumpfile)
        pulseseq.seq()
        output = pulseseq.output
#        output,msg=ps.pulses(units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, scantime, t2jumpfile)
        status = pb.pulseblaster_program(output)
        time.sleep(0.1)
        pb.pulseblaster_start()        
        return
    
    def OnStartScan(self):
        self.TriggerLog.append('Starting scan.')
        
        self.readboard_class()
        
        self.scanmode = self.ScanModeBox.currentText()
        self.delaymode = self.DelayModeBox.currentText()
        self.bgmode = self.BGModeBox.currentText()
        self.repetitions = self.ScanRepsBox.value()
        self.numberofsegments = self.AverageSBox.value()
        
        self.scanchannel = 'No'
        self.relchannel = 'No'
        
        
        if self.scanmode == '3d-scan':
            self.ScanBar(0)
            self.Scan3d()
            return
        
        if self.scanmode == 'normal':
            self.scan = ListenerThread1d(self)
            
        elif self.scanmode == 'fake':
            self.scan = ListenerThreadFake(self)
#                
        self.scan.Plot.connect(self.PlotScan)
        self.scan.Log.connect(self.fScanLog)
        self.scan.Bar.connect(self.ScanBar)
        self.scan.ScanTime.connect(self.ScanTime)
        self.scan.FailSafe.connect(self.FailSafe)
        self.scan.start()
            
        
        time.sleep(1e-1)
            
            
        return
        
        
    def OnStopTrigger(self):
        if devmode==1:
            self.TriggerLog.append('No trigger running: devmode.')
        else:
            self.TriggerLog.append('Stopping trigger.')
            try:
                pb.pulseblaster_stop()
            except RuntimeError:
                self.TriggerLog.append('Trigger already stopped.')

        self.RunningStatusBox.setStyleSheet("QTextBrowser {background-color: rgb(255,0,0);}")
        return
        
    def start_working(self):
        self.TriggerLog.append('Welcome! You are using HyT-Control v.%s.\r'\
                                'Have fun!'%version)
        if devmode == 1:
            self.TriggerLog.append('Running in developer mode: No communication to devices.')
        
        return
    
    def OnScopeUpdate(self):
        scope = lc.WaveRunner64MXsB()
        scope.connect()
        channels = []
        
        for ch in range (1,5):
            if eval('self.ScopeCB%i.isChecked()'%ch):
                channels.append(ch)
                print(eval('scope.ch%i_verticaloffset'%(ch)))

                exec('ch%i_verticaloffset = self.ScopeVOffs%i.value()'%(ch,ch))
                exec('ch%i_verticalcoupling = self.ScopeCP%i.currentText()'%(ch,ch))
                exec('ch%i_invert = self.ScopeInv%i.isChecked()'%(ch,ch))
                exec('ch%i_scale = self.ScopeVScale%i.value()'%(ch,ch))

                exec("self.ScopeVOffs%i_act.setText(str(ch%i_verticaloffset))"%(ch,ch))
                exec('self.ScopeCP%i_act.setText(str(ch%i_verticalcoupling))'%(ch,ch))
                exec('self.ScopeInv%i_act.setText(str(ch%i_invert))'%(ch,ch))
                exec('self.ScopeVScale%i_act.setText(str(ch%i_scale))'%(ch,ch))
                
                exec('scope.ch%i_verticaloffset = self.ScopeVOffs%i.value()'%(ch,ch))
                exec('scope.ch%i_verticalcoupling = self.ScopeCP%i.currentText()'%(ch,ch))
                exec('scope.ch%i_invert = self.ScopeInv%i.isChecked()'%(ch,ch))
                exec('scope.ch%i_scale = self.ScopeVScale%i.value()'%(ch,ch))
                
                print(eval('scope.ch%i_verticaloffset'%(ch)))
                scope.numberofsegments = self.AverageSBox.value()

        scope.initialize(channels)
        scope.disconnect()

        return

        
    def OnT2Jump(self):
        t2jumpfile = QtWidgets.QFileDialog.getOpenFileName(self,'Open file', './', 'Data files (*.dat *.txt *out)')[0]
        self.TriggerLog.append('T2Jump file set to: %s'%t2jumpfile)
        self.T2JumpLine.setText(t2jumpfile)
        return
    
    def OnScan3dFile(self):
        scan3dfile = QtWidgets.QFileDialog.getOpenFileName(self,'Open file', './', 'Data files (*.dat *.txt *out)')[0]
        self.ScanLog.append('3d Scan file set to: %s'%scan3dfile)
        self.Scan3dText.setText(scan3dfile)
        return
    
    def OnScan3dListFile(self):
        scan3dlistfile = QtWidgets.QFileDialog.getOpenFileName(self,'Open file', './', 'Data files (*.dat *.txt *out)')[0]
        self.ScanLog.append('3d Scan List file set to: %s'%scan3dlistfile)
        self.Scan3dListText.setText(scan3dlistfile)
        return
    
    
    def OnListScanFile(self):
        listscanfile = QtWidgets.QFileDialog.getOpenFileName(self,'Open file', './', 'Data files (*.dat *.txt *out)')[0]
        self.ScanLog.append('1d listscan file set to: %s'%listscanfile)
        self.Scan3dText.setText(listscanfile)
        return
        
    def OnConfig(self):
        configfile = QtWidgets.QFileDialog.getOpenFileName(self,'Open file', './', 'Data files (*.dat *.txt)')[0]
        self.TriggerLog.append('Config file set to: %s'%configfile)
        self.ConfigLine.setText(configfile)
        return
        
    
    def OnWavelengthUpdate(self):
        self.wavelength = float(self.WavelengthSB.value())
        self.TrigLog('New wavelength: %s 1/cm'%self.wavelength)
        self.SetWaveLength(self.wavelength)
        return
    
    def SetWaveLength(self, wavelength):
        puls = pulsare.pulsare()
        puls.pulsare_ip = self.pulsare_ip
        puls.pulsare_port = self.pulsare_port
        puls.connect()
        time.sleep(1e-1)
        puls.setWavelength(wavelength)
        puls.disconnect()
        return
    
    
    def OnStopScan(self):
        try:   
            self.scan.f.close()
            self.scan.scan.terminate()

            self.scan.terminate()
            self.ScanLog.append('Scan terminated.')
        except:
            self.OnStopTrigger()
            self.ScanLog.append('no scan running, trigger stopped.')
            self.scan.terminate()

        return

        

    def OnPlotSeq(self):
        self.TriggerLog.append('Plot pulse sequence.')
        units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, t2jumpfile = self.readboard()
        scanchannel = 'None'
        relchannel = 'None'
        scantime = 'None'
        bgchannel = 'None'
        self.TriggerLog.append('Calculating sequence...')
#        print(ps.pulses(units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, scantime, t2jumpfile))
#        output,msg=ps.pulses(units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, scantime, t2jumpfile)
        pulseseq = ps.sequence(units1, chno1, time1, units2, chno2, time2, freqs, 
                               pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, bgchannel,
                               scantime, t2jumpfile)
        try:
            pulseseq.seq()
        except:
            print('bad seq')
            return
        output = pulseseq.output
        
#        if msg == 1:
        self.TriggerLog.append('Pulse sequence ok.')
#        if msg == 0:
#            self.TriggerLog.append('Pulse sequence wrong. Probably a delay longer than 1/frequency.')
#            self.TriggerLog.append(output)
#            return
        self.TriggerLog.append('Plotting...')
        
        self.PlotT = TrigPlotThread(output,checkboxvals)
        
        self.PlotT.TrigLog.connect(self.TrigLog)
        self.PlotT.PlotSeq.connect(self.PlotSeq)
        
        self.PlotT.start()
        
        return
        
    def bug(self):
        plt.show()
        sys.exit()
    
    def OnLoadBoard(self):
        units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, t2jumpfile = self.readboard()  ###pulsenature just says which trigger is what
        scanchannel = 'None'
        relchannel = 'None'
        scantime = 'None'
        bgchannel = 'None'
        
        
#        output = ps.pulses(units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, scantime, t2jumpfile)[0]
#        
        pulseseq = ps.sequence(units1, chno1, time1, units2, chno2, time2, freqs, 
                               pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, bgchannel,
                               scantime, t2jumpfile)
        pulseseq.seq()
        output = pulseseq.output
        
        if devmode == 0:
            ans = pb.pulseblaster_program(output)
            self.TriggerLog.append('PB status: %s.'%ans)
        return
        
    def PlotSeq(self,xvals,yvals,newlabels):
        fig, ax = plt.subplots()
        ticks = [0.25+i for i in range(len(xvals))]
        ax.set_yticks(ticks)
        ax.set_yticklabels(newlabels)
        plt.ylabel('Channel no.')
        plt.xlabel('time [$\mu$s]')
        
        no=0
        for i in range(len(xvals)):
            plt.plot(np.array(xvals[i])*1e6,(np.array(yvals[i])/2)-(newlabels[i]-1)+i)
        plt.show()
        
        
    def readboard_class(self):
        ''' read values from channels board'''
        
        self.t2jumpfile = self.T2JumpLine.text()
        lastindex=12
        indexes=[['%ib'%i,'%ie'%i] for i in range(1,lastindex)]
        #~ print(indexes)
        
        units={'ns':1e-9, 'us':1e-6, 'ms':1e-3, 's':1}
        
        #~ read all channels and times from gui
        self.units1=[]
        self.units1.append(self.ch_time1_unitbox.currentText())

        #~ time 1 units in board
        for i in range(2,lastindex+1):
            self.units1.append(eval('self.ch_time1_unitbox_'+str(i)+'.currentText()'))
            
        
        #~ channel number relative
        self.chno1=[]
        self.chno1.append(self.ch_combo_chno1.currentText())
        for i in range(2,lastindex+1):
            self.chno1.append(eval('self.ch_combo_chno1_'+str(i)+'.currentText()'))

        #~ time 1
        self.time1=[]
        self.time1.append(self.ch_time1_sbox.value()*units[self.units1[0]])
        for i in range(2,lastindex+1):
            self.time1.append(eval('self.ch_time1_sbox_'+str(i)+'.value()*units[self.units1['+str(i-1)+']]'))

        #~ time 2 units
        self.units2=[]
        self.units2.append(self.ch_time2_unitbox.currentText())
        for i in range(2,lastindex+1):
            self.units2.append(eval('self.ch_time2_unitbox_'+str(i)+'.currentText()'))
        
        #~ second relative channel
        self.chno2=[]
        self.chno2.append(self.ch_combo_chno2.currentText())
        for i in range(2,lastindex+1):
            self.chno2.append(eval('self.ch_combo_chno2_'+str(i)+'.currentText()'))
        
        #~ time 2
        self.time2=[]
        self.time2.append(self.ch_time2_sbox.value()*units[self.units2[0]])
        for i in range(2,lastindex+1):
            self.time2.append(eval('self.ch_time2_sbox_'+str(i)+'.value()*units[self.units2['+str(i-1)+']]'))
        
        #~ frequency of the pulse
        self.freqs=[]
        self.freqs.append(self.ch_freqbox.value())
        for i in range(2,lastindex+1):
            self.freqs.append(eval('self.ch_freqbox_'+str(i)+'.value()'))

        #~ which kind of pulse? decelerator sequence? normal pulse?
        self.pulsenature=[]
        self.pulsenature.append(self.ch_pulsebox.currentText())
        for i in range(2,lastindex+1):
            self.pulsenature.append(eval('self.ch_pulsebox_'+str(i)+'.currentText()'))
          
        #~ checkbox : active / inactive
        self.checkboxvals=[]
        self.checkboxvals.append(self.ch_checkBox.isChecked())
        for i in range(2,lastindex+1):
            self.checkboxvals.append(eval('self.ch_checkBox_'+str(i)+'.isChecked()'))
            
        self.posnegvals=[]
        self.posnegvals.append(self.ch_posneg.currentText())
        for i in range(2,lastindex+1):
            self.posnegvals.append(eval('self.ch_posneg_'+str(i)+'.currentText()'))

        return self.units1, self.chno1, self.time1, self.units2, self.chno2, self.time2, self.freqs, self.pulsenature, self.checkboxvals, self.posnegvals, self.t2jumpfile
        
    def readboard(self):
        ''' read values from channels board'''
        
        t2jumpfile = self.T2JumpLine.text()
        lastindex=12
        indexes=[['%ib'%i,'%ie'%i] for i in range(1,lastindex)]
        #~ print(indexes)
        
        units={'ns':1e-9, 'us':1e-6, 'ms':1e-3, 's':1}
        
        #~ read all channels and times from gui
        units1=[]
        units1.append(self.ch_time1_unitbox.currentText())

        #~ time 1 units in board
        for i in range(2,lastindex+1):
            units1.append(eval('self.ch_time1_unitbox_'+str(i)+'.currentText()'))
            
        
        #~ channel number relative
        chno1=[]
        chno1.append(self.ch_combo_chno1.currentText())
        for i in range(2,lastindex+1):
            chno1.append(eval('self.ch_combo_chno1_'+str(i)+'.currentText()'))

        #~ time 1
        time1=[]
        time1.append(self.ch_time1_sbox.value()*units[units1[0]])
        for i in range(2,lastindex+1):
            time1.append(eval('self.ch_time1_sbox_'+str(i)+'.value()*units[units1['+str(i-1)+']]'))

        #~ time 2 units
        units2=[]
        units2.append(self.ch_time2_unitbox.currentText())
        for i in range(2,lastindex+1):
            units2.append(eval('self.ch_time2_unitbox_'+str(i)+'.currentText()'))
        
        #~ second relative channel
        chno2=[]
        chno2.append(self.ch_combo_chno2.currentText())
        for i in range(2,lastindex+1):
            chno2.append(eval('self.ch_combo_chno2_'+str(i)+'.currentText()'))
        
        #~ time 2
        time2=[]
        time2.append(self.ch_time2_sbox.value()*units[units2[0]])
        for i in range(2,lastindex+1):
            time2.append(eval('self.ch_time2_sbox_'+str(i)+'.value()*units[units2['+str(i-1)+']]'))
        
        #~ frequency of the pulse
        freqs=[]
        freqs.append(self.ch_freqbox.value())
        for i in range(2,lastindex+1):
            freqs.append(eval('self.ch_freqbox_'+str(i)+'.value()'))

        #~ which kind of pulse? decelerator sequence? normal pulse?
        pulsenature=[]
        pulsenature.append(self.ch_pulsebox.currentText())
        for i in range(2,lastindex+1):
            pulsenature.append(eval('self.ch_pulsebox_'+str(i)+'.currentText()'))
          
        #~ checkbox : active / inactive
        checkboxvals=[]
        checkboxvals.append(self.ch_checkBox.isChecked())
        for i in range(2,lastindex+1):
            checkboxvals.append(eval('self.ch_checkBox_'+str(i)+'.isChecked()'))
            
        posnegvals=[]
        posnegvals.append(self.ch_posneg.currentText())
        for i in range(2,lastindex+1):
            posnegvals.append(eval('self.ch_posneg_'+str(i)+'.currentText()'))

        return units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, t2jumpfile
      
    def OnSaveConfig(self):
        filename_in = './input/HyT-input.dat'
        dialog = QtWidgets.QFileDialog.getSaveFileName(self, "Save file", "", "")
        #~ dialog = QtWidgets.QFileDialog.selectFile(self,'Open file', './', 'Data files (*.dat *.txt)')
        filename_out = dialog[0]
        try:
            self.saveconfig(filename_in,filename_out)
            msg = 'Output file: %s'%filename_out
        except FileNotFoundError:
            msg = 'File not found!'
        return msg
        
    def saveconfig(self, filename_in, filename_out):
        infile = open(filename_in, "r")
        outfile = open(filename_out, "w")
        
        #~ READ OUT ALL GUI VARIABLES ! (TO DO!)---------------------------------------------------
        
        #~ triggercontrol guiparams
        trigparams=[]
        chs=19
        for i in range(1,chs):
            trigparams.append('ch%ib'%i)
            trigparams.append('ch%ie'%i)
            trigparams.append('ch%ipulse'%i)
            trigparams.append('ch%ifreq'%i)
            trigparams.append('ch%itrig'%i)
            trigparams.append('ch%icb'%i)
            for j in range(1,3):
                trigparams.append('ch%it%i'%(i,j))
                trigparams.append('ch%iu%i'%(i,j))
                 
        chu1,chb,cht1,chu2,che,cht2,chfreq,chpulse,chcb,chtrig,t2jumpfile = self.readboard()
        units={'us':1e6, 'ms':1e3, 'ns':1e9, 's': 1}
        for i in range(len(chu1)):
            cht1[i] = cht1[i]*units[chu1[i]]
            cht2[i] = cht2[i]*units[chu2[i]]
        
        # rest of the gui changeable params in dictionary saveparams
        saveparams={}
        saveparams['T2jumpfile']= self.T2JumpLine.text()
        
        while True:
            line=infile.readline()
            paramname=''
            value=''
            num=0
            
            if line[0] is not '#':
                for char in line:
                    if char == '=':
                        break
                    paramname+=char
                    num+=1
                
                paramname=paramname.strip()
                value = eval(paramname)
                if paramname in trigparams:
                    tstr=''.join(i for i in paramname if not i.isdigit())
                    if paramname[-1].isdigit():
                        tstr+=paramname[-1]
                        
                    tnum=[int(i) for i in paramname if i.isdigit()]
                    value = eval(tstr)[tnum[0]-1]
                    try:
                        float(value)
                        value = "%.3f"%value
                    except (TypeError,ValueError):
                        value = value
                        
                if paramname in saveparams:
                    value = saveparams[paramname]
                    
                    
                outfile.write('%s = %s\r\n'%(paramname,value))
                
            else:
                outfile.write(line)
                
            if 'this is the end' in line:
                break
        return 
        
        
    def OnBrowse(self):
        dialog = QtWidgets.QFileDialog.getOpenFileName(self,'Open file', './', 'Data files (*.dat *.txt)')
        return dialog
    

    
    def PlotScan(self,xdata,ydata,frame):
        xdata = np.array(xdata)
        ydata = np.array(ydata)

        exec("self.ScanPlt%i.setData(xdata, ydata)"%frame)      
        return
        
def main():
    app = QtWidgets.QApplication(sys.argv)
    form = GuiThread()
    
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()
