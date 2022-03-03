#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QThread, pyqtSignal
from PyQt5.QtWidgets import QMessageBox
import sys
import readfile as rf
import gui
import pulseblaster as pb
# ~ import lcwr_scope as scope
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
from scipy.ndimage import gaussian_filter1d  # smoothing
from scipy.signal._peak_finding import find_peaks  # peak finding algorithm
import stepper_motor
import unidaq

# start by reading the config file
version = r'0.1 α'
devmode = 0  # if 1: developer mode which means no communication to other devices
units = {'ns': 1e-9, 'us': 1e-6, 'ms': 1e-3, 's': 1, 'steps': 1}

### read default parameters from HyT-input.dat
filename = './input/HyT-input.dat'  # config file
params, values = rf.readpars(filename)
for i in range(len(params)):
    try:
        exec(params[i] + "=%s" % (values[i]))
    except (NameError, SyntaxError):
        if values[i] == '':
            print('shid')
            continue
        print(params[i], values[i])
        exec(params[i] + "='%s'" % (values[i].strip()))


class TrigPlotThread(QThread):
    TrigLog = pyqtSignal(object)
    ScanLog = pyqtSignal(object)

    PlotSeq = pyqtSignal(object, object, object)

    def __init__(self, output, checkboxvals):
        QThread.__init__(self)
        self.output = output
        self.checkboxvals = checkboxvals
        self.plotsequence = plotseq.plotsequence

    def __del__(self):
        self.wait()

    def run(self):
        output = self.output
        checkboxvals = self.checkboxvals
        xvals, yvals, newlabels = self.plotsequence(output, checkboxvals)
        self.TrigLog.emit('Done')
        self.PlotSeq.emit(xvals, yvals, newlabels)
        # ~ self.terminate()
        return


class MotorShuttleThread(QThread):
    StepperMotorLog = pyqtSignal(object)
    ShuttleProgressBar = pyqtSignal(object)
    PositionLCD = pyqtSignal(object)

    def __init__(self, parent):
        QThread.__init__(self)
        self.guithread = parent
        self.my_motor = self.guithread.my_motor
        self.shuttle_reps = self.guithread.shuttle_reps
        self.delay_time = self.guithread.shuttle_delay_time

        self.shuttle_mode = self.guithread.shuttle_mode
        # self.oneway_steps = self.guithread.shuttle_oneway_steps
        # self.stopover_time = self.guithread.shuttle_stopover_time

    def __del__(self):
        self.wait()

    def shuttle(self):
        self.oneway_steps = self.guithread.shuttle_oneway_steps
        self.stopover_time = self.guithread.shuttle_stopover_time
        if self.shuttle_mode == "Test":
            self.my_motor.dwell(self.delay_time)
            shuttle_start_t = time.perf_counter()
        self.my_motor.rotate_neg(self.oneway_steps)
        # print("Rotate %i steps" % self.oneway_steps)
        # time.sleep(self.oneway_steps / self.my_motor.run_freq)
        time.sleep(self.my_motor.rotation_time(self.oneway_steps))
        self.PositionLCD.emit(int(self.my_motor.position))
        self.my_motor.dwell(self.stopover_time)
        self.my_motor.rotate_pos(int(self.oneway_steps))
        # time.sleep(self.oneway_steps / self.my_motor.run_freq)
        time.sleep(self.my_motor.rotation_time(self.oneway_steps))
        self.PositionLCD.emit(int(self.my_motor.position))
        if self.shuttle_mode == "Test":
            print("Time elapsed for a complete shuttle:", time.perf_counter() - shuttle_start_t)
        return True

    def shuttle_with_end_check(self):
        self.oneway_steps = self.guithread.shuttle_oneway_steps
        self.stopover_time = self.guithread.shuttle_stopover_time
        if self.shuttle_mode == "Test":
            self.my_motor.dwell(self.delay_time)

        shuttle_start_t = time.perf_counter()
        self.my_motor.rotate_neg(self.oneway_steps)
        # print("Rotate %i steps" % self.oneway_steps)
        # time.sleep(self.oneway_steps / self.my_motor.run_freq)
        time.sleep(self.my_motor.rotation_time(self.oneway_steps))
        self.PositionLCD.emit(self.my_motor.position)
        self.my_motor.dwell(self.stopover_time)
        # self.my_motor.set_run_freq(4000)
        self.my_motor.rotate_pos(self.oneway_steps)
        # time.sleep(self.oneway_steps / self.my_motor.run_freq)
        time.sleep(self.my_motor.rotation_time(self.oneway_steps))
        self.guithread.OnMotorEndCheck()

        shuttle_elapsed_t = time.perf_counter() - shuttle_start_t
        if self.shuttle_mode == "Test": print("Time elapsed for a complete shuttle:", shuttle_elapsed_t)
        tolerance_time = 2
        if shuttle_elapsed_t > 2 * (
                self.oneway_steps + 900) / self.my_motor.run_freq + self.stopover_time + tolerance_time:
            self.guithread.StepperMotorLog('-Warning: motor shuttle time longer than expected, might got stuck once.',
                                           red_text=True)
            return False
        self.my_motor.position = 0
        self.PositionLCD.emit(self.my_motor.position)
        return True

    def run(self):
        self.udaq = self.guithread.udaq
        hybrid_shuttle_with_both_modes = False  # CAREFUL HERE
        if self.shuttle_mode == "Test":
            for i in range(self.shuttle_reps):
                if i % 3 == 0 and hybrid_shuttle_with_both_modes:
                    last_successful_shuttle = self.shuttle_with_end_check()
                    if last_successful_shuttle == False:
                        break
                else:
                    if self.guithread.motor_end_check:
                        last_successful_shuttle = self.shuttle_with_end_check()
                        if last_successful_shuttle == False:
                            break
                    else:
                        self.shuttle()
                self.StepperMotorLog.emit('- Motor shuttle finished %d/%d' % (i + 1, self.shuttle_reps))
                self.ShuttleProgressBar.emit(round((i + 1) / self.shuttle_reps * 100, 8))

        else:
            old_value = 0
            while self.guithread.scan_is_running:
                error, value = self.udaq.ReadDIBit(0, 0, 0)
                if value != old_value and value == 1:
                    print("Motor shuttle trigger detected")
                    if self.guithread.motor_end_check == False:
                        shuttle_success = self.shuttle()
                    else:
                        shuttle_success = self.shuttle_with_end_check()
                    if not shuttle_success:
                        self.StepperMotorLog.emit('- Fatal error: motor shuttle failed, please stop.', red_text=True)
                        break
                old_value = value
        return


class DelayScanThread(QThread):
    '''
    This thread is just to overwrite the Pulsebalster during a "normal" scan.
    '''

    def __init__(self, units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals,
                 t2jumpfile,
                 scanchannel, relchannel, bgchannel, delay, numberofsegments, scope, channel):
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
        except Exception as err:
            print(err.args[0])
            print('pulseblaster not yet running')

        pulseseq = ps.sequence(self.units1, self.chno1, self.time1, self.units2,
                               self.chno2, self.time2, self.freqs, self.pulsenature,
                               self.checkboxvals, self.posnegvals, self.scanchannel,
                               self.relchannel, self.bgchannel, self.delay, self.t2jumpfile)

        if 'b' in self.scanchannel:
            index = int(self.scanchannel[:-1]) - 1
            pulseseq.chno1[index] = self.relchannel
            pulseseq.time1[index] = self.delay

        elif 'e' in self.scanchannel:
            index = int(self.scanchannel[:-1]) - 1
            pulseseq.chno2[index] = self.relchannel
            pulseseq.time2[index] = self.delay

        pulseseq.seq()
        output = pulseseq.output
        pb.pulseblaster_program(output)
        # time.sleep(1 / np.min(self.freqs) * 1.1) #!!!!!!!!!!!!
        time.sleep(110e-3)
        pb.pulseblaster_start()

        count = 0

        try:
            while count < self.numberofsegments:
                count = self.scope.segment_counter(self.channel)[0]
                # print(count)
                time.sleep(5e-2)
        except:
            print('tofscan, numberofsegments')
        return

    def run(self):
        self.tof_scan()
        return


class ListenerThread(QThread):
    '''
    This thread is to run a scan according to different scan modes.
    '''
    Plot = pyqtSignal(object, object, object)
    Log = pyqtSignal(object)
    Bar = pyqtSignal(object)
    SetWaveLength = pyqtSignal(object)
    ScanTime = pyqtSignal(object)
    FailSafe = pyqtSignal(object)
    # Freerun_Points = pyqtSignal(object)
    MotorProgressBar = pyqtSignal(object)

    def __init__(self, parent):
        QThread.__init__(self)
        self.guithread = parent
        self.numberofsegments = self.guithread.numberofsegments

        now = datetime.datetime.now()
        self.file_pre = '%04d-%02d-%02d-' % (now.year, now.month, now.day)
        self.file_end = '.h5'

        self.directory = 'C:/Users/Gruppe Willitsch/Desktop/data/%04d-%02d-%02d/' % (now.year, now.month, now.day)
        if not os.path.exists(self.directory):
            self.Log.emit('creating directory \n %s' % self.directory)
            os.makedirs(self.directory)

        self.count = 0

        for file in os.listdir(self.directory):
            if file.endswith('.h5'):
                number = int(file[-6:-3])
                if number >= self.count:
                    self.count = number

        self.scanmode = self.guithread.ScanModeBox.currentText()
        self.delaymode = self.guithread.DelayModeBox.currentText()
        self.bgmode = self.guithread.BGModeBox.currentText()
        self.repetitions = self.guithread.ScanRepsBox.value()
        self.numberofsegments = self.guithread.AverageSBox.value()

        self.scanchannel = 'No'
        self.relchannel = 'No'

        self.bgchannel = int(self.guithread.ScanBGChannelBox.currentText())

        self.channel = 1  # signal channel scope
        self.pdchannel = 2
        self.scopewait = 1e-2

        # self.x1 = 0.195e-6  # int time start
        self.x1 = 0.3e-6  # int time start
        self.x2 = 5.0e-6  # integration time end
        self.xb1 = self.x1  # background integration
        self.xb2 = 5.0e-6  #

        self.xpd1 = 2.5e-7  # probably integration times of the detector
        self.xpd2 = 3.5e-7
        self.xpdb1 = 3.5e-7
        self.xpdb2 = 6e-7

        # self.bgchannel = int(self.guithread.ScanBGChannelBox.currentText())
        # self.freerun_points = 0 # ??
        self.number = 0

        self.unit = units[self.guithread.ScanUnitBox.currentText()]

    def __del__(self):
        self.wait()

    def filename_update(self):
        self.count += 1
        self.filename = self.file_pre + '%03d' % (self.count) + self.file_end
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
                  'self.scanchannel', 'self.relchannel', 'self.bgmode', 'self.bgchannel',
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
                self.f.create_dataset('inputs/%s' % name, (len(asciiList), 1), 'S10',
                                      asciiList)  # , compression = 'gzip')

            except:
                asciiList = eval(inputz)
                self.f.create_dataset('inputs/%s' % name, data=asciiList)  # , compression = 'gzip')

        #

        #        print(self.scanchannels,self.relchannels, self.scantimes)
        for i in range(2):
            try:
                self.f.attrs['scan%i' % i] = self.scanchannels[i]
                self.f.attrs['scanrel%i' % i] = self.relchannels[i]

            except:
                # print('no scan%i'%i)
                a = 1
        self.f.close()

    def Scan3d(self):
        self.Log.emit('Starting grid Scan.')

        inputfile = self.guithread.Scan3dText.text()
        self.data = pd.read_csv(self.guithread.Scan3dText.text(), delim_whitespace=True, header=None, comment='#')

        self.Log.emit('Scan starts with inputfile %s' % (inputfile))

        self.values0 = np.arange(self.data[2][0], self.data[3][0] + 0.1 * self.data[4][0], self.data[4][0])
        self.values1 = np.arange(self.data[2][1], self.data[3][1] + 0.1 * self.data[4][1], self.data[4][1])

        if self.delaymode == 'list':
            inputfile = self.guithread.Scan3dListText.text()
            self.Log.emit('Using list file for delays: %s' % inputfile)
            self.values2 = np.array(pd.read_csv(inputfile, delim_whitespace=True, header=None, comment='#')[0])


        elif self.delaymode == 'normal':
            self.values2 = np.arange(self.data[2][2], self.data[3][2] + 0.1 * self.data[4][2], self.data[4][2])
            self.Log.emit('Using start, stop and step values from %s' % inputfile)

        self.vals1, self.vals2, self.vals3 = np.array(np.meshgrid(self.values1, self.values2, self.values0))

        self.scanchannels = np.array(self.data[0])
        self.relchannels = np.array(self.data[1])

        self.scanchannel = self.scanchannels[-1]
        self.relchannel = self.relchannels[-1]

    def Scan1d(self):
        self.Log.emit('Starting 1d Scan. %s' % self.scanmode)

        self.scanchannel = self.guithread.ScanChannelBox.currentText()
        self.relchannel = self.guithread.ScanRelChannelBox.currentText()

        if self.delaymode == 'list':
            inputfile = self.guithread.Scan3dListText.text()
            self.Log.emit('Using list file for delays: %s' % inputfile)

            self.values2 = np.array(pd.read_csv(inputfile, delim_whitespace=True, header=None, comment='#')[0])

        elif self.delaymode == 'normal':
            '''
            For normal scan:
            When scanning "motor stopover" or "motor steps", if the Q-switch (Channel 3) trigger is relative to 
            "12e", then the Q-switch timing will also be scanned with a delay wrt the end of every shuttle; otherwise 
            it is fixed.
            '''
            self.values2 = np.arange(self.guithread.ScanStartBox.value() * self.unit,
                                         self.guithread.ScanStopBox.value() * self.unit + 0.1 * self.guithread.ScanStepBox.value() * self.unit,
                                         self.guithread.ScanStepBox.value() * self.unit)

            if self.scanchannel == 'motor stopover' or self.scanchannel == 'motor steps':
                self.scanchannel = '12e'
                self.relchannel = '12b'

            self.Log.emit('Using start, stop and step values from GUI')

        if self.scanmode == 'wavelength':
            self.scanchannel = 'No'
            self.relchannel = 'No'
            self.values2 = np.arange(self.guithread.ScanStartBoxLaser.value(), self.guithread.ScanStopBoxLaser.value(),
                                     self.guithread.ScanStepBoxLaser.value())
            self.Log.emit('Wavelength Scan')

        self.values0 = [1e-6]  # some dummy values, not used in 1d scan!
        self.values1 = [1e-6]
        self.vals1, self.vals2, self.vals3 = np.array(np.meshgrid(self.values1, self.values2, self.values0))

        self.scanchannels = np.array(['No', 'No', self.scanchannel])
        self.relchannels = np.array(['No', 'No', self.relchannel])

    # def FreeRun(self):
    #     self.Log.emit('Start continuous mode')
    #     datax = []
    #     datay = []
    #
    #     delay = 0
    #
    #     datax_bg = []
    #     datay_bg = []
    #
    #     dataypd = []
    #     dataypd_bg = []
    #     count = 0
    #     while True:
    #         self.Freerun_Points.emit(self)
    #
    #         self.scope.clear_sweeps()
    #
    #         while count < self.numberofsegments:
    #             count = self.scope.segment_counter(self.channel)[0]
    #             time.sleep(5e-2)
    #         count = 0
    #         self.data = self.scope.get_waveform(self.channel)
    #         self.photodiode = self.scope.get_waveform(self.pdchannel)  # photodiode signal
    #
    #         ypdb = np.mean(
    #             self.photodiode[1][np.logical_and(self.photodiode[0] > self.xpdb1, self.photodiode[0] < self.xpdb2)])
    #         self.photodiode[1] = self.photodiode[1] - ypdb
    #
    #         y = np.mean(self.data[1][np.logical_and(self.data[0] > self.x1, self.data[0] < self.x2)])
    #         yb = np.mean(self.data[1][np.logical_and(self.data[0] > self.xb1, self.data[0] < self.xb2)])
    #         #                    datay.append(y-yb)
    #
    #         if len(datax) >= self.freerun_points:
    #             datax.pop(0)
    #             datay.pop(0)
    #             dataypd.pop(0)
    #         datax.append(delay)
    #
    #         datay.append(y)  # here might be what u are looking for
    #
    #         # photodiode
    #         ypd = np.sum(
    #             self.photodiode[1][np.logical_and((self.photodiode[0]) > self.xpd1, (self.photodiode[0]) < self.xpd2)])
    #         dataypd.append(ypd)
    #
    #         # plot
    #         self.Plot.emit(self.data[0], self.data[1], 1)
    #         self.Plot.emit(datax, datay, 2)
    #         self.Plot.emit(self.photodiode[0], self.photodiode[1], 5)
    #         delay += 1

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
        except Exception as err:
            print(err.args[0])
            print('pulseblaster not yet running')

        pointreps = self.guithread.ScanFakePointRepBox.value()
        totalreps = self.guithread.ScanFakeTotalRepBox.value()
        photon_counting_in_fake_mode = self.guithread.PhotonCountingCheckBox.isChecked()
        print("Photon counting mode:", photon_counting_in_fake_mode)

        count = 0

        self.scanchannel = self.guithread.ScanChannelBox.currentText()
        self.relchannel = self.guithread.ScanRelChannelBox.currentText()

        self.scanchannels = np.array(['No', 'No', self.scanchannel])
        self.relchannels = np.array(['No', 'No', self.relchannel])

        self.scantimes = ['No', 'No']

        self.filename_update()
        self.file_head3d()

        for total in range(totalreps):
            self.Fail = 'None'
            for point in range(pointreps):
                self.scope.clear_sweeps()
                if self.guithread.shuttle_mode == 'Trigger':
                    self.Log.emit('Motor Scanning Fake: point %i' % point + ' in position %i' % total)
                else:
                    self.Log.emit('Scanning Fake: point %i' % point + ' in position %i' % total)
                self.scan = DelayScanThread(self.guithread.units1, self.guithread.chno1, self.guithread.time1,
                                            self.guithread.units2,
                                            self.guithread.chno2, self.guithread.time2, self.guithread.freqs,
                                            self.guithread.pulsenature,
                                            self.guithread.checkboxvals, self.guithread.posnegvals,
                                            self.guithread.t2jumpfile,
                                            'No', 'No', 'No',
                                            delay, self.numberofsegments,
                                            self.scope, self.channel)
                self.scan.start()

                while self.scan.isRunning():
                    time.sleep(self.scopewait)

                datax.append(count)
                self.data = self.scope.get_waveform(self.channel)
                self.photodiode = self.scope.get_waveform(self.pdchannel)  # photodiode signal

                ypdb = np.mean(self.photodiode[1][
                                   np.logical_and(self.photodiode[0] > self.xpdb1, self.photodiode[0] < self.xpdb2)])
                self.photodiode[1] = self.photodiode[1] - ypdb

                y = np.mean(self.data[1][np.logical_and(self.data[0] > self.x1, self.data[0] < self.x2)])
                yb = np.mean(self.data[1][np.logical_and(self.data[0] > self.xb1, self.data[0] < self.xb2)])
                #                    datay.append(y-yb)

                # y will contain at the end the number of photons per laser shot detected.
                # data[0] contains the timebase
                # data[1] contains the signal from the oscilloscope
                # x1, x2 are the gating times for the signal, currently 1.95e-7, 5e-6
                # xb1, xb2 are the gating times from the background, currently 1.95e-7, 5e-6
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
                    y = len(signal_peaks) / self.numberofsegments
                    yb = len(background_peaks) / self.numberofsegments
                    # print('number of photons y :', y)
                    # print('number of photons yb:', yb)

                # print('signal peaks', signal_peaks)
                datay.append(y)  # here might be what u are looking for

                # photodiode
                ypd = np.sum(self.photodiode[1][
                                 np.logical_and((self.photodiode[0]) > self.xpd1, (self.photodiode[0]) < self.xpd2)])
                dataypd.append(ypd)

                # plot
                self.Plot.emit(self.data[0], self.data[1], 1)
                self.Plot.emit(datax, datay, 2)
                self.Plot.emit(self.photodiode[0], self.photodiode[1], 5)
                #                    self.Plot.emit(datax,dataypd,6)

                self.f = h5py.File(self.directory + self.filename, 'a')
                dset_s = self.f.create_dataset('signal/' + '%03d-' % total + '%03d' % point,
                                               data=np.array(self.data[1]), compression='gzip')
                dset_s.attrs['x1'] = self.x1
                dset_s.attrs['x2'] = self.x2
                dset_s.attrs['y'] = y
                dset_s.attrs['xb1'] = self.xb1
                dset_s.attrs['xb2'] = self.xb2
                dset_s.attrs['y_oc'] = y - yb

                pd_bd = np.histogram(self.photodiode[0], weights=np.array(self.photodiode[1]), bins=200)

                # photodiode
                dset_s = self.f.create_dataset('signal_pd/' + '%03d-' % total + '%03d' % point, data=np.array(pd_bd[0]),
                                               compression='gzip')
                dset_s.attrs['xpd1'] = self.xpd1
                dset_s.attrs['xpd2'] = self.xpd2
                dset_s.attrs['ypd'] = ypd
                dset_s.attrs['xpdb1'] = self.xpdb1
                dset_s.attrs['xpdb2'] = self.xpdb2

                # background scan
                if self.bgmode != "None":
                    self.scan = DelayScanThread(self.guithread.units1, self.guithread.chno1, self.guithread.time1,
                                                self.guithread.units2,
                                                self.guithread.chno2, self.guithread.time2, self.guithread.freqs,
                                                self.guithread.pulsenature,
                                                self.guithread.checkboxvals, self.guithread.posnegvals,
                                                self.guithread.t2jumpfile,
                                                'No', 'No', self.bgchannel,
                                                delay, self.numberofsegments,
                                                self.scope, self.channel)

                    self.scope.clear_sweeps()
                    self.scan.start()

                    while self.scan.isRunning():
                        time.sleep(self.scopewait)

                    self.data = self.scope.get_waveform(self.channel)
                    self.photodiode = self.scope.get_waveform(self.pdchannel)  # photodiode signal
                    ypdb = np.mean(self.photodiode[1][
                                       np.logical_and(self.photodiode[0] > self.xpdb1,
                                                      self.photodiode[0] < self.xpdb2)])
                    self.photodiode[1] = self.photodiode[1] - ypdb

                    datax_bg.append(count)
                    y = np.mean(self.data[1][np.logical_and(self.data[0] > self.x1, self.data[0] < self.x2)])
                    yb = np.mean(self.data[1][np.logical_and(self.data[0] > self.xb1, self.data[0] < self.xb2)])
                    #                    datay_bg.append(y-yb)
                    datay_bg.append(y)  # super XXX

                    ypd = np.sum(self.photodiode[1][
                                     np.logical_and((self.photodiode[0]) > self.xpd1,
                                                    (self.photodiode[0]) < self.xpd2)])
                    dataypd_bg.append(ypd)

                    self.Plot.emit(self.data[0], self.data[1], 1)
                    self.Plot.emit(datax, datay_bg, 3)
                    self.Plot.emit(np.array(datax), np.array(datay) - np.array(datay_bg), 4)

                    self.Plot.emit(self.photodiode[0], self.photodiode[1], 5)
                    self.Plot.emit(datax, dataypd_bg, 6)

                    dset_bg = self.f.create_dataset('background/' + '%03d-' % total + '%03d' % point,
                                                    data=np.array(self.data[1]), compression='gzip')
                    dset_bg.attrs['x1'] = self.x1
                    dset_bg.attrs['x2'] = self.x2
                    dset_bg.attrs['y'] = y
                    dset_bg.attrs['xb1'] = self.xb1
                    dset_bg.attrs['xb2'] = self.xb2
                    dset_bg.attrs['y_oc'] = y - yb

                    # photodiode
                    pd_bd = np.histogram(self.photodiode[0], weights=np.array(self.photodiode[1]), bins=200)

                    dset_s = self.f.create_dataset('background_pd/' + '%03d-' % total + '%03d' % point,
                                                   data=np.array(pd_bd[0]), compression='gzip')
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
        letime = 2 * np.shape(self.vals2.T[0][0])[0] / np.min(self.guithread.freqs) * self.numberofsegments * (
                self.totallength - self.number)
        day = letime // (24 * 3600)
        letime = letime % (24 * 3600)
        hour = letime // 3600
        letime %= 3600
        minutes = letime // 60
        letime %= 60
        seconds = letime
        letime = "dd:hh:mm:ss-> %02d:%02d:%02d:%02d" % (day, hour, minutes, seconds)

        self.ScanTime.emit('Approx. time left: %s' % (letime))

        for repetition in range(self.repetitions):
            for index0 in range(self.vals3.T.shape[0]):
                for index1 in range(self.vals3.T.shape[1]):
                    scanval1 = self.vals3.T[index0][index1][0]
                    scanval2 = self.vals1.T[index0][index1][0]
                    self.scantimes = [scanval1, scanval2]

                    delays = self.vals2.T[index0][index1]

                    if self.guithread.ScanModeBox.currentText() == '3d-scan':
                        for index in range(2):
                            chno1r = int(self.scanchannels[index][:-1])

                            if self.scanchannels[index][-1] == 'b':
                                self.guithread.chno1[chno1r - 1] = self.relchannels[index]
                                self.guithread.time1[chno1r - 1] = self.scantimes[index]

                            elif self.scanchannels[index][-1] == 'e':
                                self.guithread.chno2[chno1r - 1] = self.relchannels[index]
                                self.guithread.time2[chno1r - 1] = self.scantimes[index]

                    for file in os.listdir(self.directory):
                        if file.endswith('.h5'):
                            number = int(file[-6:-3])
                            if number >= self.count:
                                self.count = number

                    self.filename_update()
                    self.file_head3d()
                    self.f = h5py.File(self.directory + self.filename, 'a')

                    for i in range(2):
                        self.f.attrs['scan%i_time' % i] = self.scantimes[i]
                    self.f.close()

                    self.Log.emit(
                        'Scanning %s rel to %s with bg %s' % (self.scanchannel, self.relchannel, self.bgchannel))
                    self.Log.emit('%s' % self.scanchannels)
                    self.Log.emit('%s' % self.relchannels)
                    self.Log.emit('%s' % self.scantimes)
                    datax = []
                    datay = []

                    datax_bg = []
                    datay_bg = []

                    dataypd = []
                    dataypd_bg = []

                    self.number += 1

                    for delay_counter, delay in enumerate(delays):
                        now = datetime.datetime.now()
                        tstamp = '%04d-%02d-%02d-%02d:%02d%02d' % (
                            now.year, now.month, now.day, now.hour, now.minute, now.second)

                        #signal scan
                        # time0 = time.clock()
                        self.scope.clear_sweeps()
                        if self.scanmode == 'wavelength':
                            self.SetWaveLength.emit(delay)
                            time.sleep(0.5)

                        if self.guithread.ScanChannelBox.currentText() == 'motor stopover':
                            scanchannel_delay = self.guithread.my_motor.shuttle_total_time(self.guithread.shuttle_oneway_steps, delay)
                            self.guithread.shuttle_stopover_time = delay

                        elif self.guithread.ScanChannelBox.currentText() =='motor steps':
                            scanchannel_delay = self.guithread.my_motor.shuttle_total_time(int(delay), self.guithread.shuttle_stopover_time)
                            self.guithread.shuttle_oneway_steps = int(delay)
                        else:
                            scanchannel_delay = delay

                        self.scan = DelayScanThread(self.guithread.units1, self.guithread.chno1, self.guithread.time1,
                                                    self.guithread.units2,
                                                    self.guithread.chno2, self.guithread.time2, self.guithread.freqs,
                                                    self.guithread.pulsenature,
                                                    self.guithread.checkboxvals, self.guithread.posnegvals,
                                                    self.guithread.t2jumpfile,
                                                    self.scanchannel, self.relchannel, 'None',
                                                    scanchannel_delay, self.numberofsegments,
                                                    self.scope, self.channel)
                        self.scan.start()

                        # Wait for the Pulseblaster to output pulses according to the number of averages
                        while self.scan.isRunning():
                            time.sleep(self.scopewait)
                            # pass

                        # if 'motor' in self.guithread.ScanChannelBox.currentText():
                        #     self.motor_shuttle_thread.terminate()

                        datax.append(delay)

                        self.data = self.scope.get_waveform(self.channel)
                        self.photodiode = self.scope.get_waveform(self.pdchannel)  # photodiode signal

                        ypdb = np.mean(self.photodiode[1][np.logical_and(self.photodiode[0] > self.xpdb1,
                                                                         self.photodiode[0] < self.xpdb2)])
                        self.photodiode[1] = self.photodiode[1] - ypdb

                        y = np.mean(self.data[1][np.logical_and(self.data[0] > self.x1, self.data[0] < self.x2)])
                        yb = np.mean(self.data[1][np.logical_and(self.data[0] > self.xb1, self.data[0] < self.xb2)])
                        # datay.append(y-yb)
                        datay.append(y)  # here might be what u are looking for

                        # photodiode
                        ypd = np.sum(self.photodiode[1][np.logical_and((self.photodiode[0]) > self.xpd1,
                                                                       (self.photodiode[0]) < self.xpd2)])
                        dataypd.append(ypd)

                        # plot
                        self.Plot.emit(self.data[0], self.data[1], 1)
                        self.Plot.emit(datax, datay, 2)
                        self.Plot.emit(self.photodiode[0], self.photodiode[1], 5)
                        #                    self.Plot.emit(datax,dataypd,6)

                        self.f = h5py.File(self.directory + self.filename, 'a')
                        dset_s = self.f.create_dataset('signal/' + '%s' % round(delay, 8), data=np.array(self.data[1]),
                                                       compression='gzip')
                        dset_s.attrs['x1'] = self.x1
                        dset_s.attrs['x2'] = self.x2
                        dset_s.attrs['y'] = y
                        dset_s.attrs['xb1'] = self.xb1
                        dset_s.attrs['xb2'] = self.xb2
                        dset_s.attrs['y_oc'] = y - yb
                        dset_s.attrs['tstamp'] = tstamp

                        if delay == delays[0]:
                            self.f.attrs['osci_tstart'] = self.data[0][0]
                            self.f.attrs['osci_tstop'] = self.data[0][-1]
                            self.f.attrs['osci_tbins'] = len(self.data[0])

                        pd_bd = np.histogram(self.photodiode[0], weights=np.array(self.photodiode[1]), bins=200)

                        # photodiode
                        dset_s = self.f.create_dataset('signal_pd/' + '%s' % round(delay, 8), data=np.array(pd_bd[0]),
                                                       compression='gzip')
                        dset_s.attrs['xpd1'] = self.xpd1
                        dset_s.attrs['xpd2'] = self.xpd2
                        dset_s.attrs['ypd'] = ypd
                        dset_s.attrs['xpdb1'] = self.xpdb1
                        dset_s.attrs['xpdb2'] = self.xpdb2

                        # begin background scan
                        if self.bgmode != "None":
                            now = datetime.datetime.now()
                            tstamp = '%04d-%02d-%02d-%02d:%02d%02d' % (
                                now.year, now.month, now.day, now.hour, now.minute, now.second)
                            self.scope.clear_sweeps()
                            self.scan = DelayScanThread(self.guithread.units1, self.guithread.chno1,
                                                        self.guithread.time1,
                                                        self.guithread.units2,
                                                        self.guithread.chno2, self.guithread.time2,
                                                        self.guithread.freqs,
                                                        self.guithread.pulsenature,
                                                        self.guithread.checkboxvals, self.guithread.posnegvals,
                                                        self.guithread.t2jumpfile,
                                                        self.scanchannel, self.relchannel, self.bgchannel,
                                                        scanchannel_delay, self.numberofsegments,
                                                        self.scope, self.channel)

                            self.scan.start()

                            while self.scan.isRunning():
                                time.sleep(self.scopewait)

                            self.data = self.scope.get_waveform(self.channel)
                            self.photodiode = self.scope.get_waveform(self.pdchannel)  # photodiode signal
                            ypdb = np.mean(self.photodiode[1][np.logical_and(self.photodiode[0] > self.xpdb1,
                                                                             self.photodiode[0] < self.xpdb2)])
                            self.photodiode[1] = self.photodiode[1] - ypdb

                            datax_bg.append(delay)
                            y = np.mean(self.data[1][np.logical_and(self.data[0] > self.x1, self.data[0] < self.x2)])
                            yb = np.mean(self.data[1][np.logical_and(self.data[0] > self.xb1, self.data[0] < self.xb2)])
                            #                    datay_bg.append(y-yb)
                            datay_bg.append(y)  # super XXX

                            ypd = np.sum(self.photodiode[1][np.logical_and((self.photodiode[0]) > self.xpd1,
                                                                           (self.photodiode[0]) < self.xpd2)])
                            dataypd_bg.append(ypd)

                            self.Plot.emit(self.data[0], self.data[1], 1)
                            self.Plot.emit(datax_bg, datay_bg, 3)
                            self.Plot.emit(np.array(datax), np.array(datay) - np.array(datay_bg), 4)

                            self.Plot.emit(self.photodiode[0], self.photodiode[1], 5)
                            self.Plot.emit(datax, dataypd_bg, 6)

                            dset_bg = self.f.create_dataset('background/' + '%s' % round(delay, 8),
                                                            data=np.array(self.data[1]), compression='gzip')
                            dset_bg.attrs['x1'] = self.x1
                            dset_bg.attrs['x2'] = self.x2
                            dset_bg.attrs['y'] = y
                            dset_bg.attrs['xb1'] = self.xb1
                            dset_bg.attrs['xb2'] = self.xb2
                            dset_bg.attrs['y_oc'] = y - yb
                            dset_s.attrs['tstamp'] = tstamp

                            # photodiode
                            pd_bd = np.histogram(self.photodiode[0], weights=np.array(self.photodiode[1]), bins=200)

                            dset_s = self.f.create_dataset('background_pd/' + '%s' % round(delay, 8),
                                                           data=np.array(pd_bd[0]), compression='gzip')
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

                    # time1 = time.clock()
                    #                letime = (time1-time0)*(self.totallength-self.number)
                    letime = 2 * np.shape(delays)[0] / np.min(self.guithread.freqs) * self.numberofsegments * (
                            self.totallength - self.number)
                    day = letime // (24 * 3600)
                    letime = letime % (24 * 3600)
                    hour = letime // 3600
                    letime %= 3600
                    minutes = letime // 60
                    letime %= 60
                    seconds = letime

                    letime = "dd:hh:mm:ss-> %02d:%02d:%02d:%02d" % (day, hour, minutes, seconds)

                    self.ScanTime.emit('Approx. time left: %s' % (letime))
                    self.Bar.emit(round(self.number / self.totallength * 100, 8))

    def run(self):
        self.scope = lc.WaveRunner64MXsB()
        self.scope.numberofsegments = self.numberofsegments
        self.scope.connect()
        msg = self.scope.initialize([self.channel, self.pdchannel])
        self.Log.emit('%s' % msg)

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

        # if self.scanmode == 'motor fake':
        #     self.ScanFake()

        self.guithread.scan_is_running = False
        return


class GuiThread(QtWidgets.QMainWindow, gui.Ui_MainWindow):
    plotsequence = pyqtSignal(object, object)
    tof_scan = pyqtSignal(object, object, object, object, object, object, object, object, object, object,
                          object, object, object, object)
    testfunc = pyqtSignal()

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)

        self.start_working()

        self.sigchannel_label.setVisible(False)
        self.Scan3dStartButton.setVisible(False)
        self.ScanFakeStartButton.setVisible(False)
        self.FreeRunPoints.setVisible(False)
        self.label_12.setVisible(False)

        # new : read in parameters as class attributes

        filename = './input/HyT-input.dat'  # config file
        params, values = rf.readpars_class(filename)
        for i in range(len(params)):
            try:
                exec(params[i] + "=%s" % (values[i]))
            except (NameError, SyntaxError):
                if values[i] == '':
                    print('shid')
                    continue
                print(params[i], values[i])
                exec(params[i] + "='%s'" % (values[i].strip()))

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
        self.WavelengthButton.clicked.connect(self.OnWavelengthUpdate)
        self.WavelengthUpButton.clicked.connect(self.OnWavelengthUp)
        self.WavelengthDownButton.clicked.connect(self.OnWavelengthDown)
        #        self.ScanListStartButton.clicked.connect(self.OnListScan)
        self.LinkedFreqRadioButton.setChecked(True)
        self.ch_freqbox.valueChanged.connect(self.OnChannelFreqChanged)


        # initialize the trigger control tabs as defined in config file
        channels = 12

        index = self.ch_combo_chno1.findText(ch1b)
        self.ch_combo_chno1.setCurrentIndex(index)

        index = self.ch_combo_chno2.findText(ch1e)
        self.ch_combo_chno2.setCurrentIndex(index)

        self.ch_time1_sbox.setValue(ch1t1)
        self.ch_time2_sbox.setValue(ch1t2)

        index = self.ch_time1_unitbox.findText(ch1u1)
        self.ch_time1_unitbox.setCurrentIndex(index)

        index = self.ch_time2_unitbox.findText(ch1u2)
        self.ch_time2_unitbox.setCurrentIndex(index)

        index = self.ch_pulsebox.findText(ch1pulse)
        self.ch_pulsebox.setCurrentIndex(index)

        self.ch_freqbox.setValue(ch1freq)

        index = self.ch_posneg.findText(ch1trig)
        self.ch_posneg.setCurrentIndex(index)

        self.ch_checkBox.setCheckState(ch1cb)

        self.T2jumpFile = T2jumpfile

        for chno in range(2, channels + 1):

            try:
                index = eval("self.ch_combo_chno1_%i.findText(ch%ib)" % (chno, chno))
                exec('self.ch_combo_chno1_%i.setCurrentIndex(%i)' % (chno, index))

                index = eval("self.ch_combo_chno2_%i.findText(ch%ie)" % (chno, chno))
                exec('self.ch_combo_chno2_%i.setCurrentIndex(%i)' % (chno, index))

                exec("self.ch_time1_sbox_%i.setValue(ch%it1)" % (chno, chno))
                exec("self.ch_time2_sbox_%i.setValue(ch%it2)" % (chno, chno))

                # ~ print(ch1u1)
                index = eval("self.ch_time1_unitbox_%i.findText(ch%iu1)" % (chno, chno))
                # ~ print(index)
                exec("self.ch_time1_unitbox_%i.setCurrentIndex(%i)" % (chno, index))
                index = eval("self.ch_time2_unitbox_%i.findText(ch%iu2)" % (chno, chno))
                exec("self.ch_time2_unitbox_%i.setCurrentIndex(%i)" % (chno, index))

                index = eval("self.ch_pulsebox_%i.findText(ch%ipulse)" % (chno, chno))
                exec("self.ch_pulsebox_%i.setCurrentIndex(%i)" % (chno, index))

                exec("self.ch_freqbox_%i.setValue(ch%ifreq)" % (chno, chno))

                index = eval("self.ch_posneg_%i.findText(ch%itrig)" % (chno, chno))
                exec("self.ch_posneg_%i.setCurrentIndex(%i)" % (chno, index))

                exec("self.ch_checkBox_%i.setCheckState(ch%icb)" % (chno, chno))

            except NameError:
                msg = 0

        self.InitMotor()
        self.T2JumpLine.setText(self.T2jumpFile)

        self.configfile = './input/HyT-input.dat'
        self.ConfigLine.setText(self.configfile)

        self.units = {'ns': 1e-9, 'us': 1e-6, 'ms': 1e-3, 's': 1}

        index = self.ScanChannelBox.findText(scan_ch)
        self.ScanChannelBox.setCurrentIndex(index)

        index = self.ScanRelChannelBox.findText(rel_ch)
        self.ScanRelChannelBox.setCurrentIndex(index)

        index = self.ScanBGChannelBox.findText('%s' % bg_ch)
        self.ScanBGChannelBox.setCurrentIndex(index)

        self.ScanStartBox.setValue(scan_start)
        self.ScanStopBox.setValue(scan_stop)
        self.ScanStepBox.setValue(scan_step)
        self.AverageSBox.setValue(numberofsegments)

        index = self.ScanUnitBox.findText(scan_start_unit)
        self.ScanUnitBox.setCurrentIndex(index)

        self.Scan3dText.setText(scan3dfile)
        self.Scan3dListText.setText(listscanfile)

        self.PhotonCountingCheckBox.setChecked(photon_counting)
        self.ScanFakePointRepBox.setValue(fakescan_pointreps)
        self.ScanFakeTotalRepBox.setValue(fakescan_totalreps)

        #        XXX todo scope settings
        index = self.ScopeTrigSource.findText(str(trig_source))
        self.ScopeTrigSource.setCurrentIndex(index)

        index = self.ScopeTrigEdge.findText(str(trig_slope))
        self.ScopeTrigEdge.setCurrentIndex(index)

        self.ScopeTrigLevel.setValue(trig_level)

        self.WavelengthSB.setValue(wavelength)
        self.WavelengthIncrementBox.setValue(wavelength_increment)

        self.define_pulsare = False  # set to True if u're ready to do wavelength scans, otherwise False
        if self.define_pulsare:
            # Pietro here, summer 2021
            # I move pulsare to the init.
            # this may crash eveything if not TCT connection can be established, but avoids
            # to init the connection at every scan.
            self.my_pulsare = pulsare.pulsare()  # in the init there is also the connect method and the setWavelength method,
            # that connect and initialize the laser to the default wavelength
            self.my_pulsare.getActualPosition()

        channels = 4
        for ch in range(1, channels + 1):
            exec('self.ScopeCB%i.setCheckState(ch%i_cb)' % (ch, ch))
            exec('index = self.ScopeCP%i.findText(ch%i_verticalcoupling)' % (ch, ch))
            exec('self.ScopeCP%i.setCurrentIndex(index)' % (ch))
            exec('self.ScopeCP%i_act.setText(str(ch%i_verticalcoupling))' % (ch, ch))
            exec('self.ScopeVScale%i.setValue(ch%i_scale)' % (ch, ch))
            exec('self.ScopeVScale%i_act.setText(str(ch%i_scale))' % (ch, ch))
            exec('self.ScopeVOffs%i.setValue(ch%i_verticaloffset)' % (ch, ch))
            exec('self.ScopeVOffs%i_act.setText(str(ch%i_verticaloffset))' % (ch, ch))
            exec('self.ScopeInv%i.setCheckState(ch%i_invert)' % (ch, ch))
            exec('self.ScopeInv%i_act.setText(str(ch%i_invert))' % (ch, ch))

        for frames in range(1, 7):
            exec("self.ScanPlot%i = pg.PlotWidget()" % frames)
            exec("self.ScanPlot%i.setObjectName('ScanPlot%i')" % (frames, frames))
            exec("self.ScanPlotFrame%i.addWidget(self.ScanPlot%i)" % (frames, frames))
            exec("self.ScanPlt%i = self.ScanPlot%i.plot()" % (frames, frames))
            exec("self.ScanPlot%i.setLabel('bottom', text='time', units='s')" % (frames))

    def FailSafe(self, thread):
        choice = QtGui.QMessageBox.question(self, 'Failsafe',
                                            'Continue ?',
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)

        thread.Fail = choice == QtGui.QMessageBox.Yes
        return choice == QtGui.QMessageBox.Yes

    def TrigLog(self, msg, red_text=False):
        if red_text:
            self.TriggerLog.setTextColor(QtGui.QColor("red"))
            self.TriggerLog.append(msg)
            self.TriggerLog.setTextColor(QtGui.QColor("black"))
        else:
            self.TriggerLog.setTextColor(QtGui.QColor("black"))
            self.TriggerLog.append(msg)
        return

    def fScanLog(self, msg, red_text=False):
        if red_text:
            self.ScanLog.setTextColor(QtGui.QColor("red"))
            self.ScanLog.append(msg)
            self.ScanLog.setTextColor(QtGui.QColor("black"))
        else:
            self.TriggerLog.setTextColor(QtGui.QColor("black"))
            self.ScanLog.append(msg)
        return

    def StepperMotorLog(self, msg, red_text=False):
        if red_text:
            self.MotorLog.setTextColor(QtGui.QColor("red"))
            self.MotorLog.append(msg)
            self.MotorLog.setTextColor(QtGui.QColor("black"))
        else:
            self.TriggerLog.setTextColor(QtGui.QColor("black"))
            self.MotorLog.append(msg)
        return

    def OnChannelFreqChanged(self):
        channels = 12
        if self.LinkedFreqRadioButton.isChecked():
            for chno in range(3, channels + 1):
                exec("self.ch_freqbox_%i.setValue(%s)" % (chno, self.ch_freqbox.value()))

    def ScanBar(self, msg):
        self.scan_progressBar.setValue(msg)
        return

    def ScanTime(self, msg):
        self.scan_timeleftText.setText(msg)
        return

    def Freerun_Points(self, thread):
        thread.freerun_points = self.FreeRunPoints.value()
        return

    def OnScan(self):
        # if not pb.pulseblaster_status()['running'] or pb.pulseblaster_status()['stopped']:
        #     self.fScanLog('Error: pulseblaster not yet running!', red_text=True)
        #     return
        self.TriggerTab.setDisabled(True)
        self.OnStopTrigger()
        self.readboard_class()

        time_max = max([max(self.time1[i], self.time2[i]) for i in range(len(self.chno1)) if self.checkboxvals[i]])
        if self.ScanChannelBox.currentText() == 'motor steps':
            time_max=max(time_max, self.my_motor.shuttle_total_time(oneway_steps=self.ScanStopBox.value(),
                                                                    stopover_time=self.shuttle_stopover_time))
        if self.ScanChannelBox.currentText() == 'motor stopover':
            time_max=max(time_max, self.my_motor.shuttle_total_time(oneway_steps=self.shuttle_oneway_steps,
                                                                    stopover_time=self.ScanStopBox.value() * units[self.ScanUnitBox.currentText()]))
        if max(self.freqs[2:]) > 1.0/(time_max + 0.005):
            self.fScanLog('Error: Trigger frequency too high!\r - switch motor shuttle mode to Test if not needed\r - or set to recommended value: < %.3f Hz' % (1.0/(time_max + 0.005)),
                red_text=True)
            return

        # if self.ScanModeBox.currentText() == 'motor fake' or self.ScanChannelBox.currentText()=='motor steps' or self.ScanChannelBox.currentText()=='motor stopover':
        #     if self.shuttle_mode != 'Trigger':
        #         self.fScanLog("Error: Motor shuttle mode should be Trigger!", red_text=True)
        #         return

        self.scan = ListenerThread(self)
        self.scan.Plot.connect(self.PlotScan)
        self.scan.Log.connect(self.fScanLog)
        self.scan.Bar.connect(self.ScanBar)
        self.scan.ScanTime.connect(self.ScanTime)
        self.scan.FailSafe.connect(self.FailSafe)
        self.scan.SetWaveLength.connect(self.SetWaveLength)
        # self.scan.Freerun_Points.connect(self.Freerun_Points)
        self.scan.MotorProgressBar.connect(self.OnMotorShuttleProgressBar)
        self.scan.start()
        self.scan_is_running = True

        # if self.ScanModeBox.currentText() == 'motor fake' or self.ScanChannelBox.currentText() == 'motor steps' or self.ScanChannelBox.currentText() == 'motor stopover':
        if self.shuttle_mode == 'Trigger':
            self.fScanLog('Attention: Motor shuttle is running (Trigger mode)!')
            self.shuttle_thread = MotorShuttleThread(self)
            self.shuttle_thread.PositionLCD.connect(self.OnMotorPositionLCD)
            self.shuttle_thread.StepperMotorLog.connect(self.StepperMotorLog)
            self.shuttle_thread.ShuttleProgressBar.connect(self.OnMotorShuttleProgressBar)
            self.shuttle_thread.start()
        return


    def OnStartTrigger(self):
        self.TriggerLog.append('Starting trigger...')
        units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, t2jumpfile = self.readboard()
        self.TriggerLog.append('Starting Pulse Blaster...')
        scantime = 'None'
        scanchannel = 'None'
        relchannel = 'None'
        bgchannel = 'None'

        pulseseq = ps.sequence(units1, chno1, time1, units2, chno2, time2, freqs,
                               pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, bgchannel,
                               scantime, t2jumpfile)
        try:
            pulseseq.seq()
        except Exception as err:
            self.TrigLog(str(err.args[0]), red_text=True)
            return False

        output = pulseseq.output
        # output,msg=ps.pulses(units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, scantime, t2jumpfile)
        status = pb.pulseblaster_program(output)
        time.sleep(0.1)
        pb.pulseblaster_start()
        # print(pb.pulseblaster_status())
        self.RunningStatusBox.setStyleSheet("QTextBrowser {background-color: rgb(0,255,0);}")
        self.TriggerLog.append('Done.')
        return

    def OnStopTrigger(self):
        if devmode == 1:
            self.TriggerLog.append('No trigger running: devmode.')
        else:
            self.TriggerLog.append('Stopping trigger...')
            try:
                pb.pulseblaster_stop()
            except RuntimeError:
                self.TriggerLog.append('Trigger already stopped')
        # print(pb.pulseblaster_status())
        self.RunningStatusBox.setStyleSheet("QTextBrowser {background-color: rgb(255,0,0);}")
        self.TriggerLog.append('Done.')
        return

    def start_working(self):
        self.TriggerLog.append('Welcome! You are using HyT-Control v.%s.\r' \
                               'Have fun!' % version)

        if devmode == 1:
            self.TriggerLog.append('Running in developer mode: No communication to devices.')

        return

    def OnScopeUpdate(self):
        scope = lc.WaveRunner64MXsB()
        scope.connect()
        channels = []

        for ch in range(1, 5):
            if eval('self.ScopeCB%i.isChecked()' % ch):
                channels.append(ch)
                print(eval('scope.ch%i_verticaloffset' % (ch)))

                exec('ch%i_verticaloffset = self.ScopeVOffs%i.value()' % (ch, ch))
                exec('ch%i_verticalcoupling = self.ScopeCP%i.currentText()' % (ch, ch))
                exec('ch%i_invert = self.ScopeInv%i.isChecked()' % (ch, ch))
                exec('ch%i_scale = self.ScopeVScale%i.value()' % (ch, ch))

                exec("self.ScopeVOffs%i_act.setText(str(ch%i_verticaloffset))" % (ch, ch))
                exec('self.ScopeCP%i_act.setText(str(ch%i_verticalcoupling))' % (ch, ch))
                exec('self.ScopeInv%i_act.setText(str(ch%i_invert))' % (ch, ch))
                exec('self.ScopeVScale%i_act.setText(str(ch%i_scale))' % (ch, ch))

                exec('scope.ch%i_verticaloffset = self.ScopeVOffs%i.value()' % (ch, ch))
                exec('scope.ch%i_verticalcoupling = self.ScopeCP%i.currentText()' % (ch, ch))
                exec('scope.ch%i_invert = self.ScopeInv%i.isChecked()' % (ch, ch))
                exec('scope.ch%i_scale = self.ScopeVScale%i.value()' % (ch, ch))

                print(eval('scope.ch%i_verticaloffset' % (ch)))
                scope.numberofsegments = self.AverageSBox.value()

        scope.initialize(channels)
        scope.disconnect()

        return

    def OnT2Jump(self):
        Temp_file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', './', 'Data files (*.dat *.txt *out)')[0]
        if Temp_file:
            t2jumpfile = Temp_file
            self.t2jumpfile = t2jumpfile
            self.TrigLog('T2Jump file set to: %s' % t2jumpfile)
            self.T2JumpLine.setText(t2jumpfile)
        else:
            self.TrigLog('Warning: no T2jump.out file has been chosen!', red_text=True)
        return

    def OnScan3dFile(self):
        Temp_file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', './', 'Data files (*.dat *.txt *out)')[0]
        if Temp_file:
            scan3dfile = Temp_file
            self.fScanLog('3d Scan file set to: %s' % scan3dfile)
            self.Scan3dText.setText(scan3dfile)
        else:
            self.fScanLog('Warning: no new file has been chosen!', red_text=True)
        return

    def OnScan3dListFile(self):
        Temp_file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', './', 'Data files (*.dat *.txt *out)')[0]
        if Temp_file:
            scan3dlistfile = Temp_file
            self.fScanLog('3d Scan List file set to: %s' % scan3dlistfile)
            self.Scan3dListText.setText(scan3dlistfile)
        else:
            self.fScanLog('Warning: no new file has been chosen!', red_text=True)
        return

    def OnListScanFile(self):
        listscanfile = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', './', 'Data files (*.dat *.txt *out)')[
            0]
        self.ScanLog.append('1d listscan file set to: %s' % listscanfile)
        self.Scan3dText.setText(listscanfile)
        return

    def OnConfig(self):
        configfile = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', './', 'Data files (*.dat *.txt)')[0]
        self.TriggerLog.append('Config file set to: %s' % configfile)
        self.ConfigLine.setText(configfile)
        return

    def OnWavelengthUpdate(self):
        self.wavelength = float(self.WavelengthSB.value())
        self.SetWaveLength(self.wavelength)
        self.fScanLog('New wavelength: %s 1/cm' % self.wavelength)
        return

    def OnWavelengthUp(self):
        self.wavelength = float(self.WavelengthSB.value()) + float(self.WavelengthIncrementBox.value())
        self.SetWaveLength(self.wavelength)
        self.WavelengthSB.setValue(self.wavelength)
        self.fScanLog('New wavelength: %s 1/cm' % self.wavelength)
        return

    def OnWavelengthDown(self):
        self.wavelength = float(self.WavelengthSB.value()) - float(self.WavelengthIncrementBox.value())
        self.SetWaveLength(self.wavelength)
        self.WavelengthSB.setValue(self.wavelength)
        self.fScanLog('New wavelength: %s 1/cm' % self.wavelength)
        return

    def SetWaveLength(self, wavelength):
        print('Setting wavelength to ', wavelength, ' cm**-1')
        if self.define_pulsare:
            self.my_pulsare.setWavelength(wavelength_cmminus1=wavelength, use_nanometers=False)
            self.my_pulsare.getActualPosition()
        else:
            print('Self.define_pulsare set to False. No connection to the dye laser, change the __init__')
        return

    def OnStopScan(self):
        try:
            self.TriggerTab.setDisabled(False)
            self.scan.f.close()
            # if self.ScanModeBox.currentText() == 'motor fake' or self.ScanChannelBox.currentText()=='motor steps' or self.ScanChannelBox.currentText()=='motor stopover':
            if self.shuttle_mode == 'Trigger':
                self.shuttle_thread.terminate()
                self.OnMotorGohome()
            self.scan.terminate()
            self.scan_is_running = False
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
        pulseseq = ps.sequence(units1, chno1, time1, units2, chno2, time2, freqs,
                               pulsenature, checkboxvals, posnegvals, scanchannel, relchannel, bgchannel,
                               scantime, t2jumpfile)
        try:
            pulseseq.seq()
        except Exception as err:
            self.TrigLog(str(err.args[0]), red_text=True)
            return
        output = pulseseq.output
        # self.TriggerLog.append('Pulse sequence ok.')
        self.TriggerLog.append('Plotting...')

        self.PlotT = TrigPlotThread(output, checkboxvals)
        self.PlotT.TrigLog.connect(self.TrigLog)
        self.PlotT.PlotSeq.connect(self.PlotSeq)
        self.PlotT.start()

        return

    def bug(self):
        plt.show()
        sys.exit()

    def OnLoadBoard(self):
        self.OnStopTrigger()
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
        # pulseseq.seq()
        try:
            pulseseq.seq()
        except Exception as err:
            self.TrigLog(str(err.args[0]), red_text=True)
            return
        output = pulseseq.output

        # if devmode == 0:
        #     ans = pb.pulseblaster_program(output)
        #     self.TriggerLog.append('PB status: %s.' % ans)
        return

    def PlotSeq(self, xvals, yvals, newlabels):
        fig, ax = plt.subplots()
        ticks = [0.25 + i for i in range(len(xvals))]
        ax.set_yticks(ticks)
        ax.set_yticklabels(newlabels)
        plt.ylabel('Channel no.')
        plt.xlabel('time [$\mu$s]')

        no = 0
        for i in range(len(xvals)):
            plt.plot(np.array(xvals[i]) * 1e6, (np.array(yvals[i]) / 2) - (newlabels[i] - 1) + i)
        plt.show()

    def readboard_class(self):
        ''' read values from channels board'''

        self.t2jumpfile = self.T2JumpLine.text()
        lastindex = 12
        indexes = [['%ib' % i, '%ie' % i] for i in range(1, lastindex)]
        # ~ print(indexes)

        units = {'ns': 1e-9, 'us': 1e-6, 'ms': 1e-3, 's': 1}

        # ~ read all channels and times from gui
        self.units1 = []
        self.units1.append(self.ch_time1_unitbox.currentText())

        # ~ time 1 units in board
        for i in range(2, lastindex + 1):
            self.units1.append(eval('self.ch_time1_unitbox_' + str(i) + '.currentText()'))

        # ~ channel number relative
        self.chno1 = []
        self.chno1.append(self.ch_combo_chno1.currentText())
        for i in range(2, lastindex + 1):
            self.chno1.append(eval('self.ch_combo_chno1_' + str(i) + '.currentText()'))

        # ~ time 1
        self.time1 = []
        self.time1.append(self.ch_time1_sbox.value() * units[self.units1[0]])
        for i in range(2, lastindex + 1):
            self.time1.append(eval('self.ch_time1_sbox_' + str(i) + '.value()*units[self.units1[' + str(i - 1) + ']]'))

        # ~ time 2 units
        self.units2 = []
        self.units2.append(self.ch_time2_unitbox.currentText())
        for i in range(2, lastindex + 1):
            self.units2.append(eval('self.ch_time2_unitbox_' + str(i) + '.currentText()'))

        # ~ second relative channel
        self.chno2 = []
        self.chno2.append(self.ch_combo_chno2.currentText())
        for i in range(2, lastindex + 1):
            self.chno2.append(eval('self.ch_combo_chno2_' + str(i) + '.currentText()'))

        # ~ time 2
        self.time2 = []
        self.time2.append(self.ch_time2_sbox.value() * units[self.units2[0]])
        for i in range(2, lastindex + 1):
            self.time2.append(eval('self.ch_time2_sbox_' + str(i) + '.value()*units[self.units2[' + str(i - 1) + ']]'))

        # ~ frequency of the pulse
        self.freqs = []
        self.freqs.append(self.ch_freqbox.value())
        for i in range(2, lastindex + 1):
            self.freqs.append(eval('self.ch_freqbox_' + str(i) + '.value()'))

        # ~ which kind of pulse? decelerator sequence? normal pulse?
        self.pulsenature = []
        self.pulsenature.append(self.ch_pulsebox.currentText())
        for i in range(2, lastindex + 1):
            self.pulsenature.append(eval('self.ch_pulsebox_' + str(i) + '.currentText()'))

        # ~ checkbox : active / inactive
        self.checkboxvals = []
        self.checkboxvals.append(self.ch_checkBox.isChecked())
        for i in range(2, lastindex + 1):
            self.checkboxvals.append(eval('self.ch_checkBox_' + str(i) + '.isChecked()'))

        self.posnegvals = []
        self.posnegvals.append(self.ch_posneg.currentText())
        for i in range(2, lastindex + 1):
            self.posnegvals.append(eval('self.ch_posneg_' + str(i) + '.currentText()'))

        self.shuttle_reps = self.MotorShuttleRepetitionsBox.value()
        self.shuttle_mode = self.MotorShuttleModeCombo.currentText()
        self.shuttle_delay_time = self.MotorShuttleDelayBox.value() * units[self.MotorShuttleDelayUnitBox.currentText()]
        self.shuttle_oneway_steps = self.MotorShuttleStepsBox.value()
        self.shuttle_stopover_time = self.MotorShuttleStopoverBox.value() * units[
            self.MotorShuttleStopoverUnitBox.currentText()]

        return self.units1, self.chno1, self.time1, self.units2, self.chno2, self.time2, self.freqs, self.pulsenature, self.checkboxvals, self.posnegvals, self.t2jumpfile

    def readboard(self):
        ''' read values from channels board'''

        t2jumpfile = self.T2JumpLine.text()
        lastindex = 12
        indexes = [['%ib' % i, '%ie' % i] for i in range(1, lastindex)]
        # ~ print(indexes)

        units = {'ns': 1e-9, 'us': 1e-6, 'ms': 1e-3, 's': 1}

        # ~ read all channels and times from gui
        units1 = []
        units1.append(self.ch_time1_unitbox.currentText())

        # ~ time 1 units in board
        for i in range(2, lastindex + 1):
            units1.append(eval('self.ch_time1_unitbox_' + str(i) + '.currentText()'))

        # ~ channel number relative
        chno1 = []
        chno1.append(self.ch_combo_chno1.currentText())
        for i in range(2, lastindex + 1):
            chno1.append(eval('self.ch_combo_chno1_' + str(i) + '.currentText()'))

        # ~ time 1
        time1 = []
        time1.append(self.ch_time1_sbox.value() * units[units1[0]])
        for i in range(2, lastindex + 1):
            time1.append(eval('self.ch_time1_sbox_' + str(i) + '.value()*units[units1[' + str(i - 1) + ']]'))

        # ~ time 2 units
        units2 = []
        units2.append(self.ch_time2_unitbox.currentText())
        for i in range(2, lastindex + 1):
            units2.append(eval('self.ch_time2_unitbox_' + str(i) + '.currentText()'))

        # ~ second relative channel
        chno2 = []
        chno2.append(self.ch_combo_chno2.currentText())
        for i in range(2, lastindex + 1):
            chno2.append(eval('self.ch_combo_chno2_' + str(i) + '.currentText()'))

        # ~ time 2
        time2 = []
        time2.append(self.ch_time2_sbox.value() * units[units2[0]])
        for i in range(2, lastindex + 1):
            time2.append(eval('self.ch_time2_sbox_' + str(i) + '.value()*units[units2[' + str(i - 1) + ']]'))

        # ~ frequency of the pulse
        freqs = []
        freqs.append(self.ch_freqbox.value())
        for i in range(2, lastindex + 1):
            freqs.append(eval('self.ch_freqbox_' + str(i) + '.value()'))

        # ~ which kind of pulse? decelerator sequence? normal pulse?
        pulsenature = []
        pulsenature.append(self.ch_pulsebox.currentText())
        for i in range(2, lastindex + 1):
            pulsenature.append(eval('self.ch_pulsebox_' + str(i) + '.currentText()'))

        # ~ checkbox : active / inactive
        checkboxvals = []
        checkboxvals.append(self.ch_checkBox.isChecked())
        for i in range(2, lastindex + 1):
            checkboxvals.append(eval('self.ch_checkBox_' + str(i) + '.isChecked()'))

        posnegvals = []
        posnegvals.append(self.ch_posneg.currentText())
        for i in range(2, lastindex + 1):
            posnegvals.append(eval('self.ch_posneg_' + str(i) + '.currentText()'))

        return units1, chno1, time1, units2, chno2, time2, freqs, pulsenature, checkboxvals, posnegvals, t2jumpfile

    def OnSaveConfig(self):
        filename_in = './input/HyT-input.dat'
        dialog = QtWidgets.QFileDialog.getSaveFileName(self, "Save file", "", "")
        # ~ dialog = QtWidgets.QFileDialog.selectFile(self,'Open file', './', 'Data files (*.dat *.txt)')
        filename_out = dialog[0]
        try:
            self.saveconfig(filename_in, filename_out)
            msg = 'Output file: %s' % filename_out
        except FileNotFoundError:
            msg = 'File not found!'
        return msg

    def saveconfig(self, filename_in, filename_out):
        infile = open(filename_in, "r")
        outfile = open(filename_out, "w")

        # ~ READ OUT ALL GUI VARIABLES ! (TO DO!)---------------------------------------------------

        # ~ triggercontrol guiparams
        trigparams = []
        chs = 19
        for i in range(1, chs):
            trigparams.append('ch%ib' % i)
            trigparams.append('ch%ie' % i)
            trigparams.append('ch%ipulse' % i)
            trigparams.append('ch%ifreq' % i)
            trigparams.append('ch%itrig' % i)
            trigparams.append('ch%icb' % i)
            for j in range(1, 3):
                trigparams.append('ch%it%i' % (i, j))
                trigparams.append('ch%iu%i' % (i, j))

        chu1, chb, cht1, chu2, che, cht2, chfreq, chpulse, chcb, chtrig, t2jumpfile = self.readboard()
        units = {'us': 1e6, 'ms': 1e3, 'ns': 1e9, 's': 1}
        for i in range(len(chu1)):
            cht1[i] = cht1[i] * units[chu1[i]]
            cht2[i] = cht2[i] * units[chu2[i]]

        # rest of the gui changeable params in dictionary saveparams
        saveparams = {}
        saveparams['T2jumpfile'] = self.T2JumpLine.text()

        while True:
            line = infile.readline()
            paramname = ''
            value = ''
            num = 0

            if line[0] != '#':
                for char in line:
                    if char == '=':
                        break
                    paramname += char
                    num += 1

                paramname = paramname.strip()
                value = eval(paramname)
                if paramname in trigparams:
                    tstr = ''.join(i for i in paramname if not i.isdigit())
                    if paramname[-1].isdigit():
                        tstr += paramname[-1]

                    tnum = [int(i) for i in paramname if i.isdigit()]
                    value = eval(tstr)[tnum[0] - 1]
                    try:
                        float(value)
                        value = "%.3f" % value
                    except (TypeError, ValueError):
                        value = value

                if paramname in saveparams:
                    value = saveparams[paramname]

                outfile.write('%s = %s\r\n' % (paramname, value))

            else:
                outfile.write(line)

            if 'this is the end' in line:
                break
        return

    def OnBrowse(self):
        dialog = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', './', 'Data files (*.dat *.txt)')
        return dialog

    def PlotScan(self, xdata, ydata, frame):
        xdata = np.array(xdata)
        ydata = np.array(ydata)

        exec("self.ScanPlt%i.setData(xdata, ydata)" % frame)
        return

    def InitMotor(self):
        self.udaq = unidaq.UniDaq()
        self.udaq.initalize()
        print(self.udaq.GetCardInfo(0))

        self.my_motor = stepper_motor.motor()
        self.MotorLog.append('Notes:\r' \
                             '- Sign of relative steps: (+) towards Stark Dec, (-) towards ion trap.\r' \
                             '- Home position: when program starts, the last position motor stopped; click Set ' \
                             'home to set it to the current position.\r' \
                             '- Test mode: shuttle without trigger; Trigger mode: shuttle triggered by CH12.\r' \
                             '- Stopover: the time the motor stays before shuttling back.\n\r\nHistory:')

        self.MotorSetRunFreqButton.clicked.connect(self.OnMotorSetFreq)
        self.MotorGoRelativeStepsButton.clicked.connect(self.OnMotorGoRelativeSteps)
        self.MotorSetHomeButton.clicked.connect(self.OnMotorSethome)
        self.MotorGoHomeButton.clicked.connect(self.OnMotorGohome)
        self.MotorDisconnectButton.clicked.connect(self.OnMotorDisconnect)
        self.MotorReconnectButton.clicked.connect(self.OnMotorReconnect)
        self.MotorStartShuttleButton.clicked.connect(self.OnMotorStartShuttle)
        self.MotorStopShuttleButton.clicked.connect(self.OnMotorStopShuttle)
        self.MotorEnableCheckbox.clicked.connect(self.OnMotorEnable)
        self.MotorShuttleModeCombo.activated.connect(self.OnMotorShuttleMode)
        self.ScanChannelBox.activated.connect(self.OnScanChannelChoose)
        # self.ScanRelChannelBox.activated.connect(self.OnScanRelativeChannelChoose)

        self.MotorEnableCheckbox.setCheckState(enable_cb_motor)
        self.MotorRunFreqBox.setValue(run_freq_motor)
        self.MotorRelativeStepsBox.setValue(go_relative_steps_motor)
        self.MotorShuttleRepetitionsBox.setValue(shuttle_repetitions_motor)
        index = self.MotorShuttleModeCombo.findText(shuttle_mode_motor)
        self.MotorShuttleModeCombo.setCurrentIndex(index)
        self.MotorShuttleDelayBox.setValue(shuttle_delay_motor)
        index = self.MotorShuttleDelayUnitBox.findText(shuttle_delay_unit_motor)
        self.MotorShuttleDelayUnitBox.setCurrentIndex(index)
        self.MotorShuttleEndCheckRadioButton.setChecked(shuttle_end_check_motor)
        self.MotorShuttleEndCheckRadioButton.clicked.connect(self.OnShuttleEndCheckSwitch)
        self.OnShuttleEndCheckSwitch()
        self.MotorShuttleStepsBox.setValue(shuttle_steps_motor)
        self.MotorShuttleStopoverBox.setValue(shuttle_stopover_motor)
        index = self.MotorShuttleStopoverUnitBox.findText(shuttle_stopover_unit_motor)
        self.MotorShuttleStopoverUnitBox.setCurrentIndex(index)
        self.OnMotorShuttleMode()
        self.MotorShuttleStepsBox.valueChanged.connect(self.OnShuttleTotalTimeUpdate)
        self.MotorShuttleStopoverBox.valueChanged.connect(self.OnShuttleTotalTimeUpdate)
        self.OnShuttleTotalTimeUpdate() ## Update the duration for Channel 12 (Motor shuttle) according to the total shuttle time for one shuttle

        if enable_cb_motor:
            try:
                self.my_motor.open_serial_port()
            except Exception as err:
                self.StepperMotorLog('- Motor connection failed:' + str(err.args))
            else:
                self.StepperMotorLog('- Motor connected')

                self.my_motor.set_run_freq(run_freq_motor)
                read_run_freq = self.my_motor.read_run_freq()
                if read_run_freq and str(run_freq_motor) in read_run_freq:
                    self.StepperMotorLog("- Run frequency set to %d Hz" % run_freq_motor)
                else:
                    self.StepperMotorLog("- Run frequency setting failed. Please check connection.")

                self.my_motor.set_acc_ramp(acc_ramp_motor)
                read_acc_ramp = self.my_motor.read_acc_ramp()
                if read_acc_ramp and str(acc_ramp_motor) in read_acc_ramp:
                    self.StepperMotorLog("- Acceleration ramp set to %d Hz/s" % acc_ramp_motor)
                else:
                    self.StepperMotorLog("- Acceleration ramp setting failed. Please check connection.")

                self.MotorStatusLabel.setStyleSheet("QLabel {background-color: rgb(0,255,0);}")
        else:
            self.frame_Motor.setVisible(False)

    def OnMotorEnable(self):
        if self.MotorEnableCheckbox.isChecked():
            self.frame_Motor.show()
            self.OnMotorReconnect()
        else:
            self.frame_Motor.hide()
            self.OnMotorDisconnect()

    def OnMotorSetFreq(self):
        if self.my_motor.isConnected:
            run_freq_motor = self.MotorRunFreqBox.value()
            self.my_motor.set_run_freq(run_freq_motor)
            read_run_freq = self.my_motor.read_run_freq()
            if read_run_freq and str(run_freq_motor) in read_run_freq:
                self.StepperMotorLog("- Run frequency set to %d Hz" % run_freq_motor)
            else:
                self.StepperMotorLog("- Run frequency setting failed. Please check connection.")
        else:
            self.StepperMotorLog('- Warning: Motor is not connected', red_text=True)
        return

    def OnMotorGoRelativeSteps(self):
        if self.my_motor.isConnected:
            rel_steps = int(self.MotorRelativeStepsBox.value())
            self.my_motor.rotate(rel_steps)
            self.StepperMotorLog("- Motor moved %d steps" % rel_steps)
            # self.MotorPositionLCD.display(self.my_motor.position)
            self.OnMotorPositionLCD(self.my_motor.position)
        else:
            self.StepperMotorLog('- Warning: Motor is not connected', red_text=True)
        return

    def OnMotorSethome(self):
        if self.my_motor.isConnected:
            self.my_motor.position = 0
            self.MotorPositionLCD.display(0)
            self.StepperMotorLog('- Set home')
        else:
            self.StepperMotorLog('- Warning: Motor is not connected', red_text=True)
        return

    def OnMotorGohome(self):
        # if self.my_motor.isConnected:
        #     self.my_motor.rotate(-self.my_motor.position)
        #     self.MotorPositionLCD.display(self.my_motor.position)
        #     self.StepperMotorLog('- Motor gone home')
        # else:
        #     self.StepperMotorLog('- Warning: Motor is not connected', red_text=True)
        # return
        if self.my_motor.isConnected:
            if self.motor_end_check:
                # motor_end_check_fail = False
                #
                # # print('just before motor command', time.perf_counter())
                # self.my_motor.rotate_pos(12000)
                # t0 = time.perf_counter()
                # # print('just after motor command', time.perf_counter())
                # while self.udaq.ReadDIBit(0, 0, 1)[1] == 1:
                #     if time.perf_counter()-t0 > 10:
                #         motor_end_check_fail = True
                #         self.my_motor.stop_rotation()
                #         break
                #     elif time.perf_counter()-t0 > 12000/self.my_motor.run_freq:
                #         self.my_motor.rotate_pos(10000)
                #     pass
                # # print('just after while cycle', time.perf_counter())
                #
                # if motor_end_check_fail:
                #     self.StepperMotorLog('- Motor failed to go home, please do it manually', red_text=True)
                # else:
                #     # print('just before motor stop', time.perf_counter())
                #     self.my_motor.stop_rotation()
                #     # print('just after motor stop', time.perf_counter())
                #     print("Motor end switch")
                #     self.my_motor.rotate_neg(350)
                #     time.sleep(350 / self.my_motor.run_freq)
                if self.OnMotorEndCheck(look_for_home=True) == False:
                    self.StepperMotorLog('- Motor failed to go home, please do it manually', red_text=True)
                else:
                    self.MotorPositionLCD.display(self.my_motor.position)
                    self.StepperMotorLog('- Motor gone home')
            else:
                self.my_motor.rotate(-self.my_motor.position)
                self.MotorPositionLCD.display(self.my_motor.position)
                self.StepperMotorLog('- Motor gone home')
        else:
            self.StepperMotorLog('- Warning: Motor is not connected', red_text=True)
        return

    def OnMotorEndCheck(self, look_for_home=False):
        if look_for_home:
            self.my_motor.rotate_pos(int(10600))
        else:
            self.my_motor.rotate_pos(int(870))

        t0 = time.perf_counter()
        while self.udaq.ReadDIBit(0, 0, 1)[1] == 1:
            if time.perf_counter() - t0 > 10:
                self.my_motor.stop_rotation()
                return False
            elif time.perf_counter() - t0 > 3:  # 000 / self.my_motor.run_freq:
                self.my_motor.rotate_pos(self.MotorShuttleEndCheckStepsBox.currentValue())
                time.sleep(2000 / self.my_motor.run_freq)
                print('steps missed', self.my_motor.run_freq)
            pass
        self.my_motor.stop_rotation()
        print("Motor end switch detected")
        self.my_motor.rotate_neg(self.MotorShuttleEndCheckStepsBox.currentValue()+50)
        # time.sleep(850 / self.my_motor.run_freq)
        return True

    def OnMotorDisconnect(self):
        if self.my_motor.isConnected:
            try:
                self.my_motor.close_serial_port()
            except Exception as err:
                self.StepperMotorLog('- Motor disconnection failed:' + str(err.args))
            else:
                self.StepperMotorLog('- Motor disconnected')
                self.MotorStatusLabel.setStyleSheet("QLabel {background-color: rgb(255,0,0);}")
        return

    def OnMotorReconnect(self):
        if not self.my_motor.isConnected:
            try:
                self.my_motor.open_serial_port()
            except Exception as err:
                self.StepperMotorLog('- Motor connection failed:' + str(err.args))
            else:
                self.StepperMotorLog('- Motor connected')
                self.MotorStatusLabel.setStyleSheet("QLabel {background-color: rgb(0,255,0);}")
        return

    def OnMotorShuttleMode(self):
        if self.MotorShuttleModeCombo.currentText() == "Test":
            self.shuttle_mode = "Test"
            self.ch_checkBox_12.setCheckState(False)
            self.ch_checkBox_12.setDisabled(True)
            self.MotorShuttleRepetitionsBox.setDisabled(False)
            self.MotorShuttleDelayBox.setDisabled(False)
            self.MotorShuttleDelayUnitBox.setDisabled(False)
            self.label_18.setDisabled(False)
            self.label_55.setDisabled(False)
            self.label_62.setVisible(True)
            self.MotorStartShuttleButton.setDisabled(False)
            self.MotorStopShuttleButton.setDisabled(False)
            self.MotorShuttleProgressBar.setVisible(True)
        else:
            self.shuttle_mode = "Trigger"
            self.ch_checkBox_12.setCheckState(True)
            self.ch_checkBox_12.setDisabled(True)
            self.MotorShuttleRepetitionsBox.setDisabled(True)
            self.MotorShuttleDelayBox.setDisabled(True)
            self.MotorShuttleDelayUnitBox.setDisabled(True)
            self.label_18.setDisabled(True)
            self.label_55.setDisabled(True)
            self.label_62.setVisible(False)
            self.MotorStartShuttleButton.setDisabled(True)
            self.MotorStopShuttleButton.setDisabled(True)
            self.MotorShuttleProgressBar.setVisible(False)

    def OnScanChannelChoose(self):
        # self.OnScanRelativeChannelChoose()
        if self.ScanChannelBox.currentText()=='motor steps':
            self.ScanRelChannelBox.setDisabled(True)
            index = self.ScanUnitBox.findText("steps")
            self.ScanUnitBox.setCurrentIndex(index)
        elif self.ScanChannelBox.currentText() == 'motor stopover':
            self.ScanRelChannelBox.setDisabled(True)
            index = self.ScanUnitBox.findText("s")
            self.ScanUnitBox.setCurrentIndex(index)
        else:
            self.ScanRelChannelBox.setDisabled(False)
            index = self.ScanUnitBox.findText("us")
            self.ScanUnitBox.setCurrentIndex(index)

    # def OnScanRelativeChannelChoose(self):
    #     if self.ScanRelChannelBox.currentText() == 'MSend' and self.ScanChannelBox.currentText() != '3b':
    #         msg = QMessageBox()
    #         msg.setWindowTitle('Warning')
    #         msg.setIcon(QMessageBox.Warning)
    #         msg.setText('Only a combination of scan channel 3b with MSend make sense! \n Scan channel has been set to 3b.')
    #         msg.exec()
    #         index = self.ScanChannelBox.findText("3b")
    #         self.ScanChannelBox.setCurrentIndex(index)

    def OnMotorStartShuttle(self):
        self.shuttle_reps = self.MotorShuttleRepetitionsBox.value()
        self.shuttle_mode = self.MotorShuttleModeCombo.currentText()
        self.shuttle_delay_time = self.MotorShuttleDelayBox.value() * units[self.MotorShuttleDelayUnitBox.currentText()]
        self.motor_end_check = self.MotorShuttleEndCheckRadioButton.isChecked()
        self.shuttle_oneway_steps = self.MotorShuttleStepsBox.value()
        self.shuttle_stopover_time = self.MotorShuttleStopoverBox.value() * units[
            self.MotorShuttleStopoverUnitBox.currentText()]
        self.MotorShuttleProgressBar.setValue(0)

        self.shuttle_thread = MotorShuttleThread(self)
        self.shuttle_thread.PositionLCD.connect(self.OnMotorPositionLCD)
        self.shuttle_thread.StepperMotorLog.connect(self.StepperMotorLog)
        self.shuttle_thread.ShuttleProgressBar.connect(self.OnMotorShuttleProgressBar)

        if self.my_motor.isConnected:
            self.StepperMotorLog('- Motor shuttle started (Test mode)')
            self.shuttle_thread.start()
        else:
            self.StepperMotorLog('- Warning: Motor is not connected', red_text=True)
        return

    def OnMotorStopShuttle(self):
        if self.my_motor.isConnected:
            self.StepperMotorLog('- Motor shuttle stopping...')
            self.shuttle_thread.terminate()
            # self.OnMotorGohome()
            self.my_motor.stop_rotation()
        else:
            self.StepperMotorLog('- Warning: Motor is not connected', red_text=True)
        return

    def OnMotorShuttleProgressBar(self, msg):
        self.MotorShuttleProgressBar.setValue(msg)
        return

    def OnMotorPositionLCD(self, position):
        self.MotorPositionLCD.display(position)
        return

    def OnShuttleTotalTimeUpdate(self):
        '''
        Update the Motor Shuttle duration in Channel 12 according to the total shuttle time for one shuttle
        '''
        self.shuttle_oneway_steps = self.MotorShuttleStepsBox.value()
        self.shuttle_stopover_time = self.MotorShuttleStopoverBox.value() * units[self.MotorShuttleStopoverUnitBox.currentText()]
        self.ch_time2_sbox_12.setValue(self.my_motor.shuttle_total_time(self.shuttle_oneway_steps,self.shuttle_stopover_time))
        index = self.ch_time2_unitbox_12.findText("s")
        self.ch_time2_unitbox_12.setCurrentIndex(index)
        return

    def OnShuttleEndCheckSwitch(self):
        if self.MotorShuttleEndCheckRadioButton.isChecked():
            self.motor_end_check = True
            self.MotorShuttleEndCheckStepsBox.setDisabled(False)
        else:
            self.motor_end_check = False
            self.MotorShuttleEndCheckStepsBox.setDisabled(True)

def main():
    app = QtWidgets.QApplication(sys.argv)
    form = GuiThread()

    form.show()
    app.exec_()


if __name__ == '__main__':
    main()
