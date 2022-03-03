# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:10:02 2019

@author: Gruppe Willitsch
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QThread, pyqtSignal, Qt
import gui
import sys
import andor
import numpy as np
from PyQt5.QtWidgets import QApplication, QWidget, QLabel
from PyQt5.QtGui import QIcon, QPixmap, QImage
import pyqtgraph as pg
from matplotlib import cm
import numpy as np
import time
from ctypes import *
from ctypes.wintypes import *
from unidaq import *
#import tmcl
import h5py
import datetime
from pytz import reference

class GuiThread(QtWidgets.QMainWindow, gui.Ui_gui):
    start_cont_gui = pyqtSignal()
    shutdown_cam_gui = pyqtSignal()
    expotime_cam_gui = pyqtSignal(object)
    stop_cam_gui = pyqtSignal()
    change_bin_gui = pyqtSignal(object)
    motor_home_gui = pyqtSignal()

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        
        # exposure time initiation
        self.exptime = 0.3
        self.ExpTimeSB.setValue(self.exptime)
        self.bindeg = 1 # binning degree

        # endcap dc offset
        self.EndcapOffsetA = 0
        self.EndcapOffsetB = 0
        self.EndcapSetButton.clicked.connect(self.OnEndcaps)
        
        #Camera buttons
        self.CamStartButton.clicked.connect(self.StartCam)
        self.CamConnectButton.clicked.connect(self.CamConnect)
        self.CamStopButton.clicked.connect(self.StopCam)
        self.CamDisconnectButton.clicked.connect(self.CamDisconnect)
        self.BinComboBox.currentIndexChanged.connect(self.ChangeBin)
#       self.MotorHomeButton.clicked.connect(self.MotorHome)
        self.ClearPlotButton.clicked.connect(self.OnClear)
        self.RefreshTrapButton.clicked.connect(self.OnRefreshTrap)
        self.ExpTimeButton.clicked.connect(self.OnExpoTimeGui)
        self.BackgroundButton.clicked.connect(self.OnBackground)
        
        # Region of interest
        self.ROIx0L.editingFinished.connect(self.setROI)
        self.ROIy0L.editingFinished.connect(self.setROI)
        self.ROIxdL.editingFinished.connect(self.setROI)
        self.ROIydL.editingFinished.connect(self.setROI)
        
        # change of background checkbox
        self.BackgroundCheckBox.stateChanged.connect(self.BackGndCam)
        
        # save checkbox
        self.SaveCB.stateChanged.connect(self.SaveImages)

        # Camera image widgets
        self.pixmap_label = QLabel()
        self.imagewidget = QLabel()
        
        self.win = pg.GraphicsLayoutWidget()
        self.plot = self.win.addPlot()
        self.img = pg.ImageItem()
#        self.img = pg.ImageView()

        self.plot.addItem(self.img)
        self.plot.setAspectLocked(lock = True, ratio = 1)
        
        self.ImageLayout.addWidget(self.win)
        
#       Generate image data for welcome image
        size = 1000
        self.data = np.random.normal(size=(size, size))
        
        #H
        self.data[20:40, 20:980] += 2.
        self.data[40:200, 490:510] += 2.
        self.data[200:220, 20:980] += 2.
        #E
        self.data[240:260, 20:980] +=2.
        self.data[260:420, 20:40] += 2.
        self.data[260:420, 490:510] += 2.
        self.data[260:420, 960:980] += 2.
        #L
        self.data[440:460, 20:980] +=2.
        self.data[460:620, 20:40] +=2.
        #L
        self.data[640:660, 20:980] +=2.
        self.data[660:820, 20:40] +=2.
        #O
        self.data[840:860, 20:980] +=2.
        self.data[860:980, 20:40] += 2.
        self.data[860:980, 960:980] += 2.
        self.data[980:1000, 20:980] +=2.
#
#
        self.data = pg.gaussianFilter( self.data, (3, 3))
        self.data += np.random.normal(size=(1000, 1000)) * 0.1
        

        self.img.setImage( self.data)
#        self.draw_setup(data)

        # Isocurve drawing unused
#        isodat = pg.gaussianFilter( self.data, (2, 2))
#        self.iso = pg.IsocurveItem(level=0.8, pen='g')
#        self.iso.setParentItem(self.img)
#        self.iso.setZValue(5)
        
        self.xmean = []
        self.ymean = []
        self.meanlength = 1000
        self.meanPlot = pg.PlotWidget()
        self.meanPlot.setObjectName('meanPlot')
        self.PlotLayout_2.addWidget(self.meanPlot)
        self.meanPlt = self.meanPlot.plot()
        

        self.plotymean = []
        self.plotmeanPlot = pg.PlotWidget()
        self.plotmeanPlot.setObjectName('meanPlot')
        self.PlotLayout_3.addWidget(self.plotmeanPlot)
        self.plotmeanPlt = self.plotmeanPlot.plot()
        

        self.plotmean = []
        self.blockmeanslist = []

        self.blockPlot = pg.PlotWidget()
        self.PlotLayout_4.addWidget(self.blockPlot)
        self.colors = [(255,0,0),(0,255,0),(0,0,255),(255,255,0),(0,255,255),(255,0,255),(255,102,102),(0,0,204),(153,153,255)]
        for i in range(9):
            exec('self.blockPlt%i = self.blockPlot.plot()'%i)
            exec('self.blockPlt%i.setData([1,2],[2+%i,3+%i], pen = %s)'%(i,i,i,self.colors[i]))
        
        # Contrast/color control
        self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img)
        self.hist.setLevels( self.data.min(),  self.data.max())
        self.hist.autoHistogramRange()
        
#         add histogram
        self.hist.vb.setMouseEnabled(y=False) # makes user interaction a little easier
        self.win.addItem(self.hist)
        self.win.resize(1000,1000)
        
        # Custom ROI for selecting an image region
        self.ROIx0 = 400
        self.ROIxd = 400
        self.ROIy0 = 200
        self.ROIyd = 200
        self.roi = pg.ROI([self.ROIx0, self.ROIxd], [self.ROIy0, self.ROIyd])
        self.roi.addScaleHandle([0.5, 1], [0.5, 0.5])
        self.roi.addScaleHandle([0, 0.5], [0.5, 0.5])
        self.plot.addItem(self.roi)
        self.roi.setZValue(10)  # make sure ROI is drawn above image
        self.roi.sigRegionChangeFinished.connect(self.updateROI)
        self.roi.getArrayRegion(self.data, self.img)
        self.updateROI()
        
#         Draggable line for setting isocurve level
#        self.isoLine = pg.InfiniteLine(angle=0, movable=True, pen='g')
#        self.hist.vb.addItem(self.isoLine)
#        self.isoLine.setValue(np.mean(data))
#        self.isoLine.setZValue(1000) # bring iso line above contrast controls
#        self.iso.setData(isodat)
#        self.isoLine.sigDragged.connect(self.updateIsocurve)
        
## Set a custom color map
#        colors = [
#            (0, 0, 0),
#            (45, 5, 61),
#            (84, 42, 55),
#            (150, 87, 60),
#            (208, 171, 141),
#            (255, 255, 255)]
#        cmap = pg.ColorMap(pos=np.linspace(0.0, 1.0, 6), color=colors)
#        self.img.setColorMap(cmap)#        
        
        
        # initialize some state vaiables
        self.started = False # cam running
        self.connected_cam = False # cam connected
        self.generatePgColormap('seismic') 
        # colormap, at the moment overwritten by self.hist = pg.HistogramLUTItem()
        

        
    def start_working(self):
        self.workT = CameraThread()
        worker=self.workT
        
        worker.draw_work.connect(self.draw)
        worker.start()
        return
        
#    isocurve unused at the moment
#    def updateIsocurve(self):
#        self.iso.setLevel(self.isoLine.value())
#        return
        
    def OnBackground(self):
        if self.connected_cam and self.started:
            self.camthread.getimgthread.imagenr = 0
        return
    
    def SaveImages(self):
        if self.connected_cam and self.started:
            self.camthread.getimgthread.count = 0
            self.camthread.getimgthread.save = self.SaveCB.isChecked()
        return
    
    def updateROI(self):
#        updates the region of interest data
        self.data_roi = self.roi.getArrayRegion(self.data, self.img)
#        self.ROISumLineEdit.setText('%.5f'%np.sum(self.data_roi))
        self.ROIMeanLineEdit.setText('%.5f'%np.mean(self.data_roi))
        self.ROIStdLineEdit.setText('%.5f'%np.std(self.data_roi))
        
        if len(self.ymean)>=10:
            self.PlotMeanLineEdit.setText('%.5f'%np.mean(self.ymean[-10:]))
        else:
            self.PlotMeanLineEdit.setText('%.5f'%np.mean(self.ymean))
        
        roipos = self.roi.pos()
        roisize = self.roi.size()
        self.ROIx0L.setText('%.5f'%roipos[0])
        self.ROIy0L.setText('%.5f'%roipos[1])
        self.ROIxdL.setText('%.5f'%roisize[0])
        self.ROIydL.setText('%.5f'%roisize[1])
        return
            
    def generatePgColormap(self, cm_name):
        '''
        https://github.com/pyqtgraph/pyqtgraph/issues/561
         justengel commented on Sep 15, 2017
        '''
        cmap = cm.get_cmap(cm_name)
        cmap._init()
        lut = (cmap._lut * 255).view(np.ndarray)  # Convert matplotlib colormap from 0-1 to 0 -255 for Qt
        # Apply the colormap
        self.img.setLookupTable(lut)
#        self.img.setColorMap(cmap)
        return
    
    def OnRefreshTrap(self):
        refreshtime = 10.0
        self.refresh = TrapRefreshThread(self,refreshtime)
        self.refresh.start()
        return
    
    def OnEndcaps(self):
#   #daq card
        boardno = 0
        channelA = 1
        channelB = 2
        self.EndcapOffsetA = self.EndcapOffsetAsb.value()
        self.EndcapOffsetB = self.EndcapOffsetBsb.value()
        
        print('setting endcap offset: A %.2f \t B %.2f'%(self.EndcapOffsetA,self.EndcapOffsetB))

        udaq = UniDaq()
        udaq.initalize()
        udaq.GetCardInfo(boardno)
        udaq.ConfigAO(boardno,channelA,3)
        udaq.WriteAOVoltage(boardno,channelA,self.EndcapOffsetA)
        time.sleep(1e-3)
        udaq.ConfigAO(boardno,channelB,3)
        udaq.WriteAOVoltage(boardno,channelB,self.EndcapOffsetB)
        time.sleep(1e-3)
        udaq.close()
        return
    
    def CamConnect(self):
        if not self.connected_cam:
            print('connecting to camera...')
            self.camthread = CameraThread(self)
            self.camthread.draw_work.connect(self.draw)
            self.camthread.draw_setup_work.connect(self.draw_setup)
            self.change_bin_gui.connect(self.camthread.ChangeBin)
            self.start_cont_gui.connect(self.camthread.start_cont)
            self.shutdown_cam_gui.connect(self.camthread.OnStop_work)
            self.stop_cam_gui.connect(self.camthread.OnStop_collect)
            self.expotime_cam_gui.connect(self.camthread.CamExpTime)
            self.camthread.start()
            self.connected_cam = True
        else:
            print('already connected')
        return
    
    def ChangeBin(self):
        bins = int(self.BinComboBox.currentText())
        self.bindeg = bins
        roipos = self.roi.pos()/self.bindeg
        roisize = self.roi.size()/self.bindeg
        self.roi.setPos(roipos)
        self.roi.setSize(roisize)
        self.change_bin_gui.emit(bins)
        return
    
    def OnExpoTimeGui(self):
        self.StopCam()
        self.exptime = self.ExpTimeSB.value()
        print('new cam exp time', self.exptime)
        if self.connected_cam:
            self.expotime_cam_gui.emit(self.exptime)
        else:
            print('cam not started')
            
        if self.started:
            self.StartCam()
        return
    
    def setROI(self):
        try:
            x0 = float(self.ROIx0L.text())
            y0 = float(self.ROIy0L.text())
            xd = float(self.ROIxdL.text())
            yd = float(self.ROIydL.text())
            self.roi.setPos((x0, y0))
            self.roi.setSize((xd, yd))
        except:
            print('non float value in ROI')
            

#        self.updateROI()
        return
    
    def BackGndCam(self):
        if self.started:
            self.camthread.getimgthread.bgmode = self.BackgroundCheckBox.isChecked()
            print('background set to: %s'%(self.BackgroundCheckBox.isChecked()))
        else:
            print('camera continuous thread not running')
            pass
        return
    
    def StartCam(self):
        if not self.started:
            self.start_cont_gui.emit()
            self.started = True
        else:
            print('already started')
        return
    
    def StopCam(self):
        if self.started:
            self.stop_cam_gui.emit()
            self.started = False
            self.camthread.getimgthread.imagenr = 0
            self.camthread.getimgthread.accimagenr = 0
            self.camthread.getimgthread.accimage = np.zeros([1000,1000])
            self.camthread.getimgthread.bgimage = np.zeros([1000,1000])
            self.OnClear()
        else:
            print('not started yet')
        return
    
    def CamDisconnect(self):
        if self.connected_cam:
            while self.camthread.camera_running:
                time.sleep(1e-3)
                
            self.shutdown_cam_gui.emit()
            self.connected_cam = False
            print('cam disconnected')
        else:
            print('cam not connected')
#        udaq.close()
            
        return
    
    def MotorHome(self):
        print('motorhome gui')
        print('leckt')
        #self.MotorConnect()
        #self.motor_home_gui.emit()
        return
    
    def MotorConnect(self):
        print('motorconnect gui')
        print('Deine')
        #self.motor = MotorThread(self)
        #self.motor_home_gui.connect(self.motor.goHome)
        #self.motor.MotorDisconnect_work.connect(self.MotorDisconnect)
        
        #self.motor.start()
        time.sleep(1)
        return
    
    def MotorDisconnect(self):
        print('motordisconnect gui')
        #self.motor.terminate()
        print("Mudda")
        return
    
    def OnClear(self):
        self.xmean = []
        self.ymean = []
        self.plotymean = []
        self.blockmeanslist = []
        return
        
    def draw(self, data, bg, blockmeans):
        self.data = data-bg
        self.img.setImage(self.data, autoLevels = False, autoHistogramRange = False)
        self.plot.setXRange(0, np.shape(data)[0])
        self.plot.setYRange(0, np.shape(data)[1])
        self.updateROI()
        

        ymeanwindow = 10
        self.ymean.append(np.mean(self.data_roi))
        if len(self.ymean)>=ymeanwindow:
            self.plotymean.append(np.mean(self.ymean[-ymeanwindow:]))
        else:
            self.plotymean.append(np.mean(self.ymean))
        self.xmean = np.arange(0 , len(self.ymean), 1)
        if len(self.ymean) >= self.meanlength:
            self.ymean = self.ymean[1:]
            self.xmean = self.xmean[1:]
            self.plotymean = self.plotymean[1:]
        self.meanPlt.setData(self.xmean, np.array(self.ymean))
        self.plotmeanPlt.setData(self.xmean, np.array(self.plotymean))
#        self.iso.setData(pg.gaussianFilter(data, (2, 2)))
#        self.isoLine.setValue(np.mean(data))
        
#        self.blockmeanslist.append(blockmeans)
#        blockmeanslistnp = np.array(self.blockmeanslist).T
#        
#        for i in range(9):
#            exec('self.blockPlt%i.setData(self.xmean, blockmeanslistnp[%i], pen = %s)'%(i,i,self.colors[i]))
#        
        
        return
    
    def draw_setup(self, data):
        self.hist.setLevels(np.min(data), np.percentile(data, 99))
        return
        
class GetImageCont(QThread):
    def __init__(self, parent):
        QThread.__init__(self)
        self.parent = parent
        self.running = True
        self.imagenr = 0
        self.accimagenr = 0
        self.accimage = np.zeros([1000,1000])
        self.bgimage = np.zeros([1000,1000])
        self.bgmode = self.parent.parent.BackgroundCheckBox.isChecked()
        self.bgavgs = 10
        self.images = 200
        self.save = False
        self.count = 0
        self.filecount = 0

        
    def __del__(self):
        self.wait()
        
    def shop_bin(self,data, bindeg):
        bins = int(1000/bindeg)
        print(bins)
        x_mesh,y_mesh = np.meshgrid(np.linspace(0,999,1000),np.linspace(0,999,1000))
        x_mesh = x_mesh.ravel()
        y_mesh = y_mesh.ravel()
        data = data[:1000,:1000]
        data_hist,x,y = np.histogram2d(x_mesh,y_mesh,weights=data.ravel(),bins=(bins,bins))
        data_hist = data_hist.T
        return data_hist/(bindeg*bindeg)
        
    def run(self):
        blocks = 3
        thresh = 1000
        bmeansbg = np.zeros([blocks*blocks])
        filename = './images%i.h5'%(self.filecount)
        hf = h5py.File(filename, 'w')
        hf.close()

        while self.running:
            hf = h5py.File(filename, 'a')
            self.parent.cam.GetImageData()
            datasig = self.parent.cam.imagedata
            data = self.parent.shapeData(datasig)
            data = self.shop_bin(data,self.parent.bindeg)
            xlen, ylen = np.array(data).shape

            blockmeans = []
            bmeanss = np.zeros([blocks*blocks])
            for i in range (blocks):
                for j in range (blocks):
    #                mean = np.mean(self.data[int(np.floor(i*xlen / blocks)): int(np.floor((i+1)*xlen / blocks)), int(np.floor(j*ylen / blocks)):int(np.floor((j+1)*ylen / blocks))])
                    meanb = np.mean(data[int(np.floor(i*xlen / blocks)): int(np.floor((i+1)*xlen / blocks)), int(np.floor(j*ylen / blocks)):int(np.floor((j+1)*ylen / blocks))] > thresh)
                    val = meanb
                    blockmeans.append(val)
            bmeanss += np.array(blockmeans).reshape([blocks*blocks])
#            self.accimage += data
#            self.accimagenr += 1
#            
#            if self.accimagenr == 20:
#                np.save('name%i'%self.imagenr, self.accimage)
#                self.accimage = np.zeros_like(self.accimage)
#                self.imagenr += 1
#                self.accimagenr = 0
#                print('saved data', self.imagenr-1)
            if self.bgmode:
                if self.imagenr == 0:
                    print('taking background')
                    print('refreshing trap')
            #   #daq card
                    boardno = 0
                    channel = 0
                    voltage = 4.0
                    waittime = 4.0
                    
                    print('ramping down rf voltage')
                    udaq = UniDaq()
                    udaq.initalize()
                    udaq.GetCardInfo(boardno)
                    udaq.ConfigAO(boardno,channel,3)
                    udaq.WriteAOVoltage(boardno,channel,0.0)
                    
                    time.sleep(waittime)
            
                    
                    
                    self.bgimage = np.zeros_like(data)
                    bmeansbg = np.zeros([blocks*blocks])
                    for i in range(self.bgavgs):
                        self.parent.cam.GetImageData()
                        databg = self.parent.cam.imagedata
                        databg = self.parent.shapeData(databg)
                        databg = self.shop_bin(databg, self.parent.bindeg)
                        self.bgimage += databg
                        
                        
                        blockmeans = []
                        for i in range (blocks):
                            for j in range (blocks):
                #                mean = np.mean(self.data[int(np.floor(i*xlen / blocks)): int(np.floor((i+1)*xlen / blocks)), int(np.floor(j*ylen / blocks)):int(np.floor((j+1)*ylen / blocks))])
                                meanb = np.mean(databg[int(np.floor(i*xlen / blocks)): int(np.floor((i+1)*xlen / blocks)), int(np.floor(j*ylen / blocks)):int(np.floor((j+1)*ylen / blocks))] > thresh)
                                val = meanb
                                blockmeans.append(val)
                        bmeansbg += np.array(blockmeans).reshape([blocks*blocks])
                    bmeansbg /= self.bgavgs
                        
                        
                    self.bgimage = self.bgimage / self.bgavgs
                    print('ramping up rf voltage')
                    print(udaq.WriteAOVoltage(boardno,channel,voltage))
                    udaq.close()
                    time.sleep(waittime)

                self.imagenr += 1
#                if self.imagenr == self.images:
#                    self.imagenr = 0
                if data.shape == self.bgimage.shape:
                    self.parent.draw_work.emit(data, self.bgimage, bmeanss - bmeansbg)
            else: # no bgmode
                self.parent.draw_work.emit(data, np.zeros_like(data), bmeanss - bmeansbg)
                if self.save:
                    # new claudio changes here XXX
                    localtime = reference.LocalTimezone()
                    now = datetime.datetime.now()
                    ts = str(now.strftime("%Y-%m-%d %H:%M:%S" + localtime.tzname(now))).strip('W. Europe Daylight Time')  
                    try:
                        dset = hf.create_dataset('data%i'%self.count, data=data,compression = 'gzip')
                    except RuntimeError:
                        hf.close
                        self.filecount +=1
                        filename = './images%i.h5'%(self.filecount)
                        hf = h5py.File(filename, 'w')
                        dset = hf.create_dataset('data%i'%self.count, data=data,compression = 'gzip')
                    # new claudio changes here XXX
                    dset.attrs['timestamp'] = np.string_(ts)
                    dset.attrs['gain'] = self.parent.gain
                    dset.attrs['exposuretime'] = self.parent.exptime
                    hf.close()
                    self.count+=1

                # end if bgmode
            time.sleep(1e-1)
        print('getimg stopped')
        return
    
class CameraThread(QThread):
    draw_work = pyqtSignal(object, object, object)
    draw_setup_work = pyqtSignal(object)
    
    def __init__(self, parent):
        QThread.__init__(self)
        self.exptime = parent.exptime
        self.gain = int(10*4.066)#4066# 0 to 4066  (4.066 corresponds to 1 in andor software)
        self.state = 1
        self.trigmode = 0
        self.emccdgainmode = 1
        self.temperature=-20
        self.hbin = 1
        self.vbin = 1
        self.hstart = 1
        self.hend = 1
        self.vstart = 1
        self.vend = 1
        self.parent = parent
        self.camera_running = False
        self.binning = True
        self.bindeg = self.parent.bindeg

        
    def __del__(self):
        self.wait()

    def run(self):
        self.andor_initialize()
        return
    
    def start_cont(self):
        self.getimgthread = GetImageCont(self)
        self.getimgthread.start()
        self.camera_running = True

        print('starting')
        return
    
    def shop_bin(self,data, bindeg):
        bins = 1000/bindeg
        x_mesh,y_mesh = np.meshgrid(np.linspace(0,999,1000),np.linspace(0,999,1000))
        x_mesh = x_mesh.ravel()
        y_mesh = y_mesh.ravel()
        data = data[:1000,:1000]
        data_hist,x,y = np.histogram2d(x_mesh,y_mesh,weights=data.ravel(),bins=(bins,bins))
        data_hist = data_hist.T
        return data_hist
        
    def workerfunc_work(self,msg):
        self.Log_work.emit(msg)
        return
    
    def ChangeBin(self, bins):
        self.bindeg = bins
        return
    
    def shapeData(self, data):
        return np.array(data).reshape(int(self.cam.height/self.hbin),int(self.cam.width/self.vbin)).T
    
    def OnStop_work(self):
#        self.getimgthread.running = False
        self.cam.ShutDown()
        return
    
    def OnStop_collect(self):
        self.getimgthread.running = False
        print('stopping')
        while self.getimgthread.isRunning():
            time.sleep(1e-3)
        print('terminat0r')
        self.getimgthread.terminate()
        self.camera_running = False
        return
    
    def CamExpTime(self, val):
        self.exptime = val
        self.cam.SetExposureTime(self.exptime)
        return

    def andor_initialize(self):
        self.cam = andor.Andor()
        print('cam init', self.cam.Initialize())
        print(self.cam.GetDetector())
        print(self.cam.GetHeadModel())
        print(self.cam.HeadModel)
        self.cam.SetAcquisitionMode(1)
        self.cam.SetReadMode(4)
        self.cam.GetFastestRecommendedVSSpeed()
        self.cam.SetVSSpeed(self.cam.VSNumber)
        
        STemp = 0
        HSnumber = 0
        ADnumber = 0
        
        errorValue = self.cam.GetNumberADChannels()
        nAD = self.cam.noADChannels
        
        for iAD in range(nAD):
            index = c_int(0)
            self.cam.dll.GetNumberHSSpeeds(c_int(iAD), c_int(0), byref(index))
            
            for iSpeed in range(index.value):
                speed = c_float(0)
                self.cam.dll.GetHSSpeed(c_int(iAD), c_int(0), c_int(iSpeed), byref(speed))
                
                if speed.value > STemp:
                    STemp = speed.value
                    HSnumber = iSpeed
                    ADnumber = iAD
        
    #    print(STemp, HSnumber, ADnumber)
        self.cam.SetADChannel(ADnumber)
        self.cam.SetHSSpeed(0,HSnumber)
        
        self.cam.GetStatus()
        self.cam.GetTemperatureRange()
        print('temprange:', self.cam.tempmin, self.cam.tempmax)
        msg=self.cam.GetTemperature()
        self.cam.GetTemperatureRange()

        self.cam.StabilizeTemperature(self.temperature) #selfmade temperature set routine
        self.cam.GetTemperature()
        print('current temperature in degC:', self.cam.temperature)
    #    print(cam.SetShutter(typ, mode, closingtime, openingtime))
        
        self.hbin=1
        self.vbin=1
        self.hstart=1
        self.hend=self.cam.width
        self.vstart=1
        self.vend=self.cam.height
        self.cam.SetTriggerMode(0)
        
        self.cam.SetExposureTime(self.exptime)
        self.cam.SetTriggerMode(self.trigmode)
        self.cam.GetEMGainRange()
        print('gainrange', self.cam.gainRange)
        print('emccd mode', self.cam.SetEMAdvanced(self.state))
        print('emgainmode', self.cam.SetEMCCDGainMode(self.emccdgainmode))
        print('egain get', self.cam.GetEMCCDGain())
        self.cam.GetEMGainRange()
        print('gainrange', self.cam.gainRange)
        print('gainval', self.cam.SetEMCCDGain(self.gain))
        print(self.cam.gain)
    
        print(self.cam.GetAcquisitionTimings())
        print('Acq timings: exposure time, accumulation time, kinetic:\n %f \n %f \n %f'
          %(self.cam.exposure,self.cam.accumulate,self.cam.kinetic))
        

        print(self.cam.SetImage(self.hbin,self.vbin,self.hstart,self.hend,self.vstart,self.vend))
        print('signal')
        
#        self.cam.GetImageData()
#        datasig = self.cam.imagedata
#        data = self.shapeData(datasig)
#        data = self.shop_bin(data, self.bindeg)
#        self.draw_work.emit(data)
#        self.draw_setup_work.emit(data)
        return
    
#class MotorThread(QThread):    
#    MotorDisconnect_work = pyqtSignal()
#    def __init__(self, parent):
#        QThread.__init__(self)
##        # tmcl is the motor controller
##        self.tmcl = tmcl.TMCL()
##        self.tmcl.baudrate = 9600
##        self.tmcl.timeout = 1
##        self.tmcl.COM = 'COM28'
##        self.tmcl.connected = False
##        self.tmcl.position = 0
##        self.tmcl.MODULE_ADDRESS = 1
##        self.tmcl.current_max = 120
##        self.tmcl.current_sb = 100
##        self.tmcl.vel_max = 2047
##        self.tmcl.acc_max = 1000
##        self.tmcl.pulse_divisor = 0
##        self.tmcl.ramp_divisor = 7
##        self.tmcl.initialized = False
#        
#        self.vel_home = 100  # velocity for finding home position
#
#        self.udaq = unidaq.UniDaq()
#        
#    def connect(self):
#        print('motor connected')
#        self.tmcl.connect()
#        self.tmcl.initialize()
#        
#        self.udaq.initalize()
#        return
#    
#    def disconnect(self):
#        print('motor connected')
#
#        self.tmcl.disconnect()
#        self.udaq.close()
#        return
#    
#    def run(self):
#        print('run motorthread')
#
#        self.connect()
#        return
#    
#    def goHome(self):
#        print('goHome motorthread')
#        state = self.udaq.ReadDI(0, 0)[1]
#        if state:
#            print('truuu',state)
#        else:
#            print('false',state)
#        self.tmcl.motor.rotate_right(self.vel_home)
#        
#        time.sleep(1)
##        
#        while self.udaq.ReadDI(0, 0)[1]:
#            time.sleep(1e-2)
#            
#        print('reached home')
##        
#        self.tmcl.motor.stop()
#        self.disconnect()
#        self.MotorDisconnect_work.emit()
#        return
        
class TrapRefreshThread(QThread):    
    def __init__(self, parent, time):
        QThread.__init__(self)
        self.timeoff = time
        
    def run(self):
        print('refreshing trap')
        #   #daq card
        boardno = 0
        channel = 0
        voltage = 4.0
        waittime = self.timeoff
        
        print('ramping down rf voltage')
        udaq = UniDaq()
        udaq.initalize()
        print(udaq.GetCardInfo(boardno))
        print(udaq.ConfigAO(boardno,channel,3))
        print(udaq.WriteAOVoltage(boardno,channel,0.0))
#        udaq.close()
        
        time.sleep(waittime)

        print('ramping up rf voltage')
#        udaq = UniDaq()
#        udaq.initalize()
#        print(udaq.GetCardInfo(boardno))
#        print(udaq.ConfigAO(boardno,channel,3))
        print(udaq.WriteAOVoltage(boardno,channel,voltage))
        udaq.close()
        return
    
    

def main():
    app = QtWidgets.QApplication(sys.argv)
    form = GuiThread()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()

    
