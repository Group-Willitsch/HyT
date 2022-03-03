#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QThread, pyqtSignal
from PyQt5.QtGui import QIcon, QPixmap, QMovie
from PfeifferVacuum import MaxiGauge, MaxiGaugeError
import serial
import gui
import sys
import glob
import time
import numpy as np
import pyqtgraph as pg
import datetime
from pytz import reference
#from Users.Gruppe Willitsch.Desktop.Lab315.Stark.Software.lab_watch import send_email
sys.path.insert(1, 'C:/Users/Gruppe Willitsch/Desktop/Lab315/Lab315/Software/lab_watch')
import lab_watch


devmode = 0
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'b')

def serial_ports():
    """ Lists serial port names

        :raises EnvironmentError:
            On unsupported or unknown platforms
        :returns:
            A list of the serial ports available on the system
    """
    if sys.platform.startswith('win'):
        ports = ['COM%s' % (i + 1) for i in range(256)]
    elif sys.platform.startswith('linux') or sys.platform.startswith('cygwin'):
        # this excludes your current terminal "/dev/tty"
        ports = glob.glob('/dev/tty[A-Za-z]*')
    elif sys.platform.startswith('darwin'):
        ports = glob.glob('/dev/tty.*')
    else:
        raise EnvironmentError('Unsupported platform')

    result = []
    for port in ports:
        try:
            s = serial.Serial(port)
            s.close()
            result.append(port)
        except (OSError, serial.SerialException):
            pass
    return result

class GuiThread(QtWidgets.QMainWindow, gui.Ui_my_pressures):
#    workerfunc_gui = pyqtSignal(object)
    connectPfeiffer_gui = pyqtSignal(object)
    clearGraphs_gui = pyqtSignal()
    PlotMode_gui = pyqtSignal(object)
  
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)
        self.ports = None
        self.COM = None
        self.COM2 = None

        
        self.setWindowIcon(QtGui.QIcon('my_pressures.jpg'))

        
        self.PlotMode = 'window'
        index = self.PlotBox.findText(self.PlotMode)
        self.PlotBox.setCurrentIndex(index)
        
        self.PlotBox.currentIndexChanged.connect(self.PlotBoxFunc)
        self.COMButton.clicked.connect(self.OnCOM)
        self.ConCOMButton.clicked.connect(self.OnConCOM)
        self.ClearButton.clicked.connect(self.ClearGraphs)
        self.LogScaleCB.stateChanged.connect(self.LogScaleFunc)
        
        self.LogMode =  bool(self.LogScaleCB.checkState()/2)
        #~ self.PlotTimeWindow =
        setgif = True
        if setgif:
            label = QtWidgets.QLabel()
            pixmap = QPixmap('my_pressures.gif').scaled(self.pic_frame.width(), self.pic_frame.height())
            label.setPixmap(pixmap)
            movie = QMovie('my_pressures.gif')
            movie.setScaledSize(self.pic_frame.size()*0.9)
            label.setMovie(movie)
            movie.start()
        else:
            label = QtWidgets.QLabel()
            pixmap = QPixmap('my_pressures.jpg').scaled(self.pic_frame.width(), self.pic_frame.height())
            label.setPixmap(pixmap)
        
        self.pic_frame_layout.addWidget(label)
        
        for channel in range(1,8):
            exec('self.plot%i = pg.PlotWidget()'%channel)
            exec('self.plot%i.setObjectName("plot%i")'%(channel,channel))
            exec('self.PlotLayout%i.addWidget(self.plot%i)'%(channel,channel))
            exec('self.plt%i = self.plot%i.plot()'%(channel,channel))
            if channel<7:
                exec("self.plot%i.setLabel('left', text='p', units='mbar')"%(channel))
            else:
                exec("self.plot%i.setLabel('left', text='t', units='C')"%(channel))

        
    def start_working(self):
        self.workT = WorkThread(self.COM,self.COM2)
        worker=self.workT

        self.connectPfeiffer_gui.connect(worker.connectPfeiffer)
        self.clearGraphs_gui.connect(worker.ClearGraphs)
        self.PlotMode_gui.connect(worker.PlotModeFunc)
        worker.pressure_TB_work.connect(self.pressure_TB)
        worker.replot_work.connect(self.replot)
#        self.workerfunc_gui.connect(worker.workerfunc_work)
#        worker.Log_work.connect(self.Log)
        worker.start()
        return
        
    def OnCOM(self):
        self.COMbox.clear()
        self.ports = serial_ports()
        i=0
        for port in self.ports:
            self.COMbox.addItem("")
            self.COMbox.setItemText(i, port)
            self.COMbox2.addItem("")
            self.COMbox2.setItemText(i, port)
            i+=1
        return
    
    def OnConCOM(self):
        time.sleep(0.1)
        self.COM = self.COMbox.currentText()
        #~ self.COM2 = self.COMbox2.currentText()
        self.start_working()
        return
    
    def LogScaleFunc(self):
        self.LogMode = bool(self.LogScaleCB.checkState()/2)
        try:
            self.workT.LogMode = self.LogMode
        except:
            return
        return
    
    def pressure_TB(self,pressures,temperature):
        time.sleep(0.1)
        for chno in range(1,7):
            exec("self.Ch%iTB.setText(str(pressures[%i]))"%(chno,chno-1))
        chno = 7
        exec("self.Ch%iTB.setText(str(temperature))"%(chno))
        return
    
    
    
    def replot(self,times,timestamps,data1,data2,data3,data4,data5,data6,temperature):
        xticks = 3
        
        indexes = np.array([False]*len(times))
        trues = [int(0+i*len(times)/xticks) for i in range(xticks)]
        for index in range(indexes.shape[0]):
            if index in trues:
                indexes[index] == True
        #~ print(np.unique(trues))            
        
        
        for channel in range (1,7):
            times = np.array(times)
            exec('self.plt%i.setLogMode(False,%s)'%(channel,self.LogMode))
#            exec('data%i = np.array(data%i)'%(channel,channel))
#            exec("self.plot%i.getAxis('bottom').setTickSpacing(0.9*np.min(data%i),1.1*np.max(data%i),4)"%(channel,channel,channel))
#            exec("self.plot%i.getAxis('bottom').setTickSpacing(5, 3)"%(channel))
            exec('self.plt%i.setData(times, data%i)'%(channel,channel))
            time.sleep(1e-3)
        
        channel = 7 
        exec('self.plt%i.setData(times, temperature)'%(channel))
        return
        
    def ClearGraphs(self):
        self.clearGraphs_gui.emit()
        return
            
    
    def PlotBoxFunc(self):
        time.sleep(1e-3)
        try:
            self.PlotMode_gui.emit(self.PlotBox.currentText())
            return 1
        except:
            return 0
#    def OnClick(self):
#        self.somefunction()
#        return    
    
#    def somefunction(self):
#        self.workerfunc_gui.emit('lala')
#        return
#    def Log(self,msg):
#        self.textbox.append(msg)
#        return
        
        


class WorkThread(QThread,object):
#    Log_work = pyqtSignal(object)
    pressure_TB_work = pyqtSignal(object,object)
    replot_work = pyqtSignal(object,object,object,
                             object,object,object,
                             object,object,object)
    
    def __init__(self,COM,COM2):
        QThread.__init__(self)
        self.connected = False
        self.COM = COM
        print('com1 = %s'%COM)
        #~ print('com2 = %s'%COM2)
        #~ self.COM2 = COM2
        self.OnClear = False
        self.PlotMode = 'window'
        self.PlotWindow = 500
        self.LogMode = False

    def __del__(self):
        self.wait()

    def run(self):
        
        if devmode ==0:
            self.connectPfeiffer(self.COM)
            #~ self.connectTemperature(self.COM2)
        
        self.lab_watch  = False
        self.flag_email = False
        self.threshold  = 1e-5
        self.time_start = time.time()
        
        
        self.times = []
        self.timestamps = []
        
        lasttime = 0
        
        for channel in range(1,7):
            exec('self.data%i = []'%channel)
        self.temperatures = []
        
        # check if file for this month already exists
        now = datetime.datetime.now()-datetime.datetime(1970,1,1)
        thetime = now.total_seconds()
        day_month = time.strftime('%Y-%m-%d', time.localtime(thetime))
        filename = 'pressure_'+day_month+'.txt'
        LogDir = 'C:/Users/Gruppe Willitsch/Desktop/data/pressure_log/'
        filename = LogDir + filename
        try: 
            pressfile = open(filename, 'r')
            print('file already exists dawg')
        except FileNotFoundError:
            pressfile = open(filename, 'w+')
            pressfile.write('Time YYYY-mm-dd HH:MM:SS \t p(Sour) \t p(Dec) \t p(Trap) \t Foreline(Dec) \t Foreline(Trap) \t Nothing \t Temperature(degC)\n')
            print('new Log file: %s'%filename)
        pressfile.close()
        
            
        i = 0
        while True:
#            QtCore.QCoreApplication.processEvents()
            self.pressures = np.zeros(6)
            if self.OnClear:
                self.times = []
                self.temperatures = []
                for channel in range(1,7):
                    exec('self.data%i = []'%channel)
                self.OnClear = False


            #~ temperature = self.arduino.readline()[:-3]
            try:
                temperature = float(temperature) 
                # first read out is an empty string... replacing it by 25.0 degC !! maybe NaN would be better?
            except:
                temperature = 25.0
                
            self.temperatures.append(temperature)

            for port in range(1,7):
                if devmode == 0:
                    pressure = self.pressure(port)
                if devmode ==1:
                    pressure = 10.0
					
                self.pressures[port-1] = pressure
                exec('self.data%i.append(pressure)'%port)
                time.sleep(0.1)
                
            if self.PlotMode == 'window':
                if len(self.times)>=self.PlotWindow:
                    self.times.pop(0)
                    self.temperatures.pop(0)
                    for channel in range(1,7):
                        exec('self.data%i.pop(0)'%channel)
                    
            localtime = reference.LocalTimezone()
            now = datetime.datetime.now()
            ts = str(now.strftime("%Y-%m-%d %H:%M:%S" + localtime.tzname(now))).strip('W. Europe Daylight Time')        
            #~ now = datetime.datetime.now()-datetime.datetime(1970,1,1)
            thetime = (now-datetime.datetime(1970,1,1)).total_seconds()
            
            self.timestamps.append(now)
            timediff = 120   # time difference between print to file in s
            

            
            if thetime - lasttime >=timediff:
                if time.strftime('%Y-%m-%d', time.localtime(thetime))!=day_month:
                    day_month = time.strftime('%Y-%m-%d', time.localtime(thetime))
                    filename = LogDir + 'pressure_'+day_month+'.txt'
                    pressfile = open(filename, 'w+')
                    pressfile.write('Time YYYY-mm-dd HH:MM:SS \t p(Sour) \t p(Dec) \t p(Trap) \t Foreline(Dec) \t Foreline(Trap) \t Nothing \t Temperature(degC)\n')
                    print('new Log file: %s'%filename)
                    pressfile.close()

                pressfile = open(filename, 'a')
                lasttime = thetime
                #~ time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(thetime))
                string = ts
                
                for channel in range(0,6):
                    string += '\t'+str(self.pressures[channel])
                   
                string += '\t'+str(temperature)
                string += '\n'
                pressfile.write(string)
                pressfile.close()
                
                
                
#            timediff = time.strftime('%Y-%m-%d', time.localtime(timediff))
#            print(now)

            self.times.append(i)
                
            lent = len(self.times)
            lents = len(self.timestamps)
            len1 = len(self.data1)
            len2 = len(self.data2)
            len3 = len(self.data3)
            len4 = len(self.data4)
            len5 = len(self.data5)
            len6 = len(self.data6)
            lentemp = len(self.temperatures)
            
            lengths = [lent, lents, len1, len2, len3, len4, len5, len6, lentemp]
            
            if np.min(lengths) != np.max(lengths):
                if lengths[0] > np.min(lengths):
                    self.times = self.times[:-1]
                if lengths[1] > np.min(lengths):
                    self.timestamps = self.timestamps[:-1]
                if lengths[2] > np.min(lengths):
                    self.data1 = self.data1[:-1]
                if lengths[3] > np.min(lengths):
                    self.data2 = self.data2[:-1]
                if lengths[4] > np.min(lengths):
                    self.data3 = self.data3[:-1]
                if lengths[5] > np.min(lengths):
                    self.data4 = self.data4[:-1]
                if lengths[6] > np.min(lengths):
                    self.data5 = self.data5[:-1]
                if lengths[7] > np.min(lengths):
                    self.data6 = self.data6[:-1]
                if lengths[8] > np.min(lengths):
                    self.temperatures = self.temperatures[:-1]
                    
            self.pressure_TB_work.emit(self.pressures, temperature) 
            
            if self.LogMode:
                self.replot_work.emit(np.array(self.times),np.array(self.timestamps),np.array(self.data1),np.array(self.data2),np.array(self.data3),
                                      np.array(self.data4),np.array(self.data5),np.array(self.data6),np.array(self.temperatures))
            else:
                self.replot_work.emit(np.array(self.times),np.array(self.timestamps),np.array(self.data1),np.array(self.data2),np.array(self.data3),
                                      np.array(self.data4),np.array(self.data5),np.array(self.data6),np.array(self.temperatures)) 
            
            i+=1
            
            
#            ####Add the lab-watch here
#            if lab_watch:
#                t_act = time.time()
#                ###--- Triggered by threshold
#                if self.data3[-1]>self.threshold:
#                    if not self.flag_email:
#                        lab_watch.send_email('The pressure is below %i'%self.threshold)
#                        self.flag_email = True
#                        t_send = time.time()
#                        print('Email was send to Mr. x')
#                    else:
#                        if(t_act-t_send)>43200:
#                            self.flag_email = False
#                ###--- Daily check if everthing is fine
#                if (t_act-self.time_start)>86400:
#                    lab_watch.send_email('Daily mail: pressure is %s'%self.data3[-1])
#                    self.time_start = time.time()
#                        
            time.sleep(1)
            
            
        return
    
    def pressure(self,port):
        msg = self.mg.my_pressure(port)
        pressure = msg[0]
#        status = msg[1]
        return pressure
    
    def connectPfeiffer(self,com):
        try:
            self.mg = MaxiGauge(com)
            self.connected = True
            print('Connected to Pfeiffer')
        except:
            msg = 'No connection to port %s'%com
            print(msg)
            return msg
            
        print(self.mg.checkDevice())
        time.sleep(0.2)
        return 'connected'
    
    def ClearGraphs(self):
        self.OnClear = True
        return
    
    def connectTemperature(self,com):
        try:
            self.arduino = serial.Serial('COM18', 9600, timeout=10.0)
            print('Connected to temperature reader')
        except:
            msg = 'No connection to port %s'%com
            print(msg)
            return msg
            
        print(self.mg.checkDevice())
        time.sleep(0.2)
        return 'connected'
    
    def ClearGraphs(self):
        self.OnClear = True
        return
        
    def PlotModeFunc(self,mode):
        self.PlotMode = mode
        return
    
#    def workerfunc_work(self,msg):
#        self.Log_work.emit(msg)
#        return
    

            

def main():
    app = QtWidgets.QApplication(sys.argv)
    form = GuiThread()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()
