import numpy as np
import sys
from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import QThread, pyqtSignal, Qt

from Mainwindow import Ui_MainWindow
import locale
import time
from xmlrpc.client import ServerProxy
import numpy as np
import unidaq as ud



class GuiThread(QtWidgets.QMainWindow, Ui_MainWindow):
    
    ###pyqtSignals are used to cummunicate between threads, here GuiThread and testprint
    stopper_sig = pyqtSignal()         ###Stopping command, no object (or data) required
    connect_sig = pyqtSignal()         ###Start signal
    set_WL_sig  = pyqtSignal(object)   ###Send the wavelength (WL, aka here just a number) from one thread to another
    set_U_ec1_sig   = pyqtSignal(object)   
    set_U_ec2_sig   = pyqtSignal(object)   
    set_U_comp_v_sig   = pyqtSignal(object)   
    set_U_comp_h_sig   = pyqtSignal(object)   
    set_pid_p_sig   = pyqtSignal(object)   
    set_pid_i_sig   = pyqtSignal(object)   
    set_pid_ched_sig    = pyqtSignal()
    

    def __init__(self):
        super(self.__class__, self).__init__()       ###  Some standard lines to set up the gui
        self.setupUi(self)

        ###Variables
        self.show_nbr = 0                             ### Variable to display the number
        self.lcdnumber.display(str(self.show_nbr))    ### Connect variable to display
        self.set_WL = locale.atof(self.Set_wavelength.text())   ###Read the value. Had a problem with the "," in german language, that is why locale.atof
        
        ###Get the valvues for the voltages
        self.set_u_ec1 = locale.atof(self.ec_voltage_1.text())
        self.set_u_ec2 = locale.atof(self.ec_voltage_2.text())
        self.set_u_comp_v = locale.atof(self.dc_offset_vertical.text())
        self.set_u_comp_h = locale.atof(self.dc_offset_horizontal.text())
        self.set_pid_p = locale.atof(self.pic_control_p.text())
        self.set_pid_i = locale.atof(self.pic_control_i.text())
        self.reset  = True

        ###Buttons actions
        ###Buttons actions
        ###Buttons actions
        self.start_button.clicked.connect(self.StartPIDen)    ####Button names are defined in mainwindow.py
        self.stop_button.clicked.connect(self.EndPIDen)       ####Button names are defined in mainwindow.py
        self.Set_wavelength.valueChanged.connect(self.update_set_wl)  ### Important :)
        self.ec_voltage_1.valueChanged.connect(self.update_set_ec1)  ### Important :)
        self.ec_voltage_2.valueChanged.connect(self.update_set_ec2)  ### Important :)
        self.dc_offset_vertical.valueChanged.connect(self.update_set_comp_v)  ### Important :)
        self.dc_offset_horizontal.valueChanged.connect(self.update_set_comp_h)  ### Important :)
        self.pic_control_p.valueChanged.connect(self.update_set_pid_p)
        self.pic_control_i.valueChanged.connect(self.update_set_pid_i)
        self.PID_check_box.stateChanged.connect(self.pid_control_status)

        ###Threads
        self.lockwl = LaserLock_Thread() 

        ###Connect signals/threads with methods (fcts)
        self.stopper_sig.connect(self.lockwl.stop_locking)
        self.connect_sig.connect(self.lockwl.connect)        
        self.set_WL_sig.connect(self.lockwl.update_wl)
        self.lockwl.wl_printer.connect(self.update_display_wl)
        self.set_U_ec1_sig.connect(self.lockwl.update_u_ec1)
        self.set_U_ec2_sig.connect(self.lockwl.update_u_ec2)
        self.set_U_comp_v_sig.connect(self.lockwl.update_u_comp_v)
        self.set_U_comp_h_sig.connect(self.lockwl.update_u_comp_h)
        self.set_pid_p_sig.connect(self.lockwl.update_pid_p)
        self.set_pid_i_sig.connect(self.lockwl.update_pid_i)
        self.set_pid_ched_sig.connect(self.lockwl.stop_PID)
        

    def StartPIDen(self):                                      ###Start executing
        if self.reset:
            self.connect_sig.emit()
            self.set_WL_sig.emit(self.set_WL)
            self.set_U_ec1_sig.emit(self.set_u_ec1)
            self.lockwl.start()

            self.reset = False
        else:
            self.lockwl.start()
        
        
    def EndPIDen(self):
        self.stopper_sig.emit()
        
    def update_set_wl(self):
        self.set_WL = locale.atof(self.Set_wavelength.text())
        self.set_WL_sig.emit(self.set_WL)
        return
    
    def update_display_wl(self,current_wl):                    ###Update the display of the ui
        self.show_nbr = current_wl
        self.lcdnumber.display(str(np.round(self.show_nbr,6)))


    def update_set_ec1(self):
        self.set_u_ec1 = locale.atof(self.ec_voltage_1.text())
        self.set_U_ec1_sig.emit(self.set_u_ec1)
        return
    def update_set_ec2(self):
        self.set_u_ec2 = locale.atof(self.ec_voltage_2.text())
        self.set_U_ec2_sig.emit(self.set_u_ec2)
        return
    def update_set_comp_v(self):
        self.set_u_comp_v = locale.atof(self.dc_offset_vertical.text())
        self.set_U_comp_v_sig.emit(self.set_u_comp_v)
        return
    def update_set_comp_h(self):
        self.set_u_comp_h = locale.atof(self.dc_offset_horizontal.text())
        self.set_U_comp_h_sig.emit(self.set_u_comp_h)
        return
    
    def update_set_pid_p(self):
        self.set_pid_p = locale.atof(self.pic_control_p.text())
        self.set_pid_p_sig.emit(self.set_pid_p)
        return
    
    def update_set_pid_i(self):
        self.set_pid_i = locale.atof(self.pic_control_i.text())
        self.set_pid_i_sig.emit(self.set_pid_i)
        return
    
    def pid_control_status(self):
        if self.PID_check_box.isChecked():
            self.set_pid_ched_sig.emit()
            print('start')
        else:
            print('stop')
            self.set_pid_ched_sig.emit()
        return


class LaserLock_Thread(QThread):
    wl_printer = pyqtSignal(object)                             ###Used to send the current wl to the GUI Thread
    
    def __init__(self, parent, time):
        QThread.__init__(self)
        self.timeoff = time

    def __init__(self):
        QThread.__init__(self)
        self.kp = 0.5
        self.ki = 0.00425
        self.setwl = 396.959185
        self.piezovoltage = -0.0
        self.integral_window = 100
        self.ud_channel = 1
        self.ud_boardno = 0
        self.sleeptime = 0.2
        self.server_ip = "10.40.10.120"
        self.server_port = 8000
        self.piezo_max = 0.5
        print('Super Duper Ultra Wave Length Client 3000')
        print('ip', self.server_ip)
        print('port', self.server_port)
     
        ####Variables to optimize PI control
        self.counter = 0
        self.value = []
        self.time  = []
        self.pid_killer = 1

        ####DC Endcaps
        self.channel_ec1 = 0
        self.channel_ec2 = 2
        self.u_ec1 = 2.
        self.u_ec2 = 2.
        
        ####DC offset channels and voltages
        self.channel_l = 3
        self.channel_r = 4
        self.channel_q = 5
        self.channel_m = 6
        
        self.u_l = 0
        self.u_r = 0
        self.u_q = 0
        self.u_m = 0
        
        
        ###--- DC offset parameters in 10 mV
        range_v = [0,55.0]
        range_h = [0,-45.0]
        which = 0
        
        self.vertical = 0
        self.horizontal = 0
        


        ###--- Variables for running the script
        self.excecute = False
        
        
    def connect(self):
        link = 'http://' + self.server_ip + ':' + str(self.server_port)
        self.server = ServerProxy(link)
        ts, wl = self.server.get_wavelength()
        self.udaq = ud.UniDaq()
        self.udaq.initalize()
        self.udaq.GetCardInfo(self.ud_boardno)
        self.udaq.ConfigAO(self.ud_boardno,self.ud_channel,3)    ###### Mode defines the resolution, 1 = +/- 5 V, 2 = 0-10, 3 = +/- 10 V
                                                                  ###http://ftp.icpdas.com/pub/cd/iocard/pci/napdos/pci/unidaq/manual/unidaq%20dll%20user%20manual_eng_2.6.pdf
        self.udaq.ConfigAO(self.ud_boardno,self.channel_ec1,3)
        self.udaq.ConfigAO(self.ud_boardno,self.channel_ec2,3)

        self.udaq.ConfigAO(self.ud_boardno,self.channel_l,3)
        self.udaq.ConfigAO(self.ud_boardno,self.channel_r,3)
        self.udaq.ConfigAO(self.ud_boardno,self.channel_q,3)
        self.udaq.ConfigAO(self.ud_boardno,self.channel_m,3)

        return

    def set_endcapvoltage(self):
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_ec1,self.u_ec1)
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_ec2,self.u_ec2)
        return

    def set_piezovoltage(self, voltage):
        self.udaq.WriteAOVoltage(self.ud_boardno,self.ud_channel,voltage)
        return

    def set_dc_offset_voltage(self):
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_l,self.u_l)
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_r,self.u_r)
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_q,self.u_q)
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_m,self.u_m)
        print(self.u_m)
        return
    
    def stop_PID(self):
        if self.pid_killer == 1:
            self.pid_killer = 0
        else:
            self.pid_killer = 1
        return
    
    def update_wl(self,new_wl):
        self.setwl = new_wl
        return
    
    def update_pid_p(self,new_p):
        self.kp = new_p
        print(self.kp)
        return
    
    def update_pid_i(self,new_i):
        self.ki = new_i
        print(self.ki)
        return
        
    def update_u_ec1(self,new_voltage):
        self.u_ec1 = new_voltage
        if self.excecute:
            self.set_endcapvoltage()
        return
    
    def update_u_ec2(self,new_voltage):
        self.u_ec2 = new_voltage
        if self.excecute:
            self.set_endcapvoltage()
        return
    
    def update_u_comp_v(self,new_voltage):
        self.vertical = new_voltage
        self.get_dc_offset_voltages()
        if self.excecute:
            self.set_dc_offset_voltage()
        return
    
    def update_u_comp_h(self,new_voltage):
        self.horizontal = new_voltage
        self.get_dc_offset_voltages()
        if self.excecute:
            self.set_dc_offset_voltage()
        return
        
    def get_dc_offset_voltages(self):
        u_vertical = 0.01*self.vertical       ###### Convert volts into 10 mV
        u_horizontal = 0.01*self.horizontal   ###### Convert volts into 10 mV
        ##Vertical pairs
        self.u_l = u_vertical+u_horizontal    ###### Calculate the voltages of the pairs of electrodes
        self.u_q = u_vertical
        self.u_r = u_horizontal
        self.u_m = 0
        return
    
    def stop_locking(self):
        self.excecute = False

    def run(self):     #def lock(self):
        self.excecute = True
        integral_list = np.zeros(self.integral_window)
        
        self.set_endcapvoltage()
        self.set_dc_offset_voltage()
        
        while self.excecute == True:
            ts, wavelength = self.server.get_wavelength()
            wavelength = np.round(wavelength,6)
            print('eaz'+str(wavelength))
            dwl = wavelength - self.setwl
            if abs(dwl)<5e-6:
                dwl = 0

            tmp = np.insert(integral_list, 0, dwl)
            integral_list = np.delete(tmp,-1)
            integral_value = np.sum(integral_list)

            dv = self.kp * dwl + self.ki * integral_value
            dv = dv*self.pid_killer
            self.piezovoltage = self.piezovoltage + dv


            if np.abs(self.piezovoltage) < self.piezo_max:
                self.set_piezovoltage(self.piezovoltage)
            else:
                print('Voltage not within boundaries')
           
            print(str(np.round(wavelength,6))+'  '+str(dwl)+'  '+str(np.round(self.piezovoltage,4)))#+'  '+str(np.round(self.u_ec1,2)))

            ###Emit the wl to the GUI
            self.wl_printer.emit(wavelength)
 
           
            time.sleep(self.sleeptime)
 #          ####optimize PI
 #          if self.counter<500:
 #               self.value = np.append(self.value,np.round(wavelength,8))
 #               self.time  = np.append(self.time,time.time())
 #               self.counter +=1
 #               print(self.counter)
 
 #           if self.counter==500:
 #               np.save('C:/Users/Gruppe Willitsch/Desktop/data/2020-04-22/ki_0p00425_kp_p0p2_3',self.value)
 #               #np.save('C:/Users/Gruppe Willitsch/Desktop/data/2020-04-22/time_long',self.time)
 #               print('do noth8ing')
 #               self.counter +=1


def main():
    app = QtWidgets.QApplication(sys.argv)
    window = GuiThread()
    window.show()
    app.exec()

if __name__ == '__main__':
    main()
