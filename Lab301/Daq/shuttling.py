import time
from xmlrpc.client import ServerProxy
import numpy as np
import unidaq as ud


class owesome():
    def __init__(self):
        self.ud_boardno = 0
        self.ud_channel = 11    ####--- Channels 0-11 are working fine. 
        self.voltage_1 =  0.

    def connect(self):
        self.udaq = ud.UniDaq()
        self.udaq.initalize()
        self.udaq.GetCardInfo(self.ud_boardno)
        self.udaq.ConfigAO(self.ud_boardno,self.ud_channel,3)
        
      
    def output_volts(self):
        self.udaq.WriteAOVoltage(self.ud_boardno,self.ud_channel,self.voltage_1)
        print("test")


class LaserLock():
    def __init__(self):
        self.ud_channel = 1
        self.ud_boardno = 0
        self.sleeptime = 0.2
        self.server_ip = "131.152.105.25"
        self.server_port = 8000
        self.piezo_max = 0.5
        print('Super Duper Ultra Wave Length Client 3000')
        print('ip', self.server_ip)
        print('port', self.server_port)

        
        ####DC Endcaps
        self.channel_ec1 = 6
        self.channel_ec2 = 12
        self.u_ec1 = 2.
        self.u_ec2 = 2.
        
        ####DC offset channels and voltages
        self.channel_l = 9
        self.channel_r = 10
        self.channel_q = 13
        self.channel_m = 14
        
        
        
        ###--- DC offset parameters in 10 mV
        range_v = [120,55.0]
        range_h = [40,-45.0]
        which = 0
        
        vertical = range_v[which]
        horizontal = range_h[which]
        u_vertical = 0.01*vertical
        u_horizontal = 0.01*horizontal
        
        ##Vertical pairs
        self.u_l = u_vertical+u_horizontal
        self.u_q = u_vertical
        
        self.u_r = u_horizontal
        self.u_m = 0       


    def connect(self):
        link = 'http://' + self.server_ip + ':' + str(self.server_port)
        self.server = ServerProxy(link)
        ts, wl = self.server.get_wavelength()
        self.udaq = ud.UniDaq()
        self.udaq.initalize()
        self.udaq.GetCardInfo(self.ud_boardno)
        self.udaq.ConfigAO(self.ud_boardno,self.ud_channel,3)
        ###### Mode defines the resolution, 1 = +/- 5 V, 2 = 0-10, 3 = +/- 10 V
        ###http://ftp.icpdas.com/pub/cd/iocard/pci/napdos/pci/unidaq/manual/unidaq%20dll%20user%20manual_eng_2.6.pdf
        
        self.udaq.ConfigAO(self.ud_boardno,self.channel_ec1,3)
        self.udaq.ConfigAO(self.ud_boardno,self.channel_ec2,3)
        
        self.udaq.ConfigAO(self.ud_boardno,self.channel_l,3)
        self.udaq.ConfigAO(self.ud_boardno,self.channel_r,3)
        self.udaq.ConfigAO(self.ud_boardno,self.channel_q,3)
        self.udaq.ConfigAO(self.ud_boardno,self.channel_m,3)
        
        self.udaq.ConfigAO(self.ud_boardno,1,3)
        print(ts)
        print(wl)
        return
    
    def set_endcapvoltage(self):
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_ec1,self.u_ec1)
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_ec2,self.u_ec2)

    
    def set_dc_offset_voltage(self):
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_l,self.u_l)
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_r,self.u_r)
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_q,self.u_q)
        self.udaq.WriteAOVoltage(self.ud_boardno,self.channel_m,self.u_m)
        return


if __name__ == '__main__':


    
    shuttle = owesome()
    shuttle.connect()
    shuttle.output_volts()
