# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 14:00:50 2019

@author: Gruppe Willitsch

following 
https://pypi.org/project/pyTMCL/
"""
from serial import Serial
from time import sleep
import pyTMCL

class TMCL():
    def __init__(self):
        self.baudrate = 9600
        self.timeout = 1
        self.COM = 'COM26'
        self.connected = False
        self.position = 0
        self.MODULE_ADDRESS = 1
        self.current_max = 110
        self.current_sb = 0
        self.vel_max = 2047
        self.acc_max = 1000
        self.pulse_divisor = 0
        self.ramp_divisor = 7
        self.initialized = False
        # 1382852 steps (around)
        
    def connect(self):
        try: 
            self.serial_port  = Serial(self.COM, self.baudrate) 
                                     #timeout = self.timeout, 
                                     #parity = self.parity,
                                     #stopbits = self.stopbits)
            self.connected = True
            
            ## Create a Bus instance using the open serial port
            self.bus = pyTMCL.connect(self.serial_port)
            self.motor = self.bus.get_motor(self.MODULE_ADDRESS)
            self.initialized = True
            print('1')
            return (1, self.connected)
        except:
            print('2')
            return (0, self.connected)
        
    def initialize(self):
        if not self.initialized:
            return (0, 'not initialized')
        self.motor.axis.max_current = self.current_max
        self.motor.axis.standby_current = self.current_sb

        self.motor.axis.max_positioning_speed = self.vel_max
        self.motor.axis.max_accelleration = self.acc_max
        
        self.motor.axis.pulse_divisor = self.pulse_divisor
        self.motor.axis.ramp_divisor = self.ramp_divisor

        
        print('3')
        return (1, self.connected)
    
    def disconnect(self):
        if self.initialized:
            self.serial_port.close()
            print('4')
        
if __name__ == "__main__":
    
    # some example program that runs the motor back and forth and measures the motor tmeperature
    import time
    import serial
    import datetime
    from matplotlib import pyplot as plt
    import pandas as pd
    import numpy as np

#    arduino = serial.Serial('COM26', 9600, timeout=10.0)
    
#    temperature = arduino.readline()[:-3]
#    print(temperature)
    distance = 1382852
    runs = 20
    waittime = 10
    tmcl = TMCL()
    tmcl.connect()
    tmcl.initialize()
    
    zerotime = datetime.datetime(1970,1,1)
    times = []
    temperatures = []
    for i in range(runs):
        print('hin', i)
        tmcl.motor.move_relative(+distance)
        time.sleep(waittime)
        print('zurueck', i)
        tmcl.motor.move_relative(-distance)
        time.sleep(waittime)
        
#        temperature = float(arduino.readline()[:-3])
#        temperatures.append(temperature)
#        print(temperature)
        now = datetime.datetime.now()
        times.append((now-zerotime).total_seconds())

#        
    timesnp = np.array(times)
    timesnp = timesnp - np.min(timesnp)
    tmcl.disconnect()
#    arduino.close()

#    temperatures2 = [float(temp) for temp in temperatures]

#    plt.plot(timesnp, temperatures2)
#    plt.xlabel('time /s')
#    plt.ylabel('temperature /degC')
#    
#    datadf = pd.DataFrame([timesnp, temperatures2]).T
#    datadf.to_csv('motor_100runs_2.txt', index = False)
