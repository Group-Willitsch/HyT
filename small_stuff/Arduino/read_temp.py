# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 10:27:51 2020

@author: Claudio ganz allein !!!!! 
"""

##############
## Script listens to serial port and writes contents into a file
##############
## requires pySerial to be installed 
import serial  # sudo pip install pyserial should work
import time
import numpy as np
import pyqtgraph as pg
import datetime
from pytz import reference



# check if file for this month already exists
now = datetime.datetime.now()-datetime.datetime(1970,1,1)
thetime = now.total_seconds()
day_month = time.strftime('%Y-%m-%d', time.localtime(thetime))
filename = 'temp_'+day_month+'.txt'
LogDir = 'C:/Users/Gruppe Willitsch/Desktop/data/temp_log/'
filename = LogDir + filename
try: 
    pressfile = open(filename, 'r')
    print('file already exists dawg')
except FileNotFoundError:
    pressfile = open(filename, 'w+')
    pressfile.write('Time YYYY-mm-dd HH:MM:SS \t t(deg C)\n')
    print('new Log file: %s'%filename)
    pressfile.close()


serial_port = 'COM32'
baud_rate = 9600 #In arduino, Serial.begin(baud_rate)
#write_to_file_path = "./output.txt"

#output_file = open(write_to_file_path, "w+")
ser = serial.Serial(serial_port, baud_rate)
while True:
    line = ser.readline()
    line = line.decode("utf-8") #ser.readline returns a binary, convert to string
    line = line[:-3]
    print(line)
    time.sleep(.2)

    localtime = reference.LocalTimezone()
    now = datetime.datetime.now()
    ts = str(now.strftime("%Y-%m-%d %H:%M:%S" + localtime.tzname(now))).strip('W. Europe Daylight Time') 
    
    if time.strftime('%Y-%m-%d', time.localtime(thetime))!=day_month:
        day_month = time.strftime('%Y-%m-%d', time.localtime(thetime))
        filename = LogDir + 'temp_'+day_month+'.txt'
        pressfile = open(filename, 'w+')
        pressfile.write('Time YYYY-mm-dd HH:MM:SS \t t(deg C)\n')
        print('new Log file: %s'%filename)
        pressfile.close()

    pressfile = open(filename, 'a')
    lasttime = thetime
    #~ time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(thetime))
    string = ts+'\t'+str(line)
    string += '\n'
    pressfile.write(string)
    pressfile.close()

