# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 16:06:57 2020

@author: Gruppe Willitsch
"""

import time
from unidaq import *

def run():
    print('refreshing trap')
    #   #daq card
    boardno = 0
    channel = 0
    voltage = 4.0
    waittime = 10
    
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



if __name__ == '__main__':
    i = 0
    while(i<12):
        run()
        time.sleep(10)
        i+=1
        print(i)
