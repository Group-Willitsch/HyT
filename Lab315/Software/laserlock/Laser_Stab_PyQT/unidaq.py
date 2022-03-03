# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:47:44 2019

@author: Claudantus
"""

'''
    17.04.2019
    Unidaq dll comm

'''



from ctypes import *
from ctypes.wintypes import *
import win32event # for interrupt events


ERRORS = {0: 'Successfully',
          1: 'Open Driver Failure',
          2: 'Plug&Play Failure',
          3: 'Driver was not open',
          4: 'Get Driver Version Failure',
          5: 'Board Number Error',
          6: 'Cannot Find Board',
          7: 'Board Mapping Error',
          8: 'Configure DIO Port Failure',
          9: 'Invalid Address',
          10: 'Invalid Size',
          11: 'Invalid Port Number',
          12: 'Model Is Not Supported',
          13: 'Function Is Not Supported',
          14: 'Invalid Channel Number',
          15: 'Invalid Value',
          16: 'Invalid Mode',
          17: 'Data Not Ready',
          18: 'Timeout',
          19: 'Cannot Find Configuration Code Index',
          20: 'ADC Timeout',
          21: 'Cannot Find Board Index',
          22: 'Invaild Setting',
          23: 'Allocate Memory Space Failed',
          24: 'Install Interrupt Event Failure',
          25: 'Install Interrupt IRQ Failure',
          26: 'Remove Interrupt IRQ Failure',
          27: 'Clear Interrupt Count Failure',
          28: 'Get System Buffer Failure',
          29: 'Call CreateEvent() Failed',
          30: 'Resolution IS Not Supported',
          31: 'Call CreateThread() Failed'}

#EVENTS
IXUD_HARDWARE_INT = WORD(1)     #Device generated a Hardware interrupt
IXUD_APC_READY_INT = WORD(2)    #Interrupt generated from analog input data ready.
IXUD_ACTIVE_LOW = WORD(4)       #Interrupt generated from digital input port failing edge
IXUD_ACTIVE_HIGH = WORD(8)      #Interrupt generated from digital input port raising edge.




class UniDaq():
    def __init__(self):
        self.dwModelNo = DWORD(0x800400)
        self.wBoardNo = WORD(0)
        self.initialized = False
        self.AOconfig = False
        self.wCfgCode = 3 # for PIO-DA4/8/16: +/- 10 V, 16 for 0-20 V
        self.CardInfo = False

        
    def initalize(self):
        print('udaq init')
        libname = 'UniDAQ.dll'
        self.dll = cdll.LoadLibrary(libname)
        
        self.dwDLLVer =  DWORD(0)
        self.dll.Ixud_GetDllVersion.restype = WORD
        self.dll.Ixud_GetDllVersion.argtypes = [POINTER(DWORD)]
        error = self.dll.Ixud_GetDllVersion(byref(self.dwDLLVer))
        if error != 0:
            print(ERRORS[error])
#        print('DLL version: %s' %self.dwDLLVer.value)

        self.wTotalBoards = WORD(0)
        self.dll.Ixud_DriverInit.restype = WORD
        self.dll.Ixud_DriverInit.argtypes = [POINTER(WORD)]
        error = self.dll.Ixud_DriverInit(byref(self.wTotalBoards))
        print('Number of boards: %i' %self.wTotalBoards.value)
        if error != 0:
            print(ERRORS[error])
#
        self.initialized = True
        return (0, 'DAQ initialized')
    
    def ConfigAO(self, wBoardNo, wChannel, wCfgCode):
        if not self.initialized:
            return (1, 'DAQ not initialized')
        
        self.dll.Ixud_ConfigAO.restype = c_uint16
        self.dll.Ixud_ConfigAO.argtypes = [c_uint16, c_uint16, c_uint16]
        
        error = self.dll.Ixud_ConfigAO(c_uint16(wBoardNo), 
                                       c_uint16(wChannel), 
                                       c_uint16(wCfgCode))
        if error != 0:
            print(ERRORS[error])
            
        if error == 0:
            self.AOconfig = True
            
#        print(ERRORS[error])
        
        return (0, 'Ch %i configured with code %i'%(wChannel, wCfgCode))
        
    def WriteAOVoltage(self, wBoardNo, wChannel, fValue):
        if not self.initialized:
            return (1, 'DAQ not initialized')
        
        if not self.AOconfig:
            return (1, 'AO not configured')
        
        self.dll.Ixud_WriteAOVoltage.restype = c_uint16
        self.dll.Ixud_WriteAOVoltage.argtypes = [c_uint16, c_uint16, c_float]
        
        error = self.dll.Ixud_WriteAOVoltage(c_uint16(wBoardNo), 
                                             c_uint16(wChannel), 
                                             c_float(fValue))   
        if error != 0:
            print(ERRORS[error])
            
        return (0, 'AOV: %.2f V on Ch %i'%(fValue, wChannel))
    
    def GetCardInfo(self, wBoardNo):
        if not self.initialized:
            return (1, 'DAQ not initialized')
        
        dev_arr = [0 for i in range(20)]
        card_arr = [0 for i in range(18)]
        model_arr = [0 for i in range(19)]
        
        wBoardNo = c_uint16(wBoardNo)
        sDevInfo = (c_uint32 * len(dev_arr))(*dev_arr)
        sCardInfo = (c_uint32 * len(card_arr))(*card_arr)
        szModelName = (c_uint8 * len(model_arr))(*model_arr)
        
        self.dll.Ixud_GetCardInfo.restype = c_uint16
        self.dll.Ixud_GetCardInfo.argtypes = [c_uint16, POINTER(c_uint32 * len(dev_arr)), POINTER(c_uint32 * len(card_arr)), POINTER(c_uint8 * len(model_arr))]
        error = self.dll.Ixud_GetCardInfo(wBoardNo, pointer(sDevInfo), pointer(sCardInfo), pointer(szModelName))
#       
        if error != 0:
            print(ERRORS[error])        
            
        self.sDevInfo = sDevInfo
        self.sCardInfo = sCardInfo
        #   
        divisor = 65536
        self.AOChannels , self.AIChannels = divmod(self.sCardInfo[3], divisor)
        
        # get model name from list of integers
        self.szModelName = ''
        for char in szModelName[:]:
            self.szModelName += str(bytearray.fromhex(str(hex(char)).lstrip('0x')).decode())
        self.szModelName = self.szModelName.strip(' ')
        
#        print('Model ' + self.szModelName)
#        print('channels: %i'%self.AOChannels)
        
        if error == 0:
            self.CardInfo = True
        
        return (0, self.szModelName)
    
    def ReadDI(self, wBoardNo, wPortNo):
        dwDIVal = DWORD(0)
        self.dll.Ixud_ReadDI.restype = c_uint16
        self.dll.Ixud_ReadDI.argtypes = [WORD, WORD, POINTER(DWORD)]
        error = self.dll.Ixud_ReadDI(WORD(wBoardNo), WORD(wPortNo), byref(dwDIVal))
        if error != 0:
            print(ERRORS[error])

        return (0, dwDIVal.value)
    
    def SetEventCallback(self, wBoardNo, wEventTypeMask, wInterruptSource, CallbackFun):
        '''
        Enable the callback function on interrupt event event, 
        when stop the callback function, it must call Ixud_RemoveEventCallback function to disable it.
        Not used so far, but does not give errors...
        '''
#        hEvent = win32event.CreateEvent(None, 0, 0, None)
        hEvent = HANDLE()
        wInterruptSource = WORD(wInterruptSource)
        dwCallBackParameter = DWORD(0)
        self.dll.Ixud_SetEventCallback.restype = c_uint16
        self.dll.Ixud_SetEventCallback.argtypes = [WORD, WORD, WORD, POINTER(HANDLE), c_void_p, DWORD]
        
        error = self.dll.Ixud_SetEventCallback(wBoardNo, wEventTypeMask , wInterruptSource,
                                               byref(hEvent), CallbackFun, dwCallBackParameter)
        if error != 0:
            print(ERRORS[error])
            
        return(0, dwCallBackParameter)
        
    def InstallIrq(self, wBoardNo, dwInterruptMask):
        self.dll.Ixud_InstallIrq.restype = c_uint16
        self.dll.Ixud_InstallIrq.argtypes = [WORD, DWORD]
        error = self.dll.Ixud_InstallIrq(WORD(wBoardNo), DWORD(dwInterruptMask))
        if error != 0:
            print(ERRORS[error])
            
        return(0, 1)
        
    def callback0(self):
        print('el callobacko')
        return 0
        
    def close(self):
        self.dll.Ixud_DriverClose()
        return
    
    
if __name__ == '__main__':
    import time
    udaq = UniDaq()
    udaq.initalize()
    
    boardno = 0
    channel = 0
    voltage = 3.0
    
    print(udaq.GetCardInfo(boardno))
    print(udaq.ConfigAO(boardno,channel,3))
    print(udaq.WriteAOVoltage(boardno,channel,3.0))
#    time.sleep(1e-4)
#    print(udaq.WriteAOVoltage(boardno,channel,voltage))
    
#    value = 0
#    while value == 0:
#        error, value = udaq.ReadDI(0,0)
#        print(value)
#        time.sleep(1)
    print('finito')
#    Example for SetEventCallback
#    wEventTypeMask = IXUD_ACTIVE_LOW
#    wInterruptSource = 0
#    dwInterruptMask = 1
#    c_callback = CFUNCTYPE(None)(udaq.callback0)
#    print(udaq.SetEventCallback(boardno, wEventTypeMask, wInterruptSource, c_callback))
#    print(udaq.InstallIrq(boardno, dwInterruptMask))
#    error, value = udaq.ReadDI(0,0)
#    print(value)
#    while True:
#        time.sleep(1e-3)
#    print(udaq.WriteAOVoltage(boardno,channel,0))

    udaq.close()
    
