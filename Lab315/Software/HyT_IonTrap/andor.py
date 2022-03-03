#   pyAndor - A Python wrapper for Andor's scientific cameras
#   Copyright (C) 2009  Hamid Ohadi
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.


'''

    Adapted by Claudio von Planta

'''

import platform
from ctypes import *
from PIL import Image
import sys
import time
import numpy as np

"""Andor class which is meant to provide the Python version of the same
   functions that are defined in the Andor's SDK. Since Python does not
   have pass by reference for immutable variables, some of these variables
   are actually stored in the class instance. For example the temperature,
   gain, gainRange, status etc. are stored in the class. """

class Andor:
    def __init__(self):
        
        libname = 'atmcd64d.dll'
        self.dll = cdll.LoadLibrary(libname)
        self.width       = None
        self.height      = None
        self.temperature = None
        self.tempmin     = None
        self.tempmax     = None

        self.set_T       = None
        self.gain        = None
        self.gainRange   = None
        #~ self.status      = ERROR_CODE[error]
        self.verbosity   = True
        self.preampgain  = None
        self.channel     = None
        self.outamp      = None
        self.hsspeed     = None
        self.vsspeed     = None
        self.serial      = None
        self.exposure    = None
        self.accumulate  = None
        self.kinetic     = None
        self.ReadMode    = None
        self.AcquisitionMode = None
        self.scans       = 1
        self.hbin        = 1
        self.vbin        = 1
        self.hstart      = 1
        self.hend        = None
        self.vstart      = 1
        self.vend        = None
        self.cooler      = None
        # Check operating system and load library
        # for Windows
        
        

        
    def __del__(self):
        error = self.dll.ShutDown()
        return
    
    def verbose(self, error, function=''):
        if self.verbosity is True:
#            print ("[%s]: %s" %(function, error))
            return

    def SetVerbose(self, state=True):
        self.verbose = state

    def AbortAcquisition(self):
        error = self.dll.AbortAcquisition()
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def Initialize(self):
        tekst = c_char(0)  
        self.dll.Initialize.restype = c_uint16
        self.dll.Initialize.argtypes = [POINTER(c_char)]

        error = self.dll.Initialize(byref(tekst))
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
    
    def GetDetector(self):
        cw = c_int(0)
        ch = c_int(1)
        self.dll.GetDetector.restype = c_uint16
        self.dll.GetDetector.argtypes = [POINTER(c_int), POINTER(c_int)]
        error = self.dll.GetDetector(byref(cw), byref(ch))
        
        self.width       = cw.value
        self.height      = ch.value
        
        self.hend        = cw
        self.vend        = ch
        return ERROR_CODE[error]
        
    def ShutDown(self):
        error = self.dll.ShutDown()
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
#        libhandle = self.dll._handle
#        windll.kernel32.FreeLibrary(libhandle)
        return ERROR_CODE[error]
        
    def GetCameraSerialNumber(self):
        serial = c_int()
        error = self.dll.GetCameraSerialNumber(byref(serial))
        self.serial = serial.value
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return (self.serial,ERROR_CODE[error])

    def SetReadMode(self, mode):
        #0: Full vertical binning
        #1: multi track
        #2: random track
        #3: single track
        #4: image
        error = self.dll.SetReadMode(mode)
        self.ReadMode = mode
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def SetAcquisitionMode(self, mode):
        #1: Single scan
        #3: Kinetic scan
        error = self.dll.SetAcquisitionMode(mode)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        self.AcquisitionMode = mode
        return ERROR_CODE[error]
        
    def SetNumberKinetics(self,numKin):
        error = self.dll.SetNumberKinetics(numKin)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        self.scans = numKin
        return ERROR_CODE[error]

    def SetNumberAccumulations(self,number):
        error = self.dll.SetNumberAccumulations(number)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def SetAccumulationCycleTime(self,time):
        error = self.dll.SetAccumulationCycleTime(c_float(time))
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def SetKineticCycleTime(self,time):
        error = self.dll.SetKineticCycleTime(c_float(time))
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def SetShutter(self,typ,mode,closingtime,openingtime):
        error = self.dll.SetShutter(typ,mode,closingtime,openingtime)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def SetImage(self,hbin,vbin,hstart,hend,vstart,vend):
        self.hbin = hbin
        self.vbin = vbin
        self.hstart = hstart
        self.hend = hend
        self.vstart = vstart
        self.vend = vend
        
        error = self.dll.SetImage(hbin,vbin,hstart,hend,vstart,vend)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
    
    # cvp
    def GetHeadModel(self):
        head_arr = [0 for i in range(32)]
        HeadModel = (c_char * len(head_arr))(*head_arr)

        self.dll.GetHeadModel.restype = c_uint16
        self.dll.GetHeadModel.argtypes = [POINTER(c_char * len(head_arr))]

        error = self.dll.GetHeadModel(pointer(HeadModel))
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        self.HeadModel = ''
        
        
        for char in HeadModel[:]:
            self.HeadModel += str(bytearray.fromhex(str(hex(char)).lstrip('0x')).decode())
            
        return ERROR_CODE[error]
    
    # cvp
#    def GetHeadCapabilities(self):
#        caps = 
    
    #cvp
    def GetFastestRecommendedVSSpeed(self):
        index = c_int(0)
        speed = c_float(0)
        self.dll.GetFastestRecommendedVSSpeed.restype = c_uint16
        self.dll.GetFastestRecommendedVSSpeed.argtypes = [POINTER(c_int), POINTER(c_float)]
        error = self.dll.GetFastestRecommendedVSSpeed(byref(index), byref(speed))
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        
        self.VSNumber = index
        self.VSSpeed = speed
        return ERROR_CODE[error]
    
    # cvp not available for our camera!!
    def SetBaselineClamp(self, state):
        self.dll.SetBaselineClamp.restype = c_uint16
        self.dll.SetBaselineClamp.argtypes = [c_int]
        error = self.dll.SetBaselineClamp(c_int(state))
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
    
    #cvp
    def GetImage(self):
        self.StartAcquisition()
        msg = self.GetStatus()

        while msg!='DRV_IDLE':
            msg = self.GetStatus()
        data = []
        msg = self.GetAcquiredData(data)
        
#        self.SaveAsTxt("raw.txt")
        datanew = np.array(data).reshape(int(cam.height/hbin),int(cam.width/vbin))
        self.imagedata = datanew
        return msg
    
    def GetImageData(self):
        self.StartAcquisition()
        msg = self.GetStatus()

        while msg!='DRV_IDLE':
            msg = self.GetStatus()
        data = []
        msg = self.GetAcquiredData(data)
        self.imagedata = data
#        self.SaveAsTxt("raw.txt")
        return msg

    def StartAcquisition(self):
        error = self.dll.StartAcquisition()
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
    
    def WaitForAcquisition(self):
        error = self.dll.WaitForAcquisition()
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)

        return ERROR_CODE[error]

    def GetAcquiredData(self,imageArray):
        if (self.ReadMode==4):
            if (self.AcquisitionMode==1):
                dim = int(self.width * self.height / self.hbin / self.vbin)

            elif (self.AcquisitionMode==3):
                dim = self.width * self.height / self.hbin / self.vbin * self.scans
        elif (self.ReadMode==3 or self.ReadMode==0):
            if (self.AcquisitionMode==1):
                dim = self.width
            elif (self.AcquisitionMode==3):
                dim = self.width * self.scans

#        print('dim = %s'%dim)
        cimageArray = c_int * dim
        cimage = cimageArray()
        error = self.dll.GetAcquiredData(pointer(cimage),dim)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)

        for i in range(len(cimage)):
            imageArray.append(cimage[i])

        self.imageArray = imageArray[:]
#        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def SetExposureTime(self, time):
        error = self.dll.SetExposureTime(c_float(time))
        self.exposure = time
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
        
    def GetAcquisitionTimings(self):
        exposure   = c_float()
        accumulate = c_float()
        kinetic    = c_float()
        error = self.dll.GetAcquisitionTimings(byref(exposure),byref(accumulate),byref(kinetic))
        self.exposure = exposure.value
        self.accumulate = accumulate.value
        self.kinetic = kinetic.value
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def SetSingleScan(self):
        self.SetReadMode(4)
        self.SetAcquisitionMode(1)
        self.SetImage(1,1,1,self.width,1,self.height)

    def SetCoolerMode(self, mode):
        error = self.dll.SetCoolerMode(mode)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
        
    def SetFanMode(self, mode):
        #0: fan on full
        #1: fan on low
        #2: fna off
        error = self.dll.SetFanMode(mode)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def SaveAsBmp(self, path):
        im=Image.new("RGB",(self.width,self.height),"white")
        pix = im.load()

        for i in range(len(self.imageArray)):
            (row, col) = divmod(i,self.width)
            picvalue = int(round(self.imageArray[i]*255.0/65535))
            pix[col,row] = (picvalue,picvalue,picvalue)

        im.save(path,"BMP")

    def SaveAsTxt(self, path):
        file = open(path, 'w')

        for line in self.imageArray:
            file.write("%g\n" % line)

        file.close()

    def SetImageRotate(self, iRotate):
        error = self.dll.SetImageRotate(iRotate)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)

    def SaveAsBmpNormalised(self, path):

        im=Image.new("RGB",(self.width,self.height),"white")
        pix = im.load()

        maxIntensity = max(self.imageArray)

        for i in range(len(self.imageArray)):
            (row, col) = divmod(i,self.width)
            picvalue = int(round(self.imageArray[i]*255.0/maxIntensity))
            pix[col,row] = (picvalue,picvalue,picvalue)

        im.save(path,"BMP")
        
    def SaveAsFITS(self, filename, type):
        error = self.dll.SaveAsFITS(filename, type)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def CoolerON(self):
        error = self.dll.CoolerON()
        self.cooler = 1
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def CoolerOFF(self):
        error = self.dll.CoolerOFF()
        self.cooler = 0
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def IsCoolerOn(self):
        iCoolerStatus = c_int()
        self.cooler = iCoolerStatus
        error = self.dll.IsCoolerOn(byref(iCoolerStatus))
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return iCoolerStatus.value

    def GetTemperature(self):
        ctemperature = c_int()
        error = self.dll.GetTemperature(byref(ctemperature))
        self.temperature = ctemperature.value
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
    
    def GetTemperatureRange(self):
        ctempmin = c_int()
        ctempmax = c_int()

        error = self.dll.GetTemperatureRange(byref(ctempmin),byref(ctempmax))
        self.tempmin = ctempmin.value
        self.tempmax = ctempmax.value

        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
    
    def StabilizeTemperature(self,tset):
        self.SetTemperature(tset)
        self.CoolerON()
        msg=self.GetTemperature()
        print('Cooling started: %.2f degC'%tset)
        
        while msg!='DRV_TEMP_STABILIZED':
            msg = self.GetTemperature()
            print('...')
            time.sleep(1)
        
        print('Temperature stabilized')

    def SetTemperature(self,temperature):
        #ctemperature = c_int(temperature)
        #error = self.dll.SetTemperature(byref(ctemperature))
        error = self.dll.SetTemperature(temperature)
        self.set_T = temperature
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def GetEMCCDGain(self):
        gain = c_int()
        error = self.dll.GetEMCCDGain(byref(gain))
        self.gain = gain.value
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
     
    def SetEMCCDGainMode(self, gainMode):
#        error = self.dll.SetEMCCDGainMode(gainMode)
        error = self.dll.SetEMGainMode(gainMode)

        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]   
        
    def SetEMCCDGain(self, gain):
        error = self.dll.SetEMCCDGain(gain)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
        
    def SetEMAdvanced(self, gainAdvanced):
        error = self.dll.SetEMAdvanced(gainAdvanced)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def GetEMGainRange(self):
        low = c_int()
        high = c_int()
        error = self.dll.GetEMGainRange(byref(low),byref(high))
        self.gainRange = (low.value, high.value)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
      
    def GetNumberADChannels(self):
        noADChannels = c_int()
        error = self.dll.GetNumberADChannels(byref(noADChannels))
        self.noADChannels = noADChannels.value
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def GetBitDepth(self):
        bitDepth = c_int()

        self.bitDepths = []

        for i in range(self.noADChannels):
            self.dll.GetBitDepth(i,byref(bitDepth))
            self.bitDepths.append(bitDepth.value)

    def SetADChannel(self, index):
        error = self.dll.SetADChannel(index)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        self.channel = index
        return ERROR_CODE[error]  
        
    def SetOutputAmplifier(self, index):
        error = self.dll.SetOutputAmplifier(index)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        self.outamp = index
        return ERROR_CODE[error]
        
    def GetNumberHSSpeeds(self):
        noHSSpeeds = c_int()
        error = self.dll.GetNumberHSSpeeds(self.channel, self.outamp, byref(noHSSpeeds))
        self.noHSSpeeds = noHSSpeeds.value
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def GetHSSpeed(self):
        HSSpeed = c_float()

        self.HSSpeeds = []

        for i in range(self.noHSSpeeds):
            self.dll.GetHSSpeed(self.channel, self.outamp, i, byref(HSSpeed))
            self.HSSpeeds.append(HSSpeed.value)
            
    def SetHSSpeed(self, itype, index):
        error = self.dll.SetHSSpeed(itype,index)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        self.hsspeed = index
        return ERROR_CODE[error]
        
    def GetNumberVSSpeeds(self):
        noVSSpeeds = c_int()
        error = self.dll.GetNumberVSSpeeds(byref(noVSSpeeds))
        self.noVSSpeeds = noVSSpeeds.value
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def GetVSSpeed(self):
        VSSpeed = c_float()

        self.VSSpeeds = []

        for i in range(self.noVSSpeeds):
            self.dll.GetVSSpeed(i,byref(VSSpeed))
            self.preVSpeeds.append(VSSpeed.value)

    def SetVSSpeed(self, index):
        error = self.dll.SetVSSpeed(index)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        self.vsspeed = index
        return ERROR_CODE[error] 
    
    def GetNumberPreAmpGains(self):
        noGains = c_int()
        error = self.dll.GetNumberPreAmpGains(byref(noGains))
        self.noGains = noGains.value
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def GetPreAmpGain(self):
        gain = c_float()

        self.preAmpGain = []

        for i in range(self.noGains):
            self.dll.GetPreAmpGain(i,byref(gain))
            self.preAmpGain.append(gain.value)

    def SetPreAmpGain(self, index):
        error = self.dll.SetPreAmpGain(index)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        self.preampgain = index
        return ERROR_CODE[error]

    def SetTriggerMode(self, mode):
        error = self.dll.SetTriggerMode(mode)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def GetStatus(self):
        status = c_int()
        error = self.dll.GetStatus(byref(status))
        self.status = ERROR_CODE[status.value]
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return self.status
        
    def GetSeriesProgress(self):
        acc = c_long()
        series = c_long()
        error = self.dll.GetAcquisitionProgress(byref(acc),byref(series))
        if ERROR_CODE[error] == "DRV_SUCCESS":
            return series.value
        else:
            return None
             
    def GetAccumulationProgress(self):
        acc = c_long()
        series = c_long()
        error = self.dll.GetAcquisitionProgress(byref(acc),byref(series))
        if ERROR_CODE[error] == "DRV_SUCCESS":
            return acc.value
        else:
            return None
        
    def SetFrameTransferMode(self, frameTransfer):
        error = self.dll.SetFrameTransferMode(frameTransfer)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
        
    def SetShutterEx(self, typ, mode, closingtime, openingtime, extmode):
        error = self.dll.SetShutterEx(typ, mode, closingtime, openingtime, extmode)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
        
    def SetSpool(self, active, method, path, framebuffersize):
        error = self.dll.SetSpool(active, method, c_char_p(path), framebuffersize)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]

    def SetSingleTrack(self, centre, height):
        error = self.dll.SetSingleTrack(centre, height)
        self.verbose(ERROR_CODE[error], sys._getframe().f_code.co_name)
        return ERROR_CODE[error]
    
    def SetDemoReady(self):
        error = self.SetSingleScan()
        error = self.SetTriggerMode(0)
        error = self.SetShutter(1,0,30,30)
        error = self.SetExposureTime(0.01)
        return error
    
    def SetBinning(self,binningmode):
        if (binningmode==1):
            self.SetImage(1,1,1,self.width,1,self.height)
        elif (binningmode==2):
            self.SetImage(2,2,1,self.width,1,self.height)
        elif (binningmode==4):
            self.SetImage(4,4,1,self.width,1,self.height)
        else:
            self.verbose("Binning mode not found")

ERROR_CODE = {
    20001: "DRV_ERROR_CODES",
    20002: "DRV_SUCCESS",
    20003: "DRV_VXNOTINSTALLED",
    20006: "DRV_ERROR_FILELOAD",
    20007: "DRV_ERROR_VXD_INIT",
    20010: "DRV_ERROR_PAGELOCK",
    20011: "DRV_ERROR_PAGE_UNLOCK",
    20013: "DRV_ERROR_ACK",
    20024: "DRV_NO_NEW_DATA",
    20026: "DRV_SPOOLERROR",
    20034: "DRV_TEMP_OFF",
    20035: "DRV_TEMP_NOT_STABILIZED",
    20036: "DRV_TEMP_STABILIZED",
    20037: "DRV_TEMP_NOT_REACHED",
    20038: "DRV_TEMP_OUT_RANGE",
    20039: "DRV_TEMP_NOT_SUPPORTED",
    20040: "DRV_TEMP_DRIFT",
    20050: "DRV_COF_NOTLOADED",
    20053: "DRV_FLEXERROR",
    20066: "DRV_P1INVALID",
    20067: "DRV_P2INVALID",
    20068: "DRV_P3INVALID",
    20069: "DRV_P4INVALID",
    20070: "DRV_INIERROR",
    20071: "DRV_COERROR",
    20072: "DRV_ACQUIRING",
    20073: "DRV_IDLE",
    20074: "DRV_TEMPCYCLE",
    20075: "DRV_NOT_INITIALIZED",
    20076: "DRV_P5INVALID",
    20077: "DRV_P6INVALID",
    20083: "P7_INVALID",
    20089: "DRV_USBERROR",
    20091: "DRV_NOT_SUPPORTED",
    20095: "DRV_INVALID_TRIGGER_MODE",
    20099: "DRV_BINNING_ERROR",
    20990: "DRV_NOCAMERA",
    20991: "DRV_NOT_SUPPORTED",
    20992: "DRV_NOT_AVAILABLE"
}



if __name__ == '__main__':
    from unidaq import *
    import time
    import threading
    import numpy as np
    from matplotlib import pyplot as plt
    from scipy.ndimage import gaussian_filter1d

    def shop_bin(data):
        x_mesh,y_mesh = np.meshgrid(np.linspace(0,999,1000),np.linspace(0,999,1000))
        x_mesh = x_mesh.ravel()
        y_mesh = y_mesh.ravel()
        data = data[:1000,:1000]
        data_hist,x,y = np.histogram2d(x_mesh,y_mesh,weights=data.ravel(),bins=(200,200))
        data_hist = data_hist.T
        return(data_hist[60:150,50:160])
        
    cam = Andor()
    print(cam.Initialize())
    print(cam.GetDetector())
    print(cam.GetHeadModel())
    print(cam.HeadModel)
    print(cam.SetAcquisitionMode(1))
    print(cam.SetReadMode(4))
    print(cam.GetFastestRecommendedVSSpeed())
    print(cam.SetVSSpeed(cam.VSNumber))
    
    STemp = 0
    HSnumber = 0
    ADnumber = 0
    
    errorValue = cam.GetNumberADChannels()
    nAD = cam.noADChannels
    
    for iAD in range(nAD):
        index = c_int(0)
        cam.dll.GetNumberHSSpeeds(c_int(iAD), c_int(0), byref(index))
        
        for iSpeed in range(index.value):
            speed = c_float(0)
            cam.dll.GetHSSpeed(c_int(iAD), c_int(0), c_int(iSpeed), byref(speed))
            
            if speed.value > STemp:
                STemp = speed.value
                HSnumber = iSpeed
                ADnumber = iAD
    
#    print(STemp, HSnumber, ADnumber)
    print(cam.SetADChannel(ADnumber))
    print(cam.SetHSSpeed(0,HSnumber))
    
    print(cam.GetStatus())
    print(cam.GetTemperatureRange())
    print('temprange:', cam.tempmin, cam.tempmax)
    msg=cam.GetTemperature()
    cam.GetTemperatureRange()
    temperature=-20
    cam.StabilizeTemperature(temperature) #selfmade temperature set routine
    cam.GetTemperature()
    print('current temperature in degC:', cam.temperature)
#    print(cam.SetShutter(typ, mode, closingtime, openingtime))
    
    print(cam.SetTriggerMode(0))
    
    exptime=1.0
    gain = 4066
    state = 1
    cam.SetExposureTime(exptime)
    cam.SetTriggerMode(0)
    cam.GetEMGainRange()
    print('gainrange', cam.gainRange)
    print('emccd mode', cam.SetEMAdvanced(state))
    print('emgainmode', cam.SetEMCCDGainMode(1))
    print('egain get', cam.GetEMCCDGain())
    cam.GetEMGainRange()
    print('gainrange', cam.gainRange)
    print('gainval', cam.SetEMCCDGain(gain))
    print(cam.gain)

    print(cam.GetAcquisitionTimings())
    print('Acq timings: exposure time, accumulation time, kinetic:\n %f \n %f \n %f'
      %(cam.exposure,cam.accumulate,cam.kinetic))
    
    hbin=1
    vbin=1
    hstart=1
    hend=cam.width
    vstart=1
    vend=cam.height
    print(cam.SetImage(hbin,vbin,hstart,hend,vstart,vend))
    
    # measure reference
    
    
    thread = threading.Thread(target=cam.GetImage)
    thread.daemon = True
    thread.start()
    
    print(cam.GetImage())
    dataref = cam.imagedata
    dataref = shop_bin(dataref)
    
    
    

        
#    smoothval = 2
    
#    plt.close('all')
#    gs = plt.GridSpec(1,2)
#    fig = plt.figure(figsize=(12,5))
#    ax = fig.add_subplot(gs[0,0])
#    img1 = plt.imshow(dataref,cmap='seismic',vmin = -7e3,vmax = 7e3)
#    plt.colorbar()
#    ax = fig.add_subplot(gs[0,1])
#    
#    img2, = plt.plot(np.sum(dataref,axis=0))
#    plt.ylim([-30000,30000])
#    img3, = plt.plot(np.sum(dataref,axis=1))
#    plt.ylim([-30000,30000])
#
#    plt.ion()

#    print('Enter 1 to continue, 0  to stop')
#    cont = int(input())
#    
#    while cont:
#        plt.pause(1e-3)
#        print('background')
#        
#        #daq card
#        boardno = 0
#        channel = 0
#        voltage = 4.0
#        udaq = UniDaq()
#        udaq.initalize()
#        waittime = 4
#        
#        print(udaq.GetCardInfo(boardno))
#        print(udaq.ConfigAO(boardno,channel,3))
#        print(udaq.WriteAOVoltage(boardno,channel,0))
#        udaq.close()
#        print('ramping down rf voltage')
#        time.sleep(waittime)
#
#        print(cam.GetImage())
#        dataref = cam.imagedata
#        dataref = shop_bin(dataref)    
#
#        print('enter number of signal to continue (0 to leave)')
#        number = int(input())
#        
#        udaq = UniDaq()
#        udaq.initalize()
#        
#        print(udaq.GetCardInfo(boardno))
#        print(udaq.ConfigAO(boardno,channel,3))
#        print(udaq.WriteAOVoltage(boardno,channel,voltage))
#        udaq.close()
#        print('ramping up rf voltage')
#        
#        time.sleep(waittime)
#        
#        while number > 0.5:
#            plt.pause(1e-3)
#            for i in range(number):
#                plt.pause(1e-3)
#                print('signal')
#                print(cam.GetImage())
#                datasig = cam.imagedata
#                datasig = shop_bin(datasig)
#                
#                img1.set_data(datasig - dataref)
##                plt.draw()
#                plt.pause(1e-3)
#                
#                img2.set_xdata(np.arange(0,len(datasig[0,]),1))
#                img2.set_ydata(gaussian_filter1d(np.sum(datasig,axis=0)-np.sum(dataref,axis=0),smoothval))
#                img3.set_xdata(np.arange(0,len(datasig),1))
#                img3.set_ydata(gaussian_filter1d(np.sum(datasig,axis=1)-np.sum(dataref,axis=1),smoothval))
#                plt.pause(1e-3)
#                
#                plt.draw()
#            print('enter number of signal to continue (0 for next bg)')
#            number = int(input())
#            
#        print('Enter 1 to continue, 0  to stop')
#        cont = int(input())
#            
#            
#    print('Shutting down...')
    
    cam.ShutDown()

