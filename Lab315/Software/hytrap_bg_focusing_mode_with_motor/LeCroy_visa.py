# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 14:01:28 2018

@author: Gruppe Willitsch
"""

import visa
import sys
import struct
import numpy as np
import time
import readfile as rf



'''

    remote manual: (pdf should be in the dirname_pressure_dataset)
    MAUI Oscilloscopes Remote Control and Automation Manual
    © 2017 Teledyne LeCroy, Inc. All rights reserved.


    waveform template:
    lecroy_waveform_template.txt
'''

'''
    Some important settings
    
    # verticalcoupling 0: A1M, 1: D1M, 2: D50, 3: GND, 4: DC, 5: AC
    # samplemode 0: RealTime, 1: RIS, 2: Sequence, 3: Roll, 4: WaveStream
    # bandwidth 0: Fulll, 1: 20 MHZ, 2: 25 MHZ, 3: 200 MHZ, 4: 250 MHZ, 5: 1 GHZ
    # vertical range VDIV: v/div in V
    # vertical offset OFST: offset in V
    # Probe Attenuation ATTN: Selects the vertical attenuation factor of the probe. Values up to 10000.
    # Interpolation InterpolateType: 0: Linear, 1: Sinx/x
    # Noise Filter EnhanceResType: 0: None, 1: +0.5 bits, 2: +1.0 bits ,... 6: + 3bits
    # Trigger Source SR: 0: Channel 1, 1: Channel 2, 2: Channel 3, 3: Channel 4, 4: Ext, 
        5: Ext10, 6: Line, 7: Pattern
    # Qualifier Source QL
    # Trigger Qualifier (Do we need this?)

        When the trigger qualifier is on, the test set samples the input signal when 
        a trigger is received. It then determines if the input signal was valid by looking 
        at its power level. If the power level during sampling did not meet the requirements 
        of a valid signal, the state returns to wait for trigger without processing the samples. 
        Trigger qualifier is available for GSM/GPRS TX Power and Phase Frequency Error measurements only.
'''


class WaveRunner64MXsB():
    def __init__(self):
        filename='./input/HyT-input.dat' # inputfile
        names,values=rf.readpars_class(filename)    # define inputfile defined variables
    
        for i in range (len(names)):
            try:
                exec(names[i]+"=%s"%(values[i]))
            except (NameError, SyntaxError):
                try:
                    exec(names[i]+"='%s'"%(values[i]))
                except:
                    print('problem after %s = %s'%(names[i-1],values[i-1]))
                    a = 1
                    
        self.scope_ip = self.scope_ip.strip(' ') #strip read in ip from spaces
        
                    
        
    def connect(self):
        try:
            self.rm = visa.ResourceManager()
            self.scope = self.rm.open_resource('TCPIP0::'+self.scope_ip+'::inst0::INSTR')
            self.scope_id = self.scope.query("*IDN?")
            output = 'Connected to '+self.scope_id
            self.scope.write("COMM_HEADER OFF")
            
            self.scope.write("WFSU SP,0,NP,0,FP,0,SN,%i"%(self.numberofsegments)) # waveform setup
            self.scope.write("CFMT DEF9,WORD,BIN") # BYTE <> WORD
        except:
            output = str(sys.exc_info())+'\n connect() failed'
        return output
    
    def get_template(self):
        try:
            self.scope.write('TMPL?')
            output = self.scope.read_bytes(22245)    
        except:
            output = sys.exc_info()+'\n get_template() failed'
        return output
    
    def disconnect(self):
        try:
            self.scope.close() #always close in the end
            self.rm.close() #always close in the end
            output = 'Disconnected from '+self.scope_id
        except:
            output = str(sys.exc_info())+'\n disconnect() failed'
        return output
    
    def initialize(self,channels):
        try:
            if self.numberofsegments == 1:
                self.scope.write("VBS 'app.Acquisition.Horizontal.SampleMode=0'")
            else:
                self.scope.write("VBS 'app.Acquisition.Horizontal.SampleMode=%s'"%(self.samplemode))
            self.scope.write("VBS 'app.Acquisition.Horizontal.NumSegments=%d'"%(self.numberofsegments))
            self.scope.write("WFSU SP,0,NP,0,FP,0,SN,0")
            
            # Trigger
            self.scope.write("TRSE %s,SR,C%i"%(self.trig_type,self.trig_source))
            self.scope.write("C%i:TRLV %f"%(self.trig_source,self.trig_level))
            self.scope.write("C%i:TRSL %s"%(self.trig_source,self.trig_slope))
            self.scope.write("TRDL %s"%(str(self.timeoffset)))   # time delay (offset)
            self.scope.write("TDIV %s"%(str(self.timebase)))     # time division (time base)
            self.scope.write("MSIZ %s"%(self.max_sample))     # Maximum sample size.
            
            self.scope.write("BWL OFF") # no bandwidth limit
            
    #        set up the channels and set up acquisition:
            for channel in channels:
                chstr = "C%i:"%channel
                self.read_descriptor(channel)

                self.scope.write(chstr+"TRA %s"%(self.tracemode))
                self.scope.write(chstr+"ATTN %s"%(self.probeattenuation))  
                self.scope.write(chstr+"VDIV %f"%(eval('self.ch%i_scale'%(channel))))
                self.scope.write(chstr+"OFST %f"%(eval('self.ch%i_verticaloffset'%(channel))))
                self.scope.write(chstr+"CPL %s"%(eval('self.ch%i_verticalcoupling'%(channel))))
       
                acquisition="VBS 'app.acquisition.C%i."%channel
                self.scope.write(acquisition+"AverageSweeps=%i'"%(self.numberofsegments))
                self.scope.write(acquisition+"InterpolateType=%i'"%(self.interpolatetype))
                self.scope.write(acquisition+"EnhanceResType=%i'"%(self.enhancerestype))
                self.scope.write(acquisition+"Invert=%s'"%(eval('self.ch%i_invert'%(channel))))
            output = 'Scope initialized'
        except:
            output = str(sys.exc_info())+'\n initialize() failed'
        return output
    
    def arm(self):
        self.scope.write('*CLS')
        self.scope.write('ARM')
        return
    
    def read_descriptor(self,channel):
        try:
            self.scope.write('C%i:WF? DESC'%channel)
            desc = self.scope.read_bytes(363)
#            print(desc)
#            p 172 ff in manual
            strip = 9   # String 002000348, signaling the beginning of a binary block where nine ASCII
                        # integers are used to give the length of the block (2000348 bytes).
            todel = 7   # position of #9 in bytestring
            header = strip + todel
            desc = desc[header:] # cut away unnecessary part
#            print(desc)
            
#            define some important stuff (see template and remote control manual)
            exec("self.command_type_%i = struct.unpack('i',desc[32:36])[0]"%channel)
            exec("self.first_valid_point_%i = struct.unpack('i',desc[124:128])[0]"%channel)
            exec("self.last_valid_point_%i = struct.unpack('i',desc[128:132])[0]"%channel)
            exec("self.vertical_offset_%i = struct.unpack('f',desc[160:164])[0]"%channel)
            exec("self.vertical_gain_%i = struct.unpack('f',desc[156:160])[0]"%channel)
            exec("self.horizontal_interval_%i = struct.unpack('f',desc[176:180])[0]"%channel)
            exec("self.horizontal_offset_%i = struct.unpack('d',desc[180:188])[0]"%channel)
            exec("self.number_of_points_%i = struct.unpack('i',desc[116:120])[0]"%channel)
            exec("self.vertunit_%i = desc[196:198].decode('ascii')"%channel)
            exec("self.wave_array_count_%i = struct.unpack('l',desc[116:120])[0]"%channel)
            output = desc
        except:
            output = str(sys.exc_info())+'\n read_descriptor() failed'
 
        return output
    
    def get_waveform(self,channel):
        self.read_descriptor(channel)
        xoffset = 0
        if not "self.command_type_%i"%channel:
            return 'Channel %i: descriptor for channel not read'%channel
        
        try:
            # read out the waveform as 8bit integers
            ydata = np.array(self.scope.query_binary_values("C%i:WF? DAT1"%channel,datatype='h',is_big_endian=False))
            # scale the output to the actual value in V
            gain = eval('self.vertical_gain_%i'%channel)
            offset = eval('self.vertical_offset_%i'%channel)
            ydata = gain * ydata - offset
#            print(gain,offset)
#            get x data
            dt = eval('self.horizontal_interval_%i'%channel)
            xoffset = eval('self.horizontal_offset_%i'%channel)
#            dt = 1e-9
            xoffset = 0
            xdata = np.array([xoffset + i * dt for i in range(ydata.shape[0])])
#            print(dt,xoffset)
            output = np.array([xdata,ydata])
        except:
            output = str(sys.exc_info())+'\n get_waveform() failed'
            
        return output
    
    def segment_counter(self,channel):
        try:
            self.scope.write('C%i:WF? DESC'%channel)
            desc = self.scope.read_bytes(363)
            
#            p 172 ff in manual
            strip = 9   # String 002000348, signaling the beginning of a binary block where nine ASCII
                        # integers are used to give the length of the block (2000348 bytes).
            todel = 7   # position of #9 in bytestring
            header = strip + todel
            desc = desc[header:] # cut away unnecessary part
            exec("self.segment_count_%i = struct.unpack('l',desc[144:148])[0]"%channel)
            exec("self.segment_index_%i = struct.unpack('l',desc[140:144])[0]"%channel)
            exec("self.sweeps_per_acq_%i = struct.unpack('l',desc[148:152])[0]"%channel)
            exec("self.subarray_count_%i = struct.unpack('h',desc[174:176])[0]"%channel)


            output = [eval("self.sweeps_per_acq_%i"%channel),
                      eval("self.segment_count_%i"%channel),
                      eval("self.segment_index_%i"%channel),
                      eval("self.subarray_count_%i"%channel)]

        except:
            output = str(sys.exc_info())+'\n segment_counter() failed'
            
        return output

    def clear_sweeps(self):
        try:
            self.scope.write("VBS app.ClearSweeps")
            output = 'Sweeps cleared'
        except:
            output = str(sys.exc_info())+'\n clear_sweeps() failed'
        
        return output



if __name__ == "__main__":

    from matplotlib import pyplot as plt
    import pandas as pd
    scopo = WaveRunner64MXsB()
    print(scopo.connect())
    channel = 1
    print(scopo.initialize([channel]))
    scopo.read_descriptor(channel)
    
    print(scopo.command_type_1)
    print(scopo.first_valid_point_1)
    print(scopo.last_valid_point_1)
    print(scopo.vertical_offset_1)
    print(scopo.vertical_gain_1)
    print(scopo.horizontal_interval_1)
    print(scopo.horizontal_offset_1)
    print(scopo.number_of_points_1)
    print(scopo.vertunit_1)
    print(scopo.wave_array_count_1)
    

    data = scopo.get_waveform(channel)
    print('segments',scopo.segment_counter(channel)[0])
    scopo.clear_sweeps()

    print(scopo.disconnect())

    plt.plot(data[0],data[1])

    datadf = pd.DataFrame(data).T
            
    datadf.to_csv('data.txt')


        