import socket    # used for TCP/IP communication
import readfile as rf


'''
    16.04.2019
    Class file for pulsare dye laser (fine adjustment) remote control (minimal)
    example see if __name__ == '__main__'  
'''

class pulsare():
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
#                    print('problem after %s = %s'%(names[i-1],values[i-1]))
                    a=1
        self.BUFFER_SIZE = 4096
        self.pulsare_ip = self.pulsare_ip.lstrip(' ')

    def connect(self):
        try:
            self.pulsare = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.pulsare.settimeout(3)
        except socket.error as msg:
            self.pulsare = None
            
            
        try:
            self.pulsare.connect((self.pulsare_ip, self.pulsare_port))

        except :
            self.pulsare.close()
            self.pulsare = None
            
        
        if self.pulsare is None:  
            return (False,'Could not connect to Pulsare!')
        
        msg=b"RemoteConnect\r\n"
        try:
            self.pulsare.sendall(msg)
            resp = self.pulsare.recv(self.BUFFER_SIZE)
        except:
            return (False,'Could not connect to Pulsare!')
        
        return (True,str(resp)+' IP %s!'%self.pulsare_ip)
    
    def setWavelength(self,wavelength):
        wn_to_nm = 10000000.00
        wlval = 1/(wavelength*100)*1e9 #cm^-1 to nm
        msg = b'SetWavelength %.8f\r\n'%wlval
        try:
            self.pulsare.sendall(msg)
#            resp = self.pulsare.recv(BUFFER_SIZE)
            resp = 'dini mama'
        except:
            resp = 'pulsare not connected'
        return resp
    
    def disconnect(self):
        msg=b"RemoteDisconnect\r\n"
        try:
            self.pulsare.sendall(msg)
            resp = self.pulsare.recv(self.BUFFER_SIZE)
        except:
            return (False,'Not connected to Pulsare!')
        
        try:
            self.pulsare.close()
            resp = 'pulsare disconnected'
        except:
            resp = 'pulsare not connected'
            
        return resp
 
            
if __name__ == '__main__':

    import time
    puls = pulsare()
    print(puls.pulsare_ip, puls.pulsare_port)
    print(puls.connect())
    time.sleep(2)

    print(puls.setWavelength(17734.88))
    print(puls.disconnect())
    
