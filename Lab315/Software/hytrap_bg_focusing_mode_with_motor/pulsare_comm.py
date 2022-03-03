import socket  # used for TCP/IP communication
import readfile as rf

'''
    16.04.2019
    Class file for pulsare dye laser (fine adjustment) remote control (minimal)
    example see if __name__ == '__main__'  
Default values stored in HyT-input.dat as of 16/07/2021
####################################################
# PULSARE pulsed laser
####################################################
#
pulsare_ip = 10.40.10.113
pulsare_port = 65510
wavelength = 17734.84 (in cm **-1)

'''


class pulsare():
    def __init__(self):
        '''
        set buffer size, ip address and default wavelength
        '''
        self.input_filename = './input/HyT-input.dat'  # inputfile
        names, values = rf.readpars_class(self.input_filename)  # define inputfile defined variables
        for i in range(len(names)):
            try:
                exec(names[i] + "=%s" % (values[i]))
            except (NameError, SyntaxError):
                try:
                    exec(names[i] + "='%s'" % (values[i]))
                except:
                    a=0
#                    print("Some funny problem while reading the input file. ")
        self.BUFFER_SIZE = 4096
        self.pulsare_ip = self.pulsare_ip.lstrip(' ')
        self.wavelength_cmminus1 = self.wavelength # stores wavelength in cm**-1
        self.wavelength_nm = 1e7/self.wavelength_cmminus1 # stores wavelength in nm
        print('Ip address:', self.pulsare_ip, '\t port:', self.pulsare_port)
        print('Default input wavelength:', self.wavelength_cmminus1, ' cm**-1  or ', self.wavelength_nm, ' nm')

        # connect to the laser via TCT/IP (must be enable on the PulsareV3 program first)
        [connection_outcome, connection_reply] = self.connect()

        # set default wavelengths (maybe dangerous)
        self.setWavelength(self.wavelength_cmminus1)


    def connect(self):
        try:
            self.pulsare = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.pulsare.settimeout(10)
        except socket.error as msg:
            self.pulsare = None

        try:
            self.pulsare.connect((self.pulsare_ip, self.pulsare_port))
        except:
            self.pulsare.close()
            self.pulsare = None
            print('Could not connect to Pulsare!')
            return (False, 'Could not connect to Pulsare!')
        else:
            print('No exceptions while connecting to Pulsare')

        msg = b"GetRemoteStatus\r\n"
        try:
            self.pulsare.sendall(msg)
            resp = self.pulsare.recv(self.BUFFER_SIZE)
            print('GetRemoteStatus reply: ', resp)
        except:
            return (False, 'Could not GetRemoteStatus !')

        msg = b"RemoteConnect\r\n"
        try:
            self.pulsare.sendall(msg)
            resp = self.pulsare.recv(self.BUFFER_SIZE)
            print('RemoteConnect reply: ', resp)
        except:
            return (False, 'Could not connect to Pulsare!')

        msg = b"GetStatus\r\n"
        try:
            self.pulsare.sendall(msg)
            temp_resp = self.pulsare.recv(self.BUFFER_SIZE)
            print('GetStatus reply: ', temp_resp)
        except:
            return (False, 'Could not connect to Pulsare!')

        return (True, str(resp))

    def setWavelength(self, wavelength_cmminus1, use_nanometers=False):
        '''
        Command like “SetWavelength 540.1234\r\n" with a space between command and wavelength
        @param wavelength: must be in cm**-1! Gets converted later in nanometers
        with a dot separating the decimal aprt
        @return: None
        '''
        if use_nanometers:
            wavelength_nm = wavelength_cmminus1
            print('Setting wavelength to ', wavelength_nm, 'nm ')
        else:  # I have to convert the input wavelength from cm**-1 to nm
            wavelength_nm = 1e7/wavelength_cmminus1
            print('Setting wavelength to ', wavelength_nm, 'nm  or ', wavelength_cmminus1, ' cm**-1')

        msg = b'SetWavelength %.8f\r\n' % wavelength_nm
        try:
            self.pulsare.sendall(msg)
            resp = self.pulsare.recv(self.BUFFER_SIZE)
            print('setWavelength reply: ', resp)
        except:
            resp = 'setWavelength threw an exception'
            print(resp)
        return resp

    def disconnect(self):
        msg = b"RemoteDisconnect\r\n"
        try:
            self.pulsare.sendall(msg)
            resp = self.pulsare.recv(self.BUFFER_SIZE)
            print('RemoteDisconnect reply: ', resp)
        except:
            return (False, 'Not connected to Pulsare!')

        try:
            self.pulsare.close()
            resp = 'pulsare disconnected'
        except:
            resp = 'pulsare not connected'

        return resp

    def getStatus(self):
        msg = b"GetStatus\r\n"
        try:
            self.pulsare.sendall(msg)
            temp_resp = self.pulsare.recv(self.BUFFER_SIZE)
            print('GetStatus reply: ', temp_resp)
        except:
            return (False, 'GetStatus threw an exception')

    def getActualPosition(self):
        msg = b"GetActualPosition\r\n"
        try:
            self.pulsare.sendall(msg)
            temp_resp = self.pulsare.recv(self.BUFFER_SIZE)
            print('GetStatus reply: ', temp_resp)
        except:
            return (False, 'GetActualPosition threw an exception')

if __name__ == '__main__':
    import time

    my_pulsare = pulsare()
    names, values = rf.readpars(my_pulsare.input_filename)

    print(my_pulsare.pulsare_ip, my_pulsare.pulsare_port)
    time.sleep(2)

    my_pulsare.setWavelength(17734.60)
    my_pulsare.getActualPosition()
    # my_pulsare.disconnect()