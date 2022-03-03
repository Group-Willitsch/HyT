"""
Created on Tuesday 28/07/2021
@author: Yanning
Notes:
    - This program includes control functions for the Phytron stepper motor (controller: MCM02.1 and I1AM02.1 S1:APS01.1)

To-do:
    -
"""

import serial
import time
import numpy as np


class motor:
    def __init__(self):
        self.run_freq = 8000  # default frequency in Hz
        self.acc_ramp = 64000 # default acc/dec ramp in Hz/s
        self.position = 0  # in steps, always the position when motor object is initialized
        self.isConnected = False

        open_port_when_init = False
        if open_port_when_init:
            self.open_serial_port()
            # self.set_run_freq(self.run_freq)

    def open_serial_port(self):
        try:
            self.serial_port = serial.Serial('com12', 115200, timeout=1)
        except serial.SerialException as err:
            print('Error:', err.args)
            raise Exception(err.args[0])
        else:
            if self.serial_port.isOpen():
                print("Serial port opened:", self.serial_port)
                self.isConnected = True
            else:
                print('Error: Serial port open failed')

    def send_cmd(self, command):
        cmd = '0' + command.upper() + ':'
        byte_cmd = bytearray(cmd, encoding='utf-8')
        xor_result = 0
        for n in byte_cmd:
            xor_result = int_xor(n, xor_result)
            # print(xor_result)
        cmd = chr(2) + cmd + str(format(xor_result, '02X')) + chr(3)
        try:
            self.serial_port.write(cmd.encode())
            # time.sleep(0.25)
            # print("Sent>>>:", cmd)
        except:
            raise Exception("Error: Serial port failed to write, please check connection")

        return

    def read_response(self):
        try:
            response = self.serial_port.readline().decode()
            if True:#response[1] == chr(6):
                print("Received<<<", response[2:-4])
                return response[2:-4]
            else:
                print("Warning: Unexpected response")
                return False
        except:
            raise Exception('Error: Serial port failed to read, please check connection')

    def rotate_pos(self, num_steps):
        """
        Rotate the motor towards OH side
        Note:
            - After sending the rotating command, it will sleep for the rotating time to prevent further command.
        @param num_steps: number of steps
        """
        command = "1.1+" + str(int(num_steps))
        self.send_cmd(command)
        # rotating_time = num_steps / self.run_freq
        # time.sleep(rotating_time*1)
        self.position += int(num_steps)
        return

    def rotate_neg(self, num_steps):
        """
        Rotate the motor towards ion side
        Note:
            - After sending the rotating command, it will sleep for the rotating time to prevent further command.
        @param num_steps: number of steps
        """
        command = '1.1-' + str(int(num_steps))
        self.send_cmd(command)
        # rotating_time = num_steps / self.run_freq
        # time.sleep(rotating_time*1)
        self.position -= int(num_steps)
        return

    def rotate(self, signed_num_steps):
        """
        Rotate the motor according to the signed number of steps
        Note:
            - After sending the rotating command, it will sleep for the rotating time to prevent further command.
        @param signed_num_steps: number of steps
        """
        signed_num_steps = int(signed_num_steps)
        if signed_num_steps == 0:
            return
        else:
            if signed_num_steps > 0:
                self.rotate_pos(signed_num_steps)
            else:
                self.rotate_neg(-signed_num_steps)
        return

    def stop_rotation(self):
        self.send_cmd('1.1S')

    def dwell(self, dwell_time):
        """
        Stop the motor for a duration given by dwell_time in microseconds.
        """
        print('Dwell for', dwell_time, ' s')
        time.sleep(dwell_time)

    def rotation_time(self, num_steps, tolerance=0.):
        '''
        Estimate rotation time (one-way) according to number of steps
        @param tolerance: tolerance time added to the estimated value
        '''
        return self.run_freq/self.acc_ramp + int(num_steps)/self.run_freq + tolerance


    def set_run_freq(self, run_freq):
        """
        Set the motor run frequency in Hz (see the manual Principles of Positioning for Stepper Motor Controllers for
        more info about relation between rpm and Hz, page 7)
        """
        self.run_freq = int(run_freq)
        command = "1.1P14S" + str(self.run_freq)
        self.send_cmd(command)
        return

    def set_acc_ramp(self, acc_ramp):
        """
        Set the motor acc/dec ramp in Hz/s (see the manual Principles of Positioning for Stepper Motor Controllers for
        more info, page 10)
        """
        self.acc_ramp = int(acc_ramp)
        command = "1.1P15S" + str(self.acc_ramp)
        self.send_cmd(command)
        return

    def read_run_freq(self):
        self.serial_port.flush()
        self.send_cmd("1.1P14R")
        response = self.read_response()
        if response:
            # print("The readout run frequency is:", response, "Hz")
            return response
        else:
            return False

    def read_acc_ramp(self):
        self.serial_port.flush()
        self.send_cmd("1.1P15R")
        response = self.read_response()
        if response:
            # print("The readout acceleration ramp is:", response, "Hz/s")
            return response
        else:
            return False

    def close_serial_port(self):
        try:
            self.serial_port.close()
            print('Serial port closed')
            self.isConnected = False
        except:
            print('Error: Failed to close serial port')
            raise Exception('Error: Failed to close serial port')

    def shuttle_total_time(self, oneway_steps, stopover_time, delay_of_start = 0., repetitions = 1):
        return (self.rotation_time(num_steps=oneway_steps)*2+stopover_time+delay_of_start)*repetitions

def int_xor(n, m):
    """
    Calculate the logic "xor" of two decimal integers by converting them to binaries, doing xor operation and
    then converting back.
    """
    bool_list_xor = np.logical_xor([int(x) for x in format(n, '08b')], [int(x) for x in format(m, '08b')])
    t = ''.join([str(int(b)) for b in bool_list_xor])
    return int(t, 2)


if __name__ == "__main__":
    my_motor = motor()
    my_motor.open_serial_port()

    serial_com_test = True

    if serial_com_test:
        my_motor.send_cmd('1.1P14R')
        my_motor.read_response()
        my_motor.send_cmd('1.1P15R') ##read the acc/dec ramp
        my_motor.read_response()

        my_motor.read_run_freq()
        my_motor.read_acc_ramp()
        # my_motor.send_cmd('1.1P04S400')
        # my_motor.read_response()
        # my_motor.send_cmd('1.1P04R')
        # my_motor.read_response()

    shuttle_test = False

    if shuttle_test:
        num_steps = 3000
        dwell_time = 6  # in s
        my_motor.dwell(dwell_time)
        my_motor.rotate_neg(num_steps)
        dwell_time = 3  # in s
        my_motor.dwell(dwell_time)
        my_motor.rotate_pos(num_steps)

    my_motor.close_serial_port()
