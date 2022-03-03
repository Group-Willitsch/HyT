"""
Created on Tuesday 28/07/2021
@author: Yanning
Notes:
    - This program includes control functions for the stepper motor (controller: MCM02.1 and I1AM02.1 S1:APS01.1)

To-do:
    - Link to the HyT program
"""

import serial
import time
import numpy as np


class motor:
    def __init__(self):
        self.serial_port = serial.Serial('com12', 115200, timeout=1)
        if self.serial_port.isOpen():
            print('Serial port opened')
            print("Serial port details:", self.serial_port)
        else:
            print('Error: Failed to open serial port')

        self.run_freq = 1100  # default frequency in Hz
        # self.set_run_freq(self.run_freq)

    def send_cmd(self, command):
        cmd = '0' + command.upper() + ':'
        byte_cmd = bytearray(cmd, encoding='utf-8')
        xor_result = 0
        for n in byte_cmd:
            xor_result = int_xor(n, xor_result)
            # print(xor_result)
        cmd = chr(2) + cmd + str(format(xor_result, '02X')) + chr(3)
        print("Sent>>>:", cmd)
        self.serial_port.write(cmd.encode())
        time.sleep(0.25)

    def read_response(self):
        try:
            response = self.serial_port.readline().decode()
            if response[1] == chr(6):
                print("Received<<<:", response[2:-4])
                return response[2:-4]
            else:
                print("Warning: unexpected response")
        except:
            return False, 'Error: Failed to read response'

    def rotate_pos(self, num_steps):
        """
        Rotate the motor towards OH side
        Note:
            - After sending the rotating command, it will sleep for the rotating time to prevent further command.
        @param num_steps: number of steps
        """
        command = "1.1+" + str(num_steps)
        self.send_cmd(command)
        rotating_time = num_steps/self.run_freq
        time.sleep(rotating_time)

    def rotate_neg(self, num_steps):
        """
        Rotate the motor towards ion side
        Note:
            - After sending the rotating command, it will sleep for the rotating time to prevent further command.
        @param num_steps: number of steps
        """
        command = '1.1-' + str(num_steps)
        self.send_cmd(command)
        rotating_time = num_steps/self.run_freq
        time.sleep(rotating_time)

    def dwell(self, dwell_time):
        """
        Stop the motor for a duration given by dwell_time in microseconds.
        """
        print('Dwell for', dwell_time * 1e-3, ' ms')
        time.sleep(dwell_time / 1000000)

    def set_run_freq(self, run_freq):
        """
        Set the motor run frequency in Hz (see the manual Principles of Positioning for Stepper Motor Controllers for
        more info about relation between rpm and Hz, page 7)
        """
        self.run_freq = run_freq
        command = "1.1P14S" + str(self.run_freq)
        self.send_cmd(command)
        # print("Frequency is set to: ", run_freq, "Hz")

    def read_run_freq(self):
        command = "1.1P14R"
        self.send_cmd(command)
        print("The run frequency is:", self.read_response(), " Hz")
        return self.read_response()

    def close_serial_port(self):
        try:
            self.serial_port.close()
            print('Serial port closed')
        except:
            print('Error: Failed to close serial port')


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

    serial_com_test = False

    if serial_com_test:
        my_motor.send_cmd('1.1P04R')
        my_motor.read_response()
        my_motor.send_cmd('1.1P04S400')
        my_motor.read_response()
        my_motor.send_cmd('1.1P04R')
        my_motor.read_response()

    shuttle_test = True

    if shuttle_test:
        num_steps = 3000
        dwell_time = 10e6  # in us
        my_motor.dwell(dwell_time)
        my_motor.rotate_neg(num_steps)
        dwell_time = 10e6  # in us
        my_motor.dwell(dwell_time)
        my_motor.rotate_pos(num_steps)

    my_motor.close_serial_port()
