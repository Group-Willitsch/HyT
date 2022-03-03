import numpy as np


def send_cmd(cmd):
    cmd = '0'+cmd.upper()+':'
    byte_cmd = bytearray(cmd,encoding = 'utf-8')
    xor_result = 0
    for n in byte_cmd:
        xor_result = int_xor(n,xor_result)
        # print(xor_result)
    cmd = chr(2)+cmd+str(format(xor_result,'02X'))+chr(3)
    print(cmd)


def int_xor(n, m):
    """
    Calculate the logic "xor" of two decimal integers by converting them to binaries, doing xor operation and
    then converting back.
    """
    bool_list_xor=np.logical_xor([int(x) for x in format(n,'08b')], [int(x) for x in format(m, '08b')])
    t=''.join([str(int(b)) for b in bool_list_xor])
    return int(t,2)

cmd='1.1+1000'
send_cmd(cmd)
# print(int_xor(58,24))