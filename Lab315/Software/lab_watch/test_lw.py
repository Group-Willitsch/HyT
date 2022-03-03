import numpy as np
from lab_watch import send_email
import time

flag_email = False   ####used for pressure email

for i in range(10):
    if i > 1:
        if not flag_email:
            send_email('The pressure is too low %i'%i)
            flag_email = True
            t_send = time.time()
            print('YES')
        else:
            t_act  = time.time()
            if(t_act-t_send)>3:
                flag_email = False


    time.sleep(1)
