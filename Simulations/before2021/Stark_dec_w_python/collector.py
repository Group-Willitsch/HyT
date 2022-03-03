import numpy as np


data = 0

for i in range(1,7):
    if i == 1:
        data = np.load('part%s.npy'%i)
    else:
        data = np.append(data,np.load('part%s.npy'%i),axis=1)
