import numpy as np
from numba import jit
import time

@jit(nopython = True)
def multiply(x,y):
    #z = x*y
    z = np.empty(np.shape(x))
    for i in range(0,len(x)):
        z[i] = x[i]*y[i]
    return(z)


@jit(nopython=True)
def indexing_2(x,y,v):
    idx_x = (np.abs(x - v)).argmin()
    idx_y = (np.abs(y - v)).argmin()
    return(idx_x,idx_y)

@jit(nopython=True)
def indexing(array,value):
    index = 0
    diff  = 1e3
    for i in range(len(array)):
        diff_new = array[i]-value
        if abs(diff_new) < abs(diff):
            diff = diff_new
            index = i
    return(index)



@jit(nopython=True)
def indexing_for_interpolation(array,value):
    index = 0
    diff  = 1e3
    for i in range(len(array)):
        diff_new = array[i]-value
        if abs(diff_new) < abs(diff):
            diff = diff_new
            index = i
    return(index)



N = int(1e2)
x = np.random.random_sample(N)
y = np.random.random_sample(N)

v = 0.5

end = 1000000


###--- Test 1: Multiplication---:
test_multiply = False
###--- Test 2: Indexing---:
test_index    = True


if test_index:

    t1 = time.time()
    for i in range(end):
        idx_x = (np.abs(x - v)).argmin()
        idx_y = (np.abs(y - v)).argmin()
    t2 = time.time()
    print(str(t2-t1))

    t1 = time.time()
    for i in range(end):
        idx_x,idx_y = indexing_2(x,y,v)
    t2 = time.time()
    print(str(t2-t1))

    t1 = time.time()
    for i in range(end):
        idx_x_1 = indexing(x,v)
        idx_y_1 = indexing(y,v)
    t2 = time.time()
    print(str(t2-t1))

    @jit()
    def bla(x,y,v,end):
        for i in range(end):
            idx_x_1 = indexing(x,v)
            idx_y_1 = indexing(y,v)
        return(idx_x_1,idx_y_1)


    t1 = time.time()
    x_t,y_t = bla(x,y,v,end)
    t2 = time.time()
    print(str(t2-t1))


if test_multiply:
    print('\n')

    t1 = time.time()
    for i in range(0,end):
        z = x*y
    t2 = time.time()
    print(str(t2-t1))

    t1 = time.time()
    for i in range(0,end):
        z = multiply(x,y)

    t2 = time.time()
    print(str(t2-t1))

    t1 = time.time()
    for i in range(0,end):
        z = np.multiply(x,y)
    t2 = time.time()
    print(str(t2-t1))
