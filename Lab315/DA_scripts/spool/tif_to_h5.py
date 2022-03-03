# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:35:41 2020
@author: T Kierspel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
import sys
import h5py
from PIL import Image
from scipy.ndimage import gaussian_filter1d



plt.close('all')

location = "C:/Users/Gruppe Willitsch/Desktop/data/2020-04-21/"
file     = '010'


filename = location+file+'.tif'

im = Image.open(filename)
im.seek(0)
    
imag   = np.array(im)    
    
images = [imag]
for i in range(1,int(1e6)):
    try:
        im.seek(i)
        images.append(np.array(im))
            
    except:
        print('%i images in this tiff'%i)
        break
    

filename = location+file+'.hdf5'
f = h5py.File(filename, "w")
f.create_dataset("movie", data = images, compression = 'gzip')
f.close()
print("Done")

##                        hf.create_dataset('data%i'%self.count, data=data,compression = 'gzip')
