# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import rasterio as rio
from osgeo import gdal

# parameters and arguments

in_dir = "/home/bparmentier/c_drive/Users/bparmentier/Data/India_research/Data"
out_dir ="/home/bparmentier/c_drive/Users/bparmentier/Data/India_research/outputs"
arr = np.arange(10)
in_file_tellus = "GRCTellus.CSR.200204_201701.LND.RL05.DSTvSCS1409.nc"
x=arr
y= arr +2
plt.plot(x,y)
arr = np.zeros((3,3,3)) #this takes in a tuple
type(arr)
arr.size
arr.shape # this is 3 matrices of 3x3
mat=arr.matrix()
arr=np.zeros((3,3))
mat=arr.matrix()

r=rio.open(os.path.join(in_dir,in_file_tellus))
r_np=r.read()
r_np.shape #this is a multiband file with 159 dates
r_np=r.read(1)
r_np.shape
plt.imshow(r_np)





