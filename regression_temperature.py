# -*- coding: utf-8 -*-
"""
Spyder Editor.
"""
#################################### Regression Temperature #######################################
######################## Analyze and predict air temperature with Earth Observation data #######################################
#This script performs analyses to predict air temperature using several coveriates.
#The goal is to predict air temperature using Remotely Sensing data as well as compare measurements
# from the ground station to the remotely sensed measurements.
#
#AUTHORS: Benoit Parmentier
#DATE CREATED: 09/07/2018
#DATE MODIFIED: 03/07/2019
#Version: 1
#PROJECT: SESYNC Geospatial Course and AAG 2019 Python Geospatial Course
#TO DO:
#
#COMMIT: clean up code for workshop
#
#################################################################################################

###### Library used in this script
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import seaborn as sns
import rasterio
import subprocess
import pandas as pd
import os, glob
from rasterio import plot
import geopandas as gpd
import georasters as gr
import gdal
import rasterio
import descartes
import pysal as ps
from cartopy import crs as ccrs
from pyproj import Proj
from osgeo import osr
from shapely.geometry import Point
from collections import OrderedDict
import webcolors
import sklearn

################ NOW FUNCTIONS  ###################

##------------------
# Functions used in the script
##------------------

def create_dir_and_check_existence(path):
    #Create a new directory
    try:
        os.makedirs(path)
    except:
        print ("directory already exists")

############################################################################
#####  Parameters and argument set up ###########

#ARGS 1
in_dir = "/home/bparmentier/c_drive/Users/bparmentier/Data/python/climate_regression/data/Oregon_covariates"
#ARGS 2
out_dir = "/home/bparmentier/c_drive/Users/bparmentier/Data/python/climate_regression/outputs"

#in_dir="/nfs/bparmentier-data/Data/workshop_spatial/climate_regression/data/Oregon_covariates"
#out_dir="/nfs/bparmentier-data/Data/workshop_spatial/climate_regression/outputs"

#ARGS 3:
create_out_dir=True #create a new ouput dir if TRUE
#ARGS 7
out_suffix = "exercise4_03032019" #output suffix for the files and ouptut folder
#ARGS 8
NA_value = -9999 # NA flag balue
file_format = ".tif"

#NLCD coordinate reference system: we will use this projection rather than TX.
CRS_reg = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
method_proj_val = "bilinear" # method option for the reprojection and resampling
gdal_installed = True #if TRUE, GDAL is used to generate distance files


#epsg 2991
crs_reg = "+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

#infile = "mean_month1_rescaled.rst" # mean LST for January
infile = "lst_mean_month1_rescaled.tif" 
infile_forest_perc =""
ghcn_filename = "ghcn_or_tmax_covariates_06262012_OR83M.shp" # climate stations

prop = 0.3
random_seed= 100

################# START SCRIPT ###############################

######### PART 0: Set up the output dir ################

#set up the working directory
#Create output directory

if create_out_dir==True:
    #out_path<-"/data/project/layers/commons/data_workflow/output_data"
    out_dir = "output_data_"+out_suffix
    out_dir = os.path.join(in_dir,out_dir)
    create_dir_and_check_existence(out_dir)
    os.chdir(out_dir)        #set working directory
else:
    os.chdir(create_out_dir) #use working dir defined earlier

###########################################
### PART I: READ AND VISUALIZE DATA #######

data_gpd = gpd.read_file(os.path.join(in_dir,ghcn_filename))

data_gpd.head()  

## Extracting information from raster using raster io object
lst = rasterio.open(os.path.join(in_dir,infile))
lst.crs # explore Coordinate Reference System 
lst.shape
lst.height
plot.show(lst)

## Read raster bands directly to Numpy arrays and visualize data
r_lst = src.read(1,masked=True) #read first array with masked value, nan are assigned for NA
spatial_extent = rasterio.plot.plotting_extent(src)
type(r_lst)
r_lst.size
plot.show(r_lst)
#Can also use the regular matplotlib library function to plot images
plt.imshow(r_lst)
plt.imshow(r_lst, clim=(259.0, 287.0))

# Explore values distribution
plt.hist(r_lst.ravel(),bins=256,range=(259.0,287.0))

##### Combine raster layer and geogpanda layer

data_gpd.plot(marker="*",color="green",markersize=5)
station_or = data_gpd.to_crs({'init': 'epsg:2991'}) #reproject to  match the  raster image

##### How to combine plots with rasterio package
fig, ax = plt.subplots()
with rasterio.open(os.path.join(in_dir,infile)) as src:
        rasterio.plot.show((src,1),ax=ax,
                          clim=(259.0,287.0),)
station_or.plot(ax=ax,marker="*",
              color="red",
               markersize=10)
               
##### How to combine plots with matplotlib package
fig, ax = plt.subplots(figsize = (8,3))
lst_plot = ax.imshow(r_lst, 
                       cmap='Greys', 
                       extent=spatial_extent)
ax.set_title("Long term mean for January land surface temperature", fontsize= 20)
fig.colorbar(lst_plot)
# turn off the x and y axes for prettier plotting
#ax.set_axis_off(); #this removes coordinates on the plot

###########################################
### PART II : Extract information from raster and prepare covariates #######
#raster = './data/slope.tif'

data=gr.from_file(os.path.join(in_dir,infile))
type(data) # check that we have a georaster object
# Plot data
data.plot()
data.plot(clim=(259.0, 287.0))

#### Extract information from raster using coordinates
x_coord = station_or.geometry.x # pands.core.series.Series
y_coord = station_or.geometry.y
# Find value at point (x,y) or at vectors (X,Y)
values = data.map_pixel(x_coord,y_coord)
station_or.columns #get names of col

station_or['year'].value_counts()
station_or.groupby(['month'])['value'].mean()
     
print("number of rows:",station_or.station.count(),"number of stations:",len(station_or.station.unique()))
station_or['LST1'] = values - 273.15 #create new column

station_or_jan = station_or.loc[(station_or['month']==1) & (station_or['value']!=-9999)]
station_or_jan.head()
station_or_jan.columns

#avg_df = station_or.groupby(['station'])['value'].mean())
avg_df = station_or_jan.groupby(['station'])['value','LST1'].mean()
avg_df['value']= avg_df['value']/10
avg_df.head()
         
################################################
###  PART III : Fit model and generate prediction

### Add split training and testing!!!
### Add additionl covariates!!

#selected_covariates_names_updated = selected_continuous_var_names + names_cat 
selected_covariates_names_updated = ['LST1'] 
selected_target_names = ['value']
## Split training and testing

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(avg_df[selected_covariates_names_updated], 
                                                    avg_df[selected_target_names], 
                                                    test_size=prop, 
                                                    random_state=random_seed)

X_train.shape

from sklearn.linear_model import LinearRegression
regr = LinearRegression().fit(X_train,y_train)

#regr = linear_model.LinearRegression()
regr.fit(x, y)

plt.scatter(X_train, y_train,  color='black')
plt.plot(X_train, regr.predict(X_train), color='blue', linewidth=3)
plt.xticks(())
plt.yticks(())
plt.show()

print('reg coef',regr.coef_)
print('reg intercept',regr.intercept_)

#reg.predict(x) # Note this is a fit!
#reg.score(x, y)

regr.predict(X_train) # Note this is a fit!
regr.score(X_train, y_train)

X_test.shape
regr.predict(X_test) # Note this is a fit!
regr.score(X_test, y_test)

## As the plot shows for 2006, we have 15 land cover types. Analyzing such complex categories in terms of decreasse (loss), increase (gain),
### Do models for January,July with LST and with/without land cover % of forest
## Calculate MAE,RMSE,R2,etc. inspire yourself from paper. Save this into a CSV file.

############################# END OF SCRIPT ###################################

















