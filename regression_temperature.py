# -*- coding: utf-8 -*-
"""
Spyder Editor.
"""
#################################### Regression Temperature #######################################
############################ Analyze and predict air temperature with Earth Observation data #######################################
#This script performs analyses for the Exercise 4 of the AAG Course using aggregated NLCD values.
#The goal is to assess land cover change using two land cover maps in the Houston areas.
#Additional datasets are provided for the land cover change modeling. A model is built for Harris county.
#
#AUTHORS: Benoit Parmentier
#DATE CREATED: 01/07/2019
#DATE MODIFIED: 02/16/2019
#Version: 1
#PROJECT: AAG 2019 Geospatial Short Course
#TO DO:
#
#COMMIT: clean up code for workshop
#
#################################################################################################


###### Library used in this script
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import rasterio
import subprocess
import pandas as pd
import os, glob
from rasterio import plot
import geopandas as gpd
import georasters as gr



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

def open_image(url):
    image_data = open_http_query(url)

    if not image_data:
            return None

    mmap_name = "/vsimem/"+uuid4().get_hex()
    gdal.FileFromMemBuffer(mmap_name, image_data.read())
    gdal_dataset = gdal.Open(mmap_name)
    image = gdal_dataset.GetRasterBand(1).ReadAsArray()
    gdal_dataset = None
    gdal.Unlink(mmap_name)

    return image

############################################################################
#####  Parameters and argument set up ###########

#ARGS 1
in_dir = "/home/bparmentier/c_drive/Users/bparmentier/Data/python/Exercise_4/data"
#ARGS 2
out_dir = "/home/bparmentier/c_drive/Users/bparmentier/Data/python/Exercise_4/outputs"
#ARGS 3:
create_out_dir=True #create a new ouput dir if TRUE
#ARGS 7
out_suffix = "exercise4_02162018" #output suffix for the files and ouptut folder
#ARGS 8
NA_value = -9999 # number of cores
file_format = ".tif"

#NLCD coordinate reference system: we will use this projection rather than TX.
CRS_reg = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
method_proj_val = "bilinear" # method option for the reprojection and resampling
gdal_installed = True #if TRUE, GDAL is used to generate distance files

in_dir="/nfs/bparmentier-data/Data/workshop_spatial/climate_regression/data/Oregon_covariates"
out_dir="/nfs/bparmentier-data/Data/workshop_spatial/climate_regression/outputs"

#epsg 2991
crs_reg = "+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=400000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

infile = "mean_month1_rescaled.rst"
infile_forest_perc =""
#infile = "mean_month1_rescaled.tif"
#infile = "lst_mean_month1_rescaled.tif"

ghcn_filename = "ghcn_or_tmax_covariates_06262012_OR83M.shp"


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


#######################################
### PART 1: Read in DATA #######

###########################################
### PART I: READ AND VISUALIZE DATA #######

data_gpd = gpd.read_file(os.path.join(in_dir,ghcn_filename))

data_gpd 

# -12 layers from land cover concensus (Jetz lab)
fileglob = "*.rst"
pathglob = os.path.join(in_dir, fileglob)
l_f = glob.glob(pathglob)
l_f.sort() #order input by decade
l_dir = map(lambda x: os.path.splitext(x)[0],l_f) #remmove extension
l_dir = map(lambda x: os.path.join(out_dir,os.path.basename(x)),l_dir) #set the directory output
 

# Read raster bands directly to Numpy arrays.
with rasterio.open(os.path.join(in_dir,infile)) as src:
        r_lst = src.read(1,masked=True) #read first array with masked value, nan are assigned for NA
        spatial_extent = rasterio.plot.plotting_extent(src)

plot.show(r_lst)
#plot.show(r_lst,cmap='viridis',scheme='quantiles')

src.crs # not defined with *.rst
#Note that you can also plot the raster io data reader
type(r_lst)

r_lst.size
#r_lst.ndim #array dimension
src.height
#src.profile
type(r_lst)

#Can also use the regular matplotlib library function to plot images
#plt.imshow(subset)
plt.imshow(r_lst)
#plt.hist(r_lst)

#see: https://matplotlib.org/users/image_tutorial.html
plt.imshow(r_lst, clim=(259.0, 287.0))
plt.hist(r_lst.ravel(),bins=256,range=(259.0,287.0))

data_gpd.plot(marker="*",color="green",markersize=5)
station_or = data_gpd.to_crs({'init': 'epsg:2991'})

#https://www.earthdatascience.org/courses/earth-analytics-python/lidar-raster-data/customize-matplotlib-raster-maps/

fig, ax = plt.subplots()
with rasterio.open(os.path.join(in_dir,infile)) as src:
        rasterio.plot.show((src,1),ax=ax,
                          clim=(259.0,287.0),)

#plot.show(r_lst, clim=(259.0, 287.0),ax=ax)
#with rasterio.plot.show((src,1),ax=ax)
station_or.plot(ax=ax,marker="*",
              color="red",
               markersize=10)

fig, ax = plt.subplots(figsize = (8,3))
lst_plot = ax.imshow(r_lst, 
                       cmap='Greys', 
                       extent=spatial_extent)
ax.set_title("Long term mean for January land surface temperature", fontsize= 20)
fig.colorbar(lst_plot)
# turn off the x and y axes for prettier plotting
#ax.set_axis_off(); #this removes coordinates on the plot

#raster = './data/slope.tif'
data=gr.from_file(os.path.join(in_dir,infile))
# Plot data
data.plot()

x_coord = station_or.geometry.x # pands.core.series.Series
y_coord = station_or.geometry.y

# Find value at point (x,y) or at vectors (X,Y)
values = data.map_pixel(x_coord,y_coord)
list(station_or) #get names of col
station_or['year'].value_counts()
station_or.groupby(['month'])['value'].mean()
     
print("number of rows:",station_or.station.count(),"number of stations:",len(station_or.station.unique()))
station_or['LST1'] = value-273.15
station_or_jan = station_or.loc[(station_or['month']==1) & (station_or['value']!=-9999)]
station_or_jan.head()
#avg_df = station_or.groupby(['station'])['value'].mean())
avg_df = station_or_jan.groupby(['station'])['value','LST1'].mean()
avg_df['value']= avg_df['value']/10
avg_df.head()
         
#######
################################################
###  PART II : Analyze change and transitions

from sklearn.linear_model import LinearRegression
x=avg_df.LST1.values
y=avg_df.value.values
x = x.reshape(len(x), 1)
y = y.reshape(len(y), 1)
regr = LinearRegression().fit(x,y)

#regr = linear_model.LinearRegression()
regr.fit(x, y)

# plot it as in the example at http://scikit-learn.org/
plt.scatter(x, y,  color='black')
plt.plot(x, regr.predict(x), color='blue', linewidth=3)
plt.xticks(())
plt.yticks(())
plt.show()

print('reg coef',regr.coef_)
print('reg intercept',regr.intercept_)

reg.predict(x) # Note this is a fit!
reg.score(x, y)
## As the plot shows for 2006, we have 15 land cover types. Analyzing such complex categories in terms of decreasse (loss), increase (gain),

###################### END OF SCRIPT #####################

















