# -*- coding: utf-8 -*-
"""
Spyder Editor.
"""
#################################### Land Use and Land Cover Change #######################################
############################ Analyze Land Cover change in Houston #######################################
#This script performs analyses for the Exercise 4 of the AAG Course using aggregated NLCD values.
#The goal is to assess land cover change using two land cover maps in the Houston areas.
#Additional datasets are provided for the land cover change modeling. A model is built for Harris county.
#
#AUTHORS: Benoit Parmentier
#DATE CREATED: 01/07/2019
#DATE MODIFIED: 01/08/2019
#Version: 1
#PROJECT: AAG 2019 Geospatial Short Course
#TO DO:
#
#COMMIT: clean up code for workshop
#
#################################################################################################
	
	
###### Library used in this script

import gdal
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import rasterio
import subprocess
import pandas as pd
import os, glob
from rasterio import plot
import geopandas as gpd
import descartes
import pysal as ps
from cartopy import crs as ccrs
from pyproj import Proj
from osgeo import osr
from shapely.geometry import Point

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
out_suffix = "exercise4_01072018" #output suffix for the files and ouptut folder
#ARGS 8
NA_value = -9999 # number of cores
file_format = ".tif"

#NLCD coordinate reference system: we will use this projection rather than TX.
CRS_reg = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
method_proj_val = "bilinear" # method option for the reprojection and resampling
gdal_installed = True #if TRUE, GDAL is used to generate distance files
		
### Input data files
rastername_county_harris = "harris_county_mask.tif" #Region of interest: extent of Harris County
elevation_fname = "srtm_Houston_area_90m.tif" #SRTM elevation
roads_fname = "r_roads_Harris.tif" #Road count for Harris county
	
### Aggreagate NLCD input files
infile_land_cover_date1 = "agg_3_r_nlcd2001_Houston.tif"
infile_land_cover_date2 = "agg_3_r_nlcd2006_Houston.tif"
infile_land_cover_date3 = "agg_3_r_nlcd2011_Houston.tif"
	
infile_name_nlcd_legend = "nlcd_legend.txt"
infile_name_nlcd_classification_system = "classification_system_nlcd_legend.xlsx"
	
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
	
infile_land_cover_date1 = os.path.join(in_dir,infile_land_cover_date1) #NLCD 2001
infile_land_cover_date2 = os.path.join(in_dir,infile_land_cover_date2) #NLCD 2006
infile_land_cover_date3 = os.path.join(in_dir,infile_land_cover_date3) #NLCD 2011

lc_date1 = rasterio.open(infile_land_cover_date1) 
r_lc_date1 = src.read(1,masked=True) #read first array with masked value, nan are assigned for NA
lc_date2 = rasterio.open(infile_land_cover_date2) 
r_lc_date2 = src.read(1,masked=True) #read first array with masked value, nan are assigned for NA
lc_date3= rasterio.open(infile_land_cover_date2) 
r_lc_date3 = src.read(1,masked=True) #read first array with masked value, nan are assigned for NA

spatial_extent = rasterio.plot.plotting_extent(lc_date1)
plot.show(r_lc_date1)
#plot.show(r_lst,cmap='viridis',scheme='quantiles')
lc_date1.crs # not defined with *.rst
lc_legend_df = pd.read_table(os.path.join(in_dir,infile_name_nlcd_legend),sep=",")
	
lc_legend_df.head() # Inspect data
plot.show(r_lc_date2) # View NLCD 2006, we will need to add the legend use the appropriate palette!!
	
### Let's add legend and examine existing land cover categories
	
freq_tb_date2 <- freq(r_lc_date2)
head(freq_tb_date2) #view first 5 rows, note this is a matrix object.
	
### Let's generate a palette from the NLCD legend information to view the existing land cover for 2006.
names(lc_legend_df)
dim(lc_legend_df) #contains a lot of empty rows
	
lc_legend_df<- subset(lc_legend_df,COUNT>0) #subset the data to remove unsured rows
### Generate a palette color from the input Red, Green and Blue information using RGB encoding:
	
lc_legend_df$rgb <- paste(lc_legend_df$Red,lc_legend_df$Green,lc_legend_df$Blue,sep=",") #combine
	
### row 2 correspond to the "open water" category
color_val_water <- rgb(lc_legend_df$Red[2],lc_legend_df$Green[2],lc_legend_df$Blue[2],maxColorValue = 255)
color_val_developed_high <- rgb(lc_legend_df$Red[7],lc_legend_df$Green[7],lc_legend_df$Blue[7],maxColorValue = 255)
	
lc_col_palette <- c(color_val_water,color_val_developed_high)
	
barplot(c(1,1),
col=lc_col_palette,
main="Visualization of color palette for NLCD land cover",
names.arg=c("Open water", "Developed, High Intensity"),las=1)
	
### Let's generate a color for all the land cover categories by using lapply and function
n_cat <- nrow(lc_legend_df)
lc_col_palette <- lapply(1:n_cat,
FUN=function(i){rgb(lc_legend_df$Red[i],lc_legend_df$Green[i],lc_legend_df$Blue[i],maxColorValue = 255)})
lc_col_palette <- unlist(lc_col_palette)
	
lc_legend_df$palette <- lc_col_palette
	
r_lc_date2 <- ratify(r_lc_date2) # create a raster layer with categorical information
rat <- levels(r_lc_date2)[[1]] #This is a data.frame with the categories present in the raster
	
lc_legend_df_date2 <- subset(lc_legend_df,lc_legend_df$ID%in% (rat[,1])) #find the land cover types present in date 2 (2006)
rat$legend <- lc_legend_df_date2$NLCD.2006.Land.Cover.Class #assign it back in case it is missing
levels(r_lc_date2) <- rat #add the information to the raster layer
	
### Now generate a plot of land cover with the NLCD legend and palette
levelplot(r_lc_date2,
	col.regions = lc_legend_df_date2$palette,
	scales=list(draw=FALSE),
	main = "NLCD 2006")
	


###################### END OF SCRIPT #####################

















