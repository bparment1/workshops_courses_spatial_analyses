# -*- coding: utf-8 -*-
"""
Spyder Editor.
"""
####################################    Spatial Analyses: SYRACUSE   #######################################
#######################################  Analyse data from Census #######################################
#This script performs basic analyses for the Exercise 1 of the workshop using Census data.
# The overall goal is to explore spatial autocorrelation and aggregation of units of analyses.     
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 12/29/2018 
#DATE MODIFIED: 12/29/2018
#Version: 1
#PROJECT: AAG 2019 workshop preparation
#TO DO:
#
#COMMIT: added Moran'I and spatial regression, AAG workshop
#
##################################################################################################

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
in_dir = "/home/bparmentier/c_drive/Users/bparmentier/Data/python/Exercise_1/data"
#ARGS 2
out_dir = "/home/bparmentier/c_drive/Users/bparmentier/Data/python/Exercise_1/outputs"
#ARGS 3:
create_out_dir=True #create a new ouput dir if TRUE
#ARGS 7
out_suffix = "exercise1_12292018" #output suffix for the files and ouptut folder
#ARGS 8
NA_value = -9999 # number of cores
file_format = ".tif"

ct_2000_fname = "ct_00.shp" # CT_00: Cencus Tracts 2000
bg_2000_fname = "bg_00.shp" # BG_00: Census Blockgroups 2000
bk_2000_fname = "bk_00.shp" # BK_00: Census Blocks 2000

census_table_fname = "census.csv" #contains data from census to be linked
soil_PB_table_fname = "Soil_PB.csv" #same as census table
tgr_shp_fname = "tgr36067lkA.shp" #contains data from census to be linked

metals_table_fname = "SYR_metals.xlsx" #contains metals data to be linked

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

#ct_2000_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(ct_2000_fname)))
#bg_2000_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(bg_2000_fname)))
#bk_2000_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(bk_2000_fname)))

## Counties for Syracuse in 2000
ct_2000_filename = os.path.join(in_dir,ct_2000_fname)
## block groups for Syracuse in 2000
bg_2000_filename = os.path.join(in_dir,bg_2000_fname)
## block for Syracuse in 200
bk_2000_filename = os.path.join(in_dir,bk_2000_fname)

 
ct_2000_gpd = gpd.read_file(ct_2000_filename)
bg_2000_gpd = gpd.read_file(bg_2000_filename)
bk_2000_gpd = gpd.read_file(bk_2000_filename)

ct_2000_gpd.describe()
ct_2000_gpd.plot(column="CNTY_FIPS")
ct_2000_gpd.head()

metals_df = pd.read_excel(os.path.join(in_dir,metals_table_fname))
census_syr_df = pd.read_csv(os.path.join(in_dir,census_table_fname),sep=",",header=0)

#Soil lead samples: UTM z18 coordinates
soil_PB_df = pd.read_csv(os.path.join(in_dir,soil_PB_table_fname),sep=",",header=0) #point locations

census_syr_df.shape #47 spatial entities
ct_2000_gpd.shape #47 spatial entities
metals_df.shape #47 entities
bg_2000_gpd.shape #147 spatial entities

######## PRODUCE MAPS OF 2000 Population #########

#First need to link it.

bg_2000_gdp.columns
census_syr_df.columns
#key is "TRACT" but with a different format.
#First fix the format
bg_2000_gdp.head()
census_syr_df.BKG_KEY.head()

#as.numeric(as.character(ct_2000_sp$TRACT))
ct_2000_gpd.TRACT.dtype
bg_2000_gpd.BKG_KEY.dtypes
census_syr_df.dtypes
census_syr_df.BKG_KEY.dtypes


#bg_2000_gpd['BKG_KEY'].astype(census_syr_df.BKG_KEY.dtypes)
bg_2000_gpd['BKG_KEY']=bg_2000_gpd['BKG_KEY'].astype('int64')
3bg_2000_gpd['BKG_KEY'].astype('|S')

>>> df['purchase'].astype(int)
ct_2000_sp$TRACT <- as.numeric(as.character(ct_2000_sp$TRACT))

bg_2000_sp = merge(bg_2000_sp,census_syr_df,by="BKG_KEY")
bg_2000_gpd = bg_2000_gpd.merge(census_syr_df, on='BKG_KEY')
# country_shapes = country_shapes.merge(country_names, on='iso_a3')

spplot(bg_2000_sp,"POP2000",main="POP2000") #quick visualization of population 

























