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
#DATE MODIFIED: 01/02/2019
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
soil_PB_df = pd.read_csv(os.path.join(in_dir,soil_PB_table_fname),sep=",",header=None) #point locations

census_syr_df.shape #57 spatial entities
ct_2000_gpd.shape #57 spatial entities
metals_df.shape #57 entities
bg_2000_gpd.shape #147 spatial entities

#chekc the crs

######## PRODUCE MAPS OF 2000 Population #########

#First need to link it.

bg_2000_gpd.columns
census_syr_df.columns
#key is "TRACT" but with a different format.
#First fix the format
bg_2000_gpd.head()
census_syr_df.BKG_KEY.head()

#as.numeric(as.character(ct_2000_sp$TRACT))
ct_2000_gpd.TRACT.dtype
bg_2000_gpd.BKG_KEY.dtypes
census_syr_df.dtypes
census_syr_df.BKG_KEY.dtypes

#bg_2000_gpd['BKG_KEY'].astype(census_syr_df.BKG_KEY.dtypes)
bg_2000_gpd['BKG_KEY']=bg_2000_gpd['BKG_KEY'].astype('int64')

#ct_2000_sp$TRACT <- as.numeric(as.character(ct_2000_sp$TRACT))

#bg_2000_sp = merge(bg_2000_sp,census_syr_df,by="BKG_KEY")
bg_2000_gpd = bg_2000_gpd.merge(census_syr_df, on='BKG_KEY')
# country_shapes = country_shapes.merge(country_names, on='iso_a3')

#spplot(bg_2000_sp,"POP2000",main="POP2000") #quick visualization of population 
bg_2000_gpd.plot(column='POP2000',cmap="OrRd",
                 scheme='quantiles')
#plt.title('POP2000')

### Let's use more option with matplotlib
fig, ax = plt.subplots()
bg_2000_gpd.plot(column='POP2000',cmap="OrRd",
                 scheme='quantiles',
                 ax=ax)
ax.set_title('POP2000')
##Now change the classes!

### Summarize by census track
#census_2000_sp <- aggregate(bg_2000_sp , 
#                            by="TRACT",FUN=sum)
census_2000_df = bg_2000_gpd.groupby(['TRACT']).sum()
type(census_2000_df)
#To keep geometry, we must use dissolve method from geopanda
census_2000_gpd = bg_2000_gpd.dissolve(by='TRACT',aggfunc='sum')
type(census_2000_gpd)
## Sum is 57!!!
sum(census_2000_gpd.POP2000==census_2000_df.POP2000)
census_2000_gpd.shape == census_2000_df.shape
##compare to sp!!
#df_test <- aggregate(POP2000 ~ TRACT, bg_2000_sp , 
#                     FUN=sum)

### Check if the new geometry of entities is the same as census

fig, ax = plt.subplots(figsize=(12,8))

# set aspect to equal. This is done automatically
# when using *geopandas* plot on it's own, but not when
# working with pyplot directly.
ax.set_aspect('equal')
census_2000_gpd.plot(ax=ax,column='POP2000',cmap='OrRd')
ct_2000_gpd.plot(ax=ax,color=None,edgecolor="red")

ax.set_title("Population", fontsize= 20)
fig.colorbar(ax) #add palette later
ax.set_axis_off()
plt.show()

#plot(census_2000_sp)
#plot(ct_2000_sp,border="red",add=T)
#nrow(census_2000_sp)==nrow(ct_2000_sp)

df_summary_by_census <- aggregate(. ~ TRACT, bg_2000_sp 
                                  , FUN=sum) #aggregate all variables from the data.frame

##Join by key table id:
dim(ct_2000_sp)
ct_2000_sp <- merge(ct_2000_sp,df_summary_by_census,by="TRACT")
dim(ct_2000_sp)

#save as sp and text table
#write.table(file.path(out_dir,)

### Do another map with different class intervals: 

title_str <- "Population by census tract in 2000"
range(ct_2000_sp$POP2000)
col_palette <- matlab.like(256)

## 9 classes with fixed and constant break
break_seq <- seq(0,9000,1000)
breaks.qt <- classIntervals(ct_2000_sp$POP2000, n=length(break_seq), 
                            style="fixed", fixedBreaks=break_seq, intervalClosure='right')

## Color for each class
#colcode = findColours(breaks.qt , c('darkblue', 'blue', 'lightblue', 'palegreen','yellow','lightpink', 'pink','brown3',"red","darkred"))
p_plot_pop2000_ct <- spplot(ct_2000_sp,
                            "POP2000",
                            col="transparent", #transprent color boundaries for polygons
                            col.regions = col_palette ,
                            main=title_str,
                            at = breaks.qt$brks)
print(p_plot_pop2000_ct)

##### PART II: SPATIAL QUERY #############

## Join metals to census track 
## Join lead (pb) measurements to census tracks

#soil_PB_df <- read.table(file.path(in_dir_var,census_table_fname),sep=",",header=T)
#metals_df= pd.read_xls(os.path.join(in_dir,metals_table_fname))
#metals_df <- read.xls(file.path(in_dir_var,metals_table_fname),sep=",",header=T)

#View(soil_PB_df)
metals_df.head()

##This suggests matching to the following spatial entities
#nrow(metals_df)==nrow(ct_2000_sp)
metals_df.shape[0]== ct_2000_gpd.shape[0]
#nrow(soil_PB_df)==nrow(bg_2000_sp)

#dim(bg_2000_sp)
#census_metals_sp <- merge(ct_2000_sp,metals_df,by.x="TRACT",by.y="ID")
#Check data types before joining tables with "merge"
metals_df.dtypes
ct_2000_gpd.dtypes
ct_2000_gpd['TRACT']=ct_2000_gpd.TRACT.astype('int64')
census_metals_gpd = ct_2000_gpd.merge(metals_df,left_on='TRACT',right_on='ID')

########processing lead data
### Now let's plot lead data 
#Soil lead samples: UTM z18 coordinates
#soil_PB_df <- read.table(file.path(in_dir_var,soil_PB_table_fname),sep=",",header=T) #point locations
census_metals_gpd.crs
#proj4string(census_metals_sp) #
soil_PB_df.columns = ["x","y","ID","ppm"]
#names(soil_PB_df)
soil_PB_df.head()
#names(soil_PB_df) <- c("x","y","ID","ppm") 

soil_PB_gpd = soil_PB_df.copy()
type(soil_PB_df)
soil_PB_gpd['Coordinates']=list(zip(soil_PB_gpd.x,soil_PB_gpd.y))
#coordinates(soil_PB_sp) <- soil_PB_sp[,c("x","y")]
#coordinates(soil_PB_sp) <- soil_PB_sp[,c("x","y")]
type(soil_PB_gpd)
soil_PB_gpd['Coordinates']= soil_PB_gpd.Coordinates.apply(Point)
soil_PB_gpd = gpd.GeoDataFrame(soil_PB_gpd,geometry='Coordinates')

#### Check the coordinates reference system
type(census_metals_gpd.crs) #dictionary
epsg_code = census_metals_gpd.crs.get('init').split(':')[1]

inproj = osr.SpatialReference()
inproj.ImportFromEPSG(int(epsg_code))
inproj.ExportToProj4()

#Assign projection system
soil_PB_gpd.crs= census_metals_gpd.crs
#proj4string(soil_PB_sp) <- proj4string(census_metals_sp)
#dim(soil_PB_sp)

#
#soil_PB_sp <- soil_PB_sp[,c("ID","ppm","x","y")]
#View(soil_PB_sp)
soil_PB_gpd.head()

fig, ax = plt.subplots()

census_metals_gpd.plot(ax=ax,color='white',edgecolor='red')
soil_PB_gpd.plot(ax=ax,marker='*',
                 color='black',
                 markersize=2)
                 
#plot(census_metals_sp)
#plot(soil_PB_sp,add=T)

###### Spatial query: associate points of pb measurements to each census tract
### Get the ID and 
#Use sjoin
test=gpd.tools.sjoin(soil_PB_gpd,census_2000_gpd,how="left")
len(test.BKG_KEY.value_counts()) #associated BKG Key to points
len(test.index_right.value_counts())
test.columns
grouped = test.groupby(['index_right']).mean()
#soil_tract_id_df <- over(soil_PB_sp,census_2000_sp,fn=mean)
#soil_PB_sp <- intersect(soil_PB_sp,census_2000_sp)
#test4 <- gIntersection(soil_PB_sp,census_2000_sp,byid=T)
#names(soil_PB_sp)
soil_PB_sp <- rename(soil_PB_sp, c("d"="TRACT")) #from package plyr

census_pb_avg <- aggregate(ppm ~ TRACT,(soil_PB_sp),FUN=mean)
census_pb_avg <- rename(census_pb_avg,c("ppm"="pb_ppm"))

##Now join
census_metals_pb_sp <- merge(census_metals_sp,census_pb_avg,by="TRACT")
### write out final table and shapefile

outfile<-paste("census_metals_pb_sp","_",
               out_suffix,sep="")
writeOGR(census_metals_pb_sp,dsn= out_dir,layer= outfile, driver="ESRI Shapefile",overwrite_layer=TRUE)

outfile_df_name <- file.path(out_dir,paste0(outfile,".txt"))
write.table(as.data.frame(census_metals_pb_sp),file=outfile_df_name,sep=",")

## For kriging use scipy.interpolate
#https://stackoverflow.com/questions/45175201/how-can-i-interpolate-station-data-with-kriging-in-python


















