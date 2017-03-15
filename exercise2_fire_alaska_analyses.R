####################################    Fire Alaska Analyses   #######################################
############################  Analyse using time series and fire locations  #######################################
#This script performs analyzes for the Exercise 2 of the workshop using NDVI data.
# The overall goal is to explore the impact of fire on surface state variables observed from Remote Sensing.     
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/15/2017 
#DATE MODIFIED: 03/15/2017
#Version: 1
#PROJECT: AAG 2017 workshop preparation
#TO DO:
#
#COMMIT: initial commit for exercise 2, AAG workshop
#
#################################################################################################

###Loading R library and packages                                                      

library(sp) # spatial/geographfic objects and functions
library(rgdal) #GDAL/OGR binding for R with functionalities
library(spdep) #spatial analyses operations, functions etc.
library(gtools) # contains mixsort and other useful functions
library(maptools) # 
library(parallel) # parallel computation, part of base package no
library(rasterVis) # raster visualization operations
library(raster) # raster functionalities
library(forecast) #ARIMA forecasting
library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(lubridate) # dates functionality
library(colorRamps) #contains matlab.like color palette
library(rgeos) #contains topological operations
library(sphet) #contains spreg, spatial regression modeling
library(BMS) #contains hex2bin and bin2hex
library(bitops) #

###### Functions used in this script

function_preprocessing_and_analyses <- "fire_alaska_analyses_preprocessing_functions_03102017.R" #PARAM 1
script_path <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/R_scripts"
source(file.path(script_path,function_preprocessing_and_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

in_dir_NDVI <- "~/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Fire_Alaska_analyses/NDVI_alaska_2005"
in_dir_var <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Fire_Alaska_analyses/data_alaska_analayses"
out_dir <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Fire_Alaska_analyses/outputs"

CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"exercise_03152017" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######


## First create an output directory

if(is.null(out_dir)){
  out_dir <- dirname(in_dir) #output will be created in the input dir
}

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix_s)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

## Second list.files and create raster images stack

lf_NDVI <- list.files(path=in_dir_NDVI, pattern="*.tif",full.names=T)
r_NDVI_ts <- stack(lf_NDVI)

plot(r_NDVI_ts)
date_range <- c("2005.01.01","2005.12.31") #NDVI Alaska
  
#generate dates for 16 days product
dates_val <- generate_dates_by_step(date_range[1],date_range[2],16)$dates #NDVI Katrina

file.info(lf_NDVI[1])$size/(1024*1024) # this is in bytes, convert to mb
dataType(r_NDVI_ts)
inMemory(r_NDVI_ts)
dim(r_NDVI_ts)

##### Now examine other files from Alaska

lf_var <- list.files(path=in_dir_var,pattern="*.tif$",full.names=T)
r_var <- stack(lf_var)
dim(r_var)
plot(r_var)

#filename(r_wwf)
#plot(r_wwf)

##### 
tb_freq <- freq(subset(r_var,6)) #  count of pixels by burn scars
projection(r_var)

####### PART II: Change analyses via image differencing and ratioing.

## Studying change by image differencing
#look for NDVI 2002 and 2009 and select each average

r_NDVI_avg_2002 <- subset(r_var,4)
r_NDVI_avg_2009 <- subset(r_var,5)

r_diff_NDVI <- r_NDVI_avg_2009 - r_NDVI_avg_2002
r_ratio_NDVI <- r_NDVI_avg_2009/r_NDVI_avg_2002
  
plot(r_diff_NDVI,col=matlab.like(255))
plot(r_diff_NDVI,col=matlab.like(255),zlim=c(-0.5,0.5))

plot(r_ratio_NDVI,col=matlab.like(255),zlim=c(0.5,1.5))

### Quick histogram
hist(r_diff_NDVI)
##Adjust the bins

hist_bins <- seq(-1,1,by=0.05)
hist(r_diff_NDVI,breaks=hist_bins)
#hist_bins <- seq(-0.3,0.3,by=0.05)
hist(r_diff_NDVI,breaks=hist_bins,xlim=c(-0.3,0.3))

## Threshold your inputs
plot(r_diff_NDVI < -0.2)
plot(r_diff_NDVI < -0.1)

##Standardize images

##Extract values by fire scars

## Plot changes over the time period

### Your turn: use albedo images to define image of changes, what do you think is a good threshold

#1.Select Albedo images
#2.Do differencing, and standardization
#3.Generate Image of cahnge

######### PART III: time series analyses

#1. Extract time series
#2. Generate movies through looping
#3. Compute temporal corr
#4. Compute autocorr
#5. Performa PCA
#4. compute averages by WWF ecoregions

# Using LST time series perform similar analysis

################### End of Script #########################

#To explore later:
#http://stackoverflow.com/questions/12670972/doing-pca-on-very-large-data-set-in-r
