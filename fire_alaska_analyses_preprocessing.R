####################################    Fire Alaska Analyses   #######################################
############################  Analyse using time series and fire locations  #######################################
#This script prepares data for classes and workshops.
# The overall goal is to explore the impact of fire on surface state variables observed from Remote Sensing.     
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/09/2017 
#DATE MODIFIED: 03/11/2017
#Version: 2
#PROJECT: AAG 2017 workshop preparation
#TO DO:
#
#COMMIT: masking data and transforming to tiff
#
#################################################################################################

###Loading R library and packages                                                      

library(sp)
library(rgdal)
library(spdep)
library(gtools)
library(maptools)
library(parallel)
library(rasterVis)
library(raster)
library(forecast) #ARIMA forecasting
library(xts)
library(zoo)
library(lubridate)
library(colorRamps) #contains matlab.like color palette
library(rgeos)
library(sphet) #contains spreg
library(BMS) #contains hex2bin and bin2hex
library(bitops)

###### Functions used in this script

function_preprocessing_and_analyses <- "fire_alaska_analyses_preprocessing_functions_03102017.R" #PARAM 1
script_path <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/R_scripts"
source(file.path(script_path,function_preprocessing_and_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

in_dir_NDVI <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Fire_Alaska_analyses/NDVI_Alaska"

out_dir <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Fire_Alaska_analyses/outputs"

mask_alaska_filename <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Fire_Alaska_analyses/Lab_time_series/Data_unprocessed/mask_Alaska_11112014.rst" 

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" #CONST 1
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4


file_format <- ".rst" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"alaska_03092017" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######
#

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


r_mask <- raster(mask_alaska_filename)
NAvalue(r_mask) <- -999
plot(r_mask)

lf <- mixedsort(list.files(path=in_dir_NDVI,pattern="NDVI_2001_2009_filling6__NDVI_gap_filling2_07062011_test1_.*.rst$",full.names=T))

r_NDVI_data <- stack(lf)

plot(r_NDVI_data,y=1)

r_NDVI_data_mask <- mask(r_NDVI_data,mask=r_mask,filename=file.path(out_dir,"r_NDVI_data_alaska_mask2.tif"))

plot(r_NDVI_data_mask,y=2)
filename(r_NDVI_data_mask) 

################### End of Script #########################
