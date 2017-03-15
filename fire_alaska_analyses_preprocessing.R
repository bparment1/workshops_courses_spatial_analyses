####################################    Fire Alaska Analyses   #######################################
############################  Analyse using time series and fire locations  #######################################
#This script prepares data for classes and workshops.
# The overall goal is to explore the impact of fire on surface state variables observed from Remote Sensing.     
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/09/2017 
#DATE MODIFIED: 03/15/2017
#Version: 2
#PROJECT: AAG 2017 workshop preparation
#TO DO:
#
#COMMIT: masking data and transforming to tiff
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

in_dir_NDVI <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Fire_Alaska_analyses/NDVI_Alaska"

out_dir <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Fire_Alaska_analyses/outputs"

mask_alaska_filename <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Fire_Alaska_analyses/Lab_time_series/Data_unprocessed/mask_Alaska_11112014.rst" 

proj_modis_str <-"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" #CONST 1
#CRS_interp <-"+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0" #Station coords WGS84 # CONST 2
proj_str<- CRS_WGS84 
CRS_reg <- CRS_WGS84 # PARAM 4

masking_process <- F
file_format <- ".tif" #PARAM5
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
freq(r_mask)# Alaska mask count of pixels

lf <- mixedsort(list.files(path=in_dir_NDVI,pattern="NDVI_2001_2009_filling6__NDVI_gap_filling2_07062011_test1_.*.rst$",full.names=T))

r_NDVI_data <- stack(lf)

plot(r_NDVI_data,y=1)

if(masking_process==T){
  r_NDVI_data_mask <- mask(r_NDVI_data,mask=r_mask,filename=file.path(out_dir,"r_NDVI_data_alaska_mask3.tif"))
  
  date_range <- c("2001.01.01","2009.12.31") #NDVI Katrina
  
  #generate dates for 16 days product
  dates_val <- generate_dates_by_step(date_range[1],date_range[2],16)$dates #NDVI Katrina
  dates_val_str <- format(dates_val,format="%Y_%m_%d") #format to insert in names of raste files, convert "-" in "_"
  
  out_suffix_str <- paste0("NDVI_",dates_val_str) # this needs to be matching the number of outputs files writeRaste
  
  #note idrisi files are 19 mb compared to 6.3 after compression in the tif format
  writeRaster(r_NDVI_data_mask,filename="alaska.tif",
              bylayer=T,datatype="FLT4S",options="COMPRESS=LZW",suffix=out_suffix_str,overwrite=T)
  
}else{
  r_NDVI_data_mask <- brick(file.path(out_dir,"r_NDVI_data_alaska_mask3.tif"))
}

##
plot(r_NDVI_data_mask,y=2)
filename(r_NDVI_data_mask) 

lf_masked <- mixedsort(list.files(path= out_dir, pattern="^alaska_NDVI_.*.tif",full.names=T))

#Comparing size of files after conversion, masking and compression in tif
file.info(lf_masked[1])$size/(1024*1024) # this is in bytes, convert to mb
file.info(lf[1])$size/(1024*1024)

file.info(lf_masked[1])$size/(1024*1024)/(file.info(lf[1])$size/(1024*1024)) # this is about 35.8% the size of the original raster file 

dim(r_NDVI_data_mask)

##### Now convert other files from *.rst to *.tif

#~/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Fire_Alaska_analyses/Lab_time_series/Data_unprocessed
lf_var <- list.files(path="~/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Fire_Alaska_analyses/Lab_time_series/Data_unprocessed",pattern="*.rst$",full.names=T)

r_var <- stack(lf_var)
dim(r_var)
#plot(r_var,8)
#freq(subset(r_var,8))

## now mask and transform one by one...
#r_var_data_mask <- mask(r_var,mask=r_mask,filename=file.path(out_dir,"r_var_masked.tif"),overwrite=T)
#Error in .checkLevels(levs[[j]], value[[j]]) : 
#  new raster attributes (factor values) should be in a data.frame (inside a list)
# need to fix issue with layer 8: wwf_terr_ecos_Alaska_ECOREGIONS_ECOSYS_ALB83.rst

r_var_data_mask <- mask(subset(r_var,1:7),mask=r_mask,filename=file.path(out_dir,"r_var_masked.tif"),overwrite=T)
plot(r_var_data_mask,y=1:7)

out_suffix_str <- names(r_var)[1:7] # this needs to be matching the number of outputs files writeRaste
#Fix here problem.
#note idrisi files are 19 mb compared to 6.3 after compression in the tif format
writeRaster(r_var_data_mask,filename="r.tif",
            bylayer=T,datatype="FLT4S",options="COMPRESS=LZW",suffix=out_suffix_str,overwrite=T)

### Now deal with ecoregions
r_wwf <- subset(r_var,8)
r_wwf_mask <- mask(r_wwf,r_mask)
plot(r_wwf_mask,colNA=c("black"))

raster_name <- paste0(names(r_wwf),file_format,sep="")
r_wwf <- writeRaster(r_wwf_mask,filename=file.path(out_dir,raster_name),
            bylayer=F,datatype="FLT4S",options="COMPRESS=LZW",overwrite=T)
filename(r_wwf)
plot(r_wwf)

################### End of Script #########################

#To explore later:
#http://stackoverflow.com/questions/12670972/doing-pca-on-very-large-data-set-in-r
