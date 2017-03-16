####################################    Spatial Analyses: MA   #######################################
############################  Analyse data from Census #######################################
#This script performs basic analyses for the Exercise 1 of the workshop using Census data.
# The overall goal is to explore spatial autocorrelation and aggregation of units of analyses.     
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/16/2017 
#DATE MODIFIED: 03/16/2017
#Version: 1
#PROJECT: AAG 2017 workshop preparation
#TO DO:
#
#COMMIT: initial commit for exercise 1, AAG workshop
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

in_dir_var <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/WS_AAG2017/MAUP_Lab"
out_dir <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/WS_AAG2017/MAUP_Lab/outputs"

#CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"exercise1_03162017" #output suffix for the files and ouptu folder #PARAM 8
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

census_blocks_2010_fname <- "CENSUS2010BLOCKS_POLY.shp"
census_blockgroups_2010_fname <- "CENSUS2010BLOCKGROUPS_POLY.shp"

#stat_reg <- readOGR(dsn=dirname(census_blocks_2010_fname),sub(".shp","",basename(census_blocks_2010_fname))
#blocks2010_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(census_blocks_2010_fname))) #too large for workshop
blockgroups2010_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(census_blockgroups_2010_fname)))
spplot(blockgroups2010_sp,"POP100_RE")   # somewhat slow

## Aggregate by census tracks
#ag5 <- aggregate(x = wards2, by = census_tracts_in_wards, FUN = mean)


###################### END OF SCRIPT #####################
