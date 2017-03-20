####################################    Spatial Analyses: Suitability analyses   #######################################
############################  Analyse data from Florida #######################################
#This script performs basic analyses for the Exercise 3 of the workshop using Florida data.
# The overall goal is to perform a multi-criteria/sustainabiltiy analysis to select areas suitable for conservation.     
#
#Goal: Determine the ten (10) parcels of land within Clay County, FL most suitable for purchase
#towards conversion to conservation lands.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/17/2017 
#DATE MODIFIED: 03/17/2017
#Version: 1
#PROJECT: AAG 2017 workshop preparation
#TO DO:
#
#COMMIT: processing data for the exercise, AAG workshop
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

in_dir_var <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Exercise_3/data"
out_dir <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Exercise_3/outputs"

#CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"exercise3_03172017" #output suffix for the files and ouptu folder #PARAM 8
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

##Inputs
#1) Strategic Habitat conservation areas raster file
#2) Tiger Roads shapefile
#3) County shapefile
#4) Priority Wetlands Habitat raster file
#5) Clay County parcel shapefile
#6) General Habitat raster file
#7) Biodiversity hotspot raster file
#8) Florida managed areas shapefile

#tl_2011_12019_roads.shp
#tl_2011_12019_roads_prj.shp
#tl_2010_12019_roads.shp

strat_hab_fname <- "Strat_hab_con_areas1.tif" #1)Strategic Habitat conservation areas raster file
roads_fname <- "tl_2011_12019_roads_prj.shp" #2) Missing:Tiger Roads shapefile
regional_counties_fname <- "Regional_Counties.shp" #3) County shapefile
priority_wet_habitats_fname <- "Priority_Wet_Habitats1.tif" #4) Priority Wetlands Habitat raster file
clay_parcels_fname <- "Clay_Parcels.shp" #5) Clay County parcel shapefile
habitat_fname <- "Habitat.tif" #6) General Habitat raster file
biodiversity_hotspot_fname <- "Biodiversity_Hot_Spots1.tif" #7) Biodiversity hotspot raster file
florida_managed_areas_fname <- "flma_jun13.shp" #8) Florida managed areas shapefile

## read in
r_strat_hab <- raster(file.path(in_dir_var,strat_hab_fname))
roads_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(roads_fname))) #too large for workshop
reg_counties_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(regional_counties_fname))) #too large for workshop
r_priority_wet_hab <- raster(file.path(in_dir_var,priority_wet_habitats_fname))
clay_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(clay_parcels_fname))) #too large for workshop
r_habitat <- raster(file.path(in_dir_var,habitat_fname))
r_bio_hotspot <- raster(file.path(in_dir_var,biodiversity_hotspot_fname))
fma_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(florida_managed_areas_fname))) #too large for workshop

plot(r_strat_hab)
plot(r_priority_wet_hab)
#plot(r_habitat,add=T)
plot(roads_sp,add=T)


###################### END OF SCRIPT #####################
