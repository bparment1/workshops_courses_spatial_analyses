####################################    Fire Alaska Analyses   #######################################
############################  Analyses using time series and fire locations  #######################################
#This script contains functions used in the analyses of Exercise 2 for the workshop using NDVI data.
# The overall goal is to explore the impact of fire on surface state variables observed from Remote Sensing.     
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/23/2017 
#DATE MODIFIED: 03/23/2017
#Version: 1
#PROJECT: AAG 2017 workshop preparation
#TO DO:
#
#COMMIT: initial commit AAG workshop
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

#function_preprocessing_and_analyses <- "fire_alaska_analyses_preprocessing_functions_03102017.R" #PARAM 1
#script_path <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/R_scripts"
#source(file.path(script_path,function_preprocessing_and_analyses)) #source all functions used in this script 1.

pca_to_raster_fun<-function(pc_spdf,ref_raster,NA_val,file_format=".tif",out_suffix=""){
  #Input arguments:
  #pc_spdf: spdf with scores components and
  #         must include x,y in the last 2 columns
  #ref_raster: reference raster giving the output extent and NA if present
  #NA_val<- values assigned to NA pixels, if "NA", use default
  
  npc<-ncol(pc_spdf)-2 #number of components 
  pc_scores_lf<-vector("list",npc) 
  for (k in 1:npc){
    pc_scores<-pc_spdf[,k] #vector containing scores for components k
    pc_name<-names(pc_spdf)[k] #name of the current component
    raster_name<-paste("pc_component_",k,"_",out_suffix,file_format,sep="")
    
    if (NA_val=="NA"){
      pc_lag<-rasterize(pc_scores,ref_raster,pc_name,fun=min,overwrite=TRUE,
                        filename=raster_name)
    }
    #Use provided NA value for the background
    if (NA_val!="NA"){
      pc_lag <- rasterize(pc_scores,ref_raster,pc_name,background=NA_val,fun=min,overwrite=TRUE,
                          filename=raster_name)
      NAvalue(pc_lag) <- NA_val
    }
    pc_scores_lf[[k]]<-raster_name
  } 
  return(pc_scores_lf)
}

################### End of Script #########################

#To explore later:
#http://stackoverflow.com/questions/12670972/doing-pca-on-very-large-data-set-in-r
