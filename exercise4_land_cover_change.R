####################################   Land Use and Land Cover Change   #######################################
############################  Analyze Land Cover change in Houston  #######################################
#This script performs analyses for the Exercise 4 of the Short Course using reflectance data derived from MODIS.
#The goal is to assess land cover change using two land cover maps.
#Additional datasets are provided for the land cover change modeling. 
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/16/2018 
#DATE MODIFIED: 03/16/2018
#Version: 1
#PROJECT: SESYNC and AAG 2018 workshop/Short Course preparation
#TO DO:
#
#COMMIT: initial commit exercise 4
#
#################################################################################################

###Loading R library and packages                                                      

library(sp) # spatial/geographfic objects and functions
library(rgdal) #GDAL/OGR binding for R with functionalities
library(spdep) #spatial analyses operations, functions etc.
library(gtools) # contains mixsort and other useful functions
library(maptools) # tools to manipulate spatial data
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
library(BMS) #contains hex2bin and bin2hex, Bayesian methods
library(bitops) # function for bitwise operations
library(foreign) # import datasets from SAS, spss, stata and other sources
library(gdata) #read xls, dbf etc., not recently updated but useful
library(classInt) #methods to generate class limits
library(plyr) #data wrangling: various operations for splitting, combining data
library(gstat) #spatial interpolation and kriging methods
library(readxl) #functionalities to read in excel type data
library(psych) #pca/eigenvector decomposition functionalities
library(sf)
library(plotrix) #various graphic functions e.g. draw.circle

###### Functions used in this script

create_dir_fun <- function(outDir,out_suffix=NULL){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

#function_analyses <- "exercise2_fire_alaska_analyses_functions_03232017.R" #PARAM 1
script_path <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/R_scripts"
#source(file.path(script_path,function_preprocessing_and_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#in_dir_reflectance <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/Exercise_6/data/reflectance_RITA"
in_dir_var <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_4/data/"
out_dir <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_4/outputs"
infile_reg_outline <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_4/revised_area_Rita/new_strata_rita_10282017.shp"

#region coordinate reference system
#http://spatialreference.org/ref/epsg/nad83-texas-state-mapping-system/proj4/
CRS_reg <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 
file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"exercise4_03152018" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
date_event <- ""
#ARG4
method_proj_val <- "bilinear" # "ngb"

#ARG9
#local raster name defining resolution, extent
ref_rast_name <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_4/data/r_ref_Houston_RITA.tif"

################# START SCRIPT ###############################

### PART I: READ AND PREPARE DATA FOR ANALYSES #######


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

# agg_3_r_nlcd2006_Houston.tif
# agg_3_r_nlcd2011_Houston.tif
# nlcd_legend.txt
# r_2006_nlcd30m_Houston.tif.aux.xml
# r_2006_nlcd30m_Houston.tif
# r_2011_nlcd30m_Houston.tif.aux.xml
# r_2011_nlcd30m_Houston.tif

infile_land_cover_date1 <- "r_2006_nlcd30m_Houston.tif"
infile_land_cover_date2 <- "r_2011_nlcd30m_Houston.tif"
 
r_lc_date1 <- raster(infile_land_cover_date1) 
r_lc_date2 <- raster(infile_land_cover_date2) 

lc_legend_df <- read.table("nlcd_legend.txt",sep=",")
lc_legend_df
View(lc_legend_df)
plot(r_lc_date1==90)
plot(r_lc_date2==95)
plot(r_lc_date1==11)

plot(r_lc_date1)

freq_tb_date1 <- freq(r_lc_date1)
freq_tb_date2 <- freq(r_lc_date2)

View(freq_tb_date1)
View(freq_tb_date2)

freq_tb_date1 <- merge(freq_tb_date1,lc_legend_df,by.x="value",by.y="ID")
freq_tb_date2 <- merge(freq_tb_date2,lc_legend_df,by.x="value",by.y="ID")

freq_tb_date1$date <- 2006
freq_tb_date2$date <- 2011

### No category disappeared:
freq_tb_date1$value==freq_tb_date2$value

lc_df <- data.frame(ID=freq_tb_date1$value,
           lc2006=freq_tb_date1$COUNT,
           lc2011=freq_tb_date2$COUNT,
           name=freq_tb_date1$NLCD.2006.Land.Cover.Class)

lc_df$diff <- lc_df$lc2011 - lc_df$lc2006 

View(lc_df)

################### End of Script #########################

