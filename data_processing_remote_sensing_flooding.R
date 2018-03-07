####################################    Flood Mapping Analyses   #######################################
############################  Analyze and map flooding from RITA hurricane  #######################################
#This script performs data processing for the Exercise 4, Exercise 5 and Exercise 6 for the Geospatial Short Course 
#This is using reflectance data derived from MODIS MOD09 and NLCD data.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/07/2018 
#DATE MODIFIED: 03/07/2018
#Version: 1
#PROJECT: SESYNC and AAG 2018 workshop/Short Course preparation
#TO DO:
#
#COMMIT: initial commit
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

###### Functions used in this script

data_processing_functions <- "data_processing_remote_sensing_flooding_functions_03072018c.R" #PARAM 1
script_path <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/R_scripts"
source(file.path(script_path,data_processing_functions)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

in_dir_reflectance <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/Exercise_5/data/reflectance_RITA"
in_dir_var <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/Exercise_5/data/"
out_dir <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/Exercise_5/outputs"

#infile_ecoreg <- "wwf_terr_ecos_Alaska_ECOREGIONS_ECOSYS_ALB83.shp" #WWF ecoregions 2001 for Alaska
infile_reg_outline <- "/nfs/bparmentier-data/Data/Space_beats_time/Data/data_RITA_reflectance/revised_area_Rita/new_strata_rita_10282017.shp"

#region coordinate reference system
#http://spatialreference.org/ref/epsg/nad83-texas-state-mapping-system/proj4/
CRS_reg <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 

file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"exercise5_03052018" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
date_event <- ""

#ARG4
method_proj_val <- "bilinear" # "ngb"

#ARG9
#local raster name defining resolution, extent
ref_rast_name <- "/nfs/bparmentier-data/Data/Space_beats_time/Data/data_RITA_reflectance/revised_area_Rita/r_ref_Houston_RITA.tif"
#ref_rast_name <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Harvey_NDVI/revised_area_Rita/r_ref_Houston_RITA.tif"
#ref_rast_name <- "/home/bparmentier/Google Drive/Space_beats_time/Data/data_Harvey_NDVI/rita_outline_reg/Study_Area_Rita_New.shp"
#ARG11
date_param <- "2005.01.01;2005.12.31;8" #start date, end date, time_step
#ARG12
#ARG13
#scaling_factors <- c(1,-273.15) #set up as slope (a) and intercept (b), if NULL, no scaling done, setting for LST 
scaling_factors <- c(0.0001,0) #set up as slope (a) and intercept (b), if NULL, no scaling done, setting for NDVI 
#ARGS14
product_type = c("reflectance") #can be LST, ALBEDO etc.#this can be set from the modis product!! #param 19
#ARG15
multiband <- TRUE #This is only used for multiband products?
#ARGS16: This can be removed in the future by stating LST_Day as a product type
num_cores <- 4 #option for parallel processes

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

lf_reflectance <- list.files(path=in_dir_reflectance, pattern="*.tif",full.names=T)
r_refl_ts <- stack(lf_reflectance) #          #note this creates 46*7 bands

plot(r_NDVI_ts)
date_range <- unlist(strsplit(date_param,";")) #NDVI Alaska, year 2005 (this is a 16 days product)
  
df_dates <- as.data.frame(generate_dates_by_step(start_date=date_range[1],
                                                 end_date=date_range[2],
                                                 step_date=as.numeric(date_range[3])))
#34: 2005265 this is Sept 22
#35: 2005273 this is Sept 30

#generate dates for 16 days product
#dates_val <- generate_dates_by_step(date_range[1],date_range[2],16)$dates #NDVI Katrina

file.info(lf_reflectance[1])$size/(1024*1024) # this is in bytes, convert to mb
dataType(r_refl_ts) #Examine the data type used in the storing of data, this is float 32 signed: FLT4S
inMemory(r_refl_ts) #Is the data in memory? Raster package does not load in memory automatically.
dim(r_refl_ts) #dimension of the raster object: rows, cols, layers/bands


######################## PART I: PROCESS NLCD ##############

r_2006_nlcd30m <- raster("/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/nlcd_2006_landcover_2011_edition_2014_10_10/nlcd_2006_landcover_2011_edition_2014_10_10.img")
r_2011_nlcd30m <- raster("/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img")

#r_nlcd_2006_30m <- raster("/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/nlcd_2006_landcover_2011_edition_2014_10_10/nlcd_2006.tif")
r_2006_nlcd30m
r_2011_nlcd30m
dataType(r_2006_nlcd30m)
dataType(r_2011_nlcd30m)

reg_sf_nlcd <- st_transform(reg_sf,projection(r_2006_nlcd30m))
reg_sp_nlcd <- as(reg_sf_nlcd,"Spatial")
r_2006_nlcd30m_RITA <- crop(r_2006_nlcd30m,reg_sp_nlcd,"r_2006_nlcd30m.tif",overwrite=T)
r_2011_nlcd30m_RITA <- crop(r_2011_nlcd30m,reg_sp_nlcd,"r_2011_nlcd30m.tif",overwrite=T)

plot(r_2006_nlcd30m_RITA)
plot(r_2011_nlcd30m_RITA)

#freq(r_nlcd30m_RITA)

## input files to aggregate
l_rast <- list(r_2006_nlcd30m_RITA,r_2011_nlcd30m_RITA)
#cat_names <- NULL
names(l_rast) <- c("nlcd2006_RITA","nlcd2011_RITA")
cat_names <- c("nlcd2006_RITA","nlcd2011_RITA")

rast_ref_1km <- "/nfs/bparmentier-data/Data/Space_beats_time/Data/data_RITA_reflectance/revised_area_Rita/r_ref_Houston_RITA.tif"

### Get to 1km:
debug(aggregate_raster_fun)


agg_obj_1km <- aggregate_raster_fun(l_rast,
                                    cat_names=cat_names,
                                    use_majority=T,
                                    agg_fact=NULL, #if null will look for the ref image to determine
                                    agg_fun=mean,
                                    file_format=file_format,
                                    rast_ref=rast_ref_1km,
                                    num_cores=num_cores,
                                    out_suffix=out_suffix, 
                                    out_dir=out_dir)

################## PART II: MODIS REFLECTANCE ##############

  
#lf_var <- list.files(path=in_dir_var,pattern="*.tif$",full.names=T)
#r_var <- stack(lf_var) # create a raster stack, this is not directly stored in memory

#dim(r_var) #dimension of the stack with 
#plot(r_var)

reg_sf <- st_read(infile_reg_outline)
reg_sf <- st_transform(reg_sf,
                       crs=CRS_reg)
reg_sp <-as(reg_sf, "Spatial") 
plot(reg_sf$geometry)

#r_tmp <- subset(r_refl_ts,1)
r_tmp <- brick(lf_reflectance[1])

r_ref <- rasterize(reg_sp,
                   r_tmp,
                   field="OBJECTID_1",
                   fun="first")

r_before <- brick(lf_reflectance[34])
r_after <- brick(lf_reflectance[35])

r_before <- r_before*(1/0.0001)
r_after <- r_after*(1/0.0001)

#SWIR1 (1230–1250 nm), SWIR2 (1628–1652 nm) and SWIR3 (2105–2155 nm).
band_refl_order <- c(3,4,1,2,5,6,7)

names(r_before) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")
names(r_after) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")

#resample()


################### End of Script #########################
