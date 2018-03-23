####################################    Flood Mapping Analyses   #######################################
############################  Analyze and map flooding from RITA hurricane  #######################################
#This script performs data processing for the Exercise 4, Exercise 5 and Exercise 6 for the Geospatial Short Course 
#This is using reflectance data derived from MODIS MOD09 and NLCD data.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/07/2018 
#DATE MODIFIED: 03/23/2018
#Version: 1
#PROJECT: SESYNC and AAG 2018 workshop/Short Course preparation
#TO DO:
#
#COMMIT: generate input mod09 images
#
#################################################################################################

###Loading R library and packages                                                      
#library(gstat) #spatial interpolation and kriging methods
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
library(readxl) #functionalities to read in excel type data
library(psych) #pca/eigenvector decomposition functionalities
library(sf)
library(tigris)

###### Functions used in this script

data_processing_functions <- "data_processing_remote_sensing_flooding_functions_03162018.R" #PARAM 1
script_path <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/R_scripts"

source(file.path(script_path,data_processing_functions)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

in_dir_reflectance <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/reflectance_RITA"
in_dir_var <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/"
out_dir <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/"

#Source:https://cohgis-mycity.opendata.arcgis.com/datasets/houston-city-limit
infile_reg_outline_Houston_city_limits <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/Houston_City_Limit/Houston_City_Limit.shp"
#                                           /nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/Houston_City_Limit
infile_reg_outline_RITA <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/revised_area_Rita/new_strata_rita_10282017.shp"
infilename_2006_nlcd30m <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/nlcd_2006_landcover_2011_edition_2014_10_10/nlcd_2006_landcover_2011_edition_2014_10_10.img"
infilename_2001_nlcd30m <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/nlcd_2001_landcover_2011_edition_2014_10_10/nlcd_2001_landcover_2011_edition_2014_10_10.img"
infilename_2011_nlcd30m <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img"

#region coordinate reference system
#http://spatialreference.org/ref/epsg/nad83-texas-state-mapping-system/proj4/
CRS_reg <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 

file_format <- ".tif" #PARAM5
NA_flag_val <- -9999 #PARAM7
out_suffix <-"data_preprocessing_03212018" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
date_event <- ""

#ARG4
method_proj_val <- "bilinear" # "ngb"

#ARG9
#local raster name defining resolution, extent
ref_rast_name <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/revised_area_Rita/r_ref_Houston_RITA.tif"
#ARG11
date_param <- "2005.01.01;2005.12.31;8" #start date, end date, time_step
#ARG13
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

################## PART II: PROCESSING MODIS REFLECTANCE ##############

#lf_var <- list.files(path=in_dir_var,pattern="*.tif$",full.names=T)
#r_var <- stack(lf_var) # create a raster stack, this is not directly stored in memory

#dim(r_var) #dimension of the stack with 
#plot(r_var)

reg_sf <- st_read(infile_reg_outline_RITA)
reg_sf <- st_transform(reg_sf,
                       crs=CRS_reg)
reg_sp <- as(reg_sf, "Spatial") 
plot(reg_sf$geometry)

#this is at 926m
rast_ref <- raster(ref_rast_name)
projection(rast_ref) <- CRS_reg
#r_tmp <- subset(r_refl_ts,1)
r_tmp <- brick(lf_reflectance[1])

r_ref <- rasterize(reg_sp,
                   r_tmp,
                   field="OBJECTID_1",
                   fun="first")
res(r_ref) 
#[1] 463.3127 463.3127

r_before <- brick(lf_reflectance[34])
r_after <- brick(lf_reflectance[35])

r_before <- r_before*(1/0.0001)
r_after <- r_after*(1/0.0001)

projection(r_before) <- CRS_reg
projection(r_after) <- CRS_reg

#SWIR1 (1230–1250 nm), SWIR2 (1628–1652 nm) and SWIR3 (2105–2155 nm).
band_refl_order <- c(3,4,1,2,5,6,7)

names(r_before) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")
names(r_after) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")

#resample()
plot(r_before)

######## Let's reproject and match to the resolution ofr 926m, the NDVI data
r_before <- projectRaster(r_before,rast_ref,method="bilinear")
r_after <- projectRaster(r_after,rast_ref,method="bilinear")
plot(r_after)
plot(r_before)
res(r_before)
res(r_after)

#### now write out data at 926m

raster_name_tmp <- lf_reflectance[34]
raster_name_tmp <- "mosaiced_MOD09A1_A2005265__006_reflectance_masked_RITA_reg_1km.tif"
bylayer_val <- FALSE
out_suffix_str <- NULL
data_type_str <- dataType(r_before)

writeRaster(r_before,
            filename=file.path(out_dir,raster_name_tmp),
            bylayer=bylayer_val,
            suffix=suffix_str,
            overwrite=TRUE,
            NAflag=NA_flag_val,
            datatype=data_type_str,
            options=c("COMPRESS=LZW"))

### Now second date:

raster_name_tmp <- lf_reflectance[35]
raster_name_tmp <- "mosaiced_MOD09A1_A2005273__006_reflectance_masked_RITA_reg_1km.tif"
bylayer_val <- FALSE
out_suffix_str <- NULL
data_type_str <- dataType(r_after)

writeRaster(r_after,
            filename=file.path(out_dir,raster_name_tmp),
            bylayer=bylayer_val,
            suffix=suffix_str,
            overwrite=TRUE,
            NAflag=NA_flag_val,
            datatype=data_type_str,
            options=c("COMPRESS=LZW"))

#####################
## Creating a true color composite with streching

histogram(r_before)
histogram(r_after)

minValue(r_before)
#[1] 0.0003000000 0.0000000000 0.0001000000 0.0073000001 0.0000000000 0.0013417750 0.0002549418
maxValue(r_before)
#[1] 0.5316203 0.5905742 0.4861007 0.5169186 0.5570464 0.5135773 0.4598759

min_val <- minValue(r_before)
max_val <- maxValue(r_before)

r <- subset(r_before,1)
r_test <- round(255*(r-min_val[1])/(max_val[1]-min_val[1]))
q_val <- quantile(r,probs=seq(0,1,0.01))


i <- 1

scale_rast_fun <- function(i,r_stack,min_val=NULL,max_val=NULL){
  #
  #
  
  r <- subset(r_stack,i)
  if(is.null(min_val)){
    min_val <- minValue(r)
  }
  if(is.null(max_val)){
    max_val <- maxValue(r)
  }
  r_scaled <- round(255*(r-min_val)/(max_val-min_val))
  return(r_scaled)
}

r_red <- scale_rast_fun(1,r_before)
r_blue <- scale_rast_fun(3,r_before)
r_green <- scale_rast_fun(4,r_before)
r_nir <- scale_rast_fun(2,r_before)

r_test <- stretch(r_red,minq=0,maxq=100)
histogram(r_test)
histogram(r_red)
plot(r_test)

histogram(r_red)
plot(r_red)
plot(r_blue)
plot(r_green)

r_rgb <- stack(r_red,r_green,r_blue,r_nir)

plotRGB(r_rgb,
        r=1,
        g=2,
        b=3,
        scale=255,
        strech="hist")

### False color composite:

plotRGB(r_rgb,
        r=4,
        g=1,
        b=2,
        #scale=255,
        strech="hist")

#R = XS3 (NIR band)
#G = XS2 (red band)
#B = XS1 (green band)

writeRaster(r_red,"r_red.rst")
writeRaster(r_blue,"r_blue.rst")
writeRaster(r_green,"r_green.rst")

r_red <- scale_rast_fun(1,r_before,0,0.3)
r_blue <- scale_rast_fun(3,r_before,0,0.3)
r_green <- scale_rast_fun(4,r_before,0,0.3)
r_nir <- scale_rast_fun(2,r_before,0,0.3)

r_rgb <- stack(r_red,r_green,r_blue,r_nir)

plotRGB(r_rgb,
        r=1,
        g=2,
        b=3,
        #scale=255,
        strech="hist")
### False color composite:

plotRGB(r_rgb,
        r=4,
        g=1,
        b=2,
        #scale=255,
        strech="hist")

#R = XS3 (NIR band)
#G = XS2 (red band)
#B = XS1 (green band)

r <- subset(r_before,3)

r_blue <- round(255*(r-min_val[1])/(max_val[1]-min_val[1]))

r <- subset(r_before,4)
r_blue <- round(255*(r-min_val[1])/(max_val[1]-min_val[1]))

r_test <- (r-0)/(0.2-0)

r_test <- (r-0)*(255-0)
plot(r_test)

plot(r_test)
?quantile

quantile(r, probs = c(0.25, 0.75), type=7,names = FALSE)

plot(r_test)
histogram(r_test)
histogram(r)


###### For Houston area:
#mosaiced_MOD09A1_A2005265__006_reflectance_masked_RITA_02132018.tif
#mosaiced_MOD09A1_A2005273__006_reflectance_masked_RITA_02132018.tif

mosaic_filename <- "/nfs/bparmentier-data/Data/Space_beats_time/Data/data_RITA_reflectance/mosaic_output/mosaiced_MOD09A1_A2005265__006_reflectance_masked_RITA_02132018.tif"
r_mosaiced <- brick(mosaic_filename)

plot(r_mosaiced)
#crop()

#ref_rast_tmp <-raster(list_var_mosaiced[[1]]) 
method_proj_val <- "bilinear" 
#method_proj_val <- "ngb" 

ref_rast_prj <-projectRaster(from=r_mosaiced,
                             res=res(r_mosaiced), #set resolution to the same as input
                             crs=CRS_reg,
                             method=method_proj_val)

##This is the mosaiced and reproject tiles matching:
ref_rast_prj_name_generated <- paste("ref_mosaiced_input_rast_",out_suffix,file_format,sep="")
writeRaster( ref_rast_prj,file.path(out_dir,ref_rast_prj_name_generated))

#to define a local reference system and reproject later!!
#Assign new projection system here in the argument CRS_reg (!it is used later)

infile_reg_outline <- infile_reg_outline_Houston_city_limits
reg_sf <- st_read(infile_reg_outline)
reg_sf <- st_transform(reg_sf,crs=CRS_reg)
reg_sp <-as(reg_sf, "Spatial") 

ref_rast <- crop(ref_rast_prj,reg_sp)  
ref_rast_name_generated <- paste("ref_rast_generated_Houston_",out_suffix,file_format,sep="")
writeRaster(ref_rast,file.path(out_dir,ref_rast_name_generated))

### This is only for 2005!!!
plot(ref_rast)
names(r_before) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")

#mosaic_filename <- "/nfs/bparmentier-data/Data/Space_beats_time/Data/data_RITA_reflectance/mosaic_output/mosaiced_MOD09A1_A2005265__006_reflectance_masked_RITA_02132018.tif"
#r_mosaiced <- brick(mosaic_filename)

######################## PART II: PROCESSING NLCD ##############

#infile_reg_outline_Houston_city_limits <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/Houston_City_Limit.shp"
#infile_reg_outline_RITA <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/revised_area_Rita/new_strata_rita_10282017.shp"

###### Need to process for two areas of study!!!

r_2001_nlcd30m <- raster(infilename_2001_nlcd30m)
r_2006_nlcd30m <- raster(infilename_2006_nlcd30m)
r_2011_nlcd30m <- raster(infilename_2011_nlcd30m)

r_2006_nlcd30m
r_2001_nlcd30m

r_2011_nlcd30m
dataType(r_2001_nlcd30m)
dataType(r_2011_nlcd30m)
legend_col <- r_2006_nlcd30m@legend
class(legend_col)

lc_nlcd_legend <- r_2006_nlcd30m@data@attributes[[1]]
class(lc_nlcd_legend) #this is a data.frame
View(lc_nlcd_legend)
names(lc_nlcd_legend)

lc_nlcd_legend_filename <- "nlcd_legend.txt"
write.table(lc_nlcd_legend,file=lc_nlcd_legend_filename,sep=",")
#lc_types <- r_2006_nlcd30m@data@attributes[[1]]$Land.Cover.Class
#unique(legend_col@colortable)
#[1] "#000000" "#00F900" "#476BA0" "#D1DDF9" "#DDC9C9" "#D89382"
#[7] "#ED0000" "#AA0000" "#B2ADA3" "#F9F9F9" "#68AA63" "#1C6330"
#[13] "#B5C98E" "#A58C30" "#CCBA7C" "#E2E2C1" "#C9C977" "#99C147"
#[19] "#77AD93" "#DBD83D" "#AA7028" "#BAD8EA" "#B5D3E5" "#70A3BA"

#Need to reclass values in NLCD and plot different classes

########### Let's crop to match the area of interests: Houston and RITA  ####

#infile_reg_outline_Houston_city_limits <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/Houston_City_Limit.shp"
#infile_reg_outline_RITA <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/revised_area_Rita/new_strata_rita_10282017.shp"

### Now generate NLCD for each region:

list_reg_outline <- list(infile_reg_outline_RITA,infile_reg_outline_Houston_city_limits)
#list_ref_filename <- 
#Make a function later:
#inputs:
#ref
#this is at 926m
rast_ref_filename <- "/nfs/bparmentier-data/Data/Space_beats_time/Data/data_RITA_reflectance/revised_area_Rita/r_ref_Houston_RITA.tif"
list_out_suffix <- list("RITA","Houston")
list_ref_rast_name <- list(rast_ref_filename,NULL)
list_agg_fact_val <- list("NULL",3) 
#names(r_before) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")

i<-2

for(i in 1:length(list_reg_outline)){
  #
  #
  
  infile_reg_outline <- list_reg_outline[[i]]
  out_suffix_str <- list_out_suffix[[i]]
  ref_rast_name <- list_ref_rast_name[[i]]
  agg_fact_val <- list_agg_fact_val[[i]] #3 for Houston
  
  reg_sf <- st_read(infile_reg_outline)
  reg_sf <- st_transform(reg_sf,
                              crs=CRS_reg)
  reg_sp <- as(reg_sf, "Spatial") 
  plot(reg_sf$geometry)
  
  #this is at 926m
  if(!is.null(ref_rast_name)){
    rast_ref <- raster(ref_rast_name)
    projection(rast_ref) <- CRS_reg
    #r_tmp <- subset(r_refl_ts,1)
    
  }else{
    rast_ref <- NULL
  }
  
  reg_sf_nlcd <- st_transform(reg_sf,projection(r_2001_nlcd30m))
  reg_sp_nlcd <- as(reg_sf_nlcd,"Spatial")
  #out_suffix_str
  r_2001_nlcd30m_reg <- crop(r_2001_nlcd30m,
                         reg_sp_nlcd,
                         paste0("r_2001_nlcd30m_",out_suffix_str,file_format),
                         overwrite=T)
  r_2011_nlcd30m_reg <- crop(r_2011_nlcd30m,
                         reg_sp_nlcd,
                         filename=paste0("r_2011_nlcd30m_",out_suffix_str,file_format),
                         overwrite=T)
  
  plot(r_2001_nlcd30m_reg,main="2006")
  plot(r_2011_nlcd30m_reg,main="2011")
  
  #freq(r_nlcd30m_RITA)
  
  ## input files to aggregate
  l_rast <- list(r_2001_nlcd30m,r_2011_nlcd30m)
  #cat_names <- NULL
  names(l_rast) <- c("nlcd2001","nlcd2011")
  cat_names <- c("nlcd2001","nlcd2011")
  
  ### Get to 1km (or ~926m):
  #debug(aggregate_raster_fun)
  #agg_fact_val <- 3
  #rast_ref <- NULL
  
  agg_obj <- aggregate_raster_fun(l_rast,
                                      cat_names=cat_names,
                                      agg_method_cat="majority",
                                      agg_fact=agg_fact_val, #if null will look for the ref image to determine
                                      agg_fun=mean,
                                      file_format=file_format,
                                      rast_ref=rast_ref,
                                      num_cores=num_cores,
                                      out_suffix=out_suffix_str, 
                                      out_dir=out_dir)
  agg_obj$l_rast_cat
  names(agg_obj)
  #rast_agg_nlcd2006_aea <- raster("agg_31_r_nlcd2006_RITA_nlcd2006_RITA_data_preprocessing_03142018.tif")
  #rast_agg_nlcd2011_aea_RITA <- raster("agg_31_r_nlcd2011_RITA_nlcd2011_RITA_data_preprocessing_03142018.tif")
  
  if(!is.null(rast_ref)){
    nlcd2006_agg_reg <- projectRaster(agg_obj$l_rast_cat[[1]],
                                      rast_ref,
                                      method ="ngb")
    
    plot(nlcd2006_agg_reg)
    nlcd2011_reg_RITA <- projectRaster(agg_obj$l_rast_cat[[1]],
                                       rast_ref,
                                       method="ngb")
    plot(nlcd2011_agg_reg)
    
  }
  
  #out_filename <- filename(agg_obj$l_rast_cat[[1]])
  #writeRaster(nlcd2006_reg_RITA,filename = "nlcd_2006_RITA.tif",overwrite=T)
  #writeRaster(nlcd2011_reg_RITA,filename = "nlcd_2011_RITA.tif",overwrite=T)
  
}


########### Part III: Additional data ##############################

### Get srtm elevation data

#getData('SRTM', lon=5, lat=45)
#<- st_centroid(reg_sf)
#library(tigris)

#reg_sf
#bbx_reg <- (st_bbox(reg_sf))
#str(st_bbox(reg_sf))
#srtm <- getData('SRTM', lon=16, lat=48)
#srtm <- getData('SRTM', lon=bbx_reg[1], lat=bbx_reg[2],path=out_dir)

#Source:https://cohgis-mycity.opendata.arcgis.com/datasets/houston-city-limit
infile_reg_outline_Houston_city_limits <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/data/Houston_City_Limit/Houston_City_Limit.shp"

srtm_list <- list.files(path=in_dir_var,pattern="*.tif",full.names=T)
l_srtm <- lapply(srtm_list,function(x){raster(x)})

r_m <- mosaic(l_srtm[[2]],l_srtm[[3]],l_srtm[[4]],l_srtm[[5]],fun=mean)

reg_sf <- st_read(infile_reg_outline_Houston_city_limits)
reg_sf <- st_transform(reg_sf,
                       crs=projection(r_m))
reg_sp <- as(reg_sf, "Spatial") 

plot(r_m)
plot(reg_sf$geometry,add=T)
plot(reg_sp,add=T)

r_srtm <- crop(r_m,reg_sp)
plot(r_srtm)
plot(reg_sp,add=T)

srtm_out_rastername <- "srtm_Houston_area_90m.tif"
writeRaster(r_srtm,filename = file.path(out_dir,srtm_out_rastername))

### Disaggregate at 30m??

####### Get roads data
test<- roads("Texas",county="Harris County",year=2011)
roads2<- roads("Texas",county="Harris County",year=2011)

list_roads_Houston <- mclapply(list_counties,
                               FUN=function(x){roads("Texas",county=x,year=2011)},
         mc.preschedule = F,mc.cores = 4)
list_counties <- c("Harris County","Fort Bend County","Montgomery County","Brazoria County",
"Galveston County","Liberty County","Waller County","Chambers County","Austin County")

class(list_roads_Houston[[1]])
projection(list_roads_Houston[[1]])

plot(reg_sf$geometry)
plot(list_roads_Houston[[1]],add=T)
#Harris County – 4,092,459
#Fort Bend County – 585,375
#Montgomery County – 455,746
#Brazoria County – 313,166
#Galveston County – 291,309
#Liberty County – 75,643
#Waller County – 43,205
#Chambers County – 35,096
#Austin County – 28,417

r_2001_nlcd30m_reg

plot(r_2001_nlcd30m_reg)
plot(reg_sf$geometry,add=T)

### Neeed to make this a function!!

road_reg_sample <- spTransform(list_roads_Houston[[1]],projection(r_2001_nlcd30m_reg))
#?rasterize

road_reg_sample@data$STATEFP <- as.numeric(road_reg_sample@data$STATEFP)

r_roads_sample <- rasterize(x=road_reg_sample, 
          y=r_2001_nlcd30m_reg, 
          field="STATEFP", fun='count', background=NA,
          mask=FALSE, update=FALSE, updateValue='all', filename="",
          getCover=FALSE, silent=TRUE)

roads_rastername <- "r_roads_Harris.tif"
writeRaster(r_roads_sample,filename=roads_rastername)

names(road_reg_sample)
View(road_reg_sample@data)

####### Processing training and testing information

class1_data_sf <- st_read(file.path(in_dir_var,"class1.shp"))
class2_data_sf <- st_read(file.path(in_dir_var,"class2.shp"))
class3_data_sf <- st_read(file.path(in_dir_var,"class3.shp"))

class1_data_sf$class_ID <- 1
class2_data_sf$class_ID <- 2
class3_data_sf$class_ID <- 3

list_class_data_sf <- list(class1_data_sf,class2_data_sf,class3_data_sf)
names(class1_data_sf)
table(class1_data_sf$REC_NUM)
clean_up_data_sf <- function(i,x_sf){
  data_sf <- x_sf[[i]]
  data_sf$class_ID <- i
  data_sf$record_ID <- data_sf$REC_NUM
  data_sf <- data_sf[,c("record_ID","class_ID")]
  return(data_sf)
}

test_sf <- clean_up_data_sf(1,list_class_data_sf)
View(test_sf)
list_ground_sf <- lapply(1:length(list_class_data_sf),FUN=clean_up_data_sf,x_sf=list_class_data_sf)
View(list_ground_sf[[1]])

st_write(list_ground_sf[[1]],"class1_sites.shp",layer_options="OVERWRITE=yes")
st_write(list_ground_sf[[2]],"class2_sites.shp")
st_write(list_ground_sf[[3]],"class3_sites.shp")

################################## End of Script #########################################
