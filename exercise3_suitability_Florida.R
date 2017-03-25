####################################    Spatial Analyses: Suitability analyses   #######################################
############################################  Analyse data from Florida #######################################
#This script performs basic analyses for the Exercise 3 of the workshop using Florida data.
# The overall goal is to perform a multi-criteria/sustainabiltiy analysis to select areas suitable for conservation.     
#
#Goal: Determine the ten (10) parcels of land within Clay County, FL most suitable for purchase
#towards conversion to conservation lands.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/17/2017 
#DATE MODIFIED: 03/25/2017
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

gdal_installed <- TRUE #if true use the system/shell command else use the distance layer provided
file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"exercise3_03242017" #output suffix for the files and ouptu folder #PARAM 8
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
clay_parcels_fname <- "Clay_Parcels.shp" #5) Clay County parcel shapefile

plot(r_strat_hab, main="strategic habitat")
plot(reg_counties_sp,add=T)
plot(r_priority_wet_hab)
#plot(r_habitat,add=T)
plot(roads_sp,add=T)

##### Before starting the production of sustainability let's check the projection, resolution for each layer

#raster layers:

list_raster <- c(r_strat_hab,r_priority_wet_hab,r_habitat,r_bio_hotspot)

## Examine information on rasters using lapply
lapply(list_raster,function(x){res(x)}) #spatial resolution
lapply(list_raster,function(x){projection(x)}) #spatial projection
lapply(list_raster,function(x){extent(x)}) #extent of rasters

### PART 0: generate reference layer

## let's use the resolution 55x55 as the reference since it corresponds to finer resolution
## Select clay county
clay_county_sp <- subset(reg_counties_sp,NAME=="CLAY")
plot(r_strat_hab)
plot(clay_county_sp,border="red",add=T)

## Crop r_strat_hab
r_ref <- crop(r_strat_hab,clay_county_sp) #make a reference image for use in the processing
plot(r_ref)
plot(clay_county_sp,border="red",add=T)

r_clay <- rasterize(clay_county_sp,r_ref) #this can be used as mask for the study area
freq(r_clay)
##Use raster of Clay county definining the study area to mask pixels
plot(r_clay)
dim(r_clay)
#dim(r_priority_wet_hab_w)

# PART I: P1- PRIORITY1- IDENTIFY LANDS WITH HIGH NATIVE BIODIVERSITY

### STEP 1: Strategic Habitat conservation areas

#Input data layer: Habitat
#Criteria for value assignment: Habitat ranked by the Natural Heritage Program as having high native
#biodiversity were given a value of 9. Habitat ranked as having a moderate native biodiversity were given
#a value of 5, and all other habitat types were given a value of 1.
#Rationale for value assignment: Certain habitat types are known to have higher native biodiversity than
#others, consequently those with higher native biodiversity were given higher suitability rankings.
#Output: Habitat Biodiversity SUA (CG1O11SO114)

#strategic habitat conservation areas
plot(r_strat_hab)
r_strat_hab_w <- crop(r_strat_hab,r_clay) #did this up there
r_strat_hab_masked <- mask(r_strat_hab_w,r_clay)

## Now reclassify
m <- c(5, 1000, 9,  
       4, 5, 5,  
       1, 3, 1)  

rclmat <- matrix(m, ncol=3, byrow=TRUE)

#?raster::reclassify: ton find out information on 

rc_strat_hab_reg <- reclassify(r_strat_hab_masked, rclmat)
plot(rc_strat_hab_reg)

### STEP 2: Identify Lands With High Native Biodiversity: based on species count 

## Crop bio raster
r_bio_hotspot_w <- crop(r_bio_hotspot,clay_county_sp)
plot(r_bio_hotspot_w)
plot(clay_county_sp,border="red",add=T)

#r_bio_clay_masked <- mask(r_bio_hotspot_w,r_clay) ## Does not work!! because resolution don't match
#match resolution:
projection(r_bio_hotspot_w)==projection(r_clay)
res(r_bio_hotspot_w)==res(r_clay)

## Find about resample
#?raster::resample
r_bio_hotspot_reg <- resample(r_bio_hotspot_w,r_clay, method='bilinear') #Use resample to match resolutions

r_bio_hotspot_reg <- mask(r_bio_hotspot_reg,r_clay) ## It now works because resolutions were matched
plot(r_bio_hotspot_reg)

### Reclassify using instructions/informaiton given to us:
#Criteria for value assignment: Cells with 7 or more focal species or an actual species occurrence record
#location were assigned a value of 9; a value of 8 was assigned to cells with 5–8 focal species; 7 was
#assigned to cells with 3–4 focal species; all other cells were assigned a value of 1.

m <- c(9, 1000, 9,  
       5, 8, 8,  
       3, 4, 7,  
       1, 2,1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

rc_bio_hotspot_reg <- reclassify(r_bio_hotspot_reg, rclmat)

plot(rc_bio_hotspot_reg)

### STEP 3: Wetland priority 

#Input data layer: Priority Wetland Habitats
#Criteria for value assignment: Values were assigned based on the number of focal species present in
#each cell. The value of 9 was assigned to 10–12 wetland focal species, 8 was assigned to 7–9 wetland focal
#species, 7 was assigned to 4–6 wetland focal species and 4–6 upland focal species, 6 was assigned to
#1–3 wetland or upland focal species. The value 1 was assigned to all other cells.
#Rationale for value assignment: The better the habitat for focal wetland species, the higher the priority.
#Output: Wetland Biodiversity SUA (CG1O11SO111)

plot(r_priority_wet_hab)
#check projection
projection(r_priority_wet_hab)

## Crop Wetland priority raster
r_priority_wet_hab_w <- crop(r_priority_wet_hab,clay_county_sp)
plot(r_priority_wet_hab_w)
plot(clay_county_sp,border="red",add=T)

#r_priority_wet_hab_reg <- mask(r_priority_wet_hab_w,r_clay) ## Does not work!! because resolution don't match
#match resolution:
## Find about resample

r_priority_wet_hab_reg <- resample(r_priority_wet_hab_w,r_clay, method='bilinear') #resolution matching the study region

#r_priority_wet_hab_reg <- mask(r_priority_wet_hab_reg,r_clay) ## Does not work!! because resolution don't match
plot(r_priority_wet_hab_reg)

### Now reclass
#The value of 9 was assigned to 10–12 wetland focal species, 8 was assigned to 7–9 wetland focal
#species, 7 was assigned to 4–6 wetland focal species and 4–6 upland focal species, 6 was assigned to
#1–3 wetland or upland focal species. The value 1 was assigned to all other cells.
m <- c(10, 12, 9,  
       7, 9, 8,  
       4, 6, 7,  
       1, 3,6,
       -1, 1,1)

rclmat <- matrix(m, ncol=3, byrow=TRUE)

rc_priority_wet_hab_reg <- reclassify(r_priority_wet_hab_reg, rclmat)
plot(rc_priority_wet_hab_reg )
freq_tb <- freq(rc_priority_wet_hab_reg)


### STEP 4: Combine all the three layers with weigthed sum

f_weights <- c(1,1,1)/3
r_bio_es_factor <- (f_weights[1]*rc_strat_hab_reg + f_weights[2]*rc_bio_hotspot_reg + f_weights[3]*rc_priority_wet_hab_reg) #weighted sum

f_weights <- c(1,1.5,1.5)/3
r_bio_ws_factor <- (f_weights[1]*rc_strat_hab_reg + f_weights[2]*rc_bio_hotspot_reg + f_weights[3]*rc_priority_wet_hab_reg) #weighted sum

r_bio_factor <- stack(r_bio_es_factor,r_bio_ws_factor)
names(r_bio_factor) <- c("equal_weights","weigthed_sum")
plot(r_bio_factor)

##########################
#P2- IDENTIFY POTENTIAL CONSERVATION LANDS IN RELATION WITH ROAD DENSITY AND EXISTING MANAGED LANDS
#GOAL: Create two raster maps showing lands in Clay County, 
#Florida that have would have higher conservation potential based on 
#local road density and distance from existing managed lands using a combination 
#of the Clip, Extract by Mask, Euclidean Distance, Line Density, Project, Reclassify, and Weighted Sum tools.

#Step 1: prepare files to create a distance to road layer

#r_roads <- raster("roads.tif")
r_roads <- raster("roads_counts.tif")
### Processs roads first
#set1f <- function(x){rep(NA, x)}
#r_init <- init(r_clay, fun=set1f)
r_roads_bool <- r_roads > 0
NAvalue(r_roads_bool ) <- 0 
roads_bool_fname <- file.path(out_dir,paste0("roads_bool_",out_suffix,file_format))
r_roads_bool <- writeRaster(r_roads_bool,filename=roads_bool_fname,overwrite=T)

#setp 2: prepare files to create a distance to existing managed land

florida_managed_areas_fname <- "flma_jun13.shp" #8) Florida managed areas shapefile
flma_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(florida_managed_areas_fname ))) #too large for workshop
r_flma_clay <- rasterize(flma_sp,r_clay,"OBJECTID_1",fun="max")
r_flma_clay_bool <- r_flma_clay > 0
NAvalue(r_flma_clay_bool) <- 0 
r_flma_clay_bool_fname <- file.path(out_dir,paste0("r_flma_clay_bool_",out_suffix,file_format))
r_flma_clay_bool <- writeRaster(r_flma_clay_bool,filename=r_flma_clay_bool_fname,overwrite=T)
plot(r_flma_clay_bool)

if(gdal_installed==TRUE){
  
  ## Roads
  srcfile <- roads_bool_fname 
  dstfile_roads <- file.path(out_dir,paste("roads_distance_",out_suffix,file_format,sep=""))
  n_values <- "1"
  
  ### Note that gdal_proximity doesn't like when path is too long
  cmd_roads_str <- paste("gdal_proximity.py",basename(srcfile),basename(dstfile_roads),"-values",n_values,sep=" ")
  #cmd_str <- paste("gdal_proximity.py", srcfile, dstfile,sep=" ")
  
  ### Prepare command for FLMA
  
  srcfile <- r_flma_clay_bool_fname 
  dstfile_flma <- file.path(out_dir,paste("r_flma_clay_bool_distance_",out_suffix,file_format,sep=""))
  n_values <- "1"
  
  ### Note that gdal_proximity doesn't like when path is too long
  cmd_flma_str <- paste("gdal_proximity.py",basename(srcfile),basename(dstfile_flma),"-values",n_values,sep=" ")
  #cmd_str <- paste("gdal_proximity.py", srcfile, dstfile,sep=" ")
  
  sys_info<- as.list(Sys.info())$sysname
  
  if(sys_info$sysname=="Windows"){
    shell(cmd_roads_str)
    shell(cmd_flma_str)
  }else{
    system(cmd_roads_str)
    system(cmd_flma_str)
  }
  r_flma_distance <- raster(dstfile_flma)
  r_roads_distance <- raster(dstfile_roads)
  
}else{
  r_roads_distance <- raster(file.path(in_dir,paste("roads_bool_distance_",file_format,sep="")))
  r_flma_distance <- raster(file.path(in_dir_var,paste("r_flma_clay_bool_distance",file_format,sep="")))
}

#Now reverse the distance...

r_flma_distance <- raster(dstfile_flma)
r_roads_distance <- raster(dstfile_roads)

min_val <- cellStats(r_flma_distance,min) 
max_val <- cellStats(r_flma_distance,max)
r <- abs(r_flma_distance  - min_val)/ (max_val - min_val) #no need to inverse...

#Get distance from managed land
#b. Which parts of Clay County contain proximity-to-managed-lands characteristics that would make them more favorable to be used as conservation lands?


##########################
#P3- IDENTIFY POTENTIAL CONSERVATION LANDS IN RELATION TO PARCEL VALUES AND ACREAGES
#GOAL: Create two raster maps showing parcels in Clay County, Florida that have would 
#have higher conservation potential based on parcel value and parcel size (acreage) using 
#a combination of the Extract by Mask, Polygon to Raster, Project, Reclassify, and Weighted Sum tools.

#Now summarize by parcels!!
r_focus_zone1 <- raster("focus_zone1.tif")

#clay_sp_parcels_reg <- spTransform(clay_sp,projection(r_clay))
r_parcels_clay <- rasterize(clay_sp_parcels_reg,r_clay,"OBJECTID_1",fun="min")

r_parcels_clay_min <- rasterize(clay_sp_parcels_reg,r_clay,"OBJECTID_1",fun="min")
#r_parcels_clay <- rasterize(clay_sp_parcels_reg,r_clay,"OBJECTID_1",fun="max")
#r_parcels_clay_max <- rasterize(clay_sp_parcels_reg,r_clay,"OBJECTID_1",fun="max")

#writeRaster(subset(r_bio_factor,"equal_weights"),"bio_factor_equal_weights.tif")
writeRaster(subset(r_bio_factor,1),"bio_factor_equal_weights.tif")

parcels_clay_avg_index <- zonal(subset(r_bio_factor,1), z=r_parcels_clay, fun='mean', digits=0, na.rm=TRUE) 

#This step takes about 18 minutes on my laptop.
#r_parcels_clay_index <- extract(subset(r_bio_factor,1),clay_sp_parcels_reg,sp=T)
parcels_clay_avg_index <- extract(subset(r_bio_factor,1),clay_sp_parcels_reg,fun=mean)
class(parcels_clay_avg_index)
View(parcels_clay_avg_index)

###################### END OF SCRIPT #####################
