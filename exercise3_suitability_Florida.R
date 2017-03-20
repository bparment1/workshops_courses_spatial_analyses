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
#DATE MODIFIED: 03/19/2017
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

plot(r_strat_hab, "startegic habitat")
plot(reg_counties_sp,add=T)
plot(r_priority_wet_hab)
#plot(r_habitat,add=T)
plot(roads_sp,add=T)


##### Before starting the production of sustainability let's check the projection, resolution for each layer

#raster layers:

list_raster <- c(r_strat_hab,r_priority_wet_hab,r_habitat,r_bio_hotspot)

lapply(list_raster,function(x){res(x)})
lapply(list_raster,function(x){projection(x)})
lapply(list_raster,function(x){extent(x)})

## let's use the resolution 55x55 as the reference


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
r_strat_hab <- raster(file.path(in_dir_var,strat_hab_fname))

r_strat_hab




### STEP 2: Identify Lands With High Native Biodiversity: based on species count 

## Select clay county
clay_county_sp <- subset(reg_counties_sp,NAME=="CLAY")

plot(r_strat_hab)
plot(clay_county_sp,border="red")
r_bio_hotspot <- raster(file.path(in_dir_var,biodiversity_hotspot_fname))

## Crop bio raster
r_test <- crop(r_bio_hotspot,clay_county_sp)
plot(r_test)
plot(clay_county_sp,border="red",add=T)

r_clay <- rasterize(clay_county_sp,r_test)
freq(r_clay)
r_bio_clay <- mask(r_test,r_clay)
plot(r_bio_clay)
freq(r_bio_clay)

#Criteria for value assignment: Cells with 7 or more focal species or an actual species occurrence record
#location were assigned a value of 9; a value of 8 was assigned to cells with 5–8 focal species; 7 was
#assigned to cells with 3–4 focal species; all other cells were assigned a value of 1.

# reclassify the values into three groups 
# all values >= 0 and <= 0.25 become 1, etc.
#m <- c(0, 0.25, 1,  
#       0.25, 0.5, 2,  
#       0.5, 1, 3)

#rclmat <- matrix(m, ncol=3, byrow=TRUE)
#rc <- reclassify(r, rclmat)

m <- c(9, 1000, 9,  
       5, 8, 8,  
       3, 4, 7,  
       1, 2,1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

rc <- reclassify(r_bio_clay, rclmat)
plot(rc)

### STEP 3: Wetland priority 

#Input data layer: Priority Wetland Habitats
#Criteria for value assignment: Values were assigned based on the number of focal species present in
#each cell. The value of 9 was assigned to 10–12 wetland focal species, 8 was assigned to 7–9 wetland focal
#species, 7 was assigned to 4–6 wetland focal species and 4–6 upland focal species, 6 was assigned to
#1–3 wetland or upland focal species. The value 1 was assigned to all other cells.
#Rationale for value assignment: The better the habitat for focal wetland species, the higher the priority.
#Output: Wetland Biodiversity SUA (CG1O11SO111)

#r_priority_wet_hab <- raster(file.path(in_dir_var,priority_wet_habitats_fname))

plot(r_priority_wet_hab)
#check projection
projection(r_priority_wet_hab)

## Crop Wetland priority raster
r_priority_wet_hab_w <- crop(r_priority_wet_hab,clay_county_sp)
plot(r_priority_wet_hab_w)
plot(clay_county_sp,border="red",add=T)

##Use raster of Clay county definining the study area to mask pixels
plot(r_clay)
dim(r_clay)
dim(r_priority_wet_hab_w)

r_priority_wet_hab_clay <- mask(r_priority_wet_hab_w,r_clay) ## Does not work!! because resolution don't match
#match resolution:
#projection(r_priority_wet_hab_w)==projection(r_clay)
#res(r_priority_wet_hab_w)==res(r_clay)

## Find about resample
?raster::resample
r_priority_wet_hab_clay <- resample(r_priority_wet_hab_w,r_clay, method='bilinear')

r_priority_wet_hab_clay <- mask(r_priority_wet_hab_clay,r_clay) ## Does not work!! because resolution don't match
plot(r_priority_wet_hab_clay)

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

rc_priority_wet_hab_clay <- reclassify(r_priority_wet_hab_clay, rclmat)
plot(rc_priority_wet_hab_clay)
freq_tb <- freq(rc_priority_wet_hab_clay)


### STEP 4: Combine all the three layers with weigthed sum




###################### END OF SCRIPT #####################
