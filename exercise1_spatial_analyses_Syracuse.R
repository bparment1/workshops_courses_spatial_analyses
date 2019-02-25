####################################    Spatial Analyses: SYRACUSE   #######################################
#######################################  Analyse data from Census #######################################
#
# This script performs analyses for the Exercise 1 of the Geospatial Analysis course created for the AAG2017.
# The goal is to assess land cover change using two land cover maps in the Houston areas.
# This script performs basic analyses for the Exercise 1 using Census data for Syracuse.
# The overall goal is to explore spatial autocorrelation and aggregation of units of analyses.     

#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/21/2017 
#DATE MODIFIED: 02/15/2019
#Version: 1
#PROJECT: AAG 2019 Geospatial workshop and Sesync Geopstial Data Analyses course, Geocompuation Yale
#TO DO:
#
#COMMIT: Geospatial Data Analysis Course
#
#################################################################################################

###Loading R library and packages                                                      
library(gstat) #spatial interpolation and kriging methods
library(sp) # spatial/geographic objects and functions
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
library(sf) # spatial ojbects simple feature model implementation OGC
#library(gstat)
#library(spacetime)

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

#####  Parameters and argument set up ###########

##Change path to your location where you downlaoded the material
in_dir_var <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2019_geospatial_workshop/Exercise_1/data"
out_dir <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2019_geospatial_workshop/Exercise_1/outputs"

file_format <- ".tif" # raster file format extension
NA_flag_val <- -9999 #NA value, default 
out_suffix <-"exercise1_02122019" #output suffix for the files and output folder #PARAM 8
create_out_dir_param=TRUE #if true, a new output directory is generated

################################   START SCRIPT   ###################################

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

### PART I: EXPLORE DATA: READ AND DISPLAY #######

ct_2000_fname <- "ct_00.shp" # CT_00: Cencus Tracts 2000
bg_2000_fname <- "bg_00.shp" # BG_00: Census Blockgroups 2000
bk_2000_fname <- "bk_00.shp" # BK_00: Census Blocks 2000

census_table_fname <- "census.csv" #contains data from census to be linked
soil_PB_table_fname <- "Soil_PB.csv" #same as census table

metals_table_fname <- "SYR_metals.xlsx" #contains metals data to be linked

ct_2000_sf <- st_read(file.path(in_dir_var,ct_2000_fname)) #read in shapefile
bg_2000_sf <- st_read(file.path(in_dir_var,bg_2000_fname))
bk_2000_sf <- st_read(file.path(in_dir_var,bk_2000_fname))

census_syr_df <- read.table(file.path(in_dir_var,census_table_fname),sep=",",header=T) #read in textfile
metals_df <-read_excel( file.path(in_dir_var,metals_table_fname),1) #use function from readxl

#Soil lead samples: UTM z18 coordinates
soil_PB_df <- read.table(file.path(in_dir_var,soil_PB_table_fname),sep=",",header=T) #point locations

dim(ct_2000_sf) #57 spatial entities corresponding to census tracks
dim(metals_df) #57 entities
dim(bg_2000_sf) #147 spatial entities corresponding to blockgroups
dim(bk_2000_sf) #2025 spatial entities corresponding to blocks 

###PRODUCE MAPS OF 2000 Population #########

#First need to link it data attribute table to sf features:

#names(bg_2000_sp) #missing census data
names(bg_2000_sf) #missing census data

names(census_syr_df)
#key is "TRACT" but with a different format.
#First fix the format
#head(bg_2000_sp)
head(bg_2000_sf)

head(census_syr_df$BKG_KEY)
class(census_syr_df$BKG_KEY)
class(bg_2000_sf$BKG_KEY)
bg_2000_sf$BKG_KEY <- as.numeric(as.character(bg_2000_sf$BKG_KEY)) 

#as.numeric(as.character(ct_2000_sp$TRACT))
class(ct_2000_sf$TRACT)
class(bg_2000_sf$TRACT)

#ct_2000_sp$TRACT <- as.numeric(as.character(ct_2000_sp$TRACT)) 
ct_2000_sf$TRACT <- as.numeric(as.character(ct_2000_sf$TRACT)) 

## Join based on common key id
#bg_2000_sp <- merge(bg_2000_sp,census_syr_df,by="BKG_KEY") #Join 
bg_2000_sf <- merge(bg_2000_sf,census_syr_df,by="BKG_KEY") #Join 

#Plot the spatial object
#spplot(bg_2000_sp,"POP2000",main="POP2000") #quick visualization of population 
plot(bg_2000_sf["POP2000"],main="POP2000")

##Aggregate data from block group to census

### Summarize by census track
#census_2000_sp <- aggregate(bg_2000_sp , by="TRACT",FUN=sum)
bg_2000_sf$TRACT <- as.numeric(as.character(bg_2000_sf$TRACT)) 

## 
census_2000_sf <- aggregate(bg_2000_sf , by=list(bg_2000_sf$TRACT),FUN=sum)
dim(census_2000_sf)
class(census_2000_sf)

### Check if the new geometry of entities is the same as census
#plot(census_2000_sp)
plot(census_2000_sf$geometry)

plot(ct_2000_sf,border="red",pal=NA,add=T)
nrow(census_2000_sf)==nrow(ct_2000_sf)

#df_summary_by_census <- aggregate(. ~ TRACT, bg_2000_sf , FUN=sum) #aggregate all variables from the data.frame
#df_summary_by_census <- aggregate(. ~ TRACT, bg_2000_sf , FUN=sum) #aggregate all variables from the data.frame

##Join by key table id:
dim(ct_2000_sf)
dim(census_2000_sf)
## merge two sf not available

census_2000_df <- as.data.frame(census_2000_sf)
dim(census_2000_df)
class(census_2000_df)
head(census_2000_df)
nrow(census_2000_df)

#test <- merge(ct_2000_sf,census_2000_df,by="TRACT")

#test <- st_equals(ct_2000_sf,census_2000_sf#by="TRACT"))

#dim(ct_2000_sp)

#dim(ct_2000_sp)
#ct_2000_sp <- merge(ct_2000_sp,df_summary_by_census,by="TRACT")
#dim(ct_2000_sp)

#save as sp and text table
#write.table(file.path(out_dir,)

### Do another map with different class intervals: 

title_str <- "Population by census tract in 2000"
#range(ct_2000_sp$POP2000)
range(census_2000_sf$POP2000)


## 9 classes with fixed and constant break
break_seq <- seq(0,9000,1000)
#breaks.qt <- classIntervals(ct_2000_sp$POP2000, n=length(break_seq), 
#                            style="fixed", fixedBreaks=break_seq, intervalClosure='right')
breaks.qt <- classIntervals(census_2000_sf$POP2000, n=length(break_seq), 
                            style="fixed", fixedBreaks=break_seq, intervalClosure='right')

n_classes <- length(breaks.qt$brks) -1 
col_palette <- matlab.like(n_classes)
#col_palette <- matlab.like(256)

## generate plot using sp function:
#p_plot_pop2000_ct <- spplot(ct_2000_sp,
#                            "POP2000",
#                            col="transparent", #transprent color boundaries for polygons
#                            col.regions = col_palette ,
#                           main=title_str,
#                            at = breaks.qt$brks)
#print(p_plot_pop2000_ct)

plot(census_2000_sf["POP2000"],
                            #col=col_palette, #transprent color boundaries for polygons
                            pal = col_palette ,
                            main=title_str,
                            at = breaks.qt$brks)
### Another map with different class intervals

n_classes <- 5
#n_breaks <- n_classes+1
breaks.qt <- classIntervals(census_2000_sf$POP2000, n = n_classes, style = "quantile", intervalClosure = "right")
col_palette <- matlab.like(n_classes)
length(breaks.qt$brks)

#Fix error here
#plot(census_2000_sf["POP2000"],
#                            #col="transparent", #transprent color boundaries for polygons
#                            pal = col_palette,
#                            main=title_str,
#                            at = breaks.qt$brks)
#length(breaks.qt$brks)
#length(col_palette)
#print(p_plot_pop2000_ct)

####### PART II: SPATIAL QUERY #############

## Join metals to census track 
## Join lead (pb) measurements to census tracks

#metals_df <- read.xls(file.path(in_dir_var,metals_table_fname),sep=",",header=T)
metals_df <-read_excel( file.path(in_dir_var,metals_table_fname),1) #use function from readxl

head(soil_PB_df)
head(metals_df)

##This suggests matching to the following spatial entities
nrow(metals_df)==nrow(ct_2000_sf)
nrow(metals_df)==nrow(census_2000_sf)

census_2000_sf$TRACT <- census_2000_sf$Group.1
#nrow(soil_PB_df)==nrow(bg_2000_sp)

#dim(bg_2000_sp)
metals_df$TRACT <- metals_df$ID
#census_metals_sf <- merge(ct_2000_sf,metals_df,by="TRACT")
census_metals_sf <- merge(census_2000_sf,metals_df,by="TRACT")
#names(census_2000_sf$TRACT)
#class(census_2000_sf$TRACT)
plot(census_2000_sf$geometry)
plot(census_metals_sf$geometry,col="red")

########processing lead data
### Now let's plot lead data 
#Soil lead samples: UTM z18 coordinates

soil_PB_df <- read.table(file.path(in_dir_var,soil_PB_table_fname),sep=",",header=T) #point locations

st_crs(census_metals_sf) #
names(soil_PB_df)
names(soil_PB_df) <- c("x","y","ID","ppm") 
soil_PB_sp <- soil_PB_df
class(soil_PB_df)

###Create a sf from coordinates
epsg_code <- st_crs(census_metals_sf)$epsg
#df.SP <- st_as_sf(soil_PB_df, coords = c("LONG", "LAT"), crs = 4326)
soil_PB_sf <- st_as_sf(soil_PB_df, coords = c("x", "y"), crs = epsg_code)
plot(soil_PB_sf)

#coordinates(soil_PB_sp) <- soil_PB_sp[,c("x","y")]
class(soil_PB_sf)
#proj4string(soil_PB_sp) <- proj4string(census_metals_sp)
dim(soil_PB_sf)
#soil_PB_sp <- soil_PB_sp[,c("ID","ppm","x","y")]

plot(census_metals_sf$geometry,border='blue')
plot(soil_PB_sf$geometry,add=T,pch='+')

###### Spatial query: associate points of pb measurements to each census tract
### Get the ID and 

#soil_tract_id_df <- over(soil_PB_sp,census_2000_sp,fn=mean)
soil_PB_join_sf <- st_join(x=soil_PB_sf, y=census_2000_sf, join = st_within)
#test <- aggregate(x=soil_PB_sf,by=census_2000_sf,FUN=mean,join=st_within)

#head(test)

#test <- st_intersection(soil_PB_sf,census_2000_sf)
census_pb_avg <- aggregate(ppm ~ TRACT,soil_PB_join_sf,FUN=mean)
#census_pb_avg <- aggregate(ppm ~ TRACT,(soil_PB_sp),FUN=mean)
#census_pb_avg <- rename(census_pb_avg,c("ppm"="pb_ppm"))

#soil_PB_sp <- intersect(soil_PB_sp,census_2000_sp)
#head(soil_PB_sp$ID)==head(soil_PB_sp$ID)
#names(soil_PB_sp)
names(census_pb_avg)
#soil_PB_sp <- rename(soil_PB_sp, c("d"="TRACT")) #from package plyr

##Now join
#census_metals_pb_sp <- merge(census_metals_sp,census_pb_avg,by="TRACT")
census_metals_pb_sf <- merge(census_metals_sf,census_pb_avg,by="TRACT")

### write out final table and shapefile
plot(census_metals_pb_sf['ppm'])
outfile <-paste("census_metals_pb_sp","_",
               out_suffix,".shp",sep="")

#writeOGR(census_metals_pb_sp,dsn= out_dir,layer= outfile, driver="ESRI Shapefile",overwrite_layer=TRUE)

st_write(census_metals_pb_sf,
         file.path(out_dir,outfile),
         #driver="ESRI Shapefile",
         delete_dsn=TRUE)

outfile_df_name <- file.path(out_dir,paste0(outfile,".txt"))
write.table(as.data.frame(census_metals_pb_sf),file=outfile_df_name,sep=",")

####################  PART III: RASTER FROM KRIGING   ######################
#Generating raster lead surface from point and comparing aggregation ###################

#Now generate a raster image to create grid of cell for kriging
extent_reg <- extent(census_metals_pb_sf)
plot(extent_reg)
#plot(census_metals_pb_sp,add=T)
plot(census_metals_pb_sf$geometry,add=T)

extent_matrix <- as.matrix(extent_reg)
extent_matrix

x_length_reg <- extent_matrix[1,2] - extent_matrix[1,1] 
y_length_reg <- extent_matrix[2,2] - extent_matrix[2,1] 

print(c(x_length_reg,y_length_reg))

## Based on the size of the extent, let's set the size for our new raster layer: 
#we don't want too fine as resolution: let's do 100m, it will keep the grid small

res_val <- 100
r = raster(ext=extent_reg, res=res_val)
dim(r)
values(r) <- 1:ncell(r) # Assign values to raster, ID for each pixel
#assign projection system


projection(r) <- st_crs(census_metals_pb_sf)$proj4string

######Visualize the data first

plot(r)
#generate grid from raster as poly for visualization
r_poly<- rasterToPolygons(r)
plot(extent_reg,add=T,col="red")
plot(census_metals_pb_sf['ppm'],border="blue",add=T)
### Let's show the grid first
plot(r_poly,add=T)

## Transform the raster layer into a sp Grid object for kriging
r_sgdf <- as(r, 'SpatialGridDataFrame')
class(r_sgdf)

## Generate and plot a sample variogram from lead data
#soil_PB_sp <- 
v_ppm <- variogram(ppm ~ 1,
    locations = ~ x + y,
    data = soil_PB_df)
#v_ppm <- variogram(ppm ~ 1,soil_PB_sp)
plot(v_ppm)

## Fit a variogram model from lead data

v_ppm_fit <- fit.variogram(v_ppm,model=vgm(1,"Sph",900,1))
plot(v_ppm,v_ppm_fit)

##Generate a kriging surface using data and modeled variogram: this may take more than 3 minutes
soil_PB_sp <- as(soil_PB_sf,"Spatial") # convert to spatial object
#pred_ppm_xy <- g$krige(
#  ppm ~ 1,
#  locations = ~ x + y,
# data = lead_xy,
#  newdata = pred_ppm_xy,
#  model = v_ppm_fit)

ppm_lead_spg <- krige(ppm ~ 1, 
                      soil_PB_sp, 
                      r_sgdf, 
                      model=v_ppm_fit)

class(ppm_lead_spg)
r_lead <- raster(ppm_lead_spg)
rm(ppm_lead_spg) #remove grid object from memory
r_lead #examine new layer

col_palette <- matlab.like(256)
plot(r_lead,col=col_palette)
plot(census_metals_pb_sf$geometry,border="blue",add=T)

## Save raster layers produced from kriging
raster_name <- file.path(out_dir,paste0("r_lead",out_suffix,file_format))
writeRaster(r_lead,filename = raster_name,overwrite=T)

#### Comparison of aggregations ###
## Compare values from averages from kriging surface and averages from block groups

census_metals_pb_sp <- as(census_metals_pb_sf,'Spatial') #convert sf to sp object
census_lead_sp <- extract(r_lead,census_metals_pb_sp,sp=T,fun=mean) #extract average values by census track
spplot(census_metals_pb_sp,"ppm",col.regions=col_palette,main="Averaged from blockgroups") #
spplot(census_lead_sp,"var1.pred",col.regions=col_palette,main="Averaged from kriging ") 

census_lead_sp$diff <- census_metals_pb_sp$ppm - census_lead_sp$var1.pred #comparing the averages
hist(census_lead_sp$diff)
spplot(census_lead_sp,"diff",col.regions=col_palette,main="Difference in averages")

##### PART IV: Spatial autocorrelation and regression #############
## Examine spatial autocorrelation
#Examine the relationship between metals, Pb and vulnerable populations in Syracuse

list_nb <- poly2nb(census_lead_sp) #generate neighbours based on polygons
summary(list_nb)
plot(census_lead_sp,border="blue")
plot.nb(list_nb,coordinates(census_lead_sp),add=T)

#generate weights
#nb2listw
list_w <- nb2listw(list_nb, glist=NULL, style="W", zero.policy=NULL) #use row standardized
can.be.simmed(list_w)
summary(list_w)

## Compute Moran's I and display it
moran(census_lead_sp$ppm,list_w,n=nrow(census_lead_sp), Szero(list_w))
moran.plot(census_lead_sp$ppm, list_w,
           labels=as.character(census_lead_sp$TRACT), pch=19)

##### Now do a spatial regression

## Is there are relationship between minorities and high level of lead?
# As a way to explore use,  perc_hispa as explanatory variable

#linear model without taking into account spatial autocorrelation
mod_lm <- lm(ppm ~ perc_hispa, data=census_lead_sp)
#autoregressive model
mod_lag <- lagsarlm(ppm ~ perc_hispa, data=census_lead_sp, list_w, tol.solve=1.0e-30)

### Checking for autocorrelation in residuals
moran.test(mod_lm$residuals,list_w)
moran.test(mod_lag$residuals,list_w) #Note that Moran'sI is close to zero in the lag model

#### Compare Moran's I from raster to Moran's I from polygon sp
# Rook's case
f <- matrix(c(0,1,0,1,0,1,0,1,0), nrow=3)
Moran(r_lead, f) 
r_moran <- MoranLocal(r_lead)
plot(r_moran) # hotspots?

######################### END OF SCRIPT ##################################


