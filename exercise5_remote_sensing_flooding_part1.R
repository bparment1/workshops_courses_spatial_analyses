####################################    Flood Mapping Analyses   #######################################
############################  Analyze and map flooding from RITA hurricane  #######################################
#This script performs analyses for the Exercise 5 of the Short Course using reflectance data derived from MODIS.
#The goal is to map flooding from RITA using various reflectance bands from Remote Sensing platforms.
#Additional data is provided including FEMA flood region. 
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/05/2018 
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

generate_dates_by_step <-function(start_date,end_date,step_date){
  #library(xts) declare out of this function
  #library(zoo)
  #library(lubridate)
  
  st <- as.Date(start_date,format="%Y.%m.%d")
  en <- as.Date(end_date,format="%Y.%m.%d")
  #year_list <-seq(format(st,"%Y"),format(en,"%Y")) #extract year
  year_list <- seq(as.numeric(strftime(st,"%Y")),as.numeric(strftime(en,"%Y"))) #extract year
  
  ll_list <- vector("list",length=length(year_list))
  for (i in 1:length(year_list)){
    if(i==1){
      first_date <-st
    }else{
      first_date<-paste(year_list[[i]],"-01","-01",sep="")
    }
    if(i==length(year_list)){
      last_date <-en
    }else{
      last_date<-paste(year_list[[i]],"-12","-31",sep="")
    }
    #ll <- seq.Date(st, en, by=step)
    ll <- seq.Date(as.Date(first_date), as.Date(last_date), by=step_date)
    ll_list[[i]]<-as.character(ll)
    #paste(yday(ll,)
  }
  
  #
  dates_modis <-as.Date(unlist((ll_list))) 
  
  dates_DOY_modis <- as.character(paste(year(dates_modis),sprintf("%03d", yday(dates_modis)),sep=""))
  dates_obj <- list(dates_modis,dates_DOY_modis)
  names(dates_obj) <- c("dates","doy")  
  return(dates_obj)
}

#function_preprocessing_and_analyses <- "fire_alaska_analyses_preprocessing_functions_03102017.R" #PARAM 1
#function_analyses <- "exercise2_fire_alaska_analyses_functions_03232017.R" #PARAM 1
script_path <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/R_scripts"
#source(file.path(script_path,function_preprocessing_and_analyses)) #source all functions used in this script 1.
#source(file.path(script_path,function_analyses)) #source all functions used in this script 1.

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


##### PART I: DISPLAY AND EXPLORE DATA ##############

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

###
tb_before <- freq(r_before,merge=T)
View(tb_before)
plot(r_after,colNA="black")
plot(r_before,colNA="black")

test <- st_centroid(reg_sf)
test

df_before <- extract(r_before,test)
df_after <- extract(r_after,test)

### We need to rethink the ordering of band:
plot(df_before[2,],type="l")
lines(df_after[2,],col="red")


## Read in table info?
### Need to reorder the bands:
#500m Surface Reflectance Band 1 (620–670 nm)	Reflectance	16-bit signed integer	-28672	-100–16000	0.0001
#500m Surface Reflectance Band 2 (841–876 nm)	Reflectance	16-bit signed integer	-28672	-100–16000	0.0001
#500m Surface Reflectance Band 3 (459–479 nm)	Reflectance	16-bit signed integer	-28672	-100–16000	0.0001
#500m Surface Reflectance Band 4 (545–565 nm)	Reflectance	16-bit signed integer	-28672	-100–16000	0.0001
#500m Surface Reflectance Band 5 (1230–1250 nm)	Reflectance	16-bit signed integer	-28672	-100–16000	0.0001
#500m Surface Reflectance Band 6 (1628–1652 nm)	Reflectance	16-bit signed integer	-28672	-100–16000	0.0001
#500m Surface Reflectance Band 7 (2105–2155 nm)	Reflectance	16-bit signed integer	-28672	-100–16000	0.0001

#SWIR1 (1230–1250 nm), SWIR2 (1628–1652 nm) and SWIR3 (2105–2155 nm).
band_refl_order <- c(3,4,1,2,5,6,7)

names(r_before) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")
names(r_after) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")

plot(df_before[2,band_refl_order],type="l")
lines(df_after[2,band_refl_order],col="red")
plot(df_before[1,band_refl_order],type="l")
lines(df_after[1,band_refl_order],col="red")
plot(df_after[1,band_refl_order],col="red")

###### Now do a extraction for nlcd data


nlcd_2006_filename <- "nlcd_2006_RITA.tif"
nlcd2006_reg <- raster(nlcd_2006_filename)

avg_nlcd <- as.data.frame(zonal(r_before,nlcd2006_reg,fun="mean"))
avg_nlcd <- as.data.frame(avg_nlcd)

avg_nlcd
View(avg_nlcd)
names(avg_nlcd)
plot()

col_ordering <- band_refl_order + 1
plot(avg_nlcd[9,col_ordering],type="l") #42 evergreen forest
lines(avg_nlcd[6,col_ordering],type="l") #22 developed,High intensity
class(avg_nlcd)

plot(avg_nlcd$Red,avg_nlcd$NIR)

#############



#Feature space NIR1 and Red
plot(subset(r_before,"Red"),subset(r_before,"NIR"))

df_test <- as.data.frame(stack(r_before,nlcd2006_reg))

View(df_test)


#Forest:
plot(df_test[df_test$nlcd_2006_RITA==42,c("Green")],
     df_test[df_test$nlcd_2006_RITA==42,c("Red")],
     col="green",cex=0.15)

#Urban: dense
points(df_test[df_test$nlcd_2006_RITA==22,c("Green")],
       df_test[df_test$nlcd_2006_RITA==22,c("Red")],
       col="brown",cex=0.15)

#Water
points(df_test[df_test$nlcd_2006_RITA==11,c("Green")],
       df_test[df_test$nlcd_2006_RITA==11,c("Red")],
       col="blue",cex=0.15)

#### Feature space green and red

#Forest:
plot(df_test[df_test$nlcd_2006_RITA==42,c("Red")],
     df_test[df_test$nlcd_2006_RITA==42,c("NIR")],
     col="green",cex=0.15)

#Urban: dense
points(df_test[df_test$nlcd_2006_RITA==22,c("Red")],
     df_test[df_test$nlcd_2006_RITA==22,c("NIR")],
     col="brown",cex=0.15)

#Water
points(df_test[df_test$nlcd_2006_RITA==11,c("Red")],
       df_test[df_test$nlcd_2006_RITA==11,c("NIR")],
       col="blue",cex=0.15)


#names(r_before) <- 
#### Add by land cover here:
#Label all pixel with majority vegetation in NIR-Red
#Label all pixel with majority water in NIR-Red
#Label all pixel with majority urban in NIR-Red

## Show water and forest? on a plot

#show expected change vector from forest to water
## Show e

############## Generating indices:


####
# Generate flood index?
#The NDWI (Gao, 1996) is a satellite-derived index from the NIR and Short Wave Infrared (SWIR) channels, defined as:
#NDVI = NIR - VIS / (NIR +VIS)

#MODIS has two bands in the SWIR region: band 5 (1230 to
#                                                1250 nm) and band 6 (1.628 to 1.652 m) while band 2
#represents the NIR region; while ABI has two bands in the
#SWIR region, band 4 (1.3705 to 1.3855 m) and band 5 (1.58
#                                                    to 1.64 m).
#
# NDWI = NIR - SWIR / NIR + SWIR
#
names(r_before)
r_date1_NDVI <- subset(r_before,"NIR") - subset(r_before,"Red")/(subset(r_before,"NIR") - subset(r_before,"Red"))
plot(r_date1_NDVI,zlim=c(-1,1))
r_date2_NDVI <- subset(r_after,"NIR") - subset(r_after,"Red")/(subset(r_after,"NIR") - subset(r_after,"Red"))
plot(r_date2_NDVI,zlim=c(-1,1))

#Modified NDWI MNDWI : Green-SWIR1/(NIR+SWIR1)
#
#MMNDWI : Green-Red/(NIR+Red)
#LSWI (Land Surface Water (ndex)): LSWI: NIR -SWIR2/(NIR +SWIR2)

plot(r_NDWI)
tb_NDWI <- freq(r_NDWI)
barplot(tb_NDWI[,2])
barplot(tb_NDWI[,2],xlim=c(-100,100))
histogram(tb_NDWI[,2],xlim=c(-100,100))
hist(tb_NDWI[,2])
View(tb_NDWI)
r_date2_NDWI <- r_after[,a,] - r_after[,,6] / r_after[,,2] - r_after[,,6]  
r_date2_NDWI <- subset(r_after,2) - subset(r_after,6) / (subset(r_after,2) - subset(r_after,6))  

plot(r_date1_NDWI)
plot(r_date2_NDWI)

plot(r_date1_NDWI,zlim=c(-1,1))
plot(r_date1_NDWI,zlim=c(-100,100))

plot(r_date2_NDWI,zlim=c(-1,1))

r_test <- r_date2_NDWI - r_date1_NDWI
plot(r_test,zlim=c(-1,1))
plot(r_date1_NDVI > 1000)

#According to Gao (1996), NDWI is a good
#indicator for vegetation liquid water content and is less
#sensitive to atmospheric scattering effects than NDVI. In this
#study, MODIS band 6 is used for the NDWI calculation,
#because it is found that MODIS band 6 is sensitive to water
#types and contents (Li et al., 2011), while band 5 is sensitive
#to vegetation liquid water content (Gao, 1996).


MDWI
#According to Gao (1996), NDWI is a good
#indicator for vegetation liquid water content and is less
#sensitive to atmospheric scattering effects than NDVI. In this
#study, MODIS band 6 is used for the NDWI calculation,
#because it is found that MODIS band 6 is sensitive to water
#types and contents (Li et al., 2011), while band 5 is sensitive
#to vegetation liquid water content (Gao, 1996).

resample()

# NRT MODIS

# Other

# Do relationship with flood zone using ROC?

###############################################
######## Let's carry out a PCA in T-mode #######

#Correlate long term mean to PC!
cor_mat_layerstats <- layerStats(r_before, 'pearson', na.rm=T)
cor_matrix <- cor_mat_layerstats$`pearson correlation coefficient`
class(cor_matrix)
dim(cor_matrix)
View(cor_matrix)
image(cor_matrix)

pca_mod <-principal(cor_matrix,nfactors=3,rotate="none")
class(pca_mod$loadings)
str(pca_mod$loadings)
plot(pca_mod$loadings[,1][band_refl_order],type="b",
     xlab="time steps",
     ylab="PC loadings",
     ylim=c(-1,1),
     col="blue")
lines(-1*(pca_mod$loadings[,2][band_refl_order]),type="b",col="red")
lines(pca_mod$loadings[,3][band_refl_order],type="b",col="black")
title("Loadings for the first three components using T-mode")

##Make this a time series
loadings_df <- as.data.frame(pca_mod$loadings[,1:3])
pca_loadings_dz <- zoo(loadings_df,dates_val) #create zoo object from data.frame and date sequence object
#?plot.zoo to find out about zoo time series plotting of indexes
plot(pca_loadings_dz,
     type="b",
     plot.type="single",
     col=c("blue","red","black"),
     xlab="time steps",
     ylab="PC loadings",
     ylim=c(-1,1))
title("Loadings for the first three components using T-mode")
names_vals <- c("pc1","pc2","pc3")
legend("topright",legend=names_vals,
       pt.cex=0.8,cex=1.1,col=c("blue","red","black"),
       lty=c(1,1), # set legend symbol as lines
       pch=1, #add circle symbol to line
       lwd=c(1,1),bty="n")

## Add scree plot
plot(pca_mod$values,main="Scree plot: Variance explained",type="b")

### Generate scores from eigenvectors
## Do it two different ways:

#By converting data and working with matrix:
df_NDVI_ts <- as(r_NDVI_ts,"SpatialPointsDataFrame")
df_NDVI_ts <- as.data.frame(df_NDVI_ts) #convert to data.frame because predict works on data.frame
names(df_NDVI_ts) <- c(paste0("pc_",seq(1,23,1)),"x","y")
names(df_NDVI_ts)

pca_all <- as.data.frame(predict(pca_mod,df_NDVI_ts[,1:23])) ## Apply model object on the data.frame
#pca_all <-predict(pca_mod,df_NDVI_ts) ## Apply model object on the data.frame

coordinates(pca_all) <- df_NDVI_ts[,c("x","y")] #Assign coordinates
class(pca_all) #Check type of class

raster_name <- paste0("pc1_NDVI_2005.tif") #output raster name for component 1
r_pc1<-rasterize(pca_all,r_NDVI_ts,"PC1",fun=min,overwrite=TRUE,
                 filename=raster_name)
raster_name <- paste0("pc2_NDVI_2005.tif") #output raster name for component 2
r_pc2<-rasterize(pca_all,r_NDVI_ts,"PC2",fun=min,overwrite=TRUE,
                 filename=raster_name)

### Using predict function: this is recommended for raster imagery!!
# note the use of the 'index' argument
r_pca <- predict(r_before, pca_mod, index=1:3,filename="pc_scores.tif",overwrite=T) # fast
plot(-1*r_pca,y=2,zlim=c(-2,2))
plot(r_pca,y=1,zlim=c(-2,2))
plot(r_pca,y=3,zlim=c(-2,2))

plot(subset(r_pca,1),subset(r_pca,2))
plot(subset(r_pca,2),subset(r_pca,3))

plot(stack(r_pc1,r_pc2))
#layerStats(r_pc1,r_NDVI_mean )
cor_pc <- layerStats(stack(r_pc1,r_NDVI_mean),'pearson', na.rm=T)
cor_pc #PC1 correspond to the average mean by pixel as expected.
plot(r_pc2)


####### PART II: Change analyses via image differencing and ratioing ##########

## Studying change by image differencing
#look for NDVI 2002 and 2009 and select each average

r_NDVI_avg_2002 <- subset(r_var,4) #select layer 4 from raster stack
r_NDVI_avg_2009 <- subset(r_var,5)

r_diff_NDVI <- r_NDVI_avg_2009 - r_NDVI_avg_2002 #if negative Higher NDVI in 2002, hence decrease in NDVI over 2002-2009
r_ratio_NDVI <- r_NDVI_avg_2009/r_NDVI_avg_2002

plot(r_ratio_NDVI,col=matlab.like(255),zlim=c(0.5,1.5)) #zlim to control the range displayed
plot(r_diff_NDVI,col=matlab.like(255),zlim=c(-0.5,0.5))

### Quick histogram
hist(r_diff_NDVI)
##Adjust the bins

hist_bins <- seq(-1,1,by=0.05)
hist(r_diff_NDVI,breaks=hist_bins)
hist(r_diff_NDVI,breaks=hist_bins,xlim=c(-0.3,0.3))

plot(r_diff_NDVI,col=matlab.like(255),zlim=c(-0.2,0.2))

## Threshold your inputs
plot(r_diff_NDVI < -0.2)
plot(r_diff_NDVI < -0.1)
plot(r_diff_NDVI > 0.1)

##Standardize images and define change
#r <- as.data.frame(cellStats(x,mean))
val_diff_mean <- cellStats(r_diff_NDVI,mean)
val_diff_sd <- cellStats(r_diff_NDVI,sd)

r_diff_NDVI_standardized <- (r_diff_NDVI - val_diff_mean)/val_diff_sd
plot(r_diff_NDVI_standardized,col=matlab.like(255))#
plot(r_diff_NDVI_standardized,col=matlab.like(255),zlim=c(-10,5))#

hist_bins <- seq(-15,15,by=0.5)
hist(r_diff_NDVI_standardized,breaks=hist_bins)
hist(r_diff_NDVI_standardized,breaks=hist_bins,xlim=c(-5,5)) # zoom in

hist(r_diff_NDVI_standardized)
#shapiro.test(values(r_diff_NDVI_standardized))
#qqnorm(values(r_diff_NDVI_standardized))

r_change_NDVI_pos <-  r_diff_NDVI_standardized > 1.96
r_change_NDVI_neg <-  r_diff_NDVI_standardized < -1.96
plot(r_change_NDVI_pos)
plot(r_change_NDVI_neg)

writeRaster(r_change_NDVI_pos,"r_change_NDVI_pos_196.tif",overwrite=T)
writeRaster(r_change_NDVI_neg,"r_change_NDVI_neg_196.tif",overwrite=T)

####### PART III: Change analyses by comparing averages in fire polygons ##########

##Extract values by fire scars

r_fire_poly <- subset(r_var,6)
plot(r_fire_poly)
mean_fire_poly_tb <- zonal(stack(r_NDVI_avg_2002,r_NDVI_avg_2009),r_fire_poly,fun="mean") #mean square error
mean_fire_poly_df <- as.data.frame(t(mean_fire_poly_tb[-1,-1]))

names(mean_fire_poly_df) <- c("fire_pol1","fire_pol2","fire_pol3")
#write.table(as.data.frame(mean_wind_zones_tb),"mean_wind_zones_tb.txt",sep=",")
mean_fire_poly_df$year <- c(2002,2009)

## Plot the average NDVI by burn scars
#Note that the decrease in NDVI varies according to the burned intensity
plot(fire_pol1 ~year,data=mean_fire_poly_df,type="b",
     ylim=c(0.1,0.4),ylab="Average NDVI")
lines(fire_pol2 ~year,data=mean_fire_poly_df,type="b",col="red")
lines(fire_pol3 ~year,data=mean_fire_poly_df,type="b",col="blue")
title("Average NDVI for fire polygons")
names_vals <- c("pol1","pol2","pol3")
legend("topright",legend=names_vals,
       pt.cex=0.8,cex=1.1,col=c("black","red","blue"),
       lty=c(1,1), # set legend symbol as lines
       pch=1, #add circle symbol to line
       lwd=c(1,1),bty="n")

### Your turn: use albedo images to define image of changes, what do you think is a good threshold for change?

#1.Select Albedo images
#2.Peform differencing, and standardization
#3.Generate Image of changes
#4.Compute average by polygons of fire and compare to NDVI.

######### PART IV: time series analyses #################

#1. Extract time series from fire polygon
#2. Visualize time series using zoo object
#3. Compute ACF
#4. Perform PCA
#5. Generate movie using animate

# Using LST time series perform similar analysis
#r_NDVI_mean <- stackApply(r_NDVI_ts, indices=rep(1,23), fun=mean,na.rm=T) # works too but slower
r_NDVI_mean <- mean(r_NDVI_ts, na.rm=TRUE) # mean by pixel
projection(r_NDVI_ts) <- CRS_reg

mean_fire_NDVI_ts_tb <- zonal(stack(r_NDVI_ts),r_fire_poly,fun="mean") #mean square error
mean_fire_NDVI_ts_df <- as.data.frame(t(mean_fire_NDVI_ts_tb[-1,-1]))
names(mean_fire_NDVI_ts_df) <- c("fire_pol1","fire_pol2","fire_pol3")

## Make a time series object
NDVI_fire_dat_dz <- zoo(mean_fire_NDVI_ts_df,dates_val) #create zoo object from data.frame and date sequence object
plot(NDVI_fire_dat_dz,main="Times series of average NDVI in fire polygons for 2005",type="b")
#Note the sudden decrease in NDVI. It is related to the fire event. It is an average so the decresease may not be
#strong. Let's examine pixel level temporal profiles now.

#Read in fire plygon
fire_poly_sp <-readOGR(dsn=in_dir_var,sub(".shp","",fire_poly_shp_fname))
proj4string(fire_poly_sp) <- CRS_reg
plot(r_diff_NDVI_standardized,ext=extent(fire_poly_sp))
plot(fire_poly_sp,add=T)
spplot(fire_poly_sp)

r_stack <- stack(r_NDVI_ts,r_fire_poly)
inMemory(r_stack) #check that it is not stored in memory
poly_fire_spdf <- as(crop(r_stack,extent(fire_poly_sp)),"SpatialPointsDataFrame")
poly_fire_spdf <- rename(poly_fire_spdf,c("r_OVERLAY_ID_83_399_144_TEST_BURNT_83_144_399_reclassed"="fire_id"))
table(poly_fire_spdf$fire_id)
barplot(table(poly_fire_spdf$fire_id), main="Pixel count by fire polygon (1,2,3) in focus area")
       
## labels
labelat <- c(0, 1, 2,3)

labeltext <- c("background","fire1","fire2","fire3")

spplot(poly_fire_spdf,"fire_id",        
        main="Pixel count in small area",
        col.regions=c("white","blue","yellow","red")#,
        #colorkey = list(width=1,
        #                space="right",
        #                tick.number=5,
        #                labels = list(at = labelat,labeltext)
        #               )
)

#let's plot polygons one since it has a variety of pixels and was heavily affected

#pix_fire_poly1_df  <- as.data.frame(t(subset(poly_fire_spdf,fire_id==1)))
pix_fire_poly1_df  <- as.data.frame(t(as.data.frame(subset(poly_fire_spdf,fire_id==1))))
names(pix_fire_poly1_df) <- paste0("pix_",1:ncol(pix_fire_poly1_df))
NDVI_fire_dat_dz <- zoo(mean_fire_NDVI_ts_df,dates_val) #create zoo object from data.frame and date sequence object
dim(NDVI_fire_dat_dz)

pix_fire_poly1_dz <- zoo(pix_fire_poly1_df,dates_val) #create zoo object from data.frame and date sequence object
#colnames(pix_fire_poly1_dz) <- names(pix_fire_poly1_df)
##Explore Time series
plot(pix_fire_poly1_dz[,1000:1010])
acf(pix_fire_poly1_dz[,1000], type = "correlation") #burnt pixel
acf(pix_fire_poly1_dz[,1007], type = "correlation") #Note the difference in the acf for unburnt pixel


## Use animate function
#animate(r_NDVI_ts, pause=0.25, n=1)

#### Your turn: 
#Compute the average using the ecoregions as zonal areas for PC1 and PC2
#Plot averages in the PC1-PC2. Are these variables useful to separate the various ecoregions?
#Correlate average times series for PC2 to the loadings. What does it tell you?

################### End of Script #########################
