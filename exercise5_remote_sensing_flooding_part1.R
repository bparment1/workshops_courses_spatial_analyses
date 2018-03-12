####################################    Flood Mapping Analyses   #######################################
############################  Analyze and map flooding from RITA hurricane  #######################################
#This script performs analyses for the Exercise 5 of the Short Course using reflectance data derived from MODIS.
#The goal is to map flooding from RITA using various reflectance bands from Remote Sensing platforms.
#Additional data is provided including FEMA flood region. 
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/05/2018 
#DATE MODIFIED: 03/12/2018
#Version: 1
#PROJECT: SESYNC and AAG 2018 workshop/Short Course preparation
#TO DO:
#
#COMMIT: generating indices
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
                   r_before,
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

df_modis_band_info <- data.frame("band_name"=NA,"band_number"=NA,"start_wlength"=NA,"end_wlength"=NA)
df_modis_band_info[1,] <- c("Red", 1, 620,670)
df_modis_band_info[2,] <- c("NIR", 2,841,876) 
df_modis_band_info[3,] <- c("Blue",3,459,479)
df_modis_band_info[4,] <- c("Green",4,545,565)
df_modis_band_info[5,] <- c("SWIR1",5,1230,1250)
df_modis_band_info[6,] <- c("SWIR2",6,1628,1652)
df_modis_band_info[7,] <- c("SWIR3",7,2105,2155)

#View(df_modis_band_info)
write.table(df_modis_band_info,
            file=paste0("df_modis_band_info",".txt"),
            sep=",")

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
#View(avg_nlcd)
names(avg_nlcd)

col_ordering <- band_refl_order + 1
plot(as.numeric(avg_nlcd[9,col_ordering]),type="l") #42 evergreen forest
lines(as.numeric(avg_nlcd[6,col_ordering]),type="l") #22 developed,High intensity
class(avg_nlcd)

plot(avg_nlcd$Red,avg_nlcd$NIR)

#############

#Feature space NIR1 and Red
plot(subset(r_before,"Red"),subset(r_before,"NIR"))
plot(subset(r_before,"Green"),subset(r_before,"Red"))
plot(subset(r_before,"SWIR1"),subset(r_before,"NIR"))
plot(subset(r_before,"Red"),subset(r_before,"SWIR1"))

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

#### Feature space Red and Green

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

#### Feature space Blue and Red

#Water:
plot(df_test[df_test$nlcd_2006_RITA==11,c("Blue")],
       df_test[df_test$nlcd_2006_RITA==11,c("Red")],
       col="blue",cex=0.15)

#Forest:
points(df_test[df_test$nlcd_2006_RITA==42,c("Blue")],
     df_test[df_test$nlcd_2006_RITA==42,c("Red")],
     col="green",cex=0.15)

#Urban dense:
points(df_test[df_test$nlcd_2006_RITA==22,c("Blue")],
       df_test[df_test$nlcd_2006_RITA==22,c("Red")],
       col="brown",cex=0.15)

#### Stretch to 0-255 range:

plotRGB(r_before,
        r=1,
        g=4,
        b=3,
        scale=0.6,
        strech="hist")

plotRGB(r_after,
        r=1,
        g=4,
        b=3,
        scale=0.6,
        strech="hist")

r_test### Experiment with threshold:
?colorRamps
col_palette <- colorRampPalette(c("black","blue"))(255)
#colorRampPalette(c("red", "white", "blue"))(255)
plot(subset(r_before,"NIR") < 0.2)
plot(subset(r_before,"Blue"),col=col_palette)

plot(subset(r_before,"NIR"))
plot(subset(r_after,"NIR"))

### THis is suggesting flooding!!!
plot(subset(r_before,"NIR") < 0.2)
plot(subset(r_after,"NIR") < 0.2)

#Compare to actual flooding data


#names(r_before) <- 
#### Add by land cover here:
#Label all pixel with majority vegetation in NIR-Red
#Label all pixel with majority water in NIR-Red
#Label all pixel with majority urban in NIR-Red

## Show water and forest? on a plot

#show expected change vector from forest to water
## Show e

### Make polygons in lake 


############## Generating indices:

#compare indices with FEMA map? Use ROC.

#

####
# Generate flood index?

#1) NDVI = (NIR - Red)/(NIR+Red)
#2) NDWI = (Green - NIR)/(Green + NIR)
#3) MNDWI = Green - SWIR2 / Green + SWIR2
#4) NDWI2 (LSWIB5) =  (NIR - SWIR1)/(NIR + SWIR1)
#5) LSWI (LSWIB5) =  (NIR - SWIR2)/(NIR + SWIR2)
#6) TCWI =  0.10839 * Red+ 0.0912 * NIR +0.5065 * Blue+ 0.404 * Green 
#            - 0.241 * SWIR1- 0.4658 * SWIR2-
#           0.5306 * SWIR3
#7) TCBI = 0.3956 * Red + 0.4718 * NIR +0.3354 * Blue+ 0.3834 * Green
#           + 0.3946 * SWIR1 + 0.3434 * SWIR2+ 0.2964 * SWIR3

names(r_before)
r_date1_NDVI <- (subset(r_before,"NIR") - subset(r_before,"Red")) / (subset(r_before,"NIR") + subset(r_before,"Red"))

#r_date1_NDVI <- subset(r_before,"NIR") - subset(r_before,"Red")/(subset(r_before,"NIR") + subset(r_before,"Red"))

plot(r_date1_NDVI,zlim=c(-1,1))
plot(r_date1_NDVI,zlim=c(-1,1),col=matlab.like(255))

r_date2_NDVI <- (subset(r_after,"NIR") - subset(r_after,"Red"))/ (subset(r_after,"NIR") + subset(r_after,"Red"))

plot(r_date1_NDVI,zlim=c(-1,1),col=matlab.like(255))
plot(r_date2_NDVI,zlim=c(-1,1),col=matlab.like2(255))

### THis is suggesting flooding!!!
plot(r_date1_NDVI < -0.5)
plot(r_date2_NDVI < -0.5)
plot(r_date1_NDVI < -0.1)
plot(r_date2_NDVI < -0.1)

#2) NDWI = (Green - NIR)/(Green + NIR)
#3) MNDWI = Green - SWIR2 / Green + SWIR2

#Modified NDWI MNDWI : Green - SWIR2/(Green + SWIR2)

#According to Gao (1996), NDWI is a good
#indicator for vegetation liquid water content and is less
#sensitive to atmospheric scattering effects than NDVI. In this
#study, MODIS band 6 is used for the NDWI calculation,
#because it is found that MODIS band 6 is sensitive to water
#types and contents (Li et al., 2011), while band 5 is sensitive
#to vegetation liquid water content (Gao, 1996).

names(r_before)
r_date1_MNDWI <- (subset(r_before,"Green") - subset(r_before,"SWIR2")) / (subset(r_before,"Green") + subset(r_before,"SWIR2"))
r_date2_MNDWI <- (subset(r_after,"Green") - subset(r_after,"SWIR2")) / (subset(r_after,"Green") + subset(r_after,"SWIR2"))

plot(r_date1_MNDWI,zlim=c(-1,1))
plot(r_date1_MNDWI,zlim=c(-1,1),col=rev(matlab.like(255)))
plot(r_date2_MNDWI,zlim=c(-1,1))
plot(r_date2_MNDWI,zlim=c(-1,1),col=rev(matlab.like(255)))

### THis is suggesting flooding!!!
plot(r_date1_MNDWI > 0.5)
plot(r_date2_MNDWI > 0.5)
plot(r_date1_MNDWI > 0.1)
plot(r_date2_MNDWI > 0.1)

r_test <- r_date2_MNDWI - r_date1_MNDWI
plot(r_test)

# NRT MODIS
# Other

# Do relationship with flood zone using ROC?
### Generate a map of flooding with MNDWI and compare to FEMA map:

r_date2_flood <- r_date2_MNDWI > 0.1
plot(r_date2_flood)


#reclass in zero/1!!!

df <- data.frame(id=c(1,2), v=c(0,1))
r_ref <- subs(r_ref, df)


ref_test_tb <- crosstab(r_date2_flood,r_ref)

ref_test_tb

###############################################
######## Let's carry out a PCA in T-mode #######

#Correlate long term mean to PC!
cor_mat_layerstats <- layerStats(r_after, 'pearson', na.rm=T)
cor_matrix <- cor_mat_layerstats$`pearson correlation coefficient`
class(cor_matrix)
dim(cor_matrix)
View(cor_matrix)
image(cor_matrix)

pca_mod <-principal(cor_matrix,nfactors=7,rotate="none")
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
loadings_df <- as.data.frame(pca_mod$loadings[,1:7])
#pca_loadings_dz <- zoo(loadings_df,dates_val) #create zoo object from data.frame and date sequence object
#?plot.zoo to find out about zoo time series plotting of indexes
#plot(loadings_df ~ 1:7,
#     #type="b",
#     col=c("blue","red","black","orange","green","purple","brown"),
#     #xlab="time steps",
#     #ylab="PC loadings",
#     #ylim=c(-1,1))
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
### Using predict function: this is recommended for raster imagery!!
# note the use of the 'index' argument
r_pca <- predict(r_before, pca_mod, index=1:7,filename="pc_scores.tif",overwrite=T) # fast
plot(-1*r_pca,y=2,zlim=c(-2,2))
plot(r_pca,y=1,zlim=c(-2,2))
plot(r_pca,y=3,zlim=c(-2,2))

plot(subset(r_pca,1),subset(r_pca,2))
plot(subset(r_pca,2),subset(r_pca,3))

#plot(stack(r_pc1,r_pc2))
#layerStats(r_pc1,r_NDVI_mean )
#cor_pc <- layerStats(stack(r_pc1,r_NDVI_mean),'pearson', na.rm=T)
#cor_pc #PC1 correspond to the average mean by pixel as expected.
#plot(r_pc2)


################### End of Script #########################
