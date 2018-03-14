####################################    Flood Mapping Analyses   #######################################
############################  Analyze and map flooding from RITA hurricane  #######################################
#This script performs analyses for the Exercise 6 of the Short Course using reflectance data derived from MODIS.
#The goal is to map flooding from RITA using various reflectance bands from Remote Sensing platforms.
#Additional data is provided including FEMA flood region. 
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/13/2018 
#DATE MODIFIED: 03/14/2018
#Version: 1
#PROJECT: SESYNC and AAG 2018 workshop/Short Course preparation
#TO DO:
#
#COMMIT: setting up input data
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

#function_preprocessing_and_analyses <- "fire_alaska_analyses_preprocessing_functions_03102017.R" #PARAM 1
#function_analyses <- "exercise2_fire_alaska_analyses_functions_03232017.R" #PARAM 1
script_path <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/R_scripts"
#source(file.path(script_path,function_preprocessing_and_analyses)) #source all functions used in this script 1.
#source(file.path(script_path,function_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#in_dir_reflectance <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/Exercise_6/data/reflectance_RITA"
in_dir_var <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/Exercise_6/data/"
out_dir <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/Exercise_6/outputs"
infile_reg_outline <- "/nfs/bparmentier-data/Data/Space_beats_time/Data/data_RITA_reflectance/revised_area_Rita/new_strata_rita_10282017.shp"
#region coordinate reference system
#http://spatialreference.org/ref/epsg/nad83-texas-state-mapping-system/proj4/
CRS_reg <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 
file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"exercise6_03132018" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
date_event <- ""
#ARG4
method_proj_val <- "bilinear" # "ngb"

#ARG9
#local raster name defining resolution, extent
ref_rast_name <- "/nfs/bparmentier-data/Data/Space_beats_time/Data/data_RITA_reflectance/revised_area_Rita/r_ref_Houston_RITA.tif"

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

plot(nlcd2006_reg_RITA==90)
plot(nlcd2006_reg_RITA==95)
plot(nlcd2006_reg_RITA==11)

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

nlcd_2006_filename <- file.path(in_dir_var,"nlcd_2006_RITA.tif")
nlcd2006_reg <- raster(nlcd_2006_filename)

plot(nlcd2006_reg)

avg_nlcd <- as.data.frame(zonal(r_after,nlcd2006_reg))
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

# NRT MODIS
# Other

##### plot feature space:

df_raster_val <- as.data.frame(stack(r_after,r_pca,nlcd2006_reg))

### Now do a unsupervised
## do a supervised
## training and testing... 


# Do relationship with flood zone using ROC?
### Generate a map of flooding with MNDWI and compare to FEMA map:

r_date2_flood <- r_date2_MNDWI > 0.1
plot(r_date2_flood)

#reclass in zero/1!!!

df <- data.frame(id=c(1,2), v=c(0,1))
r_ref_test <- subs(r_ref, df)
plot(r_ref_test)
#plot(r_ref)
ref_test_tb <- crosstab(r_date2_flood,r_ref_test)

## Compute Jaccard Index:

ref_test_tb$Freq[5]/(sum(ref_test_tb[ref_test_tb$Var1==1,c("Freq")],na.rm = T)+ 
                     sum(ref_test_tb[ref_test_tb$Var2==1,c("Freq")],na.rm = T))

r_date2_flood <- mask(r_date2_flood,r_ref)

ref_test_tb <- crosstab(r_date2_flood,r_ref_test)

## Compute Jaccard Index:

ref_test_tb$Freq[5]/(sum(ref_test_tb[ref_test_tb$Var1==1,c("Freq")],na.rm = T)+ 
                       sum(ref_test_tb[ref_test_tb$Var2==1,c("Freq")],na.rm = T))

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

#### Generate a plot for PCA with loadings and compare to Tassel Cap

var_labels <- rownames(loadings_df)

plot(loadings_df[,1],loadings_df[,2],
     type="p",
     pch = 20,
     col ="blue",
     xlab=names(loadings_df)[1],
     ylab=names(loadings_df)[2],
     ylim=c(-1,1),
     xlim=c(-1,1),
     axes = FALSE,
     cex.lab = 1.2)
axis(1, at=seq(-1,1,0.2),cex=1.2)
axis(2, las=1,at=seq(-1,1,0.2),cex=1.2) # "1' for side=below, the axis is drawned  on the right at location 0 and 1

box()    #This draws a box...

title(paste0("Loadings for component ", names(loadings_df)[1]," and " ,names(loadings_df)[2] ))
draw.circle(0,0,c(1.0,1),nv=200)#,border="purple",
text(loadings_df[,1],loadings_df[,2],var_labels,pos=1,cex=1)            
grid(2,2)


################### End of Script #########################

#plot(stack(r_pc1,r_pc2))
#layerStats(r_pc1,r_NDVI_mean )
#cor_pc <- layerStats(stack(r_pc1,r_NDVI_mean),'pearson', na.rm=T)
#cor_pc #PC1 correspond to the average mean by pixel as expected.
#plot(r_pc2)

### Potential follow up questions:
#1) compare wetness index from TCP in predicted flooded areas
#2) Use ROC to compare?
#3) Use NLCD and extract other reflectance curve.
#4) Use NLCD and recombine values using the general legend
#5) Use disaggregated values from NLCD %cover from 30m and correlate with the new indices...
#6) Compare area before and after classified as water with thresholding (don't forget to mask)
#