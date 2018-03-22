####################################    Flood Mapping Analyses   #######################################
############################  Analyze and map flooding from RITA hurricane  #######################################
#This script performs analyses for the Exercise 6 of the Short Course using reflectance data derived from MODIS.
#The goal is to map flooding from RITA using various reflectance bands from Remote Sensing platforms.
#Additional data is provided including FEMA flood region. 
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/13/2018 
#DATE MODIFIED: 03/21/2018
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
library(nnet)

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
script_path <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/R_scripts"
#source(file.path(script_path,function_preprocessing_and_analyses)) #source all functions used in this script 1.
#source(file.path(script_path,function_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#in_dir_reflectance <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/Exercise_6/data/reflectance_RITA"
in_dir_var <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_6/data/"
out_dir <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop//Exercise_6/outputs"
infile_reg_outline_RITA <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_6/data/revised_area_Rita/new_strata_rita_10282017.shp"

#region coordinate reference system
#http://spatialreference.org/ref/epsg/nad83-texas-state-mapping-system/proj4/
CRS_reg <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 
file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"exercise6_03212018" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
date_event <- ""
#ARG4
method_proj_val <- "bilinear" # "ngb"

#ARG9
#local raster name defining resolution, extent
ref_rast_name <- "/nfs/bparmentier-data/Data/Space_beats_time/Data/data_RITA_reflectance/revised_area_Rita/r_ref_Houston_RITA.tif"
infile_RITA_reflectance_date1 <- "mosaiced_MOD09A1_A2005265__006_reflectance_masked_RITA_reg_1km.tif"
infile_RITA_reflectance_date2 <- "mosaiced_MOD09A1_A2005273__006_reflectance_masked_RITA_reg_1km.tif"

nlcd_2006_filename <- "nlcd_2006_RITA.tif"
#infile_land_cover_date2 <- "nlcd_2011_RITA.tif"

infilename_class1 <- "class1.shp" #
infilename_class2 <- "class2.shp" #
infilename_class3 <- "class3.shp" #

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

###
## Second list.files and create raster images stack

nlcd2006_reg_RITA <- raster(file.path(in_dir_var,nlcd_2006_filename)) 
#11: open water
#90:Woody Wetlands
#95:Emergent Herbaceuous Wetlands

plot(nlcd2006_reg_RITA==90)
plot(nlcd2006_reg_RITA==95)
plot(nlcd2006_reg_RITA==11)

lc_legend_df <- read.table(file.path(in_dir_var,"nlcd_legend.txt"),sep=",")
lc_legend_df
View(lc_legend_df)

###### Read in modis 09

r_date1 <- brick(file.path(in_dir_var,infile_RITA_reflectance_date1))
r_date2 <- brick(file.path(in_dir_var,infile_RITA_reflectance_date2))

#SWIR1 (1230–1250 nm), SWIR2 (1628–1652 nm) and SWIR3 (2105–2155 nm).
band_refl_order <- c(3,4,1,2,5,6,7)

names(r_date1) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")
names(r_date2) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")

r_date2_MNDWI <- (subset(r_date2,"Green") - subset(r_date2,"SWIR2")) / (subset(r_date2,"Green") + subset(r_date2,"SWIR2"))
plot(r_date2_MNDWI)
r_date1_MNDWI <- (subset(r_date1,"Green") - subset(r_date1,"SWIR2")) / (subset(r_date1,"Green") + subset(r_date1,"SWIR2"))
plot(r_date1_MNDWI)

r_date2_NDVI <- (subset(r_date2,"NIR") - subset(r_date2,"Red")) / (subset(r_date2,"NIR") + subset(r_date2,"Red"))
plot(r_date2_NDVI)
r_date1_NDVI <- (subset(r_date1,"NIR") - subset(r_date1,"Red")) / (subset(r_date1,"NIR") + subset(r_date1,"Red"))
plot(r_date1_NDVI)
writeRaster(r_date1_NDVI,"ndvi_date1.rst")
NAvalue(r_date1_NDVI) <- 9999
writeRaster(r_date2_NDVI,"ndvi_date2.rst")
NAvalue(r_date2_NDVI) <- 9999
plot(r_date2_MNDWI)

#training_data_sf <- st_read("training1.shp")
class1_data_sf <- st_read("class1.shp")
class2_data_sf <- st_read("class2.shp")
class3_data_sf <- st_read("class3.shp")

list_class_sf <- list(class1_data_sf,class2_data_sf,class3_data_sf)
list_class_sp <- lapply(list_class_sf,function(x){as(x,"Spatial")})

r_x <- init(r_date2,"x") #raster with coordinates x
r_y <- init(r_date2,"x") #raster with coordiates y

r_stack <- stack(r_x,r_y,r_date2,r_date2_NDVI,r_date2_MNDWI)
names(r_stack) <- c("x","y","Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3","NDVI","MNDWI")

#lapply()
#list_pixels_df_test1 <- lapply(list_class_sp,function(x){extract(r_stack,x,df=T)})
#list_pixels_df_test2 <- lapply(list_class_sp,function(x){extract(r_stack,x,sp=T,df=T)})
#list_pixels_df <- lapply(1:length(list_class_sp),function(i){pix_df <- extract(r_stack,list_class_sp[[i]],df=T);pix_df$class_ID <- i})
list_pixels_df <- lapply(list_class_sp,function(x){extract(r_stack,x,df=T)})
list_class_pixels_df <- lapply(1:length(list_pixels_df),function(i){list_pixels_df[[i]]$class_ID <- i})

list_class_pixels_df <- list_pixels_df
list_class_pixels_df[[1]]$class_ID <- 1
list_class_pixels_df[[2]]$class_ID <- 2
list_class_pixels_df[[3]]$class_ID <- 3

list_tmp <- list_pixels_df
list_tmp[[1]]$class_id <- 1
View(list_tmp[[1]])
View(list_pixels_df[[1]])
View(list_class_pixels_df[[1]])

#list_pixels_df <- list_pixels_df_test2
#class1_data_sp$class_ID <- 1
#class2_data_sp$class_ID <- 2
#class3_data_sp$class_ID <- 3

## Note that both step above can be combined but for ease of understanding we kept them separate.
#pix_df <- extract(r_stack,training_data_sf,df=T)
pixels_df <- do.call("rbind",list_pixels_df)
pixels_df <- do.call("rbind",list_class_pixels_df)

View(pixels_df)
pix_df <- as.data.frame(pix_df)
names(pix_df)
names(pix_df) <- c("poly_ID","Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3","NDVI","MNDWI")
dim(pix_df)

#Show average MNDWI and NDVI for each class

### Need to add other extract for other land cover!!

#1) vegetation abd other
#2) Flooded vegetation
#3) Flooded area, or water (lake etc)

#so the classification will have three classes!!!


## Maybe also do a unsupervised?

#############

#Water
x_range <- range(pixels_df$Green,na.rm=T)
y_range <- range(pixels_df$NIR,na.rm=T)

#plot(Green~NIR,xlim=x_range,ylim=y_range,col="blue",subset(pixels_df,class_ID==1))
plot(Green~NIR,xlim=c(0,1),ylim=c(0,1),col="blue",subset(pixels_df,class_ID==1))
points(Green~NIR,col="green",subset(pixels_df,class_ID==2))
points(Green~NIR,col="red",subset(pixels_df,class_ID==3))

plot(Green~NIR,xlim=c(0,0.5),ylim=c(0,0.2),cex=0.2,col="blue",subset(pixels_df,class_ID==1))
points(Green~NIR,col="green",cex=0.2,subset(pixels_df,class_ID==2))
points(Green~NIR,col="red",cex=0.2,subset(pixels_df,class_ID==3))

histogram(r_date2)
histogram(pix_df)

# NRT MODIS
# Other

############### Using Classification and Regression Tree model (CART) #########

library(rpart)
# grow tree 
fit <- rpart(class_ID ~ Red +NIR + Blue + Green + SWIR1 + SWIR2 + SWIR3,
             method="class", 
             data=pixels_df)
plot(fit)
text(fit,cex=0.8)

# Much cleaner way is to plot the trained classification tree
plot(model.class, uniform=TRUE, main="Classification Tree")
text(model.class, cex=.8)
plot(fit, uniform=TRUE, 
     main="Classification Tree for Kyphosis")
text(fit, use.n=TRUE, all=TRUE, cex=.8)
#https://www.statmethods.net/advstats/cart.html

# prune the tree 
pfit<- prune(fit, cp=   fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])

# plot the pruned tree 
plot(pfit, uniform=TRUE, 
     main="Pruned Classification Tree for Kyphosis")
text(pfit, use.n=TRUE, all=TRUE, cex=.8)

r_predicted <- predict(r_stack,fit)
plot(r_predicted)
test <- unique(r_predicted)

test

# Now predict the subset data based on the model; prediction for entire area takes longer time
r_predicted <- predict(r_stack, fit, type='class', progress = 'text')

plot(r_predicted)


############### Using Neural Network #########

##### plot feature space:
pixels_df <- na.omit(pixels_df)
dim(pixels_df)
selected_var <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")#,"NDVI","MNDWI")
nrow(pixels_df)

test <- nnet(class_ID ~ Red +NIR + Blue + Green + SWIR1 + SWIR2 + SWIR3,
             x=subset(pixels_df,select=selected_var),
             y=subset(pixels_df,select=c("class_ID")),weights = rep( 1,1520),
             size=7)

test <- nnet(class_ID ~ Red +NIR + Blue + Green + SWIR1 + SWIR2 + SWIR3,
             x=subset(pixels_df,select=selected_var),
             y=subset(pixels_df,select=c("class_ID")))

?nnet

             #right_side_formula <- paste(explanatory_variables,collapse = " + ")
#model_formula_str <- paste0(y_var," ~ ",right_side_formula)
             
names(pixels_df)
neuralnet
pixels_df
##let's keep 30% of data for testing for each class


### Now do a unsupervised
## do a supervised
## Split training and testing... 

## Do ROC here!!!
## Do neural net, cart, random forest,svm

#nnet()

############### Using KNN or SVM #########

library(e1071)

svm_model <- svm(Species~ ., data=trainset, method="C-classification", kernel="linear")

svm_model <- svm(Species~ ., data=trainset, method="C-classification", kernel="linear")

mod_svm <- svm(class_ID ~ Red +NIR + Blue + Green + SWIR1 + SWIR2 + SWIR3,
             data=pixels_df,
             method="C-classification",
             kernel="linear")

mod_svm
summary(mod_svm)

plot(mod_svm)

# Now predict the subset data based on the model; prediction for entire area takes longer time
r_predicted_svm <- predict(r_stack, mod_svm)

plot(r_predicted_svm)

############# Compare methods with ROC #######


################### End of Script #########################

