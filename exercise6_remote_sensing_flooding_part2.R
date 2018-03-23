####################################    Flood Mapping Analyses   #######################################
############################  Analyze and map flooding from RITA hurricane  #######################################
#This script performs analyses for the Exercise 6 of the Short Course using reflectance data derived from MODIS.
#The goal is to map flooding from RITA using various reflectance bands from Remote Sensing platforms.
#Additional data is provided including FEMA flood region. 
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/13/2018 
#DATE MODIFIED: 03/23/2018
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
#library(gstat) #spatial interpolation and kriging methods
library(readxl) #functionalities to read in excel type data
library(psych) #pca/eigenvector decomposition functionalities
library(sf)
library(plotrix) #various graphic functions e.g. draw.circle
library(nnet)
library(rpart)
library(e1071)
library(caret)

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
out_suffix <-"exercise6_03232018" #output suffix for the files and ouptu folder #PARAM 8
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


###### Read in modis 09

r_date1 <- brick(file.path(in_dir_var,infile_RITA_reflectance_date1))
r_date2 <- brick(file.path(in_dir_var,infile_RITA_reflectance_date2))

#SWIR1 (1230–1250 nm), SWIR2 (1628–1652 nm) and SWIR3 (2105–2155 nm).
band_refl_order <- c(3,4,1,2,5,6,7)

names(r_date1) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")
names(r_date2) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")

r_date2_MNDWI <- (subset(r_date2,"Green") - subset(r_date2,"SWIR2")) / (subset(r_date2,"Green") + subset(r_date2,"SWIR2"))
plot(r_date2_MNDWI,zlim=c(-1,1))
r_date1_MNDWI <- (subset(r_date1,"Green") - subset(r_date1,"SWIR2")) / (subset(r_date1,"Green") + subset(r_date1,"SWIR2"))
plot(r_date1_MNDWI,zlim=c(-1,1))

r_date2_NDVI <- (subset(r_date2,"NIR") - subset(r_date2,"Red")) / (subset(r_date2,"NIR") + subset(r_date2,"Red"))
plot(r_date2_NDVI)
r_date1_NDVI <- (subset(r_date1,"NIR") - subset(r_date1,"Red")) / (subset(r_date1,"NIR") + subset(r_date1,"Red"))
plot(r_date1_NDVI)
#writeRaster(r_date1_NDVI,"ndvi_date1.rst")
#NAvalue(r_date1_NDVI) <- 9999
#writeRaster(r_date2_NDVI,"ndvi_date2.rst")
#NAvalue(r_date2_NDVI) <- 9999
#plot(r_date2_MNDWI)

class2_sites.shp
#training_data_sf <- st_read("training1.shp")
class1_data_sf <- st_read(file.path(in_dir_var,"class1_sites.shp"))
class2_data_sf <- st_read(file.path(in_dir_var,"class2_sites.shp"))
class3_data_sf <- st_read(file.path(in_dir_var,"class3_sites.shp"))

class_data_sf <- rbind(class1_data_sf,class2_data_sf,class3_data_sf)
class_data_sf$poly_ID <- 1:nrow(class_data_sf) #unique ID for each polygon
nrow(class_data_sf)

class_data_sp <- as(class_data_sf,"Spatial")

###merg sf data
list_class_sf <- list(class1_data_sf,class2_data_sf,class3_data_sf)
list_class_sp <- lapply(list_class_sf,function(x){as(x,"Spatial")})

r_x <- init(r_date2,"x") #raster with coordinates x
r_y <- init(r_date2,"x") #raster with coordiates y

r_stack <- stack(r_x,r_y,r_date2,r_date2_NDVI,r_date2_MNDWI)
names(r_stack) <- c("x","y","Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3","NDVI","MNDWI")
pixels_df <- extract(r_stack,class_data_sp,df=T)

dim(pixels_df) #We have 1547 pixels extracted
class_data_df <- class_data_sf
st_geometry(class_data_df) <- NULL #this will coerce the sf object into a data.frame
pixels_df <- merge(pixels_extracted_df,class_data_df,by.x="ID",by.y="poly_ID")

#Show average MNDWI and NDVI for each class
### Need to add other extract for other land cover!!
#View(pixels_df)
#names(pixels_df)
#class(pixels_df)

#1) vegetation abd other
#2) Flooded vegetation
#3) Flooded area, or water (lake etc)

#so the classification will have three classes!!!

######## Examining sites data used for the classification

#Water
x_range <- range(pixels_df$Green,na.rm=T)
y_range <- range(pixels_df$NIR,na.rm=T)

###Add legend?
plot(NIR~Green,xlim=x_range,ylim=y_range,cex=0.2,col="blue",subset(pixels_df,class_ID==1))
points(NIR~Green,col="green",cex=0.2,subset(pixels_df,class_ID==2))
points(NIR~Green,col="red",cex=0.2,subset(pixels_df,class_ID==3))

plot(NDVI~MNDWI,xlim=c(-1,1),ylim=c(-1,1),cex=0.2,col="blue",subset(pixels_df,class_ID==1))
points(NDVI~MNDWI,col="green",cex=0.2,subset(pixels_df,class_ID==2))
points(NDVI~MNDWI,col="red",cex=0.2,subset(pixels_df,class_ID==3))

histogram(r_date2)

boxplot(MNDWI~class_ID,pixels_df,main="Boxplot for MNDWI per class")

############### Split training and testing #############

##let's keep 30% of data for testing for each class
## Split training and testing... 

#pixels_df
pixels_df$pix_ID <- 1:nrow(pixels_df)
prop <- 0.3
table(pixels_df$class_ID)
set.seed(100)
### This is for one class:
#Better as a function but we use a loop for clarity here:

list_data_df <- vector("list",length=3)
for(i in 1:3){
  data_df <- subset(pixels_df,class_ID==i)
  data_df$pix_id <- 1:nrow(data_df)
  indices <- as.vector(createDataPartition(data_df$pix_ID,p=0.7,list=F))
  data_df$training <-  as.numeric(data_df$pix_id %in% indices)
  #table(data_df$training)
  list_data_df[[i]] <- data_df
}

data_df <- do.call(rbind,list_data_df)

dim(data_df)
View(data_df)

data_training <- subset(data_df,training==1)

## Do neural net, cart, random forest,svm
############### Using Classification and Regression Tree model (CART) #########

# grow tree 
#mod_rpart <- rpart(class_ID ~ Red +NIR + Blue + Green + SWIR1 + SWIR2 + SWIR3,
#             method="class", 
#             data=pixels_df)
mod_rpart <- rpart(class_ID ~ Red +NIR + Blue + Green + SWIR1 + SWIR2 + SWIR3,
                   method="class", 
                   data=data_training)

# Plot the fitted  classification tree
plot(mod_rpart, uniform=TRUE, main="Classification Tree")
text(mod_rpart, cex=.8)

plot(mod_rpart, uniform=TRUE, main="Classification Tree")
text(mod_rpart, use.n=TRUE, all=TRUE, cex=.8)

# Now predict the subset data based on the model; prediction for entire area takes longer time
r_predicted_rpart <- predict(r_stack,mod_rpart, type='class', progress = 'text')

plot(r_predicted_rpart)
r_predicted_rpart <- ratify(r_predicted_rpart)
rat <- levels(r_predicted_rpart)[[1]]
rat$legend <- c("vegetation","wetland","water")
levels(r_predicted_rpart) <- rat
levelplot(r_predicted_rpart, maxpixels = 1e6,
          col.regions = c("green","blue","darkblue"),
          scales=list(draw=FALSE),
          main = "Classification Tree")

############### Using KNN or SVM #########

## set class_ID as factor to generate classification
pixels_df$class_ID <- as.factor(pixels_df$class_ID)
mod_svm <- svm(class_ID ~ Red +NIR + Blue + Green + SWIR1 + SWIR2 + SWIR3,
               data=data_training,
               method="C-classification",
               kernel="linear") # can be radial

summary(mod_svm)

#plot(mod_svm)

# Now predict the subset data based on the model; prediction for entire area takes longer time
raster_outfilename <- paste0("r_predicted_svm_",out_suffix,file_format)
r_predicted_svm <- predict(r_stack, mod_svm,
                           progress = 'text',
                           filename=raster_outfilename,
                           overwrite=T)

plot(r_predicted_svm)
histogram(r_predicted_svm)
r_predicted_svm <- ratify(r_predicted_svm)
rat <- levels(r_predicted_svm)[[1]]
rat$legend <- c("vegetation","wetland","water")
levels(r_predicted_svm) <- rat
levelplot(r_predicted_svm, maxpixels = 1e6,
          col.regions = c("green","blue","darkblue"),
          scales=list(draw=FALSE),
          main = "SVM classification")

### get confusion matrix?
#table(pred,y)
#https://rischanlab.github.io/SVM.html

################# Using Neural Network ##################

##### plot feature space:
data_training_df <- na.omit(data_training)
dim(data_training_df)
selected_var <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")#,"NDVI","MNDWI")
nrow(data_training_df)

data_training_df_subset <- subset(data_training,select=c("class_ID",selected_var))
dim(data_training_df_subset)
data_training_df_subset$class_ID

mod_nnet <- nnet(as.factor(class_ID) ~ Red +NIR + Blue + Green + SWIR1 + SWIR2 + SWIR3,
             #x=subset(pixels_df,select=selected_var),
             #y=subset(pixels_df,select=c("class_ID")),
             data=data_training_df_subset,
             #weights = 1,
             size=10,
             type="class")

?nnet

NAvalue(r_stack) <- -9999
raster_outfilename <- paste0("r_predicted_nnet_",out_suffix,file_format)
r_predicted_nnet <- predict(subset(r_stack,3:9),mod_nnet,type="raw",overwrite=T)
r_predicted_nnet <- predict(r_stack,mod_nnet,
                            type="raw",
                            filename=raster_outfilename,
                            overwrite=T)
??raster::predict
r_predicted_nnet <- predict(subset(r_stack,3:9),mod_nnet,type="class",na.rm=T,overwrite=T)
names(r_stack)
pnnet <- raster::predict(object=xn, model=results, fun=predict.NNET, na.rm=TRUE, factors=factors,
                         filename=fullname, progress='text', overwrite=TRUE)
https://rdrr.io/cran/BiodiversityR/src/R/ensemble.raster.R
mod_nnet$fitted.values
?predict.nnet
#plot(r_predicted_nnet)
#unique(r_predicted_nnet)
#model <- mlp(dat$inputsTrain, dat$targetsTrain, size=1, learnFunc="SCG", learnFuncParams=c(0, 0, 0, 0), 
#             maxit=400, inputsTest=dat$inputsTest, targetsTest=dat$targetsTest)

plot(r_predicted_nnet)
histogram(r_predicted_nnet)

#nnet()

######## Compare methods for the performance #################

dim(data_df)
data_testing <- subset(data_df,training==0)
#predict()
#Not working here:
#testing_rpart <- predict.rpart(mod_rpart, data_testing,type='class')
testing_rpart <- predict(mod_rpart, data_testing,type='class')

#class(data_testing)
#testing_svm <- e1071::predict(data_testing,mod_svm, type='class')
testing_svm <- predict(mod_svm,data_testing, type='class')

#### Generate confusion matrix!!!

tb_rpart <- table(testing_rpart,data_test$class_ID)
tb_svm <- table(testing_svm,data_test$class_ID)

sum(table(testing_rpart))

tb_rpart[1]/sum(tb_rpart[1,])
tb_svm[1]/sum(tb_svm[,1])
#overall accuracy for svm
sum(diag(tb_svm))/sum(table(testing_svm))
#overall accuracy for rpart
sum(diag(tb_rpart))/sum(table(testing_rpart))


#Generate more accuracy measurements from CARET
accuracy_info_svm <- confusionMatrix(testing_svm,data_test$class_ID, positive = NULL)
accuracy_info_rpart <- confusionMatrix(testing_rpart,data_test$class_ID, positive = NULL)

accuracy_info_rpart$overall
accuracy_info_svm$overall

#### write out the results:


################### End of Script #########################

