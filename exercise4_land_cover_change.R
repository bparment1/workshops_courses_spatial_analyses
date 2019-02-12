####################################   Land Use and Land Cover Change   #######################################
############################  Analyze Land Cover change in Houston  #######################################
#This script performs analyses for a land cover change model using NLCD aggregated values.
#The goal is to assess land cover change using two land cover maps in the Houston areas.
#Additional datasets are provided for the land cover change modeling. A model is built for Harris county.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/16/2018 
#DATE MODIFIED: 02/11/2019
#Version: 1
#PROJECT: AAG 2019 Geospatial Analysis workshop, Geospatial Data Analysis Short Course, Geo-Computation and Environmental Analysis Yale.  
#TO DO:
#
#COMMIT: adding sampling
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
library(sf) #spatial objects and functionalities
library(plotrix) #various graphic functions e.g. draw.circle
library(TOC) # TOC and ROC for raster images
library(ROCR) # ROCR general for data.frame
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

#####  Parameters and argument set up ###########

#Separate inputs and outputs directories
#in_dir_var <- "/home/user/ost4sem/exercise/Exercise_4/data/"
#out_dir <- "/home/user/ost4sem/exercise/Exercise_4/outputs"

in_dir_var <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2019_geospatial_workshop/Exercise_4/data"
out_dir <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2019_geospatial_workshop/Exercise_4/outputs"

### General parameters

#NLCD coordinate reference system: we will use this projection rather than TX.
CRS_reg <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
file_format <- ".tif" #raster output format 
NA_flag_val <- -9999 # NA value assigned to output raster
out_suffix <-"exercise4_02072019" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE # if TRUE, a output dir using output suffix will be created
method_proj_val <- "bilinear" # method option for the reprojection and resampling 
gdal_installed <- TRUE #if TRUE, GDAL is used to generate distance files

### Input data files
rastername_county_harris <- "harris_county_mask.tif" #Region of interest: extent of Harris County
elevation_fname <- "srtm_Houston_area_90m.tif" #SRTM elevation
roads_fname <- "r_roads_Harris.tif" #Road count for Harris county

### Aggreagate NLCD input files
infile_land_cover_date1 <- "agg_3_r_nlcd2001_Houston.tif"
infile_land_cover_date2 <- "agg_3_r_nlcd2006_Houston.tif"
infile_land_cover_date3 <- "agg_3_r_nlcd2011_Houston.tif"

infile_name_nlcd_legend <- "nlcd_legend.txt"
infile_name_nlcd_classification_system <- "classification_system_nlcd_legend.xlsx"

######################### START SCRIPT ###############################

## First create an output directory to separate inputs and outputs

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

###########################################
### PART I: READ AND VISUALIZE DATA #######

r_lc_date1 <- raster(file.path(in_dir_var,infile_land_cover_date1)) #NLCD 2001 
r_lc_date2 <- raster(file.path(in_dir_var,infile_land_cover_date2)) #NLCD 2006
r_lc_date3 <- raster(file.path(in_dir_var,infile_land_cover_date2)) #NLCD 2011

lc_legend_df <- read.table(file.path(in_dir_var,infile_name_nlcd_legend),
                           stringsAsFactors = F,
                           sep=",")

head(lc_legend_df) # Inspect data

plot(r_lc_date2) # View NLCD 2006, we will need to add the legend use the appropriate palette!!

### Let's add legend and examine existing land cover categories

freq_tb_date2 <- freq(r_lc_date2)
head(freq_tb_date2) #view first 5 rows, note this is a matrix object.

### Let's generate a palette from the NLCD legend information to view the existing land cover for 2006.
names(lc_legend_df)
dim(lc_legend_df) #contains a lot of empty rows

lc_legend_df<- subset(lc_legend_df,COUNT>0) #subset the data to remove unsured rows
### Generate a palette color from the input Red, Green and Blue information using RGB encoding:


lc_legend_df$rgb <- paste(lc_legend_df$Red,lc_legend_df$Green,lc_legend_df$Blue,sep=",") #combine

### row 2 correspond to the "open water" category
color_val_water <- rgb(lc_legend_df$Red[2],lc_legend_df$Green[2],lc_legend_df$Blue[2],maxColorValue = 255)
color_val_developed_high <- rgb(lc_legend_df$Red[7],lc_legend_df$Green[7],lc_legend_df$Blue[7],maxColorValue = 255)

lc_col_palette <- c(color_val_water,color_val_developed_high)

barplot(c(1,1), 
        col=lc_col_palette,
        main="Visualization of color palette for NLCD land cover",
        names.arg=c("Open water",	"Developed, High Intensity"),las=1)

### Let's generate a color for all the land cover categories by using lapply and function
n_cat <- nrow(lc_legend_df)
lc_col_palette <- lapply(1:n_cat,
                 FUN=function(i){rgb(lc_legend_df$Red[i],lc_legend_df$Green[i],lc_legend_df$Blue[i],maxColorValue = 255)})
lc_col_palette <- unlist(lc_col_palette)

lc_legend_df$palette <- lc_col_palette

r_lc_date2 <- ratify(r_lc_date2) # create a raster layer with categorical information
rat <- levels(r_lc_date2)[[1]] #This is a data.frame with the categories present in the raster

lc_legend_df_date2 <- subset(lc_legend_df,lc_legend_df$ID%in% (rat[,1])) #find the land cover types present in date 2 (2006)
rat$legend <- lc_legend_df_date2$NLCD.2006.Land.Cover.Class #assign it back in case it is missing
levels(r_lc_date2) <- rat #add the information to the raster layer

### Now generate a plot of land cover with the NLCD legend and palette
levelplot(r_lc_date2, 
          col.regions = lc_legend_df_date2$palette,
          scales=list(draw=FALSE),
          main = "NLCD 2006")

################################################
###  PART II : Analyze change and transitions

## As the plot shows for 2006, we have 15 land cover types. Analyzing such complex categories in terms of decreasse (loss), increase (gain), 
# persistence in land cover will generate a large number of transitions (potential up to 15*15=225 transitions in this case!)

## To generalize the information, let's aggregate leveraging the hierachical nature of NLCD Anderson Classification system.

lc_system_nlcd_df <- read_xlsx(file.path(in_dir_var,infile_name_nlcd_classification_system))
head(lc_system_nlcd_df) #inspect data

### Let's identify existing cover and compute change:
r_stack_nlcd <- stack(r_lc_date1,r_lc_date2)
freq_tb_nlcd <- as.data.frame(freq(r_stack_nlcd,merge=T))
head(freq_tb_nlcd)

dim(lc_system_nlcd_df) # We have categories that are not relevant to the study area and time period.
lc_system_nlcd_df <- subset(lc_system_nlcd_df,id_l2%in%freq_tb_nlcd$value ) 
dim(lc_system_nlcd_df) # Now 15 land categories instead of 20.

### Selectet relevant columns for the reclassification
rec_df <- lc_system_nlcd_df[,c(2,1)]
r_date1_rec <- subs(r_lc_date1,rec_df,by="id_l2","id_l1")
r_date2_rec <- subs(r_lc_date2,rec_df,by="id_l2","id_l1")

plot(r_date1_rec)

rec_xtab_df <- crosstab(r_date1_rec,r_date2_rec,long=T)
names(rec_xtab_df) <- c("2001","2006","freq")

head(rec_xtab_df)
dim(rec_xtab_df) #9*9 possible transitions if we include NA values
print(rec_xtab_df) # View the full table

which.max(rec_xtab_df$freq)
rec_xtab_df[11,] # Note the most important transition is persistence!!

### Let's rank the transition:
class(rec_xtab_df)
is.na(rec_xtab_df$freq)
rec_xtab_df_ranked <- rec_xtab_df[order(rec_xtab_df$freq,decreasing=T) , ]
head(rec_xtab_df_ranked) # Unsurprsingly, top transitions are persistence categories

### Let's examine the overall change in categories rather than transitions

label_legend_df <- data.frame(ID=lc_system_nlcd_df$id_l1,name=lc_system_nlcd_df$name_l1)
r_stack <- stack(r_date1_rec,r_date2_rec)

lc_df <- freq(r_stack,merge=T)
names(lc_df) <- c("value","date1","date2")
lc_df$diff <- lc_df$date2 - lc_df$date1 #difference for each land cover categories over the 2001-2011 time period
head(lc_df) # Quickly examine the output

### Add relevant categories
lc_df <- merge(lc_df,label_legend_df,by.x="value",by.y="ID",all.y=F)
lc_df <- lc_df[!duplicated(lc_df),] #remove duplictates
head(lc_df) # Note the overall cahnge

#### Now visualize the overall land cover changes
barplot(lc_df$diff,names.arg=lc_df$name,las=2)
total_val  <- sum(lc_df$date1)
lc_df$perc_change <- 100*lc_df$diff/total_val 
barplot(lc_df$perc_change,names.arg=lc_df$name,las=2)

### Create a change image to map all pixels that transitioned to the developed category:  

r_cat2 <- r_date2_rec==2 # developped on date 2
r_not_cat2 <- r_date1_rec!=2 #remove areas that were already developed in date1, we do not want persistence

r_change <- r_cat2 * r_not_cat2 #mask
plot(r_change,main="Land transitions to developed over 2001-2006")
change_tb <- freq(r_change) #Find out how many pixels transitions to developped

#####################################
############# PART III: Process and Prepare variables for land change modeling ##############

## y= 1 if change to urban over 2001-2011
### Explanatory variables:
#var1: distance to existing urban in 2001
#var2: distance to road in 2001
#var3: elevation, low slope better for new development
#var4: past land cover state that may influence future land change

## 1) Generate var1 and var2 : distance to developped and distance to roads

### Distance to existing in 2001: prepare information
r_cat2<- r_date1_rec==2 #developped in 2001
plot(r_cat2)
cat_bool_fname <- "developped_2001.tif" #input for the distance to road computation
writeRaster(r_cat2,filename = cat_bool_fname,overwrite=T)

### Read in data for road count
r_roads <- raster(file.path(in_dir_var,roads_fname))
plot(r_roads,colNA="black")
res(r_roads)

#### Aggregate to match the NLCD data resolution
r_roads_90m <- aggregate(r_roads,
                         fact=3, #factor of aggregation in x and y
                         fun=mean) #function used in aggregation values
plot(r_roads_90m)

r_roads_bool <- r_roads_90m > 0
plot(r_roads_bool)

roads_bool_fname <- "roads_bool.tif" #input for the distance to road computation
writeRaster(r_roads_bool,filename = roads_bool_fname,overwrite=T)

### This part could be transformed into a function but we keep it for clarity and learning:
if(gdal_installed==TRUE){
  
  ## Distance from developped land in 2001
  srcfile <- cat_bool_fname  
  dstfile_developped <- file.path(out_dir,paste("developped_distance_",out_suffix,file_format,sep=""))
  n_values <- "1"
  
  ### Prepare GDAL command: note that gdal_proximity doesn't like when path is too long
  cmd_developped_str <- paste("gdal_proximity.py",
                              basename(srcfile),
                              basename(dstfile_developped),
                              "-values",n_values,sep=" ")

  ### Distance from roads
  srcfile <- roads_bool_fname 
  dstfile_roads <- file.path(out_dir,paste("roads_bool_distance_",out_suffix,file_format,sep=""))
  n_values <- "1"
  
  ### Prepare GDAL command: note that gdal_proximity doesn't like when path is too long
  cmd_roads_str <- paste("gdal_proximity.py",basename(srcfile),
                         basename(dstfile_roads),
                         "-values",n_values,sep=" ")

  sys_os <- as.list(Sys.info())$sysname #Find what OS system is in use.
  
  if(sys_os=="Windows"){
    shell(cmd_developped_str)
    shell(cmd_roads_str)
  }else{
    system(cmd_developped_str)
    system(cmd_roads_str)
  }
  r_roads_distance <- raster(dstfile_roads)
  r_developped_distance <- raster(dstfile_developped)
  
}else{
  r_developped_distance <- raster(file.path(in_dir,paste("developped_distance",file_format,sep="")))
  r_roads_distance <- raster(file.path(in_dir_var,paste("roads_bool_distance",file_format,sep="")))
}

plot(r_developped_distance) #This is at 90m.
plot(r_roads_distance) #This is at 90m.

#Now rescale the distance...
min_val <- cellStats(r_roads_distance,min) 
max_val <- cellStats(r_roads_distance,max)
r_roads_dist <-  (max_val - r_roads_distance) / (max_val - min_val) #high values close to 1 for areas close to roads


min_val <- cellStats(r_developped_distance,min) 
max_val <- cellStats(r_developped_distance,max)
r_developped_dist <-  (max_val - r_developped_distance) / (max_val - min_val)

## 2) Generate var3 : slope
r_elevation <- raster(file.path(in_dir_var,elevation_fname))

projection(r_elevation) # This is not in the same projection as the study area. 
r_elevation_reg <- projectRaster(from= r_elevation, #input raster to reproject
                                 to= r_date1_rec, #raster with desired extent, resolution and projection system
                                 method= method_proj_val) #method used in the reprojection

r_slope <- terrain(r_elevation_reg,unit="degrees")

## 3) Generate var4 : past land cover state

### reclass Land cover
r_mask <- r_date1_rec==2 #Remove developed land from  
r_date1_rec_masked <- mask(r_date1_rec,r_mask,maskvalue=1)
plot(r_date1_rec_masked)

#####################################
############# PART IV: Run Model and perform assessment ##############

##############
###### Step 1: Consistent masking and generate mask removing water (1) and developped (2) in 2001

#r_mask <- (r_date1_rec!=2)*(r_date1_rec!=1)*r_county_harris
r_mask <- (r_date1_rec!=2)*(r_date1_rec!=1)
plot(r_mask)
NAvalue(r_mask) <- 0 

### Read in focus area for the modeling:
r_county_harris <- raster(file.path(in_dir_var,rastername_county_harris))
### Screen for area of interest
r_mask <- r_mask * r_county_harris
r_mask[r_mask==0]<-NA
plot(r_mask)

### Check the number of NA and pixels in the study area:
tb_study_area <- freq(r_mask)
print(tb_study_area)

## Generate dataset for Harris county
r_variables <- stack(r_change,r_date1_rec_masked,r_slope,r_roads_dist,r_developped_dist)
r_variables <- mask(r_variables,mask=r_mask) # mask to keep relevant area
names(r_variables) <- c("change","land_cover","slope","roads_dist","developped_dist")
NAvalue(r_variables) <- NA_flag_val

## Examine all the variables
plot(r_variables)

### Check for consistency in mask:
NA_freq_tb <- freq(r_variables,value=NA,merge=T)
### Notice that the number of NA is not consistent.
### Let's recombine all NA for consstencies:
plot(r_variables)
r_NA <- r_variables > -1 #There are no negative values in the raster stack

r_valid_pixels <- overlay(r_NA,fun=sum)
plot(r_valid_pixels)
dim(r_NA)
freq(r_valid_pixels)

r_mask <- r_valid_pixels > 0
r_variables <- mask(r_variables,r_mask)
#r_variables <- freq(r_test2,value=NA,merge=T)
names(r_variables) <- c("change","land_cover","slope","roads_dist","developped_dist")

###############
###### Step 2: Fit glm model and generate predictions

variables_df <- na.omit(as.data.frame(r_variables)) # convert raster stack to data.frame
dim(variables_df)
variables_df$land_cover <- as.factor(variables_df$land_cover) # convert to categorical variable
variables_df$change <- as.factor(variables_df$change) # convert to categorical variable

names(variables_df)
#names(variables_df) <- c("change","land_cover","elevation","roads_dist","developped_dist")

## Now generate the glm model using the logistic specification:
mod_glm <- glm(change ~ land_cover + slope + roads_dist + developped_dist, 
           data=variables_df , family=binomial())

print(mod_glm)
summary(mod_glm)

### Genreate probability map from fitted model:
r_p <- predict(r_variables, mod_glm, type="response")
plot(r_p)
histogram(r_p)


### save outputs:
r_out <- stack(r_variables,r_p)
names(r_out)[nlayers(r_out)] <- "prob"

out_filename <- paste0("r_variables_harris_county","_",out_suffix,file_format)
writeRaster(r_out,
            filename=file.path(out_dir,out_filename),
            bylayer=T,
            suffix=names(r_out),
            overwrite=T)
dim(r_p)
ncell(r_p)
dim(variables_df)

###
r_x <-init(r_out,v="x")
r_y <-init(r_out,v="y")

r_out <- stack(r_out,r_x,r_y)

n_name<- nlayers((r_out))
names(r_out)[(n_name-1):n_name] <- c("x","y")

#test <- st_as_sf(r_out)
variables_out_df <- na.omit(as.data.frame(r_out)) # convert raster stack to data.frame
dim(variables_out_df)
names(variables_out_df)

out_filename <- paste0("r_variables_harris_county","_",out_suffix,".txt")

write.table(variables_out_df,
            file=file.path(out_dir,out_filename),
            sep=",",
            row.names=F)

###############
###### Step 3: Model assessment with ROC

## We use the TOC package since it allows for the use of raster layers.

r_change_harris <- subset(r_variables,"change")

## These are the inputs for the assessment
plot(r_change_harris) # boolean reference variable
plot(r_mask) # mask for relevant observation
plot(r_p) # index variable to assess

roc_rast <- ROC(index=r_p, 
                  boolean=r_change_harris, 
                  mask=r_mask,
                  nthres=100)

plot(roc_rast)
slot(roc_rast,"AUC") #this is the AUC from ROC for the logistic modeling

toc_rast <- TOC(index=r_p, 
                  boolean=r_change_harris, 
                  mask=r_mask,
                  nthres=100)

plot(toc_rast)
slot(toc_rast,"AUC") #this is the AUC from TOC for the logistic modeling

##############
##### Step 4: Prepare model for predictions: Split data into training and testing

seed_number <- 100
if (seed_number>0) {
  set.seed(seed_number)                        #Using a seed number allow results based on random number to be compared...
}

#variables_df
#nel<-length(variables_df)
#dates_list<-vector("list",nel) #list of one row data.frame
prop <- 0.3

n<- nrow(variables_df)

ns<-n-round(n*prop)   #Create a sample from the data frame with 70% of the rows
nv<-n-ns              #create a sample for validation with prop of the rows
ind.training <- sample(nrow(variables_df), size=ns, replace=FALSE) #This selects the index position for 70% of the rows taken randomly
ind.testing <- setdiff(1:nrow(variables_df), ind.training)

variables_df$training[ind.training] <-  1
variables_df$training[ind.testing] <-  0

sum(variables_df$training)/length(variables_df$training)

data_training_df <- variables_df[variables_df$training==1,]
data_testing_df <- variables_df[variables_df$training==1,]

## Now generate the glm model using the logistic specification:
mod_glm <- glm(change ~ land_cover + slope + roads_dist + developped_dist, 
               data=data_training_df , family=binomial())

print(mod_glm)
summary(mod_glm)

r_p <- predict(data_testing, mod_glm, type="response")

r_training <- rasterize(variables_df,y=r_variables,field="training")

### Generate probability map from fitted model:
r_p <- predict(r_variables, mod_glm, type="response")
plot(r_p)
histogram(r_p)



###############################  End of script  #####################################

