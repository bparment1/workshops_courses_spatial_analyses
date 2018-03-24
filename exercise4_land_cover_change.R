####################################   Land Use and Land Cover Change   #######################################
############################  Analyze Land Cover change in Houston  #######################################
#This script performs analyses for the Exercise 4 of the Short Course using reflectance data derived from MODIS.
#The goal is to assess land cover change using two land cover maps.
#Additional datasets are provided for the land cover change modeling. 
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/16/2018 
#DATE MODIFIED: 03/24/2018
#Version: 1
#PROJECT: SESYNC and AAG 2018 workshop/Short Course preparation
#TO DO:
#
#COMMIT: more processing of elevevation and roads
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

#function_analyses <- "exercise2_fire_alaska_analyses_functions_03232017.R" #PARAM 1
script_path <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/R_scripts"
#source(file.path(script_path,function_preprocessing_and_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

#in_dir_reflectance <- "/nfs/bparmentier-data/Data/workshop_spatial/GIS_training/Exercise_6/data/reflectance_RITA"
in_dir_var <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_4/data/"
out_dir <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_4/outputs"
infile_reg_outline <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_4/revised_area_Rita/new_strata_rita_10282017.shp"

#region coordinate reference system
#http://spatialreference.org/ref/epsg/nad83-texas-state-mapping-system/proj4/
CRS_reg <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 
file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"exercise4_03242018" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
date_event <- ""
#ARG4
method_proj_val <- "bilinear" # "ngb"

#ARG9
#local raster name defining resolution, extent
ref_rast_name <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_4/data/r_ref_Houston_RITA.tif"

elevation_fname <- "srtm_Houston_area_90m.tif"
roads_fname <- "r_roads_Harris.tif"

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

# agg_3_r_nlcd2006_Houston.tif
# agg_3_r_nlcd2011_Houston.tif
# nlcd_legend.txt
# r_2006_nlcd30m_Houston.tif.aux.xml
# r_2006_nlcd30m_Houston.tif
# r_2011_nlcd30m_Houston.tif.aux.xml
# r_2011_nlcd30m_Houston.tif
infile_land_cover_date1 <- file.path(in_dir_var,"r_2001_nlcd30m_Houston.tif")
#infile_land_cover_date1 <- "r_2006_nlcd30m_Houston.tif"
infile_land_cover_date2 <- file.path(in_dir_var,"r_2011_nlcd30m_Houston.tif")
 
r_lc_date1 <- raster(infile_land_cover_date1) 
r_lc_date2 <- raster(infile_land_cover_date2) 

lc_legend_df <- read.table(file.path(in_dir_var,"nlcd_legend.txt"),sep=",")
lc_legend_df
View(lc_legend_df)

### get confusion matrix?
#table(pred,y)
#https://rischanlab.github.io/SVM.html

#
plot(r_lc_date1,col=lc_col) #will need to add the legend!!
plot(r_lc_date2,col=lc_col)

names(lc_legend_df)
dim(lc_legend_df)

freq_tb_date1 <- freq(r_lc_date1)
freq_tb_date2 <- freq(r_lc_date2)

View(freq_tb_date1)
View(freq_tb_date2)

freq_tb_date1 <- merge(freq_tb_date1,lc_legend_df,by.x="value",by.y="ID")
freq_tb_date2 <- merge(freq_tb_date2,lc_legend_df,by.x="value",by.y="ID")

freq_tb_date1$date <- 2001
freq_tb_date2$date <- 2011

### No category disappeared:
freq_tb_date1$value==freq_tb_date2$value
#lc_df <- data.frame(ID=freq_tb_date1$value,
#                    lc2006=freq_tb_date1$COUNT,#count is wrong!!!! that is from 2006
#                    lc2011=freq_tb_date2$COUNT)

lc_df <- data.frame(ID=freq_tb_date1$value,
           lc2001=freq_tb_date1$count,
           lc2011=freq_tb_date2$count,
           name=freq_tb_date1$NLCD.2006.Land.Cover.Class)

lc_df$diff <- lc_df$lc2011 - lc_df$lc2001 

View(lc_df)

### Let's make plot of land cover types and differences
lc_legend_df<- subset(lc_legend_df,COUNT>0)

lc_legend_df$rgb <- paste(lc_legend_df$Red,lc_legend_df$Green,lc_legend_df$Blue,sep=",")
i<-1
n_cat <- nrow(lc_legend_df)
lc_col <- lapply(1:n_cat,function(i){rgb(lc_legend_df$Red[i],lc_legend_df$Green[i],lc_legend_df$Blue[i],maxColorValue = 255)})
lc_col <- unlist(lc_col)

r_lc_date1 <- ratify(r_lc_date1)
rat <- levels(r_lc_date1)[[1]]

as.character(lc_df$name)
#rat$legend <- c("vegetation","wetland","water")
rat$legend <- as.character(lc_df$name)
levels(r_lc_date1) <- rat
levelplot(r_lc_date1, maxpixels = 1e6,
          col.regions = lc_col,
          scales=list(draw=FALSE),
          main = "NLCD 2001")
plot(r_lc_date1)
#positive means increase, negative a decrease

xtab_df <- crosstab(r_lc_date1,r_lc_date2)
dim(xtab_df)
View(xtab_df)

### This removes the zero transitions:
xtab_df_long <- crosstab(r_lc_date1,r_lc_date2,long=T)
View(xtab_df_long)

### What is the largest land transition?

xtab_df_long[which.max(xtab_df_long$Freq),]
#23:Developed, Medium Intensity

lc_transitions_df <- xtab_df_long[order(xtab_df_long$Freq,decreasing=T),]
View(lc_transitions_df)
#let's remove all the persisence classes:

#?crosstab

persistence_cat <- lc_transitions_df$r_2001_nlcd30m_Houston==lc_transitions_df$r_2011_nlcd30m_Houston

sum(persistence_cat)

lc_change_df <- lc_transitions_df[!persistence_cat,]
dim(lc_change_df)
dim(lc_transitions_df)

#View(lc_change_df)
# Hay pasture: 81
#Cultivated Crops: 82

#21: Developed, Open Space
#22: Developed, Low Intensity
#23: Developed, Medium Intensity
#24: Developed, High Intensity

#highest transition: 21 to 23
#second hightest: 81-23

# Too much information: let's aggregate and summize the info:

unique(lc_legend_df$ID)

df_reclasss <- lc_legend_df$ID

lc_df$ID
View(lc_df)

#as.character(lc_df$ID)[1]

infile_name_nlcd_legend <- list.files(path=in_dir_var,pattern="*.xlsx",full.names=T)

nlcd_legend_df <- read_xlsx(infile_name_nlcd_legend)
View(nlcd_legend_df)
names(nlcd_legend_df)

class(lc_df$ID)
nlcd_legend_df$id_l2
nlcd_legend_df <- subset(nlcd_legend_df,id_l2%in%lc_df$ID ) 
dim(nlcd_legend_df)

rec_df <- nlcd_legend_df[,c(2,1)]

#r_date1_rec <- subs(r_lc_date1,nlcd_legend_df[,1:2],by="id_l1","id_l2")
r_date1_rec <- subs(r_lc_date1,rec_df,by="id_l2","id_l1")
r_date2_rec <- subs(r_lc_date2,rec_df,by="id_l2","id_l1")

plot(r_date1_rec)

rec_xtab_df <- crosstab(r_date1_rec,r_date2_rec,long=T)
names(rec_xtab_df) <- c("2001","2011","freq")
#View(rec_xtab_df)

### plot urban growth and urban loss?

ncell(r_date1_rec)

### Make this a function:

label_legend_df <- data.frame(ID=nlcd_legend_df$id_l1,name=nlcd_legend_df$name_l1)
r_stack <- stack(r_date1_rec,r_date2_rec)

lc_df <- freq(r_stack,merge=T)
names(lc_df) <- c("value","date1","date2")
lc_df$diff <- lc_df$date2 - lc_df$date1

lc_df <- merge(lc_df,label_legend_df,by.x="value",by.y="ID",all.y=F)
lc_df <- lc_df[!duplicated(test_df),]
barplot(lc_df$diff,names.arg=lc_df$name,las=2)
total_val  <- sum(lc_df$date1)
lc_df$perc_change <- 100*lc_df$diff/total_val 
barplot(lc_df$perc_change,names.arg=lc_df$name,las=2)

View(lc_df)
## Plot the changes here by land cover classes

### reclassify:  

#devopped
r_cat2 <- r_date2_rec==2
r_not_cat2 <- r_date1_rec!=2

r_change <- r_cat2 * r_not_cat2
plot(r_change)
change_tb <- freq(r_change) #this is about 500,000 pixels!!!

# change to urban from 2001 to 2011
# compute rate of growth for a year and project in 2022
# show in plot:

## y= 1 if change to urban over 2001-2011
### Suitability:
#var1: distance to existing urban in 2001
#var2: distance to road in 2001
#var3: elevation, low slope
#var4: landcover before

#Not available yer: var5: conservation areas

#Suitable land cover:
#-meadow is cheap: easier than forest
#-forest ok
#-urban: NA (already urban cannot transition)
#

## could also do a logistic

r_cat2<- r_date1_rec==2
plot(r_cat2)

writeRaster(r_cat2,filename = "developped_2001.tif")
### distance to existing in 2001

r_roads <- raster(file.path(in_dir_var,roads_fname))
#<- "r_roads_Harris.tif"
plot(r_roads)
r_roads_bool <- r_roads >0
roads_bool_fname <- "roads_bool.tif" 
writeRaster(r_roads_bool,filename = roads_bool_fname,overwrite=T)

if(gdal_installed==TRUE){
  
  ## Distance from developped land
  srcfile <- cat_bool_fname 
  
  dstfile_developped <- file.path(out_dir,paste("developped_distance_",out_suffix,file_format,sep=""))
  n_values <- "1"
  
  ### Note that gdal_proximity doesn't like when path is too long
  cmd_developped_str <- paste("gdal_proximity.py",basename(srcfile),basename(dstfile_developped),"-values",n_values,sep=" ")

  ### Distance from roads
  
  srcfile <- roads_bool_fname 
  dstfile_roads <- file.path(out_dir,paste("roads_bool_distance_",out_suffix,file_format,sep=""))
  n_values <- "1"
  
  ### Note that gdal_proximity doesn't like when path is too long
  cmd_roads_str <- paste("gdal_proximity.py",basename(srcfile),
                         basename(dstfile_roads),
                         "-values",n_values,sep=" ")

  sys_os <- as.list(Sys.info())$sysname
  
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
  r_developped_distance <- raster(file.path(in_dir,paste("developped_distance_",file_format,sep="")))
  r_roads_distance <- raster(file.path(in_dir_var,paste("roads_bool_distance",file_format,sep="")))
}

plot(r_developped_distance)
plot(r_roads_distance)

#Now rescale the distance...
min_val <- cellStats(r_roads_distance,min) 
max_val <- cellStats(r_roads_distance,max)

#Linear rescaling:
#y = ax + b with b=0
#with 9 being new max and 0 being new min
a = (1 - 0) /(max_val - min_val)
a=1
r_roads_dist <- r_roads_distance * a
plot(r_roads_dist)

#Get distance from managed land
#b. Which parts of Clay County contain proximity-to-managed-lands characteristics that would make them more favorable to be used as conservation lands?

#min_val <- cellStats(r_developped_distance,min) 
#max_val <- cellStats(r_developped_distance,max)

a = (1 - 0) /(max_val - min_val) #linear rescaling factor
r_developped_dist <- r_developped_distance * a

############ Now deal with elevation

r_elevation <- raster(file.path(in_dir_var,elevation_fname))
#<- "srtm_Houston_area_90m.tif"
r_elevation_30m <- disaggregate(r_elevation,fact=3)
projection(r_elevation_30m)
r_elevation_reg <- projectRaster(r_elevation_30m,r_date1_rec)

### reclass Land cover
#?mask
r_mask <- r_date1_rec==2

NAvalue(r_date1_rec_masked)
r_date1_rec_masked <- mask(r_date1_rec,r_mask,maskvalue=1)
#r_date1_rec[r_date1_rec==2] <- NA

plot(r_date1_rec_masked)

###### The logistic regression comes here:

r_elevation_reg
r_date1_rec_masked
r_roads_dist
r_developped_dist

r_change

### let's split training and testing??
#removing water and developped in 2001
r_mask <- (r_date1_rec!=2)*(r_date1_rec!=1)
plot(r_mask)

plot(r_change)
r_variables <- stack(r_change,r_date1_rec_masked,r_elevation_reg,r_roads_dist,r_developped_dist)
r_variables <- mask(r_variables,mask=r_mask,maskvalue=0)
#plot(r_not_cat2)
names(r_variables) <- c("change","land_cover","elevation","roads_dist","developped_dist")

### May be useful to have x and y locations

variables_df <- na.omit(as.data.frame(r_variables))
#variables_df <- na.omit(variable_df)
dim(variables_df)

variables_df$land_cover <- as.factor(variables_df$land_cover)
variables_df$change <- as.factor(variables_df$change)

names(variables_df)
#names(variables_df) <- c("change","land_cover","elevation","roads_dist","developped_dist")

mod <- glm(change ~ land_cover + elevation + roads_dist + developped_dist, 
           data=variables_df , family=binomial())
mod
r_p <- predict(r_variables, mod, type="response")
plot(r_p)

histogram(r_p)
histogram(p,xlim=c(0,1),breaks=10)

#### Do AUC to check how good it is?
plot(r_change)

plot(r_date1_rec_masked)
tb_freq <- freq(r_date1_rec_masked)
View(tb_freq)

####################### End of script #####################################