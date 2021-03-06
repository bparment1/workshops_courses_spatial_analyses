####################################    Fire Alaska Analyses   #######################################
############################  Analyses using time series and fire locations  #######################################
#This script performs analyses for the Exercise 2 of the Geospatial Course using NDVI data derived from MODIS.
#The overall goal is to explore the impact of fire on surface state variables observed from Remote Sensing.
#Additional data is provided with location of three fire polygons containing pixels affected by fire. 
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/15/2017 
#DATE MODIFIED: 02/25/2019
#Version: 1
#PROJECT: AAG 2019 workshop and Geospatial Data Analysis course at SESYNC 
#TO DO: Update to include sf
#
#COMMIT: Fixing and adding change with reclassify
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
library(readxl) #functionalities to read in excel type data
library(psych) #pca/eigenvector decomposition functionalities
library(sf) #spatial objects and methods

###### Functions used in this script

function_preprocessing_and_analyses <- "fire_alaska_analyses_preprocessing_functions_03102017.R" #PARAM 1
function_analyses <- "exercise2_fire_alaska_analyses_functions_03232017.R" #PARAM 1
script_path <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2019_geospatial_workshop/Exercise_2/scripts" #this must be changed to local path
source(file.path(script_path,function_preprocessing_and_analyses)) #source all functions used in this script 1.
source(file.path(script_path,function_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

in_dir_NDVI <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2019_geospatial_workshop/Exercise_2/data/NDVI_alaska_2005"
in_dir_var <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2019_geospatial_workshop/Exercise_2/data"
out_dir <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2019_geospatial_workshop/Exercise_2/outputs"

infile_ecoreg <- "wwf_terr_ecos_Alaska_ECOREGIONS_ECOSYS_ALB83.shp" #WWF ecoregions 2001 for Alaska
fire_poly_shp_fname <- "OVERLAY_ID_83_399_144_TEST_BURNT_83_144_399_reclassed.shp"

#region coordinate reference system
CRS_reg <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

file_format <- ".tif" #PARAM5
NA_flag_val <- -9999
out_suffix <-"exercise2_02252019" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

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

lf_NDVI <- list.files(path=in_dir_NDVI, pattern="*.tif",full.names=T)
r_NDVI_ts <- stack(lf_NDVI)

plot(r_NDVI_ts)
date_range <- c("2005.01.01","2005.12.31") #NDVI Alaska, year 2005 (this is a 16 days product)
  
#generate dates for 16 days product
dates_val <- generate_dates_by_step(date_range[1],date_range[2],16)$dates #NDVI Katrina

file.info(lf_NDVI[1])$size/(1024*1024) # this is in bytes, convert to mb
dataType(r_NDVI_ts) #Examine the data type used in the storing of data, this is float 32 signed: FLT4S
inMemory(r_NDVI_ts) #Is the data in memory? Raster package does not load in memory automatically.
dim(r_NDVI_ts) #dimension of the raster object: rows, cols, layers/bands

##### PART I: DISPLAY AND EXPLORE DATA ##############

lf_var <- list.files(path=in_dir_var,pattern="*.tif$",full.names=T)
r_var <- stack(lf_var) # create a raster stack, this is not directly stored in memory
dim(r_var) #dimension of the stack with 
plot(r_var)

##### For now using sp, could use sf too
tb_freq <- freq(subset(r_var,6)) #  count of pixels by burn scars
projection(r_var) #note that there is no projection assigned
projection(r_var) <- CRS_reg

ecoreg_spdf <-readOGR(dsn=in_dir_var,sub(".shp","",infile_ecoreg))
proj4string(ecoreg_spdf)
proj4string(ecoreg_spdf) <- CRS_reg
plot(ecoreg_spdf)
spplot(ecoreg_spdf)

# Wrong order in terms of the categories of ecoreg so assign them
cat_names<-c("Alaska Peninsula montane taiga",
             "Alaska St Elias Range tundra",
             "Aleutian Islands tundra",
             "Arctic coastal tundra",
             "Arctic foothills tundra",
             "Beringia lowland tundra",
             "Beringia upland tundra",
             "Brooks-British Range tundra",
             "Cook Inlet taiga",
             "Copper Plateau taiga",
             "Interior Alaska-Yukon lowland taiga",
             "Interior Yukon-Alaska alpine tundra",
             "Northern Cordillera forests",
             "Northern Pacific coastal forests",
             "Ogilvie-Mackenzie alpine tundra",
             "Pacific Coastal Mountain icefields and tundra",
             "Rock and Ice")

ecoreg_spdf$ECO_NAME<-cat_names

nb_col <- length(cat_names)

#problem with extent, ecoreg_spdf is not the same extent as raster images!!
lf_eco<-list.files(pattern="ecoregion_map.rst$")

ecoreg_rast<-rasterize(ecoreg_spdf,r_var,"DATA_VALUE")
data_name<-paste("ecoregion_map",sep="")
raster_name<-paste(data_name,file_format, sep="") #Remember file_format ".tif" is set at the beginning
writeRaster(ecoreg_rast, filename=raster_name,NAflag=NA_flag_val,overwrite=TRUE)  #Writing the data in a raster file format...

###An example of plottting raster data for publication:

## Generate a color palette/color ramp
col_eco <- rainbow(nb_col)
col_eco[16]<-"brown"
col_eco[11]<-"darkgreen"
col_eco[7]<-"lightblue"
col_eco[6]<-"grey"
col_eco[12]<-"yellowgreen"

plot(ecoreg_rast)
plot(ecoreg_rast,col=col_eco)

##Figure 1: wwf ecoregion
res_pix<-960
col_mfrow<-1
row_mfrow<-1
png(filename=paste("Figure1_paper1_wwf_ecoreg_Alaska",out_suffix,".png",sep=""),
    width=col_mfrow*res_pix,height=row_mfrow*res_pix)
#par(mfrow=c(1,2))


plot(ecoreg_rast,col=col_eco,legend=FALSE,axes="FALSE")
legend("topright",legend=cat_names,title="WWF ecoregions",
       pt.cex=1.1,cex=1.1,fill=col_eco,bty="n")
scale_position<-c(450000, 600000)
arrow_position<-c(900000, 600000)

label_scalebar<-c("0","125","250")
scalebar(d=250000, xy=scale_position, type = 'bar', 
         divs=3,label=label_scalebar,below="kilometers",
         cex=1.8)
#this work on non sp plot too
SpatialPolygonsRescale(layout.north.arrow(), offset = arrow_position, 
                       scale = 150000, fill=c("transparent","black"),plot.grid=FALSE)
dev.off()

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
val_diff_mean <- cellStats(r_diff_NDVI,mean)
val_diff_sd <- cellStats(r_diff_NDVI,sd)

r_diff_NDVI_standardized <- (r_diff_NDVI - val_diff_mean)/val_diff_sd
plot(r_diff_NDVI_standardized,col=matlab.like(255),main='Standardized NDVI Difference')#
plot(r_diff_NDVI_standardized,col=matlab.like(255),zlim=c(-10,5),main='Standardized NDVI Difference')#

hist_bins <- seq(-15,15,by=0.5)
hist(r_diff_NDVI_standardized,breaks=hist_bins,main='Standardized NDVI Difference')
hist(r_diff_NDVI_standardized,breaks=hist_bins,xlim=c(-5,5),main='Standardized NDVI Difference') # zoom in

hist(r_diff_NDVI_standardized, main='Standardized NDVI Difference')

##Generate change map using thesholding of standardized NDVI
r_change_NDVI_pos <-  r_diff_NDVI_standardized > 1.96
r_change_NDVI_neg <-  r_diff_NDVI_standardized < -1.96

cat_names=c("No Change","Change")
col_palette <- c("black","red")

plot(r_change_NDVI_pos,main='Change: increases in NDVI',
     legend=FALSE,axes="FALSE",col=col_palette)
legend("topright",legend=cat_names,title="NDVI change",
            pt.cex=1.1,cex=1.1,fill=col_palette,bty="n")

plot(r_change_NDVI_neg,main='Change: decreases in NDVI',
     legend=FALSE,axes="FALSE",col=col_palette)
legend("topright",legend=cat_names,title="NDVI change",
       pt.cex=1.1,cex=1.1,fill=col_palette,bty="n")

#Use reclassify to combine change or on larger raster images
rec_val <- c(-10000, -1.96, 1, #if  
                -1.96, 1.96, 0,  
                 1.96, 10000, 1)
rclmat <- matrix(rec_val, ncol=3, byrow=TRUE)
r_change_NDVI <- reclassify(r_diff_NDVI_standardized, 
                            rclmat)
plot(r_change_NDVI,main='Change in NDVI',
     legend=FALSE,axes="FALSE",col=col_palette)
legend("topright",legend=cat_names,title="NDVI change",
       pt.cex=1.1,cex=1.1,fill=col_palette,bty="n")

writeRaster(r_change_NDVI_pos,"r_change_NDVI_pos_196.tif",overwrite=T)
writeRaster(r_change_NDVI_neg,"r_change_NDVI_neg_196.tif",overwrite=T)

# Better way: use names with output suffix to take changes into account
out_filename <- paste0("r_change_NDVI_",out_suffix,file_format)
writeRaster(r_change_NDVI,out_filename,overwrite=T)

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

###############################################
######## Let's carry out a PCA in T-mode #######

#Correlate long term mean to PC!
cor_mat_layerstats <- layerStats(r_NDVI_ts, 'pearson', na.rm=T)
cor_matrix <- cor_mat_layerstats$`pearson correlation coefficient`
class(cor_matrix)
dim(cor_matrix)
print(cor_matrix) #use View(cor_matrix in Rstudio)
image(cor_matrix)

#### Carry out a PCA using principal from psych package
pca_mod <-principal(cor_matrix,nfactors=3,rotate="none")
class(pca_mod$loadings)
str(pca_mod$loadings)
plot(pca_mod$loadings[,1],type="b",
     xlab="time steps",
     ylab="PC loadings",
     ylim=c(-1,1),
     col="blue")
lines(pca_mod$loadings[,2],type="b",col="red")
lines(pca_mod$loadings[,3],type="b",col="black")
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
plot(pca_mod$values,
     main="Scree plot: Variance explained",
     type="b",
     xlab="Principal Components",
     ylab="Eigen values")

### Generate scores from eigenvectors
## Do it two different ways:

#By converting data and working with matrix:
df_NDVI_ts <- as(r_NDVI_ts,"SpatialPointsDataFrame")
df_NDVI_ts <- as.data.frame(df_NDVI_ts) #convert to data.frame because predict works on data.frame
names(df_NDVI_ts) <- c(paste0("pc_",seq(1,23,1)),"x","y")
names(df_NDVI_ts)

pca_all <- as.data.frame(predict(pca_mod,df_NDVI_ts[,1:23])) ## Apply model object on the data.frame

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
r_pca <- predict(r_NDVI_ts, pca_mod, index=1:3,filename="pc_scores.tif",overwrite=T) # fast
plot(r_pca)

plot(stack(r_pc1,r_pc2))
#layerStats(r_pc1,r_NDVI_mean )
cor_pc <- layerStats(stack(r_pc1,r_NDVI_mean),'pearson', na.rm=T)
cor_pc #PC1 correspond to the average mean by pixel as expected.
plot(r_pc2)

## Use animate function
#animate(r_NDVI_ts, pause=0.25, n=1)

#### Your turn: 
#Compute the average using the ecoregions as zonal areas for PC1 and PC2
#Plot averages in the PC1-PC2. Are these variables useful to separate the various ecoregions?
#Correlate average times series for PC2 to the loadings. What does it tell you?

################################ End of Script ######################################
