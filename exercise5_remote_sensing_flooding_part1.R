####################################    Flood Mapping Analyses   #######################################
############################  Analyze and map flooding from RITA hurricane  #######################################
#This script performs analyses for the Exercise 5 of the Short Course using reflectance data derived from MODIS.
#The goal is to map flooding from RITA using various reflectance bands from Remote Sensing platforms.
#Additional data is provided including FEMA flood region. 
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/05/2018 
#DATE MODIFIED: 03/29/2018
#Version: 1
#PROJECT: SESYNC and AAG 2018 workshop/Short Course preparation
#TO DO:
#
#COMMIT: PCA loading space
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

in_dir_reflectance <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_5/data/reflectance_RITA"
in_dir_var <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_5/data"
out_dir <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2018_workshop/Exercise_5/outputs"

#region coordinate reference system
#http://spatialreference.org/ref/epsg/nad83-texas-state-mapping-system/proj4/
CRS_reg <- "+proj=lcc +lat_1=27.41666666666667 +lat_2=34.91666666666666 +lat_0=31.16666666666667 +lon_0=-100 +x_0=1000000 +y_0=1000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" 
file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"exercise5_03292018" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9
method_proj_val <- "bilinear" # "ngb"
multiband <- TRUE #This is only used for multiband products?

### Input data files used:
infile_reg_outline <- "new_strata_rita_10282017.shp"
ref_rast_name <- "r_ref_Houston_RITA.tif"
infile_modis_bands_information <- "df_modis_band_info.txt"
nlcd_2006_filename <- "nlcd_2006_RITA.tif"
infile_name_nlcd_legend <- "nlcd_legend.txt"
#34: 2005-09-22 2005265
infile_reflectance_date1 <- "mosaiced_MOD09A1_A2005265__006_reflectance_masked_RITA_reg_1km.tif"
#35: 2005-09-30 2005273
infile_reflectance_date2 <- "mosaiced_MOD09A1_A2005273__006_reflectance_masked_RITA_reg_1km.tif"

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

##### PART I: DISPLAY AND EXPLORE DATA ##############

#lf_var <- list.files(path=in_dir_var,pattern="*.tif$",full.names=T)
#r_var <- stack(lf_var) # create a raster stack, this is not directly stored in memory

#dim(r_var) #dimension of the stack with 
#plot(r_var)

## #34: 2005265 this is Sept 22
## #35: 2005273 this is Sept 30
r_before <- brick(file.path(in_dir_var,infile_reflectance_date1)) # <- "mosaiced_MOD09A1_A2005265__006_reflectance_masked_RITA_reg_1km.tif"
r_after <- brick(file.path(in_dir_var,infile_reflectance_date2)) # <- "mosaiced_MOD09A1_A2005265__006_reflectance_masked_RITA_reg_1km.tif"

plot(r_before)

reg_sf <- st_read(file.path(in_dir_var,infile_reg_outline))
reg_sf <- st_transform(reg_sf,
                       crs=CRS_reg)
reg_sp <-as(reg_sf, "Spatial") 
plot(reg_sf$geometry)

r_ref <- rasterize(reg_sp,
                   r_before,
                   field="OBJECTID_1",
                   fun="first")
plot(r_ref) # zone 2 is flooded and zone 1 is not flooded

#### Let's examine a location within FEMA flooded zone and outside: use centroids
centroids_sf <- st_centroid(reg_sf)

df_before <- extract(r_before,centroids_sf)
df_after <- extract(r_after,centroids_sf)

### Plot values of bands before and after for flooded region:
plot(df_before[2,],type="l")
lines(df_after[2,],col="red")

## Read band information since it is more informative!!
df_modis_band_info <- read.table(file.path(in_dir_var,infile_modis_bands_information),
                                 sep=",",
                                 stringsAsFactors = F)
print(df_modis_band_info)

names(r_before) <- df_modis_band_info$band_name
names(r_after) <- df_modis_band_info$band_name

### Order band in terms of wavelenth:

df_modis_band_info <- df_modis_band_info[order(df_modis_band_info$start_wlength),]
#SWIR1 (1230–1250 nm), SWIR2 (1628–1652 nm) and SWIR3 (2105–2155 nm).
band_refl_order <- df_modis_band_info$band_number

plot(df_before[2,band_refl_order],type="l",main="Reflectance profile for centroid of flooded area")
lines(df_after[2,band_refl_order],col="red")
# Add legend

###### Now do a extraction for nlcd data
nlcd2006_reg <- raster(file.path(in_dir_var,nlcd_2006_filename))

avg_nlcd <- as.data.frame(zonal(r_before,nlcd2006_reg,fun="mean"))

lc_legend_df <- read.table(file.path(in_dir_var,infile_name_nlcd_legend),
                           stringsAsFactors = F,
                           sep=",")
print(avg_nlcd)

names(lc_legend_df)
### Add relevant categories
lc_legend_df_subset <- subset(lc_legend_df,select=c("ID","NLCD.2006.Land.Cover.Class")]
avg_nlcd_test <- merge(avg_nlcd,lc_legend_df_subset,by.x="zone",by.y="ID",all.y=F)
View(avg_nlcd_test)
lc_df <- lc_df[!duplicated(lc_df),] #remove duplictates
head(lc_df) # Note the overall cahnge

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

##### plot feature space:

df_raster_val <- as.data.frame(stack(r_after,r_pca,nlcd2006_reg))

View(df_raster_val)

#Water
plot(df_raster_val[df_raster_val$nlcd_2006_RITA==11,c("pc_scores.1")],
       df_raster_val[df_raster_val$nlcd_2006_RITA==11,c("pc_scores.2")],
       col="blue",cex=0.15)

#Urban: dense
points(df_raster_val[df_raster_val$nlcd_2006_RITA==22,c("pc_scores.1")],
       df_raster_val[df_raster_val$nlcd_2006_RITA==22,c("pc_scores.2")],
       col="brown",cex=0.15)

#Forest:
points(df_raster_val[df_raster_val$nlcd_2006_RITA==42,c("pc_scores.1")],
     df_raster_val[df_raster_val$nlcd_2006_RITA==42,c("pc_scores.2")],
     col="green",cex=0.15)

### Note the negative values are related to Forest on PC2

#Water
plot(df_raster_val[df_raster_val$nlcd_2006_RITA==11,c("pc_scores.1")],
     df_raster_val[df_raster_val$nlcd_2006_RITA==11,c("pc_scores.2")],
     col="blue",cex=0.15,
     ylim=c(-2,2),
     xlim=c(-2,2))

#Urban: dense
points(df_raster_val[df_raster_val$nlcd_2006_RITA==22,c("pc_scores.1")],
       df_raster_val[df_raster_val$nlcd_2006_RITA==22,c("pc_scores.2")],
       col="brown",cex=0.15)

#Forest:
points(df_raster_val[df_raster_val$nlcd_2006_RITA==42,c("pc_scores.1")],
       df_raster_val[df_raster_val$nlcd_2006_RITA==42,c("pc_scores.2")],
       col="green",cex=0.15)


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