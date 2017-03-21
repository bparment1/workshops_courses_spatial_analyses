####################################    Spatial Analyses: SYRACUSE   #######################################
#######################################  Analyse data from Census #######################################
#This script performs basic analyses for the Exercise 1 of the workshop using Census data.
# The overall goal is to explore spatial autocorrelation and aggregation of units of analyses.     
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/21/2017 
#DATE MODIFIED: 03/22/2017
#Version: 1
#PROJECT: AAG 2017 workshop preparation
#TO DO:
#
#COMMIT: initial commit for exercise 1, AAG workshop
#
#################################################################################################

###Loading R library and packages                                                      

library(sp) # spatial/geographfic objects and functions
library(rgdal) #GDAL/OGR binding for R with functionalities
library(spdep) #spatial analyses operations, functions etc.
library(gtools) # contains mixsort and other useful functions
library(maptools) # 
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
library(BMS) #contains hex2bin and bin2hex
library(bitops) #
#
library(foreign) # import datasets from SAS, spss, stata and other sources
library(gdata) #read xls, dbf etc.
library(classInt) #methods to generate class limits
#library(sqldf) #Not available for 3.3.3
library(plyr)

###### Functions used in this script

function_preprocessing_and_analyses <- "fire_alaska_analyses_preprocessing_functions_03102017.R" #PARAM 1
script_path <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/R_scripts"
source(file.path(script_path,function_preprocessing_and_analyses)) #source all functions used in this script 1.

#####  Parameters and argument set up ###########

in_dir_var <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Exercise_1/data"
out_dir <- "/home/bparmentier/Google Drive/Data/Seminars_talks_workshops/workshops/AAG2017_spatial_temporal_analysis_R/Exercise_1/outputs"

#CRS_reg <- CRS_WGS84 # PARAM 4

file_format <- ".tif" #PARAM5
NA_value <- -9999 #PARAM6
NA_flag_val <- NA_value #PARAM7
out_suffix <-"exercise1_03212017" #output suffix for the files and ouptu folder #PARAM 8
create_out_dir_param=TRUE #PARAM9

################# START SCRIPT ###############################

### PART I READ AND PREPARE DATA FOR REGRESSIONS #######

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

ct_2000_fname <- "ct_00.shp" # CT_00: Cencus Tracts 2000
bg_2000_fname <- "bg_00.shp" # BG_00: Census Blockgroups 2000
bk_2000_fname <- "bk_00.shp" # BK_00: Census Blocks 2000

census_table_fname <- "census.csv" #contains data from census to be linked
soil_PB_table_fname <- "Soil_PB.csv" #same as census table
tgr_shp_fname <- "tgr36067lkA.shp" #contains data from census to be linked

metals_table_fname <- "SYR_metals.xlsx" #contains metals data to be linked

ct_2000_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(ct_2000_fname)))
bg_2000_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(bg_2000_fname)))
bk_2000_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(bk_2000_fname)))

#roads
#tgr_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(tgr_shp_fname)))
#View(tgr_sp)
#tgr_shp_fname

census_syr_df <- read.table(file.path(in_dir_var,census_table_fname),sep=",",header=T)
metals_df <- read.xls(file.path(in_dir_var,metals_table_fname),sep=",",header=T)

#Soil lead samples: UTM z18 coordinates
soil_PB_df <- read.table(file.path(in_dir_var,soil_PB_table_fname),sep=",",header=T) #point locations

dim(census_syr_df) #47 spatial entities
dim(ct_2000_sp) #47 spatial entities
dim(metals_df) #47 entities
dim(bg_2000_sp) #147 spatial entities

######## PRODUCE MAPS OF 2000 Population #########

#First need to link it.

names(bg_2000_sp)
names(census_syr_df)
#key is "TRACT" but with a different format.
#First fix the format
head(bg_2000_sp)
head(census_syr_df$BKG_KEY)
#as.numeric(as.character(ct_2000_sp$TRACT))
ct_2000_sp$TRACT <- as.numeric(as.character(ct_2000_sp$TRACT))

bg_2000_sp <- merge(bg_2000_sp,census_syr_df,by="BKG_KEY")

spplot(bg_2000_sp,"POP2000",main="POP2000") #quick visualization of population 

##Now change the classes!

### Summarize by census track
census_2000_sp <- aggregate(bg_2000_sp , by="TRACT",FUN=sum)
##compare to sp!!
df_test <- aggregate(POP2000 ~ TRACT, bg_2000_sp , FUN=sum)

### Check if the new geometry of entities is the same as census
plot(census_2000_sp)
plot(ct_2000_sp,border="red",add=T)
nrow(census_2000_sp)==nrow(ct_2000_sp)

df_summary_by_census <- aggregate(. ~ TRACT, bg_2000_sp , FUN=sum) #aggregate all variables from the data.frame

##Join by key table id:
dim(ct_2000_sp)
ct_2000_sp <- merge(ct_2000_sp,df_summary_by_census,by="TRACT")
dim(ct_2000_sp)

#save as sp and text table
#write.table(file.path(out_dir,)

### Do another map with different class intervals: 

title_str <- "Population by census tract in 2000"
range(ct_2000_sp$POP2000)
col_palette <- matlab.like(256)

## 9 classes with fixed and constant break
break_seq <- seq(0,9000,1000)
breaks.qt <- classIntervals(ct_2000_sp$POP2000, n=length(break_seq), 
                            style="fixed", fixedBreaks=break_seq, intervalClosure='right')

## Color for each class
#colcode = findColours(breaks.qt , c('darkblue', 'blue', 'lightblue', 'palegreen','yellow','lightpink', 'pink','brown3',"red","darkred"))
p_plot_pop2000_ct <- spplot(ct_2000_sp,
                            "POP2000",
                            col="transparent", #transprent color boundaries for polygons
                            col.regions = col_palette ,
                            main=title_str,
                            at = breaks.qt$brks)
print(p_plot_pop2000_ct)

### Another map:

breaks.qt <- classIntervals(ct_2000_sp$POP2000, n = 6, style = "quantile", intervalClosure = "right")

p_plot_pop2000_ct <- spplot(ct_2000_sp,
                            "POP2000",
                            col="transparent", #transprent color boundaries for polygons
                            col.regions = col_palette,
                            main=title_str,
                            at = breaks.qt$brks)
print(p_plot_pop2000_ct)

## plot polygons with colors and legend
#plot(ct_2000_sp, col=colcode,border="red")

#title(title_str)
#legend('topleft', legend=c(names(attr(colcode, 'table')),'no data'), 
#       fill=c(attr(colcode, 'palette'),'white'), title=title_str)

##### PART II: Vulnerability to metals #############

#Examine the relationship between metals, Pb and vulnerable populations in Syracuse

#soil_PB_df <- read.table(file.path(in_dir_var,census_table_fname),sep=",",header=T)
metals_df <- read.xls(file.path(in_dir_var,metals_table_fname),sep=",",header=T)

#View(soil_PB_df)
View(metals_df)

##This suggests matching to the following spatial entities
nrow(metals_df)==nrow(ct_2000_sp)
#nrow(soil_PB_df)==nrow(bg_2000_sp)

#dim(bg_2000_sp)
census_metals_sp <- merge(ct_2000_sp,metals_df,by.x="TRACT",by.y="ID")

########processing lead data
### Now let's plot lead data 
#Soil lead samples: UTM z18 coordinates
soil_PB_df <- read.table(file.path(in_dir_var,soil_PB_table_fname),sep=",",header=T) #point locations

proj4string(census_metals_sp) #
names(soil_PB_df)
names(soil_PB_df) <- c("x","y","ID","ppm") 
soil_PB_sp <- soil_PB_df
class(soil_PB_df)
coordinates(soil_PB_sp) <- soil_PB_sp[,c("x","y")]
class(soil_PB_sp)
proj4string(soil_PB_sp) <- proj4string(census_metals_sp)
dim(soil_PB_sp)
soil_PB_sp <- soil_PB_sp[,c("ID","ppm","x","y")]
View(soil_PB_sp)

plot(census_metals_sp)
plot(soil_PB_sp,add=T)


###### Spatial query: associate points of pb measurements to each census tract
### Get the ID and 
soil_tract_id_df <- over(soil_PB_sp,census_2000_sp,fn=mean)
soil_PB_sp <- intersect(soil_PB_sp,census_2000_sp)
#test4 <- gIntersection(soil_PB_sp,census_2000_sp,byid=T)
head(soil_PB_sp$ID)==head(soil_PB_sp$ID)
names(soil_PB_sp)
soil_PB_sp <- rename(soil_PB_sp, c("d"="TRACT")) #from package plyr

census_pb_avg <- aggregate(ppm ~ TRACT,(soil_PB_sp),FUN=mean)

##

### Aggregate by county
#http://gis.stackexchange.com/questions/137621/join-spatial-point-data-to-polygons-in-r
#test_df2 <- over(census_2000_sp,soil_PB_sp,fn=mean)


statesAg2 <- aggregate(tweets["actor.friendsCount"], by = states, mean)



#P2- SPATIAL AND NON SPATIAL QUERIES (cannot use spatial join)
#GOAL: Answer a set of questions using spatial and attribute queries and their combinations

#Produce:
#  a) two different maps based on two different definitions that answer the question:  which areas have high levels of children and are predominantly minority AND are at risk of heavy metal exposure using at least three variables. Use only tabular operations
#b) Same question as a) but using both spatial and tabular operations

#Note: In both cases include the method, variables used and your definition of risk areas in each 4 maps. The definition of risk is your own, you can also follow an established standard that would make sense or is official.  
#From these products, the layman should be able to answer the following questions:
#  a. Where are the areas of high heavy metal exposure that also have high levels of children population that belong to a demographic minority(s)? 
#b. Is there a different outcome in using tabular methods only vs combining tabular and spatial query methods?

########### PART III: Spatial Moran's I and regression ###################



###################### END OF SCRIPT #####################


