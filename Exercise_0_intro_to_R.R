#################################### An Introduction to R #######################################
#This script provides a gentle introduction to some basic commands and tools for geospatial analyses, geospatial processing
#and geovisualization.
#
#Original script from spatial filtering workshop substantially modified.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/15/2017 
#DATE MODIFIED: 02/26/2019
#
#Version: 2
#PROJECT: AAG 2019 workshop  
#TO DO: 
#
#COMMIT: 
#
#################################################################################################

########################### PART 1: R BASICS ################

####################
## Numeric vectors

x <- c(7,13,20) #assigning values to x (numeric vector)
x #displaying the value of x
mean(x) #mean function works on vectors
y <- x-mean(x)
y
z <- y > 0
z             # logical

####################
## Character vectors

landuse <- c("urban","forest","water")
landuse
landuse[1]
lu.char <- sample(landuse,7,replace=T) #use this to draw randomly values
lu.char
lu.fact <- factor(lu.char) # create a factor/categorical variable
class(lu.fact)
lu.fact

####################
## matricies

mat <- matrix(c(1,2,3,4,5,6), 
              nrow=3, 
              ncol=2)

mat <- matrix(1:6, 
              nrow=2, 
              ncol=3)
mdat <- matrix(seq(1,9,1), 
               nrow=3)
mdat
mdat %*% mdat #matrix multiplication
mdat[1:2,]
diag(mdat)
image(mdat) #visualize matrix as an image

####################
## lists

my.cats <- c("Loulou","Lily","Charlie") 
age.cats <- c(5,15,4)
my.list <- list(catname=my.cats, age=age.cats)
my.list
my.list$catname
my.list$age[1]

####################
## data frame

my.data <- data.frame(name=my.cats,age=age.cats)
my.data
dim(my.data) #show the dimension of the data.frame (columns and rows)
ls()       # you may have more objects
my.data$age #select a specific column to view data "$"
my.data[,2] #select column 2
my.data[,"age"] #select by name

# attach/detach
#age #not found
attach(my.data) #not advised!! 
age
name
detach(my.data)
#age #not found
searchpaths()      # Search paths will be displayed
rm(my.data)

####################
## read data and dealing with paths

in_dir <- "/nfs/bparmentier-data/Data/workshop_spatial/sesync2019_geospatial_workshop/Exercise_0/data"

#setwd("c:\\workspace") #Windows style...
setwd(in_dir) 
getwd() #Get the current workind directory
library(spdep)
#if not installed:
#install.packages("spdep")

data(columbus) #load columbus dataset from spdep package
dim(columbus) # size of columbus in rows x columns
head(columbus) #view first six rows of data.frame

write.table(columbus,"columbus.csv",sep=",",row.names=F)
columbus <- read.csv("columbus.csv", header=T)
colnames(columbus)
library(foreign) # to read a variety of data types
concord_data <- read.spss("Concord1.sav", to.data.frame = TRUE)
concord_data[1:2,c(3,8)]

####################
## data management
tr <- read.dbf("turnout.dbf") #read in turnout data from Italian elections

tr.sub <- tr[,c("FID_1","GDPCAP","TURNOUT01")] #subset a data.frame
nrow(tr)
summary(tr) #summary of the data.frame
tr.sub$log_GDPCAP <- log(tr.sub[,2]) #creating column and assigning new values using log of existing column
colnames(tr.sub)
tr.sub <- tr.sub[,-2] #drop column 2
colnames(tr.sub)
tr.sub$log_GDPCAP <- NULL #remove column

## Other common solutions: select and filter with dplyr
library(dplyr)
tr.sub2 <- select(tr,c("FID_1","GDPCAP","TURNOUT01","REGIONE")) #subset a data.frame
class(tr.sub2)
tr.sub2 <- filter(tr.sub2,REGIONE=='Veneto')

##########################
## descriptive statistics

mean(tr.sub$TURNOUT01)
cor(tr.sub) #correlation matrix!!
x <- tr.sub[,2]; sd(x)
sqrt(sum((x - mean(x))^2)/(length(x)-1))
summary(tr.sub)     # summary of each variable

##########################
## graphics

# histogram
attach(tr)
hist(GDPCAP)	# histogram

plot(GDPCAP, TURNOUT01, main="Scatterplot")  #scatterplot
pairs(tr.sub)            # scatterplot matrix

qqnorm(GDPCAP)            # normal quantile plot: qqplot
qqline(GDPCAP)            # add normal quantile line

boxplot(GDPCAP~SOUTH,main="GDP PER REGION")    #box plot
detach(tr)

##########################
## Generate a linear model

# linear regression
mod_lm <- lm(TURNOUT01 ~ GDPCAP, data=tr)
summary(mod_lm)

# graphical diagnostics
plot(fitted(mod_lm),residuals(mod_lm))
abline(h=0,col="4")
title("Scatterplot of fitted values vs. residuals")

qqnorm(resid(mod_lm)) #plotting residuals qqplot
qqline(resid(mod_lm)) # 1 to 1 line for qqplot

#prepare for plotting the model object
par(mfrow=c(2,2))
plot(mod_lm)

##########################
## save and load results

save(mod_lm, file="mod_lm.RData") #save R object 
load("mod_lm.RData") #reload R object

# export data
tr$ols_pred <- fitted(mod_lm) 
write.table(tr, file="tr_mod_lm.txt", sep=",")
write.csv(tr, file="tr_mod_lm.csv")

##########################
## getting help

help(rnorm)
?lm

args(rnorm)
example(lm)

################### PART 2: Spatial R ##########

#to install a new package use:
#install.packages("package_name") eg install.packages("colorRamps")

##########################
#create a raster image:
library(raster)

r <- raster() #create a raster image, default is latitude/longitude degrees 
values(r) <- rnorm(ncell(r)) #assign random values to the raster
plot(r) #quick plot of the raster image
projection(r) #display the current projection system
writeRaster(r,"test.tif",overwrite=TRUE) #save ther aster layer
r  #check the layer

#read a raster
r1 <- raster("moore_reg_mosaiced_MOD13A2_A2005209__005_1_km_16_days_NDVI_09242013_09242013_04062014.tif")
plot(r1)
r_mean <- cellStats(r1,stat="mean")
r_sd <- cellStats(r1,stat="sd")

r2 <- (r1 - r_mean)/r_sd #standardized raster using mean and sd
library(colorRamps) #contains matlab.like color palette)
plot(r2,col=matlab.like(25))
r_outliers_pos <- r2 > 1.96 #determine positive outliers
r_outliers_neg <- r2 < -1.96 #determine positive outliers
freq(r_outliers_pos)
freq(r_outliers_neg)

#########################
## Read in Polygon data eg ESRI shapefile

#shapefile of Italy regions
in_dir <- getwd() #this should be ../data
infile_shp <- file.path(in_dir,"turnout.shp") #joint path and name of input file

## Use sp spatial objects
library(rgdal)
layer_path <- dirname(infile_shp) #extract path from string
layer_name <- gsub(".shp","",basename(infile_shp)) #extract layer name, remove extension ".shp"
tr_sp <- readOGR(layer_path,layer_name) # read in the shapefile

## Use sf spatial objects
library(sf)
## use sf functionality:
tr_sf <- st_read(infile_shp) # read directly in the shapefile
tr_sf

#Plotting with sp: quick plot of vote turnout by Italian regions
spplot(tr_sp,"TURNOUT01",col.regions=matlab.like(25),
       main="TURNOUT") #plot turnout using specific color ramp
### Plotting with sf:
plot(tr_sf["TURNOUT01"],
       main="TURNOUT") #plot turnout using specific color ramp

### save figure in png format
out_suffix <- "workshop_04032019"
png(filename=paste("Figure_turnout_Italy","_",out_suffix,".png",sep=""))
plot(tr_sf["TURNOUT01"],
     main="TURNOUT") #plot turnout using specific color ramp
dev.off()


################ END OF SCRIPT ###################