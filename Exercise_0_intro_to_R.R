#################################### An Introduction to R #######################################

#This script was distributed during the Spatial Filtering Workshop in 2008
#Modified substantially by Benoit Parmentier on April 4, 2017 
#for AAG2017 workshop.

########################### PART 1: R BASICS ################
 
####################
## Numeric vectors

x <- c(1,2,3) #assigning values to x (numeric vector)
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
lu.char <- sample(landuse,7,replace=T)
lu.char
lu.fact <- factor(lu.char)
class(lu.fact)
lu.fact

####################
## matricies

mat <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=2)
mat <- matrix(1:6, nrow=2, ncol=3)
mdat <- matrix(seq(1,9,1), nrow=3)
mdat
mdat %*% mdat #matrix multiplication
mdat[1:2,]
diag(mdat)

####################
## lists

my.cats <- c("Henry","Lily","Charlie")
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
my.data$age
my.data[,2]
my.data[,"age"]

# attach/detach
age
attach(my.data) #not advised!! Benoit 07/15
age
name
detach(my.data)
age
searchpaths()      # Search paths will be displayed
rm(my.data)

####################
## read data
in_dir <- "/Users/benoitparmentier/Dropbox/GeospatialPatternCourse/Module3/Intro_to_R/data"
#setwd("c:\\workspace") #Windows style...
setwd(in_dir) 
getwd() #Get the current workind directory
library(spdep)
data(columbus)
write.table(columbus,"columbus.csv",sep=",",row.names=F)
columbus <- read.csv("columbus.csv", header=T)
colnames(columbus)
library(foreign)
concord_data <- read.spss("Concord1.sav", to.data.frame = TRUE)
concord_data[1:2,c(3,8)]

####################
## data management
tr <- read.dbf("turnout.dbf")
tr.sub <- tr[,c("FID_1","GDPCAP","TURNOUT01")] #subset a data.frame
nrow(tr)
summary(tr) #summary of the data.frame
tr.sub$log_GDPCAP <- log(tr.sub[,2]) #creating column and assigning new values using log of existing column
colnames(tr.sub)
tr.sub <- tr.sub[,-2] #drop column 2
colnames(tr.sub)
tr.sub$log_GDPCAP <- NULL #remove column

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
## linear model

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

history()
savehistory(file="IntroR.Rhistory")

###################PART 2: Spatial R ##########

#to install a new package use:
#install.package("package_name") eg install.package("colorRamps")

##########################
#create a raster image:
library(raster)

r <- raster()
values(r) <- rnorm(ncell(r))
plot(r)
projection(r) #display the current projection system
writeRaster(r,"test.tif",overwrite=TRUE)

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
infile_shp <- "/Users/benoitparmentier/Dropbox/GeospatialPatternCourse/Module3/srm_data_v2/turnout.shp"

layer_path <- dirname(infile_shp) #extract path from string
layer_name <- gsub(".shp","",basename(infile_shp)) #extract layer name
tr <- readOGR(layer_path,layer_name)

spplot(tr,"TURNOUT01",col.regions=matlab.like(25),
       main="TURNOUT") #plot turnout using specific color ramp

### save figure
out_suffix <- "wm_07152014"
png(filename=paste("Figure_turnout_Italy","_",out_suffix,".png",sep=""))
spplot(tr,"TURNOUT01",col.regions=matlab.like(25),
       main="TURNOUT") #plot turnout using specific color ramp
dev.off()

##### 

library(dismo)
mymap <- gmap("Belgium")  # choose whatever country
plot(mymap)

################ END OF SCRIPT ###################