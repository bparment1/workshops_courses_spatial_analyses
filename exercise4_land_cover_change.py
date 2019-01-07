# -*- coding: utf-8 -*-
"""
Spyder Editor.
"""
#################################### Land Use and Land Cover Change #######################################
############################ Analyze Land Cover change in Houston #######################################
#This script performs analyses for the Exercise 4 of the AAG Course using aggregated NLCD values.
#The goal is to assess land cover change using two land cover maps in the Houston areas.
#Additional datasets are provided for the land cover change modeling. A model is built for Harris county.
#
#AUTHORS: Benoit Parmentier
#DATE CREATED: 01/07/2019
#DATE MODIFIED: 01/07/2019
#Version: 1
#PROJECT: AAG 2019 Geospatial Short Course
#TO DO:
#
#COMMIT: clean up code for workshop
#
#################################################################################################
	
	
###### Library used in this script

import gdal
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import rasterio
import subprocess
import pandas as pd
import os, glob
from rasterio import plot
import geopandas as gpd
import descartes
import pysal as ps
from cartopy import crs as ccrs
from pyproj import Proj
from osgeo import osr
from shapely.geometry import Point

################ NOW FUNCTIONS  ###################

##------------------
# Functions used in the script 
##------------------

def create_dir_and_check_existence(path):
    #Create a new directory
    try:
        os.makedirs(path)
    except:
        print ("directory already exists")

def open_image(url):
    image_data = open_http_query(url)
    
    if not image_data:
            return None
            
    mmap_name = "/vsimem/"+uuid4().get_hex()
    gdal.FileFromMemBuffer(mmap_name, image_data.read())
    gdal_dataset = gdal.Open(mmap_name)
    image = gdal_dataset.GetRasterBand(1).ReadAsArray()
    gdal_dataset = None
    gdal.Unlink(mmap_name)
    
    return image

############################################################################
#####  Parameters and argument set up ########### 

#ARGS 1
in_dir = "/home/bparmentier/c_drive/Users/bparmentier/Data/python/Exercise_4/data"
#ARGS 2
out_dir = "/home/bparmentier/c_drive/Users/bparmentier/Data/python/Exercise_4/outputs"
#ARGS 3:
create_out_dir=True #create a new ouput dir if TRUE
#ARGS 7
out_suffix = "exercise4_01072018" #output suffix for the files and ouptut folder
#ARGS 8
NA_value = -9999 # number of cores
file_format = ".tif"

#NLCD coordinate reference system: we will use this projection rather than TX.
CRS_reg <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
	
ct_2000_fname = "ct_00.shp" # CT_00: Cencus Tracts 2000
bg_2000_fname = "bg_00.shp" # BG_00: Census Blockgroups 2000
bk_2000_fname = "bk_00.shp" # BK_00: Census Blocks 2000

census_table_fname = "census.csv" #contains data from census to be linked
soil_PB_table_fname = "Soil_PB.csv" #same as census table
tgr_shp_fname = "tgr36067lkA.shp" #contains data from census to be linked

metals_table_fname = "SYR_metals.xlsx" #contains metals data to be linked

################# START SCRIPT ###############################

######### PART 0: Set up the output dir ################

#set up the working directory
#Create output directory

if create_out_dir==True:
    #out_path<-"/data/project/layers/commons/data_workflow/output_data"
    out_dir = "output_data_"+out_suffix
    out_dir = os.path.join(in_dir,out_dir)
    create_dir_and_check_existence(out_dir)
    os.chdir(out_dir)        #set working directory
else:
    os.chdir(create_out_dir) #use working dir defined earlier
    
    
#######################################
### PART 1: Read in DATA #######

#ct_2000_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(ct_2000_fname)))
#bg_2000_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(bg_2000_fname)))
#bk_2000_sp <- readOGR(dsn=in_dir_var,sub(".shp","",basename(bk_2000_fname)))

## Counties for Syracuse in 2000
ct_2000_filename = os.path.join(in_dir,ct_2000_fname)
## block groups for Syracuse in 2000
bg_2000_filename = os.path.join(in_dir,bg_2000_fname)
## block for Syracuse in 200
bk_2000_filename = os.path.join(in_dir,bk_2000_fname)

 
ct_2000_gpd = gpd.read_file(ct_2000_filename)
bg_2000_gpd = gpd.read_file(bg_2000_filename)
bk_2000_gpd = gpd.read_file(bk_2000_filename)

ct_2000_gpd.describe()
ct_2000_gpd.plot(column="CNTY_FIPS")
ct_2000_gpd.head()

metals_df = pd.read_excel(os.path.join(in_dir,metals_table_fname))
census_syr_df = pd.read_csv(os.path.join(in_dir,census_table_fname),sep=",",header=0)

#Soil lead samples: UTM z18 coordinates
soil_PB_df = pd.read_csv(os.path.join(in_dir,soil_PB_table_fname),sep=",",header=None) #point locations

census_syr_df.shape #57 spatial entities
ct_2000_gpd.shape #57 spatial entities
metals_df.shape #57 entities
bg_2000_gpd.shape #147 spatial entities
bk_2000_gpd.shape #2025 spatial entities


######## PRODUCE MAPS OF 2000 Population #########

#First need to link it.

bg_2000_gpd.columns # check columns' name for the data frame
census_syr_df.columns #contains census variables
#key is "TRACT" but with a different format.
#First fix the format
bg_2000_gpd.head()
census_syr_df.BKG_KEY.head()

#as.numeric(as.character(ct_2000_sp$TRACT))
ct_2000_gpd.TRACT.dtype
bg_2000_gpd.BKG_KEY.dtypes
census_syr_df.dtypes
census_syr_df.BKG_KEY.dtypes

#bg_2000_gpd['BKG_KEY'].astype(census_syr_df.BKG_KEY.dtypes)
bg_2000_gpd['BKG_KEY']=bg_2000_gpd['BKG_KEY'].astype('int64')

#ct_2000_sp$TRACT <- as.numeric(as.character(ct_2000_sp$TRACT))

#bg_2000_sp = merge(bg_2000_sp,census_syr_df,by="BKG_KEY")
bg_2000_gpd = bg_2000_gpd.merge(census_syr_df, on='BKG_KEY')
# country_shapes = country_shapes.merge(country_names, on='iso_a3')

#spplot(bg_2000_sp,"POP2000",main="POP2000") #quick visualization of population 
bg_2000_gpd.plot(column='POP2000',cmap="OrRd",
                 scheme='quantiles')
plt.title('POPULATION 2000')
#plt.title('POP2000')

### Let's use more option with matplotlib
fig, ax = plt.subplots()
bg_2000_gpd.plot(column='POP2000',cmap="OrRd",
                 scheme='quantiles',
                 ax=ax)
ax.set_title('POP2000')
##Now change the classes!

### Summarize by census track
#census_2000_sp <- aggregate(bg_2000_sp , 
#                            by="TRACT",FUN=sum)
##Note losing TRACT field
census_2000_df = bg_2000_gpd.groupby('TRACT',as_index=False).sum()
type(census_2000_df)
#To keep geometry, we must use dissolve method from geopanda
census_2000_gpd = bg_2000_gpd.dissolve(by='TRACT',
                                       aggfunc='sum')
type(census_2000_gpd)
census_2000_gpd.index
#note that the TRACT field has become the index
#census_2000_gpd['TRACT']=census_2000_gpd.index
#census_2000_gpd.set_index(list(np.arange(0,census_2000_df.shape[0])))
census_2000_gpd=census_2000_gpd.reset_index()

## Sum is 57!!!
sum(census_2000_gpd.POP2000==census_2000_df.POP2000)
census_2000_df.POP2000[0:2]
census_2000_gpd.POP2000[0:2]

##Index don't match
#census_2000_df.index['TRACT']
census_2000_gpd.shape[0] == census_2000_df.shape[0]
## Add back the tract ID Using sjoin??
#test3=gpd.tools.sjoin(census_2000_gpd,ct_2000_gpd,
#                      how='left')

##compare to sp!!
#df_test <- aggregate(POP2000 ~ TRACT, bg_2000_sp , 
#                     FUN=sum)

### Check if the new geometry of entities is the same as census

fig, ax = plt.subplots(figsize=(12,8))

# set aspect to equal. This is done automatically
# when using *geopandas* plot on it's own, but not when
# working with pyplot directly.
ax.set_aspect('equal')
census_2000_gpd.plot(ax=ax,column='POP2000',cmap='OrRd')
ct_2000_gpd.plot(ax=ax,color='white',edgecolor="red",alpha=0.7)

ax.set_title("Population", fontsize= 20)
#fig.colorbar(ax) #add palette later
#ax.set_axis_off()
#plt.show()

#plot(census_2000_sp)
#plot(ct_2000_sp,border="red",add=T)
#nrow(census_2000_sp)==nrow(ct_2000_sp)

#df_summary_by_census <- aggregate(. ~ TRACT, 
#                                  bg_2000_sp 
#                                  ,FUN=sum) #aggregate all variables from the data.frame

##Join by key table id:
#dim(ct_2000_sp)
#ct_2000_sp <- merge(ct_2000_sp,df_summary_by_census,by="TRACT")
#dim(ct_2000_sp)
metals_df.dtypes
ct_2000_gpd.dtypes
ct_2000_gpd['TRACT']=ct_2000_gpd.TRACT.astype('int64')

ct_2000_gpd.shape
ct_2000_gpd = ct_2000_gpd.merge(census_2000_df, on='TRACT')
ct_2000_gpd.shape
#save as sp and text table
#write.table(file.path(out_dir,)

### Do another map with different class intervals: 

title_str = "Population by census tract in 2000"
#range(ct_2000_sp$POP2000)
#col_palette <- matlab.like(256)

## 9 classes with fixed and constant break
#break_seq <- seq(0,9000,1000)
#breaks.qt <- classIntervals(ct_2000_sp$POP2000, n=length(break_seq), 
#                            style="fixed", fixedBreaks=break_seq, intervalClosure='right')

## Color for each class
#colcode = findColours(breaks.qt , c('darkblue', 'blue', 'lightblue', 'palegreen','yellow','lightpink', 'pink','brown3',"red","darkred"))
#p_plot_pop2000_ct <- spplot(ct_2000_sp,
#                            "POP2000",
#                            col="transparent", #transprent color boundaries for polygons
#                            col.regions = col_palette ,
#                            main=title_str,
#                            at = breaks.qt$brks)
#print(p_plot_pop2000_ct)

##### PART II: SPATIAL QUERY #############

## Join metals to census track 
## Join lead (pb) measurements to census tracks

#soil_PB_df <- read.table(file.path(in_dir_var,census_table_fname),sep=",",header=T)
#metals_df= pd.read_xls(os.path.join(in_dir,metals_table_fname))
#metals_df <- read.xls(file.path(in_dir_var,metals_table_fname),sep=",",header=T)

#View(soil_PB_df)
metals_df.head()

##This suggests matching to the following spatial entities
#nrow(metals_df)==nrow(ct_2000_sp)
metals_df.shape[0]== ct_2000_gpd.shape[0]
#nrow(soil_PB_df)==nrow(bg_2000_sp)

#dim(bg_2000_sp)
#census_metals_sp <- merge(ct_2000_sp,metals_df,by.x="TRACT",by.y="ID")
#Check data types before joining tables with "merge"
#metals_df.dtypes
#ct_2000_gpd.dtypes
#ct_2000_gpd['TRACT']=ct_2000_gpd.TRACT.astype('int64')
census_metals_gpd = ct_2000_gpd.merge(metals_df,left_on='TRACT',right_on='ID')

########processing lead data
### Now let's plot lead data 
#Soil lead samples: UTM z18 coordinates
#soil_PB_df <- read.table(file.path(in_dir_var,soil_PB_table_fname),sep=",",header=T) #point locations
census_metals_gpd.crs
#proj4string(census_metals_sp) #
soil_PB_df.columns = ["x","y","ID","ppm"]
#names(soil_PB_df)
soil_PB_df.head()
#names(soil_PB_df) <- c("x","y","ID","ppm") 

soil_PB_gpd = soil_PB_df.copy()
type(soil_PB_df)
soil_PB_gpd['Coordinates']=list(zip(soil_PB_gpd.x,soil_PB_gpd.y))
#coordinates(soil_PB_sp) <- soil_PB_sp[,c("x","y")]
#coordinates(soil_PB_sp) <- soil_PB_sp[,c("x","y")]
type(soil_PB_gpd)
soil_PB_gpd['Coordinates']= soil_PB_gpd.Coordinates.apply(Point)
soil_PB_gpd = gpd.GeoDataFrame(soil_PB_gpd,geometry='Coordinates')

#### Check the coordinates reference system
type(census_metals_gpd.crs) #dictionary
epsg_code = census_metals_gpd.crs.get('init').split(':')[1]

inproj = osr.SpatialReference()
inproj.ImportFromEPSG(int(epsg_code))
inproj.ExportToProj4()

#Assign projection system
soil_PB_gpd.crs= census_metals_gpd.crs
#proj4string(soil_PB_sp) <- proj4string(census_metals_sp)
#dim(soil_PB_sp)

#
#soil_PB_sp <- soil_PB_sp[,c("ID","ppm","x","y")]
#View(soil_PB_sp)
soil_PB_gpd.head()

fig, ax = plt.subplots()

census_metals_gpd.plot(ax=ax,color='white',edgecolor='red')
soil_PB_gpd.plot(ax=ax,marker='*',
                 color='black',
                 markersize=0.8)
                 
#plot(census_metals_sp)
#plot(soil_PB_sp,add=T)

###### Spatial query: associate points of pb measurements to each census tract
### Get the ID and 
#Use sjoin
##### Problem here ********************
test=gpd.tools.sjoin(soil_PB_gpd,census_2000_gpd,
                     how="left")
#test2=gpd.tools.sjoin(test,ct_2000_gpd,"left")
len(test.BKG_KEY.value_counts()) #associated BKG Key to points
len(test.index_right.value_counts())
test.columns
grouped = test.groupby(['index_right']).mean()
grouped = grouped.reset_index()
#soil_tract_id_df <- over(soil_PB_sp,census_2000_sp,fn=mean)
#soil_PB_sp <- intersect(soil_PB_sp,census_2000_sp)
#test4 <- gIntersection(soil_PB_sp,census_2000_sp,byid=T)
#names(soil_PB_sp)
#grouped = grouped.rename(columns={'index_right': 'TRACT',
#                            'ppm': 'pb_ppm' })
grouped = grouped.rename(columns={'ppm': 'pb_ppm' })

#soil_PB_sp <- rename(soil_PB_sp, c("d"="TRACT")) #from package plyr

#census_pb_avg <- aggregate(ppm ~ TRACT,(soil_PB_sp),
#                           FUN=mean)
#census_pb_avg <- rename(census_pb_avg,c("ppm"="pb_ppm"))

##Now join
#census_metals_pb_sp <- merge(census_metals_sp,
#                             census_pb_avg,by="TRACT")

census_metals_gpd = census_metals_gpd.merge(grouped,on="TRACT")

### write out final table and shapefile

#outfile<-paste("census_metals_pb_sp","_",
#               out_suffix,sep="")
#writeOGR(census_metals_pb_sp,dsn= out_dir,layer= outfile, driver="ESRI Shapefile",overwrite_layer=TRUE)

outfile = "census_metals_pb_"+'_'+out_suffix+'.shp'
census_metals_gpd.to_file(outfile)

#outfile_df_name <- file.path(out_dir,paste0(outfile,".txt"))
#write.table(as.data.frame(census_metals_pb_sp),file=outfile_df_name,sep=",")

## For kriging use scipy.interpolate
#https://stackoverflow.com/questions/45175201/how-can-i-interpolate-station-data-with-kriging-in-python


##### PART IV: Vulnerability to metals #############
#Examine the relationship between metals, Pb and vulnerable populations in Syracuse

#P2- SPATIAL AND NON SPATIAL QUERIES (cannot use spatial join)
#GOAL: Answer a set of questions using spatial and attribute queries and their combinations

#Produce:
#  a) two different maps based on two different definitions that answer the question:  which areas have high levels of children and are predominantly minority AND are at risk of heavy metal exposure using at least three variables. Use only tabular operations
#b) Same question as a) but using both spatial and tabular operations

#Note: In both cases include the method, variables used and your definition of risk areas in each 4 maps. The definition of risk is your own, you can also follow an established standard that would make sense or is official.  
#From these products, the layman should be able to answer the following questions:
#  a. Where are the areas of high heavy metal exposure that also have high levels of children population that belong to a demographic minority(s)? 
#b. Is there a different outcome in using tabular methods only vs combining tabular and spatial query methods?

#lm(,data=census_metals_pb_sp)
#moran(x, listw, n, S0, zero.policy=NULL, NAOK=FALSE)

census_metals_gpd.index
#census_metals_gpd.reset_index(drop=True)
census_metals_gpd = census_metals_gpd.set_index('TRACT')

w_queen = ps.weights.queen_from_shapefile(outfile,idVariable='TRACT')

#q_weights = ps.weights.queen_from_shapefile(census_metals_gpd,idVariable='TRACT')
w_queen.transform = 'R'
w_queen.neighbors
w_queen.n # number of observations (spatia features)
w_queen.mean_neighbors

y = census_metals_gpd.pb_ppm_x

m_I = ps.Moran(y,w_queen)

m_I.I
m_I.EI

y_lag = ps.lag_spatial(w_queen,y) #this is a numpy array
census_metals_gpd['y'] = census_metals_gpd.pb_ppm_x
census_metals_gpd['y_lag'] = y_lag

sns.regplot(x=y,y=y_lag,data=census_metals_gpd)

#list_nb <- poly2nb(census_lead_sp) #generate neighbours based on polygons
#summary(list_nb)
#plot(census_lead_sp,border="blue")
#plot.nb(list_nb,coordinates(census_lead_sp),add=T)

#generate weights
#nb2listw
#list_w <- nb2listw(list_nb, glist=NULL, style="W", zero.policy=NULL) #use row standardized
#can.be.simmed(list_w)
#summary(list_w)
#plot(as.matrix(list_w))
#moran(census_lead_sp$pb_ppm,list_w,n=nrow(census_lead_sp), Szero(list_w))
#moran.plot(census_lead_sp$pb_ppm, list_w,
#           labels=as.character(census_lead_sp$TRACT), pch=19)

##### Now do a spatial regression
#Use numpy array so convert with values
y.values.shape #not the right dimension
y = y.values.reshape(len(y),1)
y_lag = y_lag.reshape(len(y_lag),1)

x = census_metals_gpd['perc_hispa']
x = x.values.reshape(len(x),1)

mod_ols = ps.spreg.OLS(y,x)
mod_ols.u 
m_I_residuals = ps.Moran(mod_ols.u,w_queen)
#take into account autocorr in spreg
mod_ols.summary
mod_ols_test = ps.spreg.OLS(y,x,w_queen)
mod_ols_test.summary

mod_ml_lag = ps.spreg.ML_Lag(y,x,w_queen)
mod_ml_lag.summary
#replace explicative variable later! 

#mod_lm <- lm(pb_ppm ~ perc_hispa, data=census_lead_sp)
#mod_lag <- lagsarlm(pb_ppm ~ perc_hispa, data=census_lead_sp, list_w, tol.solve=1.0e-30)

#moran.test(mod_lm$residuals,list_w)
#moran.test(mod_lag$residuals,list_w)

#### Compare Moran's I from raster to Moran's I from polygon sp
# Rook's case
f <- matrix(c(0,1,0,1,0,1,0,1,0), nrow=3)
Moran(r_lead, f)

#http://rspatial.org/analysis/rst/7-spregression.html

###################### END OF SCRIPT #####################

















