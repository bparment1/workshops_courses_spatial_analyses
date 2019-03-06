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
#DATE MODIFIED: 03/06/2019
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
import matplotlib.cm as cm
import matplotlib.colors as colors
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
from collections import OrderedDict
import webcolors

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

############################################################################
#####  Parameters and argument set up ########### 

#ARGS 1
in_dir = "/home/bparmentier/c_drive/Users/bparmentier/Data/python/Exercise_4/data"
#ARGS 2
out_dir = "/home/bparmentier/c_drive/Users/bparmentier/Data/python/Exercise_4/outputs"
#ARGS 3:
create_out_dir=True #create a new ouput dir if TRUE
#ARGS 7
out_suffix = "exercise4_03062018" #output suffix for the files and ouptut folder
#ARGS 8
NA_value = -9999 # number of cores
file_format = ".tif"

#NLCD coordinate reference system: we will use this projection rather than TX.
CRS_reg = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
method_proj_val = "bilinear" # method option for the reprojection and resampling
gdal_installed = True #if TRUE, GDAL is used to generate distance files
		
### Input data files
rastername_county_harris = "harris_county_mask.tif" #Region of interest: extent of Harris County
elevation_fname = "srtm_Houston_area_90m.tif" #SRTM elevation
roads_fname = "r_roads_Harris.tif" #Road count for Harris county
	
### Aggreagate NLCD input files
infile_land_cover_date1 = "agg_3_r_nlcd2001_Houston.tif"
infile_land_cover_date2 = "agg_3_r_nlcd2006_Houston.tif"
infile_land_cover_date3 = "agg_3_r_nlcd2011_Houston.tif"
	
infile_name_nlcd_legend = "nlcd_legend.txt"
infile_name_nlcd_classification_system = "classification_system_nlcd_legend.xlsx"
	
data_fname = 'r_variables_harris_county_exercise4_02072019.txt'

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
    
    
###########################################
### PART I: READ AND VISUALIZE DATA #######
	
infile_land_cover_date1 = os.path.join(in_dir,infile_land_cover_date1) #NLCD 2001
infile_land_cover_date2 = os.path.join(in_dir,infile_land_cover_date2) #NLCD 2006
infile_land_cover_date3 = os.path.join(in_dir,infile_land_cover_date3) #NLCD 2011

lc_date1 = rasterio.open(infile_land_cover_date1) 
r_lc_date1 = lc_date1.read(1,masked=True) #read first array with masked value, nan are assigned for NA
lc_date2 = rasterio.open(infile_land_cover_date2) 
r_lc_date2 = lc_date2.read(1,masked=True) #read first array with masked value, nan are assigned for NA
lc_date3= rasterio.open(infile_land_cover_date2) 
r_lc_date3 = lc_date3.read(1,masked=True) #read first array with masked value, nan are assigned for NA

spatial_extent = rasterio.plot.plotting_extent(lc_date1)
plot.show(r_lc_date1)

#Note that you can also plot the raster io data reader
type(lc_date2)
plot.show(lc_date2) #no need to read in memory

lc_date1.crs # not defined with *.rst
lc_legend_df = pd.read_table(os.path.join(in_dir,infile_name_nlcd_legend),sep=",")
	
lc_legend_df.head() # Inspect data
plot.show(lc_date2) # View NLCD 2006, we will need to add the legend use the appropriate palette!!
plot.show(lc_date2,cmap=plt.cm.get_cmap('cubehelix',16 ))	
### Let's generate a palette from the NLCD legend information to view the existing land cover for 2006.
#names(lc_legend_df)
lc_legend_df.columns
lc_legend_df.shape
#subset the data to remove unsured rows
lc_legend_df = lc_legend_df[lc_legend_df['COUNT']>0] 


# Generate palette

colors_val = ['linen', 'lightgreen', 'darkgreen', 'maroon']

cmap = colors.ListedColormap(colors_val) # can be used directly
webcolors.rgb_to_name
webcolors.rgb_to_name((0, 0, 0)) #default is css3 convention
webcolors.rgb_to_name((255,255,255))
webcolors.name_to_rgb('navy')

### Generate a palette color from the input Red, Green and Blue information using RGB encoding:
rgb_col=list(zip(lc_legend_df.Red,lc_legend_df.Green,lc_legend_df.Blue))
len(rgb_col)
rgb_col[0]
#lc_legend_df$rgb <- paste(lc_legend_df$Red,lc_legend_df$Green,lc_legend_df$Blue,sep=",") #combine
#','.join([lc_legend_df.Red,lc_legend_df.Green, lc_legend_df.Blue]) 	

#lc_legend_df['rgb'] = lc_legend_df[['Red','Green','Blue']].apply[lambda x:]
lc_legend_df['rgb'] = rgb_col
### row 2 correspond to the "open water" category
#webcolors.rgb_to_name(rgb_col[1])

color_val_water = rgb_col[1]
color_val_developed_high = rgb_col[7]

#color_val_water <- rgb(lc_legend_df$Red[2],lc_legend_df$Green[2],lc_legend_df$Blue[2],maxColorValue = 255)
#color_val_developed_high <- rgb(lc_legend_df$Red[7],lc_legend_df$Green[7],lc_legend_df$Blue[7],maxColorValue = 255)
	
#lc_col_palette <- c(color_val_water,color_val_developed_high)
	
#barplot(c(1,1),
#col=lc_col_palette,
#main="Visualization of color palette for NLCD land cover",
#names.arg=c("Open water", "Developed, High Intensity"),las=1)
	
### Let's generate a color for all the land cover categories by using lapply and function
#n_cat <- nrow(lc_legend_df)
#lc_col_palette <- lapply(1:n_cat,
#FUN=function(i){rgb(lc_legend_df$Red[i],lc_legend_df$Green[i],lc_legend_df$Blue[i],maxColorValue = 255)})
#lc_col_palette <- unlist(lc_col_palette)
	
	
### Now generate a plot of land cover with the NLCD legend and palette
#levelplot(r_lc_date2,
#	col.regions = lc_legend_df_date2$palette,
#	scales=list(draw=FALSE),
#	main = "NLCD 2006")
	

#######
################################################
###  PART II : Analyze change and transitions

## As the plot shows for 2006, we have 15 land cover types. Analyzing such complex categories in terms of decreasse (loss), increase (gain), 
# persistence in land cover will generate a large number of transitions (potential up to 15*15=225 transitions in this case!)

## To generalize the information, let's aggregate leveraging the hierachical nature of NLCD Anderson Classification system.

lc_system_nlcd_df = pd.read_excel(os.path.join(in_dir,infile_name_nlcd_classification_system))
lc_system_nlcd_df.head #inspect data

val, cnts =np.unique(r_lc_date1,return_counts=True)

df = pd.DataFrame(np.ma.filled(val))
pd.set_option('display.float_format', lambda x: '%.3f' % x)
df_date1 = pd.DataFrame(val,cnts)
df_date1 = df_date1.reset_index()
df_date1.columns = ['y_2001','value']

val, cnts =np.unique(r_lc_date2,return_counts=True)
df = pd.DataFrame(np.ma.filled(val))
pd.set_option('display.float_format', lambda x: '%.3f' % x)
df_date2 = pd.DataFrame(val,cnts)
df_date2 = df_date2.reset_index()
df_date2.columns = ['y_2006','value']

### Let's identify existing cover and compute change:
freq_tb_nlcd = pd.merge(df_date1,df_date2,on='value')
#reorder columns 
freq_tb_nlcd = freq_tb_nlcd[['value','y_2001','y_2006']]
freq_tb_nlcd.head()
#dim(lc_system_nlcd_df) # We have categories that are not relevant to the study area and time period.
#lc_system_nlcd_df <- subset(lc_system_nlcd_df,id_l2%in%freq_tb_nlcd$value ) 
#dim(lc_system_nlcd_df) # Now 15 land categories instead of 20.

lc_system_nlcd_df.shape
selected_cat = lc_system_nlcd_df.id_l2.isin(freq_tb_nlcd.value)
lc_system_nlcd_df = lc_system_nlcd_df[selected_cat]

### Select relevant columns for the reclassification
rec_df <- lc_system_nlcd_df[,c(2,1)]
r_date1_rec <- subs(r_lc_date1,rec_df,by="id_l2","id_l1")
r_date2_rec <- subs(r_lc_date2,rec_df,by="id_l2","id_l1")

plot(r_date1_rec)

rec_xtab_df <- crosstab(r_date1_rec,r_date2_rec,long=T)
names(rec_xtab_df) <- c("2001","2011","freq")

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

###########################################
############# PART III: Process and prepare for land change modeling ####################
## add this later


###########################################
### PART IV: Run model and perform assessment ###########################

################
##### Step 1: consistent masking and generate mask removing water (1) and developped (2) in 2001
	
infile_land_cover_date1 = os.path.join(in_dir,infile_land_cover_date1) #NLCD 2001

#data_df = pd.read_table(os.path.join(in_dir,data_fname))
data_df = pd.read_csv(os.path.join(in_dir,data_fname))
data_df.columns

################
##### Step 2: Prepare model for predictions: Split data into training and testing and rescaling

from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
from numpy import array

## Relevant variables used:
selected_covariates_names = ['land_cover', 'slope', 'roads_dist', 'developped_dist']
selected_target_names = ['change'] #also called dependent variable

selected_categorical_var_names=['land_cover']

selected_continuous_var_names=list(set(selected_covariates_names) - set(selected_categorical_var_names))
##Find frequency of unique values:
freq_val_df = data_df[selected_categorical_var_names].apply(pd.value_counts)
print(freq_val_df.head())

values_cat = array(data_df[selected_categorical_var_names].values) #note this is assuming only one cat val here

label_encoder = LabelEncoder() 
one_hot_encoder = OneHotEncoder(sparse=False)

### First integer encode:
integer_encoded = label_encoder.fit_transform(values_cat)
print(integer_encoded)

# Binary encode:

integer_encoded = integer_encoded.reshape(len(integer_encoded),1)
print(integer_encoded)

onehot_encoded = one_hot_encoder.fit_transform(integer_encoded)
print(onehot_encoded)
onehot_encoded.shape
type(onehot_encoded)

#invert to check value?
onehot_encoded[0:5,]
values_cat[0:5,]

inverted = label_encoder.inverse_transform([np.argmax(onehot_encoded[0,:])])
inverted = label_encoder.inverse_transform([np.argmax(onehot_encoded[1,:])])
print(inverted)

#assign back to the data.frame

unique_val = np.array(freq_val_df.index)
unique_val = np.sort(unique_val)

print(unique_val)
names_cat = ['lc_' + str(i) for i in unique_val]

print(names_cat)
onehot_encoded_df = pd.DataFrame(onehot_encoded,columns=names_cat)
onehot_encoded_df.columns
onehot_encoded_df.head()
onehot_encoded_df.shape
data_df.shape
## Combine back!!

data_df= pd.concat([data_df,onehot_encoded_df],sort=False,axis=1)
data_df.shape
data_df.head()

selected_covariates_names_updated = selected_continuous_var_names + names_cat 

## Split training and testing

from sklearn.model_selection import train_test_split

prop = 0.3
random_seed = 100

X_train, X_test, y_train, y_test = train_test_split(data_df[selected_covariates_names_updated], 
                                                    data_df[selected_target_names], 
                                                    test_size=prop, 
                                                    random_state=random_seed)

X_train.shape

#### Scaling between 0-1 for continuous variables

from sklearn.preprocessing import MinMaxScaler

# Data needs to be scaled to a small range like 0 to 1 for the neural
# network to work well.
scaler = MinMaxScaler(feature_range=(0, 1))

##### need to use one hot encoding or text embedding to normalize categorical variables
### need to select only the continuous var:
scaled_training = scaler.fit_transform(X_train[selected_continuous_var_names])
scaled_testing = scaler.transform(X_test[selected_continuous_var_names])

type(scaled_training) # array
scaled_training.shape

## Concatenate column-wise
X_testing_df = pd.DataFrame(np.concatenate((X_test[names_cat].values,scaled_testing),axis=1),
                                            columns=names_cat+selected_continuous_var_names)

X_training_df = pd.DataFrame(np.concatenate((X_train[names_cat].values,scaled_training),axis=1),
                                            columns=names_cat+selected_continuous_var_names)

#scaled_training_df.to_csv("sales_data_training_scaled.csv", index=False)
#scaled_testing_df.to_csv("sales_data_testing_scaled.csv", index=False)

####################
###### Step 3: fit glm logistic model and generate predictions

##################################
### logistic model
from sklearn.linear_model import LogisticRegression

model_logistic = LogisticRegression()

model_logistic = model_logistic.fit(X_train.values,y_train.values.ravel())

model_logistic.coef_
selected_covariates_names_updated

pred_test = model_logistic.predict(X_test.values)
pred_test_prob = model_logistic.predict_proba(X_test.values)

pred_test_prob[:,1] # this is the prob for 1
y_test[0:5]
pred_test_prob[0:5,:]

predicted_classes = model.predict(X)
accuracy = accuracy_score(y.flatten(),pred_test)
parameters = model.coef_
#pred_test = model_logistic.predict(X_test)

model_logistic.score(pred_test,y_test)

pred_test_prob = model_logistic.predict_proba(X_test.values)
y_scores_test = pred_test_prob[:,1]

pred_train_prob = model_logistic.predict_proba(X_train.values)
y_scores_train = pred_test_prob[:,1]

### Note that we only have about 10% change in the dataset so setting 50% does not make sense!!
sum(data_df.change)/data_df.shape[0]
sum(y_train.change)/y_train.shape[0]

data_df['change'].value_counts()

sns.set(color_codes=True) #improves layout with bar and background grid
sns.countplot(x='change',data=data_df)
#,palette='hls')
plt.show()
plt.savefig('count_plot')

plt.figure()
sns.distplot(y_scores_test) #histogram for the predicted probablities

sns.distplot(y_scores_train) #histogram for the predicted probablities

####################
###### Step 5: Model assessment with ROC and AUC


#https://towardsdatascience.com/building-a-logistic-regression-in-python-301d27367c24
#This is for ROC curve
#https://towardsdatascience.com/building-a-logistic-regression-in-python-step-by-step-becd4d56c9c8

from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

y_true = y_test
y_scores = pred_test_prob[:,1]
auc_val_test =roc_auc_score(y_true,y_scores)
auc_val_test =roc_auc_score(y_true,y_scores)

#fpr, tpr, thresholds = roc_curve(y_test, 
#                                 logreg.predict_proba(X_test)[:,1])

fpr, tpr, thresholds = roc_curve(y_test, 
                                 y_scores)
plt.figure()
plt.plot(fpr, tpr, 
         label='Logistic Regression (area = %0.2f)' % auc_val)
plt.plot([0, 1], [0, 1],'r--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.savefig('Log_ROC')
plt.show()

###################### END OF SCRIPT #####################










                   






