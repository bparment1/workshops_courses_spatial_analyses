# -*- coding: utf-8 -*-
"""
Spyder Editor.
"""
#################################### Regression Temperature #######################################
######################## Analyze and predict air temperature with Earth Observation data #######################################
#This script performs analyses to predict air temperature using several coveriates.
#The goal is to predict air temperature using Remotely Sensing data as well as compare measurements
# from the ground station to the remotely sensed measurements.
#
#AUTHORS: Benoit Parmentier
#DATE CREATED: 09/07/2018
#DATE MODIFIED: 03/04/2019
#Version: 1
#PROJECT: SESYNC Geospatial Course and AAG 2019 Python Geospatial Course
#TO DO:
#
#COMMIT: clean up code for workshop
#
#################################################################################################


###### Libr  
# Read raster bands directly to Numpy arrays.
with rasterio.open(os.path.join(in_dir,infile)) as src:
        r_lst = src.read(1,masked=True) #read first array with masked value, nan are assigned for NA
        spatial_extent = rasterio.plot.plotting_extent(src)

#test = rasterio.open(os.path.join(in_dir,infile))
plot.show(r_lst)
#plot.show(r_lst,cmap='viridis',scheme='quantiles')

src.crs # explore Coordinate Reference System 
plot.show(src)
type(r_lst)
r_lst.size
src.shape
src.height

#Can also use the regular matplotlib library function to plot images
plt.imshow(r_lst)

#see: https://matplotlib.org/users/image_tutorial.html
plt.imshow(r_lst, clim=(259.0, 287.0))
plt.hist(r_lst.ravel(),bins=256,range=(259.0,287.0))

data_gpd.plot(marker="*",color="green",markersize=5)
station_or = data_gpd.to_crs({'init': 'epsg:2991'}) #reproject to  match the  raster image

#https://www.earthdatascience.org/courses/earth-analytics-python/lidar-raster-data/customize-matplotlib-raster-maps/

##### How to combine plots:
fig, ax = plt.subplots()
with rasterio.open(os.path.join(in_dir,infile)) as src:
        rasterio.plot.show((src,1),ax=ax,
                          clim=(259.0,287.0),)

#plot.show(r_lst, clim=(259.0, 287.0),ax=ax)
#with rasterio.plot.show((src,1),ax=ax)
station_or.plot(ax=ax,marker="*",
              color="red",
               markersize=10)

fig, ax = plt.subplots(figsize = (8,3))
lst_plot = ax.imshow(r_lst, 
                       cmap='Greys', 
                       extent=spatial_extent)
ax.set_title("Long term mean for January land surface temperature", fontsize= 20)
fig.colorbar(lst_plot)
# turn off the x and y axes for prettier plotting
#ax.set_axis_off(); #this removes coordinates on the plot

###########################################
### PART II : Extract information from raster and prepare covariates #######
#raster = './data/slope.tif'

data=gr.from_file(os.path.join(in_dir,infile))
type(data) # check that we have a georaster object
# Plot data
data.plot()
data.plot(clim=(259.0, 287.0))

#### Extract information from raster using coordinates
x_coord = station_or.geometry.x # pands.core.series.Series
y_coord = station_or.geometry.y
# Find value at point (x,y) or at vectors (X,Y)
values = data.map_pixel(x_coord,y_coord)
station_or.columns #get names of col

station_or['year'].value_counts()
station_or.groupby(['month'])['value'].mean()
     
print("number of rows:",station_or.station.count(),"number of stations:",len(station_or.station.unique()))
station_or['LST1'] = values - 273.15 #create new column

station_or_jan = station_or.loc[(station_or['month']==1) & (station_or['value']!=-9999)]
station_or_jan.head()
station_or_jan.columns

#avg_df = station_or.groupby(['station'])['value'].mean())
avg_df = station_or_jan.groupby(['station'])['value','LST1'].mean()
avg_df['value']= avg_df['value']/10
avg_df.head()
         
################################################
###  PART III : Fit model and generate prediction

### Add split training and testing!!!
### Add additionl covariates!!

#selected_covariates_names_updated = selected_continuous_var_names + names_cat 
selected_covariates_names_updated = ['LST1'] 
selected_target_names = ['value']
## Split training and testing

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(avg_df[selected_covariates_names_updated], 
                                                    avg_df[selected_target_names], 
                                                    test_size=prop, 
                                                    random_state=random_seed)

X_train.shape

from sklearn.linear_model import LinearRegression
regr = LinearRegression().fit(X_train,y_train)

#regr = linear_model.LinearRegression()
regr.fit(x, y)

plt.scatter(X_train, y_train,  color='black')
plt.plot(X_train, regr.predict(X_train), color='blue', linewidth=3)
plt.xticks(())
plt.yticks(())
plt.show()

print('reg coef',regr.coef_)
print('reg intercept',regr.intercept_)

#reg.predict(x) # Note this is a fit!
#reg.score(x, y)

regr.predict(X_train) # Note this is a fit!
regr.score(X_train, y_train)

X_test.shape
regr.predict(X_test) # Note this is a fit!
regr.score(X_test, y_test)

## As the plot shows for 2006, we have 15 land cover types. Analyzing such complex categories in terms of decreasse (loss), increase (gain),
### Do models for January,July with LST and with/without land cover % of forest
## Calculate MAE,RMSE,R2,etc. inspire yourself from paper. Save this into a CSV file.

############################# END OF SCRIPT ###################################

















