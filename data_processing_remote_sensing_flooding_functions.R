####################################    Flood Mapping Analyses   #######################################
############################  Analyze and map flooding from RITA hurricane  #######################################
#This script performs data processing for the Exercise 4, Exercise 5 and Exercise 6 for the Geospatial Short Course 
#This is using reflectance data derived from MODIS MOD09 and NLCD data.
#
#AUTHORS: Benoit Parmentier                                             
#DATE CREATED: 03/07/2018 
#DATE MODIFIED: 03/07/2018
#Version: 1
#PROJECT: SESYNC and AAG 2018 workshop/Short Course preparation
#TO DO:
#
#COMMIT: initial commit
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
library(gstat) #spatial interpolation and kriging methods
library(readxl) #functionalities to read in excel type data
library(psych) #pca/eigenvector decomposition functionalities
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

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

generate_dates_by_step <-function(start_date,end_date,step_date){
  #library(xts) declare out of this function
  #library(zoo)
  #library(lubridate)
  
  st <- as.Date(start_date,format="%Y.%m.%d")
  en <- as.Date(end_date,format="%Y.%m.%d")
  #year_list <-seq(format(st,"%Y"),format(en,"%Y")) #extract year
  year_list <- seq(as.numeric(strftime(st,"%Y")),as.numeric(strftime(en,"%Y"))) #extract year
  
  ll_list <- vector("list",length=length(year_list))
  for (i in 1:length(year_list)){
    if(i==1){
      first_date <-st
    }else{
      first_date<-paste(year_list[[i]],"-01","-01",sep="")
    }
    if(i==length(year_list)){
      last_date <-en
    }else{
      last_date<-paste(year_list[[i]],"-12","-31",sep="")
    }
    #ll <- seq.Date(st, en, by=step)
    ll <- seq.Date(as.Date(first_date), as.Date(last_date), by=step_date)
    ll_list[[i]]<-as.character(ll)
    #paste(yday(ll,)
  }
  
  #
  dates_modis <-as.Date(unlist((ll_list))) 
  
  dates_DOY_modis <- as.character(paste(year(dates_modis),sprintf("%03d", yday(dates_modis)),sep=""))
  dates_obj <- list(dates_modis,dates_DOY_modis)
  names(dates_obj) <- c("dates","doy")  
  return(dates_obj)
}

aggregate_raster_fun <- function(l_rast,cat_names,use_majority,agg_fact,agg_fun,file_format=".tif",rast_ref=NULL,num_cores=1,out_suffix=NULL, out_dir=NULL){
  #
  #Aggregate raster from raster input and reference file
  #INPUT arguments:
  #1) l_rast : set of input raster layers as list
  #2) cat_names: within the list, give names of raster that contain categorical variables  
  #3) use_majority: if TRUE, will use this rule for cat variables
  #4) agg_fact: factor to aggregate
  #5) agg_fun: default is mean
  #6) file_format: e.g. ".tif"
  #7) rast_ref: reference raster to match in resolution, if NULL then send a message
  #6) file_Format: raster format used e.g. .tif
  #5) num_cores: number of cores to use:
  #5) out_suffix: output suffix
  #7) out_dir: output directory
  #8) out_rast_name: output raster name if null it is created from the input file name
  #OUTPUT:
  # out_raster_name: name of the file containing the aggregated raster
  #
  # Authors: Benoit Parmentier
  # Created: 03/02/2017
  # Modified: 03/08/2018
  # To Do: 
  # - Add option to disaggregate
  #
  ################################
  
  
  #
  #Function to aggregate input raster stack
  #if use majority then the zonal layer is aggregated and then reclassfied based by the majority rule
  
  #debug(aggregate_raster)
  #lf_agg_test <- aggregate_raster(l_rast[[1]],
  #                   #r_in=raster(lf_layerized_bool[1]),
  #                   agg_fact=agg_fact,
  #                   reg_ref_rast=NULL,
  #                   #agg_fun="mean",
  #                   agg_fun=agg_fun,
  #                   out_suffix=NULL,
  #                   file_format=file_format,
  #                   out_dir=out_dir,
  #                   out_rast_name = NULL) 
  
  
  if(!is.null(l_rast)){
    lf_agg <- mclapply(l_rast,
                       FUN=aggregate_raster,
                       #r_in=raster(lf_layerized_bool[1]),
                       agg_fact=agg_fact,
                       reg_ref_rast=reg_ref_rast,
                       #agg_fun="mean",
                       agg_fun=agg_fun,
                       out_suffix=NULL,
                       file_format=file_format,
                       out_dir=out_dir,
                       out_rast_name = NULL,
                       mc.preschedule=FALSE,
                       mc.cores = num_cores) 
    
    l_rast_original <- l_rast
    l_rast <- unlist(lf_agg) 
    
  }
  
  ###Break out and get mean per class and do majority rule!
  
  if(use_majority==TRUE){
    
    #l_rast_original
    #r_r_srtm_Katrina_rec2_NDVI_Katrina_03162017.rst"
    #r <- raster(paste0("r_",zonal_colnames,"_",out_suffix,file_format,sep=""))
    
    selected_cat_layers <- names(l_rast)==cat_names
    #raster_name <- (paste0("r_",cat_names,"_",out_suffix,file_format,sep="")) #can be a list of names
    #out_suffix_str <- paste0("agg_zonal","_",out_suffix)
    out_suffix_str <- out_suffix #may change this to included "majority" in the name
    #debug(generate_soft_cat_aggregated_raster_fun)
    #l_rast[selected_cat_layers]
    
    ## Use loop because we already have a num_cores
    for(i in 1:length(selected_cat_layers)){
      
      #debug(generate_soft_cat_aggregated_raster_fun)
      raster_name_cat <- l_rast[selected_cat_layers][[i]]
      lf_agg_soft <- generate_soft_cat_aggregated_raster_fun(raster_name_cat,
                                                             reg_ref_rast=reg_ref_rast,
                                                             agg_fact,
                                                             agg_fun,
                                                             num_cores,
                                                             NA_flag_val=NA_flag_val,
                                                             file_format,
                                                             out_dir,
                                                             out_suffix_str)
      
      if(raster_name_cat=="character"){
        raster_name_cat <- raster(raster_name_cat)
      }
      
      reclass_val <- unique(raster_name_cat) #unique zonal values to reassign
      #reclass_val <- c(0,1,2) # value for the elevation reclassified
      
      #debug(reclass_in_majority)
      r_stack <- stack(lf_agg_soft)
      r_reclass_obj <- reclass_in_majority(r_stack= r_stack,
                                           threshold_val=NULL,
                                           max_aggregation = TRUE,
                                           reclass_val = reclass_val)
      
      plot(r_reclass_obj$r_rec)
      rast_zonal <- r_reclass_obj$r_rec
      #zonal_colnames
      raster_name <- paste0("agg_",agg_fact,"_","r_",cat_names[i],"_",out_suffix,file_format)
      
      writeRaster(rast_zonal,
                  filename=file.path(out_dir,raster_name),
                  overwrite=TRUE)  
      
      
    }
    
    
  }
  
  if(use_majority==FALSE){
    #sure SRTM and reclass based on threshold values?
    
  }
  
  #r_srtm_Katrina_rec2
  #-rw-rw-r-- 1 bparmentier bparmentier 1894 Apr  7 12:33 r_r_srtm_Katrina_rec2_NDVI_Katrina_04062017.tif
  #-rw-rw-r-- 1 bparmentier bparmentier 1016 Apr  7 12:34 agg_5_r_r_srtm_Katrina_rec2_NDVI_Katrina_04062017.tif
  
  ###
  zonal_colnames <- gsub(extension(raster_name),"",raster_name)
  ##
  
  ##########################
  #### prepare return object
  
  obj <- list(zonal_colnames,l_rast,l_rast_original)
  names(obj) <- c("zonal_colnames","l_rast","l_rast_original")
  
  return(obj)
}

reclass_in_majority <- function(r_stack,threshold_val=0.5,max_aggregation=FALSE,reclass_val){
  ##
  #This function reclassify a set of soft layers using the majority or maximum value rule.
  #When max_aggregation is TRUE, the max value rule is used in the aggregation.
  #
  #INPUTS
  #1) r_stack
  #2) threshold_val
  #3) max_aggregation
  #4) reclass_val
  #
  
  ## Reclass
  if(!is.null(threshold_val) & (max_aggregation==FALSE)){
    r_rec_threshold <- r_stack > threshold_val
    #use calc later to save directly the file
    #
    
    r_rec_val_s <- lapply(1:nlayers(r_rec_threshold),
                          function(i,r_stack){df_subs <- data.frame(id=c(0,1),v=c(0,reclass_val[i]));
                          x <- subs(subset(r_stack,i), df_subs)},r_stack=r_rec_threshold)
    r_rec_val_s <- stack(r_rec_val_s) #this contains pixel above 0.5 with re-assigned values
    r_rec <- calc(r_test,function(x){sum(x)})
    
    ### prepare return object
    reclass_obj <- list(r_rec,r_rec_val_s)
    names(reclass_obj) <- c("r_rec","r_rec_val_s")
    
  }
  
  if(max_aggregation==TRUE){
    #r_zonal_agg_soft <- stack(lf_agg_soft)
    #Find the max, in stack of pixels (can be used for maximum compositing)
    r_max_s <- calc(r_stack, function(x) max(x, na.rm = TRUE))
    #maxStack <- stackApply(r_zonal_agg_soft, indices=rep(1,nlayers(r_zonal_agg_soft)), fun = max, na.rm=TRUE)
    r_max_rec_s <- overlay(r_stack,r_max_s, fun=function(x,y){as.numeric(x==y)})
    r_ties <- sum(r_max_rec_s) #find out ties
    #this may be long
    #freq_r_rec_df <- freq(r_rec_max,merge=T)
    
    r_rec_val_s <- lapply(1:nlayers(r_max_rec_s),
                          function(i,r_stack){df_subs <- data.frame(id=c(0,1),v=c(0,reclass_val[i]));
                          x <- subs(subset(r_stack,i), df_subs)},r_stack=r_max_rec_s)
    r_rec_val_s <- stack(r_rec_val_s)
    r_rec <- calc(r_rec_val_s,function(x){sum(x)})#overlays the layer with sum, 
    #x2 <- subs(r, df, subsWithNA=FALSE)
    
    ### prepare return object
    reclass_obj <- list(r_rec,r_rec_val_s,r_max_rec_s,r_ties)
    names(reclass_obj) <- c("r_rec","r_rec_val_s","r_max_rec_s","r_ties")
  }
  
  ###
  return(reclass_obj)
}

generate_soft_cat_aggregated_raster_fun <- function(r,reg_ref_rast,agg_fact,agg_fun,num_cores,NA_flag_val,file_format,out_dir,out_suffix){
  ## Function to aggregate categories
  
  ##INPUTS
  #1) r: raster to aggregate and crop/project if necessary
  #2) reg_ref_rast: reference raster, it must have a coordinate system defined
  #3) agg_fact:
  #4) agg_fun:
  #5) num_cores: number of cores used in the parallel processsing
  #6) NA_flag_val: flag value used for no data
  #7) file_format: raster file format e.g. .tif, .rst
  #8) out_dir: output directory
  #9) out_suffix: output suffix added to files
  #
  #OUTPUTS
  #
  #
  #
  
  
  ##### STEP 1: Check input
  
  if(class(r)!="RasterLayer"){
    r <- raster(r)
  }


  NAvalue(r) <- NA_flag_val #make sure we have a flag value assigned
  
  ###### STEP 1: BREAK OUT
  ## Breakout layers: separate categories in individual layers
  
  freq_tb <- as.data.frame(freq(r)) #zero is NA?
  out_filename_freq <- paste("freq_tb_",out_suffix,".txt",sep="")
  write.table(freq_tb,file= file.path(out_dir,out_filename_freq))
  
  ## get the names
  names_val <- freq_tb$value
  names_val <- names_val[!is.na(names_val)] #remove NA
  
  ## Make a brick composed of multiple layers: one layer per category (in one unique multiband file)
  out_raster_name <- file.path(out_dir,paste("r_layerized_bool_",out_suffix,file_format,sep=""))
  r_layerized <- layerize(r,
                          classes=names_val,
                          filename= out_raster_name,
                          overwrite=T)
  list_out_raster_name <- file.path(out_dir,paste("r_layerized_bool_",names_val,"_",out_suffix,file_format,sep=""))
  
  ## Now write out separate layers in one unique file (one per categories)
  #notice some issue here, looping through
  #for(i in 1:nlayers){
  #  writeRaster(r_layerized,
  #              bylayer=T,
  #              suffix=paste0(names_val,"_",out_suffix),
  #              filename=paste("r_layerized_bool",
  #                             file_format,sep="")
  #             ,overwrite=T)
  #}
  
  writeRaster(r_layerized,
              bylayer=T,
              suffix=paste0(names_val,"_",out_suffix),
              filename=paste("r_layerized_bool",
                             file_format,sep="")
              ,overwrite=T)
  
  #browser()
  
  lf_layerized_bool <- paste("r_layerized_bool_",names_val,"_",out_suffix,file_format,sep="")
  #names(r_layerized) <- 
  #inMemory(r_date_layerized)
  #filename(r_date_layerized)
  #r_test <- raster(lf_layerized_bool[1])
  
  ### STEP 2: aggregate
  ## Aggregate
  
  ##To do
  ##Set correct names based on categories of input
  ##Set output name
  
  #r_agg_fname <-aggregate_raster(agg_fact=agg_fact,
  #                               r_in=r_date_layerized,reg_ref_rast=NULL,agg_fun="mean",out_suffix=NULL,file_format=".tif",out_dir=NULL)
  #debug(aggregate_raster)
  #r_agg_fname <-aggregate_raster(r_in=raster(lf_layerized_bool[1]),
  #                               agg_fact=agg_fact,
  #                               reg_ref_rast=NULL,
  #                               #agg_fun="mean",
  #                               agg_fun=agg_fun,
  #                               out_suffix=NULL,
  #                               file_format=file_format,
  #                               out_dir=out_dir,
  #                               out_rast_name = NULL)
  #r_var_s <- mclapply(1:length(infile_var),FUN=import_list_modis_layers_fun,list_param=list_param_import_modis,mc.preschedule=FALSE,mc.cores = num_cores) #This is the end bracket from mclapply(...) statement
  
  #debug(aggregate_raster)
  #lf_agg <- aggregate_raster(lf_layerized_bool[1],
  #                           #r_in=raster(lf_layerized_bool[1]),
  #                           agg_fact=agg_fact,
  #                           reg_ref_rast=reg_ref_rast,
  #                           #agg_fun="mean",
  #                           agg_fun=agg_fun,
  #                           out_suffix=NULL,
  #                           file_format=file_format,
  #                           out_dir=out_dir,
  #                           out_rast_name = NULL)

  lf_agg <- mclapply(lf_layerized_bool,
                     FUN=aggregate_raster,
                     #r_in=raster(lf_layerized_bool[1]),
                     agg_fact=agg_fact,
                     reg_ref_rast=reg_ref_rast,
                     #agg_fun="mean",
                     agg_fun=agg_fun,
                     out_suffix=NULL,
                     file_format=file_format,
                     out_dir=out_dir,
                     out_rast_name = NULL,
                     mc.preschedule=FALSE,
                     mc.cores = num_cores) 
  #r_agg <- stack(unlist(lf_agg))
  raster_outname <- unlist(lf_agg)
  #apply function by layer using lapply.
  
  ## Reclassify by labels
  return(raster_outname)
}



#Function to aggregate from fine to coarse resolution, this will change accordingly once the input raster ref is given..
#
aggregate_raster <- function(r_in, agg_fact, reg_ref_rast=NULL,agg_fun="mean",out_suffix=NULL,file_format=".tif",out_dir=NULL,out_rast_name=NULL){
  #Aggregate raster from raster input and reference file
  #INPUT arguments:
  #1) r_in: input raster layer
  #2) agg_fact: factor to aggregate
  #3) reg_ref_rast: reference raster to match in resolution, if NULL then send a message
  #4) agg_fun: default is mean
  #5) out_suffix: output suffix
  #6) file_Format: raster format used e.g. .tif
  #7) out_dir: output directory
  #8) out_rast_name: output raster name if null it is created from the input file name
  #OUTPUT:
  # out_raster_name: name of the file containing the aggregated raster
  #
  # Authors: Benoit Parmentier
  # Created: 10/15/2015
  # Modified: 03/08/2017
  # To Do: 
  # - Add option to disaggregate
  #
  ################################
  
  ##If file is provided rather than RasterLayer
  if(class(r_in)!="RasterLayer"){
    r_in <- raster(r_in)
  }

  ##If file is provided rather than RasterLayer
  if(class(reg_ref_rast)!="RasterLayer"){
    reg_ref_rast <- raster(reg_ref_rast)
  }
  
  if(is.null(agg_fact)){
    res_ref <- res(reg_ref_rast)[1] #assumes square cells, and decimal degrees from WGS84 for now...
    res_in <- res(r_in)[1] #input resolution, assumes decimal degrees
    agg_fact <-round(res_ref/res_in) #find the factor needed..
    #fix this to add other otpions e.g. aggregating down
  }
  
  #Default values...
  if(is.null(out_suffix)){
    out_suffix <- ""
  }
  
  if(is.null(out_dir)){
    out_dir <- "."
  }
  
  ## Create output raster name if out_rast_name is null
  if(is.null(out_rast_name)){
    raster_name <- filename(r_in)
    extension_str <- extension(raster_name)
    raster_name_tmp <- gsub(extension_str,"",basename(raster_name))
    out_rast_name <- file.path(out_dir,paste("agg_",agg_fact,"_",raster_name_tmp,out_suffix,file_format,sep="")) #output name for aster file
    #out_rast_name <- raster_name #for use in function later...
  }
  
  r_agg <- aggregate(r_in, 
                     fact=agg_fact,
                     FUN=agg_fun,
                     filename=out_rast_name,
                     overwrite=TRUE)
  
  return(out_rast_name)
  
}


################################# End of Script #######################################
