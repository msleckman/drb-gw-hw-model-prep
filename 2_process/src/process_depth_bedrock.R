library(sf)
library(dplyr)
library(terra)
library(targets)
library(raster)
library(mapview)

targets::tar_load(p1_nhd_reaches)
targets::tar_load(p1_nhd_reaches_along_NHM)

## target #1 
## buffer reaches
buffer_in_m <- 100
buff_reaches <- st_buffer(p1_nhd_reaches, dist = buffer_in_m)
buff_reaches_along_nhm <- st_buffer(p1_nhd_reaches_along_NHM, dist = buffer_in_m)


## target #2
## catchments
start <- Sys.time()
nhd_catchments_along_nhm <- nhdplusTools::get_nhdplus(comid = p1_nhd_reaches_along_NHM$comid,
                                                      realization = 'catchment') 
end <- Sys.time()
print(end - start)

start <- Sys.time()
nhd_catchments <- nhdplusTools::get_nhdplus(comid = p1_nhd_reaches$comid, realization = 'catchment') 
end <- Sys.time()
print(end - start)

## function 
process_depth_to_bedrock <- function(adf_file_path,
                                    reach_catchments,
                                    comid_col = 'comid'){

  ## Vars -- to be removed
  #adf_file_path <- '1_fetch/in/Shangguan_dtb_cm_250m_clip 2/w001001.adf'
  # reach_catchments <- nhd_catchments_along_nhm %>% head()
  #reach_catchments = buff_reaches %>% head()
  # ## Vis reach_catchments ##
  # mapview(st_as_sf(reach_catchments))+
  #   mapview(p1_nhd_reaches_along_NHM[p1_nhd_reaches_along_NHM$comid %in% reach_catchments$featureid,],
  #           color = 'red')
  # mapview(reach_catchments)+mapview(raster(raster_cm))
  
  raster_cm <- rast(adf_file_path)
  
  ## checks
  ## check vector geometries 
  if(any(!st_is_valid(reach_catchments))){
    vector_sf <- st_make_valid(reach_catchments)
    message('shp geometries fixed')
  } 
  
  print(st_crs(reach_catchments))
  print(st_crs(raster_cm))
  
  ## match crs
  if(!st_crs(raster_cm) == st_crs(reach_catchments)){
    message('crs are different. Transforming ...')
    vector_sf <- st_transform(reach_catchments, crs = st_crs(raster_cm))
    if(st_crs(raster_cm) == st_crs(reach_catchments)){
      message('crs now aligned')}
  }else{
    message('crs are already aligned')
  }

## extract raster values
  raster_val_per_polygon <- terra::extract(raster_cm, vect(reach_catchments),
                                           fun = mean, weighted = TRUE,
                                           na.rm = TRUE) %>%
    as.data.frame()

  reach_catchments$weighted_mean_raster <- round(raster_val_per_polygon[,2], 3)

 #mapview(st_as_sf(reach_catchments), zcol = 'weighted_mean_raster')

  return(reach_catchments[c(comid_col,"weighted_mean_raster")])
  
}

tt <- process_depth_to_bedrock(adf_file_path = '1_fetch/in/Shangguan_dtb_cm_250m_clip 2/w001001.adf',
                             reach_catchments = nhd_catchments_along_nhm, 
                             comid_col = 'featureid')

ttp <- process_depth_to_bedrock(adf_file_path = '1_fetch/in/Shangguan_dtb_cm_250m_clip 2/w001001.adf',
                               reach_catchments = buff_reaches, 
                               comid_col = 'comid')

## Viewing catchments + polygons
mapview(tt, zcol = 'weighted_mean_raster', color = 'red')+mapview(ttp, zcol = 'weighted_mean_raster')

## TEMP WORK BELOW ###


## TERRA APPROACH
## clip depth to bedrock data to the boundary of reaches in drb
raster_mask_cm <- terra::mask(raster_cm, terra::vect(reach_catchments))
raster_mask_cm <- terra::mask(raster(raster_cm), reach_catchments)

## mapview of clip
mapview(raster(raster_mask_cm))+mapview(reach_catchments)

raster_val_mean <- terra::extract(rast(raster_mask_cm), terra::vect(reach_catchments), fun=mean) %>% as.data.frame()
raster_val_mean <- terra::extract(raster_mask_cm, reach_catchments, fun=mean) %>% as.data.frame()

reach_catchments$raster_val_mean <- raster_val_mean[,2]

install.packages('exactextractr')
library(exactextractr)



## extract values for each comid buffer


## Compare raster vs terra processing ##
## Terra
raster_cm <- rast(adf_file_path)
reach_catchments <- vect(nhd_catchments_along_nhm %>% head())
raster_mask_cm <- terra::mask(raster_cm, reach_catchments)
raster_val_mean <- terra::extract(raster_mask_cm, reach_catchments, fun = mean, weights = TRUE) %>% as.data.frame()
raster_val_mean
# > raster_val_mean
# ID  w001001
# 1  1 1463.756
# 2  2 1340.902
# 3  3 1376.819
# 4  4 1338.603
# 5  5 1463.000
# 6  6 1416.091

## raster
raster_cm <- raster(adf_file_path)
reach_catchments <- nhd_catchments_along_nhm %>% head()
raster_mask_cm <- terra::mask(raster_cm, reach_catchments)
raster_val_mean <- terra::extract(raster_mask_cm, reach_catchments, fun=mean, weights=TRUE) %>% as.data.frame()
raster_val_mean
# > raster_val_mean
# V1
# 1 1463.756
# 2 1340.902
# 3 1376.819
# 4 1338.603
# 5 1463.000
# 6 1416.091

reach_catchments %>% mutate(mean = exact_extract(raster_mask_cm, reach_catchments,'weighted_mean', weights = reach_catchments$areasqkm))
