library(sf)
library(dplyr)
library(terra)


process_depth_to_bedrock(adf_file_path){}

adf_file_path <- '1_fetch/in/Shangguan_dtb_cm_250m_clip 2/w001001.adf'

#huc12 <- nhdplusTools::get_huc12(AOI = p1_nhd_reaches)
polygon_sf <- huc12 %>% head()
raster <- raster(adf_file_path)
vector_sf <- polygon_sf

## check vector geometries 
if(any(!st_is_valid(vector_sf))){
  vector_sf <- st_make_valid(vector_sf)
  message('shp geometries fixed')
} 

## match crs
if(!st_crs(raster) == st_crs(vector_sf)){
  message('crs are different. Transforming ...')
  vector_sf <- st_transform(vector_sf, crs = st_crs(raster))
  if(st_crs(raster) == st_crs(vector_sf)){
    message('crs now aligned')}
}else{
  message('crs are already aligned')
}

raster_mask <- terra::mask(raster, vector_sf)
