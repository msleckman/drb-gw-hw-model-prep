source("2_process/src/subset_closest_nhd.R")
source('2_process/src/raster_in_polygon_weighted_mean.R')

p2_targets_list <- list(

  # Match temperature monitoring locations to "mainstem" NHDPlusv2 flowline
  # reaches, i.e. those that NHD reaches that intersect the NHM river network.
  # Note that this site-to-segment matching procedure emulates the process
  # used in delaware-model-prep for modeling water temperature in the Delaware
  # River Basin. Here, match each site to an NHDPlusv2 reach, preferring reaches 
  # for which the downstream vertex (endpoint) is close to the site point.
  tar_target(
    p2_drb_temp_sites_w_segs,
    subset_closest_nhd(nhd_lines = p1_nhd_reaches_along_NHM,
                       sites = p1_drb_temp_sites_sf)
  ),
  
  tar_target(p2_buffered_nhd_reaches,
              st_buffer(p1_nhd_reaches, dist = 250)
  ),

  tar_target(p2_buffered_nhd_reaches_along_nhm,
              st_buffer(p1_nhd_reaches_along_NHM, dist = 250)
  ),

  ## Catchment 
  # tar_target(p2_depth_to_bedrock_catchments,
  #            raster_in_polygon_weighted_mean(raster = '1_fetch/in/Shangguan_dtb_cm_250m_clip 2/w001001.adf',
  #                                            nhd_polygon_layer =  p1_nhd_catchments,
  #                                            comid_col = 'COMID')
  # ),
   
  tar_target(p2_depth_to_bedrock_catchments_along_nhm,
             raster_in_polygon_weighted_mean(raster = '1_fetch/in/Shangguan_dtb_cm_250m_clip 2/w001001.adf',
                                             nhd_polygon_layer =  p1_nhd_catchments_along_nhm,
                                             comid_col = 'COMID')
  ),

  tar_target(p2_depth_to_bedrock_reaches,
             raster_in_polygon_weighted_mean(raster = '1_fetch/in/Shangguan_dtb_cm_250m_clip 2/w001001.adf',
                                             nhd_polygon_layer =  p2_buffered_nhd_reaches,
                                             comid_col = 'comid')
  ),
  
  tar_target(p2_depth_to_bedrock_reaches_along_nhm,
             raster_in_polygon_weighted_mean(raster = '1_fetch/in/Shangguan_dtb_cm_250m_clip 2/w001001.adf',
                                             nhd_polygon_layer =  p2_buffered_nhd_reaches_along_nhm,
                                             comid_col = 'comid')
  )
)
