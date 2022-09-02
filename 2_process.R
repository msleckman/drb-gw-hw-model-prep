source("2_process/src/subset_closest_nhd.R")
source('2_process/src/raster_in_polygon_weighted_mean.R')
source("2_process/src/estimate_mean_width.R")
source("2_process/src/write_data.R")

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
  
  # ## Dissolving all reaches to prms scale (joining with xwalk table to do this)
  tar_target(p2_buffered_nhd_reaches_along_nhm,
             p1_nhd_reaches_along_NHM %>% 
               mutate(COMID = as.character(comid)) %>%
               left_join(.,
                         p1_drb_comids_all_tribs %>%
                           mutate(COMID = as.character(COMID)), 
                         by = 'COMID') %>%
               st_make_valid() %>% 
               # Dissolving by PRMS segid
               group_by(PRMS_segid) %>%
               dplyr::summarize(geometry = sf::st_union(geometry)) %>% 
               st_buffer(., dist = units::set_units(250, m))
  ),
  
  # Depth to bedrock processing
  ## Note: If you do not have Shangguan_dtb_cm_250m_clip_path data, you must grab it from the caldera project folder. 
  ## Dataset accessible on caldera in project folder sub-dir: 1_fetch/in. scp to your local 1_fetch/in folder in this repo in order to run this piece of pipeline
  ## original source: http://globalchange.bnu.edu.cn/research/dtb.jsp. Data was clipped to drb before getting added to caldera.
  
  # Reach -- depth_to_bedrock data for each nhm reach buffered at 250m  
  tar_target(p2_depth_to_bedrock_reaches_along_nhm,
             raster_in_polygon_weighted_mean(raster = p1_depth_to_bedrock_tif,
                                             nhd_polygon_layer =  p2_buffered_nhd_reaches_along_nhm,
                                             feature_id = 'PRMS_segid', 
                                             weighted_mean_col_name = 'dtb_weighted_mean')
  ),
  
  # Catchment -- depth_to_bedrock data for each nhm upstream catchment 
  tar_target(p2_depth_to_bedrock_catchments_along_nhm_dissolved,
             raster_in_polygon_weighted_mean(raster = p1_depth_to_bedrock_tif,
                                             nhd_polygon_layer =  p1_nhm_catchments_dissolved,
                                             feature_id = 'PRMS_segid',
                                             weighted_mean_col_name  = 'dtb_weighted_mean') %>% 
               ## tacking on 287_1 dtb value for reach because it 287_1 doesn't have a catchment 
               rbind(.,
                     p2_depth_to_bedrock_reaches_along_nhm[p2_depth_to_bedrock_reaches_along_nhm$PRMS_segid == '287_1',])
  ),
  
  
  # Estimate mean width for each "mainstem" NHDv2 reach 
  tar_target(
    p2_nhd_mainstem_reaches_w_width,
    estimate_mean_width(p1_nhd_reaches_along_NHM, 
                        estimation_method = 'nwis',
                        network_pos_variable = 'arbolate_sum',
                        ref_gages = p1_ref_gages_sf)
  ),
  
  # Pull static segment attributes from PRMS SNTemp model driver data
  tar_target(
    p2_input_drivers_prms,
    p1_sntemp_input_output %>%
      group_by(seg_id_nat) %>%
      summarize(seg_elev = unique(seg_elev),
                seg_slope = unique(seg_slope),
                seg_width = mean(seg_width, na.rm = TRUE)) %>%
      mutate(segidnat = as.character(seg_id_nat)) %>%
      select(segidnat, seg_elev, seg_slope, seg_width)
  ),
  
  # Compile river-dl input drivers at NHDv2 resolution, including river width
  # (meters), slope (unitless), and min/max elevation (transformed to meters
  # from cm). Note that these input drivers represent "mainstem" NHDv2 reaches 
  # only (i.e., those reaches that intersect the NHM fabric).
  tar_target(
    p2_input_drivers_nhd,
    p2_nhd_mainstem_reaches_w_width %>%
      sf::st_drop_geometry() %>%
      mutate(COMID = as.character(comid),
             min_elev_m = minelevsmo/100, 
             max_elev_m = maxelevsmo/100,
             slope = if_else(slope == -9998, NA_real_, slope)) %>%
      # add NHM segment identifier segidnat
      left_join(y = p1_drb_comids_all_tribs, by = "COMID") %>%
      rename(subsegid = PRMS_segid) %>%
      # calculate length-weighted average slope for NHDv2 reaches associated
      # with each NHM reach. For simplicity, weight by the reach length rather
      # than another value-added attribute, slopelenkm, which represents the
      # length over which the NHDv2 attribute slope was computed.
      group_by(subsegid) %>%
      mutate(slope_len_wtd_mean = weighted.mean(x = slope, w = lengthkm, na.rm = TRUE),
             seg_width_max = max(est_width_m, na.rm = TRUE), 
             seg_elev_min = min(min_elev_m, na.rm = TRUE)) %>%
      ungroup() %>%
      # join select attributes from PRMS-SNTemp
      left_join(y = p2_input_drivers_prms, by = "segidnat") %>%
      select(COMID, segidnat, subsegid, est_width_m, slope, slopelenkm, slope_len_wtd_mean, 
             lengthkm, min_elev_m, max_elev_m, seg_elev, seg_slope, seg_width, seg_width_max,
             seg_elev_min) 
  ),
  
  # Save river-dl input drivers at NHDv2 resolution as a feather file
  tar_target(
    p2_input_drivers_nhd_zarr,
    write_df_to_zarr(p2_input_drivers_nhd, 
                     index_cols = c("COMID"), 
                     "2_process/out/nhdv2_inputs_io.zarr"),
    format = "file"
  )
  
)
