source("2_process/src/subset_closest_nhd.R")
source("2_process/src/estimate_mean_width.R")

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
  
  # Estimate mean width for each NHDv2 reach 
  tar_target(
    p2_nhd_reaches_w_width,
    estimate_mean_width(p1_nhd_reaches, 
                        estimation_method = 'nwis',
                        network_pos_variable = 'arbolate_sum',
                        ref_gages = p1_ref_gages_sf)
  ),
  
  # Compile river-dl input drivers at NHDv2 resolution, including river 
  # width (meters), slope (unitless), and min/max elevation (transformed
  # to meters from cm)
  tar_target(
    p2_input_drivers_nhd,
    p2_nhd_reaches_w_width %>%
      sf::st_drop_geometry() %>%
      mutate(COMID = as.character(comid),
             min_elev_m = minelevsmo/100, 
             max_elev_m = maxelevsmo/100,
             slope = if_else(slope == -9998, NA_real_,slope)) %>%
      left_join(p1_drb_comids_all_tribs, by = "COMID") %>%
      # calculate length-weighted average slope for NHDv2 reaches associated
      # with each NHM reach. For simplicity, weight by the reach length rather
      # than another value-added attribute, slopelenkm, which represents the
      # length over which the NHDv2 attribute slope was computed.
      group_by(segidnat) %>%
      mutate(slope_len_wtd_mean = weighted.mean(x = slope, w = lengthkm, na.rm = TRUE)) %>%
      ungroup() %>%
      select(COMID, segidnat, PRMS_segid, est_width_m, slope, slopelenkm, 
             slope_len_wtd_mean, lengthkm, min_elev_m, max_elev_m) 
  )
  
)
