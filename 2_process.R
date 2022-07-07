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
                        network_pos_variable = 'arbolate_sum')
  )
  
)
