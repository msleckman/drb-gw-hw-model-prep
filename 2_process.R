source("2_process/src/match_sites_to_reaches.R")

p2_targets_list <- list(
  
  # Match temperature monitoring locations to NHDPlusv2 flowline reaches.
  # Use the previously-matched NHM segment when matching sites to 
  # NHDv2 reaches; note that this will result in sites matched to NHDv2
  # reaches that are up to 10 km away from the site. 
  tar_target(
    p2_drb_temp_sites_w_segs,
    match_sites_to_reaches(nhd_lines = p1_nhd_reaches, 
                           sites = p1_drb_temp_sites_sf, 
                           comids_segs = p1_drb_comids_segs,
                           use_NHM_match = TRUE)
  )
  
)
