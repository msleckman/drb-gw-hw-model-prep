source("1_fetch/src/download_nhdplus_flowlines.R")

p1_targets_list <- list(
  
  # Read in the NHM - NHDv2 crosswalk file
  tar_target(
    p1_GFv1_NHDv2_xwalk_csv,
    # For now, download the NHM/PRMS - NHDv2 crosswalk table with divergent
    # reaches omitted from the GW/HW Quest "Data" folder on Sharepoint and
    # place in a directory called 1_fetch/in. The crosswalk table will 
    # ultimately be pulled directly from the drb-network-prep github repo.
    "1_fetch/in/GFv1_NHDv2_xwalk.csv",
    format = "file"
  ),
  
  tar_target(
    p1_GFv1_NHDv2_xwalk,
    read_csv(p1_GFv1_NHDv2_xwalk_csv,col_types = cols(.default = "c"))
  ),
  
  # Reshape crosswalk table to return all COMIDs that overlap each NHM segment
  tar_target(
    p1_drb_comids_segs, 
    p1_GFv1_NHDv2_xwalk %>%
      select(PRMS_segid, comid_seg) %>% 
      tidyr::separate_rows(comid_seg,sep=";") %>% 
      rename(COMID = comid_seg)
  ),
  
  # Reshape crosswalk table to return all COMIDs that drain to each NHM segment
  tar_target(
    p1_drb_comids_all_tribs, 
    p1_GFv1_NHDv2_xwalk %>%
      select(PRMS_segid, comid_cat) %>% 
      tidyr::separate_rows(comid_cat,sep=";") %>% 
      rename(COMID = comid_cat)
  ),
  
  # Use crosswalk table to fetch just the NHDv2 reaches that overlap the NHM network
  tar_target(
    p1_nhd_reaches_along_NHM,
    download_nhdplus_flowlines(p1_drb_comids_segs$COMID)
  ),
  
  # Use crosswalk table to fetch all NHDv2 reaches in the DRB 
  tar_target(
    p1_nhd_reaches,
    download_nhdplus_flowlines(p1_drb_comids_all_tribs$COMID)
  ),
  
  # Manually download temperature site locations from ScienceBase using the
  # commented-out code below and place the downloaded zip file in 1_fetch/in. 
  # Note that you'll be prompted for your username and password and will need 
  # authorization to download the temperature files while the data release is 
  # still in process:
  # sbtools::authenticate_sb()
  # download_sb_file(sb_id = "623e54c4d34e915b67d83580",
  #                 file_name = "study_monitoring_sites.zip",
  #                 out_dir = "1_fetch/in")
  tar_target(
    p1_drb_temp_sites_shp,
    {
      file_names <- unzip(zipfile = "1_fetch/in/study_monitoring_sites.zip", 
                          exdir = "1_fetch/out", 
                          overwrite = TRUE)
      grep(".shp",file_names, value = TRUE, ignore.case = TRUE)
    },
    format = "file"
  ),
  
  # Read in temperature site locations
  tar_target(
    p1_drb_temp_sites_sf,
    sf::read_sf(p1_drb_temp_sites_shp, crs = 4326)
  ),
  
  # Manually download unaggregated temperature observations from ScienceBase using
  # the commented-out code below and place the downloaded zip file in 1_fetch/in. 
  # Note that you'll be prompted for your username and password and will need 
  # authorization to download the temperature files while the data release is 
  # still in process:
  # sbtools::authenticate_sb()
  # download_sb_file(sb_id = "623e550ad34e915b67d8366e",
  #                 file_name = "unaggregated_temperature_observations_drb.zip",
  #                 out_dir = "1_fetch/in")
  tar_target(
    p1_drb_temp_obs_csv,
    unzip(zipfile = "1_fetch/in/unaggregated_temperature_observations_drb.zip", 
          overwrite = TRUE, 
          exdir = "1_fetch/out"),
    format = "file"
  ),
  
  # Read in temperature observations
  tar_target(
    p1_drb_temp_obs,
    read_csv(p1_drb_temp_obs_csv, col_types = list(seg_id_nat = "c"))
  )

  
)

  