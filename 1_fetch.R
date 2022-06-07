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
  )

  
)

  