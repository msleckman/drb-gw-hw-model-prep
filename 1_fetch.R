source("1_fetch/src/download_nhdplus_flowlines.R")
source('1_fetch/src/sb_read_filter_by_comids.R')

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
  
  # STATSGO SOIL Characteristics
  ## get selected child items nhdv2 STATSGO Soil Characteristics
  ## 1) Text attributes and 2) Layer attributes
  tar_target(
    p1_selected_statsgo_sbid_children,
    sbtools::item_list_children(sb_id = nhd_statsgo_parent_sbid) %>% 
      Filter(function(x){str_detect(x[['title']],'Text|Layer')},
             .)
    ),
  
  ## download selected CONUS STATSGO datasets from Science base
  tar_target(
    p1_download_statsgo_text_layer_zip,
    lapply(p1_selected_statsgo_sbid_children,
           function(x){sbtools::item_file_download(x$id,
                                                   dest_dir = '1_fetch/out/statsgo',
                                                   overwrite_file = TRUE)}
           ) %>%
      unlist(),
    format = 'file'
    ),
  
  ## Combine statsgo TEXT and Layer Attributes for CAT and TOT and filter to drb
  tar_target(
    p1_statsgo_soil_df,
    sb_read_filter_by_comids(data_path = '1_fetch/out/statsgo',
                             comid = p1_drb_comids_all_tribs$COMID,
                             sb_comid_col = 'COMID',
                             selected_cols_contains = c("KFACT","KFACT_UP","NO10AVE",
                                                        "NO4AVE","SILTAVE","CLAYAVE",
                                                        "SANDAVE",'WTDEP'),
                             cbind = TRUE)
  )
)
