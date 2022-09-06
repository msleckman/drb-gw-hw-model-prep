source("1_fetch/src/download_nhdplus_flowlines.R")
source('1_fetch/src/sb_read_filter_by_comids.R')
source("1_fetch/src/download_nhdplus_catchments.R")
source("1_fetch/src/download_file.R")
source("1_fetch/src/read_netcdf.R")

p1_targets_list <- list(
  
  # Important! This pipeline uses two versions of the NHM-NHDv2 crosswalk table
  # for different purposes. All targets pertaining to the "dendritic" version
  # of the crosswalk table will include the string "dendritic" in the name.
  # Crosswalk #1: Read in the NHM - NHDv2 crosswalk file that contains all 
  # NHDPlusv2 COMIDs in the DRB (including divergent reaches and reaches 
  # without catchments)
  tar_target(
    p1_GFv1_NHDv2_xwalk,
    read_csv(GFv1_NHDv2_xwalk_url,col_types = cols(.default = "c"))
  ),
  
  # Reshape crosswalk table to return all COMIDs that drain to each NHM segment.
  tar_target(
    p1_drb_comids_all_tribs, 
    p1_GFv1_NHDv2_xwalk %>%
      select(PRMS_segid, segidnat, comid_cat) %>% 
      tidyr::separate_rows(comid_cat,sep=";") %>% 
      rename(COMID = comid_cat)
  ),
  
  # Reshape crosswalk table to return all COMIDs that overlap each NHM segment.
  tar_target(
    p1_drb_comids_segs, 
    p1_GFv1_NHDv2_xwalk %>%
      select(PRMS_segid, segidnat, comid_seg) %>% 
      tidyr::separate_rows(comid_seg,sep=";") %>% 
      rename(COMID = comid_seg)
  ),
  
  # Use crosswalk table to fetch all NHDv2 reaches in the DRB. These COMIDs 
  # should be used for preparing feature data, including aggregating feature 
  # values from the NHD-scale to the NHM-scale and/or for deriving feature 
  # values from raster data.
  tar_target(
    p1_nhd_reaches,
    download_nhdplus_flowlines(p1_drb_comids_all_tribs$COMID)
  ),
  
  # Use crosswalk table to fetch just the NHDv2 reaches that overlap the NHM network.
  tar_target(
    p1_nhd_reaches_along_NHM,
    download_nhdplus_flowlines(p1_drb_comids_segs$COMID)
  ),
  
  # Download all NHDPlusv2 catchments in the DRB (may take awhile to run).
  tar_target(
    p1_nhd_catchments,
    get_nhdplusv2_catchments(comid = p1_nhd_reaches$comid)
  ),
  
  # Dissolve NHDPlusv2 catchments to create a single catchment polygon for
  # each NHM segment (analogous to HRU).
  tar_target(
    p1_nhm_catchments_dissolved,
    {sf_use_s2(FALSE)
      left_join(p1_nhd_catchments %>% mutate(COMID = as.character(COMID)),
                p1_drb_comids_all_tribs %>% mutate(COMID = as.character(COMID)),
                by = 'COMID') %>%
        group_by(PRMS_segid) %>%
        dplyr::summarize(geometry = sf::st_union(geometry))
    }
  ),
  
  # Crosswalk #2: Read in the NHM - NHDv2 crosswalk file that corresponds to 
  # the *dendritic* network only (i.e., divergent reaches have been omitted). 
  # This version of the crosswalk table was used to build the network distance
  # matrix in drb-network-prep, and so includes the COMIDs that we will make 
  # predictions on in the NHD-downscaling set of experiments. 
  tar_target(
    p1_GFv1_NHDv2_xwalk_dendritic,
    read_csv(GFv1_NHDv2_xwalk_dendritic_url,col_types = cols(.default = "c"))
  ),
  
  # Reshape dendritic crosswalk table to return all COMIDs that overlap each 
  # NHM segment
  tar_target(
    p1_drb_comids_dendritic_segs, 
    p1_GFv1_NHDv2_xwalk_dendritic %>%
      select(PRMS_segid, segidnat, comid_seg) %>% 
      tidyr::separate_rows(comid_seg,sep=";") %>% 
      rename(COMID = comid_seg)
  ),
  
  # Use crosswalk table to fetch just the dendritic NHDv2 reaches that overlap
  # the NHM network
  tar_target(
    p1_dendritic_nhd_reaches_along_NHM,
    download_nhdplus_flowlines(p1_drb_comids_dendritic_segs$COMID)
  ),
  
  # Download temperature site locations from ScienceBase
  tar_target(
    p1_drb_temp_sites_zip,
    download_sb_file(sb_id = "623e54c4d34e915b67d83580",
                     file_name = "study_monitoring_sites.zip",
                     out_dir = "1_fetch/out"),
    format = "file"
  ),

  # Unzip downloaded temperature site locations and save shp file in 1_fetch/out
  # TODO: Getting the following error message right now: error 1 in extracting from
  # zip file (I cannot manually unzip the downloaded file, either)
  tar_target(
    p1_drb_temp_sites_shp,
    {
      file_names <- unzip(zipfile = p1_drb_temp_sites_zip, 
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
  
  # Download unaggregated temperature observations from ScienceBase 
  tar_target(
    p1_drb_temp_obs_zip,
    download_sb_file(sb_id = "623e550ad34e915b67d8366e",
                     file_name = "unaggregated_temperature_observations_drb.zip",
                     out_dir = "1_fetch/out"),
    format = "file"
  ),

  # Unzip unaggregated temperature observations and save csv file in 1_fetch/out
  tar_target(
    p1_drb_temp_obs_csv,
    unzip(zipfile = p1_drb_temp_obs_zip, 
          overwrite = TRUE, 
          exdir = "1_fetch/out"),
    format = "file"
  ),
  
  # Read in temperature observations
  tar_target(
    p1_drb_temp_obs,
    read_csv(p1_drb_temp_obs_csv, col_types = list(seg_id_nat = "c"))
  ),
  
  # Manually download PRMS-SNTemp model driver data from ScienceBase 
  # using the commented-out code below and place the download zip file in 
  # 1_fetch/in. Note that you'll be prompted for your username and password
  # and will need authorization to download the model driver data while the 
  # data release is still in process:
  #sbtools::authenticate_sb()
  #download_sb_file(sb_id = "623e5587d34e915b67d83806",
  #                 file_name = "uncal_sntemp_input_output.nc.zip",
  #                 out_dir = "1_fetch/in")
  tar_target(
    p1_sntemp_input_output_zip,
    "1_fetch/in/uncal_sntemp_input_output.nc.zip",
    format = "file"
  ),
  tar_target(
    p1_sntemp_input_output_nc,
    unzip(zipfile = p1_sntemp_input_output_zip, 
          overwrite = TRUE, 
          exdir = "1_fetch/out"),
    format = "file"
  ),
  
  # Read in PRMS-SNTemp model driver data
  tar_target(
    p1_sntemp_input_output,
    read_netcdf(p1_sntemp_input_output_nc)
  ),
  
  # Read in meteorological data aggregated to NHDPlusV2 catchments for the 
  # DRB (prepped in https://github.com/USGS-R/drb_gridmet_tools). Note that
  # the DRB met data file must be stored in 1_fetch/in. If working outside
  # of tallgrass/caldera, this file will need to be downloaded from the
  # PGDL-DO project's S3 bucket and manually placed in 1_fetch/in.
  tar_target(
    p1_drb_nhd_gridmet,
    "1_fetch/in/drb_climate_2022_06_14.nc",
    format = "file"
  ),
  
  # Download ref-gages v0.6 to help QC matching NWIS sites to NHDv2 flowlines
  tar_target(
    p1_ref_gages_geojson,
    download_file(url = 'https://github.com/internetofwater/ref_gages/releases/download/v0.6/usgs_nldi_gages.geojson',
                  fileout = "1_fetch/out/ref_gages.geojson",
                  mode = "wb", quiet = TRUE)
  ),
  
  # Read in ref-gages v0.6 as an sf object
  tar_target(
    p1_ref_gages_sf,
    sf::st_read(p1_ref_gages_geojson, quiet = TRUE) %>%
      mutate(COMID_refgages = as.character(nhdpv2_COMID))
  ),
  
  # STATSGO SOIL Characteristics
  # get selected child items nhdv2 STATSGO Soil Characteristics,
  # including 1) texture and 2) layer attributes.
  tar_target(
    p1_selected_statsgo_sbid_children,
    sbtools::item_list_children(sb_id = nhd_statsgo_parent_sbid) %>% 
      Filter(function(x){str_detect(x[['title']],'Text|Layer')},
             .)
  ),
  
  # Download selected CONUS STATSGO datasets from Science base
  tar_target(
    p1_download_statsgo_text_layer_zip,
    lapply(p1_selected_statsgo_sbid_children,
           function(x){download_sb_file(sb_id = x$id,
                                        out_dir = '1_fetch/out/statsgo',
                                        file_name = NULL,
                                        overwrite_file = TRUE)}
           ) %>%
      unlist(),
    format = 'file'
  ),

  # Combine statsgo TEXT and Layer Attributes for CAT and TOT and filter to drb
  tar_target(
    p1_statsgo_soil_df,
    sb_read_filter_by_comids(data_path = '1_fetch/out/statsgo',
                             comid = p1_drb_comids_all_tribs$COMID,
                             sb_comid_col = 'COMID',
                             selected_cols_contains = c("KFACT","KFACT_UP","NO10AVE",
                                                        "NO4AVE","SILTAVE","CLAYAVE",
                                                        "SANDAVE",'WTDEP'),
                             cbind = TRUE)
  ),
  
  # Track depth to bedrock raster dataset in 1_fetch/in
  tar_target(
    p1_depth_to_bedrock_tif,
    Shangguan_dtb_cm_250m_clip_path,
    format = "file"
  )
  
)
