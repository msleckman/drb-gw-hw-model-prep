source("2a_process_nhd_downscaling/src/write_data.R")
source('2b_process_nhm_groundwater/src/raster_in_polygon_weighted_mean.R')
source("2b_process_nhm_groundwater/src/coarse_stratified_sediment_processing.R")
source("2b_process_nhm_groundwater/src/process_channel_confinement.R")
source("2b_process_nhm_groundwater/src/process_nhdv2_attr.R")

p2b_targets_list <- list(

  # Create buffer sf object of nhm reaches
  # Use xwalk nhd reaches along nhm and dissolve all reaches to NHM scale.
  tar_target(
    p2b_buffered_nhm_reaches,
    # Join with xwalk table to get NHM segment identifiers
    p1_nhd_reaches_along_NHM %>% 
      mutate(COMID = as.character(comid)) %>%
      left_join(y = p1_drb_comids_all_tribs,
                by = 'COMID') %>%
      sf::st_make_valid() %>% 
      # dissolve reaches by seg_id_nat - nrow, old: = 3229; nrow, new: = 456 
      group_by(seg_id_nat) %>%
      dplyr::summarize(geometry = sf::st_union(geometry)) %>% 
      # create 250 meter buffer around each reach
      sf::st_buffer(., dist = units::set_units(250, m)) %>% 
      # create new col with area of buffer - useful for downstream targets that use buffered reaches
      mutate(total_reach_buffer_area_km2 = units::set_units(sf::st_area(.), km^2)) %>% 
      relocate(geometry, .after = last_col())
  ),
  
  # Depth to bedrock processing
  # Reach -- depth_to_bedrock data for each buffered NHM reach.  
  # Note: In function, we transform the proj of vector to the raster (4326) to 
  # perform weighted average. Retransform to 5070 after computation at end of code chunk.  
  tar_target(
    p2b_depth_to_bedrock_reaches_along_nhm,
    raster_in_polygon_weighted_mean(raster = p1_depth_to_bedrock_tif,
                                    nhd_polygon_layer =  p2b_buffered_nhm_reaches,
                                    feature_id = 'seg_id_nat', 
                                    weighted_mean_col_name = 'dtb_weighted_mean')
  ),
  
  # Catchment -- depth_to_bedrock data for each NHM catchment. 
  # Note: In function, we transform the proj of vector to the raster (4326) to 
  # perform weighted average. Retransform to 5070 after computation at end of code chunk.  
  tar_target(
    p2b_depth_to_bedrock_catchments_along_nhm_dissolved,
    raster_in_polygon_weighted_mean(raster = p1_depth_to_bedrock_tif,
                                    nhd_polygon_layer =  p1_nhm_catchments_dissolved,
                                    feature_id = 'seg_id_nat',
                                    weighted_mean_col_name  = 'dtb_weighted_mean') %>% 
      # append dtb value subsegid = 287_1/segidnat 1721. Because this reach doesn't have 
      # an NHDPlusv2 catchment (and the catchment is presumably very small for this short
      # reach), we assume that the buffered reach value also applies to the catchment.
      rbind(.,
            p2b_depth_to_bedrock_reaches_along_nhm[p2b_depth_to_bedrock_reaches_along_nhm$seg_id_nat == '1721',]) 
  ),
  
  # Process Soller et al. coarse stratified sediments to the scale of the buffered NHM segments.
  tar_target(
    p2b_soller_coarse_sediment_reaches_nhm,
    coarse_sediment_area_calc(buffered_reaches_sf = p2b_buffered_nhm_reaches,
                              buffered_reaches_area_col = 'total_reach_buffer_area_km2',
                              coarse_sediments_area_sf = p1_soller_coarse_sediment_drb_sf,
                              prms_col = 'seg_id_nat')
  ),
  
  # Process McManamay channel confinement dataset, including reaggregating from
  # NHDPlusv2 to NHM. Replace any width values in the McManamay and DeRolph (2018)
  # dataset that are less than 1m with 1m (the minimum width allowed). See further
  # discussion in https://github.com/USGS-R/drb-gw-hw-model-prep/issues/41.
  tar_target(
    p2b_confinement_mcmanamay,
    aggregate_mcmanamay_confinement(confinement_data = p1_confinement_mcmanamay, 
                                    nhd_nhm_xwalk = p1_drb_comids_segs, 
                                    force_min_width_m = 1,
                                    network = "nhm",
                                    nhm_identifier_col = "seg_id_nat")
  ),
  
  # The processed McManamay channel confinement data are missing confinement
  # values for 17 NHM segments. The objective of this target is to impute
  # the NA values with confinement estimates from the nearest upstream 
  # segment with a non-NA confinement value. 
  tar_target(
    p2b_confinement_mcmanamay_filled,
    refine_from_neighbors(attr_df = p2b_confinement_mcmanamay, 
                          nhm_identifier_col = "seg_id_nat",
                          attr_name = "confinement_calc_mcmanamay",
                          reach_distances = p1_nhm_distance_matrix,
                          neighbors = "nearest")
  ),
  
  # Pull centroid of each reach within FACET stream network.
  # [Lauren] This target took ~1.3 hours to build on my local machine. It
  # takes a long time because the FACET stream network is so dense. We 
  # may look for ways to optimize this function in the future.
  tar_target(
    p2b_facet_network_centroid,
    get_reach_centroids(p1_facet_network)
  ),
  
  # Process FACET DRB geomorphometry dataset by first Spatially joining the FACET
  # stream reaches that have their center within an NHDPlusv2 catchment. For each
  # NHDPlusv2 catchment, subset the FACET segment with the largest shreve magnitude
  # (if multiple with the same magnitude, break a tie with upstream area). Select 
  # columns for mean channel width (between 5th and 95th percentiles within reach) 
  # and floodplain width as described in https://doi.org/10.1088/1748-9326/ac6e47. 
  # If aggregation to NHM segments is requested, FACET floodplain width and channel
  # width values for each NHDPlusv2 COMID are summarized as a length-weighted mean 
  # before calculating channel confinement.
  tar_target(
    p2b_confinement_facet,
    calculate_facet_confinement(facet_network = p2b_facet_network_centroid, 
                                facet_width_col = "CW955mean_1D",
                                facet_floodplain_width_col = "FWmean_1D_FP",
                                nhd_catchment_polygons = p1_nhd_catchments,
                                nhd_nhm_xwalk = p1_drb_comids_segs,
                                network = "nhm",
                                nhm_identifier_col = "seg_id_nat")
  ),
  
  # The processed FACET channel confinement data are missing confinement
  # values for 93 NHM segments. The objective of this target is to impute
  # the NA values with confinement estimates from the nearest upstream 
  # segment with a non-NA confinement value. 
  tar_target(
    p2b_confinement_facet_filled,
    refine_from_neighbors(attr_df = p2b_confinement_facet, 
                          nhm_identifier_col = "seg_id_nat",
                          attr_name = "confinement_calc_facet",
                          reach_distances = p1_nhm_distance_matrix,
                          neighbors = "nearest")
  ),
  
  # Process Wieczorek NHDPlusv2 attributes referenced to cumulative upstream
  # area; returns object target of class "list". 
  # We are using the "TOT" (total cumulative drainage area) columns in the 
  # Wieczorek attribute data files to represent the cumulative upstream
  # attribute values rather than "ACC" (divergence-routed accumulate values). 
  # Note that if "TOT" is selected below, list elements for CAT_PPT 
  # (catchment-scale precip) and ACC_PPT (watershed-accumulated precip) will 
  # only contain the PRMS_segid column and so will functionally be omitted 
  # when creating the `p2b_nhdv2_attr` target below. 
  tar_target(
    p2b_nhdv2_attr_upstream,
    process_cumulative_nhdv2_attr(file_path = p1_sb_attributes_downloaded_csvs,
                                  segs_w_comids = p1_drb_comids_down %>%
                                    # for the split segments, just take the more downstream segment
                                    filter(!PRMS_segid %in% c("3_1", "8_1", "51_1")) %>%
                                    select(seg_id_nat, COMID),
                                  cols = c("TOT")),
    pattern = map(p1_sb_attributes_downloaded_csvs),
    iteration = "list"
  ),
  
  # Process Wieczorek NHDPlusv2 attributes scaled to the catchment that directly
  # drains into each NHM segment; returns object target of class "list" that is 
  # nested and contains the aggregated data as well as a separate NA diagnostics
  # data table for each NHDv2 attribute.
  tar_target(
    p2b_nhdv2_attr_catchment,
    process_catchment_nhdv2_attr(file_path = p1_sb_attributes_downloaded_csvs,
                                 vars_table = p1_sb_attributes,
                                 segs_w_comids = select(p1_drb_comids_all_tribs, seg_id_nat, COMID),
                                 nhm_identifier_col = "seg_id_nat",
                                 nhd_lines = p1_nhd_reaches),
    pattern = map(p1_sb_attributes_downloaded_csvs),
    iteration = "list"
  ),
  
  # Create combined NHDv2 attribute data frame that includes both the cumulative 
  # upstream and catchment-scale values that have been aggregated to the NHM scale.
  tar_target(
    p2b_nhdv2_attr,
    create_nhdv2_attr_table(attr_data_upstream = p2b_nhdv2_attr_upstream, 
                            attr_data_catchment = p2b_nhdv2_attr_catchment, 
                            nhm_identifier_col = "seg_id_nat")
  ),
  
  # Format NHM-scale attributes, including empirical width
  # seg_id_nat's 1437, 1442, and 1485 each have two subsegid's attached to one 
  # seg_id_nat (e.g. 3_1 and 3_2 in the case of seg_id_nat 1437). Since we only
  # consider the seg_id_nat values when running river-dl, here we assume that 
  # the sub-segments all have the same width value, which is equal to the maximum
  # width between the two contributing subsegid's. 
  tar_target(
    p2b_static_inputs_nhm_formatted,
    p2a_static_inputs_nhd_formatted %>%
      group_by(seg_id_nat) %>%
      summarize(seg_width_empirical = max(seg_width_empirical),
                .groups = "drop")
  ),
  
  # Combine NHM-scale river and catchment attributes into a single data frame
  tar_target(
    p2b_static_inputs_nhm_combined,
    p2b_static_inputs_nhm_formatted %>%
      left_join(y = p2a_static_inputs_prms, by = "seg_id_nat") %>%
      left_join(y = p2b_soller_coarse_sediment_reaches_nhm  %>%
                  sf::st_drop_geometry() %>% 
                  select(seg_id_nat, cs_area_proportion) %>% 
                  rename(prop_reach_w_coarse_sediment = cs_area_proportion),
                by = "seg_id_nat") %>%
      left_join(y = p2b_depth_to_bedrock_reaches_along_nhm %>%
                  sf::st_drop_geometry() %>%
                  rename(dtb_weighted_mean_reach = dtb_weighted_mean),
                by = "seg_id_nat") %>%
      left_join(y = p2b_depth_to_bedrock_catchments_along_nhm_dissolved %>%
                  sf::st_drop_geometry() %>%
                  rename(dtb_weighted_mean_catchment = dtb_weighted_mean),
                by = "seg_id_nat") %>%
      left_join(y = rename(p2b_confinement_mcmanamay_filled, flag_gaps_mcmanamay = flag_gaps), 
                by = "seg_id_nat") %>%
      left_join(y = rename(p2b_confinement_facet_filled, flag_gaps_facet = flag_gaps), 
                by = "seg_id_nat") %>%
      left_join(y = p2b_nhdv2_attr, by = "seg_id_nat") %>%
      mutate(seg_id_nat = as.integer(seg_id_nat)) %>%
      left_join(y = select(p1_modflow_params, -COMID), by = "seg_id_nat") %>%
      left_join(y = select(p1_modflow_discharge, -COMID), by = "seg_id_nat")
  ),
  
  # Save a feather file that contains the formatted NHM-scale attributes
  tar_target(
    p2b_static_inputs_nhm_formatted_feather,
    write_feather(p2b_static_inputs_nhm_combined, sprintf("2b_process_nhm_groundwater/out/nhm_attributes_%s.feather", format(Sys.Date(), "%Y%m%d"))),
    format = "file"
  )
  
)
