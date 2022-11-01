source("2_process/src/subset_closest_nhd.R")
source('2_process/src/raster_in_polygon_weighted_mean.R')
source("2_process/src/estimate_mean_width.R")
source("2_process/src/write_data.R")
source("2_process/src/combine_nhd_input_drivers.R")
source("2_process/src/munge_split_temp_dat.R")
source("2_process/src/coarse_stratified_sediment_processing.R")
source("2_process/src/process_channel_confinement.R")
source("2_process/src/process_nhdv2_attr.R")

p2_targets_list <- list(

  # Subset NHDv2 reaches that overlap the NHM network to only include those 
  # that have a corresponding catchment (and meteorological data)
  tar_target(
    p2_dendritic_nhd_reaches_along_NHM_w_cats,
    p1_dendritic_nhd_reaches_along_NHM %>%
      filter(areasqkm > 0)
  ),
  
  # Match temperature monitoring locations to "mainstem" NHDPlusv2 flowline
  # reaches, i.e. those that NHD reaches that intersect the NHM river network.
  # Note that this site-to-segment matching procedure emulates the process
  # used in delaware-model-prep for modeling water temperature in the Delaware
  # River Basin. Here, match each site to an NHDPlusv2 reach, preferring reaches 
  # for which the downstream vertex (endpoint) is close to the site point.
  tar_target(
    p2_drb_temp_sites_w_segs,
    subset_closest_nhd(nhd_lines = p2_dendritic_nhd_reaches_along_NHM_w_cats,
                       sites = p1_drb_temp_sites_sf)
  ),
  
  # Match temperature observational time series to "mainstem" NHDPlusv2 flowline
  # reaches by COMID. Rename columns to match the column names used in the aggregated
  # temps data file in the temperature forecasting data release. 
  tar_target(
    p2_drb_temp_obs_w_segs,
    p1_drb_temp_obs %>%
      left_join(y = p2_drb_temp_sites_w_segs[,c("site_id","comid")], by = "site_id")
  ),
  
  # Resolve duplicate temperature observations and summarize the temperature data
  # to return one value for each COMID-date. By setting prioritize_nwis_sites to
  # FALSE we include all observations from non-NWIS sites in the summarized data.
  tar_target(
    p2_drb_temp_obs_by_comid,
    p2_drb_temp_obs_w_segs %>%
      munge_split_temp_dat(., prioritize_nwis_sites = FALSE) %>%
      rename(COMID = comid)
  ),
  
  # Save temperature observations mapped to NHDPlusv2 COMIDs as a zarr data store
  tar_target(
    p2_drb_temp_obs_w_segs_zarr,
    write_df_to_zarr(p2_drb_temp_obs_by_comid, 
                     index_cols = c("date", "COMID"), 
                     "2_process/out/drb_temp_observations_nhdv2.zarr"),
    format = "file"
  ),
  
  # Create buffer sf object of nhm reaches
  # Use xwalk nhd reaches along nhm and Dissolve all reaches to NHM scale
  tar_target(
    p2_buffered_nhm_reaches,
    ## Join with xwalk table to get PRMS_segids for nhm network
    p1_nhd_reaches_along_NHM %>% 
      mutate(COMID = as.character(comid)) %>%
      left_join(.,
                p1_drb_comids_all_tribs %>%
                  mutate(COMID = as.character(COMID)),
                by = 'COMID') %>%
      sf::st_make_valid() %>% 
      ## Dissolving by PRMS segid - old nrow = 3229, new nrow = 459 
      group_by(PRMS_segid) %>%
      dplyr::summarize(geometry = sf::st_union(geometry)) %>% 
      ## Buffer reach segments to 250 
      sf::st_buffer(.,dist = units::set_units(250, m)) %>% 
      ## creating new col with area of buffer - useful for downstream targets that uses buffered reaches
      mutate(total_reach_buffer_area_km2 = units::set_units(st_area(.), km^2)) %>% 
      relocate(geometry, .after = last_col())
  ),
  
  # Depth to bedrock processing
  ## Note: If you do not have Shangguan_dtb_cm_250m_clip_path data, you must grab 
  ## it from the caldera project folder. Dataset accessible on caldera in project 
  ## folder sub-dir: 1_fetch/in. scp to your local 1_fetch/in folder in this repo 
  ## in order to run this piece of pipeline original source: 
  ## http://globalchange.bnu.edu.cn/research/dtb.jsp. 
  ## Data was clipped to drb before getting added to caldera.
  
  # Reach -- depth_to_bedrock data for each nhm reach buffered at 250m  
  ## Note: In function, we transform the proj of vector to the raster (4326) to 
  ## perform weighted average. Retransform to 5070 after computation at end of code chunk.  
  tar_target(
    p2_depth_to_bedrock_reaches_along_nhm,
    raster_in_polygon_weighted_mean(raster = p1_depth_to_bedrock_tif,
                                    nhd_polygon_layer =  p2_buffered_nhm_reaches,
                                    feature_id = 'PRMS_segid', 
                                    weighted_mean_col_name = 'dtb_weighted_mean')
  ),
  
  # Catchment -- depth_to_bedrock data for each nhm upstream catchment 
  ## Note: In function, we transform the proj of vector to the raster (4326) to perform weighted average. Retransform to 5070 after computation at end of code chunk.  
  tar_target(
    p2_depth_to_bedrock_catchments_along_nhm_dissolved,
    raster_in_polygon_weighted_mean(raster = p1_depth_to_bedrock_tif,
                                    nhd_polygon_layer =  p1_nhm_catchments_dissolved,
                                    feature_id = 'PRMS_segid',
                                    weighted_mean_col_name  = 'dtb_weighted_mean') %>% 
      # append dtb value subsegid = 287_1 because this reach doesn't have an nhd catchment
      rbind(.,
            p2_depth_to_bedrock_reaches_along_nhm[p2_depth_to_bedrock_reaches_along_nhm$PRMS_segid == '287_1',]) 
  ),
  
  ## Soller coarse stratified Sediment processing to buffered-reach scale
  tar_target(
    p2_soller_coarse_sediment_reaches_nhm,
    coarse_sediment_area_calc(buffered_reaches_sf = p2_buffered_nhm_reaches,
                              buffered_reaches_area_col = 'total_reach_buffer_area_km2',
                              coarse_sediments_area_sf = p1_soller_coarse_sediment_drb_sf,
                              prms_col = 'PRMS_segid')
    ),
  
  # Process McManamay channel confinement dataset, including reaggregating
  # from NHDPlusv2 to NHM.
  tar_target(
    p2_confinement_mcmanamay,
    aggregate_mcmanamay_confinement(confinement_data = p1_confinement_mcmanamay, 
                                    nhd_nhm_xwalk = p1_drb_comids_segs, 
                                    force_min_width_m = 1,
                                    network = "nhm",
                                    nhm_identifier_col = "seg_id_nat")
  ),
  
  # Pull centroid of each reach within FACET stream network.
  # [Lauren] This target took ~1.3 hours to build on my local machine. It
  # takes a long time because the FACET stream network is so dense. We 
  # may look for ways to optimize this function in the future.
  tar_target(
    p2_facet_network_centroid,
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
    p2_confinement_facet,
    calculate_facet_confinement(facet_network = p2_facet_network_centroid, 
                                facet_width_col = "CW955mean_1D",
                                facet_floodplain_width_col = "FWmean_1D_FP",
                                nhd_catchment_polygons = p1_nhd_catchments,
                                nhd_nhm_xwalk = p1_drb_comids_segs,
                                network = "nhm",
                                nhm_identifier_col = "seg_id_nat")
  ),
  
  # Process Wieczorek NHDPlusv2 attributes referenced to cumulative upstream
  # area; returns object target of class "list". 
  # We are using the "TOT" (total cumulative drainage area) columns in the 
  # Wieczorek attribute data files to represent the cumulative upstream
  # attribute values rather than "ACC" (divergence-routed accumulate values). 
  # Note that if "TOT" is selected below, list elements for CAT_PPT 
  # (catchment-scale precip) and ACC_PPT (watershed-accumulated precip) will 
  # only contain the PRMS_segid column and so will functionally be omitted 
  # when creating the `p2_nhdv2_attr` target below. 
  tar_target(
    p2_nhdv2_attr_upstream,
    process_cumulative_nhdv2_attr(file_path = p1_sb_attributes_downloaded_csvs,
                                  segs_w_comids = p1_drb_comids_down,
                                  cols = c("TOT")),
    pattern = map(p1_sb_attributes_downloaded_csvs),
    iteration = "list"
  ),
  
  # Process Wieczorek NHDPlusv2 attributes scaled to the catchment that directly
  # drains into each NHM segment; returns object target of class "list" that is 
  # nested and contains the aggregated data as well as a separate NA diagnostics
  # data table for each NHDv2 attribute.
  tar_target(
    p2_nhdv2_attr_catchment,
    process_catchment_nhdv2_attr(file_path = p1_sb_attributes_downloaded_csvs,
                                 vars_table = p1_sb_attributes,
                                 segs_w_comids = p1_drb_comids_all_tribs,
                                 nhd_lines = p1_nhd_reaches),
    pattern = map(p1_sb_attributes_downloaded_csvs),
    iteration = "list"
  ),
  
  # Create combined NHDv2 attribute data frame that includes both the cumulative 
  # upstream and catchment-scale values that have been aggregated to the NHM scale.
  tar_target(
    p2_nhdv2_attr,
    create_nhdv2_attr_table(p2_nhdv2_attr_upstream, p2_nhdv2_attr_catchment)
  ),
  
  # Estimate mean width for each "mainstem" NHDv2 reach. 
  # Note that one NHM segment, segidnat 1721 (subsegid 287_1) is not included
  # in the dendritic nhd reaches w cats data frame because the only COMID that
  # intersects the segment (COMID 4188275) does not have an NHD catchment area. 
  # So in addition to estimating width for the COMIDs represented in 
  # p2_dendritic_nhd_reaches_along_NHM_w_cats, we also want to estimate width 
  # for COMID 4188275.
  tar_target(
    p2_nhd_mainstem_reaches_w_width,
    {
      nhd_lines <- bind_rows(p2_dendritic_nhd_reaches_along_NHM_w_cats,
                             filter(p1_dendritic_nhd_reaches_along_NHM,
                                    comid == "4188275"))
      estimate_mean_width(nhd_lines, 
                          estimation_method = 'nwis',
                          network_pos_variable = 'arbolate_sum',
                          ref_gages = p1_ref_gages_sf)
    }
  ),

  # Pull static segment attributes from PRMS SNTemp model driver data
  tar_target(
    p2_static_inputs_prms,
    p1_sntemp_input_output %>%
      group_by(seg_id_nat) %>%
      summarize(seg_elev = unique(seg_elev),
                seg_slope = unique(seg_slope),
                seg_width = mean(seg_width, na.rm = TRUE)) %>%
      mutate(seg_id_nat = as.character(seg_id_nat)) %>%
      select(seg_id_nat, seg_elev, seg_slope, seg_width)
  ),
  
  # Pull dynamic segment attributes from PRMS SNTemp model driver data
  tar_target(
    p2_dynamic_inputs_prms,
    p1_sntemp_input_output %>%
      mutate(seg_id_nat = as.character(seg_id_nat)) %>%
      select(seg_id_nat, date, seginc_potet)
  ),
  
  # Subset the DRB meteorological data to only include the NHDPlusv2 catchments 
  # (COMID) that intersect the NHM segments. `subset_nc_to_comid()` originally
  # developed by Jeff Sadler as part of the PGDL-DO project:
  # https://github.com/USGS-R/drb-do-ml/blob/main/2_process/src/subset_nc_to_comid.py
  # The resulting target is ~0.6 GB.
  tar_target(
    p2_met_data_nhd_mainstem_reaches,
    {
      reticulate::source_python("2_process/src/subset_nc_to_comid.py")
      subset_nc_to_comids(p1_drb_nhd_gridmet, 
                          p2_dendritic_nhd_reaches_along_NHM_w_cats$comid) %>%
        as_tibble() %>%
        relocate(c(COMID,time), .before = "tmmx") %>%
        # format dates
        mutate(date = lubridate::as_date(time, tz = "UTC")) %>%
        # convert gridmet precip units from inches to meters, and temperature
        # units from degrees Farenheit to degrees Celsius. 
        # Note that we average the daily minimum and maximum temperatures from
        # gridmet and call that "seg_tave_air", which we assume corresponds 
        # approximately to the PRMS-SNTemp variable with the same name. 
        mutate(tmmean = rowMeans(select(., c(tmmn,tmmx))),
               seg_tave_air = ((tmmean - 32) * (5/9)),
               seg_rain = pr * 0.0254) %>%
        # rename gridmet columns to conform to PRMS-SNTemp names used in river-dl
        select(COMID, date, seg_tave_air, srad, seg_rain) %>%
        rename(seginc_swrad = srad)
    }
  ),
  
  # Compile river-dl static input drivers at NHDv2 resolution, including river 
  # width (m) slope (unitless), and min/max elevation (m). Note that these input 
  # drivers represent "mainstem" NHDv2 reaches only (i.e., those reaches that 
  # intersect the NHM fabric). Some of the columns in p2_static_input_drivers_nhd 
  # are meant to facilitate comparison/EDA between segment attributes at the NHM and 
  # NHDPlusv2 scales (i.e., seg_slope ~ slope_len_wtd_mean; seg_elev ~ seg_elev_min; 
  # seg_width ~ seg_width_max). 
  tar_target(
    p2_static_inputs_nhd,
    prepare_nhd_static_inputs(nhd_flowlines = p2_nhd_mainstem_reaches_w_width,
                              prms_inputs = p2_static_inputs_prms,
                              nhd_nhm_xwalk = p1_drb_comids_all_tribs)
  ),
  
  # Format NHD-scale static input drivers
  tar_target(
    p2_static_inputs_nhd_formatted,
    p2_static_inputs_nhd %>%
      select(COMID, seg_id_nat, subsegid, 
             est_width_m, min_elev_m, slope) %>%
      rename(seg_width_empirical = est_width_m,
             seg_elev = min_elev_m,
             seg_slope = slope,
             seg_id_nat = seg_id_nat)
  ),
  
  # Combine NHD-scale static input drivers with dynamic climate drivers. 
  # Note that there is currently no pot ET analog variable at the NHD-scale,
  # so we are using seginc_potet from PRMS-SNTemp and assuming that there is 
  # not much intra-segment variation in potential ET among NHD reaches that 
  # contribute to a given NHM segment. The resulting target contains 15,800 
  # days of climate data across each of 3,173 COMIDs = 50,133,400 total rows.
  tar_target(
    p2_input_drivers_nhd,
    combine_nhd_input_drivers(nhd_static_inputs = p2_static_inputs_nhd_formatted, 
                              climate_inputs = p2_met_data_nhd_mainstem_reaches,
                              prms_dynamic_inputs = p2_dynamic_inputs_prms,
                              earliest_date = "1979-01-01",
                              latest_date = "2022-04-04")
  ),
  
  # Save river-dl input drivers at NHDv2 resolution as a zarr data store
  tar_target(
    p2_input_drivers_nhd_zarr,
    write_df_to_zarr(p2_input_drivers_nhd, 
                     index_cols = c("date", "COMID"), 
                     "2_process/out/nhdv2_inputs_io.zarr"),
    format = "file"
  ),
  
  # Format NHM-scale attributes, including empirical width
  # seg_id_nat's 1437, 1442, and 1485 each have two subsegid's attached to one 
  # seg_id_nat (e.g. 3_1 and 3_2 in the case of seg_id_nat 1437). Since we only
  # consider the seg_id_nat values when running river-dl, here we assume that 
  # the sub-segments all have the same width value, which is equal to the maximum
  # width between the two contributing subsegid's. 
  tar_target(
    p2_static_inputs_nhm_formatted,
    p2_static_inputs_nhd_formatted %>%
      group_by(seg_id_nat) %>%
      summarize(seg_width_empirical = max(seg_width_empirical),
                .groups = "drop") %>%
      mutate(seg_id_nat = as.integer(seg_id_nat))
  ),
  
  # Save a feather file that contains the formatted NHM-scale attributes
  tar_target(
    p2_static_inputs_nhm_formatted_feather,
    write_feather(p2_static_inputs_nhm_formatted, "2_process/out/nhm_attributes.feather"),
    format = "file"
  )
  
)
