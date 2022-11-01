#' @title Process McManamay channel confinement data
#' 
#' @description 
#' Function to process McManamay and DeRolph (2018) channel confinement dataset,
#' including to aggregate the values from NHDv2 reaches to NHM segments.
#' 
#' @param confinement_data data frame containing McManamay channel confinement
#' data. Must include columns "COMID", "VBL", "RL", "RWA", and "VBA".
#' @param nhd_nhm_xwalk data frame that specifies how NHDPlusv2 COMIDs map
#' onto NHM segment identifiers. Must contain columns "COMID", "PRMS_segid",
#' and "seg_id_nat".
#' @param force_min_width_m numeric value that indicates what the minimum 
#' width value used in channel confinement calculations should be. All width
#' values below this minimum width will be set to this value. Optional;
#' defaults to zero, which amounts to this setting being ignored.
#' @param network character string indicating the requested resolution of the
#' channel confinement values. Options include "nhdv2" or "nhm". If "nhm", 
#' the river and valley bottom lengths and areas will be aggregated to NHM
#' flowlines and the channel confinement categories will be assigned using
#' the same criteria as described in the McManamay and DeRolph (2018) Scientific
#' Data paper: https://doi.org/10.1038/sdata.2019.17.
#' @param nhm_identifier_col if `network` is "nhm" identify which column in
#' `confinement_data` contains the unique identifier to use for aggregation (e.g.
#' "PRMS_segid" or "seg_id_nat"). Defaults to "seg_id_nat".
#' 
#' @returns 
#' Returns a data frame with one row per reach/segment and columns representing
#' the valley bottom length: river length ratio, the valley bottom area: river
#' area ratio, and the confinement category that was assigned ("Confinement_calc").
#' 
aggregate_mcmanamay_confinement <- function(confinement_data, 
                                            nhd_nhm_xwalk, 
                                            force_min_width_m = 0,
                                            network = "nhdv2", 
                                            nhm_identifier_col = "seg_id_nat"){
  
  # Check that the value for `network` matches one of two options we expect
  if(!network %in% c("nhdv2","nhm")){
    stop(paste0("The network argument accepts 'nhdv2' or 'nhm'. Please check ",
                "that the requested network matches one of these two options."))
  }
  
  # Subset full confinement dataset to the NHDPlusv2 COMID's of interest.
  confinement_subset <- confinement_data %>%
    filter(COMID %in% nhd_nhm_xwalk$COMID)
  
  # Note that McManamay and DeRolph develop a categorization scheme 
  # (described in https://doi.org/10.1038/sdata.2019.17) to assign
  # confinement categories: "unconfined," "moderately confined," or "confined."
  # Even though they only present confinement classes, we could also calculate
  # a numeric value for confinement based on the information given.
  confinement_nhd <- confinement_subset %>%
    rename(reach_length_km = RL,
           valley_bottom_length = VBL,
           reach_area = RWA,
           valley_bottom_area = VBA,
           vbl_rl_ratio = VBL_RL_R,
           vba_ra_ratio = VBA_RWA_R) %>%
    # Back-calculate river width and floodplain width from the values given.
    mutate(river_width_m = reach_area/(reach_length_km*1000),
           floodplain_width_m = if_else(valley_bottom_length == 0, 
                                        NA_real_,
                                        valley_bottom_area/(valley_bottom_length*1000))) %>%
    # Replace all values of width below the user-specified value in `force_min_width_m`.
    # If `force_min_width_m` equals zero, this step is functionally omitted. 
    mutate(river_width_m = if_else(river_width_m < force_min_width_m, force_min_width_m, river_width_m)) %>%
    # Note that the calculation below is equal to vba_ra_ratio/vbl_rl_ratio. If
    # denominator is equal to zero, just assign NA for confinement.
    mutate(confinement_calc_mcmanamay = if_else(river_width_m == 0 | floodplain_width_m == 0, 
                                                NA_real_,
                                                floodplain_width_m/river_width_m))
  
  if(network == "nhdv2"){
    confinement_nhd_out <- confinement_nhd %>%
      select(COMID, reach_length_km, river_width_m, floodplain_width_m, confinement_calc_mcmanamay)
    return(confinement_nhd_out)
  }
  
  if(network == "nhm"){
    # Join subsetted confinement dataset to NHM segment identifiers and format
    # column names in preparation for aggregation step.
    confinement_w_nhm_segs <- confinement_nhd %>%
      rename(reach_length_km_comid = reach_length_km,
             confinement_calc_comid = confinement_calc_mcmanamay,
             river_width_comid = river_width_m,
             floodplain_width_comid = floodplain_width_m) %>%
      mutate(COMID = as.character(COMID)) %>%
      left_join(nhd_nhm_xwalk, by = "COMID")
    
    # Group data by NHM segment and calculate a length-weighted average of the
    # confinement values among individual COMIDs that make up the NHM segment.
    confinement_nhm <- confinement_w_nhm_segs %>%
      group_by(.data[[nhm_identifier_col]]) %>%
      summarize(
        reach_length_km = sum(reach_length_km_comid),
        lengthkm_mcmanamay_is_na = sum(reach_length_km_comid[is.na(confinement_calc_comid)]),
        prop_reach_w_mcmanamay = 1-(lengthkm_mcmanamay_is_na/reach_length_km),
        confinement_calc_mcmanamay = if_else(prop_reach_w_mcmanamay > 0,
                                             weighted.mean(x = confinement_calc_comid,
                                                           w = reach_length_km_comid,
                                                           na.rm = TRUE),
                                             NA_real_),
        .groups = "drop") %>%
      mutate(
        flag_mcmanamay = if_else(prop_reach_w_mcmanamay < 0.7, 
                             paste0("Note that <70% of the NHM segment is ",
                                    "covered by a COMID with McManamay confinement data."),
                             NA_character_)
      )

    return(confinement_nhm)
  }
  
}



#' @title Estimate channel confinement from FACET geomorphic data
#' 
#' @description
#' Function to aggregate geomorphic data derived from 3-m LiDAR from the
#' DRB FACET dataset to estimate channel confinement for individual 
#' NHDPlusv2 or NHM segments.
#' 
#' @details 
#' See FACET metadata for further details about the geomorphic dataset
#' obtained from the FACET model and the geomorphic metric columns used
#' for `facet_width_col` and `facet_floodplain_width_col` here.
#' https://doi.org/10.5066/P9RQJPT1.
#' 
#' @param facet_network sf linestring object containing the FACET stream network.
#' Must contain columns "UniqueID", Magnitude", and "USContArea" in addition to 
#' the columns defined in `facet_width_col` and `facet_floodplain_width_col`.
#' @param facet_width_col character string indicating which column from the FACET
#' dataset should be used to represent channel width. Defaults to "CW955mean_1D", 
#' channel width, mean, for values <95th and >5th percentile (within the reach). 
#' @param facet_floodplain_width_col character string indicating which column from
#' the FACET dataset should be used to represent floodplain width. Defaults to
#' "FWmean_1D_FP", total floodplain width (fpwid_1d), mean. 
#' @param nhd_catchment_polygons sf polygon object containing the NHDPlusv2
#' catchments. Must contain column "COMID".
#' @param nhd_nhm_xwalk data frame that specifies how NHDPlusv2 COMIDs map
#' onto NHM segment identifiers. Must contain columns "COMID", "PRMS_segid",
#' and "seg_id_nat".
#' @param network character string indicating the requested resolution of the
#' channel confinement values. Options include "nhdv2" or "nhm".
#' @param nhm_identifier_col if `network` is "nhm" identify which column in
#' `confinement_data` contains the unique identifier to use for aggregation (e.g.
#' "PRMS_segid" or "seg_id_nat"). Defaults to "seg_id_nat".
#' @param show_warnings logical; should any warnings that arise during the 
#' spatial join be printed to the console? Defaults to FALSE.
#' 
#' @returns 
#' Returns a data frame containing estimates of channel confinement, defined
#' as floodplain width/channel width (unitless). 
#' 
calculate_facet_confinement <- function(facet_network, 
                                        facet_width_col = "CW955mean_1D", 
                                        facet_floodplain_width_col = "FWmean_1D_FP", 
                                        nhd_catchment_polygons, 
                                        nhd_nhm_xwalk,
                                        network = "nhdv2",
                                        nhm_identifier_col = "seg_id_nat",
                                        show_warnings = FALSE){
  
  # Check that the value for `network` matches one of two options we expect.
  if(!network %in% c("nhdv2","nhm")){
    stop(paste0("The network argument accepts 'nhdv2' or 'nhm'. Please check ",
                "that the requested network matches one of these two options."))
  }
  
  # Spatially join the FACET stream network with the NHDPlusv2 catchments.
  # Then, for each NHDPlusv2 catchment, subset the FACET segment with the 
  # largest Shreve magnitude (if multiple segments with the same magnitude, 
  # break a tie using the upstream area). Capture all warning messages (e.g.
  # "attribute variables assumed to be spatially constant throughout all 
  # geometries" during the st_intersection step). If `show_warnings` is TRUE, 
  # print any warning messages to the console. Otherwise, hide them.
  facet_nhd <- withCallingHandlers({
    sf::st_intersection(x = facet_network, y = nhd_catchment_polygons) %>%
      group_by(COMID) %>%
      arrange(desc(Magnitude), desc(USContArea)) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(COMID = as.character(COMID))
  }, warning = function(w) {
    if(!show_warnings) invokeRestart("muffleWarning")
  })


  # Retain selected columns in the joined FACET dataset and subset
  # the data to only include the requested COMIDs.
  cols_to_keep <- c("COMID", "UniqueID", "HUC4", "Magnitude", "USContArea", 
                    facet_width_col, facet_floodplain_width_col)
  
  facet_nhd_out <- nhd_nhm_xwalk %>%
    left_join(y = facet_nhd, by = "COMID") %>%
    select(any_of(cols_to_keep)) %>%
    sf::st_drop_geometry() %>%
    rename(channel_width := {{facet_width_col}}) %>%
    rename(floodplain_width := {{facet_floodplain_width_col}}) %>%
    mutate(confinement_calc_facet = floodplain_width/channel_width)
  
  # If network = "nhdv2", return estimated confinement values.
  if(network == "nhdv2"){
    return(facet_nhd_out)
  }
  
  
  if(network == "nhm"){
    # If network = "nhm", download requested COMIDs and subset reach
    # length column. Then join NHD-scale metrics to NHM segment IDs.
    nhd_lengths <- nhdplusTools::get_nhdplus(comid = unique(nhd_nhm_xwalk$COMID),
                                             realization = "flowline") %>%
      mutate(COMID = as.character(comid)) %>%
      sf::st_drop_geometry() %>%
      select(COMID, lengthkm) %>%
      rename(lengthkm_comid = lengthkm)
    
    # Bind COMID, NHDPlusv2 reach lengths, and NHDPlusv2-scale FACET
    # metrics together.
    facet_nhm <- nhd_nhm_xwalk %>%
      left_join(y = nhd_lengths, by = "COMID") %>%
      left_join(y = facet_nhd_out, by = "COMID") %>%
      rename(confinement_calc_comid = confinement_calc_facet)
    
    # For each NHM segment, calculate a length-weighted average of the NHDPlusv2-scale
    # channel confinement values. In addition, add a flag for any NHM segments where 
    # less than 70% of the segment is covered by a COMID with FACET values.
    facet_nhm_out <- facet_nhm %>%
      group_by(.data[[nhm_identifier_col]]) %>%
      summarize(
        lengthkm = sum(lengthkm_comid),
        lengthkm_facet_is_na = sum(lengthkm_comid[is.na(confinement_calc_comid)]),
        prop_reach_w_facet = 1-(lengthkm_facet_is_na/lengthkm),
        confinement_calc_facet = if_else(prop_reach_w_facet > 0,
                                         weighted.mean(x = confinement_calc_comid, 
                                                       w = lengthkm_comid, 
                                                       na.rm = TRUE),
                                         NA_real_)) %>%
      mutate(
        flag_facet = if_else(prop_reach_w_facet < 0.7, 
                             paste0("Note that <70% of the NHM segment is ",
                                    "covered by a COMID with FACET data."),
                             NA_character_)
      )
    
    return(facet_nhm_out)
  }
  
}



#' @title Get centroid of stream reaches
#' 
#' @description 
#' Function to find the centroid of each linestring object representing
#' one stream reach within the network.
#' 
#' @param network sf linestring object containing the stream network. Must
#' contain column "UniqueID".
#' @param show_warnings logical; should any warnings that arise during the 
#' spatial join be printed to the console? Defaults to FALSE.
#' 
#' @returns 
#' Returns sf point object containing one centroid point per linestring in `network`.
#' 
get_reach_centroids <- function(network, show_warnings = FALSE){
  
  message("Finding the centroid for each reach within the network. This may take awhile...")
  
  network_pts_at_centroids <- network %>%
    split(., .$UniqueID) %>%
    lapply(., function(x){
      # 1) cast the reach to points:
      reach_pts <- withCallingHandlers({
        sf::st_cast(x, "POINT")
      }, warning = function(w){
        if(!show_warnings) invokeRestart("muffleWarning")
      })
      
      # 2) Find the centroid of each reach:
      reach_centroid <- withCallingHandlers({
        sf::st_centroid(x)
      }, warning = function(w){
        if(!show_warnings) invokeRestart("muffleWarning")
      })
      
      # 3) Snap reach centroid to the reach:
      pt_at_centroid <- reach_pts[which.min(sf::st_distance(reach_centroid, reach_pts)),]
      return(pt_at_centroid)
    }) %>%
    bind_rows()

  return(network_pts_at_centroids)
}





