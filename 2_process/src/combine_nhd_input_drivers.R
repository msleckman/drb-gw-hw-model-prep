#' @title Combine NHD input drivers
#' 
#' @description 
#' Function to combine river-dl input drivers for NHDPlusv2 reaches, including
#' mean width, reach slope, reach elevation, and meteorological data.
#' 
#' @param nhd_flowlines sf object containing NHDPlusv2 flowline reaches. Must
#' contain columns "comid", "minelevsmo", "maxelevsmo", "slope"
#' @param prms_inputs data frame containing input drivers for each NHM segment.
#' Must include columns "segidnat", "seg_elev", "seg_slope", and "seg_width"
#' @param nhd_nhm_xwalk data frame that specifies how NHDPlusv2 COMIDs map
#' onto NHM segment identifiers. Must contain columns "COMID", "PRMS_segid",
#' and "segidnat".
#' @param climate_inputs data frame containing daily meteorological data to join
#' with the NHDPlusv2 static attributes. Must contain column "COMID".
#'
#' @returns 
#' Returns a data frame with one row per COMID in `nhd_flowlines`. Columns
#' indicate the NHM segment identifiers, the mean width for the NHD reach 
#' ("est_width_m"), the slope of the NHD reach ("slope"), the length-weighted
#' NHD slope for the COMIDs that make up each NHM segment ("slope_len_wtd_mean"), 
#' the length of the NHD reach ("lengthkm"), the min and max NHD reach elevation
#' ("min_elev_m" and "max_elev_m"), the elevation, slope, and mean width of the 
#' corresponding NHM segment ("seg_elev", "seg_slope", and "seg_width"), the
#' maximum NHD reach width among the COMIDs that make up each NHM segment
#' ("seg_width_max"), and the minimum NHD slope among the COMIDs that make up
#' each NHM segment ("seg_elev_min").
#' 
combine_nhd_input_drivers <- function(nhd_flowlines, prms_inputs, nhd_nhm_xwalk,
                                      climate_inputs){
  
  # Subset NHDPlusv2 flowlines to return the desired attributes, 
  # including elevation and slope
  nhd_attributes <- nhd_flowlines %>%
    sf::st_drop_geometry() %>%
    # transform elevation from cm to meters
    mutate(COMID = as.character(comid),
           min_elev_m = minelevsmo/100, 
           max_elev_m = maxelevsmo/100,
           slope = if_else(slope == -9998, NA_real_, slope)) 
  
  # Add NHM segment identifier to NHDv2 attributes table
  nhd_attributes_w_nhm <- nhd_attributes %>%
    left_join(y = nhd_nhm_xwalk, by = "COMID") %>%
    rename(subsegid = PRMS_segid)
  
  # Format NHDPlusv2 input data.
  # calculate length-weighted average slope for NHDv2 reaches associated
  # with each NHM reach. For simplicity, weight by the reach length rather
  # than another value-added attribute, slopelenkm, which represents the
  # length over which the NHDv2 attribute slope was computed.
  nhd_static_inputs <- nhd_attributes_w_nhm %>%
    group_by(subsegid) %>%
    mutate(slope_len_wtd_mean = weighted.mean(x = slope, w = lengthkm, na.rm = TRUE),
           seg_width_max = max(est_width_m, na.rm = TRUE), 
           seg_elev_min = min(min_elev_m, na.rm = TRUE)) %>%
    ungroup() %>%
    # join select attributes from PRMS-SNTemp
    left_join(y = prms_inputs, by = "segidnat") %>%
    select(COMID, segidnat, subsegid, est_width_m, slope, slope_len_wtd_mean, 
           lengthkm, min_elev_m, max_elev_m, seg_elev, seg_slope, seg_width, seg_width_max,
           seg_elev_min)
  
  message("Combining NHD static inputs with dynamic climate drivers...")
  
  # Combine dynamic meteorological inputs with static input data
  nhd_all_inputs <- nhd_static_inputs %>%
    left_join(y = mutate(climate_inputs, COMID = as.character(COMID)), 
                         by = "COMID") %>%
    relocate(segidnat, .after = COMID) %>%
    relocate(subsegid, .after = segidnat)
  
  return(nhd_all_inputs)

}

