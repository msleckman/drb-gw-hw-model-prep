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
#' @param network character string indicating the requested resolution of the
#' channel confinement values. Options include "nhdv2" or "nhm". If "nhm", 
#' the river and valley bottom lengths and areas will be aggregated to NHM
#' flowlines and the channel confinement categories will be assigned using
#' the same criteria as described in the McManamay and DeRolph (2018) Scientific
#' Data paper: https://doi.org/10.1038/sdata.2019.17.
#' 
#' @returns 
#' Returns a data frame with one row per reach/segment and columns representing
#' the valley bottom length: river length ratio, the valley bottom area: river
#' area ratio, and the confinement category that was assigned ("Confinement_calc").
#' 
aggregate_mcmanamay_confinement <- function(confinement_data, nhd_nhm_xwalk, network = "nhdv2"){
  
  # Check that the value for `network` matches one of two options we expect
  if(!network %in% c("nhdv2","nhm")){
    stop(paste0("The network argument accepts 'nhdv2' or 'nhm'. Please check ",
                "that the requested network matches one of these two options."))
  }
  
  # Subset full confinement dataset to the NHDPlusv2 COMID's of interest
  confinement_subset <- confinement_data %>%
    filter(COMID %in% nhd_nhm_xwalk$COMID)
  
  # Join subsetted confinement dataset to NHM segment identifiers
  confinement_w_nhm_segs <- confinement_subset %>%
    mutate(COMID = as.character(COMID)) %>%
    left_join(nhd_nhm_xwalk, by = "COMID")
  
  # Group data by NHM segment and sum reach length, valley bottom length, 
  # reach area, and valley bottom area from individual COMIDs that comprise
  # each NHM segment. 
  confinement_nhm <- confinement_w_nhm_segs %>%
    group_by(seg_id_nat) %>%
    summarize(reach_length = sum(RL),
              valley_bottom_length = sum(VBL),
              reach_area = sum(RWA),
              valley_bottom_area = sum(VBA)) %>%
    ungroup() %>%
    mutate(vbl_rl_ratio = valley_bottom_length/reach_length,
           vba_ra_ratio = ifelse((valley_bottom_area == 0 & reach_area == 0), 
                                 0, (valley_bottom_area/reach_area))) %>%
    # Now use McManamay categorization scheme (described in https://doi.org/10.1038/sdata.2019.17)
    # to assign confinement categories, "unconfined," "moderately confined," or "confined."
    mutate(Confinement_calc = case_when(
      vbl_rl_ratio > 0.50 & vba_ra_ratio > 4 ~ "Unconfined",
      vbl_rl_ratio > 0.50 & vba_ra_ratio < 4 & vba_ra_ratio > 2 ~ "Mod Confined",
      vba_ra_ratio > 4 & vbl_rl_ratio > 0.25 & vbl_rl_ratio < 0.5 ~ "Mod Confined",
      TRUE ~ "Confined"
    ))
  
  # Return processed confinement dataset depending on the network requested
  if(network == "nhdv2"){
    return(confinement_subset)
  }else{
    return(confinement_nhm)
  }
  
}
