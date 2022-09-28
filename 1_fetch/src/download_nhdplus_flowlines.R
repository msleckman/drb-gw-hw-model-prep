#' @title Download NHDPlusv2 flowlines.
#'
#' @description
#' Function to download NHDPlusv2 flowlines given a set of COMIDs
#' 
#' @details 
#' This function was pulled from drb-network-prep/get_nhdv2_flowlines.R:
#' https://github.com/USGS-R/drb-network-prep/blob/main/1_fetch/src/get_nhdplusv2.R 
#' 
#' @param comid character string or vector of character strings containing
#' the common identifier (COMID) of the desired flowline(s)
#' 
download_nhdplus_flowlines <- function(comid, crs = 4326){
  
  # Chunk desired COMIDs into groups, where each group has no more than
  # 50 COMID's to avoid timeout errors when downloading nhdplus subsets
  # using helper functions from nhdplusTools
  comid_df <- tibble(COMID = comid) %>%
    mutate(comid_n = row_number(),
           download_grp = ((comid_n -1) %/% 50) + 1)
  
  # Download flowlines associated with desired COMIDs and return an sf
  # data frame containing the linestrings and value-added attributes 
  flowlines <- comid_df %>%
    split(., .$download_grp) %>%
    lapply(., function(x){
      flines_sub <- nhdplusTools::get_nhdplus(comid = x$COMID, 
                                              realization = "flowline",
                                              t_srs = crs)
      # format certain columns to allow merging chunked flowlines into a single
      # data frame
      flines_sub_out <- flines_sub %>%
        mutate(across(c(lakefract, surfarea, rareahload,hwnodesqkm), as.character))
    }) %>%
    bind_rows() 
  
  return(flowlines)
  
}

