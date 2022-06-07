#' Function to download NHDPlus flowlines given a set of COMIDs
#' 
#' @param comid character string or vector of character strings containing
#' the common identifier (COMID) of the desired flowline(s)
#' 
download_nhdplus_flowlines <- function(comid){
  
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
                                              realization = "flowline")
      # format certain columns to allow merging chunked flowlines into a single
      # data frame
      flines_sub_out <- flines_sub %>%
        mutate(across(c(lakefract, surfarea, rareahload,hwnodesqkm), as.character))
    }) %>%
    bind_rows()
  
  return(flowlines)
  
}

