## functions taken from drb-network-prep/get_nhdv2_flowlines.R 

get_nhdplusv2_catchments <- function(comid, crs = 4326){
  
  #' @description Function to download NHDPlusv2 catchments given a set of COMIDs
  #' @param comid character string or vector of character strings containing
  #' the common identifier (COMID or featureid) of the desired catchment(s).
  #' @param crs specified crs. default is 4326 because that is the projection of the output from nhdplusv2

  # Chunk desired COMIDs into groups, where each group has no more than
  # 50 COMID's to avoid timeout errors when downloading nhdplus subsets
  # using helper functions from nhdplusTools
  comid_df <- tibble(COMID = comid) %>%
    mutate(comid_n = row_number(),
           download_grp = ((comid_n -1) %/% 50) + 1)
  
  # Download catchments associated with desired COMIDs and return an
  # sf data frame containing the polygons and value-added attributes.
  # suppress messages "Found invalid geometry, attempting to fix" from
  # printing to the console.
  catchments <- comid_df %>%
    split(., .$download_grp) %>%
    lapply(., function(x){
      cats_sub <- suppressMessages(nhdplusTools::get_nhdplus(comid = x$COMID, 
                                                             realization = "catchment",
                                                             t_srs = crs))
    }) %>%
    bind_rows() 
  
  # Format variable names
  catchments_out <- catchments %>%
    rename_with(., toupper, id:shape_area) %>%
    rename(COMID = FEATUREID)
  
  # Inform user how many catchments were returned
  message(sprintf("Returning %s catchments of %s total COMID's requested",
                  format(length(catchments_out$ID), big.mark=","), 
                  format(length(comid), big.mark=",")))
  
  return(catchments_out)
  
}

