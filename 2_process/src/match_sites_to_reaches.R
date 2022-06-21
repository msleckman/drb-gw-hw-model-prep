#' Function to match site locations to NHDPlusV2 reaches
#' 
#' @param nhd_lines sf data frame containing NHDPlusv2 flowlines
#' @param sites sf data frame containing site locations
#' @param comids_segs data frame indicating the NHDPlusv2 COMIDs that intersect
#' each PRMS/NHM segment. Contains columns `segidnat` and `COMID`
#' @param search_radius the maximum radius in meters to match a point to a segment;
#' segments outside of search_radius will not match
#' @param use_NHM_match logical; For the temperature forecasting project, NHM
#' segments were preferentially matched to minimize the "fish distance" between 
#' the downstream end of a segment and the site location, as opposed to minimizing 
#' the "bird distance" to identify the absolute closest segment to a site. When 
#' matching site locations to NHD reaches, should we only consider the NHDPlusv2 
#' COMIDs that overlap the previously-matched NHM segment? Defaults to TRUE.
#'
#' @value A data frame that contains the same columns as the `sites` data,
#' with additional columns `COMID`, `bird_dist_to_comid_m` where 
#' `bird_dist_to_comid_m` is the distance (in meters) between the site and the
#' matched NHDv2 flowline.
#' 
match_sites_to_reaches <- function(nhd_lines, sites, comids_segs, search_radius = 500, 
                                   use_NHM_match = TRUE){
  
  message('Matching site locations to NHDPlusv2 flowlines...')
  
  if(use_NHM_match){
    
    # The sites are not matched to NHM segments using a simple spatial join;
    # rather, sites are matched to an NHM segment, preferring segments for 
    # which the downstream vertex/endpoint is close to the site location. As
    # a result, we may choose to only consider NHD reaches that overlay
    # the previously-matched segidnat when matching sites to NHDPlusV2 reaches. 
    # See subset_closest in delaware-model-prep for further details: 
    # https://github.com/USGS-R/delaware-model-prep/blob/main/2_observations/src/subset_closest.R
    
    sites_matched <- sites %>%
      split(.,.$segidnat) %>%
      lapply(., function(x){
        
        # find comids that overlap the segidnat and subset nhd_lines to only those comids
        comids <- comids_segs %>%
          filter(segidnat == unique(x$segidnat)) %>%
          pull(COMID)
        
        nhd_lines_subset <- nhd_lines %>%
          filter(comid %in% comids)
        
        # Match sites to subset flowlines (use a high value for search_radius, e.g.
        # 20 km, since some sites will be matched to NHM segments as far as 10 km away)
        sites_w_reach_ids <- get_site_nhd_flowlines(sites = x, 
                                                    nhd_lines = nhd_lines_subset,
                                                    search_radius = 20000)
      }) %>%
      bind_rows()
    
  } else {
    sites_matched <- get_site_nhd_flowlines(nhd_lines, sites, search_radius)
  }
  
  return(sites_matched)
  
}




#' @description Function to match site locations to NHDv2 flowline reaches 
#' 
#' @param nhd_lines sf data frame containing NHDPlusV2 flowlines
#' @param sites sf data frame containing site locations
#' @param search_radius the maximum radius in meters to match a point to a segment;
#' segments outside of search_radius will not match
#'
#' @value A data frame that contains the same columns as the `sites` data,
#' with additional columns `COMID` and `bird_dist_to_comid_m`, where 
#' `bird_dist_to_comid_m` is the distance (in meters) between the site and the
#' matched NHDv2 flowline.
#'
get_site_nhd_flowlines <- function(nhd_lines, sites, search_radius = 500){
  
  # Project reaches to Albers Equal Area Conic so that offsets returned by 
  # get_flowline_index are in meters rather than degrees
  nhd_lines_proj <- st_transform(nhd_lines, 5070)
  sites_proj <- st_transform(sites, 5070)
  
  # Below, precision indicates the resolution of measure precision (in meters)
  # in the output; since we are interested in a more accurate estimate of the 
  # `offset` distance between a point and the matched reach, set precision to 1 m.
  # Conduct initial search using a larger radius (search_radius*5) than specified 
  # to account for any uncertainty in the RANN::nn2 nearest neighbor search. 
  flowline_indices <- nhdplusTools::get_flowline_index(flines = nhd_lines_proj,
                                                       points = sites_proj,
                                                       max_matches = 3,
                                                       search_radius = units::set_units(search_radius*5,"m"),
                                                       precision = 1)
  
  # For each site, select one 'best' match by filtering the matched segments to
  # the segment with the smallest `offset`, or distance between the site and the
  # flowline. Then filter sites to include only those sites with segment matches
  # that are within the specified `search_radius`.
  flowline_indices_filtered <- flowline_indices %>%
    group_by(id) %>%
    arrange(offset) %>%
    slice(1) %>%
    select(id, COMID, offset) %>%
    rename(bird_dist_to_comid_m = offset) %>%
    filter(bird_dist_to_comid_m <= search_radius)
  
  # nhdplusTools returns an "id" column which is just an index from 1 to 
  # the number of sites. To later join to the site-ids, we need to add
  # a matching index column.
  sites <- rowid_to_column(sites, "id")
  
  # Rejoin flowline indices to original sites df
  sites_w_reach_ids <- sites %>%
    # only retain sites that got matched to flowlines and are 
    # within specified `search_radius`
    right_join(flowline_indices_filtered, by = "id") %>%
    mutate(across(c(segidnat, COMID), as.character)) %>%
    select(-id) %>%
    relocate(any_of(c("COMID", "bird_dist_to_comid_m")), .before = "geometry")
  
  return(sites_w_reach_ids)
  
}

