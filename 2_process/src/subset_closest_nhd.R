#' Return the upstream/downstream vertices of NHDPlusv2 linestrings
#' 
#' @param nhd_lines sf linestring containing the NHDPlusv2 flowlines. Must
#' contain column COMID indicating the unique reach identifier
#' 
#' @return returns sf multipoint object containing the upstream and downstream
#' vertices for each COMID. The upstream vertex is denoted by reach_location = 0
#' and the downstream vertex is denoted by reach_location = 1.
#' 
return_nhd_vertices <- function(nhd_lines){
  
  # Reformat variable names to lowercase so that
  # code steps below will work regardless of case
  nhd_lines_proj <- nhd_lines %>%
    rename_with(.,tolower,everything()) %>%
    sf::st_transform(5070)
  
  message(sprintf("Returning flowline vertices for %s requested COMIDs...",
                  format(length(nhd_lines_proj$comid), big.mark=",")))
  
  vertices <- nhd_lines_proj %>%
    split(.,.$comid) %>%
    lapply(.,function(comid){
      # Return the upstream vertex associated with each NHDv2 reach
      vertex_up <- comid %>%
        sf::st_line_sample(type = 'regular', sample = 0) %>%
        sf::st_as_sf() %>%
        mutate(COMID = comid$comid,
               reach_location = 0,
               FROMNODE = comid$fromnode,
               TONODE = comid$tonode) %>%
        rename(geometry = x)
      
      # Return the downstream vertex associated with each NHDv2 reach
      vertex_down <- comid %>%
        sf::st_line_sample(type = 'regular', sample = 1) %>%
        sf::st_as_sf() %>%
        mutate(COMID = comid$comid,
               reach_location = 1,
               FROMNODE = comid$fromnode,
               TONODE = comid$tonode) %>%
        rename(geometry = x)
      
      # Return two vertices (up/down) for each NHDv2 reach
      net_vertices <- bind_rows(vertex_up, vertex_down) %>%
        sf::st_transform(.,st_crs(nhd_lines))
    }) %>%
    bind_rows()
  
  return(vertices)
}





#' Match each site point to an NHDPlusv2 reach, preferring reaches for which the
#' downstream vertex (endpoint) is close to the site point
#' 
#' @details This function originally created as part of the delaware-model-prep
#' pipeline and modified for use with NHDPlusv2 by L. Koenig. See:
#' https://github.com/USGS-R/delaware-model-prep/blob/main/2_observations/src/subset_closest.R
#'
#' @description Algorithm: Locate the nearest reach. If the nearest point is 
#' downstream, then return a match to that nearest reach. If the nearest reach 
#' endpoint (vertex) is upstream AND there is no branching at the upstream 
#' vertex, then return a match to the nearest upstream reach. Still report 
#' bird_dist_to_subseg_m as the distance to the initially identified nearest reach.
#'
#' @param nhd_lines sf linestring containing the NHDPlusv2 flowlines. Must
#' contain column COMID indicating the unique reach identifier.
#' @param sites sf POINT object with site identifier information in column siteid.
#' 
#' @return returns a data frame of site information, including the site_id and
#' the COMID of the most appropriate matched reach based on the algorithm details
#' described above. Data frame also contains the bird distance (in meters) to the 
#' originally matched reach and the fish distance (in meters) to the most 
#' appropriate matched reach. 
#' 
subset_closest_nhd <- function(nhd_lines, sites){

  # Get vertices associated with each reach in nhd_lines
  vertices <- return_nhd_vertices(nhd_lines)
  
  message('Matching site locations to NHDPlusv2 flowlines...')
  
  # Find the nearest NHDv2 segment to each point and return the 
  # 'bird distance' between the site and the nearest segment
  # [Alison] for guidance used for this implementation, see 
  # https://gis.stackexchange.com/questions/288570/find-nearest-point-along-polyline-using-sf-package-in-r
  nearest_segment <- nhd_lines[sf::st_nearest_feature(sites, nhd_lines), ]
  bird_dist_to_seg_m <- sf::st_length(
    sf::st_nearest_points(sites, nearest_segment, pairwise = TRUE))
  
  # Find the nearest NHDv2 vertex to each point and return the
  # 'bird distance' between the site and the nearest vertex
  nearest_vertex <- vertices[sf::st_nearest_feature(sites, vertices), ]
  bird_dist_to_vertex_m <- sf::st_length(
    sf::st_nearest_points(sites, nearest_vertex, pairwise = TRUE))
  
  nearest <- tibble(
    site_id = sites$siteid,
    # creates a cool horizontally nested tibble structure with multiple geometries:
    nearest_segment = nearest_segment,
    bird_dist_to_seg_m = units::drop_units(bird_dist_to_seg_m),
    # creates a cool horizontally nested tibble structure with multiple geometries:
    nearest_vertex = nearest_vertex,
    bird_dist_to_vertex_m = units::drop_units(bird_dist_to_vertex_m)) %>%
    sf::st_set_geometry(st_geometry(sites))
  
  # For each site, use nested conditionals to navigate a set of possibilities
  # for what the best site-reach match is. See notes below.
  crosswalk <- nearest %>%
    split(., .$site_id) %>%
    purrr::map_dfr(function(site_sf){
      # handle duplicated site_id
      if(nrow(site_sf) > 1) {
        bbox_diagonal <- sf::st_length(sf::st_sfc(
          sf::st_linestring(matrix(sf::st_bbox(site_sf), ncol=2, byrow=TRUE)), 
          crs=sf::st_crs(site_sf)))
        if(bbox_diagonal > units::set_units(1, m)) {
          warning(sprintf('site %s has diverse coordinates across databases, with bbox diagonal = %0.03f m', 
                          site_sf$site_id[1], bbox_diagonal))
        }
        # if a site id has different coordinates across databases, just keep the first one
        site_sf <- site_sf[1,] 
        }
      # start site processing by initializing a new non-sf tibble to return
      site <- tibble(site_id = site_sf$site_id) 
      
      # Calculate fish distances to the downstream vertex of the matched reach
      # and the downstream vertex of the nearest upstream reach 
      vertex_downstream <- vertices %>%
        filter(COMID == site_sf$nearest_segment$comid,
               reach_location == 1)
      vertex_upstream <- vertices %>%
        filter(COMID == site_sf$nearest_segment$comid,
               reach_location == 0)
      
      # cast matched segment to points (would work poorly if there was a big,
      # straight reach with no intermediate points)
      segment_as_points <- site_sf$nearest_segment %>%
        sf::st_geometry() %>%
        sf::st_cast(.,"POINT")
      point_loc_in_segment <- sf::st_nearest_feature(site_sf, segment_as_points)
      
      # confirm that the points are listed upstream to downstream
      stopifnot(sf::st_nearest_feature(vertex_upstream, segment_as_points) == 1) 
      
      # Calculate 'fish distance' between the site and the top/bottom of the 
      # matched reach
      fish_dist_upstream_m <- sf::st_combine(
        segment_as_points[1:point_loc_in_segment]) %>%
        sf::st_cast("LINESTRING") %>%
        sf::st_length() %>%
        units::drop_units()
      
      fish_dist_downstream_m <- sf::st_combine(
        segment_as_points[point_loc_in_segment:length(segment_as_points)]) %>%
        sf::st_cast("LINESTRING") %>%
        sf::st_length() %>%
        units::drop_units()
      
      # Decide which reach to use.
      # Because the model predicts values for the downstream point of each stream 
      # reach, we will sometimes want to use the reach upstream of the matched 
      # reach (if the site point was very close to the upstream point). So we 
      # need some conditionals.
      if(fish_dist_downstream_m < fish_dist_upstream_m) {
        # The nearest point (by fish distance) is downstream, so use the current reach
        site$comid <- site_sf$nearest_segment$comid
        site$bird_dist_to_subseg_m <- site_sf$bird_dist_to_seg_m
        site$fish_dist_to_outlet_m <- fish_dist_downstream_m
        } else {
          # The nearest point is upstream, so count the reaches immediately
          # upstream to decide what to do
          upstream_segs <- nhd_lines %>% 
            filter(tonode %in% site_sf$nearest_segment$fromnode)
          if(nrow(upstream_segs) == 0) {
            # The current reach is a headwater, so use the downstream point
            # and this reach
            site$comid <- site_sf$nearest_segment$comid
            site$fish_dist_to_outlet_m <- fish_dist_downstream_m
            } else if(nrow(upstream_segs) == 1) {
              # The upstream reach exists and does not fork, so match to the
              # upstream point and the reach it drains
              site$comid <- upstream_segs$comid
              # negative distance to indicate upstream:
              site$fish_dist_to_outlet_m <- -fish_dist_upstream_m 
              } else if(nrow(upstream_segs) > 1) {
                # The upstream reach does fork, so just stick to matching to the current reach
                site$comid <- site_sf$nearest_segment$comid
                site$fish_dist_to_outlet_m <- fish_dist_downstream_m
              }
        }
      
      # Regardless of whether we return the current or upstream reach as the
      # best match, the as-a-bird-flies distance to the river should still be
      # the distance to the initial nearest reach, so that we can use this
      # measure to decide whether we were able to snap the site to the river
      # network successfully
      site$bird_dist_to_subseg_m <- site_sf$bird_dist_to_seg_m
      
      return(site)
    })

return(crosswalk)
  
}
  
      

    