coarse_sediment_area_calc <- function(buffered_reaches_sf, 
                                 coarse_sediments_area_sf,
                                 prms_col = 'PRMS_segid'){
  #' @description function takes the a vector polygon sf object describing areas of coarse stratified sediment surface material
  #' and calculates the proportion of this area within the NHM reach buffer zones (250 m).
  #' @param buffered_reaches_sf
  #' @param coarse_sediments_area_sf
  #' @param prms_col
  #' @value The output is a sf object of the buffered reaches and three calculated columns: 
  #' total_reach_buffer_area_km2 = total area of buffered reach
  #' cs_area_km2 = area of coarse stratified sediment within given buffered reach  
  #' cs_area_proportion = proportion of coarse stratified sediment area within reach area 
  
  # Grab area of coarse sediment within buffer
  cs_df <- sf::st_intersection(buffered_reaches_sf,
                               coarse_sediments_area_sf) %>% 
    group_by(.data[[prms_col]]) %>% 
    dplyr::summarize(geometry = sf::st_union(geometry)) %>% 
    mutate(cs_area_km2 = units::set_units(st_area(.), km^2)) %>% 
    ungroup() %>% 
    st_drop_geometry()
 
  # Join to final buffered reaches df 
  buffered_reaches_w_cs_df <- left_join(buffered_reaches_sf, 
                             cs_df, 
                             by = prms_col) %>%
  ## creating area cols for final df with all PRMS reaches  
  mutate(cs_area_km2 = units::set_units(ifelse(is.na(cs_area_km2), 0, cs_area_km2), km^2),
           total_reach_buffer_area_km2 = units::set_units(st_area(.), km^2),
          ## calculating proportion of coarse stratified sediment within buffered reach 
           cs_area_proportion = round(cs_area_km2 / total_reach_buffer_area_km2, 3)) %>%
    ## re order cols 
    relocate(geometry, .after = cs_area_proportion)
  
  ## removing unit on proportion
  buffered_reaches_w_cs_df$cs_area_proportion <- units::drop_units(buffered_reaches_w_cs_df$cs_area_proportion)
  
  return(buffered_reaches_w_cs_df)
  
}


