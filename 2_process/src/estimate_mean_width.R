#' Function to estimate mean width for an NHDPlusv2 reach
#'
#' @param nhd_lines sf linestring containing the NHDPlusv2 flowlines. Must
#' contain column COMID indicating the unique reach identifier
#' @param buffer_dist_m integer; look for NWIS sites within this distance
#' in meters. Defaults to 500 m.
#' @param estimation_method character string; what method should be used 
#' to estimate NHDflowline widths? valid options include "nwis".
#' @param network_pos_variable character string; what NHDPlusv2 attribute
#' should be used to build an empirical regression that predicts mean width?
#' Defaults to "arbolate_sum". Other valid entries are "upstream_area" and
#' "ma_flow", which is the gage-adjusted mean annual flow for a given COMID
#' as indicated by the value-added attribute QE_MA.
#' 
estimate_mean_width <- function(nhd_lines, 
                                buffer_dist_m = 500, 
                                estimation_method = 'nwis', 
                                network_pos_variable = 'arbolate_sum'){
  
  # Check that user entries match expected inputs
  stopifnot("Value provided for estimation_method is not a valid entry" =
              estimation_method %in% c("nwis"))
  stopifnot("Value provided for network_pos_variable is not a valid entry" =
              network_pos_variable %in% c("arbolate_sum","upstream_area","ma_flow"))
  
  # 1) One approach is to estimate width using an empirical relationship
  # between measured width and upstream area that is developed using
  # USGS field measurements of width. 
  if(estimation_method == 'nwis'){
    
    message("Downloading width data from NWIS...")
    
    # return a list of nwis sites within the bounding box of nhd_lines
    nwis_sites <- dataRetrieval::whatNWISsites(bBox = sf::st_bbox(nhd_lines),
                                               parameterCd = "00060")
    
    nwis_sites_sf <- nwis_sites %>% 
      sf::st_as_sf(coords = c('dec_long_va','dec_lat_va'), crs = 4269) %>%
      sf::st_transform(st_crs(nhd_lines))
    
    # The bounding box approach above returns more sites than we're
    # interested in, omit sites that aren't within some distance of
    # the given nhd network 
    nwis_sites_along_nhd <- nwis_sites_sf %>%
      sf::st_filter(y = nhd_lines, 
                    .predicate = st_is_within_distance,
                    dist = units::set_units(buffer_dist_m, m))

    # Join nwis sites to nearest nhd flowline
    nwis_sites_snapped <- nwis_sites_along_nhd %>%
      sf::st_transform(st_crs(nhd_lines)) %>%
      sf::st_join(y = nhd_lines, join = st_nearest_feature)
    
    # Download width measurements from NWIS
    nwis_widths <- nwis_sites_snapped %>%
      sf::st_drop_geometry() %>%
      split(.,.$site_no) %>%
      lapply(., function(x){
        tryCatch(summarize_nwis_widths(x$site_no), 
                 error = function(cond){
                   return(NULL)
                 },
                 warning = function(cond){
                   return(NULL)
                 })
      }) %>%
      bind_rows()
    
    # create empirical relationship between NWIS widths and arbolate sum, 
    # as in this USGS report, https://pubs.usgs.gov/sir/2021/5116/sir20215116.pdf,
    # or by total upstream area, or by the gage-adjusted estimated mean annual 
    # discharge (NHDPlus attribute qe_ma)
    message(sprintf(
      paste0("Estimating width for NHDPlusv2 flowlines based on empirical fit ",
             "between NWIS measurements and attribute %s..."), network_pos_variable))
    
    nwis_sites_w_width <- nwis_widths %>%
      left_join(nwis_sites_snapped, by = "site_no")
      
    fit_arbsum <- summary(lm(log10(width_m) ~ log10(arbolatesu), data = nwis_sites_w_width))
    fit_totda <- summary(lm(log10(width_m) ~ log10(totdasqkm), data = nwis_sites_w_width))
    fit_qe_ma <- summary(lm(log10(width_m) ~ log10(qe_ma), data = nwis_sites_w_width))
    
    width_fits <- tibble(network_pos_variable = c('arbolate_sum','upstream_area','ma_flow'),
                         slope = c(fit_arbsum$coefficients[2],
                                   fit_totda$coefficients[2],
                                   fit_qe_ma$coefficients[2]),
                         intercept = c(fit_arbsum$coefficients[1],
                                       fit_totda$coefficients[1],
                                       fit_qe_ma$coefficients[1]),
                         r.squared = c(fit_arbsum$r.squared, 
                                       fit_totda$r.squared, 
                                       fit_qe_ma$r.squared))
    
    if(network_pos_variable == "arbolate_sum"){
      slope_arbsum <- width_fits$slope[width_fits$network_pos_variable == "arbolate_sum"]
      int_arbsum <- width_fits$intercept[width_fits$network_pos_variable == "arbolate_sum"]
      
      nhd_lines_out <- nhd_lines %>%
        mutate(est_mean_width_m = arbolatesu * slope_arbsum + int_arbsum)
    }
    if(network_pos_variable == "upstream_area"){
      slope_totda <- width_fits$slope[width_fits$network_pos_variable == "upstream_area"]
      int_totda <- width_fits$intercept[width_fits$network_pos_variable == "upstream_area"]
      
      nhd_lines_out <- nhd_lines %>%
        mutate(est_mean_width_m = totdasqkm * slope_totda + int_totda)
    }
    if(network_pos_variable == "ma_flow"){
      slope_qe_ma <- width_fits$slope[width_fits$network_pos_variable == "ma_flow"]
      int_qe_ma <- width_fits$intercept[width_fits$network_pos_variable == "ma_flow"]
      
      nhd_lines_out <- nhd_lines %>%
        mutate(est_mean_width_m = qe_ma * slope_qe_ma + int_qe_ma)
    }
  }    
    
    return(nhd_lines_out)

}


#' Function to download USGS in situ field measurements of river width
#' 
#' @param site_no character string indicating the USGS site number
#' 
#' @examples download_nwis_widths("01481000")
#' 
download_nwis_widths <- function(site_no){
  
  # download width measurements from NWIS
  width_data <- dataRetrieval::readNWISmeas(siteNumbers = site_no, expanded = TRUE)
    
  # Format in situ width measurements
  width_data_out <- width_data %>%
    filter(!is.na(chan_width)) %>%
    mutate(width_m = chan_width * 0.3048) 
  
  return(width_data_out)
  
}



#' Function to download USGS in situ field measurements of river width
#' and return the median width in meters.
#' 
#' @param site_no character string indicating the USGS site number
#' @param earliest_date character string formatted as "YYYY-MM-DD", will 
#' return width measurements made on or after this date. Defaults to NULL, 
#' which returns all available width measurements.
#' @param latest_date character string formatted as "YYYY-MM-DD", will 
#' return width measurements made on or before this date. Defaults to NULL,
#' which returns all available width measurements.
#' 
#' @examples summarize_nwis_widths("01481000", earliest_date = "1980-10-01")
#'
summarize_nwis_widths <- function(site_no, earliest_date = NULL, latest_date = NULL){
  
  # Download width data from NWIS
  width_data <- download_nwis_widths(site_no)
  
  # Summarize in situ field widths
  width_data_out <- width_data %>%
    # For days with multiple width measurements, summarize
    # data into a single width value for that date:
    group_by(measurement_dt) %>%
    summarize(site_no = unique(site_no),
              width_median_dt = median(width_m),
              .groups = 'drop') %>%
    # If `earliest_date` is not NULL, return only width measurements
    # made on or after `earliest_date`:
    {
      if(!is.null(earliest_date)){
        filter(.,measurement_dt >= earliest_date)
      } else {.}
    } %>%
    # If 'latest_date' is not NULL, return only width measurements
    # made on or before `latest_date`:
    {
      if(!is.null(latest_date)){
        filter(.,measurement_dt <= latest_date)
      } else {.}
    } %>%
    summarize(site_no = unique(site_no),
              width_m = median(width_median_dt),
              width_measures_n = n(),
              .groups = 'drop')
  
  return(width_data_out)
  
}

  



