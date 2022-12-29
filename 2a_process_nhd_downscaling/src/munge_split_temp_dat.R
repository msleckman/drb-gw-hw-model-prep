#' @title Clean and summarize water temperature data
#'
#' @description 
#' This function reduces duplicate temperature values among rows that share
#' the same site_id and date, and summarizes multiple temperature values along
#' each NHDPlusv2 COMID. This function comes from USGS-R/delaware-model-prep
#' and has been lightly modified here to return one temperature value per COMID-date.
#' https://github.com/USGS-R/delaware-model-prep/blob/main/2_observations.yml
#'
#' @param df a data frame of data. Must contain columns "date", "site_id", 
#' "n_obs", and "sub_location".
#' @prioritize_nwis_sites logical; indicates whether segment summaries should
#' prioritize NWIS observations when available. Defaults to TRUE. If TRUE, any
#' observations from non-NWIS sites will be omitted from the summarized data. 
#'
#' @returns 
#' Returns a data frame with one row per COMID-date.
#'
munge_split_temp_dat <- function(data, prioritize_nwis_sites = TRUE){
  
  # first handle multiple sub_locations that are causing >1 obs per site_id-date.
  # slice_max on n_sub selects the sub-location with the most observations
  sub_location_res <- data %>%
    filter(!grepl('piezometer', sub_location, ignore.case = TRUE)) %>%
    filter(!is.na(sub_location)) %>%
    group_by(site_id, sub_location) %>%
    mutate(n_sub = n()) %>% 
    ungroup() %>%
    group_by(site_id, date) %>%
    slice_max(order_by = n_sub, n = 1, with_ties = FALSE) %>% 
    ungroup()
  
  # bind back with data without sub_location data
  # resolve remaining site_ids with multiple obs
  data2 <- bind_rows(sub_location_res,
                     filter(data, is.na(sub_location)))
  
  # resolve remaining site_ids with >1 obs per date
  # take the site_id with the most n_obs
  dat_dup_resolved <- data2 %>%
    group_by(site_id, date) %>%
    slice_max(order_by = n_obs, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # The munge_split_temp_dat function from USGS-R/delaware-model-prep goes on
  # to summarize temperature by flowline segment. We want to retain the 
  # unaggregated temperature sites so we will return dat_dup_resolved.
  dat_by_seg <- dat_dup_resolved %>%
    group_by(comid, date) %>%
    # If prioritize_nwis_sites is TRUE, check whether data for that segment
    # comes from multiple sources. If multiple distinct site id's, retain only
    # NWIS sites for that segment-date; otherwise, retain all samples
    {if(prioritize_nwis_sites){
      filter(., if(n_distinct(site_id) > 1 & any(grepl("nwis", source, ignore.case = TRUE))) grepl("nwis", source, ignore.case = TRUE) else TRUE)
    } else {.}
    } %>%
    summarize(site_id = paste0(site_id, collapse = ', '),
              source = paste0(unique(source), collapse = ', '),
              time = ifelse(n() > 1, NA, time), # keep timestamp if represents a single value
              mean_temp_c = round(mean(mean_temp_degC), 1),
              min_temp_c = min(min_temp_degC),
              max_temp_c = max(max_temp_degC),
              sd_mean_temp_c = round(sd(mean_temp_degC), 1), # provide an indicator of variability across sites
              flag = paste(unique(flag)[!is.na(unique(flag))], collapse = '; '),
              .groups = 'drop') %>%
    rowwise() %>%
    mutate(flag = paste(unique(unlist(strsplit(flag, '; '))), collapse = '; ')) 
  
  return(dat_by_seg)
  
}
