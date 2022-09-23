#' @title Process NHDPlusv2 attributes for cumulative upstream watershed
#' 
#' @description 
#' Function to read in downloaded NHDv2 attribute data and join with river segment ID's.
#' 
#' @details This function was pulled and modified from the inland salinity ml project:
#' https://github.com/USGS-R/drb-inland-salinity-ml/blob/main/2_process/src/process_nhdv2_attr.R
#'
#' @param file_path file path of downloaded NHDv2 attribute data table, including file extension
#' @param segs_w_comids data frame containing the PRMS segment ids and the comids of interest
#' segs_w_comids must contain variables PRMS_segid and COMID
#' @param cols character string indicating which columns to retain from downloaded attribute data; 
#' cols can take values "ACC" or "TOT"
#'
#' @value A data frame containing PRMS_id and columns representing the NHDv2 attribute data referenced to the 
#' cumulative upstream watershed.
#' 
process_cumulative_nhdv2_attr <- function(file_path,segs_w_comids,cols){

  message(file_path)
  # Read in downloaded data 
  # only specify col_type for COMID since cols will differ for each downloaded data file
  dat <- read_csv(file_path, col_types = cols(COMID = "c"), show_col_types = FALSE)
  

  # LAUREN: Commenting out the two lines below that are used in inland salinity but
  # not for our groundwater data prep.
  # For PPT data we want to return the long-term (1971-2000) monthly averages 
  # instead of the monthly values for each year
  #if(grepl("PPT_TOT",file_path)|grepl("PPT_ACC",file_path)){
  #  message("Calculating long-term monthly average precipitation from annual data")
  #  dat <- calc_monthly_avg_ppt(dat)
  #}
  
  # LAUREN: Commenting out the two lines below that are used in inland salinity
  # but not for our groundwater data prep.
  # For NADP data we want to return the long-term (1984-2014) average 
  # instead of annual values for each constituent
  #if(grepl("NADP",file_path)){
  #  message("Calculating long-term average NADP from annual data")
  #  dat <- calc_avg_NADP(dat)
  #}
    
  # Process downloaded data
  dat_proc <- dat %>%
    # retain desired columns ('ACC' or 'TOT')
    select(c(COMID, starts_with(cols))) %>%
    # join data to {segs_w_comids} data frame by COMID
    right_join(., segs_w_comids,by = "COMID") %>%
    relocate("PRMS_segid", .before = "COMID") %>%
    relocate("seg_id_nat", .before = "COMID") %>%
    select(-COMID)
  
  # Flag columns with undesired flag values (e.g. -9999)
  flag_cols <- dat_proc %>%
    select(where(function(x) -9999 %in% x)) %>% 
    names()
  
  # For columns with undesired flag values, replace -9999 with NA, else use existing value
  dat_proc_out <- dat_proc %>%
    mutate(across(all_of(flag_cols), 
                  ~case_when(. == -9999 ~ NA_real_, 
                             TRUE ~ as.numeric(.))))

  return(dat_proc_out)
  
}



