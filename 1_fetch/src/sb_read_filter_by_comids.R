#' Function to filter sb datasets to given comids in AOI. Generates unzip if needed.  

sb_read_filter_by_comids <- function(comid, data_path, sb_comid_col = 'COMID',
                                     selected_cols_contains = NULL, cbind = TRUE){
  
  #' @description function to unzip, fetch and munge datasets from Sciencebase. Built for STATSGO but generalized for ther SB datasets
  #' @param comid character string or vector of character strings containing the common identifier (COMID) of the desired flowline(s)
  #' @param data_path folder path to files downloaded from sb
  #' @param sb_comid_col str of the the column name from the comids. Default is 'COMID' as this is typically the name of this column from ds from SB0
  #' @param selected_cols_contains char vector of columns of interest from original dataset e.g. c('KFACT','KFACT_UP','NO10AVE','NO4AVE','SILTAVE','CLAYAVE','SANDAVE'). Do not include the comid id colname
  #' @param cbind if true read in files are merged by column COMID
  
  # 1) unzip - if no zip this will pass 
  data_zip_lst <- list.files(path = data_path, pattern = '*.zip', full.names = TRUE)
  
  if(!is_empty(data_zip_lst)){
  ##unzipping in same location as zipped. No subfolder created 
  lapply(data_zip_lst, function(x) unzip(x, exdir = data_path, overwrite = TRUE))
  }
  
  # 2) read
  ## Pull downloaded files
  data_lst <- list.files(path = data_path, pattern = '*.TXT|*.csv', full.names = TRUE)
  if(is_empty(data_lst)){
    stop('No txt or csv file in this data path')
  }
  
  ## read in + clip
  data <- lapply(data_lst, function(x) read.csv(x, sep = ',') %>% 
                   ### selecting only specified columns of param selected_cols_contains
                   {if(!is.null(selected_cols_contains))
                     select(.,contains(append(sb_comid_col,
                                              selected_cols_contains))) else .} %>% 
                   ### filter to comid in AOI
                   filter(.data[[sb_comid_col]] %in% comid))
  
  # 3) cbind if true
  if(cbind == TRUE){
    full_df <- data %>% reduce(full_join, by = sb_comid_col)
  } else{full_df <- data}
  
  return(full_df)  
} 
