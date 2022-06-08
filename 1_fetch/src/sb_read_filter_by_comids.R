#' Function to filter sb datasets to given comids in AOI. Generates unzip if needed.  

sb_read_filter_by_comids <- function(data_path = '1_fetch/out/statsgo', sb_comid_col = 'COMID', comid, cbind = TRUE){
  
  #' @param comid character string or vector of character strings containing
  #' the common identifier (COMID) of the desired flowline(s)
  #' @param data_path folder path to files downloaded from sb
  #' @param cbind if true read in files are merged by column COMID
  
  ##unzip - if no zip this will 
  data_zip_lst <- list.files(path = data_path, pattern = '*.zip', full.names = TRUE)
  
  if(!is_empty(data_zip_lst)){
  ##unzipping in sale location as zipped. No subfolder created 
  lapply(data_zip_lst, function(x) unzip(x, exdir = data_path))
  }

  ## Pull downloaded files
  data_lst <- list.files(path = data_path, pattern = '*.TXT|*.csv', full.names = TRUE)
  
  ## read in + clip
  data <- lapply(data_lst[1:4], function(x) read.delim(x, sep = ',') %>% 
                   filter(.data[[sb_comid_col]] %in% comid))
                 
  ## cbind if true
  if(cbind == TRUE){
    full_df <- data %>% reduce(full_join, by = sb_comid_col)
  } else{full_df <- data}
  
  return(full_df)  
} 
