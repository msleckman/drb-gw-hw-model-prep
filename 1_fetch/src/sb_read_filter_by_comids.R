#' Function to filter sb datasets to given comids in AOI. Generates unzip if needed.  

sb_read_filter_by_comids <- function(data_path = '1_fetch/out/statsgo', comid, cbind = TRUE){
  
  #' @param comid character string or vector of character strings containing
  #' the common identifier (COMID) of the desired flowline(s)
  #' @param data_path folder path to files downloaded from sb
  #' @param cbind if true read in files are merged by column COMID
  
  ##unzip  
  data_zip_lst <- list.files(path = data_path, pattern = '*.zip', full.names = TRUE)
  lapply(data_zip_lst, function(x) unzip(x, exdir = '1_fetch/out/statsgo',overwrite = FALSE))
  
  ## identify unzipped files
  data_lst <- list.files(path = data_path, pattern = '*.TXT|*.csv', full.names = TRUE)
  
  ## read in + clip
  data <- lapply(data_lst, function(x) read.delim(x, sep = ',') %>% 
                   filter(COMID %in% comid)
                 )
  ## cbind if true
  if(cbind == TRUE){
    full_df <- data %>% reduce(full_join, by = 'COMID')
  } else{full_df <- data}
  
  return(full_df)  
} 
