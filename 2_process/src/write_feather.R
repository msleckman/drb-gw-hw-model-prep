#' Function to save data frames as a feather file
#' 
#' @param data data frame object to save as a feather file
#' @param file_out character string indicating the name of the output file,
#' including file path and .feather extension
#' 
write_feather <- function(data, file_out){
  
  # save data as a feather file
  arrow::write_feather(data, file_out)
  
  # return the file name
  return(file_out)
  
}

