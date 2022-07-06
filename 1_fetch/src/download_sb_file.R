download_sb_file <- function(sb_id, out_dir, file_name = NULL, overwrite_file = TRUE){
  #'
  #' @description Function to download file from ScienceBase 
  #'
  #' @param sb_id string - the id of the science base item
  #' @param file_name string - the name of the file in the science base item to download
  #' @param out_dir string - the directory where you want the file downloaded to
  #'
  #' @value string the out_path
  
  if(!is.null(file_name)){
    
    out_path = file.path(out_dir, file_name)
    # Get the data from ScienceBase:
    sbtools::item_file_download(sb_id = sb_id,
                                names = file_name,
                                destinations = out_path,
                                overwrite_file = overwrite_file)
  } else {
    # Get the data from ScienceBase:
    out_path <- sbtools::item_file_download(sb_id = sb_id,
                                            dest_dir = out_dir, 
                                            overwrite_file = overwrite_file)
    }
  
  return(out_path)
}
