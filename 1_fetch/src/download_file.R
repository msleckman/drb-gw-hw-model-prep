#'
#' @title Download file from url
#' 
#' @description 
#' Download a file from the internet and return a character string indicating
#' the file path.
#' 
#' @param url character string indicating the URL of the file to be downloaded
#' @param fileout character string with the name and file path where the downloaded file
#' should be saved
#' @param mode character string indicating the mode used to write the file; see 
#' ??utils::download.file for details
#' @param quiet logical (defaults to FALSE); if TRUE, suppresses any status messages
#' 
download_file <- function(url,fileout,mode,quiet = FALSE){

  
  # download the file
  utils::download.file(url, destfile = fileout, mode = mode, quiet = quiet)

  # return the file name including file path and extension
  return(fileout)
  
}


#' @title Download file from ScienceBase
#'
#' @description 
#' Function to download files from ScienceBase and return a character string 
#' indicating the file path.
#'
#' @param sb_id string - the id of the science base item
#' @param file_name string - the name of the file in the science base item to download
#' @param out_dir string - the directory where you want the file downloaded to
#' @param overwrite_file whether to overwrite download if file already present in out_dir. Default TRUE. see sbtools::item_file_download() for further details 
#'
#' @value string the out_path
#' 
download_sb_file <- function(sb_id, out_dir, file_name = NULL, overwrite_file = TRUE){

  
  if(!is.null(file_name)){
    out_path = file.path(out_dir, file_name)
    tryCatch(
      # Get the data from ScienceBase:
      {sbtools::item_file_download(sb_id = sb_id,
                                   names = file_name,
                                   destinations = out_path,
                                   overwrite_file = overwrite_file)},
      # Catching error if file_name does not output from item_file_download
      error = function(e){
        sb_str_format_url <- 'https://www.sciencebase.gov/catalogMaps/mapping/ows/%s?service=wfs&request=GetFeature&typeName=sb:%s&outputFormat=shape-zip&version=1.0.0'
        ## sb url with facet defined - optional alt. commenting out for now. Consider using at a later point
        # sb_str_format_url_fct <- 'https://www.sciencebase.gov/catalog/file/get/%s?facet=%s'
        
        # if else depending on file name. Query cannot have *zip ext.
        query <- ifelse(
          grepl('.zip',file_name),
          sprintf(sb_str_format_url,
                  sb_id,
                  substr(file_name, start = 1, stop = nchar(file_name))),
          sprintf(sb_str_format_url, sb_id, file_name))
        
        message(paste('File not directly accessible with input sb_id using sbtools.\nDownloading', file_name, 'directly from:', query))
        
        httr::GET(query,httr::write_disk(out_path, overwrite=TRUE))
      }
    )
    
  } else {
    # Get the data from ScienceBase:
    out_path <- sbtools::item_file_download(sb_id = sb_id,
                                            dest_dir = out_dir, 
                                            overwrite_file = overwrite_file)
  }
  
  return(out_path)
  
}


