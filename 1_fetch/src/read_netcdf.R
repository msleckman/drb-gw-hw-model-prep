#' Read netcdf file and format as an R data frame
#' 
#' @details This function contains three different methods for
#' reading in and formatting netcdf files. Two methods (using 
#' functions from raster and ncdf4) are commented out and we're 
#' currently using tidync::hyper_tibble() to format the input/output
#' data as a data frame.
#' 
#' @param nc_file character string indicating the netcdf file name, 
#' including file path and .nc extension
#' 
read_netcdf <- function(nc_file){
  
  # Find variables and get netcdf attributes
  nc <- ncdf4::nc_open(nc_file)
  variables <- names(nc[['var']])
  seg_id_nat <- ncdf4::ncvar_get(nc, "seg_id_nat")
  tunits <- ncdf4::ncatt_get(nc,"time","units")
  start_date <- sub( "^\\D+", "", tunits$value)
  time <- ncdf4::ncvar_get(nc, "time")
  
  # 1) Read in individual variables from netcdf using {raster}
  # [Lauren] this implementation is clean and seems to work, but there are 
  # a couple of issues that I don't know how to address right now: 1) xy is 
  # set up for lat/lon data and so is reading in seg_id_nat as numeric and 
  # reformatting all of the values into floats. I haven't been able to figure
  # out a workaround for that; 2) more importantly, this raster implementation
  # returns an output that differs from the ncdf4 output(!). Comparing static
  # attribute values from the 2021 forecasting data release leads me to believe
  # the ncdf4 code returns the correct values, so what's going on with raster?
  # Leaving this code here but commented out in case we want to return to this.
  #ncvars_ls <- lapply(variables, function(i){
  #  r <- raster::raster(nc_file,varname = i, stopIfNotEqualSpaced = FALSE)
  #  r_df <- raster::as.data.frame(r, xy = TRUE) %>%
  #    rename(time = x, seg_id_nat = y, !!i := 3) %>%
  #    mutate(date = as.Date(as.numeric(time), origin = start_date))
  #})
  #names(ncvars_ls) <- variables
  
  # 2) Read in individual variables from netcdf using {ncdf4}
  #ncvars_ls <- lapply(variables, function(i){
  #  # load variable data
  #  r <- ncdf4::ncvar_get(nc, i)
  #  
  #  # format as a data frame
  #  r_df <- as.data.frame(r)
  #  colnames(r_df) <- seg_id_nat
  #  r_df$time <- seq(time)-1 # seq starts at 1 instead of zero, so subtract 1
  #  
  #  # reformat data frame from wide to long 
  #  r_df_out <- r_df %>%
  #    pivot_longer(!time, names_to = 'seg_id_nat', values_to = i) %>%
  #    mutate(date = as.Date(time, origin = start_date)) 
  #}) 
  
  # convert list into a data frame
  #ncvars_df <- purrr::reduce(ncvars_ls, full_join, by = c('date','time','seg_id_nat')) %>%
  #  arrange(seg_id_nat, time)
  
  # close netcdf file
  ncdf4::nc_close(nc)
  
  # 3) Use tidync package instead: https://ropensci.org/blog/2019/11/05/tidync/
  ncvars_df <- nc_file %>% 
    tidync::hyper_tibble() %>% 
    mutate(date = as.Date(time, origin = start_date)) %>%
    relocate(seg_id_nat, .before = seg_rain) %>% 
    relocate(time, .after = seg_id_nat) %>%
    relocate(date, .after = time)
  
  return(ncvars_df)
  
}

