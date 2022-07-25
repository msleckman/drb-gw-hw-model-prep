library(targets)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse","sbtools","sf","nhdplusTools", 'purrr', 'terra','raster', "arrow", "tidync", "ncdf4", "reticulate")) 

# dir for selected datasets soil characteristics from nhd statsgo
dir.create('1_fetch/out/statsgo', showWarnings = FALSE)

source("1_fetch.R")
source("2_process.R")

# Define the url for the NHGFv1 to NHDv2 crosswalk 
# Note that there are different versions of the crosswalk table that might be best suited
# to different use cases (e.g. retain divergences? retain zero-area flowlines?)
# see https://github.com/USGS-R/drb-network-prep/commit/3637931f5a17469a4234eaed3d20ed44ba45958d
GFv1_NHDv2_xwalk_url <- paste0("https://raw.githubusercontent.com/USGS-R/drb-network-prep/",
                               "3637931f5a17469a4234eaed3d20ed44ba45958d",
                               "/2_process/out/GFv1_NHDv2_xwalk.csv")


## nhd parent id 
# pulled from https://www.sciencebase.gov/catalog/item/5728d6ace4b0b13d3918a992
nhd_statsgo_parent_sbid <- '5728d6ace4b0b13d3918a992'

drb_huc8s <- c("02040101","02040102","02040103","02040104","02040105","02040106",
               "02040201","02040202","02040203","02040204","02040205","02040206","02040207")

# Return the complete list of targets
c(p1_targets_list, p2_targets_list)


