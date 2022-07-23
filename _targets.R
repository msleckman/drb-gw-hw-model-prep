library(targets)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse","sbtools","sf","nhdplusTools", 'purrr', 'terra','raster')) 

# dir for selected datasets soil characteristics from nhd statsgo
dir.create('1_fetch/out/statsgo', showWarnings = FALSE)

source("1_fetch.R")
source("2_process.R")

# Define the url for the NHGFv1 to NHDv2 crosswalk 
# see https://github.com/USGS-R/drb-network-prep/commit/614177cc96ed955db364dc0a1dc7d03adffbd33b
GFv1_NHDv2_xwalk_url <- paste0("https://raw.githubusercontent.com/USGS-R/drb-network-prep/",
                               "614177cc96ed955db364dc0a1dc7d03adffbd33b",
                               "/2_process/out/GFv1_NHDv2_xwalk.csv")


## nhd parent id 
# pulled from https://www.sciencebase.gov/catalog/item/5728d6ace4b0b13d3918a992
nhd_statsgo_parent_sbid <- '5728d6ace4b0b13d3918a992'

drb_huc8s <- c("02040101","02040102","02040103","02040104","02040105","02040106",
               "02040201","02040202","02040203","02040204","02040205","02040206","02040207")

# Return the complete list of targets
c(p1_targets_list, p2_targets_list)


