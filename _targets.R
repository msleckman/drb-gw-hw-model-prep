library(targets)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse","sbtools","sf","nhdplusTools", 'purrr')) 

# dir for selected datasets from nhd statsgo
dir.create('1_fetch/out/statsgo', showWarnings = FALSE)

source("1_fetch.R")

# Define minor HUCs (hydrologic unit codes) that make up the DRB to use in calls to dataRetrieval functions
# Lower Delaware: 0204 subregion (for now, exclude New Jersey Coastal (https://water.usgs.gov/GIS/huc_name.html)
drb_huc8s <- c("02040101","02040102","02040103","02040104","02040105","02040106",
               "02040201","02040202","02040203","02040204","02040205",
               "02040206","02040207")

## nhd parent id 
# pulled from https://www.sciencebase.gov/catalog/item/5728d6ace4b0b13d3918a992
nhd_statsgo_parent_sbid <- '5728d6ace4b0b13d3918a992'

# Return the complete list of targets
c(p1_targets_list)


