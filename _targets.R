library(targets)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse","sbtools","sf","nhdplusTools", 'purrr')) 

# dir for selected datasets soil characteristics from nhd statsgo
dir.create('1_fetch/out/statsgo', showWarnings = FALSE)

source("1_fetch.R")

## nhd parent id 
# pulled from https://www.sciencebase.gov/catalog/item/5728d6ace4b0b13d3918a992
nhd_statsgo_parent_sbid <- '5728d6ace4b0b13d3918a992'

# Return the complete list of targets
c(p1_targets_list)


