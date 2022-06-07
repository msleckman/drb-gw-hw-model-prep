library(targets)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse","sbtools","sf","nhdplusTools")) 

source("1_fetch.R")

# Return the complete list of targets
c(p1_targets_list)


