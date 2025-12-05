#' -------------------------------------
#' Look at BRIV predation success and rate
#' across berms versus controls and
#' -------------------------------------

source("scripts/helpers.R")


##
# Set-up
##

#load libraries
library(tidyverse)
library(broom)

# load data
briv_dat <- read_csv("data/brivs.csv")

## Verify same amount of brivs per site-treatment
# briv_dat |>
#   group_by(site, treatment) |>
#   reframe(briv_count = n()) |>
#   View()
# Not quite even, Coughlin has 5 in each, bayside is 9/8, dux is 10/10


## Next steps: track down inconsistencies in the brivs data, where is notes column and why are there all these NAs!