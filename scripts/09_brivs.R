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


# Clean up data
briv_scrub_dat <- briv_dat |> 
  select(52, everything ()) |> 
  filter(notes != "No visibility") |> 
  filter(notes != "video data not usable") |> 
  filter(notes != "No data, camera did not record") |> 
  filter(notes != "No Visibility") |> 
  filter(notes != "No visibility for first 2 videos, count started after bait was eaten") |>
  filter(notes != "No Visibility, cannot see bait for much of video, benthos cannot be seen at all. CAMA at 1098 was persistent and potentially consumer, however visibility prevents conclusive determination.") |>
  filter(notes != "No Visibility, cannot see bait or benthos") |>
  filter(notes != "No Visibility, cannot see bait or benthos; bait quickly targeted by CAMA, which dissapears and reappears frequently due to poor visibility. Likely only one individual, and this islikely the consumer of bait. Impossible to confirm the moment that the bait is actualy removed entirely.") |> 
  filter(notes != "No visibility. A crab appears at somepoint clasped to the BRIV and feeding on squid. Unknown when it first appears, whether it eats the bait completly, and what it is.") |> 
  filter(notes != "No visibility for much of video.Unknown what ultimatly consumes the bait. Likely CAMA but cannot be confirmed from the video")
  
briv_scrub_dat |> 
  group_by(site, treatment) |> 
  summarise(n = n())
# When scrubbed of poor quality video, no replication remains for predation analysis
# will just look at some descriptive stats for this

briv_dat |> 
  group_by(site, treatment, first_predator) |>
  filter(!is.na(first_predator)) |> 
  filter(first_predator != "NOSP") |> 
  summarise(n = n())
  
briv_dat |> 
  group_by(site, treatment, bait_consumer) |>
  filter(!is.na(bait_consumer)) |> 
  filter(bait_consumer != "NOSP") |> 
  summarise(n = n())

# For usable videos, save for two instances where an unknown crab could not be identified, the first predator was always CAMA, across all sites and treatments.
# Likewise, bait consumer was always CAMA.