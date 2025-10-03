#' -------------------------------------
#' Look at total species richness
#' across berms versus controls and
#' quads (all), seines, traps, and brivs combined
#' -------------------------------------
source("scripts/helpers.R")

## load libraries
library(tidyverse)


## load data, bind them together, and add taxonomic data
quad_dat <- read_csv("data/quads.csv")
seine_dat <- read_csv("data/seine.csv")
trap_dat <- read_csv("data/traps.csv")
farm_dat <- read_csv("data/farms.csv")
briv_dat <- read_csv("data/brivs.csv")
msl <- read_csv("data/master_sp_list.csv")

# Drop uneeded cols in msl and rename col
msl <- select(msl, 1:10) |> 
  rename(species_code = code)

# bind data together, and drop the two NAs in Bayside Berm
richness_dat <- rbind(quad_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)), 
                      seine_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)),
                      trap_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)),
                      briv_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(first_predator)),
                      briv_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(bait_consumer))
) |> 
  filter(species_code != is.na(species_code)) |> 
  filter(species_code != "NOSP")
  

# attach taxonomic data
richness_dat <- left_join(richness_dat, 
                          msl, 
                          "species_code")

## Write a function that references the MSL, checks if a given site-treatment 
## has a lower order (i.e. more specific) spp entry, if yes then drop the
## higher order code and dont count it in richness, if no then keep and count, 
## function should spit out a summarised df with richness counts by level

## Calculate various richnesses

kingdom_richness_dat <- richness_dat |> 
  group_by(site, treatment) |>
  select(site, treatment, kingdom) |> 
  unique() |> 
  reframe(kingdom_r = n())
# Additional group is from VERR

# phylum level richness
phylum_richness_dat <- richness_dat |> 
  select(site, treatment, phylum) |> 
  group_by(site, treatment) |>
  unique() |> 
  reframe(richness = n())

richness_bylevel_dat <- 


### These lower levels need to keep or drop higher levels with no finer level assigned based on whether there is already one (see above notes on a function) 
# class level richness
# order level richness
# family level richness
# genus level richness
# species level richness



## Initial Viz
richness_plot <- ggplot(data = richness_dat,
       mapping = aes(x = treatment,
                     y = richness,
                     fill = treatment)) +
  geom_boxplot(show.legend = FALSE) +
  ylab("Species Richness") +
  xlab("Treatment") +
  scale_fill_manual(values = c("Berm" = "#af8dc3", "Control" = "#7fbf7b"))

# write out richness figure
ggsave("figures/richness_plot.jpg", 
       richness_plot,
       dpi = 600)

## Analysis of Richness


## The plan
# sumarise and calculate richness by site and taxonomic level 
# (should be one df w/ header:site, phy rich, cl rich, etc

## Then ttest to determine if the richnesses are different by treatment and at what levels