#' -------------------------------------
#' Look at functional group cover
#' across berms versus controls and
#' quads cover only
#' -------------------------------------
source("scripts/helpers.R")

## Planning
#' Determine what func groups we are interested in: Habitat formers (ASCO, FUSP, CHCR/MAST), 
#' sessile inverts (SEBA, MYED); Predators, i.e. crabs
#' Pull out the relevant taxa by selecting those codes
#' Compare abundance (pc or count) by quad for each species
#' lm brown_macrophytes ~ site + treatment + zone
#' 

## load libraries
library(tidyverse)
library(broom)
library(vegan)
library(patchwork)

## load data
quad_dat <- read_csv("data/quads.csv")

# Sum across squares for percent data
percent_dat <- quad_dat |> 
  filter(measurement_type != "Count") |> 
  group_by(site, treatment, height, quadrat, measurement_type, species_code) |> 
  reframe(sum = sum(measurement)) |> 
  filter(species_code != "NOSP") |> 
  ungroup()

# Sum across squares for count data
count_dat <- quad_dat |> 
  filter(measurement_type == "Count") |> 
  group_by(site, treatment, height, quadrat, measurement_type, species_code) |> 
  reframe(sum = sum(measurement)) |> 
  filter(species_code != "NOSP") |> 
  ungroup()

##
# Brown Macrophytes
# Analysis of change in ASNO, FUCU, FUSP, FUDI, FUVE
##

phaeo_dat <- percent_dat |> 
  filter(species_code %in% c("ASNO", "FUCU", "FUSP", "FUVE", "FUDI")) |> 
  group_by(site, treatment, height) |> 
  summarise(total_percent = sum(sum))
# Basically only Bayside berm mid and high have any amount of brown macrophytes.

##
# Red Macrophytes
# Analysis of change in CHCR, MAST
##

rhodo_dat <- percent_dat |> 
  filter(species_code %in% c("CHCR", "MAST")) |> 
  group_by(site, treatment, height) |> 
  summarise(total_percent = sum(sum))

# Only found in low zone at bay berm, and cough berm and control. Essentially the same

##
# Sessile Inverts
# Analysis of change in SEBA, MYED
##

invert_dat <- percent_dat |> 
  filter(species_code %in% c("SEBA", "MYED")) |> 
  group_by(site, treatment, height, species_code) |> 
  summarise(total_percent = sum(sum))

seba_plot <- ggplot(percent_dat |> filter(species_code == "SEBA"),
                    mapping = 
                          aes(x = treatment, 
                              y = sum, 
                              fill = treatment)) +
  facet_grid(rows = vars(height), cols = vars(site)) +
  geom_point(alpha = .6, position = "jitter") +
  geom_boxplot(alpha = .6)

# Stats
percent_dat |> 
  filter(species_code == "SEBA") |>
  pivot_wider(names_from = treatment)
  t.test(x, y, paired = TRUE)

## Results - no difference in SEBA btw sites. Better at bayside, worse at coughlin, 
## the same at duxbury, and non at powderpoint.

# Analysis
# Take code from above, but dont combine through the quads. Then do a t-test comparing sites
