#' -------------------------------------
#' Look at functional group cover
#' across berms versus controls and
#' quads cover only
#' -------------------------------------
source("scripts/helpers.R")

## Planning
#' Determine what func groups we are interested in: Habitat formers (ASCO, FUSP, CHCR/MAST), 
#' sessile inverts (SEBA, MYED)
#' lm brown_macrophytes ~ site + treatment + zone
#' EXCLUDE MYED from year 1 due to count-PC issue

## load libraries
library(tidyverse)
library(broom)
library(vegan)
library(patchwork)

## load data
quad_s23_dat <- read_csv("data/quads_s23.csv")
quad_f23_dat <- read_csv("data/quads_f23.csv")
quad_s24_dat <- read_csv("data/quads_s24.csv")
quad_s25_dat <- read_csv("data/quads_s25.csv")

# Bind data
func_dat <- rbind(quad_s23_dat |> 
                        group_by(season, site, treatment) |> 
                        reframe(species_code = unique(species_code)),
                      quad_f23_dat |> 
                        group_by(season, site, treatment) |> 
                        reframe(species_code = unique(species_code)), 
                      quad_s24_dat |> 
                        group_by(season, site, treatment) |> 
                        reframe(species_code = unique(species_code)), 
                      quad_s25_dat |> 
                        group_by(season, site, treatment) |> 
                        reframe(species_code = unique(species_code))
) |> 
  filter(species_code != is.na(species_code)) |> 
  filter(!species_code %in% c("none_persent", "NOSP"))

# Drop data that are no longer needed
rm(quad_s23_dat)
rm(quad_f23_dat)
rm(quad_s24_dat)
rm(quad_s25_dat)

# Select relevant spp, and sum across squares
func_dat <- func_dat |>
  filter(species_code %in% c("ASNO", "FUCU", "FUSP", "FUVE", "FUDI", "CHCR", "MAST", "SEBA", "MYED")) |> 
  group_by(site, treatment, height, quadrat, measurement_type, species_code) |> 
  reframe(sum = sum(measurement)) |> 
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
