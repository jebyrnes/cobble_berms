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
  summarise(phaeo_pc = sum(sum))

##
# Red Macrophytes
# Analysis of change in CHCR, MAST
##

phaeo_dat <- percent_dat |> 
  filter(species_code %in% c("CHCR", "MAST")) |> 
  group_by(site, treatment, height) |> 
  summarise(phaeo_pc = sum(sum))

