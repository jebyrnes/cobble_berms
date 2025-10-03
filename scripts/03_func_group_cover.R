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