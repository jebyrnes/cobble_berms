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



## Calculate Richness by taxonomic level

# Write a function to compute richness for a given taxonomic level
compute_richness <- function(df, level) {
  level_sym <- sym(level)
  
  df |> 
    distinct(site, treatment, !!level_sym) |> 
    count(site, treatment, name = paste0(level, "_richness"))
}

# List of taxonomic levels to include
tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

# Compute richness for each level
richness_list <- lapply(tax_levels, function(level) {
  compute_richness(richness_dat, level)
})

# Reduce the list of data frames into one by joining on site and treatment
richness_summary <- Reduce(function(x, y) full_join(x, y, by = c("site", "treatment")), richness_list)

# Rearrange
richness_summary <- richness_summary |> arrange(site, treatment)



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