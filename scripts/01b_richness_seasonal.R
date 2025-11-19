#' -------------------------------------
#' Look at total species richness
#' across berms versus controls and
#' quads (all), seines, traps, and brivs combined
#' -------------------------------------
source("scripts/helpers.R")

## load libraries
library(tidyverse)
library(broom)

## load data, bind them together, and add taxonomic data
quad_s23_dat <- read_csv("data/quads_s23.csv")
quad_f23_dat <- read_csv("data/quads_f23.csv")
quad_s24_dat <- read_csv("data/quads_s24.csv")
quad_s25_dat <- read_csv("data/quads_s25.csv")
msl <- read_csv("data/master_sp_list.csv")

# Drop uneeded cols in msl and rename col
msl <- select(msl, 1:10) |> 
  rename(species_code = code)

# bind data together, and drop the two NAs in Bayside Berm
seasonal_dat <- rbind(quad_s23_dat |> 
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
  filter(species_code != "NOSP")

# attach taxonomic data
seasonal_dat <- left_join(seasonal_dat, 
                                   msl, 
                                   "species_code")

## Write a function to calculate richness by season
seasonal_richness <- function(df, level) {
  level_sym <- sym(level)
  
  df |> 
    distinct(season, site, treatment, !!level_sym) |> 
    count(season, site, treatment, name = paste0(level, "_richness"))
}

# List of taxonomic levels to include
tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

## Calculate Richness by taxonomic level and season

# Compute richness for each level
seasonal_list <- lapply(tax_levels, function(level) {
  seasonal_richness(seasonal_dat, level)
})

# Reduce the list of df into one
seasonal_richness <- Reduce(function(x, y) 
  full_join(x, y, by = c("season", "site", "treatment")), seasonal_list)

# Rearrange
seasonal_richness <- seasonal_richness |> 
  arrange(site, treatment)

# Convert wide format to long format for plotting
seasonal_richness <- seasonal_richness |> 
  pivot_longer(
    cols = ends_with("_richness"),
    names_to = "taxonomic_level",
    values_to = "richness"
  ) |>
  mutate(
    taxonomic_level = gsub("_richness", "", taxonomic_level), 
    taxonomic_level = factor(taxonomic_level, 
                             levels = c("kingdom", "phylum", "class", "order", 
                                        "family", "genus", "species")) 
  )

##
# Visualize
##

# time series of genus richness faceted by site
seasonal_plot <- ggplot(seasonal_richness |> 
                          filter(taxonomic_level == "genus"), 
                        mapping = 
                          aes(x = season, 
                              y = richness, 
                              shape = treatment)) +
  geom_point() +
  facet_wrap(vars(site), scales = "free_y")

# write out richness figure
ggsave("figures/richness_plot.jpg", 
       richness_plot,
       dpi = 600)
