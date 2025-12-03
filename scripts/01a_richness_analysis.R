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
seine_dat <- read_csv("data/seine.csv")
trap_dat <- read_csv("data/traps.csv")
farm_dat <- read_csv("data/farms.csv")
briv_dat <- read_csv("data/brivs.csv")
msl <- read_csv("data/master_sp_list.csv")

# Drop uneeded cols in msl and rename col
msl <- select(msl, 1:10) |> 
  rename(species_code = code)

# bind data together, and drop the two NAs in Bayside Berm
richness_dat <- rbind(quad_s23_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)),
                      quad_f23_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)), 
                      quad_s24_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)), 
                      quad_s25_dat |> 
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
  filter(species_code != "NOSP") |> 
  filter(species_code != "none_present") |> 
  filter(species_code != "UCRA")

# attach taxonomic data
richness_dat <- left_join(richness_dat, 
                          msl, 
                          "species_code")


##
# Write function to compute richness for a given taxonomic level
##
compute_richness <- function(df, level) {
  level_sym <- sym(level)
  
  df |> 
    distinct(site, treatment, !!level_sym) |> 
    count(site, treatment, name = paste0(level, "_richness"))
}

# List of taxonomic levels to include
tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

## 
# Calculate Richness by taxonomic level
## 

# Compute richness for each level
richness_list <- lapply(tax_levels, function(level) {
  compute_richness(richness_dat, level)
})

# Reduce the list of df into one
richness_summary <- Reduce(function(x, y) 
  full_join(x, y, by = c("site", "treatment")), richness_list)

# Rearrange
richness_summary <- richness_summary |> 
  arrange(site, treatment)

# Convert wide format to long format for plotting
richness_summary <- richness_summary |> 
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

# Create boxplots
richness_plot <- ggplot(richness_summary |> 
                          filter(taxonomic_level == "species"), 
                        mapping = 
         aes(x = treatment, 
             y = richness, 
             fill = treatment)) +
  geom_boxplot() +
  facet_wrap(vars(taxonomic_level), scales = "free_y") +
  labs(title = "Species Richness by Treatment",
       x = "Treatment",
       y = "Species Richness",
       fill = "Treatment") +
  scale_fill_manual(values = c("Berm" = "#af8dc3", "Control" = "#7fbf7b"))

# write out richness figure
ggsave("figures/richness_plot.jpg", 
       richness_plot,
       dpi = 600)

## Analysis of Species Richness

# Set up data
species_ttest_data <- richness_dat |> 
  select()
  filter(taxonomic_level == "species") |> 
  select(-"taxonomic_level") |> 
  pivot_wider(names_from = "treatment",
              values_from = "richness")
# Run ttest 
richness_ttest <- t.test(species_ttest_data$Berm, species_ttest_data$Control, paired = TRUE)


    
    
    
# Helper function that safely runs t.test and returns NA on error
safe_t_test <- function(x, y) {
  tryCatch(
    t.test(x, y, paired = TRUE),
    error = function(e) NULL
  )
}

# Run paired t-tests by level
richness_ttest <- richness_summary_modified |> 
  group_by(level) |> 
  summarise(
    t_test = list(safe_t_test(Berm, Control)),
    .groups = "drop"
  ) |> 
  # Remove levels where t-test failed (returned NULL)
  filter(!sapply(t_test, is.null)) |> 
  mutate(results = map(t_test, tidy)) |> 
  unnest(results) |> 
  select(level, estimate, statistic, p.value, conf.low, conf.high)


## Archive
# Pivot longer
richness_summary_modified <- pivot_longer(data = richness_summary,
             cols = "taxonomic_level",
             names_to = "level") |> 
  pivot_wider(names_from = "treatment")

# Helper function that safely runs t.test and returns NA on error
safe_t_test <- function(x, y) {
  tryCatch(
    t.test(x, y, paired = TRUE),
    error = function(e) NULL
  )
}

# Run paired t-tests by level
richness_ttest <- richness_summary_modified |> 
  group_by(level) |> 
  summarise(
    t_test = list(safe_t_test(Berm, Control)),
    .groups = "drop"
  ) |> 
  # Remove levels where t-test failed (returned NULL)
  filter(!sapply(t_test, is.null)) |> 
  mutate(results = map(t_test, tidy)) |> 
  unnest(results) |> 
  select(level, estimate, statistic, p.value, conf.low, conf.high)
