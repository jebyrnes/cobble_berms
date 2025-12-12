#' -------------------------------------
#' Look at trap count data
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
trap_dat <- read_csv("data/traps.csv")

## Verify same amount of traps per site
# trap_dat |>
#   group_by(site, treatment, trap_name) |>
#   reframe(trap_count = n_distinct(trap_number)) |>
#   View()
# All there!

# Drop the NOSP
trap_dat <- trap_dat |> 
  filter(species_code != "NOSP")

##
# Exploring the data
##

# Abundance of each spp per trap
abundance_dat <- trap_dat |>
    group_by(site, treatment, trap_name, trap_number, species_code) |>
    reframe(abundance = n()) |> 
  pivot_wider(names_from = species_code, values_from = abundance, values_fill = NA)

# CAMA only
CAMA_dat <- trap_dat |>
  filter(species_code == "CAMA") |> 
  group_by(site, treatment, trap_name, trap_number) |>
  reframe(abundance = n())

##
# Visualize
##

# CAMA abundance by site as a function of treatment
CAMA_plot <- ggplot(CAMA_dat, 
                    mapping = aes(x = treatment, 
                              y = abundance)) +
  geom_boxplot() +
  geom_point(position = "jitter",
              alpha = .6)+
  labs(title = "CAMA abundance in traps",
       x = "Treatment",
       y = "CAMA Abundance") +
    facet_wrap(vars(site))

##
# Stats
##

# ttest for CAMA abundnace
  CAMA_ttest <- t.test(abundance ~ treatment, data = CAMA_dat, var.equal = FALSE) |> 
    tidy()

#  estimate estimate1 estimate2 statistic p.value parameter conf.low conf.high 
#  0.0814      6.24      6.15    0.0283   0.978      27.8    -5.82      5.98 
  
# Conclusions: no difference in CAMA abundance in a trap given a treatment
  