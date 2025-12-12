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

##
# Exploring the data
##

# Abundance of each spp per trap
abundance_dat <- trap_dat |>
    group_by(site, treatment, trap_name, trap_number, species_code) |>
    reframe(abundance = n()) |> 
  pivot_wider(names_from = species_code, values_from = abundance, values_fill = 0) |> 
  select(-NOSP)

# CAMA only
CAMA_dat <- abundance_dat |>
  select(1:5) |> 
  filter(trap_name == "Crab") |> 
  select(-trap_name)

# All spp for table
abundance_dat <- abundance_dat |>
  group_by(site, treatment, trap_name) |> 
  summarise(
    n = n_distinct(trap_number),
    across(CAMA:FUNS, ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )


##
# Visualize
##

# CAMA abundance by site as a function of treatment
CAMA_plot <- ggplot(CAMA_dat, 
                    mapping = aes(x = treatment, 
                              y = CAMA)) +
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
  CAMA_ttest <- t.test(CAMA ~ treatment, data = CAMA_dat, var.equal = FALSE) |> 
    tidy()

# estimate estimate1 estimate2 statistic p.value parameter conf.low   conf.high 
# 1.73      6.87      5.13     0.573   0.572      24.5    -4.51      7.97
#   
# Conclusions: no difference in CAMA abundance in a trap given a treatment
  