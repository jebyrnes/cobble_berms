#' -------------------------------------
#' Look at seine count data
#' across berms versus controls and
#' -------------------------------------


source("scripts/helpers.R")

##
# Set-up
##

#load libraries
library(tidyverse)
library(broom)
library(lubridate)

# load data
seine_dat <- read_csv("data/seine.csv")

# Add date column
seine_dat <- seine_dat |> 
  mutate(date = make_date(year, month, day))

# Estimate Volume sampled
seine_dat <- seine_dat|>
  group_by(site, treatment, date) |>
  mutate(vol_sampled = ((estimated_depth_pole_1_m/estimated_depth_pole_2_m)/2)*width_net_aperture_m*distance_walked_m)

# Drop cols not interested in anymore
seine_dat <- seine_dat |> 
  select(-c(3:17,19:22))

# Reframe for abundance per m3 sampled
seine_dat <- seine_dat|>
  group_by(site, treatment, date, vol_sampled, species_code) |>
  reframe(abundance = n()) |> 
  mutate(concentration = abundance/vol_sampled)

# Make a table showing concentration by species for each tow
seine_table_dat <- seine_dat |> 
  select(-abundance) |> 
  mutate(concentration_100m3 = concentration*100) |>
  select(-concentration) |> 
  pivot_wider(names_from = species_code, 
              values_from = concentration_100m3,
              values_fill = 0) |> 
  mutate(across(where(is.numeric), ~ ifelse(.x == 0, "0", format(round(.x, 2), nsmall = 2))))

seine_table <- seine_table_dat |> 
  kable("html", caption = "Seine Haul Species Abundance per 100 m3", digits = 2) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) %>%
  row_spec(0, bold = TRUE, background = "#D3D3D3")

##
# Visualize
##

# CAMA abundance by site as a function of treatment
#seine_plot <- 
  ggplot(seine_dat, 
                    mapping = aes(x = treatment, 
                                  y = concentration)) +
  geom_boxplot() +
  geom_point(position = "jitter",
             alpha = .6) +
  labs(title = "Individual concentration per m3",
       x = "Treatment",
       y = "Concentration (indiv/m3)") +
  facet_wrap(vars(site))

##
# Stats
##

# ttest for CAMA abundnace
CAMA_seine_dat <- seine_dat |> 
  filter(species_code == "CAMA")

CAMA_seine_ttest <- t.test(data = CAMA_seine_dat, 
                     concentration ~ treatment, var.equal = FALSE) 

#t = -1.1274, df = 7.4981, p-value = 0.2944
# Conclusions: no difference in CAMA abundance in a seine given a treatment

## NEXT: Does probability of finding a different species vary with treatment?
