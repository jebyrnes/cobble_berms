#' -------------------------------------
#' Look at FARM counts
#' across berms versus controls and
#' -------------------------------------

source("scripts/helpers.R")


#load libraries
library(tidyverse)
library(broom)
library(lubridate)
library(glmmTMB)
library(performance)
library(car)
library(broom.mixed)

# load data
farm_dat <- read_csv("data/farms.csv")

# Work it up
farm_dat <- farm_dat |> 
  mutate(date_deployed = make_date(year_deployed, month_deployed, day_deployed)) |>
  mutate(date_retrieved = make_date(year_retrieved, month_retrieved, day_retrieved)) |>
  mutate(soak_duration = date_retrieved - date_deployed) |> 
  select(-c(4:14,16:20)) |> 
  select(1,2,3,5,6,7,4) |> 
  group_by(site, treatment, farms_number, date_deployed, date_retrieved, soak_duration, species_code) |>
  reframe(abundance = n()) |> 
  pivot_wider(names_from = species_code, values_from = abundance, values_fill = 0) |> 
  filter(!is.na(date_retrieved)) |> 
  select(-c(NOSP, LESP, CRFO, ULSL)) |> # get rid of nothing placeholder, ULSL is supposed to be ULSP which is an algae, not sure why it was recorded, CRFO only recoreded on one (i.e. not consistently recorded), ditto for lesp
  mutate(soak_class = if_else(condition = soak_duration > 100, "long", "short")) |> 
  rename(FUNS = FUSP) # correcting this code. FUSP is an algae, FUNS is a fish


# Calc abundance found per trap by retrieval season/soak duration
farm_soak_dat <- farm_dat |>
  select(1,2,20, 7:19) 

farm_soak_dat <- farm_soak_dat |> 
  group_by(site, treatment, soak_class) |> 
  summarise(
    n_farms = n(),
    across(
      where(is.numeric),
      \(x) sum(x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  mutate(
    across(
      -c(1:4),
      ~ round(.x / n_farms,2)))

  
  
# FPEN, CAMA, LILI, PHGU, PAHE are the codes that occur in abundance
# LILI exclusivly in the short/summer soak
# FPEN most abundant in the long soak, but still present in the short soak
# PAHE only in winter/long soak. Seasonality or duration thing?
# PHGU mostly in winter/long soak. Seasonality or duration thing?

# For this analysis, I will focus only on those spp which occur at some abundance (FPEN, CAMA, LILI, PHGU, PAHE)
# Also excluding coughlin due to not enough farms recovered

# CAMA
CAMA_mod <- glmmTMB(CAMA ~ treatment + soak_class, 
                     data = farm_dat |> filter(site == "Bayside"), 
                     family = poisson(link = "log"))

check_model(CAMA_mod)

# summarizing tests
Anova(CAMA_mod) |> tidy() 
tidy(CAMA_mod)

# FPEN
FPEN_mod <- glmmTMB(FPEN ~ treatment + soak_class, 
                    data = farm_dat |> filter(site == "Bayside"), 
                    family = poisson(link = "log"))

check_model(FPEN_mod)

# summarizing tests
Anova(FPEN_mod) |> tidy() 
tidy(FPEN_mod)

# LILI - Although exclusivly in the summer/short, so maybe should just do a ttest? 0 inflation for remainder
LILI_mod <- glmmTMB(LILI ~ treatment + soak_class, 
                    data = farm_dat |> filter(site == "Bayside"), 
                    family = poisson(link = "log"))

check_model(LILI_mod)

# summarizing tests
Anova(LILI_mod) |> tidy() 
tidy(LILI_mod)

# PHGU
PHGU_mod <- glmmTMB(PHGU ~ treatment + soak_class, 
                    data = farm_dat |> filter(site == "Bayside"), 
                    family = poisson(link = "log"))

check_model(PHGU_mod)

# summarizing tests
Anova(PHGU_mod) |> tidy() 
tidy(PHGU_mod)

# PAHE
PAHE_mod <- glmmTMB(PAHE ~ treatment + soak_class,
                    data = farm_dat |> filter(site == "Bayside"),
                    family = poisson(link = "log"))

check_model(PAHE_mod)

# summarizing tests
Anova(PAHE_mod) |> tidy() 
tidy(PAHE_mod)
