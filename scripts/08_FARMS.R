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
library(DHARMa)
library(emmeans)

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
  select(-c(NOSP, LESP, CRFO, ULSL, LILI)) |> # get rid of nothing placeholder, ULSL is supposed to be ULSP which is an algae, not sure why it was recorded, CRFO only recoreded on one (i.e. not consistently recorded), ditto for lesp. Also droppuing LILI, not assessed in winter I dont think
  mutate(soak_class = if_else(condition = soak_duration > 100, "long", "short")) |> 
  rename(FUNS = FUSP) # correcting this code. FUSP is an algae, FUNS is a fish


# Calc abundance found per farm by retrieval season/soak duration
farm_soak_dat <- farm_dat |>
  select(1,2,19, 7:18) 

farm_soak_dat <- farm_soak_dat |> 
  group_by(site, treatment, soak_class) |> 
  summarise(
    n_farms = n(),
    across(
      where(is.numeric),
      \(x) sum(x, na.rm = TRUE)), # sum up across farms
    .groups = "drop"
  ) |>
  mutate(
    across(
      -c(1:4),
      ~ round(.x / n_farms,2))) # standardize by count per farm

# Excluding coughlin due to not enough farms recovered
bayside_dat <- farm_dat |>
  filter(site == "Bayside",
         !is.na(FPEN),
         !is.na(PHGU),
         !is.na(treatment),
         !is.na(soak_class)) |>
  mutate(treatment = factor(treatment),
         soak_class = factor(soak_class))

# Focusing on PHGU and FPEN
# FPEN most abundant in the long soak, but still present in the short soak
# PHGU mostly in winter/long soak. Seasonality or duration thing?


##
# FPEN -------------------------------------------
##

# FPEN_mod <- lm(FPEN ~ treatment + soak_class,
#                     data = farm_dat |> filter(site == "Bayside"))

FPEN_nb <- glmmTMB(FPEN ~ treatment + soak_class + treatment *soak_class, 
                    data = bayside_dat, 
                    family = nbinom2(link = "log"))

# FPEN_pois <- glmmTMB(FPEN ~ treatment + soak_class,
#                      data = bayside_dat,
#                      family = poisson(link = "log"))

AIC(FPEN_mod, FPEN_pois, FPEN_nb) # FPEN_pois  3df,  218.2968aic; FPEN_nb    4 df, 129.9723aic
# nb is best fit

check_model(FPEN_nb)

sim_nb <- simulateResiduals(
  fittedModel = FPEN_nb,
  n = 1000)

plot(sim_nb) # Problems with residuals v predicted quantile deviations, fixed by adding interaction
# 
# testDispersion(sim_nb)
# testOutliers(sim_nb)

summary(FPEN_nb)
tidy(FPEN_nb)

FPEN_emm <- emmeans(FPEN_nb, 
        ~ treatment | soak_class,
        type = "response")

pairs(FPEN_emm, adjust = "holm")

# Visualie Predicted shrimp abundance on bayside

# Create all df with all treatment and soak class combinations
pred_dat <- expand.grid(
  treatment = levels(bayside_dat$treatment),
  soak_class = levels(bayside_dat$soak_class))

# Predict (on link/log scale)
pred <- predict(FPEN_nb,
  newdata = pred_dat,
  type = "link",
  se.fit = TRUE,
  re.form = NA)

# Convert predictions/CI to count
pred_dat <- pred_dat |>
  mutate(
    predicted = exp(pred$fit),
    lower = exp(pred$fit - 1.96 * pred$se.fit),
    upper = exp(pred$fit + 1.96 * pred$se.fit))

# Plot predicted counts
ggplot(
  pred_dat,
  aes(
    x = treatment,
    y = predicted,
    color = soak_class,
    group = soak_class
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.1
  ) +
  labs(
    x = "Treatment",
    y = "Predicted shrimp abundance",
    color = "Soak Duration")


##
# PHGU ------------------
##

# PHGU - Bayside only, long only (short has a count of 1 on the berm, all others 0)
# PHGU_lm <- lm(PHGU ~ treatment, 
#                     data = bayside_dat |> filter(soak_class == "long"))

# PHGU_nb <- glmmTMB(PHGU ~ treatment + soak_class, 
#                    data = bayside_dat, 
#                    family = nbinom2(link = "log"))
# 
PHGU_pois <- glmmTMB(PHGU ~ treatment,
                      data = bayside_dat |> filter(soak_class == "long"),
                      family = poisson(link = "log"))
# 
# AIC(PHGU_lm, PHGU_pois)


# sim_PHGU_lm <- simulateResiduals(
#   fittedModel = PHGU_lm,
#   n = 1000)

# plot(sim_PHGU_lm)
# testDispersion(sim_PHGU_lm)
# testOutliers(sim_PHGU_lm)

summary(PHGU_pois)
tidy(PHGU_pois)



# check_model(PHGU_mod)
# 
# # summarizing tests
# Anova(PHGU_mod) |> tidy() 
# tidy(PHGU_mod)
# Lose PHGU in the control and the short


ggplot(farm_dat |> filter(site == "Bayside"), 
       mapping = aes(x = treatment, 
                     y = PHGU,
                     fill = treatment)) +
  geom_boxplot() +
  geom_jitter(aes(shape = soak_class),  
              width = 0.2,
              height = 0, 
              alpha = 0.6, 
              size = 2) +
  labs(title = "PHGU abundance in farms",
       x = "Treatment",
       y = "PHGU Abundance") +
  facet_wrap(vars(soak_class))
