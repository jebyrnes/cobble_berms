#' -------------------------------------
#' Look at snail counts
#' across berms versus controls and
#' quads count only
#' -------------------------------------
source("scripts/helpers.R") # has ggplot2

## load libraries
pacman::p_load(dplyr,
               tidyr,
               readr,
               glmmTMB, # zero inflation!
               car,
               broom.mixed,
               emmeans
               )

## load the data and filter to Littorina littorea
## then sum for each quad
snail_dat <- read_csv("data/quads.csv") |>
  filter(measurement_type == "Count") |>
  filter(species_code == "LILI") |>
  group_by(site, treatment, quadrat, height) |>
  summarize(snails = sum(measurement, na.rm = TRUE)) |>
  ungroup() |> fix_data_for_plots()

## Fix up snail data to fill in 0 for quads where there
## were no snails
snail_dat <- snail_dat |>
  complete(site, treatment, nesting(quadrat, height),
           fill = list(snails = 0)) |>
  mutate(are_there_snails = as.numeric(snails !=0))


## an initial plot
snail_plot <- ggplot(snail_dat,
       aes(x = site, color = treatment,
           y = snails)) +
  geom_point(position = position_dodge(width = 0.5),
             alpha = 0.5,
             size = 2) +
  scale_color_manual(values = berm_control_colors) +
  facet_wrap(vars(height)) +
  labs(x="", 
       y = "Littorina per\n0.25 sq m",
       color = "") 

# write out snail figure
ggsave("figures/snail_plot.jpg", 
       snail_plot,
       dpi = 600)

# what's the distribution of snails?
ggplot(snail_dat, 
       aes(x = snails)) +
  geom_histogram() +
  facet_wrap(vars(height)) +
  labs(x = "Littorina per\n0.25 sq m")

##
# Analysis of where are there even snails
##
snail_pres_mod <- glm(are_there_snails ~
                        height * treatment*site,
                      data = snail_dat,
                      family = binomial)

Anova(snail_pres_mod) |> tidy() |>
  mutate(p.value = round(p.value, digits = 10))

## Complete separation problem - write about things
## without stats in the text!

##
# Analysis of Low tide height
##
## an analysis using zero inflated negative binomial
## just for the low where we don't have site:treatments
## that are completely 0
snail_mod <- glmmTMB(snails ~ treatment*site,
                     zi = ~  site,
                     data = snail_dat |> filter(height == "Low"),
                     family = nbinom2)

# model diagnostics
performance::check_model(snail_mod)


# summarizing tests
Anova(snail_mod) |> tidy()
tidy(snail_mod)

# posthoc
snail_em <- emmeans(snail_mod,
                    ~treatment | site)

contrast(snail_em, method = "pairwise")





