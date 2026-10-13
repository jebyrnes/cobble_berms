#' -------------------------------------
#' Look at quadrat level richness, shannon diversity, and eveness
#' across berms versus controls and
#' for percent cover and count seperatly
#' -------------------------------------
source("scripts/helpers.R")

##
#  Set-up ----------------------------------
##


# load libraries
library(tidyverse)
library(easystats)
library(broom)
library(broom.mixed)
library(vegan)
library(car)
# library(patchwork)
library(emmeans)
library(brglm2)
library(glmmTMB)

# load data
quad_s23_dat <- read_csv("data/quads_s23.csv")
quad_f23_dat <- read_csv("data/quads_f23.csv")
quad_s24_dat <- read_csv("data/quads_s24.csv")
quad_s25_dat <- read_csv("data/quads_s25.csv")

# Bind together
quad_dat <- rbind(quad_s23_dat,
                  quad_f23_dat,
                  quad_s24_dat,
                  quad_s25_dat)

# rm/ spaces in col titles and change none_present to NOSP, and change measurement type = NA to PErcent, also drop the duplicate entry of a quad
quad_dat <- quad_dat |> 
  mutate(measurement_type = if_else(measurement_type == "Percent Cover", "Percent", measurement_type))

quad_dat <- quad_dat |> 
  mutate(species_code = if_else(species_code == "none_present", "NOSP", species_code),
         measurement_type = if_else(species_code == "NOSP" & is.na(measurement_type),"Percent", measurement_type),
         height = factor(height, levels = c("Low", "Mid", "High"))) |> 
  slice(-808)


# Drop data that are no longer needed
rm(quad_s23_dat)
rm(quad_f23_dat)
rm(quad_s24_dat)
rm(quad_s25_dat)

## Verify correct amount of quads done
quad_dat |>
  group_by(site, treatment, season, height) |>
  reframe(quad_count = n_distinct(quadrat)) |>
  View()
  
# Problem data: 
# Coughlin-Control-fall_23-High-2 quads only
# Coughlin-Control-fall_23-Mid-4
# Bayside-Berm-summer_25-High-5
# Coughlin-Berm-summer_25-Low-5

# Store total quad count for later use
n_total <- quad_dat |>
group_by(site, treatment, season, height) |>
  summarise(quad_count = n_distinct(quadrat), .groups = "drop") |>
  summarise(n_total = sum(quad_count))

# Make sure every quad with 0 has at a NOSP place holder for percent and for count
quad_dat <- quad_dat |> 
  filter(species_code == "NOSP") |> 
  group_by(season, site, treatment, height, quadrat) |> 
  complete(measurement_type = c("Percent", "Count"),
    fill = list(value = 0, species_code = "NOSP")) |> 
  ungroup() |>
  bind_rows(quad_dat %>% filter(species_code != "NOSP"))

# quad_dat |>
#   filter(species_code == "NOSP") |>
#   count(season, site, treatment, height, quadrat) |> 
#   View()
# All there!

# Sum to the quad level
quad_dat <- quad_dat |> 
  group_by(season, site, treatment, height, quadrat, measurement_type, species_code) |> 
  reframe(sum = sum(measurement)) |> 
  ungroup()

# Pivot data for use in vegan::diversity() 
quad_dat <- quad_dat |> 
  pivot_wider(names_from = "species_code",
              values_from = "sum",
              values_fill = 0)

# Drop placeholder spp for quads w/ nothing and remove na
quad_dat <- quad_dat |> 
  select(-c("NOSP"))

quad_dat[is.na(quad_dat)] <- 0

##
#  Richness and Shannon Biodiversity set up ----------------------------------
##

# Calculate Shannon, InvSimpson, richness, evenness, and rearrange
div_dat <- quad_dat |> 
  mutate(shannon = diversity(quad_dat[,-c(1:6)], index = "shannon")) # measures the uncertainty in predicting the species identity of an individual chosen at random from the community. A higher value indicates higher diversity (more species and/or more even distribution), and a lower value indicates lower diversity. It equals 0 only when there is no uncertainty—meaning there is only ONE species present.

div_dat <- div_dat |>  # not going to analyze. Shannon and richness should suffice
  mutate(invsimpson = diversity(div_dat[,-c(1:6,85)], index = "invsimpson")) 

div_dat <- div_dat |> 
  mutate(richness = specnumber(div_dat[,-c(1:6,85,86)])) 

div_dat <- div_dat |>  # Going to ignore this in anlysis since they're similar and I think this will be better captured by the analysis in the cover composition section
  mutate(evenness = shannon/log(richness)) # Pielous eveness how evenly individuals are distributed across the different species in a community. It ranges from 0 (completely uneven) to 1 (perfectly even).# Evenness of 0 occurs when richness is 1 (i.e. log 1 = 0, dividing by 0 = NaN)
  
div_dat <- div_dat |> 
select(richness, shannon, evenness, invsimpson, everything())
  
# Set up data for plotting
div_plot_dat <- div_dat[,-c(11:88)] |> 
  pivot_longer(cols = richness:invsimpson, 
               values_to = "value", 
               names_to = "index")

div_plot_dat$height <- factor(div_plot_dat$height,
                              levels = c("Low", "Mid", "High"))


##
# Visualize
##


# ggplot(div_plot_dat |> filter(measurement_type == "Count"),
#        aes(x = treatment,
#            y = value)) +
#   geom_violin() +
#   facet_grid(rows = vars(index), 
#              cols = vars(height),
#              scales = "free_y")
# 
# ggplot(div_plot_dat |> filter(measurement_type == "Percent"),
#        aes(x = treatment,
#            y = value)) +
#   geom_violin() +
#   facet_grid(rows = vars(index), 
#              cols = vars(height),
#              scales = "free_y")

##
#  Shannon analysis ----------------------------------
##

## Shannon
# Hurdle Model - Odds of being over 0
shannon_percent_hurdle_mod <- glm(shannon > 0 ~ treatment + height + site + treatment:height,
                                  family = binomial(link = "logit"),
                                  data = div_dat |> filter(measurement_type == "Percent"))

shannon_count_hurdle_mod <- glm(shannon > 0 ~ treatment + height + site + treatment:height,
                                  family = binomial(link = "logit"),
                                  data = div_dat |> filter(measurement_type == "Count"))

# Models for those over 0
shannon_percent_posthurdle_mod <- lm(sqrt(shannon) ~ treatment + height + site + treatment:height, # note that shannon at richness of 1 is 0, no diversity, only 1 sp present which is the same as div of 0
                          data = div_dat |> filter(measurement_type == "Percent", shannon > 0))

shannon_count_posthurdle_mod <- lm(sqrt(shannon) ~ treatment + height + site + treatment:height,
                        data = div_dat |> filter(measurement_type == "Count", shannon > 0))

# Get model stats
# check_model(shannon_percent_hurdle_mod) # Homogeneity of variance seems to be the only problem
# check_model(shannon_count_hurdle_mod) # Homogeneity of variance, also high collinearity for interaction parameter
# check_model(shannon_percent_posthurdle_mod) # PPC could be better, Homogeneity of variance not great
# check_model(shannon_count_posthurdle_mod) # treatment and treatment:height has high colineaerity

## So what does this tell us
# Hurdle (0,1) for shannon percent
Anova(shannon_percent_hurdle_mod) |> tidy() # height plys biggest roll, but site and treatment also play a roll
tidy(shannon_percent_hurdle_mod)

shannon_percent_hurdle_em <- emmeans(shannon_percent_hurdle_mod,
                    ~treatment | height)

contrast(shannon_percent_hurdle_em, method = "pairwise") # Significant effect of treatment in the low and mid


# Hurdle (0,1) for shannon count
Anova(shannon_count_hurdle_mod) |> tidy() # all important, by again height the most so
tidy(shannon_count_hurdle_mod) # need to back transform the estimate?

shannon_count_hurdle_em <- emmeans(shannon_count_hurdle_mod,
                                     ~treatment | height)

contrast(shannon_count_hurdle_em, method = "pairwise") # Significant effect of treatment in the low and mid

emmeans(
  shannon_count_hurdle_mod,
  ~ treatment * height
) |>
  contrast(
    interaction = "pairwise"
  )

# Post-hurdle (>0) for shannon percent
Anova(shannon_percent_posthurdle_mod) |> tidy() # Percent Shannon  differs significantly across heights and sites, but is not strongly affected by treatment, and the interaction between treatment and height is also not significant.
tidy(shannon_percent_posthurdle_mod) |>  
  mutate(estimate_original = sign(estimate) * estimate^2,
         std.error_original = 2 * abs(estimate) * std.error) |> 
  mutate(lower.CL = ((estimate_original - std.error_original)^2),
         upper.CL = ((estimate_original + std.error_original)^2))

summary(shannon_percent_posthurdle_mod)
summary(shannon_percent_posthurdle_mod)$r.squared
summary(shannon_percent_posthurdle_mod)$adj.r.squared

shannon_percent_posthurdle_em <- emmeans(shannon_percent_posthurdle_mod,
                                   ~treatment | height)

contrast(shannon_percent_posthurdle_em, method = "pairwise") # low and mid treatment has a significant effect

# Post-hurdle (>0) for shannon count
Anova(shannon_count_posthurdle_mod) |> tidy() # essentially the same as percent
tidy(shannon_count_posthurdle_mod) |>  # height has largest effect, while there is some decrease between treatment
  mutate(estimate_original = sign(estimate) * estimate^2, # note that the negative estimates are lost after the back transformation, so i've multiplioed them by the sign to return the original sign
         std.error_original = 2 * abs(estimate) * std.error) |> 
  mutate(lower.CL = ((estimate_original - std.error_original)^2),
         upper.CL = ((estimate_original + std.error_original)^2))


shannon_count_posthurdle_em <- emmeans(shannon_count_posthurdle_mod,
                                   ~treatment | height)

contrast(shannon_count_posthurdle_em, method = "pairwise") # no significance

##
#  Richness analysis ----------------------------------
##

# Just realized that this doesn't need to be done seperatly for count and percent. Redoing the analysis.
rich_dat <- div_dat |> 
  group_by(season, site, treatment, height, quadrat) |> 
  summarise(richness = sum(richness))

#Hurdle mod
richness_hurdle_mod <- glm(richness > 0 ~ treatment + height + site + treatment:height, 
                                   family = binomial,
                                   data = rich_dat,
                                   method = "brglmFit")

# Check models predicted probabilities
emmeans(richness_hurdle_mod,
        ~ treatment | height,
        type = "response")

# And compare to the observed probabilities -- looks good!
rich_dat |>
  group_by(treatment, height) |>
  summarise(
    n = n(),
    nonzero = sum(richness > 0),
    prop_nonzero = mean(richness > 0),
    .groups = "drop"
  )


# Models for those over 0
richness_posthurdle_mod <- lm(sqrt(richness) ~ treatment + height + site + treatment:height,
                                      data = rich_dat |> 
                                filter(richness > 0))
# RESUME HERE ----------------------------------
# Potentially worth considering:
richness_posthurdle_nb <- glmmTMB(
  richness ~ treatment + height + site + treatment:height,
  family = truncated_nbinom2(link = "log"),
  data = rich_dat |> 
    filter(richness > 0)
)

# summary(richness_posthurdle_nb)
diagnose(richness_posthurdle_nb)
richness_posthurdle_nb_em <- emmeans(
  richness_posthurdle_nb,
  ~ treatment | height,
  type = "response"
)

summary(richness_posthurdle_nb_em)
pairs(richness_posthurdle_nb_em)

library(DHARMa)

sim_nb <- simulateResiduals(
  richness_posthurdle_nb,
  n = 1000
)

plot(sim_nb)
  
  testDispersion(sim_nb)
testZeroInflation(sim_nb)
testUniformity(sim_nb)

# Get model stats
# check_model(richness_hurdle_mod)
# check_model(richness_posthurdle_mod)

# Interpretation

## Hurdle
summary(richness_hurdle_mod)
richness_hurdle_em <- emmeans(richness_hurdle_mod,
                                     ~treatment | height)

richness_hurdle_em_rg <- regrid(richness_hurdle_em,
                                        transform = "response")

pairs(richness_hurdle_em_rg)
summary(richness_hurdle_em_rg)

# Richness Post-hurdle
Anova(richness_posthurdle_mod) |> tidy() # effect of height and site, not treatment 
tidy(richness_posthurdle_mod) |>  
  mutate(estimate_original = sign(estimate) * estimate^2,
         std.error_original = 2 * abs(estimate) * std.error)|> 
  mutate(lower.CL = ((estimate_original - std.error_original)^2),
         upper.CL = ((estimate_original + std.error_original)^2))

summary(richness_posthurdle_mod)
richness_posthurdle_em <- emmeans(richness_posthurdle_mod,
                                         ~treatment | height)

contrast(richness_posthurdle_em, method = "pairwise")

richness_posthurdle_em_rg <- regrid(
  richness_posthurdle_em,
  transform = "response"
)

summary(richness_posthurdle_em_rg)

pairs(richness_posthurdle_em_rg)

# ##
# # Analysis of non-selected indices. Not including this in report and thus did not finish these
# ##
# 
# ## Evenness
# # Hurdle Model - Odds of being over 0
# evenness_percent_hurdle_mod <- glm(evenness > 0 ~ treatment + height + site,
#                                   family = binomial(link = "logit"),
#                                   data = div_dat |> filter(measurement_type == "Percent"))
# 
# evenness_count_hurdle_mod <- glm(evenness > 0 ~ treatment + height + site,
#                                 family = binomial(link = "logit"),
#                                 data = div_dat |> filter(measurement_type == "Count"))
# 
# # Models for those over 0
# evenness_percent_posthurdle_mod <- lm(sqrt(evenness) ~ treatment + height + site,
#                                      data = div_dat |> filter(measurement_type == "Percent"))
# 
# evenness_count_posthurdle_mod <- lm(sqrt(evenness) ~ treatment + height + site,
#                                    data = div_dat |> filter(measurement_type == "Count"))
# 
# # Get model stats
# check_model(evenness_percent_hurdle_mod)
# check_model(evenness_count_hurdle_mod)
# check_model(evenness_percent_posthurdle_mod)
# check_model(evenness_count_posthurdle_mod)
# 
# ## invsimpson
# # Hurdle Model - Odds of being over 0
# invsimpson_percent_hurdle_mod <- glm(invsimpson > 0 ~ treatment + height + site,
#                                    family = binomial(link = "logit"),
#                                    data = div_dat |> filter(measurement_type == "Percent"))
# 
# invsimpson_count_hurdle_mod <- glm(invsimpson > 0 ~ treatment + height + site,
#                                  family = binomial(link = "logit"),
#                                  data = div_dat |> filter(measurement_type == "Count"))
# 
# # Models for those over 0
# invsimpson_percent_posthurdle_mod <- lm(sqrt(invsimpson) ~ treatment + height + site,
#                                       data = div_dat |> filter(measurement_type == "Percent"))
# 
# invsimpson_count_posthurdle_mod <- lm(sqrt(invsimpson) ~ treatment + height + site,
#                                     data = div_dat |> filter(measurement_type == "Count"))
# 
# # Get model stats
# check_model(invsimpson_percent_hurdle_mod)
# check_model(invsimpson_count_hurdle_mod)
# check_model(invsimpson_percent_posthurdle_mod)
# check_model(invsimpson_count_posthurdle_mod)

##
# Viz
##

# Shannon, post-hurdle
## by percent cover

# Summarize mean and standard error for plotting
plot_data <- div_dat %>%
  filter(measurement_type == "Percent") %>%
  group_by(height, site, treatment) %>%
  summarise(
    mean_shannon = mean(shannon, na.rm = TRUE),
    se_shannon = sd(shannon, na.rm = TRUE)/sqrt(n()),
    .groups = "drop") |> 
  mutate(height = factor(height, levels = c("Low", "Mid", "High")))

ggplot(plot_data, aes(x = treatment, y = mean_shannon, color = height, group = height)) +
  geom_line(position = position_dodge(width = 0.8)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = mean_shannon - se_shannon,
        ymax = mean_shannon + se_shannon),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~site) +
  labs(
    y = "Shannon Diversity",
    x = "Treatment",
    title = "Effect of Height and Site on Shannon Diversity for Percent Cover"
  )



# Plot that shows shannon decreasing by zone and site, but relatively unchanged by treatment
ggplot(plot_data, aes(x = height, y = mean_shannon, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon),
                width = 0.2, position = position_dodge(width = 0.8)) +
  facet_wrap(~treatment) +
  labs(
    y = "Shannon Diversity",
    x = "Height",
    title = "Effect of Height and Site on Shannon Diversity for Percent Cover")

# Plot that shows shannon unchanged by treatment
ggplot(plot_data, aes(x = height, y = mean_shannon, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon),
                width = 0.2, position = position_dodge(width = 0.8)) +
  facet_wrap(~site) +
  labs(
    y = "Shannon Diversity",
    x = "Height",
    title = "Effect of Treatment on Shannon Diversity for Percent Cover")

## by count

# Summarize mean and standard error for plotting
plot_data <- div_dat %>%
  filter(measurement_type == "Count") %>%
  group_by(height, site, treatment) %>%
  summarise(
    mean_shannon = mean(shannon, na.rm = TRUE),
    se_shannon = sd(shannon, na.rm = TRUE)/sqrt(n()),
    .groups = "drop") |> 
  mutate(height = factor(height, levels = c("Low", "Mid", "High")))

# Plot that shows shannon decreasing by zone and site, but relatively unchanged by treatment
ggplot(plot_data, aes(x = height, y = mean_shannon, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon),
                width = 0.2, position = position_dodge(width = 0.8)) +
  facet_wrap(~treatment) +
  labs(
    y = "Shannon Diversity",
    x = "Height",
    title = "Effect of Height and Site on Shannon Diversity for Count")

# Plot that shows shannon unchanged by treatment
ggplot(plot_data, aes(x = height, y = mean_shannon, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon),
                width = 0.2, position = position_dodge(width = 0.8)) +
  facet_wrap(~site) +
  labs(
    y = "Shannon Diversity",
    x = "Height",
    title = "Effect of Treatment on Shannon Diversity for Count")

# Richness, post-hurdle
## by percent cover

# Summarize mean and standard error for plotting
plot_data <- div_dat %>%
  filter(measurement_type == "Percent") %>%
  group_by(height, site, treatment) %>%
  summarise(
    mean_richness = mean(richness, na.rm = TRUE),
    se_richness = sd(richness, na.rm = TRUE)/sqrt(n()),
    .groups = "drop") |> 
  mutate(height = factor(height, levels = c("Low", "Mid", "High")))

# Plot that shows richness decreasing by zone and site, but relatively unchanged by treatment
ggplot(plot_data, aes(x = height, y = mean_richness, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_richness - se_richness, ymax = mean_richness + se_richness),
                width = 0.2, position = position_dodge(width = 0.8)) +
  facet_wrap(~treatment) +
  labs(
    y = "Species Richness",
    x = "Height",
    title = "Effect of Height and Site on Species Richness for Percent Cover")

# Plot that shows richness unchanged by treatment
ggplot(plot_data, aes(x = height, y = mean_richness, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_richness - se_richness, ymax = mean_richness + se_richness),
                width = 0.2, position = position_dodge(width = 0.8)) +
  facet_wrap(~site) +
  labs(
    y = "Species Richness",
    x = "Height",
    title = "Effect of Treatment on Species Richness for Percent Cover")

## by count

# Summarize mean and standard error for plotting
plot_data <- div_dat %>%
  filter(measurement_type == "Count") %>%
  group_by(height, site, treatment) %>%
  summarise(
    mean_richness = mean(richness, na.rm = TRUE),
    se_richness = sd(richness, na.rm = TRUE)/sqrt(n()),
    .groups = "drop") |> 
  mutate(height = factor(height, levels = c("Low", "Mid", "High")))

# Plot that shows richness decreasing by zone and site, but relatively unchanged by treatment
ggplot(plot_data, aes(x = height, y = mean_richness, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_richness - se_richness, ymax = mean_richness + se_richness),
                width = 0.2, position = position_dodge(width = 0.8)) +
  facet_wrap(~treatment) +
  labs(
    y = "Species Richness",
    x = "Height",
    title = "Effect of Height and Site on Species Richness for Count")

# Plot that shows richness unchanged by treatment
ggplot(plot_data, aes(x = height, y = mean_richness, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_richness - se_richness, ymax = mean_richness + se_richness),
                width = 0.2, position = position_dodge(width = 0.8)) +
  facet_wrap(~site) +
  labs(
    y = "Species Richness",
    x = "Height",
    title = "Effect of Treatment on Species Richness for Count")

