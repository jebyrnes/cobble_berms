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

# Drop placeholder spp for quads w/ nothing 
quad_dat <- quad_dat |> 
  select(-c("NOSP"))

##
#  Richness and Shannon Biodiversity analysis ----------------------------------
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


ggplot(div_plot_dat |> filter(measurement_type == "Count"),
       aes(x = treatment,
           y = value)) +
  geom_violin() +
  facet_grid(rows = vars(index), 
             cols = vars(height),
             scales = "free_y")

ggplot(div_plot_dat |> filter(measurement_type == "Percent"),
       aes(x = treatment,
           y = value)) +
  geom_violin() +
  facet_grid(rows = vars(index), 
             cols = vars(height),
             scales = "free_y")

##
# Analysis
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
shannon_percent_posthurdle_mod <- lm(sqrt(shannon) ~ treatment + height + site + treatment:height, # note that shannon at richness of 1 is 0, no diversity, only 1 sp present which is the same as richness of 0
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

contrast(shannon_percent_hurdle_em, method = "pairwise") # Significant effect of treatment in the low and mid

# Post-hurdle (>0) for shannon percent
Anova(shannon_percent_posthurdle_mod) |> tidy() # Percent Shannon  differs significantly across heights and sites, but is not strongly affected by treatment, and the interaction between treatment and height is also not significant.
tidy(shannon_percent_posthurdle_mod) |>  
  mutate(estimate_original = sign(estimate) * estimate^2,
         std.error_original = 2 * abs(estimate) * std.error) |> 
  mutate(lower.CL = ((estimate_original - std.error_original)^2),
         upper.CL = ((estimate_original + std.error_original)^2))

# Height matters a lot: Low  plots have the highest Shannon diversity (0.768 and 0.535).Mid heights are intermediate (0.511 and 0.375). High height plots have the lowest diversity (0.070 and 0.049).
# Treatment effect is small: Within each height category, Berm has slightly higher Shannon diversity than Control.
# So these differences are statistically robust (all p values are << .001.)
# Differences between heights are much larger than differences between treatments. Matches ANOVA: height is the main driver, treatment is minor.

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

## richness
# Hurdle Model - Odds of being over 0
richness_percent_hurdle_mod <- glm(richness > 0 ~ treatment + height + site + treatment:height,
                                   family = binomial(link = "logit"),
                                   data = div_dat |> filter(measurement_type == "Percent"))

richness_count_hurdle_mod <- glm(richness > 0 ~ treatment + height + site + treatment:height,
                                 family = binomial(link = "logit"),
                                 data = div_dat |> filter(measurement_type == "Count"))

# Models for those over 0
richness_percent_posthurdle_mod <- lm(sqrt(richness) ~ treatment + height + site + treatment:height,
                                      data = div_dat |> filter(measurement_type == "Percent", richness > 0))

richness_count_posthurdle_mod <- lm(sqrt(richness) ~ treatment + height + site + treatment:height,
                                    data = div_dat |> filter(measurement_type == "Count", richness > 0))

# Get model stats
# check_model(richness_percent_hurdle_mod)
# check_model(richness_count_hurdle_mod)
# check_model(richness_percent_posthurdle_mod)
# check_model(richness_count_posthurdle_mod) 

## So what does this tell us
# Hurdle (0,1) for richness percent
Anova(richness_percent_hurdle_mod) |> tidy() # all sig, height > treatment > site
tidy(richness_percent_hurdle_mod)

richness_percent_hurdle_em <- emmeans(richness_percent_hurdle_mod,
                                     ~treatment | height)

contrast(richness_percent_hurdle_em, method = "pairwise") # Significant effect of treatment in all (results identical?)


# Hurdle (0,1) for richness count
Anova(richness_count_hurdle_mod) |> tidy() # all sig, height > treatment > site
tidy(richness_count_hurdle_mod) # effect of treatment and height evident

richness_count_hurdle_em <- emmeans(richness_count_hurdle_mod,
                                   ~treatment | height)

contrast(richness_count_hurdle_em, method = "pairwise") #Significant effect of treatment in all (results identical?)

# Post-hurdle (>0) for richness percent
Anova(richness_percent_posthurdle_mod) |> tidy() # effect of height and site, not treatment 
tidy(richness_percent_posthurdle_mod) |>  
  mutate(estimate_original = sign(estimate) * estimate^2,
         std.error_original = 2 * abs(estimate) * std.error)|> 
  mutate(lower.CL = ((estimate_original - std.error_original)^2),
         upper.CL = ((estimate_original + std.error_original)^2))

richness_percent_posthurdle_em <- emmeans(richness_percent_posthurdle_mod,
                                         ~treatment | height)

contrast(richness_percent_posthurdle_em, method = "pairwise") #no sig

# Post-hurdle (>0) for richness count
Anova(richness_count_posthurdle_mod) |> tidy() # height and site

tidy(richness_count_posthurdle_mod) |>  # low zone and duxbury are particularlly impactful
  mutate(estimate_original = sign(estimate) * estimate^2, # note that the negative estimates are lost after the back transformation, so i've multiplioed them by the sign to return the original sign
         std.error_original = 2 * abs(estimate) * std.error) |> 
  mutate(lower.CL = ((estimate_original - std.error_original)^2),
         upper.CL = ((estimate_original + std.error_original)^2))


richness_count_posthurdle_em <- emmeans(richness_count_posthurdle_mod,
                                       ~treatment + height)

contrast(richness_count_posthurdle_em, method = "pairwise") # no significance




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


##
#  Functional Group Analysis ----------------------------------------------
##

#' -------------------------------------
#' Look at functional group cover
#' across berms versus controls and
#' quads cover only
#' -------------------------------------


# Seperate out functional groups
func_dat <- quad_dat |> 
  select(c(1:6,"ASNO", "FUCU", "FUSP", "FUVE", "FUDI", "CHCR", "MAST", "SEBA", "MYED")) 

# Code to drop rows w/o any of the selected groups. Removed to retain 0 quads
  # mutate(temp = rowSums(across(c(ASNO, FUCU, FUSP, FUVE, FUDI, CHCR, MAST, SEBA, MYED)))) |> 
  # filter(temp > 0) |> 
  # select(-temp)


# Where has what
func_dat <- func_dat |> 
  group_by(season, site, treatment, height) |> 
  mutate(phaeo = rowSums(across(c(ASNO, FUCU, FUSP, FUVE, FUDI)))) |> 
  mutate(rhodo = rowSums(across(c(CHCR, MAST)))) |> 
  select(-c(ASNO, FUCU, FUSP, FUVE, FUDI, CHCR, MAST)) |> 
  mutate(MYED_pres_abs = ifelse(test = MYED > 0, yes = 1, no = 0)) |> 
  select(-MYED) # MYED has some erroneous count data, and will be analyzed as presence absence instead of abundance

func_dat |> 
  group_by(season, site, treatment, height) |> 
  summarise(across(c(SEBA, phaeo, rhodo, MYED_pres_abs), sum, na.rm = TRUE))

# Phaeo

n <- sum(func_dat$phaeo > 0, na.rm = TRUE)
ggplot(func_dat |> filter(phaeo > 0), aes(phaeo)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "black") +
  labs(x = "Value", y = "Count", title = paste("Phaeo abundance > 0 (n =", n, " / ", n_total, ")"))

n <- sum(func_dat$rhodo > 0, na.rm = TRUE)
ggplot(func_dat |> filter(rhodo > 0), aes(rhodo)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "black") +
  labs(x = "Value", y = "Count", title = paste("Rhodo abundance > 0 (n =", n, " / ", n_total, ")"))

n <- sum(func_dat$SEBA > 0, na.rm = TRUE)
ggplot(func_dat |> filter(SEBA > 0), aes(SEBA)) +
  geom_histogram(binwidth = 1, boundary = 0, color = "black") +
  labs(x = "Value", y = "Count", title = paste("SEBA abundance > 0 (n =", n, " / ", n_total, ")"))

n <- sum(func_dat$MYED_pres_abs, na.rm = TRUE)
func_dat |> 
  group_by(treatment) |> 
  summarise(sum(MYED_pres_abs))
# MYED present in 58/424 quads. 30 on the berm, 28 on the control


# Model of Phaeo abundance
phaeo_mod <- lm(log(phaeo) ~ site + treatment + height,
                data = func_dat |> filter(phaeo > 0))

# check_model(phaeo_mod)



##
# Brown Macrophytes
# Analysis of change in ASNO, FUCU, FUSP, FUDI, FUVE
##


##
# Red Macrophytes
# Analysis of change in CHCR, MAST
##

rhodo_dat <- percent_dat |> 
  filter(species_code %in% c("CHCR", "MAST")) |> 
  group_by(site, treatment, height) |> 
  summarise(total_percent = sum(sum))

# Only found in low zone at bay berm, and cough berm and control. Essentially the same

##
# Sessile Inverts
# Analysis of change in SEBA, MYED
##

invert_dat <- percent_dat |> 
  filter(species_code %in% c("SEBA", "MYED")) |> 
  group_by(site, treatment, height, species_code) |> 
  summarise(total_percent = sum(sum))

seba_plot <- ggplot(percent_dat |> filter(species_code == "SEBA"),
                    mapping = 
                      aes(x = treatment, 
                          y = sum, 
                          fill = treatment)) +
  facet_grid(rows = vars(height), cols = vars(site)) +
  geom_point(alpha = .6, position = "jitter") +
  geom_boxplot(alpha = .6)

# Stats
percent_dat |> 
  filter(species_code == "SEBA") |>
  pivot_wider(names_from = treatment)
t.test(x, y, paired = TRUE)

## Results - no difference in SEBA btw sites. Better at bayside, worse at coughlin, 
## the same at duxbury, and non at powderpoint.

# Analysis
# Take code from above, but dont combine through the quads. Then do a t-test comparing sites

# Option for completing percent cover and coutns for all sites and treatments to have 0
## Fix up snail data to fill in 0 for quads where there
## were no snails
# snail_dat <- snail_dat |>
#   complete(site, treatment, nesting(quadrat, height),
#            fill = list(snails = 0)) |>
#   mutate(are_there_snails = as.numeric(snails !=0))

##
#  Cover Composition ----------------------------------------------
##

#' Cover Composition
#' Look at cover composition with gllvm
#' across berms versus controls and
#' quads cover only

