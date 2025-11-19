#' -------------------------------------
#' Look at shannon diversity
#' across berms versus controls and
#' quads %, quads #, seines, traps, and brivs 
#' separately
#' -------------------------------------
source("scripts/helpers.R")

## load libraries
library(tidyverse)
library(broom)
library(vegan)
library(patchwork)

## load data
quad_s23_dat <- read_csv("data/quads.csv")

## Planning
#' Process all the below into this new rproject(its unedited since I was working w/o git)
#' Give it essentially the same treatment as the quads (i.e. dont double count lower levels)
#' Calc div by taxonomic levels
#' Analyze the same as richness, essentially ttest


# Sum across squares for percent data
percent_dat <- quad_s23_dat |> 
  filter(measurement_type != "Count") |> 
  group_by(site, treatment, height, quadrat, measurement_type, species_code) |> 
  reframe(sum = sum(measurement)) |> 
  filter(species_code != "NOSP") |> 
  ungroup()

# Sum across squares for count data
count_dat <- quad_s23_dat |> 
  filter(measurement_type == "Count") |> 
  group_by(site, treatment, height, quadrat, measurement_type, species_code) |> 
  reframe(sum = sum(measurement)) |> 
  filter(species_code != "NOSP") |> 
  ungroup()

# Pivot data for use in vegan::diversity() 
percent_dat <- percent_dat |> 
  pivot_wider(names_from = "species_code",
              values_from = "sum",
              values_fill = 0)

count_dat <- count_dat |> 
  pivot_wider(names_from = "species_code",
              values_from = "sum",
              values_fill = 0)

# Calc Diversity Indicies

# Shannon and InvSimpson
percent_dat <- percent_dat |> 
  mutate(shannon = diversity(percent_dat[,-c(1:5)], index = "shannon")) |> 
  mutate(invsimpson = diversity(percent_dat[,-c(1:5)], index = "invsimpson"))

count_dat <- count_dat |> 
  mutate(shannon = diversity(count_dat[,-c(1:5)], index = "shannon")) |> 
  mutate(invsimpson = diversity(count_dat[,-c(1:5)], index = "invsimpson"))

# Hurdle Mod Analysis

# Odds of being over 0
percent_dat_hurdle_mod <- glm(shannon > 0 ~ treatment + height + site,
                             family = binomial(link = "logit"),
                             data = percent_dat)

count_dat_hurdle_mod <- glm(shannon > 0 ~ treatment + height + site,
                             family = binomial(link = "logit"),
                             data = count_dat)

# Model for those over 0
percent_dat_mod <- lm(sqrt(shannon) ~ treatment + height + site,
                     data = percent_dat |> 
                       filter(shannon > 0))

count_dat_mod <- lm(sqrt(shannon) ~ treatment + height + site,
                     data = count_dat |> 
                       filter(shannon > 0))

# viz

## by count

jitter_width <- 0.1
jitter_positions <- position_dodge(width = jitter_width)

count_em <- tidy(emmeans(count_dat_mod, ~height + treatment)) |> 
  mutate(shannon = estimate^2,
         lower.CL = ((estimate - std.error)^2),
         upper.CL = ((estimate + std.error)^2))

shannon_count_plot <- ggplot(count_em,
                             mapping = aes(x = treatment, y = shannon,
                                           ymin = lower.CL, ymax = upper.CL,
                                           color = height,
                                           group = height)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), 
                  position = jitter_positions) +
  geom_line(linewidth = 1, position = jitter_positions) + 
  ggtitle("Count")  +
  ylab("Shannon Diversity") +
  xlab(element_blank()) +
  theme(legend.position="none") +
  scale_color_manual(values = c("Low" = "#af8dc3", 
                                "Mid" = "#bdbdbd", 
                                "High" = "#7fbf7b"))


## by percent cover

percent_em <- tidy(emmeans(percent_dat_mod, ~height + treatment)) |> 
  mutate(shannon = estimate^2,
         lower.CL = ((estimate - std.error)^2),
         upper.CL = ((estimate + std.error)^2))

shannon_cover_plot <- ggplot(percent_em,
                             mapping = aes(x = treatment, 
                                           y = shannon,
                                           ymin = lower.CL, 
                                           ymax = upper.CL,
                                           color = height,
                                           group = height)) +
  geom_pointrange(aes(ymin = lower.CL, 
                      ymax = upper.CL), 
                  position = jitter_positions) +
  geom_line(linewidth = 1, position = jitter_positions) + 
  ggtitle("Percent Cover") +
  ylab("Shannon Diversity") +
  xlab(element_blank()) +
  scale_color_manual(values = c("Low" = "#af8dc3", 
                                "Mid" = "#bdbdbd", 
                                "High" = "#7fbf7b")) +
  theme(legend.position="none")

shannon_plot <- shannon_count_plot + shannon_cover_plot +
  plot_annotation(tag_levels = 'I') 

ggsave("figures/shannon_plot.jpg", 
       shannon_plot,
       dpi = 600)

# Get model stats

summary(percent_dat_hurdle_mod)
summary(percent_dat_mod)
summary(count_dat_hurdle_mod)
summary(count_dat_mod)
