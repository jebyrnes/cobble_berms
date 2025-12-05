#' -------------------------------------
#' Look at quadrat level richness, shannon diversity, and eveness
#' across berms versus controls and
#' for percent cover and count seperatly
#' -------------------------------------
source("scripts/helpers.R")

##
# Set-up
##

# load libraries
library(tidyverse)
library(easystats)
# library(broom)
library(vegan)
# library(patchwork)
# library(emmeans)

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

# Drop data that are no longer needed
rm(quad_s23_dat)
rm(quad_f23_dat)
rm(quad_s24_dat)
rm(quad_s25_dat)

## Verify correct amount of quads done
# quad_dat |> 
#   group_by(site, treatment, season, height) |> 
#   reframe(quad_count = n_distinct(quadrat)) |> 
#   View()
  
# Problem data: 
# Coughlin-Control-fall_23-High-2 quads only
# Coughlin-Control-fall_23-Mid-4
# Bayside-Berm-summer_25-High-5
# Coughlin-Berm-summer_25-Low-5


##
# Calculate Richness, biodiversity, eveness
##

# Sum to the quad level
quad_dat <- quad_dat |> 
  group_by(season, site, treatment, height, quadrat, measurement_type, species_code) |> 
  reframe(sum = sum(measurement)) |> 
  filter(species_code != "NOSP") |> 
  filter(species_code != "none_present") |> 
  ungroup()

# Pull out quad richness w/o measurement_type
quad_richness <- quad_dat |> 
  group_by(season, site, treatment, height, quadrat) |> 
  reframe(richness = n_distinct(species_code))

# Pivot data for use in vegan::diversity() 
quad_dat <- quad_dat |> 
  pivot_wider(names_from = "species_code",
              values_from = "sum",
              values_fill = 0)

# Calculate Shannon, InvSimpson, richness, evenness, and rearrange
quad_dat <- quad_dat |> 
  mutate(shannon = diversity(quad_dat[,-c(1:6)], index = "shannon")) # measures the uncertainty in predicting the species identity of an individual chosen at random from the community. A higher value indicates higher diversity (more species and/or more even distribution), and a lower value indicates lower diversity. It equals 0 only when there is no uncertainty—meaning there is only ONE species present.

quad_dat <- quad_dat |> 
  mutate(invsimpson = diversity(quad_dat[,-c(1:6,85)], index = "invsimpson")) # The Inverse Simpson Index is an alternative diversity measure that focuses more on the dominance of species. The lower the value, the greater the diversity. Unlike the Shannon-Wiener index, which treats all species equally, the Inverse Simpson Index weighs more heavily on rare species, but it gives more weight to dominant species.

quad_dat <- quad_dat |> mutate(richness = specnumber(quad_dat[,-c(1:6),85,86])) |> 
  mutate(evenness = shannon/log(richness))  # Pielous eveness how evenly individuals are distributed across the different species in a community. It ranges from 0 (completely uneven) to 1 (perfectly even).
  
quad_dat <- quad_dat |> 
select(richness, shannon, evenness, invsimpson, everything())
  
## Something is wrong here. Its detecting richness of 2 when there is only a single species present

  
##
# Visualize
##

# Set up data for plotting
div_plot_dat <- quad_dat[,-c(11:88)] |> 
  pivot_longer(cols = richness:invsimpson, 
               values_to = "value", 
               names_to = "index")


ggplot(div_plot_dat,
       aes(x = treatment,
           y = value)) +
  geom_boxplot() +
  facet_grid(rows = vars(height), cols = vars(index))



##
# Analysis
##

# Models for those over 0
percent_shannon_mod <- lm(sqrt(shannon) ~ treatment + height + site,
                      data = quad_dat |> 
                        filter(measurement_type != "Count"))

count_shannon_mod <- lm(sqrt(shannon) ~ treatment + height + site,
                    data = quad_dat |> 
                      filter(measurement_type == "Count"))



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








##
# Analysis
##


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



# Get model stats

summary(percent_dat_hurdle_mod)
summary(percent_dat_mod)
summary(count_dat_hurdle_mod)
summary(count_dat_mod)
