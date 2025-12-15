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

# rm/ spaces in col titles
quad_dat <- quad_dat |> 
  mutate(measurement_type = recode(measurement_type, "Percent Cover" = "Percent"))

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
  ungroup()

# Pivot data for use in vegan::diversity() 
quad_dat <- quad_dat |> 
  pivot_wider(names_from = "species_code",
              values_from = "sum",
              values_fill = 0)

# Drop placeholder spp for quads w/ nothing 
quad_dat <- quad_dat |> 
  select(-c("NOSP", "none_present"))

# Seperate out functional groups
func_dat <- quad_dat |> 
  select(c(1:6,"ASNO", "FUCU", "FUSP", "FUVE", "FUDI", "CHCR", "MAST", "SEBA", "MYED")) |> 
  mutate(temp = rowSums(across(c(ASNO, FUCU, FUSP, FUVE, FUDI, CHCR, MAST, SEBA, MYED)))) |> 
  filter(temp > 0) |> 
  select(-temp)

# Calculate Shannon, InvSimpson, richness, evenness, and rearrange
quad_dat <- quad_dat |> 
  mutate(shannon = diversity(quad_dat[,-c(1:6)], index = "shannon")) # measures the uncertainty in predicting the species identity of an individual chosen at random from the community. A higher value indicates higher diversity (more species and/or more even distribution), and a lower value indicates lower diversity. It equals 0 only when there is no uncertainty—meaning there is only ONE species present.

quad_dat <- quad_dat |> 
  mutate(invsimpson = diversity(quad_dat[,-c(1:6,85)], index = "invsimpson")) # The Inverse Simpson Index is an alternative diversity measure that focuses more on the dominance of species. The lower the value, the greater the diversity. Unlike the Shannon-Wiener index, which treats all species equally, the Inverse Simpson Index weighs more heavily on rare species, but it gives more weight to dominant species.

quad_dat <- quad_dat |> 
  mutate(richness = specnumber(quad_dat[,-c(1:6,85,86)])) 

quad_dat <- quad_dat |> 
  mutate(evenness = shannon/log(richness)) # Pielous eveness how evenly individuals are distributed across the different species in a community. It ranges from 0 (completely uneven) to 1 (perfectly even).# Evenness of 0 occurs when richness is 1 (i.e. log 1 = 0, dividing by 0 = NaN)
  
quad_dat <- quad_dat |> 
select(richness, shannon, evenness, invsimpson, everything())
  
# Set up data for plotting
div_plot_dat <- quad_dat[,-c(11:88)] |> 
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

# Hurdle Model - Odds of being over 0
shannon_percent_hurdle_mod <- glm(shannon > 0 ~ treatment + height + site,
                                  family = binomial(link = "logit"),
                                  data = quad_dat |> filter(measurement_type == "Percent"))

shannon_count_hurdle_mod <- glm(shannon > 0 ~ treatment + height + site,
                                  family = binomial(link = "logit"),
                                  data = quad_dat |> filter(measurement_type == "Count"))

# Models for those over 0
shannon_percent_posthurdle_mod <- lm(sqrt(shannon) ~ treatment + height + site,
                          data = quad_dat |> filter(measurement_type == "Percent"))

shannon_count_posthurdle_mod <- lm(sqrt(shannon) ~ treatment + height + site,
                        data = quad_dat |> filter(measurement_type == "Count"))

# Get model stats
check_model(shannon_percent_hurdle_mod)
check_model(shannon_count_hurdle_mod)
check_model(shannon_percent_posthurdle_mod)
check_model(shannon_count_posthurdle_mod)

# Quad richness
percent_shannon_mod <- lm(sqrt(shannon) ~ treatment + height + site,
                          data = quad_dat |> 
                            filter(measurement_type != "Count"))



















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


# Make sure that for every 0 quad that there is a 0 placeholder for both percent and count 


##
# Functional Group Analysis
##

# Where has what
func_dat <- func_dat |> 
  group_by(season, site, treatment, height) |> 
  mutate(phaeo = rowSums(across(c(ASNO, FUCU, FUSP, FUVE, FUDI)))) |> 
  select(-c(ASNO, FUCU, FUSP, FUVE, FUDI)) |> 
  mutate(rhodo = rowSums(across(c(CHCR, MAST)))) |> 
  select(-c(CHCR, MAST))

# Model of Phaeo abundance
phaeo_mod <- lm(log(phaeo) ~ site + treatment + height,
                data = func_dat |> filter(phaeo > 0))

check_model(phaeo_mod)

hist(func_dat$phaeo[func_dat$phaeo > 0])

# Start old script from here

#' -------------------------------------
#' Look at functional group cover
#' across berms versus controls and
#' quads cover only
#' -------------------------------------
source("scripts/helpers.R")

## Planning
#' Determine what func groups we are interested in: Habitat formers (ASCO, FUSP, CHCR/MAST), 
#' sessile inverts (SEBA, MYED)
#' lm brown_macrophytes ~ site + treatment + zone
#' EXCLUDE MYED from year 1 due to count-PC issue

## load libraries
library(tidyverse)
library(broom)
library(vegan)
library(patchwork)

## load data
quad_s23_dat <- read_csv("data/quads_s23.csv")
quad_f23_dat <- read_csv("data/quads_f23.csv")
quad_s24_dat <- read_csv("data/quads_s24.csv")
quad_s25_dat <- read_csv("data/quads_s25.csv")

# Bind data
func_dat <- rbind(quad_s23_dat |> 
                    group_by(season, site, treatment) |> 
                    reframe(species_code = unique(species_code)),
                  quad_f23_dat |> 
                    group_by(season, site, treatment) |> 
                    reframe(species_code = unique(species_code)), 
                  quad_s24_dat |> 
                    group_by(season, site, treatment) |> 
                    reframe(species_code = unique(species_code)), 
                  quad_s25_dat |> 
                    group_by(season, site, treatment) |> 
                    reframe(species_code = unique(species_code))
) |> 
  filter(species_code != is.na(species_code)) |> 
  filter(!species_code %in% c("none_persent", "NOSP"))

# Drop data that are no longer needed
rm(quad_s23_dat)
rm(quad_f23_dat)
rm(quad_s24_dat)
rm(quad_s25_dat)

# Select relevant spp, and sum across squares
func_dat <- func_dat |>
  filter(species_code %in% c("ASNO", "FUCU", "FUSP", "FUVE", "FUDI", "CHCR", "MAST", "SEBA", "MYED")) |> 
  group_by(site, treatment, height, quadrat, measurement_type, species_code) |> 
  reframe(sum = sum(measurement)) |> 
  ungroup()





##
# Brown Macrophytes
# Analysis of change in ASNO, FUCU, FUSP, FUDI, FUVE
##

phaeo_dat <- percent_dat |> 
  filter(species_code %in% c("ASNO", "FUCU", "FUSP", "FUVE", "FUDI")) |> 
  group_by(site, treatment, height) |> 
  summarise(total_percent = sum(sum))
# Basically only Bayside berm mid and high have any amount of brown macrophytes.

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


#' -------------------------------------
#' Cover Composition
#' Look at cover composition with gllvm
#' across berms versus controls and
#' quads cover only
#' -------------------------------------

