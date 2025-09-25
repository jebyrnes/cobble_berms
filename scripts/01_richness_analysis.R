#' -------------------------------------
#' Look at total species richness
#' across berms versus controls and
#' quads (all), seines, traps, and brivs combined
#' -------------------------------------
source("scripts/helpers.R")

## load libraries

snail_dat <- read_csv("data/quads.csv") |>
  filter(measurement_type == "Count") |>
  filter(species_code == "LILI") |>
  group_by(site, treatment, quadrat, height) |>
  summarize(snails = sum(measurement, na.rm = TRUE)) |>
  ungroup() |> fix_data_for_plots()

## load data then bind them together
quad_dat <- read_csv("data/quads.csv")
seine_dat <- read_csv("data/seine.csv")
trap_dat <- read_csv("data/traps.csv")
farm_dat <- read_csv("data/farms.csv")
briv_dat <- read_csv("data/brivs.csv")

richness_dat <- rbind(quad_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)), 
                      seine_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)),
                      trap_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)),
                      briv_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(first_predator)),
                      briv_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(bait_consumer))
)

## Calculate Richness
richness_dat <- richness_dat |>
  unique() |>
  group_by(site, treatment) |>
  reframe(richness = n())

## Visualize
ggplot(data = richness_dat,
       mapping = aes(x = treatment,
                     y = richness,
                     fill = treatment)) +
  geom_boxplot() +
  ylab("Species Richness") +
  xlab("Treatment") +
  scale_fill_manual(values = c("Berm" = "#af8dc3", "Control" = "#7fbf7b"))
