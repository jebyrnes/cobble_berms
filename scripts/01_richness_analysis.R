#' -------------------------------------
#' Look at total species richness
#' across berms versus controls and
#' quads (all), seines, traps, and brivs combined
#' -------------------------------------
source("scripts/helpers.R")

## load libraries


## bind all datasets together
richness_dat <- rbind(quads |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)), 
                      seine |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)),
                      traps |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)),
                      brivs |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(first_predator)),
                      brivs |> 
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
