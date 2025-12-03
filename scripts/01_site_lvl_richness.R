#' -------------------------------------
#' Calculate richness at the site level.
#' Analyze with a t.test
#' -------------------------------------
source("scripts/helpers.R")

##
# Set-up
##

#load libraries
library(tidyverse)
library(broom)
library(knitr)
library(kableExtra)

# load data, bind them together, and add taxonomic data
quad_s23_dat <- read_csv("data/quads_s23.csv")
quad_f23_dat <- read_csv("data/quads_f23.csv")
quad_s24_dat <- read_csv("data/quads_s24.csv")
quad_s25_dat <- read_csv("data/quads_s25.csv")
seine_dat <- read_csv("data/seine.csv")
trap_dat <- read_csv("data/traps.csv")
farm_dat <- read_csv("data/farms.csv")
briv_dat <- read_csv("data/brivs.csv")
msl <- read_csv("data/master_sp_list.csv")

# Drop uneeded cols in msl and rename col
msl <- select(msl, 1:10) |> 
  rename(species_code = code)

# bind data together, and drop the two NAs in Bayside Berm
richness_dat <- rbind(quad_s23_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)),
                      quad_f23_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)), 
                      quad_s24_dat |> 
                        group_by(site, treatment) |> 
                        reframe(species_code = unique(species_code)), 
                      quad_s25_dat |> 
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
) |> 
  filter(species_code != is.na(species_code)) |> 
  filter(species_code != "NOSP") |> 
  filter(species_code != "none_present") |> 
  filter(species_code != "UCRA")

# attach taxonomic data
richness_dat <- left_join(richness_dat, 
                          msl, 
                          "species_code")

## Drop data that are no longer needed
rm(quad_s23_dat)
rm(quad_f23_dat)
rm(quad_s24_dat)
rm(quad_s25_dat)
rm(seine_dat)
rm(trap_dat)
rm(farm_dat)
rm(briv_dat)
rm(msl)

# Remove duplicate rows
richness_dat <- unique(richness_dat)

##
# Verification: Drop higher order codes
##

#Bayside Berm
# richness_dat |> 
#   filter(site == "Bayside") |> 
#   filter(treatment == "Berm") |> 
#   View()
# Codes to drop: UFAC, FUCU, UBBL, URAL

# Bayside Control
# richness_dat |> 
#   filter(site == "Bayside") |> 
#   filter(treatment == "Control") |> 
#   View()
# Codes to drop: UFAC, UFIS, FUCU, UBBL, URFI, URBL, UCCA

# Coughlin Berm
# richness_dat |> 
#   filter(site == "Coughlin") |> 
#   filter(treatment == "Berm") |> 
#   View()
# Codes to drop: UNAL, UFAC, FUCU, UBAL, UBBL, UBFI, URCB, UCL, URFC, UCCA, UBFC

# Coughlin Control
# richness_dat |> 
#   filter(site == "Coughlin") |> 
#   filter(treatment == "Control") |> 
#   View()
# Codes to drop: UNAL, UFAC, PHSP, UCCA

# Duxbury_Beach Berm
# richness_dat |> 
#   filter(site == "Duxbury_Beach") |> 
#   filter(treatment == "Berm") |> 
#   View()
# Codes to drop: UFAC, PASP, URFC

# Duxbury_Beach Control
# richness_dat |> 
#   filter(site == "Duxbury_Beach") |> 
#   filter(treatment == "Control") |> 
#   View()
# Codes to drop: UFAC, URFC

# Filter our codes to drop
richness_dat <- richness_dat |> 
  filter(
    !(site == "Bayside" & treatment == "Berm" & species_code %in% c("UFAC", "FUCU", "UBBL", "URAL")) &
    !(site == "Bayside" & treatment == "Control" & species_code %in% c("UFAC", "UFIS", "FUCU", "UBBL", "URFI", "URBL", "UCCA")) &
    !(site == "Coughlin" & treatment == "Berm" & species_code %in% c("UNAL", "UFAC", "FUCU", "UBAL", "UBBL", "UBFI", "URCB", "UCL", "URFC", "UCCA", "UBFC")) &
    !(site == "Coughlin" & treatment == "Control" & species_code %in% c("UNAL", "UFAC", "PHSP", "UCCA")) &
    !(site == "Duxbury_Beach" & treatment == "Berm" & species_code %in% c("UFAC", "PASP", "URFC")) &
    !(site == "Duxbury_Beach" & treatment == "Control" & species_code %in% c("UFAC", "URFC"))
  )

##
# Calculate Richness
##

# Calc richness per site_treatment
richness_dat <- richness_dat |>
  group_by(site, treatment) |>
  reframe(richness = n())

##
# Analyze Richness
## 


## Analysis of Species Richness

# Set up data
ttest_dat <- richness_dat |> 
  pivot_wider(names_from = "treatment",
              values_from = "richness")
# Run ttest 
richness_ttest <- t.test(ttest_dat$Berm, ttest_dat$Control, paired = TRUE) |> 
tidy()

# estimate statistic p.value parameter conf.low conf.high method        alternative
# <dbl>     <dbl>   <dbl>     <dbl>    <dbl>     <dbl> <chr>         <chr>      
#   1     4.33      2.98  0.0964         2    -1.92      10.6 Paired t-test two.sided  

# Unsurprisingly not significant w/ large p value. Tiny sample! 


##
# Visualize
##

# Create boxplots
richness_plot <- ggplot(richness_dat, 
                        mapping = 
         aes(x = treatment, 
             y = richness, 
             group = site,
             fill = treatment,
             shape = site)) +
  geom_line(color = "black") +
  geom_point(size = 3, 
             color = "black", 
             stroke = 1.2) +
  labs(title = "Site Level Species Richness by Treatment",
       x = "Treatment",
       y = "Species Richness") +
  scale_fill_manual(values = c("Berm" = "#af8dc3", "Control" = "#7fbf7b")) + 
  scale_shape_manual(values = c(21, 22, 24)) +
  guides(fill = "none")

# write out richness figure
ggsave("figures/richness_sitelvl_plot.jpg", 
       richness_plot,
       dpi = 600)


richness_table <- ttest_dat |> 
  kable("html", caption = "Site Level Species Richness by Treatment") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) %>%
  row_spec(0, bold = TRUE, background = "#D3D3D3")

save_kable(richness_table, file = "figures/richness_sitelvl_table.html")
