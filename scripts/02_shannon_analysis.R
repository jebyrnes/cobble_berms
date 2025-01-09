#' -------------------------------------
#' Look at shannon diversity
#' across berms versus controls and
#' quads %, quads #, seines, traps, and brivs 
#' separately
#' -------------------------------------
source("scripts/helpers.R")


# Sum across squares
div_qp_dat <- quads |> 
  group_by(site, treatment, height, quadrat, measurement_type, species_code) |> 
  reframe(sum = sum(measurement)) |> 
  ungroup()

# Update spp codes and ditch unecessary ones
diversity_quad_L23 <- diversity_quad_L23 |> 
  mutate(species_code = toupper(species_code)) |> 
  left_join(y = replacement_codes, 
            by = join_by(species_code == old_code)) |> 
  mutate(species_code = ifelse(test = !is.na(updated_code), 
                               yes = updated_code, 
                               no = species_code)) |> 
  select(-updated_code)

diversity_quad_L23 <- diversity_quad_L23 |> 
  filter(!(species_code %in% c("DIAT", # dont consistently assess diatoms or micro orgs
                               "MCRO", # see above
                               "UES", # don't need unknown egg sack
                               "NONE_PRESENT"))) |>  # Not needed for diversity
  filter(!is.na(site))


# Pivot bio data for use in vegan::diversity() 
diversity_quad_L23 <- diversity_quad_L23 |> 
  pivot_wider(names_from = "species_code",
              values_from = "sum",
              values_fill = 0)

##
# Calc Diversity Indicies
##

# Shannon and InvSimpson
diversity_quad_L23 <- diversity_quad_L23 |> 
  mutate(shannon = diversity(diversity_quad_L23[,-c(1:6)], index = "shannon")) |> 
  mutate(invsimpson = diversity(diversity_quad_L23[,-c(1:6, 41)], index = "invsimpson")) |> 
  select(-c(season, 7:40))