#' -------------------------------------
#' Script for other exploratory biodiversity analyses
#' Keep scripts 1 and 2 pure
#' -------------------------------------
source("scripts/helpers.R")

##
#  Set-up ----------------------------------
##

# Run after script 2


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

##
# NMDS ----------------------------------------------
##
# 
# 
# comm_percent <- quad_dat |>
#   filter(measurement_type == "Percent") |>
#   filter(height == "Low") |>
#   select(7:84) |>
#   mutate(across(everything(), ~ as.numeric(.x))) |>
#   mutate(across(everything(), ~ replace_na(.x, 0)))
# 
# meta_percent <- quad_dat |>
#   filter(measurement_type == "Percent")|>
#   filter(height == "Low") |>
#   select(site, treatment, height)
# 
# meta_percent <- meta_percent[
#   rowSums(comm_percent) > 0,
# ]
# 
# comm_percent <- comm_percent[
#   rowSums(comm_percent) > 0,
# ] #remove all zeros
# 
# nmds_percent <- metaMDS(
#   comm_percent,
#   distance = "bray",
#   k = 2,
#   trymax = 100,
#   autotransform = FALSE
# )
# nmds_plot <- scores(nmds_percent, display = "sites") |>
#   as.data.frame() |>
#   bind_cols(meta_percent)
# 
# ggplot(nmds_plot, aes(NMDS1, NMDS2)) +
#   geom_point(
#     aes(color = site, shape = treatment),
#     size = 4
#   ) +
#   stat_ellipse(
#     aes(
#       color = site,
#       group = interaction(site, treatment)
#     ),
#     linetype = 2,
#     linewidth = 0.8
#   ) +
#   theme_classic() +
#   labs(
#     x = "NMDS1",
#     y = "NMDS2",
#     color = "Site",
#     shape = "Treatment",
#     title = "NMDS of Percent Cover"
#   )
# 
