#' -------------------------------------
#' Download data files from the google drive
#' -------------------------------------

# load libraries for accessing google drive
# writing out CSVs
# and doing some light data cleaning
pacman::p_load(googlesheets4,
               readr,
               dplyr,
               janitor,
               tidyr)

# authorize
gs4_auth()

###
# Species lists
###

###
# Quad Data
###
quads <- read_sheet("https://docs.google.com/spreadsheets/d/1eOic8zhjAKyBWss7iENRJUrdz52BS-GOp260UR6Pnno/edit?gid=0#gid=0",
                    sheet = "data")
quads <- quads |>
  remove_empty("cols") |> 
  clean_names() |> 
  separate_wider_regex(cols = site, 
                       patterns = c(site = ".*", "_", 
                                    treatment = ".*")) |> 
  mutate(species_code = ifelse(is.na(species_code), "NOSP", species_code))

visdat::vis_dat(quads)
skimr::skim(quads)

write_csv(quads, "data/quads.csv")



###
# Seine Data
###

###
# Trap Data
###

###
# BRIV Data
###

###
# FARM Data
###
