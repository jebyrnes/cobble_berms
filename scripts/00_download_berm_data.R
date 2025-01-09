#' -------------------------------------
#' Download data files from the google drive
#' -------------------------------------

# load libraries for accessing google drive
# writing out CSVs
# and doing some light data cleaning
pacman::p_load(googlesheets4,
               googledrive,
               readr,
               readxl,
               dplyr,
               janitor,
               tidyr)

# authorize
gs4_auth()

###
# Helper functions
###

# a function that takes a data frame
# removes empty columns, cleans up names
# and then splits the site into site and treatment
# as well as filling in a marker for no species
clean_berm_file <- function(dat) {
  dat |>
    remove_empty("cols") |>
    clean_names() |>
    separate_wider_regex(cols = site,
                         patterns = c(site = ".*", "_", treatment = ".*")) |>
    mutate(species_code = ifelse(is.na(species_code), "NOSP", species_code))
  
}

update_spp_codes <- function(dat) {
  dat |> 
      filter(!(species_code %in% c("DIAT", # dont consistently assess diatoms or micro orgs
                                   "MCRO", # see above
                                   "UES"))) |>  # don't need unknown egg sack
      left_join(y = sp_replace_codes, 
                by = join_by(species_code == old_codes)) |> 
      mutate(species_code = ifelse(test = !is.na(updated_codes), 
                                   yes = updated_codes, 
                                   no = species_code)) |> 
      select(-updated_codes)
}

###
# Species lists
###
#drive_download("https://docs.google.com/spreadsheets/d/1xHQB7VxzH33HgSl9hZfqV9YAtrns48Rq/edit?usp=sharing&ouid=105030083689633832675&rtpof=true&sd=true",
#                                 path = "data/species_list.xlsx")

master_sp_list <- read_excel("data/species_list.xlsx",
                             sheet = "Master_Species_List") |> 
  clean_names()

sp_replace_codes <- read_excel("data/species_list.xlsx",
                               sheet = "Replaced_Codes") |> 
  clean_names()

write_csv(master_sp_list, "data/master_sp_list.csv")

###
# Quad Data
###
quads <- read_sheet("https://docs.google.com/spreadsheets/d/1eOic8zhjAKyBWss7iENRJUrdz52BS-GOp260UR6Pnno/edit?gid=0#gid=0",
                    sheet = "data")

quads <- quads |>
  clean_berm_file() |> 
  update_spp_codes()


# check
visdat::vis_dat(quads)
skimr::skim(quads)

# write out
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

brivs <- read_sheet("https://docs.google.com/spreadsheets/d/1Th0rDY5iQckVB4O9QrHg2y0xfj-e-3mz/edit?gid=1022311473#gid=1022311473",
                    sheet = "SLL_BRUV_BRIV_Data")

###
# FARM Data
###
