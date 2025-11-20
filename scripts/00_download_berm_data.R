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
# Duplicate function but w/o the site-treatment separating aspect. Data from summer 24 and 25 were entered w/ that already being done
clean_berm_file2 <- function(dat) {
  dat |>
    remove_empty("cols") |>
    clean_names() |>
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
      select(-c(updated_codes,
             notes.y,
             notes.x,
             old_name,
             updated_name))
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
# Quad Data (summer & fall 23; summer 24; summer 25)
###

#summer 23
quads_s23 <- read_sheet("https://docs.google.com/spreadsheets/d/1eOic8zhjAKyBWss7iENRJUrdz52BS-GOp260UR6Pnno/edit?gid=0#gid=0",
                    sheet = "data",
                    na = c("", "NA"))

quads_s23 <- quads_s23 |>
  mutate(season = "summer_23") |>
  clean_berm_file() |> 
  update_spp_codes()|>
  filter(!c(site %in% c("Thompson_Island", #bad install
                        "Trunk",
                        "Powder_Point_Bridge"))) #not yet a berm

# fall 23
quads_f23 <- read_sheet("https://docs.google.com/spreadsheets/d/182VCWn9GilEg4azCROsCkW-00jX9Vh0m2vD0MaJSiyI/edit?usp=sharing",
                        sheet = "Data",
                        na = c("", "NA"))

quads_f23 <- quads_f23 |>
  mutate(season = "fall_23") |>
  clean_berm_file() |> 
  update_spp_codes()|>
  filter(!c(site %in% c("Thompson_Island", #bad install
                        "Trunk",
                        "Powder_Point_Bridge"))) #not yet a berm

# summer 24
quads_s24 <- read_sheet("https://docs.google.com/spreadsheets/d/1kTPiddF99S1HdPHgCMqbkuihAMIQtLUUdf7wVVuIK0g/edit?usp=sharing",
                        sheet = "data",
                        na = c("", "NA"))

quads_s24 <- quads_s24 |>
  mutate(season = "summer_24") |> 
  clean_berm_file2() |> 
  update_spp_codes()|>
  filter(!c(site %in% c("Thompson_Island", #bad install
                        "Trunk",
                        "Powder_Point_Bridge"))) #not yet a berm

# summer 25
quads_s25 <- read_sheet("https://docs.google.com/spreadsheets/d/1qpe9UG-28v-oLscKI0Qw3N-6jHnes01so5xNrJuElQQ/edit?usp=sharing",
                        sheet = "data",
                        na = c("", "NA"))

quads_s25 <- quads_s25 |>
  mutate(season = "summer_25") |>
  clean_berm_file2() |> 
  update_spp_codes()|>
  filter(!c(site %in% c("Thompson_Island", #bad install
                        "Trunk",
                        "Powder_Point_Bridge"))) #not yet a berm


# check
visdat::vis_dat(quads_s23)
skimr::skim(quads_s23)

visdat::vis_dat(quads_f23)
skimr::skim(quads_f23)

visdat::vis_dat(quads_s24)
skimr::skim(quads_s24)

visdat::vis_dat(quads_s25)
skimr::skim(quads_s25)

# write out
write_csv(quads_s23, "data/quads_s23.csv")
write_csv(quads_f23, "data/quads_f23.csv")
write_csv(quads_s24, "data/quads_s24.csv")
write_csv(quads_s25, "data/quads_s25.csv")

###
# Seine Data
###

seine <- read_sheet("https://docs.google.com/spreadsheets/d/1g552Pl-xXLSsi5g-H22jOiRkE4x1T2I-tZuuus_bvZk/edit?gid=0#gid=0",
                    sheet = "Data")

seine <- seine |>
  clean_berm_file() |> 
  update_spp_codes() |>
  filter(!c(site %in% c("Thompson_Island", #bad install
                        "Trunk",
                        "Powder_Point_Bridge"))) #not yet a berm

# write out
write_csv(seine, "data/seine.csv")

###
# Trap Data
###

traps <- read_sheet("https://docs.google.com/spreadsheets/d/1ox_JR305ZaYIgTVGWsVf4pn01OXQ5F7EJbdWPZvxVyk/edit?gid=0#gid=0",
                    sheet = "Data")

traps <- traps |>
  clean_berm_file() |> 
  update_spp_codes() |>
  separate(trap_name_number,
           c("trap_name", "trap_number"),
           sep = "_") |>
  filter(!c(site %in% c("Thompson_Island", #bad install
                        "Trunk",
                        "Powder_Point_Bridge"))) #not yet a berm


# write out
write_csv(traps, "data/traps.csv")

###
# BRIV Data
###
# drive_download("https://docs.google.com/spreadsheets/d/1Th0rDY5iQckVB4O9QrHg2y0xfj-e-3mz/edit?usp=sharing&ouid=105030083689633832675&rtpof=true&sd=true",
#                "data/briv_raw.xlsx")

brivs <- read_excel("data/briv_raw.xlsx",
                    sheet = "SLL_BRUV_BRIV_Data") |>
  remove_empty("cols") |>
  clean_names() |>
  select(-c(bruv_or_briv,
            depth_ft, 
            tidal_cycle,
            bait)) |>
  separate_wider_regex(cols = site,
                       patterns = 
                         c(site = ".*", "_", 
                           treatment = ".*")) |>
  mutate(first_predator = ifelse(first_predator == "CRNK",
                                 "UCRA",
                                 first_predator)) |>
  filter(!c(site %in% c("Thompson_Island", #bad install
                        "Trunk",
                        "Powder_Point_Bridge"))) #not yet a berm


write_csv(brivs, "data/brivs.csv")

###
# FARM Data
###

farms <- read_sheet("https://docs.google.com/spreadsheets/d/11H5kCxaMysnp7sXp_dN-k1DT-MCO58loKfJrMLUGeh4/edit?gid=0#gid=0",
                    sheet = "Data") |>
  clean_berm_file() |> 
  update_spp_codes()

write_csv(farms, "data/farms.csv")

