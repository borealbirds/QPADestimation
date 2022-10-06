# ---
# title: "QPAD estimation - clean bird data"
# author: "Elly Knight"
# created: "September 22, 2022"
# ---

library(tidyverse) #basic data wrangling
library(lubridate) #date and time wrangling

#Load previous dataset to look at it----
#load into separate environment to avoid overwriting things
# e <- new.env()
# load("data/new_offset_data_package_2017-03-01.Rdata", envir = e)
#
# #look at data structure
# str(e$pc)

#1. Load in raw  & visit dataset----
load("data/wildtrax_data_2022-07-24.Rdata")
load("data/visit_data_2022-07-24.Rdata")

#2. Filter, tidy, create foreign key for visit table----
#Remove surveys with no location
#Remove outliers for day of year (use 99% quantile)
#Ensure there's visit data
bird <- raw %>%
    mutate(datetime = ymd_hms(date),
           julian = yday(datetime),
           abundance = as.numeric(abundance),
           id = paste(location, observer, datetime)) %>%
    rename(lat = latitude, lon = longitude) %>%
    dplyr::filter(!is.na(lat),
                  lat > 0,
                  !is.na(date)) %>%
    dplyr::filter(julian > quantile(julian, 0.005),
                  julian < quantile(julian, 0.995)) %>%
    dplyr::filter(id %in% visit$id) %>%
    dplyr::select(id, organization, project, location, lat, lon, observer, datetime, distanceMethod, durationMethod, speciesCode, distanceBand, durationInterval, abundance, vocalization)

#3. Add species list----
species <- read.csv("data/singing-species.csv") %>%
    rename(speciesCode = Species_ID)

#4. Save----
save(visit, bird, species,  file="data/cleaned_data_2022-07-24.Rdata")
