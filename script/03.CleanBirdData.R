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

#2. Filter----
#Remove surveys with no location
#Remove surveys with no temporal or spatial structure within a visit
#Remove outliers for day of year
use <- raw %>%
    mutate(datetime = ymd_hms(date),
           julian = yday(datetime),
           abundance = as.numeric(abundance)) %>%
    dplyr::filter(!is.na(latitude),
                  latitude > 0,
                  !(distanceMethod %in% c("0m-INF", NA, "0m-INF-ARU", "UNKNOWN") &
                        durationMethod %in% c("0-3min", "0-5min", "0-10min", "0-2min", "0-20min", "UNKNOWN", NA))) %>%
    dplyr::filter(julian > quantile(julian, 0.01),
                  julian < quantile(julian, 0.99))

#3. Tidy, & create foreign key for visit table----
bird <- use %>%
    rename(lat = latitude, lon = longitude) %>%
    dplyr::select(organization, project, location, lat, lon, observer, datetime, distanceMethod, durationMethod, speciesCode, distanceBand, durationInterval, abundance, vocalization) %>%
    mutate(id = paste(location, observer, datetime))

#4. Add species list----
species <- read.csv("data/singing-species.csv") %>%
    rename(speciesCode = Species_ID)

#5. Save----
visit <- visit %>%
    dplyr::filter(!is.na(tsg),
                  !is.na(tssr))  %>%
    mutate(id = paste(location, observer, datetime)) %>%
    mutate(lat = round(lat, 5),
           lon = round(lon, 5)) %>%
    mutate(tsg = tsg/365,
           tsg2 = tsg^2) %>%
    unique()

save(visit, bird, species,  file="data/cleaned_data_2022-07-24.Rdata")
