# ---
# title: "QPAD estimation - clean bird data"
# author: "Elly Knight"
# created: "September 22, 2022"
# ---

#NOTE: NEED TO DEAL WITH DISTANCE AND DURATION BINS####

library(tidyverse) #basic data wrangling

#Load previous dataset to look at it----
#load into separate environment to avoid overwriting things
e <- new.env()
load("data/new_offset_data_package_2017-03-01.Rdata", envir = e)
#look at data structure
str(e$pc)

#1. Load in raw  & visit dataset----
load("data/wildtrax_data_2022-07-24.Rdata")
load("data/cleaned_data_2022-07-24.Rdata")

#2. Filter & wrangle----
#Remove surveys with no location
#Remove surveys with no temporal or spatial structure within a visit
bird <- raw %>%
    dplyr::filter(!is.na(latitude),
                  latitude > 0,
                  !(distanceMethod %in% c("0m-INF", NA, "0m-INF-ARU", "UNKNOWN") &
                        durationMethod %in% c("0-3min", "0-5min", "0-10min", "0-2min", "0-20min", "UNKNOWN", NA))) %>%
    dplyr::select(organization, project, location, latitude, longitude, observer, date, recording_date, distanceMethod, durationMethod, speciesCode, distanceBand, durationInterval, abundance)

#3. Save----
save(visit, bird,  file="data/cleaned_data_2022-07-24.Rdata")
