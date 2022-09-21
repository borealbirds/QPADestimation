# ---
# title: "QPAD estimation - clean visit data"
# author: "Elly Knight"
# created: "July 24, 2022"
# ---

library(tidyverse) #basic data wrangling
library(lubridate) #date manipulation
library(suncalc) #sunrise time retrieval
#library(rnaturalearth) #get state/province data
library(sf) #spatial manipulation
library(downloader) #download spatial data

#NOTE: INVESTIGATE MISSING LAT LONS####
#NOTE: NEED TO GET ARU MODEL####

#1. Load previous dataset----
#load into separate environment to avoid overwriting things
e <- new.env()
load("data/new_offset_data_package_2017-03-01.Rdata", envir = e)

#2. Load in new dataset----
load("data/wildtrax_data_2022-07-24.Rdata")

#3. Subset to visits----
visit <- raw %>%
    dplyr::filter(!is.na(latitude)) %>%
    dplyr::select(organization, project, location, latitude, longitude, observer, date, recording_date, distanceMethod, durationMethod) %>%
    unique() %>%
    sample_n(1000)
#Remove surveys with no location
#Take subsample for testing

#4. Determine survey type---
type <- visit %>%
    mutate(type = case_when(!is.na(date) ~ "human",
                            !is.na(recording_date) ~ "aru"))

#5. Wrangle temporal variables----
temporal <- type %>%
    mutate(datetime = case_when(type=="human" ~ ymd_hms(date),
                                type=="aru" ~ ymd_hms(recording_date)),
           year = year(datetime),
           julian = yday(datetime),
           start_time = hour(datetime) + minute(datetime)/60)

#6. Get sunrise time----

#NOTE: HOW TO DEAL WITH SURVEYS TOO FAR NORTH FOR SUNRISE????####
#Peter just had NAs in his data, so presumably removed them
summary(e$dat %>% dplyr::filter(Y>66) %>% dplyr::select(TSSR, srise, Y))

sun <- temporal %>%
    mutate(date = ymd(str_sub(datetime, 1, 10))) %>%
    rename(lat = latitude, lon = longitude)

sun$sunrise <- getSunlightTimes(data=sun, keep="sunrise")$sunrise
sun$tssr <- as.numeric(difftime(sun$datetime, sun$sunrise), units="hours")

#7. Lookup method codes----

#NOTE: CHANGE TO DOWNLOADING DIRECTLY FROM GIT REPO####
#NOTE: CREATE ALPHA CODE OR USE HISTORIC?####

dist_codes <- read.csv("data/distance_band_codes.csv")
dur_codes <- read.csv("data/durint.csv")

#8. Get region----

#Download data
# temp <- tempfile()
# download("https://birdscanada.org/download/gislab/bcr_terrestrial_shape.zip", temp)
# unzip(zipfile=temp, exdir="gis")
# unlink(temp)

#Read in & wrangle shapefile
bcrshp <- read_sf("gis/BCR_Terrestrial/BCR_Terrestrial_master.shp") %>%
    dplyr::select(BCR, BCRNAME, PROVINCE_S, COUNTRY) %>%
    rename(bcr = BCR, bcrname = BCRNAME, province = PROVINCE_S, country = COUNTRY) %>%
    st_transform(crs=4326) %>%
    st_make_valid()

#Intersect with visits
bcr <- sun %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    st_intersection(bcrshp)

#9. Get covariates----
#road
#tree, tree3
#hab_NALC1
#hab_NALC2
#SPRNG
#TSLS
#DD5

#6. Get landcover covariates----

#7. Get climate covariates----



#10. Tidy & write out----


