# ---
# title: "QPAD estimation - clean visit data"
# author: "Elly Knight"
# created: "July 24, 2022"
# updated: "September 22, 2022"
# ---

library(tidyverse) #basic data wrangling
library(lubridate) #date manipulation
library(suncalc) #sunrise time retrieval
library(sf) #spatial manipulation
library(terra) #raster handling
library(downloader) #download region file

#Load previous dataset to look at it----
#load into separate environment to avoid overwriting things
# e <- new.env()
# load("data/new_offset_data_package_2017-03-01.Rdata", envir = e)
# #look at data structure
# str(e$dat)

#1. Load in new dataset----
load("data/wildtrax_data_2022-10-06.Rdata")

#2. Subset to visits & filter----
#Remove surveys with no location
#Standardize sig figs for location to remove duplicates
#Remove ARU surveys (have recording_date instead of date)
dat <- raw %>%
    dplyr::select(organization, project, location, latitude, longitude, observer, date, distanceMethod, durationMethod) %>%
    unique() %>%
    mutate(latitude = round(latitude, 5),
           longitude = round(longitude, 5)) %>%
    dplyr::filter(!is.na(latitude),
                  latitude > 0,
                  !is.na(date)) %>%
        unique()

#3. Wrangle temporal variables----
temporal <- dat %>%
    mutate(datetime = ymd_hms(date),
           year = year(datetime),
           julian = yday(datetime),
           start_time = hour(datetime) + minute(datetime)/60) %>%
    dplyr::filter(year > 1900)

#4. Get sunrise time----
sun <- temporal %>%
    mutate(date = ymd(str_sub(datetime, 1, 10))) %>%
    rename(lat = latitude, lon = longitude)

sun$sunrise <- getSunlightTimes(data=sun, keep="sunrise")$sunrise
sun$hssr <- as.numeric(difftime(sun$datetime, sun$sunrise), units="hours")

#5. Get region----

#5a. Download data
# temp <- tempfile()
# download("https://birdscanada.org/download/gislab/bcr_terrestrial_shape.zip", temp)
# unzip(zipfile=temp, exdir="gis")
# unlink(temp)

#5b.Read in & wrangle shapefile
bcrshp <- read_sf("gis/BCR_Terrestrial/BCR_Terrestrial_master.shp") %>%
    dplyr::select(BCR, PROVINCE_S, COUNTRY) %>%
    rename(bcr = BCR, province = PROVINCE_S, country=COUNTRY) %>%
    dplyr::filter(country %in% c("USA", "CANADA")) %>%
    mutate(countryid = as.numeric(as.factor(country)),
           provinceid = as.numeric(as.factor(province))) %>%
    st_transform(crs=3857) %>%
    st_make_valid() %>%
    vect()

#5c. Create rasters (much faster than from polygon)
r <- rast(ext(bcrshp), resolution=1000)

bcr <- rasterize(x=bcrshp, y=r, field="bcr")
province <- rasterize(x=bcrshp, y=r, field="provinceid")
country <- rasterize(x=bcrshp, y=r, field="countryid")

regionstack <- rast(list(bcr, province, country))

#5d. Extract values
regionids <- sun %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    st_transform(crs=3857) %>%
    vect() %>%
    terra::extract(x=regionstack)

#5e. Create lookup table to join back ids
bcrtbl <- st_as_sf(bcrshp) %>%
    as.data.frame() %>%
    dplyr::select(-geometry) %>%
    unique()

#5f. Put together
region <- sun %>%
    cbind(regionids) %>%
    dplyr::select(-ID) %>%
    left_join(bcrtbl)

#6. Get covariates----

#6a. Read in rasters----
lcc <- rast("gis/lcc.tif")
sg <- rast("gis/seedgrow.tif")
tree <- rast("gis/tree.tif")

#6b. Get values----
covsf <- region %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    st_transform(crs=crs(lcc)) %>%
    vect()

lccval <- terra::extract(covsf, x=lcc) %>%
    dplyr::select(-ID)
sgval <- terra::extract(covsf, x=sg) %>%
    dplyr::select(-ID)
treeval <- terra::extract(covsf, x=tree) %>%
    dplyr::select(-ID)

#6c.Create lookup table for lcc
# 0: No data (NA/NA)
# 1: Temperate or sub-polar needleleaf forest (Conif/Forest)
# 2: Sub-polar taiga needleleaf forest (Conif/Forest)
# 5: Temperate or sub-polar broadleaf deciduous (DecidMixed/Forest)
# 6:  Mixed Forest (DecidMixed/Forest)
# 8: Temperate or sub-polar shrubland (Open/OpenWet)
# 10: Temperate or sub-polar grassland (Open/OpenWet)
# 11: Sub-polar or polar shrubland-lichen-moss (Open/OpenWet)
# 12: Sub-polar or polar grassland-lichen-moss (Open/OpenWet)
# 13: Sub-polar or polar barren-lichen-moss (Open/OpenWet)
# 14: Wetland (Wet/OpenWet)
# 15: Cropland (Open/OpenWet)
# 16: Barren Lands (Open/OpenWet)
# 17: Urban and Built-up (Open/OpenWet)
# 18: Water (NA/NA)
# 19: Snow and Ice (NA/NA)
lcctbl <- data.frame(lcc=c(0:19),
                     lcc4=c("", "Conif", "Conif", "", "", "DecidMixed", "DecidMixed", "", "Open", "", "Open", "Open", "Open", "Open", "Wet", "Open", "Open", "Open", "", "")) %>%
    mutate(lcc2 = case_when(lcc4 %in% c("Conif", "DecidMixed") ~ "Forest",
                            lcc4 %in% c("Open", "Wet") ~ "OpenWet",
                            !is.na(lcc4) ~ lcc4))
#6d. Put together----
covariates <- region %>%
    cbind(lccval, sgval, treeval) %>%
    left_join(lcctbl)

#7. Tidy, create primary key, standardize & create polynomial variables----
visit <- covariates %>%
    mutate(tsg = (as.numeric(date)-seedgrow)/365,
           jday = julian/365,
           tssr = hssr/24,
           tsg = tsg/365,
           jday2 = jday^2,
           tssr2 = tssr^2,
           tsg2 = tsg^2) %>%
    dplyr::select(organization, project,  observer, distanceMethod, durationMethod, location, lon, lat, datetime, year, julian, jday, jday2, hssr, tssr, tssr2, seedgrow, tsg, tsg2, bcr, province, country, tree, lcc2, lcc4) %>%
    mutate(id = paste(location, observer, datetime))

#8. Save----
save(visit, file="data/visit_data_2022-10-06.Rdata")
