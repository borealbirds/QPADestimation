# ---
# title: "QPAD estimation - clean visit data"
# author: "Elly Knight"
# created: "July 24, 2022"
# updated: "November 6, 2022"
# ---

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(lubridate) #date manipulation
library(suncalc) #sunrise time retrieval
library(sf) #spatial manipulation
library(terra) #raster handling
library(downloader) #download region file

#2. Set root path for data on google drive----
root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/"

#1. Load in new dataset----
load(file.path(root, "/qpadv4_raw.Rdata"))

#2. Wrangle duration method----
#Remove surveys with "none" method
#Remove surveys that are not a factor of 60s
#Remove tags that are after the indicated duration
method <- use %>%
  dplyr::filter(sensor=="ARU") %>% 
  separate(ARUMethod, into=c("duration", "tagMethod"), sep=" ", remove=FALSE) %>%
  mutate(minutes = as.numeric(str_sub(duration, -100, -2))/60) %>% 
  dplyr::filter(minutes %in% c(1:10),
                tagStart <= minutes*60) %>%
  mutate(durationMethod = case_when(minutes==1 ~ "0-0.5-1min",
                                    minutes==2 ~ "0-1-2min",
                                    minutes==3 ~ "0-1-2-3min",
                                    minutes==4 ~ "0-1-2-3-4min",
                                    minutes==5 ~ "0-1-2-3-4-5min",
                                    minutes==6 ~ "0-1-2-3-4-5-6min",
                                    minutes==7 ~ "0-1-2-3-4-5-6-7min",
                                    minutes==8 ~ "0-1-2-3-4-5-6-7-8min",
                                    minutes==9 ~ "0-1-2-3-4-5-6-7-8-9min",
                                    minutes==10 ~ "0-1-2-3-4-5-6-7-8-9-10min")) %>% 
  dplyr::select(colnames(use)) %>% 
  rbind(use %>% 
          dplyr::filter(sensor=="PC"))

#3. Subset to visits & filter----
#Remove surveys with no location
#Standardize sig figs for location to remove duplicates
#Remove outliers for day of year (use 99% quantile)
#Take out BBS because isn't useful for removal or distance sampling
dat <- method %>%
  dplyr::filter(!is.na(date),
                project!="BAM-BBS") %>% 
    dplyr::select(id, source, project, sensor, singlesp, location, buffer, lat, lon, year, date, observer, distanceMethod, durationMethod) %>%
    mutate(lat = round(lat, 5),
           lon = round(lon, 5),
           buffer = ifelse(is.na(buffer), 0, buffer),
           julian = yday(date)) %>% 
    unique() %>% 
  dplyr::filter(!is.na(lat),
                lat > 0, 
                lon < 0,
                !is.na(date),
                year(date) > 1900,
                julian > quantile(julian, 0.005),
                julian < quantile(julian, 0.995))

#4. Wrangle temporal variables----
temporal <- dat %>%
    mutate(datetime = ymd_hms(date),
           year = year(datetime),
           julian = yday(datetime),
           start_time = hour(datetime) + minute(datetime)/60) %>%
    dplyr::filter(year > 1900)

#check BAM dataset distribution against WT for time zone issues
ggplot(temporal) +
  geom_histogram(aes(x=start_time, fill=source))
#looks OK

#5. Get sunrise time----
sun <- temporal %>%
    mutate(date = ymd(str_sub(datetime, 1, 10))) %>%
    rename(lat = lat, lon = lon)

sun$sunrise <- getSunlightTimes(data=sun, keep="sunrise")$sunrise
sun$hssr <- as.numeric(difftime(sun$datetime, sun$sunrise), units="hours")

#6. Get region----

#6a. Download data
temp <- tempfile()
download("https://birdscanada.org/download/gislab/bcr_terrestrial_shape.zip", temp)
unzip(zipfile=temp, exdir=root)
unlink(temp)

#6b.Read in & wrangle shapefile
bcrshp <- read_sf("gis/BCR_Terrestrial/BCR_Terrestrial_master.shp") %>%
    dplyr::select(BCR, PROVINCE_S, COUNTRY) %>%
    rename(bcr = BCR, province = PROVINCE_S, country=COUNTRY) %>%
    dplyr::filter(country %in% c("USA", "CANADA")) %>%
    mutate(countryid = as.numeric(as.factor(country)),
           provinceid = as.numeric(as.factor(province))) %>%
    st_transform(crs=3857) %>%
    st_make_valid() %>%
    vect()

#6c. Create rasters (much faster than from polygon)
r <- rast(ext(bcrshp), resolution=1000)

bcr <- rasterize(x=bcrshp, y=r, field="bcr")
province <- rasterize(x=bcrshp, y=r, field="provinceid")
country <- rasterize(x=bcrshp, y=r, field="countryid")

regionstack <- rast(list(bcr, province, country))

#6d. Extract values
regionids <- sun %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    st_transform(crs=3857) %>%
    vect() %>%
    terra::extract(x=regionstack)

#6e. Create lookup table to join back ids
bcrtbl <- st_as_sf(bcrshp) %>%
    as.data.frame() %>%
    dplyr::select(-geometry) %>%
    unique()

#6f. Put together
region <- sun %>%
    cbind(regionids) %>%
    dplyr::select(-ID) %>%
    left_join(bcrtbl)

#7. Get covariates----

#7a. Download data
download("https://github.com/borealbirds/qpad-offsets/tree/main/data/lcc.tif", destfile=file.path(root, "lcc.tif"))
download("https://github.com/borealbirds/qpad-offsets/tree/main/data/seedgrow.tif", destfile=file.path(root, "seedgrow.tif"))
download("https://github.com/borealbirds/qpad-offsets/tree/main/data/tree.tif", destfile=file.path(root, "tree.tif"))

#7b. Read in rasters----
lcc <- rast("gis/lcc.tif")
sg <- rast("gis/seedgrow.tif")
tree <- rast("gis/tree.tif")

#7c. Get values----
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

#7c.Create lookup table for lcc
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
#7d. Put together----
covariates <- region %>%
    cbind(lccval, sgval, treeval) %>%
    left_join(lcctbl) %>%
    mutate(tree = tree/100,
           tree = ifelse(tree > 1, 0, tree))

#8. Tidy, create primary key, standardize----
#select processing method with higher resolution for recordings that are processed twice
visit <- covariates %>%
    mutate(tsg = (as.numeric(date)-seedgrow)/365,
           jday = julian/365,
           tssr = hssr/24,
           tsg = tsg/365,
           methodlength = nchar(durationMethod)) %>%
    group_by(id) %>% 
    dplyr::filter(methodlength == max(methodlength)) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    dplyr::select(id, source, project, sensor, singlesp, location, buffer, lat, lon, year, date, observer, distanceMethod, durationMethod, julian, jday, hssr, tssr, seedgrow, tsg, bcr, province, country, tree, lcc2, lcc4)

#9. Save----
save(visit, file="G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/qpadv4_visit.Rdata")
