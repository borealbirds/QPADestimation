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

#A. CLEAN VISIT DATA####################

#1. Load in new dataset----
load(file.path(root, "wildtrax_raw_2023-01-20.Rdata"))

#2. Subset to columns of interest----

colnms <- c("project", "sensor", "location", "buffer", "lat", "lon", "year", "date", "observer", "ARUMethod", "durationMethod", "distanceMethod", "species", "abundance", "individual", "durationInterval", "distanceBand", "tagStart", "isHeard", "isSeen")

use <- raw.wt %>% 
  dplyr::select(-observer) %>% 
  rename(species = speciesCode, lat = latitude, lon = longitude, individual = individual_appearance_order, ARUMethod = method, tagStart = tag_start_s, observer = observer_id) %>% 
  full_join(projects %>% 
              rename(sensor = sensorId) %>% 
              dplyr::select(project_id, project, sensor)) %>% 
  mutate(date = ymd_hms(date),
         year = year(date)) %>% 
  dplyr::select(all_of(colnms)) %>% 
  mutate(id = paste(project, location, lat, lon, observer, date, ARUMethod)) %>% 
  dplyr::filter(!is.na(date))

#2. Wrangle duration method----
#Remove surveys that are not a factor of 60s
#Remove tags that are after the indicated duration
#create combination of sensor and tagmethod
#identify ARU surveys in the PC sensor and treat accordingly
method <- use %>%
  dplyr::filter(sensor=="ARU") %>% 
  separate(ARUMethod, into=c("duration", "tagMethod"), sep=" ", remove=FALSE) %>%
  mutate(minutes = as.numeric(str_sub(duration, -100, -2))/60,
         tagMethod = paste0("ARU-", tagMethod)) %>% 
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
  dplyr::select(c(colnames(use), tagMethod)) %>% 
  rbind(use %>% 
          dplyr::filter(sensor=="PC" & distanceMethod != "0m-INF-ARU") %>% 
          mutate(tagMethod="PC")) %>% 
  rbind(use %>% 
          dplyr::filter(sensor=="PC" & distanceMethod=="0m-INF-ARU") %>% 
          mutate(tagMethod="ARU-1SPT",
                 sensor="ARU"))

#3. Subset to visits & filter----
#Remove surveys with no location
#Standardize sig figs for location to remove duplicates
#Remove outliers for day of year (use 99% quantile)
#Take out BBS because isn't useful for removal or distance sampling
#Remove none tag methods
dat <- method %>%
  dplyr::filter(!is.na(date),
                project!="BAM-BBS",
                tagMethod!="None") %>% 
    dplyr::select(id, project, sensor, location, buffer, lat, lon, year, date, observer, distanceMethod, durationMethod, tagMethod) %>%
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

#check distribution against for time zone issues
ggplot(temporal) +
  geom_histogram(aes(x=start_time, fill=sensor))
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
unzip(zipfile=temp, exdir=file.path(root, "gis"))
unlink(temp)

#6b.Read in & wrangle shapefile
bcrshp <- read_sf(file.path(root, "gis", "BCR_Terrestrial/BCR_Terrestrial_master.shp")) %>%
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
    mutate(tsg = (yday(date)-seedgrow)/365,
           jday = julian/365,
           tssr = hssr/24,
           methodlength = nchar(durationMethod)) %>%
    group_by(id) %>% 
    dplyr::filter(methodlength == max(methodlength)) %>% 
    sample_n(1) %>% 
  ungroup() %>% 
  dplyr::select(id, project, sensor, location, buffer, lat, lon, year, date, observer, distanceMethod, durationMethod, tagMethod, jday, hssr, tssr, seedgrow, tsg, bcr, province, country, tree, lcc2, lcc4)

#B. CLEAN BIRD DATA################

#1. Filter, tidy, create foreign key for visit table----
#Ensure there's visit data
#Remove UNSPs
#Remove records without abundance
#Replace GRAJ with CAJA
#Remove records with auditory detections
dat <- use %>% 
  dplyr::filter(str_sub(species, 1, 2)!="UN",
                id %in% visit$id,
                !is.na(abundance),
                !abundance %in% c("CI 1", "CI 2", "CI 3", "N/A", "0", ""),
                !isHeard %in% c("f", "no", "No")) %>% 
  mutate(species = ifelse(species=="GRAJ", "CAJA", species))

#2. Filter to just first detection per individual and bin in minutes----
#bin in 30 s bins for 60 s surveys
#Fill in distance & duration method from visit object
first <- dat %>%
  dplyr::filter(sensor=="ARU") %>% 
  dplyr::select(-distanceMethod, -durationMethod) %>% 
  left_join(visit %>% 
              dplyr::select(id, distanceMethod, durationMethod, tagMethod)) %>%
  group_by(id, project, sensor, location, buffer, lat, lon, year, date,  observer, species, abundance, individual, isSeen, isHeard) %>% 
  mutate(firstTag = min(tagStart)) %>%
  ungroup() %>%
  dplyr::filter(tagStart == firstTag) %>% 
  mutate(end=ifelse(durationMethod=="0-0.5-1min", ceiling(tagStart/30), ceiling(tagStart/60)),
         end = ifelse(end==0, 1, end),
         durationInterval = case_when(durationMethod=="0-0.5-1min" & end==1 ~ "0-0.5min",
                                      durationMethod=="0-0.5-1min" & end==2 ~ "0.5-1min",
                                      !is.na(durationMethod) ~ paste0(end-1, "-", end, "min")),
         distanceBand = "UNKNOWN") %>% 
  dplyr::select(c(colnames(dat), tagMethod)) %>% 
  rbind(dat %>% 
          dplyr::filter(sensor=="PC" & distanceMethod!="0m-INF-ARU") %>% 
          mutate(tagMethod="PC")) %>%
  rbind(dat %>% 
          dplyr::filter(sensor=="PC" & distanceMethod=="0m-INF-ARU") %>% 
          mutate(tagMethod="ARU-1SPT")) %>% 
  dplyr::select(id, project, sensor, location, buffer, lat, lon, year, date, observer, distanceMethod, durationMethod, tagMethod, distanceBand, durationInterval, species, abundance, isSeen, isHeard)

#3. Replace TMTTs with predicted abundance----
tmtt <- read.csv("C:/Users/Elly Knight/Documents/ABMI/Projects/Wildtrax/TMTT/data/tmtt_predictions_mean.csv") %>% 
  rename(species = species_code, observer = observer_id)

bird <- first %>% 
  dplyr::filter(abundance=="TMTT") %>% 
  mutate(species = ifelse(species %in% tmtt$species, species, "species"),
         observer = as.integer(ifelse(observer %in% tmtt$observer, observer, 0))) %>% 
  data.frame() %>% 
  left_join(tmtt) %>% 
  mutate(abundance = round(pred)) %>% 
  dplyr::select(colnames(first)) %>% 
  rbind(first %>% 
          dplyr::filter(abundance!="TMTT") %>% 
          mutate(abundance = as.numeric(abundance))) 

#B. PACKAGE AND SAVE###########################

#1. Rename to match QPAD V3-----
visit <- visit %>% 
  rename(TSSR = tssr, JDAY = jday, DSLS = tsg, LCC2 = lcc2, LCC4 = lcc4, TREE = tree,
         TM = tagMethod)
bird <- bird %>% 
  rename(TM = tagMethod)

#2. Add species list----
species <- read.csv(file.path(root, "lookups", "singing-species.csv")) %>%
  rename(species = Species_ID)

#3. Save----
save(visit, bird, species,  file="G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/qpadv4_clean.Rdata")

