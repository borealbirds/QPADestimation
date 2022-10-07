# ---
# title: "QPAD estimation - test inclusion of ARU data in removal models"
# author: "Elly Knight"
# created: "Oct 6, 2022"
# ---

library(tidyverse) #basic data wrangling
library(wildRtrax) #download data from wildtrax api
library(data.table) #collapse list to dataframe
library(lubridate) #date manipulation
library(suncalc) #sunrise time retrieval
library(sf) #spatial manipulation
library(terra) #raster handling
library(downloader) #download region file
library(detect) #removal models

#A. GET DATA####

#1. Get list of projects from WildTrax----
wt_auth()

#sensor = PC gives all ARU and point count projects
project.list <- wt_get_download_summary(sensor_id = 'ARU')

#Convert to a plain dataframe
projects.aru <- data.frame(project = as.character(project.list$project),
                       project_id = as.numeric(project.list$project_id),
                       sensorId = as.character(project.list$sensorId),
                       tasks = as.numeric(project.list$tasks),
                       status = as.character(project.list$status)) %>%
    dplyr::filter(sensorId=="ARU")

dat.list <- list()
for(i in 1:nrow(projects.aru)){

    try(dat.list[[i]] <- wt_download_report(project_id = projects.aru$project_id[i], sensor_id = "ARU", cols_def = F))

    print(paste0("Finished dataset ", projects.aru$project[i], " : ", i, " of ", nrow(projects.aru), " projects"))

}

#4. Collapse list----
raw.aru <- rbindlist(dat.list, fill=TRUE)

#5. Save date stamped data & project list----
save(raw.aru, projects.aru, file=paste0("data/wildtrax_data_aru_", Sys.Date(), ".Rdata"))
load("data/wildtrax_data_aru_2022-10-07.Rdata")

#B. CLEAN VISIT DATA####

#1. Subset to visits & filter----
#Remove surveys with no location
#Standardize sig figs for location to remove duplicates
#Remove human surveys (have date instead of recording_date)
#Remove surveys with "none" method
dat <- raw.aru %>%
    dplyr::select(organization, project_name, location, latitude, longitude, observer, recording_date, method) %>%
    rename(project = project_name, durationMethod = method) %>%
    unique() %>%
    mutate(latitude = round(latitude, 5),
           longitude = round(longitude, 5)) %>%
    dplyr::filter(!is.na(latitude),
                  latitude > 0,
                  !is.na(recording_date)) %>%
    dplyr::filter(durationMethod!="None") %>%
    unique()

#2. Wrangle temporal variables----
temporal <- dat %>%
    mutate(datetime = ymd_hms(recording_date),
           year = year(datetime),
           julian = yday(datetime),
           start_time = hour(datetime) + minute(datetime)/60) %>%
    dplyr::filter(year > 1900)

#3. Get sunrise time----
sun <- temporal %>%
    mutate(date = ymd(str_sub(datetime, 1, 10))) %>%
    rename(lat = latitude, lon = longitude)

sun$sunrise <- getSunlightTimes(data=sun, keep="sunrise")$sunrise
sun$hssr <- as.numeric(difftime(sun$datetime, sun$sunrise), units="hours")

#4. Get region----

#4a. Download data
# temp <- tempfile()
# download("https://birdscanada.org/download/gislab/bcr_terrestrial_shape.zip", temp)
# unzip(zipfile=temp, exdir="gis")
# unlink(temp)

#4b.Read in & wrangle shapefile
bcrshp <- read_sf("gis/BCR_Terrestrial/BCR_Terrestrial_master.shp") %>%
    dplyr::select(BCR, PROVINCE_S, COUNTRY) %>%
    rename(bcr = BCR, province = PROVINCE_S, country=COUNTRY) %>%
    dplyr::filter(country %in% c("USA", "CANADA")) %>%
    mutate(countryid = as.numeric(as.factor(country)),
           provinceid = as.numeric(as.factor(province))) %>%
    st_transform(crs=3857) %>%
    st_make_valid() %>%
    vect()

#4c. Create rasters (much faster than from polygon)
r <- rast(ext(bcrshp), resolution=1000)

bcr <- rasterize(x=bcrshp, y=r, field="bcr")
province <- rasterize(x=bcrshp, y=r, field="provinceid")
country <- rasterize(x=bcrshp, y=r, field="countryid")

regionstack <- rast(list(bcr, province, country))

#4d. Extract values
regionids <- sun %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    st_transform(crs=3857) %>%
    vect() %>%
    terra::extract(x=regionstack)

#4e. Create lookup table to join back ids
bcrtbl <- st_as_sf(bcrshp) %>%
    as.data.frame() %>%
    dplyr::select(-geometry) %>%
    unique()

#4f. Put together
region <- sun %>%
    cbind(regionids) %>%
    dplyr::select(-ID) %>%
    left_join(bcrtbl)

#5. Get covariates----

#5a. Read in rasters----
lcc <- rast("gis/lcc.tif")
sg <- rast("gis/seedgrow.tif")
tree <- rast("gis/tree.tif")

#5b. Get values----
covsf <- region %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    st_transform(crs=crs(lcc)) %>%
    vect()

#Can't stack because extents don't match :(
lccval <- terra::extract(covsf, x=lcc) %>%
    dplyr::select(-ID)
sgval <- terra::extract(covsf, x=sg) %>%
    dplyr::select(-ID)
treeval <- terra::extract(covsf, x=tree) %>%
    dplyr::select(-ID)

#5c.Create lookup table for lcc
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
#5d. Put together----
covariates <- region %>%
    cbind(lccval, sgval, treeval) %>%
    left_join(lcctbl)

#6. Tidy, create primary key, standardize & create polynomial variables----
visit.aru <- covariates %>%
    mutate(tsg = (as.numeric(date)-seedgrow)/365,
           jday = julian/365,
           tssr = hssr/24,
           tsg = tsg/365,
           jday2 = jday^2,
           tssr2 = tssr^2,
           tsg2 = tsg^2) %>%
    dplyr::select(organization, project,  observer, durationMethod, location, lon, lat, datetime, year, julian, jday, jday2, hssr, tssr, tssr2, seedgrow, tsg, tsg2, bcr, province, country, tree, lcc2, lcc4) %>%
    mutate(id = paste(project, location, observer, datetime),
           distanceMethod = NA)

#C. CLEAN BIRD DATA####

#1. Filter, tidy, create foreign key for visit table----
#Remove surveys with no location
#Remove outliers for day of year (use 99% quantile)
#Ensure there's visit data
dat.bird <- raw.aru %>%
    rename(lat = latitude, lon = longitude, project = project_name, speciesCode = species_code, durationMethod = method) %>%
    mutate(datetime = ymd_hms(recording_date),
           julian = yday(datetime),
           id = paste(project, location, observer, datetime)) %>%
    dplyr::filter(!is.na(lat),
                  lat > 0,
                  !is.na(recording_date)) %>%
    dplyr::filter(julian > quantile(julian, 0.005),
                  julian < quantile(julian, 0.995)) %>%
    dplyr::filter(id %in% visit.aru$id) %>%
    dplyr::select(id, organization, project, location, lat, lon, observer, datetime, durationMethod, speciesCode, tag_start_s, abundance, individual_appearance_order)

#2. Filter to just first detection per individual and bin in minutes----
first <- dat.bird %>%
    group_by(id, organization, project, location, lat, lon, observer, datetime, durationMethod, speciesCode, abundance, individual_appearance_order) %>%
    mutate(first_tag = min(tag_start_s)) %>%
    ungroup() %>%
    dplyr::filter(tag_start_s == first_tag) %>%
    dplyr::select(-first_tag) %>%
    dplyr::filter(!abundance %in% c("CI 1", "CI 2", "CI 3", "N/A")) %>%
    mutate(end=ceiling(tag_start_s/60),
           end = ifelse(end==0, 1, end),
           durationInterval = paste0(end-1, "-", end, "min"))

#3. Replace TMTTs----
#99 percentile for combination of species & observer
tmtttbl <- first %>%
    dplyr::filter(abundance!="TMTT") %>%
    group_by(id, observer, speciesCode) %>%
    summarize(abundance = sum(as.numeric(abundance))) %>%
    group_by(speciesCode, observer) %>%
    summarize(q = quantile(abundance, 0.99)) %>%
    ungroup() %>%
    mutate(val = ceiling(q))

bird.aru <- first %>%
    left_join(tmtttbl) %>%
    mutate(abundance = ifelse(abundance=="TMTT", val, as.numeric(abundance)),
           distanceMethod = NA,
           distanceBand = NA) %>%
    dplyr::filter(!is.na(abundance)) %>%
    dplyr::select(id, organization, project, location, lat, lon, observer, datetime, durationMethod, distanceMethod, speciesCode, durationInterval, distanceBand, abundance)

#4. Save----
save(visit.aru, bird.aru, species,  file="data/cleaned_data_aru_2022-10-07.Rdata")

#D. REMOVAL MODELLING####

#1. Create list of models----
#jday = day of year as a decimal between 0 and 1
#tssr = time since sunrise as a decimal between 0 and 1
#tsg = days since start of seedgrowth from seedgrow layers

mods <- list(
    ~ 1,
    ~ jday,
    ~ tssr,
    ~ jday + jday2,
    ~ tssr + tssr2,
    ~ jday + tssr,
    ~ jday + jday2 + tssr,
    ~ jday + tssr + tssr2,
    ~ jday + jday2 + tssr + tssr2,
    ~ tsg,
    ~ tsg + tsg2,
    ~ tsg + tssr,
    ~ tsg + tsg2 + tssr,
    ~ tsg + tssr + tssr2,
    ~ tsg + tsg2 + tssr + tssr2)
names(mods) <- 0:14
modnames <- c("(Intercept)",
              "(Intercept) + jday",
              "(Intercept) + tssr",
              "(Intercept) + jday + jday2",
              "(Intercept) + tssr + tssr2",
              "(Intercept) + jday + tssr",
              "(Intercept) + jday + jday2 + tssr",
              "(Intercept) + jday + tssr + tssr2",
              "(Intercept) + jday + jday2 + tssr + tssr2",
              "(Intercept) + tsg",
              "(Intercept) + tsg + tsg2",
              "(Intercept) + tsg + tssr",
              "(Intercept) + tsg + tsg2 + tssr",
              "(Intercept) + tsg + tssr + tssr2",
              "(Intercept) + tsg + tsg2 + tssr + tssr2")

#2. Load point count data and put together----
load("data/cleaned_data_2022-10-06.Rdata")

visit.all <- rbind(visit, visit.aru)
bird.all <- rbind(bird, bird.aru)

#3. Create design lookup table that describes duration method for each protocol----
#filter out duration methods that aren't appropriate for removal modelling (only have 1 time bin)
design <- visit %>%
    dplyr::select(durationMethod) %>%
    unique() %>%
    dplyr::filter(!durationMethod %in% c("UNKNOWN", "0-10min", "0-20min", "0-5min", "0-3min", "0-2min")) %>%
    mutate(dm = str_sub(durationMethod, -100, -4)) %>%
    separate(dm, into=c("t00", "t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), remove=TRUE, sep="-") %>%
    dplyr::select(-t00) %>%
    mutate_at(c("t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), ~as.numeric(.))

#4. Get list of species to process----
spp <- species %>%
    dplyr::filter(Singing_birds==TRUE) %>%
    left_join(bird %>%
                  dplyr::select(speciesCode) %>%
                  unique()) %>%
    arrange(speciesCode)

#5. Set up loop for species----
species.list <- list()
for(i in 1:nrow(spp)){

    #6. Filter abundance data for species---
    # filter out observations with unknown duration method or interval
    # filter to observations with covariates
    bird.i <- bird %>%
        dplyr::filter(speciesCode==spp$speciesCode[i],
                      durationMethod %in% design$durationMethod,
                      durationInterval!="UNKNOWN") %>%
        group_by(id, durationMethod, durationInterval) %>%
        summarize(abundance = sum(abundance)) %>%
        ungroup()

    #only model if there is data
    if(nrow(bird.i) > 0){

        #7. Replace abundance outliers (defined as 99% quantile across species) with 99% quantile of abundance----
        if(max(bird.i$abundance) > quantile(bird$abundance, 0.99)){
            bird.i <- bird.i %>%
                mutate(abundance = ifelse(abundance > quantile(abundance, 0.99), round(quantile(abundance, 0.99)), abundance))
        }

        #8. Filter visit covariates----
        #Remove records with nas in covariates
        x <- visit %>%
            dplyr::filter(id %in% unique(bird.i$id),
                          !is.na(tssr),
                          !is.na(tsg),
                          !is.na(jday)) %>%
            arrange(id) %>%
            dplyr::select(id, durationMethod, jday, jday2, tssr, tssr2, tsg, tsg2)

        #9. Create design matrix----
        d <- x %>%
            dplyr::select(durationMethod) %>%
            left_join(design, by="durationMethod") %>%
            dplyr::select(-durationMethod) %>%
            as.matrix()

        #10. Format abundance matrix----
        #add dummy variables for each of the columns to make sure the matrix is the right width
        y <- bird.i %>%
            dplyr::filter(id %in% x$id) %>%
            separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE, extra="drop", fill="right") %>%
            mutate(start = as.numeric(start),
                   end = as.numeric(str_sub(end, -100, -4))) %>%
            left_join(design %>%
                          pivot_longer(t01:t10, values_to="end", names_to="position"),
                      by=c("durationMethod", "end")) %>%
            dplyr::select(id, position, abundance) %>%
            arrange(position) %>%
            rbind(data.frame(id="dummy", position=colnames(d), abundance=NA)) %>%
            pivot_wider(id_cols=id, names_from=position, values_from=abundance, values_fill=0) %>%
            dplyr::filter(id!="dummy") %>%
            arrange(id) %>%
            dplyr::select(-id) %>%
            as.matrix()

        #11. Change zeros to NAs in the abundance matrix to match the design matrix----
        for (j in 1:nrow(y)){
            indices <- which(is.na(d[j,]))
            y[j, indices] <- NA
        }

        #12. Fit models----
        #Save a bunch of metadata like sample size and aic value
        mod.list <- list()
        for (j in 1:length(mods)) {
            f <- as.formula(paste0("y | d ", paste(as.character(mods[[j]]), collapse=" ")))
            mod <- try(cmulti(f, x, type="rem"))
            if (!inherits(mod, "try-error")) {
                rmvl <- data.frame(t(data.frame(mod["coefficients"]))) %>%
                    mutate(nobs=mod["nobs"]$nobs,
                           loglik = mod["loglik"]$loglik,
                           df = length(coef(mod)),
                           aic = AIC(mod),
                           aicc = aic + (2*df^2+2*df) / (nrow(y)-df-1),
                           model = modnames[j],
                           species = spp$speciesCode[i])
            } else {
                rmvl <- data.frame(nobs="try-error",
                                   species=spp$speciesCode[i])
            }
            mod.list[[j]] <- rmvl
        }

        #13. Save model results to species list---
        species.list[[i]] <- rbindlist(mod.list, fill=TRUE)

    }

    print(paste0("Finished modelling species ", spp$English_Name[i], ": ", i, " of ", nrow(spp), " species"))

}

species.out <- rbindlist(species.list, fill=TRUE)

#14. Identify species that failed----
species.fail <- species.out %>%
    dplyr::filter(nobs=="try-error") %>%
    dplyr::select(species) %>%
    unique()

species.use <- species.out %>%
    dplyr::filter(!species %in% species.fail$species)

#15. Save
save(species.out, species.use, file="results/availability_results_aru_2022-10-07.Rdata")
