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
load("data/wildtrax_data_aru_2022-10-15.Rdata")

#B. CLEAN VISIT DATA####

#1. Wrangle duration method----
#Remove surveys with "none" method
#Remove surveys that are not a factor of 60s
#Remove tags that are after the indicated duration
method <- raw.aru %>%
    separate(method, into=c("duration", "tagMethod"), sep=" ", remove=FALSE) %>%
    mutate(minutes = as.numeric(str_sub(duration, -100, -2))/60) %>%
    dplyr::filter(tagMethod!="None",
                  minutes %in% c(1:10),
                  tag_start_s <= minutes*60) %>%
    mutate(durationMethod = case_when(minutes==1 ~ "0-0.5-1min",
                                      minutes==2 ~ "0-1-2min",
                                      minutes==3 ~ "0-1-2-3min",
                                      minutes==4 ~ "0-1-2-3-4min",
                                      minutes==5 ~ "0-1-2-3-4-5min",
                                      minutes==6 ~ "0-1-2-3-4-5-6min",
                                      minutes==7 ~ "0-1-2-3-4-5-6-7min",
                                      minutes==8 ~ "0-1-2-3-4-5-6-7-8min",
                                      minutes==9 ~ "0-1-2-3-4-5-6-7-8-9min",
                                      minutes==10 ~ "0-1-2-3-4-5-6-7-8-9-10min"))

#1. Subset to visits & filter----
#Remove surveys with no location
#Standardize sig figs for location to remove duplicates
#Remove human surveys (have date instead of recording_date)
dat <- method %>%
    dplyr::select(organization, project_name, location, latitude, longitude, observer, recording_date, durationMethod) %>%
    rename(project = project_name) %>%
    unique() %>%
    mutate(latitude = round(latitude, 5),
           longitude = round(longitude, 5)) %>%
    dplyr::filter(!is.na(latitude),
                  latitude > 0,
                  !is.na(recording_date)) %>%
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
    mutate(id = paste(project, location, observer, durationMethod, datetime),
           distanceMethod = NA)

#C. CLEAN BIRD DATA####

#1. Filter, tidy, create foreign key for visit table----
#Remove surveys with no location
#Remove outliers for day of year (use 99% quantile)
#Ensure there's visit data
#Remove UN sps
dat.bird <- method %>%
    rename(lat = latitude, lon = longitude, project = project_name, speciesCode = species_code) %>%
    mutate(datetime = ymd_hms(recording_date),
           julian = yday(datetime),
           id = paste(project, location, observer, durationMethod, datetime)) %>%
    dplyr::filter(!is.na(lat),
                  lat > 0,
                  !is.na(recording_date)) %>%
    dplyr::filter(julian > quantile(julian, 0.005),
                  julian < quantile(julian, 0.995)) %>%
    dplyr::filter(id %in% visit.aru$id) %>%
    dplyr::filter(str_sub(speciesCode, 1, 2)!="UN") %>%
    dplyr::select(id, organization, project, location, lat, lon, observer, datetime, durationMethod, speciesCode, tag_start_s, abundance, individual_appearance_order, species_comments)

#2. Filter to just first detection per individual and bin in minutes----
#bin in 30 s bins for 60 s surveys
first <- dat.bird %>%
    group_by(id, organization, project, location, lat, lon, observer, datetime, durationMethod, speciesCode, abundance, individual_appearance_order) %>%
    mutate(first_tag = min(tag_start_s)) %>%
    ungroup() %>%
    dplyr::filter(tag_start_s == first_tag) %>%
    dplyr::select(-first_tag) %>%
    dplyr::filter(!abundance %in% c("CI 1", "CI 2", "CI 3", "N/A")) %>%
    mutate(end=ifelse(durationMethod=="0-0.5-1min", ceiling(tag_start_s/30), ceiling(tag_start_s/60)),
           end = ifelse(end==0, 1, end),
           durationInterval = case_when(durationMethod=="0-0.5-1min" & end==1 ~ "0-0.5min",
                                        durationMethod=="0-0.5-1min" & end==2 ~ "0.5-1min",
                                        !is.na(durationMethod) ~ paste0(end-1, "-", end, "min")))

#3. Replace TMTTs----
#99 percentile for combination of species & observer
tmtttbl <- first %>%
    dplyr::filter(abundance!="TMTT") %>%
    group_by(id, observer, speciesCode) %>%
    summarize(abundance = sum(as.numeric(abundance))) %>%
    group_by(speciesCode, observer) %>%
    summarize(q = quantile(abundance, 0.99),
              nobs = n()) %>%
    ungroup() %>%
    mutate(val = ceiling(q))

#Replace observer species combos with less than 20 with generic "observer"
tmttn <- first %>%
    dplyr::filter(abundance=="TMTT") %>%
    group_by(observer, speciesCode) %>%
    summarize(ntmtt=n()) %>%
    ungroup() %>%
    left_join(tmtttbl) %>%
    dplyr::filter(!is.na(val)) %>%
    mutate(ptmtt = ntmtt/nobs) %>%
    mutate(obs = ifelse(ntmtt < 20, "observer", observer))

#Model 99% quantile with random effects for species and observer
lm.tmtt <- lme4::lmer(val ~ nobs + (1|speciesCode) + (1|obs), data=tmttn)

#Check variance explained by REs
summary(lm.tmtt)

#Predict to get values to replace TMTTs with
tmttpred <- data.frame(pred = predict(lm.tmtt)) %>%
    cbind(data.frame(tmttn)) %>%
    mutate(predabun = round(pred))

#Plot quantiles vs predictions
ggplot(tmttpred) +
    geom_jitter(aes(x=val, y=predabun, colour=log(nobs))) +
    xlab("Raw quantile") +
    ylab("Predicted value") +
    geom_abline(aes(intercept=0, slope=1)) +
    scale_colour_viridis_c()

ggplot(tmttpred) +
    geom_smooth(aes(x=nobs, y=val))

ggplot(tmttpred) +
    geom_jitter(aes(x=nobs, y=predabun)) +
    geom_smooth(aes(x=nobs, y=predabun))

#Replace TMTTs
bird.aru <- first %>%
    left_join(tmttpred) %>%
    mutate(abundance = as.numeric(ifelse(abundance=="TMTT", as.numeric(predabun), abundance)),
           distanceMethod = NA,
           distanceBand = NA) %>%
    dplyr::filter(!is.na(abundance)) %>%
    dplyr::select(id, organization, project, location, lat, lon, observer, datetime, durationMethod, distanceMethod, speciesCode, durationInterval, distanceBand, abundance)

#4. Save----
save(visit.aru, bird.aru, file="data/cleaned_data_aru_2022-10-15.Rdata")

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
modnames <- list(
    "0"="(Intercept)",
    "1"=c("(Intercept)", "jday"),
    "2"=c("(Intercept)", "tssr"),
    "3"=c("(Intercept)", "jday", "jday2"),
    "4"=c("(Intercept)", "tssr", "tssr2"),
    "5"=c("(Intercept)", "jday", "tssr"),
    "6"=c("(Intercept)", "jday", "jday2", "tssr"),
    "7"=c("(Intercept)", "jday", "tssr", "tssr2"),
    "8"=c("(Intercept)", "jday", "jday2", "tssr", "tssr2"),
    "9"=c("(Intercept)", "tsg"),
    "10"=c("(Intercept)", "tsg", "tsg2"),
    "11"=c("(Intercept)", "tsg", "tssr"),
    "12"=c("(Intercept)", "tsg", "tsg2", "tssr"),
    "13"=c("(Intercept)", "tsg", "tssr", "tssr2"),
    "14"=c("(Intercept)", "tsg", "tsg2", "tssr", "tssr2"))

#2. Load data----
load("data/cleaned_data_2022-10-06.Rdata")
load("data/cleaned_data_aru_2022-10-15.Rdata")

bird.all <- rbind(bird, bird.aru)
visit.all <- rbind(visit, visit.aru)

#3. Create design lookup table that describes duration method for each protocol----
#filter out duration methods that aren't appropriate for removal modelling (only have 1 time bin)
durdesign <- visit.all %>%
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
    left_join(bird.all %>%
                  dplyr::select(speciesCode) %>%
                  unique()) %>%
    arrange(speciesCode)

#5. Set up loop for species----
avail <- list()
for(i in 1:nrow(spp)){

    #6. Filter abundance data for species---
    # filter out observations with unknown duration method or interval
    # filter to observations with covariates
    bird.i <- bird.all %>%
        dplyr::filter(speciesCode==spp$speciesCode[i],
                      durationMethod %in% durdesign$durationMethod,
                      durationInterval!="UNKNOWN") %>%
        group_by(id, durationMethod, durationInterval) %>%
        summarize(abundance = sum(abundance)) %>%
        ungroup()

    #only model if there is data
    if(nrow(bird.i) > 0){

        #7. Replace abundance outliers (defined as 99% quantile across species) with 99% quantile of abundance----
        if(max(bird.i$abundance) > quantile(bird$abundance, 0.99)){
            bird.i <- bird.i %>%
                mutate(abundance = ifelse(abundance > quantile(abundance, 0.99), ceiling(quantile(abundance, 0.99)), abundance))
        }

        #8. Filter visit covariates----
        #Remove records with nas in covariates
        x <- visit.all %>%
            dplyr::filter(id %in% unique(bird.i$id),
                          !is.na(tssr),
                          !is.na(tsg),
                          !is.na(jday)) %>%
            arrange(id) %>%
            dplyr::select(id, durationMethod, jday, jday2, tssr, tssr2, tsg, tsg2)

        #9. Create design matrix----
        d <- x %>%
            dplyr::select(durationMethod) %>%
            left_join(durdesign, by="durationMethod") %>%
            dplyr::select(-durationMethod) %>%
            as.matrix()

        #10. Format abundance matrix----
        #add dummy variables for each of the columns to make sure the matrix is the right width
        y <- bird.i %>%
            dplyr::filter(id %in% x$id) %>%
            separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE, extra="drop", fill="right") %>%
            mutate(start = as.numeric(start),
                   end = as.numeric(str_sub(end, -100, -4))) %>%
            left_join(durdesign %>%
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
                rmvl <- mod[c("coefficients","vcov","nobs","loglik")]
                rmvl$p <- length(coef(mod))
                rmvl$names <- modnames[[j]]
            } else {
                rmvl <- mod
            }
            mod.list[[names(modnames)[j]]] <- rmvl
        }

        #13. Save model results to species list---
        avail[[i]] <- mod.list

    }

    print(paste0("Finished modelling species ", spp$English_Name[i], ": ", i, " of ", nrow(spp), " species"))

}

names(avail) <- spp$speciesCode

#14. Save out results----
save(durdesign, avail, file="results/availability_results_aru_2022-10-15.Rdata")

#E. PACKAGE####
load("results/availability_results_aru_2022-10-15.Rdata")
load("results/BAMCOEFS_QPAD_v4.rda")

#1. Remove species that didn't have enough data----
resDurOK <- avail[!sapply(avail, length)==0]
c(OK=length(resDurOK), failed=length(avail)-length(resDurOK), all=length(avail))

#2. Remove species where null model failed---
resDur <- resDurOK[!t(sapply(resDurOK, function(z) ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))[,1]==0]

#3. Create 0/1 table for model fit----
sra_mod <- t(sapply(resDur, function(z)
    ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))

#3. Adjust lists to include all species there's an estimate for----
tmp <- rownames(sra_mod)

sra_models <- matrix(0L, length(tmp), ncol(sra_mod))
dimnames(sra_models) <- list(tmp, colnames(sra_mod))
sra_models[rownames(sra_mod),] <- sra_mod

#4. Check for dropped factor levels in removal models----
for (spp in rownames(sra_mod)) {
    for (mid in colnames(sra_models)) {
        if (!inherits(resDur[[spp]][[mid]], "try-error")) {
            lcf <- length(resDur[[spp]][[mid]]$coefficients)
            lnm <- length(resDur[[spp]][[mid]]$names)
            if (lcf != lnm) {
                cat("SRA conflict for", spp, "model", mid,
                    "( len.coef =", lcf, ", len.name =", lnm, ")\n")
                sra_models[spp,mid] <- 0
            }
        } else {
            resDur[[spp]][[mid]] <- structure("Error", class = "try-error")
        }
        flush.console()
    }
}

#6. Exclude species with no null model----
sra_models[sra_models[,1]==0,] <- 0L

#7. Exclude species with phi estimate < 0.01----
phi0 <- sapply(resDur, function(z) exp(z[["0"]]$coefficients))
names(phi0) <- names(resDur)
sra_models[names(phi0)[phi0 < 0.01],] <- 0L

#8. Get number of models----
sra_nmod <- ncol(sra_mod)

#9. Get sample sizes----
sra_n <- numeric(length(tmp))
names(sra_n) <- tmp
sra_nn <- sapply(resDur, function(z) ifelse(inherits(z[["0"]], "try-error"),
                                            NA, z[["0"]]$nobs))
sra_n[names(sra_nn)] <- sra_nn

#10. exclude all models for species with < n.con observations----
n.con <- 25
sra_models[sra_n < n.con, ] <- 0L

#11. Exclude everything but null for species with n.con < observations < n.min----
n.min <- 75
sra_models[sra_n < n.min & sra_n >= n.con, 2:ncol(sra_models)] <- 0L

#12. ID species to keep----
spp <- tmp[rowSums(sra_models) > 0]
length(spp)

sra_models <- sra_models[spp,]
sra_models <- sra_models[.BAMCOEFS4$spp,]

sra_n <- sra_n[spp]

#13. Get number of parameters----
#use OVEN as template
sra_df <- sapply(resDur[["OVEN"]][1:sra_nmod], "[[", "p")

#14. Get estimates----
sra_estimates <- resDur[spp]

#15. Get species table----
tax <- read.csv("data/taxonomytable.csv")
tax <- tax[!duplicated(tax$Species_ID),]
rownames(tax) <- tax$Species_ID
spp_table <- data.frame(spp=spp,
                        scientific_name=tax[spp, "Scientific_Name"],
                        common_name=tax[spp, "English_Name"])
rownames(spp_table) <- spp
spp_table <- droplevels(spp_table)

#16. Get variable names for different models----
#use OVEN as template
sra_list <- sapply(sra_estimates[["OVEN"]], function(z) paste(z$names, collapse=" + "))

#16. Get loglik values---
sra_loglik <- sra_models
sra_loglik[] <- -Inf
for (i in .BAMCOEFS4$spp) { # species
    for (j in 1:sra_nmod) { # models
        if (sra_models[i,j] > 0)
            sra_loglik[i,j] <- resDur[[i]][[j]]$loglik
    }
}

#17. Get AIC values----
sra_aic <- sra_aicc <- sra_bic <- sra_loglik
sra_aic[] <- Inf
sra_aicc[] <- Inf
sra_bic[] <- Inf
for (i in .BAMCOEFS4$spp) {
    sra_aic[i,] <- -2*sra_loglik[i,] + 2*sra_df
    sra_aicc[i,] <- sra_aic[i,] + (2*sra_df*(sra_df+1)) / (sra_n[i]-sra_df-1)
    sra_bic[i,] <- -2*sra_loglik[i,] + log(sra_n[i])*sra_df
}

#18. Rank models----
sra_aicrank <- t(apply(sra_aic, 1, rank))*sra_models
sra_aicrank[sra_aicrank==0] <- NA

sra_aiccrank <- t(apply(sra_aicc, 1, rank))*sra_models
sra_aiccrank[sra_aiccrank==0] <- NA

sra_bicrank <- t(apply(sra_bic, 1, rank))*sra_models
sra_bicrank[sra_bicrank==0] <- NA

sra_aicbest <- apply(sra_aicrank, 1, function(z) colnames(sra_models)[which.min(z)])
sra_aiccbest <- apply(sra_aiccrank, 1, function(z) colnames(sra_models)[which.min(z)])
sra_bicbest <- apply(sra_bicrank, 1, function(z) colnames(sra_models)[which.min(z)])

#19. Set version----
version <- "4aru"

#20. Bundle----
bamcoefs <- list(spp=spp,
                 spp_table=spp_table,
                 edr_list=.BAMCOEFS4$edr_list,
                 sra_list=sra_list,
                 edr_models=.BAMCOEFS4$edr_models,
                 sra_models=sra_models,
                 edr_n=.BAMCOEFS4$edr_n,
                 sra_n=sra_n,
                 edr_df=.BAMCOEFS4$edr_df,
                 sra_df=sra_df,
                 edr_loglik=.BAMCOEFS4$edr_loglik,
                 sra_loglik=sra_loglik,
                 edr_aic=.BAMCOEFS4$edr_aic,
                 sra_aic=sra_aic,
                 edr_aicc=.BAMCOEFS4$edr_aicc,
                 sra_aicc=sra_aicc,
                 edr_bic=.BAMCOEFS4$edr_bic,
                 sra_bic=sra_bic,
                 edr_aicrank=.BAMCOEFS4$edr_aicrank,
                 sra_aicrank=sra_aicrank,
                 edr_aiccrank=.BAMCOEFS4$edr_aiccrank,
                 sra_aiccrank=sra_aiccrank,
                 edr_bicrank=.BAMCOEFS4$edr_bicrank,
                 sra_bicrank=sra_bicrank,
                 edr_aicbest=.BAMCOEFS4$edr_aicbest,
                 sra_aicbest=sra_aicbest,
                 edr_aiccbest=.BAMCOEFS4$edr_aiccbest,
                 sra_aiccbest=sra_aiccbest,
                 edr_bicbest=.BAMCOEFS4$edr_bicbest,
                 sra_bicbest=sra_bicbest,
                 edr_estimates=.BAMCOEFS4$edr_estimates,
                 sra_estimates=sra_estimates,
                 version=version)
.BAMCOEFS4aru <- list2env(bamcoefs)

save(.BAMCOEFS4aru, file="results/BAMCOEFS_QPAD_v4_aru.rda")

#F. COMPARE####
load("results/BAMCOEFS_QPAD_v4.rda")
load("results/BAMCOEFS_QPAD_v4_aru.rda")

my.theme <- theme_classic() +
    theme(text=element_text(size=12, family="Arial"),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(margin=margin(10,0,0,0)),
          axis.title.y=element_text(margin=margin(0,10,0,0)),
          axis.line.x=element_line(linetype=1),
          axis.line.y=element_line(linetype=1),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          plot.title=element_text(size=12, hjust = 0.5))

#Recall list of models
names4 <- data.frame(name = .BAMCOEFS4$sra_list,
                     mod = names(.BAMCOEFS4$sra_list))

#3. Sample size----
n4 <- data.frame(sra.n4=.BAMCOEFS4$sra_n,
                 species=.BAMCOEFS4$spp)
naru <- data.frame(sra.naru=.BAMCOEFS4aru$sra_n,
                 species=.BAMCOEFS4aru$spp)

n4aru <- full_join(n4, naru) %>%
    mutate(sra.n4 = as.numeric(sra.n4),
           sra.naru = ifelse(is.na(sra.naru), 0, sra.naru),
           sra.n4 = ifelse(is.na(sra.n4), 0, sra.n4),
           sra.n = sra.naru/sra.n4) %>%
    mutate(version = case_when(sra.n4==0 ~ "V4",
                               sra.naru==0 ~ "V4 + ARU",
                               !is.na(sra.n) ~ "Both"))

sum(n4aru$sra.n4)
sum(n4aru$sra.naru)

ggplot(n4aru) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x=sra.n4, y=sra.naru, fill=version), pch=21, alpha = 0.5, size=4) +
    xlab("V4 sample size") +
    ylab("V4 + ARU sample size") +
    scale_fill_manual(values=c("grey80", "blue", "orange"), name="") +
    my.theme

ggsave(filename="figures/ARU_samplesize.jpeg", width =7, height=6)

#4. Null estimates----
est4 <-data.frame()
for(i in 1:length(.BAMCOEFS4$spp)){
    sra4 <- exp(.BAMCOEFS4$sra_estimates[[i]]$`0`$coefficients)
    edr4 <- exp(.BAMCOEFS4$edr_estimates[[i]]$`0`$coefficients)
    est4 <- rbind(est4,
                  data.frame(sra4=sra4, edr4=edr4, species=.BAMCOEFS4$spp[i]))
}
rownames(est4) <- NULL

estaru <-data.frame()
for(i in 1:length(.BAMCOEFS4$spp)){
    sraaru <- exp(.BAMCOEFS4aru$sra_estimates[[i]]$`0`$coefficients)
    edraru <- exp(.BAMCOEFS4aru$edr_estimates[[i]]$`0`$coefficients)
    estaru <- rbind(estaru,
                  data.frame(sraaru=sraaru, edraru=edraru, species=.BAMCOEFS4$spp[i]))
}
rownames(estaru) <- NULL

est4aru <- full_join(est4, estaru) %>%
    mutate(edraru = ifelse(is.na(edraru), 0, edraru),
           sraaru = ifelse(is.na(sraaru), 0, sraaru),
           edr4 = ifelse(is.na(edr4), 0, edr4),
           sra4 = ifelse(is.na(sra4), 0, sra4)) %>%
    mutate(version = case_when(edraru==0 ~ "V4",
                               edr4==0 ~ "V4 + ARU",
                               !is.na(edr4) ~ "Both")) %>%
    left_join(n4aru)

est4aru$SampleSizeRatio <- est4aru$sra.n
ggplot(est4aru) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x=sra4, y=sraaru, size = SampleSizeRatio), pch=21, alpha = 0.5, fill="grey80") +
    xlab("V4 availability estimate (phi)") +
    ylab("V4 + ARU availability estimate (phi)") +
    my.theme

ggsave(filename="figures/ARU_phi.jpeg", width =7, height=6)

ggplot(est4aru) +
    geom_abline(intercept = 0, slope = 1) +
    geom_text(aes(x=log(sra4), y=log(sraaru), label=species)) +
    xlab("V4 availability estimate (phi)") +
    ylab("V4 + ARU availability estimate (phi)") +
    my.theme

ggsave(filename="figures/ARU_phi_species.jpeg", width =7, height=6)
