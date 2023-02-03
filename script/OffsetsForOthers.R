library(tidyverse)
library(lubridate)
library(data.table)

#A. GET DATA#####

#1. Load wildRtrax package----
library(wildRtrax)

#2. Login to WildTrax----
config <- "script/login.R"
source(config)
wt_auth()

#3. List of available projects----
project.list <- wt_get_download_summary(sensor_id = 'PC')
projects <- data.frame(project = as.character(project.list$project),
                       project_id = as.numeric(project.list$project_id),
                       sensorId = as.character(project.list$sensorId),
                       tasks = as.numeric(project.list$tasks),
                       status = as.character(project.list$status))

#4. Pick projects to download----
dl <- projects %>% 
  dplyr::filter(project_id %in% c(801, 1095, 1384, 1385))

#5. Download data----
dat.list <- list()
error.log <- data.frame()
for(i in 1:nrow(dl)){
  
  #Do each sensor type separately because the reports have different columns and we need different things for each sensor type
  if(dl$sensorId[i]=="ARU"){
    
    #Get summary report
    report.try <- try(wt_download_report(project_id = dl$project_id[i], sensor_id = dl$sensorId[i], weather_cols = F, report = "summary"))
    
    #Get task report for ARU model
    task.try <- try(wt_download_report(project_id = dl$project_id[i], sensor_id = dl$sensorId[i], weather_cols = F, report = "task"))
    
    if(class(report.try)=="data.frame"){
      dat.try <- report.try %>% 
        left_join(task.try  %>% 
                    dplyr::select("organization", "project_name", "location", "recording_date", "longitude", "latitude", "method", "status", "observer", "observer_id", "equipment_used", "buffer"))
    }
  }
  
  if(dl$sensorId[i]=="PC"){
    
    dat.try <- try(wt_download_report(project_id = dl$project_id[i], sensor_id = dl$sensorId[i], weather_cols = F, report="report"))
    
  }
  
  if(class(dat.try)=="data.frame"){
    dat.list[[i]] <- dat.try
  }
  
  #Log projects that error
  if(class(dat.try)!="data.frame"){
    error.log <- rbind(error.log, 
                       dl[i,])
    
  }
  
  print(paste0("Finished dataset ", dl$project[i], " : ", i, " of ", nrow(dl), " projects"))
  
}

#6. Collapse list----
#standardize column names between sensor types
dat.raw <- rbindlist(dat.list, fill=TRUE) %>% 
  left_join(dl %>% 
              dplyr::rename(project_status = status))

#B. PREPARE QPAD###########

#1. Load basic packages----
library(maptools)
library(intrval)
library(raster)

#2. Set file path for data & results----
root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/OffsetCalculation/VenierFeb2023"

#3. Install most recent version of QPAD----
#library(devtools)
#install_github("borealbirds/QPAD@dev-v4", force=TRUE)

#4. Load package & offsets----
library(QPAD)
load_BAM_QPAD(4)

#5. Set WD to qpad-offsets package----
setwd("C:/Users/Elly Knight/Documents/BAM/Projects/QPAD/qpad-offsets")

#6. Read raster data----
rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
rd1 <- raster("./data/seedgrow.tif")
crs <- proj4string(rtree)

#7. Source functions----
source("functions.R")

#C. WRANGLE DATA####

#1. Wrangle----
spp <- getBAMspecieslist()

dat.use <- dat.raw %>% 
  dplyr::select(organization, project, sensorId, location, latitude, longitude, observer, date, distanceMethod, durationMethod) %>% 
  unique() %>% 
  left_join(dat.raw %>% 
              mutate(abundance = as.numeric(abundance)) %>% 
              dplyr::filter(speciesCode %in% spp,
                            !is.na(abundance)) %>% 
              dplyr::select(organization, project, sensorId, location, latitude, longitude, observer, date, distanceMethod, durationMethod, speciesCode, abundance)) %>% 
  mutate(speciesCode = ifelse(is.na(speciesCode), "NONE", speciesCode),
         abundance = ifelse(is.na(abundance), 0, abundance),
         distanceMethod = case_when(distanceMethod=="0m-100m-INF" ~ Inf,
                         distanceMethod=="0m-400m" ~ 400),
         durationMethod = case_when(durationMethod=="0-3min" ~ 3,
                         durationMethod=="0-5min" ~ 5),
         time = str_sub(date, 12, 16),
         date = as.character(as.Date(date))) %>% 
  rename(species = speciesCode, tagmeth = sensorId, lat = latitude, lon = longitude, dur = durationMethod, dis = distanceMethod) %>% 
  dplyr::filter(lat >=39, lat <= 69) %>% 
  pivot_wider(names_from=species, values_from=abundance, values_fn=sum, values_fill=0, names_sort=TRUE) %>% 
  dplyr::select(-NONE) %>% 
  mutate(id = row_number())

#2. Make prediction df----
x <- make_x(dat.use, tz="local")

#3. Filter out NAs----
data <- dat.use[!(is.na(x$TSSR) | is.na(x$LCC4) | is.na(x$DSLS)),] %>% 
  mutate(id = row_number())
data <- data[,c(197, 1:196)]
x.use <- x %>% 
  dplyr::filter(!is.na(TSSR), !is.na(LCC4), !is.na(DSLS))
nrow(data)==nrow(x.use)

#C. CREATE OFFSETS#######

#1. Get list of species----
SPP <- colnames(data)[-c(1:12)]

#2. Set up output----
OFF <- matrix(0, nrow(x.use), length(SPP))
rownames(OFF) <- data$id
colnames(OFF) <- SPP

#3. Make OFF----
for (spp in SPP) {
  cat(spp, "\n")
  flush.console()
  o <- make_off(spp, x.use, useMethod="y")
  OFF[,spp] <- round(o$offset, 4)
}
offsets <- data.frame(id=data$id) %>% 
  cbind(data.frame(OFF))

#D. PACKAGE####
save(data, offsets, file=file.path(root, paste0("Data&Offsets_Venier_", Sys.Date(), ".Rdata")))
