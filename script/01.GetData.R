# ---
# title: "QPAD estimation - get data"
# author: "Elly Knight"
# created: "July 24, 2022"
# updated: "November 6, 2022"
# ---

#TO DO: UPDATE EVERYTHING TO COV NAMES IN PREVIOUS QPAD VERSIONS####

#NOTES################################

#QPAD V4 and the National models V4.1 were updated concurrently and use the same base dataset. The following code is similar to the code for the national models and the data is stored in the national models folder on Google Drive

#The "BAMProjects_WildTrax.csv" file is a list of all projects currently in WildTrax that BAM can use in the national models. This file should be updated for each iteration of the national models in collaboration with Erin Bayne. Future versions of this spreadsheet can hopefully be derived by a combination of organization and a google form poll for consent from other organizations.

#The "BAMProjects_WildTrax.csv" file also contains information on which ARU projects are processed for a single species or taxa and therefore those visits should only be used for models for the appropriate taxa. This file should be updated for each iteration of the national models in collaboration with Erin Bayne.

#There are a handful of projects that are not downloading properly via wildRtrax. An issue is open on this. These projects are listed in the error.log object.

#BAM patch was compiled by Melina Houle to include datasets that are not yet loaded into WildTrax. Future iterations of the national models should not need this patch, with the exception of the BBS data which is likely too large to ever be uploaded in WT.

#The column "sensor" currently only differentiates between ARU & human point count data types. Future versions should consider differentiating between SM2 ARU data types and other ARU data types due to differences in the perceptibility of these two approaches, either via QPAD or a correction factor.

#The replace TMTTs script will be replaced by a wildRtrax function in the near future.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(wildRtrax) #to download data from wildtrax
library(data.table) #for binding lists into dataframes
library(lubridate) #date wrangling
library(auk) #eBird wrangling

#2. Set root path for data on google drive----
root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/NationalModelsV4.1/PointCount/"

#A. DOWNLOAD DATA FROM WILDTRAX#######################

#1. Login to WildTrax----
config <- "script/login.R"
source(config)

#1. Get list of projects from WildTrax----
wt_auth()

#sensor = PC gives all ARU and point count projects
project.list <- wt_get_download_summary(sensor_id = 'PC')

#2. Convert to a plain dataframe----
projects <- data.frame(project = as.character(project.list$project),
                       project_id = as.numeric(project.list$project_id),
                       sensorId = as.character(project.list$sensorId),
                       tasks = as.numeric(project.list$tasks),
                       status = as.character(project.list$status))

#3. Filter by list of projects we can use----
#Fix names to match report
project.names <- data.frame(project_id = c(325, 432, 856, 977, 978),
                            newname = c("Tłı̨chǫ Winter Road CWS Northern Region 2019",
                                        "Ts’udé Nilįné Tuyeta PA CWS Northern Region 2020",
                                        "Ts’udé Nilįné Tuyeta PA CWS Northern Region 2021",
                                        "Nááts'įhch'oh NPR 2018 CWS Northern Region",
                                        "Nááts'įhch'oh NPR 2019 CWS Northern Region"))

use <- read.csv(file.path(root, "BAMProjects_WildTrax.csv"))

projects.wt <- projects %>% 
  dplyr::filter(project_id %in% use$project_id) %>% 
  left_join(project.names) %>% 
  mutate(project = ifelse(is.na(newname), project, newname)) %>% 
  dplyr::select(-newname)

#4. Loop through projects to download data----
dat.list <- list()
error.log <- data.frame()
for(i in 1:nrow(projects.wt)){
  
  #Do each sensor type separately because the reports have different columns and we need different things for each sensor type
  if(projects.wt$sensorId[i]=="ARU"){
    dat.try <- try(wt_download_report(project_id = projects.wt$project_id[i], sensor_id = projects.wt$sensorId[i], weather_cols = F, report = "summary"))
  }
  
  if(projects.wt$sensorId[i]=="PC"){
    dat.try <- try(wt_download_report(project_id = projects.wt$project_id[i], sensor_id = projects.wt$sensorId[i], weather_cols = F, report="report"))
  }
  
  if(class(dat.try)=="data.frame"){
    dat.list[[i]] <- dat.try
  }
  
  #Log projects that error
  if(class(dat.try)!="data.frame"){
    error.log <- rbind(error.log, 
                       projects.wt[i,])
    
  }
  
  print(paste0("Finished dataset ", projects.wt$project[i], " : ", i, " of ", nrow(projects.wt), " projects"))
  
}

#5. Collapse list----
#standardize column names between sensor types
#fix special character project names
report.names <- data.frame(project=c("EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2016",
                                     "EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2019",
                                     "TÅ‚Ä±Ì¨chÇ« Winter Road CWS Northern Region 2019",
                                     "Tsâ€™udÃ© NilÄ¯nÃ© Tuyeta PA CWS Northern Region 2020",
                                     "ForÃªt Montmorency long-term bird survey 2000",
                                     "RÃ©gularisation du Lac KÃ©nogami EIA 2001",
                                     "Tsâ€™udÃ© NilÄ¯nÃ© Tuyeta PA CWS Northern Region 2021",
                                     "NÃ¡Ã¡ts'Ä¯hch'oh NPR 2018 CWS Northern Region",
                                     "NÃ¡Ã¡ts'Ä¯hch'oh NPR 2019 CWS Northern Region",
                                     "EdÃ©hzhÃ­e National Wildlife Area CWS Northern Region 2021",
                                     "Thaidene NÃ«nÃ© PA CWS Northern Region 2022"),
                           newname=c("Edéhzhíe National Wildlife Area CWS Northern Region 2016",
                                     "Edéhzhíe National Wildlife Area CWS Northern Region 2019",
                                     "Tchicho Winter Road CWS Northern Region 2019",
                                     "Ts’udé Niliné Tuyeta PA CWS Northern Region 2020",
                                     "Forêt Montmorency long-term bird survey 2000",
                                     "Régularisation du Lac Kénogami EIA 2001",
                                     "Ts’udé Niliné Tuyeta PA CWS Northern Region 2021",
                                     "Nááts'ihch'oh NPR 2018 CWS Northern Region",
                                     "Nááts'ihch'oh NPR 2019 CWS Northern Region",
                                     "Edéhzhíe National Wildlife Area CWS Northern Region 2021",
                                     "Thaidene Nëné PA CWS Northern Region 2022"))

raw.wt <- rbindlist(dat.list, fill=TRUE) %>%
  mutate(project = ifelse(is.na(project), project_name, project),
         speciesCode = ifelse(is.na(speciesCode), species_code, speciesCode),
         date = ifelse(is.na(date), recording_date, date)) %>%
  dplyr::select(-project_name, -recording_date, -species_code) %>% 
  left_join(report.names) %>% 
  mutate(project = ifelse(is.na(newname), project, newname)) %>% 
  dplyr::select(-newname)

#7. Save date stamped data & project list----
save(raw.wt, projects.wt, error.log, file=paste0(root, "/wildtrax_raw_", Sys.Date(), ".Rdata"))
load(file.path(root, "wildtrax_raw_2022-11-22.Rdata"))

#B. GET PATCH DATA###############################

#1. BAM patch----
raw.bam <- readRDS(file.path(root, "pc_patch.rds"))

#C. HARMONIZE###############################

#1. Set desired columns----
colnms <- c("source", "project", "sensor", "singlesp", "location", "buffer", "lat", "lon", "year", "date", "observer", "ARUMethod", "durationMethod", "distanceMethod", "species", "abundance", "individual", "durationInterval", "distanceBand", "tagStart", "isHeard", "isSeen")

#2. Wrangle wildtrax data-----
#identify single species projects
ssp <- read.csv(file.path(root, "BAMProjects_WildTrax.csv")) %>% 
  dplyr::filter(single.species=="y")

use.wt <- raw.wt %>% 
  rename(species = speciesCode, buffer = bufferRadius.m., lat = latitude, lon = longitude, individual = individual_appearance_order, ARUMethod = method, tagStart = tag_start_s) %>% 
  full_join(projects.wt %>% 
              rename(sensor = sensorId) %>% 
              dplyr::select(project_id, project, sensor)) %>% 
  mutate(source = "WildTrax", 
         date = ymd_hms(date),
         year = year(date),
         singlesp = ifelse(project_id %in% ssp$project_id, "y", "n")) %>% 
  dplyr::select(all_of(colnms))


#3. Wrangle BAM patch data----
#wrangle distance and duration maximums
#remove counts with odd duration method entries
#replace all unknown dates (MN-BBATLAS, NEFBMP2012-19) with June 15 of 2012
use.bam <- raw.bam %>% 
  dplyr::filter(!durationMethod %in% c("", " during the 10 minutes.", " seemed to bring food then brooded; assume nestlings stil", " timeperiod C")) %>% 
  mutate(source = "BAM",
         sensor = "PC",
         singlesp = "n",
         date = ymd_hms(date),
         year = year(date),
         observer = NA,
         ARUMethod= NA,
         individual = NA,
         tagStart = NA) %>% 
  rename(buffer = 'bufferRadius(m)', lat = latitude, lon = longitude, species=speciesCode) %>% 
  dplyr::select(all_of(colnms))

#D. PUT DATASETS TOGETHER & SAVE######################

#1. Put together----
use <- rbind(use.wt, use.bam) %>% 
  mutate(id = paste(project, location, lat, lon, observer, date, ARUMethod)) %>% 
  dplyr::filter(!is.na(date))

#2. Save----
save(use, file="G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/qpadv4_raw.Rdata")
