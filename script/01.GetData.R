# ---
# title: "QPAD estimation - get data"
# author: "Elly Knight"
# created: "July 24, 2022"
# updated: "November 6, 2022"
# ---

#TO DO: UPDATE EVERYTHING TO COV NAMES IN PREVIOUS QPAD VERSIONS####

#NOTES################################

#QPAD V4 and the National models V4.1 were updated concurrently and use the same base dataset. The following code is similar to the code for the national models and the data is stored in the national models folder on Google Drive

#The "projectInstructions.csv" file is a list of all projects currently in WildTrax should not be used in ABMI models (instructions=="DO NOT USE"). This file should be updated for each iteration of in collaboration with Erin Bayne. Future versions of this spreadsheet can hopefully be derived by a combination of organization and a google form poll for consent from other organizations. Note this category also includes all ABMI projects with inaccurate (i.e., buffered) coordinates.

#The "projectInstructions.csv" file also contains information on which ARU projects are processed for a single species or taxa (instructions=="DO NOT USE") and therefore those visits should only be used for models for the appropriate taxa. This file should be updated for each iteration of the national models in collaboration with Erin Bayne. These projects are currently not included in the models.

#There are a handful of projects that are not downloading properly via wildRtrax. An issue is open on this. These projects are listed in the error.log object. These files should be downloaded manually.

#The column "sensor" currently only differentiates between ARU & human point count data types. Future versions should consider differentiating between SM2 ARU data types and other ARU data types due to differences in the perceptibility of these two approaches, either via QPAD or a correction factor.

#The replace TMTTs script will be replaced by a wildRtrax function in the near future.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(wildRtrax) #to download data from wildtrax
library(data.table) #for binding lists into dataframes
library(lubridate) #date wrangling

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_RshProjs/PopnStatus/QPAD"

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

#3. Loop through projects to download data----
dat.list <- list()
error.log <- data.frame()
for(i in 1:nrow(projects)){
  
  #Do each sensor type separately because the reports have different columns and we need different things for each sensor type
  if(projects$sensorId[i]=="ARU"){
    
    #Get summary report
    report.try <- try(wt_download_report(project_id = projects$project_id[i], sensor_id = projects$sensorId[i], weather_cols = F, report = "summary"))
    
    #Get task report for ARU model
    task.try <- try(wt_download_report(project_id = projects$project_id[i], sensor_id = projects$sensorId[i], weather_cols = F, report = "task"))
    
    if(class(report.try)=="data.frame"){
      dat.try <- report.try %>% 
        left_join(task.try  %>% 
                    dplyr::select("organization", "project_name", "location", "recording_date", "longitude", "latitude", "method", "status", "observer", "observer_id", "equipment_used", "buffer"))
    }
  }
  
  if(projects$sensorId[i]=="PC"){
    
    dat.try <- try(wt_download_report(project_id = projects$project_id[i], sensor_id = projects$sensorId[i], weather_cols = F, report="report"))
    
  }
  
  if(class(dat.try)=="data.frame"){
    dat.list[[i]] <- dat.try
  }
  
  #Log projects that error
  if(class(dat.try)!="data.frame"){
    error.log <- rbind(error.log, 
                       projects[i,])
    
  }
  
  print(paste0("Finished dataset ", projects$project[i], " : ", i, " of ", nrow(projects), " projects"))
  
}

#4. Go download error log projects from wildtrax.ca----

#5. Read in error projects----
error.files <- list.files(file.path(root, "errorFiles"), full.names = TRUE)

raw.error <- data.frame()
for(i in 1:length(error.files)){
  raw.error <- read.csv(error.files[i]) %>% 
    rbind(raw.error)
}

#6. Collapse list----
#standardize column names between sensor types
all.wt <- rbindlist(dat.list, fill=TRUE)  %>% 
  rbind(raw.error, fill=TRUE) %>% 
  mutate(project = ifelse(is.na(project), project_name, project),
         speciesCode = ifelse(is.na(speciesCode), species_code, speciesCode),
         date = ifelse(is.na(date), recording_date, date),
         buffer = ifelse(is.na(buffer), bufferRadius.m., buffer)) %>%
  dplyr::select(-project_name, -recording_date, -species_code, -bufferRadius.m.) %>% 
  left_join(projects %>% 
              dplyr::rename(project_status = status))

#7. Filter out projects that shouldn't be used----
#nothing in BU training & all "DO NOT USE" projects in projectInventory file
#filter out 'NONE' method ARU projects later after this field is parsed out

instructions <- read.csv(file.path(root, "projectInventory", "projectInstructions.csv"))

raw.wt <- all.wt %>% 
  anti_join(instructions) %>%
  dplyr::filter(organization!="BU-TRAINING")

#9. Save date stamped data & project list----
save(raw.wt, projects, error.log, file=paste0(root, "/wildtrax_raw_", Sys.Date(), ".Rdata"))
