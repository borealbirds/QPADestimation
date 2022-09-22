# ---
# title: "QPAD estimation - get data"
# author: "Elly Knight"
# created: "July 24, 2022"
# ---

library(tidyverse)
library(wildRtrax)
library(data.table)

#NOTE: DOCUMENT TAXONOMY####

#1. Get list of projects from WildTrax----
wt_auth()

#sensor = PC gives all ARU and point count projects
project.list <- wt_get_download_summary(sensor_id = 'PC')

#2. Convert to a plain dataframe (not sure why it's a named list)
projects <- data.frame(project = as.character(project.list$project),
                       project_id = as.numeric(project.list$project_id),
                       sensorId = as.character(project.list$sensorId),
                       tasks = as.numeric(project.list$tasks),
                       status = as.character(project.list$status))

#3. Loop through projects to download data----
dat.list <- list()
for(i in 1:nrow(projects)){

    try(dat.list[[i]] <- wt_download_report(project_id = projects$project_id[i], sensor_id = projects$sensorId[i], cols_def = F, weather_cols = T))

    print(paste0("Finished dataset ", projects$project[i], " : ", i, " of ", nrow(projects), " projects"))

}

#4. Collapse list----
raw <- rbindlist(dat.list, fill=TRUE)

#5. Save date stamped data & project list----
save(raw, projects, file=paste0("data/wildtrax_data_", Sys.Date(), ".Rdata"))
