library(tidyverse) #basic data wrangling

root <- "G:/Shared drives/BAM_RshProjs/PopnStatus/QPAD"

#1. Load data----
load(file.path(root, "Data", "qpadv4_clean.Rdata"))

#2. Create anonmyous project & location id----
project <- data.frame(project = unique(visit$project)) %>% 
  mutate(projectid = row_number()) %>% 
  left_join(visit %>% 
              dplyr::select(project, location) %>% 
              unique()) %>% 
  group_by(project) %>% 
  mutate(locationid = row_number()) %>% 
  ungroup()

#3. Strip to just data for detectability modelling----
#remove ARU tag method 'NONE'
#fill in distance method for ARU data
det <- bird %>% 
  left_join(project) %>% 
  dplyr::select(projectid, TM, locationid, date, distanceMethod, durationMethod, distanceBand, durationInterval, species, abundance, isSeen, isHeard) %>% 
  rename(sensorMethod = TM) %>% 
  dplyr::filter(sensorMethod!="ARU-NONE") %>% 
  mutate(distanceMethod = ifelse(sensorMethod!="PC", Inf, distanceMethod))

#4. Save-----
saveRDS(det, file.path(root, "Data", "BAM_QPADData_DavidIles.rds"))
saveRDS(det, "G:/Shared Drives/BAM/Data/StuffBAM_QPADData_DavidIles.rds")