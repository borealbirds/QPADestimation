# ---
# title: "QPAD estimation - clean bird data"
# author: "Elly Knight"
# created: "September 22, 2022"
# updated: "November 6, 2022"
# ---

library(tidyverse) #basic data wrangling
library(lubridate) #date and time wrangling

root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/"

#1. Load in use  & visit dataset----
load(file.path(root,"qpadv4_raw.Rdata"))
load(file.path(root,"qpadv4_visit.Rdata"))

#2. Filter, tidy, create foreign key for visit table----
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
  group_by(id, source, project, sensor, singlesp, location, buffer, lat, lon, year, date,  observer, species, abundance, individual, isSeen, isHeard) %>% 
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
  dplyr::select(id, source, project, sensor, singlesp, location, buffer, lat, lon, year, date, observer, distanceMethod, durationMethod, tagMethod, distanceBand, durationInterval, species, abundance, isSeen, isHeard)

#3. Replace TMTTs with predicted abundance----
tmtt <- read.csv("C:/Users/Elly Knight/Documents/ABMI/Projects/TMTT/data/tmtt_predictions.csv") %>% 
  rename(species = species_code)

user <- read.csv("C:/Users/Elly Knight/Documents/ABMI/Projects/TMTT/data/app_user.csv") %>% 
  rename(observer = user_name) %>% 
  dplyr::select(observer, user_id)

bird <- first %>% 
  dplyr::filter(abundance=="TMTT") %>% 
  left_join(user) %>% 
  mutate(user_id = ifelse(is.na(user_id), observer, user_id))%>% 
  mutate(boot = round(runif(max(row_number()), 1, 100)),
         species = ifelse(species %in% tmtt$species, species, "species"),
         user_id = as.integer(ifelse(user_id %in% tmtt$user_id, user_id, 0))) %>% 
  data.frame() %>% 
  left_join(tmtt) %>% 
  mutate(abundance = round(pred)) %>% 
  dplyr::select(colnames(first)) %>% 
  rbind(first %>% 
          dplyr::filter(abundance!="TMTT") %>% 
          mutate(abundance = as.numeric(abundance))) 


#3. Add species list----
species <- read.csv(file.path(root, "singing-species.csv")) %>%
    rename(species = Species_ID)

#4. Save----
save(visit, bird, species,  file="G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/qpadv4_clean.Rdata")
