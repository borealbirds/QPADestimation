# ---
# title: "QPAD estimation - clean bird data"
# author: "Elly Knight"
# created: "September 22, 2022"
# updated: "November 6, 2022"
# ---

library(tidyverse) #basic data wrangling
library(lubridate) #date and time wrangling

#1. Load in use  & visit dataset----
load("G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/qpadv4_dat_2022-11-06.Rdata")
load("G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/qpadv4_visit_2022-11-06.Rdata")

#2. Filter, tidy, create foreign key for visit table----
#Ensure there's visit data
#Remove UNSPs
dat <- use %>% 
    mutate(datetime = ymd_hms(date),
           julian = yday(datetime),
           id = paste(location, observer, datetime)) %>%
    rename(lat = latitude, lon = longitude) %>% 
    dplyr::filter(str_sub(speciesCode, 1, 2)!="UN",
                  id %in% visit$id)

#2. Filter to just first detection per individual and bin in minutes----
#bin in 30 s bins for 60 s surveys
first <- dat %>%
  dplyr::filter(method=="ARU") %>% 
  group_by(id, organization, project, location, lat, lon, observer, datetime, durationMethod, speciesCode, abundance, individual_appearance_order) %>%
  mutate(first_tag = min(tag_start_s)) %>%
  ungroup() %>%
  dplyr::filter(tag_start_s == first_tag) %>%
  dplyr::select(-first_tag) %>%
  dplyr::filter(!abundance %in% c("CI 1", "CI 2", "CI 3", "N/A", "0")) %>%
  mutate(end=ifelse(durationMethod=="0-0.5-1min", ceiling(tag_start_s/30), ceiling(tag_start_s/60)),
         end = ifelse(end==0, 1, end),
         durationInterval = case_when(durationMethod=="0-0.5-1min" & end==1 ~ "0-0.5min",
                                      durationMethod=="0-0.5-1min" & end==2 ~ "0.5-1min",
                                      !is.na(durationMethod) ~ paste0(end-1, "-", end, "min"))) %>% 
    dplyr::select(colnames(dat)) %>% 
    rbind(dat %>% 
            dplyr::filter(method=="PC")) %>%
  dplyr::select(id, dataset, organization, project, method, location, lat, lon, observer, datetime, distanceMethod, durationMethod, speciesCode, distanceBand, durationInterval, abundance)

#3. Replace TMTTs----
#99 percentile for combination of species & observer
tmtttbl <- first %>%
  dplyr::filter(abundance!="TMTT",
                method=="ARU") %>%
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

#Replace TMTTs
bird <- first %>%
  left_join(tmttpred) %>%
  mutate(abundance = as.numeric(ifelse(abundance=="TMTT", as.numeric(predabun), abundance)),
         distanceMethod = NA,
         distanceBand = NA) %>%
  dplyr::filter(!is.na(abundance)) %>%
  dplyr::select(id, dataset, organization, method, project, location, lat, lon, observer, datetime, durationMethod, distanceMethod, speciesCode, durationInterval, distanceBand, abundance)

#3. Add species list----
species <- read.csv("G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/singing-species.csv") %>%
    rename(speciesCode = Species_ID)

#4. Save----
save(visit, bird, species,  file="G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/Data/qpadv4_clean_2022-11-06.Rdata")
