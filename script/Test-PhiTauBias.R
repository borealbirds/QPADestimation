# PREAMBLE####

#1. Load packages----
library(tidyverse) #basic data wrangling
library(detect) #removal models
library(data.table) #collapse list to dataframe
library(gridExtra)

options(dplyr.summarise.inform = FALSE, scipen=9999)

#2. Data root----
root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/"

#3. Load data----
load(file.path(root, "Data", "qpadv4_clean.Rdata"))

#Set factor levels
visit$TM <- factor(visit$TM, levels=c("PC", "ARU-1SPT", "ARU-1SPM"))
#visit$sensor <- factor(visit$sensor, levels=c("PC", "ARU"))

#4. Load V4 null estimates----
est34 <- read.csv(file.path(root, "Results/QPADV3V4NullEstimates.csv"))

#A. AVAILABILITY####

#1. Create list of models----
#JDAY = day of year as a decimal between 0 and 1
#TSSR = time since sunrise as a decimal between 0 and 1
#DSLS = ("datys since local spring") days since start of seedgrowth from seedgrow layers
#TM = TM ("PC" = point count, "ARU-1SPM" = 1 tag per minute, "ARU-1SPT" = 1 tag per task - ie.e., recording)

#1a. Lists for species that only have data for one method----
mods <- list(
  ~ 1,
  ~ JDAY,
  ~ TSSR,
  ~ poly(JDAY, 2),
  ~ poly(TSSR, 2),
  ~ JDAY + TSSR,
  ~ poly(JDAY, 2) + TSSR,
  ~ JDAY + poly(TSSR, 2),
  ~ poly(JDAY, 2) + poly(TSSR, 2),
  ~ DSLS,
  ~ poly(DSLS, 2),
  ~ DSLS + TSSR,
  ~ poly(DSLS, 2) + TSSR,
  ~ DSLS + poly(TSSR, 2),
  ~ poly(DSLS, 2) + poly(TSSR, 2),
  ~ TM,
  ~ TM + JDAY,
  ~ TM + TSSR,
  ~ TM + poly(JDAY, 2),
  ~ TM + poly(TSSR, 2),
  ~ TM + JDAY + TSSR,
  ~ TM + poly(JDAY, 2) + TSSR,
  ~ TM + JDAY + poly(TSSR, 2),
  ~ TM + poly(JDAY, 2) + poly(TSSR, 2),
  ~ TM + DSLS,
  ~ TM + poly(DSLS, 2),
  ~ TM + DSLS + TSSR,
  ~ TM + poly(DSLS, 2) + TSSR,
  ~ TM + DSLS + poly(TSSR, 2),
  ~ TM + poly(DSLS, 2) + poly(TSSR, 2))
names(mods) <- 0:29
modnames <- list(
  "0"=c("(Intercept)"),
  "1"=c("(Intercept)", "JDAY"),
  "2"=c("(Intercept)", "TSSR"),
  "3"=c("(Intercept)", "JDAY", "JDAY2"),
  "4"=c("(Intercept)", "TSSR", "TSSR2"),
  "5"=c("(Intercept)", "JDAY", "TSSR"),
  "6"=c("(Intercept)", "JDAY", "JDAY2", "TSSR"),
  "7"=c("(Intercept)", "JDAY", "TSSR", "TSSR2"),
  "8"=c("(Intercept)", "JDAY", "JDAY2", "TSSR", "TSSR2"),
  "9"=c("(Intercept)", "DSLS"),
  "10"=c("(Intercept)", "DSLS", "DSLS2"),
  "11"=c("(Intercept)", "DSLS", "TSSR"),
  "12"=c("(Intercept)", "DSLS", "DSLS2", "TSSR"),
  "13"=c("(Intercept)", "DSLS", "TSSR", "TSSR2"),
  "14"=c("(Intercept)", "DSLS", "DSLS2", "TSSR", "TSSR2"),
  "15"=c("(Intercept)", "TM1SPT", "TM1SPM"),
  "16"=c("(Intercept)", "TM1SPT", "TM1SPM", "JDAY"),
  "17"=c("(Intercept)", "TM1SPT", "TM1SPM", "TSSR"),
  "18"=c("(Intercept)", "TM1SPT", "TM1SPM", "JDAY", "JDAY2"),
  "19"=c("(Intercept)", "TM1SPT", "TM1SPM", "TSSR", "TSSR2"),
  "20"=c("(Intercept)", "TM1SPT", "TM1SPM", "JDAY", "TSSR"),
  "21"=c("(Intercept)", "TM1SPT", "TM1SPM", "JDAY", "JDAY2", "TSSR"),
  "22"=c("(Intercept)", "TM1SPT", "TM1SPM", "JDAY", "TSSR", "TSSR2"),
  "23"=c("(Intercept)", "TM1SPT", "TM1SPM", "JDAY", "JDAY2", "TSSR", "TSSR2"),
  "24"=c("(Intercept)", "TM1SPT", "TM1SPM", "DSLS"),
  "25"=c("(Intercept)", "TM1SPT", "TM1SPM", "DSLS", "DSLS2"),
  "26"=c("(Intercept)", "TM1SPT", "TM1SPM", "DSLS", "TSSR"),
  "27"=c("(Intercept)", "TM1SPT", "TM1SPM", "DSLS", "DSLS2", "TSSR"),
  "28"=c("(Intercept)", "TM1SPT", "TM1SPM", "DSLS", "TSSR", "TSSR2"),
  "29"=c("(Intercept)", "TM1SPT", "TM1SPM", "DSLS", "DSLS2", "TSSR", "TSSR2"))

#2. Create design lookup table that describes duration method for each protocol----
#filter out duration methods that aren't appropriate for removal modelling (only have 1 time bin)
durdesign <- visit %>%
  dplyr::select(durationMethod) %>%
  unique() %>%
  dplyr::filter(!durationMethod %in% c("UNKNOWN", "0-10min", "0-20min", "0-5min", "0-3min", "0-2min"),
                !is.na(durationMethod)) %>%
  mutate(dm = str_sub(durationMethod, -100, -4)) %>%
  separate(dm, into=c("t00", "t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), remove=TRUE, sep="-") %>%
  dplyr::select(-t00) %>%
  mutate_at(c("t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), ~as.numeric(.))

durdesign.long <- durdesign %>% 
  pivot_longer(t01:t10, values_to="end", names_to="position")

distdesign <- visit %>%
  dplyr::select(distanceMethod) %>%
  unique() %>%
  dplyr::filter(!distanceMethod %in% c("UNKNOWN", "0m-INF", "0m-100m", "0m-400m", "0m-80m", "0m-50m", "0m-INF-ARU")) %>%
  separate(distanceMethod, into=c("d00", "d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), remove=FALSE, sep="-") %>%
  dplyr::select(-d00) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~str_sub(., -100, -2)) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~ifelse(.=="IN", Inf, .)) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10",  "d11", "d12", "d13"), ~as.numeric(.)/100)

distdesign.long <- distdesign %>% 
  pivot_longer(d01:d13, values_to="far", names_to="position")

#3. Get list of species to process----
set.seed(1234)
spp.avail <- data.frame(phi = seq(round(min(est34$sra4), 1), round(max(est34$sra4)), 0.1)) %>% left_join(est34 %>% mutate(phi = round(sra4, 1)), multiple="all") %>% 
  dplyr::filter(sra.n4 > 5000) %>% 
  group_by(phi) %>% 
  sample_n(5, replace=TRUE) %>% 
  ungroup() %>% 
  unique()

#4. Set up loop for species----
avail <- list()
for(i in 1:nrow(spp.avail)){
  
  #6. Filter abundance data for species---
  # filter out observations with unknown duration method or interval
  # filter out observations with "none" tag method
  # filter out observations where interval does not match method
  # filter to observations with covariates
  # filter human PC to observations in less than 100 m bin
  bird.i <- bird %>%
    dplyr::filter(species==spp.avail$species[i],
                  durationMethod %in% durdesign$durationMethod,
                  !durationInterval %in% c("UNKNOWN", "before or after/incidental"),
                  !is.na(TM)) %>% 
    separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE, extra="drop", fill="right") %>%
    mutate(start = as.numeric(start),
           end = as.numeric(str_sub(end, -100, -4))) %>% 
    separate(distanceBand, into=c("close", "far"), sep="-", remove=FALSE) %>%
    mutate(far = as.numeric(str_sub(far, -100, -2))/100,
           far = ifelse(is.na(far), Inf, far))  %>%
    left_join(durdesign.long, by=c("durationMethod", "end")) %>% 
    dplyr::filter(!is.na(position),
                  (far < 1 | TM!="PC")) %>% 
    group_by(id, durationMethod, durationInterval, position) %>%
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
    x <- visit %>%
      dplyr::filter(id %in% unique(bird.i$id),
                    !is.na(TSSR),
                    !is.na(DSLS),
                    !is.na(JDAY),
                    !is.na(TM),
                    TM!="ARU-None",
                    durationMethod %in% durdesign$durationMethod) %>%
      arrange(id) %>%
      dplyr::select(id, durationMethod, TM, year, JDAY, TSSR, DSLS)
    
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
      dplyr::select(id, position, abundance) %>%
      arrange(position)  %>% 
      rbind(data.frame(id="dummy", position=colnames(d), abundance=NA)) %>% 
      pivot_wider(id_cols=id, names_from=position, values_from=abundance, values_fill=0) %>%
      dplyr::filter(id!="dummy") %>%
      arrange(id) %>%
      dplyr::select(-id) %>%
      as.matrix()
    if(nrow(y) > 0){
      
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
    
  }
  
  print(paste0("Finished modelling species ", spp.avail$English_Name[i], ": ", i, " of ", nrow(spp.avail), " species"))
  
}

names(avail) <- spp.avail$species

save(avail, file=file.path(root, "Results/availability_PhiTauBias.Rdata")) 

#B. PERCEPTIBILITY####

#1. Create list of models----
mods <- list(
  ~ 1,
  ~ TREE,
  ~ LCC2,
  ~ LCC4,
  ~ LCC2 + TREE,
  ~ LCC4 + TREE)
names(mods) <- 0:5
modnames <- list(
  "0"="(Intercept)",
  "1"=c("(Intercept)", "TREE"),
  "2"=c("(Intercept)", "LCC2OpenWet"),
  "3"=c("(Intercept)", "LCC4Conif", "LCC4Open", "LCC4Wet"),
  "4"=c("(Intercept)", "LCC2OpenWet", "TREE"),
  "5"=c("(Intercept)", "LCC4Conif", "LCC4Open", "LCC4Wet", "TREE"))

#4. Get list of species to process----
set.seed(1234)
spp.percep <- data.frame(tau = seq(round(min(est34$edr4), 1), round(max(est34$edr4)), 0.1)) %>% left_join(est34 %>% mutate(tau = round(edr4, 1))) %>% 
  dplyr::filter(edr.n4 > 5000) %>% 
  group_by(tau) %>% 
  sample_n(5, replace=TRUE) %>% 
  ungroup() %>% 
  unique()

#5. Set up loop for species----
percep <- list()
for(i in 1:nrow(spp.percep)){
  
  #6. Filter abundance data for species---
  # filter out observations with unknown duration method or interval
  # filter to observations with covariates
  bird.i <- bird %>%
    dplyr::filter(species==spp.percep$species[i],
                  distanceMethod %in% distdesign$distanceMethod,
                  distanceBand!="UNKNOWN") %>% 
    separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE, extra="drop", fill="right") %>%
    mutate(start = as.numeric(start),
           end = as.numeric(str_sub(end, -100, -4))) %>%
    separate(distanceBand, into=c("close", "far"), sep="-", remove=FALSE) %>%
    mutate(far = as.numeric(str_sub(far, -100, -2))/100,
           far = ifelse(is.na(far), Inf, far)) %>%
    left_join(distdesign.long, by=c("distanceMethod", "far"), multiple="all") %>% 
    left_join(durdesign.long  %>% rename(durposition = position), by=c("durationMethod", "end"), multiple="all") %>% 
    dplyr::filter(!is.na(position),
                  (durposition=="t01")) %>% 
    group_by(id, durationMethod, durationInterval, position) %>%
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
    #Remove nas
    x <- visit %>%
      dplyr::filter(id %in% unique(bird.i$id),
                    !is.na(TREE),
                    !is.na(LCC2),
                    !is.na(LCC4),
                    LCC2!="",
                    LCC4!="",
                    buffer < 1000) %>% 
      arrange(id) %>%
      dplyr::select(id, distanceMethod, TREE, LCC2, LCC4)
    
    #9. Create design matrix----
    d <- x %>%
      dplyr::select(distanceMethod) %>%
      left_join(distdesign, by="distanceMethod") %>%
      dplyr::select(-distanceMethod) %>%
      as.matrix()
    
    #10. Format abundance matrix----
    #add dummy variables for each of the columns to make sure the matrix is the right width
    y <- bird.i %>%
      dplyr::filter(id %in% x$id) %>% 
      dplyr::select(id, position, abundance) %>%
      arrange(position) %>%
      rbind(data.frame(id="dummy", position=colnames(d), abundance=NA)) %>%
      pivot_wider(id_cols=id, names_from=position, values_from=abundance, values_fill=0) %>%
      dplyr::filter(id!="dummy") %>%
      arrange(id) %>%
      dplyr::select(-id) %>%
      as.matrix()
    
    if(nrow(y) > 0){
      
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
        mod <- try(cmulti(f, x, type="dis"))
        if (!inherits(mod, "try-error")) {
          dista <- mod[c("coefficients","vcov","nobs","loglik")]
          dista$p <- length(coef(mod))
          dista$names <- modnames[[j]]
        } else {
          dista <- mod
        }
        mod.list[[names(modnames)[j]]] <- dista
      }
      
      #13. Save model results to species list---
      percep[[i]] <- mod.list
      
    }
    
  }
  
  print(paste0("Finished modelling species ", spp.percep$English_Name[i], ": ", i, " of ", nrow(spp.percep), " species"))
  
}

names(percep) <- spp.percep$species

save(percep, file=file.path(root, "Results/perceptibility_PhiTauBias.Rdata")) 

#C. COMPARISON####

load(file.path(root, "Results/availability_PhiTauBias.Rdata"))
load(file.path(root, "Results/perceptibility_PhiTauBias.Rdata"))

#1. Wrangle availability results----
est.avail <- est34 %>% 
  dplyr::select(sra4, sra4.pc, sra4.spt, sra4.spm, edr4, species) %>% 
  mutate(method="all") %>% 
  dplyr::filter(species %in% spp.avail$species)

for(i in 1:nrow(spp.avail)){
  
  sra4 <- exp(avail[[i]]$"0"$coefficients[1])
  sra4.pc <- exp(avail[[i]]$"15"$coefficients[[1]])
  sra4.spt <- exp(sum(avail[[i]]$"15"$coefficients[c(1,2)]))
  sra4.spm <- exp(sum(avail[[i]]$"15"$coefficients[c(1,3)]))
  est.avail <- rbind(est.avail,
                data.frame(sra4=sra4, sra4.pc=sra4.pc, sra4.spt=sra4.spt, sra4.spm=sra4.spm, edr4=est.avail$edr4[est.avail$species==spp.avail$species[i]],species=spp.avail$species[i], method = "first"))
}

est.avail.wide <- est.avail %>% 
  pivot_wider(id_cols=species, names_from=method, values_from=sra4:edr4) %>% 
  data.frame()

est.avail.wide$density_all <- 10/((pi*est.avail.wide$edr4_all^2)*(1-exp(-5*est.avail.wide$sra4.pc_all)))
est.avail.wide$density_first <- 10/((pi*est.avail.wide$edr4_first^2)*(1-exp(-5*est.avail.wide$sra4.pc_first)))
                                  

#2. Wrangle perceptibility results----
est.percep <- est34 %>% 
  dplyr::select(edr4, sra4.pc, species) %>% 
  mutate(method="all") %>% 
  dplyr::filter(species %in% names(percep))

for(i in 1:length(percep)){
  
  edr4 <- exp(percep[[i]]$"0"$coefficients[1])
  est.percep <- rbind(est.percep,
                     data.frame(edr4=edr4, sra4.pc=est.percep$sra4.pc[est.percep$species==names(percep)[i]], species=names(percep)[i], method = "first"))
}

est.percep.wide <- est.percep %>% 
  pivot_wider(id_cols=species, names_from=method, values_from=edr4:sra4.pc) %>% 
  data.frame()

est.percep.wide$density_all <- 10/((pi*est.percep.wide$edr4_all^2)*(1-exp(-5*est.percep.wide$sra4.pc_all)))
est.percep.wide$density_first <- 10/((pi*est.percep.wide$edr4_first^2)*(1-exp(-5*est.percep.wide$sra4.pc_first)))

#3. Calculate density----

#3. Plot----
plot.avial1 <- ggplot(est.avail.wide) +
  geom_text(aes(x=sra4_all, y=sra4_first, label=species)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("all data") +
  ylab("distance < 100 m") +
  ggtitle("Null phi estimate")
plot.avial1

plot.avail2 <- ggplot(est.avail.wide) +
  geom_text(aes(x=sra4.pc_all, y=sra4.pc_first, label=species)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("all data") +
  ylab("distance < 100 m") +
  ggtitle("Human point count phi estimate")
plot.avail2

plot.avail.density <- ggplot(est.avail.wide) +
  geom_text(aes(x=density_all, y=density_first, label=species)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("all data") +
  ylab("distance < 100 m") +
  ggtitle("Human point count phi estimate")
plot.avail.density

plot.percep <- ggplot(est.percep.wide) +
  geom_text(aes(x=all, y=first, label=species)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("all data") +
  ylab("first time bin") +
  ggtitle("Null tau estimate")
plot.percep

plot.percep.density <- ggplot(est.percep.wide) +
  geom_text(aes(x=density_all, y=density_first, label=species)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("all data") +
  ylab("distance < 100 m") +
  ggtitle("Human point count phi estimate")
plot.percep.density

ggsave(grid.arrange(plot.avial1, plot.avail2, plot.percep, nrow=3),
       filename=file.path(root, "Figures", "Tau&PhiTruncation.jpeg"), 
       height=18, width=6, device="jpeg")
