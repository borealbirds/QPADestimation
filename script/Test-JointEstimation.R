library(tidyverse) #basic data wrangling
library(detect) #traditional models
library(data.table) #collapse list to dataframe
library(CmultiJoint.dev) #joint models

options(dplyr.summarise.inform = FALSE, scipen=9999)

root <- "G:/Shared drives/BAM_RshProjs/PopnStatus/QPAD"

#A. ESTIMATES####

#1. Load data----
load(file.path(root, "Data", "qpadv4_clean.Rdata"))

#Set factor levels - turn NONE to NA
visit$TM <- factor(visit$TM, levels=c("PC", "ARU-1SPT", "ARU-1SPM"))
bird$TM <- factor(bird$TM, levels=c("PC", "ARU-1SPT", "ARU-1SPM"))

#fix ARU distance band
bird$distanceBand <- ifelse(bird$TM!="PC", "0m-INF", bird$distanceBand)

#fix ARU distance design
visit$distanceMethod <- ifelse(visit$distanceMethod=="0m-INF-ARU", "0m-INF", visit$distanceMethod)
visit$distanceMethod <- ifelse(visit$TM %in% c("ARU-1SPT", "ARU-1SPM"), "0m-INF", visit$distanceMethod)
bird$distanceMethod <- ifelse(bird$distanceMethod=="0m-INF-ARU", "0m-INF", bird$distanceMethod)
bird$distanceMethod <- ifelse(bird$TM %in% c("ARU-1SPT", "ARU-1SPM"), "0m-INF", bird$distanceMethod)

#collapse 0-3-5-10-10min+ and 0-3-5-10min durationMethods (no 10+ detections in dataset)
visit$durationMethod <- ifelse(visit$durationMethod=="0-3-5-10-10min+", "0-3-5-10min", visit$durationMethod)
bird$durationMethod <- ifelse(bird$durationMethod=="0-3-5-10-10min+", "0-3-5-10min", bird$durationMethod)

#2. Filter bird data----
#no NONE tag method
#no na Tag method
#no unknown methods, bands, intervals
bird.use <- dplyr::filter(bird,
                      TM!="ARU-None",
                      !is.na(TM),
                      durationMethod!="UNKNOWN",
                      distanceMethod!="UNKNOWN",
                      !durationInterval %in% c("UNKNOWN", "before or after/incidental"),
                      distanceBand!="UNKNOWN")

#3. Create design lookup tables----
durdesign <- visit %>%
  dplyr::select(durationMethod) %>%
  unique() %>%
  dplyr::filter(!durationMethod %in% c("UNKNOWN"),
                !is.na(durationMethod)) %>%
  mutate(dm = str_sub(durationMethod, -100, -4)) %>%
  separate(dm, into=c("t00", "t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), remove=TRUE, sep="-", extra="drop", fill="right") %>%
  dplyr::select(-t00) %>% 
  mutate_at(c("t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), ~as.numeric(.))

durdesign.long <- durdesign %>% 
  pivot_longer(t01:t10, values_to="end", names_to="durposition") %>% 
  dplyr::filter(!is.na(end))

distdesign <- visit %>%
  dplyr::select(distanceMethod) %>%
  unique() %>%
  dplyr::filter(!distanceMethod %in% c("UNKNOWN"),
                !is.na(distanceMethod)) %>%
  separate(distanceMethod, into=c("d00", "d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), remove=FALSE, sep="-", extra="drop", fill="right") %>%
  dplyr::select(-d00) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~str_sub(., -100, -2)) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~ifelse(.=="IN", Inf, .)) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10",  "d11", "d12", "d13"), ~as.numeric(.)/100)

distdesign.long <- distdesign %>% 
  pivot_longer(d01:d13, values_to="far", names_to="distposition") %>% 
  dplyr::filter(!is.na(far))

#4. Get list of species to process----
#only species with at least 25 obs per method
specieslist <- species %>% 
  dplyr::filter(Singing_birds==TRUE)

spp <- bird.use %>% 
  dplyr::filter(species %in% specieslist$species) %>% 
  dplyr::select(species, project, location, date, TM, distanceBand, durationInterval) %>% 
  unique() %>% 
  group_by(species, TM) %>% 
  summarize(obs=n()) %>% 
  ungroup() %>% 
  dplyr::filter(obs >= 25) %>% 
  group_by(species) %>% 
  summarize(obs=sum(obs),
            methods = n()) %>% 
  ungroup() %>% 
  dplyr::filter(methods==3)

#5. Set up loop for species----
mod.list <- list()
for(i in 100:nrow(spp)){
  
  start <- Sys.time()
  
  #6. Filter abundance data for species---
  # filter out observations with unknown duration method or interval
  # filter out observations where interval does not match method
  # sum abundance for each bin in each visit
  bird.i <- bird.use %>%
    dplyr::filter(species==spp$species[i],
                  durationMethod %in% durdesign$durationMethod,
                  distanceMethod %in% distdesign$distanceMethod) %>% 
    separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE, extra="drop", fill="right") %>%
    mutate(start = as.numeric(start),
           end = as.numeric(str_sub(end, -100, -4))) %>% 
    left_join(durdesign.long, by=c("durationMethod", "end")) %>% 
    separate(distanceBand, into=c("near", "far"), sep="-", remove=FALSE, extra="drop", fill="right") %>% 
    mutate(near = as.numeric(str_sub(near, -100, -2)),
           far = ifelse(far=="INF", Inf, suppressWarnings(as.numeric(str_sub(far, -100, -2))/100))) %>% 
    left_join(distdesign.long, by=c("distanceMethod", "far")) %>% 
    dplyr::filter(!is.na(durposition),
                  !is.na(distposition)) %>% 
    group_by(id, TM, durationMethod, distanceMethod, durationInterval, distanceBand, end, far, durposition, distposition) %>%
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
      dplyr::filter(id %in% unique(bird.i$id)) %>%
      arrange(id) %>%
      dplyr::select(id, durationMethod, distanceMethod, TM)
    
    #9. Create design matrices----
    tarray <- x %>%
      dplyr::select(durationMethod) %>%
      left_join(durdesign, by="durationMethod") %>%
      dplyr::select(-durationMethod) %>%
      as.matrix()
    dimnames(tarray) <- NULL
    
    rarray <- x %>%
      dplyr::select(distanceMethod) %>%
      left_join(distdesign, by="distanceMethod") %>%
      dplyr::select(-distanceMethod) %>%
      as.matrix()
    dimnames(rarray) <- NULL
    
    #10. Create visit matrix-----
    Yarray <- array(NA,dim=c(nrow(x),ncol(distdesign)-1,ncol(durdesign)-1))
    for(j in 1:nrow(x)){
      
      bird.j <- bird.i %>% 
        dplyr::filter(id==x$id[j])
      bird.j$durint <- cut(bird.j$end,c(0,tarray[j,]))
      bird.j$distint <- cut(bird.j$far, c(0,rarray[j,]))
      
      Y <- table(bird.j[,c("distint", "durint")])
      Yarray[j,
             1:length(levels(bird.j$distint)),
             1:length(levels(bird.j$durint))] <- Y
    }
    
    #11. Create covariate matrices----
    X1.TM <- model.matrix(~x$TM)
    colnames(X1.TM) <- c("tau_int","tau_SPT", "tau_SPM")
    
    X2.TM <- model.matrix(~x$TM)
    colnames(X2.TM) <- c("phi_int","phi_SPT", "phi_SPM")
    
    X1 <- list(NULL, X1.TM, NULL, X1.TM)
    X2 <- list(NULL, NULL, X2.TM, X2.TM)
    
    names <- list(c("tau_int", "phi_int"),
                  c("tau_int", "tau_SPT", "tau_SPM", "phi_int"),
                  c("tau_int", "phi_int", "phi_SPT", "phi_SPM"),
                  c("tau_int", "tau_SPT", "tau_SPM", "phi_int", "phi_SPT", "phi_SPM"))
    
    fit <- list()
    for(j in 1:4){
      fit[[j]] <- try(cmulti_fit_joint(Yarray,
                                       rarray,
                                       tarray,
                                       X1 = X1[[j]],
                                       X2 = X2[[j]]))
    }
    names(fit) <- c("tau-null_phi-null", "tau-TM_phi-null", "tau-null_phi-TM", "tau-TM_phi-TM")
    mod.list[[i]] <- fit
    
    names(mod.list) <- spp$species[1:i]
  
  }
  
  #14. Save results to local----
  
  save(mod.list, file=file.path(root, "Results/jointmodels.Rdata")) 
  
  end <- Sys.time()
  
  print(paste0("Finished modelling species ", spp$species[i], ": ", i, " of ", nrow(spp), " species in ", round(end-start, 2), " minutes"))
}

