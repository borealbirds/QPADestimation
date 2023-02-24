library(tidyverse) #basic data wrangling
library(detect) #traditional models
library(data.table) #collapse list to dataframe
library(CmultiJoint.dev) #joint models

options(dplyr.summarise.inform = FALSE, scipen=9999)

root <- "G:/Shared drives/BAM_RshProjs/PopnStatus/QPAD"

#1. Load data----
load(file.path(root, "Data", "qpadv4_clean.Rdata"))

#Set factor levels
visit$TM <- factor(visit$TM, levels=c("PC", "ARU-1SPT", "ARU-1SPM"))

#fix ARU distance band
bird$distanceBand <- ifelse(bird$TM!="PC", "0m-INF", bird$distanceBand)

#fix ARU distance design
visit$distanceMethod <- ifelse(visit$distanceMethod=="0m-INF-ARU", "0m-INF", visit$distanceMethod)

#2. Filter bird data----
#no NONE tag method
#no na TM
bird.use <- dplyr::filter(bird,
                      TM!="ARU-NONE",
                      !is.na(TM))

#3. Create design lookup tables----
durdesign <- visit %>%
  dplyr::select(durationMethod) %>%
  unique() %>%
  dplyr::filter(!durationMethod %in% c("UNKNOWN"),
                !is.na(durationMethod)) %>%
  mutate(dm = str_sub(durationMethod, -100, -4)) %>%
  separate(dm, into=c("t00", "t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), remove=TRUE, sep="-") %>%
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
  separate(distanceMethod, into=c("d00", "d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), remove=FALSE, sep="-") %>%
  dplyr::select(-d00) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~str_sub(., -100, -2)) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~ifelse(.=="IN", Inf, .)) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10",  "d11", "d12", "d13"), ~as.numeric(.)/100)

distdesign.long <- distdesign %>% 
  pivot_longer(d01:d13, values_to="far", names_to="distposition") %>% 
  dplyr::filter(!is.na(far))

#4. Get list of species to process----
spp <- species %>%
      dplyr::filter(Singing_birds==TRUE) %>%
  left_join(bird %>%
              dplyr::select(species) %>%
              unique()) %>%
  arrange(species)

#5. Set up loop for species----
for(i in 1:nrow(spp)){
  
  start <- Sys.time()
  
  #6. Filter abundance data for species---
  # filter out observations with unknown duration method or interval
  # filter out observations with "none" tag method
  # filter out observations where interval does not match method
  bird.i <- bird.use %>%
    dplyr::filter(species==spp$species[i],
                  durationMethod %in% durdesign$durationMethod,
                  distanceMethod %in% distdesign$distanceMethod,
                  !durationInterval %in% c("UNKNOWN", "before or after/incidental"),
                  !is.na(TM)) %>% 
    separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE, extra="drop", fill="right") %>%
    mutate(start = as.numeric(start),
           end = as.numeric(str_sub(end, -100, -4))) %>%
    left_join(durdesign.long, by=c("durationMethod", "end")) %>% 
    separate(distanceBand, into=c("near", "far"), sep="-", remove=FALSE, extra="drop", fill="right") %>% 
    mutate(near = as.numeric(str_sub(near, -100, -2)),
           far = ifelse(far=="INF", Inf, as.numeric(str_sub(far, -100, -2)))) %>% 
    left_join(distdesign.long, by=c("distanceMethod", "far")) %>% 
    dplyr::filter(!is.na(durposition),
                  !is.na(distposition)) %>% 
    group_by(id, durationMethod, distanceMethod, durationInterval, distanceBand, end, far, durposition, distposition) %>%
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
    
    rarray <- x %>%
      dplyr::select(distanceMethod) %>%
      left_join(distdesign, by="distanceMethod") %>%
      dplyr::select(-distanceMethod) %>%
      as.matrix()
    
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
    
    #11. Create covariate matrix----
    X1.null <- model.matrix(~rep(1, nrow(x)))
    colnames(X1.null) <- c("tau_int")
    
    # Design matrix for phi (scaled covariate for better model convergence)
    X2.null <- model.matrix(~rep(1, nrow(x)))
    colnames(X2.null) <- c("phi_int")
    X2.TM <- model.matrix(~x$TM)
    colnames(X2.TM) <- c("phi_int","phi_SPT", "phi_SPM")
    
    #12. Fit models----
    #Save a bunch of metadata like sample size and aic value
    mod.list <- list()
    fit.null <- try(cmulti_fit_joint(Yarray,
                            rarray,
                            tarray,
                            X1 = X1.null,
                            X2 = X2.null))
    
    if (!inherits(fit.null, "try-error")) {
      out <- fit.null[c("coefficients","vcov","loglik")]
      out$nobs <- nrow(x)
      out$p <- length(coef(fit.null))
      out$names <- "TM"
    } else {
      out <- fit.null
    }
    mod.list[[1]] <- out
    
    fit.TM <- try(cmulti_fit_joint(Yarray,
                                 rarray,
                                 tarray,
                                 X1 = X1.null,
                                 X2 = X2.TM))
    
    if (!inherits(fit.TM, "try-error")) {
      out <- fit.TM[c("coefficients","vcov","loglik")]
      out$nobs <- nrow(x)
      out$p <- length(coef(fit.TM))
      out$names <- "TM"
    } else {
      out <- fit.TM
    }
    mod.list[[2]] <- out

    #13. Save model results to species list---
    joint[[i]] <- mod.list
  
  }
  
  #14. Save results to local----
  names(joint) <- spp$species[1:i]
  
  save(joint, file=file.path(root, "Results/jointestimates.Rdata")) 
  
  end <- Sys.time()
  
  print(paste0("Finished modelling species ", spp$English_Name[i], ": ", i, " of ", nrow(spp), " species in ", end-start, " minutes"))
}
