# ---
# title: "QPAD estimation - estimate availability"
# author: "Elly Knight"
# created: "September 22, 2022"
# adapted from: "QPAD version 3 documentation" by Peter Solymos https://github.com/borealbirds/bamanalytics/blob/master/projects/qpad_v3/QPAD-v3-report.Rmd
# ---

library(tidyverse) #basic data wrangling
library(detect) #removal models
library(data.table) #collapse list to dataframe

options(dplyr.summarise.inform = FALSE, scipen=9999)

root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/"

#1. Create list of models----
#jday = day of year as a decimal between 0 and 1
#tssr = time since sunrise as a decimal between 0 and 1
#tsg = days since start of seedgrowth from seedgrow layers

#1a. Lists for species that have data for both sensors----
mods.2 <- list(
    ~ sensor,
    ~ sensor + jday,
    ~ sensor + tssr,
    ~ sensor + poly(jday, 2),
    ~ sensor + poly(tssr, 2),
    ~ sensor + jday + tssr,
    ~ sensor + poly(jday, 2) + tssr,
    ~ sensor + jday + poly(tssr, 2),
    ~ sensor + poly(jday, 2) + poly(tssr, 2),
    ~ sensor + tsg,
    ~ sensor + poly(tsg, 2),
    ~ sensor + tsg + tssr,
    ~ sensor + poly(tsg, 2) + tssr,
    ~ sensor + tsg + poly(tssr, 2),
    ~ sensor + poly(tsg, 2) + poly(tssr, 2))
names(mods.2) <- 0:14
modnames.2 <- list(
    "0"=c("(Intercept)", "sensor"),
    "1"=c("(Intercept)", "sensor", "jday"),
    "2"=c("(Intercept)", "sensor", "tssr"),
    "3"=c("(Intercept)", "sensor", "jday", "jday2"),
    "4"=c("(Intercept)", "sensor", "tssr", "tssr2"),
    "5"=c("(Intercept)", "sensor", "jday", "tssr"),
    "6"=c("(Intercept)", "sensor", "jday", "jday2", "tssr"),
    "7"=c("(Intercept)", "sensor", "jday", "tssr", "tssr2"),
    "8"=c("(Intercept)", "sensor", "jday", "jday2", "tssr", "tssr2"),
    "9"=c("(Intercept)", "sensor", "tsg"),
    "10"=c("(Intercept)", "sensor", "tsg", "tsg2"),
    "11"=c("(Intercept)", "sensor", "tsg", "tssr"),
    "12"=c("(Intercept)", "sensor", "tsg", "tsg2", "tssr"),
    "13"=c("(Intercept)", "sensor", "tsg", "tssr", "tssr2"),
    "14"=c("(Intercept)", "sensor", "tsg", "tsg2", "tssr", "tssr2"))

#1b. Lists for species that only have data for one sensor----
mods.1 <- list(
  ~ 1,
  ~ jday,
  ~ tssr,
  ~ poly(jday, 2),
  ~ poly(tssr, 2),
  ~ jday + tssr,
  ~ poly(jday, 2) + tssr,
  ~ jday + poly(tssr, 2),
  ~ poly(jday, 2) + poly(tssr, 2),
  ~ tsg,
  ~ poly(tsg, 2),
  ~ tsg + tssr,
  ~ poly(tsg, 2) + tssr,
  ~ tsg + poly(tssr, 2),
  ~ poly(tsg, 2) + poly(tssr, 2))
names(mods.1) <- 0:14
modnames.1 <- list(
  "0"=c("(Intercept)"),
  "1"=c("(Intercept)", "jday"),
  "2"=c("(Intercept)", "tssr"),
  "3"=c("(Intercept)", "jday", "jday2"),
  "4"=c("(Intercept)", "tssr", "tssr2"),
  "5"=c("(Intercept)", "jday", "tssr"),
  "6"=c("(Intercept)", "jday", "jday2", "tssr"),
  "7"=c("(Intercept)", "jday", "tssr", "tssr2"),
  "8"=c("(Intercept)", "jday", "jday2", "tssr", "tssr2"),
  "9"=c("(Intercept)", "tsg"),
  "10"=c("(Intercept)", "tsg", "tsg2"),
  "11"=c("(Intercept)", "tsg", "tssr"),
  "12"=c("(Intercept)", "tsg", "tsg2", "tssr"),
  "13"=c("(Intercept)", "tsg", "tssr", "tssr2"),
  "14"=c("(Intercept)", "tsg", "tsg2", "tssr", "tssr2"))


#2. Load data----
load(file.path(root, "Data/qpadv4_clean.Rdata"))

#set intercept as point count sensor
visit$sensor <- factor(visit$sensor, levels=c("PC", "ARU"))

#3. Create design lookup table that describes duration method for each protocol----
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

#4. Get list of species to process----
spp <- species %>%
    dplyr::filter(Singing_birds==TRUE) %>%
    left_join(bird %>%
                  dplyr::select(species) %>%
                  unique()) %>%
    arrange(species)

#5. Set up loop for species----
avail <- list()
for(i in 1:nrow(spp)){

    #6. Filter abundance data for species---
    # filter out observations with unknown duration method or interval
    # filter out observations where interval does not match method
    # filter to observations with covariates
    bird.i <- bird %>%
        dplyr::filter(species==spp$species[i],
                      durationMethod %in% durdesign$durationMethod,
                      !durationInterval %in% c("UNKNOWN", "before or after/incidental")) %>% 
      separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE, extra="drop", fill="right") %>%
      mutate(start = as.numeric(start),
             end = as.numeric(str_sub(end, -100, -4))) %>%
      left_join(durdesign.long, by=c("durationMethod", "end")) %>% 
      dplyr::filter(!is.na(position)) %>% 
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
                          !is.na(tssr),
                          !is.na(tsg),
                          !is.na(jday),
                          durationMethod %in% durdesign$durationMethod) %>%
            arrange(id) %>%
            dplyr::select(id, durationMethod, sensor, year, jday, tssr, tsg)

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

        #11. Change zeros to NAs in the abundance matrix to match the design matrix----
        for (j in 1:nrow(y)){
            indices <- which(is.na(d[j,]))
            y[j, indices] <- NA
        }
        
        #12. Determine model set based on # of sensors in data and comparison with sensor only model----
        if(length(unique(x$sensor))> 1){
          
          mod.null <- try(cmulti(y | d ~ 1, x, type="rem"))
          mod.sensor <- try(cmulti(y | d ~ sensor, x, type="rem"))
          
          if(class(mod.null)!="try-error" & class(mod.null)!="try-error"){
            
            if(AIC(mod.sensor) < AIC(mod.null)){
              mods <- mods.2
              modnames <- modnames.2
            }
            
            if(AIC(mod.sensor) >= AIC(mod.null)){
              mods <- mods.1
              modnames <- modnames.1
            }
            
            #species for which the PC model gives an unrealistically low estimate
            if(exp(mod.sensor$coefficients[[1]] < 0.01)){
              mods <- mods.1
              modnames <- modnames.1
            }
            
          }
          
          if(class(mod.null)=="try-error" | class(mod.sensor)=="try-error"){
            mods <- mods.1
            modnames <- modnames.1
          }
          
        }
        
        if(length(unique(x$sensor))==1){
          mods <- mods.1
          modnames <- modnames.1
        }
        
        #13. Fit models----
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
        
        #14. Save model results to species list---
        avail[[i]] <- mod.list
        
    }
    
    print(paste0("Finished modelling species ", spp$English_Name[i], ": ", i, " of ", nrow(spp), " species"))
    
}

names(avail) <- spp$species

#15. Save out results----
save(durdesign, avail, file=file.path(root, "Results/availability.Rdata")) 
