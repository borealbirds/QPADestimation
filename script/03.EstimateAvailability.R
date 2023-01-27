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

#TO DO: TAKE OUT RECORDS WITH DUMMY ENTRY FOR TIME AND DATE

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

#2. Load data----
load(file.path(root, "Data", "qpadv4_clean.Rdata"))

 #Set factor levels
visit$TM <- factor(visit$TM, levels=c("PC", "ARU-1SPT", "ARU-1SPM"))
#visit$sensor <- factor(visit$sensor, levels=c("PC", "ARU"))

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
    # filter out observations with "none" tag method
    # filter out observations where interval does not match method
    # filter to observations with covariates
    bird.i <- bird %>%
        dplyr::filter(species==spp$species[i],
                      durationMethod %in% durdesign$durationMethod,
                      !durationInterval %in% c("UNKNOWN", "before or after/incidental"),
                      !is.na(TM)) %>% 
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
    
    print(paste0("Finished modelling species ", spp$English_Name[i], ": ", i, " of ", nrow(spp), " species"))
    
}

names(avail) <- spp$species

#14. Save out results----
save(durdesign, avail, file=file.path(root, "Results/availability.Rdata")) 
