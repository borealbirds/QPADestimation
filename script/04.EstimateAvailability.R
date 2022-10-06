# ---
# title: "QPAD estimation - estimate availability"
# author: "Elly Knight"
# created: "September 22, 2022"
# adapted from: "QPAD estimation by Peter Solymos"
# ---

library(tidyverse) #basic data wrangling
library(detect) #removal models

options(dplyr.summarise.inform = FALSE, scipen=9999)

#1. Create list of models----
#jday = day of year as a decimal between 0 and 1
#tssr = time since sunrise as a decimal between 0 and 1
#tsg = days since start of seedgrowth from seedgrow layers

mods <- list(
    ~ 1,
    ~ jday,
    ~ tssr,
    ~ jday + jday2,
    ~ tssr + tssr2,
    ~ jday + tssr,
    ~ jday + jday2 + tssr,
    ~ jday + tssr + tssr2,
    ~ jday + jday2 + tssr + tssr2,
    ~ tsg,
    ~ tsg + tsg2,
    ~ tsg + tssr,
    ~ tsg + tsg2 + tssr,
    ~ tsg + tssr + tssr2,
    ~ tsg + tsg2 + tssr + tssr2)
names(mods) <- 0:14
modnames <- list(
    "0"="(Intercept)",
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
load("data/cleaned_data_2022-07-24.Rdata")

#3. Create design lookup table that describes duration method for each protocol----
#filter out duration methods that aren't appropriate for removal modelling (only have 1 time bin)
design <- visit %>%
    dplyr::select(durationMethod) %>%
    unique() %>%
    dplyr::filter(!durationMethod %in% c("UNKNOWN", "0-10min", "0-20min", "0-5min", "0-3min", "0-2min")) %>%
    mutate(dm = str_sub(durationMethod, -100, -4)) %>%
    separate(dm, into=c("t00", "t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), remove=TRUE, sep="-") %>%
    dplyr::select(-t00) %>%
    mutate_at(c("t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), ~as.numeric(.))

#4. Get list of species to process----
spp <- species %>%
    dplyr::filter(Singing_birds==TRUE) %>%
    left_join(bird %>%
                  dplyr::select(speciesCode) %>%
                  unique()) %>%
    arrange(speciesCode)

#5. Set up loop for species----
species.list <- list()
for(i in 1:nrow(spp)){

    #6. Filter abundance data for species---
    # filter out observations with unknown duration method or interval
    # filter to observations with covariates
    bird.i <- bird %>%
        dplyr::filter(speciesCode==spp$speciesCode[i],
                      durationMethod %in% design$durationMethod) %>%
        group_by(id, durationMethod, durationInterval) %>%
        summarize(abundance = sum(abundance)) %>%
        ungroup()

    #7. Replace abundance outliers (defined as 99% quantile across species) with 99% quantile of abundance----
    if(max(bird.i$abundance) > quantile(bird$abundance, 0.99)){
        bird.i <- bird.i %>%
            mutate(abundance = ifelse(abundance > quantile(abundance, 0.99), quantile(abundance, 0.99), abundance))
    }

    #8. Filter visit covariates----
    #Remove nas
    x <- visit %>%
        dplyr::filter(id %in% unique(bird.i$id),
                      !is.na(tssr),
                      !is.na(tsg),
                      !is.na(jday)) %>%
        arrange(id) %>%
        dplyr::select(id, durationMethod, jday, jday2, tssr, tssr2, tsg, tsg2)

    #9. Create design matrix----
    d <- x %>%
        dplyr::select(durationMethod) %>%
        left_join(design, by="durationMethod") %>%
        dplyr::select(-durationMethod) %>%
        as.matrix()

    #10. Format abundance matrix----
    #add dummy variables for each of the columns to make sure the matrix is the right width
    y <- bird.i %>%
        dplyr::filter(id %in% x$id) %>%
        separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE) %>%
        mutate(start = as.numeric(start),
               end = as.numeric(str_sub(end, -100, -4))) %>%
        left_join(design %>%
                      pivot_longer(t01:t10, values_to="end", names_to="position"),
                  by=c("durationMethod", "end")) %>%
        dplyr::select(id, position, abundance) %>%
        arrange(position) %>%
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
    mod.list <- list()
    for (j in 1:length(mods)) {
        f <- as.formula(paste0("y | d ", paste(as.character(mods[[j]]), collapse=" ")))
        mod <- try(cmulti(f, x, type="rem"))
        if (!inherits(mod, "try-error")) {
            rval <- mod[c("coefficients","vcov","nobs","loglik")]
            rval$p <- length(coef(mod))
            rval$names <- modnames[[j]]
            rval$aic <- AIC(mod)
            rval$aicc <- rval$aic + (2*rval$p^2+2*rval$p) / (nrow(y)-rval$p-1)
        } else {
            rval <- mod
        }
        mod.list[[names(modnames)[[j]]]] <- rval
    }

    #13. Save model results to species list---
    species.list[[spp$speciesCode[[i]]]] <- mod.list

    print(paste0("Finished modelling species ", spp$English_Name[i], ": ", i, " of ", nrow(spp), " species"))

}

#14. Remove species that failed----

#15. Save
save(species.list, file="results/availability_results_2022-07-24.Rdata")
