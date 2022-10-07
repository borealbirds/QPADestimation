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
modnames <- c("(Intercept)",
              "(Intercept) + jday",
              "(Intercept) + tssr",
              "(Intercept) + jday + jday2",
              "(Intercept) + tssr + tssr2",
              "(Intercept) + jday + tssr",
              "(Intercept) + jday + jday2 + tssr",
              "(Intercept) + jday + tssr + tssr2",
              "(Intercept) + jday + jday2 + tssr + tssr2",
              "(Intercept) + tsg",
              "(Intercept) + tsg + tsg2",
              "(Intercept) + tsg + tssr",
              "(Intercept) + tsg + tsg2 + tssr",
              "(Intercept) + tsg + tssr + tssr2",
              "(Intercept) + tsg + tsg2 + tssr + tssr2")

#2. Load data----
load("data/cleaned_data_2022-10-06.Rdata")

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
                      durationMethod %in% design$durationMethod,
                      durationInterval!="UNKNOWN") %>%
        group_by(id, durationMethod, durationInterval) %>%
        summarize(abundance = sum(abundance)) %>%
        ungroup()

    #only model if there is data
    if(nrow(bird.i) > 0){

        #7. Replace abundance outliers (defined as 99% quantile across species) with 99% quantile of abundance----
        if(max(bird.i$abundance) > quantile(bird$abundance, 0.99)){
            bird.i <- bird.i %>%
                mutate(abundance = ifelse(abundance > quantile(abundance, 0.99), round(quantile(abundance, 0.99)), abundance))
        }

        #8. Filter visit covariates----
        #Remove records with nas in covariates
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
            separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE, extra="drop", fill="right") %>%
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
        #Save a bunch of metadata like sample size and aic value
        mod.list <- list()
        for (j in 1:length(mods)) {
            f <- as.formula(paste0("y | d ", paste(as.character(mods[[j]]), collapse=" ")))
            mod <- try(cmulti(f, x, type="rem"))
            if (!inherits(mod, "try-error")) {
                rmvl <- data.frame(t(data.frame(mod["coefficients"]))) %>%
                    mutate(nobs=mod["nobs"]$nobs,
                           loglik = mod["loglik"]$loglik,
                           df = length(coef(mod)),
                           aic = AIC(mod),
                           aicc = aic + (2*df^2+2*df) / (nrow(y)-df-1),
                           model = modnames[j],
                           species = spp$speciesCode[i])
            } else {
                rmvl <- data.frame(nobs="try-error",
                                   species=spp$speciesCode[i])
            }
            mod.list[[j]] <- rmvl
        }

        #13. Save model results to species list---
        species.list[[i]] <- rbindlist(mod.list, fill=TRUE)

    }

    print(paste0("Finished modelling species ", spp$English_Name[i], ": ", i, " of ", nrow(spp), " species"))

}

species.out <- rbindlist(species.list, fill=TRUE)

#14. Identify species that failed----
species.fail <- species.out %>%
    dplyr::filter(nobs=="try-error") %>%
    dplyr::select(species) %>%
    unique()

species.use <- species.out %>%
    dplyr::filter(!species %in% species.fail$species)

#15. Save
save(species.out, species.use, file="results/availability_results_2022-10-06.Rdata")
