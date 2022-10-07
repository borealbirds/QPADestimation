# ---
# title: "QPAD estimation - estimate perceptability"
# author: "Elly Knight"
# created: "October 6, 2022"
# adapted from: "QPAD version 3 documentation" by Peter Solymos https://github.com/borealbirds/bamanalytics/blob/master/projects/qpad_v3/QPAD-v3-report.Rmd
# ---

library(tidyverse) #basic data wrangling
library(detect) #removal models
library(data.table) #collapse list to dataframe

#1. Create list of models----
mods <- list(
    ~ 1,
    ~ tree,
    ~ lcc2,
    ~ lcc4,
    ~ lcc2 + tree,
    ~ lcc4 + tree)
names(mods) <- 0:5
modnames <- list(
    "0"="(Intercept)",
    "1"=c("(Intercept)", "tree"),
    "2"=c("(Intercept)", "lcc2openwet"),
    "3"=c("(Intercept)", "lcc4conif", "lcc4open", "lcc4wet"),
    "4"=c("(Intercept)", "lcc2openwet", "tree"),
    "5"=c("(Intercept)", "lcc4conif", "lcc4open", "lcc4wet", "tree"))

#2. Load data----
load("data/cleaned_data_2022-10-06.Rdata")

#3. Create design lookup table that describes duration method for each protocol----
#filter out distance methods that aren't appropriate for removal modelling (only have 1 bin)
design <- visit %>%
    dplyr::select(distanceMethod) %>%
    unique() %>%
    dplyr::filter(!distanceMethod %in% c("UNKNOWN", "0m-INF", "0m-100m", "0m-400m", "0m-80m", "0m-50m", "0m-INF-ARU")) %>%
    separate(distanceMethod, into=c("d00", "d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), remove=FALSE, sep="-") %>%
    dplyr::select(-d00) %>%
    mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~str_sub(., -100, -2)) %>%
    mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~ifelse(.=="IN", Inf, .)) %>%
    mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10",  "d11", "d12", "d13"), ~as.numeric(.)/100)

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
                      distanceMethod %in% design$distanceMethod,
                      distanceBand!="UNKNOWN") %>%
        group_by(id, distanceMethod, distanceBand) %>%
        summarize(abundance = sum(abundance)) %>%
        ungroup()

    #only model if there is data
    if(nrow(bird.i) > 0){

        #7. Replace abundance outliers (defined as 99% quantile across species) with 99% quantile of abundance----
        if(max(bird.i$abundance) > quantile(bird$abundance, 0.99)){
            bird.i <- bird.i %>%
                mutate(abundance = ifelse(abundance > quantile(abundance, 0.99), quantile(abundance, 0.99), abundance))
        }

        #8. Filter visit covariates----
        #Remove nas
        x <- visit %>%
            dplyr::filter(id %in% unique(bird.i$id),
                          !is.na(tree),
                          !is.na(lcc2),
                          !is.na(lcc4)) %>%
            arrange(id) %>%
            dplyr::select(id, distanceMethod, tree, lcc2, lcc4)

        #9. Create design matrix----
        d <- x %>%
            dplyr::select(distanceMethod) %>%
            left_join(design, by="distanceMethod") %>%
            dplyr::select(-distanceMethod) %>%
            as.matrix()

        #10. Format abundance matrix----
        #add dummy variables for each of the columns to make sure the matrix is the right width
        y <- bird.i %>%
            dplyr::filter(id %in% x$id) %>%
            separate(distanceBand, into=c("start", "end"), sep="-", remove=FALSE) %>%
            mutate(end = as.numeric(str_sub(end, -100, -2))/100,
                   end = ifelse(is.na(end), Inf, end)) %>%
            left_join(design %>%
                          pivot_longer(d01:d13, values_to="end", names_to="position"),
                      by=c("distanceMethod", "end")) %>%
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
            mod <- try(cmulti(f, x, type="dis"))
            if (!inherits(mod, "try-error")) {
                dist <- data.frame(t(data.frame(mod["coefficients"]))) %>%
                    mutate(nobs=mod["nobs"]$nobs,
                           loglik = mod["loglik"]$loglik,
                           df = length(coef(mod)),
                           aic = AIC(mod),
                           aicc = aic + (2*df^2+2*df) / (nrow(y)-df-1),
                           model = modnames[j],
                           species = spp$speciesCode[i])
            } else {
                dist <- data.frame(nobs="try-error",
                                   species=spp$speciesCode[i])
            }
            mod.list[[j]] <- dist
        }

        #13. Save model results to species list---
        species.list[[i]] <- rbindlist(mod.list, fill=TRUE)

    }

    print(paste0("Finished modelling species ", spp$English_Name[i], ": ", i, " of ", nrow(spp), " species"))

}
