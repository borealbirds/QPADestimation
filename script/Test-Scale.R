# ---
# title: "QPAD estimation - test use of annual landcover"
# author: "Elly Knight"
# created: "Oct 6, 2022"
# ---

library(tidyverse) #basic data wrangling
library(detect) #removal models
library(data.table) #collapse list to dataframe
library(sf) #spatial manipulation
library(rgee) #GEE

#A. GET COVARIATES####

#1. Load data----
load("data/cleaned_data_2022-10-06.Rdata")

#2. Initialize rgee----
ee_Initialize()
ee_check()

#3. Set up to loop through data years----
years <- sort(unique(visit$year))

out.list <- list()
for(i in 1:length(years)){

    #4. Filter to year and subset into chunks of 5000----
    dat <- visit %>%
        st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
        dplyr::filter(year==years[i]) %>%
        mutate(loop = ceiling(row_number()/5000))

    dat.out <- data.frame()
    #5. Set up loop----
    for(j in 1:max(dat$loop)){

        #6. Format for rgee----
        dat.j <- dat %>%
            dplyr::filter(loop==j)

        dat.ee <- dat.j %>%
            dplyr::select(geometry) %>%
            sf_as_ee()

        #7. Get Hansen dataset for canopy cover----
        tree <- ee$Image('UMD/hansen/global_forest_change_2021_v1_9')

        dat.tree <- ee_extract(
            x=tree,
            y=dat.ee,
            scale=100,
            sf=FALSE
        ) %>%
            dplyr::select(-first_b30, -first_b40, -first_b50, -first_b70, -last_b30, -last_b40, -last_b50, -last_b70, -datamask)

        #8. Get Hermosilla landcover----
        yeartbl <- data.frame(year = years,
                              year.i = c(1988, 1990, 1990, 1993, 1993, 1995, 1995, 1997, 1997, 1999, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2007, 2007, 2008, 2010, 2010, 2010, 2010, 2015, 2015, 2016, 2018, 2018, 2019, 2019))
        year.i <- yeartbl$year.i[i]

        start<-paste0(year.i, "-01-01")
        end<-paste0(year.i,"-12-31")

        lc <- ee$ImageCollection('projects/sat-io/open-datasets/CA_FOREST_LC_VLCE2')$filterDate(start, end)

        dat.lc <- ee_extract(
            x=lc,
            y=dat.ee,
            scale=100,
            sf=FALSE
        )
        colnames(dat.lc) <- "hermosilla"

        #9. Get 2015 copernicus to fill gaps----
        start <- "2015-01-01"
        end <- "2015-12-31"

        cop <- ee$ImageCollection('COPERNICUS/Landcover/100m/Proba-V-C3/Global')$filterDate(start, end)

        dat.cop <- ee_extract(
            x=cop,
            y=dat.ee,
            scale=100,
            sf=FALSE
        )
        colnames(dat.cop) <- c("bare", "crop", "density", "copernicus", "prob", "forest", "grass", "moss", "shrub", "snow", "tree", "urban", "permanentwater", "seasonalwater")

        #10.Put everything together----
        dat.out <- data.frame(st_coordinates(dat.j)) %>%
            rename(lat = Y, lon = X) %>%
            cbind(data.frame(dat.j) %>%
                      dplyr::select(-geometry)) %>%
            cbind(dat.tree, dat.lc, dat.cop) %>%
            rbind(dat.out)

    }

    #11. Save to list----
    out.list[[i]] <- dat.out

    print(paste0("Finished year ", years[i], ": ", i, " of ", length(years), " years"))

}

#12. Create lookup tables for landcover classes----
#hermosilla codes
hermosillacodes <- data.frame(hermosilla = c(0, 20, 31, 32, 33, 40, 50, 80, 81, 100, 210, 220, 230),
                              hermosilladesc = c("unclassified", "water", "snow", "rock", "barren", "bryoid", "shrub", "wetland", "wetlandtreed", "herb", "conifer", "deciduous", "mixedwood"),
                              hermosilla2 = c(NA, "open", "open", "open", "open", "open", "open", "open", "treed", "open", "treed", "treed", "treed"),
                              hermosilla4 = c(NA, "open", "open", "open", "open", "open", "open", "wetland", "conifer", "open", "conifer", "decidmixed", "decidmixed"))

#copernicus codes
copernicuscodes <- data.frame(copernicus = c(0, 20, 30, 40, 50, 60, 70, 80, 90, 100, 111, 112, 113, 114, 115, 116, 121, 122, 123, 124, 125, 126, 200),
                              copernicusdesc = c("unclassified", "shrub", "herb", "crop", "urban", "barren", "snow", "water", "wetland", "bryoid", "coniferousclosed", "deciduousclosed", "coniferousclosed", "deciduousclosed", "mixedclosed", "treedclosed", "coniferousopen", "deciduousopen", "coniferousopen", "deciduousopen", "mixedopen", "treedopen", "marine"),
                              copernicus2 = c(NA, "open", "open", "open", "open", "open", "open", "open", "open", "open", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "open"),
                              copernicus4 = c(NA, "open", "open", "open", "open", "open", "open", "open", "wetland", "open", "conifer", "decidmixed", "conifer", "decidmixed", "decidmixed", "decidmixed", "conifer", "decidmixed", "conifer", "decidmixed", "decidmixed", "decidmixed", "open"))

#13. Tidy & create primary key----
tidy <- rbindlist(out.list, fill=TRUE) %>%
    data.frame() %>%
    mutate(lossyear = ifelse(is.na(lossyear), 0, lossyear+2000),
           cover = ifelse(lossyear >= year, 0, treecover2000),
           cover = ifelse(year >= 2012 & gain==1, 1, cover)) %>%
    left_join(hermosillacodes) %>%
    left_join(copernicuscodes) %>%
    mutate(lc2 = ifelse(is.na(hermosilla2), copernicus2, hermosilla2),
           lc4 = ifelse(is.na(hermosilla4), copernicus4, hermosilla4)) %>%
    dplyr::select(colnames(visit), cover, lc2, lc4)

visit <- tidy

#14. Save----
save(visit, bird, species,  file="data/cleaned_data_scale_2022-10-06.Rdata")
load("data/cleaned_data_scale_2022-10-06.Rdata")

#B. MODEL PERCEPTABILITY####

#1. Create list of models----
mods <- list(
    ~ 1,
    ~ cover,
    ~ lc2,
    ~ lc4,
    ~ lc2 + cover,
    ~ lc4 + cover)
names(mods) <- 0:5
modnames <- list(
    "0"="(Intercept)",
    "1"=c("(Intercept)", "cover"),
    "2"=c("(Intercept)", "lc2OpenWet"),
    "3"=c("(Intercept)", "lc4Conif", "lc4Open", "lc4Wet"),
    "4"=c("(Intercept)", "lc2OpenWet", "cover"),
    "5"=c("(Intercept)", "lc4Conif", "lc4Open", "lc4Wet", "cover"))

#2. Create design lookup table that describes duration method for each protocol----
#filter out distance methods that aren't appropriate for removal modelling (only have 1 bin)
distdesign <- visit %>%
    dplyr::select(distanceMethod) %>%
    unique() %>%
    dplyr::filter(!distanceMethod %in% c("UNKNOWN", "0m-INF", "0m-100m", "0m-400m", "0m-80m", "0m-50m", "0m-INF-ARU")) %>%
    separate(distanceMethod, into=c("d00", "d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), remove=FALSE, sep="-") %>%
    dplyr::select(-d00) %>%
    mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~str_sub(., -100, -2)) %>%
    mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~ifelse(.=="IN", Inf, .)) %>%
    mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10",  "d11", "d12", "d13"), ~as.numeric(.)/100)

#3. Get list of species to process----
spp <- species %>%
    dplyr::filter(Singing_birds==TRUE) %>%
    left_join(bird %>%
                  dplyr::select(speciesCode) %>%
                  unique()) %>%
    arrange(speciesCode)

#4. Set up loop for species----
percep <- list()
for(i in 1:nrow(spp)){

    #5. Filter abundance data for species---
    # filter out observations with unknown duration method or interval
    # filter to observations with covariates
    bird.i <- bird %>%
        dplyr::filter(speciesCode==spp$speciesCode[i],
                      distanceMethod %in% distdesign$distanceMethod,
                      distanceBand!="UNKNOWN") %>%
        group_by(id, distanceMethod, distanceBand) %>%
        summarize(abundance = sum(abundance)) %>%
        ungroup()

    #only model if there is data
    if(nrow(bird.i) > 0){

        #6. Replace abundance outliers (defined as 99% quantile across species) with 99% quantile of abundance----
        if(max(bird.i$abundance) > quantile(bird$abundance, 0.99)){
            bird.i <- bird.i %>%
                mutate(abundance = ifelse(abundance > quantile(abundance, 0.99), ceiling(quantile(abundance, 0.99)), abundance))
        }

        #7. Filter visit covariates----
        #Remove nas
        x <- visit %>%
            dplyr::filter(id %in% unique(bird.i$id),
                          !is.na(cover),
                          !is.na(lc2),
                          !is.na(lc4)) %>%
            arrange(id) %>%
            dplyr::select(id, distanceMethod, cover, lc2, lc4)

        #8. Create design matrix----
        d <- x %>%
            dplyr::select(distanceMethod) %>%
            left_join(distdesign, by="distanceMethod") %>%
            dplyr::select(-distanceMethod) %>%
            as.matrix()

        #9. Format abundance matrix----
        #add dummy variables for each of the columns to make sure the matrix is the right width
        y <- bird.i %>%
            dplyr::filter(id %in% x$id) %>%
            separate(distanceBand, into=c("start", "end"), sep="-", remove=FALSE) %>%
            mutate(end = as.numeric(str_sub(end, -100, -2))/100,
                   end = ifelse(is.na(end), Inf, end)) %>%
            left_join(distdesign %>%
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

        #10. Change zeros to NAs in the abundance matrix to match the design matrix----
        for (j in 1:nrow(y)){
            indices <- which(is.na(d[j,]))
            y[j, indices] <- NA
        }

        #11. Fit models----
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

        #12. Save model results to species list---
        percep[[i]] <- mod.list

    }

    print(paste0("Finished modelling species ", spp$English_Name[i], ": ", i, " of ", nrow(spp), " species"))

}

names(percep) <- spp$speciesCode

#13. Save out results----
save(percep, distdesign, file="results/perceptability_results_scale_2022-10-06.Rdata")

