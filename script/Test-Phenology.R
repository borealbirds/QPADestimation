# ---
# title: "QPAD estimation - test use of other phenology variables"
# author: "Elly Knight"
# created: "Oct 6, 2022"
# ---

library(tidyverse) #basic data wrangling
library(lubridate) #date wrangling
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

        #7. Get modis greenup day----
        if(years[i] < 2001){year.i <- 2001}
        if(years[i] >= 2001 & years[i] <= 2019){year.i <- years[i]}
        if(years[i] > 2019){year.i <-2019}

        start<-paste0(year.i, "-01-01")
        end<-paste0(year.i,"-12-31")

        green <- ee$ImageCollection('MODIS/006/MCD12Q2')$select(c('Greenup_1', 'MidGreenup_1', 'Peak_1', 'Maturity_1'))$filterDate(start, end)

        dat.green <- ee_extract(
            x=green,
            y=dat.ee,
            scale=100,
            sf=FALSE
        )
        colnames(dat.green) <- c("green15", "green50", "greenmax", "green90")

        #8. Get mean annual temperature----
        temp <- ee$ImageCollection('ECMWF/ERA5/MONTHLY')$select('mean_2m_air_temperature')$filterDate(start, end)$mean()

        dat.temp <- ee_extract(
            x=temp,
            y=dat.ee,
            scale=100,
            sf=FALSE
        )
        colnames(dat.temp) <- c("annualtemp")

        #9. Get mean temperature----
        tempall <- ee$ImageCollection('ECMWF/ERA5/MONTHLY')$select('mean_2m_air_temperature')$mean()

        dat.tempall <- ee_extract(
            x=tempall,
            y=dat.ee,
            scale=100,
            sf=FALSE
        )
        colnames(dat.tempall) <- c("temp")

        #10.Put everything together----
        dat.out <- data.frame(st_coordinates(dat.j)) %>%
            rename(lat = Y, lon = X) %>%
            cbind(data.frame(dat.j) %>%
                      dplyr::select(-geometry)) %>%
            cbind(dat.green, dat.temp, dat.tempall) %>%
            rbind(dat.out)

    }

    #11. Save to list----
    out.list[[i]] <- dat.out

    print(paste0("Finished year ", years[i], ": ", i, " of ", length(years), " years"))

}

#12. Tidy, standardize covs, create polynomials----
tidy <- rbindlist(out.list, fill=TRUE) %>%
    data.frame() %>%
    mutate(dsgyear = case_when(year < 2001 ~ 2001,
                               year >= 2001 & year <= 2019 ~ year,
                               year > 2019 ~ 2019)) %>%
    mutate(dsg15 = green15 - as.numeric(ymd(paste(dsgyear, "-01-01"))),
           dsg50 = green50 - as.numeric(ymd(paste(dsgyear, "-01-01"))),
           dsg90 = green90 - as.numeric(ymd(paste(dsgyear, "-01-01"))),
           dsgmax = greenmax - as.numeric(ymd(paste(dsgyear, "-01-01")))) %>%
    mutate(tsg15 = dsg15/365,
           tsg50 = dsg15/365,
           tsg90 = dsg15/365,
           tsgmax = dsg15/365,
           jday = julian/365,
           tssr = hssr/24,
           tsg = tsg/365,
           jday2 = jday^2,
           tssr2 = tssr^2,
           tsg152 = tsg15^2,
           tsg502 = tsg50^2,
           tsg902 = tsg90^2,
           tsgmax2 = tsgmax^2,
           temps = (temp-223.15)/100,
           tempannuals = (annualtemp-223.15)/100)

visit <- tidy

#13. Save----
save(visit, bird, species,  file="data/cleaned_data_phenology_2022-10-06.Rdata")
load("data/cleaned_data_phenology_2022-10-06.Rdata")

#B. MODEL AVAILABILITY####

#1. Create list of models for day of year----
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
    ~ tsg + tsg2 + tssr + tssr2,
    ~ jday + jday:temps,
    ~ jday + jday2 + jday:temps + jday2:temps,
    ~ jday + jday:temps + tssr,
    ~ jday + jday2 + jday:temps + jday2:temps + tssr,
    ~ jday + jday:temps + tssr + tssr2,
    ~ jday + jday2 + jday:temps + jday2:temps + tssr + tssr2,
    ~ jday + jday:tempannuals,
    ~ jday + jday2 + jday:tempannuals + jday2:tempannuals,
    ~ jday + jday:tempannuals + tssr,
    ~ jday + jday2 + jday:tempannuals + jday2:tempannuals + tssr,
    ~ jday + jday:tempannuals + tssr + tssr2,
    ~ jday + jday2 + jday:tempannuals + jday2:tempannuals + tssr + tssr2,
    ~ jday + jday:lat,
    ~ jday + jday2 + jday:lat + jday2:lat,
    ~ jday + jday:lat + tssr,
    ~ jday + jday2 + jday:lat + jday2:lat + tssr,
    ~ jday + jday:lat + tssr + tssr2,
    ~ jday + jday2 + jday:lat + jday2:lat + tssr + tssr2,
    ~ jday + jday:lat:lon,
    ~ jday + jday2 + jday:lat:lon + jday2:lat:lon,
    ~ jday + jday:lat:lon + tssr,
    ~ jday + jday2 + jday:lat:lon + jday2:lat:lon + tssr,
    ~ jday + jday:lat:lon + tssr + tssr2,
    ~ jday + jday2 + jday:lat:lon + jday2:lat:lon + tssr + tssr2,
    ~ tsg15,
    ~ tsg15 + tsg152,
    ~ tsg15 + tssr,
    ~ tsg15 + tsg152 + tssr,
    ~ tsg15 + tssr + tssr2,
    ~ tsg15 + tsg152 + tssr + tssr2,
    ~ tsg50,
    ~ tsg50 + tsg502,
    ~ tsg50 + tssr,
    ~ tsg50 + tsg502 + tssr,
    ~ tsg50 + tssr + tssr2,
    ~ tsg50 + tsg502 + tssr + tssr2)
names(mods) <- 0:50
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
    "14"=c("(Intercept)", "tsg", "tsg2", "tssr", "tssr2"),
    "15"=c("(Intercept)", "jday", "jday:temps"),
    "16"=c("(Intercept)", "jday", "jday:temps", "jday2", "jday2:temps"),
    "17"=c("(Intercept)", "jday", "jday:temps", "tssr"),
    "18"=c("(Intercept)", "jday", "jday:temps", "jday2", "jday2:temps", "tssr"),
    "19"=c("(Intercept)", "jday", "jday:temps", "tssr", "tssr2"),
    "20"=c("(Intercept)", "jday", "jday:temps", "jday2", "jday2:temps", "tssr", "tssr2"),
    "21"=c("(Intercept)", "jday", "jday:tempannuals"),
    "22"=c("(Intercept)", "jday", "jday:tempannuals", "jday2", "jday2:tempannuals"),
    "23"=c("(Intercept)", "jday", "jday:tempannuals", "tssr"),
    "24"=c("(Intercept)", "jday", "jday:tempannuals", "jday2", "jday2:tempannuals", "tssr"),
    "25"=c("(Intercept)", "jday", "jday:tempannuals", "tssr", "tssr2"),
    "26"=c("(Intercept)", "jday", "jday:tempannuals", "jday2", "jday2:tempannuals", "tssr", "tssr2"),
    "27"=c("(Intercept)", "jday", "jday:lat"),
    "28"=c("(Intercept)", "jday", "jday:lat", "jday2", "jday2:lat"),
    "29"=c("(Intercept)", "jday", "jday:lat", "tssr"),
    "30"=c("(Intercept)", "jday", "jday:lat", "jday2", "jday2:lat", "tssr"),
    "31"=c("(Intercept)", "jday", "jday:lat", "tssr", "tssr2"),
    "32"=c("(Intercept)", "jday", "jday:lat", "jday2", "jday2:lat", "tssr", "tssr2"),
    "33"=c("(Intercept)", "jday", "jday:lat:lon"),
    "34"=c("(Intercept)", "jday", "jday:lat:lon", "jday2", "jday2:lat:lon"),
    "35"=c("(Intercept)", "jday", "jday:lat:lon", "tssr"),
    "36"=c("(Intercept)", "jday", "jday:lat:lon", "jday2", "jday2:lat:lon", "tssr"),
    "37"=c("(Intercept)", "jday", "jday:lat:lon", "tssr", "tssr2"),
    "38"=c("(Intercept)", "jday", "jday:lat:lon", "jday2", "jday2:lat:lon", "tssr", "tssr2"),
    "39"=c("(Intercept)", "tsg15"),
    "40"=c("(Intercept)", "tsg15", "tsg152"),
    "41"=c("(Intercept)", "tsg15", "tssr"),
    "42"=c("(Intercept)", "tsg15", "tsg152", "tssr"),
    "43"=c("(Intercept)", "tsg15", "tssr", "tssr2"),
    "44"=c("(Intercept)", "tsg15", "tsg152", "tssr", "tssr2"),
    "45"=c("(Intercept)", "tsg50"),
    "46"=c("(Intercept)", "tsg50", "tsg502"),
    "47"=c("(Intercept)", "tsg50", "tssr"),
    "48"=c("(Intercept)", "tsg50", "tsg502", "tssr"),
    "49"=c("(Intercept)", "tsg50", "tssr", "tssr2"),
    "50"=c("(Intercept)", "tsg50", "tsg502", "tssr", "tssr2"))

#3. Create design lookup table that describes duration method for each protocol----
#filter out duration methods that aren't appropriate for removal modelling (only have 1 time bin)
durdesign <- visit %>%
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
avail <- list()
for(i in 1:nrow(spp)){

    #6. Filter abundance data for species---
    # filter out observations with unknown duration method or interval
    # filter to observations with covariates
    bird.i <- bird %>%
        dplyr::filter(speciesCode==spp$speciesCode[i],
                      durationMethod %in% durdesign$durationMethod,
                      durationInterval!="UNKNOWN") %>%
        group_by(id, durationMethod, durationInterval) %>%
        summarize(abundance = sum(abundance)) %>%
        ungroup()

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
                      !is.na(tsg15),
                      !is.na(temps),
                      !is.na(tempannuals)) %>%
        arrange(id) %>%
        dplyr::select(id, durationMethod, lat, lon, jday, jday2, tssr, tssr2, tsg, tsg2, tsg15, tsg152, tsg50, tsg502, tsg90, tsg902, tsgmax, tsgmax2, temps, tempannuals)

    #only model if there is data

        if(nrow(x) > 0){

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
            separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE, extra="drop", fill="right") %>%
            mutate(start = as.numeric(start),
                   end = as.numeric(str_sub(end, -100, -4))) %>%
            left_join(durdesign %>%
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

names(avail) <- spp$speciesCode

#14. Save out results----
save(durdesign, avail, file="results/availability_results_phenology_2022-10-06.Rdata")

#C. PACKAGE####

load("results/availability_results_phenology_2022-10-06.Rdata")
load("results/BAMCOEFS_QPAD_v4.rda")

#1. Remove species that didn't have enough data----
resDurOK <- avail[!sapply(avail, length)==0]
c(OK=length(resDurOK), failed=length(avail)-length(resDurOK), all=length(avail))

#2. Remove species where null model failed---
resDur <- resDurOK[!t(sapply(resDurOK, function(z) ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))[,1]==0]

#3. Create 0/1 table for model fit----
sra_mod <- t(sapply(resDur, function(z)
    ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))

#3. Adjust lists to include all species there's an estimate for----
tmp <- rownames(sra_mod)

sra_models <- matrix(0L, length(tmp), ncol(sra_mod))
dimnames(sra_models) <- list(tmp, colnames(sra_mod))
sra_models[rownames(sra_mod),] <- sra_mod

#4. Check for dropped factor levels in removal models----
for (spp in rownames(sra_mod)) {
    for (mid in colnames(sra_models)) {
        if (!inherits(resDur[[spp]][[mid]], "try-error")) {
            lcf <- length(resDur[[spp]][[mid]]$coefficients)
            lnm <- length(resDur[[spp]][[mid]]$names)
            if (lcf != lnm) {
                cat("SRA conflict for", spp, "model", mid,
                    "( len.coef =", lcf, ", len.name =", lnm, ")\n")
                sra_models[spp,mid] <- 0
            }
        } else {
            resDur[[spp]][[mid]] <- structure("Error", class = "try-error")
        }
        flush.console()
    }
}

#6. Exclude species with no null model----
sra_models[sra_models[,1]==0,] <- 0L

#7. Exclude species with phi estimate < 0.01----
phi0 <- sapply(resDur, function(z) exp(z[["0"]]$coefficients))
names(phi0) <- names(resDur)
sra_models[names(phi0)[phi0 < 0.01],] <- 0L

#8. Get number of models----
sra_nmod <- ncol(sra_mod)

#9. Get sample sizes----
sra_n <- numeric(length(tmp))
names(sra_n) <- tmp
sra_nn <- sapply(resDur, function(z) ifelse(inherits(z[["0"]], "try-error"),
                                            NA, z[["0"]]$nobs))
sra_n[names(sra_nn)] <- sra_nn

#10. exclude all models for species with < n.con observations----
n.con <- 25
sra_models[sra_n < n.con, ] <- 0L

#11. Exclude everything but null for species with n.con < observations < n.min----
n.min <- 75
sra_models[sra_n < n.min & sra_n >= n.con, 2:ncol(sra_models)] <- 0L

#12. ID species to keep----
spp <- tmp[rowSums(sra_models) > 0]
length(spp)

sra_models <- sra_models[spp,]
sra_models <- sra_models[.BAMCOEFS4$spp,]

sra_n <- sra_n[spp]

#13. Get number of parameters----
#use OVEN as template
sra_df <- sapply(resDur[["OVEN"]][1:sra_nmod], "[[", "p")

#14. Get estimates----
sra_estimates <- resDur[spp]

#15. Get species table----
tax <- read.csv("data/taxonomytable.csv")
tax <- tax[!duplicated(tax$Species_ID),]
rownames(tax) <- tax$Species_ID
spp_table <- data.frame(spp=spp,
                        scientific_name=tax[spp, "Scientific_Name"],
                        common_name=tax[spp, "English_Name"])
rownames(spp_table) <- spp
spp_table <- droplevels(spp_table)

#16. Get variable names for different models----
#use OVEN as template
sra_list <- sapply(sra_estimates[["OVEN"]], function(z) paste(z$names, collapse=" + "))

#16. Get loglik values---
sra_loglik <- sra_models
sra_loglik[] <- -Inf
for (i in .BAMCOEFS4$spp) { # species
    for (j in 1:sra_nmod) { # models
        if (sra_models[i,j] > 0)
            sra_loglik[i,j] <- resDur[[i]][[j]]$loglik
    }
}

#17. Get AIC values----
sra_aic <- sra_aicc <- sra_bic <- sra_loglik
sra_aic[] <- Inf
sra_aicc[] <- Inf
sra_bic[] <- Inf
for (i in .BAMCOEFS4$spp) {
    sra_aic[i,] <- -2*sra_loglik[i,] + 2*sra_df
    sra_aicc[i,] <- sra_aic[i,] + (2*sra_df*(sra_df+1)) / (sra_n[i]-sra_df-1)
    sra_bic[i,] <- -2*sra_loglik[i,] + log(sra_n[i])*sra_df
}

#18. Rank models----
sra_aicrank <- t(apply(sra_aic, 1, rank))*sra_models
sra_aicrank[sra_aicrank==0] <- NA

sra_aiccrank <- t(apply(sra_aicc, 1, rank))*sra_models
sra_aiccrank[sra_aiccrank==0] <- NA

sra_bicrank <- t(apply(sra_bic, 1, rank))*sra_models
sra_bicrank[sra_bicrank==0] <- NA

sra_aicbest <- apply(sra_aicrank, 1, function(z) colnames(sra_models)[which.min(z)])
sra_aiccbest <- apply(sra_aiccrank, 1, function(z) colnames(sra_models)[which.min(z)])
sra_bicbest <- apply(sra_bicrank, 1, function(z) colnames(sra_models)[which.min(z)])

#19. Set version----
version <- "4phen"

#20. Bundle----
bamcoefs <- list(spp=spp,
                 spp_table=spp_table,
                 edr_list=.BAMCOEFS4$edr_list,
                 sra_list=sra_list,
                 edr_models=.BAMCOEFS4$edr_models,
                 sra_models=sra_models,
                 edr_n=.BAMCOEFS4$edr_n,
                 sra_n=sra_n,
                 edr_df=.BAMCOEFS4$edr_df,
                 sra_df=sra_df,
                 edr_loglik=.BAMCOEFS4$edr_loglik,
                 sra_loglik=sra_loglik,
                 edr_aic=.BAMCOEFS4$edr_aic,
                 sra_aic=sra_aic,
                 edr_aicc=.BAMCOEFS4$edr_aicc,
                 sra_aicc=sra_aicc,
                 edr_bic=.BAMCOEFS4$edr_bic,
                 sra_bic=sra_bic,
                 edr_aicrank=.BAMCOEFS4$edr_aicrank,
                 sra_aicrank=sra_aicrank,
                 edr_aiccrank=.BAMCOEFS4$edr_aiccrank,
                 sra_aiccrank=sra_aiccrank,
                 edr_bicrank=.BAMCOEFS4$edr_bicrank,
                 sra_bicrank=sra_bicrank,
                 edr_aicbest=.BAMCOEFS4$edr_aicbest,
                 sra_aicbest=sra_aicbest,
                 edr_aiccbest=.BAMCOEFS4$edr_aiccbest,
                 sra_aiccbest=sra_aiccbest,
                 edr_bicbest=.BAMCOEFS4$edr_bicbest,
                 sra_bicbest=sra_bicbest,
                 edr_estimates=.BAMCOEFS4$edr_estimates,
                 sra_estimates=sra_estimates,
                 version=version)
.BAMCOEFS4phen <- list2env(bamcoefs)

save(.BAMCOEFS4phen, file="results/BAMCOEFS_QPAD_v4_phenology.rda")

#D. COMPARE####

load("results/BAMCOEFS_QPAD_v4.rda")
load("results/BAMCOEFS_QPAD_v4_phenology.rda")

my.theme <- theme_classic() +
    theme(text=element_text(size=12, family="Arial"),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(margin=margin(10,0,0,0)),
          axis.title.y=element_text(margin=margin(0,10,0,0)),
          axis.line.x=element_line(linetype=1),
          axis.line.y=element_line(linetype=1),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          plot.title=element_text(size=12, hjust = 0.5))

#Recall list of models
names4 <- data.frame(name = .BAMCOEFS4$sra_list,
                    mod = names(.BAMCOEFS4$sra_list))
namesphen <- data.frame(name = .BAMCOEFS4phen$sra_list,
                     mod = names(.BAMCOEFS4phen$sra_list)) %>%
    mutate(group = c("Intercept", "day", "tssr only", "day", "tssr only", rep("day", 4), rep("seedgrow", 6), rep("annual temp", 6), rep("mean temp", 6), rep("lat", 6), rep("lat:lon", 6), rep("15% evi", 6), rep("50% evi", 6)))

#1. Best model----
mod <- data.frame(mod4 = .BAMCOEFS4$sra_aicbest,
                  modphen = .BAMCOEFS4phen$sra_aicbest,
                  spp = .BAMCOEFS4$spp) %>%
    mutate(same = ifelse(mod4==modphen, 1, 0)) %>%
    left_join(namesphen %>%
                  rename(mod4 = mod, name4 = name, group4 = group)) %>%
    left_join(namesphen %>%
                  rename(modphen = mod, namephen = name, groupphen = group))

modn.mod <- mod %>%
    group_by(mod4) %>%
    summarize(n4=n()) %>%
    rename(mod = mod4) %>%
    full_join(mod %>%
                  group_by(modphen) %>%
                  summarize(nphen=n()) %>%
                  rename(mod = modphen)) %>%
    full_join(namesphen) %>%
    mutate(mod = as.numeric(mod)) %>%
    arrange(mod)

modn.group <- mod %>%
    group_by(group4) %>%
    summarize(n4=n()) %>%
    rename(group = group4) %>%
    full_join(mod %>%
                  group_by(groupphen) %>%
                  summarize(nphen=n()) %>%
                  rename(group = groupphen))
