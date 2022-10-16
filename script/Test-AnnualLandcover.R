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
        cover <- ee$Image('UMD/hansen/global_forest_change_2021_v1_9')

        dat.cover <- ee_extract(
            x=cover,
            y=dat.ee,
            scale=1000,
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
            scale=1000,
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
            scale=1000,
            sf=FALSE
        )
        colnames(dat.cop) <- c("bare", "crop", "density", "copernicus", "prob", "forest", "grass", "moss", "shrub", "snow", "cover", "urban", "permanentwater", "seasonalwater")

        #10.Put everything together----
        dat.out <- data.frame(st_coordinates(dat.j)) %>%
            rename(lat = Y, lon = X) %>%
            cbind(data.frame(dat.j) %>%
                      dplyr::select(-geometry)) %>%
            cbind(dat.cover, dat.lc, dat.cop) %>%
            rbind(dat.out)

    }

    #11. Save to list----
    out.list[[i]] <- dat.out

    print(paste0("Finished year ", years[i], ": ", i, " of ", length(years), " years"))

}

#12. Create lookup tables for landcover classes----
#hermosilla codes
hermosillacodes <- data.frame(hermosilla = c(0, 20, 31, 32, 33, 40, 50, 80, 81, 100, 210, 220, 230),
                              hermosilladesc = c("unclassified", "water", "snow", "rock", "barren", "bryoid", "shrub", "wetland", "wetlandcoverd", "herb", "conifer", "deciduous", "mixedwood"),
                              hermosilla2 = c(NA, "open", "open", "open", "open", "open", "open", "open", "coverd", "open", "coverd", "coverd", "coverd"),
                              hermosilla4 = c(NA, "open", "open", "open", "open", "open", "open", "wetland", "conifer", "open", "conifer", "decidmixed", "decidmixed"))

#copernicus codes
copernicuscodes <- data.frame(copernicus = c(0, 20, 30, 40, 50, 60, 70, 80, 90, 100, 111, 112, 113, 114, 115, 116, 121, 122, 123, 124, 125, 126, 200),
                              copernicusdesc = c("unclassified", "shrub", "herb", "crop", "urban", "barren", "snow", "water", "wetland", "bryoid", "coniferousclosed", "deciduousclosed", "coniferousclosed", "deciduousclosed", "mixedclosed", "coverdclosed", "coniferousopen", "deciduousopen", "coniferousopen", "deciduousopen", "mixedopen", "coverdopen", "marine"),
                              copernicus2 = c(NA, "open", "open", "open", "open", "open", "open", "open", "open", "open", "coverd", "coverd", "coverd", "coverd", "coverd", "coverd", "coverd", "coverd", "coverd", "coverd", "coverd", "coverd", "open"),
                              copernicus4 = c(NA, "open", "open", "open", "open", "open", "open", "open", "wetland", "open", "conifer", "decidmixed", "conifer", "decidmixed", "decidmixed", "decidmixed", "conifer", "decidmixed", "conifer", "decidmixed", "decidmixed", "decidmixed", "open"))

#13. Tidy & create primary key----
tidy <- rbindlist(out.list, fill=TRUE) %>%
    data.frame() %>%
    mutate(lossyear = ifelse(is.na(lossyear), 0, lossyear+2000),
           cover = ifelse(lossyear >= year, 0, covercover2000),
           cover = ifelse(year >= 2012 & gain==1, 1, cover)) %>%
    left_join(hermosillacodes) %>%
    left_join(copernicuscodes) %>%
    mutate(lc2 = ifelse(is.na(hermosilla2), copernicus2, hermosilla2),
           lc4 = ifelse(is.na(hermosilla4), copernicus4, hermosilla4)) %>%
    dplyr::select(colnames(visit), cover, lc2, lc4)

visit <- tidy

#14. Save----
save(visit, bird, species,  file="data/cleaned_data_annualLC_2022-10-06.Rdata")
load("data/cleaned_data_annualLC_2022-10-06.Rdata")

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
save(percep, distdesign, file="results/perceptability_results_annualLC_2022-10-06.Rdata")

#C. PACKAGE####

load("results/perceptability_results_annualLC_2022-10-06.Rdata")
load("results/BAMCOEFS_QPAD_v4.rda")

#1. Remove species that didn't have enough data----
resDisOK <- percep[!sapply(percep, length)==0]
c(OK=length(resDisOK), failed=length(percep)-length(resDisOK), all=length(percep))

#2. Remove species where null model failed----
resDis <- resDisOK[!t(sapply(resDisOK, function(z) ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))[,1]==0]

#3. Create 0/1 table for model fit----
edr_mod <- t(sapply(resDis, function(z)
    ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))

#4. Adjust lists to include all species there's an estimate for----
tmp <- rownames(edr_mod)
edr_models <- matrix(0L, length(tmp), ncol(edr_mod))
dimnames(edr_models) <- list(tmp, colnames(edr_mod))
edr_models[rownames(edr_mod),] <- edr_mod

#4. Check for dropped factor levels & min # of detections within each class for detectability models----
n.min.class <- 5 # min number of detections within each class
for (spp in rownames(edr_mod)) {
    ## data for checking detections in classes
    bird.i <- bird %>%
        dplyr::filter(speciesCode==spp,
                      distanceMethod %in% distdesign$distanceMethod,
                      distanceBand!="UNKNOWN") %>%
        dplyr::select(id) %>%
        unique()
    Dat <- visit %>%
        dplyr::filter(id %in% unique(bird.i$id),
                      !is.na(tree),
                      !is.na(lcc2),
                      !is.na(lcc4)) %>%
        arrange(id) %>%
        dplyr::select(lcc2, lcc4)
    for (mid in colnames(edr_models)) {
        if (!inherits(resDis[[spp]][[mid]], "try-error")) {
            lcf <- length(resDis[[spp]][[mid]]$coefficients)
            lnm <- length(resDis[[spp]][[mid]]$names)
            if (lcf != lnm) {
                cat("EDR conflict for", spp, "model", mid,
                    "( len.coef =", lcf, ", len.name =", lnm, ")\n")
                edr_models[spp,mid] <- 0
            } else {
                if (mid %in% c("2", "4") && min(table(Dat$lcc2)) < n.min.class) {
                    cat("EDR LCC2 min issue for", spp, "model", mid, "\n")
                    edr_models[spp,mid] <- 0
                }
                if (mid %in% c("3", "5") && min(table(Dat$lcc4)) < n.min.class) {
                    cat("EDR LCC4 min issue for", spp, "model", mid, "\n")
                    edr_models[spp,mid] <- 0
                }
                if (mid %in% c("1", "4", "5") &&
                    resDis[[spp]][[mid]]$coefficients["log.tau_cover"] > 0) {
                    cat("EDR TREE > 0 issue for", spp, "model", mid, "\n")
                    edr_models[spp,mid] <- 0
                }
            }
        } else {
            resDis[[spp]][[mid]] <- structure("Error", class = "try-error")
            attributes(resDis[[spp]][[mid]]) <- NULL
            class(resDis[[spp]][[mid]]) <- "try-error"
        }
        flush.console()
    }
}

#6. Exclude species with no null model----
edr_models[edr_models[,1]==0,] <- 0L

#8. Get number of models----
edr_nmod <- ncol(edr_mod)

#9. Get sample sizes----
edr_n <- numeric(length(tmp))
names(edr_n) <- tmp
edr_nn <- sapply(resDis, function(z) ifelse(inherits(z[["0"]], "try-error"),
                                            NA, z[["0"]]$nobs))
edr_n[names(edr_nn)] <- edr_nn

#10. exclude all models for species with < n.con observations----
n.con <- 25
edr_models[edr_n < n.con, ] <- 0L

#11. Exclude everything but null for species with n.con < observations < n.min----
n.min <- 75
edr_models[edr_n < n.min & edr_n >= n.con, 2:ncol(edr_models)] <- 0L

#12. ID species to keep----
spp <- tmp[rowSums(edr_models) > 0]
length(spp)

edr_models <- edr_models[spp,]
edr_models <- edr_models[.BAMCOEFS4$spp,]

edr_n <- edr_n[.BAMCOEFS4$spp]

#13. Get number of parameters----
#use OVEN as template
edr_df <- sapply(resDis[["OVEN"]][1:edr_nmod], "[[", "p")

#14. Get estimates----
edr_estimates <- resDis[spp]

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
edr_list <- sapply(edr_estimates[["OVEN"]], function(z) paste(z$names, collapse=" + "))

#16. Get loglik values---
edr_loglik <- edr_models
edr_loglik[] <- -Inf
for (i in .BAMCOEFS4$spp) { # species
    for (j in 1:edr_nmod) { # models
        if (edr_models[i,j] > 0)
            edr_loglik[i,j] <- resDis[[i]][[j]]$loglik
    }
}

#17. Get AIC values----
edr_aic <- edr_aicc <- edr_bic <- edr_loglik
edr_aic[] <- Inf
edr_aicc[] <- Inf
edr_bic[] <- Inf
for (i in .BAMCOEFS4$spp) {
    edr_aic[i,] <- -2*edr_loglik[i,] + 2*edr_df
    edr_aicc[i,] <- edr_aic[i,] + (2*edr_df*(edr_df+1)) / (edr_n[i]-edr_df-1)
    edr_bic[i,] <- -2*edr_loglik[i,] + log(edr_n[i])*edr_df
}

#18. Rank models----
edr_aicrank <- t(apply(edr_aic, 1, rank))*edr_models
edr_aicrank[edr_aicrank==0] <- NA

edr_aiccrank <- t(apply(edr_aicc, 1, rank))*edr_models
edr_aiccrank[edr_aiccrank==0] <- NA

edr_bicrank <- t(apply(edr_bic, 1, rank))*edr_models
edr_bicrank[edr_bicrank==0] <- NA

edr_aicbest <- apply(edr_aicrank, 1, function(z) colnames(edr_models)[which.min(z)])
edr_aiccbest <- apply(edr_aiccrank, 1, function(z) colnames(edr_models)[which.min(z)])
edr_bicbest <- apply(edr_bicrank, 1, function(z) colnames(edr_models)[which.min(z)])

#19. Set version----
version <- "4lc"

#20. Bundle----
bamcoefs <- list(spp=.BAMCOEFS4$spp,
                 spp_table=spp_table,
                 edr_list=edr_list,
                 sra_list=.BAMCOEFS4$sra_list,
                 edr_models=edr_models,
                 sra_models=.BAMCOEFS4$sra_models,
                 edr_n=edr_n,
                 sra_n=.BAMCOEFS4$sra_n,
                 edr_df=edr_df,
                 sra_df=.BAMCOEFS4$sra_df,
                 edr_loglik=edr_loglik,
                 sra_loglik=.BAMCOEFS4$sra_loglik,
                 edr_aic=edr_aic,
                 sra_aic=.BAMCOEFS4$sra_aic,
                 edr_aicc=edr_aicc,
                 sra_aicc=.BAMCOEFS4$sra_aicc,
                 edr_bic=edr_bic,
                 sra_bic=.BAMCOEFS4$sra_bic,
                 edr_aicrank=edr_aicrank,
                 sra_aicrank=.BAMCOEFS4$sra_aicrank,
                 edr_aiccrank=edr_aiccrank,
                 sra_aiccrank=.BAMCOEFS4$sra_aiccrank,
                 edr_bicrank=edr_bicrank,
                 sra_bicrank=.BAMCOEFS4$sra_bicrank,
                 edr_aicbest=edr_aicbest,
                 sra_aicbest=.BAMCOEFS4$sra_aicbest,
                 edr_aiccbest=edr_aiccbest,
                 sra_aiccbest=.BAMCOEFS4$sra_aiccbest,
                 edr_bicbest=edr_bicbest,
                 sra_bicbest=.BAMCOEFS4$sra_bicbest,
                 edr_estimates=edr_estimates,
                 sra_estimates=.BAMCOEFS4$sra_estimates,
                 version=version)
.BAMCOEFS4lc <- list2env(bamcoefs)

save(.BAMCOEFS4lc, file="results/BAMCOEFS_QPAD_v4_lc.rda")

#D. COMPARE####

load("results/BAMCOEFS_QPAD_v4.rda")
load("results/BAMCOEFS_QPAD_v4_lc.rda")

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
names <- data.frame(name = .BAMCOEFS4$edr_list,
                    mod = names(.BAMCOEFS4$edr_list))

#1. Best model----
mod <- data.frame(mod4 = .BAMCOEFS4$edr_aicbest,
                    mod4lc = .BAMCOEFS4lc$edr_aicbest,
                    spp = .BAMCOEFS4$spp) %>%
    mutate(same = ifelse(mod4==mod4lc, 1, 0)) %>%
    left_join(names %>%
                  rename(mod4=mod, name4=name)) %>%
    left_join(names %>%
                  rename(mod4lc = mod, name4lc = name))
table(mod$same)
table(mod$mod4)
table(mod$mod4lc)

#2. Tree effect----
#Multiply hermosilla by 100 to put on same scale
treedf <-data.frame()
for(i in 1:length(.BAMCOEFS4$spp)){
    tree <- .BAMCOEFS4$edr_estimates[[i]]$`1`$coefficients[2]
    treelc <- .BAMCOEFS4lc$edr_estimates[[i]]$`1`$coefficients[2]*100
    treedf <- rbind(treedf,
                  data.frame(tree=tree, treelc=treelc, .BAMCOEFS4$spp[[i]]))
}
rownames(treedf) <- NULL
treedf$conflict <- ifelse(treedf$tree > 0, 1, 0)
treedf$conflictlc <- ifelse(treedf$treelc > 0, 1, 0)
treedf$conflictv <- case_when((treedf$conflict==1 & treedf$conflictlc==1) ~ "Both",
                              (treedf$conflict==1 & treedf$conflictlc==0) ~ "V4",
                              (treedf$conflict==0 & treedf$conflictlc==1) ~ "Annual LC",
                              !is.na(treedf$conflict) ~ "Neither")


treen <- data.frame(table(treedf$conflictv)) %>%
    rename(conflictv = Var1) %>%
    mutate(text = paste0("n = ", Freq),
           x = c(-1,1,-1,1),
           y=c(1,1,-1,-1))

ggplot(treedf) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(aes(x=tree, y=treelc, fill=conflictv), pch=21, alpha = 0.5, size=4) +
    geom_text(data=treen, aes(x=x, y=y, label=text)) +
    xlab("Existing cover dataset effect on perceptibility") +
    ylab("Hansen annual cover effect on perceptibility") +
    scale_fill_viridis_d(name="Beta > 0") +
    ylim(c(-1.5,1.5)) +
    xlim(c(-1.5,1.5)) +
    my.theme

ggsave(filename="figures/AnnualLC_tree.jpeg", width =7, height=6)


#3. Open effect----
opendf <-data.frame()
for(i in 1:length(.BAMCOEFS4$spp)){
    rm(open, openlc)
    if(class(.BAMCOEFS4$edr_estimates[[i]]$`2`)!="try-error"){
        open <- .BAMCOEFS4$edr_estimates[[i]]$`2`$coefficients[2]
    }
    if(class(.BAMCOEFS4lc$edr_estimates[[i]]$`2`)!="try-error"){
        openlc <- .BAMCOEFS4lc$edr_estimates[[i]]$`2`$coefficients[2]
    }
    opendf <- rbind(opendf,
                    data.frame(open=ifelse(exists("open"), open, NA),
                               openlc=ifelse(exists("openlc"), openlc, NA),
                               .BAMCOEFS4$spp[[i]]))
}
rownames(opendf) <- NULL

opendf$conflict <- ifelse(opendf$open < 0, 1, 0)
opendf$conflictlc <- ifelse(opendf$openlc < 0, 1, 0)
opendf$conflictv <- case_when((opendf$conflict==1 & opendf$conflictlc==1) ~ "Both",
                              (opendf$conflict==1 & opendf$conflictlc==0) ~ "V4",
                              (opendf$conflict==0 & opendf$conflictlc==1) ~ "Annual LC",
                              !is.na(opendf$conflict) ~ "Neither")


openn <- data.frame(table(opendf$conflictv)) %>%
    rename(conflictv = Var1) %>%
    mutate(text = paste0("n = ", Freq),
           x = c(-1,1,-1,1),
           y=c(1,1,-1,-1))

ggplot(opendf) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(aes(x=open, y=openlc, fill=conflictv), pch=21, alpha = 0.5, size=4) +
    geom_text(data=openn, aes(x=x, y=y, label=text)) +
    xlab("LCC open effect on perceptibility") +
    ylab("Hermosilla annual open effect on perceptibility") +
    scale_fill_viridis_d(name="Beta < 0") +
    ylim(c(-1.5,1.5)) +
    xlim(c(-1.5,1.5)) +
    my.theme

ggsave(filename="figures/AnnualLC_open.jpeg", width =7, height=6)
