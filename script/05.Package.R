# ---
# title: "QPAD estimation - packaging"
# author: "Elly Knight"
# created: "October 11, 2022"
# adapted from: "QPAD version 3 documentation" by Peter Solymos https://github.com/borealbirds/bamanalytics/blob/master/projects/qpad_v3/QPAD-v3-report.Rmd
# ---

#NOTE: CONSIDER REMOVING TM OPTION FOR MODELS THAT HAVE A 0 ESTIMATE FOR PC DESPITE METHOD BEING BEST MODEL

library(tidyverse) #basic data wrangling

root <- "G:/Shared drives/BAM_RshProjs/PopnStatus/QPAD"

load(file.path(root, "Results/availability.Rdata"))
load(file.path(root, "Results/perceptability.Rdata"))
load(file.path(root, "Data/qpadv4_clean.Rdata"))
tax <- read.csv(file.path(root, "Data/lookups/taxonomytable.csv"))

#1. Remove species that didn't have enough data----
resDurOK <- avail[!sapply(avail, length)==0]
c(OK=length(resDurOK), failed=length(avail)-length(resDurOK), all=length(avail))

resDisOK <- percep[!sapply(percep, length)==0]
c(OK=length(resDisOK), failed=length(percep)-length(resDisOK), all=length(percep))

#2. Remove species where null model failed---
resDur <- resDurOK[!t(sapply(resDurOK, function(z) ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))[,1]==0]
resDis <- resDisOK[!t(sapply(resDisOK, function(z) ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))[,1]==0]

#3. Create 0/1 table for model fit----
sra_mod <- t(sapply(resDur, function(z)
    ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))
edr_mod <- t(sapply(resDis, function(z)
    ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))

#3. Adjust lists to include all species there's an estimate for----
tmp <- union(rownames(edr_mod), rownames(sra_mod))

sra_models <- matrix(0L, length(tmp), ncol(sra_mod))
dimnames(sra_models) <- list(tmp, colnames(sra_mod))
sra_models[rownames(sra_mod),] <- sra_mod

edr_models <- matrix(0L, length(tmp), ncol(edr_mod))
dimnames(edr_models) <- list(tmp, colnames(edr_mod))
edr_models[rownames(edr_mod),] <- edr_mod

#4. Check for dropped factor levels & min # of detections within each class for removal models----
n.min.class <- 5 # min number of detections within each class
for (i in 1:length(rownames(sra_mod))) {
  spp <- rownames(sra_mod)[i]
  ## data for checking detections in classes
  bird.i <- bird %>%
    dplyr::filter(species==spp,
                  durationMethod %in% durdesign$durationMethod,
                  !durationInterval %in% c("UNKNOWN", "before or after/incidental"),
                  !is.na(TM)) %>%
    dplyr::select(id) %>%
    unique()
  Dat <- visit %>%
    dplyr::filter(id %in% unique(bird.i$id),
                  !is.na(TSSR),
                  !is.na(DSLS),
                  !is.na(JDAY),
                  TM!="ARU-None",
                  !is.na(TM)) %>%
    group_by(TM) %>% 
    summarize(n=n()) %>% 
    ungroup()
  for(j in 1:length(colnames(sra_mod))){
      mid <- colnames(sra_mod)[j]
        if (!inherits(resDur[[spp]][[mid]], "try-error")) {
            lcf <- length(resDur[[spp]][[mid]]$coefficients)
            lnm <- length(resDur[[spp]][[mid]]$names)
            if (lcf != lnm) {
                cat("SRA conflict for", spp, "model", mid,
                    "( len.coef =", lcf, ", len.name =", lnm, ")\n")
                sra_models[spp,mid] <- 0
            }
            else {
              if (mid %in% colnames(sra_mod)[16:30] && min(Dat$n) < n.min.class) {
              cat("Tag method min issue for", spp, "model", mid, "\n")
              sra_models[spp,mid] <- 0
              }
              if (mid %in% colnames(sra_mod)[16] && exp(sum(resDur[[i]][[j]]$coefficients[c(1,2)])) > 2){
                cat("Tag method SPT estimate issue for", spp, "model", mid, "\n")
                sra_models[spp,c(16:30)] <- 0
              }
              if (mid %in% colnames(sra_mod)[16] && exp(sum(resDur[[i]][[j]]$coefficients[c(1,3)])) > 3){
                cat("Tag method SPM estimate issue for", spp, "model", mid, "\n")
                sra_models[spp,c(16:30)] <- 0
              }
        } 
      } else {
            resDur[[spp]][[mid]] <- structure("Error", class = "try-error")
        }
        flush.console()
  }
}

#5. Check for dropped factor levels & min # of detections within each class for detectability models----
n.min.class <- 5 # min number of detections within each class
for (spp in rownames(edr_mod)) {
    ## data for checking detections in classes
    bird.i <- bird %>%
        dplyr::filter(species==spp,
                      distanceMethod %in% distdesign$distanceMethod,
                      distanceBand!="UNKNOWN") %>%
        dplyr::select(id) %>%
        unique()
    Dat <- visit %>%
        dplyr::filter(id %in% unique(bird.i$id),
                      !is.na(TREE),
                      !is.na(LCC2),
                      !is.na(LCC4)) %>%
        arrange(id) %>%
        dplyr::select(LCC2, LCC4)
    for (mid in colnames(edr_models)) {
        if (!inherits(resDis[[spp]][[mid]], "try-error")) {
            lcf <- length(resDis[[spp]][[mid]]$coefficients)
            lnm <- length(resDis[[spp]][[mid]]$names)
            if (lcf != lnm) {
                cat("EDR conflict for", spp, "model", mid,
                    "( len.coef =", lcf, ", len.name =", lnm, ")\n")
                edr_models[spp,mid] <- 0
            } else {
                if (mid %in% c("2", "4") && min(table(Dat$LCC2)) < n.min.class) {
                    cat("EDR LCC2 min issue for", spp, "model", mid, "\n")
                    edr_models[spp,mid] <- 0
                }
                if (mid %in% c("3", "5") && min(table(Dat$LCC4)) < n.min.class) {
                    cat("EDR LCC4 min issue for", spp, "model", mid, "\n")
                    edr_models[spp,mid] <- 0
                }
                if (mid %in% c("1", "4", "5") &&
                    resDis[[spp]][[mid]]$coefficients["log.tau_TREE"] > 0) {
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
sra_models[sra_models[,1]==0,] <- 0L
edr_models[edr_models[,1]==0,] <- 0L

#7. Exclude species with phi estimate < 0.01----
phi0 <- sapply(resDur, function(z) exp(z[["0"]]$coefficients[[1]]))
names(phi0) <- names(resDur)
sra_models[names(phi0)[phi0 < 0.01],] <- 0L

#8. Get number of models----
edr_nmod <- ncol(edr_mod)
sra_nmod <- ncol(sra_mod)

#9. Get sample sizes----
edr_n <- sra_n <- numeric(length(tmp))
names(edr_n) <- names(sra_n) <- tmp
edr_nn <- sapply(resDis, function(z) ifelse(inherits(z[["0"]], "try-error"),
                                            NA, z[["0"]]$nobs))
edr_n[names(edr_nn)] <- edr_nn
sra_nn <- sapply(resDur, function(z) ifelse(inherits(z[["0"]], "try-error"),
                                            NA, z[["0"]]$nobs))
sra_n[names(sra_nn)] <- sra_nn

#10. exclude all models for species with < n.con observations----
#Except NESP because has several points in the great lakes due to the use of the low resolution rasters from the package and so n < 25
n.con <- 25
sra_models[sra_n < n.con, ] <- 0L
edr_n["NESP"] <- 25
edr_models[edr_n < n.con, ] <- 0L

#11. Exclude everything but null for species with n.con < observations < n.min----
n.min <- 75
sra_models[sra_n < n.min & sra_n >= n.con, 2:ncol(sra_models)] <- 0L
edr_models[edr_n < n.min & edr_n >= n.con, 2:ncol(edr_models)] <- 0L

#12. ID species to keep----
spp <- tmp[rowSums(edr_models) > 0 & rowSums(sra_models) > 0]
length(spp)

edr_models <- edr_models[spp,]
sra_models <- sra_models[spp,]
edr_n <- edr_n[spp]
sra_n <- sra_n[spp]

#13. Get number of parameters----
#use OVEN as template
edr_df <- sapply(resDis[["OVEN"]][1:edr_nmod], "[[", "p")
sra_df <- sapply(resDur[["OVEN"]][1:sra_nmod], "[[", "p")

#14. Get estimates----
edr_estimates <- resDis[spp]
sra_estimates <- resDur[spp]

#15. Get species table----
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
edr_list <- sapply(edr_estimates[["OVEN"]], function(z) paste(z$names, collapse=" + "))

#16. Get loglik values---
sra_loglik <- sra_models
sra_loglik[] <- -Inf
for (i in spp) { # species
    for (j in 1:sra_nmod) { # models
        if (sra_models[i,j] > 0)
            sra_loglik[i,j] <- resDur[[i]][[j]]$loglik
    }
}
edr_loglik <- edr_models
edr_loglik[] <- -Inf
for (i in spp) { # species
    for (j in 1:edr_nmod) { # models
        if (edr_models[i,j] > 0)
            edr_loglik[i,j] <- resDis[[i]][[j]]$loglik
    }
}

#17. Get AIC values----
sra_aic <- sra_aicc <- sra_bic <- sra_loglik
sra_aic[] <- Inf
sra_aicc[] <- Inf
sra_bic[] <- Inf
edr_aic <- edr_aicc <- edr_bic <- edr_loglik
edr_aic[] <- Inf
edr_aicc[] <- Inf
edr_bic[] <- Inf
for (i in spp) {
  sra_aic[i,] <- -2*sra_loglik[i,] + 2*sra_df
  sra_aicc[i,] <- sra_aic[i,] + (2*sra_df*(sra_df+1)) / (sra_n[i]-sra_df-1)
  sra_bic[i,] <- -2*sra_loglik[i,] + log(sra_n[i])*sra_df
  edr_aic[i,] <- -2*edr_loglik[i,] + 2*edr_df
  edr_aicc[i,] <- edr_aic[i,] + (2*edr_df*(edr_df+1)) / (edr_n[i]-edr_df-1)
  edr_bic[i,] <- -2*edr_loglik[i,] + log(edr_n[i])*edr_df
}

#18. Rank models----
sra_aicrank <- t(apply(sra_aic, 1, rank))*sra_models
sra_aicrank[sra_aicrank==0] <- NA
edr_aicrank <- t(apply(edr_aic, 1, rank))*edr_models
edr_aicrank[edr_aicrank==0] <- NA

sra_aiccrank <- t(apply(sra_aicc, 1, rank))*sra_models
sra_aiccrank[sra_aiccrank==0] <- NA
edr_aiccrank <- t(apply(edr_aicc, 1, rank))*edr_models
edr_aiccrank[edr_aiccrank==0] <- NA

sra_bicrank <- t(apply(sra_bic, 1, rank))*sra_models
sra_bicrank[sra_bicrank==0] <- NA
edr_bicrank <- t(apply(edr_bic, 1, rank))*edr_models
edr_bicrank[edr_bicrank==0] <- NA

sra_aicbest <- apply(sra_aicrank, 1, function(z) colnames(sra_models)[which.min(z)])
sra_aiccbest <- apply(sra_aiccrank, 1, function(z) colnames(sra_models)[which.min(z)])
sra_bicbest <- apply(sra_bicrank, 1, function(z) colnames(sra_models)[which.min(z)])
edr_aicbest <- apply(edr_aicrank, 1, function(z) colnames(edr_models)[which.min(z)])
edr_aiccbest <- apply(edr_aiccrank, 1, function(z) colnames(edr_models)[which.min(z)])
edr_bicbest <- apply(edr_bicrank, 1, function(z) colnames(edr_models)[which.min(z)])

#19. Replace EDR estimates with V3 for NESP----
library(QPAD)
load_BAM_QPAD(3)

edr_models["NESP",] <- .BAMCOEFS$edr_models["NESP",]
edr_loglik["NESP",] <- .BAMCOEFS$edr_loglik["NESP",]
edr_aic["NESP",] <- .BAMCOEFS$edr_aic["NESP",]
edr_aicc["NESP",] <- .BAMCOEFS$edr_aicc["NESP",]
edr_bic["NESP",] <- .BAMCOEFS$edr_bic["NESP",]
edr_aicrank["NESP",] <- .BAMCOEFS$edr_aicrank["NESP",]
edr_aiccrank["NESP",] <- .BAMCOEFS$edr_aiccrank["NESP",]
edr_bicrank["NESP",] <- .BAMCOEFS$edr_bicrank["NESP",]
edr_aicbest["NESP"] <- .BAMCOEFS$edr_aicbest["NESP"]
edr_aiccbest["NESP"] <- .BAMCOEFS$edr_aiccbest["NESP"]
edr_bicbest["NESP"] <- .BAMCOEFS$edr_bicbest["NESP"]
edr_estimates["NESP"] <- .BAMCOEFS$edr_estimates["NESP"]

#19. Set version----
version <- "4"

#20. Bundle----
bamcoefs <- list(spp=spp,
                 spp_table=spp_table,
                 edr_list=edr_list,
                 sra_list=sra_list,
                 edr_models=edr_models,
                 sra_models=sra_models,
                 edr_n=edr_n,
                 sra_n=sra_n,
                 edr_df=edr_df,
                 sra_df=sra_df,
                 edr_loglik=edr_loglik,
                 sra_loglik=sra_loglik,
                 edr_aic=edr_aic,
                 sra_aic=sra_aic,
                 edr_aicc=edr_aicc,
                 sra_aicc=sra_aicc,
                 edr_bic=edr_bic,
                 sra_bic=sra_bic,
                 edr_aicrank=edr_aicrank,
                 sra_aicrank=sra_aicrank,
                 edr_aiccrank=edr_aiccrank,
                 sra_aiccrank=sra_aiccrank,
                 edr_bicrank=edr_bicrank,
                 sra_bicrank=sra_bicrank,
                 edr_aicbest=edr_aicbest,
                 sra_aicbest=sra_aicbest,
                 edr_aiccbest=edr_aiccbest,
                 sra_aiccbest=sra_aiccbest,
                 edr_bicbest=edr_bicbest,
                 sra_bicbest=sra_bicbest,
                 edr_estimates=edr_estimates,
                 sra_estimates=sra_estimates,
                 version=version)
.BAMCOEFS4 <- list2env(bamcoefs)

#21. Save for internal evaluation----
#for internal evaluation
save(.BAMCOEFS4, file=file.path(root, "Results/BAMCOEFS_QPAD.rda"))

#to QPAD package
saveRDS(bamcoefs, file="C:/Users/elly/Documents/BAM/QPAD/QPAD/inst/estimates/QPAD_v4.rds")
