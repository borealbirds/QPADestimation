library(tidyverse) #basic data wrangling
library(detect) #traditional models
library(data.table) #collapse list to dataframe
library(CmultiJoint.dev) #joint models

options(dplyr.summarise.inform = FALSE, scipen=9999)

root <- "G:/Shared drives/BAM_RshProjs/PopnStatus/QPAD"

#A. ESTIMATES####

#1. Load data----
load(file.path(root, "Data", "qpadv4_clean.Rdata"))

#Set factor levels
visit$TM <- factor(visit$TM, levels=c("PC", "ARU-1SPT", "ARU-1SPM"))

#fix ARU distance band
bird$distanceBand <- ifelse(bird$TM!="PC", "0m-INF", bird$distanceBand)

#fix ARU distance design
visit$distanceMethod <- ifelse(visit$distanceMethod=="0m-INF-ARU", "0m-INF", visit$distanceMethod)
bird$distanceMethod <- ifelse(bird$distanceMethod=="0m-INF-ARU", "0m-INF", bird$distanceMethod)

#collapse 0-3-5-10-10min+ and 0-3-5-10min durationMethods (no 10+ detections in dataset)
visit$durationMethod <- ifelse(visit$durationMethod=="0-3-5-10-10min+", "0-3-5-10min", visit$durationMethod)
bird$durationMethod <- ifelse(bird$durationMethod=="0-3-5-10-10min+", "0-3-5-10min", bird$durationMethod)

#2. Filter bird data----
#no NONE tag method
#no na Tag method
#no unknown methods, bands, intervals
bird.use <- dplyr::filter(bird,
                      TM!="ARU-NONE",
                      !is.na(TM),
                      durationMethod!="UNKNOWN",
                      distanceMethod!="UNKNOWN",
                      durationInterval!="UNKNOWN",
                      distanceBand!="UNKNOWN")

#3. Create design lookup tables----
durdesign <- visit %>%
  dplyr::select(durationMethod) %>%
  unique() %>%
  dplyr::filter(!durationMethod %in% c("UNKNOWN"),
                !is.na(durationMethod)) %>%
  mutate(dm = str_sub(durationMethod, -100, -4)) %>%
  separate(dm, into=c("t00", "t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), remove=TRUE, sep="-", extra="drop", fill="right") %>%
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
  separate(distanceMethod, into=c("d00", "d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), remove=FALSE, sep="-", extra="drop", fill="right") %>%
  dplyr::select(-d00) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~str_sub(., -100, -2)) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10", "d11", "d12", "d13"), ~ifelse(.=="IN", Inf, .)) %>%
  mutate_at(c("d01", "d02", "d03", "d04", "d05", "d06", "d07", "d08", "d09", "d10",  "d11", "d12", "d13"), ~as.numeric(.)/100)

distdesign.long <- distdesign %>% 
  pivot_longer(d01:d13, values_to="far", names_to="distposition") %>% 
  dplyr::filter(!is.na(far))

#4. Get list of species to process----
# spp <- species %>%
#       dplyr::filter(Singing_birds==TRUE) %>%
#   left_join(bird %>%
#               dplyr::select(species) %>%
#               unique()) %>%
#   arrange(species)
spp <- data.frame(species=c("WETA", "AMGO", "PAWA", "LEFL", "ALFL"))

#5. Set up loop for species----
#joint <- list()
for(i in 1:nrow(spp)){
  
  start <- Sys.time()
  
  #6. Filter abundance data for species---
  # filter out observations with unknown duration method or interval
  # filter out observations where interval does not match method
  # sum abundance for each bin in each visit
  bird.i <- bird.use %>%
    dplyr::filter(species==spp$species[i],
                  durationMethod %in% durdesign$durationMethod,
                  distanceMethod %in% distdesign$distanceMethod,
                  !durationInterval %in% c("UNKNOWN", "before or after/incidental")) %>% 
    separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE, extra="drop", fill="right") %>%
    mutate(start = as.numeric(start),
           end = as.numeric(str_sub(end, -100, -4))) %>% 
    left_join(durdesign.long, by=c("durationMethod", "end")) %>% 
    separate(distanceBand, into=c("near", "far"), sep="-", remove=FALSE, extra="drop", fill="right") %>% 
    mutate(near = as.numeric(str_sub(near, -100, -2)),
           far = ifelse(far=="INF", Inf, suppressWarnings(as.numeric(str_sub(far, -100, -2))/100))) %>% 
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
    dimnames(tarray) <- NULL
    
    rarray <- x %>%
      dplyr::select(distanceMethod) %>%
      left_join(distdesign, by="distanceMethod") %>%
      dplyr::select(-distanceMethod) %>%
      as.matrix()
    dimnames(rarray) <- NULL
    
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
    
    #11. Create covariate matrices----
    X1.TM <- model.matrix(~x$TM)
    colnames(X1.TM) <- c("tau_int","tau_SPT", "tau_SPM")
    
    X2.TM <- model.matrix(~x$TM)
    colnames(X2.TM) <- c("phi_int","phi_SPT", "phi_SPM")
    
    X1 <- list(NULL, X1.TM, NULL, X1.TM)
    X2 <- list(NULL, NULL, X2.TM, X2.TM)
    
    names <- list(c("tau_int", "phi_int"),
                  c("tau_int", "tau_SPT", "tau_SPM", "phi_int"),
                  c("tau_int", "phi_int", "phi_SPT", "phi_SPM"),
                  c("tau_int", "tau_SPT", "tau_SPM", "phi_int", "phi_SPT", "phi_SPM"))
    
    #Save a bunch of metadata like sample size and aic value
    fit <- list()
    for(j in 1:4){
      fit[[j]] <- try(cmulti_fit_joint(Yarray,
                                       rarray,
                                       tarray,
                                       X1 = X1[[j]],
                                       X2 = X2[[j]]))
    }
    
    mod.list <- list()
    for(j in 1:length(fit)){
      if (!inherits(fit[[j]], "try-error")) {
        mod.list[[j]] <- fit[[j]][c("coefficients","vcov","loglik")]
        mod.list[[j]]$nobs <- nrow(x)
        mod.list[[j]]$p <- length(coef(fit[[j]]))
        mod.list[[j]]$names <- names[[j]]
      } else {
        mod.list[[j]] <- fit[[j]]
      }
    }
    names(mod.list) <- c("tau-null_phi-null", "tau-TM_phi-null", "tau-null_phi-TM", "tau-TM_phi-TM")
    
    #13. Save model results to species list---
    joint[[i]] <- mod.list
    
    names(joint) <- spp$species[1:i]
  
  }
  
  #14. Save results to local----
  
  save(joint, file=file.path(root, "Results/jointestimates.Rdata")) 
  
  end <- Sys.time()
  
  print(paste0("Finished modelling species ", spp$English_Name[i], ": ", i, " of ", nrow(spp), " species in ", round(end-start, 2), " minutes"))
}

#B. PACKAGE####

load(file.path(root, "Results/jointestimates.Rdata"))

#1. Remove species that didn't have enough data----
resOK <- joint[!sapply(joint, length)==0]
c(OK=length(resOK), failed=length(joint)-length(resOK), all=length(joint))

#2. Remove species where null model failed---
res <- resOK[!t(sapply(resOK, function(z) ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))[,1]==0]

#3. Create 0/1 table for model fit----
mod <- t(sapply(res, function(z)
  ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))

models <- matrix(0L, nrow(mod), ncol(mod))
dimnames(models) <- list(rownames(mod), colnames(mod))
models[rownames(mod),] <- mod

#4. Check for dropped factor levels & min # of detections within each class for removal models----
n.min.class <- 5 # min number of detections within each class
for (i in 1:length(rownames(mod))) {
  spp <- rownames(mod)[i]
  ## data for checking detections in classes
  bird.i <- bird.use %>%
    dplyr::filter(species==spp,
                  durationMethod %in% durdesign$durationMethod,
                  distanceMethod %in% distdesign$distanceMethod,
                  !durationInterval %in% c("UNKNOWN", "before or after/incidental"),
                  !is.na(TM)) %>% 
    dplyr::select(id) %>%
    unique()
  Dat <- visit %>%
    dplyr::filter(id %in% unique(bird.i$id),
                  TM!="ARU-None",
                  !is.na(TM)) %>%
    group_by(TM) %>% 
    summarize(n=n()) %>% 
    ungroup()
  for(j in 1:length(colnames(mod))){
    mid <- colnames(mod)[j]
    if (!inherits(res[[spp]][[mid]], "try-error")) {
      lcf <- length(res[[spp]][[mid]]$coefficients)
      lnm <- length(res[[spp]][[mid]]$names)
      if (lcf != lnm) {
        cat("conflict for", spp, "model", mid,
            "( len.coef =", lcf, ", len.name =", lnm, ")\n")
        models[spp,mid] <- 0
      }
      else {
        if (mid %in% colnames(mod)[2:4] && min(Dat$n) < n.min.class) {
          cat("Tag method min issue for", spp, "model", mid, "\n")
          models[spp,mid] <- 0
        }
      } 
    } else {
      res[[spp]][[mid]] <- structure("Error", class = "try-error")
    }
    flush.console()
  }
}

#6. Exclude species with no null model----
models[models[,1]==0,] <- 0L

#7. Exclude species with phi estimate < 0.01----
phi0 <- sapply(res, function(z) exp(z[["tau-null_phi-null"]]$coefficients[[2]]))
names(phi0) <- names(res)
models[names(phi0)[phi0 < 0.01],] <- 0L

#8. Get number of models----
nmod <- ncol(mod)

#9. Get sample sizes----
n <- numeric(length(models))
names(n) <- rownames(models)
nn <- sapply(res, function(z) ifelse(inherits(z[["tau-null_phi-null"]], "try-error"),
                                            NA, z[["tau-null_phi-null"]]$nobs))

#10. exclude all models for species with < n.con observations----
#Except NESP because has several points in the great lakes due to the use of the low resolution rasters from the package and so n < 25
n.con <- 25
models[nn < n.con, ] <- 0L

#11. Exclude everything but null for species with n.con < observations < n.min----
n.min <- 75
models[nn < n.min & nn >= n.con, 2:ncol(models)] <- 0L

#12. ID species to keep----
spp <- rownames(models)[rowSums(models) > 0]
length(spp)

models <- models[spp,]
n <- n[spp]

#13. Get number of parameters----
#use OVEN as template
df <- sapply(res[["OVEN"]][1:nmod], "[[", "p")

#14. Get estimates----
estimates <- res[spp]

#15. Get species table----
tax <- read.csv(file.path(root, "Data/lookups/taxonomytable.csv"))
tax <- tax[!duplicated(tax$Species_ID),]
rownames(tax) <- tax$Species_ID
spp_table <- data.frame(spp=spp,
                        scientific_name=tax[spp, "Scientific_Name"],
                        common_name=tax[spp, "English_Name"])
rownames(spp_table) <- spp
spp_table <- droplevels(spp_table)

#16. Get variable names for different models----
#use OVEN as template
mod_list <- sapply(estimates[["OVEN"]], function(z) paste(z$names, collapse=" + "))

#16. Get loglik values---
loglik <- models
loglik[] <- -Inf
for (i in spp) { # species
  for (j in 1:nmod) { # models
    if (models[i,j] > 0)
      loglik[i,j] <- res[[i]][[j]]$loglik
  }
}

#17. Get AIC values----
aic <- aicc <- bic <- loglik
aic[] <- Inf
aicc[] <- Inf
bic[] <- Inf
for (i in spp) {
  aic[i,] <- -2*loglik[i,] + 2*df
  aicc[i,] <- aic[i,] + (2*df*(df+1)) / (n[i]-df-1)
  bic[i,] <- -2*loglik[i,] + log(n[i])*df
}

#18. Rank models----
aicrank <- t(apply(aic, 1, rank))*models
aicrank[aicrank==0] <- NA

aiccrank <- t(apply(aicc, 1, rank))*models
aiccrank[aiccrank==0] <- NA

bicrank <- t(apply(bic, 1, rank))*models
bicrank[bicrank==0] <- NA

aicbest <- apply(aicrank, 1, function(z) colnames(models)[which.min(z)])
aiccbest <- apply(aiccrank, 1, function(z) colnames(models)[which.min(z)])
bicbest <- apply(bicrank, 1, function(z) colnames(models)[which.min(z)])

#19. Set version----
version <- "joint"

#20. Bundle----
bamcoefs <- list(spp=spp,
                 spp_table=spp_table,
                 list=list,
                 models=models,
                 n=n,
                 df=df,
                 loglik=loglik,
                 aic=aic,
                 aicc=aicc,
                 bic=bic,
                 aicrank=aicrank,
                 aiccrank=aiccrank,
                 bicrank=bicrank,
                 aicbest=aicbest,
                 aiccbest=aiccbest,
                 bicbest=bicbest,
                 estimates=estimates,
                 version=version)
.BAMCOEFSjoint <- list2env(bamcoefs)

#21. Save for internal evaluation----
#for internal evaluation
save(.BAMCOEFSjoint, file=file.path(root, "Results/BAMCOEFS_QPAD_joint.rda"))
saveRDS(bamcoefs, file=file.path(root, "Results/BAMCOEFS_QPAD_joint.rds"))

#C. COMPARE####

.BAMCOEFSjoint <- readRDS(file=file.path(root, "Results/BAMCOEFS_QPAD_joint.rds"))
spp <- .BAMCOEFSjoint$spp

#1. Read in QPADV4 nulll estimates----
est34 <- read.csv(file.path(root, "Results/QPADV3V4NullEstimates.csv"))

#4. Null estimates----
estj <-data.frame()
for(i in 1:length(spp)){
  sraj <- exp(.BAMCOEFSjoint$estimates[[i]]$"tau-null_phi-null"$coefficients[2])
  edrj <- exp(.BAMCOEFSjoint$estimates[[i]]$"tau-null_phi-null"$coefficients[1])
  estj <- rbind(estj,
                data.frame(sraj=sraj, edrj=edrj, species=spp[i]))
}

estnull <- est34 %>% 
  dplyr::select(species, sra4, edr4) %>% 
  left_join(estj) %>% 
  mutate(offj = )

#Null estimates----
ggplot(estnull) +
  geom_point(aes(x=sra4, y=sraj, label=species)) +
  geom_smooth(aes(x=sra4, y=sraj), method="lm") +
  geom_abline(intercept = 0, slope = 1)+
  xlab("Phi - QPAD V4") +
  ylab("Phi - joint model") +
  theme_bw()

ggsave(file.path(root, "Figures", "JointModelComparison - phi - null.jpeg"), width = 6, height = 4, units="in")

t.test(estnull$sra4, estnull$sraj, paired=TRUE)

mod1 <- lm(sraj~ sra4, data=estnull)
summary(mod1)

ggplot(estnull) +
  geom_point(aes(x=edr4, edrj)) +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Tau - QPAD V4") +
  ylab("Tau - joint model") +
  theme_bw()

ggsave(file.path(root, "Figures", "JointModelComparison - EDR - null.jpeg"), width = 6, height = 4, units="in")

#5. TM estimates----
estj <-data.frame(edrj.tauTM.PC=rep(NA, length(spp)), edrj.tauTM.SPT=rep(NA, length(spp)), edrj.tauTM.SPM=rep(NA, length(spp)), sraj.tauTM.all=rep(NA, length(spp)), edrj.phiTM.all=rep(NA, length(spp)), sraj.phiTM.PC=rep(NA, length(spp)), sraj.phiTM.SPT=rep(NA, length(spp)), sraj.phiTM.SPM=rep(NA, length(spp)), edrj.bothTM.PC=rep(NA, length(spp)), edrj.bothTM.SPT=rep(NA, length(spp)), edrj.bothTM.SPM=rep(NA, length(spp)), sraj.bothTM.PC=rep(NA, length(spp)), sraj.bothTM.SPT=rep(NA, length(spp)), sraj.bothTM.SPM=rep(NA, length(spp)))
for(i in 1:length(spp)){
  
  estj$edrj.tauTM.PC[i] <- exp(.BAMCOEFSjoint$estimates[[i]]$"tau-TM_phi-null"$coefficients[c(1),])
  estj$edrj.tauTM.SPT[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-TM_phi-null"$coefficients[c(1,2),]))
  estj$edrj.tauTM.SPM[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-TM_phi-null"$coefficients[c(1,3),]))
  estj$sraj.tauTM.all[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-TM_phi-null"$coefficients[c(4),]))
  
  estj$edrj.phiTM.all[i] <- exp(.BAMCOEFSjoint$estimates[[i]]$"tau-null_phi-TM"$coefficients[1,])
  estj$sraj.phiTM.PC[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-null_phi-TM"$coefficients[2,]))
  estj$sraj.phiTM.SPT[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-null_phi-TM"$coefficients[c(2,3),]))
  estj$sraj.phiTM.SPM[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-null_phi-TM"$coefficients[c(2,4),]))
  
  estj$edrj.bothTM.PC[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-TM_phi-TM"$coefficients[c(1),]))
  estj$edrj.bothTM.SPT[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-TM_phi-TM"$coefficients[c(1,2),]))
  estj$edrj.bothTM.SPM[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-TM_phi-TM"$coefficients[c(1,3),]))
  estj$sraj.bothTM.PC[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-TM_phi-TM"$coefficients[c(4),]))
  estj$sraj.bothTM.SPT[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-TM_phi-TM"$coefficients[c(4,5),]))
  estj$sraj.bothTM.SPM[i] <- exp(sum(.BAMCOEFSjoint$estimates[[i]]$"tau-TM_phi-TM"$coefficients[c(4,6),]))
  
  estj$species[i] <- spp[i]

}

estTM <- est34 %>% 
  dplyr::select(species, edr4, sra4, sra4.pc, sra4.spt, sra4.spm) %>% 
  left_join(estj) %>% 
  left_join(estnull)

#PHI - vsQPAD4 - with method on phi side----
ggplot(estTM) +
  geom_point(aes(x=sra4.pc, sraj.phiTM.PC)) +
  geom_abline(intercept = 0, slope = 1)

ggplot(estTM) +
  geom_point(aes(x=sra4.spt, sraj.phiTM.SPT)) +
  geom_abline(intercept = 0, slope = 1) +
  ylim(c(0,5))

ggplot(estTM) +
  geom_point(aes(x=sra4.spm, sraj.phiTM.SPM)) +
  geom_abline(intercept = 0, slope = 1) +
  ylim(c(0,5))

#PHI - vs QPAD4 - with method on both sides----
ggplot(estTM) +
  geom_point(aes(x=sra4.pc, sraj.bothTM.PC)) +
  geom_abline(intercept = 0, slope = 1)

ggplot(estTM) +
  geom_point(aes(x=sra4.spt, sraj.bothTM.SPT)) +
  geom_abline(intercept = 0, slope = 1) +
  ylim(c(0,5))

ggplot(estTM) +
  geom_point(aes(x=sra4.spm, sraj.bothTM.SPM)) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0,5)) +
  ylim(c(0,5))

#PHI - method on phi side vs method on both sides----
ggplot(estTM) +
  geom_point(aes(x=sraj.phiTM.PC, sraj.bothTM.PC)) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0,5)) +
  ylim(c(0,5)) +
  xlab("Phi - point count - joint model - method on phi side") +
  ylab("Phi - point count - joint model - method on both sides") +
  theme_bw()

ggsave(file.path(root, "Figures", "JointModelComparison - PCPhi - MethodModels.jpeg"), width=6, height = 4, units="in")

ggplot(estTM) +
  geom_point(aes(x=sraj, sraj.bothTM.PC)) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0,5)) +
  ylim(c(0,5)) +
  xlab("Phi - point count - joint model - null model") +
  ylab("Phi - point count - joint model - method on both sides") +
  theme_bw()

ggsave(file.path(root, "Figures", "JointModelComparison - PCPhi - NullVMethodModels.jpeg"), width=6, height = 4, units="in")

ggplot(estTM) +
  geom_point(aes(x=sraj.phiTM.SPT, sraj.bothTM.SPT)) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0,5)) +
  ylim(c(0,5))

ggplot(estTM) +
  geom_point(aes(x=sraj.phiTM.SPM, sraj.bothTM.SPM)) +
  geom_abline(intercept = 0, slope = 1) +
  xlim(c(0,5)) +
  ylim(c(0,5))

#PHI point count vs null----
ggplot(estTM) +
  geom_point(aes(x=sraj, y=sraj.phiTM.PC)) + 
  geom_abline(intercept = 0, slope = 1)

ggplot(estTM) +
  geom_point(aes(x=sraj, y=sraj.bothTM.PC)) + 
  geom_abline(intercept = 0, slope = 1)
  
#EDR - vsQPAD4 - with method on EDR side----
ggplot(estTM) +
  geom_point(aes(x=edr4, edrj.tauTM.PC)) +
  geom_abline(intercept = 0, slope = 1)

ggplot(estTM) +
  geom_point(aes(x=edr4, edrj.tauTM.SPM)) +
  geom_abline(intercept = 0, slope = 1)

ggplot(estTM) +
  geom_point(aes(x=edr4, edrj.tauTM.SPT)) +
  geom_abline(intercept = 0, slope = 1)

#EDR - vs QPAD4 - with method on both sides----
ggplot(estTM) +
  geom_point(aes(x=edr4, edrj.bothTM.PC)) +
  geom_abline(intercept = 0, slope = 1)

ggplot(estTM) +
  geom_point(aes(x=edr4, edrj.bothTM.SPM)) +
  geom_abline(intercept = 0, slope = 1)

ggplot(estTM) +
  geom_point(aes(x=edr4, edrj.bothTM.SPT)) +
  geom_abline(intercept = 0, slope = 1)

#EDR - method on EDR side vs method on both sides----
ggplot(estTM) +
  geom_point(aes(x=edrj.tauTM.PC, edrj.bothTM.PC)) +
  geom_abline(intercept = 0, slope = 1)

ggplot(estTM) +
  geom_point(aes(x=edrj.tauTM.SPT, edrj.bothTM.SPT)) +
  geom_abline(intercept = 0, slope = 1)

ggplot(estTM) +
  geom_point(aes(x=edrj.tauTM.SPM, edrj.bothTM.SPM)) +
  geom_abline(intercept = 0, slope = 1)

#6. Compare method model formulations----

#PHI----
estj.phi <- estTM %>% 
  dplyr::select(species, sraj.tauTM.all, sraj.phiTM.PC, sraj.phiTM.SPT, sraj.phiTM.SPM, sraj.bothTM.PC, sraj.bothTM.SPT, sraj.bothTM.SPM) %>% 
  full_join(estnull %>% 
              dplyr::select(species, sraj) %>% 
              rename(sraj.null.all = sraj)) %>% 
  pivot_longer(cols=sraj.tauTM.all:sraj.null.all, names_to="estimate", values_to="phi") %>% 
  separate(estimate, into=c("sra", "model", "method")) %>% 
  dplyr::select(-sra) %>% 
  dplyr::filter(!is.na(phi))
  
ggplot(estj.phi) +
  geom_violin(aes(x=method, y=log(phi), colour=model), position = position_dodge(width = 0.7)) +
  xlab("Survey method type") +
  ylab("log(phi)")

ggsave(file.path(root, "Figures", "JointModelComparison - Method.jpeg"), width = 6, height = 4, units="in")

#TAU----
estj.tau <- estTM %>% 
  dplyr::select(species, edrj.phiTM.all, edrj.tauTM.PC, edrj.tauTM.SPT, edrj.tauTM.SPM, edrj.bothTM.PC, edrj.bothTM.SPT, edrj.bothTM.SPM) %>% 
  full_join(estnull %>% 
              dplyr::select(species, edrj) %>% 
              rename(edrj.null.all = edrj)) %>% 
  pivot_longer(cols=edrj.phiTM.all:edrj.null.all, names_to="estimate", values_to="tau") %>% 
  separate(estimate, into=c("edr", "model", "method")) %>% 
  dplyr::select(-edr) %>% 
  dplyr::filter(!is.na(tau))

ggplot(estj.tau) +
  geom_violin(aes(x=method, y=log(tau), colour=model))
