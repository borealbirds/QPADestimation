library(tidyverse) #basic data wrangling
library(detect) #traditional models
library(data.table) #collapse list to dataframe
library(CmultiJoint.dev) #joint models

#A. PREAMBLE####

#1. System options----
options(dplyr.summarise.inform = FALSE, scipen=9999)

#2. Working google drive path----
root <- "G:/Shared drives/BAM_RshProjs/PopnStatus/QPAD"

#3. Load model output----
load(file.path(root, "Results/jointmodels.Rdata"))

#4. Fix offset calculation function bug----
calculate_offsets <- function (fit, rarray = rarray, tarray = tarray, X1 = NULL, 
                               X2 = NULL) {
  rarray_fit = fit$input_data$rarray
  tarray_fit = fit$input_data$tarray
  X1_fit = fit$input_data$X1
  X2_fit = fit$input_data$X2
  maxdistint = fit$input_data$maxdistint
  nsurvey <- dim(rarray)[1]
  nrint <- apply(rarray, 1, function(x) length(na.omit(x)))
  ntint <- apply(tarray, 1, function(x) length(na.omit(x)))
  if (!is.null(X1)) {
    tau_params <- colnames(X1)
  }
  else {
    X1 <- matrix(1, nrow = nsurvey, ncol = 1)
    tau_params <- colnames(X1)[1] <- "log_tau"
  }
  if (!is.null(X2)) {
    phi_params <- colnames(X2)
  }
  else {
    X2 <- matrix(1, nrow = nsurvey, ncol = 1)
    phi_params <- colnames(X2)[1] <- "log_phi"
  }
  max_r <- apply(rarray, 1, max, na.rm = TRUE)
  max_r[max_r == Inf] <- maxdistint
  tau_params <- fit$coefficients[1:length(tau_params)]
  phi_params <- fit$coefficients[(length(tau_params) + 1):length(fit$coefficients)]
  tau <- poisson("log")$linkinv(drop(X1 %*% tau_params))
  phi <- poisson("log")$linkinv(drop(X2 %*% phi_params))
  p <- A <- rep(NA, nsurvey)
  for (k in 1:nsurvey) {
    tau_k <- tau[k]
    phi_k <- phi[k]
    f_d = function(dmax) {
      integrand = substitute(2 * pi * dmax * (1 - exp(-phi * 
                                                        tmax * exp(-dmax^2/tau^2))), list(phi = phi_k, 
                                                                                          tau = tau_k, tmax = tmax))
      eval(integrand)
    }
    #    Y <- Yarray[k, 1:nrint[k], 1:ntint[k]]
    CDF_binned <- matrix(NA, nrow = nrint[k], ncol = ntint[k])
    for (j in 1:ntint[k]) {
      tmax = max(tarray[k, j])
      for (i in 1:nrint[k]) {
        upper_r = rarray[k, i]
        if (upper_r == Inf) 
          upper_r = max_r[k]
        CDF_binned[i, j] = integrate(f_d, lower = 0.01, 
                                     upper = upper_r, subdivisions = 500)$value
      }
    }
    tmp1 = CDF_binned
    if (nrow(tmp1) > 1) {
      for (i in 2:nrint[k]) {
        tmp1[i, ] <- CDF_binned[i, ] - CDF_binned[i - 
                                                    1, ]
      }
    }
    p_matrix = tmp1
    if (ncol(p_matrix) > 1) {
      for (j in 2:ntint[k]) {
        p_matrix[, j] <- tmp1[, j] - tmp1[, j - 1]
      }
    }
    p[k] <- sum(p_matrix)
    A[k] <- pi * max_r[k]^2
  }
  log_offset <- log(p)
  log_offset
}

#B. NULL ESTIMATES####

#1. Make arrays for offset calculation----
rarray <- array(Inf,dim=c(1,1))
tarray <- array(10,dim=c(1,1))

#2. Get estimates----
estj <- data.frame()
for(i in 1:length(mod.list)){
  
  offj <- calculate_offsets(mod.list[[i]][[1]],
                           rarray = rarray,
                           tarray = tarray)
  
  countj <- exp(offj) # Expected count if density is 1
  
  phij <- exp(mod.list[[i]][[1]]$coefficients[2])
  
  tauj <- exp(mod.list[[i]][[1]]$coefficients[1])

  estj <- rbind(estj, data.frame(offj=offj, countj=countj, phij=phij, tauj=tauj, species = names(mod.list)[i]))
}

#3. Read in QPADV3 & V4 estimates----
#Wrangle
#Calculate density est
est34 <- read.csv(file.path(root, "Results/QPADV3V4NullEstimates.csv")) %>% 
  rename(phi3=sra3, tau3=edr3, phi4=sra4, tau4=edr4) %>% 
  dplyr::select(species, phi3, tau3, phi4, tau4) %>% 
  mutate(count3 = (pi*tau3^2)*(1-exp(-10*phi3)),
         count4 = (pi*tau4^2)*(1-exp(-10*phi4)))

#4. Put together with joint estimates----
est <- inner_join(estj, est34) %>% 
  mutate(diff3 = 100*(countj/count3-1),
         diff4 = 100*(countj/count4-1)) %>%
  arrange(diff4)

#rearrange by percent diff
est$species <- factor(est$species, levels=est$species)

#5. Plot----

#5a. Phi----
ggplot(est) +
  geom_abline(intercept=0, slope=1) +
  geom_text(aes(x=phi3, y=phi4, label=species)) +
  theme_bw()

ggplot(est) +
  geom_abline(intercept=0, slope=1) +
  geom_text(aes(x=phi3, y=phij, label=species)) +
  theme_bw()

plot.phi <- ggplot(est %>% dplyr::filter(species!="LCSP")) +
    geom_abline(intercept=0, slope=1) +
    geom_text(aes(x=phi4, y=phij, label=species)) +
    theme_bw() +
  xlab("Phi - QPADV4") +
  ylab("Phi - joint estimation")
plot.phi

#5b. Tau----
ggplot(est) +
  geom_abline(intercept=0, slope=1) +
  geom_text(aes(x=tau3, y=tau4, label=species)) +
  theme_bw()

ggplot(est) +
  geom_abline(intercept=0, slope=1) +
  geom_text(aes(x=tau3, y=tauj, label=species)) +
  theme_bw()

plot.tau <- ggplot(est) +
  geom_abline(intercept=0, slope=1) +
  geom_text(aes(x=tau4, y=tauj, label=species)) +
  theme_bw() +
  xlab("Tau - QPADV4") +
  ylab("Tau - joint estimation")
plot.tau

#5c. Offset----
ggplot(est) +
  geom_abline(intercept=0, slope=1) +
  geom_text(aes(x=count3, y=count4, label=species)) +
  theme_bw()

ggplot(est) +
  geom_abline(intercept=0, slope=1) +
  geom_text(aes(x=count3, y=countj, label=species)) +
  theme_bw()

plot.off <- ggplot(est) +
  geom_abline(intercept=0, slope=1) +
  geom_text(aes(x=count4, y=countj, label=species)) +
  theme_bw() +
  xlab("Offset - QPADV4") +
  ylab("Offset - joint estimation")
plot.off

#5d. Population estimation difference----
plot.diff <- ggplot(arrange(est, diff4)) +
  geom_bar(aes(x = species, y = diff4), stat = "identity", fill = "grey70")+
  xlab("Species")+
  ylab("Percent difference in density estimate") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
plot.diff

ggsave(gridExtra::grid.arrange(plot.phi, plot.tau, plot.off, plot.diff, ncol=2, nrow=2), filename=file.path(root, "Figures", "JointModelComparison.jpeg"), width=15, height=15)

ggsave(gridExtra::grid.arrange(plot.phi, plot.tau, ncol=2, nrow=1), filename=file.path(root, "Figures", "JointModelComparison - phi&tau.jpeg"), width=15, height=8)

ggsave(gridExtra::grid.arrange(plot.off, plot.diff, ncol=2, nrow=1), filename=file.path(root, "Figures", "JointModelComparison - off&diff.jpeg"), width=15, height=8)

#B. SENSOR ESTIMATES####

#1. Make arrays for offset calculation----
rarray <- array(Inf,dim=c(1,1))
tarray <- array(10,dim=c(1,1))
x <- data.frame(TM=factor(c("PC", "ARU-1SPT", "ARU-1SPM"), levels=c("PC", "ARU-1SPT", "ARU-1SPM")))
XPC <- model.matrix(~x$TM[1])
XSPT <- model.matrix(~x$TM[2])
XSPM <- model.matrix(~x$TM[3])
X <- list(XPC, XSPT, XSPM)
names(X) <- c("Point count", "ARU - 1SPT", "ARU - 1SPM")

#2. Get estimates----
estj <- data.frame()

#2a. tau-TM_phi-null----
for(i in 1:length(mod.list)){
  
  for(j in 1:length(X)){
    
    offj <- calculate_offsets(mod.list[[i]][[2]],
                              rarray = rarray,
                              tarray = tarray,
                              X1=X[[j]])
    
    countj <- exp(offj) # Expected count if density is 1
    
    phij <- exp(mod.list[[i]][[2]]$coefficients[4])
    
    if(j==1){
      tauj <- exp(mod.list[[i]][[2]]$coefficients[1])
    }
    if(j>1){
      tauj <- exp(sum(mod.list[[i]][[2]]$coefficients[c(1,j)]))
    }

    
    estj <- rbind(estj, data.frame(offj=offj, countj=countj, phij=phij, tauj=tauj, species = names(mod.list)[i], phisensor="null", tausensor = names(X)[[j]], phimodel="null", taumodel="TM"))
  }

}

#2b. tau-null_phi-TM----
for(i in 1:length(mod.list)){
  
  for(j in 1:length(X)){
    
    offj <- calculate_offsets(mod.list[[i]][[3]],
                              rarray = rarray,
                              tarray = tarray,
                              X2=X[[j]])
    
    countj <- exp(offj) # Expected count if density is 1
    
    tauj <- exp(mod.list[[i]][[3]]$coefficients[1])
    
    if(j==1){
      phij <- exp(mod.list[[i]][[3]]$coefficients[2])
    }
    if(j>1){
      phij <- exp(sum(mod.list[[i]][[3]]$coefficients[c(2,j+1)]))
    }
    
    
    estj <- rbind(estj, data.frame(offj=offj, countj=countj, phij=phij, tauj=tauj, species = names(mod.list)[i], phisensor = names(X)[[j]], tausensor="null", phimodel="TM", taumodel="null"))
  }
  
}

#2c. tau-TM_phi-TM----
for(i in 1:length(mod.list)){
  
  for(j in 1:length(X)){
    
    offj <- calculate_offsets(mod.list[[i]][[4]],
                              rarray = rarray,
                              tarray = tarray,
                              X1=X[[j]],
                              X2=X[[j]])
    
    countj <- exp(offj) # Expected count if density is 1
    
    if(j==1){
      tauj <- exp(mod.list[[i]][[4]]$coefficients[1])
      phij <- exp(mod.list[[i]][[4]]$coefficients[4])
    }
    if(j>1){
      tauj <- exp(sum(mod.list[[i]][[4]]$coefficients[c(1,j)]))
      phij <- exp(sum(mod.list[[i]][[4]]$coefficients[c(4,j+3)]))
    }
    
    
    estj <- rbind(estj, data.frame(offj=offj, countj=countj, phij=phij, tauj=tauj, species = names(mod.list)[i], phisensor = names(X)[[j]], tausensor = names(X)[[j]], phimodel="TM", taumodel="TM"))
  }
  
}

#3. Read in QPADV3 & V4 estimates----
#Wrangle
#Calculate density est
est34 <- read.csv(file.path(root, "Results/QPADV3V4NullEstimates.csv")) %>% 
  dplyr::select(species, sra4, sra4.pc, sra4.spt, sra4.spm, edr4) %>% 
  pivot_longer(cols=sra4:sra4.spm, names_to="phisensor", values_to="phi4") %>% 
  rename(tau4=edr4) %>% 
  mutate(phimodel = ifelse(phisensor=="sra4", "null", "TM"),
         phisensor = case_when(phisensor=="sra4" ~ "null",
                               phisensor=="sra4.pc" ~ "Point count",
                               phisensor=="sra4.spt" ~ "ARU - 1SPT",
                               phisensor=="sra4.spm" ~ "ARU - 1SPM"),
         tausensor = "null",
         taumodel = "null",
         count4 = (pi*tau4^2)*(1-exp(-10*phi4))) %>% 
  dplyr::filter(species!="CCLO")

#4. Put together with joint estimates----
est <- left_join(estj, est34) %>% 
  mutate(diff = 100*(countj/count4-1),
         sensor = ifelse(tausensor=="null", phisensor, tausensor)) %>% 
  arrange(phisensor, diff) 

#rearrange by percent diff
est$species <- factor(est$species, levels=unique(est$species))

#5. Plot V4 vs joint----

#5a. Phi----
plot.phi <- ggplot(est %>% dplyr::filter(phimodel=="TM", taumodel=="null", species!="CCLO")) +
  geom_abline(intercept=0, slope=1) +
  geom_text(aes(x=phi4, y=phij, label=species)) +
  theme_bw() +
  facet_wrap(~phisensor, scales="free") +
  xlab("Phi - QPADV4") +
  ylab("Phi - joint estimation")
plot.phi

ggsave(filename=file.path(root, "Figures", "JointModelComparison - PhiMethod - Phi.jpeg"), width=15, height=8)

#5b. Tau----
plot.tau <- ggplot(est %>% dplyr::filter(phimodel=="TM", taumodel=="null", species!="CCLO")) +
  geom_abline(intercept=0, slope=1) +
  geom_text(aes(x=tau4, y=tauj, label=species)) +
  theme_bw() +
  facet_wrap(~phisensor, scales="free") +
  xlab("Tau - QPADV4") +
  ylab("Tau - joint estimation")
plot.tau

ggsave(filename=file.path(root, "Figures", "JointModelComparison - PhiMethod - Tau.jpeg"), width=15, height=8)


#5c. Offset----
plot.off <- ggplot(est %>% dplyr::filter(phimodel=="TM", taumodel=="null", species!="CCLO")) +
  geom_abline(intercept=0, slope=1) +
  geom_text(aes(x=count4, y=countj, label=species)) +
  theme_bw() +
  facet_wrap(~phisensor, scales="free") +
  xlab("Offset - QPADV4") +
  ylab("Offset - joint estimation")
plot.off

ggsave(filename=file.path(root, "Figures", "JointModelComparison - PhiMethod - offset.jpeg"), width=15, height=8)

#5d. Percent difference----
plot.diff <- ggplot(est %>% dplyr::filter(phimodel=="TM", taumodel=="null", species!="CCLO")) +
  geom_bar(aes(x = species, y = diff), stat = "identity", fill = "grey70")+
  xlab("Species")+
  ylab("Percent difference in density estimate") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
  facet_wrap(~phisensor, scales="free")
plot.diff

ggsave(filename=file.path(root, "Figures", "JointModelComparison - PhiMethod - difference.jpeg"), width=15, height=8)

ggsave(gridExtra::grid.arrange(plot.phi, plot.tau, plot.off, plot.diff, ncol=1, nrow=4), filename=file.path(root, "Figures", "JointModelComparison - PhiMethod.jpeg"), width=15, height=20)

#6. Plot different joint model formulations----

#6a. Phi----
ggplot(est) +
  geom_violin(aes(y=log(phij), x=paste(phimodel, taumodel), colour=phisensor))

#6b. Tau----
ggplot(est) +
  geom_violin(aes(y=log(tauj), x=paste(phimodel, taumodel), colour=tausensor))

#6c. Offset----
ggplot(est) +
  geom_violin(aes(y=offj, x=paste(phimodel, taumodel), colour=sensor)) 
  

#9. Plot phi----
phi.sensor <- est %>% 
  mutate(model = paste0("phi", phimodel, "_tau", taumodel)) %>% 
  dplyr::select(species, phij, model, sensor) %>% 
  pivot_wider(names_from=model, values_from=phij)

ggplot(phi.sensor) +
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(x=phiTM_taunull, y=phiTM_tauTM)) +
  facet_wrap(~sensor, scales="free")

ggplot(phi.sensor) +
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(x=phinull_tauTM, y=phiTM_tauTM)) +
  facet_wrap(~sensor, scales="free")

ggplot(phi.sensor) +
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(x=phiTM_taunull, y=phinull_tauTM)) +
  facet_wrap(~sensor, scales="free")

#9. Plot tau----
tau.sensor <- est %>% 
  mutate(model = paste0("phi", phimodel, "_tau", taumodel)) %>% 
  dplyr::select(species, tauj, model, sensor)  %>% 
  dplyr::filter(tauj < 100) %>% 
  pivot_wider(names_from=model, values_from=tauj)

ggplot(tau.sensor) +
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(x=phiTM_taunull, y=phiTM_tauTM)) +
  facet_wrap(~sensor, scales="free")

ggplot(tau.sensor) +
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(x=phinull_tauTM, y=phiTM_tauTM)) +
  facet_wrap(~sensor, scales="free")

ggplot(tau.sensor) +
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(x=phiTM_taunull, y=phinull_tauTM)) +
  facet_wrap(~sensor, scales="free")

#9. Plot offsets----
off.sensor <- est %>% 
  mutate(model = paste0("phi", phimodel, "_tau", taumodel)) %>% 
  dplyr::select(species, offj, model, sensor) %>% 
  pivot_wider(names_from=model, values_from=offj)

ggplot(off.sensor) +
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(x=phiTM_taunull, y=phiTM_tauTM)) +
  facet_wrap(~sensor, scales="free")

ggplot(off.sensor) +
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(x=phinull_tauTM, y=phiTM_tauTM)) +
  facet_wrap(~sensor, scales="free")

ggplot(off.sensor) +
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(x=phiTM_taunull, y=phinull_tauTM)) +
  facet_wrap(~sensor, scales="free")
