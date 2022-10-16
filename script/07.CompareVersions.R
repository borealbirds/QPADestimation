# ---
# title: "QPAD estimation - compare versions"
# author: "Elly Knight"
# created: "October 11, 2022"
# ---

library(tidyverse)
library(QPAD)

options(dplyr.summarise.inform = FALSE, scipen=9999)

load_BAM_QPAD(3)

load("results/availability_results_2022-10-06.Rdata")
load("results/perceptability_results_2022-10-06.Rdata")

#1. Number of species----
#1a. Version 3
spp3 <- getBAMspecieslist()

#1b. Version 4
spp4 <- species.use.avail %>%
    dplyr::select(species) %>%
    unique() %>%
    inner_join(species.use.percep %>%
                   dplyr::select(species) %>%
                   unique())

#1c. New species
sppnew <- spp4 %>%
    dplyr::filter(!species %in% spp3)

#2. Samples size----

#2a. Version 3
n3 <- data.frame(sra.n3=.BAMCOEFS$sra_n, edr.n3=.BAMCOEFS$edr_n,
                 species=.BAMCOEFS$spp)

#2b. Version 4
n4 <- data.frame(species = spp4$species) %>%
    left_join(species.use.avail %>%
                  dplyr::select(species, nobs) %>%
                  unique() %>%
                  rename(sra.n4=nobs)) %>%
    left_join(species.use.percep %>%
                  dplyr::select(species, nobs) %>%
                  unique() %>%
                  rename(edr.n4=nobs))

#2c. Compare
n34 <- full_join(n3, n4) %>%
    mutate(edr.n4 = as.numeric(edr.n4),
           sra.n4 = as.numeric(sra.n4),
           edr.n3 = ifelse(is.na(edr.n3), 0, edr.n3),
           sra.n3 = ifelse(is.na(sra.n3), 0, sra.n3),
           edr.n4 = ifelse(is.na(edr.n4), 0, edr.n4),
           sra.n4 = ifelse(is.na(sra.n4), 0, sra.n4),
           edr.n = edr.n3/edr.n4,
           sra.n = sra.n3/sra.n4) %>%
    mutate(version = case_when(edr.n3==0 ~ "V4",
                           edr.n4==0 ~ "V3",
                           !is.na(edr.n) ~ "Both"))

ggplot(n34) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x=sra.n3, y=sra.n4, colour=version))

sum(n34$edr.n4/sum(n34$edr.n3))
sum(n34$sra.n4/sum(n34$sra.n3))

#3. Top model for availability----

#3c. Version 3
mod3 <- data.frame(sra.mod3 = .BAMCOEFS$sra_aiccbest,
                   edr.mod3 = .BAMCOEFS$edr_aicbest)

#4. Top model for perceptibility----

#5. Estimates----
est3 <- data.frame(exp(t(sapply(spp3, function(i) unlist(coefBAMspecies(i))))),
                   species = .BAMCOEFS$spp)
colnames(est3) <- c("sra3", "edr3", "species")

est4 <- species.use.avail %>%
    group_by(species) %>%
    arrange(aicc) %>%
    dplyr::filter(row_number()==1) %>%
    ungroup() %>%
    mutate(sra4 = exp(log.phi_.Intercept.)) %>%
    dplyr::select(species, sra4) %>%
    inner_join(species.use.percep %>%
                   group_by(species) %>%
                   arrange(aicc) %>%
                   dplyr::filter(row_number()==1) %>%
                   ungroup() %>%
                   mutate(edr4 = exp(log.tau_.Intercept.)) %>%
                   dplyr::select(species, edr4))

est34 <- full_join(est3, est4) %>%
    mutate(edr4 = as.numeric(edr4),
           sra4 = as.numeric(sra4),
           edr3 = ifelse(is.na(edr3), 0, edr3),
           sra3 = ifelse(is.na(sra3), 0, sra3),
           edr4 = ifelse(is.na(edr4), 0, edr4),
           sra4 = ifelse(is.na(sra4), 0, sra4)) %>%
    mutate(version = case_when(edr3==0 ~ "V4",
                               edr4==0 ~ "V3",
                               !is.na(edr4) ~ "Both")) %>%
    left_join(n34)

ggplot(est34) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x=sra3, y=sra4, colour=version)) +
    ylim(c(0, 2))

ggplot(est34) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x=edr3, y=edr4, colour=version)) +
    ylim(c(0, 2))

#6. Effect of tree----
tree3 <- data.frame(exp(t(sapply(spp3, function(i) unlist(coefBAMspecies(i))))),
                   species = .BAMCOEFS$spp)
