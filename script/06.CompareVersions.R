# ---
# title: "QPAD estimation - compare versions"
# author: "Elly Knight"
# created: "October 11, 2022"
# ---

library(tidyverse)
library(QPAD)

root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/"

options(scipen=9999)

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

#1. Load results----
load_BAM_QPAD(3)
load(file.path(root, "Results/BAMCOEFS_QPAD_v4.rda"))

#2. Number of species----
spp3 <- getBAMspecieslist()
spp4 <- .BAMCOEFS4$spp

spp34 <- data.frame(species = union(spp3, spp4)) %>%
    left_join(data.frame(species = spp3, version3 = 1)) %>%
    left_join(data.frame(species = spp4, version4 = 1))

sppmissing <- spp34 %>%
    dplyr::filter(is.na(version4))

sppnew <- spp34 %>%
    dplyr::filter(is.na(version3))

#3. Sample size----
n3 <- data.frame(sra.n3=.BAMCOEFS$sra_n, edr.n3=.BAMCOEFS$edr_n,
                 species=.BAMCOEFS$spp) %>% 
  mutate(species = ifelse(species=="GRAJ", "CAJA", species))
n4 <- data.frame(sra.n4=.BAMCOEFS4$sra_n, edr.n4=.BAMCOEFS4$edr_n,
                 species=.BAMCOEFS4$spp)

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

sum(n3$sra.n3)
sum(n4$sra.n4)
sum(n3$edr.n3)
sum(n4$edr.n4)

ggplot(n34) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x=sra.n3, y=sra.n4, fill=version), pch=21, alpha = 0.5, size=4) +
    xlab("V3 sample size") +
    ylab("V4 sample size") +
    scale_fill_manual(values=c("grey80", "blue", "orange"), name="") +
    my.theme

#ggsave(filename="figures/V3V4_samplesize_log.jpeg", width =7, height=6)

ggplot(n34) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x=edr.n3, y=edr.n4, fill=version), pch=21, alpha = 0.5, size=4) +
    scale_fill_manual(values=c("grey80", "blue", "orange"), name="") +
    my.theme

#ggsave(filename="figures/V3V4_samplesize_log.jpeg", width =7, height=6)

ggplot(n34) +
    geom_abline(intercept = 0, slope = 1) +
    geom_text(aes(x=log(sra.n3), y=log(sra.n4), label=species)) +
    xlab("log(V3 sample size)") +
    ylab("log(V4 sample size)") +
    my.theme

#ggsave(filename="figures/V3V4_samplesize_species.jpeg", width =7, height=6)

#4. Null estimates----
est3 <-data.frame()
for(i in 1:length(spp3)){
    sra3 <- exp(.BAMCOEFS$sra_estimates[[i]]$`0`$coefficients)
    edr3 <- exp(.BAMCOEFS$edr_estimates[[i]]$`0`$coefficients)
    est3 <- rbind(est3,
                  data.frame(sra3=sra3, edr3=edr3, species=spp3[i]))
}
rownames(est3) <- NULL
est3 <- est3 %>% 
  mutate(species = ifelse(species=="GRAJ", "CAJA", species))

est4 <-data.frame()
for(i in 1:length(spp4)){
  
    sra4 <- exp(.BAMCOEFS4$sra_estimates[[i]]$`0`$coefficients[1])
    if(.BAMCOEFS4$sra_models[[i,"15"]]==1){
      sra4.pc <- exp(.BAMCOEFS4$sra_estimates[[i]]$'15'$coefficients[[1]])
      sra4.aru <- exp(sum(.BAMCOEFS4$sra_estimates[[i]]$`15`$coefficients[c(1,2)]))
    }
    else{
      sra4.pc <- NA
      sra4.aru <- NA
    }
    edr4 <- exp(.BAMCOEFS4$edr_estimates[[i]]$`0`$coefficients)
    est4 <- rbind(est4,
                  data.frame(sra4=sra4, sra4.aru=sra4.aru, sra4.pc=sra4.pc, edr4=edr4, species=spp4[i]))
}
rownames(est4) <- NULL

est34 <- full_join(est3, est4) %>%
    mutate(edr3 = ifelse(is.na(edr3), 0, edr3),
           sra3 = ifelse(is.na(sra3), 0, sra3),
           edr4 = ifelse(is.na(edr4), 0, edr4),
           sra4 = ifelse(is.na(sra4), 0, sra4),
           sra4.aru = ifelse(is.na(sra4.aru), 0, sra4.aru),
           sra4.pc = ifelse(is.na(sra4.pc), 0, sra4.pc)) %>%
    mutate(version = case_when(edr3==0 ~ "V4",
                               edr4==0 ~ "V3",
                               !is.na(edr4) ~ "Both")) %>%
    left_join(n34) %>% 
  mutate(sra34 = sra3/sra4,
         edr34 = edr3/edr4,
         sra34.abs = abs(1-sra34),
         edr34.abs = abs(1-edr34))

write.csv(est34, file.path(root, "Results/QPADV3V4NullEstimates.csv"), row.names=FALSE)

ggplot(est34) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x=sra3, y=sra4, fill=version), pch=21, alpha = 0.5, size=4) +
#    geom_text(aes(x=sra3, y=sra4, label=species)) +
    xlab("V3 availability estimate (phi)") +
    ylab("V4 availability estimate (phi)") +
    scale_fill_manual(values=c("grey80", "blue", "orange"), name="") +
    my.theme

#ggsave(filename="figures/V3V4_phi_species.jpeg", width =7, height=6)

ggplot(est34) +
    geom_abline(intercept = 0, slope = 1) +
    geom_point(aes(x=edr3, y=edr4, fill=version), pch=21, alpha = 0.5, size=4) +
#    geom_text(aes(x=edr3, y=edr4, label=species)) +
    xlab("V3 perceptability estimate (tau)") +
    ylab("V4 perceptability estimate (tau)") +
    scale_fill_manual(values=c("grey80", "blue", "orange"), name="") +
    my.theme

#ggsave(filename="figures/V3V4_tau_species.jpeg", width =7, height=6)

ggplot(est34) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x=sra4, y=sra4.aru, fill=version), pch=21, alpha = 0.5, size=4) +
  geom_smooth(aes(x=sra4, y=sra4.aru), method="lm") +
#  geom_text(aes(x=sra4, y=sra4.aru, label=species)) +
  xlab("V4 point count availability estimate (phi)") +
  ylab("V4 SPT availability estimate (phi)") +
  scale_fill_manual(values=c("grey80", "blue", "orange"), name="") +
  my.theme

ggsave(filename="figs/ARUvsPC.jpeg", width=8, height=6, units="in")

lm1 <- lm(sra4.aru ~ sra4, data=est34)
summary(lm1)

#5. Sample size vs estimate----
ggplot(est34 %>% 
         dplyr::filter(version=="Both")) +
  geom_point(aes(x=sra.n, y=sra34.abs, colour=sra.n4)) +
  geom_smooth(aes(x=sra.n, y=sra34.abs)) +
  scale_colour_viridis_c()+
  ylim(c(0,1))

ggplot(est34 %>% 
         dplyr::filter(version=="Both")) +
  geom_point(aes(x=sra.n4, y=sra34.abs, colour=sra.n4)) +
  geom_smooth(aes(x=sra.n4, y=sra34.abs)) +
  scale_colour_viridis_c() +
  ylim(c(0,1))

ggplot(est34 %>% 
         dplyr::filter(version=="Both")) +
  geom_point(aes(x=edr.n, y=edr34.abs, colour=edr.n4)) +
  geom_smooth(aes(x=edr.n, y=edr34.abs)) +
  scale_colour_viridis_c()

ggplot(est34 %>% 
         dplyr::filter(version=="Both")) +
  geom_point(aes(x=edr.n4, y=edr34.abs, colour=edr.n4)) +
  geom_smooth(aes(x=edr.n4, y=edr34.abs)) +
  scale_colour_viridis_c()