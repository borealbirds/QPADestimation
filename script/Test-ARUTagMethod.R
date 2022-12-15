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

root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/"

#1. Create list of models----
#jday = day of year as a decimal between 0 and 1
#tssr = time since sunrise as a decimal between 0 and 1
#tsg = days since start of seedgrowth from seedgrow layers

mods <- list(
  ~ 1,
  ~ tagMethod)
names(mods) <- c(0:1)
modnames <- list(
  "0"=c("(Intercept)"),
  "1"=c("(Intercept)", "tagMethod"))

#2. Load data----
load(file.path(root, "Data/qpadv4_clean.Rdata"))
visit <- dplyr::filter(visit, tagMethod!="None", sensor=="ARU")

#set intercept as 1spt-----
visit$tagMethod <- factor(visit$tagMethod, levels=c("1SPT", "1SPM"))

#3. Create design lookup table that describes duration method for each protocol----
#filter out duration methods that aren't appropriate for removal modelling (only have 1 time bin)
durdesign <- visit %>%
  dplyr::select(durationMethod) %>%
  unique() %>%
  dplyr::filter(!durationMethod %in% c("UNKNOWN", "0-10min", "0-20min", "0-5min", "0-3min", "0-2min"),
                !is.na(durationMethod)) %>%
  mutate(dm = str_sub(durationMethod, -100, -4)) %>%
  separate(dm, into=c("t00", "t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), remove=TRUE, sep="-") %>%
  dplyr::select(-t00) %>%
  mutate_at(c("t01", "t02", "t03", "t04", "t05", "t06", "t07", "t08", "t09", "t10"), ~as.numeric(.))

durdesign.long <- durdesign %>% 
  pivot_longer(t01:t10, values_to="end", names_to="position")

#4. Get list of species to process----
spp <- species %>%
  dplyr::filter(Singing_birds==TRUE) %>%
  left_join(bird %>%
              dplyr::select(species) %>%
              unique()) %>%
  arrange(species)

#5. Set up loop for species----
avail <- list()
for(i in 1:nrow(spp)){
  
  #6. Filter abundance data for species---
  # filter out observations with unknown duration method or interval
  # filter out observations where interval does not match method
  # filter to observations with covariates
  bird.i <- bird %>%
    dplyr::filter(species==spp$species[i],
                  durationMethod %in% durdesign$durationMethod,
                  !durationInterval %in% c("UNKNOWN", "before or after/incidental")) %>% 
    separate(durationInterval, into=c("start", "end"), sep="-", remove=FALSE, extra="drop", fill="right") %>%
    mutate(start = as.numeric(start),
           end = as.numeric(str_sub(end, -100, -4))) %>%
    left_join(durdesign.long, by=c("durationMethod", "end")) %>% 
    dplyr::filter(!is.na(position)) %>% 
    group_by(id, durationMethod, durationInterval, position) %>%
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
      dplyr::filter(id %in% unique(bird.i$id),
                    !is.na(tssr),
                    !is.na(tsg),
                    !is.na(jday),
                    durationMethod %in% durdesign$durationMethod) %>%
      arrange(id) %>%
      dplyr::select(id, durationMethod, tagMethod, year, jday, tssr, tsg)
    
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
      dplyr::select(id, position, abundance) %>%
      arrange(position)  %>% 
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
spp <- spp %>% 
  dplyr::filter(species!="WATA")
names(avail) <- spp$species

#14. Save out results----
save(durdesign, avail, file=file.path(root, "Results/availability-tagMethod.Rdata")) 

#15. Summarize----

#remove species that didn't run
avail.fit <- avail[sapply(avail, length) > 0]
avail.fit <- avail.fit[!t(sapply(avail.fit, function(z) ifelse(sapply(z, inherits, what="try-error"), 0L, 1L)))[,1]==0]

#16. Get loglik values---
out <- data.frame()
for(i in 1:length(avail.fit)){
  
   out <- data.frame(loglik.null = avail.fit[[i]][[1]]$loglik,
                       loglik.method = avail.fit[[i]][[2]]$loglik,
                       est.null = exp(avail.fit[[i]][[1]]$coefficients),
                       est.method = exp(sum(avail.fit[[i]][[2]]$coefficients)),
                       species = names(avail.fit)[[i]]) %>% 
     rbind(out)
}

out.aic <- out %>% 
  mutate(aic.null = -2*loglik.null + 2,
         aic.method = -2*loglik.method + 4,
         deltaaic = aic.null - aic.method,
         better = ifelse(deltaaic > 2, "Tagmethod model selected", "Null model selected"))
table(out.aic$better)

ggplot(out.aic) +
  geom_point(aes(x=est.null, y=est.method)) + 
  facet_wrap(~better, scales="free") +
  geom_abline(aes(intercept=0, slope=1)) +
  xlim(0,3) +
  ylim(0,3) +
  xlab("Null model estimate") +
  ylab("Tag method estimate")

ggsave(filename="figs/TagMethod.jpeg", width=8, height = 4)