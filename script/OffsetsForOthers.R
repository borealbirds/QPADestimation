#A. PREPARATION###########

#1. Load basic packages----
library(tidyverse)
library(maptools)
library(intrval)
library(raster)
library(lubridate)

#2. Set file path for data & results----
root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/QPAD/OffsetCalculation/VenierFeb2023"

#3. Install most recent version of QPAD----
library(devtools)
install_github("borealbirds/QPAD@dev-v4", force=TRUE)

#4. Load package & offsets----
library(QPAD)
load_BAM_QPAD(4)

#5. Set WD to qpad-offsets package----
setwd("C:/Users/Elly Knight/Documents/BAM/Projects/QPAD/qpad-offsets")

#6. Read raster data----
rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
rd1 <- raster("./data/seedgrow.tif")
crs <- proj4string(rtree)

#7. Source functions----
source("functions.R")

#B. WRANGLE VISIT DATA#####

#1. Read in data----
visit <- read.csv(file.path(root, "PC_offset_run.csv"))

#2. Wrangle----
dat <- visit %>% 
  mutate(date = as.character(ymd(paste0(YEAR, "-", MONTH, "-", DAY))),
         time = as.character(hm(paste0(HOUR, ":", MIN))),
         lon = NA,
         lat = NA,
         dur = NA,
         dist = NA,
         tagmeth="PC") %>% 
  dplyr::select(id, date, time, lon, lat, dur, dist, tagmeth)

#3. Make prediction df----
x <- make_x(dat)
str(x)

#C. CREATE OFFSETS#######

#1. Get list of species----

#2. Make offsets----
SPP <- getBAMspecieslist()
OFF <- matrix(0, nrow(x), length(SPP))
rownames(OFF) <- rownames(x) # your survey IDs here
colnames(OFF) <- SPP

for (spp in SPP) {
  cat(spp, "\n")
  flush.console()
  o <- make_off(spp, x, useMeth)
  OFF[,spp] <- o$offset
}
str(OFF)

#D. PACKAGE####