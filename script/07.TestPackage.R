# ---
# title: "QPAD estimation - test QPAD & qpad-offsets"
# author: "Elly Knight"
# created: "January 4, 2023"
# ---

#NOTES################################

#Use of this script requires updating the QPAD R package:
#1. Save the output as an RDS to the inst/estimates folder of the QPAD package directory on your local (this is part of the 06.Package script)
#2. Create a new R script in the inst/estimates folder that loads that rds and writes it to a hidden .BAMCOEFS environment. Use the V4 version as a template
#3. Add the new version to the load_BAM_QPAD function in the R folder.
#4. Push your changes to the QPAD repo (use a development branch...)

#If model structure is changed, use of this script requires updating the qpad-offsets repo:
#1. Update the make_x() and make_off() functions in the functions.R script.

#If model structure is changed, use of this script may also require updating the QPAD package
#1. Update the 

#A. TEST QPAD PACKAGE####

#1. Download test branch of QPAD package----
library(devtools)
install_github("borealbirds/QPAD@dev-v4", force=TRUE)

#2. Load QPAD package----
library(QPAD)

#3. Check that all versions of estimates load----
load_BAM_QPAD(2)
length(.BAMCOEFS$spp) #75
rm(.BAMCOEFS)

load_BAM_QPAD(3)
length(.BAMCOEFS$spp) #151
rm(.BAMCOEFS)

load_BAM_QPAD(4)
length(.BAMCOEFS$spp) #186

#B. TEST qpad-offsets REPO####

#1. Set WD to qpad-offsets package----
setwd("C:/Users/elly/Documents/BAM/QPAD/qpad-offsets")

#2. Load packages----
library(maptools)
library(intrval)
library(raster)

#3. Read raster data----
rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
rd1 <- raster("./data/seedgrow.tif")
crs <- proj4string(rtree)

#4. Source functions----
source("functions.R")

#5. Make up some data----
dat <- data.frame(date = c("2019-06-07", "2019-06-17", "2019-06-27"),
                  time = rep("05:20", 3),
                  lon = seq(-115, -113, 1),
                  lat = seq(53, 55, 1),
                  dur = rep(10, 3),
                  dist = rep(100, 3),
                  tagmeth = c("PC", "1SPT", "1SPM"))

#6. Make prediction df----
x <- make_x(dat)
str(x)

#7. Make offsets----
spp <- "OVEN"
useMeth <- "y"
o <- make_off(spp, x, useMeth)
str(o)

#C. TEST README SCRIPT####

setwd("C:/Users/elly/Documents/BAM/QPAD/qpad-offsets")

## load packages
library(QPAD)
library(maptools)
library(intrval)
library(raster)

## load v4 estimates
load_BAM_QPAD(version = 4)
if (!getBAMversion() %in% c("3", "4"))
  stop("This script requires BAM version 3 or version 4")

## read raster data
rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
rd1 <- raster("./data/seedgrow.tif")
crs <- proj4string(rtree)

## source functions
source("functions.R")

## dataframe
dat <- data.frame(date = c("2019-06-07", "2019-06-17", "2019-06-27"),
                  time = rep("05:20", 3), 
                  lon = rep(-115, 3),
                  lat = rep(53, 3),
                  dur = rep(10, 3), 
                  dist = rep(100, 3),
                  tagmeth = rep("PC", 3)) 

## timezone argument
tz <- "local"

## organize predictors
x <- make_x(dat, tz)
str(x)
# 'data.frame':	3 obs. of  9 variables:
#   $ TSSR  : num  0.0024 0.0044 0.0026
# $ JDAY  : num  0.43 0.458 0.485
# $ DSLS  : num  0.11 0.137 0.164
# $ LCC2  : Factor w/ 2 levels "Forest","OpenWet": 2 2 2
# $ LCC4  : Factor w/ 4 levels "DecidMixed","Conif",..: 3 3 3
# $ TREE  : num  0.3 0.3 0.3
# $ MAXDUR: num  10 10 10
# $ MAXDIS: num  1 1 1
# $ TM    : chr  "PC" "PC" "PC"

## species of interest
spp <- "OVEN"

## take method into acount
useMeth <- "y"

o <- make_off(spp, x, useMeth)
str(o)
# 'data.frame':	3 obs. of  5 variables:
#   $ p         : num  0.982 0.975 0.966
# $ q         : num  0.58 0.58 0.58
# $ A         : num  3.14 3.14 3.14
# $ correction: num  1.79 1.78 1.76
# $ offset    : num  0.582 0.575 0.566

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
# num [1:3, 1:187] -0.0683 -0.0683 -0.0683 -0.01 -0.0111 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:3] "1" "2" "3"
# ..$ : chr [1:187] "ACFL" "ALFL" "AMCR" "AMGO" ...