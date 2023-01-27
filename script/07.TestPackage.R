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
#1. 

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
length(.BAMCOEFS$spp) #187

#B. TEST qpad-offsets REPO####

#1. Set WD to qpad-offsets package----
setwd("C:/Users/Elly Knight/Documents/BAM/Projects/QPAD/qpad-offsets")

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
                  tagmeth = rep("PC", 3))

#6. Make prediction df----
x <- make_x(dat)
x

#7. Make offsets----
spp <- "OVEN"
meth <- 1
off <- make_off(spp, x, meth)
off
