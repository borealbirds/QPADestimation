# ---
# title: "QPAD estimation - clean visit data"
# author: "Elly Knight"
# created: "July 24, 2022"
# updated: "September 22, 2022"
# ---

#NOTE: INVESTIGATE MISSING LAT LONS - email sent####
#NOTE: NEED TO GET ARU MODEL - waiting for next WT update Sept 28####
#NOTE: NEED TO DECIDE HOW TO CREATE ALPHA CODE OR USE HISTORIC? IF HISTORIC, DOWNLOAD FROM REPO####
#NOTE: HOW TO DEAL WITH SURVEYS TOO FAR NORTH FOR SUNRISE????####
#Peter just had NAs in his data, so presumably removed them
#NOTE: DATA PRIOR TO 2000 USES 2000 VALUE FOR CANOPY COVER. PONDER####
#NOTE: DATA SOUTH OF HERMOSILLA JUST USES COPERNICUS 2015 CLASS####
#NOTE: COME BACK TO FILTERING DATASETS FOR ARU DATA####
#NOTE: FILL GREENDAY NAs WITH VALUE FROM ANOTHER YEAR####

library(tidyverse) #basic data wrangling
library(lubridate) #date manipulation
library(suncalc) #sunrise time retrieval
library(sf) #spatial manipulation
library(terra) #raster creation
library(downloader) #download spatial data
library(data.table) #to collapse rgee output

#Load previous dataset to look at it----
#load into separate environment to avoid overwriting things
e <- new.env()
load("data/new_offset_data_package_2017-03-01.Rdata", envir = e)
#look at data structure
str(e$dat)

#1. Load in new dataset----
load("data/wildtrax_data_2022-07-24.Rdata")

#2. Subset to visits & filter----
#Remove surveys with no location
#Remove surveys with no temporal or spatial structure within a visit
dat <- raw %>%
    dplyr::select(organization, project, location, latitude, longitude, observer, date, recording_date, distanceMethod, durationMethod) %>%
    unique() %>%
    dplyr::filter(!is.na(latitude),
                  latitude > 0,
                  !(distanceMethod %in% c("0m-INF", NA, "0m-INF-ARU", "UNKNOWN") &
                        durationMethod %in% c("0-3min", "0-5min", "0-10min", "0-2min", "0-20min", "UNKNOWN", NA)))

#Take subsample for testing
#dat <- sample_n(dat, 1000)

#3. Determine survey type---
type <- dat %>%
    mutate(type = case_when(!is.na(date) ~ "human",
                            !is.na(recording_date) ~ "aru"))

#3. Wrangle temporal variables----
temporal <- type %>%
    mutate(datetime = case_when(type=="human" ~ ymd_hms(date),
                                type=="aru" ~ ymd_hms(recording_date)),
           year = year(datetime),
           julian = yday(datetime),
           start_time = hour(datetime) + minute(datetime)/60) %>%
    dplyr::filter(year > 1900)

#4. Get sunrise time----
sun <- temporal %>%
    mutate(date = ymd(str_sub(datetime, 1, 10))) %>%
    rename(lat = latitude, lon = longitude)

sun$sunrise <- getSunlightTimes(data=sun, keep="sunrise")$sunrise
sun$tssr <- as.numeric(difftime(sun$datetime, sun$sunrise), units="hours")

#5. Lookup method codes----


dist_codes <- read.csv("data/distance_band_codes.csv")
dur_codes <- read.csv("data/durint.csv")

#6. Get region----

#6a. Download data
# temp <- tempfile()
# download("https://birdscanada.org/download/gislab/bcr_terrestrial_shape.zip", temp)
# unzip(zipfile=temp, exdir="gis")
# unlink(temp)

#6b.Read in & wrangle shapefile
bcrshp <- read_sf("gis/BCR_Terrestrial/BCR_Terrestrial_master.shp") %>%
    dplyr::select(BCR, PROVINCE_S, COUNTRY) %>%
    rename(bcr = BCR, province = PROVINCE_S, country=COUNTRY) %>%
    dplyr::filter(country %in% c("USA", "CANADA")) %>%
    mutate(countryid = as.numeric(as.factor(country)),
           provinceid = as.numeric(as.factor(province))) %>%
    st_transform(crs=3857) %>%
    st_make_valid() %>%
    vect()

#6c. Create rasters (much faster than from polygon)
r <- rast(ext(bcrshp), resolution=1000)

bcr <- rasterize(x=bcrshp, y=r, field="bcr")
province <- rasterize(x=bcrshp, y=r, field="provinceid")
country <- rasterize(x=bcrshp, y=r, field="countryid")

#6d. Extract values
bcrid <- sun %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    st_transform(crs=3857) %>%
    vect() %>%
    terra::extract(x=bcr)

provinceid <- sun %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    st_transform(crs=3857) %>%
    vect() %>%
    terra::extract(x=province)

countryid <- sun %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    st_transform(crs=3857) %>%
    vect() %>%
    terra::extract(x=country)

#6e. Create lookup table to join back ids
bcrtbl <- st_as_sf(bcrshp) %>%
    as.data.frame() %>%
    dplyr::select(-geometry) %>%
    unique()

#6f. Put together
region <- sun %>%
    cbind(bcrid %>% dplyr::select(bcr)) %>%
    cbind(provinceid %>% dplyr::select(provinceid)) %>%
    cbind(countryid %>% dplyr::select(countryid)) %>%
    left_join(bcrtbl) %>%
    dplyr::select(-provinceid, -countryid)

#7. Get covariates----

#7a. Initialize rgee
library(rgee)
ee_Initialize()
ee_check()

#7b. Set up to loop through data years
years <- sort(unique(region$year))

out.list <- list()
for(i in 1:length(years)){

    #7c. Filter to year and subset into chunks of 5000
    dat <- region %>%
        st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
        dplyr::filter(year==years[i]) %>%
        mutate(loop = ceiling(row_number()/5000))

    dat.out <- data.frame()
    #7d. Set up loop
    for(j in 1:max(dat$loop)){

        #7e. Format for rgee
        dat.j <- dat %>%
            dplyr::filter(loop==j)

        dat.ee <- dat.j %>%
            dplyr::select(geometry) %>%
            sf_as_ee()

        #7f. Get Hansen dataset for canopy cover
        tree <- ee$Image('UMD/hansen/global_forest_change_2021_v1_9')

        dat.tree <- ee_extract(
            x=tree,
            y=dat.ee,
            scale=100,
            sf=FALSE
        ) %>%
            dplyr::select(-first_b30, -first_b40, -first_b50, -first_b70, -last_b30, -last_b40, -last_b50, -last_b70, -datamask)

        #7g. Get Hermosilla landcover
        yeartbl <- data.frame(year = years,
                              year.i = c(1988, 1990, 1990, 1993, 1993, 1995, 1995, 1997, 1997, 1999, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2007, 2007, 2008, 2010, 2010, 2010, 2010, 2015, 2015, 2016, 2018, 2018, 2019, 2019))
        year.i <- yeartbl$year.i[i]

        start<-paste0(year.i, "-01-01")
        end<-paste0(year.i,"-12-31")

        lc <- ee$ImageCollection('projects/sat-io/open-datasets/CA_FOREST_LC_VLCE2')$filterDate(start, end)

        dat.lc <- ee_extract(
            x=lc,
            y=dat.ee,
            scale=100,
            sf=FALSE
        )
        colnames(dat.lc) <- "hermosilla"

        #7h. Get 2015 copernicus to fill gaps
        start <- "2015-01-01"
        end <- "2015-12-31"

        cop <- ee$ImageCollection('COPERNICUS/Landcover/100m/Proba-V-C3/Global')$filterDate(start, end)

        dat.cop <- ee_extract(
            x=cop,
            y=dat.ee,
            scale=100,
            sf=FALSE
        )
        colnames(dat.cop) <- c("bare", "crop", "density", "copernicus", "prob", "forest", "grass", "moss", "shrub", "snow", "tree", "urban", "permanentwater", "seasonalwater")


        #7i. Get modis greenup day
        if(years[i] < 2001){year.i <- 2001}
        if(years[i] >= 2001 & years[i] <= 2019){year.i <- years[i]}
        if(years[i] > 2019){year.i <-2019}

        start<-paste0(year.i, "-01-01")
        end<-paste0(year.i,"-12-31")

        green <- ee$ImageCollection('MODIS/006/MCD12Q2')$select('Greenup_1')$filterDate(start, end)

        dat.green <- ee_extract(
            x=green,
            y=dat.ee,
            scale=100,
            sf=FALSE
        )
        colnames(dat.green) <- "green"
        dat.green$greenday <- dat.green$green - as.numeric(ymd(start))

        #7j.Put everything together
        dat.out <- data.frame(st_coordinates(dat.j)) %>%
            rename(lat = Y, lon = X) %>%
            cbind(data.frame(dat.j) %>%
                      dplyr::select(-geometry)) %>%
            cbind(dat.tree, dat.lc, dat.cop, dat.green) %>%
            rbind(dat.out)

    }

    #7k. Save to list
    out.list[[i]] <- dat.out

    print(paste0("Finished year ", years[i], ": ", i, " of ", length(years), " years"))

}

#8. Create lookup tables for landcover classes----
#hermosilla codes
hermosillacodes <- data.frame(hermosilla = c(0, 20, 31, 32, 33, 40, 50, 80, 81, 100, 210, 220, 230),
                              hermosilladesc = c("unclassified", "water", "snow", "rock", "barren", "bryoid", "shrub", "wetland", "wetlandtreed", "herb", "conifer", "deciduous", "mixedwood"),
                              hermosilla1 = c(NA, "open", "open", "open", "open", "open", "open", "open", "treed", "open", "treed", "treed", "treed"),
                              hermosilla2 = c(NA, "open", "open", "open", "open", "open", "open", "wetland", "conifer", "open", "conifer", "deciduous", "mixed"))

#copernicus codes
copernicuscodes <- data.frame(copernicus = c(0, 20, 30, 40, 50, 60, 70, 80, 90, 100, 111, 112, 113, 114, 115, 116, 121, 122, 123, 124, 125, 126, 200),
                              copernicusdesc = c("unclassified", "shrub", "herb", "crop", "urban", "barren", "snow", "water", "wetland", "bryoid", "coniferousclosed", "deciduousclosed", "coniferousclosed", "deciduousclosed", "mixedclosed", "treedclosed", "coniferousopen", "deciduousopen", "coniferousopen", "deciduousopen", "mixedopen", "treedopen", "marine"),
                              copernicus1 = c(NA, "open", "open", "open", "open", "open", "open", "open", "open", "open", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "treed", "open"),
                              copernicus2 = c(NA, "open", "open", "open", "open", "open", "open", "open", "wetland", "open", "conifer", "deciduous", "conifer", "deciduous", "mixed", "mixed", "conifer", "deciduous", "conifer", "deciduous", "mixed", "mixed", "open"))

#9. Tidy and save----
visit <- rbindlist(out.list, fill=TRUE) %>%
    data.frame() %>%
    mutate(lossyear = ifelse(is.na(lossyear), 0, lossyear+2000),
           cover = ifelse(lossyear >= year, 0, treecover2000),
           cover = ifelse(year >= 2012 & gain==1, 1, cover)) %>%
    left_join(hermosillacodes) %>%
    left_join(copernicuscodes) %>%
    mutate(lc1 = ifelse(is.na(hermosilla1), copernicus1, hermosilla1),
           lc2 = ifelse(is.na(hermosilla2), copernicus2, hermosilla2)) %>%
    mutate(tsg = (as.numeric(date)-greenday)/365) %>%
    dplyr::select(organization, project, distanceMethod, durationMethod, location, lon, lat, date, datetime, year, julian, tssr, type, observer, bcr, province, country, cover, lc1, lc2, greenday, tsg)

save(visit, file="data/visit_data_2022-07-24.Rdata")
