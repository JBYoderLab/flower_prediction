# Build a static species distribution model to project suitable habitat changes
# Assumes local environment
# jby 2025.05.19

# setwd("~/Documents/Active_projects/flower_prediction")

# Clear the environment and load key packages
rm(list=ls())

library("tidyverse") 

library("terra")
library("geodata")
library("rgbif")
library("prism")
library("rnaturalearth")
library("rnaturalearthdata")
library("CoordinateCleaner")

library("sf")
library("embarcadero")

# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia: 57250

prism_set_dl_dir("../data/PRISM") # set PRISM data directory

# get the USFS range polygon for toyon
usfs.Range <- read_sf("../data/spatial/wpetry-USTreeAtlas-4999258/shp/photarbu/", layer="photarbu", crs=4326)
usfs.buff <- st_transform(st_buffer(st_transform(usfs.Range, crs=3857), 10000), crs=4326) %>% st_intersection(filter(ne_countries(scale=10, continent="north america", returnclass="sf"), name_en=="United States of America")) # 10km buffer?

plot(usfs.Range)
plot(usfs.buff)



# new composite
broad.Range <- read_sf(paste0("output/broad_range_polygon_", taxon, ".shp"))
broad.Range

ggplot() + geom_sf(data=broad.Range)

#-----------------------------------------------------------
# Load and manage input data: occurrences from GBIF

taxonKey <- name_backbone("Heteromeles arbutifolia")$usageKey
taxonKey

gbif_download <- occ_download(pred("taxonKey", taxonKey), pred_lt("coordinateUncertaintyInMeters", 500), format = "SIMPLE_CSV") 

# This doesn't immediately download your data --- it puts the request in a queue
# on the GBIF servers. The text output to the screen has the details, and code
# to use to actually download the data and load it directly into memory in R
# when it's ready.

# This function will pause until there's a download waiting and then let you
# know when it's ready 
occ_download_wait(gbif_download)

# And this code will complete the download and read it into memory
sppOcc_data <- occ_download_get(gbif_download) %>% occ_download_import()

# Take a look at your data
glimpse(sppOcc_data)

# And now take a look at the actual spatial distribution of the occurrences.
# Because GBIF draws from a lot of sources, erroneous records are very possible.

ggplot() + geom_sf(data=ne_countries(continent = "north america", returnclass = "sf")) + 
	geom_point(data=sppOcc_data, aes(x=decimalLongitude, y=decimalLatitude))

# If you see occurrences wildly outside the range you expect, you should filter
# them out, using latitude and longitude. GBIF has a whole pipeline of suggested
# filters to improve data quality:

sppOcc_clean <- sppOcc_data %>%
	setNames(tolower(names(.))) %>% # set lowercase column names to work with CoordinateCleaner
	filter(occurrencestatus  == "PRESENT") %>%
	filter(!basisofrecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN")) %>%
	filter(coordinateprecision < 0.01 | is.na(coordinateprecision)) %>% 
	filter(coordinateuncertaintyinmeters < 5000 | is.na(coordinateuncertaintyinmeters)) %>%
	filter(!coordinateuncertaintyinmeters %in% c(301,3036,999,9999)) %>% 
	filter(!decimallatitude == 0 | !decimallongitude == 0) %>%
	cc_cen(lon="decimallongitude", lat="decimallatitude", buffer = 2000) %>% # remove country centroids within 2km 
	cc_cap(lon="decimallongitude", lat="decimallatitude", buffer = 2000) %>% # remove capitals centroids within 2km
	cc_inst(lon="decimallongitude", lat="decimallatitude", buffer = 2000) %>% # remove zoo and herbaria within 2km 
	cc_sea(lon="decimallongitude", lat="decimallatitude") %>% # remove from ocean 
	distinct(decimallongitude,decimallatitude,specieskey,datasetkey, .keep_all = TRUE) %>% 
	filter(decimallongitude < -110, decimallatitude  > 32) %>% 
	filter(!(decimallongitude < -120 & decimallatitude < 35)) # a priori bad records?

# Now, check the range of your cleaned occurrences and restrict your map to that
# range:
range(sppOcc_clean$decimallongitude)
range(sppOcc_clean$decimallatitude)

ggplot() + geom_sf(data=ne_countries(continent = "north america", returnclass = "sf")) + 
	geom_sf(data=spp.Range, color=NA, fill="forestgreen", alpha=0.5) + 
	geom_point(data=sppOcc_clean, aes(x=decimallongitude, y=decimallatitude), size=0.5, alpha=0.5) +
	coord_sf(xlim = c(-125,-114), ylim = c(28,42), expand = FALSE)


# SAVE the data at this point
write.table(sppOcc_clean, paste0("data/GBIF_cleaned_", taxon, ".csv"), sep=",", col.names=TRUE, row.names=FALSE)

# RELOAD it (if you skipped down to here)
sppOcc_clean <- read.csv(paste0("data/GBIF_cleaned_", taxon, ".csv"))


# what is the TIMING of these records?
glimpse(sppOcc_clean)
hist(sppOcc_clean$year)
table(sppOcc_clean$year) # yeah okay ... not a lot of early-20th century records!
table(sppOcc_clean$year<1950) # LOL a whole 27
table(sppOcc_clean$year<1980) # ... 78 ...
table(sppOcc_clean$year<1990) # ... 103 ...
table(sppOcc_clean$year<2000) # ... 120 ...

# let's plot them ...
sppOcc_clean$late <- sppOcc_clean$year>1990

ggplot() + geom_sf(data=ne_countries(continent = "north america", returnclass = "sf")) + 
	geom_sf(data=spp.Range, color=NA, fill="forestgreen", alpha=0.5) + 
	geom_point(data=sppOcc_clean, aes(x=decimallongitude, y=decimallatitude), size=0.5, alpha=0.1) +
	facet_wrap("late") + 
	coord_sf(xlim = c(-125,-114), ylim = c(28,42), expand = FALSE) +
	theme(panel.background=element_rect(fill="slategray3", color=NA), panel.grid=element_line(color="slategray2"))

#-------------------------------------------------------------------------
# Building a better distribution polygon

countries <- ne_countries(scale=10, continent="north america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

sppOcc_clean_recent <- sppOcc_clean %>% filter(year>1990)
glimpse(sppOcc_clean_recent) # 17k, whee

usfs.Range <- read_sf("../data/spatial/wpetry-USTreeAtlas-4999258/shp/photarbu/", layer="photarbu", crs=4326)
usfs.Range
plot(usfs.Range)

sppOcc_composite <- sppOcc_clean_recent %>% # early records
	st_as_sf(coords=c("decimallongitude", "decimallatitude"), crs=4326) %>% # coverted to sf, scaled in deg
	st_transform(crs=3857) %>% # transformed to units in meters
	st_buffer(15000) %>% st_union() %>% # buffer by ... 15km?
	st_union(., st_transform(usfs.Range, crs=3857)) %>% st_union() %>% # merge with USFS polygon
	st_concave_hull(ratio=0.125) %>% 
	st_simplify(preserveTopology=TRUE, dTolerance=5000) %>% st_buffer(50000) %>% 
	st_transform(crs=4326) %>% st_intersection(filter(countries, name_en=="United States of America")) %>% st_union() %>% st_as_sf() # back to lat-lon

plot(sppOcc_composite) # okay that might be enough to work with?

st_write(sppOcc_composite, paste0("output/broad_range_polygon_", taxon, ".shp"), append=FALSE)

ggplot() + geom_sf(data=ne_countries(scale=10, continent = "north america", returnclass = "sf")) + 
	geom_sf(data=sppOcc_composite, color=NA, fill="forestgreen", alpha=0.5) + 
	geom_sf(data=usfs.Range, color="white", fill=NA) + 
	geom_point(data=sppOcc_clean_recent, aes(x=decimallongitude, y=decimallatitude), size=0.5, alpha=0.1) +
	coord_sf(xlim = c(-125,-114), ylim = c(28,42), expand = FALSE) +
	theme(panel.background=element_rect(fill="slategray3", color=NA), panel.grid=element_line(color="slategray2"))



#-------------------------------------------------------------------------
# early- and recent-period climate averages

files_early <- list.files("../data/spatial/custom_BIOCLIM/", pattern="1901_1930", full=TRUE)

BIOCLIM_1901_1930 <- rast(files_early)
names(BIOCLIM_1901_1930) <- gsub(".+BIO_\\d\\d_(.+)_1901_1930\\.tiff", "\\1", files_early)

BIOCLIM_1901_1930

files_recent <- list.files("../data/spatial/custom_BIOCLIM/", pattern="1991_2020", full=TRUE)

BIOCLIM_1991_2020 <- rast(files_recent)
names(BIOCLIM_1991_2020) <- gsub(".+BIO_\\d\\d_(.+)_1991_2020\\.tiff", "\\1", files_recent)

BIOCLIM_1991_2020


#-------------------------------------------------------------------------
# create occurrence and non-occurrence records, hitch them to climate data

plot(sppOcc_composite) # starting from this as our "presence" polygon

# mask a BIOCLIM layer to the broad distribution, take the cell coordinates, pull 2000 at random
presence <- rasterize(st_as_sf(sppOcc_clean_recent, coords=c("decimallongitude", "decimallatitude"), crs=4326), BIOCLIM_1991_2020[[1]], fun=function(x) count(x)>0) %>% crds(df=TRUE) %>% rename(lon=x, lat=y)
glimpse(presence) # nice, there we go

glimpse(presence) # okay cool

# how do the points look?
ggplot() + geom_sf(data=ne_countries(scale=10, continent = "north america", returnclass = "sf")) + 
	geom_sf(data=sppOcc_composite, color=NA, fill="forestgreen", alpha=0.5) + 
	geom_point(data=presence, aes(x=lon, y=lat), size=0.5, alpha=0.1) +
	coord_sf(xlim = c(-125,-114), ylim = c(28,42), expand = FALSE) +
	theme(panel.background=element_rect(fill="slategray3", color=NA), panel.grid=element_line(color="slategray2"))

pseudabs <- mask(BIOCLIM_1991_2020[[1]], st_transform(sppOcc_composite, crs=crs(BIOCLIM_1991_2020[[1]]))) %>% 
	crds(df=TRUE) %>% rename(lon=x, lat=y) %>% 
	filter(!(paste(lon,lat) %in% paste(presence$lon, presence$lat))) %>% # remove cells represented in presence
	slice_sample(n=nrow(presence)) # pull a sample the same size as presence

glimpse(pseudabs)

# how do the points look?
ggplot() + geom_sf(data=ne_countries(scale=10, continent = "north america", returnclass = "sf")) + 
	geom_sf(data=sppOcc_composite, color=NA, fill="forestgreen", alpha=0.5) + 
	geom_point(data=presence, aes(x=lon, y=lat), size=0.5, alpha=0.1) +
	geom_point(data=pseudabs, aes(x=lon, y=lat), size=0.5, alpha=0.1, color="red") +
	coord_sf(xlim = c(-125,-114), ylim = c(28,42), expand = FALSE) +
	theme(panel.background=element_rect(fill="slategray3", color=NA), panel.grid=element_line(color="slategray2"))


#-------------------------------------------------------------------------
# link custom predictor values to presence/pseudoabsence points

PA_BIOCLIM_recent <- data.frame(pres=rep(1:0,each=nrow(presence)), rbind(presence, pseudabs), terra::extract(BIOCLIM_1991_2020, rbind(presence, pseudabs)))

glimpse(PA_BIOCLIM_recent)

write.table(PA_BIOCLIM_recent, paste0("output/presence_pseudoabsence_BIOCLIM_1991_2020_", taxon, ".csv"), sep=",", col.names=TRUE, row.names=FALSE)

# PA_BIOCLIM_recent <- read.csv(paste0("output/presence_pseudoabsence_BIOCLIM_1991_2020_", taxon, ".csv"))

#-------------------------------------------------------------------------
# Fit a species distribution model with BARTs

# Now we have presence and (pseudo)absence data, and values for a bunch of 
# climate variables at the corresponding locations. This is what we need to fit
# a species distribution model, predicting binary presence/absence with one or
# more environmental variables. Just as we model flowering activity with BARTs,
# we can model species presence with BARTs

library("embarcadero")

# Start with a vector of all our predictor names, for convenience
xnames <- colnames(PA_BIOCLIM_recent)[-1:-4]

# Then, predictor selection using varimp()
predSel <- varimp.diag(y.data=PA_BIOCLIM_recent[,"pres"], x.data=PA_BIOCLIM_recent[,xnames]) 
# [This may take a while. WAIT HERE to check in before moving on]

write_rds(predSel, file=paste0("output/BART/bart.SDM.varimp.", taxon, ".rds")) # save varimp() results
# predSel <- read_rds(paste0("output/BART/bart.SDM.varimp.", taxon, ".rds"))

# generate a better-organized varimp() figure
predSel$data <- predSel$data |> mutate(trees = factor(trees, c(10,20,50,100,150,200)))

predSel$labels$group <- "Trees"
predSel$labels$colour <- "Trees"

label_parse <- function(breaks){ parse(text=breaks) } # need this, for reasons

{cairo_pdf(paste("output/figures/SDM_varimp_", taxon, ".pdf", sep=""), width=6, height=5)

predSel + scale_x_discrete(label=label_parse) + 
theme_bw(base_size=12) +
theme(legend.position="inside", legend.position.inside=c(0.85, 0.7), axis.text.x=element_text(angle=45, hjust=1))  # okay nice

}
dev.off()

# Stepwise fitting is also an option if varimp() results aren't crystal-clear
sdm.step <- bart.step(y.data=as.numeric(PA[,"pres"]), x.data=PA[,xnames], full=FALSE, quiet=TRUE)


# okay, let's spell out the top predictors
topx <- c("TS", "ITH", "PS", "TAR", "PCoQ", "PWeQ", "PWeM")

# and train our working model
bart_sdm <- bart(y.train=as.numeric(PA_BIOCLIM_recent[,"pres"]), x.train=PA_BIOCLIM_recent[,topx], keeptrees=TRUE)

summary(bart_sdm) # LOL, LMAO

invisible(bart_sdm$fit$state)
write_rds(bart_sdm, file=paste0("output/BART/bart.SDM.", taxon, ".rds")) # save model
# bart_sdm <- read_rds(paste0("output/BART/bart.SDM.", taxon, ".rds"))



#-------------------------------------------------------------------------
# Now, predict species presence

# We don't want to do this worldwide, so mask the BioClim data to a region defined by presence data:
BIOCLIM_1991_2020_masked <- crop(mask(BIOCLIM_1991_2020, st_transform(broad.Range, crs=crs(BIOCLIM_1991_2020)), touches=TRUE), extent(st_transform(broad.Range, crs=crs(BIOCLIM_1991_2020))))
plot(BIOCLIM_1991_2020_masked[[1]])

bc_recent.df <- as.data.frame(BIOCLIM_1991_2020_masked)

BIOCLIM_1901_1930_masked <- crop(mask(BIOCLIM_1901_1930, st_transform(broad.Range, crs=crs(BIOCLIM_1991_2020)), touches=TRUE), extent(st_transform(broad.Range, crs=crs(BIOCLIM_1901_1930))))
plot(BIOCLIM_1901_1930_masked[[1]])

# for the recent period ...
# make a NaN vector of length ncell(stanRas)
pred_recent.vect <- rep(NaN, ncell(BIOCLIM_1991_2020_masked))
# insert prediction values at positions where stanRas has non-NaNs
pred_recent.vect[which(!is.na(values(BIOCLIM_1991_2020_masked[[1]])))] <- apply(stats::predict(object=bart_sdm, newdata=as.data.frame(BIOCLIM_1991_2020_masked)), 2, mean)
# convert it into a SpatRaster
pred_recent <- rast(ncol=ncol(BIOCLIM_1991_2020_masked), nrow=nrow(BIOCLIM_1991_2020_masked), crs=crs(BIOCLIM_1991_2020_masked), extent=ext(BIOCLIM_1991_2020_masked), vals=pred_recent.vect)

plot(pred_recent)

writeRaster(pred_recent, paste0("output/BART/SDM_", taxon, "_prediction_1991_2020.tiff"), overwrite=TRUE) 


# for the early period ...
# make a NaN vector of length ncell(stanRas)
pred_early.vect <- rep(NaN, ncell(BIOCLIM_1901_1930_masked))
# insert prediction values at positions where stanRas has non-NaNs
pred_early.vect[which(!is.na(values(BIOCLIM_1901_1930_masked[[1]])))] <- apply(stats::predict(object=bart_sdm, newdata=as.data.frame(BIOCLIM_1901_1930_masked)), 2, mean)
# convert it into a SpatRaster
pred_early <- rast(ncol=ncol(BIOCLIM_1901_1930_masked), nrow=nrow(BIOCLIM_1901_1930_masked), crs=crs(BIOCLIM_1901_1930_masked), extent=ext(BIOCLIM_1901_1930_masked), vals=pred_early.vect)

plot(pred_early)

writeRaster(pred_early, paste0("output/BART/SDM_", taxon, "_prediction_1901_1930.tiff"), overwrite=TRUE) 


# and, binary classification:
summary(bart_sdm)
cutoff <- 0.4489313
plot(pred_early >= cutoff)
plot(pred_recent >= cutoff) # hmm

writeRaster(pred_recent>=cutoff, paste0("output/BART/SDM_", taxon, "_prediction_binary_1991_2020.tiff"), overwrite=TRUE) 
writeRaster(pred_early>=cutoff, paste0("output/BART/SDM_", taxon, "_prediction_binary_1901_1930.tiff"), overwrite=TRUE) 

length(which(values(pred_early)>=cutoff))
length(which(values(pred_recent)>=cutoff))

hist(values(pred_recent-pred_early))
quantile(values(pred_recent-pred_early), c(0.025, 0.5, 0.975), na.rm=TRUE) # median = 0.11, 95%CI -0.30 to 0.28
t.test(values(pred_recent-pred_early)) # mean > 0, p < 2.2e-16


