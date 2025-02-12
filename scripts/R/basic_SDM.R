# Build a species distribution polygon for downstream analyses
# Assumes local environment
# jby 2025.02.08

# setwd("~/Documents/Active_projects/flower_prediction")

# Clear the environment and load key packages
rm(list=ls())

library("tidyverse") 
library("raster")
library("terra")
library("geodata")
library("sf")
library("embarcadero")

# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia: 57250

#-----------------------------------------------------------
# Load and manage input data: occurrences and environmental values

# rasterized iNat records to use as occurrence data
inat <- read.csv(paste("data/inat_phenology_data_", taxon,"_cleaned.csv", sep=""))
glimpse(inat)

# Species area crop extent (deliberately generous)
# note the padding factors assume we're NW of 0,0
SppExt <- round(c(range(inat$longitude), range(inat$latitude)) * c(1.01,0.99,0.9,1.1),0) # for toyon; need to adjust accordingly


# BioClim normals for 19 standard "biologically relevant" climate variables
BClim <- worldclim_tile(var="bio", res=0.5, path="../data/spatial", lon=)

# BClim is a raster data object, a RasterStack, with a "layer" for each of the
# 19 BioClim variables.
BClim

# We can rename the raster layers for improved understanding:
names(BClim) <- c("MAT", "MDR", "ITH", "TS", "MWaMT", "MCMT", "TAR", "MWeQT", "MDQT", "MWaQT", "MCQT", "AP", "PWeM", "PDM", "PS", "PWeQ", "PDQ", "PWaQ", "PCQ")

# These are explained in detail here:
# https://www.worldclim.org/data/bioclim.html

# And we can plot individual layers of the RasterStack thus:
plot(BClim, "MAT", cex=0.5, legend=T, mar=par("mar"), 
     xaxt="n", yaxt="n", main="Mean annual temperature (ºC x 10)")


#-----------------------------------------------------------
# Create non-occurrence records

# SDMs need observations of places where a species is present AND places where 
# it isn't present --- the contrast between conditions at these places is the
# key focus of an SDM. Often, though, we don't have formal records of places
# where a species is absent! So we create pseudo-absences to contrast with our
# presences. How we create these is a bit tricky.

# First, your occurrence records are probably unevenly distributed in space.
# We can fix this (a bit) by converting them to a raster at the same resoluiton
# as our climate data, then converting the raster back to points.
occ_rast <- rasterize(as.matrix(inat[,c("longitude", "latitude")]), BClim[[1]], fun=sum, background=0)
occ_thin <- as.data.frame(rasterToPoints(occ_rast, function(x) x > 0)[,1:2])
colnames(occ_thin) <- c("lon", "lat")
glimpse(occ_thin)
# This is one presence point for each grid cell containing ANY points in the 
# original data set.

# Now, we need to generate pseudo-absences. First, we define regions that are
# not so close to our presence points that they probably are effectively also
# habitat for our focal species, but not so far that they're beyond where the 
# species could reasonably migrate to. This is a VERY subjective judgement, and
# often based on general estimates of dispersal distance or on-the-ground 
# experience. 
occ_sf <- st_as_sf(occ_thin, coords=c("lon", "lat"), crs=4326) # coordinates are in degrees
occ_sf <- st_transform(occ_sf, crs=3857) # units in meters

x.inner <- st_union(st_buffer(occ_sf, 10000)) # polygons based on 10km radii around occ_thin
x.outer <- st_union(st_buffer(occ_sf, 50000)) # polygons based on 50km radii
x.donuts <- st_difference(x.outer, x.inner) # difference between above

ggplot() + geom_sf(data=x.donuts) 
# A conservative region for pseudoabsences: Not within 5km of a presence point,
# but not more than 50km from a presence point. You may want to adjust this!

# Now, we draw random points within the "donut" regions --- the same number as 
# we have presence points, for balance
pseudabs <- st_sample(x.donuts, nrow(occ_thin))

ggplot() + 
  geom_sf(data=occ_sf, color="blue", size=1) + 
  geom_sf(data=pseudabs, color="red", size=0.5)
# There's our presence points and pseudo-absences, visualized in space.

#-------------------------------------------------------------------------
# 4. Collect environmental values for our presence and absence locations

pres_envs <- data.frame(raster::extract(BClim, occ_thin[,c("lon","lat")]))
glimpse(pres_envs) # BioClim values for the presence locations

pa_sp <- as_Spatial(st_transform(pseudabs, crs=4326))
pa_envs <- data.frame(raster::extract(BClim, pa_sp))
glimpse(pa_envs) # And ditto for the pseudoabsences

# Now, assemble a single data frame with presence/absence coordinates and the
# associated BioClim values:

PA <- rbind(data.frame(lon=occ_thin$lon, lat=occ_thin$lat, pres=1, pres_envs), 
          data.frame(lon=coordinates(pa_sp)[,1], lat=coordinates(pa_sp)[,2], pres=0,pa_envs)) %>% 
          mutate(id = row_number()) %>% 
          filter(!is.na(MAT)) # if your taxon is on the coast, you may have presences OR pseudo-absences without data

glimpse(PA)
dim(PA) # should be twice the size of the presence data set (barring removals for missing BioClim data)

# You may want to write out this data frame for later /read back in 
write.table(PA, paste("data/presence-pseudoabsence_",taxon, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

# PA <- rpaste("data/presence-pseudoabsence_",taxon, ".csv")
glimpse(PA)

#-------------------------------------------------------------------------
# Fit a species distribution model with BARTs

# Now we have presence and (pseudo)absence data, and values for a bunch of 
# climate variables at the corresponding locations. This is what we need to fit
# a species distribution model, predicting binary presence/absence with one or
# more environmental variables. Just as we model flowering activity with BARTs,
# we can model species presence with BARTs

library("embarcadero")

# Start with a vector of all our predictor names, for convenience
xnames <- c("MAT", "MDR", "ITH", "TS", "MWaMT", "MCMT", "TAR", "MWeQT", "MDQT", "MWaQT", "MCQT", "AP", "PWeM", "PDM", "PS", "PWeQ", "PDQ", "PWaQ", "PCQ")

# Then, predictor selection using varimp()
predSel <- varimp.diag(y.data=PA[,"pres"], x.data=PA[,xnames]) 
# [This may take a while. WAIT HERE to check in before moving on]

# generate a better-organized varimp() figure
predSel$data <- predSel$data |> mutate(trees = factor(trees, c(10,20,50,100,150,200)))

predSel$labels$group <- "Trees"
predSel$labels$colour <- "Trees"

label_parse <- function(breaks){ parse(text=breaks) } # need this, for reasons

{cairo_pdf(paste("output/figures/SDM_varimp_", taxon, ".pdf", sep=""), width=6, height=5)

predSel + scale_x_discrete(label=label_parse) + 
theme_bw(base_size=12) +
theme(legend.position="inside", legend.position.inside=c(0.8, 0.7), axis.text.x=element_text(angle=45, hjust=1))  # okay nice

}
dev.off()

# Stepwise fitting is also an option if varimp() results aren't crystal-clear
sdm.step <- bart.step(y.data=as.numeric(PA[,"pres"]), x.data=PA[,xnames], full=FALSE, quiet=TRUE)


#-------------------------------------------------------------------------
# Now, predict species presence

# We don't want to do this worldwide, so mask the BioClim data to a region defined by presence data:
bcMask <- st_union(st_buffer(occ_sf, 100000)) # polygons based on 100km radii

BC.mask <- mask(crop(BClim, SppExt), st_as_sf(st_transform(bcMask, crs=4269)), touches=TRUE)
plot(BC.mask[[1]]) # see what we have, to check

# predict species presence from the SDM of your choice in the masked layers
pred <- predict(sdm.step, BC.mask) 

# The output is another raster layer, with probability of your species' 
# occurrence in each cell:
plot(pred) # How does that look?

# from here, we need to threshold the predicted presence to create a range polygon
# we get the best-power threshold from the summary() of our model
summary(sdm.step)
cutoff <- 0.5967603

plot(pred>cutoff) # how's that look?

writeRaster(pred, paste("output/BART_SDM_",taxon,".grd", sep="")) # write out the modeled presence probability


# building new distribution polys ...
predPts <- rasterToPoints(pred, function(x) x>=cutoff)
# predPts
plot(predPts, col="black", pch=16, cex=0.3) # okay

class(predPts)
colnames(predPts) <- c("lon", "lat", "prob")

presPts.sf <- st_as_sf(data.frame(predPts), coords=c("lon", "lat"), crs=4326) # coordinates are in degrees

presPts.sf <- st_transform(presPts.sf, crs=3857) # units in meters

presPts.10k <- st_union(st_buffer(presPts.sf, 10000)) # polygons based on 10km radii around bkgd obs

plot(presPts.10k) # I don't hate this?

pres.poly <- st_as_sf(st_transform(presPts.10k, crs=4326)) # this seems right?

plot(pres.poly) # and this works

if(!dir.exists(paste("output/BART_SDM_", taxon, sep=""))) dir.create(paste("output/BART_SDM_", taxon, sep=""))
st_write(pres.poly, paste("output/BART_SDM_", taxon, "/BART_SDM_range_poly_", taxon, ".shp", sep=""), append=FALSE) # and write it out


