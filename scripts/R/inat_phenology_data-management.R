# working with phenology-annotated iNat observations
# Assumes MAJEL environment 
# jby 2024.03.26

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flower_prediction")

library("tidyverse")
library("lubridate")

library("raster")
library("sf")

# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia = 57250

# Species area crop extent (deliberately generous)
SppExt <- extent(-125, -109, 23, 42) # for toyon; need to adjust accordingly


#-------------------------------------------------------------------------
# read in iNat observations compiled using inat_phenology_obs.R

inat <- read.csv(paste("data/inat_phenology_data_", taxon, ".csv", sep=""), h=TRUE) %>% mutate(observed_on = ymd(observed_on))

glimpse(inat) # how many raw observations?
table(inat$phenology) # by phenophase
table(inat$phenology, inat$year) # by phenophase and year

# is fruit ever observed in Jan, Feb, or March? That's really (maybe?) the PREVIOUS flowering year
# Take (again) the toyon example: per iNat, peak flowering is in ~June, peak fruiting in Nov
# To account for this we need to see how often fruit is logged BEFORE peak flowering
filter(inat, phenology=="Fruiting", month(observed_on)<6) # if this is not zero, we need to deal with these!

# create a "flowering year" variable that wraps early-year observations of fruit into the previous year
inat$flr_yr <- inat$year
inat$flr_yr[inat$phenology=="Fruiting" & month(inat$observed_on)<6] <- inat$year[inat$phenology=="Fruiting" & month(inat$observed_on)<6]-1

#-------------------------------------------------------------------------
# organize iNat observations for extraction of summarized PRISM data

# data structure setup
flowering <- data.frame(matrix(0,0,4))
names(flowering) <- c("lon","lat","year", "flr")

prism_temp_rast <- raster(paste("data/PRISM/annual.", taxon, "/tmax_cropped_2010Q1.bil", sep="")) # raster grid base

# then LOOP over years in the raw data ....
for(yr in unique(inat$flr_yr)){

# yr <- 2012

if(length(which(inat$flr_yr==yr & inat$phenology!="No Evidence of Flowering"))>0){
	yes <- rasterize(dplyr::filter(inat, flr_yr==yr, phenology!="No Evidence of Flowering")[,c("longitude","latitude")], prism_temp_rast, fun=sum, background=0)
	yearyes <- rasterToPoints(yes, fun=function(x){x>=1})
	outyes <-data.frame(lon=yearyes[,"x"], lat=yearyes[,"y"], year=yr, flr=TRUE)
}else{
	outyes <- NULL
}

if(length(which(inat$flr_yr==yr & inat$phenology=="No Evidence of Flowering"))>0){
	no <- rasterize(dplyr::filter(inat, flr_yr==yr, phenology=="No Evidence of Flowering")[,c("longitude","latitude")], prism_temp_rast, fun=sum, background=0)
	yearno <- rasterToPoints(no, fun=function(x){x>=1})
	outno <- data.frame(lon=yearno[,"x"], lat=yearno[,"y"], year=yr, flr=FALSE)	
}else{
	outno <- NULL
}

# put it all together
flowering <- rbind(flowering, outyes, outno)
cat("Done with year", yr, "\n")

} # END LOOP over years

head(flowering)
glimpse(flowering) # okay okay okay!
table(flowering$year, flowering$flr) # smiling serenely

write.table(flowering, paste("output/flowering_obs_rasterized_", taxon,".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE)


#-------------------------------------------------------------------------
# attach PRISM data to flowering/not flowering observations

# we've got variables aggregated quarterly, assume we want Year 0 (year of obs) and Year -1 (year before obs)
varnames <- paste(rep(c("ppt", "tmax", "tmin", "vpdmax", "vpdmin"), each=4), paste(rep(c("y0", "y1"), each=20), paste("q", 1:4, sep=""), sep=""), sep=".")

# new data structure
flr.clim <- data.frame(matrix(0,0,ncol(flowering)+length(varnames)))
names(flr.clim) <- c(colnames(flowering), varnames)

# LOOP over years, because of the current year previous year thing ...
for(yr in sort(unique(flowering$year))){

# yr <- 2020 # test condition

# specific years to target
y0 <- yr
y1 <- yr-1

# list of files ... this assumes known variables, I should consider rewriting to use list.files()
prismfiles <- paste("data/PRISM/annual.", taxon, "/", rep(c("ppt", "tmax", "tmin", "vpdmax", "vpdmin"), each=4), "_cropped_", rep(c(y0, y1), each=20), "Q", 1:4, ".bil", sep="")

# assemble weather data predictors for a given year
preds <- stack(lapply(prismfiles, function(x) crop(raster(x), SppExt)))
names(preds) <- varnames

# pull subset of flowering observations for year
flsub <- subset(flowering, year==yr)
flsub <- cbind(flsub, raster::extract(preds, flsub[,c("lon","lat")], df=FALSE))

flr.clim <- rbind(flr.clim,flsub) 

write.table(flr.clim, paste("output/flowering_obs_climate_", taxon, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE)

cat("Done with year", yr, "\n")


} # END LOOP over years

glimpse(flr.clim)


# and that's generated a data file we can feed into embarcadero ... in the next script!








