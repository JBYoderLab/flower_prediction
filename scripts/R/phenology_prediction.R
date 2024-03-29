# Using BARTs to predict historical flowering in a target species
# last used/modified jby, 2024.03.28

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/flower_prediction")

library("tidyverse")

library("embarcadero")


#-----------------------------------------------------------
# initial file loading


# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia = 57250

flow <- read.csv(paste("output/flowering_obs_climate_", taxon, ".csv", sep="")) %>% filter(!is.na(ppt.y0q1)) # flowering/not flowering, biologically-informed candidate predictors, subspecies id'd

dim(flow)
glimpse(flow)


#-----------------------------------------------------------
# build species predictions from historical PRISM layers
# computation time is determined by the time span covered and the area of prediction; 
# this is one reason to try to restrict the size of the cropped area of consideration!

if(!dir.exists(paste("output/BART/predictions.", taxon, sep=""))) dir.create(paste("output/BART/predictions.", taxon, sep=""))

# load the saved model developed in `phenology_modeling.R`
flr.mod <- read_rds(file=paste("output/BART/bart.model.", taxon, ".rds", sep="")) # swap in RI if needed
summary(flr.mod)

flower.preds <- attr(flr.mod$fit$data@x, "term.labels")
flower.preds

# LOOP over years
for(yr in 1900:2023){ # adjust year range based on data available

# yr <- 2015

stanRas <- raster(paste("data/PRISM/annual.", taxon, "/ppt_cropped_",yr,"Q1.bil", sep=""))

# Building a brick with all possible predictors
preds <- brick(
	resample(raster(paste("data/PRISM/annual.", taxon, "/ppt_cropped_",yr,"Q1.bil", sep="")), stanRas, method="ngb"),
	resample(raster(paste("data/PRISM/annual.", taxon, "/tmax_cropped_",yr,"Q1.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual.", taxon, "/tmin_cropped_",yr,"Q1.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual.", taxon, "/vpdmax_cropped_",yr,"Q1.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual.", taxon, "/vpdmin_cropped_",yr,"Q1.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual.", taxon, "/ppt_cropped_",yr-1,"Q4.bil", sep="")), stanRas, method="ngb"),
	resample(raster(paste("data/PRISM/annual.", taxon, "/tmax_cropped_",yr-1,"Q4.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual.", taxon, "/tmin_cropped_",yr-1,"Q4.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual.", taxon, "/vpdmax_cropped_",yr-1,"Q4.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual.", taxon, "/vpdmin_cropped_",yr-1,"Q4.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual.", taxon, "/ppt_cropped_",yr-1,"Q3.bil", sep="")), stanRas, method="ngb"),
	resample(raster(paste("data/PRISM/annual.", taxon, "/tmax_cropped_",yr-1,"Q3.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual.", taxon, "/tmin_cropped_",yr-1,"Q3.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual.", taxon, "/vpdmax_cropped_",yr-1,"Q3.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual.", taxon, "/vpdmin_cropped_",yr-1,"Q3.bil", sep="")), stanRas, method="ngb") 	
	)
names(preds) <- c("ppt.y0q1", "tmax.y0q1", "tmin.y0q1", "vpdmax.y0q1", "vpdmin.y0q1", "ppt.y1q4", "tmax.y1q4", "tmin.y1q4", "vpdmax.y1q4", "vpdmin.y1q4", "ppt.y1q3", "tmax.y1q3", "tmin.y1q3", "vpdmax.y1q3", "vpdmin.y1q3")


# prediction with the RI predictor (year) removed
pred.yr <- predict(flr.mod, preds[[flower.preds]], splitby=20)

# pred.yr # useful confirmation if you're doing a single test year

writeRaster(pred.yr, paste("output/BART/predictions.", taxon, "/BART_predicted_flowering_", taxon, "_", yr, ".bil", sep=""), overwrite=TRUE)

} # END loop over years



