# Using Continuous-response DARTs to predict historical flowering in a target species
# last used/modified jby, 2025.02.05

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/flower_prediction")

library("tidyverse")

library("embarcadero")
library("SoftBart")

library("terra")

#-----------------------------------------------------------
# initial file loading

# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia = 57250

flow <- read.csv(paste("output/flowering_freq_climate_", taxon, ".csv", sep="")) %>% filter(!is.na(ppt.y0q1)) # flowering/not flowering, biologically-informed candidate predictors

dim(flow)
glimpse(flow)

table(flow$year, flow$prop_flr>0)
hist(flow$prop_flr)

SppExt <- round(c(range(flow$lon), range(flow$lat)) * c(1.05,0.95,0.95,1.05),0) # for toyon; need to adjust accordingly



#-----------------------------------------------------------
# build species predictions from historical PRISM layers
# computation time is determined by the time span covered and the area of prediction; 
# this is one reason to try to restrict the size of the cropped area of consideration!

if(!dir.exists(paste("output/models/DART_predictions.", taxon, sep=""))) dir.create(paste("output/models/DART_predictions.", taxon, sep=""), recursive = TRUE)

# load the saved model developed in `phenology_modeling.R`
flr.mod <- read_rds(file=paste0("output/models/flrDARTtop_", taxon, ".rds")) # swap in RI if needed
summary(flr.mod)

flower.preds <- strsplit(gsub(".+~ (.+)", "\\1", flr.mod$formula), split="\\+")[[1]]
flower.preds

# LOOP over years ---------------------------------------------------
for(yr in 1900:2024){ # adjust year range based on data available

# yr <- 2015

stanRas <- rast(paste0("../data/PRISM/quarterlies/ppt_cropped_",yr,"Q1.bil"))

# Building a brick with all possible predictors
preds <- c(
	resample(rast(paste0("../data/PRISM/quarterlies/ppt_cropped_",yr,"Q1.bil")), stanRas, method="near"),
	resample(rast(paste0("../data/PRISM/quarterlies/tmax_cropped_",yr,"Q1.bil")), stanRas, method="near"), 	
	resample(rast(paste0("../data/PRISM/quarterlies/tmin_cropped_",yr,"Q1.bil")), stanRas, method="near"), 	
	resample(rast(paste0("../data/PRISM/quarterlies/vpdmax_cropped_",yr,"Q1.bil")), stanRas, method="near"), 	
	resample(rast(paste0("../data/PRISM/quarterlies/vpdmin_cropped_",yr,"Q1.bil")), stanRas, method="near"), 	
	resample(rast(paste0("../data/PRISM/quarterlies/ppt_cropped_",yr-1,"Q4.bil")), stanRas, method="near"),
	resample(rast(paste0("../data/PRISM/quarterlies/tmax_cropped_",yr-1,"Q4.bil")), stanRas, method="near"), 	
	resample(rast(paste0("../data/PRISM/quarterlies/tmin_cropped_",yr-1,"Q4.bil")), stanRas, method="near"), 	
	resample(rast(paste0("../data/PRISM/quarterlies/vpdmax_cropped_",yr-1,"Q4.bil")), stanRas, method="near"), 	
	resample(rast(paste0("../data/PRISM/quarterlies/vpdmin_cropped_",yr-1,"Q4.bil")), stanRas, method="near"), 	
	resample(rast(paste0("../data/PRISM/quarterlies/ppt_cropped_",yr-1,"Q3.bil")), stanRas, method="near"),
	resample(rast(paste0("../data/PRISM/quarterlies/tmax_cropped_",yr-1,"Q3.bil")), stanRas, method="near"), 	
	resample(rast(paste0("../data/PRISM/quarterlies/tmin_cropped_",yr-1,"Q3.bil")), stanRas, method="near"), 	
	resample(rast(paste0("../data/PRISM/quarterlies/vpdmax_cropped_",yr-1,"Q3.bil")), stanRas, method="near"), 	
	resample(rast(paste0("../data/PRISM/quarterlies/vpdmin_cropped_",yr-1,"Q3.bil")), stanRas, method="near") 	
	)
names(preds) <- c("ppt.y0q1", "tmax.y0q1", "tmin.y0q1", "vpdmax.y0q1", "vpdmin.y0q1", "ppt.y1q4", "tmax.y1q4", "tmin.y1q4", "vpdmax.y1q4", "vpdmin.y1q4", "ppt.y1q3", "tmax.y1q3", "tmin.y1q3", "vpdmax.y1q3", "vpdmin.y1q3")

preds.crop <- crop(preds, SppExt)

# predictor dataframe from raster layers
preds.df <- as.data.frame(preds.crop[[flower.preds]])

# make a NaN vector of length ncell(stanRas)
pred.vect.mu_mean <- rep(NaN, ncell(preds.crop[[1]]))
# insert prediction values at positions where stanRas has non-NaNs
pred.vect.mu_mean[which(!is.nan(values(preds.crop[[1]])))] <- predict(object=flr.mod, newdata=data.frame(prop_flr=0, preds.df))$mu_mean

# convert it into a SpatRaster
pred.rast.mu_mean <- rast(ncol=ncol(preds.crop[[1]]), nrow=nrow(preds.crop[[1]]), crs=crs(preds.crop[[1]]), extent=ext(preds.crop[[1]]), vals=pred.vect.mu_mean)

# pred.rast.mu_mean
# plot(pred.rast.mu_mean) # helllll yeah

# and write out the layer
writeRaster(pred.rast.mu_mean, paste0("output/models/DART_predictions.", taxon, "/DART_predicted_flowering_", taxon, "_", yr, ".tiff"), overwrite=TRUE) # got to figure this out

cat("DONE with predictions for", yr, "\n")

} # END loop over years



