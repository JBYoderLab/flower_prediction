# Using BARTs to predict historical flowering in a target species
# last used/modified jby, 2023.11.22

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/flowering_prediction")

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
# build species predictions from historical PRISM (let's do 1900-2022)


# flower, with RI --------------------------------
if(!dir.exists("output/BART/RI.predictions")) dir.create("output/BART/RI.predictions")

# load the saved model developed in `phenology_modeling.R`
flower.RImod <- read_rds(file=paste("output/BART/bart.RImodel.", taxon, ".rds", sep="")) # this now WORKS

flower.preds <- attr(flower.RImod$fit[[1]]$data@x, "term.labels")

# LOOP over years
for(yr in 2010:2022){ # adjust year range based on data available

# yr <- 2015

stanRas <- raster(paste("data/PRISM/annual/ppt_cropped_",yr,"Q1.bil", sep=""))

# Building a brick with all possible predictors
preds <- brick(
	resample(raster(paste("data/PRISM/annual/ppt_cropped_",yr,"Q1.bil", sep="")), stanRas, method="ngb"),
	resample(raster(paste("data/PRISM/annual/tmax_cropped_",yr,"Q1.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual/tmin_cropped_",yr,"Q1.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual/vpdmax_cropped_",yr,"Q1.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual/vpdmin_cropped_",yr,"Q1.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual/ppt_cropped_",yr-1,"Q4.bil", sep="")), stanRas, method="ngb"),
	resample(raster(paste("data/PRISM/annual/tmax_cropped_",yr-1,"Q4.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual/tmin_cropped_",yr-1,"Q4.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual/vpdmax_cropped_",yr-1,"Q4.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual/vpdmin_cropped_",yr-1,"Q4.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual/ppt_cropped_",yr-1,"Q3.bil", sep="")), stanRas, method="ngb"),
	resample(raster(paste("data/PRISM/annual/tmax_cropped_",yr-1,"Q3.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual/tmin_cropped_",yr-1,"Q3.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual/vpdmax_cropped_",yr-1,"Q3.bil", sep="")), stanRas, method="ngb"), 	
	resample(raster(paste("data/PRISM/annual/vpdmin_cropped_",yr-1,"Q3.bil", sep="")), stanRas, method="ngb") 	
	)
names(preds) <- c("ppt.y0q1", "tmax.y0q1", "tmin.y0q1", "vpdmax.y0q1", "vpdmin.y0q1", "ppt.y1q4", "tmax.y1q4", "tmin.y1q4", "vpdmax.y1q4", "vpdmin.y1q4", "ppt.y1q3", "tmax.y1q3", "tmin.y1q3", "vpdmax.y1q3", "vpdmin.y1q3")


# prediction with the RI predictor (year) removed
pred.ri0 <- predict(flower.RImod, preds[[attr(flower.RImod$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/RI.predictions/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

} # END loop over years



