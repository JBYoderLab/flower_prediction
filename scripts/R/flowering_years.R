# Analyzing predicted historical flowering 
# jby 2024.03.28

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flower_prediction")

# Clear the environment and load key packages
rm(list=ls())

library("tidyverse")

library("raster")
library("sp")
library("sf")
library("hexbin")
library("embarcadero")

#-------------------------------------------------------------------------
# Load and prep data

# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia = 57250

sdm.pres <- read_sf(paste("output/BART_SDM_", taxon, "/BART_SDM_range_poly_", taxon, ".shp", sep="")) 

# flowering observation data
obs <- read.csv(paste("output/flowering_obs_climate_", taxon, ".csv", sep=""))
# raster files of predicted prFL
pred.files <- list.files(paste("output/BART/predictions.", taxon, sep=""), pattern=".bil", full=TRUE)

# useful bits and bobs
# Species area crop extent (deliberately generous)
# note the padding factors assume we're NW of 0,0
SppExt <- round(c(range(obs$lon), range(obs$lat)) * c(1.01,0.99,0.9,1.1),0) # for toyon; need to adjust accordingly

sdm.buff <- st_buffer(st_transform(sdm.pres[,2], crs=3857), 1000) # put a 1km buffer on the range polygons

#-------------------------------------------------------------------------
# Process prediction layers

pred.histStack <- raster::stack(sapply(pred.files, function(x) crop(raster::raster(x), SppExt))) # read in prediction layers, crop to the species extent to reduce data volume later
names(pred.histStack) <- paste("prFL", 1900:2023, sep=".") # name the layers by year
projection(pred.histStack) <- CRS("+init=epsg:4269") # add CRS

pred.histStack # what's this look like

writeRaster(pred.histStack, paste("output/BART/BART_predicted_flowering_", taxon, "_1900-2023_nomask.grd", sep=""), overwrite=TRUE) # write out the cropped stack as a single raster file

# mask to the SDM polygon
pred.maskHist <- mask(pred.histStack, st_transform(sdm.buff, crs=4269), touches=TRUE)

writeRaster(pred.maskHist, paste("output/BART/BART_predicted_flowering_", taxon, "_1900-2023.grd", sep=""), overwrite=TRUE) # write out the masked layers as a single raster file

# build a data frame for analysis
hist.flowering <- cbind(coordinates(pred.maskHist), as.data.frame(pred.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(hist.flowering)


#-------------------------------------------------------------------------
# Reconstructed flowering years for 1900-2023

# get the classification cutoff that maximizes TSS
flr.mod <- read_rds(paste("output/BART/BART.model.", taxon, ".rds", sep=""))
summary(flr.mod) # for toyon example, this is 0.79

cutoff <- 0.7929546

flyrs <- hist.flowering %>% group_by(lon,lat) %>% 
summarize(
	flyrs_all=length(which(prFL>=cutoff)), 
	flyrs_1900_1929=length(which(prFL>=cutoff & year<=1929)), # earliest three decades
	flyrs_1990_2019=length(which(prFL>=cutoff & year>=1990 & year<=2019)) # most recent three-decade period
	) |> mutate(
	flyrs_change=flyrs_1990_2019-flyrs_1900_1929
	) # rangewide model
range(hist.flowering$year) # should be: 1900-2023
glimpse(flyrs)

write.table(flyrs, paste("output/reconstructed_flowering_years_", taxon, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE)
# flyrs <- read.csv(paste("output/reconstructed_flowering_years_", taxon, ".csv", sep=""))

# Inspect the flowering years .. does this make biological sense?
quantile(flyrs$flyrs_all, c(0.025,0.5,0.975)) 
quantile(flyrs$flyrs_all, c(0.025,0.5,0.975))/124 

quantile(flyrs$flyrs_1900_1929, c(0.025,0.5,0.975))/30
quantile(flyrs$flyrs_1990_2019, c(0.025,0.5,0.975))/30

quantile(flyrs$flyrs_change, c(0.025,0.5,0.975)) 
mean(flyrs$flyrs_change) 


#-------------------------------------------------------------------------
# visualizations


# histograms
ggplot(flyrs, aes(x=flyrs_all)) + geom_histogram(bins=30)
ggplot(flyrs, aes(x=flyrs_1900_1929)) + geom_histogram(bins=30)
ggplot(flyrs, aes(x=flyrs_1990_2019)) + geom_histogram(bins=30)

# distributions of flowering years
{cairo_pdf(paste("output/figures/flowering_years_", taxon, ".pdf", sep=""), width=3.5, height=2.5)

ggplot(flyrs, aes(x=flyrs_all)) + geom_histogram(fill="#b2df8a") + labs(x="Projected flowering years, 1900-2023", y="Grid cells") + 

geom_vline(xintercept=median(flyrs$flyrs_all), color="#33a02c") +

theme_bw() + theme(legend.position="none", plot.margin=margin(0.1,0.15,0.1,0.15, "inches"))

}
dev.off()


# MAPS ----------------------------------------------------

library("rnaturalearth")
library("rnaturalearthdata")

# map elements
sdm.pres <- read_sf(paste("output/BART_SDM_", taxon, "/BART_SDM_range_poly_", taxon, ".shp", sep="")) 

states <- ne_states(country="united states of america", returnclass="sf")
countries <- ne_countries(scale=10, continent="north america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

flrfrq_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
	geom_sf(data=countries, fill="antiquewhite3", color="antiquewhite4") + 
	geom_sf(data=states, fill="antiquewhite2", color="antiquewhite4") + 
	
	geom_tile(data=flyrs, aes(x=lon, y=lat, fill=flyrs_all/124)) + 
	
	scale_fill_distiller(type="seq", palette="Greens", direction=1, name="Flowering frequency,\n1900-2023", breaks=c(0,0.25,0.5)) + labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = SppExt[1:2], ylim = SppExt[3:4], expand = FALSE) +
	
	theme_minimal(base_size=9) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.1, "in"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=8), legend.title=element_text(size=9), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())


{cairo_pdf(paste("output/figures/flfrq_map_", taxon, ".pdf", sep=""), width=5, height=5)

flrfrq_map

}
dev.off()

flrfrq_change <- 




