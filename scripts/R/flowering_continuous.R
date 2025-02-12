# Analyzing predicted historical flowering 
# jby 2025.02.11

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flower_prediction")

# Clear the environment and load key packages
rm(list=ls())

library("tidyverse")

library("terra")
library("geodata")
library("sp")
library("sf")

library("embarcadero")

#-------------------------------------------------------------------------
# Load and prep data

# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia = 57250

# flowering observation data
obs <- read.csv(paste("output/flowering_freq_climate_", taxon, ".csv", sep=""))
# raster files of predicted prFL
pred.files <- list.files(paste("output/models/DART_predictions.", taxon, sep=""), pattern=".tiff", full=TRUE)

# useful bits and bobs
# Species area crop extent (deliberately generous)
# note the padding factors assume we're NW of 0,0
SppExt <- round(c(range(obs$lon), range(obs$lat)) * c(1.01,0.99,0.9,1.1),0) # for toyon; need to adjust accordingly

# species distribution polygon (derive from SDM or use prior from USFS)
spp.Range <- read_sf("../data/spatial/wpetry-USTreeAtlas-4999258/shp/photarbu/photarbu.shp", crs=4267)
spp.Range
spp.Range.buff <- st_transform(st_buffer(st_transform(spp.Range, crs=3857), 10000), crs=4267) # 100km buffer?

#-------------------------------------------------------------------------
# Process prediction layers

historic.preds <- sapply(pred.files, simplify=TRUE, function(x) crop(rast(x), SppExt)) # read in prediction layers, crop to the species extent to reduce data volume later
historic.stack <- rast(historic.preds)
names(historic.stack) <- paste("prFL", 1900:2023, sep=".") # name the layers by year

historic.stack # what's this look like

writeRaster(historic.stack, paste("output/models/DART_predicted_flowering_", taxon, "_1900-2023_nomask.tiff", sep=""), overwrite=TRUE) # write out the cropped stack as a single raster file

# mask to the SDM polygon
historic.stack.masked <- mask(historic.stack, st_transform(spp.Range.buff, crs=4326), touches=TRUE)

writeRaster(historic.stack.masked, paste0("output/models/DART_predicted_flowering_", taxon, "_1900-2023_masked.tiff"), overwrite=TRUE) # write out the masked layers as a single raster file

# read back in
historic.stack.masked <- rast(paste0("output/models/DART_predicted_flowering_", taxon, "_1900-2023.tiff"))

#-------------------------------------------------------------------------
# pair with relevant info, reformat for analysis

elev <- rast("../data/spatial/SR_50M/SR_50M.tif")
elev

# build a data frame for analysis? 
historic.flowering <- cbind(crds(historic.stack.masked, df=TRUE), as.data.frame(historic.stack.masked)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(historic.flowering) # oh that's a lot

write.table(historic.flowering, paste0("output/DART_predicted_flowering_", taxon, "_1900-2023.csv"), sep=",")

# historic.flowering <- read.csv(paste0("output/DART_predicted_flowering_", taxon, "_1900-2023.csv"))

#-------------------------------------------------------------------------
# Trend in flowering frequency, 1900-2023

# correlation coefficients
timecor <- historic.flowering %>% group_by(lat, lon) %>% do(broom::tidy(cor.test(~prFL+year, data=., method="spearman"))) %>% ungroup()

timecor$elev_m <- terra::extract(elev, timecor[,c("lon", "lat")])$SR_50M

glimpse(timecor)

hist(timecor$estimate) # interesting
quantile(timecor$estimate, c(0.025, 0.5, 0.975))
t.test(timecor$estimate) # mean > 0, p < 2.2e-16

cor.test(~estimate+lat, data=timecor, method="pearson")
# cor = 0.07, p = 9.1e-10 --- greater positive trend farther north

cor.test(~estimate+elev_m, data=timecor, method="pearson")
# cor = 0.07, p = 4.5e-09 --- greater positive trend farther north


# map this ------------------------------------------------
library("rnaturalearth")
library("rnaturalearthdata")

# map elements
states <- ne_states(country="united states of america", returnclass="sf")
countries <- ne_countries(scale=10, continent="north america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")


flrcor_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=timecor, aes(x=lon, y=lat, fill=estimate)) + 
	
	geom_sf(data=spp.Range, fill=NA, color="black", linewidth=0.5, linetype=2) + 
	
	scale_fill_gradient2(low="#762a83", mid="white", high="#1b7837", name="Trend in flowering\nintensity, 1900-2023", breaks=c(-0.25,0,0.25)) + labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125,-114), ylim = c(32,42), expand = FALSE) +
	
	theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.4, "inches"), legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())


{cairo_pdf(paste("output/figures/DART_flowering_trend_map_", taxon, ".pdf", sep=""), width=4, height=5.25)

flrcor_map

}
dev.off()



# barbar
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

flrfrq_change <- <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
	geom_sf(data=countries, fill="antiquewhite3", color="antiquewhite4") + 
	geom_sf(data=states, fill="antiquewhite2", color="antiquewhite4") + 
	
	geom_tile(data=flyrs, aes(x=lon, y=lat, fill=flyrs_change)) + 
	
	scale_fill_distiller(type="seq", palette="Greens", direction=1, name="Flowering frequency,\n1900-2023", breaks=c(0,0.25,0.5)) + labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = SppExt[1:2], ylim = SppExt[3:4], expand = FALSE) +
	
	theme_minimal(base_size=9) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.1, "in"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=8), legend.title=element_text(size=9), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())


{cairo_pdf(paste("output/figures/flfrq_change_", taxon, ".pdf", sep=""), width=5, height=5)

flrfrq_change

}
dev.off()





