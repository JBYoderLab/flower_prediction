# Analyzing predicted historical flowering 
# jby 2025.05.19

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flower_prediction")

# Clear the environment and load key packages
rm(list=ls())

library("tidyverse")

library("terra")
library("geodata")
library("sp")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

library("cowplot")

library("embarcadero")


#-------------------------------------------------------------------------
# Load and prep data

# set parameters as variables
taxon <- 53405 # toyon
# Prunus ilicifolia = 57250

# flowering observation data
obs <- read.csv(paste("output/flowering_freq_climate_", taxon, ".csv", sep=""))
# raster files of predicted prFL
pred.files <- list.files(paste("output/models/DART_predictions.", taxon, sep=""), pattern=".tiff", full=TRUE)

# useful bits and bobs
# Species area crop extent (deliberately generous)
# note the padding factors assume we're NW of 0,0
SppExt <- round(c(range(obs$lon), range(obs$lat)) * c(1.01,0.99,0.9,1.1),0) # for toyon; need to adjust accordingly

# get the USFS range polygon for toyon
usfs.Range <- read_sf("../data/spatial/wpetry-USTreeAtlas-4999258/shp/photarbu/", layer="photarbu", crs=4326)
usfs.buff <- st_transform(st_buffer(st_transform(usfs.Range, crs=3857), 10000), crs=4326) %>% st_intersection(filter(ne_countries(scale=10, continent="north america", returnclass="sf"), name_en=="United States of America")) # 10km buffer?

plot(usfs.Range)
plot(usfs.buff)

# new composite for "broad" study area
broad.Range <- read_sf(paste0("output/broad_range_polygon_", taxon, ".shp"))
broad.Range

plot(broad.Range)

#-------------------------------------------------------------------------
# Process prediction layers

historic.preds <- sapply(pred.files, simplify=TRUE, function(x) crop(rast(x), SppExt)) # read in prediction layers, crop to the species extent to reduce data volume later
historic.stack <- rast(historic.preds)
names(historic.stack) <- paste("prFL", 1900:2024, sep=".") # name the layers by year

historic.stack # what's this look like

writeRaster(historic.stack, paste("output/models/DART_predicted_flowering_", taxon, "_1900-2024_nomask.tiff", sep=""), overwrite=TRUE) # write out the cropped stack as a single raster file

# mask to the SDM polygon
historic.stack.masked <- mask(historic.stack, st_transform(broad.Range, crs=4326), touches=TRUE)
historic.stack.masked.usfs <- mask(historic.stack, st_transform(usfs.buff, crs=4326), touches=TRUE)

writeRaster(historic.stack.masked, paste0("output/models/DART_predicted_flowering_", taxon, "_1900-2024_masked.tiff"), overwrite=TRUE) # write out the masked layers as a single raster file

writeRaster(historic.stack.masked.usfs, paste0("output/models/DART_predicted_flowering_", taxon, "_1900-2024_masked_to_range.tiff"), overwrite=TRUE) # write out the masked layers as a single raster file


# read back in
historic.stack.masked <- rast(paste0("output/models/DART_predicted_flowering_", taxon, "_1900-2024_masked.tiff"))

historic.stack.masked.usfs <- rast(paste0("output/models/DART_predicted_flowering_", taxon, "_1900-2024_masked_to_range.tiff"))

#-------------------------------------------------------------------------
# pair with relevant info, reformat for analysis

elev <- rast("../data/spatial/SR_50M/SR_50M.tif")
elev

# build a data frame for analysis? 
historic.flowering.usfs <- cbind(crds(historic.stack.masked.usfs, df=TRUE), as.data.frame(historic.stack.masked.usfs)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

historic.flowering <- cbind(crds(historic.stack.masked, df=TRUE), as.data.frame(historic.stack.masked)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(historic.flowering.usfs) # oh that's a lot
glimpse(historic.flowering) # oh that's a LOT

historic.flowering$usfs <- paste(historic.flowering$lat, historic.flowering$lon) %in% paste(historic.flowering.usfs$lat, historic.flowering.usfs$lon)

table(historic.flowering$usfs)

write.table(historic.flowering, paste0("output/DART_predicted_flowering_broad_", taxon, "_1900-2024.csv"), sep=",")

# historic.flowering <- read.csv(paste0("output/DART_predicted_flowering_broad_", taxon, "_1900-2024.csv"))

#-------------------------------------------------------------------------
# Trend in flowering frequency, 1900-2023

# correlation coefficients
timecor <- historic.flowering %>% group_by(lat, lon, usfs) %>% do(broom::tidy(cor.test(~prFL+year, data=., method="spearman"))) %>% ungroup() %>% rename(rho=estimate)

timecor$elev_m <- terra::extract(elev, timecor[,c("lon", "lat")])$SR_50M 

sumFlr <- historic.flowering %>% group_by(lat, lon, usfs) %>% summarize(mnPrFlr = mean(prFL), sdPrFlr = sd(prFL), CVPrFlr = sdPrFlr/(mnPrFlr)) %>% ungroup() %>% left_join(timecor)

glimpse(sumFlr)


# Summary and analysis of average .........................
hist(sumFlr$mnPrFlr[sumFlr$usfs])
quantile(sumFlr$mnPrFlr[sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.81; 95% CI 0.69, 0.89
t.test(sumFlr$mnPrFlr[sumFlr$usfs]) # mean > 0, p < 2.2e-16

cor.test(~mnPrFlr+lat, data=filter(sumFlr, usfs), method="pearson")
# cor = 0.83, p = 2.2e-16 --- bigger average PrFlr farther north

cor.test(~mnPrFlr+elev_m, data=filter(sumFlr, usfs), method="pearson")
# cor = 0.07, p = 5.9e-08 --- slightly greater average PrFlr farther uphill

# Summary and analysis of variation (sd) ..................
hist(sumFlr$sdPrFlr[sumFlr$usfs])
quantile(sumFlr$sdPrFlr[sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.05, 95% CI 0.03, 0.08

cor.test(~sdPrFlr+lat, data=filter(sumFlr, usfs), method="pearson")
# cor = -0.72, p = 2.2e-16 --- more variation in the south

cor.test(~sdPrFlr+elev_m, data=filter(sumFlr, usfs), method="pearson")
# cor = -0.17, p < 2.2e-16 --- less variation in PrFlr farther uphill

# Summary and analysis of variation (CV) ..................
hist(sumFlr$CVPrFlr[sumFlr$usfs])
quantile(sumFlr$CVPrFlr[sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.07, 95% CI 0.03, 0.11

cor.test(~CVPrFlr+lat, data=filter(sumFlr, usfs), method="pearson")
# cor = -0.78, p < 2.2e-16 --- more variation in the south

cor.test(~CVPrFlr+elev_m, data=filter(sumFlr), method="pearson")
# cor = -0.10, p < 2.2e-16 --- less variation farther uphill



# Summary and analysis of the trend .......................
hist(sumFlr$rho[sumFlr$usfs]) # interesting
quantile(sumFlr$rho[sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.13; 95% CI -0.15, 0.46
t.test(sumFlr$rho) # mean > 0, p < 2.2e-16

cor.test(~rho+lat, data=filter(sumFlr, usfs), method="pearson")
# n.s.

cor.test(~rho+elev_m, data=filter(sumFlr, usfs), method="pearson")
# cor = 0.07, p = 9.0e-09 --- more positive trend at higher elevation


# some comparisons for paranoia, maybe? ...................
cor.test(~rho+mnPrFlr, data=filter(sumFlr, usfs), method="spearman")
# cor = -0.04, p = 0.001 ... hm yeah that's higher correlations at lower mean intensity?

# okay this is interestingly messy, not sure what I think of it
ggplot(filter(sumFlr, usfs), aes(x=mnPrFlr, y=rho)) + geom_point(alpha=0.1) + geom_smooth(method="lm", color="white", linewidth=0.5) + theme_bw()

ggplot(filter(sumFlr, usfs), aes(x=mnPrFlr, y=sdPrFlr)) + geom_point(alpha=0.1) + geom_smooth(method="lm", color="white", linewidth=0.5) + theme_bw()

# map this ------------------------------------------------
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")

# map elements
states <- ne_states(country="united states of america", returnclass="sf")
countries <- ne_countries(scale=10, continent="north america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

	
mnFlr_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=sumFlr, aes(x=lon, y=lat, fill=mnPrFlr)) + 

	geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
	geom_sf(data=usfs.buff, fill=NA, color="black", linewidth=0.3, linetype=2) + 

	annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
	annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
	
	scale_fill_gradient(low="#e5f5f9", high="#2ca25f", name="Mean Pr(flowers),\n1900-2024", breaks=c(0.7,0.8,0.9)) + 
	labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +
	annotation_scale(location = "bl", width_hint = 0.3) + 
	annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.15, "in"), pad_y = unit(0.25, "in"), style = north_arrow_fancy_orienteering, height=unit(0.75, "in"), width=unit(0.5, "in")) +
	
	theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.25, "inches"), legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12, margin=margin(0, 10, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())
	
cvFlr_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=sumFlr, aes(x=lon, y=lat, fill=CVPrFlr)) + 

	geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
	geom_sf(data=usfs.buff, fill=NA, color="black", linewidth=0.3, linetype=2) + 

	annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
	annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
	
	scale_fill_gradient(low="#fff7bc", high="#d95f0e", name="CV of Pr(flowers),\n1900-2024", breaks=c(0.05,0.1)) + labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +
	annotation_scale(location = "bl", width_hint = 0.3) + 
	annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.15, "in"), pad_y = unit(0.25, "in"), style = north_arrow_fancy_orienteering, height=unit(0.75, "in"), width=unit(0.5, "in")) +
	
	theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.25, "inches"),
	legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(),
	axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12, margin=margin(0, 10, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())

flrcor_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=sumFlr, aes(x=lon, y=lat, fill=rho)) + 

	geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
	geom_sf(data=usfs.buff, fill=NA, color="black", linewidth=0.3, linetype=2) + 
	
	annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
	annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
	
	scale_fill_gradient2(low="#f1a340", mid="white", high="#998ec3", name="Spearman's rho,\nPr(flowers) vs year", breaks=c(-0.25,0,0.25,0.5)) + labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +
	annotation_scale(location = "bl", width_hint = 0.3) + 
	annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.15, "in"), pad_y = unit(0.25, "in"), style = north_arrow_fancy_orienteering, height=unit(0.75, "in"), width=unit(0.5, "in")) +
	
	theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.25, "inches"), legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12, margin=margin(0, 20, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())

{cairo_pdf(paste("output/figures/DART_flowering_summaries_map_", taxon, ".pdf", sep=""), width=10, height=4.5)

ggdraw() + draw_plot(mnFlr_map, 0, 0, 0.33, 1) + draw_plot(cvFlr_map, 0.33, 0, 0.33, 1) + draw_plot(flrcor_map, 0.66, 0, 0.33, 1) + draw_plot_label(x=c(0.01, 0.34, 0.67), y=0.99, label=c("A", "B", "C"), size=20)

}
dev.off()

#-------------------------------------------------------------------------
# map training data

library("rnaturalearth")
library("rnaturalearthdata")

# map elements
states <- ne_states(country="united states of america", returnclass="sf")
countries <- ne_countries(scale=10, continent="north america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")


glimpse(obs) # confirm I've got that

propFlr_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 

	geom_sf(data=broad.Range, fill="darkseagreen", color=NA) + 
	
	geom_tile(data=filter(obs, year>=2021), aes(x=lon, y=lat, fill=prop_flr)) + 	
	
	facet_wrap("year", nrow=1) +
	
	scale_fill_gradient(low="#f7f7f7", high="#762a83", name="Observed prop. flowering", breaks=c(0.25,0.5,0.75)) + labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-124.6,-116), ylim = c(32,42), expand = FALSE) +
	
	theme_minimal(base_size=20) + theme(legend.position="bottom", legend.key.width=unit(0.35, "inches"),
	legend.key.height=unit(0.15, "in"), legend.direction="horizontal", axis.text=element_blank(),
	axis.title=element_blank(), plot.margin=unit(c(0.05,0.01,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=12), legend.title=element_text(size=16, margin=margin(0, 10, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank(), strip_text=element_text(size=22))

{cairo_pdf(paste("output/figures/FlrFrq_2021_2024_", taxon, ".pdf", sep=""), width=10, height=4.5)

propFlr_map

}
dev.off()

#-------------------------------------------------------------------------
# Comparison to static SDM

early_SDM <- rast(paste0("output/BART/SDM_", taxon, "_prediction_1901_1930.tiff"))
early_SDM_binary <- rast(paste0("output/BART/SDM_", taxon, "_prediction_binary_1901_1930.tiff"))
recent_SDM <- rast(paste0("output/BART/SDM_", taxon, "_prediction_1991_2020.tiff"))
recent_SDM_binary <- rast(paste0("output/BART/SDM_", taxon, "_prediction_binary_1991_2020.tiff"))

glimpse(sumFlr)

sumFlr$SDM_1901_1930 <- terra::extract(early_SDM, sumFlr[,c("lon", "lat")])$lyr.1
sumFlr$SDM_pres_1901_1930 <- terra::extract(early_SDM_binary, sumFlr[,c("lon", "lat")])$lyr.1

sumFlr$SDM_1991_2020 <- terra::extract(recent_SDM, sumFlr[,c("lon", "lat")])$lyr.1
sumFlr$SDM_pres_1991_2020 <- terra::extract(recent_SDM_binary, sumFlr[,c("lon", "lat")])$lyr.1

# comparisons
sumFlr$SDM_change <- sumFlr$SDM_1991_2020 - sumFlr$SDM_1901_1930 # more positive, greater gain
sumFlr$SDM_gain <- mapply(function(x,y) x==0 && y==1, sumFlr$SDM_pres_1901_1930, sumFlr$SDM_pres_1991_2020)
sumFlr$SDM_loss <- mapply(function(x,y) x==1 && y==0, sumFlr$SDM_pres_1901_1930, sumFlr$SDM_pres_1991_2020)

sumFlr$SDM_gl <- "No change"
sumFlr$SDM_gl[sumFlr$SDM_gain] <- "Gain"
sumFlr$SDM_gl[sumFlr$SDM_loss] <- "Loss"

glimpse(sumFlr)
length(which(sumFlr$SDM_gain)) # 629 cells
length(which(sumFlr$SDM_loss)) # 1063 cells

write.table(sumFlr, paste0("output/DART_and_SDM_summary_", taxon, ".csv"), sep=",", col.names=TRUE, row.names=FALSE)

# And the stats!
cor.test(~SDM_change+rho, data=sumFlr, method="spearman") 
# cor = 0.12, p < 2.2e-16 ... THERE we go
plot(recent_SDM - early_SDM) # hm

### OH AND
cor.test(~SDM_change+elev_m, data=sumFlr, method="spearman") # n.s.
cor.test(~SDM_change+elev_m, data=filter(sumFlr, usfs), method="spearman") # n.s.

cor.test(~SDM_change+lat, data=sumFlr, method="spearman") # rho = -0.13, p < 2.2e-16 whaaa
cor.test(~SDM_change+lat, data=filter(sumFlr, usfs), method="spearman") # rho = -0.08, p < 1.2e-11 whaaa




t.test(sumFlr$rho[sumFlr$SDM_gain], sumFlr$rho[sumFlr$SDM_loss & sumFlr$usfs], alt="less")
#  p-value = 0.52 welp

t.test(sumFlr$rho[sumFlr$SDM_gain]) # mean > 0, p < 2.2e-16
t.test(sumFlr$rho[sumFlr$SDM_loss]) # mean > 0, p < 2.2e-16 LMAO


ggplot(sumFlr, aes(x=rho)) + geom_histogram() + facet_wrap("SDM_gl", ncol=1, scale="free")

ggplot(filter(sumFlr, SDM_loss), aes(x=rho, fill=SDM_gain)) + geom_histogram() 
# not super-clear that's a divergent subset, hmm
t.test(sumFlr$rho[sumFlr$SDM_gain], sumFlr$rho[!sumFlr$SDM_loss]) # OKAY well

