# Analyzing predicted historical flowering 
# jby 2025.08.07

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
library("ggspatial")

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
# Trend in flowering frequency, 1900-2024

# correlation coefficients across all years
timecor <- historic.flowering %>% group_by(lat, lon, usfs) %>% do(broom::tidy(cor.test(~prFL+year, data=., method="spearman"))) %>% ungroup() %>% rename(rho=estimate)

timecor$elev_m <- terra::extract(elev, timecor[,c("lon", "lat")])$SR_50M 

sumFlr <- historic.flowering %>% group_by(lat, lon, usfs) %>% summarize(mnPrFlr = mean(prFL), sdPrFlr = sd(prFL), CVPrFlr = sdPrFlr/(mnPrFlr)) %>% ungroup() %>% left_join(timecor)

glimpse(sumFlr)

# change in early vs recent 30-year periods ... 
earlyFlr <- historic.flowering %>% filter(year %in% 1901:1930) %>% group_by(lat, lon, usfs) %>% summarize(eMnPrFlr = mean(prFL), eSdPrFlr = sd(prFL), eCVPrFlr = eSdPrFlr/(eMnPrFlr))
glimpse(earlyFlr)

recentFlr <- historic.flowering %>% filter(year %in% 1991:2020) %>% group_by(lat, lon, usfs) %>% summarize(rMnPrFlr = mean(prFL), rSdPrFlr = sd(prFL), rCVPrFlr = rSdPrFlr/(rMnPrFlr))
glimpse(recentFlr)

sumFlr <- sumFlr %>% left_join(earlyFlr) %>% left_join(recentFlr) %>% mutate(delMnPrFlr = rMnPrFlr-eMnPrFlr, delSdPrFlr = rSdPrFlr-eSdPrFlr, delCVPrFlr = rCVPrFlr - eCVPrFlr)

glimpse(sumFlr)

# 1901-1930 MnPrFlr .............................
quantile(sumFlr$eMnPrFlr, c(0.025, 0.5, 0.975)) # median 0.80, 95% CI 0.66, 0.89
quantile(sumFlr$eMnPrFlr[sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.79, 95% CI 0.67, 0.88
quantile(sumFlr$eMnPrFlr[!sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.81, 95% CI 0.66, 0.88
wilcox.test(eMnPrFlr~usfs, data=sumFlr, alt="g") # false > true, p < 2.2e-16

cor.test(~eMnPrFlr+lat, data=sumFlr, m="spearman")
# cor = 0.84, p < 2.2e-16 --- bigger PrFlr farther north
cor.test(~eMnPrFlr+lat, data=filter(sumFlr, usfs), m="spearman")
# cor = 0.81 p < 2.2e-16 --- bigger PrFlr farther north
cor.test(~eMnPrFlr+lat, data=filter(sumFlr, !usfs), m="spearman")
# cor = 0.83, p < 2.2e-16 --- bigger PrFlr farther north

cor.test(~eMnPrFlr+elev_m, data=sumFlr, m="spearman")
# cor = 0.02, n.s. --- bigger PrFlr farther north
cor.test(~eMnPrFlr+elev_m, data=filter(sumFlr, usfs), m="spearman")
# cor = 0.05, p = 5.7e-05 --- bigger PrFlr farther north
cor.test(~eMnPrFlr+elev_m, data=filter(sumFlr, !usfs), m="spearman")
# cor = -0.01, n.s. --- bigger PrFlr farther north


# 1991-2020 MnPrFlr .............................
quantile(sumFlr$rMnPrFlr, c(0.025, 0.5, 0.975)) # median 0.84, 95% CI 0.70, 0.90
quantile(sumFlr$rMnPrFlr[sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.82, 95% CI 0.70, 0.91
quantile(sumFlr$rMnPrFlr[!sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.84, 95% CI 0.69, 0.90
wilcox.test(rMnPrFlr~usfs, data=sumFlr, alt="g") # false > true, p < 2.2e-16

cor.test(~rMnPrFlr+lat, data=sumFlr, m="spearman")
# cor = 0.84, p < 2.2e-16 --- bigger PrFlr farther north
cor.test(~rMnPrFlr+lat, data=filter(sumFlr, usfs), m="spearman")
# cor = 0.82, p < 2.2e-16 --- bigger PrFlr farther north
cor.test(~rMnPrFlr+lat, data=filter(sumFlr, !usfs), m="spearman")
# cor = 0.84, p < 2.2e-16 --- bigger PrFlr farther north

cor.test(~rMnPrFlr+elev_m, data=sumFlr, m="spearman")
# cor = 0.03, p = 1.5e-05 --- bigger PrFlr farther north
cor.test(~rMnPrFlr+elev_m, data=filter(sumFlr, usfs), m="spearman")
# cor = 0.05, p = 1.5e-05 --- bigger PrFlr farther north
cor.test(~rMnPrFlr+elev_m, data=filter(sumFlr, !usfs), m="spearman")
# cor = 0.01, p = 0.22 --- n.s.

# change in MnPrFlr .............................
ggplot(sumFlr, aes(x=delMnPrFlr, fill=usfs)) + geom_histogram(position="dodge") + theme_bw() 
quantile(sumFlr$delMnPrFlr, c(0.025, 0.5, 0.975)) # median 0.02; 95% CI -0.02, 0.08
quantile(sumFlr$delMnPrFlr[sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.02; 95% CI -0.03, 0.08
quantile(sumFlr$delMnPrFlr[!sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.03; 95% CI -0.02, 0.07
wilcox.test(delMnPrFlr~usfs, data=sumFlr, alt="g") # false > true, p = 5.6e-09

cor.test(~delMnPrFlr+lat, data=sumFlr, m="spearman")
# cor = -0.22, p < 2.2e-16 --- bigger increase in PrFlr farther south?
cor.test(~delMnPrFlr+lat, data=filter(sumFlr, usfs), m="spearman")
# cor = -0.08, p = 4.9e-10 --- slightly bigger increase in PrFlr farther south?
cor.test(~delMnPrFlr+lat, data=filter(sumFlr, !usfs), m="spearman")
# cor = -0.34, p < 2.2e-16 --- bigger increase in PrFlr farther south?

cor.test(~delMnPrFlr+elev_m, data=sumFlr, m="spearman")
# cor = 0.04, p = 1.8e-08 --- greater increase in PrFlr farther uphill
cor.test(~delMnPrFlr+elev_m, data=filter(sumFlr, usfs), m="spearman")
# cor = 0.01, p = 0.26 --- n.s.
cor.test(~delMnPrFlr+elev_m, data=filter(sumFlr, !usfs), m="spearman")
# cor = 0.07, p =5.4e-12 --- greater increase in PrFlr farther uphill

# 1901-1930 CVPrFlr .............................
quantile(sumFlr$eCVPrFlr, c(0.025, 0.5, 0.975)) # median 0.06, 95% CI 0.02, 0.10
quantile(sumFlr$eCVPrFlr[sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.06, 95% CI 0.02, 0.10
quantile(sumFlr$eCVPrFlr[!sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.05, 95% CI 0.01, 0.10
wilcox.test(eCVPrFlr~usfs, data=sumFlr, alt="l") # true > false, p < 2.2e-16

cor.test(~eCVPrFlr+lat, data=sumFlr, m="spearman")
# cor = -0.61, p < 2.2e-16 --- bigger CVPrFlr farther south
cor.test(~eCVPrFlr+lat, data=filter(sumFlr, usfs), m="spearman")
# cor = -0.64, p < 2.2e-16 --- bigger CVPrFlr farther south
cor.test(~eCVPrFlr+lat, data=filter(sumFlr, !usfs), m="spearman")
# cor = -0.56, p < 2.2e-16 --- bigger CVPrFlr farther south

cor.test(~eCVPrFlr+elev_m, data=sumFlr, m="spearman")
# cor = -0.09, p < 2.2e-16 --- bigger CVPrFlr downhill
cor.test(~eCVPrFlr+elev_m, data=filter(sumFlr, usfs), m="spearman")
# cor = -0.09, p = 4.0e-13 --- bigger CVPrFlr downhill
cor.test(~eCVPrFlr+elev_m, data=filter(sumFlr, !usfs), m="spearman")
# cor = -0.10, p < 2.2e-16 --- (barely) bigger CVPrFlr downhill


# 1991-2020 CVPrFlr .............................
quantile(sumFlr$rCVPrFlr, c(0.025, 0.5, 0.975)) # median 0.06, 95% CI 0.03, 0.13
quantile(sumFlr$rCVPrFlr[sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.07, 95% CI 0.03, 0.14
quantile(sumFlr$rCVPrFlr[!sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.05, 95% CI 0.03, 0.11
wilcox.test(rCVPrFlr~usfs, data=sumFlr, alt="l") # true > false, p < 2.2e-16

cor.test(~rCVPrFlr+lat, data=sumFlr, m="spearman")
# cor = -0.63, p < 2.2e-16 --- bigger CVPrFlr farther south
cor.test(~rCVPrFlr+lat, data=filter(sumFlr, usfs), m="spearman")
# cor = -0.72, p < 2.2e-16 --- bigger CVPrFlr farther south
cor.test(~rCVPrFlr+lat, data=filter(sumFlr, !usfs), m="spearman")
# cor = -0.52, p < 2.2e-16 --- bigger CVPrFlr farther south

cor.test(~rCVPrFlr+elev_m, data=sumFlr, m="spearman")
# cor = -0.10, p < 2.2e-16 --- bigger CVPrFlr downhill
cor.test(~rCVPrFlr+elev_m, data=filter(sumFlr, usfs), m="spearman")
# cor = -0.09, p = 1.4e-12 --- bigger CVPrFlr downhill
cor.test(~rCVPrFlr+elev_m, data=filter(sumFlr, !usfs), m="spearman")
# cor = -0.12, p < 2.2e-16 --- bigger CVPrFlr downhill


# change in CVPrFlr .............................
ggplot(sumFlr, aes(x=delCVPrFlr, fill=usfs)) + geom_histogram(position="dodge") + theme_bw() 
quantile(sumFlr$delCVPrFlr, c(0.025, 0.5, 0.975)) # median 0.01; 95% CI -0.02, 0.05
quantile(sumFlr$delCVPrFlr[sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.01; 95% CI -0.01, 0.06
quantile(sumFlr$delCVPrFlr[!sumFlr$usfs], c(0.025, 0.5, 0.975)) # median 0.005; 95% CI -0.02, 0.04
wilcox.test(delCVPrFlr~usfs, data=sumFlr, alt="l") # false < true, p < 2.2e-16 HMM

cor.test(~delCVPrFlr+lat, data=sumFlr, m="spearman")
# cor = -0.10, p < 2.2e-16 --- much bigger increase farther south
cor.test(~delCVPrFlr+lat, data=filter(sumFlr, usfs), m="spearman")
# cor = -0.26, p < 2.2e-16 --- much bigger increase in CVFlr farther south
cor.test(~delCVPrFlr+lat, data=filter(sumFlr, !usfs), m="spearman")
# cor = 0.09, p < 2.2e-16 --- slightly bigger increase in CVFlr farther NORTH

cor.test(~delCVPrFlr+elev_m, data=sumFlr, m="spearman")
# cor = -0.04, p = 2.5e-07 --- greater increase in CVFlr at lower elevations
cor.test(~delCVPrFlr+elev_m, data=filter(sumFlr, usfs), m="spearman")
# cor = -0.03, p = 0.06 --- n.s.
cor.test(~delCVPrFlr+elev_m, data=filter(sumFlr, !usfs), m="spearman")
# cor = -0.04, p = 3.3e-05 --- greater increase in CVFlr at lower elevations




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
# Comparisons to static SDM

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
length(which(sumFlr$SDM_gain))/nrow(sumFlr) # 4% gain

length(which(sumFlr$SDM_gain & sumFlr$usfs)) # 367 cells
length(which(sumFlr$SDM_gain & sumFlr$usfs))/length(which(sumFlr$usfs)) # 4% gain


length(which(sumFlr$SDM_loss)) # 1063 cells
length(which(sumFlr$SDM_loss))/nrow(sumFlr) # 7% loss

length(which(sumFlr$SDM_loss & sumFlr$usfs)) # 550 cells
length(which(sumFlr$SDM_loss & sumFlr$usfs))/length(which(sumFlr$usfs)) # 9% loss


write.table(sumFlr, paste0("output/DART_and_SDM_summary_", taxon, ".csv"), sep=",", col.names=TRUE, row.names=FALSE)

# sumFlr <- read.csv(paste0("output/DART_and_SDM_summary_", taxon, ".csv"))

# And the stats!
# SDM suitability in early period vs mean flowering over early period
cor.test(~SDM_1901_1930+eMnPrFlr, data=sumFlr, m="spearman")
# rho = 0.04, p = 7.6e-07
cor.test(~SDM_1901_1930+eMnPrFlr, data=filter(sumFlr, usfs), m="spearman")
# rho = -0.19, p < 2.2e-16
cor.test(~SDM_1901_1930+eMnPrFlr, data=filter(sumFlr, !usfs), m="spearman")
# rho = 0.36, p < 2.2e-16


# SDM suitability in early period vs CV of flowering over early period
cor.test(~SDM_1901_1930+eCVPrFlr, data=sumFlr, m="spearman")
# rho = 0.11, p < 2.2e-16
cor.test(~SDM_1901_1930+eCVPrFlr, data=filter(sumFlr, usfs), m="spearman")
# rho = 0.23, p < 2.2e-16
cor.test(~SDM_1901_1930+eCVPrFlr, data=filter(sumFlr, !usfs), m="spearman")
# rho = -0.07, p = 9.0e-14


# SDM suitability in recent period vs mean flowering over recent period
cor.test(~SDM_1991_2020+rMnPrFlr, data=sumFlr, m="spearman")
# rho = -0.05, p = 4.4e-11
cor.test(~SDM_1991_2020+rMnPrFlr, data=filter(sumFlr, usfs), m="spearman")
# rho = -0.20, p < 2.2e-16
cor.test(~SDM_1991_2020+rMnPrFlr, data=filter(sumFlr, !usfs), m="spearman")
# rho = 0.19, p < 2.2e-16


# SDM suitability in recent period vs CV of flowering over recent period
cor.test(~SDM_1991_2020+rCVPrFlr, data=sumFlr, m="spearman")
# rho = 0.21, p < 2.2e-16
cor.test(~SDM_1991_2020+rCVPrFlr, data=filter(sumFlr, usfs), m="spearman")
# rho = 0.001, p = 0.89
cor.test(~SDM_1991_2020+rCVPrFlr, data=filter(sumFlr, !usfs), m="spearman")
# rho = 0.12, p < 2.2e-16


# change in SDM suitability vs change in mean flowering intensity
cor.test(~SDM_change+delMnPrFlr, data=sumFlr, m="spearman") 
# cor = 0.19, p < 2.2e-16 ... THERE we go
cor.test(~SDM_change+delMnPrFlr, data=filter(sumFlr, usfs), m="spearman") 
# cor = 0.06, p = 1.8e-06
cor.test(~SDM_change+delMnPrFlr, data=filter(sumFlr, !usfs), m="spearman") 
# cor = 0.28, p < 2.2e-16


# change in SDM suitability vs change in CV of flowering intensity
cor.test(~SDM_change+delCVPrFlr, data=sumFlr, m="spearman") 
# cor = -0.24, p < 2.2e-16 ... OHO
cor.test(~SDM_change+delCVPrFlr, data=filter(sumFlr, usfs), m="spearman") 
# cor = -0.21, p < 2.2e-16
cor.test(~SDM_change+delCVPrFlr, data=filter(sumFlr, !usfs), m="spearman") 
# cor = -0.30, p < 2.2e-16


### OH AND
cor.test(~SDM_change+elev_m, data=sumFlr, method="spearman") # n.s.
cor.test(~SDM_change+elev_m, data=filter(sumFlr, usfs), method="spearman") # n.s.

cor.test(~SDM_change+lat, data=sumFlr, method="spearman") # rho = -0.13, p < 2.2e-16 whaaa
cor.test(~SDM_change+lat, data=filter(sumFlr, usfs), method="spearman") # rho = -0.08, p < 1.2e-11 whaaa


t.test(sumFlr$delMnPrFlr[sumFlr$SDM_gain], sumFlr$delMnPrFlr[sumFlr$SDM_loss])
#  p-value = 0.13 welp
t.test(sumFlr$delCVPrFlr[sumFlr$SDM_gain], sumFlr$delCVPrFlr[sumFlr$SDM_loss])
#  p-value < 2.2e-16 OH indeed

# does SDM predicted gain/loss align with change in flowering intensity?
mod1 <- lm(delMnPrFlr~SDM_gl, data=sumFlr)
summary(mod1)
TukeyHSD(aov(mod1))
#                        diff          lwr         upr     p adj
# Loss-Gain      -0.001665956 -0.004682841 0.001350928 0.3983628
# No change-Gain  0.011116058  0.008672968 0.013559147 0.0000000
# No change-Loss  0.012782014  0.010875678 0.014688349 0.0000000

# does SDM predicted gain/loss align with change in CV of flowering intensity?
mod2 <- lm(delCVPrFlr~SDM_gl, data=sumFlr)
summary(mod2)
TukeyHSD(aov(mod2))
#                        diff          lwr         upr     p adj
# Loss-Gain       0.008334204  0.006262368  0.010406040 0.0e+00
# No change-Gain -0.003637466 -0.005315250 -0.001959682 1.2e-06
# No change-Loss -0.011971670 -0.013280840 -0.010662500 0.0e+00


ggplot(sumFlr, aes(x=delMnPrFlr)) + geom_histogram() + facet_wrap("SDM_gl", ncol=1, scale="free")

ggplot(filter(sumFlr, SDM_loss), aes(x=rho, fill=SDM_gain)) + geom_histogram() 
# not super-clear that's a divergent subset, hmm
t.test(sumFlr$rho[sumFlr$SDM_gain], sumFlr$rho[!sumFlr$SDM_gain]) # OKAY well

# map SDM results and trends in different regions ...

sdm_vs_delMean <- ggplot() +
	geom_point(data=filter(sumFlr, SDM_gl=="No change"), aes(x=SDM_change, y=delMnPrFlr, color=SDM_gl), alpha=0.2) +
	geom_point(data=filter(sumFlr, SDM_gl!="No change"), aes(x=SDM_change, y=delMnPrFlr, color=SDM_gl), alpha=0.5) +
	geom_smooth(data=sumFlr, aes(x=SDM_change, y=delMnPrFlr), method="lm", linewidth=0.5, color="black", se=FALSE) +
	scale_color_manual(values=c('#66c2a5','#fc8d62','#8da0cb'), name="Change in\nsuitability") + 
	labs(x = "SDM predicted change in Pr(present)", y = "Change in mean flowering intensity") + 
	theme_bw(base_size=10) + 
	theme(legend.position="none", axis.title=element_text(size=11))

sdm_vs_delMean

sdm_vs_delCV <- ggplot() +
	geom_point(data=filter(sumFlr, SDM_gl=="No change"), aes(x=SDM_change, y=delCVPrFlr, color=SDM_gl), alpha=0.2) +
	geom_point(data=filter(sumFlr, SDM_gl!="No change"), aes(x=SDM_change, y=delCVPrFlr, color=SDM_gl), alpha=0.5) +
	geom_smooth(data=sumFlr, aes(x=SDM_change, y=delCVPrFlr), method="lm", linewidth=0.5, color="black", se=FALSE) +
	scale_color_manual(values=c('#66c2a5','#fc8d62','#8da0cb'), name="Change in\nsuitability") + 
	labs(x = "SDM predicted change in Pr(present)", y = "Change in CV of flowering intensity") + 
	theme_bw(base_size=10) + 
	theme(legend.position="none", axis.title=element_text(size=11))

sdm_vs_delCV

sdm_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=sumFlr, aes(x=lon, y=lat, fill=SDM_gl)) + 

	geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
	geom_sf(data=usfs.buff, fill=NA, color="black", linewidth=0.3, linetype=2) + 
	
	annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
	annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
	
	scale_fill_manual(values=c('#66c2a5','#fc8d62','#8da0cb'), name="SDM-predicted change in\nhabitat suitability,\n1991-2020 vs 1901-1930") + 
	labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +
	annotation_scale(location = "bl", width_hint = 0.3) + 
	annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.15, "in"), pad_y = unit(0.25, "in"),
	style = north_arrow_fancy_orienteering, height=unit(0.75, "in"), width=unit(0.5, "in")) +
	
	theme_minimal(base_size=10) + theme(legend.position="bottom", legend.key.width=unit(0.1, "inches"),
	legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=8), legend.title=element_text(size=10, margin=margin(0, 10, 5, 0, unit="pt")), legend.title.position="top", panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())

sdm_map

sdm_trend_mean <- ggplot(sumFlr, aes(x=delMnPrFlr, fill=SDM_gl, color=usfs)) + geom_histogram(bins=20) + 
	facet_wrap("SDM_gl", ncol=1, scale="free_y") + 
	scale_fill_manual(values=c('#66c2a5','#fc8d62','#8da0cb'), name="Change in suitability,\n1991-2020 vs 1901-1930") + 
	scale_color_manual(values=c('white', 'black')) +
	geom_text(data=data.frame(SDM_gl="Gain", x = 0.08, y = 47, label="Cells in\ncore range", usfs=NA), aes(x=x, y=y, label=label), color="black", lineheight=0.75, size=3) +
	labs(x = "Change in mean flowering intensity", y = "Raster grid cells") + 
	theme_bw(base_size=12) + 
	theme(legend.position="none", axis.title=element_text(size=10))

sdm_trend_mean

sdm_trend_cv <- ggplot(sumFlr, aes(x=delCVPrFlr, fill=SDM_gl, color=usfs)) + geom_histogram(bins=20) + 
	facet_wrap("SDM_gl", ncol=1, scale="free_y") + 
	scale_fill_manual(values=c('#66c2a5','#fc8d62','#8da0cb'), name="Change in suitability,\n1991-2020 vs 1901-1930") + 
	scale_color_manual(values=c('white', 'black')) +
	geom_text(data=data.frame(SDM_gl="Gain", x = 0.05, y = 47, label="Cells in\ncore range", usfs=NA), aes(x=x, y=y, label=label), color="black", lineheight=0.75, size=3) +
	labs(x = "Change in CV of flowering intensity", y = "Raster grid cells") + 
	theme_bw(base_size=12) + 
	theme(legend.position="none", axis.title=element_text(size=10))

sdm_trend_cv


{cairo_pdf(paste("output/figures/SDM_vs_hindcast_", taxon, ".pdf", sep=""), width=9, height=3.5)

ggdraw() + draw_plot(sdm_map, 0, 0, 0.25, 1) + draw_plot(sdm_vs_delMean, 0.25, 0, 0.375, 1) + draw_plot(sdm_vs_delCV, 0.625, 0, 0.375, 1) + draw_plot_label(x=c(0, 0.25, 0.625), y=1, label=c("A", "B", "C"), size=18)

}
dev.off()


# map this ------------------------------------------------
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")

# map elements
states <- ne_states(country="united states of america", returnclass="sf")
countries <- ne_countries(scale=10, continent="north america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

MnFlr <- sumFlr %>% dplyr::select(lat, lon, eMnPrFlr, rMnPrFlr) %>% rename(`1901-1930`=eMnPrFlr, `1991-2020`=rMnPrFlr) %>% pivot_longer(3:4, names_to="period", values_to="MnPrFlr")

facet_labs <- data.frame(period=c("1901-1930", "1991-2020"), x=-116, y=41.75)

MnFlr_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=MnFlr, aes(x=lon, y=lat, fill=MnPrFlr)) + 

	geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
	geom_sf(data=usfs.buff, fill=NA, color="black", linewidth=0.3, linetype=2) + 

	geom_text(data=facet_labs, aes(label=period, x=x, y=y), size=6, hjust=1, vjust=1) +
	
	facet_wrap("period") +

	annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
	annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
	
	scale_fill_gradient(low="#e5f5f9", high="#2ca25f", name="Mean flowering intensity\n(hindcast flowering frequency)", breaks=c(0.7,0.8,0.9)) + 
	labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +

	theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.25, "inches"), legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12, margin=margin(0, 10, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank(), strip.text=element_blank(), strip.background=element_blank())
	
MnFlr_map


ChangeMnFlr_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=sumFlr, aes(x=lon, y=lat, fill=delMnPrFlr)) + 

	geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
	geom_sf(data=usfs.buff, fill=NA, color="black", linewidth=0.3, linetype=2) + 

	annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
	annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
	
	scale_fill_gradient2(low="#f1a340", mid="white", high="#998ec3", name="Change in mean\nflowering intensity", breaks=c(-0.05, 0, 0.05, 0.1)) + 
	labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +
	
	theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.25, "inches"), legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12, margin=margin(0, 10, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())
	
ChangeMnFlr_map

CVFlr <- sumFlr %>% dplyr::select(lat, lon, eCVPrFlr, rCVPrFlr) %>% rename(`1901-1930`=eCVPrFlr, `1991-2020`=rCVPrFlr) %>% pivot_longer(3:4, names_to="period", values_to="CVPrFlr")


CVFlr_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=CVFlr, aes(x=lon, y=lat, fill=CVPrFlr)) + 

	geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
	geom_sf(data=usfs.buff, fill=NA, color="black", linewidth=0.3, linetype=2) + 

	geom_text(data=facet_labs, aes(label=period, x=x, y=y), size=6, hjust=1, vjust=1) +
	
	facet_wrap("period") +

	annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
	annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
	
	scale_fill_gradient(low="#feedde", high="#d94801", name="CV of flowering intensity\n(hindcast flowering frequency)", breaks=c(0,0.05,0.1,0.15)) + 
	labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +

	theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.25, "inches"), legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12, margin=margin(0, 10, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank(), strip.text=element_blank(), strip.background=element_blank())
	
CVFlr_map


ChangeCVFlr_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=sumFlr, aes(x=lon, y=lat, fill=delCVPrFlr)) + 

	geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
	geom_sf(data=usfs.buff, fill=NA, color="black", linewidth=0.3, linetype=2) + 

	annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
	annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
	
	scale_fill_gradient2(low="#4d9221", mid="white", high="#c51b7d", name="Change in CV of\nflowering intensity", breaks=c(-0.05, 0, 0.05, 0.1)) + 
	labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +
	
	theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.25, "inches"), legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12, margin=margin(0, 10, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())
	
ChangeCVFlr_map

SDM <- sumFlr %>% dplyr::select(lat, lon, SDM_1901_1930, SDM_1991_2020) %>% rename(`1901-1930`=SDM_1901_1930, `1991-2020`=SDM_1991_2020) %>% pivot_longer(3:4, names_to="period", values_to="SDMPrPresent")

SDMs_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=SDM, aes(x=lon, y=lat, fill=SDMPrPresent)) + 

	geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
	geom_sf(data=usfs.buff, fill=NA, color="black", linewidth=0.3, linetype=2) + 

	geom_text(data=facet_labs, aes(label=period, x=x, y=y), size=6, hjust=1, vjust=1) +
	
	facet_wrap("period") +

	annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
	annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
	
	scale_fill_gradient(low="#f7fcfd", high="#810f7c", name="SDM-predicted\nPr(present)", breaks=c(0,0.25,0.5,0.75,1)) + 
	labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +

	theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.25, "inches"), legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12, margin=margin(0, 10, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank(), strip.text=element_blank(), strip.background=element_blank())
	
SDMs_map


ChangeSDM_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=sumFlr, aes(x=lon, y=lat, fill=SDM_change)) + 

	geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
	geom_sf(data=usfs.buff, fill=NA, color="black", linewidth=0.3, linetype=2) + 

	annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
	annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
	
	scale_fill_gradient2(low="#bf812d", mid="white", high="#35978f", name="Change in SDM-\npredicted Pr(present)", breaks=c(-0.5, 0, 0.5)) + 
	labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +
	
	theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.25, "inches"), legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12, margin=margin(0, 10, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())
	
ChangeSDM_map



{png(paste("output/figures/DART_flowering_summaries_map_v3_", taxon, ".png", sep=""), width=750, height=1000)

ggdraw() + draw_plot(MnFlr_map, 0, 0.66, 0.66, 0.33) + draw_plot(ChangeMnFlr_map, 0.66, 0.66, 0.33, 0.33)  + draw_plot(CVFlr_map, 0, 0.33, 0.66, 0.33) + draw_plot(ChangeCVFlr_map, 0.66, 0.33, 0.33, 0.33) + draw_plot(SDMs_map, 0, 0, 0.66, 0.33) + draw_plot(ChangeSDM_map, 0.66, 0, 0.33, 0.33) + draw_plot_label(label=LETTERS[1:9], x=rep(c(0.02, 0.33, 0.675), 3), y=rep(c(0.99, 0.66, 0.33), each=3), size=20)

}
dev.off()



#-------------------------------------------------------------------------
# lat/elev correlations … over time?

# new correlation setup needed ...
historic.flowering$elev_m <- terra::extract(elev, historic.flowering[,c("lon", "lat")])$SR_50M 

elev_cor <- historic.flowering %>% group_by(year) %>% do(broom::tidy(cor.test(~prFL+elev_m, data=., method="spearman"))) %>% ungroup() %>% rename(rho_elev=estimate)

elev_cor_usfs <- historic.flowering %>% filter(usfs) %>% group_by(year) %>% do(broom::tidy(cor.test(~prFL+elev_m, data=., method="spearman"))) %>% ungroup() %>% rename(rho_elev_core=estimate)

lat_cor <- historic.flowering %>% group_by(year) %>% do(broom::tidy(cor.test(~prFL+lat, data=., method="spearman"))) %>% ungroup() %>% rename(rho_lat=estimate)

lat_cor_usfs <- historic.flowering %>% filter(usfs) %>% group_by(year) %>% do(broom::tidy(cor.test(~prFL+lat, data=., method="spearman"))) %>% ungroup() %>% rename(rho_lat_core=estimate)

# merge it all ....
space_cors <- elev_cor %>% dplyr::select(year, rho_elev) %>% 
	left_join(dplyr::select(elev_cor_usfs, year, rho_elev_core)) %>% 
	left_join(dplyr::select(lat_cor, year, rho_lat)) %>% 
	left_join(dplyr::select(lat_cor_usfs, year, rho_lat_core))

glimpse(space_cors)

ggplot(space_cors, aes(x = year, y=rho_elev)) + geom_smooth() + geom_point(alpha=0.1)
ggplot(space_cors, aes(x = year, y=rho_elev_core)) + geom_smooth() + geom_point(alpha=0.1)
ggplot(space_cors, aes(x = year, y=rho_lat)) + geom_smooth() + geom_point(alpha=0.1)
ggplot(space_cors, aes(x = year, y=rho_lat_core)) + geom_smooth() + geom_point(alpha=0.1)

# HUH

cor.test(~rho_elev+year, data=space_cors, method="sp") # n.s.
cor.test(~rho_elev_core+year, data=space_cors, method="sp") # n.s.
cor.test(~rho_lat+year, data=space_cors, method="sp") # p = 0.04
cor.test(~rho_lat_core+year, data=space_cors, method="sp") # n.s.

# welp?


#-------------------------------------------------------------------------
# flowering activity changes vs predictor changes ...

# change in predictors .................................
sample_sites <- sumFlr %>% dplyr::select(lat, lon)

# port over from model training ...
preds <- c("tmin.y1q3", "ppt.y0q1", "ppt.y1q4", "tmin.y1q4", "ppt.y1q3", "vpdmin.y0q1")

sample_predictor_history <- data.frame(matrix(0,0,3+length(preds)))
names(sample_predictor_history) <- c("lat", "lon", "year", preds)

for(yr in 1900:2024){ # loop over years

# yr = 1900
y0 = yr
y1 = yr - 1

YrPredFiles <- c(list.files("../data/PRISM/quarterlies/", pattern=paste0("_",y0,"Q\\d\\.bil"), full.names = TRUE), list.files("../data/PRISM/quarterlies/", pattern=paste0("_",y1,"Q\\d\\.bil"), full.names = TRUE))

YrPredStack <- rast(YrPredFiles)
names(YrPredStack) <- paste0(rep(c("ppt", "tmax", "tmin", "vpdmax", "vpdmin"), each=4), ".", rep(c("y0","y1"), each=20), rep(paste0("q",1:4), 10))

sample_predictor_history <- rbind(sample_predictor_history, data.frame(sample_sites, year=yr, terra::extract(YrPredStack[[preds]], sample_sites[,c("lon","lat")])[,preds]))

}

# now compile trends ...
pred_history_long <- sample_predictor_history %>% pivot_longer(all_of(preds), names_to="predictor", values_to="value") %>% mutate(predictor = factor(predictor, preds))

glimpse(pred_history_long)

levels(pred_history_long$predictor) <- c("Tmin Y1Q3", "PPT Y0Q1", "PPT Y1Q4", "Tmin Y1Q4", "PPT Y1Q3", "VPDmin Y0Q1")

# change between the early-recent timeframe ...
predictor_mean_early <- pred_history_long %>% filter(year%in%1901:1930) %>% group_by(lat, lon, predictor) %>% summarize(mean_1901_1930=mean(value))

predictor_mean_late <- pred_history_long %>% filter(year%in%1991:2020) %>% group_by(lat, lon, predictor) %>% summarize(mean_1991_2020=mean(value))

predictor_mean_change <- full_join(predictor_mean_early, predictor_mean_late) %>% mutate(mean_change=mean_1991_2020-mean_1901_1930) %>% left_join(sumFlr)

glimpse(predictor_mean_change)

eMnPrFlr, rMnPrFlr, eCVPrFlr, rCVPrFlr


pred_changes <- predictor_mean_change %>% dplyr::select(lat, lon, predictor, mean_1901_1930, mean_1991_2020) %>% pivot_longer(starts_with("mean"), names_to="timeframe", values_to="pred_value")
pred_changes$timeframe[pred_changes$timeframe=="mean_1901_1930"] <- "1901_1930"
pred_changes$timeframe[pred_changes$timeframe=="mean_1991_2020"] <- "1991_2020"

MnPrFlr_changes <- predictor_mean_change %>% dplyr::select(lat, lon, predictor, eMnPrFlr, rMnPrFlr) %>% pivot_longer(ends_with("PrFlr"), names_to="timeframe", values_to="MnPrFlr")
MnPrFlr_changes$timeframe[MnPrFlr_changes$timeframe=="eMnPrFlr"] <- "1901_1930"
MnPrFlr_changes$timeframe[MnPrFlr_changes$timeframe=="rMnPrFlr"] <- "1991_2020"

CVPrFlr_changes <- predictor_mean_change %>% dplyr::select(lat, lon, predictor, eCVPrFlr, rCVPrFlr) %>% pivot_longer(ends_with("PrFlr"), names_to="timeframe", values_to="CVPrFlr")
CVPrFlr_changes$timeframe[CVPrFlr_changes$timeframe=="eCVPrFlr"] <- "1901_1930"
CVPrFlr_changes$timeframe[CVPrFlr_changes$timeframe=="rCVPrFlr"] <- "1991_2020"

pred_flr_changes <- pred_changes %>% left_join(MnPrFlr_changes) %>% left_join(CVPrFlr_changes)
glimpse(pred_flr_changes)

plot_slice <- pred_flr_changes %>% group_by(timeframe, predictor) %>% slice_sample(n=300)

flr_change_summ <- pred_flr_changes %>% group_by(timeframe, predictor) %>% summarize(med_pred = median(pred_value), med_MnPrFlr = median(MnPrFlr)) %>% pivot_wider(names_from=timeframe, values_from=c(med_pred, med_MnPrFlr))


MnPredChange <- ggplot(predictor_mean_change, aes(x=mean_change, y=delMnPrFlr, color=predclass)) + 
	geom_point(alpha=0.1) +
	geom_smooth(method="lm", color="black", se=FALSE, linewidth=0.25) + 
	facet_wrap("predictor", nrow=2, scale="free_x") +
	scale_color_manual(values = c("#a6cee3", "#fdbf6f", "#cab2d6")) + 
	labs(x="Change in predictor, 1991-2020 vs 1901-1930", y="Change in flowering intensity") +
	theme_bw(base_size=18) +
	theme(legend.position="none")

MnPredChange

CVPredChange <- ggplot(predictor_mean_change, aes(x=mean_change, y=delCVPrFlr, color=predclass)) + 
	geom_point(alpha=0.1) +
	geom_smooth(method="lm", color="black", se=FALSE, linewidth=0.25) + 
	facet_wrap("predictor", nrow=2, scale="free_x") +
	scale_color_manual(values = c("#a6cee3", "#fdbf6f", "#cab2d6")) + 
	labs(x="Change in predictor, 1991-2020 vs 1901-1930", y="Change in CV of flowering intensity") +
	theme_bw(base_size=18) +
	theme(legend.position="none")



{png(paste0("output/figures/pred_flr_change_", taxon, ".png"), width=750, height=1000)

ggdraw() + draw_plot(MnPredChange, 0, 0.5, 1, 0.5) + draw_plot(CVPredChange, 0, 0, 1, 0.5) + draw_plot_label(label=c("A", "B"), x=0, y=c(1, 0.5), size=24)

}
dev.off()

