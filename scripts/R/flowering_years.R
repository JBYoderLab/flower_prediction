# Analyzing predicted historical flowering 
# jby 2024.03.27

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flower_prediction")
# setwd("~/Documents/Academic/Active_projects/flower_prediction")

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

sdm.pres <- read_sf("../data/Yucca/Jotr_SDM2023_range/Jotr_SDM2023_range.shp") # this is for Joshua tree

# flowering observation data
obs <- read.csv(paste("output/flowering_obs_climate_", taxon, ".csv", sep=""))
# raster files of predicted prFL
pred.files <- list.files(paste("output/BART/predictions_", taxon, sep=""), pattern=".bil", full=TRUE)

# useful bits and bobs
MojExt <- extent(-119, -112, 33, 40) # Mojave extent, maybe useful

jotr.preds <- c("pptY1Y2", "pptY0Y1", "vpdmaxY0", "vpdminY0Y1", "tminY0", "tmaxY0Y1")

sdm.buff <- st_buffer(st_transform(sdm.pres[,2], crs=3857), 1000) # put a 1km buffer on the range polygons
yubr.buff <- st_buffer(st_transform(yubr.pres, crs=3857), 1000)
yuja.buff <- st_buffer(st_transform(yuja.pres, crs=3857), 1000)

#-------------------------------------------------------------------------
# Process prediction layers

# without RI --------------------------------
jotr.histStack <- raster::stack(sapply(jotr.files, function(x) crop(raster::raster(x), MojExt)))
names(jotr.histStack) <- paste("prFL",1900:2023,sep=".")
projection(jotr.histStack)<-CRS("+init=epsg:4269")

jotr.histStack
writeRaster(jotr.histStack, "output/BART/jotr_BART_predicted_flowering_1900-2023_nomask.grd", overwrite=TRUE)

jotr.maskHist <- mask(jotr.histStack, st_transform(sdm.buff, crs=4269), touches=TRUE)

writeRaster(jotr.maskHist, "output/BART/jotr_BART_predicted_flowering_1900-2023.grd", overwrite=TRUE)

# build a data frame
hist.flowering <- cbind(coordinates(jotr.maskHist), as.data.frame(jotr.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(jotr.hist.flowering)

# with RI ------------------------------------
jotr.ri.histStack <- raster::stack(sapply(jotr.ri.files, function(x) crop(raster::raster(x), MojExt)))
names(jotr.ri.histStack) <- paste("prFL",1900:2023,sep=".")
projection(jotr.ri.histStack)<-CRS("+init=epsg:4269")

jotr.ri.histStack
writeRaster(jotr.ri.histStack, "output/BART/jotr_BART_RI_predicted_flowering_1900-2023_nomask.grd", overwrite=TRUE)

jotr.ri.maskHist <- mask(jotr.ri.histStack, st_transform(sdm.buff, crs=4269), touches=TRUE)

writeRaster(jotr.ri.maskHist, "output/BART/jotr_BART_RI_predicted_flowering_1900-2023.grd", overwrite=TRUE)

# build a data frame
ri.hist.flowering <- cbind(coordinates(jotr.ri.maskHist), as.data.frame(jotr.ri.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="RI.prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(jotr.ri.hist.flowering)


# build a data frame with base and RI versions
all.hist.flowering <- jotr.hist.flowering |> left_join(jotr.ri.hist.flowering)

glimpse(all.hist.flowering)

write.table(all.hist.flowering, "output/historic_flowering_reconst.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE) # write out and read in

all.hist.flowering <- read.csv("output/historic_flowering_reconst.csv") %>% filter(year !=2023)


#-------------------------------------------------------------------------
# Flowering years

# best-power prediction thresholds from the original models
# jotr: 0.26
# jotr with RI: 0.25

# YUBR: 0.30
# YUJA: 0.23

flyrs <- all.hist.flowering %>% group_by(lon,lat) %>% 
summarize(
	flyrs_all=length(which(prFL>=0.26)), 
	flyrs_1900_1929=length(which(prFL>=0.26 & year<=1929)), 
	flyrs_1990_2019=length(which(prFL>=0.26 & year>=1990 & year<=2019)),
	flyrs_RI_all=length(which(RI.prFL>=0.25)), 
	flyrs_RI_1900_1929=length(which(RI.prFL>=0.25 & year<=1929)), 
	flyrs_RI_1990_2019=length(which(RI.prFL>=0.25 & year>=1990 & year<=2019)),
	) |> mutate(
	flyrs_change=flyrs_1990_2019-flyrs_1900_1929,
	flyrs_RI_change=flyrs_RI_1990_2019-flyrs_RI_1900_1929
	) # rangewide model
range(all.hist.flowering$year) # remember: 1900-2022
glimpse(flyrs)

write.table(flyrs, "output/jotr_reconstructed_flowering_years.csv", sep=",", col.names=TRUE, row.names=FALSE)
# flyrs <- read.csv("output/jotr_reconstructed_flowering_years.csv")

# Inspect the flowering years .. does this make biological sense?
# without RI
quantile(flyrs$flyrs_all, c(0.025,0.5,0.975))
quantile(flyrs$flyrs_all, c(0.025,0.5,0.975))/123 # median 0.24

quantile(flyrs$flyrs_1900_1929, c(0.025,0.5,0.975))/30 # median 0.20
quantile(flyrs$flyrs_1990_2019, c(0.025,0.5,0.975))/30 # median 0.27

quantile(flyrs$flyrs_change, c(0.025,0.5,0.975)) # median 2; -3 to +6
mean(flyrs$flyrs_change) # mean 1.63

# with RI
quantile(flyrs$flyrs_RI_all, c(0.025,0.5,0.975))
quantile(flyrs$flyrs_RI_all, c(0.025,0.5,0.975))/123 # median 0.16

quantile(flyrs$flyrs_RI_1900_1929, c(0.025,0.5,0.975))/30 # median 0.17
quantile(flyrs$flyrs_RI_1990_2019, c(0.025,0.5,0.975))/30 #  median 0.13

quantile(flyrs$flyrs_RI_change, c(0.025,0.5,0.975)) # median 0; -6 to +8
mean(flyrs$flyrs_RI_change) # mean 0.22

# histograms
ggplot(flyrs, aes(x=flyrs_1900_1929)) + geom_histogram(bins=30)
ggplot(flyrs, aes(x=flyrs_1990_2019)) + geom_histogram(bins=30)

# distributions of flowering years
{cairo_pdf("output/figures/flowering-years_jotr.pdf", width=3.5, height=2.5)

ggplot(flyrs.jotr, aes(x=flyrs_all)) + geom_histogram(fill="#ccece6") + labs(x="Projected flowering years, 1900-2022", y="Grid cells") + 

geom_vline(xintercept=median(flyrs.jotr$flyrs_all), color="#006d2c") +

theme_minimal() + theme(legend.position="none", plot.margin=margin(0.1,0.15,0.1,0.15, "inches"))

}
dev.off()


# maps
# map elements
sdm.pres <- read_sf("../data/Yucca/Jotr_SDM2023_range/Jotr_SDM2023_range.shp")

states <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_admin_1_states_provinces")
coast <- read_sf("../data/spatial/10m_physical/ne_10m_coastline", "ne_10m_coastline")


flrfrq_map <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=flyrs.jotr, aes(x=lon, y=lat, fill=flyrs_all/123)) + 
	
	scale_fill_distiller(type="seq", palette="Greens", direction=1, name="Flowering frequency,\n1900-2022", breaks=c(0,0.25,0.5)) + labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-119.5, -112), ylim = c(33.5, 38.3), expand = FALSE) +
	
	theme_minimal(base_size=9) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.1, "in"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=8), legend.title=element_text(size=9), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())


{cairo_pdf("output/figures/flfrq_map_jotr.pdf", width=5, height=5)

flrfrq_map

}
dev.off()

flrfrq_change <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
		
	geom_tile(data=flyrs.jotr, aes(x=lon, y=lat, fill=flyrs_change)) + 
		
	scale_fill_gradient2(low="#762a83", mid="white", high="#1b7837", name="Change in flowering years,\n1990-2019 vs 1900-1929", breaks=seq(-12.5,12.5,by=2.5), labels=c("",-10,"",-5,"",0,"",5,"",10,"")) + labs(x="Longitude", y="Latitude") + 
			
	coord_sf(xlim = c(-119.5, -112), ylim = c(33.5, 38.3), expand = FALSE) +
	
	theme_minimal(base_size=9) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.1, "in"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=8), legend.title=element_text(size=9), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())


{cairo_pdf("output/figures/base_change_map_jotr.pdf", width=5, height=5)

flrfrq_change

} 
dev.off()

