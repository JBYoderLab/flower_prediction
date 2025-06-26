# working with phenology-annotated iNat observations
# Assumes local environment 
# jby 2025.05.21

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flower_prediction")

library("tidyverse")
library("terra")
library("sf")

# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia: 57250
# E nauseosa: 57394

#-------------------------------------------------------------------------
# read in iNat observations compiled using inat_phenology_obs.R

inat <- read.csv(paste("data/inat_phenology_data_", taxon, "_cleaned.csv", sep=""), h=TRUE) %>% mutate(observed_on = ymd(observed_on))

glimpse(inat) # how many raw observations?
table(inat$phenology) # by phenophase
table(inat$phenology, inat$year) # by phenophase and year

# Species area crop extent (deliberately generous)
# note the padding factors assume we're NW of 0,0
SppExt <- round(c(range(inat$longitude), range(inat$latitude)) * c(1.1,0.9,0.9,1.1),0) # for toyon; need to adjust accordingly

# is fruit ever observed in Jan, Feb, or March? That's really (maybe?) the PREVIOUS flowering year
# Take (again) the toyon example: per iNat, peak flowering is in ~June, peak fruiting in Nov
# To account for this we need to see how often fruit is logged BEFORE peak flowering
filter(inat, phenology=="Fruiting", month(observed_on)<6) # if this is not zero, we need to deal with these!

# create a "flowering year" variable that wraps early-year observations of fruit into the previous year
inat$flr_yr <- inat$year
inat$flr_yr[inat$phenology=="Fruiting" & month(inat$observed_on)<6] <- inat$year[inat$phenology=="Fruiting" & month(inat$observed_on)<6]-1

glimpse(inat)


#-------------------------------------------------------------------------
# organize iNat observations for extraction of summarized PRISM data

# data structure setup
flowering <- data.frame(matrix(0,0,5))
names(flowering) <- c("lon","lat","year", "prop_flr", "n_obs")

prism_temp_rast <- rast(paste("../data/PRISM/quarterlies/tmax_cropped_2010Q1.bil", sep="")) # raster grid base

# then LOOP over years in the raw data ....
for(yr in unique(inat$flr_yr)){

# yr <- 2020

obs <- rasterize(as.matrix(dplyr::filter(inat, flr_yr==yr)[,c("longitude","latitude")]), prism_temp_rast, fun=length, background=NA)

flr <- rasterize(as.matrix(dplyr::filter(inat, flr_yr==yr, phenology!="No Evidence of Flowering")[,c("longitude","latitude")]), prism_temp_rast, fun=length, background=0)

flrfrq <- flr/obs

outp <- data.frame(lon=crds(flrfrq)[,"x"], lat=crds(flrfrq)[,"y"], year=yr, prop_flr=as.data.frame(flrfrq)$values, n_obs=as.data.frame(obs)$values)

# put it all together
flowering <- rbind(flowering, outp)

cat("Done with year", yr, "\n")

} # END LOOP over years

head(flowering)
glimpse(flowering) # okay okay okay!
table(flowering$year) # smiling serenely
hist(flowering$prop_flr)
hist(flowering$n_obs)


write.table(flowering, paste("output/flowering_freq_rasterized_", taxon,".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE)
# flowering <- read.csv(paste("output/flowering_freq_rasterized_", taxon,".csv", sep=""))

#-------------------------------------------------------------------------
# attach PRISM data to flowering/not flowering observations

# we've got variables aggregated quarterly, assume we want Year 0 (year of obs) and Year -1 (year before obs)
varnames <- paste(rep(c("ppt", "tmax", "tmin", "vpdmax", "vpdmin"), each=4), paste(rep(c("y0", "y1"), each=20), paste("q", 1:4, sep=""), sep=""), sep=".")

prism_temp_rast <- rast(paste("../data/PRISM/quarterlies/tmax_cropped_2010Q1.bil", sep="")) # raster grid base

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
prismfiles <- paste("../data/PRISM/quarterlies/", rep(c("ppt", "tmax", "tmin", "vpdmax", "vpdmin"), each=4), "_cropped_", rep(c(y0, y1), each=20), "Q", 1:4, ".bil", sep="")

# assemble weather data predictors for a given year
preds <- rast(lapply(prismfiles, function(x) resample(rast(x), prism_temp_rast)))
names(preds) <- varnames

# pull subset of flowering observations for year
flsub <- subset(flowering, year==yr)
flsub <- cbind(flsub, terra::extract(preds, flsub[,c("lon","lat")], df=FALSE))

flr.clim <- rbind(flr.clim,flsub) 

write.table(flr.clim, paste0("output/flowering_freq_climate_", taxon, ".csv"), sep=",", col.names=TRUE, row.names=FALSE)

cat("Done with year", yr, "\n")


} # END LOOP over years

glimpse(flr.clim)


# and that's generated a data file we can feed into embarcadero ... in the next script!

#-------------------------------------------------------------------------
# map gridded records

flr.clim <- read.csv(paste0("output/flowering_freq_climate_", taxon, ".csv"))
glimpse(flr.clim)

flr.clim.summed <- flr.clim %>% group_by(lat, lon) %>% summarize(tot_obs = sum(n_obs))
glimpse(flr.clim.summed)

library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library("sf")

# map elements
states <- ne_states(country="united states of america", returnclass="sf")
countries <- ne_countries(scale=10, continent="north america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

# get the USFS range polygon for toyon
usfs.Range <- read_sf("../data/spatial/wpetry-USTreeAtlas-4999258/shp/photarbu/", layer="photarbu", crs=4326)
usfs.buff <- st_transform(st_buffer(st_transform(usfs.Range, crs=3857), 10000), crs=4326) %>% st_intersection(filter(ne_countries(scale=10, continent="north america", returnclass="sf"), name_en=="United States of America")) # 10km buffer?


{cairo_pdf(paste("output/figures/record_distribution_map_", taxon, ".pdf", sep=""), width=3.2, height=4.5)

ggplot() + 

geom_sf(data=coast, color="slategray2", linewidth=3) + 
geom_sf(data=countries, fill="antiquewhite4", color="antiquewhite4") + 
geom_sf(data=states, fill="darkseagreen3", color="antiquewhite4") + 
#geom_sf(data=filter(states, name=="California"), fill="cornsilk3", color="antiquewhite4") + 

geom_sf(data=usfs.buff, fill="darkseagreen4", color=NA, linewidth=0.3, linetype=2) + 

geom_tile(data=flr.clim.summed, aes(x=lon, y=lat, fill=log10(tot_obs))) + 

geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
annotate("text", x=-117.9, y=35.1, label="Core range", size=4, fontface="bold", color="darkseagreen4") + 

	
scale_fill_gradient(low="#fde0dd", high="#49006a", name=expression(log[10]("iNat records")), breaks=c(0,1,2)) + 
labs(x="Longitude", y="Latitude") + 
		
coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +
annotation_scale(location = "bl", width_hint = 0.3) + 
annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.15, "in"), pad_y = unit(0.25, "in"), style = north_arrow_fancy_orienteering, height=unit(0.75, "in"), width=unit(0.5, "in")) +
	
theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.25, "inches"), legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.01,0.1,0.01,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12, margin=margin(0, 10, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())

}
dev.off()




