# Analyzing predicted historical flowering 
# jby 2025.07.08

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

# Read in prediction layers
historic.stack.masked <- rast(paste0("output/models/DART_predicted_flowering_", taxon, "_1900-2024_masked.tiff"))

historic.stack.masked.usfs <- rast(paste0("output/models/DART_predicted_flowering_", taxon, "_1900-2024_masked_to_range.tiff"))



#-------------------------------------------------------------------------
# USA-NPN data

npn <- read.csv("data/datasheet_1752018515618/status_intensity_observation_data.csv") %>% mutate(Year = year(ymd(npn$Observation_Date)), Month=month(ymd(npn$Observation_Date)))

glimpse(npn)

# what do we have?
table(npn$Site_ID, npn$Year) 
# okay a LOT of sites, actually; decent temporal cover also
table(npn$Phenophase_Description) 

npn %>% filter(Site_ID==8957, Year==2019) %>% group_by(Phenophase_Description, Month) %>% summarize(Np = n())

npn_by_year_loc <- npn %>% group_by(Latitude, Longitude, Year) %>% summarize(Nobs = n())
npn_by_year_loc

npn %>% filter(Phenophase_Status==1) %>% group_by(Phenophase_Description, Phenophase_Status) %>% summarize(Np = n())
# oookay, this is workable: Phenophase_Status is -1 if absent, 1 if present. So I can get a proxy for intensity by dividing n(1) for a given year/location by (n(-1)+n(1)), I think?


# map npn sites ------------------------------------------------
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")

# map elements
states <- ne_states(country="united states of america", returnclass="sf")
countries <- ne_countries(scale=10, continent="north america", returnclass="sf")
coast <- ne_coastline(scale=10, returnclass="sf")

	
ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=3) + 
	geom_sf(data=countries, fill="cornsilk3", color="antiquewhite4") + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 

	geom_sf(data=states, fill=NA, color="antiquewhite4") + 
	
	geom_sf(data=usfs.buff, fill="darkseagreen", color=NA) + 

	geom_point(data=npn_by_year_loc, aes(x=Longitude, y=Latitude, color=Year)) + 
	
	annotate("text", x=-119, y=36, label="CA", size=12, color="white", alpha=0.35) + 
	annotate("text", x=-117.5, y=40, label="NV", size=12, color="white", alpha=0.35) + 
	
#	scale_fill_discrete(low="#e5f5f9", high="#2ca25f", name="Mean Pr(flowers),\n1900-2024", breaks=c(0.7,0.8,0.9)) + 
	labs(x="Longitude", y="Latitude") + 
		
	coord_sf(xlim = c(-125.5,-115.5), ylim = c(32,42), expand = FALSE) +
	annotation_scale(location = "bl", width_hint = 0.3) + 
	annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.15, "in"), pad_y = unit(0.25, "in"), style = north_arrow_fancy_orienteering, height=unit(0.75, "in"), width=unit(0.5, "in")) +
	
	theme_minimal(base_size=12) + theme(legend.position="bottom", legend.key.width=unit(0.25, "inches"), legend.key.height=unit(0.1, "in"), legend.direction="horizontal", axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.05,0.1,0.05,0.01), "inches"), legend.box.spacing=unit(0.001,"inches"), legend.box="horizontal", legend.text=element_text(size=10), legend.title=element_text(size=12, margin=margin(0, 10, 0, 10, unit="pt")), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())


#-------------------------------------------------------------------------
# aggregate into flowering intensity per year/site

# create a "flowering year" variable that wraps early-year observations of fruit into the previous year
npn$flr_yr <- npn$Year
npn$flr_yr[npn$Phenophase_Description%in%c("Fruits", "Recent fruit or seed drop", "Ripe fruits") & npn$Month<6] <- npn$Year[npn$Phenophase_Description%in%c("Fruits", "Recent fruit or seed drop", "Ripe fruits") & npn$Month<6]-1

glimpse(npn)

npn_intensity <- npn %>% filter(Phenophase_Status==1) %>% group_by(Site_ID, Latitude, Longitude, flr_yr) %>% 
				summarize(Ntot=n(), Nflr=length(which(Phenophase_Description%in%c("Fruits", "Recent fruit or seed drop", "Ripe fruits")))) %>% mutate(Prop_flr = Nflr/Ntot)
npn_intensity # I think this does it?


#-------------------------------------------------------------------------
# pair with model predictions

historic.stack.masked # read in at the top

npn_intensity_predictions <- npn_intensity %>% filter(flr_yr!=2025) %>% mutate(prFL=NA)

for(yr in unique(npn_intensity_predictions$flr_yr)){

# yr <- 2020

subs <- npn_intensity_predictions %>% filter(flr_yr==yr)

subs_pred <- terra::extract(historic.stack.masked[paste0("prFL.", yr)], subs[,c("Longitude", "Latitude")])

npn_intensity_predictions$prFL[npn_intensity_predictions$flr_yr==yr] <- subs_pred[,2]

} # this will throw an error for 2025 because we don't have predictions that year!

npn_intensity_predictions <- npn_intensity_predictions %>% filter(!is.na(prFL))

glimpse(npn_intensity_predictions) # 380 year-location records

length(table(npn_intensity_predictions$Site_ID)) # 197 unique locations (!)
range(table(npn_intensity_predictions$Site_ID))
median(table(npn_intensity_predictions$Site_ID), 0.5)

table(npn_intensity_predictions$flr_yr)

# moment of truth
cor.test(~Prop_flr+prFL, data=npn_intensity_predictions, method="spearman")
# rho = 0.23, p = 7.165e-06 


{cairo_pdf(paste("output/figures/DART_predictions_validation_", taxon, ".pdf", sep=""), width=4.5, height=4)

ggplot(npn_intensity_predictions, aes(x=Prop_flr, y=prFL)) + 
	geom_smooth(method="lm", color="white") + geom_point(alpha=0.5) + 
	labs(x="USA-NPN flowering intensity", y="Modeled Pr(flowering)") + 
	theme_bw(base_size=14)

}
dev.off()



# hey kid, wanna see something freaky
# (the correlation varies between years, I think, because within years there's only geographic variation, and
# the USA-NPN sites don't cover that wide a range of geographic variation ...)
ggplot(npn_intensity_predictions, aes(x=Prop_flr, y=prFL)) + 
	geom_smooth(method="lm") + geom_point() + 
	facet_wrap("flr_yr") +
	labs(x="USA-NPN flowering intensity", y="Modeled Pr(flowering)") + 
	theme_bw(base_size=14)





