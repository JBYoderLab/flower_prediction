# Inventorying iNat observations for candidate taxa
# Assumes local environment 
# jby 2023.07.18

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flowering_prediction")
# setwd("~/Documents/Academic/Active_projects/flowering_prediction")

library("tidyverse")

library("sf")
library("ggspatial")

library("patchwork")

source("scripts/R/get_inat.R") # attaches rinat and hacks the key function
source("../shared/Rscripts/base_graphics.R")

#-------------------------------------------------------------------------
# Pull down iNat observations of target taxon with specific phenology code
# Joshua tree, e.g.: taxon_id=47785
# actually chosen for example figure ...
# Prunus ilicifolia = 57250
# Quercus lobata = 49011
# Lupinus bicolor = 50614
# Impatiens capensis = 47888
# Trillium erectum = 50855
# Liriodendron tulipifera = 53582
# Cypripedium reginae = 51434
# Carnegia gigantea = 54449
# Lythrum salicaria = 61321

taxon <- 61321 

# per https://forum.inaturalist.org/t/how-to-use-inaturalists-search-urls-wiki-part-2-of-2/18792
# term id: 12 for Plant Phenology then term_id_value: 13 =Flowering, 14 =Fruiting, 15 =Flower Budding, 21 =No evidence of flowering
# place id: place_id=53170 for Mojave Desert; generally don't add a place restriction to start, but expect to do some geo-fencing later


# trial run, to make sure it works as expected
test <- get_inat_obs(quality="research", taxon_id=taxon, term_id=12, term_value_id=14, year=2021, maxresults=1e4)

glimpse(test)


# ACTUALLY RUN THE THING
inat_pheno_data <- data.frame(matrix(0,0,7))
names(inat_pheno_data) <- c("scientific_name", "latitude", "longitude", "url", "image_url", "observed_on", "phenology")


# set parameters as variables
years <- 2008:2022 # to run everything

# to read back in and continue
# inat_pheno_data <- read.csv(paste("data/inat_phenology_data_", taxon, ".csv", sep=""), h=TRUE)

# okay let's pull this stuff down already
# n.b. for-looping this borks up in a way that makes me suspect it's overloading the API
for(y in years){

# y <- 2016

bud.y <- try(get_inat_obs(quality="research", taxon_id=taxon, term_id=12, term_value_id=15, year=y, maxresults=1e4))
Sys.sleep(5) # throttling under the API limit, maybe?
flo.y <- try(get_inat_obs(quality="research", taxon_id=taxon, term_id=12, term_value_id=13, year=y, maxresults=1e4))
Sys.sleep(5)  
fru.y <- try(get_inat_obs(quality="research", taxon_id=taxon, term_id=12, term_value_id=14, year=y, maxresults=1e4))
Sys.sleep(5)  
non.y <- try(get_inat_obs(quality="research", taxon_id=taxon, term_id=12, term_value_id=21, year=y, maxresults=1e4))
Sys.sleep(5)  


if(class(bud.y)=="data.frame") bud.o <- bud.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="Flower Budding", year=gsub("(\\d{4})-.+","\\1", observed_on)) else bud.o <- NULL

if(class(flo.y)=="data.frame") flo.o <- flo.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="Flowering", year=gsub("(\\d{4})-.+","\\1", observed_on)) else flo.o <- NULL

if(class(fru.y)=="data.frame") fru.o <- fru.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="Fruiting", year=gsub("(\\d{4})-.+","\\1", observed_on)) else fru.o <- NULL

if(class(non.y)=="data.frame") non.o <- non.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="No Evidence of Flowering", year=gsub("(\\d{4})-.+","\\1", observed_on)) else non.o <- NULL


inat_pheno_data <- rbind(inat_pheno_data, bud.o, flo.o, fru.o, non.o)


if(!file.exists("data")) dir.create("data") # make sure there's a folder to write to!

write.table(inat_pheno_data, paste("data/inat_phenology_data_", taxon, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

} # END loop over years

# expect error messages if searches return zero obs with a given phenology status; this may not be a problem, but see what the final data table looks like
glimpse(inat_pheno_data) # okay

inat_pheno_data$id <- gsub(".+/(\\d+)", "\\1", inat_pheno_data$url)

write.table(inat_pheno_data, paste("data/inat_phenology_data_", taxon, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

table(inat_pheno_data$year, inat_pheno_data$phenology) # and okay

# ALL data, even without phenology annotation ............................
inat_all_data <- data.frame(matrix(0,0,7))
names(inat_all_data) <- c("scientific_name", "latitude", "longitude", "url", "image_url", "observed_on", "phenology")

for(y in years){

# y <- 2016

all.y <- try(get_inat_obs(quality="research", taxon_id=taxon, year=y, maxresults=1e4))
Sys.sleep(5) # throttling under the API limit, maybe?

if(class(all.y)=="data.frame") all.o <- all.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="Flower Budding", year=gsub("(\\d{4})-.+","\\1", observed_on)) else all.o <- NULL

inat_all_data <- rbind(inat_all_data, all.o)

write.table(inat_all_data, paste("data/inat_all_data_", taxon, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

} # END loop over years

glimpse(inat_all_data) # okay

inat_all_data$id <- gsub(".+/(\\d+)", "\\1", inat_all_data$url)

write.table(inat_all_data, paste("data/inat_all_data_", taxon, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

table(inat_all_data$year)


# put annotated and non-annotated records together ...
inat_data <- rbind(
	data.frame(inat_pheno_data[,c("scientific_name", "id", "longitude", "latitude", "year")], pheno=TRUE),
	data.frame(inat_all_data[!(inat_all_data$id%in%inat_pheno_data$id),c("scientific_name", "id", "longitude", "latitude", "year")], pheno=FALSE)
	) |> filter(!is.na(latitude))
glimpse(inat_data)
table(inat_data$pheno)

hist(inat_data$latitude)
hist(inat_data$longitude)

inat_data <- filter(inat_data, latitude>=20) # remove geographic outliers
inat_data <- filter(inat_data, longitude<=-40) # remove geographic outliers

write.table(inat_data, paste("data/inat_proposal_data_", taxon, ".csv", sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)


#-------------------------------------------------------------------------
# visualize for panels of a larger figure
inat_data <- read.csv(paste("data/inat_proposal_data_", taxon, ".csv", sep=""), h=TRUE) |> filter(year>=2008)

glimpse(inat_data)

# barplot of record count by phenology and year -----------
phenocols <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c')

flr.raw.ln <- table(inat_data$year, inat_data$pheno) %>% as.data.frame() %>% rename(year=Var1, phenology=Var2, observations=Freq) 

sampleN <- ggplot(flr.raw.ln, aes(x=year, y=observations, fill=phenology)) + 
geom_bar(stat="identity", position="stack") + 
scale_fill_manual(values=c("gray80", phenocols[4]), labels=c("Not annotated", "Phenology annotated")) + labs(x="Year", y="Records") + 
scale_x_discrete(labels=c("","",2010,"","","","",2015,"","","","",2020,"","")) + 
theme_minimal(base_size=24) + theme(legend.position=c(0.285, 0.9), legend.title=element_blank(), legend.box.background=element_rect(fill="white", color=NA), panel.grid.major.x=element_blank(), axis.title.x=element_blank(), legend.key.size=unit(0.75, "cm"))

sampleN

# map making ----------------------------------------------

# load some useful stuff ...................
states <- read_sf(dsn = "../data/spatial/10m_cultural/", lay="ne_10m_admin_1_states_provinces")
coast <- read_sf("../data/spatial/10m_physical/ne_10m_coastline", "ne_10m_coastline")
lakes <- read_sf("../data/spatial/10m_physical/ne_10m_lakes", "ne_10m_lakes")

mean(range(inat_data$latitude))
mean(range(inat_data$longitude))

# convert to a pretty polyconic projection
#zoom_to <- c(-112, 31)  # ~ center of the species range, west coastish
#zoom_to <- c(-80, 42)  # ~ center of the species range, east of the Mississippi
zoom_to <- c(-89, 43)  # ~ center of the species range, NAmish
zoom_level <- 2.5
# Lambert azimuthal equal-area projection around center of interest
target_crs <- sprintf('+proj=laea +lon_0=%f +lat_0=%f', zoom_to[1], zoom_to[2])
C <- 40075016.686   # ~ circumference of Earth in meters
x_span <- C / 2^zoom_level
y_span <- C / 2^(zoom_level+1)
zoom_to_xy <- st_transform(st_sfc(st_point(zoom_to), crs = 4326), crs = target_crs)
zoom_to_xy
disp_window <- st_sfc(
    st_point(st_coordinates(zoom_to_xy - c(x_span / 2, y_span / 2))),
    st_point(st_coordinates(zoom_to_xy + c(x_span / 2, y_span / 2))),
    crs = target_crs
)

# Records, illustrated ...
exmap <- ggplot() + 

geom_sf(data=coast, color="white", linewidth=1) + 
geom_sf(data=states, fill="gray80", color=NA) + 
geom_sf(data=states, fill=NA, color="white") + # turn off for NAm focus

geom_sf(data=lakes, fill="slategray3", color=NA) + 

geom_sf(data=st_as_sf(inat_data, coords=c("longitude", "latitude"), crs=4267), aes(shape=pheno, color=pheno), alpha=0.7, size=2) + 

#annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +

coord_sf(xlim = c(-3e6, 2.6e6), ylim = c(-2.5e6, 2.5e6), crs = target_crs, datum = target_crs) +
# nb don't change this if you don't want to f*ck up the aspect ratio
# originally 12 x 14; so to ... 
# Desert SW: xlim = c(-1e6, 8e5), ylim = c(-9e5, 7e5)
# Eastern Nam: xlim = c(-1.7e6, 2.1e6), ylim = c(-1.9e6, 1.5e6)
# most of NAm: xlim = c(-3e6, 2.6e6), ylim = c(-2.5e6, 2.5e6)

# legend
scale_color_manual(labels=c("Not annotated", "Phenology annotated"), values=c("gray50", phenocols[4])) + 
scale_shape_manual(labels=c("Not annotated", "Phenology annotated"), values=c(21,20)) + 

theme_minimal() + theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = "slategray3"), axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(), legend.position="none", legend.title=element_blank(), legend.box.background=element_rect(fill="white"))

# end figure-building ...

exmap


{png(paste("output/figures/example_records_", taxon, ".png", sep=""), height=720, width=480)

ggdraw() + draw_plot(exmap, 0, 0, 1, 0.6) + draw_plot(sampleN, 0, 0.6, 1, 0.4)

}
dev.off()

#-------------------------------------------------------------------------
# summaries for the proposal text

# Carnegia gigantea = 54449
# Liriodendron tulipifera = 53582
# Cypripedium reginae = 51434
# Lythrum salicaria = 61321

CAGI <- read.csv("data/inat_proposal_data_54449.csv", h=TRUE)
LITU <- read.csv("data/inat_proposal_data_53582.csv", h=TRUE)
CYRE <- read.csv("data/inat_proposal_data_51434.csv", h=TRUE)
LYSA <- read.csv("data/inat_proposal_data_61321.csv", h=TRUE)

glimpse(CAGI)
table(CAGI$pheno)[2]/nrow(CAGI)

glimpse(LITU)
table(LITU$pheno)[2]/nrow(LITU)

glimpse(CYRE)
table(CYRE$pheno)[2]/nrow(CYRE)

glimpse(LYSA)
table(LYSA$pheno)[2]/nrow(LYSA)

