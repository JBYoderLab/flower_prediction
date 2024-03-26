# Scraping phenology-annotated iNat observations
# Assumes local environment 
# jby 2023.11.02

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flowering_prediction")
# setwd("~/Documents/Academic/Active_projects/flowering_prediction")

library("tidyverse")

source("scripts/R/get_inat.R") # attaches rinat and hacks the key function

#-------------------------------------------------------------------------
# Pull down iNat observations of target taxon with specific phenology code
# Joshua tree, e.g.: taxon_id=47785
# Toyon: 53405
# Carnegia gigantea: 54449
# Prunus ilicifolia = 57250
# Ericameria nauseosa = 57934

taxon <- 57934 

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
years <- 2010:2022 # to run everything

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
glimpse(inat_pheno_data) # ooh not bad
table(inat_pheno_data$year, inat_pheno_data$phenology) # uh wow!


#-------------------------------------------------------------------------
# visualize, if you like
inat_pheno_data <- read.csv(paste("data/inat_phenology_data_", taxon, ".csv", sep=""), h=TRUE)

flr.raw.ln <- table(inat_pheno_data$year, inat_pheno_data$phenology) %>% as.data.frame() %>% rename(year=Var1, phenology=Var2, observations=Freq)

if(!file.exists("output")) dir.create("output") # make sure there's a folder to write to!
if(!file.exists("output/figures")) dir.create("output/figures") # make sure there's a folder to write to!


{cairo_pdf(paste("output/figures/iNat_obs_raw_", taxon, ".pdf", sep=""), width=11, height=5)

ggplot(flr.raw.ln, aes(x=year, y=observations, fill=phenology)) + geom_bar(stat="identity", position="dodge") + labs(x="Year of observation", y="iNat records (research grade)") + theme_bw()

}
dev.off()

