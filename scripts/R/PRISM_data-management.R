# working with PRISM historical monthlys
# Assumes local environment 
# jby 2024.03.26

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flower_prediction")

library("tidyverse")
library("lubridate")

library("raster")

library("prism")

prism_set_dl_dir("../data/PRISM") # system-specific; this goes to a directory shared among projects

# set parameters as variables
taxon <- 53405 # toyon!
# taxon <- 57250 # Prunus ilicifolia

# Crop extent --- CHANGE THIS TO MATCH YOUR TARGET TAXON
SppExt <- extent(-125, -109, 23, 42) # for toyon; need to adjust accordingly


#-------------------------------------------------------------------------
# Download the PRISM monthlies, if you haven't already (otherwise, skip)

for(yr in 1895:2023){ # note that years are up to you, 1895 is earliest avialable

# yr <- 2021

get_prism_monthlys(type="tmax", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="tmin", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="ppt", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmax", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmin", mon=1:12, year=yr, keepZip=FALSE)

}


#-------------------------------------------------------------------------
# process PRISM data layers into cropped quarterly values for analysis

# make a place to stash files --- note taxon specificity, because of extent crop
if(!dir.exists(paste("data/PRISM/annual.", taxon, sep=""))) dir.create(paste("data/PRISM/annual.", taxon, sep=""))

# parse monthly values into quarterlies
for(yr in 1895:2023){

# yr <- 2021

# FOR LOOP over quarters to do the thing ...
	for(q in 1:4){
	
	# q <- 1
	
	mos <- list(1:3,4:6,7:9,10:12)[[q]]
	
	# max temp in each quarter
	tmaxQ <- crop(max(raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[3])))), SppExt)
	
	writeRaster(tmaxQ, paste("data/PRISM/annual.", taxon, "/tmax_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# min temp in each quarter
	tminQ <- crop(min(raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[3])))), SppExt)

	writeRaster(tminQ, paste("data/PRISM/annual.", taxon, "/tmin_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# sum of precip in each quarter
	pptQ <- crop(sum(raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[3])))), SppExt)

	writeRaster(pptQ, paste("data/PRISM/annual.", taxon, "/ppt_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	vpdmaxQ <- crop(sum(raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[3])))), SppExt)

	writeRaster(vpdmaxQ, paste("data/PRISM/annual.", taxon, "/vpdmax_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	vpdminQ <- crop(sum(raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[3])))), SppExt)

	writeRaster(vpdminQ, paste("data/PRISM/annual.", taxon, "/vpdmin_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	
	} # END loop over quarters

cat("Done with data from", yr, "\n\n")

} # END loop over years


# and now I have quarterly historical data cropped for the species-specific extent for all downstream analysis




