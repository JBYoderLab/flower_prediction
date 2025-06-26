# working with PRISM historical monthlys
# Assumes local environment 
# jby 2024.10.03

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flower_prediction")

library("tidyverse")

library("terra")

library("prism")

prism_set_dl_dir("../data/PRISM") # system-specific; this goes to a directory shared among projects

# set parameters as variables
taxon <- 53405 # toyon!
# taxon <- 57250 # Prunus ilicifolia


#-------------------------------------------------------------------------
# Download the PRISM monthlies, if you haven't already (otherwise, skip)

for(yr in 1895:2023){ # note that years are up to you, 1895 is earliest avialable

# yr <- 2024

get_prism_monthlys(type="tmax", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="tmin", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="ppt", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmax", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmin", mon=1:12, year=yr, keepZip=FALSE)

}

#-------------------------------------------------------------------------
# process PRISM data layers into quarterly values for analysis

if(!dir.exists("../data/PRISM/quarterlies")) dir.create("../data/PRISM/quarterlies")

# parse monthly values into quarterlies
for(yr in 1895:2023){

# yr <- 2024

# FOR LOOP over quarters to do the thing ...
	for(q in 1:4){
	
	# q <- 1
	
	mos <- list(1:3,4:6,7:9,10:12)[[q]]
	
	# max temp in each quarter
	tmaxQ <- max(raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[3]))))
	
	writeRaster(tmaxQ, paste("../data/PRISM/quarterlies/tmax_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# min temp in each quarter
	tminQ <- min(raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[3]))))

	writeRaster(tminQ, paste("../data/PRISM/quarterlies/tmin_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# sum of precip in each quarter
	pptQ <- sum(raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[3]))))

	writeRaster(pptQ, paste("../data/PRISM/quarterlies/ppt_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	vpdmaxQ <- sum(raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[3]))))

	writeRaster(vpdmaxQ, paste("../data/PRISM/quarterlies/vpdmax_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	vpdminQ <- sum(raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[3]))))

	writeRaster(vpdminQ, paste("../data/PRISM/quarterlies/vpdmin_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	
	} # END loop over quarters

cat("Done with data from", yr, "\n\n")

} # END loop over years



# and now I have quarterly historical data cropped for the species-specific extent for all downstream analysis




