# working with PRISM historical monthlys
# Assumes local environment
# jby 2026.08.03

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flower_modeling")

library("tidyverse")

library("terra")

library("prism")

prism_set_dl_dir("../data/PRISM") # system-specific; this goes to a directory shared among projects

#-------------------------------------------------------------------------
# Download the PRISM monthlies, if you haven't already (otherwise, skip)

for(yr in 1895:2025){ # note that years are up to you, 1895 is earliest available

# yr <- 2024

get_prism_monthlys(type="tmax", mon=1:12, resolution = "4km", year=yr, keepZip=FALSE)
get_prism_monthlys(type="tmin", mon=1:12, resolution = "4km", year=yr, keepZip=FALSE)
get_prism_monthlys(type="ppt", mon=1:12, resolution = "4km", year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmax", mon=1:12, resolution = "4km", year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmin", mon=1:12, resolution = "4km", year=yr, keepZip=FALSE)

}

#-------------------------------------------------------------------------
# process PRISM data layers into quarterly values for analysis

if(!dir.exists("../data/PRISM/quarterlies")) dir.create("../data/PRISM/quarterlies")

# parse monthly values into quarterlies
for(yr in 1895:2025){

# yr <- 2024

# FOR LOOP over quarters to do the thing ...
	for(q in 1:4){

	# q <- 1

	mos <- list(1:3,4:6,7:9,10:12)[[q]]

	# max temp in each quarter
	tmaxQ <- max(rast(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[1], resolution="4km"))), rast(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[2], resolution="4km"))), rast(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[3], resolution="4km"))))

	writeRaster(tmaxQ, paste0("../data/PRISM/quarterlies/tmax_", yr, "Q", q, ".tiff"), overwrite=TRUE)

	# min temp in each quarter
	tminQ <- min(rast(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[1], resolution = "4km"))), rast(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[2], resolution = "4km"))), rast(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[3], resolution = "4km"))))

	writeRaster(tminQ, paste0("../data/PRISM/quarterlies/tmin_", yr, "Q", q, ".tiff"), overwrite=TRUE)

	# sum of precip in each quarter
	pptQ <- sum(rast(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[1], resolution = "4km"))), rast(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[2], resolution = "4km"))), rast(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[3], resolution = "4km"))))

	writeRaster(pptQ, paste0("../data/PRISM/quarterlies/ppt_", yr, "Q", q, ".tiff"), overwrite=TRUE)

	vpdmaxQ <- sum(rast(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[1], resolution = "4km"))), rast(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[2], resolution = "4km"))), rast(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[3], resolution = "4km"))))

	writeRaster(vpdmaxQ, paste0("../data/PRISM/quarterlies/vpdmax_", yr, "Q", q, ".tiff"), overwrite=TRUE)

	vpdminQ <- sum(rast(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[1], resolution = "4km"))), rast(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[2], resolution = "4km"))), rast(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[3], resolution = "4km"))))

	writeRaster(vpdminQ, paste0("../data/PRISM/quarterlies/vpdmin_", yr, "Q", q, ".tiff"), overwrite=TRUE)


	} # END loop over quarters

cat("Done with data from", yr, "\n\n")

} # END loop over years



# and now I have quarterly historical data for all downstream analysis




