# working with PRISM historical monthlys
# Assumes MAJEL environment 
# jby 2022.12.20

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/flowering_prediction")
# setwd("~/Documents/Academic/Active_projects/flowering_prediction")

library("tidyverse")
library("lubridate")

library("raster")

library("prism")

if(!dir.exists("data/PRISM")) dir.create("data/PRISM")

prism_set_dl_dir("data/PRISM")

# Crop extent --- CHANGE THIS TO MATCH YOUR TARGET TAXON
SppExt <- extent(-125, -109, 23, 42) # for toyon; need to adjust accordingly


#-------------------------------------------------------------------------
# process PRISM data layers into a timeline archive I can use later

# make places to stash actual values and anomaly values
if(!dir.exists("data/PRISM/annual")) dir.create("data/PRISM/annual")


# FOR LOOP each year in PRISM
# start with years of observation to play nice with the db
for(yr in 2008:2022){

# yr=2021

get_prism_monthlys(type="tmax", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="tmin", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="ppt", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmax", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmin", mon=1:12, year=yr, keepZip=FALSE)


# FOR LOOP over quarters to do the thing ...
	for(q in 1:4){
	
	# q <- 1
	
	mos <- list(1:3,4:6,7:9,10:12)[[q]]
	
	# max temp in each quarter
	tmaxQ <- crop(max(raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[3])))), SppExt)
	
	writeRaster(tmaxQ, paste("data/PRISM/annual/tmax_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# min temp in each quarter
	tminQ <- crop(min(raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[3])))), SppExt)

	writeRaster(tminQ, paste("data/PRISM/annual/tmin_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# sum of precip in each quarter
	pptQ <- crop(sum(raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[3])))), SppExt)

	writeRaster(pptQ, paste("data/PRISM/annual/ppt_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	vpdmaxQ <- crop(sum(raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[3])))), SppExt)

	writeRaster(vpdmaxQ, paste("data/PRISM/annual/vpdmax_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	vpdminQ <- crop(sum(raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[3])))), SppExt)

	writeRaster(vpdminQ, paste("data/PRISM/annual/vpdmin_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	
	} # END loop over quarters

# clean out raw data
sapply(prism_archive_subset("tmax", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))
sapply(prism_archive_subset("tmin", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))
sapply(prism_archive_subset("ppt", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))
sapply(prism_archive_subset("vpdmax", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))
sapply(prism_archive_subset("vpdmin", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))


} # END loop over years


#-------------------------------------------------------------------------
# normalize PRISM data to a 1981-2010 baseline

# make places to stash actual values and anomaly values
if(!dir.exists("data/PRISM/annual_normed_1981-2010")) dir.create("data/PRISM/annual_normed_1981-2010")
if(!dir.exists("data/PRISM/norms_1981-2010")) dir.create("data/PRISM/norms_1981-2010")

quarterlies <- list.files("data/PRISM/annual", pattern=".bil")

# first, our normals, 1981-2010
for(q in 1:4){

# q = 1

tmax_mn <- calc(brick(lapply(paste("data/PRISM/annual/tmax_cropped_",1981:2010,"Q",q,".bil", sep=""), raster)), mean)
tmax_sd <- calc(brick(lapply(paste("data/PRISM/annual/tmax_cropped_",1981:2010,"Q",q,".bil", sep=""), raster)), sd)

writeRaster(tmax_mn, paste("data/PRISM/norms_1981-2010/tmax_mean_cropped_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)
writeRaster(tmax_sd, paste("data/PRISM/norms_1981-2010/tmax_sd_cropped_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)


tmin_mn <- calc(brick(lapply(paste("data/PRISM/annual/tmin_cropped_",1981:2010,"Q",q,".bil", sep=""), raster)), mean)
tmin_sd <- calc(brick(lapply(paste("data/PRISM/annual/tmin_cropped_",1981:2010,"Q",q,".bil", sep=""), raster)), sd)

writeRaster(tmin_mn, paste("data/PRISM/norms_1981-2010/tmin_mean_cropped_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)
writeRaster(tmin_sd, paste("data/PRISM/norms_1981-2010/tmin_sd_cropped_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)


ppt_mn <- calc(brick(lapply(paste("data/PRISM/annual/ppt_cropped_",1981:2010,"Q",q,".bil", sep=""), raster)), mean)
ppt_sd <- calc(brick(lapply(paste("data/PRISM/annual/ppt_cropped_",1981:2010,"Q",q,".bil", sep=""), raster)), sd)

writeRaster(ppt_mn, paste("data/PRISM/norms_1981-2010/ppt_mean_cropped_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)
writeRaster(ppt_sd, paste("data/PRISM/norms_1981-2010/ppt_sd_cropped_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)


vpdmax_mn <- calc(brick(lapply(paste("data/PRISM/annual/vpdmax_cropped_",1981:2010,"Q",q,".bil", sep=""), raster)), mean)
vpdmax_sd <- calc(brick(lapply(paste("data/PRISM/annual/vpdmax_cropped_",1981:2010,"Q",q,".bil", sep=""), raster)), sd)

writeRaster(vpdmax_mn, paste("data/PRISM/norms_1981-2010/vpdmax_mean_cropped_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)
writeRaster(vpdmax_sd, paste("data/PRISM/norms_1981-2010/vpdmax_sd_cropped_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)


vpdmin_mn <- calc(brick(lapply(paste("data/PRISM/annual/vpdmin_cropped_",1981:2010,"Q",q,".bil", sep=""), raster)), mean)
vpdmin_sd <- calc(brick(lapply(paste("data/PRISM/annual/vpdmin_cropped_",1981:2010,"Q",q,".bil", sep=""), raster)), sd)

writeRaster(vpdmin_mn, paste("data/PRISM/norms_1981-2010/vpdmin_mean_cropped_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)
writeRaster(vpdmin_sd, paste("data/PRISM/norms_1981-2010/vpdmin_sd_cropped_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)

}


# FOR LOOP to norm the historical data ...
for(yr in 1895:2022){

# yr=1895

# FOR LOOP over quarters to do the thing ...
	for(q in 1:4){
	
	# q <- 1
	
	# max temp in each quarter
	tmaxQnorm <- (raster(paste("data/PRISM/annual/tmax_cropped_", yr, "Q", q, ".bil", sep="")) - raster(paste("data/PRISM/norms_1981-2010/tmax_mean_cropped_1981-2010_Q", q, ".bil", sep="")))/ raster(paste("data/PRISM/norms_1981-2010/tmax_sd_cropped_1981-2010_Q", q, ".bil", sep=""))
	
	writeRaster(tmaxQnorm, paste("data/PRISM/annual_normed_1981-2010/tmax_normed_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# min temp in each quarter
	tminQnorm <- (raster(paste("data/PRISM/annual/tmin_cropped_", yr, "Q", q, ".bil", sep="")) - raster(paste("data/PRISM/norms_1981-2010/tmin_mean_cropped_1981-2010_Q", q, ".bil", sep="")))/ raster(paste("data/PRISM/norms_1981-2010/tmin_sd_cropped_1981-2010_Q", q, ".bil", sep=""))
	
	writeRaster(tminQnorm, paste("data/PRISM/annual_normed_1981-2010/tmin_normed_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# sum of precip in each quarter
	pptQnorm <- crop((raster(paste("data/PRISM/annual/ppt_cropped_", yr, "Q", q, ".bil", sep="")) - raster(paste("data/PRISM/norms_1981-2010/ppt_mean_cropped_1981-2010_Q", q, ".bil", sep="")))/ raster(paste("data/PRISM/norms_1981-2010/ppt_sd_cropped_1981-2010_Q", q, ".bil", sep="")), MojExt)
	
	writeRaster(pptQnorm, paste("data/PRISM/annual_normed_1981-2010/ppt_normed_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# max VPD in each quarter
	vpdmaxQnorm <- (raster(paste("data/PRISM/annual/vpdmax_cropped_", yr, "Q", q, ".bil", sep="")) - raster(paste("data/PRISM/norms_1981-2010/vpdmax_mean_cropped_1981-2010_Q", q, ".bil", sep="")))/ raster(paste("data/PRISM/norms_1981-2010/vpdmax_sd_cropped_1981-2010_Q", q, ".bil", sep=""))
	
	writeRaster(vpdmaxQnorm, paste("data/PRISM/annual_normed_1981-2010/vpdmax_normed_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)	

	# min VPD in each quarter
	vpdminQnorm <- (raster(paste("data/PRISM/annual/vpdmin_cropped_", yr, "Q", q, ".bil", sep="")) - raster(paste("data/PRISM/norms_1981-2010/vpdmin_mean_cropped_1981-2010_Q", q, ".bil", sep="")))/ raster(paste("data/PRISM/norms_1981-2010/vpdmin_sd_cropped_1981-2010_Q", q, ".bil", sep=""))
	
	writeRaster(vpdminQnorm, paste("data/PRISM/annual_normed_1981-2010/vpdmin_normed_cropped_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)	

	
	} # END loop over quarters

}

# and now I have quarterly historical data for the Mojave to work with for all downstream analysis




