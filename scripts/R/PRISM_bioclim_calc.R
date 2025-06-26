# working with PRISM historical monthlys
# Assumes local environment 
# jby 2025.02.27

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
# Download the PRISM monthlies for 1901-1930, 1991-2020

for(yr in c(1901:1930, 1991:2020)){ # note that years are up to you, 1895 is earliest avialable

# yr <- 1902
get_prism_monthlys(type="tmean", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="tmax", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="tmin", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="ppt", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmax", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmin", mon=1:12, year=yr, keepZip=FALSE)

}

# oh yes nice I have these already


#-------------------------------------------------------------------------
# process PRISM data layers into the BIOCLIM variables
# Following O’Donnell, M.S., and Ignizio, D.A., 2012, Bioclimatic predictors for supporting ecological applications in the conterminous United States: U.S. Geological Survey Data Series 691, 10 

stdrast <- rast(pd_to_file(prism_archive_subset("tmean", "monthly", years=1901))) # might want this

# save them as we go
if(!dir.exists("../data/spatial/custom_BIOCLIM")) dir.create("../data/spatial/custom_BIOCLIM")

# I want

# BIO1 = Annual Mean Temperature
AMT_1901_1930 <- app(rast(pd_to_file(prism_archive_subset("tmean", "monthly", years=1901:1930))), fun=mean, na.rm = TRUE)

AMT_1901_1930
plot(AMT_1901_1930)
writeRaster(AMT_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_01_AMT_1901_1930.tiff", overwrite=TRUE) 

AMT_1991_2020 <- app(rast(pd_to_file(prism_archive_subset("tmean", "monthly", years=1991:2020))), fun=mean, na.rm = TRUE)

AMT_1991_2020
plot(AMT_1991_2020)
writeRaster(AMT_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_01_AMT_1991_2020.tiff", overwrite=TRUE) 


# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
MDR_1901_1930 <- app(rast(pd_to_file(prism_archive_subset("tmax", "monthly", years=1901:1930))) - rast(pd_to_file(prism_archive_subset("tmin", "monthly", years=1901:1930))), fun=mean, na.rm = TRUE)

MDR_1901_1930
plot(MDR_1901_1930)
writeRaster(MDR_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_02_MDR_1901_1930.tiff", overwrite=TRUE) 


MDR_1991_2020 <- app(rast(pd_to_file(prism_archive_subset("tmax", "monthly", years=1991:2020))) - rast(pd_to_file(prism_archive_subset("tmin", "monthly", years=1991:2020))), fun=mean, na.rm = TRUE)

MDR_1991_2020
plot(MDR_1991_2020)
writeRaster(MDR_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_02_MDR_1991_2020.tiff", overwrite=TRUE) 



# BIO3 = Isothermality (BIO2/BIO7) (×100)
ITH_1901_1930 <- MDR_1901_1930/TAR_1901_1930*100

ITH_1901_1930
plot(ITH_1901_1930)
writeRaster(ITH_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_03_ITH_1901_1930.tiff", overwrite=TRUE) 


ITH_1991_2020 <- MDR_1991_2020/TAR_1991_2020*100

ITH_1991_2020
plot(ITH_1991_2020)
writeRaster(ITH_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_03_ITH_1991_2020.tiff", overwrite=TRUE) 


# BIO4 = Temperature Seasonality (standard deviation ×100)
TS_1901_1930 <- app(tapp(rast(pd_to_file(prism_archive_subset("tmean", "monthly", years=1901:1930))), index=rep(1:30, each=12), fun=sd, na.rm=TRUE), fun=mean, na.rm=TRUE)

TS_1901_1930
plot(TS_1901_1930)
writeRaster(TS_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_04_TS_1901_1930.tiff", overwrite=TRUE) 


TS_1991_2020 <- app(tapp(rast(pd_to_file(prism_archive_subset("tmean", "monthly", years=1991:2020))), index=rep(1:30, each=12), fun=sd, na.rm=TRUE), fun=mean, na.rm=TRUE)

TS_1991_2020
plot(TS_1991_2020)
writeRaster(TS_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_04_TS_1991_2020.tiff", overwrite=TRUE) 


# BIO5 = Max Temperature of Warmest Month
MTWaM_1901_1930 <- app(tapp(rast(pd_to_file(prism_archive_subset("tmax", "monthly", years=1901:1930))), index=rep(1:30, each=12), fun=max, na.rm=TRUE), fun=mean, na.rm=TRUE)

MTWaM_1901_1930
plot(MTWaM_1901_1930)
writeRaster(MTWaM_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_05_MTWaM_1901_1930.tiff", overwrite=TRUE) 


MTWaM_1991_2020 <- app(tapp(rast(pd_to_file(prism_archive_subset("tmax", "monthly", years=1991:2020))), index=rep(1:30, each=12), fun=max, na.rm=TRUE), fun=mean, na.rm=TRUE)

MTWaM_1991_2020
plot(MTWaM_1991_2020)
writeRaster(MTWaM_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_05_MTWaM_1991_2020.tiff", overwrite=TRUE) 


# BIO6 = Min Temperature of Coldest Month
MTCoM_1901_1930 <- app(tapp(rast(pd_to_file(prism_archive_subset("tmin", "monthly", years=1901:1930))), index=rep(1:30, each=12), fun=min, na.rm=TRUE), fun=mean, na.rm=TRUE)

MTCoM_1901_1930
plot(MTCoM_1901_1930)
writeRaster(MTCoM_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_06_MTCoM_1901_1930.tiff", overwrite=TRUE) 


MTCoM_1991_2020 <- app(tapp(rast(pd_to_file(prism_archive_subset("tmin", "monthly", years=1991:2020))), index=rep(1:30, each=12), fun=min, na.rm=TRUE), fun=mean, na.rm=TRUE)

MTCoM_1991_2020
plot(MTCoM_1991_2020)
writeRaster(MTCoM_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_06_MTCoM_1991_2020.tiff", overwrite=TRUE) 


# BIO7 = Temperature Annual Range (BIO5-BIO6)
TAR_1901_1930 <- MTWaM_1901_1930 - MTCoM_1901_1930

TAR_1901_1930
plot(TAR_1901_1930)
writeRaster(TAR_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_07_TAR_1901_1930.tiff", overwrite=TRUE) 


TAR_1991_2020 <- MTWaM_1991_2020 - MTCoM_1991_2020

TAR_1991_2020
plot(TAR_1991_2020)
writeRaster(TAR_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_07_TAR_1991_2020.tiff", overwrite=TRUE) 


# BIO8 = Mean Temperature of Wettest Quarter
# sum quarterly precip for 1901-1930, take means for each quarter across all 30 years, identify which quarter max
WettestQs_1901_1930 <- app(tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1901:1930))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE), fun=which.max)
plot(WettestQs_1901_1930, main="Wettest quarter, 1901-1930") # cooooool

MTQs_1901_1930 <- tapp(tapp(rast(pd_to_file(prism_archive_subset("tmean", "monthly", years=1901:1930))), index=rep(1:120, each=3), fun=mean, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE)
MTQs_1901_1930

MTWeQ_1901_1930.df <- as.data.frame(c(WettestQs_1901_1930, MTQs_1901_1930)) %>% rename(WQ = `which.max`)
MTWeQ_1901_1930.df$MTWeQ <- apply(as.matrix(MTWeQ_1901_1930.df), 1, function(x) c(x[2:5])[x[1]])

MTWeQ_1901_1930.vect <- rep(NaN, ncell(WettestQs_1901_1930)) # make a NaN vector of length ncell(WettestQs_1901_1930)
MTWeQ_1901_1930.vect[which(!is.na(values(WettestQs_1901_1930)))] <- MTWeQ_1901_1930.df$MTWeQ # values at positions where WettestQs_1901_1930 has non-NaNs

MTWeQ_1901_1930 <- rast(ncol=ncol(WettestQs_1901_1930), nrow=nrow(WettestQs_1901_1930), crs=crs(WettestQs_1901_1930), extent=ext(WettestQs_1901_1930), vals=MTWeQ_1901_1930.vect) # FINALLY convert it into a SpatRaster

MTWeQ_1901_1930
plot(MTWeQ_1901_1930) # glory be
writeRaster(MTWeQ_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_08_MTWeQ_1901_1930.tiff", overwrite=TRUE) 


# same but for 1991-2020
WettestQs_1991_2020 <- app(tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1991:2020))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE), fun=which.max)
plot(WettestQs_1991_2020, main="Wettest quarter, 1991-2020") # cooooool

MTQs_1991_2020 <- tapp(tapp(rast(pd_to_file(prism_archive_subset("tmean", "monthly", years=1991:2020))), index=rep(1:120, each=3), fun=mean, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE)
MTQs_1991_2020

MTWeQ_1991_2020.df <- as.data.frame(c(WettestQs_1991_2020, MTQs_1991_2020)) %>% rename(WQ = `which.max`)
MTWeQ_1991_2020.df$MTWeQ <- apply(as.matrix(MTWeQ_1991_2020.df), 1, function(x) c(x[2:5])[x[1]])

MTWeQ_1991_2020.vect <- rep(NaN, ncell(WettestQs_1991_2020)) # make a NaN vector of length ncell(WettestQs_1991_2020)
MTWeQ_1991_2020.vect[which(!is.na(values(WettestQs_1991_2020)))] <- MTWeQ_1991_2020.df$MTWeQ # values at positions where WettestQs_1991_2020 has non-NaNs

MTWeQ_1991_2020 <- rast(ncol=ncol(WettestQs_1991_2020), nrow=nrow(WettestQs_1991_2020), crs=crs(WettestQs_1991_2020), extent=ext(WettestQs_1991_2020), vals=MTWeQ_1991_2020.vect) # FINALLY convert it into a SpatRaster

MTWeQ_1991_2020
plot(MTWeQ_1991_2020) # glory be
writeRaster(MTWeQ_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_08_MTWeQ_1991_2020.tiff", overwrite=TRUE) 


# BIO9 = Mean Temperature of Driest Quarter
# sum quarterly precip for 1901-1930, take means for each quarter across all 30 years, identify which quarter min
DriestQs_1901_1930 <- app(tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1901:1930))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE), fun=which.min)
plot(DriestQs_1901_1930, main="Driest quarter, 1901-1930") # cooooool

MTQs_1901_1930 # confirm we still have this from above

MTDrQ_1901_1930.df <- as.data.frame(c(DriestQs_1901_1930, MTQs_1901_1930)) %>% rename(WQ = `which.min`)
MTDrQ_1901_1930.df$MTDrQ <- apply(as.matrix(MTDrQ_1901_1930.df), 1, function(x) c(x[2:5])[x[1]])

MTDrQ_1901_1930.vect <- rep(NaN, ncell(DriestQs_1901_1930)) # make a NaN vector of length ncell(WettestQs_1901_1930)
MTDrQ_1901_1930.vect[which(!is.na(values(DriestQs_1901_1930)))] <- MTDrQ_1901_1930.df$MTDrQ # values at positions where WettestQs_1901_1930 has non-NaNs

MTDrQ_1901_1930 <- rast(ncol=ncol(DriestQs_1901_1930), nrow=nrow(DriestQs_1901_1930), crs=crs(DriestQs_1901_1930), extent=ext(DriestQs_1901_1930), vals=MTDrQ_1901_1930.vect) # FINALLY convert it into a SpatRaster

MTDrQ_1901_1930
plot(MTDrQ_1901_1930) # glory be
writeRaster(MTDrQ_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_09_MTDrQ_1901_1930.tiff", overwrite=TRUE) 


# and now same for 1991-2020
DriestQs_1991_2020 <- app(tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1991:2020))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE), fun=which.min)
plot(DriestQs_1991_2020, main="Driest quarter, 1991-2020") # cooooool

MTQs_1991_2020 # confirm we still have this from above

MTDrQ_1991_2020.df <- as.data.frame(c(DriestQs_1991_2020, MTQs_1991_2020)) %>% rename(WQ = `which.min`)
MTDrQ_1991_2020.df$MTDrQ <- apply(as.matrix(MTDrQ_1991_2020.df), 1, function(x) c(x[2:5])[x[1]])

MTDrQ_1991_2020.vect <- rep(NaN, ncell(DriestQs_1991_2020)) # make a NaN vector of length ncell(WettestQs_1991_2020)
MTDrQ_1991_2020.vect[which(!is.na(values(DriestQs_1991_2020)))] <- MTDrQ_1991_2020.df$MTDrQ # values at positions where WettestQs_1991_2020 has non-NaNs

MTDrQ_1991_2020 <- rast(ncol=ncol(DriestQs_1991_2020), nrow=nrow(DriestQs_1991_2020), crs=crs(DriestQs_1991_2020), extent=ext(DriestQs_1991_2020), vals=MTDrQ_1991_2020.vect) # FINALLY convert it into a SpatRaster

MTDrQ_1991_2020
plot(MTDrQ_1991_2020) # glory be
writeRaster(MTDrQ_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_09_MTDrQ_1991_2020.tiff", overwrite=TRUE) 


# BIO10 = Mean Temperature of Warmest Quarter
MTWQ_1901_1930 <- app(MTQs_1901_1930, fun=max)
MTWQ_1901_1930
plot(MTWQ_1901_1930)
writeRaster(MTWQ_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_10_MTWQ_1901_1930.tiff", overwrite=TRUE) 

MTWQ_1991_2020 <- app(MTQs_1991_2020, fun=max)
MTWQ_1991_2020
plot(MTWQ_1991_2020)
writeRaster(MTWQ_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_10_MTWQ_1991_2020.tiff", overwrite=TRUE) 


# BIO11 = Mean Temperature of Coldest Quarter
MTCQ_1901_1930 <- app(MTQs_1901_1930, fun=min)
MTCQ_1901_1930
plot(MTCQ_1901_1930)
writeRaster(MTCQ_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_11_MTCQ_1901_1930.tiff", overwrite=TRUE) 

MTCQ_1991_2020 <- app(MTQs_1991_2020, fun=min)
MTCQ_1991_2020
plot(MTCQ_1991_2020)
writeRaster(MTCQ_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_11_MTCQ_1991_2020.tiff", overwrite=TRUE) 


# BIO12 = Annual Precipitation
AP_1901_1930 <- app(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1901:1930))), index=rep(1:30, each=12), fun=sum, na.rm=TRUE), fun=mean, na.rm=TRUE)

AP_1901_1930
plot(AP_1901_1930)
writeRaster(AP_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_12_AP_1901_1930.tiff", overwrite=TRUE) 


AP_1991_2020 <- app(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1991:2020))), index=rep(1:30, each=12), fun=sum, na.rm=TRUE), fun=mean, na.rm=TRUE)

AP_1991_2020
plot(AP_1991_2020)
writeRaster(AP_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_12_AP_1991_2020.tiff", overwrite=TRUE) 


# BIO13 = Precipitation of Wettest Month
PWeM_1901_1930 <- app(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1901:1930))), index=rep(1:30, each=12), fun=max, na.rm=TRUE), fun=mean, na.rm=TRUE)

PWeM_1901_1930
plot(PWeM_1901_1930)
writeRaster(PWeM_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_13_PWeM_1901_1930.tiff", overwrite=TRUE) 


PWeM_1991_2020 <- app(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1991:2020))), index=rep(1:30, each=12), fun=max, na.rm=TRUE), fun=mean, na.rm=TRUE)

PWeM_1991_2020
plot(PWeM_1991_2020)
writeRaster(PWeM_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_13_PWeM_1991_2020.tiff", overwrite=TRUE) 


# BIO14 = Precipitation of Driest Month
PDrM_1901_1930 <- app(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1901:1930))), index=rep(1:30, each=12), fun=min, na.rm=TRUE), fun=mean, na.rm=TRUE)

PDrM_1901_1930
plot(PDrM_1901_1930)
writeRaster(PDrM_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_14_PDrM_1901_1930.tiff", overwrite=TRUE) 


PDrM_1991_2020 <- app(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1991:2020))), index=rep(1:30, each=12), fun=min, na.rm=TRUE), fun=mean, na.rm=TRUE)

PDrM_1991_2020
plot(PDrM_1991_2020)
writeRaster(PDrM_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_14_PDrM_1991_2020.tiff", overwrite=TRUE) 

# BIO15 = Precipitation Seasonality (Coefficient of Variation)
PS_1901_1930 <- app(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1901:1930))), index=rep(1:30, each=12), fun=sd, na.rm=TRUE), fun=mean, na.rm=TRUE)/(1+AP_1901_1930/12)

PS_1901_1930
plot(PS_1901_1930)
writeRaster(PS_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_15_PS_1901_1930.tiff", overwrite=TRUE) 


PS_1991_2020 <- app(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1991:2020))), index=rep(1:30, each=12), fun=sd, na.rm=TRUE), fun=mean, na.rm=TRUE)/(1+AP_1991_2020/12)

PS_1991_2020
plot(PS_1991_2020)
writeRaster(PS_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_15_PS_1991_2020.tiff", overwrite=TRUE) 


# BIO16 = Precipitation of Wettest Quarter
PWeQ_1901_1930 <- app(tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1901:1930))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE), fun=max)

PWeQ_1901_1930
plot(PWeQ_1901_1930) 
writeRaster(PWeQ_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_16_PWeQ_1901_1930.tiff", overwrite=TRUE) 

PWeQ_1991_2020 <- app(tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1991:2020))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE), fun=max)

PWeQ_1991_2020
plot(PWeQ_1991_2020) 
writeRaster(PWeQ_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_16_PWeQ_1991_2020.tiff", overwrite=TRUE) 


# BIO17 = Precipitation of Driest Quarter
PDrQ_1901_1930 <- app(tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1901:1930))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE), fun=min)

PDrQ_1901_1930
plot(PDrQ_1901_1930) 
writeRaster(PDrQ_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_17_PDrQ_1901_1930.tiff", overwrite=TRUE) 

PDrQ_1991_2020 <- app(tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1991:2020))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE), fun=min)

PDrQ_1991_2020
plot(PDrQ_1991_2020) 
writeRaster(PDrQ_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_17_PDrQ_1991_2020.tiff", overwrite=TRUE) 


# BIO18 = Precipitation of Warmest Quarter
MTQs_1901_1930 # confirm we still have this from above

WarmestQs_1901_1930 <- app(MTQs_1901_1930, fun=which.max)
plot(WarmestQs_1901_1930, main="Warmest quarter, 1901-1930") # cooooool

PQs_1901_1930 <- tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1901:1930))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE)

PWaQ_1901_1930.df <- as.data.frame(c(WarmestQs_1901_1930, PQs_1901_1930)) %>% rename(WQ = `which.max`)
PWaQ_1901_1930.df$PWaQ <- apply(as.matrix(PWaQ_1901_1930.df), 1, function(x) c(x[2:5])[x[1]])

PWaQ_1901_1930.vect <- rep(NaN, ncell(WarmestQs_1901_1930)) # make a NaN vector of length ncell(WettestQs_1901_1930)
PWaQ_1901_1930.vect[which(!is.na(values(WarmestQs_1901_1930)))] <- PWaQ_1901_1930.df$PWaQ # values at positions where WettestQs_1901_1930 has non-NaNs

PWaQ_1901_1930 <- rast(ncol=ncol(DriestQs_1901_1930), nrow=nrow(WarmestQs_1901_1930), crs=crs(WarmestQs_1901_1930), extent=ext(WarmestQs_1901_1930), vals=PWaQ_1901_1930.vect) # FINALLY convert it into a SpatRaster

PWaQ_1901_1930
plot(PWaQ_1901_1930) # glory be
writeRaster(PWaQ_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_18_PWaQ_1901_1930.tiff", overwrite=TRUE) 

# and for 1991-2020
MTQs_1991_2020 # confirm we still have this from above

WarmestQs_1991_2020 <- app(MTQs_1991_2020, fun=which.max)
plot(WarmestQs_1991_2020, main="Warmest quarter, 1901-1930") # cooooool

PQs_1991_2020 <- tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1991:2020))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE)

PWaQ_1991_2020.df <- as.data.frame(c(WarmestQs_1991_2020, PQs_1991_2020)) %>% rename(WQ = `which.max`)
PWaQ_1991_2020.df$PWaQ <- apply(as.matrix(PWaQ_1991_2020.df), 1, function(x) c(x[2:5])[x[1]])

PWaQ_1991_2020.vect <- rep(NaN, ncell(WarmestQs_1991_2020)) # make a NaN vector of length ncell(WettestQs_1991_2020)
PWaQ_1991_2020.vect[which(!is.na(values(WarmestQs_1991_2020)))] <- PWaQ_1991_2020.df$PWaQ # values at positions where WettestQs_1991_2020 has non-NaNs

PWaQ_1991_2020 <- rast(ncol=ncol(DriestQs_1991_2020), nrow=nrow(WarmestQs_1991_2020), crs=crs(WarmestQs_1991_2020), extent=ext(WarmestQs_1991_2020), vals=PWaQ_1991_2020.vect) # FINALLY convert it into a SpatRaster

PWaQ_1991_2020
plot(PWaQ_1991_2020) # glory be
writeRaster(PWaQ_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_18_PWaQ_1991_2020.tiff", overwrite=TRUE) 


# BIO19 = Precipitation of Coldest Quarter
MTQs_1901_1930 # confirm we still have this from above

ColdestQs_1901_1930 <- app(MTQs_1901_1930, fun=which.min)
plot(ColdestQs_1901_1930, main="Coldest quarter, 1901-1930") # cooooool

PQs_1901_1930 <- tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1901:1930))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE)

PCoQ_1901_1930.df <- as.data.frame(c(ColdestQs_1901_1930, PQs_1901_1930)) %>% rename(WQ = `which.min`)
PCoQ_1901_1930.df$PCoQ <- apply(as.matrix(PCoQ_1901_1930.df), 1, function(x) c(x[2:5])[x[1]])

PCoQ_1901_1930.vect <- rep(NaN, ncell(ColdestQs_1901_1930)) # make a NaN vector of length ncell(WettestQs_1901_1930)
PCoQ_1901_1930.vect[which(!is.na(values(ColdestQs_1901_1930)))] <- PCoQ_1901_1930.df$PCoQ # values at positions where WettestQs_1901_1930 has non-NaNs

PCoQ_1901_1930 <- rast(ncol=ncol(DriestQs_1901_1930), nrow=nrow(ColdestQs_1901_1930), crs=crs(ColdestQs_1901_1930), extent=ext(ColdestQs_1901_1930), vals=PCoQ_1901_1930.vect) # FINALLY convert it into a SpatRaster

PCoQ_1901_1930
plot(PCoQ_1901_1930) # glory be
writeRaster(PCoQ_1901_1930, "../data/spatial/custom_BIOCLIM/BIO_19_PCoQ_1901_1930.tiff", overwrite=TRUE) 

# and for 1991-2020
MTQs_1991_2020 # confirm we still have this from above

ColdestQs_1991_2020 <- app(MTQs_1991_2020, fun=which.min)
plot(ColdestQs_1991_2020, main="Coldest quarter, 1901-1930") # cooooool

PQs_1991_2020 <- tapp(tapp(rast(pd_to_file(prism_archive_subset("ppt", "monthly", years=1991:2020))), index=rep(1:120, each=3), fun=sum, na.rm=TRUE), index=rep(1:4, 30), fun=mean, na.rm=TRUE)

PCoQ_1991_2020.df <- as.data.frame(c(ColdestQs_1991_2020, PQs_1991_2020)) %>% rename(WQ = `which.min`)
PCoQ_1991_2020.df$PCoQ <- apply(as.matrix(PCoQ_1991_2020.df), 1, function(x) c(x[2:5])[x[1]])

PCoQ_1991_2020.vect <- rep(NaN, ncell(ColdestQs_1991_2020)) # make a NaN vector of length ncell(WettestQs_1991_2020)
PCoQ_1991_2020.vect[which(!is.na(values(ColdestQs_1991_2020)))] <- PCoQ_1991_2020.df$PCoQ # values at positions where WettestQs_1991_2020 has non-NaNs

PCoQ_1991_2020 <- rast(ncol=ncol(DriestQs_1991_2020), nrow=nrow(ColdestQs_1991_2020), crs=crs(ColdestQs_1991_2020), extent=ext(ColdestQs_1991_2020), vals=PCoQ_1991_2020.vect) # FINALLY convert it into a SpatRaster

PCoQ_1991_2020
plot(PCoQ_1991_2020) # glory be
writeRaster(PCoQ_1991_2020, "../data/spatial/custom_BIOCLIM/BIO_19_PCoQ_1991_2020.tiff", overwrite=TRUE) 





# and now I have quarterly historical data cropped for the species-specific extent for all downstream analysis




