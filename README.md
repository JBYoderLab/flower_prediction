Studying flowering demography with crowd-sourced observations 
==============================================================

Readme updated 2 Sep 2025. 


Project description
-------------------

This repo contains code to (1) use the [iNaturalist](https://www.inaturalist.org) API to download species observations based on phenology annotations, by modifying code from the [`rinat`](https://cran.r-project.org/web/packages/rinat/index.html) package, and (2) model the relationship between flowering and weather using spatially interpolated records from [PRISM](https://prism.oregonstate.edu) and Bayesian additive regression tree methods implemented in [`dbarts`](https://cran.r-project.org/web/packages/dbarts) with utilities from [`embarcadero`](https://github.com/cjcarlson/embarcadero).


![A Joshua tree with, a conical cluster of white-green flowers on one branch and a cluster of green, golf-ball-sized fruits on another](protocol_manual/Joshua_tree_flowering_fruiting.jpeg "A Joshua tree bearing open flowers and mature fruit, in Walker Pass, California")


Contents
--------

Subfolders in the repository and their contents:
 
- `protocol_manual` --- a PDF document describing our protocol for adding phenology annotations to the iNaturalist database, with supporting Markdown and image files
- `scripts` --- all project scripts
	- `R` --- all R scripts
		- `get_inat.R` --- script to load [`rinat`](https://cran.r-project.org/web/packages/rinat/index.html) and modify the `get_inat` function to allow searches on phenology state annotation.
		- `inat_phenology_download.R` --- script to use `get_inat.R` to download phenology-annotated observations using the iNaturalist API, clean them up for downstream analysis, and visualize this initial data. Creates:
			- `data/inat_phenology_data_[taxon code].csv` --- iNat records for the taxon with numeric id `[taxon code]`
			- `data/inat_phenology_data_[taxon code]_cleaned.csv` --- the above, cleaned using utilities in the [`CoordinateCleaner`](https://ropensci.github.io/CoordinateCleaner/)
		- `PRISM_data-management.R` --- downloads monthly climate data by year, crops it to an extent defined for the  taxon with numeric id [taxon code], and summarizes it to quarterlies for downstream work. Creates files in `data/PRISM/annual.[taxon code]`, using functionality from the [`prism`](https://cran.r-project.org/web/packages/prism/index.html) package
		- `inat_phenology_data-management.R` --- organization of data output from `inat_phenology_download.R` into rasterized, binary flowering occurrences, and pairing of these records with PRISM data by location and year. Creates: 
			- `output/flowering_obs_rasterized_[taxon code].csv` (only rasterized flowering frequency records) and 
			- `output/flowering_obs_climate_[taxon code].csv` (the above with PRISM data)
		- Scripts for using a "static" SDM to predict habitat suitability change from the early 20th century (1901--1930) to recent times (1991--2020)
			- `PRISM_bioclim_calc.R` --- use PRISM data to calculate the [19 Bioclim variables](https://pubs.usgs.gov/ds/691/) as averages for the early and recent periods
			- `time_shift_SDM.R` --- script to use GBIF records for the target species, and the early-recent Bioclim averages calculated in `PRISM_bioclim_calc.R`, to train a SDM on recent-period records and predict habitat suitability in the early period
		- Scripts for *binary-response* modeling of flowering activity
			- `phenology_modeling_binary.R` --- modeling annualized, rasterized observations of binary flowering (or no flowering) predicted with weather data using Bayesian additive regression tree (BART) methods. Includes code to perform predictor selection and diagnostic plot(s) and trains a final model with top predictors. Saves modeling objects to `output/models` as `rds` files with the `[taxon code]` numeric id tag; and creates predictor partial-effect plots in `output/figures`
			- `phenology_prediction_binary.R` --- uses the BART model trained in  `phenology_modeling_binary.R` to predict flowering from PRISM data for 1900-present, output as annual spatial layers to `output/models/predictions.[taxon code]`
			- `flowering_years.R` --- analyses of the historic flowering activity predictions from `phenology_prediction_binary.R`, with figures created in `output/figures`
		- Scripts for *continuous-response* modeling of flowering activity
			- `phenology_modeling_continous.R` --- modeling annualized, rasterized observations of flowering frequency predicted with weather data using Bayesian additive regression tree (BART) methods. Includes code to perform predictor selection and diagnostic plot(s) and trains a final model with top predictors. Saves modeling objects to `output/models` as `rds` files with the `[taxon code]` numeric id tag; and creates predictor partial-effect plots in `output/figures`
			- `phenology_prediction_continous.R` --- uses the BART model trained in  `phenology_modeling_continous.R` to predict flowering activity from PRISM data for 1900-present, output as annual spatial layers to `output/models/predictions.[taxon code]`
			- `flowering_continouos.R` --- analyses of the historic flowering activity predictions from `phenology_prediction_continuous.R`, with figures created in `output/figures`

The following subfolders are generated in the course of running the pipeline, but they are not part of the version-controlled repository:

- `data` --- for collecting data in various stages of processing/cleaning
	- `PRISM` --- folder for PRISM data, downloaded and modified with `PRISM_data-management.R`
	- `BClim` --- folder for Bioclim variable normals, calculated from PRISM data in `PRISM_bioclim_calc.R` and used as predictors in `time_shift_SDM.R`
	- iNat observations with phenology annotation, output from `inat_phenology_download.R`, stored as delimited text files with numeric taxon codes in the file names. 
- `output` --- transitional data products, modeling results, analysis and figures
	- `models` --- saved BART models and prediction layers
	- `figures` --- what it says on the tin
	
