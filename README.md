Studying flowering demography with crowd-sourced observations 
==============================================================

Readme updated 27 Mar 2024. 


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
		- `basic_SDM.R` --- script to use iNat records to build a species distribution polygon, which is necessary for later-stage analyses. This is a simplified, rough approach, which may not produce publication-quality results --- you may way to use an independently derived SDM polygon if you have it!
		- `inat_phenology_download.R` --- script to use `get_inat.R` to download phenology-annotated observations using the iNaturalist API, clean them up for downstream analysis, and visualize this initial data. Creates:
			- `data/inat_phenology_data_[taxon code].csv` --- iNat records for the taxon with numeric id `[taxon code]`
			- `data/inat_phenology_data_[taxon code]_cleaned.csv` --- the above, cleaned using utilities in the [`CoordinateCleaner`](https://ropensci.github.io/CoordinateCleaner/)
		- `PRISM_data-management.R` --- downloads monthly climate data by year, crops it to an extent defined for the  taxon with numeric id [taxon code], and summarizes it to quarterlies for downstream work. Creates files in `data/PRISM/annual.[taxon code]`, using functionality from the [`prism`](https://cran.r-project.org/web/packages/prism/index.html) package
		- `inat_phenology_data-management.R` --- organization of data output from `inat_phenology_download.R` into rasterized, binary flowering occurrences, and pairing of these records with PRISM data by location and year. Creates: 
			- `output/flowering_obs_rasterized_[taxon code].csv` (only rasterized binary flowering records) and 
			- `output/flowering_obs_climate_[taxon code].csv` (the above with PRISM data)
		- `phenology_modeling.R` --- modeling annualized, rasterized observations of flowering (or no flowering) predicted with weather data using Bayesian additive regression tree (BART) methods. Saves modeling objects to `output/BART` as `rds` files with the `[taxon code]` numeric id tag
		- `phenology_prediction.R` --- uses the BART model trained in  `phenology_modeling.R` to predict flowering from PRISM data for 1900-present, output as annual spatial layers to `output/BART/predictions.[taxon code]`
		- `flowering_years.R` --- analysis of the historic flowering predictions from `phenology_modeling.R`

The following subfolders are generated in the course of running the pipeline, but they are not part of the version-controlled repository:

- `data` --- for collecting data in various stages of processing/cleaning
	- `PRISM` --- folder for PRISM data, downloaded and modified with `PRISM_data-management.R`
	- iNat observations with phenology annotation, output from `inat_phenology_download.R`, stored as delimited text files with numeric taxon codes in the file names. 
- `output` --- transitional data products, modeling results, analysis and figures
	- `BART` --- saved BART models and prediction layers
	- `figures` --- what it says on the tin
	
	
Usage
-----

Some of these steps can be skipped once data is downloaded/organized, but to perform the full analysis for the first time, use the scripts in this order:

1. First use `inat_phenology_download.R` to download iNaturalist observations with phenophase annotated;
2. Use `PRISM_data-management.R` to download spatially interpolated weather data at 4km resolution and process it into quarterly aggregates;
3. Use `inat_phenology_data-management.R` to match the iNat observations to weather results for the years leading up to each observation;
3.5 (optional) You will want a species range polygon to do some of the downstream analysis --- a very simplified version of this is possible with the data you've already got on hand, and `basic_SDM.R` provides code for that;
4. Use `phenology_modeling.R` to evaluate and train a BART model predicting flowering status with weather data; 
5. Use `phenology_prediction.R` with the resulting model to predict what flowering was like in years when we have weather data but no iNaturalist observations;
6. Finally, use `flowering_years.R` to visualize and analyze the reconstructed flowering history

