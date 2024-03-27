Studying flowering demography with crowd-sourced observations 
==============================================================

Readme updated 26 Mar 2024. 


Project description
-------------------

This repo contains code to (1) use the [iNaturalist](https://inaturalist.org) API to download species observations based on phenology annotations, by modifying code from the [`rinat`](https://cran.r-project.org/web/packages/rinat/index.html) package, and (2) model the relationship between flowering and weather using spatially interpolated records from [PRISM](https://prism.oregonstate.edu) and Bayesian additive regression tree methods implemented in [`embarcadero`](https://github.com/cjcarlson/embarcadero).


![A Joshua tree with, a conical cluster of white-green flowers on one branch and a cluster of green, golf-ball-sized fruits on another](protocol_manual/Joshua_tree_flowering_fruiting.jpeg "A Joshua tree bearing open flowers and mature fruit, in Walker Pass, California")


Contents
--------

Subfolder names and contents:
 
- `data` --- raw data, as output by the `inat_phenology_download.R` script
- `output` --- processed data, figures, etc
- `protocol_manual` --- a PDF document describing our protocol for adding phenology annotations to the iNaturalist database, with supporting Markdown and image files
- `scripts` --- all project scripts
	- `R` --- all R scripts, notably
		- `get_inat.R` --- script to load [`rinat`](https://cran.r-project.org/web/packages/rinat/index.html) and modify the `get_inat` function to allow searches on phenology state annotation.
		- `inat_phenology_download.R` --- script to use `get_inat.R` to download phenology-annotated observations using the iNaturalist API
		- `PRISM_data-management.R` --- downloads monthly climate data by year, summarizes to quarterly, normalizes to 1981-2010, and crops to the Mojave extent to build a repository for downstream work
		- `inat_phenology_data-management.R` --- organization of data output from `inat_phenology_download.R`, and pairing with PRISM data with functionality from the [`prism`](https://cran.r-project.org/web/packages/prism/index.html) package.
		- `phenology_modeling.R` --- modeling annualized, rasterized observations of flowering (or no flowering) predicted with weather data using Bayesian additive regression tree (BART) methods
		- `phenology_prediction.R` --- fits random intercept models based on the best-fit from `phenology_modeling.R`, with year as the random effect, uses that to predict flowering from PRISM data for 1900-present
		- `historic_flowering_analysis.R` --- analysis of the historic flowering prediction from `phenology_modeling.R`
- `data` --- raw data
	- `inat_phenology_data.csv` --- cleaned iNat observations with phenology annotation, output from `inat_phenology_download.R`
	- `PRISM` --- folder for PRISM data, downloaded and modified with `PRISM_data-management.R`
- `output` --- transitional data products, modeling results, analysis and figures
	- `flowering_obs_rasterized.csv` --- flowering observations, rendered binary (flowering/no) and rasterized to the 4km PRISM resolution, output from `inat_phenology_data-management.R`
	- `flowering_obs_climate.csv` --- rasterized binary flowering observations merged with predictors generated from PRISM weather data, output from `inat_phenology_data-management.R`
	- `BART` --- saved BART models, prediction layers
	- `figures` --- what it says on the tin
	
	
Usage
-----

Some of these steps can be skipped once data is downloaded/organized, but to perform the full analysis for the first time, use the scripts in this order:

1. First use `inat_phenology_download.R` to download all of the iNaturalist observations with flowering status annotated;
2. Then use `PRISM_data-management.R` to download spatially interpolated weather data at 4km resolution and process it into quarterly aggregates, normalized to 1981-2010;
3. Then, use `inat_phenology_data-management.R` to match iNat observations to weather results for the years leading up to each observation;
4. Finally, use `phenology_modeling.R` to builds a BART model predicting flowering status with weather data, which are stored in `output/flowering_obs_climate_normed.csv`; and use `phenology_prediction.R` with the resulting model to predict what flowering was like in years when we have weather data but no iNaturalist observations.

