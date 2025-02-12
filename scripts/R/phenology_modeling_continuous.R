# Using BARTs to model flowering activity as a continuous response
# last used/modified jby, 2025.02.05

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/flower_prediction")

library("tidyverse")
library("embarcadero")
library("SoftBart")
# devtools::install_github("theodds/SoftBART")

#-----------------------------------------------------------
# initial file loading

# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia = 57250

flow <- read.csv(paste("output/flowering_freq_climate_", taxon, ".csv", sep="")) %>% filter(!is.na(ppt.y0q1)) # flowering/not flowering, biologically-informed candidate predictors

dim(flow)
glimpse(flow)

table(flow$year, flow$prop_flr>0)
hist(flow$prop_flr)

#-------------------------------------------------------------------------
# fit a first DART model

if(!dir.exists("output/models")) dir.create("output/models")

# predictors
xvars <- c("n_obs", "ppt.y0q1", "tmax.y0q1", "tmin.y0q1", "vpdmax.y0q1", "vpdmin.y0q1", "ppt.y1q4", "tmax.y1q4", "tmin.y1q4", "vpdmax.y1q4", "vpdmin.y1q4", "ppt.y1q3", "tmax.y1q3", "tmin.y1q3", "vpdmax.y1q3", "vpdmin.y1q3") # weather data, curated --- we don't want weather from later in y0 than flowering season, for instance!

# softBART with all candidate predictors
flrDART <- softbart_regression(paste("prop_flr ~", paste(xvars, collapse="+")), data=flow, test_data=flow, k=2, opts=Opts(num_burn=2000, num_save=1000, num_thin=100)) # this may be slow

invisible(flrDART$forest)
write_rds(flrDART, paste0("output/models/flrDART_", taxon, ".rds"))
flrDART <- read_rds(paste0("output/models/flrDART_", taxon, ".rds"))

plot(flrDART$sigma_mu) # examine MCMC sampling

#-------------------------------------------------------------------------
# variable importance

variable_selection <- data.frame(varimp=posterior_probs(flrDART)$varimp, post_prob=posterior_probs(flrDART)$post_probs, predictor=xvars) %>% dplyr::arrange(desc(post_prob)) %>% mutate(predictor=factor(predictor, predictor))

# which predictors have posterior inclusion probability > 0.55? (Arbitrary threshold)
variable_selection
variable_selection %>% filter(post_prob > 0.55) # arbitrary threshold --- inspect visually to tweak

{png(paste0("output/figures/flrDART_", taxon, "_variable_selection.png"), width=1000, height=500)

ggplot() +
	geom_point(data=variable_selection, aes(x=predictor, y=post_prob, color=post_prob<0.8), size=4) +
	geom_hline(yintercept=0.8, linetype=2) + 
	labs(x = "Predictor", y = "Posterior inclusion prob") +
	theme_bw(base_size=24) + theme(axis.text.x=element_text(angle=75, hjust=1), legend.position="none")

}
dev.off()

topXvars <- variable_selection %>% filter(post_prob>0.9) %>% .$predictor %>% as.character
topXvars

#-------------------------------------------------------------------------
# fit a model with top-candidate predictors

flrDARTtop <- softbart_regression(paste("prop_flr ~", paste(topXvars, collapse="+")), data=flow, test_data=flow, k=2, opts=Opts(num_burn=4000, num_save=1000, num_thin=100)) # this may be slow

plot(flrDARTtop$sigma_mu) # examine MCMC sampling

invisible(flrDARTtop$forest)
write_rds(flrDARTtop, paste0("output/models/flrDARTtop_", taxon, ".rds"))
flrDARTtop <- read_rds(paste0("output/models/flrDARTtop_", taxon, ".rds"))


#-------------------------------------------------------------------------
#  visualize partial effects of individual predictors 

# do partial regressions for each top predictor ...
partials <- matrix(0,0,4)
colnames(partials) <- c("pred", "value", "mu", "sample")

for(p in topXvars){

# p <- ppt.y1q4

grid_pred <- seq(from = min(flow[,p]), to = max(flow[,p]), length = 50)
pdf_pred <- partial_dependence_regression(flrDARTtop, flow, p, grid_pred)

pdf_out <- data.frame(pred=p, value=pdf_pred$pred_df[,3], mu=pdf_pred$pred_df$mu, sample=pdf_pred$pred_df$sample)

partials <- rbind(partials, pdf_out)

}
# this takes long enough we should save results
write.table(partials, "output/flrDARTtop_partials.csv", sep=",", col.names=TRUE, row.names=FALSE)

# construct a multi-panel figure
{png("output/figures/flrDARTtop_partial_ppt.y1q4.png", width=750, height=500)

ggplot(partials, aes(x = value, y = mu)) +
	geom_line(stat = "summary", fun = mean) +
	geom_ribbon(stat = "summary", alpha = 0.3, fun.min = function(x) quantile(x, 0.025), fun.max = function(x) quantile(x, 0.975)) + xlab("Predictor value") + ylab("Marginal flowering frequency") +
	facet_wrap("pred", nrow=2, scale="free_x") +
	theme_bw(base_size=18)

}
dev.off()



# SPARTIALS ----------------------------------------------- ???

if(!dir.exists("output/models/BART_spartials")) dir.create("output/BART/BART_spartials", recursive=TRUE) # make sure there's a folder to write to!

# read back in, if necessary
flr.mod <- read_rds(paste("output/BART/bart.model.", taxon, ".rds", sep=""))
preds <- attr(flr.mod$fit$data@x, "term.labels")

prism_temp_rast <- raster(paste("data/PRISM/annual.", taxon, "/tmax_cropped_2010Q1.bil", sep="")) # raster grid base

# LOOP over years in the training data
for(y0 in sort(unique(flow$year))){

# y0 <- 2021

# parse the year into the corresponding quarterly predictor layers
vars <- gsub("(\\w+)\\..+", "\\1", preds)
yrs <- c(y0=y0, y1=y0-1)[gsub("\\w+\\.(y\\d).+", "\\1", preds)]
qs <- as.numeric(gsub("\\w+\\.y\\dq(\\d)", "\\1", preds))

predfiles <- paste("data/PRISM/annual.", taxon, "/", vars, "_cropped_", yrs, "Q", qs, ".bil", sep="")
predbrick <- brick(lapply(predfiles, function(x) resample(raster(x), prism_temp_rast)))
names(predbrick) <- preds

# estimate the spatial partial effects as a new raster layer
spYr <- spartial(flr.mod, predbrick, x.vars=preds)

# write out a plot (may need to pilot and adjust page dimensions)
{cairo_pdf(paste("output/figures/BART_spartials_", taxon, "_", y0, ".pdf", sep=""), width=9, height=6.5)
	plot(spYr)
}
dev.off()

# write out a raster file
writeRaster(spYr, paste("output/BART/BART_spartials/BART_spartials_", taxon, "_", y0,".grd", sep=""), overwrite=TRUE) 

# status updates to the terminal
cat("Done with spartials for", y0, "\n")

}
# END loop over years

