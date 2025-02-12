# Using BARTs to model flowering activity
# best run on MAJEL
# last used/modified jby, 2024.03.27

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/flower_prediction")

library("tidyverse")
library("embarcadero")

#-----------------------------------------------------------
# initial file loading

# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia = 57250

flow <- read.csv(paste("output/flowering_obs_climate_", taxon, ".csv", sep="")) %>% filter(!is.na(ppt.y0q1)) # flowering/not flowering, biologically-informed candidate predictors

dim(flow)
glimpse(flow)

table(flow$year, flow$flr)

#-------------------------------------------------------------------------
# fit candidate BART models

if(!dir.exists("output/BART")) dir.create("output/BART")

# predictors
xnames <- c("ppt.y0q1", "tmax.y0q1", "tmin.y0q1", "vpdmax.y0q1", "vpdmin.y0q1", "ppt.y1q4", "tmax.y1q4", "tmin.y1q4", "vpdmax.y1q4", "vpdmin.y1q4", "ppt.y1q3", "tmax.y1q3", "tmin.y1q3", "vpdmax.y1q3", "vpdmin.y1q3") # weather data, curated --- we don't want weather from later in y0 than flowering season, for instance!

# VARIMP variable importance across the whole predictor set .........
flow.varimp <- varimp.diag(y.data=as.numeric(flow[,"flr"]), x.data=flow[,xnames], ri.data=flow[,"year"])

write_rds(flow.varimp, file=paste("output/BART/bart.varimp.", taxon, ".rds", sep="")) # save varimp() results
# flow.varimp <- read_rds(file=paste("output/BART/bart.varimp.", taxon, ".rds", sep=""))

# generate a better-organized varimp() figure
flow.varimp$data <- flow.varimp$data |> mutate(trees = factor(trees, c(10,20,50,100,150,200)))

flow.varimp$labels$group <- "Trees"
flow.varimp$labels$colour <- "Trees"

label_parse <- function(breaks){ parse(text=breaks) } # need this, for reasons

{cairo_pdf(paste("output/figures/varimp_", taxon, ".pdf", sep=""), width=6, height=5)

flow.varimp + scale_x_discrete(label=label_parse) + 
theme_bw(base_size=12) +
theme(legend.position="inside", legend.position.inside=c(0.8, 0.7), axis.text.x=element_text(angle=45, hjust=1))  # okay nice

}
dev.off()

# fill this in based on results of varimp.diag()
preds <- c("ppt.y1q3", "vpdmax.y0q1", "tmin.y1q4", "vpdmin.y1q4", "tmax.y0q1", "ppt.y1q4")

# STEPWISE model training to confirm selection ......................
flr.mod.step <- bart.step(y.data=as.numeric(flow[,"flr"]), x.data=flow[,xnames], ri.data=flow[,"year"], full=FALSE, quiet=TRUE)

invisible(flr.mod.step$fit$state)
write_rds(flr.mod.step, file=paste("output/BART/bart.step.models.", taxon, ".rds", sep="")) # save stepwise model
# flr.mod.step <- read_rds(file=paste("output/BART/bart.step.models.", taxon, ".rds", sep=""))

summary(flr.mod.step) # note predictors selected this way, compare to varimp() output


# STANDARD model with varimp() selection ..................
flr.mod <- bart(y.train=as.numeric(flow[,"flr"]), x.train=flow[,preds], keeptrees=TRUE)

invisible(flr.mod$fit$state)
write_rds(flr.mod, file=paste("output/BART/bart.model.", taxon, ".rds", sep="")) # save model
# flr.mod <- read_rds(paste("output/BART/bart.model.", taxon, ".rds", sep=""))

summary(flr.mod) # AUC reflects classification accuracy, how's that look?

p <- partial(flr.mod, preds, trace=FALSE, smooth=5) # visualize partials
varimp(flr.mod)

write_rds(p, file=paste("output/BART/bart.model.partials.", taxon, ".rds", sep=""))
# p <- read_rds(paste("output/BART/bart.model.partials.", taxon, ".rds", sep=""))


# random intercept model ..................................
flr.RImod <- rbart_vi(as.formula(paste(paste('flr', paste(preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = flow,
	group.by = flow[,'year'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE)

summary(flr.RImod)

invisible(flr.RImod$fit[[1]]$state) # MUST do this to save
write_rds(flr.RImod, file=paste("output/BART/bart.RImodel.", taxon, ".rds", sep="")) # save model
# flr.RImod <- read_rds(paste("output/BART/bart.RImodel.", taxon, ".rds", sep=""))

# visualize RI estimates --- is there a temporal trend? If so, keep RI, otherwise stick with non-RI model

{cairo_pdf(paste("output/figures/BART_RI_estimates_", taxon, ".pdf", sep=""), width=5, height=3.5)
plot.ri(flr.RImod, temporal=TRUE) + labs(title="Random intercept effects of observation year", x="Observation year") + theme_bw(base_size=12) + theme(axis.text.x=element_text(angle=0, size=12))
}
dev.off()

#-------------------------------------------------------------------------
# plot the predictors in raw observations ...

pred.plot <- flow %>% dplyr::select(year, flr, all_of(preds)) %>% pivot_longer(all_of(preds), names_to="Predictor", values_to="Value")


{cairo_pdf(file=paste("output/figures/mod_best_predictors_", taxon, ".pdf", sep=""), width=6, height=2.5)

ggplot(pred.plot, aes(x=flr, y=Value)) + geom_jitter(alpha=0.25, size=0.25, color="#a6cee3") + geom_boxplot(alpha=0.5, width=0.5) + 

#scale_color_manual(values=park_palette("JoshuaTree")[c(1,7,7)], guide=FALSE) +
#scale_fill_manual(values=park_palette("JoshuaTree")[c(1,7,7)], guide=FALSE) +

facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value") + theme_bw(base_size=8) + theme(legend.position="none")

}
dev.off()

#-------------------------------------------------------------------------
# Partials and spartials in example years

# read back in, if necessary
flr.mod <- read_rds(paste("output/BART/bart.model.", taxon, ".rds", sep=""))
preds <- attr(flr.mod$fit$data@x, "term.labels")


# PARTIALS ------------------------------------------------
p <- read_rds(paste("output/BART/bart.model.partials.", taxon, ".rds", sep=""))

partvals <- NULL # empty object to hold data

# reorganize raw data underlying the partials
for(part in 1:length(preds)){
	partvals <- rbind(partvals, data.frame(predictor=preds[part], p[[part]]$data))
}

glimpse(partvals)

# generate a figure .............................

# This may need tweaking, to 
# - provide more readable predictor labels
# - reorient or reorganize if the predictor count is not a multiple of three
{cairo_pdf(paste("output/figures/predictor_partials_", taxon, ".pdf", sep=""), width=5, height=6)

ggplot(partvals) + 
	geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#a6cee3") + 
	geom_line(aes(x=x, y=med), color="white") + 
	facet_wrap("predictor", nrow=3, labeller="label_parsed", scale="free") + 
	labs(y="Marginal Pr(Flowers)", x="Predictor value") + 
	theme_bw() + theme(panel.spacing=unit(0.2,"in"))

}
dev.off()

# SPARTIALS -----------------------------------------------

if(!dir.exists("output/BART/BART_spartials")) dir.create("output/BART/BART_spartials", recursive=TRUE) # make sure there's a folder to write to!

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

