# Using BARTs to model flowering activity
# best run on MAJEL
# last used/modified jby, 2024.02.27

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/flowering_prediction")

library("tidyverse")
library("embarcadero")
library("ggdark")

source("../shared/Rscripts/base_graphics.R")

#-----------------------------------------------------------
# initial file loading

# flow <- read.csv("output/flowering_obs_climate.csv") # flowering/not flowering, gridded and annualized
# flow <- read.csv("output/flowering_obs_climate_normed.csv") # flowering/not flowering, gridded and annualized and normed

# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia = 57250

flow <- read.csv(paste("output/flowering_obs_climate_", taxon, ".csv", sep="")) %>% filter(!is.na(ppt.y0q1)) # flowering/not flowering, biologically-informed candidate predictors, subspecies id'd

dim(flow)
glimpse(flow)


#-------------------------------------------------------------------------
# fit candidate BART models, stepwise

if(!dir.exists("output/BART")) dir.create("output/BART")

# predictors
xnames <- c("ppt.y0q1", "tmax.y0q1", "tmin.y0q1", "vpdmax.y0q1", "vpdmin.y0q1", "ppt.y1q4", "tmax.y1q4", "tmin.y1q4", "vpdmax.y1q4", "vpdmin.y1q4", "ppt.y1q3", "tmax.y1q3", "tmin.y1q3", "vpdmax.y1q3", "vpdmin.y1q3") # weather data, curated

# variable importance across the whole predictor set
flow.varimp <- varimp.diag(y.data=as.numeric(flow[,"flr"]), x.data=flow[,xnames], ri.data=flow[,"year"])

# enter based on results of varimp.diag()
preds <- c("ppt.y1q3", "vpdmin.y1q3", "tmin.y1q4", "ppt.y1q4", "tmax.y0q1", "vpdmin.y1q4", "vpdmax.y1q3")

write_rds(flow.varimp, file=paste("output/BART/bart.varimp.", taxon, ".rds", sep=""))

# variable selection, stepwise
flr.mod.step <- bart.step(y.data=as.numeric(flow[,"flr"]), x.data=flow[,xnames], ri.data=flow[,"year"], full=FALSE, quiet=TRUE)

invisible(flr.mod.step$fit$state)
write_rds(flr.mod.step, file=paste("output/BART/bart.step.models.", taxon, ".rds", sep=""))

summary(flr.mod.step) # note predictors selected this way, compare to varimp() output


# standard model with varimp() selection ..................
flr.mod <- bart(y.train=as.numeric(flow[,"flr"]), x.train=flow[,preds], keeptrees=TRUE)

summary(flr.mod)

invisible(flr.mod$fit$state)
write_rds(flr.mod, file=file=paste("output/BART/bart.model.", taxon, ".rds", sep=""))
# jotr.mod <- read_rds(paste("output/BART/bart.model.", taxon, ".rds", sep=""))

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
write_rds(flr.RImod, file=paste("output/BART/bart.RImodel.", taxon, ".rds", sep="")) # write out for downstream use

# flr.mod <- read_rds(paste("output/BART/bart.RImodel.", taxon, ".rds", sep=""))


#-------------------------------------------------------------------------
# plot the predictors in raw observations ...

pred.plot <- flow %>% dplyr::select(year, flr, all_of(preds)) %>% pivot_longer(all_of(preds), names_to="Predictor", values_to="Value")


{cairo_pdf(file=paste("output/figures/mod_best_predictors_", taxon, ".pdf", sep=""), width=6, height=2.5)

ggplot(pred.plot, aes(x=flr, y=Value)) + geom_jitter(alpha=0.25, size=0.25) + geom_boxplot(alpha=0.5, aes(color=Predictor, fill=Predictor), width=0.5) + 

#scale_color_manual(values=park_palette("JoshuaTree")[c(1,7,7)], guide=FALSE) +
#scale_fill_manual(values=park_palette("JoshuaTree")[c(1,7,7)], guide=FALSE) +

facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value") + theme_minimal(base_size=8) + theme(legend.position="none")

}
dev.off()


