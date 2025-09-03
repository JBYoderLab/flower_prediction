# Using BARTs to model flowering activity as a continuous response
# last used/modified jby, 2025.05.21

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/flower_prediction")

library("tidyverse")
library("cowplot")
library("embarcadero")
library("SoftBart")
# devtools::install_github("theodds/SoftBART")

set.seed(19820604)

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
xvars <- c("ppt.y0q1", "tmax.y0q1", "tmin.y0q1", "vpdmax.y0q1", "vpdmin.y0q1", "ppt.y1q4", "tmax.y1q4", "tmin.y1q4", "vpdmax.y1q4", "vpdmin.y1q4", "ppt.y1q3", "tmax.y1q3", "tmin.y1q3", "vpdmax.y1q3", "vpdmin.y1q3") # weather data, curated --- we don't want weather from later in y0 than flowering season, for instance!

# softBART with all candidate predictors
flrDART <- softbart_regression(paste("prop_flr ~", paste(xvars, collapse="+")), data=flow, test_data=flow, k=2, opts=Opts(num_burn=2000, num_save=1000, num_thin=100)) # this may be slow

invisible(flrDART$forest)
write_rds(flrDART, paste0("output/models/flrDART_", taxon, ".rds"))
flrDART <- read_rds(paste0("output/models/flrDART_", taxon, ".rds"))

plot(flrDART$sigma_mu) # examine MCMC sampling

#-------------------------------------------------------------------------
# variable importance/selection --- single-model PIP

# extract varimp info from model fitted above
variable_selection <- data.frame(varimp=posterior_probs(flrDART)$varimp, post_prob=posterior_probs(flrDART)$post_probs, predictor=xvars) %>% dplyr::arrange(desc(post_prob)) %>% mutate(predictor=factor(predictor, predictor))

# which predictors have posterior inclusion probability > 0.55? (Arbitrary threshold)
variable_selection
variable_selection %>% filter(post_prob > 0.55) # arbitrary threshold --- inspect visually to tweak


# save the top-PIP set in this vector
topXvars <- c("ppt.y1q3", "ppt.y0q1", "tmin.y1q3", "ppt.y1q3", "tmin.y1q4")


#-------------------------------------------------------------------------
# variable importance/selection --- varying model complexity

# fit models with varying complexity
flrDART_compare <- lapply(c(200,100,50,20,10), function(trees) softbart_regression(paste("prop_flr ~", paste(xvars, collapse="+")), data=flow, test_data=flow, num_tree = trees, k=2, opts=Opts(num_burn=2000, num_save=1000, num_thin=100))) # nb training a 200-tree model with these parameters has eta 3h on my M1 MacBook

invisible(flrDART_compare[[1]]$forest)
invisible(flrDART_compare[[2]]$forest)
invisible(flrDART_compare[[3]]$forest)
invisible(flrDART_compare[[4]]$forest)
invisible(flrDART_compare[[5]]$forest)

write_rds(flrDART_compare, paste0("output/models/flrDART_compare_", taxon, ".rds"))
# flrDART_compare <- read_rds(paste0("output/models/flrDART_compare_", taxon, ".rds"))


# doing it for the multi-complexity approach
xvars_ord <- xvars[order(posterior_probs(flrDART_compare[[4]])$post_probs, decreasing=TRUE)] # may tweak this

var_sel_compare <- data.frame(trees=rep(c(200,100,50,20,10), each=15), 
	rbind(data.frame(varimp=posterior_probs(flrDART_compare[[1]])$varimp, post_prob=posterior_probs(flrDART_compare[[1]])$post_probs, predictor=xvars),
	data.frame(varimp=posterior_probs(flrDART_compare[[2]])$varimp, post_prob=posterior_probs(flrDART_compare[[2]])$post_probs, predictor=xvars), data.frame(varimp=posterior_probs(flrDART_compare[[3]])$varimp, post_prob=posterior_probs(flrDART_compare[[3]])$post_probs, predictor=xvars),
	data.frame(varimp=posterior_probs(flrDART_compare[[4]])$varimp, post_prob=posterior_probs(flrDART_compare[[4]])$post_probs, predictor=xvars),
	data.frame(varimp=posterior_probs(flrDART_compare[[5]])$varimp, post_prob=posterior_probs(flrDART_compare[[5]])$post_probs, predictor=xvars))) %>% 
	mutate(predictor=factor(predictor, xvars_ord), trees=factor(trees, c(200,100,50,20,10)))
glimpse(var_sel_compare)

levels(var_sel_compare$predictor) 

levels(var_sel_compare$predictor) <- c("Tmin Y1Q3", "PPT Y0Q1", "PPT Y1Q4", "Tmin Y1Q4", "PPT Y1Q3", "VPDmin Y0Q1", "VPDmin Y1Q4", "VPDmax Y0Q1", "Tmax Y0Q1", "VPDmax Y1Q4", "VPDmax Y1Q3", "Tmax Y1Q3", "VPDmin Y1Q3", "Tmax Y1Q4", "Tmin Y0Q1")

predsel <- ggplot(data=filter(var_sel_compare, trees%in%c(10,20,50,100,200)), aes(x=predictor, y=post_prob, color=trees, group=trees)) +
	geom_line(linewidth=0.5) + geom_point(size=2) +
	labs(x = "Predictor", y = "Posterior inclusion prob.") +
	scale_color_manual(values=c('#ccebc5','#7bccc4','#4eb3d3','#2b8cbe','#08589e'), name="N trees") +
	theme_bw(base_size=12) + theme(axis.text.x=element_text(angle=75, hjust=1), legend.position="inside", legend.position.inside=c(0.15,0.35), legend.key.spacing.y=unit(0.01, "in"))


{cairo_pdf(paste0("output/figures/DART_varimp_", taxon, ".pdf"), width=5, height=4)

predsel

}
dev.off()

# okay that worked GREAT

topXvars <- xvars_ord[1:6]

topXvars <- c("tmin.y1q3", "ppt.y0q1", "ppt.y1q4", "tmin.y1q4", "ppt.y1q3", "vpdmin.y0q1")

#-------------------------------------------------------------------------
# fit a model with top-candidate predictors

flrDARTtop <- softbart_regression(paste("prop_flr ~", paste(topXvars, collapse="+")), data=flow, test_data=flow, num_tree=50, k=2, opts=Opts(num_burn=4000, num_save=1000, num_thin=100)) # this may be slow

plot(flrDARTtop$sigma_mu) # examine MCMC sampling

invisible(flrDARTtop$forest)
write_rds(flrDARTtop, paste0("output/models/flrDARTtop_", taxon, ".rds"))
# flrDARTtop <- read_rds(paste0("output/models/flrDARTtop_", taxon, ".rds"))

# RMSE of the model real quick
plot(colMeans(flrDARTtop$mu_test), flow$prop_flr)
abline(a = 0, b = 1)

rmse(colMeans(flrDARTtop$mu_test), flow$prop_flr) # 0.3470725

#-------------------------------------------------------------------------
# visualize partial effects of individual predictors 

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

glimpse(partials)
partials$pred[partials$pred=="tmin.y1q3"] <- "Tmin Y1Q3"
partials$pred[partials$pred=="ppt.y0q1"] <- "PPT Y0Q1"
partials$pred[partials$pred=="ppt.y1q4"] <- "PPT Y1Q4"
partials$pred[partials$pred=="tmin.y1q4"] <- "Tmin Y1Q4"
partials$pred[partials$pred=="ppt.y1q3"] <- "PPT Y1Q3"
partials$pred[partials$pred=="vpdmin.y0q1"] <- "VPDmin Y0Q1"

partials$predclass <- "Precip"
partials$predclass[grep("Tm", partials$pred)] <- "Temp"
partials$predclass[grep("VPD", partials$pred)] <- "VPD"

# this takes long enough we should save results
write.table(partials, paste0("output/flrDARTtop_partials_", taxon,".csv"), sep=",", col.names=TRUE, row.names=FALSE)
# partials <- read.csv(paste0("output/flrDARTtop_partials_", taxon,".csv"))

partials$pred <- factor(partials$pred, c("Tmin Y1Q3", "PPT Y0Q1", "PPT Y1Q4", "Tmin Y1Q4", "PPT Y1Q3", "VPDmin Y0Q1"))
glimpse(partials)

# scale_fill_manual(values=c('#1f78b4','#a6cee3','#33a02c','#b2df8a'), name="Phenology annotated") + 


# construct a multi-panel figure
library("RColorBrewer")

partplots <- ggplot(partials, aes(x = value, y = mu, fill = predclass)) +
	geom_ribbon(stat = "summary", alpha = 1, fun.min = function(x) quantile(x, 0.025), fun.max = function(x) quantile(x, 0.975)) + 
	geom_line(stat = "summary", fun = mean, color="white") +
	xlab("Predictor value") + ylab("Marginal flowering frequency") +
	scale_fill_manual(values = c("#a6cee3", "#fdbf6f", "#cab2d6")) + 
	facet_wrap("pred", nrow=2, scale="free") +
	theme_bw(base_size=16) + theme(panel.spacing=unit(0.2,"in"), legend.position="none")


{cairo_pdf(paste0("output/figures/flrDARTtop_partials_", taxon, ".pdf"), width=10, height=5)

partplots

}
dev.off()

{cairo_pdf(paste("output/figures/DART_predsel_partials_", taxon, ".pdf", sep=""), width=9, height=4)

ggdraw() + draw_plot(predsel, 0, 0, 0.38, 1) + draw_plot(partplots, 0.38, 0, 0.61, 1) + draw_plot_label(label=c("A", "B"), x=c(0, 0.39), y=1)

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

