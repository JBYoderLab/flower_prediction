# Using BARTs to model flowering activity
# best run on MAJEL
# last used/modified jby, 2025.07.09

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/flower_prediction")

library("tidyverse")
library("embarcadero")
library("cowplot")

set.seed(19820604)

#-----------------------------------------------------------
# initial file loading

# set parameters as variables
taxon <- 53405 # toyon!
# Prunus ilicifolia = 57250

flow <- read.csv(paste("output/flowering_freq_climate_", taxon, ".csv", sep="")) %>% filter(!is.na(ppt.y0q1)) %>% mutate(flr = prop_flr == 0) # flowering/not flowering, biologically-informed candidate predictors

dim(flow)
glimpse(flow)

hist(flow$prop_flr) # check the cutoff I've set for binary flowering

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
var_sel_compare <- flow.varimp$data |> mutate(trees = factor(trees, c(200, 150, 100, 50, 20, 10)), variable = factor(variable, xnames[order(filter(flow.varimp$data, trees==10)$imp, decreasing=TRUE)]))

levels(var_sel_compare$variable) <- c("PPT Y1Q3", "VPDmin Y1Q4", "VPDmax Y0Q1", "Tmin Y1Q4", "PPT Y1Q4", "Tmin Y1Q3", "PPT Y0Q1", "VPDmin Y0Q1", "Tmax Y0Q1", "VPDmin Y1Q3", "VPDmax Y1Q3", "Tmax Y1Q4", "VPDmax Y1Q4", "Tmax Y1Q3", "Tmin Y0Q1")

# colors for this: '#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#08589e'

predsel <- ggplot(data=filter(var_sel_compare, trees%in%c(10,20,50,100,200)), aes(x=variable, y=imp, color=trees, group=trees)) +
	geom_line(linewidth=0.5) + geom_point(size=2) +
	labs(x = "Predictor", y = "Importance (prop. splits)") +
	scale_color_manual(values=c('#ccebc5','#7bccc4','#4eb3d3','#2b8cbe','#08589e'), name="N trees") +
	theme_bw(base_size=12) + theme(axis.text.x=element_text(angle=75, hjust=1), legend.position="inside", legend.position.inside=c(0.85,0.65), legend.key.spacing.y=unit(0.01, "in"))


{cairo_pdf(paste0("output/figures/varimp_", taxon, ".pdf"), width=5, height=4)

predsel

}
dev.off()

# fill this in based on results of varimp.diag()
topX <- c("ppt.y1q3", "vpdmin.y1q4", "vpdmax.y0q1", "tmin.y1q4", "ppt.y1q4", "tmin.y1q3")


# STANDARD model with varimp() selection ..................
flr.mod <- bart(y.train=as.numeric(flow[,"flr"]), x.train=flow[,topX], keeptrees=TRUE)

invisible(flr.mod$fit$state)
write_rds(flr.mod, file=paste("output/BART/bart.model.", taxon, ".rds", sep="")) # save model
# flr.mod <- read_rds(paste("output/BART/bart.model.", taxon, ".rds", sep=""))

summary(flr.mod) # AUC reflects classification accuracy, how's that look?

mod_valid <- summary(flr.mod)$data %>% dplyr::select(fitted, observed) %>% mutate(type="Training data", classified=fitted>0.1352385)

rmse(mod_valid$classified, mod_valid$observed) # RMSE = 0.5807869


p <- partial(flr.mod, topX, trace=FALSE, smooth=5) # visualize partials
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
# Partials and spartials in example years

# read back in, if necessary
flr.mod <- read_rds(paste("output/BART/bart.model.", taxon, ".rds", sep=""))
topX <- attr(flr.mod$fit$data@x, "term.labels")


# PARTIALS ------------------------------------------------
p <- read_rds(paste("output/BART/bart.model.partials.", taxon, ".rds", sep=""))

partvals <- NULL # empty object to hold data

# reorganize raw data underlying the partials
for(part in 1:length(topX)){
	partvals <- rbind(partvals, data.frame(predictor=topX[part], p[[part]]$data))
}

glimpse(partvals)

partvals$predictor <- factor(partvals$predictor, topX)
levels(partvals$predictor) <- c("PPT Y1Q3", "VPDmin Y1Q4", "VPDmax Y0Q1", "Tmin Y1Q4", "PPT Y1Q4", "Tmin Y1Q3")
partvals$predtype <- NA
partvals$predtype[grepl("PPT", partvals$predictor)] <- "precip"
partvals$predtype[grepl("Tm", partvals$predictor)] <- "temp"
partvals$predtype[grepl("VPD", partvals$predictor)] <- "vpd"

table(partvals$predtype)

# generate a figure .............................

# This may need tweaking, to 
# - provide more readable predictor labels
# - reorient or reorganize if the predictor count is not a multiple of three

# color set: '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'

partplots <- ggplot(partvals) + 
	geom_ribbon(aes(x=x, ymin=q05, ymax=q95, fill=predtype)) + 
	geom_line(aes(x=x, y=med), color="white") + 
	facet_wrap("predictor", nrow=2, scale="free") + 
	scale_fill_manual(values = c("#a6cee3", "#fdbf6f", "#cab2d6")) +
	labs(y="Marginal Pr(Flowers)", x="Predictor value") + 
	theme_bw() + theme(panel.spacing=unit(0.2,"in"), legend.position="none")


{cairo_pdf(paste("output/figures/predictor_partials_", taxon, ".pdf", sep=""), width=8, height=5)

partplots

}
dev.off()

# predictor selection and partials together

{cairo_pdf(paste("output/figures/BART_predsel_partials_", taxon, ".pdf", sep=""), width=9, height=4)

ggdraw() + draw_plot(predsel, 0, 0, 0.39, 1) + draw_plot(partplots, 0.39, 0, 0.6, 1) + draw_plot_label(label=c("A", "B"), x=c(0, 0.39), y=1)

}
dev.off()



