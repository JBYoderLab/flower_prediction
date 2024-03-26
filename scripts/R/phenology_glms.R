# Using glms to model Joshua tree flowering
# best run on MAJEL
# last used/modified jby, 2022.07.30

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("/Volumes/GoogleDrive/Other computers/My MacBook Pro 2020/Documents/Academic/Active_projects/Jotr_phenology")
# setwd("~/Jotr_phenology-main")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")
library("lme4")
library("MuMIn")

library("embarcadero")

#-----------------------------------------------------------
# initial file loading

# flow <- read.csv("output/flowering_obs_climate.csv") # flowering/not flowering, gridded and annualized
# flow <- read.csv("output/flowering_obs_climate_normed.csv") # flowering/not flowering, gridded and annualized and normed

flow <- read.csv("output/flowering_obs_climate_v2_subsp.csv") # flowering/not flowering, biologically-informed candidate predictors, subspecies id'd

dim(flow)
glimpse(flow)

# variant datasets -- dealing with the second flowering in 2019
flow2 <- flow %>% filter(year!=2019.5) # drop the weird observations
flow3 <- flow
flow3$year[flow3$year==2019.5] <- 2019 # or merge 2019.5 into 2019?

glimpse(flow2)

# split by subspecies
# swap input datasets to change --- current most trustworthy is flow2, ignoring 2019.5
yuja <- filter(flow2, type=="Eastern") 
yubr <- filter(flow2, type=="Western")

glimpse(yuja) # 1,048 obs
glimpse(yubr) # 1,280 obs

#-------------------------------------------------------------------------
# fit candidate glms, in various combinations of predictors

# predictors (climate only)
pnames <- c("pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1") 
tnames <- c("tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1") 
vnames <- c("vpdmaxW0", "vpdminW0", "vpdmaxW0vW1", "vpdminW0vW1")

# going to compare all individual weather variables; all P+T, P+V, and T+V combos ...
PT <- c(sapply(pnames, function(x) sapply(tnames, function(y) paste(c(x,y), collapse=' + '))))
PT # wheeeeeew
PV <- c(sapply(pnames, function(x) sapply(vnames, function(y) paste(c(x,y), collapse=' + '))))
TV <- c(sapply(tnames, function(x) sapply(vnames, function(y) paste(c(x,y), collapse=' + '))))
# ... and all P+T+V combos
PTV <- c(sapply(pnames, function(x) sapply(tnames, function(y) sapply(vnames, function(z) paste(c(x,y,z), collapse=' + ')))))
PTV # ahahaha oh god

# put it allllll together
predsets <- c(pnames, tnames, vnames, PT, PV, TV, PTV)
length(predsets)
names(predsets) <- predsets


# YUBR --------------------------------

# make sure this works ...
test <- glmer(flr ~ pptY0Y1 + (1 | year), data = yubr, family=binomial(link="logit"))

yubr.std <- cbind(yubr[,c("flr","year")], apply(yubr[,c(pnames,tnames,vnames)], 2, function(x) (x-mean(x))/sd(x))) # 0-mean 1-sd scale
glimpse(yubr.std)

# batch-fit this mfer
yubr.mods <- lapply(predsets, function(x){
	glmer(as.formula(paste(paste('flr', paste(x,  collapse=' + '), sep = ' ~ '), '(1 | year)', sep=' + ')), data = yubr.std, family = binomial(link="logit"))
	}
)

save(yubr.mods, file="output/GLMs_YUBR.Rdata")

# incredibly, that doesn't take so long
# now, to extract AUCs ...
yubr.aics <- unlist(lapply(yubr.mods, AICc))

# pull out the best-fit to inspect ...
yubr.best <- yubr.mods[[which(yubr.aics==min(yubr.aics))]]
summary(yubr.best)
plot(yubr.best)


# format and usefully save the output ...
yubr.mods.aics <- data.frame(predictors=names(yubr.aics), RI="year", AICc=yubr.aics) %>% mutate(dAICc = AICc - min(AICc))
glimpse(yubr.mods.aics)

filter(yubr.mods.aics, dAICc <=4) # grrrr

write.table(yubr.mods.aics, "output/GLM_AICs_YUBR.txt", sep="\t", col.names=TRUE, row.names=FALSE)

# plot the predictors in raw observations ...
yubr.preds <- unlist(strsplit(names(yubr.aucs)[which(yubr.aucs==max(yubr.aucs))], split=" \\+ "))
yubr.pred.plot <- yubr %>% dplyr::select(year, flr, all_of(yubr.preds)) %>% pivot_longer(all_of(yubr.preds), names_to="Predictor", values_to="Value")

{cairo_pdf(file="output/figures/RImod_best_predictors_YUBR.pdf", width=6, height=2.5)
ggplot(yubr.pred.plot, aes(x=flr, y=Value)) + geom_jitter(alpha=0.5) + geom_boxplot() + facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value")
}
dev.off()


# for YUJA
yuja.RImods <- lapply(predsets, function(x){

rbart_vi(as.formula(paste(paste('flr', paste(c("year", x),  collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = yuja,
	group.by = yuja[,'year'],
	n.chains = 1,
	power = 2,
	base = 0.95,
	keepTrees = TRUE)
}

)

save(yuja.RImods, file="output/BART/apriori_BART_mods_YUJA.Rdata")

# incredibly, that doesn't take so long
# now, to extract AUCs ...
yuja.aucs <- unlist(lapply(yuja.RImods, function(x) auc(summary(x)$data$observed, summary(x)$data$fitted)))

# pull out the best-fit to inspect ...
yuja.best <- yuja.RImods[[which(yuja.aucs==max(yuja.aucs))]]
summary(yuja.best)
varimp(yuja.best)

# format and usefully save the output ...
yuja.RImod.aucs <- data.frame(predictors=names(yuja.aucs), RI="year", AUC=yuja.aucs)
glimpse(yuja.RImod.aucs)
write.table(yuja.RImod.aucs, "output/BART/RImod_candidates_AUCs_YUJA.txt", sep="\t", col.names=TRUE, row.names=FALSE)

# plot the predictors in raw observations ...
yuja.preds <- unlist(strsplit(names(yuja.aucs)[which(yuja.aucs==max(yuja.aucs))], split=" \\+ "))
yuja.pred.plot <- yuja %>% dplyr::select(year, flr, all_of(yuja.preds)) %>% pivot_longer(all_of(yuja.preds), names_to="Predictor", values_to="Value")

{cairo_pdf(file="output/figures/RImod_best_predictors_YUJA.pdf", width=6, height=2.5)
ggplot(yuja.pred.plot, aes(x=flr, y=Value)) + geom_jitter(alpha=0.5) + geom_boxplot() + facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value")
}
dev.off()


