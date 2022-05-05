# R scripts for the course "Species distribution models and fuzzy logic: Combining model predictions into threshold-free estimates of diversity and change"
# A. Marcia Barbosa (https://modtools.wordpress.com)
# Alba Estrada (https://www.researchgate.net/profile/Alba-Estrada)


# SETUP ####

# set the working directory to the folder containing this script:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# load required packages:
library(terra)
library(fuzzySim)
library(modEvA)

# import results of previous script:
models <- readRDS("../outputs/models.rds")
dat <- read.csv("../outputs/dat.csv")

names(dat)
spc_cols <- 1:3  # species codes are in these columns IN THIS DATASET (change as necessary!) JUST column #s

species <- names(dat)[spc_cols]


# CONVERT PROBABILITY TO FAVOURABILITY ####

?Fav
# this function computes favourability from either:
# - a model object of an implemented class (e.g. GLM, GAM, GBM, Random Forest, BART), or:
# - a numeric vector or a raster map of predicted probability values, and either 1) the vector of observed presence-(pseudo)absence data that produced those probabilities, i.e. the model TRAINING zeros and ones; or 2) the prevalence (proportion of presences) in the model TRAINING sample
# fav = probability without the effect of sample prevalence

# for our example we are using the WHOLE dataset for training
# at this phase you would have to specify the training prevalence versus the testing prevalence, this is what will be done when extrapolating into different regions or time periods

species
spc <- species[1]  # or any of the other 'species' available

Fav(model = models[[spc]][["GLM"]])
# or:
Fav(obs = dat[ , spc], pred = dat[ , paste0(spc, "_GLM_P")])
# or:
Fav(pred = dat[ , paste0(spc, "_GLM_P")], sample.preval = prevalence(dat[ , spc]))


# GET FAVOURABILITY FOR SEVERAL SPECIES IN A LOOP ####

for (spc in species) {
  fav_GLM <- Fav(obs = dat[ , spc], pred = dat[ , paste0(spc, "_GLM_P")])
  fav_GAM <- Fav(obs = dat[ , spc], pred = dat[ , paste0(spc, "_GAM_P")])
  # add as named columns to the table:
  dat[ , paste0(spc, "_GLM_F")] <- fav_GLM
  dat[ , paste0(spc, "_GAM_F")] <- fav_GAM
}
# note that Maxent predictions are not converted to favourability, because they are not presence probability and they do not incorporate sample prevalence
# once again expression generates a score 0 to 1

head(dat)

# note that favourability is higher than probability for low-prevalence (rare) species, and lower than probability for high-prevalence (common) species
# differences between probability and favourability are larger for data with less balanced prevalence (i.e., farther from 0.5)



prevalence(dat$Lutlut) # closest to .5
prevalence(dat$Muslut) # closer to 0
prevalence(dat$Neovis) # closer to 0


# plot sorted predictions for each model:

species
spc <- species[3]  # or any of the other 'species' available

par(mfrow = c(2, 1), mar = c(2, 2, 1, 1))  # 2 x 1 plots per window

plot(sort(dat[, paste0(spc, "_GLM_P")]), pch = ".", main = "probability")
plot(sort(dat[, paste0(spc, "_GLM_F")]), pch = ".", main = "favourability")

plot(sort(dat[, paste0(spc, "_GAM_P")]), pch = ".", main = "probability")
plot(sort(dat[, paste0(spc, "_GAM_F")]), pch = ".", main = "favourability")

# change 'spc' above and rerun the last 5 commands
# shape changes, likely based on the balance between presences and absences in the study area

# map probability and favourability:

# first, convert the table to a spatial object:
dat_map <- vect(dat, geom = c("x", "y"))

for (spc in rev(species)) {
  plot(dat_map, paste0(spc, c("_GLM_P", "_GLM_F", "_GAM_P", "_GAM_F", "_MXT_P")), col = hcl.colors(100), type = "continuous", range = c(0, 1), cex = 0.5, nr = 3, nc = 2)
}  # press the 'back' arrow in the plotting pane to see maps for other species

# notice that Maxent predictions resemble favourability more than probability (especially when prevalence is far from 0.5)
# Lutlut is closest to .5, smallest MaxEnt difference
# Muslut large delta between 

# save the current version of the table, now with probability and favourability:

write.csv(dat, "../outputs/dat.csv", row.names = FALSE)


# if you have presence probability models of your own, convert them to favourability and notice the changes

# remember that favourability requires TRAINING observations or TRAINING prevalence

# if you have any difficulties, during the course days we're here to help!
