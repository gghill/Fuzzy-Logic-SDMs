# R scripts for the course "Species distribution models and fuzzy logic: Combining model predictions into threshold-free estimates of diversity and change"
# A. Marcia Barbosa (https://modtools.wordpress.com)
# Alba Estrada (https://www.researchgate.net/profile/Alba-Estrada)


# SETUP ####

# the following command sets the working directory to the folder containing this script (similar to RStudio menu "Session - Set Working Directory - To Source File Location"):
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# the next command creates a folder for the outputs (if it doesn't already exist):
if (!file.exists("../outputs")) dir.create("../outputs")  # '../' goes up one level from the current working directory, so this creates the 'outputs' folder just outside the 'scripts' folder

# load required packages:

# first, run script 0 to install the latest versions of these packages!
# as we will need functions not yet available in their CRAN versions

library(terra)
library(fuzzySim)
library(gam)
library(maxnet)


# IMPORT SPECIES OCCURRENCE DATA AND RASTER VARIABLES ####

occurrences <- read.csv("../data/species_occurrences/species_occurences.csv")

head(occurrences)

unique(occurrences$species) # 3 species, European mink, American mink, Eurasian otter 


variables <- rast(list.files("../data/variables", pattern = "\\.tif$", full.names = TRUE))

variables # 20 bio variables (worldclim, https://www.worldclim.org/data/bioclim.html)

plot(variables, maxnl = 20) # maxnl = max # layers


# GRID OCCURRENCE DATA TO THE VARIABLE PIXELS ####

# to get environmental values at pixels with and without presence records for each species

# note that you normally need to clean the species occurrences and define an appropriate extent and spatial resolution for your models!
# we won't go into that here, because it is outside the scope of this course
# you can check e.g. ?scrubr::scrubr or ?CoordinateCleaner::CoordinateCleaner for cleaning occurrence data before doing any real work on them

# first, create a 3+3 letter code for each species, to use as column names below:
# spCodes is designed specifically for this purpose - obtain unique abbreviations
occurrences$spcode <- spCodes(occurrences$species, 3, 3)
unique(occurrences[ , c("species", "spcode")])

?gridRecords

dat <- gridRecords(rst = variables, pres.coords = occurrences[ , c("lon", "lat")], species = occurrences$spcode)
# note you can use also the 'abs.coords' argument if you have absence or background point coordinates

head(dat)
head(dat[ , 1:7]) # don't worry about row names/indices
nrow(dat)

plot(variables[[1]]) # presence absence over enviro variable
points(subset(dat, Muslut == 1, select = c("x", "y")), pch = 20, cex = 0.5, col = "blue")
points(subset(dat, Muslut == 0, select = c("x", "y")), pch = 20, cex = 0.1, col = "red")


# COMPUTE SOME DISTRIBUTION MODELS FOR ONE SPECIES ####

# DOUBLE CHECK WHEN ADAPTING #
names(dat)
spc_cols <- 1:3  # species codes are in these columns IN THIS DATASET (change as necessary!)
var_cols <- 7:25  # variables are in these columns IN THIS DATASET (change as necessary!)
names(dat)[spc_cols]  # check if OK!
names(dat)[var_cols]  # check if OK!


## choose an example species to model:

species <- names(dat)[spc_cols]
species
spc <- species[1]


## remove redundant variables for this species:

?corSelect
vars_sel <- corSelect(data = dat, sp.cols = spc, var.cols = var_cols, cor.thresh = 0.7) # what if results not identical based on variable selection params?
vars_sel # VIF = variance inflation factor
vars_sel <- vars_sel$selected.vars
vars_sel

# note this is just one possible way of selecting variables for modelling!
# you should normally use additional selection criteria, depending on the modelling method
# we won't go deeper here, because it is outside the scope of this course


## make a GLM for this species ####

?glm

form_GLM <- as.formula(paste(spc, "~", paste(vars_sel, collapse = "+")))
form_GLM # creating the correct formula combining species and variables

mod_GLM <- glm(formula = form_GLM, family = binomial, data = dat)

summary(mod_GLM) # all high sig


## make a GAM for this species ####

?gam

# useful for forming future models
form_GAM <- as.formula(paste(spc, "~", paste0("s(", vars_sel, ")", collapse = "+")))  # GAM with smoothing splines ('s')
form_GAM

mod_GAM <- gam(formula = form_GAM, family = binomial, data = dat)

summary(mod_GAM) # some sig


## make a Maxent model for this species ####

?maxnet # make maxent model using glmnet for fitting

mod_MXT <- maxnet(p = dat[ , spc], data = dat[ , vars_sel], f = maxnet.formula(dat[ , spc], dat[ , vars_sel], classes = "lq"))  # using linear ('l') and quadratic ('q') features

mod_MXT


# GET MODEL PREDICTIONS FROM THE TABLE OF VARIABLES ####

# (variables must have exact same names as in the models):

# moving from model to prediction (quantitative, continuous)
?predict.glm
predict(mod_GLM, newdata = dat, type = "response")  # "response" transforms the linear predictor with the appropriate link function, yielding results as probability values

?predict.Gam
predict(mod_GAM, dat, type = "response")

?predict.maxnet
predict(mod_MXT, newdata = dat, type = "cloglog")


# COMPUTE THE SAME MODELS FOR SEVERAL SPECIES IN A LOOP ####

# create an empty list to receive the results for each species:
models <- vector("list", length(species))
names(models) <- species # ordering important

for (spc in species) {
  vars_sel <- corSelect(data = dat, sp.cols = spc, var.cols = var_cols, cor.thresh = 0.7)$selected.vars
  
  form_GLM <- as.formula(paste(spc, "~", paste(vars_sel, collapse = "+")))
  mod_GLM <- glm(formula = form_GLM, family = binomial, data = dat)
  
  form_GAM <- as.formula(paste(spc, "~", paste0("s(", vars_sel, ")", collapse = "+")))  # GAM with smoothing splines ('s')
  mod_GAM <- gam(formula = form_GAM, family = binomial, data = dat)
  
  mod_MXT <- maxnet(p = dat[ , spc], data = dat[ , vars_sel], f = maxnet.formula(dat[ , spc], dat[ , vars_sel], classes = "lq"))  # using linear ('l') and quadratic ('q') features

  models[[spc]][["GLM"]] <- mod_GLM
  models[[spc]][["GAM"]] <- mod_GAM
  models[[spc]][["MXT"]] <- mod_MXT
}

models # large list of models which are each the same as above but all at once
names(models)
names(models[[1]])
models[["Muslut"]][["GLM"]]

# save results to disk:
saveRDS(models, "../outputs/models.rds")


# GET MODEL PREDICTIONS FOR SEVERAL SPECIES IN A LOOP ####

# and add each prediction as a new named column in the table

for (spc in species) {
  pred_GLM <- predict(models[[spc]][["GLM"]], newdata = dat, type = "response")
  pred_GAM <- predict(models[[spc]][["GAM"]], newdata = dat, type = "response")
  pred_MXT <- as.vector(predict(models[[spc]][["MXT"]], newdata = dat, type = "cloglog"))
  
  dat[ , paste0(spc, "_GLM_P")] <- pred_GLM
  dat[ , paste0(spc, "_GAM_P")] <- pred_GAM
  dat[ , paste0(spc, "_MXT_P")] <- pred_MXT
}

# see the table with the newly added columns:
head(dat) # now we have the predictions which are using the model algorithm to generate predicted values that have coordinates in space for each relevant bio variable (MAPPABLE)


# MAP MODEL PREDICTIONS ####

# first we need to convert the table to a spatial object, using the coordinate columns as geometry:
dat_map <- vect(dat, geom = c("x", "y"))

par(mfrow = c(3, 3))

for (spc in species) {
  plot(dat_map, paste0(spc, "_GLM_P"), col = hcl.colors(100), type = "continuous", range = c(0, 1), cex = 0.5, axes = FALSE)
  plot(dat_map, paste0(spc, "_GAM_P"), col = hcl.colors(100), type = "continuous", range = c(0, 1), cex = 0.5, axes = FALSE)
  plot(dat_map, paste0(spc, "_MXT_P"), col = hcl.colors(100), type = "continuous", range = c(0, 1), cex = 0.5, axes = FALSE)
}


# observe the differences between GLM and GAM probability vs. Maxent suitability predictions
# for species with lower and higher prevalence (proportion of presences)
# similar patterns, but difference between Maxent and others
names(dat)
# presence absence, 0s and 1s, GLM and GAM are both presence absence, GLM simpler relying on linear relationships, GAM can fit curves
# prediction columns are predicting occurence from 0 (impossible) to 1 (certain)
table(dat$Muslut)
table(dat$Neovis)
table(dat$Lutlut)


# save table with predictions to disk:
write.csv(dat, "../outputs/dat.csv", row.names = FALSE)


# CHECK CORRELATIONS AMONG PREDICTIONS ####

species
spc <- species[1]  # or any of the other 'species' available

names(dat)
# only the model columns are relevant to the actual predictions at this point
pred_cols <- paste0(spc, c("_GLM_P", "_GAM_P" , "_MXT_P"))
pred_cols

cor(dat[ , pred_cols], use = "pairwise.complete.obs")
# comparing the different models via scatterplot where each axis is a prediction by one model
# analogous to the prediction map
pairs(dat[ , pred_cols], pch = 20, cex = 0.5)

# change 'spc' above and rerun the last 5 commands


# note that you normally need to thoroughly evaluate and/or cross-validate your models before moving on!
# we won't do that here, because it is outside the scope of this course
# but you can check e.g. the 'modEvA' package for evaluating models before basing any real work on them
