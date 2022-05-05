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
library(gam)
library(maxnet)

# import results of previous scripts:
dat <- read.csv("../outputs/dat.csv")
models <- readRDS("../outputs/models.rds")

# define occurrence and favourability columns for the species of interest:
names(dat)
spc_cols <- 1:3  # species are in these columns IN THIS DATASET (change as necessary!)
fav_cols <- grep("_GLM_F", names(dat))  # using GLM predictions for the examples
names(dat)[spc_cols]  # check if OK
names(dat)[fav_cols]  # check if OK


# DISTRIBUTIONAL SIMILARITY between different species ####
# (binary similarity if you use the presence-absence data; fuzzy similarity if you use the favourability values)

# map presence and favourability for 2 related species - e.g. Lutra lutra and Neovisn vison, which are potential competitors for habitat and food:

dat_map <- vect(dat, geom = c("x", "y"))
plot(dat_map, c("Lutlut", "Neovis", "Lutlut_GLM_F", "Neovis_GLM_F"), col = hcl.colors(100), type = "continuous", range = c(0, 1), cex = 0.4)

?fuzSim
with(dat, fuzSim(Lutlut, Neovis, method = "Jaccard"))
with(dat, fuzSim(Lutlut_GLM_F, Neovis_GLM_F, method = "Jaccard"))
with(dat, fuzSim(Lutlut, Neovis, method = "Baroni"))
with(dat, fuzSim(Lutlut_GLM_F, Neovis_GLM_F, method = "Baroni"))
# notice that distributions tend to be more similar with fuzzy than with crisp (observed presence/absence) values, as crisp values overlap less

# calculate binary and fuzzy similarity among SEVERAL species:
bin_sim_mat <- simMat(dat[ , spc_cols], method = "Jaccard")
fuz_sim_mat <- simMat(dat[ , fav_cols], method = "Jaccard")
bin_sim_mat
fuz_sim_mat

# represent similarity matrix with colours:

par(mfrow = c(1, 1), mar = c(6, 6, 2, 1))
image(x = 1:ncol(bin_sim_mat), y = 1:nrow(bin_sim_mat), z = bin_sim_mat, col = rev(heat.colors(10)), xlab = "", ylab = "", axes = FALSE, main = "Binary similarity")
axis(side = 1, at = 1:ncol(bin_sim_mat), tick = FALSE, labels = colnames(bin_sim_mat), las = 2, cex.axis = 0.8)
axis(side = 2, at = 1:nrow(bin_sim_mat), tick = FALSE, labels = rownames(bin_sim_mat), las = 2, cex.axis = 0.8)

image(x = 1:ncol(fuz_sim_mat), y = 1:nrow(fuz_sim_mat), z = fuz_sim_mat, col = rev(heat.colors(10)), xlab = "", ylab = "", axes = FALSE, main = "Fuzzy similarity")
axis(side = 1, at = 1:ncol(fuz_sim_mat), tick = FALSE, labels = colnames(fuz_sim_mat), las = 2, cex.axis = 0.8)
axis(side = 2, at = 1:nrow(fuz_sim_mat), tick = FALSE, labels = rownames(fuz_sim_mat), las = 2, cex.axis = 0.8)

# plot a dendrogram based on the distributional similarities:
plot(hclust(as.dist(1 - bin_sim_mat), method = "average"), main = "Binary dendrogram")
plot(hclust(as.dist(1 - fuz_sim_mat), method = "average"), main = "Fuzzy dendrogram")

# see ?simMat for additional options


# FUTURE FAVOURABILITY for different species ####

# first, we need to get the predictions of our models for a future scenario:

names(models)
names(models[[1]])

species <- names(models)
species

mods <- names(models[[1]])
mods

future <- rast(list.files("../data/variables/future", pattern = "\\.tif$", full.names = TRUE))
future
plot(future, maxnl = 20)

# note this is just one of many future climate scenarios, used here as an example!
# you should normally analyse several scenarios available for your set of variables and study region


# extract future variables from the rasters to a table:
head(dat)
dat_fut <- extract(future, dat[ , c("x", "y")], cells = TRUE, xy = TRUE)
head(dat_fut)


# get future probability and future favourability for each species and each model:

for (spc in species) {
  mods <- models[[spc]][c("GLM", "GAM")]  # Maxent left out because it does not produce probability
  
  for (m in names(mods)) {
    # get predicted probability for the future:
    prob <- predict(mods[[m]], newdata = dat_fut, type = "response")
  
    # convert to favourability for the future (using modelled prevalence):
    fav <- Fav(pred = prob, sample.preval = prevalence(model = mods[[m]]))
    
    # add results as new named columns to the table:
    dat_fut[ , paste(spc, m, "P_fut", sep = "_")] <- prob
    dat_fut[ , paste(spc, m, "F_fut", sep = "_")] <- fav
  }
}

head(dat_fut)
names(dat_fut)


# map the results:

dat_fut_map <- vect(dat_fut, geom = c("x", "y"))

fav_cols <- grep("_F$", names(dat_map))
fav_cols_fut <- grep("_F_fut", names(dat_fut_map)) # grep out column #s

names(dat)

# present
plot(dat_map, fav_cols, cex = 0.4, col = hcl.colors(100), type = "continuous", range = c(0, 1), nr = 3, nc = 2)
# future
plot(dat_fut_map, fav_cols_fut, cex = 0.4, col = hcl.colors(100), type = "continuous", range = c(0, 1), nr = 3, nc = 2)


spc <- "Neovis"
mod <- "GAM"

par(mfrow = c(2, 1))
plot(dat_map, paste(spc, mod, "F", sep = "_"), cex = 0.8, col = hcl.colors(100), type = "continuous", range = c(0, 1), nc = 1)
plot(dat_fut_map, paste(spc, mod, "F_fut", sep = "_"), cex = 0.8, col = hcl.colors(100), type = "continuous", range = c(0, 1), nc = 1)

# repeat the last two commands after changing 'spc' and/or 'mod' above
# to see present and future favourability for each species and model
# huge contraction of favorability in the future for Neovis based on GAM


# (FUZZY) RANGE CHANGE BETWEEN TIME PERIODS ####

?fuzzyRangeChange

names(dat)
names(dat_fut)

par(mfrow = c(1, 1), mar = c(5, 4, 2, 1))
fuzzyRangeChange(dat[ , "Muslut_GLM_F"], dat_fut[ , "Muslut_GLM_F_fut"], las = 2, main = "Muslut")
fuzzyRangeChange(dat[ , "Neovis_GLM_F"], dat_fut[ , "Neovis_GLM_F_fut"], las = 2, main = "Neovis")
fuzzyRangeChange(dat[ , "Lutlut_GLM_F"], dat_fut[ , "Lutlut_GLM_F_fut"], las = 2, main = "Lutlut")


# try fuzzyRangeChange with your own probability predictions for different time periods
# don't forget to convert to favourability first, using training (present) prevalence!


# SHARED FAVOURABILITY FOR TWO INTERACTING SPECIES ####

?sharedFav  # read the help file for detailed info and references!
# circles and a continuous line representing favourability for the stronger species, and squares and a dashed line representing favourability for the weaker species. The height of the bars at the bottom represents the proportional sample size of each bin

par(mfrow = c(1, 1))
names(dat)
with(dat, sharedFav(Neovis_GLM_F, Muslut_GLM_F, las = 2))  # this may generate an error, see "Note" in ?sharedFav
with(dat, sharedFav(Neovis_GLM_F, Muslut_GLM_F, bin_interval = "quantiles", las = 2))

with(dat, sharedFav(Neovis_GLM_F, Lutlut_GLM_F, las = 2))
with(dat, sharedFav(Neovis_GLM_F, Lutlut_GLM_F, bin_interval = "quantiles", las = 2))


# BIOTIC THREAT FOR TWO INTERACTING SPECIES ####

?bioThreat  # read the help file for detailed info and references!

dat$threat_NvMl <- with(dat, bioThreat(Neovis_GLM_F, Muslut_GLM_F))
dat$threat_NvLl <- with(dat, bioThreat(Neovis_GLM_F, Lutlut_GLM_F))
dat$threat_LlNv <- with(dat, bioThreat(Lutlut_GLM_F, Neovis_GLM_F))

dat_map <- vect(dat, geom = c("x", "y"))

# grey is areas of complete non-overal (no threat)
par(mfrow = c(2, 2))
plot(dat_map, "Neovis_GLM_F", cex = 0.4, col = hcl.colors(100), type = "continuous", range = c(0, 1), main = "Neovis favourability")
plot(dat_map, "Muslut_GLM_F", cex = 0.4, col = hcl.colors(100), type = "continuous", range = c(0, 1), main = "Muslut favourability")
plot(dat_map, "threat_NvMl", col = c("lightgrey", "green", "yellow", "orange", "red"), cex = 0.4, main = "Biotic threat\nNeovis over Muslut")

par(mfrow = c(2, 2))
plot(dat_map, "Neovis_GLM_F", cex = 0.4, col = hcl.colors(100), type = "continuous", range = c(0, 1), main = "Neovis favourability")
plot(dat_map, "Lutlut_GLM_F", cex = 0.4, col = hcl.colors(100), type = "continuous", range = c(0, 1), main = "Lutlut favourability")
plot(dat_map, "threat_NvLl", col = c("lightgrey", "green", "yellow", "orange", "red"), cex = 0.4, main = "Biotic threat\nNeovis over Lutlut")
plot(dat_map, "threat_LlNv", col = c("lightgrey", "green", "yellow", "orange", "red"), cex = 0.4, main = "Biotic threat\nLutlut over Neovis")
