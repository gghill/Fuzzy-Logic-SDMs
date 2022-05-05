# R scripts for the course "Species distribution models and fuzzy logic: Combining model predictions into threshold-free estimates of diversity and change"
# A. Marcia Barbosa (https://modtools.wordpress.com)
# Alba Estrada (https://www.researchgate.net/profile/Alba-Estrada)


# SETUP ####

# set the working directory to the folder containing this script:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# load required packages:
library(fuzzySim)
library(terra)

# import results of previous scripts:
dat <- read.csv("../outputs/dat.csv")

# define occurrence and favourability columns for all species of interest:
names(dat)
spc_cols <- 1:3  # species are in these columns IN THIS DATASET (change as necessary!)
fav_cols <- grep("_GLM_F$", names(dat))  # using GLM predictions for the examples
names(dat)[spc_cols]  # check if OK
names(dat)[fav_cols]  # check if OK


# FUZZY SET OPERATIONS BETWEEN SPECIES DISTRIBUTIONS ####

?fuzzyOverlay


# (fuzzy) intersection ####
# areas of (favourability for) co-occurrence:
# using pure occurence columns for crips, favorability columns for fuzzy

dat$intersection <- fuzzyOverlay(data = dat, overlay.cols = spc_cols, op = "intersection")
dat$intersection_fuzzy <- fuzzyOverlay(data = dat, overlay.cols = fav_cols, op = "intersection")
# you can also intersect e.g. favourability for the same species in different time periods


# (fuzzy) union ####
# areas with (or favourable for) either species

dat$union <- fuzzyOverlay(data = dat, overlay.cols = spc_cols, op = "union")
dat$union_fuzzy <- fuzzyOverlay(data = dat, overlay.cols = fav_cols, op = "union")


# (fuzzy) segregation ####
# areas with (or favourable for) one species and not another
# pairwise comparisons, 2 rounds to accomodate 3 species.

dat$LLnotNV <- fuzzyOverlay(data = dat, overlay.cols = c("Lutlut", "Neovis"), op = "AnotB")
dat$LLnotNV_fuzzy <- fuzzyOverlay(data = dat, overlay.cols = c("Lutlut_GLM_F", "Neovis_GLM_F"), op = "AnotB")

dat$NVnotLL <- fuzzyOverlay(data = dat, overlay.cols = c("Neovis", "Lutlut"), op = "AnotB")
dat$NVnotLL_fuzzy <- fuzzyOverlay(data = dat, overlay.cols = c("Neovis_GLM_F", "Lutlut_GLM_F"), op = "AnotB")


# map the results:

dat_map <- vect(dat, geom = c("x", "y"))

# intersection
plot(dat_map, c("intersection", "intersection_fuzzy", "union", "union_fuzzy"), cex = 0.5, type = "continuous", col = hcl.colors(100), range = c(0, 1))

# segregation
plot(dat_map, c("LLnotNV", "LLnotNV_fuzzy", "NVnotLL", "NVnotLL_fuzzy"), cex = 0.5, type = "continuous", col = hcl.colors(100), range = c(0, 1))


# CRISP AND FUZZY AREA OF OCCUPANCY of a species ####

# area of occupancy = sum of occupied cells in study area (cell values 0/1)
# fuzzy / potential area of occupancy = sum of favourability values in study area (cell values 0-1)

plot(dat_map, c("Muslut", "Muslut_GLM_F"), cex = 0.5, col = hcl.colors(100), type = "continuous", range = c(0, 1), nc = 1)
sum(dat[ , "Muslut"], na.rm = TRUE)
sum(dat[ , "Muslut_GLM_F"], na.rm = TRUE)

plot(dat_map, c("Neovis", "Neovis_GLM_F"), cex = 0.5, col = hcl.colors(100), type = "continuous", range = c(0, 1), nc = 1)
sum(dat[ , "Neovis"], na.rm = TRUE)
sum(dat[ , "Neovis_GLM_F"], na.rm = TRUE)

plot(dat_map, c("Lutlut", "Lutlut_GLM_F"), cex = 0.5, col = hcl.colors(100), type = "continuous", range = c(0, 1), nc = 1)
sum(dat[ , "Lutlut"], na.rm = TRUE)
sum(dat[ , "Lutlut_GLM_F"], na.rm = TRUE)
# Lutlut most widely distributed

# (CRISP AND) FUZZY ENTROPY ####

?entropy
# entropy (simply) is distance from 0 or 1, values close to .5 (uncertain whether present/absent, favorable/unfavorable)

# fuzzy entropy (requires a fuzzy system, such as favourability):
par(mfrow = c(1, 1), mar = c(5, 4, 2, 1))
entropy(dat, spc_cols)  # fuzzy entropy is zero because all cells are classified as either presence or absence, i.e. values are not fuzzy
entropy(dat, fav_cols)  # entropy is higher because values are fuzzy, rather than clear zeros and ones
# the fuzzier the values (i.e., the farther from 0 or 1), the higher the entropy
# so in our example, lutlut has the highest fuzzy area of occupancy, but also the highest entropy


# CRISP AND FUZZY SPECIES RICHNESS (potential biodiversity) ####

# species richness = sum of species in each cell
# fuzzy / potential species richness = sum of favourability values in each cell

names(dat)
dat$SR <- rowSums(dat[ , spc_cols]) # richness = sums
dat$SR_fuzzy <- rowSums(dat[ , fav_cols])
head(dat)


# map the results:

dat_map <- vect(dat, geom = c("x", "y"))

plot(dat_map, c("SR", "SR_fuzzy"), cex = 0.5, col = hcl.colors(100), type = "continuous", range = range(dat_map$SR, na.rm = TRUE), nc = 1)


# CRISP AND FUZZY RARITY ####

?rarity
# Rarity is like a (potential) richness index in which rarer species have higher weight.
# fuzzy rarity: species with a smaller sum of favorabilities get weighted more

# crisp and fuzzy rarity for each individual species:
ra <- sapply(dat[ , spc_cols], rarity)
ra
ra_fuzzy <- sapply(dat[ , fav_cols], rarity)
ra_fuzzy

# visualize as bar plots:
par(mar = c(4, 7, 2, 1), mfrow = c(2, 1))
barplot(sort(ra), horiz = TRUE, las = 2, main = "Rarity")
barplot(sort(ra_fuzzy), horiz = TRUE, las = 2, main = "Fuzzy rarity")
# trend is the same, but Lutlut jumps up in fuzzy rarity.

# combined crisp and fuzzy rarity for all species across the study area:
dat$rarity <- rarity(dat, sp.cols = spc_cols)
dat$rarity_fuzzy <- rarity(dat, sp.cols = fav_cols)
head(dat)

# rarity range and histograms:
range(dat$rarity, na.rm = TRUE)
range(dat$rarity_fuzzy, na.rm = TRUE)
hist(dat$rarity)
hist(dat$rarity_fuzzy)

# map rarity and fuzzy rarity:
dat_map <- vect(dat, geom = c("x", "y"))
plot(dat_map, c("rarity", "rarity_fuzzy"), cex = 0.5, col = hcl.colors(100), type = "continuous", nc = 1)
# you can try plotting both maps with the same colour scale by adding range = range(dat[ , c("rarity", "rarity_fuzzy")], na.rm = TRUE)


# CRISP AND FUZZY VULNERABILITY ####

?vulnerability

# set the conservation status of each species:
Muslut_cat <- 16  # CR, as per https://www.iucnredlist.org
Lutlut_cat <- 2  # NT, as per https://www.iucnredlist.org
# we don't add Neovis because it's an invasive species here

names(dat)
dat$vulnerability <- vulnerability(dat, c("Muslut", "Lutlut"), categories = c(Muslut_cat, Lutlut_cat))
dat$vulnerability_fuzzy <- vulnerability(dat, c("Muslut_GLM_F", "Lutlut_GLM_F"), categories = c(Muslut_cat, Lutlut_cat))
head(dat)

# vulnerability range and histograms:
range(dat$vulnerability, na.rm = TRUE)
range(dat$vulnerability_fuzzy, na.rm = TRUE)
hist(dat$vulnerability)
hist(dat$vulnerability_fuzzy)

# map vulnerability and fuzzy vulnerability:
dat_map <- vect(dat, geom = c("x", "y"))
names(dat_map)
plot(dat_map, c("vulnerability", "vulnerability_fuzzy"), col = hcl.colors(100), type = "continuous", range = range(dat[ , c("vulnerability", "vulnerability_fuzzy")], na.rm = TRUE), nc = 1)

# so the rarer species are concentrated in the North and center of iberian peninsula while the more prevalent Lutlut is widespread. Rarity and vulnerability concentrate in the range of the less prevalent 

# save the current version of the table, now with the conservation metrics:

write.csv(dat, "../outputs/dat.csv", row.names = FALSE)
