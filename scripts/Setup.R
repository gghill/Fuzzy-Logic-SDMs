## SETUP ##
require('terra')
library(spocc)
library(rgbif)
require('leaflet')
require('raster')
require('modEvA')
require(gam)
require(fuzzySim)
require(maxnet)

my_species <- c('Caligus elongatus', 'Lepeophtheirus salmonis')

spocc_data <- occ(my_species, from = c("gbif", "obis"), has_coords = TRUE, limit = 5000)

names(spocc_data)
spocc_data # 5000 limit seems to have captured everything
# spocc_data <- spocc_data[c("gbif","obis")] # slice out only the used databases

# save this object to disk:
saveRDS(spocc_data, paste0("./outputs/spocc_", my_species[1], "_raw.rds"))  # we name it "raw" because the data haven't been cleaned yet

spocc_data <- readRDS(paste0("./outputs/spocc_", my_species[1], "_raw.rds"))

# combine the 'spocc_data' list into one data frame:
?occ2df
spocc_df <- occ2df(spocc_data)
spocc_df

table(spocc_df$prov)  # providers with records for my_species

str(spocc_df)  # longitude and latitude are character...
# convert them to numeric for mapping:
spocc_df$longitude <- as.numeric(spocc_df$longitude)
spocc_df$latitude <- as.numeric(spocc_df$latitude)

spocc_df <- spocc_df[spocc_df$latitude>0,]
write.table(spocc_df, file = 'lice_df.txt', sep = "\t", row.names = FALSE, col.names = TRUE)
spocc_df <- read.delim(file = 'lice_df_clean_names.txt', header = TRUE)
unique(spocc_df$name)
C_elongatus <- c("Caligus elongatus Nordmann, 1832", # accommodating alt. names in database
                 "BOLD:AAE8404",
                 "BOLD:AAE8403",
                 "Caligus elongatus")

L_salmonis <- c("Lepeophtheirus salmonis (Krøyer, 1838)", # accommodating alt. names in database
                "Lepeophtheirus salmonis oncorhynchi Skern-Mauritzen, Torrissen & Glover, 2014",
                "BOLD:AAA1252",
                "BOLD:AAA1251",
                "Lepeophtheirus salmonis salmonis",
                "Caligus salmonis Krøyer, 1837")


# TEST PLOT ---------------------------------------------------------------
# legend function to work with circle markers
addLegendCustom <- function(map, colors, labels, sizes = 20, opacity = 0.5){
  colorAdditions <- paste0(colors, "; border-radius: 50%; width:", sizes, "px; height:", sizes, "px")
  labelAdditions <- paste0("<div style='display: inline-block;height: ", 
                           sizes, "px;margin-top: 4px;line-height: ", sizes, "px;'>", 
                           labels, "</div>")
  
  return(addLegend(map, colors = colorAdditions, 
                   labels = labelAdditions, opacity = opacity))
}
leaflet() %>%
  addTiles() %>%
  # add global raster + legend
  # add occurrence points (circles), defining which TABLE and COLUMNS contain their longitude and latitude coordinates:
  addCircles(data = subset(spocc_df, (name) %in% (C_elongatus)), lng = ~ longitude, lat = ~ latitude, color = 'blue') %>%
  addCircles(data = subset(spocc_df, (name) %in% (L_salmonis)), lng = ~ longitude, lat = ~ latitude, color = 'red') %>%
  addLegendCustom(colors = c('blue', 'red'),
                  labels = my_species)

