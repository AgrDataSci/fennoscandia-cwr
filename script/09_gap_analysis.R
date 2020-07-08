##Load package
library("raster")
library("data.table")
library("GapAnalysis")
library("rgdal")
library("tidyverse")
library("magrittr")

sessioninfo::session_info()
# write session info
capture.output(sessioninfo::session_info(),
               file = "script/session_info/09_gap_analysis.txt")

# presence rasters
load("processing/presence_raster/presence_raster.rda")

# observed points
dt <- fread("data/passport_data.csv")

species <- fread("data/species_names.csv")

# ecoregions
eco <- readOGR(dsn = "data/gadm/ecoregions", layer = "ecoregions_eur")

# protected areas
pa <- raster("data/gadm/wdpa/wdpa_raster.tif")

# dissolve the raster
pa[pa[] > 1 ] <- 1


# ................................
# .................................
# configure observed points 
dt$source <- as.factor(ifelse(dt$source == "genesys", "G", 
                              "H"))

names(dt)[names(dt) == "lon"] <- "longitude"
names(dt)[names(dt) == "lat"] <- "latitude"
names(dt)[names(dt) == "acronym"] <- "taxon"
names(dt)[names(dt) == "source"] <- "type"

dt <- dt[, .(taxon, latitude, longitude, type)]

dt$taxon <- as.factor(dt$taxon)

dt <- as.data.frame(dt)

# get species names from the data
spp <- unique(dt$taxon)

# ..................................
# rename raster layers
sp_p <- sp_p[as.character(spp)]

for(i in seq_along(spp)) {
  names(sp_p[[i]]) <- as.character(spp)[[i]]
}

sp_p

output <- "output/gap_analysis/"
dir.create(output, showWarnings = FALSE, recursive = TRUE)


#Running all three ex situ gap analysis steps using FCSex function
FCSex_df <- FCSex(Species_list = spp,
                  Occurrence_data = dt,
                  Raster_list = sp_p,
                  Buffer_distance = 50000,
                  Ecoregions_shp = eco)

write.csv(FCSex_df, paste0(output, "fcsex.csv"), row.names = FALSE)

#Running all three in situ gap analysis steps using FCSin function
FCSin_df <- FCSin(Species_list = spp,
                  Occurrence_data = dt,
                  Raster_list = sp_p,
                  Ecoregions_shp = eco,
                  Pro_areas = pa)

write.csv(FCSin_df, paste0(output, "fcsin.csv"), row.names = FALSE)

## Combine gap analysis metrics
FCSc_mean_df <- FCSc_mean(FCSex_df = FCSex_df,FCSin_df = FCSin_df)

gap <- merge(FCSex_df, FCSin_df, by = "species")

FCSc_mean_df <- FCSc_mean_df[, c(1, which(!names(FCSc_mean_df) %in% names(gap)))]

gap <- merge(gap, FCSc_mean_df, by = "species")

names(gap)[names(gap)=="species"] <- "acronym"

write.csv(gap, paste0(output, "gap.csv"), row.names = FALSE)

# remove cols with gap class and min max values
gap <- gap[,!grepl("_class$", names(gap))]
gap <- gap[,!grepl("_min$|_max$", names(gap))]


gap %<>% 
  inner_join(species[,c("genus","species","use","acronym")], ., by = "acronym")

gap$taxa <- paste0(substr(gap$genus, start = 1, stop = 1), ". ", gap$species)

# sort values from least to high priority
gap %<>% 
  arrange(FCSc_mean)

gap$taxa <- factor(gap$taxa, levels = gap$taxa)

gap %<>% 
  select(-genus, -species, -use, -acronym)

# reshape the data.frame
gap %<>% 
  pivot_longer(-taxa, names_to = "index", values_to = "value")

gap$index <- gsub("ex$"," ex-situ", gap$index)
gap$index <- gsub("in$"," in-situ", gap$index)
gap$index <- gsub("FCSc_mean","Final Conservation Score", gap$index)

unique(gap$index)

gap$index <- factor(gap$index, levels = c("ERS ex-situ", "FCS ex-situ", "GRS ex-situ", "SRS ex-situ",
                                          "ERS in-situ", "FCS in-situ", "GRS in-situ", "SRS in-situ", 
                                          "Final Conservation Score"))

spp

ggplot() +
  geom_point(gap[!grepl("Final", gap$index), ], 
             mapping = aes(x = value, y = taxa, 
                           fill = index, color = index), 
             shape = 21) +
  annotate("rect", 
           xmin = 0, xmax = 25, 
           ymin = 0, ymax = length(spp) + 1,
           alpha = 0.3, fill = "#fb6a4a") +
  annotate("rect", 
           xmin = 25, xmax = 50, 
           ymin = 0, ymax = length(spp) + 1,
           alpha = 0.3, fill = "#fc9272") +
  annotate("rect", 
           xmin = 50, xmax = 75, 
           ymin = 0, ymax = length(spp) + 1,
           alpha = 0.3, fill = "#fcbba1") +
  annotate("rect", 
           xmin = 75, xmax = 100, 
           ymin = 0, ymax = length(spp) + 1,
           alpha = 0.3, fill = "#fee0d2") +
  annotate("text",
           x = c(12.5, 37.5, 62.5, 87.5),
           y = length(spp) + 0.7,
           label = c("High priority","Medium priority","Low priority","Sufficiently conserved"),
           colour = "black",
           size = 2.5) +
  labs(y = "", x = "Conservation score") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#f0f0f0", linetype = "dashed"))

  geom_point(gap[grepl("in-", gap$index), ], 
             mapping = aes(x = value, y = taxa,
                           fill = index, color = index), 
             shape = 25) +
  geom_point(gap[grepl("Final", gap$index), ], 
             mapping = aes(x = value, y = taxa), 
             shape = 23, col = "red", fill = "red", size = 3)


rects <- data.frame(xstart = seq(0, 75, 25), 
                    xend = seq(25, 100, 25), 
                    col = letters[1:4])


