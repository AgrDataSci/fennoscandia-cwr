##Load package
library("raster")
library("data.table")
library("GapAnalysis")
library("rgdal")

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

gap <- merge(species[,c("genus","species","use","acronym")], gap, by = "acronym")

write.csv(gap, paste0(output, "gap.csv"), row.names = FALSE)

