##Load package
library("raster")
library("data.table")
library("GapAnalysis")
library("rgdal")

sessioninfo::session_info()
# write session info
capture.output(sessioninfo::session_info(),
               file = "script/session_info/09_gap_analysis.txt")

# obtaining occurrences from GENESYS
# data(CucurbitaData)

# presence rasters
load("processing/presence_raster/presence_raster.rda")

# observed points
dt <- fread("data/passport_data.csv")

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

##Obtaining raster_list
#data(CucurbitaRasters)
#CucurbitaRasters <- raster::unstack(CucurbitaRasters)



##Obtaining protected areas raster
data(ProtectedAreas)

##Obtaining ecoregions shapefile
data(ecoregions)

#Running all three ex situ gap analysis steps using FCSex function
FCSex_df <- FCSex(Species_list=speciesList,
                  Occurrence_data=CucurbitaData,
                  Raster_list=CucurbitaRasters,
                  Buffer_distance=50000,
                  Ecoregions_shp=ecoregions
)

#Running all three in situ gap analysis steps using FCSin function
FCSin_df <- FCSin(Species_list=speciesList,
                  Occurrence_data=CucurbitaData,
                  Raster_list=CucurbitaRasters,
                  Ecoregions_shp=ecoregions,
                  Pro_areas=ProtectedAreas)

## Combine gap analysis metrics
FCSc_mean_df <- FCSc_mean(FCSex_df = FCSex_df,FCSin_df = FCSin_df)

##Running Conservation indicator across taxa
indicator_df  <- indicator(FCSc_mean_df)

