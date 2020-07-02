##Load package
library(raster)
library(GapAnalysis)

##Obtaining occurrences from example
data(CucurbitaData)

##Obtaining species names from the data
speciesList <- unique(CucurbitaData$taxon)

##Obtaining raster_list
data(CucurbitaRasters)
CucurbitaRasters <- raster::unstack(CucurbitaRasters)

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

## Generate summary HTML file with all result
GetDatasets()
summaryHTML_file <- SummaryHTML(Species_list=speciesList,
                                Occurrence_data = CucurbitaData,
                                Raster_list=CucurbitaRasters,
                                Buffer_distance=50000,
                                Ecoregions_shp=ecoregions,
                                Pro_areas=ProtectedAreas,
                                Output_Folder=".",
                                writeRasters=F)

