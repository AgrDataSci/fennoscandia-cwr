# clip bioclim rasters within 
# the target research area
# source data is available in an external hard drive or 
# http://worldclim.org/
library("raster")
library("rgdal")

source("script/helper01_functions.R")

gadm <- readOGR("data/gadm/europe", "europe")

ext <- floor(extent(bbox(gadm)))

r <- raster(ext, res = 0.04166666)

gadm <- rasterize(gadm, r, field = 1)

# ...........................................
# ...........................................
# worldclim 1.0
ph <- "/Users/kauedesousa/Large-raster-data/wc2.1_2.5m_bio"
# read rasters from an external hard drive
r <- stack(list.files(ph, 
                      pattern = ".tif",
                      full.names = TRUE))
# crop it
r <- crop(r, ext)

# mask it
r <- mask(r, gadm)

# check the first layer
plot(r[[1]])

# rename rasters
names(r)
names(r) <- gsub("wc2.1_2.5m_","",names(r))
names(r)
# save it on the project folder
output <- "data/bioclim/"

dir.create(output,
           showWarnings = FALSE,
           recursive = TRUE)

r <- stack(r)

names_r <- paste0(output, names(r))

list.files(output)
file.remove(list.files("data/bioclim/",
                       pattern = ".bil|.hdr|.tif", 
                       full.names = TRUE))

list.files(output)


writeRaster(r, 
            filename = names_r, 
            format = "GTiff",
            bylayer = TRUE,
            overwrite = TRUE)

# ..................................
# ..................................
# evapotranspiration
ph <- "/Users/kauedesousa/Large-raster-data/et0_yr/"
r <- raster(paste0(ph, "et0_yr.tif"))

r <- crop(r, ext)

# reduce resolution
r <- aggregate(r, fact = 5)

r <- mask(r, gadm)

plot(r)

writeRaster(r,
            filename = paste0(output, "eto.tif"),
            format = "GTiff",
            overwrite = TRUE)

# ..................................
# ..................................
# solar radiation
ph <- "/Users/kauedesousa/Large-raster-data/wc2.1_2.5m_srad"
ph2 <- "/Users/kauedesousa/Large-raster-data/wc2.1_2.5m_vapr"
ph3 <- "/Users/kauedesousa/Large-raster-data/wc2.1_2.5m_wind"
# read rasters from an external hard drive
r <- stack(list.files(ph, 
                      pattern = ".tif",
                      full.names = TRUE))

r2 <- stack(list.files(ph2, 
                       pattern = ".tif",
                       full.names = TRUE))

r3 <- stack(list.files(ph3, 
                       pattern = ".tif",
                       full.names = TRUE))


r <- stack(r, r2, r3)

# crop it
r <- crop(r, ext)

# mask it
r <- mask(r, gadm)

# check the first layer
plot(r[[1]])

# rename rasters
names(r)
names(r) <- gsub("wc2.1_2.5m_","",names(r))
names(r)

r <- stack(r)

names_r <- paste0(output, names(r))

list.files(output)

writeRaster(r, 
            filename = names_r, 
            format = "GTiff",
            bylayer = TRUE,
            overwrite = TRUE)
