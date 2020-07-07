library("sf")
library("rgdal")
library("rgeos")
library("raster")
library("maptools")
library("tidyverse")
library("magrittr")


# read file with countries iso codes
iso <- read_csv("data/country_iso.csv")

iso %<>% 
  filter(region == "Europe") %>% 
  select(alpha3) %>% 
  t()

# .............................
# .............................
# read world polygon and keep European countries
# data is kept in an external hard drive
eur <- readOGR(dsn = "/Volumes/BioversityInt/shapefiles/world_borders/", layer = "world_borders_adm0")

eur <- eur[eur$ADM0_A3 %in% iso, ]

# pick focal columns 
eur@data <- eur@data[,c("ADMIN","ADM0_A3", "GEOUNIT")]

# set a new bbox to remove countries areas outside Europe
e <- as(raster::extent(-25, 60, 34, 71), "SpatialPolygons")
proj4string(e) <- proj4string(eur)

# intersect
eur <- raster::intersect(eur, e)

plot(eur)

writeOGR(eur, 
         dsn = "data/gadm/europe",
         layer = "europe_adm0",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)

# ................................
# ................................
# dissolve eur polygon 

# now dissolve
eur@data$eur <- 1

eur <- aggregate(eur, c("eur"))

plot(eur)

writeOGR(eur, 
         dsn = "data/gadm/europe",
         layer = "europe",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)

# ......................................
# ......................................
# create buffer area
buff <- SpatialPolygons(eur@polygons, 
                        proj4string = CRS(proj4string(eur)) )

buff <- gBuffer(buff, width = 0.5)

# buffer into polygons dataframe
ids <- sapply(slot(buff, "polygons"), function(x) slot(x, "ID"))
df <- data.frame(region="Europe", 
                 adm = "buffer_010arcmin", 
                 row.names=ids)

# Create a spatial polygon data frame (includes shp attributes)
spdf <- SpatialPolygonsDataFrame(buff, df)

plot(spdf)

writeOGR(spdf, 
         dsn = "data/gadm/europe",
         layer = "europe_buffer_010arcmin",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)

# .........................................
# ........................................
# create points around coast

buff <- SpatialPolygons(eur@polygons, 
                        proj4string = CRS(proj4string(eur)))

buff <- gBuffer(buff, width = -0.55)

e <- raster::intersect(buff, eur)

# set eur as raster
ext <- floor(extent(eur))
proj <- proj4string(buff)

r <- raster(ext, res = 0.05)

eur <- rasterize(eur, r, field = 1)

# and the inverse buffer
e <- rasterize(e, r, field = 1)

# now mask eur using the inverse raster
m <- mask(eur, e, inverse = TRUE)

p <- as.data.frame(m, xy = TRUE)

p <- p[!is.na(p[,3]),]

# remove russia east
p <- p[p[[1]] <= xmax(ext)-1, ]

p <- p[,-3]

p <- SpatialPoints(p, proj4string = CRS(proj))

df <- as.data.frame(p@coords)

pdf <- SpatialPointsDataFrame(p, df)

plot(pdf)

writeOGR(pdf, 
         dsn = "data/gadm/europe",
         layer = "europe_coast_points",
         driver = "ESRI Shapefile",
         overwrite_layer = TRUE)

# .....................................
# .....................................
# Ecorregions ####

eco <- st_read("/Volumes/BioversityInt/shapefiles/tnc_terr_ecoregions/tnc_terr_ecoregions.shp")
eco <- st_as_sf(eco)

eur <- st_read("data/gadm/europe/europe.shp")
eur <- st_as_sf(eur)

st_crop(eco, eur)


