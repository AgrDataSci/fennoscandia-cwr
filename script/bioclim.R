#...................................................
#...................................................
# Model species distribution using bioclimatic variables
# 
# Kaue de Sousa
# Inland Norway University
# 
# ...................................................
# ...................................................
# Packages ####
library("data.table")
library("blockCV")
library("rgdal")
library("raster")
library("dismo")
library("maptools")
library("sf")

#...................................................
#...................................................
# Data ####
# the BiodiversityR saves outputs on the current working directory
# get the parent wd to return here if needed
parentwd <- getwd()

outputwd <- "processing/enm/"
dir.create(outputwd, showWarnings = FALSE, recursive = TRUE)

# bioclimatic variables
bio <- list.files("data/bioclim",
                  pattern = ".tif$",
                  full.names = TRUE)

bio <- stack(bio)

# define projection and extension
myproj <- proj4string(bio)
myext  <- extent(bio) 
myres  <- res(bio)

myext[1] <- -10.5
myext[2] <- 48

bio <- crop(bio, myext)

# passport data
dt <- fread("data/passport_data.csv")

sp <- sort(unique(dt$acronym))
sp

i=1


sp_i <- sp[[i]]

coords <- dt[dt$acronym == sp_i, ]

coords <- coords[, .(lon, lat)]

# create a buffer with presence points to avoid background
# points too close to presence
lonlatb <- st_as_sf(coords, 
                    coords = c("lon", "lat"), 
                    crs = 4326)

# set the buffer around the points
lonlatb <- suppressWarnings(st_buffer(lonlatb,
                                      dist = 1))

lonlatb <- st_union(lonlatb)

lonlatb <- as_Spatial(lonlatb)

lonlatb <- rasterize(lonlatb, bio[[1]], field = 1, background = NA)

bg_mask <- mask(bio[[1]], lonlatb, inverse = TRUE)

# create background points over the area
nbg <- ceiling(nrow(coords))
bg <- randomPoints(bg_mask, 
                   n = nbg, 
                   p = coords, 
                   ext = myext, 
                   extf = 1.25)
bg <- as.data.frame(bg)
names(bg) <- c("lon", "lat")

# reduce sampling bias removing points within the same grid cell
r <- raster(myext)
# set the resolution of the cells to ~ 6 arc-min
res(r) <- c(1, 1)

coords <- as.data.frame(coords)

coords <- gridSample(coords, r, n = 1)

k <- 15

set.seed(123)
group <- kfold(coords, k)

set.seed(321)
group2 <- kfold(bg, k)


layers <- list()
aucs <- rep(NA, k)
thres <- rep(NA, k)

for(j in seq_len(k)) {
  print(j)
  
  pres_train <- coords[group != j , c("lon","lat")]
  back_train <- bg[group2 != j, c("lon","lat")]
  
  pres_test <- coords[group == j, c("lon","lat")]
  back_test <- bg[group2 == j, c("lon","lat")]
  
  bc <- bioclim(bio, pres_train)
  
  e <- evaluate(pres_test, back_test, bc, bio)
  
  aucs[j] <- e@auc
  
  tr <- threshold(e, 'equal_sens_spec', sensitivity = 0.9)

  thres[j] <- tr
    
  pb <- predict(bio, bc, ext = myext)
  
  pb[pb[] < tr ] <- NA
  
  layers[[j]] <- pb
  
}

x <- stack(layers)

x <- calc(x, mean)

r <- raster(bio, 1)
plot(!is.na(r), col=c('white', 'light grey'), legend=FALSE)
plot(x, add = TRUE)

x

