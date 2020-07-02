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
library("rgdal")
library("raster")
library("dismo")
library("rgeos")
library("gbm")
library("glmnet")
library("maxnet")
library("sf")
library("BiodiversityR")
library("PresenceAbsence")

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

# passport data
dt <- fread("data/passport_data.csv")

sp <- sort(unique(dt$acronym))
sp

#...................................................
#...................................................
# In case you need to stop the process here is the 
# point where the algorithm looks for species already 
# processed
filepattern <- paste(c("current"), collapse = "|")
nfiles <- length(c("current")) * 2

for (i in seq_along(sp)) {
  
  pres <- paste0(outputwd, sp[i], "/ensembles/presence/")
  pres <- grepl(filepattern, list.files(pres))
  
  suit <- paste0(outputwd, sp[i], "/ensembles/suitability/")
  suit <- grepl(filepattern, list.files(suit))
  
  # if folder doesnt contain these files then it should be removed
  if (sum(suit) != nfiles | sum(pres) != nfiles) {
    
    unlink(paste0(outputwd, sp[i]), recursive = TRUE)
    
  }
  
}

done <- list.dirs(outputwd)[-1]
done <- strsplit(done, "/")
done <- suppressWarnings(do.call("rbind", done)[,4])
done <- unique(done)


# filter and run the model for the remnant 
sp <- sp[!sp %in% done] 

#...................................................
#...................................................
# Run ensemble modelling ####
for (i in seq_along(sp)) {
  
  # create a dir for the species and work in that dir
  output <- paste0(outputwd, sp[i])
  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  
  setwd(output)
  
  k <- dt$acronym == sp[i]
  
  # subset data
  coord <- dt[k, ]
  
  message("\n Ensemble modelling for ", 
          unique(coord$acronym), "\n Time: ", 
          date(), "\n")
  
  coord <- coord[, .(lon, lat)]
  
  coord <- as.data.frame(coord)
  
  # create a convexHull to limit the model to the 
  # area where the species is actually present
  # calculate largest distance
  largedist <- pointDistance(coord, longlat = FALSE)
  largedist <- max(largedist, na.rm = TRUE)

  # make a convex hull 
  hull <- convHull(coord, lonlat = TRUE)
  
  # extent convex hull
  hull <- gBuffer(hull@polygons, 
                  width = 0.1 * largedist)
  
  crs(hull) <- myproj
  
  # convert into a raster
  r <- raster(hull, res = myres, ext = myext)

  hull <- rasterize(hull, r, field = 1,
                    background = NA)
  
  hull[hull == 0] <- NA
  
  # crop bioclim layers to fit this extention
  bio_i <- mask(bio, hull)
  
  bio_i <- crop(bio_i, hull)
  
  bio_i <- stack(bio_i)
  
  ext <- extent(bio_i)
  
  # create a buffer with presence points to avoid background
  # points too close to presence
  lonlatb <- st_as_sf(coord, 
                      coords = c("lon", "lat"), 
                      crs = 4326)
  
  # set the buffer around the points
  lonlatb <- suppressWarnings(st_buffer(lonlatb,
                                        dist = 1))
  
  lonlatb <- st_union(lonlatb)
  
  lonlatb <- as_Spatial(lonlatb)
  
  lonlatb <- rasterize(lonlatb, bio_i[[1]], field = 1, background = NA)

  bg_mask <- mask(bio_i[[1]], lonlatb, inverse = TRUE)
  
  # create background points over the area
  nbg <- ceiling(nrow(coord))
  bg <- randomPoints(bg_mask, 
                     n = nbg, 
                     p = coord, 
                     ext = ext, 
                     extf = 1.25)
  bg <- as.data.frame(bg)
  names(bg) <- c("lon", "lat")
  
  # reduce sampling bias removing points within the same grid cell
  r <- raster(ext)
  # set the resolution of the cells to ~ 6 arc-min
  res(r) <- c(1, 1)

  coord <- as.data.frame(coord)
  
  coord <- gridSample(coord, r, n = 1)
  
  # Run ensemble modelling
  # step 1: model calibration
  # here the function tests for the best algorithms
  # since the algorithms were previously selected,
  # a 3-fold cross-validation is performed to make sure that all 
  # pass the output.weights threshold
  message("\n Step 1: Calibrating algorithms \n", "Time: ", date(), "\n")
  set.seed(9999)
  enm_step1 <- ensemble.calibrate.weights(x = bio_i, 
                                          p = coord, 
                                          a = bg,
                                          k = 3,
                                          layer.drops = NULL,
                                          SINK = TRUE, 
                                          species.name = sp[[i]],
                                          BIOCLIM = 1,
                                          MAHAL = 0,
                                          DOMAIN = 1, 
                                          MAXNET = 1,
                                          GBM = 0,
                                          GAM = 0,
                                          GLM = 0,
                                          SVM = 0,
                                          RPART = 0, 
                                          GBMSTEP = 0,
                                          MAXENT = 0, 
                                          NNET = 0, 
                                          RF = 0, 
                                          EARTH = 0,
                                          GLMSTEP = 0, 
                                          GAMSTEP = 0, 
                                          MGCV = 0, 
                                          MGCVFIX = 0, 
                                          CF = 0, 
                                          FDA = 0,
                                          SVME = 0,
                                          ENSEMBLE.tune = TRUE, 
                                          PROBIT = TRUE,
                                          # see Liu et al (2013) doi:10.1111/jbi.12058
                                          threshold.method = "threshold2013.mean", 
                                          threshold.PresenceAbsence = TRUE, 
                                          ENSEMBLE.min = 0.7)
  
  # step 2: create models that will be used for the raster predictions
  # models with output.weights <0.05 are excluded
  output_weights <- enm_step1$output.weights
  output_weights[output_weights < 0.05] <- 0
  
  message("Step 2: Model species distribution with selected ENM algorithms \n")
  
  set.seed(9999)
  enm_step2 <- ensemble.calibrate.models(x = bio_i, 
                                         p = coord, 
                                         a = bg,
                                         k = 10, 
                                         layer.drops = NULL,
                                         SINK = TRUE, 
                                         species.name = sp[[i]],
                                         models.keep = TRUE,
                                         input.weights = output_weights,
                                         # see Liu et al (2013) doi:10.1111/jbi.12058
                                         threshold.method = "threshold2013.mean", 
                                         threshold.PresenceAbsence = TRUE, 
                                         ENSEMBLE.tune = FALSE, 
                                         ENSEMBLE.min = 0.7,
                                         PROBIT = TRUE,
                                         Yweights = "BIOMOD", 
                                         models.save = FALSE)
  
  
  # save AUCs
  auc <- data.frame(auc = enm_step2$AUC.testing)
  auc$model <- rownames(auc)
  auc <- auc[!grepl("MAHAL|ENSEMBLE", rownames(auc)), ]
  write.csv(auc, file = "outputs/auc_testing.csv")
  
  message("Step 3.1: Generate map of current distribution \n")
  #step3: use previously calibrated models to construct consensus layers
  ensemble_current <- ensemble.raster(xn = bio_i,
                                      models.list = enm_step2$models,
                                      input.weights = output_weights,
                                      thresholds = enm_step2$models$thresholds,
                                      SINK = TRUE,
                                      RASTER.species.name = sp[[i]], 
                                      RASTER.stack.name = "current")
  
  # remove working files created in the third step
  unlink("models", recursive = TRUE)
  unlink("ensembles/count", recursive = TRUE)
  file.remove(list.files(pattern = "working", full.names = TRUE))
  
  
  # return to parent wd
  setwd(parentwd)

}

message("Done at ", Sys.time())
