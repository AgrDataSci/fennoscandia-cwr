#...................................................
# Select bioclimatic variables for modelling
# using VIF analysis

#...................................................
#...................................................
# Packages ####
library("data.table")
library("raster")
library("dismo")
library("BiodiversityR")
library("car")

sessioninfo::session_info()
# write session info
capture.output(sessioninfo::session_info(),
               file = "script/session_info/04_vif_analysis_variable_selection.txt")

#...................................................
#...................................................
# Data ####

# bioclimatic variables
bio <- list.files("data/bioclim",
                  pattern = ".tif$",
                  full.names = TRUE)

bio <- stack(bio)
names(bio)
# define projection and extension
myproj <- proj4string(bio)
myext  <- extent(bio)
myres  <- res(bio)

# species acronyms
list.files("data")

# passport data
df <- fread("data/passport_data.csv")

# .......................................
# .......................................
# Set background points ####
xy <- df[, c("lon", "lat")]

xy <- unique(xy, by = c("lon", "lat"))

set.seed(123)
bg <-randomPoints(bio[[1]], 
                  n = 10000, 
                  ext = myext, 
                  extf = 1.25)

plot(bio[[1]])
points(bg)
#...................................................
#...................................................
# Variable selection with VIF ####
vif <- ensemble.VIF(
  x = bio,
  a = xy,
  an = bg,
  VIF.max = 10,
  keep = NULL,
  factors = NULL,
  dummy.vars = NULL
)

# save outputs
output <- "processing/vif/"

dir.create(output,
           recursive = TRUE,
           showWarnings = FALSE)

save(vif, file = paste0(output, "vif_results.rda"))

# remove files not selected by vif analysis
out <- vif$var.drops

file.remove(paste0("data/bioclim/", out, ".tif"))
