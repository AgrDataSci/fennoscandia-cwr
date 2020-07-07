# ...............................................
# ...............................................
# Process layers from species distribution models
# ...............................................
# ...............................................
# Packages
library("data.table")
library("ggplot2")
library("patchwork")
library("raster")
library("sp")
library("sf")

sessioninfo::session_info()
# write session info
capture.output(sessioninfo::session_info(),
               file = "script/session_info/06_process_layers.txt")

# ...............................................
# ...............................................
# Files and parameters
path <- "processing/enm/"

# file with spp acronyms and uses
spnames <- fread("data/species_names.csv")

# read the passport data
pass <- fread("data/passport_data.csv")

# shape with adm units
adm <- st_read("data/gadm/europe/europe_adm0.shp")
adm <- st_as_sf(adm)
adm <- adm[-c(2:3)]
adm

# ...............................................
# ...............................................
# get species names from processed models
sp <- spnames$acronym
n  <- max(seq_along(sp))

# filter spnames and passport data
keep <- spnames$acronym %in% sp

spnames <- spnames[keep, ]

keep <- pass$acronym %in% sp

pass <- pass[keep, ]

# each species layer has its own extent based on the max hull for the 
# presence points used take a raster that includes all the regional
# area used here to create a new layer with the same extent for all species
# use one of the bioclim layers and set all values as zero
f <- list.files("data/bioclim", pattern = ".tif", full.names = TRUE)[[2]]
eur <- raster(f)
eur[eur[] != 0] <- 0

myext <- extent(eur)
#myext[1] <- -10.5
myext[2] <- 48

eur <- crop(eur, myext)

adm <- st_crop(adm, myext)

myproj <- as.character(crs(eur))

plot(adm)

# ...............................................
# ...............................................
# Run over current presence ####
# list with presence layers
sp_p <- list()

pb <- txtProgressBar(min = 1, max = n, initial = 0)

for(i in seq_along(sp)) {
  
  # layer with presence
  p_i <- paste0(path, "/", 
                sp[i], 
                "/ensembles/presence/")
  
  p_i <- stack(list.files(p_i,
                          pattern = "current.grd",
                          full.names = TRUE))

  # presence is defined as 1,
  # set all values different tha 1 as NA
  p_i[p_i[] != 1 ] <- NA

  crs(p_i) <- myproj

  p_i <- crop(p_i, myext)
  
  # reconstruct layer using the regional layer as baseline
  p_i <- mosaic(eur, p_i, fun = sum)
  # set the regional layer as a mask
  p_i <- mask(p_i, eur)

  sp_p[[sp[i]]] <- p_i

  setTxtProgressBar(pb, i)
  
}


# ...............................................
# ...............................................
# Run over current distribution ####
# list with distribution layers
sp_d <- list()

pb <- txtProgressBar(min = 1, max = n, initial = 0)

for(i in seq_along(sp)) {
  
  d_i <- paste0(path, "/",
                sp[i],
                "/ensembles/suitability/")
  
  d_i <- stack(list.files(d_i,
                          pattern = "current.grd",
                          full.names = TRUE))
  
  crs(d_i) <- myproj
  
  d_i <- crop(d_i, myext)
  
  # reconstruct layer using the regional layer as baseline
  d_i <- mosaic(eur, d_i, fun = sum)
  # set the regional layer as a mask
  d_i <- mask(d_i, eur)
  
  # get the presence layers and set it as a mask
  # for the distribution threshold
  p <- sp_p[[i]]
  p[p[] != 1] <- NA
  
  d_i <- mask(d_i, p)
  
  sp_d[[sp[i]]]  <- d_i
  
  setTxtProgressBar(pb, i)
  
}

# ...............................................
# ...............................................
# Create maps ###

# define colours for map of gradient of distribution
colpall <- colorRampPalette(c("#FFFFFF", "#FFFF80", 
                              "#38E009","#1A93AB", 
                              "#0C1078"))

# order spp names by uses
spnames <- spnames[order(spnames$acronym), ]

# apply the same order across the list of maps
sp_d <- sp_d[spnames$acronym]

rm(sp)

# remove points outside bbox 
keep <- pass$lon > myext[1]

keep <- pass$lon < myext[2] & keep

pass <- pass[keep, ]

# Plot individual maps
output <- "output/individual_maps/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

for(i in seq_along(sp_d)) {
  
  print(spnames[i,"taxa"])
  
  r <- sp_d[[i]]
  r <- as.data.frame(r, xy = TRUE)
  r <- r[!is.na(r[, "layer"]), ]
  r[, "layer"] <- round(r[, "layer"] / 1000, 2)
  
  gen <- pass[pass$acronym == spnames[[i,"acronym"]], ]
  gen <- gen[gen$source == "genesys"]
  gen <- gen[, .(lon, lat)]
  
  # if(nrow(gen) > 5){
  #   x <- extract(sp_d[[i]], as.data.frame(gen))
  #   gen <- gen[!is.na(x), ]
  # }
  
  acronym <- names(sp_d[i])
  acronym <- which(spnames$acronym %in% acronym)
  species <- paste(spnames[acronym,.(genus, species)], collapse = " ")
  
  p <- ggplot() +
    geom_tile(r, mapping = aes(x = x, y = y, fill = layer)) +
    geom_sf(adm$geometry, mapping = aes(), colour = "black", fill = NA) +
    geom_point(gen, mapping = aes(x = lon, y = lat), size = 1.2, col = "red", pch = 18) +
    scale_fill_gradientn(name = NULL, 
                         colours = colpall(10),
                         limits = c(0, 1)) +
    theme_void() + 
    labs(title = species) +
    theme(legend.position = "right",
          plot.title = element_text(size = 14, 
                                    colour = "black", 
                                    face = "italic"),
          legend.text = element_text(size = 14),
          plot.margin = unit(c(1,5,1,1), "mm"))
  p
  ggsave(paste0(output, gsub(" ","_", species),  ".png"),
         plot = p,
         width = 20,
         height = 20,
         units = "cm")
  
}

# .........................................
# .........................................
# Plot all presences to see the areas with higher diversity ####

output <- "output/diversity/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

div <- stack(sp_p)
div <- calc(div, sum)

r <- as.data.frame(div, xy = TRUE)
r <- r[!is.na(r$layer) & r$layer > 0, ]

colpall <- colorRampPalette(c('#ffffcc','#ffeda0','#fed976',
                              '#feb24c','#fd8d3c','#fc4e2a',
                              '#e31a1c','#bd0026','#800026'))

p <- 
ggplot() +
  geom_tile(r, mapping = aes(x = x, y = y, fill = layer)) +
  geom_sf(adm$geometry, mapping = aes(), colour = "black", fill = NA) +
  scale_fill_gradientn(name = NULL,
                       colours = colpall(18),
                       labels = c(5,15,25,35),
                       breaks = c(5,15,25,35)) +
  theme_void() +
  theme(legend.text = element_text(size = 14),
        plot.margin = unit(c(1,5,1,1), "mm"))

ggsave(paste0(output, "diversity.png"),
       plot = p,
       width = 20,
       height = 20,
       dpi = 500,
       units = "cm")

# .........................................
# .........................................
# Plot by uses and all together ####
colpall <- colorRampPalette(c("#FFFFFF", "#FFFF80", 
                              "#38E009","#1A93AB", 
                              "#0C1078"))
plots <- list()

for(i in seq_along(sp_d)) {
  
  print(spnames[i,"taxa"])
  
  r <- sp_d[[i]]
  r <- as.data.frame(r, xy = TRUE)
  r <- r[!is.na(r[, "layer"]), ]
  r[, "layer"] <- round(r[, "layer"] / 1000, 2)
  
  gen <- pass[pass$acronym == spnames[[i,"acronym"]], ]
  gen <- gen[gen$source == "genesys"]
  gen <- gen[, .(lon, lat)]
  
  # if(nrow(gen) > 5){
  #   x <- extract(sp_d[[i]], as.data.frame(gen))
  #   gen <- gen[!is.na(x), ]
  # }
  
  legend <- "none"
  
  
  acronym <- names(sp_d[i])
  acronym <- which(spnames$acronym %in% acronym)
  spp <- spnames[acronym,.(genus)]
  spp <- substr(spp, start = 1, stop = 1)
  
  species <- paste0(spp, ". ", spnames[acronym,.(species)])
  
  
  p <- ggplot() +
    geom_tile(r, mapping = aes(x = x, y = y, fill = layer)) +
    geom_sf(adm$geometry, mapping = aes(), colour = "black", fill = NA) +
    geom_point(gen, mapping = aes(x = lon, y = lat), size = 0.6, col = "red", pch = 18) +
    scale_fill_gradientn(name = NULL, 
                         colours = colpall(10),
                         limits = c(0, 1)) +
    theme_void() + 
    labs(title = species) +
    theme(legend.position = legend,
          plot.title = element_text(size = 10, 
                                    colour = "black", 
                                    face = "italic"))
  
  
  
  plots[[i]] <- p
  
}

names(plots) <- names(sp_d)

table(spnames$use)

# leafy spp
leafy <- sort(spnames$acronym[spnames$use == "Leafy"])

leafy <- plots[leafy]

last <- length(leafy)

leafy[[last]] <-
  leafy[[last]] +
  theme(legend.position = "right")

p <- 
  leafy[[1]] + leafy[[2]] + leafy[[3]] + leafy[[4]] + 
  leafy[[5]] + leafy[[6]] + leafy[[7]] + leafy[[8]] + 
  leafy[[9]] + leafy[[10]] + leafy[[11]] + leafy[[12]] +
  leafy[[13]] + leafy[[14]] + leafy[[15]] + leafy[[16]] +
  leafy[[17]] + leafy[[18]] + leafy[[19]] + leafy[[20]] +
  plot_layout(ncol = 4, nrow = 5)

output <- "output/combined_map/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

ggsave(paste0(output, "leafy_map.png"),
       plot = p,
       width = 20,
       height = 25,
       dpi = 950,
       units = "cm")

# .................................
# .................................
# the other uses
other <- sort(spnames$acronym[spnames$use != "Leafy"])

other <- plots[other]

last <- length(other)

other[[last]] <-
  other[[last]] +
  theme(legend.position = "right")

p <- 
  other[[1]] + other[[2]] + other[[3]] + 
  other[[4]] + other[[5]] + other[[6]] + 
  other[[7]] + other[[8]] + other[[9]] + 
  other[[10]] + other[[11]] + other[[12]] +
  other[[13]] + other[[14]] + other[[15]] +
  plot_layout(ncol = 3, nrow = 5)

ggsave(paste0(output, "all_other_uses_map.png"),
       plot = p,
       width = 18,
       height = 25,
       dpi = 950,
       units = "cm")


