# ...............................................
# ...............................................
# Make table with available coordinates per species
# ...............................................
# ...............................................
# Packages
library("tidyverse")
library("magrittr")
library("raster")
library("sf")

# # Present in target area (GBIF)	GBIF cluster			
# # Genesys total wild 
# # NOR wild Genesys	
# # SWE wild Genesys	
# # FIN wild Genesys	
# # DENM wild Genesys	
# # Main provenance of material 

# Genesys data
gen <- load("data/raw/genesys_raw.rda")

# GBIF data
gbif <- read_csv("data/passport_data.csv")
gbif %<>% 
  filter(source == "gbif")


# ADM data
# shape with adm units
adm <- st_read("data/gadm/europe/europe_adm0.shp")
adm <- st_as_sf(adm)
adm <- adm[-c(2:3)]
adm

myext <- extent(adm)
myext[1] <- -10.5
myext[2] <- 48

adm <- st_crop(adm, myext)

# Euro buffer
buf <- st_read("data/gadm/europe/europe_buffer_010arcmin.shp")
buf <- st_as_sf(buf)
buf <- buf[-c(1)]
buf <- st_crop(buf, myext)

# .......................................
# .......................................
# Add country data to GBIF
# remove points outside bbox 
keep <- gbif$lon > myext[1]

keep <- gbif$lon < myext[2] & keep

gbif <- gbif[keep, ]

coord <- st_as_sf(gbif, coords = c("lon", "lat"), crs = 4326)

x <- st_intersects(coord, adm)

gbif$adm0 <- as.character(adm$ADMIN[unlist(x)])

# .......................................
# .......................................
# Same for Genesys
gen <- lapply(gen_df, function(x){
  x[,c("acronym","countryOfOrigin.name","geo.latitude", "geo.longitude","inSvalbard")]
})

gen <- do.call(rbind, gen)

gen %<>% 
  as_tibble() %>% 
  rename(country = countryOfOrigin.name,
         lon = geo.longitude,
         lat = geo.latitude) %>% 
  mutate(country = ifelse(country == "", "Unknown", country))

sort(unique(gen$country))

gen$country[gen$country=="Czechoslovakia, Czechoslovak Socialist Republic"] <- "Czechia"
gen$country[gen$country=="Czechia"] <- "Czechia"
gen$country[gen$country=="German Democratic Republic"] <- "Germany"
gen$country[gen$country=="Romania, Socialist Republic of"] <- "Romania"
gen$country[gen$country=="Serbia and Montenegro"] <- "Serbia"
gen$country[gen$country=="Yugoslavia, Socialist Federal Republic of"] <- "Serbia"
gen$country[gen$country=="USSR, Union of Soviet Socialist Republics"] <- "Russian Federation"

# get list of iso country names
ct <- read_csv("data/country_iso.csv")
ct %<>% 
  filter(region == "Europe") %>% 
  dplyr::select(name:alpha3)

keep <- gen$country %in% union(ct$name, "Unknown")

gen <- gen[keep, ]

# remove points out the adm area
# first create an index for the entries without coordinates
naindex <- is.na(gen$lat) | is.na(gen$lon)

# coordinates into sf object
coord <- st_as_sf(gen[!naindex,], coords = c("lon","lat"), crs = 4326)

# interset the coordinates with the adm area
x <- st_intersects(coord, buf)
# check for the empty values (outside the area)
x <- lapply(x, function(y){
  ifelse(length(y)==0,FALSE,TRUE)
})
# vector with TRUE for coordinates inside the area
x <- unlist(x)

# send it back to the main df 
# combine the coordinates within the adm area and the entries
# with Unknown country
gen <- gen[sort(union(which(!naindex)[x], which(naindex))),]

# remove Unknown country without coordinates
out <- (is.na(gen$lon) + is.na(gen$lat)) > 0 & gen$country == "Unknown"
gen <- gen[!out, ]

out <- (is.na(gen$lon) + is.na(gen$lat)) == 1

gen <- gen[!out, ]

# check country names and add names to Unknown with coordinates 
naindex <- is.na(gen$lat)

coord <- st_as_sf(gen[!naindex, ], coords = c("lon", "lat"), crs = 4326)

x <- st_intersects(coord, adm)

x <- lapply(x, function(y){
  ifelse(length(y)==0, NA, as.character(adm$ADMIN[y]))
})

x <- unlist(x)

gen$adm <- NA

gen[!naindex, "adm"] <- x

gen$country <- ifelse(gen$country == "Unknown", gen$adm, gen$country)

gen$country[gen$country=="Russian Federation"] <- "Russia"
gen$country[gen$country=="Republic of Serbia"] <- "Serbia"

x <- cbind(acronym = c("RIBNIG", "RIBSPI", "ALLSCO","RIBPET","RIBPET","RIBPET"),
           country = c("Sweden", "Poland", "Russia", "Georgia","Georgia","Georgia"),
           lat = c(66.1, 53.1936,	41.316667,42.56333333,NA,NA),
           lon = c(22.1, 17.8902, 47.916667,45.03722222,NA,NA),
           inSvalbard = 0,
           adm = NA)

gen <- rbind(gen, x)

gen <- gen[!is.na(gen$country), ]

#keep <- is.na(gen$lat) + is.na(gen$lon) == 0
#gen <- gen[keep, ]

# .............................................
# .............................................
# ............................................
# Make the summary table

main <- 
  gen %>% 
  group_by(acronym) %>% 
  mutate(inSvalbard = as.integer(inSvalbard)) %>% 
  summarise(Genesys = length(acronym),
            Georeferenced = sum(!is.na(lon)),
            inSvalbard = sum(inSvalbard == 1 & !is.na(lon)),
            Norway = sum(country == "Norway" & !is.na(lon)),
            Sweden = sum(country == "Sweden" & !is.na(lon)),
            Finland = sum(country == "Finland" & !is.na(lon)),
            Denmark = sum(country == "Denmark" & !is.na(lon)),
            Russia = sum(country == "Russia"  & !is.na(lon)))
  
main

keep <- is.na(gen$lat) + is.na(gen$lon) == 0
gen <- gen[keep, ]

max_country <-
  gen %>% 
  group_by(acronym) %>% 
  count(country) %>% 
  filter(n == max(n)) %>% 
  distinct(acronym, .keep_all = TRUE) %>% 
  mutate(`Main provenance` = paste0(country, " (n=", n, ")")) %>% 
  dplyr::select(-country, -n)

sum_gbif <-
gbif %>% 
  group_by(acronym) %>% 
  count(acronym) %>% 
  rename(GBIF = n)

main %<>% 
  inner_join(., max_country, by = "acronym") %>% 
  inner_join(., sum_gbif, by = "acronym")


main <- main[union(c("acronym","GBIF", "Genesys"), names(main))]

spp <- read_csv("data/species_names.csv")

spp %<>% 
  mutate(taxa = paste(genus, species, authority)) %>% 
  dplyr::select(family, taxa, acronym)


spp %<>% 
  inner_join(., main, by = "acronym")

output <- "output/summary_table/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

write_csv(spp, paste0(output, "summary_table.csv"))
