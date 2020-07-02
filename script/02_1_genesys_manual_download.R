# Combine data downloaded manually from Genesys
# K. de Sousa 
# Inland Norway University
#..........................................
#..........................................
# Packages ####
library("data.table")

sessioninfo::session_info()
# write session info
capture.output(sessioninfo::session_info(),
               file = "script/session_info/02_1_genesys_manual_download.txt")


dt <- fread("data/species_names.csv")

taxa <- dt[, .(genus,species,acronym)]

a <- taxa$acronym

f <- list.dirs("data/raw/genesys_manual_download")

gen <- data.frame()

for(i in seq_along(a)){
  
  f <- paste0("data/raw/genesys_manual_download/",
              a[i], "/",
              "geo.csv")
  
  f <- fread(f)
  
  f <- f[,.(longitude, latitude)]

  names(f) <- c("lon", "lat")
  
  f[, acronym := a[i]]
  
  f[, genus := taxa[i, "genus"]]
  
  f[, species := taxa[i, "species"]]
  
  f <- na.omit(f)
  
  f[, country := NA]
  
  cat(a[i], " ", nrow(f), "\n")
  
  gen <- rbind(gen, f)
    
}


table(gen$acronym)

length(table(gen$acronym))

write.csv(gen, file = "data/raw/genesys_occurrences_manual_download.csv", row.names = FALSE)
