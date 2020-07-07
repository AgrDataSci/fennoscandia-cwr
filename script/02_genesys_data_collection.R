# Get Genesys occurrence data for distribution analysis 
# updated by K. de Sousa 
# Inland Norway University
#..........................................
#..........................................
# Packages ####
library("genesysr")
library("data.table")

sessioninfo::session_info()
# write session info
capture.output(sessioninfo::session_info(),
               file = "script/session_info/02_genesys_data_collection.txt")

#..........................................
#..........................................
# Data ####

df <- fread("data/species_names.csv")

taxa <- df[ ,c("genus","species", "acronym")]

taxa <- split(taxa, 1:nrow(taxa))

#..........................................
#..........................................
# Set Genesys query ####

setup_sandbox()

user_login()

gen_df <- list()

for(i in seq_along(taxa)) {

  cat(i, " ",  as.vector(unlist(taxa[[i]])), "\n")

  genus <- as.character(as.matrix(taxa[[i]][,"genus"]))
  sp <- as.character(as.matrix(taxa[[i]][,"species"]))

  tax <- list(genus = genus, species = sp)

  call <- get_accessions(filters = mcpd_filter(#SAMPSTAT = c(100, 110, 120, 130, 200),
                                               GENUS = tax$genus,
                                               SPECIES = tax$species))

  if(dim(call)[[1]] > 0) {
    call$acronym <- taxa[[i]]$acronym
  }

  gen_df[[i]] <- call

}

save(gen_df, file = "data/raw/genesys_raw.rda")

#..........................................
#..........................................
# Combine data from genesys query ####

vars <- c("acronym","accessionNumber","acquisitionDate","geo.latitude","geo.longitude",
          "coll.collSite","countryOfOrigin.code2","id","inSvalbard","institute.acronym",
          "taxonomy.genus","taxonomy.species","taxonomy.subtaxa","institute.fullName",
          "institute.country.name","institute.country.code2")

gen_sub <- data.frame(matrix(NA, 
                             ncol = length(vars),
                             nrow = 0,
                             dimnames = list(NULL, vars)))


for (i in seq_along(gen_df)) {
  
  if (dim(gen_df[[i]])[1] == 0) { 
    next
  }
  
  # take the data
  x <- gen_df[[i]]
  
  # select the variable available in the data
  in_x <- vars %in% names(x)
  
  # if any missing variable, then add it as NAs
  if (any(!in_x)) {
    
    miss <- vars[!in_x]
    
    miss <- data.frame(matrix(NA, 
                              ncol = length(miss),
                              nrow = nrow(x),
                              dimnames = list(1:nrow(x), miss)))
    
    x <- cbind(x, miss)
    
  }
  
  x <- x[, vars]
  
  # bind with the main data
  gen_sub <- rbind(gen_sub, x)
  
  
}

gen_sub <- as.data.table(gen_sub)

sum(is.na(gen_sub$acronym))

table(gen_sub$acronym)

length(table(gen_sub$acronym))

write.csv(gen_sub, file = "data/raw/genesys_occurrences.csv", row.names = FALSE)
