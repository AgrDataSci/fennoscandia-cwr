# return number of decimals 
# https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r/5173906
.decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}


substrRight <- function(x, n){
  substr(x, nchar(x) - n + 1, nchar(x))
}


#' Plot map using mapview
#' @param data a data frame
#' @param coords index of data for the lonlat coordinates
#' @param add any additional index for colunms in data to add to the map
plot_map <- function(data, coords = NULL, add = NULL, ...) {
  
  lonlat <- data[, c(coords, add)]
  
  lonlat <- stats::na.omit(lonlat)
  
  lonlat <- suppressWarnings(
    sf::st_as_sf(lonlat, coords = c("lon","lat"), crs = 4326)
  )
  
  suppressWarnings(
    mapview::mapview(lonlat, ...)
  )
}



