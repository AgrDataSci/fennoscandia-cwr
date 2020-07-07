# GADM

This is Global Administrative Areas data downloaded from <https://gadm.org/index.html>. Country polygons were downloaded and combined into the mainland Europe. Overseas Territories were removed. Details in `script/helper02_combine_polygons_eur.R`

Citation:
>Global Administrative Areas (2012). GADM database of Global Administrative Areas, version 2.0. [online] URL: www.gadm.org.

# Ecorregions

Terrestrial Ecoregions shapefile downloaded from <http://maps.tnc.org/gis_data.html>. The shapefile was clipped in ArcGIS using the file `gadm/europe/europe.shp` as baseline.


Citation:
>David M. Olson, Eric Dinerstein, Eric D. Wikramanayake, Neil D. Burgess, George V. N. Powell, Emma C. Underwood, Jennifer A. D'amico, Illanga Itoua, Holly E. Strand, John C. Morrison, Colby J. Loucks, Thomas F. Allnutt, Taylor H. Ricketts, Yumiko Kura, John F. Lamoreux, Wesley W. Wettengel, Prashant Hedao, Kenneth R. Kassem, Terrestrial Ecoregions of the World: A New Map of Life on Earth: A new global map of terrestrial ecoregions provides an innovative tool for conserving biodiversity, BioScience 51(11):933-938. https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2


# Protected areas

World Database on Protected Areas downloaded from <https://www.protectedplanet.net/>. The shapefile was clipped in R using the file `gadm/europe/europe.shp` as baseline. However the files are too large and were not uploaded to GitHub.

Citation:
>UNEP-WCMC (2020), The World Database on Protected Areas (WDPA) statistics. Cambridge, UK: UNEP- WCMC. Accessed on: [05/07/2020]. Available at: <https://www.protectedplanet.net/>

