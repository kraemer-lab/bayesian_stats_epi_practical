library(sf)

### Example on how to load and plot shp files
### Improving plots' labels and saving them is part of the task


gdf <- st_read("./shapefiles/sim_country.shp")

cases_map <- gdf[c("geometry", "cases")]

plot(cases_map)


temp_map <- gdf[c("geometry", "temperatur")]

plot(temp_map)