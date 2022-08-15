source("factored.R")

# Coordinates of the edges of the New Forest
longmin <- -1.8
longmax <- -1.4
latmin <- 50.75
latmax <- 51

# Set patch cutoff, search radius, and max neighbours
patch_cutoff <- 0.9
radius <- Inf
k <- 100

# Fetch and wrangle data.
data <- fetch_data(latmin, latmax, longmin, longmax)
gis_data <- fetch_raster(latmin, latmax, longmin, longmax)
wrangled <- wrangle_data(data)

# Get the raw connectivity metrics.
conn <- get_conn_metrics(
  gis_data, wrangled,
  radius=radius, k=k, patch_cutoff=patch_cutoff
)

# Calculate patches and lines.
coords <- data.frame(
  strsplit(unique(paste(
    wrangled$Longitude, wrangled$Latitude,
    sep = "#"
  )), split = "#")
)
site_locs <- data.frame(
  Longitude = as.numeric(coords[1, ]),
  Latitude = as.numeric(coords[2, ])
)
# Transform sites into points object.
sites <- terra::vect(
  site_locs, # source data frame
  geom = c("Longitude", "Latitude"), # columns for coordinates
  crs = "WGS84"
)
# Get the patch polygons.
polys <- get_patch_polys(gis_data, patch_cutoff)
# Calculate the distances to nearest habitat patch.
neighbours <- nngeo::st_nn(
  st_as_sf(sites), st_as_sf(polys),
  k = min(k, length(polys)),
  maxdist = radius,
  parallel = n_cores,
  returnDist = TRUE
)
# Calculate the lines between each site and its neighbor.
lines_between <- as(nngeo::st_connect(
  st_as_sf(sites), st_as_sf(polys),
  ids = neighbours[["nn"]]
), "SpatVector")

######################
## Make a nice plot ##
######################

pdf(
  "../results/metrics-figure.pdf",
  width = 10, height = 8 # inches
)

# Plot the habitat
plot(
  gis_data$primary + gis_data$secondary,
  axes = FALSE, mar = c(3.1, 3.1, 3.1, 4.1)
)

# Plot all patches
plot(polys, add=T, col='black')

# Highlight patches within 5km of site of interest
nearby <- conn[conn$site == 50 & conn$dist < 5000,]$id
plot(polys[nearby], add=T, col='blue')

# # Plot all sites
# points(sites, pch=4)

# Highlight the site of interest
points(sites[50], col='red', pch=16, cex=2)

# Plot a 5km buffer
plot(terra::buffer(sites[50], 5000), add=T)

# Add the connecting lines
line_ids <- as.numeric(rownames(conn[conn$site == 50,]))
connecting_lines <- lines_between[line_ids]
lines(connecting_lines)

# Highlight the line to the nearest neighbour
lines(connecting_lines[1], col='red', lwd=2)

dev.off()
