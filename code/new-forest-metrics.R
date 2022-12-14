# Include the analysis code I wrote elsewhere.
source("factored.R")

# Set up the suffix for files produced.
file_suffix <- "_background_new-forest"

# Coordinates of the edges of the New Forest
longmin <- -1.8
longmax <- -1.4
latmin <- 50.75
latmax <- 51

# Set patch cutoff, search radius, and max neighbours
patch_cutoff <- 0.9
radius <- Inf
k <- 100

# Get the raster.
gis_data <- fetch_raster(latmin, latmax, longmin, longmax)

# Extract points for aggregated land cells.
sites <- terra::as.points(
  gis_data, values = FALSE,
  na.rm = TRUE
)
site_locs <- data.frame(
  Longitude = crds(sites)[, 1],
  Latitude = crds(sites)[, 2]
)

# Extract the raw connectivity metrics -- this will take some time.
conn <- get_conn_metrics(
  gis_data, site_locs, patch_cutoff,
  radius, k, file_suffix
)

## MAKE THE FIGURE

####################################
## Calculate connectivity metrics ##
####################################

# Convert area from sq. m to sq. km, and distance to km
conn$area <- conn$area / 1e6
conn$dist <- conn$dist / 1e3

# Add the penetrability index along each path (using the BII).
conn$pen_index <- calculate_index2(conn[7:12], "forests")

# Add the weighted path length.
conn$weighted_dist <- conn$dist * conn$pen_index

# Group the connectivity metric by site, to calculate the aggregate metrics.
conn2 <- conn %>%
  group_by(site, Latitude, Longitude) %>%
  dplyr::summarise(
    # Calculate the area of the patch the site is in
    site_area = ifelse(
      !is.na(min(dist)) &    # at least one nearby patch
        min(dist) == 0,      # within a patch
      patch_cutoff * area[dist == min(dist)],
      ifelse(
        all(is.na(pri_site)),   # no site level data
        NA,
        mean(pri_site[!is.na(pri_site)]) + mean(sec_site[!is.na(sec_site)])
      )
    ),
    # Unweighted distance to nearest neighbour
    d_nn_uw = min(dist),
    # Weighted distance to nearest neighbour
    d_nn_w = min(weighted_dist),
    # Unweighted buffer, radius 2km
    b2_uw = patch_cutoff * sum(area[which(dist < 2)]),
    # Unweighted buffer, radius 5km
    b5_uw = patch_cutoff * sum(area[which(dist < 5)]),
    # Unweighted buffer, radius 7km
    b7_uw = patch_cutoff * sum(area[which(dist < 7)]),
    # Unweighted buffer, radius 10km
    b10_uw = patch_cutoff * sum(area[which(dist < 10)]),
    # Weighted buffer, radius 3km
    b2_w = patch_cutoff * sum(area[which(weighted_dist < 2)]),
    # Weighted buffer, radius 5km
    b5_w = patch_cutoff * sum(area[which(weighted_dist < 5)]),
    # Weighted buffer, radius 7km
    b7_w = patch_cutoff * sum(area[which(weighted_dist < 7)]),
    # Weighted buffer, radius 10km
    b10_w = patch_cutoff * sum(area[which(weighted_dist < 10)]),
    # Unweighted and weighted incidence function model. We want to take account
    # of the area of the patch in which the site is situated -- and for that I
    # will use the coverage at the site, or the patch area if the site is within
    # a patch. This will be a bit different in the case of the stretch goal
    # prioritisation map, since there we will be interested in the _potential_
    # patch size. Use a conservative estimate for the coverage within patches
    # (the cutoff) to scale the area.
    ifm_uw = ifelse(   # if there are no neighbours, this is 0.
      is.na(id),
      0,
      patch_cutoff * sum(
        exp(-dist[which(dist >= 0)]) * area[which(dist >= 0)]
      )
    ),
    ifm_w = ifelse(
      is.na(id),
      0,
      patch_cutoff * sum(
        exp(-weighted_dist[which(dist >= 0)]) * area[which(dist >= 0)]
      )
    )
  ) %>%
  ungroup() %>%
  unique()   # remove the duplicate rows so we have a single row per site

# Log the connectivity metrics. Add a small increment to avoid 0s.
conn2 <- mutate(conn2, across(d_nn_uw:ifm_w, \(x) log(x + 1e-5)))

## RASTERS
# Create a points vector object
points <- terra::vect(
  conn2[-1], geom = c("Longitude", "Latitude"),
  crs = "WGS84"
)
# Get the patch polygons.
polys <- get_patch_polys(gis_data, patch_cutoff)

make_map <- function(points, polys, metric, flip = FALSE) {
  for_map <- terra::rasterize(
    points, gis_data, field = metric
  )
  # Site
  slg <- -1.653148
  slt <- 50.9034
  # Plot
  pdf(
    paste(
      "../results/new-forest-", metric, ".pdf",
      sep = ""
    )
  )
  plot(
    if (flip) -1 * for_map$lyr1 else for_map$lyr1,
    axes = FALSE, mar = c(3.1, 3.1, 3.1, 4.1),
    legend = FALSE
  )
  plot(polys, col = 'black', add = TRUE)
  points(slg, slt, col = 'red', pch=16, cex=2)
  dev.off()
}

make_map(points, polys, "b5_uw")
make_map(points, polys, "d_nn_uw", flip = TRUE)
make_map(points, polys, "ifm_uw")
