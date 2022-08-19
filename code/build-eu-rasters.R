source("factored.R")

### RASTERS AND MAPS
conn <- readRDS("../data/europe_background_conn_2022-08-19.rds")
patch_cutoff <- 0.9

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
        exp(-dist[which(dist > 0)]) * area[which(dist > 0)]
      )
    ),
    ifm_w = ifelse(
      is.na(id),
      0,
      patch_cutoff * sum(
        exp(-weighted_dist[which(dist > 0)]) * area[which(dist > 0)]
      )
    ),
    ifm_w_05 = ifelse(
      is.na(id),
      0,
      patch_cutoff * sum(
        exp(-0.5 * weighted_dist[which(dist > 0)]) * area[which(dist > 0)]
      )
    )
  ) %>%
  ungroup() %>%
  unique()   # remove the duplicate rows so we have a single row per site

# Log the connectivity metrics. Add a small increment to avoid 0s.
conn2 <- mutate(conn2, across(d_nn_uw:ifm_w_05, \(x) log(x + 1e-5)))

## RASTERS
# Create a points vector object
points <- terra::vect(
  conn2[-1], geom = c("Longitude", "Latitude"),
  crs = "WGS84"
)

# Coordinates of the edges of Europe.
longmin <- -11
longmax <- 30
latmin <- 35
latmax <- 70
# Get the original raster.
gis_data <- fetch_raster(latmin, latmax, longmin, longmax)
e <- terra::ext(-10.5, 29.5, 35.5, 69.5)
gis_data2 <- terra::aggregate(
  terra::crop(gis_data, e),
  fact = 30,  # 30km points
  na.rm = TRUE
)

for (metric in names(points)[-1]) {
  terra::rasterize(
    points, gis_data2, field = metric, filename = paste(
      "../results/rasters/background-eu-", metric, ".tif",
      sep = ""
    ), overwrite = TRUE
  )
}

## Final rasters of interest
for_maps <- c(
  terra::rast("../results/rasters/background-eu-ifm_w.tif"),
  terra::rast("../results/rasters/background-eu-b7_w.tif")
)     # applies WGS84 by itself
names(for_maps) <- c("ifm_w", "b7_w")

par(mfrow=c(1,2))
plot(for_maps$ifm_w)
plot(gis_data2$primary + gis_data2$secondary)

par(mfrow=c(1,2))
plot(for_maps$b7_w)
plot(gis_data2$primary + gis_data2$secondary)
