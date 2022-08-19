# Include the analysis code I wrote elsewhere.
source("factored.R")

# Set up the suffix for files produced.
file_suffix <- "_background"

# Coordinates of the edges of Europe.
longmin <- -11
longmax <- 30
latmin <- 35
latmax <- 70

patch_cutoff <- 0.9

# Get the raster.
gis_data <- fetch_raster(latmin, latmax, longmin, longmax)

# Extract points for aggregated land cells. Crop to try and avoid memory leaks
# around the edge? Cropping cuts ~500 30km sites
e <- terra::ext(-10.5, 29.5, 35.5, 69.5)
gis_data2 <- terra::aggregate(
  terra::crop(gis_data, e),
  fact = 30,  # 30km points   # TODO: also remove points that fall in sea?
  na.rm = TRUE
)
sites <- terra::as.points(
  gis_data2, values = FALSE,
  na.rm = TRUE
)
saveRDS(
  sites, file = paste(
    "../data/sites", file_suffix, ".rds",
    sep=""
  )
)

# # Retrieve
# sites <- readRDS(
#   paste(
#     "../data/sites", file_suffix, ".rds",
#     sep=""
#   )
# )
# sites <- terra::vect(sites)

# Extract the raw connectivity metrics -- this will take some time.
logger::log_info(
  "Found ", as.character(length(sites)), " sites for analysis"
)

# # Get the patch polygons.
# logger::log_info("Extracting polygons...")
# polys <- get_patch_polys(gis_data, patch_cutoff)
# saveRDS(
#   polys, file = paste(
#     "../data/polys", file_suffix, ".rds",
#     sep=""
#   )
# )

# Retrieve
polys <- readRDS(
  paste(
    "../data/polys", file_suffix, ".rds",
    sep=""
  )
)
polys <- terra::vect(polys)

# # Get the land-use cover in the cell of the site itself.
# logger::log_info("Collecting land cover data at sites...")
# lu_sites <- terra::extract(
#   gis_data, sites
# )
# saveRDS(
#   lu_sites, file = paste(
#     "../data/lu_sites", file_suffix, ".rds",
#     sep=""
#   )
# )

# Retrieve
lu_sites <- readRDS(
  paste(
    "../data/lu_sites", file_suffix, ".rds",
    sep=""
  )
)


radius <- 20000
k <- 20

# # Run garbage collection to avoid issues.
# logger::log_info("Running garbage collection")
# gc()
# # Calculate the distances to nearest habitat patch.
# logger::log_info("Finding neighbours...")
# neighbours <- nngeo::st_nn(
#   st_as_sf(sites), st_as_sf(polys),
#   k = min(k, length(polys)),
#   maxdist = radius,
#   returnDist = T
# )
# saveRDS(
#   neighbours, file = paste(
#     "../data/neighbours", file_suffix, ".rds",
#     sep=""
#   )
# )

# Retrieve
neighbours <- readRDS(
  paste(
    "../data/neighbours", file_suffix, ".rds",
    sep=""
  )
)

# Batching to dodge random nngeo memory leaks.
nns <- neighbours$nn
nn_len <- length(nns)
n_batches <- 20
batch_groups <- sample(letters[1:n_batches], nn_len, replace = TRUE)
batches <- split(nns, batch_groups)
batched_sites <- split(st_as_sf(sites), batch_groups)

# for (batch in letters[1:n_batches]) {
#   # Calculate the lines between each site and its neighbor.
#   lines_between <- as(nngeo::st_connect(
#     batched_sites[[batch]], st_as_sf(polys),
#     ids = batches[[batch]], k = k, maxdist = radius
#   ), "SpatVector")
#   # Extract the land-use values under the line, and summarise these.
#   logger::log_info("Extracting path data...")
#   vals <- terra::extract(
#     gis_data, lines_between
#   )
#   logger::log_info("Modifying vals...")
#   vals <- vals %>%
#     # Make sure we keep track of how many patches of water we cover.
#     mutate(
#       water = is.nan(primary),
#       across(, function(x) ifelse(is.nan(x), 0, x))
#     ) %>%
#     group_by(ID) %>%
#     # And sum the proportions of each land use along the lines.
#     summarise(across(, sum)) %>%
#     ungroup()
#   # Save progress in case of failure
#   logger::log_info("Saving data...")
#   saveRDS(
#     vals, file = paste(
#       "../data/vals", "_batch_", batch, file_suffix, ".rds",
#       sep=""
#     )
#   )
# }

# We have a bunch of batched data now, but the ordering of the paths is all over
# the place. We can handle this properly later.
batched_vals <- list()
for (batch in letters[1:n_batches]) {
  batched_vals[[batch]] <- readRDS(
    paste(
      "../data/vals", "_batch_", batch, file_suffix, ".rds",
      sep=""
    )
  )
}

# Get coordinates of the points
site_locs <- terra::crds(sites)
colnames(site_locs) <- c("Longitude", "Latitude")
# Put all of this data together.
res_cols <- c(
  # Site
  "site", "Longitude", "Latitude",
  # Id, area, and distance to the neighbour
  "id", "area", "dist",
  # Land use cover of the path
  "pri_path", "sec_path", "past_path", "crop_path", "urb_path",
  # Number of intervening water cells
  "n_water",
  # Land use of the site
  "pri_site", "sec_site", "past_site", "crop_site", "urb_site"
)
n_neighbors <- sapply(neighbours$nn, length)
n_rows <- sum(n_neighbors + (n_neighbors == 0))
res <- matrix(NA, nrow=n_rows, ncol=length(res_cols))
colnames(res) <- res_cols
# track which connecting line we're dealing with for each batch
i_line <- as.list(rep(1,n_batches))
names(i_line) <- letters[1:n_batches]
# track which row of the results matrix we're filling in
i_row <- 1
# Progress bar to track the sites
logger::log_info("Processing sites...")
conn_metrics <- {
  pb <- txtProgressBar(
    min = 0, max = length(sites), initial = 0, style = 3
  )
  for (i in 1:length(sites)) {
    setTxtProgressBar(pb, i)
    nns <- neighbours$nn[[i]]
    # To get the right path val, we need to find the batch
    site_batch <- batch_groups[i]
    if (length(nns) == 0) {
      # no neighbours
      res[i_row,1:3] <- c(i, as.numeric(site_locs[i,]))
      i_row <- i_row + 1
    } else {
      # some neighbours
      for (j in 1:length(nns)) {
        nn <- nns[j]
        res[i_row,] <- c(
          # Site
          i, as.numeric(site_locs[i,]),
          # Id, area, and distance to the neighbour
          nn, polys$area[nn], neighbours$dist[[i]][j],
          # Land use cover of the path and intervening water
          as.numeric(
            batched_vals[[site_batch]][i_line[[site_batch]],-1]
          ),
          # Land use of the site
          as.numeric(lu_sites[i,-1])
        )
        i_line[[site_batch]] <- i_line[[site_batch]] + 1
        i_row <- i_row + 1
      }
    }
  }
  close(pb)  # finish the progress bar
  as.data.frame(res)
}

# Store this.
saveRDS(
  conn_metrics,
  paste(
    "../data/europe_background_conn_",
    format(Sys.time(), "%Y-%m-%d"),
    ".rds",
    sep = ""
  )
)

