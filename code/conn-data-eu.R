source("factored.R")

# Don't round the coordinates of the sites.
dp <- NA

# Patch cutoff.
# With patch cutoff of 0.9, get 99289 polys.
# With patch cutoff of 0.7, get 81688 polys (they are more easily joined up).
patch_cutoff <- 0.9

# Coordinates of the edges of Europe.
longmin <- -11
longmax <- 30
latmin <- 35
latmax <- 70
site_latmin <- 37  # to avoid sites on the border

# Get the data and the rasters, and wrangle.
data <- fetch_data(site_latmin, latmax, longmin, longmax)
gis_data <- fetch_raster(latmin, latmax, longmin, longmax)
wrangled <- wrangle_data(data, dp = dp)

# Config properties of the metrics
k <- 20
radius <- 20000

# Set up file names for saving intermediate data
file_suffix <- paste(
  "_eu_",
  dp, "_",
  patch_cutoff * 100, "_",
  "r", radius, "_k", k, "_",
  # format(Sys.time(), "%Y-%m-%d"),
  "2022-07-28",
  sep = ""
)

# Extract the raw connectivity metrics -- this will take some time.
conn_metrics <- get_conn_metrics(
  gis_data, wrangled, patch_cutoff = patch_cutoff,
  radius = radius, k = k, file_suffix = file_suffix
)

# Finish it off with the saved data -- still don't know why this OOMs inside the
# function call above.
site_locs <- readRDS(paste("../data/site_locs", file_suffix, ".rds", sep = ""))
polys_area <- readRDS(
  paste("../data/polys_area", file_suffix, ".rds", sep = "")
)
neighbours <- readRDS(
  paste("../data/neighbours", file_suffix, ".rds", sep = "")
)
vals <- readRDS(paste("../data/vals", file_suffix, ".rds", sep = ""))
lu_sites <- readRDS(paste("../data/lu_sites", file_suffix, ".rds", sep = ""))

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
res <- matrix(NA, nrow = n_rows, ncol = length(res_cols))
colnames(res) <- res_cols
i_line <- 1   # track which connecting line we're dealing with
i_row <- 1    # track which row of the results matrix we're filling in
# Progress bar to track the sites
logger::log_info("Processing sites...")
for (i in 1:nrow(site_locs)) {
  nns <- neighbours$nn[[i]]
  if (length(nns) == 0) {
    # no neighbours
    res[i_row, 1:3] <- c(i, as.numeric(site_locs[i, ]))
    res[i_row, 13:17] <- as.numeric(lu_sites[i, -1])
    i_row <- i_row + 1
  } else {
    # some neighbours
    for (j in 1:length(nns)) {
      nn <- nns[j]
      res[i_row, ] <- c(
        # Site
        i, as.numeric(site_locs[i, ]),
        # Id, area, and distance to the neighbour
        nn, polys_area[nn], neighbours$dist[[i]][j],
        # Land use cover of the path and intervening water
        as.numeric(vals[i_line, -1]),
        # Land use of the site
        as.numeric(lu_sites[i, -1])
      )
      i_line <- i_line + 1
      i_row <- i_row + 1
    }
  }
}
conn_metrics <- as.data.frame(res)


# Store this.
saveRDS(
  conn_metrics,
  paste(
    "../data/europe_conn_",
    file_suffix,
    ".rds",
    sep = ""
  )
)
