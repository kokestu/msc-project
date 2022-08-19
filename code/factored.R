## Library functions for all of the end-to-end analysis of the diversity data.
## It should be possible to call the top function on some spatial extent +
## objects defining metrics to include, and produce images and model output.
library(logger)  # for logging stuff (duh)

# Import the required libraries.
library(dplyr)
library(tidyr)
library(terra)
library(nngeo)
library(lwgeom)
library(sf)
library(units)
library(gamm4)
library(lme4)
library(blme)   # for fitting me models with bayesian methods

# Disable the use of s2 for coordinate geometry. Otherwise, the vector
# geometries made by terra break things (annoyingly).
sf::sf_use_s2(FALSE)

# Choose the number of cores to use for parallel computation.
n_cores <- 1

# Fetch the PREDICTS data for the defined range.
fetch_data <- function(latmin, latmax, longmin, longmax) {
  logger::log_info("Loading PREDICTS data")
  diversity <- readRDS("../data/diversity-2022-04-05-02-33-06.rds")
  # Crop to the appropriate range.
  diversity[
    diversity$Latitude < latmax & diversity$Latitude > latmin &
      diversity$Longitude < longmax & diversity$Longitude > longmin,
  ]
}


# Fetch the land use rasters for the defined range.
fetch_raster <- function(latmin, latmax, longmin, longmax) {
  logger::log_info("Loading rasters")
  # Load rasters
  gis_data <- c(
    # We don't need the ice data
    # "../data/landuse/ICE/ICE_1km_2005.bil",            # ice cover
    # primary vegetation
    terra::rast("../data/landuse/PRI_2005/PRI_1km_2005_0ice.bil"),
    # secondary vegetation
    terra::rast("../data/landuse/SEC_2005/SEC_1km_2005_0ice.bil"),
    # pasture
    terra::rast("../data/landuse/PAS_2005/PAS_1km_2005_0ice.bil"),
    # cropland
    terra::rast("../data/landuse/CRP_2005/CRP_1km_2005_0ice.bil"),
    # urban cover
    terra::rast("../data/landuse/URB_2005/URB_1km_2005_0ice.bil")
  )     # applies WGS84 by itself
  names(gis_data) <- c(
    "primary", # layer of primary vegetation
    "secondary", # layer of secondary vegetation
    "pasture", # layer of pasture
    "cropland", # layer of cropland
    "urban" # layer of urban
  )
  # Crop to the appropriate range
  so_extent <- terra::ext(longmin, longmax, latmin, latmax)
  terra::crop(gis_data, so_extent)
}


# Wrangle the PREDICTS data to:
#   - Mark all secondary habitat in the same category: we will use the age
#     values for secondary vegetation in the model.
#   - Add effort-corrected diversity metrics
#   - Remove duplicate sites data
#   - Sqrt effort-corrected measurement to make less non-normal
#   - Round the Longitude/Latitude to a similar resolution to the raster.
wrangle_data <- function(data, dp=NA) {
  logger::log_info("Wrangling data")
  # Mark all secondary vegetation as the same category.
  all_secondary <- grep(
    "secondary", data$Predominant_habitat, ignore.case = TRUE
  )
  data$land_use <- as.character(data$Predominant_habitat)
  data[all_secondary, ]$land_use <- "Secondary vegetation"
  data$land_use <- factor(data$land_use)
  # Build effort-corrected measurements.
  data <- data %>%
    dplyr::mutate(x = tidyr::replace_na(Sampling_effort, 1)) %>%
    dplyr::mutate(
      Corrected_sampling_effort = ifelse(
        Diversity_metric_is_effort_sensitive,
        1, Sampling_effort
      ),
      Effort_corrected_measurement = ifelse(
        Diversity_metric_is_effort_sensitive,
        Measurement / Sampling_effort, Measurement
      )
    )
  # Now we want to merge data from repeated measurements of the same site.
  data <- data %>%
    # Group by aspects of the sites that should be identical if we need to merge
    # the abundances.
    # I only want to merge abundances if they are within the same study and
    # block, as I'm assuming that even if the locations and sampling times are
    # the same, if the blocks or studies are different, then there is some good
    # reason for this.
    dplyr::group_by(
      Source_ID, Study_number, Study_name, Block,
      # diversity metric type
      Diversity_metric, Diversity_metric_type, Diversity_metric_unit,
      Diversity_metric_is_effort_sensitive,
      # details of the sites
      Predominant_habitat, Use_intensity,
      Years_since_fragmentation_or_conversion,
      # details of the sampling method
      Sampling_method, Sampling_effort_unit,
      # species identity
      Study_common_taxon, Rank_of_study_common_taxon,
      Taxon_number, Taxon_name_entered,
      Indication, Parsed_name,
      Best_guess_binomial, COL_ID, Taxon, Name_status,
      Rank, Kingdom, Phylum, Class, Order, Family, Genus, Species,
      Higher_taxon,
      # site location
      Longitude, Latitude,
      # sampling time
      Sample_start_earliest, Sample_end_latest, Sample_date_resolution
    ) %>%
    # If the diversity metric is occurrence then if it is present at all, give
    # it a 1, if it is always absent, give it a 0.
    # Otherwise (if the metric is either abundance or species richness),
    # calculate the weighted abundance/richness for each taxonomic group,
    # weighted by sampling effort.
    dplyr::mutate(
      merged_diversity = ifelse(
        Diversity_metric_type == "Occurrence",
        # If any of the occurrence values are 1, `any` will return
        # TRUE.
        as.integer(any(Effort_corrected_measurement > 0)),
        # note that since we've already corrected the sampling effort, this is
        # essentially a mean rather than a weighted mean for abundance
        # measurements. It's a weighted mean for species richness though where
        # sampling efforts vary.
        stats::weighted.mean(
          x = Effort_corrected_measurement,
          w = Corrected_sampling_effort
        )
      )
    )
  # Extract the group data so that we can create metrics about the merging.
  group_dat <- data %>%
    dplyr::group_data() %>%
    dplyr::mutate(
      nvals_merged = lengths(.rows),
      merge_id = dplyr::row_number()
    )
  # And ungroup the data now that we have the group summaries.
  data <- dplyr::ungroup(data)
  # Include the grouping data in the main merged table.
  data <- dplyr::left_join(
    data, group_dat,
    by = c(
      "Source_ID", "Study_number", "Study_name", "Diversity_metric",
      "Diversity_metric_unit", "Diversity_metric_type",
      "Diversity_metric_is_effort_sensitive", "Sampling_method",
      "Sampling_effort_unit", "Block", "Sample_start_earliest",
      "Sample_end_latest", "Sample_date_resolution", "Predominant_habitat",
      "Use_intensity", "Years_since_fragmentation_or_conversion",
      "Longitude", "Latitude", "Taxon_number", "Taxon_name_entered",
      "Indication", "Parsed_name", "COL_ID", "Taxon", "Name_status",
      "Rank", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
      "Species", "Higher_taxon", "Study_common_taxon",
      "Rank_of_study_common_taxon", "Best_guess_binomial"
    )
  )
  # Round the latitude and longitude to approx. 1km
  # 1 deg ~ 110km
  if(!is.na(dp)) {
    data$Latitude <- round(data$Latitude, dp)
    data$Longitude <- round(data$Longitude, dp)
  }
  # Square-root transform the abundance data to try and eliminate issues with
  # non-normal error structure. Also scale per-study so that the max.
  # abundance is 1, before sqrt.
  data %>%
    group_by(SSBS, Class) %>%
    dplyr::mutate(
      class_abund = sum(Effort_corrected_measurement),
    ) %>%
    ungroup() %>%
    distinct(SSBS, Class, .keep_all = TRUE) %>%
    group_by(SS, Class) %>%
    dplyr::mutate(
      study_max = max(class_abund)
    ) %>%
    ungroup() %>%
    dplyr::mutate(
      scaled_abund = ifelse(
        study_max == 0,
        0,
        class_abund / study_max
      ),
      sqrt_abund = sqrt(scaled_abund)
    )
}


# Calculate the metrics provided in the list and add the results as columns to
# the dataframe.
calculate_metrics <- function(
  wrangled,   # wrangled data to extract context features for
  raster,     # raster containing land use information
  metrics,    # list of functions that generate metrics
  cex,        # size of the points in output plots
  file_suffix # suffix to differentiate file output, file_suffix)
) {
  # Get hold of the unique sites and build a dataframe.
  coords <- data.frame(
    strsplit(unique(paste(
      wrangled$Longitude, wrangled$Latitude, wrangled$Ecoregion,
      sep = "#"
    )), split = "#")
  )
  sites <- data.frame(
    Longitude = as.numeric(coords[1, ]),
    Latitude = as.numeric(coords[2, ]),
    Ecoregion = as.character(coords[3, ])
  )
  # Transform sites into points object.
  site_points <- terra::vect(
    sites, # source data frame
    geom = c("Longitude", "Latitude"), # columns for coordinates
    crs = "WGS84"
  )
  logger::log_info(
    "Found ", as.character(length(site_points)), " sites for analysis"
  )
  # Extract metrics and add them to the dataframe.
  for (metric in metrics) {
    logger::log_info(
      "Extracting metrics for ", paste(metric$outputs, collapse = ",")
    )
    res <- metric$f(site_points, raster, cex, file_suffix)
    sites <- left_join(
      sites, res,
      by = c("Longitude", "Latitude")
    )
  }
  # Return the final results of all the metrics joined to the original data.
  left_join(
    wrangled, sites,
    by = c("Longitude", "Latitude", "Ecoregion")
  )
}


# Fetch connectivity-relevant metrics from the GIS and site data. This is an
# alternative to the separate metrics approached I used above, and should get us
# the data we need to calculate all the relevant metrics of connectivity, that
# is, for each site-neighbour pair:
#   - Area of the neighbour
#   - direct distance to the neighbour (unweighted)
#   - proportion of path per land use type
#   - land use cover of the site, from the raster
get_conn_metrics <- function(
  gis_data, wrangled, patch_cutoff = 0.7,
  radius = 10000, k = 10, file_suffix = ""
) {
  # Get hold of the unique sites and build a dataframe.
  logger::log_info("Building sites...")
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
  logger::log_info(
    "Found ", as.character(length(sites)), " sites for analysis"
  )
  # Get the patch polygons.
  logger::log_info("Extracting polygons...")
  polys <- get_patch_polys(gis_data, patch_cutoff)
  # Get the land-use cover in the cell of the site itself.
  logger::log_info("Collecting land cover data at sites...")
  lu_sites <- terra::extract(
    gis_data, sites
  )
  # Run garbage collection to avoid issues.
  logger::log_info("Running garbage collection")
  gc()
  # Calculate the distances to nearest habitat patch.
  logger::log_info("Finding neighbours...")
  neighbours <- nngeo::st_nn(
    st_as_sf(sites), st_as_sf(polys),
    k = min(k, length(polys)),
    maxdist = radius,
    parallel = n_cores,
    returnDist = T
  )
  # Calculate the lines between each site and its neighbor.
  lines_between <- as(nngeo::st_connect(
    st_as_sf(sites), st_as_sf(polys),
    ids = neighbours[["nn"]]
  ), "SpatVector")
  # Extract the land-use values under the line, and summarise these.
  logger::log_info("Extracting path data...")
  vals <- terra::extract(
    gis_data, lines_between,
    along = TRUE
  )
  vals <- vals %>%
    # Make sure we keep track of how many patches of water we cover.
    mutate(
      water = is.nan(primary),
      across(, function(x) ifelse(is.nan(x), 0, x))
    ) %>%
    group_by(ID) %>%
    # And sum the proportions of each land use along the lines.
    summarise(across(, sum)) %>%
    ungroup()
  # Save progress in case of failure
  logger::log_info("Saving data...")
  saveRDS(
    site_locs, file = paste(
      "../data/site_locs", file_suffix, ".rds",
      sep = ""
    )
  )
  saveRDS(
    polys$area, file = paste(
      "../data/polys_area", file_suffix, ".rds",
      sep = ""
    )
  )
  saveRDS(
    neighbours, file = paste(
      "../data/neighbours", file_suffix, ".rds",
      sep = ""
    )
  )
  saveRDS(
    vals, file = paste(
      "../data/vals", file_suffix, ".rds",
      sep = ""
    )
  )
  saveRDS(
    lu_sites, file = paste(
      "../data/lu_sites", file_suffix, ".rds",
      sep = ""
    )
  )
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
  i_line <- 1   # track which connecting line we're dealing with
  i_row <- 1    # track which row of the results matrix we're filling in
  # Progress bar to track the sites
  logger::log_info("Processing sites...")
  pb <- txtProgressBar(
    min = 0, max = length(sites), initial = 0, style = 3
  )
  for (i in 1:length(sites)) {
    setTxtProgressBar(pb, i)
    nns <- neighbours$nn[[i]]
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
          as.numeric(vals[i_line,-1]),
          # Land use of the site
          as.numeric(lu_sites[i,-1])
        )
        i_line <- i_line + 1
        i_row <- i_row + 1
      }
    }
  }
  close(pb)  # finish the progress bar
  # Return the result
  as.data.frame(res)
}



# TODO: Other metrics that can be interesting:
#   - Habitat coverage of adjacent cells

# A metric that counts the number of habitat patches within a given search
# radius (default 10km). Might be possible to build this out into a proper
# habitat network and do some analyses.
near_patches_metric <- function(radius) {
  list(
    outputs = paste(
      c("n_nearby_patches"),
      radius,
      sep=""
    ),
    f = function(sites, gis_data, cex, file_suffix) {
      # Get the patches
      polys <- get_patch_polys(gis_data)
      # Find the nearest patches
      nns <- nngeo::st_nn(
        st_as_sf(sites), st_as_sf(polys),
        k = length(polys),
        maxdist = radius,  # max search radius
        parallel = n_cores,
      )
      # Get number of nearby patches
      sites$n_nns <- sapply(nns, length)
      # Plot the sites by their metric value.
      pdf(
        paste(
          "../results/sites_with_nns", radius,
          file_suffix, ".pdf",
          sep = ""
        ),
        width = 10, height = 8 # inches
      )
      palette(terrain.colors(15))
      plot(gis_data$primary, legend = F)
      palette(hcl.colors(15))
      plot(
        sites, "n_nns",
        add = T,
        cex = cex, pch = 16, legend = T
      )
      palette("default")
      dev.off()
      # Return results
      coords <- crds(sites)
      res <- data.frame(
        Longitude = coords[, 1],
        Latitude = coords[, 2]
      )
      res[,paste("n_nearby_patches",radius,sep="")] <- sites$n_nns
      res
    }
  )
}

# A metric that produces a weighted path length to the nearest habitat patch,
# according to the BII of the intervening habitat.
w_path_metric <- list(
  outputs = c("weighted_path"),
  f = function(sites, gis_data, cex = 2, file_suffix = "") {
    # Get the patches
    polys <- get_patch_polys(gis_data)
    # Extract landscape context data for the sites only.
    logger::log_info("  Extracting spatial data")
    res_sites <- get_distance_data(sites, polys, gis_data)
    # Calculate the metric
    logger::log_info("  Calculating metric")
    res_sites$points$index <- mapply(
      calculate_index,
      split(res_sites$vals, res_sites$vals$ID),
      res_sites$points$Ecoregion
    )
    res_sites$points$weighted_path <-
      res_sites$points$index * res_sites$points$mindist
    # Plot the sites by their connectivity metric value.
    pdf(
      paste(
        "../results/sites_with_context", file_suffix, ".pdf",
        sep = ""
      ),
      width = 10, height = 8 # inches
    )
    palette(terrain.colors(15))
    plot(gis_data$primary, legend = F)
    palette("default")
    plot(polys, add = T, col = "lightblue")
    palette(hcl.colors(15))
    plot(
      res_sites$points, "weighted_path",
      add = T,
      cex = cex, pch = 16, legend = T
    )
    palette("default")
    plot(res_sites$edges, add = T)
    dev.off()
    # Return a dataframe as result.
    coords <- crds(res_sites$points)
    data.frame(
      Longitude = coords[, 1],
      Latitude = coords[, 2],
      weighted_path = res_sites$points$weighted_path,
      mindist = res_sites$points$mindist
    )
  }
)


# Function to build polygons of habitat cover from the provided raster.
# Patch cutoff is the proportion of habitat cover to consider it a "patch". In
# the experiment (connectivity-analysis), it's at approximately this point
# where we get convergence to a single patch.
get_patch_polys <- function(gis_data, patch_cutoff=0.7) {
  # Merge the primary and secondary vegetation rasters to get habitat coverage
  merged <- sum(gis_data[[c("primary", "secondary")]])
  # Potential alternative idea: maybe we can use BII to quantify something
  # like habitat quality, instead of just coverage. Then we can use a cutoff
  # in that for patches instead.
  # Classify according to habitat cover -- we need to do this to turn these
  # into vector features.
  classified <- merged > patch_cutoff
  classified[!classified] <- NA
  # Run garbage collection to avoid issues.
  logger::log_info("Running garbage collection")
  gc()
  # Create polygons out of the classified data
  polys <- terra::as.polygons(classified)
  # need to turn the single feature into many
  polys <- terra::disagg(polys)
  # Add some information about the patches to them and return.
  polys$area <- terra::expanse(polys)
  polys$perimeter <- terra::perim(polys)
  # # Sadly, this OOMs, so we'll have to just use 0.7 instead :(
  # polys$coverage <- terra::extract(merged, polys, fun=mean)
  polys
}


# Function that will extract connectivity features of points to polygons
# (representing habitat patches).
# Expects:
#     points: A set of points whose connectivity we're intersted in.
#     polys: A set of polyons representing habitat patches
#     raster: A raster containing data we want to extract along the paths.
# Returns:
#     points: The original points object with a new feature for the distances
#       from each point to the nearest patch.
#     vals: The values along the paths to the nearest patch, extracted from
#       the raster.
#     edges: The vector geometries of the connecting paths.
#     n_water_cells: The number of cells of water each path passes through.
get_distance_data <- function(sites, polys, gis_data) {
  # Run garbage collection to avoid issues.
  logger::log_info("Running garbage collection")
  gc()
  # Calculate the distances to nearest habitat patch.
  mindist <- nngeo::st_nn(
    st_as_sf(sites), st_as_sf(polys),
    k = 1,
    maxdist = Inf,    # unlimited search radius
    parallel = n_cores, returnDist = T,
  )
  sites$mindist <- as.numeric(mindist$dist)
  sites$nearest <- as.numeric(mindist[["nn"]])
  # Don't use terra for this atm.
  # mindist <- terra::nearby(sites, polys, centroids = FALSE)
  # sites$nearest <- mindist$k1
  #
  # Instead of just getting the distances between points, we can get the lines
  # that connect them too, with `st_nearest_points(poly1, poly2)`. It returns
  # a linestring with two points (the endpoints).
  # Then we can extract land use percentages along this line and weight the
  # result according to some metric of penetrability.
  #
  # TODO: This only does direct nearest neighbour, which should be okay at
  # this scale, I think. We need to use geosphere::greatCircle to do distances
  # properly, but I suspect this will end up being kind of tough...
  #
  # Run garbage collection to avoid issues.
  logger::log_info("Running garbage collection")
  gc()
  lines_between <- as(nngeo::st_connect(
    st_as_sf(sites), st_as_sf(polys),
    ids = mindist[["nn"]]
  ), "SpatVector")
  # Don't use terra for this (even though it's faster) since it doesn't notice
  # when sites are already within a polygon.
  # lines_between <- terra::nearest(
  #   sites, polys,
  #   centroids = FALSE,  # don't use the patch centroids, just the edges
  #   lines = TRUE   # fetch the lines themselves
  # )
  # Extract values along the line. This gets us all land use types -- probably
  # we will need to aggregate this a bit to make sure that it doesn't take
  # forever.
  #
  # This data can be used to weight the distance according to penetrability --
  # and we also have water marked as NA in the results too.
  vals <- terra::extract(
    gis_data, lines_between,
    along = TRUE
  )
  # Count instances of intervening water and report.
  n_water_cells <- vals %>%
    group_by(ID) %>%
    summarise(sum(is.nan(primary))) %>%
    `[[`(2)
  # Return this data
  list(
    points = sites,
    vals = vals,
    edges = lines_between,
    n_water_cells = n_water_cells
  )
}

# Set up a function for calculating a weighted path length for a point, based on
# the path to the nearest patch. NOTE: These paths are straight lines, not
# land-based paths. For now, we will treat intervening water as 100% urban.
# Use the ecoregion to fetch BII for forested or non-forested.
calculate_index <- function(path_data, ecoregion) {
  # Handle points very close to patches
  if (nrow(path_data) == 0) {
    return(0)
  }
  # Extract main types of landscape
  prim <- path_data[, "primary"]
  sec <- path_data[, "secondary"]
  past <- path_data[, "pasture"]
  crop <- path_data[, "cropland"]
  urb <- path_data[, "urban"]
  # Treat water as 100% urban for now.
  prim[is.na(prim)] <- 0
  sec[is.na(sec)] <- 0
  past[is.na(past)] <- 0
  crop[is.na(crop)] <- 0
  urb[is.na(urb)] <- 1
  # Weight according to BII, depending on whether or not the background
  # potential ecosystem is forested or non-forested. Data from:
  # https://www.nature.com/articles/s41586-020-2705-y/figures/3
  prim_bii <- 1 # primary
  mat_sec_f <- 1
  mat_sec_nf <- 0.775 # mature secondary
  yng_sec_f <- 0.875
  yng_sec_nf <- 0.775 # young secondary
  mat_sec_mu_f <- 1
  mat_sec_mu_nf <- 0.7 # mature secondary min. use
  past_f <- 0.7
  past_nf <- 0.37 # managed pasture
  pcrop_f <- 0.525
  pcrop_nf <- 0.7 # perennial cropland
  acrop_f <- 0.5
  acrop_nf <- 0.65 # annual cropland
  urb_f <- 0.55
  urb_nf <- 0.7 # urban
  # Find the appropriate value given the ecoregion.
  is_forest <- grepl("woodland|forest", ecoregion, ignore.case = T)
  sec_bii <- if (is_forest) {
    mean(c(mat_sec_f, yng_sec_f, mat_sec_mu_f))
  } else {
    mean(c(mat_sec_nf, yng_sec_nf, mat_sec_mu_nf))
  }
  crop_bii <- if (is_forest) {
    mean(c(pcrop_f, acrop_f))
  } else {
    mean(c(pcrop_nf, acrop_nf))
  }
  past_bii <- if (is_forest) past_f else past_nf
  urb_bii <- if (is_forest) urb_f else urb_nf
  # Weight the path according to the length travelled through each land use
  # type.
  return(
    sum(
      prim / prim_bii + sec / sec_bii + past / past_bii +
        crop / crop_bii + urb / urb_bii
    ) / nrow(path_data)
  )
}

# Calculate the BII penetrability index for the new data format.
calculate_index2 <- function(path_data, ecoregion) {
  # Extract main types of landscape
  prim <- path_data$pri_path
  sec <- path_data$sec_path
  past <- path_data$past_path
  crop <- path_data$crop_path
  # Treat water as 100% urban for now.
  urb <- path_data$urb_path + path_data$n_water
  # Weight according to BII, depending on whether or not the background
  # potential ecosystem is forested or non-forested. Data from:
  # https://www.nature.com/articles/s41586-020-2705-y/figures/3
  prim_bii <- 1 # primary
  mat_sec_f <- 1
  mat_sec_nf <- 0.775 # mature secondary
  yng_sec_f <- 0.875
  yng_sec_nf <- 0.775 # young secondary
  mat_sec_mu_f <- 1
  mat_sec_mu_nf <- 0.7 # mature secondary min. use
  past_f <- 0.7
  past_nf <- 0.37 # managed pasture
  pcrop_f <- 0.525
  pcrop_nf <- 0.7 # perennial cropland
  acrop_f <- 0.5
  acrop_nf <- 0.65 # annual cropland
  urb_f <- 0.55
  urb_nf <- 0.7 # urban
  # Find the appropriate value given the ecoregion.
  is_forest <- grepl("woodland|forest|taiga", ecoregion, ignore.case = TRUE)
  sec_bii <- ifelse(
    is_forest,
    mean(c(mat_sec_f, yng_sec_f, mat_sec_mu_f)),
    mean(c(mat_sec_nf, yng_sec_nf, mat_sec_mu_nf))
  )
  crop_bii <- ifelse(
    is_forest,
    mean(c(pcrop_f, acrop_f)),
    mean(c(pcrop_nf, acrop_nf))
  )
  past_bii <- ifelse(is_forest, past_f, past_nf)
  urb_bii <- ifelse(is_forest, urb_f, urb_nf)
  # Weight the path according to the length travelled through each land use
  # type.
  return(
    (prim / prim_bii + sec / sec_bii + past / past_bii +
      crop / crop_bii + urb / urb_bii)
    / (prim + sec + past + crop + urb)
  )
}

# Fit a GAMM and return the model object.
fit_gamm <- function(with_metrics, eqn) {
  logger::log_info("Fitting the GAMM")
  # Fit the GAMM.
  gamm4(
    formula(eqn), random = ~ (1 | SSB) + (1 | SS),
    data = with_metrics
  )
}

# Fit a GLMM and return the model output.
fit_glmm <- function(with_metrics, eqn, REML, optim = NA) {
  logger::log_info("Fitting the GLMM")
  if (is.na(optim)) optim <- "bobyqa"
  lme4::lmer(
    formula(eqn),
    data = with_metrics,
    REML = REML,   # whether to adjust fitting procedure with REML
    control = lme4::lmerControl(
      optimizer = optim,
      optCtrl = list(maxfun = 9999999)  # To give time for convergence
    )
  )
}

# Fit a mixed model with bayesian methods.
fit_blme <- function(with_metrics, eqn) {
  logger::log_info("Fitting the GLMM")
  blme::blmer(
    formula(eqn),
    data = with_metrics
  )
}

