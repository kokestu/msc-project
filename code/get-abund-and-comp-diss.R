# Adapted from https://adrianadepalma.github.io/BII_tutorial/bii_example.html
# for easy data manipulation
library(dplyr)
library(tidyr)
library(tibble)

# Running loops (and in parallel)
library(purrr)
library(furrr)

library(magrittr) # for piping
library(lme4) # for mixed effects models
library(car) # for logit transformation with adjustment
library(betapart) # for calculating balanced bray-curtis dissimilarity
library(raster) # for working with raster data
library(geosphere) # calculating geographic distance between sites
library(viridis) # figure colours for colour blindness

# Coordinates of the edges of Europe.
longmin <- -11
longmax <- 30
latmin <- 3
latmax <- 70

# Get the data for europe
source('factored.R')
diversity <- fetch_data(latmin, latmax, longmin, longmax)

# Let's look at the data -- how much do we have for different land uses and
# intensities? We see a pretty good balance. "minimal use primary vegetation" is
# the closest we have to a "pristine" baseline -- so it's often good to keep
# this separate for comparison.
table(diversity$Predominant_habitat, diversity$Use_intensity)

# Let's simplify the data to make this example easier.
diversity <- diversity %>%
  mutate(
    # Separate the "pristine" (ish) data
    LandUse = ifelse(
      grepl("primary", tolower(Predominant_habitat)) &
        Use_intensity == "Minimal use",
      "Primary minimal",
      as.character(Predominant_habitat)
    ),
    # Make all undecided cases into NA
    LandUse = ifelse(Predominant_habitat == "Cannot decide", NA, LandUse),
    # Relevel the factor so that Primary minimal is the first level -- this
    # makes it the intercept term in models.
    LandUse = relevel(factor(LandUse), ref='Primary minimal')
  )


#### Build effort-corrected measurements.
diversity <- diversity %>%
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

### TOTAL ABUNDANCE ###
abundance_data <- diversity %>%
  # Extract just the abundance measures
  filter(Diversity_metric_type == "Abundance") %>%
  # Group by SSBS -- so that each row is a unique site
  group_by(SSBS) %>%
  # Sum the abundance measurements per site
  mutate(TotalAbundance = sum(Effort_corrected_measurement)) %>%
  ungroup() %>%
  # Remove duplicate site measurements to simplify
  distinct(SSBS, .keep_all=T) %>%
  # Now find the maximum abundance per study
  group_by(SS) %>%
  mutate(MaxAbundance = max(TotalAbundance)) %>%
  ungroup() %>%
  # Now rescale values within each study
  mutate(RescaledAbundance = TotalAbundance/MaxAbundance)
# Store this.
saveRDS(
  abundance_data,
  paste(
    "../data/site_abund_europe_",
    format(Sys.time(), "%Y-%m-%d"),
    ".rds",
    sep = ""
  )
)

# Get the information about the taxa in each site
site_data <- diversity %>%
  filter(Class != "") %>%
  group_by(SSBS) %>%   # for each site
  summarise({
    row <- as.numeric(levels(Class) %in% Class)
    names(row) <- levels(Class)
    as.data.frame(t(row))
  })
saveRDS(
  site_data[-2],
  paste(
    "../data/site_taxa_europe_",
    format(Sys.time(), "%Y-%m-%d"),
    ".rds",
    sep = ""
  )
)

# An important thing to notice here is that we use the "effort corrected
# measurement". This is because in some studies the "sampling effort" varies
# between sites -- e.g. the number of samples taken might vary per site. For
# abundance this is okay, since we know that abundance scales linearly with
# sampling effort. However, for species richness, this doesn't hold so below
# we'll use the raw measurements, and only include studies with invariable
# sampling effort.

### COMPOSITIONAL SIMILARITY ###
# Originally, the PREDICTS project used the "asymmetric Jacard Index" but this
# has changed recently to the "balanced Bray-Curtis". This is because the Jacard
# Index doesn't take into account changes in community structure, only in the
# species present. The Bray-Curtis (on the other hand) will show increased
# dissimilarity for sites where dominance is different, even when the species
# pools are identical.

# First, let's filter the data to make sure everything is relevant
bray_curtis_input <- diversity %>%
  # Remove unknown land use
  filter(!is.na(LandUse)) %>%
  # Only look at abundance data
  filter(Diversity_metric_type =="Abundance") %>%
  # Calculate some aggregate stuff per study
  group_by(SS) %>%
    # Number of unique sampling efforts
    mutate(n_sample_effort = n_distinct(Sampling_effort)) %>%
    # Number of unique species
    mutate(n_species = n_distinct(Taxon_name_entered)) %>%
    # Number of primary minimal sites
    mutate(n_primin_records = sum(LandUse == "Primary minimal")) %>%
  ungroup() %>%
  # Now we filter:
  filter(
    n_sample_effort == 1 &   # invariable sampling effort
    n_species > 1 &  # otherwise likely not looking at community-level diversity
    n_primin_records > 0  # at least some primary minimal data
  ) %>%
  # Get rid of empty factor levels
  droplevels()

# Function for calculating the balanced Bray-Curtis compositional similarity
# between a single pair of sites.
get_bray_curtis <- function(s1, s2, data) {
  require(betapart)
  sp_data <- data %>%
    # Extract the data for these sites
    filter(SSBS %in% c(s1, s2)) %>%
    # Get the columns we're interested in
    dplyr::select(
      SSBS, Taxon_name_entered, Measurement
    ) %>%
    # Pivot so that the columns are species
    pivot_wider(
      names_from = Taxon_name_entered,
      values_from = Measurement
    ) %>%
    # Make the SSBS the rownames
    column_to_rownames("SSBS")
  # Calculate the similarity
  n_empty_sites <- sum(rowSums(sp_data) == 0, na.rm=T)
  if (n_empty_sites == 1) {
    # One site empty, we have similarity of 0
    bray <- 0
  } else if (n_empty_sites == 2) {
    # Both sites are empty, have NA similarity
    bray <- NA
  } else {
    # Otherwise, calculate the similarity as 1 - bray
    bray <- 1 - pluck(bray.part(sp_data), 'bray.bal', 1)
  }
  return(bray)
}

# To calculate the index, first let's get all of the study IDs that we want to
# compare.
studies <- bray_curtis_input %>%
  distinct(SS) %>%
  pull()   # to get the result as a vector

# Now use purrr to do a loop to get the pairings we need for calculating the
# similarity. Using map_dfr so that we get a dataframe at the end.
site_comparisons <- map_dfr(.x=studies, .f=function(x) {
  # Get the study we're interested in
  site_data <- bray_curtis_input[bray_curtis_input$SS == x,] %>%
    # Extract SSBS and LandUse information
    dplyr::select(SSBS, LandUse) %>%
    # Remove duplicate rows per site, to simplify the data
    distinct(SSBS, .keep_all=T)
  # Extract the primary minimal sites: we're only interested in comparisons with
  # this baseline.
  baseline_sites <- site_data %>%
    filter(LandUse == "Primary minimal") %>%
    pull(SSBS)
  # Get all baseline site x site comparisons
  site_comparisons <- expand.grid(baseline_sites, pull(site_data, SSBS)) %>%
    # Rename the columns for our pairwise function
    rename(s1=Var1, s2=Var2) %>%
    # Remove the self-comparisons
    filter(s1 != s2) %>%
    mutate(
      # Factor -> string
      s1 = as.character(s1),
      s2 = as.character(s2),
      # Add the full name
      contrast = paste(s1, 'vs', s2, sep='_'),
      # And the study id
      SS = as.character(x)
    )
  return(site_comparisons)
})


####### BEWARE, this takes ages #######
# Set up the `future` plan -- used for parallel computation.
future::plan(
  multisession,   # do the computation asynchronously on separate R sessions
  workers = 3     # avoid killing my 4-core machine
)
# Do the map -- use map2_dbl for a two parameter function that returns a double.
bray <- future_map2_dbl(
  .x = site_comparisons$s1,
  .y = site_comparisons$s2,
  ~get_bray_curtis(s1=.x, s2=.y, data=bray_curtis_input)
)
# Revert to sequential execution
future::plan(sequential)
# Store this.
saveRDS(
  bray,
  paste(
    "../data/bray_europe_",
    format(Sys.time(), "%Y-%m-%d"),
    ".rds",
    sep = ""
  )
)
#####################################

# Now, we'll pick up the other data we need that doesn't need a loop.
latlongs <- bray_curtis_input %>%
  group_by(SSBS) %>%   # for each site
  summarise(
    Lat = unique(Latitude),
    Long = unique(Longitude)
  )

# Get the land uses
lus <- bray_curtis_input %>%
  group_by(SSBS) %>%
  summarise(lu = unique(LandUse))

# Put all the data together
cd_data <- site_comparisons %>%
  # Bray-Curtis data is already in the right order
  mutate(bray=bray) %>%
  # Latitude and longitude for both sites
  left_join(latlongs, by=c("s1"="SSBS")) %>%
  rename(s1_lat = Lat, s1_long = Long) %>%
  left_join(latlongs, by=c("s2"="SSBS")) %>%
  rename(s2_lat = Lat, s2_long = Long) %>%
  # Geographic distance between them
  mutate(geog_dist=distHaversine(  # default unit is metres
    cbind(s1_long, s1_lat), cbind(s2_long, s2_lat)
  )) %>%
  # Land uses for both sites
  left_join(lus, by=c("s1"="SSBS")) %>%
  rename(s1_lu = lu) %>%
  left_join(lus, by=c("s2"="SSBS")) %>%
  rename(s2_lu = lu) %>%
  # Contrast column -- for the modelling below
  mutate(lu_contrast=paste(s1_lu, s2_lu, sep='_vs_'))


# Store this.
saveRDS(
  cd_data,
  paste(
    "../data/comp_diss_europe_",
    format(Sys.time(), "%Y-%m-%d"),
    ".rds",
    sep = ""
  )
)
