library(dplyr)

# Get the connectivity metrics
# file_suffix <- "_70"
# patch_cutoff <- 0.7
# conn <- as_tibble(readRDS(
#   "../data/europe_conn_2_70_r10000_k10_2022-06-22.rds"
# ))
# ... and for 90% cutoff
file_suffix <- "_90"
patch_cutoff <- 0.9
conn <- as_tibble(readRDS(
  "../data/europe_conn__eu_NA_90_r20000_k20_2022-07-28.rds"
))

# dp <- 3   # rounding of the long/lat for site matching
dp <- NA

# Get the abundance data for each site
abundance <- as_tibble(readRDS("../data/site_abund_europe_2022-07-14.rds"))

# Get the compositional dissimilarity data for the sites
comp_diss <- as_tibble(readRDS("../data/comp_diss_europe_2022-07-14.rds"))

######################################
## Process abundance and Bray index ##
######################################

if (!is.na(dp)) {
  comp_diss <- comp_diss %>%
    mutate(   # add the rounded coordinates
      Longitude = round(s2_long, digits = dp),
      Latitude = round(s2_lat, digits = dp)
    )
  abundance <- abundance %>%
    mutate(   # add the rounded coordinates
      Longitude = round(Longitude, digits = dp),
      Latitude = round(Latitude, digits = dp)
    )
  conn <- conn %>%
    mutate(   # add the rounded coordinates
      Longitude = round(Longitude, digits = dp),
      Latitude = round(Latitude, digits = dp)
    )
} else {
  comp_diss <- comp_diss %>%
    mutate(
      Longitude = s2_long,
      Latitude = s2_lat
    )
}


####################################
## Calculate connectivity metrics ##
####################################

# Convert area from sq. m to sq. km, and distance to km
conn$area <- conn$area / 1e6
conn$dist <- conn$dist / 1e3

# # Add the Ecoregion to the sites data.
# conn <- dplyr::left_join(
#   conn, distinct(abundance, Longitude, Latitude, Ecoregion),
#   by = c("Longitude", "Latitude")
# )

source("factored.R")
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
        exp(-dist[dist > 0]) * area[dist > 0]
      )
    ),
    ifm_w = ifelse(
      is.na(id),
      0,
      patch_cutoff * sum(
        exp(-weighted_dist[dist > 0]) * area[dist > 0]
      )
    ),
    ifm_w_05 = ifelse(
      is.na(id),
      0,
      patch_cutoff * sum(
        exp(-0.5 * weighted_dist[dist > 0]) * area[dist > 0]
      )
    )
  ) %>%
  ungroup() %>%
  unique()   # remove the duplicate rows so we have a single row per site

# Log the connectivity metrics. Add a small increment to avoid 0s.
conn2 <- mutate(conn2, across(d_nn_uw:ifm_w_05, \(x) log(x + 1e-5)))

# Finally, join the metrics to the data
abund_data <- abundance %>%
  dplyr::left_join(
    conn2,
    by = c("Longitude", "Latitude")
  ) %>%
  dplyr::filter(!is.na(site))  # remove the sites with no connectivity data


comp_diss_data <- comp_diss %>%
  dplyr::left_join(
    conn2,
    by = c("Longitude", "Latitude")
  ) %>%
  dplyr::filter(!is.na(site))  # remove the sites with no connectivity data

# Remove data with indeterminate land use (already gone from the comp diss
# data).
abund_data <- filter(
  abund_data,
  !is.na(LandUse)
)

# Count how many sites are not in a forested ecoregion
ecoregs <- table(abund_data$Ecoregion)
ecoregs[ecoregs != 0]

# Count how many patches had paths with intervening water
sum(ifelse(is.na(conn$n_water), FALSE, conn$n_water > 0)) / nrow(conn)
sum(ifelse(is.na(conn$n_water), FALSE, conn$n_water < 4)) / nrow(conn)
# hist(conn$n_water[conn$n_water > 0])

# Add a reduced habitat classification.
reduced_habitat <- function(predominant_habitat) {
  mod_hab <- factor(as.character(Map(
    function(x) {
      if (x == "Primary forest" || x == "Primary non-forest") {
        "Primary"
      } else if (x == "Young secondary vegetation") {
        "Young"
      } else if (x == "Intermediate secondary vegetation") {
        "Intermediate"
      } else if (x == "Mature secondary vegetation") {
        "Mature"
      } else if (x == "Secondary vegetation (indeterminate age)") {
        "Indeterminate age"
      } else if (x == "Primary minimal") {
        "Primary minimal"
      } else {
        "In use"
      }
    },
    predominant_habitat
  )))
  # make sure Primary minimal is the reference level.
  relevel(mod_hab, "Primary minimal")
}
abund_data$mod_hab <- reduced_habitat(abund_data$LandUse)
comp_diss_data$mod_hab <- reduced_habitat(comp_diss_data$s2_lu)


# Check for co-linearity in the explanatory variables. We see some here, since
# obviously the buffer measures are related to each other, as are the IFM
# measures. If we remove these kinds of duplicates though, the values are good.
source("extern/HighstatLibV10.R")
corvif(abund_data[c("LandUse", "d_nn_w")])
corvif(abund_data[c("LandUse", "b2_w")])
corvif(abund_data[c("LandUse", "b5_w")])
corvif(abund_data[c("LandUse", "b7_w")])
corvif(abund_data[c("LandUse", "b10_w")])
corvif(abund_data[c("LandUse", "ifm_w")])

corvif(comp_diss_data[c("geog_dist", "mod_hab", "d_nn_w")])
corvif(comp_diss_data[c("geog_dist", "mod_hab", "b2_w")])
corvif(comp_diss_data[c("geog_dist", "mod_hab", "b5_w")])
corvif(comp_diss_data[c("geog_dist", "mod_hab", "b7_w")])
corvif(comp_diss_data[c("geog_dist", "mod_hab", "b10_w")])
corvif(comp_diss_data[c("geog_dist", "mod_hab", "ifm_w")])

# This is interesting, because if we increase the area of our nearby patches by
# getting further away from a small patch that we were close too, this might be
# bad.
corvif(abund_data[c("LandUse", "d_nn_w", "b2_w")])

# But it looks like there's a negative association -- further from nearest
# patch, lower buffer measure.
m <- lm(d_nn_w ~ b2_w, abund_data)
summary(m)
# plot(abund_data$d_nn_w, abund_data$b2_w)


# Continue with the analysis
library(car)
comp_diss_data <- comp_diss_data %>%
  # Since the similarity index is bounded between 0 and 1, the researchers have
  # found that a logit transfomation goes some way towards normalising the
  # errors -- as well as being conceptually more appropriate than a log
  # transform since it recognises the boundedness of the data.
  mutate(
    logitCS = car::logit(bray, adjust=0.001, percents=F)
  ) %>%
  # Also log-transform the distance between sites.
  mutate(
    log10geo = log10(geog_dist + 1)  # add 1 to remove zeros
  ) %>%
  # Set primary-minimal as the baseline
  mutate(
    lu_contrast = relevel(
      factor(lu_contrast), ref='Primary minimal_vs_Primary minimal'
    )
  )

# Doing this stops us having rank deficiency. Before, the indeterminate age is
# only associated with single values for the connectivity.
table(comp_diss_data$mod_hab, round(comp_diss_data$ifm_w, 0))
comp_diss_data <- comp_diss_data[
  comp_diss_data$mod_hab != "Indeterminate age",
]

# Count how many sites we have now.
sites_info <- data.frame(
  name = c("abundance", "cs"),
  sites = c(
    length(unique(abund_data$SSBS)),
    length(unique(comp_diss_data$s2))
  ),
  studies = c(
    length(unique(abund_data$SS)),
    length(unique(comp_diss_data$SS))
  )
)
write.csv(
  sites_info,
  paste(
    "../results/sites", file_suffix, ".csv",
    sep = ""
  ),
  row.names = FALSE, na = "", quote = FALSE
)

# Make a nice illustration of the sites.
abund_sites_all <- dplyr::distinct(abund_data, Latitude, Longitude)
cs_sites_all <- dplyr::distinct(comp_diss_data, Latitude, Longitude)
both_sites <- dplyr::inner_join(abund_sites_all, cs_sites_all)
abund_sites <- setdiff(abund_sites_all, both_sites)
cs_sites <- setdiff(cs_sites_all, both_sites)

longmin <- -11
longmax <- 30
latmin <- 35
latmax <- 70
gis_data <- fetch_raster(latmin, latmax, longmin, longmax)

pdf(
  paste(
    "../results/europe-figure", file_suffix, ".pdf",
    sep = ""
  )
  # width = 10, height = 8 # inches
)
plot(gis_data$primary + gis_data$secondary, axes = FALSE)
points(abund_sites, col = "blue", pch = 16, cex = 0.7)
# points(cs_sites, col = "red", pch = 16)    # there are none
points(both_sites, col = "black", pch = 16, cex = 0.7)
dev.off()

# Add the abundance data to the cs data -- to compare AICs for abundance on the
# reduced dataset.
cs_dat_ab <- comp_diss_data %>%
  dplyr::left_join(
    abund_data, by = c("s2" = "SSBS")
  )


#####################
## Fit some models ##
#####################
# This function is Adriana's.
model_plot <-function(mod.for.plot){
  par(mfrow = c(1,3))
  qqnorm(resid(mod.for.plot))
  qqline(resid(mod.for.plot), col = 2)
  plot(
    fitted(mod.for.plot),
    resid(mod.for.plot),
    xlab = "Fitted Values", ylab = "Residuals",
    main = "Residuals vs fitted"
  )
  abline(h=0, lty=2)
  lines(smooth.spline(fitted(mod.for.plot), resid(mod.for.plot)), col = "red")
  hist(resid(mod.for.plot))
  par(ask = TRUE)    # let me scroll through
  lme4::dotplot.ranef.mer(ranef(mod.for.plot,condVar = TRUE))
}

# Abundance
get_abund_eqn <- function(metric, random_slope = NA) {
  paste(
    "sqrt(RescaledAbundance) ~ ",
      # "LandUse + ",
      # "mod_hab : ", metric,
      if (!is.na(metric)) paste("LandUse * ", metric, sep="") else "LandUse",
      " + (1|SS)",
      " + (1|SSB)",
      if (!is.na(random_slope)) paste(" + ", random_slope, sep=""),
    sep = ""
  )
}

fit_models <- function(data, eqn) {
  list(
    null = fit_glmm(
      data, eqn(NA)
    ),
    ifm_w = fit_glmm(
      data, eqn("ifm_w")
    ),
    ifm_uw = fit_glmm(
      data, eqn("ifm_uw")
    ),
    b2_w = fit_glmm(
      data, eqn("b2_w")
    ),
    b2_uw = fit_glmm(
      data, eqn("b2_uw"), optim = "Nelder_Mead"
    ),
    b5_w = fit_glmm(
      data, eqn("b5_w")
    ),
    b5_uw = fit_glmm(
      data, eqn("b5_uw")
    ),
    b7_w = fit_glmm(
      data, eqn("b7_w")
    ),
    b7_uw = fit_glmm(   # to help the CS model converge.
      data, eqn("b7_uw"), optim = "Nelder_Mead"
    ),
    b10_w = fit_glmm(
      data, eqn("b10_w")
    ),
    b10_uw = fit_glmm(
      data, eqn("b10_uw")
    ),
    d_nn_w = fit_glmm(
      data, eqn("d_nn_w")
    ),
    d_nn_uw = fit_glmm(
      data, eqn("d_nn_uw")
    )
  )
}

model_abund <- fit_models(abund_data, get_abund_eqn)

# See if mod_hab gives improvement over null -- turns out it doesn't. The gap is
# too big to be closed by removing this few parameters.
model_abund_mod_hab <- fit_models(
  abund_data, \(metric) {
    paste(
      "sqrt(RescaledAbundance) ~ ",
      if (!is.na(metric)) {
        paste("mod_hab * ", metric, sep="")
      } else {
        "mod_hab"
      },
      " + (1|SS)",
      " + (1|SSB)", sep=""
    )
  }
)

# Examine the results for data reduced to those where we have nearest-neighbour
# values.
model_abund_red <- fit_models(
  abund_data[!is.na(abund_data$d_nn_w),],
  get_abund_eqn
)

# Examine how AICs compare on the smaller dataset alone. TODO: make this work
model_abund_cs_data <- fit_models(
  cs_dat_ab, \(metric) {
    paste(
      "sqrt(RescaledAbundance) ~ ",
      if (!is.na(metric)) {
        paste("mod_hab.y * ", metric, ".y", sep="")
      } else {
        "mod_hab.y"
      },
      " + (1|SS.y)",
      " + (1|SSB)", sep=""
    )
  }
)

# Doesn't seem to be a problem with only 0/1:
table(round(cs_dat_ab$RescaledAbundance, 1), cs_dat_ab$mod_hab.y)


# Compositional dissimilarity
get_cs_eqn <- function(metric) {
  paste(
    "logitCS ~ ",
      # "s2_lu * ", metric,   # rank deficiency
      if (!is.na(metric)) paste("mod_hab * ", metric, sep="") else "mod_hab",
      " + log10geo",
      " + (1|SS)",
      " + (1|s2)",
    sep = ""
  )
}

model_cs <- fit_models(comp_diss_data, get_cs_eqn)

# all <- lme4::allFit(model_cs$ifm_w)
# ss <- summary(all)
# ss$which.OK

# Examine the results for data reduced to those where we have nearest-neighbour
# values.
model_cs_red <- fit_models(
  comp_diss_data[!is.na(comp_diss_data$d_nn_w),],
  get_cs_eqn
)

## AIC
save_aic_data <- function(models, file_suffix, inc_nn = FALSE) {
  aic <- data.frame(
    model = c(
      "Null",
      "IFM", NA,
      "2km buffer", NA,
      "5km buffer", NA,
      "7km buffer", NA,
      "10km buffer", NA,
      if(inc_nn) c("Nearest-neighbour", NA) else c()
    ),
    weighting = c("NA", rep(c("Weighted", "Unweighted"), if(inc_nn) 6 else 5)),
    aic = if(inc_nn) {
      AIC(
        models$null,
        models$ifm_w,
        models$ifm_uw,
        models$b2_w,
        models$b2_uw,
        models$b5_w,
        models$b5_uw,
        models$b7_w,
        models$b7_uw,
        models$b10_w,
        models$b10_uw,
        models$d_nn_w,
        models$d_nn_uw
      )$AIC
    } else {
      AIC(
        models$null,
        models$ifm_w,
        models$ifm_uw,
        models$b2_w,
        models$b2_uw,
        models$b5_w,
        models$b5_uw,
        models$b7_w,
        models$b7_uw,
        models$b10_w,
        models$b10_uw    # can't include, due to missing values
        # model_abund$dnn_w,
        # model_abund$dnn_uw
      )$AIC
    }
  )
  aic$aicDelta <- aic$aic - min(aic$aic)
  x <- exp(-0.5 * aic$aicDelta)
  aic$weight <- x / sum(x)
  aic$aicDelta <- round(aic$aicDelta, 3)
  aic$weight <- round(aic$weight, 3)
  write.csv(
    aic[c("model", "weighting", "aicDelta", "weight")],
    paste(
      "../results/aic", file_suffix, ".csv",
      sep = ""
    ),
    row.names = FALSE, na = "", quote = FALSE
  )
  return(aic)
}

abund_aic <- save_aic_data(model_abund, paste("_abund", file_suffix, sep=""))
abund_red_aic <- save_aic_data(
  model_abund_red,
  paste("_abund_red", file_suffix, sep=""),
  inc_nn = TRUE
)
cs_aic <- save_aic_data(model_cs, paste("_cs", file_suffix, sep=""))
cs_red_aic <- save_aic_data(
  model_cs_red,
  paste("_cs_red", file_suffix, sep=""),
  inc_nn = TRUE
)

### BEST MODELS
## using `emmeans` here to look at the contrasts in the interacting variables
library(emmeans)

best_ab <- model_abund$ifm_w
ab_metric <- "ifm_w"
best_cs <- model_cs$ifm_uw
cs_metric <- "ifm_uw"

# Slopes
emm <- emtrends(best_ab, "LandUse", var = ab_metric)
ss_abund <- summary(emm, infer = c(TRUE, TRUE), null = 0)
emm <- emtrends(best_cs, "mod_hab", var = cs_metric)
ss_cs <- summary(emm, infer = c(TRUE, TRUE), null = 0)

# # Intercepts
# emm <- emtrends(model_abund$b2_w, "LandUse", var = "b2_w")
# ss_abund <- summary(emm, infer = c(TRUE, TRUE), null = 0)
# emm <- emtrends(model_cs$ifm_uw, "mod_hab", var = "ifm_uw")
# ss_cs <- summary(emm, infer = c(TRUE, TRUE), null = 0)

write_model_results <- function(ss, name) {
  ss <- as.data.frame(ss)
  res_numeric <- ss[c(
    2,     # coefficient
    3,     # standard error
    5, 6,  # lower and upper confidence
    7, 8   # z-ratio and p-value
  )]
  colnames(res_numeric) <- c(
    "estimate", "stderr", "lower", "upper", "z", "p"
  )
  res <- cbind(ss[1], signif(res_numeric, 3))
  colnames(res)[1] <- "landuse"
  write.csv(
    res,
    paste(
      "../results/model_", name, file_suffix, ".csv",
      sep = ""
    ),
    row.names = FALSE, na = "", quote = FALSE
  )
  return(res)
}

ab_res <- write_model_results(ss_abund, "abund")
cs_res <- write_model_results(ss_cs, "cs")


#### ANOVA
library(car)

Anova(model_abund$ifm_w)

anova(model_abund$ifm_w, model_abund$null)