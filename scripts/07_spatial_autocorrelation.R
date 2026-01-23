#!/usr/bin/env Rscript
# ============================================================================
# 07_spatial_autocorrelation.R - Spatial Autocorrelation Analysis
# ============================================================================
#
# PURPOSE: Test for spatial non-independence in CAFI data and assess
#          implications for statistical models
#
# ANALYSES:
#   1. Global Spatial Autocorrelation (Moran's I)
#      - CAFI abundance, species richness, Shannon diversity
#      - Inverse distance weights
#      - 999 permutation tests
#
#   2. Local Spatial Patterns (LISA)
#      - Local Moran's I for each coral
#      - Cluster identification (High-High, Low-Low, High-Low, Low-High)
#
#   3. Distance-Decay Analysis
#      - Mantel test: community dissimilarity vs geographic distance
#
#   4. Site-Specific Patterns
#      - Separate analyses for HAU, MAT, MRB
#
# DEPENDENCIES: 00_setup.R, coral_master.rds
#
# OUTPUTS:
#   Figures:
#     - spatial_autocorrelation_map.png
#     - distance_decay.png
#     - lisa_clusters.png
#   Tables:
#     - morans_i_results.csv
#     - mantel_test.csv
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-21
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI SURVEY - Spatial Autocorrelation Analysis\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Additional spatial packages - auto-install if needed
cran_mirror <- "https://cran.rstudio.com"
suppressPackageStartupMessages({
  if (!requireNamespace("spdep", quietly = TRUE)) {
    install.packages("spdep", repos = cran_mirror)
  }
  library(spdep)

  if (!requireNamespace("sf", quietly = TRUE)) {
    install.packages("sf", repos = cran_mirror)
  }
  library(sf)
})

# vegan already loaded in setup for Mantel test

cat("[OK] Spatial packages loaded\n\n")

# Output directory for this script
FIG_DIR <- file.path(PATHS$figures, "07_spatial")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")

coral_master <- load_object("coral_master")
community_matrix <- load_object("community_matrix")

cat("  Corals:", nrow(coral_master), "\n")
cat("  Sites:", paste(unique(coral_master$site), collapse = ", "), "\n\n")

# ============================================================================
# PREPARE SPATIAL DATA
# ============================================================================

cat("Preparing spatial data...\n")

# Check for coordinates
has_coords <- all(c("lat", "long") %in% names(coral_master))

if (!has_coords) {
  stop("GPS coordinates (lat, long) not found in coral_master. Cannot perform spatial analysis.")
}

# Filter to corals with valid coordinates
spatial_data <- coral_master %>%
  mutate(
    # Ensure coordinates are numeric
    lat = as.numeric(lat),
    long = as.numeric(long)
  ) %>%
  filter(!is.na(lat), !is.na(long),
         is.finite(lat), is.finite(long)) %>%
  mutate(
    # Ensure diversity metrics are available
    species_richness = replace_na(otu_richness, 0),
    shannon = replace_na(shannon_diversity, 0)
  )

cat("  Corals with GPS coordinates:", nrow(spatial_data), "\n")
cat("  Coordinate ranges:\n")
cat("    Latitude:", round(min(spatial_data$lat), 5), "to", round(max(spatial_data$lat), 5), "\n")
cat("    Longitude:", round(min(spatial_data$long), 5), "to", round(max(spatial_data$long), 5), "\n\n")

# Create sf object for mapping
spatial_sf <- st_as_sf(spatial_data, coords = c("long", "lat"), crs = 4326)

# ============================================================================
# HELPER FUNCTION: Create Inverse Distance Weights
# ============================================================================

create_idw_weights <- function(coords, max_dist = NULL) {
  # Calculate distance matrix
  dist_mat <- as.matrix(dist(coords))

  # Convert to km (approximate at Mo'orea latitude ~17.5S)
  # 1 degree lat ~ 111 km, 1 degree long ~ 106 km at this latitude
  dist_km <- dist_mat * 111

  # Determine bandwidth if not specified
  if (is.null(max_dist)) {
    max_dist <- quantile(dist_km[upper.tri(dist_km)], 0.25)
  }

  # Create inverse distance weights
  n <- nrow(coords)
  nb_list <- vector("list", n)
  weights_list <- vector("list", n)

  for (i in 1:n) {
    # Find neighbors within bandwidth
    neighbors <- which(dist_km[i, ] > 0 & dist_km[i, ] <= max_dist)

    if (length(neighbors) == 0) {
      # Use k-nearest if no neighbors within bandwidth
      neighbors <- order(dist_km[i, ])[2:min(6, n)]
    }

    nb_list[[i]] <- as.integer(neighbors)

    # Inverse distance weights (add small constant to prevent division by zero)
    dists <- dist_km[i, neighbors]
    dists[dists == 0] <- 0.001  # Minimum distance of 1 meter
    w <- 1 / dists
    w[!is.finite(w)] <- 0  # Handle any remaining Inf/NaN
    if (sum(w) > 0) {
      weights_list[[i]] <- w / sum(w)  # Row standardize
    } else {
      weights_list[[i]] <- rep(1/length(neighbors), length(neighbors))  # Equal weights fallback
    }
  }

  # Convert to nb object
  class(nb_list) <- "nb"
  attr(nb_list, "region.id") <- as.character(1:n)

  # Convert to listw
  listw <- nb2listw(nb_list, glist = weights_list, style = "W", zero.policy = TRUE)

  return(list(listw = listw, nb = nb_list, bandwidth = max_dist))
}

# ============================================================================
# GLOBAL SPATIAL AUTOCORRELATION (MORAN'S I)
# ============================================================================

cat("Testing global spatial autocorrelation (Moran's I)...\n")

# Extract coordinates
coords <- as.matrix(spatial_data[, c("long", "lat")])

# Create inverse distance weights
idw <- create_idw_weights(coords)
weights <- idw$listw

cat("  Bandwidth:", round(idw$bandwidth, 2), "km\n")

# Variables to test
test_vars <- c("total_cafi", "species_richness", "shannon")
var_labels <- c("CAFI Abundance", "Species Richness", "Shannon Diversity")

# Store results
moran_results <- data.frame(
  variable = character(),
  morans_i = numeric(),
  expected_i = numeric(),
  variance = numeric(),
  z_score = numeric(),
  p_value = numeric(),
  interpretation = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  label <- var_labels[i]

  # Extract values
  values <- spatial_data[[var]]

  # Moran's I test with permutation
  moran_test <- moran.test(values, weights, randomisation = TRUE)

  # Monte Carlo test for robust p-value
  moran_mc <- moran.mc(values, weights, nsim = 999, zero.policy = TRUE)

  # Calculate z-score
  z_score <- (moran_test$estimate["Moran I statistic"] - moran_test$estimate["Expectation"]) /
              sqrt(moran_test$estimate["Variance"])

  # Interpretation
  if (moran_mc$p.value < 0.05) {
    if (moran_test$estimate["Moran I statistic"] > 0) {
      interp <- "Significant positive autocorrelation (clustering)"
    } else {
      interp <- "Significant negative autocorrelation (dispersion)"
    }
  } else {
    interp <- "No significant spatial autocorrelation"
  }

  moran_results <- rbind(moran_results, data.frame(
    variable = label,
    morans_i = round(moran_test$estimate["Moran I statistic"], 4),
    expected_i = round(moran_test$estimate["Expectation"], 4),
    variance = round(moran_test$estimate["Variance"], 6),
    z_score = round(z_score, 3),
    p_value = round(moran_mc$p.value, 4),
    interpretation = interp,
    stringsAsFactors = FALSE
  ))

  cat("  ", label, ": I =", round(moran_test$estimate["Moran I statistic"], 3),
      "(p =", round(moran_mc$p.value, 4), ")\n")
}

cat("\n")

# Save Moran's I results
save_table(moran_results, "morans_i_results")
cat("  Saved: morans_i_results.csv\n\n")

# ============================================================================
# LOCAL INDICATORS OF SPATIAL ASSOCIATION (LISA)
# ============================================================================

cat("Calculating local spatial autocorrelation (LISA)...\n")

# Local Moran's I for CAFI abundance
local_moran <- localmoran(spatial_data$total_cafi, weights, zero.policy = TRUE)

# Add to spatial data
spatial_data$local_i <- local_moran[, "Ii"]
spatial_data$local_z <- local_moran[, "Z.Ii"]
spatial_data$local_p <- local_moran[, "Pr(z != E(Ii))"]
spatial_data$local_sig <- spatial_data$local_p < 0.05

# Classify cluster types
z_values <- scale(spatial_data$total_cafi)[, 1]
lag_values <- lag.listw(weights, z_values, zero.policy = TRUE)

spatial_data$cluster_type <- case_when(
  !spatial_data$local_sig ~ "Not Significant",
  z_values > 0 & lag_values > 0 ~ "High-High",
  z_values < 0 & lag_values < 0 ~ "Low-Low",
  z_values > 0 & lag_values < 0 ~ "High-Low",
  z_values < 0 & lag_values > 0 ~ "Low-High",
  TRUE ~ "Not Significant"
)

# Summarize clusters
cluster_summary <- spatial_data %>%
  group_by(cluster_type) %>%
  summarise(
    n = n(),
    mean_abundance = mean(total_cafi),
    mean_richness = mean(species_richness),
    .groups = "drop"
  )

cat("  LISA cluster summary:\n")
print(cluster_summary)
cat("\n")

# ============================================================================
# LISA CLUSTER MAP
# ============================================================================

cat("Creating LISA cluster map...\n")

# Define cluster colors
cluster_colors <- c(
  "High-High" = "#D55E00",      # Orange-red (hot spots)
  "Low-Low" = "#0072B2",        # Blue (cold spots)
  "High-Low" = "#E69F00",       # Orange (spatial outliers)
  "Low-High" = "#56B4E9",       # Light blue (spatial outliers)
  "Not Significant" = "gray70"
)

p_lisa <- ggplot(spatial_data, aes(x = long, y = lat)) +
  geom_point(aes(color = cluster_type, size = total_cafi), alpha = 0.7) +
  scale_color_manual(
    values = cluster_colors,
    name = "Cluster Type"
  ) +
  scale_size_continuous(range = c(2, 8), name = "CAFI\nAbundance") +
  labs(
    title = "Local Indicators of Spatial Association (LISA)",
    subtitle = "Spatial clusters of CAFI abundance",
    x = "Longitude",
    y = "Latitude"
  ) +
  coord_quickmap() +
  theme_publication() +
  theme(legend.position = "right")

save_figure(p_lisa, "lisa_clusters", width = 10, height = 8, script_dir = "07_spatial")

# ============================================================================
# SPATIAL AUTOCORRELATION MAP (Combined Overview)
# ============================================================================

cat("Creating spatial autocorrelation overview map...\n")

# Map showing local Moran's I values
p_spatial_overview <- ggplot(spatial_data, aes(x = long, y = lat)) +
  geom_point(aes(color = local_i, size = total_cafi), alpha = 0.7) +
  scale_color_gradient2(
    low = "#0072B2", mid = "gray90", high = "#D55E00",
    midpoint = 0,
    name = "Local\nMoran's I"
  ) +
  scale_size_continuous(range = c(2, 8), name = "CAFI\nAbundance") +
  facet_wrap(~site, scales = "free") +
  labs(
    title = "Spatial Autocorrelation in CAFI Abundance",
    subtitle = paste("Global Moran's I =", round(moran_results$morans_i[1], 3),
                    "| p =", round(moran_results$p_value[1], 4)),
    x = "Longitude",
    y = "Latitude"
  ) +
  coord_quickmap() +
  theme_publication() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

save_figure(p_spatial_overview, "spatial_autocorrelation_map", width = 12, height = 8, script_dir = "07_spatial")

# ============================================================================
# DISTANCE-DECAY ANALYSIS (MANTEL TEST)
# ============================================================================

cat("Performing Mantel test (distance decay)...\n")

# Match community matrix to spatial data
common_corals <- intersect(rownames(community_matrix), spatial_data$coral_id)

if (length(common_corals) > 20) {

  # Subset data
  comm_subset <- community_matrix[common_corals, ]
  spatial_subset <- spatial_data %>%
    filter(coral_id %in% common_corals) %>%
    arrange(match(coral_id, common_corals))

  # Remove zero-sum columns (species not present in subset)
  comm_subset <- comm_subset[, colSums(comm_subset) > 0]

  # Calculate distance matrices
  comm_dist <- vegdist(comm_subset, method = "bray")
  geo_coords <- as.matrix(spatial_subset[, c("long", "lat")])

  # Geographic distance in km
  geo_dist_raw <- as.matrix(dist(geo_coords))
  geo_dist_km <- geo_dist_raw * 111  # Convert to km
  geo_dist <- as.dist(geo_dist_km)

  # Mantel test
  mantel_result <- mantel(comm_dist, geo_dist, method = "pearson", permutations = 999)

  # Save results
  mantel_summary <- data.frame(
    test = "Mantel (Bray-Curtis vs Geographic Distance)",
    mantel_r = round(mantel_result$statistic, 4),
    p_value = round(mantel_result$signif, 4),
    permutations = 999,
    n_samples = length(common_corals),
    interpretation = ifelse(mantel_result$signif < 0.05,
                           "Significant distance decay of community similarity",
                           "No significant distance decay")
  )

  save_table(mantel_summary, "mantel_test")

  cat("  Mantel r:", round(mantel_result$statistic, 4), "\n")
  cat("  P-value:", round(mantel_result$signif, 4), "\n")
  cat("  Interpretation:", mantel_summary$interpretation, "\n\n")

  # Distance-decay plot
  cat("Creating distance-decay plot...\n")

  tryCatch({
    decay_data <- data.frame(
      geographic_km = as.vector(geo_dist),
      community_dissim = as.vector(comm_dist)
    )

    # Check if we have valid data
    if (nrow(decay_data) < 10 || all(is.na(decay_data$geographic_km))) {
      stop("Insufficient valid data for distance-decay plot")
    }

    # Bin data for cleaner visualization
    max_dist <- max(decay_data$geographic_km, na.rm = TRUE)
    if (!is.finite(max_dist) || max_dist <= 0) {
      stop("Invalid distance range for binning")
    }

    decay_data$dist_bin <- cut(decay_data$geographic_km,
                               breaks = seq(0, max_dist + 0.1, length.out = 20))

    decay_binned <- decay_data %>%
      filter(!is.na(dist_bin)) %>%
      group_by(dist_bin) %>%
      summarise(
        mean_dist = mean(geographic_km, na.rm = TRUE),
        mean_dissim = mean(community_dissim, na.rm = TRUE),
        se_dissim = sd(community_dissim, na.rm = TRUE) / sqrt(n()),
        n = n(),
        .groups = "drop"
      ) %>%
      filter(n >= 3)  # Reduced threshold to allow more bins

    # Create the plot
    p_decay <- ggplot(decay_data, aes(x = geographic_km, y = community_dissim)) +
      geom_point(alpha = 0.2, size = 1.5, color = "gray50") +
      geom_smooth(method = "lm", se = TRUE, color = "#D55E00", linewidth = 1.2) +
      labs(
        title = "Distance Decay of Community Similarity",
        subtitle = paste("Mantel r =", round(mantel_result$statistic, 3),
                        "| p =", round(mantel_result$signif, 4)),
        x = "Geographic Distance (km)",
        y = "Community Dissimilarity (Bray-Curtis)"
      ) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_publication()

    # Add binned means if we have enough bins
    if (nrow(decay_binned) >= 3) {
      p_decay <- p_decay +
        geom_point(data = decay_binned, aes(x = mean_dist, y = mean_dissim),
                   color = "#0072B2", size = 3) +
        geom_errorbar(data = decay_binned,
                      aes(x = mean_dist, ymin = mean_dissim - se_dissim,
                          ymax = mean_dissim + se_dissim),
                      width = 0, color = "#0072B2")
    }

    save_figure(p_decay, "distance_decay", width = 8, height = 6, script_dir = "07_spatial")
    cat("  Saved: distance_decay.png\n")

  }, error = function(e) {
    cat("  Warning: Could not create distance-decay plot:", conditionMessage(e), "\n")
    # Create informative placeholder figure
    p_placeholder <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("Distance-decay analysis unavailable\n\n",
                            "Reason:", conditionMessage(e), "\n\n",
                            "This may occur when:\n",
                            "- Insufficient GPS data overlap with community data\n",
                            "- All corals are from the same location"),
               size = 4, hjust = 0.5) +
      theme_void() +
      labs(title = "Distance Decay of Community Similarity")
    save_figure(p_placeholder, "distance_decay", width = 8, height = 6, script_dir = "07_spatial")
    cat("  Saved: distance_decay.png (placeholder with explanation)\n")
  })

} else {
  cat("  Insufficient overlapping samples for Mantel test\n\n")
  mantel_summary <- data.frame(
    test = "Mantel test",
    mantel_r = NA,
    p_value = NA,
    permutations = NA,
    n_samples = length(common_corals),
    interpretation = "Insufficient data"
  )
  save_table(mantel_summary, "mantel_test")
}

# ============================================================================
# SITE-SPECIFIC SPATIAL PATTERNS
# ============================================================================

cat("Analyzing site-specific spatial patterns...\n\n")

sites <- unique(spatial_data$site)
site_moran_results <- data.frame()

for (site_name in sites) {
  cat("  Site:", site_name, "\n")

  site_data <- spatial_data %>% filter(site == site_name)

  if (nrow(site_data) < 20) {
    cat("    Insufficient samples (n =", nrow(site_data), ")\n\n")
    next
  }

  # Create site-specific weights
  site_coords <- as.matrix(site_data[, c("long", "lat")])

  tryCatch({
    site_idw <- create_idw_weights(site_coords)
    site_weights <- site_idw$listw

    # Test each variable
    for (i in seq_along(test_vars)) {
      var <- test_vars[i]
      label <- var_labels[i]

      values <- site_data[[var]]

      # Skip if no variance
      if (sd(values, na.rm = TRUE) == 0) {
        cat("    ", label, ": No variance\n")
        next
      }

      # Moran's I test
      moran_mc <- moran.mc(values, site_weights, nsim = 999, zero.policy = TRUE)

      site_moran_results <- rbind(site_moran_results, data.frame(
        site = site_name,
        variable = label,
        n = nrow(site_data),
        morans_i = round(moran_mc$statistic, 4),
        p_value = round(moran_mc$p.value, 4),
        significant = moran_mc$p.value < 0.05,
        stringsAsFactors = FALSE
      ))

      cat("    ", label, ": I =", round(moran_mc$statistic, 3),
          "(p =", round(moran_mc$p.value, 3), ")\n")
    }
    cat("\n")

  }, error = function(e) {
    cat("    Error:", conditionMessage(e), "\n\n")
  })
}

# Save site-specific results
if (nrow(site_moran_results) > 0) {
  save_table(site_moran_results, "morans_i_by_site")
  cat("  Saved: morans_i_by_site.csv\n\n")
}

# ============================================================================
# MODEL IMPLICATIONS
# ============================================================================

cat("============================================================\n")
cat("MODEL IMPLICATIONS\n")
cat("============================================================\n\n")

# Determine overall autocorrelation status
sig_autocorr <- moran_results %>% filter(p_value < 0.05)

if (nrow(sig_autocorr) > 0) {
  cat("SIGNIFICANT SPATIAL AUTOCORRELATION DETECTED\n\n")

  cat("Variables with significant autocorrelation:\n")
  for (i in 1:nrow(sig_autocorr)) {
    cat("  -", sig_autocorr$variable[i], "(I =", sig_autocorr$morans_i[i],
        ", p =", sig_autocorr$p_value[i], ")\n")
  }
  cat("\n")

  # Estimate variance inflation
  # Using formula: VIF_approx = (1 + rho) / (1 - rho) where rho ~ Moran's I
  max_moran <- max(sig_autocorr$morans_i)
  if (max_moran > 0 && max_moran < 1) {
    vif_approx <- (1 + max_moran) / (1 - max_moran)
    effective_n <- nrow(spatial_data) / vif_approx

    cat("VARIANCE INFLATION:\n")
    cat("  Maximum Moran's I:", round(max_moran, 3), "\n")
    cat("  Approximate variance inflation factor:", round(vif_approx, 2), "\n")
    cat("  Effective sample size:", round(effective_n), "(nominal:", nrow(spatial_data), ")\n\n")
  }

  cat("RECOMMENDATIONS:\n")
  cat("  1. Consider adding spatial random effects to mixed models\n")
  cat("     - Use site as random effect (minimum: random intercept)\n")
  cat("     - For GLMMs: library(glmmTMB) with spatial covariance\n")
  cat("     - For Bayesian: library(brms) with gp() for Gaussian process\n")
  cat("  2. For hypothesis tests:\n")
  cat("     - Use permutation tests that respect spatial structure\n")
  cat("     - Consider Dutilleul's modified t-test for correlations\n")
  cat("  3. Standard errors may be underestimated\n")
  cat("     - Bootstrap confidence intervals recommended\n")
  cat("     - Cluster-robust standard errors by site\n\n")

  # Write recommendations to file
  recommendations <- data.frame(
    category = c("Model Structure", "Model Structure", "Hypothesis Testing",
                 "Standard Errors", "Standard Errors"),
    recommendation = c(
      "Add spatial random effects (site + spatial covariance)",
      "Use glmmTMB or brms for spatial mixed models",
      "Use permutation tests or Dutilleul's modified t-test",
      "Use bootstrap confidence intervals",
      "Consider cluster-robust SEs by site"
    ),
    rationale = c(
      paste("Moran's I =", round(max_moran, 3), "indicates clustering"),
      "These packages support spatial correlation structures",
      "Standard tests assume independence",
      paste("Effective N =", round(effective_n), "vs nominal N =", nrow(spatial_data)),
      "Accounts for within-site dependence"
    )
  )

  save_table(recommendations, "spatial_modeling_recommendations")

} else {
  cat("NO SIGNIFICANT SPATIAL AUTOCORRELATION DETECTED\n\n")
  cat("Standard statistical approaches are appropriate.\n")
  cat("No spatial random effects required, but site as random effect\n")
  cat("is still recommended for hierarchical structure.\n\n")

  recommendations <- data.frame(
    category = "Model Structure",
    recommendation = "Standard mixed models with site random effect",
    rationale = "No spatial autocorrelation detected"
  )

  save_table(recommendations, "spatial_modeling_recommendations")
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("OUTPUTS:\n")
cat("  Tables:\n")
cat("    - morans_i_results.csv (global autocorrelation)\n")
cat("    - mantel_test.csv (distance decay)\n")
cat("    - morans_i_by_site.csv (site-specific patterns)\n")
cat("    - spatial_modeling_recommendations.csv\n")
cat("  Figures:\n")
cat("    - spatial_autocorrelation_map.png\n")
cat("    - lisa_clusters.png\n")
cat("    - distance_decay.png\n\n")

cat("KEY FINDINGS:\n")
cat("  Global Moran's I (CAFI abundance):", moran_results$morans_i[1], "\n")
cat("  P-value:", moran_results$p_value[1], "\n")
if (exists("mantel_summary") && !is.na(mantel_summary$mantel_r)) {
  cat("  Mantel r (distance decay):", mantel_summary$mantel_r, "\n")
  cat("  Mantel p-value:", mantel_summary$p_value, "\n")
}
cat("\n")
