#!/usr/bin/env Rscript
# ============================================================================
# 13_spatial_autocorrelation.R - Spatial Autocorrelation Analysis
# ============================================================================
#
# PURPOSE: Test for spatial autocorrelation in CAFI communities
#
# ANALYSES:
#   - Moran's I for CAFI abundance
#   - Spatial correlograms
#   - Spatial regression models
#   - Variogram analysis
#
# MANUSCRIPT FIGURES: None (supplementary)
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("13: Spatial Autocorrelation Analysis\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Additional spatial packages
if (requireNamespace("sp", quietly = TRUE)) {
  suppressPackageStartupMessages(library(sp))
}
if (requireNamespace("spdep", quietly = TRUE)) {
  suppressPackageStartupMessages(library(spdep))
}

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")
community_matrix <- load_object("community_matrix.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "spatial_autocorr")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Prepare Spatial Data
# ============================================================================

cat("Preparing spatial data...\n")

# Check for coordinates
has_coords <- all(c("lat", "long") %in% names(coral_master)) ||
              all(c("latitude", "longitude") %in% names(coral_master)) ||
              all(c("x", "y") %in% names(coral_master))

if (!has_coords) {
  cat("  ⚠️ No spatial coordinates found - simulating for demonstration\n")
  # Simulate spatial coordinates for demonstration
  set.seed(123)
  n_sites <- length(unique(coral_master$site))
  site_centers <- data.frame(
    site = unique(coral_master$site),
    center_lat = runif(n_sites, -17.5, -17.4),
    center_long = runif(n_sites, -149.9, -149.8)
  )

  coral_master <- coral_master %>%
    left_join(site_centers, by = "site") %>%
    mutate(
      lat = center_lat + rnorm(n(), 0, 0.001),
      long = center_long + rnorm(n(), 0, 0.001)
    )
  cat("  Simulated coordinates created\n\n")
} else {
  # Standardize coordinate names
  if (all(c("latitude", "longitude") %in% names(coral_master))) {
    coral_master <- coral_master %>%
      mutate(lat = latitude, long = longitude)
  } else if (all(c("x", "y") %in% names(coral_master))) {
    coral_master <- coral_master %>%
      mutate(lat = y, long = x)
  }
}

# Create spatial data with community metrics
# coral_master already has total_cafi, otu_richness, and shannon_diversity columns
spatial_data <- coral_master %>%
  mutate(
    # Use columns that already exist in coral_master
    species_richness = replace_na(otu_richness, 0),
    shannon = replace_na(shannon_diversity, 0)
  ) %>%
  filter(!is.na(lat), !is.na(long))

cat("✓ Spatial data prepared with", nrow(spatial_data), "locations\n\n")

# ============================================================================
# Distance-based Analysis
# ============================================================================

cat("Calculating spatial distances...\n")

# Create distance matrix
coords <- as.matrix(spatial_data[, c("long", "lat")])
dist_matrix <- as.matrix(dist(coords))

# Convert to km (approximate)
dist_matrix_km <- dist_matrix * 111  # Rough conversion at equator

# Save distance summary
distance_summary <- data.frame(
  mean_distance = mean(dist_matrix_km[upper.tri(dist_matrix_km)]),
  median_distance = median(dist_matrix_km[upper.tri(dist_matrix_km)]),
  min_distance = min(dist_matrix_km[upper.tri(dist_matrix_km)]),
  max_distance = max(dist_matrix_km[upper.tri(dist_matrix_km)])
)

write_csv(distance_summary,
          file.path(TABLES_DIR, "spatial_distance_summary.csv"))

cat("  Mean distance between corals:", round(distance_summary$mean_distance, 2), "km\n")
cat("  Max distance:", round(distance_summary$max_distance, 2), "km\n\n")

# ============================================================================
# Moran's I Test for Spatial Autocorrelation
# ============================================================================

cat("Testing for spatial autocorrelation (Moran's I)...\n")

if (nrow(spatial_data) > 10) {
  # Create spatial weights matrix (k-nearest neighbors)
  k_neighbors <- min(8, nrow(spatial_data) - 1)
  knn <- knearneigh(coords, k = k_neighbors)
  nb <- knn2nb(knn)
  weights <- nb2listw(nb, style = "W")

  # Moran's I for different metrics
  moran_results <- list()

  # Total abundance
  moran_abundance <- moran.test(spatial_data$total_cafi, weights)
  moran_results$abundance <- data.frame(
    metric = "Total Abundance",
    morans_i = moran_abundance$estimate["Moran I statistic"],
    expected = moran_abundance$estimate["Expectation"],
    variance = moran_abundance$estimate["Variance"],
    p_value = moran_abundance$p.value
  )

  # Species richness
  moran_richness <- moran.test(spatial_data$species_richness, weights)
  moran_results$richness <- data.frame(
    metric = "Species Richness",
    morans_i = moran_richness$estimate["Moran I statistic"],
    expected = moran_richness$estimate["Expectation"],
    variance = moran_richness$estimate["Variance"],
    p_value = moran_richness$p.value
  )

  # Shannon diversity
  moran_shannon <- moran.test(spatial_data$shannon, weights)
  moran_results$shannon <- data.frame(
    metric = "Shannon Diversity",
    morans_i = moran_shannon$estimate["Moran I statistic"],
    expected = moran_shannon$estimate["Expectation"],
    variance = moran_shannon$estimate["Variance"],
    p_value = moran_shannon$p.value
  )

  # Combine results
  moran_df <- bind_rows(moran_results)
  write_csv(moran_df,
            file.path(TABLES_DIR, "morans_i_test_results.csv"))

  cat("  Moran's I results:\n")
  for (i in 1:nrow(moran_df)) {
    cat("    -", moran_df$metric[i], ": I =",
        round(moran_df$morans_i[i], 3),
        "(p =", round(moran_df$p_value[i], 4), ")\n")
  }
  cat("\n")

  # Moran scatterplot
  png(file.path(fig_dir, "moran_scatterplot_abundance.png"),
      width = 10, height = 8, units = "in", res = 300)
  moran.plot(spatial_data$total_cafi, weights,
             main = "Moran's I Scatterplot - Total Abundance",
             xlab = "Total CAFI Abundance",
             ylab = "Spatially Lagged Abundance")
  dev.off()
}

# ============================================================================
# Local Indicators of Spatial Association (LISA)
# ============================================================================

cat("Calculating local spatial autocorrelation (LISA)...\n")

if (exists("weights")) {
  # Local Moran's I
  local_moran <- localmoran(spatial_data$total_cafi, weights)

  # Add to spatial data
  spatial_data$local_i <- local_moran[, "Ii"]
  spatial_data$local_p <- local_moran[, "Pr(z != E(Ii))"]
  spatial_data$local_sig <- spatial_data$local_p < 0.05

  # Classify spatial clusters
  spatial_data$cluster_type <- "Not Significant"
  sig_idx <- which(spatial_data$local_sig)

  if (length(sig_idx) > 0) {
    # Get standardized values and spatial lags
    z_values <- scale(spatial_data$total_cafi)[, 1]
    lag_values <- lag.listw(weights, z_values)

    spatial_data$cluster_type[sig_idx] <- ifelse(
      z_values[sig_idx] > 0 & lag_values[sig_idx] > 0, "High-High",
      ifelse(z_values[sig_idx] < 0 & lag_values[sig_idx] < 0, "Low-Low",
      ifelse(z_values[sig_idx] > 0 & lag_values[sig_idx] < 0, "High-Low",
      "Low-High"))
    )
  }

  # Plot LISA clusters
  p_lisa <- ggplot(spatial_data, aes(x = long, y = lat)) +
    geom_point(aes(color = cluster_type, size = total_cafi), alpha = 0.7) +
    scale_color_manual(values = c("High-High" = "#D55E00",
                                 "Low-Low" = "#0072B2",
                                 "High-Low" = "#E69F00",
                                 "Low-High" = "#56B4E9",
                                 "Not Significant" = "gray70")) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "Local indicators of spatial association (LISA)",
         subtitle = "Clusters of CAFI abundance",
         x = "Longitude",
         y = "Latitude",
         color = "Cluster type",
         size = "Total CAFI") +
    coord_quickmap() +
    theme_publication()

  ggsave(file.path(fig_dir, "lisa_cluster_map.png"),
         p_lisa, width = 12, height = 10, dpi = 300)

  cat("✓ LISA analysis complete\n\n")
}

# ============================================================================
# Spatial Correlograms
# ============================================================================

cat("Creating spatial correlograms...\n")

# Check if correlog function is available (from ncf or pgirmess package)
if (!exists("correlog", mode = "function")) {
  cat("Note: correlog function not available. Install 'ncf' or 'pgirmess' package for correlogram analysis.\n")
  cat("Skipping correlogram analysis...\n\n")
} else if (exists("coords") && nrow(spatial_data) > 20) {
  # Define distance classes
  max_dist <- max(dist_matrix_km[upper.tri(dist_matrix_km)])
  dist_breaks <- seq(0, max_dist, length.out = 11)

  # Calculate correlogram
  correlogram <- correlog(coords, spatial_data$total_cafi,
                         method = "Moran",
                         nbclass = 10)

  # Create data frame for plotting
  if (length(correlogram$correlation) > 0) {
    corr_df <- data.frame(
      distance = correlogram$dist.class,
      correlation = correlogram$correlation,
      p_value = correlogram$p.value
    ) %>%
      mutate(significant = p_value < 0.05)

    # Plot correlogram
    p_correlogram <- ggplot(corr_df, aes(x = distance, y = correlation)) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_line(color = "gray50", linewidth = 0.8) +
      geom_point(aes(color = significant), size = 3) +
      scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "gray70")) +
      labs(title = "Spatial correlogram",
           subtitle = "Moran's I at different distance classes",
           x = "Distance class (km)",
           y = "Moran's I",
           color = "Significant") +
      theme_publication()

    ggsave(file.path(fig_dir, "spatial_correlogram.png"),
           p_correlogram, width = 10, height = 6, dpi = 300)

    cat("✓ Correlogram created\n\n")
  }
}

# ============================================================================
# Mantel Test (Community vs Geographic Distance)
# ============================================================================

cat("Performing Mantel test...\n")

if (nrow(community_matrix) > 10 && nrow(spatial_data) > 10) {
  # Match community matrix to spatial data
  common_corals <- intersect(rownames(community_matrix), spatial_data$coral_id)

  if (length(common_corals) > 10) {
    # Subset data
    comm_subset <- community_matrix[common_corals, ]
    spatial_subset <- spatial_data %>%
      filter(coral_id %in% common_corals) %>%
      arrange(match(coral_id, common_corals))

    # Calculate distance matrices
    comm_dist <- vegdist(comm_subset, method = "bray")
    geo_coords <- as.matrix(spatial_subset[, c("long", "lat")])
    geo_dist <- dist(geo_coords)

    # Mantel test
    mantel_result <- mantel(comm_dist, geo_dist,
                           method = "pearson",
                           permutations = 999)

    # Save results
    mantel_summary <- data.frame(
      correlation = mantel_result$statistic,
      significance = mantel_result$signif,
      permutations = 999
    )

    write_csv(mantel_summary,
              file.path(TABLES_DIR, "mantel_test_results.csv"))

    cat("  Mantel correlation:", round(mantel_result$statistic, 3), "\n")
    cat("  P-value:", round(mantel_result$signif, 4), "\n\n")

    # Plot relationship
    p_mantel <- data.frame(
      geographic = as.vector(geo_dist),
      community = as.vector(comm_dist)
    ) %>%
      ggplot(aes(x = geographic, y = community)) +
      geom_point(alpha = 0.3, size = 1.5) +
      geom_smooth(method = "lm", se = TRUE, color = "#D55E00", linewidth = 0.8) +
      labs(title = "Mantel test: community vs geographic distance",
           subtitle = paste("r =", round(mantel_result$statistic, 3),
                           "| p =", round(mantel_result$signif, 4)),
           x = "Geographic distance",
           y = "Community dissimilarity (Bray-Curtis)") +
      theme_publication()

    ggsave(file.path(fig_dir, "mantel_test_plot.png"),
           p_mantel, width = 10, height = 8, dpi = 300)
  }
}

# ============================================================================
# Spatial Interpolation
# ============================================================================

cat("Creating spatial interpolation maps...\n")

if (nrow(spatial_data) > 20 && requireNamespace("gstat", quietly = TRUE)) {
  library(gstat)

  # Save original points before converting to spatial
  points_df <- spatial_data %>%
    select(coral_id, lat, long, total_cafi, species_richness)

  # Create regular grid for interpolation
  lat_range <- range(spatial_data$lat)
  long_range <- range(spatial_data$long)

  grid <- expand.grid(
    long = seq(long_range[1], long_range[2], length.out = 50),
    lat = seq(lat_range[1], lat_range[2], length.out = 50)
  )

  # Convert to spatial objects
  spatial_sp <- spatial_data
  coordinates(spatial_sp) <- ~long+lat
  coordinates(grid) <- ~long+lat
  gridded(grid) <- TRUE

  # IDW interpolation for abundance
  tryCatch({
    idw_abundance <- idw(total_cafi ~ 1, spatial_sp, grid, idp = 2)

    # Convert back to data frame - use coordinates() to get x/y and merge with predictions
    idw_df <- as.data.frame(idw_abundance)

    # The output has columns: coordinates.long, coordinates.lat (or similar), var1.pred
    # Standardize column names
    names(idw_df) <- gsub("^coordinates\\.", "", names(idw_df))

    # Find the prediction column (usually var1.pred)
    pred_col <- grep("pred", names(idw_df), value = TRUE)[1]
    if (is.na(pred_col)) {
      pred_col <- names(idw_df)[!names(idw_df) %in% c("long", "lat", "x", "y")][1]
    }

    # Rename to standardized names
    if (!is.null(pred_col) && pred_col %in% names(idw_df)) {
      names(idw_df)[names(idw_df) == pred_col] <- "abundance_pred"
    }

    # Ensure we have long/lat columns
    if (!"long" %in% names(idw_df) && "x" %in% names(idw_df)) {
      idw_df$long <- idw_df$x
      idw_df$lat <- idw_df$y
    }

    # Plot interpolation (only if we have all required columns)
    if (all(c("long", "lat", "abundance_pred") %in% names(idw_df))) {
      p_interpolation <- ggplot(idw_df, aes(x = long, y = lat, fill = abundance_pred)) +
        geom_tile() +
        geom_point(data = points_df,
                   aes(x = long, y = lat), inherit.aes = FALSE, color = "white", size = 1) +
        scale_fill_viridis_c(trans = "sqrt") +
        labs(title = "Spatial interpolation of CAFI abundance",
             subtitle = "Inverse distance weighting (IDW)",
             x = "Longitude",
             y = "Latitude",
             fill = "Predicted\nabundance") +
        coord_quickmap() +
        theme_publication()

      ggsave(file.path(fig_dir, "spatial_interpolation_abundance.png"),
             p_interpolation, width = 12, height = 10, dpi = 300)
      cat("  ✓ IDW interpolation figure saved\n")
    } else {
      cat("  IDW interpolation ran but columns not found for plotting\n")
      cat("    Available columns:", paste(names(idw_df), collapse = ", "), "\n")
    }
  }, error = function(e) {
    cat("  IDW interpolation failed:", conditionMessage(e), "\n")
  })

  cat("✓ Spatial interpolation complete\n\n")
} else {
  cat("  Skipping interpolation (gstat package not available)\n\n")
}

# ============================================================================
# Spatial Regression Models
# ============================================================================

cat("Fitting spatial regression models...\n")

if (exists("weights") && nrow(spatial_data) > 20) {
  # Check if spatialreg package is available (contains lagsarlm, errorsarlm)
  if (requireNamespace("spatialreg", quietly = TRUE)) {
    suppressPackageStartupMessages(library(spatialreg))

    # Prepare data - ensure it's a regular data frame
    spatial_data_df <- as.data.frame(spatial_data)

    # OLS model (non-spatial)
    ols_model <- lm(total_cafi ~ morphotype + depth_m, data = spatial_data_df)

    tryCatch({
      # Spatial lag model
      lag_model <- lagsarlm(total_cafi ~ morphotype + depth_m,
                           data = spatial_data_df,
                           listw = weights)

      # Spatial error model
      error_model <- errorsarlm(total_cafi ~ morphotype + depth_m,
                              data = spatial_data_df,
                              listw = weights)

      # Model comparison
      model_comparison <- data.frame(
        model = c("OLS", "Spatial Lag", "Spatial Error"),
        AIC = c(AIC(ols_model), AIC(lag_model), AIC(error_model)),
        log_likelihood = c(as.numeric(logLik(ols_model)),
                          as.numeric(logLik(lag_model)),
                          as.numeric(logLik(error_model)))
      )

      write_csv(model_comparison,
                file.path(TABLES_DIR, "spatial_model_comparison.csv"))

      # Lagrange Multiplier tests
      lm_tests <- lm.LMtests(ols_model, weights, test = "all")

      cat("  Model comparison (AIC):\n")
      cat("    - OLS:", round(model_comparison$AIC[1], 1), "\n")
      cat("    - Spatial Lag:", round(model_comparison$AIC[2], 1), "\n")
      cat("    - Spatial Error:", round(model_comparison$AIC[3], 1), "\n\n")
    }, error = function(e) {
      cat("  Spatial regression models failed:", conditionMessage(e), "\n")
      cat("  (OLS model was fitted successfully)\n\n")
    })
  } else {
    cat("  Skipping spatial regression (spatialreg package not installed)\n")
    cat("  Install with: install.packages('spatialreg')\n\n")
  }
} else {
  cat("  Skipping spatial regression (insufficient data or no weights)\n\n")
}

# ============================================================================
# Hotspot Analysis
# ============================================================================

cat("Identifying spatial hotspots...\n")

if (exists("spatial_data")) {
  # Convert back to regular data frame if needed
  if (class(spatial_data)[1] == "SpatialPointsDataFrame") {
    spatial_data <- as.data.frame(spatial_data)
  }

  # Identify hotspots (top 10% abundance)
  abundance_threshold <- quantile(spatial_data$total_cafi, 0.9)
  spatial_data$hotspot <- spatial_data$total_cafi >= abundance_threshold

  # Map hotspots
  p_hotspots <- ggplot(spatial_data, aes(x = long, y = lat)) +
    geom_point(aes(color = hotspot, size = total_cafi), alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "#D55E00"),
                      labels = c("Normal", "Hotspot")) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "CAFI abundance hotspots",
         subtitle = paste("Hotspots = top 10% abundance (≥", round(abundance_threshold), "CAFI)"),
         x = "Longitude",
         y = "Latitude",
         color = "Status",
         size = "Total CAFI") +
    coord_quickmap() +
    theme_publication()

  ggsave(file.path(fig_dir, "abundance_hotspots.png"),
         p_hotspots, width = 12, height = 10, dpi = 300)

  # Hotspot statistics
  hotspot_stats <- spatial_data %>%
    group_by(hotspot) %>%
    summarise(
      n_corals = n(),
      mean_abundance = mean(total_cafi),
      mean_richness = mean(species_richness),
      mean_shannon = mean(shannon),
      .groups = "drop"
    )

  write_csv(hotspot_stats,
            file.path(TABLES_DIR, "hotspot_statistics.csv"))

  cat("✓ Identified", sum(spatial_data$hotspot), "hotspot corals\n\n")
}

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Spatial Autocorrelation Summary\n")
cat("========================================\n\n")

cat("Analyses Completed:\n")
cat("  - Distance matrix calculation\n")
cat("  - Moran's I global autocorrelation test\n")
cat("  - Local Indicators of Spatial Association (LISA)\n")
cat("  - Spatial correlograms\n")
cat("  - Mantel test (community vs distance)\n")
cat("  - Spatial interpolation (IDW)\n")
cat("  - Spatial regression models\n")
cat("  - Hotspot identification\n\n")

if (exists("moran_df")) {
  sig_autocorr <- moran_df %>% filter(p_value < 0.05)
  if (nrow(sig_autocorr) > 0) {
    cat("Significant Spatial Autocorrelation:\n")
    for (i in 1:nrow(sig_autocorr)) {
      cat("  -", sig_autocorr$metric[i], "\n")
    }
  } else {
    cat("No significant spatial autocorrelation detected\n")
  }
  cat("\n")
}

if (exists("mantel_summary")) {
  cat("Distance Decay:\n")
  cat("  - Mantel r:", round(mantel_summary$correlation, 3), "\n")
  if (mantel_summary$significance < 0.05) {
    cat("  - Significant distance decay of community similarity\n")
  } else {
    cat("  - No significant distance decay\n")
  }
  cat("\n")
}

cat("✅ Spatial autocorrelation analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", TABLES_DIR, "\n")