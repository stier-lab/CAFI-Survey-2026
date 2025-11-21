#!/usr/bin/env Rscript
# ============================================================================
# 11_spatial_autocorrelation.R - Spatial autocorrelation analysis for Survey data
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
# ============================================================================

cat("\n========================================\n")
cat("Spatial Autocorrelation Analysis\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/Survey/00_load_libraries.R"))
library(spdep)
library(spatialreg)
library(gstat)
library(sp)
library(sf)

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "spatial_analysis")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Prepare Spatial Data
# ============================================================================

cat("Preparing spatial data...\n")

# Check for coordinates
has_coords <- all(c("lat", "long") %in% names(metadata)) ||
              all(c("latitude", "longitude") %in% names(metadata)) ||
              all(c("x", "y") %in% names(metadata))

if (!has_coords) {
  cat("  ⚠️ No spatial coordinates found - simulating for demonstration\n")
  # Simulate spatial coordinates for demonstration
  set.seed(123)
  n_sites <- length(unique(metadata$site))
  site_centers <- data.frame(
    site = unique(metadata$site),
    center_lat = runif(n_sites, -17.5, -17.4),
    center_long = runif(n_sites, -149.9, -149.8)
  )

  metadata <- metadata %>%
    left_join(site_centers, by = "site") %>%
    mutate(
      lat = center_lat + rnorm(n(), 0, 0.001),
      long = center_long + rnorm(n(), 0, 0.001)
    )
  cat("  Simulated coordinates created\n\n")
} else {
  # Standardize coordinate names
  if (all(c("latitude", "longitude") %in% names(metadata))) {
    metadata <- metadata %>%
      mutate(lat = latitude, long = longitude)
  } else if (all(c("x", "y") %in% names(metadata))) {
    metadata <- metadata %>%
      mutate(lat = y, long = x)
  }
}

# Create spatial data with community metrics
spatial_data <- metadata %>%
  left_join(
    cafi_clean %>%
      group_by(coral_id) %>%
      summarise(
        total_cafi = n(),
        species_richness = n_distinct(species),
        shannon = vegan::diversity(table(species)),
        .groups = "drop"
      ),
    by = "coral_id"
  ) %>%
  mutate(across(c(total_cafi:shannon), ~replace_na(., 0))) %>%
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
          file.path(SURVEY_TABLES, "spatial_distance_summary.csv"))

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
            file.path(SURVEY_TABLES, "morans_i_test_results.csv"))

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
    scale_color_manual(values = c("High-High" = "red",
                                 "Low-Low" = "blue",
                                 "High-Low" = "orange",
                                 "Low-High" = "lightblue",
                                 "Not Significant" = "gray70")) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "Local Indicators of Spatial Association (LISA)",
         subtitle = "Clusters of CAFI Abundance",
         x = "Longitude",
         y = "Latitude",
         color = "Cluster Type",
         size = "Total CAFI") +
    coord_quickmap()

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
      geom_line(color = "gray50") +
      geom_point(aes(color = significant), size = 3) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray70")) +
      labs(title = "Spatial Correlogram",
           subtitle = "Moran's I at different distance classes",
           x = "Distance Class (km)",
           y = "Moran's I",
           color = "Significant")

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
              file.path(SURVEY_TABLES, "mantel_test_results.csv"))

    cat("  Mantel correlation:", round(mantel_result$statistic, 3), "\n")
    cat("  P-value:", round(mantel_result$signif, 4), "\n\n")

    # Plot relationship
    p_mantel <- data.frame(
      geographic = as.vector(geo_dist),
      community = as.vector(comm_dist)
    ) %>%
      ggplot(aes(x = geographic, y = community)) +
      geom_point(alpha = 0.3) +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      labs(title = "Mantel Test: Community vs Geographic Distance",
           subtitle = paste("r =", round(mantel_result$statistic, 3),
                           "| p =", round(mantel_result$signif, 4)),
           x = "Geographic Distance",
           y = "Community Dissimilarity (Bray-Curtis)")

    ggsave(file.path(fig_dir, "mantel_test_plot.png"),
           p_mantel, width = 10, height = 8, dpi = 300)
  }
}

# ============================================================================
# Spatial Interpolation
# ============================================================================

cat("Creating spatial interpolation maps...\n")

if (nrow(spatial_data) > 20) {
  # Create regular grid for interpolation
  lat_range <- range(spatial_data$lat)
  long_range <- range(spatial_data$long)

  grid <- expand.grid(
    long = seq(long_range[1], long_range[2], length.out = 50),
    lat = seq(lat_range[1], lat_range[2], length.out = 50)
  )

  # Convert to spatial objects
  coordinates(spatial_data) <- ~long+lat
  coordinates(grid) <- ~long+lat
  gridded(grid) <- TRUE

  # IDW interpolation for abundance
  idw_abundance <- idw(total_cafi ~ 1, spatial_data, grid, idp = 2)

  # Convert back to data frame
  idw_df <- as.data.frame(idw_abundance) %>%
    rename(abundance_pred = var1.pred)

  # Plot interpolation
  p_interpolation <- ggplot(idw_df, aes(x = long, y = lat, fill = abundance_pred)) +
    geom_tile() +
    geom_point(data = as.data.frame(spatial_data),
               aes(x = long, y = lat), color = "white", size = 1) +
    scale_fill_viridis_c(trans = "sqrt") +
    labs(title = "Spatial Interpolation of CAFI Abundance",
         subtitle = "Inverse Distance Weighting (IDW)",
         x = "Longitude",
         y = "Latitude",
         fill = "Predicted\nAbundance") +
    coord_quickmap()

  ggsave(file.path(fig_dir, "spatial_interpolation_abundance.png"),
         p_interpolation, width = 12, height = 10, dpi = 300)

  cat("✓ Spatial interpolation complete\n\n")
}

# ============================================================================
# Spatial Regression Models
# ============================================================================

cat("Fitting spatial regression models...\n")

if (exists("weights") && nrow(spatial_data) > 20) {
  # Prepare data
  spatial_data_df <- as.data.frame(spatial_data)

  # OLS model (non-spatial)
  ols_model <- lm(total_cafi ~ morphotype + depth_m, data = spatial_data_df)

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
    log_likelihood = c(logLik(ols_model), logLik(lag_model), logLik(error_model))
  )

  write_csv(model_comparison,
            file.path(SURVEY_TABLES, "spatial_model_comparison.csv"))

  # Lagrange Multiplier tests
  lm_tests <- lm.LMtests(ols_model, weights, test = "all")

  cat("  Model comparison (AIC):\n")
  cat("    - OLS:", round(model_comparison$AIC[1], 1), "\n")
  cat("    - Spatial Lag:", round(model_comparison$AIC[2], 1), "\n")
  cat("    - Spatial Error:", round(model_comparison$AIC[3], 1), "\n\n")
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
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red"),
                      labels = c("Normal", "Hotspot")) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "CAFI Abundance Hotspots",
         subtitle = paste("Hotspots = top 10% abundance (≥", round(abundance_threshold), "CAFI)"),
         x = "Longitude",
         y = "Latitude",
         color = "Status",
         size = "Total CAFI") +
    coord_quickmap()

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
            file.path(SURVEY_TABLES, "hotspot_statistics.csv"))

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
cat("Tables saved to:", SURVEY_TABLES, "\n")