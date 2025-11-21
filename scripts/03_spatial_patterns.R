#!/usr/bin/env Rscript
# ============================================================================
# 03_spatial_patterns.R - Analyze spatial patterns in Survey data
#
# Purpose: Examines geographic distribution, depth patterns, and spatial
#          autocorrelation of CAFI communities across survey sites
#
# Outputs: - Survey location maps
#          - Depth distribution plots
#          - Spatial autocorrelation metrics
#          - Site boundary visualizations
#
# Author: CAFI Analysis Pipeline
# Date: 2025-10-31
# ============================================================================

cat("\n========================================\n")
cat("Spatial Pattern Analysis\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/Survey/00_load_libraries.R"))

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "spatial_patterns")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Spatial Distribution of Survey Sites
# ============================================================================

cat("Analyzing spatial distribution...\n")

# Get unique coral locations
# Each coral has GPS coordinates representing its surveyed location
coral_locations <- metadata %>%
  filter(!is.na(lat) & !is.na(long)) %>%
  distinct(coral_id, .keep_all = TRUE)

# Calculate spatial extent of survey
spatial_summary <- coral_locations %>%
  summarise(
    n_locations = n(),
    lat_min = min(lat),
    lat_max = max(lat),
    long_min = min(long),
    long_max = max(long),
    lat_range = lat_max - lat_min,
    long_range = long_max - long_min
  )

write_csv(spatial_summary,
          file.path(SURVEY_TABLES, "spatial_extent_summary.csv"))

# Plot survey locations with branch width (the actual measurable trait)
# Shows geographic spread of sampling effort across sites
# NOTE: Branch width (tight/wide) is the real morphological trait
# "Morphotype" names (meandrina/eudoxi/verucosa) are NOT confirmed species

# Check if branch_width exists, otherwise use morphotype
if ("branch_width" %in% colnames(coral_locations)) {
  p_survey_map <- ggplot(coral_locations, aes(x = long, y = lat)) +
    geom_point(aes(color = site, shape = branch_width), size = 3, alpha = 0.7) +
    scale_color_viridis_d() +
    scale_shape_manual(values = c("tight" = 16, "wide" = 17),
                       na.value = 4) +
    labs(
      title = "Spatial Distribution of Survey Corals",
      subtitle = paste("Summer 2019 Survey -", nrow(coral_locations), "corals"),
      x = "Longitude",
      y = "Latitude",
      color = "Site",
      shape = "Branch Width"
    ) +
    coord_quickmap() +
    theme_bw()
} else {
  # Fallback to morphotype if branch_width not available
  p_survey_map <- ggplot(coral_locations, aes(x = long, y = lat)) +
    geom_point(aes(color = site, shape = morphotype), size = 3, alpha = 0.7) +
    scale_color_viridis_d() +
    labs(
      title = "Spatial Distribution of Survey Corals",
      subtitle = paste("Summer 2019 Survey -", nrow(coral_locations), "corals"),
      x = "Longitude",
      y = "Latitude",
      color = "Site",
      shape = "Morphotype\n(not species)"
    ) +
    coord_quickmap() +
    theme_bw()
}

ggsave(file.path(fig_dir, "survey_locations_map.png"),
       p_survey_map, width = 10, height = 8, dpi = 300)

cat("  ✓ Survey map created\n")

# ============================================================================
# Distance-Based Analysis
# ============================================================================

cat("Calculating spatial distances...\n")

# Calculate pairwise distances between all surveyed corals
# Uses Haversine formula for accurate distances on Earth's surface
if (nrow(coral_locations) > 1) {
  # Create coordinate matrix for distance calculation
  coords <- coral_locations %>%
    select(long, lat) %>%
    as.matrix()

  # Calculate great circle distances in meters
  dist_matrix <- distm(coords, fun = distHaversine)
  rownames(dist_matrix) <- coral_locations$coral_id
  colnames(dist_matrix) <- coral_locations$coral_id

  # Convert to distance object for downstream analyses
  spatial_dist <- as.dist(dist_matrix)

  # Summarize distance metrics
  dist_summary <- data.frame(
    min_distance_m = min(spatial_dist),
    max_distance_m = max(spatial_dist),
    mean_distance_m = mean(spatial_dist),
    median_distance_m = median(spatial_dist)
  )

  write_csv(dist_summary,
            file.path(SURVEY_TABLES, "spatial_distance_summary.csv"))

  cat("  - Min distance:", round(dist_summary$min_distance_m, 1), "m\n")
  cat("  - Max distance:", round(dist_summary$max_distance_m, 1), "m\n")
  cat("  - Mean distance:", round(dist_summary$mean_distance_m, 1), "m\n")
  cat("  ✓ Distance matrix calculated\n\n")
}

# ============================================================================
# Depth Patterns
# ============================================================================

cat("Analyzing depth patterns...\n")

# Depth is a key environmental gradient affecting CAFI communities
depth_summary <- metadata %>%
  filter(!is.na(depth_m)) %>%
  group_by(site) %>%
  summarise(
    n_corals = n(),
    mean_depth = mean(depth_m, na.rm = TRUE),
    sd_depth = sd(depth_m, na.rm = TRUE),
    min_depth = min(depth_m, na.rm = TRUE),
    max_depth = max(depth_m, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(depth_summary,
          file.path(SURVEY_TABLES, "depth_distribution_by_site.csv"))

# Visualize depth distribution by site
p_depth_dist <- metadata %>%
  filter(!is.na(depth_m)) %>%
  ggplot(aes(x = depth_m, fill = site)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  facet_wrap(~site, scales = "free_y") +
  scale_fill_viridis_d() +
  labs(
    title = "Depth Distribution by Site",
    subtitle = "Survey depth range reflects habitat heterogeneity",
    x = "Depth (m)",
    y = "Number of Corals",
    fill = "Site"
  ) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "depth_distribution.png"),
       p_depth_dist, width = 10, height = 6, dpi = 300)

cat("  ✓ Depth analysis complete\n\n")

# ============================================================================
# Community Patterns Along Environmental Gradients
# ============================================================================

cat("Analyzing community patterns along gradients...\n")

# Merge community metrics with spatial data
# Tests if diversity varies along depth gradient
spatial_community <- coral_locations %>%
  left_join(
    cafi_clean %>%
      group_by(coral_id) %>%
      summarise(
        species_richness = n_distinct(species),
        total_abundance = n(),
        .groups = "drop"
      ),
    by = "coral_id"
  ) %>%
  mutate(
    species_richness = replace_na(species_richness, 0),
    total_abundance = replace_na(total_abundance, 0)
  )

# Plot species richness vs depth with smoothed trend
# LOESS regression reveals nonlinear depth-diversity relationships
if (sum(!is.na(spatial_community$depth_m)) > 10) {
  p_richness_depth <- ggplot(spatial_community,
                             aes(x = depth_m, y = species_richness)) +
    geom_point(aes(color = site), size = 3, alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, color = "black", alpha = 0.2) +
    scale_color_viridis_d() +
    labs(
      title = "OTU Richness Along Depth Gradient",
      subtitle = "LOESS smoothing shows depth-diversity relationship",
      x = "Depth (m)",
      y = "OTU Richness",
      color = "Site"
    )

  ggsave(file.path(fig_dir, "richness_vs_depth.png"),
         p_richness_depth, width = 10, height = 6, dpi = 300)

  cat("  ✓ Depth gradient analysis complete\n\n")
}

# ============================================================================
# Spatial Autocorrelation
# ============================================================================

cat("Testing spatial autocorrelation...\n")

# Moran's I tests if nearby corals have similar community composition
# Significant positive autocorrelation suggests spatial structure
spatial_data <- coral_locations %>%
  left_join(
    cafi_clean %>%
      group_by(coral_id) %>%
      summarise(
        total_cafi = n(),
        species_richness = n_distinct(species),
        .groups = "drop"
      ),
    by = "coral_id"
  ) %>%
  filter(!is.na(species_richness) & !is.na(lat) & !is.na(long))

if (nrow(spatial_data) > 20) {
  # Create spatial points object
  coordinates(spatial_data) <- ~long+lat

  # Calculate inverse distance weights for Moran's I
  # Closer corals get higher weights in spatial correlation
  dists <- as.matrix(dist(coordinates(spatial_data)))
  dists.inv <- 1/dists
  diag(dists.inv) <- 0  # No self-correlation
  dists.inv[is.infinite(dists.inv)] <- 0  # Handle zero distances

  # Compute Moran's I statistic
  # I > 0: positive spatial autocorrelation (clustering)
  # I ≈ 0: random spatial pattern
  # I < 0: negative spatial autocorrelation (dispersion)
  moran_test <- Moran.I(spatial_data$total_cafi, dists.inv)

  moran_summary <- data.frame(
    metric = "total_cafi",
    observed = moran_test$observed,
    expected = moran_test$expected,
    sd = moran_test$sd,
    p_value = moran_test$p.value,
    interpretation = ifelse(
      moran_test$p.value < 0.05,
      ifelse(moran_test$observed > 0, "Clustered", "Dispersed"),
      "Random"
    )
  )

  write_csv(moran_summary,
            file.path(SURVEY_TABLES, "morans_i_spatial_autocorrelation.csv"))

  cat("  - Moran's I:", round(moran_test$observed, 3), "\n")
  cat("  - P-value:", format(moran_test$p.value, scientific = TRUE), "\n")
  cat("  - Pattern:", moran_summary$interpretation, "\n")
  cat("  ✓ Spatial autocorrelation tested\n\n")
}

# ============================================================================
# Site-Level Spatial Patterns
# ============================================================================

cat("Analyzing site-level patterns...\n")

# Calculate site centroids and spatial spread
site_centroids <- coral_locations %>%
  group_by(site) %>%
  summarise(
    n_corals = n(),
    lat_center = mean(lat),
    long_center = mean(long),
    lat_spread = sd(lat),
    long_spread = sd(long),
    .groups = "drop"
  )

write_csv(site_centroids,
          file.path(SURVEY_TABLES, "site_spatial_centroids.csv"))

# Plot site boundaries with 95% confidence ellipses
# Shows spatial clustering of survey effort within each site
# Wrap stat_ellipse in tryCatch to handle sites with too few points
p_site_boundaries <- ggplot(coral_locations, aes(x = long, y = lat)) +
  geom_point(aes(color = site), size = 2, alpha = 0.5) +
  geom_point(data = site_centroids,
             aes(x = long_center, y = lat_center),
             size = 5, shape = 21, fill = "white", color = "black") +
  scale_color_viridis_d() +
  labs(
    title = "Site Boundaries and Centroids",
    subtitle = "Large points show site centers; dots show individual corals",
    x = "Longitude",
    y = "Latitude",
    color = "Site"
  ) +
  coord_quickmap()

# Add ellipses only for sites with sufficient points
# stat_ellipse requires at least 3 points per group
coral_counts_by_site <- coral_locations %>%
  group_by(site) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= 3)

if (nrow(coral_counts_by_site) > 0) {
  sites_with_ellipses <- coral_counts_by_site$site
  data_for_ellipses <- coral_locations %>%
    filter(site %in% sites_with_ellipses)

  p_site_boundaries <- p_site_boundaries +
    stat_ellipse(data = data_for_ellipses,
                 aes(color = site),
                 level = 0.95,
                 linewidth = 1.5)  # Use linewidth instead of deprecated size
}

ggsave(file.path(fig_dir, "site_boundaries.png"),
       p_site_boundaries, width = 10, height = 8, dpi = 300)

cat("  ✓ Site boundaries mapped\n\n")

# ============================================================================
# Branch Width Spatial Distribution
# ============================================================================

cat("Analyzing branch width spatial distribution...\n")

# Test if tight vs wide branching corals partition space differently
# This is the real measurable trait (not "morphotype" which are not confirmed species)
if ("branch_width" %in% colnames(metadata)) {
  branch_spatial <- metadata %>%
    filter(!is.na(branch_width)) %>%
    group_by(site, branch_width) %>%
    summarise(
      n_corals = n(),
      mean_depth = mean(depth_m, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(site) %>%
    mutate(proportion = n_corals / sum(n_corals))
} else {
  # Skip branch width analysis if column not available
  branch_spatial <- NULL
}

# Visualize branch width proportions across sites
if (!is.null(branch_spatial)) {
  p_branch_site <- ggplot(branch_spatial,
                          aes(x = site, y = proportion, fill = branch_width)) +
    geom_col(position = "stack") +
    scale_fill_viridis_d() +
    labs(
      title = "Branch Width Distribution by Site",
      subtitle = "Relative abundance of tight vs wide branching corals",
      x = "Site",
      y = "Proportion",
      fill = "Branch Width"
    )

  ggsave(file.path(fig_dir, "branch_width_by_site.png"),
         p_branch_site, width = 10, height = 6, dpi = 300)

  cat("  ✓ Branch width spatial analysis complete\n\n")
} else {
  cat("  ⚠ Branch width spatial analysis skipped (column not available)\n\n")
}

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Spatial Patterns Summary\n")
cat("========================================\n\n")

cat("Spatial Coverage:\n")
cat("  - Survey locations:", nrow(coral_locations), "\n")
cat("  - Latitude range:", round(spatial_summary$lat_range, 4), "degrees\n")
cat("  - Longitude range:", round(spatial_summary$long_range, 4), "degrees\n")
cat("  - Major sites:", length(unique(coral_locations$site)), "\n\n")

cat("Depth Distribution:\n")
cat("  - Mean depth:", round(mean(metadata$depth_m, na.rm = TRUE), 1), "m\n")
cat("  - Depth range:", round(min(metadata$depth_m, na.rm = TRUE), 1),
    "-", round(max(metadata$depth_m, na.rm = TRUE), 1), "m\n\n")

cat("✅ Spatial pattern analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", SURVEY_TABLES, "\n")
