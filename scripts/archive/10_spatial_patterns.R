#!/usr/bin/env Rscript
# ============================================================================
# 10_spatial_patterns.R - Spatial Pattern Analysis
# ============================================================================
#
# PURPOSE: Examine geographic distribution and spatial patterns of corals/CAFI
#
# ANALYSES:
#   - Survey location maps (Mo'orea)
#   - Depth distribution by site
#   - Spatial extent of sampling
#   - Site boundary visualizations
#
# MANUSCRIPT FIGURES: None (supplementary maps)
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("10: Spatial Pattern Analysis\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Additional spatial packages
if (requireNamespace("sp", quietly = TRUE)) {
  suppressPackageStartupMessages(library(sp))
}

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "spatial")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Spatial Distribution of Survey Sites
# ============================================================================

cat("Analyzing spatial distribution...\n")

# Get unique coral locations
# Each coral has GPS coordinates representing its surveyed location
coral_locations <- coral_master %>%
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
          file.path(TABLES_DIR, "spatial_extent_summary.csv"))

# Plot survey locations with branch width (the actual measurable trait)
# Shows geographic spread of sampling effort across sites
# NOTE: Branch width (tight/wide) is the real morphological trait
# "Morphotype" names (meandrina/eudoxi/verucosa) are NOT confirmed species

# Check if branch_width exists, otherwise use morphotype
if ("branch_width" %in% colnames(coral_locations)) {
  p_survey_map <- ggplot(coral_locations, aes(x = long, y = lat)) +
    geom_point(aes(color = site, shape = branch_width), size = 3, alpha = 0.7) +
    scale_color_site() +
    scale_shape_manual(values = c("tight" = 16, "wide" = 17),
                       na.value = 4) +
    labs(
      title = "Spatial distribution of survey corals",
      subtitle = paste("Summer 2019 survey -", nrow(coral_locations), "corals"),
      x = "Longitude",
      y = "Latitude",
      color = "Site",
      shape = "Branch width"
    ) +
    coord_quickmap() +
    theme_publication()
} else {
  # Fallback to morphotype if branch_width not available
  p_survey_map <- ggplot(coral_locations, aes(x = long, y = lat)) +
    geom_point(aes(color = site, shape = morphotype), size = 3, alpha = 0.7) +
    scale_color_site() +
    labs(
      title = "Spatial distribution of survey corals",
      subtitle = paste("Summer 2019 survey -", nrow(coral_locations), "corals"),
      x = "Longitude",
      y = "Latitude",
      color = "Site",
      shape = "Morphotype\n(not species)"
    ) +
    coord_quickmap() +
    theme_publication()
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
  # Use geosphere package if available, otherwise use simple Euclidean
  if (requireNamespace("geosphere", quietly = TRUE)) {
    dist_matrix <- geosphere::distm(coords, fun = geosphere::distHaversine)
  } else {
    # Fallback: Euclidean distance (approximate for small areas)
    cat("  Note: geosphere package not available, using Euclidean distance\n")
    dist_matrix <- as.matrix(dist(coords))
  }
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
            file.path(TABLES_DIR, "spatial_distance_summary.csv"))

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
depth_summary <- coral_master %>%
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
          file.path(TABLES_DIR, "depth_distribution_by_site.csv"))

# Visualize depth distribution by site
p_depth_dist <- coral_master %>%
  filter(!is.na(depth_m)) %>%
  ggplot(aes(x = depth_m, fill = site)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  facet_wrap(~site, scales = "free_y") +
  scale_fill_site() +
  labs(
    title = "Depth distribution by site",
    subtitle = "Survey depth range reflects habitat heterogeneity",
    x = "Depth (m)",
    y = "Number of corals",
    fill = "Site"
  ) +
  theme_publication() +
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
# coral_master already has otu_richness and total_cafi columns
spatial_community <- coral_locations %>%
  mutate(
    species_richness = replace_na(otu_richness, 0),
    total_abundance = replace_na(total_cafi, 0)
  )

# Plot species richness vs depth with smoothed trend
# LOESS regression reveals nonlinear depth-diversity relationships
if (sum(!is.na(spatial_community$depth_m)) > 10) {
  p_richness_depth <- ggplot(spatial_community,
                             aes(x = depth_m, y = species_richness)) +
    geom_point(aes(color = site, shape = site), size = 3, alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, color = "black", alpha = 0.2,
                linewidth = 0.8) +
    scale_color_site() +
    scale_shape_site() +
    labs(
      title = "OTU richness along depth gradient",
      subtitle = "LOESS smoothing shows depth-diversity relationship",
      x = "Depth (m)",
      y = "OTU richness",
      color = "Site",
      shape = "Site"
    ) +
    theme_publication()

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
# Use spatial_community which has species_richness from earlier processing
spatial_data <- spatial_community %>%
  filter(!is.na(lat) & !is.na(long))

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
            file.path(TABLES_DIR, "morans_i_spatial_autocorrelation.csv"))

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
          file.path(TABLES_DIR, "site_spatial_centroids.csv"))

# Plot site boundaries with 95% confidence ellipses
# Shows spatial clustering of survey effort within each site
# Wrap stat_ellipse in tryCatch to handle sites with too few points
p_site_boundaries <- ggplot(coral_locations, aes(x = long, y = lat)) +
  geom_point(aes(color = site, shape = site), size = 2, alpha = 0.5) +
  geom_point(data = site_centroids,
             aes(x = long_center, y = lat_center),
             size = 5, shape = 21, fill = "white", color = "black") +
  scale_color_site() +
  scale_shape_site() +
  labs(
    title = "Site boundaries and centroids",
    subtitle = "Large points show site centers; dots show individual corals",
    x = "Longitude",
    y = "Latitude",
    color = "Site",
    shape = "Site"
  ) +
  coord_quickmap() +
  theme_publication()

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
if ("branch_width" %in% colnames(coral_master)) {
  branch_spatial <- coral_master %>%
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
      title = "Branch width distribution by site",
      subtitle = "Relative abundance of tight vs wide branching corals",
      x = "Site",
      y = "Proportion",
      fill = "Branch width"
    ) +
    theme_publication()

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
cat("  - Mean depth:", round(mean(coral_master$depth_m, na.rm = TRUE), 1), "m\n")
cat("  - Depth range:", round(min(coral_master$depth_m, na.rm = TRUE), 1),
    "-", round(max(coral_master$depth_m, na.rm = TRUE), 1), "m\n\n")

# ============================================================================
# STATISTICAL RESULTS SUMMARY
# ============================================================================

cat("Compiling statistical results...\n")

stats_results <- init_results_df()

# 1. Moran's I spatial autocorrelation
if(exists("moran_test")) {
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H_spatial",
      question = "Is CAFI abundance spatially autocorrelated?",
      test_name = "Moran's I",
      test_statistic = moran_test$observed,
      df = "Inverse distance weights",
      p_value = moran_test$p.value,
      effect_size = moran_test$observed - moran_test$expected,
      effect_type = "I - E[I]",
      n = nrow(spatial_data),
      notes = moran_summary$interpretation
    )
  )
}

# 2. Depth-richness correlation (if richness_vs_depth was created)
if(sum(!is.na(spatial_community$depth_m)) > 10) {
  depth_richness_data <- spatial_community %>%
    filter(!is.na(depth_m), !is.na(species_richness))

  if(nrow(depth_richness_data) > 10) {
    cor_depth_rich <- cor.test(depth_richness_data$depth_m,
                                depth_richness_data$species_richness,
                                method = "spearman")

    stats_results <- bind_rows(stats_results,
      create_result_row(
        hypothesis = "H_depth_rich",
        question = "Does CAFI richness vary with depth?",
        test_name = "Spearman correlation",
        test_statistic = cor_depth_rich$statistic,
        df = nrow(depth_richness_data) - 2,
        p_value = cor_depth_rich$p.value,
        effect_size = cor_depth_rich$estimate,
        effect_type = "rho",
        n = nrow(depth_richness_data),
        notes = "Correlation: depth vs OTU richness"
      )
    )
  }
}

# 3. Kruskal-Wallis for depth by site
if(sum(!is.na(coral_master$depth_m) & !is.na(coral_master$site)) > 10) {
  depth_site_data <- coral_master %>%
    filter(!is.na(depth_m), !is.na(site))

  kw_depth_site <- kruskal.test(depth_m ~ site, data = depth_site_data)
  n_total <- nrow(depth_site_data)
  k <- n_distinct(depth_site_data$site)
  epsilon_sq <- kw_depth_site$statistic / (n_total - 1)

  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H_depth_site",
      question = "Does survey depth differ among sites?",
      test_name = "Kruskal-Wallis",
      test_statistic = kw_depth_site$statistic,
      df = k - 1,
      p_value = kw_depth_site$p.value,
      effect_size = as.numeric(epsilon_sq),
      effect_type = "ε²",
      n = n_total,
      notes = "Tests for depth differences across sites"
    )
  )
}

# 4. Chi-square for branch width by site (if available)
if(!is.null(branch_spatial) && nrow(branch_spatial) > 0) {
  branch_contingency <- coral_master %>%
    filter(!is.na(branch_width), !is.na(site)) %>%
    count(site, branch_width) %>%
    pivot_wider(names_from = branch_width, values_from = n, values_fill = 0)

  if(nrow(branch_contingency) > 1 && ncol(branch_contingency) > 2) {
    chi_branch <- chisq.test(branch_contingency %>% select(-site) %>% as.matrix())

    n <- sum(branch_contingency %>% select(-site))
    cramers_v <- sqrt(chi_branch$statistic / (n * (min(nrow(branch_contingency), ncol(branch_contingency) - 1) - 1)))

    stats_results <- bind_rows(stats_results,
      create_result_row(
        hypothesis = "H_branch_site",
        question = "Does branch width distribution differ among sites?",
        test_name = "Chi-square test",
        test_statistic = chi_branch$statistic,
        df = chi_branch$parameter,
        p_value = chi_branch$p.value,
        effect_size = as.numeric(cramers_v),
        effect_type = "Cramér's V",
        n = n,
        notes = "Tests for non-random branch width distribution"
      )
    )
  }
}

# 5. Spatial extent statistics (descriptive)
stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "Desc",
    question = "What is the spatial extent of the survey?",
    test_name = "Descriptive statistics",
    test_statistic = spatial_summary$n_locations,
    df = NA,
    p_value = NA,
    effect_size = spatial_summary$lat_range,
    effect_type = "Lat range (deg)",
    n = nrow(coral_locations),
    notes = sprintf("Long range: %.4f deg", spatial_summary$long_range)
  )
)

# Save statistical results
save_stats_summary(stats_results, "10_spatial_patterns", "Spatial Pattern Analysis")

cat("✅ Spatial pattern analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", TABLES_DIR, "\n")
