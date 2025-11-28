#!/usr/bin/env Rscript
# ============================================================================
# investigate_nmds_outliers.R - Investigate NMDS outliers for anomalies
#
# Purpose: Identify corals that are outliers in NMDS ordination space and
# examine whether they have partial mortality, unusual CAFI communities,
# or other characteristics that explain their extreme positions.
#
# Author: CAFI Analysis Pipeline
# Date: 2025-11-22
# ============================================================================

cat("\n========================================\n")
cat("NMDS Outlier Investigation\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/00_load_libraries.R"))

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))

# Also load raw coral data to get field_obs and notes
coral_raw <- read_csv(here::here("data/1. survey_coral_characteristics_merged_v2.csv"),
                      show_col_types = FALSE)

# Create output directory
fig_dir <- file.path(SURVEY_FIGURES, "nmds_outliers")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# STEP 1: Perform NMDS and identify outliers
# ============================================================================

cat("STEP 1: Performing NMDS ordination...\n")

# Remove empty rows/columns
community_matrix_clean <- community_matrix[rowSums(community_matrix) > 0, ]
community_matrix_clean <- community_matrix_clean[, colSums(community_matrix_clean) > 0]

# Perform NMDS
set.seed(123)
nmds_bray <- metaMDS(community_matrix_clean, distance = "bray", k = 2,
                     trymax = 100, autotransform = FALSE)

cat("  NMDS stress:", round(nmds_bray$stress, 3), "\n")

# Extract scores
nmds_scores <- as.data.frame(scores(nmds_bray, display = "sites"))
nmds_scores$coral_id <- rownames(nmds_scores)

# Calculate distance from centroid (multivariate outlier measure)
nmds_scores <- nmds_scores %>%
  mutate(
    centroid_dist = sqrt((NMDS1 - mean(NMDS1))^2 + (NMDS2 - mean(NMDS2))^2)
  )

# Identify outliers (top 10% by distance from centroid)
outlier_threshold <- quantile(nmds_scores$centroid_dist, 0.90)
nmds_scores <- nmds_scores %>%
  mutate(
    is_outlier = centroid_dist > outlier_threshold,
    outlier_rank = rank(-centroid_dist)
  )

n_outliers <- sum(nmds_scores$is_outlier)
cat("  Identified", n_outliers, "outliers (top 10% by distance from centroid)\n\n")

# ============================================================================
# STEP 2: Merge with coral metadata
# ============================================================================

cat("STEP 2: Merging with coral characteristics...\n")

# Merge with metadata and survey_master (for volume and other data)
nmds_full <- nmds_scores %>%
  left_join(metadata, by = "coral_id") %>%
  left_join(
    survey_master %>%
      select(coral_id, volume_lab, notes,
             mean_live_volume_of_neighbors, mean_total_volume_of_neighbors),
    by = "coral_id"
  )

# Calculate CAFI metrics for each coral
cafi_summary <- cafi_clean %>%
  group_by(coral_id) %>%
  summarise(
    cafi_abundance = n(),
    cafi_richness = n_distinct(species),
    n_crabs = sum(type == "crab"),
    n_shrimp = sum(type == "shrimp"),
    n_snails = sum(type == "snail"),
    n_fish = sum(type == "fish"),
    prop_crabs = n_crabs / cafi_abundance,
    prop_shrimp = n_shrimp / cafi_abundance,
    .groups = "drop"
  )

nmds_full <- nmds_full %>%
  left_join(cafi_summary, by = "coral_id")

# ============================================================================
# STEP 3: Compare outliers vs non-outliers
# ============================================================================

cat("STEP 3: Comparing outlier characteristics...\n\n")

# Summary comparison
comparison <- nmds_full %>%
  group_by(is_outlier) %>%
  summarise(
    n = n(),
    mean_volume = mean(volume_lab, na.rm = TRUE),
    sd_volume = sd(volume_lab, na.rm = TRUE),
    mean_depth = mean(depth_m, na.rm = TRUE),
    mean_abundance = mean(cafi_abundance, na.rm = TRUE),
    mean_richness = mean(cafi_richness, na.rm = TRUE),
    pct_wide_branch = mean(branch_width == "wide", na.rm = TRUE) * 100,
    mean_prop_crabs = mean(prop_crabs, na.rm = TRUE),
    mean_prop_shrimp = mean(prop_shrimp, na.rm = TRUE),
    .groups = "drop"
  )

cat("COMPARISON: Outliers vs Non-outliers\n")
cat("=====================================\n")
print(comparison)
cat("\n")

# Statistical tests
cat("Statistical tests:\n")

# Volume
vol_test <- wilcox.test(volume_lab ~ is_outlier, data = nmds_full)
cat("  Volume: Wilcoxon p =", round(vol_test$p.value, 3), "\n")

# Abundance
abund_test <- wilcox.test(cafi_abundance ~ is_outlier, data = nmds_full)
cat("  Abundance: Wilcoxon p =", round(abund_test$p.value, 3), "\n")

# Richness
rich_test <- wilcox.test(cafi_richness ~ is_outlier, data = nmds_full)
cat("  Richness: Wilcoxon p =", round(rich_test$p.value, 3), "\n")

# Site distribution
site_test <- chisq.test(table(nmds_full$is_outlier, nmds_full$site))
cat("  Site distribution: Chi-squared p =", round(site_test$p.value, 3), "\n\n")

# ============================================================================
# STEP 4: Examine individual outliers
# ============================================================================

cat("STEP 4: Examining top outliers...\n\n")

top_outliers <- nmds_full %>%
  filter(is_outlier) %>%
  arrange(outlier_rank) %>%
  select(coral_id, outlier_rank, centroid_dist, NMDS1, NMDS2,
         site, morphotype, branch_width, volume_lab, depth_m,
         cafi_abundance, cafi_richness, prop_crabs, prop_shrimp,
         field_obs, notes)

cat("TOP 10 OUTLIERS:\n")
cat("================\n\n")

for (i in 1:min(10, nrow(top_outliers))) {
  coral <- top_outliers[i, ]
  cat(sprintf("Rank %d: %s\n", coral$outlier_rank, coral$coral_id))
  cat(sprintf("  Site: %s | Morphotype: %s | Branch: %s\n",
              coral$site, coral$morphotype, coral$branch_width))
  cat(sprintf("  Volume: %.0f cm³ | Depth: %.1f m\n",
              coral$volume_lab, coral$depth_m))
  cat(sprintf("  CAFI abundance: %d | Richness: %d\n",
              coral$cafi_abundance, coral$cafi_richness))
  cat(sprintf("  Composition: %.0f%% crabs, %.0f%% shrimp\n",
              coral$prop_crabs * 100, coral$prop_shrimp * 100))
  cat(sprintf("  NMDS position: (%.2f, %.2f) | Distance: %.2f\n",
              coral$NMDS1, coral$NMDS2, coral$centroid_dist))

  if (!is.na(coral$field_obs) && coral$field_obs != "") {
    cat(sprintf("  Field obs: %s\n", coral$field_obs))
  }
  if (!is.na(coral$notes) && coral$notes != "") {
    cat(sprintf("  Notes: %s\n", coral$notes))
  }
  cat("\n")
}

# Save outlier details
write_csv(top_outliers,
          file.path(SURVEY_TABLES, "nmds_outlier_details.csv"))

# ============================================================================
# STEP 5: Visualizations
# ============================================================================

cat("STEP 5: Creating visualizations...\n")

# 5A: NMDS plot with outliers highlighted
p_outliers <- ggplot(nmds_full, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site, size = is_outlier, alpha = is_outlier)) +
  scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4)) +
  scale_alpha_manual(values = c("FALSE" = 0.5, "TRUE" = 1)) +
  scale_color_viridis_d() +
  geom_text(data = filter(nmds_full, outlier_rank <= 5),
            aes(label = coral_id), size = 3, nudge_y = 0.05) +
  labs(
    title = "NMDS Ordination with Outliers Highlighted",
    subtitle = paste("Top 10% (n =", n_outliers, ") by distance from centroid"),
    x = "NMDS1", y = "NMDS2"
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "01_nmds_outliers_highlighted.png"),
       p_outliers, width = 10, height = 8, dpi = 300, bg = "white")

# 5B: Outliers colored by CAFI characteristics
p_abundance <- ggplot(nmds_full, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = cafi_abundance, shape = is_outlier),
             size = 3, alpha = 0.7) +
  scale_color_viridis_c(option = "plasma") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +
  labs(
    title = "NMDS Colored by CAFI Abundance",
    subtitle = "Triangles = outliers",
    x = "NMDS1", y = "NMDS2",
    color = "Abundance"
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "02_nmds_by_abundance.png"),
       p_abundance, width = 10, height = 8, dpi = 300, bg = "white")

# 5C: Outliers colored by richness
p_richness <- ggplot(nmds_full, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = cafi_richness, shape = is_outlier),
             size = 3, alpha = 0.7) +
  scale_color_viridis_c(option = "viridis") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17)) +
  labs(
    title = "NMDS Colored by CAFI Richness",
    subtitle = "Triangles = outliers",
    x = "NMDS1", y = "NMDS2",
    color = "Richness"
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "03_nmds_by_richness.png"),
       p_richness, width = 10, height = 8, dpi = 300, bg = "white")

# 5D: Box plots comparing outliers
p_compare <- nmds_full %>%
  select(is_outlier, volume_lab, cafi_abundance, cafi_richness, depth_m) %>%
  pivot_longer(cols = -is_outlier, names_to = "variable", values_to = "value") %>%
  mutate(variable = case_when(
    variable == "volume_lab" ~ "Volume (cm³)",
    variable == "cafi_abundance" ~ "CAFI Abundance",
    variable == "cafi_richness" ~ "CAFI Richness",
    variable == "depth_m" ~ "Depth (m)"
  )) %>%
  ggplot(aes(x = is_outlier, y = value, fill = is_outlier)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~variable, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "red")) +
  labs(
    title = "Characteristics of NMDS Outliers vs Non-outliers",
    x = "Is Outlier",
    y = "Value"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "04_outlier_comparison_boxplots.png"),
       p_compare, width = 12, height = 8, dpi = 300, bg = "white")

# 5E: Composition comparison
composition_compare <- nmds_full %>%
  select(is_outlier, n_crabs, n_shrimp, n_snails, n_fish) %>%
  pivot_longer(cols = -is_outlier, names_to = "taxon", values_to = "count") %>%
  mutate(taxon = gsub("n_", "", taxon)) %>%
  group_by(is_outlier, taxon) %>%
  summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop")

p_composition <- ggplot(composition_compare,
                        aes(x = taxon, y = mean_count, fill = is_outlier)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "red"),
                    labels = c("Non-outlier", "Outlier")) +
  labs(
    title = "Mean CAFI Composition: Outliers vs Non-outliers",
    x = "Taxonomic Group",
    y = "Mean Count per Coral",
    fill = ""
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "05_composition_comparison.png"),
       p_composition, width = 10, height = 6, dpi = 300, bg = "white")

# ============================================================================
# STEP 6: Check for species driving outlier positions
# ============================================================================

cat("\nSTEP 6: Identifying species driving outlier positions...\n\n")

# Get species scores from NMDS
species_scores <- as.data.frame(scores(nmds_bray, display = "species"))
species_scores$species <- rownames(species_scores)
species_scores <- species_scores %>%
  mutate(score_magnitude = sqrt(NMDS1^2 + NMDS2^2)) %>%
  arrange(desc(score_magnitude))

# For each outlier, check which species are overrepresented
outlier_ids <- nmds_full %>% filter(is_outlier) %>% pull(coral_id)

# Get community matrix for outliers
outlier_comm <- community_matrix_clean[rownames(community_matrix_clean) %in% outlier_ids, ]
nonoutlier_comm <- community_matrix_clean[!rownames(community_matrix_clean) %in% outlier_ids, ]

# Compare mean abundances
outlier_means <- colMeans(outlier_comm)
nonoutlier_means <- colMeans(nonoutlier_comm)

species_comparison <- data.frame(
  species = names(outlier_means),
  outlier_mean = outlier_means,
  nonoutlier_mean = nonoutlier_means
) %>%
  mutate(
    fold_change = (outlier_mean + 0.01) / (nonoutlier_mean + 0.01),
    diff = outlier_mean - nonoutlier_mean
  ) %>%
  arrange(desc(abs(diff)))

cat("Species most different between outliers and non-outliers:\n")
print(head(species_comparison, 10))

write_csv(species_comparison,
          file.path(SURVEY_TABLES, "nmds_outlier_species_comparison.csv"))

# ============================================================================
# Summary
# ============================================================================

cat("\n========================================\n")
cat("NMDS Outlier Investigation Summary\n")
cat("========================================\n\n")

cat("Total corals:", nrow(nmds_full), "\n")
cat("Outliers (top 10%):", n_outliers, "\n")
cat("Outlier threshold (centroid distance):", round(outlier_threshold, 3), "\n\n")

cat("Key findings:\n")
cat("  - Volume difference: Wilcoxon p =", round(vol_test$p.value, 3), "\n")
cat("  - Abundance difference: Wilcoxon p =", round(abund_test$p.value, 3), "\n")
cat("  - Richness difference: Wilcoxon p =", round(rich_test$p.value, 3), "\n")
cat("  - Site distribution: Chi-squared p =", round(site_test$p.value, 3), "\n\n")

cat("Outputs saved to:\n")
cat("  - Figures:", fig_dir, "\n")
cat("  - Tables:", SURVEY_TABLES, "\n")

cat("\n✓ NMDS outlier investigation complete\n")
