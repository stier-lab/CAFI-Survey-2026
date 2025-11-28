#!/usr/bin/env Rscript
# ============================================================================
# create_nmds_figures.R - Publication NMDS figures with species vectors
#
# Purpose: Create main text NMDS figure excluding depauperate corals,
# with species vectors showing the top species distinguishing sites.
# Also create supplementary figure documenting excluded outliers.
#
# Author: CAFI Analysis Pipeline
# Date: 2025-11-22
# ============================================================================

cat("\n========================================\n")
cat("Creating Publication NMDS Figures\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/00_load_libraries.R"))

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))

# Create output directories
fig_dir <- file.path(SURVEY_FIGURES, "diversity")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# STEP 1: Identify and exclude depauperate corals
# ============================================================================

cat("STEP 1: Identifying depauperate corals to exclude...\n")

# Calculate CAFI metrics for each coral
cafi_summary <- cafi_clean %>%
  group_by(coral_id) %>%
  summarise(
    cafi_abundance = n(),
    cafi_richness = n_distinct(species),
    .groups = "drop"
  )

# Define exclusion criteria: corals with very low richness (<=3 species)
# These represent depauperate communities that obscure main patterns
min_richness_threshold <- 3

excluded_corals <- cafi_summary %>%
  filter(cafi_richness <= min_richness_threshold) %>%
  pull(coral_id)

n_excluded <- length(excluded_corals)
n_total <- nrow(cafi_summary)

cat(sprintf("  Total corals: %d\n", n_total))
cat(sprintf("  Excluded (≤%d species): %d (%.1f%%)\n",
            min_richness_threshold, n_excluded, 100 * n_excluded / n_total))
cat(sprintf("  Retained: %d\n\n", n_total - n_excluded))

# Get details on excluded corals
excluded_details <- cafi_summary %>%
  filter(coral_id %in% excluded_corals) %>%
  left_join(metadata, by = "coral_id") %>%
  left_join(survey_master %>% select(coral_id, volume_lab), by = "coral_id") %>%
  arrange(cafi_richness, cafi_abundance)

cat("Excluded coral summary:\n")
print(excluded_details %>%
        group_by(site) %>%
        summarise(n = n(),
                  mean_richness = mean(cafi_richness),
                  mean_abundance = mean(cafi_abundance),
                  .groups = "drop"))

# Save excluded coral details
write_csv(excluded_details,
          file.path(SURVEY_TABLES, "nmds_excluded_depauperate_corals.csv"))

# ============================================================================
# STEP 2: Create filtered community matrix and run NMDS
# ============================================================================

cat("\nSTEP 2: Running NMDS on filtered data...\n")

# Filter community matrix
community_filtered <- community_matrix[!rownames(community_matrix) %in% excluded_corals, ]
community_filtered <- community_filtered[rowSums(community_filtered) > 0, ]
community_filtered <- community_filtered[, colSums(community_filtered) > 0]

cat(sprintf("  Filtered matrix: %d corals × %d species\n",
            nrow(community_filtered), ncol(community_filtered)))

# Run NMDS
set.seed(123)
nmds_filtered <- metaMDS(community_filtered, distance = "bray", k = 2,
                         trymax = 100, autotransform = FALSE)

cat(sprintf("  NMDS stress: %.3f\n", nmds_filtered$stress))

# Extract site scores
site_scores <- as.data.frame(scores(nmds_filtered, display = "sites"))
site_scores$coral_id <- rownames(site_scores)
site_scores <- site_scores %>%
  left_join(metadata, by = "coral_id") %>%
  left_join(cafi_summary, by = "coral_id")

# Extract species scores
species_scores <- as.data.frame(scores(nmds_filtered, display = "species"))
species_scores$species <- rownames(species_scores)

# ============================================================================
# STEP 3: Identify top species distinguishing sites
# ============================================================================

cat("\nSTEP 3: Identifying species that distinguish sites...\n")

# Calculate species importance (distance from origin in NMDS space)
species_scores <- species_scores %>%
  mutate(
    importance = sqrt(NMDS1^2 + NMDS2^2),
    # Clean species names for display
    species_label = gsub("sp_", "", species),
    species_label = gsub("_", " ", species_label)
  ) %>%
  arrange(desc(importance))

# Select top 15 species
top_species <- species_scores %>%
  slice_head(n = 15)

cat("Top 15 species by NMDS importance:\n")
print(top_species %>% select(species_label, NMDS1, NMDS2, importance))

# ============================================================================
# STEP 4: Create main text NMDS figure
# ============================================================================

cat("\nSTEP 4: Creating main text NMDS figure...\n")

# Define site colors
site_colors <- c("HAU" = "#440154", "MAT" = "#21918c", "MRB" = "#fde725")

# Calculate centroids for site labels
site_centroids <- site_scores %>%
  group_by(site) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  )

# Main NMDS plot with species vectors
p_nmds_main <- ggplot() +
  # Site ellipses (95% confidence)
  stat_ellipse(data = site_scores,
               aes(x = NMDS1, y = NMDS2, color = site),
               level = 0.95, linewidth = 1, linetype = "dashed") +
  # Points colored by site
  geom_point(data = site_scores,
             aes(x = NMDS1, y = NMDS2, color = site, size = cafi_richness),
             alpha = 0.6) +
  # Species vectors
  geom_segment(data = top_species,
               aes(x = 0, y = 0, xend = NMDS1 * 1.5, yend = NMDS2 * 1.5),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray40", linewidth = 0.5, alpha = 0.7) +
  # Species labels
  geom_text(data = top_species,
            aes(x = NMDS1 * 1.7, y = NMDS2 * 1.7, label = species_label),
            size = 2.5, color = "gray20", fontface = "italic") +
  # Site labels at centroids
  geom_label(data = site_centroids,
             aes(x = NMDS1, y = NMDS2, label = site, fill = site),
             color = "white", fontface = "bold", size = 4,
             label.padding = unit(0.3, "lines")) +
  # Styling
  scale_color_manual(values = site_colors, guide = "none") +
  scale_fill_manual(values = site_colors, guide = "none") +
  scale_size_continuous(range = c(1, 5), name = "Species\nrichness") +
  labs(
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.9, 0.15),
    legend.background = element_rect(fill = "white", color = "gray80"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  # Add stress annotation

  annotate("text", x = Inf, y = -Inf,
           label = sprintf("Stress = %.3f", nmds_filtered$stress),
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "gray40")

ggsave(file.path(fig_dir, "nmds_ordination_cleaned.png"),
       p_nmds_main, width = 10, height = 8, dpi = 300, bg = "white")

cat("  ✓ Saved: nmds_ordination_cleaned.png\n")

# ============================================================================
# STEP 5: Create supplementary figure showing excluded corals
# ============================================================================

cat("\nSTEP 5: Creating supplementary figure for excluded corals...\n")

# Run NMDS on full dataset for comparison
community_full <- community_matrix[rowSums(community_matrix) > 0, ]
community_full <- community_full[, colSums(community_full) > 0]

set.seed(123)
nmds_full <- metaMDS(community_full, distance = "bray", k = 2,
                     trymax = 100, autotransform = FALSE)

# Extract scores
full_scores <- as.data.frame(scores(nmds_full, display = "sites"))
full_scores$coral_id <- rownames(full_scores)
full_scores <- full_scores %>%
  left_join(metadata, by = "coral_id") %>%
  left_join(cafi_summary, by = "coral_id") %>%
  mutate(excluded = coral_id %in% excluded_corals)

# Supplementary figure showing what was excluded
p_supp <- ggplot(full_scores, aes(x = NMDS1, y = NMDS2)) +
  # All points, highlight excluded
  geom_point(aes(color = excluded, size = cafi_richness, shape = excluded),
             alpha = 0.7) +
  # Labels for excluded points
  geom_text(data = filter(full_scores, excluded),
            aes(label = coral_id),
            size = 2, nudge_y = 0.08, color = "red") +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red"),
                     labels = c("Retained", "Excluded"),
                     name = "") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("Retained", "Excluded"),
                     name = "") +
  scale_size_continuous(range = c(1, 5), name = "Species\nrichness") +
  labs(
    title = "NMDS with Depauperate Corals Highlighted",
    subtitle = sprintf("Excluded: %d corals with ≤%d species",
                       n_excluded, min_richness_threshold),
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("Stress = %.3f", nmds_full$stress),
           hjust = 1.1, vjust = -0.5, size = 3.5, color = "gray40")

ggsave(file.path(fig_dir, "nmds_excluded_corals_supplement.png"),
       p_supp, width = 10, height = 8, dpi = 300, bg = "white")

cat("  ✓ Saved: nmds_excluded_corals_supplement.png\n")

# ============================================================================
# STEP 6: Create comparison of excluded vs retained
# ============================================================================

cat("\nSTEP 6: Creating comparison statistics...\n")

# Comparison statistics
comparison_stats <- full_scores %>%
  group_by(excluded) %>%
  summarise(
    n = n(),
    mean_richness = mean(cafi_richness),
    sd_richness = sd(cafi_richness),
    mean_abundance = mean(cafi_abundance),
    sd_abundance = sd(cafi_abundance),
    .groups = "drop"
  ) %>%
  mutate(
    group = ifelse(excluded, "Excluded", "Retained")
  )

cat("\nComparison of excluded vs retained corals:\n")
print(comparison_stats)

# Boxplot comparison
p_compare <- full_scores %>%
  mutate(group = ifelse(excluded, "Excluded", "Retained")) %>%
  select(group, cafi_richness, cafi_abundance) %>%
  pivot_longer(cols = c(cafi_richness, cafi_abundance),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = case_when(
    metric == "cafi_richness" ~ "Species Richness",
    metric == "cafi_abundance" ~ "Total Abundance"
  )) %>%
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c("Excluded" = "red", "Retained" = "steelblue")) +
  labs(
    title = "CAFI Community Characteristics",
    subtitle = "Excluded corals have depauperate communities",
    x = "",
    y = "Value"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "nmds_excluded_comparison.png"),
       p_compare, width = 8, height = 5, dpi = 300, bg = "white")

cat("  ✓ Saved: nmds_excluded_comparison.png\n")

# Save comparison statistics
write_csv(comparison_stats,
          file.path(SURVEY_TABLES, "nmds_excluded_vs_retained_stats.csv"))

# ============================================================================
# Summary
# ============================================================================

cat("\n========================================\n")
cat("NMDS Figure Generation Complete\n")
cat("========================================\n\n")

cat("Main text figure:\n")
cat(sprintf("  - %d corals (excluding %d with ≤%d species)\n",
            n_total - n_excluded, n_excluded, min_richness_threshold))
cat(sprintf("  - Stress = %.3f\n", nmds_filtered$stress))
cat("  - Top 15 species vectors shown\n\n")

cat("Outputs:\n")
cat("  - Main: nmds_ordination_cleaned.png\n")
cat("  - Supplement: nmds_excluded_corals_supplement.png\n")
cat("  - Supplement: nmds_excluded_comparison.png\n")
cat("  - Tables: nmds_excluded_depauperate_corals.csv\n")
cat("           nmds_excluded_vs_retained_stats.csv\n\n")

cat("✓ NMDS figure generation complete\n")
