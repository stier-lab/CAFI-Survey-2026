#!/usr/bin/env Rscript
# ============================================================================
# Fig6_comprehensive_neighborhood_effects.R - PRD Figure 6
#
# Comprehensive visualization of LOCAL (meter-scale) neighborhood effects
# on CAFI communities, testing Hypothesis H3: Propagule Redirection
#
# Key metrics from PRD:
#   1. Neighbor density (count within 5m)
#   2. Neighbor volume (total habitat nearby)
#   3. Isolation index (distance to neighbors)
#   4. Relative size (focal vs. neighbor volumes)
#   5. Spillover potential (larvae from nearby corals)
#
# Expected Results:
#   - Positive isolation effect (propagule redirection)
#   - Possible positive density effect (spillover/facilitation)
#
# Author: CAFI Analysis Pipeline
# Date: 2025-11-23
# ============================================================================

cat("\n========================================\n")
cat("FIGURE 6: Comprehensive Neighborhood Effects\n")
cat("Testing meter-scale local effects (NOT spatial autocorrelation)\n")
cat("========================================\n\n")

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(scales)
})

# Set working directory
setwd("/Users/adrianstiermbp2023/CAFI-Survey-2026")

# Load data
cat("Loading and preparing data...\n")
cafi_data <- read.csv("data/1. survey_cafi_data_w_taxonomy_summer2019_v5.csv")
coral_data <- read.csv("data/1. survey_coral_characteristics_merged_v2.csv")

# Aggregate CAFI
cafi_summary <- cafi_data %>%
  filter(!is.na(genus) & genus != "" & genus != "NA") %>%
  group_by(coral_id) %>%
  summarise(
    cafi_abundance = n(),
    cafi_richness = n_distinct(paste(genus, species)),
    shannon = vegan::diversity(table(paste(genus, species))),
    .groups = "drop"
  )

# Merge and calculate neighborhood metrics
neighborhood_data <- coral_data %>%
  left_join(cafi_summary, by = "coral_id") %>%
  mutate(
    # Replace NA with 0
    cafi_abundance = replace_na(cafi_abundance, 0),
    cafi_richness = replace_na(cafi_richness, 0),
    shannon = replace_na(shannon, 0),

    # Coral volume
    volume = coalesce(volume_field, volume_lab,
                     length_field * width_field * height_field),

    # ==========================================
    # KEY NEIGHBORHOOD METRICS FROM PRD
    # ==========================================

    # 1. NEIGHBOR DENSITY (count within 5m)
    neighbor_density = number_of_neighbors,

    # 2. TOTAL NEIGHBOR VOLUME (habitat amount nearby)
    neighbor_volume = coalesce(combined_total_volume_of_neighbors, 0),

    # 3. ISOLATION INDEX (normalized by coral size)
    mean_dist = coalesce(mean_neighbor_distance, 500) / 100,  # Convert to meters
    isolation_index = mean_dist / (volume^(1/3) + 1),

    # 4. RELATIVE SIZE (focal vs neighbors)
    mean_neighbor_vol = coalesce(mean_total_volume_of_neighbors, 1),
    relative_size = volume / (mean_neighbor_vol + 1),

    # 5. SPILLOVER POTENTIAL (volume/distance)
    spillover_potential = neighbor_volume / (mean_dist * 100 + 1),

    # Additional metrics
    crowding_index = neighbor_volume / (mean_dist + 0.1),

    # Log transforms
    log_volume = log(volume + 1),
    log_abundance = log(cafi_abundance + 1),

    # Categories
    density_category = cut(neighbor_density,
                          breaks = c(-0.1, 0, 2, 4, 100),
                          labels = c("Isolated", "Low", "Medium", "High")),

    isolation_category = cut(isolation_index,
                            breaks = quantile(isolation_index, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                            labels = c("Clustered", "Intermediate", "Isolated"),
                            include.lowest = TRUE),

    # Clean site
    site = str_extract(site, "^[A-Z]+")
  ) %>%
  filter(!is.na(volume), volume > 0)

cat(sprintf("✓ Data prepared: %d corals with neighborhood metrics\n\n", nrow(neighborhood_data)))

# ============================================================================
# CREATE COMPREHENSIVE FIGURE 6
# ============================================================================

cat("Creating comprehensive neighborhood effects figure...\n")

# Theme
theme_publication <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "grey30"),
      axis.title = element_text(face = "bold", size = 10)
    )
}

# Colors
site_colors <- c("HAU" = "#E69F00", "MAT" = "#56B4E9", "MRB" = "#009E73")

# ============================================================================
# Panel A: NEIGHBOR DENSITY EFFECT
# ============================================================================

# Fit GAM for smooth relationship
gam_density <- gam(cafi_abundance ~ s(neighbor_density, k = 5) + s(log_volume),
                  data = neighborhood_data, family = nb())

# Create prediction data
pred_density <- data.frame(
  neighbor_density = seq(0, max(neighborhood_data$neighbor_density, na.rm = TRUE), length.out = 100),
  log_volume = mean(neighborhood_data$log_volume, na.rm = TRUE)
)
pred_density$predicted <- predict(gam_density, newdata = pred_density, type = "response")

p_density <- ggplot(neighborhood_data, aes(x = neighbor_density)) +
  # Raw data
  geom_jitter(aes(y = cafi_abundance, color = site),
             alpha = 0.4, width = 0.2, size = 2) +
  # GAM smooth
  geom_line(data = pred_density, aes(y = predicted),
           color = "black", linewidth = 1.2) +
  # Confidence band
  geom_smooth(aes(y = cafi_abundance), method = "gam",
             formula = y ~ s(x, k = 5),
             se = TRUE, color = "black", fill = "grey30", alpha = 0.2) +
  scale_color_manual(values = site_colors) +
  scale_y_sqrt(breaks = c(0, 10, 25, 50, 100)) +
  labs(
    title = "A. Neighbor Density Effect",
    subtitle = "Count of corals within 5 meters",
    x = "Number of Neighboring Corals",
    y = "CAFI Abundance (sqrt scale)",
    color = "Site"
  ) +
  theme_publication()

# ============================================================================
# Panel B: NEIGHBOR VOLUME EFFECT
# ============================================================================

p_volume <- neighborhood_data %>%
  filter(neighbor_volume > 0) %>%
  ggplot(aes(x = neighbor_volume, y = cafi_abundance)) +
  geom_point(aes(color = site), alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x,
             color = "black", linewidth = 1.2, se = TRUE) +
  scale_x_log10(breaks = c(100, 1000, 10000, 100000)) +
  scale_y_sqrt(breaks = c(0, 10, 25, 50, 100)) +
  scale_color_manual(values = site_colors) +
  annotation_logticks(sides = "b") +
  labs(
    title = "B. Neighbor Volume Effect",
    subtitle = "Total habitat amount nearby",
    x = "Total Neighbor Volume (cm³, log scale)",
    y = "CAFI Abundance (sqrt scale)",
    color = "Site"
  ) +
  theme_publication()

# ============================================================================
# Panel C: ISOLATION INDEX EFFECT
# ============================================================================

p_isolation <- ggplot(neighborhood_data, aes(x = isolation_index, y = cafi_abundance)) +
  geom_point(aes(color = site, size = volume), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
             color = "black", linewidth = 1.2, se = TRUE) +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  scale_color_manual(values = site_colors) +
  scale_y_sqrt(breaks = c(0, 10, 25, 50, 100)) +
  labs(
    title = "C. Isolation Index Effect",
    subtitle = "Distance normalized by coral size",
    x = "Isolation Index (higher = more isolated)",
    y = "CAFI Abundance (sqrt scale)",
    color = "Site"
  ) +
  theme_publication()

# ============================================================================
# Panel D: RELATIVE SIZE EFFECT
# ============================================================================

p_relative <- neighborhood_data %>%
  filter(relative_size < 10) %>%  # Remove extreme outliers
  ggplot(aes(x = relative_size, y = cafi_abundance)) +
  geom_point(aes(color = site), alpha = 0.5, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
             color = "black", linewidth = 1.2, se = TRUE) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5)) +
  scale_y_sqrt(breaks = c(0, 10, 25, 50, 100)) +
  scale_color_manual(values = site_colors) +
  annotation_logticks(sides = "b") +
  labs(
    title = "D. Relative Size Effect",
    subtitle = "Focal coral size vs. neighbors",
    x = "Relative Size (focal/neighbor, log scale)",
    y = "CAFI Abundance (sqrt scale)",
    color = "Site"
  ) +
  theme_publication()

# ============================================================================
# Panel E: SPILLOVER POTENTIAL
# ============================================================================

p_spillover <- neighborhood_data %>%
  filter(spillover_potential > 0) %>%
  ggplot(aes(x = spillover_potential, y = cafi_abundance)) +
  geom_point(aes(color = site), alpha = 0.5, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
             color = "black", linewidth = 1.2, se = TRUE) +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) +
  scale_y_sqrt(breaks = c(0, 10, 25, 50, 100)) +
  scale_color_manual(values = site_colors) +
  annotation_logticks(sides = "b") +
  labs(
    title = "E. Spillover Potential",
    subtitle = "Neighbor volume / distance",
    x = "Spillover Potential (log scale)",
    y = "CAFI Abundance (sqrt scale)",
    color = "Site"
  ) +
  theme_publication()

# ============================================================================
# Panel F: SYNTHESIS - Categories
# ============================================================================

# Summary by categories
category_summary <- neighborhood_data %>%
  group_by(density_category) %>%
  summarise(
    n = n(),
    mean_cafi = mean(cafi_abundance),
    se_cafi = sd(cafi_abundance) / sqrt(n()),
    mean_richness = mean(cafi_richness),
    se_richness = sd(cafi_richness) / sqrt(n()),
    .groups = "drop"
  )

p_synthesis <- ggplot(category_summary, aes(x = density_category)) +
  geom_col(aes(y = mean_cafi), fill = "steelblue", alpha = 0.7, width = 0.7) +
  geom_errorbar(aes(ymin = mean_cafi - se_cafi,
                   ymax = mean_cafi + se_cafi),
               width = 0.2, linewidth = 0.8) +
  geom_text(aes(y = mean_cafi + se_cafi + 2, label = paste0("n=", n)),
           size = 3.5) +
  labs(
    title = "F. Synthesis: Density Categories",
    subtitle = "Mean CAFI abundance by neighbor density",
    x = "Neighborhood Density Category",
    y = "Mean CAFI Abundance (±SE)"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# ============================================================================
# COMBINE ALL PANELS
# ============================================================================

final_figure <- (p_density | p_volume | p_isolation) /
                (p_relative | p_spillover | p_synthesis) +
  plot_annotation(
    title = "Figure 6: Comprehensive LOCAL Neighborhood Effects (Meter-Scale)",
    subtitle = "Testing H3: Propagule redirection and spillover effects within 5 meters of focal corals",
    caption = sprintf("n = %d corals | GAM smoothers with 95%% CI | Local effects only (NOT spatial autocorrelation)",
                     nrow(neighborhood_data)),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "grey30"),
      plot.caption = element_text(size = 10, hjust = 1)
    )
  )

# Save figure
dir.create("output/figures/hypothesis_tests", showWarnings = FALSE, recursive = TRUE)
ggsave("output/figures/hypothesis_tests/Fig6_neighborhood_effects_comprehensive.png",
       final_figure,
       width = 18, height = 12, dpi = 300, bg = "white")

cat("   ✓ Figure saved: Fig6_neighborhood_effects_comprehensive.png\n\n")

# ============================================================================
# STATISTICAL TESTS FOR MANUSCRIPT
# ============================================================================

cat("========================================\n")
cat("STATISTICAL RESULTS FOR MANUSCRIPT\n")
cat("========================================\n\n")

# Test each metric
cat("1. Neighbor Density Effect:\n")
m_density <- glm(cafi_abundance ~ neighbor_density + log_volume,
                data = neighborhood_data, family = "poisson")
density_coef <- coef(m_density)[2]
density_p <- summary(m_density)$coefficients[2, 4]
cat(sprintf("   β = %.3f, p = %.4f\n", density_coef, density_p))
if (density_coef > 0 && density_p < 0.05) {
  cat("   ✓ Positive density effect (possible spillover)\n\n")
} else {
  cat("   ✗ No significant positive density effect\n\n")
}

cat("2. Neighbor Volume Effect:\n")
volume_data <- neighborhood_data %>% filter(neighbor_volume > 0)
m_volume <- glm(cafi_abundance ~ log(neighbor_volume + 1) + log_volume,
               data = volume_data, family = "poisson")
volume_coef <- coef(m_volume)[2]
volume_p <- summary(m_volume)$coefficients[2, 4]
cat(sprintf("   β = %.3f, p = %.4f\n", volume_coef, volume_p))

cat("\n3. Isolation Index Effect:\n")
m_isolation <- glm(cafi_abundance ~ isolation_index + log_volume,
                  data = neighborhood_data, family = "poisson")
isolation_coef <- coef(m_isolation)[2]
isolation_p <- summary(m_isolation)$coefficients[2, 4]
cat(sprintf("   β = %.3f, p = %.4f\n", isolation_coef, isolation_p))
if (isolation_coef > 0 && isolation_p < 0.05) {
  cat("   ✓ Positive isolation effect (propagule redirection)\n\n")
} else {
  cat("   ✗ No significant isolation effect\n\n")
}

cat("4. Relative Size Effect:\n")
m_relative <- glm(cafi_abundance ~ log(relative_size + 0.1) + log_volume,
                 data = neighborhood_data, family = "poisson")
relative_coef <- coef(m_relative)[2]
relative_p <- summary(m_relative)$coefficients[2, 4]
cat(sprintf("   β = %.3f, p = %.4f\n", relative_coef, relative_p))

cat("\n5. Spillover Potential:\n")
spillover_data <- neighborhood_data %>% filter(spillover_potential > 0)
m_spillover <- glm(cafi_abundance ~ log(spillover_potential + 1) + log_volume,
                  data = spillover_data, family = "poisson")
spillover_coef <- coef(m_spillover)[2]
spillover_p <- summary(m_spillover)$coefficients[2, 4]
cat(sprintf("   β = %.3f, p = %.4f\n", spillover_coef, spillover_p))

# Summary for H3
cat("\n========================================\n")
cat("HYPOTHESIS H3 SUMMARY:\n")
cat("========================================\n")

if ((density_coef > 0 && density_p < 0.05) ||
    (isolation_coef > 0 && isolation_p < 0.05) ||
    (spillover_coef > 0 && spillover_p < 0.05)) {
  cat("✓ H3 PARTIALLY SUPPORTED: Local neighborhood effects detected\n")
} else {
  cat("✗ H3 NOT STRONGLY SUPPORTED: Weak local neighborhood effects\n")
}

cat("\nKey Finding: Effects are LOCAL (within meters), not broad spatial patterns\n")
cat("This is about propagule interception and spillover, NOT autocorrelation\n")

# Save results
results <- list(
  density_effect = density_coef,
  volume_effect = volume_coef,
  isolation_effect = isolation_coef,
  relative_size_effect = relative_coef,
  spillover_effect = spillover_coef,
  p_values = c(density_p, volume_p, isolation_p, relative_p, spillover_p)
)

saveRDS(results, "output/objects/H3_neighborhood_results.rds")
cat("\n✓ Analysis complete. Results saved.\n")