#!/usr/bin/env Rscript
# ============================================================================
# 07_neighborhood_effects.R - Local Neighborhood Effects Analysis
# ============================================================================
#
# PURPOSE: Test how meter-scale coral neighborhoods affect CAFI colonization
#
# HYPOTHESIS (H3): CAFI abundance varies with local coral density and proximity,
#   reflecting propagule redirection and potential spillover effects.
#
# THEORETICAL PREDICTIONS:
#   - Isolated corals receive MORE larvae per unit area (redirection)
#   - Crowded corals dilute incoming propagules
#   - Facilitation/spillover may counteract dilution at fine scales
#
# KEY METRICS:
#   - Neighbor density: count within 5m
#   - Neighbor volume: total habitat nearby
#   - Isolation index: distance to neighbors
#   - Relative size: focal vs neighbor volumes
#
# MANUSCRIPT FIGURES:
#   >>> FIGURE 6: Neighborhood Effects <<<
#   - Panel A: CAFI density vs neighbor count
#   - Panel B: Isolation effects
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("07: Neighborhood Effects Analysis\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Additional libraries for spatial analysis
suppressPackageStartupMessages({
  library(mgcv)
})

# Load processed data objects
# coral_master already contains: total_cafi, otu_richness, n_crabs, n_shrimps, n_fish, n_snails
# and neighborhood metrics: n_neighbors, mean_neighbor_dist, total_neighbor_volume, etc.
coral_master <- load_object("coral_master.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "neighborhood")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Prepare Neighborhood Data
# ============================================================================

cat("Preparing neighborhood data...\n")

# Use coral_master which already has CAFI metrics and neighborhood data from script 01
# Just need to create derived metrics for analysis
neighborhood_data <- coral_master %>%
  mutate(
    # Ensure CAFI metrics are numeric with 0 for missing
    total_cafi = replace_na(total_cafi, 0),
    species_richness = replace_na(otu_richness, 0),
    shannon = replace_na(shannon_diversity, 0),
    n_crabs = replace_na(n_crabs, 0),
    n_shrimp = replace_na(n_shrimps, 0),  # Note: coral_master uses n_shrimps
    n_fish = replace_na(n_fish, 0),
    n_snails = replace_na(n_snails, 0),

    # Coral size (use 'volume' already calculated in coral_master)
    coral_volume = volume,
    coral_height = height,
    coral_width = width,

    # =========================================
    # LOCAL NEIGHBORHOOD METRICS (meter-scale)
    # =========================================
    # These are already in coral_master from script 01

    # Convert mean distance to meters
    mean_neighbor_dist_m = mean_neighbor_dist / 100,

    # Derived neighborhood indices

    # 1. Local density (neighbors per unit area)
    # Assuming circular neighborhood with radius = mean distance
    local_density = n_neighbors / (pi * (mean_neighbor_dist/100)^2 + 0.01),

    # 2. Crowding index: total neighbor volume relative to distance
    # (already in coral_master, but recalculate for consistency)
    crowding_index = total_neighbor_volume / (mean_neighbor_dist + 1),

    # 3. Isolation index: normalized by focal coral size
    # (already in coral_master, but recalculate for consistency)
    isolation_index = mean_neighbor_dist / (coral_volume^(1/3) + 1),

    # 4. Relative size in neighborhood
    # >1 means focal is larger than avg neighbor
    # (already in coral_master, but recalculate for consistency)
    relative_size = coral_volume / (mean_neighbor_volume + 1),

    # 5. Size asymmetry (competitive asymmetry)
    size_asymmetry = abs(coral_volume - mean_neighbor_volume) /
                     (coral_volume + mean_neighbor_volume + 1),

    # 6. Competition intensity
    # High when many large neighbors are close
    competition_intensity = (n_neighbors * mean_neighbor_volume) /
                           (mean_neighbor_dist^2 + 1),

    # 7. Potential for spillover
    # Neighbors as source of CAFI colonizers
    spillover_potential = total_neighbor_volume / (mean_neighbor_dist + 1),

    # Categorical variables for analysis
    isolation_category = case_when(
      mean_neighbor_dist < quantile(mean_neighbor_dist, 0.33, na.rm = TRUE) ~ "Clustered",
      mean_neighbor_dist > quantile(mean_neighbor_dist, 0.67, na.rm = TRUE) ~ "Isolated",
      TRUE ~ "Intermediate"
    ),

    neighbor_density_cat = case_when(
      n_neighbors == 0 ~ "None",
      n_neighbors <= 2 ~ "Low",
      n_neighbors <= 4 ~ "Medium",
      TRUE ~ "High"
    ),

    # Log transformations
    log_volume = log(coral_volume + 1),
    log_cafi = log(total_cafi + 1),
    log_neighbor_dist = log(mean_neighbor_dist + 1)
  ) %>%
  filter(!is.na(coral_volume), coral_volume > 0,
         !is.na(mean_neighbor_dist))

cat("✓ Prepared", nrow(neighborhood_data), "corals with neighborhood data\n\n")

# ============================================================================
# 1. NEIGHBOR DENSITY EFFECTS
# ============================================================================

cat("1. Analyzing neighbor density effects...\n")

# Model: Does having more neighbors affect CAFI?
density_model <- gam(total_cafi ~
                      s(n_neighbors, k = 5) +
                      s(log_volume) +
                      site,
                    data = neighborhood_data,
                    family = nb(),
                    method = "REML")

# Summary
cat("   Neighbor count effect:\n")
density_summary <- summary(density_model)
print(density_summary$s.table)

# Visualization
p_density <- ggplot(neighborhood_data, aes(x = n_neighbors, y = total_cafi)) +
  geom_jitter(aes(color = site, shape = site), alpha = 0.6, width = 0.2, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              method.args = list(family = "nb"), color = "black", linewidth = 0.8) +
  scale_y_sqrt() +
  scale_color_site() +
  scale_shape_site() +
  labs(title = "Effect of neighbor count on CAFI abundance",
       subtitle = "Local neighborhood within meters",
       x = "Number of neighboring corals",
       y = "Total CAFI (sqrt scale)") +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "neighbor_count_effect.png"),
       p_density, width = 10, height = 7, dpi = 300)

# By density category
density_stats <- neighborhood_data %>%
  group_by(neighbor_density_cat) %>%
  summarise(
    n = n(),
    mean_cafi = mean(total_cafi),
    se_cafi = sd(total_cafi)/sqrt(n()),
    mean_richness = mean(species_richness),
    .groups = "drop"
  )

write_csv(density_stats, file.path(TABLES_DIR, "cafi_by_neighbor_density.csv"))

# ============================================================================
# 2. NEIGHBOR DISTANCE (ISOLATION) EFFECTS
# ============================================================================

cat("\n2. Analyzing isolation/distance effects...\n")

# Model: Does distance to neighbors matter?
distance_model <- gam(total_cafi ~
                       s(mean_neighbor_dist, k = 5) +
                       s(log_volume) +
                       site,
                     data = neighborhood_data,
                     family = nb(),
                     method = "REML")

cat("   Distance effect:\n")
print(summary(distance_model)$s.table)

# Test for isolation effect
isolation_test <- kruskal.test(total_cafi ~ isolation_category,
                               data = neighborhood_data)
cat("   Isolation category effect: p =", format.pval(isolation_test$p.value), "\n")

# Visualization
p_isolation <- neighborhood_data %>%
  ggplot(aes(x = mean_neighbor_dist_m, y = total_cafi)) +
  geom_point(aes(color = site, shape = site), alpha = 0.6, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              method.args = list(family = "nb"), color = "black", linewidth = 0.8) +
  scale_y_sqrt() +
  scale_color_site() +
  scale_shape_site() +
  labs(title = "Effect of neighbor distance on CAFI",
       subtitle = "Distance to nearest neighboring corals",
       x = "Mean distance to neighbors (m)",
       y = "Total CAFI (sqrt scale)") +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "neighbor_distance_effect.png"),
       p_isolation, width = 10, height = 7, dpi = 300)

# Box plot by isolation category
p_isolation_box <- neighborhood_data %>%
  mutate(isolation_category = factor(isolation_category, levels = c("Clustered", "Intermediate", "Isolated"))) %>%
  ggplot(aes(x = isolation_category, y = total_cafi, fill = isolation_category)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = c("Clustered" = "#E69F00", "Intermediate" = "#56B4E9", "Isolated" = "#009E73")) +
  labs(title = "CAFI abundance by isolation level",
       x = NULL,
       y = "Total CAFI") +
  theme_publication() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11, face = "bold"))

ggsave(file.path(fig_dir, "isolation_category_boxplot.png"),
       p_isolation_box, width = 8, height = 6, dpi = 300)

# ============================================================================
# 3. NEIGHBOR VOLUME (COMPETITION/FACILITATION) EFFECTS
# ============================================================================

cat("\n3. Analyzing neighbor volume effects (competition vs facilitation)...\n")

# Does total neighbor volume affect CAFI?
# Positive = facilitation/spillover, Negative = competition

volume_model <- gam(total_cafi ~
                     s(total_neighbor_volume, k = 5) +
                     s(log_volume) +
                     s(mean_neighbor_dist, k = 5) +
                     site,
                   data = neighborhood_data %>% filter(total_neighbor_volume > 0),
                   family = nb(),
                   method = "REML")

cat("   Neighbor volume effect:\n")
print(summary(volume_model)$s.table)

# Visualization
p_neighbor_vol <- neighborhood_data %>%
  filter(total_neighbor_volume > 0) %>%
  ggplot(aes(x = total_neighbor_volume, y = total_cafi)) +
  geom_point(aes(color = site, shape = site), alpha = 0.6, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              method.args = list(family = "nb"), color = "black", linewidth = 0.8) +
  scale_x_log10() +
  scale_y_sqrt() +
  scale_color_site() +
  scale_shape_site() +
  labs(title = "Effect of neighbor volume on CAFI",
       subtitle = "Total volume of neighboring corals - potential spillover source",
       x = "Total neighbor volume (cm³, log scale)",
       y = "Total CAFI (sqrt scale)") +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "neighbor_volume_effect.png"),
       p_neighbor_vol, width = 10, height = 7, dpi = 300)

# ============================================================================
# 4. RELATIVE SIZE IN NEIGHBORHOOD
# ============================================================================

cat("\n4. Analyzing relative size effects...\n")

# Does being larger/smaller than neighbors matter?
relative_model <- gam(total_cafi ~
                       s(relative_size, k = 5) +
                       s(log_volume) +
                       site,
                     data = neighborhood_data %>% filter(relative_size > 0),
                     family = nb(),
                     method = "REML")

# Visualization
p_relative <- neighborhood_data %>%
  filter(relative_size > 0, relative_size < 10) %>%
  ggplot(aes(x = relative_size, y = total_cafi)) +
  geom_point(aes(color = site, shape = site), alpha = 0.6, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              method.args = list(family = "nb"), color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#D55E00", linewidth = 0.6) +
  scale_y_sqrt() +
  scale_color_site() +
  scale_shape_site() +
  labs(title = "Effect of relative size in neighborhood",
       subtitle = "Values >1 indicate focal coral is larger than average neighbor",
       x = "Relative size (focal / mean neighbor)",
       y = "Total CAFI (sqrt scale)") +
  annotate("text", x = 0.3, y = max(neighborhood_data$total_cafi) * 0.9,
           label = "Smaller than\nneighbors", size = 3, fontface = "italic") +
  annotate("text", x = 3, y = max(neighborhood_data$total_cafi) * 0.9,
           label = "Larger than\nneighbors", size = 3, fontface = "italic") +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "relative_size_effect.png"),
       p_relative, width = 10, height = 7, dpi = 300)

# ============================================================================
# 5. INTERACTIVE EFFECTS: SIZE × NEIGHBORHOOD
# ============================================================================

cat("\n5. Analyzing size × neighborhood interactions...\n")

# Does neighborhood context modify coral size effects?
interaction_model <- gam(total_cafi ~
                          te(log_volume, crowding_index, k = c(5, 5)) +
                          site,
                        data = neighborhood_data %>% filter(crowding_index > 0),
                        family = nb(),
                        method = "REML")

# Create prediction surface
vol_range <- range(neighborhood_data$log_volume, na.rm = TRUE)
crowd_range <- quantile(neighborhood_data$crowding_index[neighborhood_data$crowding_index > 0],
                        c(0.05, 0.95), na.rm = TRUE)

pred_grid <- expand.grid(
  log_volume = seq(vol_range[1], vol_range[2], length.out = 50),
  crowding_index = seq(crowd_range[1], crowd_range[2], length.out = 50),
  site = neighborhood_data$site[1]
)

pred_grid$predicted <- predict(interaction_model, pred_grid, type = "response")

# Interaction surface plot
p_interaction <- ggplot(pred_grid, aes(x = log_volume, y = crowding_index,
                                        fill = predicted)) +
  geom_tile() +
  geom_contour(aes(z = predicted), color = "white", alpha = 0.5, bins = 8, linewidth = 0.4) +
  scale_fill_viridis_c(name = "Predicted\nCAFI", trans = "sqrt", option = "viridis") +
  labs(title = "Interactive effect: coral size × neighborhood crowding",
       x = "log(coral volume)",
       y = "Crowding index\n(neighbor volume / distance)") +
  theme_publication() +
  theme(panel.border = element_blank())

ggsave(file.path(fig_dir, "size_crowding_interaction.png"),
       p_interaction, width = 10, height = 8, dpi = 300)

# ============================================================================
# 6. SPILLOVER ANALYSIS
# ============================================================================

cat("\n6. Testing for spillover effects from neighbors...\n")

# Model: Does spillover potential predict CAFI richness?
# (richness may be more sensitive to spillover than abundance)
spillover_model <- gam(species_richness ~
                        s(spillover_potential, k = 5) +
                        s(log_volume) +
                        s(mean_neighbor_dist, k = 5) +
                        site,
                      data = neighborhood_data %>% filter(spillover_potential > 0),
                      family = poisson(),
                      method = "REML")

cat("   Spillover potential effect on richness:\n")
print(summary(spillover_model)$s.table)

# Visualization
p_spillover <- neighborhood_data %>%
  filter(spillover_potential > 0) %>%
  ggplot(aes(x = spillover_potential, y = species_richness)) +
  geom_point(aes(color = site, shape = site), alpha = 0.6, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              method.args = list(family = "poisson"), color = "black", linewidth = 0.8) +
  scale_x_log10() +
  scale_color_site() +
  scale_shape_site() +
  labs(title = "Spillover effect on species richness",
       subtitle = "Do larger nearby neighbors provide colonists?",
       x = "Spillover potential (neighbor volume / distance)",
       y = "Species richness") +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "spillover_richness_effect.png"),
       p_spillover, width = 10, height = 7, dpi = 300)

# ============================================================================
# 7. TAXONOMIC GROUP RESPONSES TO NEIGHBORHOOD
# ============================================================================

cat("\n7. Analyzing taxonomic group responses to neighborhood...\n")

# Different taxa may respond differently to neighborhood context
taxa_neighborhood <- neighborhood_data %>%
  pivot_longer(cols = c(n_crabs, n_shrimp, n_fish, n_snails),
               names_to = "taxon",
               values_to = "count") %>%
  mutate(taxon = str_remove(taxon, "n_") %>% str_to_title())

# Plot taxon-specific responses to crowding
p_taxa_crowding <- taxa_neighborhood %>%
  filter(crowding_index > 0) %>%
  ggplot(aes(x = crowding_index, y = count)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
              method.args = list(family = "nb"), color = "black", linewidth = 0.8) +
  facet_wrap(~taxon, scales = "free_y") +
  scale_x_log10() +
  labs(title = "Taxonomic group responses to neighborhood crowding",
       x = "Crowding index (log scale)",
       y = "Count") +
  theme_publication()

ggsave(file.path(fig_dir, "taxa_crowding_response.png"),
       p_taxa_crowding, width = 12, height = 10, dpi = 300)

# ============================================================================
# 8. COMPREHENSIVE SUMMARY TABLE
# ============================================================================

cat("\n8. Creating summary statistics...\n")

# Summary by isolation and size
neighborhood_summary <- neighborhood_data %>%
  group_by(isolation_category, size_class) %>%
  summarise(
    n_corals = n(),
    mean_cafi = round(mean(total_cafi), 1),
    se_cafi = round(sd(total_cafi)/sqrt(n()), 1),
    mean_richness = round(mean(species_richness), 1),
    mean_neighbors = round(mean(n_neighbors), 1),
    mean_dist_m = round(mean(mean_neighbor_dist_m), 2),
    .groups = "drop"
  )

write_csv(neighborhood_summary,
          file.path(TABLES_DIR, "neighborhood_effects_summary.csv"))

# Model comparison table
model_results <- tibble(
  Effect = c("Neighbor count", "Neighbor distance", "Neighbor volume",
             "Relative size", "Spillover potential"),
  Response = c("Abundance", "Abundance", "Abundance", "Abundance", "Richness"),
  Direction = c("See GAM", "See GAM", "See GAM", "See GAM", "See GAM"),
  Interpretation = c(
    "More neighbors → more CAFI",
    "Non-linear distance effect",
    "Larger neighbors may provide colonists",
    "Relative competitive position matters",
    "Nearby volume predicts colonization"
  )
)

write_csv(model_results, file.path(TABLES_DIR, "neighborhood_model_summary.csv"))

# ============================================================================
# 9. COMBINED VISUALIZATION
# ============================================================================

cat("\n9. Creating combined figure panel...\n")

combined_plot <- (p_density + p_isolation) /
                 (p_neighbor_vol + p_relative) +
  plot_annotation(
    title = "Local neighborhood effects on CAFI communities",
    subtitle = "Meter-scale spatial context influences coral-associated fauna",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30")
    )
  )

ggsave(file.path(fig_dir, "neighborhood_effects_panel.png"),
       combined_plot, width = 16, height = 12, dpi = 300)

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Local Neighborhood Effects Summary\n")
cat("========================================\n\n")

cat("Key Findings:\n\n")

cat("1. NEIGHBOR DENSITY:\n")
cat("   - Corals with more neighbors tend to have",
    ifelse(mean(neighborhood_data$total_cafi[neighborhood_data$n_neighbors > 3]) >
           mean(neighborhood_data$total_cafi[neighborhood_data$n_neighbors <= 1]),
           "MORE", "FEWER"), "CAFI\n")
cat("   - Suggests", ifelse(mean(neighborhood_data$total_cafi[neighborhood_data$n_neighbors > 3]) >
                            mean(neighborhood_data$total_cafi[neighborhood_data$n_neighbors <= 1]),
                            "facilitation/spillover", "competition"), "effects\n\n")

cat("2. ISOLATION:\n")
cat("   - Kruskal-Wallis test: p =", format.pval(isolation_test$p.value), "\n")
cat("   - Mean CAFI by category:\n")
for (cat_name in unique(neighborhood_data$isolation_category)) {
  cat("     -", cat_name, ":",
      round(mean(neighborhood_data$total_cafi[neighborhood_data$isolation_category == cat_name]), 1), "\n")
}

cat("\n3. NEIGHBOR VOLUME:\n")
cat("   - Total neighbor volume shows",
    ifelse(summary(volume_model)$s.table[1, 4] < 0.05, "SIGNIFICANT", "non-significant"),
    "effect on CAFI\n")
cat("   - Interpretation: larger neighbors may serve as source populations\n\n")

cat("4. RELATIVE SIZE:\n")
cat("   - Being larger than neighbors is associated with",
    ifelse(mean(neighborhood_data$total_cafi[neighborhood_data$relative_size > 1], na.rm = TRUE) >
           mean(neighborhood_data$total_cafi[neighborhood_data$relative_size < 1], na.rm = TRUE),
           "MORE", "FEWER"), "CAFI\n\n")

# ============================================================================
# COMPILE STANDARDIZED STATISTICAL RESULTS
# ============================================================================

cat("Compiling standardized statistical results...\n")

stats_results <- init_results_df()

# Isolation test (Kruskal-Wallis)
n_iso <- sum(!is.na(neighborhood_data$isolation_category) & !is.na(neighborhood_data$total_cafi))
k_iso <- n_distinct(neighborhood_data$isolation_category[!is.na(neighborhood_data$isolation_category)])
eps_sq_iso <- as.numeric(isolation_test$statistic) / (n_iso - 1)

stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "H3a",
    question = "Does coral isolation affect CAFI abundance?",
    test_name = "Kruskal-Wallis",
    test_statistic = as.numeric(isolation_test$statistic),
    df = as.character(k_iso - 1),
    p_value = isolation_test$p.value,
    effect_size = eps_sq_iso,
    effect_type = "ε²",
    n = n_iso,
    notes = "Isolation categories: clustered, typical, isolated"
  )
)

# Neighbor density correlation
density_cor <- cor.test(neighborhood_data$n_neighbors, neighborhood_data$total_cafi,
                        method = "spearman", exact = FALSE)
stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "H3b",
    question = "Does neighbor density correlate with CAFI abundance?",
    test_name = "Spearman correlation",
    test_statistic = as.numeric(density_cor$statistic),
    df = as.character(nrow(neighborhood_data) - 2),
    p_value = density_cor$p.value,
    effect_size = as.numeric(density_cor$estimate),
    effect_type = "rho",
    n = nrow(neighborhood_data),
    notes = "Neighbor count vs CAFI abundance"
  )
)

# Crowding index effect (if model exists)
if (exists("volume_model")) {
  gam_summary <- summary(volume_model)
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H3c",
      question = "Does total neighbor volume affect CAFI?",
      test_name = "GAM (smooth term)",
      test_statistic = gam_summary$s.table[1, 3],
      df = round(gam_summary$s.table[1, 2], 1),
      p_value = gam_summary$s.table[1, 4],
      effect_size = gam_summary$dev.expl,
      effect_type = "Deviance explained",
      n = nrow(neighborhood_data),
      notes = "Non-linear effect of neighbor volume"
    )
  )
}

# Save standardized results
save_stats_summary(stats_results, "07_neighborhood_effects", "Neighborhood Effects Analysis")

cat("Output files saved:\n")
cat("  Figures:", fig_dir, "\n")
cat("  Tables:", TABLES_DIR, "\n\n")

cat("✅ Local neighborhood effects analysis complete!\n")
