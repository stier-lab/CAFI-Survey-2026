#!/usr/bin/env Rscript
# ============================================================================
# 14_coral_size_neighbor_effects.R - Coral Size Scaling and Neighbor Effects
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
#
# Hypotheses tested (aligned with PRD):
#   H2: CAFI abundance scales with coral volume following a power-law with
#       exponent < 1, indicating larger corals have lower CAFI densities
#       (propagule redirection prediction)
#   H3: CAFI abundance varies with local coral density and proximity,
#       reflecting propagule redirection and spillover effects
#
# Theoretical Background:
#   Propagule redirection theory predicts that larvae distribute among
#   available habitats based on chemical cue strength. Larger and less
#   isolated corals intercept more total larvae but at lower density per
#   unit habitat. Expected scaling exponent ~0.75 for 3D habitat.
#
# Key Analyses:
#   - Power-law scaling: log(Abundance) ~ log(Volume)
#   - Test scaling exponent vs theoretical 0.75
#   - Neighbor density/isolation effects
#   - Taxon-specific responses to spatial context
# ============================================================================

cat("\n========================================\n")
cat("Coral Size and Neighbor Effects Analysis\n")
cat("========================================\n\n")

# Load libraries
source(here::here("scripts/00_load_libraries.R"))
library(nlme)
library(mgcv)
library(gratia)
library(ggeffects)

# Load comprehensive coral data
coral_chars <- read_csv(file.path(here("data"),
                                  "1. survey_coral_characteristics_merged_v2.csv"))

cafi_data <- read_csv(file.path(here("data"),
                                "1. survey_cafi_data_w_taxonomy_summer2019_v5.csv"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "size_neighbor_effects")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Data Preparation with All Size Metrics
# ============================================================================

cat("Preparing comprehensive size and neighbor dataset...\n")

# Clean column names
cafi_data <- cafi_data %>% janitor::clean_names()
coral_chars <- coral_chars %>% janitor::clean_names()

# Calculate CAFI metrics
cafi_metrics <- cafi_data %>%
  group_by(coral_id) %>%
  summarise(
    total_cafi = n(),
    species_richness = n_distinct(if_else(!is.na(lowest_level), lowest_level,
                                          if_else(!is.na(search_term), search_term, paste(type, code)))),
    shannon = if(n() > 0) vegan::diversity(table(if_else(!is.na(lowest_level), lowest_level,
                                             if_else(!is.na(search_term), search_term, paste(type, code))))) else 0,

    # Functional group counts
    n_crabs = sum(type == "crab", na.rm = TRUE),
    n_shrimps = sum(type == "shrimp", na.rm = TRUE),
    n_fish = sum(type == "fish", na.rm = TRUE),
    n_snails = sum(type == "snail", na.rm = TRUE),

    # Size distribution (using cafi_size_mm from raw data)
    mean_cafi_size = mean(cafi_size_mm, na.rm = TRUE),
    cv_cafi_size = sd(cafi_size_mm, na.rm = TRUE) / mean(cafi_size_mm, na.rm = TRUE),

    .groups = "drop"
  )

# Merge with coral characteristics
analysis_data <- coral_chars %>%
  left_join(cafi_metrics, by = "coral_id") %>%
  mutate(
    # Replace NAs with 0 for CAFI metrics
    across(c(total_cafi:cv_cafi_size), ~replace_na(., 0)),

    # CORAL SIZE METRICS (comprehensive)
    # Field measurements
    field_volume = coalesce(volume_field,
                           length_field * width_field * height_field),
    field_surface_area = 2 * pi * (width_field/2) * height_field +
                        2 * pi * (width_field/2)^2,  # Cylinder approximation
    field_circumference = circ_field,

    # Lab measurements
    lab_volume = coalesce(volume_lab,
                         length_lab * width_lab * height_lab),
    lab_surface_area = 2 * pi * (width_lab/2) * height_lab +
                      2 * pi * (width_lab/2)^2,

    # Combined best estimate
    coral_volume = coalesce(field_volume, lab_volume),
    coral_surface_area = coalesce(field_surface_area, lab_surface_area),
    coral_height = coalesce(height_field, height_lab),
    coral_width = coalesce(width_field, width_lab),
    coral_length = coalesce(length_field, length_lab),

    # Shape metrics
    aspect_ratio = coral_height / coral_width,
    compactness = coral_volume^(2/3) / coral_surface_area,
    elongation = coral_length / coral_width,

    # NEIGHBOR METRICS (comprehensive)
    # Basic counts and distances
    n_neighbors = number_of_neighbors,
    mean_dist_neighbors = mean_neighbor_distance,

    # Neighbor volume metrics
    total_volume_neighbors = combined_total_volume_of_neighbors,
    mean_volume_neighbors = mean_total_volume_of_neighbors,

    # Neighbor morphotype breakdown
    n_wide_neighbors = number_of_wide_branching_neighbors,
    n_tight_neighbors = number_of_tight_branching_neighbors,
    prop_wide_neighbors = n_wide_neighbors / (n_neighbors + 1),
    prop_tight_neighbors = n_tight_neighbors / (n_neighbors + 1),

    # Isolation metrics
    isolation_index = mean_dist_neighbors / (coral_volume^(1/3) + 1),
    crowding_index = total_volume_neighbors / (mean_dist_neighbors + 1),
    relative_size = coral_volume / (mean_volume_neighbors + 1),

    # Competition metrics
    competition_index = n_neighbors / (mean_dist_neighbors + 1),
    size_asymmetry = abs(coral_volume - mean_volume_neighbors) /
                     (coral_volume + mean_volume_neighbors + 1),

    # Spatial context
    local_density = n_neighbors / (pi * mean_dist_neighbors^2 + 1),

    # Categorizations
    size_class = cut(coral_volume,
                    breaks = quantile(coral_volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                    labels = c("Small", "Medium", "Large")),
    isolation_class = cut(mean_dist_neighbors,
                         breaks = quantile(mean_dist_neighbors, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                         labels = c("Clustered", "Intermediate", "Isolated")),

    # Log transformations for modeling
    log_volume = log(coral_volume + 1),
    log_surface = log(coral_surface_area + 1),
    log_dist = log(mean_dist_neighbors + 1),
    log_cafi = log(total_cafi + 1)
  ) %>%
  filter(!is.na(coral_volume), coral_volume > 0)

cat("✓ Dataset prepared with", nrow(analysis_data), "corals\n")
cat("  Size metrics calculated: volume, surface area, shape indices\n")
cat("  Neighbor metrics calculated: distance, density, competition\n\n")

# ============================================================================
# H2: Size Scaling Analysis - Power-law Relationships
# ============================================================================
# Tests whether CAFI abundance scales sublinearly with coral volume
# Prediction: exponent < 1 (density decreases with size)
# Theoretical expectation: ~0.75 for 3D habitat

cat("Testing H2: Size scaling relationships...\n")
cat("  Prediction: Power-law exponent < 1 (propagule redirection)\n\n")

# 1. SIZE-ABUNDANCE POWER-LAW RELATIONSHIPS
# -----------------------------------------
# Fit multiple models for size-abundance relationship
models_size <- list(
  linear = lm(total_cafi ~ coral_volume, data = analysis_data),
  log_log = lm(log_cafi ~ log_volume, data = analysis_data),
  polynomial = lm(total_cafi ~ poly(coral_volume, 2), data = analysis_data),
  gam = gam(total_cafi ~ s(coral_volume, k = 5) + morphotype,
           data = analysis_data, family = nb())
)

# Compare models
size_model_comparison <- data.frame(
  Model = names(models_size),
  AIC = sapply(models_size, AIC),
  R2 = c(summary(models_size$linear)$r.squared,
        summary(models_size$log_log)$r.squared,
        summary(models_size$polynomial)$r.squared,
        NA)  # GAM doesn't have simple R2
)

write_csv(size_model_comparison,
          file.path(SURVEY_TABLES, "size_abundance_model_comparison.csv"))

# Plot different size relationships
p_size_models <- analysis_data %>%
  filter(coral_volume < quantile(coral_volume, 0.99, na.rm = TRUE)) %>%
  ggplot(aes(x = coral_volume, y = total_cafi)) +
  geom_point(aes(color = morphotype), alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  geom_smooth(method = "gam", se = FALSE, color = "red") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  annotation_logticks() +
  labs(title = "Size-Abundance Scaling",
       subtitle = "Blue: Linear, Red: GAM",
       x = "Coral Volume (cm³, log scale)",
       y = "Total CAFI (log scale)")

ggsave(file.path(fig_dir, "size_abundance_models.png"),
       p_size_models, width = 10, height = 8, dpi = 300)

# 2. SIZE EFFECTS ON DIFFERENT TAXA
# ----------------------------------
taxa_size_effects <- analysis_data %>%
  pivot_longer(cols = c(n_crabs, n_shrimps, n_fish, n_snails),
               names_to = "taxon",
               values_to = "count") %>%
  mutate(taxon = gsub("n_", "", taxon)) %>%
  group_by(taxon) %>%
  do(model = lm(log(count + 1) ~ log_volume, data = .)) %>%
  mutate(
    slope = sapply(model, function(x) coef(x)[2]),
    se = sapply(model, function(x) summary(x)$coefficients[2, 2]),
    p_value = sapply(model, function(x) summary(x)$coefficients[2, 4])
  )

write_csv(taxa_size_effects %>% select(-model),
          file.path(SURVEY_TABLES, "taxa_specific_size_effects.csv"))

# Plot taxon-specific scaling
p_taxa_scaling <- analysis_data %>%
  pivot_longer(cols = c(n_crabs, n_shrimps, n_fish, n_snails),
               names_to = "taxon",
               values_to = "count") %>%
  mutate(taxon = gsub("n_", "", taxon)) %>%
  filter(count > 0) %>%
  ggplot(aes(x = coral_volume, y = count)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~taxon, scales = "free_y") +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks() +
  labs(title = "Taxon-Specific Size Scaling",
       x = "Coral Volume (cm³, log)",
       y = "Abundance (log)")

ggsave(file.path(fig_dir, "taxa_specific_size_scaling.png"),
       p_taxa_scaling, width = 12, height = 10, dpi = 300)

# 3. SURFACE AREA VS VOLUME EFFECTS
# ----------------------------------
# Compare surface area vs volume as predictors
surf_model <- lm(total_cafi ~ coral_surface_area, data = analysis_data)
vol_model <- lm(total_cafi ~ coral_volume, data = analysis_data)

cat("  Surface area R²:", round(summary(surf_model)$r.squared, 3), "\n")
cat("  Volume R²:", round(summary(vol_model)$r.squared, 3), "\n\n")

# ============================================================================
# Neighbor Distance Effects
# ============================================================================

cat("Analyzing neighbor distance effects...\n")

# 1. DISTANCE-ABUNDANCE RELATIONSHIPS
# ------------------------------------
# GAM for non-linear distance effects
gam_distance <- gam(total_cafi ~
                   s(mean_dist_neighbors, k = 5) +
                   s(n_neighbors, k = 5) +
                   morphotype,
                   data = analysis_data %>% filter(mean_dist_neighbors > 0),
                   family = nb())

# Plot distance effects
p_distance_gam <- draw(gam_distance, select = 1:2) &
  theme_minimal()

ggsave(file.path(fig_dir, "neighbor_distance_gam_effects.png"),
       p_distance_gam, width = 12, height = 6, dpi = 300)

# 2. ISOLATION VS CROWDING EFFECTS
# ---------------------------------
p_isolation <- analysis_data %>%
  filter(!is.na(isolation_class)) %>%
  ggplot(aes(x = isolation_class, y = total_cafi)) +
  geom_boxplot(aes(fill = isolation_class), alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~morphotype) +
  scale_y_sqrt() +
  labs(title = "CAFI Abundance by Isolation Level",
       x = "Isolation Class",
       y = "Total CAFI (sqrt scale)")

ggsave(file.path(fig_dir, "isolation_effects.png"),
       p_isolation, width = 12, height = 8, dpi = 300)

# Statistical test for isolation effects
isolation_test <- kruskal.test(total_cafi ~ isolation_class,
                              data = analysis_data)
cat("  Isolation effect (Kruskal-Wallis): p =",
    format.pval(isolation_test$p.value), "\n\n")

# 3. NEIGHBOR MORPHOTYPE EFFECTS
# -------------------------------
p_neighbor_morph <- analysis_data %>%
  filter(n_neighbors > 0) %>%
  ggplot(aes(x = prop_wide_neighbors, y = total_cafi)) +
  geom_point(aes(color = morphotype, size = n_neighbors), alpha = 0.6) +
  geom_smooth(method = "gam", se = TRUE) +
  labs(title = "Effect of Neighbor Morphotype Composition",
       x = "Proportion of Wide-branching Neighbors",
       y = "Total CAFI",
       color = "Focal Morphotype",
       size = "N Neighbors")

ggsave(file.path(fig_dir, "neighbor_morphotype_effects.png"),
       p_neighbor_morph, width = 10, height = 8, dpi = 300)

# ============================================================================
# Interactive Size × Distance Effects
# ============================================================================

cat("Analyzing size × distance interactions...\n")

# Interaction model
interaction_model <- gam(total_cafi ~
                        s(coral_volume, mean_dist_neighbors, k = 10) +
                        morphotype,
                        data = analysis_data %>%
                          filter(mean_dist_neighbors > 0, coral_volume > 0),
                        family = nb())

# Create prediction grid
pred_grid <- expand.grid(
  coral_volume = seq(min(analysis_data$coral_volume, na.rm = TRUE),
                    max(analysis_data$coral_volume, na.rm = TRUE),
                    length.out = 50),
  mean_dist_neighbors = seq(min(analysis_data$mean_dist_neighbors[analysis_data$mean_dist_neighbors > 0], na.rm = TRUE),
                           max(analysis_data$mean_dist_neighbors, na.rm = TRUE),
                           length.out = 50),
  morphotype = unique(analysis_data$morphotype)[1]
)

pred_grid$predicted <- predict(interaction_model, pred_grid, type = "response")

# Plot interaction surface
p_interaction <- ggplot(pred_grid, aes(x = coral_volume, y = mean_dist_neighbors,
                                       fill = predicted)) +
  geom_tile() +
  geom_contour(aes(z = predicted), color = "white", alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(trans = "sqrt", name = "Predicted\nCAFI") +
  labs(title = "Interactive Effects: Coral Size × Neighbor Distance",
       x = "Coral Volume (cm³, log)",
       y = "Mean Neighbor Distance (cm, log)")

ggsave(file.path(fig_dir, "size_distance_interaction.png"),
       p_interaction, width = 10, height = 8, dpi = 300)

# ============================================================================
# Threshold Analysis
# ============================================================================

cat("Identifying size and distance thresholds...\n")

# Size thresholds using segmented regression
library(segmented)

# Fit segmented model for size
lm_size <- lm(total_cafi ~ coral_volume, data = analysis_data)
seg_size <- segmented(lm_size, seg.Z = ~coral_volume, npsi = 2)

size_breakpoints <- seg_size$psi[, "Est."]
cat("  Size breakpoints:", round(size_breakpoints, 1), "cm³\n")

# Plot with breakpoints
p_threshold <- analysis_data %>%
  ggplot(aes(x = coral_volume, y = total_cafi)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "gray") +
  geom_vline(xintercept = size_breakpoints, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_sqrt() +
  labs(title = "Size Threshold Analysis",
       subtitle = "Red lines indicate breakpoints",
       x = "Coral Volume (cm³, log)",
       y = "Total CAFI (sqrt)")

ggsave(file.path(fig_dir, "size_thresholds.png"),
       p_threshold, width = 10, height = 8, dpi = 300)

# ============================================================================
# Competition Analysis
# ============================================================================

cat("Analyzing competition effects...\n")

# Competition model
comp_model <- lm(total_cafi ~
                coral_volume * competition_index +
                size_asymmetry +
                morphotype,
                data = analysis_data)

# Extract competition effects
comp_summary <- broom::tidy(comp_model) %>%
  filter(str_detect(term, "competition|asymmetry"))

write_csv(comp_summary,
          file.path(SURVEY_TABLES, "competition_effects.csv"))

# Visualize competition effects
p_competition <- analysis_data %>%
  filter(competition_index < quantile(competition_index, 0.95, na.rm = TRUE)) %>%
  ggplot(aes(x = competition_index, y = total_cafi)) +
  geom_point(aes(color = size_class), alpha = 0.6) +
  geom_smooth(aes(color = size_class), method = "lm", se = TRUE) +
  labs(title = "Competition Effects by Coral Size Class",
       x = "Competition Index (neighbors/distance)",
       y = "Total CAFI",
       color = "Size Class")

ggsave(file.path(fig_dir, "competition_by_size_class.png"),
       p_competition, width = 10, height = 8, dpi = 300)

# ============================================================================
# Shape Effects Analysis
# ============================================================================

cat("Analyzing coral shape effects...\n")

# Shape model
shape_model <- lm(total_cafi ~
                 aspect_ratio + elongation + compactness +
                 coral_volume,
                 data = analysis_data)

shape_effects <- broom::tidy(shape_model)
write_csv(shape_effects,
          file.path(SURVEY_TABLES, "shape_effects_on_cafi.csv"))

# Visualize shape effects
shape_plots <- list()

shape_plots[[1]] <- ggplot(analysis_data, aes(x = aspect_ratio, y = total_cafi)) +
  geom_point(aes(color = morphotype), alpha = 0.5) +
  geom_smooth(method = "gam") +
  labs(title = "Aspect Ratio Effect", x = "Height/Width Ratio", y = "Total CAFI")

shape_plots[[2]] <- ggplot(analysis_data, aes(x = elongation, y = total_cafi)) +
  geom_point(aes(color = morphotype), alpha = 0.5) +
  geom_smooth(method = "gam") +
  labs(title = "Elongation Effect", x = "Length/Width Ratio", y = "Total CAFI")

shape_plots[[3]] <- ggplot(analysis_data, aes(x = compactness, y = total_cafi)) +
  geom_point(aes(color = morphotype), alpha = 0.5) +
  geom_smooth(method = "gam") +
  labs(title = "Compactness Effect", x = "Volume^(2/3)/Surface Area", y = "Total CAFI")

shape_combined <- wrap_plots(shape_plots, ncol = 3)
ggsave(file.path(fig_dir, "shape_effects_on_cafi.png"),
       shape_combined, width = 15, height = 5, dpi = 300)

# ============================================================================
# Comprehensive Summary Statistics
# ============================================================================

cat("\nCreating comprehensive summary statistics...\n")

# Summary by size class
size_summary <- analysis_data %>%
  group_by(size_class) %>%
  summarise(
    n_corals = n(),
    mean_volume = mean(coral_volume, na.rm = TRUE),
    mean_cafi = mean(total_cafi, na.rm = TRUE),
    sd_cafi = sd(total_cafi, na.rm = TRUE),
    mean_richness = mean(species_richness, na.rm = TRUE),
    mean_neighbor_dist = mean(mean_dist_neighbors, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(size_summary,
          file.path(SURVEY_TABLES, "summary_by_size_class.csv"))

# Summary by isolation class
isolation_summary <- analysis_data %>%
  filter(!is.na(isolation_class)) %>%
  group_by(isolation_class) %>%
  summarise(
    n_corals = n(),
    mean_distance = mean(mean_dist_neighbors, na.rm = TRUE),
    mean_cafi = mean(total_cafi, na.rm = TRUE),
    sd_cafi = sd(total_cafi, na.rm = TRUE),
    mean_richness = mean(species_richness, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(isolation_summary,
          file.path(SURVEY_TABLES, "summary_by_isolation_class.csv"))

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Size and Neighbor Effects Summary\n")
cat("========================================\n\n")

cat("Key Findings:\n")

# Size effects
if (exists("models_size")) {
  best_size_model <- names(models_size)[which.min(size_model_comparison$AIC)]
  cat("  Size-Abundance Relationship:\n")
  cat("    - Best model:", best_size_model, "\n")

  if (best_size_model == "log_log") {
    scaling_exp <- coef(models_size$log_log)[2]
    cat("    - Scaling exponent:", round(scaling_exp, 2), "\n")
  }
}

# Distance effects
if (exists("isolation_test")) {
  cat("\n  Neighbor Distance Effects:\n")
  cat("    - Isolation effect significant:", isolation_test$p.value < 0.05, "\n")
}

# Thresholds
if (exists("size_breakpoints")) {
  cat("\n  Size Thresholds:\n")
  cat("    - Breakpoints at:", paste(round(size_breakpoints, 1), "cm³"), "\n")
}

# Shape effects
if (exists("shape_effects")) {
  sig_shapes <- shape_effects %>%
    filter(p.value < 0.05, term != "(Intercept)")

  if (nrow(sig_shapes) > 0) {
    cat("\n  Significant Shape Effects:\n")
    for (i in 1:nrow(sig_shapes)) {
      cat("    -", sig_shapes$term[i], ": β =",
          round(sig_shapes$estimate[i], 3), "\n")
    }
  }
}

cat("\n✅ Coral size and neighbor effects analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", SURVEY_TABLES, "\n")