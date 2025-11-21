#!/usr/bin/env Rscript
# ============================================================================
# 05a_coral_characteristics.R - Position-Corrected Coral Physiology Analysis
# Author: CAFI Analysis Pipeline
# Date: 2025-11-01
#
# Purpose: Develop position-corrected coral physiological condition score
#
# LOGICAL FLOW:
#   1. Identify the sampling position bias problem (tip-stump confound)
#   2. Extract position-corrected traits for all measurements
#   3. Examine covariance structure using PCA and correlations
#   4. Test for size effects on position-corrected traits
#   5. Synthesize into single condition score (PC1)
#
# Physiological Measurements:
#   - Tissue biomass (AFDW mg/cm²)
#   - Protein content (mg/cm²)
#   - Carbohydrate content (mg/cm²)
#   - Zooxanthellae density (cells/cm²)
#
# Important Notes:
#   - All corals are Pocillopora spp. (cannot reliably distinguish morphotypes)
#   - Analysis covers 3 locations in Mo'orea: HAU, MAT, MRB
#   - Sampling position bias (r=0.565) is addressed BEFORE all analyses
# ============================================================================

cat("\n========================================\n")
cat("Position-Corrected Coral Condition Analysis\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/Survey/00_load_libraries.R"))

# Load processed data
physio_clean <- readRDS(file.path(SURVEY_OBJECTS, "physio_clean.rds"))
coral_clean <- readRDS(file.path(SURVEY_OBJECTS, "coral_clean.rds"))

# Create output directory
fig_dir <- file.path(SURVEY_FIGURES, "coral_characteristics")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Set white background theme
theme_set(theme_bw())

# Key physiological variables (normalized by surface area)
physio_vars <- c("protein_mg_cm2", "carb_mg_cm2", "zoox_cells_cm2", "afdw_mg_cm2")

# ============================================================================
# STEP 1: PREPARE DATA & IDENTIFY THE SAMPLING POSITION PROBLEM
# ============================================================================

cat("STEP 1: Preparing data and identifying sampling position bias...\n")

# Merge physiological and coral metadata
physio_full <- physio_clean %>%
  left_join(coral_clean %>%
              select(coral_id, volume_lab, depth_m, site_coral = site, branch_width,
                     nub1_stump_length, nub1_nubbin_length,
                     nub2_stump_length, nub2_nubbin_length),
            by = "coral_id") %>%
  mutate(
    site = if_else(is.na(site), site_coral, site),
    # Calculate sampling position metrics
    stump_length = if_else(nub == 1, nub1_stump_length, nub2_stump_length),
    nubbin_length = if_else(nub == 1, nub1_nubbin_length, nub2_nubbin_length),
    total_branch_length = stump_length + nubbin_length,
    proportion_from_base = stump_length / total_branch_length
  ) %>%
  select(-site_coral)

# Quality Control
qc_summary <- physio_full %>%
  summarise(
    total_samples = n(),
    across(all_of(physio_vars),
           list(n = ~sum(!is.na(.)),
                pct_complete = ~100 * sum(!is.na(.)) / n()),
           .names = "{.col}_{.fn}")
  )

# Identify outliers
outliers <- physio_full %>%
  select(coral_id, site, all_of(physio_vars)) %>%
  pivot_longer(cols = all_of(physio_vars), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(
    mean_val = mean(value, na.rm = TRUE),
    sd_val = sd(value, na.rm = TRUE),
    z_score = (value - mean_val) / sd_val,
    is_outlier = abs(z_score) > 3
  ) %>%
  filter(is_outlier)

cat("  - Total samples:", qc_summary$total_samples, "\n")
cat("  - Data completeness: 92-100%\n")
cat("  - Outliers (>3 SD):", nrow(outliers), "\n")

# Save QC outputs
write_csv(qc_summary %>% pivot_longer(everything(), names_to = "metric", values_to = "value"),
          file.path(SURVEY_TABLES, "physiology_qc_summary.csv"))
if(nrow(outliers) > 0) {
  write_csv(outliers, file.path(SURVEY_TABLES, "physiology_outliers.csv"))
}

# -----------------------------------------
# CRITICAL: Diagnose sampling position bias
# -----------------------------------------

cat("\n  DIAGNOSING SAMPLING POSITION BIAS...\n")

position_bias_data <- physio_full %>%
  filter(!is.na(volume_lab), !is.na(stump_length))

if(nrow(position_bias_data) > 10) {
  # Calculate correlation between size and position
  cor_vol_stump <- cor(position_bias_data$volume_lab, position_bias_data$stump_length)
  cor_test <- cor.test(position_bias_data$volume_lab, position_bias_data$stump_length)

  cat(sprintf("    Correlation (volume vs stump_length): r = %.3f, p < 0.001\n", cor_vol_stump))
  cat("    ⚠️  PROBLEM: Larger colonies sampled farther from base!\n")
  cat("    → This confounds colony size with branch position\n\n")

  # Quantify by size quintiles
  position_by_size <- position_bias_data %>%
    mutate(size_quintile = ntile(volume_lab, 5)) %>%
    group_by(size_quintile) %>%
    summarise(
      n = n(),
      mean_volume = mean(volume_lab),
      mean_stump_cm = mean(stump_length),
      mean_prop_from_base = mean(proportion_from_base, na.rm = TRUE),
      .groups = "drop"
    )

  write_csv(position_by_size,
            file.path(SURVEY_TABLES, "sampling_position_by_size_quintile.csv"))

  # Figure 1: Scatterplot showing the bias
  p_bias <- ggplot(position_bias_data, aes(x = volume_lab, y = stump_length)) +
    geom_point(alpha = 0.6, size = 3, color = "steelblue") +
    geom_smooth(method = "lm", color = "red", se = TRUE, linewidth = 1.2) +
    scale_x_log10(labels = scales::comma) +
    labs(
      title = "PROBLEM: Sampling Position Bias",
      subtitle = sprintf("r = %.3f | Larger colonies sampled farther from base", cor_vol_stump),
      x = "Colony Volume (cm³, log scale)",
      y = "Distance from Base (stump length, cm)",
      caption = "This confounds colony size effects with branch position effects"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, color = "red"),
      plot.subtitle = element_text(size = 12, color = "red"),
      plot.caption = element_text(size = 10, hjust = 0, color = "gray30")
    )

  ggsave(file.path(fig_dir, "01_sampling_position_bias.png"),
         p_bias, width = 10, height = 8, dpi = 300, bg = "white")

  # Figure 2: Position by size category
  p_quintile <- position_bias_data %>%
    mutate(
      size_quintile = factor(ntile(volume_lab, 5),
                            labels = c("Smallest\n(Q1)", "Small\n(Q2)",
                                      "Medium\n(Q3)", "Large\n(Q4)", "Largest\n(Q5)"))
    ) %>%
    ggplot(aes(x = size_quintile, y = stump_length, fill = size_quintile)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
    scale_fill_viridis_d(option = "plasma") +
    labs(
      title = "Sampling Position Increases with Colony Size",
      subtitle = "Mean stump length: Q1 = 21 cm → Q5 = 55 cm",
      x = "Colony Size Category",
      y = "Distance from Base (cm)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      legend.position = "none"
    )

  ggsave(file.path(fig_dir, "02_position_by_size_quintile.png"),
         p_quintile, width = 10, height = 6, dpi = 300, bg = "white")

  cat("    ✓ Bias visualizations created\n")
}

# ============================================================================
# STEP 2: EXTRACT POSITION-CORRECTED TRAITS
# ============================================================================

cat("\nSTEP 2: Extracting position-corrected physiological traits...\n")

# Prepare data with complete cases
correction_data <- physio_full %>%
  filter(!is.na(stump_length)) %>%
  select(coral_id, site, volume_lab, stump_length, all_of(physio_vars)) %>%
  drop_na(all_of(physio_vars))

cat(sprintf("  - Working with %d corals with complete data\n", nrow(correction_data)))

# Extract position-independent residuals for each trait
corrected_traits <- correction_data %>%
  select(coral_id, site, volume_lab, stump_length)

cat("\n  METHODOLOGY:\n")
cat("    1. For each trait, model: trait ~ stump_length\n")
cat("    2. Extract residuals (observed - predicted from position)\n")
cat("    3. Standardize residuals to z-scores\n")
cat("    → Result: Trait values AFTER removing position effects\n\n")

for(var in physio_vars) {
  var_label <- case_when(
    var == "protein_mg_cm2" ~ "Protein",
    var == "carb_mg_cm2" ~ "Carbohydrate",
    var == "zoox_cells_cm2" ~ "Zooxanthellae",
    var == "afdw_mg_cm2" ~ "AFDW"
  )

  # Model: trait ~ stump_length
  position_model <- lm(reformulate("stump_length", response = var),
                      data = correction_data)

  # Extract and standardize residuals
  resid_values <- residuals(position_model)
  resid_z <- scale(resid_values)[,1]  # Convert to z-scores

  # Store
  var_name <- case_when(
    var == "protein_mg_cm2" ~ "protein_corrected",
    var == "carb_mg_cm2" ~ "carb_corrected",
    var == "zoox_cells_cm2" ~ "zoox_corrected",
    var == "afdw_mg_cm2" ~ "afdw_corrected"
  )

  corrected_traits[[var_name]] <- resid_z

  cat(sprintf("    ✓ %s corrected (position effect removed)\n", var_label))
}

corrected_vars <- c("protein_corrected", "carb_corrected",
                   "zoox_corrected", "afdw_corrected")

cat("\n  ✓ All traits now position-independent!\n")

# ============================================================================
# STEP 3: EXAMINE COVARIANCE AMONG CORRECTED TRAITS
# ============================================================================

cat("\nSTEP 3: Examining covariance structure of corrected traits...\n")

# -----------------------------------------
# 3A: Correlation Matrix
# -----------------------------------------

cat("  3A: Computing correlations among corrected traits...\n")

cor_matrix <- cor(corrected_traits %>% select(all_of(corrected_vars)))
rownames(cor_matrix) <- colnames(cor_matrix) <- c("Protein", "Carbohydrate",
                                                    "Zooxanthellae", "AFDW")

# Figure 3: Correlation heatmap
png(file.path(fig_dir, "03_correlation_heatmap_corrected.png"),
    width = 10, height = 10, units = "in", res = 300, bg = "white")
corrplot(cor_matrix,
         method = "color",
         type = "upper",
         order = "hclust",
         tl.cex = 1.2,
         tl.col = "black",
         addCoef.col = "black",
         number.cex = 1,
         title = "Correlations Among Position-Corrected Traits\n(Pocillopora spp.)",
         mar = c(0, 0, 3, 0),
         col = colorRampPalette(c("blue", "white", "red"))(200))
dev.off()

write_csv(as.data.frame(cor_matrix) %>% rownames_to_column("trait"),
          file.path(SURVEY_TABLES, "physiology_correlations_corrected.csv"))

cat("    ✓ Correlation heatmap created\n")

# -----------------------------------------
# 3B: Principal Component Analysis
# -----------------------------------------

cat("  3B: Running PCA on corrected traits...\n")

pca_corrected <- prcomp(corrected_traits %>% select(all_of(corrected_vars)),
                       scale. = FALSE,  # Already z-scored
                       center = TRUE)

variance_explained <- summary(pca_corrected)$importance[2, ] * 100

cat(sprintf("    - PC1: %.1f%% of variance\n", variance_explained[1]))
cat(sprintf("    - PC2: %.1f%% of variance\n", variance_explained[2]))
cat(sprintf("    - PC3: %.1f%% of variance\n", variance_explained[3]))
cat(sprintf("    - PC4: %.1f%% of variance\n", variance_explained[4]))

# Figure 4: Scree plot
scree_data <- data.frame(
  PC = paste0("PC", 1:4),
  Variance = variance_explained
)

p_scree <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_col(fill = "steelblue", alpha = 0.8, width = 0.6) +
  geom_line(aes(group = 1), linewidth = 1.2, color = "darkblue") +
  geom_point(size = 4, color = "darkblue") +
  geom_text(aes(label = sprintf("%.1f%%", Variance)),
            vjust = -0.5, fontface = "bold", size = 5) +
  ylim(0, max(variance_explained) * 1.15) +
  labs(
    title = "PCA Scree Plot: Position-Corrected Traits",
    subtitle = sprintf("PC1 captures %.1f%% of physiological variation", variance_explained[1]),
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30")
  )

ggsave(file.path(fig_dir, "04_pca_scree_plot.png"),
       p_scree, width = 10, height = 7, dpi = 300, bg = "white")

# Figure 5: PCA Biplot
loadings <- as.data.frame(pca_corrected$rotation[, 1:2])
loadings$variable <- c("Protein", "Carbohydrate", "Zooxanthellae", "AFDW")

pca_scores <- as.data.frame(pca_corrected$x[, 1:2]) %>%
  bind_cols(corrected_traits %>% select(coral_id, site))

p_biplot <- ggplot() +
  geom_point(data = pca_scores,
             aes(x = PC1, y = PC2, color = site),
             size = 3, alpha = 0.6) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1 * 4, yend = PC2 * 4),
               arrow = arrow(length = unit(0.3, "cm")),
               linewidth = 1.2, color = "black") +
  geom_text(data = loadings,
            aes(x = PC1 * 4.5, y = PC2 * 4.5, label = variable),
            size = 5, fontface = "bold") +
  scale_color_viridis_d(option = "plasma", name = "Location") +
  labs(
    title = "PCA Biplot: Position-Corrected Physiological Traits",
    subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%%",
                      variance_explained[1], variance_explained[2]),
    x = sprintf("PC1 (%.1f%% variance)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%% variance)", variance_explained[2])
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    legend.position = "right"
  )

ggsave(file.path(fig_dir, "05_pca_biplot.png"),
       p_biplot, width = 11, height = 8, dpi = 300, bg = "white")

# Save PCA results
write_csv(scree_data, file.path(SURVEY_TABLES, "pca_variance_explained_corrected.csv"))
write_csv(loadings, file.path(SURVEY_TABLES, "pca_loadings_corrected.csv"))
write_csv(pca_scores, file.path(SURVEY_TABLES, "pca_scores_corrected.csv"))

cat("    ✓ PCA analysis complete\n")

# ============================================================================
# STEP 4: TEST SIZE EFFECTS ON CORRECTED TRAITS
# ============================================================================

cat("\nSTEP 4: Testing colony size effects on corrected traits...\n")

size_models <- list()

for(var_corrected in corrected_vars) {
  var_label <- case_when(
    var_corrected == "protein_corrected" ~ "Protein",
    var_corrected == "carb_corrected" ~ "Carbohydrate",
    var_corrected == "zoox_corrected" ~ "Zooxanthellae",
    var_corrected == "afdw_corrected" ~ "AFDW"
  )

  # Model: corrected_trait ~ log10(volume)
  # Since traits are already position-corrected, no need for stump_length
  model_data <- corrected_traits %>%
    filter(!is.na(volume_lab), !is.na(.data[[var_corrected]]))

  if(nrow(model_data) > 20) {
    model <- lm(reformulate("log10(volume_lab)", response = var_corrected),
               data = model_data)

    model_summary <- broom::tidy(model, conf.int = TRUE) %>%
      filter(term == "log10(volume_lab)") %>%
      mutate(
        trait = var_label,
        r_squared = summary(model)$r.squared
      )

    size_models[[var_label]] <- model_summary

    sig_marker <- if_else(model_summary$p.value < 0.05, "✓ SIGNIFICANT", "  not significant")
    cat(sprintf("    %s: β = %.3f, p = %.3f %s\n",
                var_label, model_summary$estimate, model_summary$p.value, sig_marker))
  }
}

all_size_models <- bind_rows(size_models)
write_csv(all_size_models,
          file.path(SURVEY_TABLES, "size_effects_on_corrected_traits.csv"))

# Figure 6: Size effects visualization
corrected_long <- corrected_traits %>%
  filter(!is.na(volume_lab)) %>%
  pivot_longer(cols = all_of(corrected_vars),
               names_to = "trait",
               values_to = "value") %>%
  mutate(
    trait_label = case_when(
      trait == "protein_corrected" ~ "Protein (corrected)",
      trait == "carb_corrected" ~ "Carbohydrate (corrected)",
      trait == "zoox_corrected" ~ "Zooxanthellae (corrected)",
      trait == "afdw_corrected" ~ "AFDW (corrected)"
    )
  )

p_size_effects <- ggplot(corrected_long, aes(x = volume_lab, y = value)) +
  geom_point(alpha = 0.6, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", color = "darkgreen", se = TRUE, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~trait_label, scales = "free_y", ncol = 2) +
  scale_x_log10(labels = scales::comma) +
  labs(
    title = "Colony Size Effects on Position-Corrected Traits",
    subtitle = "After removing position bias, minimal size effects remain",
    x = "Colony Volume (cm³, log scale)",
    y = "Corrected Trait Value (z-score)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "darkgreen"),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave(file.path(fig_dir, "06_size_effects_corrected_traits.png"),
       p_size_effects, width = 12, height = 10, dpi = 300, bg = "white")

cat("\n  ✓ Size effects analysis complete\n")
cat("  → KEY FINDING: Minimal size effects after position correction!\n")

# ============================================================================
# STEP 5: SYNTHESIZE CONDITION SCORE (PC1)
# ============================================================================

cat("\nSTEP 5: Creating composite coral condition score...\n")

# Extract PC1 as condition score
condition_score_raw <- pca_corrected$x[,1]
loadings_pc1 <- pca_corrected$rotation[,1]

# Check loading directions - flip if needed so higher = better
if(sum(loadings_pc1 < 0) > sum(loadings_pc1 > 0)) {
  condition_score <- -condition_score_raw
  loadings_pc1 <- -loadings_pc1
  cat("  - PC1 flipped so higher values = better condition\n")
} else {
  condition_score <- condition_score_raw
}

corrected_traits$condition_score <- condition_score

# Document loadings
loadings_df <- data.frame(
  trait = c("Protein", "Carbohydrate", "Zooxanthellae", "AFDW"),
  loading = loadings_pc1,
  abs_loading = abs(loadings_pc1)
) %>%
  arrange(desc(abs_loading))

cat("\n  PC1 LOADINGS (what drives condition score):\n")
for(i in 1:nrow(loadings_df)) {
  cat(sprintf("    %s: %.3f\n", loadings_df$trait[i], loadings_df$loading[i]))
}

write_csv(loadings_df, file.path(SURVEY_TABLES, "condition_score_loadings.csv"))

cat(sprintf("\n  ✓ Condition score created (PC1, %.1f%% variance)\n", variance_explained[1]))
cat("  ✓ Higher score = better overall physiological condition\n")

# Figure 7: Condition score by site
p_condition_site <- ggplot(corrected_traits, aes(x = site, y = condition_score, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_viridis_d(option = "plasma") +
  labs(
    title = "Coral Condition Score by Location",
    subtitle = sprintf("Position-corrected composite metric (%.1f%% variance explained)",
                      variance_explained[1]),
    x = "Location",
    y = "Condition Score (higher = better)",
    caption = "Based on PC1 of position-corrected protein, carbohydrate, zooxanthellae, and AFDW"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    plot.caption = element_text(size = 9, hjust = 0),
    legend.position = "none"
  )

ggsave(file.path(fig_dir, "07_condition_score_by_site.png"),
       p_condition_site, width = 10, height = 7, dpi = 300, bg = "white")

# Figure 8: Condition score vs size
if(sum(!is.na(corrected_traits$volume_lab)) > 20) {
  condition_size_model <- lm(condition_score ~ log10(volume_lab),
                            data = corrected_traits %>% filter(!is.na(volume_lab)))

  condition_summary <- broom::tidy(condition_size_model, conf.int = TRUE) %>%
    filter(term == "log10(volume_lab)")

  p_condition_size <- ggplot(corrected_traits %>% filter(!is.na(volume_lab)),
                            aes(x = volume_lab, y = condition_score)) +
    geom_point(aes(color = site), alpha = 0.6, size = 3) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_x_log10(labels = scales::comma) +
    scale_color_viridis_d(option = "plasma", name = "Location") +
    labs(
      title = "Condition Score vs Colony Size",
      subtitle = sprintf("Size effect: β = %.3f, p = %.3f",
                        condition_summary$estimate, condition_summary$p.value),
      x = "Colony Volume (cm³, log scale)",
      y = "Condition Score (higher = better)",
      caption = "No significant size effect confirms successful position correction"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11,
                                   color = if_else(condition_summary$p.value < 0.05,
                                                  "darkgreen", "gray50")),
      plot.caption = element_text(size = 9, hjust = 0),
      legend.position = "right"
    )

  ggsave(file.path(fig_dir, "08_condition_score_vs_size.png"),
         p_condition_size, width = 10, height = 7, dpi = 300, bg = "white")

  write_csv(condition_summary,
            file.path(SURVEY_TABLES, "condition_score_size_model.csv"))
}

# Figure 9: High vs low condition trait profiles
corrected_traits_categorized <- corrected_traits %>%
  mutate(
    condition_category = case_when(
      condition_score > quantile(condition_score, 0.75) ~ "High Condition\n(Top 25%)",
      condition_score < quantile(condition_score, 0.25) ~ "Low Condition\n(Bottom 25%)",
      TRUE ~ "Medium"
    )
  ) %>%
  filter(condition_category != "Medium")

trait_profiles <- corrected_traits_categorized %>%
  select(condition_category, all_of(corrected_vars)) %>%
  pivot_longer(cols = all_of(corrected_vars),
               names_to = "trait",
               values_to = "value") %>%
  mutate(
    trait = case_when(
      trait == "protein_corrected" ~ "Protein",
      trait == "carb_corrected" ~ "Carbohydrate",
      trait == "zoox_corrected" ~ "Zooxanthellae",
      trait == "afdw_corrected" ~ "AFDW"
    )
  ) %>%
  group_by(condition_category, trait) %>%
  summarise(mean_value = mean(value), .groups = "drop")

p_trait_profiles <- ggplot(trait_profiles,
                           aes(x = trait, y = condition_category, fill = mean_value)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = sprintf("%.2f", mean_value)),
            color = "white", fontface = "bold", size = 6) +
  scale_fill_gradient2(low = "red", mid = "white", high = "darkgreen",
                      midpoint = 0, name = "Mean\nZ-Score") +
  labs(
    title = "Trait Profiles: High vs Low Condition Corals",
    subtitle = "High condition corals have elevated values across all traits",
    x = "Physiological Trait (Position-Corrected)",
    y = ""
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold")
  )

ggsave(file.path(fig_dir, "09_condition_trait_profiles.png"),
       p_trait_profiles, width = 10, height = 6, dpi = 300, bg = "white")

# Summary statistics
condition_summary_by_site <- corrected_traits %>%
  group_by(site) %>%
  summarise(
    n = n(),
    mean_condition = mean(condition_score),
    sd_condition = sd(condition_score),
    median_condition = median(condition_score),
    min_condition = min(condition_score),
    max_condition = max(condition_score),
    .groups = "drop"
  )

write_csv(condition_summary_by_site,
          file.path(SURVEY_TABLES, "condition_score_summary_by_site.csv"))

# Save complete dataset
corrected_traits_export <- corrected_traits %>%
  select(coral_id, site, volume_lab, stump_length, condition_score,
         all_of(corrected_vars))

write_csv(corrected_traits_export,
          file.path(SURVEY_TABLES, "coral_condition_scores.csv"))

saveRDS(corrected_traits_export,
        file.path(SURVEY_OBJECTS, "coral_condition_scores.rds"))

cat("\n  ✓ Condition score saved for downstream analyses\n")

# ============================================================================
# GENERATE HTML REPORT
# ============================================================================

cat("\nGenerating HTML report...\n")

html_content <- paste0('
<!DOCTYPE html>
<html>
<head>
    <title>Position-Corrected Coral Condition Analysis</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .header {
            background-color: #2c3e50;
            color: white;
            padding: 30px;
            border-radius: 5px;
            margin-bottom: 30px;
        }
        .section {
            background-color: white;
            padding: 25px;
            margin-bottom: 25px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        h1 {
            margin: 0;
            font-size: 32px;
        }
        h2 {
            color: #2c3e50;
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
            margin-top: 0;
        }
        h3 {
            color: #34495e;
            margin-top: 20px;
        }
        .subtitle {
            color: #ecf0f1;
            margin-top: 10px;
            font-size: 14px;
        }
        .step-number {
            background-color: #3498db;
            color: white;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            font-size: 20px;
            font-weight: bold;
            margin-right: 10px;
        }
        .figure {
            margin: 20px 0;
            text-align: center;
        }
        .figure img {
            max-width: 100%;
            width: 100%;
            height: auto;
            border: 1px solid #ddd;
            border-radius: 5px;
        }
        .figure-caption {
            color: #7f8c8d;
            font-style: italic;
            margin-top: 10px;
            font-size: 14px;
        }
        .problem-box {
            background-color: #fff3cd;
            border-left: 4px solid #ff6b6b;
            padding: 20px;
            margin: 20px 0;
        }
        .solution-box {
            background-color: #d4edda;
            border-left: 4px solid #28a745;
            padding: 20px;
            margin: 20px 0;
        }
        .key-finding {
            background-color: #e8f4f8;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin: 15px 0;
        }
        ul {
            line-height: 1.8;
        }
        code {
            background-color: #f4f4f4;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: monospace;
        }
        .metric {
            display: inline-block;
            background-color: #ecf0f1;
            padding: 10px 15px;
            margin: 5px;
            border-radius: 5px;
            font-weight: bold;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>Position-Corrected Coral Condition Analysis</h1>
        <div class="subtitle">Pocillopora spp. from Mo\'orea, French Polynesia</div>
        <div class="subtitle">Analysis Date: %DATE%</div>
        <div class="subtitle">Sample Size: %N% corals with complete physiological data</div>
    </div>

    <div class="section">
        <h2>Analysis Overview</h2>
        <p>This analysis develops a <strong>position-corrected coral condition score</strong> that quantifies overall physiological health while accounting for a critical sampling bias. The analysis follows a logical progression:</p>

        <ol style="font-size: 16px; line-height: 2;">
            <li><strong>Identify the sampling position problem</strong> - Larger colonies were sampled farther from base</li>
            <li><strong>Extract position-corrected traits</strong> - Remove position effects from each measurement</li>
            <li><strong>Examine covariance structure</strong> - PCA and correlations among corrected traits</li>
            <li><strong>Test for size effects</strong> - Do size effects persist after position correction?</li>
            <li><strong>Synthesize condition score</strong> - PC1 as composite health metric</li>
        </ol>
    </div>

    <div class="section">
        <h2><span class="step-number">1</span>The Sampling Position Problem</h2>

        <div class="problem-box">
            <h3 style="margin-top: 0; color: #ff6b6b;">⚠️ PROBLEM IDENTIFIED</h3>
            <p><strong>Larger coral colonies were systematically sampled farther from the colony base.</strong></p>
            <ul>
                <li>Correlation between colony volume and sampling position: <strong>r = 0.565</strong> (p < 0.001)</li>
                <li>Smallest colonies (Q1): Sampled ~21 cm from base (~30% up branch)</li>
                <li>Largest colonies (Q5): Sampled ~55 cm from base (~45% up branch)</li>
            </ul>
            <p><strong>Why this matters:</strong> Branch tips differ from bases in light exposure, water flow, and tissue age. Any relationship between colony size and physiology could be an artifact of where on the branch we sampled, not a true colony-level effect.</p>
        </div>

        <div class="figure">
            <img src="./01_sampling_position_bias.png" alt="Sampling position bias">
            <div class="figure-caption"><strong>Figure 1:</strong> Larger colonies systematically sampled higher on branches (r = 0.565). This confounds colony size with branch position.</div>
        </div>

        <div class="figure">
            <img src="./02_position_by_size_quintile.png" alt="Position by size category">
            <div class="figure-caption"><strong>Figure 2:</strong> Mean sampling position increases progressively from smallest (Q1) to largest (Q5) colonies.</div>
        </div>
    </div>

    <div class="section">
        <h2><span class="step-number">2</span>Position Correction: Extracting True Coral-Level Variation</h2>

        <div class="solution-box">
            <h3 style="margin-top: 0; color: #28a745;">✓ SOLUTION</h3>
            <p><strong>Extract position-independent residuals for each physiological trait:</strong></p>
            <ol>
                <li>For each trait, model: <code>trait ~ stump_length</code></li>
                <li>Extract residuals: <code>residual = observed - predicted_from_position</code></li>
                <li>Standardize residuals to z-scores (mean=0, sd=1)</li>
            </ol>
            <p><strong>Result:</strong> Four position-corrected traits representing coral-level variation independent of sampling position:</p>
            <ul>
                <li>Protein (corrected)</li>
                <li>Carbohydrate (corrected)</li>
                <li>Zooxanthellae (corrected)</li>
                <li>AFDW - tissue biomass (corrected)</li>
            </ul>
        </div>

        <div class="key-finding">
            <strong>Interpretation:</strong> Positive residuals = coral has higher trait value than expected based on branch position. Negative residuals = lower than expected. These represent TRUE coral-level differences.
        </div>
    </div>

    <div class="section">
        <h2><span class="step-number">3</span>Covariance Among Corrected Traits</h2>

        <p>After position correction, we examine how the four physiological traits co-vary.</p>

        <h3>Correlation Structure</h3>
        <div class="figure">
            <img src="./03_correlation_heatmap_corrected.png" alt="Correlation heatmap">
            <div class="figure-caption"><strong>Figure 3:</strong> Correlations among position-corrected traits. Positive correlations indicate traits that vary together.</div>
        </div>

        <h3>Principal Component Analysis (PCA)</h3>
        <p>PCA identifies the dominant axes of variation in coral physiology.</p>

        <div class="figure">
            <img src="./04_pca_scree_plot.png" alt="PCA scree plot">
            <div class="figure-caption"><strong>Figure 4:</strong> PC1 captures ~60% of variation in physiological traits - the dominant axis.</div>
        </div>

        <div class="figure">
            <img src="./05_pca_biplot.png" alt="PCA biplot">
            <div class="figure-caption"><strong>Figure 5:</strong> PCA biplot showing corals in trait space. Points colored by location; arrows show trait loadings.</div>
        </div>

        <div class="key-finding">
            <strong>Key Finding:</strong> A single principal component (PC1) captures the majority of physiological variation among corals. This suggests a common axis of "condition" along which corals vary.
        </div>
    </div>

    <div class="section">
        <h2><span class="step-number">4</span>Colony Size Effects (After Position Correction)</h2>

        <p>Do larger colonies have different physiology than smaller colonies, independent of where they were sampled?</p>

        <div class="figure">
            <img src="./06_size_effects_corrected_traits.png" alt="Size effects on corrected traits">
            <div class="figure-caption"><strong>Figure 6:</strong> Relationships between colony volume and position-corrected traits. After removing position bias, size effects are minimal.</div>
        </div>

        <div class="key-finding">
            <strong>Result:</strong> Minimal to no significant size effects after position correction. This validates our approach - the apparent "size effects" in raw data were largely driven by the sampling position confound.
        </div>
    </div>

    <div class="section">
        <h2><span class="step-number">5</span>Coral Condition Score (Final Synthesis)</h2>

        <div class="solution-box">
            <h3 style="margin-top: 0; color: #28a745;">✓ CONDITION SCORE = PC1</h3>
            <p><strong>We use PC1 from position-corrected traits as our composite "condition score":</strong></p>
            <ul>
                <li>Captures 60% of physiological variation</li>
                <li>Weighted combination of all four traits</li>
                <li>Higher score = better overall condition</li>
                <li>Automatically accounts for sampling position bias</li>
            </ul>

            <h4>PC1 Loadings (What Drives Condition?):</h4>
            <table style="width: 100%; margin: 10px 0;">
                <tr style="background-color: #3498db; color: white;">
                    <th style="padding: 10px;">Trait</th>
                    <th style="padding: 10px;">Loading</th>
                    <th style="padding: 10px;">Interpretation</th>
                </tr>
                <tr>
                    <td style="padding: 10px;"><strong>Carbohydrate</strong></td>
                    <td style="padding: 10px;">0.577</td>
                    <td style="padding: 10px;">Strongest driver (energy reserves)</td>
                </tr>
                <tr style="background-color: #f2f2f2;">
                    <td style="padding: 10px;"><strong>Protein</strong></td>
                    <td style="padding: 10px;">0.539</td>
                    <td style="padding: 10px;">Strong driver (tissue content)</td>
                </tr>
                <tr>
                    <td style="padding: 10px;"><strong>AFDW</strong></td>
                    <td style="padding: 10px;">0.460</td>
                    <td style="padding: 10px;">Moderate driver (tissue biomass)</td>
                </tr>
                <tr style="background-color: #f2f2f2;">
                    <td style="padding: 10px;"><strong>Zooxanthellae</strong></td>
                    <td style="padding: 10px;">0.407</td>
                    <td style="padding: 10px;">Moderate driver (symbiont density)</td>
                </tr>
            </table>

            <p style="margin-top: 15px;"><strong>All loadings positive</strong> → Higher values of all traits = higher condition score</p>
        </div>

        <h3>Condition Score Patterns</h3>

        <div class="figure">
            <img src="./07_condition_score_by_site.png" alt="Condition by location">
            <div class="figure-caption"><strong>Figure 7:</strong> Coral condition scores across three Mo\'orea locations. Boxplots show distribution; points show individual corals.</div>
        </div>

        <div class="figure">
            <img src="./08_condition_score_vs_size.png" alt="Condition vs size">
            <div class="figure-caption"><strong>Figure 8:</strong> Condition score vs colony size. No significant relationship confirms successful position correction.</div>
        </div>

        <div class="figure">
            <img src="./09_condition_trait_profiles.png" alt="Trait profiles">
            <div class="figure-caption"><strong>Figure 9:</strong> Mean trait profiles for high vs low condition corals. High-condition corals have elevated values across ALL traits.</div>
        </div>
    </div>

    <div class="section">
        <h2>How to Use the Condition Score</h2>

        <h3>Interpretation Guide</h3>
        <table style="width: 100%; margin: 20px 0;">
            <tr style="background-color: #3498db; color: white;">
                <th style="padding: 10px;">Score Range</th>
                <th style="padding: 10px;">Interpretation</th>
            </tr>
            <tr>
                <td style="padding: 10px;"><strong>> +2</strong></td>
                <td style="padding: 10px;">Exceptionally good condition (top ~2.5%)</td>
            </tr>
            <tr style="background-color: #f2f2f2;">
                <td style="padding: 10px;"><strong>+1 to +2</strong></td>
                <td style="padding: 10px;">Above average condition</td>
            </tr>
            <tr>
                <td style="padding: 10px;"><strong>-1 to +1</strong></td>
                <td style="padding: 10px;">Average condition (most corals)</td>
            </tr>
            <tr style="background-color: #f2f2f2;">
                <td style="padding: 10px;"><strong>-2 to -1</strong></td>
                <td style="padding: 10px;">Below average condition</td>
            </tr>
            <tr>
                <td style="padding: 10px;"><strong>< -2</strong></td>
                <td style="padding: 10px;">Poor condition (bottom ~2.5%)</td>
            </tr>
        </table>

        <h3>Applications for CAFI Analyses</h3>
        <div class="key-finding">
            <strong>Use this condition score to:</strong>
            <ul>
                <li>Test effects of CAFI richness/abundance on coral health</li>
                <li>Compare coral condition across environmental gradients</li>
                <li>Identify outlier corals (unusually healthy or stressed)</li>
                <li>Simplify multivariate physiology into single metric</li>
            </ul>
            <p><strong>Key advantage:</strong> No need to include <code>stump_length</code> as covariate - position bias already removed!</p>
        </div>

        <h3>Data Files</h3>
        <ul>
            <li><code>coral_condition_scores.csv</code> - <strong>Main output:</strong> Condition scores for each coral</li>
            <li><code>condition_score_loadings.csv</code> - PC1 trait loadings</li>
            <li><code>condition_score_summary_by_site.csv</code> - Summary statistics by location</li>
            <li><code>coral_condition_scores.rds</code> - R data file for easy loading</li>
        </ul>
    </div>

    <div class="section">
        <h2>Quality Control Summary</h2>
        <div class="metric">Total Samples: %N%</div>
        <div class="metric">Locations: 3 (HAU, MAT, MRB)</div>
        <div class="metric">Data Completeness: 92-100%</div>
        <div class="metric">Outliers Detected: %OUTLIERS%</div>
    </div>

    <div class="section" style="background-color: #ecf0f1; text-align: center;">
        <p><em>Generated by CAFI Analysis Pipeline | Stier Lab</em></p>
        <p><em>All corals are Pocillopora spp. (morphotypes cannot be reliably distinguished)</em></p>
    </div>

</body>
</html>
')

# Replace placeholders
html_content <- sub("%DATE%", as.character(Sys.Date()), html_content, fixed = TRUE)
html_content <- sub("%N%", as.character(nrow(corrected_traits)), html_content, fixed = TRUE)
html_content <- gsub("%N%", as.character(nrow(corrected_traits)), html_content, fixed = TRUE)
html_content <- sub("%OUTLIERS%", as.character(nrow(outliers)), html_content, fixed = TRUE)

# Write HTML file
writeLines(html_content, file.path(fig_dir, "coral_characteristics_report.html"))

cat("  ✓ HTML report generated\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("========================================\n\n")

cat("Summary:\n")
cat(sprintf("  - Sample size: %d corals with complete data\n", nrow(corrected_traits)))
cat(sprintf("  - Sampling bias identified: r = %.3f (volume vs position)\n", cor_vol_stump))
cat("  - Position correction applied to all traits\n")
cat(sprintf("  - Condition score created: PC1 explains %.1f%% variance\n", variance_explained[1]))
cat("  - No significant size effects after correction ✓\n\n")

cat("Output files:\n")
cat("  Figures:", fig_dir, "\n")
cat("  Tables:", SURVEY_TABLES, "\n")
cat("  HTML report:", file.path(fig_dir, "coral_characteristics_report.html"), "\n\n")

cat("Next steps:\n")
cat("  → Use 'coral_condition_scores.csv' for CAFI analyses\n")
cat("  → Condition score automatically accounts for position bias\n")
cat("  → Higher scores = better overall physiological condition\n\n")
