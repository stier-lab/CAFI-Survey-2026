#!/usr/bin/env Rscript
# ============================================================================
# 05_fig5_condition_landscape.R - FIGURE 5: Coral Condition vs Landscape
# ============================================================================
#
# PURPOSE: Test whether coral physiological condition varies with size or
# landscape position, INDEPENDENT of CAFI effects. This establishes the
# baseline landscape patterns in coral condition before examining CAFI
# feedbacks in Figure 6.
#
# HYPOTHESES:
#   H1: Larger corals may show higher biomass but potentially lower symbiont
#       density (common density-dependent pattern)
#   H2: Isolated corals may experience different condition depending on
#       CAFI-mediated processes
#   H3: PCA summarizing condition metrics reveals overall quality gradients
#
# PHYSIOLOGICAL METRICS (position-corrected):
#   - Protein content (mg/cm2)
#   - Carbohydrate content (mg/cm2)
#   - Zooxanthellae density (cells/cm2)
#   - AFDW tissue biomass (mg/cm2)
#   - Condition Score = PC1 of above (higher = better)
#
# STATISTICAL APPROACHES:
#   - Linear models for condition metrics vs landscape predictors
#   - All analyses use position-corrected traits (sampling bias removed)
#   - PCA biplot to visualize trait covariance
#
# Author: CAFI Analysis Pipeline
# Date: 2025-12-03
# ============================================================================

cat("\n========================================\n")
cat("FIGURE 5: Coral Condition vs Landscape\n")
cat("========================================\n\n")

# Load setup and data
source(here::here("scripts/publication/00_setup.R"))

# Load pre-processed data
# coral_master already contains condition_score from 01_load_data.R
coral_master <- load_object("coral_master.rds")

# Also load condition_scores for the corrected trait components (for PCA biplot)
condition_scores <- load_object("condition_scores.rds")

# Output directory
FIG_DIR <- FIGURE_DIRS$fig5
cat("Outputs will be saved to:", FIG_DIR, "\n\n")

# ============================================================================
# PREPARE ANALYSIS DATA
# ============================================================================

cat("Preparing analysis data...\n")

# Use coral_master which already has condition_score
# Add corrected traits from condition_scores for PCA visualization
analysis_data <- coral_master %>%
  left_join(
    condition_scores %>% select(coral_id, protein_corr, carb_corr, zoox_corr, afdw_corr),
    by = "coral_id"
  ) %>%
  filter(!is.na(condition_score), !is.na(volume)) %>%
  mutate(
    log_volume = log10(volume),
    proximity_m = mean_neighbor_dist / 100,

    # Standardize for models
    volume_z = scale(log_volume)[, 1],
    proximity_z = scale(proximity_m)[, 1],
    neighbor_count_z = scale(n_neighbors)[, 1],

    # Size class for visualization
    size_class = cut(volume,
                     breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE)
  )

# Handle any potential column conflicts from joins
if ("site.x" %in% names(analysis_data)) {
  analysis_data <- analysis_data %>%
    mutate(site = coalesce(site.x, site.y)) %>%
    select(-site.x, -site.y)
}

cat("  - Sample size:", nrow(analysis_data), "corals with condition data\n")
cat("  - Sites:", paste(unique(analysis_data$site), collapse = ", "), "\n\n")

# ============================================================================
# ANALYSIS 1: CONDITION vs CORAL SIZE
# ============================================================================

cat("1. Testing coral condition ~ coral size...\n")

model_cond_size <- lm(condition_score ~ log_volume + site, data = analysis_data)
coef_size <- extract_model_stats(model_cond_size, "log_volume")

cat(sprintf("  Size effect on condition: %.3f (SE = %.3f)\n",
            coef_size$estimate, coef_size$se))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n", coef_size$conf.low, coef_size$conf.high))
cat(sprintf("  P-value %s\n", format_p(coef_size$p.value)))

# Figure 5A: Condition vs Size
show_line_5a <- coef_size$p.value < 0.10

p_5a <- ggplot(analysis_data, aes(x = volume, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  {if (show_line_5a)
    geom_smooth(method = "lm", formula = y ~ log10(x),
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_x_log10(labels = scales::comma) +
  scale_color_site() +
  labs(
    title = "A. Coral Condition vs Colony Size",
    subtitle = sprintf("Size effect = %.2f, p %s | Position-corrected condition score",
                       coef_size$estimate, format_p(coef_size$p.value)),
    x = expression(Coral~Volume~(cm^3~","~log~scale)),
    y = "Condition Score (higher = better)",
    color = "Site"
  ) +
  theme_publication()

save_pub_fig(p_5a, "fig5a_condition_vs_size.png", dir = FIG_DIR)

# ============================================================================
# ANALYSIS 2: CONDITION vs NEIGHBOR PROXIMITY
# ============================================================================

cat("\n2. Testing coral condition ~ neighbor proximity...\n")

model_cond_prox <- lm(condition_score ~ proximity_m + log_volume + site,
                      data = analysis_data)
coef_prox <- extract_model_stats(model_cond_prox, "proximity_m")

cat(sprintf("  Proximity effect on condition: %.3f (SE = %.3f)\n",
            coef_prox$estimate, coef_prox$se))
cat(sprintf("  P-value %s\n", format_p(coef_prox$p.value)))

# Figure 5B: Condition vs Proximity
show_line_5b <- coef_prox$p.value < 0.10

p_5b <- ggplot(analysis_data, aes(x = proximity_m, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  {if (show_line_5b)
    geom_smooth(method = "lm", formula = y ~ x,
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_site() +
  labs(
    title = "B. Coral Condition vs Neighbor Proximity",
    subtitle = sprintf("Proximity effect = %.3f, p %s | Controlling for size",
                       coef_prox$estimate, format_p(coef_prox$p.value)),
    x = "Mean Distance to Neighbors (m)",
    y = "Condition Score (higher = better)",
    color = "Site"
  ) +
  theme_publication()

save_pub_fig(p_5b, "fig5b_condition_vs_proximity.png", dir = FIG_DIR)

# ============================================================================
# ANALYSIS 3: CONDITION vs NEIGHBOR COUNT
# ============================================================================

cat("\n3. Testing coral condition ~ neighbor count...\n")

model_cond_neighbors <- lm(condition_score ~ n_neighbors + log_volume + site,
                           data = analysis_data)
coef_neighbors <- extract_model_stats(model_cond_neighbors, "n_neighbors")

cat(sprintf("  Neighbor count effect: %.3f (SE = %.3f)\n",
            coef_neighbors$estimate, coef_neighbors$se))
cat(sprintf("  P-value %s\n", format_p(coef_neighbors$p.value)))

# Figure 5C: Condition vs Neighbor Count
show_line_5c <- coef_neighbors$p.value < 0.10

p_5c <- ggplot(analysis_data, aes(x = n_neighbors, y = condition_score)) +
  geom_jitter(aes(color = site), alpha = 0.6, size = 3, width = 0.2) +
  {if (show_line_5c)
    geom_smooth(method = "lm", formula = y ~ x,
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_site() +
  labs(
    title = "C. Coral Condition vs Neighbor Count",
    subtitle = sprintf("Neighbor effect = %.3f, p %s | Controlling for size",
                       coef_neighbors$estimate, format_p(coef_neighbors$p.value)),
    x = "Number of Neighboring Corals",
    y = "Condition Score (higher = better)",
    color = "Site"
  ) +
  theme_publication()

save_pub_fig(p_5c, "fig5c_condition_vs_neighbor_count.png", dir = FIG_DIR)

# ============================================================================
# ANALYSIS 4: CONDITION BY SITE
# ============================================================================

cat("\n4. Testing site differences in coral condition...\n")

# ANOVA for site effect
model_site <- lm(condition_score ~ site, data = analysis_data)
anova_site <- anova(model_site)

cat(sprintf("  Site effect: F = %.2f, p %s\n",
            anova_site$`F value`[1], format_p(anova_site$`Pr(>F)`[1])))

# Post-hoc comparisons
if (anova_site$`Pr(>F)`[1] < 0.05) {
  emm <- emmeans(model_site, "site")
  pairs_result <- pairs(emm)
  cat("  Post-hoc pairwise comparisons:\n")
  print(pairs_result)
}

# Figure 5D: Condition by Site
p_5d <- ggplot(analysis_data, aes(x = site, y = condition_score, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_fill_site() +
  labs(
    title = "D. Coral Condition by Site",
    subtitle = sprintf("ANOVA: F = %.2f, p %s",
                       anova_site$`F value`[1], format_p(anova_site$`Pr(>F)`[1])),
    x = "Site",
    y = "Condition Score (higher = better)"
  ) +
  theme_publication() +
  theme(legend.position = "none")

save_pub_fig(p_5d, "fig5d_condition_by_site.png", dir = FIG_DIR)

# ============================================================================
# ANALYSIS 5: INDIVIDUAL PHYSIOLOGICAL TRAITS
# ============================================================================

cat("\n5. Testing individual traits vs landscape...\n")

# If individual corrected traits are available
trait_vars <- c("protein_corr", "carb_corr", "zoox_corr", "afdw_corr")
trait_labels <- c("Protein", "Carbohydrate", "Zooxanthellae", "AFDW")

trait_results <- list()

for (i in seq_along(trait_vars)) {
  var <- trait_vars[i]
  label <- trait_labels[i]

  # Check if trait is in analysis_data (already joined earlier)
  if (var %in% names(analysis_data)) {
    # Filter to non-NA values for this trait
    trait_data <- analysis_data %>%
      filter(!is.na(.data[[var]]))

    if (nrow(trait_data) > 20) {
      # Model: trait ~ size + proximity + site
      model <- lm(reformulate(c("log_volume", "proximity_m", "site"), response = var),
                  data = trait_data)

      coef_size_trait <- extract_model_stats(model, "log_volume")
      coef_prox_trait <- extract_model_stats(model, "proximity_m")

      trait_results[[label]] <- tibble(
        trait = label,
        size_effect = coef_size_trait$estimate,
        size_p = coef_size_trait$p.value,
        prox_effect = coef_prox_trait$estimate,
        prox_p = coef_prox_trait$p.value
      )

      cat(sprintf("  %s: size = %.3f (p %s), prox = %.3f (p %s)\n",
                  label,
                  coef_size_trait$estimate, format_p(coef_size_trait$p.value),
                  coef_prox_trait$estimate, format_p(coef_prox_trait$p.value)))
    }
  }
}

trait_results_df <- bind_rows(trait_results)

# Figure 5E: Individual trait effects
if (nrow(trait_results_df) > 0) {

  trait_plot_data <- trait_results_df %>%
    pivot_longer(cols = c(size_effect, prox_effect),
                 names_to = "predictor",
                 values_to = "effect") %>%
    mutate(
      p_value = ifelse(predictor == "size_effect", size_p, prox_p),
      significant = p_value < 0.05,
      predictor = case_when(
        predictor == "size_effect" ~ "Coral Size",
        predictor == "prox_effect" ~ "Neighbor Proximity"
      )
    )

  p_5e <- ggplot(trait_plot_data, aes(x = effect, y = trait, color = predictor)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(shape = significant), size = 4,
               position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("Coral Size" = "#0072B2",
                                  "Neighbor Proximity" = "#D55E00")) +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16),
                       labels = c("Not significant", "p < 0.05")) +
    labs(
      title = "E. Landscape Effects on Individual Traits",
      subtitle = "Position-corrected physiological traits",
      x = "Effect Size",
      y = "Physiological Trait",
      color = "Predictor",
      shape = "Significance"
    ) +
    theme_publication()

  save_pub_fig(p_5e, "fig5e_individual_traits.png", dir = FIG_DIR)
}

# ============================================================================
# ANALYSIS 6: FULL LANDSCAPE MODEL
# ============================================================================

cat("\n6. Fitting full landscape model...\n")

# Full model with all predictors
model_full <- lm(condition_score ~ log_volume + proximity_m + n_neighbors +
                   total_neighbor_volume + site,
                 data = analysis_data %>%
                   filter(!is.na(total_neighbor_volume)))

# Check VIF for multicollinearity
vif_values <- car::vif(model_full)
cat("  Variance Inflation Factors:\n")
print(round(vif_values, 2))

# Model summary
full_summary <- summary(model_full)
cat(sprintf("\n  Model R-squared: %.3f\n", full_summary$r.squared))
cat(sprintf("  Adjusted R-squared: %.3f\n", full_summary$adj.r.squared))

# ============================================================================
# COMBINED FIGURE 5 PANEL
# ============================================================================

cat("\n7. Creating combined Figure 5 panel...\n")

fig5_combined <- (p_5a + p_5b) /
  (p_5c + p_5d) +
  plot_annotation(
    title = "Figure 5. Coral Condition Across Landscape Characteristics",
    subtitle = "Position-corrected physiological condition vs coral size, proximity, and site",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_pub_fig(fig5_combined, "fig5_condition_landscape_combined.png",
             dir = FIG_DIR, width = 14, height = 12)

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n8. Saving results...\n")

# Compile model results
fig5_results <- tibble(
  predictor = c("Coral Size (log)", "Proximity", "Neighbor Count"),
  estimate = c(coef_size$estimate, coef_prox$estimate, coef_neighbors$estimate),
  std_error = c(coef_size$se, coef_prox$se, coef_neighbors$se),
  conf_low = c(coef_size$conf.low, coef_prox$conf.low, coef_neighbors$conf.low),
  conf_high = c(coef_size$conf.high, coef_prox$conf.high, coef_neighbors$conf.high),
  p_value = c(coef_size$p.value, coef_prox$p.value, coef_neighbors$p.value),
  significant = c(coef_size$significant, coef_prox$significant, coef_neighbors$significant)
)

save_table(fig5_results, "fig5_landscape_condition_results.csv")
if (nrow(trait_results_df) > 0) {
  save_table(trait_results_df, "fig5_individual_trait_results.csv")
}

save_object(list(
  model_size = model_cond_size,
  model_prox = model_cond_prox,
  model_neighbors = model_cond_neighbors,
  model_full = model_full
), "fig5_models.rds")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("Figure 5 Analysis Summary\n")
cat("========================================\n\n")

cat("LANDSCAPE EFFECTS ON CORAL CONDITION:\n\n")

cat("1. CORAL SIZE:\n")
if (coef_size$significant) {
  if (coef_size$estimate > 0) {
    cat("   -> Larger corals have BETTER condition\n")
  } else {
    cat("   -> Larger corals have WORSE condition\n")
  }
} else {
  cat("   -> No significant size-condition relationship\n")
}
cat(sprintf("   Effect = %.3f, p %s\n\n", coef_size$estimate, format_p(coef_size$p.value)))

cat("2. NEIGHBOR PROXIMITY:\n")
if (coef_prox$significant) {
  if (coef_prox$estimate > 0) {
    cat("   -> Isolated corals have BETTER condition\n")
  } else {
    cat("   -> Clustered corals have BETTER condition\n")
  }
} else {
  cat("   -> No significant proximity-condition relationship\n")
}
cat(sprintf("   Effect = %.3f, p %s\n\n", coef_prox$estimate, format_p(coef_prox$p.value)))

cat("3. NEIGHBOR COUNT:\n")
if (coef_neighbors$significant) {
  cat(sprintf("   -> Each additional neighbor associated with %.3f change in condition\n",
              coef_neighbors$estimate))
} else {
  cat("   -> No significant neighbor count effect\n")
}
cat(sprintf("   Effect = %.3f, p %s\n\n", coef_neighbors$estimate, format_p(coef_neighbors$p.value)))

cat("INTERPRETATION:\n")
cat("   These landscape effects on condition are INDEPENDENT of CAFI.\n")
cat("   Figure 6 tests whether CAFI mediate additional condition effects.\n\n")

cat("Figures saved to:", FIG_DIR, "\n")
cat("Tables saved to:", TABLES_DIR, "\n\n")
cat("Figure 5 analysis complete!\n\n")
