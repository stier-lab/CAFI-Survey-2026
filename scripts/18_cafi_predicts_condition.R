#!/usr/bin/env Rscript
# ============================================================================
# 18_cafi_predicts_condition.R - Do CAFI Predict Coral Condition?
#
# Purpose: Test whether CAFI community attributes predict coral physiological
#          condition (reverse causation from H4)
#
# Hypotheses tested:
#   - Does total CAFI abundance predict coral condition?
#   - Do specific taxonomic groups (crabs, shrimps, fish, snails) predict condition?
#   - Does CAFI community composition (PCA1) predict condition (PCA1)?
#   - Are there synergistic/antagonistic effects of different groups?
#
# Theoretical Background:
#   While H4 asks if healthy corals attract more CAFI, this script asks the
#   reverse: do CAFI benefit coral health? Possible mechanisms:
#   - Mutualistic crabs remove sediment and defend against predators
#   - Fish provide nutrient subsidies via waste
#   - Some associates may be parasitic or neutral
#
# Author: CAFI Analysis Pipeline
# Date: 2025-11-23
# ============================================================================

cat("============================================================\n")
cat("Script 18: Do CAFI Predict Coral Condition?\n")
cat("============================================================\n\n")

# ============================================================================
# Setup
# ============================================================================

# Load libraries and theme
source(here::here("scripts/00_load_libraries.R"))

# Load data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))
condition_scores <- readRDS(file.path(SURVEY_OBJECTS, "coral_condition_scores.rds"))

# Create output directory
fig_dir <- file.path(SURVEY_FIGURES, "cafi_predicts_condition")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("Data loaded successfully\n")
cat("  - Condition scores:", nrow(condition_scores), "corals\n")
cat("  - Community matrix:", nrow(community_matrix), "corals x", ncol(community_matrix), "OTUs\n\n")

# ============================================================================
# 1. PREPARE PREDICTOR VARIABLES
# ============================================================================

cat("1. Preparing CAFI predictor variables...\n")

# Calculate taxonomic group abundances per coral
# Note: using lowest_level as species identifier since otu_id doesn't exist
taxonomic_predictors <- cafi_clean %>%
  group_by(coral_id) %>%
  summarise(
    total_cafi = n(),
    n_crabs = sum(type == "crab", na.rm = TRUE),
    n_shrimps = sum(type == "shrimp", na.rm = TRUE),
    n_fish = sum(type == "fish", na.rm = TRUE),
    n_snails = sum(type == "snail", na.rm = TRUE),
    otu_richness = n_distinct(lowest_level),
    shannon = vegan::diversity(table(lowest_level), index = "shannon"),
    .groups = "drop"
  )

cat("  ✓ Calculated abundances for", nrow(taxonomic_predictors), "corals\n")
cat("    - Crabs range:", min(taxonomic_predictors$n_crabs), "-", max(taxonomic_predictors$n_crabs), "\n")
cat("    - Shrimps range:", min(taxonomic_predictors$n_shrimps), "-", max(taxonomic_predictors$n_shrimps), "\n")
cat("    - Fish range:", min(taxonomic_predictors$n_fish), "-", max(taxonomic_predictors$n_fish), "\n")
cat("    - Snails range:", min(taxonomic_predictors$n_snails), "-", max(taxonomic_predictors$n_snails), "\n\n")

# ============================================================================
# 2. PCA OF CAFI COMMUNITY
# ============================================================================

cat("2. Running PCA on CAFI community composition...\n")

# Hellinger transformation for PCA (recommended for community data)
community_hell <- decostand(community_matrix, method = "hellinger")

# Run PCA
pca_cafi <- prcomp(community_hell, scale. = TRUE)

# Extract scores
cafi_pc_scores <- data.frame(
  coral_id = rownames(community_matrix),
  cafi_PC1 = pca_cafi$x[, 1],
  cafi_PC2 = pca_cafi$x[, 2]
)

# Variance explained
var_explained_cafi <- summary(pca_cafi)$importance[2, 1:2] * 100
cat("  ✓ Community PCA complete\n")
cat("    - PC1:", round(var_explained_cafi[1], 1), "% variance\n")
cat("    - PC2:", round(var_explained_cafi[2], 1), "% variance\n\n")

# ============================================================================
# 3. MERGE ALL DATA
# ============================================================================

cat("3. Merging predictor and response data...\n")

# Merge all datasets
# Note: volume_lab is coral volume from survey_master
analysis_data <- condition_scores %>%
  select(coral_id, condition_score, site) %>%
  left_join(taxonomic_predictors, by = "coral_id") %>%
  left_join(cafi_pc_scores, by = "coral_id") %>%
  left_join(survey_master %>% select(coral_id, coral_volume = volume_lab,
                                      morphotype, depth_m), by = "coral_id") %>%
  filter(!is.na(condition_score), !is.na(total_cafi))

cat("  ✓ Merged data:", nrow(analysis_data), "corals with complete data\n\n")

# ============================================================================
# 4. UNIVARIATE ANALYSES: TAXONOMIC GROUPS
# ============================================================================

cat("4. Testing taxonomic group predictors of coral condition...\n\n")

# Function to run and report regression
run_condition_model <- function(data, predictor_name, predictor_col) {
  formula <- as.formula(paste("condition_score ~", predictor_col))
  model <- lm(formula, data = data)
  summary_model <- summary(model)

  result <- data.frame(
    predictor = predictor_name,
    estimate = coef(model)[2],
    se = summary_model$coefficients[2, 2],
    t_value = summary_model$coefficients[2, 3],
    p_value = summary_model$coefficients[2, 4],
    r_squared = summary_model$r.squared,
    adj_r_squared = summary_model$adj.r.squared,
    f_statistic = summary_model$fstatistic[1],
    n = nrow(data)
  )

  return(list(model = model, result = result))
}

# Test each predictor
predictors <- list(
  c("Total CAFI", "total_cafi"),
  c("Crab abundance", "n_crabs"),
  c("Shrimp abundance", "n_shrimps"),
  c("Fish abundance", "n_fish"),
  c("Snail abundance", "n_snails"),
  c("OTU richness", "otu_richness"),
  c("Shannon diversity", "shannon"),
  c("CAFI PC1", "cafi_PC1")
)

univariate_results <- list()
univariate_models <- list()

for (pred in predictors) {
  result <- run_condition_model(analysis_data, pred[1], pred[2])
  univariate_results[[pred[1]]] <- result$result
  univariate_models[[pred[1]]] <- result$model
}

# Combine results
univariate_df <- bind_rows(univariate_results)

# Print results
cat("  Univariate regression results (Condition ~ Predictor):\n\n")
univariate_df %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ ".",
      TRUE ~ " "
    )
  ) %>%
  select(predictor, estimate, t_value, p_value, r_squared, sig)

print(as.data.frame(univariate_df %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ ".",
      TRUE ~ " "
    )
  ) %>%
  select(predictor, estimate, t_value, p_value, r_squared, sig)))

cat("\n")

# Save results
write_csv(univariate_df, file.path(SURVEY_TABLES, "cafi_predicts_condition_univariate.csv"))

# ============================================================================
# 5. VISUALIZATION: UNIVARIATE RELATIONSHIPS
# ============================================================================

cat("5. Creating univariate relationship plots...\n")

# Panel plot for all taxonomic groups
p_crabs <- ggplot(analysis_data, aes(x = n_crabs, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  scale_color_site() +
  labs(title = "Crabs", x = "Number of crabs", y = "Condition score") +
  theme_publication() +
  theme(legend.position = "none")

p_shrimps <- ggplot(analysis_data, aes(x = n_shrimps, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  scale_color_site() +
  labs(title = "Shrimps", x = "Number of shrimps", y = "Condition score") +
  theme_publication() +
  theme(legend.position = "none")

p_fish <- ggplot(analysis_data, aes(x = n_fish, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  scale_color_site() +
  labs(title = "Fish", x = "Number of fish", y = "Condition score") +
  theme_publication() +
  theme(legend.position = "none")

p_snails <- ggplot(analysis_data, aes(x = n_snails, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  scale_color_site() +
  labs(title = "Snails", x = "Number of snails", y = "Condition score") +
  theme_publication() +
  theme(legend.position = "none")

# Combine with patchwork
p_taxonomic_panel <- (p_crabs | p_shrimps) / (p_fish | p_snails) +
  plot_annotation(
    title = "Taxonomic Group Abundances as Predictors of Coral Condition",
    subtitle = "Linear regressions: condition score ~ group abundance",
    caption = "Condition = position-corrected PC1 of protein, carbs, zoox, AFDW"
  )

ggsave(file.path(fig_dir, "taxonomic_groups_vs_condition.png"),
       p_taxonomic_panel, width = 10, height = 8, dpi = 300, bg = "white")

cat("  ✓ Saved taxonomic groups panel\n")

# Total CAFI vs condition
p_total <- ggplot(analysis_data, aes(x = total_cafi, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_site() +
  scale_x_sqrt(breaks = c(0, 25, 100, 225)) +
  labs(
    title = "Total CAFI Abundance as Predictor of Coral Condition",
    subtitle = sprintf("R² = %.3f, p = %.3f",
                       univariate_df$r_squared[univariate_df$predictor == "Total CAFI"],
                       univariate_df$p_value[univariate_df$predictor == "Total CAFI"]),
    x = "Total CAFI abundance (sqrt scale)",
    y = "Coral condition score"
  ) +
  theme_publication()

ggsave(file.path(fig_dir, "total_cafi_vs_condition.png"),
       p_total, width = 10, height = 7, dpi = 300, bg = "white")

cat("  ✓ Saved total CAFI plot\n")

# ============================================================================
# 6. PCA1 vs PCA1: COMMUNITY vs CONDITION
# ============================================================================

cat("6. Testing community composition (PC1) vs condition (PC1)...\n")

# Plot PCA1 vs PCA1
p_pca <- ggplot(analysis_data, aes(x = cafi_PC1, y = condition_score)) +
  geom_point(aes(color = site, size = total_cafi), alpha = 0.6) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_site() +
  scale_size_continuous(range = c(2, 8), name = "CAFI\nabundance") +
  labs(
    title = "CAFI Community Composition vs Coral Condition",
    subtitle = sprintf("Community PC1 (%.1f%% var) vs Condition PC1 | R² = %.3f, p = %.3f",
                       var_explained_cafi[1],
                       univariate_df$r_squared[univariate_df$predictor == "CAFI PC1"],
                       univariate_df$p_value[univariate_df$predictor == "CAFI PC1"]),
    x = paste0("CAFI Community PC1 (", round(var_explained_cafi[1], 1), "%)"),
    y = "Coral Condition Score (PC1)"
  ) +
  theme_publication()

ggsave(file.path(fig_dir, "community_pc1_vs_condition_pc1.png"),
       p_pca, width = 10, height = 8, dpi = 300, bg = "white")

cat("  ✓ Saved PCA1 vs PCA1 plot\n")

# ============================================================================
# 7. MULTIPLE REGRESSION: ALL GROUPS TOGETHER
# ============================================================================

cat("7. Running multiple regression with all taxonomic groups...\n")

# Full model with all groups
model_full <- lm(condition_score ~ n_crabs + n_shrimps + n_fish + n_snails +
                   coral_volume + site,
                 data = analysis_data)

cat("\n  Full model: Condition ~ Crabs + Shrimps + Fish + Snails + Volume + Site\n")
print(summary(model_full))

# Model without coral volume (to see pure CAFI effects)
model_cafi_only <- lm(condition_score ~ n_crabs + n_shrimps + n_fish + n_snails,
                      data = analysis_data)

cat("\n  CAFI-only model: Condition ~ Crabs + Shrimps + Fish + Snails\n")
print(summary(model_cafi_only))

# Save model summaries
capture.output(
  {
    cat("=== FULL MODEL ===\n")
    print(summary(model_full))
    cat("\n\n=== CAFI-ONLY MODEL ===\n")
    print(summary(model_cafi_only))
  },
  file = file.path(SURVEY_TABLES, "cafi_predicts_condition_multiple_regression.txt")
)

# ============================================================================
# 8. COEFFICIENT PLOT
# ============================================================================

cat("8. Creating coefficient comparison plot...\n")

# Extract coefficients for visualization
coef_df <- univariate_df %>%
  filter(predictor %in% c("Crab abundance", "Shrimp abundance",
                          "Fish abundance", "Snail abundance")) %>%
  mutate(
    ci_lower = estimate - 1.96 * se,
    ci_upper = estimate + 1.96 * se,
    significant = p_value < 0.05,
    predictor = factor(predictor,
                       levels = c("Crab abundance", "Shrimp abundance",
                                  "Fish abundance", "Snail abundance"))
  )

p_coef <- ggplot(coef_df, aes(x = predictor, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, color = significant),
                  size = 1, linewidth = 1) +
  geom_text(aes(label = sprintf("p = %.3f", p_value)),
            hjust = -0.2, vjust = 0.5, size = 3.5) +
  scale_color_manual(values = c("TRUE" = "#2e7d32", "FALSE" = "#757575"),
                     name = "Significant\n(p < 0.05)") +
  coord_flip() +
  labs(
    title = "Effect of Taxonomic Groups on Coral Condition",
    subtitle = "Univariate regression coefficients with 95% CI",
    x = "",
    y = "Regression coefficient (effect on condition score)"
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "coefficient_comparison.png"),
       p_coef, width = 10, height = 6, dpi = 300, bg = "white")

cat("  ✓ Saved coefficient plot\n")

# ============================================================================
# 9. INTERACTIONS: GROUP RATIOS
# ============================================================================

cat("9. Testing group composition ratios...\n")

# Calculate ratios (with small constant to avoid division by zero)
analysis_data <- analysis_data %>%
  mutate(
    prop_crabs = n_crabs / (total_cafi + 1),
    prop_shrimps = n_shrimps / (total_cafi + 1),
    prop_fish = n_fish / (total_cafi + 1),
    prop_snails = n_snails / (total_cafi + 1),
    crab_to_shrimp = (n_crabs + 1) / (n_shrimps + 1)
  )

# Test if proportions matter
model_proportions <- lm(condition_score ~ prop_crabs + prop_shrimps +
                          prop_fish + prop_snails + total_cafi,
                        data = analysis_data)

cat("\n  Proportion model: Condition ~ prop_crabs + prop_shrimps + prop_fish + prop_snails + total\n")
print(summary(model_proportions))

# Crab-to-shrimp ratio (potential competition/complementarity)
model_ratio <- lm(condition_score ~ crab_to_shrimp + total_cafi,
                  data = analysis_data)

cat("\n  Ratio model: Condition ~ crab:shrimp ratio + total_cafi\n")
print(summary(model_ratio))

# ============================================================================
# 10. SITE-SPECIFIC EFFECTS
# ============================================================================

cat("10. Testing site-specific effects...\n")

# Run separate models for each site
site_results <- analysis_data %>%
  group_by(site) %>%
  group_modify(~ {
    if (nrow(.x) >= 10) {
      model <- lm(condition_score ~ n_crabs + n_shrimps + n_fish + n_snails, data = .x)
      tibble(
        n = nrow(.x),
        r_squared = summary(model)$r.squared,
        f_stat = summary(model)$fstatistic[1],
        p_value = pf(summary(model)$fstatistic[1],
                     summary(model)$fstatistic[2],
                     summary(model)$fstatistic[3],
                     lower.tail = FALSE)
      )
    } else {
      tibble(n = nrow(.x), r_squared = NA, f_stat = NA, p_value = NA)
    }
  }) %>%
  ungroup()

cat("\n  Site-specific model results (Condition ~ all groups):\n")
print(site_results)

# Faceted plot by site
p_site_facet <- ggplot(analysis_data, aes(x = total_cafi, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  facet_wrap(~site, scales = "free_x") +
  scale_color_site() +
  scale_x_sqrt() +
  labs(
    title = "CAFI-Condition Relationship by Site",
    subtitle = "Site-specific patterns may reveal environmental moderation",
    x = "Total CAFI abundance (sqrt scale)",
    y = "Condition score"
  ) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "cafi_condition_by_site.png"),
       p_site_facet, width = 12, height = 5, dpi = 300, bg = "white")

cat("  ✓ Saved site-specific plot\n")

# ============================================================================
# 11. SUMMARY STATISTICS TABLE
# ============================================================================

cat("11. Creating summary statistics for manuscript...\n")

# Prepare statistics for extraction
stats_summary <- univariate_df %>%
  mutate(
    analysis = "CAFI Predicts Condition",
    hypothesis = "Reverse H4",
    test_type = "Linear regression",
    response = "Coral condition score",
    stat_type = "F",
    ci_95 = sprintf("[%.3f, %.3f]", estimate - 1.96 * se, estimate + 1.96 * se)
  ) %>%
  select(
    Analysis = analysis,
    Hypothesis = hypothesis,
    Test = test_type,
    Predictor = predictor,
    Response = response,
    `Stat Type` = stat_type,
    Statistic = f_statistic,
    `P-value` = p_value,
    `Effect Type` = predictor,
    `Effect Size` = r_squared,
    `95% CI` = ci_95,
    N = n
  )

write_csv(stats_summary, file.path(SURVEY_TABLES, "cafi_predicts_condition_stats.csv"))

cat("  ✓ Saved statistics summary\n")

# ============================================================================
# 12. COMPREHENSIVE SUMMARY FIGURE
# ============================================================================

cat("12. Creating comprehensive summary figure...\n")

# Combine key plots
p_summary <- (p_total + p_pca) / (p_taxonomic_panel) +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(
    title = "Do CAFI Predict Coral Condition?",
    subtitle = "Testing reverse causation: CAFI abundance/composition as predictors of host physiology",
    caption = "All analyses control for coral volume and site effects where appropriate"
  )

ggsave(file.path(fig_dir, "comprehensive_summary.png"),
       p_summary, width = 14, height = 16, dpi = 300, bg = "white")

cat("  ✓ Saved comprehensive summary figure\n")

# ============================================================================
# 13. FINAL SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE: Do CAFI Predict Coral Condition?\n")
cat("============================================================\n\n")

# Identify significant predictors
sig_predictors <- univariate_df %>%
  filter(p_value < 0.05) %>%
  pull(predictor)

if (length(sig_predictors) > 0) {
  cat("SIGNIFICANT PREDICTORS (p < 0.05):\n")
  for (pred in sig_predictors) {
    row <- univariate_df[univariate_df$predictor == pred, ]
    direction <- ifelse(row$estimate > 0, "positive", "negative")
    cat(sprintf("  - %s: %s effect (β = %.3f, R² = %.3f, p = %.3f)\n",
                pred, direction, row$estimate, row$r_squared, row$p_value))
  }
} else {
  cat("NO SIGNIFICANT PREDICTORS FOUND (p < 0.05)\n\n")
  cat("Key findings:\n")
  # Show best predictor even if not significant
  best_pred <- univariate_df[which.min(univariate_df$p_value), ]
  cat(sprintf("  - Best predictor: %s (p = %.3f, R² = %.3f)\n",
              best_pred$predictor, best_pred$p_value, best_pred$r_squared))
}

cat("\nInterpretation:\n")
cat("  - If no significant predictors: CAFI do not predict coral condition\n")
cat("  - This suggests the relationship is unidirectional (if any)\n")
cat("  - Combined with H4 null result, CAFI-condition relationship may be truly absent\n")

cat("\nOutput saved to:\n")
cat("  - Figures:", fig_dir, "\n")
cat("  - Tables:", SURVEY_TABLES, "\n\n")

cat("✅ Script 18 complete!\n")
