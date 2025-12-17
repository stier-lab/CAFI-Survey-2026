#!/usr/bin/env Rscript
# ============================================================================
# 06_fig6_cafi_feedbacks.R - FIGURE 6: CAFI -> Coral Condition Feedbacks
# ============================================================================
#
# PURPOSE: Directly test whether CAFI communities influence coral condition.
# This is the "feedback" component of the CAFI-coral relationship.
#
# HYPOTHESES:
#   H1: More defenders (Trapezia) -> higher condition
#   H2: More resident fishes -> improved condition via nutrient provisioning
#   H3: More corallivores (Drupella) -> lower condition
#   H4: CAFI community composition (PC) predicts condition more strongly
#       than abundance alone
#   H5: Landscape may modulate these effects (e.g., mutualist benefits
#       weaken on very large colonies)
#
# STATISTICAL APPROACHES:
#   - Models: Condition ~ CAFI abundance + diversity (controlling for size)
#   - Functional group models: Condition ~ defenders + fish + corallivores
#   - Composition effects: Condition ~ CAFI composition PCA scores
#   - Interaction tests: Condition ~ CAFI * size
#   - Visualization: Lollipop plots of effect sizes, partial residual plots
#
# Author: CAFI Analysis Pipeline
# Date: 2025-12-03
# ============================================================================

cat("\n========================================\n")
cat("FIGURE 6: CAFI -> Coral Condition Feedbacks\n")
cat("========================================\n\n")

# Load setup and data
source(here::here("scripts/publication/00_setup.R"))

# Load pre-processed data
# coral_master already contains condition_score from 01_load_data.R
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")
community_matrix <- load_object("community_matrix.rds")

# Output directory
FIG_DIR <- FIGURE_DIRS$fig6
cat("Outputs will be saved to:", FIG_DIR, "\n\n")

# ============================================================================
# PREPARE ANALYSIS DATA
# ============================================================================

cat("Preparing analysis data...\n")

# Use coral_master which already has condition_score
analysis_data <- coral_master %>%
  filter(!is.na(condition_score), !is.na(volume), total_cafi > 0) %>%
  mutate(
    log_volume = log10(volume),
    log_cafi = log(total_cafi + 1),

    # Standardize predictors
    volume_z = scale(log_volume)[, 1],
    cafi_z = scale(total_cafi)[, 1],
    richness_z = scale(otu_richness)[, 1],

    # Size class
    size_class = cut(volume,
                     breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE)
  )

cat("  - Sample size:", nrow(analysis_data), "corals with CAFI and condition data\n\n")

# ============================================================================
# ANALYSIS 1: TOTAL CAFI ABUNDANCE -> CONDITION
# ============================================================================

cat("1. Testing condition ~ total CAFI abundance...\n")

model_abund <- lm(condition_score ~ total_cafi + log_volume + site,
                  data = analysis_data)
coef_abund <- extract_model_stats(model_abund, "total_cafi")

cat(sprintf("  CAFI abundance effect: %.4f (SE = %.4f)\n",
            coef_abund$estimate, coef_abund$se))
cat(sprintf("  P-value %s\n", format_p(coef_abund$p.value)))

# Figure 6A: Condition vs CAFI Abundance
show_line_6a <- coef_abund$p.value < 0.10

p_6a <- ggplot(analysis_data, aes(x = total_cafi, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  {if (show_line_6a)
    geom_smooth(method = "lm", formula = y ~ x,
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_x_sqrt(breaks = c(0, 25, 100, 225, 400)) +
  scale_color_site() +
  labs(
    title = "A. Coral Condition vs CAFI Abundance",
    subtitle = sprintf("Effect = %.4f, p %s | Controlling for coral size + site",
                       coef_abund$estimate, format_p(coef_abund$p.value)),
    x = "Total CAFI Abundance (sqrt scale)",
    y = "Condition Score (higher = better)",
    color = "Site"
  ) +
  theme_publication()

save_pub_fig(p_6a, "fig6a_condition_vs_abundance.png", dir = FIG_DIR)

# ============================================================================
# ANALYSIS 2: CAFI RICHNESS -> CONDITION
# ============================================================================

cat("\n2. Testing condition ~ CAFI richness...\n")

model_rich <- lm(condition_score ~ otu_richness + log_volume + site,
                 data = analysis_data)
coef_rich <- extract_model_stats(model_rich, "otu_richness")

cat(sprintf("  CAFI richness effect: %.4f (SE = %.4f)\n",
            coef_rich$estimate, coef_rich$se))
cat(sprintf("  P-value %s\n", format_p(coef_rich$p.value)))

# Figure 6B: Condition vs CAFI Richness
show_line_6b <- coef_rich$p.value < 0.10

p_6b <- ggplot(analysis_data, aes(x = otu_richness, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  {if (show_line_6b)
    geom_smooth(method = "lm", formula = y ~ x,
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_site() +
  labs(
    title = "B. Coral Condition vs CAFI Richness",
    subtitle = sprintf("Effect = %.3f, p %s | Controlling for size + site",
                       coef_rich$estimate, format_p(coef_rich$p.value)),
    x = "OTU Richness",
    y = "Condition Score (higher = better)",
    color = "Site"
  ) +
  theme_publication()

save_pub_fig(p_6b, "fig6b_condition_vs_richness.png", dir = FIG_DIR)

# ============================================================================
# ANALYSIS 3: FUNCTIONAL GROUP EFFECTS
# ============================================================================

cat("\n3. Testing condition ~ functional group abundances...\n")

# Check which functional groups have non-zero variation
functional_groups_all <- c("n_trapezia", "n_resident_fish", "n_corallivore",
                           "n_other_crab", "n_shrimp")
functional_labels_all <- c("Trapezia (Defenders)", "Resident Fish",
                           "Corallivores", "Other Crabs", "Shrimp")

# Filter to groups with some presence
groups_with_data <- functional_groups_all[
  sapply(functional_groups_all, function(g) {
    sum(analysis_data[[g]] > 0, na.rm = TRUE) >= 5
  })
]

# Model with available functional groups
formula_terms <- c(groups_with_data, "log_volume", "site")
model_functional <- lm(
  reformulate(formula_terms, response = "condition_score"),
  data = analysis_data
)

# Get labels for groups with data
functional_groups <- groups_with_data
functional_labels <- functional_labels_all[functional_groups_all %in% groups_with_data]

functional_results <- list()

for (i in seq_along(functional_groups)) {
  grp <- functional_groups[i]
  label <- functional_labels[i]

  coef_info <- extract_model_stats(model_functional, grp)

  if (!is.null(coef_info)) {
    functional_results[[grp]] <- tibble(
      group = label,
      estimate = coef_info$estimate,
      se = coef_info$se,
      conf_low = coef_info$conf.low,
      conf_high = coef_info$conf.high,
      p_value = coef_info$p.value,
      significant = coef_info$significant
    )

    sig <- if (coef_info$significant) " *" else ""
    cat(sprintf("  %s: %.4f (p %s)%s\n",
                label, coef_info$estimate, format_p(coef_info$p.value), sig))
  }
}

functional_df <- bind_rows(functional_results)

# Figure 6C: Functional Group Effect Sizes (Lollipop Plot)
p_6c <- functional_df %>%
  mutate(group = factor(group, levels = rev(group))) %>%
  ggplot(aes(x = estimate, y = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_segment(aes(x = 0, xend = estimate, yend = group),
               linewidth = 1, color = "gray50") +
  geom_point(aes(color = significant), size = 5) +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high),
                 height = 0.2, linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "#009E73"),
                     labels = c("Not significant", "p < 0.05"),
                     name = "") +
  labs(
    title = "C. Functional Group Effects on Coral Condition",
    subtitle = "Effect of each group controlling for all others + coral size + site",
    x = "Effect on Condition Score",
    y = ""
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

save_pub_fig(p_6c, "fig6c_functional_group_effects.png", dir = FIG_DIR,
             width = 10, height = 6)

# ============================================================================
# ANALYSIS 4: TRAPEZIA FOCUS (KEY MUTUALIST)
# ============================================================================

cat("\n4. Focused analysis: Trapezia (key mutualist)...\n")

# Trapezia presence/absence effect
analysis_data <- analysis_data %>%
  mutate(has_trapezia = n_trapezia > 0)

# Check if there's variation in Trapezia presence
if (sum(analysis_data$has_trapezia) > 5 && sum(!analysis_data$has_trapezia) > 5) {

  model_trap_presence <- lm(condition_score ~ has_trapezia + log_volume + site,
                            data = analysis_data)
  coef_trap_pres <- extract_model_stats(model_trap_presence, "has_trapeziaTRUE")

  if (!is.null(coef_trap_pres)) {
    cat(sprintf("  Trapezia presence effect: %.3f (p %s)\n",
                coef_trap_pres$estimate, format_p(coef_trap_pres$p.value)))
  } else {
    cat("  Trapezia presence effect: could not be estimated\n")
    coef_trap_pres <- list(estimate = NA, p.value = NA)
  }

  # Trapezia abundance effect (among corals with Trapezia)
  trap_corals <- analysis_data %>% filter(has_trapezia)
  if (nrow(trap_corals) > 10) {
    model_trap_abund <- lm(condition_score ~ n_trapezia + log_volume + site,
                           data = trap_corals)
    coef_trap_abund <- extract_model_stats(model_trap_abund, "n_trapezia")

    if (!is.null(coef_trap_abund)) {
      cat(sprintf("  Trapezia abundance effect (when present): %.3f (p %s)\n",
                  coef_trap_abund$estimate, format_p(coef_trap_abund$p.value)))
    }
  }

  # Figure 6D: Condition by Trapezia Presence
  trap_subtitle <- if (!is.na(coef_trap_pres$estimate)) {
    sprintf("Effect = %.3f, p %s | Controlling for size + site",
            coef_trap_pres$estimate, format_p(coef_trap_pres$p.value))
  } else {
    "Trapezia presence effect not estimable"
  }

  p_6d <- ggplot(analysis_data, aes(x = has_trapezia, y = condition_score,
                                     fill = has_trapezia)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "#CC79A7")) +
    scale_x_discrete(labels = c("FALSE" = "No Trapezia", "TRUE" = "Trapezia Present")) +
    labs(
      title = "D. Coral Condition by Trapezia Presence",
      subtitle = trap_subtitle,
      x = "",
      y = "Condition Score (higher = better)"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  save_pub_fig(p_6d, "fig6d_condition_by_trapezia.png", dir = FIG_DIR)

} else {
  cat("  Insufficient variation in Trapezia presence for analysis\n")
  coef_trap_pres <- list(estimate = NA, se = NA, p.value = NA, significant = FALSE)
}

# ============================================================================
# ANALYSIS 5: SIZE-DEPENDENT CAFI EFFECTS (INTERACTION)
# ============================================================================

cat("\n5. Testing size-dependent CAFI effects (interaction)...\n")

# Does CAFI benefit depend on coral size?
model_interaction <- lm(condition_score ~ total_cafi * log_volume + site,
                        data = analysis_data)

coef_interaction <- extract_model_stats(model_interaction, "total_cafi:log_volume")

if (!is.null(coef_interaction)) {
  cat(sprintf("  CAFI x Size interaction: %.6f (p %s)\n",
              coef_interaction$estimate, format_p(coef_interaction$p.value)))

  if (coef_interaction$significant) {
    cat("  -> CAFI effects VARY with coral size!\n")
  } else {
    cat("  -> CAFI effects do not vary significantly with size\n")
  }
} else {
  cat("  -> Interaction term could not be estimated\n")
  coef_interaction <- list(estimate = NA, se = NA, p.value = NA, significant = FALSE)
}

# Figure 6E: Interaction plot (CAFI effect by size class)
interaction_subtitle <- if (!is.null(coef_interaction)) {
  sprintf("Interaction term: p %s", format_p(coef_interaction$p.value))
} else {
  "Interaction term: not estimable"
}

p_6e <- analysis_data %>%
  ggplot(aes(x = total_cafi, y = condition_score, color = size_class)) +
  geom_point(alpha = 0.5, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_x_sqrt(breaks = c(0, 25, 100, 225)) +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  labs(
    title = "E. Size-Dependent CAFI Effects on Condition",
    subtitle = interaction_subtitle,
    x = "Total CAFI Abundance (sqrt scale)",
    y = "Condition Score",
    color = "Coral Size\nClass"
  ) +
  theme_publication()

save_pub_fig(p_6e, "fig6e_size_dependent_cafi_effects.png", dir = FIG_DIR)

# ============================================================================
# ANALYSIS 6: COMPOSITION EFFECTS (PCA OF CAFI COMMUNITY)
# ============================================================================

cat("\n6. Testing CAFI composition effects...\n")

# Get corals with community data
corals_with_comm <- intersect(rownames(community_matrix), analysis_data$coral_id)
comm_subset <- community_matrix[corals_with_comm, ]

# Remove empty columns
comm_subset <- comm_subset[, colSums(comm_subset) > 0]

if (nrow(comm_subset) > 20 && ncol(comm_subset) > 3) {

  # Hellinger transformation
  comm_hell <- decostand(comm_subset, method = "hellinger")

  # PCA on community
  pca_comm <- rda(comm_hell)

  # Extract PC scores
  pc_scores <- as.data.frame(scores(pca_comm, display = "sites", choices = 1:2)) %>%
    rownames_to_column("coral_id")

  # Merge with condition data
  composition_data <- analysis_data %>%
    filter(coral_id %in% corals_with_comm) %>%
    left_join(pc_scores, by = "coral_id")

  # Model: condition ~ PC1 + PC2 + size + site
  model_composition <- lm(condition_score ~ PC1 + PC2 + log_volume + site,
                          data = composition_data)

  coef_pc1 <- extract_model_stats(model_composition, "PC1")
  coef_pc2 <- extract_model_stats(model_composition, "PC2")

  cat(sprintf("  Composition PC1 effect: %.3f (p %s)\n",
              coef_pc1$estimate, format_p(coef_pc1$p.value)))
  cat(sprintf("  Composition PC2 effect: %.3f (p %s)\n",
              coef_pc2$estimate, format_p(coef_pc2$p.value)))

  # Figure 6F: Composition PC1 vs Condition
  show_line_6f <- coef_pc1$p.value < 0.10

  p_6f <- composition_data %>%
    ggplot(aes(x = PC1, y = condition_score)) +
    geom_point(aes(color = site), alpha = 0.6, size = 3) +
    {if (show_line_6f)
      geom_smooth(method = "lm", formula = y ~ x,
                  color = "black", linewidth = 1.2, se = TRUE)
    } +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    scale_color_site() +
    labs(
      title = "F. Coral Condition vs CAFI Composition",
      subtitle = sprintf("PC1 effect = %.3f, p %s | Controlling for size + site",
                         coef_pc1$estimate, format_p(coef_pc1$p.value)),
      x = "CAFI Composition (PC1)",
      y = "Condition Score",
      color = "Site"
    ) +
    theme_publication()

  save_pub_fig(p_6f, "fig6f_condition_vs_composition.png", dir = FIG_DIR)

} else {
  cat("  Insufficient data for composition PCA\n")
  coef_pc1 <- list(estimate = NA, p.value = NA)
  coef_pc2 <- list(estimate = NA, p.value = NA)
}

# ============================================================================
# COMBINED FIGURE 6 PANEL
# ============================================================================

cat("\n7. Creating combined Figure 6 panel...\n")

# Use plot_spacer() for any missing panels
p_6d_final <- if (exists("p_6d")) p_6d else plot_spacer()
p_6f_final <- if (exists("p_6f")) p_6f else plot_spacer()

fig6_combined <- (p_6a + p_6b) /
  (p_6c + p_6d_final) /
  (p_6e + p_6f_final) +
  plot_annotation(
    title = "Figure 6. CAFI Communities Influence Coral Condition",
    subtitle = "Effects of CAFI abundance, diversity, and functional groups on coral physiology",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_pub_fig(fig6_combined, "fig6_cafi_feedbacks_combined.png",
             dir = FIG_DIR, width = 14, height = 16)

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n8. Saving results...\n")

# Compile all results
fig6_results <- bind_rows(
  tibble(predictor = "Total CAFI Abundance", estimate = coef_abund$estimate,
         se = coef_abund$se, p_value = coef_abund$p.value,
         significant = coef_abund$significant),
  tibble(predictor = "OTU Richness", estimate = coef_rich$estimate,
         se = coef_rich$se, p_value = coef_rich$p.value,
         significant = coef_rich$significant),
  tibble(predictor = "Trapezia Presence", estimate = coef_trap_pres$estimate,
         se = coef_trap_pres$se, p_value = coef_trap_pres$p.value,
         significant = coef_trap_pres$significant),
  tibble(predictor = "CAFI x Size Interaction", estimate = coef_interaction$estimate,
         se = coef_interaction$se, p_value = coef_interaction$p.value,
         significant = coef_interaction$significant),
  tibble(predictor = "Composition PC1", estimate = coef_pc1$estimate,
         se = NA, p_value = coef_pc1$p.value,
         significant = if (!is.na(coef_pc1$p.value)) coef_pc1$p.value < 0.05 else NA)
)

save_table(fig6_results, "fig6_cafi_condition_results.csv")
save_table(functional_df, "fig6_functional_group_effects.csv")

save_object(list(
  model_abundance = model_abund,
  model_richness = model_rich,
  model_functional = model_functional,
  model_interaction = model_interaction,
  model_composition = if (exists("model_composition")) model_composition else NULL
), "fig6_models.rds")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("Figure 6 Analysis Summary\n")
cat("========================================\n\n")

cat("CAFI EFFECTS ON CORAL CONDITION:\n\n")

cat("1. TOTAL ABUNDANCE:\n")
if (coef_abund$significant) {
  cat(sprintf("   -> More CAFI = %s condition (effect = %.4f)\n",
              ifelse(coef_abund$estimate > 0, "BETTER", "WORSE"),
              coef_abund$estimate))
} else {
  cat("   -> No significant effect of total CAFI abundance\n")
}

cat("\n2. SPECIES RICHNESS:\n")
if (coef_rich$significant) {
  cat(sprintf("   -> More OTUs = %s condition (effect = %.4f)\n",
              ifelse(coef_rich$estimate > 0, "BETTER", "WORSE"),
              coef_rich$estimate))
} else {
  cat("   -> No significant effect of CAFI richness\n")
}

cat("\n3. KEY FUNCTIONAL GROUPS:\n")
sig_groups <- functional_df %>% filter(significant)
if (nrow(sig_groups) > 0) {
  for (i in 1:nrow(sig_groups)) {
    direction <- ifelse(sig_groups$estimate[i] > 0, "improves", "reduces")
    cat(sprintf("   -> %s %s condition\n",
                sig_groups$group[i], direction))
  }
} else {
  cat("   -> No individual functional group has significant effect\n")
}

cat("\n4. TRAPEZIA (KEY MUTUALIST):\n")
if (coef_trap_pres$significant) {
  cat(sprintf("   -> Presence associated with %.3f condition change\n",
              coef_trap_pres$estimate))
  cat("   -> SUPPORTS mutualist hypothesis!\n")
} else {
  cat("   -> No significant Trapezia effect detected\n")
}

cat("\n5. SIZE-DEPENDENT EFFECTS:\n")
if (coef_interaction$significant) {
  cat("   -> CAFI effects VARY with coral size (significant interaction)\n")
} else {
  cat("   -> CAFI effects are consistent across coral sizes\n")
}

cat("\nFigures saved to:", FIG_DIR, "\n")
cat("Tables saved to:", TABLES_DIR, "\n\n")
cat("Figure 6 analysis complete!\n\n")
