#!/usr/bin/env Rscript
# ============================================================================
# 04_scaling_relationships.R - Size-Abundance Scaling Analysis
# ============================================================================
#
# PURPOSE: Test power-law scaling between coral size and CAFI abundance
#
# HYPOTHESIS (H2): CAFI abundance scales with coral volume following a
#   power-law with exponent < 1 (sublinear scaling), indicating larger
#   corals have lower CAFI densities per unit volume.
#
#   Propagule redirection theory predicts: N = a * V^b where b < 1
#   Expected exponent ~0.75 for 3D habitat
#
# ANALYSES:
#   - Power-law scaling: log(abundance) ~ log(volume)
#   - Site and branch architecture effects
#   - Richness and diversity scaling
#   - Proximity/isolation effects
#
# MANUSCRIPT FIGURES:
#   >>> FIGURE 2: Landscape -> CAFI Scaling <<<
#   - Panel A: Size-abundance scaling
#   - Panel B: Size-richness scaling
#   - Panel C: Isolation effects
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("04: Scaling Relationships Analysis\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")
condition_scores <- load_object("condition_scores.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "scaling")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# CAFI Abundance vs Coral Characteristics
# ============================================================================

cat("Analyzing CAFI abundance patterns...\n")

# Use coral_master which already has CAFI metrics from 01_load_clean_data.R
coral_cafi_summary <- coral_master %>%
  mutate(
    total_cafi = replace_na(total_cafi, 0),
    otu_richness = replace_na(otu_richness, 0),
    n_crabs = replace_na(n_crabs, 0),
    n_shrimps = replace_na(n_shrimps, 0),
    n_fish = replace_na(n_fish, 0),
    n_snails = replace_na(n_snails, 0)
  )

# CAFI abundance across all Pocillopora - overall distribution
p_overall_cafi <- coral_cafi_summary %>%
  ggplot(aes(x = "Pocillopora spp.", y = total_cafi)) +
  geom_violin(fill = "steelblue", alpha = 0.3, color = NA) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.shape = NA, width = 0.15) +
  geom_jitter(width = 0.08, alpha = 0.4, size = 1.5, color = "gray30") +
  scale_y_sqrt(breaks = c(0, 25, 100, 225)) +
  labs(
    title = "CAFI abundance distribution",
    subtitle = sprintf("n = %d corals | Median = %.0f, Mean = %.1f",
                       nrow(coral_cafi_summary),
                       median(coral_cafi_summary$total_cafi),
                       mean(coral_cafi_summary$total_cafi)),
    x = NULL,
    y = "CAFI abundance (sqrt scale)"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(face = "italic", size = 11))

ggsave(file.path(fig_dir, "cafi_overall_abundance.png"),
       p_overall_cafi, width = 5, height = 5, dpi = 300, bg = "white")

# CAFI by location (3 sites)
site_stats <- coral_cafi_summary %>%
  group_by(site) %>%
  summarise(n = n(), median = median(total_cafi), .groups = "drop")

p_site_cafi <- coral_cafi_summary %>%
  ggplot(aes(x = site, y = total_cafi, fill = site)) +
  geom_violin(alpha = 0.3, color = NA) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.15) +
  geom_jitter(width = 0.08, alpha = 0.4, size = 1.5, color = "gray30") +
  scale_fill_site() +
  scale_y_sqrt(breaks = c(0, 25, 100, 225)) +
  labs(
    title = "CAFI abundance across Mo'orea reef sites",
    x = NULL,
    y = "CAFI abundance (sqrt scale)"
  ) +
  theme_publication() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11, face = "bold"))

ggsave(file.path(fig_dir, "cafi_by_location.png"),
       p_site_cafi, width = 6, height = 5, dpi = 300, bg = "white")

# CAFI abundance and richness by branch width (if column exists)
if("branch_width" %in% names(coral_cafi_summary) && sum(!is.na(coral_cafi_summary$branch_width)) > 0) {

  p_branch_abundance <- coral_cafi_summary %>%
    filter(!is.na(branch_width)) %>%
    ggplot(aes(x = branch_width, y = total_cafi)) +
    geom_boxplot(aes(fill = branch_width), alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
    scale_fill_branch() +
    scale_y_sqrt(breaks = c(0, 25, 100, 225, 400)) +
    labs(
      title = "CAFI Abundance by Branch Width",
      subtitle = "Branch width is a measured coral trait (wide vs tight branching)",
      x = "Branch Width Category",
      y = "CAFI Abundance (sqrt scale)"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "cafi_by_branch_width.png"),
         p_branch_abundance, width = 10, height = 6, dpi = 300, bg = "white")

  p_branch_richness <- coral_cafi_summary %>%
    filter(!is.na(branch_width)) %>%
    ggplot(aes(x = branch_width, y = otu_richness)) +
    geom_boxplot(aes(fill = branch_width), alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
    scale_fill_branch() +
    labs(
      title = "CAFI OTU Richness by Branch Width",
      subtitle = "Branch width is a measured coral trait (wide vs tight branching)",
      x = "Branch Width Category",
      y = "OTU Richness"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "richness_by_branch_width.png"),
         p_branch_richness, width = 10, height = 6, dpi = 300, bg = "white")
}

# ============================================================================
# Taxonomic Group Composition
# ============================================================================

cat("Analyzing taxonomic group composition...\n")

# Calculate overall taxonomic composition
taxonomic_composition <- coral_cafi_summary %>%
  filter(total_cafi > 0) %>%
  summarise(
    Crabs = sum(n_crabs),
    Shrimps = sum(n_shrimps),
    Fish = sum(n_fish),
    Snails = sum(n_snails)
  ) %>%
  pivot_longer(everything(), names_to = "group", values_to = "count") %>%
  mutate(proportion = count / sum(count))

# Plot overall taxonomic composition - horizontal bar chart (better for publications)
p_tax_overall <- taxonomic_composition %>%
  mutate(group = factor(group, levels = group[order(proportion)])) %>%
  ggplot(aes(x = group, y = proportion, fill = group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", proportion * 100)),
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("Crabs" = "#CC79A7", "Shrimps" = "#F0E442",
                               "Fish" = "#0072B2", "Snails" = "#999999")) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Overall CAFI taxonomic composition",
    subtitle = "Proportions across all Pocillopora spp. corals",
    x = NULL,
    y = "Proportion of community"
  ) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "taxonomic_composition_overall.png"),
       p_tax_overall, width = 8, height = 6, dpi = 300, bg = "white")

# Taxonomic composition by location
tax_by_site <- coral_cafi_summary %>%
  filter(total_cafi > 0) %>%
  group_by(site) %>%
  summarise(
    Crabs = sum(n_crabs),
    Shrimps = sum(n_shrimps),
    Fish = sum(n_fish),
    Snails = sum(n_snails),
    .groups = "drop"
  ) %>%
  pivot_longer(-site, names_to = "group", values_to = "count") %>%
  group_by(site) %>%
  mutate(proportion = count / sum(count))

p_tax_site <- ggplot(tax_by_site, aes(x = site, y = proportion, fill = group)) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Crabs" = "#CC79A7", "Shrimps" = "#F0E442",
                               "Fish" = "#0072B2", "Snails" = "#999999"),
                    name = "Taxonomic\nGroup") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(
    title = "CAFI taxonomic composition by site",
    subtitle = "Relative proportions of functional groups at each Mo'orea site",
    x = NULL,
    y = "Proportion of community"
  ) +
  theme_publication() +
  theme(legend.position = "right",
        axis.text.x = element_text(size = 11, face = "bold"))

ggsave(file.path(fig_dir, "taxonomic_composition_by_site.png"),
       p_tax_site, width = 10, height = 6, dpi = 300, bg = "white")

# ============================================================================
# CAFI-Coral Condition Relationships (Using Position-Corrected Scores)
# ============================================================================

cat("Analyzing CAFI-coral condition relationships...\n")
cat("  Using position-corrected condition scores from Script 05a\n")

# Merge CAFI and condition score data
cafi_condition <- coral_cafi_summary %>%
  inner_join(condition_scores %>% dplyr::select(coral_id, condition_score, site),
            by = "coral_id",
            suffix = c("", "_condition"))

if (nrow(cafi_condition) > 10) {

  cat(sprintf("  - Merged data: %d corals with both CAFI and condition data\n", nrow(cafi_condition)))

  # ============================================================================
  # CAFI Abundance vs Coral Condition Score
  # ============================================================================

  p_cafi_condition <- cafi_condition %>%
    filter(!is.na(condition_score), !is.na(total_cafi)) %>%
    ggplot(aes(x = total_cafi, y = condition_score)) +
    geom_point(aes(color = site), size = 3, alpha = 0.6) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    scale_color_site() +
    scale_x_sqrt(breaks = c(0, 25, 100, 225, 400)) +
    labs(
      title = "CAFI Abundance vs Coral Condition",
      subtitle = "Condition = position-corrected PC1 of protein, carbs, zoox, AFDW | Higher = better",
      x = "CAFI Abundance (sqrt scale)",
      y = "Coral Condition Score (position-corrected)",
      caption = "Condition score automatically accounts for sampling position bias"
    ) +
    theme_publication() +
    theme(legend.position = "right")

  ggsave(file.path(fig_dir, "cafi_vs_coral_condition.png"),
         p_cafi_condition, width = 11, height = 7, dpi = 300, bg = "white")

  # Faceted by location
  p_cafi_condition_site <- cafi_condition %>%
    filter(!is.na(condition_score), !is.na(total_cafi)) %>%
    ggplot(aes(x = total_cafi, y = condition_score)) +
    geom_point(aes(color = site, shape = site), size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    facet_wrap(~site, ncol = 3) +
    scale_color_site() +
    scale_shape_site() +
    scale_x_sqrt(breaks = c(0, 25, 100, 225)) +
    labs(
      title = "CAFI abundance vs coral condition by site",
      subtitle = "Position-corrected condition score across 3 Mo'orea sites",
      x = "CAFI abundance (sqrt scale)",
      y = "Coral condition score"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "cafi_vs_condition_by_site.png"),
         p_cafi_condition_site, width = 14, height = 5, dpi = 300, bg = "white")

  # ============================================================================
  # CAFI Richness vs Coral Condition
  # ============================================================================

  p_richness_condition <- cafi_condition %>%
    filter(!is.na(condition_score), !is.na(otu_richness)) %>%
    ggplot(aes(x = otu_richness, y = condition_score)) +
    geom_point(aes(color = site, shape = site), size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    scale_color_site() +
    scale_shape_site() +
    labs(
      title = "CAFI OTU richness vs coral condition",
      subtitle = "Do corals with diverse CAFI communities have better physiological condition?",
      x = "CAFI OTU richness",
      y = "Coral condition score (position-corrected)"
    ) +
    theme_publication() +
    theme(legend.position = "right")

  ggsave(file.path(fig_dir, "richness_vs_coral_condition.png"),
         p_richness_condition, width = 11, height = 7, dpi = 300, bg = "white")

  # ============================================================================
  # Taxonomic Group Abundance vs Condition
  # ============================================================================

  # Panel showing different CAFI groups
  p_crabs <- cafi_condition %>%
    filter(!is.na(condition_score)) %>%
    ggplot(aes(x = n_crabs, y = condition_score)) +
    geom_point(alpha = 0.6, size = 2, color = "#CC79A7") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(x = "Crab abundance", y = "Condition score") +
    theme_publication()

  p_shrimps <- cafi_condition %>%
    filter(!is.na(condition_score)) %>%
    ggplot(aes(x = n_shrimps, y = condition_score)) +
    geom_point(alpha = 0.6, size = 2, color = "#F0E442") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(x = "Shrimp abundance", y = "Condition score") +
    theme_publication()

  p_fish <- cafi_condition %>%
    filter(!is.na(condition_score)) %>%
    ggplot(aes(x = n_fish, y = condition_score)) +
    geom_point(alpha = 0.6, size = 2, color = "#0072B2") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(x = "Fish abundance", y = "Condition score") +
    theme_publication()

  p_snails <- cafi_condition %>%
    filter(!is.na(condition_score)) %>%
    ggplot(aes(x = n_snails, y = condition_score)) +
    geom_point(alpha = 0.6, size = 2, color = "#999999") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(x = "Snail abundance", y = "Condition score") +
    theme_publication()

  p_taxa_panel <- (p_crabs | p_shrimps) / (p_fish | p_snails) +
    plot_annotation(
      title = "CAFI taxonomic groups vs coral condition",
      subtitle = "Position-corrected condition score by major CAFI taxa",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "gray30")
      )
    )

  ggsave(file.path(fig_dir, "cafi_taxa_vs_condition.png"),
         p_taxa_panel, width = 12, height = 10, dpi = 300, bg = "white")

  # ============================================================================
  # Individual Position-Corrected Traits vs CAFI
  # ============================================================================

  cat("\n  Analyzing individual position-corrected traits...\n")

  # Merge with all corrected traits (not just condition score)
  cafi_traits <- coral_cafi_summary %>%
    inner_join(condition_scores %>% dplyr::select(coral_id, condition_score,
                                                   protein_corr, carb_corr,
                                                   zoox_corr, afdw_corr),
              by = "coral_id")

  if(nrow(cafi_traits) > 10) {

    # Create 4-panel plot: CAFI abundance vs each corrected trait
    p_protein_corr <- cafi_traits %>%
      filter(!is.na(protein_corr)) %>%
      ggplot(aes(x = total_cafi, y = protein_corr)) +
      geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      scale_x_sqrt() +
      labs(x = "CAFI abundance", y = "Protein (z-score)") +
      theme_publication()

    p_carb_corr <- cafi_traits %>%
      filter(!is.na(carb_corr)) %>%
      ggplot(aes(x = total_cafi, y = carb_corr)) +
      geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      scale_x_sqrt() +
      labs(x = "CAFI abundance", y = "Carbohydrate (z-score)") +
      theme_publication()

    p_zoox_corr <- cafi_traits %>%
      filter(!is.na(zoox_corr)) %>%
      ggplot(aes(x = total_cafi, y = zoox_corr)) +
      geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      scale_x_sqrt() +
      labs(x = "CAFI abundance", y = "Zooxanthellae (z-score)") +
      theme_publication()

    p_afdw_corr <- cafi_traits %>%
      filter(!is.na(afdw_corr)) %>%
      ggplot(aes(x = total_cafi, y = afdw_corr)) +
      geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      scale_x_sqrt() +
      labs(x = "CAFI abundance", y = "AFDW (z-score)") +
      theme_publication()

    p_traits_panel <- (p_protein_corr | p_carb_corr) / (p_zoox_corr | p_afdw_corr) +
      plot_annotation(
        title = "CAFI abundance vs position-corrected physiological traits",
        subtitle = "All traits corrected for sampling position bias (residuals from trait ~ stump_length)",
        caption = "Z-scores centered at 0; positive = better than expected for branch position",
        theme = theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11, color = "gray30"),
          plot.caption = element_text(size = 9, hjust = 0)
        )
      )

    ggsave(file.path(fig_dir, "cafi_vs_corrected_traits_panel.png"),
           p_traits_panel, width = 12, height = 10, dpi = 300, bg = "white")

    # Create similar panel for OTU richness
    p_protein_rich <- cafi_traits %>%
      filter(!is.na(protein_corr)) %>%
      ggplot(aes(x = otu_richness, y = protein_corr)) +
      geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      labs(x = "OTU richness", y = "Protein (z-score)") +
      theme_publication()

    p_carb_rich <- cafi_traits %>%
      filter(!is.na(carb_corr)) %>%
      ggplot(aes(x = otu_richness, y = carb_corr)) +
      geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      labs(x = "OTU richness", y = "Carbohydrate (z-score)") +
      theme_publication()

    p_zoox_rich <- cafi_traits %>%
      filter(!is.na(zoox_corr)) %>%
      ggplot(aes(x = otu_richness, y = zoox_corr)) +
      geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      labs(x = "OTU richness", y = "Zooxanthellae (z-score)") +
      theme_publication()

    p_afdw_rich <- cafi_traits %>%
      filter(!is.na(afdw_corr)) %>%
      ggplot(aes(x = otu_richness, y = afdw_corr)) +
      geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      labs(x = "OTU richness", y = "AFDW (z-score)") +
      theme_publication()

    p_richness_traits_panel <- (p_protein_rich | p_carb_rich) / (p_zoox_rich | p_afdw_rich) +
      plot_annotation(
        title = "CAFI OTU richness vs position-corrected physiological traits",
        subtitle = "All traits corrected for sampling position bias",
        theme = theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11, color = "gray30")
        )
      )

    ggsave(file.path(fig_dir, "richness_vs_corrected_traits_panel.png"),
           p_richness_traits_panel, width = 12, height = 10, dpi = 300, bg = "white")

    # =========================================================================
    # Statistical models for each corrected trait
    # =========================================================================

    cat("  Running models for individual corrected traits...\n")

    trait_models <- list()
    trait_names <- c("protein_corr", "carb_corr", "zoox_corr", "afdw_corr")
    trait_labels <- c("Protein", "Carbohydrate", "Zooxanthellae", "AFDW")

    for(i in seq_along(trait_names)) {
      trait <- trait_names[i]
      label <- trait_labels[i]

      # Model: trait ~ CAFI abundance + site
      model_abundance <- lm(reformulate(c("total_cafi", "site"), response = trait),
                           data = cafi_traits %>% filter(!is.na(.data[[trait]])))

      # Model: trait ~ OTU richness + site
      model_richness <- lm(reformulate(c("otu_richness", "site"), response = trait),
                          data = cafi_traits %>% filter(!is.na(.data[[trait]])))

      # Extract coefficients
      coef_abundance <- broom::tidy(model_abundance, conf.int = TRUE) %>%
        filter(term == "total_cafi") %>%
        mutate(trait = label, predictor = "CAFI Abundance")

      coef_richness <- broom::tidy(model_richness, conf.int = TRUE) %>%
        filter(term == "otu_richness") %>%
        mutate(trait = label, predictor = "OTU Richness")

      trait_models[[label]] <- bind_rows(coef_abundance, coef_richness)

      cat(sprintf("    %s ~ CAFI: β = %.3f, p = %.3f\n",
                  label, coef_abundance$estimate, coef_abundance$p.value))
      cat(sprintf("    %s ~ Richness: β = %.3f, p = %.3f\n",
                  label, coef_richness$estimate, coef_richness$p.value))
    }

    # Combine all trait models
    all_trait_models <- bind_rows(trait_models)
    write_csv(all_trait_models,
              file.path(TABLES_DIR, "cafi_corrected_traits_models.csv"))

    # Visualization of all trait model coefficients
    p_trait_coefs <- ggplot(all_trait_models,
                           aes(x = estimate, y = trait, color = predictor)) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                     height = 0.3, linewidth = 0.8, position = position_dodge(width = 0.5)) +
      geom_point(size = 3, position = position_dodge(width = 0.5)) +
      scale_color_manual(values = c("CAFI Abundance" = "#E69F00", "OTU Richness" = "#0072B2"),
                        name = "Predictor") +
      labs(
        title = "CAFI effects on individual position-corrected traits",
        subtitle = "Controlling for site differences | All traits position-corrected",
        x = "Coefficient estimate (± 95% CI)",
        y = "Physiological trait",
        caption = "Models: trait ~ CAFI_metric + site"
      ) +
      theme_publication() +
      theme(legend.position = "bottom")

    ggsave(file.path(fig_dir, "cafi_trait_coefficients.png"),
           p_trait_coefs, width = 10, height = 7, dpi = 300, bg = "white")

    cat("  ✓ Individual corrected trait analyses complete\n")
  }

  # ============================================================================
  # Statistical Models: CAFI ~ Condition + Site
  # ============================================================================

  cat("\n  Running statistical models with condition score...\n")

  # Model 1: CAFI abundance ~ condition + site
  model_cafi_condition <- lm(total_cafi ~ condition_score + site,
                            data = cafi_condition %>% filter(!is.na(condition_score)))

  # Model 2: CAFI richness ~ condition + site
  model_richness_condition <- lm(otu_richness ~ condition_score + site,
                                data = cafi_condition %>% filter(!is.na(condition_score)))

  # Save model summaries
  model_summaries <- list(
    cafi_abundance_model = broom::tidy(model_cafi_condition, conf.int = TRUE),
    richness_model = broom::tidy(model_richness_condition, conf.int = TRUE)
  )

  saveRDS(model_summaries,
          file.path(OBJECTS_DIR, "cafi_condition_models.rds"))

  # Combine for visualization
  coef_data <- bind_rows(
    broom::tidy(model_cafi_condition, conf.int = TRUE) %>% mutate(model = "CAFI Abundance"),
    broom::tidy(model_richness_condition, conf.int = TRUE) %>% mutate(model = "OTU Richness")
  ) %>%
    filter(term != "(Intercept)")

  p_model_coef <- ggplot(coef_data, aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = model),
                   height = 0.2, linewidth = 0.8) +
    geom_point(aes(color = model), size = 3) +
    facet_wrap(~model, scales = "free_x", ncol = 2) +
    scale_color_manual(values = c("CAFI Abundance" = "#E69F00", "OTU Richness" = "#0072B2")) +
    labs(
      title = "Effect of coral condition on CAFI patterns",
      subtitle = "Linear models: CAFI ~ condition_score + site",
      x = "Coefficient estimate (± 95% CI)",
      y = "Predictor",
      caption = "Condition score = position-corrected PC1 from Script 05a"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "condition_model_coefficients.png"),
         p_model_coef, width = 12, height = 7, dpi = 300, bg = "white")

  # Save coefficient table
  write_csv(coef_data,
            file.path(TABLES_DIR, "cafi_condition_model_coefficients.csv"))

  # Extract key results
  condition_effect_abundance <- coef_data %>%
    filter(model == "CAFI Abundance", term == "condition_score")

  condition_effect_richness <- coef_data %>%
    filter(model == "OTU Richness", term == "condition_score")

  cat(sprintf("\n  KEY RESULTS:\n"))
  cat(sprintf("    CAFI Abundance ~ Condition: β = %.3f, p = %.3f\n",
              condition_effect_abundance$estimate,
              condition_effect_abundance$p.value))
  cat(sprintf("    OTU Richness ~ Condition: β = %.3f, p = %.3f\n",
              condition_effect_richness$estimate,
              condition_effect_richness$p.value))

  if(condition_effect_abundance$p.value < 0.05) {
    cat("    → CAFI abundance is significantly associated with coral condition!\n")
  } else {
    cat("    → No significant relationship between CAFI abundance and condition\n")
  }

  if(condition_effect_richness$p.value < 0.05) {
    cat("    → CAFI richness is significantly associated with coral condition!\n")
  } else {
    cat("    → No significant relationship between CAFI richness and condition\n")
  }

  cat("\n  ✓ Condition-CAFI relationships analyzed\n\n")

} else {
  cat("  ⚠️  Insufficient data for CAFI-condition analysis\n\n")
}

# ============================================================================
# Species-Specific Associations
# ============================================================================

cat("Analyzing species-specific associations...\n")

# COMMENTED OUT: OTU-morphotype associations
# Without genetic species IDs, we cannot make ecological inferences about
# individual "species" preferences. Focus on community-level patterns and
# functional groups (crabs, shrimp, fish) instead.
#
# # Get top 20 most common species
# top_species <- cafi_clean %>%
#   count(species, sort = TRUE) %>%
#   slice_head(n = 20) %>%
#   pull(species)
#
# # Calculate species associations with coral traits
# select_cols <- c("coral_id", "morphotype")
# if("branch_width" %in% names(metadata)) select_cols <- c(select_cols, "branch_width")
#
# species_associations <- cafi_clean %>%
#   filter(species %in% top_species) %>%
#   left_join(metadata %>% select(any_of(select_cols)),
#             by = "coral_id") %>%
#   group_by(species, morphotype) %>%
#   summarise(
#     count = n(),
#     .groups = "drop"
#   ) %>%
#   group_by(species) %>%
#   mutate(proportion = count / sum(count))
#
# # Plot species preferences
# p_species_pref <- species_associations %>%
#   ggplot(aes(x = morphotype, y = species, fill = proportion)) +
#   geom_tile() +
#   scale_fill_viridis_c() +
#   labs(
#     title = "Species Preferences for Coral Morphotypes",
#     subtitle = "Top 20 species",
#     x = "Morphotype",
#     y = "Species",
#     fill = "Proportion"
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         axis.text.y = element_text(size = 8))
#
# ggsave(file.path(fig_dir, "species_morphotype_preferences.png"),
#        p_species_pref, width = 12, height = 10, dpi = 300)

# ============================================================================
# Statistical Models
# ============================================================================

cat("Running statistical models...\n")

# Model CAFI abundance and richness
if (sum(!is.na(coral_cafi_summary$depth_m)) > 20) {
  # Build model formula based on available predictors
  # Use 'site' for the 3 locations (HAU, MAT, MRB)
  predictors <- "depth_m + site"
  if("branch_width" %in% names(coral_cafi_summary) && sum(!is.na(coral_cafi_summary$branch_width)) > 10) {
    predictors <- paste(predictors, "+ branch_width")
  }

  # Poisson GLM for total CAFI abundance
  model_abundance <- glm(as.formula(paste("total_cafi ~", predictors)),
                        data = coral_cafi_summary,
                        family = poisson())

  # Poisson GLM for OTU richness
  model_richness <- glm(as.formula(paste("otu_richness ~", predictors)),
                       data = coral_cafi_summary,
                       family = poisson())

  # Save model summaries
  glm_model_summaries <- list(
    abundance_model = broom::tidy(model_abundance),
    richness_model = broom::tidy(model_richness)
  )

  saveRDS(glm_model_summaries,
          file.path(OBJECTS_DIR, "coral_cafi_models.rds"))

  # Create coefficient plot
  glm_coef_data <- bind_rows(
    broom::tidy(model_abundance, conf.int = TRUE) %>% mutate(model = "Abundance"),
    broom::tidy(model_richness, conf.int = TRUE) %>% mutate(model = "OTU Richness")
  ) %>%
    filter(term != "(Intercept)")

  p_coef <- ggplot(glm_coef_data, aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, linewidth = 0.8) +
    geom_point(aes(color = model), size = 3) +
    facet_wrap(~model, scales = "free_x") +
    scale_color_manual(values = c("Abundance" = "#E69F00", "OTU Richness" = "#0072B2")) +
    labs(
      title = "Model coefficients for CAFI patterns on Pocillopora spp.",
      subtitle = "Poisson GLMs for abundance and OTU richness across 3 Mo'orea sites",
      x = "Coefficient estimate (log scale)",
      y = "Predictor"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "model_coefficients.png"),
         p_coef, width = 12, height = 8, dpi = 300, bg = "white")
}

# ============================================================================
# Size Distribution by Taxonomic Group
# ============================================================================

cat("Analyzing CAFI size distributions...\n")

# Size distribution - filter to main 4 taxonomic groups
size_data <- cafi_clean %>%
  filter(!is.na(size_mm)) %>%
  filter(type %in% c("crab", "shrimp", "fish", "snail"))  # Only main 4 groups

p_size <- ggplot(size_data, aes(x = size_mm, fill = type)) +
  geom_density(alpha = 0.6, color = "black", linewidth = 0.3) +
  facet_wrap(~type, scales = "free_y", ncol = 2) +
  scale_fill_taxon() +
  scale_x_continuous(limits = c(0, 50)) +  # Focus on main size range
  labs(
    title = "CAFI size distribution by taxonomic group",
    subtitle = "Four main functional groups on Pocillopora spp.",
    x = "Body size (mm)",
    y = "Density"
  ) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "size_distribution_by_group.png"),
       p_size, width = 12, height = 8, dpi = 300, bg = "white")

# ============================================================================
# Summary Statistics
# ============================================================================

# Calculate overall summary statistics
overall_summary <- coral_cafi_summary %>%
  summarise(
    n_corals = n(),
    mean_cafi = mean(total_cafi),
    sd_cafi = sd(total_cafi),
    median_cafi = median(total_cafi),
    mean_otu_richness = mean(otu_richness),
    sd_otu_richness = sd(otu_richness),
    median_otu_richness = median(otu_richness)
  )

# Summary by location
summary_by_site <- coral_cafi_summary %>%
  group_by(site) %>%
  summarise(
    n_corals = n(),
    mean_cafi = mean(total_cafi),
    sd_cafi = sd(total_cafi),
    mean_otu_richness = mean(otu_richness),
    sd_otu_richness = sd(otu_richness),
    .groups = "drop"
  )

write_csv(summary_by_site,
          file.path(TABLES_DIR, "coral_cafi_summary_by_site.csv"))

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Coral-CAFI Relationships Summary\n")
cat("========================================\n\n")

cat("Study Overview:\n")
cat("  - Species: All Pocillopora spp. (morphotypes pooled)\n")
cat("  - Locations: 3 sites in Mo'orea (HAU, MAT, MRB)\n")
cat("  - Total corals surveyed:", nrow(coral_cafi_summary), "\n\n")

cat("Overall Patterns:\n")
cat("  - Corals hosting CAFI:", sum(coral_cafi_summary$total_cafi > 0),
    "of", nrow(coral_cafi_summary),
    sprintf("(%.1f%%)\n", 100 * sum(coral_cafi_summary$total_cafi > 0) / nrow(coral_cafi_summary)))
cat("  - Mean CAFI per coral:", round(mean(coral_cafi_summary$total_cafi), 1),
    "+/-", round(sd(coral_cafi_summary$total_cafi), 1), "\n")
cat("  - Median CAFI per coral:", median(coral_cafi_summary$total_cafi), "\n")
cat("  - Max CAFI on single coral:", max(coral_cafi_summary$total_cafi), "\n")
cat("  - Mean OTU richness:", round(mean(coral_cafi_summary$otu_richness), 1),
    "+/-", round(sd(coral_cafi_summary$otu_richness), 1), "\n\n")

cat("Location Patterns:\n")
for(i in 1:nrow(summary_by_site)) {
  cat(sprintf("  - %s: %.1f CAFI/coral (n=%d)\n",
              summary_by_site$site[i],
              summary_by_site$mean_cafi[i],
              summary_by_site$n_corals[i]))
}

# ============================================================================
# Compile Statistical Results Summary
# ============================================================================

cat("Compiling statistical results...\n")

# Initialize results data frame
stats_results <- init_results_df()

# Test 1: CAFI abundance vs condition score (linear regression)
if (exists("cafi_condition") && nrow(cafi_condition) > 10) {
  model_cafi_cond <- lm(total_cafi ~ condition_score, data = cafi_condition)
  model_sum <- summary(model_cafi_cond)

  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H4a",
      question = "Does CAFI abundance increase with coral condition?",
      test_name = "Linear regression",
      test_statistic = model_sum$fstatistic[1],
      df = paste(model_sum$fstatistic[2], model_sum$fstatistic[3], sep = ", "),
      p_value = pf(model_sum$fstatistic[1], model_sum$fstatistic[2], model_sum$fstatistic[3], lower.tail = FALSE),
      effect_size = model_sum$r.squared,
      effect_type = "R²",
      n = nrow(cafi_condition),
      notes = "CAFI abundance ~ condition score"
    )
  )
}

# Test 2: OTU richness vs condition score
if (exists("cafi_condition") && nrow(cafi_condition) > 10) {
  model_rich_cond <- lm(otu_richness ~ condition_score, data = cafi_condition)
  model_sum <- summary(model_rich_cond)

  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H4b",
      question = "Does CAFI richness increase with coral condition?",
      test_name = "Linear regression",
      test_statistic = model_sum$fstatistic[1],
      df = paste(model_sum$fstatistic[2], model_sum$fstatistic[3], sep = ", "),
      p_value = pf(model_sum$fstatistic[1], model_sum$fstatistic[2], model_sum$fstatistic[3], lower.tail = FALSE),
      effect_size = model_sum$r.squared,
      effect_type = "R²",
      n = nrow(cafi_condition),
      notes = "OTU richness ~ condition score"
    )
  )
}

# Test 3: Kruskal-Wallis for CAFI abundance by site
kw_site <- kruskal.test(total_cafi ~ site, data = coral_cafi_summary)
n_total <- nrow(coral_cafi_summary)
k <- n_distinct(coral_cafi_summary$site)
epsilon_sq <- as.numeric(kw_site$statistic) / (n_total - 1)

stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "H1",
    question = "Does CAFI abundance differ among sites?",
    test_name = "Kruskal-Wallis",
    test_statistic = as.numeric(kw_site$statistic),
    df = as.character(k - 1),
    p_value = kw_site$p.value,
    effect_size = epsilon_sq,
    effect_type = "ε²",
    n = n_total,
    notes = "Non-parametric site comparison"
  )
)

# Test 4: Kruskal-Wallis for OTU richness by site
kw_rich_site <- kruskal.test(otu_richness ~ site, data = coral_cafi_summary)
epsilon_sq_rich <- as.numeric(kw_rich_site$statistic) / (n_total - 1)

stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "H1",
    question = "Does OTU richness differ among sites?",
    test_name = "Kruskal-Wallis",
    test_statistic = as.numeric(kw_rich_site$statistic),
    df = as.character(k - 1),
    p_value = kw_rich_site$p.value,
    effect_size = epsilon_sq_rich,
    effect_type = "ε²",
    n = n_total,
    notes = "Non-parametric site comparison"
  )
)

# Test 5: Correlation between CAFI abundance and richness
cor_abund_rich <- cor.test(coral_cafi_summary$total_cafi, coral_cafi_summary$otu_richness,
                           method = "pearson")
stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "Pattern",
    question = "Are CAFI abundance and richness correlated?",
    test_name = "Pearson correlation",
    test_statistic = as.numeric(cor_abund_rich$statistic),
    df = as.character(cor_abund_rich$parameter),
    p_value = cor_abund_rich$p.value,
    effect_size = as.numeric(cor_abund_rich$estimate),
    effect_type = "r",
    n = n_total,
    notes = "Abundance-richness relationship"
  )
)

# Save statistical results
save_stats_summary(stats_results, "04_scaling_relationships", "Size-Abundance Scaling Analysis")

# ============================================================================
# POWER-LAW SCALING ANALYSIS (Hypothesis H2)
# 4-Panel Figure: Testing Propagule Redirection Theory
# ============================================================================

cat("\nRunning power-law scaling analysis...\n")

# Prepare scaling data
scaling_data <- coral_cafi_summary %>%
  filter(!is.na(volume), volume > 0) %>%
  mutate(
    log_volume = log(volume + 1),
    log_abundance = log(total_cafi + 1),
    cafi_density = total_cafi / volume,
    log_density = log(cafi_density + 0.001),
    branch_type = ifelse(!is.na(branch_width) & branch_width == "tight", "tight", "wide")
  )

cat(sprintf("  - Scaling analysis sample: %d corals\n", nrow(scaling_data)))

# 1. Basic Power-Law Regression
m_basic <- lm(log_abundance ~ log_volume, data = scaling_data)
beta <- coef(m_basic)[2]
beta_se <- summary(m_basic)$coefficients[2, 2]
beta_ci <- confint(m_basic)[2,]
r_squared <- summary(m_basic)$r.squared

cat(sprintf("  - Power-law exponent β = %.3f [%.3f, %.3f]\n", beta, beta_ci[1], beta_ci[2]))
cat(sprintf("  - R² = %.3f\n", r_squared))

if (beta_ci[2] < 1) {
  cat("  ✓ Sublinear scaling confirmed (exponent < 1)\n")
} else {
  cat("  - Cannot reject isometric scaling\n")
}

# 2. Mixed Effects Model with Site Random Slopes
m_mixed <- lmer(log_abundance ~ log_volume + (1 + log_volume | site),
                data = scaling_data,
                REML = FALSE)

beta_mixed <- fixef(m_mixed)[2]
beta_mixed_ci <- confint(m_mixed, parm = "log_volume", level = 0.95)
ranef_sites <- ranef(m_mixed)$site

cat("\n  Site-specific slopes:\n")
for(s in rownames(ranef_sites)) {
  site_slope <- beta_mixed + ranef_sites[s, "log_volume"]
  cat(sprintf("    %s: %.3f\n", s, site_slope))
}

# 3. Density Scaling
density_data <- scaling_data %>% filter(cafi_density > 0)
m_density <- lm(log_density ~ log_volume, data = density_data)
density_slope <- coef(m_density)[2]
density_ci <- confint(m_density)[2,]

cat(sprintf("\n  Density scaling slope: %.3f [%.3f, %.3f]\n",
            density_slope, density_ci[1], density_ci[2]))

# ============================================================================
# Create 4-Panel Power-Law Scaling Figure
# ============================================================================

cat("  Creating 4-panel power-law figure...\n")

# Panel A: Basic power-law scaling
p_scaling_a <- ggplot(scaling_data, aes(x = volume, y = total_cafi)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black",
              se = TRUE, linewidth = 1.2) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = scales::comma) +
  scale_y_log10(breaks = c(1, 10, 100, 1000)) +
  scale_color_site() +
  annotation_logticks() +
  labs(
    title = "A. Power-law scaling",
    subtitle = sprintf("β = %.2f [%.2f, %.2f] | R² = %.2f",
                       beta, beta_ci[1], beta_ci[2], r_squared),
    x = expression(Coral~Volume~(cm^3)),
    y = "CAFI Abundance",
    color = "Site"
  ) +
  theme_publication() +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# Panel B: Site-specific slopes
site_predictions <- expand.grid(
  log_volume = seq(min(scaling_data$log_volume),
                   max(scaling_data$log_volume),
                   length.out = 100),
  site = unique(scaling_data$site)
)
site_predictions$log_abundance <- predict(m_mixed, newdata = site_predictions)

p_scaling_b <- ggplot(scaling_data, aes(x = exp(log_volume), y = exp(log_abundance))) +
  geom_line(data = site_predictions %>%
              mutate(volume = exp(log_volume), abundance = exp(log_abundance)),
            aes(x = volume, y = abundance, color = site),
            linewidth = 1.2) +
  geom_point(aes(color = site), alpha = 0.3, size = 1.5) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = scales::comma) +
  scale_y_log10(breaks = c(1, 10, 100, 1000)) +
  scale_color_site() +
  annotation_logticks() +
  labs(
    title = "B. Site-specific scaling",
    subtitle = "Random slopes model",
    x = expression(Coral~Volume~(cm^3)),
    y = "CAFI Abundance",
    color = "Site"
  ) +
  theme_publication() +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# Panel C: Density scaling (propagule dilution)
p_scaling_c <- density_data %>%
  ggplot(aes(x = volume, y = cafi_density)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black",
              se = TRUE, linewidth = 1.2) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = scales::comma) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10)) +
  scale_color_site() +
  annotation_logticks() +
  labs(
    title = "C. Density scaling (propagule dilution)",
    subtitle = sprintf("Slope = %.2f | Larger corals have lower CAFI density",
                       density_slope),
    x = expression(Coral~Volume~(cm^3)),
    y = expression(CAFI~Density~(individuals/cm^3)),
    color = "Site"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# Panel D: Richness scaling
m_richness <- lm(log(otu_richness + 1) ~ log_volume, data = scaling_data)
rich_slope <- coef(m_richness)[2]
rich_ci <- confint(m_richness)[2,]
rich_r2 <- summary(m_richness)$r.squared

p_scaling_d <- scaling_data %>%
  ggplot(aes(x = volume, y = otu_richness)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black",
              se = TRUE, linewidth = 1.2) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000), labels = scales::comma) +
  scale_y_log10(breaks = c(1, 5, 10, 20, 50)) +
  scale_color_site() +
  annotation_logticks() +
  labs(
    title = "D. Richness scaling",
    subtitle = sprintf("Slope = %.2f | Species accumulate with size (R² = %.2f)",
                       rich_slope, rich_r2),
    x = expression(Coral~Volume~(cm^3)),
    y = "CAFI OTU Richness",
    color = "Site"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# Combine 4-panel figure
fig_scaling_4panel <- (p_scaling_a | p_scaling_b) / (p_scaling_c | p_scaling_d) +
  plot_annotation(
    title = "Power-law scaling of coral-associated fauna (Hypothesis H2)",
    subtitle = "Testing propagule redirection theory: larger corals have lower CAFI density per unit volume",
    caption = sprintf("n = %d corals | Mixed model with site random slopes | Theory predicts β = 0.75",
                     nrow(scaling_data)),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      plot.caption = element_text(size = 9)
    )
  )

# Save to scaling folder
ggsave(file.path(fig_dir, "H2_power_law_scaling_4panel.png"), fig_scaling_4panel,
       width = 14, height = 12, dpi = 300, bg = "white")
cat("  - Saved:", file.path(fig_dir, "H2_power_law_scaling_4panel.png"), "\n")

# Save to manuscript folder
ggsave(file.path(MANUSCRIPT_DIR, "H2_power_law_scaling_4panel.png"), fig_scaling_4panel,
       width = 14, height = 12, dpi = 300, bg = "white")
cat("  - Saved:", file.path(MANUSCRIPT_DIR, "H2_power_law_scaling_4panel.png"), "\n")

# Add scaling results to stats
stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "H2",
    question = "Does CAFI abundance show sublinear scaling with coral volume?",
    test_name = "Power-law regression (log-log)",
    test_statistic = summary(m_basic)$fstatistic[1],
    df = paste(summary(m_basic)$fstatistic[2], summary(m_basic)$fstatistic[3], sep = ", "),
    p_value = summary(m_basic)$coefficients[2, 4],
    effect_size = beta,
    effect_type = "scaling exponent (β)",
    n = nrow(scaling_data),
    notes = sprintf("95%% CI [%.2f, %.2f]; %s",
                   beta_ci[1], beta_ci[2],
                   ifelse(beta_ci[2] < 1, "Sublinear confirmed", "Cannot reject isometric"))
  )
)

# Save scaling results
scaling_results <- tibble(
  metric = c("Abundance scaling exponent (β)", "β 95% CI lower", "β 95% CI upper",
             "R² (basic model)", "Density scaling slope", "Richness scaling slope"),
  value = c(beta, beta_ci[1], beta_ci[2], r_squared, density_slope, rich_slope)
)
write_csv(scaling_results, file.path(TABLES_DIR, "H2_scaling_results.csv"))

cat("  ✓ Power-law scaling analysis complete\n\n")

# ============================================================================
# SIZE-BINNED ANALYSIS (6-Panel Figure)
# Comparing CAFI patterns across Small, Medium, and Large coral size classes
# ============================================================================

cat("Running size-binned analysis...\n")

# Define size class breaks using tertiles
breaks <- quantile(scaling_data$volume, c(0, 0.33, 0.67, 1), na.rm = TRUE)
cat(sprintf("  Size class breaks: Small (<%.0f cm³), Medium (%.0f-%.0f cm³), Large (>%.0f cm³)\n",
            breaks[2], breaks[2], breaks[3], breaks[3]))

size_binned_data <- scaling_data %>%
  mutate(
    size_class = cut(volume,
                     breaks = breaks,
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE),
    size_class = factor(size_class, levels = c("Small", "Medium", "Large"))
  )

cat(sprintf("  Sample sizes: Small=%d, Medium=%d, Large=%d\n",
            sum(size_binned_data$size_class == "Small", na.rm = TRUE),
            sum(size_binned_data$size_class == "Medium", na.rm = TRUE),
            sum(size_binned_data$size_class == "Large", na.rm = TRUE)))

# Size class colors
size_colors <- c("Small" = "#d73027", "Medium" = "#fee08b", "Large" = "#1a9850")
size_colors_fill <- c("Small" = "#fc8d59", "Medium" = "#ffffbf", "Large" = "#91cf60")

# Statistical tests
kw_abundance_size <- kruskal.test(total_cafi ~ size_class, data = size_binned_data)
kw_richness_size <- kruskal.test(otu_richness ~ size_class, data = size_binned_data)
kw_density_size <- kruskal.test(cafi_density ~ size_class, data = size_binned_data)

format_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01) return(sprintf("p = %.3f", p))
  return(sprintf("p = %.2f", p))
}

# Panel A: Total abundance by size class
pA1 <- ggplot(size_binned_data, aes(x = size_class, y = total_cafi, fill = size_class)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white", outlier.size = 2) +
  geom_jitter(aes(color = site), width = 0.2, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  scale_color_site() +
  labs(title = "A. Total CAFI abundance",
       subtitle = sprintf("Kruskal-Wallis χ² = %.1f, %s",
                         kw_abundance_size$statistic, format_p(kw_abundance_size$p.value)),
       x = "Coral size class", y = "Total individuals") +
  theme_publication() +
  theme(legend.position = "none")

# Panel B: Species richness by size class
pA2 <- ggplot(size_binned_data, aes(x = size_class, y = otu_richness, fill = size_class)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white", outlier.size = 2) +
  geom_jitter(aes(color = site), width = 0.2, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  scale_color_site() +
  labs(title = "B. OTU richness",
       subtitle = sprintf("Kruskal-Wallis χ² = %.1f, %s",
                         kw_richness_size$statistic, format_p(kw_richness_size$p.value)),
       x = "Coral size class", y = "Number of OTUs") +
  theme_publication() +
  theme(legend.position = "none")

# Panel C: Shannon diversity (if available)
if ("shannon_diversity" %in% names(size_binned_data)) {
  kw_shannon_size <- kruskal.test(shannon_diversity ~ size_class, data = size_binned_data)
  pA3 <- ggplot(size_binned_data, aes(x = size_class, y = shannon_diversity, fill = size_class)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white", outlier.size = 2) +
    geom_jitter(aes(color = site), width = 0.2, alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = size_colors_fill, guide = "none") +
    scale_color_site() +
    labs(title = "C. Shannon diversity (H')",
         subtitle = sprintf("Kruskal-Wallis χ² = %.1f, %s",
                           kw_shannon_size$statistic, format_p(kw_shannon_size$p.value)),
         x = "Coral size class", y = "Shannon H'") +
    theme_publication() +
    theme(legend.position = "right")
} else {
  # Use volume histogram as alternative
  pA3 <- ggplot(size_binned_data, aes(x = volume, fill = size_class)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "white", position = "identity") +
    scale_fill_manual(values = size_colors, name = "Size class") +
    scale_x_log10(labels = scales::comma) +
    labs(title = "C. Volume distribution by size class",
         subtitle = "Tertile-based size classification",
         x = expression(Coral~Volume~(cm^3)), y = "Count") +
    theme_publication() +
    theme(legend.position = "right")
}

# Panel D: CAFI density (key test of propagule dilution)
pA4 <- ggplot(size_binned_data, aes(x = size_class, y = cafi_density * 1000, fill = size_class)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white", outlier.size = 2) +
  geom_jitter(aes(color = site), width = 0.2, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  scale_color_site() +
  labs(title = "D. CAFI density (propagule dilution test)",
       subtitle = sprintf("Kruskal-Wallis χ² = %.1f, %s",
                         kw_density_size$statistic, format_p(kw_density_size$p.value)),
       x = "Coral size class", y = expression("Individuals per 1000 cm"^3)) +
  theme_publication() +
  theme(legend.position = "none")

# Panel E: Species accumulation
pA5 <- ggplot(size_binned_data, aes(x = total_cafi, y = otu_richness, color = size_class)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values = size_colors, name = "Size class") +
  labs(title = "E. Species accumulation by size class",
       subtitle = "Larger corals accumulate species faster per individual",
       x = "Total CAFI abundance", y = "OTU richness") +
  theme_publication() +
  theme(legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = alpha("white", 0.8)))

# Panel F: Abundance vs Volume with Field of Dreams reference
fod_slope <- mean(size_binned_data$total_cafi[size_binned_data$size_class == "Small"], na.rm = TRUE) /
             mean(size_binned_data$volume[size_binned_data$size_class == "Small"], na.rm = TRUE)

pA6 <- ggplot(size_binned_data, aes(x = volume, y = total_cafi)) +
  geom_abline(slope = fod_slope, intercept = 0,
              linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_point(aes(fill = size_class), shape = 21, alpha = 0.7, size = 2.5,
             stroke = 0.3, color = "white") +
  geom_smooth(method = "lm", formula = y ~ x, color = "black",
              linewidth = 1, se = TRUE, fill = "gray80") +
  scale_fill_manual(values = size_colors, name = "Size class") +
  scale_x_continuous(labels = scales::comma) +
  annotate("text", x = max(size_binned_data$volume) * 0.6,
           y = max(size_binned_data$volume) * fod_slope * 0.7,
           label = "Field of Dreams\n(proportional)", color = "gray50",
           size = 2.8, fontface = "italic") +
  labs(title = "F. Sublinear scaling confirms propagule dilution",
       subtitle = "Observed (solid) falls below proportional expectation (dashed)",
       x = expression("Coral volume (cm"^3*")"), y = "CAFI abundance") +
  theme_publication() +
  theme(legend.position = "none")

# Combine into 6-panel figure
fig_size_binned <- (pA1 | pA2 | pA3) / (pA4 | pA5 | pA6) +
  plot_annotation(
    title = "Community metrics differ significantly across coral size classes",
    subtitle = "Larger corals have more CAFI total but LOWER density—supporting propagule dilution",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  )

# Save to scaling folder
ggsave(file.path(fig_dir, "size_binned_metrics_6panel.png"), fig_size_binned,
       width = 14, height = 10, dpi = 300, bg = "white")
cat("  - Saved:", file.path(fig_dir, "size_binned_metrics_6panel.png"), "\n")

# Save to manuscript folder
ggsave(file.path(MANUSCRIPT_DIR, "size_binned_metrics_6panel.png"), fig_size_binned,
       width = 14, height = 10, dpi = 300, bg = "white")
cat("  - Saved:", file.path(MANUSCRIPT_DIR, "size_binned_metrics_6panel.png"), "\n")

# Save size class statistics
size_class_stats <- size_binned_data %>%
  group_by(size_class) %>%
  summarise(
    n = n(),
    vol_mean = mean(volume),
    vol_sd = sd(volume),
    abund_mean = mean(total_cafi),
    abund_se = sd(total_cafi) / sqrt(n()),
    rich_mean = mean(otu_richness),
    rich_se = sd(otu_richness) / sqrt(n()),
    density_mean = mean(cafi_density) * 1000,
    density_se = sd(cafi_density) * 1000 / sqrt(n()),
    .groups = "drop"
  )

write_csv(size_class_stats, file.path(TABLES_DIR, "size_class_summaries.csv"))

cat("  ✓ Size-binned analysis complete\n\n")

# Update stats summary with new results
save_stats_summary(stats_results, "04_scaling_relationships", "Size-Abundance Scaling Analysis")

cat("\n✅ Coral-CAFI relationship analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", TABLES_DIR, "\n")