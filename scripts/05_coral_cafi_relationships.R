#!/usr/bin/env Rscript
# ============================================================================
# 05_coral_cafi_relationships.R - Coral-CAFI Relationships (H2, H4)
# Author: CAFI Analysis Pipeline
# Date: 2025-11-01
#
# Hypotheses tested (aligned with PRD):
#   H2: CAFI abundance scales with coral volume following a power-law with
#       exponent < 1, indicating larger corals have lower CAFI densities
#   H4: Coral physiological condition positively predicts CAFI diversity
#
# Theoretical Background:
#   Propagule redirection theory predicts that occupant abundance scales
#   nonlinearly with habitat amount, causing occupant density (per unit
#   habitat) to decrease as habitat increases. Larvae distribute among
#   available habitats based on chemical cue strength, so larger corals
#   intercept more propagules but at lower density per unit volume.
#   Expected scaling exponent ~0.75 for 3D habitat.
#
# Key Analyses:
#   - CAFI abundance patterns across sites
#   - Size-abundance relationships (power-law scaling)
#   - Branch architecture effects
#   - Condition-diversity relationships
#   - Taxonomic group composition
#
# Important Notes:
#   - All corals are Pocillopora spp. (cannot reliably distinguish
#     morphotypes/species in field)
#   - branch_width is a real measured trait (wide vs tight branching)
#   - "OTU" = Operational Taxonomic Unit (genetic cluster, not necessarily
#     a true species without further validation)
#   - Analysis covers 3 locations in Mo'orea: Hauru (HAU), Maatea (MAT),
#     Moorea Barrier Reef (MRB)
# ============================================================================

cat("\n========================================\n")
cat("Coral-CAFI Relationship Analysis\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/00_load_libraries.R"))

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
physio_clean <- readRDS(file.path(SURVEY_OBJECTS, "physio_clean.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))

# Load position-corrected condition scores from Script 05a
condition_scores <- readRDS(file.path(SURVEY_OBJECTS, "coral_condition_scores.rds"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "coral_cafi_relationships")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Use publication theme (loaded via 00_load_libraries.R)

# ============================================================================
# CAFI Abundance vs Coral Characteristics
# ============================================================================

cat("Analyzing CAFI abundance patterns...\n")

# Prepare coral-level summary
coral_cafi_summary <- cafi_clean %>%
  group_by(coral_id) %>%
  summarise(
    total_cafi = n(),
    otu_richness = n_distinct(species),
    n_crabs = sum(type == "crab"),
    n_shrimps = sum(type == "shrimp"),
    n_fish = sum(type == "fish"),
    n_snails = sum(type == "snail"),
    mean_size = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(metadata, by = "coral_id") %>%
  mutate(across(c(total_cafi:n_snails), ~replace_na(., 0)))

# CAFI abundance across all Pocillopora - overall distribution
p_overall_cafi <- coral_cafi_summary %>%
  ggplot(aes(x = "Pocillopora spp.", y = total_cafi)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 2, color = "gray30") +
  scale_y_sqrt(breaks = c(0, 25, 100, 225, 400)) +
  labs(
    title = "CAFI Abundance on Pocillopora spp.",
    subtitle = "All corals pooled (morphotypes cannot be reliably distinguished)",
    x = "",
    y = "CAFI Abundance (sqrt scale)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.text.x = element_text(face = "bold", size = 12)
  )

ggsave(file.path(fig_dir, "cafi_overall_abundance.png"),
       p_overall_cafi, width = 8, height = 6, dpi = 300, bg = "white")

# CAFI by location (3 sites)
p_site_cafi <- coral_cafi_summary %>%
  ggplot(aes(x = site, y = total_cafi)) +
  geom_boxplot(aes(fill = site), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
  scale_fill_site() +
  scale_y_sqrt(breaks = c(0, 25, 100, 225, 400)) +
  labs(
    title = "CAFI Abundance by Location",
    subtitle = "Three locations in Mo'orea on Pocillopora spp.",
    x = "Location",
    y = "CAFI Abundance (sqrt scale)"
  ) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "cafi_by_location.png"),
       p_site_cafi, width = 10, height = 6, dpi = 300, bg = "white")

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

# Plot overall taxonomic composition
p_tax_overall <- ggplot(taxonomic_composition, aes(x = "", y = proportion, fill = group)) +
  geom_col(width = 1, color = "white", linewidth = 1) +
  coord_polar(theta = "y") +
  scale_fill_viridis_d(option = "plasma") +
  labs(
    title = "Overall CAFI Taxonomic Composition",
    subtitle = "Proportions across all Pocillopora spp. corals",
    fill = "Taxonomic Group"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),
    legend.position = "right"
  )

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
  geom_col(position = "stack", alpha = 0.9) +
  scale_fill_manual(values = c("Crabs" = "#E69F00", "Shrimps" = "#56B4E9",
                                "Fish" = "#009E73", "Snails" = "#CC79A7")) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "CAFI Taxonomic Composition by Location",
    subtitle = "Relative proportions of functional groups at each Mo'orea location",
    x = "Location",
    y = "Proportion of CAFI Community",
    fill = "Taxonomic Group"
  ) +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "taxonomic_composition_by_site.png"),
       p_tax_site, width = 10, height = 6, dpi = 300, bg = "white")

# ============================================================================
# CAFI-Coral Condition Relationships (Using Position-Corrected Scores)
# ============================================================================

cat("Analyzing CAFI-coral condition relationships...\n")
cat("  Using position-corrected condition scores from Script 05a\n")

# Merge CAFI and condition score data
cafi_condition <- coral_cafi_summary %>%
  inner_join(condition_scores %>% select(coral_id, condition_score, site),
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
    geom_point(aes(color = site), size = 3, alpha = 0.6) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    facet_wrap(~site, ncol = 3) +
    scale_color_viridis_d(option = "plasma") +
    scale_x_sqrt(breaks = c(0, 25, 100, 225)) +
    labs(
      title = "CAFI Abundance vs Coral Condition by Location",
      subtitle = "Position-corrected condition score across 3 Mo'orea sites",
      x = "CAFI Abundance (sqrt scale)",
      y = "Coral Condition Score"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

  ggsave(file.path(fig_dir, "cafi_vs_condition_by_site.png"),
         p_cafi_condition_site, width = 14, height = 5, dpi = 300, bg = "white")

  # ============================================================================
  # CAFI Richness vs Coral Condition
  # ============================================================================

  p_richness_condition <- cafi_condition %>%
    filter(!is.na(condition_score), !is.na(otu_richness)) %>%
    ggplot(aes(x = otu_richness, y = condition_score)) +
    geom_point(aes(color = site), size = 3, alpha = 0.6) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    scale_color_viridis_d(option = "plasma", name = "Location") +
    labs(
      title = "CAFI OTU Richness vs Coral Condition",
      subtitle = "Do corals with diverse CAFI communities have better physiological condition?",
      x = "CAFI OTU Richness",
      y = "Coral Condition Score (position-corrected)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      legend.position = "right"
    )

  ggsave(file.path(fig_dir, "richness_vs_coral_condition.png"),
         p_richness_condition, width = 11, height = 7, dpi = 300, bg = "white")

  # ============================================================================
  # Taxonomic Group Abundance vs Condition
  # ============================================================================

  # Panel showing different CAFI groups
  p_crabs <- cafi_condition %>%
    filter(!is.na(condition_score)) %>%
    ggplot(aes(x = n_crabs, y = condition_score)) +
    geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(x = "Crab Abundance", y = "Condition Score")

  p_shrimps <- cafi_condition %>%
    filter(!is.na(condition_score)) %>%
    ggplot(aes(x = n_shrimps, y = condition_score)) +
    geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(x = "Shrimp Abundance", y = "Condition Score")

  p_fish <- cafi_condition %>%
    filter(!is.na(condition_score)) %>%
    ggplot(aes(x = n_fish, y = condition_score)) +
    geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(x = "Fish Abundance", y = "Condition Score")

  p_snails <- cafi_condition %>%
    filter(!is.na(condition_score)) %>%
    ggplot(aes(x = n_snails, y = condition_score)) +
    geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    labs(x = "Snail Abundance", y = "Condition Score")

  p_taxa_panel <- (p_crabs | p_shrimps) / (p_fish | p_snails) +
    plot_annotation(
      title = "CAFI Taxonomic Groups vs Coral Condition",
      subtitle = "Position-corrected condition score by major CAFI taxa",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray30")
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
    inner_join(condition_scores %>% select(coral_id, condition_score,
                                           protein_corrected, carb_corrected,
                                           zoox_corrected, afdw_corrected),
              by = "coral_id")

  if(nrow(cafi_traits) > 10) {

    # Create 4-panel plot: CAFI abundance vs each corrected trait
    p_protein_corr <- cafi_traits %>%
      filter(!is.na(protein_corrected)) %>%
      ggplot(aes(x = total_cafi, y = protein_corrected)) +
      geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      scale_x_sqrt() +
      labs(x = "CAFI Abundance", y = "Protein (pos-corrected z-score)") +
      theme_minimal()

    p_carb_corr <- cafi_traits %>%
      filter(!is.na(carb_corrected)) %>%
      ggplot(aes(x = total_cafi, y = carb_corrected)) +
      geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      scale_x_sqrt() +
      labs(x = "CAFI Abundance", y = "Carbohydrate (pos-corrected z-score)") +
      theme_minimal()

    p_zoox_corr <- cafi_traits %>%
      filter(!is.na(zoox_corrected)) %>%
      ggplot(aes(x = total_cafi, y = zoox_corrected)) +
      geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      scale_x_sqrt() +
      labs(x = "CAFI Abundance", y = "Zooxanthellae (pos-corrected z-score)") +
      theme_minimal()

    p_afdw_corr <- cafi_traits %>%
      filter(!is.na(afdw_corrected)) %>%
      ggplot(aes(x = total_cafi, y = afdw_corrected)) +
      geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      scale_x_sqrt() +
      labs(x = "CAFI Abundance", y = "AFDW (pos-corrected z-score)") +
      theme_minimal()

    p_traits_panel <- (p_protein_corr | p_carb_corr) / (p_zoox_corr | p_afdw_corr) +
      plot_annotation(
        title = "CAFI Abundance vs Position-Corrected Physiological Traits",
        subtitle = "All traits corrected for sampling position bias (residuals from trait ~ stump_length)",
        caption = "Z-scores centered at 0; positive = better than expected for branch position",
        theme = theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "gray30"),
          plot.caption = element_text(size = 9, hjust = 0)
        )
      )

    ggsave(file.path(fig_dir, "cafi_vs_corrected_traits_panel.png"),
           p_traits_panel, width = 12, height = 10, dpi = 300, bg = "white")

    # Create similar panel for OTU richness
    p_protein_rich <- cafi_traits %>%
      filter(!is.na(protein_corrected)) %>%
      ggplot(aes(x = otu_richness, y = protein_corrected)) +
      geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      labs(x = "OTU Richness", y = "Protein (z-score)") +
      theme_minimal()

    p_carb_rich <- cafi_traits %>%
      filter(!is.na(carb_corrected)) %>%
      ggplot(aes(x = otu_richness, y = carb_corrected)) +
      geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      labs(x = "OTU Richness", y = "Carbohydrate (z-score)") +
      theme_minimal()

    p_zoox_rich <- cafi_traits %>%
      filter(!is.na(zoox_corrected)) %>%
      ggplot(aes(x = otu_richness, y = zoox_corrected)) +
      geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      labs(x = "OTU Richness", y = "Zooxanthellae (z-score)") +
      theme_minimal()

    p_afdw_rich <- cafi_traits %>%
      filter(!is.na(afdw_corrected)) %>%
      ggplot(aes(x = otu_richness, y = afdw_corrected)) +
      geom_point(alpha = 0.6, size = 2.5, color = "steelblue") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      labs(x = "OTU Richness", y = "AFDW (z-score)") +
      theme_minimal()

    p_richness_traits_panel <- (p_protein_rich | p_carb_rich) / (p_zoox_rich | p_afdw_rich) +
      plot_annotation(
        title = "CAFI OTU Richness vs Position-Corrected Physiological Traits",
        subtitle = "All traits corrected for sampling position bias",
        theme = theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "gray30")
        )
      )

    ggsave(file.path(fig_dir, "richness_vs_corrected_traits_panel.png"),
           p_richness_traits_panel, width = 12, height = 10, dpi = 300, bg = "white")

    # =========================================================================
    # Statistical models for each corrected trait
    # =========================================================================

    cat("  Running models for individual corrected traits...\n")

    trait_models <- list()
    trait_names <- c("protein_corrected", "carb_corrected", "zoox_corrected", "afdw_corrected")
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
              file.path(SURVEY_TABLES, "cafi_corrected_traits_models.csv"))

    # Visualization of all trait model coefficients
    p_trait_coefs <- ggplot(all_trait_models,
                           aes(x = estimate, y = trait, color = predictor)) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                     height = 0.3, linewidth = 1, position = position_dodge(width = 0.5)) +
      geom_point(size = 4, position = position_dodge(width = 0.5)) +
      scale_color_viridis_d(option = "plasma", name = "Predictor") +
      labs(
        title = "CAFI Effects on Individual Position-Corrected Traits",
        subtitle = "Controlling for site differences | All traits position-corrected",
        x = "Coefficient Estimate (± 95% CI)",
        y = "Physiological Trait",
        caption = "Models: trait ~ CAFI_metric + site"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray30"),
        plot.caption = element_text(size = 9, hjust = 0),
        legend.position = "bottom"
      )

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
          file.path(SURVEY_OBJECTS, "cafi_condition_models.rds"))

  # Combine for visualization
  coef_data <- bind_rows(
    broom::tidy(model_cafi_condition, conf.int = TRUE) %>% mutate(model = "CAFI Abundance"),
    broom::tidy(model_richness_condition, conf.int = TRUE) %>% mutate(model = "OTU Richness")
  ) %>%
    filter(term != "(Intercept)")

  p_model_coef <- ggplot(coef_data, aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = model),
                   height = 0.2, linewidth = 1) +
    geom_point(aes(color = model), size = 4) +
    facet_wrap(~model, scales = "free_x", ncol = 2) +
    scale_color_viridis_d(option = "plasma") +
    labs(
      title = "Effect of Coral Condition on CAFI Patterns",
      subtitle = "Linear models: CAFI ~ condition_score + site",
      x = "Coefficient Estimate (± 95% CI)",
      y = "Predictor",
      caption = "Condition score = position-corrected PC1 from Script 05a"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      plot.caption = element_text(size = 9, hjust = 0),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "none"
    )

  ggsave(file.path(fig_dir, "condition_model_coefficients.png"),
         p_model_coef, width = 12, height = 7, dpi = 300, bg = "white")

  # Save coefficient table
  write_csv(coef_data,
            file.path(SURVEY_TABLES, "cafi_condition_model_coefficients.csv"))

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
  model_summaries <- list(
    abundance_model = broom::tidy(model_abundance),
    richness_model = broom::tidy(model_richness)
  )

  saveRDS(model_summaries,
          file.path(SURVEY_OBJECTS, "coral_cafi_models.rds"))

  # Create coefficient plot
  coef_data <- bind_rows(
    broom::tidy(model_abundance, conf.int = TRUE) %>% mutate(model = "Abundance"),
    broom::tidy(model_richness, conf.int = TRUE) %>% mutate(model = "OTU Richness")
  ) %>%
    filter(term != "(Intercept)")

  p_coef <- ggplot(coef_data, aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, linewidth = 0.8) +
    geom_point(aes(color = model), size = 3) +
    facet_wrap(~model, scales = "free_x") +
    scale_color_viridis_d(option = "plasma") +
    labs(
      title = "Model Coefficients for CAFI Patterns on Pocillopora spp.",
      subtitle = "Poisson GLMs for abundance and OTU richness across 3 Mo'orea locations",
      x = "Coefficient Estimate (log scale)",
      y = "Predictor"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

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
  geom_density(alpha = 0.6) +
  facet_wrap(~type, scales = "free_y", ncol = 2) +
  scale_fill_viridis_d(option = "plasma") +
  scale_x_continuous(limits = c(0, 50)) +  # Focus on main size range
  labs(
    title = "CAFI Size Distribution by Taxonomic Group",
    subtitle = "Four main functional groups on Pocillopora spp.",
    x = "Body Size (mm)",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

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
          file.path(SURVEY_TABLES, "coral_cafi_summary_by_site.csv"))

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

cat("\n✅ Coral-CAFI relationship analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", SURVEY_TABLES, "\n")