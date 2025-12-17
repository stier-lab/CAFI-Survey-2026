#!/usr/bin/env Rscript
# ============================================================================
# 09_functional_groups.R - Functional Group Analysis
# ============================================================================
#
# PURPOSE: Analyze CAFI functional group patterns, especially Trapezia defenders
#
# FUNCTIONAL GROUPS:
#   - Trapezia (guardian crabs): Defend coral from predators, remove sediment
#   - Resident fish: Nutrient subsidies via waste
#   - Corallivores: Tissue consumption (negative effects)
#   - Other cryptofauna: Mixed/unknown effects
#
# ANALYSES:
#   - Trapezia abundance patterns by coral size/site
#   - Functional group scaling relationships
#   - Group-specific responses to landscape
#
# MANUSCRIPT FIGURES:
#   >>> FIGURE 3: Functional Group Responses <<<
#   - Panel A: Trapezia scaling with size
#   - Panel B: Fish scaling with size
#   - Panel C: Group-specific slope coefficients
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("09: Functional Group Analysis\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")
functional_summary <- load_object("functional_summary.rds")

# Load physio_clean if available (for trapezid-physiology analysis)
physio_file <- file.path(OBJECTS_DIR, "physio_clean.rds")
if (file.exists(physio_file)) {
  physio_clean <- load_object("physio_clean.rds")
} else {
  physio_clean <- data.frame()  # Empty placeholder
}

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "functional_groups")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Identify Trapezids
# ============================================================================

cat("Identifying Trapezid crabs...\n")

# Trapezids: Family Trapeziidae - guardian crabs of Pocillopora
# In the CAFI data, these typically have code starting with "TR" or are
# identified as Trapeziidae in the family/order columns

# Check what identification methods we have
if("family" %in% names(cafi_clean)) {
  # Use family column if available
  trapezids <- cafi_clean %>%
    filter(grepl("trapez", family, ignore.case = TRUE) |
           grepl("^TR", code, ignore.case = TRUE))

  cat("  Identified using family column\n")
} else if("code" %in% names(cafi_clean)) {
  # Fall back to code prefix
  trapezids <- cafi_clean %>%
    filter(grepl("^TR", code, ignore.case = TRUE))

  cat("  Identified using code prefix (TR)\n")
} else {
  cat("  Warning: Cannot identify Trapezids without family or code column\n")
  trapezids <- data.frame()
}

if(nrow(trapezids) > 0) {
  cat("  Total Trapezid individuals:", nrow(trapezids), "\n")
  cat("  Trapezid species/OTUs:", n_distinct(trapezids$species), "\n\n")

  # ============================================================================
  # Trapezid Metrics per Coral
  # ============================================================================

  cat("Calculating per-coral Trapezid metrics...\n")

  # Summarize Trapezids per coral
  trapezid_summary <- trapezids %>%
    group_by(coral_id) %>%
    summarise(
      n_trapezids = n(),
      trapezid_richness = n_distinct(species),
      trapezid_species = paste(sort(unique(species)), collapse = "; "),
      mean_trap_size = mean(size_mm, na.rm = TRUE),
      .groups = "drop"
    )

  # Join with all corals (add 0s for corals without Trapezids)
  trap_data <- coral_master %>%
    left_join(trapezid_summary, by = "coral_id") %>%
    mutate(
      n_trapezids = replace_na(n_trapezids, 0),
      trapezid_richness = replace_na(trapezid_richness, 0)
    )

  # Add coral size metrics
  if(all(c("volume_field", "volume_lab") %in% names(trap_data))) {
    trap_data <- trap_data %>%
      mutate(coral_volume = coalesce(volume_field, volume_lab))
  }

  cat("✓ Trapezid metrics calculated for", nrow(trap_data), "corals\n\n")

  # ============================================================================
  # Trapezid Abundance Patterns
  # ============================================================================

  cat("Analyzing Trapezid abundance patterns...\n")

  # Trapezid abundance by site
  p_trap_site <- trap_data %>%
    filter(!is.na(site)) %>%
    ggplot(aes(x = site, y = n_trapezids)) +
    geom_boxplot(aes(fill = site), alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
    scale_fill_site() +
    labs(
      title = "Trapezid crab abundance by site",
      subtitle = "Family Trapeziidae - Guardian crabs of Pocillopora",
      x = NULL,
      y = "Trapezid abundance"
    ) +
    theme_publication() +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "trapezids_by_site.png"),
         p_trap_site, width = 10, height = 6, dpi = 300, bg = "white")

  # Trapezid abundance by morphotype
  if("morphotype" %in% names(trap_data) && sum(!is.na(trap_data$morphotype)) > 0) {
    p_trap_morph <- trap_data %>%
      filter(!is.na(morphotype)) %>%
      ggplot(aes(x = morphotype, y = n_trapezids)) +
      geom_boxplot(aes(fill = morphotype), alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
      scale_fill_viridis_d() +
      labs(
        title = "Trapezid crab abundance by coral morphotype",
        x = "Morphotype",
        y = "Trapezid abundance"
      ) +
      theme_publication() +
      theme(legend.position = "none")

    ggsave(file.path(fig_dir, "trapezids_by_morphotype.png"),
           p_trap_morph, width = 10, height = 6, dpi = 300, bg = "white")
  }

  # ============================================================================
  # Trapezid vs Coral Volume (inspired by Curtis 2019)
  # ============================================================================

  if("coral_volume" %in% names(trap_data)) {
    cat("Analyzing Trapezid-coral size relationships...\n")

    # Trapezid abundance vs coral volume
    p_trap_volume <- trap_data %>%
      filter(!is.na(coral_volume), coral_volume > 0) %>%
      ggplot(aes(x = coral_volume, y = n_trapezids)) +
      geom_point(aes(color = site, shape = site), size = 2.5, alpha = 0.7) +
      geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 0.8) +
      scale_color_site() +
      scale_shape_site() +
      scale_x_log10() +
      labs(
        title = "Trapezid abundance vs coral volume",
        subtitle = "Inspired by Curtis 2019 analysis",
        x = expression("Coral volume (cm"^3*", log scale)"),
        y = "Trapezid abundance"
      ) +
      theme_publication() +
      theme(legend.position = "right")

    ggsave(file.path(fig_dir, "trapezids_vs_coral_volume.png"),
           p_trap_volume, width = 10, height = 8, dpi = 300, bg = "white")

    # Faceted by branch width (if available)
    if("branch_width" %in% names(trap_data) && sum(!is.na(trap_data$branch_width)) > 5) {
      p_trap_volume_branch <- trap_data %>%
        filter(!is.na(coral_volume), !is.na(branch_width), coral_volume > 0) %>%
        ggplot(aes(x = coral_volume, y = n_trapezids)) +
        geom_point(aes(color = morphotype), size = 3, alpha = 0.6) +
        geom_smooth(method = "lm", color = "black", se = TRUE) +
        facet_wrap(~branch_width) +
        scale_color_viridis_d() +
        scale_x_log10() +
        labs(
          title = "Trapezid Abundance vs Coral Volume by Branch Width",
          x = expression("Coral Volume (cm"^3*", log scale)"),
          y = "Trapezid Abundance"
        )

      ggsave(file.path(fig_dir, "trapezids_vs_volume_by_branch.png"),
             p_trap_volume_branch, width = 12, height = 6, dpi = 300)
    }

    # Trapezid richness vs coral volume
    p_traprich_volume <- trap_data %>%
      filter(!is.na(coral_volume), coral_volume > 0) %>%
      ggplot(aes(x = coral_volume, y = trapezid_richness)) +
      geom_point(aes(color = site), size = 3, alpha = 0.6) +
      geom_smooth(method = "lm", color = "black", se = TRUE) +
      scale_color_viridis_d() +
      scale_x_log10() +
      labs(
        title = "Trapezid Richness vs Coral Volume",
        x = expression("Coral Volume (cm"^3*", log scale)"),
        y = "Trapezid Species Richness"
      )

    ggsave(file.path(fig_dir, "trapezid_richness_vs_volume.png"),
           p_traprich_volume, width = 10, height = 8, dpi = 300)
  }

  # ============================================================================
  # Trapezid vs Physiology (if available)
  # ============================================================================

  trap_physio <- trap_data %>%
    inner_join(physio_clean, by = "coral_id")

  if(nrow(trap_physio) > 10 && "afdw_mg_cm2" %in% names(trap_physio)) {
    cat("Analyzing Trapezid-physiology relationships...\n")

    # Trapezid abundance vs tissue biomass
    p_trap_biomass <- trap_physio %>%
      filter(!is.na(afdw_mg_cm2), !is.na(n_trapezids)) %>%
      ggplot(aes(x = n_trapezids, y = afdw_mg_cm2)) +
      geom_point(aes(color = morphotype), size = 3, alpha = 0.6) +
      geom_smooth(method = "lm", color = "black", se = TRUE) +
      scale_color_viridis_d() +
      labs(
        title = "Trapezid Abundance vs Coral Tissue Biomass",
        subtitle = "Testing for mutualist effects on coral health",
        x = "Trapezid Abundance",
        y = expression("Tissue Biomass (AFDW mg/cm"^2*")")
      )

    ggsave(file.path(fig_dir, "trapezids_vs_tissue_biomass.png"),
           p_trap_biomass, width = 10, height = 8, dpi = 300)
  }

  # ============================================================================
  # Summary Statistics
  # ============================================================================

  cat("Generating summary statistics...\n")

  # Overall Trapezid summary
  trap_stats <- data.frame(
    metric = c(
      "Total Trapezid individuals",
      "Trapezid species/OTUs",
      "Corals with Trapezids",
      "Percent corals with Trapezids",
      "Mean Trapezids per coral (all)",
      "Mean Trapezids per coral (occupied)",
      "Max Trapezids on single coral"
    ),
    value = c(
      nrow(trapezids),
      n_distinct(trapezids$species),
      sum(trap_data$n_trapezids > 0),
      round(sum(trap_data$n_trapezids > 0) / nrow(trap_data) * 100, 1),
      round(mean(trap_data$n_trapezids), 2),
      round(mean(trap_data$n_trapezids[trap_data$n_trapezids > 0]), 2),
      max(trap_data$n_trapezids)
    )
  )

  write_csv(trap_stats,
            file.path(TABLES_DIR, "trapezid_summary_statistics.csv"))

  # Trapezid species list
  if("species" %in% names(trapezids)) {
    trap_species <- trapezids %>%
      group_by(species, type, code) %>%
      summarise(
        n_individuals = n(),
        n_corals = n_distinct(coral_id),
        mean_size_mm = mean(size_mm, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n_individuals))

    write_csv(trap_species,
              file.path(TABLES_DIR, "trapezid_species_list.csv"))
  }

  cat("✓ Trapezid analysis complete\n\n")

  # ============================================================================
  # Summary Report
  # ============================================================================

  cat("========================================\n")
  cat("Trapezid Analysis Summary\n")
  cat("========================================\n\n")

  cat("Trapeziidae Family (Guardian Crabs):\n")
  cat("  Total individuals:", nrow(trapezids), "\n")
  cat("  Species/OTUs:", n_distinct(trapezids$species), "\n")
  cat("  Corals occupied:", sum(trap_data$n_trapezids > 0), "of", nrow(trap_data),
      sprintf("(%.1f%%)", sum(trap_data$n_trapezids > 0) / nrow(trap_data) * 100), "\n")
  cat("  Mean per coral:", sprintf("%.2f", mean(trap_data$n_trapezids)), "\n")
  cat("  Max on one coral:", max(trap_data$n_trapezids), "\n\n")

  cat("Ecological Notes:\n")
  cat("  - Trapezids are obligate symbionts of Pocillopora corals\n")
  cat("  - Act as 'guardian crabs' defending coral from predators\n")
  cat("  - Presence may indicate coral health and habitat quality\n")
  cat("  - Larger corals typically host more Trapezid individuals\n\n")

  cat("Figures generated in:", fig_dir, "\n")
  cat("Tables generated in:", TABLES_DIR, "\n\n")

} else {
  cat("No Trapezid crabs identified in dataset\n")
  cat("  Check that CAFI data includes family classification or TR codes\n\n")
}

# ============================================================================
# STATISTICAL RESULTS SUMMARY
# ============================================================================

cat("Compiling statistical results...\n")

stats_results <- init_results_df()

if(nrow(trapezids) > 0 && exists("trap_data")) {

  # 1. Kruskal-Wallis test for trapezid abundance by site
  if(sum(!is.na(trap_data$site)) > 10) {
    kw_trap_site <- kruskal.test(n_trapezids ~ site,
                                  data = trap_data %>% filter(!is.na(site)))
    n_total <- nrow(trap_data %>% filter(!is.na(site)))
    k <- n_distinct(trap_data$site[!is.na(trap_data$site)])
    epsilon_sq <- kw_trap_site$statistic / (n_total - 1)

    stats_results <- bind_rows(stats_results,
      create_result_row(
        hypothesis = "H_trap_site",
        question = "Does Trapezid abundance differ among sites?",
        test_name = "Kruskal-Wallis",
        test_statistic = kw_trap_site$statistic,
        df = k - 1,
        p_value = kw_trap_site$p.value,
        effect_size = as.numeric(epsilon_sq),
        effect_type = "ε²",
        n = n_total,
        notes = "Non-parametric test for Trapezid abundance"
      )
    )
  }

  # 2. Trapezid-volume correlation (if coral_volume exists)
  if("coral_volume" %in% names(trap_data)) {
    trap_vol_data <- trap_data %>%
      filter(!is.na(coral_volume), coral_volume > 0)

    if(nrow(trap_vol_data) > 10) {
      cor_trap_vol <- cor.test(log10(trap_vol_data$coral_volume),
                                trap_vol_data$n_trapezids,
                                method = "spearman")

      stats_results <- bind_rows(stats_results,
        create_result_row(
          hypothesis = "H_trap_size",
          question = "Does Trapezid abundance scale with coral volume?",
          test_name = "Spearman correlation",
          test_statistic = cor_trap_vol$statistic,
          df = nrow(trap_vol_data) - 2,
          p_value = cor_trap_vol$p.value,
          effect_size = cor_trap_vol$estimate,
          effect_type = "rho",
          n = nrow(trap_vol_data),
          notes = "Correlation: log(volume) vs Trapezid count"
        )
      )

      # Linear regression for scaling
      trap_model <- lm(n_trapezids ~ log10(coral_volume),
                       data = trap_vol_data)
      trap_summary <- summary(trap_model)

      stats_results <- bind_rows(stats_results,
        create_result_row(
          hypothesis = "H_trap_scaling",
          question = "What is the scaling relationship between coral size and Trapezid abundance?",
          test_name = "Linear regression",
          test_statistic = trap_summary$coefficients["log10(coral_volume)", "t value"],
          df = trap_model$df.residual,
          p_value = trap_summary$coefficients["log10(coral_volume)", "Pr(>|t|)"],
          effect_size = trap_summary$r.squared,
          effect_type = "R²",
          n = nrow(trap_vol_data),
          notes = sprintf("β=%.3f", trap_summary$coefficients["log10(coral_volume)", "Estimate"])
        )
      )
    }
  }

  # 3. Trapezid-physiology relationship (if data exists)
  if(exists("trap_physio") && nrow(trap_physio) > 10 && "afdw_mg_cm2" %in% names(trap_physio)) {
    trap_physio_data <- trap_physio %>%
      filter(!is.na(afdw_mg_cm2), !is.na(n_trapezids))

    if(nrow(trap_physio_data) > 10) {
      cor_trap_biomass <- cor.test(trap_physio_data$n_trapezids,
                                    trap_physio_data$afdw_mg_cm2,
                                    method = "spearman")

      stats_results <- bind_rows(stats_results,
        create_result_row(
          hypothesis = "H_trap_health",
          question = "Is Trapezid abundance associated with coral tissue biomass?",
          test_name = "Spearman correlation",
          test_statistic = cor_trap_biomass$statistic,
          df = nrow(trap_physio_data) - 2,
          p_value = cor_trap_biomass$p.value,
          effect_size = cor_trap_biomass$estimate,
          effect_type = "rho",
          n = nrow(trap_physio_data),
          notes = "Trapezid count vs AFDW (tissue biomass)"
        )
      )
    }
  }

  # 4. Chi-square for trapezid presence/absence by site
  if(sum(!is.na(trap_data$site)) > 10) {
    trap_presence <- trap_data %>%
      filter(!is.na(site)) %>%
      mutate(has_trapezid = n_trapezids > 0)

    if(n_distinct(trap_presence$site) > 1 &&
       sum(trap_presence$has_trapezid) > 0 &&
       sum(!trap_presence$has_trapezid) > 0) {

      chi_trap <- chisq.test(table(trap_presence$site, trap_presence$has_trapezid))

      # Cramér's V
      n <- nrow(trap_presence)
      cramers_v <- sqrt(chi_trap$statistic / (n * (min(2, n_distinct(trap_presence$site)) - 1)))

      stats_results <- bind_rows(stats_results,
        create_result_row(
          hypothesis = "H_trap_occupancy",
          question = "Does Trapezid occupancy rate differ among sites?",
          test_name = "Chi-square test",
          test_statistic = chi_trap$statistic,
          df = chi_trap$parameter,
          p_value = chi_trap$p.value,
          effect_size = as.numeric(cramers_v),
          effect_type = "Cramér's V",
          n = nrow(trap_presence),
          notes = sprintf("%.1f%% overall occupancy", mean(trap_presence$has_trapezid) * 100)
        )
      )
    }
  }
}

# Save statistical results
save_stats_summary(stats_results, "09_functional_groups", "Functional Group Analysis")

cat("✅ Trapezid guild analysis complete!\n\n")
