#!/usr/bin/env Rscript
# ============================================================================
# 17_trapezid_guild_analysis.R - Trapezid crab guild-specific analyses
# Author: CAFI Analysis Pipeline
# Date: 2025-10-31
# Purpose: Analyze Trapeziidae crab family (guardian crabs) patterns
# Inspired by: Curtis 2019 original Survey analysis
# ============================================================================

cat("\n========================================\n")
cat("Trapezid Guild Analysis\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/00_load_libraries.R"))

# Load processed data
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))
physio_clean <- readRDS(file.path(SURVEY_OBJECTS, "physio_clean.rds"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "trapezid_analysis")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

table_dir <- SURVEY_TABLES

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
  trap_data <- metadata %>%
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
    geom_boxplot(aes(fill = site), alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    scale_fill_viridis_d() +
    labs(
      title = "Trapezid Crab Abundance by Site",
      subtitle = "Family Trapeziidae - Guardian crabs of Pocillopora",
      x = "Site",
      y = "Trapezid Abundance"
    ) +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "trapezids_by_site.png"),
         p_trap_site, width = 10, height = 6, dpi = 300)

  # Trapezid abundance by morphotype
  if("morphotype" %in% names(trap_data) && sum(!is.na(trap_data$morphotype)) > 0) {
    p_trap_morph <- trap_data %>%
      filter(!is.na(morphotype)) %>%
      ggplot(aes(x = morphotype, y = n_trapezids)) +
      geom_boxplot(aes(fill = morphotype), alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.3) +
      scale_fill_viridis_d() +
      labs(
        title = "Trapezid Crab Abundance by Coral Morphotype",
        x = "Morphotype",
        y = "Trapezid Abundance"
      ) +
      theme(legend.position = "none")

    ggsave(file.path(fig_dir, "trapezids_by_morphotype.png"),
           p_trap_morph, width = 10, height = 6, dpi = 300)
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
      geom_point(aes(color = site), size = 3, alpha = 0.6) +
      geom_smooth(method = "lm", color = "black", se = TRUE) +
      scale_color_viridis_d() +
      scale_x_log10() +
      labs(
        title = "Trapezid Abundance vs Coral Volume",
        subtitle = "Inspired by Curtis 2019 analysis",
        x = expression("Coral Volume (cm"^3*", log scale)"),
        y = "Trapezid Abundance"
      )

    ggsave(file.path(fig_dir, "trapezids_vs_coral_volume.png"),
           p_trap_volume, width = 10, height = 8, dpi = 300)

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
            file.path(table_dir, "trapezid_summary_statistics.csv"))

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
              file.path(table_dir, "trapezid_species_list.csv"))
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
  cat("Tables generated in:", table_dir, "\n\n")

} else {
  cat("No Trapezid crabs identified in dataset\n")
  cat("  Check that CAFI data includes family classification or TR codes\n\n")
}

cat("✅ Trapezid guild analysis complete!\n\n")
