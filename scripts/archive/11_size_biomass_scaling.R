#!/usr/bin/env Rscript
# ============================================================================
# 11_size_biomass_scaling.R - CAFI Size Structure & Biomass Scaling
# ============================================================================
#
# PURPOSE: Analyze CAFI body size distributions and biomass estimates
#
# ANALYSES:
#   - Size frequency distributions by taxa
#   - Length-weight conversions
#   - Biomass scaling relationships
#   - Size spectra analysis
#
# MANUSCRIPT FIGURES: None (supplementary)
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("11: Size-Biomass Scaling Analysis\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "biomass")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Size-Biomass Conversions
# ============================================================================

cat("Converting size to biomass estimates...\n")

# Length-weight relationships from literature (example coefficients)
# W = a * L^b where W is weight in mg, L is length in mm
lw_relationships <- data.frame(
  type = c("crab", "shrimp", "fish", "snail"),
  a = c(0.0012, 0.0008, 0.0095, 0.0023),  # These are example values
  b = c(2.89, 2.95, 3.02, 2.67)
)

# Calculate biomass for each individual
cafi_biomass <- cafi_clean %>%
  left_join(lw_relationships, by = "type") %>%
  mutate(
    biomass_mg = a * (size_mm ^ b),
    biomass_g = biomass_mg / 1000,
    log_size = log10(size_mm),
    log_biomass = log10(biomass_mg)
  )

# Aggregate biomass by coral
coral_biomass <- cafi_biomass %>%
  group_by(coral_id) %>%
  summarise(
    total_biomass_g = sum(biomass_g, na.rm = TRUE),
    mean_ind_biomass_mg = mean(biomass_mg, na.rm = TRUE),
    total_abundance = n(),
    biomass_per_capita = total_biomass_g / total_abundance,
    .groups = "drop"
  )

write_csv(coral_biomass,
          file.path(TABLES_DIR, "coral_cafi_biomass_estimates.csv"))

cat("✓ Biomass calculations complete\n\n")

# ============================================================================
# Size Spectra Analysis
# ============================================================================

cat("Analyzing size spectra...\n")

# Create size bins (log2 scale)
size_bins <- 2^seq(log2(min(cafi_biomass$size_mm, na.rm = TRUE)),
                   log2(max(cafi_biomass$size_mm, na.rm = TRUE)),
                   length.out = 20)

# Calculate normalized size spectra
size_spectra <- cafi_biomass %>%
  mutate(size_class = cut(size_mm, breaks = size_bins, include.lowest = TRUE)) %>%
  group_by(size_class) %>%
  summarise(
    abundance = n(),
    total_biomass = sum(biomass_mg, na.rm = TRUE),
    mean_size = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(size_class)) %>%
  mutate(
    bin_width = diff(size_bins)[as.numeric(size_class)],
    normalized_abundance = abundance / bin_width,
    normalized_biomass = total_biomass / bin_width,
    log_size = log10(mean_size),
    log_norm_abundance = log10(normalized_abundance + 1),
    log_norm_biomass = log10(normalized_biomass + 1)
  )

# Fit power law to size spectrum
if (nrow(size_spectra) > 5) {
  spectrum_model <- lm(log_norm_abundance ~ log_size, data = size_spectra)
  spectrum_slope <- coef(spectrum_model)[2]
  spectrum_r2 <- summary(spectrum_model)$r.squared

  cat("  Size spectrum slope:", round(spectrum_slope, 2), "\n")
  cat("  R-squared:", round(spectrum_r2, 3), "\n\n")

  # Plot size spectrum
  p_spectrum <- ggplot(size_spectra, aes(x = mean_size, y = normalized_abundance)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.8) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(
      title = "Community size spectrum",
      subtitle = paste("Slope =", round(spectrum_slope, 2), ", R² =", round(spectrum_r2, 3)),
      x = "Size (mm)",
      y = "Normalized abundance (N/mm)"
    ) +
    theme_publication()

  ggsave(file.path(fig_dir, "community_size_spectrum.png"),
         p_spectrum, width = 10, height = 6, dpi = 300)
}

# ============================================================================
# Taxonomic Size Spectra
# ============================================================================

cat("Analyzing taxonomic-specific size spectra...\n")

# Calculate size spectra by taxonomic group
taxonomic_spectra <- cafi_biomass %>%
  filter(!is.na(size_mm)) %>%
  group_by(type) %>%
  filter(n_distinct(size_mm) > 1) %>%  # Need at least 2 unique values for binning
  mutate(size_class = cut(size_mm, breaks = 10)) %>%
  group_by(type, size_class) %>%
  summarise(
    abundance = n(),
    mean_size = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(size_class))

# Fit models for each group
taxonomic_slopes <- cafi_biomass %>%
  filter(!is.na(size_mm)) %>%
  group_by(type) %>%
  filter(n() > 20) %>%
  do({
    df_with_rank <- mutate(., rank = row_number())
    mod <- lm(log10(size_mm) ~ log10(rank), data = df_with_rank)
    data.frame(
      slope = coef(mod)[2],
      intercept = coef(mod)[1],
      r2 = summary(mod)$r.squared,
      n = nrow(.)
    )
  }) %>%
  ungroup()

write_csv(taxonomic_slopes,
          file.path(TABLES_DIR, "taxonomic_size_spectrum_slopes.csv"))

# Plot taxonomic size distributions
p_tax_size <- cafi_biomass %>%
  ggplot(aes(x = size_mm, fill = type)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~type, scales = "free_y") +
  scale_x_log10() +
  scale_fill_taxon() +
  labs(
    title = "Size distributions by taxonomic group",
    x = "Size (mm, log scale)",
    y = "Count",
    fill = "Type"
  ) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "taxonomic_size_distributions.png"),
       p_tax_size, width = 12, height = 8, dpi = 300)

# ============================================================================
# Biomass-Abundance Relationships
# ============================================================================

cat("Analyzing biomass-abundance relationships...\n")

# Calculate B-N relationships
bn_data <- cafi_biomass %>%
  group_by(coral_id, type) %>%
  summarise(
    abundance = n(),
    total_biomass = sum(biomass_g, na.rm = TRUE),
    mean_ind_mass = mean(biomass_g, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(abundance > 0, total_biomass > 0)

# Fit B-N scaling
bn_model <- lm(log10(total_biomass) ~ log10(abundance), data = bn_data)
bn_slope <- coef(bn_model)[2]
bn_r2 <- summary(bn_model)$r.squared

# Energy equivalence test (slope should be ~1 if energy equivalence holds)
energy_equiv_test <- summary(bn_model)$coefficients[2, ]
is_energy_equiv <- abs(bn_slope - 1) < 2 * energy_equiv_test[2]  # Within 2 SE of 1

cat("  B-N scaling exponent:", round(bn_slope, 3), "\n")
cat("  Energy equivalence:", ifelse(is_energy_equiv, "Supported", "Not supported"), "\n\n")

# Plot B-N relationship
p_bn <- ggplot(bn_data, aes(x = abundance, y = total_biomass)) +
  geom_point(aes(color = type, shape = type), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks() +
  scale_color_taxon() +
  scale_shape_taxon() +
  labs(
    title = "Biomass-abundance scaling",
    subtitle = paste("Slope =", round(bn_slope, 2), ", R² =", round(bn_r2, 3)),
    x = "Abundance (N)",
    y = "Total biomass (g)",
    color = "Type",
    shape = "Type"
  ) +
  theme_publication()

ggsave(file.path(fig_dir, "biomass_abundance_scaling.png"),
       p_bn, width = 10, height = 6, dpi = 300)

# ============================================================================
# Size Structure by Habitat
# ============================================================================

cat("Analyzing size structure across habitats...\n")

# Size distributions by morphotype
size_by_habitat <- cafi_biomass %>%
  left_join(coral_master %>% dplyr::select(coral_id, morphotype, depth_m), by = "coral_id") %>%
  filter(!is.na(morphotype))

# Calculate size metrics by morphotype
size_metrics_habitat <- size_by_habitat %>%
  group_by(morphotype, type) %>%
  summarise(
    n_individuals = n(),
    mean_size = mean(size_mm, na.rm = TRUE),
    median_size = median(size_mm, na.rm = TRUE),
    size_cv = sd(size_mm, na.rm = TRUE) / mean(size_mm, na.rm = TRUE),
    skewness = moments::skewness(size_mm, na.rm = TRUE),
    total_biomass = sum(biomass_g, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(size_metrics_habitat,
          file.path(TABLES_DIR, "size_structure_by_morphotype.csv"))

# Violin plot of size distributions
p_size_habitat <- size_by_habitat %>%
  ggplot(aes(x = morphotype, y = size_mm, fill = morphotype)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.1, outlier.alpha = 0.3) +
  facet_wrap(~type, scales = "free_y") +
  scale_y_log10() +
  scale_fill_viridis_d() +
  labs(
    title = "Size structure across morphotypes",
    x = "Morphotype",
    y = "Size (mm, log scale)"
  ) +
  theme_publication() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "size_structure_by_morphotype.png"),
       p_size_habitat, width = 12, height = 8, dpi = 300)

# ============================================================================
# Allometric Scaling Relationships
# ============================================================================

cat("Analyzing allometric relationships...\n")

# Test for allometric scaling within species
# Select species with sufficient data
species_for_allometry <- cafi_biomass %>%
  group_by(species) %>%
  filter(n() >= 20) %>%
  pull(species) %>%
  unique()

allometry_results <- list()

for (sp in species_for_allometry[1:min(10, length(species_for_allometry))]) {
  sp_data <- cafi_biomass %>%
    filter(species == sp, !is.na(size_mm), size_mm > 0)

  # Need at least 10 observations with valid size AND variation in size
  sp_data <- sp_data %>%
    filter(!is.na(log_biomass), !is.na(log_size),
           is.finite(log_biomass), is.finite(log_size))

  if (nrow(sp_data) > 10 && sd(sp_data$log_size, na.rm = TRUE) > 0.01) {
    # Fit SMA regression (appropriate for allometry)
    # Use smatr package if available, otherwise use OLS as approximation
    tryCatch({
      if (requireNamespace("smatr", quietly = TRUE)) {
        sma_model <- smatr::sma(log_biomass ~ log_size, data = sp_data)

        allometry_results[[sp]] <- data.frame(
          species = sp,
          n = nrow(sp_data),
          slope = sma_model$coef[[1]][2, 1],
          slope_lower = sma_model$coef[[1]][2, 2],
          slope_upper = sma_model$coef[[1]][2, 3],
          r2 = sma_model$r2[[1]],
          p_value = sma_model$pval[[1]],
          isometric = sma_model$coef[[1]][2, 1] > 2.5 & sma_model$coef[[1]][2, 1] < 3.5
        )
      } else {
        # Fallback: OLS regression
        ols_model <- lm(log_biomass ~ log_size, data = sp_data)
        ols_summ <- summary(ols_model)

        allometry_results[[sp]] <- data.frame(
          species = sp,
          n = nrow(sp_data),
          slope = coef(ols_model)[2],
          slope_lower = confint(ols_model)[2, 1],
          slope_upper = confint(ols_model)[2, 2],
          r2 = ols_summ$r.squared,
          p_value = coef(ols_summ)[2, 4],
          isometric = coef(ols_model)[2] > 2.5 & coef(ols_model)[2] < 3.5
        )
      }
    }, error = function(e) {
      cat("  Skipping", sp, "- regression failed:", conditionMessage(e), "\n")
    })
  }
}

if (length(allometry_results) > 0) {
  allometry_summary <- bind_rows(allometry_results)
  write_csv(allometry_summary,
            file.path(TABLES_DIR, "species_allometric_scaling.csv"))

  # Plot allometric slopes
  p_allometry <- allometry_summary %>%
    ggplot(aes(x = reorder(species, slope), y = slope)) +
    geom_hline(yintercept = 3, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_errorbar(aes(ymin = slope_lower, ymax = slope_upper), width = 0.2) +
    geom_point(aes(color = isometric), size = 3) +
    coord_flip() +
    scale_color_manual(values = c("FALSE" = "#0072B2", "TRUE" = "#009E73")) +
    labs(
      title = "Allometric scaling exponents by species",
      subtitle = "Red line = isometric scaling (b = 3)",
      x = "Species",
      y = "Scaling exponent (b)",
      color = "Isometric"
    ) +
    theme_publication()

  ggsave(file.path(fig_dir, "species_allometric_exponents.png"),
         p_allometry, width = 10, height = 8, dpi = 300)
}

# ============================================================================
# Biomass Dominance Patterns
# ============================================================================

cat("Analyzing biomass dominance...\n")

# Calculate biomass dominance by species
biomass_dominance <- cafi_biomass %>%
  group_by(species, type) %>%
  summarise(
    total_biomass = sum(biomass_g, na.rm = TRUE),
    total_abundance = n(),
    mean_ind_biomass = mean(biomass_g, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_biomass)) %>%
  mutate(
    biomass_rank = row_number(),
    cumulative_biomass = cumsum(total_biomass),
    cumulative_proportion = cumulative_biomass / sum(total_biomass)
  )

# Identify biomass dominants (80% of total biomass)
biomass_dominants <- biomass_dominance %>%
  filter(cumulative_proportion <= 0.8)

write_csv(biomass_dominants,
          file.path(TABLES_DIR, "biomass_dominant_species.csv"))

# Lorenz curve for biomass inequality
p_lorenz <- biomass_dominance %>%
  mutate(
    species_proportion = row_number() / n()
  ) %>%
  ggplot(aes(x = species_proportion, y = cumulative_proportion)) +
  geom_line(linewidth = 1.5, color = "#0072B2") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Biomass Lorenz curve",
    subtitle = paste(nrow(biomass_dominants), "species account for 80% of biomass"),
    x = "Cumulative proportion of species",
    y = "Cumulative proportion of biomass"
  ) +
  coord_equal() +
  theme_publication()

ggsave(file.path(fig_dir, "biomass_lorenz_curve.png"),
       p_lorenz, width = 8, height = 8, dpi = 300)

# ============================================================================
# Size-Based Trophic Structure
# ============================================================================

cat("Inferring size-based trophic structure...\n")

# Assign trophic levels based on size and type (simplified)
trophic_assignments <- cafi_biomass %>%
  mutate(
    trophic_level = case_when(
      type == "snail" & size_mm < 5 ~ 2.0,  # Small grazers
      type == "snail" ~ 2.2,
      type == "shrimp" & size_mm < 10 ~ 2.5,  # Detritivores
      type == "shrimp" ~ 2.8,
      type == "crab" & size_mm < 15 ~ 2.8,  # Small omnivores
      type == "crab" ~ 3.2,
      type == "fish" & size_mm < 30 ~ 3.0,  # Small planktivores
      type == "fish" ~ 3.5,  # Larger predators
      TRUE ~ 2.5
    )
  )

# Trophic pyramid
trophic_pyramid <- trophic_assignments %>%
  mutate(trophic_bin = cut(trophic_level, breaks = c(1.5, 2.5, 3.0, 3.5, 4.0),
                           labels = c("Primary", "Secondary", "Tertiary", "Quaternary"))) %>%
  group_by(trophic_bin) %>%
  summarise(
    total_biomass = sum(biomass_g, na.rm = TRUE),
    total_abundance = n(),
    mean_size = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(trophic_bin))

# Plot trophic pyramid
p_trophic <- ggplot(trophic_pyramid, aes(x = trophic_bin, y = total_biomass)) +
  geom_col(aes(fill = trophic_bin), width = 0.8) +
  scale_fill_viridis_d() +
  coord_flip() +
  labs(
    title = "Size-based trophic pyramid",
    x = "Trophic level",
    y = "Total biomass (g)",
    fill = "Trophic level"
  ) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "trophic_pyramid_biomass.png"),
       p_trophic, width = 10, height = 6, dpi = 300)

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Size and Biomass Scaling Summary\n")
cat("========================================\n\n")

cat("Biomass Statistics:\n")
cat("  - Total CAFI biomass:", round(sum(coral_biomass$total_biomass_g), 2), "g\n")
cat("  - Mean biomass per coral:", round(mean(coral_biomass$total_biomass_g), 3), "g\n")
cat("  - Biomass per capita:", round(mean(coral_biomass$biomass_per_capita), 4), "g\n\n")

if (exists("spectrum_slope")) {
  cat("Size Spectrum:\n")
  cat("  - Community slope:", round(spectrum_slope, 2), "\n")
  cat("  - Spectrum R²:", round(spectrum_r2, 3), "\n\n")
}

cat("Biomass-Abundance Scaling:\n")
cat("  - B-N exponent:", round(bn_slope, 3), "\n")
cat("  - Energy equivalence:", ifelse(is_energy_equiv, "Supported", "Not supported"), "\n\n")

cat("Biomass Dominance:\n")
cat("  - Species for 80% biomass:", nrow(biomass_dominants), "of",
    n_distinct(cafi_biomass$species), "\n")
cat("  - Top biomass species:", biomass_dominance$species[1],
    "(", round(biomass_dominance$total_biomass[1], 2), "g)\n\n")

# ============================================================================
# STATISTICAL RESULTS SUMMARY
# ============================================================================

cat("Compiling statistical results...\n")

stats_results <- init_results_df()

# 1. Size spectrum slope
if(exists("spectrum_model")) {
  spec_summary <- summary(spectrum_model)
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H_spectrum",
      question = "What is the slope of the community size spectrum?",
      test_name = "Linear regression (log abundance ~ log size)",
      test_statistic = spec_summary$coefficients["log_size", "t value"],
      df = spectrum_model$df.residual,
      p_value = spec_summary$coefficients["log_size", "Pr(>|t|)"],
      effect_size = spectrum_r2,
      effect_type = "R²",
      n = nrow(size_spectra),
      notes = sprintf("Slope=%.3f", spectrum_slope)
    )
  )
}

# 2. Biomass-abundance scaling
if(exists("bn_model")) {
  bn_summary <- summary(bn_model)
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H_BN",
      question = "Does biomass scale with abundance (energy equivalence)?",
      test_name = "Linear regression (log B ~ log N)",
      test_statistic = bn_summary$coefficients["log10(abundance)", "t value"],
      df = bn_model$df.residual,
      p_value = bn_summary$coefficients["log10(abundance)", "Pr(>|t|)"],
      effect_size = bn_r2,
      effect_type = "R²",
      n = nrow(bn_data),
      notes = sprintf("Slope=%.3f, Energy equiv=%s", bn_slope,
                      ifelse(is_energy_equiv, "Yes", "No"))
    )
  )
}

# 3. Allometric scaling results
if(exists("allometry_summary") && nrow(allometry_summary) > 0) {
  for(i in 1:min(5, nrow(allometry_summary))) {
    row <- allometry_summary[i, ]
    stats_results <- bind_rows(stats_results,
      create_result_row(
        hypothesis = paste0("H_allom_", i),
        question = paste("Is", row$species, "allometry isometric (b≈3)?"),
        test_name = "SMA/OLS regression",
        test_statistic = row$slope,
        df = row$n - 2,
        p_value = row$p_value,
        effect_size = row$r2,
        effect_type = "R²",
        n = row$n,
        notes = sprintf("Slope=%.3f [%.3f-%.3f], Isometric=%s",
                        row$slope, row$slope_lower, row$slope_upper,
                        ifelse(row$isometric, "Yes", "No"))
      )
    )
  }
}

# 4. Kruskal-Wallis for size by taxonomic group
if(nrow(cafi_biomass) > 10) {
  kw_size_type <- kruskal.test(size_mm ~ type, data = cafi_biomass %>% filter(!is.na(size_mm)))
  n_total <- nrow(cafi_biomass %>% filter(!is.na(size_mm)))
  k <- n_distinct(cafi_biomass$type)
  epsilon_sq <- kw_size_type$statistic / (n_total - 1)

  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H_size_type",
      question = "Does body size differ among taxonomic groups?",
      test_name = "Kruskal-Wallis",
      test_statistic = kw_size_type$statistic,
      df = k - 1,
      p_value = kw_size_type$p.value,
      effect_size = as.numeric(epsilon_sq),
      effect_type = "ε²",
      n = n_total,
      notes = "Size differences across crab, shrimp, fish, snail"
    )
  )
}

# 5. Kruskal-Wallis for biomass by taxonomic group
if(nrow(cafi_biomass) > 10 && sum(!is.na(cafi_biomass$biomass_g)) > 10) {
  kw_biomass_type <- kruskal.test(biomass_g ~ type, data = cafi_biomass %>% filter(!is.na(biomass_g)))
  n_total <- nrow(cafi_biomass %>% filter(!is.na(biomass_g)))
  k <- n_distinct(cafi_biomass$type)
  epsilon_sq <- kw_biomass_type$statistic / (n_total - 1)

  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H_biomass_type",
      question = "Does individual biomass differ among taxonomic groups?",
      test_name = "Kruskal-Wallis",
      test_statistic = kw_biomass_type$statistic,
      df = k - 1,
      p_value = kw_biomass_type$p.value,
      effect_size = as.numeric(epsilon_sq),
      effect_type = "ε²",
      n = n_total,
      notes = "Biomass differences across taxonomic groups"
    )
  )
}

# 6. Biomass dominance (descriptive)
if(exists("biomass_dominants") && nrow(biomass_dominants) > 0) {
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "Desc_dominance",
      question = "How many species account for 80% of biomass?",
      test_name = "Cumulative biomass analysis",
      test_statistic = nrow(biomass_dominants),
      df = NA,
      p_value = NA,
      effect_size = nrow(biomass_dominants) / n_distinct(cafi_biomass$species),
      effect_type = "Proportion of species",
      n = n_distinct(cafi_biomass$species),
      notes = sprintf("Top species: %s (%.2f g)", biomass_dominance$species[1], biomass_dominance$total_biomass[1])
    )
  )
}

# Save statistical results
save_stats_summary(stats_results, "11_size_biomass_scaling", "Size-Biomass Scaling Analysis")

cat("✅ Size and biomass scaling analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", TABLES_DIR, "\n")