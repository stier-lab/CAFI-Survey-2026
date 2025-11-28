#!/usr/bin/env Rscript
# ============================================================================
# 08_size_biomass_scaling.R - Size structure and biomass scaling analysis
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
# ============================================================================

cat("\n========================================\n")
cat("Size Structure and Biomass Scaling Analysis\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/00_load_libraries.R"))
# smatr is optional - load if available
if (requireNamespace("smatr", quietly = TRUE)) {
  library(smatr)
} else {
  cat("Note: smatr package not installed. Some scaling analyses will be skipped.\n")
}
library(segmented)

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "size_biomass_scaling")
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
          file.path(SURVEY_TABLES, "coral_cafi_biomass_estimates.csv"))

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
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(
      title = "Community Size Spectrum",
      subtitle = paste("Slope =", round(spectrum_slope, 2), ", R² =", round(spectrum_r2, 3)),
      x = "Size (mm)",
      y = "Normalized Abundance (N/mm)"
    )

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
          file.path(SURVEY_TABLES, "taxonomic_size_spectrum_slopes.csv"))

# Plot taxonomic size distributions
p_tax_size <- cafi_biomass %>%
  ggplot(aes(x = size_mm, fill = type)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~type, scales = "free_y") +
  scale_x_log10() +
  scale_fill_viridis_d() +
  labs(
    title = "Size Distributions by Taxonomic Group",
    x = "Size (mm, log scale)",
    y = "Count",
    fill = "Type"
  ) +
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
  geom_point(aes(color = type), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks() +
  scale_color_viridis_d() +
  labs(
    title = "Biomass-Abundance Scaling",
    subtitle = paste("Slope =", round(bn_slope, 2), ", R² =", round(bn_r2, 3)),
    x = "Abundance (N)",
    y = "Total Biomass (g)",
    color = "Type"
  )

ggsave(file.path(fig_dir, "biomass_abundance_scaling.png"),
       p_bn, width = 10, height = 6, dpi = 300)

# ============================================================================
# Size Structure by Habitat
# ============================================================================

cat("Analyzing size structure across habitats...\n")

# Size distributions by morphotype
size_by_habitat <- cafi_biomass %>%
  left_join(metadata %>% dplyr::select(coral_id, morphotype, depth_m), by = "coral_id") %>%
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
          file.path(SURVEY_TABLES, "size_structure_by_morphotype.csv"))

# Violin plot of size distributions
p_size_habitat <- size_by_habitat %>%
  ggplot(aes(x = morphotype, y = size_mm, fill = morphotype)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.1, outlier.alpha = 0.3) +
  facet_wrap(~type, scales = "free_y") +
  scale_y_log10() +
  scale_fill_viridis_d() +
  labs(
    title = "Size Structure Across Morphotypes",
    x = "Morphotype",
    y = "Size (mm, log scale)"
  ) +
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

  if (nrow(sp_data) > 10) {
    # Fit SMA regression (appropriate for allometry)
    sma_model <- sma(log_biomass ~ log_size, data = sp_data)

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
  }
}

if (length(allometry_results) > 0) {
  allometry_summary <- bind_rows(allometry_results)
  write_csv(allometry_summary,
            file.path(SURVEY_TABLES, "species_allometric_scaling.csv"))

  # Plot allometric slopes
  p_allometry <- allometry_summary %>%
    ggplot(aes(x = reorder(species, slope), y = slope)) +
    geom_hline(yintercept = 3, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_errorbar(aes(ymin = slope_lower, ymax = slope_upper), width = 0.2) +
    geom_point(aes(color = isometric), size = 3) +
    coord_flip() +
    scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "green")) +
    labs(
      title = "Allometric Scaling Exponents by Species",
      subtitle = "Red line = isometric scaling (b = 3)",
      x = "Species",
      y = "Scaling Exponent (b)",
      color = "Isometric"
    )

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
          file.path(SURVEY_TABLES, "biomass_dominant_species.csv"))

# Lorenz curve for biomass inequality
p_lorenz <- biomass_dominance %>%
  mutate(
    species_proportion = row_number() / n()
  ) %>%
  ggplot(aes(x = species_proportion, y = cumulative_proportion)) +
  geom_line(size = 1.5, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Biomass Lorenz Curve",
    subtitle = paste(nrow(biomass_dominants), "species account for 80% of biomass"),
    x = "Cumulative Proportion of Species",
    y = "Cumulative Proportion of Biomass"
  ) +
  coord_equal()

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
    title = "Size-Based Trophic Pyramid",
    x = "Trophic Level",
    y = "Total Biomass (g)",
    fill = "Trophic Level"
  ) +
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

cat("✅ Size and biomass scaling analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", SURVEY_TABLES, "\n")