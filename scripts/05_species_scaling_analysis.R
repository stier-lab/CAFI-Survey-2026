# ============================================================================
# 05_species_scaling_analysis.R - Comprehensive CAFI-Coral Scaling Analysis
# ============================================================================
#
# PURPOSE: Test how CAFI abundance scales with coral colony volume
#
# THEORETICAL FRAMEWORK:
# ----------------------
# The relationship between habitat (coral) size and inhabitant (CAFI) abundance
# follows a power law: N = a × V^β, where:
#   - N = abundance (total CAFI, species count, etc.)
#   - V = coral volume (cm³)
#   - a = intercept (baseline density)
#   - β = scaling exponent (the key parameter of interest)
#
# In log-log space: log(N) = log(a) + β × log(V)
#
# THREE COMPETING HYPOTHESES:
# ---------------------------
#
# H1: FIELD OF DREAMS (β = 1.0)
#     "If you build it, they will come"
#     - Abundance scales PROPORTIONALLY with volume
#     - Per-capita density constant across coral sizes
#     - Mechanism: Passive habitat filling
#
# H2: PROPAGULE REDISTRIBUTION (β < 1.0)
#     Following Stier & Osenberg (2010) Ecology
#     - Abundance scales SUBLINEARLY with volume
#     - Per-capita density DECREASES in larger corals (dilution effect)
#     - Mechanism: Propagule supply limited; larger corals "dilute" settlers
#     - Also called "dilution effect" or "propagule limitation"
#
# H3: SUPER-LINEAR SCALING (β > 1.0)
#     - Abundance scales FASTER than volume
#     - Per-capita density INCREASES in larger corals
#     - Mechanism: Larger corals disproportionately attractive
#     - Could indicate aggregation, facilitation, or preferred habitat
#
# ANALYSES IN THIS SCRIPT:
# ------------------------
# PART 1: Total CAFI Abundance (community-level)
# PART 2: Species Richness
# PART 3: Shannon Diversity
# PART 4: Individual Species (abundant taxa with sufficient prevalence)
# PART 5: Taxonomic Groups (crabs, shrimps, fishes, snails)
# PART 6: Taxonomic Groups (Trapezia, Gastropods, Fish, etc.)
#
# MODEL SPECIFICATION:
#   log(Y) ~ β × log(V) + site
#   - Negative binomial GLM for counts (handles overdispersion)
#   - Site as fixed effect (accounts for among-site variation)
#   - log10 transformation for interpretability
#
# KEY TEST:
#   H0: β = 1 (Field of Dreams)
#   HA: β ≠ 1
#   Using z-test: z = (β - 1) / SE(β)
#
# REFERENCES:
#   Stier AC, Osenberg CW (2010) Propagule redirection: habitat availability
#     reduces colonization and increases recruitment in reef fishes.
#     Ecology 91:2884-2892
#
#   Stier AC, Osenberg CW (2024) Field of dreams or propagule redistribution?
#     Current Biology 34:R1186-R1189
#
# OUTPUTS:
#   Figures: output/figures/05_scaling/
#   Tables: output/tables/
#   Report: output/reports/SCALING_ANALYSIS.html
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-10
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI-CORAL SCALING ANALYSIS\n")
cat("    Field of Dreams vs. Propagule Redistribution\n")
cat("============================================================\n\n")

# Load setup and data
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

# Use script-specific output directory
FIG_DIR <- PATHS$fig_05_scaling

# Load CAFI data
cafi_clean <- readRDS(file.path(PATHS$objects, "cafi_clean.rds"))

# ============================================================================
# HELPER FUNCTION: Fit Scaling Model with Bootstrap CI
# ============================================================================
# Following Stier et al. (2024) Ecology Letters methodology:
# Use bootstrap resampling to generate 95% CI around the scaling exponent
# If 95% CI includes 1 → Field of Dreams
# If 95% CI excludes 1 and β < 1 → Propagule Redirection
# ============================================================================

#' Bootstrap function to estimate scaling exponent
#' @param data Data with abundance and volume
#' @param indices Bootstrap sample indices
#' @return Scaling exponent beta
boot_scaling <- function(data, indices) {
  d <- data[indices, ]
  tryCatch({
    # Include site as fixed effect to match main model specification
    if ("site" %in% names(d) && length(unique(d$site)) > 1) {
      model <- MASS::glm.nb(abundance ~ log10(volume) + site, data = d)
    } else {
      model <- MASS::glm.nb(abundance ~ log10(volume), data = d)
    }
    coef(model)["log10(volume)"]
  }, error = function(e) NA)
}

#' Fit a scaling model and test against Field of Dreams (β = 1)
#'
#' @param data Data frame with 'abundance' and 'volume' columns
#' @param response_name Character string for labeling
#' @param min_nonzero Minimum non-zero observations required (default 15)
#' @param n_boot Number of bootstrap replicates (default 2000)
#' @return List with model results and interpretation
fit_scaling_model <- function(data, response_name, min_nonzero = 15, n_boot = 2000) {

  # Calculate summary stats
  n_total <- nrow(data)
  n_nonzero <- sum(data$abundance > 0)
  n_zero <- sum(data$abundance == 0)
  total_abundance <- sum(data$abundance)

  # Check for sufficient data
  if (n_nonzero < min_nonzero) {
    return(list(
      response = response_name,
      n_corals = n_total,
      n_nonzero = n_nonzero,
      n_zero = n_zero,
      total_abundance = total_abundance,
      beta = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      boot_ci_lower = NA_real_,
      boot_ci_upper = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      z_vs_1 = NA_real_,
      p_vs_1 = NA_real_,
      pseudo_r2 = NA_real_,
      interpretation = "Insufficient data",
      model = NULL,
      converged = FALSE
    ))
  }

  # Fit negative binomial GLM
  tryCatch({
    model <- MASS::glm.nb(abundance ~ log10(volume) + site, data = data)

    # Extract coefficients
    coefs <- summary(model)$coefficients
    beta <- coefs["log10(volume)", "Estimate"]
    se <- coefs["log10(volume)", "Std. Error"]
    z_val <- coefs["log10(volume)", "z value"]
    p_val <- coefs["log10(volume)", "Pr(>|z|)"]

    # Profile-based 95% CI
    ci <- confint(model, "log10(volume)", level = 0.95)

    # Bootstrap 95% CI (following Stier et al. 2024 methodology)
    boot_ci_lower <- NA_real_
    boot_ci_upper <- NA_real_

    if (requireNamespace("boot", quietly = TRUE) && n_nonzero >= 20) {
      suppressWarnings({
        # Stratify by site to maintain site proportions in each bootstrap replicate
        boot_result <- boot::boot(data, boot_scaling, R = n_boot,
                                  strata = factor(data$site))
        valid_betas <- boot_result$t[!is.na(boot_result$t)]

        if (length(valid_betas) >= n_boot * 0.8) {  # Require 80% successful fits
          boot_ci <- quantile(valid_betas, probs = c(0.025, 0.975))
          boot_ci_lower <- boot_ci[1]
          boot_ci_upper <- boot_ci[2]
        }
      })
    }

    # Use bootstrap CI if available, otherwise profile CI
    test_ci_lower <- if (!is.na(boot_ci_lower)) boot_ci_lower else ci[1]
    test_ci_upper <- if (!is.na(boot_ci_upper)) boot_ci_upper else ci[2]

    # Test vs β = 1 using z-test
    z_vs_1 <- (beta - 1) / se
    p_vs_1 <- 2 * pnorm(-abs(z_vs_1))

    # Pseudo-R² (McFadden's)
    pseudo_r2 <- calc_pseudo_r2(model)

    # Interpret based on CI (whether 95% CI includes 1)
    # Following Stier et al. (2024) approach
    if (test_ci_upper < 1) {
      interpretation <- "Redistribution (β < 1)"
    } else if (test_ci_lower > 1) {
      interpretation <- "Super-linear (β > 1)"
    } else {
      interpretation <- "Field of Dreams (β ≈ 1)"
    }

    list(
      response = response_name,
      n_corals = n_total,
      n_nonzero = n_nonzero,
      n_zero = n_zero,
      total_abundance = total_abundance,
      beta = beta,
      se = se,
      ci_lower = ci[1],
      ci_upper = ci[2],
      boot_ci_lower = boot_ci_lower,
      boot_ci_upper = boot_ci_upper,
      z_value = z_val,
      p_value = p_val,
      z_vs_1 = z_vs_1,
      p_vs_1 = p_vs_1,
      pseudo_r2 = pseudo_r2,
      interpretation = interpretation,
      model = model,
      converged = TRUE
    )

  }, error = function(e) {
    list(
      response = response_name,
      n_corals = n_total,
      n_nonzero = n_nonzero,
      n_zero = n_zero,
      total_abundance = total_abundance,
      beta = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      boot_ci_lower = NA_real_,
      boot_ci_upper = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      z_vs_1 = NA_real_,
      p_vs_1 = NA_real_,
      interpretation = paste("Model error:", conditionMessage(e)),
      model = NULL,
      converged = FALSE
    )
  })
}

#' Convert scaling result list to single-row tibble
result_to_row <- function(result) {
  tibble(
    response = result$response,
    n_corals = result$n_corals,
    n_nonzero = result$n_nonzero,
    n_zero = result$n_zero,
    total_abundance = result$total_abundance,
    beta = result$beta,
    se = result$se,
    ci_lower = result$ci_lower,
    ci_upper = result$ci_upper,
    boot_ci_lower = if (!is.null(result$boot_ci_lower)) result$boot_ci_lower else NA_real_,
    boot_ci_upper = if (!is.null(result$boot_ci_upper)) result$boot_ci_upper else NA_real_,
    z_value = result$z_value,
    p_value = result$p_value,
    z_vs_1 = result$z_vs_1,
    p_vs_1 = result$p_vs_1,
    interpretation = result$interpretation
  )
}

# ============================================================================
# PREPARE BASE DATA
# ============================================================================

cat("Preparing data for scaling analysis...\n\n")

# Base coral data with volume
coral_data <- coral_master %>%
  filter(!is.na(volume), volume > 0) %>%
  dplyr::select(coral_id, site, volume, log_volume, total_cafi, otu_richness, shannon,
                n_trapezia, n_resident_fish, n_corallivore, n_other_crab, n_shrimp, n_other,
                n_crabs, n_shrimps, n_fish, n_snails)

cat("Sample size:", nrow(coral_data), "corals with volume data\n")
cat("Sites:", paste(unique(coral_data$site), collapse = ", "), "\n")
cat("Volume range:", round(min(coral_data$volume)), "-", round(max(coral_data$volume)), "cm³\n")
cat("Total CAFI:", sum(coral_data$total_cafi), "individuals\n\n")

# ############################################################################
#                    PART 1: TOTAL CAFI ABUNDANCE
# ############################################################################

cat("############################################################\n")
cat("PART 1: TOTAL CAFI ABUNDANCE SCALING\n")
cat("############################################################\n\n")

# Prepare data
abundance_data <- coral_data %>%
  dplyr::select(coral_id, site, volume, abundance = total_cafi)

# Fit model
total_result <- fit_scaling_model(abundance_data, "Total CAFI Abundance")

cat("Model: Total CAFI ~ log10(Volume) + site (Negative Binomial GLM)\n\n")
cat("Results:\n")
cat("  Scaling exponent β =", round(total_result$beta, 3), "\n")
cat("  Standard error =", round(total_result$se, 3), "\n")
cat("  Profile 95% CI: [", round(total_result$ci_lower, 3), ",", round(total_result$ci_upper, 3), "]\n")
if (!is.na(total_result$boot_ci_lower)) {
  cat("  Bootstrap 95% CI: [", round(total_result$boot_ci_lower, 3), ",", round(total_result$boot_ci_upper, 3), "]\n")
}
cat("  z-value =", round(total_result$z_value, 2), ", p =", format.pval(total_result$p_value, 3), "\n\n")

cat("Test vs Field of Dreams (H0: β = 1):\n")
cat("  z =", round(total_result$z_vs_1, 2), ", p =", format.pval(total_result$p_vs_1, 3), "\n")
if (!is.na(total_result$boot_ci_lower)) {
  cat("  Bootstrap CI includes 1:", ifelse(total_result$boot_ci_lower <= 1 & total_result$boot_ci_upper >= 1, "YES", "NO"), "\n")
}
cat("  Interpretation:", total_result$interpretation, "\n\n")

# Store for combined results
all_results <- list(total_abundance = total_result)

# ############################################################################
#                    PART 2: SPECIES RICHNESS SCALING
# ############################################################################

cat("############################################################\n")
cat("PART 2: SPECIES RICHNESS SCALING\n")
cat("############################################################\n\n")

richness_data <- coral_data %>%
  dplyr::select(coral_id, site, volume, abundance = otu_richness)

richness_result <- fit_scaling_model(richness_data, "Species Richness")

cat("Model: OTU Richness ~ log10(Volume) + site (Negative Binomial GLM)\n\n")
cat("Results:\n")
cat("  Scaling exponent β =", round(richness_result$beta, 3), "\n")
cat("  95% CI: [", round(richness_result$ci_lower, 3), ",", round(richness_result$ci_upper, 3), "]\n")
cat("  Test vs β = 1: z =", round(richness_result$z_vs_1, 2), ", p =", format.pval(richness_result$p_vs_1, 3), "\n")
cat("  Interpretation:", richness_result$interpretation, "\n\n")

all_results$species_richness <- richness_result

# ---- Model Diagnostics (DHARMa simulated residuals) ----
cat("--- Model Diagnostics ---\n")
if (requireNamespace("DHARMa", quietly = TRUE) && !is.null(total_result$model)) {
  # Total abundance (NB GLM)
  dharma_total <- DHARMa::simulateResiduals(total_result$model, n = 1000, plot = FALSE)
  dharma_tests <- DHARMa::testResiduals(dharma_total, plot = FALSE)
  cat("  TOTAL ABUNDANCE (NB GLM):\n")
  cat("    Uniformity (KS): p =", format.pval(dharma_tests$uniformity$p.value, 3), "\n")
  cat("    Dispersion: p =", format.pval(dharma_tests$dispersion$p.value, 3), "\n")
  cat("    Zero-inflation: p =", format.pval(
    DHARMa::testZeroInflation(dharma_total, plot = FALSE)$p.value, 3), "\n")

  # Species richness (NB GLM)
  if (!is.null(richness_result$model)) {
    dharma_rich <- DHARMa::simulateResiduals(richness_result$model, n = 1000, plot = FALSE)
    dharma_tests_rich <- DHARMa::testResiduals(dharma_rich, plot = FALSE)
    cat("  SPECIES RICHNESS (NB GLM):\n")
    cat("    Uniformity (KS): p =", format.pval(dharma_tests_rich$uniformity$p.value, 3), "\n")
    cat("    Dispersion: p =", format.pval(dharma_tests_rich$dispersion$p.value, 3), "\n")
    cat("    Zero-inflation: p =", format.pval(
      DHARMa::testZeroInflation(dharma_rich, plot = FALSE)$p.value, 3), "\n")
  }

  # Save diagnostic plots
  fig_dir <- file.path(PATHS$figures, "05_scaling")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  png(file.path(fig_dir, "dharma_diagnostics.png"), width = 10, height = 5, units = "in", res = 200)
  par(mfrow = c(1, 2))
  plot(dharma_total, main = "Total CAFI (NB)")
  if (!is.null(richness_result$model)) plot(dharma_rich, main = "Richness (NB)")
  dev.off()
  cat("  Saved: dharma_diagnostics.png\n")
} else {
  cat("  DHARMa not available or model not retained - skipping\n")
}
cat("\n")

# ############################################################################
#                    PART 3: SHANNON DIVERSITY SCALING
# ############################################################################

cat("############################################################\n")
cat("PART 3: SHANNON DIVERSITY SCALING\n")
cat("############################################################\n\n")

# Shannon requires linear model (continuous, not count)
shannon_data <- coral_data %>%
  filter(shannon > 0) %>%
  dplyr::select(coral_id, site, volume, shannon)

cat("Model: Shannon H' ~ log10(Volume) + site (Linear Model)\n\n")

if (nrow(shannon_data) >= 15) {
  shannon_model <- lm(shannon ~ log10(volume) + site, data = shannon_data)
  shannon_summary <- summary(shannon_model)

  shannon_beta <- coef(shannon_model)["log10(volume)"]
  shannon_se <- shannon_summary$coefficients["log10(volume)", "Std. Error"]
  shannon_ci <- confint(shannon_model, "log10(volume)", level = 0.95)
  shannon_t <- shannon_summary$coefficients["log10(volume)", "t value"]
  shannon_p <- shannon_summary$coefficients["log10(volume)", "Pr(>|t|)"]

  # Test vs β = 0 (does Shannon diversity scale with volume at all?)
  # Note: Testing vs β=1 is meaningless for Shannon H' because it is already

  # log-transformed (H' = -Σ pi*ln(pi)). The relevant null is β=0 (no scaling).
  shannon_t_vs_0 <- shannon_t  # standard t-test from lm() is already vs 0
  shannon_p_vs_0 <- shannon_p

  cat("Results:\n")
  cat("  Scaling exponent β =", round(shannon_beta, 3), "\n")
  cat("  95% CI: [", round(shannon_ci[1], 3), ",", round(shannon_ci[2], 3), "]\n")
  cat("  R² =", round(shannon_summary$r.squared, 3), "\n")
  cat("  Test vs β = 0: t =", round(shannon_t_vs_0, 2), ", p =", format.pval(shannon_p_vs_0, 3), "\n\n")

  shannon_result <- list(
    response = "Shannon Diversity",
    n_corals = nrow(shannon_data),
    n_nonzero = nrow(shannon_data),
    n_zero = 0,
    total_abundance = NA,
    beta = shannon_beta,
    se = shannon_se,
    ci_lower = shannon_ci[1],
    ci_upper = shannon_ci[2],
    z_value = shannon_t,
    p_value = shannon_p,
    z_vs_1 = shannon_t_vs_0,
    p_vs_1 = shannon_p_vs_0,
    interpretation = ifelse(shannon_p_vs_0 < 0.05 && shannon_beta > 0,
                           "Positive scaling (β > 0)",
                           ifelse(shannon_p_vs_0 >= 0.05, "No significant scaling",
                                  "Negative scaling (β < 0)")),
    model = shannon_model,
    converged = TRUE
  )
} else {
  shannon_result <- list(
    response = "Shannon Diversity",
    n_corals = nrow(shannon_data), n_nonzero = NA, n_zero = NA, total_abundance = NA,
    beta = NA, se = NA, ci_lower = NA, ci_upper = NA, z_value = NA, p_value = NA,
    z_vs_1 = NA, p_vs_1 = NA, interpretation = "Insufficient data", model = NULL, converged = FALSE
  )
  cat("Insufficient data for Shannon scaling analysis\n\n")
}

all_results$shannon_diversity <- shannon_result

# ############################################################################
#                    PART 4: TAXONOMIC GROUP SCALING
# ############################################################################

cat("############################################################\n")
cat("PART 4: TAXONOMIC GROUP SCALING\n")
cat("############################################################\n\n")

cat("Testing scaling for major taxonomic groups...\n\n")

taxonomic_groups <- list(
  "Crabs" = "n_crabs",
  "Shrimps" = "n_shrimps",
  "Fishes" = "n_fish",
  "Snails" = "n_snails"
)

taxonomic_results <- list()

for (group_name in names(taxonomic_groups)) {
  col_name <- taxonomic_groups[[group_name]]

  group_data <- coral_data %>%
    dplyr::select(coral_id, site, volume, abundance = all_of(col_name))

  result <- fit_scaling_model(group_data, group_name)
  taxonomic_results[[group_name]] <- result

  if (result$converged && !is.na(result$beta)) {
    cat(sprintf("%-10s β = %6.3f [%5.2f, %5.2f], p vs 1: %s - %s\n",
                group_name, result$beta, result$ci_lower, result$ci_upper,
                format.pval(result$p_vs_1, 2), result$interpretation))
  } else {
    cat(sprintf("%-10s %s\n", group_name, result$interpretation))
  }
}

cat("\n")
all_results$taxonomic_groups <- taxonomic_results

# ############################################################################
#                    PART 5: FUNCTIONAL GROUP SCALING
# ############################################################################

cat("############################################################\n")
cat("PART 5: FUNCTIONAL GROUP SCALING\n")
cat("############################################################\n\n")

cat("Testing scaling for functional groups...\n\n")

functional_groups <- list(
  "Trapezia crabs" = "n_trapezia",
  "Fish" = "n_resident_fish",
  "Gastropods" = "n_corallivore",  # 73% Galeropsis monodonta (Coralliophilinae)
  "Other crabs" = "n_other_crab",
  "Shrimps" = "n_shrimp",
  "Other invertebrates" = "n_other"
)

functional_results <- list()

for (group_name in names(functional_groups)) {
  col_name <- functional_groups[[group_name]]

  group_data <- coral_data %>%
    dplyr::select(coral_id, site, volume, abundance = all_of(col_name))

  result <- fit_scaling_model(group_data, group_name)
  functional_results[[group_name]] <- result

  if (result$converged && !is.na(result$beta)) {
    cat(sprintf("%-25s β = %6.3f [%5.2f, %5.2f], p vs 1: %s - %s\n",
                group_name, result$beta, result$ci_lower, result$ci_upper,
                format.pval(result$p_vs_1, 2), result$interpretation))
  } else {
    cat(sprintf("%-25s %s\n", group_name, result$interpretation))
  }
}

cat("\n")
all_results$functional_groups <- functional_results

# ############################################################################
#                    PART 6: INDIVIDUAL SPECIES SCALING
# ############################################################################

cat("############################################################\n")
cat("PART 6: INDIVIDUAL SPECIES SCALING\n")
cat("############################################################\n\n")

# Calculate species prevalence and abundance
species_summary <- cafi_clean %>%
  group_by(otu) %>%
  summarise(
    total_individuals = n(),
    n_corals_present = n_distinct(coral_id),
    prevalence = n_distinct(coral_id) / n_distinct(cafi_clean$coral_id) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(total_individuals))

# Get functional group info
species_functional <- cafi_clean %>%
  dplyr::select(otu, functional_group, type) %>%
  distinct()

species_summary <- species_summary %>%
  left_join(species_functional, by = "otu")

# Selection criteria: >= 30 individuals AND >= 15% prevalence
key_species <- species_summary %>%
  filter(total_individuals >= 30, prevalence >= 15) %>%
  arrange(desc(total_individuals))

cat("Selection criteria:\n")
cat("  - Minimum 30 total individuals\n")
cat("  - Minimum 15% prevalence (present on >= 17 corals)\n")
cat("  - Result:", nrow(key_species), "species selected\n\n")

# Create species × coral abundance matrix
species_by_coral <- cafi_clean %>%
  filter(otu %in% key_species$otu) %>%
  group_by(coral_id, otu) %>%
  summarise(abundance = n(), .groups = "drop") %>%
  complete(coral_id = unique(cafi_clean$coral_id), otu, fill = list(abundance = 0))

# Merge with coral data
species_scaling_data <- species_by_coral %>%
  left_join(coral_data %>% dplyr::select(coral_id, site, volume), by = "coral_id") %>%
  filter(!is.na(volume))

cat("Fitting models for", nrow(key_species), "species...\n\n")

species_results <- map(key_species$otu, function(sp) {
  sp_data <- species_scaling_data %>%
    filter(otu == sp) %>%
    dplyr::select(coral_id, site, volume, abundance)

  result <- fit_scaling_model(sp_data, sp)

  if (result$converged && !is.na(result$beta)) {
    cat(sprintf("  %-30s β = %6.3f, p vs 1: %s - %s\n",
                sp, result$beta, format.pval(result$p_vs_1, 2), result$interpretation))
  } else {
    cat(sprintf("  %-30s %s\n", sp, result$interpretation))
  }

  result
})

names(species_results) <- key_species$otu
all_results$individual_species <- species_results

# ############################################################################
#                    PART 6B: KEY SPECIES FROM EXPERIMENTAL PAPER
# ############################################################################
# Test scaling patterns for species highlighted in Stier et al. experimental work
# These species showed significant effects on coral condition:
#   POSITIVE effects: Caracanthus maculatus, Alpheus lottini (and Alpheus spp.)
#   NEGATIVE effects: Cymo quadrilobatus, Luniella pugil
#
# Key question: Do these functionally important species show heterogeneous
# scaling responses (as predicted by "Field of Dreams" heterogeneity)?
# ############################################################################

cat("\n############################################################\n")
cat("PART 6B: KEY SPECIES FROM EXPERIMENTAL PAPER\n")
cat("############################################################\n\n")

cat("Testing scaling patterns for species with documented effects on coral condition\n")
cat("from Stier et al. experimental studies.\n\n")

# Define key species from experimental paper
key_exp_species <- list(
  list(
    name = "Caracanthus maculatus",
    pattern = "Caracanthus",
    effect = "positive",
    description = "Coral croucher (guards against corallivores)"
  ),
  list(
    name = "Alpheus lottini",
    pattern = "^Alpheus lottini$",
    effect = "positive",
    description = "Snapping shrimp (mutualist)"
  ),
  list(
    name = "All Alpheus spp.",
    pattern = "^Alpheus",
    effect = "positive",
    description = "All snapping shrimp species"
  ),
  list(
    name = "Cymo + Luniella (xanthids)",
    pattern = "Cymo|Luniella",
    effect = "negative",
    description = "Harmful xanthid crabs (tissue consumers)"
  )
)

# Calculate abundances for each key species/group
cat("Calculating abundances for key species...\n")

key_species_abundances <- lapply(key_exp_species, function(sp) {
  counts <- cafi_clean %>%
    filter(grepl(sp$pattern, otu, ignore.case = TRUE)) %>%
    group_by(coral_id) %>%
    summarise(abundance = n(), .groups = "drop")

  # Complete with zeros for corals without this species
  complete_counts <- coral_data %>%
    dplyr::select(coral_id, site, volume) %>%
    left_join(counts, by = "coral_id") %>%
    mutate(abundance = replace_na(abundance, 0))

  list(
    name = sp$name,
    effect = sp$effect,
    description = sp$description,
    data = complete_counts,
    total_individuals = sum(complete_counts$abundance),
    n_corals_present = sum(complete_counts$abundance > 0)
  )
})

names(key_species_abundances) <- sapply(key_exp_species, function(x) x$name)

# Report occurrence summary
cat("\nKey species occurrence summary:\n")
cat(sprintf("  %-30s %8s %12s\n", "Species/Group", "Total N", "Corals Present"))
cat("  ", strrep("-", 52), "\n", sep = "")

for (sp in key_species_abundances) {
  cat(sprintf("  %-30s %8d %12d (%.1f%%)\n",
              sp$name,
              sp$total_individuals,
              sp$n_corals_present,
              100 * sp$n_corals_present / nrow(coral_data)))
}

# Fit scaling models with RELAXED thresholds for key species
# Using min_nonzero = 5 to include rarer but important species
cat("\nFitting scaling models (relaxed threshold: min 5 corals with species present)...\n\n")

key_species_results <- lapply(key_species_abundances, function(sp) {
  result <- fit_scaling_model(
    sp$data,
    sp$name,
    min_nonzero = 5,  # Relaxed threshold for key species
    n_boot = 2000
  )

  # Add experimental paper metadata
  result$exp_effect <- sp$effect
  result$description <- sp$description

  # Report
  if (result$converged && !is.na(result$beta)) {
    cat(sprintf("  %-30s β = %6.3f [%5.3f, %5.3f], p vs 1: %s\n",
                sp$name,
                result$beta,
                result$boot_ci_lower,
                result$boot_ci_upper,
                format.pval(result$p_vs_1, 3)))
    cat(sprintf("    → %s (experimental effect: %s on coral condition)\n",
                result$interpretation, sp$effect))
  } else {
    cat(sprintf("  %-30s %s\n", sp$name, result$interpretation))
  }

  result
})

names(key_species_results) <- names(key_species_abundances)
all_results$key_exp_species <- key_species_results

# Create summary table for key species
cat("\n--- KEY SPECIES SCALING SUMMARY ---\n\n")

key_species_summary <- map_df(key_species_results, function(r) {
  tibble(
    species = r$response,
    exp_effect = r$exp_effect,
    description = r$description,
    n_total = r$n_corals,
    n_present = r$n_nonzero,
    total_abundance = r$total_abundance,
    beta = r$beta,
    ci_lower = r$boot_ci_lower,
    ci_upper = r$boot_ci_upper,
    p_vs_1 = r$p_vs_1,
    interpretation = r$interpretation
  )
})

# Add Field of Dreams / Redistribution classification
key_species_summary <- key_species_summary %>%
  mutate(
    scaling_pattern = case_when(
      is.na(beta) ~ "Insufficient data",
      interpretation == "Field of Dreams" ~ "Field of Dreams (β ≈ 1)",
      grepl("Redistribution", interpretation) ~ "Propagule Redistribution (β < 1)",
      TRUE ~ interpretation
    )
  )

print(key_species_summary %>%
        dplyr::select(species, exp_effect, n_present, beta, ci_lower, ci_upper, scaling_pattern) %>%
        mutate(across(where(is.numeric) & !c(n_present), ~round(., 3))))

# Save key species scaling table
save_table(key_species_summary, "key_species_scaling_experimental")
cat("\n  Saved: key_species_scaling_experimental.csv\n")

# Create visualization: Key species forest plot
cat("\nCreating key species scaling forest plot...\n")

key_plot_data <- key_species_summary %>%
  filter(!is.na(beta)) %>%
  mutate(
    species = factor(species, levels = rev(species)),
    exp_effect_label = ifelse(exp_effect == "positive",
                               "Beneficial (positive effect on coral)",
                               "Harmful (negative effect on coral)")
  )

if (nrow(key_plot_data) > 0) {
  p_key_forest <- ggplot(key_plot_data, aes(x = beta, y = species, color = exp_effect_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3, linewidth = 0.8, orientation = "y") +
    geom_point(size = 4) +
    scale_color_manual(
      values = c("Beneficial (positive effect on coral)" = "#0072B2",
                 "Harmful (negative effect on coral)" = "#D55E00"),
      name = "Effect on coral condition\n(from experimental paper)"
    ) +
    annotate("text", x = 1.15, y = 0.5, label = "Field of Dreams\n(β = 1)",
             hjust = 0, size = 3, color = "gray40") +
    annotate("text", x = 0.5, y = 0.5, label = "Redistribution\n(β < 1)",
             hjust = 0.5, size = 3, color = "gray40") +
    labs(
      title = "Scaling Patterns for Key Species from Experimental Paper",
      subtitle = "Comparing Field of Dreams (β = 1) vs. Propagule Redistribution (β < 1)",
      x = "Scaling exponent (β) with 95% bootstrap CI",
      y = NULL,
      caption = "Dashed line = Field of Dreams prediction (β = 1)\nSpecies colored by documented effect on coral condition"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      plot.caption = element_text(hjust = 0, size = 9),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 11),
      plot.title = element_text(face = "bold", size = 13)
    )

  ggsave(file.path(FIG_DIR, "key_species_scaling_forest.png"), p_key_forest,
         width = 10, height = 6, dpi = 300, bg = "white")
  cat("  Saved: key_species_scaling_forest.png\n")
}

# Summary interpretation
cat("\n--- INTERPRETATION: HETEROGENEITY IN KEY SPECIES SCALING ---\n\n")

n_analyzed <- sum(!is.na(key_species_summary$beta))
n_fod <- sum(key_species_summary$interpretation == "Field of Dreams", na.rm = TRUE)
n_redist <- sum(grepl("Redistribution", key_species_summary$interpretation), na.rm = TRUE)

cat("Of", n_analyzed, "key species/groups analyzed:\n")
cat("  Field of Dreams (β ≈ 1):", n_fod, "\n")
cat("  Propagule Redistribution (β < 1):", n_redist, "\n\n")

cat("Comparison with experimental paper predictions:\n")
cat("  The experimental work (Stier et al.) predicts heterogeneous scaling responses\n")
cat("  across species, with some following Field of Dreams and others showing\n")
cat("  Propagule Redistribution. This survey tests whether the same pattern holds\n")
cat("  for observational data across natural coral size variation.\n\n")

# Check if already in individual species results
already_in_main <- key_species_summary %>%
  filter(species %in% names(species_results))

if (nrow(already_in_main) > 0) {
  cat("Note: Some key species were also included in main species analysis (met occurrence thresholds).\n")
}

# ############################################################################
#                    PART 7: SUMMARY & SYNTHESIS
# ############################################################################

cat("\n############################################################\n")
cat("PART 7: SUMMARY & SYNTHESIS\n")
cat("############################################################\n\n")

# Compile all results into a single table
compile_results <- function(results_list) {
  rows <- list()

  # Community-level metrics
  rows[[1]] <- result_to_row(results_list$total_abundance) %>% mutate(category = "Community")
  rows[[2]] <- result_to_row(results_list$species_richness) %>% mutate(category = "Community")
  rows[[3]] <- result_to_row(results_list$shannon_diversity) %>% mutate(category = "Community")

  # Taxonomic groups
  for (name in names(results_list$taxonomic_groups)) {
    rows[[length(rows) + 1]] <- result_to_row(results_list$taxonomic_groups[[name]]) %>%
      mutate(category = "Taxonomic Group")
  }

  # Functional groups
  for (name in names(results_list$functional_groups)) {
    rows[[length(rows) + 1]] <- result_to_row(results_list$functional_groups[[name]]) %>%
      mutate(category = "Functional Group")
  }

  # Individual species
  for (name in names(results_list$individual_species)) {
    rows[[length(rows) + 1]] <- result_to_row(results_list$individual_species[[name]]) %>%
      mutate(category = "Species")
  }

  bind_rows(rows)
}

all_results_df <- compile_results(all_results)

# Apply FDR correction within species-level tests (multiple testing)
all_results_df <- all_results_df %>%
  group_by(category) %>%
  mutate(
    p_vs_1_fdr = p.adjust(p_vs_1, method = "BH"),
    significant_fdr = !is.na(p_vs_1_fdr) & p_vs_1_fdr < 0.05
  ) %>%
  ungroup()

cat("FDR CORRECTION (Benjamini-Hochberg within category):\n")
fdr_summary <- all_results_df %>%
  filter(!is.na(p_vs_1)) %>%
  group_by(category) %>%
  summarise(
    n_tests = n(),
    n_sig_raw = sum(p_vs_1 < 0.05, na.rm = TRUE),
    n_sig_fdr = sum(p_vs_1_fdr < 0.05, na.rm = TRUE),
    .groups = "drop"
  )
print(fdr_summary)
cat("\n")

# Summary statistics
valid_results <- all_results_df %>% filter(!is.na(beta))

cat("OVERALL SUMMARY:\n")
cat("  Total analyses:", nrow(all_results_df), "\n")
cat("  Successful fits:", nrow(valid_results), "\n\n")

# Summary by interpretation
interpretation_summary <- valid_results %>%
  group_by(interpretation) %>%
  summarise(
    n = n(),
    mean_beta = mean(beta),
    .groups = "drop"
  )

cat("INTERPRETATION DISTRIBUTION:\n")
print(interpretation_summary)
cat("\n")

# Summary by category
category_summary <- valid_results %>%
  group_by(category) %>%
  summarise(
    n = n(),
    mean_beta = mean(beta),
    median_beta = median(beta),
    .groups = "drop"
  )

cat("\nBY CATEGORY:\n")
print(category_summary)
cat("\n")

# Meta-analytic approach: weighted mean of species betas
# Note: Simple t-test is inappropriate because species betas have different
# precision (SEs) and share the same coral volume predictor (non-independence)
species_data <- valid_results %>%
  filter(category == "Species", !is.na(se), se > 0)

if (nrow(species_data) >= 3) {
  # Inverse-variance weighting (more precise estimates get more weight)
  weights <- 1 / (species_data$se^2)
  weighted_mean_beta <- sum(species_data$beta * weights) / sum(weights)
  weighted_se <- sqrt(1 / sum(weights))

  # 95% CI for weighted mean
  weighted_ci_lower <- weighted_mean_beta - 1.96 * weighted_se
  weighted_ci_upper <- weighted_mean_beta + 1.96 * weighted_se

  # Z-test: does weighted mean differ from 1?
  z_stat <- (weighted_mean_beta - 1) / weighted_se
  p_val_weighted <- 2 * pnorm(-abs(z_stat))

  cat("\nSPECIES-LEVEL WEIGHTED META-ANALYSIS (H0: mean β = 1):\n")
  cat("  Weighted mean β =", round(weighted_mean_beta, 3), "\n")
  cat("  95% CI: [", round(weighted_ci_lower, 3), ",", round(weighted_ci_upper, 3), "]\n")
  cat("  z =", round(z_stat, 2), ", p =", format.pval(p_val_weighted, 3), "\n")
  cat("  n species:", nrow(species_data), "\n")

  # Heterogeneity: I² statistic
  Q_stat <- sum(weights * (species_data$beta - weighted_mean_beta)^2)
  df_Q <- nrow(species_data) - 1
  I2 <- max(0, (Q_stat - df_Q) / Q_stat) * 100
  cat("  Heterogeneity: I² =", round(I2, 1), "% (Q =", round(Q_stat, 1), ", df =", df_Q, ")\n")

  if (p_val_weighted < 0.05 && weighted_mean_beta < 1) {
    cat("  Conclusion: Species on average show REDISTRIBUTION (β < 1)\n")
  } else if (p_val_weighted < 0.05 && weighted_mean_beta > 1) {
    cat("  Conclusion: Species on average show SUPER-LINEAR scaling (β > 1)\n")
  } else {
    cat("  Conclusion: Species on average support FIELD OF DREAMS (β ≈ 1)\n")
  }

  # Also report unweighted t-test for comparison (with caveat)
  species_betas <- species_data$beta
  t_test <- t.test(species_betas, mu = 1)
  cat("\n  [Unweighted t-test for comparison - interpret cautiously due to precision heterogeneity]\n")
  cat("  Unweighted mean β =", round(mean(species_betas), 3), "±", round(sd(species_betas)/sqrt(length(species_betas)), 3), "\n")
  cat("  t =", round(t_test$statistic, 2), ", df =", round(t_test$parameter, 1),
      ", p =", format.pval(t_test$p.value, 3), "\n")

  # Store for summary
  meta_results <- list(
    weighted_mean = weighted_mean_beta,
    weighted_se = weighted_se,
    weighted_ci = c(weighted_ci_lower, weighted_ci_upper),
    z_stat = z_stat,
    p_val = p_val_weighted,
    I2 = I2,
    n_species = nrow(species_data)
  )
} else {
  t_test <- NULL
  meta_results <- NULL
  species_betas <- NULL
}

# ############################################################################
#                    FIGURES
# ############################################################################

cat("\n############################################################\n")
cat("CREATING FIGURES\n")
cat("############################################################\n\n")

# Define interpretation colors (colorblind-safe)
interpretation_colors <- c(
  "Redistribution (β < 1)" = "#D55E00",
  "Field of Dreams (β ≈ 1)" = "#009E73",
  "Super-linear (β > 1)" = "#0072B2"
)

# ============================================================================
# FIGURE 1: Community-Level Scaling (Total, Richness, Shannon)
# ============================================================================

cat("Creating Figure 1: Community-level scaling...\n")

# Total abundance plot
p_total <- ggplot(coral_data, aes(x = volume, y = total_cafi, color = site)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(aes(group = 1), method = "glm.nb", formula = y ~ log10(x),
              se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  scale_color_manual(values = SITE_COLORS) +
  labs(
    x = expression("Coral Volume (cm"^3*")"),
    y = "Total CAFI",
    title = "A. Total Abundance",
    subtitle = sprintf("β = %.2f [%.2f, %.2f]",
                       total_result$beta, total_result$ci_lower, total_result$ci_upper)
  ) +
  theme(legend.position = c(0.15, 0.85))

# Richness plot
p_richness <- ggplot(coral_data, aes(x = volume, y = otu_richness, color = site)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(aes(group = 1), method = "glm.nb", formula = y ~ log10(x),
              se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  scale_color_manual(values = SITE_COLORS) +
  labs(
    x = expression("Coral Volume (cm"^3*")"),
    y = "Species Richness",
    title = "B. Species Richness",
    subtitle = sprintf("β = %.2f [%.2f, %.2f]",
                       richness_result$beta, richness_result$ci_lower, richness_result$ci_upper)
  ) +
  theme(legend.position = "none")

# Shannon plot
p_shannon <- coral_data %>%
  filter(shannon > 0) %>%
  ggplot(aes(x = volume, y = shannon, color = site)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ log10(x),
              se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = SITE_COLORS) +
  labs(
    x = expression("Coral Volume (cm"^3*")"),
    y = "Shannon H'",
    title = "C. Shannon Diversity",
    subtitle = if(!is.na(shannon_result$beta)) sprintf("β = %.2f [%.2f, %.2f]",
                       shannon_result$beta, shannon_result$ci_lower, shannon_result$ci_upper) else "Insufficient data"
  ) +
  theme(legend.position = "none")

# Combine community plots
fig_community <- (p_total + p_richness + p_shannon) +
  plot_annotation(
    title = "Community-Level Scaling: CAFI vs Coral Volume",
    subtitle = "Negative binomial GLM | Dashed line = β = 1 (Field of Dreams)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(FIG_DIR, "community_scaling.png"), fig_community,
       width = 14, height = 5, dpi = 300, bg = "white")
cat("  Saved: community_scaling.png\n")

# ============================================================================
# FIGURE 2: Taxonomic & Functional Group Forest Plot
# ============================================================================

cat("Creating Figure 2: Taxonomic & functional group scaling...\n")

# Combine taxonomic and functional results
group_results <- bind_rows(
  map_df(names(all_results$taxonomic_groups), ~result_to_row(all_results$taxonomic_groups[[.x]])) %>%
    mutate(group_type = "Taxonomic"),
  map_df(names(all_results$functional_groups), ~result_to_row(all_results$functional_groups[[.x]])) %>%
    mutate(group_type = "Functional")
) %>%
  filter(!is.na(beta))

if (nrow(group_results) > 0) {
  # Create unique labels for plotting (append group type to avoid duplicates)
  group_results <- group_results %>%
    mutate(plot_label = paste0(response, " (", substr(group_type, 1, 1), ")"))

  p_groups <- group_results %>%
    mutate(plot_label = factor(plot_label, levels = unique(plot_label[order(beta)]))) %>%
    ggplot(aes(x = beta, y = plot_label, color = interpretation, shape = group_type)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.3, linewidth = 0.6) +
    geom_point(size = 4) +
    scale_color_manual(values = interpretation_colors, name = "Interpretation") +
    scale_shape_manual(values = c("Taxonomic" = 16, "Functional" = 17), name = "Group Type") +
    labs(
      x = expression("Scaling Exponent (" * beta * ")"),
      y = NULL,
      title = "Scaling by Taxonomic and Functional Groups",
      subtitle = "95% CI shown | Vertical line = Field of Dreams (β = 1)"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 10)
    )

  ggsave(file.path(FIG_DIR, "group_scaling_forest.png"), p_groups,
         width = 10, height = 7, dpi = 300, bg = "white")
  cat("  Saved: group_scaling_forest.png\n")
}

# ============================================================================
# FIGURE 3: Species-Level Forest Plot
# ============================================================================

cat("Creating Figure 3: Species-level scaling forest plot...\n")

species_plot_data <- valid_results %>%
  filter(category == "Species") %>%
  left_join(species_functional %>% dplyr::select(otu, functional_group, type),
            by = c("response" = "otu")) %>%
  mutate(response = factor(response, levels = response[order(beta)]))

if (nrow(species_plot_data) > 0) {
  p_species_forest <- ggplot(species_plot_data, aes(x = beta, y = response, color = interpretation)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray80") +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.3, linewidth = 0.6) +
    geom_point(size = 3) +
    scale_color_manual(values = interpretation_colors, name = "Scaling Pattern") +
    labs(
      x = expression("Scaling Exponent (" * beta * ")"),
      y = NULL,
      title = "Species-Level Scaling: Abundance vs Coral Volume",
      subtitle = sprintf("N = %d species | Mean β = %.2f | Test vs 1: p = %s",
                         nrow(species_plot_data),
                         mean(species_plot_data$beta),
                         if(!is.null(t_test)) format.pval(t_test$p.value, 2) else "NA")
    ) +
    theme(
      legend.position = c(0.85, 0.15),
      legend.background = element_rect(fill = "white", color = "gray80"),
      axis.text.y = element_text(face = "italic", size = 9)
    ) +
    coord_cartesian(xlim = c(-0.5, 3))

  ggsave(file.path(FIG_DIR, "species_scaling_forest.png"), p_species_forest,
         width = 10, height = 8, dpi = 300, bg = "white")
  cat("  Saved: species_scaling_forest.png\n")
}

# ============================================================================
# FIGURE 4: Species Scatter Panels
# ============================================================================

cat("Creating Figure 4: Species scatter panels...\n")

# Select top species by sample size
top_species <- valid_results %>%
  filter(category == "Species") %>%
  arrange(desc(n_nonzero)) %>%
  head(12)

make_species_panel <- function(sp_name, sp_result) {
  sp_data <- species_scaling_data %>%
    filter(otu == sp_name, abundance > 0)

  beta_val <- sp_result$beta
  interp <- sp_result$interpretation

  title_color <- case_when(
    interp == "Redistribution (β < 1)" ~ "#D55E00",
    interp == "Super-linear (β > 1)" ~ "#0072B2",
    TRUE ~ "#009E73"
  )

  ggplot(sp_data, aes(x = volume, y = abundance, color = site)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(aes(group = 1), method = "glm.nb", formula = y ~ log10(x),
                se = TRUE, color = "black", linewidth = 0.8) +
    scale_x_log10(labels = scales::comma) +
    scale_y_continuous() +
    scale_color_manual(values = SITE_COLORS) +
    labs(
      x = expression("Volume (cm"^3*")"),
      y = "Abundance",
      title = str_replace(sp_name, "_", " ")
    ) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("β = %.2f", beta_val),
             hjust = 1.1, vjust = 1.3, size = 3, fontface = "italic") +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold.italic", size = 9, color = title_color)
    )
}

if (nrow(top_species) > 0) {
  species_panels <- map2(top_species$response,
                          map(top_species$response, ~all_results$individual_species[[.x]]),
                          make_species_panel)

  fig_species <- wrap_plots(species_panels, ncol = 4) +
    plot_annotation(
      title = "Species-Specific Scaling Relationships",
      subtitle = "Top 12 species by sample size | Title color indicates interpretation",
      caption = paste0("Red = Redistribution (β < 1), Green = Field of Dreams (β ≈ 1), Blue = Super-linear (β > 1)"),
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11)
      )
    )

  ggsave(file.path(FIG_DIR, "species_scaling_panels.png"), fig_species,
         width = 14, height = 10, dpi = 300, bg = "white")
  cat("  Saved: species_scaling_panels.png\n")
}

# ============================================================================
# FIGURE 5: Summary Visualization
# ============================================================================

cat("Creating Figure 5: Summary visualization...\n")

# Beta distribution by category
p_beta_dist <- valid_results %>%
  filter(category %in% c("Species", "Taxonomic Group", "Functional Group")) %>%
  ggplot(aes(x = beta, fill = category)) +
  geom_histogram(bins = 15, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  scale_fill_viridis_d(option = "D", name = "Category") +
  labs(
    x = expression("Scaling Exponent (" * beta * ")"),
    y = "Count",
    title = "Distribution of Scaling Exponents",
    subtitle = "Dashed line = Field of Dreams (β = 1)"
  ) +
  theme(legend.position = "right")

# Beta by interpretation
p_interp <- valid_results %>%
  ggplot(aes(x = interpretation, y = beta, fill = interpretation)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = interpretation_colors) +
  labs(
    x = NULL,
    y = expression("Scaling Exponent (" * beta * ")"),
    title = "Scaling by Interpretation"
  ) +
  theme(legend.position = "none") +
  coord_flip()

fig_summary <- (p_beta_dist / p_interp) +
  plot_annotation(
    title = "Scaling Analysis Summary",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(FIG_DIR, "scaling_summary.png"), fig_summary,
       width = 10, height = 10, dpi = 300, bg = "white")
cat("  Saved: scaling_summary.png\n")

# ############################################################################
#                    SAVE TABLES
# ############################################################################

cat("\n############################################################\n")
cat("SAVING TABLES\n")
cat("############################################################\n\n")

# Full results table
full_table <- all_results_df %>%
  mutate(across(where(is.numeric) & !contains("n_"), ~round(., 3))) %>%
  arrange(category, desc(n_nonzero))

save_table(full_table, "scaling_results_all")
cat("  Saved: scaling_results_all.csv\n")

# Summary by category
save_table(category_summary, "scaling_summary_by_category")
cat("  Saved: scaling_summary_by_category.csv\n")

# Interpretation summary
save_table(interpretation_summary, "scaling_interpretation_summary")
cat("  Saved: scaling_interpretation_summary.csv\n")

# ############################################################################
#                    SAVE RESULTS OBJECT
# ############################################################################

scaling_analysis_results <- list(
  # Sample info
  n_corals = nrow(coral_data),
  total_cafi = sum(coral_data$total_cafi),

  # Full results table
  all_results = all_results_df,

  # Model objects
  models = all_results,

  # Summary statistics
  category_summary = category_summary,
  interpretation_summary = interpretation_summary,

  # T-test for species
  species_t_test = t_test,

  # Species metadata
  species_summary = species_summary,

  # Key species from experimental paper
  key_species_scaling = key_species_summary
)

save_object(scaling_analysis_results, "scaling_analysis_results")
cat("  Saved: scaling_analysis_results.rds\n")

# ############################################################################
#                    FINAL SUMMARY
# ############################################################################

cat("\n")
cat("============================================================\n")
cat("    SCALING ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("SAMPLE:\n")
cat("  Corals:", nrow(coral_data), "\n")
cat("  Total CAFI:", sum(coral_data$total_cafi), "\n")
cat("  Species analyzed:", sum(valid_results$category == "Species"), "\n\n")

cat("KEY RESULTS:\n\n")

cat("1. TOTAL ABUNDANCE:\n")
cat("   β =", round(total_result$beta, 3), "[", round(total_result$ci_lower, 3), ",",
    round(total_result$ci_upper, 3), "]\n")
cat("   Interpretation:", total_result$interpretation, "\n\n")

cat("2. SPECIES RICHNESS:\n")
cat("   β =", round(richness_result$beta, 3), "[", round(richness_result$ci_lower, 3), ",",
    round(richness_result$ci_upper, 3), "]\n")
cat("   Interpretation:", richness_result$interpretation, "\n\n")

cat("3. SPECIES-LEVEL PATTERNS:\n")
if (!is.null(meta_results)) {
  cat("   Weighted mean β =", round(meta_results$weighted_mean, 3), "\n")
  cat("   95% CI: [", round(meta_results$weighted_ci[1], 3), ",", round(meta_results$weighted_ci[2], 3), "]\n")
  cat("   Meta-analysis (H0: β = 1): p =", format.pval(meta_results$p_val, 3), "\n")
  cat("   Heterogeneity I² =", round(meta_results$I2, 1), "%\n")
}
cat("   Distribution:\n")
for (i in 1:nrow(interpretation_summary)) {
  cat("     ", interpretation_summary$interpretation[i], ":", interpretation_summary$n[i], "\n")
}

cat("\n4. KEY SPECIES FROM EXPERIMENTAL PAPER:\n")
if (exists("key_species_summary") && nrow(key_species_summary) > 0) {
  for (i in 1:nrow(key_species_summary)) {
    row <- key_species_summary[i, ]
    if (!is.na(row$beta)) {
      cat(sprintf("   %-25s β = %.3f [%.3f, %.3f] → %s\n",
                  row$species, row$beta, row$ci_lower, row$ci_upper, row$scaling_pattern))
    } else {
      cat(sprintf("   %-25s %s\n", row$species, row$interpretation))
    }
  }
  cat("   (See output/tables/key_species_scaling_experimental.csv for full details)\n")
}

cat("\nOUTPUT FILES:\n")
cat("  Figures: output/figures/05_scaling/\n")
cat("    - community_scaling.png\n")
cat("    - group_scaling_forest.png\n")
cat("    - species_scaling_forest.png\n")
cat("    - species_scaling_panels.png\n")
cat("    - scaling_summary.png\n")
cat("    - key_species_scaling_forest.png (NEW)\n")
cat("  Tables: output/tables/\n")
cat("    - scaling_results_all.csv\n")
cat("    - scaling_summary_by_category.csv\n")
cat("    - scaling_interpretation_summary.csv\n")
cat("    - key_species_scaling_experimental.csv (NEW)\n\n")
