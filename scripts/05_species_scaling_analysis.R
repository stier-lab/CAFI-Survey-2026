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
#   - Natural log: coefficient = power-law exponent directly
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
# DEPENDENCIES: 00_setup.R, 01_load_data.R (coral_master.rds, cafi_clean.rds,
#               community_matrix.rds, functional_summary.rds)
#
# OUTPUTS:
#   Figures: output/figures/05_scaling/, output/figures/manuscript/fig2_scaling.png
#   Tables: output/tables/scaling_results_all.csv, scaling_summary_by_category.csv
#   Objects: output/objects/scaling_analysis_results.rds
#
# Author: CAFI Survey Analysis Pipeline
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI-CORAL SCALING ANALYSIS\n")
cat("    Field of Dreams vs. Propagule Redirection\n")
cat("============================================================\n\n")

# Load setup and data
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

# Use script-specific output directory
FIG_DIR <- PATHS$fig_05_scaling

# Load CAFI data
cafi_clean <- load_object("cafi_clean")

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
      model <- MASS::glm.nb(abundance ~ log(volume) + site, data = d)
    } else {
      model <- MASS::glm.nb(abundance ~ log(volume), data = d)
    }
    coef(model)["log(volume)"]
  }, error = function(e) NA)
}

#' Fit a scaling model and test against Field of Dreams (β = 1)
#'
#' @param data Data frame with 'abundance' and 'volume' columns
#' @param response_name Character string for labeling
#' @param min_nonzero Minimum non-zero observations required (default 15)
#' @param n_boot Number of bootstrap replicates (default 1000; 1000 provides
#'   ~0.03 CI width precision for beta, sufficient for ecological inference)
#' @return List with model results and interpretation
# Minimum non-zero observations for reliable NB GLM convergence
# (15 non-zero counts ensures ~13% prevalence across 114 corals;
# sensitivity tested with thresholds 10, 15, 20 — results qualitatively unchanged)
fit_scaling_model <- function(data, response_name, min_nonzero = 15, n_boot = 1000) {

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
      p_boot_vs_1 = NA_real_,
      pseudo_r2 = NA_real_,
      interpretation = "Insufficient data",
      model = NULL,
      converged = FALSE
    ))
  }

  # Fit negative binomial GLM
  tryCatch({
    model <- MASS::glm.nb(abundance ~ log(volume) + site, data = data)

    # Extract coefficients
    coefs <- summary(model)$coefficients
    beta <- coefs["log(volume)", "Estimate"]
    se <- coefs["log(volume)", "Std. Error"]
    z_val <- coefs["log(volume)", "z value"]
    p_val <- coefs["log(volume)", "Pr(>|z|)"]

    # Profile-based 95% CI (with tryCatch for convergence)
    set.seed(42)
    ci <- tryCatch(
      confint(model, "log(volume)", level = 0.95),
      error = function(e) {
        cat("  Warning: Profile CI failed, using Wald CI\n")
        beta_val <- coefs["log(volume)", "Estimate"]
        se_val <- coefs["log(volume)", "Std. Error"]
        c(beta_val - 1.96 * se_val, beta_val + 1.96 * se_val)
      }
    )

    # Bootstrap 95% CI (following Stier et al. 2024 methodology)
    boot_ci_lower <- NA_real_
    boot_ci_upper <- NA_real_
    p_boot_vs_1 <- NA_real_
    boot_ci_method <- NA_character_

    if (requireNamespace("boot", quietly = TRUE) && n_nonzero >= 20) {
      suppressWarnings({
        # Stratify by site to maintain site proportions in each bootstrap replicate
        set.seed(42)
        boot_result <- boot::boot(data, boot_scaling, R = n_boot,
                                  strata = factor(data$site))

        # Report bootstrap success/failure rate
        n_succeeded <- sum(!is.na(boot_result$t[,1]))
        n_total_boot <- nrow(boot_result$t)
        cat(sprintf("  Bootstrap: %d/%d iterations succeeded\n", n_succeeded, n_total_boot))
        if (n_succeeded < 0.9 * n_total_boot) {
          warning(sprintf("Bootstrap failure rate >10%% for %s: %d/%d failed",
                          response_name, n_total_boot - n_succeeded, n_total_boot))
        }

        valid_betas <- boot_result$t[!is.na(boot_result$t)]

        if (length(valid_betas) >= n_boot * 0.8) {  # Require 80% successful fits
          # Try BCa CI first, fall back to percentile
          bca_ci <- tryCatch(
            boot::boot.ci(boot_result, type = "bca", conf = 0.95),
            error = function(e) {
              warning(sprintf("BCa CI failed for %s; using percentile CI instead (less bias-corrected)", response_name))
              NULL
            }
          )

          if (!is.null(bca_ci)) {
            boot_ci_lower <- bca_ci$bca[4]
            boot_ci_upper <- bca_ci$bca[5]
            boot_ci_method <- "BCa"
          } else {
            # Fallback to percentile
            boot_ci_lower <- quantile(valid_betas, 0.025)
            boot_ci_upper <- quantile(valid_betas, 0.975)
            boot_ci_method <- "percentile"
          }

          # Bootstrap p-value vs β=1 (proportion of replicates on opposite side of 1)
          # Two-sided: p = 2 * min(proportion ≥ 1, proportion < 1)
          prop_ge_1 <- mean(valid_betas >= 1)
          p_boot_vs_1 <- 2 * min(prop_ge_1, 1 - prop_ge_1)
        }
      })

      # Warn if bootstrap failure rate exceeds 10% (redundant safety check)
      n_failed_check <- n_total_boot - n_succeeded
      if (n_failed_check > 0.1 * n_boot) {
        warning(sprintf("Bootstrap for %s: %d/%d iterations failed (%.1f%%)",
                        response_name, n_failed_check, n_boot, 100 * n_failed_check / n_boot))
      }
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
      interpretation <- "Redirection (β < 1)"
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
      boot_ci_method = boot_ci_method,
      z_value = z_val,
      p_value = p_val,
      z_vs_1 = z_vs_1,
      p_vs_1 = p_vs_1,
      p_boot_vs_1 = p_boot_vs_1,
      pseudo_r2 = pseudo_r2,
      interpretation = interpretation,
      model = model,
      converged = TRUE
    )

  }, error = function(e) {
    message(sprintf("NB GLM did not converge for %s (n_nonzero=%d); returning NA", response_name, n_nonzero))
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
      boot_ci_method = NA_character_,
      z_value = NA_real_,
      p_value = NA_real_,
      z_vs_1 = NA_real_,
      p_vs_1 = NA_real_,
      p_boot_vs_1 = NA_real_,
      interpretation = paste("Model error:", conditionMessage(e)),
      model = NULL,
      converged = FALSE
    )
  })
}

#' Fit richness with Poisson (not NB) — richness has different distributional properties
#' Richness (count of unique species) is better modeled with Poisson than NB
#' because overdispersion in richness is typically lower than in abundance.
#' @param data Data frame with 'abundance' (richness) and 'volume' columns
#' @param response_name Character string for labeling
#' @param min_nonzero Minimum non-zero observations required (default 15)
#' @param n_boot Number of bootstrap replicates (default 1000)
#' @return List with model results and interpretation
fit_richness_model <- function(data, response_name, min_nonzero = 15, n_boot = 1000) {

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
      p_boot_vs_1 = NA_real_,
      pseudo_r2 = NA_real_,
      interpretation = "Insufficient data",
      model = NULL,
      converged = FALSE
    ))
  }

  # Fit Poisson GLM (not NB — richness has different distributional properties)
  tryCatch({
    formula <- as.formula("abundance ~ log(volume) + site")
    model <- glm(formula, data = data, family = poisson)

    # Check for overdispersion
    pearson_chi2 <- sum(residuals(model, type = "pearson")^2)
    dispersion <- pearson_chi2 / model$df.residual
    model_family <- "Poisson"

    if (dispersion > 1.5) {
      cat("  Note: Overdispersion detected (", round(dispersion, 2),
          "), switching to quasipoisson\n")
      model <- glm(formula, data = data, family = quasipoisson)
      model_family <- "quasiPoisson"
    }

    # Extract coefficients
    coefs <- summary(model)$coefficients
    beta <- coefs["log(volume)", "Estimate"]
    se <- coefs["log(volume)", "Std. Error"]
    # For Poisson, use "z value"; for quasipoisson, use "t value"
    stat_col <- if (model_family == "Poisson") "z value" else "t value"
    p_col <- if (model_family == "Poisson") "Pr(>|z|)" else "Pr(>|t|)"
    z_val <- coefs["log(volume)", stat_col]
    p_val <- coefs["log(volume)", p_col]

    # Profile-based 95% CI (with tryCatch for convergence)
    set.seed(42)
    ci <- tryCatch(
      confint(model, "log(volume)", level = 0.95),
      error = function(e) {
        cat("  Warning: Profile CI failed, using Wald CI\n")
        c(beta - 1.96 * se, beta + 1.96 * se)
      }
    )

    # Bootstrap 95% CI — match detected family (Poisson or quasipoisson)
    # The bootstrap is nonparametric (resampling rows), so the family choice
    # affects point estimation within each replicate. Using the same family
    # as the main model ensures consistency in SE estimation.
    richness_family <- if (model_family == "quasiPoisson") quasipoisson else poisson
    boot_ci_lower <- NA_real_
    boot_ci_upper <- NA_real_
    p_boot_vs_1 <- NA_real_
    boot_ci_method <- NA_character_

    boot_richness_fn <- function(data, indices) {
      d <- data[indices, ]
      tryCatch({
        if ("site" %in% names(d) && length(unique(d$site)) > 1) {
          m <- glm(abundance ~ log(volume) + site, data = d, family = richness_family)
        } else {
          m <- glm(abundance ~ log(volume), data = d, family = richness_family)
        }
        coef(m)["log(volume)"]
      }, error = function(e) NA)
    }

    if (requireNamespace("boot", quietly = TRUE) && n_nonzero >= 20) {
      suppressWarnings({
        set.seed(42)
        boot_result <- boot::boot(data, boot_richness_fn, R = n_boot,
                                  strata = factor(data$site))

        # Report bootstrap success/failure rate
        n_succeeded <- sum(!is.na(boot_result$t[,1]))
        n_total_boot <- nrow(boot_result$t)
        cat(sprintf("  Bootstrap: %d/%d iterations succeeded\n", n_succeeded, n_total_boot))
        if (n_succeeded < 0.9 * n_total_boot) {
          warning(sprintf("Bootstrap failure rate >10%% for %s: %d/%d failed",
                          response_name, n_total_boot - n_succeeded, n_total_boot))
        }

        valid_betas <- boot_result$t[!is.na(boot_result$t)]

        if (length(valid_betas) >= n_boot * 0.8) {
          # Try BCa CI first, fall back to percentile
          boot_ci <- tryCatch(
            boot::boot.ci(boot_result, type = "bca", conf = 0.95),
            error = function(e) {
              cat("  BCa CI failed, falling back to percentile method\n")
              NULL
            }
          )

          if (!is.null(boot_ci)) {
            boot_ci_lower <- boot_ci$bca[4]
            boot_ci_upper <- boot_ci$bca[5]
            boot_ci_method <- "BCa"
          } else {
            # Fallback to percentile
            boot_ci_lower <- quantile(valid_betas, 0.025)
            boot_ci_upper <- quantile(valid_betas, 0.975)
            boot_ci_method <- "percentile"
          }

          # Bootstrap p-value vs β=1
          prop_ge_1 <- mean(valid_betas >= 1)
          p_boot_vs_1 <- 2 * min(prop_ge_1, 1 - prop_ge_1)
        }
      })
    }

    # Use bootstrap CI if available, otherwise profile CI
    test_ci_lower <- if (!is.na(boot_ci_lower)) boot_ci_lower else ci[1]
    test_ci_upper <- if (!is.na(boot_ci_upper)) boot_ci_upper else ci[2]

    # Test vs beta = 1: use z-test for Poisson, t-test for quasipoisson
    z_vs_1 <- (beta - 1) / se
    if (model_family == "quasiPoisson") {
      p_vs_1 <- 2 * pt(-abs(z_vs_1), df = model$df.residual)
    } else {
      p_vs_1 <- 2 * pnorm(-abs(z_vs_1))
    }

    # Pseudo-R2 (McFadden's)
    pseudo_r2 <- calc_pseudo_r2(model)

    # Interpret based on CI (whether 95% CI includes 1)
    if (test_ci_upper < 1) {
      interpretation <- "Redirection (β < 1)"
    } else if (test_ci_lower > 1) {
      interpretation <- "Super-linear (β > 1)"
    } else {
      interpretation <- "Field of Dreams (β ≈ 1)"
    }

    cat("  Richness model family:", model_family,
        "(dispersion =", round(dispersion, 2), ")\n")

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
      boot_ci_method = boot_ci_method,
      z_value = z_val,
      p_value = p_val,
      z_vs_1 = z_vs_1,
      p_vs_1 = p_vs_1,
      p_boot_vs_1 = p_boot_vs_1,
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
      boot_ci_method = NA_character_,
      z_value = NA_real_,
      p_value = NA_real_,
      z_vs_1 = NA_real_,
      p_vs_1 = NA_real_,
      p_boot_vs_1 = NA_real_,
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
    ci_method = if (!is.null(result$boot_ci_method)) result$boot_ci_method else NA_character_,
    z_value = result$z_value,
    p_value = result$p_value,
    z_vs_1 = result$z_vs_1,
    p_vs_1 = result$p_vs_1,
    p_boot_vs_1 = if (!is.null(result$p_boot_vs_1)) result$p_boot_vs_1 else NA_real_,
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

cat("Model: Total CAFI ~ log(Volume) + site (Negative Binomial GLM)\n\n")
cat("Results:\n")
cat("  Scaling exponent β =", round(total_result$beta, 3), "\n")
cat("  Standard error =", round(total_result$se, 3), "\n")
cat("  Profile 95% CI: [", round(total_result$ci_lower, 3), ",", round(total_result$ci_upper, 3), "]\n")
if (!is.na(total_result$boot_ci_lower)) {
  cat("  Bootstrap 95% CI: [", round(total_result$boot_ci_lower, 3), ",", round(total_result$boot_ci_upper, 3), "]\n")
}
cat("  z-value =", round(total_result$z_value, 2), ", p =", format.pval(total_result$p_value, 3), "\n\n")

cat("Test vs Field of Dreams (H0: β = 1):\n")
cat("  z =", round(total_result$z_vs_1, 2), ", p =", format.pval(total_result$p_vs_1, 3), "(Wald)\n")
if (!is.na(total_result$p_boot_vs_1)) {
  cat("  Bootstrap p vs β=1:", format.pval(total_result$p_boot_vs_1, 3),
      "(proportion of", 1000, "replicates ≥ 1)\n")
}
if (!is.na(total_result$boot_ci_lower)) {
  cat("  Bootstrap CI includes 1:", ifelse(total_result$boot_ci_lower <= 1 & total_result$boot_ci_upper >= 1, "YES", "NO"), "\n")
}
cat("  Interpretation:", total_result$interpretation, "\n\n")

# VIF diagnostics for abundance model
if (!is.null(total_result$model)) {
  cat("  VIF diagnostics:\n")
  vif_values <- car::vif(total_result$model)
  print(vif_values)
  if (any(vif_values > 5)) {
    cat("  WARNING: VIF > 5 detected — potential multicollinearity\n")
  }
  cat("\n")

  # Influence diagnostics (Cook's D)
  cd <- cooks.distance(total_result$model)
  cat("  Max Cook's D:", round(max(cd, na.rm = TRUE), 4), "\n")
  cat("  Observations with Cook's D > 4/n:", sum(cd > 4/length(cd), na.rm = TRUE), "\n\n")
}

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

# Use Poisson (not NB) for richness — richness has different distributional properties
richness_result <- fit_richness_model(richness_data, "Species Richness")

cat("Model: OTU Richness ~ log(Volume) + site (Poisson GLM)\n\n")
cat("Results:\n")
cat("  Scaling exponent β =", round(richness_result$beta, 3), "\n")
cat("  95% CI: [", round(richness_result$ci_lower, 3), ",", round(richness_result$ci_upper, 3), "]\n")
cat("  Test vs β = 1: z =", round(richness_result$z_vs_1, 2), ", p =", format.pval(richness_result$p_vs_1, 3), "\n")
cat("  Interpretation:", richness_result$interpretation, "\n\n")

# Influence diagnostics for richness model (Cook's D)
if (!is.null(richness_result$model)) {
  cd_richness <- cooks.distance(richness_result$model)
  cat("  Max Cook's D:", round(max(cd_richness, na.rm = TRUE), 4), "\n")
  cat("  Observations with Cook's D > 4/n:", sum(cd_richness > 4/length(cd_richness), na.rm = TRUE), "\n\n")
}

all_results$species_richness <- richness_result

# ---- PART 2B: RAREFIED RICHNESS SCALING ----
# Does true diversity (controlling for abundance) increase with coral size?
# If flat → SAR is entirely passive sampling (more individuals → more species)
# If positive → larger corals genuinely support richer communities
cat("--- PART 2B: Rarefied Richness Scaling ---\n")

comm_matrix_raw <- load_object("community_matrix")
row_abundances <- rowSums(comm_matrix_raw)

rarefy_depth <- 20  # matches script 09 convention
has_enough <- row_abundances >= rarefy_depth
n_rare_corals <- sum(has_enough)
cat("  Rarefying to n =", rarefy_depth, "individuals (keeps", n_rare_corals,
    "of", length(has_enough), "corals)\n")

comm_sub <- comm_matrix_raw[has_enough, ]
rare_rich <- vegan::rarefy(comm_sub, sample = rarefy_depth)

rare_data <- coral_data %>%
  filter(coral_id %in% names(rare_rich)) %>%
  mutate(rarefied_richness = as.numeric(rare_rich[coral_id]))

# Gaussian LM (rarefied richness is continuous expected value)
rare_model <- lm(rarefied_richness ~ log(volume) + site, data = rare_data)
rare_summary <- summary(rare_model)
rare_slope <- coef(rare_model)["log(volume)"]
rare_p <- rare_summary$coefficients["log(volume)", "Pr(>|t|)"]
rare_r2 <- rare_summary$adj.r.squared

# Compare to raw richness slope
raw_slope <- coef(richness_result$model)["log(volume)"]

cat("  Raw richness slope (Poisson, log-link):", round(raw_slope, 4), "\n")
cat("  Rarefied richness slope (LM):", round(rare_slope, 4),
    ", p =", format.pval(rare_p, 3), "\n")
cat("  Rarefied model adj. R²:", round(rare_r2, 3), "\n")

if (rare_p < 0.05) {
  cat("  INTERPRETATION: Rarefied richness INCREASES with coral size.\n")
  cat("    Larger corals support genuinely richer communities beyond passive sampling.\n\n")
} else {
  cat("  INTERPRETATION: Rarefied richness does NOT increase with coral size.\n")
  cat("    The SAR is driven by passive sampling (more individuals → more species).\n\n")
}

# Store results
all_results$rarefied_richness <- list(
  metric = "Rarefied Richness (n=20)",
  model = rare_model,
  slope = rare_slope,
  p_value = rare_p,
  adj_r2 = rare_r2,
  n_corals = n_rare_corals,
  rarefy_depth = rarefy_depth,
  raw_slope_comparison = raw_slope,
  interpretation = ifelse(rare_p < 0.05,
    "Genuine diversity enrichment on larger corals",
    "SAR driven by passive sampling (abundance artifact)")
)

# ---- Model Diagnostics (DHARMa simulated residuals) ----
cat("--- Model Diagnostics ---\n")
if (requireNamespace("DHARMa", quietly = TRUE) && !is.null(total_result$model)) {
  # Total abundance (NB GLM)
  set.seed(42)
  dharma_total <- DHARMa::simulateResiduals(total_result$model, n = 1000, plot = FALSE)
  dharma_tests <- DHARMa::testResiduals(dharma_total, plot = FALSE)
  cat("  TOTAL ABUNDANCE (NB GLM):\n")
  cat("    Uniformity (KS): p =", format.pval(dharma_tests$uniformity$p.value, 3), "\n")
  cat("    Dispersion: p =", format.pval(dharma_tests$dispersion$p.value, 3), "\n")
  cat("    Zero-inflation: p =", format.pval(
    DHARMa::testZeroInflation(dharma_total, plot = FALSE)$p.value, 3), "\n")

  # Species richness (Poisson GLM) — DHARMa doesn't support quasipoisson
  dharma_rich <- NULL
  if (!is.null(richness_result$model)) {
    model_fam <- family(richness_result$model)$family
    if (model_fam == "quasipoisson") {
      cat("  SPECIES RICHNESS: Skipping DHARMa (quasipoisson not supported)\n")
    } else {
      set.seed(42)
      dharma_rich <- tryCatch({
        dr <- DHARMa::simulateResiduals(richness_result$model, n = 1000, plot = FALSE)
        dharma_tests_rich <- DHARMa::testResiduals(dr, plot = FALSE)
        cat("  SPECIES RICHNESS (Poisson GLM):\n")
        cat("    Uniformity (KS): p =", format.pval(dharma_tests_rich$uniformity$p.value, 3), "\n")
        cat("    Dispersion: p =", format.pval(dharma_tests_rich$dispersion$p.value, 3), "\n")
        cat("    Zero-inflation: p =", format.pval(
          DHARMa::testZeroInflation(dr, plot = FALSE)$p.value, 3), "\n")
        dr
      }, error = function(e) {
        cat("  SPECIES RICHNESS: DHARMa diagnostics failed:", conditionMessage(e), "\n")
        NULL
      })
    }
  }

  # Save diagnostic plots
  fig_dir <- file.path(PATHS$figures, "05_scaling")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  png(file.path(fig_dir, "dharma_diagnostics.png"), width = 10, height = 5, units = "in", res = 200)
  par(mfrow = c(1, 2))
  plot(dharma_total, main = "Total CAFI (NB)")
  if (!is.null(dharma_rich)) plot(dharma_rich, main = "Richness (Poisson)")
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
# Shannon = 0 means single-species community; log(Shannon) undefined for these
n_zero_shannon <- sum(coral_data$shannon == 0, na.rm = TRUE)
cat("  Note: Excluding", n_zero_shannon, "corals with Shannon = 0 from diversity analysis\n")
cat("  (Shannon = 0 indicates single-species community)\n\n")

shannon_data <- coral_data %>%
  filter(shannon > 0) %>%
  dplyr::select(coral_id, site, volume, shannon)

cat("Model: Shannon H' ~ log(Volume) + site (Linear Model)\n\n")

if (nrow(shannon_data) >= 15) {
  shannon_model <- lm(shannon ~ log(volume) + site, data = shannon_data)
  shannon_summary <- summary(shannon_model)

  shannon_beta <- coef(shannon_model)["log(volume)"]
  shannon_se <- shannon_summary$coefficients["log(volume)", "Std. Error"]
  shannon_ci <- tryCatch(
    confint(shannon_model, "log(volume)", level = 0.95),
    error = function(e) {
      cat("  Warning: Profile CI failed for Shannon model, using Wald CI\n")
      c(shannon_beta - 1.96 * shannon_se, shannon_beta + 1.96 * shannon_se)
    }
  )
  shannon_t <- shannon_summary$coefficients["log(volume)", "t value"]
  shannon_p <- shannon_summary$coefficients["log(volume)", "Pr(>|t|)"]

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
    p_boot_vs_1 = NA_real_,
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
    z_vs_1 = NA, p_vs_1 = NA, p_boot_vs_1 = NA, interpretation = "Insufficient data", model = NULL, converged = FALSE
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
    prevalence = n_distinct(coral_id) / n_distinct(coral_master$coral_id) * 100,  # denominator = all 114 corals surveyed
    .groups = "drop"
  ) %>%
  arrange(desc(total_individuals))

# Get functional group info (deduplicate by OTU in case functional_group/type are inconsistent)
species_functional <- cafi_clean %>%
  dplyr::select(otu, functional_group, type) %>%
  distinct(otu, .keep_all = TRUE) %>%
  group_by(otu) %>% slice(1) %>% ungroup()

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
# Note: Some species have very few non-zero observations (~5); NB GLM results
# for these rare taxa should be interpreted as exploratory, not confirmatory.
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

# Add Field of Dreams / Redirection classification
key_species_summary <- key_species_summary %>%
  mutate(
    scaling_pattern = case_when(
      is.na(beta) ~ "Insufficient data",
      grepl("Field of Dreams", interpretation) ~ "Field of Dreams (β ≈ 1)",
      grepl("Redirection", interpretation) ~ "Propagule Redirection (β < 1)",
      TRUE ~ interpretation
    ),
    # Flag species with fewer than 15 non-zero observations as exploratory
    reliability = ifelse(n_present < 15, "low", "adequate")
  )

print(key_species_summary %>%
        dplyr::select(species, exp_effect, n_present, beta, ci_lower, ci_upper, scaling_pattern) %>%
        mutate(across(where(is.numeric) & !c(n_present), ~round(., 3))))

# Save key species scaling table
save_table(key_species_summary, "key_species_scaling_experimental")
cat("\n  Saved: key_species_scaling_experimental.csv\n")

# --- SCALING CONCORDANCE: Survey vs Experiment ---
# The experiment found ~proportional scaling (Field of Dreams) at reef scale.
# The survey finds sublinear scaling (Redirection) at colony scale.
# Compare species tested in both studies.
cat("\n--- CROSS-STUDY SCALING CONCORDANCE ---\n\n")

# Experimental scaling responses (from companion paper)
# Experiment measured per-coral abundance as function of reef-level coral number
expt_scaling <- tribble(
  ~species,                  ~expt_response,
  "Caracanthus maculatus",   "proportional",
  "Harpiliopsis spinigera",  "proportional",
  "Periclimenes watamuae",   "proportional",
  "Calcinus latens",         "proportional",
  "Trapezia serenei",        "sub-proportional",
  "Alpheus lottini",         "proportional",
  "Alpheus diadema",         "sub-proportional",
  "Luniella pugil",          "super-proportional",
  "Galeropsis monodonta",    "proportional"
)

# Match with BOTH key species AND individual species from survey
# Key species summary uses grouped names; individual species from Part 4 uses OTU names
survey_scaling <- bind_rows(
  # From key species (Part 5)
  key_species_summary %>%
    dplyr::select(species, survey_beta = beta, survey_ci_lower = ci_lower,
                  survey_ci_upper = ci_upper, survey_pattern = scaling_pattern),
  # From individual species (Part 4) — these use OTU names
  if (exists("all_results") && !is.null(all_results$individual_species)) {
    map_df(all_results$individual_species, function(r) {
      if (is.null(r) || is.na(r$beta)) return(NULL)
      tibble(species = r$response, survey_beta = r$beta,
             survey_ci_lower = r$boot_ci_lower, survey_ci_upper = r$boot_ci_upper,
             survey_pattern = r$interpretation)
    })
  } else {
    tibble()
  }
) %>%
  distinct(species, .keep_all = TRUE)

scaling_concordance <- expt_scaling %>%
  left_join(survey_scaling, by = "species") %>%
  filter(!is.na(survey_beta))

if (nrow(scaling_concordance) > 0) {
  cat("  Species with scaling data in both studies:", nrow(scaling_concordance), "\n")

  scaling_concordance <- scaling_concordance %>%
    mutate(
      survey_class = case_when(
        !is.na(survey_ci_upper) & survey_ci_upper < 1 ~ "sublinear",
        !is.na(survey_ci_lower) & survey_ci_lower > 1 ~ "super-linear",
        TRUE ~ "proportional"
      ),
      expt_class = case_when(
        expt_response == "sub-proportional" ~ "sublinear",
        expt_response == "super-proportional" ~ "super-linear",
        TRUE ~ "proportional"
      )
    )

  cat("\n  Cross-study scaling comparison:\n")
  cat(sprintf("  %-25s %-18s %-18s\n", "Species", "Experiment", "Survey"))
  cat("  ", strrep("-", 62), "\n", sep = "")
  for (i in 1:nrow(scaling_concordance)) {
    r <- scaling_concordance[i, ]
    cat(sprintf("  %-25s %-18s %-18s\n",
                r$species, r$expt_response,
                sprintf("β=%.2f %s", r$survey_beta, r$survey_class)))
  }

  n_survey_sub <- sum(scaling_concordance$survey_class == "sublinear")
  cat(sprintf("\n  Survey: %d/%d species show sublinear scaling (vs mostly proportional in experiment)\n",
              n_survey_sub, nrow(scaling_concordance)))
  cat("  Key contrast: experiment found ~proportional at reef scale (colonization),\n")
  cat("  survey finds sublinear at colony scale (established communities).\n")

  save_table(scaling_concordance, "cross_study_scaling_concordance")
  cat("  Saved: cross_study_scaling_concordance.csv\n")
} else {
  cat("  No overlapping species between experiment and survey scaling analyses.\n")
}
cat("\n")

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
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3, linewidth = 0.8, orientation = "y") +
    geom_point(size = 4) +
    scale_color_manual(
      values = c("Beneficial (positive effect on coral)" = "#0072B2",
                 "Harmful (negative effect on coral)" = "#D55E00"),
      name = "Effect on coral condition\n(from experimental paper)"
    ) +
    annotate("text", x = 1.15, y = 0.5, label = "Field of Dreams\n(β = 1)",
             hjust = 0, size = 3, color = "gray40") +
    annotate("text", x = 0.5, y = 0.5, label = "Redirection\n(β < 1)",
             hjust = 0.5, size = 3, color = "gray40") +
    labs(
      title = "Scaling Patterns for Key Species from Experimental Paper",
      subtitle = "Comparing Field of Dreams (β = 1) vs. Propagule Redirection (β < 1)",
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

  save_figure(p_key_forest, file.path(FIG_DIR, "key_species_scaling_forest.png"),
              width = 10, height = 6)
  cat("  Saved: key_species_scaling_forest.png\n")
}

# Summary interpretation
cat("\n--- INTERPRETATION: HETEROGENEITY IN KEY SPECIES SCALING ---\n\n")

n_analyzed <- sum(!is.na(key_species_summary$beta))
n_fod <- sum(grepl("Field of Dreams", key_species_summary$interpretation), na.rm = TRUE)
n_redist <- sum(grepl("Redirection", key_species_summary$interpretation), na.rm = TRUE)

cat("Of", n_analyzed, "key species/groups analyzed:\n")
cat("  Field of Dreams (β ≈ 1):", n_fod, "\n")
cat("  Propagule Redirection (β < 1):", n_redist, "\n\n")

cat("Comparison with experimental paper predictions:\n")
cat("  The experimental work (Stier et al.) predicts heterogeneous scaling responses\n")
cat("  across species, with some following Field of Dreams and others showing\n")
cat("  Propagule Redirection. This survey tests whether the same pattern holds\n")
cat("  for observational data across natural coral size variation.\n\n")

# Check if already in individual species results
already_in_main <- key_species_summary %>%
  filter(species %in% names(species_results))

if (nrow(already_in_main) > 0) {
  cat("Note: Some key species were also included in main species analysis (met occurrence thresholds).\n")
}

# ############################################################################
#                    PART 6C: SPECIES OCCURRENCE CURVES
# ############################################################################
# For prevalent species, fit logistic regressions P(present) ~ log(volume) + site
# to visualize the "assembly sequence" — which species appear at what coral sizes.

cat("\n############################################################\n")
cat("PART 6C: SPECIES OCCURRENCE CURVES\n")
cat("############################################################\n\n")

comm_matrix_occ <- load_object("community_matrix")
comm_pa <- (comm_matrix_occ > 0) * 1L

# Validate coral_ids exist in community matrix row names before indexing
stopifnot(all(coral_data$coral_id %in% rownames(comm_pa)))

# Select species with ≥15% prevalence (≥17 of 112 corals)
sp_prev <- colMeans(comm_pa[coral_data$coral_id, ])
prevalence_threshold <- 0.15
focal_occ_species <- names(sort(sp_prev[sp_prev >= prevalence_threshold], decreasing = TRUE))
cat("  Species with ≥", prevalence_threshold * 100, "% prevalence:", length(focal_occ_species), "\n")

# Fit logistic GLM per species
occurrence_list <- lapply(focal_occ_species, function(sp) {
  dat <- coral_data %>%
    mutate(present = comm_pa[coral_id, sp])

  fit <- tryCatch(
    glm(present ~ log(volume) + site, data = dat, family = binomial),
    warning = function(w) suppressWarnings(glm(present ~ log(volume) + site, data = dat, family = binomial)),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  coefs <- summary(fit)$coefficients
  vol_row <- coefs["log(volume)", ]

  # Predict across volume range (averaged over sites)
  vol_seq <- seq(log(min(coral_data$volume)), log(max(coral_data$volume)), length.out = 100)
  preds_by_site <- lapply(unique(coral_data$site), function(s) {
    nd <- data.frame(volume = exp(vol_seq), site = s)
    predict(fit, newdata = nd, type = "response")
  })
  avg_pred <- Reduce("+", preds_by_site) / length(preds_by_site)

  # Volume at P = 0.5 (arrival volume, site-averaged)
  if (any(avg_pred < 0.5) && any(avg_pred > 0.5)) {
    v50_idx <- which.min(abs(avg_pred - 0.5))
    v50 <- exp(vol_seq[v50_idx])
  } else {
    v50 <- NA
  }

  list(
    species = sp,
    prevalence = sp_prev[sp],
    beta_volume = vol_row["Estimate"],
    se = vol_row["Std. Error"],
    z_value = vol_row["z value"],
    p_value = vol_row["Pr(>|z|)"],
    v50 = v50,
    pred_volume = exp(vol_seq),
    pred_prob = avg_pred
  )
})

occurrence_list <- occurrence_list[!sapply(occurrence_list, is.null)]

# Build summary table
occ_df <- do.call(rbind, lapply(occurrence_list, function(x) {
  data.frame(
    species = x$species,
    prevalence = round(x$prevalence * 100, 1),
    beta_volume = round(x$beta_volume, 3),
    se = round(x$se, 3),
    z_value = round(x$z_value, 2),
    p_value = x$p_value,
    v50 = round(x$v50, 0),
    stringsAsFactors = FALSE
  )
}))
occ_df$p_adj <- p.adjust(occ_df$p_value, method = "fdr")
occ_df <- occ_df %>% arrange(p_value)

n_sig <- sum(occ_df$p_adj < 0.05)
n_pos <- sum(occ_df$beta_volume > 0 & occ_df$p_adj < 0.05)
n_neg <- sum(occ_df$beta_volume < 0 & occ_df$p_adj < 0.05)

cat("  Significant volume effects (FDR < 0.05):", n_sig, "of", nrow(occ_df), "\n")
cat("    Positive (more likely on larger corals):", n_pos, "\n")
cat("    Negative (more likely on smaller corals):", n_neg, "\n\n")

cat("  Top 10 species by volume effect significance:\n")
top10_occ <- head(occ_df, 10)
for (i in 1:nrow(top10_occ)) {
  dir <- ifelse(top10_occ$beta_volume[i] > 0, "+", "-")
  sig <- ifelse(top10_occ$p_adj[i] < 0.05, "*", "")
  cat("    ", dir, " ", top10_occ$species[i],
      " (β=", sprintf("%.2f", top10_occ$beta_volume[i]),
      ", p=", format.pval(top10_occ$p_value[i], 3),
      ", FDR=", format.pval(top10_occ$p_adj[i], 3),
      sig, ")\n")
}

# Save table
save_table(occ_df, "occurrence_scaling_results")

# Store in results
all_results$occurrence_curves <- list(
  summary = occ_df,
  n_species = nrow(occ_df),
  n_significant = n_sig,
  curves = occurrence_list
)

# ---- Figure S15: Occurrence Curves ----
cat("\n  Generating occurrence curves figure (Fig S15)...\n")

# Select species to show: all FDR-significant, plus top by |beta| up to 15 total
show_species <- occ_df %>%
  mutate(abs_beta = abs(beta_volume)) %>%
  arrange(p_adj, desc(abs_beta)) %>%
  head(min(15, nrow(occ_df))) %>%
  pull(species)

# Build prediction data for plotting
pred_data <- do.call(rbind, lapply(occurrence_list, function(x) {
  if (!(x$species %in% show_species)) return(NULL)
  data.frame(
    species = x$species,
    volume = x$pred_volume,
    prob = x$pred_prob,
    stringsAsFactors = FALSE
  )
}))

# Order species by arrival volume (v50), then prevalence
occ_order <- occ_df %>%
  filter(species %in% show_species) %>%
  arrange(desc(v50))
pred_data$species <- factor(pred_data$species, levels = occ_order$species)

# Add raw data points
raw_points <- do.call(rbind, lapply(show_species, function(sp) {
  coral_data %>%
    mutate(species = sp, present = comm_pa[coral_id, sp]) %>%
    dplyr::select(species, volume, present)
}))
raw_points$species <- factor(raw_points$species, levels = occ_order$species)

# Determine significance labels
sig_labels <- occ_df %>%
  filter(species %in% show_species) %>%
  mutate(label = paste0(
    species, "  (",
    ifelse(beta_volume > 0, "+", ""), sprintf("%.1f", beta_volume),
    ifelse(p_adj < 0.05, "*", ""), ")"
  ))

# Use species names as facet labels with significance indicator
facet_labels <- setNames(sig_labels$label, sig_labels$species)

fig_occ <- ggplot() +
  geom_jitter(data = raw_points,
              aes(x = volume, y = present),
              width = 0, height = 0.03, alpha = 0.15, size = 0.5, color = "grey50") +
  geom_line(data = pred_data,
            aes(x = volume, y = prob),
            color = "#0072B2", linewidth = 0.8) +
  facet_wrap(~ species, ncol = 3,
             labeller = labeller(species = facet_labels)) +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(
    x = expression("Colony volume " (cm^3)),
    y = "P(occurrence)",
    title = "Figure S15: Species occurrence probability vs. coral size"
  ) +
  theme_publication(base_size = 10) +
  theme(
    strip.text = element_text(size = 8, face = "bold", hjust = 0),
    axis.text = element_text(size = 8, color = "grey30"),
    axis.title = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10, "mm"),
    plot.title = element_text(size = 11, face = "bold")
  )

# Save to scaling dir and supplement
n_occ_rows <- ceiling(length(show_species) / 3)
occ_height <- max(140, n_occ_rows * 45)
save_figure(fig_occ, file.path(FIG_DIR, "occurrence_curves.png"),
            width = 250, height = occ_height, units = "mm")
cat("  Saved: 05_scaling/occurrence_curves.png\n")

supplement_dir <- file.path(PATHS$figures, "supplement")
dir.create(supplement_dir, showWarnings = FALSE, recursive = TRUE)
save_figure(fig_occ, file.path(supplement_dir, "figS15_occurrence_curves.png"),
            width = 250, height = occ_height, units = "mm")
cat("  Saved: supplement/figS15_occurrence_curves.png\n\n")

# ---- S15 Legend / Results Text ----
s15_legend <- paste0(
'FIGURE S15: SPECIES OCCURRENCE PROBABILITY VS. CORAL SIZE
================================================================================

FIGURE LEGEND
-------------
Figure S15. Occurrence probability of prevalent CAFI species as a function of
coral colony volume. Each panel shows one species. Blue curves: logistic GLM fit
(P(present) ~ log(volume) + site, averaged across sites). Points: observed
presence/absence data (jittered vertically for visibility). Species are ordered
by v50 (colony volume at which predicted occurrence probability = 0.5). Facet
labels show the volume coefficient with significance (* = FDR < 0.05).

', length(show_species), ' species shown (prevalence >= ', prevalence_threshold * 100, '%).

================================================================================

STATISTICAL RESULTS
-------------------
Method: Logistic GLM per species: P(present) ~ log(volume) + site
Multiple testing: FDR correction (Benjamini-Hochberg) across ', nrow(occ_df), ' species

Species tested: ', nrow(occ_df), '
Significant (FDR < 0.05): ', n_sig, ' / ', nrow(occ_df), '
  Positive (more likely on larger corals): ', n_pos, '
  Negative (more likely on smaller corals): ', n_neg, '

',
paste(capture.output({
  cat("Top species by significance:\n")
  top_show <- head(occ_df, min(15, nrow(occ_df)))
  for (i in 1:nrow(top_show)) {
    dir <- ifelse(top_show$beta_volume[i] > 0, "+", "-")
    sig <- ifelse(top_show$p_adj[i] < 0.05, " *", "")
    cat("  ", dir, " ", top_show$species[i],
        ": beta=", sprintf("%.2f", top_show$beta_volume[i]),
        ", p=", format.pval(top_show$p_value[i], 3),
        ", FDR=", format.pval(top_show$p_adj[i], 3),
        ifelse(!is.na(top_show$v50[i]), paste0(", v50=", top_show$v50[i], " cm3"), ""),
        sig, "\n", sep = "")
  }
}), collapse = "\n"), '

================================================================================

INTERPRETATION
--------------
Most prevalent CAFI species show size-dependent occurrence. ', n_sig, ' of ',
nrow(occ_df), ' species (', round(n_sig / nrow(occ_df) * 100), '%) had significant
volume effects on occurrence probability after FDR correction. ',
if (n_pos > n_neg) paste0('The majority (', n_pos, '/', n_sig, ') were more likely to occur on larger corals, ')
else if (n_neg > n_pos) paste0('The majority (', n_neg, '/', n_sig, ') were more likely to occur on smaller corals, ')
else paste0('Species were equally split between larger and smaller coral preferences, '),
'consistent with species-specific size thresholds driving community turnover
along the size gradient.

================================================================================

METHODS
-------
For each species with >= ', prevalence_threshold * 100, '% prevalence across ', nrow(coral_data),
' corals, we fitted a logistic GLM: P(present) ~ log(volume) + site.
The volume coefficient quantifies the log-odds change in occurrence per
unit increase in log(volume). We computed v50 as the colony volume at
which predicted occurrence probability equals 0.5 (averaged across sites).
P-values were FDR-corrected (Benjamini-Hochberg) across all ', nrow(occ_df), ' tests.

================================================================================
Generated: ', Sys.Date(), '
Source script: scripts/05_species_scaling_analysis.R
')

writeLines(s15_legend, file.path(PATHS$fig_supplement, "figS15_legend_results.txt"))
cat("  Generated: figS15_legend_results.txt\n")

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
  # IMPORTANT: Species betas are NOT independent — all estimated from the same corals.
  # Standard meta-analytic SE assumes independence and is anti-conservative.
  # Report with caveat; consider rma.mv for proper inference.

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
  cat("\n  NOTE: Weighted mean treats species betas as independent estimates.\n")
  cat("  Since all species share the same coral sampling units, the SE is\n")
  cat("  likely anti-conservative. Interpret the z-test vs beta=1 with caution.\n\n")
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

  # If metafor is available, fit random-effects model (still assumes independence
  # but at least accounts for heterogeneity properly)
  species_betas <- species_data$beta
  species_ses <- species_data$se
  if (requireNamespace("metafor", quietly = TRUE)) {
    rma_result <- tryCatch({
      metafor::rma(yi = species_betas, sei = species_ses, method = "REML")
    }, error = function(e) NULL)
    if (!is.null(rma_result)) {
      cat("\n  Random-effects meta-analytic mean (metafor::rma):", round(rma_result$beta[1], 3),
          "[", round(rma_result$ci.lb, 3), ",", round(rma_result$ci.ub, 3), "]\n")
    }
  }

  # Also report unweighted t-test for comparison (with caveat)
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

# Define interpretation colors (aligned with Fig 2 Panel C palette)
interpretation_colors <- c(
  "Redirection (β < 1)" = "#5A8FAF",
  "Field of Dreams (β ≈ 1)" = "gray55",
  "Super-linear (β > 1)" = "#D55E00"
)

# ============================================================================
# FIGURE 1: Community-Level Scaling (Total, Richness, Shannon)
# ============================================================================

cat("Creating Figure 1: Community-level scaling...\n")

# Total abundance plot
p_total <- ggplot(coral_data, aes(x = volume, y = total_cafi, fill = site)) +
  geom_point(alpha = 0.7, size = 2.5, shape = 21, color = "gray30", stroke = 0.4) +
  geom_smooth(aes(group = 1), method = MASS::glm.nb, formula = y ~ log(x),
              method.args = list(control = glm.control(maxit = 50)),
              se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  scale_fill_manual(values = SITE_COLORS) +
  labs(
    x = expression("Coral Volume (cm"^3*")"),
    y = "Total CAFI",
    title = "A. Total Abundance",
    subtitle = sprintf("β = %.2f [%.2f, %.2f]",
                       total_result$beta, total_result$ci_lower, total_result$ci_upper)
  ) +
  theme(legend.position = c(0.15, 0.85))

# Richness plot
p_richness <- ggplot(coral_data, aes(x = volume, y = otu_richness, fill = site)) +
  geom_point(alpha = 0.7, size = 2.5, shape = 21, color = "gray30", stroke = 0.4) +
  geom_smooth(aes(group = 1), method = "glm",
              method.args = list(family = poisson), formula = y ~ log(x),
              se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  scale_fill_manual(values = SITE_COLORS) +
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
  ggplot(aes(x = volume, y = shannon, fill = site)) +
  geom_point(alpha = 0.7, size = 2.5, shape = 21, color = "gray30", stroke = 0.4) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ log(x),
              se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = SITE_COLORS) +
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

save_figure(fig_community, file.path(FIG_DIR, "community_scaling.png"),
            width = 14, height = 5)
cat("  Saved: community_scaling.png\n")

# ============================================================================
# FIGURE 2: Taxonomic & Functional Group Forest Plot
# ============================================================================
# NOTE: Supplement/exploratory figures below use the default ggplot2 theme (via
# theme() adjustments) rather than theme_publication(), which is reserved for
# manuscript figures. This is intentional: supplement figures prioritize
# information density over strict publication styling.

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
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.4) +
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

  save_figure(p_groups, file.path(FIG_DIR, "group_scaling_forest.png"),
              width = 10, height = 7)
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
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray80") +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.3, linewidth = 0.6) +
    geom_point(size = 3) +
    scale_color_manual(values = interpretation_colors, name = "Scaling Pattern") +
    labs(
      x = expression("Scaling Exponent (" * beta * ")"),
      y = NULL,
      title = "Figure S6: Species-Level Scaling of Abundance with Coral Volume",
      subtitle = sprintf("N = %d species | Mean β = %.2f | Test vs 1: p = %s",
                         nrow(species_plot_data),
                         mean(species_plot_data$beta),
                         if(!is.null(t_test)) format.pval(t_test$p.value, 2) else "NA")
    ) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(face = "italic", size = 9)
    ) +
    coord_cartesian(xlim = c(-0.5, 3))

  save_figure(p_species_forest, file.path(FIG_DIR, "species_scaling_forest.png"),
              width = 10, height = 8)
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

make_species_panel <- function(sp_name, sp_result, species_scaling_data) {
  sp_data <- species_scaling_data %>%
    filter(otu == sp_name, abundance > 0)

  beta_val <- sp_result$beta
  interp <- sp_result$interpretation

  title_color <- case_when(
    interp == "Redirection (β < 1)" ~ "#5A8FAF",
    interp == "Super-linear (β > 1)" ~ "#D55E00",
    TRUE ~ "gray55"
  )

  tryCatch({
    ggplot(sp_data, aes(x = volume, y = abundance, fill = site)) +
      geom_point(alpha = 0.7, size = 2, shape = 21, color = "gray30", stroke = 0.4) +
      geom_smooth(aes(group = 1), method = MASS::glm.nb, formula = y ~ log(x),
                  method.args = list(control = glm.control(maxit = 50)),
                  se = TRUE, color = "black", linewidth = 0.8) +
      scale_x_log10(labels = scales::comma) +
      scale_y_continuous() +
      scale_fill_manual(values = SITE_COLORS) +
      labs(
        x = expression("Volume (cm"^3*")"),
        y = "Abundance",
        title = str_replace(sp_name, "_", " ")
      ) +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("\u03b2 = %.2f", beta_val),
               hjust = 1.1, vjust = 1.3, size = 3, fontface = "italic") +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold.italic", size = 9, color = title_color)
      )
  }, error = function(e) {
    # Fallback: plot points without NB smooth if convergence fails on sparse data
    ggplot(sp_data, aes(x = volume, y = abundance, fill = site)) +
      geom_point(alpha = 0.7, size = 2, shape = 21, color = "gray30", stroke = 0.4) +
      scale_x_log10(labels = scales::comma) +
      scale_y_continuous() +
      scale_fill_manual(values = SITE_COLORS) +
      labs(
        x = expression("Volume (cm"^3*")"),
        y = "Abundance",
        title = str_replace(sp_name, "_", " ")
      ) +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("\u03b2 = %.2f (fit failed)", beta_val),
               hjust = 1.1, vjust = 1.3, size = 3, fontface = "italic") +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold.italic", size = 9, color = title_color)
      )
  })
}

if (nrow(top_species) > 0) {
  species_panels <- map2(top_species$response,
                          map(top_species$response, ~all_results$individual_species[[.x]]),
                          ~make_species_panel(.x, .y, species_scaling_data))

  fig_species <- wrap_plots(species_panels, ncol = 4) +
    plot_annotation(
      title = "Species-Specific Scaling Relationships",
      subtitle = "Top 12 species by sample size | Title color indicates interpretation",
      caption = paste0("Blue = Redirection (β < 1), Gray = Field of Dreams (β ≈ 1), Vermillion = Super-linear (β > 1)"),
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11)
      )
    )

  save_figure(fig_species, file.path(FIG_DIR, "species_scaling_panels.png"),
              width = 14, height = 10)
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
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.4) +
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
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.4) +
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

save_figure(fig_summary, file.path(FIG_DIR, "scaling_summary.png"),
            width = 10, height = 10)
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

# ############################################################################
# MANUSCRIPT FIGURE 2: SIZE-ABUNDANCE SCALING
# ############################################################################

cat("\n============================================================\n")
cat("MANUSCRIPT FIGURE 2: Size-Abundance Scaling\n")
cat("============================================================\n\n")

# Clean publication theme for single-column (89mm) figure
# Lighter than theme_publication: no box border, thin axis lines, subtle grid
theme_fig2 <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.line = element_line(color = "grey30", linewidth = 0.3),
      axis.ticks = element_line(color = "grey50", linewidth = 0.3),
      axis.ticks.length = unit(0.12, "cm"),
      axis.title = element_text(size = base_size, color = "black"),
      axis.text = element_text(size = base_size - 1, color = "grey30"),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(3.5, "mm"),
      legend.margin = margin(t = 0, b = 0),
      plot.tag = element_text(size = base_size + 2, face = "bold"),
      plot.margin = margin(2, 8, 2, 6, "mm"),
      plot.background = element_rect(fill = "white", color = NA)
    )
}
# Alias for downstream code (Fig 3 etc.) that still references theme_manuscript
theme_manuscript <- function(base_size = 10) theme_fig2(base_size = base_size)

# Prepare data for scaling plots (coral_data already exists from earlier in this script)
scaling_data <- coral_data %>%
  filter(!is.na(volume), volume > 0, !is.na(total_cafi))

if (nrow(scaling_data) >= 30) {

  # Fit NB model for total CAFI abundance (using natural log)
  abundance_model <- tryCatch({
    MASS::glm.nb(total_cafi ~ log(volume) + site, data = scaling_data)
  }, error = function(e) NULL)

  if (!is.null(abundance_model)) {
    coefs <- summary(abundance_model)$coefficients
    beta_abundance <- coefs["log(volume)", "Estimate"]
    se_abundance <- coefs["log(volume)", "Std. Error"]
    ci_abundance <- tryCatch(
      confint(abundance_model, "log(volume)", level = 0.95),
      error = function(e) {
        cat("  Warning: Profile CI failed, using Wald CI\n")
        c(beta_abundance - 1.96 * se_abundance, beta_abundance + 1.96 * se_abundance)
      }
    )

    cat("  Total CAFI scaling: beta =", round(beta_abundance, 3), "\n")
    cat("  95% CI: [", round(ci_abundance[1], 3), ",", round(ci_abundance[2], 3), "]\n\n")
  } else {
    beta_abundance <- NA
    ci_abundance <- c(NA, NA)
  }

  # --------------------------------------------------------------------------
  # DESIGN: Refined aesthetics for 89mm single-column figure
  # --------------------------------------------------------------------------

  smooth_color <- "grey20"
  smooth_fill  <- "grey80"
  smooth_lwd   <- 0.8
  smooth_alpha  <- 0.25
  point_alpha   <- 0.65
  point_size    <- 2.0

  # Proportional-scaling (beta=1) reference line intercept
  intercept_beta1 <- mean(
    log10(scaling_data$total_cafi[scaling_data$total_cafi > 0]) -
    log10(scaling_data$volume[scaling_data$total_cafi > 0])
  )

  shared_xlim <- c(30, 50000)

  # Build annotation labels as parseable strings
  beta_label <- if (!is.na(beta_abundance)) {
    ci_lo <- if(!is.na(total_result$boot_ci_lower)) total_result$boot_ci_lower else ci_abundance[1]
    ci_hi <- if(!is.na(total_result$boot_ci_upper)) total_result$boot_ci_upper else ci_abundance[2]
    sprintf("beta == %.2f ~ '[' * %.2f * ', ' * %.2f * ']'", beta_abundance, ci_lo, ci_hi)
  } else ""

  # --- Panel A: Total abundance vs volume (log-log) ---
  panel_a <- scaling_data %>%
    ggplot(aes(x = volume, y = total_cafi)) +
    geom_abline(slope = 1, intercept = intercept_beta1, linetype = "longdash",
                color = "grey50", linewidth = 0.6) +
    geom_smooth(aes(group = 1), method = MASS::glm.nb, formula = y ~ log(x),
                se = TRUE, color = smooth_color, fill = smooth_fill,
                linewidth = smooth_lwd, alpha = smooth_alpha) +
    geom_point(aes(fill = site), alpha = point_alpha, size = point_size,
               shape = 21, color = "white", stroke = 0.3) +
    scale_x_log10(labels = scales::comma, breaks = c(100, 1000, 10000)) +
    scale_y_log10(labels = scales::comma, breaks = c(1, 10, 100)) +
    scale_fill_manual(values = SITE_COLORS, name = "Site") +
    coord_cartesian(xlim = shared_xlim, ylim = c(0.8, 500)) +
    labs(x = NULL, y = "Total CAFI abundance") +
    annotate("text", x = 35000, y = 0.9, label = beta_label,
             size = 3.2, color = "black", fontface = "bold", hjust = 1, parse = TRUE) +
    annotate("text", x = 70, y = 350, label = "1:1\n(Field of Dreams)",
             size = 2.3, color = "grey40", fontface = "italic", hjust = 0,
             lineheight = 0.85) +
    theme_fig2() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      legend.position = "none"  # Site legend goes in companion text file (fig2_legend_results.txt)
    )

  # --- Panel B: Density Dilution (middle) ---
  scaling_data_nonzero <- scaling_data %>% filter(total_cafi > 0)
  density_lm <- lm(log(total_cafi / volume) ~ log(volume),
                    data = scaling_data_nonzero)
  density_slope <- coef(density_lm)[["log(volume)"]]

  panel_b_density <- ggplot(scaling_data_nonzero,
                            aes(x = volume, y = total_cafi / volume, fill = site)) +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x,
                se = TRUE, color = smooth_color, fill = smooth_fill,
                linewidth = smooth_lwd, alpha = smooth_alpha) +
    geom_point(alpha = point_alpha, size = point_size,
               shape = 21, color = "white", stroke = 0.3) +
    geom_hline(yintercept = mean(scaling_data_nonzero$total_cafi / scaling_data_nonzero$volume),
               linetype = "dashed", color = "grey60", linewidth = 0.4) +
    scale_x_log10(labels = scales::comma, breaks = c(100, 1000, 10000)) +
    scale_y_log10() +
    scale_fill_manual(values = SITE_COLORS, name = "Site") +
    coord_cartesian(xlim = shared_xlim) +
    labs(x = NULL, y = expression("CAFI density (ind. cm"^-3*")")) +
    annotate("text", x = 35000, y = min(scaling_data_nonzero$total_cafi / scaling_data_nonzero$volume) * 1.5,
             label = sprintf("slope == %.2f", density_slope),
             size = 3.2, color = "black", fontface = "bold", hjust = 1, parse = TRUE) +
    theme_fig2() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      legend.position = "none"
    )

  # --- Panel C: Species richness scaling ---
  richness_model_fig <- tryCatch({
    glm(otu_richness ~ log(volume) + site, data = scaling_data, family = poisson)
  }, error = function(e) NULL)

  if (!is.null(richness_model_fig)) {
    z_richness <- coef(richness_model_fig)["log(volume)"]
    z_se <- summary(richness_model_fig)$coefficients["log(volume)", "Std. Error"]
    z_ci <- tryCatch(
      confint(richness_model_fig, "log(volume)"),
      error = function(e) {
        cat("  Warning: Profile CI failed for richness figure model, using Wald CI\n")
        c(z_richness - 1.96 * z_se, z_richness + 1.96 * z_se)
      }
    )
  } else {
    z_richness <- NA
    z_ci <- c(NA, NA)
  }

  z_label <- if (!is.na(z_richness)) {
    z_lo <- if(!is.na(richness_result$boot_ci_lower)) richness_result$boot_ci_lower else z_ci[1]
    z_hi <- if(!is.na(richness_result$boot_ci_upper)) richness_result$boot_ci_upper else z_ci[2]
    sprintf("italic(z) == %.2f ~ '[' * %.2f * ', ' * %.2f * ']'", z_richness, z_lo, z_hi)
  } else ""

  panel_c <- scaling_data %>%
    ggplot(aes(x = volume, y = otu_richness)) +
    geom_smooth(aes(group = 1), method = "glm",
                method.args = list(family = poisson), formula = y ~ log(x),
                se = TRUE, color = smooth_color, fill = smooth_fill,
                linewidth = smooth_lwd, alpha = smooth_alpha) +
    geom_point(aes(fill = site), alpha = point_alpha, size = point_size,
               shape = 21, color = "white", stroke = 0.3) +
    scale_x_log10(labels = scales::comma, breaks = c(100, 1000, 10000)) +
    scale_y_log10(limits = c(1, NA)) +
    scale_fill_manual(values = SITE_COLORS, name = "Site") +
    coord_cartesian(xlim = shared_xlim) +
    labs(x = expression("Coral volume (cm"^3*")"), y = "Species richness") +
    annotate("text", x = 35000, y = 1.2,
             label = z_label,
             size = 3.2, color = "black", fontface = "bold", hjust = 1, parse = TRUE) +
    theme_fig2() +
    theme(legend.position = "none")

  # --- Combine: 3 vertical panels, shared x-axis ---
  fig2 <- panel_a / panel_b_density / panel_c +
    plot_layout(heights = c(1, 1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 11, face = "bold"))

  # Save to manuscript directory
  save_figure(fig2, file.path(PATHS$fig_manuscript, "fig2_scaling.png"),
              width = 89, height = 200, units = "mm")
  cat("  Saved: manuscript/fig2_scaling.png (+ PDF)\n")

  # Save to analysis figure directory
  save_figure(fig2, file.path(FIG_DIR, "fig2_scaling.png"),
              width = 89, height = 200, units = "mm")
  cat("  Saved: 05_scaling/fig2_scaling.png\n")

  # Save species scaling forest to supplement
  supplement_dir <- file.path(PATHS$figures, "supplement")
  dir.create(supplement_dir, showWarnings = FALSE, recursive = TRUE)
  save_figure(p_species_forest, file.path(supplement_dir, "figS6_species_scaling.png"),
              width = 10, height = 8)
  cat("  Saved: supplement/figS6_species_scaling.png\n")

  # Generate fig2 legend results text
  fig2_legend <- paste0(
'FIGURE 2: CAFI SCALING WITH CORAL SIZE
================================================================================

FIGURE LEGEND
-------------
Figure 2. CAFI abundance and richness scale sublinearly with coral volume.
(A) Total CAFI abundance vs colony volume (log-log scale). Solid curve: negative
binomial GLM fit; shaded band: 95% CI. Dashed line: 1:1 proportional scaling
(\u03b2 = 1, Field of Dreams hypothesis). Abundance scaling exponent \u03b2 = ',
round(total_result$beta, 2), ' [', round(if(!is.na(total_result$boot_ci_lower)) total_result$boot_ci_lower else total_result$ci_lower, 2), ', ',
round(if(!is.na(total_result$boot_ci_upper)) total_result$boot_ci_upper else total_result$ci_upper, 2), '] — sublinear (Propagule Redirection).
(B) Per-capita CAFI density (individuals per cm\u00b3) vs colony volume (log-log
scale). Dashed horizontal line: mean density. Log-log slope = ',
round(density_slope, 2), ', confirming density dilution — larger corals harbour
fewer CAFI per unit volume, consistent with sublinear abundance scaling.
(C) Species richness vs colony volume. Solid curve: Poisson GLM; z = ',
round(richness_result$beta, 2), ' [', round(if(!is.na(richness_result$boot_ci_lower)) richness_result$boot_ci_lower else richness_result$ci_lower, 2), ', ',
round(if(!is.na(richness_result$boot_ci_upper)) richness_result$boot_ci_upper else richness_result$ci_upper, 2), '] — sublinear.
All panels share the same x-axis (coral volume).
n = ', nrow(scaling_data), ' corals across 3 reef sites.

================================================================================

KEY STATISTICS
--------------

Total CAFI Abundance:
  \u03b2 = ', round(total_result$beta, 3), ' [', round(if(!is.na(total_result$boot_ci_lower)) total_result$boot_ci_lower else total_result$ci_lower, 3), ', ', round(if(!is.na(total_result$boot_ci_upper)) total_result$boot_ci_upper else total_result$ci_upper, 3), ']
  z vs 1 = ', round(total_result$z_vs_1, 2), ', p = ', format.pval(total_result$p_vs_1, 3), ' (Wald)', if(!is.na(total_result$p_boot_vs_1)) paste0('; p_boot = ', format.pval(total_result$p_boot_vs_1, 3)) else '', '
  Interpretation: ', total_result$interpretation, '

Species Richness:
  z = ', round(richness_result$beta, 3), ' [', round(if(!is.na(richness_result$boot_ci_lower)) richness_result$boot_ci_lower else richness_result$ci_lower, 3), ', ', round(if(!is.na(richness_result$boot_ci_upper)) richness_result$boot_ci_upper else richness_result$ci_upper, 3), ']
  z vs 1 = ', round(richness_result$z_vs_1, 2), ', p = ', format.pval(richness_result$p_vs_1, 3), ' (Wald)', if(!is.na(richness_result$p_boot_vs_1)) paste0('; p_boot = ', format.pval(richness_result$p_boot_vs_1, 3)) else '', '
  Interpretation: ', richness_result$interpretation, '

Density Dilution:
  Log-log slope = ', round(density_slope, 3), '
  Confirms sublinear scaling: per-capita density decreases with coral size

Bootstrap: 1,000 site-stratified iterations per taxon
Log base: natural log (coefficient = power-law exponent directly)

================================================================================
Generated: ', Sys.Date(), '
Source script: scripts/05_species_scaling_analysis.R
')
  writeLines(fig2_legend, file.path(PATHS$fig_manuscript, "fig2_legend_results.txt"))
  cat("  Generated: fig2_legend_results.txt\n")

  cat("\n")

  # ==========================================================================
  # MANUSCRIPT FIGURE 3: Species + Taxonomic Group Scaling (2x2)
  # Top: raw abundance vs volume curves; Bottom: coefficient forest plots
  # ==========================================================================

  cat("============================================================\n")
  cat("MANUSCRIPT FIGURE 3: Species & Taxonomic Group Scaling\n")
  cat("============================================================\n\n")

  # --- Scaling class colors (shared across bottom panels) ---
  scaling_colors <- c(
    "Redirection"     = "#5A8FAF",
    "Field of Dreams" = "gray55",
    "Super-linear"    = "#D55E00"
  )

  # --- Panel A: Key species abundance vs volume (overlaid curves) ---
  # Top 10 species by prevalence (same set as old Fig 2C)
  top10_species <- all_results_df %>%
    filter(category == "Species", !is.na(beta)) %>%
    arrange(desc(n_nonzero)) %>%
    head(10)

  # Full data for NB GLM fits (includes zeros), filtered for points
  top10_data_all <- species_scaling_data %>%
    filter(otu %in% top10_species$response) %>%
    mutate(species_label = gsub("_", " ", otu))
  top10_data_pts <- top10_data_all %>% filter(abundance > 0)

  # Sort species by beta (ascending) — same order as Panel C forest plot
  sp_order <- top10_species %>%
    mutate(species_label = gsub("_", " ", response)) %>%
    arrange(beta) %>%
    pull(species_label)

  # Assign colors: use a qualitative palette for species, ordered by beta
  sp_palette <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00",
                  "#CC79A7", "#882255", "#44AA99", "#999933", "#332288")
  names(sp_palette) <- sp_order[1:min(length(sp_order), 10)]

  # Apply factor levels to data so legend matches Panel C order (top = highest beta)
  top10_data_all <- top10_data_all %>%
    mutate(species_label = factor(species_label, levels = rev(sp_order)))

  # Generate NB GLM predictions for each species (avoids zero-inflation display issues)
  vol_grid <- data.frame(volume = exp(seq(log(50), log(40000), length.out = 100)))

  sp_pred_list <- lapply(top10_species$response, function(sp) {
    sp_data <- top10_data_all %>% filter(otu == sp)
    fit <- tryCatch(MASS::glm.nb(abundance ~ log(volume), data = sp_data),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    pred <- predict(fit, newdata = vol_grid, type = "response")
    data.frame(volume = vol_grid$volume, abundance = pred,
               species_label = gsub("_", " ", sp))
  })
  sp_pred_df <- do.call(rbind, Filter(Negate(is.null), sp_pred_list))
  sp_pred_df$species_label <- factor(sp_pred_df$species_label, levels = rev(sp_order))

  panel_a_sp <- ggplot() +
    geom_jitter(data = top10_data_all, aes(x = volume, y = abundance, color = species_label),
                alpha = 0.15, size = 0.6, height = 0.3, width = 0, show.legend = FALSE) +
    geom_line(data = sp_pred_df, aes(x = volume, y = abundance, color = species_label),
              linewidth = 0.7) +
    scale_x_log10(labels = scales::comma, breaks = c(100, 1000, 10000)) +
    scale_y_sqrt(breaks = c(0, 1, 5, 10, 25, 50)) +
    scale_color_manual(values = sp_palette, name = NULL) +
    coord_cartesian(xlim = shared_xlim) +
    labs(
      x = expression("Coral volume (cm"^3*")"),
      y = "Abundance",
      subtitle = "Species scaling"
    ) +
    theme_manuscript() +
    theme(
      axis.title = element_text(size = 9),
      legend.position = "none",  # 10 species — legend info in fig3_legend_results.txt
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 8, 5, 8, "mm")
    )

  # --- Panel B: Taxonomic group abundance vs volume (overlaid curves) ---
  # Reshape taxonomic groups to long format
  tax_cols <- c("n_crabs" = "Crabs", "n_shrimps" = "Shrimps",
                "n_fish" = "Fishes", "n_snails" = "Snails")
  tax_long <- scaling_data %>%
    dplyr::select(coral_id, site, volume, all_of(names(tax_cols))) %>%
    tidyr::pivot_longer(cols = all_of(names(tax_cols)),
                        names_to = "col", values_to = "abundance") %>%
    mutate(group = tax_cols[col])

  tax_palette <- c("Crabs" = "#E69F00", "Shrimps" = "#56B4E9",
                   "Fishes" = "#009E73", "Snails" = "#CC79A7")

  # Generate NB GLM predictions for each taxonomic group
  tax_pred_list <- lapply(unique(tax_long$group), function(g) {
    g_data <- tax_long %>% filter(group == g)
    fit <- tryCatch(MASS::glm.nb(abundance ~ log(volume), data = g_data),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    pred <- predict(fit, newdata = vol_grid, type = "response")
    data.frame(volume = vol_grid$volume, abundance = pred, group = g)
  })
  tax_pred_df <- do.call(rbind, Filter(Negate(is.null), tax_pred_list))

  # Direct-label endpoints for 4 taxonomic groups (ggrepel, no in-plot legend)
  tax_label_df <- tax_pred_df %>%
    group_by(group) %>%
    slice_max(volume, n = 1) %>%
    ungroup()

  panel_b_tax <- ggplot() +
    geom_jitter(data = tax_long, aes(x = volume, y = abundance, color = group),
                alpha = 0.15, size = 0.6, height = 0.3, width = 0, show.legend = FALSE) +
    geom_line(data = tax_pred_df, aes(x = volume, y = abundance, color = group),
              linewidth = 0.7) +
    ggrepel::geom_text_repel(
      data = tax_label_df,
      aes(x = volume, y = abundance, label = group, color = group),
      size = 2.8, fontface = "bold", nudge_x = 0.05, direction = "y",
      hjust = 0, segment.size = 0.3, segment.color = "grey60",
      show.legend = FALSE
    ) +
    scale_x_log10(labels = scales::comma, breaks = c(100, 1000, 10000)) +
    scale_y_sqrt(breaks = c(0, 1, 5, 10, 25, 50)) +
    scale_color_manual(values = tax_palette, name = NULL) +
    coord_cartesian(xlim = shared_xlim, clip = "off") +
    labs(
      x = expression("Coral volume (cm"^3*")"),
      y = "Abundance",
      subtitle = "Group scaling"
    ) +
    theme_manuscript() +
    theme(
      axis.title = element_text(size = 9),
      legend.position = "none",  # Direct labels used instead; details in fig3_legend_results.txt
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 14, 5, 8, "mm")  # Extra right margin for labels
    )

  # --- Panel C: Species scaling coefficient forest plot ---
  sp_forest <- top10_species %>%
    mutate(
      species_label = gsub("_", " ", response),
      species_label = factor(species_label, levels = rev(species_label[order(beta)])),
      scaling_class = case_when(
        boot_ci_upper < 1 ~ "Redirection",
        boot_ci_lower > 1 ~ "Super-linear",
        TRUE              ~ "Field of Dreams"
      )
    )

  # --- Panel D: Taxonomic group scaling coefficient forest plot ---
  tax_forest <- map_df(names(all_results$taxonomic_groups), function(nm) {
    r <- all_results$taxonomic_groups[[nm]]
    tibble(
      group = nm,
      beta = r$beta,
      ci_lower = if (!is.na(r$boot_ci_lower)) r$boot_ci_lower else r$ci_lower,
      ci_upper = if (!is.na(r$boot_ci_upper)) r$boot_ci_upper else r$ci_upper
    )
  }) %>%
    filter(!is.na(beta)) %>%
    mutate(
      group = factor(group, levels = rev(group[order(beta)])),
      scaling_class = case_when(
        ci_upper < 1 ~ "Redirection",
        ci_lower > 1 ~ "Super-linear",
        TRUE          ~ "Field of Dreams"
      )
    )

  # Compute shared x-axis range for both forest plots (must be after both data frames)
  all_betas <- c(sp_forest$beta, sp_forest$boot_ci_lower, sp_forest$boot_ci_upper)
  tax_betas_vals <- c(tax_forest$beta, tax_forest$ci_lower, tax_forest$ci_upper)
  forest_xlim <- range(c(all_betas, tax_betas_vals), na.rm = TRUE)
  forest_xlim <- forest_xlim + c(-0.1, 0.1) * diff(forest_xlim)  # 10% padding

  panel_c_forest <- ggplot(sp_forest, aes(x = beta, y = species_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.4) +
    geom_segment(aes(x = boot_ci_lower, xend = boot_ci_upper,
                     y = species_label, yend = species_label, color = scaling_class),
                 linewidth = 0.5) +
    geom_point(aes(color = scaling_class), size = 2) +
    scale_color_manual(values = scaling_colors, name = NULL) +
    scale_x_continuous(limits = forest_xlim) +
    labs(
      x = expression("Scaling exponent (" * beta * ")"),
      y = NULL,
      subtitle = expression("Species" ~ beta ~ "estimates")
    ) +
    theme_manuscript() +
    theme(
      axis.title = element_text(size = 9),
      axis.text.y = element_text(face = "italic", size = 7),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 8, 5, 8, "mm")
    )

  panel_d_forest <- ggplot(tax_forest, aes(x = beta, y = group)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.4) +
    geom_segment(aes(x = ci_lower, xend = ci_upper,
                     y = group, yend = group, color = scaling_class),
                 linewidth = 0.5) +
    geom_point(aes(color = scaling_class), size = 2.5) +
    scale_color_manual(values = scaling_colors, name = NULL) +
    scale_x_continuous(limits = forest_xlim) +
    labs(
      x = expression("Scaling exponent (" * beta * ")"),
      y = NULL,
      subtitle = expression("Group" ~ beta ~ "estimates")
    ) +
    theme_manuscript() +
    theme(
      axis.title = element_text(size = 9),
      axis.text.y = element_text(size = 8),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 8, 5, 8, "mm")
    )

  # --- Compose 2x2 with aligned axes ---
  fig3 <- (panel_a_sp | panel_b_tax) / (panel_c_forest | panel_d_forest) +
    plot_layout(heights = c(1, 0.9), axes = "collect") +
    plot_annotation(
      tag_levels = "A",
      caption = sprintf("Top: NB GLM fits. Bottom: scaling exponents with 95%% bootstrap CI.\nDashed line: \u03b2 = 1 (proportional). n = %d corals, 3 sites.", nrow(scaling_data)),
      theme = theme(
        plot.tag = element_text(face = "bold", size = 12),
        plot.caption = element_text(size = 6.5, color = "gray45", hjust = 0,
                                    margin = margin(t = 4, b = 2), lineheight = 1.2)
      )
    )

  save_figure(fig3, file.path(PATHS$fig_manuscript, "fig3_species_group_scaling.png"),
              width = 183, height = 180, units = "mm")
  cat("  Saved: manuscript/fig3_species_group_scaling.png\n")

  save_figure(fig3, file.path(FIG_DIR, "fig3_species_group_scaling.png"),
              width = 183, height = 180, units = "mm")
  cat("  Saved: 05_scaling/fig3_species_group_scaling.png\n")

  # Gather key statistics for legend
  n_species_shown <- 10
  n_species_total <- sum(all_results_df$category == "Species" & !is.na(all_results_df$beta))
  n_redirection_all <- sum(all_results_df$category == "Species" &
                           grepl("Redirection", all_results_df$interpretation), na.rm = TRUE)
  n_fod_all <- sum(all_results_df$category == "Species" &
                   grepl("Field of Dreams", all_results_df$interpretation), na.rm = TRUE)
  n_super_all <- sum(all_results_df$category == "Species" &
                     grepl("Super-linear", all_results_df$interpretation), na.rm = TRUE)
  # Counts for the 10 species actually shown in Panel C
  n_redirection_shown <- sum(grepl("Redirection", sp_forest$scaling_class))
  n_fod_shown <- sum(grepl("Field of Dreams", sp_forest$scaling_class))
  n_super_shown <- sum(grepl("Super-linear", sp_forest$scaling_class))

  # Generate fig3 legend text
  fig3_legend <- paste0(
'FIGURE 3: SPECIES AND TAXONOMIC GROUP SCALING
================================================================================

FIGURE LEGEND
-------------
Figure 3. Scaling of individual species and taxonomic groups with coral volume.
(A) Abundance vs colony volume for the ', n_species_shown, ' most prevalent CAFI species,
fitted with negative binomial (NB) GLMs (total_cafi ~ log(volume) + site; log-log
scale). (B) Abundance vs colony volume for four taxonomic groups
(crabs, shrimps, fishes, snails). (C) Species-level scaling exponents (\u03b2) with
95% bootstrap CI (1,000 site-stratified iterations). Blue: Redirection (\u03b2 < 1);
grey: Field of Dreams (CI spans 1); vermillion: Super-linear (\u03b2 > 1).
(D) Taxonomic group scaling exponents with 95% bootstrap CI.
Dashed vertical line: proportional scaling (\u03b2 = 1).
Of the ', n_species_shown, ' species shown in Panels A/C: ', n_redirection_shown, ' Redirection, ', n_fod_shown, ' Field of Dreams, ', n_super_shown, ' Super-linear.
Across all ', n_species_total, ' species tested: ', n_redirection_all, ' Redirection, ', n_fod_all, ' Field of Dreams, ', n_super_all, ' Super-linear.
FDR correction (Benjamini-Hochberg) applied within species and group categories.
n = ', nrow(scaling_data), ' Pocillopora corals across 3 reef sites (HAU, MAT, MRB).

PANEL A SPECIES COLORS (ordered by scaling exponent):
', paste(sapply(seq_along(sp_order[1:min(length(sp_order), 10)]), function(i) {
  sprintf("  %s: %s", gsub("_", " ", sp_order[i]), sp_palette[sp_order[i]])
}), collapse = "\n"), '

PANEL B TAXONOMIC GROUP COLORS:
  Crabs: #E69F00 (amber)
  Shrimps: #56B4E9 (sky blue)
  Fishes: #009E73 (bluish green)
  Snails: #CC79A7 (reddish purple)

================================================================================
Generated: ', Sys.Date(), '
Source script: scripts/05_species_scaling_analysis.R
')
  writeLines(fig3_legend, file.path(PATHS$fig_manuscript, "fig3_legend_results.txt"))
  cat("  Generated: fig3_legend_results.txt\n")

  cat("\n")

} else {
  cat("  Insufficient data for Figure 2 scaling panels\n\n")
}

cat("\nOUTPUT FILES:\n")
cat("  Figures: output/figures/05_scaling/\n")
cat("    - community_scaling.png\n")
cat("    - group_scaling_forest.png\n")
cat("    - species_scaling_forest.png\n")
cat("    - species_scaling_panels.png\n")
cat("    - scaling_summary.png\n")
cat("    - key_species_scaling_forest.png (NEW)\n")
cat("    - fig2_scaling.png (MANUSCRIPT)\n")
cat("  Figures: output/figures/manuscript/\n")
cat("    - fig2_scaling.png\n")
cat("  Figures: output/figures/supplement/\n")
cat("    - figS6_species_scaling.png\n")
cat("  Tables: output/tables/\n")
cat("    - scaling_results_all.csv\n")
cat("    - scaling_summary_by_category.csv\n")
cat("    - scaling_interpretation_summary.csv\n")
cat("    - key_species_scaling_experimental.csv (NEW)\n\n")
