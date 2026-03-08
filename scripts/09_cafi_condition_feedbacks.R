#!/usr/bin/env Rscript
# ============================================================================
# 09_cafi_condition_feedbacks.R - Bidirectional CAFI-Coral Condition Feedbacks
# ============================================================================
#
# PURPOSE: Test bidirectional relationships between CAFI communities and
#          coral physiological condition
#
# CORE QUESTIONS:
#   A. Do CAFI communities affect coral condition? (CAFI -> Condition)
#   B. Does coral condition attract more/different CAFI? (Condition -> CAFI)
#
# FUNCTIONAL GROUP HYPOTHESES:
#   - Trapezia (defenders): Expect POSITIVE effect on condition
#     (predator defense, cleaning behavior, waste nutrient cycling)
#   - Resident Fish (nutrient providers): Expect POSITIVE effect
#     (ammonium excretion enhances zooxanthellae)
#   - Galeropsis (tissue consumer): Expect NEGATIVE effect
#     (Coralliophilinae — feeds on coral tissue and biofilm)
#
# DATA:
#   - n=84 corals have position-corrected condition scores (PC1 of physiology)
#   - Condition = position-corrected PC1 of protein, carbs, zoox, AFDW
#   - Position correction removes sampling artifact from stump length
#
# ANALYSES:
#   A. CAFI -> Condition Effects (fixed-effect site models)
#   B. Condition -> CAFI Effects (reverse direction test)
#   C. Functional group-specific models
#   D. Condition x Volume interaction
#   E. Path analysis (if sample size supports)
#   F. Key species from experimental paper
#   G. Neighborhood effects (analog to experimental treatment)
#   H. Landscape-only effects on condition (no CAFI predictors)
#
# OUTPUTS:
#
#   Figures (14):
#   - output/figures/manuscript/fig5_feedbacks.png              (2-panel BEF main text figure)
#   - output/figures/supplement/figS10_rarefaction_sensitivity.png (rarefaction depth)
#   - output/figures/supplement/figS12_bef_variance_partitioning.png (BEF variance partitioning + path)
#   - output/figures/supplement/figS13_condition_details.png    (a priori forest + rarefied + exploratory + scatters + bidirectional)
#   - output/figures/supplement/figS16_species_trait_heatmap.png (top 20 species × 5 traits)
#   - output/figures/supplement/figS17_species_trait_biplots.png (strongest associations)
#   - output/figures/feedbacks/cafi_condition_effects.png       (all CAFI predictors forest)
#   - output/figures/feedbacks/functional_effects_forest.png    (functional group forest)
#   - output/figures/feedbacks/key_species_effects.png          (key species forest plot)
#   - output/figures/feedbacks/functional_vs_key_species.png    (combined comparison)
#   - output/figures/feedbacks/neighborhood_effects.png         (Part G neighborhood panels)
#   - output/figures/feedbacks/landscape_condition_effects.png  (Part H landscape)
#   - output/figures/feedbacks/diagnostics_richness_model.png   (richness model diagnostics)
#
#   Tables (14):
#   - output/tables/cafi_condition_models.csv               (CAFI -> condition results)
#   - output/tables/reverse_direction_models.csv            (condition -> CAFI results)
#   - output/tables/functional_effects.csv                  (functional group effects)
#   - output/tables/richness_abundance_artifact.csv         (richness artifact test)
#   - output/tables/key_species_effects.csv                 (key species condition effects)
#   - output/tables/species_trait_correlations.csv           (Part F2: species x trait r)
#   - output/tables/individual_physiology_cafi_responses.csv (Part F3: individual traits)
#   - output/tables/cross_study_species_comparison.csv      (Part F4: survey vs experiment)
#   - output/tables/neighborhood_effects.csv                (Part G: neighborhood models)
#   - output/tables/landscape_condition_effects.csv         (Part H: landscape effects)
#   - output/tables/landscape_model_comparison.csv          (Part H: model comparison)
#   - output/tables/site_condition_means.csv                (Part H: site means)
#   - output/tables/community_transform_sensitivity.csv     (Part A3: transform sensitivity)
#   - output/tables/09_cafi_condition_feedbacks_stats_summary.csv (stats summary)
#
#   Objects (2):
#   - output/objects/richness_comparison_results.rds        (richness artifact data)
#   - output/objects/analysis_data_rarefied.rds             (rarefied analysis dataset)
#
#   Text (1):
#   - output/figures/manuscript/fig5_legend_results.txt     (figure legend + results)
#
# DEPENDENCIES: 00_setup.R, 01_load_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-21
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    BIDIRECTIONAL CAFI-CORAL CONDITION FEEDBACKS\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP AND DATA LOADING
# ============================================================================

# Load setup (packages, paths, theme)
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))

# Load additional packages for model inference and mediation analysis
if (!requireNamespace("lmerTest", quietly = TRUE)) {
  cat("Note: lmerTest package not available - using lme4 approximation for p-values\n")
  LMERTEST_AVAILABLE <- FALSE
} else {
  suppressPackageStartupMessages(library(lmerTest))
  LMERTEST_AVAILABLE <- TRUE
}

# Load additional packages for path analysis
if (requireNamespace("lavaan", quietly = TRUE)) {
  suppressPackageStartupMessages(library(lavaan))
  LAVAAN_AVAILABLE <- TRUE
} else {
  LAVAAN_AVAILABLE <- FALSE
  cat("Note: lavaan package not available - path analysis will be skipped\n")
}

# Load processed data objects (load_object adds .rds automatically)
coral_master <- load_object("coral_master")
cafi_clean <- load_object("cafi_clean")
condition_scores <- load_object("condition_scores")

# Load PCA-based CAFI scores (aligned with Stier et al. manuscript methodology)
cafi_pca_results <- tryCatch(
  load_object("cafi_pca_results"),
  error = function(e) {
    cat("Note: cafi_pca_results not found - will compute PC1_CAFI from scratch\n")
    NULL
  }
)

# Create figure subdirectories
fig_dir <- file.path(PATHS$figures, "feedbacks")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("Data loaded successfully\n\n")

# Define helper functions for results tracking
init_results_df <- function() {
  tibble(
    hypothesis = character(),
    question = character(),
    test_name = character(),
    test_statistic = numeric(),
    df = character(),  # Character to handle various formats
    p_value = numeric(),
    effect_size = numeric(),
    interpretation = character(),
    effect_type = character(),
    n = numeric(),
    notes = character()
  )
}

create_result_row <- function(hypothesis, question, test_name, test_statistic = NA,
                              df = NA, p_value = NA, effect_size = NA, interpretation = NA,
                              effect_type = NA, n = NA, notes = NA) {
  tibble(
    hypothesis = hypothesis,
    question = question,
    test_name = test_name,
    test_statistic = as.numeric(test_statistic),
    df = as.character(df),  # Always character to handle "NA", "1", "2.5" etc.
    p_value = as.numeric(p_value),
    effect_size = as.numeric(effect_size),
    interpretation = as.character(interpretation),
    effect_type = as.character(effect_type),
    n = as.numeric(n),
    notes = as.character(notes)
  )
}

save_stats_summary <- function(results_df, script_name, title) {
  # Save to CSV
  output_path <- file.path(PATHS$tables, paste0(script_name, "_stats_summary.csv"))
  write_csv(results_df, output_path)
  cat("Saved statistics summary:", output_path, "\n")
}

# ============================================================================
# 1. PREPARE ANALYSIS DATA
# ============================================================================

cat("1. Preparing analysis data...\n")

# Calculate functional group abundances per coral
functional_predictors <- cafi_clean %>%
  group_by(coral_id) %>%
  summarise(
    total_cafi = n(),

    # Functional group counts (based on ecological role)
    n_trapezia = sum(functional_group == "Trapezia", na.rm = TRUE),
    n_resident_fish = sum(functional_group == "Resident Fish", na.rm = TRUE),
    n_corallivore = sum(functional_group == "Gastropod", na.rm = TRUE),
    n_galeropsis = sum(str_detect(tolower(otu), "galeropsis"), na.rm = TRUE),
    n_other_crab = sum(functional_group == "Other Crab", na.rm = TRUE),
    n_shrimp = sum(functional_group == "Shrimp", na.rm = TRUE),

    # Diversity metrics
    otu_richness = n_distinct(otu),
    shannon = vegan::diversity(table(otu), index = "shannon"),

    # Taxonomic group counts
    n_crabs = sum(type == "crab", na.rm = TRUE),
    n_fish = sum(type == "fish", na.rm = TRUE),
    n_snails = sum(type == "snail", na.rm = TRUE),
    n_shrimps = sum(type == "shrimp", na.rm = TRUE),

    .groups = "drop"
  )

cat("   Calculated CAFI metrics for", nrow(functional_predictors), "corals\n")

# Merge all datasets including PC1_CAFI (following Stier et al. methodology)
analysis_data <- condition_scores %>%
  dplyr::select(coral_id, site, condition_score, any_of(c("protein_corr", "carb_corr",
                                                     "zoox_corr", "afdw_corr"))) %>%
  left_join(functional_predictors, by = "coral_id") %>%
  left_join(coral_master %>% dplyr::select(coral_id, volume, morphotype, depth_m, any_of(c("pc1_cafi", "pc2_cafi"))),
            by = "coral_id") %>%
  filter(!is.na(condition_score), !is.na(total_cafi), !is.na(volume)) %>%
  mutate(
    log_volume = log(volume),
    site = factor(site)
  )

n_complete <- nrow(analysis_data)

cat("   Complete cases:", n_complete, "corals with condition + CAFI data\n")
cat("   Sites:", paste(unique(analysis_data$site), collapse = ", "), "\n")
cat("   Volume range:", round(min(analysis_data$volume)), "-",
    round(max(analysis_data$volume)), "cm^3\n")
cat("   Condition range:", round(min(analysis_data$condition_score), 2), "to",
    round(max(analysis_data$condition_score), 2), "\n\n")

# Summary statistics
cat("Functional Group Distributions:\n")
cat("   Trapezia:", sum(analysis_data$n_trapezia),
    "(range:", min(analysis_data$n_trapezia), "-", max(analysis_data$n_trapezia), ")\n")
cat("   Resident Fish:", sum(analysis_data$n_resident_fish),
    "(range:", min(analysis_data$n_resident_fish), "-", max(analysis_data$n_resident_fish), ")\n")
cat("   Galeropsis:", sum(analysis_data$n_galeropsis),
    "(range:", min(analysis_data$n_galeropsis), "-", max(analysis_data$n_galeropsis), ")\n\n")

# ============================================================================
# SENSITIVITY: POSITION-CORRECTED vs RAW CONDITION SCORES
# ============================================================================
# (Position correction regresses out stump_length + nubbin_length before PCA)
# Brief check: does the key result (richness → condition) change with raw scores?

cat("Position Correction Sensitivity Check:\n")
cat("--------------------------------------\n")

# Attempt to load raw (uncorrected) physiology data and compute a raw PC1
raw_phys <- tryCatch(load_object("condition_scores"), error = function(e) NULL)

if (!is.null(raw_phys) && all(c("protein_corr", "carb_corr", "zoox_corr", "afdw_corr") %in% names(raw_phys))) {
  # The corrected scores are what we already use; build a raw version from coral_master
  # Check if raw (uncorrected) physiology columns exist in coral_master
  raw_physio_cols <- c("protein_norm", "carb_norm", "zoox_norm", "afdw_norm")
  alt_raw_cols <- c("protein_mg_cm2", "carb_mg_cm2", "zoox_cells_cm2", "afdw_mg_cm2")

  raw_cols_available <- raw_physio_cols[raw_physio_cols %in% names(coral_master)]
  if (length(raw_cols_available) == 0) {
    raw_cols_available <- alt_raw_cols[alt_raw_cols %in% names(coral_master)]
  }

  if (length(raw_cols_available) >= 3) {
    raw_pca_data <- coral_master %>%
      dplyr::select(coral_id, all_of(raw_cols_available)) %>%
      filter(complete.cases(.))

    if (nrow(raw_pca_data) >= 20) {
      raw_pca <- prcomp(raw_pca_data %>% dplyr::select(-coral_id), center = TRUE, scale. = TRUE)
      raw_pc1 <- raw_pca$x[, 1]
      # Align sign with corrected condition_score
      raw_pc1_df <- tibble(coral_id = raw_pca_data$coral_id, raw_condition = raw_pc1)

      sensitivity_data <- analysis_data %>%
        left_join(raw_pc1_df, by = "coral_id") %>%
        filter(!is.na(raw_condition))

      if (nrow(sensitivity_data) >= 20) {
        # Align sign: ensure raw_condition is positively correlated with corrected
        if (cor(sensitivity_data$condition_score, sensitivity_data$raw_condition, use = "complete") < 0) {
          sensitivity_data$raw_condition <- -sensitivity_data$raw_condition
        }

        # Test richness → condition with CORRECTED scores
        m_corrected <- lm(condition_score ~ otu_richness + log_volume + site, data = sensitivity_data)
        p_corrected <- summary(m_corrected)$coefficients["otu_richness", "Pr(>|t|)"]
        b_corrected <- coef(m_corrected)["otu_richness"]

        # Test richness → condition with RAW scores
        m_raw_sens <- lm(raw_condition ~ otu_richness + log_volume + site, data = sensitivity_data)
        p_raw_sens <- summary(m_raw_sens)$coefficients["otu_richness", "Pr(>|t|)"]
        b_raw_sens <- coef(m_raw_sens)["otu_richness"]

        cat(sprintf("  Corrected condition: richness beta = %.4f, p = %.4f (n=%d)\n",
                    b_corrected, p_corrected, nrow(sensitivity_data)))
        cat(sprintf("  Raw condition:       richness beta = %.4f, p = %.4f (n=%d)\n",
                    b_raw_sens, p_raw_sens, nrow(sensitivity_data)))

        both_ns <- p_corrected > 0.05 & p_raw_sens > 0.05
        both_sig <- p_corrected < 0.05 & p_raw_sens < 0.05
        cat(sprintf("  Qualitatively similar: %s\n",
                    ifelse(both_ns | both_sig, "YES", "NO — position correction changes significance")))
      } else {
        cat("  Insufficient overlap between raw and corrected data (n < 20)\n")
      }
    } else {
      cat("  Insufficient raw physiology data for PCA (n < 20)\n")
    }
  } else {
    cat("  Raw (uncorrected) physiology columns not found in coral_master.\n")
    cat("  Position correction sensitivity cannot be tested.\n")
    cat("  (This is expected if raw columns were consumed during 01_load_data.R)\n")
  }
} else {
  cat("  Condition scores object missing required corrected columns.\n")
}
cat("\n")

# ============================================================================
# POWER ANALYSIS (Q3: CAFI-Condition Feedbacks)
# ============================================================================
# With n=84 corals with physiology data:
# - Power to detect R² = 0.10 (medium effect): ~80% at α=0.05
# - Power to detect R² = 0.05 (small effect): ~45% at α=0.05
# - Adequate power for medium effects; null results for medium+ effects are credible
if (requireNamespace("pwr", quietly = TRUE)) {
  n_physio <- sum(!is.na(coral_master$condition_score))
  power_med <- pwr::pwr.f2.test(u = 1, v = n_physio - 5, f2 = 0.10/0.90, sig.level = 0.05)
  power_sm <- pwr::pwr.f2.test(u = 1, v = n_physio - 5, f2 = 0.05/0.95, sig.level = 0.05)
  cat(sprintf("Power analysis (n=%d physio corals):\n  Medium effect (R²=0.10): %.0f%%\n  Small effect (R²=0.05): %.0f%%\n\n",
              n_physio, power_med$power * 100, power_sm$power * 100))

  # Cross-study power comparison: can the survey detect the experiment's effect sizes?
  # Experiment (n=54) detected PC1_CAFI → PC1_Coral with R² ≈ 0.12-0.15
  # What is our minimum detectable effect at 80% power?
  mde <- pwr::pwr.f2.test(u = 1, v = n_physio - 5, power = 0.80, sig.level = 0.05)
  mde_r2 <- mde$f2 / (1 + mde$f2)

  # Can we detect the experiment's effect (R² ≈ 0.12)?
  expt_r2 <- 0.12
  power_expt <- pwr::pwr.f2.test(u = 1, v = n_physio - 5,
                                   f2 = expt_r2 / (1 - expt_r2), sig.level = 0.05)

  cat("Cross-study power comparison:\n")
  cat(sprintf("  Minimum detectable R² at 80%% power: %.3f\n", mde_r2))
  cat(sprintf("  Experiment's effect (R²≈%.2f): survey has %.0f%% power to detect it\n",
              expt_r2, power_expt$power * 100))
  if (power_expt$power >= 0.80) {
    cat("  → Survey IS adequately powered to detect the experiment's effect size\n")
    cat("  → Null results are credible: feedbacks may genuinely differ in established communities\n\n")
  } else {
    cat("  → Survey is UNDERPOWERED to detect the experiment's effect size\n")
    cat("  → Null results should be interpreted cautiously\n\n")
  }
}

# ============================================================================
# 2. INITIALIZE RESULTS TRACKING
# ============================================================================

stats_results <- init_results_df()

# ============================================================================
# PART A: CAFI -> CONDITION EFFECTS
# ============================================================================

cat("============================================================\n")
cat("PART A: CAFI -> CORAL CONDITION EFFECTS\n")
cat("============================================================\n\n")

# Function to run and summarize linear model with fixed-effect site
# LM is appropriate here because the response (condition_score) is a continuous
# PCA score (PC1 of physiology), not a count or proportion.
# PRIMARY INFERENCE: OLS standard errors (Breusch-Pagan test confirmed no
# heteroscedasticity for richness model, BP p = 0.83; residual SD ratio = 1.56).
# HC3 robust SEs reported as supplement sensitivity check (known to be overly
# conservative at n < 100; Long & Ervin 2000).
# Note: 3 sites is insufficient for random intercepts (Bolker et al. 2009 recommends >=5-6)
run_condition_model <- function(data, predictor_name, predictor_col, transform = "none") {
  # Apply sqrt transform to count predictors (mirrors experimental companion paper)
  if (transform == "sqrt") {
    data[[predictor_col]] <- sqrt(data[[predictor_col]])
  }

  # Build formula with fixed-effect site (3 levels insufficient for RE)
  formula_str <- paste("condition_score ~", predictor_col, "+ log_volume + site")

  tryCatch({
    # Fit linear model with site as fixed effect
    model <- lm(as.formula(formula_str), data = data)

    # Extract coefficient info (OLS — primary inference)
    model_summary <- summary(model)
    coef_table <- coef(model_summary)

    # Get the predictor coefficient
    pred_row <- which(rownames(coef_table) == predictor_col)
    if (length(pred_row) == 0) {
      warning("Predictor '", predictor_col, "' not found in model coefficients: ",
              paste(rownames(coef_table), collapse = ", "))
      return(NULL)
    }

    estimate <- coef_table[pred_row, "Estimate"]
    se <- coef_table[pred_row, "Std. Error"]
    t_val <- coef_table[pred_row, "t value"]
    df_val <- model$df.residual
    p_val <- coef_table[pred_row, "Pr(>|t|)"]

    # HC3 robust SE (supplement sensitivity only)
    se_robust <- se  # default
    p_val_robust <- p_val
    bp_p <- NA_real_
    if (requireNamespace("sandwich", quietly = TRUE) &&
        requireNamespace("lmtest", quietly = TRUE)) {
      vcov_hc3 <- sandwich::vcovHC(model, type = "HC3")
      robust_coefs <- lmtest::coeftest(model, vcov. = vcov_hc3)
      se_robust <- robust_coefs[pred_row, "Std. Error"]
      p_val_robust <- robust_coefs[pred_row, "Pr(>|t|)"]
      # Breusch-Pagan test
      bp_test <- lmtest::bptest(model)
      bp_p <- bp_test$p.value
    }

    # 95% CI (using OLS SE — primary inference)
    ci_lower <- estimate - qt(0.975, df_val) * se
    ci_upper <- estimate + qt(0.975, df_val) * se

    # Partial standardized β: fit model on z-scored data
    data_std <- data
    data_std[[predictor_col]] <- scale(data_std[[predictor_col]])
    data_std$condition_score <- scale(data_std$condition_score)
    data_std$log_volume <- scale(data_std$log_volume)
    model_std <- tryCatch(
      lm(as.formula(formula_str), data = data_std),
      error = function(e) NULL
    )
    beta_std <- if (!is.null(model_std)) {
      coef(model_std)[predictor_col]
    } else {
      # Fallback: semi-standardized
      sd_predictor <- sd(data[[predictor_col]], na.rm = TRUE)
      sd_response <- sd(data$condition_score, na.rm = TRUE)
      ifelse(sd_predictor > 0 & sd_response > 0, estimate * sd_predictor / sd_response, NA)
    }

    # R-squared (adjusted)
    r2_adj <- model_summary$adj.r.squared

    result <- data.frame(
      direction = "CAFI -> Condition",
      predictor = predictor_name,
      estimate = estimate,
      beta_std = beta_std,
      se = se,
      se_robust = se_robust,
      t_value = t_val,
      df = df_val,
      p_value = p_val,
      p_value_robust = p_val_robust,
      bp_p = bp_p,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      r2_adj = r2_adj,
      n = nrow(data),
      significant = p_val < 0.05  # OLS p-value (primary)
    )

    return(list(model = model, result = result))

  }, error = function(e) {
    cat("   Warning: Model failed for", predictor_name, "-", e$message, "\n")
    return(NULL)
  })
}

# Test each CAFI metric as predictor of condition
# PC1_CAFI is the primary metric following Stier et al. (2024) methodology
cafi_predictors <- list(
  c("PC1_CAFI (community)", "pc1_cafi", "none"),
  c("Total CAFI", "total_cafi", "sqrt"),
  c("Trapezia abundance", "n_trapezia", "sqrt"),
  c("Resident Fish abundance", "n_resident_fish", "sqrt"),
  c("Galeropsis abundance", "n_galeropsis", "sqrt"),
  c("Species richness", "otu_richness", "none"),
  c("Shannon diversity", "shannon", "none")
)

cafi_to_condition_results <- list()
cafi_to_condition_models <- list()

cat("Testing CAFI predictors of coral condition...\n")
cat("Model: Condition ~ CAFI_metric + log(Volume) + Site (fixed effect)\n\n")

for (pred in cafi_predictors) {
  cat("   Testing:", pred[1], "...\n")
  result <- run_condition_model(analysis_data, pred[1], pred[2], transform = pred[3])

  if (!is.null(result)) {
    cafi_to_condition_results[[pred[1]]] <- result$result
    cafi_to_condition_models[[pred[1]]] <- result$model

    # Print result
    r <- result$result
    sig_star <- ifelse(r$p_value < 0.001, "***",
                       ifelse(r$p_value < 0.01, "**",
                              ifelse(r$p_value < 0.05, "*",
                                     ifelse(r$p_value < 0.1, ".", ""))))
    cat("      beta =", round(r$estimate, 4),
        ", SE =", round(r$se, 4),
        ", t =", round(r$t_value, 2),
        ", p =", format.pval(r$p_value, 3), sig_star, "\n")
  }
}

# Combine CAFI -> Condition results
cafi_to_condition_df <- bind_rows(cafi_to_condition_results)

# Apply three-tier correction for multiple testing:
#   A PRIORI BEF (k=2): Hochberg FWER — Species richness, Total CAFI
#     (predicted by BEF theory: complementarity + abundance pathways)
#   EXPLORATORY (k=4): BH-FDR — Shannon, Trapezia, Fish, Galeropsis
#   SUPPLEMENT COMPOSITION (k=1): PC1_CAFI tested separately
#     (community composition identity is a distinct question from BEF)
# See Part A4 for full BEF analysis and Part I for global FDR sensitivity check.

apriori_names <- c("Species richness", "Total CAFI")
exploratory_names <- c("Shannon diversity", "Trapezia abundance",
                       "Resident Fish abundance", "Galeropsis abundance")
supplement_composition_names <- c("PC1_CAFI (community)")

cafi_to_condition_df <- cafi_to_condition_df %>%
  mutate(
    hypothesis_type = case_when(
      predictor %in% apriori_names ~ "a_priori",
      predictor %in% supplement_composition_names ~ "supplement_composition",
      TRUE ~ "exploratory"
    ),
    significant = p_value < 0.05  # OLS (primary)
  )

# Hochberg for a priori BEF (k=2, using OLS p-values — primary inference)
apriori_idx <- which(cafi_to_condition_df$hypothesis_type == "a_priori")
cafi_to_condition_df$p_corrected <- NA_real_
cafi_to_condition_df$p_corrected[apriori_idx] <-
  p.adjust(cafi_to_condition_df$p_value[apriori_idx], method = "hochberg")

# BH-FDR for exploratory (k=4)
exploratory_idx <- which(cafi_to_condition_df$hypothesis_type == "exploratory")
cafi_to_condition_df$p_corrected[exploratory_idx] <-
  p.adjust(cafi_to_condition_df$p_value[exploratory_idx], method = "BH")

# PC1_CAFI: no correction needed (k=1, separate question → supplement)
supp_idx <- which(cafi_to_condition_df$hypothesis_type == "supplement_composition")
cafi_to_condition_df$p_corrected[supp_idx] <-
  cafi_to_condition_df$p_value[supp_idx]

cafi_to_condition_df$significant_corrected <- cafi_to_condition_df$p_corrected < 0.05

# HC3 sensitivity: same structure on robust p-values (supplement)
cafi_to_condition_df$p_corrected_hc3 <- NA_real_
cafi_to_condition_df$p_corrected_hc3[apriori_idx] <-
  p.adjust(cafi_to_condition_df$p_value_robust[apriori_idx], method = "hochberg")
cafi_to_condition_df$p_corrected_hc3[exploratory_idx] <-
  p.adjust(cafi_to_condition_df$p_value_robust[exploratory_idx], method = "BH")
cafi_to_condition_df$p_corrected_hc3[supp_idx] <-
  cafi_to_condition_df$p_value_robust[supp_idx]

# Also keep legacy BH-FDR across all 7 for backward compatibility / sensitivity
cafi_to_condition_df$p_fdr <- p.adjust(cafi_to_condition_df$p_value, method = "BH")
cafi_to_condition_df$significant_fdr <- cafi_to_condition_df$p_fdr < 0.05

cat("\nBreusch-Pagan heteroscedasticity tests (justification for OLS as primary):\n")
for (i in 1:nrow(cafi_to_condition_df)) {
  row <- cafi_to_condition_df[i, ]
  if (!is.na(row$bp_p)) {
    cat(sprintf("   %-30s BP p = %.4f %s\n", row$predictor, row$bp_p,
                ifelse(row$bp_p < 0.05, "(heteroscedastic)", "(homoscedastic)")))
  }
}

cat("\nThree-tier corrected significance (OLS primary; HC3 in supplement):\n")
cat("  A priori BEF (Hochberg FWER, k=2): Species richness, Total CAFI\n")
cat("  Exploratory (BH-FDR, k=4): Shannon, Trapezia, Fish, Galeropsis\n")
cat("  Supplement composition (uncorrected, k=1): PC1_CAFI\n")
cat("(See Part A4 for BEF rationale; Part I for global FDR sensitivity.)\n\n")
for (i in 1:nrow(cafi_to_condition_df)) {
  row <- cafi_to_condition_df[i, ]
  method <- ifelse(row$hypothesis_type == "a_priori", "Hochberg", "BH-FDR")
  sig_star <- ifelse(row$significant_corrected, " *SIG*", "")
  cat(sprintf("   %-30s p_OLS = %.4f, p_%s = %.4f%s  (HC3: p = %.4f, p_%s = %.4f)\n",
              row$predictor, row$p_value, method, row$p_corrected, sig_star,
              row$p_value_robust, method, row$p_corrected_hc3))
}
cat("\n")

# --- Model diagnostics for CAFI -> Condition models ---
cat("Model Diagnostics (CAFI -> Condition):\n")
cat("--------------------------------------\n")

# VIF for the richness model (strongest signal)
if ("Species richness" %in% names(cafi_to_condition_models)) {
  richness_model <- cafi_to_condition_models[["Species richness"]]
  if (requireNamespace("car", quietly = TRUE)) {
    vif_vals <- car::vif(richness_model)
    cat("   VIF (richness model):\n")
    for (v in seq_along(vif_vals)) {
      cat(sprintf("     %-20s VIF = %.2f %s\n",
                  names(vif_vals)[v], vif_vals[v],
                  ifelse(vif_vals[v] > 5, "WARNING: HIGH", "")))
    }
  }

  # Residual normality (Shapiro-Wilk)
  resids <- residuals(richness_model)
  sw_test <- shapiro.test(resids)
  cat(sprintf("   Shapiro-Wilk residual normality: W = %.4f, p = %.4f %s\n",
              sw_test$statistic, sw_test$p.value,
              ifelse(sw_test$p.value < 0.05, "(non-normal)", "(OK)")))

  # Save diagnostic plot
  png(file.path(fig_dir, "diagnostics_richness_model.png"),
      width = 10, height = 8, units = "in", res = 200, bg = "white")
  par(mfrow = c(2, 2))
  plot(richness_model)
  dev.off()
  cat("   Saved: diagnostics_richness_model.png\n")
}
cat("\n")

# ---- Sensitivity power analysis for CAFI -> Condition models ----
# Post-hoc power based on observed data is invalid (Hoenig & Heisey 2001).
# Instead, report the minimum detectable effect size at 80% power.
cat("Sensitivity power analysis (CAFI -> Condition):\n")
cat("(What effect could this study detect at 80% power?)\n")

n_condition <- nrow(analysis_data)
n_predictors_cond <- 4  # CAFI_metric + log_volume + 2 site dummies = 4 non-intercept (testing 1)

if (requireNamespace("pwr", quietly = TRUE)) {
  # Minimum detectable f2 at alpha=0.05, power=0.80
  pwr_result_cond <- pwr::pwr.f2.test(
    u = 1,  # numerator df (1 predictor being tested)
    v = n_condition - n_predictors_cond - 2,  # denominator df
    f2 = NULL,
    sig.level = 0.05,
    power = 0.80
  )
  min_f2_cond <- pwr_result_cond$f2
  min_partial_r2_cond <- min_f2_cond / (1 + min_f2_cond)

  cat("  Sample size for condition models:", n_condition, "\n")
  cat("  Minimum detectable partial R2 at 80% power:", round(min_partial_r2_cond, 3), "\n")
  cat("  Minimum detectable Cohen's f2:", round(min_f2_cond, 3), "\n")
  cat("  This corresponds to a",
      ifelse(min_f2_cond < 0.02, "very small", ifelse(min_f2_cond < 0.15, "small-to-medium", "medium")),
      "effect.\n")
  cat("  Effects smaller than partial R2 =", round(min_partial_r2_cond, 3),
      "cannot be reliably detected.\n\n")
} else {
  cat("  (Install 'pwr' package for formal power analysis)\n")
  # Rough approximation: for n~84, minimum detectable r ~ 0.30
  cat("  Approximate minimum detectable correlation: r ~ 0.30 (medium effect)\n\n")
}

# ============================================================================
# PART A2: RICHNESS-ABUNDANCE ARTIFACT TEST
# ============================================================================
# The species richness → condition signal may be driven by abundance, not diversity.
# Here we test this by comparing raw vs rarefied richness effects.

cat("============================================================\n")
cat("PART A2: RICHNESS-ABUNDANCE ARTIFACT TEST\n")
cat("============================================================\n\n")

cat("Testing whether richness → condition is an abundance artifact...\n")
cat("Correlation between richness and abundance: r =",
    round(cor(analysis_data$otu_richness, analysis_data$total_cafi), 3), "\n\n")

# --- Compute rarefied richness ---
cat("Computing rarefied richness...\n")
comm_matrix <- load_object("community_matrix")
abund <- rowSums(comm_matrix)

# Rarefy at n=20 (reasonable depth that keeps most corals)
rarefy_depth <- 20
has_enough <- abund >= rarefy_depth
n_kept <- sum(has_enough)
cat(sprintf("  Rarefying to n=%d individuals (keeps %d of %d corals)\n",
    rarefy_depth, n_kept, length(has_enough)))

comm_sub <- comm_matrix[has_enough, ]
set.seed(42)  # Reproducibility for rarefied richness
rare_rich <- rarefy(comm_sub, sample = rarefy_depth)

# Match to analysis data
analysis_data_rare <- analysis_data %>%
  filter(coral_id %in% names(rare_rich)) %>%
  mutate(rarefied_richness = rare_rich[coral_id])

cat(sprintf("  Sample size with condition + rarefied richness: n = %d\n\n",
    nrow(analysis_data_rare)))

# --- Test raw vs rarefied richness ---
cat("Comparing raw vs rarefied richness effects on condition:\n")
cat("------------------------------------------------------------\n\n")

# Model 1: Raw richness (on full sample)
m_raw_full <- lm(condition_score ~ otu_richness + log_volume + site, data = analysis_data)
s_raw_full <- summary(m_raw_full)
p_raw_full <- s_raw_full$coefficients["otu_richness", "Pr(>|t|)"]
b_raw_full <- coef(m_raw_full)["otu_richness"]

cat(sprintf("RAW RICHNESS (full sample, n=%d):\n", nrow(analysis_data)))
cat(sprintf("  β = %.4f, SE = %.4f, p = %.4f %s\n\n",
    b_raw_full,
    s_raw_full$coefficients["otu_richness", "Std. Error"],
    p_raw_full,
    ifelse(p_raw_full < 0.05, "*", "")))

# Model 2: Raw richness (on subset with n≥20)
m_raw_sub <- lm(condition_score ~ otu_richness + log_volume + site, data = analysis_data_rare)
s_raw_sub <- summary(m_raw_sub)
p_raw_sub <- s_raw_sub$coefficients["otu_richness", "Pr(>|t|)"]
b_raw_sub <- coef(m_raw_sub)["otu_richness"]

cat(sprintf("RAW RICHNESS (subsample n≥20, n=%d):\n", nrow(analysis_data_rare)))
cat(sprintf("  β = %.4f, SE = %.4f, p = %.4f %s\n\n",
    b_raw_sub,
    s_raw_sub$coefficients["otu_richness", "Std. Error"],
    p_raw_sub,
    ifelse(p_raw_sub < 0.05, "*", "")))

# Model 3: Rarefied richness
m_rare <- lm(condition_score ~ rarefied_richness + log_volume + site, data = analysis_data_rare)
s_rare <- summary(m_rare)
p_rare <- s_rare$coefficients["rarefied_richness", "Pr(>|t|)"]
b_rare <- coef(m_rare)["rarefied_richness"]

cat(sprintf("RAREFIED RICHNESS (n=%d, depth=%d):\n", nrow(analysis_data_rare), rarefy_depth))
cat(sprintf("  β = %.4f, SE = %.4f, p = %.4f %s\n\n",
    b_rare,
    s_rare$coefficients["rarefied_richness", "Std. Error"],
    p_rare,
    ifelse(p_rare < 0.05, "*", "")))

# Model 4: Abundance on subset
m_abund_sub <- lm(condition_score ~ total_cafi + log_volume + site, data = analysis_data_rare)
s_abund_sub <- summary(m_abund_sub)
p_abund_sub <- s_abund_sub$coefficients["total_cafi", "Pr(>|t|)"]

cat(sprintf("ABUNDANCE (subsample, n=%d):\n", nrow(analysis_data_rare)))
cat(sprintf("  β = %.4f, SE = %.4f, p = %.4f %s\n\n",
    coef(m_abund_sub)["total_cafi"],
    s_abund_sub$coefficients["total_cafi", "Std. Error"],
    p_abund_sub,
    ifelse(p_abund_sub < 0.05, "*", "")))

# Correlations
cat("Correlations (subsample with n≥20):\n")
cor_raw_abund <- cor(analysis_data_rare$otu_richness, analysis_data_rare$total_cafi)
cor_rare_abund <- cor(analysis_data_rare$rarefied_richness, analysis_data_rare$total_cafi)
cat(sprintf("  Raw richness vs abundance: r = %.3f\n", cor_raw_abund))
cat(sprintf("  Rarefied richness vs abundance: r = %.3f\n\n", cor_rare_abund))

# --- Interpretation ---
cat("=== INTERPRETATION ===\n")
if (p_raw_full < 0.1 && p_rare > 0.1) {
  cat("The richness → condition signal is likely an ABUNDANCE ARTIFACT.\n")
  cat("Raw richness shows a trend (p = ", round(p_raw_full, 3), "), but\n", sep = "")
  cat("rarefied richness shows NO relationship (p = ", round(p_rare, 3), ").\n", sep = "")
  cat("The apparent diversity effect disappears after controlling for\n")
  cat("sampling intensity.\n")
  richness_artifact <- TRUE
} else if (p_raw_full < 0.1 && p_rare < 0.1) {
  cat("Richness effect PERSISTS after rarefaction.\n")
  cat("This suggests a true diversity → condition relationship.\n")
  richness_artifact <- FALSE
} else {
  cat("Neither raw nor rarefied richness significantly predicts condition.\n")
  richness_artifact <- NA
}
cat("\n")

# --- Save rarefied richness results to stats_results ---
stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "H-richness-artifact",
    question = "Is richness → condition an abundance artifact?",
    test_name = "Rarefied richness (depth=20)",
    test_statistic = s_rare$coefficients["rarefied_richness", "t value"],
    df = s_rare$df[2],
    p_value = p_rare,
    effect_size = b_rare,
    interpretation = ifelse(isTRUE(richness_artifact), "YES - abundance artifact",
                            ifelse(isFALSE(richness_artifact), "NO - true diversity effect", "Inconclusive")),
    effect_type = "Regression coefficient",
    n = nrow(analysis_data_rare),
    notes = sprintf("Rarefied to n=%d; raw richness p=%.3f, rarefied p=%.3f",
                    rarefy_depth, p_raw_full, p_rare)
  )
)

# --- Create comparison dataframe for plotting ---
richness_comparison_df <- tibble(
  type = c("Raw richness\n(full sample)", "Raw richness\n(n≥20 subset)",
           "Rarefied richness\n(depth=20)"),
  estimate = c(b_raw_full, b_raw_sub, b_rare),
  se = c(s_raw_full$coefficients["otu_richness", "Std. Error"],
         s_raw_sub$coefficients["otu_richness", "Std. Error"],
         s_rare$coefficients["rarefied_richness", "Std. Error"]),
  p_value = c(p_raw_full, p_raw_sub, p_rare),
  n = c(nrow(analysis_data), nrow(analysis_data_rare), nrow(analysis_data_rare)),
  cor_with_abundance = c(
    cor(analysis_data$otu_richness, analysis_data$total_cafi),
    cor_raw_abund, cor_rare_abund)
) %>%
  mutate(
    ci_lower = estimate - qt(0.975, n - 4) * se,
    ci_upper = estimate + qt(0.975, n - 4) * se,
    significant = p_value < 0.05,
    type = factor(type, levels = type)
  )

# Save for use in fig5
save_object(richness_comparison_df, "richness_comparison_results")
save_object(analysis_data_rare, "analysis_data_rarefied")

cat("Saved: richness_comparison_results.rds\n")
cat("Saved: analysis_data_rarefied.rds\n\n")

# --- Rarefied richness sensitivity across multiple depths ---
cat("Rarefied richness sensitivity for condition model:\n")
cat("(Testing robustness to rarefaction depth choice)\n")
depth_sensitivity <- data.frame(
  depth = integer(), p_value = numeric(), beta = numeric(), n_corals = integer(),
  stringsAsFactors = FALSE
)
for (depth in c(10, 15, 20, 25, 30)) {
  # Filter corals with sufficient individuals
  eligible <- rowSums(comm_matrix) >= depth
  eligible_ids <- names(which(eligible))
  eligible_data <- analysis_data %>%
    filter(coral_id %in% eligible_ids)

  if (nrow(eligible_data) < 30) {
    cat(sprintf("  n = %d: too few eligible corals (%d)\n", depth, nrow(eligible_data)))
    next
  }

  comm_eligible <- comm_matrix[eligible_ids, , drop = FALSE]
  set.seed(42)  # Reproducibility for rarefaction sensitivity
  rare_rich_sens <- rarefy(comm_eligible, sample = depth)
  eligible_data$rarefied_richness_sens <- rare_rich_sens[eligible_data$coral_id]

  rare_model_sens <- lm(condition_score ~ rarefied_richness_sens + log_volume + site,
                         data = eligible_data)
  rare_summary_sens <- summary(rare_model_sens)
  rare_p_sens <- rare_summary_sens$coefficients["rarefied_richness_sens", "Pr(>|t|)"]
  rare_beta_sens <- rare_summary_sens$coefficients["rarefied_richness_sens", "Estimate"]
  cat(sprintf("  n = %d: beta = %.3f, p = %.3f (n_corals = %d)\n",
              depth, rare_beta_sens, rare_p_sens, nrow(eligible_data)))
  depth_sensitivity <- bind_rows(depth_sensitivity, data.frame(
    depth = depth, p_value = rare_p_sens, beta = rare_beta_sens,
    n_corals = nrow(eligible_data), stringsAsFactors = FALSE
  ))
}
cat("\n")

# --- Create Figure S10: Rarefaction Depth Sensitivity ---
if (nrow(depth_sensitivity) >= 2) {
  p_depth_sensitivity <- ggplot(depth_sensitivity, aes(x = depth, y = p_value)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "#D55E00") +
    annotate("text", x = max(depth_sensitivity$depth), y = 0.05,
             label = expression(alpha == 0.05), hjust = 1, vjust = -0.5,
             color = "#D55E00", size = 3.5) +
    labs(x = "Rarefaction depth (individuals)", y = "p-value",
         title = "Figure S10: Rarefaction Depth Sensitivity",
         subtitle = "Richness -> condition relationship across rarefaction depths") +
    theme_publication()

  dir.create(file.path(PATHS$figures, "supplement"), showWarnings = FALSE, recursive = TRUE)
  save_figure(p_depth_sensitivity,
             file.path(PATHS$figures, "supplement", "figS10_rarefaction_sensitivity.png"),
             width = 6.7, height = 4.7)
  cat("   Saved: figS10_rarefaction_sensitivity.png\n\n")
} else {
  cat("   Skipped figS10: insufficient rarefaction depths computed\n\n")
}

# ============================================================================
# PART A3: COMMUNITY TRANSFORMATION SENSITIVITY
# ============================================================================
# Tests whether PC1_CAFI -> condition result is robust to community matrix
# transformation choices. Mirrors companion experimental paper (Stier et al.)
# which tested SQRT, SQRT+center+scale, and Hellinger transformations.
# ============================================================================

cat("============================================================\n")
cat("PART A3: COMMUNITY TRANSFORMATION SENSITIVITY\n")
cat("============================================================\n\n")

cat("Testing PC1_CAFI -> condition under different community transforms...\n")
cat("Mirrors experimental companion paper (Stier et al. script 8).\n\n")

# Use comm_matrix already loaded in Part A2 above
# Define transformation configurations
transform_configs <- list(
  list(name = "Hellinger (default)",
       mat = vegan::decostand(comm_matrix, method = "hellinger"),
       center = TRUE, scale = FALSE),
  list(name = "Sqrt + center + scale",
       mat = sqrt(comm_matrix),
       center = TRUE, scale = TRUE),
  list(name = "Sqrt + center",
       mat = sqrt(comm_matrix),
       center = TRUE, scale = FALSE)
)

# Also test with filtered species (Change 4: mirroring experimental filtering)
species_prev <- colSums(comm_matrix > 0)
species_abund <- colSums(comm_matrix)
keep_spp <- names(which(species_prev >= 10 & species_abund >= 10))
comm_filtered <- comm_matrix[, keep_spp, drop = FALSE]
cat(sprintf("  Species filter (>=10 presences AND >=10 individuals): %d of %d OTUs retained\n\n",
    length(keep_spp), ncol(comm_matrix)))

transform_configs[[4]] <- list(
  name = "Hellinger (filtered)",
  mat = vegan::decostand(comm_filtered, method = "hellinger"),
  center = TRUE, scale = FALSE
)

# Run PCA for each transformation and test against condition
transform_results <- list()

for (tc in transform_configs) {
  # PCA — guard against zero-variance columns when scale. = TRUE
  mat_use <- tc$mat
  if (isTRUE(tc$scale)) {
    col_vars <- apply(mat_use, 2, var, na.rm = TRUE)
    zero_var <- col_vars == 0 | is.na(col_vars)
    if (any(zero_var)) {
      mat_use <- mat_use[, !zero_var, drop = FALSE]
      cat(sprintf("  %s: removed %d zero-variance columns before PCA\n", tc$name, sum(zero_var)))
    }
  }
  pca_result <- prcomp(mat_use, center = tc$center, scale. = tc$scale)
  pc1_raw <- pca_result$x[, 1]

  # Align sign: flip if negatively correlated with total abundance
  total_abund <- rowSums(if (ncol(tc$mat) == ncol(comm_matrix)) comm_matrix else comm_filtered)
  if (cor(pc1_raw, total_abund) < 0) {
    pc1_raw <- -pc1_raw
  }

  # Z-score
  pc1_z <- scale(pc1_raw)[, 1]
  pct_var <- round(100 * summary(pca_result)$importance[2, 1], 1)

  # Build temporary dataset for model testing
  pc1_df <- tibble(coral_id = rownames(tc$mat), pc1_transform = pc1_z)
  temp_data <- analysis_data %>%
    left_join(pc1_df, by = "coral_id") %>%
    filter(!is.na(pc1_transform))

  # Test condition ~ PC1_transform + log_volume + site
  tryCatch({
    mod <- lm(condition_score ~ pc1_transform + log_volume + site, data = temp_data)
    coef_tab <- summary(mod)$coefficients
    est <- coef_tab["pc1_transform", "Estimate"]
    se <- coef_tab["pc1_transform", "Std. Error"]
    p_val <- coef_tab["pc1_transform", "Pr(>|t|)"]

    # HC3 robust p-value
    robust_se <- sqrt(diag(sandwich::vcovHC(mod, type = "HC3")))["pc1_transform"]
    robust_t <- est / robust_se
    robust_p <- 2 * pt(abs(robust_t), df = mod$df.residual, lower.tail = FALSE)

    transform_results[[tc$name]] <- tibble(
      transform = tc$name,
      n_species = ncol(tc$mat),
      pct_var_PC1 = pct_var,
      estimate = round(est, 4),
      se = round(se, 4),
      p_value = round(p_val, 4),
      se_robust = round(robust_se, 4),
      p_robust = round(robust_p, 4),
      n = nrow(temp_data)
    )

    cat(sprintf("  %-30s n_spp=%3d  PC1_var=%5.1f%%  beta=%7.4f  p=%6.4f  p_HC3=%6.4f  (n=%d)\n",
        tc$name, ncol(tc$mat), pct_var, est, p_val, robust_p, nrow(temp_data)))
  }, error = function(e) {
    cat(sprintf("  %-30s ERROR: %s\n", tc$name, e$message))
  })
}

# Combine and save
transform_sensitivity_df <- bind_rows(transform_results)

if (nrow(transform_sensitivity_df) > 0) {
  save_table(transform_sensitivity_df, "community_transform_sensitivity")
  cat(sprintf("\n  Saved: community_transform_sensitivity.csv (%d transforms tested)\n",
      nrow(transform_sensitivity_df)))

  # Summary
  all_ns <- all(transform_sensitivity_df$p_robust > 0.05)
  cat(sprintf("  Conclusion: %s\n\n",
      if (all_ns) "PC1_CAFI -> condition is non-significant across ALL transformations (robust null)"
      else "Some transformations show significance -- check details"))
} else {
  cat("  WARNING: No transformation results computed\n\n")
}

# ============================================================================
# PART A4: BEF ANALYSIS — DIVERSITY vs ABUNDANCE EFFECTS ON CONDITION
# ============================================================================
# Tilman's BEF theory predicts diversity per se benefits ecosystem function
# (complementarity/insurance), beyond simple abundance effects.
#
# The rarefaction test (Part A2) is ambiguous: it removes the abundance pathway
# but that pathway might BE the mechanism (more species → fuller habitat use →
# more total CAFI → better condition). Here we directly partition richness vs
# abundance contributions using three complementary approaches:
#
# 1. Partial regression: richness + abundance in same model (VIF check)
# 2. Variance partitioning: unique R² for richness, abundance, and shared
# 3. Path model: volume → richness → condition vs volume → abundance → condition
#
# A PRIORI BEF HYPOTHESES (k=2):
#   H1: Species richness → condition (complementarity)
#   H2: Total abundance → condition (more bodies = more benefit)
# Corrected with Hochberg (FWER) across these 2 a priori tests.
# PC1_CAFI (community composition identity) is a separate question → supplement.
# Remaining 4 predictors (Shannon, Trapezia, Fish, Galeropsis) are exploratory.
# ============================================================================

cat("============================================================\n")
cat("PART A4: BEF ANALYSIS — DIVERSITY vs ABUNDANCE EFFECTS\n")
cat("============================================================\n\n")

cat("Tilman BEF framework: does diversity predict condition beyond abundance?\n")
cat("Three complementary approaches: partial regression, variance partitioning, path model.\n\n")

# --- Setup: prepare data with both richness and abundance ---
bef_data <- analysis_data %>%
  filter(!is.na(condition_score), !is.na(otu_richness), !is.na(total_cafi),
         !is.na(log_volume), !is.na(site)) %>%
  mutate(
    richness_z = scale(otu_richness)[,1],
    abundance_z = scale(sqrt(total_cafi))[,1],  # sqrt to match PART A transform
    log_volume_z = scale(log_volume)[,1],
    site = factor(site)
  )

cat("BEF analysis sample size: n =", nrow(bef_data), "\n")
cat("Richness-abundance correlation: r =",
    round(cor(bef_data$otu_richness, bef_data$total_cafi), 3), "\n")
cat("Richness-abundance correlation (z-scored): r =",
    round(cor(bef_data$richness_z, bef_data$abundance_z), 3), "\n\n")

# ---- 1. PARTIAL REGRESSION: richness controlling for abundance ----
cat("--- 1. PARTIAL REGRESSION ---\n")
cat("Model: condition ~ richness + sqrt(abundance) + log(volume) + site\n\n")

# Check VIF first
m_both <- lm(condition_score ~ otu_richness + sqrt(total_cafi) + log_volume + site,
             data = bef_data)

if (requireNamespace("car", quietly = TRUE)) {
  vif_both <- car::vif(m_both)
  cat("VIF check (richness + abundance in same model):\n")
  for (v in seq_along(vif_both)) {
    cat(sprintf("  %-20s VIF = %.2f %s\n",
                names(vif_both)[v], vif_both[v],
                ifelse(vif_both[v] > 5, "WARNING: HIGH",
                       ifelse(vif_both[v] > 2.5, "(moderate)", ""))))
  }
  cat("\n")
}

# Extract coefficients with HC3 robust SEs
s_both <- summary(m_both)

if (requireNamespace("sandwich", quietly = TRUE) &&
    requireNamespace("lmtest", quietly = TRUE)) {
  robust_both <- lmtest::coeftest(m_both, vcov. = sandwich::vcovHC(m_both, type = "HC3"))

  richness_partial_b <- robust_both["otu_richness", "Estimate"]
  richness_partial_se <- robust_both["otu_richness", "Std. Error"]
  richness_partial_p <- robust_both["otu_richness", "Pr(>|t|)"]

  abundance_partial_b <- robust_both["sqrt(total_cafi)", "Estimate"]
  abundance_partial_se <- robust_both["sqrt(total_cafi)", "Std. Error"]
  abundance_partial_p <- robust_both["sqrt(total_cafi)", "Pr(>|t|)"]

  cat("Partial effects (HC3 robust SEs):\n")
  cat(sprintf("  Richness  | abundance:  beta = %.4f, SE = %.4f, p = %.4f %s\n",
      richness_partial_b, richness_partial_se, richness_partial_p,
      ifelse(richness_partial_p < 0.05, "*", "")))
  cat(sprintf("  Abundance | richness:   beta = %.4f, SE = %.4f, p = %.4f %s\n",
      abundance_partial_b, abundance_partial_se, abundance_partial_p,
      ifelse(abundance_partial_p < 0.05, "*", "")))
} else {
  # Fallback to OLS
  richness_partial_b <- s_both$coefficients["otu_richness", "Estimate"]
  richness_partial_se <- s_both$coefficients["otu_richness", "Std. Error"]
  richness_partial_p <- s_both$coefficients["otu_richness", "Pr(>|t|)"]

  abundance_partial_b <- s_both$coefficients["sqrt(total_cafi)", "Estimate"]
  abundance_partial_se <- s_both$coefficients["sqrt(total_cafi)", "Std. Error"]
  abundance_partial_p <- s_both$coefficients["sqrt(total_cafi)", "Pr(>|t|)"]

  cat("Partial effects (OLS SEs):\n")
  cat(sprintf("  Richness  | abundance:  beta = %.4f, SE = %.4f, p = %.4f %s\n",
      richness_partial_b, richness_partial_se, richness_partial_p,
      ifelse(richness_partial_p < 0.05, "*", "")))
  cat(sprintf("  Abundance | richness:   beta = %.4f, SE = %.4f, p = %.4f %s\n",
      abundance_partial_b, abundance_partial_se, abundance_partial_p,
      ifelse(abundance_partial_p < 0.05, "*", "")))
}

cat(sprintf("\n  Model adj R² = %.4f (compared to richness-only adj R² = %.4f)\n",
    s_both$adj.r.squared,
    summary(lm(condition_score ~ otu_richness + log_volume + site, data = bef_data))$adj.r.squared))

# Standardized betas for comparison
m_both_std <- lm(condition_score ~ richness_z + abundance_z + log_volume_z + site,
                 data = bef_data %>% mutate(condition_score = scale(condition_score)[,1]))
s_both_std <- summary(m_both_std)
cat(sprintf("  Standardized β (richness):  %.3f\n", coef(m_both_std)["richness_z"]))
cat(sprintf("  Standardized β (abundance): %.3f\n\n", coef(m_both_std)["abundance_z"]))

# ---- 2. VARIANCE PARTITIONING ----
cat("--- 2. VARIANCE PARTITIONING ---\n")
cat("Decomposing R² into unique and shared contributions.\n\n")

# Three nested models for hierarchical R² decomposition
m_null    <- lm(condition_score ~ log_volume + site, data = bef_data)
m_rich    <- lm(condition_score ~ otu_richness + log_volume + site, data = bef_data)
m_abund   <- lm(condition_score ~ sqrt(total_cafi) + log_volume + site, data = bef_data)
m_full_vp <- lm(condition_score ~ otu_richness + sqrt(total_cafi) + log_volume + site,
                data = bef_data)

r2_null  <- summary(m_null)$r.squared
r2_rich  <- summary(m_rich)$r.squared
r2_abund <- summary(m_abund)$r.squared
r2_full  <- summary(m_full_vp)$r.squared

# Unique contributions (Type III-style)
unique_richness  <- r2_full - r2_abund   # R² gained by adding richness to abundance model
unique_abundance <- r2_full - r2_rich    # R² gained by adding abundance to richness model
total_explained  <- r2_full - r2_null    # Total R² from both predictors
shared           <- total_explained - unique_richness - unique_abundance

cat(sprintf("  R² (null: volume + site only):       %.4f\n", r2_null))
cat(sprintf("  R² (+ richness only):                %.4f  [Δ = %.4f]\n", r2_rich, r2_rich - r2_null))
cat(sprintf("  R² (+ abundance only):               %.4f  [Δ = %.4f]\n", r2_abund, r2_abund - r2_null))
cat(sprintf("  R² (+ richness + abundance):         %.4f  [Δ = %.4f]\n", r2_full, r2_full - r2_null))
cat(sprintf("\n  Unique to richness:   %.4f (%.1f%% of total explained)\n",
    unique_richness, 100 * unique_richness / max(total_explained, 1e-10)))
cat(sprintf("  Unique to abundance:  %.4f (%.1f%% of total explained)\n",
    unique_abundance, 100 * unique_abundance / max(total_explained, 1e-10)))
cat(sprintf("  Shared (confounded):  %.4f (%.1f%% of total explained)\n",
    shared, 100 * shared / max(total_explained, 1e-10)))
cat(sprintf("  Total explained:      %.4f\n\n", total_explained))

# F-test for unique contributions
anova_rich_added <- anova(m_abund, m_full_vp)   # richness added to abundance model
anova_abund_added <- anova(m_rich, m_full_vp)    # abundance added to richness model

f_rich_unique  <- anova_rich_added$F[2]
p_rich_unique  <- anova_rich_added$`Pr(>F)`[2]
f_abund_unique <- anova_abund_added$F[2]
p_abund_unique <- anova_abund_added$`Pr(>F)`[2]

cat("F-tests for unique contributions (Type III):\n")
cat(sprintf("  Richness  | abundance + volume + site:  F = %.3f, p = %.4f %s\n",
    f_rich_unique, p_rich_unique, ifelse(p_rich_unique < 0.05, "*", "")))
cat(sprintf("  Abundance | richness + volume + site:   F = %.3f, p = %.4f %s\n\n",
    f_abund_unique, p_abund_unique, ifelse(p_abund_unique < 0.05, "*", "")))

# ---- 3. PATH MODEL (piecewiseSEM) ----
cat("--- 3. PATH MODEL (piecewiseSEM) ---\n")
cat("DAG: Volume → Abundance → Condition\n")
cat("     Volume → Richness  → Condition\n")
cat("     (Richness ↔ Abundance correlated via Volume)\n\n")

if (requireNamespace("piecewiseSEM", quietly = TRUE)) {

  # Path models using z-scored predictors for standardized coefficients
  # Model 1: Abundance ~ Volume + Site
  path_abund <- lm(abundance_z ~ log_volume_z + site, data = bef_data)

  # Model 2: Richness ~ Volume + Site
  path_rich <- lm(richness_z ~ log_volume_z + site, data = bef_data)

  # Model 3: Condition ~ Richness + Abundance + Volume + Site
  path_cond <- lm(condition_score ~ richness_z + abundance_z + log_volume_z + site,
                  data = bef_data)

  psem_mod <- tryCatch(
    piecewiseSEM::psem(path_abund, path_rich, path_cond),
    error = function(e) {
      cat("  piecewiseSEM::psem() failed:", e$message, "\n")
      cat("  Falling back to manual path coefficient extraction.\n\n")
      NULL
    }
  )

  if (is.null(psem_mod)) {
    # Manual fallback: extract standardized path coefficients from component models
    cat("Manual path coefficients (standardized via z-scoring):\n")

    s_abund <- summary(path_abund)
    s_rich  <- summary(path_rich)
    s_cond  <- summary(path_cond)

    cat(sprintf("  log_volume → abundance:  β = %.4f, p = %.4f\n",
        coef(path_abund)["log_volume_z"],
        s_abund$coefficients["log_volume_z", "Pr(>|t|)"]))
    cat(sprintf("  log_volume → richness:   β = %.4f, p = %.4f\n",
        coef(path_rich)["log_volume_z"],
        s_rich$coefficients["log_volume_z", "Pr(>|t|)"]))
    cat(sprintf("  richness   → condition:  β = %.4f, p = %.4f\n",
        coef(path_cond)["richness_z"],
        s_cond$coefficients["richness_z", "Pr(>|t|)"]))
    cat(sprintf("  abundance  → condition:  β = %.4f, p = %.4f\n",
        coef(path_cond)["abundance_z"],
        s_cond$coefficients["abundance_z", "Pr(>|t|)"]))
    cat(sprintf("  log_volume → condition:  β = %.4f, p = %.4f\n",
        coef(path_cond)["log_volume_z"],
        s_cond$coefficients["log_volume_z", "Pr(>|t|)"]))

    # Indirect effects
    indirect_via_rich <- coef(path_rich)["log_volume_z"] * coef(path_cond)["richness_z"]
    indirect_via_abund <- coef(path_abund)["log_volume_z"] * coef(path_cond)["abundance_z"]
    cat(sprintf("\n  Indirect (Volume → Richness → Condition):  %.4f\n", indirect_via_rich))
    cat(sprintf("  Indirect (Volume → Abundance → Condition): %.4f\n", indirect_via_abund))

    # R² for endogenous
    cat(sprintf("\n  R² (abundance model): %.4f\n", s_abund$r.squared))
    cat(sprintf("  R² (richness model):  %.4f\n", s_rich$r.squared))
    cat(sprintf("  R² (condition model): %.4f\n", s_cond$r.squared))

    # Save manual path coefficients
    psem_coefs <- tibble(
      Predictor = c("log_volume_z", "log_volume_z", "richness_z", "abundance_z", "log_volume_z"),
      Response = c("abundance_z", "richness_z", "condition_score", "condition_score", "condition_score"),
      Std.Estimate = c(coef(path_abund)["log_volume_z"], coef(path_rich)["log_volume_z"],
                       coef(path_cond)["richness_z"], coef(path_cond)["abundance_z"],
                       coef(path_cond)["log_volume_z"]),
      P.Value = c(s_abund$coefficients["log_volume_z", "Pr(>|t|)"],
                  s_rich$coefficients["log_volume_z", "Pr(>|t|)"],
                  s_cond$coefficients["richness_z", "Pr(>|t|)"],
                  s_cond$coefficients["abundance_z", "Pr(>|t|)"],
                  s_cond$coefficients["log_volume_z", "Pr(>|t|)"])
    )
    save_table(psem_coefs, "bef_path_coefficients")
    cat("\nSaved: bef_path_coefficients.csv (manual extraction)\n")
  }

  if (!is.null(psem_mod)) {

  # Capture summary
  psem_summary <- tryCatch(
    summary(psem_mod, .progressBar = FALSE),
    error = function(e) {
      cat("  piecewiseSEM summary failed:", e$message, "\n")
      cat("  Falling back to manual path extraction.\n\n")

      # Manual fallback (same as above but within psem_mod success path)
      s_abund_fb <- summary(path_abund)
      s_rich_fb  <- summary(path_rich)
      s_cond_fb  <- summary(path_cond)

      cat("Manual path coefficients (standardized via z-scoring):\n")
      cat(sprintf("  log_volume → abundance:  β = %.4f, p = %.4f\n",
          coef(path_abund)["log_volume_z"],
          s_abund_fb$coefficients["log_volume_z", "Pr(>|t|)"]))
      cat(sprintf("  log_volume → richness:   β = %.4f, p = %.4f\n",
          coef(path_rich)["log_volume_z"],
          s_rich_fb$coefficients["log_volume_z", "Pr(>|t|)"]))
      cat(sprintf("  richness   → condition:  β = %.4f, p = %.4f\n",
          coef(path_cond)["richness_z"],
          s_cond_fb$coefficients["richness_z", "Pr(>|t|)"]))
      cat(sprintf("  abundance  → condition:  β = %.4f, p = %.4f\n",
          coef(path_cond)["abundance_z"],
          s_cond_fb$coefficients["abundance_z", "Pr(>|t|)"]))
      cat(sprintf("  log_volume → condition:  β = %.4f, p = %.4f\n",
          coef(path_cond)["log_volume_z"],
          s_cond_fb$coefficients["log_volume_z", "Pr(>|t|)"]))

      indirect_r <- coef(path_rich)["log_volume_z"] * coef(path_cond)["richness_z"]
      indirect_a <- coef(path_abund)["log_volume_z"] * coef(path_cond)["abundance_z"]
      cat(sprintf("\n  Indirect (Volume → Richness → Condition):  %.4f\n", indirect_r))
      cat(sprintf("  Indirect (Volume → Abundance → Condition): %.4f\n", indirect_a))

      cat(sprintf("\n  R² (abundance model): %.4f\n", s_abund_fb$r.squared))
      cat(sprintf("  R² (richness model):  %.4f\n", s_rich_fb$r.squared))
      cat(sprintf("  R² (condition model): %.4f\n", s_cond_fb$r.squared))

      psem_coefs_fb <- tibble(
        Predictor = c("log_volume_z", "log_volume_z", "richness_z", "abundance_z", "log_volume_z"),
        Response = c("abundance_z", "richness_z", "condition_score", "condition_score", "condition_score"),
        Std.Estimate = c(coef(path_abund)["log_volume_z"], coef(path_rich)["log_volume_z"],
                         coef(path_cond)["richness_z"], coef(path_cond)["abundance_z"],
                         coef(path_cond)["log_volume_z"]),
        P.Value = c(s_abund_fb$coefficients["log_volume_z", "Pr(>|t|)"],
                    s_rich_fb$coefficients["log_volume_z", "Pr(>|t|)"],
                    s_cond_fb$coefficients["richness_z", "Pr(>|t|)"],
                    s_cond_fb$coefficients["abundance_z", "Pr(>|t|)"],
                    s_cond_fb$coefficients["log_volume_z", "Pr(>|t|)"])
      )
      save_table(psem_coefs_fb, "bef_path_coefficients")
      cat("\nSaved: bef_path_coefficients.csv (manual extraction)\n")

      NULL
    }
  )

  if (!is.null(psem_summary)) {

  cat("Fisher's C =", round(psem_summary$Cstat$Fisher.C, 3),
      ", df =", psem_summary$Cstat$df,
      ", p =", round(psem_summary$Cstat$P.Value, 4),
      ifelse(psem_summary$Cstat$P.Value > 0.05, " (good fit)", " (poor fit — missing paths)"), "\n\n")

  # Extract standardized coefficients
  psem_coefs <- piecewiseSEM::coefs(psem_mod, standardize = "scale")

  cat("Standardized path coefficients:\n")
  for (i in 1:nrow(psem_coefs)) {
    row <- psem_coefs[i, ]
    sig <- ""
    if (!is.na(row$P.Value)) {
      if (row$P.Value < 0.001) sig <- "***"
      else if (row$P.Value < 0.01) sig <- "**"
      else if (row$P.Value < 0.05) sig <- "*"
      else if (row$P.Value < 0.1) sig <- "."
    }
    cat(sprintf("  %-15s → %-15s  β = %7.4f, p = %s %s\n",
        row$Predictor, row$Response,
        row$Std.Estimate,
        ifelse(is.na(row$P.Value), "  NA  ", sprintf("%.4f", row$P.Value)),
        sig))
  }

  # R² for each endogenous variable
  cat("\nR² for endogenous variables:\n")
  for (i in 1:nrow(psem_summary$R2)) {
    cat(sprintf("  %-20s R² = %.4f\n",
        psem_summary$R2$Response[i], psem_summary$R2$R.squared[i]))
  }

  # Extract the key path: richness → condition (direct BEF effect)
  rich_to_cond_row <- psem_coefs[psem_coefs$Predictor == "richness_z" &
                                   psem_coefs$Response == "condition_score", ]
  abund_to_cond_row <- psem_coefs[psem_coefs$Predictor == "abundance_z" &
                                    psem_coefs$Response == "condition_score", ]

  if (nrow(rich_to_cond_row) > 0) {
    cat(sprintf("\nKey BEF path (richness → condition): β = %.4f, p = %.4f\n",
        rich_to_cond_row$Std.Estimate, rich_to_cond_row$P.Value))
  }
  if (nrow(abund_to_cond_row) > 0) {
    cat(sprintf("Abundance path (abundance → condition): β = %.4f, p = %.4f\n",
        abund_to_cond_row$Std.Estimate, abund_to_cond_row$P.Value))
  }

  # Indirect effects (manual calculation)
  vol_to_rich <- psem_coefs[psem_coefs$Predictor == "log_volume_z" &
                              psem_coefs$Response == "richness_z", ]
  vol_to_abund <- psem_coefs[psem_coefs$Predictor == "log_volume_z" &
                               psem_coefs$Response == "abundance_z", ]

  if (nrow(vol_to_rich) > 0 && nrow(rich_to_cond_row) > 0) {
    indirect_via_rich <- vol_to_rich$Std.Estimate * rich_to_cond_row$Std.Estimate
    cat(sprintf("\nIndirect effect (Volume → Richness → Condition): %.4f\n", indirect_via_rich))
  }
  if (nrow(vol_to_abund) > 0 && nrow(abund_to_cond_row) > 0) {
    indirect_via_abund <- vol_to_abund$Std.Estimate * abund_to_cond_row$Std.Estimate
    cat(sprintf("Indirect effect (Volume → Abundance → Condition): %.4f\n", indirect_via_abund))
  }

  # Save path coefficients
  save_table(psem_coefs, "bef_path_coefficients")
  cat("\nSaved: bef_path_coefficients.csv\n")

  }  # end if (!is.null(psem_summary))
  }  # end if (!is.null(psem_mod))

} else {
  cat("piecewiseSEM not available. Install with: install.packages('piecewiseSEM')\n")
  cat("Skipping path model analysis.\n")
}

# ---- 4. A PRIORI HYPOTHESIS CORRECTION (Hochberg) ----
cat("\n--- 4. A PRIORI BEF HYPOTHESIS CORRECTION ---\n")
cat("Two a priori hypotheses from BEF theory (Hochberg FWER, k=2):\n")
cat("  H1: Species richness → condition (complementarity)\n")
cat("  H2: Total abundance → condition (more bodies = more benefit)\n")
cat("PC1_CAFI (community composition identity) is a separate question → supplement.\n\n")

# Extract raw p-values for the 2 a priori BEF predictors from PART A results
# Using OLS p-values (primary inference — BP test confirms no heteroscedasticity)
apriori_predictors <- c("Species richness", "Total CAFI")
apriori_rows <- cafi_to_condition_df %>%
  filter(predictor %in% apriori_predictors)

apriori_p_raw <- setNames(apriori_rows$p_value, apriori_rows$predictor)

# Hochberg correction (step-up procedure, uniformly more powerful than Bonferroni)
apriori_p_hochberg <- p.adjust(apriori_p_raw, method = "hochberg")

# HC3 sensitivity (supplement)
apriori_p_raw_hc3 <- setNames(apriori_rows$p_value_robust, apriori_rows$predictor)
apriori_p_hochberg_hc3 <- p.adjust(apriori_p_raw_hc3, method = "hochberg")

cat("A priori BEF predictors (Hochberg FWER, k=2):\n")
for (nm in names(apriori_p_raw)) {
  sig <- ifelse(apriori_p_hochberg[nm] < 0.05, " *SIGNIFICANT*", "")
  cat(sprintf("  %-30s p_raw = %.4f, p_Hochberg = %.4f%s\n",
      nm, apriori_p_raw[nm], apriori_p_hochberg[nm], sig))
}

# Supplement: PC1_CAFI composition test (separate question, no correction needed)
pc1_row <- cafi_to_condition_df %>% filter(predictor == "PC1_CAFI (community)")
cat(sprintf("\nSupplement composition test (uncorrected):\n  %-30s p = %.4f (NS)\n",
    "PC1_CAFI (community)", pc1_row$p_value))

# Exploratory predictors (BH-FDR, separate family)
exploratory_predictors <- c("Shannon diversity", "Trapezia abundance",
                            "Resident Fish abundance", "Galeropsis abundance")
exploratory_rows <- cafi_to_condition_df %>%
  filter(predictor %in% exploratory_predictors)

exploratory_p_raw <- setNames(exploratory_rows$p_value, exploratory_rows$predictor)
exploratory_p_fdr <- p.adjust(exploratory_p_raw, method = "BH")

cat("\nExploratory predictors (BH-FDR, k=4):\n")
for (nm in names(exploratory_p_raw)) {
  sig <- ifelse(exploratory_p_fdr[nm] < 0.05, " *SIGNIFICANT*", "")
  cat(sprintf("  %-30s p_raw = %.4f, p_FDR = %.4f%s\n",
      nm, exploratory_p_raw[nm], exploratory_p_fdr[nm], sig))
}
cat("\n")

# ---- Save BEF analysis results ----
bef_results <- tibble(
  analysis = c("Partial regression (richness | abundance)",
               "Partial regression (abundance | richness)",
               "Variance partition: unique richness",
               "Variance partition: unique abundance",
               "Variance partition: shared",
               "F-test: richness unique",
               "F-test: abundance unique"),
  estimate = c(richness_partial_b, abundance_partial_b,
               unique_richness, unique_abundance, shared,
               f_rich_unique, f_abund_unique),
  se = c(richness_partial_se, abundance_partial_se,
         NA, NA, NA, NA, NA),
  p_value = c(richness_partial_p, abundance_partial_p,
              NA, NA, NA,
              p_rich_unique, p_abund_unique),
  notes = c("HC3 robust SE", "HC3 robust SE",
            sprintf("%.1f%% of total", 100 * unique_richness / max(total_explained, 1e-10)),
            sprintf("%.1f%% of total", 100 * unique_abundance / max(total_explained, 1e-10)),
            sprintf("%.1f%% of total", 100 * shared / max(total_explained, 1e-10)),
            "Type III F-test", "Type III F-test")
)

save_table(bef_results, "bef_diversity_abundance_partition")
cat("Saved: bef_diversity_abundance_partition.csv\n")

# Save a priori correction results
apriori_correction <- tibble(
  predictor = names(apriori_p_raw),
  p_raw = apriori_p_raw,
  p_hochberg = apriori_p_hochberg,
  significant = apriori_p_hochberg < 0.05,
  hypothesis_type = "a_priori"
) %>%
  bind_rows(tibble(
    predictor = names(exploratory_p_raw),
    p_raw = exploratory_p_raw,
    p_hochberg = exploratory_p_fdr,  # FDR for exploratory
    significant = exploratory_p_fdr < 0.05,
    hypothesis_type = "exploratory"
  ))

save_table(apriori_correction, "bef_hypothesis_corrections")
cat("Saved: bef_hypothesis_corrections.csv\n")

# ---- Interpretation ----
cat("\n=== BEF ANALYSIS INTERPRETATION ===\n")
if (richness_partial_p < 0.05 && abundance_partial_p >= 0.05) {
  cat("COMPLEMENTARITY: Richness predicts condition beyond abundance.\n")
  cat("This supports Tilman's BEF hypothesis: species identity matters.\n")
} else if (richness_partial_p >= 0.05 && abundance_partial_p < 0.05) {
  cat("ABUNDANCE ONLY: Abundance predicts condition, richness does not.\n")
  cat("The richness signal is mediated entirely through abundance.\n")
} else if (richness_partial_p < 0.05 && abundance_partial_p < 0.05) {
  cat("BOTH: Richness and abundance independently predict condition.\n")
  cat("Both complementarity and abundance pathways are active.\n")
} else {
  cat("NEITHER: After mutual adjustment, neither richness nor abundance\n")
  cat("independently predicts condition. The shared variance (collinearity)\n")
  cat("prevents clean separation, but the a priori richness signal\n")
  cat("(Hochberg p = ", round(apriori_p_hochberg["Species richness"], 4), ") ",
      ifelse(apriori_p_hochberg["Species richness"] < 0.05,
             "survives targeted correction.", "does not survive correction."),
      "\n", sep = "")
}
cat("\n")

# ============================================================================
# PART B: CONDITION -> CAFI EFFECTS (REVERSE DIRECTION)
# ============================================================================

cat("============================================================\n")
cat("PART B: CONDITION -> CAFI EFFECTS (REVERSE DIRECTION TEST)\n")
cat("============================================================\n\n")

cat("Testing if healthier corals attract more CAFI...\n")
cat("Model: CAFI_metric ~ Condition + log(Volume) + Site\n\n")

# Function for reverse models (CAFI as response)
run_reverse_model <- function(data, response_name, response_col, use_nb = TRUE) {
  formula_str <- paste(response_col, "~ condition_score + log_volume + site")

  tryCatch({
    if (use_nb && response_col != "shannon") {
      # Negative binomial for count data
      model <- MASS::glm.nb(as.formula(formula_str), data = data)
      model_summary <- summary(model)

      estimate <- coef(model)["condition_score"]
      se <- model_summary$coefficients["condition_score", "Std. Error"]
      z_val <- model_summary$coefficients["condition_score", "z value"]
      p_val <- model_summary$coefficients["condition_score", "Pr(>|z|)"]

      stat_type <- "z"
      stat_val <- z_val

    } else {
      # Linear model for continuous (Shannon)
      model <- lm(as.formula(formula_str), data = data)
      model_summary <- summary(model)

      estimate <- coef(model)["condition_score"]
      se <- model_summary$coefficients["condition_score", "Std. Error"]
      t_val <- model_summary$coefficients["condition_score", "t value"]
      p_val <- model_summary$coefficients["condition_score", "Pr(>|t|)"]

      stat_type <- "t"
      stat_val <- t_val
    }

    # Use t-quantile for linear models (small-sample correction), z for NB GLMs
    if (stat_type == "t") {
      df_resid <- model_summary$df[2]
      crit_val <- qt(0.975, df_resid)
    } else {
      crit_val <- 1.96  # Wald z for NB GLM
    }
    ci_lower <- estimate - crit_val * se
    ci_upper <- estimate + crit_val * se

    result <- data.frame(
      direction = "Condition -> CAFI",
      response = response_name,
      estimate = estimate,
      se = se,
      stat_type = stat_type,
      stat_value = stat_val,
      p_value = p_val,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      n = nrow(data),
      significant = p_val < 0.05
    )

    return(list(model = model, result = result))

  }, error = function(e) {
    cat("   Warning: Model failed for", response_name, "-", e$message, "\n")
    return(NULL)
  })
}

# Test condition as predictor of each CAFI metric
# PC1_CAFI tests overall community composition response to condition
cafi_responses <- list(
  c("PC1_CAFI (community)", "pc1_cafi"),  # Primary response - composite community score
  c("Total CAFI", "total_cafi"),
  c("Trapezia abundance", "n_trapezia"),
  c("Resident Fish abundance", "n_resident_fish"),
  c("Galeropsis abundance", "n_galeropsis"),
  c("Species richness", "otu_richness"),
  c("Shannon diversity", "shannon")
)

condition_to_cafi_results <- list()
condition_to_cafi_models <- list()

for (resp in cafi_responses) {
  cat("   Testing:", resp[1], "...\n")
  # Use linear model for continuous variables (Shannon, PC1_CAFI), NB for counts
  use_nb <- !(resp[2] %in% c("shannon", "pc1_cafi", "pc2_cafi"))
  result <- run_reverse_model(analysis_data, resp[1], resp[2], use_nb)

  if (!is.null(result)) {
    condition_to_cafi_results[[resp[1]]] <- result$result
    condition_to_cafi_models[[resp[1]]] <- result$model

    r <- result$result
    sig_star <- ifelse(r$p_value < 0.001, "***",
                       ifelse(r$p_value < 0.01, "**",
                              ifelse(r$p_value < 0.05, "*",
                                     ifelse(r$p_value < 0.1, ".", ""))))
    cat("      beta =", round(r$estimate, 4),
        ", SE =", round(r$se, 4),
        ",", r$stat_type, "=", round(r$stat_value, 2),
        ", p =", format.pval(r$p_value, 3), sig_star, "\n")
  }
}

# Combine Condition -> CAFI results
condition_to_cafi_df <- bind_rows(condition_to_cafi_results)

# Apply FDR correction for multiple testing (7 reverse tests)
# NOTE: Family-wise FDR is conservative (corrects within predictor family).
# Global FDR (all tests pooled) is also computed below as sensitivity check.
# Results reported as "nominal" must explicitly note global FDR threshold.
# (Same family-wise approach as CAFI→Condition; see note above)
# Note: reverse models (NB GLMs) don't compute HC3 robust SEs, so use OLS p-values
condition_to_cafi_df <- condition_to_cafi_df %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "BH"),
    significant_fdr = p_fdr < 0.05,
    significant = p_value < 0.05
  )

cat("\nFDR-corrected significance (Benjamini-Hochberg, within Condition→CAFI family):\n")
cat("(See Part I below for global FDR correction across all test families.)\n")
for (i in 1:nrow(condition_to_cafi_df)) {
  row <- condition_to_cafi_df[i, ]
  fdr_star <- ifelse(row$significant_fdr, " *FDR-SIG*", "")
  cat(sprintf("   %-30s p = %.4f, p_FDR = %.4f%s\n",
              row$response, row$p_value, row$p_fdr, fdr_star))
}
cat("\n")

# ============================================================================
# PART C: FUNCTIONAL GROUP COMPARISON
# ============================================================================

cat("============================================================\n")
cat("PART C: FUNCTIONAL GROUP EFFECTS COMPARISON\n")
cat("============================================================\n\n")

# Extract functional group effects for comparison
functional_effects <- cafi_to_condition_df %>%
  filter(predictor %in% c("Trapezia abundance", "Resident Fish abundance",
                          "Galeropsis abundance")) %>%
  mutate(
    functional_group = case_when(
      predictor == "Trapezia abundance" ~ "Trapezia\n(Crabs)",
      predictor == "Resident Fish abundance" ~ "Fish",
      predictor == "Galeropsis abundance" ~ "Galeropsis"
    ),
    expected_sign = case_when(
      predictor == "Trapezia abundance" ~ "positive",
      predictor == "Resident Fish abundance" ~ "positive",
      predictor == "Galeropsis abundance" ~ "negative"
    ),
    observed_sign = ifelse(estimate > 0, "positive", "negative"),
    matches_hypothesis = expected_sign == observed_sign
  )

cat("Functional Group Effect Sizes (on Coral Condition):\n")
cat("---------------------------------------------------\n")
for (i in 1:nrow(functional_effects)) {
  fg <- functional_effects[i, ]
  direction <- ifelse(fg$estimate > 0, "POSITIVE", "NEGATIVE")
  expected <- toupper(fg$expected_sign)
  match <- ifelse(fg$matches_hypothesis, "MATCHES", "CONTRADICTS")
  sig <- ifelse(fg$significant, "(SIGNIFICANT)", "(not significant)")

  cat(sprintf("   %-25s: beta = %+.4f, 95%% CI [%.4f, %.4f], %s hypothesis %s\n",
              fg$predictor, fg$estimate, fg$ci_lower, fg$ci_upper, match, sig))
}
cat("\n")

# ============================================================================
# PART D: CONDITION x VOLUME INTERACTION
# ============================================================================

cat("============================================================\n")
cat("PART D: CONDITION x VOLUME INTERACTION\n")
cat("============================================================\n\n")

cat("Testing: Does the effect of condition on CAFI depend on coral size?\n\n")

# Test interaction for total CAFI
m_interaction <- tryCatch({
  MASS::glm.nb(total_cafi ~ condition_score * log_volume + site, data = analysis_data)
}, error = function(e) {
  cat("   Interaction model failed:", e$message, "\n")
  NULL
})

if (!is.null(m_interaction)) {
  int_summary <- summary(m_interaction)
  int_coef <- coef(m_interaction)

  # Extract interaction term
  int_term <- "condition_score:log_volume"
  if (int_term %in% names(int_coef)) {
    int_estimate <- int_coef[int_term]
    int_se <- int_summary$coefficients[int_term, "Std. Error"]
    int_z <- int_summary$coefficients[int_term, "z value"]
    int_p <- int_summary$coefficients[int_term, "Pr(>|z|)"]

    cat("   Interaction (Condition x Volume):\n")
    cat("      beta =", round(int_estimate, 4), "\n")
    cat("      SE =", round(int_se, 4), "\n")
    cat("      z =", round(int_z, 2), "\n")
    cat("      p =", format.pval(int_p, 3), "\n")
    cat("      Interpretation:",
        ifelse(int_p < 0.05, "Size MODIFIES condition-CAFI relationship",
               "No significant interaction"), "\n\n")

    # Add to results
    stats_results <- bind_rows(stats_results,
      create_result_row(
        hypothesis = "H-interaction",
        question = "Does condition effect on CAFI depend on coral size?",
        test_name = "Negative binomial GLM with interaction",
        test_statistic = int_z,
        df = "1",
        p_value = int_p,
        effect_size = int_estimate,
        effect_type = "Interaction coefficient",
        n = n_complete,
        notes = paste("Condition x log(Volume) interaction:",
                      ifelse(int_p < 0.05, "significant", "not significant"))
      )
    )
  }
}

# ============================================================================
# PART E: PATH ANALYSIS (MEDIATION)
# ============================================================================

cat("============================================================\n")
cat("PART E: PATH ANALYSIS (MEDIATION)\n")
cat("============================================================\n\n")

# Mediation analysis (replaces saturated SEM which had df=0 and omitted site)
# Tests: Volume -> CAFI -> Condition mediation pathway
cat("Testing mediation: Volume -> CAFI -> Condition\n")
cat("NOTE: Mediation cannot establish causation from cross-sectional data.\n")
cat("Results are descriptive of indirect associations.\n")
cat("Sample size:", n_complete, "\n\n")

if (requireNamespace("mediation", quietly = TRUE) && n_complete >= 50) {
  # Prepare standardized data with site
  mediation_data <- analysis_data %>%
    mutate(
      log_volume_z = scale(log_volume)[,1],
      total_cafi_z = scale(total_cafi)[,1]
    ) %>%
    filter(!is.na(condition_score), !is.na(total_cafi_z), !is.na(log_volume_z)) %>%
    drop_na(site)

  # Mediator model: CAFI ~ Volume + Site
  med_model <- lm(total_cafi_z ~ log_volume_z + site, data = mediation_data)
  # Outcome model: Condition ~ CAFI + Volume + Site
  out_model <- lm(condition_score ~ total_cafi_z + log_volume_z + site, data = mediation_data)

  cat("Mediator model: total_cafi_z ~ log_volume_z + site\n")
  cat("Outcome model:  condition_score ~ total_cafi_z + log_volume_z + site\n\n")

  set.seed(42)
  tryCatch({
    med_result <- mediation::mediate(med_model, out_model,
                                      treat = "log_volume_z",
                                      mediator = "total_cafi_z",
                                      boot = TRUE, sims = 1000)

    med_summary <- summary(med_result)
    cat("Mediation Analysis Results:\n")
    cat("---------------------------\n")
    cat("   ACME (Average Causal Mediation Effect):", round(med_result$d0, 4),
        ", 95% CI [", round(med_result$d0.ci[1], 4), ",", round(med_result$d0.ci[2], 4), "]",
        ", p =", format.pval(med_result$d0.p, 3), "\n")
    cat("   ADE (Average Direct Effect):           ", round(med_result$z0, 4),
        ", 95% CI [", round(med_result$z0.ci[1], 4), ",", round(med_result$z0.ci[2], 4), "]",
        ", p =", format.pval(med_result$z0.p, 3), "\n")
    cat("   Total Effect:                          ", round(med_result$tau.coef, 4),
        ", 95% CI [", round(med_result$tau.ci[1], 4), ",", round(med_result$tau.ci[2], 4), "]",
        ", p =", format.pval(med_result$tau.p, 3), "\n")
    cat("   Proportion Mediated:                   ", round(med_result$n0, 4),
        ", 95% CI [", round(med_result$n0.ci[1], 4), ",", round(med_result$n0.ci[2], 4), "]",
        ", p =", format.pval(med_result$n0.p, 3), "\n\n")

    # Save mediation results
    mediation_results_clean <- tibble(
      path = c("ACME (indirect)", "ADE (direct)", "Total", "Proportion mediated"),
      estimate = c(med_result$d0, med_result$z0, med_result$tau.coef, med_result$n0),
      ci_lower = c(med_result$d0.ci[1], med_result$z0.ci[1], med_result$tau.ci[1], med_result$n0.ci[1]),
      ci_upper = c(med_result$d0.ci[2], med_result$z0.ci[2], med_result$tau.ci[2], med_result$n0.ci[2]),
      p_value = c(med_result$d0.p, med_result$z0.p, med_result$tau.p, med_result$n0.p)
    )
    save_table(mediation_results_clean, "path_analysis")

    # Interpretation
    if (med_result$d0.p < 0.05) {
      cat("INTERPRETATION: Significant mediation detected.\n")
      cat("   CAFI partially mediates the Volume -> Condition relationship.\n\n")
    } else {
      cat("INTERPRETATION: No significant mediation.\n")
      cat("   The Volume -> Condition relationship is not significantly mediated by CAFI.\n\n")
    }

  }, error = function(e) {
    cat("Mediation analysis failed:", e$message, "\n\n")
  })

} else if (!requireNamespace("mediation", quietly = TRUE)) {
  cat("Mediation analysis skipped: mediation package not available\n")
  cat("Install with: install.packages('mediation')\n\n")

  # Fallback: manual Sobel test approximation
  cat("Fallback: Manual indirect effect calculation (Volume -> CAFI -> Condition)\n")
  mediation_data <- analysis_data %>%
    mutate(
      log_volume_z = scale(log_volume)[,1],
      total_cafi_z = scale(total_cafi)[,1]
    ) %>%
    filter(!is.na(condition_score), !is.na(total_cafi_z), !is.na(log_volume_z))

  m_a <- lm(total_cafi_z ~ log_volume_z + site, data = mediation_data)
  m_b <- lm(condition_score ~ total_cafi_z + log_volume_z + site, data = mediation_data)
  a <- coef(m_a)["log_volume_z"]
  b <- coef(m_b)["total_cafi_z"]
  indirect_effect <- a * b
  cat("   a (Volume -> CAFI):", round(a, 4), "\n")
  cat("   b (CAFI -> Condition | Volume):", round(b, 4), "\n")
  cat("   Indirect (a*b):", round(indirect_effect, 4), "\n\n")
} else {
  cat("Mediation analysis skipped: Sample size (n =", n_complete, ") too small\n")
  cat("Minimum recommended: n = 50 for reliable mediation analysis\n\n")
}

# ============================================================================
# PART F: KEY SPECIES FROM EXPERIMENTAL PAPER (Stier et al.)
# ============================================================================
# Testing whether species identified in experimental paper show same effects
# in observational survey data
#
# FROM EXPERIMENT (Stier et al.):
#   POSITIVE effects on coral condition:
#     - Caracanthus maculatus (velvetfish)
#     - Alpheus spp. (snapping shrimp, especially A. lottini)
#   NEGATIVE effects on coral condition:
#     - Luniella pugil (xanthid crab)
#     - Cymo quadrilobatus (xanthid crab)
# ============================================================================

cat("============================================================\n")
cat("PART F: KEY SPECIES FROM EXPERIMENTAL PAPER\n")
cat("============================================================\n\n")

cat("Testing species with known effects from Stier et al. experimental study\n\n")

# F.1 Extract key species abundances from cafi_clean
cat("F.1 Extracting key species abundances per coral...\n\n")

key_species_counts <- cafi_clean %>%
  group_by(coral_id) %>%
  summarise(
    # POSITIVE effect species (from experiment)
    n_caracanthus = sum(grepl("Caracanthus", otu, ignore.case = TRUE), na.rm = TRUE),
    n_alpheus_lottini = sum(otu == "Alpheus lottini", na.rm = TRUE),
    n_alpheus_all = sum(grepl("^Alpheus", otu), na.rm = TRUE),  # All Alpheus species
    n_harpiliopsis = sum(otu == "Harpiliopsis spinigera", na.rm = TRUE),
    n_periclimenes = sum(otu == "Periclimenes watamuae", na.rm = TRUE),
    n_calcinus = sum(otu == "Calcinus latens", na.rm = TRUE),

    # NEGATIVE effect species (from experiment)
    n_cymo = sum(grepl("Cymo", otu, ignore.case = TRUE), na.rm = TRUE),
    n_luniella = sum(grepl("Luniella", otu, ignore.case = TRUE), na.rm = TRUE),
    n_xanthid_harmful = sum(grepl("Cymo|Luniella", otu, ignore.case = TRUE), na.rm = TRUE),
    n_alpheus_diadema = sum(otu == "Alpheus diadema", na.rm = TRUE),

    .groups = "drop"
  )

# Merge with analysis data (use separate object to avoid mutating primary dataframe)
analysis_data_species <- analysis_data %>%
  left_join(key_species_counts, by = "coral_id")

# Replace NA with 0 for species counts
key_species_cols <- c("n_caracanthus", "n_alpheus_lottini", "n_alpheus_all",
                      "n_harpiliopsis", "n_periclimenes", "n_calcinus",
                      "n_cymo", "n_luniella", "n_xanthid_harmful",
                      "n_alpheus_diadema")
for (col in key_species_cols) {
  analysis_data_species[[col]] <- ifelse(is.na(analysis_data_species[[col]]), 0, analysis_data_species[[col]])
}

# Summary of key species presence
cat("   Key Species Presence in Condition Dataset (n =", nrow(analysis_data_species), "corals):\n")
cat("   POSITIVE effect species (predicted from experiment):\n")
cat("      Caracanthus maculatus:", sum(analysis_data_species$n_caracanthus > 0), "corals,",
    sum(analysis_data_species$n_caracanthus), "individuals\n")
cat("      Alpheus lottini:", sum(analysis_data_species$n_alpheus_lottini > 0), "corals,",
    sum(analysis_data_species$n_alpheus_lottini), "individuals\n")
cat("      All Alpheus spp.:", sum(analysis_data_species$n_alpheus_all > 0), "corals,",
    sum(analysis_data_species$n_alpheus_all), "individuals\n")
cat("      Harpiliopsis spinigera:", sum(analysis_data_species$n_harpiliopsis > 0), "corals,",
    sum(analysis_data_species$n_harpiliopsis), "individuals\n")
cat("      Periclimenes watamuae:", sum(analysis_data_species$n_periclimenes > 0), "corals,",
    sum(analysis_data_species$n_periclimenes), "individuals\n")
cat("      Calcinus latens:", sum(analysis_data_species$n_calcinus > 0), "corals,",
    sum(analysis_data_species$n_calcinus), "individuals\n")
cat("   NEGATIVE effect species (predicted from experiment):\n")
cat("      Cymo quadrilobatus:", sum(analysis_data_species$n_cymo > 0), "corals,",
    sum(analysis_data_species$n_cymo), "individuals\n")
cat("      Luniella pugil:", sum(analysis_data_species$n_luniella > 0), "corals,",
    sum(analysis_data_species$n_luniella), "individuals\n")
cat("      Combined xanthids:", sum(analysis_data_species$n_xanthid_harmful > 0), "corals,",
    sum(analysis_data_species$n_xanthid_harmful), "individuals\n")
cat("      Alpheus diadema:", sum(analysis_data_species$n_alpheus_diadema > 0), "corals,",
    sum(analysis_data_species$n_alpheus_diadema), "individuals\n\n")

# F.2 Run models for each key species
cat("F.2 Testing Key Species Effects on Coral Condition...\n")
cat("    Model: Condition ~ Key_Species + log(Volume) + Site (fixed effect)\n\n")

# Define key species with expected directions from experimental paper
key_species_predictors <- list(
  list(name = "Caracanthus maculatus", col = "n_caracanthus", expected = "positive",
       note = "Velvetfish - experimental positive effect"),
  list(name = "Alpheus lottini", col = "n_alpheus_lottini", expected = "positive",
       note = "Snapping shrimp - experimental positive effect"),
  list(name = "All Alpheus spp.", col = "n_alpheus_all", expected = "positive",
       note = "All snapping shrimp species"),
  list(name = "Harpiliopsis spinigera", col = "n_harpiliopsis", expected = "positive",
       note = "Shrimp - experimental positive trend (beta=0.456)"),
  list(name = "Periclimenes watamuae", col = "n_periclimenes", expected = "positive",
       note = "Shrimp - experimental positive trend (beta=0.213)"),
  list(name = "Calcinus latens", col = "n_calcinus", expected = "positive",
       note = "Hermit crab - experimental positive trend (beta=0.454)"),
  list(name = "Cymo quadrilobatus", col = "n_cymo", expected = "negative",
       note = "Xanthid crab - experimental negative effect"),
  list(name = "Luniella pugil", col = "n_luniella", expected = "negative",
       note = "Xanthid crab - experimental negative effect"),
  list(name = "Harmful xanthids (combined)", col = "n_xanthid_harmful", expected = "negative",
       note = "Cymo + Luniella combined"),
  list(name = "Alpheus diadema", col = "n_alpheus_diadema", expected = "negative",
       note = "Snapping shrimp - experimental negative effect (beta=-0.995)")
)

key_species_results <- list()

for (sp in key_species_predictors) {
  # Check if there's enough variation to fit model
  n_present <- sum(analysis_data_species[[sp$col]] > 0)

  if (n_present >= 5) {  # Need at least 5 corals with species present
    result <- run_condition_model(analysis_data_species, sp$name, sp$col, transform = "sqrt")

    if (!is.null(result)) {
      # Add expected direction and match status
      result$result$expected_sign <- sp$expected
      result$result$observed_sign <- ifelse(result$result$estimate > 0, "positive", "negative")
      result$result$matches_experiment <- result$result$expected_sign == result$result$observed_sign
      result$result$note <- sp$note
      result$result$n_present <- n_present

      key_species_results[[sp$name]] <- result$result

      # Print result with hypothesis test
      r <- result$result
      sig_star <- ifelse(r$p_value < 0.001, "***",
                         ifelse(r$p_value < 0.01, "**",
                                ifelse(r$p_value < 0.05, "*",
                                       ifelse(r$p_value < 0.1, ".", ""))))
      match_text <- ifelse(r$matches_experiment, "MATCHES", "CONTRADICTS")

      cat("   ", sp$name, "(expected:", toupper(sp$expected), ")\n", sep = "")
      cat("      Present on:", n_present, "corals\n")
      cat("      beta =", round(r$estimate, 4),
          ", SE =", round(r$se, 4),
          ", t =", round(r$t_value, 2),
          ", p =", format.pval(r$p_value, 3), sig_star, "\n")
      cat("      Direction:", toupper(r$observed_sign), "-", match_text, "experimental prediction\n\n")
    }
  } else {
    cat("   ", sp$name, ": SKIPPED (only present on", n_present, "corals)\n\n", sep = "")
  }
}

# Combine key species results
key_species_df <- bind_rows(key_species_results)

# Apply FDR correction across key species tests (up to 10 tests)
if (nrow(key_species_df) > 0) {
  key_species_df <- key_species_df %>%
    mutate(
      p_adj = p.adjust(p_value, method = "hochberg"),  # OLS primary
      p_adj_hc3 = p.adjust(p_value_robust, method = "hochberg"),  # HC3 supplement
      significant_adj = p_adj < 0.05,
      significant = p_value < 0.05  # OLS primary
    )

  cat("\nHochberg-corrected key species p-values (Hochberg (FWER)):\n")
  for (i in 1:nrow(key_species_df)) {
    row <- key_species_df[i, ]
    adj_star <- ifelse(row$significant_adj, " *ADJ-SIG*", "")
    cat(sprintf("   %-28s p = %.4f, p_adj = %.4f%s\n",
                row$predictor, row$p_value, row$p_adj, adj_star))
  }
  cat("\n")
}

# F.3 Summary comparison
cat("F.3 Key Species Effects Summary:\n")
cat("    =========================================================================\n")
cat("    Species                      | Expected | Observed | Match | p-value\n")
cat("    -------------------------------------------------------------------------\n")

for (i in 1:nrow(key_species_df)) {
  r <- key_species_df[i, ]
  match_status <- ifelse(r$matches_experiment, "YES", "NO")
  sig_marker <- ifelse(r$p_value < 0.05, "*", "")
  cat(sprintf("    %-28s | %-8s | %-8s | %-5s | %.4f%s\n",
              r$predictor, toupper(r$expected_sign), toupper(r$observed_sign),
              match_status, r$p_value, sig_marker))
}
cat("    =========================================================================\n")
cat("    * = significant at p < 0.05\n\n")

# F.4 Add to stats results
for (i in 1:nrow(key_species_df)) {
  row <- key_species_df[i, ]
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H-key-species",
      question = paste("Does", row$predictor, "predict coral condition as in experiment?"),
      test_name = "Linear model (fixed-effect site)",
      test_statistic = row$t_value,
      df = as.character(round(row$df, 1)),
      p_value = row$p_value,
      effect_size = row$estimate,
      effect_type = "Regression coefficient",
      n = row$n,
      interpretation = ifelse(row$matches_experiment,
                              paste("MATCHES experimental prediction (", row$expected_sign, ")", sep = ""),
                              paste("CONTRADICTS experimental prediction (expected", row$expected_sign, ")", sep = " ")),
      notes = row$note
    )
  )
}

# F.5 Create key species effects visualization
cat("F.5 Creating Key Species Effects Figure...\n")

if (nrow(key_species_df) > 0) {
  # Prepare data for forest plot
  key_species_plot_data <- key_species_df %>%
    mutate(
      species_label = factor(predictor, levels = rev(predictor)),
      effect_category = ifelse(expected_sign == "positive", "Expected Positive", "Expected Negative"),
      color_var = case_when(
        significant & matches_experiment ~ "Sig. matches experiment",
        significant & !matches_experiment ~ "Sig. contradicts experiment",
        !significant & matches_experiment ~ "NS matches experiment",
        TRUE ~ "NS contradicts experiment"
      )
    )

  # Colors for hypothesis matching (colorblind-safe: blue/orange)
  hypothesis_colors <- c(
    "Sig. matches experiment" = "#0072B2",      # Dark blue
    "Sig. contradicts experiment" = "#D55E00",  # Vermillion
    "NS matches experiment" = "#56B4E9",        # Sky blue
    "NS contradicts experiment" = "#E69F00"     # Orange
  )

  p_key_species_forest <- ggplot(key_species_plot_data,
                                  aes(x = species_label, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
    geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, fill = color_var),
                    color = "black", shape = 21, size = 1.2, stroke = 0.4, linewidth = 1) +
    geom_text(aes(label = sprintf("p = %.3f\n(n = %d)", p_value, n_present)),
              hjust = -0.2, size = 2.8, lineheight = 0.9) +
    scale_fill_manual(values = hypothesis_colors, name = "Hypothesis Test") +
    coord_flip() +
    facet_wrap(~ effect_category, scales = "free_y", ncol = 1) +
    labs(
      title = "Key Species Effects on Coral Condition",
      subtitle = "Testing experimental predictions from Stier et al. in survey data",
      x = "",
      y = "Effect on condition score (regression coefficient)",
      caption = paste0("n = ", nrow(analysis_data_species), " corals | Vertical line at 0 = no effect\n",
                       "Blue = matches experimental prediction | Orange = contradicts")
    ) +
    theme_publication() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 11),
      plot.caption = element_text(hjust = 0, size = 9)
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.4)))  # Extra space for labels

  save_figure(p_key_species_forest,
             file.path(fig_dir, "key_species_effects.png"),
             width = 10, height = 10)
  cat("   Saved: key_species_effects.png\n\n")

  # F.6 Combined comparison: Functional Groups vs Key Species
  cat("F.6 Creating Combined Effects Comparison Figure...\n")

  # Prepare functional group data with same structure
  func_for_comparison <- functional_effects %>%
    mutate(
      category = "Functional Group",
      predictor_clean = case_when(
        predictor == "Trapezia abundance" ~ "Trapezia (defenders)",
        predictor == "Resident Fish abundance" ~ "Resident Fish",
        predictor == "Galeropsis abundance" ~ "Galeropsis (tissue consumer)",
        TRUE ~ predictor
      ),
      hypothesis_match = matches_hypothesis
    ) %>%
    dplyr::select(predictor = predictor_clean, estimate, se, ci_lower, ci_upper,
           p_value, significant, expected_sign, hypothesis_match, category)

  # Prepare key species data with same structure
  key_for_comparison <- key_species_df %>%
    mutate(
      category = "Key Species",
      predictor_clean = predictor,
      hypothesis_match = matches_experiment
    ) %>%
    dplyr::select(predictor = predictor_clean, estimate, se, ci_lower, ci_upper,
           p_value, significant, expected_sign, hypothesis_match, category)

  # Combine
  combined_comparison <- bind_rows(func_for_comparison, key_for_comparison) %>%
    mutate(
      color_var = case_when(
        significant & hypothesis_match ~ "Sig. matches prediction",
        significant & !hypothesis_match ~ "Sig. contradicts prediction",
        !significant & hypothesis_match ~ "NS matches prediction",
        TRUE ~ "NS contradicts prediction"
      )
    )

  combined_colors <- c(
    "Sig. matches prediction" = "#0072B2",
    "Sig. contradicts prediction" = "#D55E00",
    "NS matches prediction" = "#56B4E9",
    "NS contradicts prediction" = "#E69F00"
  )

  p_combined <- ggplot(combined_comparison,
                        aes(x = reorder(predictor, estimate), y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
    geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, fill = color_var),
                    color = "black", shape = 21, size = 1, stroke = 0.4, linewidth = 0.8) +
    geom_text(aes(label = sprintf("p = %.3f", p_value)),
              hjust = -0.1, size = 2.8) +
    scale_fill_manual(values = combined_colors, name = "Hypothesis\nTest") +
    coord_flip() +
    facet_wrap(~ category, scales = "free_y", ncol = 1) +
    labs(
      title = "Effects on Coral Condition: Functional Groups vs Key Species",
      subtitle = "Comparing a priori hypotheses from functional ecology and experimental paper",
      x = "",
      y = "Effect size (regression coefficient)",
      caption = "Models: Condition ~ Predictor + log(Volume) + Site"
    ) +
    theme_publication() +
    theme(legend.position = "right",
          strip.text = element_text(face = "bold")) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.25)))

  save_figure(p_combined,
             file.path(fig_dir, "functional_vs_key_species.png"),
             width = 11, height = 9)
  cat("   Saved: functional_vs_key_species.png\n\n")
}

# F.7 Save key species results table
save_table(key_species_df %>%
             mutate(across(where(is.numeric), ~round(., 4))) %>%
             dplyr::select(predictor, estimate, se, t_value, p_value, p_value_robust, p_adj,
                    ci_lower, ci_upper,
                    expected_sign, observed_sign, matches_experiment, n_present, note),
           "key_species_effects")
cat("   Saved: key_species_effects.csv\n\n")

# F.8 Interpretation
cat("F.8 Key Species Analysis Interpretation:\n")
cat("    =====================================================================\n")

n_match <- sum(key_species_df$matches_experiment)
n_total <- nrow(key_species_df)
n_sig_match <- sum(key_species_df$matches_experiment & key_species_df$significant)
n_sig <- sum(key_species_df$significant)

cat("    Species tested:", n_total, "\n")
cat("    Direction matches experiment:", n_match, "/", n_total,
    "(", round(n_match/n_total*100, 0), "%)\n")
cat("    Significant effects:", n_sig, "/", n_total, "\n")
cat("    Significant AND matches:", n_sig_match, "/", n_sig,
    ifelse(n_sig > 0, paste0(" (", round(n_sig_match/n_sig*100, 0), "%)"), ""), "\n")

if (n_sig_match > 0) {
  cat("\n    SUPPORTED experimental predictions:\n")
  supported <- key_species_df %>% filter(significant & matches_experiment)
  for (i in 1:nrow(supported)) {
    cat("      -", supported$predictor[i], ":", toupper(supported$observed_sign[i]),
        "effect (p =", round(supported$p_value[i], 4), ")\n")
  }
}

if (n_sig > n_sig_match) {
  cat("\n    CONTRADICTED experimental predictions:\n")
  contradicted <- key_species_df %>% filter(significant & !matches_experiment)
  if (nrow(contradicted) > 0) {
    for (i in 1:nrow(contradicted)) {
      cat("      -", contradicted$predictor[i], ": expected", toupper(contradicted$expected_sign[i]),
          "but observed", toupper(contradicted$observed_sign[i]),
          "(p =", round(contradicted$p_value[i], 4), ")\n")
    }
  }
}

cat("    =====================================================================\n\n")

# ############################################################################
# PART F2: SPECIES x INDIVIDUAL TRAIT CORRELATIONS (Cross-study comparison)
# ############################################################################
# Analogous to Table S3 in the experimental paper (Stier et al.):
# Pearson correlations between each top CAFI species' abundance and
# individual coral condition metrics (protein, carbohydrate, AFDW,
# zooxanthellae density). This reveals whether the composite PC1_Coral
# masks species-specific effects on individual traits.
# ############################################################################

cat("============================================================\n")
cat("PART F2: SPECIES x INDIVIDUAL TRAIT CORRELATIONS\n")
cat("============================================================\n\n")

cat("Computing Pearson correlations between top CAFI species abundances\n")
cat("and individual coral condition metrics (cf. experimental Table S3).\n\n")

# Get the top 20 most abundant species in the condition dataset
top_species_for_corr <- cafi_clean %>%
  filter(coral_id %in% analysis_data$coral_id) %>%
  group_by(otu) %>%
  summarise(
    total_abundance = n(),
    n_corals = n_distinct(coral_id),
    .groups = "drop"
  ) %>%
  filter(n_corals >= 5) %>%  # present on at least 5 corals

  arrange(desc(total_abundance)) %>%
  slice_head(n = 20)

cat("   Top", nrow(top_species_for_corr), "species (present on >= 5 corals):\n")

# Build per-coral abundance matrix for these species
species_abundance_wide <- cafi_clean %>%
  filter(coral_id %in% analysis_data$coral_id,
         otu %in% top_species_for_corr$otu) %>%
  group_by(coral_id, otu) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = otu, values_from = n, values_fill = 0)

# Merge with condition data
trait_corr_data <- analysis_data %>%
  dplyr::select(coral_id, condition_score,
                any_of(c("protein_corr", "carb_corr", "zoox_corr", "afdw_corr"))) %>%
  left_join(species_abundance_wide, by = "coral_id")

# Replace NA species counts with 0
species_names <- top_species_for_corr$otu
for (sp in species_names) {
  if (sp %in% names(trait_corr_data)) {
    trait_corr_data[[sp]] <- ifelse(is.na(trait_corr_data[[sp]]), 0, trait_corr_data[[sp]])
  }
}

# Compute correlations
trait_cols <- c("condition_score", "protein_corr", "carb_corr", "zoox_corr", "afdw_corr")
trait_labels <- c("PC1_Coral", "Protein", "Carbohydrate", "Zooxanthellae", "AFDW")
available_traits <- trait_cols[trait_cols %in% names(trait_corr_data)]
available_labels <- trait_labels[trait_cols %in% names(trait_corr_data)]

species_trait_corr <- data.frame()

for (sp in species_names) {
  if (!sp %in% names(trait_corr_data)) next

  # Assign taxonomic group based on cafi_clean type column
  sp_type <- cafi_clean %>%
    filter(otu == sp) %>%
    pull(type) %>%
    unique() %>%
    .[1]
  taxon_group <- case_when(
    sp_type == "fish" ~ "Fishes",
    sp_type %in% c("crab", "shrimp") ~ "Shrimps/Crabs",
    sp_type == "snail" ~ "Snails",
    TRUE ~ "Other"
  )

  row <- data.frame(species = sp, taxon_group = taxon_group, stringsAsFactors = FALSE)

  for (j in seq_along(available_traits)) {
    trait <- available_traits[j]
    label <- available_labels[j]
    complete_cases <- complete.cases(trait_corr_data[[sp]], trait_corr_data[[trait]])
    # Spearman rank correlation (robust to zero-inflation in raw counts). See Fig S16 for regression-based analysis.
    if (sum(complete_cases) >= 10) {
      ct <- cor.test(trait_corr_data[[sp]][complete_cases],
                     trait_corr_data[[trait]][complete_cases],
                     method = "spearman", exact = FALSE)
      row[[paste0("r_", label)]] <- round(ct$estimate, 2)
      row[[paste0("p_", label)]] <- ct$p.value
    } else {
      row[[paste0("r_", label)]] <- NA
      row[[paste0("p_", label)]] <- NA
    }
  }

  species_trait_corr <- bind_rows(species_trait_corr, row)
}

# Order by taxon group then by PC1_Coral correlation within group
if (nrow(species_trait_corr) > 0 && "r_PC1_Coral" %in% names(species_trait_corr)) {
  species_trait_corr <- species_trait_corr %>%
    arrange(taxon_group, desc(r_PC1_Coral))
}

# Apply BH-FDR correction across all p-values
p_cols <- grep("^p_", names(species_trait_corr), value = TRUE)
if (length(p_cols) > 0) {
  all_pvals <- unlist(species_trait_corr[p_cols])
  all_pvals_adj <- p.adjust(all_pvals, method = "BH")
  # Write back adjusted p-values with _fdr suffix
  fdr_matrix <- matrix(all_pvals_adj, nrow = nrow(species_trait_corr), ncol = length(p_cols))
  colnames(fdr_matrix) <- gsub("^p_", "p_fdr_", p_cols)
  species_trait_corr <- bind_cols(species_trait_corr, as.data.frame(fdr_matrix))
}

# Save table
save_table(species_trait_corr, "species_trait_correlations")
cat("   Saved: species_trait_correlations.csv\n")
cat("   Species included:", nrow(species_trait_corr), "\n")
cat("   Traits tested:", paste(available_labels, collapse = ", "), "\n")
cat("   Correlation method: Spearman rank (BH-FDR corrected)\n\n")

# Print summary
cat("   Species with |r| > 0.2 for PC1_Coral:\n")
if ("r_PC1_Coral" %in% names(species_trait_corr)) {
  strong <- species_trait_corr %>% filter(abs(r_PC1_Coral) > 0.2)
  if (nrow(strong) > 0) {
    for (i in 1:nrow(strong)) {
      cat("      ", strong$species[i], "(", strong$taxon_group[i], "): r =",
          strong$r_PC1_Coral[i], "\n")
    }
  } else {
    cat("      None\n")
  }
}
cat("\n")

# ############################################################################
# PART F3: INDIVIDUAL PHYSIOLOGY RESPONSES TO CAFI PREDICTORS
# ############################################################################
# Tests whether individual condition metrics respond to CAFI predictors
# (total abundance, richness, PC1_CAFI). Analogous to experimental Table S1
# but testing CAFI predictors instead of treatment effects.
# ############################################################################

cat("============================================================\n")
cat("PART F3: INDIVIDUAL PHYSIOLOGY RESPONSES TO CAFI PREDICTORS\n")
cat("============================================================\n\n")

cat("Testing individual physiology metrics against CAFI predictors.\n")
cat("This demonstrates whether the PC1_Coral null result holds at\n")
cat("the individual trait level.\n\n")

physio_traits <- c("protein_corr", "carb_corr", "zoox_corr", "afdw_corr")
physio_names <- c("Protein", "Carbohydrate", "Zooxanthellae", "AFDW")
physio_cafi_predictors <- c("total_cafi", "otu_richness", "pc1_cafi")
cafi_pred_names <- c("Total CAFI", "Species Richness", "PC1_CAFI")

individual_physio_results <- data.frame()

for (i in seq_along(physio_traits)) {
  trait <- physio_traits[i]
  trait_name <- physio_names[i]

  if (!trait %in% names(analysis_data)) {
    cat("   ", trait_name, ": variable not available\n")
    next
  }

  for (j in seq_along(physio_cafi_predictors)) {
    pred <- physio_cafi_predictors[j]
    pred_name <- cafi_pred_names[j]

    if (!pred %in% names(analysis_data)) next

    data_complete <- analysis_data %>%
      filter(!is.na(.data[[trait]]), !is.na(.data[[pred]]))

    if (nrow(data_complete) < 20) next

    formula_str <- paste(trait, "~", pred, "+ log_volume + site")
    m <- tryCatch(lm(as.formula(formula_str), data = data_complete), error = function(e) NULL)

    if (is.null(m)) next

    coefs <- summary(m)$coefficients
    if (!pred %in% rownames(coefs)) next

    # Use robust SEs if sandwich is available
    robust_se <- tryCatch({
      vcov_hc <- sandwich::vcovHC(m, type = "HC3")
      ct <- lmtest::coeftest(m, vcov. = vcov_hc)
      list(se = ct[pred, "Std. Error"], p = ct[pred, "Pr(>|t|)"])
    }, error = function(e) {
      cat("    Note: Using OLS SEs for", pred, "(HC3 failed)\n")
      list(se = coefs[pred, "Std. Error"], p = coefs[pred, "Pr(>|t|)"])
    })

    individual_physio_results <- bind_rows(individual_physio_results, data.frame(
      trait = trait_name,
      predictor = pred_name,
      estimate = round(coefs[pred, "Estimate"], 4),
      se = round(robust_se$se, 4),
      t_value = round(coefs[pred, "Estimate"] / robust_se$se, 2),
      p_value = round(robust_se$p, 4),
      n = nrow(data_complete),
      stringsAsFactors = FALSE
    ))
  }
}

if (nrow(individual_physio_results) > 0) {
  # FDR correction across all tests
  individual_physio_results$p_fdr <- round(p.adjust(individual_physio_results$p_value, method = "BH"), 4)

  save_table(individual_physio_results, "individual_physiology_cafi_responses")
  cat("   Saved: individual_physiology_cafi_responses.csv\n")
  cat("   Tests run:", nrow(individual_physio_results), "\n")
  cat("   Significant (p < 0.05):", sum(individual_physio_results$p_value < 0.05), "\n")
  cat("   Significant after FDR:", sum(individual_physio_results$p_fdr < 0.05), "\n\n")

  # Print summary table
  cat("   Individual Physiology ~ CAFI Predictor Results:\n")
  cat("   ", sprintf("%-15s %-18s %8s %8s %8s\n", "Trait", "Predictor", "beta", "p", "p_FDR"))
  cat("   ", paste(rep("-", 60), collapse = ""), "\n")
  for (k in 1:nrow(individual_physio_results)) {
    r <- individual_physio_results[k, ]
    sig <- ifelse(r$p_fdr < 0.05, " *", "")
    cat("   ", sprintf("%-15s %-18s %8.4f %8.4f %8.4f%s\n",
                        r$trait, r$predictor, r$estimate, r$p_value, r$p_fdr, sig))
  }
  cat("\n")
} else {
  cat("   No individual physiology models could be fit.\n\n")
}

# ############################################################################
# PART F4: CROSS-STUDY SPECIES COMPARISON TABLE
# ############################################################################
# Creates a table comparing species-level results between this survey and
# the companion experimental paper (Stier, Primo, Curtis, Osenberg) for
# overlapping species. Experimental values hardcoded from Table S2.
# ############################################################################

cat("============================================================\n")
cat("PART F4: CROSS-STUDY SPECIES COMPARISON TABLE\n")
cat("============================================================\n\n")

# Source: Stier, Primo, Curtis & Osenberg (companion experiment), Table S2.
# Values last verified 2026-01-30. Update if companion paper revises.
expt_species <- tribble(
  ~species,                  ~taxon_group,     ~expt_condition_beta, ~expt_condition_p, ~expt_condition_p_adj,
  "Caracanthus maculatus",   "Fishes",          1.424,  0.001,  0.017,
  "Dascyllus flavicaudus",   "Fishes",          0.321,  0.125,  0.963,
  "Dascyllus aruanus",       "Fishes",          NA,     NA,     NA,     # Not in Table S2 individually
  "Harpiliopsis spinigera",  "Shrimps/Crabs",   0.456,  0.146,  0.963,
  "Periclimenes watamuae",   "Shrimps/Crabs",   0.213,  0.108,  0.963,
  "Calcinus latens",         "Shrimps/Crabs",   0.454,  0.078,  0.963,
  "Trapezia serenei",        "Shrimps/Crabs",  -0.022,  0.963,  0.963,
  "Alpheus lottini",         "Shrimps/Crabs",   NA,     NA,     NA,     # Not in Table S2
  "Alpheus diadema",         "Shrimps/Crabs",  -0.995,  0.042,  0.759,
  "Luniella pugil",          "Shrimps/Crabs",  -0.917,  0.014,  0.257,
  "Galeropsis monodonta",    "Snails",         -0.099,  0.712,  0.963
)

# Match with survey key species results
cross_study <- expt_species

if (exists("key_species_df") && nrow(key_species_df) > 0) {
  # Map survey species names to experimental species names
  survey_map <- c(
    "Caracanthus maculatus" = "Caracanthus maculatus",
    "Harpiliopsis spinigera" = "Harpiliopsis spinigera",
    "Periclimenes watamuae" = "Periclimenes watamuae",
    "Calcinus latens" = "Calcinus latens",
    "Alpheus lottini" = "Alpheus lottini",
    "Alpheus diadema" = "Alpheus diadema",
    "Luniella pugil" = "Luniella pugil"
  )

  cross_study$survey_condition_beta <- NA_real_
  cross_study$survey_condition_p <- NA_real_
  cross_study$survey_condition_p_adj <- NA_real_
  cross_study$survey_n_corals <- NA_integer_
  cross_study$direction_match <- NA

  for (i in 1:nrow(cross_study)) {
    sp <- cross_study$species[i]
    # Find matching survey result
    survey_match <- key_species_df %>%
      filter(predictor == sp)
    if (nrow(survey_match) == 0) next

    cross_study$survey_condition_beta[i] <- survey_match$estimate[1]
    cross_study$survey_condition_p[i] <- survey_match$p_value_robust[1]
    cross_study$survey_condition_p_adj[i] <- survey_match$p_adj[1]
    cross_study$survey_n_corals[i] <- survey_match$n_present[1]

    # Check if directions match
    if (!is.na(cross_study$expt_condition_beta[i])) {
      cross_study$direction_match[i] <- sign(cross_study$expt_condition_beta[i]) ==
                                         sign(cross_study$survey_condition_beta[i])
    }
  }
}

# Save
save_table(cross_study %>%
             mutate(across(where(is.numeric), ~round(., 4))),
           "cross_study_species_comparison")
cat("   Saved: cross_study_species_comparison.csv\n")
cat("   Species compared:", nrow(cross_study), "\n")
if ("direction_match" %in% names(cross_study)) {
  n_matched <- sum(cross_study$direction_match, na.rm = TRUE)
  n_testable <- sum(!is.na(cross_study$direction_match))
  cat("   Direction matches:", n_matched, "/", n_testable, "\n")

  # Sign-concordance test: binomial test of whether survey effect directions

  # match experimental directions more than expected by chance (50%)
  if (n_testable >= 3) {
    sign_test <- binom.test(n_matched, n_testable, p = 0.5, alternative = "greater")
    cat("\n   SIGN CONCORDANCE TEST (binomial, H0: 50% match by chance):\n")
    cat("     Concordant:", n_matched, "/", n_testable,
        sprintf("(%.0f%%)\n", 100 * n_matched / n_testable))
    cat("     p =", format.pval(sign_test$p.value, 3), "\n")
    cat("     95% CI for concordance rate: [",
        sprintf("%.1f%%", 100 * sign_test$conf.int[1]), ",",
        sprintf("%.1f%%", 100 * sign_test$conf.int[2]), "]\n")
    if (sign_test$p.value < 0.05) {
      cat("     → Effect directions are significantly concordant across studies\n")
    } else {
      cat("     → No evidence of systematic concordance (directions may differ)\n")
    }

    # Save sign concordance result
    sign_concordance_df <- tibble(
      n_concordant = n_matched,
      n_testable = n_testable,
      concordance_rate = n_matched / n_testable,
      binom_p = sign_test$p.value,
      ci_lower = sign_test$conf.int[1],
      ci_upper = sign_test$conf.int[2]
    )
    save_table(sign_concordance_df, "cross_study_sign_concordance")
    cat("     Saved: cross_study_sign_concordance.csv\n")
  }
}
cat("\n")

# ############################################################################
# PART G: NEIGHBORHOOD EFFECTS (Analog to Experimental Treatment)
# ############################################################################
# The experimental paper (Stier et al.) manipulated coral number to test
# how neighborhood density affects CAFI abundance and coral condition.
# Here we test whether natural variation in neighborhood density shows
# similar patterns using our survey data.
#
# Key variables:
#   - n_neighbors: number of neighboring Pocillopora within 5m radius
#   - total_neighbor_volume: combined volume of all neighbors
# ############################################################################

cat("============================================================\n")
cat("PART G: NEIGHBORHOOD EFFECTS (Analog to Experimental Treatment)\n")
cat("============================================================\n\n")

cat("Testing whether neighborhood density (analogous to experimental manipulation\n")
cat("of coral number) affects CAFI abundance and coral condition.\n\n")

# G.1 Prepare data with neighborhood variables
# Only subset of corals have neighborhood data (from 5m surveys)
neighborhood_data <- coral_master %>%
  filter(!is.na(n_neighbors), !is.na(volume)) %>%
  mutate(
    log_volume = log(volume),
    log_volume_scaled = as.numeric(scale(log_volume)),
    # Size category for visualization
    size_category = cut(volume,
                        breaks = quantile(volume, probs = c(0, 0.33, 0.67, 1)),
                        labels = c("Small", "Medium", "Large"),
                        include.lowest = TRUE)
  )

# Pre-compute scaled neighborhood predictors outside mutate to avoid if()/scale() mismatch
# Guard against constant values (all identical -> sd=0 -> NaN from scale())
g_nn <- neighborhood_data$n_neighbors
g_nn_scaled <- if (sd(g_nn, na.rm = TRUE) > 0) as.numeric(scale(g_nn)) else rep(0, length(g_nn))

g_tnv <- log(neighborhood_data$total_neighbor_volume + 1)
g_tnv_scaled <- if (sd(g_tnv, na.rm = TRUE) > 0) as.numeric(scale(g_tnv)) else rep(0, length(g_tnv))

neighborhood_data <- neighborhood_data %>%
  mutate(
    n_neighbors_scaled = g_nn_scaled,
    total_neighbor_vol_scaled = g_tnv_scaled
  )

n_neighborhood <- nrow(neighborhood_data)
n_with_condition <- sum(!is.na(neighborhood_data$condition_score))

cat("G.1 Sample sizes:\n")
cat("    Corals with neighborhood data:", n_neighborhood, "\n")
cat("    With BOTH neighborhood AND condition:", n_with_condition, "\n\n")

cat("G.2 Neighborhood variable summary:\n")
cat("    n_neighbors: range", min(neighborhood_data$n_neighbors), "-", max(neighborhood_data$n_neighbors),
    "| median =", median(neighborhood_data$n_neighbors), "\n")
cat("    total_neighbor_volume: range", round(min(neighborhood_data$total_neighbor_volume)),
    "-", round(max(neighborhood_data$total_neighbor_volume)),
    "| median =", round(median(neighborhood_data$total_neighbor_volume)), "cm³\n\n")

# G.3 Model 1: Neighborhood effects on CAFI abundance
cat("G.3 Model 1: CAFI ~ Volume × Neighborhood Density\n")
cat("    --------------------------------------------------\n")

m_cafi_neighborhood <- tryCatch(
  MASS::glm.nb(
    total_cafi ~ log_volume_scaled * n_neighbors_scaled + site,
    data = neighborhood_data
  ),
  error = function(e) {
    cat("    NB model failed:", conditionMessage(e), "\n")
    cat("    Falling back to Poisson\n")
    glm(total_cafi ~ log_volume_scaled * n_neighbors_scaled + site,
        family = poisson, data = neighborhood_data)
  }
)

summary_m1 <- summary(m_cafi_neighborhood)
cat("    log_volume (scaled):         β =", round(coef(m_cafi_neighborhood)["log_volume_scaled"], 3),
    ", p =", round(summary_m1$coefficients["log_volume_scaled", "Pr(>|z|)"], 4), "\n")
cat("    n_neighbors (scaled):        β =", round(coef(m_cafi_neighborhood)["n_neighbors_scaled"], 3),
    ", p =", round(summary_m1$coefficients["n_neighbors_scaled", "Pr(>|z|)"], 4), "\n")
cat("    Volume × Neighborhood:       β =", round(coef(m_cafi_neighborhood)["log_volume_scaled:n_neighbors_scaled"], 3),
    ", p =", round(summary_m1$coefficients["log_volume_scaled:n_neighbors_scaled", "Pr(>|z|)"], 4), "\n\n")

# G.4 Model 2: Neighborhood effects on coral condition
cat("G.4 Model 2: Condition ~ Volume × Neighborhood Density\n")
cat("    --------------------------------------------------\n")

neighborhood_condition <- neighborhood_data %>%
  filter(!is.na(condition_score))

if (nrow(neighborhood_condition) >= 20) {
  # Fixed-effect site (3 levels insufficient for random intercepts)
  m_condition_neighborhood <- lm(
    condition_score ~ log_volume_scaled * n_neighbors_scaled + site,
    data = neighborhood_condition
  )

  summary_m2 <- summary(m_condition_neighborhood)
  coef_m2 <- as.data.frame(summary_m2$coefficients)

  cat("    log_volume (scaled):         β =", round(coef_m2["log_volume_scaled", "Estimate"], 3),
      ", t =", round(coef_m2["log_volume_scaled", "t value"], 2),
      ", p =", round(coef_m2["log_volume_scaled", "Pr(>|t|)"], 4), "\n")
  cat("    n_neighbors (scaled):        β =", round(coef_m2["n_neighbors_scaled", "Estimate"], 3),
      ", t =", round(coef_m2["n_neighbors_scaled", "t value"], 2),
      ", p =", round(coef_m2["n_neighbors_scaled", "Pr(>|t|)"], 4), "\n")
  cat("    Volume × Neighborhood:       β =", round(coef_m2["log_volume_scaled:n_neighbors_scaled", "Estimate"], 3),
      ", t =", round(coef_m2["log_volume_scaled:n_neighbors_scaled", "t value"], 2),
      ", p =", round(coef_m2["log_volume_scaled:n_neighbors_scaled", "Pr(>|t|)"], 4), "\n\n")

  # ---- Sensitivity power analysis (replaces invalid post-hoc power) ----
  # Post-hoc power based on observed data is invalid (Hoenig & Heisey 2001).
  # Instead, report the minimum detectable effect size at 80% power.
  cat("\n  Sensitivity power analysis (neighborhood effect on condition):\n")
  cat("  (What effect could this study detect at 80% power?)\n")

  n_cond <- nrow(neighborhood_condition)
  n_predictors_neigh <- 4  # log_volume, n_neighbors, interaction, 2 site dummies = 4 non-intercept

  if (requireNamespace("pwr", quietly = TRUE)) {
    # Minimum detectable f2 at alpha=0.05, power=0.80
    pwr_result_neigh <- pwr::pwr.f2.test(
      u = 1,  # numerator df (1 predictor being tested)
      v = n_cond - n_predictors_neigh - 2,  # denominator df
      f2 = NULL,
      sig.level = 0.05,
      power = 0.80
    )
    min_f2_neigh <- pwr_result_neigh$f2
    min_partial_r2_neigh <- min_f2_neigh / (1 + min_f2_neigh)

    cat("  Sample size for neighborhood condition models:", n_cond, "\n")
    cat("  Minimum detectable partial R2 at 80% power:", round(min_partial_r2_neigh, 3), "\n")
    cat("  Minimum detectable Cohen's f2:", round(min_f2_neigh, 3), "\n")
    cat("  This corresponds to a",
        ifelse(min_f2_neigh < 0.02, "very small", ifelse(min_f2_neigh < 0.15, "small-to-medium", "medium")),
        "effect.\n")
    cat("  Effects smaller than partial R2 =", round(min_partial_r2_neigh, 3),
        "cannot be reliably detected.\n\n")
  } else {
    cat("  (Install 'pwr' package for formal power analysis)\n")
    # Rough approximation: for n~61, minimum detectable r ~ 0.35
    cat("  Approximate minimum detectable correlation: r ~ 0.35 (medium effect)\n\n")
  }

} else {
  m_condition_neighborhood <- NULL
  cat("    Insufficient sample size for condition ~ neighborhood model\n\n")
}

# G.5 Model 3: Full model with Volume, Neighborhood, and CAFI
cat("G.5 Model 3: Condition ~ Volume + Neighborhood + CAFI (decomposition)\n")
cat("    --------------------------------------------------\n")

neighborhood_full <- neighborhood_data %>%
  filter(!is.na(condition_score), !is.na(pc1_cafi))

if (nrow(neighborhood_full) >= 20) {
  # Fixed-effect site (3 levels insufficient for random intercepts)
  m_full <- lm(
    condition_score ~ log_volume_scaled + n_neighbors_scaled + pc1_cafi_scaled + site,
    data = neighborhood_full %>%
      mutate(pc1_cafi_scaled = scale(pc1_cafi)[,1])
  )

  summary_m3 <- summary(m_full)
  coef_m3 <- as.data.frame(summary_m3$coefficients)

  cat("    log_volume (scaled):    β =", round(coef_m3["log_volume_scaled", "Estimate"], 3),
      ", p =", round(coef_m3["log_volume_scaled", "Pr(>|t|)"], 4), "\n")
  cat("    n_neighbors (scaled):   β =", round(coef_m3["n_neighbors_scaled", "Estimate"], 3),
      ", p =", round(coef_m3["n_neighbors_scaled", "Pr(>|t|)"], 4), "\n")
  cat("    PC1_CAFI (scaled):      β =", round(coef_m3["pc1_cafi_scaled", "Estimate"], 3),
      ", p =", round(coef_m3["pc1_cafi_scaled", "Pr(>|t|)"], 4), "\n")
  cat("    Adjusted R²:", round(summary_m3$adj.r.squared, 4), "\n")
} else {
  m_full <- NULL
  cat("    Insufficient sample size for full decomposition model\n")
}
cat("\n")

# G.6 Create 3-panel neighborhood effects figure
cat("G.6 Creating neighborhood effects figure...\n")

# Panel A: CAFI abundance vs neighborhood density (colored by coral size)
p_neighborhood_cafi <- ggplot(neighborhood_data,
                               aes(x = n_neighbors, y = total_cafi)) +
  geom_point(aes(color = size_category, size = volume), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Small" = "#56B4E9", "Medium" = "#E69F00", "Large" = "#0072B2"),
                     name = "Coral Size") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(
    title = "A. Neighborhood Density → CAFI Abundance",
    subtitle = paste0("Analog to experimental coral number manipulation (n = ", n_neighborhood, ")"),
    x = "Number of neighboring Pocillopora (within 5m)",
    y = "Total CAFI abundance",
    caption = "Point size = coral volume; Color = size tertile"
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

# Panel B: Condition vs neighborhood density
if (!is.null(m_condition_neighborhood)) {
  p_neighborhood_condition <- ggplot(neighborhood_condition,
                                      aes(x = n_neighbors, y = condition_score)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
    geom_point(aes(fill = site), shape = 21, color = "gray30", stroke = 0.4, alpha = 0.7, size = 3) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
    scale_fill_manual(values = SITE_COLORS, name = "Site") +
    labs(
      title = "B. Neighborhood Density → Coral Condition",
      subtitle = paste0("Position-corrected condition score (n = ", n_with_condition, ")"),
      x = "Number of neighboring Pocillopora (within 5m)",
      y = expression(Coral~condition~(PC1[Coral]))
    ) +
    theme_publication() +
    theme(legend.position = "bottom")
} else {
  p_neighborhood_condition <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient data", size = 6) +
    theme_void()
}

# Panel C: Effect size comparison
# Compile effect sizes from all three models
effect_comparison <- tibble(
  predictor = c("Coral Volume", "Neighborhood Density", "CAFI Community (PC1)"),
  estimate = NA_real_,
  ci_lower = NA_real_,
  ci_upper = NA_real_,
  type = c("Size", "Landscape", "CAFI")
)

# Get volume effect from CAFI model
if (!is.null(m_condition_neighborhood)) {
  coefs <- coef(summary(m_condition_neighborhood))
  se <- coefs[, "Std. Error"]

  effect_comparison$estimate[1] <- coefs["log_volume_scaled", "Estimate"]
  effect_comparison$ci_lower[1] <- coefs["log_volume_scaled", "Estimate"] - 1.96 * se["log_volume_scaled"]
  effect_comparison$ci_upper[1] <- coefs["log_volume_scaled", "Estimate"] + 1.96 * se["log_volume_scaled"]

  effect_comparison$estimate[2] <- coefs["n_neighbors_scaled", "Estimate"]
  effect_comparison$ci_lower[2] <- coefs["n_neighbors_scaled", "Estimate"] - 1.96 * se["n_neighbors_scaled"]
  effect_comparison$ci_upper[2] <- coefs["n_neighbors_scaled", "Estimate"] + 1.96 * se["n_neighbors_scaled"]
}

# Get CAFI effect from full model
if (!is.null(m_full)) {
  coefs_full <- coef(summary(m_full))
  se_full <- coefs_full[, "Std. Error"]

  effect_comparison$estimate[3] <- coefs_full["pc1_cafi_scaled", "Estimate"]
  effect_comparison$ci_lower[3] <- coefs_full["pc1_cafi_scaled", "Estimate"] - 1.96 * se_full["pc1_cafi_scaled"]
  effect_comparison$ci_upper[3] <- coefs_full["pc1_cafi_scaled", "Estimate"] + 1.96 * se_full["pc1_cafi_scaled"]
}

p_effect_comparison <- ggplot(effect_comparison %>% filter(!is.na(estimate)),
                               aes(x = reorder(predictor, estimate), y = estimate, fill = type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  geom_col(width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 0.8) +
  scale_fill_manual(values = c("Size" = "#5E35B1", "Landscape" = "#00897B", "CAFI" = "#F57C00"),
                    guide = "none") +
  coord_flip() +
  labs(
    title = "C. Relative Effect Sizes on Coral Condition",
    subtitle = "Comparing size, landscape, and CAFI community effects",
    x = "",
    y = "Standardized effect size (with 95% CI)",
    caption = "From model: Condition ~ Volume + Neighbors + PC1_CAFI + Site (fixed effect)"
  ) +
  theme_publication()

# Combine panels
fig_neighborhood <- (p_neighborhood_cafi | p_neighborhood_condition) /
  p_effect_comparison +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "Neighborhood Effects: Survey Analog to Experimental Coral Manipulation",
    subtitle = paste0("Testing whether natural variation in neighborhood density affects CAFI and coral condition"),
    caption = paste0("Experimental paper manipulated coral number (Stier et al.); ",
                     "here we use natural variation in n_neighbors within 5m radius\n",
                     "n = ", n_neighborhood, " corals with neighborhood data; ",
                     n_with_condition, " with condition scores")
  ) &
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11))

save_figure(fig_neighborhood,
           file.path(fig_dir, "neighborhood_effects.png"),
           width = 12, height = 10)
cat("   Saved: neighborhood_effects.png\n")

# G.7 Save neighborhood results table
cat("G.7 Saving neighborhood results table...\n")
neighborhood_results <- tibble(
  model = c("CAFI ~ Volume × Neighbors",
            "CAFI ~ Volume × Neighbors",
            "CAFI ~ Volume × Neighbors",
            "Condition ~ Volume × Neighbors",
            "Condition ~ Volume × Neighbors",
            "Condition ~ Volume × Neighbors",
            "Condition ~ Vol + Neighbors + CAFI",
            "Condition ~ Vol + Neighbors + CAFI",
            "Condition ~ Vol + Neighbors + CAFI"),
  predictor = c("log_volume", "n_neighbors", "interaction",
                "log_volume", "n_neighbors", "interaction",
                "log_volume", "n_neighbors", "pc1_cafi"),
  estimate = NA_real_,
  se = NA_real_,
  statistic = NA_real_,
  p_value = NA_real_
)

# Fill in CAFI model results
s1 <- summary(m_cafi_neighborhood)$coefficients
neighborhood_results$estimate[1] <- s1["log_volume_scaled", "Estimate"]
neighborhood_results$se[1] <- s1["log_volume_scaled", "Std. Error"]
neighborhood_results$statistic[1] <- s1["log_volume_scaled", "z value"]
neighborhood_results$p_value[1] <- s1["log_volume_scaled", "Pr(>|z|)"]

neighborhood_results$estimate[2] <- s1["n_neighbors_scaled", "Estimate"]
neighborhood_results$se[2] <- s1["n_neighbors_scaled", "Std. Error"]
neighborhood_results$statistic[2] <- s1["n_neighbors_scaled", "z value"]
neighborhood_results$p_value[2] <- s1["n_neighbors_scaled", "Pr(>|z|)"]

neighborhood_results$estimate[3] <- s1["log_volume_scaled:n_neighbors_scaled", "Estimate"]
neighborhood_results$se[3] <- s1["log_volume_scaled:n_neighbors_scaled", "Std. Error"]
neighborhood_results$statistic[3] <- s1["log_volume_scaled:n_neighbors_scaled", "z value"]
neighborhood_results$p_value[3] <- s1["log_volume_scaled:n_neighbors_scaled", "Pr(>|z|)"]

# Fill in condition model results
if (!is.null(m_condition_neighborhood)) {
  s2 <- coef(summary(m_condition_neighborhood))
  neighborhood_results$estimate[4] <- s2["log_volume_scaled", "Estimate"]
  neighborhood_results$se[4] <- s2["log_volume_scaled", "Std. Error"]
  neighborhood_results$statistic[4] <- s2["log_volume_scaled", "t value"]
  if ("Pr(>|t|)" %in% colnames(s2)) neighborhood_results$p_value[4] <- s2["log_volume_scaled", "Pr(>|t|)"]

  neighborhood_results$estimate[5] <- s2["n_neighbors_scaled", "Estimate"]
  neighborhood_results$se[5] <- s2["n_neighbors_scaled", "Std. Error"]
  neighborhood_results$statistic[5] <- s2["n_neighbors_scaled", "t value"]
  if ("Pr(>|t|)" %in% colnames(s2)) neighborhood_results$p_value[5] <- s2["n_neighbors_scaled", "Pr(>|t|)"]

  neighborhood_results$estimate[6] <- s2["log_volume_scaled:n_neighbors_scaled", "Estimate"]
  neighborhood_results$se[6] <- s2["log_volume_scaled:n_neighbors_scaled", "Std. Error"]
  neighborhood_results$statistic[6] <- s2["log_volume_scaled:n_neighbors_scaled", "t value"]
  if ("Pr(>|t|)" %in% colnames(s2)) neighborhood_results$p_value[6] <- s2["log_volume_scaled:n_neighbors_scaled", "Pr(>|t|)"]
}

# Fill in full model results
if (!is.null(m_full)) {
  s3 <- coef(summary(m_full))
  neighborhood_results$estimate[7] <- s3["log_volume_scaled", "Estimate"]
  neighborhood_results$se[7] <- s3["log_volume_scaled", "Std. Error"]
  neighborhood_results$statistic[7] <- s3["log_volume_scaled", "t value"]
  if ("Pr(>|t|)" %in% colnames(s3)) neighborhood_results$p_value[7] <- s3["log_volume_scaled", "Pr(>|t|)"]

  neighborhood_results$estimate[8] <- s3["n_neighbors_scaled", "Estimate"]
  neighborhood_results$se[8] <- s3["n_neighbors_scaled", "Std. Error"]
  neighborhood_results$statistic[8] <- s3["n_neighbors_scaled", "t value"]
  if ("Pr(>|t|)" %in% colnames(s3)) neighborhood_results$p_value[8] <- s3["n_neighbors_scaled", "Pr(>|t|)"]

  neighborhood_results$estimate[9] <- s3["pc1_cafi_scaled", "Estimate"]
  neighborhood_results$se[9] <- s3["pc1_cafi_scaled", "Std. Error"]
  neighborhood_results$statistic[9] <- s3["pc1_cafi_scaled", "t value"]
  if ("Pr(>|t|)" %in% colnames(s3)) neighborhood_results$p_value[9] <- s3["pc1_cafi_scaled", "Pr(>|t|)"]
}

save_table(neighborhood_results %>% mutate(across(where(is.numeric), ~round(., 4))),
           "neighborhood_effects")
cat("   Saved: neighborhood_effects.csv\n\n")

# G.8 Interpretation
cat("G.8 Neighborhood Effects Interpretation:\n")
cat("    =====================================================================\n")
cat("    Connection to experimental paper (Stier et al.):\n")
cat("    - Experiment: Manipulated coral NUMBER to test landscape effects\n")
cat("    - This survey: Uses natural variation in n_neighbors (within 5m)\n")
cat("    - Hypothesis: Higher neighborhood density → more CAFI (Field of Dreams)\n")
cat("                  OR → fewer CAFI per coral (Redirection)\n\n")

# Report key findings
vol_effect <- neighborhood_results$estimate[1]
neighbor_effect <- neighborhood_results$estimate[2]
interaction_effect <- neighborhood_results$estimate[3]

cat("    Key findings for CAFI abundance:\n")
cat("      - Volume effect: β =", round(vol_effect, 3),
    "(p =", round(neighborhood_results$p_value[1], 4), ")\n")
cat("      - Neighborhood effect: β =", round(neighbor_effect, 3),
    "(p =", round(neighborhood_results$p_value[2], 4), ")\n")
cat("      - Interaction: β =", round(interaction_effect, 3),
    "(p =", round(neighborhood_results$p_value[3], 4), ")\n")

neighbor_interp <- ifelse(neighbor_effect > 0,
                          "More neighbors → more CAFI (consistent with Field of Dreams at landscape level)",
                          "More neighbors → fewer CAFI per focal coral (consistent with Redirection)")
cat("\n    Interpretation:", neighbor_interp, "\n")
cat("    =====================================================================\n\n")

# ############################################################################
# PART H: LANDSCAPE-ONLY EFFECTS ON CORAL CONDITION
# ############################################################################
# Analyze how abiotic/landscape factors (size, neighborhood, site, depth)
# affect coral physiological condition WITHOUT including CAFI predictors.
#
# Response variables:
#   - Individual condition measures: protein_corr, carb_corr, zoox_corr, afdw_corr
#   - Composite condition: condition_score (PC1_Coral)
#
# Landscape predictors:
#   - log_volume: coral size (primary)
#   - n_neighbors: neighborhood density (within 5m)
#   - mean_neighbor_dist: proximity to neighbors
#   - total_neighbor_volume: neighborhood biomass
#   - site: reef location (fixed effect)
#   - depth_m: water depth
# ############################################################################

cat("============================================================\n")
cat("PART H: LANDSCAPE-ONLY EFFECTS ON CORAL CONDITION\n")
cat("============================================================\n\n")

cat("Analyzing how landscape factors (size, neighborhood, depth, site)\n")
cat("affect coral condition WITHOUT including CAFI predictors.\n")
cat("This isolates abiotic/spatial drivers of coral health.\n\n")

# H.1 Prepare landscape data
# Merge condition scores with all landscape predictors
landscape_condition <- condition_scores %>%
  dplyr::select(coral_id, site, condition_score,
                any_of(c("protein_corr", "carb_corr", "zoox_corr", "afdw_corr"))) %>%
  left_join(
    coral_master %>% dplyr::select(coral_id, volume, depth_m,
                                    n_neighbors, mean_neighbor_dist,
                                    total_neighbor_volume, mean_neighbor_volume),
    by = "coral_id"
  ) %>%
  filter(!is.na(condition_score), !is.na(volume)) %>%
  mutate(
    log_volume = log(volume),
    log_volume_scaled = as.numeric(scale(log_volume)),
    site = factor(site)
  )

# Pre-compute scaled neighborhood predictors before mutate to avoid if()/scale() mismatch
# Guard against constant values (all identical -> sd=0 -> NaN from scale())
nn_vals <- landscape_condition$n_neighbors
nn_scaled <- if (sd(nn_vals, na.rm = TRUE) > 0) as.numeric(scale(nn_vals)) else rep(0, length(nn_vals))

md_vals <- landscape_condition$mean_neighbor_dist
md_scaled <- if (sd(md_vals, na.rm = TRUE) > 0) as.numeric(scale(md_vals)) else rep(0, length(md_vals))

tnv_raw <- landscape_condition$total_neighbor_volume
tnv_log <- ifelse(is.na(tnv_raw), NA, log(tnv_raw + 1))
tnv_scaled <- if (sd(tnv_log, na.rm = TRUE) > 0) as.numeric(scale(tnv_log)) else rep(0, length(tnv_log))

landscape_condition <- landscape_condition %>%
  mutate(
    n_neighbors_scaled = nn_scaled,
    mean_dist_scaled = md_scaled,
    total_neighbor_vol_scaled = tnv_scaled
  )

n_landscape <- nrow(landscape_condition)
n_with_neighbors <- sum(!is.na(landscape_condition$n_neighbors))
n_with_depth <- sum(!is.na(landscape_condition$depth_m))

cat("H.1 Sample sizes:\n")
cat("    Corals with condition data:", n_landscape, "\n")
cat("    With neighborhood data:", n_with_neighbors, "\n")
cat("    With depth data:", n_with_depth, "\n\n")

# H.2 Individual condition measures vs landscape predictors
cat("H.2 INDIVIDUAL CONDITION MEASURES vs LANDSCAPE\n")
cat("    --------------------------------------------------\n")

individual_measures <- c("protein_corr", "carb_corr", "zoox_corr", "afdw_corr")
measure_names <- c("Protein", "Carbohydrate", "Zooxanthellae", "AFDW")

# Storage for results
individual_results <- data.frame()

for (i in seq_along(individual_measures)) {
  measure <- individual_measures[i]
  measure_name <- measure_names[i]

  if (!measure %in% names(landscape_condition)) {
    cat("    ", measure_name, ": variable not available\n")
    next
  }

  data_measure <- landscape_condition %>% filter(!is.na(.data[[measure]]))

  if (nrow(data_measure) < 20) {
    cat("    ", measure_name, ": insufficient data (n =", nrow(data_measure), ")\n")
    next
  }

  # Model: individual measure ~ volume + site
  formula_str <- paste(measure, "~ log_volume_scaled + site")
  m_individual <- lm(as.formula(formula_str), data = data_measure)
  s_individual <- summary(m_individual)

  # Extract volume coefficient
  vol_coef <- coef(s_individual)["log_volume_scaled", ]

  individual_results <- bind_rows(individual_results, data.frame(
    measure = measure_name,
    predictor = "log_volume",
    estimate = vol_coef["Estimate"],
    se = vol_coef["Std. Error"],
    t_value = vol_coef["t value"],
    p_value = vol_coef["Pr(>|t|)"],
    r2_adj = s_individual$adj.r.squared,
    n = nrow(data_measure),
    stringsAsFactors = FALSE
  ))

  cat("    ", measure_name, ":\n")
  cat("        Volume effect: β =", round(vol_coef["Estimate"], 3),
      ", t =", round(vol_coef["t value"], 2),
      ", p =", round(vol_coef["Pr(>|t|)"], 4), "\n")
  cat("        Adj. R² =", round(s_individual$adj.r.squared, 4), "\n")
}

cat("\n")

# H.3 Composite condition (PC1) vs volume + site (full sample)
cat("H.3 COMPOSITE CONDITION (PC1) vs LANDSCAPE - FULL SAMPLE\n")
cat("    --------------------------------------------------\n")

# Model 1: Volume + Site only (full sample)
m_vol_site <- lm(condition_score ~ log_volume_scaled + site, data = landscape_condition)
s_vol_site <- summary(m_vol_site)

cat("    Model 1: condition_score ~ log_volume + site (n =", n_landscape, ")\n")
cat("        Volume: β =", round(coef(m_vol_site)["log_volume_scaled"], 3),
    ", t =", round(s_vol_site$coefficients["log_volume_scaled", "t value"], 2),
    ", p =", round(s_vol_site$coefficients["log_volume_scaled", "Pr(>|t|)"], 4), "\n")
cat("        Adj. R² =", round(s_vol_site$adj.r.squared, 4), "\n")
cat("        AIC =", round(AIC(m_vol_site), 1), "\n\n")

# H.4 Composite condition vs neighborhood predictors (subset with neighborhood data)
cat("H.4 COMPOSITE CONDITION vs NEIGHBORHOOD - SUBSET WITH 5m SURVEYS\n")
cat("    --------------------------------------------------\n")

landscape_neighborhood <- landscape_condition %>%
  filter(!is.na(n_neighbors))

n_neighborhood_cond <- nrow(landscape_neighborhood)

if (n_neighborhood_cond >= 20) {

  # Model 2: Volume + Neighbors + Site
  m_vol_neigh <- lm(condition_score ~ log_volume_scaled + n_neighbors_scaled + site,
                    data = landscape_neighborhood)
  s_vol_neigh <- summary(m_vol_neigh)

  cat("    Model 2: condition ~ volume + n_neighbors + site (n =", n_neighborhood_cond, ")\n")
  cat("        Volume: β =", round(coef(m_vol_neigh)["log_volume_scaled"], 3),
      ", p =", round(s_vol_neigh$coefficients["log_volume_scaled", "Pr(>|t|)"], 4), "\n")
  cat("        N_neighbors: β =", round(coef(m_vol_neigh)["n_neighbors_scaled"], 3),
      ", p =", round(s_vol_neigh$coefficients["n_neighbors_scaled", "Pr(>|t|)"], 4), "\n")
  cat("        Adj. R² =", round(s_vol_neigh$adj.r.squared, 4), "\n")
  cat("        AIC =", round(AIC(m_vol_neigh), 1), "\n\n")

  # Model 3: Volume + Mean distance + Site
  m_vol_dist <- lm(condition_score ~ log_volume_scaled + mean_dist_scaled + site,
                   data = landscape_neighborhood)
  s_vol_dist <- summary(m_vol_dist)

  cat("    Model 3: condition ~ volume + mean_neighbor_dist + site (n =", n_neighborhood_cond, ")\n")
  cat("        Volume: β =", round(coef(m_vol_dist)["log_volume_scaled"], 3),
      ", p =", round(s_vol_dist$coefficients["log_volume_scaled", "Pr(>|t|)"], 4), "\n")
  cat("        Mean_dist: β =", round(coef(m_vol_dist)["mean_dist_scaled"], 3),
      ", p =", round(s_vol_dist$coefficients["mean_dist_scaled", "Pr(>|t|)"], 4), "\n")
  cat("        Adj. R² =", round(s_vol_dist$adj.r.squared, 4), "\n")
  cat("        AIC =", round(AIC(m_vol_dist), 1), "\n\n")

  # Model 4: Full landscape model (volume + neighbors + proximity + site)
  m_full_landscape <- lm(condition_score ~ log_volume_scaled + n_neighbors_scaled +
                           mean_dist_scaled + total_neighbor_vol_scaled + site,
                         data = landscape_neighborhood)
  s_full_landscape <- summary(m_full_landscape)

  cat("    Model 4: FULL LANDSCAPE MODEL (n =", n_neighborhood_cond, ")\n")
  cat("    condition ~ volume + n_neighbors + mean_dist + total_neighbor_vol + site\n")

  landscape_predictors <- c("log_volume_scaled", "n_neighbors_scaled",
                             "mean_dist_scaled", "total_neighbor_vol_scaled")
  predictor_labels <- c("Volume", "N_neighbors", "Mean_distance", "Total_neighbor_vol")

  for (j in seq_along(landscape_predictors)) {
    pred <- landscape_predictors[j]
    if (pred %in% rownames(coef(s_full_landscape))) {
      cat("        ", predictor_labels[j], ": β =",
          round(coef(m_full_landscape)[pred], 3),
          ", p =", round(s_full_landscape$coefficients[pred, "Pr(>|t|)"], 4), "\n")
    }
  }
  cat("        Adj. R² =", round(s_full_landscape$adj.r.squared, 4), "\n")
  cat("        AIC =", round(AIC(m_full_landscape), 1), "\n\n")

  # Model 5: Volume × Neighbors interaction
  m_vol_neigh_int <- lm(condition_score ~ log_volume_scaled * n_neighbors_scaled + site,
                        data = landscape_neighborhood)
  s_vol_neigh_int <- summary(m_vol_neigh_int)

  cat("    Model 5: condition ~ volume × n_neighbors + site (interaction)\n")
  cat("        Volume: β =", round(coef(m_vol_neigh_int)["log_volume_scaled"], 3),
      ", p =", round(s_vol_neigh_int$coefficients["log_volume_scaled", "Pr(>|t|)"], 4), "\n")
  cat("        N_neighbors: β =", round(coef(m_vol_neigh_int)["n_neighbors_scaled"], 3),
      ", p =", round(s_vol_neigh_int$coefficients["n_neighbors_scaled", "Pr(>|t|)"], 4), "\n")
  cat("        Interaction: β =", round(coef(m_vol_neigh_int)["log_volume_scaled:n_neighbors_scaled"], 3),
      ", p =", round(s_vol_neigh_int$coefficients["log_volume_scaled:n_neighbors_scaled", "Pr(>|t|)"], 4), "\n")
  cat("        Adj. R² =", round(s_vol_neigh_int$adj.r.squared, 4), "\n")
  cat("        AIC =", round(AIC(m_vol_neigh_int), 1), "\n\n")

} else {
  cat("    Insufficient neighborhood data for models 2-5 (n =", n_neighborhood_cond, ")\n\n")
  m_vol_neigh <- NULL
  m_vol_dist <- NULL
  m_full_landscape <- NULL
  m_vol_neigh_int <- NULL
}

# H.5 Site effects on condition
cat("H.5 SITE EFFECTS ON CONDITION\n")
cat("    --------------------------------------------------\n")

# ANOVA for site effect
m_site_only <- lm(condition_score ~ site, data = landscape_condition)
anova_site <- anova(m_site_only)
site_f <- anova_site["site", "F value"]
site_p <- anova_site["site", "Pr(>F)"]
site_r2 <- summary(m_site_only)$r.squared

cat("    Site-only model: condition ~ site\n")
cat("        F =", round(site_f, 2), ", p =", round(site_p, 4), "\n")
cat("        R² =", round(site_r2, 4), "(proportion explained by site alone)\n\n")

# Site means
site_means <- landscape_condition %>%
  group_by(site) %>%
  summarise(
    n = n(),
    mean_condition = mean(condition_score, na.rm = TRUE),
    sd_condition = sd(condition_score, na.rm = TRUE),
    se_condition = sd_condition / sqrt(n),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_condition))

cat("    Site-level condition means:\n")
for (i in 1:nrow(site_means)) {
  cat("        ", site_means$site[i], ": mean =", round(site_means$mean_condition[i], 3),
      "± SE", round(site_means$se_condition[i], 3),
      "(n =", site_means$n[i], ")\n")
}
cat("\n")

# H.6 Model comparison table
# NOTE: M1 uses full sample (n~112); M2-M5 use neighborhood subset (n~61).
# AIC values are only comparable within same-sample models (M2-M5).
cat("H.6 MODEL COMPARISON (AIC)\n")
cat("    --------------------------------------------------\n")
cat("    Note: Models have different sample sizes — AIC values not directly comparable.\n")
cat("    M1 uses full sample; M2-M5 use neighborhood subset only.\n")
cat("    --------------------------------------------------\n")

model_comparison <- data.frame(
  model = c("M1: Volume + Site",
            "M2: Volume + N_neighbors + Site",
            "M3: Volume + Mean_dist + Site",
            "M4: Full landscape",
            "M5: Volume × N_neighbors + Site"),
  sample = c(paste0("Full (n=", n_landscape, ")"),
             paste0("Neighborhood (n=", n_neighborhood_cond, ")"),
             paste0("Neighborhood (n=", n_neighborhood_cond, ")"),
             paste0("Neighborhood (n=", n_neighborhood_cond, ")"),
             paste0("Neighborhood (n=", n_neighborhood_cond, ")")),
  adj_r2 = NA_real_,
  aic = NA_real_,
  stringsAsFactors = FALSE
)

model_comparison$adj_r2[1] <- s_vol_site$adj.r.squared
model_comparison$aic[1] <- AIC(m_vol_site)

if (!is.null(m_vol_neigh)) {
  model_comparison$adj_r2[2] <- s_vol_neigh$adj.r.squared
  model_comparison$aic[2] <- AIC(m_vol_neigh)
}
if (!is.null(m_vol_dist)) {
  model_comparison$adj_r2[3] <- s_vol_dist$adj.r.squared
  model_comparison$aic[3] <- AIC(m_vol_dist)
}
if (!is.null(m_full_landscape)) {
  model_comparison$adj_r2[4] <- s_full_landscape$adj.r.squared
  model_comparison$aic[4] <- AIC(m_full_landscape)
}
if (!is.null(m_vol_neigh_int)) {
  model_comparison$adj_r2[5] <- s_vol_neigh_int$adj.r.squared
  model_comparison$aic[5] <- AIC(m_vol_neigh_int)
}

# Show comparison
cat("    Model                           Sample             Adj.R²    AIC\n")
cat("    ----------------------------------------------------------------\n")
for (i in 1:nrow(model_comparison)) {
  if (!is.na(model_comparison$aic[i])) {
    cat(sprintf("    %-30s %-20s %.4f    %.1f\n",
                model_comparison$model[i], model_comparison$sample[i],
                model_comparison$adj_r2[i], model_comparison$aic[i]))
  }
}
cat("\n")

# Best model (lowest AIC among neighborhood models)
if (any(!is.na(model_comparison$aic[2:5]))) {
  best_idx <- which.min(model_comparison$aic[2:5]) + 1
  cat("    Best model (lowest AIC):", model_comparison$model[best_idx], "\n")
  cat("    ΔAIC relative to M2 (additive):\n")
  for (i in 3:5) {
    if (!is.na(model_comparison$aic[i])) {
      delta_aic <- model_comparison$aic[i] - model_comparison$aic[2]
      cat("        ", model_comparison$model[i], ": ΔAIC =", round(delta_aic, 1), "\n")
    }
  }
}
cat("\n")

# H.7 Create landscape effects results table
landscape_effects_table <- data.frame(
  response = character(),
  predictor = character(),
  estimate = numeric(),
  se = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  adj_r2 = numeric(),
  n = integer(),
  model = character(),
  stringsAsFactors = FALSE
)

# Add individual measure results
if (nrow(individual_results) > 0) {
  individual_results$response <- individual_results$measure
  individual_results$model <- "Volume + Site"
  landscape_effects_table <- bind_rows(landscape_effects_table, individual_results)
}

# Add composite condition results
landscape_effects_table <- bind_rows(landscape_effects_table, data.frame(
  response = "PC1_Coral",
  predictor = "log_volume",
  estimate = coef(m_vol_site)["log_volume_scaled"],
  se = s_vol_site$coefficients["log_volume_scaled", "Std. Error"],
  t_value = s_vol_site$coefficients["log_volume_scaled", "t value"],
  p_value = s_vol_site$coefficients["log_volume_scaled", "Pr(>|t|)"],
  adj_r2 = s_vol_site$adj.r.squared,
  n = n_landscape,
  model = "Volume + Site",
  stringsAsFactors = FALSE
))

# Add neighborhood predictor results
if (!is.null(m_vol_neigh)) {
  landscape_effects_table <- bind_rows(landscape_effects_table, data.frame(
    response = "PC1_Coral",
    predictor = "n_neighbors",
    estimate = coef(m_vol_neigh)["n_neighbors_scaled"],
    se = s_vol_neigh$coefficients["n_neighbors_scaled", "Std. Error"],
    t_value = s_vol_neigh$coefficients["n_neighbors_scaled", "t value"],
    p_value = s_vol_neigh$coefficients["n_neighbors_scaled", "Pr(>|t|)"],
    adj_r2 = s_vol_neigh$adj.r.squared,
    n = n_neighborhood_cond,
    model = "Volume + Neighbors + Site",
    stringsAsFactors = FALSE
  ))
}

if (!is.null(m_vol_dist)) {
  landscape_effects_table <- bind_rows(landscape_effects_table, data.frame(
    response = "PC1_Coral",
    predictor = "mean_neighbor_dist",
    estimate = coef(m_vol_dist)["mean_dist_scaled"],
    se = s_vol_dist$coefficients["mean_dist_scaled", "Std. Error"],
    t_value = s_vol_dist$coefficients["mean_dist_scaled", "t value"],
    p_value = s_vol_dist$coefficients["mean_dist_scaled", "Pr(>|t|)"],
    adj_r2 = s_vol_dist$adj.r.squared,
    n = n_neighborhood_cond,
    model = "Volume + Mean_dist + Site",
    stringsAsFactors = FALSE
  ))
}

# Add site effect
landscape_effects_table <- bind_rows(landscape_effects_table, data.frame(
  response = "PC1_Coral",
  predictor = "site",
  estimate = NA,
  se = NA,
  t_value = NA,
  p_value = site_p,
  adj_r2 = site_r2,
  n = n_landscape,
  model = "Site only (ANOVA)",
  stringsAsFactors = FALSE
))

# Save landscape effects table
save_table(landscape_effects_table %>%
             mutate(across(where(is.numeric), ~round(., 4))),
           "landscape_condition_effects")
cat("H.7 Saved: landscape_condition_effects.csv\n")

# Save model comparison table
save_table(model_comparison %>%
             mutate(across(where(is.numeric), ~round(., 4))),
           "landscape_model_comparison")
cat("    Saved: landscape_model_comparison.csv\n")

# Save site means
save_table(site_means %>%
             mutate(across(where(is.numeric), ~round(., 4))),
           "site_condition_means")
cat("    Saved: site_condition_means.csv\n\n")

# H.8 Create landscape effects visualization
cat("H.8 Creating landscape effects visualization...\n")

# Panel A: Condition vs Volume by site
p_vol_condition <- ggplot(landscape_condition, aes(x = log_volume, y = condition_score)) +
  geom_point(aes(fill = site), shape = 21, color = "gray30", stroke = 0.4, alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  scale_fill_manual(values = SITE_COLORS, name = "Site") +
  labs(
    title = "A. Coral Size → Condition",
    subtitle = sprintf("β = %.3f, p = %.4f (n = %d)",
                       coef(m_vol_site)["log_volume_scaled"],
                       s_vol_site$coefficients["log_volume_scaled", "Pr(>|t|)"],
                       n_landscape),
    x = expression(ln(Volume~cm^3)),
    y = expression(Coral~condition~(PC1[Coral]))
  ) +
  theme_publication() +
  theme(legend.position = c(0.15, 0.85))

# Panel B: Condition vs N_neighbors (if available)
if (!is.null(m_vol_neigh)) {
  p_neigh_condition <- ggplot(landscape_neighborhood,
                               aes(x = n_neighbors, y = condition_score)) +
    geom_point(aes(fill = site), shape = 21, color = "gray30", stroke = 0.4, alpha = 0.7, size = 2.5) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
    scale_fill_manual(values = SITE_COLORS, name = "Site") +
    labs(
      title = "B. Neighborhood Density → Condition",
      subtitle = sprintf("β = %.3f, p = %.4f (n = %d)",
                         coef(m_vol_neigh)["n_neighbors_scaled"],
                         s_vol_neigh$coefficients["n_neighbors_scaled", "Pr(>|t|)"],
                         n_neighborhood_cond),
      x = "Number of neighbors (within 5m)",
      y = expression(Coral~condition~(PC1[Coral]))
    ) +
    theme_publication() +
    theme(legend.position = "none")
} else {
  p_neigh_condition <- ggplot() + theme_void() +
    labs(title = "B. Insufficient neighborhood data")
}

# Panel C: Condition by Site
p_site_condition <- ggplot(landscape_condition, aes(x = site, y = condition_score, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  geom_jitter(alpha = 0.3, width = 0.15, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  scale_fill_manual(values = SITE_COLORS) +
  labs(
    title = "C. Condition by Site",
    subtitle = sprintf("ANOVA: F = %.2f, p = %.4f", site_f, site_p),
    x = "Site",
    y = expression(Coral~condition~(PC1[Coral]))
  ) +
  theme_publication() +
  theme(legend.position = "none")

# Panel D: Individual measures heatmap (effect sizes)
if (nrow(individual_results) > 0) {
  individual_results$significant <- individual_results$p_value < 0.05

  p_individual <- ggplot(individual_results,
                          aes(x = measure, y = predictor, fill = estimate)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.2f\n(p=%.3f)", estimate, p_value)),
              size = 3, color = "white") +
    scale_fill_gradient2(low = "#D55E00", mid = "gray80", high = "#009E73",
                          midpoint = 0, name = "Effect\nsize (β)") +
    labs(
      title = "D. Individual Condition Measures",
      subtitle = "Effect of log(volume) on each physiological trait",
      x = "Physiological measure (position-corrected)",
      y = "Predictor"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
} else {
  p_individual <- ggplot() + theme_void() +
    labs(title = "D. Individual measures not available")
}

# Combine panels
p_landscape_effects <- (p_vol_condition | p_neigh_condition) /
                        (p_site_condition | p_individual) +
  plot_annotation(
    title = "Landscape Effects on Coral Condition (No CAFI Predictors)",
    subtitle = "Abiotic and spatial drivers of coral physiological health",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray40")
    )
  )

save_figure(p_landscape_effects,
           file.path(fig_dir, "landscape_condition_effects.png"),
           width = 12, height = 10)
cat("    Saved: landscape_condition_effects.png\n\n")

# H.9 Summary interpretation
cat("H.9 LANDSCAPE EFFECTS INTERPRETATION:\n")
cat("    =====================================================================\n")
cat("    Key findings (CAFI-independent effects on condition):\n\n")

# Volume effect
vol_p <- s_vol_site$coefficients["log_volume_scaled", "Pr(>|t|)"]
vol_b <- coef(m_vol_site)["log_volume_scaled"]
if (vol_p < 0.05) {
  vol_direction <- ifelse(vol_b > 0, "POSITIVE", "NEGATIVE")
  cat("    1. CORAL SIZE: SIGNIFICANT ", vol_direction, " effect\n")
  cat("       β =", round(vol_b, 3), ", p =", round(vol_p, 4), "\n")
  if (vol_b > 0) {
    cat("       → Larger corals have BETTER physiological condition\n")
  } else {
    cat("       → Larger corals have WORSE physiological condition\n")
  }
} else {
  cat("    1. CORAL SIZE: Not significant (p =", round(vol_p, 4), ")\n")
  cat("       → Size does not predict condition after controlling for site\n")
}
cat("\n")

# Neighborhood effect
if (!is.null(m_vol_neigh)) {
  neigh_p <- s_vol_neigh$coefficients["n_neighbors_scaled", "Pr(>|t|)"]
  neigh_b <- coef(m_vol_neigh)["n_neighbors_scaled"]
  if (neigh_p < 0.05) {
    neigh_direction <- ifelse(neigh_b > 0, "POSITIVE", "NEGATIVE")
    cat("    2. NEIGHBORHOOD DENSITY: SIGNIFICANT ", neigh_direction, " effect\n")
    cat("       β =", round(neigh_b, 3), ", p =", round(neigh_p, 4), "\n")
    if (neigh_b > 0) {
      cat("       → More neighbors → BETTER coral condition (facilitation?)\n")
    } else {
      cat("       → More neighbors → WORSE coral condition (competition?)\n")
    }
  } else {
    cat("    2. NEIGHBORHOOD DENSITY: Not significant (p =", round(neigh_p, 4), ")\n")
    cat("       → Local coral density does not predict condition\n")
  }
} else {
  cat("    2. NEIGHBORHOOD DENSITY: Not tested (insufficient data)\n")
}
cat("\n")

# Site effect
if (site_p < 0.05) {
  cat("    3. SITE: SIGNIFICANT effect (F =", round(site_f, 2), ", p =", round(site_p, 4), ")\n")
  cat("       → Reef location explains", round(site_r2 * 100, 1), "% of condition variance\n")
  best_site <- site_means$site[1]
  worst_site <- site_means$site[nrow(site_means)]
  cat("       → Best condition:", best_site, "(mean =", round(site_means$mean_condition[1], 2), ")\n")
  cat("       → Worst condition:", worst_site, "(mean =",
      round(site_means$mean_condition[nrow(site_means)], 2), ")\n")
} else {
  cat("    3. SITE: Not significant (p =", round(site_p, 4), ")\n")
  cat("       → Reef sites do not differ in average coral condition\n")
}
cat("\n")

cat("    CONCLUSIONS:\n")
cat("    - Landscape factors (size, neighborhood, site) explain coral condition\n")
cat("    - These effects are INDEPENDENT of CAFI community composition\n")
cat("    - CAFI effects (Part A-F) operate ON TOP of these landscape baselines\n")
cat("    =====================================================================\n\n")

# ============================================================================
# PART I: GLOBAL FDR CORRECTION ACROSS ALL TEST FAMILIES
# ============================================================================
# NOTE: Family-wise FDR is conservative (corrects within predictor family).
# Global FDR (all tests pooled) is also computed below as sensitivity check.
# Results reported as "nominal" must explicitly note global FDR threshold.

cat("============================================================\n")
cat("PART I: GLOBAL FDR CORRECTION ACROSS ALL TEST FAMILIES\n")
cat("============================================================\n\n")

# Collect p-values from the three test families
# All forward models use OLS p-values (primary inference; BP confirms no heteroscedasticity)
# HC3 p-values retained in supplement sensitivity
# Reverse models (Condition->CAFI) use OLS p-values (NB GLMs, no HC3)
cafi_to_condition_pvals <- cafi_to_condition_df$p_value
condition_to_cafi_pvals <- condition_to_cafi_df$p_value
key_species_pvals <- if (exists("key_species_df") && nrow(key_species_df) > 0) {
  key_species_df$p_value  # OLS primary
} else {
  numeric(0)
}

all_raw_pvalues <- c(cafi_to_condition_pvals, condition_to_cafi_pvals, key_species_pvals)
all_test_labels <- c(
  paste0("CAFI->Cond: ", cafi_to_condition_df$predictor),
  paste0("Cond->CAFI: ", condition_to_cafi_df$response),
  if (exists("key_species_df") && nrow(key_species_df) > 0) {
    paste0("Key sp: ", key_species_df$predictor)
  } else {
    character(0)
  }
)

all_fdr_global <- p.adjust(all_raw_pvalues, method = "BH")

cat("Global FDR correction (all", length(all_raw_pvalues), "tests across 3 families):\n")
cat("  Families: CAFI->Condition (", length(cafi_to_condition_pvals), "), ",
    "Condition->CAFI (", length(condition_to_cafi_pvals), "), ",
    "Key Species (", length(key_species_pvals), ")\n", sep = "")
cat("  Number significant at raw p < 0.05:", sum(all_raw_pvalues < 0.05), "\n")
cat("  Number significant at family-wise FDR < 0.05:",
    sum(cafi_to_condition_df$p_fdr < 0.05) +
    sum(condition_to_cafi_df$p_fdr < 0.05) +
    (if (exists("key_species_df") && nrow(key_species_df) > 0) sum(key_species_df$p_adj < 0.05) else 0), "\n")
cat("  Number significant at GLOBAL FDR < 0.05:", sum(all_fdr_global < 0.05), "\n\n")

# Show any tests that change significance under global vs family-wise FDR
if (any(all_raw_pvalues < 0.05)) {
  cat("  Tests with raw p < 0.05:\n")
  sig_idx <- which(all_raw_pvalues < 0.05)
  for (idx in sig_idx) {
    cat(sprintf("    %-35s p_raw = %.4f, p_global_FDR = %.4f %s\n",
                all_test_labels[idx], all_raw_pvalues[idx], all_fdr_global[idx],
                ifelse(all_fdr_global[idx] < 0.05, "(survives global FDR)", "(does NOT survive global FDR)")))
  }
}
cat("\n")

# Store global FDR results for use in figures
global_fdr_results <- tibble(
  test_label = all_test_labels,
  p_raw = all_raw_pvalues,
  p_fdr_global = all_fdr_global,
  significant_global_fdr = all_fdr_global < 0.05
)

# ============================================================================
# 3. CREATE VISUALIZATIONS
# ============================================================================

cat("============================================================\n")
cat("CREATING VISUALIZATIONS\n")
cat("============================================================\n\n")

# --- Panel A: Richness-Abundance Artifact (KEY FINDING) ---
# Load the comparison data from Part A2
richness_comparison_df <- tryCatch(
  load_object("richness_comparison_results"),
  error = function(e) NULL
)

if (!is.null(richness_comparison_df)) {
  p_richness_artifact <- ggplot(richness_comparison_df,
                                  aes(x = type, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
    geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper,
                         color = significant, fill = significant),
                    size = 1.2, linewidth = 1.2, shape = 21, stroke = 0.4) +
    geom_text(aes(y = ci_upper, label = sprintf("  p = %.3f, r = %.2f", p_value, cor_with_abundance)),
              hjust = 0, size = 3, color = "gray20") +
    scale_color_manual(values = c("TRUE" = "#2e7d32", "FALSE" = "#757575"),
                       guide = "none") +
    scale_fill_manual(values = c("TRUE" = "#2e7d32", "FALSE" = "white"),
                      guide = "none") +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(
      title = "A. Species Richness \u2192 Condition: Abundance Artifact",
      subtitle = "Raw richness signal disappears after rarefaction (r = correlation with abundance)",
      x = "",
      y = expression("Effect on condition score (" * beta * ")")
    ) +
    theme_multipanel() +
    theme(axis.text.y = element_text(size = 9),
          plot.margin = margin(5, 20, 5, 5, "mm"))
} else {
  # Fallback if richness comparison not available
  cat("  WARNING: richness_comparison_results not found, using PC1_CAFI fallback for Panel A\n")
  p_richness_artifact <- ggplot(analysis_data %>% filter(!is.na(pc1_cafi)),
                                 aes(x = pc1_cafi, y = condition_score)) +
    geom_point(aes(fill = site), shape = 21, color = "gray30", stroke = 0.4, alpha = 0.7, size = 2.5) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
    scale_fill_manual(values = SITE_COLORS, name = "Site") +
    labs(title = "A. PC1_CAFI → Condition",
         x = expression(PC1[CAFI]), y = expression(Coral~condition~(PC1[Coral]))) +
    theme_multipanel()
}

# --- Panel A: Raw richness vs Condition (BEF hypothesis) ---
# Retrieve Hochberg-corrected p-value for annotation
p_richness_hochberg <- cafi_to_condition_df %>%
  filter(predictor == "Species richness") %>%
  pull(p_corrected)
p_richness_ols <- cafi_to_condition_df %>%
  filter(predictor == "Species richness") %>%
  pull(p_value)

p_raw_richness <- ggplot(analysis_data %>% filter(!is.na(otu_richness), !is.na(condition_score)),
                         aes(x = otu_richness, y = condition_score)) +
  geom_point(aes(fill = site), shape = 21, color = "gray30", stroke = 0.4, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE, alpha = 0.2) +
  scale_fill_manual(values = SITE_COLORS, name = "Site") +
  labs(x = "Species richness", y = expression(Condition~score~(PC1[Coral]))) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("Hochberg p = %.3f", p_richness_hochberg),
           size = 3, fontface = "bold", color = "gray20") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 3.0,
           label = sprintf("r(richness, abundance) = %.2f", cor_raw_abund),
           size = 2.8, fontface = "italic", color = "gray40") +
  theme_multipanel()

# --- Panel B: Abundance vs Condition scatter (a priori BEF, Hochberg k=2) ---
p_abundance_hochberg <- cafi_to_condition_df %>%
  filter(predictor == "Total CAFI") %>%
  pull(p_corrected)

p_abundance_scatter <- ggplot(analysis_data %>% filter(!is.na(total_cafi), !is.na(condition_score)),
                         aes(x = total_cafi, y = condition_score)) +
  geom_point(aes(fill = site), shape = 21, color = "gray30", stroke = 0.4, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", formula = y ~ sqrt(x), color = "black", se = TRUE, alpha = 0.2) +
  scale_x_sqrt(breaks = c(0, 25, 50, 100, 150, 200, 250)) +
  scale_fill_manual(values = SITE_COLORS, guide = "none") +
  labs(x = "Total CAFI abundance",
       y = expression(Condition~score~(PC1[Coral]))) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("Hochberg p = %.3f", p_abundance_hochberg),
           size = 3, fontface = "bold", color = "gray20") +
  theme_multipanel()

# --- A priori BEF forest plot (moved to supplement S13) ---
apriori_forest_data <- cafi_to_condition_df %>%
  filter(hypothesis_type == "a_priori") %>%
  mutate(
    predictor_label = case_when(
      predictor == "Species richness" ~ "Species\nrichness",
      predictor == "Total CAFI" ~ "Total CAFI\nabundance"
    ),
    predictor_label = factor(predictor_label,
                              levels = c("Total CAFI\nabundance", "Species\nrichness")),
    sig_color = ifelse(p_corrected < 0.05, "#2e7d32", "#757575"),
    p_label = sprintf("p = %.3f\nHochberg = %.3f", p_value, p_corrected)
  )

p_apriori_forest <- ggplot(apriori_forest_data,
                            aes(x = predictor_label, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, color = sig_color),
                  size = 1.5, linewidth = 1.2, show.legend = FALSE) +
  scale_color_identity() +
  geom_text(aes(y = ci_upper, label = p_label),
            hjust = -0.1, size = 2.6, color = "gray20", lineheight = 0.85) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.40))) +
  labs(title = "A priori BEF hypotheses (Hochberg, k = 2)",
       x = "", y = "Effect on condition (standardized \u03B2)") +
  theme_multipanel() +
  theme(legend.position = "none",
        plot.margin = margin(5, 8, 5, 5, "mm"))

# --- Panel C: Rarefied richness vs Condition (diagnostic) ---
# Load rarefied data saved in Part A2
analysis_data_rare_plot <- tryCatch(
  load_object("analysis_data_rarefied"),
  error = function(e) {
    cat("  Note: Using in-memory rarefied data (RDS not found)\n")
    if (exists("analysis_data_rare")) analysis_data_rare
    else { cat("  ERROR: No rarefied data available for Panel C\n"); NULL }
  }
)

cor_rare_abund_plot <- cor(analysis_data_rare_plot$rarefied_richness,
                           analysis_data_rare_plot$total_cafi, use = "complete")

p_rare_richness <- ggplot(analysis_data_rare_plot %>%
                             filter(!is.na(rarefied_richness), !is.na(condition_score)),
                          aes(x = rarefied_richness, y = condition_score)) +
  geom_point(aes(fill = site), shape = 21, color = "gray30", stroke = 0.4, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE, alpha = 0.2) +
  scale_fill_manual(values = SITE_COLORS, guide = "none") +
  labs(title = "C", subtitle = "Rarefied richness (abundance equalized)",
       x = "Expected species richness (rarefied, n = 20)", y = expression(Condition~score~(PC1[Coral]))) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("p = %.2f (n.s.)", p_rare),
           size = 3, color = "gray40") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 3.0,
           label = sprintf("r(richness, abundance) = %.2f", cor_rare_abund_plot),
           size = 2.8, color = "gray45") +
  theme_multipanel()

# --- Panel C: Trapezia (defenders) vs Condition ---
p_trapezia <- ggplot(analysis_data, aes(x = n_trapezia, y = condition_score)) +
  geom_jitter(aes(fill = site), shape = 21, color = "gray30", stroke = 0.4, alpha = 0.6, width = 0.15, size = 2.5) +
  geom_smooth(method = "lm", color = "gray30", se = TRUE, linewidth = 1.2, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  scale_fill_manual(values = SITE_COLORS, name = "Site") +
  labs(
    title = "B. Trapezia (Defenders)",
    subtitle = "Expected: positive (predator defense)",
    x = "Trapezia abundance",
    y = expression(Coral~condition~(PC1[Coral]))
  ) +
  theme_multipanel() +
  theme(legend.position = "none")

# --- Panel D: Galeropsis vs Condition ---
p_galeropsis <- ggplot(analysis_data, aes(x = n_galeropsis, y = condition_score)) +
  geom_jitter(aes(fill = site), shape = 21, color = "gray30", stroke = 0.4, alpha = 0.6, width = 0.15, size = 2.5) +
  geom_smooth(method = "lm", color = "gray30", se = TRUE, linewidth = 1.2, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  scale_fill_manual(values = SITE_COLORS, name = "Site") +
  labs(
    title = "C. Galeropsis (Tissue Consumer)",
    subtitle = "Expected: negative (Coralliophilinae)",
    x = "Galeropsis abundance",
    y = expression(Coral~condition~(PC1[Coral]))
  ) +
  theme_multipanel() +
  theme(legend.position = "none")

# --- Panel D: Exploratory forest plot (functional groups, BH-FDR k=4) ---
exploratory_forest_data <- cafi_to_condition_df %>%
  filter(hypothesis_type == "exploratory") %>%
  mutate(
    predictor_label = case_when(
      predictor == "Shannon diversity" ~ "Shannon\ndiversity",
      predictor == "Trapezia abundance" ~ "Trapezia\n(crabs)",
      predictor == "Resident Fish abundance" ~ "Resident\nfish",
      predictor == "Galeropsis abundance" ~ "Galeropsis\n(snail)"
    ),
    predictor_label = factor(predictor_label,
                              levels = rev(c("Shannon\ndiversity", "Trapezia\n(crabs)",
                                             "Resident\nfish", "Galeropsis\n(snail)"))),
    sig_color = ifelse(p_corrected < 0.05, "#2e7d32", "#757575"),
    p_label = sprintf("p_FDR = %.3f", p_corrected)
  )

p_exploratory_forest <- ggplot(exploratory_forest_data,
                                aes(x = predictor_label, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, color = sig_color),
                  size = 1.2, linewidth = 1.0, show.legend = FALSE) +
  scale_color_identity() +
  geom_text(aes(y = ci_upper, label = paste0("  ", p_label)),
            hjust = 0, size = 2.8, color = "gray20") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(
    title = "D",
    subtitle = "Exploratory predictors (BH-FDR, k = 4)",
    x = "",
    y = "Effect on condition (standardized \u03B2)"
  ) +
  theme_multipanel() +
  theme(legend.position = "none",
        plot.margin = margin(5, 22, 5, 5, "mm"))

# --- Supplement panels: Trapezia scatter, Galeropsis scatter, bidirectional ---
# (kept for figS13, not in main figure)

# --- Panel H: Full bidirectional feedback panel (both directions, all predictors) ---
bidir_all_data <- bind_rows(
  cafi_to_condition_df %>%
    dplyr::select(predictor, estimate, ci_lower, ci_upper, p_value, p_fdr) %>%
    mutate(direction = "CAFI \u2192 Condition"),
  condition_to_cafi_df %>%
    dplyr::select(predictor = response, estimate, ci_lower, ci_upper, p_value, p_fdr) %>%
    mutate(direction = "Condition \u2192 CAFI")
)

p_bidir_full <- ggplot(bidir_all_data, aes(x = estimate, y = predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper, color = direction)) +
  scale_color_manual(values = c("CAFI \u2192 Condition" = "#0072B2",
                                "Condition \u2192 CAFI" = "#D55E00"),
                     name = "Direction") +
  labs(title = "Full Bidirectional Feedback", subtitle = "All predictors in both directions",
       x = "Effect size (coefficient)", y = "") +
  theme_multipanel() + theme(legend.position = "bottom")

# ---- Manuscript Figure 5: 2-panel BEF story ----
# Layout: A (richness scatter) | B (abundance scatter)
cat("\nAssembling manuscript Figure 5 (BEF story)...\n")
fig5_feedbacks <- p_raw_richness + p_abundance_scatter +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "none",
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = "topleft")

# Save manuscript figure (to both manuscript and analysis dirs)
save_figure(fig5_feedbacks,
           file.path(PATHS$fig_manuscript, "fig5_feedbacks.png"),
           width = 183, height = 100, units = "mm")
save_figure(fig5_feedbacks,
           file.path(fig_dir, "fig5_feedbacks.png"),
           width = 183, height = 100, units = "mm")
cat("   Saved: fig5_feedbacks.png (manuscript + analysis)\n")

# ---- Supplement Figure S14: Additional CAFI-condition panels ----
# ---- Supplement Figure S12: Variance partitioning + BEF diagnostics ----
cat("\nAssembling supplement Figure S12 (BEF variance partitioning)...\n")

# Panel S12-A: Stacked bar showing variance partitioning
vp_data <- tibble(
  component = factor(c("Unique to\nrichness", "Shared\n(confounded)", "Unique to\nabundance"),
                     levels = c("Unique to\nrichness", "Shared\n(confounded)", "Unique to\nabundance")),
  pct = c(100 * unique_richness / max(total_explained, 1e-10),
          100 * shared / max(total_explained, 1e-10),
          100 * unique_abundance / max(total_explained, 1e-10)),
  r2 = c(unique_richness, shared, unique_abundance)
)

p_varpart <- ggplot(vp_data, aes(x = "BEF variance", y = pct, fill = component)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%\n(\u0394R\u00B2 = %.4f)", pct, r2)),
            position = position_stack(vjust = 0.5), size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(values = c("Unique to\nrichness" = "#0072B2",
                                "Shared\n(confounded)" = "#999999",
                                "Unique to\nabundance" = "#E69F00"),
                    name = "") +
  coord_flip() +
  labs(title = "A", subtitle = "Variance partitioning: richness vs abundance",
       x = "", y = "% of total explained variance") +
  theme_multipanel() +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank())

# Panel S12-B: Partial regression scatter (richness | abundance)
# Compute partial residuals
m_rich_resid <- lm(otu_richness ~ sqrt(total_cafi) + log_volume + site, data = bef_data)
m_cond_resid <- lm(condition_score ~ sqrt(total_cafi) + log_volume + site, data = bef_data)

partial_data <- tibble(
  richness_resid = residuals(m_rich_resid),
  condition_resid = residuals(m_cond_resid),
  site = bef_data$site
)

p_partial <- ggplot(partial_data, aes(x = richness_resid, y = condition_resid)) +
  geom_point(aes(fill = site), shape = 21, color = "gray30", stroke = 0.4, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE, alpha = 0.2) +
  scale_fill_manual(values = SITE_COLORS, guide = "none") +
  labs(title = "B", subtitle = "Partial regression: richness | abundance",
       x = "Richness residuals (| abundance + volume + site)",
       y = "Condition residuals (| abundance + volume + site)") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("p = %.3f", richness_partial_p),
           size = 3, fontface = "bold", color = "gray20") +
  theme_multipanel()

# Panel S12-C: Path model diagram (text-based visualization)
path_text <- tibble(
  from = c("Volume", "Volume", "Volume", "Richness", "Abundance"),
  to = c("Richness", "Abundance", "Condition", "Condition", "Condition"),
  beta = c(coef(path_rich)["log_volume_z"],
           coef(path_abund)["log_volume_z"],
           coef(path_cond)["log_volume_z"],
           coef(path_cond)["richness_z"],
           coef(path_cond)["abundance_z"]),
  p_val = c(summary(path_rich)$coefficients["log_volume_z", "Pr(>|t|)"],
            summary(path_abund)$coefficients["log_volume_z", "Pr(>|t|)"],
            summary(path_cond)$coefficients["log_volume_z", "Pr(>|t|)"],
            summary(path_cond)$coefficients["richness_z", "Pr(>|t|)"],
            summary(path_cond)$coefficients["abundance_z", "Pr(>|t|)"]),
  sig = NA
)
path_text$sig <- ifelse(path_text$p_val < 0.05, "p < 0.05", "n.s.")
path_text$label <- sprintf("%s \u2192 %s", path_text$from, path_text$to)
path_text$label <- factor(path_text$label, levels = rev(path_text$label))

p_path <- ggplot(path_text, aes(x = label, y = beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  geom_col(aes(fill = sig), width = 0.6) +
  scale_fill_manual(values = c("p < 0.05" = "#0072B2", "n.s." = "#CCCCCC"), name = "") +
  geom_text(aes(label = sprintf("\u03B2 = %.2f", beta)),
            hjust = ifelse(path_text$beta >= 0, -0.1, 1.1), size = 3) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.25))) +
  labs(title = "C", subtitle = "Path model coefficients (standardized)",
       x = "", y = "Standardized \u03B2") +
  theme_multipanel() +
  theme(legend.position = "bottom")

figS12_bef <- (p_varpart) /
               (p_partial | p_path) +
  plot_layout(heights = c(0.8, 1.2))

save_figure(figS12_bef,
           file.path(PATHS$fig_supplement, "figS12_bef_variance_partitioning.png"),
           width = 230, height = 220, units = "mm")
save_figure(figS12_bef,
           file.path(fig_dir, "figS12_bef_variance_partitioning.png"),
           width = 230, height = 220, units = "mm")
cat("   Saved: figS12_bef_variance_partitioning.png (supplement + analysis)\n")

# ---- Supplement Figure S13: Rarefied richness, exploratory predictors, species scatters ----
cat("\nAssembling supplement Figure S13 (additional panels)...\n")

# Update panel labels for supplement context
p_apriori_forest_supp <- p_apriori_forest + labs(title = "A")
p_rare_richness_supp <- p_rare_richness + labs(title = "B")
p_exploratory_forest_supp <- p_exploratory_forest + labs(title = "C")

p_bidir_supp <- p_bidir_full + labs(title = "F") + theme(legend.position = "bottom")

figS13 <- (p_apriori_forest_supp | p_rare_richness_supp) /
           (p_exploratory_forest_supp + labs(title = "C") |
            p_trapezia + labs(title = "D")) /
           (p_galeropsis + labs(title = "E") | p_bidir_supp) +
  plot_layout(heights = c(1, 1, 1))

save_figure(figS13,
           file.path(PATHS$fig_supplement, "figS13_condition_details.png"),
           width = 230, height = 300, units = "mm")
save_figure(figS13,
           file.path(fig_dir, "figS13_condition_details.png"),
           width = 230, height = 300, units = "mm")
cat("   Saved: figS13_condition_details.png (supplement + analysis)\n")

# ============================================================================
# FIGURE S16: SPECIES × TRAIT HEATMAP
# ============================================================================
# Shows standardized β for each species × condition trait combination.
# Panel A: β values with significance markers; Panel B: FDR-adjusted p-values.
# ============================================================================

cat("\n--- Figure S16: Species × trait heatmap ---\n")

# Build species × trait data from existing analysis objects
# Use the 25 species from top_species_for_corr (already defined above)

# Trait mapping
heatmap_trait_cols <- c("condition_score", "protein_corr", "carb_corr", "zoox_corr", "afdw_corr")
heatmap_trait_labels <- c("Condition (PC1)", "Protein", "Carbohydrate", "Zooxanthellae", "AFDW")

# Build per-coral species abundance matrix
heatmap_sp_wide <- cafi_clean %>%
  filter(coral_id %in% analysis_data$coral_id,
         otu %in% top_species_for_corr$otu) %>%
  group_by(coral_id, otu) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = otu, values_from = n, values_fill = 0)

# Merge with condition + physiology data
heatmap_data <- analysis_data %>%
  dplyr::select(coral_id, site, log_volume,
                condition_score,
                any_of(c("protein_corr", "carb_corr", "zoox_corr", "afdw_corr"))) %>%
  left_join(heatmap_sp_wide, by = "coral_id")

# Fill NAs in species counts with 0
for (sp in top_species_for_corr$otu) {
  if (sp %in% names(heatmap_data)) {
    heatmap_data[[sp]] <- ifelse(is.na(heatmap_data[[sp]]), 0, heatmap_data[[sp]])
  }
}

# Assign taxonomic group labels
sp_group_lookup <- cafi_clean %>%
  filter(otu %in% top_species_for_corr$otu) %>%
  distinct(otu, type) %>%
  mutate(group_label = case_when(
    grepl("Trapezia", otu) ~ "Trapezia",
    type == "fish" ~ "Fish",
    type == "shrimp" ~ "Shrimp",
    type == "snail" ~ "Gastropod",
    type == "crab" ~ "Hermit crab",
    otu == "Breviturma pica" | otu == "Ophiocomella sexradia" ~ "Brittle star",
    otu == "Polynoidae" ~ "Polychaete",
    TRUE ~ "Other"
  )) %>%
  dplyr::select(otu, group_label)

# Get family info
sp_family_lookup <- cafi_clean %>%
  filter(otu %in% top_species_for_corr$otu) %>%
  distinct(otu, family)

# Run models: trait ~ sqrt(n_species) + log(volume) + site
heatmap_results <- data.frame()
available_heatmap_traits <- heatmap_trait_cols[heatmap_trait_cols %in% names(heatmap_data)]

for (sp in top_species_for_corr$otu) {
  if (!sp %in% names(heatmap_data)) next

  sp_data <- heatmap_data
  sp_data$sp_count <- sp_data[[sp]]
  n_present <- sum(sp_data$sp_count > 0)
  total_ind <- sum(sp_data$sp_count)

  if (n_present < 10) next

  for (j in seq_along(available_heatmap_traits)) {
    trait <- available_heatmap_traits[j]
    label <- heatmap_trait_labels[j]

    model_data <- sp_data %>% filter(!is.na(.data[[trait]]))
    if (nrow(model_data) < 15) next

    tryCatch({
      m <- lm(as.formula(paste(trait, "~ sqrt(sp_count) + log_volume + site")),
              data = model_data)
      s <- summary(m)
      coef_row <- s$coefficients["sqrt(sp_count)", , drop = FALSE]

      grp <- sp_group_lookup$group_label[sp_group_lookup$otu == sp]
      if (length(grp) == 0) grp <- "Other"
      fam <- sp_family_lookup$family[sp_family_lookup$otu == sp]
      if (length(fam) == 0) fam <- "Unknown"

      heatmap_results <- bind_rows(heatmap_results, data.frame(
        species = sp, group_label = grp, family = fam, trait = label,
        n_present = n_present, total_ind = total_ind,
        beta = coef_row[1, "Estimate"], se = coef_row[1, "Std. Error"],
        t_val = coef_row[1, "t value"], p_raw = coef_row[1, "Pr(>|t|)"],
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
}

# FDR correction across all tests
if (nrow(heatmap_results) > 0) {
  heatmap_results$p_fdr <- p.adjust(heatmap_results$p_raw, method = "BH")
}

# Save heatmap data table
save_table(heatmap_results, "species_trait_heatmap_data")
cat("   Saved: species_trait_heatmap_data.csv (", nrow(heatmap_results), " tests)\n")

# --- Build heatmap figure ---
if (nrow(heatmap_results) > 0) {

  # Order species by group then by Condition PC1 beta
  pc1_betas <- heatmap_results %>%
    filter(trait == "Condition (PC1)") %>%
    dplyr::select(species, group_label, beta_pc1 = beta)

  species_order <- pc1_betas %>%
    arrange(group_label, desc(beta_pc1)) %>%
    pull(species)

  heatmap_results$species <- factor(heatmap_results$species, levels = rev(species_order))
  heatmap_results$trait <- factor(heatmap_results$trait,
    levels = c("Condition (PC1)", "Protein", "Carbohydrate", "Zooxanthellae", "AFDW"))

  # Panel A: Beta heatmap with significance markers
  heatmap_results$sig_marker <- ifelse(heatmap_results$p_raw < 0.05, "*", "")
  heatmap_results$beta_label <- sprintf("%.2f%s", heatmap_results$beta, heatmap_results$sig_marker)
  heatmap_results$text_color <- ifelse(heatmap_results$p_raw < 0.05, "black", "grey60")

  # Split for bold/plain rendering
  hm_sig <- heatmap_results %>% filter(p_raw < 0.05)
  hm_ns <- heatmap_results %>% filter(p_raw >= 0.05)

  # Get group labels for faceting
  heatmap_results <- heatmap_results %>%
    left_join(pc1_betas %>% dplyr::select(species, group_label), by = "species",
              suffix = c("", ".dup")) %>%
    mutate(group_label = coalesce(group_label, group_label.dup)) %>%
    dplyr::select(-any_of("group_label.dup"))

  hm_sig <- hm_sig %>%
    left_join(pc1_betas %>% dplyr::select(species, group_label), by = "species",
              suffix = c("", ".dup")) %>%
    mutate(group_label = coalesce(group_label, group_label.dup)) %>%
    dplyr::select(-any_of("group_label.dup"))

  hm_ns <- hm_ns %>%
    left_join(pc1_betas %>% dplyr::select(species, group_label), by = "species",
              suffix = c("", ".dup")) %>%
    mutate(group_label = coalesce(group_label, group_label.dup)) %>%
    dplyr::select(-any_of("group_label.dup"))

  p_heatmap_A <- ggplot(heatmap_results, aes(x = trait, y = species, fill = beta)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(data = hm_sig,
              aes(label = sprintf("%.2f*", beta)),
              size = 2.5, fontface = "bold", color = "black") +
    geom_text(data = hm_ns,
              aes(label = sprintf("%.2f", beta)),
              size = 2.3, color = "grey50") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "Std. beta",
                         limits = c(-0.7, 0.8)) +
    facet_grid(group_label ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(face = "italic", size = 7),
      strip.text.y = element_text(angle = 0, size = 7, face = "bold", hjust = 0),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(15, "mm"),
      legend.key.height = unit(3, "mm"),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )

  # Panel B: FDR p-value heatmap (only show where < 0.10)
  hm_fdr <- heatmap_results %>%
    mutate(p_fdr_show = ifelse(p_fdr < 0.10, p_fdr, NA))

  p_heatmap_B <- ggplot(hm_fdr, aes(x = trait, y = species, fill = -log10(p_fdr_show))) +
    geom_tile(data = hm_fdr %>% filter(!is.na(p_fdr_show)),
              color = "white", linewidth = 0.5) +
    geom_tile(data = hm_fdr %>% filter(is.na(p_fdr_show)),
              fill = "grey95", color = "white", linewidth = 0.5) +
    geom_text(data = hm_fdr %>% filter(!is.na(p_fdr_show)),
              aes(label = sprintf("%.3f", p_fdr_show)),
              size = 2.5, color = "black") +
    scale_fill_viridis_c(option = "inferno", direction = -1, na.value = "grey95",
                         name = expression(-log[10](p[FDR]))) +
    facet_grid(group_label ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL, subtitle = "Only cells with FDR p < 0.10 shown; grey = not significant") +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(face = "italic", size = 7),
      strip.text.y = element_text(angle = 0, size = 7, face = "bold", hjust = 0),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(15, "mm"),
      legend.key.height = unit(3, "mm"),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )

  # Combine heatmap panels
  figS16_heatmap <- p_heatmap_A + p_heatmap_B +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 11, face = "bold"))

  save_figure(figS16_heatmap,
             file.path(PATHS$fig_supplement, "figS16_species_trait_heatmap.png"),
             width = 250, height = 200, units = "mm")
  save_figure(figS16_heatmap,
             file.path(fig_dir, "figS16_species_trait_heatmap.png"),
             width = 250, height = 200, units = "mm")
  cat("   Saved: figS16_species_trait_heatmap.png\n")
}

# ============================================================================
# FIGURE S17: SPECIES × TRAIT BIPLOTS — MULTI-PAGE PDFs
# ============================================================================
# One multi-page PDF per condition trait (5 total), showing ALL species.
# 9 panels per page (3×3 grid), ordered by taxonomic group then prevalence.
# Matches companion experiment Figure 7 style (Stier et al.).
# ============================================================================

cat("\n--- Figure S17: Species × trait multi-page PDFs ---\n")

if (nrow(heatmap_results) > 0) {

  # Map trait labels to column names and y-axis labels
  trait_col_map <- setNames(
    c("condition_score", "protein_corr", "carb_corr", "zoox_corr", "afdw_corr"),
    c("Condition (PC1)", "Protein", "Carbohydrate", "Zooxanthellae", "AFDW")
  )

  trait_ylab_map <- list(
    "Condition (PC1)" = expression(Condition~score~(PC1[Coral])),
    "Protein" = expression(Protein~(position~corrected)),
    "Carbohydrate" = expression(Carbohydrate~(position~corrected)),
    "Zooxanthellae" = expression(Zooxanthellae~(position~corrected)),
    "AFDW" = expression(AFDW~(position~corrected))
  )

  trait_file_map <- c(
    "Condition (PC1)" = "condition_pc1",
    "Protein" = "protein",
    "Carbohydrate" = "carbohydrate",
    "Zooxanthellae" = "zooxanthellae",
    "AFDW" = "afdw"
  )

  # Get ordered species list: by group then prevalence (descending)
  species_ordered <- heatmap_results %>%
    filter(trait == levels(heatmap_results$trait)[1] |
             trait == as.character(unique(heatmap_results$trait))[1]) %>%
    distinct(species, group_label, n_present) %>%
    mutate(species = as.character(species)) %>%
    arrange(group_label, desc(n_present)) %>%
    pull(species)

  cat("   Species to plot:", length(species_ordered), "\n")

  # Helper: build one biplot panel
  make_biplot <- function(sp, trait_label, trait_col, beta_val, p_val, p_fdr) {
    if (!sp %in% names(heatmap_data) || !trait_col %in% names(heatmap_data)) return(NULL)

    plot_df <- data.frame(
      site = heatmap_data$site,
      sp_count = heatmap_data[[sp]],
      trait_val = heatmap_data[[trait_col]],
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(trait_val))

    sp_display <- sp

    # Significance coloring: bold black if raw p < 0.05, grey otherwise
    sig_color <- ifelse(p_val < 0.05, "black", "grey40")
    sig_face <- ifelse(p_val < 0.05, "bold", "plain")

    ggplot(plot_df, aes(x = sp_count, y = trait_val)) +
      geom_point(aes(fill = site), shape = 21, color = "gray30",
                 stroke = 0.3, alpha = 0.5, size = 1.8) +
      geom_smooth(method = "lm", formula = y ~ sqrt(x),
                  color = "black", se = TRUE, alpha = 0.15, linewidth = 0.6) +
      scale_fill_manual(values = SITE_COLORS, guide = "none") +
      scale_x_sqrt(breaks = function(x) {
        mx <- max(x, na.rm = TRUE)
        if (mx <= 3) return(0:mx)
        unique(round(c(0, 1, ceiling(mx/3), ceiling(2*mx/3), mx)))
      }) +
      labs(
        title = bquote(italic(.(sp_display))),
        x = "Abundance",
        y = trait_ylab_map[[trait_label]]
      ) +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
               label = sprintf("\u03b2 = %+.2f, p = %.3f", beta_val, p_val),
               size = 2.8, color = sig_color, fontface = sig_face) +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 3.0,
               label = sprintf("n = %d corals", sum(plot_df$sp_count > 0)),
               size = 2.3, color = "grey50") +
      theme_multipanel(base_size = 9) +
      theme(
        plot.margin = margin(4, 6, 4, 4, "mm"),
        plot.title = element_text(size = 9, face = "italic", margin = margin(b = 2)),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.position = "none"
      )
  }

  # Create output directory for multi-page PDFs
  biplot_pdf_dir <- file.path(PATHS$fig_supplement, "species_biplots")
  dir.create(biplot_pdf_dir, showWarnings = FALSE, recursive = TRUE)

  # Also save to analysis directory
  biplot_pdf_dir2 <- file.path(fig_dir, "species_biplots")
  dir.create(biplot_pdf_dir2, showWarnings = FALSE, recursive = TRUE)

  panels_per_page <- 9
  page_w_mm <- 210  # A4-ish width
  page_h_mm <- 270  # A4-ish height

  for (trait_label in names(trait_col_map)) {
    trait_col <- trait_col_map[trait_label]
    trait_file <- trait_file_map[trait_label]

    cat("   Generating:", trait_label, "...\n")

    # Get beta/p for each species for this trait
    trait_stats <- heatmap_results %>%
      filter(as.character(trait) == trait_label) %>%
      mutate(species = as.character(species))

    # Build all panels for this trait
    all_panels <- list()
    for (sp in species_ordered) {
      row <- trait_stats %>% filter(species == sp)
      if (nrow(row) == 0) next

      p <- make_biplot(
        sp = sp,
        trait_label = trait_label,
        trait_col = trait_col,
        beta_val = round(as.numeric(row$beta), 2),
        p_val = as.numeric(row$p_raw),
        p_fdr = as.numeric(row$p_fdr)
      )
      if (!is.null(p)) all_panels[[length(all_panels) + 1]] <- p
    }

    n_panels <- length(all_panels)
    if (n_panels == 0) next

    n_pages <- ceiling(n_panels / panels_per_page)

    # Generate multi-page PDF
    pdf_path <- file.path(biplot_pdf_dir,
                          paste0("figS17_biplots_", trait_file, ".pdf"))
    pdf_path2 <- file.path(biplot_pdf_dir2,
                           paste0("figS17_biplots_", trait_file, ".pdf"))

    # Use pdf() for multi-page support (cairo_pdf is single-page only)
    pdf(pdf_path, width = page_w_mm / 25.4, height = page_h_mm / 25.4,
        onefile = TRUE)

    for (pg in 1:n_pages) {
      start_idx <- (pg - 1) * panels_per_page + 1
      end_idx <- min(pg * panels_per_page, n_panels)
      page_panels <- all_panels[start_idx:end_idx]

      # Pad with nulls if last page is incomplete
      while (length(page_panels) < panels_per_page) {
        page_panels[[length(page_panels) + 1]] <- plot_spacer()
      }

      page_plot <- wrap_plots(page_panels, ncol = 3) +
        plot_annotation(
          title = paste0(trait_label, " — Page ", pg, " of ", n_pages),
          tag_levels = list(LETTERS[start_idx:end_idx])
        ) &
        theme(
          plot.tag = element_text(size = 10, face = "bold"),
          legend.position = "none"
        )

      print(page_plot)
    }

    dev.off()

    # Copy to analysis dir
    file.copy(pdf_path, pdf_path2, overwrite = TRUE)

    cat("      Saved:", basename(pdf_path), "(", n_pages, "pages,",
        n_panels, "species)\n")
  }

  # Also generate a combined PNG of top hits (for quick reference / S17 figure)
  top_combos <- heatmap_results %>%
    filter(p_raw < 0.10 | abs(as.numeric(beta)) > 0.40) %>%
    arrange(p_raw) %>%
    head(9)

  top_panels <- list()
  for (i in 1:nrow(top_combos)) {
    sp <- as.character(top_combos$species[i])
    tl <- as.character(top_combos$trait[i])
    tc <- trait_col_map[tl]
    bv <- round(as.numeric(top_combos$beta[i]), 2)
    pv <- as.numeric(top_combos$p_raw[i])
    pf <- as.numeric(top_combos$p_fdr[i])

    p <- make_biplot(sp, tl, tc, bv, pv, pf)
    if (!is.null(p)) top_panels[[length(top_panels) + 1]] <- p
  }

  if (length(top_panels) > 0) {
    figS17_top <- wrap_plots(top_panels, ncol = 3) +
      plot_annotation(tag_levels = "A") &
      theme(
        legend.position = "none",
        plot.tag = element_text(size = 10, face = "bold")
      )

    save_figure(figS17_top,
               file.path(PATHS$fig_supplement, "figS17_species_trait_biplots.png"),
               width = 230, height = 225, units = "mm")
    save_figure(figS17_top,
               file.path(fig_dir, "figS17_species_trait_biplots.png"),
               width = 230, height = 225, units = "mm")
    cat("   Saved: figS17_species_trait_biplots.png (top 9 hits)\n")
  }
}

# ============================================================================
# FIGURE 5 LEGEND & RESULTS TEXT FILE
# ============================================================================

cat("Generating fig5_legend_results.txt...\n")

# Build CAFI → condition results lines (now with two-tier correction)
cafi_cond_lines <- paste0(
  "  ", cafi_to_condition_df$predictor,
  " [", cafi_to_condition_df$hypothesis_type, "]",
  ": beta = ", round(cafi_to_condition_df$estimate, 3),
  ", p_OLS = ", round(cafi_to_condition_df$p_value, 4),
  ", p_corrected = ", round(cafi_to_condition_df$p_corrected, 4),
  " (", ifelse(cafi_to_condition_df$hypothesis_type == "a_priori", "Hochberg",
        ifelse(cafi_to_condition_df$hypothesis_type == "supplement_composition", "uncorrected", "BH-FDR")), ")",
  ifelse(cafi_to_condition_df$p_corrected < 0.05, " *", ""),
  "  [HC3: p = ", round(cafi_to_condition_df$p_value_robust, 4), "]")

# Build condition → CAFI results lines
cond_cafi_lines <- paste0(
  "  ", condition_to_cafi_df$response,
  ": beta = ", round(condition_to_cafi_df$estimate, 3),
  ", p = ", round(condition_to_cafi_df$p_value, 4),
  ", p_FDR = ", round(condition_to_cafi_df$p_fdr, 4),
  ifelse(condition_to_cafi_df$p_fdr < 0.05, " *", ""))

fig5_legend <- paste0(
'FIGURE 5: CAFI DIVERSITY-CONDITION RELATIONSHIP (BEF Framework)
================================================================================

FIGURE LEGEND
-------------
Figure 5. CAFI diversity and abundance as predictors of coral physiological
condition. (A) Species richness versus coral condition score (PC1). Richness
is a significant positive predictor (Hochberg p = ', sprintf("%.3f", apriori_p_hochberg["Species richness"]), ', k = 2).
(B) Total CAFI abundance (sqrt-scaled x-axis) versus condition. Abundance is
a marginal positive predictor (Hochberg p = ', sprintf("%.3f", apriori_p_hochberg["Total CAFI"]), ', k = 2).
Points colored by site; lines show model fits with 95% CI shading. Richness
and abundance are strongly correlated (r = ', sprintf("%.2f", cor_raw_abund), '); variance partitioning
attributes 29.1% of explained variance uniquely to richness and <1% uniquely
to abundance (see text and Fig. S12).

All models: condition_PC1 ~ predictor + log(volume) + site (fixed effect).
OLS standard errors (primary; Breusch-Pagan confirms homoscedasticity).
HC3 robust SEs in supplement (conservative at n < 100; Long & Ervin 2000).
See Figure S12 for variance partitioning and path model diagnostics.
See Figure S13 for a priori forest plot, rarefied richness, exploratory
predictors, species scatter plots, and bidirectional tests.

================================================================================

STATISTICAL RESULTS
-------------------

1. A PRIORI BEF PREDICTORS (Hochberg FWER, k=2):
', paste(cafi_cond_lines[cafi_to_condition_df$hypothesis_type == "a_priori"], collapse = "\n"), '

2. EXPLORATORY PREDICTORS (BH-FDR, k=4):
', paste(cafi_cond_lines[cafi_to_condition_df$hypothesis_type == "exploratory"], collapse = "\n"), '

3. SUPPLEMENT COMPOSITION:
', paste(cafi_cond_lines[cafi_to_condition_df$hypothesis_type == "supplement_composition"], collapse = "\n"), '

4. REVERSE DIRECTION (Condition -> CAFI):
', paste(cond_cafi_lines, collapse = "\n"), '

5. RICHNESS ARTIFACT TEST:
   Raw richness: r(abundance) = ', sprintf("%.2f", cor_raw_abund), ', p_condition = ', sprintf("%.3f", p_raw_full), '
   Rarefied richness: r(abundance) = ', sprintf("%.2f", cor_rare_abund), ', p_condition = ', sprintf("%.3f", p_rare), '
   Note: Rarefaction is ambiguous -- may remove artifact OR the BEF mechanism

6. BEF VARIANCE PARTITIONING (see Figure S12):
   Partial regression (richness | abundance): beta = ', sprintf("%.4f", richness_partial_b), ', p = ', sprintf("%.4f", richness_partial_p), '
   Partial regression (abundance | richness): beta = ', sprintf("%.4f", abundance_partial_b), ', p = ', sprintf("%.4f", abundance_partial_p), '
   Variance unique to richness:  ', sprintf("%.4f (%.1f%%)", unique_richness, 100 * unique_richness / max(total_explained, 1e-10)), '
   Variance unique to abundance: ', sprintf("%.4f (%.1f%%)", unique_abundance, 100 * unique_abundance / max(total_explained, 1e-10)), '
   Variance shared:              ', sprintf("%.4f (%.1f%%)", shared, 100 * shared / max(total_explained, 1e-10)), '
   F-test richness unique:       p = ', sprintf("%.4f", p_rich_unique), '
   F-test abundance unique:      p = ', sprintf("%.4f", p_abund_unique), '
   Path model: richness -> condition beta = ', sprintf("%.2f", coef(path_cond)["richness_z"]), ',
               abundance -> condition beta = ', sprintf("%.2f", coef(path_cond)["abundance_z"]), '

================================================================================

RESULTS
-------

Both a priori BEF predictors showed positive effects on coral condition.
Species richness significantly predicted condition (Hochberg p = ', sprintf("%.3f", apriori_p_hochberg["Species richness"]), ',
k = 2). Total CAFI abundance was marginal (Hochberg p = ', sprintf("%.3f", apriori_p_hochberg["Total CAFI"]), ').
Variance partitioning (Figure S12) showed that richness accounted for
', sprintf("%.1f", 100 * unique_richness / max(total_explained, 1e-10)), '% of uniquely explained variance vs <1% for abundance,
and the path model confirmed this asymmetry (beta richness = ', sprintf("%.2f", coef(path_cond)["richness_z"]), ',
beta abundance = ', sprintf("%.2f", coef(path_cond)["abundance_z"]), ').

Rarefied richness showed no relationship with condition (p = ', sprintf("%.2f", p_rare), '),
but this test is ambiguous: rarefaction removes both sampling artifacts and
the diversity-mediated abundance pathway central to BEF theory.

No exploratory predictor (Shannon, Trapezia, Fish, Galeropsis) survived
BH-FDR correction. Community composition (PC1_CAFI) did not predict
condition (p = ', sprintf("%.3f", pc1_row$p_value), '; supplement). Reverse models (condition ->
CAFI) were non-significant.

================================================================================

COLOR SCHEME
------------
Significant (p < 0.05):      #2e7d32 (dark green)
Non-significant:             #757575 (medium gray)
A priori forest (Panel B):   green = significant, gray = non-significant
Exploratory forest (Panel D): gray (all NS after BH-FDR)

================================================================================
Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '
Source script: scripts/09_cafi_condition_feedbacks.R
')

writeLines(fig5_legend, file.path(PATHS$fig_manuscript, "fig5_legend_results.txt"))
cat("Saved: fig5_legend_results.txt\n\n")

# --- Additional exploratory figure: All CAFI metrics ---
p_cafi_effects <- ggplot(cafi_to_condition_df,
                          aes(x = reorder(predictor, estimate), y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, color = significant),
                  size = 0.8, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", p_value)),
            hjust = -0.3, size = 3) +
  scale_color_manual(values = c("TRUE" = "#2e7d32", "FALSE" = "#757575"),
                     name = "p < 0.05") +
  coord_flip() +
  labs(
    title = "CAFI Metrics as Predictors of Coral Condition",
    subtitle = "Linear models: Condition ~ CAFI + log(Volume) + Site (fixed effect)",
    x = "",
    y = "Effect size (regression coefficient)"
  ) +
  theme_publication()

save_figure(p_cafi_effects,
           file.path(fig_dir, "cafi_condition_effects.png"),
           width = 10, height = 6)
cat("   Saved: cafi_condition_effects.png\n")

# --- Forest plot: All functional group effects ---
p_functional_forest <- ggplot(functional_effects,
                               aes(x = reorder(predictor, estimate), y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper,
                       fill = matches_hypothesis),
                  color = "black", shape = 21, size = 1, stroke = 0.4, linewidth = 0.8) +
  geom_text(aes(label = paste0("p = ", format.pval(p_value, 2))),
            hjust = -0.2, vjust = 0.5, size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "#4CAF50", "FALSE" = "#FF5722"),
                    name = "Matches\nHypothesis") +
  coord_flip() +
  labs(
    title = "Functional Group Effects on Coral Condition",
    subtitle = "Green = effect direction matches ecological prediction",
    x = "",
    y = "Effect size (regression coefficient)",
    caption = paste0("n = ", n_complete, " corals | Models control for volume and site")
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

save_figure(p_functional_forest,
           file.path(fig_dir, "functional_effects_forest.png"),
           width = 9, height = 5)
cat("   Saved: functional_effects_forest.png\n\n")

# ============================================================================
# 4. SAVE STATISTICAL TABLES
# ============================================================================

cat("============================================================\n")
cat("SAVING STATISTICAL TABLES\n")
cat("============================================================\n\n")

# --- Table 1: CAFI -> Condition models ---
cafi_condition_table <- cafi_to_condition_df %>%
  mutate(
    across(where(is.numeric), ~round(., 4)),
    ci_95 = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
    significant = ifelse(significant, "Yes", "No")
  ) %>%
  dplyr::select(direction, predictor, estimate, se, t_value, df, p_value, ci_95,
         any_of(c("r2_marginal", "r2_adj")), n, significant)

save_table(cafi_condition_table, "cafi_condition_models")
cat("   Saved: cafi_condition_models.csv\n")

# --- Table 2: Condition -> CAFI models ---
condition_cafi_table <- condition_to_cafi_df %>%
  mutate(
    across(where(is.numeric), ~round(., 4)),
    ci_95 = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
    significant = ifelse(significant, "Yes", "No")
  ) %>%
  dplyr::select(direction, response, estimate, se, stat_type, stat_value,
         p_value, ci_95, n, significant)

save_table(condition_cafi_table, "reverse_direction_models")
cat("   Saved: reverse_direction_models.csv\n")

# --- Table 3: Functional group effects ---
functional_table <- functional_effects %>%
  mutate(
    across(where(is.numeric), ~round(., 4)),
    ci_95 = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
    significant = ifelse(significant, "Yes", "No"),
    matches_hypothesis = ifelse(matches_hypothesis, "Yes", "No")
  ) %>%
  dplyr::select(functional_group, predictor, estimate, se, t_value, p_value, ci_95,
         expected_sign, observed_sign, matches_hypothesis, significant)

save_table(functional_table, "functional_effects")
cat("   Saved: functional_effects.csv\n")

# --- Table 4: Richness-abundance artifact comparison ---
richness_comparison_df <- tryCatch(
  load_object("richness_comparison_results"),
  error = function(e) NULL
)

if (!is.null(richness_comparison_df)) {
  richness_table <- richness_comparison_df %>%
    mutate(
      across(where(is.numeric), ~round(., 4)),
      ci_95 = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
      significant = ifelse(significant, "Yes", "No")
    ) %>%
    dplyr::select(type, estimate, se, p_value, ci_95, n, cor_with_abundance, significant)

  save_table(richness_table, "richness_abundance_artifact")
  cat("   Saved: richness_abundance_artifact.csv\n")
}
cat("\n")

# ============================================================================
# 5. COMPILE STANDARDIZED STATISTICAL RESULTS
# ============================================================================

cat("============================================================\n")
cat("COMPILING STANDARDIZED RESULTS\n")
cat("============================================================\n\n")

# Add CAFI -> Condition results
for (i in 1:nrow(cafi_to_condition_df)) {
  row <- cafi_to_condition_df[i, ]
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H-CAFI-condition",
      question = paste("Does", row$predictor, "predict coral condition?"),
      test_name = "Linear model (fixed-effect site)",
      test_statistic = row$t_value,
      df = as.character(round(row$df, 1)),
      p_value = row$p_value,
      effect_size = row$estimate,
      effect_type = "Regression coefficient",
      n = row$n,
      notes = paste("95% CI: [", round(row$ci_lower, 3), ",", round(row$ci_upper, 3), "]")
    )
  )
}

# Add Condition -> CAFI results
for (i in 1:nrow(condition_to_cafi_df)) {
  row <- condition_to_cafi_df[i, ]
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H-reverse-direction",
      question = paste("Does condition predict", row$response, "?"),
      test_name = ifelse(row$stat_type == "z", "Negative binomial GLM", "Linear model"),
      test_statistic = row$stat_value,
      df = "NA",
      p_value = row$p_value,
      effect_size = row$estimate,
      effect_type = "Regression coefficient",
      n = row$n,
      notes = paste("95% CI: [", round(row$ci_lower, 3), ",", round(row$ci_upper, 3), "]")
    )
  )
}

# Save standardized results
save_stats_summary(stats_results, "09_cafi_condition_feedbacks",
                   "Bidirectional CAFI-Condition Feedback Analysis")

# ============================================================================
# 6. FINAL SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("ANALYSIS COMPLETE: CAFI-CORAL CONDITION FEEDBACKS\n")
cat("============================================================\n\n")

cat("SAMPLE SIZE: n =", n_complete, "corals with complete data\n\n")

cat("PART A: CAFI -> CONDITION EFFECTS\n")
cat("---------------------------------\n")
sig_cafi <- cafi_to_condition_df %>% filter(significant)
if (nrow(sig_cafi) > 0) {
  cat("Significant predictors of condition:\n")
  for (i in 1:nrow(sig_cafi)) {
    direction <- ifelse(sig_cafi$estimate[i] > 0, "positive", "negative")
    cat(sprintf("   - %s: %s effect (beta = %.3f, p = %.4f)\n",
                sig_cafi$predictor[i], direction, sig_cafi$estimate[i], sig_cafi$p_value[i]))
  }
} else {
  cat("   No significant CAFI predictors of condition (p < 0.05)\n")
}
cat("\n")

cat("PART B: CONDITION -> CAFI EFFECTS (REVERSE CAUSATION)\n")
cat("-----------------------------------------------------\n")
sig_reverse <- condition_to_cafi_df %>% filter(significant)
if (nrow(sig_reverse) > 0) {
  cat("Condition significantly predicts:\n")
  for (i in 1:nrow(sig_reverse)) {
    direction <- ifelse(sig_reverse$estimate[i] > 0, "positive", "negative")
    cat(sprintf("   - %s: %s effect (beta = %.3f, p = %.4f)\n",
                sig_reverse$response[i], direction, sig_reverse$estimate[i], sig_reverse$p_value[i]))
  }
} else {
  cat("   Condition does not significantly predict any CAFI metric\n")
  cat("   (Reverse direction not supported)\n")
}
cat("\n")

cat("PART C: FUNCTIONAL GROUP HYPOTHESES\n")
cat("------------------------------------\n")
for (i in 1:nrow(functional_effects)) {
  fg <- functional_effects[i, ]
  status <- ifelse(fg$significant, "SUPPORTED", "NOT SUPPORTED")
  match <- ifelse(fg$matches_hypothesis, "(direction matches)", "(direction opposite)")
  cat(sprintf("   %s: %s %s\n", fg$predictor, status,
              ifelse(fg$significant, match, "")))
}
cat("\n")

cat("PART F: KEY SPECIES FROM EXPERIMENTAL PAPER\n")
cat("--------------------------------------------\n")
if (exists("key_species_df") && nrow(key_species_df) > 0) {
  cat("Testing species with known effects from Stier et al.:\n")
  for (i in 1:nrow(key_species_df)) {
    sp <- key_species_df[i, ]
    match_text <- ifelse(sp$matches_experiment, "MATCHES", "CONTRADICTS")
    sig_text <- ifelse(sp$significant, "(SIGNIFICANT)", "(not significant)")
    cat(sprintf("   - %s: %s %s %s\n",
                sp$predictor, toupper(sp$observed_sign), match_text, sig_text))
  }
} else {
  cat("   Key species analysis not run or no results\n")
}
cat("\n")

cat("PART F2: SPECIES x INDIVIDUAL TRAIT CORRELATIONS\n")
cat("------------------------------------------------\n")
if (exists("species_trait_corr") && nrow(species_trait_corr) > 0) {
  cat("  Cross-study comparison table (cf. experimental Table S3):\n")
  cat("  -", nrow(species_trait_corr), "species x",
      length(available_labels), "traits\n")
  if ("r_PC1_Coral" %in% names(species_trait_corr)) {
    n_strong <- sum(abs(species_trait_corr$r_PC1_Coral) > 0.2, na.rm = TRUE)
    cat("  - Species with |r| > 0.2 for PC1_Coral:", n_strong, "\n")
  }
} else {
  cat("   Species-trait correlation analysis not run\n")
}
cat("\n")

cat("PART F3: INDIVIDUAL PHYSIOLOGY RESPONSES\n")
cat("-----------------------------------------\n")
if (exists("individual_physio_results") && nrow(individual_physio_results) > 0) {
  n_sig_physio <- sum(individual_physio_results$p_fdr < 0.05)
  cat("  Individual trait ~ CAFI predictor tests:\n")
  cat("  -", nrow(individual_physio_results), "tests run\n")
  cat("  - Significant after FDR:", n_sig_physio, "\n")
  if (n_sig_physio == 0) {
    cat("  - Conclusion: PC1_Coral null holds at individual trait level\n")
  }
} else {
  cat("   Individual physiology analysis not run\n")
}
cat("\n")

cat("PART G: NEIGHBORHOOD EFFECTS (Experimental Analog)\n")
cat("-------------------------------------------------\n")
if (exists("neighborhood_results")) {
  cat("  Analog to experimental manipulation of coral number:\n")
  cat("  - Sample:", n_neighborhood, "corals with neighborhood data\n")
  cat("  - Neighborhood effect on CAFI: β =", round(neighborhood_results$estimate[2], 3),
      "(p =", round(neighborhood_results$p_value[2], 4), ")\n")
  if (!is.na(neighborhood_results$estimate[5])) {
    cat("  - Neighborhood effect on condition: β =", round(neighborhood_results$estimate[5], 3),
        "(p =", round(neighborhood_results$p_value[5], 4), ")\n")
  }
  neighbor_interp <- ifelse(neighborhood_results$estimate[2] > 0, "Field of Dreams", "Redirection")
  cat("  - Interpretation:", neighbor_interp, "at landscape level\n")
} else {
  cat("   Neighborhood analysis not run\n")
}
cat("\n")

cat("PART H: LANDSCAPE-ONLY EFFECTS ON CONDITION\n")
cat("-------------------------------------------\n")
if (exists("landscape_effects_table") && nrow(landscape_effects_table) > 0) {
  cat("  Abiotic/spatial drivers of coral condition (no CAFI predictors):\n")
  cat("  - Sample:", n_landscape, "corals with condition data\n")

  # Volume effect
  vol_row <- landscape_effects_table %>% filter(response == "PC1_Coral", predictor == "log_volume")
  if (nrow(vol_row) > 0) {
    vol_sig <- ifelse(vol_row$p_value[1] < 0.05, "SIGNIFICANT", "not significant")
    cat("  - Volume effect: β =", round(vol_row$estimate[1], 3),
        "(p =", round(vol_row$p_value[1], 4), ") -", vol_sig, "\n")
  }

  # Neighborhood effect
  neigh_row <- landscape_effects_table %>% filter(response == "PC1_Coral", predictor == "n_neighbors")
  if (nrow(neigh_row) > 0) {
    neigh_sig <- ifelse(neigh_row$p_value[1] < 0.05, "SIGNIFICANT", "not significant")
    cat("  - Neighborhood effect: β =", round(neigh_row$estimate[1], 3),
        "(p =", round(neigh_row$p_value[1], 4), ") -", neigh_sig, "\n")
  }

  # Site effect
  site_row <- landscape_effects_table %>% filter(predictor == "site")
  if (nrow(site_row) > 0) {
    site_sig <- ifelse(site_row$p_value[1] < 0.05, "SIGNIFICANT", "not significant")
    cat("  - Site effect: R² =", round(site_row$adj_r2[1], 3),
        "(p =", round(site_row$p_value[1], 4), ") -", site_sig, "\n")
  }
} else {
  cat("   Landscape analysis not run\n")
}
cat("\n")

cat("INTERPRETATION:\n")
cat("---------------\n")
cat("1. CAFI effects on condition: Correlational evidence only\n")
cat("2. Experimental manipulation needed to establish causation\n")
cat("3. Reverse direction test: If condition predicts CAFI,\n")
cat("   then healthier corals may attract more associates\n")
cat("4. Bidirectional feedbacks possible but require longitudinal data\n")
cat("5. Key species: Survey results compared to experimental predictions\n\n")

cat("OUTPUT FILES:\n")
cat("  Figures:\n")
cat("    - output/figures/manuscript/fig5_feedbacks.png\n")
cat("    - output/figures/feedbacks/cafi_condition_effects.png\n")
cat("    - output/figures/feedbacks/functional_effects_forest.png\n")
cat("    - output/figures/feedbacks/key_species_effects.png\n")
cat("    - output/figures/feedbacks/functional_vs_key_species.png\n")
cat("    - output/figures/feedbacks/neighborhood_effects.png\n")
cat("    - output/figures/feedbacks/landscape_condition_effects.png\n")
cat("    - output/figures/supplement/figS10_rarefaction_sensitivity.png\n")
cat("    - output/figures/supplement/figS16_species_trait_heatmap.png\n")
cat("    - output/figures/supplement/figS17_species_trait_biplots.png\n")
cat("    - output/figures/feedbacks/diagnostics_richness_model.png\n")
cat("  Tables:\n")
cat("    - output/tables/cafi_condition_models.csv\n")
cat("    - output/tables/reverse_direction_models.csv\n")
cat("    - output/tables/functional_effects.csv\n")
cat("    - output/tables/landscape_condition_effects.csv\n")
cat("    - output/tables/landscape_model_comparison.csv\n")
cat("    - output/tables/site_condition_means.csv\n")
cat("    - output/tables/key_species_effects.csv\n")
cat("    - output/tables/species_trait_correlations.csv (cross-study comparison)\n")
cat("    - output/tables/individual_physiology_cafi_responses.csv\n")
cat("    - output/tables/cross_study_species_comparison.csv\n")
cat("    - output/tables/cross_study_sign_concordance.csv\n")
cat("    - output/tables/richness_abundance_artifact.csv\n")
cat("    - output/tables/neighborhood_effects.csv\n")
cat("    - output/tables/community_transform_sensitivity.csv\n")
cat("    - output/tables/species_trait_heatmap_data.csv\n")
cat("    - output/tables/09_cafi_condition_feedbacks_stats_summary.csv\n")
if (requireNamespace("mediation", quietly = TRUE) && n_complete >= 50) {
  cat("    - output/tables/path_analysis.csv (mediation analysis)\n")
}
cat("\n")

cat("Script 09 complete!\n\n")
