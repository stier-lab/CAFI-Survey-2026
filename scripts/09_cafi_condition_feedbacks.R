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
#   A. CAFI -> Condition Effects (mixed models)
#   B. Condition -> CAFI Effects (reverse direction test)
#   C. Functional group-specific models
#   D. Condition x Volume interaction
#   E. Path analysis (if sample size supports)
#
# OUTPUTS:
#   - output/figures/manuscript/fig5_feedbacks.png
#   - output/figures/cafi_condition_effects.png
#   - output/figures/functional_effects_forest.png
#   - output/tables/cafi_condition_models.csv
#   - output/tables/functional_effects.csv
#   - output/tables/path_analysis.csv (if implemented)
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
source(here::here("scripts/00_setup.R"))

# Load additional packages for mixed models and path analysis
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
  left_join(coral_master %>% dplyr::select(coral_id, volume, morphotype, depth_m, pc1_cafi, pc2_cafi),
            by = "coral_id") %>%
  filter(!is.na(condition_score), !is.na(total_cafi), !is.na(volume)) %>%
  mutate(
    log_volume = log10(volume),
    sqrt_volume = sqrt(volume),
    site = factor(site),
    # Standardize condition score for comparable effect sizes
    condition_z = scale(condition_score)[,1]
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
# Note: 3 sites is insufficient for random intercepts (Bolker et al. 2009 recommends >=5-6)
run_condition_model <- function(data, predictor_name, predictor_col) {
  # Build formula with fixed-effect site (3 levels insufficient for RE)
  formula_str <- paste("condition_score ~", predictor_col, "+ log_volume + site")

  tryCatch({
    # Fit linear model with site as fixed effect
    model <- lm(as.formula(formula_str), data = data)

    # Extract coefficient info
    model_summary <- summary(model)
    coef_table <- coef(model_summary)

    # Get the predictor coefficient
    pred_row <- which(rownames(coef_table) == predictor_col)
    if (length(pred_row) == 0) {
      pred_row <- 2  # Fallback to second row
    }

    estimate <- coef_table[pred_row, "Estimate"]
    se <- coef_table[pred_row, "Std. Error"]
    t_val <- coef_table[pred_row, "t value"]
    df_val <- model$df.residual
    p_val <- coef_table[pred_row, "Pr(>|t|)"]

    # Heteroscedasticity-robust SE (HC3 sandwich estimator)
    se_robust <- se  # default
    p_val_robust <- p_val
    if (requireNamespace("sandwich", quietly = TRUE) &&
        requireNamespace("lmtest", quietly = TRUE)) {
      vcov_hc3 <- sandwich::vcovHC(model, type = "HC3")
      robust_coefs <- lmtest::coeftest(model, vcov. = vcov_hc3)
      se_robust <- robust_coefs[pred_row, "Std. Error"]
      p_val_robust <- robust_coefs[pred_row, "Pr(>|t|)"]
    }

    # 95% CI (using robust SE)
    ci_lower <- estimate - qt(0.975, df_val) * se_robust
    ci_upper <- estimate + qt(0.975, df_val) * se_robust

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
      coef(model_std)[pred_row]
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
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      r2_adj = r2_adj,
      n = nrow(data),
      significant = p_val_robust < 0.05  # Use robust p-value
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
  c("PC1_CAFI (community)", "pc1_cafi"),  # Primary predictor - composite community score
  c("Total CAFI", "total_cafi"),
  c("Trapezia abundance", "n_trapezia"),
  c("Resident Fish abundance", "n_resident_fish"),
  c("Galeropsis abundance", "n_galeropsis"),
  c("Species richness", "otu_richness"),
  c("Shannon diversity", "shannon")
)

cafi_to_condition_results <- list()
cafi_to_condition_models <- list()

cat("Testing CAFI predictors of coral condition...\n")
cat("Model: Condition ~ CAFI_metric + log(Volume) + Site (fixed effect)\n\n")

for (pred in cafi_predictors) {
  cat("   Testing:", pred[1], "...\n")
  result <- run_condition_model(analysis_data, pred[1], pred[2])

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

# Apply FDR correction for multiple testing (7 CAFI metrics tested)
cafi_to_condition_df <- cafi_to_condition_df %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "BH"),
    significant_fdr = p_fdr < 0.05,
    significant = p_value < 0.05  # Keep unadjusted for reference
  )

cat("\nFDR-corrected significance (Benjamini-Hochberg):\n")
for (i in 1:nrow(cafi_to_condition_df)) {
  row <- cafi_to_condition_df[i, ]
  fdr_star <- ifelse(row$significant_fdr, " *FDR-SIG*", "")
  cat(sprintf("   %-30s p = %.4f, p_FDR = %.4f%s\n",
              row$predictor, row$p_value, row$p_fdr, fdr_star))
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
      model <- glm.nb(as.formula(formula_str), data = data)
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

    ci_lower <- estimate - 1.96 * se
    ci_upper <- estimate + 1.96 * se

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
condition_to_cafi_df <- condition_to_cafi_df %>%
  mutate(
    p_fdr = p.adjust(p_value, method = "BH"),
    significant_fdr = p_fdr < 0.05,
    significant = p_value < 0.05
  )

cat("\nFDR-corrected significance (Benjamini-Hochberg):\n")
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
  glm.nb(total_cafi ~ condition_score * log_volume + site, data = analysis_data)
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

if (LAVAAN_AVAILABLE && n_complete >= 50) {
  cat("Testing mediation: Volume -> CAFI -> Condition\n")
  cat("Sample size:", n_complete, "(minimum 50 recommended)\n\n")

  # Prepare standardized data
  path_data <- analysis_data %>%
    mutate(
      volume_z = scale(log_volume)[,1],
      cafi_z = scale(total_cafi)[,1],
      condition_z = scale(condition_score)[,1]
    ) %>%
    dplyr::select(volume_z, cafi_z, condition_z) %>%
    drop_na()

  # Specify path model
  path_model <- '
    # Direct effects
    cafi_z ~ a*volume_z        # Volume -> CAFI
    condition_z ~ b*cafi_z     # CAFI -> Condition
    condition_z ~ c*volume_z   # Volume -> Condition (direct)

    # Indirect effect (mediation)
    indirect := a*b

    # Total effect
    total := c + a*b

    # Proportion mediated
    prop_mediated := (a*b) / (c + a*b)
  '

  # Fit the model
  tryCatch({
    path_fit <- sem(path_model, data = path_data, se = "bootstrap", bootstrap = 1000)

    # Extract results
    path_results <- parameterEstimates(path_fit, standardized = TRUE,
                                        boot.ci.type = "perc")

    cat("Path Analysis Results:\n")
    cat("----------------------\n")

    # Key paths
    a_path <- path_results %>% filter(label == "a")
    b_path <- path_results %>% filter(label == "b")
    c_path <- path_results %>% filter(label == "c")
    indirect <- path_results %>% filter(label == "indirect")
    total <- path_results %>% filter(label == "total")
    prop_med <- path_results %>% filter(label == "prop_mediated")

    cat("\n   a: Volume -> CAFI:        beta =", round(a_path$est, 3),
        ", 95% CI [", round(a_path$ci.lower, 3), ",", round(a_path$ci.upper, 3), "]",
        ", p =", format.pval(a_path$pvalue, 3), "\n")
    cat("   b: CAFI -> Condition:     beta =", round(b_path$est, 3),
        ", 95% CI [", round(b_path$ci.lower, 3), ",", round(b_path$ci.upper, 3), "]",
        ", p =", format.pval(b_path$pvalue, 3), "\n")
    cat("   c: Volume -> Condition:   beta =", round(c_path$est, 3),
        "(direct)",
        ", 95% CI [", round(c_path$ci.lower, 3), ",", round(c_path$ci.upper, 3), "]",
        ", p =", format.pval(c_path$pvalue, 3), "\n")
    cat("\n   Indirect (a*b):           beta =", round(indirect$est, 3),
        ", 95% CI [", round(indirect$ci.lower, 3), ",", round(indirect$ci.upper, 3), "]",
        ", p =", format.pval(indirect$pvalue, 3), "\n")
    cat("   Total effect:             beta =", round(total$est, 3), "\n")

    if (!is.na(prop_med$est)) {
      cat("   Proportion mediated:     ", round(prop_med$est * 100, 1), "%\n")
    }

    # Model fit
    # NOTE: This path model is just-identified (saturated, df=0) with 3 observed variables
    # and all possible paths. Fit indices are trivially perfect and uninformative.
    model_df <- fitMeasures(path_fit, "df")
    fit_measures <- fitMeasures(path_fit, c("cfi", "rmsea", "srmr"))
    cat("\n   Model df =", model_df,
        ifelse(model_df == 0, " [SATURATED - fit indices uninformative]", ""), "\n")
    if (model_df > 0) {
      cat("   Model Fit: CFI =", round(fit_measures["cfi"], 3),
          ", RMSEA =", round(fit_measures["rmsea"], 3),
          ", SRMR =", round(fit_measures["srmr"], 3), "\n\n")
    } else {
      cat("   (CFI=1, RMSEA=0, SRMR=0 by definition for saturated model)\n")
      cat("   Path coefficients and bootstrap CIs remain valid.\n\n")
    }

    # Save path analysis results
    path_results_clean <- path_results %>%
      dplyr::select(label, lhs, op, rhs, est, se, ci.lower, ci.upper, pvalue) %>%
      rename(
        path = label,
        from = lhs,
        operator = op,
        to = rhs,
        estimate = est,
        p_value = pvalue
      )

    save_table(path_results_clean, "path_analysis")

    # Mediation interpretation
    if (indirect$pvalue < 0.05) {
      cat("INTERPRETATION: Significant mediation detected.\n")
      cat("   CAFI partially mediates the Volume -> Condition relationship.\n\n")
    } else {
      cat("INTERPRETATION: No significant mediation.\n")
      cat("   The Volume -> Condition relationship is not significantly mediated by CAFI.\n\n")
    }

  }, error = function(e) {
    cat("Path analysis failed:", e$message, "\n\n")
  })

} else if (!LAVAAN_AVAILABLE) {
  cat("Path analysis skipped: lavaan package not available\n")
  cat("Install with: install.packages('lavaan')\n\n")
} else {
  cat("Path analysis skipped: Sample size (n =", n_complete, ") too small\n")
  cat("Minimum recommended: n = 50 for reliable SEM\n\n")
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

    # NEGATIVE effect species (from experiment)
    n_cymo = sum(grepl("Cymo", otu, ignore.case = TRUE), na.rm = TRUE),
    n_luniella = sum(grepl("Luniella", otu, ignore.case = TRUE), na.rm = TRUE),
    n_xanthid_harmful = sum(grepl("Cymo|Luniella", otu, ignore.case = TRUE), na.rm = TRUE),

    .groups = "drop"
  )

# Merge with analysis data
analysis_data <- analysis_data %>%
  left_join(key_species_counts, by = "coral_id")

# Replace NA with 0 for species counts
key_species_cols <- c("n_caracanthus", "n_alpheus_lottini", "n_alpheus_all",
                      "n_cymo", "n_luniella", "n_xanthid_harmful")
for (col in key_species_cols) {
  analysis_data[[col]] <- ifelse(is.na(analysis_data[[col]]), 0, analysis_data[[col]])
}

# Summary of key species presence
cat("   Key Species Presence in Condition Dataset (n =", nrow(analysis_data), "corals):\n")
cat("   POSITIVE effect species (predicted from experiment):\n")
cat("      Caracanthus maculatus:", sum(analysis_data$n_caracanthus > 0), "corals,",
    sum(analysis_data$n_caracanthus), "individuals\n")
cat("      Alpheus lottini:", sum(analysis_data$n_alpheus_lottini > 0), "corals,",
    sum(analysis_data$n_alpheus_lottini), "individuals\n")
cat("      All Alpheus spp.:", sum(analysis_data$n_alpheus_all > 0), "corals,",
    sum(analysis_data$n_alpheus_all), "individuals\n")
cat("   NEGATIVE effect species (predicted from experiment):\n")
cat("      Cymo quadrilobatus:", sum(analysis_data$n_cymo > 0), "corals,",
    sum(analysis_data$n_cymo), "individuals\n")
cat("      Luniella pugil:", sum(analysis_data$n_luniella > 0), "corals,",
    sum(analysis_data$n_luniella), "individuals\n")
cat("      Combined xanthids:", sum(analysis_data$n_xanthid_harmful > 0), "corals,",
    sum(analysis_data$n_xanthid_harmful), "individuals\n\n")

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
  list(name = "Cymo quadrilobatus", col = "n_cymo", expected = "negative",
       note = "Xanthid crab - experimental negative effect"),
  list(name = "Luniella pugil", col = "n_luniella", expected = "negative",
       note = "Xanthid crab - experimental negative effect"),
  list(name = "Harmful xanthids (combined)", col = "n_xanthid_harmful", expected = "negative",
       note = "Cymo + Luniella combined")
)

key_species_results <- list()

for (sp in key_species_predictors) {
  # Check if there's enough variation to fit model
  n_present <- sum(analysis_data[[sp$col]] > 0)

  if (n_present >= 5) {  # Need at least 5 corals with species present
    result <- run_condition_model(analysis_data, sp$name, sp$col)

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

# Apply FDR correction across key species tests (6 tests)
if (nrow(key_species_df) > 0) {
  key_species_df <- key_species_df %>%
    mutate(
      p_fdr = p.adjust(p_value, method = "BH"),
      significant_fdr = p_fdr < 0.05,
      significant = p_value < 0.05  # Keep unadjusted for reference
    )

  cat("\nFDR-corrected key species p-values (Benjamini-Hochberg):\n")
  for (i in 1:nrow(key_species_df)) {
    row <- key_species_df[i, ]
    fdr_star <- ifelse(row$significant_fdr, " *FDR-SIG*", "")
    cat(sprintf("   %-28s p = %.4f, p_FDR = %.4f%s\n",
                row$predictor, row$p_value, row$p_fdr, fdr_star))
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
cat("F.4 Creating Key Species Effects Figure...\n")

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
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, fill = color_var),
                    color = "black", shape = 21, size = 1.2, linewidth = 1) +
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
      caption = paste0("n = ", nrow(analysis_data), " corals | Vertical line at 0 = no effect\n",
                       "Green = matches experimental prediction | Red = contradicts")
    ) +
    theme_publication() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 11),
      plot.caption = element_text(hjust = 0, size = 9)
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.4)))  # Extra space for labels

  ggsave(file.path(fig_dir, "key_species_effects.png"),
         p_key_species_forest, width = 10, height = 8, dpi = 300, bg = "white")
  cat("   Saved: key_species_effects.png\n\n")

  # F.6 Combined comparison: Functional Groups vs Key Species
  cat("F.5 Creating Combined Effects Comparison Figure...\n")

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
    "Sig. matches prediction" = "#2e7d32",
    "Sig. contradicts prediction" = "#c62828",
    "NS matches prediction" = "#81c784",
    "NS contradicts prediction" = "#ef9a9a"
  )

  p_combined <- ggplot(combined_comparison,
                        aes(x = reorder(predictor, estimate), y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, fill = color_var),
                    color = "black", shape = 21, size = 1, linewidth = 0.8) +
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

  ggsave(file.path(fig_dir, "functional_vs_key_species.png"),
         p_combined, width = 11, height = 9, dpi = 300, bg = "white")
  cat("   Saved: functional_vs_key_species.png\n\n")
}

# F.7 Save key species results table
save_table(key_species_df %>%
             mutate(across(where(is.numeric), ~round(., 4))) %>%
             dplyr::select(predictor, estimate, se, t_value, p_value, ci_lower, ci_upper,
                    expected_sign, observed_sign, matches_experiment, n_present, note),
           "key_species_effects")
cat("   Saved: key_species_effects.csv\n\n")

# F.8 Interpretation
cat("F.6 Key Species Analysis Interpretation:\n")
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

cat("############################################################\n")
cat("PART G: NEIGHBORHOOD EFFECTS (Analog to Experimental Treatment)\n")
cat("############################################################\n\n")

cat("Testing whether neighborhood density (analogous to experimental manipulation\n")
cat("of coral number) affects CAFI abundance and coral condition.\n\n")

# G.1 Prepare data with neighborhood variables
# Only subset of corals have neighborhood data (from 5m surveys)
neighborhood_data <- coral_master %>%
  filter(!is.na(n_neighbors), !is.na(volume_field)) %>%
  mutate(
    log_volume = log10(volume_field),
    log_volume_scaled = scale(log_volume)[,1],
    n_neighbors_scaled = scale(n_neighbors)[,1],
    total_neighbor_vol_scaled = scale(log10(total_neighbor_volume + 1))[,1],
    # Size category for visualization
    size_category = cut(volume_field,
                        breaks = quantile(volume_field, probs = c(0, 0.33, 0.67, 1)),
                        labels = c("Small", "Medium", "Large"),
                        include.lowest = TRUE)
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

m_cafi_neighborhood <- MASS::glm.nb(
  total_cafi ~ log_volume_scaled * n_neighbors_scaled + site,
  data = neighborhood_data
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
  geom_point(aes(color = size_category, size = volume_field), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Small" = "#2196F3", "Medium" = "#FF9800", "Large" = "#E91E63"),
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
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", color = NA))

# Panel B: Condition vs neighborhood density
if (!is.null(m_condition_neighborhood)) {
  p_neighborhood_condition <- ggplot(neighborhood_condition,
                                      aes(x = n_neighbors, y = condition_score)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_point(aes(color = site), alpha = 0.7, size = 3) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    labs(
      title = "B. Neighborhood Density → Coral Condition",
      subtitle = paste0("Position-corrected condition score (n = ", n_with_condition, ")"),
      x = "Number of neighboring Pocillopora (within 5m)",
      y = "Coral condition (PC1_Coral)"
    ) +
    theme_publication() +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
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
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
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
    caption = "From model: Condition ~ Volume + Neighbors + PC1_CAFI + (1|Site)"
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

ggsave(file.path(fig_dir, "neighborhood_effects.png"),
       fig_neighborhood, width = 12, height = 10, dpi = 300, bg = "white")
cat("   Saved: neighborhood_effects.png\n")

# G.7 Save neighborhood results table
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
cat("G.7 Neighborhood Effects Interpretation:\n")
cat("    =====================================================================\n")
cat("    Connection to experimental paper (Stier et al.):\n")
cat("    - Experiment: Manipulated coral NUMBER to test landscape effects\n")
cat("    - This survey: Uses natural variation in n_neighbors (within 5m)\n")
cat("    - Hypothesis: Higher neighborhood density → more CAFI (Field of Dreams)\n")
cat("                  OR → fewer CAFI per coral (Redistribution)\n\n")

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
                          "More neighbors → fewer CAFI per focal coral (consistent with Redistribution)")
cat("\n    Interpretation:", neighbor_interp, "\n")
cat("    =====================================================================\n\n")

# ============================================================================
# 3. CREATE VISUALIZATIONS
# ============================================================================

cat("============================================================\n")
cat("CREATING VISUALIZATIONS\n")
cat("============================================================\n\n")

# --- Panel A: PC1_CAFI (primary metric) vs Condition ---
# This is the key bidirectional test following Stier et al. (2024)
p_pc1_cafi <- ggplot(analysis_data %>% filter(!is.na(pc1_cafi)),
                      aes(x = pc1_cafi, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_manual(values = SITE_COLORS, name = "Site") +
  labs(
    title = "A. PC1_CAFI → Condition",
    subtitle = "Community composition predicts coral health",
    x = "PC1_CAFI (community score)",
    y = "Coral condition (PC1_Coral)"
  ) +
  theme_publication() +
  theme(legend.position = c(0.15, 0.85))

# --- Panel Aa: Trapezia (defenders) vs Condition ---
p_trapezia <- ggplot(analysis_data, aes(x = n_trapezia, y = condition_score)) +
  geom_jitter(aes(color = site), alpha = 0.6, width = 0.15, size = 2.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_manual(values = SITE_COLORS, name = "Site") +
  labs(
    title = "B. Trapezia (Defenders)",
    subtitle = "Expected: positive (predator defense)",
    x = "Trapezia abundance",
    y = "Coral condition score"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# --- Panel C: Galeropsis vs Condition ---
p_galeropsis <- ggplot(analysis_data, aes(x = n_galeropsis, y = condition_score)) +
  geom_jitter(aes(color = site), alpha = 0.6, width = 0.15, size = 2.5) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_manual(values = SITE_COLORS, name = "Site") +
  labs(
    title = "C. Galeropsis (Tissue Consumer)",
    subtitle = "Expected: negative (Coralliophilinae)",
    x = "Galeropsis abundance",
    y = "Coral condition score"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# --- Panel C: Forest plot of functional group effects ---
forest_data <- functional_effects %>%
  mutate(
    functional_group = factor(functional_group,
                               levels = c("Trapezia\n(Crabs)",
                                          "Fish",
                                          "Galeropsis"))
  )

# Add expected direction indicator
forest_data$expected_color <- ifelse(forest_data$expected_sign == "positive",
                                      "#009E73", "#D55E00")

p_forest <- ggplot(forest_data, aes(x = functional_group, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, color = significant),
                  size = 1.2, linewidth = 1.2) +
  geom_text(aes(label = sprintf("p = %.3f", p_value)),
            vjust = -1.5, size = 3.5) +
  scale_color_manual(values = c("TRUE" = "#2e7d32", "FALSE" = "#757575"),
                     name = "Significant\n(p < 0.05)",
                     labels = c("TRUE" = "Yes", "FALSE" = "No")) +
  coord_flip() +
  labs(
    title = "D. Functional Group Effects on Condition",
    subtitle = "Effect sizes with 95% CI (mixed model estimates)",
    x = "",
    y = "Effect on condition score (standardized)"
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

# --- Panel D: Bidirectionality comparison ---
# Combine both directions for visualization
bidirectional_data <- bind_rows(
  cafi_to_condition_df %>%
    filter(predictor == "Total CAFI") %>%
    mutate(direction = "CAFI -> Condition",
           label = "CAFI predicts\ncondition"),
  condition_to_cafi_df %>%
    filter(response == "Total CAFI") %>%
    rename(predictor = response) %>%
    mutate(direction = "Condition -> CAFI",
           label = "Condition predicts\nCAFI")
)

p_bidirectional <- ggplot(bidirectional_data,
                           aes(x = label, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, color = significant),
                  size = 1.2, linewidth = 1.2) +
  geom_text(aes(label = sprintf("p = %.3f", p_value)),
            vjust = -1.5, size = 3.5) +
  scale_color_manual(values = c("TRUE" = "#2e7d32", "FALSE" = "#757575"),
                     name = "Significant") +
  coord_flip() +
  labs(
    title = "E. Bidirectional Effects",
    subtitle = "Testing both directions of causation",
    x = "",
    y = "Standardized effect size"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# --- Combine into manuscript figure ---
# Following Stier et al. (2024) layout with PC1_CAFI as primary metric
fig5_feedbacks <- (p_pc1_cafi | p_trapezia | p_galeropsis) /
                   (p_forest | p_bidirectional) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Figure 5: CAFI-Coral Condition Feedbacks",
    subtitle = paste0("n = ", n_complete, " corals | PC1_CAFI = PCA-based community composition score"),
    caption = "PC1_Coral (condition) = position-corrected PC1 of protein, carbs, zoox, AFDW\nModels control for coral volume and site | Methodology follows Stier et al. (2024)"
  ) &
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11))

# Save manuscript figure
ggsave(file.path(PATHS$fig_manuscript, "fig5_feedbacks.png"),
       fig5_feedbacks, width = 12, height = 10, dpi = 300, bg = "white")
cat("   Saved: fig5_feedbacks.png (manuscript)\n")

# --- Additional exploratory figure: All CAFI metrics ---
p_cafi_effects <- ggplot(cafi_to_condition_df,
                          aes(x = reorder(predictor, estimate), y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper, color = significant),
                  size = 0.8, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", p_value)),
            hjust = -0.3, size = 3) +
  scale_color_manual(values = c("TRUE" = "#2e7d32", "FALSE" = "#757575"),
                     name = "p < 0.05") +
  coord_flip() +
  labs(
    title = "CAFI Metrics as Predictors of Coral Condition",
    subtitle = "Mixed effects models: Condition ~ CAFI + log(Volume) + (1|Site)",
    x = "",
    y = "Effect size (regression coefficient)"
  ) +
  theme_publication()

ggsave(file.path(fig_dir, "cafi_condition_effects.png"),
       p_cafi_effects, width = 10, height = 6, dpi = 300, bg = "white")
cat("   Saved: cafi_condition_effects.png\n")

# --- Forest plot: All functional group effects ---
p_functional_forest <- ggplot(functional_effects,
                               aes(x = reorder(predictor, estimate), y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper,
                       fill = matches_hypothesis),
                  color = "black", shape = 21, size = 1, linewidth = 0.8) +
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

ggsave(file.path(fig_dir, "functional_effects_forest.png"),
       p_functional_forest, width = 9, height = 5, dpi = 300, bg = "white")
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
         r2_marginal, n, significant)

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
cat("   Saved: functional_effects.csv\n\n")

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
  cat("   (Reverse causation not supported)\n")
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
  neighbor_interp <- ifelse(neighborhood_results$estimate[2] > 0, "Field of Dreams", "Redistribution")
  cat("  - Interpretation:", neighbor_interp, "at landscape level\n")
} else {
  cat("   Neighborhood analysis not run\n")
}
cat("\n")

cat("INTERPRETATION:\n")
cat("---------------\n")
cat("1. CAFI effects on condition: Correlational evidence only\n")
cat("2. Experimental manipulation needed to establish causation\n")
cat("3. Reverse causation test: If condition predicts CAFI,\n")
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
cat("    - output/figures/feedbacks/neighborhood_effects.png (NEW)\n")
cat("  Tables:\n")
cat("    - output/tables/cafi_condition_models.csv\n")
cat("    - output/tables/reverse_direction_models.csv\n")
cat("    - output/tables/functional_effects.csv\n")
cat("    - output/tables/key_species_effects.csv\n")
cat("    - output/tables/neighborhood_effects.csv (NEW)\n")
if (LAVAAN_AVAILABLE && n_complete >= 50) {
  cat("    - output/tables/path_analysis.csv\n")
}
cat("\n")

cat("Script 09 complete!\n\n")
