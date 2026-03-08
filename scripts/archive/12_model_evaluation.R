# ============================================================================
# 12_model_evaluation.R - Rigorous Model Evaluation and Diagnostics
# ============================================================================
#
# PURPOSE: Comprehensive evaluation and diagnostics for machine learning models
#
# RUNS AFTER: 11_machine_learning.R
#
# EXPECTED INPUTS (from 11_machine_learning.R):
#   - output/objects/machine_learning_models.rds  : Main ML results
#   - output/objects/rf_model_abundance.rds       : Random Forest abundance model
#   - output/objects/rf_model_richness.rds        : Random Forest richness model
#   - output/objects/xgb_model_abundance.rds      : XGBoost abundance model
#
# EVALUATION COMPONENTS:
#   A. Cross-Validation Analysis
#   B. Prediction Accuracy (Observed vs Predicted)
#   C. Residual Diagnostics
#   D. Site-Specific Performance
#   E. Ecological Validity
#   F. Uncertainty Quantification
#   G. Learning Curves
#   H. Model Selection Recommendations
#
# OUTPUTS:
#   Figures:
#     - output/figures/observed_vs_predicted.png
#     - output/figures/residual_diagnostics.png
#     - output/figures/learning_curves.png
#     - output/figures/site_performance.png
#   Tables:
#     - output/tables/cv_results_detailed.csv
#     - output/tables/site_metrics.csv
#     - output/tables/model_recommendations.csv
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-21
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    MODEL EVALUATION AND DIAGNOSTICS\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

# Load setup
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))

# Additional packages for model evaluation
suppressPackageStartupMessages({
  if (!require("caret", quietly = TRUE)) {
    cat("Installing caret package...\n")
    install.packages("caret", quiet = TRUE)
    library(caret)
  } else {
    library(caret)
  }

  if (!require("ranger", quietly = TRUE)) {
    cat("Note: ranger package not installed - some RF diagnostics may be limited\n")
  }
})

# Create evaluation figures directory
FIG_DIR <- file.path(PATHS$figures, "12_evaluation")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Calculate comprehensive regression metrics
#' @param observed Vector of observed values
#' @param predicted Vector of predicted values
#' @return Named list of metrics
calc_metrics <- function(observed, predicted) {
  # Remove NAs
  valid <- complete.cases(observed, predicted)
  obs <- observed[valid]
  pred <- predicted[valid]

  if (length(obs) < 3) {
    return(list(
      n = length(obs),
      r_squared = NA_real_,
      rmse = NA_real_,
      mae = NA_real_,
      mape = NA_real_,
      bias = NA_real_,
      cor = NA_real_
    ))
  }

  # Calculate metrics
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r_squared <- 1 - (ss_res / ss_tot)
  rmse <- sqrt(mean((obs - pred)^2))
  mae <- mean(abs(obs - pred))

  # MAPE (avoid division by zero)
  mape <- if (any(obs != 0)) {
    mean(abs((obs - pred) / obs)[obs != 0]) * 100
  } else NA_real_

  # Bias (mean prediction error)
  bias <- mean(pred - obs)

  # Correlation
  cor_val <- cor(obs, pred)

  list(
    n = length(obs),
    r_squared = r_squared,
    rmse = rmse,
    mae = mae,
    mape = mape,
    bias = bias,
    cor = cor_val
  )
}

#' Create observed vs predicted plot with diagnostics
#' @param observed Vector of observed values
#' @param predicted Vector of predicted values
#' @param title Plot title
#' @param site Optional site vector for coloring
#' @return ggplot object
plot_obs_vs_pred <- function(observed, predicted, title = "Observed vs Predicted",
                              site = NULL, model_name = "") {
  df <- data.frame(
    observed = observed,
    predicted = predicted
  )

  if (!is.null(site)) {
    df$site <- site
  }

  # Calculate metrics for subtitle
  metrics <- calc_metrics(observed, predicted)

  subtitle <- sprintf("R² = %.3f, RMSE = %.2f, MAE = %.2f",
                      metrics$r_squared, metrics$rmse, metrics$mae)

  # Base plot
  p <- ggplot(df, aes(x = observed, y = predicted))

  if (!is.null(site)) {
    p <- p + geom_point(aes(color = site), alpha = 0.7, size = 2.5) +
      scale_color_manual(values = SITE_COLORS, name = "Site")
  } else {
    p <- p + geom_point(alpha = 0.7, size = 2.5, color = "#4477AA")
  }

  # Add 1:1 line and smoothed trend
  p <- p +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "black", linewidth = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#CC6677",
                linewidth = 0.6, alpha = 0.2) +
    labs(
      x = "Observed",
      y = "Predicted",
      title = title,
      subtitle = subtitle
    ) +
    coord_equal() +
    theme(legend.position = if(!is.null(site)) "right" else "none")

  p
}

#' Create residual diagnostic plots
#' @param observed Vector of observed values
#' @param predicted Vector of predicted values
#' @param predictors Optional data frame of predictors
#' @return patchwork object with multiple diagnostic plots
plot_residual_diagnostics <- function(observed, predicted, predictors = NULL,
                                       model_name = "") {
  residuals <- observed - predicted

  df <- data.frame(
    observed = observed,
    predicted = predicted,
    residuals = residuals,
    std_residuals = scale(residuals)[, 1]
  )

  # 1. Residuals vs Fitted
  p1 <- ggplot(df, aes(x = predicted, y = residuals)) +
    geom_point(alpha = 0.6, size = 2, color = "#4477AA") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = TRUE, color = "#CC6677",
                linewidth = 0.6, alpha = 0.2) +
    labs(
      x = "Fitted Values",
      y = "Residuals",
      title = "A. Residuals vs Fitted"
    )

  # 2. Q-Q Plot for normality
  p2 <- ggplot(df, aes(sample = std_residuals)) +
    stat_qq(alpha = 0.6, color = "#4477AA") +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(
      x = "Theoretical Quantiles",
      y = "Standardized Residuals",
      title = "B. Q-Q Plot (Normality)"
    )

  # 3. Scale-Location (for heteroscedasticity)
  df$sqrt_abs_std_resid <- sqrt(abs(df$std_residuals))

  p3 <- ggplot(df, aes(x = predicted, y = sqrt_abs_std_resid)) +
    geom_point(alpha = 0.6, size = 2, color = "#4477AA") +
    geom_smooth(method = "loess", se = TRUE, color = "#CC6677",
                linewidth = 0.6, alpha = 0.2) +
    labs(
      x = "Fitted Values",
      y = expression(sqrt("|Standardized Residuals|")),
      title = "C. Scale-Location (Heteroscedasticity)"
    )

  # 4. Histogram of residuals
  p4 <- ggplot(df, aes(x = residuals)) +
    geom_histogram(bins = 20, fill = "#4477AA", color = "white", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x = "Residuals",
      y = "Count",
      title = "D. Residual Distribution"
    )

  # Combine plots
  combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = paste("Residual Diagnostics:", model_name),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )

  combined
}

#' Calculate site-specific metrics
#' @param observed Vector of observed values
#' @param predicted Vector of predicted values
#' @param site Vector of site labels
#' @return Data frame with metrics by site
calc_site_metrics <- function(observed, predicted, site) {
  df <- data.frame(
    observed = observed,
    predicted = predicted,
    site = site
  )

  site_metrics <- df %>%
    group_by(site) %>%
    summarise(
      n = n(),
      mean_observed = mean(observed, na.rm = TRUE),
      mean_predicted = mean(predicted, na.rm = TRUE),
      r_squared = {
        m <- calc_metrics(observed, predicted)
        m$r_squared
      },
      rmse = {
        m <- calc_metrics(observed, predicted)
        m$rmse
      },
      mae = {
        m <- calc_metrics(observed, predicted)
        m$mae
      },
      bias = {
        m <- calc_metrics(observed, predicted)
        m$bias
      },
      .groups = "drop"
    )

  site_metrics
}

#' Generate learning curve data
#' @param model_fun Function that fits model and returns predictions
#' @param data Training data
#' @param response Response variable name
#' @param fractions Vector of training set fractions to use
#' @param n_iter Number of iterations per fraction
#' @return Data frame with learning curve data
generate_learning_curve <- function(data, response_col, predictor_cols,
                                    fractions = seq(0.1, 1, by = 0.1),
                                    n_iter = 5, model_type = "lm") {
  n <- nrow(data)
  results <- list()

  for (frac in fractions) {
    train_size <- round(n * frac)

    for (iter in 1:n_iter) {
      # Sample training indices
      train_idx <- sample(1:n, train_size, replace = FALSE)
      train_data <- data[train_idx, ]
      test_data <- data[-train_idx, ]

      # Skip if test set is too small
      if (nrow(test_data) < 5) {
        test_data <- train_data  # Use training data for evaluation
      }

      # Fit model based on type
      tryCatch({
        formula <- as.formula(paste(response_col, "~",
                                    paste(predictor_cols, collapse = " + ")))

        if (model_type == "lm") {
          model <- lm(formula, data = train_data)
        } else if (model_type == "glm.nb") {
          model <- MASS::glm.nb(formula, data = train_data)
        } else if (model_type == "rf" && require("ranger", quietly = TRUE)) {
          model <- ranger::ranger(formula, data = train_data, num.trees = 100)
        } else {
          model <- lm(formula, data = train_data)
        }

        # Predictions
        if (model_type == "rf" && require("ranger", quietly = TRUE)) {
          train_pred <- predict(model, data = train_data)$predictions
          test_pred <- predict(model, data = test_data)$predictions
        } else {
          train_pred <- predict(model, newdata = train_data, type = "response")
          test_pred <- predict(model, newdata = test_data, type = "response")
        }

        # Metrics
        train_metrics <- calc_metrics(train_data[[response_col]], train_pred)
        test_metrics <- calc_metrics(test_data[[response_col]], test_pred)

        results[[length(results) + 1]] <- data.frame(
          fraction = frac,
          train_size = train_size,
          iteration = iter,
          train_rmse = train_metrics$rmse,
          test_rmse = test_metrics$rmse,
          train_r2 = train_metrics$r_squared,
          test_r2 = test_metrics$r_squared
        )

      }, error = function(e) {
        # Skip failed iterations
        NULL
      })
    }
  }

  if (length(results) > 0) {
    bind_rows(results)
  } else {
    data.frame()
  }
}

#' Plot learning curves
#' @param learning_data Data frame from generate_learning_curve
#' @param title Plot title
#' @return ggplot object
plot_learning_curves <- function(learning_data, title = "Learning Curves") {
  if (nrow(learning_data) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "Insufficient data for learning curves") +
             theme_void())
  }

  # Summarize across iterations
  summary_data <- learning_data %>%
    group_by(fraction, train_size) %>%
    summarise(
      train_rmse_mean = mean(train_rmse, na.rm = TRUE),
      train_rmse_se = sd(train_rmse, na.rm = TRUE) / sqrt(n()),
      test_rmse_mean = mean(test_rmse, na.rm = TRUE),
      test_rmse_se = sd(test_rmse, na.rm = TRUE) / sqrt(n()),
      train_r2_mean = mean(train_r2, na.rm = TRUE),
      test_r2_mean = mean(test_r2, na.rm = TRUE),
      .groups = "drop"
    )

  # RMSE plot
  p_rmse <- summary_data %>%
    pivot_longer(
      cols = c(train_rmse_mean, test_rmse_mean),
      names_to = "set",
      values_to = "rmse"
    ) %>%
    mutate(set = ifelse(set == "train_rmse_mean", "Training", "Validation")) %>%
    ggplot(aes(x = train_size, y = rmse, color = set)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Training" = "#4477AA", "Validation" = "#CC6677")) +
    labs(
      x = "Training Set Size",
      y = "RMSE",
      title = "A. RMSE vs Training Size",
      color = "Dataset"
    )

  # R² plot
  p_r2 <- summary_data %>%
    pivot_longer(
      cols = c(train_r2_mean, test_r2_mean),
      names_to = "set",
      values_to = "r2"
    ) %>%
    mutate(set = ifelse(set == "train_r2_mean", "Training", "Validation")) %>%
    ggplot(aes(x = train_size, y = r2, color = set)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Training" = "#4477AA", "Validation" = "#CC6677")) +
    labs(
      x = "Training Set Size",
      y = expression(R^2),
      title = expression("B. "*R^2*" vs Training Size"),
      color = "Dataset"
    )

  # Combine
  (p_rmse | p_r2) +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
}

# ============================================================================
# LOAD DATA AND MODELS
# ============================================================================

cat("Loading data and models...\n\n")

# Load coral master data
if (!file.exists(file.path(PATHS$objects, "coral_master.rds"))) {
  cat("Running 01_load_data.R to create coral_master...\n")
  source(here::here("scripts/01_load_data.R"))
}
coral_master <- readRDS(file.path(PATHS$objects, "coral_master.rds"))

# Check for ML results from 11_machine_learning.R
# Actual output files from 11_machine_learning.R:
ml_files <- list(
  machine_learning_models = file.path(PATHS$objects, "machine_learning_models.rds"),
  rf_model_abundance = file.path(PATHS$objects, "rf_model_abundance.rds"),
  rf_model_richness = file.path(PATHS$objects, "rf_model_richness.rds"),
  xgb_model_abundance = file.path(PATHS$objects, "xgb_model_abundance.rds")
)

# Check which ML models are available
ml_available <- file.exists(ml_files$machine_learning_models)
rf_available <- file.exists(ml_files$rf_model_abundance)
xgb_available <- file.exists(ml_files$xgb_model_abundance)

if (ml_available || rf_available) {
  cat("[OK] ML results found from 11_machine_learning.R\n")

  # Load main ML models object
  if (ml_available) {
    ml_models <- readRDS(ml_files$machine_learning_models)
    cat("  Loaded: machine_learning_models.rds\n")
  }

  # Load individual models
  if (rf_available) {
    rf_abundance <- readRDS(ml_files$rf_model_abundance)
    cat("  Loaded: rf_model_abundance.rds\n")
  }
  if (file.exists(ml_files$rf_model_richness)) {
    rf_richness <- readRDS(ml_files$rf_model_richness)
    cat("  Loaded: rf_model_richness.rds\n")
  }
  if (xgb_available) {
    xgb_abundance <- readRDS(ml_files$xgb_model_abundance)
    cat("  Loaded: xgb_model_abundance.rds\n")
  }

  # Load model comparison table if available
  model_comparison_file <- file.path(PATHS$tables, "model_comparison.csv")
  if (file.exists(model_comparison_file)) {
    model_comparison <- read_csv(model_comparison_file, show_col_types = FALSE)
    cat("  Loaded: model_comparison.csv\n")
  }

  # Load variable importance if available
  var_imp_file <- file.path(PATHS$tables, "variable_importance.csv")
  if (file.exists(var_imp_file)) {
    variable_importance <- read_csv(var_imp_file, show_col_types = FALSE)
    cat("  Loaded: variable_importance.csv\n")
  }

} else {
  cat("[INFO] ML results not found. Running evaluation on existing statistical models.\n")
  cat("       Run 11_machine_learning.R first for full ML evaluation.\n\n")

  # Load scaling analysis results if available
  if (file.exists(file.path(PATHS$objects, "scaling_analysis_results.rds"))) {
    cat("[OK] Loading scaling analysis results for evaluation\n")
    scaling_results <- readRDS(file.path(PATHS$objects, "scaling_analysis_results.rds"))
  } else {
    scaling_results <- NULL
  }
}

# ============================================================================
# PART A: CROSS-VALIDATION ANALYSIS
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PART A: CROSS-VALIDATION ANALYSIS\n")
cat("============================================================\n\n")

# Prepare analysis data
analysis_data <- coral_master %>%
  filter(!is.na(volume), volume > 0, !is.na(total_cafi)) %>%
  mutate(
    # Natural log, consistent with core pipeline (scripts 01-09)
    log_volume = log(volume),
    log_cafi = log10(total_cafi + 1)
  )

cat("Analysis dataset:", nrow(analysis_data), "corals\n\n")

# Define target variables and predictors
targets <- c("total_cafi", "otu_richness", "shannon")
predictors <- c("log_volume", "site")

# Additional predictors if neighborhood data available
if (sum(!is.na(analysis_data$n_neighbors)) > 30) {
  landscape_data <- analysis_data %>%
    filter(!is.na(n_neighbors), !is.na(total_neighbor_volume), !is.na(mean_neighbor_dist))

  extended_predictors <- c("log_volume", "n_neighbors",
                           "total_neighbor_volume", "mean_neighbor_dist", "site")
  cat("Neighborhood data available:", nrow(landscape_data), "corals\n")
} else {
  landscape_data <- NULL
  extended_predictors <- predictors
}

# Perform k-fold cross-validation
k <- 5
set.seed(42)

cv_results_list <- list()

for (target in targets) {
  cat("\nCross-validating:", target, "\n")

  # Remove NAs for this target
  cv_data <- analysis_data %>%
    filter(!is.na(.data[[target]]))

  if (nrow(cv_data) < 20) {
    cat("  Insufficient data, skipping\n")
    next
  }

  # Create folds
  folds <- createFolds(cv_data[[target]], k = k, list = TRUE)

  fold_results <- list()

  for (i in 1:k) {
    test_idx <- folds[[i]]
    train_data <- cv_data[-test_idx, ]
    test_data <- cv_data[test_idx, ]

    # Fit model based on target type
    tryCatch({
      if (target %in% c("total_cafi", "otu_richness")) {
        # Negative binomial for counts
        formula <- as.formula(paste(target, "~ log_volume + site"))
        model <- MASS::glm.nb(formula, data = train_data)
        predictions <- predict(model, newdata = test_data, type = "response")
      } else {
        # Linear model for continuous (Shannon)
        formula <- as.formula(paste(target, "~ log_volume + site"))
        model <- lm(formula, data = train_data)
        predictions <- predict(model, newdata = test_data)
      }

      # Calculate metrics
      metrics <- calc_metrics(test_data[[target]], predictions)

      fold_results[[i]] <- data.frame(
        target = target,
        fold = i,
        n_train = nrow(train_data),
        n_test = nrow(test_data),
        r_squared = metrics$r_squared,
        rmse = metrics$rmse,
        mae = metrics$mae,
        bias = metrics$bias
      )

    }, error = function(e) {
      cat("    Fold", i, "failed:", conditionMessage(e), "\n")
    })
  }

  if (length(fold_results) > 0) {
    cv_results_list[[target]] <- bind_rows(fold_results)
  }
}

# Combine CV results
cv_results_df <- bind_rows(cv_results_list)

if (nrow(cv_results_df) > 0) {
  # Summary statistics
  cv_summary <- cv_results_df %>%
    group_by(target) %>%
    summarise(
      n_folds = n(),
      mean_r2 = mean(r_squared, na.rm = TRUE),
      sd_r2 = sd(r_squared, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      sd_rmse = sd(rmse, na.rm = TRUE),
      mean_mae = mean(mae, na.rm = TRUE),
      sd_mae = sd(mae, na.rm = TRUE),
      cv_r2 = sd(r_squared, na.rm = TRUE) / mean(r_squared, na.rm = TRUE),
      .groups = "drop"
    )

  cat("\nCross-Validation Summary:\n")
  print(cv_summary)

  # Save CV results
  save_table(cv_results_df, "cv_results_detailed")
  cat("\n  Saved: cv_results_detailed.csv\n")

  # Plot CV fold performance
  p_cv_folds <- cv_results_df %>%
    pivot_longer(cols = c(r_squared, rmse, mae),
                 names_to = "metric", values_to = "value") %>%
    ggplot(aes(x = factor(fold), y = value, fill = target)) +
    geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
    facet_wrap(~metric, scales = "free_y") +
    scale_fill_viridis_d(option = "D", name = "Target") +
    labs(
      x = "CV Fold",
      y = "Value",
      title = "Cross-Validation Performance by Fold",
      subtitle = paste(k, "-fold CV")
    )

  ggsave(file.path(FIG_DIR, "cv_fold_performance.png"), p_cv_folds,
         width = 12, height = 5, dpi = 300, bg = "white")
  cat("  Saved: cv_fold_performance.png\n")
}

# ============================================================================
# PART A.2: LEAVE-ONE-SITE-OUT VALIDATION
# ============================================================================

cat("\n--- Leave-One-Site-Out Validation ---\n\n")

loso_results <- list()

for (target in targets) {
  cat("LOSO for:", target, "\n")

  cv_data <- analysis_data %>%
    filter(!is.na(.data[[target]]))

  sites <- unique(cv_data$site)

  for (held_out_site in sites) {
    train_data <- cv_data %>% filter(site != held_out_site)
    test_data <- cv_data %>% filter(site == held_out_site)

    if (nrow(test_data) < 3) next

    tryCatch({
      if (target %in% c("total_cafi", "otu_richness")) {
        formula <- as.formula(paste(target, "~ log_volume + site"))
        model <- MASS::glm.nb(formula, data = train_data)
        predictions <- predict(model, newdata = test_data, type = "response")
      } else {
        formula <- as.formula(paste(target, "~ log_volume + site"))
        model <- lm(formula, data = train_data)
        predictions <- predict(model, newdata = test_data)
      }

      metrics <- calc_metrics(test_data[[target]], predictions)

      loso_results[[length(loso_results) + 1]] <- data.frame(
        target = target,
        held_out_site = held_out_site,
        n_train = nrow(train_data),
        n_test = nrow(test_data),
        r_squared = metrics$r_squared,
        rmse = metrics$rmse,
        mae = metrics$mae,
        bias = metrics$bias
      )

      cat("  Held out", held_out_site, ": R² =", round(metrics$r_squared, 3),
          ", RMSE =", round(metrics$rmse, 2), "\n")

    }, error = function(e) {
      cat("  Error for", held_out_site, ":", conditionMessage(e), "\n")
    })
  }
}

loso_df <- bind_rows(loso_results)

if (nrow(loso_df) > 0) {
  cat("\nLOSO Summary:\n")
  loso_summary <- loso_df %>%
    group_by(target) %>%
    summarise(
      mean_r2 = mean(r_squared, na.rm = TRUE),
      range_r2 = paste(round(min(r_squared, na.rm = TRUE), 2), "-",
                       round(max(r_squared, na.rm = TRUE), 2)),
      .groups = "drop"
    )
  print(loso_summary)
}

# ============================================================================
# PART A.3: ML MODEL CV EVALUATION (if available)
# ============================================================================

if (exists("ml_models") || exists("rf_abundance")) {
  cat("\n--- Machine Learning Model Evaluation ---\n\n")

  # Summarize ML model comparison if available
  if (exists("model_comparison") && is.data.frame(model_comparison)) {
    cat("ML Model Comparison (from 11_machine_learning.R):\n")
    print(model_comparison)
    cat("\n")

    # Find best model by RMSE
    best_model_idx <- which.min(model_comparison$rmse_cv)
    if (length(best_model_idx) > 0) {
      cat("Best model by CV RMSE:", model_comparison$model[best_model_idx], "\n")
      cat("  RMSE:", round(model_comparison$rmse_cv[best_model_idx], 3), "\n")
      if ("r_squared_cv" %in% names(model_comparison)) {
        cat("  R²:", round(model_comparison$r_squared_cv[best_model_idx], 3), "\n")
      }
    }
  }

  # Evaluate Random Forest model if available
  if (exists("rf_abundance") && !is.null(rf_abundance)) {
    cat("\nRandom Forest (Abundance) Details:\n")

    # Check if it's a ranger object
    if ("ranger" %in% class(rf_abundance)) {
      cat("  Number of trees:", rf_abundance$num.trees, "\n")
      cat("  Min node size:", rf_abundance$min.node.size, "\n")
      cat("  OOB prediction error (MSE):", round(rf_abundance$prediction.error, 3), "\n")
      cat("  OOB R²:", round(rf_abundance$r.squared, 3), "\n")

      # Variable importance if available
      if (!is.null(rf_abundance$variable.importance)) {
        cat("\n  Top 5 Variable Importance:\n")
        imp_sorted <- sort(rf_abundance$variable.importance, decreasing = TRUE)
        for (i in 1:min(5, length(imp_sorted))) {
          cat("    ", names(imp_sorted)[i], ":", round(imp_sorted[i], 3), "\n")
        }
      }
    } else if ("caret" %in% class(rf_abundance) || "train" %in% class(rf_abundance)) {
      # caret train object
      cat("  Final model details from caret:\n")
      if (!is.null(rf_abundance$results)) {
        best_tune <- rf_abundance$bestTune
        cat("  Best hyperparameters:\n")
        print(best_tune)
      }
    }
  }

  # Evaluate XGBoost model if available
  if (exists("xgb_abundance") && !is.null(xgb_abundance)) {
    cat("\nXGBoost (Abundance) Details:\n")

    if ("xgb.Booster" %in% class(xgb_abundance)) {
      # Get parameters - use tryCatch as API changed in newer xgboost versions
      tryCatch({
        # Try newer xgboost API (xgb.attr for individual attributes)
        nrounds <- xgb.attr(xgb_abundance, "niter")
        if (!is.null(nrounds)) cat("  Number of rounds:", nrounds, "\n")
      }, error = function(e) {
        cat("  (Parameter extraction not available in this xgboost version)\n")
      })
    }
  }

  cat("\n")
}

# ============================================================================
# PART B: PREDICTION ACCURACY (OBSERVED VS PREDICTED)
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PART B: PREDICTION ACCURACY\n")
cat("============================================================\n\n")

# Fit models and generate predictions for each target
obs_pred_plots <- list()
all_predictions <- list()

for (target in targets) {
  cat("Evaluating:", target, "\n")

  eval_data <- analysis_data %>%
    filter(!is.na(.data[[target]]))

  if (nrow(eval_data) < 20) {
    cat("  Insufficient data\n")
    next
  }

  tryCatch({
    if (target %in% c("total_cafi", "otu_richness")) {
      formula <- as.formula(paste(target, "~ log_volume + site"))
      model <- MASS::glm.nb(formula, data = eval_data)
      predictions <- predict(model, type = "response")
    } else {
      formula <- as.formula(paste(target, "~ log_volume + site"))
      model <- lm(formula, data = eval_data)
      predictions <- predict(model)
    }

    observed <- eval_data[[target]]

    # Store predictions
    all_predictions[[target]] <- data.frame(
      coral_id = eval_data$coral_id,
      site = eval_data$site,
      observed = observed,
      predicted = predictions,
      residual = observed - predictions,
      target = target
    )

    # Calculate overall metrics
    metrics <- calc_metrics(observed, predictions)

    cat("  n =", metrics$n, "\n")
    cat("  R² =", round(metrics$r_squared, 3), "\n")
    cat("  RMSE =", round(metrics$rmse, 2), "\n")
    cat("  MAE =", round(metrics$mae, 2), "\n")
    cat("  Bias =", round(metrics$bias, 3), "\n\n")

    # Create obs vs pred plot
    title_map <- c(
      "total_cafi" = "Total CAFI Abundance",
      "otu_richness" = "Species Richness",
      "shannon" = "Shannon Diversity"
    )

    obs_pred_plots[[target]] <- plot_obs_vs_pred(
      observed = observed,
      predicted = predictions,
      title = title_map[target],
      site = eval_data$site,
      model_name = target
    )

  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
  })
}

# Combine and save observed vs predicted plots
if (length(obs_pred_plots) > 0) {
  n_plots <- length(obs_pred_plots)

  if (n_plots == 1) {
    combined_obs_pred <- obs_pred_plots[[1]]
  } else if (n_plots == 2) {
    combined_obs_pred <- obs_pred_plots[[1]] + obs_pred_plots[[2]]
  } else {
    combined_obs_pred <- wrap_plots(obs_pred_plots, ncol = 2)
  }

  combined_obs_pred <- combined_obs_pred +
    plot_annotation(
      title = "Observed vs Predicted: Model Performance",
      subtitle = "Dashed line = 1:1 | Red line = actual fit",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )

  ggsave(file.path(FIG_DIR, "observed_vs_predicted.png"), combined_obs_pred,
         width = 12, height = 8, dpi = 300, bg = "white")
  cat("  Saved: observed_vs_predicted.png\n")
}

# ============================================================================
# PART C: RESIDUAL DIAGNOSTICS
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PART C: RESIDUAL DIAGNOSTICS\n")
cat("============================================================\n\n")

residual_plots <- list()

for (target in names(all_predictions)) {
  pred_df <- all_predictions[[target]]

  cat("Residual diagnostics for:", target, "\n")

  # Test for normality
  shapiro_test <- if (nrow(pred_df) >= 5 && nrow(pred_df) <= 5000) {
    shapiro.test(pred_df$residual)
  } else {
    list(statistic = NA, p.value = NA)
  }

  cat("  Shapiro-Wilk W =", round(shapiro_test$statistic, 3),
      ", p =", format.pval(shapiro_test$p.value, 3), "\n")

  # Test for heteroscedasticity (Breusch-Pagan approximation)
  cor_resid_fitted <- cor(abs(pred_df$residual), pred_df$predicted)
  cat("  Cor(|residuals|, fitted) =", round(cor_resid_fitted, 3), "\n")

  if (abs(cor_resid_fitted) > 0.3) {
    cat("  Warning: Potential heteroscedasticity detected\n")
  }

  # Mean and variance of residuals
  cat("  Mean residual =", round(mean(pred_df$residual), 4), "\n")
  cat("  SD residuals =", round(sd(pred_df$residual), 3), "\n\n")

  # Create diagnostic plot
  title_map <- c(
    "total_cafi" = "Total CAFI",
    "otu_richness" = "Species Richness",
    "shannon" = "Shannon H'"
  )

  residual_plots[[target]] <- plot_residual_diagnostics(
    observed = pred_df$observed,
    predicted = pred_df$predicted,
    model_name = title_map[target]
  )
}

# Save residual diagnostics (for first target as main figure)
if (length(residual_plots) > 0) {
  main_resid_plot <- residual_plots[[1]]

  ggsave(file.path(FIG_DIR, "residual_diagnostics.png"), main_resid_plot,
         width = 12, height = 10, dpi = 300, bg = "white")
  cat("  Saved: residual_diagnostics.png\n")

  # Save all individual plots
  for (target in names(residual_plots)) {
    ggsave(file.path(FIG_DIR, paste0("residuals_", target, ".png")),
           residual_plots[[target]],
           width = 12, height = 10, dpi = 300, bg = "white")
  }
  cat("Saved: Individual residual plots to", FIG_DIR, "\n")
}

# ============================================================================
# PART D: SITE-SPECIFIC PERFORMANCE
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PART D: SITE-SPECIFIC PERFORMANCE\n")
cat("============================================================\n\n")

site_metrics_list <- list()

for (target in names(all_predictions)) {
  pred_df <- all_predictions[[target]]

  site_metrics <- calc_site_metrics(
    observed = pred_df$observed,
    predicted = pred_df$predicted,
    site = pred_df$site
  ) %>%
    mutate(target = target)

  site_metrics_list[[target]] <- site_metrics

  cat(target, "by site:\n")
  print(site_metrics %>% dplyr::select(site, n, r_squared, rmse, bias))
  cat("\n")
}

# Combine site metrics
all_site_metrics <- bind_rows(site_metrics_list)

# Identify potential biases
cat("--- Site-Specific Bias Analysis ---\n\n")

bias_flags <- all_site_metrics %>%
  group_by(target) %>%
  mutate(
    overall_rmse = mean(rmse, na.rm = TRUE),
    rel_rmse = rmse / overall_rmse,
    bias_flag = case_when(
      abs(bias) > overall_rmse * 0.5 ~ "High bias",
      rel_rmse > 1.5 ~ "High error",
      r_squared < 0 ~ "Poor fit",
      TRUE ~ "OK"
    )
  ) %>%
  ungroup()

flagged <- bias_flags %>% filter(bias_flag != "OK")

if (nrow(flagged) > 0) {
  cat("Potential issues detected:\n")
  print(flagged %>% dplyr::select(target, site, n, r_squared, bias, bias_flag))
} else {
  cat("No major site-specific issues detected.\n")
}

# Save site metrics
save_table(all_site_metrics, "site_metrics")
cat("\nSaved: site_metrics.csv\n")

# Plot site performance
p_site_r2 <- all_site_metrics %>%
  ggplot(aes(x = site, y = r_squared, fill = site)) +
  geom_col(alpha = 0.8) +
  facet_wrap(~target, scales = "free_y") +
  scale_fill_manual(values = SITE_COLORS) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    x = "Site",
    y = expression(R^2),
    title = "A. Model Performance by Site"
  ) +
  theme(legend.position = "none")

p_site_bias <- all_site_metrics %>%
  ggplot(aes(x = site, y = bias, fill = site)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~target, scales = "free_y") +
  scale_fill_manual(values = SITE_COLORS) +
  labs(
    x = "Site",
    y = "Bias (Predicted - Observed)",
    title = "B. Prediction Bias by Site"
  ) +
  theme(legend.position = "none")

p_site_combined <- (p_site_r2 / p_site_bias) +
  plot_annotation(
    title = "Site-Specific Model Performance",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(FIG_DIR, "site_performance.png"), p_site_combined,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("  Saved: site_performance.png\n")

# ============================================================================
# PART E: ECOLOGICAL VALIDITY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PART E: ECOLOGICAL VALIDITY\n")
cat("============================================================\n\n")

for (target in names(all_predictions)) {
  pred_df <- all_predictions[[target]]

  cat("Ecological validity check for:", target, "\n")

  # Observed data range
  obs_range <- range(pred_df$observed, na.rm = TRUE)
  pred_range <- range(pred_df$predicted, na.rm = TRUE)

  cat("  Observed range: [", round(obs_range[1], 2), ",", round(obs_range[2], 2), "]\n")
  cat("  Predicted range: [", round(pred_range[1], 2), ",", round(pred_range[2], 2), "]\n")

  # Check for unrealistic predictions
  if (target %in% c("total_cafi", "otu_richness")) {
    n_negative <- sum(pred_df$predicted < 0)
    n_extreme <- sum(pred_df$predicted > obs_range[2] * 2)

    cat("  Negative predictions:", n_negative,
        "(", round(n_negative/nrow(pred_df)*100, 1), "%)\n")
    cat("  Extreme predictions (>2x max observed):", n_extreme,
        "(", round(n_extreme/nrow(pred_df)*100, 1), "%)\n")

    if (n_negative > 0) {
      cat("  WARNING: Count predictions should not be negative!\n")
    }
  }

  if (target == "shannon") {
    n_negative <- sum(pred_df$predicted < 0)
    n_extreme <- sum(pred_df$predicted > 3)  # Shannon rarely exceeds 3-4

    cat("  Negative predictions:", n_negative, "\n")
    cat("  Predictions > 3:", n_extreme, "\n")
  }

  # Check prediction coverage
  within_range <- sum(pred_df$predicted >= obs_range[1] &
                       pred_df$predicted <= obs_range[2])
  cat("  Predictions within observed range:", within_range,
      "(", round(within_range/nrow(pred_df)*100, 1), "%)\n\n")
}

# ============================================================================
# PART F: UNCERTAINTY QUANTIFICATION
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PART F: UNCERTAINTY QUANTIFICATION\n")
cat("============================================================\n\n")

# Bootstrap confidence intervals for metrics
cat("--- Bootstrap Confidence Intervals ---\n\n")

n_boot <- 1000
set.seed(42)

boot_results <- list()

for (target in names(all_predictions)) {
  pred_df <- all_predictions[[target]]
  n <- nrow(pred_df)

  boot_metrics <- replicate(n_boot, {
    idx <- sample(1:n, n, replace = TRUE)
    m <- calc_metrics(pred_df$observed[idx], pred_df$predicted[idx])
    c(r2 = m$r_squared, rmse = m$rmse, mae = m$mae)
  })

  boot_ci <- data.frame(
    target = target,
    metric = c("R²", "RMSE", "MAE"),
    estimate = c(
      calc_metrics(pred_df$observed, pred_df$predicted)$r_squared,
      calc_metrics(pred_df$observed, pred_df$predicted)$rmse,
      calc_metrics(pred_df$observed, pred_df$predicted)$mae
    ),
    ci_lower = c(
      quantile(boot_metrics["r2", ], 0.025, na.rm = TRUE),
      quantile(boot_metrics["rmse", ], 0.025, na.rm = TRUE),
      quantile(boot_metrics["mae", ], 0.025, na.rm = TRUE)
    ),
    ci_upper = c(
      quantile(boot_metrics["r2", ], 0.975, na.rm = TRUE),
      quantile(boot_metrics["rmse", ], 0.975, na.rm = TRUE),
      quantile(boot_metrics["mae", ], 0.975, na.rm = TRUE)
    )
  )

  boot_results[[target]] <- boot_ci

  cat(target, ":\n")
  for (i in 1:nrow(boot_ci)) {
    cat("  ", boot_ci$metric[i], "=", round(boot_ci$estimate[i], 3),
        "[", round(boot_ci$ci_lower[i], 3), ",", round(boot_ci$ci_upper[i], 3), "]\n")
  }
  cat("\n")
}

boot_df <- bind_rows(boot_results)

# --- Prediction Intervals ---
cat("--- Prediction Intervals (Individual Predictions) ---\n\n")

prediction_intervals_list <- list()

for (target in names(all_predictions)) {
  pred_df <- all_predictions[[target]]

  cat(target, ":\n")

  # Calculate residual standard error
  residual_se <- sd(pred_df$residual, na.rm = TRUE)

  # Generate 95% prediction intervals
  # For GLM counts, use asymmetric intervals based on log-scale
  if (target %in% c("total_cafi", "otu_richness")) {
    # Use log-normal approximation for count data
    log_pred <- log(pred_df$predicted + 1)
    log_se <- sd(log(pred_df$observed + 1) - log(pred_df$predicted + 1), na.rm = TRUE)

    pred_df$pi_lower <- pmax(0, exp(log_pred - 1.96 * log_se) - 1)
    pred_df$pi_upper <- exp(log_pred + 1.96 * log_se) - 1
  } else {
    # Normal intervals for continuous variables
    pred_df$pi_lower <- pred_df$predicted - 1.96 * residual_se
    pred_df$pi_upper <- pred_df$predicted + 1.96 * residual_se
  }

  # Calculate coverage
  coverage <- mean(pred_df$observed >= pred_df$pi_lower &
                   pred_df$observed <= pred_df$pi_upper, na.rm = TRUE)

  # Calculate mean interval width
  mean_width <- mean(pred_df$pi_upper - pred_df$pi_lower, na.rm = TRUE)

  cat("  Residual SE:", round(residual_se, 3), "\n")
  cat("  95% PI coverage:", round(coverage * 100, 1), "%\n")
  cat("  Mean PI width:", round(mean_width, 2), "\n\n")

  # Store results
  pred_df$target <- target
  prediction_intervals_list[[target]] <- pred_df
}

# Combine all prediction intervals
all_prediction_intervals <- bind_rows(prediction_intervals_list)

# Save prediction intervals table
prediction_intervals_summary <- all_prediction_intervals %>%
  group_by(target) %>%
  summarise(
    n = n(),
    mean_predicted = mean(predicted, na.rm = TRUE),
    mean_pi_lower = mean(pi_lower, na.rm = TRUE),
    mean_pi_upper = mean(pi_upper, na.rm = TRUE),
    mean_pi_width = mean(pi_upper - pi_lower, na.rm = TRUE),
    coverage = mean(observed >= pi_lower & observed <= pi_upper, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# Save detailed prediction intervals
prediction_intervals_detail <- all_prediction_intervals %>%
  dplyr::select(target, coral_id, site, observed, predicted, pi_lower, pi_upper, residual)

write_csv(prediction_intervals_detail, file.path(PATHS$tables, "prediction_intervals.csv"))
cat("Saved: prediction_intervals.csv\n")

# Create prediction interval visualization
pi_plots <- list()

for (target in unique(all_prediction_intervals$target)) {
  target_data <- all_prediction_intervals %>%
    filter(target == !!target) %>%
    arrange(predicted)

  # Add row index for plotting
  target_data$index <- 1:nrow(target_data)

  title_map <- c(
    "total_cafi" = "Total CAFI Abundance",
    "otu_richness" = "Species Richness",
    "shannon" = "Shannon Diversity"
  )

  p <- ggplot(target_data, aes(x = index)) +
    geom_ribbon(aes(ymin = pi_lower, ymax = pi_upper),
                fill = "#377EB8", alpha = 0.3) +
    geom_line(aes(y = predicted), color = "#377EB8", linewidth = 1) +
    geom_point(aes(y = observed, color = site), size = 1.5, alpha = 0.7) +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    labs(
      x = "Coral (ordered by prediction)",
      y = title_map[target],
      title = paste(title_map[target], "- 95% Prediction Intervals"),
      subtitle = sprintf("Coverage: %.1f%% | Mean Width: %.2f",
                        prediction_intervals_summary$coverage[prediction_intervals_summary$target == target],
                        prediction_intervals_summary$mean_pi_width[prediction_intervals_summary$target == target])
    ) +
    theme_publication()

  pi_plots[[target]] <- p
}

# Combine and save
if (length(pi_plots) > 0) {
  combined_pi <- wrap_plots(pi_plots, ncol = 1)

  ggsave(file.path(FIG_DIR, "prediction_intervals.png"), combined_pi,
         width = 10, height = 4 * length(pi_plots), dpi = 300, bg = "white")
  cat("Saved: prediction_intervals.png to", FIG_DIR, "\n\n")
}

# ============================================================================
# PART G: LEARNING CURVES
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PART G: LEARNING CURVES\n")
cat("============================================================\n\n")

# Generate learning curves for primary target (total_cafi)
cat("Generating learning curves for total_cafi...\n")

lc_data <- analysis_data %>%
  filter(!is.na(total_cafi), !is.na(log_volume)) %>%
  dplyr::select(total_cafi, log_volume, site) %>%
  mutate(site = as.factor(site))

# Add dummy coding for site
lc_data_model <- lc_data %>%
  mutate(
    site_HAU = as.numeric(site == "HAU"),
    site_MAT = as.numeric(site == "MAT")
  )

learning_curve_data <- generate_learning_curve(
  data = lc_data_model,
  response_col = "total_cafi",
  predictor_cols = c("log_volume", "site_HAU", "site_MAT"),
  fractions = seq(0.2, 1, by = 0.1),
  n_iter = 10,
  model_type = "glm.nb"
)

if (nrow(learning_curve_data) > 0) {
  p_learning <- plot_learning_curves(
    learning_curve_data,
    title = "Learning Curves: Total CAFI Abundance Model"
  )

  ggsave(file.path(FIG_DIR, "learning_curves.png"), p_learning,
         width = 12, height = 5, dpi = 300, bg = "white")
  cat("  Saved: learning_curves.png\n")

  # Analyze learning curve
  lc_summary <- learning_curve_data %>%
    group_by(fraction) %>%
    summarise(
      train_rmse = mean(train_rmse, na.rm = TRUE),
      test_rmse = mean(test_rmse, na.rm = TRUE),
      gap = mean(train_rmse, na.rm = TRUE) - mean(test_rmse, na.rm = TRUE),
      .groups = "drop"
    )

  # Check for overfitting
  if (nrow(lc_summary) > 3) {
    final_gap <- tail(lc_summary$gap, 1)
    early_gap <- head(lc_summary$gap, 1)

    if (abs(final_gap) < 0.5 * abs(early_gap)) {
      cat("\nLearning curve indicates good generalization.\n")
    } else if (final_gap > 0) {
      cat("\nLearning curve suggests possible overfitting.\n")
    }

    # Check if performance is still improving
    last_improvement <- diff(tail(lc_summary$test_rmse, 2))
    if (last_improvement < -0.5) {
      cat("Model performance may benefit from more data.\n")
    } else {
      cat("Model appears to have saturated with current features.\n")
    }
  }
} else {
  cat("Could not generate learning curves.\n")
}

# ============================================================================
# PART H: MODEL SELECTION RECOMMENDATIONS
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PART H: MODEL SELECTION RECOMMENDATIONS\n")
cat("============================================================\n\n")

# Compile recommendations
recommendations <- list()

for (target in names(all_predictions)) {
  pred_df <- all_predictions[[target]]
  metrics <- calc_metrics(pred_df$observed, pred_df$predicted)

  # Get CV results if available
  cv_info <- if (exists("cv_summary") && target %in% cv_summary$target) {
    cv_summary %>% filter(target == !!target)
  } else {
    data.frame(mean_r2 = NA, cv_r2 = NA)
  }

  # Get site bias info
  site_bias_range <- all_site_metrics %>%
    filter(target == !!target) %>%
    summarise(max_bias = max(abs(bias), na.rm = TRUE)) %>%
    pull(max_bias)

  # Generate recommendation
  quality <- case_when(
    metrics$r_squared > 0.7 ~ "Excellent",
    metrics$r_squared > 0.5 ~ "Good",
    metrics$r_squared > 0.3 ~ "Moderate",
    metrics$r_squared > 0.1 ~ "Weak",
    TRUE ~ "Poor"
  )

  interpretability <- case_when(
    target %in% c("total_cafi", "otu_richness") ~ "High (GLM with clear coefficients)",
    TRUE ~ "High (Linear model)"
  )

  recommendation <- data.frame(
    target = target,
    model_type = ifelse(target %in% c("total_cafi", "otu_richness"),
                        "Negative Binomial GLM", "Linear Regression"),
    r_squared = round(metrics$r_squared, 3),
    rmse = round(metrics$rmse, 2),
    cv_r2_mean = round(cv_info$mean_r2[1], 3),
    cv_r2_cv = round(cv_info$cv_r2[1], 3),
    max_site_bias = round(site_bias_range, 2),
    quality = quality,
    interpretability = interpretability,
    recommendation = case_when(
      quality == "Excellent" ~ "Recommended for prediction and inference",
      quality == "Good" ~ "Suitable for most applications",
      quality == "Moderate" ~ "Use with caution, consider additional predictors",
      quality == "Weak" ~ "Not recommended for prediction",
      TRUE ~ "Model inadequate, needs revision"
    )
  )

  recommendations[[target]] <- recommendation
}

recommendations_df <- bind_rows(recommendations)

cat("MODEL RECOMMENDATIONS:\n")
cat("----------------------\n\n")

for (i in 1:nrow(recommendations_df)) {
  r <- recommendations_df[i, ]
  cat(r$target, ":\n")
  cat("  Model type:", r$model_type, "\n")
  cat("  R² =", r$r_squared, ", RMSE =", r$rmse, "\n")
  cat("  CV mean R² =", r$cv_r2_mean, "(CV =", r$cv_r2_cv, ")\n")
  cat("  Quality:", r$quality, "\n")
  cat("  Interpretability:", r$interpretability, "\n")
  cat("  Recommendation:", r$recommendation, "\n\n")
}

# Save recommendations
save_table(recommendations_df, "model_recommendations")
cat("Saved: model_recommendations.csv\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    MODEL EVALUATION COMPLETE\n")
cat("============================================================\n\n")

cat("EVALUATION SUMMARY:\n")
cat("-------------------\n\n")

cat("Sample size:", nrow(analysis_data), "corals\n")
cat("Targets evaluated:", paste(targets, collapse = ", "), "\n\n")

cat("Cross-Validation:\n")
if (exists("cv_summary") && nrow(cv_summary) > 0) {
  for (i in 1:nrow(cv_summary)) {
    cat("  ", cv_summary$target[i], ": R² =", round(cv_summary$mean_r2[i], 3),
        "+/-", round(cv_summary$sd_r2[i], 3), "\n")
  }
}

cat("\nBest performing model:",
    recommendations_df$target[which.max(recommendations_df$r_squared)],
    "(R² =", max(recommendations_df$r_squared), ")\n")

cat("\nOUTPUT FILES:\n")
cat("  Figures:\n")
cat("    - output/figures/observed_vs_predicted.png\n")
cat("    - output/figures/residual_diagnostics.png\n")
cat("    - output/figures/learning_curves.png\n")
cat("    - output/figures/site_performance.png\n")
cat("    - output/figures/12_evaluation/ (additional diagnostics)\n")
cat("  Tables:\n")
cat("    - output/tables/cv_results_detailed.csv\n")
cat("    - output/tables/site_metrics.csv\n")
cat("    - output/tables/model_recommendations.csv\n")
cat("    - output/tables/prediction_intervals.csv\n\n")

cat("============================================================\n")
cat("    Evaluation complete. See outputs for detailed results.\n")
cat("============================================================\n\n")
