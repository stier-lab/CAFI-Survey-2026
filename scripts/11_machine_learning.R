#!/usr/bin/env Rscript
# ============================================================================
# 11_machine_learning.R - Machine Learning Models for CAFI Prediction
# ============================================================================
#
# PURPOSE: Build and optimize predictive models for CAFI abundance and diversity
#
# TARGET VARIABLES:
#   - cafi_abundance: Total CAFI count per coral
#   - cafi_richness: Species count per coral
#   - shannon: Shannon diversity index
#
# MODELS IMPLEMENTED:
#   A. Baseline: GLM (negative binomial for counts, Gaussian for Shannon)
#   B. Random Forest: ranger package with hyperparameter tuning
#   C. Gradient Boosting: XGBoost with cross-validation
#
# CROSS-VALIDATION:
#   - 5-fold CV stratified by site
#   - Leave-one-site-out validation
#
# OUTPUTS:
#   - output/objects/rf_model_abundance.rds
#   - output/objects/rf_model_richness.rds
#   - output/objects/xgb_model_abundance.rds
#   - output/tables/model_comparison.csv
#   - output/tables/variable_importance.csv
#   - output/figures/variable_importance.png
#   - output/figures/partial_dependence.png
#   - output/figures/model_comparison.png
#
# DEPENDENCIES: 00_setup.R, 01_load_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-21
# ============================================================================

cat("\n========================================\n")
cat("11: Machine Learning Models\n")
cat("========================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Additional ML packages
suppressPackageStartupMessages({
  if (!requireNamespace("ranger", quietly = TRUE)) {
    cat("Installing ranger...\n")
    install.packages("ranger", repos = "https://cran.rstudio.com/")
  }
  library(ranger)

  if (!requireNamespace("xgboost", quietly = TRUE)) {
    cat("Installing xgboost...\n")
    install.packages("xgboost", repos = "https://cran.rstudio.com/")
  }
  library(xgboost)

  if (!requireNamespace("caret", quietly = TRUE)) {
    cat("Installing caret...\n")
    install.packages("caret", repos = "https://cran.rstudio.com/")
  }
  library(caret)
})

# Install and load optional packages for advanced visualizations
has_pdp <- FALSE
has_shapviz <- FALSE

tryCatch({
  if (!requireNamespace("pdp", quietly = TRUE)) {
    cat("Installing pdp package for partial dependence plots...\n")
    install.packages("pdp", repos = "https://cran.rstudio.com/", quiet = TRUE)
  }
  library(pdp)
  has_pdp <- TRUE
}, error = function(e) {
  cat("Note: pdp package not available\n")
})

tryCatch({
  if (!requireNamespace("shapviz", quietly = TRUE)) {
    cat("Installing shapviz package for SHAP values...\n")
    install.packages("shapviz", repos = "https://cran.rstudio.com/", quiet = TRUE)
  }
  library(shapviz)
  has_shapviz <- TRUE
}, error = function(e) {
  cat("Note: shapviz package not available\n")
})

cat("[OK] ML packages loaded\n")
cat("   - ranger: Random Forest\n")
cat("   - xgboost: Gradient Boosting\n")
cat("   - caret: Model training/CV\n")
if (has_pdp) cat("   - pdp: Partial dependence plots\n")
if (has_shapviz) cat("   - shapviz: SHAP values\n")
cat("\n")

# Create figure subdirectory
fig_dir <- file.path(PATHS$figures, "machine_learning")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")

# Check if objects exist, otherwise run data loading script
if (!file.exists(file.path(PATHS$objects, "coral_master.rds"))) {
  cat("  Running 01_load_data.R to generate data objects...\n")
  source(here::here("scripts/01_load_data.R"))
}

coral_master <- readRDS(file.path(PATHS$objects, "coral_master.rds"))
cafi_clean <- readRDS(file.path(PATHS$objects, "cafi_clean.rds"))

cat("  Coral master:", nrow(coral_master), "corals,", ncol(coral_master), "variables\n")
cat("  CAFI clean:", nrow(cafi_clean), "records\n\n")

# ============================================================================
# FEATURE ENGINEERING
# ============================================================================

cat("Engineering features for ML models...\n")

# Create additional CAFI features per coral
cafi_features <- cafi_clean %>%
  group_by(coral_id) %>%
  summarise(
    # Size metrics
    mean_cafi_size = mean(size_mm, na.rm = TRUE),
    sd_cafi_size = sd(size_mm, na.rm = TRUE),
    max_cafi_size = if(any(!is.na(size_mm))) max(size_mm, na.rm = TRUE) else NA_real_,
    min_cafi_size = if(any(!is.na(size_mm))) min(size_mm, na.rm = TRUE) else NA_real_,

    # Type proportions
    prop_crabs = sum(type == "crab") / n(),
    prop_shrimps = sum(type == "shrimp") / n(),
    prop_fish = sum(type == "fish") / n(),
    prop_snails = sum(type == "snail") / n(),

    # Dominant species proportion
    dominant_prop = if(n() > 0) max(table(species)) / n() else NA_real_,

    .groups = "drop"
  )

# Create ML dataset
ml_data <- coral_master %>%
  left_join(cafi_features, by = "coral_id") %>%
  mutate(
    # Ensure target variables exist
    cafi_abundance = total_cafi,
    cafi_richness = otu_richness,

    # Log transformations for predictors
    log_volume = log10(volume + 1),
    log_total_neighbor_volume = log10(coalesce(total_neighbor_volume, 0) + 1),

    # Derived features
    cafi_density = total_cafi / (volume + 1) * 1000,

    # Encode factors for ML models
    morphotype_factor = as.factor(morphotype),
    site_factor = as.factor(site),

    # Numeric encoding for tree models
    morphotype_num = as.numeric(morphotype_factor),
    site_num = as.numeric(site_factor)
  ) %>%
  # Handle missing values - use any_of() to handle columns that may not exist
  mutate(
    across(any_of(c("mean_cafi_size", "sd_cafi_size", "max_cafi_size", "min_cafi_size",
             "prop_crabs", "prop_shrimps", "prop_fish", "prop_snails", "dominant_prop")),
           ~replace_na(., 0)),
    across(any_of(c("n_neighbors", "total_neighbor_volume", "mean_neighbor_dist")),
           ~replace_na(., 0))
  ) %>%
  filter(!is.na(volume), volume > 0, !is.na(depth_m))

cat("  ML dataset:", nrow(ml_data), "samples,", ncol(ml_data), "features\n")

# Define feature sets
feature_cols_base <- c("log_volume", "depth_m", "site_num", "morphotype_num")

# Extended features (includes neighborhood if available)
feature_cols_extended <- c(
  feature_cols_base,
  "n_neighbors",
  "log_total_neighbor_volume",
  "mean_neighbor_dist"
)

# Filter to valid features
feature_cols_base <- feature_cols_base[feature_cols_base %in% names(ml_data)]
feature_cols_extended <- feature_cols_extended[feature_cols_extended %in% names(ml_data)]

cat("  Base features:", length(feature_cols_base), "\n")
cat("  Extended features:", length(feature_cols_extended), "\n\n")

# ============================================================================
# TRAIN/TEST SPLIT
# ============================================================================

cat("Creating train/test split...\n")

set.seed(123)

# Stratified split by site
train_indices <- createDataPartition(ml_data$site_factor, p = 0.8, list = FALSE)
train_data <- ml_data[train_indices, ]
test_data <- ml_data[-train_indices, ]

cat("  Training set:", nrow(train_data), "samples\n")
cat("  Test set:", nrow(test_data), "samples\n")
cat("  Sites in train:", paste(unique(train_data$site), collapse = ", "), "\n")
cat("  Sites in test:", paste(unique(test_data$site), collapse = ", "), "\n\n")

# ============================================================================
# CROSS-VALIDATION SETUP
# ============================================================================

cat("Setting up cross-validation...\n")

# 5-fold CV stratified by site
folds_5 <- createFolds(train_data$site_factor, k = 5, returnTrain = TRUE)

# Leave-one-site-out CV
sites <- unique(train_data$site)
folds_loso <- lapply(sites, function(s) {
  which(train_data$site != s)
})
names(folds_loso) <- paste0("Leave_", sites)

cat("  5-fold CV created\n")
cat("  Leave-one-site-out CV:", length(sites), "folds\n\n")

# ============================================================================
# PART 1: BASELINE GLM MODELS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 1: BASELINE GLM MODELS\n")
cat("------------------------------------------------------------\n\n")

results_list <- list()

# --- 1A. GLM for Abundance (Negative Binomial) ---
cat("1A. GLM: CAFI Abundance (Negative Binomial)...\n")

glm_abundance_formula <- as.formula(paste(
  "cafi_abundance ~", paste(feature_cols_base, collapse = " + ")
))

glm_abundance <- tryCatch({
  glm.nb(glm_abundance_formula, data = train_data)
}, error = function(e) {
  cat("  Warning: glm.nb failed, using Poisson\n")
  glm(glm_abundance_formula, data = train_data, family = poisson)
})

# Predictions
train_data$glm_pred_abundance <- predict(glm_abundance, train_data, type = "response")
test_data$glm_pred_abundance <- predict(glm_abundance, test_data, type = "response")

# Metrics
glm_abund_rmse_train <- sqrt(mean((train_data$cafi_abundance - train_data$glm_pred_abundance)^2))
glm_abund_rmse_test <- sqrt(mean((test_data$cafi_abundance - test_data$glm_pred_abundance)^2))
glm_abund_r2_test <- cor(test_data$cafi_abundance, test_data$glm_pred_abundance, use = "complete")^2
glm_abund_mae_test <- mean(abs(test_data$cafi_abundance - test_data$glm_pred_abundance))

results_list$glm_abundance <- data.frame(
  model = "GLM (NB)", target = "abundance",
  rmse_train = glm_abund_rmse_train, rmse_test = glm_abund_rmse_test,
  r2_test = glm_abund_r2_test, mae_test = glm_abund_mae_test
)

cat("  RMSE (test):", round(glm_abund_rmse_test, 2), "\n")
cat("  R2 (test):", round(glm_abund_r2_test, 3), "\n")
cat("  MAE (test):", round(glm_abund_mae_test, 2), "\n\n")

# --- 1B. GLM for Richness (Poisson) ---
cat("1B. GLM: CAFI Richness (Poisson)...\n")

glm_richness_formula <- as.formula(paste(
  "cafi_richness ~", paste(feature_cols_base, collapse = " + ")
))

glm_richness <- glm(glm_richness_formula, data = train_data, family = poisson)

train_data$glm_pred_richness <- predict(glm_richness, train_data, type = "response")
test_data$glm_pred_richness <- predict(glm_richness, test_data, type = "response")

glm_rich_rmse_test <- sqrt(mean((test_data$cafi_richness - test_data$glm_pred_richness)^2))
glm_rich_r2_test <- cor(test_data$cafi_richness, test_data$glm_pred_richness, use = "complete")^2
glm_rich_mae_test <- mean(abs(test_data$cafi_richness - test_data$glm_pred_richness))

results_list$glm_richness <- data.frame(
  model = "GLM (Poisson)", target = "richness",
  rmse_train = sqrt(mean((train_data$cafi_richness - train_data$glm_pred_richness)^2)),
  rmse_test = glm_rich_rmse_test,
  r2_test = glm_rich_r2_test, mae_test = glm_rich_mae_test
)

cat("  RMSE (test):", round(glm_rich_rmse_test, 2), "\n")
cat("  R2 (test):", round(glm_rich_r2_test, 3), "\n\n")

# --- 1C. GLM for Shannon (Gaussian) ---
cat("1C. GLM: Shannon Diversity (Gaussian)...\n")

glm_shannon_formula <- as.formula(paste(
  "shannon ~", paste(feature_cols_base, collapse = " + ")
))

glm_shannon <- lm(glm_shannon_formula, data = train_data)

train_data$glm_pred_shannon <- predict(glm_shannon, train_data)
test_data$glm_pred_shannon <- predict(glm_shannon, test_data)

glm_shan_rmse_test <- sqrt(mean((test_data$shannon - test_data$glm_pred_shannon)^2))
glm_shan_r2_test <- cor(test_data$shannon, test_data$glm_pred_shannon, use = "complete")^2
glm_shan_mae_test <- mean(abs(test_data$shannon - test_data$glm_pred_shannon))

results_list$glm_shannon <- data.frame(
  model = "GLM (Gaussian)", target = "shannon",
  rmse_train = sqrt(mean((train_data$shannon - train_data$glm_pred_shannon)^2)),
  rmse_test = glm_shan_rmse_test,
  r2_test = glm_shan_r2_test, mae_test = glm_shan_mae_test
)

cat("  RMSE (test):", round(glm_shan_rmse_test, 3), "\n")
cat("  R2 (test):", round(glm_shan_r2_test, 3), "\n\n")

# ============================================================================
# PART 2: RANDOM FOREST MODELS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 2: RANDOM FOREST MODELS\n")
cat("------------------------------------------------------------\n\n")

# --- 2A. Random Forest for Abundance ---
cat("2A. Random Forest: CAFI Abundance...\n")

# Hyperparameter grid
rf_grid <- expand.grid(
  mtry = c(2, 3, floor(sqrt(length(feature_cols_extended)))),
  min.node.size = c(3, 5, 10),
  num.trees = c(500, 1000)
)

# Use OOB error to select best hyperparameters
rf_results <- list()
for (i in 1:nrow(rf_grid)) {
  rf_temp <- ranger(
    cafi_abundance ~ .,
    data = train_data[, c("cafi_abundance", feature_cols_extended)],
    num.trees = rf_grid$num.trees[i],
    mtry = rf_grid$mtry[i],
    min.node.size = rf_grid$min.node.size[i],
    importance = "impurity"
  )
  rf_results[[i]] <- data.frame(
    rf_grid[i, ],
    oob_rmse = sqrt(rf_temp$prediction.error)
  )
}

rf_tuning <- bind_rows(rf_results)
best_rf_params <- rf_tuning[which.min(rf_tuning$oob_rmse), ]

cat("  Best params: mtry=", best_rf_params$mtry,
    ", min.node.size=", best_rf_params$min.node.size,
    ", num.trees=", best_rf_params$num.trees, "\n")
cat("  OOB RMSE:", round(best_rf_params$oob_rmse, 2), "\n")

# Train final RF abundance model
rf_abundance <- ranger(
  cafi_abundance ~ .,
  data = train_data[, c("cafi_abundance", feature_cols_extended)],
  num.trees = best_rf_params$num.trees,
  mtry = best_rf_params$mtry,
  min.node.size = best_rf_params$min.node.size,
  importance = "impurity"
)

train_data$rf_pred_abundance <- predict(rf_abundance, train_data)$predictions
test_data$rf_pred_abundance <- predict(rf_abundance, test_data)$predictions

rf_abund_rmse_test <- sqrt(mean((test_data$cafi_abundance - test_data$rf_pred_abundance)^2))
rf_abund_r2_test <- cor(test_data$cafi_abundance, test_data$rf_pred_abundance, use = "complete")^2
rf_abund_mae_test <- mean(abs(test_data$cafi_abundance - test_data$rf_pred_abundance))

results_list$rf_abundance <- data.frame(
  model = "Random Forest", target = "abundance",
  rmse_train = sqrt(rf_abundance$prediction.error),
  rmse_test = rf_abund_rmse_test,
  r2_test = rf_abund_r2_test, mae_test = rf_abund_mae_test
)

cat("  RMSE (test):", round(rf_abund_rmse_test, 2), "\n")
cat("  R2 (test):", round(rf_abund_r2_test, 3), "\n\n")

# --- 2B. Random Forest for Richness ---
cat("2B. Random Forest: CAFI Richness...\n")

rf_richness <- ranger(
  cafi_richness ~ .,
  data = train_data[, c("cafi_richness", feature_cols_extended)],
  num.trees = best_rf_params$num.trees,
  mtry = best_rf_params$mtry,
  min.node.size = best_rf_params$min.node.size,
  importance = "impurity"
)

train_data$rf_pred_richness <- predict(rf_richness, train_data)$predictions
test_data$rf_pred_richness <- predict(rf_richness, test_data)$predictions

rf_rich_rmse_test <- sqrt(mean((test_data$cafi_richness - test_data$rf_pred_richness)^2))
rf_rich_r2_test <- cor(test_data$cafi_richness, test_data$rf_pred_richness, use = "complete")^2

results_list$rf_richness <- data.frame(
  model = "Random Forest", target = "richness",
  rmse_train = sqrt(rf_richness$prediction.error),
  rmse_test = rf_rich_rmse_test,
  r2_test = rf_rich_r2_test,
  mae_test = mean(abs(test_data$cafi_richness - test_data$rf_pred_richness))
)

cat("  RMSE (test):", round(rf_rich_rmse_test, 2), "\n")
cat("  R2 (test):", round(rf_rich_r2_test, 3), "\n\n")

# --- 2C. Random Forest for Shannon ---
cat("2C. Random Forest: Shannon Diversity...\n")

rf_shannon <- ranger(
  shannon ~ .,
  data = train_data[, c("shannon", feature_cols_extended)],
  num.trees = best_rf_params$num.trees,
  mtry = best_rf_params$mtry,
  min.node.size = best_rf_params$min.node.size,
  importance = "impurity"
)

train_data$rf_pred_shannon <- predict(rf_shannon, train_data)$predictions
test_data$rf_pred_shannon <- predict(rf_shannon, test_data)$predictions

rf_shan_rmse_test <- sqrt(mean((test_data$shannon - test_data$rf_pred_shannon)^2))
rf_shan_r2_test <- cor(test_data$shannon, test_data$rf_pred_shannon, use = "complete")^2

results_list$rf_shannon <- data.frame(
  model = "Random Forest", target = "shannon",
  rmse_train = sqrt(rf_shannon$prediction.error),
  rmse_test = rf_shan_rmse_test,
  r2_test = rf_shan_r2_test,
  mae_test = mean(abs(test_data$shannon - test_data$rf_pred_shannon))
)

cat("  RMSE (test):", round(rf_shan_rmse_test, 3), "\n")
cat("  R2 (test):", round(rf_shan_r2_test, 3), "\n\n")

# ============================================================================
# PART 3: XGBOOST MODELS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 3: XGBOOST MODELS\n")
cat("------------------------------------------------------------\n\n")

# Prepare XGBoost matrices
xgb_train_matrix <- as.matrix(train_data[, feature_cols_extended])
xgb_test_matrix <- as.matrix(test_data[, feature_cols_extended])

# --- 3A. XGBoost for Abundance ---
cat("3A. XGBoost: CAFI Abundance...\n")

dtrain_abund <- xgb.DMatrix(data = xgb_train_matrix, label = train_data$cafi_abundance)
dtest_abund <- xgb.DMatrix(data = xgb_test_matrix, label = test_data$cafi_abundance)

# Cross-validation for hyperparameter tuning
xgb_params_abund <- list(
  objective = "count:poisson",
  eval_metric = "rmse",
  max_depth = 6,
  eta = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8
)

set.seed(123)
xgb_cv_abund <- xgb.cv(
  params = xgb_params_abund,
  data = dtrain_abund,
  nrounds = 500,
  nfold = 5,
  early_stopping_rounds = 20,
  verbose = 0
)

# Use best_iteration if available, otherwise default to 100
best_nrounds_abund <- if (!is.null(xgb_cv_abund$best_iteration) &&
                          length(xgb_cv_abund$best_iteration) > 0 &&
                          xgb_cv_abund$best_iteration > 0) {
  xgb_cv_abund$best_iteration
} else {
  100
}
cat("  Best nrounds:", best_nrounds_abund, "\n")

# Train final XGBoost abundance model
xgb_abundance <- xgb.train(
  params = xgb_params_abund,
  data = dtrain_abund,
  nrounds = best_nrounds_abund,
  verbose = 0
)

train_data$xgb_pred_abundance <- predict(xgb_abundance, dtrain_abund)
test_data$xgb_pred_abundance <- predict(xgb_abundance, dtest_abund)

xgb_abund_rmse_test <- sqrt(mean((test_data$cafi_abundance - test_data$xgb_pred_abundance)^2))
xgb_abund_r2_test <- cor(test_data$cafi_abundance, test_data$xgb_pred_abundance, use = "complete")^2
xgb_abund_mae_test <- mean(abs(test_data$cafi_abundance - test_data$xgb_pred_abundance))

results_list$xgb_abundance <- data.frame(
  model = "XGBoost", target = "abundance",
  rmse_train = sqrt(mean((train_data$cafi_abundance - train_data$xgb_pred_abundance)^2)),
  rmse_test = xgb_abund_rmse_test,
  r2_test = xgb_abund_r2_test, mae_test = xgb_abund_mae_test
)

cat("  RMSE (test):", round(xgb_abund_rmse_test, 2), "\n")
cat("  R2 (test):", round(xgb_abund_r2_test, 3), "\n\n")

# --- 3B. XGBoost for Richness ---
cat("3B. XGBoost: CAFI Richness...\n")

dtrain_rich <- xgb.DMatrix(data = xgb_train_matrix, label = train_data$cafi_richness)
dtest_rich <- xgb.DMatrix(data = xgb_test_matrix, label = test_data$cafi_richness)

xgb_params_rich <- list(
  objective = "count:poisson",
  eval_metric = "rmse",
  max_depth = 5,
  eta = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8
)

set.seed(123)
xgb_cv_rich <- xgb.cv(
  params = xgb_params_rich,
  data = dtrain_rich,
  nrounds = 500,
  nfold = 5,
  early_stopping_rounds = 20,
  verbose = 0
)

best_nrounds_rich <- if (!is.null(xgb_cv_rich$best_iteration) &&
                         length(xgb_cv_rich$best_iteration) > 0 &&
                         xgb_cv_rich$best_iteration > 0) {
  xgb_cv_rich$best_iteration
} else {
  100
}

xgb_richness <- xgb.train(
  params = xgb_params_rich,
  data = dtrain_rich,
  nrounds = best_nrounds_rich,
  verbose = 0
)

train_data$xgb_pred_richness <- predict(xgb_richness, dtrain_rich)
test_data$xgb_pred_richness <- predict(xgb_richness, dtest_rich)

xgb_rich_rmse_test <- sqrt(mean((test_data$cafi_richness - test_data$xgb_pred_richness)^2))
xgb_rich_r2_test <- cor(test_data$cafi_richness, test_data$xgb_pred_richness, use = "complete")^2

results_list$xgb_richness <- data.frame(
  model = "XGBoost", target = "richness",
  rmse_train = sqrt(mean((train_data$cafi_richness - train_data$xgb_pred_richness)^2)),
  rmse_test = xgb_rich_rmse_test,
  r2_test = xgb_rich_r2_test,
  mae_test = mean(abs(test_data$cafi_richness - test_data$xgb_pred_richness))
)

cat("  RMSE (test):", round(xgb_rich_rmse_test, 2), "\n")
cat("  R2 (test):", round(xgb_rich_r2_test, 3), "\n\n")

# --- 3C. XGBoost for Shannon ---
cat("3C. XGBoost: Shannon Diversity...\n")

dtrain_shan <- xgb.DMatrix(data = xgb_train_matrix, label = train_data$shannon)
dtest_shan <- xgb.DMatrix(data = xgb_test_matrix, label = test_data$shannon)

xgb_params_shan <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  max_depth = 5,
  eta = 0.05,
  subsample = 0.8,
  colsample_bytree = 0.8
)

set.seed(123)
xgb_cv_shan <- xgb.cv(
  params = xgb_params_shan,
  data = dtrain_shan,
  nrounds = 500,
  nfold = 5,
  early_stopping_rounds = 20,
  verbose = 0
)

best_nrounds_shan <- if (!is.null(xgb_cv_shan$best_iteration) &&
                         length(xgb_cv_shan$best_iteration) > 0 &&
                         xgb_cv_shan$best_iteration > 0) {
  xgb_cv_shan$best_iteration
} else {
  100
}

xgb_shannon <- xgb.train(
  params = xgb_params_shan,
  data = dtrain_shan,
  nrounds = best_nrounds_shan,
  verbose = 0
)

train_data$xgb_pred_shannon <- predict(xgb_shannon, dtrain_shan)
test_data$xgb_pred_shannon <- predict(xgb_shannon, dtest_shan)

xgb_shan_rmse_test <- sqrt(mean((test_data$shannon - test_data$xgb_pred_shannon)^2))
xgb_shan_r2_test <- cor(test_data$shannon, test_data$xgb_pred_shannon, use = "complete")^2

results_list$xgb_shannon <- data.frame(
  model = "XGBoost", target = "shannon",
  rmse_train = sqrt(mean((train_data$shannon - train_data$xgb_pred_shannon)^2)),
  rmse_test = xgb_shan_rmse_test,
  r2_test = xgb_shan_r2_test,
  mae_test = mean(abs(test_data$shannon - test_data$xgb_pred_shannon))
)

cat("  RMSE (test):", round(xgb_shan_rmse_test, 3), "\n")
cat("  R2 (test):", round(xgb_shan_r2_test, 3), "\n\n")

# ============================================================================
# PART 4: LEAVE-ONE-SITE-OUT CROSS-VALIDATION
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 4: LEAVE-ONE-SITE-OUT CROSS-VALIDATION\n")
cat("------------------------------------------------------------\n\n")

loso_results <- list()

for (test_site in sites) {
  cat("  Testing on site:", test_site, "\n")

  train_loso <- ml_data %>% filter(site != test_site)
  test_loso <- ml_data %>% filter(site == test_site)

  # RF model
  rf_loso <- ranger(
    cafi_abundance ~ .,
    data = train_loso[, c("cafi_abundance", feature_cols_extended)],
    num.trees = 500,
    mtry = best_rf_params$mtry,
    min.node.size = best_rf_params$min.node.size
  )

  pred_loso <- predict(rf_loso, test_loso)$predictions

  loso_results[[test_site]] <- data.frame(
    site = test_site,
    n_test = nrow(test_loso),
    rmse = sqrt(mean((test_loso$cafi_abundance - pred_loso)^2)),
    r2 = cor(test_loso$cafi_abundance, pred_loso, use = "complete")^2,
    mae = mean(abs(test_loso$cafi_abundance - pred_loso))
  )
}

loso_summary <- bind_rows(loso_results)
cat("\nLeave-One-Site-Out Results:\n")
print(loso_summary)
cat("\nMean RMSE across sites:", round(mean(loso_summary$rmse), 2), "\n")
cat("Mean R2 across sites:", round(mean(loso_summary$r2), 3), "\n\n")

# ============================================================================
# PART 5: MODEL COMPARISON
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 5: MODEL COMPARISON\n")
cat("------------------------------------------------------------\n\n")

# Combine all results
model_comparison <- bind_rows(results_list) %>%
  arrange(target, rmse_test)

cat("Model Comparison Table:\n")
print(model_comparison)
cat("\n")

# Save comparison table
write_csv(model_comparison, file.path(PATHS$tables, "model_comparison.csv"))
cat("Saved: model_comparison.csv\n\n")

# Model comparison figure - show R² directly (easier to interpret)
p_comparison <- model_comparison %>%
  mutate(model_target = paste(model, "-", target)) %>%
  ggplot(aes(x = reorder(model_target, r2_test), y = r2_test, fill = target)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = round(r2_test, 2)), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2", name = "Target") +
  scale_y_continuous(limits = c(0, max(model_comparison$r2_test, na.rm = TRUE) * 1.15)) +
  labs(
    title = "Model Comparison: Test Set Performance",
    subtitle = "Higher R² indicates better predictive performance",
    x = "Model - Target",
    y = expression("R"^2~"(Test Set)")
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "model_comparison.png"),
       p_comparison, width = 10, height = 8, dpi = 300)
cat("Saved: model_comparison.png\n\n")

# ============================================================================
# PART 6: VARIABLE IMPORTANCE
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 6: VARIABLE IMPORTANCE\n")
cat("------------------------------------------------------------\n\n")

# RF Variable Importance (impurity-based)
rf_importance <- data.frame(
  variable = names(rf_abundance$variable.importance),
  importance_rf = rf_abundance$variable.importance
) %>%
  mutate(importance_rf_scaled = importance_rf / max(importance_rf) * 100)

# XGBoost Variable Importance
xgb_importance_raw <- xgb.importance(
  feature_names = feature_cols_extended,
  model = xgb_abundance
)

xgb_importance <- data.frame(
  variable = xgb_importance_raw$Feature,
  importance_xgb = xgb_importance_raw$Gain * 100
)

# Combine
variable_importance <- rf_importance %>%
  left_join(xgb_importance, by = "variable") %>%
  mutate(
    importance_xgb = replace_na(importance_xgb, 0),
    importance_mean = (importance_rf_scaled + importance_xgb) / 2
  ) %>%
  arrange(desc(importance_mean))

cat("Variable Importance (Top 10):\n")
print(head(variable_importance, 10))
cat("\n")

# Save importance table
write_csv(variable_importance, file.path(PATHS$tables, "variable_importance.csv"))
cat("Saved: variable_importance.csv\n\n")

# Variable importance figure
importance_long <- variable_importance %>%
  pivot_longer(
    cols = c(importance_rf_scaled, importance_xgb),
    names_to = "model",
    values_to = "importance"
  ) %>%
  mutate(
    model = case_when(
      model == "importance_rf_scaled" ~ "Random Forest",
      model == "importance_xgb" ~ "XGBoost"
    )
  )

p_importance <- importance_long %>%
  ggplot(aes(x = reorder(variable, importance_mean), y = importance, fill = model)) +
  geom_col(position = "dodge", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Random Forest" = "#0072B2", "XGBoost" = "#009E73")) +
  labs(
    title = "Variable Importance for CAFI Abundance Prediction",
    x = "Variable",
    y = "Relative Importance (%)",
    fill = "Model"
  ) +
  theme_publication() +
  theme(legend.position = "top")

ggsave(file.path(fig_dir, "variable_importance.png"),
       p_importance, width = 10, height = 8, dpi = 300)
cat("Saved: variable_importance.png\n\n")

# ============================================================================
# PART 7: PARTIAL DEPENDENCE PLOTS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 7: PARTIAL DEPENDENCE PLOTS\n")
cat("------------------------------------------------------------\n\n")

if (has_pdp) {
  cat("Creating partial dependence plots...\n")

  # Top 5 variables by importance
  top_vars <- head(variable_importance$variable, 5)

  # Wrapper function for ranger predictions
  pred_fun <- function(object, newdata) {
    predict(object, newdata)$predictions
  }

  pd_plots <- list()

  for (var in top_vars) {
    if (var %in% names(train_data)) {
      pd_data <- partial(
        rf_abundance,
        pred.var = var,
        train = train_data[, feature_cols_extended],
        pred.fun = pred_fun,
        grid.resolution = 20
      )

      pd_plots[[var]] <- autoplot(pd_data) +
        labs(title = paste("PD:", var)) +
        theme_publication()
    }
  }

  # Combine plots
  if (length(pd_plots) > 0) {
    p_pd_combined <- wrap_plots(pd_plots, ncol = 3)

    ggsave(file.path(fig_dir, "partial_dependence.png"),
           p_pd_combined, width = 14, height = 10, dpi = 300)
    cat("Saved: partial_dependence.png\n\n")
  }
} else {
  cat("Note: pdp package not available, skipping partial dependence plots\n\n")
}

# ============================================================================
# PART 8: SHAP VALUES (if available)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 8: SHAP VALUES\n")
cat("------------------------------------------------------------\n\n")

if (has_shapviz) {
  cat("Calculating SHAP values for XGBoost...\n")

  shap_obj <- shapviz(xgb_abundance, X_pred = xgb_train_matrix)

  # SHAP summary plot
  p_shap <- sv_importance(shap_obj, kind = "bee") +
    labs(title = "SHAP Values: XGBoost Abundance Model") +
    theme_publication()

  ggsave(file.path(fig_dir, "shap_summary.png"),
         p_shap, width = 12, height = 8, dpi = 300)
  cat("Saved: shap_summary.png\n\n")

  # SHAP interaction plots for top 2 variables
  top2 <- head(variable_importance$variable, 2)
  if (length(top2) == 2) {
    p_shap_interact <- sv_dependence(shap_obj, v = top2[1], color_var = top2[2]) +
      labs(
        title = paste("SHAP Dependence:", top2[1]),
        subtitle = paste("Colored by", top2[2]),
        x = top2[1],
        y = "SHAP Value"
      ) +
      theme_publication()

    ggsave(file.path(fig_dir, "shap_interaction.png"),
           p_shap_interact, width = 10, height = 8, dpi = 300)
    cat("Saved: shap_interaction.png\n\n")
  }
} else {
  cat("Note: shapviz package not available, skipping SHAP analysis\n")
  cat("Install with: install.packages('shapviz')\n\n")
}

# ============================================================================
# PART 9: PREDICTIONS VS ACTUAL PLOTS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 9: PREDICTIONS VS ACTUAL\n")
cat("------------------------------------------------------------\n\n")

# Abundance predictions
p_pred_abund <- test_data %>%
  pivot_longer(
    cols = c(glm_pred_abundance, rf_pred_abundance, xgb_pred_abundance),
    names_to = "model",
    values_to = "predicted"
  ) %>%
  mutate(
    model = case_when(
      model == "glm_pred_abundance" ~ "GLM (NB)",
      model == "rf_pred_abundance" ~ "Random Forest",
      model == "xgb_pred_abundance" ~ "XGBoost"
    )
  ) %>%
  ggplot(aes(x = cafi_abundance, y = predicted)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~model, scales = "free") +
  labs(
    title = "Predicted vs Actual: CAFI Abundance",
    x = "Actual Abundance",
    y = "Predicted Abundance"
  ) +
  theme_publication()

ggsave(file.path(fig_dir, "predictions_vs_actual_abundance.png"),
       p_pred_abund, width = 14, height = 5, dpi = 300)

# Richness predictions
p_pred_rich <- test_data %>%
  pivot_longer(
    cols = c(glm_pred_richness, rf_pred_richness, xgb_pred_richness),
    names_to = "model",
    values_to = "predicted"
  ) %>%
  mutate(
    model = case_when(
      model == "glm_pred_richness" ~ "GLM (Poisson)",
      model == "rf_pred_richness" ~ "Random Forest",
      model == "xgb_pred_richness" ~ "XGBoost"
    )
  ) %>%
  ggplot(aes(x = cafi_richness, y = predicted)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~model, scales = "free") +
  labs(
    title = "Predicted vs Actual: CAFI Richness",
    x = "Actual Richness",
    y = "Predicted Richness"
  ) +
  theme_publication()

ggsave(file.path(fig_dir, "predictions_vs_actual_richness.png"),
       p_pred_rich, width = 14, height = 5, dpi = 300)

cat("Saved: predictions_vs_actual_*.png\n\n")

# ============================================================================
# SAVE MODELS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("SAVING MODELS\n")
cat("------------------------------------------------------------\n\n")

# Save Random Forest models
saveRDS(rf_abundance, file.path(PATHS$objects, "rf_model_abundance.rds"))
saveRDS(rf_richness, file.path(PATHS$objects, "rf_model_richness.rds"))
cat("Saved: rf_model_abundance.rds\n")
cat("Saved: rf_model_richness.rds\n")

# Save XGBoost models
xgb.save(xgb_abundance, file.path(PATHS$objects, "xgb_model_abundance.xgb"))
saveRDS(xgb_abundance, file.path(PATHS$objects, "xgb_model_abundance.rds"))
cat("Saved: xgb_model_abundance.rds\n")

# Save all models in one object
ml_models <- list(
  glm_abundance = glm_abundance,
  glm_richness = glm_richness,
  glm_shannon = glm_shannon,
  rf_abundance = rf_abundance,
  rf_richness = rf_richness,
  rf_shannon = rf_shannon,
  xgb_abundance = xgb_abundance,
  xgb_richness = xgb_richness,
  xgb_shannon = xgb_shannon,
  feature_cols = feature_cols_extended,
  best_rf_params = best_rf_params,
  model_comparison = model_comparison,
  variable_importance = variable_importance,
  loso_results = loso_summary
)

saveRDS(ml_models, file.path(PATHS$objects, "machine_learning_models.rds"))
cat("Saved: machine_learning_models.rds\n\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("MACHINE LEARNING SUMMARY\n")
cat("========================================\n\n")

cat("Models Trained:\n")
cat("  - GLM: Negative Binomial (abundance), Poisson (richness), Gaussian (shannon)\n")
cat("  - Random Forest: ranger with tuned hyperparameters\n")
cat("  - XGBoost: with cross-validation tuning\n\n")

cat("Best Models by Target (test R2):\n")
best_by_target <- model_comparison %>%
  group_by(target) %>%
  slice_max(r2_test, n = 1) %>%
  ungroup()

for (i in 1:nrow(best_by_target)) {
  cat("  ", best_by_target$target[i], ": ", best_by_target$model[i],
      " (R2 = ", round(best_by_target$r2_test[i], 3), ")\n", sep = "")
}

cat("\nTop Variable Importance (abundance):\n")
for (i in 1:min(5, nrow(variable_importance))) {
  cat("  ", i, ". ", variable_importance$variable[i],
      " (", round(variable_importance$importance_mean[i], 1), "%)\n", sep = "")
}

cat("\nLeave-One-Site-Out CV (RF abundance):\n")
cat("  Mean RMSE:", round(mean(loso_summary$rmse), 2), "\n")
cat("  Mean R2:", round(mean(loso_summary$r2), 3), "\n")

cat("\nOutputs saved to:\n")
cat("  Models:", PATHS$objects, "\n")
cat("  Tables:", PATHS$tables, "\n")
cat("  Figures:", fig_dir, "\n\n")

cat("Machine learning analysis complete!\n")
