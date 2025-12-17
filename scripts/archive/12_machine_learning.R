#!/usr/bin/env Rscript
# ============================================================================
# 12_machine_learning.R - Machine Learning Predictions
# ============================================================================
#
# PURPOSE: Use ML models to predict CAFI abundance from coral/landscape features
#
# ANALYSES:
#   - Random Forest variable importance
#   - Gradient boosting (XGBoost)
#   - Cross-validation performance
#   - Partial dependence plots
#
# MANUSCRIPT FIGURES: None (supplementary)
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("12: Machine Learning Analysis\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Additional ML packages (required for this script)
if (requireNamespace("caret", quietly = TRUE)) {
  suppressPackageStartupMessages(library(caret))
} else {
  cat("  Note: caret package not installed, ML analysis will be limited\n")
}
if (requireNamespace("randomForest", quietly = TRUE)) {
  suppressPackageStartupMessages(library(randomForest))
}
if (requireNamespace("ranger", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ranger))
}

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")
community_matrix <- load_object("community_matrix.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "machine_learning")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Feature Engineering
# ============================================================================

cat("Engineering features for ML models...\n")

# Create ML dataset with engineered features
# coral_master already has total_cafi, otu_richness, etc., so add only new features

cafi_features <- cafi_clean %>%
  group_by(coral_id) %>%
  summarise(
    # Diversity metrics (may differ from coral_master)
    shannon_div = vegan::diversity(table(species)),
    simpson = vegan::diversity(table(species), index = "simpson"),

    # Size metrics
    mean_size = mean(size_mm, na.rm = TRUE),
    sd_size = sd(size_mm, na.rm = TRUE),
    max_size = if(any(!is.na(size_mm))) max(size_mm, na.rm = TRUE) else NA_real_,
    min_size = if(any(!is.na(size_mm))) min(size_mm, na.rm = TRUE) else NA_real_,

    # Type proportions
    prop_crabs = sum(type == "crab") / n(),
    prop_shrimps = sum(type == "shrimp") / n(),
    prop_fish = sum(type == "fish") / n(),
    prop_snails = sum(type == "snail") / n(),

    # Dominant species
    dominant_prop = if(n() > 0) max(table(species)) / n() else NA_real_,

    .groups = "drop"
  )

ml_data <- coral_master %>%
  left_join(cafi_features, by = "coral_id") %>%
  mutate(
    # Handle missing values
    across(where(is.numeric), ~replace_na(., 0)),

    # Ensure key columns exist
    total_cafi = replace_na(total_cafi, 0),
    shannon_div = replace_na(shannon_div, 0),

    # Create categorical targets
    high_diversity = factor(ifelse(shannon_div > median(shannon_div, na.rm = TRUE),
                                   "High", "Low")),
    high_abundance = factor(ifelse(total_cafi > median(total_cafi, na.rm = TRUE),
                                  "High", "Low")),

    # Encode factors as numeric for some models
    morphotype_num = as.numeric(factor(morphotype)),
    site_num = as.numeric(factor(site))
  ) %>%
  # Add branch_width_num only if branch_width exists
  {if("branch_width" %in% names(.)) mutate(., branch_width_num = as.numeric(factor(branch_width))) else .} %>%
  filter(!is.na(morphotype), !is.na(depth_m))

# Split data
set.seed(123)
train_idx <- createDataPartition(ml_data$high_diversity, p = 0.7, list = FALSE)
train_data <- ml_data[train_idx, ]
test_data <- ml_data[-train_idx, ]

cat("✓ Features engineered:", ncol(ml_data), "features\n")
cat("  - Training set:", nrow(train_data), "samples\n")
cat("  - Test set:", nrow(test_data), "samples\n")
cat("  - Available columns:", paste(head(names(ml_data), 15), collapse = ", "), "...\n\n")

# ============================================================================
# Random Forest Classification
# ============================================================================

cat("Training Random Forest models...\n")

# Prepare features (only include columns that exist in ml_data)
potential_features <- c("depth_m", "morphotype_num", "branch_width_num", "lat", "long", "site_num")
feature_cols <- potential_features[potential_features %in% names(ml_data)]
cat("  Using features:", paste(feature_cols, collapse = ", "), "\n")

# Build dynamic formulas
diversity_formula <- as.formula(paste("high_diversity ~", paste(feature_cols, collapse = " + ")))
abundance_formula <- as.formula(paste("total_cafi ~", paste(feature_cols, collapse = " + ")))

# Classification for high diversity
rf_diversity <- randomForest(
  diversity_formula,
  data = train_data,
  ntree = 500,
  importance = TRUE,
  na.action = na.omit
)

# Variable importance
importance_df <- as.data.frame(randomForest::importance(rf_diversity)) %>%
  rownames_to_column("variable") %>%
  arrange(desc(MeanDecreaseGini))

write_csv(importance_df,
          file.path(TABLES_DIR, "rf_variable_importance.csv"))

# Plot importance
p_importance <- importance_df %>%
  ggplot(aes(x = reorder(variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "#0072B2") +
  coord_flip() +
  labs(title = "Random forest variable importance",
       x = "Variable",
       y = "Mean decrease in Gini index") +
  theme_publication()

ggsave(file.path(fig_dir, "rf_variable_importance.png"),
       p_importance, width = 10, height = 6, dpi = 300)

# Predictions
test_data$rf_pred_diversity <- predict(rf_diversity, test_data)

# Confusion matrix
conf_matrix <- confusionMatrix(test_data$rf_pred_diversity, test_data$high_diversity)
write.csv(as.matrix(conf_matrix$table),
          file.path(TABLES_DIR, "rf_confusion_matrix.csv"))

cat("  RF Accuracy:", round(conf_matrix$overall["Accuracy"], 3), "\n\n")

# ============================================================================
# Regression: Predicting CAFI Abundance
# ============================================================================

cat("Training regression models for abundance...\n")

# Random Forest regression
rf_abundance <- ranger(
  abundance_formula,
  data = train_data,
  num.trees = 500,
  importance = "impurity"
)

# Predictions
test_data$rf_pred_abundance <- predict(rf_abundance, test_data)$predictions

# Calculate performance metrics
rf_rmse <- sqrt(mean((test_data$total_cafi - test_data$rf_pred_abundance)^2, na.rm = TRUE))
rf_mae <- mean(abs(test_data$total_cafi - test_data$rf_pred_abundance), na.rm = TRUE)
rf_r2 <- cor(test_data$total_cafi, test_data$rf_pred_abundance, use = "complete.obs")^2

cat("  RF Regression - RMSE:", round(rf_rmse, 2), "\n")
cat("  RF Regression - MAE:", round(rf_mae, 2), "\n")
cat("  RF Regression - R²:", round(rf_r2, 3), "\n\n")

# Plot predictions vs actual
p_pred_actual <- ggplot(test_data, aes(x = total_cafi, y = rf_pred_abundance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, linewidth = 0.8) +
  labs(title = "Random forest predictions vs actual",
       subtitle = paste("R² =", round(rf_r2, 3), "| RMSE =", round(rf_rmse, 2)),
       x = "Actual CAFI abundance",
       y = "Predicted CAFI abundance") +
  theme_publication()

ggsave(file.path(fig_dir, "rf_predictions_vs_actual.png"),
       p_pred_actual, width = 10, height = 8, dpi = 300)

# ============================================================================
# XGBoost Models
# ============================================================================

cat("Training XGBoost models...\n")

# XGBoost models (if package available)
if (requireNamespace("xgboost", quietly = TRUE)) {
  library(xgboost)

  # Prepare data for XGBoost
  xgb_features <- train_data %>%
    select(any_of(feature_cols)) %>%
    as.matrix()

  xgb_labels <- as.numeric(train_data$high_diversity) - 1  # Convert to 0/1

  # Train XGBoost
  xgb_model <- xgboost(
    data = xgb_features,
    label = xgb_labels,
    nrounds = 100,
    objective = "binary:logistic",
    eval_metric = "logloss",
    verbose = 0
  )

  # Test predictions
  xgb_test_features <- test_data %>%
    select(any_of(feature_cols)) %>%
    as.matrix()

  test_data$xgb_prob <- predict(xgb_model, xgb_test_features)
  test_data$xgb_pred <- factor(ifelse(test_data$xgb_prob > 0.5, "High", "Low"),
                               levels = c("High", "Low"))

  # XGBoost accuracy
  xgb_accuracy <- mean(test_data$xgb_pred == test_data$high_diversity, na.rm = TRUE)
  cat("  XGBoost Accuracy:", round(xgb_accuracy, 3), "\n\n")

  # Feature importance
  xgb_importance <- xgb.importance(
    feature_names = colnames(xgb_features),
    model = xgb_model
  )

  # Plot XGBoost importance
  p_xgb_imp <- xgb_importance %>%
    ggplot(aes(x = reorder(Feature, Gain), y = Gain)) +
    geom_col(fill = "#009E73") +
    coord_flip() +
    labs(title = "XGBoost feature importance",
         x = "Feature",
         y = "Gain") +
    theme_publication()

  ggsave(file.path(fig_dir, "xgboost_feature_importance.png"),
         p_xgb_imp, width = 10, height = 6, dpi = 300)
} else {
  cat("  Skipping XGBoost (package not installed)\n\n")
}

# ============================================================================
# Species Distribution Models
# ============================================================================

cat("Building species distribution models...\n")

# Select top species for modeling
top_species <- cafi_clean %>%
  count(species) %>%
  top_n(10, n) %>%
  pull(species)

sdm_results <- list()

for (sp in top_species) {
  # Create presence/absence data
  sp_data <- ml_data %>%
    left_join(
      cafi_clean %>%
        filter(species == sp) %>%
        distinct(coral_id) %>%
        mutate(present = 1),
      by = "coral_id"
    ) %>%
    mutate(present = replace_na(present, 0))

  # Split data
  sp_train <- sp_data[train_idx, ]
  sp_test <- sp_data[-train_idx, ]

  # Train model (use branch_width only if it exists)
  if("branch_width" %in% names(sp_train)) {
    sp_model <- glm(present ~ depth_m + morphotype + branch_width,
                    data = sp_train,
                    family = binomial())
  } else {
    sp_model <- glm(present ~ depth_m + morphotype,
                    data = sp_train,
                    family = binomial())
  }

  # Predictions
  sp_test$predicted_prob <- predict(sp_model, sp_test, type = "response")

  # Calculate AUC (if ROCR package available)
  if (requireNamespace("ROCR", quietly = TRUE) &&
      sum(sp_test$present) > 0 && sum(sp_test$present) < nrow(sp_test)) {
    pred_obj <- ROCR::prediction(sp_test$predicted_prob, sp_test$present)
    auc <- ROCR::performance(pred_obj, "auc")@y.values[[1]]

    sdm_results[[sp]] <- list(
      species = sp,
      auc = auc,
      prevalence = mean(sp_data$present),
      model = sp_model
    )
  }
}

# Summarize SDM results (only if ROCR was available and models were fitted)
if (length(sdm_results) > 0) {
  sdm_summary <- bind_rows(lapply(sdm_results, function(x) {
    data.frame(species = x$species, auc = x$auc, prevalence = x$prevalence)
  })) %>%
    arrange(desc(auc))

  write_csv(sdm_summary,
            file.path(TABLES_DIR, "species_distribution_model_performance.csv"))

  # Plot SDM performance
  p_sdm <- sdm_summary %>%
    ggplot(aes(x = reorder(species, auc), y = auc)) +
    geom_col(aes(fill = prevalence)) +
    geom_hline(yintercept = 0.7, linetype = "dashed", alpha = 0.5) +
    scale_fill_viridis_c() +
    coord_flip() +
    labs(title = "Species distribution model performance",
         x = "Species",
         y = "AUC",
         fill = "Prevalence") +
    theme_publication()

  ggsave(file.path(fig_dir, "sdm_performance.png"),
         p_sdm, width = 10, height = 8, dpi = 300)

  cat("✓ SDMs fitted for", length(sdm_results), "species\n\n")
} else {
  cat("  Note: SDMs skipped (ROCR package not installed)\n\n")
}

# ============================================================================
# Cross-validation and Model Comparison
# ============================================================================

cat("Performing cross-validation...\n")

# Set up cross-validation
ctrl <- trainControl(
  method = "cv",
  number = 10,
  savePredictions = TRUE,
  classProbs = TRUE
)

# Train multiple models
models <- list()

# Build CV formula based on available columns
cv_features <- intersect(c("depth_m", "morphotype_num", "branch_width_num", "site_num"),
                         names(train_data))
cv_formula <- as.formula(paste("high_diversity ~", paste(cv_features, collapse = " + ")))

# Logistic regression
models$logistic <- train(
  cv_formula,
  data = train_data,
  method = "glm",
  family = "binomial",
  trControl = ctrl
)

# Random Forest
models$rf <- train(
  cv_formula,
  data = train_data,
  method = "rf",
  trControl = ctrl,
  ntree = 100
)

# Compare models
model_comparison <- resamples(models)
comparison_summary <- summary(model_comparison)

# Plot model comparison
png(file.path(fig_dir, "model_comparison_boxplot.png"),
    width = 10, height = 8, units = "in", res = 300)
bwplot(model_comparison)
dev.off()

cat("✓ Cross-validation complete\n\n")

# ============================================================================
# Partial Dependence Plots
# ============================================================================

cat("Creating partial dependence plots...\n")

if (requireNamespace("pdp", quietly = TRUE) && exists("rf_abundance")) {
  # Partial dependence for depth
  pd_depth <- pdp::partial(rf_abundance, pred.var = "depth_m",
                     train = train_data[, feature_cols])

  p_pd_depth <- pdp::autoplot(pd_depth) +
    labs(title = "Partial dependence: depth effect on CAFI abundance",
         x = "Depth (m)",
         y = "Partial dependence") +
    theme_publication()

  ggsave(file.path(fig_dir, "partial_dependence_depth.png"),
         p_pd_depth, width = 10, height = 6, dpi = 300)
  cat("✓ Partial dependence plots created\n\n")
} else {
  cat("  Note: Partial dependence plots skipped (pdp package not installed)\n\n")
}

# ============================================================================
# Ensemble Predictions
# ============================================================================

cat("Creating ensemble predictions...\n")

# Combine predictions from multiple models
test_data$ensemble_prob <- (test_data$rf_pred_diversity == "High") * 0.5 +
                           test_data$xgb_prob * 0.5

test_data$ensemble_pred <- factor(ifelse(test_data$ensemble_prob > 0.5,
                                         "High", "Low"),
                                  levels = c("High", "Low"))

# Ensemble accuracy
ensemble_accuracy <- mean(test_data$ensemble_pred == test_data$high_diversity,
                         na.rm = TRUE)

cat("  Ensemble Accuracy:", round(ensemble_accuracy, 3), "\n\n")

# ============================================================================
# Save ML Models
# ============================================================================

ml_models <- list(
  rf_diversity = rf_diversity,
  rf_abundance = rf_abundance,
  xgb_model = xgb_model,
  sdm_models = sdm_results
)

saveRDS(ml_models, file.path(OBJECTS_DIR, "machine_learning_models.rds"))

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Machine Learning Summary\n")
cat("========================================\n\n")

cat("Models Trained:\n")
cat("  - Random Forest (classification & regression)\n")
cat("  - XGBoost (classification)\n")
cat("  - Species Distribution Models (", length(sdm_results), "species)\n")
cat("  - Ensemble model\n\n")

cat("Performance Metrics:\n")
cat("  Classification (Diversity):\n")
cat("    - Random Forest:", round(conf_matrix$overall["Accuracy"], 3), "\n")
cat("    - XGBoost:", round(xgb_accuracy, 3), "\n")
cat("    - Ensemble:", round(ensemble_accuracy, 3), "\n\n")

cat("  Regression (Abundance):\n")
cat("    - RMSE:", round(rf_rmse, 2), "\n")
cat("    - R²:", round(rf_r2, 3), "\n\n")

if (exists("sdm_summary")) {
  cat("  Best SDM Performance:\n")
  cat("    - Species:", sdm_summary$species[1], "\n")
  cat("    - AUC:", round(sdm_summary$auc[1], 3), "\n\n")
}

# ============================================================================
# STATISTICAL RESULTS SUMMARY
# ============================================================================

cat("Compiling statistical results...\n")

stats_results <- init_results_df()

# 1. Random Forest Classification (Diversity)
if(exists("conf_matrix")) {
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "ML_RF_class",
      question = "Can Random Forest classify high/low diversity corals?",
      test_name = "Random Forest Classification",
      test_statistic = conf_matrix$overall["Accuracy"],
      df = paste(nrow(train_data), "train,", nrow(test_data), "test"),
      p_value = NA,
      effect_size = conf_matrix$overall["Kappa"],
      effect_type = "Kappa",
      n = nrow(ml_data),
      notes = sprintf("Sensitivity=%.3f, Specificity=%.3f",
                      conf_matrix$byClass["Sensitivity"],
                      conf_matrix$byClass["Specificity"])
    )
  )
}

# 2. Random Forest Regression (Abundance)
if(exists("rf_rmse")) {
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "ML_RF_reg",
      question = "Can Random Forest predict CAFI abundance?",
      test_name = "Random Forest Regression",
      test_statistic = rf_rmse,
      df = paste(nrow(train_data), "train,", nrow(test_data), "test"),
      p_value = NA,
      effect_size = rf_r2,
      effect_type = "R²",
      n = nrow(ml_data),
      notes = sprintf("RMSE=%.2f, MAE=%.2f", rf_rmse, rf_mae)
    )
  )
}

# 3. XGBoost Classification
if(exists("xgb_accuracy")) {
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "ML_XGB",
      question = "Can XGBoost classify high/low diversity corals?",
      test_name = "XGBoost Classification",
      test_statistic = xgb_accuracy,
      df = "100 rounds",
      p_value = NA,
      effect_size = xgb_accuracy,
      effect_type = "Accuracy",
      n = nrow(ml_data),
      notes = "Binary classification for diversity"
    )
  )
}

# 4. Ensemble Model
if(exists("ensemble_accuracy")) {
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "ML_Ensemble",
      question = "Does ensemble improve over individual models?",
      test_name = "Ensemble (RF + XGBoost)",
      test_statistic = ensemble_accuracy,
      df = "50/50 weight",
      p_value = NA,
      effect_size = ensemble_accuracy,
      effect_type = "Accuracy",
      n = nrow(test_data),
      notes = "Average of RF and XGBoost predictions"
    )
  )
}

# 5. SDM results for top species
if(exists("sdm_summary") && nrow(sdm_summary) > 0) {
  for(i in 1:min(5, nrow(sdm_summary))) {
    row <- sdm_summary[i, ]
    stats_results <- bind_rows(stats_results,
      create_result_row(
        hypothesis = paste0("SDM_", i),
        question = paste("Can we predict", row$species, "distribution?"),
        test_name = "Species Distribution Model (GLM)",
        test_statistic = row$auc,
        df = NA,
        p_value = NA,
        effect_size = row$auc,
        effect_type = "AUC",
        n = nrow(ml_data),
        notes = sprintf("Prevalence=%.3f", row$prevalence)
      )
    )
  }
}

# 6. Variable importance (descriptive)
if(exists("importance_df") && nrow(importance_df) > 0) {
  top_var <- importance_df$variable[1]
  top_gini <- importance_df$MeanDecreaseGini[1]

  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "ML_VarImp",
      question = "What is the most important predictor?",
      test_name = "RF Variable Importance (Gini)",
      test_statistic = top_gini,
      df = ncol(xgb_features),
      p_value = NA,
      effect_size = top_gini / sum(importance_df$MeanDecreaseGini),
      effect_type = "Prop. of total Gini",
      n = nrow(train_data),
      notes = sprintf("Top variable: %s", top_var)
    )
  )
}

# Save statistical results
save_stats_summary(stats_results, "12_machine_learning", "Machine Learning Analysis")

cat("✅ Machine learning analysis complete!\n")
cat("Models saved to:", file.path(OBJECTS_DIR, "machine_learning_models.rds"), "\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", TABLES_DIR, "\n")