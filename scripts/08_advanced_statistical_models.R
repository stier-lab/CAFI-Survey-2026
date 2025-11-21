#!/usr/bin/env Rscript
# ============================================================================
# 08_advanced_statistical_models.R - Advanced statistical modeling of Survey data
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
# ============================================================================

cat("\n========================================\n")
cat("Advanced Statistical Modeling\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/Survey/00_load_libraries.R"))
library(glmmTMB)
library(DHARMa)
library(performance)
library(MuMIn)
library(effects)

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "advanced_models")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Data Preparation for Modeling
# ============================================================================

cat("Preparing data for modeling...\n")

# Create comprehensive modeling dataset
model_data <- metadata %>%
  left_join(
    cafi_clean %>%
      group_by(coral_id) %>%
      summarise(
        total_cafi = n(),
        species_richness = n_distinct(species),
        shannon = vegan::diversity(table(species)),
        mean_size = mean(size_mm, na.rm = TRUE),
        n_crabs = sum(type == "crab"),
        n_shrimps = sum(type == "shrimp"),
        n_fish = sum(type == "fish"),
        n_snails = sum(type == "snail"),
        .groups = "drop"
      ),
    by = "coral_id"
  ) %>%
  mutate(
    across(c(total_cafi:n_snails), ~replace_na(., 0)),
    log_cafi = log(total_cafi + 1),
    sqrt_cafi = sqrt(total_cafi),
    morphotype = factor(morphotype),
    site = factor(site),
    depth_scaled = scale(depth_m)[, 1]
  ) %>%
  {if("branch_width" %in% names(.)) mutate(., branch_width = factor(branch_width)) else .} %>%
  filter(!is.na(morphotype), !is.na(depth_m))

cat("✓ Model data prepared with", nrow(model_data), "observations\n\n")

# ============================================================================
# Generalized Linear Mixed Models (GLMMs)
# ============================================================================

cat("Fitting GLMMs for CAFI abundance...\n")

# Build formula based on available columns
if("branch_width" %in% names(model_data)) {
  base_formula <- total_cafi ~ morphotype + branch_width + depth_scaled + (1|site)
} else {
  base_formula <- total_cafi ~ morphotype + depth_scaled + (1|site)
}

# Negative Binomial GLMM for total CAFI
m1_nb <- glmmTMB(base_formula,
                 data = model_data,
                 family = nbinom2())

# Zero-inflated Negative Binomial
m2_zinb <- glmmTMB(base_formula,
                   data = model_data,
                   ziformula = ~1,
                   family = nbinom2())

# Poisson GLMM for comparison
m3_pois <- glmmTMB(base_formula,
                   data = model_data,
                   family = poisson())

# Model comparison
model_comparison <- compare_performance(m1_nb, m2_zinb, m3_pois,
                                       metrics = c("AIC", "BIC", "RMSE"))

write_csv(as.data.frame(model_comparison),
          file.path(SURVEY_TABLES, "glmm_model_comparison.csv"))

# Select best model
best_model <- m1_nb  # Based on AIC

# Model diagnostics
png(file.path(fig_dir, "model_diagnostics_abundance.png"),
    width = 12, height = 10, units = "in", res = 300)
par(mfrow = c(2, 2))
simulateResiduals(best_model) %>% plot()
dev.off()

# Extract coefficients
model_coefs <- broom.mixed::tidy(best_model, conf.int = TRUE)
write_csv(model_coefs,
          file.path(SURVEY_TABLES, "abundance_model_coefficients.csv"))

# Plot coefficients
p_coefs <- model_coefs %>%
  filter(term != "(Intercept)") %>%
  ggplot(aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(size = 3) +
  labs(title = "GLMM Coefficients for CAFI Abundance",
       x = "Coefficient Estimate",
       y = "Predictor")

ggsave(file.path(fig_dir, "glmm_coefficients_abundance.png"),
       p_coefs, width = 10, height = 6, dpi = 300)

cat("✓ Abundance models fitted\n\n")

# ============================================================================
# Species Richness Models
# ============================================================================

cat("Modeling species richness...\n")

# Build richness formula based on available columns
if("branch_width" %in% names(model_data)) {
  richness_formula <- species_richness ~ morphotype + branch_width + depth_scaled + offset(log(total_cafi + 1)) + (1|site)
} else {
  richness_formula <- species_richness ~ morphotype + depth_scaled + offset(log(total_cafi + 1)) + (1|site)
}

# GLM for species richness
m_richness <- glmmTMB(richness_formula,
                      data = model_data %>% filter(total_cafi > 0),
                      family = poisson())

# Model summary
richness_summary <- broom.mixed::tidy(m_richness, conf.int = TRUE)
write_csv(richness_summary,
          file.path(SURVEY_TABLES, "richness_model_summary.csv"))

# Predicted values
model_data$predicted_richness <- NA
model_data$predicted_richness[model_data$total_cafi > 0] <- predict(m_richness, type = "response")

# Plot predictions
p_richness_pred <- ggplot(model_data %>% filter(total_cafi > 0),
                          aes(x = depth_m, y = species_richness)) +
  geom_point(aes(color = morphotype), alpha = 0.6) +
  geom_line(aes(y = predicted_richness, color = morphotype), size = 1) +
  scale_color_viridis_d() +
  labs(title = "Species Richness Model Predictions",
       x = "Depth (m)",
       y = "Species Richness",
       color = "Morphotype") +
  facet_wrap(~morphotype)

ggsave(file.path(fig_dir, "richness_model_predictions.png"),
       p_richness_pred, width = 12, height = 8, dpi = 300)

cat("✓ Richness models fitted\n\n")

# ============================================================================
# Multivariate Models (Community Composition)
# ============================================================================

cat("Fitting multivariate community models...\n")

# Prepare community data for modeling
comm_hellinger <- decostand(community_matrix, method = "hellinger")

# db-RDA with environmental constraints
env_subset <- model_data %>%
  filter(coral_id %in% rownames(comm_hellinger)) %>%
  filter(complete.cases(morphotype, depth_scaled)) %>%  # Remove rows with NAs in predictors
  column_to_rownames("coral_id")

# Match order and ensure no NAs
comm_subset <- comm_hellinger[rownames(env_subset), ]
comm_subset <- comm_subset[complete.cases(comm_subset), ]  # Remove any rows with NAs in community data
env_subset <- env_subset[rownames(comm_subset), ]  # Match again after removing NAs

# Fit db-RDA (only if we have enough data left)
if(nrow(comm_subset) > 5) {
  if('branch_width' %in% names(env_subset)) {
    dbrda_model <- dbrda(comm_subset ~ morphotype + branch_width + depth_scaled, data = env_subset)
  } else {
    dbrda_model <- dbrda(comm_subset ~ morphotype + depth_scaled, data = env_subset)
  }
} else {
  cat("  Warning: Not enough complete cases for dbRDA\n")
  dbrda_model <- NULL
}

# ANOVA-like permutation test (only if model exists)
if(!is.null(dbrda_model)) {
  dbrda_anova <- anova(dbrda_model, by = "terms", permutations = 999)
  write.csv(as.data.frame(dbrda_anova),
            file.path(SURVEY_TABLES, "dbrda_permutation_test.csv"))

  # Variance partitioning
  if('branch_width' %in% names(env_subset)) {
    varpart_result <- varpart(comm_subset,
                              ~ morphotype,
                              ~ branch_width,
                              ~ depth_scaled,
                              data = env_subset)
  } else {
    varpart_result <- varpart(comm_subset,
                              ~ morphotype,
                              ~ depth_scaled,
                              data = env_subset)
  }
}

# Plot variance partitioning
png(file.path(fig_dir, "variance_partitioning.png"),
    width = 10, height = 8, units = "in", res = 300)
plot(varpart_result, Xnames = c("Morphotype", "Branch Width", "Depth"),
     main = "Variance Partitioning of Community Composition")
dev.off()

# Extract db-RDA scores
dbrda_scores <- scores(dbrda_model, display = c("sites", "species"))

# Plot db-RDA
p_dbrda <- as.data.frame(dbrda_scores$sites) %>%
  rownames_to_column("coral_id") %>%
  left_join(model_data, by = "coral_id") %>%
  ggplot(aes(x = dbRDA1, y = dbRDA2)) +
  geom_point(aes(color = morphotype, shape = site), size = 3, alpha = 0.7) +
  stat_ellipse(aes(color = morphotype), level = 0.95) +
  scale_color_viridis_d() +
  labs(title = "db-RDA Ordination",
       subtitle = paste("Constrained variance:",
                       round(sum(dbrda_model$CCA$eig) / sum(dbrda_model$tot.chi) * 100, 1), "%"),
       x = paste0("dbRDA1 (", round(dbrda_model$CCA$eig[1] / sum(dbrda_model$CCA$eig) * 100, 1), "%)"),
       y = paste0("dbRDA2 (", round(dbrda_model$CCA$eig[2] / sum(dbrda_model$CCA$eig) * 100, 1), "%)"))

ggsave(file.path(fig_dir, "dbrda_ordination.png"),
       p_dbrda, width = 10, height = 8, dpi = 300)

cat("✓ Multivariate models fitted\n\n")

# ============================================================================
# Hierarchical Models for Nested Data
# ============================================================================

cat("Fitting hierarchical models...\n")

# Three-level model: CAFI within corals within sites
if('branch_width' %in% names(model_data)) {
  m_hierarchical <- glmmTMB(total_cafi ~ morphotype * depth_scaled + branch_width +
                           (depth_scaled|site) + (1|coral_id),
                           data = model_data,
                           family = nbinom2())
} else {
  m_hierarchical <- glmmTMB(total_cafi ~ morphotype * depth_scaled +
                           (depth_scaled|site) + (1|coral_id),
                           data = model_data,
                           family = nbinom2())
}

# Extract random effects
random_effects <- ranef(m_hierarchical)

# Plot random effects
p_ranef <- as.data.frame(random_effects$cond$site) %>%
  rownames_to_column("site") %>%
  pivot_longer(cols = -site, names_to = "effect", values_to = "value") %>%
  ggplot(aes(x = site, y = value, fill = effect)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  labs(title = "Site-level Random Effects",
       x = "Site",
       y = "Random Effect Value",
       fill = "Effect Type")

ggsave(file.path(fig_dir, "hierarchical_random_effects.png"),
       p_ranef, width = 10, height = 6, dpi = 300)

# ICC calculation
icc_result <- performance::icc(m_hierarchical)
cat("  Intraclass Correlation Coefficient:", round(icc_result$ICC_adjusted, 3), "\n\n")

# ============================================================================
# Model Selection and Averaging
# ============================================================================

cat("Performing model selection...\n")

# Create candidate model set
if('branch_width' %in% names(model_data)) {
  m_full <- glmmTMB(total_cafi ~ morphotype + branch_width + depth_scaled +
                    morphotype:depth_scaled + (1|site),
                    data = model_data,
                    family = nbinom2(),
                    na.action = "na.fail")
} else {
  m_full <- glmmTMB(total_cafi ~ morphotype + depth_scaled +
                    morphotype:depth_scaled + (1|site),
                    data = model_data,
                    family = nbinom2(),
                    na.action = "na.fail")
}

# Dredge for best models
model_set <- dredge(m_full, rank = "AICc")

# Get top models (delta AICc < 2)
top_models <- subset(model_set, delta < 2)
write.csv(as.data.frame(top_models),
          file.path(SURVEY_TABLES, "model_selection_results.csv"))

# Model averaging
if (nrow(top_models) > 1) {
  avg_model <- model.avg(top_models)
  avg_coefs <- summary(avg_model)$coefmat.subset

  write.csv(avg_coefs,
            file.path(SURVEY_TABLES, "model_averaged_coefficients.csv"))

  cat("✓ Model averaging complete (", nrow(top_models), "models)\n\n")
}

# ============================================================================
# Non-linear Relationships (GAMs)
# ============================================================================

cat("Fitting GAMs for non-linear relationships...\n")

library(mgcv)

# GAM for abundance
if('branch_width' %in% names(model_data)) {
  gam_abundance <- gam(total_cafi ~ s(depth_m, k = 5) + morphotype + branch_width +
                      s(lat, long, k = 10) + s(site, bs = "re"),
                      data = model_data,
                      family = nb())
} else {
  gam_abundance <- gam(total_cafi ~ s(depth_m, k = 5) + morphotype +
                      s(lat, long, k = 10) + s(site, bs = "re"),
                      data = model_data,
                      family = nb())
}

# GAM for richness
gam_richness <- gam(species_richness ~ s(depth_m, k = 5) + morphotype +
                   s(total_cafi, k = 5) + s(site, bs = "re"),
                   data = model_data %>% filter(total_cafi > 0),
                   family = poisson())

# Model summaries
gam_summary <- list(
  abundance = broom::tidy(gam_abundance),
  richness = broom::tidy(gam_richness)
)

saveRDS(gam_summary, file.path(SURVEY_OBJECTS, "gam_model_summaries.rds"))

# Plot smooth terms
png(file.path(fig_dir, "gam_smooth_terms.png"),
    width = 12, height = 10, units = "in", res = 300)
par(mfrow = c(2, 3))
plot(gam_abundance, pages = 1, main = "GAM: Abundance")
dev.off()

# Prediction surfaces
if (sum(!is.na(model_data$lat) & !is.na(model_data$long)) > 20) {
  # Create prediction grid
  lat_range <- range(model_data$lat, na.rm = TRUE)
  long_range <- range(model_data$long, na.rm = TRUE)

  pred_grid <- expand.grid(
    lat = seq(lat_range[1], lat_range[2], length.out = 50),
    long = seq(long_range[1], long_range[2], length.out = 50),
    depth_m = mean(model_data$depth_m, na.rm = TRUE),
    morphotype = levels(model_data$morphotype)[1],
    # branch_width omitted (not available in data),
    site = levels(model_data$site)[1]
  )

  pred_grid$predicted <- predict(gam_abundance, pred_grid, type = "response")

  # Plot prediction surface
  p_spatial_pred <- ggplot(pred_grid, aes(x = long, y = lat, fill = predicted)) +
    geom_tile() +
    geom_point(data = model_data, aes(x = long, y = lat), size = 1, alpha = 0.5) +
    scale_fill_viridis_c(trans = "sqrt") +
    labs(title = "Spatial Prediction of CAFI Abundance",
         x = "Longitude",
         y = "Latitude",
         fill = "Predicted\nAbundance") +
    coord_quickmap()

  ggsave(file.path(fig_dir, "spatial_abundance_prediction.png"),
         p_spatial_pred, width = 10, height = 8, dpi = 300)
}

cat("✓ GAM models fitted\n\n")

# ============================================================================
# Effect Sizes and Predictions
# ============================================================================

cat("Calculating effect sizes...\n")

# Calculate marginal effects
if (exists("best_model")) {
  # Create prediction dataset
  pred_data <- expand.grid(
    morphotype = levels(model_data$morphotype),
    # branch_width omitted (not available in data),
    depth_scaled = seq(-2, 2, length.out = 50),
    site = levels(model_data$site)[1]
  )

  # Predictions with confidence intervals
  preds <- predict(best_model, pred_data,
                   type = "response",
                   se.fit = TRUE,
                   re.form = NA)  # Population-level predictions

  pred_data$predicted <- preds$fit
  pred_data$se <- preds$se.fit
  pred_data$lower <- pred_data$predicted - 1.96 * pred_data$se
  pred_data$upper <- pred_data$predicted + 1.96 * pred_data$se

  # Unscale depth
  pred_data$depth_m <- pred_data$depth_scaled * sd(model_data$depth_m, na.rm = TRUE) +
                       mean(model_data$depth_m, na.rm = TRUE)

  # Plot marginal effects
  p_marginal <- ggplot(pred_data, aes(x = depth_m, y = predicted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = morphotype), alpha = 0.3) +
    geom_line(aes(color = morphotype), size = 1.5) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    labs(title = "Marginal Effects of Depth on CAFI Abundance",
         x = "Depth (m)",
         y = "Predicted CAFI Abundance",
         color = "Morphotype",
         fill = "Morphotype")

  ggsave(file.path(fig_dir, "marginal_effects_depth.png"),
         p_marginal, width = 10, height = 6, dpi = 300)
}

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Advanced Statistical Modeling Summary\n")
cat("========================================\n\n")

cat("Models Fitted:\n")
cat("  - GLMMs for abundance (3 variants)\n")
cat("  - GLMMs for species richness\n")
cat("  - db-RDA for community composition\n")
cat("  - Hierarchical models with random slopes\n")
cat("  - GAMs for non-linear relationships\n\n")

if (exists("model_comparison")) {
  cat("Best Model (by AIC):\n")
  best_idx <- which.min(model_comparison$AIC)
  cat("  - Model:", rownames(model_comparison)[best_idx], "\n")
  cat("  - AIC:", model_comparison$AIC[best_idx], "\n\n")
}

if (exists("dbrda_model")) {
  cat("Community Composition:\n")
  cat("  - Constrained variance:",
      round(sum(dbrda_model$CCA$eig) / sum(dbrda_model$tot.chi) * 100, 1), "%\n")
}

if (exists("icc_result")) {
  cat("  - ICC (site effect):", round(icc_result$ICC_adjusted, 3), "\n\n")
}

cat("✅ Advanced statistical modeling complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", SURVEY_TABLES, "\n")
