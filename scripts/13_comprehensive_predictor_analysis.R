#!/usr/bin/env Rscript
# ============================================================================
# 13_comprehensive_predictor_analysis.R - Analysis of all key predictors including
# coral size, neighbor distance, spatial location, and morphology
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
# ============================================================================

cat("\n========================================\n")
cat("Comprehensive Predictor Analysis\n")
cat("========================================\n\n")

# Load libraries
source(here::here("scripts/Survey/00_load_libraries.R"))
library(glmmTMB)
library(MuMIn)
library(performance)
library(car)
library(corrplot)

# Load data
cat("Loading comprehensive coral characteristics data...\n")

# Load coral characteristics with all predictors
coral_chars <- read_csv(file.path(here("data/Survey"),
                                  "1. survey_coral_characteristics_merged_v2.csv"))

# Load CAFI data
cafi_data <- read_csv(file.path(here("data/Survey"),
                                "1. survey_cafi_data_w_taxonomy_summer2019_v5.csv"))

# Load physiology data if available
physio_data <- read_csv(file.path(here("data/Survey"),
                                  "1. survey_master_phys_data_v3.csv"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "comprehensive_predictors")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Data Preparation - All Predictors
# ============================================================================

cat("Preparing comprehensive predictor dataset...\n")

# Calculate CAFI metrics per coral
# First clean column names
cafi_data <- cafi_data %>% janitor::clean_names()

cafi_metrics <- cafi_data %>%
  group_by(coral_id) %>%
  summarise(
    total_cafi = n(),
    species_richness = n_distinct(if_else(!is.na(lowest_level), lowest_level,
                                          if_else(!is.na(search_term), search_term, paste(type, code)))),
    shannon = if(n() > 0) vegan::diversity(table(if_else(!is.na(lowest_level), lowest_level,
                                             if_else(!is.na(search_term), search_term, paste(type, code))))) else 0,
    simpson = if(n() > 0) vegan::diversity(table(if_else(!is.na(lowest_level), lowest_level,
                                             if_else(!is.na(search_term), search_term, paste(type, code)))), index = "simpson") else 0,

    # Type-specific counts
    n_crabs = sum(type == "crab", na.rm = TRUE),
    n_shrimps = sum(type == "shrimp", na.rm = TRUE),
    n_fish = sum(type == "fish", na.rm = TRUE),
    n_snails = sum(type == "snail", na.rm = TRUE),

    # Size metrics (using cafi_size_mm from raw data)
    mean_cafi_size = mean(cafi_size_mm, na.rm = TRUE),
    max_cafi_size = max(cafi_size_mm, na.rm = TRUE),

    .groups = "drop"
  )

# Merge all data
comprehensive_data <- coral_chars %>%
  left_join(cafi_metrics, by = "coral_id") %>%
  left_join(physio_data %>% select(coral_id, any_of(c(names(physio_data)[grepl("^zoox|^chl", names(physio_data))]))),
            by = "coral_id") %>%
  mutate(
    # Replace NAs in CAFI metrics with 0
    across(c(total_cafi:max_cafi_size), ~replace_na(., 0)),

    # Calculate coral size metrics
    coral_volume = coalesce(volume_field, volume_lab),
    coral_surface_area = 2 * pi * (width_field/2) * height_field,  # Cylinder approximation
    coral_complexity = height_field / width_field,  # Height:width ratio

    # Calculate neighbor metrics
    neighbor_density = number_of_neighbors / mean_neighbor_distance,
    isolation_index = mean_neighbor_distance / coral_volume^(1/3),

    # Categorize branch width
    branch_category = factor(branch_width, levels = c("tight", "wide")),

    # Scale continuous predictors
    depth_scaled = scale(depth)[,1],
    volume_scaled = scale(log(coral_volume + 1))[,1],
    neighbor_dist_scaled = scale(log(mean_neighbor_distance + 1))[,1],

    # Create interaction terms
    volume_x_morphotype = coral_volume * as.numeric(factor(morphotype)),
    depth_x_morphotype = depth * as.numeric(factor(morphotype)),

    # Spatial coordinates
    lat_scaled = scale(lat)[,1],
    long_scaled = scale(long)[,1]
  ) %>%
  filter(!is.na(coral_volume), !is.na(total_cafi))

cat("✓ Dataset prepared with", nrow(comprehensive_data), "corals and",
    ncol(comprehensive_data), "variables\n\n")

# ============================================================================
# Exploratory Analysis - Key Predictors
# ============================================================================

cat("Analyzing key predictor relationships...\n")

# 1. CORAL SIZE EFFECTS
# --------------------
p_size_abundance <- ggplot(comprehensive_data,
                           aes(x = coral_volume, y = total_cafi)) +
  geom_point(aes(color = morphotype), alpha = 0.6) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE) +
  scale_x_log10() +
  scale_y_sqrt() +
  labs(title = "CAFI Abundance vs Coral Volume",
       x = "Coral Volume (cm³, log scale)",
       y = "Total CAFI (sqrt scale)",
       color = "Morphotype")

p_size_richness <- ggplot(comprehensive_data,
                          aes(x = coral_volume, y = species_richness)) +
  geom_point(aes(color = morphotype), alpha = 0.6) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE) +
  scale_x_log10() +
  labs(title = "Species Richness vs Coral Volume",
       x = "Coral Volume (cm³, log scale)",
       y = "Species Richness",
       color = "Morphotype")

# Combine size plots
size_plots <- p_size_abundance / p_size_richness +
  plot_annotation(title = "Coral Size Effects on CAFI Communities")

ggsave(file.path(fig_dir, "coral_size_effects.png"),
       size_plots, width = 12, height = 10, dpi = 300)

# 2. NEIGHBOR DISTANCE EFFECTS
# ----------------------------
p_neighbor_dist <- ggplot(comprehensive_data %>% filter(mean_neighbor_distance > 0),
                         aes(x = mean_neighbor_distance, y = total_cafi)) +
  geom_point(aes(color = number_of_neighbors), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_log10() +
  scale_color_viridis_c() +
  labs(title = "CAFI Abundance vs Neighbor Distance",
       x = "Mean Distance to Neighbors (cm, log scale)",
       y = "Total CAFI",
       color = "Number of\nNeighbors")

p_neighbor_density <- ggplot(comprehensive_data %>% filter(number_of_neighbors > 0),
                            aes(x = number_of_neighbors, y = total_cafi)) +
  geom_boxplot(aes(group = cut_number(number_of_neighbors, 5)), alpha = 0.7) +
  geom_smooth(method = "gam", se = TRUE) +
  labs(title = "CAFI vs Neighbor Density",
       x = "Number of Neighboring Corals",
       y = "Total CAFI")

neighbor_plots <- p_neighbor_dist | p_neighbor_density

ggsave(file.path(fig_dir, "neighbor_effects.png"),
       neighbor_plots, width = 14, height = 6, dpi = 300)

# 3. DEPTH GRADIENT
# -----------------
p_depth <- ggplot(comprehensive_data, aes(x = depth, y = total_cafi)) +
  geom_point(aes(color = morphotype), alpha = 0.6) +
  geom_smooth(aes(color = morphotype), method = "gam", se = TRUE) +
  facet_wrap(~site, scales = "free_x") +
  labs(title = "Depth Effects on CAFI Abundance by Site",
       x = "Depth (m)",
       y = "Total CAFI",
       color = "Morphotype")

ggsave(file.path(fig_dir, "depth_effects_by_site.png"),
       p_depth, width = 14, height = 8, dpi = 300)

cat("✓ Key predictor visualizations created\n\n")

# ============================================================================
# Correlation Analysis - All Predictors
# ============================================================================

cat("Performing correlation analysis of predictors...\n")

# Select numeric predictors
numeric_predictors <- comprehensive_data %>%
  select(coral_volume, coral_surface_area, coral_complexity,
         depth, number_of_neighbors, mean_neighbor_distance,
         neighbor_density, isolation_index,
         mean_total_volume_of_neighbors,
         total_cafi, species_richness, shannon) %>%
  na.omit()

# Calculate correlation matrix
cor_matrix <- cor(numeric_predictors, use = "complete.obs")

# Plot correlation matrix
png(file.path(fig_dir, "predictor_correlation_matrix.png"),
    width = 12, height = 10, units = "in", res = 300)
corrplot(cor_matrix, method = "color", type = "upper",
         order = "hclust", addrect = 3,
         col = colorRampPalette(c("blue", "white", "red"))(100),
         tl.col = "black", tl.srt = 45)
dev.off()

# Identify highly correlated predictors
if (requireNamespace("caret", quietly = TRUE)) {
  high_cor <- caret::findCorrelation(cor_matrix, cutoff = 0.7)
  if (length(high_cor) > 0) {
    cat("  Highly correlated predictors (r > 0.7):\n")
    cat("   ", names(numeric_predictors)[high_cor], "\n\n")
  }
} else {
  cat("  Note: Correlation filtering skipped (caret package not installed)\n\n")
}

# ============================================================================
# Comprehensive Statistical Models
# ============================================================================

cat("Fitting comprehensive statistical models...\n")

# 1. FULL MODEL WITH ALL KEY PREDICTORS
# --------------------------------------
cat("  Fitting full model with all predictors...\n")

# Full model for abundance
m_full <- glmmTMB(
  total_cafi ~
    # Coral characteristics
    scale(log(coral_volume + 1)) +
    morphotype +
    branch_category +

    # Spatial/environmental
    poly(depth, 2) +
    scale(lat) + scale(long) +

    # Neighbor effects
    scale(log(mean_neighbor_distance + 1)) +
    scale(number_of_neighbors) +

    # Interactions
    scale(log(coral_volume + 1)):morphotype +
    poly(depth, 2):morphotype +

    # Random effects
    (1|site),

  data = comprehensive_data,
  family = nbinom2(),
  na.action = "na.fail"
)

# Model diagnostics
png(file.path(fig_dir, "full_model_diagnostics.png"),
    width = 12, height = 10, units = "in", res = 300)
par(mfrow = c(2, 2))
plot(simulateResiduals(m_full))
dev.off()

# Extract coefficients
full_model_coefs <- broom.mixed::tidy(m_full, conf.int = TRUE) %>%
  filter(effect == "fixed", term != "(Intercept)")

write_csv(full_model_coefs,
          file.path(SURVEY_TABLES, "comprehensive_model_coefficients.csv"))

# Plot coefficients
p_coefs <- full_model_coefs %>%
  mutate(term = gsub("scale\\(|\\)", "", term)) %>%
  mutate(term = gsub("log\\(|\\+ 1", "", term)) %>%
  ggplot(aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(size = 3) +
  labs(title = "Comprehensive Model Coefficients",
       subtitle = "Negative Binomial GLMM for CAFI Abundance",
       x = "Coefficient Estimate",
       y = "Predictor")

ggsave(file.path(fig_dir, "comprehensive_model_coefficients.png"),
       p_coefs, width = 10, height = 8, dpi = 300)

# 2. MODEL SELECTION
# ------------------
cat("  Performing model selection...\n")

# Dredge for best model combinations
model_set <- dredge(m_full, rank = "AICc",
                   subset = dc(morphotype, branch_category))  # Keep morphology

# Get top models
top_models <- subset(model_set, delta < 4)
write.csv(as.data.frame(top_models),
          file.path(SURVEY_TABLES, "comprehensive_model_selection.csv"))

# Model averaging
if (nrow(top_models) > 1) {
  avg_model <- model.avg(top_models)
  avg_coefs <- summary(avg_model)$coefmat.subset

  write.csv(avg_coefs,
            file.path(SURVEY_TABLES, "model_averaged_coefficients_comprehensive.csv"))

  cat("  Model averaging complete (", nrow(top_models), "models)\n")
}

# ============================================================================
# Variance Partitioning Analysis
# ============================================================================

cat("Performing variance partitioning...\n")

# Prepare community matrix
comm_matrix <- comprehensive_data %>%
  select(coral_id, starts_with("n_")) %>%
  column_to_rownames("coral_id") %>%
  as.matrix()

# Hellinger transform
comm_hell <- decostand(comm_matrix, method = "hellinger")

# Variance partitioning between major predictor groups
varpart_result <- varpart(
  comm_hell,
  # Coral size and morphology
  ~ coral_volume + morphotype + branch_category,
  # Spatial and environmental
  ~ depth + lat + long,
  # Neighbor effects
  ~ mean_neighbor_distance + number_of_neighbors,
  data = comprehensive_data
)

# Plot variance partitioning
png(file.path(fig_dir, "variance_partitioning_comprehensive.png"),
    width = 10, height = 8, units = "in", res = 300)
plot(varpart_result,
     Xnames = c("Coral Characteristics", "Spatial/Environmental", "Neighbor Effects"),
     main = "Variance Partitioning of CAFI Community")
dev.off()

# Save variance partitioning results
varpart_df <- data.frame(
  Component = c("Coral Characteristics", "Spatial/Environmental", "Neighbor Effects",
                "Coral + Spatial", "Coral + Neighbor", "Spatial + Neighbor",
                "All Three", "Residuals"),
  Variance_Explained = c(varpart_result$part$indfract$Adj.R.squared,
                        varpart_result$part$fract$Adj.R.squared[7])
)

write_csv(varpart_df,
          file.path(SURVEY_TABLES, "variance_partitioning_results.csv"))

cat("✓ Variance partitioning complete\n\n")

# ============================================================================
# Distance Decay Analysis
# ============================================================================

cat("Analyzing distance decay patterns...\n")

# Calculate pairwise distances
geo_dist <- dist(comprehensive_data[, c("lat", "long")])
size_dist <- dist(log(comprehensive_data$coral_volume + 1))
neighbor_dist <- dist(comprehensive_data$mean_neighbor_distance)

# Calculate community dissimilarity
cafi_comm <- comprehensive_data %>%
  select(coral_id, starts_with("n_")) %>%
  column_to_rownames("coral_id")

comm_dist <- vegdist(cafi_comm, method = "bray")

# Multiple regression on distance matrices (MRM)
library(ecodist)

# Prepare distance matrices
dist_data <- data.frame(
  comm_dist = as.vector(comm_dist),
  geo_dist = as.vector(geo_dist),
  size_dist = as.vector(size_dist),
  neighbor_dist = as.vector(neighbor_dist)
)

# MRM analysis
mrm_result <- MRM(comm_dist ~ geo_dist + size_dist + neighbor_dist,
                  data = dist_data, nperm = 999)

# Save MRM results
write.csv(mrm_result$coef,
          file.path(SURVEY_TABLES, "distance_decay_mrm_results.csv"))

# Plot distance decay
p_decay <- dist_data %>%
  sample_n(min(5000, nrow(.))) %>%  # Subsample for plotting
  pivot_longer(cols = -comm_dist, names_to = "distance_type", values_to = "distance") %>%
  mutate(distance_type = factor(distance_type,
                                labels = c("Geographic", "Size", "Neighbor"))) %>%
  ggplot(aes(x = distance, y = comm_dist)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam", color = "red") +
  facet_wrap(~distance_type, scales = "free_x") +
  labs(title = "Distance Decay of Community Similarity",
       x = "Distance",
       y = "Community Dissimilarity (Bray-Curtis)")

ggsave(file.path(fig_dir, "distance_decay_patterns.png"),
       p_decay, width = 14, height = 6, dpi = 300)

cat("✓ Distance decay analysis complete\n\n")

# ============================================================================
# Hierarchical Analysis of Predictors
# ============================================================================

cat("Performing hierarchical analysis...\n")

# Random forest for variable importance
library(randomForest)

# Prepare data for RF
rf_data <- comprehensive_data %>%
  select(total_cafi, coral_volume, depth, morphotype, branch_category,
         number_of_neighbors, mean_neighbor_distance, lat, long) %>%
  na.omit() %>%
  mutate(across(where(is.character), as.factor))

# Fit random forest
rf_model <- randomForest(total_cafi ~ ., data = rf_data,
                         ntree = 500, importance = TRUE)

# Extract importance
importance_df <- importance(rf_model) %>%
  as.data.frame() %>%
  rownames_to_column("predictor") %>%
  arrange(desc(`%IncMSE`))

write_csv(importance_df,
          file.path(SURVEY_TABLES, "predictor_importance_rf.csv"))

# Plot importance
p_importance <- importance_df %>%
  ggplot(aes(x = reorder(predictor, `%IncMSE`), y = `%IncMSE`)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Predictor Importance (Random Forest)",
       x = "Predictor",
       y = "% Increase in MSE")

ggsave(file.path(fig_dir, "predictor_importance.png"),
       p_importance, width = 10, height = 8, dpi = 300)

cat("✓ Hierarchical analysis complete\n\n")

# ============================================================================
# Size-Scaling Relationships
# ============================================================================

cat("Analyzing size-scaling relationships...\n")

# Power law scaling
scaling_models <- comprehensive_data %>%
  filter(coral_volume > 0, total_cafi > 0) %>%
  group_by(morphotype) %>%
  do(
    abundance_model = lm(log(total_cafi) ~ log(coral_volume), data = .),
    richness_model = lm(log(species_richness + 1) ~ log(coral_volume), data = .)
  )

# Extract scaling exponents
scaling_exponents <- scaling_models %>%
  summarise(
    abundance_slope = coef(abundance_model)[2],
    abundance_se = summary(abundance_model)$coefficients[2, 2],
    richness_slope = coef(richness_model)[2],
    richness_se = summary(richness_model)$coefficients[2, 2]
  )

write_csv(scaling_exponents,
          file.path(SURVEY_TABLES, "size_scaling_exponents.csv"))

# Plot scaling relationships
p_scaling <- comprehensive_data %>%
  filter(coral_volume > 0, total_cafi > 0) %>%
  ggplot(aes(x = coral_volume, y = total_cafi)) +
  geom_point(aes(color = morphotype), alpha = 0.5) +
  geom_smooth(aes(color = morphotype), method = "lm", se = TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Size-Abundance Scaling Relationships",
       subtitle = paste("Scaling exponents by morphotype"),
       x = "Coral Volume (cm³, log scale)",
       y = "Total CAFI (log scale)",
       color = "Morphotype") +
  annotation_logticks()

ggsave(file.path(fig_dir, "size_scaling_relationships.png"),
       p_scaling, width = 10, height = 8, dpi = 300)

cat("✓ Size-scaling analysis complete\n\n")

# ============================================================================
# Interactive Effects Analysis
# ============================================================================

cat("Analyzing interactive effects...\n")

# Three-way interaction: Volume × Morphotype × Depth
interaction_model <- glmmTMB(
  total_cafi ~
    scale(log(coral_volume + 1)) * morphotype * poly(depth, 2) +
    (1|site),
  data = comprehensive_data,
  family = nbinom2()
)

# Predict for visualization
pred_grid <- expand.grid(
  coral_volume = exp(seq(log(min(comprehensive_data$coral_volume, na.rm = TRUE)),
                         log(max(comprehensive_data$coral_volume, na.rm = TRUE)),
                         length.out = 50)),
  morphotype = unique(comprehensive_data$morphotype),
  depth = seq(min(comprehensive_data$depth, na.rm = TRUE),
             max(comprehensive_data$depth, na.rm = TRUE),
             length.out = 50),
  site = levels(factor(comprehensive_data$site))[1]
)

pred_grid$predicted <- predict(interaction_model, pred_grid,
                               type = "response", re.form = NA)

# Plot interaction effects
p_interaction <- pred_grid %>%
  filter(!is.na(morphotype)) %>%
  ggplot(aes(x = coral_volume, y = depth, fill = predicted)) +
  geom_tile() +
  scale_x_log10() +
  scale_fill_viridis_c(trans = "sqrt", name = "Predicted\nCAFI") +
  facet_wrap(~morphotype) +
  labs(title = "Interactive Effects: Volume × Depth × Morphotype",
       x = "Coral Volume (cm³, log scale)",
       y = "Depth (m)")

ggsave(file.path(fig_dir, "interaction_effects_3way.png"),
       p_interaction, width = 14, height = 6, dpi = 300)

cat("✓ Interactive effects analysis complete\n\n")

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Comprehensive Predictor Analysis Summary\n")
cat("========================================\n\n")

cat("Key Predictors Analyzed:\n")
cat("  - Coral size metrics (volume, surface area, complexity)\n")
cat("  - Neighbor effects (distance, density, isolation)\n")
cat("  - Spatial factors (depth, latitude, longitude)\n")
cat("  - Morphological traits (morphotype, branch width)\n")
cat("  - Interactive effects (3-way interactions)\n\n")

cat("Statistical Analyses:\n")
cat("  - Comprehensive GLMMs with all predictors\n")
cat("  - Model selection and averaging\n")
cat("  - Variance partitioning\n")
cat("  - Distance decay (MRM)\n")
cat("  - Random forest importance\n")
cat("  - Size-scaling relationships\n\n")

if (exists("avg_model")) {
  cat("Top Predictors (by model averaging):\n")
  top_preds <- as.data.frame(avg_coefs) %>%
    rownames_to_column("predictor") %>%
    arrange(desc(abs(Estimate))) %>%
    head(5)

  for (i in 1:nrow(top_preds)) {
    cat("  ", i, ". ", top_preds$predictor[i],
        " (β = ", round(top_preds$Estimate[i], 3), ")\n", sep = "")
  }
  cat("\n")
}

if (exists("varpart_df")) {
  cat("Variance Explained:\n")
  cat("  - Coral characteristics: ",
      round(varpart_df$Variance_Explained[1] * 100, 1), "%\n", sep = "")
  cat("  - Spatial/Environmental: ",
      round(varpart_df$Variance_Explained[2] * 100, 1), "%\n", sep = "")
  cat("  - Neighbor effects: ",
      round(varpart_df$Variance_Explained[3] * 100, 1), "%\n", sep = "")
  cat("\n")
}

if (exists("scaling_exponents")) {
  cat("Size-Scaling Exponents:\n")
  for (i in 1:nrow(scaling_exponents)) {
    cat("  - ", scaling_exponents$morphotype[i], ": ",
        round(scaling_exponents$abundance_slope[i], 2),
        " ± ", round(scaling_exponents$abundance_se[i], 2), "\n", sep = "")
  }
}

cat("\n✅ Comprehensive predictor analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", SURVEY_TABLES, "\n")