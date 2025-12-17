#!/usr/bin/env Rscript
# ============================================================================
# 19_condition_vs_landscape.R - Coral Condition vs Landscape Characteristics
#
# Purpose: Test whether coral physiological condition varies with landscape
#          characteristics (size, neighborhood, isolation)
#
# Questions:
#   1. Does coral condition vary with coral size?
#   2. Does condition vary with neighborhood density/composition?
#   3. Does isolation affect coral condition?
#   4. Do these patterns hold for individual physiological metrics?
#
# Response variables:
#   - Condition PC1 (composite score)
#   - Protein (mg/cm²)
#   - Carbohydrate (mg/cm²)
#   - Zooxanthellae density (cells/cm²)
#   - AFDW (mg/cm²)
#
# Predictor variables:
#   - Coral volume (focal coral size)
#   - Neighbor count (within 5m)
#   - Neighbor volume (total habitat nearby)
#   - Isolation index (distance / size^1/3)
#   - Relative size (focal / neighbor volume)
#   - Crowding index (neighbor_volume / mean_distance)
#   - Site (environmental context)
#
# Author: CAFI Analysis Pipeline
# Date: 2025-12-02
# ============================================================================

cat("============================================================\n")
cat("Script 19: Coral Condition vs Landscape Characteristics\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(car)
  library(broom)
  library(patchwork)
  library(corrplot)
})

# Set paths
setwd("/Users/adrianstiermbp2023/CAFI-Survey-2026")

# Create output directories
fig_dir <- "output/figures/condition_vs_landscape"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")

# Load condition scores
condition_file <- "output/tables/coral_condition_scores.csv"
if (file.exists(condition_file)) {
  condition_data <- read.csv(condition_file)
} else {
  condition_data <- readRDS("output/objects/coral_condition_scores.rds")
}

# Load coral characteristics
coral_data <- read.csv("data/survey_coral_characteristics_merged_v2.csv")

# Load neighborhood data if available
neighbor_file <- "output/tables/neighborhood_metrics.csv"
if (file.exists(neighbor_file)) {
  neighbor_data <- read.csv(neighbor_file)
  has_neighbors <- TRUE
} else {
  # Try to calculate from coral data
  has_neighbors <- FALSE
  cat("  Note: Pre-calculated neighborhood metrics not found\n")
  cat("        Will calculate from GPS coordinates if available\n")
}

cat("  - Condition data:", nrow(condition_data), "corals\n")
cat("  - Coral characteristics:", nrow(coral_data), "corals\n\n")

# ============================================================================
# DATA PREPARATION
# ============================================================================

cat("Preparing analysis data...\n")

# Process coral data - neighborhood variables already exist in raw data
coral_processed <- coral_data %>%
  mutate(
    site = str_extract(site, "^[A-Z]+"),
    volume = coalesce(volume_field, volume_lab),
    log_volume = log10(volume + 1),
    depth_m = depth,
    # Rename neighborhood variables from raw data
    n_neighbors = number_of_neighbors,
    neighbor_volume = combined_total_volume_of_neighbors,
    mean_neighbor_distance_raw = mean_neighbor_distance,
    mean_neighbor_volume = mean_total_volume_of_neighbors
  ) %>%
  filter(!is.na(volume), volume > 0, site %in% c("HAU", "MAT", "MRB")) %>%
  # Calculate derived neighborhood metrics
  mutate(
    isolation_index = ifelse(n_neighbors > 0 & !is.na(mean_neighbor_distance_raw),
                             mean_neighbor_distance_raw / (volume^(1/3)), NA),
    relative_size = ifelse(!is.na(mean_neighbor_volume) & mean_neighbor_volume > 0,
                           volume / mean_neighbor_volume, NA),
    crowding_index = ifelse(!is.na(mean_neighbor_distance_raw) & mean_neighbor_distance_raw > 0,
                            neighbor_volume / mean_neighbor_distance_raw, NA)
  )

# Check if we have neighborhood data
has_neighbors <- sum(!is.na(coral_processed$n_neighbors)) > 20
cat("  - Neighborhood data available:", has_neighbors, "\n")
if (has_neighbors) {
  cat("    Corals with neighbor data:", sum(!is.na(coral_processed$n_neighbors)), "\n")
}

# Merge condition with coral data
# Note: condition_data already has site, so we use suffix to avoid conflict
analysis_data <- condition_data %>%
  left_join(coral_processed %>%
              select(coral_id, volume, log_volume,
                     depth_m, lat, long, morphotype, branch_width,
                     n_neighbors, neighbor_volume, mean_neighbor_distance_raw,
                     mean_neighbor_volume, isolation_index, relative_size, crowding_index),
            by = "coral_id") %>%
  filter(!is.na(volume), !is.na(condition_score))

cat("  - Analysis dataset:", nrow(analysis_data), "corals with condition data\n")
if (has_neighbors) {
  n_with_neighbors <- sum(!is.na(analysis_data$n_neighbors))
  cat("  - With neighborhood data:", n_with_neighbors, "corals\n")
}
cat("\n")

# ============================================================================
# PART 1: CONDITION PC1 VS LANDSCAPE
# ============================================================================

cat("============================================================\n")
cat("PART 1: Condition Score (PC1) vs Landscape Characteristics\n")
cat("============================================================\n\n")

# 1a. Condition vs Coral Size
cat("1a. Condition vs Coral Size:\n")

m_size <- lm(condition_score ~ log_volume, data = analysis_data)
size_summary <- summary(m_size)

cat(sprintf("    Condition ~ log(Volume):\n"))
cat(sprintf("      β = %.4f (SE = %.4f)\n", coef(m_size)[2], size_summary$coefficients[2, 2]))
cat(sprintf("      R² = %.4f, p = %.4f\n", size_summary$r.squared, size_summary$coefficients[2, 4]))

# Plot
p_size <- ggplot(analysis_data, aes(x = volume, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_log10() +
  scale_color_viridis_d() +
  labs(
    title = "Coral Condition vs Size",
    subtitle = sprintf("β = %.3f, R² = %.3f, p = %.4f",
                       coef(m_size)[2], size_summary$r.squared,
                       size_summary$coefficients[2, 4]),
    x = "Coral Volume (cm³, log scale)",
    y = "Condition Score (PC1)",
    color = "Site"
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(fig_dir, "condition_vs_size.png"), p_size,
       width = 10, height = 7, dpi = 300, bg = "white")

# 1b. Condition vs Site
cat("\n1b. Condition vs Site:\n")

m_site <- lm(condition_score ~ site, data = analysis_data)
anova_site <- anova(m_site)
cat(sprintf("    ANOVA: F = %.2f, p = %.4f\n",
            anova_site$`F value`[1], anova_site$`Pr(>F)`[1]))

# Site means
site_means <- analysis_data %>%
  group_by(site) %>%
  summarise(
    n = n(),
    mean_condition = mean(condition_score, na.rm = TRUE),
    se = sd(condition_score, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
print(site_means)

# Plot
p_site <- ggplot(analysis_data, aes(x = site, y = condition_score, fill = site)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_fill_viridis_d() +
  labs(
    title = "Coral Condition by Site",
    subtitle = sprintf("ANOVA: F = %.2f, p = %.4f",
                       anova_site$`F value`[1], anova_site$`Pr(>F)`[1]),
    x = "Site",
    y = "Condition Score (PC1)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "condition_by_site.png"), p_site,
       width = 8, height = 6, dpi = 300, bg = "white")

# 1c. Condition vs Neighborhood (if available)
if (has_neighbors && sum(!is.na(analysis_data$n_neighbors)) > 20) {
  cat("\n1c. Condition vs Neighborhood Characteristics:\n")

  neighbor_vars <- c("n_neighbors", "neighbor_volume", "isolation_index",
                     "relative_size", "crowding_index")

  neighbor_results <- data.frame()

  for (var in neighbor_vars) {
    if (var %in% names(analysis_data) && sum(!is.na(analysis_data[[var]])) > 15) {
      # Simple correlation
      valid <- !is.na(analysis_data[[var]]) & !is.na(analysis_data$condition_score)
      if (sum(valid) > 15) {
        cor_test <- cor.test(analysis_data[[var]][valid],
                             analysis_data$condition_score[valid])

        # Linear model
        formula <- as.formula(paste("condition_score ~", var))
        m <- lm(formula, data = analysis_data)
        m_sum <- summary(m)

        neighbor_results <- rbind(neighbor_results, data.frame(
          variable = var,
          n = sum(valid),
          correlation = cor_test$estimate,
          cor_p = cor_test$p.value,
          beta = coef(m)[2],
          se = m_sum$coefficients[2, 2],
          r_squared = m_sum$r.squared,
          lm_p = m_sum$coefficients[2, 4]
        ))

        cat(sprintf("    %s: r = %.3f (p = %.4f), β = %.4f, R² = %.4f\n",
                    var, cor_test$estimate, cor_test$p.value,
                    coef(m)[2], m_sum$r.squared))
      }
    }
  }

  if (nrow(neighbor_results) > 0) {
    write.csv(neighbor_results,
              file.path("output/tables", "condition_vs_neighborhood_results.csv"),
              row.names = FALSE)

    # Multi-panel plot
    neighbor_plots <- list()

    for (var in neighbor_vars[neighbor_vars %in% names(analysis_data)]) {
      if (sum(!is.na(analysis_data[[var]])) > 15) {
        p <- ggplot(analysis_data, aes_string(x = var, y = "condition_score")) +
          geom_point(aes(color = site), alpha = 0.6, size = 2) +
          geom_smooth(method = "lm", color = "black", se = TRUE) +
          geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
          scale_color_viridis_d() +
          labs(title = var, x = var, y = "Condition PC1") +
          theme_bw(base_size = 10) +
          theme(legend.position = "none")

        neighbor_plots[[var]] <- p
      }
    }

    if (length(neighbor_plots) > 0) {
      p_neighbors <- wrap_plots(neighbor_plots, ncol = 3) +
        plot_annotation(
          title = "Coral Condition vs Neighborhood Characteristics",
          subtitle = "Each panel shows condition PC1 vs a different neighborhood metric"
        )

      ggsave(file.path(fig_dir, "condition_vs_neighborhood_panel.png"),
             p_neighbors, width = 14, height = 10, dpi = 300, bg = "white")
    }
  }
}

# 1d. Combined model (size + site + neighbors)
cat("\n1d. Combined Model:\n")

if (has_neighbors && sum(!is.na(analysis_data$n_neighbors)) > 20) {
  m_combined <- lm(condition_score ~ log_volume + site + n_neighbors,
                   data = analysis_data)
} else {
  m_combined <- lm(condition_score ~ log_volume + site, data = analysis_data)
}

cat("    Model:", deparse(formula(m_combined)), "\n")
print(summary(m_combined)$coefficients)
cat(sprintf("    Overall R² = %.4f\n", summary(m_combined)$r.squared))

# ============================================================================
# PART 2: INDIVIDUAL PHYSIOLOGICAL METRICS VS LANDSCAPE
# ============================================================================

cat("\n============================================================\n")
cat("PART 2: Individual Physiological Metrics vs Landscape\n")
cat("============================================================\n\n")

# Get physiological variables (position-corrected if available)
physio_vars <- c("protein_corrected", "carb_corrected", "zoox_corrected", "afdw_corrected")
physio_names <- c("Protein", "Carbohydrate", "Zooxanthellae", "AFDW")

# Check which exist
available_physio <- physio_vars[physio_vars %in% names(analysis_data)]

if (length(available_physio) > 0) {

  # 2a. Each physio metric vs size
  cat("2a. Individual Metrics vs Coral Size:\n\n")

  physio_size_results <- data.frame()
  physio_plots <- list()

  for (i in seq_along(available_physio)) {
    var <- available_physio[i]
    name <- physio_names[which(physio_vars == var)]

    valid <- !is.na(analysis_data[[var]]) & !is.na(analysis_data$log_volume)

    if (sum(valid) > 15) {
      formula <- as.formula(paste(var, "~ log_volume"))
      m <- lm(formula, data = analysis_data)
      m_sum <- summary(m)

      physio_size_results <- rbind(physio_size_results, data.frame(
        metric = name,
        variable = var,
        n = sum(valid),
        beta = coef(m)[2],
        se = m_sum$coefficients[2, 2],
        r_squared = m_sum$r.squared,
        p_value = m_sum$coefficients[2, 4]
      ))

      cat(sprintf("    %s ~ log(Volume): β = %.4f, R² = %.4f, p = %.4f\n",
                  name, coef(m)[2], m_sum$r.squared, m_sum$coefficients[2, 4]))

      # Plot
      p <- ggplot(analysis_data, aes_string(x = "volume", y = var)) +
        geom_point(aes(color = site), alpha = 0.6, size = 2) +
        geom_smooth(method = "lm", color = "black", se = TRUE) +
        scale_x_log10() +
        scale_color_viridis_d() +
        labs(title = name, x = "Coral Volume (cm³)", y = paste(name, "(corrected)")) +
        theme_bw(base_size = 10) +
        theme(legend.position = "none")

      physio_plots[[name]] <- p
    }
  }

  # Combined physio vs size plot
  if (length(physio_plots) > 0) {
    p_physio_size <- wrap_plots(physio_plots, ncol = 2) +
      plot_annotation(
        title = "Individual Physiological Metrics vs Coral Size",
        subtitle = "Position-corrected metrics (residuals from stump length regression)"
      )

    ggsave(file.path(fig_dir, "individual_physio_vs_size.png"),
           p_physio_size, width = 12, height = 10, dpi = 300, bg = "white")
  }

  # 2b. Each physio metric vs site
  cat("\n2b. Individual Metrics by Site:\n\n")

  physio_site_results <- data.frame()

  for (i in seq_along(available_physio)) {
    var <- available_physio[i]
    name <- physio_names[which(physio_vars == var)]

    valid <- !is.na(analysis_data[[var]]) & !is.na(analysis_data$site)

    if (sum(valid) > 15) {
      formula <- as.formula(paste(var, "~ site"))
      m <- lm(formula, data = analysis_data)
      anova_result <- anova(m)

      physio_site_results <- rbind(physio_site_results, data.frame(
        metric = name,
        variable = var,
        n = sum(valid),
        f_value = anova_result$`F value`[1],
        p_value = anova_result$`Pr(>F)`[1]
      ))

      cat(sprintf("    %s ~ Site: F = %.2f, p = %.4f\n",
                  name, anova_result$`F value`[1], anova_result$`Pr(>F)`[1]))
    }
  }

  # 2c. Each physio metric vs neighborhood
  if (has_neighbors && sum(!is.na(analysis_data$n_neighbors)) > 15) {
    cat("\n2c. Individual Metrics vs Neighborhood:\n\n")

    physio_neighbor_results <- data.frame()

    for (var in available_physio) {
      name <- physio_names[which(physio_vars == var)]

      for (nvar in c("n_neighbors", "neighbor_volume", "isolation_index")) {
        if (nvar %in% names(analysis_data)) {
          valid <- !is.na(analysis_data[[var]]) & !is.na(analysis_data[[nvar]])

          if (sum(valid) > 15) {
            cor_test <- cor.test(analysis_data[[var]][valid],
                                 analysis_data[[nvar]][valid])

            physio_neighbor_results <- rbind(physio_neighbor_results, data.frame(
              physio_metric = name,
              neighbor_metric = nvar,
              n = sum(valid),
              correlation = cor_test$estimate,
              p_value = cor_test$p.value
            ))
          }
        }
      }
    }

    if (nrow(physio_neighbor_results) > 0) {
      # Print summary
      wide_results <- physio_neighbor_results %>%
        select(physio_metric, neighbor_metric, correlation) %>%
        pivot_wider(names_from = neighbor_metric, values_from = correlation)

      cat("    Correlation matrix (physio vs neighbor metrics):\n")
      print(wide_results)

      write.csv(physio_neighbor_results,
                file.path("output/tables", "physio_vs_neighborhood_correlations.csv"),
                row.names = FALSE)
    }
  }

  # Save results
  if (nrow(physio_size_results) > 0) {
    write.csv(physio_size_results,
              file.path("output/tables", "physio_vs_size_results.csv"),
              row.names = FALSE)
  }

  if (nrow(physio_site_results) > 0) {
    write.csv(physio_site_results,
              file.path("output/tables", "physio_vs_site_results.csv"),
              row.names = FALSE)
  }

} else {
  cat("  Note: Individual physiological metrics not found in condition data\n")
  cat("        Only condition PC1 available\n")
}

# ============================================================================
# PART 3: COMPREHENSIVE CORRELATION MATRIX
# ============================================================================

cat("\n============================================================\n")
cat("PART 3: Comprehensive Correlation Matrix\n")
cat("============================================================\n\n")

# Select numeric variables for correlation
corr_vars <- c("condition_score", "log_volume", "depth_m")
if (length(available_physio) > 0) {
  corr_vars <- c(corr_vars, available_physio)
}
if (has_neighbors) {
  neighbor_vars_available <- intersect(
    c("n_neighbors", "neighbor_volume", "isolation_index", "relative_size", "crowding_index"),
    names(analysis_data)
  )
  corr_vars <- c(corr_vars, neighbor_vars_available)
}

# Get only variables that exist
corr_vars <- corr_vars[corr_vars %in% names(analysis_data)]

# Create correlation matrix
corr_data <- analysis_data[, corr_vars, drop = FALSE]
corr_data <- corr_data[complete.cases(corr_data), ]

if (nrow(corr_data) > 10 && ncol(corr_data) > 2) {
  corr_matrix <- cor(corr_data, use = "pairwise.complete.obs")

  # Save correlation matrix
  write.csv(as.data.frame(corr_matrix),
            file.path("output/tables", "condition_landscape_correlation_matrix.csv"))

  # Plot correlation matrix
  png(file.path(fig_dir, "condition_landscape_correlation_matrix.png"),
      width = 12, height = 10, units = "in", res = 300, bg = "white")

  corrplot(corr_matrix,
           method = "color",
           type = "upper",
           order = "hclust",
           addCoef.col = "black",
           number.cex = 0.7,
           tl.col = "black",
           tl.srt = 45,
           title = "Correlation Matrix: Condition & Landscape Variables",
           mar = c(0, 0, 2, 0))

  dev.off()

  cat("  Correlation matrix saved with", ncol(corr_matrix), "variables\n")
}

# ============================================================================
# PART 4: SUMMARY TABLE FOR MANUSCRIPT
# ============================================================================

cat("\n============================================================\n")
cat("PART 4: Summary for Manuscript\n")
cat("============================================================\n\n")

summary_results <- data.frame(
  Response = character(),
  Predictor = character(),
  N = integer(),
  Effect = numeric(),
  SE = numeric(),
  Statistic = character(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# Condition PC1 vs size
summary_results <- rbind(summary_results, data.frame(
  Response = "Condition PC1",
  Predictor = "log(Volume)",
  N = nrow(analysis_data),
  Effect = coef(m_size)[2],
  SE = summary(m_size)$coefficients[2, 2],
  Statistic = sprintf("R² = %.3f", summary(m_size)$r.squared),
  P_value = summary(m_size)$coefficients[2, 4]
))

# Condition PC1 vs site
summary_results <- rbind(summary_results, data.frame(
  Response = "Condition PC1",
  Predictor = "Site",
  N = nrow(analysis_data),
  Effect = NA,
  SE = NA,
  Statistic = sprintf("F = %.2f", anova_site$`F value`[1]),
  P_value = anova_site$`Pr(>F)`[1]
))

# Add neighbor results if available
if (exists("neighbor_results") && nrow(neighbor_results) > 0) {
  for (i in 1:nrow(neighbor_results)) {
    summary_results <- rbind(summary_results, data.frame(
      Response = "Condition PC1",
      Predictor = neighbor_results$variable[i],
      N = neighbor_results$n[i],
      Effect = neighbor_results$beta[i],
      SE = neighbor_results$se[i],
      Statistic = sprintf("R² = %.3f", neighbor_results$r_squared[i]),
      P_value = neighbor_results$lm_p[i]
    ))
  }
}

# Add individual physio results
if (exists("physio_size_results") && nrow(physio_size_results) > 0) {
  for (i in 1:nrow(physio_size_results)) {
    summary_results <- rbind(summary_results, data.frame(
      Response = physio_size_results$metric[i],
      Predictor = "log(Volume)",
      N = physio_size_results$n[i],
      Effect = physio_size_results$beta[i],
      SE = physio_size_results$se[i],
      Statistic = sprintf("R² = %.3f", physio_size_results$r_squared[i]),
      P_value = physio_size_results$p_value[i]
    ))
  }
}

# Save summary
write.csv(summary_results,
          file.path("output/tables", "condition_vs_landscape_summary.csv"),
          row.names = FALSE)

# Print summary
cat("Summary of Condition vs Landscape Results:\n")
cat("─────────────────────────────────────────────────────────────────────\n")
print(summary_results %>%
        mutate(P_value = format.pval(P_value, digits = 3)))
cat("─────────────────────────────────────────────────────────────────────\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("Key Findings:\n")
cat(sprintf("  1. Condition vs Size: β = %.3f, p = %.4f\n",
            coef(m_size)[2], summary(m_size)$coefficients[2, 4]))
cat(sprintf("  2. Condition vs Site: F = %.2f, p = %.4f\n",
            anova_site$`F value`[1], anova_site$`Pr(>F)`[1]))

if (has_neighbors && exists("neighbor_results") && nrow(neighbor_results) > 0) {
  best_neighbor <- neighbor_results[which.min(neighbor_results$lm_p), ]
  cat(sprintf("  3. Best neighborhood predictor: %s (R² = %.3f, p = %.4f)\n",
              best_neighbor$variable, best_neighbor$r_squared, best_neighbor$lm_p))
}

cat("\nOutput files:\n")
cat("  - Figures:", fig_dir, "\n")
cat("  - Tables: output/tables/condition_vs_*.csv\n\n")

cat("✅ Script 19 complete!\n")
