#!/usr/bin/env Rscript
# ============================================================================
# 19_figure2_abundance_richness_landscape.R - Landscape Effects on CAFI
# ============================================================================
#
# PURPOSE: Test how coral size and neighborhood context shape CAFI abundance,
# species richness, and Shannon diversity.
#
# HYPOTHESES:
#   H1: Larger corals host more CAFI (positive size-abundance slope)
#   H2: Isolated corals host more CAFI per unit size (propagule redirection)
#   H3: Diversity increases with coral size but may vary with neighborhood
#
# THEORETICAL BACKGROUND:
#   Propagule redirection theory predicts that occupant abundance scales
#   nonlinearly with habitat amount (power-law with exponent < 1). This causes
#   occupant density to decrease as habitat increases - larger corals intercept
#   more propagules but at lower density per unit volume.
#
# OUTPUTS:
#   - 6-panel figure: Abundance, Richness, Diversity vs Size and Proximity
#   - Statistical model results table
#   - Saved to neighborhood figures and manuscript folders
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("19: Landscape Effects on CAFI Communities\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Load MASS for negative binomial
library(MASS)

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "neighborhood")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# PREPARE ANALYSIS DATA
# ============================================================================

cat("Preparing analysis data...\n")

# Filter to corals with complete data for landscape analyses
landscape_data <- coral_master %>%
  filter(
    !is.na(volume), volume > 0,
    !is.na(total_cafi),
    !is.na(mean_neighbor_dist)
  ) %>%
  mutate(
    # Log transformations
    log_volume = log10(volume),
    log_cafi = log(total_cafi + 1),

    # Proximity metric: distance to nearest neighbors (cm -> m)
    proximity_m = mean_neighbor_dist / 100,
    log_proximity = log10(proximity_m + 0.1),

    # Standardize predictors for model stability
    volume_z = scale(log_volume)[, 1],
    proximity_z = scale(proximity_m)[, 1],

    # Create isolation categories for visualization
    isolation_cat = cut(proximity_m,
                        breaks = quantile(proximity_m, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                        labels = c("Clustered", "Intermediate", "Isolated"),
                        include.lowest = TRUE)
  )

cat(sprintf("  - Analysis sample: %d corals\n", nrow(landscape_data)))
cat(sprintf("  - Volume range: %.0f - %.0f cm³\n",
            min(landscape_data$volume), max(landscape_data$volume)))
cat(sprintf("  - Proximity range: %.2f - %.2f m\n",
            min(landscape_data$proximity_m), max(landscape_data$proximity_m)))

# Helper function to extract model statistics
extract_model_stats <- function(model, term) {
  coefs <- summary(model)$coefficients
  if (!(term %in% rownames(coefs))) {
    return(list(estimate = NA, se = NA, p.value = NA, conf.low = NA, conf.high = NA, significant = FALSE))
  }

  est <- coefs[term, "Estimate"]
  se <- coefs[term, "Std. Error"]

  # Get p-value (different column names for different model types)
  if ("Pr(>|z|)" %in% colnames(coefs)) {
    p <- coefs[term, "Pr(>|z|)"]
  } else if ("Pr(>|t|)" %in% colnames(coefs)) {
    p <- coefs[term, "Pr(>|t|)"]
  } else {
    p <- NA
  }

  # 95% CI
  ci <- tryCatch({
    confint(model, parm = term, level = 0.95)
  }, error = function(e) {
    c(est - 1.96 * se, est + 1.96 * se)
  })

  if (is.matrix(ci)) {
    conf.low <- ci[1]
    conf.high <- ci[2]
  } else {
    conf.low <- ci[1]
    conf.high <- ci[2]
  }

  list(
    estimate = est,
    se = se,
    p.value = p,
    conf.low = conf.low,
    conf.high = conf.high,
    significant = !is.na(p) && p < 0.05
  )
}

format_p <- function(p) {
  if (p < 0.001) return("< 0.001")
  if (p < 0.01) return(sprintf("= %.3f", p))
  sprintf("= %.2f", p)
}

# ============================================================================
# ANALYSIS 1: CAFI ABUNDANCE vs CORAL SIZE (Power-law Scaling)
# ============================================================================

cat("\n1. Testing CAFI abundance ~ coral size (power-law scaling)...\n")

# Negative binomial GLM (appropriate for overdispersed counts)
model_abundance_size <- glm.nb(
  total_cafi ~ log_volume + site,
  data = landscape_data
)

coef_size_abund <- extract_model_stats(model_abundance_size, "log_volume")

cat(sprintf("  Scaling exponent: %.3f (SE = %.3f)\n",
            coef_size_abund$estimate, coef_size_abund$se))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n",
            coef_size_abund$conf.low, coef_size_abund$conf.high))
cat(sprintf("  P-value %s\n", format_p(coef_size_abund$p.value)))

if (coef_size_abund$conf.high < 1) {
  cat("  -> SUBLINEAR scaling confirmed (exponent < 1)\n")
} else if (coef_size_abund$conf.low > 1) {
  cat("  -> SUPERLINEAR scaling (exponent > 1)\n")
} else {
  cat("  -> Linear scaling (CI includes 1)\n")
}

# Panel A: Size-Abundance Scaling
show_line_2a <- coef_size_abund$p.value < 0.10

p_2a <- ggplot(landscape_data, aes(x = volume, y = total_cafi)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  {if (show_line_2a)
    geom_smooth(method = "glm.nb", formula = y ~ log10(x),
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_x_log10(labels = scales::comma) +
  scale_y_sqrt(breaks = c(0, 25, 100, 225, 400)) +
  scale_color_site() +
  labs(
    title = "A. CAFI abundance scales with coral size",
    subtitle = sprintf("Exponent = %.2f [%.2f, %.2f], p %s",
                       coef_size_abund$estimate,
                       coef_size_abund$conf.low,
                       coef_size_abund$conf.high,
                       format_p(coef_size_abund$p.value)),
    x = expression(Coral~Volume~(cm^3~","~log~scale)),
    y = "CAFI Abundance (sqrt scale)",
    color = "Site"
  ) +
  theme_publication() +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# ============================================================================
# ANALYSIS 2: CAFI ABUNDANCE vs PROXIMITY (Propagule Redirection)
# ============================================================================

cat("\n2. Testing CAFI abundance ~ proximity (propagule redirection)...\n")

model_abundance_prox <- glm.nb(
  total_cafi ~ proximity_m + log_volume + site,
  data = landscape_data
)

coef_prox_abund <- extract_model_stats(model_abundance_prox, "proximity_m")

cat(sprintf("  Proximity effect: %.3f (SE = %.3f)\n",
            coef_prox_abund$estimate, coef_prox_abund$se))
cat(sprintf("  P-value %s\n", format_p(coef_prox_abund$p.value)))

if (coef_prox_abund$significant && coef_prox_abund$estimate > 0) {
  cat("  -> Isolated corals host MORE CAFI (propagule redirection supported)\n")
} else if (coef_prox_abund$significant && coef_prox_abund$estimate < 0) {
  cat("  -> Clustered corals host MORE CAFI (spillover/facilitation)\n")
} else {
  cat("  -> No significant proximity effect\n")
}

# Panel B: Proximity effect
show_line_2b <- coef_prox_abund$p.value < 0.10

p_2b <- ggplot(landscape_data, aes(x = proximity_m, y = total_cafi)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  {if (show_line_2b)
    geom_smooth(method = "glm.nb", formula = y ~ x,
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_y_sqrt(breaks = c(0, 25, 100, 225, 400)) +
  scale_color_site() +
  labs(
    title = "B. CAFI abundance vs neighbor proximity",
    subtitle = sprintf("Proximity effect = %.3f, p %s | Controlling for coral size",
                       coef_prox_abund$estimate,
                       format_p(coef_prox_abund$p.value)),
    x = "Mean Distance to Neighbors (m)",
    y = "CAFI Abundance (sqrt scale)",
    color = "Site"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# ============================================================================
# ANALYSIS 3: SPECIES RICHNESS vs CORAL SIZE
# ============================================================================

cat("\n3. Testing OTU richness ~ coral size...\n")

# Start with Poisson, check for overdispersion
model_richness_size <- glm(
  otu_richness ~ log_volume + site,
  data = landscape_data,
  family = poisson()
)

overdispersion <- sum(residuals(model_richness_size, type = "pearson")^2) /
  model_richness_size$df.residual

if (overdispersion > 1.5) {
  cat("  Overdispersion detected, using negative binomial...\n")
  model_richness_size <- glm.nb(
    otu_richness ~ log_volume + site,
    data = landscape_data
  )
}

coef_size_rich <- extract_model_stats(model_richness_size, "log_volume")

cat(sprintf("  Size effect on richness: %.3f (SE = %.3f)\n",
            coef_size_rich$estimate, coef_size_rich$se))
cat(sprintf("  P-value %s\n", format_p(coef_size_rich$p.value)))

# Panel C: Richness vs Size
show_line_2c <- coef_size_rich$p.value < 0.10

p_2c <- ggplot(landscape_data, aes(x = volume, y = otu_richness)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  {if (show_line_2c)
    geom_smooth(method = "glm", formula = y ~ log10(x),
                method.args = list(family = poisson()),
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_x_log10(labels = scales::comma) +
  scale_color_site() +
  labs(
    title = "C. OTU richness scales with coral size",
    subtitle = sprintf("Size effect = %.2f, p %s",
                       coef_size_rich$estimate,
                       format_p(coef_size_rich$p.value)),
    x = expression(Coral~Volume~(cm^3~","~log~scale)),
    y = "OTU Richness",
    color = "Site"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# ============================================================================
# ANALYSIS 4: SPECIES RICHNESS vs PROXIMITY
# ============================================================================

cat("\n4. Testing OTU richness ~ proximity...\n")

model_richness_prox <- glm.nb(
  otu_richness ~ proximity_m + log_volume + site,
  data = landscape_data
)

coef_prox_rich <- extract_model_stats(model_richness_prox, "proximity_m")

cat(sprintf("  Proximity effect on richness: %.3f (SE = %.3f)\n",
            coef_prox_rich$estimate, coef_prox_rich$se))
cat(sprintf("  P-value %s\n", format_p(coef_prox_rich$p.value)))

# Panel D: Richness vs Proximity
show_line_2d <- coef_prox_rich$p.value < 0.10

p_2d <- ggplot(landscape_data, aes(x = proximity_m, y = otu_richness)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  {if (show_line_2d)
    geom_smooth(method = "glm.nb", formula = y ~ x,
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_color_site() +
  labs(
    title = "D. OTU richness vs neighbor proximity",
    subtitle = sprintf("Proximity effect = %.3f, p %s | Controlling for size",
                       coef_prox_rich$estimate,
                       format_p(coef_prox_rich$p.value)),
    x = "Mean Distance to Neighbors (m)",
    y = "OTU Richness",
    color = "Site"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# ============================================================================
# ANALYSIS 5: SHANNON DIVERSITY vs CORAL SIZE
# ============================================================================

cat("\n5. Testing Shannon diversity ~ coral size...\n")

# Linear model for diversity (continuous response)
model_shannon_size <- lm(
  shannon ~ log_volume + site,
  data = landscape_data %>% filter(shannon > 0)
)

coef_size_div <- extract_model_stats(model_shannon_size, "log_volume")

cat(sprintf("  Size effect on diversity: %.3f, p %s\n",
            coef_size_div$estimate, format_p(coef_size_div$p.value)))

# Panel E: Diversity vs Size
show_line_2e <- coef_size_div$p.value < 0.10

p_2e <- landscape_data %>%
  filter(shannon > 0) %>%
  ggplot(aes(x = volume, y = shannon)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  {if (show_line_2e)
    geom_smooth(method = "lm", formula = y ~ log10(x),
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_x_log10(labels = scales::comma) +
  scale_color_site() +
  labs(
    title = "E. Shannon diversity vs coral size",
    subtitle = sprintf("Size effect = %.2f, p %s",
                       coef_size_div$estimate,
                       format_p(coef_size_div$p.value)),
    x = expression(Coral~Volume~(cm^3~","~log~scale)),
    y = "Shannon Diversity (H')",
    color = "Site"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# ============================================================================
# ANALYSIS 6: SHANNON DIVERSITY vs PROXIMITY
# ============================================================================

cat("\n6. Testing Shannon diversity ~ proximity...\n")

model_shannon_prox <- lm(
  shannon ~ proximity_m + log_volume + site,
  data = landscape_data %>% filter(shannon > 0)
)

coef_prox_div <- extract_model_stats(model_shannon_prox, "proximity_m")

cat(sprintf("  Proximity effect on diversity: %.3f, p %s\n",
            coef_prox_div$estimate, format_p(coef_prox_div$p.value)))

# Panel F: Diversity vs Proximity
show_line_2f <- coef_prox_div$p.value < 0.10

p_2f <- landscape_data %>%
  filter(shannon > 0) %>%
  ggplot(aes(x = proximity_m, y = shannon)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  {if (show_line_2f)
    geom_smooth(method = "lm", formula = y ~ x,
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_color_site() +
  labs(
    title = "F. Shannon diversity vs neighbor proximity",
    subtitle = sprintf("Proximity effect = %.3f, p %s",
                       coef_prox_div$estimate,
                       format_p(coef_prox_div$p.value)),
    x = "Mean Distance to Neighbors (m)",
    y = "Shannon Diversity (H')",
    color = "Site"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# ============================================================================
# COMBINED 6-PANEL FIGURE
# ============================================================================

cat("\nCreating combined 6-panel figure...\n")

fig_landscape <- (p_2a | p_2b) /
  (p_2c | p_2d) /
  (p_2e | p_2f) +
  plot_annotation(
    title = "Landscape characteristics shape CAFI communities",
    subtitle = "Effects of coral size (left) and neighbor proximity (right) on abundance, richness, and diversity",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray30")
    )
  )

# Save to neighborhood folder
ggsave(file.path(fig_dir, "landscape_effects_6panel.png"), fig_landscape,
       width = 14, height = 14, dpi = 300, bg = "white")
cat("  - Saved:", file.path(fig_dir, "landscape_effects_6panel.png"), "\n")

# Save to manuscript folder
ggsave(file.path(MANUSCRIPT_DIR, "landscape_effects_6panel.png"), fig_landscape,
       width = 14, height = 14, dpi = 300, bg = "white")
cat("  - Saved:", file.path(MANUSCRIPT_DIR, "landscape_effects_6panel.png"), "\n")

# ============================================================================
# SAVE MODEL RESULTS
# ============================================================================

cat("\nSaving model results...\n")

# Compile all model results
landscape_results <- tibble(
  response = c("Abundance", "Abundance", "Richness", "Richness",
               "Shannon", "Shannon"),
  predictor = c("Coral Size (log)", "Proximity", "Coral Size (log)", "Proximity",
                "Coral Size (log)", "Proximity"),
  estimate = c(coef_size_abund$estimate, coef_prox_abund$estimate,
               coef_size_rich$estimate, coef_prox_rich$estimate,
               coef_size_div$estimate, coef_prox_div$estimate),
  std_error = c(coef_size_abund$se, coef_prox_abund$se,
                coef_size_rich$se, coef_prox_rich$se,
                coef_size_div$se, coef_prox_div$se),
  conf_low = c(coef_size_abund$conf.low, coef_prox_abund$conf.low,
               coef_size_rich$conf.low, coef_prox_rich$conf.low,
               coef_size_div$conf.low, coef_prox_div$conf.low),
  conf_high = c(coef_size_abund$conf.high, coef_prox_abund$conf.high,
                coef_size_rich$conf.high, coef_prox_rich$conf.high,
                coef_size_div$conf.high, coef_prox_div$conf.high),
  p_value = c(coef_size_abund$p.value, coef_prox_abund$p.value,
              coef_size_rich$p.value, coef_prox_rich$p.value,
              coef_size_div$p.value, coef_prox_div$p.value),
  significant = c(coef_size_abund$significant, coef_prox_abund$significant,
                  coef_size_rich$significant, coef_prox_rich$significant,
                  coef_size_div$significant, coef_prox_div$significant)
)

write_csv(landscape_results, file.path(TABLES_DIR, "landscape_model_results.csv"))

# Save model objects
saveRDS(list(
  abundance_size = model_abundance_size,
  abundance_prox = model_abundance_prox,
  richness_size = model_richness_size,
  richness_prox = model_richness_prox,
  shannon_size = model_shannon_size,
  shannon_prox = model_shannon_prox
), file.path(OBJECTS_DIR, "landscape_models.rds"))

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("Landscape Effects Analysis Summary\n")
cat("========================================\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. SIZE-ABUNDANCE SCALING:\n")
cat(sprintf("   Exponent = %.2f [%.2f, %.2f]\n",
            coef_size_abund$estimate, coef_size_abund$conf.low, coef_size_abund$conf.high))
if (coef_size_abund$significant) {
  if (coef_size_abund$conf.high < 1) {
    cat("   -> SUBLINEAR: Larger corals have lower CAFI density\n")
  } else {
    cat("   -> Significant positive size effect\n")
  }
}

cat("\n2. SIZE-RICHNESS:\n")
cat(sprintf("   Effect = %.2f, p %s\n",
            coef_size_rich$estimate, format_p(coef_size_rich$p.value)))
if (coef_size_rich$significant) {
  cat("   -> Larger corals support more species\n")
}

cat("\n3. SIZE-DIVERSITY:\n")
cat(sprintf("   Effect = %.2f, p %s\n",
            coef_size_div$estimate, format_p(coef_size_div$p.value)))

cat("\n4. PROXIMITY EFFECTS:\n")
if (coef_prox_abund$significant) {
  if (coef_prox_abund$estimate > 0) {
    cat("   Abundance: Isolated corals host more CAFI\n")
  } else {
    cat("   Abundance: Clustered corals host more CAFI\n")
  }
} else {
  cat("   Abundance: No significant proximity effect\n")
}
if (coef_prox_rich$significant) {
  cat(sprintf("   Richness: Effect = %.2f\n", coef_prox_rich$estimate))
} else {
  cat("   Richness: No significant proximity effect\n")
}
if (coef_prox_div$significant) {
  cat(sprintf("   Diversity: Effect = %.2f\n", coef_prox_div$estimate))
} else {
  cat("   Diversity: No significant proximity effect\n")
}

cat("\nFigures saved to:\n")
cat("  -", file.path(fig_dir, "landscape_effects_6panel.png"), "\n")
cat("  -", file.path(MANUSCRIPT_DIR, "landscape_effects_6panel.png"), "\n")

cat("\n\u2705 Script 19 complete!\n")
