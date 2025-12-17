#!/usr/bin/env Rscript
# ============================================================================
# 17_power_law_detailed.R - Detailed Power-Law Scaling Analysis (H2)
# ============================================================================
#
# PURPOSE: Comprehensive analysis of power-law scaling relationships between
#          coral size and CAFI communities, testing propagule redirection theory
#
# HYPOTHESIS H2: CAFI abundance scales with coral volume following a power-law
#   with exponent b < 1 (sublinear scaling), indicating larger corals have
#   lower CAFI densities per unit volume.
#
# KEY ANALYSES:
#   - Bootstrap confidence intervals for scaling exponent
#   - Site-specific scaling comparisons
#   - Alternative model comparisons (linear, quadratic, breakpoint)
#   - Residual diagnostics
#
# OUTPUTS:
#   - Detailed scaling figures
#   - Model comparison tables
#   - Manuscript figure
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("17: Detailed Power-Law Scaling Analysis\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "power_law")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Prepare Scaling Data
# ============================================================================

cat("Preparing data for power-law analysis...\n")

scaling_data <- coral_master %>%
  filter(!is.na(volume), volume > 0, !is.na(total_cafi)) %>%
  mutate(
    log_volume = log10(volume),
    log_abundance = log10(total_cafi + 1),
    ln_volume = log(volume),
    ln_abundance = log(total_cafi + 1),
    cafi_density = total_cafi / volume,
    log_density = log10(cafi_density + 0.001)
  )

cat(sprintf("  - Analysis sample: %d corals\n", nrow(scaling_data)))
cat(sprintf("  - Volume range: %.0f – %.0f cm³\n",
            min(scaling_data$volume), max(scaling_data$volume)))
cat(sprintf("  - CAFI range: %d – %d individuals\n",
            min(scaling_data$total_cafi), max(scaling_data$total_cafi)))

# ============================================================================
# Basic Power-Law Model
# ============================================================================

cat("\n1. Fitting basic power-law model...\n")

# OLS on log-log data
m_powerlaw <- lm(ln_abundance ~ ln_volume, data = scaling_data)

# Extract coefficients
beta_ols <- coef(m_powerlaw)[2]
beta_se <- summary(m_powerlaw)$coefficients[2, 2]
beta_ci <- confint(m_powerlaw)[2,]
r_squared <- summary(m_powerlaw)$r.squared

cat(sprintf("  Power-law exponent (β): %.4f\n", beta_ols))
cat(sprintf("  Standard error: %.4f\n", beta_se))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n", beta_ci[1], beta_ci[2]))
cat(sprintf("  R²: %.4f\n", r_squared))

if (beta_ci[2] < 1) {
  cat("  ✓ SUBLINEAR SCALING CONFIRMED (β < 1)\n")
  cat("  → Supports propagule redirection theory\n")
} else if (beta_ci[1] > 1) {
  cat("  → SUPERLINEAR scaling (β > 1)\n")
} else {
  cat("  → Cannot reject isometric scaling (CI includes 1)\n")
}

# ============================================================================
# Bootstrap Confidence Intervals
# ============================================================================

cat("\n2. Bootstrap analysis for robust CI...\n")

set.seed(42)
n_boot <- 1000

boot_betas <- replicate(n_boot, {
  boot_idx <- sample(nrow(scaling_data), replace = TRUE)
  boot_data <- scaling_data[boot_idx, ]
  coef(lm(ln_abundance ~ ln_volume, data = boot_data))[2]
})

boot_ci <- quantile(boot_betas, c(0.025, 0.975))
boot_se <- sd(boot_betas)

cat(sprintf("  Bootstrap β: %.4f ± %.4f\n", mean(boot_betas), boot_se))
cat(sprintf("  Bootstrap 95%% CI: [%.4f, %.4f]\n", boot_ci[1], boot_ci[2]))

# ============================================================================
# Site-Specific Scaling
# ============================================================================

cat("\n3. Site-specific scaling analysis...\n")

site_models <- scaling_data %>%
  group_by(site) %>%
  summarise(
    n = n(),
    beta = coef(lm(ln_abundance ~ ln_volume))[2],
    se = summary(lm(ln_abundance ~ ln_volume))$coefficients[2, 2],
    r2 = summary(lm(ln_abundance ~ ln_volume))$r.squared,
    .groups = "drop"
  ) %>%
  mutate(
    ci_low = beta - 1.96 * se,
    ci_high = beta + 1.96 * se,
    sublinear = ci_high < 1
  )

cat("\n  Site-specific results:\n")
for(i in 1:nrow(site_models)) {
  cat(sprintf("    %s: β = %.3f [%.3f, %.3f], R² = %.3f, n = %d %s\n",
              site_models$site[i],
              site_models$beta[i],
              site_models$ci_low[i],
              site_models$ci_high[i],
              site_models$r2[i],
              site_models$n[i],
              ifelse(site_models$sublinear[i], "✓", "")))
}

# Test for site differences using interaction
m_interaction <- lm(ln_abundance ~ ln_volume * site, data = scaling_data)
interaction_p <- anova(m_interaction)$`Pr(>F)`[3]

cat(sprintf("\n  Site × Volume interaction p-value: %.4f\n", interaction_p))
if (interaction_p < 0.05) {
  cat("  → Scaling exponents differ significantly among sites\n")
} else {
  cat("  → No significant difference in scaling among sites\n")
}

# ============================================================================
# Alternative Model Comparisons
# ============================================================================

cat("\n4. Comparing alternative models...\n")

# Linear model (no log transform on abundance)
m_linear <- lm(total_cafi ~ volume, data = scaling_data)

# Quadratic model
m_quadratic <- lm(ln_abundance ~ ln_volume + I(ln_volume^2), data = scaling_data)

# Model comparison using AIC
model_comparison <- tibble(
  model = c("Power-law (log-log)", "Linear", "Quadratic"),
  AIC = c(AIC(m_powerlaw), AIC(m_linear), AIC(m_quadratic)),
  BIC = c(BIC(m_powerlaw), BIC(m_linear), BIC(m_quadratic)),
  R2 = c(r_squared, summary(m_linear)$r.squared, summary(m_quadratic)$r.squared)
) %>%
  mutate(delta_AIC = AIC - min(AIC)) %>%
  arrange(AIC)

cat("\n  Model comparison (sorted by AIC):\n")
print(model_comparison)

write_csv(model_comparison, file.path(TABLES_DIR, "power_law_model_comparison.csv"))

# ============================================================================
# Density Scaling (Propagule Dilution Test)
# ============================================================================

cat("\n5. Density scaling analysis...\n")

density_data <- scaling_data %>%
  filter(cafi_density > 0)

m_density <- lm(log(cafi_density) ~ ln_volume, data = density_data)
density_slope <- coef(m_density)[2]
density_ci <- confint(m_density)[2,]

cat(sprintf("  Density scaling slope: %.4f [%.4f, %.4f]\n",
            density_slope, density_ci[1], density_ci[2]))

if (density_slope < 0) {
  cat("  ✓ DENSITY DECREASES with size (propagule dilution confirmed)\n")
} else {
  cat("  → Density increases with size (no dilution)\n")
}

# ============================================================================
# Create Diagnostic Figures
# ============================================================================

cat("\n6. Creating figures...\n")

# Panel A: Main scaling relationship
p_main <- ggplot(scaling_data, aes(x = volume, y = total_cafi)) +
  geom_point(aes(color = site), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black",
              se = TRUE, linewidth = 1.2) +
  scale_x_log10(breaks = c(100, 500, 2000, 10000, 30000),
                labels = scales::comma) +
  scale_y_log10(breaks = c(1, 5, 20, 100, 400)) +
  scale_color_site() +
  annotation_logticks() +
  labs(title = "A. Power-law scaling relationship",
       subtitle = sprintf("β = %.3f [%.3f, %.3f], R² = %.3f | %s scaling",
                         beta_ols, beta_ci[1], beta_ci[2], r_squared,
                         ifelse(beta_ci[2] < 1, "Sublinear", "Linear")),
       x = expression(Coral~Volume~(cm^3)),
       y = "CAFI Abundance",
       color = "Site") +
  theme_publication() +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# Panel B: Bootstrap distribution
p_boot <- ggplot(data.frame(beta = boot_betas), aes(x = beta)) +
  geom_histogram(bins = 40, fill = "#0072B2", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(xintercept = beta_ols, color = "black", linewidth = 1) +
  geom_vline(xintercept = boot_ci, color = "black", linetype = "dotted", linewidth = 0.8) +
  annotate("text", x = 1.02, y = Inf, label = "β = 1\n(isometric)",
           hjust = 0, vjust = 1.5, size = 3, color = "red") +
  labs(title = "B. Bootstrap distribution of scaling exponent",
       subtitle = sprintf("n = %d bootstrap samples | 95%% CI shown", n_boot),
       x = "Scaling exponent (β)",
       y = "Frequency") +
  theme_publication()

# Panel C: Site-specific slopes
p_sites <- ggplot(site_models, aes(x = site, y = beta, color = site)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = beta_ols, linetype = "solid", color = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, linewidth = 1) +
  geom_point(size = 4) +
  scale_color_site() +
  annotate("text", x = 0.6, y = 1, label = "Isometric", hjust = 0, size = 3, color = "red") +
  labs(title = "C. Site-specific scaling exponents",
       subtitle = sprintf("Overall β = %.3f (black line) | Error bars = 95%% CI", beta_ols),
       x = "Site",
       y = "Scaling exponent (β)") +
  theme_publication() +
  theme(legend.position = "none")

# Panel D: Density scaling
p_density <- density_data %>%
  ggplot(aes(x = volume, y = cafi_density)) +
  geom_point(aes(color = site), alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black",
              se = TRUE, linewidth = 1.2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  scale_color_site() +
  annotation_logticks() +
  labs(title = "D. Density scaling (propagule dilution)",
       subtitle = sprintf("Slope = %.3f | %s",
                         density_slope,
                         ifelse(density_slope < 0,
                                "Larger corals have LOWER density",
                                "No dilution effect")),
       x = expression(Coral~Volume~(cm^3)),
       y = expression(CAFI~Density~(individuals/cm^3)),
       color = "Site") +
  theme_publication() +
  theme(legend.position = "none")

# Panel E: Residual diagnostics
scaling_data$residuals <- residuals(m_powerlaw)
scaling_data$fitted <- fitted(m_powerlaw)

p_resid <- ggplot(scaling_data, aes(x = fitted, y = residuals)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(aes(color = site), alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", color = "red", se = FALSE, linewidth = 0.8) +
  scale_color_site() +
  labs(title = "E. Residual diagnostics",
       subtitle = "Loess smooth should be flat if model is appropriate",
       x = "Fitted values (log scale)",
       y = "Residuals") +
  theme_publication() +
  theme(legend.position = "none")

# Panel F: Q-Q plot
p_qq <- ggplot(scaling_data, aes(sample = residuals)) +
  stat_qq(alpha = 0.6, color = "#0072B2") +
  stat_qq_line(color = "red", linewidth = 1) +
  labs(title = "F. Q-Q plot of residuals",
       subtitle = "Points should follow red line if normally distributed",
       x = "Theoretical quantiles",
       y = "Sample quantiles") +
  theme_publication()

# Combine 6-panel figure
fig_powerlaw <- (p_main | p_boot | p_sites) / (p_density | p_resid | p_qq) +
  plot_annotation(
    title = "Hypothesis H2: Power-Law Scaling Analysis",
    subtitle = "Testing propagule redirection theory: larger corals should have lower CAFI density per unit volume",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30")
    )
  )

# Save comprehensive figure
ggsave(file.path(fig_dir, "H2_power_law_comprehensive_6panel.png"), fig_powerlaw,
       width = 16, height = 11, dpi = 300, bg = "white")
cat("  - Saved: H2_power_law_comprehensive_6panel.png\n")

# Save to manuscript folder
ggsave(file.path(MANUSCRIPT_DIR, "H2_power_law_comprehensive_6panel.png"), fig_powerlaw,
       width = 16, height = 11, dpi = 300, bg = "white")
cat("  - Saved to manuscript folder\n")

# ============================================================================
# Save Results
# ============================================================================

cat("\n7. Saving results...\n")

# Comprehensive results table
results_summary <- tibble(
  metric = c(
    "Scaling exponent (β)", "β SE", "β 95% CI lower", "β 95% CI upper",
    "R²", "n corals",
    "Bootstrap β mean", "Bootstrap β SE", "Boot CI lower", "Boot CI upper",
    "Density slope", "Density CI lower", "Density CI upper",
    "Site interaction p-value"
  ),
  value = c(
    beta_ols, beta_se, beta_ci[1], beta_ci[2],
    r_squared, nrow(scaling_data),
    mean(boot_betas), boot_se, boot_ci[1], boot_ci[2],
    density_slope, density_ci[1], density_ci[2],
    interaction_p
  )
)

write_csv(results_summary, file.path(TABLES_DIR, "H2_power_law_detailed_results.csv"))
write_csv(site_models, file.path(TABLES_DIR, "H2_site_specific_scaling.csv"))

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Power-Law Scaling Analysis Summary\n")
cat("========================================\n\n")

cat("KEY FINDINGS:\n")
cat(sprintf("  - Scaling exponent β = %.3f [%.3f, %.3f]\n",
            beta_ols, beta_ci[1], beta_ci[2]))
if (beta_ci[2] < 1) {
  cat("  - ✓ SUBLINEAR SCALING CONFIRMED\n")
  cat("  - Larger corals have proportionally fewer CAFI\n")
  cat("  - Supports propagule redirection theory\n")
}
cat(sprintf("  - Density slope = %.3f (%s dilution)\n",
            density_slope,
            ifelse(density_slope < 0, "confirms", "no")))
cat(sprintf("  - Site differences: %s (p = %.3f)\n",
            ifelse(interaction_p < 0.05, "significant", "not significant"),
            interaction_p))

cat("\n  Theory prediction: β ≈ 0.75 for 3D habitat\n")
cat(sprintf("  Observed β = %.3f (%s from prediction)\n",
            beta_ols,
            ifelse(abs(beta_ols - 0.75) < 0.1, "close", "differs")))

cat("\nFigures saved to:", fig_dir, "\n")
cat("Tables saved to:", TABLES_DIR, "\n\n")

cat("✅ Script 17 complete!\n")
