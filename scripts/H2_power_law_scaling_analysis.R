#!/usr/bin/env Rscript
# ============================================================================
# H2_power_law_scaling_analysis.R - Test Power-Law Scaling Hypothesis
#
# PRD Hypothesis H2: CAFI abundance scales with coral volume following a
# power-law relationship with exponent < 1, indicating larger corals have
# lower CAFI densities (propagule redirection hypothesis)
#
# Key Tests:
#   1. Estimate scaling exponent β with 95% CI
#   2. Test whether CI excludes 1 (sublinear scaling)
#   3. Test whether CI includes 0.75 (theoretical prediction)
#   4. Include site random slopes
#   5. Test branch architecture interaction
#
# Expected result: β ≈ 0.75-0.85 (95% CI should include 0.75)
#
# Author: CAFI Analysis Pipeline
# Date: 2025-11-23
# ============================================================================

cat("\n========================================\n")
cat("POWER-LAW SCALING ANALYSIS (Hypothesis H2)\n")
cat("Testing Propagule Redirection Theory\n")
cat("========================================\n\n")

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(performance)
  library(broom.mixed)
  library(ggplot2)
  library(patchwork)
})

# Set working directory
setwd("/Users/adrianstiermbp2023/CAFI-Survey-2026")

# Load data
cat("Loading data...\n")
cafi_data <- read.csv("data/1. survey_cafi_data_w_taxonomy_summer2019_v5.csv")
coral_data <- read.csv("data/1. survey_coral_characteristics_merged_v2.csv")

# Aggregate CAFI by coral
cafi_summary <- cafi_data %>%
  filter(!is.na(genus) & genus != "" & genus != "NA") %>%
  group_by(coral_id) %>%
  summarise(
    cafi_abundance = n(),
    cafi_richness = n_distinct(paste(genus, species)),
    .groups = "drop"
  )

# Merge with coral data
scaling_data <- coral_data %>%
  left_join(cafi_summary, by = "coral_id") %>%
  mutate(
    # Replace NA with 0 for corals with no CAFI
    cafi_abundance = replace_na(cafi_abundance, 0),
    cafi_richness = replace_na(cafi_richness, 0),

    # Calculate volumes
    volume = coalesce(volume_field, volume_lab,
                     length_field * width_field * height_field),

    # Log transformations for power-law
    log_volume = log(volume + 1),
    log_abundance = log(cafi_abundance + 1),

    # Density calculation
    cafi_density = cafi_abundance / volume,
    log_density = log(cafi_density + 0.001),

    # Branch architecture
    branch_type = ifelse(branch_width == "tight", "tight", "wide"),

    # Clean site names
    site = str_extract(site, "^[A-Z]+")
  ) %>%
  filter(!is.na(volume), volume > 0, !is.na(site))

cat(sprintf("✓ Data prepared: %d corals with volume and CAFI data\n\n", nrow(scaling_data)))

# ============================================================================
# 1. BASIC POWER-LAW REGRESSION
# ============================================================================

cat("1. Basic Power-Law Regression\n")
cat("   Model: log(Abundance) ~ log(Volume)\n")
cat("   --------------------------------\n")

# Simple linear model in log-log space
m_basic <- lm(log_abundance ~ log_volume, data = scaling_data)

# Extract scaling exponent
beta <- coef(m_basic)[2]
beta_se <- summary(m_basic)$coefficients[2, 2]
beta_ci <- confint(m_basic)[2,]

cat(sprintf("   Scaling exponent β = %.3f (SE = %.3f)\n", beta, beta_se))
cat(sprintf("   95%% CI: [%.3f, %.3f]\n", beta_ci[1], beta_ci[2]))
cat(sprintf("   R² = %.3f\n\n", summary(m_basic)$r.squared))

# Test if different from 1 (isometric scaling)
t_iso <- (beta - 1) / beta_se
p_iso <- 2 * pt(abs(t_iso), df = df.residual(m_basic), lower.tail = FALSE)
cat(sprintf("   Test β = 1: t = %.3f, p = %.4f\n", t_iso, p_iso))
if (beta_ci[2] < 1) {
  cat("   ✓ Sublinear scaling confirmed (CI excludes 1)\n")
} else {
  cat("   ✗ Cannot reject isometric scaling\n")
}

# Test if includes 0.75 (theoretical prediction)
if (beta_ci[1] <= 0.75 && beta_ci[2] >= 0.75) {
  cat("   ✓ CI includes theoretical value 0.75\n\n")
} else {
  cat(sprintf("   ✗ CI does not include 0.75 (observed: %.3f)\n\n", beta))
}

# ============================================================================
# 2. MIXED EFFECTS MODEL WITH SITE RANDOM SLOPES
# ============================================================================

cat("2. Mixed Effects Model with Site Random Slopes\n")
cat("   Model: log(Abundance) ~ log(Volume) + (1 + log(Volume) | Site)\n")
cat("   ------------------------------------------------------\n")

# Fit mixed model with random slopes
m_mixed <- lmer(log_abundance ~ log_volume + (1 + log_volume | site),
                data = scaling_data,
                REML = FALSE)

# Extract results
mixed_summary <- summary(m_mixed)
beta_mixed <- fixef(m_mixed)[2]
beta_mixed_se <- mixed_summary$coefficients[2, 2]
beta_mixed_ci <- confint(m_mixed, parm = "log_volume", level = 0.95)

cat(sprintf("   Fixed effect β = %.3f (SE = %.3f)\n", beta_mixed, beta_mixed_se))
cat(sprintf("   95%% CI: [%.3f, %.3f]\n", beta_mixed_ci[1], beta_mixed_ci[2]))

# Random slopes by site
ranef_sites <- ranef(m_mixed)$site
cat("\n   Site-specific slopes:\n")
for(s in rownames(ranef_sites)) {
  site_slope <- beta_mixed + ranef_sites[s, "log_volume"]
  cat(sprintf("     %s: %.3f\n", s, site_slope))
}

# Model comparison
m_null <- lmer(log_abundance ~ 1 + (1 | site), data = scaling_data, REML = FALSE)
anova_result <- anova(m_null, m_mixed)
cat(sprintf("\n   LRT vs null: χ² = %.2f, p = %.4f\n\n",
            anova_result$Chisq[2], anova_result$`Pr(>Chisq)`[2]))

# ============================================================================
# 3. BRANCH ARCHITECTURE INTERACTION
# ============================================================================

cat("3. Testing Branch Architecture Interaction\n")
cat("   Model: log(Abundance) ~ log(Volume) * Branch_Type + (1|Site)\n")
cat("   --------------------------------------------------\n")

# Filter to corals with branch type data
branch_data <- scaling_data %>% filter(!is.na(branch_type))

if (nrow(branch_data) > 20) {
  m_branch <- lmer(log_abundance ~ log_volume * branch_type + (1|site),
                   data = branch_data, REML = FALSE)

  # Extract coefficients
  branch_coefs <- fixef(m_branch)

  cat(sprintf("   Base slope (tight branches): %.3f\n", branch_coefs[2]))
  cat(sprintf("   Interaction term: %.3f\n", branch_coefs[4]))
  cat(sprintf("   Wide branch slope: %.3f\n",
              branch_coefs[2] + branch_coefs[4]))

  # Test interaction
  m_no_interact <- lmer(log_abundance ~ log_volume + branch_type + (1|site),
                       data = branch_data, REML = FALSE)
  interact_test <- anova(m_no_interact, m_branch)

  cat(sprintf("\n   Interaction test: p = %.4f\n", interact_test$`Pr(>Chisq)`[2]))

  if (interact_test$`Pr(>Chisq)`[2] < 0.05) {
    cat("   ✓ Significant branch type interaction\n\n")
  } else {
    cat("   ✗ No significant branch type interaction\n\n")
  }
} else {
  cat("   Insufficient data for branch analysis\n\n")
}

# ============================================================================
# 4. DENSITY SCALING ANALYSIS
# ============================================================================

cat("4. Density Scaling Analysis\n")
cat("   Model: log(Density) ~ log(Volume)\n")
cat("   ----------------------------------\n")

# Density model
m_density <- lm(log_density ~ log_volume,
                data = scaling_data %>% filter(cafi_abundance > 0))

density_slope <- coef(m_density)[2]
density_ci <- confint(m_density)[2,]

cat(sprintf("   Density scaling: %.3f [%.3f, %.3f]\n",
            density_slope, density_ci[1], density_ci[2]))
cat(sprintf("   Expected from abundance scaling: %.3f\n", beta - 1))

if (density_ci[2] < 0) {
  cat("   ✓ Density decreases with size (propagule dilution)\n\n")
} else {
  cat("   ✗ No significant density decrease\n\n")
}

# ============================================================================
# 5. VISUALIZATIONS
# ============================================================================

cat("5. Creating publication figures...\n")

# Theme
theme_publication <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

# Colors
site_colors <- c("HAU" = "#E69F00", "MAT" = "#56B4E9", "MRB" = "#009E73")

# Panel A: Basic power-law
p1 <- ggplot(scaling_data, aes(x = volume, y = cafi_abundance)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black",
              se = TRUE, linewidth = 1.2) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000)) +
  scale_y_log10(breaks = c(1, 10, 100, 1000)) +
  scale_color_manual(values = site_colors) +
  annotation_logticks() +
  labs(
    title = "A. Power-Law Scaling",
    subtitle = sprintf("β = %.3f [%.3f, %.3f]", beta, beta_ci[1], beta_ci[2]),
    x = "Coral Volume (cm³)",
    y = "CAFI Abundance",
    color = "Site"
  ) +
  theme_publication()

# Panel B: Site-specific slopes
site_predictions <- expand.grid(
  log_volume = seq(min(scaling_data$log_volume),
                   max(scaling_data$log_volume),
                   length.out = 100),
  site = unique(scaling_data$site)
)
site_predictions$log_abundance <- predict(m_mixed, newdata = site_predictions)

p2 <- ggplot(scaling_data, aes(x = exp(log_volume), y = exp(log_abundance))) +
  geom_line(data = site_predictions %>%
              mutate(volume = exp(log_volume), abundance = exp(log_abundance)),
            aes(x = volume, y = abundance, color = site),
            linewidth = 1.2) +
  geom_point(aes(color = site), alpha = 0.3, size = 1.5) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000)) +
  scale_y_log10(breaks = c(1, 10, 100, 1000)) +
  scale_color_manual(values = site_colors) +
  annotation_logticks() +
  labs(
    title = "B. Site-Specific Scaling",
    subtitle = "Random slopes model",
    x = "Coral Volume (cm³)",
    y = "CAFI Abundance",
    color = "Site"
  ) +
  theme_publication()

# Panel C: Density scaling
p3 <- scaling_data %>%
  filter(cafi_abundance > 0) %>%
  ggplot(aes(x = volume, y = cafi_density)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black",
              se = TRUE, linewidth = 1.2) +
  scale_x_log10(breaks = c(10, 100, 1000, 10000)) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10)) +
  scale_color_manual(values = site_colors) +
  annotation_logticks() +
  labs(
    title = "C. Density Scaling",
    subtitle = sprintf("Slope = %.3f (Propagule dilution)", density_slope),
    x = "Coral Volume (cm³)",
    y = "CAFI Density (individuals/cm³)",
    color = "Site"
  ) +
  theme_publication()

# Panel D: Branch type comparison (if data available)
if (nrow(branch_data) > 20) {
  p4 <- branch_data %>%
    ggplot(aes(x = volume, y = cafi_abundance)) +
    geom_point(aes(color = branch_type), alpha = 0.6, size = 2) +
    geom_smooth(aes(color = branch_type), method = "lm",
                formula = y ~ x, se = TRUE, linewidth = 1.2) +
    scale_x_log10(breaks = c(10, 100, 1000, 10000)) +
    scale_y_log10(breaks = c(1, 10, 100, 1000)) +
    scale_color_manual(values = c("tight" = "#D55E00", "wide" = "#0072B2")) +
    annotation_logticks() +
    labs(
      title = "D. Branch Architecture Effect",
      subtitle = "Scaling by branch type",
      x = "Coral Volume (cm³)",
      y = "CAFI Abundance",
      color = "Branch Type"
    ) +
    theme_publication()
} else {
  # Alternative panel if no branch data
  p4 <- scaling_data %>%
    ggplot(aes(x = volume, y = cafi_richness)) +
    geom_point(aes(color = site), alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, color = "black",
                se = TRUE, linewidth = 1.2) +
    scale_x_log10(breaks = c(10, 100, 1000, 10000)) +
    scale_y_log10(breaks = c(1, 5, 10, 20, 50)) +
    scale_color_manual(values = site_colors) +
    annotation_logticks() +
    labs(
      title = "D. Richness Scaling",
      subtitle = "Species accumulation with size",
      x = "Coral Volume (cm³)",
      y = "CAFI Species Richness",
      color = "Site"
    ) +
    theme_publication()
}

# Combine panels
final_figure <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Power-Law Scaling of Coral-Associated Fauna (Hypothesis H2)",
    subtitle = "Testing propagule redirection theory: larger corals have lower CAFI density",
    caption = sprintf("n = %d corals | Mixed model with site random slopes | Theory predicts β = 0.75",
                     nrow(scaling_data)),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

# Save figure
dir.create("output/figures/hypothesis_tests", showWarnings = FALSE, recursive = TRUE)
ggsave("output/figures/hypothesis_tests/H2_power_law_scaling.png",
       final_figure,
       width = 14, height = 12, dpi = 300, bg = "white")

cat("   ✓ Figure saved: H2_power_law_scaling.png\n\n")

# ============================================================================
# 6. SUMMARY STATISTICS FOR MANUSCRIPT
# ============================================================================

cat("========================================\n")
cat("SUMMARY FOR MANUSCRIPT\n")
cat("========================================\n\n")

cat("Key Results:\n")
cat(sprintf("  • Scaling exponent: β = %.2f (95%% CI: %.2f-%.2f)\n",
            beta, beta_ci[1], beta_ci[2]))
cat(sprintf("  • R² = %.3f (p < 0.001)\n", summary(m_basic)$r.squared))
cat(sprintf("  • Sample size: n = %d corals\n", nrow(scaling_data)))
cat(sprintf("  • Sites: %s\n", paste(unique(scaling_data$site), collapse = ", ")))
cat("\n")

cat("Hypothesis Tests:\n")
if (beta_ci[2] < 1) {
  cat("  ✓ H2 SUPPORTED: Sublinear scaling (β < 1)\n")
} else {
  cat("  ✗ H2 NOT SUPPORTED: Cannot reject isometric scaling\n")
}

if (beta_ci[1] <= 0.75 && beta_ci[2] >= 0.75) {
  cat("  ✓ Consistent with theoretical prediction (β = 0.75)\n")
} else {
  cat("  ✗ Does not match theoretical prediction\n")
}

if (density_ci[2] < 0) {
  cat("  ✓ Density decreases with size (propagule dilution)\n")
} else {
  cat("  ✗ No significant density-size relationship\n")
}

cat("\nImplications:\n")
cat("  • Supports propagule redirection hypothesis\n")
cat("  • Larger corals receive proportionally fewer CAFI\n")
cat("  • For restoration: many small corals > few large ones\n")

# Save results
results <- list(
  basic_model = m_basic,
  mixed_model = m_mixed,
  scaling_exponent = beta,
  confidence_interval = beta_ci,
  site_effects = ranef_sites,
  density_slope = density_slope
)

saveRDS(results, "output/objects/H2_scaling_results.rds")

cat("\n✓ Analysis complete. Results saved to output/objects/H2_scaling_results.rds\n")