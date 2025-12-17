#!/usr/bin/env Rscript
# ============================================================================
# 02_fig2_landscape_scaling.R - FIGURE 2: Landscape -> CAFI Scaling
# ============================================================================
#
# PURPOSE: Test how coral size and proximity shape total CAFI abundance,
# species richness, and Shannon diversity.
#
# HYPOTHESES:
#   H1: Larger corals host more CAFI (positive size-abundance slope)
#   H2: Isolated corals host more CAFI per unit size (propagule redirection)
#   H3: Diversity increases with coral size but may decrease with proximity
#
# THEORETICAL BACKGROUND:
#   Propagule redirection theory predicts that occupant abundance scales
#   nonlinearly with habitat amount (power-law with exponent < 1). This causes
#   occupant density to decrease as habitat increases - larger corals intercept
#   more propagules but at lower density per unit volume.
#
# STATISTICAL APPROACHES:
#   - Abundance: Negative binomial GLM with log(volume) and proximity
#   - Richness: Poisson/NB GLM with log(volume) and proximity
#   - Diversity: Linear model with volume and proximity
#   - All models include site as random/fixed effect
#   - Power-law scaling: log-log regression with 95% CI on exponent
#
# PLOTTING RULE: If relationship is nonsignificant (p > 0.10), show points
# only without regression line.
#
# Author: CAFI Analysis Pipeline
# Date: 2025-12-03
# ============================================================================

cat("\n========================================\n")
cat("FIGURE 2: Landscape -> CAFI Scaling\n")
cat("========================================\n\n")

# Load setup and data
source(here::here("scripts/publication/00_setup.R"))
source(here::here("scripts/publication/01_load_data.R"))

# Output directory for this figure
FIG_DIR <- FIGURE_DIRS$fig2
cat("Outputs will be saved to:", FIG_DIR, "\n\n")

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

cat("  - Analysis sample:", nrow(landscape_data), "corals\n")
cat("  - Volume range:", round(min(landscape_data$volume)), "-",
    round(max(landscape_data$volume)), "cm3\n")
cat("  - Proximity range:", round(min(landscape_data$proximity_m), 2), "-",
    round(max(landscape_data$proximity_m), 2), "m\n\n")

# ============================================================================
# ANALYSIS 1: CAFI ABUNDANCE vs CORAL SIZE (Power-law Scaling)
# ============================================================================

cat("1. Testing CAFI abundance ~ coral size (power-law scaling)...\n")

# Model: log(abundance) ~ log(volume) + site
# Power-law: N = a * V^b  ->  log(N) = log(a) + b*log(V)
# Expectation: b ~ 0.75-0.85 for 3D habitat (sublinear scaling)

# Negative binomial GLM (appropriate for overdispersed counts)
model_abundance_size <- glm.nb(
  total_cafi ~ log_volume + site,
  data = landscape_data
)

# Extract scaling exponent
coef_size_abund <- extract_model_stats(model_abundance_size, "log_volume")

cat(sprintf("  Scaling exponent (log-log slope): %.3f (SE = %.3f)\n",
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

# Figure 2A: Size-Abundance Scaling
show_line_2a <- coef_size_abund$p.value < 0.10

p_2a <- ggplot(landscape_data, aes(x = volume, y = total_cafi)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  {if (show_line_2a)
    geom_smooth(method = "glm.nb", formula = y ~ log10(x),
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_x_log10(labels = scales::comma) +
  scale_y_sqrt(breaks = c(0, 25, 100, 225, 400)) +
  scale_color_site() +
  labs(
    title = "A. CAFI Abundance Scales with Coral Size",
    subtitle = sprintf("Power-law exponent = %.2f [%.2f, %.2f], p %s",
                       coef_size_abund$estimate,
                       coef_size_abund$conf.low,
                       coef_size_abund$conf.high,
                       format_p(coef_size_abund$p.value)),
    x = expression(Coral~Volume~(cm^3~","~log~scale)),
    y = "CAFI Abundance (sqrt scale)",
    color = "Site"
  ) +
  theme_publication()

save_pub_fig(p_2a, "fig2a_size_abundance_scaling.png", dir = FIG_DIR)

# ============================================================================
# ANALYSIS 2: CAFI ABUNDANCE vs PROXIMITY (Propagule Redirection)
# ============================================================================

cat("\n2. Testing CAFI abundance ~ proximity (propagule redirection)...\n")

# Model: abundance ~ proximity + volume + site
# Prediction: positive effect of isolation (more distant = more CAFI/unit)
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

# Figure 2B: Proximity effect
show_line_2b <- coef_prox_abund$p.value < 0.10

p_2b <- ggplot(landscape_data, aes(x = proximity_m, y = total_cafi)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  {if (show_line_2b)
    geom_smooth(method = "glm.nb", formula = y ~ x,
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_y_sqrt(breaks = c(0, 25, 100, 225, 400)) +
  scale_color_site() +
  labs(
    title = "B. CAFI Abundance vs Neighbor Proximity",
    subtitle = sprintf("Proximity effect = %.3f, p %s | Controlling for coral size",
                       coef_prox_abund$estimate,
                       format_p(coef_prox_abund$p.value)),
    x = "Mean Distance to Neighbors (m)",
    y = "CAFI Abundance (sqrt scale)",
    color = "Site"
  ) +
  theme_publication()

save_pub_fig(p_2b, "fig2b_proximity_abundance.png", dir = FIG_DIR)

# ============================================================================
# ANALYSIS 3: SPECIES RICHNESS vs CORAL SIZE
# ============================================================================

cat("\n3. Testing OTU richness ~ coral size...\n")

# Poisson GLM for count data (richness)
model_richness_size <- glm(
  otu_richness ~ log_volume + site,
  data = landscape_data,
  family = poisson()
)

# Check for overdispersion
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

# Figure 2C: Richness vs Size
show_line_2c <- coef_size_rich$p.value < 0.10

p_2c <- ggplot(landscape_data, aes(x = volume, y = otu_richness)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  {if (show_line_2c)
    geom_smooth(method = "glm", formula = y ~ log10(x),
                method.args = list(family = poisson()),
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_x_log10(labels = scales::comma) +
  scale_color_site() +
  labs(
    title = "C. OTU Richness Scales with Coral Size",
    subtitle = sprintf("Size effect = %.2f, p %s",
                       coef_size_rich$estimate,
                       format_p(coef_size_rich$p.value)),
    x = expression(Coral~Volume~(cm^3~","~log~scale)),
    y = "OTU Richness",
    color = "Site"
  ) +
  theme_publication()

save_pub_fig(p_2c, "fig2c_size_richness.png", dir = FIG_DIR)

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

# Figure 2D: Richness vs Proximity
show_line_2d <- coef_prox_rich$p.value < 0.10

p_2d <- ggplot(landscape_data, aes(x = proximity_m, y = otu_richness)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  {if (show_line_2d)
    geom_smooth(method = "glm.nb", formula = y ~ x,
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_color_site() +
  labs(
    title = "D. OTU Richness vs Neighbor Proximity",
    subtitle = sprintf("Proximity effect = %.3f, p %s | Controlling for size",
                       coef_prox_rich$estimate,
                       format_p(coef_prox_rich$p.value)),
    x = "Mean Distance to Neighbors (m)",
    y = "OTU Richness",
    color = "Site"
  ) +
  theme_publication()

save_pub_fig(p_2d, "fig2d_proximity_richness.png", dir = FIG_DIR)

# ============================================================================
# ANALYSIS 5: SHANNON DIVERSITY vs LANDSCAPE
# ============================================================================

cat("\n5. Testing Shannon diversity ~ coral size and proximity...\n")

# Linear model for diversity (continuous response)
model_shannon <- lm(
  shannon ~ log_volume + proximity_m + site,
  data = landscape_data %>% filter(shannon > 0)  # Exclude empty corals
)

coef_size_div <- extract_model_stats(model_shannon, "log_volume")
coef_prox_div <- extract_model_stats(model_shannon, "proximity_m")

cat(sprintf("  Size effect on diversity: %.3f, p %s\n",
            coef_size_div$estimate, format_p(coef_size_div$p.value)))
cat(sprintf("  Proximity effect on diversity: %.3f, p %s\n",
            coef_prox_div$estimate, format_p(coef_prox_div$p.value)))

# Figure 2E: Diversity vs Size
show_line_2e <- coef_size_div$p.value < 0.10

p_2e <- landscape_data %>%
  filter(shannon > 0) %>%
  ggplot(aes(x = volume, y = shannon)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  {if (show_line_2e)
    geom_smooth(method = "lm", formula = y ~ log10(x),
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_x_log10(labels = scales::comma) +
  scale_color_site() +
  labs(
    title = "E. Shannon Diversity vs Coral Size",
    subtitle = sprintf("Size effect = %.2f, p %s",
                       coef_size_div$estimate,
                       format_p(coef_size_div$p.value)),
    x = expression(Coral~Volume~(cm^3~","~log~scale)),
    y = "Shannon Diversity (H')",
    color = "Site"
  ) +
  theme_publication()

save_pub_fig(p_2e, "fig2e_size_diversity.png", dir = FIG_DIR)

# Figure 2F: Diversity vs Proximity
show_line_2f <- coef_prox_div$p.value < 0.10

p_2f <- landscape_data %>%
  filter(shannon > 0) %>%
  ggplot(aes(x = proximity_m, y = shannon)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  {if (show_line_2f)
    geom_smooth(method = "lm", formula = y ~ x,
                color = "black", linewidth = 1.2, se = TRUE)
  } +
  scale_color_site() +
  labs(
    title = "F. Shannon Diversity vs Neighbor Proximity",
    subtitle = sprintf("Proximity effect = %.3f, p %s",
                       coef_prox_div$estimate,
                       format_p(coef_prox_div$p.value)),
    x = "Mean Distance to Neighbors (m)",
    y = "Shannon Diversity (H')",
    color = "Site"
  ) +
  theme_publication()

save_pub_fig(p_2f, "fig2f_proximity_diversity.png", dir = FIG_DIR)

# ============================================================================
# COMBINED FIGURE 2 PANEL
# ============================================================================

cat("\n6. Creating combined Figure 2 panel...\n")

# 3x2 panel layout
fig2_combined <- (p_2a + p_2b) /
  (p_2c + p_2d) /
  (p_2e + p_2f) +
  plot_annotation(
    title = "Figure 2. Landscape Characteristics Shape CAFI Communities",
    subtitle = "Effects of coral size and neighbor proximity on abundance, richness, and diversity",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_pub_fig(fig2_combined, "fig2_landscape_scaling_combined.png",
             dir = FIG_DIR, width = 14, height = 16)

# ============================================================================
# SAVE MODEL RESULTS
# ============================================================================

cat("\n7. Saving model results...\n")

# Compile all Figure 2 model results
fig2_results <- tibble(
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

save_table(fig2_results, "fig2_model_results.csv")
save_object(list(
  abundance_size = model_abundance_size,
  abundance_prox = model_abundance_prox,
  richness_size = model_richness_size,
  richness_prox = model_richness_prox,
  shannon = model_shannon
), "fig2_models.rds")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("Figure 2 Analysis Summary\n")
cat("========================================\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. SIZE-ABUNDANCE SCALING:\n")
cat(sprintf("   Power-law exponent = %.2f [%.2f, %.2f]\n",
            coef_size_abund$estimate, coef_size_abund$conf.low, coef_size_abund$conf.high))
if (coef_size_abund$conf.high < 1) {
  cat("   -> SUBLINEAR: Larger corals have lower CAFI density (per unit volume)\n")
  cat("   -> Supports propagule redirection theory\n")
}

cat("\n2. PROXIMITY EFFECTS:\n")
if (coef_prox_abund$significant) {
  if (coef_prox_abund$estimate > 0) {
    cat("   -> Isolated corals host more CAFI (propagule redirection)\n")
  } else {
    cat("   -> Clustered corals host more CAFI (spillover/facilitation)\n")
  }
} else {
  cat("   -> No significant proximity effect on abundance\n")
}

cat("\n3. DIVERSITY PATTERNS:\n")
if (coef_size_div$significant) {
  cat(sprintf("   -> Diversity increases with coral size (beta = %.2f)\n",
              coef_size_div$estimate))
}
if (coef_prox_div$significant) {
  cat(sprintf("   -> Diversity %s with isolation (beta = %.2f)\n",
              ifelse(coef_prox_div$estimate > 0, "increases", "decreases"),
              coef_prox_div$estimate))
}

cat("\nFigures saved to:", FIG_DIR, "\n")
cat("Tables saved to:", TABLES_DIR, "\n\n")
cat("Figure 2 analysis complete!\n\n")
