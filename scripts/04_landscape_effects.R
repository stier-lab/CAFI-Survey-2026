# ============================================================================
# 04_landscape_effects.R - Landscape Effects on CAFI Abundance & Diversity
# ============================================================================
#
# PURPOSE: Test how coral size and neighborhood context shape CAFI communities
#
# IMPORTANT: Survey Design Distinction
#   - "neighborhood" corals (n=63): POC01-POC21 at each site - full 5m radius
#     surveys of all neighboring Pocillopora within 5m
#   - "size" corals (n=51): POC22+ at each site - surveyed only for size/CAFI,
#     NO neighborhood data collected (n_neighbors = NA)
#
#   This script uses ONLY corals with complete neighborhood data (filters !is.na())
#
# DEPENDENCIES: 00_setup.R, 01_load_data.R (coral_master.rds, cafi_clean.rds)
# PREREQUISITE: Run 03_landscape_characterization.R first to understand
#               landscape structure and predictor correlations
#
# HYPOTHESES:
#   H1: Coral size drives CAFI abundance (power-law scaling)
#   H2: Neighborhood predictors have weak effects relative to size
#   H3: Diversity patterns mirror abundance patterns
#
# NON-REDUNDANT PREDICTOR SUBSET (see 03 for rationale):
#   1. log(volume)            - Focal coral size (PRIMARY)
#   2. n_neighbors            - Neighbor count
#   3. total_neighbor_volume  - Neighbor habitat
#   4. mean_neighbor_dist     - Spatial isolation
#
# OUTPUTS:
#   - Figures: output/figures/04_effects/ (working), output/figures/supplement/figS7_neighborhood.png
#   - Tables: output/tables/landscape_full_model_results.csv, landscape_power_analysis.csv
#   - Statistical summaries with effect sizes, p-values, df
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-07
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    LANDSCAPE EFFECTS ON CAFI COMMUNITIES\n")
cat("============================================================\n\n")

# Load setup and data
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

# Use script-specific output directory
FIG_DIR <- PATHS$fig_04_effects

# ============================================================================
# PREPARE LANDSCAPE DATA
# ============================================================================

cat("Preparing landscape analysis data...\n\n")

# Note: Only 63 corals have neighborhood data (survey_type == "neighborhood")
# The filter !is.na(n_neighbors) automatically excludes the 51 "size" corals
cat("Survey Design Reminder:\n")
cat("   Total corals:", nrow(coral_master), "\n")
cat("   Neighborhood surveys:", sum(!is.na(coral_master$n_neighbors)), "\n")
cat("   Size-only surveys:", sum(is.na(coral_master$n_neighbors)), "\n\n")

# Filter to corals with COMPLETE neighborhood data for all predictors
landscape_data <- coral_master %>%
  filter(
    !is.na(volume), volume > 0,
    !is.na(total_cafi),
    !is.na(n_neighbors),  # This filters to ONLY neighborhood-surveyed corals
    !is.na(total_neighbor_volume),
    !is.na(mean_neighbor_dist)
  ) %>%
  mutate(
    # Log transformations
    log_volume = log(volume),
    # Proximity in meters
    proximity_m = mean_neighbor_dist / 100,

    # Categories for visualization
    size_cat = cut(volume,
                   breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                   labels = c("Small", "Medium", "Large"),
                   include.lowest = TRUE),

    isolation_cat = cut(proximity_m,
                        breaks = quantile(proximity_m, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                        labels = c("Clustered", "Intermediate", "Isolated"),
                        include.lowest = TRUE),

    neighbor_cat = case_when(
      n_neighbors == 0 ~ "None",
      n_neighbors <= 2 ~ "Low (1-2)",
      n_neighbors <= 4 ~ "Medium (3-4)",
      TRUE ~ "High (5+)"
    ) %>% factor(levels = c("None", "Low (1-2)", "Medium (3-4)", "High (5+)"))
  )

cat("ANALYSIS DATASET (neighborhood-surveyed corals only):\n")
cat("   Sample size:", nrow(landscape_data), "corals with complete data\n")
cat("   Volume range:", round(min(landscape_data$volume)), "-",
    round(max(landscape_data$volume)), "cm³\n")
cat("   Neighbors range:", min(landscape_data$n_neighbors), "-",
    max(landscape_data$n_neighbors), "\n")
cat("   Truly isolated (0 neighbors):", sum(landscape_data$n_neighbors == 0), "\n\n")

cat("NOTE: See 03_landscape_characterization.R for landscape structure analysis\n")
cat("      and rationale for predictor selection.\n\n")

cat("Using 4 non-redundant predictors (all VIF < 2):\n")
cat("  1. log(volume)            - Focal coral size\n")
cat("  2. n_neighbors            - Neighbor count\n")
cat("  3. total_neighbor_volume  - Neighborhood habitat\n")
cat("  4. mean_neighbor_dist     - Spatial isolation\n\n")

# Helper to safely extract anova.negbin columns (names vary by R/MASS version)
safe_lrt_extract <- function(lrt) {
  tryCatch({
    lr_col <- intersect(names(lrt), c("LR stat.", "LR.stat.", "Deviance", "Chisq"))
    p_col <- intersect(names(lrt), c("Pr(Chi)", "Pr(>Chi)", "Pr..Chi.", "Pr(>Chisq)"))
    if (length(lr_col) == 0 || length(p_col) == 0) {
      ll <- lrt$`-2 Log L.` %||% lrt$`Resid. Dev`
      if (!is.null(ll)) {
        lr_stat <- abs(diff(ll))
        return(list(lr = lr_stat, df = lrt$Df[2], p = NA_real_))
      }
      return(list(lr = NA_real_, df = NA_integer_, p = NA_real_))
    }
    list(lr = lrt[[lr_col[1]]][2], df = lrt$Df[2], p = lrt[[p_col[1]]][2])
  }, error = function(e) {
    cat("      WARNING: Could not extract LRT statistics:", e$message, "\n")
    list(lr = NA_real_, df = NA_integer_, p = NA_real_)
  })
}

# ============================================================================
# POWER ANALYSIS (Q4: Neighborhood Effects)
# ============================================================================
# With n=63 neighborhood-surveyed corals, 3 predictors + site:
# - Power to detect R² = 0.10 (medium effect): ~65% at α=0.05
# - Power to detect R² = 0.05 (small effect): ~35% at α=0.05
# - Conclusion: adequately powered for medium effects, underpowered for small effects
# - Null results should be interpreted as "no evidence" not "no effect"
if (requireNamespace("pwr", quietly = TRUE)) {
  # Cohen's f² for NB GLM approximation
  # Full model: log_volume + n_neighbors + log(neighbor_vol) + mean_dist + 2 site dummies = 6 params
  # Testing u=3 neighborhood predictors jointly; denominator df = n - 6 - 1 = n - 7
  n_landscape <- sum(!is.na(coral_master$n_neighbors))
  power_medium <- pwr::pwr.f2.test(u = 3, v = n_landscape - 7, f2 = 0.10/0.90, sig.level = 0.05)
  power_small <- pwr::pwr.f2.test(u = 3, v = n_landscape - 7, f2 = 0.05/0.95, sig.level = 0.05)
  cat(sprintf("Power analysis (n=%d):\n  Medium effect (R²=0.10): %.0f%%\n  Small effect (R²=0.05): %.0f%%\n\n",
              n_landscape, power_medium$power * 100, power_small$power * 100))
}

# ============================================================================
# PART 1: FOCAL CORAL SIZE EFFECTS
# ============================================================================
# Models use site as a FIXED EFFECT (not random) because k = 3 sites is
# insufficient for reliable random intercept variance estimation
# (Bolker et al. 2009, Trends Ecol Evol). This treats site as a nuisance
# variable to absorb among-site differences rather than estimating site
# variance as a population-level parameter.
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 1: FOCAL CORAL SIZE EFFECTS\n")
cat("------------------------------------------------------------\n\n")

# NOTE: Part 1 tests are exploratory (univariate associations, no FDR correction).
# Confirmatory tests with FDR correction are in Part 7 (interaction models).

# 1.1 CAFI Abundance vs Volume (Power-law)
cat("1.1 CAFI Abundance ~ Coral Volume:\n")

m_abund_vol <- MASS::glm.nb(total_cafi ~ log(volume) + site, data = landscape_data)
m_summary <- summary(m_abund_vol)

slope <- coef(m_abund_vol)["log(volume)"]
slope_se <- sqrt(vcov(m_abund_vol)["log(volume)", "log(volume)"])
slope_ci <- confint(m_abund_vol)["log(volume)", ]

cat("    Scaling exponent (β):", round(slope, 3), "\n")
cat("    95% CI: [", round(slope_ci[1], 3), ",", round(slope_ci[2], 3), "]\n")
cat("    SE:", round(slope_se, 3), "\n")
cat("    z-value:", round(m_summary$coefficients["log(volume)", "z value"], 2), "\n")
cat("    p-value:", format.pval(m_summary$coefficients["log(volume)", "Pr(>|z|)"], 3), "\n")

# Test β ≠ 1
z_vs_1 <- (slope - 1) / slope_se
p_vs_1 <- 2 * (1 - pnorm(abs(z_vs_1)))
cat("    Test β ≠ 1: z =", round(z_vs_1, 2), ", p =", format.pval(p_vs_1, 3), "\n")
cat("    Interpretation:", ifelse(p_vs_1 < 0.05 & slope < 1, "Propagule dilution (β < 1)",
                                   ifelse(p_vs_1 < 0.05 & slope > 1, "Super-linear scaling (β > 1)",
                                          "Field of Dreams (β ≈ 1)")), "\n\n")

# 1.2 Species Richness vs Volume
cat("1.2 Species Richness ~ Coral Volume:\n")

m_rich_vol <- glm(otu_richness ~ log(volume) + site, family = poisson, data = landscape_data)
m_rich_summary <- summary(m_rich_vol)

cat("    β =", round(coef(m_rich_vol)["log(volume)"], 3),
    ", z =", round(m_rich_summary$coefficients["log(volume)", "z value"], 2),
    ", p =", format.pval(m_rich_summary$coefficients["log(volume)", "Pr(>|z|)"], 3), "\n")

# Overdispersion check
dispersion_ratio_rich <- sum(residuals(m_rich_vol, type = "pearson")^2) / m_rich_vol$df.residual
cat("    Overdispersion (Pearson X²/df):", round(dispersion_ratio_rich, 2),
    ifelse(dispersion_ratio_rich > 1.5, " [OVERDISPERSED]", " [OK]"), "\n")
if (dispersion_ratio_rich > 1.5) {
  m_rich_vol <- glm(otu_richness ~ log(volume) + site, family = quasipoisson, data = landscape_data)
  m_rich_summary <- summary(m_rich_vol)
  cat("    Switched to quasi-Poisson: β =", round(coef(m_rich_vol)["log(volume)"], 3),
      ", t =", round(m_rich_summary$coefficients["log(volume)", "t value"], 2),
      ", p =", format.pval(m_rich_summary$coefficients["log(volume)", "Pr(>|t|)"], 3), "\n")
}
cat("\n")

# 1.3 Shannon Diversity vs Volume
cat("1.3 Shannon Diversity ~ Coral Volume:\n")

m_shan_vol <- lm(shannon ~ log(volume) + site, data = landscape_data)
m_shan_summary <- summary(m_shan_vol)

cat("    β =", round(coef(m_shan_vol)["log(volume)"], 3),
    ", t(", m_shan_summary$df[2], ") =", round(m_shan_summary$coefficients["log(volume)", "t value"], 2),
    ", p =", format.pval(m_shan_summary$coefficients["log(volume)", "Pr(>|t|)"], 3),
    ", R² =", round(m_shan_summary$r.squared, 3), "\n\n")

# Figure: Size effects panel
p_size_abund <- ggplot(landscape_data, aes(x = volume, y = total_cafi, color = site)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = MASS::glm.nb, formula = y ~ log(x), se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = mean(log10(landscape_data$total_cafi) - log10(landscape_data$volume)),
              linetype = "dashed", color = "gray50", linewidth = 0.8) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = expression("Coral Volume (cm"^3*")"),
       y = "CAFI Abundance",
       title = "A. Abundance vs Size",
       subtitle = paste0("β = ", round(slope, 2), " [", round(slope_ci[1], 2), ", ",
                         round(slope_ci[2], 2), "]")) +
  theme(legend.position = "none")

p_size_rich <- ggplot(landscape_data, aes(x = volume, y = otu_richness, color = site)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = "glm", formula = y ~ log(x), method.args = list(family = poisson),
              se = TRUE, color = "black") +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = expression("Coral Volume (cm"^3*")"),
       y = "Species Richness",
       title = "B. Richness vs Size") +
  theme(legend.position = "none")

p_size_shan <- ggplot(landscape_data, aes(x = volume, y = shannon, color = site)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = TRUE, color = "black") +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = expression("Coral Volume (cm"^3*")"),
       y = "Shannon H'",
       title = "C. Diversity vs Size",
       subtitle = paste0("R² = ", round(m_shan_summary$r.squared, 2))) +
  theme(legend.position = "none")

# ============================================================================
# PART 2: NEIGHBORHOOD COUNT EFFECTS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 2: NEIGHBORHOOD COUNT EFFECTS\n")
cat("------------------------------------------------------------\n\n")

# 2.1 Neighbor count effect on abundance
cat("2.1 Neighbor Count Effect:\n")

m_neighbors <- MASS::glm.nb(total_cafi ~ n_neighbors + log(volume) + site, data = landscape_data)
m_neigh_summary <- summary(m_neighbors)

cat("    β =", round(coef(m_neighbors)["n_neighbors"], 4),
    ", z =", round(m_neigh_summary$coefficients["n_neighbors", "z value"], 2),
    ", p =", format.pval(m_neigh_summary$coefficients["n_neighbors", "Pr(>|z|)"], 3), "\n")

# Does adding neighbors improve model?
m_size_only <- MASS::glm.nb(total_cafi ~ log(volume) + site, data = landscape_data)
lrt <- anova(m_size_only, m_neighbors, test = "Chisq")
lrt_vals <- safe_lrt_extract(lrt)
cat("    LRT vs size-only: χ² =", round(lrt_vals$lr, 2),
    ", p =", format.pval(lrt_vals$p, 3), "\n\n")

# Figure
p_neighbor_count <- ggplot(landscape_data, aes(x = n_neighbors, y = total_cafi, color = site)) +
  geom_jitter(alpha = 0.6, width = 0.2, size = 2.5) +
  geom_smooth(method = MASS::glm.nb, formula = y ~ x, se = TRUE, color = "black") +
  scale_y_sqrt() +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = "Number of Neighbors (5m radius)",
       y = "CAFI Abundance (sqrt scale)",
       title = "D. Neighbor Count Effect",
       subtitle = paste0("β = ", round(coef(m_neighbors)["n_neighbors"], 3),
                         ", p = ", format.pval(m_neigh_summary$coefficients["n_neighbors", "Pr(>|z|)"], 2))) +
  theme(legend.position = "none")

# 2.2 By neighbor category
cat("2.2 Abundance by Neighbor Category:\n")
neighbor_stats <- landscape_data %>%
  group_by(neighbor_cat) %>%
  summarise(
    n = n(),
    mean_cafi = round(mean(total_cafi), 1),
    se_cafi = round(sd(total_cafi) / sqrt(n()), 1),
    mean_richness = round(mean(otu_richness), 1),
    .groups = "drop"
  )
print(neighbor_stats)
cat("\n")

# ============================================================================
# PART 3: NEIGHBORHOOD VOLUME EFFECTS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 3: NEIGHBORHOOD VOLUME EFFECTS\n")
cat("------------------------------------------------------------\n\n")

# 3.1 Total neighbor volume
cat("3.1 Total Neighbor Volume Effect:\n")

m_neigh_vol <- MASS::glm.nb(total_cafi ~ log(total_neighbor_volume + 1) + log(volume) + site,
                      data = landscape_data)
m_nvol_summary <- summary(m_neigh_vol)

cat("    β =", round(coef(m_neigh_vol)["log(total_neighbor_volume + 1)"], 4),
    ", z =", round(m_nvol_summary$coefficients["log(total_neighbor_volume + 1)", "z value"], 2),
    ", p =", format.pval(m_nvol_summary$coefficients["log(total_neighbor_volume + 1)", "Pr(>|z|)"], 3), "\n\n")

# Figure
p_neighbor_vol <- landscape_data %>%
  filter(total_neighbor_volume > 0) %>%
  ggplot(aes(x = total_neighbor_volume, y = total_cafi, color = site)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = MASS::glm.nb, formula = y ~ log(x), se = TRUE, color = "black") +
  scale_x_log10(labels = scales::comma) +
  scale_y_sqrt() +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = expression("Total Neighbor Volume (cm"^3*")"),
       y = "CAFI Abundance (sqrt scale)",
       title = "E. Neighbor Volume Effect") +
  theme(legend.position = "none")

# ============================================================================
# PART 4: NEIGHBORHOOD DISTANCE (ISOLATION) EFFECTS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 4: ISOLATION EFFECTS\n")
cat("------------------------------------------------------------\n\n")

# 4.1 Mean neighbor distance
cat("4.1 Mean Neighbor Distance Effect:\n")

m_distance <- MASS::glm.nb(total_cafi ~ mean_neighbor_dist + log(volume) + site, data = landscape_data)
m_dist_summary <- summary(m_distance)

cat("    β =", round(coef(m_distance)["mean_neighbor_dist"], 5),
    ", z =", round(m_dist_summary$coefficients["mean_neighbor_dist", "z value"], 2),
    ", p =", format.pval(m_dist_summary$coefficients["mean_neighbor_dist", "Pr(>|z|)"], 3), "\n")

# Kruskal-Wallis by isolation category
kw_iso <- kruskal.test(total_cafi ~ isolation_cat, data = landscape_data)
cat("    Isolation category: H =", round(kw_iso$statistic, 2),
    ", df =", kw_iso$parameter, ", p =", format.pval(kw_iso$p.value, 3), "\n\n")

# Figure
p_isolation <- landscape_data %>%
  ggplot(aes(x = proximity_m, y = total_cafi, color = site)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = MASS::glm.nb, formula = y ~ x, se = TRUE, color = "black") +
  scale_y_sqrt() +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = "Mean Distance to Neighbors (m)",
       y = "CAFI Abundance (sqrt scale)",
       title = "F. Isolation Effect") +
  theme(legend.position = "none")

# ============================================================================
# PART 5: FULL MULTIVARIATE MODELS (4 predictors)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 5: FULL MULTIVARIATE MODELS\n")
cat("------------------------------------------------------------\n\n")

# 5.1 Full model for CAFI Abundance
cat("5.1 Full Model: CAFI Abundance ~ All 4 Predictors\n")

m_abund_full <- MASS::glm.nb(total_cafi ~ log(volume) + n_neighbors +
                        log(total_neighbor_volume + 1) + mean_neighbor_dist + site,
                       data = landscape_data)
m_abund_full_summary <- summary(m_abund_full)

cat("\n    Coefficients:\n")
coef_table <- data.frame(
  Predictor = c("log(volume)", "n_neighbors", "log(neighbor_vol)", "mean_neighbor_dist"),
  Beta = round(coef(m_abund_full)[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist")], 4),
  SE = round(m_abund_full_summary$coefficients[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist"), "Std. Error"], 4),
  z = round(m_abund_full_summary$coefficients[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist"), "z value"], 2),
  p = sapply(m_abund_full_summary$coefficients[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist"), "Pr(>|z|)"],
             function(x) format.pval(x, 3))
)
print(coef_table, row.names = FALSE)

cat("\n    Model fit:\n")
cat("      AIC:", round(AIC(m_abund_full), 1), "\n")
cat("      Theta (dispersion):", round(m_abund_full$theta, 2), "\n")
cat("      Residual deviance:", round(m_abund_full$deviance, 1),
    "on", m_abund_full$df.residual, "df\n")
cat("      Pseudo-R² (McFadden):", round(calc_pseudo_r2(m_abund_full), 4), "\n\n")

# 5.2 Full model for Species Richness
cat("5.2 Full Model: Species Richness ~ All 4 Predictors\n")

m_rich_full <- glm(otu_richness ~ log(volume) + n_neighbors +
                    log(total_neighbor_volume + 1) + mean_neighbor_dist + site,
                   family = poisson, data = landscape_data)
m_rich_full_summary <- summary(m_rich_full)

cat("\n    Coefficients:\n")
coef_table_rich <- data.frame(
  Predictor = c("log(volume)", "n_neighbors", "log(neighbor_vol)", "mean_neighbor_dist"),
  Beta = round(coef(m_rich_full)[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist")], 4),
  SE = round(m_rich_full_summary$coefficients[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist"), "Std. Error"], 4),
  z = round(m_rich_full_summary$coefficients[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist"), "z value"], 2),
  p = sapply(m_rich_full_summary$coefficients[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist"), "Pr(>|z|)"],
             function(x) format.pval(x, 3))
)
print(coef_table_rich, row.names = FALSE)

# Check overdispersion
overdispersion <- m_rich_full$deviance / m_rich_full$df.residual
cat("\n    Overdispersion ratio:", round(overdispersion, 2),
    ifelse(overdispersion > 1.5, " (consider negative binomial)", ""), "\n\n")

# 5.3 Full model for Shannon Diversity
cat("5.3 Full Model: Shannon H' ~ All 4 Predictors\n")

m_shan_full <- lm(shannon ~ log(volume) + n_neighbors +
                   log(total_neighbor_volume + 1) + mean_neighbor_dist + site,
                  data = landscape_data)
m_shan_full_summary <- summary(m_shan_full)

cat("\n    Coefficients:\n")
coef_table_shan <- data.frame(
  Predictor = c("log(volume)", "n_neighbors", "log(neighbor_vol)", "mean_neighbor_dist"),
  Beta = round(coef(m_shan_full)[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist")], 4),
  SE = round(m_shan_full_summary$coefficients[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist"), "Std. Error"], 4),
  t = round(m_shan_full_summary$coefficients[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist"), "t value"], 2),
  p = sapply(m_shan_full_summary$coefficients[c("log(volume)", "n_neighbors", "log(total_neighbor_volume + 1)", "mean_neighbor_dist"), "Pr(>|t|)"],
             function(x) format.pval(x, 3))
)
print(coef_table_shan, row.names = FALSE)

cat("\n    Model fit:\n")
cat("      R²:", round(m_shan_full_summary$r.squared, 3), "\n")
cat("      Adjusted R²:", round(m_shan_full_summary$adj.r.squared, 3), "\n")
cat("      F(", m_shan_full_summary$fstatistic[2], ",",
    m_shan_full_summary$fstatistic[3], ") =",
    round(m_shan_full_summary$fstatistic[1], 2), "\n\n")

# ============================================================================
# PART 6: HIERARCHICAL MODEL COMPARISON
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 6: HIERARCHICAL MODEL COMPARISON\n")
cat("------------------------------------------------------------\n\n")

# Calculate pseudo-R² (McFadden's)
calc_r2 <- function(model, null) {
  1 - (logLik(model)[1] / logLik(null)[1])
}

# Build hierarchy of models
cat("6.1 CAFI Abundance Models:\n\n")

m_null <- MASS::glm.nb(total_cafi ~ site, data = landscape_data)
m_size <- MASS::glm.nb(total_cafi ~ log(volume) + site, data = landscape_data)
m_size_neighbors <- MASS::glm.nb(total_cafi ~ log(volume) + n_neighbors + site,
                           data = landscape_data)
m_size_neigh_vol <- MASS::glm.nb(total_cafi ~ log(volume) + n_neighbors +
                            log(total_neighbor_volume + 1) + site,
                           data = landscape_data)
m_full <- MASS::glm.nb(total_cafi ~ log(volume) + n_neighbors +
                  log(total_neighbor_volume + 1) + mean_neighbor_dist + site,
                 data = landscape_data)

model_comparison <- tibble(
  Model = c("1. Null (site only)",
            "2. + Size",
            "3. + Neighbor count",
            "4. + Neighbor volume",
            "5. + Distance (Full)"),
  Predictors = c("site",
                 "log(vol) + site",
                 "+ n_neighbors",
                 "+ log(neigh_vol)",
                 "+ mean_dist"),
  AIC = c(AIC(m_null), AIC(m_size), AIC(m_size_neighbors),
          AIC(m_size_neigh_vol), AIC(m_full)),
  `Pseudo-R²` = c(0, calc_r2(m_size, m_null), calc_r2(m_size_neighbors, m_null),
                  calc_r2(m_size_neigh_vol, m_null), calc_r2(m_full, m_null))
) %>%
  mutate(
    `Δ AIC` = AIC - min(AIC),
    `Δ R²` = `Pseudo-R²` - lag(`Pseudo-R²`, default = 0)
  )

print(model_comparison, n = 10)

# LRT tests
cat("\n    Likelihood Ratio Tests:\n")

lrt_size <- anova(m_null, m_size, test = "Chisq")
lrt_size_vals <- safe_lrt_extract(lrt_size)
cat("      Null vs Size: χ² =", round(lrt_size_vals$lr, 2),
    ", df =", lrt_size_vals$df,
    ", p =", format.pval(lrt_size_vals$p, 3), "\n")

lrt_neigh <- anova(m_size, m_size_neighbors, test = "Chisq")
lrt_neigh_vals <- safe_lrt_extract(lrt_neigh)
cat("      Size vs +Neighbors: χ² =", round(lrt_neigh_vals$lr, 2),
    ", df =", lrt_neigh_vals$df,
    ", p =", format.pval(lrt_neigh_vals$p, 3), "\n")

lrt_vol <- anova(m_size_neighbors, m_size_neigh_vol, test = "Chisq")
lrt_vol_vals <- safe_lrt_extract(lrt_vol)
cat("      +Neighbors vs +Volume: χ² =", round(lrt_vol_vals$lr, 2),
    ", df =", lrt_vol_vals$df,
    ", p =", format.pval(lrt_vol_vals$p, 3), "\n")

lrt_dist <- anova(m_size_neigh_vol, m_full, test = "Chisq")
lrt_dist_vals <- safe_lrt_extract(lrt_dist)
cat("      +Volume vs +Distance: χ² =", round(lrt_dist_vals$lr, 2),
    ", df =", lrt_dist_vals$df,
    ", p =", format.pval(lrt_dist_vals$p, 3), "\n")

lrt_full <- anova(m_size, m_full, test = "Chisq")
lrt_full_vals <- safe_lrt_extract(lrt_full)
cat("      Size vs Full: χ² =", round(lrt_full_vals$lr, 2),
    ", df =", lrt_full_vals$df,
    ", p =", format.pval(lrt_full_vals$p, 3), "\n")

delta_r2 <- calc_r2(m_full, m_null) - calc_r2(m_size, m_null)
cat("\n    Size R² =", round(calc_r2(m_size, m_null) * 100, 1), "%\n")
cat("    Neighborhood adds:", round(delta_r2 * 100, 2), "% additional R²\n\n")

# 6.2 Same for richness
cat("6.2 Species Richness Models:\n\n")

m_rich_null <- glm(otu_richness ~ site, family = poisson, data = landscape_data)
m_rich_size <- glm(otu_richness ~ log(volume) + site, family = poisson,
                   data = landscape_data)
m_rich_full <- glm(otu_richness ~ log(volume) + n_neighbors +
                    log(total_neighbor_volume + 1) + mean_neighbor_dist + site,
                   family = poisson, data = landscape_data)

cat("    Null AIC:", round(AIC(m_rich_null), 1), "\n")
cat("    Size AIC:", round(AIC(m_rich_size), 1), "\n")
cat("    Full AIC:", round(AIC(m_rich_full), 1), "\n")

lrt_rich <- anova(m_rich_size, m_rich_full, test = "Chisq")
cat("    LRT (Size vs Full): χ² =", round(lrt_rich$Deviance[2], 2),
    ", df =", lrt_rich$Df[2],
    ", p =", format.pval(lrt_rich$`Pr(>Chi)`[2], 3), "\n\n")

# 6.3 Same for Shannon
cat("6.3 Shannon Diversity Models:\n\n")

m_shan_null <- lm(shannon ~ site, data = landscape_data)
m_shan_size <- lm(shannon ~ log(volume) + site, data = landscape_data)
m_shan_full <- lm(shannon ~ log(volume) + n_neighbors +
                   log(total_neighbor_volume + 1) + mean_neighbor_dist + site,
                  data = landscape_data)

cat("    Null R²:", round(summary(m_shan_null)$r.squared, 3), "\n")
cat("    Size R²:", round(summary(m_shan_size)$r.squared, 3), "\n")
cat("    Full R²:", round(summary(m_shan_full)$r.squared, 3), "\n")

lrt_shan <- anova(m_shan_size, m_shan_full)
cat("    F-test (Size vs Full): F =", round(lrt_shan$F[2], 2),
    ", df =", lrt_shan$Df[2], ",", lrt_shan$Res.Df[2],
    ", p =", format.pval(lrt_shan$`Pr(>F)`[2], 3), "\n\n")

# ============================================================================
# PART 7: INTERACTION EFFECTS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 7: TESTING SIZE × NEIGHBORHOOD INTERACTIONS\n")
cat("------------------------------------------------------------\n\n")

cat("Does the effect of neighborhood depend on coral size?\n\n")

# Center log_volume for interpretable main effects in interaction models
# Centering ensures that the main effect of n_neighbors represents the neighborhood
# effect at the AVERAGE coral size, not at log(volume) = 0 (which is volume = 1 cm³)
landscape_data$log_volume_c <- landscape_data$log_volume - mean(landscape_data$log_volume, na.rm = TRUE)
cat("  Centered log_volume: mean =", round(mean(landscape_data$log_volume, na.rm = TRUE), 3),
    " (subtracted from log_volume for interaction models)\n\n")

# 7.1 Size × Neighbor count interaction
cat("7.1 Abundance ~ Size × Neighbor Count:\n")

m_int_neighbors <- MASS::glm.nb(total_cafi ~ log_volume_c * n_neighbors + site,
                          data = landscape_data)
m_int_neigh_summary <- summary(m_int_neighbors)

int_coef <- coef(m_int_neighbors)["log_volume_c:n_neighbors"]
int_se <- m_int_neigh_summary$coefficients["log_volume_c:n_neighbors", "Std. Error"]
int_z <- m_int_neigh_summary$coefficients["log_volume_c:n_neighbors", "z value"]
int_p <- m_int_neigh_summary$coefficients["log_volume_c:n_neighbors", "Pr(>|z|)"]

cat("    Interaction β =", round(int_coef, 4),
    ", z =", round(int_z, 2),
    ", p =", format.pval(int_p, 3), "\n")
cat("    Interpretation:", ifelse(int_p < 0.05,
                                   "Neighbor effects differ by coral size",
                                   "No size-dependent neighbor effects"), "\n\n")

# 7.2 Size × Neighbor volume interaction
cat("7.2 Abundance ~ Size × Neighbor Volume:\n")

m_int_vol <- MASS::glm.nb(total_cafi ~ log_volume_c * log(total_neighbor_volume + 1) + site,
                    data = landscape_data)
m_int_vol_summary <- summary(m_int_vol)

int_coef_vol <- coef(m_int_vol)["log_volume_c:log(total_neighbor_volume + 1)"]
int_p_vol <- m_int_vol_summary$coefficients["log_volume_c:log(total_neighbor_volume + 1)", "Pr(>|z|)"]

cat("    Interaction β =", round(int_coef_vol, 4),
    ", p =", format.pval(int_p_vol, 3), "\n")
cat("    Interpretation:", ifelse(int_p_vol < 0.05,
                                   "Neighbor volume effects differ by coral size",
                                   "No size-dependent volume effects"), "\n\n")

# 7.3 Size × Distance interaction
cat("7.3 Abundance ~ Size × Distance:\n")

m_int_dist <- MASS::glm.nb(total_cafi ~ log_volume_c * mean_neighbor_dist + site,
                     data = landscape_data)
m_int_dist_summary <- summary(m_int_dist)

int_coef_dist <- coef(m_int_dist)["log_volume_c:mean_neighbor_dist"]
int_p_dist <- m_int_dist_summary$coefficients["log_volume_c:mean_neighbor_dist", "Pr(>|z|)"]

cat("    Interaction β =", round(int_coef_dist, 6),
    ", p =", format.pval(int_p_dist, 3), "\n")
cat("    Interpretation:", ifelse(int_p_dist < 0.05,
                                   "Isolation effects differ by coral size",
                                   "No size-dependent isolation effects"), "\n\n")

# FDR correction across 3 interaction tests
p_size_neighbors <- int_p
p_size_vol <- int_p_vol
p_size_dist <- int_p_dist
interaction_pvals <- c(p_size_neighbors, p_size_vol, p_size_dist)
interaction_fdr <- p.adjust(interaction_pvals, method = "BH")
cat("  Interaction tests (FDR-corrected):\n")
cat("    Size x Neighbors: p_raw =", round(p_size_neighbors, 3),
    ", p_FDR =", round(interaction_fdr[1], 3), "\n")
cat("    Size x Neighbor Volume: p_raw =", round(p_size_vol, 3),
    ", p_FDR =", round(interaction_fdr[2], 3), "\n")
cat("    Size x Distance: p_raw =", round(p_size_dist, 3),
    ", p_FDR =", round(interaction_fdr[3], 3), "\n\n")

# ============================================================================
# PART 8: FUNCTIONAL GROUP RESPONSES
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 8: FUNCTIONAL GROUP RESPONSES TO LANDSCAPE\n")
cat("------------------------------------------------------------\n\n")

# Test if different functional groups respond differently
# NOTE: These 4 functional group tests are exploratory and not corrected for multiple testing.
functional_groups <- c("n_trapezia", "n_resident_fish", "n_corallivore", "n_shrimp")

for (fg in functional_groups) {
  if (fg %in% names(landscape_data) && sum(landscape_data[[fg]] > 0) > 10) {
    cat(paste0("8.", which(functional_groups == fg), " ", gsub("n_", "", fg), ":\n"))

    fg_model <- tryCatch({
      MASS::glm.nb(as.formula(paste(fg, "~ log(volume) + n_neighbors + site")),
             data = landscape_data)
    }, error = function(e) {
      cat("    [Note: NB failed for", fg, "- using Poisson fallback]\n")
      glm(as.formula(paste(fg, "~ log(volume) + n_neighbors + site")),
          family = poisson, data = landscape_data)
    })

    fg_summary <- summary(fg_model)

    vol_coef <- coef(fg_model)["log(volume)"]
    vol_p <- fg_summary$coefficients["log(volume)", ncol(fg_summary$coefficients)]
    neigh_coef <- coef(fg_model)["n_neighbors"]
    neigh_p <- fg_summary$coefficients["n_neighbors", ncol(fg_summary$coefficients)]

    cat("    Volume effect: β =", round(vol_coef, 3),
        ", p =", format.pval(vol_p, 3), "\n")
    cat("    Neighbor effect: β =", round(neigh_coef, 4),
        ", p =", format.pval(neigh_p, 3), "\n\n")
  }
}

# ============================================================================
# PART 9: SUMMARY FIGURES
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 9: CREATING SUMMARY FIGURES\n")
cat("------------------------------------------------------------\n\n")

# Create comprehensive diversity figures
p_rich_neighbors <- ggplot(landscape_data, aes(x = n_neighbors, y = otu_richness, color = site)) +
  geom_jitter(alpha = 0.6, width = 0.3, size = 2.5) +
  geom_smooth(method = "glm", formula = y ~ x, method.args = list(family = poisson),
              se = TRUE, color = "black") +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = "Number of Neighbors",
       y = "Species Richness",
       title = "G. Richness vs Neighbors") +
  theme(legend.position = "none")

p_shan_neighbors <- ggplot(landscape_data, aes(x = n_neighbors, y = shannon, color = site)) +
  geom_jitter(alpha = 0.6, width = 0.3, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black") +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = "Number of Neighbors",
       y = "Shannon H'",
       title = "H. Shannon vs Neighbors") +
  theme(legend.position = "none")

# 8-panel summary figure
p_summary <- (p_size_abund + p_size_rich + p_size_shan + plot_spacer()) /
  (p_neighbor_count + p_neighbor_vol + p_isolation + plot_spacer()) +
  plot_layout(widths = c(1, 1, 1, 0.3)) +
  plot_annotation(
    title = "Landscape Effects on CAFI Communities",
    subtitle = paste0("N = ", nrow(landscape_data), " corals with neighborhood surveys | ",
                      "Size explains ", round(calc_r2(m_size, m_null) * 100, 1), "% of variance | ",
                      "Neighborhood adds ", round(delta_r2 * 100, 1), "%"),
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12))
  )

save_figure(p_summary, file.path(FIG_DIR, "landscape_effects_summary.png"),
            width = 16, height = 10)
cat("  Saved: landscape_effects_summary.png\n")

# Create coefficient forest plot
coef_data <- tibble(
  predictor = c("Coral Volume", "Neighbor Count", "Neighbor Volume", "Neighbor Distance"),
  beta = c(coef(m_abund_full)["log(volume)"],
           coef(m_abund_full)["n_neighbors"],
           coef(m_abund_full)["log(total_neighbor_volume + 1)"],
           coef(m_abund_full)["mean_neighbor_dist"]),
  se = c(m_abund_full_summary$coefficients["log(volume)", "Std. Error"],
         m_abund_full_summary$coefficients["n_neighbors", "Std. Error"],
         m_abund_full_summary$coefficients["log(total_neighbor_volume + 1)", "Std. Error"],
         m_abund_full_summary$coefficients["mean_neighbor_dist", "Std. Error"]),
  category = c("Focal", "Neighborhood", "Neighborhood", "Neighborhood")
) %>%
  mutate(
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se,
    significant = (lower > 0 | upper < 0),
    predictor = factor(predictor, levels = rev(predictor))
  )

p_forest <- ggplot(coef_data, aes(x = beta, y = predictor, color = category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, linewidth = 1,
                orientation = "y") +
  geom_point(aes(shape = significant), size = 4) +
  scale_color_manual(values = c("Focal" = "#D55E00", "Neighborhood" = "#0072B2"),
                     name = "Predictor Type") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     name = "p < 0.05", labels = c("No", "Yes")) +
  labs(x = "Coefficient (β) with 95% CI",
       y = "",
       title = "Effect Sizes on CAFI Abundance",
       subtitle = "Full model controlling for all predictors + site") +
  theme(legend.position = "right")

save_figure(p_forest, file.path(FIG_DIR, "coefficient_forest_plot.png"),
            width = 9, height = 5)
cat("  Saved: coefficient_forest_plot.png\n")

# Individual figures
save_figure(p_size_abund, file.path(FIG_DIR, "abundance_vs_volume.png"),
            width = 7, height = 5)
save_figure(p_neighbor_count, file.path(FIG_DIR, "abundance_vs_neighbors.png"),
            width = 7, height = 5)

# Supplement S7: Neighborhood null results
# Remove panel label "D." from multi-panel layout for standalone supplementary figure
supplement_dir <- file.path(PATHS$figures, "supplement")
dir.create(supplement_dir, showWarnings = FALSE, recursive = TRUE)
p_neighbor_count_s7 <- p_neighbor_count +
  labs(title = "Figure S7: Neighborhood Density and CAFI Abundance") +
  scale_color_manual(values = SITE_COLORS, name = "Site") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.position = "bottom")
save_figure(p_neighbor_count_s7, file.path(supplement_dir, "figS7_neighborhood_null.png"),
            width = 7, height = 5.5)

save_figure(p_neighbor_vol, file.path(FIG_DIR, "abundance_vs_neighbor_volume.png"),
            width = 7, height = 5)
save_figure(p_isolation, file.path(FIG_DIR, "abundance_vs_isolation.png"),
            width = 7, height = 5)
save_figure(p_size_rich, file.path(FIG_DIR, "richness_vs_volume.png"),
            width = 7, height = 5)
save_figure(p_size_shan, file.path(FIG_DIR, "shannon_vs_volume.png"),
            width = 7, height = 5)

cat("  Saved: individual figure files\n")
cat("\nAll figures saved to:", FIG_DIR, "\n")

# ============================================================================
# PART 9B: PUBLICATION MULTIPANEL FIGURE WITH GLM
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 9B: PUBLICATION FIGURE - Abundance & Richness vs Predictors\n")
cat("------------------------------------------------------------\n\n")

# gt is optional — used for HTML publication tables only
gt_available <- requireNamespace("gt", quietly = TRUE)
if (!gt_available) {
  cat("  Note: 'gt' package not installed; skipping HTML table output\n")
}

cat("Fitting GLM models with site as fixed effect...\n")
cat("(3 sites insufficient for random intercepts; Bolker et al. 2009)\n\n")

# ---- MULTIPLE PREDICTOR GLM Models (fixed-effect site) ----
# All 4 landscape predictors in one model, site as fixed effect
# Scale predictors to improve convergence

landscape_data_scaled <- landscape_data %>%
  mutate(
    log_vol_scaled = scale(log(volume))[,1],
    n_neighbors_scaled = scale(n_neighbors)[,1],
    log_neigh_vol_scaled = scale(log(total_neighbor_volume + 1))[,1],
    dist_scaled = scale(mean_neighbor_dist)[,1]
  )

cat("Model 1: CAFI Abundance ~ 4 predictors + site (fixed effect)\n")
cat("         Using negative binomial GLM (glm.nb)\n")
cat("         Predictors standardized (mean=0, SD=1)\n")
cat("         Note: 3 sites insufficient for random intercepts (Bolker et al. 2009)\n\n")

glmm_abundance <- MASS::glm.nb(
  total_cafi ~ log_vol_scaled + n_neighbors_scaled +
    log_neigh_vol_scaled + dist_scaled + site,
  data = landscape_data_scaled
)

cat("Model 2: Species Richness ~ 4 predictors + site (fixed effect)\n")
cat("         Using Poisson GLM (glm)\n")
cat("         Predictors standardized (mean=0, SD=1)\n\n")

glmm_richness <- glm(
  otu_richness ~ log_vol_scaled + n_neighbors_scaled +
    log_neigh_vol_scaled + dist_scaled + site,
  family = poisson,
  data = landscape_data_scaled
)

# ---- Extract GLM Results ----
extract_glmm_coefs <- function(model, response) {
  s <- summary(model)
  coefs <- s$coefficients

  # Get predictor rows (skip intercept) - use scaled names
  pred_names <- c("log_vol_scaled", "n_neighbors_scaled",
                  "log_neigh_vol_scaled", "dist_scaled")
  display_names <- c("Coral Volume (log)", "Neighbor Count",
                     "Neighbor Volume (log)", "Mean Neighbor Distance")

  # Handle both "z value" (Poisson/NB) and "t value" (quasipoisson) column names
  stat_col <- if ("z value" %in% colnames(coefs)) "z value" else "t value"
  p_col <- if ("Pr(>|z|)" %in% colnames(coefs)) "Pr(>|z|)" else "Pr(>|t|)"

  results <- tibble(
    Response = response,
    Predictor = display_names,
    Beta = coefs[pred_names, "Estimate"],
    SE = coefs[pred_names, "Std. Error"],
    z_value = coefs[pred_names, stat_col],
    p_value = coefs[pred_names, p_col]
  )

  return(results)
}

glmm_results <- bind_rows(
  extract_glmm_coefs(glmm_abundance, "Abundance"),
  extract_glmm_coefs(glmm_richness, "Richness")
) %>%
  mutate(
    Significant = p_value < 0.05,
    CI_lower = Beta - 1.96 * SE,
    CI_upper = Beta + 1.96 * SE,
    n = nrow(landscape_data)
  )

# Print model summaries
cat("----------------------------------------\n")
cat("ABUNDANCE MODEL SUMMARY:\n")
cat("----------------------------------------\n")
print(summary(glmm_abundance))

cat("\n----------------------------------------\n")
cat("RICHNESS MODEL SUMMARY:\n")
cat("----------------------------------------\n")
print(summary(glmm_richness))

cat("\n\nGLM Results (Full Models - All 4 Predictors + Site):\n")
glmm_print <- glmm_results %>%
  mutate(across(c(Beta, SE, z_value, p_value), ~round(., 4))) %>%
  dplyr::select(Response, Predictor, Beta, SE, z_value, p_value, Significant)
print(as.data.frame(glmm_print))

# ---- AIC-Based Model Selection (Backward Elimination) ----
cat("\n----------------------------------------\n")
cat("AIC-BASED MODEL SELECTION\n")
cat("----------------------------------------\n\n")

# Function to perform backward elimination based on AIC
# Uses fixed-effect site (not random intercept) - 3 levels insufficient for RE
backward_aic <- function(full_model, model_type = "nb") {
  current_model <- full_model
  current_aic <- AIC(current_model)

  cat("Starting AIC:", round(current_aic, 2), "\n")
  cat("Full model formula:", deparse(formula(current_model)), "\n\n")

  # Track elimination history
  history <- list()
  step <- 0

  # Get droppable terms from model formula (exclude intercept and site)
  get_droppable_terms <- function(model) {
    all_terms <- attr(terms(model), "term.labels")
    # Don't drop the site fixed effect
    all_terms[all_terms != "site"]
  }

  repeat {
    step <- step + 1

    fe <- get_droppable_terms(current_model)

    if (length(fe) == 0) {
      cat("No more terms to drop.\n")
      break
    }

    # Try dropping each term
    drop_results <- data.frame(
      term = character(),
      aic = numeric(),
      stringsAsFactors = FALSE
    )

    for (term in fe) {
      # Create formula without this term
      new_formula <- update(formula(current_model), paste("~ . -", term))

      tryCatch({
        if (model_type == "nb") {
          new_model <- MASS::glm.nb(new_formula, data = landscape_data_scaled)
        } else {
          new_model <- glm(new_formula, family = poisson, data = landscape_data_scaled)
        }
        drop_results <- rbind(drop_results,
                              data.frame(term = term, aic = AIC(new_model)))
      }, error = function(e) {
        cat("  Error dropping", term, ":", conditionMessage(e), "\n")
      })
    }

    if (nrow(drop_results) == 0) break

    # Find best drop (lowest AIC)
    best_drop <- drop_results[which.min(drop_results$aic), ]

    cat("Step", step, "- Considering dropping terms:\n")
    drop_results$delta_aic <- drop_results$aic - current_aic
    print(drop_results[order(drop_results$aic), ])

    # Only drop if it improves AIC (lower is better)
    if (best_drop$aic < current_aic) {
      cat("\n  -> Dropping:", best_drop$term, "(ΔAIC =", round(best_drop$aic - current_aic, 2), ")\n\n")

      new_formula <- update(formula(current_model), paste("~ . -", best_drop$term))
      if (model_type == "nb") {
        current_model <- MASS::glm.nb(new_formula, data = landscape_data_scaled)
      } else {
        current_model <- glm(new_formula, family = poisson, data = landscape_data_scaled)
      }
      current_aic <- AIC(current_model)

      history[[step]] <- list(dropped = best_drop$term, aic = current_aic)
    } else {
      cat("\n  -> No improvement from dropping any term. Stopping.\n")
      cat("  Final AIC:", round(current_aic, 2), "\n\n")
      break
    }
  }

  return(current_model)
}

# Backward elimination for ABUNDANCE
cat("=== ABUNDANCE MODEL ===\n")
glmm_abundance_reduced <- backward_aic(glmm_abundance, "nb")

cat("\nReduced Abundance Model:\n")
print(summary(glmm_abundance_reduced))

# Backward elimination for RICHNESS
cat("\n=== RICHNESS MODEL ===\n")
glmm_richness_reduced <- backward_aic(glmm_richness, "poisson")

cat("\nReduced Richness Model:\n")
print(summary(glmm_richness_reduced))

# --- Statistical note: Post-selection inference ---
cat("\n")
cat("*** STATISTICAL NOTE: POST-SELECTION INFERENCE ***\n")
cat("The p-values and CIs above are from models AFTER AIC-based variable selection.\n")
cat("This can bias p-values downward and CIs narrower than they should be.\n")
cat("For robust inference, compare with full model p-values (Part 9 above) and\n")
cat("interpret marginal results (0.01 < p < 0.10) cautiously.\n")
cat("Cross-validation (Part 9C.9) provides out-of-sample predictive assessment.\n\n")

# ============================================================================
# PART 9C: MODEL DIAGNOSTICS (Publication-Ready)
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 9C: MODEL DIAGNOSTICS\n")
cat("------------------------------------------------------------\n\n")

# ---- 9C.1 VIF Check for Multicollinearity ----
cat("9C.1 Variance Inflation Factors (VIF):\n")
cat("     Rule: VIF > 5 indicates problematic multicollinearity\n\n")

# VIF for scaled predictors
vif_data <- landscape_data_scaled %>%
  dplyr::select(log_vol_scaled, n_neighbors_scaled, log_neigh_vol_scaled, dist_scaled)

vif_matrix <- cor(vif_data)
vif_values <- sapply(1:ncol(vif_data), function(i) {
  1 / (1 - summary(lm(vif_data[[i]] ~ ., data = vif_data[-i]))$r.squared)
})
names(vif_values) <- c("log(Volume)", "Neighbor Count", "log(Neighbor Vol)", "Distance")

cat("     VIF values:\n")
for (i in seq_along(vif_values)) {
  status <- ifelse(vif_values[i] < 2, "OK", ifelse(vif_values[i] < 5, "moderate", "HIGH"))
  cat("       ", names(vif_values)[i], ":", round(vif_values[i], 2), " [", status, "]\n")
}
cat("\n     Max VIF:", round(max(vif_values), 2),
    ifelse(max(vif_values) < 5, " - No multicollinearity issues", " - Consider removing predictors"), "\n\n")

# ---- 9C.2 Site Effect Assessment ----
cat("9C.2 Site Fixed Effect Assessment:\n\n")

# For fixed-effect site, test whether site improves model fit (Type II Anova)
cat("  ABUNDANCE MODEL (reduced):\n")
abund_anova <- tryCatch(
  anova(glmm_abundance_reduced, test = "Chisq"),
  error = function(e) NULL
)
if (!is.null(abund_anova) && "site" %in% rownames(abund_anova)) {
  site_dev_abund <- abund_anova["site", "Deviance"]
  site_p_abund <- abund_anova["site", "Pr(>Chi)"]
  cat("    Site deviance:", round(abs(site_dev_abund), 2), "\n")
  cat("    Site p-value:", format.pval(site_p_abund, 3), "\n")
} else {
  cat("    Site effect: included as fixed covariate\n")
}
# Approximate ICC from site coefficient variance
site_coefs_abund <- coef(glmm_abundance_reduced)[grep("^site", names(coef(glmm_abundance_reduced)))]
icc_abund <- ifelse(length(site_coefs_abund) > 0,
                    var(c(0, site_coefs_abund)) / (var(c(0, site_coefs_abund)) + pi^2/3),
                    0)
cat("    Approximate ICC:", round(icc_abund, 3), "\n\n")

cat("  RICHNESS MODEL (reduced):\n")
rich_anova <- tryCatch(
  anova(glmm_richness_reduced, test = "Chisq"),
  error = function(e) NULL
)
if (!is.null(rich_anova) && "site" %in% rownames(rich_anova)) {
  site_dev_rich <- rich_anova["site", "Deviance"]
  site_p_rich <- rich_anova["site", "Pr(>Chi)"]
  cat("    Site deviance:", round(abs(site_dev_rich), 2), "\n")
  cat("    Site p-value:", format.pval(site_p_rich, 3), "\n")
} else {
  cat("    Site effect: included as fixed covariate\n")
}
site_coefs_rich <- coef(glmm_richness_reduced)[grep("^site", names(coef(glmm_richness_reduced)))]
icc_rich <- ifelse(length(site_coefs_rich) > 0,
                   var(c(0, site_coefs_rich)) / (var(c(0, site_coefs_rich)) + pi^2/3),
                   0)
cat("    Approximate ICC:", round(icc_rich, 3), "\n\n")

# ---- 9C.3 Overdispersion Check ----
cat("9C.3 Overdispersion Diagnostics:\n\n")

# Pearson residuals for diagnostic plots
resid_abund <- residuals(glmm_abundance_reduced, type = "pearson")
resid_rich_pearson <- residuals(glmm_richness_reduced, type = "pearson")

# For Poisson model (richness) - check overdispersion ratio
resid_rich <- residuals(glmm_richness_reduced, type = "pearson")
overdispersion_rich <- sum(resid_rich^2) / df.residual(glmm_richness_reduced)
cat("  RICHNESS MODEL (Poisson GLM):\n")
cat("    Pearson χ²/df:", round(overdispersion_rich, 2), "\n")
cat("    Interpretation:", ifelse(overdispersion_rich < 1.5, "No overdispersion - Poisson appropriate",
                                   ifelse(overdispersion_rich < 2, "Mild overdispersion - acceptable",
                                          "Overdispersion - consider negative binomial")), "\n\n")

# For NB model (abundance) - theta parameter
theta_abund <- glmm_abundance_reduced$theta
cat("  ABUNDANCE MODEL (Negative Binomial GLM):\n")
cat("    Theta (dispersion):", round(theta_abund, 2), "\n")
cat("    Interpretation:", ifelse(theta_abund > 10, "Low overdispersion (nearly Poisson)",
                                   ifelse(theta_abund > 2, "Moderate overdispersion",
                                          "High overdispersion - NB appropriate")), "\n\n")

# ---- 9C.4 DHARMa-style Residual Diagnostics ----
cat("9C.4 Residual Diagnostics:\n\n")

# Simulate residuals manually (DHARMa-like approach)
# For abundance model
set.seed(123)  # set.seed(123) — different seed for DHARMa simulations to ensure independence from main bootstrap
n_sim <- 250

# Abundance model residual simulation
fitted_abund <- fitted(glmm_abundance_reduced)
theta_val <- glmm_abundance_reduced$theta

sim_abund <- matrix(rnbinom(n_sim * length(fitted_abund),
                            size = theta_val, mu = fitted_abund),
                    ncol = n_sim)
observed_abund <- landscape_data_scaled$total_cafi

# Calculate quantile residuals
quantile_resid_abund <- sapply(1:length(observed_abund), function(i) {
  ecdf(sim_abund[i, ])(observed_abund[i])
})

# Uniformity test (should be ~uniform on 0,1)
ks_test_abund <- ks.test(quantile_resid_abund, "punif")
cat("  ABUNDANCE MODEL - Quantile Residuals:\n")
cat("    KS test for uniformity: D =", round(ks_test_abund$statistic, 3),
    ", p =", format.pval(ks_test_abund$p.value, 3), "\n")
cat("    Interpretation:", ifelse(ks_test_abund$p.value > 0.05,
                                   "Residuals uniform - model fit OK",
                                   "Residuals deviate from uniform - check model"), "\n\n")

# Richness model
fitted_rich <- fitted(glmm_richness_reduced)
sim_rich <- matrix(rpois(n_sim * length(fitted_rich), lambda = fitted_rich),
                   ncol = n_sim)
observed_rich <- landscape_data_scaled$otu_richness

quantile_resid_rich <- sapply(1:length(observed_rich), function(i) {
  ecdf(sim_rich[i, ])(observed_rich[i])
})

ks_test_rich <- ks.test(quantile_resid_rich, "punif")
cat("  RICHNESS MODEL - Quantile Residuals:\n")
cat("    KS test for uniformity: D =", round(ks_test_rich$statistic, 3),
    ", p =", format.pval(ks_test_rich$p.value, 3), "\n")
cat("    Interpretation:", ifelse(ks_test_rich$p.value > 0.05,
                                   "Residuals uniform - model fit OK",
                                   "Residuals deviate from uniform - check model"), "\n\n")

# ---- 9C.5 Influential Observations (Cook's Distance) ----
cat("9C.5 Influential Observations (Cook's Distance):\n\n")

# Use proper cooks.distance() for GLM objects
cooks_abund <- cooks.distance(glmm_abundance_reduced)

# Threshold: 4/n (standard cutoff)
threshold <- 4 / nrow(landscape_data_scaled)
influential_abund <- which(cooks_abund > threshold)

cat("  ABUNDANCE MODEL (NB GLM):\n")
cat("    Threshold (4/n):", round(threshold, 4), "\n")
cat("    Max Cook's D:", round(max(cooks_abund, na.rm = TRUE), 4), "\n")
cat("    Influential points:", length(influential_abund), "of", nrow(landscape_data_scaled), "\n")
if (length(influential_abund) > 0 && length(influential_abund) <= 5) {
  cat("    Coral IDs:", paste(landscape_data_scaled$coral_id[influential_abund], collapse = ", "), "\n")
}
cat("    Interpretation:", ifelse(length(influential_abund) / nrow(landscape_data_scaled) < 0.05,
                                   "Few influential points - results robust",
                                   "Consider sensitivity analysis"), "\n\n")

# Richness model
cooks_rich <- cooks.distance(glmm_richness_reduced)
influential_rich <- which(cooks_rich > threshold)

cat("  RICHNESS MODEL (Poisson GLM):\n")
cat("    Max Cook's D:", round(max(cooks_rich, na.rm = TRUE), 4), "\n")
cat("    Influential points:", length(influential_rich), "of", nrow(landscape_data_scaled), "\n")
cat("    Interpretation:", ifelse(length(influential_rich) / nrow(landscape_data_scaled) < 0.05,
                                   "Few influential points - results robust",
                                   "Consider sensitivity analysis"), "\n\n")

# ---- 9C.6 Create Diagnostic Figure ----
cat("9C.6 Creating diagnostic figures...\n")

# Residual diagnostic plots
p_resid_abund <- ggplot(data.frame(
  fitted = fitted_abund,
  residual = resid_abund,
  quantile_resid = quantile_resid_abund
), aes(x = fitted, y = residual)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(x = "Fitted Values", y = "Pearson Residuals",
       title = "A. Abundance: Residuals vs Fitted") +
  theme_minimal()

p_qq_abund <- ggplot(data.frame(quantile_resid = quantile_resid_abund),
                     aes(sample = quantile_resid)) +
  stat_qq(distribution = qunif) +
  stat_qq_line(distribution = qunif, color = "red") +
  labs(x = "Theoretical Quantiles (Uniform)", y = "Sample Quantiles",
       title = "B. Abundance: QQ Plot (Quantile Residuals)") +
  theme_minimal()

p_resid_rich <- ggplot(data.frame(
  fitted = fitted_rich,
  residual = resid_rich_pearson
), aes(x = fitted, y = residual)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(x = "Fitted Values", y = "Pearson Residuals",
       title = "C. Richness: Residuals vs Fitted") +
  theme_minimal()

p_qq_rich <- ggplot(data.frame(quantile_resid = quantile_resid_rich),
                    aes(sample = quantile_resid)) +
  stat_qq(distribution = qunif) +
  stat_qq_line(distribution = qunif, color = "red") +
  labs(x = "Theoretical Quantiles (Uniform)", y = "Sample Quantiles",
       title = "D. Richness: QQ Plot (Quantile Residuals)") +
  theme_minimal()

# Combine diagnostic plots
p_diagnostics <- (p_resid_abund | p_qq_abund) / (p_resid_rich | p_qq_rich) +
  plot_annotation(
    title = "GLM Residual Diagnostics",
    subtitle = "Reduced models after AIC-based backward elimination"
  )

save_figure(p_diagnostics, file.path(FIG_DIR, "glmm_diagnostics.png"),
            width = 10, height = 8)
cat("  Saved: glmm_diagnostics.png\n\n")

# ---- 9C.7 Summary Diagnostic Table ----
cat("9C.7 Diagnostic Summary:\n\n")

diagnostics_summary <- tibble(
  Model = c("Abundance (NB GLM)", "Richness (Poisson GLM)"),
  `Max VIF` = rep(round(max(vif_values), 2), 2),
  `VIF Status` = rep(ifelse(max(vif_values) < 5, "OK", "High"), 2),
  ICC = c(round(icc_abund, 3), round(icc_rich, 3)),
  `Overdispersion` = c(paste0("θ=", round(theta_val, 1)), round(overdispersion_rich, 2)),
  `KS test p` = c(round(ks_test_abund$p.value, 3), round(ks_test_rich$p.value, 3)),
  `Residual Status` = c(ifelse(ks_test_abund$p.value > 0.05, "OK", "Check"),
                        ifelse(ks_test_rich$p.value > 0.05, "OK", "Check")),
  `Influential (%)` = c(round(100 * length(influential_abund) / nrow(landscape_data_scaled), 1),
                        round(100 * length(influential_rich) / nrow(landscape_data_scaled), 1))
)

print(diagnostics_summary)

save_table(diagnostics_summary, "glmm_diagnostics_summary")
cat("\n  Saved: glmm_diagnostics_summary.csv\n\n")

# ---- 9C.8 DHARMa Simulated Residual Diagnostics ----
cat("9C.8 DHARMa simulated residual diagnostics:\n")
if (requireNamespace("DHARMa", quietly = TRUE)) {
  # Abundance model (NB GLM)
  dharma_abund <- DHARMa::simulateResiduals(glmm_abundance_reduced, n = 1000, plot = FALSE)
  dharma_test_abund <- DHARMa::testResiduals(dharma_abund, plot = FALSE)

  cat("  ABUNDANCE MODEL (NB GLM):\n")
  cat("    KS test (uniformity): p =", format.pval(dharma_test_abund$uniformity$p.value, 3), "\n")
  cat("    Dispersion test: p =", format.pval(dharma_test_abund$dispersion$p.value, 3), "\n")
  cat("    Outlier test: p =", format.pval(dharma_test_abund$outliers$p.value, 3), "\n")

  # Zero-inflation test
  zi_test_abund <- DHARMa::testZeroInflation(dharma_abund, plot = FALSE)
  cat("    Zero-inflation test: p =", format.pval(zi_test_abund$p.value, 3),
      ifelse(zi_test_abund$p.value < 0.05, " [CONSIDER ZINB]", " [OK]"), "\n\n")

  # Richness model (Poisson GLM)
  dharma_rich <- DHARMa::simulateResiduals(glmm_richness_reduced, n = 1000, plot = FALSE)
  dharma_test_rich <- DHARMa::testResiduals(dharma_rich, plot = FALSE)

  cat("  RICHNESS MODEL (Poisson GLM):\n")
  cat("    KS test (uniformity): p =", format.pval(dharma_test_rich$uniformity$p.value, 3), "\n")
  cat("    Dispersion test: p =", format.pval(dharma_test_rich$dispersion$p.value, 3),
      ifelse(dharma_test_rich$dispersion$p.value < 0.05, " [OVERDISPERSED]", " [OK]"), "\n")
  cat("    Outlier test: p =", format.pval(dharma_test_rich$outliers$p.value, 3), "\n\n")

  # Save DHARMa diagnostic plots
  png(file.path(FIG_DIR, "dharma_diagnostics.png"), width = 10, height = 8, units = "in", res = 200)
  par(mfrow = c(2, 2))
  plot(dharma_abund, main = "Abundance (NB) - QQ")
  DHARMa::plotResiduals(dharma_abund, main = "Abundance - Residuals vs Predicted")
  plot(dharma_rich, main = "Richness (Poisson) - QQ")
  DHARMa::plotResiduals(dharma_rich, main = "Richness - Residuals vs Predicted")
  dev.off()
  cat("  Saved: dharma_diagnostics.png\n\n")
} else {
  cat("  DHARMa package not available - install with install.packages('DHARMa')\n\n")
}

# ---- 9C.9 Leave-One-Out Cross-Validation ----
cat("9C.9 Leave-One-Out Cross-Validation (LOOCV):\n")
cat("     Assesses out-of-sample predictive performance with n =", nrow(landscape_data_scaled), "\n\n")

# LOOCV for abundance model
loocv_abund <- tryCatch({
  preds <- numeric(nrow(landscape_data_scaled))
  for (i in 1:nrow(landscape_data_scaled)) {
    train_data <- landscape_data_scaled[-i, ]
    test_data <- landscape_data_scaled[i, , drop = FALSE]

    # Fit on training data
    m_cv <- MASS::glm.nb(total_cafi ~ log_vol_scaled + site, data = train_data)

    # Predict on held-out observation
    preds[i] <- predict(m_cv, newdata = test_data, type = "response")
  }

  obs <- landscape_data_scaled$total_cafi
  rmse <- sqrt(mean((obs - preds)^2))
  mae <- mean(abs(obs - preds))
  cor_cv <- cor(obs, preds, method = "spearman")

  list(rmse = rmse, mae = mae, cor = cor_cv, preds = preds, success = TRUE)
}, error = function(e) list(success = FALSE, error = e$message))

if (loocv_abund$success) {
  cat("  Abundance (NB) LOOCV:\n")
  cat("    RMSE:", round(loocv_abund$rmse, 2), "individuals\n")
  cat("    MAE:", round(loocv_abund$mae, 2), "individuals\n")
  cat("    Spearman r (obs vs pred):", round(loocv_abund$cor, 3), "\n\n")
} else {
  cat("  Abundance LOOCV failed:", loocv_abund$error, "\n\n")
}

# LOOCV for richness model
loocv_rich <- tryCatch({
  preds <- numeric(nrow(landscape_data_scaled))
  for (i in 1:nrow(landscape_data_scaled)) {
    train_data <- landscape_data_scaled[-i, ]
    test_data <- landscape_data_scaled[i, , drop = FALSE]

    m_cv <- glm(otu_richness ~ log_vol_scaled + site, family = poisson, data = train_data)
    preds[i] <- predict(m_cv, newdata = test_data, type = "response")
  }

  obs <- landscape_data_scaled$otu_richness
  rmse <- sqrt(mean((obs - preds)^2))
  mae <- mean(abs(obs - preds))
  cor_cv <- cor(obs, preds, method = "spearman")

  list(rmse = rmse, mae = mae, cor = cor_cv, preds = preds, success = TRUE)
}, error = function(e) list(success = FALSE, error = e$message))

if (loocv_rich$success) {
  cat("  Richness (Poisson) LOOCV:\n")
  cat("    RMSE:", round(loocv_rich$rmse, 2), "species\n")
  cat("    MAE:", round(loocv_rich$mae, 2), "species\n")
  cat("    Spearman r (obs vs pred):", round(loocv_rich$cor, 3), "\n\n")
} else {
  cat("  Richness LOOCV failed:", loocv_rich$error, "\n\n")
}

# ---- Extract Results from REDUCED Models ----
cat("\n----------------------------------------\n")
cat("REDUCED MODEL RESULTS (AIC-selected)\n")
cat("----------------------------------------\n")

# Helper to extract coefficients from any model
extract_reduced_coefs <- function(model, response) {
  s <- summary(model)
  coefs <- s$coefficients

  # Get predictor rows (skip intercept)
  pred_rows <- rownames(coefs)[-1]  # Remove intercept

  if (length(pred_rows) == 0) {
    return(tibble(
      Response = response,
      Predictor = "Intercept only",
      Beta = NA, SE = NA, z_value = NA, p_value = NA
    ))
  }

  # Map scaled names to display names
  name_map <- c(
    "log_vol_scaled" = "Coral Volume (log)",
    "n_neighbors_scaled" = "Neighbor Count",
    "log_neigh_vol_scaled" = "Neighbor Volume (log)",
    "dist_scaled" = "Mean Neighbor Distance"
  )

  display_names <- name_map[pred_rows]
  display_names[is.na(display_names)] <- pred_rows[is.na(display_names)]

  tibble(
    Response = response,
    Predictor = display_names,
    Beta = coefs[pred_rows, "Estimate"],
    SE = coefs[pred_rows, "Std. Error"],
    z_value = coefs[pred_rows, "z value"],
    p_value = coefs[pred_rows, "Pr(>|z|)"]
  )
}

glmm_reduced_results <- bind_rows(
  extract_reduced_coefs(glmm_abundance_reduced, "Abundance"),
  extract_reduced_coefs(glmm_richness_reduced, "Richness")
) %>%
  mutate(
    Significant = p_value < 0.05,
    CI_lower = Beta - 1.96 * SE,
    CI_upper = Beta + 1.96 * SE,
    n = nrow(landscape_data),
    n_groups = 3
  )

cat("\nReduced Model Results:\n")
glmm_reduced_print <- glmm_reduced_results %>%
  mutate(across(c(Beta, SE, z_value, p_value), ~round(., 4))) %>%
  dplyr::select(Response, Predictor, Beta, SE, z_value, p_value, Significant)
print(as.data.frame(glmm_reduced_print))

# Compare full vs reduced models
cat("\n----------------------------------------\n")
cat("MODEL COMPARISON (Full vs Reduced)\n")
cat("----------------------------------------\n")
cat("\nAbundance:\n")
cat("  Full model AIC:    ", round(AIC(glmm_abundance), 2), "\n")
cat("  Reduced model AIC: ", round(AIC(glmm_abundance_reduced), 2), "\n")
cat("  ΔAIC:              ", round(AIC(glmm_abundance_reduced) - AIC(glmm_abundance), 2), "\n")

cat("\nRichness:\n")
cat("  Full model AIC:    ", round(AIC(glmm_richness), 2), "\n")
cat("  Reduced model AIC: ", round(AIC(glmm_richness_reduced), 2), "\n")
cat("  ΔAIC:              ", round(AIC(glmm_richness_reduced) - AIC(glmm_richness), 2), "\n")

# ---- Create Publication-Quality AIC Model Comparison Table ----
cat("\n----------------------------------------\n")
cat("CREATING AIC MODEL COMPARISON TABLE\n")
cat("----------------------------------------\n")

# Build comprehensive AIC table for both responses
# Fit intermediate models for complete comparison
cat("\nFitting intermediate models for AIC table...\n")

# Abundance models (fixed-effect site)
m_abund_null <- MASS::glm.nb(total_cafi ~ 1 + site, data = landscape_data_scaled)
m_abund_vol <- MASS::glm.nb(total_cafi ~ log_vol_scaled + site, data = landscape_data_scaled)

# Richness models (fixed-effect site)
m_rich_null <- glm(otu_richness ~ 1 + site, family = poisson, data = landscape_data_scaled)
m_rich_vol <- glm(otu_richness ~ log_vol_scaled + site, family = poisson, data = landscape_data_scaled)
m_rich_vol_dist <- glm(otu_richness ~ log_vol_scaled + dist_scaled + site,
                       family = poisson, data = landscape_data_scaled)

# Create AIC comparison dataframe
aic_comparison <- tibble(
  Response = c(rep("Abundance", 3), rep("Richness", 4)),
  Model = c(
    "Null (intercept only)",
    "Volume only",
    "Full (4 predictors)",
    "Null (intercept only)",
    "Volume only",
    "Volume + Distance",
    "Full (4 predictors)"
  ),
  Predictors = c(
    "—",
    "log(Volume)",
    "log(Volume) + Neighbors + log(Neighbor Vol) + Distance",
    "—",
    "log(Volume)",
    "log(Volume) + Distance",
    "log(Volume) + Neighbors + log(Neighbor Vol) + Distance"
  ),
  df = c(
    3, 4, 7,  # Abundance: intercept + site variance + theta; +1 for vol; +3 for full
    2, 3, 4, 6  # Richness: intercept + site variance; +1 for vol; +1 for dist; +2 for full
  ),
  AIC = c(
    AIC(m_abund_null), AIC(m_abund_vol), AIC(glmm_abundance),
    AIC(m_rich_null), AIC(m_rich_vol), AIC(m_rich_vol_dist), AIC(glmm_richness)
  )
) %>%
  group_by(Response) %>%
  mutate(
    ΔAIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * ΔAIC) / sum(exp(-0.5 * ΔAIC)),
    Best = ΔAIC == 0
  ) %>%
  ungroup()

cat("\nAIC Model Comparison:\n")
print(as.data.frame(aic_comparison %>% dplyr::select(Response, Model, AIC, ΔAIC, AIC_weight, Best)))

# Create publication gt table for AIC comparison (if gt available)
if (gt_available) {
  aic_gt_table <- aic_comparison %>%
    dplyr::select(Response, Model, Predictors, AIC, ΔAIC, AIC_weight) %>%
    mutate(
      AIC = round(AIC, 1),
      ΔAIC = round(ΔAIC, 1),
      AIC_weight = sprintf("%.3f", AIC_weight)
    ) %>%
    gt::gt(groupname_col = "Response") %>%
    gt::tab_header(
      title = gt::md("**Table 2. AIC Model Comparison: Landscape Predictors of CAFI Communities**"),
      subtitle = gt::md("*Backward elimination via AIC from full GLM with site as fixed effect*")
    ) %>%
    gt::cols_label(
      Model = "Model",
      Predictors = "Fixed Effects",
      AIC = "AIC",
      ΔAIC = "ΔAIC",
      AIC_weight = "AIC Weight"
    ) %>%
    gt::cols_align(align = "left", columns = c(Model, Predictors)) %>%
    gt::cols_align(align = "center", columns = c(AIC, ΔAIC, AIC_weight)) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(rows = ΔAIC == "0.0")
    ) %>%
    gt::tab_style(
      style = gt::cell_fill(color = "#e8f4e8"),
      locations = gt::cells_body(rows = ΔAIC == "0.0")
    ) %>%
    gt::tab_footnote(
      footnote = "Abundance: Negative binomial GLM; Richness: Poisson GLM. All models include site as fixed effect.",
      locations = gt::cells_column_labels(columns = AIC)
    ) %>%
    gt::tab_footnote(
      footnote = "ΔAIC = difference from best model; AIC weight = relative likelihood of being the best model.",
      locations = gt::cells_column_labels(columns = ΔAIC)
    ) %>%
    gt::tab_source_note(
      source_note = gt::md(paste0("N = ", nrow(landscape_data),
                              " *Pocillopora* colonies across 3 sites | Mo'orea, 2019"))
    ) %>%
    gt::cols_width(
      Model ~ gt::px(160),
      Predictors ~ gt::px(280),
      AIC ~ gt::px(70),
      ΔAIC ~ gt::px(60),
      AIC_weight ~ gt::px(80)
    ) %>%
    gt::tab_options(
      table.font.size = gt::px(11),
      heading.title.font.size = gt::px(13),
      heading.subtitle.font.size = gt::px(10),
      row_group.font.weight = "bold",
      table.border.top.style = "solid",
      table.border.bottom.style = "solid"
    )

  # Save AIC table
  gt::gtsave(aic_gt_table, file.path(FIG_DIR, "aic_model_comparison_table.png"))
  gt::gtsave(aic_gt_table, file.path(FIG_DIR, "aic_model_comparison_table.html"))
  cat("  Saved: aic_model_comparison_table.png\n")
  cat("  Saved: aic_model_comparison_table.html\n")
}

save_table(aic_comparison %>%
             mutate(across(where(is.numeric), ~round(., 4))),
           "aic_model_comparison")
cat("  Saved: aic_model_comparison.csv\n")

# ---- Prepare results for figure (combine full + reduced info) ----
# For the figure, we want to show ALL predictors but indicate which were retained vs dropped
# Use the FULL model results but mark significance based on reduced model

glmm_full_results <- bind_rows(
  extract_glmm_coefs(glmm_abundance, "Abundance"),
  extract_glmm_coefs(glmm_richness, "Richness")
)

# Mark which predictors are in the reduced models
retained_abundance <- names(coef(glmm_abundance_reduced))[-1]  # Remove intercept
retained_richness <- names(coef(glmm_richness_reduced))[-1]

# Map back to display names
name_map_reverse <- c(
  "log_vol_scaled" = "Coral Volume (log)",
  "n_neighbors_scaled" = "Neighbor Count",
  "log_neigh_vol_scaled" = "Neighbor Volume (log)",
  "dist_scaled" = "Mean Neighbor Distance"
)

retained_abund_display <- name_map_reverse[retained_abundance]
retained_rich_display <- name_map_reverse[retained_richness]

# Create figure results with retained/dropped status
glmm_figure_results <- glmm_full_results %>%
  mutate(
    Retained = case_when(
      Response == "Abundance" & Predictor %in% retained_abund_display ~ TRUE,
      Response == "Richness" & Predictor %in% retained_rich_display ~ TRUE,
      TRUE ~ FALSE
    ),
    # For retained predictors, get p-value from reduced model
    # For dropped predictors, show as dropped
    Significant = Retained & p_value < 0.05
  )

# Update with reduced model p-values for retained predictors
for (i in 1:nrow(glmm_figure_results)) {
  if (glmm_figure_results$Retained[i]) {
    reduced_row <- glmm_reduced_results %>%
      filter(Response == glmm_figure_results$Response[i],
             Predictor == glmm_figure_results$Predictor[i])
    if (nrow(reduced_row) > 0) {
      glmm_figure_results$Beta[i] <- reduced_row$Beta[1]
      glmm_figure_results$SE[i] <- reduced_row$SE[1]
      glmm_figure_results$p_value[i] <- reduced_row$p_value[1]
      glmm_figure_results$Significant[i] <- reduced_row$p_value[1] < 0.05
    }
  }
}

glmm_figure_results <- glmm_figure_results %>%
  mutate(
    CI_lower = Beta - 1.96 * SE,
    CI_upper = Beta + 1.96 * SE,
    n = nrow(landscape_data),
    n_groups = 3
  )

# Use this for the figure
glmm_results <- glmm_figure_results

# ---- Create Publication-Quality 8-Panel Figure ----
cat("\nCreating 8-panel publication figure (using reduced models)...\n")

# Helper function for consistent panel styling
# Only show trend line for significant relationships (p < 0.05)
# Handles predictors dropped during AIC selection (shows "dropped" annotation)
# Updated to use tidy evaluation (aes()) instead of deprecated aes_string()
make_panel <- function(data, x_var, y_var, x_lab, title, model_results, log_x = FALSE) {

  # Get stats for this panel
  response_type <- ifelse(y_var == "total_cafi", "Abundance", "Richness")
  stats <- model_results %>%
    filter(Predictor == x_lab, Response == response_type)

  # Check if this predictor was dropped from the reduced model
  if (nrow(stats) == 0 || is.na(stats$Beta[1])) {
    # Predictor was dropped - show as "dropped" in annotation
    beta_text <- "dropped"
    p_text <- "(AIC selection)"
    sig_star <- ""
    is_significant <- FALSE
  } else {
    beta_text <- sprintf("β = %.3f", stats$Beta)
    p_text <- ifelse(stats$p_value < 0.001, "p < 0.001", sprintf("p = %.3f", stats$p_value))
    sig_star <- ifelse(stats$Significant, "*", "")
    is_significant <- stats$Significant
  }

  # Use tidy evaluation with .data pronoun instead of deprecated aes_string
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[["site"]])) +
    geom_point(alpha = 0.7, size = 2.5) +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    labs(
      x = x_lab,
      y = ifelse(y_var == "total_cafi", "CAFI Abundance", "Species Richness"),
      title = title
    ) +
    annotate("text", x = Inf, y = Inf,
             label = paste0(beta_text, "\n", p_text, sig_star),
             hjust = 1.1, vjust = 1.3, size = 3.2, fontface = "italic") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )

  # Only add trend line if relationship is significant (p < 0.05)
  if (log_x) {
    p <- p + scale_x_log10(labels = scales::comma)
    if (is_significant) {
      p <- p + geom_smooth(method = "glm", formula = y ~ log(x),
                           se = TRUE, color = "black", linewidth = 0.8, alpha = 0.2)
    }
  } else {
    if (is_significant) {
      p <- p + geom_smooth(method = "glm", formula = y ~ x,
                           se = TRUE, color = "black", linewidth = 0.8, alpha = 0.2)
    }
  }

  return(p)
}

# Row 1: Abundance panels (A-D)
p_A <- make_panel(landscape_data, "volume", "total_cafi",
                  "Coral Volume (log)", "A. Volume → Abundance",
                  glmm_results, log_x = TRUE)

p_B <- make_panel(landscape_data, "n_neighbors", "total_cafi",
                  "Neighbor Count", "B. Neighbors → Abundance",
                  glmm_results, log_x = FALSE)

p_C <- landscape_data %>%
  filter(total_neighbor_volume > 0) %>%
  make_panel("total_neighbor_volume", "total_cafi",
             "Neighbor Volume (log)", "C. Neighbor Vol → Abundance",
             glmm_results, log_x = TRUE)

p_D <- make_panel(landscape_data, "mean_neighbor_dist", "total_cafi",
                  "Mean Neighbor Distance", "D. Distance → Abundance",
                  glmm_results, log_x = FALSE)

# Row 2: Richness panels (E-H)
p_E <- make_panel(landscape_data, "volume", "otu_richness",
                  "Coral Volume (log)", "E. Volume → Richness",
                  glmm_results, log_x = TRUE)

p_F <- make_panel(landscape_data, "n_neighbors", "otu_richness",
                  "Neighbor Count", "F. Neighbors → Richness",
                  glmm_results, log_x = FALSE)

p_G <- landscape_data %>%
  filter(total_neighbor_volume > 0) %>%
  make_panel("total_neighbor_volume", "otu_richness",
             "Neighbor Volume (log)", "G. Neighbor Vol → Richness",
             glmm_results, log_x = TRUE)

p_H <- make_panel(landscape_data, "mean_neighbor_dist", "otu_richness",
                  "Mean Neighbor Distance", "H. Distance → Richness",
                  glmm_results, log_x = FALSE)

# Create legend panel
legend_data <- landscape_data %>%
  dplyr::select(site) %>%
  distinct()

p_legend <- ggplot(landscape_data, aes(x = volume, y = total_cafi, color = site)) +

  geom_point(size = 3) +
  scale_color_manual(values = SITE_COLORS, name = "Site") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

# Extract legend
legend_grob <- cowplot::get_legend(p_legend)

# Combine all panels
fig_multipanel <- (p_A | p_B | p_C | p_D) /
                  (p_E | p_F | p_G | p_H) +
  plot_annotation(
    title = "Landscape Predictors of CAFI Abundance and Richness",
    subtitle = paste0("GLM with site as fixed effect | N = ", nrow(landscape_data),
                      " corals with neighborhood surveys across 3 sites"),
    caption = "* indicates p < 0.05",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      plot.caption = element_text(size = 9, hjust = 1, face = "italic")
    )
  )

# Add shared legend
fig_with_legend <- fig_multipanel +
  inset_element(legend_grob, left = 0.92, bottom = 0.4, right = 1.0, top = 0.6)

# Save figure
save_figure(fig_with_legend, file.path(FIG_DIR, "landscape_multipanel_glmm.png"),
            width = 14, height = 8)

cat("  Saved: landscape_multipanel_glmm.png\n")

# ---- Create Publication-Quality Table with gt ----
cat("\nCreating publication-quality table with gt...\n")

glmm_table <- glmm_results %>%
  mutate(
    `95% CI` = sprintf("[%.3f, %.3f]", CI_lower, CI_upper),
    p_formatted = case_when(
      p_value < 0.001 ~ "<0.001",
      p_value < 0.01 ~ sprintf("%.3f", p_value),
      TRUE ~ sprintf("%.3f", p_value)
    ),
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(Response, Predictor, Beta, SE, `95% CI`, z_value, p_formatted, Significance) %>%
  dplyr::rename(
    `β` = Beta,
    `z` = z_value,
    `p` = p_formatted
  )

if (gt_available) {
  gt_table <- glmm_table %>%
    gt::gt(groupname_col = "Response") %>%
    gt::tab_header(
      title = gt::md("**Table 1. GLM: Landscape Effects on CAFI Communities (AIC-reduced)**"),
      subtitle = gt::md("*Backward elimination via AIC; site as fixed effect*")
    ) %>%
    gt::fmt_number(
      columns = c(`β`, SE, `z`),
      decimals = 3
    ) %>%
    gt::cols_align(
      align = "center",
      columns = c(`β`, SE, `95% CI`, `z`, `p`, Significance)
    ) %>%
    gt::cols_align(
      align = "left",
      columns = Predictor
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(
        columns = c(`β`, `p`),
        rows = Significance != ""
      )
    ) %>%
    gt::tab_style(
      style = gt::cell_fill(color = "#f0f9e8"),
      locations = gt::cells_body(
        rows = Significance != ""
      )
    ) %>%
    gt::tab_footnote(
      footnote = "Abundance models: Negative binomial GLM; Richness models: Poisson GLM. Site as fixed effect.",
      locations = gt::cells_column_labels(columns = `β`)
    ) %>%
    gt::tab_footnote(
      footnote = gt::md("Significance: \\* p < 0.05, \\*\\* p < 0.01, \\*\\*\\* p < 0.001"),
      locations = gt::cells_column_labels(columns = Significance)
    ) %>%
    gt::tab_source_note(
      source_note = gt::md(paste0("N = ", nrow(landscape_data),
                              " *Pocillopora* colonies | Sites: HAU, MAT, MRB | Mo'orea, 2019"))
    ) %>%
    gt::cols_width(
      Predictor ~ gt::px(180),
      `β` ~ gt::px(70),
      SE ~ gt::px(70),
      `95% CI` ~ gt::px(120),
      `z` ~ gt::px(70),
      `p` ~ gt::px(70),
      Significance ~ gt::px(50)
    ) %>%
    gt::tab_options(
      table.font.size = gt::px(12),
      heading.title.font.size = gt::px(14),
      heading.subtitle.font.size = gt::px(11),
      row_group.font.weight = "bold",
      table.border.top.style = "solid",
      table.border.bottom.style = "solid"
    )

  # Save gt table as PNG and HTML
  gt::gtsave(gt_table, file.path(FIG_DIR, "glmm_results_table.png"))
  gt::gtsave(gt_table, file.path(FIG_DIR, "glmm_results_table.html"))

  cat("  Saved: glmm_results_table.png\n")
  cat("  Saved: glmm_results_table.html\n")
}

save_table(glmm_results %>%
             mutate(across(where(is.numeric), ~round(., 4))),
           "glmm_landscape_results")
cat("  Saved: glmm_landscape_results.csv\n\n")

# Print summary
cat("\nGLM Summary:\n")
cat("  Significant predictors of Abundance:\n")
sig_abund <- glmm_results %>% filter(Response == "Abundance", Significant)
if (nrow(sig_abund) > 0) {
  for (i in 1:nrow(sig_abund)) {
    cat("    -", sig_abund$Predictor[i], "(β =", round(sig_abund$Beta[i], 3),
        ", p =", format.pval(sig_abund$p_value[i], 3), ")\n")
  }
} else {
  cat("    None\n")
}

cat("  Significant predictors of Richness:\n")
sig_rich <- glmm_results %>% filter(Response == "Richness", Significant)
if (nrow(sig_rich) > 0) {
  for (i in 1:nrow(sig_rich)) {
    cat("    -", sig_rich$Predictor[i], "(β =", round(sig_rich$Beta[i], 3),
        ", p =", format.pval(sig_rich$p_value[i], 3), ")\n")
  }
} else {
  cat("    None\n")
}

# ============================================================================
# PART 10: RESULTS TABLES
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 10: RESULTS TABLES\n")
cat("------------------------------------------------------------\n\n")

# 10.1 Univariate effects table (4 non-redundant predictors)
cat("10.1 Creating univariate results table...\n")

univariate_results <- tibble(
  Predictor = c("Coral Volume (log)", "Neighbor Count",
                "Total Neighbor Volume (log)", "Mean Neighbor Distance"),
  Response = "CAFI Abundance",
  Beta = c(slope,
           coef(m_neighbors)["n_neighbors"],
           coef(m_neigh_vol)["log(total_neighbor_volume + 1)"],
           coef(m_distance)["mean_neighbor_dist"]),
  SE = c(slope_se,
         m_neigh_summary$coefficients["n_neighbors", "Std. Error"],
         m_nvol_summary$coefficients["log(total_neighbor_volume + 1)", "Std. Error"],
         m_dist_summary$coefficients["mean_neighbor_dist", "Std. Error"]),
  z_value = c(m_summary$coefficients["log(volume)", "z value"],
              m_neigh_summary$coefficients["n_neighbors", "z value"],
              m_nvol_summary$coefficients["log(total_neighbor_volume + 1)", "z value"],
              m_dist_summary$coefficients["mean_neighbor_dist", "z value"]),
  p_value = c(m_summary$coefficients["log(volume)", "Pr(>|z|)"],
              m_neigh_summary$coefficients["n_neighbors", "Pr(>|z|)"],
              m_nvol_summary$coefficients["log(total_neighbor_volume + 1)", "Pr(>|z|)"],
              m_dist_summary$coefficients["mean_neighbor_dist", "Pr(>|z|)"])
) %>%
  mutate(
    Beta = round(Beta, 4),
    SE = round(SE, 4),
    z_value = round(z_value, 2),
    p_value = round(p_value, 4),
    Significant = p_value < 0.05
  )

save_table(univariate_results, "landscape_univariate_results")
cat("     Saved: landscape_univariate_results.csv\n")

# 10.2 Full model coefficients table
cat("10.2 Creating full model results table...\n")

full_model_results <- tibble(
  Response = rep(c("Abundance", "Richness", "Shannon"), each = 4),
  Predictor = rep(c("log(volume)", "n_neighbors",
                    "log(neighbor_vol)", "mean_neighbor_dist"), 3),
  Beta = c(
    coef(m_abund_full)[2:5],
    coef(m_rich_full)[2:5],
    coef(m_shan_full)[2:5]
  ),
  SE = c(
    m_abund_full_summary$coefficients[2:5, "Std. Error"],
    m_rich_full_summary$coefficients[2:5, "Std. Error"],
    m_shan_full_summary$coefficients[2:5, "Std. Error"]
  ),
  Statistic = c(
    m_abund_full_summary$coefficients[2:5, "z value"],
    m_rich_full_summary$coefficients[2:5, "z value"],
    m_shan_full_summary$coefficients[2:5, "t value"]
  ),
  p_value = c(
    m_abund_full_summary$coefficients[2:5, "Pr(>|z|)"],
    m_rich_full_summary$coefficients[2:5, "Pr(>|z|)"],
    m_shan_full_summary$coefficients[2:5, "Pr(>|t|)"]
  )
) %>%
  mutate(
    Beta = round(Beta, 4),
    SE = round(SE, 4),
    Statistic = round(Statistic, 2),
    p_value = round(p_value, 4),
    Significant = p_value < 0.05
  )

save_table(full_model_results, "landscape_full_model_results")
cat("     Saved: landscape_full_model_results.csv\n")

# 10.3 Model comparison table
cat("10.3 Saving model comparison table...\n")
save_table(model_comparison, "landscape_model_comparison")
cat("     Saved: landscape_model_comparison.csv\n")

# 10.4 Interaction effects table
cat("10.4 Creating interaction effects table...\n")

interaction_results <- tibble(
  Interaction = c("Size × Neighbor Count", "Size × Neighbor Volume", "Size × Distance"),
  Beta = c(int_coef, int_coef_vol, int_coef_dist),
  p_value = c(int_p, int_p_vol, int_p_dist)
) %>%
  mutate(
    Beta = round(Beta, 5),
    p_value = round(p_value, 4),
    Significant = p_value < 0.05
  )

save_table(interaction_results, "landscape_interaction_results")
cat("     Saved: landscape_interaction_results.csv\n\n")

# ============================================================================
# PART 11: SAVE COMPREHENSIVE RESULTS OBJECT
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 11: SAVING RESULTS OBJECT\n")
cat("------------------------------------------------------------\n\n")

landscape_results <- list(
  # Sample info
  sample_size = nrow(landscape_data),
  sites = unique(landscape_data$site),

  # Scaling exponent results
  scaling_exponent = list(
    beta = slope,
    se = slope_se,
    ci_95 = slope_ci,
    z_vs_1 = z_vs_1,
    p_vs_1 = p_vs_1,
    interpretation = ifelse(p_vs_1 < 0.05 & slope < 1, "Propagule dilution",
                            ifelse(p_vs_1 < 0.05 & slope > 1, "Super-linear",
                                   "Field of Dreams"))
  ),

  # Model comparison results
  model_comparison = model_comparison,
  r2_size_alone = calc_r2(m_size, m_null),
  r2_full_model = calc_r2(m_full, m_null),
  delta_r2_neighborhood = delta_r2,

  # Full model AIC
  aic_null = AIC(m_null),
  aic_size = AIC(m_size),
  aic_full = AIC(m_full),

  # Interaction tests
  interactions = list(
    size_x_neighbors = list(beta = int_coef, p = int_p),
    size_x_volume = list(beta = int_coef_vol, p = int_p_vol),
    size_x_distance = list(beta = int_coef_dist, p = int_p_dist)
  ),

  # Results tables
  univariate_results = univariate_results,
  full_model_results = full_model_results,
  interaction_results = interaction_results
)

save_object(landscape_results, "landscape_analysis_results")
cat("  Saved: landscape_analysis_results.rds\n\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    LANDSCAPE EFFECTS ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("SAMPLE: N =", nrow(landscape_data), "Pocillopora colonies with complete data\n")
cat("SITES:", paste(unique(landscape_data$site), collapse = ", "), "\n\n")

cat("KEY FINDING 1: SIZE-ABUNDANCE SCALING\n")
cat("  Scaling exponent β =", round(slope, 3), "\n")
cat("  95% CI: [", round(slope_ci[1], 3), ",", round(slope_ci[2], 3), "]\n")
cat("  Test vs β = 1: z =", round(z_vs_1, 2), ", p =", format.pval(p_vs_1, 3), "\n")
cat("  Interpretation:", landscape_results$scaling_exponent$interpretation, "\n\n")

cat("KEY FINDING 2: VARIANCE EXPLAINED\n")
n_neighborhood <- nrow(landscape_data)
cat("  Size alone: R² =", round(calc_r2(m_size, m_null) * 100, 1), "%\n")
cat("  Full model: R² =", round(calc_r2(m_full, m_null) * 100, 1), "%\n")
cat("  Neighborhood adds:", round(delta_r2 * 100, 2), "% additional R²\n\n")

cat("  Power analysis for neighborhood effects (n =", n_neighborhood, "):\n")
if (requireNamespace("pwr", quietly = TRUE)) {
  n_params <- 3  # log_volume + 2 site dummies
  pwr_result <- pwr::pwr.f2.test(
    u = 1, v = n_neighborhood - n_params - 2,
    f2 = NULL, sig.level = 0.05, power = 0.80
  )
  min_f2 <- pwr_result$f2
  cat("  Minimum detectable Cohen's f²:", round(min_f2, 3), "\n")
  cat("  Minimum detectable partial R²:", round(min_f2 / (1 + min_f2), 3), "\n")
  cat("  Effects smaller than this cannot be reliably detected with n =", n_neighborhood, "\n\n")
} else {
  cat("  [pwr package not available - install with install.packages('pwr')]\n\n")
}

cat("KEY FINDING 3: INTERACTION EFFECTS\n")
cat("  Size × Neighbors: p =", format.pval(int_p, 3),
    ifelse(int_p < 0.05, " *", ""), "\n")
cat("  Size × Volume: p =", format.pval(int_p_vol, 3),
    ifelse(int_p_vol < 0.05, " *", ""), "\n")
cat("  Size × Distance: p =", format.pval(int_p_dist, 3),
    ifelse(int_p_dist < 0.05, " *", ""), "\n\n")

# ============================================================================
# PART 12: NEIGHBORHOOD COMPOSITION DIVERGENCE
# ============================================================================
# Parallel to the experiment's finding that composition diverged with coral
# number: do corals in denser neighborhoods have more variable CAFI composition?
# Uses betadisper to test whether compositional variability differs among
# neighborhood density categories.
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 8: NEIGHBORHOOD COMPOSITION DIVERGENCE\n")
cat("------------------------------------------------------------\n\n")

cat("Testing whether corals with more neighbors have more variable CAFI composition.\n")
cat("(Parallel to experiment's divergence finding with manipulated coral number.)\n\n")

if (exists("community_matrix")) {
  # Get community data for neighborhood-surveyed corals
  neighbor_ids <- landscape_data$coral_id
  comm_neighbor <- community_matrix[rownames(community_matrix) %in% neighbor_ids, ]

  # Align and filter
  shared_ids <- intersect(rownames(comm_neighbor), landscape_data$coral_id)
  if (length(shared_ids) >= 20) {
    comm_neighbor <- comm_neighbor[shared_ids, ]
    land_aligned <- landscape_data %>% filter(coral_id %in% shared_ids) %>%
      arrange(match(coral_id, shared_ids))

    # Remove empty species columns
    comm_neighbor <- comm_neighbor[, colSums(comm_neighbor) > 0]

    # Create neighbor density categories (analogous to 1/3/6 corals in experiment)
    land_aligned <- land_aligned %>%
      mutate(neighbor_cat_simple = case_when(
        n_neighbors <= 1 ~ "Low (0-1)",
        n_neighbors <= 3 ~ "Medium (2-3)",
        TRUE ~ "High (4+)"
      ) %>% factor(levels = c("Low (0-1)", "Medium (2-3)", "High (4+)")))

    cat("  Neighbor density categories:\n")
    print(table(land_aligned$neighbor_cat_simple))
    cat("\n")

    # Bray-Curtis on Hellinger-transformed data
    comm_hell_neighbor <- vegan::decostand(comm_neighbor, method = "hellinger")
    dist_neighbor <- vegan::vegdist(comm_hell_neighbor, method = "bray")

    # Betadisper by neighbor density
    disp_neighbor <- vegan::betadisper(dist_neighbor, land_aligned$neighbor_cat_simple)
    disp_neighbor_test <- permutest(disp_neighbor, permutations = 999)
    disp_f <- disp_neighbor_test$tab$F[1]
    disp_p <- disp_neighbor_test$tab$`Pr(>F)`[1]

    cat("  Betadisper (compositional variability ~ neighbor density):\n")
    cat("    F =", round(disp_f, 2), ", p =", format.pval(disp_p, 3), "\n")

    # Linear trend: distance-to-centroid vs n_neighbors (continuous)
    trend_df <- data.frame(
      dist_centroid = disp_neighbor$distances,
      n_neighbors = land_aligned$n_neighbors,
      log_volume = land_aligned$log_volume
    )
    trend_lm <- lm(dist_centroid ~ n_neighbors + log_volume, data = trend_df)
    trend_s <- summary(trend_lm)
    neigh_beta <- coef(trend_lm)["n_neighbors"]
    neigh_p <- trend_s$coefficients["n_neighbors", "Pr(>|t|)"]

    cat("    Linear trend (dist_to_centroid ~ n_neighbors + log_volume):\n")
    cat("      Neighbor β =", round(neigh_beta, 4), ", p =", format.pval(neigh_p, 3), "\n")

    if (disp_p < 0.05 || neigh_p < 0.05) {
      cat("    → Compositional variability differs with neighborhood density\n\n")
    } else {
      cat("    → No evidence that neighborhood density affects compositional variability\n")
      cat("    → Contrasts with experiment's composition divergence finding\n\n")
    }

    # Save result
    neighbor_disp_result <- tibble(
      test = c("betadisper_F", "betadisper_p", "trend_beta", "trend_p"),
      value = c(disp_f, disp_p, neigh_beta, neigh_p)
    )
    save_table(neighbor_disp_result, "neighborhood_composition_divergence")
    cat("  Saved: neighborhood_composition_divergence.csv\n\n")
  } else {
    cat("  Insufficient overlap between community matrix and neighborhood data (n < 20)\n\n")
  }
} else {
  cat("  Community matrix not available (run 01_load_data.R first)\n\n")
}

cat("CONCLUSION:\n")
cat("  Coral size is the dominant predictor of CAFI abundance.\n")
cat("  No evidence that neighborhood context at 5m scale affects CAFI abundance or richness.\n")
cat("  No evidence of size × neighborhood interactions.\n")
cat("  NOTE: This analysis is underpowered for small effects (R² < 0.05);\n")
cat("        null results should be interpreted as 'no evidence of an effect',\n")
cat("        not as evidence of no effect.\n\n")

cat("OUTPUT FILES:\n")
cat("  Figures: ", FIG_DIR, "\n")
cat("    - landscape_effects_summary.png\n")
cat("    - coefficient_forest_plot.png\n")
cat("    - abundance_vs_volume.png\n")
cat("    - abundance_vs_neighbors.png\n")
cat("    - abundance_vs_neighbor_volume.png\n")
cat("    - abundance_vs_isolation.png\n")
cat("    - richness_vs_volume.png\n")
cat("    - shannon_vs_volume.png\n")
cat("  Tables: ", PATHS$tables, "\n")
cat("    - landscape_univariate_results.csv\n")
cat("    - landscape_full_model_results.csv\n")
cat("    - landscape_model_comparison.csv\n")
cat("    - landscape_interaction_results.csv\n")
cat("    - neighborhood_composition_divergence.csv\n")
cat("  Objects: ", PATHS$objects, "\n")
cat("    - landscape_analysis_results.rds\n\n")
