# ============================================================================
# 05_species_scaling_analysis.R - Comprehensive CAFI-Coral Scaling Analysis
# ============================================================================
#
# PURPOSE: Test how CAFI abundance scales with coral colony volume
#
# THEORETICAL FRAMEWORK:
# ----------------------
# The relationship between habitat (coral) size and inhabitant (CAFI) abundance
# follows a power law: N = a × V^β, where:
#   - N = abundance (total CAFI, species count, etc.)
#   - V = coral volume (cm³)
#   - a = intercept (baseline density)
#   - β = scaling exponent (the key parameter of interest)
#
# In log-log space: log(N) = log(a) + β × log(V)
#
# THREE COMPETING HYPOTHESES:
# ---------------------------
#
# H1: FIELD OF DREAMS (β = 1.0)
#     "If you build it, they will come"
#     - Abundance scales PROPORTIONALLY with volume
#     - Per-capita density constant across coral sizes
#     - Mechanism: Passive habitat filling
#
# H2: PROPAGULE REDISTRIBUTION (β < 1.0)
#     Following Stier & Osenberg (2010) Ecology
#     - Abundance scales SUBLINEARLY with volume
#     - Per-capita density DECREASES in larger corals (dilution effect)
#     - Mechanism: Propagule supply limited; larger corals "dilute" settlers
#     - Also called "dilution effect" or "propagule limitation"
#
# H3: SUPER-LINEAR SCALING (β > 1.0)
#     - Abundance scales FASTER than volume
#     - Per-capita density INCREASES in larger corals
#     - Mechanism: Larger corals disproportionately attractive
#     - Could indicate aggregation, facilitation, or preferred habitat
#
# ANALYSES IN THIS SCRIPT:
# ------------------------
# PART 1: Total CAFI Abundance (community-level)
# PART 2: Species Richness
# PART 3: Shannon Diversity
# PART 4: Individual Species (abundant taxa with sufficient prevalence)
# PART 5: Taxonomic Groups (crabs, shrimps, fishes, snails)
# PART 6: Functional Groups (Trapezia, corallivores, etc.)
#
# MODEL SPECIFICATION:
#   log(Y) ~ β × log(V) + site
#   - Negative binomial GLM for counts (handles overdispersion)
#   - Site as fixed effect (accounts for among-site variation)
#   - log10 transformation for interpretability
#
# KEY TEST:
#   H0: β = 1 (Field of Dreams)
#   HA: β ≠ 1
#   Using z-test: z = (β - 1) / SE(β)
#
# REFERENCES:
#   Stier AC, Osenberg CW (2010) Propagule redirection: habitat availability
#     reduces colonization and increases recruitment in reef fishes.
#     Ecology 91:2884-2892
#
#   Stier AC, Osenberg CW (2024) Field of dreams or propagule redistribution?
#     Current Biology 34:R1186-R1189
#
# OUTPUTS:
#   Figures: output/figures/05_scaling/
#   Tables: output/tables/
#   Report: output/reports/SCALING_ANALYSIS.html
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-10
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI-CORAL SCALING ANALYSIS\n")
cat("    Field of Dreams vs. Propagule Redistribution\n")
cat("============================================================\n\n")

# Load setup and data
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

# Use script-specific output directory
FIG_DIR <- PATHS$fig_05_scaling

# Load CAFI data
cafi_clean <- readRDS(file.path(PATHS$objects, "cafi_clean.rds"))

# ============================================================================
# HELPER FUNCTION: Fit Scaling Model
# ============================================================================

#' Fit a scaling model and test against Field of Dreams (β = 1)
#'
#' @param data Data frame with 'abundance' and 'volume' columns
#' @param response_name Character string for labeling
#' @param min_nonzero Minimum non-zero observations required (default 15)
#' @return List with model results and interpretation
fit_scaling_model <- function(data, response_name, min_nonzero = 15) {

  # Calculate summary stats
  n_total <- nrow(data)
  n_nonzero <- sum(data$abundance > 0)
  n_zero <- sum(data$abundance == 0)
  total_abundance <- sum(data$abundance)

  # Check for sufficient data
  if (n_nonzero < min_nonzero) {
    return(list(
      response = response_name,
      n_corals = n_total,
      n_nonzero = n_nonzero,
      n_zero = n_zero,
      total_abundance = total_abundance,
      beta = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      z_vs_1 = NA_real_,
      p_vs_1 = NA_real_,
      interpretation = "Insufficient data",
      model = NULL,
      converged = FALSE
    ))
  }

  # Fit negative binomial GLM
  tryCatch({
    model <- MASS::glm.nb(abundance ~ log10(volume) + site, data = data)

    # Extract coefficients
    coefs <- summary(model)$coefficients
    beta <- coefs["log10(volume)", "Estimate"]
    se <- coefs["log10(volume)", "Std. Error"]
    z_val <- coefs["log10(volume)", "z value"]
    p_val <- coefs["log10(volume)", "Pr(>|z|)"]

    # 95% CI
    ci <- confint(model, "log10(volume)", level = 0.95)

    # Test vs β = 1 (Field of Dreams hypothesis)
    z_vs_1 <- (beta - 1) / se
    p_vs_1 <- 2 * pnorm(-abs(z_vs_1))

    # Interpret based on CI and statistical test
    if (p_vs_1 < 0.05 && beta < 1) {
      interpretation <- "Redistribution (β < 1)"
    } else if (p_vs_1 < 0.05 && beta > 1) {
      interpretation <- "Super-linear (β > 1)"
    } else {
      interpretation <- "Field of Dreams (β ≈ 1)"
    }

    list(
      response = response_name,
      n_corals = n_total,
      n_nonzero = n_nonzero,
      n_zero = n_zero,
      total_abundance = total_abundance,
      beta = beta,
      se = se,
      ci_lower = ci[1],
      ci_upper = ci[2],
      z_value = z_val,
      p_value = p_val,
      z_vs_1 = z_vs_1,
      p_vs_1 = p_vs_1,
      interpretation = interpretation,
      model = model,
      converged = TRUE
    )

  }, error = function(e) {
    list(
      response = response_name,
      n_corals = n_total,
      n_nonzero = n_nonzero,
      n_zero = n_zero,
      total_abundance = total_abundance,
      beta = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      z_vs_1 = NA_real_,
      p_vs_1 = NA_real_,
      interpretation = paste("Model error:", conditionMessage(e)),
      model = NULL,
      converged = FALSE
    )
  })
}

#' Convert scaling result list to single-row tibble
result_to_row <- function(result) {
  tibble(
    response = result$response,
    n_corals = result$n_corals,
    n_nonzero = result$n_nonzero,
    n_zero = result$n_zero,
    total_abundance = result$total_abundance,
    beta = result$beta,
    se = result$se,
    ci_lower = result$ci_lower,
    ci_upper = result$ci_upper,
    z_value = result$z_value,
    p_value = result$p_value,
    z_vs_1 = result$z_vs_1,
    p_vs_1 = result$p_vs_1,
    interpretation = result$interpretation
  )
}

# ============================================================================
# PREPARE BASE DATA
# ============================================================================

cat("Preparing data for scaling analysis...\n\n")

# Base coral data with volume
coral_data <- coral_master %>%
  filter(!is.na(volume), volume > 0) %>%
  dplyr::select(coral_id, site, volume, log_volume, total_cafi, otu_richness, shannon,
                n_trapezia, n_resident_fish, n_corallivore, n_other_crab, n_shrimp, n_other,
                n_crabs, n_shrimps, n_fish, n_snails)

cat("Sample size:", nrow(coral_data), "corals with volume data\n")
cat("Sites:", paste(unique(coral_data$site), collapse = ", "), "\n")
cat("Volume range:", round(min(coral_data$volume)), "-", round(max(coral_data$volume)), "cm³\n")
cat("Total CAFI:", sum(coral_data$total_cafi), "individuals\n\n")

# ############################################################################
#                    PART 1: TOTAL CAFI ABUNDANCE
# ############################################################################

cat("############################################################\n")
cat("PART 1: TOTAL CAFI ABUNDANCE SCALING\n")
cat("############################################################\n\n")

# Prepare data
abundance_data <- coral_data %>%
  dplyr::select(coral_id, site, volume, abundance = total_cafi)

# Fit model
total_result <- fit_scaling_model(abundance_data, "Total CAFI Abundance")

cat("Model: Total CAFI ~ log10(Volume) + site (Negative Binomial GLM)\n\n")
cat("Results:\n")
cat("  Scaling exponent β =", round(total_result$beta, 3), "\n")
cat("  Standard error =", round(total_result$se, 3), "\n")
cat("  95% CI: [", round(total_result$ci_lower, 3), ",", round(total_result$ci_upper, 3), "]\n")
cat("  z-value =", round(total_result$z_value, 2), ", p =", format.pval(total_result$p_value, 3), "\n\n")

cat("Test vs Field of Dreams (H0: β = 1):\n")
cat("  z =", round(total_result$z_vs_1, 2), ", p =", format.pval(total_result$p_vs_1, 3), "\n")
cat("  Interpretation:", total_result$interpretation, "\n\n")

# Store for combined results
all_results <- list(total_abundance = total_result)

# ############################################################################
#                    PART 2: SPECIES RICHNESS SCALING
# ############################################################################

cat("############################################################\n")
cat("PART 2: SPECIES RICHNESS SCALING\n")
cat("############################################################\n\n")

richness_data <- coral_data %>%
  dplyr::select(coral_id, site, volume, abundance = otu_richness)

richness_result <- fit_scaling_model(richness_data, "Species Richness")

cat("Model: OTU Richness ~ log10(Volume) + site (Negative Binomial GLM)\n\n")
cat("Results:\n")
cat("  Scaling exponent β =", round(richness_result$beta, 3), "\n")
cat("  95% CI: [", round(richness_result$ci_lower, 3), ",", round(richness_result$ci_upper, 3), "]\n")
cat("  Test vs β = 1: z =", round(richness_result$z_vs_1, 2), ", p =", format.pval(richness_result$p_vs_1, 3), "\n")
cat("  Interpretation:", richness_result$interpretation, "\n\n")

all_results$species_richness <- richness_result

# ############################################################################
#                    PART 3: SHANNON DIVERSITY SCALING
# ############################################################################

cat("############################################################\n")
cat("PART 3: SHANNON DIVERSITY SCALING\n")
cat("############################################################\n\n")

# Shannon requires linear model (continuous, not count)
shannon_data <- coral_data %>%
  filter(shannon > 0) %>%
  dplyr::select(coral_id, site, volume, shannon)

cat("Model: Shannon H' ~ log10(Volume) + site (Linear Model)\n\n")

if (nrow(shannon_data) >= 15) {
  shannon_model <- lm(shannon ~ log10(volume) + site, data = shannon_data)
  shannon_summary <- summary(shannon_model)

  shannon_beta <- coef(shannon_model)["log10(volume)"]
  shannon_se <- shannon_summary$coefficients["log10(volume)", "Std. Error"]
  shannon_ci <- confint(shannon_model, "log10(volume)", level = 0.95)
  shannon_t <- shannon_summary$coefficients["log10(volume)", "t value"]
  shannon_p <- shannon_summary$coefficients["log10(volume)", "Pr(>|t|)"]

  # Test vs β = 1
  shannon_t_vs_1 <- (shannon_beta - 1) / shannon_se
  shannon_p_vs_1 <- 2 * pt(-abs(shannon_t_vs_1), df = shannon_summary$df[2])

  cat("Results:\n")
  cat("  Scaling exponent β =", round(shannon_beta, 3), "\n")
  cat("  95% CI: [", round(shannon_ci[1], 3), ",", round(shannon_ci[2], 3), "]\n")
  cat("  R² =", round(shannon_summary$r.squared, 3), "\n")
  cat("  Test vs β = 1: t =", round(shannon_t_vs_1, 2), ", p =", format.pval(shannon_p_vs_1, 3), "\n\n")

  shannon_result <- list(
    response = "Shannon Diversity",
    n_corals = nrow(shannon_data),
    n_nonzero = nrow(shannon_data),
    n_zero = 0,
    total_abundance = NA,
    beta = shannon_beta,
    se = shannon_se,
    ci_lower = shannon_ci[1],
    ci_upper = shannon_ci[2],
    z_value = shannon_t,
    p_value = shannon_p,
    z_vs_1 = shannon_t_vs_1,
    p_vs_1 = shannon_p_vs_1,
    interpretation = ifelse(shannon_p_vs_1 < 0.05 && shannon_beta < 1, "Redistribution (β < 1)",
                           ifelse(shannon_p_vs_1 < 0.05 && shannon_beta > 1, "Super-linear (β > 1)",
                                  "Field of Dreams (β ≈ 1)")),
    model = shannon_model,
    converged = TRUE
  )
} else {
  shannon_result <- list(
    response = "Shannon Diversity",
    n_corals = nrow(shannon_data), n_nonzero = NA, n_zero = NA, total_abundance = NA,
    beta = NA, se = NA, ci_lower = NA, ci_upper = NA, z_value = NA, p_value = NA,
    z_vs_1 = NA, p_vs_1 = NA, interpretation = "Insufficient data", model = NULL, converged = FALSE
  )
  cat("Insufficient data for Shannon scaling analysis\n\n")
}

all_results$shannon_diversity <- shannon_result

# ############################################################################
#                    PART 4: TAXONOMIC GROUP SCALING
# ############################################################################

cat("############################################################\n")
cat("PART 4: TAXONOMIC GROUP SCALING\n")
cat("############################################################\n\n")

cat("Testing scaling for major taxonomic groups...\n\n")

taxonomic_groups <- list(
  "Crabs" = "n_crabs",
  "Shrimps" = "n_shrimps",
  "Fishes" = "n_fish",
  "Snails" = "n_snails"
)

taxonomic_results <- list()

for (group_name in names(taxonomic_groups)) {
  col_name <- taxonomic_groups[[group_name]]

  group_data <- coral_data %>%
    dplyr::select(coral_id, site, volume, abundance = all_of(col_name))

  result <- fit_scaling_model(group_data, group_name)
  taxonomic_results[[group_name]] <- result

  if (result$converged && !is.na(result$beta)) {
    cat(sprintf("%-10s β = %6.3f [%5.2f, %5.2f], p vs 1: %s - %s\n",
                group_name, result$beta, result$ci_lower, result$ci_upper,
                format.pval(result$p_vs_1, 2), result$interpretation))
  } else {
    cat(sprintf("%-10s %s\n", group_name, result$interpretation))
  }
}

cat("\n")
all_results$taxonomic_groups <- taxonomic_results

# ############################################################################
#                    PART 5: FUNCTIONAL GROUP SCALING
# ############################################################################

cat("############################################################\n")
cat("PART 5: FUNCTIONAL GROUP SCALING\n")
cat("############################################################\n\n")

cat("Testing scaling for functional groups...\n\n")

functional_groups <- list(
  "Trapezia (mutualists)" = "n_trapezia",
  "Resident Fish" = "n_resident_fish",
  "Corallivores" = "n_corallivore",
  "Other Crabs" = "n_other_crab",
  "Shrimps" = "n_shrimp",
  "Other" = "n_other"
)

functional_results <- list()

for (group_name in names(functional_groups)) {
  col_name <- functional_groups[[group_name]]

  group_data <- coral_data %>%
    dplyr::select(coral_id, site, volume, abundance = all_of(col_name))

  result <- fit_scaling_model(group_data, group_name)
  functional_results[[group_name]] <- result

  if (result$converged && !is.na(result$beta)) {
    cat(sprintf("%-25s β = %6.3f [%5.2f, %5.2f], p vs 1: %s - %s\n",
                group_name, result$beta, result$ci_lower, result$ci_upper,
                format.pval(result$p_vs_1, 2), result$interpretation))
  } else {
    cat(sprintf("%-25s %s\n", group_name, result$interpretation))
  }
}

cat("\n")
all_results$functional_groups <- functional_results

# ############################################################################
#                    PART 6: INDIVIDUAL SPECIES SCALING
# ############################################################################

cat("############################################################\n")
cat("PART 6: INDIVIDUAL SPECIES SCALING\n")
cat("############################################################\n\n")

# Calculate species prevalence and abundance
species_summary <- cafi_clean %>%
  group_by(otu) %>%
  summarise(
    total_individuals = n(),
    n_corals_present = n_distinct(coral_id),
    prevalence = n_distinct(coral_id) / n_distinct(cafi_clean$coral_id) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(total_individuals))

# Get functional group info
species_functional <- cafi_clean %>%
  dplyr::select(otu, functional_group, type) %>%
  distinct()

species_summary <- species_summary %>%
  left_join(species_functional, by = "otu")

# Selection criteria: >= 30 individuals AND >= 15% prevalence
key_species <- species_summary %>%
  filter(total_individuals >= 30, prevalence >= 15) %>%
  arrange(desc(total_individuals))

cat("Selection criteria:\n")
cat("  - Minimum 30 total individuals\n")
cat("  - Minimum 15% prevalence (present on >= 17 corals)\n")
cat("  - Result:", nrow(key_species), "species selected\n\n")

# Create species × coral abundance matrix
species_by_coral <- cafi_clean %>%
  filter(otu %in% key_species$otu) %>%
  group_by(coral_id, otu) %>%
  summarise(abundance = n(), .groups = "drop") %>%
  complete(coral_id = unique(cafi_clean$coral_id), otu, fill = list(abundance = 0))

# Merge with coral data
species_scaling_data <- species_by_coral %>%
  left_join(coral_data %>% dplyr::select(coral_id, site, volume), by = "coral_id") %>%
  filter(!is.na(volume))

cat("Fitting models for", nrow(key_species), "species...\n\n")

species_results <- map(key_species$otu, function(sp) {
  sp_data <- species_scaling_data %>%
    filter(otu == sp) %>%
    dplyr::select(coral_id, site, volume, abundance)

  result <- fit_scaling_model(sp_data, sp)

  if (result$converged && !is.na(result$beta)) {
    cat(sprintf("  %-30s β = %6.3f, p vs 1: %s - %s\n",
                sp, result$beta, format.pval(result$p_vs_1, 2), result$interpretation))
  } else {
    cat(sprintf("  %-30s %s\n", sp, result$interpretation))
  }

  result
})

names(species_results) <- key_species$otu
all_results$individual_species <- species_results

# ############################################################################
#                    PART 7: SUMMARY & SYNTHESIS
# ############################################################################

cat("\n############################################################\n")
cat("PART 7: SUMMARY & SYNTHESIS\n")
cat("############################################################\n\n")

# Compile all results into a single table
compile_results <- function(results_list) {
  rows <- list()

  # Community-level metrics
  rows[[1]] <- result_to_row(results_list$total_abundance) %>% mutate(category = "Community")
  rows[[2]] <- result_to_row(results_list$species_richness) %>% mutate(category = "Community")
  rows[[3]] <- result_to_row(results_list$shannon_diversity) %>% mutate(category = "Community")

  # Taxonomic groups
  for (name in names(results_list$taxonomic_groups)) {
    rows[[length(rows) + 1]] <- result_to_row(results_list$taxonomic_groups[[name]]) %>%
      mutate(category = "Taxonomic Group")
  }

  # Functional groups
  for (name in names(results_list$functional_groups)) {
    rows[[length(rows) + 1]] <- result_to_row(results_list$functional_groups[[name]]) %>%
      mutate(category = "Functional Group")
  }

  # Individual species
  for (name in names(results_list$individual_species)) {
    rows[[length(rows) + 1]] <- result_to_row(results_list$individual_species[[name]]) %>%
      mutate(category = "Species")
  }

  bind_rows(rows)
}

all_results_df <- compile_results(all_results)

# Summary statistics
valid_results <- all_results_df %>% filter(!is.na(beta))

cat("OVERALL SUMMARY:\n")
cat("  Total analyses:", nrow(all_results_df), "\n")
cat("  Successful fits:", nrow(valid_results), "\n\n")

# Summary by interpretation
interpretation_summary <- valid_results %>%
  group_by(interpretation) %>%
  summarise(
    n = n(),
    mean_beta = mean(beta),
    .groups = "drop"
  )

cat("INTERPRETATION DISTRIBUTION:\n")
print(interpretation_summary)
cat("\n")

# Summary by category
category_summary <- valid_results %>%
  group_by(category) %>%
  summarise(
    n = n(),
    mean_beta = mean(beta),
    median_beta = median(beta),
    .groups = "drop"
  )

cat("\nBY CATEGORY:\n")
print(category_summary)
cat("\n")

# T-test for species-level results
species_betas <- valid_results %>%
  filter(category == "Species") %>%
  pull(beta)

if (length(species_betas) >= 3) {
  t_test <- t.test(species_betas, mu = 1)
  cat("\nSPECIES-LEVEL T-TEST (H0: mean β = 1):\n")
  cat("  Mean β =", round(mean(species_betas), 3), "±", round(sd(species_betas)/sqrt(length(species_betas)), 3), "\n")
  cat("  t =", round(t_test$statistic, 2), ", df =", round(t_test$parameter, 1),
      ", p =", format.pval(t_test$p.value, 3), "\n")

  if (t_test$p.value < 0.05 && mean(species_betas) < 1) {
    cat("  Conclusion: Species on average show REDISTRIBUTION (β < 1)\n")
  } else if (t_test$p.value < 0.05 && mean(species_betas) > 1) {
    cat("  Conclusion: Species on average show SUPER-LINEAR scaling (β > 1)\n")
  } else {
    cat("  Conclusion: Species on average support FIELD OF DREAMS (β ≈ 1)\n")
  }
} else {
  t_test <- NULL
}

# ############################################################################
#                    FIGURES
# ############################################################################

cat("\n############################################################\n")
cat("CREATING FIGURES\n")
cat("############################################################\n\n")

# Define interpretation colors
interpretation_colors <- c(
  "Redistribution (β < 1)" = "#E41A1C",
  "Field of Dreams (β ≈ 1)" = "#4DAF4A",
  "Super-linear (β > 1)" = "#377EB8"
)

# ============================================================================
# FIGURE 1: Community-Level Scaling (Total, Richness, Shannon)
# ============================================================================

cat("Creating Figure 1: Community-level scaling...\n")

# Total abundance plot
p_total <- ggplot(coral_data, aes(x = volume, y = total_cafi, color = site)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(aes(group = 1), method = "glm.nb", formula = y ~ log10(x),
              se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  scale_color_manual(values = SITE_COLORS) +
  labs(
    x = expression("Coral Volume (cm"^3*")"),
    y = "Total CAFI",
    title = "A. Total Abundance",
    subtitle = sprintf("β = %.2f [%.2f, %.2f]",
                       total_result$beta, total_result$ci_lower, total_result$ci_upper)
  ) +
  theme(legend.position = c(0.15, 0.85))

# Richness plot
p_richness <- ggplot(coral_data, aes(x = volume, y = otu_richness, color = site)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(aes(group = 1), method = "glm.nb", formula = y ~ log10(x),
              se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  scale_color_manual(values = SITE_COLORS) +
  labs(
    x = expression("Coral Volume (cm"^3*")"),
    y = "Species Richness",
    title = "B. Species Richness",
    subtitle = sprintf("β = %.2f [%.2f, %.2f]",
                       richness_result$beta, richness_result$ci_lower, richness_result$ci_upper)
  ) +
  theme(legend.position = "none")

# Shannon plot
p_shannon <- coral_data %>%
  filter(shannon > 0) %>%
  ggplot(aes(x = volume, y = shannon, color = site)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ log10(x),
              se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = SITE_COLORS) +
  labs(
    x = expression("Coral Volume (cm"^3*")"),
    y = "Shannon H'",
    title = "C. Shannon Diversity",
    subtitle = if(!is.na(shannon_result$beta)) sprintf("β = %.2f [%.2f, %.2f]",
                       shannon_result$beta, shannon_result$ci_lower, shannon_result$ci_upper) else "Insufficient data"
  ) +
  theme(legend.position = "none")

# Combine community plots
fig_community <- (p_total + p_richness + p_shannon) +
  plot_annotation(
    title = "Community-Level Scaling: CAFI vs Coral Volume",
    subtitle = "Negative binomial GLM | Dashed line = β = 1 (Field of Dreams)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(FIG_DIR, "community_scaling.png"), fig_community,
       width = 14, height = 5, dpi = 300, bg = "white")
cat("  Saved: community_scaling.png\n")

# ============================================================================
# FIGURE 2: Taxonomic & Functional Group Forest Plot
# ============================================================================

cat("Creating Figure 2: Taxonomic & functional group scaling...\n")

# Combine taxonomic and functional results
group_results <- bind_rows(
  map_df(names(all_results$taxonomic_groups), ~result_to_row(all_results$taxonomic_groups[[.x]])) %>%
    mutate(group_type = "Taxonomic"),
  map_df(names(all_results$functional_groups), ~result_to_row(all_results$functional_groups[[.x]])) %>%
    mutate(group_type = "Functional")
) %>%
  filter(!is.na(beta))

if (nrow(group_results) > 0) {
  # Create unique labels for plotting (append group type to avoid duplicates)
  group_results <- group_results %>%
    mutate(plot_label = paste0(response, " (", substr(group_type, 1, 1), ")"))

  p_groups <- group_results %>%
    mutate(plot_label = factor(plot_label, levels = unique(plot_label[order(beta)]))) %>%
    ggplot(aes(x = beta, y = plot_label, color = interpretation, shape = group_type)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.3, linewidth = 0.6) +
    geom_point(size = 4) +
    scale_color_manual(values = interpretation_colors, name = "Interpretation") +
    scale_shape_manual(values = c("Taxonomic" = 16, "Functional" = 17), name = "Group Type") +
    labs(
      x = expression("Scaling Exponent (" * beta * ")"),
      y = NULL,
      title = "Scaling by Taxonomic and Functional Groups",
      subtitle = "95% CI shown | Vertical line = Field of Dreams (β = 1)"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 10)
    )

  ggsave(file.path(FIG_DIR, "group_scaling_forest.png"), p_groups,
         width = 10, height = 7, dpi = 300, bg = "white")
  cat("  Saved: group_scaling_forest.png\n")
}

# ============================================================================
# FIGURE 3: Species-Level Forest Plot
# ============================================================================

cat("Creating Figure 3: Species-level scaling forest plot...\n")

species_plot_data <- valid_results %>%
  filter(category == "Species") %>%
  left_join(species_functional %>% dplyr::select(otu, functional_group, type),
            by = c("response" = "otu")) %>%
  mutate(response = factor(response, levels = response[order(beta)]))

if (nrow(species_plot_data) > 0) {
  p_species_forest <- ggplot(species_plot_data, aes(x = beta, y = response, color = interpretation)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray80") +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.3, linewidth = 0.6) +
    geom_point(size = 3) +
    scale_color_manual(values = interpretation_colors, name = "Scaling Pattern") +
    labs(
      x = expression("Scaling Exponent (" * beta * ")"),
      y = NULL,
      title = "Species-Level Scaling: Abundance vs Coral Volume",
      subtitle = sprintf("N = %d species | Mean β = %.2f | Test vs 1: p = %s",
                         nrow(species_plot_data),
                         mean(species_plot_data$beta),
                         if(!is.null(t_test)) format.pval(t_test$p.value, 2) else "NA")
    ) +
    theme(
      legend.position = c(0.85, 0.15),
      legend.background = element_rect(fill = "white", color = "gray80"),
      axis.text.y = element_text(face = "italic", size = 9)
    ) +
    coord_cartesian(xlim = c(-0.5, 3))

  ggsave(file.path(FIG_DIR, "species_scaling_forest.png"), p_species_forest,
         width = 10, height = 8, dpi = 300, bg = "white")
  cat("  Saved: species_scaling_forest.png\n")
}

# ============================================================================
# FIGURE 4: Species Scatter Panels
# ============================================================================

cat("Creating Figure 4: Species scatter panels...\n")

# Select top species by sample size
top_species <- valid_results %>%
  filter(category == "Species") %>%
  arrange(desc(n_nonzero)) %>%
  head(12)

make_species_panel <- function(sp_name, sp_result) {
  sp_data <- species_scaling_data %>%
    filter(otu == sp_name, abundance > 0)

  beta_val <- sp_result$beta
  interp <- sp_result$interpretation

  title_color <- case_when(
    interp == "Redistribution (β < 1)" ~ "#E41A1C",
    interp == "Super-linear (β > 1)" ~ "#377EB8",
    TRUE ~ "#4DAF4A"
  )

  ggplot(sp_data, aes(x = volume, y = abundance, color = site)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(aes(group = 1), method = "glm.nb", formula = y ~ log10(x),
                se = TRUE, color = "black", linewidth = 0.8) +
    scale_x_log10(labels = scales::comma) +
    scale_y_continuous() +
    scale_color_manual(values = SITE_COLORS) +
    labs(
      x = expression("Volume (cm"^3*")"),
      y = "Abundance",
      title = str_replace(sp_name, "_", " ")
    ) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("β = %.2f", beta_val),
             hjust = 1.1, vjust = 1.3, size = 3, fontface = "italic") +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold.italic", size = 9, color = title_color)
    )
}

if (nrow(top_species) > 0) {
  species_panels <- map2(top_species$response,
                          map(top_species$response, ~all_results$individual_species[[.x]]),
                          make_species_panel)

  fig_species <- wrap_plots(species_panels, ncol = 4) +
    plot_annotation(
      title = "Species-Specific Scaling Relationships",
      subtitle = "Top 12 species by sample size | Title color indicates interpretation",
      caption = paste0("Red = Redistribution (β < 1), Green = Field of Dreams (β ≈ 1), Blue = Super-linear (β > 1)"),
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11)
      )
    )

  ggsave(file.path(FIG_DIR, "species_scaling_panels.png"), fig_species,
         width = 14, height = 10, dpi = 300, bg = "white")
  cat("  Saved: species_scaling_panels.png\n")
}

# ============================================================================
# FIGURE 5: Summary Visualization
# ============================================================================

cat("Creating Figure 5: Summary visualization...\n")

# Beta distribution by category
p_beta_dist <- valid_results %>%
  filter(category %in% c("Species", "Taxonomic Group", "Functional Group")) %>%
  ggplot(aes(x = beta, fill = category)) +
  geom_histogram(bins = 15, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
  scale_fill_viridis_d(option = "D", name = "Category") +
  labs(
    x = expression("Scaling Exponent (" * beta * ")"),
    y = "Count",
    title = "Distribution of Scaling Exponents",
    subtitle = "Dashed line = Field of Dreams (β = 1)"
  ) +
  theme(legend.position = "right")

# Beta by interpretation
p_interp <- valid_results %>%
  ggplot(aes(x = interpretation, y = beta, fill = interpretation)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = interpretation_colors) +
  labs(
    x = NULL,
    y = expression("Scaling Exponent (" * beta * ")"),
    title = "Scaling by Interpretation"
  ) +
  theme(legend.position = "none") +
  coord_flip()

fig_summary <- (p_beta_dist / p_interp) +
  plot_annotation(
    title = "Scaling Analysis Summary",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(FIG_DIR, "scaling_summary.png"), fig_summary,
       width = 10, height = 10, dpi = 300, bg = "white")
cat("  Saved: scaling_summary.png\n")

# ############################################################################
#                    SAVE TABLES
# ############################################################################

cat("\n############################################################\n")
cat("SAVING TABLES\n")
cat("############################################################\n\n")

# Full results table
full_table <- all_results_df %>%
  mutate(across(where(is.numeric) & !contains("n_"), ~round(., 3))) %>%
  arrange(category, desc(n_nonzero))

save_table(full_table, "scaling_results_all")
cat("  Saved: scaling_results_all.csv\n")

# Summary by category
save_table(category_summary, "scaling_summary_by_category")
cat("  Saved: scaling_summary_by_category.csv\n")

# Interpretation summary
save_table(interpretation_summary, "scaling_interpretation_summary")
cat("  Saved: scaling_interpretation_summary.csv\n")

# ############################################################################
#                    SAVE RESULTS OBJECT
# ############################################################################

scaling_analysis_results <- list(
  # Sample info
  n_corals = nrow(coral_data),
  total_cafi = sum(coral_data$total_cafi),

  # Full results table
  all_results = all_results_df,

  # Model objects
  models = all_results,

  # Summary statistics
  category_summary = category_summary,
  interpretation_summary = interpretation_summary,

  # T-test for species
  species_t_test = t_test,

  # Species metadata
  species_summary = species_summary
)

save_object(scaling_analysis_results, "scaling_analysis_results")
cat("  Saved: scaling_analysis_results.rds\n")

# ############################################################################
#                    FINAL SUMMARY
# ############################################################################

cat("\n")
cat("============================================================\n")
cat("    SCALING ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("SAMPLE:\n")
cat("  Corals:", nrow(coral_data), "\n")
cat("  Total CAFI:", sum(coral_data$total_cafi), "\n")
cat("  Species analyzed:", sum(valid_results$category == "Species"), "\n\n")

cat("KEY RESULTS:\n\n")

cat("1. TOTAL ABUNDANCE:\n")
cat("   β =", round(total_result$beta, 3), "[", round(total_result$ci_lower, 3), ",",
    round(total_result$ci_upper, 3), "]\n")
cat("   Interpretation:", total_result$interpretation, "\n\n")

cat("2. SPECIES RICHNESS:\n")
cat("   β =", round(richness_result$beta, 3), "[", round(richness_result$ci_lower, 3), ",",
    round(richness_result$ci_upper, 3), "]\n")
cat("   Interpretation:", richness_result$interpretation, "\n\n")

cat("3. SPECIES-LEVEL PATTERNS:\n")
if (!is.null(t_test)) {
  cat("   Mean β =", round(mean(species_betas), 3), "\n")
  cat("   T-test (H0: β = 1): p =", format.pval(t_test$p.value, 3), "\n")
}
cat("   Distribution:\n")
for (i in 1:nrow(interpretation_summary)) {
  cat("     ", interpretation_summary$interpretation[i], ":", interpretation_summary$n[i], "\n")
}

cat("\nOUTPUT FILES:\n")
cat("  Figures: output/figures/05_scaling/\n")
cat("    - community_scaling.png\n")
cat("    - group_scaling_forest.png\n")
cat("    - species_scaling_forest.png\n")
cat("    - species_scaling_panels.png\n")
cat("    - scaling_summary.png\n")
cat("  Tables: output/tables/\n")
cat("    - scaling_results_all.csv\n")
cat("    - scaling_summary_by_category.csv\n")
cat("    - scaling_interpretation_summary.csv\n\n")
