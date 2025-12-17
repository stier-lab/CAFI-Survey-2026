#!/usr/bin/env Rscript
# ============================================================================
# DIAGNOSTIC: Is the Richness-Condition Relationship Real or a Sampling Artifact?
#
# This script rigorously tests whether any observed richness-condition
# relationship is driven by:
#   1. SAMPLING EFFECTS: More individuals → more species detected (rarefaction)
#   2. KEY SPECIES: Specific taxa driving the pattern (leave-one-out analysis)
#   3. VOLUME CONFOUNDS: Larger corals have both more CAFI and different condition
#
# Author: CAFI Analysis Pipeline
# Date: 2025-12-01
# ============================================================================

cat("============================================================\n")
cat("DIAGNOSTIC: Richness-Condition Relationship Analysis\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(lme4)
  library(car)
})

# Load data
coral_data <- read.csv("data/survey_coral_characteristics_merged_v2.csv")
cafi_data <- read.csv("data/survey_cafi_data_w_taxonomy_summer2019_v5.csv")

# Check if condition scores exist
condition_file <- "output/tables/coral_condition_scores.csv"
if (!file.exists(condition_file)) {
  # Try RDS
  condition_rds <- "output/objects/coral_condition_scores.rds"
  if (file.exists(condition_rds)) {
    condition_data <- readRDS(condition_rds)
  } else {
    stop("Cannot find condition scores file!")
  }
} else {
  condition_data <- read.csv(condition_file)
}

cat("Data loaded:\n")
cat("  - Coral characteristics:", nrow(coral_data), "rows\n")
cat("  - CAFI observations:", nrow(cafi_data), "rows\n")
cat("  - Condition scores:", nrow(condition_data), "rows\n\n")

# ============================================================================
# DATA PREPARATION
# ============================================================================

cat("--- PREPARING DATA ---\n\n")

# Process coral data
coral_processed <- coral_data %>%
  mutate(
    site = str_extract(site, "^[A-Z]+"),
    volume = coalesce(volume_field, volume_lab),
    log_volume = log10(volume + 1)
  ) %>%
  filter(!is.na(volume), volume > 0, site %in% c("HAU", "MAT", "MRB"))

# Process CAFI data - create species identifier
cafi_processed <- cafi_data %>%
  filter(!is.na(genus) & genus != "" & genus != "NA") %>%
  mutate(
    species_id = paste(genus, species, sep = "_"),
    species_id = ifelse(is.na(species) | species == "" | species == "NA",
                        genus, species_id)
  )

# Create community matrix
comm_matrix <- cafi_processed %>%
  group_by(coral_id, species_id) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_id, values_from = count, values_fill = 0)

comm_mat <- as.matrix(comm_matrix[,-1])
rownames(comm_mat) <- comm_matrix$coral_id

cat("Community matrix:", nrow(comm_mat), "corals x", ncol(comm_mat), "species\n\n")

# Calculate diversity metrics
div_data <- data.frame(
  coral_id = comm_matrix$coral_id,
  abundance = rowSums(comm_mat),
  richness = rowSums(comm_mat > 0),
  shannon = diversity(comm_mat, "shannon"),
  simpson = diversity(comm_mat, "simpson")
)

# Calculate evenness (Pielou's J)
div_data$evenness <- div_data$shannon / log(div_data$richness)
div_data$evenness[is.nan(div_data$evenness) | is.infinite(div_data$evenness)] <- NA

# Merge all data
analysis_data <- div_data %>%
  left_join(coral_processed %>% select(coral_id, site, volume, log_volume), by = "coral_id") %>%
  left_join(condition_data %>% select(coral_id, condition_score), by = "coral_id") %>%
  filter(!is.na(volume), !is.na(condition_score))

cat("Analysis dataset:", nrow(analysis_data), "corals with complete data\n\n")

# ============================================================================
# PART 1: IS THIS A SAMPLING EFFECT?
# ============================================================================

cat("============================================================\n")
cat("PART 1: TESTING FOR SAMPLING EFFECTS\n")
cat("============================================================\n\n")

# 1a. Correlation between abundance and richness
cat("1a. Species-Abundance Relationship:\n")
abund_rich_cor <- cor.test(div_data$abundance, div_data$richness)
cat(sprintf("    Abundance-Richness correlation: r = %.3f, p < 0.001\n",
            abund_rich_cor$estimate))
cat("    → Strong positive correlation indicates sampling effect potential\n\n")

# 1b. How much of richness variance is explained by abundance?
rich_abund_model <- lm(richness ~ log10(abundance + 1), data = analysis_data)
cat("1b. Richness ~ log(Abundance) regression:\n")
cat(sprintf("    R² = %.3f (%.1f%% of richness variance explained by abundance)\n",
            summary(rich_abund_model)$r.squared,
            summary(rich_abund_model)$r.squared * 100))
cat("\n")

# 1c. Rarefied richness
cat("1c. Rarefied Richness (standardized sampling):\n")

# Find minimum sample size for rarefaction
min_n <- min(analysis_data$abundance)
cat(sprintf("    Minimum abundance: %d\n", min_n))

# Use a reasonable rarefaction depth (10 is standard for many studies)
rare_depth <- 10
valid_for_rare <- analysis_data$abundance >= rare_depth

cat(sprintf("    Using rarefaction depth: %d\n", rare_depth))
cat(sprintf("    Corals with >= %d individuals: %d (%.1f%%)\n",
            rare_depth, sum(valid_for_rare), 100 * mean(valid_for_rare)))

# Calculate rarefied richness
analysis_data$rarefied_richness <- NA
if (sum(valid_for_rare) > 0) {
  comm_subset <- comm_mat[analysis_data$coral_id[valid_for_rare], ]
  analysis_data$rarefied_richness[valid_for_rare] <- rarefy(comm_subset, rare_depth)
}

# Compare raw vs rarefied correlations with condition
cat("\n1d. Comparing Raw vs Rarefied Richness relationships with Condition:\n")

raw_cor <- cor.test(analysis_data$richness, analysis_data$condition_score)
rare_cor <- cor.test(analysis_data$rarefied_richness, analysis_data$condition_score,
                     use = "complete.obs")

cat(sprintf("    Raw richness vs condition:      r = %.3f, p = %.4f\n",
            raw_cor$estimate, raw_cor$p.value))
cat(sprintf("    Rarefied richness vs condition: r = %.3f, p = %.4f\n",
            rare_cor$estimate, rare_cor$p.value))

if (abs(raw_cor$estimate) > abs(rare_cor$estimate) * 1.5) {
  cat("    → SAMPLING EFFECT DETECTED: Raw correlation much stronger than rarefied\n")
} else if (rare_cor$p.value < 0.05) {
  cat("    → Real effect persists after rarefaction\n")
} else {
  cat("    → Effect disappears after controlling for sampling\n")
}
cat("\n")

# ============================================================================
# PART 2: KEY SPECIES ANALYSIS
# ============================================================================

cat("============================================================\n")
cat("PART 2: KEY SPECIES DRIVING THE RELATIONSHIP\n")
cat("============================================================\n\n")

# 2a. Which species correlate with condition?
cat("2a. Individual Species-Condition Correlations:\n")
cat("    (Testing which species abundances predict coral condition)\n\n")

# Get species that are present in enough corals
species_prevalence <- colSums(comm_mat > 0)
common_species <- names(species_prevalence[species_prevalence >= 10])

cat(sprintf("    Testing %d species present in >= 10 corals\n\n", length(common_species)))

# Test each species
species_condition_cors <- data.frame(
  species = character(),
  prevalence = numeric(),
  mean_abundance = numeric(),
  correlation = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

comm_df <- as.data.frame(comm_mat)
comm_df$coral_id <- rownames(comm_mat)

merged_for_species <- merge(comm_df,
                            analysis_data[, c("coral_id", "condition_score", "volume")],
                            by = "coral_id")

for (sp in common_species) {
  sp_abund <- merged_for_species[[sp]]
  cond <- merged_for_species$condition_score

  # Only test if there's variance
  if (sd(sp_abund) > 0) {
    test <- cor.test(sp_abund, cond)
    species_condition_cors <- rbind(species_condition_cors, data.frame(
      species = sp,
      prevalence = sum(sp_abund > 0),
      mean_abundance = mean(sp_abund),
      correlation = test$estimate,
      p_value = test$p.value
    ))
  }
}

# Sort by absolute correlation
species_condition_cors <- species_condition_cors %>%
  arrange(p_value) %>%
  mutate(
    significant = p_value < 0.05,
    fdr_p = p.adjust(p_value, method = "fdr")
  )

cat("    Top 10 species by condition correlation:\n")
print(head(species_condition_cors %>%
             select(species, prevalence, correlation, p_value, fdr_p), 10))

n_sig_raw <- sum(species_condition_cors$p_value < 0.05)
n_sig_fdr <- sum(species_condition_cors$fdr_p < 0.05)
cat(sprintf("\n    Significant correlations: %d (raw p < 0.05), %d (FDR < 0.05)\n",
            n_sig_raw, n_sig_fdr))

# 2b. Leave-one-species-out analysis
cat("\n2b. Leave-One-Species-Out Analysis:\n")
cat("    Testing if removing any single species eliminates the richness-condition effect\n\n")

# Base correlation
base_cor <- cor(analysis_data$richness, analysis_data$condition_score)

# Test each species removal
loo_results <- data.frame(
  species_removed = character(),
  new_correlation = numeric(),
  change = numeric(),
  stringsAsFactors = FALSE
)

for (sp in common_species) {
  # Recalculate richness without this species
  reduced_mat <- comm_mat[analysis_data$coral_id, colnames(comm_mat) != sp, drop = FALSE]
  new_richness <- rowSums(reduced_mat > 0)

  new_cor <- cor(new_richness, analysis_data$condition_score)

  loo_results <- rbind(loo_results, data.frame(
    species_removed = sp,
    new_correlation = new_cor,
    change = new_cor - base_cor
  ))
}

# Find species whose removal has biggest impact
loo_results <- loo_results %>% arrange(change)

cat("    Species whose removal most DECREASES the correlation (potential drivers):\n")
print(head(loo_results %>% select(species_removed, new_correlation, change), 5))

cat("\n    Species whose removal most INCREASES the correlation (potential suppressors):\n")
print(tail(loo_results %>% select(species_removed, new_correlation, change), 5))

# Check if any single species is critical
max_change <- max(abs(loo_results$change))
if (max_change > 0.05) {
  cat(sprintf("\n    → Maximum change from removing a species: %.3f\n", max_change))
  cat("    → Some species have notable influence on the relationship\n")
} else {
  cat("\n    → No single species dramatically changes the correlation\n")
  cat("    → Effect is distributed across multiple species\n")
}

# ============================================================================
# PART 3: TAXONOMIC GROUP ANALYSIS
# ============================================================================

cat("\n============================================================\n")
cat("PART 3: TAXONOMIC GROUP CONTRIBUTIONS\n")
cat("============================================================\n\n")

# Calculate group-specific richness
group_richness <- cafi_processed %>%
  group_by(coral_id, type) %>%
  summarise(
    group_richness = n_distinct(species_id),
    group_abundance = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = c(group_richness, group_abundance),
    values_fill = 0
  )

# Merge with analysis data
analysis_with_groups <- analysis_data %>%
  left_join(group_richness, by = "coral_id")

# Test each group
cat("Group-specific richness correlations with condition:\n\n")

group_types <- c("crab", "shrimp", "fish", "snail")
for (grp in group_types) {
  rich_col <- paste0("group_richness_", grp)
  if (rich_col %in% names(analysis_with_groups)) {
    test <- cor.test(analysis_with_groups[[rich_col]],
                     analysis_with_groups$condition_score)
    cat(sprintf("    %s richness:  r = %.3f, p = %.4f\n",
                str_to_title(grp), test$estimate, test$p.value))
  }
}

# ============================================================================
# PART 4: VOLUME-CONTROLLED ANALYSIS
# ============================================================================

cat("\n============================================================\n")
cat("PART 4: CONTROLLING FOR CORAL VOLUME\n")
cat("============================================================\n\n")

cat("Volume could confound richness-condition because:\n")
cat("  - Larger corals have more CAFI → more species detected\n")
cat("  - Larger corals may have different physiological status\n\n")

# 4a. Simple correlations
cat("4a. Correlations with log(volume):\n")
vol_abund_cor <- cor.test(analysis_data$abundance, analysis_data$log_volume)
vol_rich_cor <- cor.test(analysis_data$richness, analysis_data$log_volume)
vol_cond_cor <- cor.test(analysis_data$condition_score, analysis_data$log_volume)

cat(sprintf("    Abundance ~ log(volume): r = %.3f, p = %.4f\n",
            vol_abund_cor$estimate, vol_abund_cor$p.value))
cat(sprintf("    Richness ~ log(volume):  r = %.3f, p = %.4f\n",
            vol_rich_cor$estimate, vol_rich_cor$p.value))
cat(sprintf("    Condition ~ log(volume): r = %.3f, p = %.4f\n",
            vol_cond_cor$estimate, vol_cond_cor$p.value))

# 4b. Partial correlation (richness-condition controlling for volume)
cat("\n4b. Partial Correlations (controlling for volume):\n")

# Residualize both variables on volume
rich_resid <- residuals(lm(richness ~ log_volume, data = analysis_data))
cond_resid <- residuals(lm(condition_score ~ log_volume, data = analysis_data))

partial_cor <- cor.test(rich_resid, cond_resid)
cat(sprintf("    Richness-Condition partial r (controlling volume): %.3f, p = %.4f\n",
            partial_cor$estimate, partial_cor$p.value))

# 4c. Multiple regression
cat("\n4c. Multiple Regression: Condition ~ Richness + log(Volume)\n")
m_multi <- lm(condition_score ~ richness + log_volume, data = analysis_data)
cat("    Coefficients:\n")
print(summary(m_multi)$coefficients[, c(1, 2, 4)])

# ============================================================================
# PART 5: FINAL DIAGNOSTIC SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("DIAGNOSTIC SUMMARY\n")
cat("============================================================\n\n")

# Collect key results
results_summary <- data.frame(
  Test = c(
    "Raw richness-condition correlation",
    "Rarefied richness-condition correlation",
    "Partial correlation (controlling volume)",
    "Richness effect in multiple regression"
  ),
  Estimate = c(
    raw_cor$estimate,
    rare_cor$estimate,
    partial_cor$estimate,
    coef(m_multi)["richness"]
  ),
  P_value = c(
    raw_cor$p.value,
    rare_cor$p.value,
    partial_cor$p.value,
    summary(m_multi)$coefficients["richness", 4]
  )
)

results_summary$Significant <- results_summary$P_value < 0.05

print(results_summary)

cat("\n--- INTERPRETATION ---\n\n")

# Determine pattern
if (raw_cor$p.value < 0.05 && rare_cor$p.value >= 0.05) {
  cat("CONCLUSION: SAMPLING ARTIFACT\n")
  cat("  The apparent richness-condition relationship disappears when\n")
  cat("  richness is standardized for sampling effort (rarefaction).\n")
  cat("  This indicates the raw correlation is driven by: more individuals\n")
  cat("  → more species detected, not a true diversity-condition link.\n")
} else if (raw_cor$p.value < 0.05 && rare_cor$p.value < 0.05) {
  cat("CONCLUSION: REAL EFFECT (persists after controlling sampling)\n")
  cat("  The richness-condition relationship remains significant even\n")
  cat("  after rarefaction, suggesting a genuine ecological pattern.\n")

  # Check for key species
  if (n_sig_fdr > 0) {
    cat("\n  KEY SPECIES IDENTIFIED:\n")
    key_species <- species_condition_cors %>% filter(fdr_p < 0.05)
    for (i in 1:min(5, nrow(key_species))) {
      cat(sprintf("    - %s (r = %.3f)\n",
                  key_species$species[i], key_species$correlation[i]))
    }
  }
} else {
  cat("CONCLUSION: NO SIGNIFICANT RELATIONSHIP\n")
  cat("  Neither raw nor rarefied richness shows a significant\n")
  cat("  correlation with coral condition.\n")
}

# Save results
output_dir <- "output/tables"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(results_summary,
          file.path(output_dir, "richness_condition_diagnostic.csv"),
          row.names = FALSE)
write.csv(species_condition_cors,
          file.path(output_dir, "species_condition_correlations.csv"),
          row.names = FALSE)
write.csv(loo_results,
          file.path(output_dir, "leave_one_species_out.csv"),
          row.names = FALSE)

cat("\n\nResults saved to output/tables/:\n")
cat("  - richness_condition_diagnostic.csv\n")
cat("  - species_condition_correlations.csv\n")
cat("  - leave_one_species_out.csv\n")

cat("\n✅ Diagnostic analysis complete!\n")
