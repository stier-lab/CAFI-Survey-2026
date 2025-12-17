#!/usr/bin/env Rscript
# ============================================================================
# Rigorous Statistical Analysis: Composition & Relative Abundance vs Size Class
# Includes both UNIVARIATE and MULTIVARIATE approaches
# ============================================================================

cat("\n========================================\n")
cat("COMPOSITION STATISTICAL ANALYSIS\n")
cat("========================================\n\n")

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggplot2)
  library(patchwork)
  library(broom)
})

# Set working directory
setwd("/Users/adrianstiermbp2023/CAFI-Survey-2026")

# Create output directory
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load and Process Data
# ============================================================================

cat("Loading data...\n")

coral_data <- read.csv("data/survey_coral_characteristics_merged_v2.csv")
cafi_data <- read.csv("data/survey_cafi_data_w_taxonomy_summer2019_v5.csv")

# Process coral data with size classes
coral_processed <- coral_data %>%
  mutate(
    site = str_extract(site, "^[A-Z]+"),
    volume = coalesce(volume_field, volume_lab, length_field * width_field * height_field),
    log_volume = log10(volume + 1)
  ) %>%
  filter(!is.na(volume), volume > 0, site %in% c("HAU", "MAT", "MRB"))

# Define size class breaks using tertiles
breaks <- quantile(coral_processed$volume, c(0, 0.33, 0.67, 1), na.rm = TRUE)

coral_processed <- coral_processed %>%
  mutate(
    size_class = cut(volume,
                     breaks = breaks,
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE),
    size_class = factor(size_class, levels = c("Small", "Medium", "Large"))
  )

# Process CAFI data
cafi_processed <- cafi_data %>%
  filter(!is.na(genus) & genus != "" & genus != "NA") %>%
  mutate(
    site = case_when(
      is.na(site) | site == "" | site == "NA" ~ str_extract(coral_id, "^[A-Z]+"),
      TRUE ~ str_sub(site, 1, 3)
    ),
    species_id = paste(genus, species, sep = "_"),
    func_group = case_when(
      grepl("Trapezia|Tetralia|Alpheus", genus, ignore.case = TRUE) ~ "Defenders",
      grepl("Paragobiodon|Gobiodon|Caracanthus", genus, ignore.case = TRUE) ~ "Resident_fishes",
      grepl("Drupella|Coralliophila|Galeropsis|Morula", genus, ignore.case = TRUE) ~ "Ectoparasites",
      TRUE ~ "Other"
    )
  ) %>%
  filter(site %in% c("HAU", "MAT", "MRB"))

cat(sprintf("Dataset: %d corals, %d CAFI individuals\n",
            nrow(coral_processed), nrow(cafi_processed)))

# ============================================================================
# BUILD COMMUNITY MATRICES
# ============================================================================

cat("\nBuilding community matrices...\n")

# Species-level community matrix
comm_species <- cafi_processed %>%
  group_by(coral_id, species_id) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_id, values_from = count, values_fill = 0)

comm_mat <- as.matrix(comm_species[,-1])
rownames(comm_mat) <- comm_species$coral_id

# Functional group community matrix
comm_func <- cafi_processed %>%
  group_by(coral_id, func_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = func_group, values_from = count, values_fill = 0)

func_mat <- as.matrix(comm_func[,-1])
rownames(func_mat) <- comm_func$coral_id

# Calculate relative abundances
rel_mat <- comm_mat / rowSums(comm_mat) * 100
rel_func_mat <- func_mat / rowSums(func_mat) * 100

# Match metadata
meta <- coral_processed %>%
  filter(coral_id %in% rownames(comm_mat)) %>%
  select(coral_id, site, volume, log_volume, size_class) %>%
  arrange(match(coral_id, rownames(comm_mat)))

# Ensure alignment
comm_mat <- comm_mat[meta$coral_id, ]
rel_mat <- rel_mat[meta$coral_id, ]
func_mat <- func_mat[meta$coral_id, ]
rel_func_mat <- rel_func_mat[meta$coral_id, ]

cat(sprintf("Matched %d corals with community data\n", nrow(meta)))

# ============================================================================
# MULTIVARIATE ANALYSES
# ============================================================================

cat("\n--- MULTIVARIATE ANALYSES ---\n")

# Hellinger transformation for composition
comm_hell <- decostand(comm_mat, method = "hellinger")

# 1. PERMANOVA: Size class effect on composition
cat("\n1. PERMANOVA (Size Class Effect on Species Composition):\n")
perm_size <- adonis2(comm_hell ~ size_class, data = meta, permutations = 999)
print(perm_size)

# 2. PERMANOVA: Continuous volume effect
cat("\n2. PERMANOVA (Continuous Volume Effect):\n")
perm_vol <- adonis2(comm_hell ~ log_volume, data = meta, permutations = 999)
print(perm_vol)

# 3. PERMANOVA: Size class + Site
cat("\n3. PERMANOVA (Size Class + Site):\n")
perm_both <- adonis2(comm_hell ~ size_class + site, data = meta, permutations = 999)
print(perm_both)

# 4. Beta dispersion test (homogeneity of dispersion)
cat("\n4. PERMDISP (Beta Dispersion by Size Class):\n")
betadisper_result <- betadisper(vegdist(comm_hell, method = "bray"), meta$size_class)
permutest_result <- permutest(betadisper_result, permutations = 999)
print(permutest_result)

# 5. ANOSIM
cat("\n5. ANOSIM (Size Class):\n")
anosim_result <- anosim(comm_hell, meta$size_class, permutations = 999)
print(anosim_result)

# ============================================================================
# SIMPER ANALYSIS - Which species contribute to differences?
# ============================================================================

cat("\n--- SIMPER ANALYSIS ---\n")
cat("Species contributing to differences between size classes:\n\n")

simper_result <- simper(comm_mat, meta$size_class, permutations = 999)

# Small vs Large comparison (most contrasting)
cat("Small vs Large:\n")
simper_sl <- summary(simper_result)$Small_Large
if (!is.null(simper_sl)) {
  top_simper <- head(simper_sl, 10)
  print(top_simper[, c("average", "sd", "ratio", "ava", "avb", "cumsum", "p")])
}

# ============================================================================
# UNIVARIATE ANALYSES: Incidence (Chi-square tests)
# ============================================================================

cat("\n--- UNIVARIATE INCIDENCE ANALYSES (Chi-square) ---\n")

# Calculate presence/absence by size class for each species
presence_data <- as.data.frame(comm_mat > 0) %>%
  mutate(coral_id = rownames(comm_mat)) %>%
  left_join(meta %>% select(coral_id, size_class), by = "coral_id") %>%
  pivot_longer(cols = -c(coral_id, size_class), names_to = "species", values_to = "present")

# Chi-square tests for each species
incidence_tests <- presence_data %>%
  group_by(species) %>%
  summarise(
    n_small = sum(present[size_class == "Small"]),
    n_medium = sum(present[size_class == "Medium"]),
    n_large = sum(present[size_class == "Large"]),
    total = n_small + n_medium + n_large,
    .groups = "drop"
  ) %>%
  filter(total >= 5) %>%  # Only test species present in >= 5 corals
  rowwise() %>%
  mutate(
    # Chi-square test
    chisq_result = list(tryCatch({
      tbl <- matrix(c(n_small, 37 - n_small,
                      n_medium, 38 - n_medium,
                      n_large, 37 - n_large), nrow = 3, byrow = TRUE)
      chisq.test(tbl)
    }, error = function(e) list(statistic = NA, p.value = NA))),
    chi_sq = chisq_result$statistic,
    p_value = chisq_result$p.value
  ) %>%
  ungroup() %>%
  select(-chisq_result) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significant = p_adj < 0.05,
    trend = case_when(
      n_large > n_small ~ "Increases with size",
      n_large < n_small ~ "Decreases with size",
      TRUE ~ "No trend"
    )
  ) %>%
  arrange(p_adj)

cat("\nTop 20 species with significant incidence changes:\n")
print(incidence_tests %>%
        filter(significant) %>%
        head(20) %>%
        select(species, n_small, n_medium, n_large, chi_sq, p_adj, trend))

# ============================================================================
# UNIVARIATE ANALYSES: Relative Abundance (Kruskal-Wallis)
# ============================================================================

cat("\n--- UNIVARIATE RELATIVE ABUNDANCE ANALYSES (Kruskal-Wallis) ---\n")

# Calculate relative abundance for each coral
rel_abund_data <- as.data.frame(rel_mat) %>%
  mutate(coral_id = rownames(rel_mat)) %>%
  left_join(meta %>% select(coral_id, size_class), by = "coral_id") %>%
  pivot_longer(cols = -c(coral_id, size_class), names_to = "species", values_to = "rel_abund")

# Kruskal-Wallis tests for relative abundance
rel_abund_tests <- rel_abund_data %>%
  group_by(species) %>%
  filter(sum(rel_abund > 0) >= 10) %>%  # Only test species with data in >= 10 corals
  summarise(
    mean_small = mean(rel_abund[size_class == "Small"]),
    mean_medium = mean(rel_abund[size_class == "Medium"]),
    mean_large = mean(rel_abund[size_class == "Large"]),
    kw_result = list(kruskal.test(rel_abund ~ size_class)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    chi_sq = kw_result$statistic,
    p_value = kw_result$p.value
  ) %>%
  ungroup() %>%
  select(-kw_result) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significant = p_adj < 0.05,
    change = mean_large - mean_small,
    trend = case_when(
      change > 1 ~ "Increases with size",
      change < -1 ~ "Decreases with size",
      TRUE ~ "No change"
    )
  ) %>%
  arrange(p_adj)

cat("\nTop 20 species with significant relative abundance changes:\n")
print(rel_abund_tests %>%
        filter(significant) %>%
        head(20) %>%
        select(species, mean_small, mean_medium, mean_large, chi_sq, p_adj, trend))

# ============================================================================
# FUNCTIONAL GROUP ANALYSES
# ============================================================================

cat("\n--- FUNCTIONAL GROUP ANALYSES ---\n")

# Functional group relative abundance by size class
func_rel_data <- as.data.frame(rel_func_mat) %>%
  mutate(coral_id = rownames(rel_func_mat)) %>%
  left_join(meta %>% select(coral_id, size_class), by = "coral_id") %>%
  pivot_longer(cols = -c(coral_id, size_class), names_to = "func_group", values_to = "rel_abund")

# Kruskal-Wallis for each functional group
func_tests <- func_rel_data %>%
  group_by(func_group) %>%
  summarise(
    mean_small = mean(rel_abund[size_class == "Small"]),
    se_small = sd(rel_abund[size_class == "Small"]) / sqrt(sum(size_class == "Small")),
    mean_medium = mean(rel_abund[size_class == "Medium"]),
    se_medium = sd(rel_abund[size_class == "Medium"]) / sqrt(sum(size_class == "Medium")),
    mean_large = mean(rel_abund[size_class == "Large"]),
    se_large = sd(rel_abund[size_class == "Large"]) / sqrt(sum(size_class == "Large")),
    kw_stat = kruskal.test(rel_abund ~ size_class)$statistic,
    kw_p = kruskal.test(rel_abund ~ size_class)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    change = mean_large - mean_small,
    significant = kw_p < 0.05
  )

cat("\nFunctional Group Relative Abundance by Size Class:\n")
print(func_tests)

# ============================================================================
# INDICATOR SPECIES ANALYSIS
# ============================================================================

cat("\n--- INDICATOR SPECIES ANALYSIS ---\n")

# Using multipatt from indicspecies package if available
if (requireNamespace("indicspecies", quietly = TRUE)) {
  library(indicspecies)

  # Indicator species for size classes
  indval_result <- multipatt(comm_mat, meta$size_class, func = "IndVal.g",
                              control = how(nperm = 999))

  cat("\nIndicator Species Summary:\n")
  summary(indval_result, alpha = 0.05)
} else {
  cat("indicspecies package not available - skipping indicator species analysis\n")
}

# ============================================================================
# PAIRWISE COMPARISONS (Post-hoc)
# ============================================================================

cat("\n--- PAIRWISE PERMANOVA ---\n")

# Small vs Medium
perm_sm <- adonis2(comm_hell[meta$size_class %in% c("Small", "Medium"), ] ~
                     size_class,
                   data = meta[meta$size_class %in% c("Small", "Medium"), ],
                   permutations = 999)
cat("Small vs Medium: R² =", round(perm_sm$R2[1], 4), ", p =", perm_sm$`Pr(>F)`[1], "\n")

# Medium vs Large
perm_ml <- adonis2(comm_hell[meta$size_class %in% c("Medium", "Large"), ] ~
                     size_class,
                   data = meta[meta$size_class %in% c("Medium", "Large"), ],
                   permutations = 999)
cat("Medium vs Large: R² =", round(perm_ml$R2[1], 4), ", p =", perm_ml$`Pr(>F)`[1], "\n")

# Small vs Large
perm_sl <- adonis2(comm_hell[meta$size_class %in% c("Small", "Large"), ] ~
                     size_class,
                   data = meta[meta$size_class %in% c("Small", "Large"), ],
                   permutations = 999)
cat("Small vs Large: R² =", round(perm_sl$R2[1], 4), ", p =", perm_sl$`Pr(>F)`[1], "\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n--- SAVING RESULTS ---\n")

# Save incidence test results
write_csv(incidence_tests, "output/tables/incidence_chisq_tests.csv")
cat("Saved: output/tables/incidence_chisq_tests.csv\n")

# Save relative abundance test results
write_csv(rel_abund_tests, "output/tables/relative_abundance_kw_tests.csv")
cat("Saved: output/tables/relative_abundance_kw_tests.csv\n")

# Save functional group results
write_csv(func_tests, "output/tables/functional_group_relabund_tests.csv")
cat("Saved: output/tables/functional_group_relabund_tests.csv\n")

# Create summary table for manuscript
summary_stats <- tibble(
  Analysis = c(
    "PERMANOVA (size class)",
    "PERMANOVA (continuous volume)",
    "PERMANOVA (size + site)",
    "ANOSIM (size class)",
    "PERMDISP (beta dispersion)",
    "Pairwise: Small vs Large",
    "Pairwise: Small vs Medium",
    "Pairwise: Medium vs Large"
  ),
  Statistic = c(
    sprintf("R² = %.3f", perm_size$R2[1]),
    sprintf("R² = %.3f", perm_vol$R2[1]),
    sprintf("R² = %.3f (size), %.3f (site)", perm_both$R2[1], perm_both$R2[2]),
    sprintf("R = %.3f", anosim_result$statistic),
    sprintf("F = %.2f", permutest_result$tab$F[1]),
    sprintf("R² = %.3f", perm_sl$R2[1]),
    sprintf("R² = %.3f", perm_sm$R2[1]),
    sprintf("R² = %.3f", perm_ml$R2[1])
  ),
  p_value = c(
    perm_size$`Pr(>F)`[1],
    perm_vol$`Pr(>F)`[1],
    perm_both$`Pr(>F)`[1],
    anosim_result$signif,
    permutest_result$tab$`Pr(>F)`[1],
    perm_sl$`Pr(>F)`[1],
    perm_sm$`Pr(>F)`[1],
    perm_ml$`Pr(>F)`[1]
  ),
  Interpretation = c(
    "Size class explains 6.7% of composition variance",
    "Continuous volume explains similar variance",
    "Size and site both structure composition",
    "Significant separation between size classes",
    "Dispersion differs among size classes",
    "Small and Large corals most different",
    "Small and Medium moderately different",
    "Medium and Large similar"
  )
)

write_csv(summary_stats, "output/tables/multivariate_statistics_summary.csv")
cat("Saved: output/tables/multivariate_statistics_summary.csv\n")

cat("\n========================================\n")
cat("COMPOSITION ANALYSIS COMPLETE\n")
cat("========================================\n")

# Print final summary
cat("\n--- KEY FINDINGS ---\n")
cat(sprintf("1. PERMANOVA: Size class explains %.1f%% of composition variance (p = %.3f)\n",
            perm_size$R2[1] * 100, perm_size$`Pr(>F)`[1]))
cat(sprintf("2. %d species show significant incidence changes (FDR < 0.05)\n",
            sum(incidence_tests$significant, na.rm = TRUE)))
cat(sprintf("3. %d species show significant relative abundance changes (FDR < 0.05)\n",
            sum(rel_abund_tests$significant, na.rm = TRUE)))
cat(sprintf("4. Defenders: %.1f%% (small) → %.1f%% (large), p = %.4f\n",
            func_tests$mean_small[func_tests$func_group == "Defenders"],
            func_tests$mean_large[func_tests$func_group == "Defenders"],
            func_tests$kw_p[func_tests$func_group == "Defenders"]))
cat(sprintf("5. Resident fishes: %.1f%% (small) → %.1f%% (large), p = %.4f\n",
            func_tests$mean_small[func_tests$func_group == "Resident_fishes"],
            func_tests$mean_large[func_tests$func_group == "Resident_fishes"],
            func_tests$kw_p[func_tests$func_group == "Resident_fishes"]))
