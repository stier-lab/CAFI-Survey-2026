#!/usr/bin/env Rscript
# ============================================================================
# 99_extract_all_statistics.R - Extract All Statistical Results for Manuscript
#
# Purpose: Run all analyses and compile a comprehensive table of all
# statistical results including:
#   - Test name and hypothesis
#   - Degrees of freedom
#   - Test statistic (F, t, chi-sq, etc.)
#   - P-value
#   - Effect size (R², Cohen's d, eta², etc.)
#   - 95% CI where applicable
#
# This table serves as the master reference for manuscript statistics
# to ensure full reproducibility.
#
# Author: CAFI Analysis Pipeline
# Date: 2025-11-23
# ============================================================================

cat("\n========================================\n")
cat("Extracting All Statistical Results\n")
cat("========================================\n\n")

# Load libraries
source(here::here("scripts/00_load_libraries.R"))

# Initialize results table
stats_results <- data.frame(

  analysis = character(),
  hypothesis = character(),
  test_type = character(),
  predictor = character(),
  response = character(),
  df1 = numeric(),
  df2 = numeric(),
  statistic = numeric(),
  statistic_type = character(),
  p_value = numeric(),
  effect_size = numeric(),

effect_type = character(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  n = integer(),
  notes = character(),
  stringsAsFactors = FALSE
)

# Helper function to add results
add_result <- function(analysis, hypothesis, test_type, predictor, response,
                       df1 = NA, df2 = NA, statistic = NA, statistic_type = "",
                       p_value = NA, effect_size = NA, effect_type = "",
                       ci_lower = NA, ci_upper = NA, n = NA, notes = "") {
  new_row <- data.frame(
    analysis = analysis,
    hypothesis = hypothesis,
    test_type = test_type,
    predictor = predictor,
    response = response,
    df1 = df1,
    df2 = df2,
    statistic = statistic,
    statistic_type = statistic_type,
    p_value = p_value,
    effect_size = effect_size,
    effect_type = effect_type,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    n = n,
    notes = notes,
    stringsAsFactors = FALSE
  )
  stats_results <<- rbind(stats_results, new_row)
}

# ============================================================================
# Load all processed data
# ============================================================================

cat("Loading processed data...\n")

survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))

# Load condition scores if available
condition_scores <- tryCatch(
  readRDS(file.path(SURVEY_OBJECTS, "coral_condition_scores.rds")),
  error = function(e) NULL
)

cat("  Data loaded successfully\n\n")

# ============================================================================
# 1. COMMUNITY COMPOSITION ANALYSIS (H1)
# ============================================================================

cat("1. Analyzing community composition (H1)...\n")

# Prepare community matrix
community_matrix_clean <- community_matrix[rowSums(community_matrix, na.rm = TRUE) > 0, ]
community_matrix_clean[is.na(community_matrix_clean)] <- 0

# Align metadata
metadata_aligned <- metadata %>%
  filter(coral_id %in% rownames(community_matrix_clean)) %>%
  arrange(match(coral_id, rownames(community_matrix_clean)))

# Beta diversity
beta_bray <- vegdist(community_matrix_clean, method = "bray")

# PERMANOVA - Site effect
set.seed(123)
permanova_site <- adonis2(beta_bray ~ site,
                          data = metadata_aligned,
                          permutations = 999)

add_result(
  analysis = "Community Composition",
  hypothesis = "H1",
  test_type = "PERMANOVA",
  predictor = "Site",
  response = "Bray-Curtis dissimilarity",
  df1 = permanova_site$Df[1],
  df2 = permanova_site$Df[3],
  statistic = permanova_site$F[1],
  statistic_type = "pseudo-F",
  p_value = permanova_site$`Pr(>F)`[1],
  effect_size = permanova_site$R2[1],
  effect_type = "R2",
  n = nrow(community_matrix_clean),
  notes = "Bray-Curtis distance, 999 permutations"
)

# PERMANOVA - Morphotype effect (filter NAs first)
morph_complete <- !is.na(metadata_aligned$morphotype)
if (sum(morph_complete) > 10) {
  beta_bray_morph <- vegdist(community_matrix_clean[morph_complete, ], method = "bray")
  permanova_morph <- adonis2(beta_bray_morph ~ morphotype,
                             data = metadata_aligned[morph_complete, ],
                             permutations = 999)

  add_result(
    analysis = "Community Composition",
    hypothesis = "H1",
    test_type = "PERMANOVA",
    predictor = "Branch architecture",
    response = "Bray-Curtis dissimilarity",
    df1 = permanova_morph$Df[1],
    df2 = permanova_morph$Df[3],
    statistic = permanova_morph$F[1],
    statistic_type = "pseudo-F",
    p_value = permanova_morph$`Pr(>F)`[1],
    effect_size = permanova_morph$R2[1],
    effect_type = "R2",
    n = sum(morph_complete),
    notes = "Bray-Curtis distance, 999 permutations, NAs excluded"
  )
}

# PERMDISP - Homogeneity of dispersions
betadisp_site <- betadisper(beta_bray, metadata_aligned$site)
permdisp_result <- permutest(betadisp_site, permutations = 999)

add_result(
  analysis = "Community Composition",
  hypothesis = "H1",
  test_type = "PERMDISP",
  predictor = "Site",
  response = "Multivariate dispersion",
  df1 = permdisp_result$tab$Df[1],
  df2 = permdisp_result$tab$Df[2],
  statistic = permdisp_result$tab$F[1],
  statistic_type = "F",
  p_value = permdisp_result$tab$`Pr(>F)`[1],
  n = nrow(community_matrix_clean),
  notes = "Tests homogeneity of dispersions"
)

# NMDS stress
set.seed(123)
nmds_bray <- metaMDS(community_matrix_clean, distance = "bray", k = 2,
                     trymax = 100, autotransform = FALSE)

add_result(
  analysis = "Community Composition",
  hypothesis = "H1",
  test_type = "NMDS",
  predictor = "Community matrix",
  response = "2D ordination",
  statistic = nmds_bray$stress,
  statistic_type = "stress",
  n = nrow(community_matrix_clean),
  notes = paste("Stress <0.1 excellent, <0.2 good;",
                ifelse(nmds_bray$stress < 0.1, "excellent",
                       ifelse(nmds_bray$stress < 0.2, "good", "acceptable")))
)

cat("  PERMANOVA site: F =", round(permanova_site$F[1], 2),
    ", p =", format.pval(permanova_site$`Pr(>F)`[1], digits = 3), "\n")

# ============================================================================
# 2. DIVERSITY ANALYSIS (H4)
# ============================================================================

cat("2. Analyzing diversity patterns (H4)...\n")

# Alpha diversity metrics
alpha_diversity <- data.frame(
  coral_id = rownames(community_matrix),
  species_richness = specnumber(community_matrix),
  shannon = vegan::diversity(community_matrix, index = "shannon"),
  simpson = vegan::diversity(community_matrix, index = "simpson")
) %>%
  left_join(metadata, by = "coral_id")

# Site comparison - Kruskal-Wallis for richness
kw_richness <- kruskal.test(species_richness ~ site, data = alpha_diversity)

add_result(
  analysis = "Alpha Diversity",
  hypothesis = "H4",
  test_type = "Kruskal-Wallis",
  predictor = "Site",
  response = "Species richness",
  df1 = kw_richness$parameter,
  statistic = kw_richness$statistic,
  statistic_type = "chi-squared",
  p_value = kw_richness$p.value,
  n = nrow(alpha_diversity),
  notes = "Non-parametric test for site differences"
)

# Site comparison - Kruskal-Wallis for Shannon
kw_shannon <- kruskal.test(shannon ~ site, data = alpha_diversity)

add_result(
  analysis = "Alpha Diversity",
  hypothesis = "H4",
  test_type = "Kruskal-Wallis",
  predictor = "Site",
  response = "Shannon diversity",
  df1 = kw_shannon$parameter,
  statistic = kw_shannon$statistic,
  statistic_type = "chi-squared",
  p_value = kw_shannon$p.value,
  n = nrow(alpha_diversity),
  notes = "Non-parametric test for site differences"
)

# Gamma diversity
gamma_div <- vegan::diversity(colSums(community_matrix), index = "shannon")

add_result(
  analysis = "Regional Diversity",
  hypothesis = "H4",
  test_type = "Shannon Index",
  predictor = "All sites pooled",
  response = "Gamma diversity",
  effect_size = gamma_div,
  effect_type = "H'",
  n = sum(community_matrix),
  notes = "Total regional diversity"
)

cat("  Gamma diversity H' =", round(gamma_div, 2), "\n")

# ============================================================================
# 3. CORAL-CAFI SCALING (H2)
# ============================================================================

cat("3. Analyzing coral-CAFI scaling (H2)...\n")

# Prepare coral-level data
coral_cafi <- cafi_clean %>%
  group_by(coral_id) %>%
  summarise(
    total_cafi = n(),
    otu_richness = n_distinct(species),
    .groups = "drop"
  ) %>%
  right_join(metadata, by = "coral_id") %>%
  left_join(
    survey_master %>% select(coral_id, volume_lab),
    by = "coral_id"
  ) %>%
  mutate(
    total_cafi = replace_na(total_cafi, 0),
    otu_richness = replace_na(otu_richness, 0),
    log_volume = log10(volume_lab + 1),
    log_cafi = log10(total_cafi + 1)
  ) %>%
  filter(!is.na(volume_lab), volume_lab > 0)

# Power-law scaling (log-log regression)
scaling_model <- lm(log_cafi ~ log_volume, data = coral_cafi)
scaling_summary <- summary(scaling_model)
scaling_confint <- confint(scaling_model)

add_result(
  analysis = "Volume-Abundance Scaling",
  hypothesis = "H2",
  test_type = "Linear regression (log-log)",
  predictor = "log10(Coral volume)",
  response = "log10(CAFI abundance)",
  df1 = 1,
  df2 = scaling_summary$df[2],
  statistic = scaling_summary$fstatistic[1],
  statistic_type = "F",
  p_value = pf(scaling_summary$fstatistic[1],
               scaling_summary$fstatistic[2],
               scaling_summary$fstatistic[3],
               lower.tail = FALSE),
  effect_size = scaling_summary$r.squared,
  effect_type = "R2",
  ci_lower = scaling_confint[2, 1],
  ci_upper = scaling_confint[2, 2],
  n = nrow(coral_cafi),
  notes = paste("Scaling exponent =", round(coef(scaling_model)[2], 3),
                "; expect <1 for propagule redirection")
)

# Extract slope (scaling exponent)
add_result(
  analysis = "Scaling Exponent",
  hypothesis = "H2",
  test_type = "Regression coefficient",
  predictor = "log10(Coral volume)",
  response = "log10(CAFI abundance)",
  df1 = scaling_summary$df[2],
  statistic = scaling_summary$coefficients[2, 3],
  statistic_type = "t",
  p_value = scaling_summary$coefficients[2, 4],
  effect_size = coef(scaling_model)[2],
  effect_type = "slope (b)",
  ci_lower = scaling_confint[2, 1],
  ci_upper = scaling_confint[2, 2],
  n = nrow(coral_cafi),
  notes = "Power-law exponent; <1 indicates density decrease with size"
)

cat("  Scaling exponent =", round(coef(scaling_model)[2], 3),
    ", R2 =", round(scaling_summary$r.squared, 3), "\n")

# ============================================================================
# 4. CONDITION-DIVERSITY RELATIONSHIPS (H4)
# ============================================================================

cat("4. Analyzing condition-diversity relationships (H4)...\n")

if (!is.null(condition_scores)) {
  # Merge condition with CAFI data
  cafi_condition <- coral_cafi %>%
    inner_join(condition_scores %>% select(coral_id, condition_score),
               by = "coral_id") %>%
    filter(!is.na(condition_score))

  # CAFI abundance ~ Condition
  cond_abund_model <- lm(total_cafi ~ condition_score, data = cafi_condition)
  cond_abund_summary <- summary(cond_abund_model)
  cond_abund_confint <- confint(cond_abund_model)

  add_result(
    analysis = "Condition-Abundance",
    hypothesis = "H4",
    test_type = "Linear regression",
    predictor = "Coral condition score",
    response = "CAFI abundance",
    df1 = 1,
    df2 = cond_abund_summary$df[2],
    statistic = cond_abund_summary$fstatistic[1],
    statistic_type = "F",
    p_value = pf(cond_abund_summary$fstatistic[1],
                 cond_abund_summary$fstatistic[2],
                 cond_abund_summary$fstatistic[3],
                 lower.tail = FALSE),
    effect_size = cond_abund_summary$r.squared,
    effect_type = "R2",
    ci_lower = cond_abund_confint[2, 1],
    ci_upper = cond_abund_confint[2, 2],
    n = nrow(cafi_condition),
    notes = "Position-corrected condition score"
  )

  # OTU richness ~ Condition
  cond_rich_model <- lm(otu_richness ~ condition_score, data = cafi_condition)
  cond_rich_summary <- summary(cond_rich_model)
  cond_rich_confint <- confint(cond_rich_model)

  add_result(
    analysis = "Condition-Richness",
    hypothesis = "H4",
    test_type = "Linear regression",
    predictor = "Coral condition score",
    response = "OTU richness",
    df1 = 1,
    df2 = cond_rich_summary$df[2],
    statistic = cond_rich_summary$fstatistic[1],
    statistic_type = "F",
    p_value = pf(cond_rich_summary$fstatistic[1],
                 cond_rich_summary$fstatistic[2],
                 cond_rich_summary$fstatistic[3],
                 lower.tail = FALSE),
    effect_size = cond_rich_summary$r.squared,
    effect_type = "R2",
    ci_lower = cond_rich_confint[2, 1],
    ci_upper = cond_rich_confint[2, 2],
    n = nrow(cafi_condition),
    notes = "Position-corrected condition score"
  )

  cat("  Condition-Abundance: R2 =", round(cond_abund_summary$r.squared, 3),
      ", p =", format.pval(pf(cond_abund_summary$fstatistic[1],
                              cond_abund_summary$fstatistic[2],
                              cond_abund_summary$fstatistic[3],
                              lower.tail = FALSE), digits = 3), "\n")
}

# ============================================================================
# 5. LOCAL NEIGHBORHOOD EFFECTS (H3)
# ============================================================================

cat("5. Analyzing local neighborhood effects (H3)...\n")

# Check if neighborhood data exists
if ("mean_live_volume_of_neighbors" %in% names(survey_master)) {

  # Prepare neighborhood data
  neighborhood_data <- coral_cafi %>%
    left_join(
      survey_master %>%
        select(coral_id,
               starts_with("mean_"),
               starts_with("total_"),
               starts_with("n_")),
      by = "coral_id"
    ) %>%
    filter(!is.na(volume_lab))

  # Count neighbor columns
  neighbor_cols <- names(neighborhood_data)[grepl("neighbor", names(neighborhood_data), ignore.case = TRUE)]

  if (length(neighbor_cols) > 0) {
    # Example: Effect of neighbor count on CAFI (if available)
    if ("n_neighbors" %in% names(neighborhood_data)) {
      neighbor_model <- glm(total_cafi ~ n_neighbors,
                            data = neighborhood_data,
                            family = "poisson")
      neighbor_summary <- summary(neighbor_model)

      # Calculate pseudo-R2 (McFadden's)
      null_model <- glm(total_cafi ~ 1, data = neighborhood_data, family = "poisson")
      pseudo_r2 <- 1 - (logLik(neighbor_model) / logLik(null_model))

      add_result(
        analysis = "Neighbor Count Effect",
        hypothesis = "H3",
        test_type = "GLM (Poisson)",
        predictor = "Number of neighbors",
        response = "CAFI abundance",
        df1 = 1,
        df2 = neighbor_summary$df.residual,
        statistic = neighbor_summary$coefficients[2, 3],
        statistic_type = "z",
        p_value = neighbor_summary$coefficients[2, 4],
        effect_size = as.numeric(pseudo_r2),
        effect_type = "pseudo-R2",
        n = nrow(neighborhood_data),
        notes = "McFadden's pseudo-R2"
      )
    }
  }
}

# ============================================================================
# 6. INDICATOR SPECIES ANALYSIS
# ============================================================================

cat("6. Running indicator species analysis...\n")

# Indicator species for sites
indval_site <- multipatt(community_matrix_clean,
                         metadata_aligned$site,
                         control = how(nperm = 999))

# Count significant indicators
sig_indicators <- sum(indval_site$sign$p.value < 0.05, na.rm = TRUE)

add_result(
  analysis = "Indicator Species",
  hypothesis = "H1",
  test_type = "IndVal",
  predictor = "Site",
  response = "Species associations",
  statistic = sig_indicators,
  statistic_type = "n significant",
  n = ncol(community_matrix_clean),
  notes = paste(sig_indicators, "of", ncol(community_matrix_clean),
                "species are significant indicators (p < 0.05)")
)

cat("  Found", sig_indicators, "significant indicator species\n")

# ============================================================================
# 7. BRANCH ARCHITECTURE EFFECTS
# ============================================================================

cat("7. Analyzing branch architecture effects...\n")

# Skip branch architecture analysis - morphotype has >2 levels
# Would need Kruskal-Wallis instead of Wilcoxon
# This can be added later if needed
cat("  Skipping - morphotype has multiple levels, use Kruskal-Wallis instead\n")

if (FALSE) {  # Disabled for now

  if (nrow(branch_data) > 10) {
    # Wilcoxon test for CAFI abundance by branch width
    wilcox_branch <- wilcox.test(total_cafi ~ branch_width, data = branch_data)

    # Calculate effect size (r = Z / sqrt(N))
    z_stat <- qnorm(wilcox_branch$p.value / 2)
    effect_r <- abs(z_stat) / sqrt(nrow(branch_data))

    add_result(
      analysis = "Branch Architecture",
      hypothesis = "H2",
      test_type = "Wilcoxon rank-sum",
      predictor = "Branch width (tight vs wide)",
      response = "CAFI abundance",
      statistic = wilcox_branch$statistic,
      statistic_type = "W",
      p_value = wilcox_branch$p.value,
      effect_size = effect_r,
      effect_type = "r",
      n = nrow(branch_data),
      notes = "Non-parametric comparison of branch architectures"
    )

    # Wilcoxon test for richness
    wilcox_rich <- wilcox.test(otu_richness ~ branch_width, data = branch_data)
    z_rich <- qnorm(wilcox_rich$p.value / 2)
    effect_r_rich <- abs(z_rich) / sqrt(nrow(branch_data))

    add_result(
      analysis = "Branch Architecture",
      hypothesis = "H2",
      test_type = "Wilcoxon rank-sum",
      predictor = "Branch width (tight vs wide)",
      response = "OTU richness",
      statistic = wilcox_rich$statistic,
      statistic_type = "W",
      p_value = wilcox_rich$p.value,
      effect_size = effect_r_rich,
      effect_type = "r",
      n = nrow(branch_data),
      notes = "Non-parametric comparison of branch architectures"
    )

    cat("  Branch effect on abundance: W =", wilcox_branch$statistic,
        ", p =", format.pval(wilcox_branch$p.value, digits = 3), "\n")
  }
}

# ============================================================================
# 8. DESCRIPTIVE STATISTICS
# ============================================================================

cat("8. Calculating descriptive statistics...\n")

# Overall summary statistics
add_result(
  analysis = "Descriptive",
  hypothesis = "-",
  test_type = "Summary",
  predictor = "-",
  response = "Total corals surveyed",
  statistic = n_distinct(cafi_clean$coral_id),
  statistic_type = "n",
  notes = "Number of coral colonies sampled"
)

add_result(
  analysis = "Descriptive",
  hypothesis = "-",
  test_type = "Summary",
  predictor = "-",
  response = "Total CAFI observed",
  statistic = nrow(cafi_clean),
  statistic_type = "n",
  notes = "Total individual CAFI counted"
)

add_result(
  analysis = "Descriptive",
  hypothesis = "-",
  test_type = "Summary",
  predictor = "-",
  response = "OTU richness (total)",
  statistic = n_distinct(cafi_clean$species),
  statistic_type = "n",
  notes = "Total number of morphological OTUs"
)

# Mean CAFI per coral
cafi_per_coral <- cafi_clean %>%
  group_by(coral_id) %>%
  summarise(n = n(), .groups = "drop")

add_result(
  analysis = "Descriptive",
  hypothesis = "-",
  test_type = "Mean +/- SD",
  predictor = "-",
  response = "CAFI per coral",
  effect_size = mean(cafi_per_coral$n),
  effect_type = "mean",
  ci_lower = mean(cafi_per_coral$n) - sd(cafi_per_coral$n),
  ci_upper = mean(cafi_per_coral$n) + sd(cafi_per_coral$n),
  n = nrow(cafi_per_coral),
  notes = paste("SD =", round(sd(cafi_per_coral$n), 1))
)

# Mean richness per coral
richness_per_coral <- cafi_clean %>%
  group_by(coral_id) %>%
  summarise(richness = n_distinct(species), .groups = "drop")

add_result(
  analysis = "Descriptive",
  hypothesis = "-",
  test_type = "Mean +/- SD",
  predictor = "-",
  response = "OTU richness per coral",
  effect_size = mean(richness_per_coral$richness),
  effect_type = "mean",
  ci_lower = mean(richness_per_coral$richness) - sd(richness_per_coral$richness),
  ci_upper = mean(richness_per_coral$richness) + sd(richness_per_coral$richness),
  n = nrow(richness_per_coral),
  notes = paste("SD =", round(sd(richness_per_coral$richness), 1))
)

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\nSaving statistical results...\n")

# Format p-values for display
stats_results <- stats_results %>%
  mutate(
    p_formatted = case_when(
      is.na(p_value) ~ "-",
      p_value < 0.001 ~ "<0.001",
      p_value < 0.01 ~ sprintf("%.3f", p_value),
      TRUE ~ sprintf("%.3f", p_value)
    ),
    effect_formatted = case_when(
      is.na(effect_size) ~ "-",
      TRUE ~ sprintf("%.3f", effect_size)
    ),
    ci_formatted = case_when(
      is.na(ci_lower) | is.na(ci_upper) ~ "-",
      TRUE ~ sprintf("[%.3f, %.3f]", ci_lower, ci_upper)
    )
  )

# Save as CSV
write_csv(stats_results,
          file.path(SURVEY_TABLES, "manuscript_statistics_reference.csv"))

# Save as formatted table for easy reference
stats_display <- stats_results %>%
  select(analysis, hypothesis, test_type, predictor, response,
         statistic_type, statistic, p_formatted,
         effect_type, effect_formatted, ci_formatted, n, notes) %>%
  rename(
    Analysis = analysis,
    Hypothesis = hypothesis,
    Test = test_type,
    Predictor = predictor,
    Response = response,
    `Stat Type` = statistic_type,
    Statistic = statistic,
    `P-value` = p_formatted,
    `Effect Type` = effect_type,
    `Effect Size` = effect_formatted,
    `95% CI` = ci_formatted,
    N = n,
    Notes = notes
  )

write_csv(stats_display,
          file.path(SURVEY_TABLES, "manuscript_statistics_display.csv"))

# Also save as RDS for programmatic access
saveRDS(stats_results,
        file.path(SURVEY_OBJECTS, "manuscript_statistics.rds"))

cat("\n========================================\n")
cat("Statistical Results Summary\n")
cat("========================================\n\n")

cat("Total statistical tests recorded:", nrow(stats_results), "\n\n")

cat("Results by hypothesis:\n")
print(table(stats_results$hypothesis))

cat("\nResults by test type:\n")
print(table(stats_results$test_type))

cat("\n\nKey findings:\n")

# PERMANOVA site
perm_site <- stats_results %>%
  filter(analysis == "Community Composition", predictor == "Site", test_type == "PERMANOVA")
if (nrow(perm_site) > 0) {
  cat("  H1 - Community composition differs by site:\n")
  cat("       PERMANOVA pseudo-F =", round(perm_site$statistic, 2),
      ", R2 =", round(perm_site$effect_size, 3),
      ", p =", perm_site$p_formatted, "\n")
}

# Scaling exponent
scaling <- stats_results %>%
  filter(analysis == "Scaling Exponent")
if (nrow(scaling) > 0) {
  cat("  H2 - Abundance-volume scaling:\n")
  cat("       Exponent b =", round(scaling$effect_size, 3),
      " [", round(scaling$ci_lower, 3), ",", round(scaling$ci_upper, 3), "]",
      ", p =", scaling$p_formatted, "\n")
}

# Condition relationships
cond_abund <- stats_results %>%
  filter(analysis == "Condition-Abundance")
if (nrow(cond_abund) > 0) {
  cat("  H4 - Condition-abundance relationship:\n")
  cat("       R2 =", round(cond_abund$effect_size, 3),
      ", p =", cond_abund$p_formatted, "\n")
}

cat("\nFiles saved:\n")
cat("  - Full results:", file.path(SURVEY_TABLES, "manuscript_statistics_reference.csv"), "\n")
cat("  - Display table:", file.path(SURVEY_TABLES, "manuscript_statistics_display.csv"), "\n")
cat("  - R object:", file.path(SURVEY_OBJECTS, "manuscript_statistics.rds"), "\n")

cat("\n Use this table as the master reference for all manuscript statistics.\n")
cat("========================================\n")
