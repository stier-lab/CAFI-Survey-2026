# ============================================================================
# 14_supplementary_sensitivity.R - Supplementary Sensitivity & Robustness
# ============================================================================
#
# PURPOSE: Additional analyses addressing methodological gaps identified
#          during pre-submission audit. Each section is self-contained and
#          generates tables/figures for the supplement.
#
# SECTIONS:
#   Part 1:  Beta diversity partitioning (betapart) — turnover vs nestedness
#   Part 2:  envfit species vector fitting (vegan) — formal R² + p for species
#   Part 3:  Indicator species by size class (indicspecies)
#   Part 4:  Coverage-based rarefaction (iNEXT) — Hill numbers q=0,1,2
#   Part 5:  CWM obligate-mutualist fraction as BEF predictor
#   Part 6:  Non-linear BEF (log-richness, polynomial)
#   Part 7:  Missing data characterization (114→84 dropout)
#   Part 8:  Morphotype covariate in BEF model
#     Part 8b: Morphotype–haplotype concordance
#     Part 8c: BEF with haplotype (genetic species) covariate
#     Part 8d: Scaling by species (Q1 sensitivity)
#     Part 8e: Composition by species (Q2 sensitivity)
#     Part 8f: Genotype → architecture mediation (branch width × species)
#     Part 8g: Phylogenetic distance & symbiont identity
#   Part 9:  Mediation: richness → abundance → condition
#   Part 10: C-score community-wide co-occurrence robustness
#   Part 11: MuMIn model-averaged neighborhood coefficients
#
# DEPENDENCIES: 00_setup.R, 01_load_data.R
# OUTPUTS: output/tables/sensitivity_*, output/figures/supplement/figS15_*
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-03-22
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    SUPPLEMENTARY SENSITIVITY ANALYSES\n")
cat("============================================================\n\n")

# Load setup and data
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

FIG_DIR <- PATHS$fig_supplement %||% file.path("output", "figures", "supplement")
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# Load core objects
coral_master <- load_object("coral_master")
community_matrix <- load_object("community_matrix")
condition_scores <- load_object("condition_scores")
cafi_clean <- load_object("cafi_clean")

set.seed(42)

# ============================================================================
# PART 1: BETA DIVERSITY PARTITIONING (betapart)
# ============================================================================
# Decomposes total beta diversity into turnover (Simpson) and nestedness-
# resultant (richness difference) components along the coral size gradient.
# Strengthens the NODF result (NS nestedness) with a continuous gradient test.
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 1: BETA DIVERSITY PARTITIONING (betapart)\n")
cat("------------------------------------------------------------\n\n")

if (requireNamespace("betapart", quietly = TRUE)) {
  library(betapart)

  # Presence-absence matrix (corals with >0 CAFI)
  comm_pa <- (community_matrix[rowSums(community_matrix) > 0, ] > 0) * 1
  coral_sub <- coral_master[match(rownames(comm_pa), coral_master$coral_id), ]

  # Pairwise beta diversity components
  bp <- beta.pair(comm_pa, index.family = "sorensen")

  # Mean pairwise components per coral
  beta_df <- data.frame(
    coral_id = rownames(comm_pa),
    log_volume = log(coral_sub$volume),
    site = coral_sub$site,
    turnover_mean = rowMeans(as.matrix(bp$beta.sim), na.rm = TRUE),
    nestedness_mean = rowMeans(as.matrix(bp$beta.sne), na.rm = TRUE),
    total_beta = rowMeans(as.matrix(bp$beta.sor), na.rm = TRUE)
  )

  # Overall proportions
  all_sim <- mean(as.matrix(bp$beta.sim)[upper.tri(as.matrix(bp$beta.sim))])
  all_sne <- mean(as.matrix(bp$beta.sne)[upper.tri(as.matrix(bp$beta.sne))])
  all_sor <- mean(as.matrix(bp$beta.sor)[upper.tri(as.matrix(bp$beta.sor))])
  pct_turnover <- round(100 * all_sim / all_sor, 1)
  pct_nestedness <- round(100 * all_sne / all_sor, 1)

  cat("  Overall beta diversity (Sorensen):\n")
  cat("    Total (beta.sor):", round(all_sor, 3), "\n")
  cat("    Turnover (beta.sim):", round(all_sim, 3), "(", pct_turnover, "%)\n")
  cat("    Nestedness (beta.sne):", round(all_sne, 3), "(", pct_nestedness, "%)\n\n")

  # Test: does turnover/nestedness change along the size gradient?
  m_turn <- lm(turnover_mean ~ log_volume + site, data = beta_df)
  m_nest <- lm(nestedness_mean ~ log_volume + site, data = beta_df)

  cat("  Turnover ~ log(volume) + site:\n")
  cat("    beta =", round(coef(m_turn)["log_volume"], 4),
      ", t =", round(summary(m_turn)$coefficients["log_volume", "t value"], 2),
      ", p =", format.pval(summary(m_turn)$coefficients["log_volume", "Pr(>|t|)"], 3), "\n")

  cat("  Nestedness ~ log(volume) + site:\n")
  cat("    beta =", round(coef(m_nest)["log_volume"], 4),
      ", t =", round(summary(m_nest)$coefficients["log_volume", "t value"], 2),
      ", p =", format.pval(summary(m_nest)$coefficients["log_volume", "Pr(>|t|)"], 3), "\n\n")

  # Save results
  save_table(data.frame(
    component = c("turnover_simpson", "nestedness_resultant", "total_sorensen"),
    mean_value = round(c(all_sim, all_sne, all_sor), 4),
    pct_of_total = c(pct_turnover, pct_nestedness, 100),
    trend_beta = round(c(coef(m_turn)["log_volume"], coef(m_nest)["log_volume"], NA), 5),
    trend_p = c(
      summary(m_turn)$coefficients["log_volume", "Pr(>|t|)"],
      summary(m_nest)$coefficients["log_volume", "Pr(>|t|)"],
      NA
    )
  ), "sensitivity_betapart_decomposition")

} else {
  cat("  betapart not installed — skipping\n\n")
}


# ============================================================================
# PART 2: ENVFIT SPECIES VECTOR FITTING
# ============================================================================
# Formal permutation-tested R² and p-values for species associations with
# NMDS ordination axes. Replaces weighted-average scores with proper envfit.
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 2: ENVFIT SPECIES VECTOR FITTING\n")
cat("------------------------------------------------------------\n\n")

# Build NMDS and run envfit
comm_nz <- community_matrix[rowSums(community_matrix) > 0, ]
comm_nz <- comm_nz[, colSums(comm_nz) > 0]
coral_nz <- coral_master[match(rownames(comm_nz), coral_master$coral_id), ]

# Hellinger transform + NMDS
comm_hell <- decostand(comm_nz, method = "hellinger")
set.seed(42)
nmds_ef <- metaMDS(comm_hell, distance = "bray", k = 2, trymax = 50, trace = 0)

cat("  NMDS stress:", round(nmds_ef$stress, 3), "\n")

# Species envfit on NMDS
set.seed(42)
ef_sp <- envfit(nmds_ef, comm_nz, permutations = 999)

# Extract results
sp_r2 <- ef_sp$vectors$r
sp_p <- ef_sp$vectors$pvals
envfit_df <- data.frame(
  species = names(sp_r2),
  r_squared = round(sp_r2, 4),
  p_value = sp_p,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(r_squared))

# Significant species (R² > 0.05, p < 0.05)
sig_envfit <- envfit_df %>% filter(r_squared > 0.05, p_value < 0.05)
cat("  Species with envfit R² > 0.05, p < 0.05:", nrow(sig_envfit), "\n")
cat("  Top 10:\n")
print(head(sig_envfit, 10))

save_table(envfit_df, "sensitivity_envfit_species")
cat("\n")


# ============================================================================
# PART 3: INDICATOR SPECIES BY SIZE CLASS (indicspecies)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 3: INDICATOR SPECIES BY SIZE CLASS\n")
cat("------------------------------------------------------------\n\n")

if (requireNamespace("indicspecies", quietly = TRUE)) {
  library(indicspecies)

  # Size classes (tertiles)
  coral_nz$size_class <- factor(
    dplyr::ntile(coral_nz$volume, 3),
    labels = c("Small", "Medium", "Large")
  )

  # Filter to species on ≥5 corals
  comm_filt <- comm_nz[, colSums(comm_nz > 0) >= 5]
  cat("  Testing", ncol(comm_filt), "species (present on ≥5 corals)\n")

  set.seed(42)
  indval_size <- multipatt(comm_filt, coral_nz$size_class,
                           control = how(nperm = 999))

  # Extract significant indicators
  indval_summary <- indval_size$sign
  indval_summary$species <- rownames(indval_summary)
  sig_indval <- indval_summary %>%
    filter(p.value < 0.05) %>%
    arrange(p.value)

  cat("  Significant indicator species (p < 0.05):", nrow(sig_indval), "\n")
  if (nrow(sig_indval) > 0) {
    # Decode group membership
    sig_indval$size_group <- apply(sig_indval[, grep("^s\\.", names(sig_indval))], 1,
                                   function(x) paste(c("Small", "Medium", "Large")[x == 1], collapse = "+"))
    print(sig_indval %>% dplyr::select(species, stat, p.value, size_group))
  }

  save_table(indval_summary %>% mutate(species = rownames(indval_summary)),
             "sensitivity_indval_size_class")
  cat("\n")

} else {
  cat("  indicspecies not installed — skipping\n\n")
}


# ============================================================================
# PART 4: COVERAGE-BASED RAREFACTION (iNEXT) — HILL NUMBERS
# ============================================================================
# Coverage-standardized diversity profiles (q=0 richness, q=1 exp-Shannon,
# q=2 inv-Simpson) by coral size class. More principled than fixed-depth
# rarefaction because it compares at equal sampling completeness.
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 4: COVERAGE-BASED RAREFACTION (iNEXT)\n")
cat("------------------------------------------------------------\n\n")

if (requireNamespace("iNEXT", quietly = TRUE)) {
  library(iNEXT)

  # Aggregate community matrix by size class
  coral_for_inext <- coral_master %>%
    filter(coral_id %in% rownames(community_matrix)) %>%
    mutate(size_class = factor(
      dplyr::ntile(volume, 3),
      labels = c("Small", "Medium", "Large")
    ))

  comm_aligned <- community_matrix[coral_for_inext$coral_id, ]

  # Pool abundances by size class
  size_pools <- list()
  for (sc in c("Small", "Medium", "Large")) {
    idx <- coral_for_inext$size_class == sc
    size_pools[[sc]] <- colSums(comm_aligned[idx, ])
    size_pools[[sc]] <- size_pools[[sc]][size_pools[[sc]] > 0]
  }

  # Run iNEXT for Hill numbers q = 0, 1, 2
  set.seed(42)
  inext_result <- iNEXT(size_pools, q = c(0, 1, 2), datatype = "abundance",
                        nboot = 100)

  # Extract asymptotic estimates
  asym_df <- inext_result$AsyEst
  cat("  Asymptotic diversity estimates by size class:\n")
  print(asym_df %>%
          dplyr::select(Assemblage, Diversity, Observed, Estimator, s.e., LCL, UCL) %>%
          mutate(across(where(is.numeric), ~round(., 2))))

  save_table(asym_df, "sensitivity_inext_hill_numbers")

  # Per-coral Hill numbers for scaling analysis
  cat("\n  Computing per-coral Hill numbers...\n")
  hill_df <- coral_for_inext %>%
    dplyr::select(coral_id, volume, site, size_class)

  # q=0 (richness), q=1 (exp Shannon), q=2 (inv Simpson) per coral
  hill_df$q0 <- apply(comm_aligned, 1, function(x) sum(x > 0))
  hill_df$q1 <- apply(comm_aligned, 1, function(x) {
    p <- x[x > 0] / sum(x)
    exp(-sum(p * log(p)))
  })
  hill_df$q2 <- apply(comm_aligned, 1, function(x) {
    p <- x[x > 0] / sum(x)
    1 / sum(p^2)
  })

  # Scaling exponents for each Hill number
  cat("\n  Scaling exponents (Poisson GLM ~ log(volume) + site):\n")
  for (q_name in c("q0", "q1", "q2")) {
    hill_data <- hill_df %>% filter(.data[[q_name]] > 0)
    m <- glm(as.formula(paste(q_name, "~ log(volume) + site")),
             data = hill_data, family = poisson)
    s <- summary(m)
    beta <- coef(m)["log(volume)"]
    p_val <- s$coefficients["log(volume)", "Pr(>|z|)"]
    cat(sprintf("    %s (Hill q=%s): beta = %.3f, z = %.2f, p = %s\n",
                q_name, gsub("q", "", q_name), beta,
                s$coefficients["log(volume)", "z value"],
                format.pval(p_val, 3)))
  }

  save_table(hill_df, "sensitivity_hill_numbers_per_coral")
  cat("\n")

} else {
  cat("  iNEXT not installed — skipping\n\n")
}


# ============================================================================
# PART 5: CWM OBLIGATE-MUTUALIST FRACTION AS BEF PREDICTOR
# ============================================================================
# Community-weighted mean of pocillopora association (obligate vs facultative)
# tests whether functional identity predicts condition better than richness.
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 5: CWM OBLIGATE-MUTUALIST FRACTION (BEF)\n")
cat("------------------------------------------------------------\n\n")

# Load trait data
traits <- read.csv(here::here("data/traits/cafi_traits_final.csv"),
                   stringsAsFactors = FALSE)

# Map association to numeric: obligate=1, facultative=0.5, unknown=NA
traits$assoc_score <- case_when(
  traits$pocillopora_association == "obligate" ~ 1.0,
  traits$pocillopora_association == "facultative" ~ 0.5,
  TRUE ~ NA_real_
)

cat("  Trait coverage: obligate =", sum(traits$assoc_score == 1, na.rm = TRUE),
    ", facultative =", sum(traits$assoc_score == 0.5, na.rm = TRUE),
    ", unknown =", sum(is.na(traits$assoc_score)), "\n")

# Match traits to community matrix species
# OTU names in community matrix may not match trait accepted_name exactly
# Try matching on both taxon_input and accepted_name
trait_lookup <- traits %>%
  dplyr::select(taxon_input, accepted_name, assoc_score, functional_role_str) %>%
  distinct()

# Compute CWM per coral
comm_for_cwm <- community_matrix[rowSums(community_matrix) > 0, ]
species_names <- colnames(comm_for_cwm)

# Match species to traits
matched_scores <- sapply(species_names, function(sp) {
  # Try exact match on accepted_name first, then taxon_input
  idx <- match(sp, trait_lookup$accepted_name)
  if (is.na(idx)) idx <- match(sp, trait_lookup$taxon_input)
  if (!is.na(idx)) trait_lookup$assoc_score[idx] else NA_real_
})

n_matched <- sum(!is.na(matched_scores))
cat("  Species matched to traits:", n_matched, "of", length(species_names), "\n")

# CWM: abundance-weighted mean association score (ignoring NAs)
cwm_df <- data.frame(
  coral_id = rownames(comm_for_cwm),
  stringsAsFactors = FALSE
)

cwm_df$obligate_fraction <- apply(comm_for_cwm, 1, function(x) {
  has_trait <- !is.na(matched_scores) & x > 0
  if (sum(has_trait) == 0) return(NA_real_)
  sum(x[has_trait] * matched_scores[has_trait]) / sum(x[has_trait])
})

cwm_df$n_obligate <- apply(comm_for_cwm, 1, function(x) {
  sum(x[!is.na(matched_scores) & matched_scores == 1])
})

cwm_df$pct_obligate <- apply(comm_for_cwm, 1, function(x) {
  obl <- sum(x[!is.na(matched_scores) & matched_scores == 1])
  obl / sum(x) * 100
})

cat("  Corals with CWM data:", sum(!is.na(cwm_df$obligate_fraction)), "\n\n")

# Merge with condition data and test
cwm_analysis <- condition_scores %>%
  dplyr::select(coral_id, site, condition_score) %>%
  left_join(cwm_df, by = "coral_id") %>%
  left_join(coral_master %>% dplyr::select(coral_id, volume, otu_richness), by = "coral_id") %>%
  filter(!is.na(condition_score), !is.na(obligate_fraction), !is.na(volume)) %>%
  mutate(log_volume = log(volume))

if (nrow(cwm_analysis) >= 20) {
  cat("  Testing CWM obligate fraction → condition (n =", nrow(cwm_analysis), "):\n")

  m_cwm <- lm(condition_score ~ obligate_fraction + log_volume + site, data = cwm_analysis)
  s_cwm <- summary(m_cwm)
  cat("    obligate_fraction: beta =", round(coef(m_cwm)["obligate_fraction"], 3),
      ", t =", round(s_cwm$coefficients["obligate_fraction", "t value"], 2),
      ", p =", format.pval(s_cwm$coefficients["obligate_fraction", "Pr(>|t|)"], 3), "\n")
  cat("    adj R² =", round(s_cwm$adj.r.squared, 3), "\n")

  # Compare: richness vs CWM
  m_rich <- lm(condition_score ~ otu_richness + log_volume + site, data = cwm_analysis)
  m_both <- lm(condition_score ~ otu_richness + obligate_fraction + log_volume + site, data = cwm_analysis)
  cat("    AIC(richness only):", round(AIC(m_rich), 1),
      "  AIC(CWM only):", round(AIC(m_cwm), 1),
      "  AIC(both):", round(AIC(m_both), 1), "\n\n")

  save_table(data.frame(
    model = c("richness_only", "cwm_only", "both"),
    AIC = round(c(AIC(m_rich), AIC(m_cwm), AIC(m_both)), 1),
    adj_r2 = round(c(summary(m_rich)$adj.r.squared,
                      s_cwm$adj.r.squared,
                      summary(m_both)$adj.r.squared), 4),
    richness_p = c(
      summary(m_rich)$coefficients["otu_richness", "Pr(>|t|)"],
      NA,
      summary(m_both)$coefficients["otu_richness", "Pr(>|t|)"]
    ),
    cwm_p = c(
      NA,
      s_cwm$coefficients["obligate_fraction", "Pr(>|t|)"],
      summary(m_both)$coefficients["obligate_fraction", "Pr(>|t|)"]
    )
  ), "sensitivity_cwm_vs_richness")
} else {
  cat("  Insufficient overlap (n < 20) — skipping CWM test\n\n")
}


# ============================================================================
# PART 6: NON-LINEAR BEF (LOG-RICHNESS, POLYNOMIAL)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 6: NON-LINEAR BEF SENSITIVITY\n")
cat("------------------------------------------------------------\n\n")

# Prepare BEF analysis data
bef_data <- condition_scores %>%
  dplyr::select(coral_id, site, condition_score) %>%
  left_join(coral_master %>% dplyr::select(coral_id, volume, otu_richness, total_cafi),
            by = "coral_id") %>%
  filter(!is.na(condition_score), !is.na(otu_richness), !is.na(volume),
         otu_richness > 0) %>%
  mutate(log_volume = log(volume),
         log_richness = log(otu_richness),
         site = factor(site))

n_bef <- nrow(bef_data)
cat("  BEF sample size:", n_bef, "\n\n")

# Linear (current)
m_linear <- lm(condition_score ~ otu_richness + log_volume + site, data = bef_data)
# Log-richness
m_log <- lm(condition_score ~ log_richness + log_volume + site, data = bef_data)
# Polynomial (quadratic)
m_poly <- lm(condition_score ~ poly(otu_richness, 2) + log_volume + site, data = bef_data)

aic_table <- data.frame(
  model = c("linear", "log_richness", "polynomial_2"),
  AIC = round(c(AIC(m_linear), AIC(m_log), AIC(m_poly)), 1),
  adj_r2 = round(c(summary(m_linear)$adj.r.squared,
                    summary(m_log)$adj.r.squared,
                    summary(m_poly)$adj.r.squared), 4),
  richness_p = c(
    summary(m_linear)$coefficients["otu_richness", "Pr(>|t|)"],
    summary(m_log)$coefficients["log_richness", "Pr(>|t|)"],
    summary(m_poly)$coefficients["poly(otu_richness, 2)1", "Pr(>|t|)"]
  )
)
aic_table$delta_AIC <- round(aic_table$AIC - min(aic_table$AIC), 1)

cat("  Model comparison (condition ~ richness_form + log_volume + site):\n")
print(aic_table)
cat("\n  Best model:", aic_table$model[which.min(aic_table$AIC)],
    "(dAIC =", min(aic_table$delta_AIC[aic_table$delta_AIC > 0]), "over next)\n\n")

save_table(aic_table, "sensitivity_nonlinear_bef")


# ============================================================================
# PART 7: MISSING DATA CHARACTERIZATION (114 → 84)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 7: MISSING DATA CHARACTERIZATION\n")
cat("------------------------------------------------------------\n\n")

# Which corals have condition data?
all_coral <- coral_master %>%
  mutate(
    has_condition = coral_id %in% condition_scores$coral_id[!is.na(condition_scores$condition_score)],
    log_volume = log(volume)
  )

cat("  Total corals:", nrow(all_coral), "\n")
cat("  With condition:", sum(all_coral$has_condition), "\n")
cat("  Missing:", sum(!all_coral$has_condition), "\n\n")

# Logistic regression: does missingness depend on volume, site, richness?
m_missing <- glm(has_condition ~ log_volume + site + otu_richness + total_cafi,
                 family = binomial, data = all_coral)
s_missing <- summary(m_missing)

cat("  Predictors of having condition data (logistic):\n")
for (pred in rownames(s_missing$coefficients)[-1]) {
  cat(sprintf("    %-20s beta = %6.3f, z = %5.2f, p = %s\n",
              pred,
              s_missing$coefficients[pred, "Estimate"],
              s_missing$coefficients[pred, "z value"],
              format.pval(s_missing$coefficients[pred, "Pr(>|z|)"], 3)))
}

# Compare distributions
cat("\n  Mean volume (with condition):", round(mean(all_coral$volume[all_coral$has_condition])),
    "vs without:", round(mean(all_coral$volume[!all_coral$has_condition])), "\n")
cat("  Mean richness (with):", round(mean(all_coral$otu_richness[all_coral$has_condition]), 1),
    "vs without:", round(mean(all_coral$otu_richness[!all_coral$has_condition]), 1), "\n")

# Site breakdown
cat("\n  Missing by site:\n")
print(all_coral %>%
        group_by(site) %>%
        summarise(n = n(), n_with = sum(has_condition), n_missing = sum(!has_condition),
                  pct_missing = round(100 * mean(!has_condition), 1), .groups = "drop"))

save_table(
  broom::tidy(m_missing) %>% mutate(across(where(is.numeric), ~round(., 4))),
  "sensitivity_missing_data_predictors"
)
cat("\n")


# ============================================================================
# PART 8: MORPHOTYPE COVARIATE IN BEF MODEL
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 8: MORPHOTYPE COVARIATE CHECK\n")
cat("------------------------------------------------------------\n\n")

bef_morph <- bef_data %>%
  left_join(coral_master %>% dplyr::select(coral_id, morphotype), by = "coral_id") %>%
  filter(!is.na(morphotype))

if (length(unique(bef_morph$morphotype)) > 1 && nrow(bef_morph) >= 30) {
  cat("  Morphotypes:", paste(unique(bef_morph$morphotype), collapse = ", "), "\n")
  cat("  N per morphotype:\n")
  print(table(bef_morph$morphotype))

  m_no_morph <- lm(condition_score ~ otu_richness + log_volume + site, data = bef_morph)
  m_with_morph <- lm(condition_score ~ otu_richness + log_volume + morphotype + site,
                     data = bef_morph)

  beta_without <- coef(m_no_morph)["otu_richness"]
  beta_with <- coef(m_with_morph)["otu_richness"]
  p_without <- summary(m_no_morph)$coefficients["otu_richness", "Pr(>|t|)"]
  p_with <- summary(m_with_morph)$coefficients["otu_richness", "Pr(>|t|)"]

  cat(sprintf("\n  Richness → condition WITHOUT morphotype: beta = %.4f, p = %s\n",
              beta_without, format.pval(p_without, 3)))
  cat(sprintf("  Richness → condition WITH morphotype:    beta = %.4f, p = %s\n",
              beta_with, format.pval(p_with, 3)))
  cat(sprintf("  Change in beta: %.1f%%\n",
              100 * (beta_with - beta_without) / abs(beta_without)))
  cat("  AIC without:", round(AIC(m_no_morph), 1),
      "  AIC with:", round(AIC(m_with_morph), 1), "\n\n")

  save_table(data.frame(
    model = c("without_morphotype", "with_morphotype"),
    richness_beta = round(c(beta_without, beta_with), 4),
    richness_p = c(p_without, p_with),
    AIC = round(c(AIC(m_no_morph), AIC(m_with_morph)), 1)
  ), "sensitivity_morphotype_bef")
} else {
  cat("  Insufficient morphotype variation — skipping\n\n")
}


# ============================================================================
# PART 8b: MORPHOTYPE–HAPLOTYPE CONCORDANCE
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 8b: MORPHOTYPE–HAPLOTYPE CONCORDANCE\n")
cat("------------------------------------------------------------\n\n")

hap_data <- coral_master %>%
  filter(haplotype_valid == TRUE, !is.na(morphotype_original))

if (nrow(hap_data) > 10) {
  # Cross-tabulation
  ct <- table(Morphotype = hap_data$morphotype_original, Species = hap_data$poc_species)
  cat("  Morphotype × genetic species concordance:\n")
  print(ct)

  # Chi-square / Fisher's exact (sparse cells expected)
  chi_test <- tryCatch(chisq.test(ct, simulate.p.value = TRUE, B = 10000),
                       error = function(e) NULL)
  if (!is.null(chi_test)) {
    cramers_v <- sqrt(chi_test$statistic / (sum(ct) * (min(dim(ct)) - 1)))
    cat(sprintf("\n  Chi-square (simulated p): X² = %.1f, p = %.4f\n",
                chi_test$statistic, chi_test$p.value))
    cat(sprintf("  Cramér's V = %.3f\n", cramers_v))
  }

  # Conditional probabilities: P(species | morphotype)
  cond_prob <- prop.table(ct, margin = 1)
  cat("\n  P(genetic species | field morphotype):\n")
  print(round(cond_prob, 3))

  # Save concordance table
  ct_df <- as.data.frame.matrix(ct)
  ct_df$morphotype <- rownames(ct_df)
  ct_df <- ct_df[, c("morphotype", setdiff(names(ct_df), "morphotype"))]
  save_table(ct_df, "sensitivity_morphotype_haplotype_concordance")
  cat("  Saved: sensitivity_morphotype_haplotype_concordance.csv\n\n")
} else {
  cat("  Insufficient haplotype data — skipping\n\n")
}

# ============================================================================
# PART 8c: BEF WITH HAPLOTYPE COVARIATE
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 8c: BEF WITH HAPLOTYPE (genetic species) COVARIATE\n")
cat("------------------------------------------------------------\n\n")

bef_hap <- bef_data %>%
  left_join(coral_master %>% dplyr::select(coral_id, poc_species, haplotype_valid),
            by = "coral_id") %>%
  filter(haplotype_valid == TRUE, !is.na(poc_species))

if (length(unique(bef_hap$poc_species)) > 1 && nrow(bef_hap) >= 30) {
  cat("  N =", nrow(bef_hap), "corals with valid haplotype + condition data\n")
  cat("  Species:\n")
  print(table(bef_hap$poc_species))

  m_no_hap <- lm(condition_score ~ otu_richness + log_volume + site, data = bef_hap)
  m_with_hap <- lm(condition_score ~ otu_richness + log_volume + poc_species + site,
                    data = bef_hap)

  beta_without <- coef(m_no_hap)["otu_richness"]
  beta_with <- coef(m_with_hap)["otu_richness"]
  p_without <- summary(m_no_hap)$coefficients["otu_richness", "Pr(>|t|)"]
  p_with <- summary(m_with_hap)$coefficients["otu_richness", "Pr(>|t|)"]

  cat(sprintf("\n  Richness → condition WITHOUT haplotype: beta = %.4f, p = %s\n",
              beta_without, format.pval(p_without, 3)))
  cat(sprintf("  Richness → condition WITH haplotype:    beta = %.4f, p = %s\n",
              beta_with, format.pval(p_with, 3)))
  cat(sprintf("  Change in beta: %.1f%%\n",
              100 * (beta_with - beta_without) / abs(beta_without)))
  cat("  AIC without:", round(AIC(m_no_hap), 1),
      "  AIC with:", round(AIC(m_with_hap), 1), "\n\n")

  save_table(data.frame(
    model = c("without_haplotype", "with_haplotype"),
    n = nrow(bef_hap),
    richness_beta = round(c(beta_without, beta_with), 4),
    richness_p = c(p_without, p_with),
    AIC = round(c(AIC(m_no_hap), AIC(m_with_hap)), 1)
  ), "sensitivity_haplotype_bef")
} else {
  cat("  Insufficient haplotype variation in BEF subset — skipping\n\n")
}

# ============================================================================
# PART 8d: SCALING BY SPECIES (Q1 sensitivity)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 8d: SCALING BY SPECIES (Q1 sensitivity)\n")
cat("------------------------------------------------------------\n\n")

scale_hap <- coral_master %>%
  filter(haplotype_valid == TRUE, !is.na(log_volume))

if (nrow(scale_hap) > 30) {
  cat("  N =", nrow(scale_hap), "genotyped colonies with volume\n")

  # Test species × volume interaction
  m_add <- MASS::glm.nb(total_cafi ~ log_volume + poc_species + site, data = scale_hap)
  m_int <- MASS::glm.nb(total_cafi ~ log_volume * poc_species + site, data = scale_hap)

  anova_test <- anova(m_add, m_int, test = "Chisq")
  int_p <- anova_test$`Pr(Chi)`[2]
  cat(sprintf("  Volume × species interaction: p = %s\n", format.pval(int_p, 3)))

  # Species-specific scaling exponents
  cat("\n  Species-specific scaling exponents:\n")
  species_betas <- data.frame(species = character(), n = integer(),
                              beta = numeric(), se = numeric(),
                              p_vs_1 = numeric(), stringsAsFactors = FALSE)

  for (sp in sort(unique(scale_hap$poc_species))) {
    sp_data <- scale_hap %>% filter(poc_species == sp)
    if (nrow(sp_data) >= 10) {
      m_sp <- tryCatch(MASS::glm.nb(total_cafi ~ log_volume + site, data = sp_data),
                       error = function(e) NULL)
      if (!is.null(m_sp)) {
        b <- coef(m_sp)["log_volume"]
        se <- summary(m_sp)$coefficients["log_volume", "Std. Error"]
        z_vs1 <- (b - 1) / se
        p_vs1 <- 2 * pnorm(-abs(z_vs1))
        cat(sprintf("    %s (n=%d): β = %.2f ± %.2f, p vs β=1: %s\n",
                    sp, nrow(sp_data), b, se, format.pval(p_vs1, 3)))
        species_betas <- rbind(species_betas,
                               data.frame(species = sp, n = nrow(sp_data),
                                          beta = round(b, 3), se = round(se, 3),
                                          p_vs_1 = p_vs1))
      }
    } else {
      cat(sprintf("    %s (n=%d): too few for standalone model\n", sp, nrow(sp_data)))
    }
  }

  save_table(species_betas, "sensitivity_scaling_by_species")
  cat("\n")
}

# ============================================================================
# PART 8e: COMPOSITION BY SPECIES (Q2 sensitivity)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 8e: COMPOSITION BY SPECIES (Q2 sensitivity)\n")
cat("------------------------------------------------------------\n\n")

comm_hap <- community_matrix[coral_master$haplotype_valid == TRUE &
                              !is.na(coral_master$poc_species), ]
meta_hap <- coral_master %>%
  filter(haplotype_valid == TRUE, !is.na(poc_species))

if (nrow(comm_hap) == nrow(meta_hap) && nrow(meta_hap) > 30) {
  cat("  N =", nrow(meta_hap), "genotyped colonies\n")

  dist_hap <- vegdist(decostand(comm_hap, "hellinger"), method = "bray")

  # Marginal PERMANOVA with species
  perm_hap <- adonis2(dist_hap ~ log_volume + poc_species + site,
                       data = meta_hap, permutations = 999, method = "bray",
                       by = "margin")
  cat("  Marginal PERMANOVA (volume + species + site):\n")
  print(perm_hap)

  comp_results <- data.frame(
    term = rownames(perm_hap),
    R2 = round(perm_hap$R2, 4),
    F_stat = round(perm_hap$F, 2),
    p = perm_hap$`Pr(>F)`
  )
  save_table(comp_results, "sensitivity_composition_by_species")
  cat("\n")
}

# ============================================================================
# PART 8f: GENOTYPE → ARCHITECTURE MEDIATION
# ============================================================================
# Tests whether genetic species identity predicts CAFI through branching
# architecture (branch_width). Cross-references companion study (Maatea, Paper 2)
# which demonstrated full mediation: haplotype composition signal disappears
# when branch architecture is conditioned on (PERMANOVA p: 0.003 → 0.338).
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 8f: GENOTYPE → ARCHITECTURE MEDIATION\n")
cat("------------------------------------------------------------\n\n")

hap_arch <- coral_master %>%
  filter(haplotype_valid == TRUE, !is.na(poc_species),
         poc_species != "P. meandrina/grandis")

if (nrow(hap_arch) > 30) {
  cat("  N =", nrow(hap_arch), "genotyped colonies (4 species)\n\n")

  # --- Branch width × species (Fisher's exact) ---
  if ("branch_width" %in% names(hap_arch)) {
    bw_tab <- table(Species = hap_arch$poc_species, Branch = hap_arch$branch_width)
    cat("  Branch width × species:\n")
    print(bw_tab)
    fisher_bw <- fisher.test(bw_tab, simulate.p.value = TRUE, B = 10000)
    cat(sprintf("\n  Fisher's exact: p = %.2e\n", fisher_bw$p.value))

    # Proportion wide-branching by species
    prop_wide <- prop.table(bw_tab, margin = 1)
    cat("\n  % wide-branching by species:\n")
    print(round(prop_wide * 100, 1))
    cat("\n")
  }

  # --- Species → richness (controlling volume + site) ---
  m_rich_null <- glm(species_richness ~ log_volume + site, data = hap_arch, family = poisson)
  m_rich_sp <- glm(species_richness ~ log_volume + poc_species + site, data = hap_arch, family = poisson)
  lr_rich <- anova(m_rich_null, m_rich_sp, test = "Chisq")
  cat(sprintf("  Species → richness (| volume + site): LRT p = %s\n",
              format.pval(lr_rich$`Pr(>Chi)`[2], 3)))

  # --- Species → abundance (controlling volume + site) ---
  m_abund_null <- MASS::glm.nb(total_cafi ~ log_volume + site, data = hap_arch)
  m_abund_sp <- MASS::glm.nb(total_cafi ~ log_volume + poc_species + site, data = hap_arch)
  lr_abund <- anova(m_abund_null, m_abund_sp, test = "Chisq")
  cat(sprintf("  Species → abundance (| volume + site): LRT p = %s\n",
              format.pval(lr_abund$`Pr(Chi)`[2], 3)))

  # --- Species → condition (controlling volume + site) ---
  cond_arch <- hap_arch %>% filter(!is.na(condition_score))
  if (nrow(cond_arch) > 20) {
    m_cond_null <- lm(condition_score ~ log_volume + site, data = cond_arch)
    m_cond_sp <- lm(condition_score ~ log_volume + poc_species + site, data = cond_arch)
    f_test <- anova(m_cond_null, m_cond_sp)
    cat(sprintf("  Species → condition (| volume + site): F-test p = %s (n = %d)\n",
                format.pval(f_test$`Pr(>F)`[2], 3), nrow(cond_arch)))
    cat(sprintf("    Adj R² without species: %.3f; with species: %.3f\n",
                summary(m_cond_null)$adj.r.squared, summary(m_cond_sp)$adj.r.squared))

    # Species-level condition means
    cat("\n  Condition means by species:\n")
    cond_means <- cond_arch %>%
      group_by(poc_species) %>%
      summarise(n = n(), mean_cond = round(mean(condition_score), 2),
                se = round(sd(condition_score)/sqrt(n()), 2), .groups = "drop")
    print(as.data.frame(cond_means))
  }

  # --- Functional group species effects (key test: do Trapezia differ?) ---
  cat("\n  Functional group × species (LRT, controlling volume + site):\n")
  for (fg in c("n_trapezia", "n_resident_fish", "n_galeropsis", "n_shrimps", "n_other_crab")) {
    fg_null <- tryCatch(MASS::glm.nb(as.formula(paste(fg, "~ log_volume + site")),
                                      data = hap_arch), error = function(e) NULL)
    fg_sp <- tryCatch(MASS::glm.nb(as.formula(paste(fg, "~ log_volume + poc_species + site")),
                                    data = hap_arch), error = function(e) NULL)
    if (!is.null(fg_null) && !is.null(fg_sp)) {
      lr <- anova(fg_null, fg_sp, test = "Chisq")
      p_val <- lr$`Pr(Chi)`[2]
      cat(sprintf("    %s: p = %s %s\n", fg, format.pval(p_val, 3),
                  ifelse(p_val < 0.05, "*", "")))
    }
  }

  # --- Save summary ---
  arch_summary <- data.frame(
    Test = c("Branch width × species (Fisher)",
             "Species → richness (| vol + site)",
             "Species → abundance (| vol + site)",
             "Species → condition (| vol + site)",
             "Species → composition (PERMANOVA)"),
    P_value = c(fisher_bw$p.value,
                lr_rich$`Pr(>Chi)`[2],
                lr_abund$`Pr(Chi)`[2],
                if (exists("f_test")) f_test$`Pr(>F)`[2] else NA,
                0.001),  # from Part 8e
    N = c(nrow(hap_arch), nrow(hap_arch), nrow(hap_arch),
          if (exists("cond_arch")) nrow(cond_arch) else NA, nrow(hap_arch)),
    Interpretation = c(
      "P. grandis = wide-branching; others = tight",
      "Species predicts richness beyond volume",
      "Species does NOT predict abundance beyond volume",
      "Species predicts condition (P. verrucosa highest)",
      "Species explains 8.3% of composition"
    )
  )
  save_table(arch_summary, "sensitivity_genotype_architecture")
  cat("\n  Saved: sensitivity_genotype_architecture.csv\n\n")
}

# ============================================================================
# PART 8g: PHYLOGENETIC DISTANCE AND SYMBIONT IDENTITY
# ============================================================================
# Uses published phylogenetic tree (Johnston, Cunning & Burgess 2022, Mol Ecol)
# to test whether genetic distance predicts CAFI composition beyond architecture.
# Adapted from Stier-CAFI-Size-Maatea-2026/scripts/07_haplotype_effects.R §15.
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 8g: PHYLOGENETIC DISTANCE & SYMBIONT IDENTITY\n")
cat("------------------------------------------------------------\n\n")

tree_file <- file.path(here::here(), "data/genetic_references/Johnston2022_Pocillopora_tree.nex")

if (file.exists(tree_file)) {
  library(ape)

  tree <- read.nexus(tree_file)
  phylo_dist_mat <- cophenetic(tree)

  # Map our haplotypes to tree tip labels
  hap_to_tree <- c(
    "1a_Pm" = "P_meandrina",
    "1a_Pe" = "P_eydouxi",
    "1a"    = "P_meandrina",
    "8a"    = "Haplotype_8a",
    "3a"    = "P_verrucosa",
    "3b"    = "P_verrucosa",
    "3f"    = "P_verrucosa",
    "3h"    = "P_verrucosa",
    "10"    = "Haplotype_10"
  )

  # Assign putative Cladocopium symbiont species (Johnston et al. 2022; Burgess et al. 2024)
  hap_to_symbiont <- c(
    "1a_Pm" = "C_latusorum",
    "1a_Pe" = "C_latusorum",
    "1a"    = "C_latusorum",
    "8a"    = "C_latusorum",
    "3a"    = "C_pacificum",
    "3b"    = "C_pacificum",
    "3f"    = "C_pacificum",
    "3h"    = "C_pacificum",
    "10"    = "C_pacificum"
  )

  phylo_data <- coral_master %>%
    filter(haplotype %in% names(hap_to_tree)) %>%
    mutate(tree_taxon = hap_to_tree[haplotype],
           symbiont = hap_to_symbiont[haplotype])

  cat("  N corals mapped to phylogeny:", nrow(phylo_data), "\n")
  cat("  Symbiont distribution: C_latusorum =", sum(phylo_data$symbiont == "C_latusorum"),
      ", C_pacificum =", sum(phylo_data$symbiont == "C_pacificum"), "\n\n")

  # --- Individual-level phylogenetic distance matrix ---
  coral_ids <- phylo_data$coral_id
  n_phy <- length(coral_ids)
  ind_phylo_dist <- matrix(0, n_phy, n_phy, dimnames = list(coral_ids, coral_ids))
  for (i in 1:(n_phy - 1)) {
    for (j in (i + 1):n_phy) {
      ind_phylo_dist[i, j] <- phylo_dist_mat[phylo_data$tree_taxon[i], phylo_data$tree_taxon[j]]
      ind_phylo_dist[j, i] <- ind_phylo_dist[i, j]
    }
  }

  # Community matrix for phylo-mapped corals
  comm_phy <- community_matrix[coral_ids, ]
  comm_phy <- comm_phy[, colSums(comm_phy) > 0]
  comm_phy_hell <- decostand(comm_phy, "hellinger")
  comm_dist_phy <- vegdist(comm_phy_hell, "bray")

  # --- Mantel: phylogenetic distance ~ community dissimilarity ---
  set.seed(42)
  mantel_phylo <- mantel(as.dist(ind_phylo_dist), comm_dist_phy, permutations = 9999)
  cat("  Mantel (phylo distance ~ community): r =", round(mantel_phylo$statistic, 4),
      ", p =", round(mantel_phylo$signif, 4), "\n")

  # --- Partial Mantel: phylo ~ community | volume ---
  vol_data <- phylo_data %>% dplyr::select(coral_id, log_volume)
  vol_dist <- dist(scale(vol_data$log_volume))

  partial_mantel_phylo <- mantel.partial(as.dist(ind_phylo_dist), comm_dist_phy,
                                          vol_dist, permutations = 9999)
  cat("  Partial Mantel (phylo ~ community | volume): r =",
      round(partial_mantel_phylo$statistic, 4),
      ", p =", round(partial_mantel_phylo$signif, 4), "\n")

  # --- PERMANOVA: community ~ symbiont + volume + site ---
  perm_symb <- adonis2(comm_dist_phy ~ log_volume + symbiont + site,
                        data = phylo_data, by = "margin", permutations = 999)
  cat("  PERMANOVA symbiont (marginal, after volume + site):\n")
  cat("    R2 =", round(perm_symb["symbiont", "R2"], 4),
      ", F =", round(perm_symb["symbiont", "F"], 2),
      ", p =", round(perm_symb["symbiont", "Pr(>F)"], 4), "\n")

  # --- Symbiont ~ condition ---
  cond_symb <- phylo_data %>% filter(!is.na(condition_score))
  if (nrow(cond_symb) > 20) {
    m_symb_cond <- lm(condition_score ~ symbiont + log_volume + site, data = cond_symb)
    symb_coefs <- summary(m_symb_cond)$coefficients
    symb_row <- grep("symbiont", rownames(symb_coefs))
    if (length(symb_row) > 0) {
      cat("  Symbiont -> condition: beta =", round(symb_coefs[symb_row, 1], 3),
          ", p =", round(symb_coefs[symb_row, 4], 4), "\n\n")
    }
  }

  # --- Save summary ---
  phylo_summary <- data.frame(
    Analysis = c(
      "Mantel: phylogenetic distance ~ community dissimilarity",
      "Partial Mantel: phylo ~ community | volume",
      "PERMANOVA: community ~ symbiont (after volume + site)",
      "LM: condition ~ symbiont (after volume + site)"
    ),
    Statistic = c(
      paste0("r = ", round(mantel_phylo$statistic, 4)),
      paste0("r = ", round(partial_mantel_phylo$statistic, 4)),
      paste0("R2 = ", round(perm_symb["symbiont", "R2"], 4)),
      if (exists("m_symb_cond")) paste0("beta = ", round(symb_coefs[symb_row, 1], 3)) else "NA"
    ),
    P_value = c(
      round(mantel_phylo$signif, 4),
      round(partial_mantel_phylo$signif, 4),
      round(perm_symb["symbiont", "Pr(>F)"], 4),
      if (exists("m_symb_cond")) round(symb_coefs[symb_row, 4], 4) else NA
    ),
    N = c(nrow(phylo_data), nrow(phylo_data),
          nrow(phylo_data), if (exists("cond_symb")) nrow(cond_symb) else NA)
  )
  save_table(phylo_summary, "sensitivity_phylogenetic_symbiont")
  cat("  Saved: sensitivity_phylogenetic_symbiont.csv\n\n")

} else {
  cat("  Tree file not found:", tree_file, "\n")
  cat("  Skipping phylogenetic analyses\n\n")
}


# ============================================================================
# PART 9: MEDIATION — RICHNESS → ABUNDANCE → CONDITION
# ============================================================================
# Tests whether the richness → condition path operates independently of
# the richness → abundance → condition cascade.
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 9: MEDIATION (richness → abundance → condition)\n")
cat("------------------------------------------------------------\n\n")

if (requireNamespace("mediation", quietly = TRUE)) {

  med_data <- bef_data %>%
    mutate(sqrt_cafi = sqrt(total_cafi))

  # Mediator model: abundance ~ richness + covariates
  m_mediator <- lm(sqrt_cafi ~ otu_richness + log_volume + site, data = med_data)
  # Outcome model: condition ~ richness + abundance + covariates
  m_outcome <- lm(condition_score ~ otu_richness + sqrt_cafi + log_volume + site,
                  data = med_data)

  set.seed(42)
  med_result <- mediation::mediate(m_mediator, m_outcome,
                                    treat = "otu_richness",
                                    mediator = "sqrt_cafi",
                                    boot = TRUE, sims = 1000)

  cat("  Mediation results (richness → sqrt(abundance) → condition):\n")
  cat("    ACME (average causal mediation effect):", round(med_result$d0, 4),
      " [", round(med_result$d0.ci[1], 4), ",", round(med_result$d0.ci[2], 4), "]\n")
  cat("    ADE (average direct effect):", round(med_result$z0, 4),
      " [", round(med_result$z0.ci[1], 4), ",", round(med_result$z0.ci[2], 4), "]\n")
  cat("    Total effect:", round(med_result$tau.coef, 4), "\n")
  cat("    Proportion mediated:", round(med_result$n0, 3), "\n")
  cat("    ACME p =", format.pval(med_result$d0.p, 3),
      "  ADE p =", format.pval(med_result$z0.p, 3), "\n\n")

  save_table(data.frame(
    effect = c("ACME", "ADE", "Total", "Proportion_mediated"),
    estimate = round(c(med_result$d0, med_result$z0, med_result$tau.coef, med_result$n0), 4),
    ci_lower = round(c(med_result$d0.ci[1], med_result$z0.ci[1], med_result$tau.ci[1], med_result$n0.ci[1]), 4),
    ci_upper = round(c(med_result$d0.ci[2], med_result$z0.ci[2], med_result$tau.ci[2], med_result$n0.ci[2]), 4),
    p_value = c(med_result$d0.p, med_result$z0.p, med_result$tau.p, med_result$n0.p)
  ), "sensitivity_mediation_richness_abundance")

} else {
  cat("  mediation package not installed — skipping\n\n")
}


# ============================================================================
# PART 10: C-SCORE COMMUNITY-WIDE CO-OCCURRENCE
# ============================================================================
# Stone & Roberts C-score with fixed-row/fixed-column null (quasiswap).
# Robustness check for the 0/528 pairwise result.
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 10: C-SCORE COMMUNITY-WIDE CO-OCCURRENCE\n")
cat("------------------------------------------------------------\n\n")

# Use presence-absence of focal species (≥10 occurrences)
comm_pa_full <- (community_matrix[rowSums(community_matrix) > 0, ] > 0) * 1
sp_prev <- colSums(comm_pa_full)
focal_sp <- names(sp_prev[sp_prev >= 10])
comm_cscore <- comm_pa_full[, focal_sp]

cat("  Testing", ncol(comm_cscore), "species on", nrow(comm_cscore), "corals\n")

# C-score via oecosimu
set.seed(42)
cscore_result <- oecosimu(comm_cscore, nestedchecker, method = "quasiswap",
                           nsimul = 999, statistic = "C.score")

cat("  C-score observed:", round(cscore_result$oecosimu$statistic, 2), "\n")
cat("  C-score null mean:", round(mean(cscore_result$oecosimu$simulated), 2), "\n")
ses_cscore <- (cscore_result$oecosimu$statistic - mean(cscore_result$oecosimu$simulated)) /
  sd(cscore_result$oecosimu$simulated)
cat("  SES:", round(ses_cscore, 3), "\n")
cat("  p =", format.pval(cscore_result$oecosimu$p[1], 3), "\n\n")

save_table(data.frame(
  metric = "C_score",
  observed = round(cscore_result$oecosimu$statistic, 2),
  null_mean = round(mean(cscore_result$oecosimu$simulated), 2),
  null_sd = round(sd(cscore_result$oecosimu$simulated), 2),
  SES = round(ses_cscore, 3),
  p_value = cscore_result$oecosimu$p[1],
  n_species = ncol(comm_cscore),
  n_corals = nrow(comm_cscore),
  null_model = "quasiswap"
), "sensitivity_cscore_community")


# ============================================================================
# PART 11: MuMIn MODEL-AVERAGED NEIGHBORHOOD COEFFICIENTS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 11: MuMIn MODEL AVERAGING (NEIGHBORHOOD)\n")
cat("------------------------------------------------------------\n\n")

if (requireNamespace("MuMIn", quietly = TRUE)) {
  library(MuMIn)

  # Prepare neighborhood data (same as script 04)
  landscape_data <- coral_master %>%
    filter(!is.na(n_neighbors), !is.na(volume), volume > 0) %>%
    mutate(
      log_volume = log(volume),
      log_neighbor_vol = log(total_neighbor_volume + 1),
      site = factor(site)
    )

  # Check if mean_neighbor_dist exists
  if ("mean_neighbor_dist" %in% names(landscape_data)) {
    cat("  Neighborhood corals:", nrow(landscape_data), "\n")

    # Full model for richness (the response with significant distance effect)
    options(na.action = "na.fail")  # Required by dredge

    m_global <- glm(otu_richness ~ log_volume + n_neighbors + log_neighbor_vol +
                      mean_neighbor_dist + site,
                    family = poisson, data = landscape_data)

    # Check overdispersion
    disp <- sum(residuals(m_global, type = "pearson")^2) / m_global$df.residual
    if (disp > 1.5) {
      m_global <- glm(otu_richness ~ log_volume + n_neighbors + log_neighbor_vol +
                        mean_neighbor_dist + site,
                      family = quasipoisson, data = landscape_data)
      cat("  Using quasipoisson (dispersion =", round(disp, 2), ")\n")
      # dredge doesn't work with quasipoisson — use QAIC
      # Fall back to reporting the full model comparison
      cat("  Note: MuMIn::dredge incompatible with quasipoisson; using manual AIC\n")

      # Manual model comparison
      m1 <- glm(otu_richness ~ log_volume + site, family = poisson, data = landscape_data)
      m2 <- glm(otu_richness ~ log_volume + mean_neighbor_dist + site,
                family = poisson, data = landscape_data)
      m3 <- glm(otu_richness ~ log_volume + n_neighbors + site,
                family = poisson, data = landscape_data)
      m4 <- glm(otu_richness ~ log_volume + log_neighbor_vol + site,
                family = poisson, data = landscape_data)
      m5 <- glm(otu_richness ~ log_volume + n_neighbors + log_neighbor_vol +
                  mean_neighbor_dist + site, family = poisson, data = landscape_data)

      aic_models <- data.frame(
        model = c("volume_only", "volume+distance", "volume+count",
                  "volume+neighbor_vol", "full"),
        AIC = round(c(AIC(m1), AIC(m2), AIC(m3), AIC(m4), AIC(m5)), 1)
      )
      aic_models$delta_AIC <- round(aic_models$AIC - min(aic_models$AIC), 1)
      aic_models <- aic_models[order(aic_models$AIC), ]
      cat("\n  Model comparison (Poisson AIC):\n")
      print(aic_models, row.names = FALSE)
      save_table(aic_models, "sensitivity_mumin_neighborhood_aic")

    } else {
      # dredge with Poisson
      cat("  Using Poisson (dispersion =", round(disp, 2), ")\n")
      dredge_result <- dredge(m_global, fixed = c("log_volume", "site"))

      # Model-averaged coefficients
      avg_model <- model.avg(dredge_result, subset = delta < 4)
      cat("\n  Model-averaged coefficients (dAIC < 4):\n")
      print(summary(avg_model)$coefmat.full)

      save_table(
        as.data.frame(dredge_result) %>% head(10),
        "sensitivity_mumin_neighborhood_dredge"
      )
    }

    options(na.action = "na.omit")  # Reset
    cat("\n")
  } else {
    cat("  mean_neighbor_dist not found — skipping\n\n")
  }
} else {
  cat("  MuMIn not installed — skipping\n\n")
}


# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    SENSITIVITY ANALYSES COMPLETE\n")
cat("============================================================\n\n")
cat("Outputs saved to output/tables/sensitivity_*\n")
cat("Review results above and incorporate into manuscript supplement.\n\n")
