# ============================================================================
# 15_reviewer_response_analyses.R - Pre-Built Analyses for Reviewer Comments
# ============================================================================
#
# PURPOSE: Code ready to run if reviewers request specific analyses.
#          NOT part of the manuscript pipeline — run ad hoc as needed.
#
# SECTIONS:
#   A. mvabund model-based MANOVA (alternative to PERMANOVA)
#   B. GDM nonlinear compositional turnover
#   C. TITAN2 community breakpoints along size gradient
#   D. Bayesian SEM for BEF path model (brms)
#   E. LOOCV cross-validation of richness-condition model
#   F. Bayesian multilevel scaling (brms, partial pooling across species)
#   G. Cross-study meta-analysis of scaling exponents (metafor)
#
# DEPENDENCIES: 00_setup.R, 01_load_data.R
# INSTALL: Some require additional packages — install.packages() noted inline
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-03-22
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    REVIEWER RESPONSE ANALYSES (run as needed)\n")
cat("============================================================\n\n")

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

coral_master <- load_object("coral_master")
community_matrix <- load_object("community_matrix")
condition_scores <- load_object("condition_scores")

set.seed(42)

# ============================================================================
# A. mvabund MODEL-BASED MANOVA
# ============================================================================
# Alternative to distance-based PERMANOVA. Fits per-species NB GLMs and uses
# resampling for inference. Better handles mean-variance relationships.
#
# install.packages("mvabund")
# ============================================================================

if (FALSE) {  # Set to TRUE to run

cat("--- A. mvabund ---\n\n")

if (requireNamespace("mvabund", quietly = TRUE)) {
  library(mvabund)

  comm_nz <- community_matrix[rowSums(community_matrix) > 0, ]
  # Filter to species on ≥5 corals (reduces fitting issues)
  comm_filt <- comm_nz[, colSums(comm_nz > 0) >= 5]
  coral_sub <- coral_master[match(rownames(comm_filt), coral_master$coral_id), ]

  mv_data <- mvabund(comm_filt)
  m_mv <- manyglm(mv_data ~ log(volume) + site, data = coral_sub,
                   family = "negative.binomial")

  # This takes 5-15 minutes
  cat("  Running mvabund anova (may take 10+ minutes)...\n")
  mv_anova <- anova(m_mv, p.uni = "adjusted", nBoot = 999)
  print(mv_anova)

  # Per-species p-values for volume effect
  sp_pvals <- mv_anova$uni.p["log(volume)", ]
  sig_sp <- names(sp_pvals[sp_pvals < 0.05])
  cat("\n  Species with significant volume effect (p < 0.05):", length(sig_sp), "\n")
  cat("  ", paste(sig_sp, collapse = ", "), "\n")

  save_table(data.frame(
    species = names(sp_pvals),
    p_value = sp_pvals,
    significant = sp_pvals < 0.05
  ), "reviewer_mvabund_species_pvalues")
}

}  # end if(FALSE)


# ============================================================================
# B. GENERALIZED DISSIMILARITY MODELING (GDM)
# ============================================================================
# Models nonlinear compositional turnover using I-splines. Tests whether
# the size-composition gradient saturates at large colony sizes.
#
# install.packages("gdm")
# ============================================================================

if (FALSE) {

cat("--- B. GDM ---\n\n")

if (requireNamespace("gdm", quietly = TRUE)) {
  library(gdm)

  comm_nz <- community_matrix[rowSums(community_matrix) > 0, ]
  coral_sub <- coral_master[match(rownames(comm_nz), coral_master$coral_id), ]

  # GDM requires site-pair table format
  # Environmental predictors: log_volume, site (as dummy)
  env_data <- data.frame(
    site = rownames(comm_nz),
    log_volume = log(coral_sub$volume),
    site_HAU = as.numeric(coral_sub$site == "HAU"),
    site_MAT = as.numeric(coral_sub$site == "MAT")
  )

  gdm_data <- formatsitepair(
    bioData = comm_nz,
    bioFormat = 1,  # community matrix
    predData = env_data,
    siteColumn = "site"
  )

  gdm_model <- gdm(gdm_data, geo = FALSE)

  if (!is.null(gdm_model)) {
    cat("  GDM deviance explained:", round(gdm_model$explained, 1), "%\n")
    cat("  Variable importance (I-spline max height):\n")
    print(gdm_model$sum.importance)

    # Plot I-splines to see nonlinearity
    plot(gdm_model, plot.layout = c(1, 3))

    save_table(data.frame(
      variable = names(gdm_model$sum.importance),
      importance = round(gdm_model$sum.importance, 4)
    ), "reviewer_gdm_variable_importance")
  }
}

}  # end if(FALSE)


# ============================================================================
# C. TITAN2 COMMUNITY BREAKPOINTS
# ============================================================================
# Change-point analysis to identify coral size thresholds where community
# composition shifts abruptly.
#
# install.packages("TITAN2")
# ============================================================================

if (FALSE) {

cat("--- C. TITAN2 ---\n\n")

if (requireNamespace("TITAN2", quietly = TRUE)) {
  library(TITAN2)

  comm_nz <- community_matrix[rowSums(community_matrix) > 0, ]
  comm_filt <- comm_nz[, colSums(comm_nz > 0) >= 5]
  coral_sub <- coral_master[match(rownames(comm_filt), coral_master$coral_id), ]

  # TITAN requires environmental gradient as vector
  set.seed(42)
  titan_result <- titan(coral_sub$volume, comm_filt,
                        minSplt = 5, numPerm = 250, boot = TRUE, nBoot = 500)

  cat("  Community change points detected:\n")
  cat("    Negative taxa (decline with size):", titan_result$sppmax[1], "\n")
  cat("    Positive taxa (increase with size):", titan_result$sppmax[2], "\n")

  plot(titan_result)

  save_table(titan_result$sppmax, "reviewer_titan2_breakpoints")
}

}  # end if(FALSE)


# ============================================================================
# D. BAYESIAN SEM (brms)
# ============================================================================
# Bayesian path model: volume → richness → condition, volume → abundance → condition
# Better uncertainty propagation than piecewiseSEM at n=84.
#
# install.packages("brms")
# ============================================================================

if (FALSE) {

cat("--- D. Bayesian SEM ---\n\n")

if (requireNamespace("brms", quietly = TRUE)) {
  library(brms)

  bef_data <- condition_scores %>%
    dplyr::select(coral_id, site, condition_score) %>%
    left_join(coral_master %>% dplyr::select(coral_id, volume, otu_richness, total_cafi),
              by = "coral_id") %>%
    filter(!is.na(condition_score), !is.na(otu_richness), !is.na(volume)) %>%
    mutate(
      log_volume_z = scale(log(volume))[,1],
      richness_z = scale(otu_richness)[,1],
      abundance_z = scale(sqrt(total_cafi))[,1],
      condition_z = scale(condition_score)[,1],
      site = factor(site)
    )

  # Path model as multivariate brms
  bf_rich <- bf(richness_z ~ log_volume_z + site)
  bf_abund <- bf(abundance_z ~ log_volume_z + site)
  bf_cond <- bf(condition_z ~ richness_z + abundance_z + log_volume_z + site)

  bayes_sem <- brm(
    bf_rich + bf_abund + bf_cond + set_rescor(FALSE),
    data = bef_data,
    chains = 4, iter = 4000, warmup = 1000, cores = 4,
    seed = 42
  )

  print(summary(bayes_sem))

  # Extract path coefficients
  post <- as_draws_df(bayes_sem)
  cat("\n  Posterior medians:\n")
  cat("    richness → condition:", round(median(post$b_conditionz_richness_z), 3), "\n")
  cat("    abundance → condition:", round(median(post$b_conditionz_abundance_z), 3), "\n")
}

}  # end if(FALSE)


# ============================================================================
# E. LOOCV OF RICHNESS-CONDITION MODEL
# ============================================================================
# Leave-one-out cross-validation to assess predictive stability of the
# richness → condition relationship (adj R² = 0.04).
# ============================================================================

if (FALSE) {

cat("--- E. LOOCV ---\n\n")

bef_data <- condition_scores %>%
  dplyr::select(coral_id, site, condition_score) %>%
  left_join(coral_master %>% dplyr::select(coral_id, volume, otu_richness),
            by = "coral_id") %>%
  filter(!is.na(condition_score), !is.na(otu_richness), !is.na(volume)) %>%
  mutate(log_volume = log(volume), site = factor(site))

n <- nrow(bef_data)
cv_preds <- numeric(n)

for (i in 1:n) {
  m_cv <- lm(condition_score ~ otu_richness + log_volume + site, data = bef_data[-i, ])
  cv_preds[i] <- predict(m_cv, newdata = bef_data[i, ])
}

# Cross-validated R²
ss_res <- sum((bef_data$condition_score - cv_preds)^2)
ss_tot <- sum((bef_data$condition_score - mean(bef_data$condition_score))^2)
cv_r2 <- 1 - ss_res / ss_tot

cat("  LOOCV R²:", round(cv_r2, 4), "\n")
cat("  Model adj R²:", round(summary(lm(condition_score ~ otu_richness + log_volume + site,
                                         data = bef_data))$adj.r.squared, 4), "\n")
cat("  Difference:", round(cv_r2 - summary(lm(condition_score ~ otu_richness + log_volume + site,
                                               data = bef_data))$adj.r.squared, 4), "\n")

if (cv_r2 > 0) {
  cat("  → Model has positive predictive power (not overfit)\n")
} else {
  cat("  → WARNING: Negative CV R² — model may be overfit\n")
}

save_table(data.frame(
  metric = c("LOOCV_R2", "adj_R2", "difference"),
  value = round(c(cv_r2,
                  summary(lm(condition_score ~ otu_richness + log_volume + site,
                             data = bef_data))$adj.r.squared,
                  cv_r2 - summary(lm(condition_score ~ otu_richness + log_volume + site,
                                     data = bef_data))$adj.r.squared), 4)
), "reviewer_loocv_bef")

}  # end if(FALSE)


# ============================================================================
# F. BAYESIAN MULTILEVEL SCALING (brms)
# ============================================================================
# Partial pooling across species for species-level scaling exponents.
# Regularizes estimates for rare species.
#
# install.packages("brms")
# ============================================================================

if (FALSE) {

cat("--- F. Bayesian multilevel scaling ---\n\n")

if (requireNamespace("brms", quietly = TRUE)) {
  library(brms)

  # Long-format species data
  comm_nz <- community_matrix[rowSums(community_matrix) > 0, ]
  # Species with ≥30 individuals and ≥15% prevalence
  sp_keep <- names(which(colSums(comm_nz > 0) >= 17 & colSums(comm_nz) >= 30))

  long_data <- lapply(sp_keep, function(sp) {
    coral_sub <- coral_master[match(rownames(comm_nz), coral_master$coral_id), ]
    data.frame(
      species = sp,
      abundance = comm_nz[, sp],
      log_volume = log(coral_sub$volume),
      site = coral_sub$site
    )
  }) %>% bind_rows()

  # Multilevel NB: random slopes for log_volume by species
  bayes_scaling <- brm(
    abundance ~ log_volume + site + (1 + log_volume | species),
    data = long_data,
    family = negbinomial(),
    chains = 4, iter = 4000, warmup = 1000, cores = 4,
    seed = 42
  )

  # Extract species-level exponents (fixed + random)
  ranef_df <- ranef(bayes_scaling)$species
  cat("  Species-level β estimates (shrinkage-regularized):\n")
  print(ranef_df[, , "log_volume"])
}

}  # end if(FALSE)


# ============================================================================
# G. CROSS-STUDY META-ANALYSIS OF SCALING EXPONENTS
# ============================================================================
# Fixed-effects meta-analysis pooling β estimates across the 3 CAFI studies.
# ============================================================================

if (FALSE) {

cat("--- G. Cross-study meta-analysis ---\n\n")

if (requireNamespace("metafor", quietly = TRUE)) {
  library(metafor)

  # Data from the three studies (update when Maatea paper finalizes)
  meta_data <- data.frame(
    study = c("Survey (this paper)", "Experiment (StierInReview)"),
    beta = c(0.523, 1.0),  # Update with actual experimental β
    se = c(0.047, 0.15),   # Update with actual SE
    n = c(114, 54),
    design = c("observational", "experimental")
  )

  # Random-effects meta-analysis
  rma_result <- rma(yi = beta, sei = se, data = meta_data, method = "REML")
  cat("  Pooled β:", round(rma_result$beta[1], 3),
      "[", round(rma_result$ci.lb, 3), ",", round(rma_result$ci.ub, 3), "]\n")
  cat("  Heterogeneity I²:", round(rma_result$I2, 1), "%\n")
  cat("  Q-test p:", format.pval(rma_result$QEp, 3), "\n\n")

  # Forest plot
  forest(rma_result, slab = meta_data$study,
         xlab = "Scaling exponent (β)",
         refline = 1, header = TRUE)

  save_table(data.frame(
    pooled_beta = round(rma_result$beta[1], 3),
    ci_lower = round(rma_result$ci.lb, 3),
    ci_upper = round(rma_result$ci.ub, 3),
    I2 = round(rma_result$I2, 1),
    Q_p = rma_result$QEp
  ), "reviewer_meta_analysis_scaling")
}

}  # end if(FALSE)


cat("\n")
cat("============================================================\n")
cat("    REVIEWER RESPONSE ANALYSES — TEMPLATE COMPLETE\n")
cat("============================================================\n")
cat("Set individual sections to TRUE to run as needed.\n\n")
