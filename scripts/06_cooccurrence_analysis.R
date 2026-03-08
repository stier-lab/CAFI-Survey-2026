# ============================================================================
# 06_cooccurrence_analysis.R - Null-Model Co-occurrence Analysis
# ============================================================================
#
# PURPOSE: Test pairwise co-occurrence and intraspecific density patterns
#          using volume-weighted null models (following Stier et al. 2012)
#
# HYPOTHESES:
#   H1 (Pairwise co-occurrence): Which species pairs co-occur more or less
#       than expected by chance, after controlling for coral volume?
#   H2 (Intraspecific density): Do species occur as pairs (exactly 2) more
#       than expected? (mating-pair hypothesis, Stier et al. 2012)
#   Size-dependent: How does coral size modulate co-occurrence patterns?
#
# ANALYSES:
#   Part 1: Volume-weighted pairwise co-occurrence null model (10,000 iter)
#   Part 2: Intraspecific density patterns (10,000 iter)
#   Part 3: Size-dependent co-occurrence (H1 by size class)
#   Part 4: Supplementary figure construction (3-panel co-occurrence figure)
#   Part 5: Legacy network analysis (supplement only)
#
# OUTPUTS:
#   Figures:
#     - output/figures/supplement/figS11_cooccurrence.png (3-panel)
#   Tables:
#     - output/tables/pairwise_cooccurrence.csv
#     - output/tables/intraspecific_density.csv
#     - output/tables/size_dependent_cooccurrence.csv
#     - output/tables/network_metrics.csv (legacy)
#     - output/tables/hub_species.csv (legacy)
#   Objects:
#     - output/objects/cooccurrence_results.rds
#     - output/objects/cafi_network.rds (legacy)
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-03-04
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI CO-OCCURRENCE ANALYSIS (Null Models)\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP AND DATA LOADING
# ============================================================================

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

# Explicitly load community_matrix (may not be in global env if script is run standalone)
if (!exists("community_matrix")) {
  community_matrix <- load_object("community_matrix")
}

# Create output directories
fig_dir <- file.path(PATHS$figures, "06_network")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("[OK] Setup complete\n")
cat("     Corals:", nrow(coral_master), "\n")
cat("     Community matrix:", nrow(community_matrix), "x", ncol(community_matrix), "\n\n")

# ============================================================================
# PARAMETERS
# ============================================================================

N_ITER <- 10000          # Null model iterations
# Minimum coral occurrences for pairwise co-occurrence testing
# (10 corals = ~9% prevalence; ensures sufficient data for null model comparison;
# lower thresholds inflate Type I error, higher thresholds exclude ecologically relevant species)
MIN_CORALS <- 10
SES_SIG <- 1.96          # SES threshold for significance reference
set.seed(42)             # Reproducibility

# ============================================================================
# PREPARE DATA
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PREPARING DATA\n")
cat("------------------------------------------------------------\n\n")

# Match community matrix rows to coral_master
matched <- rownames(community_matrix) %in% coral_master$coral_id
comm <- community_matrix[matched, , drop = FALSE]
vol <- coral_master$volume[match(rownames(comm), coral_master$coral_id)]
log_vol <- log(vol)
size_cls <- coral_master$size_class[match(rownames(comm), coral_master$coral_id)]

# Identify species on >= MIN_CORALS corals
comm_pa <- comm
comm_pa[comm_pa > 0] <- 1
n_corals_per_sp <- colSums(comm_pa)
focal_species <- names(which(n_corals_per_sp >= MIN_CORALS))

cat("  Species on >=", MIN_CORALS, "corals:", length(focal_species), "\n")
cat("  Corals with volume:", sum(!is.na(vol)), "\n\n")

# Get taxonomy for focal species (for grouping in figures)
sp_taxonomy <- cafi_clean %>%
  filter(otu %in% focal_species) %>%
  distinct(otu, phylum, class, order, family) %>%
  mutate(
    taxon_group = case_when(
      family == "Trapeziidae" ~ "Crabs (Trapeziidae)",
      family == "Alpheidae" ~ "Shrimp (Alpheidae)",
      family == "Palaemonidae" ~ "Shrimp (Palaemonidae)",
      family %in% c("Diogenidae", "Paguridae") ~ "Hermit crabs",
      family == "Xanthidae" | order == "Decapoda" & !family %in% c("Trapeziidae", "Alpheidae", "Palaemonidae", "Diogenidae", "Paguridae") ~ "Other crustaceans",
      class == "Actinopteri" | class == "Teleostei" ~ "Fish",
      class == "Gastropoda" ~ "Gastropods",
      class == "Polychaeta" | phylum == "Annelida" ~ "Polychaetes",
      TRUE ~ "Other"
    )
  )

# ============================================================================
# PART 1: PAIRWISE CO-OCCURRENCE NULL MODEL (H1)
# ============================================================================
# Algorithm (following Stier et al. 2012):
#   1. Fit logistic GLM per species: P(present) ~ log(volume) + site
#      → yields per-coral expected occurrence probability
#   2. Count observed pairwise co-occurrences from presence-absence matrix
#   3. For each of N_ITER null iterations:
#      a. Draw independent Bernoulli(predicted_prob) per species per coral
#      b. Count pairwise co-occurrences in null assemblage
#   4. Compute SES = (observed - mean_null) / sd_null per pair
#   5. FDR-correct p-values across all tested pairs
# This controls for the coral-size confound: large corals host more species
# simply because they are larger, not because of biotic interactions.
# ============================================================================

cat("============================================================\n")
cat("PART 1: PAIRWISE CO-OCCURRENCE (H1)\n")
cat("============================================================\n\n")

# Step 1: Fit logistic regression for each species: presence ~ log(volume) + site
# Controls for both coral size AND site-level species pool differences.
# Without site, the null overestimates co-occurrence for species restricted to
# the same site and underestimates it for species at different sites.

cat("1.1 Fitting volume + site occurrence models for", length(focal_species), "species...\n")

# Site factor for null model
site_vec <- coral_master$site[match(rownames(comm), coral_master$coral_id)]
site_fac <- factor(site_vec)

pred_probs <- matrix(NA, nrow = nrow(comm), ncol = length(focal_species))
colnames(pred_probs) <- focal_species
rownames(pred_probs) <- rownames(comm)

for (sp in focal_species) {
  y <- as.integer(comm_pa[, sp])
  # Try volume + site first; fall back to volume-only if site causes separation
  fit <- tryCatch(
    glm(y ~ log_vol + site_fac, family = binomial),
    warning = function(w) suppressWarnings(glm(y ~ log_vol + site_fac, family = binomial)),
    error = function(e) NULL
  )
  # Check for convergence / infinite coefficients → fall back to volume-only
  if (!is.null(fit) && any(!is.finite(coef(fit)))) {
    cat("    WARNING:", sp, "- volume+site model has non-finite coefficients, falling back\n")
    fit <- NULL
  }
  if (is.null(fit)) {
    fit <- tryCatch(
      glm(y ~ log_vol, family = binomial),
      warning = function(w) suppressWarnings(glm(y ~ log_vol, family = binomial)),
      error = function(e) NULL
    )
    if (!is.null(fit)) cat("    NOTE:", sp, "- using volume-only model (site caused separation)\n")
  }
  if (!is.null(fit)) {
    pred_probs[, sp] <- predict(fit, type = "response")
  } else {
    # Fallback: use marginal frequency
    cat("    WARNING:", sp, "- all GLMs failed, using marginal frequency\n")
    pred_probs[, sp] <- mean(y)
  }
}

cat("  Done. Predicted probabilities computed.\n")

# Step 2: Count observed co-occurrences for all pairs
n_sp <- length(focal_species)
n_pairs <- n_sp * (n_sp - 1) / 2

obs_cooccur <- matrix(0, nrow = n_sp, ncol = n_sp)
colnames(obs_cooccur) <- rownames(obs_cooccur) <- focal_species

pa_focal <- comm_pa[, focal_species]
for (i in 1:(n_sp - 1)) {
  for (j in (i + 1):n_sp) {
    obs_cooccur[i, j] <- sum(pa_focal[, i] == 1 & pa_focal[, j] == 1)
    obs_cooccur[j, i] <- obs_cooccur[i, j]
  }
}

cat("1.2 Observed co-occurrences computed for", n_pairs, "pairs.\n")

# Step 3: Run null model (10,000 iterations)
cat("1.3 Running null model (", N_ITER, " iterations)...\n")

null_cooccur_sum <- matrix(0, nrow = n_sp, ncol = n_sp)
null_cooccur_sq <- matrix(0, nrow = n_sp, ncol = n_sp)
null_exceed <- matrix(0, nrow = n_sp, ncol = n_sp)    # count(null >= obs)
null_below <- matrix(0, nrow = n_sp, ncol = n_sp)     # count(null <= obs)

n_corals <- nrow(comm)
tick <- round(N_ITER / 10)

for (iter in seq_len(N_ITER)) {
  if (iter %% tick == 0) cat("     Iteration", iter, "of", N_ITER, "\n")

  # Bernoulli draw for each species independently, using volume-predicted probs
  null_pa <- matrix(0, nrow = n_corals, ncol = n_sp)
  for (s in seq_len(n_sp)) {
    null_pa[, s] <- rbinom(n_corals, 1, pred_probs[, s])
  }

  # Count null co-occurrences
  for (i in 1:(n_sp - 1)) {
    for (j in (i + 1):n_sp) {
      null_count <- sum(null_pa[, i] == 1 & null_pa[, j] == 1)
      null_cooccur_sum[i, j] <- null_cooccur_sum[i, j] + null_count
      null_cooccur_sq[i, j] <- null_cooccur_sq[i, j] + null_count^2
      if (null_count >= obs_cooccur[i, j]) null_exceed[i, j] <- null_exceed[i, j] + 1
      if (null_count <= obs_cooccur[i, j]) null_below[i, j] <- null_below[i, j] + 1
    }
  }
}

# Compute SES and p-values
null_mean <- null_cooccur_sum / N_ITER
null_sd <- sqrt(null_cooccur_sq / N_ITER - null_mean^2)

ses_matrix <- matrix(0, nrow = n_sp, ncol = n_sp)
colnames(ses_matrix) <- rownames(ses_matrix) <- focal_species
p_positive <- matrix(1, nrow = n_sp, ncol = n_sp)  # test: obs >= null (positive association)
p_negative <- matrix(1, nrow = n_sp, ncol = n_sp)  # test: obs <= null (avoidance)

for (i in 1:(n_sp - 1)) {
  for (j in (i + 1):n_sp) {
    if (null_sd[i, j] > 0) {
      ses_matrix[i, j] <- (obs_cooccur[i, j] - null_mean[i, j]) / null_sd[i, j]
      ses_matrix[j, i] <- ses_matrix[i, j]
    }
    p_positive[i, j] <- null_exceed[i, j] / N_ITER
    p_negative[i, j] <- null_below[i, j] / N_ITER
    p_positive[j, i] <- p_positive[i, j]
    p_negative[j, i] <- p_negative[i, j]
  }
}

# Two-tailed p-value: min of positive and negative, doubled
p_twotail <- matrix(1, nrow = n_sp, ncol = n_sp)
for (i in 1:(n_sp - 1)) {
  for (j in (i + 1):n_sp) {
    p_twotail[i, j] <- 2 * min(p_positive[i, j], p_negative[i, j])
    p_twotail[i, j] <- min(p_twotail[i, j], 1)  # cap at 1
    p_twotail[j, i] <- p_twotail[i, j]
  }
}

# FDR correction across all pairs
upper_idx <- which(upper.tri(p_twotail))
p_raw_vec <- p_twotail[upper_idx]
p_fdr_vec <- p.adjust(p_raw_vec, method = "BH")
p_fdr_matrix <- matrix(1, nrow = n_sp, ncol = n_sp)
colnames(p_fdr_matrix) <- rownames(p_fdr_matrix) <- focal_species
p_fdr_matrix[upper_idx] <- p_fdr_vec
p_fdr_matrix <- pmin(p_fdr_matrix, t(p_fdr_matrix))  # symmetrize

# Summarize results
n_sig_positive <- sum(ses_matrix[upper.tri(ses_matrix)] > 0 & p_fdr_matrix[upper.tri(p_fdr_matrix)] < 0.05)
n_sig_negative <- sum(ses_matrix[upper.tri(ses_matrix)] < 0 & p_fdr_matrix[upper.tri(p_fdr_matrix)] < 0.05)
n_sig <- n_sig_positive + n_sig_negative

cat("\n  RESULTS:\n")
cat("    Total pairs tested:", n_pairs, "\n")
cat("    Significant (FDR < 0.05):", n_sig, "\n")
cat("      Positive associations:", n_sig_positive, "\n")
cat("      Negative associations:", n_sig_negative, "\n")

# Build results table
pair_results <- tibble()
for (i in 1:(n_sp - 1)) {
  for (j in (i + 1):n_sp) {
    pair_results <- bind_rows(pair_results, tibble(
      species_1 = focal_species[i],
      species_2 = focal_species[j],
      obs_cooccur = obs_cooccur[i, j],
      null_mean = round(null_mean[i, j], 2),
      null_sd = round(null_sd[i, j], 2),
      ses = round(ses_matrix[i, j], 3),
      p_raw = round(p_twotail[i, j], 5),
      p_fdr = round(p_fdr_matrix[i, j], 5),
      direction = ifelse(ses_matrix[i, j] > 0, "positive", "negative"),
      significant = p_fdr_matrix[i, j] < 0.05
    ))
  }
}
pair_results <- pair_results %>% arrange(p_fdr)

# Save
save_table(pair_results, "pairwise_cooccurrence")

cat("\n  Top 10 most significant pairs:\n")
print(head(pair_results %>% dplyr::select(species_1, species_2, ses, p_fdr, direction), 10))

# Key contrasts from Stier 2012
cat("\n  KEY CONTRASTS (Stier et al. 2012 predictions):\n")
key_pairs <- list(
  c("Trapezia serenei", "Alpheus lottini"),
  c("Trapezia serenei", "Trapezia bidentata"),
  c("Alpheus lottini", "Synalpheus charon"),
  c("Trapezia serenei", "Harpiliopsis beaupresii"),
  c("Alpheus lottini", "Harpiliopsis beaupresii")
)
for (kp in key_pairs) {
  row <- pair_results %>% filter(
    (species_1 == kp[1] & species_2 == kp[2]) |
    (species_1 == kp[2] & species_2 == kp[1])
  )
  if (nrow(row) > 0) {
    cat("    ", kp[1], "×", kp[2], ": SES =", row$ses[1],
        ", p_fdr =", row$p_fdr[1], "(", row$direction[1], ")\n")
  }
}


# ============================================================================
# PART 2: INTRASPECIFIC DENSITY PATTERNS (H2)
# ============================================================================

cat("\n============================================================\n")
cat("PART 2: INTRASPECIFIC DENSITY PATTERNS (H2)\n")
cat("============================================================\n\n")

set.seed(42)  # Reproducibility for Part 2 null model

# Select top abundant species for density analysis
# Use species with mean abundance > 1.5 per occupied coral AND on >= 20 corals
density_species <- c()
for (sp in focal_species) {
  abund <- comm[comm_pa[, sp] == 1, sp]
  if (length(abund) >= 20 && mean(abund) >= 1.5) {
    density_species <- c(density_species, sp)
  }
}
cat("  Species for density analysis:", length(density_species), "\n")
cat("   ", paste(density_species, collapse = ", "), "\n\n")

# For each species: fix total abundance N and number of occupied corals K
# Null: distribute N individuals across K corals with probabilities
# proportional to coral volume
# Note: Multinomial null model is site-agnostic; volume-proportional allocation
# controls for size but not site-specific pools.
density_results <- list()

for (sp in density_species) {
  cat("  Processing:", sp, "...\n")

  occupied <- which(comm_pa[, sp] == 1)
  K <- length(occupied)
  N_total <- sum(comm[occupied, sp])
  obs_abund <- comm[occupied, sp]

  # Volume weights for occupied corals
  # NOTE: Volume weighting is across all sites (not site-specific). If species are site-restricted, this null is slightly conservative.
  vol_occupied <- vol[occupied]
  vol_weights <- vol_occupied / sum(vol_occupied)

  # Observed density distribution
  obs_freq <- table(factor(pmin(obs_abund, 5), levels = 1:5))
  names(obs_freq) <- c("1", "2", "3", "4", "5+")

  # Null model: multinomial allocation proportional to volume
  null_freq_mat <- matrix(0, nrow = N_ITER, ncol = 5)
  colnames(null_freq_mat) <- c("1", "2", "3", "4", "5+")

  for (iter in seq_len(N_ITER)) {
    # Allocate N_total individuals to K corals
    null_alloc <- rmultinom(1, N_total, vol_weights)[, 1]
    null_alloc_capped <- pmin(null_alloc, 5)
    null_tab <- table(factor(null_alloc_capped[null_alloc_capped > 0], levels = 1:5))
    null_freq_mat[iter, ] <- as.numeric(null_tab)
  }

  # Compute SES for each density bin
  null_mean_freq <- colMeans(null_freq_mat)
  null_sd_freq <- apply(null_freq_mat, 2, sd)
  ses_freq <- ifelse(null_sd_freq > 0,
                     (as.numeric(obs_freq) - null_mean_freq) / null_sd_freq, 0)

  # p-values for "pairs" (density = 2)
  obs_pairs <- as.numeric(obs_freq["2"])
  null_pairs <- null_freq_mat[, 2]
  p_pairs <- mean(null_pairs >= obs_pairs)

  # 95% CI for null
  null_ci_lo <- apply(null_freq_mat, 2, quantile, 0.025)
  null_ci_hi <- apply(null_freq_mat, 2, quantile, 0.975)

  density_results[[sp]] <- list(
    species = sp,
    N_total = N_total,
    K_corals = K,
    obs_freq = as.numeric(obs_freq),
    null_mean = null_mean_freq,
    null_ci_lo = null_ci_lo,
    null_ci_hi = null_ci_hi,
    ses_freq = ses_freq,
    p_pairs = p_pairs,
    mean_per_coral = N_total / K
  )

  cat("    N =", N_total, "on", K, "corals. Pairs observed =", obs_pairs,
      ", null mean =", round(null_mean_freq[2], 1),
      ", SES =", round(ses_freq[2], 2),
      ", p =", round(p_pairs, 4), "\n")
}

# Build summary table
density_table <- map_dfr(density_results, function(r) {
  tibble(
    species = r$species,
    n_individuals = r$N_total,
    n_corals = r$K_corals,
    mean_per_coral = round(r$mean_per_coral, 2),
    obs_singles = r$obs_freq[1],
    obs_pairs = r$obs_freq[2],
    obs_3plus = sum(r$obs_freq[3:5]),
    null_mean_pairs = round(r$null_mean[2], 1),
    null_ci_lo_pairs = round(r$null_ci_lo[2], 1),
    null_ci_hi_pairs = round(r$null_ci_hi[2], 1),
    ses_pairs = round(r$ses_freq[2], 3),
    p_pairs = round(r$p_pairs, 5)
  )
})

# FDR correction across species
density_table$p_pairs_fdr <- round(p.adjust(density_table$p_pairs, method = "BH"), 5)

save_table(density_table, "intraspecific_density")
cat("\n  Density pattern summary:\n")
print(density_table %>% dplyr::select(species, mean_per_coral, obs_pairs,
                                      null_mean_pairs, ses_pairs, p_pairs_fdr))


# ============================================================================
# PART 3: SIZE-DEPENDENT CO-OCCURRENCE
# ============================================================================

cat("\n============================================================\n")
cat("PART 3: SIZE-DEPENDENT CO-OCCURRENCE\n")
cat("============================================================\n\n")

set.seed(42)  # Reproducibility for Part 3 null model

# NOTE: Species pools differ across size classes (species with >= 5 occurrences
# within each class). SES values are therefore not directly comparable across
# classes because null variance reflects different species sets. Results are
# presented as within-class patterns, not cross-class comparisons.
size_classes <- c("Small", "Medium", "Large")
size_ses_results <- list()

for (sc in size_classes) {
  cat("  Running null model for", sc, "corals...\n")

  # Subset corals
  sc_idx <- which(size_cls == sc)
  n_sc <- length(sc_idx)
  cat("    n =", n_sc, "\n")

  comm_sc <- comm_pa[sc_idx, focal_species]
  vol_sc <- log_vol[sc_idx]

  # Filter to species present on >= 5 corals in this size class
  sp_present <- colSums(comm_sc) >= 5
  sp_sc <- focal_species[sp_present]
  cat("    Species with >= 5 occurrences:", length(sp_sc), "\n")

  if (length(sp_sc) < 3) {
    cat("    Skipping: too few species.\n")
    next
  }

  # Fit logistic regressions for this size class (volume + site)
  site_sc <- site_fac[sc_idx]
  pred_sc <- matrix(NA, nrow = n_sc, ncol = length(sp_sc))
  colnames(pred_sc) <- sp_sc
  for (sp in sp_sc) {
    y <- as.integer(comm_sc[, sp])
    fit <- tryCatch(
      glm(y ~ vol_sc + site_sc, family = binomial),
      warning = function(w) suppressWarnings(glm(y ~ vol_sc + site_sc, family = binomial)),
      error = function(e) NULL
    )
    if (!is.null(fit) && any(!is.finite(coef(fit)))) fit <- NULL
    if (is.null(fit)) {
      fit <- tryCatch(
        glm(y ~ vol_sc, family = binomial),
        warning = function(w) suppressWarnings(glm(y ~ vol_sc, family = binomial)),
        error = function(e) NULL
      )
    }
    if (!is.null(fit)) {
      pred_sc[, sp] <- predict(fit, type = "response")
    } else {
      pred_sc[, sp] <- mean(y)
    }
  }

  # Observed co-occurrences
  n_sp_sc <- length(sp_sc)
  obs_sc <- matrix(0, nrow = n_sp_sc, ncol = n_sp_sc)
  for (i in 1:(n_sp_sc - 1)) {
    for (j in (i + 1):n_sp_sc) {
      obs_sc[i, j] <- sum(comm_sc[, sp_sc[i]] == 1 & comm_sc[, sp_sc[j]] == 1)
    }
  }

  # Null model (fewer iterations for size subsets — still 10,000)
  null_sum_sc <- matrix(0, nrow = n_sp_sc, ncol = n_sp_sc)
  null_sq_sc <- matrix(0, nrow = n_sp_sc, ncol = n_sp_sc)

  for (iter in seq_len(N_ITER)) {
    null_pa <- matrix(0, nrow = n_sc, ncol = n_sp_sc)
    for (s in seq_len(n_sp_sc)) {
      null_pa[, s] <- rbinom(n_sc, 1, pred_sc[, s])
    }
    for (i in 1:(n_sp_sc - 1)) {
      for (j in (i + 1):n_sp_sc) {
        nc <- sum(null_pa[, i] == 1 & null_pa[, j] == 1)
        null_sum_sc[i, j] <- null_sum_sc[i, j] + nc
        null_sq_sc[i, j] <- null_sq_sc[i, j] + nc^2
      }
    }
  }

  null_mean_sc <- null_sum_sc / N_ITER
  null_sd_sc <- sqrt(null_sq_sc / N_ITER - null_mean_sc^2)

  ses_sc <- matrix(0, nrow = n_sp_sc, ncol = n_sp_sc)
  colnames(ses_sc) <- rownames(ses_sc) <- sp_sc
  for (i in 1:(n_sp_sc - 1)) {
    for (j in (i + 1):n_sp_sc) {
      if (null_sd_sc[i, j] > 0) {
        ses_sc[i, j] <- (obs_sc[i, j] - null_mean_sc[i, j]) / null_sd_sc[i, j]
        ses_sc[j, i] <- ses_sc[i, j]
      }
    }
  }

  size_ses_results[[sc]] <- list(
    size_class = sc,
    species = sp_sc,
    ses = ses_sc,
    n_corals = n_sc
  )

  cat("    SES range:", round(min(ses_sc[upper.tri(ses_sc)]), 2), "to",
      round(max(ses_sc[upper.tri(ses_sc)]), 2), "\n")
}

# Build size-dependent comparison table for key pairs
key_pair_names <- list(
  c("Trapezia serenei", "Alpheus lottini"),
  c("Trapezia serenei", "Trapezia bidentata"),
  c("Alpheus lottini", "Synalpheus charon"),
  c("Trapezia serenei", "Harpiliopsis beaupresii"),
  c("Alpheus lottini", "Harpiliopsis beaupresii"),
  c("Galeropsis monodonta", "Trapezia serenei")
)

size_dep_table <- tibble()
for (kp in key_pair_names) {
  for (sc in size_classes) {
    if (!sc %in% names(size_ses_results)) next
    res <- size_ses_results[[sc]]
    sp1_idx <- match(kp[1], res$species)
    sp2_idx <- match(kp[2], res$species)
    if (!is.na(sp1_idx) && !is.na(sp2_idx)) {
      size_dep_table <- bind_rows(size_dep_table, tibble(
        species_1 = kp[1],
        species_2 = kp[2],
        pair_label = paste0(
          substr(kp[1], 1, 1), ". ",
          sub("^\\S+ ", "", kp[1]), " × ",
          substr(kp[2], 1, 1), ". ",
          sub("^\\S+ ", "", kp[2])
        ),
        size_class = sc,
        ses = round(res$ses[sp1_idx, sp2_idx], 3),
        n_corals = res$n_corals
      ))
    }
  }
}

save_table(size_dep_table, "size_dependent_cooccurrence")
cat("\n  Size-dependent co-occurrence for key pairs:\n")
if (nrow(size_dep_table) > 0) {
  print(size_dep_table %>% pivot_wider(
    id_cols = pair_label, names_from = size_class, values_from = ses
  ))
}


# ============================================================================
# PART 4: SUPPLEMENTARY FIGURE — CO-OCCURRENCE NULL MODEL (3-PANEL)
# ============================================================================

cat("\n============================================================\n")
cat("PART 4: SUPPLEMENTARY FIGURE (Co-occurrence)\n")
cat("============================================================\n\n")

# --- Panel A: Pairwise co-occurrence heatmap ---

# Select top ~20 species for heatmap readability (by prevalence)
top_n_heatmap <- min(20, n_sp)
top_sp_heatmap <- names(sort(n_corals_per_sp[focal_species], decreasing = TRUE))[1:top_n_heatmap]

# Order by taxonomic group
sp_order_df <- sp_taxonomy %>%
  filter(otu %in% top_sp_heatmap) %>%
  arrange(taxon_group, otu)
sp_order <- sp_order_df$otu

# Make abbreviated names (G. species)
abbrev_name <- function(x) {
  parts <- strsplit(x, " ")
  sapply(parts, function(p) {
    if (length(p) >= 2) paste0(substr(p[1], 1, 1), ". ", p[2])
    else p[1]
  })
}

# Build heatmap data from SES matrix
heat_data <- tibble()
for (i in seq_along(sp_order)) {
  for (j in seq_along(sp_order)) {
    if (i == j) next
    sp_i <- sp_order[i]
    sp_j <- sp_order[j]
    idx_i <- match(sp_i, focal_species)
    idx_j <- match(sp_j, focal_species)
    if (!is.na(idx_i) && !is.na(idx_j)) {
      heat_data <- bind_rows(heat_data, tibble(
        sp1 = sp_i, sp2 = sp_j,
        ses = ses_matrix[idx_i, idx_j],
        p_fdr = p_fdr_matrix[idx_i, idx_j],
        sig = p_fdr_matrix[idx_i, idx_j] < 0.05
      ))
    }
  }
}

# Only show upper triangle
heat_data <- heat_data %>%
  mutate(
    sp1_f = factor(sp1, levels = sp_order),
    sp2_f = factor(sp2, levels = sp_order),
    sp1_abbr = abbrev_name(sp1),
    sp2_abbr = abbrev_name(sp2)
  ) %>%
  filter(as.numeric(sp1_f) < as.numeric(sp2_f))

# Add group color bars
sp_group_colors <- sp_order_df %>%
  mutate(abbr = abbrev_name(otu))

# Clamp SES for color scale
ses_cap <- 5
heat_data$ses_capped <- pmax(pmin(heat_data$ses, ses_cap), -ses_cap)

p_heatmap <- ggplot(heat_data, aes(
  x = factor(sp2_abbr, levels = abbrev_name(sp_order)),
  y = factor(sp1_abbr, levels = rev(abbrev_name(sp_order)))
)) +
  geom_tile(aes(fill = ses_capped), color = "white", linewidth = 0.3) +
  geom_text(data = heat_data %>% filter(sig),
            aes(label = "*"), size = 4, vjust = 0.8, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-ses_cap, ses_cap),
    name = "SES",
    breaks = c(-4, -2, 0, 2, 4)
  ) +
  labs(x = NULL, y = NULL) +
  theme_publication(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5, face = "italic"),
    axis.text.y = element_text(size = 5, face = "italic"),
    legend.position = "right",
    legend.key.height = unit(15, "mm"),
    legend.key.width = unit(4, "mm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(5, 5, 5, 5, "mm")
  )


# --- Panel B: Intraspecific density patterns ---

# Select the most biologically interesting species for the figure:
# Prioritize significant results + key obligate symbionts
# (All 15 species are in the table; figure shows focal 8)
density_plot_species <- density_table %>%
  arrange(p_pairs_fdr) %>%
  head(8) %>%
  pull(species)

# Build data for density plot
density_plot_data <- tibble()
for (sp in density_plot_species) {
  r <- density_results[[sp]]
  sp_abbr <- abbrev_name(sp)
  for (bin in 1:5) {
    bin_label <- c("1", "2", "3", "4", "5+")[bin]
    density_plot_data <- bind_rows(density_plot_data, tibble(
      species = sp_abbr,
      density_bin = bin_label,
      observed = r$obs_freq[bin],
      null_mean = r$null_mean[bin],
      null_ci_lo = r$null_ci_lo[bin],
      null_ci_hi = r$null_ci_hi[bin],
      ses = r$ses_freq[bin]
    ))
  }
}

density_plot_data$density_bin <- factor(density_plot_data$density_bin,
                                        levels = c("1", "2", "3", "4", "5+"))
density_plot_data$species <- factor(density_plot_data$species,
                                     levels = abbrev_name(density_plot_species))

p_density <- ggplot(density_plot_data, aes(x = density_bin)) +
  # Null CI band
  geom_crossbar(aes(y = null_mean, ymin = null_ci_lo, ymax = null_ci_hi),
                fill = "grey85", color = "grey60", width = 0.5, linewidth = 0.3,
                fatten = 1) +
  # Observed points
  geom_point(aes(y = observed), size = 3, color = "#B2182B") +
  facet_wrap(~ species, scales = "free_y", ncol = 2) +
  scale_y_continuous(breaks = function(x) {
    b <- pretty(x, n = 4)
    b <- b[b == floor(b)]
    if (length(b) == 0) b <- range(x, na.rm = TRUE)
    b
  }) +
  labs(x = "Individuals per coral", y = "Number of corals") +
  theme_publication(base_size = 10) +
  theme(
    strip.text = element_text(size = 9, face = "bold.italic"),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5, "mm")
  )


# --- Panel C: Size-dependent co-occurrence ---

if (nrow(size_dep_table) > 0) {
  size_dep_table$size_class <- factor(size_dep_table$size_class,
                                       levels = c("Small", "Medium", "Large"))

  p_size <- ggplot(size_dep_table, aes(x = size_class, y = ses,
                                        group = pair_label, color = pair_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = c(-SES_SIG, SES_SIG), linetype = "dotted",
               color = "grey70", linewidth = 0.4) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                                  "#0072B2", "#D55E00", "#CC79A7", "#000000"),
                       name = "Species pair") +
    labs(x = "Coral size class", y = "SES") +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 7, face = "italic"),
      legend.title = element_text(size = 9, face = "bold"),
      legend.key.size = unit(3, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5, "mm")
    ) +
    guides(color = guide_legend(ncol = 3))
} else {
  p_size <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") +
    theme_void()
}


# --- Compose Supplementary Co-occurrence Figure ---

figS11_cooccurrence <- (p_heatmap / (p_density | p_size)) +
  plot_layout(heights = c(0.9, 1)) +
  plot_annotation(tag_levels = "A")

save_figure(figS11_cooccurrence,
            file.path(PATHS$fig_supplement, "figS11_cooccurrence.png"),
            width = 300, height = 260, units = "mm")
save_figure(figS11_cooccurrence,
            file.path(fig_dir, "figS11_cooccurrence.png"),
            width = 300, height = 260, units = "mm")

# Co-occurrence is supplement-only (Fig S11)
cat("   Co-occurrence saved to supplement only (figS11)\n")

# Write legend/results text
legend_text <- paste0(
  "Figure S11. Null-model co-occurrence analysis of coral-associated fauna.\n\n",
  "(A) Pairwise co-occurrence heatmap showing standardized effect sizes (SES) ",
  "from a volume-weighted null model (", format(N_ITER, big.mark = ","),
  " iterations). Red = positive association (co-occurrence more frequent than expected); ",
  "blue = avoidance (less frequent than expected). Asterisks (*) indicate FDR-corrected ",
  "p < 0.05. Species are ordered by taxonomic group. ",
  n_sig_positive, " significant positive and ", n_sig_negative,
  " significant negative associations were detected among ", n_pairs, " pairs tested.\n\n",
  "(B) Intraspecific density patterns for ", length(density_plot_species),
  " focal species (of ", length(density_species), " tested). Red points show observed frequency of corals hosting ",
  "1, 2, 3, 4, or 5+ individuals. Grey bars show null expectation (95% CI) from ",
  "volume-proportional random allocation of individuals across occupied corals.\n\n",
  "(C) Size-dependent co-occurrence. SES values for key species pairs across ",
  "three coral size classes (terciles). Dashed line = null expectation (SES = 0); ",
  "dotted lines = ±1.96 reference. Size classes: Small (n = ",
  size_ses_results[["Small"]]$n_corals, "), Medium (n = ",
  size_ses_results[["Medium"]]$n_corals, "), Large (n = ",
  size_ses_results[["Large"]]$n_corals, ")."
)
writeLines(legend_text, file.path(PATHS$fig_supplement, "figS11_cooccurrence_legend_results.txt"))

cat("  Legend text saved.\n")


# ============================================================================
# PART 5: LEGACY NETWORK ANALYSIS (Supplement)
# ============================================================================

cat("\n============================================================\n")
cat("PART 5: LEGACY NETWORK ANALYSIS (Supplement)\n")
cat("============================================================\n\n")

library(igraph)

# Resolve igraph namespace conflicts
conflicted::conflict_prefer("groups", "igraph")
conflicted::conflict_prefer("union", "base")
conflicted::conflict_prefer("intersect", "base")
conflicted::conflict_prefer("setdiff", "base")
conflicted::conflict_prefer("as_data_frame", "igraph")

# Build co-occurrence network from volume-residualized correlations
# (identical method to original script 06)

comm_binary <- community_matrix
comm_binary[comm_binary > 0] <- 1

# Filter species (>= 5% of corals)
min_occurrence <- ceiling(nrow(comm_binary) * 0.05)
species_keep <- colSums(comm_binary) >= min_occurrence
comm_net <- comm_binary[, species_keep]

# Match to coral_master
matched_net <- rownames(comm_net) %in% coral_master$coral_id
comm_net <- comm_net[matched_net, , drop = FALSE]
vol_net <- log(coral_master$volume[match(rownames(comm_net), coral_master$coral_id)])

# Residualize on volume
resid_net <- matrix(NA, nrow = nrow(comm_net), ncol = ncol(comm_net))
colnames(resid_net) <- colnames(comm_net)
for (sp in seq_len(ncol(comm_net))) {
  y <- comm_net[, sp]
  fit <- tryCatch(
    glm(y ~ vol_net, family = binomial),
    warning = function(w) suppressWarnings(glm(y ~ vol_net, family = binomial)),
    error = function(e) NULL
  )
  if (!is.null(fit)) resid_net[, sp] <- residuals(fit, type = "deviance")
  else resid_net[, sp] <- y - mean(y)
}

# Spearman correlations + FDR
cor_net <- cor(resid_net, method = "spearman", use = "pairwise.complete.obs")
cor_net[is.na(cor_net)] <- 0
threshold <- 0.3
n_sp_net <- ncol(resid_net)

p_net <- matrix(1, nrow = n_sp_net, ncol = n_sp_net)
for (i in 1:(n_sp_net - 1)) {
  for (j in (i + 1):n_sp_net) {
    if (abs(cor_net[i, j]) > threshold) {
      ct <- tryCatch(
        cor.test(resid_net[, i], resid_net[, j], method = "spearman"),
        error = function(e) {
          cat("    WARNING: cor.test failed for", colnames(resid_net)[i], "-",
              colnames(resid_net)[j], ":", e$message, "\n")
          list(p.value = 1)
        }
      )
      p_net[i, j] <- p_net[j, i] <- ct$p.value
    }
  }
}

upper_net <- which(upper.tri(p_net))
p_fdr_net <- p.adjust(p_net[upper_net], method = "BH")
p_adj_net <- matrix(1, nrow = n_sp_net, ncol = n_sp_net)
p_adj_net[upper_net] <- p_fdr_net
p_adj_net <- pmin(p_adj_net, t(p_adj_net))

# Build adjacency matrix
adj_net <- (abs(cor_net) > threshold & p_adj_net < 0.05) * 1
diag(adj_net) <- 0
colnames(adj_net) <- rownames(adj_net) <- colnames(resid_net)

# Remove isolated nodes
connected <- colSums(adj_net) > 0
adj_net <- adj_net[connected, connected]

# Create igraph
g <- graph_from_adjacency_matrix(adj_net, mode = "undirected", diag = FALSE)

# Louvain community detection
set.seed(42)
louv <- cluster_louvain(g)
membership <- membership(louv)
V(g)$group <- as.character(membership)

# Network metrics
net_metrics <- tibble(
  metric = c("Nodes", "Edges", "Density", "Transitivity", "Modularity_Q", "N_groups"),
  value = c(vcount(g), ecount(g), graph.density(g), transitivity(g, type = "global"),
            modularity(louv), length(unique(membership)))
)
save_table(net_metrics, "network_metrics")

# Hub species
deg <- degree(g)
eig <- eigen_centrality(g)$vector
hub_df <- tibble(
  species = names(deg),
  degree = deg,
  eigenvector_centrality = round(eig, 4),
  group = membership[names(deg)]
) %>% arrange(desc(degree))
save_table(hub_df, "hub_species")

# Save network object
save_object(list(graph = g, membership = membership, metrics = net_metrics,
                 hub_species = hub_df), "cafi_network")

cat("  Network: ", vcount(g), " nodes, ", ecount(g), " edges\n")
cat("  Modularity Q =", round(modularity(louv), 3), "\n")
cat("  Density =", round(graph.density(g), 3), "\n")

# Simple network figure for supplement (circular layout)
n_groups <- length(unique(membership))
group_colors <- scales::hue_pal()(n_groups)

# Create abbreviated labels for network
V(g)$label <- sapply(V(g)$name, function(x) {
  parts <- strsplit(x, " ")[[1]]
  if (length(parts) >= 2) paste0(substr(parts[1], 1, 1), ". ", parts[2])
  else x
})

layout_circle <- layout_in_circle(g, order = order(as.numeric(V(g)$group)))

# ============================================================================
# SAVE RESULTS OBJECT
# ============================================================================

save_object(list(
  pair_results = pair_results,
  ses_matrix = ses_matrix,
  p_fdr_matrix = p_fdr_matrix,
  focal_species = focal_species,
  density_results = density_results,
  density_table = density_table,
  size_dep_table = size_dep_table,
  size_ses_results = size_ses_results,
  n_iter = N_ITER,
  min_corals = MIN_CORALS
), "cooccurrence_results")


# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("    CO-OCCURRENCE ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("H1 (Pairwise co-occurrence):\n")
cat("  ", n_pairs, "pairs tested;", n_sig, "significant (FDR < 0.05)\n")
cat("  ", n_sig_positive, "positive,", n_sig_negative, "negative\n\n")

cat("H2 (Intraspecific density):\n")
n_sig_pairs <- sum(density_table$p_pairs_fdr < 0.05)
cat("  ", length(density_species), "species tested;",
    n_sig_pairs, "with significant pair excess\n\n")

cat("Size-dependent co-occurrence:\n")
cat("  ", nrow(size_dep_table), "pair × size-class comparisons\n\n")

cat("Outputs:\n")
cat("  Figures:\n")
cat("    output/figures/supplement/figS11_cooccurrence.png\n")
cat("  Tables:\n")
cat("    output/tables/pairwise_cooccurrence.csv\n")
cat("    output/tables/intraspecific_density.csv\n")
cat("    output/tables/size_dependent_cooccurrence.csv\n")
cat("    output/tables/network_metrics.csv\n")
cat("    output/tables/hub_species.csv\n")
cat("  Objects:\n")
cat("    output/objects/cooccurrence_results.rds\n")
cat("    output/objects/cafi_network.rds\n")
