# ============================================================================
# network_abundance_artifact_test.R - Test if Network Transitivity is an
#                                     Abundance Artifact
# ============================================================================
#
# PURPOSE: Investigate whether the strong network transitivity (z=36.1, 1.66x
#          null) in CAFI co-occurrence networks is driven by species abundance
#          patterns or reflects true ecological associations.
#
# BACKGROUND:
#   The original network analysis found highly non-random clustering
#   (transitivity = 0.94 vs null = 0.56). However, abundant species may
#   appear to co-occur simply because they're everywhere, creating
#   spurious correlations.
#
#   The original script (06_network_analysis.R) controls for coral VOLUME
#   using deviance residuals from logistic GLM. But this does not control
#   for species-level ABUNDANCE confounds.
#
# ANALYSES:
#   1. Test correlation between species abundance and node centrality
#   2. Create presence/absence network (remove abundance entirely)
#   3. Create partial correlation network (regress out log-abundance)
#   4. Re-test transitivity with appropriate null models
#   5. Identify which species pairs drive clustering patterns
#
# OUTPUTS:
#   - output/tables/network_abundance_artifact.csv
#   - output/network_transitivity_artifact_analysis.txt
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-28
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    NETWORK TRANSITIVITY ABUNDANCE ARTIFACT TEST\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

# Load setup (packages, paths, theme)
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

# Load igraph for network analysis
library(igraph)

# Load original network results
network_results <- load_object("cafi_network")

cat("[OK] Setup complete\n")
cat("     Corals:", nrow(coral_master), "\n")
cat("     Community matrix:", nrow(community_matrix), "x", ncol(community_matrix), "\n\n")

# ============================================================================
# PART 1: EXAMINE CURRENT NETWORK CONSTRUCTION
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 1: CURRENT NETWORK STRUCTURE\n")
cat("------------------------------------------------------------\n\n")

# Extract original network
g_original <- network_results$graph
n_nodes <- vcount(g_original)
n_edges <- ecount(g_original)
transitivity_original <- transitivity(g_original, type = "global")

cat("Original Network (from script 06):\n")
cat("  - Nodes (species):", n_nodes, "\n")
cat("  - Edges:", n_edges, "\n")
cat("  - Network density:", round(edge_density(g_original), 4), "\n")
cat("  - Transitivity (global):", round(transitivity_original, 4), "\n")
cat("  - Original null comparison: z = 36.1, transitivity 1.66x null\n\n")

cat("IMPORTANT: The original script (06_network_analysis.R) controls for CORAL VOLUME\n")
cat("using deviance residuals from logistic GLM: P(species_i present) ~ log(volume)\n")
cat("This removes the 'larger corals host more species' confound.\n")
cat("BUT it does NOT control for species-level abundance effects:\n")
cat("  - Common species appear on many corals -> high co-occurrence with everything\n")
cat("  - This creates spurious positive correlations\n\n")

# ============================================================================
# PART 2: TEST ABUNDANCE CONFOUNDING
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 2: ABUNDANCE CONFOUNDING ANALYSIS\n")
cat("------------------------------------------------------------\n\n")

# 2.1 Calculate species-level abundance metrics
cat("2.1 Calculating species abundance metrics...\n\n")

# Get species in the network
network_species <- V(g_original)$name

# Calculate abundance metrics for network species
species_metrics <- data.frame(
  species = network_species,
  total_abundance = sapply(network_species, function(sp) {
    if (sp %in% colnames(community_matrix)) sum(community_matrix[, sp], na.rm = TRUE) else 0
  }),
  occurrence = sapply(network_species, function(sp) {
    if (sp %in% colnames(community_matrix)) sum(community_matrix[, sp] > 0, na.rm = TRUE) else 0
  }),
  stringsAsFactors = FALSE
)

species_metrics$log_abundance <- log10(species_metrics$total_abundance + 1)
species_metrics$occurrence_pct <- species_metrics$occurrence / nrow(community_matrix) * 100

# Add network centrality metrics
species_metrics$degree <- degree(g_original)[match(species_metrics$species, V(g_original)$name)]
species_metrics$eigenvector <- eigen_centrality(g_original)$vector[match(species_metrics$species, V(g_original)$name)]
species_metrics$local_transitivity <- transitivity(g_original, type = "local")[match(species_metrics$species, V(g_original)$name)]

cat("  Summary of species in network:\n")
cat("    - Median total abundance:", median(species_metrics$total_abundance), "\n")
cat("    - Median occurrence:", median(species_metrics$occurrence), "corals\n")
cat("    - Mean degree:", round(mean(species_metrics$degree), 2), "\n\n")

# 2.2 Test correlation: Abundance vs Node Degree
cat("2.2 Testing abundance-centrality correlations...\n\n")

cor_abund_degree <- cor.test(species_metrics$log_abundance, species_metrics$degree,
                              method = "spearman")
cor_occur_degree <- cor.test(species_metrics$occurrence, species_metrics$degree,
                              method = "spearman")
cor_abund_transitivity <- cor.test(species_metrics$log_abundance,
                                    species_metrics$local_transitivity,
                                    method = "spearman", use = "complete.obs")

cat("  Abundance vs Degree:\n")
cat("    - Spearman rho =", round(cor_abund_degree$estimate, 3), "\n")
cat("    - p-value =", format(cor_abund_degree$p.value, scientific = TRUE, digits = 3), "\n\n")

cat("  Occurrence vs Degree:\n")
cat("    - Spearman rho =", round(cor_occur_degree$estimate, 3), "\n")
cat("    - p-value =", format(cor_occur_degree$p.value, scientific = TRUE, digits = 3), "\n\n")

cat("  Abundance vs Local Transitivity:\n")
cat("    - Spearman rho =", round(cor_abund_transitivity$estimate, 3), "\n")
cat("    - p-value =", format(cor_abund_transitivity$p.value, scientific = TRUE, digits = 3), "\n\n")

# Interpretation
if (cor_abund_degree$p.value < 0.05 && cor_abund_degree$estimate > 0.3) {
  cat("  FINDING: Strong positive correlation between abundance and degree.\n")
  cat("           This suggests abundant species may have artificially high connectivity.\n\n")
} else {
  cat("  FINDING: Weak or no correlation between abundance and degree.\n")
  cat("           Abundance may not be driving network structure.\n\n")
}

# ============================================================================
# PART 3: ABUNDANCE-CONTROLLED NETWORK - PRESENCE/ABSENCE
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 3: PRESENCE/ABSENCE NETWORK (No Abundance Information)\n")
cat("------------------------------------------------------------\n\n")

# The original network is already built from presence/absence after volume
# residualization. Let's re-build using raw presence/absence without the
# volume correction to see how much the volume correction matters.

cat("3.1 Building raw presence/absence co-occurrence network...\n")

# Convert to binary presence/absence
comm_binary <- community_matrix
comm_binary[comm_binary > 0] <- 1

# Filter to same species as original network (>= 5% occurrence)
min_occurrence <- ceiling(nrow(comm_binary) * 0.05)
species_keep <- colSums(comm_binary) >= min_occurrence
comm_filtered_pa <- comm_binary[, species_keep]

cat("  Species retained (>= 5% occurrence):", ncol(comm_filtered_pa), "\n")

# Calculate Jaccard similarity for presence/absence
# Jaccard = intersection / union (proportion of shared presences)
n_sp <- ncol(comm_filtered_pa)
jaccard_matrix <- matrix(0, nrow = n_sp, ncol = n_sp)
colnames(jaccard_matrix) <- rownames(jaccard_matrix) <- colnames(comm_filtered_pa)

for (i in 1:(n_sp - 1)) {
  for (j in (i + 1):n_sp) {
    sp_i <- comm_filtered_pa[, i]
    sp_j <- comm_filtered_pa[, j]
    intersection <- sum(sp_i == 1 & sp_j == 1)
    union <- sum(sp_i == 1 | sp_j == 1)
    if (union > 0) {
      jaccard_matrix[i, j] <- intersection / union
      jaccard_matrix[j, i] <- jaccard_matrix[i, j]
    }
  }
}

# Build network from Jaccard > 0.3 (same threshold as original)
threshold <- 0.3
edge_indices_pa <- which(upper.tri(jaccard_matrix) & jaccard_matrix > threshold, arr.ind = TRUE)

cat("  Candidate edges (Jaccard > ", threshold, "):", nrow(edge_indices_pa), "\n")

if (nrow(edge_indices_pa) > 0) {
  edge_list_pa <- data.frame(
    sp1 = colnames(jaccard_matrix)[edge_indices_pa[, 1]],
    sp2 = colnames(jaccard_matrix)[edge_indices_pa[, 2]],
    weight = jaccard_matrix[edge_indices_pa],
    stringsAsFactors = FALSE
  )

  g_pa <- graph_from_data_frame(edge_list_pa[, c("sp1", "sp2")], directed = FALSE)
  E(g_pa)$weight <- edge_list_pa$weight

  transitivity_pa <- transitivity(g_pa, type = "global")

  cat("\n  Presence/Absence Network (raw Jaccard):\n")
  cat("    - Nodes:", vcount(g_pa), "\n")
  cat("    - Edges:", ecount(g_pa), "\n")
  cat("    - Transitivity:", round(transitivity_pa, 4), "\n\n")
} else {
  transitivity_pa <- NA
  g_pa <- NULL
  cat("  No edges above threshold - cannot build network\n\n")
}

# ============================================================================
# PART 4: PARTIAL CORRELATION NETWORK (REGRESS OUT ABUNDANCE)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 4: PARTIAL CORRELATION NETWORK (Abundance-Corrected)\n")
cat("------------------------------------------------------------\n\n")

cat("4.1 Computing species occurrence rates...\n")

# For each species pair, compute partial correlation controlling for
# the expected co-occurrence based on individual occurrence rates.
# Under null: P(both present) = P(A) * P(B) [independence]

species_occurrence_rate <- colMeans(comm_filtered_pa)

cat("  Mean occurrence rate:", round(mean(species_occurrence_rate), 3), "\n")
cat("  Range:", round(min(species_occurrence_rate), 3), "-",
    round(max(species_occurrence_rate), 3), "\n\n")

cat("4.2 Computing abundance-corrected associations...\n")
cat("  Method: Standardized residuals from expected co-occurrence\n\n")

# Compute observed vs expected co-occurrence
# Expected under independence: n_samples * P(A) * P(B)
# Standardized residual: (Observed - Expected) / sqrt(Expected * (1 - p_both))

n_samples <- nrow(comm_filtered_pa)
corrected_assoc <- matrix(0, nrow = n_sp, ncol = n_sp)
colnames(corrected_assoc) <- rownames(corrected_assoc) <- colnames(comm_filtered_pa)

for (i in 1:(n_sp - 1)) {
  for (j in (i + 1):n_sp) {
    p_i <- species_occurrence_rate[i]
    p_j <- species_occurrence_rate[j]

    # Expected co-occurrence under independence
    expected <- n_samples * p_i * p_j

    # Observed co-occurrence
    observed <- sum(comm_filtered_pa[, i] == 1 & comm_filtered_pa[, j] == 1)

    # Standardized residual (chi-square type)
    # Variance under independence: n * p_i * p_j * (1 - p_i * p_j)
    variance <- n_samples * p_i * p_j * (1 - p_i * p_j)

    if (variance > 0) {
      std_resid <- (observed - expected) / sqrt(variance)
      corrected_assoc[i, j] <- std_resid
      corrected_assoc[j, i] <- std_resid
    }
  }
}

# Convert standardized residuals to correlation-like scale [-1, 1]
# Use phi coefficient: sqrt(chi2 / n)
phi_matrix <- corrected_assoc / sqrt(n_samples)

# Build network from significant positive associations
# Use z > 1.96 as threshold (p < 0.05 one-tailed)
z_threshold <- 1.96
edge_indices_corrected <- which(upper.tri(corrected_assoc) & corrected_assoc > z_threshold, arr.ind = TRUE)

cat("  Significant positive associations (z > 1.96):", nrow(edge_indices_corrected), "\n")

if (nrow(edge_indices_corrected) > 0) {
  edge_list_corrected <- data.frame(
    sp1 = colnames(corrected_assoc)[edge_indices_corrected[, 1]],
    sp2 = colnames(corrected_assoc)[edge_indices_corrected[, 2]],
    z_score = corrected_assoc[edge_indices_corrected],
    phi = phi_matrix[edge_indices_corrected],
    stringsAsFactors = FALSE
  )

  g_corrected <- graph_from_data_frame(edge_list_corrected[, c("sp1", "sp2")], directed = FALSE)
  E(g_corrected)$weight <- edge_list_corrected$phi
  E(g_corrected)$z_score <- edge_list_corrected$z_score

  transitivity_corrected <- transitivity(g_corrected, type = "global")

  cat("\n  Abundance-Corrected Network:\n")
  cat("    - Nodes:", vcount(g_corrected), "\n")
  cat("    - Edges:", ecount(g_corrected), "\n")
  cat("    - Network density:", round(edge_density(g_corrected), 4), "\n")
  cat("    - Transitivity:", round(transitivity_corrected, 4), "\n\n")
} else {
  transitivity_corrected <- NA
  g_corrected <- NULL
  cat("  No significant positive associations found!\n\n")
}

# ============================================================================
# PART 5: NULL MODEL COMPARISON FOR CORRECTED NETWORK
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 5: NULL MODEL COMPARISON\n")
cat("------------------------------------------------------------\n\n")

# Generate null models for the abundance-corrected network

run_null_model <- function(g, n_permutations = 1000) {
  if (is.null(g) || vcount(g) < 3 || ecount(g) < 3) {
    return(list(
      obs_transitivity = NA,
      null_mean = NA,
      null_sd = NA,
      z_score = NA,
      p_value = NA,
      ratio_to_null = NA
    ))
  }

  obs_transitivity <- transitivity(g, type = "global")
  obs_degree_seq <- degree(g)

  null_trans <- numeric(n_permutations)

  set.seed(123)
  for (i in 1:n_permutations) {
    # Configuration model preserving degree sequence
    g_null <- tryCatch(
      sample_degseq(obs_degree_seq, method = "simple"),
      error = function(e) {
        tryCatch(
          sample_degseq(obs_degree_seq, method = "simple.no.multiple"),
          error = function(e2) erdos.renyi.game(vcount(g), edge_density(g), type = "gnp")
        )
      }
    )
    null_trans[i] <- transitivity(g_null, type = "global")
  }

  null_mean <- mean(null_trans, na.rm = TRUE)
  null_sd <- sd(null_trans, na.rm = TRUE)
  z_score <- (obs_transitivity - null_mean) / null_sd
  p_value <- 2 * (1 - pnorm(abs(z_score)))

  list(
    obs_transitivity = obs_transitivity,
    null_mean = null_mean,
    null_sd = null_sd,
    z_score = z_score,
    p_value = p_value,
    ratio_to_null = obs_transitivity / null_mean,
    null_distribution = null_trans
  )
}

cat("5.1 Running null models (1000 permutations each)...\n\n")

# Original network null comparison
cat("  Original network (volume-corrected only)...\n")
null_original <- run_null_model(g_original, 1000)

# Presence/absence network
if (!is.null(g_pa)) {
  cat("  Presence/absence network...\n")
  null_pa <- run_null_model(g_pa, 1000)
} else {
  null_pa <- list(obs_transitivity = NA, null_mean = NA, null_sd = NA,
                  z_score = NA, p_value = NA, ratio_to_null = NA)
}

# Abundance-corrected network
if (!is.null(g_corrected)) {
  cat("  Abundance-corrected network...\n")
  null_corrected <- run_null_model(g_corrected, 1000)
} else {
  null_corrected <- list(obs_transitivity = NA, null_mean = NA, null_sd = NA,
                         z_score = NA, p_value = NA, ratio_to_null = NA)
}

cat("\n5.2 Null Model Comparison Results:\n\n")

results_comparison <- data.frame(
  network = c("Original (volume-corrected)", "Presence/Absence (raw Jaccard)",
              "Abundance-Corrected (standardized residuals)"),
  n_nodes = c(vcount(g_original),
              ifelse(is.null(g_pa), NA, vcount(g_pa)),
              ifelse(is.null(g_corrected), NA, vcount(g_corrected))),
  n_edges = c(ecount(g_original),
              ifelse(is.null(g_pa), NA, ecount(g_pa)),
              ifelse(is.null(g_corrected), NA, ecount(g_corrected))),
  density = c(edge_density(g_original),
              ifelse(is.null(g_pa), NA, edge_density(g_pa)),
              ifelse(is.null(g_corrected), NA, edge_density(g_corrected))),
  obs_transitivity = c(null_original$obs_transitivity, null_pa$obs_transitivity,
                       null_corrected$obs_transitivity),
  null_mean = c(null_original$null_mean, null_pa$null_mean, null_corrected$null_mean),
  null_sd = c(null_original$null_sd, null_pa$null_sd, null_corrected$null_sd),
  z_score = c(null_original$z_score, null_pa$z_score, null_corrected$z_score),
  p_value = c(null_original$p_value, null_pa$p_value, null_corrected$p_value),
  ratio_to_null = c(null_original$ratio_to_null, null_pa$ratio_to_null,
                    null_corrected$ratio_to_null),
  stringsAsFactors = FALSE
)

print(results_comparison %>%
        mutate(across(where(is.numeric), ~round(., 4))))

# ============================================================================
# PART 6: IDENTIFY WHAT'S DRIVING THE PATTERN
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 6: IDENTIFYING DRIVERS OF CLUSTERING\n")
cat("------------------------------------------------------------\n\n")

if (!is.null(g_corrected) && ecount(g_corrected) > 0) {

  cat("6.1 Strongest species associations (top 20 by z-score)...\n\n")

  top_edges <- edge_list_corrected %>%
    arrange(desc(z_score)) %>%
    head(20)

  print(top_edges)

  cat("\n6.2 Taxonomic patterns in strong associations...\n\n")

  # Add taxonomic information
  get_taxon <- function(sp) {
    type_val <- cafi_clean %>%
      filter(otu == sp | species == sp) %>%
      pull(type) %>%
      unique()
    if (length(type_val) == 0) "unknown" else type_val[1]
  }

  top_edges$type_1 <- sapply(top_edges$sp1, get_taxon)
  top_edges$type_2 <- sapply(top_edges$sp2, get_taxon)
  top_edges$same_taxon <- top_edges$type_1 == top_edges$type_2

  cat("  Proportion of top 20 edges within same taxon:",
      round(mean(top_edges$same_taxon), 2), "\n")

  taxon_pairs <- table(paste(top_edges$type_1, top_edges$type_2, sep = "-"))
  cat("  Taxon pair frequencies:\n")
  print(sort(taxon_pairs, decreasing = TRUE))

  cat("\n6.3 Hub species in abundance-corrected network...\n\n")

  corrected_centrality <- data.frame(
    species = V(g_corrected)$name,
    degree = degree(g_corrected),
    eigenvector = eigen_centrality(g_corrected)$vector,
    local_trans = transitivity(g_corrected, type = "local"),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(degree))

  corrected_centrality$type <- sapply(corrected_centrality$species, get_taxon)

  cat("  Top 10 hub species by degree:\n")
  print(head(corrected_centrality, 10))

} else {
  cat("  Cannot analyze drivers - no abundance-corrected network available\n")
}

# ============================================================================
# PART 7: ALTERNATIVE TEST - CURVEBALL NULL MODEL
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 7: CURVEBALL NULL MODEL (Row-Column Marginal Preserved)\n")
cat("------------------------------------------------------------\n\n")

cat("7.1 Using curveball algorithm to preserve both:\n")
cat("    - Species occurrence totals (column sums)\n")
cat("    - Coral richness (row sums)\n")
cat("    This is a more conservative null that preserves abundance structure.\n\n")

# Implement curveball swap
curveball_swap <- function(m, n_swaps = 1000) {
  # m is a binary matrix
  m_new <- m
  n_samples <- nrow(m)
  n_species <- ncol(m)

  for (iter in 1:n_swaps) {
    # Pick two random rows
    rows <- sample(1:n_samples, 2)
    r1 <- rows[1]
    r2 <- rows[2]

    # Find species unique to each row
    sp1_only <- which(m_new[r1, ] == 1 & m_new[r2, ] == 0)
    sp2_only <- which(m_new[r1, ] == 0 & m_new[r2, ] == 1)

    if (length(sp1_only) > 0 && length(sp2_only) > 0) {
      # Randomly select one species from each set and swap
      swap1 <- sample(sp1_only, 1)
      swap2 <- sample(sp2_only, 1)

      m_new[r1, swap1] <- 0
      m_new[r1, swap2] <- 1
      m_new[r2, swap1] <- 1
      m_new[r2, swap2] <- 0
    }
  }
  return(m_new)
}

# Run curveball null model
cat("7.2 Running curveball null model (100 permutations, 10000 swaps each)...\n")

n_perm_curveball <- 100
curveball_transitivity <- numeric(n_perm_curveball)

# Use the same correlation method as original script
# (Spearman on volume-residualized presence)
# First, filter to corals that are in coral_master (to ensure volume data available)
matched_corals_pa <- rownames(comm_filtered_pa) %in% coral_master$coral_id
comm_filtered_pa_matched <- comm_filtered_pa[matched_corals_pa, , drop = FALSE]

# Get volume for these corals
volume_vec_cb <- coral_master$volume_field[match(rownames(comm_filtered_pa_matched), coral_master$coral_id)]
log_volume_cb <- log10(volume_vec_cb)

# Remove any NAs
valid_rows <- !is.na(log_volume_cb)
comm_filtered_pa_matched <- comm_filtered_pa_matched[valid_rows, , drop = FALSE]
log_volume_cb <- log_volume_cb[valid_rows]

cat("  Matrix for curveball:", nrow(comm_filtered_pa_matched), "corals x",
    ncol(comm_filtered_pa_matched), "species\n")

set.seed(456)
for (iter in 1:n_perm_curveball) {
  # Shuffle matrix preserving marginals
  comm_shuffled <- curveball_swap(comm_filtered_pa_matched, n_swaps = 10000)

  # Residualize on volume
  residual_shuffled <- matrix(NA, nrow = nrow(comm_shuffled), ncol = ncol(comm_shuffled))
  for (sp in 1:ncol(comm_shuffled)) {
    y <- comm_shuffled[, sp]
    fit <- tryCatch(
      suppressWarnings(glm(y ~ log_volume_cb, family = binomial)),
      error = function(e) NULL
    )
    if (!is.null(fit) && length(residuals(fit)) == nrow(comm_shuffled)) {
      residual_shuffled[, sp] <- residuals(fit, type = "deviance")
    } else {
      residual_shuffled[, sp] <- y - mean(y)
    }
  }

  # Compute correlation matrix
  cor_shuffled <- cor(residual_shuffled, method = "spearman", use = "pairwise.complete.obs")
  cor_shuffled[is.na(cor_shuffled)] <- 0

  # Build network with same threshold
  edge_idx <- which(upper.tri(cor_shuffled) & cor_shuffled > 0.3, arr.ind = TRUE)

  if (is.matrix(edge_idx) && nrow(edge_idx) > 0) {
    el <- data.frame(
      sp1 = colnames(cor_shuffled)[edge_idx[, 1]],
      sp2 = colnames(cor_shuffled)[edge_idx[, 2]],
      stringsAsFactors = FALSE
    )
    if (nrow(el) >= 1 && ncol(el) >= 2) {
      g_null <- tryCatch(
        graph_from_data_frame(el, directed = FALSE),
        error = function(e) NULL
      )
      if (!is.null(g_null) && ecount(g_null) > 0) {
        curveball_transitivity[iter] <- transitivity(g_null, type = "global")
      } else {
        curveball_transitivity[iter] <- NA
      }
    } else {
      curveball_transitivity[iter] <- NA
    }
  } else {
    curveball_transitivity[iter] <- NA
  }

  if (iter %% 20 == 0) cat("    ", iter, "permutations completed\n")
}

curveball_mean <- mean(curveball_transitivity, na.rm = TRUE)
curveball_sd <- sd(curveball_transitivity, na.rm = TRUE)
curveball_z <- (transitivity_original - curveball_mean) / curveball_sd
curveball_p <- 2 * (1 - pnorm(abs(curveball_z)))
curveball_ratio <- transitivity_original / curveball_mean

cat("\n7.3 Curveball Null Model Results:\n")
cat("    Observed transitivity:", round(transitivity_original, 4), "\n")
cat("    Curveball null mean:", round(curveball_mean, 4), "\n")
cat("    Curveball null SD:", round(curveball_sd, 4), "\n")
cat("    z-score:", round(curveball_z, 2), "\n")
cat("    p-value:", format(curveball_p, scientific = TRUE, digits = 3), "\n")
cat("    Ratio to null:", round(curveball_ratio, 2), "x\n\n")

# Add to results
results_comparison <- rbind(
  results_comparison,
  data.frame(
    network = "Curveball null (marginal-preserved)",
    n_nodes = n_nodes,
    n_edges = n_edges,
    density = edge_density(g_original),
    obs_transitivity = transitivity_original,
    null_mean = curveball_mean,
    null_sd = curveball_sd,
    z_score = curveball_z,
    p_value = curveball_p,
    ratio_to_null = curveball_ratio,
    stringsAsFactors = FALSE
  )
)

# ============================================================================
# PART 8: SUMMARY AND CONCLUSION
# ============================================================================

cat("============================================================\n")
cat("    SUMMARY: IS TRANSITIVITY AN ABUNDANCE ARTIFACT?\n")
cat("============================================================\n\n")

# Determine conclusion based on results
if (!is.na(null_corrected$z_score) && null_corrected$z_score > 2) {
  conclusion <- "ROBUST"
  explanation <- paste0(
    "The transitivity signal SURVIVES abundance correction:\n",
    "  - Abundance-corrected network: z = ", round(null_corrected$z_score, 1),
    " (", round(null_corrected$ratio_to_null, 2), "x null)\n",
    "  - Curveball null: z = ", round(curveball_z, 1),
    " (", round(curveball_ratio, 2), "x null)\n",
    "\nThis suggests the clustering reflects TRUE ECOLOGICAL ASSOCIATIONS,\n",
    "not an abundance artifact. Species that frequently co-occur do so more\n",
    "than expected even after controlling for their individual occurrence rates."
  )
} else if (curveball_z > 2) {
  conclusion <- "PARTIALLY ROBUST"
  explanation <- paste0(
    "The transitivity signal is REDUCED but still significant after correction:\n",
    "  - Original z = 36.1, Curveball z = ", round(curveball_z, 1), "\n",
    "\nThe pattern is partially driven by abundance but retains biological signal."
  )
} else {
  conclusion <- "ARTIFACT"
  explanation <- paste0(
    "The transitivity signal is an ABUNDANCE ARTIFACT:\n",
    "  - Original z = 36.1, but after abundance correction:\n",
    "  - Curveball z = ", round(curveball_z, 1), " (not significant)\n",
    "\nThis is analogous to the composition-divergence finding (Q2),\n",
    "where the pattern disappeared after rarefaction."
  )
}

cat("CONCLUSION:", conclusion, "\n\n")
cat(explanation, "\n\n")

# Abundance-centrality correlation summary
cat("Abundance-Centrality Correlations:\n")
cat("  - Abundance vs Degree: rho =", round(cor_abund_degree$estimate, 3),
    ", p =", format(cor_abund_degree$p.value, scientific = TRUE, digits = 2), "\n")
cat("  - Occurrence vs Degree: rho =", round(cor_occur_degree$estimate, 3),
    ", p =", format(cor_occur_degree$p.value, scientific = TRUE, digits = 2), "\n")

# ============================================================================
# PART 9: SAVE OUTPUTS
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 9: SAVING OUTPUTS\n")
cat("------------------------------------------------------------\n\n")

# 9.1 Save results table
save_table(results_comparison, "network_abundance_artifact")
cat("  Saved: output/tables/network_abundance_artifact.csv\n")

# 9.2 Save detailed analysis report
report_text <- paste0(
"# CAFI Network Transitivity: Abundance Artifact Analysis
# =======================================================
# Date: ", Sys.Date(), "

## Background
The CAFI co-occurrence network (script 06_network_analysis.R) found highly
non-random clustering:
  - Transitivity = ", round(transitivity_original, 4), " (observed)
  - Null model = ", round(null_original$null_mean, 4), " +/- ", round(null_original$null_sd, 4), "
  - z-score = ", round(null_original$z_score, 1), "
  - Ratio to null = ", round(null_original$ratio_to_null, 2), "x

This analysis tests whether this pattern is an abundance artifact.

## Methods
1. Original network controls for CORAL VOLUME (logistic GLM residuals)
   But does NOT control for species-level abundance confounds.

2. Abundance-centrality correlation tests whether common species have
   artificially high connectivity.

3. Abundance-corrected network uses standardized residuals from expected
   co-occurrence under independence: (Observed - Expected) / sqrt(Variance)

4. Curveball null model preserves both row and column marginals
   (coral richness AND species occurrence totals).

## Results

### Abundance-Centrality Correlations
  - Abundance vs Degree: rho = ", round(cor_abund_degree$estimate, 3),
    ", p = ", format(cor_abund_degree$p.value, scientific = TRUE, digits = 2), "
  - Occurrence vs Degree: rho = ", round(cor_occur_degree$estimate, 3),
    ", p = ", format(cor_occur_degree$p.value, scientific = TRUE, digits = 2), "

### Null Model Comparison
", paste(capture.output(print(results_comparison %>%
                               mutate(across(where(is.numeric), ~round(., 3))))),
         collapse = "\n"), "

## Conclusion
", conclusion, "

", explanation, "

## Implications for Q2 (Community Composition)
The pattern in Q2 (composition divergence) was shown to be an abundance artifact:
  - Before rarefaction: significant size-composition relationship
  - After rarefaction: relationship disappeared (p = 0.61)

This network analysis tests whether the transitivity finding has the same issue.

## Key Findings
1. Abundance-degree correlation: ",
  ifelse(cor_abund_degree$p.value < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"),
  " (rho = ", round(cor_abund_degree$estimate, 3), ")

2. Transitivity after abundance correction: ",
  ifelse(!is.na(null_corrected$z_score) && null_corrected$z_score > 2,
         "STILL SIGNIFICANT", "NOT SIGNIFICANT"),

"\n3. Curveball null comparison: z = ", round(curveball_z, 1), " (",
  ifelse(curveball_z > 2, "significant", "not significant"), ")

## Species Pairs Driving Clustering
",
if (!is.null(g_corrected) && exists("top_edges")) {
  paste0("Top 5 associations by z-score:\n",
         paste(apply(head(top_edges[, c("sp1", "sp2", "z_score")], 5), 1,
                     function(x) paste0("  ", x[1], " - ", x[2], " (z = ",
                                       round(as.numeric(x[3]), 2), ")")),
               collapse = "\n"))
} else {
  "  (No abundance-corrected network available)"
}
)

writeLines(report_text, file.path(PATHS$output, "network_transitivity_artifact_analysis.txt"))
cat("  Saved: output/network_transitivity_artifact_analysis.txt\n")

# 9.3 Save species metrics with abundance and centrality
if (exists("species_metrics")) {
  save_table(species_metrics, "species_abundance_centrality")
  cat("  Saved: output/tables/species_abundance_centrality.csv\n")
}

cat("\n============================================================\n")
cat("    NETWORK ABUNDANCE ARTIFACT TEST COMPLETE\n")
cat("============================================================\n\n")
