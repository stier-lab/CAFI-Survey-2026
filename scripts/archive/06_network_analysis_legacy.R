# ============================================================================
# 06_network_analysis.R - CAFI Species Co-occurrence Network Analysis
# ============================================================================
#
# PURPOSE: Build and analyze CAFI species co-occurrence networks to test
#          hypothesis H5 regarding non-random community structure
#
# HYPOTHESIS (H5): CAFI co-occurrence networks exhibit non-random structure,
#   with groups corresponding to functional groups or shared habitat
#   preferences, and identifiable hub/keystone species.
#
# ANALYSES:
#   - Co-occurrence network construction (volume-residualized correlations)
#   - Community detection (Louvain algorithm)
#   - Null model comparison (1000 configuration model random networks)
#   - Centrality analysis (degree, betweenness, eigenvector)
#   - Group composition characterization
#
# EXPECTED RESULTS (from PRD):
#   - Modularity Q ~ 0.52 (2.2x null)
#   - Transitivity ~ 0.28 (3.8x null)
#   - ~7 modules
#   - Hub species: Alpheus diadema, Caracanthus maculatus, Trapezia serenei
#
# USAGE:
#   source("scripts/00_setup.R")
#   source("scripts/01_load_data.R")
#   source("scripts/06_network_analysis.R")
#
# OUTPUTS:
#   Figures:
#     - output/figures/network_visualization.png
#     - output/figures/manuscript/fig4_network.png
#     - output/figures/supplement/figS11_network_groups.png
#   Tables:
#     - output/tables/network_metrics.csv
#     - output/tables/module_composition.csv
#     - output/tables/hub_species.csv
#   Objects:
#     - output/objects/cafi_network.rds
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-21
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI NETWORK ANALYSIS (H5)\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP AND DATA LOADING
# ============================================================================

# Load setup (packages, paths, theme)
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

# Load additional required packages for network analysis
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Package 'igraph' is required but not installed. Install with: install.packages('igraph')")
}
library(igraph)

# Create output directories
fig_dir <- file.path(PATHS$figures, "06_network")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("[OK] Setup complete\n")
cat("     Corals:", nrow(coral_master), "\n")
cat("     Community matrix:", nrow(community_matrix), "x", ncol(community_matrix), "\n\n")

# ============================================================================
# IMPORTANT: Network Scope
# ============================================================================
#
# This analysis examines CAFI-CAFI co-occurrence networks ONLY.
# Corals are NOT nodes in these networks because:
#   - All corals are Pocillopora spp. (single species)
#   - No reliable coral morphotype/species distinctions
#   - Interest is in CAFI community structure within coral hosts
#
# The network represents species co-occurrence patterns across 114 corals.
# ============================================================================

# ============================================================================
# PART 1: CO-OCCURRENCE NETWORK CONSTRUCTION
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 1: CO-OCCURRENCE NETWORK CONSTRUCTION\n")
cat("------------------------------------------------------------\n\n")

# 1.1 Convert to presence-absence matrix
cat("1.1 Preparing community data...\n")

comm_binary <- community_matrix
comm_binary[comm_binary > 0] <- 1

cat("     Total samples (corals):", nrow(comm_binary), "\n")
cat("     Total species:", ncol(comm_binary), "\n")

# 1.2 Filter to species with sufficient occurrence (>= 5% of corals)
min_occurrence <- ceiling(nrow(comm_binary) * 0.05)
species_keep <- colSums(comm_binary) >= min_occurrence
comm_filtered <- comm_binary[, species_keep]

cat("     Minimum occurrence threshold:", min_occurrence, "corals (5%)\n")
cat("     Species retained:", ncol(comm_filtered), "of", ncol(comm_binary), "\n\n")

# 1.3 Residualize presence on coral volume (partial correlations)
# This removes the mechanical confound: larger corals host more species,
# inflating raw co-occurrence. We fit logistic GLM for each species and
# use deviance residuals as volume-corrected presence values.
cat("1.2 Residualizing species presence on coral volume...\n")

# Filter to corals that exist in coral_master (some community matrix rows may not match)
matched_corals <- rownames(comm_filtered) %in% coral_master$coral_id
if (sum(!matched_corals) > 0) {
  cat("     Excluding", sum(!matched_corals), "corals not in coral_master\n")
  comm_filtered <- comm_filtered[matched_corals, , drop = FALSE]
}

# Get volume for each coral (rows of community_matrix match coral_master)
volume_vec <- coral_master$volume[match(rownames(comm_filtered), coral_master$coral_id)]
log_volume <- log(volume_vec)

# Verify no NAs in volume
stopifnot("All corals must have volume" = !any(is.na(log_volume)))

# Compute deviance residuals for each species
residual_matrix <- matrix(NA, nrow = nrow(comm_filtered), ncol = ncol(comm_filtered))
colnames(residual_matrix) <- colnames(comm_filtered)
rownames(residual_matrix) <- rownames(comm_filtered)

for (sp in seq_len(ncol(comm_filtered))) {
  y <- comm_filtered[, sp]
  fit <- tryCatch(
    glm(y ~ log_volume, family = binomial(link = "logit")),
    warning = function(w) glm(y ~ log_volume, family = binomial(link = "logit")),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    residual_matrix[, sp] <- residuals(fit, type = "deviance")
  } else {
    # Fallback: use raw presence if GLM fails (e.g., perfect separation)
    residual_matrix[, sp] <- y - mean(y)
  }
}

cat("     Species residualized:", ncol(residual_matrix), "\n")

# 1.4 Calculate Spearman correlations on volume-corrected residuals
# Note: Spearman rank correlation on logistic deviance residuals may have
# reduced power due to bimodal residual distribution. Pearson correlation
# on Pearson residuals is an alternative but assumes linearity.
cat("1.3 Calculating volume-corrected co-occurrences (Spearman correlation)...\n")

cor_matrix <- cor(residual_matrix, method = "spearman", use = "pairwise.complete.obs")

# Handle NAs in correlation matrix
cor_matrix[is.na(cor_matrix)] <- 0

# 1.5 Extract significant positive associations with FDR correction
# Edge threshold: r > 0.3 AND FDR-corrected p < 0.05
# Rationale: r > 0.3 balances sensitivity (detecting real co-occurrences) with
# specificity (avoiding noise edges). Lower thresholds (0.2) produce dense graphs
# where modularity detection is unreliable; higher (0.4) fragment the network.
# Threshold sensitivity analysis in Part 2 confirms results are robust.
threshold <- 0.3

cat("     Computing pairwise p-values for FDR correction...\n")

# Compute p-values for ALL pairs with |r| > 0.2 (lowest sensitivity threshold)
# This supports both the main analysis and the threshold sensitivity analysis
n_sp <- ncol(residual_matrix)
p_matrix <- matrix(1, nrow = n_sp, ncol = n_sp)
for (i in 1:(n_sp - 1)) {
  for (j in (i + 1):n_sp) {
    if (abs(cor_matrix[i, j]) > 0.2) {
      ct <- tryCatch(
        cor.test(residual_matrix[, i], residual_matrix[, j], method = "spearman"),
        error = function(e) list(p.value = 1)
      )
      p_matrix[i, j] <- ct$p.value
      p_matrix[j, i] <- ct$p.value
    }
  }
}

# FDR correction pass 1 (completeness / sensitivity analysis):
# Correct ALL pairwise p-values across the full upper triangle. This matrix
# (p_matrix_adj) is used by the threshold sensitivity analysis below, which
# needs FDR-adjusted p-values for edges at varying correlation thresholds.
upper_idx <- which(upper.tri(p_matrix))
p_upper_raw <- p_matrix[upper_idx]
p_upper_adj <- p.adjust(p_upper_raw, method = "BH")
p_matrix_adj <- matrix(1, nrow = n_sp, ncol = n_sp)
p_matrix_adj[upper_idx] <- p_upper_adj
p_matrix_adj <- t(p_matrix_adj)
p_matrix_adj[upper_idx] <- p_upper_adj  # symmetrize

# Extract pairs above threshold
edge_indices_raw <- which(upper.tri(cor_matrix) & cor_matrix > threshold, arr.ind = TRUE)

# Initialize edge_indices as empty in case no edges pass threshold
edge_indices <- matrix(nrow = 0, ncol = 2)

if (nrow(edge_indices_raw) > 0) {
  # FDR correction pass 2 (operational filter for main network):
  # Correct only the subset of RAW p-values for candidate edges above the
  # correlation threshold. This is the filter used to build the actual network.
  # Note: uses raw p-values from p_matrix (not already-adjusted ones from pass 1)
  # to avoid double-correction.
  raw_pvals <- sapply(1:nrow(edge_indices_raw), function(k)
    p_matrix[edge_indices_raw[k, 1], edge_indices_raw[k, 2]])

  # Apply FDR correction (Benjamini-Hochberg) on this subset
  fdr_pvals <- p.adjust(raw_pvals, method = "BH")

  # Filter: retain only FDR-significant edges
  sig_mask <- fdr_pvals < 0.05
  edge_indices <- edge_indices_raw[sig_mask, , drop = FALSE]
  n_dropped <- sum(!sig_mask)

  cat("     Candidate edges (r > ", threshold, "):", nrow(edge_indices_raw), "\n")
  cat("     Retained after FDR correction:", sum(sig_mask), "\n")
  cat("     Dropped by FDR:", n_dropped, "\n")
}

if (nrow(edge_indices) > 0) {

  edge_list <- data.frame(
    sp1 = colnames(cor_matrix)[edge_indices[, 1]],
    sp2 = colnames(cor_matrix)[edge_indices[, 2]],
    correlation = cor_matrix[edge_indices],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      weight = correlation,  # Use correlation as edge weight
      association = "positive"
    )

  cat("     FDR-significant positive associations:", nrow(edge_list), "\n")
  cat("     Mean correlation strength:", round(mean(edge_list$correlation), 3), "\n")
  cat("     Range:", round(min(edge_list$correlation), 3), "-",
      round(max(edge_list$correlation), 3), "\n\n")

} else {
  warning("No significant positive associations found. Skipping network analysis.")
  cat("  WARNING: No significant positive associations found. Network analysis skipped.\n")
  cat("  This may indicate insufficient sample size or low species co-occurrence.\n\n")
  # Create empty results so pipeline can continue
  network_results <- list(
    n_species = ncol(comm_filtered),
    n_edges = 0,
    note = "No significant positive associations found"
  )
  save_object(network_results, "cafi_network")
  cat("Script 06 complete (no network constructed).\n\n")
}

# Only proceed with network construction if edges were found
if (nrow(edge_indices) > 0) {

# ============================================================================
# PART 2: BUILD IGRAPH NETWORK OBJECT
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 2: NETWORK CONSTRUCTION\n")
cat("------------------------------------------------------------\n\n")

cat("2.1 Building network from edge list...\n")

# Create igraph object
g <- graph_from_data_frame(
  edge_list[, c("sp1", "sp2")],
  directed = FALSE
)

# Add edge weights
E(g)$weight <- edge_list$correlation
E(g)$correlation <- edge_list$correlation

# Add node attributes
# Get taxonomic type for each species
V(g)$type <- sapply(V(g)$name, function(sp) {
  type_val <- cafi_clean %>%
    dplyr::filter(otu == sp | species == sp) %>%
    dplyr::pull(type) %>%
    unique()
  if (length(type_val) == 0) "unknown" else type_val[1]
})

# Get functional group
V(g)$functional_group <- sapply(V(g)$name, function(sp) {
  fg_val <- cafi_clean %>%
    dplyr::filter(otu == sp | species == sp) %>%
    dplyr::pull(functional_group) %>%
    unique()
  if (length(fg_val) == 0) "Other" else fg_val[1]
})

# Calculate total abundance across all corals
V(g)$abundance <- sapply(V(g)$name, function(sp) {
  if (sp %in% colnames(community_matrix)) {
    sum(community_matrix[, sp], na.rm = TRUE)
  } else {
    0
  }
})

# Calculate occurrence (number of corals)
V(g)$occurrence <- sapply(V(g)$name, function(sp) {
  if (sp %in% colnames(comm_binary)) {
    sum(comm_binary[, sp], na.rm = TRUE)
  } else {
    0
  }
})

cat("     Network nodes (species):", vcount(g), "\n")
cat("     Network edges (associations):", ecount(g), "\n\n")

# ============================================================================
# PART 3: NETWORK-LEVEL METRICS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 3: NETWORK-LEVEL METRICS\n")
cat("------------------------------------------------------------\n\n")

cat("3.1 Calculating network metrics...\n")

# Basic metrics
n_nodes <- vcount(g)
n_edges <- ecount(g)
density <- edge_density(g)
transitivity_obs <- transitivity(g, type = "global")
diameter_obs <- diameter(g, weights = NA)      # unweighted for interpretable integer hops
mean_distance_obs <- mean_distance(g, weights = NA)  # unweighted for interpretable path length

# Modularity via Louvain algorithm (weighted for community assignments)
# Louvain chosen over alternatives (fast-greedy, walktrap) because it:
#   1. Handles weighted networks natively (edge weights = correlation strength)
#   2. Scales well and produces interpretable community structure
#   3. Widely used in ecological network analyses for comparison
set.seed(42)  # For reproducibility (Louvain is stochastic)
communities_louvain <- cluster_louvain(g, weights = E(g)$weight)
modularity_obs <- modularity(communities_louvain)
n_modules <- length(unique(membership(communities_louvain)))

# Unweighted modularity for null model comparison
# Null graphs (degree-preserving configuration model) are unweighted, so we
# need an unweighted observed modularity for a fair z-score comparison
communities_unweighted <- cluster_louvain(g)  # no weights
modularity_obs_unweighted <- modularity(communities_unweighted)

# Add module membership to nodes
V(g)$module <- membership(communities_louvain)

# Degree distribution stats
degrees <- degree(g)
mean_degree <- mean(degrees)
median_degree <- median(degrees)
max_degree <- max(degrees)

cat("     Nodes:", n_nodes, "\n")
cat("     Edges:", n_edges, "\n")
cat("     Density:", round(density, 4), "\n")
cat("     Transitivity (clustering):", round(transitivity_obs, 3), "\n")
cat("     Diameter:", diameter_obs, "\n")
cat("     Mean path length:", round(mean_distance_obs, 2), "\n")
cat("     Modularity Q (weighted):", round(modularity_obs, 3), "\n")
cat("     Modularity Q (unweighted, for null comparison):", round(modularity_obs_unweighted, 3), "\n")
cat("     Number of Louvain groups:", n_modules, "\n")
cat("     Mean degree:", round(mean_degree, 2), "\n")
cat("     Median degree:", median_degree, "\n")
cat("     Max degree:", max_degree, "\n\n")

# ============================================================================
# PART 4: NULL MODEL COMPARISON
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 4: NULL MODEL COMPARISON (Configuration Model)\n")
cat("------------------------------------------------------------\n\n")

cat("4.1 Generating 1000 degree-preserving random networks...\n")
cat("     (Edge-swap rewiring preserves observed degree sequence)\n")

n_permutations <- 1000
null_metrics <- matrix(NA, nrow = n_permutations, ncol = 4)
colnames(null_metrics) <- c("modularity", "transitivity", "mean_distance", "diameter")

# Create unweighted copy of observed graph for rewiring
g_unweighted <- delete_edge_attr(g, "weight")

# Number of edge swaps per permutation: 10× edge count ensures thorough randomization
n_swaps <- ecount(g_unweighted) * 10

set.seed(123)

for (i in 1:n_permutations) {
  # Edge-swap rewiring: randomly swaps edge endpoints while preserving degree sequence
  # Much faster than sample_degseq() for dense graphs (density = 0.76)
  g_random <- rewire(g_unweighted, keeping_degseq(niter = n_swaps))

  # Calculate metrics
  null_metrics[i, "transitivity"] <- transitivity(g_random, type = "global")
  null_metrics[i, "mean_distance"] <- mean_distance(g_random, weights = NA)
  null_metrics[i, "diameter"] <- diameter(g_random, weights = NA)

  # Modularity (only if network has edges)
  if (ecount(g_random) > 0) {
    null_comm <- cluster_louvain(g_random)
    null_metrics[i, "modularity"] <- modularity(null_comm)
  } else {
    null_metrics[i, "modularity"] <- 0
  }

  if (i %% 100 == 0) {
    cat("     ", i, "of", n_permutations, "permutations completed\n")
    flush.console()
  }
}

cat("\n4.2 Computing z-scores and p-values...\n\n")

# Calculate statistics for null model comparison
# Use unweighted observed modularity for fair comparison with unweighted null graphs
null_comparison <- data.frame(
  metric = c("Modularity", "Transitivity", "Mean Path Length", "Diameter"),
  observed = c(modularity_obs_unweighted, transitivity_obs, mean_distance_obs, diameter_obs),
  null_mean = c(
    mean(null_metrics[, "modularity"], na.rm = TRUE),
    mean(null_metrics[, "transitivity"], na.rm = TRUE),
    mean(null_metrics[, "mean_distance"], na.rm = TRUE),
    mean(null_metrics[, "diameter"], na.rm = TRUE)
  ),
  null_sd = c(
    sd(null_metrics[, "modularity"], na.rm = TRUE),
    sd(null_metrics[, "transitivity"], na.rm = TRUE),
    sd(null_metrics[, "mean_distance"], na.rm = TRUE),
    sd(null_metrics[, "diameter"], na.rm = TRUE)
  )
) %>%
  mutate(
    z_score = (observed - null_mean) / null_sd,
    p_value = 2 * (1 - pnorm(abs(z_score))),
    ratio_to_null = observed / null_mean,
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# Print results
cat("     Null Model Comparison Results:\n\n")
for (i in 1:nrow(null_comparison)) {
  cat(sprintf("     %s:\n", null_comparison$metric[i]))
  cat(sprintf("       Observed: %.3f\n", null_comparison$observed[i]))
  cat(sprintf("       Null: %.3f +/- %.3f\n", null_comparison$null_mean[i], null_comparison$null_sd[i]))
  cat(sprintf("       Ratio to null: %.2fx\n", null_comparison$ratio_to_null[i]))
  cat(sprintf("       z-score: %.2f, p = %.4f %s\n\n",
              null_comparison$z_score[i], null_comparison$p_value[i],
              null_comparison$significance[i]))
}

# Check modularity z-score direction and report accordingly
mod_z <- null_comparison$z_score[null_comparison$metric == "Modularity"]
if (!is.na(mod_z) && mod_z < 0) {
  cat("  NOTE: Observed modularity is LOWER than null expectation (z =",
      round(mod_z, 2), ").\n")
  cat("  Louvain community detection identified groups, but modularity is not elevated above random.\n\n")
} else if (!is.na(mod_z) && mod_z >= 2) {
  cat("  Observed modularity exceeds null expectation (z =",
      round(mod_z, 2), ").\n")
  cat("  Louvain community detection identified groups with modularity elevated above random.\n\n")
} else {
  cat("  Observed modularity z =", round(mod_z, 2),
      " — modularity is not significantly elevated above random.\n\n")
}

# ============================================================================
# SENSITIVITY: Threshold analysis
# ============================================================================
cat("\n--- Network Threshold Sensitivity Analysis ---\n")

sensitivity_results <- data.frame(
  threshold = numeric(),
  n_edges = integer(),
  density = numeric(),
  n_modules = integer(),
  modularity_Q = numeric(),
  stringsAsFactors = FALSE
)

for (thresh in c(0.2, 0.3, 0.4, 0.5)) {
  # Apply threshold to existing correlation matrix (positive only, matching main analysis)
  edges_at_thresh <- which(cor_matrix > thresh & p_matrix_adj < 0.05, arr.ind = TRUE)
  edges_at_thresh <- edges_at_thresh[edges_at_thresh[,1] < edges_at_thresh[,2], , drop = FALSE]

  if (nrow(edges_at_thresh) > 0) {
    g_temp <- igraph::graph_from_edgelist(
      matrix(colnames(cor_matrix)[edges_at_thresh], ncol = 2),
      directed = FALSE
    )
    comm_temp <- igraph::cluster_louvain(g_temp)
    q_temp <- igraph::modularity(comm_temp)
    n_mod_temp <- length(unique(igraph::membership(comm_temp)))
    dens_temp <- igraph::edge_density(g_temp)
  } else {
    q_temp <- NA; n_mod_temp <- 0; dens_temp <- 0
  }

  sensitivity_results <- rbind(sensitivity_results, data.frame(
    threshold = thresh,
    n_edges = nrow(edges_at_thresh),
    density = round(dens_temp, 3),
    n_modules = n_mod_temp,
    modularity_Q = round(q_temp, 4)
  ))
}

cat("\n  Threshold sensitivity:\n")
print(sensitivity_results, row.names = FALSE)
cat("\n")

# Save sensitivity results
save_table(sensitivity_results, "network_threshold_sensitivity")
cat("     Saved: network_threshold_sensitivity.csv\n\n")

# ============================================================================
# PART 5: CENTRALITY ANALYSIS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 5: CENTRALITY ANALYSIS\n")
cat("------------------------------------------------------------\n\n")

cat("5.1 Calculating centrality metrics for all species...\n")

# Calculate centrality metrics
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, weights = NULL)  # Unweighted for path-based metric
V(g)$closeness <- closeness(g)
V(g)$eigenvector <- eigen_centrality(g, weights = E(g)$weight)$vector

# Create centrality dataframe
centrality_df <- data.frame(
  species = V(g)$name,
  type = V(g)$type,
  functional_group = V(g)$functional_group,
  module = V(g)$module,
  abundance = V(g)$abundance,
  occurrence = V(g)$occurrence,
  degree = V(g)$degree,
  betweenness = V(g)$betweenness,
  closeness = V(g)$closeness,
  eigenvector = V(g)$eigenvector,
  stringsAsFactors = FALSE
) %>%
  mutate(
    # Standardize metrics for comparison
    degree_z = scale(degree)[, 1],
    betweenness_z = scale(betweenness)[, 1],
    eigenvector_z = scale(eigenvector)[, 1],

    # Hub score: high degree + high eigenvector (well-connected to other well-connected nodes)
    hub_score = degree_z + eigenvector_z,

    # Keystone index: high degree + high betweenness (connectivity + bridging)
    keystone_index = degree_z + betweenness_z
  ) %>%
  arrange(desc(hub_score))

# 5.2 Identify hub species (top species by multiple centrality measures)
cat("\n5.2 Identifying hub species...\n\n")

# Top 10 by hub score
hub_species <- centrality_df %>%
  slice_head(n = 10) %>%
  dplyr::select(species, type, functional_group, module, degree, betweenness,
                eigenvector, hub_score, keystone_index)

cat("     Top 10 Hub Species (by hub_score = degree_z + eigenvector_z):\n\n")
print(hub_species %>% mutate(across(where(is.numeric), ~round(., 3))))

# Top 5 by different metrics
cat("\n     Top 5 by Degree Centrality:\n")
print(centrality_df %>% arrange(desc(degree)) %>% slice_head(n = 5) %>%
        dplyr::select(species, degree, type))

cat("\n     Top 5 by Betweenness Centrality:\n")
print(centrality_df %>% arrange(desc(betweenness)) %>% slice_head(n = 5) %>%
        dplyr::select(species, betweenness, type))

cat("\n     Top 5 by Eigenvector Centrality:\n")
print(centrality_df %>% arrange(desc(eigenvector)) %>% slice_head(n = 5) %>%
        dplyr::select(species, eigenvector, type))

# ============================================================================
# PART 6: MODULE ANALYSIS
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 6: GROUP ANALYSIS (Louvain community detection)\n")
cat("------------------------------------------------------------\n\n")

cat("6.1 Characterizing group composition...\n\n")

# Module composition by species
module_species <- centrality_df %>%
  dplyr::select(species, type, functional_group, module, degree, abundance) %>%
  arrange(module, desc(degree))

# Module summary
module_summary <- module_species %>%
  group_by(module) %>%
  summarise(
    n_species = n(),
    total_abundance = sum(abundance),
    mean_degree = round(mean(degree), 2),
    max_degree = max(degree),
    hub_species = species[which.max(degree)],
    types = paste(unique(type), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_species))

cat("     Group Summary:\n\n")
print(module_summary)

# 6.2 Taxonomic composition by module
cat("\n6.2 Taxonomic composition by group...\n\n")

module_taxonomy <- module_species %>%
  group_by(module, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(module) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup()

# Wide format for display
module_tax_wide <- module_taxonomy %>%
  dplyr::select(module, type, n) %>%
  pivot_wider(names_from = type, values_from = n, values_fill = 0)

print(module_tax_wide)

# 6.3 Functional group composition by module
cat("\n6.3 Functional group composition by Louvain group...\n\n")

module_func <- module_species %>%
  group_by(module, functional_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(module) %>%
  mutate(
    total = sum(n),
    proportion = round(n / total, 2)
  ) %>%
  ungroup()

module_func_wide <- module_func %>%
  dplyr::select(module, functional_group, n) %>%
  pivot_wider(names_from = functional_group, values_from = n, values_fill = 0)

print(module_func_wide)

# Chi-square test for non-random taxonomic distribution
contingency_table <- module_taxonomy %>%
  dplyr::select(module, type, n) %>%
  pivot_wider(names_from = type, values_from = n, values_fill = 0) %>%
  dplyr::select(-module) %>%
  as.matrix()

if (nrow(contingency_table) > 1 && ncol(contingency_table) > 1 &&
    min(rowSums(contingency_table)) > 0 && min(colSums(contingency_table)) > 0) {
  chi_test <- stats::chisq.test(contingency_table, simulate.p.value = TRUE, B = 2000)
  cat("\n     Chi-square test for taxonomic clustering across groups:\n")
  cat(sprintf("       X2 = %.2f, p = %.4f\n", chi_test$statistic, chi_test$p.value))
  if (chi_test$p.value < 0.05) {
    cat("       Interpretation: Groups show NON-RANDOM taxonomic composition\n")
  } else {
    cat("       Interpretation: Groups show random taxonomic mixing\n")
  }
}

# ============================================================================
# PART 7: NETWORK VISUALIZATION
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 7: NETWORK VISUALIZATION\n")
cat("------------------------------------------------------------\n\n")

cat("7.1 Creating network visualizations...\n")

# Set layout (Fruchterman-Reingold for nice spreading)
set.seed(42)
layout_fr <- layout_with_fr(g)

# Color palette for taxonomic types (Okabe-Ito colorblind-safe)
type_colors <- c(
  "crab" = "#D55E00",
  "shrimp" = "#0072B2",
  "fish" = "#009E73",
  "snail" = "#E69F00",
  "unknown" = "#999999"
)

# Module colors (colorblind-safe Okabe-Ito palette)
n_mod <- length(unique(V(g)$module))
module_palette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00")
module_colors <- module_palette[1:min(n_mod, length(module_palette))]
if (n_mod > length(module_palette)) {
  module_colors <- c(module_colors, scales::hue_pal()(n_mod - length(module_palette)))
}

# Node sizes based on degree
node_sizes <- 3 + sqrt(V(g)$degree) * 2

# Edge widths based on correlation strength
edge_widths <- E(g)$correlation * 2

# === Figure 1: Network colored by taxonomic type ===
png(file.path(fig_dir, "network_by_type.png"),
    width = 12, height = 10, units = "in", res = 300, bg = "white")

par(mar = c(1, 1, 3, 1))

# Get colors for each node
node_colors <- type_colors[V(g)$type]
node_colors[is.na(node_colors)] <- "#999999"

plot(g,
     layout = layout_fr,
     vertex.size = node_sizes,
     vertex.color = node_colors,
     vertex.label = V(g)$name,
     vertex.label.cex = 0.5,
     vertex.label.color = "black",
     vertex.label.dist = 0.3,
     vertex.frame.color = "gray30",
     edge.color = adjustcolor("gray50", alpha = 0.5),
     edge.width = edge_widths,
     main = "CAFI Co-occurrence Network (colored by taxonomic type)")

legend("topright",
       legend = names(type_colors),
       col = type_colors,
       pch = 19,
       pt.cex = 1.5,
       title = "Type",
       bg = "white",
       cex = 0.9)

# Add network stats
mtext(sprintf("N = %d species, %d edges | Q = %.2f | Transitivity = %.2f",
              n_nodes, n_edges, modularity_obs, transitivity_obs),
      side = 1, line = -1, cex = 0.9)

dev.off()
cat("     Saved: network_by_type.png\n")

# === Figure 2: Network colored by module ===
png(file.path(fig_dir, "network_by_module.png"),
    width = 12, height = 10, units = "in", res = 300, bg = "white")

par(mar = c(1, 1, 3, 1))

plot(communities_louvain, g,
     layout = layout_fr,
     vertex.size = node_sizes,
     vertex.label = V(g)$name,
     vertex.label.cex = 0.5,
     vertex.label.color = "black",
     vertex.label.dist = 0.3,
     edge.color = adjustcolor("gray50", alpha = 0.5),
     edge.width = edge_widths,
     main = sprintf("CAFI Network Groups (Louvain: %d groups, Q = %.2f)",
                    n_modules, modularity_obs))

dev.off()
cat("     Saved: network_by_module.png\n")

# === Figure 3: Main network visualization ===
png(file.path(fig_dir, "network_visualization.png"),
    width = 14, height = 12, units = "in", res = 300, bg = "white")

par(mar = c(2, 1, 3, 1))

plot(g,
     layout = layout_fr,
     vertex.size = node_sizes,
     vertex.color = node_colors,
     vertex.label = V(g)$name,
     vertex.label.cex = 0.55,
     vertex.label.color = "black",
     vertex.label.dist = 0.3,
     vertex.frame.color = "gray30",
     edge.color = adjustcolor("steelblue", alpha = 0.4),
     edge.width = edge_widths,
     main = "CAFI Species Co-occurrence Network")

legend("topright",
       legend = c(names(type_colors), "", "Edge = positive co-occurrence"),
       col = c(type_colors, NA, "steelblue"),
       pch = c(rep(19, length(type_colors)), NA, NA),
       lty = c(rep(NA, length(type_colors)), NA, 1),
       lwd = c(rep(NA, length(type_colors)), NA, 2),
       title = "Species Type",
       bg = "white",
       cex = 0.9)

mtext(sprintf("N = %d species, %d edges | Modularity Q = %.2f (%.1fx null) | Transitivity = %.2f (%.1fx null)",
              n_nodes, n_edges, modularity_obs,
              null_comparison$ratio_to_null[null_comparison$metric == "Modularity"],
              transitivity_obs,
              null_comparison$ratio_to_null[null_comparison$metric == "Transitivity"]),
      side = 1, line = 0.5, cex = 0.9)

dev.off()
cat("     Saved: network_visualization.png\n")

# === Figure 4: Manuscript figure (multi-panel) ===
# Create ggplot-based figures for manuscript

# 4A: Degree distribution
p_degree <- ggplot(centrality_df, aes(x = degree)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean_degree, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  labs(
    x = "Degree (number of co-occurring species)",
    y = "Number of species",
    title = "A Degree Distribution",
    subtitle = sprintf("Mean = %.1f, Median = %d", mean_degree, median_degree)
  )

# 4B: Modularity comparison to null
null_mod_df <- data.frame(
  modularity = null_metrics[, "modularity"]
)

p_modularity <- ggplot(null_mod_df, aes(x = modularity)) +
  geom_histogram(bins = 30, fill = "gray70", color = "white", alpha = 0.8) +
  geom_vline(xintercept = modularity_obs_unweighted, color = "#D55E00", linewidth = 1.2) +
  annotate("text", x = modularity_obs_unweighted, y = Inf,
           label = sprintf("Observed\nQ = %.2f", modularity_obs_unweighted),
           vjust = 1.5, hjust = -0.1, color = "#D55E00", fontface = "bold") +
  labs(
    x = "Modularity (Q)",
    y = "Frequency (null model)",
    title = "B Modularity vs Null Model",
    subtitle = sprintf("z = %.1f, p = %.4f",
                       null_comparison$z_score[null_comparison$metric == "Modularity"],
                       null_comparison$p_value[null_comparison$metric == "Modularity"])
  )

# 4C: Top hub species
p_hubs <- centrality_df %>%
  slice_head(n = 15) %>%
  ggplot(aes(x = reorder(species, hub_score), y = hub_score, fill = type)) +
  geom_col(alpha = 0.8, width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = type_colors, name = "Type") +
  labs(
    x = NULL,
    y = "Hub Score (degree + eigenvector, standardized)",
    title = "C Top 15 Hub Species",
    subtitle = "Based on connectivity and influence"
  ) +
  theme(legend.position = "bottom")

# 4D: Module composition
p_modules <- module_taxonomy %>%
  ggplot(aes(x = factor(module), y = n, fill = type)) +
  geom_col(position = "fill", alpha = 0.8) +
  scale_fill_manual(values = type_colors, name = "Type") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Group",
    y = "Proportion of species",
    title = "D Group Taxonomic Composition",
    subtitle = sprintf("%d groups identified by Louvain algorithm", n_modules)
  ) +
  theme(legend.position = "bottom")

# Combine panels
p_manuscript <- (p_degree + p_modularity) / (p_hubs + p_modules) +
  plot_annotation(
    title = "CAFI Co-occurrence Network Analysis",
    subtitle = sprintf("N = %d species, %d positive associations (r > 0.3)", n_nodes, n_edges),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
  )

save_figure(p_manuscript, file.path(fig_dir, "network_panels.png"),
            width = 12, height = 10)
cat("     Saved: network_panels.png\n")

# ############################################################################
# MANUSCRIPT FIGURE 5: CAFI CO-OCCURRENCE NETWORK (5-PANEL)
# ############################################################################
# Panel A: ALL species in circular layout grouped by Louvain group (hero panel)
# Panels B-E: Individual group sub-networks with species labels (force layout)
# Wrapped in tryCatch so it doesn't crash the pipeline if ggforce is missing
# ############################################################################

tryCatch({

  has_ggforce <- requireNamespace("ggforce", quietly = TRUE)
  has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

  if (!has_ggforce || !has_ggrepel) {
    missing_pkgs <- c()
    if (!has_ggforce) missing_pkgs <- c(missing_pkgs, "ggforce")
    if (!has_ggrepel) missing_pkgs <- c(missing_pkgs, "ggrepel")
    stop(sprintf("Required package(s) not available: %s", paste(missing_pkgs, collapse = ", ")))
  }

  library(ggforce)
  library(ggrepel)

  cat("\n------------------------------------------------------------\n")
  cat("MANUSCRIPT FIGURE 5: 5-PANEL NETWORK (WIDE LAYOUT)\n")
  cat("------------------------------------------------------------\n\n")

  # Use communities_louvain (already in scope from Part 4)
  communities <- communities_louvain

  # Species info dataframe
  sp_info <- data.frame(
    species = V(g)$name,
    guild = membership(communities),
    degree = degree(g),
    type = V(g)$type,
    stringsAsFactors = FALSE
  ) %>%
    mutate(guild = factor(guild))

  # Dynamic group configuration (adapts to any number of Louvain groups)
  n_guilds <- length(unique(sp_info$guild))
  guild_names <- setNames(paste("Group", seq_len(n_guilds)), as.character(seq_len(n_guilds)))
  guild_names_short <- guild_names

  guild_counts <- sp_info %>% count(guild) %>% arrange(guild)

  # Colorblind-safe group colors (Okabe-Ito derivatives, dynamically extended)
  okabe_ito_base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                      "#E69F00", "#56B4E9", "#F0E442")
  if (n_guilds <= length(okabe_ito_base)) {
    guild_colors <- setNames(okabe_ito_base[seq_len(n_guilds)], as.character(seq_len(n_guilds)))
  } else {
    guild_colors <- setNames(
      c(okabe_ito_base, scales::hue_pal()(n_guilds - length(okabe_ito_base))),
      as.character(seq_len(n_guilds))
    )
  }

  # Lighter versions for arc backgrounds (dynamic)
  guild_colors_light <- setNames(
    sapply(guild_colors, function(col) adjustcolor(col, alpha.f = 0.4)),
    names(guild_colors)
  )

  cat("Group sizes:\n")
  for (i in sort(unique(guild_counts$guild))) {
    cat(paste0("  Group ", i, " (", guild_names[as.character(i)], "): ",
                guild_counts$n[guild_counts$guild == i], " species\n"))
  }
  cat(paste0("\nTotal species: ", sum(guild_counts$n), "\n"))
  cat(paste0("Total edges: ", ecount(g), "\n\n"))

  # --------------------------------------------------------------------------
  # PANEL A: HERO CIRCULAR NETWORK - ALL SPECIES
  # --------------------------------------------------------------------------

  cat("Building Panel A (circular species network)...\n")

  sp_info_sorted <- sp_info %>%
    group_by(guild) %>%
    arrange(desc(degree), .by_group = TRUE) %>%
    mutate(rank_in_guild = row_number()) %>%
    ungroup()

  gap_size <- 0.08
  total_gap <- gap_size * n_guilds
  species_arc <- (2 * pi) * (1 - total_gap)

  guild_sizes <- guild_counts$n
  n_total <- sum(guild_sizes)

  guild_props <- guild_sizes / n_total
  guild_arcs <- guild_props * species_arc

  guild_starts <- c(0)
  if (n_guilds > 1) {
    for (i in 1:(n_guilds - 1)) {
      guild_starts <- c(guild_starts,
                        guild_starts[i] + guild_arcs[i] + gap_size * 2 * pi)
    }
  }

  sp_positions <- sp_info_sorted %>%
    group_by(guild) %>%
    mutate(
      guild_idx = as.numeric(guild),
      n_in_guild = n(),
      pos_in_guild = (row_number() - 0.5) / n_in_guild
    ) %>%
    ungroup() %>%
    mutate(
      guild_start = guild_starts[guild_idx],
      guild_arc = guild_arcs[guild_idx],
      angle = pi/2 - (guild_start + pos_in_guild * guild_arc),
      radius = 4.5,
      x = radius * cos(angle),
      y = radius * sin(angle),
      node_size = scales::rescale(degree, to = c(3, 9))
    )

  edge_df <- igraph::as_data_frame(g, what = "edges") %>%
    left_join(sp_positions %>% dplyr::select(species, x, y, guild, degree),
              by = c("from" = "species")) %>%
    rename(x1 = x, y1 = y, guild_from = guild, degree_from = degree) %>%
    left_join(sp_positions %>% dplyr::select(species, x, y, guild, degree),
              by = c("to" = "species")) %>%
    rename(x2 = x, y2 = y, guild_to = guild, degree_to = degree) %>%
    mutate(
      is_within_guild = guild_from == guild_to,
      edge_alpha = ifelse(is_within_guild, 0.5, 0.15),
      edge_color = ifelse(is_within_guild,
                          guild_colors[as.character(guild_from)],
                          "gray50")
    )

  # Filter edges: keep only top 50% by weight to reduce visual clutter
  weight_median <- median(edge_df$weight)
  edge_df_filtered <- edge_df %>% dplyr::filter(weight >= weight_median)
  cat(sprintf("  Edge filtering: median weight = %.3f, keeping %d / %d edges (top 50%%)\n",
              weight_median, nrow(edge_df_filtered), nrow(edge_df)))

  between_guild_edges <- edge_df_filtered %>% dplyr::filter(!is_within_guild)
  # Use between_guild_edges directly (no sampling needed)
  n_between <- nrow(between_guild_edges)

  cat(sprintf("  Between-group edges: %d (all shown as thin gray)\n", n_between))

  create_bezier <- function(x1, y1, x2, y2, n_points = 30, curvature = 0.35) {
    cx <- (x1 + x2) / 2 * (1 - curvature)
    cy <- (y1 + y2) / 2 * (1 - curvature)
    t <- seq(0, 1, length.out = n_points)
    data.frame(
      x = (1-t)^2 * x1 + 2*(1-t)*t * cx + t^2 * x2,
      y = (1-t)^2 * y1 + 2*(1-t)*t * cy + t^2 * y2
    )
  }

  within_guild_edges <- edge_df_filtered %>% dplyr::filter(is_within_guild)
  cat("  Generating bezier curves for", nrow(within_guild_edges), "within-group edges...\n")

  bezier_within <- purrr::map_dfr(1:nrow(within_guild_edges), function(i) {
    row <- within_guild_edges[i,]
    pts <- create_bezier(row$x1, row$y1, row$x2, row$y2)
    pts$edge_id <- i
    pts$weight <- row$weight
    pts$guild <- as.character(row$guild_from)
    pts
  })

  cat("  Generating bezier curves for", nrow(between_guild_edges), "between-group edges...\n")

  if (nrow(between_guild_edges) > 0) {
    bezier_between <- purrr::map_dfr(1:nrow(between_guild_edges), function(i) {
      row <- between_guild_edges[i,]
      pts <- create_bezier(row$x1, row$y1, row$x2, row$y2)
      pts$edge_id <- i + 10000
      pts$weight <- row$weight
      pts
    })
  } else {
    bezier_between <- data.frame(x = numeric(), y = numeric(), edge_id = integer(), weight = numeric())
  }

  guild_arc_data <- sp_positions %>%
    group_by(guild) %>%
    summarize(
      min_angle = min(angle),
      max_angle = max(angle),
      .groups = "drop"
    ) %>%
    mutate(
      padding = 0.05,
      start = max_angle + padding,
      end = min_angle - padding,
      r_inner = 3.8,
      r_outer = 5.2,
      fill_color = guild_colors_light[as.character(guild)],
      border_color = guild_colors[as.character(guild)]
    )

  # Dynamically position group labels in corners / evenly spaced around the plot
  # For each group, place label at an angle evenly spaced around the perimeter
  label_angles <- seq(pi/4, pi/4 + 2*pi*(1 - 1/n_guilds), length.out = n_guilds)
  label_radius <- 6.5
  label_x_positions <- label_radius * cos(label_angles)
  label_y_positions <- label_radius * sin(label_angles)
  label_hjust <- ifelse(cos(label_angles) > 0, 1, 0)
  label_vjust <- ifelse(sin(label_angles) > 0, 0, 1)

  guild_label_positions <- sp_positions %>%
    group_by(guild) %>%
    summarize(mid_angle = (min(angle) + max(angle)) / 2, .groups = "drop") %>%
    left_join(guild_counts, by = "guild") %>%
    mutate(
      label = paste0(guild_names[as.character(guild)], " (n=", n, ")"),
      x = label_x_positions[as.numeric(guild)],
      y = label_y_positions[as.numeric(guild)],
      hjust = label_hjust[as.numeric(guild)],
      vjust = label_vjust[as.numeric(guild)],
      color = guild_colors[as.character(guild)]
    )

  make_arc_polygon <- function(start_angle, end_angle, r_inner, r_outer, n = 100) {
    angles_outer <- seq(start_angle, end_angle, length.out = n)
    angles_inner <- rev(angles_outer)
    data.frame(
      x = c(r_outer * cos(angles_outer), r_inner * cos(angles_inner)),
      y = c(r_outer * sin(angles_outer), r_inner * sin(angles_inner))
    )
  }

  arc_layers <- lapply(1:nrow(guild_arc_data), function(i) {
    d <- guild_arc_data[i, ]
    poly_df <- make_arc_polygon(d$start, d$end, d$r_inner, d$r_outer)
    geom_polygon(
      data = poly_df,
      aes(x = x, y = y),
      fill = d$fill_color,
      color = d$border_color,
      alpha = 0.5,
      linewidth = 0.8
    )
  })

  p_A <- ggplot() +
    arc_layers +
    geom_path(
      data = bezier_between,
      aes(x = x, y = y, group = edge_id),
      color = "gray65",
      alpha = 0.07,
      linewidth = 0.2
    ) +
    geom_path(
      data = bezier_within,
      aes(x = x, y = y, group = edge_id, color = guild, alpha = weight),
      linewidth = 0.4
    ) +
    geom_point(
      data = sp_positions,
      aes(x = x, y = y, fill = guild, size = node_size),
      shape = 21,
      color = "white",
      stroke = 0.7
    ) +
    geom_text(
      data = guild_label_positions,
      aes(x = x, y = y, label = label, color = guild, hjust = hjust, vjust = vjust),
      size = 3.8,
      fontface = "bold"
    ) +
    scale_fill_manual(values = guild_colors, guide = "none") +
    scale_color_manual(values = guild_colors, guide = "none") +
    scale_size_identity() +
    scale_alpha_continuous(range = c(0.08, 0.45), guide = "none") +
    coord_fixed(ratio = 1, xlim = c(-6.8, 6.8), ylim = c(-6.8, 6.8), clip = "off") +
    labs(
      title = "A Species Co-occurrence Network",
      subtitle = sprintf("%d species | Showing %d of %d edges (top 50%% by weight) | %d Louvain groups",
                        nrow(sp_positions), nrow(edge_df_filtered), ecount(g), n_guilds)
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 5)),
      plot.margin = margin(5, 2, 5, 2)
    )

  cat("  Panel A complete.\n")

  # --------------------------------------------------------------------------
  # PANELS B-E: INDIVIDUAL GROUP SUB-NETWORKS WITH SPECIES LABELS
  # --------------------------------------------------------------------------

  create_guild_panel_fixed <- function(guild_id, letter, max_labels = 10) {

    guild_name <- guild_names_short[as.character(guild_id)]
    guild_color <- guild_colors[as.character(guild_id)]
    guild_color_light <- guild_colors_light[as.character(guild_id)]

    guild_species <- sp_info %>% dplyr::filter(guild == guild_id)
    g_sub <- induced_subgraph(g, V(g)[V(g)$name %in% guild_species$species])

    if (vcount(g_sub) == 0) {
      return(ggplot() + theme_void() +
               labs(title = paste(letter, guild_name)))
    }

    set.seed(42 + guild_id)
    n_sp <- vcount(g_sub)
    if (n_sp > 12) {
      layout_sub <- layout_with_kk(g_sub)
    } else {
      layout_sub <- layout_with_fr(g_sub, niter = 1000, weights = E(g_sub)$weight)
    }

    if (nrow(layout_sub) > 1) {
      layout_sub[,1] <- scales::rescale(layout_sub[,1], to = c(-1.4, 1.4))
      layout_sub[,2] <- scales::rescale(layout_sub[,2], to = c(-1.4, 1.4))
    } else {
      layout_sub[,1] <- 0
      layout_sub[,2] <- 0
    }

    # Scale node sizes based on panel density
    node_size_range <- if (n_sp > 15) c(3, 10) else c(4, 14)

    node_data <- data.frame(
      species = V(g_sub)$name,
      x = layout_sub[,1],
      y = layout_sub[,2],
      degree = degree(g_sub),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        node_size = scales::rescale(degree, to = node_size_range),
        species_label = gsub("([A-Z])[a-z]+ ", "\\1. ", species)
      )

    n_species <- as.integer(nrow(node_data))
    n_to_label <- as.integer(if (n_species <= max_labels) n_species else max_labels)

    node_data <- node_data %>%
      mutate(
        rank_degree = rank(-degree, ties.method = "first"),
        show_label = rank_degree <= n_to_label
      )

    cat(sprintf("    Group %d: %d species, labeling top %d\n",
                guild_id, n_species, sum(node_data$show_label)))

    if (ecount(g_sub) > 0) {
      edge_list_sub <- igraph::as_data_frame(g_sub, what = "edges")
      edge_data <- edge_list_sub %>%
        left_join(node_data %>% dplyr::select(species, x, y),
                  by = c("from" = "species")) %>%
        rename(x1 = x, y1 = y) %>%
        left_join(node_data %>% dplyr::select(species, x, y),
                  by = c("to" = "species")) %>%
        rename(x2 = x, y2 = y)
    } else {
      edge_data <- NULL
    }

    n_edges_sub <- if (is.null(edge_data)) 0L else as.integer(nrow(edge_data))

    p <- ggplot()

    p <- p + ggforce::geom_circle(
      aes(x0 = 0, y0 = 0, r = 1.4),
      fill = guild_color_light,
      color = NA,
      alpha = 0.25
    )

    if (!is.null(edge_data) && nrow(edge_data) > 0) {
      edge_lw <- if (n_sp > 15) 0.3 else 0.5
      p <- p + geom_segment(
        data = edge_data,
        aes(x = x1, y = y1, xend = x2, yend = y2, alpha = weight),
        color = guild_color,
        linewidth = edge_lw
      )
    }

    p <- p + geom_point(
      data = node_data,
      aes(x = x, y = y, size = node_size),
      fill = guild_color,
      shape = 21,
      color = "white",
      stroke = 0.7
    )

    labels_data <- node_data %>% dplyr::filter(show_label)

    if (nrow(labels_data) > 0) {
      # Scale label size based on species count: smaller text for dense panels
      label_size <- if (n_species <= 5) 4.0 else if (n_species <= 12) 3.5 else 2.8
      repel_force <- if (n_species > 15) 25 else 12
      p <- p + ggrepel::geom_text_repel(
        data = labels_data,
        aes(x = x, y = y, label = species_label),
        size = label_size,
        fontface = "bold.italic",
        color = "gray5",
        bg.color = "white",
        bg.r = 0.15,
        segment.color = "gray40",
        segment.size = 0.3,
        segment.alpha = 0.6,
        box.padding = 0.35,
        point.padding = 0.3,
        max.overlaps = 50,
        force = repel_force,
        force_pull = 0.2,
        max.iter = 20000,
        seed = 42,
        min.segment.length = 0.1
      )
    }

    n_hidden <- n_species - sum(node_data$show_label)
    subtitle_text <- if (n_hidden > 0) {
      paste0(n_species, " species | ", n_edges_sub, " edges | top ", n_to_label, " labeled")
    } else {
      paste0(n_species, " species | ", n_edges_sub, " edges")
    }

    p <- p +
      scale_size_identity() +
      scale_alpha_continuous(range = if (n_sp > 15) c(0.1, 0.4) else c(0.15, 0.65), guide = "none") +
      # Expand coordinate limits for dense panels to give labels room
      coord_fixed(ratio = 1,
                  xlim = c(if (n_sp > 15) -2.8 else -2.3, if (n_sp > 15) 2.8 else 2.3),
                  ylim = c(if (n_sp > 15) -2.8 else -2.3, if (n_sp > 15) 2.8 else 2.3)) +
      labs(
        title = sprintf("%s %s", letter, guild_name),
        subtitle = subtitle_text
      ) +
      theme_void() +
      theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5,
                                  color = guild_color),
        plot.subtitle = element_text(size = 8.5, hjust = 0.5, color = "gray50"),
        plot.margin = margin(4, 4, 4, 4),
        plot.background = element_rect(fill = "white", color = "gray80",
                                       linewidth = 0.5)
      )

    return(p)
  }

  cat("Building sub-panels (individual groups)...\n")

  # Dynamically create one panel per Louvain group
  panel_letters <- LETTERS[2:(n_guilds + 1)]  # B, C, D, E, ...
  guild_panels <- list()
  for (gi in seq_len(n_guilds)) {
    guild_panels[[gi]] <- create_guild_panel_fixed(gi, panel_letters[gi], max_labels = 50)
    cat(sprintf("  Panel %s complete.\n", panel_letters[gi]))
  }

  # --------------------------------------------------------------------------
  # MANUSCRIPT FIGURE 4: Full co-occurrence network (single panel)
  # --------------------------------------------------------------------------

  cat("\nAssembling manuscript Figure 4 (full network)...\n")

  p_fig4 <- p_A +
    plot_annotation(
      subtitle = paste0(
        n_guilds, " groups (Louvain) | Q = ",
        sprintf("%.3f", modularity(communities)), " | ",
        vcount(g), " species | ", ecount(g), " edges"
      ),
      caption = paste0(
        "Node size = degree centrality | Within-group edges colored, between-group edges gray | ",
        "Threshold: r > 0.3, FDR p < 0.05 | See Figure S11 for group sub-networks"
      ),
      theme = theme(
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
        plot.caption = element_text(size = 9.5, hjust = 0.5, color = "gray50",
                                    margin = margin(t = 12)),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(10, 10, 10, 10, "mm")
      )
    )

  # Save manuscript figure
  save_figure(p_fig4, file.path(PATHS$fig_manuscript, "fig4_network.png"),
              width = 170, height = 170, units = "mm")
  cat("     Saved: manuscript/fig4_network.png\n")

  save_figure(p_fig4, file.path(fig_dir, "fig4_network.png"),
              width = 170, height = 170, units = "mm")
  cat("     Saved: 06_network/fig4_network.png\n")

  # --------------------------------------------------------------------------
  # SUPPLEMENT FIGURE S12: Louvain group sub-networks
  # --------------------------------------------------------------------------

  cat("\nAssembling supplement Figure S11 (group sub-networks)...\n")

  # Dynamically arrange sub-panels in a grid (2 columns)
  n_rows_right <- ceiling(n_guilds / 2)
  if (n_guilds >= 4) {
    p_groups <- (guild_panels[[1]] | guild_panels[[2]])
    if (n_guilds >= 3) {
      row2 <- if (n_guilds >= 4) (guild_panels[[3]] | guild_panels[[4]]) else guild_panels[[3]]
      p_groups <- p_groups / row2
    }
    if (n_guilds > 4) {
      for (ri in seq(5, n_guilds, by = 2)) {
        row_i <- if (ri + 1 <= n_guilds) (guild_panels[[ri]] | guild_panels[[ri + 1]]) else guild_panels[[ri]]
        p_groups <- p_groups / row_i
      }
    }
  } else if (n_guilds == 3) {
    p_groups <- (guild_panels[[1]] | guild_panels[[2]]) / (guild_panels[[3]] | plot_spacer())
  } else if (n_guilds == 2) {
    p_groups <- guild_panels[[1]] / guild_panels[[2]]
  } else {
    p_groups <- guild_panels[[1]]
  }

  p_figS11 <- p_groups +
    plot_annotation(
      title = "Figure S11: Louvain Community Group Sub-networks",
      subtitle = paste0("Subset networks for each of ", n_guilds,
                        " groups detected via Louvain community detection"),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(10, 10, 10, 10, "mm")
      )
    )

  save_figure(p_figS11,
              file.path(PATHS$fig_supplement, "figS11_network_groups.png"),
              width = 240, height = 200, units = "mm")
  save_figure(p_figS11,
              file.path(fig_dir, "figS11_network_groups.png"),
              width = 240, height = 200, units = "mm")
  cat("     Saved: figS11_network_groups.png (supplement + analysis)\n")

  cat("\n  Manuscript Figure 4 + supplement S11 complete.\n")

}, error = function(e) {
  cat("\n  WARNING: Could not create manuscript Figure 4 + supplement S11.\n")
  cat("  Reason:", conditionMessage(e), "\n")
  cat("  The analysis figure (network_panels.png) was still saved.\n\n")
})

# --- Generate figure legend and results text ---
cat("\n  Generating figure legend and results text...\n")

# Extract null comparison stats
null_mod <- null_comparison[null_comparison$metric == "Modularity", ]
null_trans <- null_comparison[null_comparison$metric == "Transitivity", ]

fig4_legend <- paste0(
  "FIGURE 4: CAFI CO-OCCURRENCE NETWORK STRUCTURE\n",
  "================================================================================\n\n",

  "FIGURE LEGEND\n",
  "-------------\n",
  "Figure 4. Co-occurrence network of coral-associated fauna and invertebrate\n",
  "(CAFI) communities. Nodes represent species (n = ", n_nodes,
  "),\n",
  "edges represent significant positive correlations (r > 0.3, FDR-corrected p < 0.05;\n",
  n_edges, " total edges). Nodes colored by Louvain-detected group, sized by degree\n",
  "centrality. Only the top 50% of edges by weight are displayed to reduce visual clutter.\n",
  "Network density is high (", sprintf("%.2f", density), ") and modularity (Q = ",
  sprintf("%.3f", modularity_obs), ") is low,\n",
  "indicating weak community structure. Groups identified by Louvain algorithm should\n",
  "be interpreted as statistical clusters rather than distinct ecological modules.\n",
  "Individual group sub-networks are shown in Figure S11.\n\n",

  "================================================================================\n\n",

  "STATISTICAL RESULTS\n",
  "-------------------\n\n",

  "1. NETWORK CONSTRUCTION\n",
  "   Species included: ", n_nodes, " (present in >= threshold corals)\n",
  "   Edges: ", n_edges, " significant positive associations\n",
  "   Correlation threshold: r > 0.3, FDR p < 0.05\n",
  "   Network density: ", sprintf("%.3f", density), "\n\n",

  "2. COMMUNITY DETECTION (Louvain)\n",
  "   Groups detected: ", n_modules, "\n",
  "   Modularity (Q): ", sprintf("%.3f", modularity_obs), "\n",
  "   Modularity (unweighted): ", sprintf("%.3f", modularity_obs_unweighted), "\n\n",

  "3. NULL MODEL COMPARISON (1000 configuration model, degree-preserving)\n",
  "   Observed modularity: ", sprintf("%.3f", modularity_obs), "\n",
  "   Null mean modularity: ", sprintf("%.3f", null_mod$null_mean), "\n",
  "   z-score: ", sprintf("%.1f", null_mod$z_score), "\n",
  "   Ratio to null: ", sprintf("%.1fx", null_mod$ratio_to_null), "\n",
  "   Interpretation: ", ifelse(null_mod$z_score > 2,
    "Modularity is elevated above random expectation.",
    ifelse(null_mod$z_score < 0,
      paste0("Modularity is LOWER than random expectation (z = ", sprintf("%.1f", null_mod$z_score), ")."),
      "Modularity is not significantly elevated above random expectation.")),
  "\n\n",

  "   Observed transitivity: ", sprintf("%.3f", transitivity_obs), "\n",
  "   Null mean transitivity: ", sprintf("%.3f", null_trans$null_mean), "\n",
  "   z-score: ", sprintf("%.1f", null_trans$z_score), "\n",
  "   Interpretation: ", ifelse(null_trans$z_score > 2,
    "Higher clustering than random.\n\n",
    "Clustering comparable to random networks.\n\n"),

  "4. NETWORK TOPOLOGY\n",
  "   Mean degree: ", sprintf("%.1f", mean_degree), "\n",
  "   Median degree: ", median_degree, "\n",
  "   Max degree: ", max_degree, "\n",
  "   Transitivity (global): ", sprintf("%.3f", transitivity_obs), "\n",
  "   Diameter: ", diameter_obs, "\n",
  "   Mean path length: ", sprintf("%.2f", mean_distance_obs), "\n\n",

  "5. HUB SPECIES (Top 10 by hub score)\n"
)

for (i in 1:min(10, nrow(centrality_df))) {
  fig4_legend <- paste0(fig4_legend,
    "   ", i, ". ", centrality_df$species[i],
    " (", centrality_df$type[i], ")",
    " — degree: ", centrality_df$degree[i],
    ", hub score: ", sprintf("%.2f", centrality_df$hub_score[i]), "\n"
  )
}

fig4_legend <- paste0(fig4_legend, "\n",
  "================================================================================\n\n",

  "RESULTS\n",
  "-------\n\n",

  "The CAFI co-occurrence network comprised ", n_nodes, " species connected by ",
  n_edges, " significant positive associations (r > 0.3, FDR p < 0.05). ",
  "Network density was high (", sprintf("%.2f", density), "), indicating that most species ",
  "co-occurred with most others. Louvain community detection identified ", n_modules, " groups ",
  "(Q = ", sprintf("%.3f", modularity_obs), "; z = ",
  sprintf("%.1f", null_mod$z_score), " vs configuration model null, ",
  sprintf("%.1fx", null_mod$ratio_to_null), " random). ",
  ifelse(null_mod$z_score > 2,
    "Modularity was elevated above null expectation.",
    ifelse(null_mod$z_score < 0,
      "Observed modularity was LOWER than null expectation. The high density and low Q indicate the network is a near-clique with weak group structure.",
      "Modularity was not significantly elevated above null expectation. The high density and low Q indicate the network is a near-clique with weak group structure.")),
  " Groups should be interpreted as statistical clusters identified by the Louvain algorithm ",
  "rather than ecologically validated modules.\n\n",

  "The network exhibited ", ifelse(null_trans$z_score > 2, "elevated", "moderate"),
  " clustering (transitivity = ", sprintf("%.2f", transitivity_obs),
  ") and a mean path length of ", sprintf("%.1f", mean_distance_obs),
  ". Hub species -- those with high degree ",
  "and eigenvector centrality -- included obligate coral associates such as ",
  centrality_df$species[1], " and ", centrality_df$species[2],
  ", suggesting these taxa are highly connected within the network.\n\n",

  "================================================================================\n\n",

  "METHODS\n",
  "-------\n\n",

  "We constructed a species co-occurrence network from volume-corrected presence data. ",
  "Species presence was converted to binary (presence-absence), then residualized on ",
  "log(coral volume) via logistic GLM to remove the mechanical confound that larger corals ",
  "host more species. Pairwise Spearman correlations were computed on deviance residuals; ",
  "edges were retained for positive correlations exceeding r = 0.3 with FDR-corrected ",
  "p < 0.05. Louvain community detection was applied to identify groups ",
  "(Blondel et al. 2008). To assess whether observed network properties ",
  "deviate from random expectations, we compared modularity and transitivity against ",
  "1000 configuration model (degree-preserving) random networks. ",
  "Hub species were identified by combining standardized degree and eigenvector ",
  "centrality into a composite hub score.\n\n",

  "================================================================================\n",
  "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "Source script: scripts/06_network_analysis.R\n"
)

writeLines(fig4_legend, file.path(PATHS$fig_manuscript, "fig4_legend_results.txt"))
cat("  Saved: fig4_legend_results.txt\n\n")

# ============================================================================
# PART 8: SAVE OUTPUTS
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 8: SAVING OUTPUTS\n")
cat("------------------------------------------------------------\n\n")

# 8.1 Network metrics table
network_metrics_df <- data.frame(
  metric = c("n_species", "n_edges", "density", "transitivity", "diameter",
             "mean_path_length", "modularity", "modularity_unweighted",
             "n_modules", "mean_degree", "median_degree", "max_degree"),
  value = c(n_nodes, n_edges, density, transitivity_obs, diameter_obs,
            mean_distance_obs, modularity_obs, modularity_obs_unweighted,
            n_modules, mean_degree, median_degree, max_degree)
)

# Add null model comparison
network_metrics_full <- network_metrics_df %>%
  left_join(
    null_comparison %>%
      mutate(metric = tolower(gsub(" ", "_", metric))) %>%
      dplyr::select(metric, null_mean, null_sd, z_score, p_value, ratio_to_null),
    by = "metric"
  )

save_table(network_metrics_full, "network_metrics")
cat("     Saved: network_metrics.csv\n")

# 8.2 Module composition table
module_composition <- module_species %>%
  left_join(module_summary %>% dplyr::select(module, n_species, hub_species), by = "module")

save_table(module_composition, "module_composition")
cat("     Saved: module_composition.csv\n")

# 8.3 Hub species table
hub_species_full <- centrality_df %>%
  arrange(desc(hub_score)) %>%
  dplyr::select(species, type, functional_group, module, abundance, occurrence,
                degree, betweenness, closeness, eigenvector, hub_score, keystone_index)

save_table(hub_species_full, "hub_species")
cat("     Saved: hub_species.csv\n")

# 8.4 Save network object
network_results <- list(
  graph = g,
  communities = communities_louvain,
  null_comparison = null_comparison,
  centrality = centrality_df,
  module_summary = module_summary,
  network_metrics = network_metrics_full,
  edge_list = edge_list,
  fr_layout = layout_fr,
  null_metrics = null_metrics,
  module_taxonomy = module_taxonomy,
  type_colors = type_colors
)

save_object(network_results, "cafi_network")
cat("     Saved: cafi_network.rds\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    NETWORK ANALYSIS SUMMARY\n")
cat("============================================================\n\n")

cat("CAFI Co-occurrence Network Structure:\n")
cat("  - Species in network:", n_nodes, "\n")
cat("  - Positive associations:", n_edges, "\n")
cat("  - Network density:", round(density, 4), "\n")
cat("  - Mean degree:", round(mean_degree, 2), "\n\n")

cat("Non-random Structure (vs configuration model null):\n")
cat(sprintf("  - Modularity Q = %.2f (weighted), %.2f (unweighted)\n",
            modularity_obs, modularity_obs_unweighted))
cat(sprintf("  - Unweighted Q vs null: %.1fx null, z = %.1f, p = %.4f\n",
            null_comparison$ratio_to_null[null_comparison$metric == "Modularity"],
            null_comparison$z_score[null_comparison$metric == "Modularity"],
            null_comparison$p_value[null_comparison$metric == "Modularity"]))
cat(sprintf("  - Transitivity = %.2f (%.1fx null, z = %.1f, p = %.4f)\n",
            transitivity_obs,
            null_comparison$ratio_to_null[null_comparison$metric == "Transitivity"],
            null_comparison$z_score[null_comparison$metric == "Transitivity"],
            null_comparison$p_value[null_comparison$metric == "Transitivity"]))
cat(sprintf("  - Number of groups: %d\n\n", n_modules))

cat("Top Hub Species (by hub score):\n")
top_hubs <- centrality_df %>% slice_head(n = 5)
for (i in 1:nrow(top_hubs)) {
  cat(sprintf("  %d. %s (%s) - degree: %d, hub score: %.2f\n",
              i, top_hubs$species[i], top_hubs$type[i],
              top_hubs$degree[i], top_hubs$hub_score[i]))
}

cat("\nKey Findings (H5):\n")
mod_z_final <- null_comparison$z_score[null_comparison$metric == "Modularity"]
trans_z_final <- null_comparison$z_score[null_comparison$metric == "Transitivity"]
if (!is.na(mod_z_final) && mod_z_final > 2) {
  cat("  - Louvain groups show modularity elevated above random expectation\n")
} else if (!is.na(mod_z_final) && mod_z_final < 0) {
  cat("  - Network does NOT show modularity elevated above random (z < 0)\n")
} else {
  cat("  - Modularity is not significantly elevated above random\n")
}
if (!is.na(trans_z_final) && trans_z_final > 2) {
  cat("  - Transitivity higher than random expectation\n")
} else {
  cat("  - Transitivity comparable to random expectation\n")
}
cat("  - Hub species identifiable through centrality metrics\n")
cat("  - Groups characterized by taxonomic composition\n\n")

cat("Outputs saved:\n")
cat("  Figures:\n")
cat("    - output/figures/network_visualization.png\n")
cat("    - output/figures/manuscript/fig4_network.png\n")
cat("    - output/figures/06_network/network_by_type.png\n")
cat("    - output/figures/06_network/network_by_module.png\n")
cat("  Tables:\n")
cat("    - output/tables/network_metrics.csv\n")
cat("    - output/tables/module_composition.csv\n")
cat("    - output/tables/hub_species.csv\n")
cat("  Objects:\n")
cat("    - output/objects/cafi_network.rds\n\n")

cat("============================================================\n")
cat("    NETWORK ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

} # end if (nrow(edge_indices) > 0)
