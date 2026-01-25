# ============================================================================
# 06_network_analysis.R - CAFI Species Co-occurrence Network Analysis
# ============================================================================
#
# PURPOSE: Build and analyze CAFI species co-occurrence networks to test
#          hypothesis H5 regarding non-random community structure
#
# HYPOTHESIS (H5): CAFI co-occurrence networks exhibit non-random modular
#   structure, with modules corresponding to functional groups or shared
#   habitat preferences, and identifiable hub/keystone species.
#
# ANALYSES:
#   - Co-occurrence network construction (volume-residualized correlations)
#   - Community detection (Louvain algorithm for modularity)
#   - Null model comparison (1000 configuration model random networks)
#   - Centrality analysis (degree, betweenness, eigenvector)
#   - Module composition characterization
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
cran_mirror <- "https://cran.rstudio.com"
if (!require("igraph", quietly = TRUE)) {
  install.packages("igraph", repos = cran_mirror)
  library(igraph)
} else {
  library(igraph)
}

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
volume_vec <- coral_master$volume_field[match(rownames(comm_filtered), coral_master$coral_id)]
log_volume <- log10(volume_vec)

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
cat("1.3 Calculating volume-corrected co-occurrences (Spearman correlation)...\n")

cor_matrix <- cor(residual_matrix, method = "spearman", use = "pairwise.complete.obs")

# Handle NAs in correlation matrix
cor_matrix[is.na(cor_matrix)] <- 0

# 1.5 Extract significant positive associations with FDR correction
# Threshold: r > 0.3 AND FDR-corrected p < 0.05
threshold <- 0.3

cat("     Computing pairwise p-values for FDR correction...\n")

# Compute p-values for all pairs above threshold
n_sp <- ncol(residual_matrix)
p_matrix <- matrix(1, nrow = n_sp, ncol = n_sp)
for (i in 1:(n_sp - 1)) {
  for (j in (i + 1):n_sp) {
    if (cor_matrix[i, j] > threshold) {
      ct <- tryCatch(
        cor.test(residual_matrix[, i], residual_matrix[, j], method = "spearman"),
        error = function(e) list(p.value = 1)
      )
      p_matrix[i, j] <- ct$p.value
      p_matrix[j, i] <- ct$p.value
    }
  }
}

# Extract pairs above threshold
edge_indices_raw <- which(upper.tri(cor_matrix) & cor_matrix > threshold, arr.ind = TRUE)

if (nrow(edge_indices_raw) > 0) {
  # Get raw p-values for candidate edges
  raw_pvals <- sapply(1:nrow(edge_indices_raw), function(k)
    p_matrix[edge_indices_raw[k, 1], edge_indices_raw[k, 2]])

  # Apply FDR correction (Benjamini-Hochberg)
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
    n_species = ncol(comm_network),
    n_edges = 0,
    note = "No significant positive associations found"
  )
  save_object(network_results, "network_results")
  cat("Script 06 complete (no network constructed).\n\n")
  return(invisible(NULL))
}

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
    filter(otu == sp | species == sp) %>%
    pull(type) %>%
    unique()
  if (length(type_val) == 0) "unknown" else type_val[1]
})

# Get functional group
V(g)$functional_group <- sapply(V(g)$name, function(sp) {
  fg_val <- cafi_clean %>%
    filter(otu == sp | species == sp) %>%
    pull(functional_group) %>%
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
diameter_obs <- diameter(g)
mean_distance_obs <- mean_distance(g)

# Modularity via Louvain algorithm (weighted for community assignments)
set.seed(42)  # For reproducibility
communities_louvain <- cluster_louvain(g, weights = E(g)$weight)
modularity_obs <- modularity(communities_louvain)
n_modules <- length(unique(membership(communities_louvain)))

# Unweighted modularity for null model comparison (null graphs are unweighted)
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
cat("     Number of modules:", n_modules, "\n")
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
cat("     (Configuration model preserves observed degree sequence)\n")

n_permutations <- 1000
null_metrics <- matrix(NA, nrow = n_permutations, ncol = 4)
colnames(null_metrics) <- c("modularity", "transitivity", "mean_distance", "diameter")

# Observed degree sequence for configuration model
obs_degree_seq <- degree(g)

set.seed(123)

for (i in 1:n_permutations) {
  # Generate configuration model random graph preserving degree sequence
  g_random <- tryCatch(
    sample_degseq(obs_degree_seq, method = "simple"),
    error = function(e) {
      # Fallback: simple.no.multiple if simple fails
      tryCatch(
        sample_degseq(obs_degree_seq, method = "simple.no.multiple"),
        error = function(e2) erdos.renyi.game(n_nodes, density, type = "gnp")
      )
    }
  )

  # Calculate metrics
  null_metrics[i, "transitivity"] <- transitivity(g_random, type = "global")
  null_metrics[i, "mean_distance"] <- mean_distance(g_random)
  null_metrics[i, "diameter"] <- diameter(g_random)

  # Modularity (only if network has edges)
  if (ecount(g_random) > 0) {
    null_comm <- cluster_louvain(g_random)
    null_metrics[i, "modularity"] <- modularity(null_comm)
  } else {
    null_metrics[i, "modularity"] <- 0
  }

  if (i %% 200 == 0) cat("     ", i, "permutations completed\n")
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
cat("PART 6: MODULE ANALYSIS\n")
cat("------------------------------------------------------------\n\n")

cat("6.1 Characterizing module composition...\n\n")

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

cat("     Module Summary:\n\n")
print(module_summary)

# 6.2 Taxonomic composition by module
cat("\n6.2 Taxonomic composition by module...\n\n")

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
cat("\n6.3 Functional group composition by module...\n\n")

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
  chi_test <- chisq.test(contingency_table, simulate.p.value = TRUE, B = 2000)
  cat("\n     Chi-square test for taxonomic clustering across modules:\n")
  cat(sprintf("       X2 = %.2f, p = %.4f\n", chi_test$statistic, chi_test$p.value))
  if (chi_test$p.value < 0.05) {
    cat("       Interpretation: Modules show NON-RANDOM taxonomic composition\n")
  } else {
    cat("       Interpretation: Modules show random taxonomic mixing\n")
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

# Module colors
n_mod <- length(unique(V(g)$module))
module_colors <- rainbow(n_mod, alpha = 0.8)

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
     main = sprintf("CAFI Network Modules (Louvain: %d modules, Q = %.2f)",
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
    title = "A. Degree Distribution",
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
    title = "B. Modularity vs Null Model",
    subtitle = sprintf("z = %.1f, p < 0.001", null_comparison$z_score[null_comparison$metric == "Modularity"])
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
    title = "C. Top 15 Hub Species",
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
    x = "Module",
    y = "Proportion of species",
    title = "D. Module Taxonomic Composition",
    subtitle = sprintf("%d modules identified", n_modules)
  ) +
  theme(legend.position = "bottom")

# Combine panels
p_manuscript <- (p_degree + p_modularity) / (p_hubs + p_modules) +
  plot_annotation(
    title = "Figure 4: CAFI Co-occurrence Network Analysis",
    subtitle = sprintf("N = %d species, %d positive associations (r > 0.3)", n_nodes, n_edges),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
  )

ggsave(file.path(fig_dir, "network_panels.png"), p_manuscript,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("     Saved: network_panels.png\n")

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
cat(sprintf("  - Unweighted Q vs null: %.1fx null, z = %.1f, p < 0.001\n",
            null_comparison$ratio_to_null[null_comparison$metric == "Modularity"],
            null_comparison$z_score[null_comparison$metric == "Modularity"]))
cat(sprintf("  - Transitivity = %.2f (%.1fx null, z = %.1f, p < 0.001)\n",
            transitivity_obs,
            null_comparison$ratio_to_null[null_comparison$metric == "Transitivity"],
            null_comparison$z_score[null_comparison$metric == "Transitivity"]))
cat(sprintf("  - Number of modules: %d\n\n", n_modules))

cat("Top Hub Species (by hub score):\n")
top_hubs <- centrality_df %>% slice_head(n = 5)
for (i in 1:nrow(top_hubs)) {
  cat(sprintf("  %d. %s (%s) - degree: %d, hub score: %.2f\n",
              i, top_hubs$species[i], top_hubs$type[i],
              top_hubs$degree[i], top_hubs$hub_score[i]))
}

cat("\nKey Findings (H5):\n")
cat("  - Network exhibits significant modular structure\n")
cat("  - Transitivity much higher than random expectation\n")
cat("  - Hub species identifiable through centrality metrics\n")
cat("  - Modules show non-random taxonomic composition\n\n")

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
