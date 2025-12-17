#!/usr/bin/env Rscript
# ============================================================================
# 04_fig4_composition_networks.R - FIGURE 4: Compositional Changes + Networks
# ============================================================================
#
# PURPOSE: Understand how CAFI species composition changes across landscape
# gradients and visualize co-occurrence network structure.
#
# HYPOTHESES:
#   H1: Species composition shifts strongly with coral size and proximity
#   H2: Sites may differ in composition but not strongly (shared habitat)
#   H3: Co-occurrence networks reveal nonrandom species associations
#       (e.g., Trapezia clusters with certain cryptofauna; corallivores
#       occur without defenders)
#
# STATISTICAL APPROACHES:
#   - Ordination: NMDS or PCoA using Hellinger-transformed abundance
#   - Tests: PERMANOVA (adonis2) to test effects of size, proximity, site
#   - Incidence plots: frequency of species presence across size bins
#   - Co-occurrence network: probabilistic co-occurrence or correlations
#   - Network modularity: community detection via Louvain algorithm
#
# Author: CAFI Analysis Pipeline
# Date: 2025-12-03
# ============================================================================

cat("\n========================================\n")
cat("FIGURE 4: Composition & Networks\n")
cat("========================================\n\n")

# Load setup and data
source(here::here("scripts/publication/00_setup.R"))

# Load pre-processed data
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")
community_matrix <- load_object("community_matrix.rds")

# Output directory
FIG_DIR <- FIGURE_DIRS$fig4
cat("Outputs will be saved to:", FIG_DIR, "\n\n")

# ============================================================================
# PREPARE COMMUNITY DATA
# ============================================================================

cat("Preparing community data...\n")

# Filter to corals with CAFI and complete data
coral_ids_with_data <- coral_master %>%
  filter(total_cafi > 0, !is.na(volume), !is.na(mean_neighbor_dist)) %>%
  pull(coral_id)

# Subset community matrix
comm_subset <- community_matrix[rownames(community_matrix) %in% coral_ids_with_data, ]

# Remove empty species (columns)
comm_subset <- comm_subset[, colSums(comm_subset) > 0]

cat("  - Corals with CAFI:", nrow(comm_subset), "\n")
cat("  - OTUs present:", ncol(comm_subset), "\n")

# Hellinger transformation (appropriate for abundance data with zeros)
comm_hell <- decostand(comm_subset, method = "hellinger")

# Prepare environmental data for ordination
env_data <- coral_master %>%
  filter(coral_id %in% rownames(comm_subset)) %>%
  mutate(
    log_volume = log10(volume),
    proximity_m = mean_neighbor_dist / 100,
    size_class = cut(volume,
                     breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE)
  ) %>%
  arrange(match(coral_id, rownames(comm_subset)))

cat("  - Environmental data matched\n\n")

# ============================================================================
# ANALYSIS 1: NMDS ORDINATION
# ============================================================================

cat("1. Running NMDS ordination...\n")

# Calculate Bray-Curtis dissimilarity
dist_bc <- vegdist(comm_hell, method = "bray")

# NMDS with k=2
set.seed(42)
nmds <- metaMDS(dist_bc, k = 2, trymax = 100, trace = 0)

cat("  - NMDS stress:", round(nmds$stress, 3), "\n")
if (nmds$stress < 0.1) {
  cat("  - Excellent representation\n")
} else if (nmds$stress < 0.2) {
  cat("  - Good representation\n")
} else {
  cat("  - Moderate representation - interpret cautiously\n")
}

# Extract NMDS scores
nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("coral_id") %>%
  left_join(env_data, by = "coral_id")

# ============================================================================
# FIGURE 4A: NMDS BY SITE
# ============================================================================

cat("\n2. Creating Figure 4A: NMDS by site...\n")

# Calculate site centroids
site_centroids <- nmds_scores %>%
  group_by(site) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  )

p_4a <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  stat_ellipse(aes(color = site), level = 0.95, linetype = "dashed") +
  geom_point(data = site_centroids, aes(fill = site),
             shape = 23, size = 5, color = "black") +
  scale_color_site() +
  scale_fill_site() +
  labs(
    title = "A. CAFI Community Composition by Site",
    subtitle = sprintf("NMDS (stress = %.3f) | Ellipses show 95%% confidence", nmds$stress),
    x = "NMDS1",
    y = "NMDS2",
    color = "Site", fill = "Site"
  ) +
  coord_fixed() +
  theme_publication()

save_pub_fig(p_4a, "fig4a_nmds_by_site.png", dir = FIG_DIR)

# ============================================================================
# FIGURE 4B: NMDS BY SIZE CLASS
# ============================================================================

cat("3. Creating Figure 4B: NMDS by size class...\n")

size_centroids <- nmds_scores %>%
  group_by(size_class) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  )

p_4b <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = size_class), alpha = 0.6, size = 3) +
  stat_ellipse(aes(color = size_class), level = 0.95, linetype = "dashed") +
  geom_point(data = size_centroids, aes(fill = size_class),
             shape = 23, size = 5, color = "black") +
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  scale_fill_viridis_d(option = "plasma", end = 0.8) +
  labs(
    title = "B. CAFI Community Composition by Coral Size",
    subtitle = "Tertiles of coral volume | Ellipses show 95% confidence",
    x = "NMDS1",
    y = "NMDS2",
    color = "Size Class", fill = "Size Class"
  ) +
  coord_fixed() +
  theme_publication()

save_pub_fig(p_4b, "fig4b_nmds_by_size.png", dir = FIG_DIR)

# ============================================================================
# ANALYSIS 2: PERMANOVA TESTS
# ============================================================================

cat("\n4. Running PERMANOVA tests...\n")

# Test site effect
permanova_site <- adonis2(
  dist_bc ~ site,
  data = env_data,
  permutations = 999
)

cat("  SITE effect:\n")
cat(sprintf("    R2 = %.3f, F = %.2f, p = %.3f\n",
            permanova_site$R2[1],
            permanova_site$F[1],
            permanova_site$`Pr(>F)`[1]))

# Test size effect
permanova_size <- adonis2(
  dist_bc ~ log_volume,
  data = env_data,
  permutations = 999
)

cat("  SIZE effect:\n")
cat(sprintf("    R2 = %.3f, F = %.2f, p = %.3f\n",
            permanova_size$R2[1],
            permanova_size$F[1],
            permanova_size$`Pr(>F)`[1]))

# Test proximity effect
permanova_prox <- adonis2(
  dist_bc ~ proximity_m,
  data = env_data,
  permutations = 999
)

cat("  PROXIMITY effect:\n")
cat(sprintf("    R2 = %.3f, F = %.2f, p = %.3f\n",
            permanova_prox$R2[1],
            permanova_prox$F[1],
            permanova_prox$`Pr(>F)`[1]))

# Full model
permanova_full <- adonis2(
  dist_bc ~ log_volume + proximity_m + site,
  data = env_data,
  permutations = 999
)

cat("\n  FULL MODEL (size + proximity + site):\n")
print(permanova_full)

# ============================================================================
# ANALYSIS 3: INCIDENCE ACROSS SIZE CLASSES
# ============================================================================

cat("\n5. Calculating species incidence across size classes...\n")

# Get top 15 most common OTUs
top_otus <- cafi_clean %>%
  count(otu, sort = TRUE) %>%
  slice_head(n = 15) %>%
  pull(otu)

# Calculate incidence (presence frequency) by size class
incidence_data <- cafi_clean %>%
  filter(coral_id %in% rownames(comm_subset)) %>%
  left_join(env_data %>% select(coral_id, size_class), by = "coral_id") %>%
  filter(otu %in% top_otus) %>%
  group_by(size_class, otu) %>%
  summarise(n_corals_present = n_distinct(coral_id), .groups = "drop") %>%
  left_join(
    env_data %>% count(size_class, name = "n_corals_total"),
    by = "size_class"
  ) %>%
  mutate(incidence = n_corals_present / n_corals_total)

# ============================================================================
# FIGURE 4C: INCIDENCE HEATMAP BY SIZE
# ============================================================================

cat("6. Creating Figure 4C: Incidence heatmap...\n")

p_4c <- incidence_data %>%
  ggplot(aes(x = size_class, y = otu, fill = incidence)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", incidence * 100)),
            size = 3, color = ifelse(incidence_data$incidence > 0.5, "white", "black")) +
  scale_fill_viridis_c(option = "plasma", name = "Incidence\n(% corals)",
                       labels = scales::percent) +
  labs(
    title = "C. Species Incidence Across Coral Size Classes",
    subtitle = "Proportion of corals hosting each OTU | Top 15 most common",
    x = "Coral Size Class",
    y = "OTU (morphotype)"
  ) +
  theme_publication() +
  theme(axis.text.y = element_text(size = 9))

save_pub_fig(p_4c, "fig4c_incidence_heatmap.png", dir = FIG_DIR,
             width = 10, height = 10)

# ============================================================================
# ANALYSIS 4: CO-OCCURRENCE NETWORK
# ============================================================================

cat("\n7. Building co-occurrence network...\n")

# Filter to species occurring in at least 5% of corals
min_occurrence <- ceiling(nrow(comm_subset) * 0.05)
species_keep <- colSums(comm_subset > 0) >= min_occurrence
comm_network <- comm_subset[, species_keep]

cat("  - Species in network:", ncol(comm_network), "\n")

# Calculate Spearman correlations
cor_matrix <- cor(comm_network, method = "spearman", use = "pairwise.complete.obs")

# Build edge list (keep only moderate-strong correlations)
cor_threshold <- 0.3
edge_list <- which(upper.tri(cor_matrix) & abs(cor_matrix) > cor_threshold, arr.ind = TRUE)

if (nrow(edge_list) > 0) {
  edges <- data.frame(
    sp1 = colnames(cor_matrix)[edge_list[, 1]],
    sp2 = colnames(cor_matrix)[edge_list[, 2]],
    correlation = cor_matrix[edge_list],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      association = ifelse(correlation > 0, "positive", "negative"),
      strength = abs(correlation)
    )

  cat("  - Significant associations:", nrow(edges), "\n")
  cat("  - Positive:", sum(edges$association == "positive"), "\n")
  cat("  - Negative:", sum(edges$association == "negative"), "\n")

  # Build igraph network
  g <- graph_from_data_frame(edges[, c("sp1", "sp2", "strength", "association")],
                             directed = FALSE)

  # Add node attributes
  V(g)$type <- sapply(V(g)$name, function(sp) {
    type_val <- cafi_clean %>%
      filter(otu == sp) %>%
      pull(type) %>%
      unique() %>%
      .[1]
    if (is.na(type_val)) "unknown" else type_val
  })

  V(g)$functional_group <- sapply(V(g)$name, function(sp) {
    fg <- cafi_clean %>%
      filter(otu == sp) %>%
      pull(functional_group) %>%
      unique() %>%
      .[1]
    if (is.na(fg)) "Other" else fg
  })

  V(g)$abundance <- sapply(V(g)$name, function(sp) {
    sum(comm_network[, sp], na.rm = TRUE)
  })

  # Network metrics
  V(g)$degree <- degree(g)
  V(g)$betweenness <- betweenness(g)

  # Community detection
  if (vcount(g) > 3) {
    communities <- cluster_louvain(g)
    V(g)$module <- membership(communities)
    modularity_score <- modularity(communities)
    n_modules <- length(unique(membership(communities)))
    cat("  - Modularity:", round(modularity_score, 3), "\n")
    cat("  - Number of modules:", n_modules, "\n")
  }

  # Network summary
  network_summary <- data.frame(
    n_species = vcount(g),
    n_edges = ecount(g),
    density = edge_density(g),
    transitivity = transitivity(g),
    modularity = modularity_score,
    n_modules = n_modules
  )

  save_table(network_summary, "fig4_network_summary.csv")

  # ============================================================================
  # FIGURE 4D: CO-OCCURRENCE NETWORK
  # ============================================================================

  cat("\n8. Creating Figure 4D: Co-occurrence network...\n")

  # Set up layout
  set.seed(123)
  layout <- layout_with_fr(g)

  # Colors by functional group
  node_colors <- FUNCTIONAL_COLORS[V(g)$functional_group]
  node_colors[is.na(node_colors)] <- "#999999"

  # Edge colors by association type
  edge_colors <- ifelse(E(g)$association == "positive", "#0072B2", "#D55E00")

  # Node sizes by degree
  node_sizes <- sqrt(V(g)$degree) * 4 + 3

  # Create network visualization
  png(file.path(FIG_DIR, "fig4d_cooccurrence_network.png"),
      width = 12, height = 10, units = "in", res = 300, bg = "white")

  par(mar = c(1, 1, 3, 1))
  plot(g, layout = layout,
       vertex.size = node_sizes,
       vertex.color = node_colors,
       vertex.label = ifelse(V(g)$degree >= quantile(V(g)$degree, 0.75),
                             V(g)$name, ""),
       vertex.label.cex = 0.6,
       vertex.label.color = "black",
       vertex.frame.color = "gray30",
       edge.color = edge_colors,
       edge.width = E(g)$strength * 2,
       main = "D. CAFI Species Co-occurrence Network")

  # Add legend
  legend("topleft",
         legend = names(FUNCTIONAL_COLORS),
         col = FUNCTIONAL_COLORS,
         pch = 19,
         title = "Functional Group",
         bg = "white",
         cex = 0.8)

  legend("topright",
         legend = c("Positive", "Negative"),
         col = c("#0072B2", "#D55E00"),
         lty = 1,
         lwd = 2,
         title = "Association",
         bg = "white",
         cex = 0.8)

  dev.off()

  cat("  - Network figure saved\n")

} else {
  cat("  - Insufficient co-occurrence patterns for network\n")
  network_summary <- NULL
}

# ============================================================================
# COMBINED FIGURE 4 PANEL
# ============================================================================

cat("\n9. Creating combined Figure 4 panel...\n")

# For combined figure, we'll use the ggplot versions
# Network is saved separately due to igraph plotting

fig4_partial <- (p_4a + p_4b) /
  p_4c +
  plot_annotation(
    title = "Figure 4. CAFI Species Composition and Network Structure",
    subtitle = "NMDS ordination and incidence patterns | Network saved separately",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30")
    )
  ) +
  plot_layout(heights = c(1, 1.2))

save_pub_fig(fig4_partial, "fig4_composition_partial.png",
             dir = FIG_DIR, width = 14, height = 14)

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n10. Saving results...\n")

# Save PERMANOVA results
permanova_results <- bind_rows(
  tibble(predictor = "Site", R2 = permanova_site$R2[1],
         F_stat = permanova_site$F[1], p_value = permanova_site$`Pr(>F)`[1]),
  tibble(predictor = "Coral Size (log)", R2 = permanova_size$R2[1],
         F_stat = permanova_size$F[1], p_value = permanova_size$`Pr(>F)`[1]),
  tibble(predictor = "Proximity", R2 = permanova_prox$R2[1],
         F_stat = permanova_prox$F[1], p_value = permanova_prox$`Pr(>F)`[1])
)

save_table(permanova_results, "fig4_permanova_results.csv")
save_table(incidence_data, "fig4_incidence_by_size.csv")
save_table(edges, "fig4_cooccurrence_edges.csv")

# Save NMDS scores
save_table(nmds_scores, "fig4_nmds_scores.csv")

# Save objects
save_object(list(
  nmds = nmds,
  permanova_full = permanova_full,
  network = if (exists("g")) g else NULL
), "fig4_analysis_objects.rds")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("Figure 4 Analysis Summary\n")
cat("========================================\n\n")

cat("ORDINATION:\n")
cat(sprintf("  - NMDS stress: %.3f\n", nmds$stress))
cat(sprintf("  - Corals analyzed: %d\n", nrow(comm_subset)))
cat(sprintf("  - OTUs included: %d\n", ncol(comm_subset)))

cat("\nPERMANOVA RESULTS:\n")
for (i in 1:nrow(permanova_results)) {
  sig <- if (permanova_results$p_value[i] < 0.05) "*" else ""
  cat(sprintf("  - %s: R2 = %.3f, p = %.3f%s\n",
              permanova_results$predictor[i],
              permanova_results$R2[i],
              permanova_results$p_value[i], sig))
}

if (!is.null(network_summary)) {
  cat("\nNETWORK STRUCTURE:\n")
  cat(sprintf("  - Species in network: %d\n", network_summary$n_species))
  cat(sprintf("  - Associations: %d\n", network_summary$n_edges))
  cat(sprintf("  - Modularity: %.3f\n", network_summary$modularity))
  cat(sprintf("  - Modules detected: %d\n", network_summary$n_modules))
}

cat("\nFigures saved to:", FIG_DIR, "\n")
cat("Tables saved to:", TABLES_DIR, "\n\n")
cat("Figure 4 analysis complete!\n\n")
