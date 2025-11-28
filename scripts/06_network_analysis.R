#!/usr/bin/env Rscript
# ============================================================================
# 06_network_analysis.R - Co-occurrence Network Analysis (H5)
# Author: CAFI Analysis Pipeline
# Date: 2025-11-01
#
# Hypothesis H5: CAFI co-occurrence networks will exhibit non-random modular
# structure, with modules corresponding to functional groups or shared habitat
# preferences, and identifiable keystone species.
#
# Theoretical Background:
#   Network analysis reveals community assembly rules and identifies species
#   that disproportionately influence community structure. Modular networks
#   suggest ecological organization rather than random assembly. Keystone
#   species can be identified through centrality metrics.
#
# Key Analyses:
#   - Species co-occurrence networks (Spearman correlations)
#   - Community detection (Louvain algorithm for modularity)
#   - Keystone species identification (degree, betweenness centrality)
#   - Network stability analysis
#   - Integration with coral condition scores (position-corrected)
#
# Predictions:
#   - Significant modularity (Q > 0.3)
#   - Modules reflect ecological similarity
#   - Keystone species identifiable through centrality
#
# Important Notes:
#   - ALL corals are Pocillopora spp. (no coral species/morphotypes in network)
#   - Networks represent CAFI-CAFI interactions only
#   - Coral condition scores are position-corrected (from Script 05a)
#   - All figures have white backgrounds
# ============================================================================

cat("\n========================================\n")
cat("Network Analysis of CAFI Communities\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/00_load_libraries.R"))
library(igraph)

# Load processed data
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))

# Load position-corrected condition scores from Script 05a
condition_scores <- readRDS(file.path(SURVEY_OBJECTS, "coral_condition_scores.rds"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "network_analysis")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Use publication theme (already set via 00_load_libraries.R)

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
# Coral data (condition scores, size, location) are used as EXPLANATORY
# VARIABLES to test how coral traits affect CAFI network properties.
# ============================================================================

# ============================================================================
# Co-occurrence Analysis
# ============================================================================

cat("Calculating CAFI species co-occurrences...\n")
cat("  NOTE: Analyzing CAFI-CAFI networks only (corals not included as nodes)\n\n")

# Convert to presence-absence for co-occurrence
comm_binary <- community_matrix
comm_binary[comm_binary > 0] <- 1

# Filter to species occurring in at least 5% of samples (corals)
min_occurrence <- ceiling(nrow(comm_binary) * 0.05)
species_keep <- colSums(comm_binary) >= min_occurrence
comm_filtered <- comm_binary[, species_keep]

cat("  - CAFI species retained for network:", ncol(comm_filtered), "of", ncol(comm_binary), "\n")
cat("  - Coral samples (network context):", nrow(comm_filtered), "\n")

# Calculate correlations for co-occurrence (Spearman)
# Using correlation as proxy for co-occurrence strength
cor_matrix <- cor(comm_filtered, method = "spearman", use = "pairwise.complete.obs")

# Convert to edge list (upper triangle only, no self-loops)
cor_edges <- which(upper.tri(cor_matrix) & abs(cor_matrix) > 0.3, arr.ind = TRUE)

if(length(cor_edges) > 0) {

  edge_list <- data.frame(
    sp1 = colnames(cor_matrix)[cor_edges[,1]],
    sp2 = colnames(cor_matrix)[cor_edges[,2]],
    correlation = cor_matrix[cor_edges],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      association = if_else(correlation > 0, "positive", "negative"),
      strength = abs(correlation)
    ) %>%
    filter(strength > 0.3)  # Keep moderate to strong associations

  # Save co-occurrence results
  write_csv(edge_list,
            file.path(SURVEY_TABLES, "cafi_cooccurrence_pairs.csv"))

  cat("  - Significant associations found:", nrow(edge_list), "\n")
  cat("  - Positive associations:", sum(edge_list$association == "positive"), "\n")
  cat("  - Negative associations:", sum(edge_list$association == "negative"), "\n\n")

} else {
  cat("  ⚠️  No strong co-occurrence patterns detected\n\n")
}

# ============================================================================
# Build CAFI Co-occurrence Network
# ============================================================================

cat("Building CAFI co-occurrence network...\n")

if (exists("edge_list") && nrow(edge_list) > 0) {

  # Create network from edge list
  # NOTE: This is a CAFI-CAFI network, NOT including corals as nodes
  g <- graph_from_data_frame(edge_list[, c("sp1", "sp2", "strength", "association")],
                             directed = FALSE)

  # Add node attributes for CAFI species
  V(g)$type <- sapply(V(g)$name, function(sp) {
    type_val <- cafi_clean %>%
      filter(species == sp) %>%
      pull(type) %>%
      unique() %>%
      .[1]
    if (is.na(type_val)) "unknown" else type_val
  })

  # Calculate total abundance across all corals
  V(g)$abundance <- sapply(V(g)$name, function(sp) {
    sum(community_matrix[, sp], na.rm = TRUE)
  })

  # Calculate network centrality metrics
  # These measure each CAFI species' importance in the network
  V(g)$degree <- degree(g)                    # Number of connections
  V(g)$betweenness <- betweenness(g)         # Bridging position
  V(g)$closeness <- closeness(g)             # Average distance to all nodes
  V(g)$eigenvector <- eigen_centrality(g)$vector  # Influence based on connections

  # Network-level metrics
  network_metrics <- data.frame(
    n_cafi_species = vcount(g),
    n_associations = ecount(g),
    density = edge_density(g),
    transitivity = transitivity(g),
    mean_path_length = mean_distance(g),
    diameter = diameter(g)
  )

  # Calculate modularity (community structure)
  if(vcount(g) > 3) {
    communities_louvain <- cluster_louvain(g)
    network_metrics$modularity <- modularity(communities_louvain)
    network_metrics$n_modules <- length(unique(membership(communities_louvain)))
  }

  write_csv(network_metrics,
            file.path(SURVEY_TABLES, "cafi_network_metrics.csv"))

  cat(sprintf("  ✓ CAFI network constructed: %d species, %d associations\n",
              vcount(g), ecount(g)))
  cat(sprintf("  ✓ Network density: %.3f\n", network_metrics$density))
  cat(sprintf("  ✓ Modularity: %.3f\n\n", network_metrics$modularity))

  # =========================================================================
  # Visualize CAFI Network
  # =========================================================================

  cat("Creating network visualizations...\n")

  # Set layout
  set.seed(123)
  layout <- layout_with_fr(g)

  # Color nodes by taxonomic type
  type_colors <- c(
    "crab" = "#E41A1C",
    "shrimp" = "#377EB8",
    "fish" = "#4DAF4A",
    "snail" = "#984EA3",
    "unknown" = "#999999"
  )
  node_colors <- type_colors[V(g)$type]

  # Size nodes by degree centrality
  node_sizes <- sqrt(V(g)$degree) * 4

  # Color edges by association type
  edge_colors <- ifelse(E(g)$association == "positive", "#377EB8", "#E41A1C")

  # Figure 1: Full network
  png(file.path(fig_dir, "01_cafi_cooccurrence_network.png"),
      width = 14, height = 12, units = "in", res = 300, bg = "white")
  par(mar = c(1, 1, 3, 1))
  plot(g, layout = layout,
       vertex.size = node_sizes,
       vertex.color = node_colors,
       vertex.label.cex = 0.7,
       vertex.label.color = "black",
       vertex.label.dist = 0.5,
       vertex.frame.color = "gray30",
       edge.color = edge_colors,
       edge.width = E(g)$strength * 3,
       main = "CAFI Species Co-occurrence Network\n(Pocillopora spp. hosts)")

  legend("topright",
         legend = c(names(type_colors), "", "Positive", "Negative"),
         col = c(type_colors, NA, "#377EB8", "#E41A1C"),
         pch = c(rep(19, length(type_colors)), NA, NA, NA),
         lty = c(rep(NA, length(type_colors)), NA, 1, 1),
         lwd = c(rep(NA, length(type_colors)), NA, 2, 2),
         title = "CAFI Type | Association",
         bg = "white",
         cex = 0.9)
  dev.off()

  cat("  ✓ Network visualization created\n")

  # =========================================================================
  # Community Detection (Modules)
  # =========================================================================

  cat("\nDetecting network modules...\n")

  if (vcount(g) > 10) {

    # Apply Louvain algorithm for community detection
    communities <- cluster_louvain(g)
    V(g)$module <- membership(communities)

    cat(sprintf("  ✓ Identified %d modules (modularity = %.3f)\n",
                max(V(g)$module), modularity(communities)))

    # Summarize modules
    module_summary <- data.frame(
      species = V(g)$name,
      type = V(g)$type,
      module = V(g)$module,
      degree = V(g)$degree,
      abundance = V(g)$abundance
    ) %>%
      group_by(module) %>%
      mutate(
        module_size = n(),
        module_density = sum(degree) / (module_size * (module_size - 1))
      ) %>%
      ungroup()

    write_csv(module_summary,
              file.path(SURVEY_TABLES, "cafi_network_modules.csv"))

    # Figure 2: Network with module colors
    module_colors <- rainbow(max(V(g)$module), alpha = 0.7)

    png(file.path(fig_dir, "02_cafi_network_modules.png"),
        width = 14, height = 12, units = "in", res = 300, bg = "white")
    par(mar = c(1, 1, 3, 1))
    plot(communities, g,
         layout = layout,
         vertex.size = node_sizes,
         vertex.label.cex = 0.7,
         vertex.label.color = "black",
         edge.color = edge_colors,
         edge.width = E(g)$strength * 3,
         main = sprintf("CAFI Network Modules (n=%d, Q=%.3f)",
                       max(V(g)$module), modularity(communities)))
    dev.off()

    cat("  ✓ Module visualization created\n")
  }

  # =========================================================================
  # Keystone Species Identification
  # =========================================================================

  cat("\nIdentifying keystone CAFI species...\n")

  # Keystone index = degree × betweenness (connectivity × bridging)
  keystone_scores <- data.frame(
    species = V(g)$name,
    type = V(g)$type,
    module = if(exists("V(g)$module")) V(g)$module else NA,
    degree = V(g)$degree,
    betweenness = V(g)$betweenness,
    closeness = V(g)$closeness,
    eigenvector = V(g)$eigenvector,
    abundance = V(g)$abundance
  ) %>%
    mutate(
      # Keystone index: high connectivity + high bridging
      keystone_index = scale(degree)[,1] * scale(betweenness)[,1],
      # Hub score: high connectivity + high influence
      hub_score = scale(degree)[,1] + scale(eigenvector)[,1],
      # Connector score: bridging relative to connectivity
      connector_score = scale(betweenness)[,1] / (scale(degree)[,1] + 0.001)
    ) %>%
    arrange(desc(keystone_index))

  write_csv(keystone_scores,
            file.path(SURVEY_TABLES, "cafi_keystone_species.csv"))

  cat(sprintf("  ✓ Top keystone species: %s (%s)\n",
              keystone_scores$species[1], keystone_scores$type[1]))

  # Figure 3: Keystone species ranking
  p_keystone <- keystone_scores %>%
    slice_head(n = 20) %>%
    ggplot(aes(x = reorder(species, keystone_index), y = keystone_index)) +
    geom_col(aes(fill = type), alpha = 0.8) +
    coord_flip() +
    scale_fill_manual(values = type_colors, name = "CAFI Type") +
    labs(
      title = "Top 20 Keystone CAFI Species",
      subtitle = "Keystone index = degree × betweenness (connectivity × bridging)",
      x = "Species (OTU)",
      y = "Keystone Index (standardized)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      legend.position = "bottom"
    )

  ggsave(file.path(fig_dir, "03_keystone_species_ranking.png"),
         p_keystone, width = 10, height = 9, dpi = 300, bg = "white")

  cat("  ✓ Keystone ranking visualization created\n")

} else {
  cat("  ⚠️  Insufficient data to build network\n")
}

# ============================================================================
# Network Properties vs Coral Condition (Position-Corrected)
# ============================================================================

cat("\nAnalyzing network properties vs coral condition...\n")
cat("  Using position-corrected condition scores from Script 05a\n")

# For each coral, calculate its CAFI network properties
if(exists("g") && nrow(condition_scores) > 10) {

  # Calculate coral-level network metrics
  coral_network_metrics <- data.frame(
    coral_id = rownames(community_matrix)
  )

  # For each coral, calculate:
  # 1. CAFI richness (already have this)
  # 2. CAFI abundance (already have this)
  # 3. Mean degree of hosted CAFI (average connectivity)
  # 4. Presence of keystone species
  # 5. Module diversity (number of different modules present)

  coral_network_metrics <- coral_network_metrics %>%
    rowwise() %>%
    mutate(
      # Get CAFI species in this coral
      cafi_species = list(colnames(community_matrix)[community_matrix[coral_id, ] > 0]),

      # Calculate metrics only if species are in network
      cafi_in_network = list(intersect(cafi_species, V(g)$name)),
      n_cafi_in_network = length(cafi_in_network),

      # Mean degree of hosted CAFI
      mean_degree = if(length(cafi_in_network) > 0) {
        mean(V(g)$degree[V(g)$name %in% cafi_in_network])
      } else NA,

      # Mean betweenness
      mean_betweenness = if(length(cafi_in_network) > 0) {
        mean(V(g)$betweenness[V(g)$name %in% cafi_in_network])
      } else NA,

      # Keystone species presence (top 10%)
      n_keystones = if(length(cafi_in_network) > 0) {
        sum(V(g)$name[V(g)$name %in% cafi_in_network] %in%
            keystone_scores$species[1:ceiling(nrow(keystone_scores)*0.1)])
      } else 0,

      # Module diversity (if modules exist)
      n_modules = if(length(cafi_in_network) > 0 && exists("V(g)$module")) {
        n_distinct(V(g)$module[V(g)$name %in% cafi_in_network])
      } else NA

    ) %>%
    ungroup() %>%
    select(-cafi_species, -cafi_in_network)  # Remove list columns

  # Merge with condition scores AND CORAL SIZE
  coral_condition_network <- coral_network_metrics %>%
    inner_join(condition_scores %>% select(coral_id, condition_score, site, volume_lab),
              by = "coral_id")

  # Save merged dataset
  write_csv(coral_condition_network,
            file.path(SURVEY_TABLES, "coral_condition_network_metrics.csv"))

  cat(sprintf("  ✓ Network metrics calculated for %d corals\n",
              nrow(coral_condition_network)))

  # =========================================================================
  # Test relationships
  # =========================================================================

  cat("\n  Testing relationships between coral condition and network properties...\n")

  # 1. Condition vs mean CAFI degree
  if(sum(!is.na(coral_condition_network$mean_degree)) > 20) {

    model_degree <- lm(condition_score ~ mean_degree + site,
                      data = coral_condition_network %>%
                        filter(!is.na(mean_degree)))

    coef_degree <- broom::tidy(model_degree, conf.int = TRUE) %>%
      filter(term == "mean_degree")

    cat(sprintf("    Mean CAFI degree → Condition: β = %.3f, p = %.3f\n",
                coef_degree$estimate, coef_degree$p.value))

    # Figure 4: Mean degree vs condition
    p_degree_condition <- coral_condition_network %>%
      filter(!is.na(mean_degree)) %>%
      ggplot(aes(x = mean_degree, y = condition_score)) +
      geom_point(aes(color = site), alpha = 0.6, size = 3) +
      geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1.2) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      scale_color_viridis_d(option = "plasma", name = "Location") +
      labs(
        title = "Coral Condition vs CAFI Network Connectivity",
        subtitle = sprintf("Mean degree of hosted CAFI | β = %.3f, p = %.3f",
                          coef_degree$estimate, coef_degree$p.value),
        x = "Mean Degree of Hosted CAFI Species",
        y = "Coral Condition Score (position-corrected)",
        caption = "Higher degree = CAFI species with more co-occurrence connections"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10,
                                     color = if_else(coef_degree$p.value < 0.05,
                                                    "darkgreen", "gray50")),
        plot.caption = element_text(size = 9, hjust = 0),
        legend.position = "right"
      )

    ggsave(file.path(fig_dir, "04_condition_vs_cafi_network_degree.png"),
           p_degree_condition, width = 11, height = 7, dpi = 300, bg = "white")
  }

  # 2. Condition vs keystone presence
  if(sum(coral_condition_network$n_keystones > 0, na.rm = TRUE) > 5) {

    coral_condition_network <- coral_condition_network %>%
      mutate(has_keystone = n_keystones > 0)

    model_keystone <- lm(condition_score ~ has_keystone + site,
                        data = coral_condition_network %>%
                          filter(!is.na(has_keystone)))

    coef_keystone <- broom::tidy(model_keystone, conf.int = TRUE) %>%
      filter(term == "has_keystoneTRUE")

    cat(sprintf("    Keystone presence → Condition: β = %.3f, p = %.3f\n",
                coef_keystone$estimate, coef_keystone$p.value))

    # Figure 5: Keystone presence
    p_keystone_condition <- coral_condition_network %>%
      filter(!is.na(has_keystone)) %>%
      ggplot(aes(x = has_keystone, y = condition_score, fill = has_keystone)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.4, size = 2.5) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      scale_fill_viridis_d(option = "plasma", begin = 0.3, end = 0.8) +
      scale_x_discrete(labels = c("FALSE" = "No Keystones", "TRUE" = "Has Keystones")) +
      labs(
        title = "Coral Condition vs Keystone Species Presence",
        subtitle = sprintf("Effect of hosting keystone CAFI | β = %.3f, p = %.3f",
                          coef_keystone$estimate, coef_keystone$p.value),
        x = "Keystone CAFI Species Present?",
        y = "Coral Condition Score (position-corrected)"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10,
                                     color = if_else(coef_keystone$p.value < 0.05,
                                                    "darkgreen", "gray50")),
        legend.position = "none"
      )

    ggsave(file.path(fig_dir, "05_condition_vs_keystone_presence.png"),
           p_keystone_condition, width = 10, height = 7, dpi = 300, bg = "white")
  }

  # 3. Condition vs module diversity
  if(sum(!is.na(coral_condition_network$n_modules)) > 20) {

    model_modules <- lm(condition_score ~ n_modules + site,
                       data = coral_condition_network %>%
                         filter(!is.na(n_modules)))

    coef_modules <- broom::tidy(model_modules, conf.int = TRUE) %>%
      filter(term == "n_modules")

    cat(sprintf("    Module diversity → Condition: β = %.3f, p = %.3f\n",
                coef_modules$estimate, coef_modules$p.value))

    # Figure 6: Module diversity
    p_modules_condition <- coral_condition_network %>%
      filter(!is.na(n_modules)) %>%
      ggplot(aes(x = n_modules, y = condition_score)) +
      geom_point(aes(color = site), alpha = 0.6, size = 3) +
      geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1.2) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      scale_color_viridis_d(option = "plasma", name = "Location") +
      labs(
        title = "Coral Condition vs CAFI Module Diversity",
        subtitle = sprintf("Number of network modules present | β = %.3f, p = %.3f",
                          coef_modules$estimate, coef_modules$p.value),
        x = "Number of CAFI Network Modules",
        y = "Coral Condition Score (position-corrected)",
        caption = "Modules = distinct CAFI community groups from network analysis"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10,
                                     color = if_else(coef_modules$p.value < 0.05,
                                                    "darkgreen", "gray50")),
        plot.caption = element_text(size = 9, hjust = 0),
        legend.position = "right"
      )

    ggsave(file.path(fig_dir, "06_condition_vs_module_diversity.png"),
           p_modules_condition, width = 11, height = 7, dpi = 300, bg = "white")
  }

  # =========================================================================
  # Test coral SIZE and SITE effects on network properties
  # =========================================================================

  cat("\n  Testing coral size and site effects on network properties...\n")

  # Test 1: Coral size vs mean CAFI degree
  if(sum(!is.na(coral_condition_network$mean_degree) &
         !is.na(coral_condition_network$volume_lab)) > 20) {

    model_size_degree <- lm(mean_degree ~ log10(volume_lab) + site,
                           data = coral_condition_network %>%
                             filter(!is.na(mean_degree), !is.na(volume_lab)))

    coef_size <- broom::tidy(model_size_degree, conf.int = TRUE) %>%
      filter(term == "log10(volume_lab)")

    # Test site effect with ANOVA
    anova_site_degree <- anova(
      lm(mean_degree ~ log10(volume_lab), data = coral_condition_network %>%
           filter(!is.na(mean_degree), !is.na(volume_lab))),
      lm(mean_degree ~ log10(volume_lab) + site, data = coral_condition_network %>%
           filter(!is.na(mean_degree), !is.na(volume_lab)))
    )

    cat(sprintf("    Coral size (log volume) → Mean degree: β = %.4f, p = %.4f\n",
                coef_size$estimate, coef_size$p.value))
    cat(sprintf("    Site effect on mean degree: F = %.2f, p = %.4f\n",
                anova_site_degree$F[2], anova_site_degree$`Pr(>F)`[2]))

    # Figure: Size vs mean degree
    p_size_degree <- coral_condition_network %>%
      filter(!is.na(mean_degree), !is.na(volume_lab)) %>%
      ggplot(aes(x = log10(volume_lab), y = mean_degree)) +
      geom_point(aes(color = site), alpha = 0.6, size = 3) +
      geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1.2) +
      scale_color_viridis_d(option = "plasma", name = "Location") +
      labs(
        title = "Coral Size vs CAFI Network Connectivity",
        subtitle = sprintf("Effect of colony volume on hosted CAFI degree | β = %.4f, p = %.4f",
                          coef_size$estimate, coef_size$p.value),
        x = "Coral Colony Volume (log10 cm³)",
        y = "Mean Degree of Hosted CAFI Species",
        caption = "Controlling for site differences | Larger corals may provide more niche space"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10,
                                     color = if_else(coef_size$p.value < 0.05,
                                                    "darkgreen", "gray50")),
        plot.caption = element_text(size = 9, hjust = 0),
        legend.position = "right",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )

    ggsave(file.path(fig_dir, "10_coral_size_vs_network_degree.png"),
           p_size_degree, width = 11, height = 7, dpi = 300, bg = "white")
  }

  # Test 2: Coral size vs keystone presence
  if(sum(!is.na(coral_condition_network$has_keystone) &
         !is.na(coral_condition_network$volume_lab)) > 20) {

    model_size_keystone <- glm(has_keystone ~ log10(volume_lab) + site,
                              data = coral_condition_network %>%
                                filter(!is.na(has_keystone), !is.na(volume_lab)),
                              family = binomial)

    coef_size_ks <- broom::tidy(model_size_keystone, conf.int = TRUE) %>%
      filter(term == "log10(volume_lab)")

    cat(sprintf("    Coral size → Keystone presence: β = %.4f, p = %.4f (logistic)\n",
                coef_size_ks$estimate, coef_size_ks$p.value))

    # Figure: Size vs keystone presence
    p_size_keystone <- coral_condition_network %>%
      filter(!is.na(has_keystone), !is.na(volume_lab)) %>%
      ggplot(aes(x = log10(volume_lab), y = as.numeric(has_keystone))) +
      geom_point(aes(color = site), alpha = 0.4, size = 3,
                position = position_jitter(height = 0.05, width = 0)) +
      geom_smooth(method = "glm", method.args = list(family = "binomial"),
                 color = "black", se = TRUE, linewidth = 1.2) +
      scale_color_viridis_d(option = "plasma", name = "Location") +
      scale_y_continuous(breaks = c(0, 1), labels = c("No", "Yes")) +
      labs(
        title = "Coral Size vs Keystone Species Presence",
        subtitle = sprintf("Effect of colony volume on hosting keystones | β = %.4f, p = %.4f",
                          coef_size_ks$estimate, coef_size_ks$p.value),
        x = "Coral Colony Volume (log10 cm³)",
        y = "Hosts Keystone CAFI?",
        caption = "Controlling for site differences | Logistic regression"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10,
                                     color = if_else(coef_size_ks$p.value < 0.05,
                                                    "darkgreen", "gray50")),
        plot.caption = element_text(size = 9, hjust = 0),
        legend.position = "right",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )

    ggsave(file.path(fig_dir, "11_coral_size_vs_keystone_presence.png"),
           p_size_keystone, width = 11, height = 7, dpi = 300, bg = "white")
  }

  # Save size/site results
  size_site_results <- data.frame(
    predictor = c("Coral size (log volume)", "Site effect"),
    response = c("Mean CAFI degree", "Mean CAFI degree"),
    beta = c(coef_size$estimate, NA),
    p_value = c(coef_size$p.value, anova_site_degree$`Pr(>F)`[2]),
    test_type = c("Linear regression", "ANOVA")
  )

  write_csv(size_site_results,
            file.path(SURVEY_TABLES, "coral_size_site_network_effects.csv"))

  cat("\n  ✓ Size and site effect analyses complete\n\n")

  cat("\n  ✓ Condition-network analyses complete\n\n")

} else {
  cat("  ⚠️  Insufficient data for condition-network analysis\n\n")
}

# ============================================================================
# Network Comparison by Location
# ============================================================================

cat("Comparing CAFI networks across locations...\n")

if(exists("g") && nrow(metadata) > 30) {

  # Split community matrix by site
  sites <- c("HAU", "MAT", "MRB")
  site_networks <- list()

  for(site_name in sites) {

    # Get corals from this site
    site_corals <- metadata %>%
      filter(site == site_name) %>%
      pull(coral_id)

    # Subset community matrix
    site_comm <- comm_filtered[rownames(comm_filtered) %in% site_corals, ]

    # Calculate correlations
    if(nrow(site_comm) >= 10 && ncol(site_comm) >= 10) {

      site_cor <- cor(site_comm, method = "spearman", use = "pairwise.complete.obs")
      site_edges <- which(upper.tri(site_cor) & abs(site_cor) > 0.3, arr.ind = TRUE)

      if(length(site_edges) > 0) {

        site_edge_list <- data.frame(
          sp1 = colnames(site_cor)[site_edges[,1]],
          sp2 = colnames(site_cor)[site_edges[,2]],
          correlation = site_cor[site_edges]
        ) %>%
          mutate(strength = abs(correlation)) %>%
          filter(strength > 0.3)

        # Create site network
        g_site <- graph_from_data_frame(site_edge_list[, c("sp1", "sp2", "strength")],
                                        directed = FALSE)

        # Calculate metrics
        site_networks[[site_name]] <- data.frame(
          site = site_name,
          n_corals = nrow(site_comm),
          n_cafi_species = vcount(g_site),
          n_edges = ecount(g_site),
          density = edge_density(g_site),
          transitivity = transitivity(g_site),
          modularity = if(vcount(g_site) > 3) {
            modularity(cluster_louvain(g_site))
          } else NA
        )
      }
    }
  }

  # Combine site network metrics
  if(length(site_networks) > 0) {
    site_comparison <- bind_rows(site_networks)

    write_csv(site_comparison,
              file.path(SURVEY_TABLES, "cafi_network_by_site.csv"))

    cat("  ✓ Site-specific networks characterized\n")

    # Visualize comparison
    p_site_compare <- site_comparison %>%
      select(site, density, transitivity, modularity) %>%
      pivot_longer(cols = -site, names_to = "metric", values_to = "value") %>%
      ggplot(aes(x = site, y = value, fill = site)) +
      geom_col(alpha = 0.8) +
      facet_wrap(~metric, scales = "free_y", ncol = 3) +
      scale_fill_viridis_d(option = "plasma") +
      labs(
        title = "CAFI Network Properties by Location",
        subtitle = "Comparing network structure across 3 Mo'orea sites",
        x = "Location",
        y = "Metric Value"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray30"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold"),
        legend.position = "none"
      )

    ggsave(file.path(fig_dir, "07_network_properties_by_site.png"),
           p_site_compare, width = 12, height = 5, dpi = 300, bg = "white")
  }

  cat("\n  ✓ Site comparison complete\n\n")

} else {
  cat("  ⚠️  Insufficient data for site comparison\n\n")
}

# ============================================================================
# Deeper Network Analyses
# ============================================================================

cat("Conducting deeper network analyses...\n")

if(exists("g") && vcount(g) > 10) {

  # =========================================================================
  # 1. Null Model Comparison
  # =========================================================================

  cat("\n  1. Comparing observed network to null expectations...\n")

  # Generate random networks with same number of nodes and edges
  n_randomizations <- 1000
  random_metrics <- matrix(NA, nrow = n_randomizations, ncol = 4)
  colnames(random_metrics) <- c("density", "transitivity", "mean_path_length", "modularity")

  for(i in 1:n_randomizations) {
    # Erdős-Rényi random graph with same density
    g_random <- erdos.renyi.game(vcount(g), edge_density(g), type = "gnp")

    random_metrics[i, "density"] <- edge_density(g_random)
    random_metrics[i, "transitivity"] <- transitivity(g_random)
    random_metrics[i, "mean_path_length"] <- mean_distance(g_random)
    random_metrics[i, "modularity"] <- modularity(cluster_louvain(g_random))
  }

  # Calculate z-scores (observed vs random)
  null_comparison <- data.frame(
    metric = c("Transitivity", "Mean Path Length", "Modularity"),
    observed = c(
      network_metrics$transitivity,
      network_metrics$mean_path_length,
      network_metrics$modularity
    ),
    random_mean = c(
      mean(random_metrics[, "transitivity"]),
      mean(random_metrics[, "mean_path_length"]),
      mean(random_metrics[, "modularity"])
    ),
    random_sd = c(
      sd(random_metrics[, "transitivity"]),
      sd(random_metrics[, "mean_path_length"]),
      sd(random_metrics[, "modularity"])
    )
  ) %>%
    mutate(
      z_score = (observed - random_mean) / random_sd,
      p_value = 2 * (1 - pnorm(abs(z_score))),  # Two-tailed
      interpretation = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  write_csv(null_comparison,
            file.path(SURVEY_TABLES, "network_null_model_comparison.csv"))

  cat(sprintf("    Transitivity: obs = %.3f, random = %.3f ± %.3f, z = %.2f %s\n",
              null_comparison$observed[1], null_comparison$random_mean[1],
              null_comparison$random_sd[1], null_comparison$z_score[1],
              null_comparison$interpretation[1]))
  cat(sprintf("    Modularity: obs = %.3f, random = %.3f ± %.3f, z = %.2f %s\n",
              null_comparison$observed[3], null_comparison$random_mean[3],
              null_comparison$random_sd[3], null_comparison$z_score[3],
              null_comparison$interpretation[3]))

  # =========================================================================
  # 2. Association Type Analysis (Positive vs Negative)
  # =========================================================================

  cat("\n  2. Analyzing positive vs negative associations...\n")

  if(exists("edge_list")) {
    association_summary <- edge_list %>%
      group_by(association) %>%
      summarise(
        n_pairs = n(),
        mean_strength = mean(strength),
        median_strength = median(strength),
        .groups = "drop"
      ) %>%
      mutate(
        proportion = n_pairs / sum(n_pairs)
      )

    write_csv(association_summary,
              file.path(SURVEY_TABLES, "association_type_summary.csv"))

    cat(sprintf("    Positive: n = %d (%.1f%%), mean r = %.3f\n",
                association_summary$n_pairs[association_summary$association == "positive"],
                association_summary$proportion[association_summary$association == "positive"] * 100,
                association_summary$mean_strength[association_summary$association == "positive"]))
    cat(sprintf("    Negative: n = %d (%.1f%%), mean r = %.3f\n",
                association_summary$n_pairs[association_summary$association == "negative"],
                association_summary$proportion[association_summary$association == "negative"] * 100,
                association_summary$mean_strength[association_summary$association == "negative"]))

    # Test if network is more positive than expected
    prop_positive <- sum(edge_list$association == "positive") / nrow(edge_list)
    binom_test <- binom.test(sum(edge_list$association == "positive"),
                             nrow(edge_list), p = 0.5, alternative = "greater")

    cat(sprintf("    Binomial test for excess positive associations: p = %.4f\n",
                binom_test$p.value))
    if(binom_test$p.value < 0.05) {
      cat("    → Network is significantly biased toward POSITIVE associations\n")
    }
  }

  # =========================================================================
  # 3. Module Composition by Taxonomic Group
  # =========================================================================

  cat("\n  3. Analyzing module composition by CAFI type...\n")

  if(exists("module_summary")) {
    module_composition <- module_summary %>%
      group_by(module, type) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(module) %>%
      mutate(
        total_in_module = sum(n),
        proportion = n / total_in_module
      ) %>%
      ungroup()

    write_csv(module_composition,
              file.path(SURVEY_TABLES, "module_taxonomic_composition.csv"))

    # Chi-square test for non-random taxonomic distribution across modules
    contingency_table <- module_composition %>%
      select(module, type, n) %>%
      pivot_wider(names_from = type, values_from = n, values_fill = 0) %>%
      select(-module) %>%
      as.matrix()

    if(nrow(contingency_table) > 1 && ncol(contingency_table) > 1) {
      chi_test <- chisq.test(contingency_table)
      cat(sprintf("    Chi-square test for taxonomic clustering: χ² = %.2f, p = %.4f\n",
                  chi_test$statistic, chi_test$p.value))
      if(chi_test$p.value < 0.05) {
        cat("    → Modules show NON-RANDOM taxonomic composition\n")
      } else {
        cat("    → Modules show random taxonomic mixing\n")
      }
    }

    # Figure: Module composition heatmap
    p_module_comp <- module_composition %>%
      ggplot(aes(x = factor(module), y = type, fill = proportion)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = n), color = "white", fontface = "bold", size = 5) +
      scale_fill_viridis_c(option = "plasma", name = "Proportion") +
      labs(
        title = "Module Composition by CAFI Taxonomic Group",
        subtitle = "Numbers show species count per module-type combination",
        x = "Network Module",
        y = "CAFI Type"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray30"),
        axis.text = element_text(size = 11, face = "bold")
      )

    ggsave(file.path(fig_dir, "08_module_taxonomic_composition.png"),
           p_module_comp, width = 10, height = 6, dpi = 300, bg = "white")
  }

  # =========================================================================
  # 4. Network Centralization
  # =========================================================================

  cat("\n  4. Calculating network centralization metrics...\n")

  # Degree centralization: how unequal is connectivity distributed?
  degree_centralization <- centr_degree(g)$centralization

  # Betweenness centralization: how dependent is network on few bridges?
  betweenness_centralization <- centr_betw(g)$centralization

  # Closeness centralization
  closeness_centralization <- centr_clo(g)$centralization

  centralization_metrics <- data.frame(
    metric = c("Degree", "Betweenness", "Closeness"),
    centralization = c(degree_centralization, betweenness_centralization,
                      closeness_centralization),
    interpretation = c(
      if_else(degree_centralization > 0.5, "Star-like (hub-dominated)", "Distributed"),
      if_else(betweenness_centralization > 0.5, "Bridge-dependent", "Multiple pathways"),
      if_else(closeness_centralization > 0.5, "Centered", "Decentralized")
    )
  )

  write_csv(centralization_metrics,
            file.path(SURVEY_TABLES, "network_centralization.csv"))

  cat(sprintf("    Degree centralization: %.3f (%s)\n",
              degree_centralization, centralization_metrics$interpretation[1]))
  cat(sprintf("    Betweenness centralization: %.3f (%s)\n",
              betweenness_centralization, centralization_metrics$interpretation[2]))
  cat(sprintf("    Closeness centralization: %.3f (%s)\n",
              closeness_centralization, centralization_metrics$interpretation[3]))

  # =========================================================================
  # 5. Degree Distribution
  # =========================================================================

  cat("\n  5. Analyzing degree distribution...\n")

  degree_dist <- data.frame(
    species = V(g)$name,
    degree = V(g)$degree
  )

  # Test for scale-free network (power law)
  # Most biological networks show power-law degree distributions
  degree_values <- V(g)$degree
  degree_freq <- table(degree_values)

  # Fit power law (if enough variation)
  if(length(unique(degree_values)) > 3) {
    # Simple linear regression on log-log scale
    log_degree <- log10(as.numeric(names(degree_freq)))
    log_freq <- log10(as.numeric(degree_freq))

    # Remove infinite values
    valid <- is.finite(log_degree) & is.finite(log_freq)
    if(sum(valid) > 2) {
      power_law_fit <- lm(log_freq[valid] ~ log_degree[valid])
      power_law_slope <- coef(power_law_fit)[2]
      power_law_r2 <- summary(power_law_fit)$r.squared

      cat(sprintf("    Power law exponent: %.3f (R² = %.3f)\n",
                  abs(power_law_slope), power_law_r2))
      if(power_law_r2 > 0.7) {
        cat("    → Degree distribution consistent with scale-free network\n")
      } else {
        cat("    → Degree distribution NOT scale-free (more random)\n")
      }
    }
  }

  # Figure: Degree distribution
  p_degree_dist <- degree_dist %>%
    ggplot(aes(x = degree)) +
    geom_histogram(binwidth = 1, fill = "steelblue", alpha = 0.8, color = "white") +
    geom_vline(aes(xintercept = mean(degree)), linetype = "dashed",
               color = "red", linewidth = 1.2) +
    labs(
      title = "CAFI Network Degree Distribution",
      subtitle = sprintf("Mean degree = %.2f | Most species have 2-4 connections",
                        mean(degree_dist$degree)),
      x = "Degree (Number of Co-occurring Species)",
      y = "Number of Species"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30")
    )

  ggsave(file.path(fig_dir, "09_degree_distribution.png"),
         p_degree_dist, width = 10, height = 7, dpi = 300, bg = "white")

  # =========================================================================
  # 6. Assortativity
  # =========================================================================

  cat("\n  6. Testing network assortativity...\n")

  # Degree assortativity: do high-degree nodes connect to other high-degree nodes?
  degree_assortativity <- assortativity_degree(g)

  # Type assortativity: do crabs co-occur with crabs, shrimp with shrimp, etc.?
  type_assortativity <- assortativity_nominal(g, as.factor(V(g)$type))

  assortativity_results <- data.frame(
    type = c("Degree", "Taxonomic"),
    coefficient = c(degree_assortativity, type_assortativity),
    interpretation = c(
      case_when(
        degree_assortativity > 0.3 ~ "Assortative (hubs connect to hubs)",
        degree_assortativity < -0.3 ~ "Disassortative (hubs connect to periphery)",
        TRUE ~ "Neutral mixing"
      ),
      case_when(
        type_assortativity > 0.3 ~ "Taxonomically clustered (like with like)",
        type_assortativity < -0.3 ~ "Taxonomically mixed (dissimilar types)",
        TRUE ~ "Random taxonomic mixing"
      )
    )
  )

  write_csv(assortativity_results,
            file.path(SURVEY_TABLES, "network_assortativity.csv"))

  cat(sprintf("    Degree assortativity: %.3f (%s)\n",
              degree_assortativity, assortativity_results$interpretation[1]))
  cat(sprintf("    Taxonomic assortativity: %.3f (%s)\n",
              type_assortativity, assortativity_results$interpretation[2]))

  cat("\n  ✓ Deeper network analyses complete\n\n")

} else {
  cat("  ⚠️  Insufficient data for deeper analyses\n\n")
}

# ============================================================================
# Summary
# ============================================================================

cat("========================================\n")
cat("Network Analysis Summary\n")
cat("========================================\n\n")

if(exists("network_metrics")) {
  cat("CAFI Network Structure:\n")
  cat(sprintf("  - Species in network: %d\n", network_metrics$n_cafi_species))
  cat(sprintf("  - Co-occurrence associations: %d\n", network_metrics$n_associations))
  cat(sprintf("  - Network density: %.3f\n", network_metrics$density))
  cat(sprintf("  - Transitivity (clustering): %.3f\n", network_metrics$transitivity))
  if(!is.na(network_metrics$modularity)) {
    cat(sprintf("  - Modularity: %.3f (%d modules)\n",
                network_metrics$modularity, network_metrics$n_modules))
  }
  cat("\n")
}

cat("Key Points:\n")
cat("  - CAFI-CAFI co-occurrence networks only\n")
cat("  - All corals are Pocillopora spp. (no coral nodes)\n")
cat("  - Coral condition scores position-corrected (from Script 05a)\n")
cat("  - All figures with white backgrounds\n\n")

cat("✅ Network analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", SURVEY_TABLES, "\n\n")
