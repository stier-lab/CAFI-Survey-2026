#!/usr/bin/env Rscript
# ============================================================================
# 02_community_composition.R - Community Composition Analysis
# ============================================================================
#
# PURPOSE: Analyze CAFI community composition across sites and coral types
#
# HYPOTHESIS (H1): CAFI community composition differs among reef sites due to
#   variation in coral landscapes and environmental conditions.
#
# ANALYSES:
#   - Taxonomic summaries and OTU ranking
#   - PERMANOVA on Bray-Curtis dissimilarity (site effect)
#   - Pairwise site comparisons
#   - Community composition by functional group
#
# MANUSCRIPT FIGURES: None directly (supports Fig 4 composition analyses)
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("02: Community Composition Analysis\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")
community_matrix <- load_object("community_matrix.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "community_composition")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Overall Community Composition
# ============================================================================

cat("Analyzing overall community composition...\n")

# CAFI OTU abundance ranking
# NOTE: "species" here refers to morphological OTUs (Operational Taxonomic Units)
# These are field-identified morphotypes, NOT genetically confirmed species
# Use for community metrics (richness, diversity) but interpret cautiously
species_abundance <- cafi_clean %>%
  group_by(species, type) %>%
  summarise(
    total_abundance = n(),
    n_corals = n_distinct(coral_id),
    mean_size = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_abundance)) %>%
  mutate(
    rank = row_number(),
    proportion = total_abundance / sum(total_abundance),
    cumulative_prop = cumsum(proportion)
  )

# Save species abundance table
write_csv(species_abundance,
          file.path(TABLES_DIR, "species_abundance_ranking.csv"))

# Top 20 most abundant species
top20_species <- species_abundance %>%
  slice_head(n = 20)

# Plot CAFI OTU rank abundance
p_rank_abundance <- ggplot(species_abundance, aes(x = rank, y = total_abundance)) +
  geom_line(color = "gray50", linewidth = 0.5) +
  geom_point(aes(color = type, shape = type), size = 2.5, alpha = 0.8) +
  scale_y_log10(labels = scales::comma,
                breaks = c(1, 10, 100, 1000),
                limits = c(1, NA)) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_color_taxon() +
  scale_shape_taxon() +
  labs(
    title = "CAFI community rank-abundance distribution",
    x = "OTU rank",
    y = "Total abundance (log scale)",
    color = "Taxonomic\nGroup"
  ) +
  theme_publication() +
  theme(legend.position = c(0.85, 0.75),
        legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3))

ggsave(file.path(fig_dir, "species_rank_abundance.png"),
       p_rank_abundance, width = 8, height = 5, dpi = 300)

# ============================================================================
# Community Composition by Site
# ============================================================================

cat("Analyzing site-level patterns...\n")

# Site-level community summary
# Comparing 3 major sites:
# - HAU (Hauru): 38 corals, north shore fringing reef
# - MAT (Maatea): 39 corals, interior lagoon
# - MRB (Moorea Barrier Reef): 37 corals, outer barrier reef
site_summary <- cafi_clean %>%
  group_by(site, species, type) %>%
  summarise(
    abundance = n(),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  mutate(
    total_site_abundance = sum(abundance),
    proportion = abundance / total_site_abundance
  )

# Species richness by site
site_richness <- cafi_clean %>%
  group_by(site) %>%
  summarise(
    n_corals = n_distinct(coral_id),
    species_richness = n_distinct(species),
    total_abundance = n(),
    shannon = vegan::diversity(table(species)),
    simpson = vegan::diversity(table(species), index = "simpson"),
    .groups = "drop"
  )

write_csv(site_richness,
          file.path(TABLES_DIR, "site_diversity_metrics.csv"))

# Plot site diversity - using points with error bars based on coral-level variation
coral_diversity <- coral_master %>%
  filter(!is.na(site)) %>%
  group_by(site) %>%
  summarise(
    n_corals = n(),
    mean_richness = mean(otu_richness, na.rm = TRUE),
    se_richness = sd(otu_richness, na.rm = TRUE) / sqrt(n()),
    mean_shannon = mean(shannon_diversity, na.rm = TRUE),
    se_shannon = sd(shannon_diversity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_site_diversity <- coral_diversity %>%
  pivot_longer(cols = c(mean_richness, mean_shannon),
               names_to = "metric",
               values_to = "mean") %>%
  mutate(
    se = if_else(metric == "mean_richness", se_richness, se_shannon),
    metric = factor(metric,
                    levels = c("mean_richness", "mean_shannon"),
                    labels = c("OTU richness", "Shannon diversity (H')"))
  ) %>%
  ggplot(aes(x = site, y = mean, fill = site)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, linewidth = 0.5) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_site() +
  labs(
    title = "CAFI diversity metrics across Mo'orea reef sites",
    x = NULL,
    y = "Mean ± SE per coral"
  ) +
  theme_publication() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11, face = "bold"))

ggsave(file.path(fig_dir, "site_diversity_metrics.png"),
       p_site_diversity, width = 8, height = 4, dpi = 300)

# ============================================================================
# Taxonomic Composition
# ============================================================================

cat("Analyzing taxonomic composition...\n")

# Taxonomic breakdown
taxonomic_summary <- cafi_clean %>%
  group_by(type) %>%
  summarise(
    total_abundance = n(),
    n_species = n_distinct(species),
    n_corals = n_distinct(coral_id),
    mean_size = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(proportion = total_abundance / sum(total_abundance))

write_csv(taxonomic_summary,
          file.path(TABLES_DIR, "taxonomic_composition.csv"))

# Plot taxonomic composition - horizontal bar chart (better than pie for publications)
p_taxonomic_bar <- taxonomic_summary %>%
  mutate(
    type = factor(type, levels = type[order(total_abundance)]),
    label = paste0(round(proportion * 100, 1), "%\n(n=", total_abundance, ")")
  ) %>%
  ggplot(aes(x = type, y = total_abundance, fill = type)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(proportion * 100, 1), "%")),
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_taxon() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Taxonomic composition of CAFI community",
    subtitle = paste0("Total individuals: n = ", sum(taxonomic_summary$total_abundance)),
    x = NULL,
    y = "Total abundance"
  ) +
  theme_publication() +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "bold", size = 11))

ggsave(file.path(fig_dir, "taxonomic_composition_pie.png"),
       p_taxonomic_bar, width = 7, height = 4, dpi = 300)

# ============================================================================
# Top Species Analysis
# ============================================================================

cat("Analyzing dominant species...\n")

# Top 10 species by site
top10_by_site <- site_summary %>%
  group_by(site) %>%
  slice_max(order_by = abundance, n = 10) %>%
  ungroup()

# Heatmap of top species across sites (using ggplot for consistency)
top_species_matrix <- top20_species$species
heatmap_data <- site_summary %>%
  filter(species %in% top_species_matrix) %>%
  dplyr::select(site, species, proportion) %>%
  # Order species by total abundance
  mutate(species = factor(species, levels = rev(top20_species$species)))

p_heatmap <- ggplot(heatmap_data, aes(x = site, y = species, fill = proportion)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(proportion > 0.01,
                                sprintf("%.1f", proportion * 100), "")),
            size = 2.5, color = "white") +
  scale_fill_viridis_c(option = "plasma", name = "Proportion\n(%)",
                       labels = scales::percent_format(scale = 100),
                       limits = c(0, NA)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "Top 20 CAFI OTU distribution across sites",
    x = NULL,
    y = NULL
  ) +
  theme_publication() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 11, face = "bold"),
        legend.position = "right",
        panel.border = element_blank())

ggsave(file.path(fig_dir, "top_species_heatmap.png"),
       p_heatmap, width = 8, height = 7, dpi = 300)

# ============================================================================
# Species Accumulation Curves
# ============================================================================

cat("Calculating species accumulation curves...\n")

# Prepare data for accumulation curves
comm_by_site <- list()
sites <- unique(coral_master$site)

for (s in sites) {
  site_corals <- coral_master %>%
    filter(site == s) %>%
    pull(coral_id)

  if (length(site_corals) > 0) {
    comm_by_site[[s]] <- community_matrix[rownames(community_matrix) %in% site_corals, ]
    # Remove empty species columns
    comm_by_site[[s]] <- comm_by_site[[s]][, colSums(comm_by_site[[s]]) > 0]
  }
}

# Calculate species accumulation
spec_accum_list <- list()
for (s in names(comm_by_site)) {
  if (nrow(comm_by_site[[s]]) > 1) {
    spec_accum_list[[s]] <- specaccum(comm_by_site[[s]], method = "random", permutations = 100)
  }
}

# Plot species accumulation curves using ggplot
if (length(spec_accum_list) > 0) {
  # Convert to data frame for ggplot
  accum_df <- do.call(rbind, lapply(names(spec_accum_list), function(site_name) {
    x <- spec_accum_list[[site_name]]
    if (!is.null(x)) {
      data.frame(
        site = site_name,
        n_corals = x$sites,
        richness = x$richness,
        sd = x$sd
      )
    }
  }))

  if (nrow(accum_df) > 0) {
    p_accum <- ggplot(accum_df, aes(x = n_corals, y = richness, color = site, fill = site)) +
      geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_site() +
      scale_fill_site() +
      labs(
        title = "Species accumulation curves by site",
        subtitle = "Shaded regions show ±1 SD from 100 random permutations",
        x = "Number of corals sampled",
        y = "OTU richness"
      ) +
      theme_publication() +
      theme(legend.position = c(0.85, 0.25),
            legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3))

    ggsave(file.path(fig_dir, "species_accumulation_curves.png"),
           p_accum, width = 7, height = 5, dpi = 300)
  }
} else {
  cat("  Note: Insufficient data for species accumulation curves\n")
}

# ============================================================================
# Community Structure Visualization
# ============================================================================

cat("Visualizing community structure...\n")

# Stacked bar chart of community composition by site
# Get top 10 species overall and group others
top10_overall <- species_abundance %>% slice_head(n = 10) %>% pull(species)

community_grouped <- site_summary %>%
  mutate(species_grouped = if_else(species %in% top10_overall, species, "Other")) %>%
  group_by(site, species_grouped) %>%
  summarise(proportion = sum(proportion), .groups = "drop") %>%
  mutate(species_grouped = factor(species_grouped,
                                   levels = c(top10_overall, "Other")))

p_community_stack <- ggplot(community_grouped,
                             aes(x = site, y = proportion, fill = species_grouped)) +
  geom_col(position = "stack", color = "white", linewidth = 0.2, width = 0.7) +
  scale_fill_manual(values = c(viridis(10), "gray70"), name = "OTU") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(
    title = "Community composition by site",
    subtitle = "Top 10 OTUs shown; remaining grouped as 'Other'",
    x = NULL,
    y = "Proportion of community"
  ) +
  theme_publication() +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 11, face = "bold"))

ggsave(file.path(fig_dir, "community_composition_stacked.png"),
       p_community_stack, width = 9, height = 6, dpi = 300)

# ============================================================================
# Size Distribution Analysis
# ============================================================================

cat("Analyzing size distributions...\n")

# Size distribution by taxonomic group - density plots for better comparison
size_stats <- cafi_clean %>%
  filter(!is.na(size_mm)) %>%
  group_by(type) %>%
  summarise(
    n = n(),
    median = median(size_mm),
    mean = mean(size_mm),
    .groups = "drop"
  )

p_size_distribution <- cafi_clean %>%
  filter(!is.na(size_mm)) %>%
  mutate(type = factor(type, levels = c("crab", "shrimp", "fish", "snail"))) %>%
  ggplot(aes(x = size_mm, fill = type, color = type)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  geom_rug(alpha = 0.3, length = unit(0.02, "npc")) +
  geom_vline(data = size_stats %>% mutate(type = factor(type, levels = c("crab", "shrimp", "fish", "snail"))),
             aes(xintercept = median, color = type),
             linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  scale_fill_taxon() +
  scale_color_taxon() +
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  labs(
    title = "Size distribution of CAFI by taxonomic group",
    subtitle = "Dashed lines indicate median size; rug marks show individual observations",
    x = "Body size (mm)",
    y = "Density",
    fill = "Taxonomic\ngroup",
    color = "Taxonomic\ngroup"
  ) +
  theme_publication() +
  theme(legend.position = c(0.85, 0.7),
        legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3))

ggsave(file.path(fig_dir, "size_distribution_by_type.png"),
       p_size_distribution, width = 8, height = 5, dpi = 300)

# ============================================================================
# Statistical Tests - Kruskal-Wallis for Site Differences
# ============================================================================

cat("Running statistical tests...\n")

# Initialize results data frame
stats_results <- init_results_df()

# Test 1: Kruskal-Wallis for coral-level richness by site
# (PERMANOVA is in script 03 - this tests if individual coral richness differs)
coral_richness_by_site <- coral_master %>%
  filter(!is.na(site), !is.na(otu_richness)) %>%
  select(coral_id, site, otu_richness)

if (nrow(coral_richness_by_site) > 10 && n_distinct(coral_richness_by_site$site) >= 2) {
  kw_richness <- kruskal.test(otu_richness ~ site, data = coral_richness_by_site)

  # Calculate effect size (epsilon-squared for Kruskal-Wallis)
  n_total <- nrow(coral_richness_by_site)
  k <- n_distinct(coral_richness_by_site$site)
  epsilon_sq <- kw_richness$statistic / (n_total - 1)

  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H1a",
      question = "Does CAFI OTU richness differ among sites?",
      test_name = "Kruskal-Wallis",
      test_statistic = kw_richness$statistic,
      df = k - 1,
      p_value = kw_richness$p.value,
      effect_size = as.numeric(epsilon_sq),
      effect_type = "ε²",
      n = n_total,
      notes = "Non-parametric test for coral-level richness"
    )
  )
}

# Test 2: Kruskal-Wallis for coral-level abundance by site
coral_abundance_by_site <- coral_master %>%
  filter(!is.na(site), !is.na(total_cafi)) %>%
  select(coral_id, site, total_cafi)

if (nrow(coral_abundance_by_site) > 10 && n_distinct(coral_abundance_by_site$site) >= 2) {
  kw_abundance <- kruskal.test(total_cafi ~ site, data = coral_abundance_by_site)

  n_total <- nrow(coral_abundance_by_site)
  k <- n_distinct(coral_abundance_by_site$site)
  epsilon_sq <- kw_abundance$statistic / (n_total - 1)

  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H1b",
      question = "Does CAFI abundance differ among sites?",
      test_name = "Kruskal-Wallis",
      test_statistic = kw_abundance$statistic,
      df = k - 1,
      p_value = kw_abundance$p.value,
      effect_size = as.numeric(epsilon_sq),
      effect_type = "ε²",
      n = n_total,
      notes = "Non-parametric test for coral-level abundance"
    )
  )
}

# Test 3: Chi-square test for taxonomic group distribution
taxon_counts <- table(cafi_clean$type)
if (length(taxon_counts) >= 2) {
  chi_taxon <- chisq.test(taxon_counts)

  # Cramér's V effect size
  n_total <- sum(taxon_counts)
  k <- length(taxon_counts)
  cramers_v <- sqrt(chi_taxon$statistic / (n_total * (k - 1)))

  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H1c",
      question = "Is taxonomic group distribution non-uniform?",
      test_name = "Chi-square goodness of fit",
      test_statistic = chi_taxon$statistic,
      df = chi_taxon$parameter,
      p_value = chi_taxon$p.value,
      effect_size = as.numeric(cramers_v),
      effect_type = "Cramér's V",
      n = n_total,
      notes = "Tests if taxonomic groups are equally represented"
    )
  )
}

# Save statistical results
save_stats_summary(stats_results, "02_community_composition", "Community Composition Analysis")

# ============================================================================
# RDA Analysis: Composition vs Coral & Neighborhood Traits (6-Panel Figure)
# ============================================================================

cat("Running db-RDA analysis for composition figure...\n")

# Prepare data for RDA
# Get corals with complete neighborhood data
coral_for_rda <- coral_master %>%
  filter(
    coral_id %in% rownames(community_matrix),
    !is.na(volume),
    !is.na(number_of_neighbors),
    !is.na(mean_neighbor_dist),
    !is.na(total_neighbor_volume)
  ) %>%
  mutate(
    log_vol = log10(volume),
    log_neighbor_vol = log10(total_neighbor_volume + 1)
  )

# Filter community matrix to matching corals
comm_rda <- community_matrix[rownames(community_matrix) %in% coral_for_rda$coral_id, ]
comm_rda <- comm_rda[, colSums(comm_rda) > 0]  # Remove empty species

# Ensure matching order
coral_for_rda <- coral_for_rda[match(rownames(comm_rda), coral_for_rda$coral_id), ]

cat("  - Corals with complete data:", nrow(comm_rda), "\n")

# Only run RDA if we have sufficient data
if (nrow(comm_rda) >= 30) {

  # Convert to presence/absence for Jaccard
  comm_pa <- (comm_rda > 0) * 1

  # Run db-RDA using capscale with Jaccard distance
  rda_full <- capscale(comm_pa ~ log_vol + number_of_neighbors +
                         mean_neighbor_dist + log_neighbor_vol,
                       data = coral_for_rda,
                       distance = "jaccard")

  # Test significance
  rda_anova_full <- anova(rda_full, permutations = 999)
  rda_margin <- anova(rda_full, by = "margin", permutations = 999)

  cat("  - Full model p-value:", rda_anova_full$`Pr(>F)`[1], "\n")

  # Variance partitioning - Size vs Neighborhood
  jacc_dist <- vegdist(comm_pa, method = "jaccard", binary = TRUE)
  var_part_2way <- varpart(jacc_dist,
                            ~ log_vol,
                            ~ number_of_neighbors + mean_neighbor_dist + log_neighbor_vol,
                            data = coral_for_rda)

  # Extract scores
  site_scores_rda <- as.data.frame(scores(rda_full, display = "sites", scaling = 2))
  if ("CAP1" %in% colnames(site_scores_rda)) {
    colnames(site_scores_rda) <- gsub("CAP", "RDA", colnames(site_scores_rda))
  }
  site_scores_rda$coral_id <- rownames(site_scores_rda)
  site_scores_rda <- site_scores_rda %>% left_join(coral_for_rda, by = "coral_id")

  species_scores_rda <- as.data.frame(scores(rda_full, display = "species", scaling = 2))
  if ("CAP1" %in% colnames(species_scores_rda)) {
    colnames(species_scores_rda) <- gsub("CAP", "RDA", colnames(species_scores_rda))
  }
  species_scores_rda$taxon <- rownames(species_scores_rda)

  env_scores_rda <- as.data.frame(scores(rda_full, display = "bp", scaling = 2))
  if ("CAP1" %in% colnames(env_scores_rda)) {
    colnames(env_scores_rda) <- gsub("CAP", "RDA", colnames(env_scores_rda))
  }
  env_scores_rda$variable <- rownames(env_scores_rda)
  env_scores_rda$label <- c("log(Volume)", "# Neighbors", "Mean Distance", "log(Neighbor Vol)")

  # Variance explained
  eig <- rda_full$CCA$eig
  var_rda1 <- eig[1] / sum(eig) * 100
  var_rda2 <- eig[2] / sum(eig) * 100
  total_constrained <- sum(eig) / rda_full$tot.chi * 100

  # Get variance partitioning values
  adj_r2_size <- var_part_2way$part$indfract$Adj.R.squared[1]
  adj_r2_shared <- var_part_2way$part$indfract$Adj.R.squared[3]
  adj_r2_neighborhood <- var_part_2way$part$indfract$Adj.R.squared[2]

  # Predictor importance
  pred_importance <- data.frame(
    predictor = c("Coral Size", "# Neighbors", "Mean Distance", "Neighbor Volume"),
    F_value = rda_margin$F[1:4],
    p_value = rda_margin$`Pr(>F)`[1:4]
  )
  pred_importance$sig <- ifelse(pred_importance$p_value < 0.001, "***",
                                 ifelse(pred_importance$p_value < 0.01, "**",
                                        ifelse(pred_importance$p_value < 0.05, "*", "")))
  pred_importance <- pred_importance %>%
    arrange(desc(F_value)) %>%
    mutate(predictor = factor(predictor, levels = predictor))

  # ==========================================================================
  # Panel A: RDA biplot
  # ==========================================================================

  sp_scale <- 1.5
  env_scale <- 0.35

  p_rda_a <- ggplot() +
    geom_point(data = site_scores_rda,
               aes(x = RDA1, y = RDA2, color = log_vol, size = volume),
               alpha = 0.7) +
    geom_segment(data = species_scores_rda,
                 aes(x = 0, y = 0, xend = RDA1 * sp_scale, yend = RDA2 * sp_scale),
                 arrow = arrow(length = unit(0.15, "cm")),
                 color = "gray30", linewidth = 0.7) +
    geom_text(data = species_scores_rda,
              aes(x = RDA1 * sp_scale * 1.12, y = RDA2 * sp_scale * 1.12, label = taxon),
              size = 2.8, fontface = "bold", color = "gray20") +
    geom_segment(data = env_scores_rda,
                 aes(x = 0, y = 0, xend = RDA1 * env_scale, yend = RDA2 * env_scale),
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "#D55E00", linewidth = 1.1) +
    geom_text(data = env_scores_rda,
              aes(x = RDA1 * env_scale * 1.15, y = RDA2 * env_scale * 1.15, label = label),
              size = 3, fontface = "bold", color = "#D55E00") +
    scale_color_viridis_c(name = expression("log"[10]*"(Vol)"), option = "plasma", end = 0.9) +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    labs(title = "A. RDA biplot: composition vs. coral & neighborhood traits",
         subtitle = sprintf("Model explains %.1f%% of variance (p = %.3f)",
                           total_constrained, rda_anova_full$`Pr(>F)`[1]),
         x = sprintf("RDA1 (%.1f%%)", var_rda1),
         y = sprintf("RDA2 (%.1f%%)", var_rda2)) +
    theme_publication() +
    theme(legend.position = c(0.12, 0.82),
          legend.background = element_rect(fill = alpha("white", 0.9)))

  # ==========================================================================
  # Panel B: Predictor importance
  # ==========================================================================

  p_rda_b <- ggplot(pred_importance, aes(x = predictor, y = F_value)) +
    geom_col(aes(fill = p_value < 0.05), width = 0.7, alpha = 0.85) +
    geom_text(aes(label = sig), vjust = -0.3, size = 5, fontface = "bold") +
    geom_text(aes(label = sprintf("p=%.3f", p_value)), vjust = 1.5, size = 2.8,
              color = "white", fontface = "bold") +
    scale_fill_manual(values = c("TRUE" = "#009E73", "FALSE" = "gray70"),
                      guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "B. Predictor importance (marginal F-test)",
         subtitle = "Unique contribution controlling for other predictors",
         x = NULL,
         y = "F-statistic") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  # ==========================================================================
  # Panel C: Variance partitioning
  # ==========================================================================

  var_df <- data.frame(
    component = c("Size\n(unique)", "Shared", "Neighborhood\n(unique)", "Residual"),
    value = c(max(0, adj_r2_size),
              max(0, adj_r2_shared),
              max(0, adj_r2_neighborhood),
              1 - max(0, adj_r2_size) - max(0, adj_r2_shared) - max(0, adj_r2_neighborhood)) * 100,
    fill_col = c("#0072B2", "#CC79A7", "#D55E00", "gray70")
  )
  var_df$component <- factor(var_df$component, levels = var_df$component)

  p_rda_c <- ggplot(var_df, aes(x = component, y = value, fill = fill_col)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.1f%%", value)), vjust = -0.3, size = 3.5, fontface = "bold") +
    scale_fill_identity() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, 100)) +
    labs(title = "C. Variance partitioning",
         subtitle = "Size vs. all neighborhood metrics combined",
         x = NULL,
         y = "Variance explained (%)") +
    theme_publication()

  # ==========================================================================
  # Panel D: Species loadings
  # ==========================================================================

  # Assign taxon colors based on type (from cafi_clean)
  taxon_types <- cafi_clean %>%
    group_by(species, type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(species) %>%
    slice_max(n) %>%
    select(species, type)

  species_scores_rda <- species_scores_rda %>%
    left_join(taxon_types, by = c("taxon" = "species")) %>%
    mutate(type = if_else(is.na(type), "other", type))

  species_scores_plot <- species_scores_rda %>%
    arrange(RDA1) %>%
    mutate(taxon = factor(taxon, levels = taxon))

  p_rda_d <- ggplot(species_scores_plot, aes(x = taxon, y = RDA1, fill = type)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_taxon() +
    coord_flip() +
    labs(title = "D. Taxon associations with primary gradient",
         subtitle = "RDA1 dominated by coral size effect",
         x = NULL,
         y = "RDA1 score") +
    theme_publication() +
    theme(axis.text.y = element_text(face = "bold", size = 8),
          legend.position = "none")

  # ==========================================================================
  # Panel E: RDA1 vs predictors
  # ==========================================================================

  # Correlation of RDA1 with each predictor
  cors_rda1 <- data.frame(
    predictor = c("log_vol", "number_of_neighbors", "mean_neighbor_dist", "log_neighbor_vol"),
    label = c("log(Volume)", "# Neighbors", "Mean Distance", "log(Neighbor Vol)")
  )
  cors_rda1$r <- sapply(cors_rda1$predictor, function(x) cor(site_scores_rda$RDA1, site_scores_rda[[x]]))
  cors_rda1$p <- sapply(cors_rda1$predictor, function(x) cor.test(site_scores_rda$RDA1, site_scores_rda[[x]])$p.value)

  site_long_rda <- site_scores_rda %>%
    select(coral_id, RDA1, log_vol, number_of_neighbors, mean_neighbor_dist, log_neighbor_vol) %>%
    pivot_longer(cols = c(log_vol, number_of_neighbors, mean_neighbor_dist, log_neighbor_vol),
                 names_to = "predictor", values_to = "value") %>%
    left_join(cors_rda1, by = "predictor") %>%
    mutate(label = factor(label, levels = c("log(Volume)", "# Neighbors", "Mean Distance", "log(Neighbor Vol)")))

  p_rda_e <- ggplot(site_long_rda, aes(x = value, y = RDA1)) +
    geom_point(alpha = 0.5, size = 1.5, color = "#0072B2") +
    geom_smooth(method = "lm", se = TRUE, color = "gray30", linewidth = 0.8) +
    geom_text(data = cors_rda1 %>%
                mutate(label = factor(label, levels = c("log(Volume)", "# Neighbors", "Mean Distance", "log(Neighbor Vol)"))),
              aes(x = -Inf, y = Inf,
                  label = sprintf("r=%.2f%s", r,
                                 ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))),
              hjust = -0.1, vjust = 1.5, size = 3, fontface = "italic", inherit.aes = FALSE) +
    facet_wrap(~ label, scales = "free_x", nrow = 1) +
    labs(title = "E. RDA1 correlation with each predictor",
         subtitle = "Coral size shows strongest association with compositional gradient",
         x = "Predictor value",
         y = "RDA1 score") +
    theme_publication() +
    theme(strip.text = element_text(size = 9))

  # ==========================================================================
  # Panel F: Environmental vector strength
  # ==========================================================================

  env_scores_rda$length <- sqrt(env_scores_rda$RDA1^2 + env_scores_rda$RDA2^2)
  env_scores_rda$angle <- atan2(env_scores_rda$RDA2, env_scores_rda$RDA1) * 180 / pi

  p_rda_f <- ggplot(env_scores_rda, aes(x = reorder(label, -length), y = length)) +
    geom_col(fill = "#D55E00", alpha = 0.8, width = 0.6) +
    geom_text(aes(label = sprintf("%.0f°", angle)), vjust = -0.3, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "F. Environmental vector strength",
         subtitle = "Length = effect magnitude; angle shown in degrees",
         x = NULL,
         y = "Vector length (RDA space)") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  # ==========================================================================
  # Combine 6-panel figure
  # ==========================================================================

  fig_rda_6panel <- (p_rda_a | p_rda_b | p_rda_c) / (p_rda_d | p_rda_f) / p_rda_e +
    plot_layout(heights = c(1.2, 1, 0.8)) +
    plot_annotation(
      title = "Multivariate analysis: Coral size dominates community composition",
      subtitle = sprintf("db-RDA on Jaccard dissimilarity (presence/absence) | Variance explained: %.1f%% (p = %.3f)",
                         total_constrained, rda_anova_full$`Pr(>F)`[1]),
      theme = theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "gray30")
      )
    )

  # Save to community_composition folder
  ggsave(file.path(fig_dir, "rda_composition_6panel.png"), fig_rda_6panel,
         width = 15, height = 13, dpi = 300, bg = "white")
  cat("  - Saved:", file.path(fig_dir, "rda_composition_6panel.png"), "\n")

  # Save to manuscript figures folder (dual output per user request)
  manuscript_dir <- file.path(FIGURES_DIR, "manuscript")
  dir.create(manuscript_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(manuscript_dir, "rda_composition_6panel.png"), fig_rda_6panel,
         width = 15, height = 13, dpi = 300, bg = "white")
  cat("  - Saved:", file.path(manuscript_dir, "rda_composition_6panel.png"), "\n")

  # Print RDA summary
  cat("\n  RDA Summary:\n")
  cat(sprintf("    Total constrained variance: %.1f%%\n", total_constrained))
  cat(sprintf("    Unique to coral size: %.1f%%\n", max(0, adj_r2_size) * 100))
  cat(sprintf("    Unique to neighborhood: %.1f%%\n", max(0, adj_r2_neighborhood) * 100))
  cat(sprintf("    Shared: %.1f%%\n", max(0, adj_r2_shared) * 100))

} else {
  cat("  - Insufficient data for RDA analysis (need >= 30 corals with complete data)\n")
}

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Community Composition Summary\n")
cat("========================================\n\n")

cat("Overall Community:\n")
cat("  - Total OTUs:", n_distinct(cafi_clean$species), "\n")
cat("  - Total individuals:", nrow(cafi_clean), "\n")
cat("  - Most abundant OTU:", top20_species$species[1],
    "(", top20_species$total_abundance[1], "individuals)\n")
cat("  - Dominant taxonomic group:", taxonomic_summary$type[1],
    "(", round(taxonomic_summary$proportion[1] * 100, 1), "%)\n\n")

cat("Site Comparison:\n")
cat("  - Highest richness:", site_richness$site[which.max(site_richness$species_richness)],
    "(", max(site_richness$species_richness), "OTUs)\n")
cat("  - Highest abundance:", site_richness$site[which.max(site_richness$total_abundance)],
    "(", max(site_richness$total_abundance), "individuals)\n")
cat("  - Highest diversity:", site_richness$site[which.max(site_richness$shannon)],
    "(H' =", round(max(site_richness$shannon), 2), ")\n\n")

cat("✅ Community composition analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", TABLES_DIR, "\n")