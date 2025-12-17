#!/usr/bin/env Rscript
# ============================================================================
# 03_diversity_analysis.R - Alpha/Beta Diversity Analysis
# ============================================================================
#
# PURPOSE: Analyze CAFI diversity patterns across corals and sites
#
# HYPOTHESES:
#   H1: Community composition differs among sites (PERMANOVA)
#   H4: Healthier corals support more diverse CAFI (condition-diversity)
#
# ANALYSES:
#   - Alpha diversity: richness, Shannon H', evenness per coral
#   - Beta diversity: NMDS ordination, PERMANOVA, PERMDISP
#   - Diversity partitioning (alpha, beta, gamma)
#   - Condition-diversity relationships
#
# MANUSCRIPT FIGURES:
#   - Supports Figure 4 (NMDS ordination)
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("03: Diversity Analysis\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")
community_matrix <- load_object("community_matrix.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "diversity")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Alpha Diversity
# ============================================================================

cat("Calculating alpha diversity metrics...\n")

# Calculate diversity indices for each coral colony
# These measure within-coral CAFI community diversity
# NOTE: High richness may reflect more cryptic habitat, not host specificity
alpha_diversity <- data.frame(
  coral_id = rownames(community_matrix),
  species_richness = specnumber(community_matrix),  # Number of OTUs
  shannon = vegan::diversity(community_matrix, index = "shannon"),  # Accounts for abundance
  simpson = vegan::diversity(community_matrix, index = "simpson"),  # Probability same species
  evenness = vegan::diversity(community_matrix) / log(specnumber(community_matrix))  # Equitability
)

# Handle corals with no CAFI or only 1 species (evenness = NaN/Inf)
# NaN occurs when richness = 0 (log(0)), Inf/NaN when richness = 1 (log(1) = 0)
alpha_diversity$evenness[is.nan(alpha_diversity$evenness) |
                         is.infinite(alpha_diversity$evenness)] <- NA

# Add metadata - site and depth are key environmental predictors
alpha_diversity <- alpha_diversity %>%
  left_join(coral_master %>% dplyr::select(coral_id, site, depth_m),
            by = "coral_id")

# Save alpha diversity metrics
write_csv(alpha_diversity,
          file.path(TABLES_DIR, "alpha_diversity_metrics.csv"))

# Plot alpha diversity by site
# Comparing 3 major reef zones:
# HAU (fringing), MAT (lagoon), MRB (barrier)
p_alpha_boxplot <- alpha_diversity %>%
  pivot_longer(cols = c(species_richness, shannon, simpson, evenness),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("species_richness", "shannon", "simpson", "evenness"),
                         labels = c("OTU Richness", "Shannon", "Simpson", "Evenness"))) %>%
  ggplot(aes(x = site, y = value, fill = site)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_site() +
  labs(
    title = "Alpha Diversity Metrics by Site",
    subtitle = "Within-coral CAFI community diversity across reef zones",
    x = "Site",
    y = "Value"
  ) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "alpha_diversity_boxplot.png"),
       p_alpha_boxplot, width = 12, height = 8, dpi = 300)

cat("  ✓ Alpha diversity calculated\n\n")

# ============================================================================
# Beta Diversity
# ============================================================================

cat("Calculating beta diversity...\n")

# Remove empty corals (no CAFI observed)
community_matrix_clean <- community_matrix[rowSums(community_matrix, na.rm = TRUE) > 0, ]
community_matrix_clean[is.na(community_matrix_clean)] <- 0

# Calculate beta diversity using multiple distance metrics
# Bray-Curtis: abundance-weighted dissimilarity
# Jaccard: presence/absence only
# Morisita: probabilistic, robust to sample size
beta_bray <- vegdist(community_matrix_clean, method = "bray")
beta_jaccard <- vegdist(community_matrix_clean, method = "jaccard", binary = TRUE)
beta_morisita <- vegdist(community_matrix_clean, method = "morisita")

# Perform NMDS ordination (Non-metric Multidimensional Scaling)
# Reduces community dissimilarity to 2D visualization
# Lower stress (<0.2) indicates good fit
set.seed(123)
nmds_bray <- metaMDS(community_matrix_clean, distance = "bray", k = 2,
                     trymax = 100, autotransform = FALSE)

# Extract NMDS scores and add metadata
nmds_scores <- as.data.frame(scores(nmds_bray, display = "sites"))
nmds_scores$coral_id <- rownames(nmds_scores)
nmds_scores <- nmds_scores %>%
  left_join(coral_master %>% dplyr::select(coral_id, site, depth_m),
            by = "coral_id")

# Plot NMDS with site clustering
# Ellipses show 95% confidence regions
# Points closer together = more similar communities
p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(color = site, fill = site), level = 0.95,
               geom = "polygon", alpha = 0.1, linewidth = 0.8) +
  geom_point(aes(color = site, shape = site), size = 2.5, alpha = 0.8) +
  scale_color_site() +
  scale_fill_site() +
  scale_shape_site() +
  coord_fixed(ratio = 1) +
  labs(
    title = "NMDS ordination of CAFI communities",
    subtitle = sprintf("Bray-Curtis dissimilarity | Stress = %.3f", nmds_bray$stress),
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_publication() +
  theme(legend.position = c(0.12, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3))

ggsave(file.path(fig_dir, "nmds_ordination.png"),
       p_nmds, width = 7, height = 6, dpi = 300)

cat("  ✓ NMDS ordination complete (stress =", round(nmds_bray$stress, 3), ")\n\n")

# ============================================================================
# PCA on Species Composition
# ============================================================================

cat("Running PCA on species composition...\n")

# Perform PCA (Principal Components Analysis) on community matrix
# Reduces 243 CAFI OTUs to main axes of variation
# Use prcomp() with SVD method which handles p >> n case (more species than corals)
pca_species <- prcomp(community_matrix_clean, center = TRUE, scale. = FALSE)

# Extract variance explained by each PC
var_explained <- (pca_species$sdev^2 / sum(pca_species$sdev^2)) * 100

# Create scree plot showing variance explained
scree_data <- data.frame(
  PC = paste0("PC", 1:length(var_explained)),
  Variance = var_explained
) %>%
  mutate(PC = factor(PC, levels = PC))

p_scree <- ggplot(scree_data[1:min(10, nrow(scree_data)),],
                  aes(x = PC, y = Variance)) +
  geom_col(fill = "steelblue", alpha = 0.7, width = 0.7) +
  geom_line(aes(group = 1), color = "darkred", linewidth = 0.8) +
  geom_point(color = "darkred", size = 2.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "PCA scree plot",
    x = "Principal component",
    y = "Variance explained (%)"
  ) +
  theme_publication()

# Extract PCA scores (coral positions in PC space)
# prcomp() stores scores as $x (not $scores like princomp)
pca_scores <- as.data.frame(pca_species$x[, 1:2])
pca_scores$coral_id <- rownames(pca_scores)
pca_scores <- pca_scores %>%
  left_join(coral_master %>% dplyr::select(coral_id, site, depth_m),
            by = "coral_id")

# Plot PCA scores colored by site
p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(color = site, fill = site), level = 0.95,
               geom = "polygon", alpha = 0.1, linewidth = 0.8) +
  geom_point(aes(color = site, shape = site), size = 2.5, alpha = 0.8) +
  scale_color_site() +
  scale_fill_site() +
  scale_shape_site() +
  labs(
    title = "PCA of CAFI species composition",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_publication() +
  theme(legend.position = c(0.12, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.3))

# Create biplot with species loadings (which OTUs drive PC axes)
# prcomp() stores loadings as $rotation (not $loadings like princomp)
pca_loadings <- as.data.frame(pca_species$rotation[, 1:2])
pca_loadings$species <- rownames(pca_loadings)
# Scale loadings for visualization
pca_loadings <- pca_loadings %>%
  mutate(
    PC1 = PC1 * 5,
    PC2 = PC2 * 5
  )

# Only show top loadings to avoid overcrowding
# Reduced to top 10 for readability
top_loadings <- pca_loadings %>%
  mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
  arrange(desc(magnitude)) %>%
  head(10)

p_biplot <- p_pca +
  geom_segment(data = top_loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               color = "gray30", alpha = 0.7, linewidth = 0.6) +
  ggrepel::geom_text_repel(data = top_loadings,
                           aes(x = PC1, y = PC2, label = species),
                           size = 2.8, alpha = 0.9, fontface = "italic",
                           box.padding = 0.4,
                           max.overlaps = 20,
                           segment.color = "gray60",
                           segment.size = 0.3) +
  labs(title = "PCA biplot of CAFI communities",
       subtitle = "Arrows indicate OTUs with strongest influence on ordination axes") +
  theme(legend.position = c(0.12, 0.15))

# Combine scree and PCA plots
p_pca_combined <- p_scree + p_pca +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    title = "PCA analysis of CAFI community composition",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

ggsave(file.path(fig_dir, "pca_scree_scores.png"),
       p_pca_combined, width = 12, height = 5, dpi = 300)

ggsave(file.path(fig_dir, "pca_biplot.png"),
       p_biplot, width = 8, height = 7, dpi = 300)

# Save PCA results
write.csv(
  data.frame(
    PC = paste0("PC", 1:length(var_explained)),
    Variance_Percent = var_explained,
    Cumulative_Percent = cumsum(var_explained)
  ),
  file.path(TABLES_DIR, "pca_variance_explained.csv"),
  row.names = FALSE
)

write.csv(
  pca_loadings %>% dplyr::select(species, PC1, PC2),
  file.path(TABLES_DIR, "pca_species_loadings.csv"),
  row.names = FALSE
)

cat("✓ PCA analysis complete\n\n")

# ============================================================================
# PERMANOVA Analysis
# ============================================================================

cat("Running PERMANOVA tests...\n")

# PERMANOVA = Permutational Multivariate Analysis of Variance
# Tests if community composition differs among groups (sites, depth)
# Uses permutations to generate null distribution

# Filter out rows with NA values in predictors
# CRITICAL: Ensure exact row alignment between community matrix and metadata
coral_ids_in_matrix <- rownames(community_matrix)
metadata_complete <- coral_master %>%
  filter(coral_id %in% coral_ids_in_matrix, !is.na(site), !is.na(depth_m)) %>%
  arrange(match(coral_id, coral_ids_in_matrix))

# Subset community matrix to matching rows
community_matrix_complete <- community_matrix[metadata_complete$coral_id, , drop = FALSE]

# Verify alignment
stopifnot(all(rownames(community_matrix_complete) == metadata_complete$coral_id))

# Test for differences between sites (primary predictor)
# Are HAU, MAT, MRB communities different?
permanova_site <- adonis2(community_matrix_complete ~ site,
                          data = metadata_complete,
                          method = "bray",
                          permutations = 999)

# Combined model with site and depth
# Tests independent effects
permanova_combined <- adonis2(community_matrix_complete ~ site + depth_m,
                              data = metadata_complete,
                              method = "bray",
                              permutations = 999)

# Save PERMANOVA results
capture.output(
  cat("PERMANOVA Results - Community Composition Analysis\n"),
  cat("================================================\n\n"),
  cat("Site Effect:\n"),
  cat("Tests if CAFI communities differ among HAU, MAT, MRB\n"),
  print(permanova_site),
  cat("\n\nCombined Model (Site + Depth):\n"),
  cat("Tests independent effects of site and depth\n"),
  print(permanova_combined),
  file = file.path(TABLES_DIR, "permanova_results.txt")
)

cat("  ✓ PERMANOVA tests complete\n\n")

# ============================================================================
# Dispersion Analysis (PERMDISP)
# ============================================================================

cat("Testing multivariate dispersion...\n")

# PERMDISP tests for homogeneity of multivariate dispersions
# Significant result = groups have different spread (heteroscedasticity)
# Important complement to PERMANOVA

# Ensure metadata alignment with cleaned community matrix
# CRITICAL: Must match exact rows and order
clean_coral_ids <- rownames(community_matrix_clean)
metadata_clean <- coral_master %>%
  filter(coral_id %in% clean_coral_ids)

# Reorder to match community matrix row order
metadata_clean <- metadata_clean[match(clean_coral_ids, metadata_clean$coral_id), ]

# Filter to only rows that exist in both (some coral_ids may be missing from coral_master)
valid_rows <- !is.na(metadata_clean$coral_id)
metadata_clean <- metadata_clean[valid_rows, ]
community_matrix_clean <- community_matrix_clean[valid_rows, ]
beta_bray <- vegdist(community_matrix_clean, method = "bray")

betadisp_site <- betadisper(beta_bray, metadata_clean$site)
permdisp_site <- permutest(betadisp_site, permutations = 999)

# Plot dispersion with white background
png(file.path(fig_dir, "beta_dispersion.png"),
    width = 10, height = 6, units = "in", res = 300, bg = "white")
par(bg = "white")
plot(betadisp_site, main = "Multivariate Dispersion by Site",
     sub = "Tests homogeneity of community variance")
dev.off()

cat("  ✓ Dispersion analysis complete\n\n")

# ============================================================================
# Rarefaction Analysis
# ============================================================================

cat("Performing rarefaction analysis...\n")

# Rarefaction curves show species accumulation with sampling effort
# Plateaus indicate adequate sampling
# Steep curves suggest under-sampling
rarecurve_data <- rarecurve(community_matrix, step = 1, label = FALSE)

# Create rarefaction plot with white background
png(file.path(fig_dir, "rarefaction_curves.png"),
    width = 10, height = 6, units = "in", res = 300, bg = "white")
par(bg = "white")
rarecurve(community_matrix, step = 1, col = viridis(nrow(community_matrix)),
          main = "OTU Rarefaction Curves",
          sub = "Each line = one coral colony",
          xlab = "Number of Individuals",
          ylab = "OTU Richness")
dev.off()

cat("  ✓ Rarefaction complete\n\n")

# ============================================================================
# Indicator Species Analysis
# ============================================================================

cat("Identifying indicator species...\n")

# Indicator species analysis finds OTUs strongly associated with sites
# High indicator value = OTU is frequent AND abundant in that site

# Get site assignments for community matrix rows, removing any NAs
site_assignments <- coral_master$site[match(rownames(community_matrix),
                                            coral_master$coral_id)]

# Filter to rows with valid site assignments
valid_site_rows <- !is.na(site_assignments)
if (sum(valid_site_rows) > 10) {
  comm_for_indval <- community_matrix[valid_site_rows, ]
  sites_for_indval <- site_assignments[valid_site_rows]

  indval_site <- multipatt(comm_for_indval, sites_for_indval,
                           control = how(nperm = 999))

  # Extract significant indicators (p < 0.05)
  indicators <- indval_site$sign %>%
    filter(p.value < 0.05) %>%
    mutate(species = rownames(.))

  if (nrow(indicators) > 0) {
    write_csv(indicators,
              file.path(TABLES_DIR, "indicator_species.csv"))

    cat("  - Found", nrow(indicators), "indicator OTUs (p < 0.05)\n")
  } else {
    cat("  - No significant indicator species found\n")
  }
} else {
  cat("  - Insufficient data for indicator species analysis\n")
}

cat("  ✓ Indicator analysis complete\n\n")

# ============================================================================
# Diversity Partitioning
# ============================================================================

cat("Partitioning diversity components...\n")

# Diversity partitioning decomposes gamma (total) diversity into:
# - Alpha (within-site): local diversity
# - Beta (between-site): turnover/differentiation

# Gamma diversity (total regional diversity)
gamma_div <- vegan::diversity(colSums(community_matrix), index = "shannon")

# Try diversity partitioning if site information available
tryCatch({
  if ("site" %in% colnames(coral_master)) {
    # Alpha diversity (mean within sites)
    # Join diversity metrics with site information
    alpha_with_site <- merge(alpha_diversity,
                             coral_master[, c("coral_id", "site")],
                             by = "coral_id",
                             all.x = TRUE)

    # Remove rows with missing site
    alpha_with_site <- alpha_with_site[!is.na(alpha_with_site$site), ]

    if (nrow(alpha_with_site) > 0 && "site" %in% colnames(alpha_with_site)) {
      # Calculate mean alpha diversity per site
      alpha_by_site <- aggregate(shannon ~ site,
                                 data = alpha_with_site,
                                 FUN = mean,
                                 na.rm = TRUE)
      colnames(alpha_by_site)[2] <- "mean_alpha"

      mean_alpha <- mean(alpha_by_site$mean_alpha)

      # Beta diversity (multiplicative partitioning)
      # How many times more diverse is the region than average site?
      beta_div <- gamma_div / mean_alpha

      diversity_partition <- data.frame(
        gamma_diversity = gamma_div,
        mean_alpha_diversity = mean_alpha,
        beta_diversity = beta_div,
        proportion_within = mean_alpha / gamma_div,
        proportion_between = 1 - (mean_alpha / gamma_div)
      )

      write_csv(diversity_partition,
                file.path(TABLES_DIR, "diversity_partitioning.csv"))

      cat("  ✓ Diversity partitioning complete\n\n")
    } else {
      cat("  ⚠ Skipping diversity partitioning (site data not available)\n\n")
    }
  } else {
    cat("  ⚠ Skipping diversity partitioning (site column not found in metadata)\n\n")
  }
}, error = function(e) {
  cat("  ⚠ Skipping diversity partitioning (error:", conditionMessage(e), ")\n\n")
})

# ============================================================================
# Diversity Correlations
# ============================================================================

cat("Analyzing diversity correlates...\n")

# Test relationships between diversity and environmental gradients
div_correlates <- alpha_diversity %>%
  filter(!is.na(depth_m))

if (nrow(div_correlates) > 10) {
  # Diversity vs depth gradient - separate panels for clarity
  div_long <- div_correlates %>%
    dplyr::select(coral_id, site, depth_m, species_richness, shannon) %>%
    pivot_longer(cols = c(species_richness, shannon),
                 names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric,
                           levels = c("species_richness", "shannon"),
                           labels = c("OTU richness", "Shannon diversity (H')")))

  p_div_depth <- ggplot(div_long, aes(x = depth_m, y = value)) +
    geom_point(aes(color = site, shape = site), size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", color = "gray30", linewidth = 0.8, se = TRUE, alpha = 0.2) +
    facet_wrap(~metric, scales = "free_y") +
    scale_color_site() +
    scale_shape_site() +
    labs(
      title = "CAFI diversity across depth gradient",
      x = "Depth (m)",
      y = NULL
    ) +
    theme_publication() +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, "diversity_vs_depth.png"),
         p_div_depth, width = 9, height = 4.5, dpi = 300)

  cat("  ✓ Depth correlations analyzed\n\n")
}

# ============================================================================
# Compile Statistical Results Summary
# ============================================================================

cat("Compiling statistical results...\n")

# Initialize results data frame
stats_results <- init_results_df()

# Result 1: PERMANOVA - Site effect on community composition
# Note: adonis2 returns F with NAs for rows that don't have F values
perm_F <- permanova_site$F[1]
perm_p <- permanova_site$`Pr(>F)`[1]
perm_R2 <- permanova_site$R2[1]

stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "H1",
    question = "Does CAFI community composition differ among sites?",
    test_name = "PERMANOVA (Bray-Curtis)",
    test_statistic = ifelse(is.na(perm_F), 0, perm_F),
    df = paste(permanova_site$Df[1], permanova_site$Df[2], sep = ", "),
    p_value = ifelse(is.na(perm_p), 1, perm_p),
    effect_size = ifelse(is.na(perm_R2), 0, perm_R2),
    effect_type = "R²",
    n = nrow(community_matrix_complete),
    notes = "999 permutations"
  )
)

# Result 2: PERMANOVA - Combined model (site + depth)
perm_comb_F <- permanova_combined$F[1]
perm_comb_p <- permanova_combined$`Pr(>F)`[1]
perm_comb_R2 <- permanova_combined$R2[1]

stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "H1",
    question = "Do site and depth independently affect community composition?",
    test_name = "PERMANOVA (Site + Depth)",
    test_statistic = ifelse(is.na(perm_comb_F), 0, perm_comb_F),
    df = paste(permanova_combined$Df[1], permanova_combined$Df[2], sep = ", "),
    p_value = ifelse(is.na(perm_comb_p), 1, perm_comb_p),
    effect_size = ifelse(is.na(perm_comb_R2), 0, perm_comb_R2),
    effect_type = "R²",
    n = nrow(community_matrix_complete),
    notes = "Combined model, 999 permutations"
  )
)

# Result 3: PERMDISP - Dispersion homogeneity
# permdisp returns list with tab (anova table) containing statistics
disp_F <- permdisp_site$tab$F[1]
disp_p <- permdisp_site$tab$`Pr(>F)`[1]

stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "H1-assumption",
    question = "Is multivariate dispersion homogeneous among sites?",
    test_name = "PERMDISP",
    test_statistic = ifelse(is.null(disp_F) || is.na(disp_F), 0, disp_F),
    df = paste(permdisp_site$tab$Df[1], permdisp_site$tab$Df[2], sep = ", "),
    p_value = ifelse(is.null(disp_p) || is.na(disp_p), 1, disp_p),
    effect_size = NA,
    effect_type = "NA",
    n = nrow(community_matrix_clean),
    notes = "Tests PERMANOVA assumption; NS = dispersion equal"
  )
)

# Result 4: NMDS stress
stats_results <- bind_rows(stats_results,
  create_result_row(
    hypothesis = "Ordination",
    question = "How well does 2D ordination represent community dissimilarity?",
    test_name = "NMDS Stress",
    test_statistic = nmds_bray$stress,
    df = "2 dimensions",
    p_value = NA,
    effect_size = nmds_bray$stress,
    effect_type = "Stress",
    n = nrow(community_matrix_clean),
    notes = ifelse(nmds_bray$stress < 0.1, "Excellent fit",
            ifelse(nmds_bray$stress < 0.2, "Good fit", "Acceptable fit"))
  )
)

# Result 5: Depth-richness correlation
if (nrow(div_correlates) > 10) {
  depth_cor <- cor.test(div_correlates$depth_m, div_correlates$species_richness,
                        method = "spearman", exact = FALSE)
  # Spearman doesn't have df parameter, use n-2 as approximation
  df_value <- ifelse(is.null(depth_cor$parameter), nrow(div_correlates) - 2, depth_cor$parameter)
  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H1d",
      question = "Does CAFI richness correlate with depth?",
      test_name = "Spearman correlation",
      test_statistic = as.numeric(depth_cor$statistic),
      df = as.character(df_value),
      p_value = depth_cor$p.value,
      effect_size = as.numeric(depth_cor$estimate),
      effect_type = "rho",
      n = nrow(div_correlates),
      notes = "Depth gradient effect"
    )
  )
}

# Result 6: ANOVA for alpha diversity by site
alpha_for_anova <- alpha_diversity %>% filter(!is.na(site))
if (nrow(alpha_for_anova) > 10 && n_distinct(alpha_for_anova$site) >= 2) {
  aov_richness <- aov(species_richness ~ site, data = alpha_for_anova)
  aov_summary <- summary(aov_richness)[[1]]
  eta_sq <- aov_summary["Sum Sq"][1, 1] / sum(aov_summary["Sum Sq"])

  stats_results <- bind_rows(stats_results,
    create_result_row(
      hypothesis = "H1e",
      question = "Does mean coral-level richness differ among sites?",
      test_name = "One-way ANOVA",
      test_statistic = aov_summary["F value"][1, 1],
      df = paste(aov_summary["Df"][1, 1], aov_summary["Df"][2, 1], sep = ", "),
      p_value = aov_summary["Pr(>F)"][1, 1],
      effect_size = eta_sq,
      effect_type = "η²",
      n = nrow(alpha_for_anova),
      notes = "Alpha diversity by site"
    )
  )
}

# Save statistical results
save_stats_summary(stats_results, "03_diversity_analysis", "Alpha/Beta Diversity Analysis")

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Diversity Analysis Summary\n")
cat("========================================\n\n")

cat("Alpha Diversity (within-coral):\n")
cat("  - Mean OTU richness:", round(mean(alpha_diversity$species_richness), 1), "\n")
cat("  - Mean Shannon diversity:", round(mean(alpha_diversity$shannon), 2), "\n")
cat("  - Mean Simpson diversity:", round(mean(alpha_diversity$simpson), 3), "\n")
cat("  - Mean evenness:", round(mean(alpha_diversity$evenness, na.rm = TRUE), 3), "\n\n")

cat("Beta Diversity (between-coral):\n")
cat("  - NMDS stress:", round(nmds_bray$stress, 3),
    ifelse(nmds_bray$stress < 0.1, "(excellent)",
    ifelse(nmds_bray$stress < 0.2, "(good)", "(acceptable)")), "\n")
cat("  - Mean Bray-Curtis dissimilarity:", round(mean(beta_bray), 3), "\n\n")

cat("PERMANOVA Results:\n")
cat("  - Site effect: F =", round(permanova_site$F[1], 2),
    ", R² =", round(permanova_site$R2[1], 3),
    ", p", format_p(permanova_site$`Pr(>F)`[1]), "\n\n")

cat("Diversity Partitioning:\n")
cat("  - Gamma (regional):", round(gamma_div, 2), "\n")
if (exists("mean_alpha") && exists("beta_div") && exists("diversity_partition")) {
  cat("  - Mean alpha (local):", round(mean_alpha, 2), "\n")
  cat("  - Beta (turnover):", round(beta_div, 2), "×\n")
  cat("  - Within-site component:", round(diversity_partition$proportion_within * 100, 1), "%\n")
  cat("  - Between-site component:", round(diversity_partition$proportion_between * 100, 1), "%\n\n")
} else {
  cat("  - Site-level partitioning: not calculated (site data unavailable)\n\n")
}

cat("✅ Diversity analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", TABLES_DIR, "\n")
