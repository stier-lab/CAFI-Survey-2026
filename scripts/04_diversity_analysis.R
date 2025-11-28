#!/usr/bin/env Rscript
# ============================================================================
# 04_diversity_analysis.R - Diversity and Condition Analysis (H1, H4)
#
# Hypothesis H1: CAFI community composition differs among reef sites due to
# variation in coral landscapes and environmental conditions.
#   Tests: PERMANOVA, NMDS ordination, beta diversity partitioning
#
# Hypothesis H4: Coral physiological condition positively predicts CAFI
# diversity, consistent with bidirectional coral-CAFI interactions.
#   Model: Shannon ~ Condition + log(Volume) + Branch Width + (1|Site)
#
# Theoretical Background:
#   CAFI provide benefits to corals (predator defense, sediment removal,
#   nutrient cycling) while corals provide habitat. This bidirectional
#   relationship predicts that healthier corals support more diverse CAFI
#   communities, creating positive feedbacks.
#
# IMPORTANT TAXONOMIC NOTES:
# - CAFI "species" are morphological OTUs (Operational Taxonomic Units)
# - NO genetic verification - field identifications only
# - Use diversity metrics (richness, Shannon) but avoid species-level inferences
# - Coral "morphotypes" (meandrina/eydoxi/verucosa) are NOT confirmed species
# - branch_width (tight/wide) is the real measurable coral trait
#
# Sites: HAU = Hauru, MAT = Maatea, MRB = Moorea Barrier Reef
#
# Key Analyses:
# - Alpha diversity (within-coral richness and evenness)
# - Beta diversity (between-coral dissimilarity via NMDS, PCA)
# - PERMANOVA (test for site/trait effects on community composition)
# - Condition-diversity relationships
# - Diversity partitioning (alpha, beta, gamma components)
#
# Author: CAFI Analysis Pipeline
# Date: 2025-10-31
# ============================================================================

cat("\n========================================\n")
cat("Survey Diversity Analysis\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/00_load_libraries.R"))

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))

# Create figure and table subdirectories
fig_dir <- file.path(SURVEY_FIGURES, "diversity")
table_dir <- SURVEY_TABLES
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
  left_join(metadata %>% select(coral_id, site, depth_m),
            by = "coral_id")

# Save alpha diversity metrics
write_csv(alpha_diversity,
          file.path(SURVEY_TABLES, "alpha_diversity_metrics.csv"))

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
  left_join(metadata %>% select(coral_id, site, depth_m),
            by = "coral_id")

# Plot NMDS with site clustering
# Ellipses show 95% confidence regions
# Points closer together = more similar communities
p_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site), size = 3, alpha = 0.7) +
  stat_ellipse(aes(color = site), level = 0.95) +
  scale_color_site() +
  labs(
    title = "NMDS Ordination of CAFI Communities",
    subtitle = paste("Stress =", round(nmds_bray$stress, 3),
                     "| Lower stress = better 2D representation"),
    x = "NMDS1",
    y = "NMDS2",
    color = "Site"
  ) +
  theme_publication()

ggsave(file.path(fig_dir, "nmds_ordination.png"),
       p_nmds, width = 10, height = 8, dpi = 300)

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
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_line(aes(group = 1), color = "darkred", linewidth = 1) +
  geom_point(color = "darkred", size = 3) +
  labs(
    title = "PCA Scree Plot - Species Composition",
    subtitle = paste0("PC1: ", round(var_explained[1], 1), "%, ",
                     "PC2: ", round(var_explained[2], 1), "%"),
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_bw()

# Extract PCA scores (coral positions in PC space)
# prcomp() stores scores as $x (not $scores like princomp)
pca_scores <- as.data.frame(pca_species$x[, 1:2])
pca_scores$coral_id <- rownames(pca_scores)
pca_scores <- pca_scores %>%
  left_join(metadata %>% select(coral_id, site, depth_m),
            by = "coral_id")

# Plot PCA scores colored by site
p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = site), size = 3, alpha = 0.7) +
  stat_ellipse(aes(color = site), level = 0.95) +
  scale_color_site() +
  labs(
    title = "PCA of CAFI Species Composition",
    subtitle = paste0("PC1: ", round(var_explained[1], 1), "%, ",
                     "PC2: ", round(var_explained[2], 1), "%"),
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    color = "Site"
  ) +
  theme_publication()

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
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black", alpha = 0.5, linewidth = 0.8) +
  ggrepel::geom_text_repel(data = top_loadings,
                           aes(x = PC1, y = PC2, label = species),
                           size = 3, alpha = 0.8,
                           box.padding = 0.5,
                           max.overlaps = 20) +
  labs(title = "PCA Biplot - Top 10 OTU Loadings and Coral Scores",
       subtitle = "Arrows show OTUs with strongest influence on PC axes")

# Combine scree and PCA plots
p_pca_combined <- p_scree + p_pca +
  plot_annotation(
    title = "PCA Analysis of CAFI Community Composition",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(file.path(fig_dir, "pca_scree_scores.png"),
       p_pca_combined, width = 14, height = 6, dpi = 300)

ggsave(file.path(fig_dir, "pca_biplot.png"),
       p_biplot, width = 12, height = 10, dpi = 300)

# Save PCA results
write.csv(
  data.frame(
    PC = paste0("PC", 1:length(var_explained)),
    Variance_Percent = var_explained,
    Cumulative_Percent = cumsum(var_explained)
  ),
  file.path(table_dir, "pca_variance_explained.csv"),
  row.names = FALSE
)

write.csv(
  pca_loadings %>% select(species, PC1, PC2),
  file.path(table_dir, "pca_species_loadings.csv"),
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
metadata_complete <- metadata %>%
  filter(!is.na(site), !is.na(depth_m))

community_matrix_complete <- community_matrix[rownames(community_matrix) %in% metadata_complete$coral_id, ]

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
  file = file.path(SURVEY_TABLES, "permanova_results.txt")
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
metadata_clean <- metadata %>%
  filter(coral_id %in% rownames(community_matrix_clean)) %>%
  arrange(match(coral_id, rownames(community_matrix_clean)))
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
indval_site <- multipatt(community_matrix,
                         metadata$site,
                         control = how(nperm = 999))

# Extract significant indicators (p < 0.05)
indicators <- indval_site$sign %>%
  filter(p.value < 0.05) %>%
  mutate(species = rownames(.))

if (nrow(indicators) > 0) {
  write_csv(indicators,
            file.path(SURVEY_TABLES, "indicator_species.csv"))

  cat("  - Found", nrow(indicators), "indicator OTUs (p < 0.05)\n")
} else {
  cat("  - No significant indicator species found\n")
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
  if ("site" %in% colnames(metadata)) {
    # Alpha diversity (mean within sites)
    # Join diversity metrics with site information
    alpha_with_site <- merge(alpha_diversity,
                             metadata[, c("coral_id", "site")],
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
                file.path(SURVEY_TABLES, "diversity_partitioning.csv"))

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
  # Diversity vs depth gradient
  # Tests if deeper/shallower corals have different CAFI diversity
  p_div_depth <- ggplot(div_correlates, aes(x = depth_m)) +
    geom_point(aes(y = species_richness, color = "Richness"), size = 2, alpha = 0.6) +
    geom_smooth(aes(y = species_richness, color = "Richness"),
                method = "loess", se = TRUE) +
    geom_point(aes(y = shannon * 10, color = "Shannon (×10)"), size = 2, alpha = 0.6) +
    geom_smooth(aes(y = shannon * 10, color = "Shannon (×10)"),
                method = "loess", se = TRUE) +
    scale_color_manual(values = c("Richness" = "#0072B2", "Shannon (×10)" = "#D55E00")) +
    labs(
      title = "Diversity Metrics vs Depth",
      subtitle = "Testing depth gradient effects on CAFI diversity",
      x = "Depth (m)",
      y = "Diversity Value",
      color = "Metric"
    ) +
    theme_publication()

  ggsave(file.path(fig_dir, "diversity_vs_depth.png"),
         p_div_depth, width = 10, height = 6, dpi = 300)

  cat("  ✓ Depth correlations analyzed\n\n")
}

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
cat("  - Site effect p-value:", format.pval(permanova_site$`Pr(>F)`[1], digits = 3), "\n\n")

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
cat("Tables saved to:", SURVEY_TABLES, "\n")
