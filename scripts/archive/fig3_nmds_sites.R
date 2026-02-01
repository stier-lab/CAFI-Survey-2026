# ============================================================================
# fig3_nmds_sites.R - Figure 3: NMDS by Reef Site with Species Vectors
# ============================================================================
#
# PURPOSE: Show that reef-scale species pools drive CAFI composition
#          Excludes low-sample corals (< 5 CAFI) that create sampling artifacts
#
# OUTPUT: output/figures/manuscript/fig3_composition.png
#
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    FIGURE 3: NMDS by Reef Site\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

source(here::here("scripts/00_setup.R"))
library(ggrepel)

# Load data
coral_master <- load_object("coral_master")
community_matrix <- load_object("community_matrix")

cat("[OK] Data loaded\n")
cat("     Corals:", nrow(coral_master), "\n")
cat("     Species:", ncol(community_matrix), "\n\n")

# ============================================================================
# DATA PREPARATION
# ============================================================================

# Align community matrix with coral data
common_ids <- intersect(rownames(community_matrix), coral_master$coral_id)
comm_aligned <- community_matrix[common_ids, ]
coral_aligned <- coral_master %>%
  filter(coral_id %in% common_ids) %>%
  arrange(match(coral_id, common_ids))

# Remove empty species columns
comm_aligned <- comm_aligned[, colSums(comm_aligned) > 0]

# Calculate total CAFI per coral
coral_aligned$total_cafi_check <- rowSums(comm_aligned)

# Identify and EXCLUDE low-sample corals (< 5 CAFI individuals)
min_cafi <- 5
low_sample <- coral_aligned$total_cafi_check < min_cafi
n_excluded <- sum(low_sample)

cat("Corals with <", min_cafi, "CAFI individuals:", n_excluded, "\n")
cat("  Excluded:", coral_aligned$coral_id[low_sample], "\n")
cat("  These create sampling artifacts - excluding completely\n\n")

# Filter to well-sampled corals only
comm_aligned <- comm_aligned[!low_sample, ]
coral_aligned <- coral_aligned[!low_sample, ]

# Remove any species that are now absent
comm_aligned <- comm_aligned[, colSums(comm_aligned) > 0]
cat("After exclusion: n =", nrow(comm_aligned), "corals,", ncol(comm_aligned), "species\n\n")

# Site colors (consistent with other figures)
site_colors <- c(HAU = "#E69F00", MAT = "#0072B2", MRB = "#009E73")

# ============================================================================
# NMDS ON WELL-SAMPLED CORALS ONLY
# ============================================================================

cat("Running NMDS on well-sampled corals (n >=", min_cafi, ")...\n")

# Hellinger transformation
comm_hell <- decostand(comm_aligned, method = "hellinger")

set.seed(42)
nmds <- metaMDS(comm_hell, distance = "bray", k = 2, trymax = 100, trace = 0)

cat("NMDS stress:", round(nmds$stress, 3), "\n")
if (nmds$stress < 0.1) {
  cat("  Interpretation: Excellent representation\n\n")
} else if (nmds$stress < 0.2) {
  cat("  Interpretation: Good representation\n\n")
} else {
  cat("  Interpretation: Fair representation\n\n")
}

# Extract scores
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$coral_id <- rownames(nmds_scores)
nmds_scores <- nmds_scores %>%
  left_join(coral_aligned %>% dplyr::select(coral_id, volume_field, site, total_cafi, species_richness),
            by = "coral_id")
# ============================================================================
# SPECIES VECTORS
# ============================================================================

cat("Fitting species vectors...\n")

# Filter to species present in at least 8 corals
species_freq <- colSums(comm_aligned > 0)
common_species <- names(species_freq[species_freq >= 8])
cat("  Species present in >=8 corals:", length(common_species), "\n")

# Fit species scores using envfit
set.seed(42)
species_fit <- envfit(nmds, comm_hell[, common_species], permutations = 999)

# Extract significant species vectors
species_scores_df <- as.data.frame(scores(species_fit, display = "vectors"))
species_scores_df$species <- rownames(species_scores_df)
species_scores_df$r2 <- species_fit$vectors$r
species_scores_df$pval <- species_fit$vectors$pvals

# Filter to significant and strong vectors
sig_species <- species_scores_df %>%
  filter(pval < 0.01, r2 > 0.12) %>%
  arrange(desc(r2))

cat("  Significant vectors (p < 0.01, R² > 0.12):", nrow(sig_species), "\n\n")
print(sig_species %>% dplyr::select(species, r2, pval))
cat("\n")

# Scale vectors for plotting
arrow_scale <- 1.2
sig_species <- sig_species %>%
  mutate(
    NMDS1_scaled = NMDS1 * arrow_scale,
    NMDS2_scaled = NMDS2 * arrow_scale,
    species_label = gsub("_", " ", species)
  )

# Take top 8 vectors
top_vectors <- sig_species %>% slice_head(n = 8)

# ============================================================================
# PERMANOVA & PERMDISP: Test site effect on location AND dispersion
# Following Anderson (2005) PERMANOVA+ guide and Anderson et al. (2011) Ecology Letters
# ============================================================================
#
# Key insight from Anderson et al. (2011): PERMANOVA tests for differences in
# multivariate location (centroids), but is ALSO sensitive to differences in
# multivariate dispersion (spread). A significant PERMANOVA could indicate:
#   (a) Groups differ in location (centroid position)
#   (b) Groups differ in dispersion (within-group variation)
#   (c) Both
#
# PERMDISP (betadisper) tests specifically for dispersion differences.
# Together, these tests allow proper interpretation of site effects.
#
# Transformation choice (Hellinger) follows Legendre & Gallagher (2001):
# - Handles double-zero problem (two sites not more similar for shared absences)
# - Emphasizes compositional differences over raw abundance
# - Bray-Curtis on Hellinger-transformed data excludes joint absences (appropriate here)
# ============================================================================

cat("============================================================\n")
cat("PERMANOVA & PERMDISP Analysis (Anderson 2005, 2011)\n")
cat("============================================================\n\n")

# Calculate distance matrix for both tests
dist_hell <- vegdist(comm_hell, method = "bray")

# --- PERMANOVA: Test for differences in multivariate location ---
cat("1. PERMANOVA: Testing site effect on composition (location)...\n")
set.seed(42)
perm_result <- adonis2(dist_hell ~ site,
                       data = coral_aligned,
                       permutations = 999)
site_r2 <- perm_result$R2[1]
site_f <- perm_result$F[1]
site_p <- perm_result$`Pr(>F)`[1]

cat("   Site effect:\n")
cat("     R² =", round(site_r2, 3), "\n")
cat("     F =", round(site_f, 2), "\n")
cat("     p =", site_p, "\n\n")

# --- PERMDISP: Test for homogeneity of multivariate dispersions ---
# This is essential per Anderson (2005, 2011) to distinguish location vs dispersion effects
cat("2. PERMDISP (betadisper): Testing homogeneity of dispersions...\n")
disp_result <- betadisper(dist_hell, coral_aligned$site)
disp_test <- permutest(disp_result, permutations = 999)

# Extract F and p from permutest output (stored in tab component)
disp_f <- disp_test$tab$F[1]
disp_p <- disp_test$tab$`Pr(>F)`[1]

cat("   Dispersion test:\n")
cat("     F =", round(disp_f, 2), "\n")
cat("     p =", round(disp_p, 3), "\n\n")

# Extract distances to centroid for each site
disp_distances <- data.frame(
  site = coral_aligned$site,
  distance_to_centroid = disp_result$distances
)
site_dispersions <- disp_distances %>%
  group_by(site) %>%
  summarise(
    mean_dispersion = mean(distance_to_centroid),
    sd_dispersion = sd(distance_to_centroid),
    n = n(),
    .groups = "drop"
  )

cat("   Average distance to centroid by site:\n")
for (i in 1:nrow(site_dispersions)) {
  cat(sprintf("     %s: %.3f ± %.3f (n=%d)\n",
              site_dispersions$site[i],
              site_dispersions$mean_dispersion[i],
              site_dispersions$sd_dispersion[i],
              site_dispersions$n[i]))
}
cat("\n")

# --- Interpretation ---
cat("3. INTERPRETATION (following Anderson et al. 2011):\n")
if (site_p < 0.05 && disp_p >= 0.05) {
  interpretation <- "Location effect confirmed"
  cat("   PERMANOVA significant, PERMDISP not significant:\n")
  cat("   → Sites differ in COMPOSITION (centroid location)\n")
  cat("   → Within-site variation is homogeneous across sites\n")
  cat("   → VALID interpretation: reef-scale species pools structure composition\n")
} else if (site_p < 0.05 && disp_p < 0.05) {
  interpretation <- "Location + dispersion effects"
  cat("   BOTH PERMANOVA and PERMDISP significant:\n")
  cat("   → Sites differ in BOTH composition AND within-site variation\n")
  cat("   → Some sites have more variable communities than others\n")
  cat("   → Interpret with caution: PERMANOVA inflated by dispersion differences\n")
} else if (site_p >= 0.05 && disp_p < 0.05) {
  interpretation <- "Dispersion effect only"
  cat("   PERMANOVA not significant, PERMDISP significant:\n")
  cat("   → Sites do NOT differ in average composition\n")
  cat("   → But sites DO differ in within-site variability\n")
} else {
  interpretation <- "No significant effects"
  cat("   Neither test significant:\n")
  cat("   → No evidence for site differences in composition or dispersion\n")
}
cat("\n")

# --- Pairwise PERMANOVA (if overall is significant) ---
if (site_p < 0.05) {
  cat("4. PAIRWISE PERMANOVA (site comparisons):\n")
  sites <- unique(coral_aligned$site)
  pairwise_results <- data.frame()

  # Convert dist object to matrix for subsetting
  dist_matrix <- as.matrix(dist_hell)

  for (i in 1:(length(sites) - 1)) {
    for (j in (i + 1):length(sites)) {
      site_pair <- c(sites[i], sites[j])
      idx <- coral_aligned$site %in% site_pair

      # Subset distance matrix and convert back to dist
      dist_subset <- as.dist(dist_matrix[idx, idx])

      set.seed(42)
      pair_perm <- adonis2(dist_subset ~ site,
                           data = coral_aligned[idx, ],
                           permutations = 999)

      pairwise_results <- rbind(pairwise_results, data.frame(
        comparison = paste(site_pair, collapse = " vs "),
        R2 = pair_perm$R2[1],
        F = pair_perm$F[1],
        p = pair_perm$`Pr(>F)`[1]
      ))
    }
  }

  # Bonferroni correction for pairwise comparisons
  pairwise_results$p_adj <- p.adjust(pairwise_results$p, method = "bonferroni")

  for (k in 1:nrow(pairwise_results)) {
    sig_marker <- ifelse(pairwise_results$p_adj[k] < 0.05, "*", "")
    cat(sprintf("   %s: R² = %.3f, F = %.2f, p = %.3f (adj = %.3f)%s\n",
                pairwise_results$comparison[k],
                pairwise_results$R2[k],
                pairwise_results$F[k],
                pairwise_results$p[k],
                pairwise_results$p_adj[k],
                sig_marker))
  }
  cat("   (* = significant after Bonferroni correction)\n\n")
}

# Store results for later use
permanova_results <- list(
  permanova = perm_result,
  permdisp = disp_test,
  site_dispersions = site_dispersions,
  interpretation = interpretation,
  pairwise = if (site_p < 0.05) pairwise_results else NULL
)

# ============================================================================
# CREATE FIGURE
# ============================================================================

cat("Creating figure...\n")

# Build plot
fig3 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  # Corals as filled points, sized by abundance
  geom_point(aes(fill = site, size = log10(total_cafi)),
             shape = 21, alpha = 0.8, stroke = 0.4, color = "white") +
  # Site ellipses
  stat_ellipse(aes(color = site),
               level = 0.95, linetype = "solid", linewidth = 1.1) +
  # Species vectors
  geom_segment(data = top_vectors,
               aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               color = "gray25", linewidth = 0.7, alpha = 0.8) +
  # Species labels
  geom_text_repel(
    data = top_vectors,
    aes(x = NMDS1_scaled * 1.1, y = NMDS2_scaled * 1.1, label = species_label),
    size = 3, fontface = "italic", color = "gray15",
    max.overlaps = 20, segment.alpha = 0.4, box.padding = 0.4,
    point.padding = 0.2
  ) +
  # Scales
  scale_fill_manual(values = site_colors, name = "Reef Site") +
  scale_color_manual(values = site_colors, guide = "none") +
  scale_size_continuous(range = c(2, 6), name = "log₁₀(Abundance)",
                        breaks = c(1, 1.5, 2), labels = c("10", "30", "100")) +
  # Labels
  labs(
    title = "Figure 3: Reef-Scale Species Pools Structure CAFI Composition",
    subtitle = sprintf("n = %d corals | Stress = %.3f | PERMANOVA R² = %.3f, p = %.3f | PERMDISP p = %.3f",
                       nrow(nmds_scores), nmds$stress, site_r2, site_p, disp_p),
    x = "NMDS1",
    y = "NMDS2",
    caption = paste0(
      "Bray-Curtis on Hellinger-transformed abundances (Legendre & Gallagher 2001). ",
      "Ellipses = 95% CI. Vectors = species (envfit, p < 0.01, R² > 0.12).\n",
      "PERMANOVA tests location; PERMDISP tests dispersion homogeneity (Anderson et al. 2011). ",
      "Corals with <", min_cafi, " CAFI (n=", n_excluded, ") excluded."
    )
  ) +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    plot.caption = element_text(size = 8, color = "gray50", hjust = 0),
    legend.position = "right",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 4)),
         size = guide_legend(order = 2))

# Save
ggsave(file.path(PATHS$fig_manuscript, "fig3_composition.png"), fig3,
       width = 10, height = 8, dpi = 300, bg = "white")

cat("\nFigure saved to: output/figures/manuscript/fig3_composition.png\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("SUMMARY: Site Effects on CAFI Composition\n")
cat("============================================================\n\n")

cat("SAMPLE SIZES:\n")
cat("  Corals included:", nrow(nmds_scores), "\n")
cat("  Corals excluded (n <", min_cafi, "):", n_excluded, "\n")
cat("  Species:", ncol(comm_aligned), "\n\n")

cat("ORDINATION (NMDS):\n")
cat("  Stress:", round(nmds$stress, 3), "\n")
cat("  Distance: Bray-Curtis on Hellinger-transformed abundances\n")
cat("  Transformation rationale: Handles double-zero problem,\n")
cat("    emphasizes composition (Legendre & Gallagher 2001)\n\n")

cat("PERMANOVA (Anderson 2001):\n")
cat("  Site R² =", round(site_r2, 3), "\n")
cat("  F =", round(site_f, 2), "\n")
cat("  p =", site_p, "\n")
cat("  Interpretation: Tests differences in multivariate location (centroids)\n\n")

cat("PERMDISP/betadisper (Anderson 2006):\n")
cat("  F =", round(disp_f, 2), "\n")
cat("  p =", round(disp_p, 3), "\n")
cat("  Interpretation: Tests homogeneity of within-group dispersion\n\n")

cat("COMBINED INTERPRETATION (Anderson et al. 2011 Ecology Letters):\n")
cat(" ", interpretation, "\n")
if (site_p < 0.05 && disp_p >= 0.05) {
  cat("  → Site differences are TRUE compositional differences\n")
  cat("  → Not an artifact of varying within-site heterogeneity\n")
} else if (site_p < 0.05 && disp_p < 0.05) {
  cat("  → Caution: PERMANOVA may be inflated by dispersion differences\n")
  cat("  → Consider reporting both effects\n")
}
cat("\n")

cat("SITE DISPERSIONS (distance to centroid):\n")
print(site_dispersions)
cat("\n")

cat("TOP SPECIES DRIVING SITE SEPARATION (envfit):\n")
print(top_vectors %>% dplyr::select(species, r2, pval))

cat("\n============================================================\n")
cat("METHODOLOGICAL NOTES (for manuscript):\n")
cat("============================================================\n")
cat("Following Anderson et al. (2011) Ecology Letters, we used:\n")
cat("  1. Hellinger transformation to handle double-zero problem\n")
cat("  2. Bray-Curtis dissimilarity (excludes joint absences)\n")
cat("  3. PERMANOVA to test location differences\n")
cat("  4. PERMDISP to test dispersion homogeneity\n")
cat("  5. Both tests together to properly interpret site effects\n")
cat("\nReferences:\n")
cat("  - Legendre & Gallagher (2001) Oecologia 129:271-280\n")
cat("  - Anderson (2001) Austral Ecol 26:32-46\n")
cat("  - Anderson (2006) Biometrics 62:245-253\n")
cat("  - Anderson et al. (2011) Ecol Lett 14:19-28\n")
cat("============================================================\n")
cat("Done!\n")
cat("============================================================\n")
