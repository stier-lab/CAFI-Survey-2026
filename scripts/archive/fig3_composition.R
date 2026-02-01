# ============================================================================
# fig3_composition.R - Figure 3: Community Composition Changes with Coral Size
# ============================================================================
#
# PURPOSE: Test and visualize how CAFI community composition changes with
#          coral size using multiple analytical approaches
#
# ANALYSES:
#   Panel A: NMDS ordination colored by coral size (continuous gradient)
#   Panel B: Nestedness vs Turnover decomposition (betapart)
#   Panel C: Rarefied species accumulation by size class (iNEXT)
#   Panel D: Indicator species for small vs large corals (indicspecies)
#
# KEY QUESTION: Do larger corals have DIFFERENT species (turnover) or just
#               MORE of the same species (nestedness)?
#
# ANSWER: After controlling for sampling intensity (rarefaction), there is
#         NO significant compositional divergence by coral size (p = 0.61).
#         The apparent "turnover" signal is an artifact of detecting more
#         rare species when sampling more individuals from larger corals.
#
# OUTPUTS:
#   - output/figures/manuscript/fig3_composition.png
#   - output/tables/beta_decomposition.csv
#   - output/tables/indicator_species_size.csv
#
# Author: CAFI Survey Analysis Pipeline
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    FIGURE 3: Community Composition Changes with Coral Size\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Set CRAN mirror for package installation
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Additional packages for this analysis
if (!require("betapart", quietly = TRUE)) install.packages("betapart", quiet = TRUE)
if (!require("iNEXT", quietly = TRUE)) install.packages("iNEXT", quiet = TRUE)
if (!require("indicspecies", quietly = TRUE)) install.packages("indicspecies", quiet = TRUE)

library(betapart)
library(iNEXT)
library(indicspecies)
library(cowplot)
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
cat("Corals with both community and size data:", length(common_ids), "\n")

# Filter to common IDs
comm_aligned <- community_matrix[common_ids, ]
coral_aligned <- coral_master %>%
  filter(coral_id %in% common_ids) %>%
  arrange(match(coral_id, common_ids))

# Remove empty species columns
comm_aligned <- comm_aligned[, colSums(comm_aligned) > 0]
cat("Species with at least 1 individual:", ncol(comm_aligned), "\n")

# Create size classes (tertiles)
coral_aligned <- coral_aligned %>%
  mutate(
    log_volume = log10(volume_field),
    size_class = cut(volume_field,
                     breaks = quantile(volume_field, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE)
  )

# Site colors (consistent with other figures)
site_colors <- c(HAU = "#E69F00", MAT = "#0072B2", MRB = "#009E73")

# Size class colors (warm gradient)
size_colors <- c(Small = "#fee090", Medium = "#fc8d59", Large = "#d73027")

cat("\nSize class distribution:\n")
print(table(coral_aligned$size_class))
cat("\n")

# ============================================================================
# PANEL A: NMDS ORDINATION BY CORAL SIZE
# ============================================================================
cat("============================================================\n")
cat("Panel A: NMDS Ordination\n")
cat("============================================================\n\n")

# Hellinger transformation and NMDS
comm_hell <- decostand(comm_aligned, method = "hellinger")

set.seed(42)
nmds <- metaMDS(comm_hell, distance = "bray", k = 2, trymax = 100, trace = 0)

cat("NMDS stress:", round(nmds$stress, 3), "\n")
if (nmds$stress < 0.1) {
  cat("  Interpretation: Excellent representation\n\n")
} else if (nmds$stress < 0.2) {
  cat("  Interpretation: Good representation\n\n")
} else {
  cat("  Interpretation: Fair representation (interpret with caution)\n\n")
}

# Extract scores
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$coral_id <- rownames(nmds_scores)
nmds_scores <- nmds_scores %>%
  left_join(coral_aligned %>% dplyr::select(coral_id, volume_field, log_volume, size_class, site),
            by = "coral_id")

# Fit species vectors to NMDS (for taxa driving site differences)
# Filter to species present in at least 10 corals (not super rare)
species_freq <- colSums(comm_aligned > 0)
common_species <- names(species_freq[species_freq >= 10])
cat("Species present in >=10 corals for vectors:", length(common_species), "\n")

# Fit species scores using envfit
set.seed(42)
species_fit <- envfit(nmds, comm_hell[, common_species], permutations = 999)

# Extract significant species vectors (p < 0.05)
species_scores_df <- as.data.frame(scores(species_fit, display = "vectors"))
species_scores_df$species <- rownames(species_scores_df)
species_scores_df$r2 <- species_fit$vectors$r
species_scores_df$pval <- species_fit$vectors$pvals

# Filter to significant and strong vectors
sig_species <- species_scores_df %>%
  filter(pval < 0.05, r2 > 0.15) %>%  # Significant and R² > 0.15
  arrange(desc(r2))

cat("Significant species vectors (p < 0.05, R² > 0.15):", nrow(sig_species), "\n")
if (nrow(sig_species) > 0) {
  print(sig_species %>% dplyr::select(species, r2, pval) %>% head(15))
}
cat("\n")

# Scale vectors for plotting (multiply by a factor for visibility)
arrow_scale <- 1.5
sig_species <- sig_species %>%
  mutate(
    NMDS1_scaled = NMDS1 * arrow_scale,
    NMDS2_scaled = NMDS2 * arrow_scale,
    # Clean species names for labels
    species_label = gsub("_", " ", species)
  )

# Take top 12 vectors by R²
top_vectors <- sig_species %>% slice_head(n = 12)

# Panel A: NMDS by SITE with species vectors
p_a <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  # Site points

geom_point(aes(fill = site), shape = 21, size = 3.5, alpha = 0.8, stroke = 0.4, color = "white") +
  # Site ellipses
  stat_ellipse(aes(color = site), level = 0.95, linetype = "solid", linewidth = 1) +
  # Species vectors
  geom_segment(data = top_vectors,
               aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "gray30", linewidth = 0.6, alpha = 0.7) +
  # Species labels
ggrepel::geom_text_repel(
    data = top_vectors,
    aes(x = NMDS1_scaled, y = NMDS2_scaled, label = species_label),
    size = 2.5, fontface = "italic", color = "gray20",
    max.overlaps = 15, segment.alpha = 0.3, box.padding = 0.3
  ) +
  scale_fill_manual(values = site_colors, name = "Reef Site") +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(
    tag = "A",
    title = "NMDS by Reef Site",
    subtitle = paste0("Stress = ", round(nmds$stress, 3),
                      "; vectors = species driving site separation (R² > 0.15, p < 0.05)"),
    x = "NMDS1",
    y = "NMDS2"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    legend.position = c(0.88, 0.15),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    panel.grid.minor = element_blank()
  )

# ============================================================================
# PANEL B: NESTEDNESS VS TURNOVER DECOMPOSITION
# ============================================================================
cat("============================================================\n")
cat("Panel B: Nestedness vs Turnover Decomposition\n")
cat("============================================================\n\n")

# Convert to presence-absence
comm_pa <- decostand(comm_aligned, method = "pa")

# Calculate beta diversity components (Sorensen family)
cat("Calculating beta diversity components...\n")
beta_pair <- beta.pair(comm_pa, index.family = "sorensen")

# Extract components
turnover_dist <- as.matrix(beta_pair$beta.sim)   # Pure species replacement
nestedness_dist <- as.matrix(beta_pair$beta.sne)  # Subset relationships
total_dist <- as.matrix(beta_pair$beta.sor)       # Total beta diversity

# Create size difference matrix
size_diff_matrix <- as.matrix(dist(coral_aligned$log_volume))

# Mantel tests: correlation between beta diversity and size difference
cat("\nMantel tests (999 permutations):\n")

set.seed(42)
mantel_turnover <- mantel(as.dist(turnover_dist), as.dist(size_diff_matrix),
                          method = "spearman", permutations = 999)
mantel_nestedness <- mantel(as.dist(nestedness_dist), as.dist(size_diff_matrix),
                            method = "spearman", permutations = 999)
mantel_total <- mantel(as.dist(total_dist), as.dist(size_diff_matrix),
                       method = "spearman", permutations = 999)

cat("  Turnover × Size:    r =", sprintf("%.3f", mantel_turnover$statistic),
    ", p =", sprintf("%.3f", mantel_turnover$signif), "\n")
cat("  Nestedness × Size:  r =", sprintf("%.3f", mantel_nestedness$statistic),
    ", p =", sprintf("%.3f", mantel_nestedness$signif), "\n")
cat("  Total Beta × Size:  r =", sprintf("%.3f", mantel_total$statistic),
    ", p =", sprintf("%.3f", mantel_total$signif), "\n\n")

# Interpretation
if (mantel_turnover$statistic > mantel_nestedness$statistic) {
  cat("  Interpretation: Turnover dominates - larger corals have DIFFERENT species\n\n")
} else {
  cat("  Interpretation: Nestedness dominates - larger corals are SUPERSETS of smaller corals\n\n")
}

# Save results
beta_results <- data.frame(
  component = c("Turnover", "Nestedness", "Total"),
  mantel_r = c(mantel_turnover$statistic, mantel_nestedness$statistic, mantel_total$statistic),
  p_value = c(mantel_turnover$signif, mantel_nestedness$signif, mantel_total$signif),
  significant = c(mantel_turnover$signif < 0.05, mantel_nestedness$signif < 0.05, mantel_total$signif < 0.05)
)
write_csv(beta_results, file.path(PATHS$tables, "beta_decomposition.csv"))

# Panel B: Bar chart of Mantel correlations
beta_plot_df <- beta_results %>%
  filter(component != "Total") %>%
  mutate(
    label = sprintf("r = %.2f\np = %.3f", mantel_r, p_value),
    component = factor(component, levels = c("Turnover", "Nestedness"))
  )

p_b <- ggplot(beta_plot_df, aes(x = component, y = mantel_r, fill = component)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
  geom_text(aes(label = label), vjust = -0.3, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Turnover" = "#0072B2", "Nestedness" = "#E69F00")) +
  scale_y_continuous(limits = c(min(0, min(beta_plot_df$mantel_r) - 0.05),
                                max(beta_plot_df$mantel_r) + 0.15),
                     expand = c(0, 0)) +
  labs(
    tag = "B",
    title = "Beta Diversity Decomposition",
    subtitle = "Correlation with coral size difference (Mantel test)",
    x = NULL,
    y = "Mantel r"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold", size = 10)
  )

# ============================================================================
# PANEL C: RAREFIED SPECIES ACCUMULATION BY SIZE CLASS
# ============================================================================
cat("============================================================\n")
cat("Panel C: Rarefied Species Accumulation\n")
cat("============================================================\n\n")

# Aggregate abundances by size class
comm_by_size <- lapply(levels(coral_aligned$size_class), function(sc) {
  idx <- which(coral_aligned$size_class == sc)
  colSums(comm_aligned[idx, ])
})
names(comm_by_size) <- levels(coral_aligned$size_class)

# Run iNEXT
cat("Running iNEXT rarefaction/extrapolation...\n")
set.seed(42)
inext_result <- iNEXT(comm_by_size, q = 0, datatype = "abundance",
                      nboot = 50, conf = 0.95)

cat("\nSpecies richness by size class:\n")
for (sc in names(comm_by_size)) {
  obs <- sum(comm_by_size[[sc]] > 0)
  n_ind <- sum(comm_by_size[[sc]])
  cat("  ", sc, ": ", obs, " species from ", n_ind, " individuals\n", sep = "")
}
cat("\n")

# Extract data for plotting
inext_df <- fortify(inext_result, type = 1)
inext_df$Assemblage <- factor(inext_df$Assemblage, levels = c("Small", "Medium", "Large"))

# Panel C: Rarefaction curves
p_c <- ggplot(inext_df, aes(x = x, y = y, color = Assemblage, fill = Assemblage)) +
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr), alpha = 0.2, color = NA) +
  geom_line(aes(linetype = Method), linewidth = 1) +
  geom_point(data = inext_df %>% filter(Method == "Observed"), size = 3, shape = 21, stroke = 1) +
  scale_color_manual(values = size_colors, name = "Size Class") +
  scale_fill_manual(values = size_colors, guide = "none") +
  scale_linetype_manual(values = c("Rarefaction" = "solid", "Extrapolation" = "dashed"),
                        guide = "none") +
  labs(
    tag = "C",
    title = "Species Accumulation",
    subtitle = "Rarefaction curves by size class (95% CI)",
    x = "Number of Individuals",
    y = "Species Richness"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    legend.position = c(0.75, 0.25),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    panel.grid.minor = element_blank()
  )

# ============================================================================
# PANEL D: INDICATOR SPECIES BY SIZE CLASS
# ============================================================================
cat("============================================================\n")
cat("Panel D: Indicator Species Analysis\n")
cat("============================================================\n\n")

# Filter to species present in at least 5 corals
comm_filtered <- comm_aligned[, colSums(comm_aligned > 0) >= 5]
cat("Species present in >=5 corals:", ncol(comm_filtered), "\n")

# Run indicator species analysis
set.seed(42)
indval_size <- multipatt(
  comm_filtered,
  coral_aligned$size_class,
  control = how(nperm = 999),
  duleg = TRUE  # Allow combination groups
)

# Extract significant indicators
indval_summary <- indval_size$sign
indval_summary$species <- rownames(indval_summary)
indval_summary <- indval_summary %>%
  filter(p.value < 0.05) %>%
  arrange(p.value)

cat("\nSignificant indicator species (p < 0.05):\n")
print(indval_summary %>% dplyr::select(species, stat, p.value))
cat("\n")

# Save full results
write_csv(indval_summary, file.path(PATHS$tables, "indicator_species_size.csv"))

# Identify which size class each species indicates
indval_summary <- indval_summary %>%
  mutate(
    indicates = case_when(
      s.Small == 1 & s.Medium == 0 & s.Large == 0 ~ "Small",
      s.Small == 0 & s.Medium == 1 & s.Large == 0 ~ "Medium",
      s.Small == 0 & s.Medium == 0 & s.Large == 1 ~ "Large",
      s.Small == 1 & s.Medium == 1 & s.Large == 0 ~ "Small+Medium",
      s.Small == 0 & s.Medium == 1 & s.Large == 1 ~ "Medium+Large",
      s.Small == 1 & s.Medium == 0 & s.Large == 1 ~ "Small+Large",
      TRUE ~ "All"
    )
  )

# Count indicators by size class
indicator_counts <- indval_summary %>%
  count(indicates) %>%
  mutate(indicates = factor(indicates, levels = c("Small", "Small+Medium", "Medium",
                                                   "Medium+Large", "Large", "All")))

cat("Indicator species counts:\n")
print(indicator_counts)
cat("\n")

# Panel D: Forest plot of top indicator species
# Take top 12 by indicator value
top_indicators <- indval_summary %>%
  slice_head(n = 12) %>%
  mutate(
    species_label = gsub("_", " ", species),
    species_label = factor(species_label, levels = rev(species_label)),
    fill_color = case_when(
      indicates == "Small" ~ size_colors["Small"],
      indicates == "Medium" ~ size_colors["Medium"],
      indicates == "Large" ~ size_colors["Large"],
      TRUE ~ "gray50"
    )
  )

p_d <- ggplot(top_indicators, aes(x = stat, y = species_label, fill = indicates)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("p=%.3f", p.value)), hjust = -0.1, size = 2.8) +
  scale_fill_manual(
    values = c("Small" = "#fee090",
               "Medium" = "#fc8d59",
               "Large" = "#d73027",
               "Small+Medium" = "#fdae61",
               "Medium+Large" = "#f46d43",
               "All" = "gray60"),
    name = "Indicates"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    tag = "D",
    title = "Indicator Species",
    subtitle = "Species associated with coral size classes",
    x = "Indicator Value (IndVal)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# ============================================================================
# COMBINE PANELS
# ============================================================================
cat("============================================================\n")
cat("Combining panels into Figure 3\n")
cat("============================================================\n\n")

# Build caption
caption_text <- paste0(
  "(A) NMDS ordination colored by reef site; vectors show species driving site separation (envfit, R² > 0.15, p < 0.05).\n",
  "(B) Mantel correlation of turnover (species replacement) and nestedness (subset structure) with coral size difference.\n",
  "(C) Sample-based rarefaction curves showing species accumulation by size class.\n",
  "(D) Significant indicator species (IndVal, p < 0.05) associated with coral size classes.\n",
  "Note: Coral size affects abundance but not composition after rarefaction (p = 0.61). Site-level species pools drive compositional differences."
)

# Combine using patchwork
fig3_combined <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(heights = c(1, 1))

# Add overall title and caption with cowplot
fig3 <- cowplot::ggdraw() +
  cowplot::draw_plot(fig3_combined, x = 0, y = 0.07, width = 1, height = 0.87) +
  cowplot::draw_label(
    "Figure 3: Community Composition Changes with Coral Size",
    x = 0.02, y = 0.98, hjust = 0, vjust = 1,
    fontface = "bold", size = 13
  ) +
  cowplot::draw_label(
    "Reef-scale species pools drive composition; coral size affects abundance but not community identity",
    x = 0.02, y = 0.95, hjust = 0, vjust = 1,
    color = "gray30", size = 9
  ) +
  cowplot::draw_label(
    caption_text,
    x = 0.02, y = 0.01, hjust = 0, vjust = 0,
    color = "gray45", size = 7.5, lineheight = 1.2
  )

# Save
ggsave(file.path(PATHS$fig_manuscript, "fig3_composition.png"), fig3,
       width = 13, height = 11, dpi = 300, bg = "white")

cat("Figure saved to: output/figures/manuscript/fig3_composition.png\n")
cat("Tables saved to: output/tables/beta_decomposition.csv\n")
cat("                 output/tables/indicator_species_size.csv\n\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================
cat("============================================================\n")
cat("SUMMARY OF RESULTS\n")
cat("============================================================\n\n")

cat("NMDS:\n")
cat("  Stress:", round(nmds$stress, 3), "\n\n")

cat("Beta Diversity Decomposition:\n")
cat("  Turnover × Size:   r =", sprintf("%.3f", mantel_turnover$statistic),
    ", p =", sprintf("%.3f", mantel_turnover$signif), "\n")
cat("  Nestedness × Size: r =", sprintf("%.3f", mantel_nestedness$statistic),
    ", p =", sprintf("%.3f", mantel_nestedness$signif), "\n")

cat("\n  CAUTION: Presence-absence turnover may be inflated by sampling intensity.\n")
cat("  Larger corals have more individuals → more rare species detected.\n")
cat("  After rarefaction (equalizing sample sizes), compositional divergence\n")
cat("  is NOT significant (p = 0.61) - see script 02_community_analysis.R\n")

cat("\nIndicator Species:\n")
cat("  Total significant (p < 0.05):", nrow(indval_summary), "\n")
if (nrow(indval_summary) > 0) {
  cat("  By size class:\n")
  print(indicator_counts)
}

cat("\n============================================================\n")
cat("Done!\n")
cat("============================================================\n")
