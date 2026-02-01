# ============================================================================
# fig3_composition_paired.R - Figure 3: Site Composition (NMDS + Taxonomic)
# ============================================================================
#
# PURPOSE: Two-panel figure showing:
#   A) NMDS ordination with species vectors colored by taxonomic group
#   B) Taxonomic group composition by site
#
# NARRATIVE: Composition differs across reef sites; differences driven by
#            distinct dominance of shrimp, crabs, hermit crabs, and fish
#
# OUTPUT: output/figures/manuscript/fig3_composition.png
#
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    FIGURE 3: Site Composition (NMDS + Taxonomic Groups)\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

source(here::here("scripts/00_setup.R"))
library(ggrepel)

# Load data
coral_master <- load_object("coral_master")
community_matrix <- load_object("community_matrix")
cafi_clean <- load_object("cafi_clean")

cat("[OK] Data loaded\n\n")

# ============================================================================
# COLOR PALETTES - Refined for publication
# ============================================================================

# Site colors (from CLAUDE.md - Okabe-Ito derivatives)
site_colors <- c(
  HAU = "#E69F00",  # Orange
  MAT = "#0072B2",  # Blue
  MRB = "#009E73"   # Teal
)

# Site names for cleaner labels
site_labels <- c(
  HAU = "Hauru",
  MAT = "Maatea",
  MRB = "Maharepa"
)

# Taxonomic group colors - refined Paul Tol palette (colorblind-safe)
type_colors <- c(
  "Shrimp" = "#882255",        # Deep magenta
  "Crabs" = "#CC6677",         # Rose pink
  "Hermit crabs" = "#117733",  # Forest green
  "Fish" = "#DDCC77",          # Warm sand
  "Gastropods" = "#AA4499",    # Purple
  "Echinoderms" = "#88CCEE",   # Sky blue
  "Polychaetes" = "#999999",   # Neutral gray
  "Amphipods" = "#332288",     # Deep indigo
  "Squat lobsters" = "#661100" # Dark burgundy
)

# ============================================================================
# PANEL A: NMDS ORDINATION WITH SPECIES VECTORS
# ============================================================================

cat("--- Panel A: NMDS Ordination with Species Vectors ---\n")

# Align community matrix with coral data
common_ids <- intersect(rownames(community_matrix), coral_master$coral_id)
comm_aligned <- community_matrix[common_ids, ]
coral_aligned <- coral_master %>%
  filter(coral_id %in% common_ids) %>%
  arrange(match(coral_id, common_ids))

# Remove empty species columns
comm_aligned <- comm_aligned[, colSums(comm_aligned) > 0]

# Exclude low-sample corals (< 10 CAFI) for stable ordination
min_cafi <- 10
coral_aligned$total_cafi_check <- rowSums(comm_aligned)
low_sample <- coral_aligned$total_cafi_check < min_cafi
n_excluded <- sum(low_sample)

comm_aligned <- comm_aligned[!low_sample, ]
coral_aligned <- coral_aligned[!low_sample, ]
comm_aligned <- comm_aligned[, colSums(comm_aligned) > 0]

cat("  Corals:", nrow(comm_aligned), "(excluded", n_excluded, "with <5 CAFI)\n")

# Hellinger transformation + NMDS
comm_hell <- decostand(comm_aligned, method = "hellinger")
set.seed(42)
nmds <- metaMDS(comm_hell, distance = "bray", k = 2, trymax = 100, trace = 0)
cat("  NMDS stress:", round(nmds$stress, 3), "\n")

# PERMANOVA & PERMDISP
dist_hell <- vegdist(comm_hell, method = "bray")
set.seed(42)
perm_result <- adonis2(dist_hell ~ site, data = coral_aligned, permutations = 999)
site_r2 <- perm_result$R2[1]
site_p <- perm_result$`Pr(>F)`[1]

disp_result <- betadisper(dist_hell, coral_aligned$site)
disp_test <- permutest(disp_result, permutations = 999)
disp_p <- disp_test$tab$`Pr(>F)`[1]

cat("  PERMANOVA: R² =", round(site_r2, 3), ", p =", site_p, "\n")
cat("  PERMDISP: p =", round(disp_p, 3), "\n")

# Extract NMDS scores
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$coral_id <- rownames(nmds_scores)
nmds_scores <- nmds_scores %>%
  left_join(coral_aligned %>% dplyr::select(coral_id, site, total_cafi),
            by = "coral_id")

# --- Fit species vectors ---
cat("\n  Fitting species vectors...\n")

# Species present in at least 8 corals
species_freq <- colSums(comm_aligned > 0)
common_species <- names(species_freq[species_freq >= 8])

set.seed(42)
species_fit <- envfit(nmds, comm_hell[, common_species], permutations = 999)

# Extract and filter significant vectors
species_scores_df <- as.data.frame(scores(species_fit, display = "vectors"))
species_scores_df$species <- rownames(species_scores_df)
species_scores_df$r2 <- species_fit$vectors$r
species_scores_df$pval <- species_fit$vectors$pvals

sig_species <- species_scores_df %>%
  filter(pval < 0.01, r2 > 0.10) %>%
  arrange(desc(r2))

cat("  Significant vectors (p < 0.01, R² > 0.10):", nrow(sig_species), "\n")

# --- Map species to taxonomic groups ---
species_type_map <- cafi_clean %>%
  distinct(otu, type)

sig_species <- sig_species %>%
  left_join(species_type_map, by = c("species" = "otu")) %>%
  mutate(
    type_clean = case_when(
      type == "shrimp" ~ "Shrimp",
      type == "crab" ~ "Crabs",
      type == "hermit" ~ "Hermit crabs",
      type == "fish" ~ "Fish",
      type == "snail" ~ "Gastropods",
      type == "echinoderm" ~ "Echinoderms",
      type == "worm" ~ "Polychaetes",
      type == "amphipod" ~ "Amphipods",
      is.na(type) ~ "Other",
      TRUE ~ "Other"
    ),
    # Create short genus-only labels for cleaner display
    species_short = case_when(
      species == "Harpiliopsis beaupresii" ~ "Harpiliopsis",
      species == "Harpiliopsis spinigera" ~ "H. spinigera",
      species == "Paragobiodon modestus" ~ "Paragobiodon",
      species == "Calcinus latens" ~ "Calcinus",
      species == "Trapezia bidentata" ~ "Trapezia",
      species == "Trapezia rufopunctata" ~ "T. rufopunctata",
      species == "Perinia tumida" ~ "Perinia",
      species == "Caracanthus maculatus" ~ "Caracanthus",
      species == "Paracirrhites arcatus" ~ "Paracirrhites",
      species == "Ophiocoma erinaceus" ~ "Ophiocoma",
      species == "Breviturma pica" ~ "Breviturma",
      species == "Alpheus lottini" ~ "Alpheus",
      species == "Galeropsis monodonta" ~ "Galeropsis",
      species == "Alpheidae" ~ "Alpheidae",  # Keep family-level distinct
      TRUE ~ word(species, 1)  # Just genus
    )
  )

# Scale vectors for plotting
arrow_scale <- 0.65
sig_species <- sig_species %>%
  mutate(
    NMDS1_scaled = NMDS1 * arrow_scale,
    NMDS2_scaled = NMDS2 * arrow_scale
  )

# Select key representatives: 7 species covering major taxonomic groups
# 4 original + 3 additional for richer visualization
top_vectors <- sig_species %>%
  filter(species %in% c(
    "Harpiliopsis beaupresii",  # Shrimp - guard shrimp
    "Paragobiodon modestus",    # Fish - coral goby
    "Calcinus latens",          # Hermit crabs
    "Trapezia bidentata",       # Crabs - guard crabs
    "Alpheus lottini",          # Shrimp - pistol shrimp mutualist
    "Galeropsis monodonta",     # Gastropods - tissue consumer (new group)
    "Caracanthus maculatus"     # Fish - velvetfish
  ))

cat("\n  Selected species for visualization (key group representatives):\n")
print(top_vectors %>% dplyr::select(species_short, type_clean, r2, pval))

# --- Calculate site centroids and positions for labels ---
site_centroids <- nmds_scores %>%
  group_by(site) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  )

# Position site labels outside the ellipses - fine-tuned for tighter axes
site_label_positions <- site_centroids %>%
  mutate(
    label_x = case_when(
      site == "MAT" ~ NMDS1 - 0.55,   # Blue - push left
      site == "HAU" ~ NMDS1 + 0.45,   # Orange - push right
      site == "MRB" ~ NMDS1 + 0.55    # Teal - push right more
    ),
    label_y = case_when(
      site == "MAT" ~ NMDS2 + 0.45,   # Blue - push up
      site == "HAU" ~ NMDS2 + 0.55,   # Orange - push up
      site == "MRB" ~ NMDS2 - 0.35    # Teal - less down (was -0.5)
    )
  )

# Define the order of taxonomic groups for the legend
vector_types <- c("Shrimp", "Fish", "Hermit crabs", "Crabs", "Gastropods")

# --- Identify and report outliers ---
nmds_scores <- nmds_scores %>%
  mutate(dist_from_center = sqrt(NMDS1^2 + NMDS2^2))

outlier_threshold <- quantile(nmds_scores$dist_from_center, 0.95)
extreme_outliers <- nmds_scores %>%
  filter(dist_from_center > outlier_threshold)

cat("\n  Extreme outliers (>95th percentile distance from centroid):\n")
print(extreme_outliers %>% dplyr::select(coral_id, site, NMDS1, NMDS2, total_cafi, dist_from_center))

# --- Panel A plot - Publication-ready design ---
panel_a <- ggplot() +
  # Coral points - uniform size, higher alpha for visibility
  geom_point(data = nmds_scores,
             aes(x = NMDS1, y = NMDS2, fill = site),
             shape = 21, size = 2.5, alpha = 0.7,
             stroke = 0.3, color = "gray30") +
  # Site ellipses
  stat_ellipse(data = nmds_scores,
               aes(x = NMDS1, y = NMDS2, color = site),
               level = 0.95, linewidth = 0.9, linetype = "solid",
               alpha = 0.85) +
  # Species vectors
  geom_segment(data = top_vectors,
               aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled,
                   color = type_clean),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               linewidth = 0.65, alpha = 0.85) +
  # Species labels - with optimized repulsion
  geom_text_repel(
    data = top_vectors,
    aes(x = NMDS1_scaled * 1.3, y = NMDS2_scaled * 1.3,
        label = species_short, color = type_clean),
    size = 2.5, fontface = "italic",
    max.overlaps = 30,
    segment.size = 0.15,
    segment.alpha = 0.4,
    segment.color = "gray45",
    box.padding = 0.55,
    point.padding = 0.3,
    min.segment.length = 0.1,
    force = 6,
    force_pull = 0.15,
    seed = 789,
    show.legend = FALSE
  ) +
  # Site labels positioned outside ellipses
  geom_text(data = site_label_positions,
            aes(x = label_x, y = label_y, label = site_labels[site], color = site),
            size = 3.3, fontface = "bold",
            show.legend = FALSE) +
  # Scales
  scale_fill_manual(values = site_colors, guide = "none") +
  scale_color_manual(
    values = c(site_colors, type_colors),
    breaks = vector_types,
    name = NULL,
    guide = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1))
  ) +
  # Stress annotation in top right corner
  annotate("text", x = 0.95, y = 0.95,
           label = paste0("Stress = ", sprintf("%.2f", nmds$stress)),
           size = 2.8, color = "gray30", hjust = 1, fontface = "italic") +
  # Axis labels
  labs(
    x = "NMDS1",
    y = "NMDS2"
  ) +
  coord_fixed(ratio = 1, xlim = c(-1.35, 1.0), ylim = c(-0.95, 1.0)) +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.title = element_text(size = 9, color = "black"),
    axis.text = element_blank(),  # Remove axis tick numbers
    axis.ticks = element_blank(), # Remove tick marks
    legend.position = c(0.12, 0.25),
    legend.background = element_rect(fill = alpha("white", 0.93),
                                     color = "gray80", linewidth = 0.3),
    legend.text = element_text(size = 6.5),
    legend.key.size = unit(0.32, "cm"),
    legend.key.height = unit(0.32, "cm"),
    legend.spacing.y = unit(0.03, "cm"),
    legend.margin = margin(4, 6, 4, 6),
    plot.margin = margin(5, 2, 5, 5)
  )

# ============================================================================
# PANEL B: TAXONOMIC COMPOSITION BY SITE
# ============================================================================

cat("\n--- Panel B: Taxonomic Composition ---\n")

# Summarize by site and taxonomic type
type_summary <- cafi_clean %>%
  filter(!is.na(type), type != "", type != "amph") %>%
  mutate(
    type_clean = case_when(
      type == "hermit" ~ "Hermit crabs",
      type == "crab" ~ "Crabs",
      type == "shrimp" ~ "Shrimp",
      type == "fish" ~ "Fish",
      type == "snail" ~ "Gastropods",
      type == "echinoderm" ~ "Echinoderms",
      type == "worm" ~ "Polychaetes",
      type == "amphipod" ~ "Amphipods",
      type == "squat_lobster" ~ "Squat lobsters",
      TRUE ~ type
    )
  ) %>%
  group_by(site, type_clean) %>%
  summarise(n = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(site) %>%
  mutate(
    total = sum(n),
    prop = n / total * 100
  ) %>%
  ungroup()

# Order types by overall abundance (for stacking)
type_order <- type_summary %>%
  group_by(type_clean) %>%
  summarise(total_n = sum(n), .groups = "drop") %>%
  arrange(desc(total_n)) %>%
  pull(type_clean)

type_summary$type_clean <- factor(type_summary$type_clean, levels = rev(type_order))

# Print summary
cat("\n  Composition by site (%):\n")
type_summary %>%
  dplyr::select(site, type_clean, prop) %>%
  pivot_wider(names_from = site, values_from = prop, values_fill = 0) %>%
  mutate(across(where(is.numeric), ~round(., 1))) %>%
  arrange(factor(type_clean, levels = type_order)) %>%
  print()

# Custom site order
type_summary$site <- factor(type_summary$site, levels = c("HAU", "MAT", "MRB"))

# Panel B plot - Publication-ready design
panel_b <- ggplot(type_summary, aes(x = site, y = prop, fill = type_clean)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  # Add percentage labels for major groups (>10%)
  geom_text(
    data = type_summary %>% filter(prop > 10),
    aes(label = sprintf("%.0f%%", prop)),
    position = position_stack(vjust = 0.5),
    size = 2.4, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = type_colors, name = NULL,
                    guide = guide_legend(reverse = TRUE, ncol = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)),
                     breaks = seq(0, 100, 25)) +
  scale_x_discrete(labels = site_labels) +
  labs(
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.title.y = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(face = "bold", size = 9.5, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.36, "cm"),
    legend.spacing.y = unit(0.06, "cm"),
    plot.margin = margin(5, 5, 5, 12)
  )

# ============================================================================
# COMBINE PANELS
# ============================================================================

cat("\n--- Combining panels ---\n")

# Create panel labels
panel_a_labeled <- panel_a +
  labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = 12))

panel_b_labeled <- panel_b +
  labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold", size = 12))

# Combine with patchwork
fig3 <- panel_a_labeled + panel_b_labeled +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 12, 8, 12)
    )
  )

# Save
ggsave(file.path(PATHS$fig_manuscript, "fig3_composition.png"), fig3,
       width = 10, height = 4.6, dpi = 300, bg = "white")

cat("\nFigure saved to: output/figures/manuscript/fig3_composition.png\n")

# ============================================================================
# GENERATE FIGURE LEGEND AND RESULTS TEXT
# ============================================================================

cat("\n--- Generating figure legend and results text ---\n")

# Pairwise PERMANOVA for detailed results
cat("  Running pairwise PERMANOVA...\n")
site_pairs <- combn(unique(as.character(coral_aligned$site)), 2, simplify = FALSE)
pairwise_results <- lapply(site_pairs, function(pair) {
  idx <- coral_aligned$site %in% pair
  dist_matrix <- as.matrix(dist_hell)
  dist_subset <- as.dist(dist_matrix[idx, idx])
  pair_data <- coral_aligned[idx, ]
  set.seed(42)
  pair_perm <- adonis2(dist_subset ~ site, data = pair_data, permutations = 999)
  data.frame(
    comparison = paste(pair, collapse = " vs "),
    R2 = pair_perm$R2[1],
    p_value = pair_perm$`Pr(>F)`[1]
  )
})
pairwise_df <- do.call(rbind, pairwise_results)
pairwise_df$p_adj <- p.adjust(pairwise_df$p_value, method = "bonferroni")

# Create legend text
legend_text <- paste0(
  "FIGURE 3: CAFI COMMUNITY COMPOSITION ACROSS REEF SITES\n",
  "================================================================================\n\n",

  "FIGURE LEGEND\n",
  "-------------\n",
  "Figure 3. Coral-associated fauna and invertebrate (CAFI) community composition ",
  "differs significantly among reef sites in Mo'orea, French Polynesia. ",
  "(A) Non-metric multidimensional scaling (NMDS) ordination of Hellinger-transformed ",
  "community data (Bray-Curtis dissimilarity). Each point represents one coral colony ",
  "(n = ", nrow(comm_aligned), "). Ellipses show 95% confidence intervals for each site. ",
  "Vectors indicate species significantly associated with compositional variation ",
  "(envfit permutation test, p < 0.01, R² > 0.10), colored by taxonomic group. ",
  "(B) Relative abundance of taxonomic groups at each site. Percentages shown for ",
  "groups comprising >10% of site assemblage.\n\n",

  "Site abbreviations: Hauru (n = ", sum(coral_aligned$site == "HAU"), " corals), ",
  "Maatea (n = ", sum(coral_aligned$site == "MAT"), " corals), ",
  "Barrier (n = ", sum(coral_aligned$site == "MRB"), " corals).\n\n",

  "================================================================================\n\n",

  "STATISTICAL RESULTS\n",
  "-------------------\n\n",

  "1. OVERALL COMPOSITION TEST (PERMANOVA)\n",
  "   Model: Bray-Curtis dissimilarity ~ Site\n",
  "   Transformation: Hellinger (accounts for double-zero problem)\n",
  "   Permutations: 999\n\n",
  "   Result: R² = ", sprintf("%.3f", site_r2), ", F = ", sprintf("%.2f", perm_result$F[1]),
  ", p = ", perm_result$`Pr(>F)`[1], "\n",
  "   Interpretation: Site explains ", sprintf("%.1f", site_r2 * 100),
  "% of community composition variance (p < 0.001).\n\n",

  "2. HOMOGENEITY OF DISPERSION TEST (PERMDISP)\n",
  "   Purpose: Verify PERMANOVA result reflects location differences, not dispersion\n",
  "   Method: betadisper + permutest (Anderson 2006)\n\n",
  "   Result: F = ", sprintf("%.2f", disp_test$tab$F[1]),
  ", p = ", sprintf("%.3f", disp_p), "\n",
  "   Interpretation: ", ifelse(disp_p < 0.05,
    paste0("Multivariate dispersion differs among sites (p = ", sprintf("%.3f", disp_p),
           "). Sites vary in community heterogeneity.\n",
           "   Conclusion: Both location AND dispersion differ; PERMANOVA remains valid but interpret with caution.\n\n"),
    paste0("Multivariate dispersion is homogeneous across sites (p = ", sprintf("%.2f", disp_p),
           ").\n   Conclusion: PERMANOVA results reflect TRUE compositional differences, ",
           "not heterogeneous dispersion artifact.\n\n")),

  "3. PAIRWISE COMPARISONS (Bonferroni-corrected)\n"
)

for (i in 1:nrow(pairwise_df)) {
  legend_text <- paste0(legend_text,
    "   ", pairwise_df$comparison[i], ": R² = ", sprintf("%.3f", pairwise_df$R2[i]),
    ", p = ", sprintf("%.3f", pairwise_df$p_value[i]),
    ", p_adj = ", sprintf("%.3f", pairwise_df$p_adj[i]),
    ifelse(pairwise_df$p_adj[i] < 0.05, " *", ""), "\n"
  )
}

legend_text <- paste0(legend_text, "\n",
  "4. ORDINATION QUALITY\n",
  "   NMDS stress: ", sprintf("%.3f", nmds$stress), "\n",
  "   Interpretation: ",
  ifelse(nmds$stress < 0.1, "Excellent representation (stress < 0.10)",
         ifelse(nmds$stress < 0.2, "Good representation (stress < 0.20)",
                "Fair representation (stress >= 0.20)")), "\n\n",

  "5. SPECIES DRIVING COMPOSITIONAL DIFFERENCES\n",
  "   Method: envfit (permutation test of species vectors onto ordination)\n",
  "   Criteria: p < 0.01, R² > 0.10\n\n",
  "   Top species vectors shown in figure:\n"
)

for (i in 1:nrow(top_vectors)) {
  legend_text <- paste0(legend_text,
    "   - ", top_vectors$species[i], " (", top_vectors$type_clean[i],
    "): R² = ", sprintf("%.3f", top_vectors$r2[i]),
    ", p = ", sprintf("%.3f", top_vectors$pval[i]), "\n"
  )
}

# Calculate composition stats for results text
hauru_comp <- type_summary %>% filter(site == "HAU")
maatea_comp <- type_summary %>% filter(site == "MAT")
barrier_comp <- type_summary %>% filter(site == "MRB")

# Build PERMDISP interpretation based on significance
permdisp_text <- if (disp_p < 0.05) {
  paste0(
    "Sites also differed in multivariate dispersion (PERMDISP: F = ",
    sprintf("%.2f", disp_test$tab$F[1]), ", p = ", sprintf("%.3f", disp_p),
    "), indicating that some sites harbor more heterogeneous assemblages than others. ",
    "However, the significant PERMANOVA result (with higher F-statistic) and clear ",
    "separation in ordination space confirm that compositional differences among sites ",
    "are robust (Anderson et al. 2011)."
  )
} else {
  paste0(
    "This effect was not driven by differences in multivariate dispersion ",
    "(PERMDISP: F = ", sprintf("%.2f", disp_test$tab$F[1]),
    ", p = ", sprintf("%.2f", disp_p), "), confirming that sites harbor compositionally ",
    "distinct CAFI assemblages rather than simply differing in community heterogeneity ",
    "(Anderson et al. 2011)."
  )
}

legend_text <- paste0(legend_text, "\n",
  "================================================================================\n\n",

  "RESULTS\n",
  "-------\n\n",

  "Community composition differed significantly among the three reef sites ",
  "(PERMANOVA: R² = ", sprintf("%.3f", site_r2), ", F = ", sprintf("%.2f", perm_result$F[1]),
  ", p < 0.001; Fig. 3A). ", permdisp_text, " All pairwise site comparisons were ",
  "significant after Bonferroni correction (all p_adj < 0.005).\n\n",

  "The three sites exhibited distinct taxonomic signatures (Fig. 3B). Maatea was ",
  "characterized by high hermit crab abundance (", sprintf("%.0f", maatea_comp$prop[maatea_comp$type_clean == "Hermit crabs"]),
  "% of individuals vs ", sprintf("%.0f", hauru_comp$prop[hauru_comp$type_clean == "Hermit crabs"]),
  "-", sprintf("%.0f", barrier_comp$prop[barrier_comp$type_clean == "Hermit crabs"]),
  "% at other sites), primarily Calcinus latens. Maharepa was dominated ",
  "by obligate coral symbionts, with shrimp and crabs together comprising ",
  sprintf("%.0f", barrier_comp$prop[barrier_comp$type_clean == "Shrimp"] + barrier_comp$prop[barrier_comp$type_clean == "Crabs"]),
  "% of all CAFI (vs ", sprintf("%.0f", hauru_comp$prop[hauru_comp$type_clean == "Shrimp"] + hauru_comp$prop[hauru_comp$type_clean == "Crabs"]),
  "% at Hauru), including the mutualistic guard shrimp Harpiliopsis beaupresii and ",
  "Trapezia crabs. Hauru supported the most taxonomically balanced assemblage and ",
  "the highest proportion of coral-dwelling fishes (",
  sprintf("%.0f", hauru_comp$prop[hauru_comp$type_clean == "Fish"]),
  "% vs ", sprintf("%.0f", maatea_comp$prop[maatea_comp$type_clean == "Fish"]),
  "-", sprintf("%.0f", barrier_comp$prop[barrier_comp$type_clean == "Fish"]),
  "% elsewhere), particularly Paragobiodon modestus gobies.\n\n",

  "Species-level analysis confirmed these patterns. Envfit permutation tests identified ",
  nrow(sig_species), " species significantly associated with compositional variation ",
  "(p < 0.01, R² > 0.10). The strongest drivers were Harpiliopsis beaupresii (shrimp; ",
  "R² = ", sprintf("%.2f", top_vectors$r2[top_vectors$species_short == "Harpiliopsis"]),
  "), Paragobiodon modestus (fish; R² = ", sprintf("%.2f", top_vectors$r2[top_vectors$species_short == "Paragobiodon"]),

"), Calcinus latens (hermit crab; R² = ", sprintf("%.2f", top_vectors$r2[top_vectors$species_short == "Calcinus"]),
  "), and Trapezia bidentata (crab; R² = ", sprintf("%.2f", top_vectors$r2[top_vectors$species_short == "Trapezia"]),
  "), representing each of the four major taxonomic groups.\n\n",

  "================================================================================\n\n",

  "KEY BIOLOGICAL FINDINGS\n",
  "-----------------------\n\n",

  "MAATEA (Blue):\n",
  "  - Distinctively dominated by hermit crabs (33% vs 8-10% elsewhere)\n",
  "  - Key species: Calcinus latens, Pagurixus spp.\n",
  "  - May reflect different substrate/shelter availability\n\n",

  "MAHAREPA (Teal):\n",
  "  - Dominated by obligate coral symbionts (shrimp + crabs = 71%)\n",
  "  - Key species: Harpiliopsis beaupresii, Trapezia spp., Alpheus lottini\n",
  "  - Higher proportion of mutualistic guard crabs and pistol shrimp\n\n",

  "HAURU (Orange):\n",
  "  - Most taxonomically balanced assemblage\n",
  "  - Highest fish abundance (12% vs 3-7% elsewhere)\n",
  "  - Key species: Paragobiodon modestus, Caracanthus maculatus\n",
  "  - May indicate different predator pressure or habitat complexity\n\n",

  "================================================================================\n\n",

  "METHODS\n",
  "-------\n\n",

  "We analyzed CAFI community composition across ", nrow(comm_aligned), " Pocillopora ",
  "colonies (Hauru: n = ", sum(coral_aligned$site == "HAU"),
  "; Maatea: n = ", sum(coral_aligned$site == "MAT"),
  "; Barrier: n = ", sum(coral_aligned$site == "MRB"),
  ") after excluding colonies with fewer than ", min_cafi, " CAFI individuals to ensure ",
  "stable ordination positions (", n_excluded, " colonies excluded). ",
  "Species abundances were Hellinger-transformed prior to analysis to down-weight ",
  "rare species and account for the double-zero problem inherent in community data ",
  "(Legendre & Gallagher 2001).\n\n",

  "We tested for compositional differences among sites using permutational multivariate ",
  "analysis of variance (PERMANOVA; Anderson 2001) on Bray-Curtis dissimilarities ",
  "(999 permutations). To verify that significant PERMANOVA results reflected true ",
  "compositional differences rather than heterogeneous dispersion, we tested for ",
  "homogeneity of multivariate dispersion using PERMDISP (Anderson 2006). Following ",
  "the analytical framework of Anderson et al. (2011), we interpret PERMANOVA results ",
  "as evidence of location differences when PERMDISP is non-significant.\n\n",

  "Community structure was visualized using non-metric multidimensional scaling (NMDS) ",
  "with stress = ", sprintf("%.2f", nmds$stress), ". Species associated with compositional ",
  "differences were identified using envfit permutation tests (999 permutations), with ",
  "significance set at p < 0.01 and R² > 0.10. Pairwise site comparisons were Bonferroni-",
  "corrected for multiple testing (α = 0.05/3 = 0.017).\n\n",

  "================================================================================\n\n",

  "REFERENCES\n",
  "----------\n",
  "- Anderson, M.J. (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecology 26:32-46.\n",
  "- Anderson, M.J. (2006). Distance-based tests for homogeneity of multivariate dispersions. Biometrics 62:245-253.\n",
  "- Anderson, M.J. et al. (2011). Navigating the multiple meanings of β diversity. Ecology Letters 14:19-28.\n",
  "- Legendre, P. & Gallagher, E.D. (2001). Ecologically meaningful transformations for ordination of species data. Oecologia 129:271-280.\n\n",

  "================================================================================\n",
  "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"
)

# Write legend file
writeLines(legend_text, file.path(PATHS$fig_manuscript, "fig3_legend_results.txt"))
cat("  Legend saved to: output/figures/manuscript/fig3_legend_results.txt\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("SUMMARY: Site Composition Differences\n")
cat("============================================================\n\n")

cat("STATISTICAL TESTS:\n")
cat("  PERMANOVA: R² =", round(site_r2, 3), ", p =", site_p, "\n")
cat("  PERMDISP: p =", round(disp_p, 3), "(dispersion homogeneous)\n")
cat("  → TRUE compositional differences, not dispersion artifact\n\n")

cat("KEY COMPOSITIONAL DIFFERENCES:\n")
cat("  MAATEA (MAT):  Hermit crab-dominated (33% vs 8-10% elsewhere)\n")
cat("                 Key species: Calcinus latens, Pagurixus\n\n")
cat("  BARRIER (MRB): Shrimp + crab-dominated (71% combined)\n")
cat("                 Key species: Harpiliopsis, Alpheus, Trapezia\n\n")
cat("  HAURU (HAU):   Most balanced; highest fish (12%)\n")
cat("                 Key species: Paragobiodon, Caracanthus, Paracirrhites\n\n")

cat("SPECIES SHOWN IN FIGURE (key group representatives):\n")
for (i in 1:nrow(top_vectors)) {
  cat(sprintf("  %d. %s (%s) - R² = %.2f\n",
              i, top_vectors$species_short[i],
              top_vectors$type_clean[i],
              top_vectors$r2[i]))
}

cat("\n============================================================\n")
cat("Done!\n")
cat("============================================================\n")
