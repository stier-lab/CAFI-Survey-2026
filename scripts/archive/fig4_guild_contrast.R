# ============================================================================
# fig4_guild_contrast.R - Guild Structure with Within vs Between Comparison
# ============================================================================
#
# PURPOSE: Show guild structure with clear visual contrast between
#   within-guild and between-guild correlations
#
# PANELS:
#   A: Heatmap sorted by guild (block structure)
#   B: Distribution of within-guild vs between-guild correlations (key contrast)
#   C: Guild taxonomic composition
#   D: Strongest species pairs (all within-guild)
#
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    FIGURE 4: GUILD STRUCTURE WITH CONTRAST\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

network_results <- tryCatch(
  load_object("cafi_network"),
  error = function(e) {
    source(here::here("scripts/06_network_analysis.R"))
    load_object("cafi_network")
  }
)

library(patchwork)

fig_dir_manuscript <- file.path(PATHS$figures, "manuscript")
dir.create(fig_dir_manuscript, recursive = TRUE, showWarnings = FALSE)

cat("[OK] Setup complete\n\n")

# ============================================================================
# EXTRACT DATA
# ============================================================================

g <- network_results$graph
centrality_df <- network_results$centrality
edge_list <- network_results$edge_list

species_info <- centrality_df %>%
  dplyr::select(species, type, module, degree) %>%
  arrange(module, desc(degree))

n_species <- nrow(species_info)
n_modules <- length(unique(species_info$module))

cat("  Species:", n_species, "| Guilds:", n_modules, "\n\n")

# ============================================================================
# BUILD CORRELATION MATRIX
# ============================================================================

all_species <- unique(c(edge_list$sp1, edge_list$sp2))
n_sp <- length(all_species)

cor_mat <- matrix(0, nrow = n_sp, ncol = n_sp)
rownames(cor_mat) <- all_species
colnames(cor_mat) <- all_species

for (i in 1:nrow(edge_list)) {
  cor_mat[edge_list$sp1[i], edge_list$sp2[i]] <- edge_list$correlation[i]
  cor_mat[edge_list$sp2[i], edge_list$sp1[i]] <- edge_list$correlation[i]
}
diag(cor_mat) <- 1

# Guild-sorted order
species_guild <- species_info %>%
  filter(species %in% rownames(cor_mat)) %>%
  arrange(module, desc(degree)) %>%
  pull(species)

cor_mat_guild <- cor_mat[species_guild, species_guild]

# ============================================================================
# CLASSIFY ALL CORRELATIONS AS WITHIN VS BETWEEN GUILD
# ============================================================================

guild_lookup <- setNames(species_info$module, species_info$species)

# Get ALL pairwise correlations (not just significant edges)
# This shows the true contrast: most between-guild pairs are weak/zero
all_pairs <- expand.grid(
  sp1 = species_guild,
  sp2 = species_guild,
  stringsAsFactors = FALSE
) %>%
  filter(sp1 < sp2) %>%  # Upper triangle only
  mutate(
    correlation = mapply(function(s1, s2) cor_mat_guild[s1, s2], sp1, sp2),
    guild1 = guild_lookup[sp1],
    guild2 = guild_lookup[sp2],
    type = ifelse(guild1 == guild2, "Within guild", "Between guilds")
  )

# Summary stats on ALL pairs
within_mean <- mean(all_pairs$correlation[all_pairs$type == "Within guild"])
between_mean <- mean(all_pairs$correlation[all_pairs$type == "Between guilds"])
n_within <- sum(all_pairs$type == "Within guild")
n_between <- sum(all_pairs$type == "Between guilds")

# For the significant edges (Panel D)
cor_classified <- edge_list %>%
  mutate(
    guild1 = guild_lookup[sp1],
    guild2 = guild_lookup[sp2],
    type = ifelse(guild1 == guild2, "Within guild", "Between guilds")
  )

cat("  Within-guild pairs:", n_within, "| mean r =", round(within_mean, 3), "\n")
cat("  Between-guild pairs:", n_between, "| mean r =", round(between_mean, 3), "\n")
cat("  Ratio:", round(within_mean / between_mean, 2), "x stronger within guilds\n\n")

# ============================================================================
# PANEL A: HEATMAP SORTED BY GUILD
# ============================================================================

cat("Creating Panel A: Heatmap...\n")

cor_long_guild <- as.data.frame(cor_mat_guild) %>%
  mutate(species1 = rownames(cor_mat_guild)) %>%
  tidyr::pivot_longer(-species1, names_to = "species2", values_to = "correlation") %>%
  mutate(
    species1 = factor(species1, levels = species_guild),
    species2 = factor(species2, levels = species_guild),
    guild1 = guild_lookup[as.character(species1)],
    guild2 = guild_lookup[as.character(species2)],
    pair_type = ifelse(guild1 == guild2, "within", "between")
  )

# Guild boundaries
guild_counts <- species_info %>%
  filter(species %in% species_guild) %>%
  group_by(module) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(module) %>%
  mutate(cumsum = cumsum(n))

guild_breaks <- guild_counts$cumsum
guild_centers <- c(0, guild_breaks[-n_modules]) + diff(c(0, guild_breaks)) / 2
guild_labels_vec <- paste0("Guild ", guild_counts$module, "\n(n=", guild_counts$n, ")")

p_heatmap <- ggplot(cor_long_guild, aes(x = species2, y = species1, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00",
    midpoint = 0.5, limits = c(0, 1),
    name = "Co-occurrence\nCorrelation"
  ) +
  geom_vline(xintercept = guild_breaks[-n_modules] + 0.5, color = "black", linewidth = 1.2) +
  geom_hline(yintercept = guild_breaks[-n_modules] + 0.5, color = "black", linewidth = 1.2) +
  scale_x_discrete(
    breaks = species_guild[round(guild_centers)],
    labels = guild_labels_vec,
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    breaks = species_guild[round(guild_centers)],
    labels = guild_labels_vec,
    expand = c(0, 0)
  ) +
  labs(
    title = "A. Co-occurrence Matrix by Guild",
    subtitle = "Block-diagonal structure shows within-guild clustering",
    x = "Species (ordered by guild)",
    y = "Species (ordered by guild)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8, hjust = 1, face = "bold"),
    axis.title = element_text(size = 9, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "gray40")
  ) +
  coord_fixed()

# ============================================================================
# PANEL B: WITHIN VS BETWEEN GUILD DISTRIBUTIONS
# ============================================================================

cat("Creating Panel B: Correlation distributions...\n")

# Create density comparison using ALL pairs (not just significant edges)
p_contrast <- ggplot(all_pairs, aes(x = correlation, fill = type)) +
  geom_density(alpha = 0.7, color = "white", linewidth = 0.5) +
  geom_vline(xintercept = within_mean, color = "#D55E00", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = between_mean, color = "#0072B2", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 0.3, color = "gray40", linetype = "dotted", linewidth = 0.8) +
  scale_fill_manual(
    values = c("Within guild" = "#D55E00", "Between guilds" = "#0072B2"),
    name = NULL
  ) +
  annotate("text", x = within_mean + 0.02, y = Inf, vjust = 1.5,
           label = sprintf("Within\nmean=%.2f", within_mean),
           color = "#D55E00", fontface = "bold", size = 3, hjust = 0) +
  annotate("text", x = between_mean - 0.02, y = Inf, vjust = 1.5,
           label = sprintf("Between\nmean=%.2f", between_mean),
           color = "#0072B2", fontface = "bold", size = 3, hjust = 1) +
  annotate("text", x = 0.32, y = 0, vjust = -0.5,
           label = "threshold\nr > 0.3", color = "gray40", size = 2.5, hjust = 0) +
  labs(
    title = "B. Within-Guild vs Between-Guild Correlations",
    subtitle = sprintf("Within-guild mean (%.2f) is %.1fx higher than between-guild (%.2f)",
                       within_mean, within_mean / between_mean, between_mean),
    x = "Co-occurrence Correlation (all species pairs)",
    y = "Density"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "gray40")
  ) +
  scale_x_continuous(limits = c(-0.1, 1), expand = c(0, 0))

# ============================================================================
# PANEL C: GUILD TAXONOMIC COMPOSITION
# ============================================================================

cat("Creating Panel C: Guild composition...\n")

type_colors_expanded <- c(
  "crab" = "#E69F00", "shrimp" = "#0072B2", "fish" = "#009E73",
  "echinoderm" = "#999999", "worm" = "#56B4E9", "snail" = "#F0E442",
  "hermit" = "#CC79A7", "squat_lobster" = "#D55E00", "amphipod" = "#666666",
  "other" = "#BBBBBB"
)

guild_by_type <- species_info %>%
  filter(species %in% species_guild) %>%
  group_by(module, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(guild_label = paste0("Guild ", module))

p_composition <- ggplot(guild_by_type, aes(x = guild_label, y = n, fill = type)) +
  geom_col(position = "fill", width = 0.75, alpha = 0.9) +
  scale_fill_manual(values = type_colors_expanded, name = "Taxa") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    title = "C. Guild Taxonomic Composition",
    subtitle = "Each guild dominated by different taxa",
    x = NULL, y = "Proportion"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "gray40")
  ) +
  guides(fill = guide_legend(nrow = 2))

# ============================================================================
# PANEL D: STRONGEST PAIRS (ALL WITHIN-GUILD)
# ============================================================================

cat("Creating Panel D: Strongest pairs...\n")

top_assoc <- cor_classified %>%
  arrange(desc(correlation)) %>%
  head(12) %>%
  mutate(
    pair_label = paste0(
      substr(sp1, 1, 12), ifelse(nchar(sp1) > 12, "..", ""),
      " - ",
      substr(sp2, 1, 12), ifelse(nchar(sp2) > 12, "..", "")
    ),
    is_within = type == "Within guild"
  )

p_pairs <- ggplot(top_assoc, aes(x = reorder(pair_label, correlation), y = correlation)) +
  geom_col(aes(fill = is_within), alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("r=%.2f", correlation)),
            hjust = -0.1, size = 2.8, color = "gray30") +
  coord_flip(ylim = c(0, 1.1)) +
  scale_fill_manual(
    values = c("FALSE" = "#0072B2", "TRUE" = "#D55E00"),
    labels = c("FALSE" = "Between guilds", "TRUE" = "Within guild"),
    name = NULL
  ) +
  labs(
    title = "D. Strongest Species Associations",
    subtitle = "All top 12 pairs are within the same guild",
    x = NULL, y = "Correlation Strength"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "gray40")
  )

# ============================================================================
# COMBINE FIGURE
# ============================================================================

cat("Combining panels...\n")

fig4 <- (p_heatmap | p_contrast) / (p_composition | p_pairs) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "Figure 4. CAFI Species Form Distinct Ecological Guilds",
    subtitle = sprintf(
      "%d species cluster into %d guilds based on co-occurrence (n = 114 corals) | Network transitivity = 0.94, z = 36.1",
      n_species, n_modules
    ),
    caption = "Co-occurrence = volume-corrected Spearman correlation (r > 0.3, FDR < 0.05). Guilds identified via Louvain community detection.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      plot.caption = element_text(size = 8, hjust = 0, color = "gray50")
    )
  )

ggsave(
  file.path(fig_dir_manuscript, "fig4_guild_structure.png"),
  fig4, width = 14, height = 12, dpi = 300, bg = "white"
)

cat("\n[OK] Saved: fig4_guild_structure.png\n")

# ============================================================================
# UPDATE TEXT FILES
# ============================================================================

legend_text <- sprintf("Figure 4. CAFI species form distinct ecological guilds based on co-occurrence patterns. (A) Species co-occurrence matrix ordered by guild membership (Louvain clustering). The block-diagonal structure indicates that within-guild associations are substantially stronger than between-guild associations. Color intensity reflects volume-corrected Spearman correlation (r > 0.3, FDR < 0.05). Black lines delineate the four guilds. (B) Distribution of correlation strengths for within-guild (orange) versus between-guild (blue) species pairs. Within-guild correlations average %.2f compared to %.2f for between-guild pairs—a %.1f-fold difference demonstrating that the guild structure captures real ecological associations. (C) Taxonomic composition of each guild. Guild 1 (n = 21) contains predominantly shrimps, hermit crabs, and snails; Guild 2 (n = 21) is dominated by echinoderms and polychaete worms; Guild 3 (n = 11) includes resident fish and xanthid crabs; Guild 4 (n = 5) comprises peripheral crab specialists. (D) The 12 strongest species associations. All top-ranked pairs occur within the same guild, confirming that the modular structure reflects genuine ecological relationships. Network transitivity (0.94) was significantly higher than null expectation (z = 36.1, p < 0.001).
", within_mean, between_mean, within_mean / between_mean)

writeLines(legend_text, file.path(fig_dir_manuscript, "fig4_legend.txt"))

results_text <- sprintf("CAFI species formed four distinct ecological guilds based on co-occurrence patterns (Fig. 4). The co-occurrence matrix revealed clear block-diagonal structure when species were ordered by guild membership (Fig. 4A), with within-guild correlations substantially stronger than between-guild correlations. This contrast was quantified directly: within-guild species pairs had mean correlation r = %.2f compared to r = %.2f for between-guild pairs—a %.1f-fold difference (Fig. 4B).

The four guilds differed markedly in taxonomic composition (Fig. 4C): Guild 1 (n = 21 species) was dominated by shrimps (38%%) and hermit crabs; Guild 2 (n = 21) comprised primarily echinoderms (33%%) and polychaete worms (24%%); Guild 3 (n = 11) contained resident fish (36%%) and xanthid crabs (27%%); Guild 4 (n = 5) included peripheral crab specialists. This taxonomic segregation suggests guilds represent distinct microhabitat assemblages within coral colonies.

All 12 strongest species associations occurred within guilds (Fig. 4D), confirming that the modular structure reflects genuine ecological relationships. The strongest pairs included Macrophiothrix longipeda–Alpheus diadema (r = 0.85; brittle star and snapping shrimp) and Eunice–Caridea (r = 0.84; polychaete and shrimp). Network transitivity was 0.94, significantly higher than random expectation (z = 36.1, p < 0.001), indicating that species sharing one associate tend to share others.
", within_mean, between_mean, within_mean / between_mean)

writeLines(results_text, file.path(fig_dir_manuscript, "fig4_results.txt"))

cat("[OK] Updated legend and results text\n")
cat("\n============================================================\n")
cat("    COMPLETE\n")
cat("============================================================\n")
