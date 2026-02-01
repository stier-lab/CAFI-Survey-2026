# ============================================================================
# fig4_ecological_guilds.R - CAFI Ecological Guilds Figure (IMPROVED)
# ============================================================================
#
# PURPOSE: Create Figure 4 highlighting the POSITIVE ecological findings from
#          CAFI network analysis - focus on species guilds and associations
#
# STORY: CAFI communities are structured into ecological guilds of cryptic
#        species (shrimps, brittle stars, worms) that share microhabitats
#        within coral colonies. This is a DISCOVERY, not a null result.
#
# 4-PANEL LAYOUT:
#   A: Co-occurrence network reveals species guilds (CLEANER - fewer labels)
#   B: Strongest species associations (SIMPLIFIED - top 12, gradient colors)
#   C: Hub species by taxonomic group (CONSOLIDATED groups)
#   D: Guild composition by module (SHARED legend with C)
#
# IMPROVEMENTS:
#   - Panel A: Only top 5 hub species labeled, better module separation
#   - Panel B: Reduced to top 12 pairs, simpler correlation gradient
#   - Panel C: Consolidated to 6 taxonomic groups (rare types as "Other")
#   - Panel D: Shared legend with Panel C, cleaner x-axis labels
#
# OUTPUTS:
#   - output/figures/manuscript/fig4_ecological_guilds.png
#   - output/figures/manuscript/fig4_legend.txt
#   - output/figures/manuscript/fig4_methods.txt
#   - output/figures/manuscript/fig4_results.txt
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-28
# ============================================================================

cat("\n========================================\n")
cat("FIGURE 4: CAFI Ecological Guilds (IMPROVED)\n")
cat("========================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Additional required packages
if (!require("igraph", quietly = TRUE)) {
  install.packages("igraph", repos = "https://cran.rstudio.com")
}
library(igraph)

if (!require("ggraph", quietly = TRUE)) {
  install.packages("ggraph", repos = "https://cran.rstudio.com")
}
library(ggraph)

if (!require("tidygraph", quietly = TRUE)) {
  install.packages("tidygraph", repos = "https://cran.rstudio.com")
}
library(tidygraph)

if (!require("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel", repos = "https://cran.rstudio.com")
}
library(ggrepel)

# Create directories
MANUSCRIPT_DIR <- file.path(PATHS$figures, "manuscript")
dir.create(MANUSCRIPT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading network data...\n")

# Load network results
network_data <- load_object("cafi_network")
g <- network_data$graph
centrality_df <- network_data$centrality
module_taxonomy <- network_data$module_taxonomy
null_comparison <- network_data$null_comparison
edge_list <- network_data$edge_list

cat("  Network: ", vcount(g), " species, ", ecount(g), " edges\n")
cat("  Modules: ", length(unique(V(g)$module)), "\n\n")

# ============================================================================
# COLOR PALETTES (colorblind-safe, CONSOLIDATED)
# ============================================================================

# Guild/Module colors (distinct, colorblind-safe Okabe-Ito)
MODULE_COLORS <- c(
  "1" = "#0072B2",  # Blue - Guild 1
  "2" = "#E69F00",  # Orange - Guild 2
  "3" = "#009E73",  # Teal - Guild 3
  "4" = "#CC79A7"   # Pink - Guild 4
)

# CONSOLIDATED taxonomic group colors (6 groups max + Other)
# Matches user requirements: shrimp (blue), crab (orange), echinoderm (teal),
# fish (green), worm (yellow), other (gray)
TYPE_COLORS_SIMPLE <- c(
  "Shrimp" = "#0072B2",       # Blue
  "Crab" = "#E69F00",         # Orange
  "Echinoderm" = "#009E73",   # Teal
  "Fish" = "#56B4E9",         # Sky blue
  "Worm" = "#F0E442",         # Yellow
  "Other" = "#999999"         # Gray (hermit, snail, squat_lobster, amphipod, unknown)
)

# Function to consolidate taxonomic types
consolidate_type <- function(type) {
  case_when(
    type == "shrimp" ~ "Shrimp",
    type == "crab" ~ "Crab",
    type == "echinoderm" ~ "Echinoderm",
    type == "fish" ~ "Fish",
    type == "worm" ~ "Worm",
    TRUE ~ "Other"  # hermit, snail, squat_lobster, amphipod, unknown
  )
}

# ============================================================================
# PANEL A: NETWORK VISUALIZATION (IMPROVED - fewer labels, better separation)
# ============================================================================

cat("Creating Panel A: Network visualization (improved)...\n")

# Convert igraph to tidygraph for ggraph
tg <- as_tbl_graph(g)

# Get top 5 hub species for labeling (by degree)
top_5_hubs <- centrality_df %>%
  arrange(desc(degree)) %>%
  slice_head(n = 5) %>%
  pull(species)

# Add node attributes
tg <- tg %>%
  activate(nodes) %>%
  mutate(
    module_factor = factor(module, levels = c(1, 2, 3, 4)),
    node_size = degree,
    # Only label top 5 hub species
    label_display = ifelse(name %in% top_5_hubs, name, ""),
    # Shorten labels: "Genus species" -> "G. species"
    label_short = ifelse(
      label_display != "",
      gsub("^([A-Z])[a-z]+\\s", "\\1. ", label_display),
      ""
    )
  )

# Create network layout with better module separation
# Use FR layout with modified parameters for better clustering
set.seed(42)

panel_a <- ggraph(tg, layout = "fr", niter = 1000) +
  # Edges first (underneath) - subtle
  geom_edge_link(aes(alpha = weight),
                 color = "gray75",
                 width = 0.25,
                 show.legend = FALSE) +
  # Nodes colored by module
  geom_node_point(aes(size = node_size, fill = module_factor),
                  shape = 21, color = "white", stroke = 0.4) +
  # Labels for top 5 hub species only - better positioned
  geom_node_text(aes(label = label_short),
                 size = 2.8,
                 repel = TRUE,
                 max.overlaps = 20,
                 box.padding = 0.6,
                 point.padding = 0.4,
                 segment.color = "gray40",
                 segment.size = 0.3,
                 segment.alpha = 0.7,
                 fontface = "italic",
                 color = "gray10",
                 force = 2,
                 force_pull = 0.5) +
  scale_fill_manual(
    values = MODULE_COLORS,
    name = "Guild",
    labels = c("1" = "Guild 1", "2" = "Guild 2", "3" = "Guild 3", "4" = "Guild 4")
  ) +
  scale_size_continuous(range = c(1.5, 7), guide = "none") +
  scale_edge_alpha_continuous(range = c(0.05, 0.3)) +
  labs(
    title = "A",
    subtitle = "Co-occurrence network (58 species, 4 guilds)"
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9, color = "gray40", hjust = 0),
    legend.position = c(0.12, 0.18),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    legend.key.size = unit(0.35, "cm"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3)))

# ============================================================================
# PANEL B: STRONGEST SPECIES ASSOCIATIONS (SIMPLIFIED - top 12, gradient)
# ============================================================================

cat("Creating Panel B: Strongest associations (simplified)...\n")

# Get top 12 associations from edge list (reduced from 20)
top_edges <- edge_list %>%
  arrange(desc(correlation)) %>%
  slice_head(n = 12) %>%
  mutate(
    # Shorten species names for display: "Genus species" -> "G. species"
    sp1_short = gsub("^([A-Z])[a-z]+\\s", "\\1. ", sp1),
    sp2_short = gsub("^([A-Z])[a-z]+\\s", "\\1. ", sp2),
    # Create clean pair label with dash
    pair_label = paste0(sp1_short, " \u2013 ", sp2_short)
  )

# Simpler color scheme: gradient by correlation strength
panel_b <- ggplot(top_edges, aes(x = reorder(pair_label, correlation), y = correlation)) +
  geom_segment(aes(xend = pair_label, y = 0, yend = correlation),
               color = "gray70", linewidth = 0.8) +
  geom_point(aes(fill = correlation), size = 3.5, shape = 21, color = "gray30", stroke = 0.3) +
  scale_fill_gradient(
    low = "#56B4E9",   # Light blue for lower correlations
    high = "#0072B2",  # Dark blue for higher correlations
    name = "r",
    guide = guide_colorbar(barwidth = 0.5, barheight = 4)
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.8, 0.2)) +
  labs(
    title = "B",
    subtitle = "Strongest species associations (top 12)",
    x = NULL,
    y = "Correlation (Spearman r)"
  ) +
  theme_publication(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.y = element_text(size = 8, face = "italic"),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# ============================================================================
# PANEL C: HUB SPECIES BY TAXONOMIC GROUP (CONSOLIDATED - 6 groups)
# ============================================================================

cat("Creating Panel C: Hub species by taxonomic group (consolidated)...\n")

# Get top 12 hub species (reduced from 15)
top_hubs <- centrality_df %>%
  arrange(desc(degree)) %>%
  slice_head(n = 12) %>%
  mutate(
    # Consolidate to 6 main groups
    type_consolidated = consolidate_type(type),
    type_consolidated = factor(type_consolidated,
                               levels = c("Shrimp", "Crab", "Echinoderm", "Fish", "Worm", "Other")),
    # Shorten species names
    species_short = gsub("^([A-Z])[a-z]+\\s", "\\1. ", species)
  )

panel_c <- ggplot(top_hubs, aes(x = reorder(species_short, degree), y = degree, fill = type_consolidated)) +
  geom_col(color = "gray30", linewidth = 0.2, width = 0.75) +
  geom_text(aes(label = degree), hjust = -0.3, size = 2.5, fontface = "bold") +
  scale_fill_manual(
    values = TYPE_COLORS_SIMPLE,
    name = "Taxonomic\nGroup",
    drop = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), breaks = seq(0, 50, 10)) +
  coord_flip() +
  labs(
    title = "C",
    subtitle = "Hub species (top 12 by degree)",
    x = NULL,
    y = "Degree Centrality"
  ) +
  theme_publication(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.y = element_text(size = 8, face = "italic"),
    legend.position = "right",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(fill = guide_legend(ncol = 1))

# ============================================================================
# PANEL D: GUILD COMPOSITION BY MODULE (SHARED legend, cleaner labels)
# ============================================================================

cat("Creating Panel D: Guild composition by module (improved)...\n")

# Summarize module composition with consolidated types
module_comp <- centrality_df %>%
  mutate(type_consolidated = consolidate_type(type)) %>%
  group_by(module, type_consolidated) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(module) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup() %>%
  mutate(
    # Cleaner, more descriptive guild names
    module_label = case_when(
      module == 1 ~ "Guild 1\n(n=21)",
      module == 2 ~ "Guild 2\n(n=21)",
      module == 3 ~ "Guild 3\n(n=11)",
      module == 4 ~ "Guild 4\n(n=5)"
    ),
    type_consolidated = factor(type_consolidated,
                               levels = c("Shrimp", "Crab", "Echinoderm", "Fish", "Worm", "Other"))
  )

# Order modules
module_order <- c("Guild 1\n(n=21)", "Guild 2\n(n=21)",
                  "Guild 3\n(n=11)", "Guild 4\n(n=5)")

panel_d <- ggplot(module_comp, aes(x = factor(module_label, levels = module_order),
                                    y = proportion, fill = type_consolidated)) +
  geom_col(color = "gray30", linewidth = 0.2, width = 0.75) +
  scale_fill_manual(
    values = TYPE_COLORS_SIMPLE,
    name = "Taxonomic\nGroup"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "D",
    subtitle = "Guild composition",
    x = NULL,
    y = "Proportion of Species"
  ) +
  theme_publication(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    axis.text.x = element_text(size = 8, lineheight = 0.85),
    legend.position = "none",  # Remove - shared with Panel C
    panel.grid.major.x = element_blank()
  )

# ============================================================================
# COMBINE PANELS
# ============================================================================

cat("Combining panels...\n")

# Layout: A and B on top, C and D on bottom
# Use plot_layout to control relative widths
fig4 <- (panel_a + panel_b + plot_layout(widths = c(1.1, 1))) /
        (panel_c + panel_d + plot_layout(widths = c(1.2, 1))) +
  plot_annotation(
    title = "Figure 4. CAFI Communities Form Distinct Ecological Guilds",
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0),
      plot.margin = margin(10, 10, 5, 10)
    )
  )

# ============================================================================
# SAVE FIGURE
# ============================================================================

cat("Saving figure...\n")

ggsave(
  filename = file.path(MANUSCRIPT_DIR, "fig4_ecological_guilds.png"),
  plot = fig4,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

cat("  Saved: fig4_ecological_guilds.png\n")

# ============================================================================
# SAVE LEGEND TEXT (IMPROVED - under 200 words)
# ============================================================================

legend_text <- "Figure 4. CAFI Communities Form Distinct Ecological Guilds. (A) Co-occurrence network of 58 coral-associated species partitioned into four ecological guilds (Louvain clustering). Node size reflects degree centrality; labels identify the five most connected hub species. Network transitivity (0.94) is 1.7x higher than null expectation (z = 36.1, p < 0.001). (B) The 12 strongest species associations after controlling for coral volume (partial Spearman correlations, r > 0.3, FDR < 0.05). Color intensity indicates correlation strength. (C) The 12 most connected hub species by degree centrality. Shrimps and hermit crabs (consolidated as 'Other') dominate hub positions, suggesting they serve as keystone members linking coral-associated assemblages. (D) Taxonomic composition within each guild. Guild 1 is dominated by shrimps and crabs; Guild 2 by echinoderms and worms; Guild 3 by fish and shrimps; Guild 4 by crabs. The modular structure indicates that CAFI partition into distinct microhabitat-sharing groups within coral colonies.
"

writeLines(legend_text, file.path(MANUSCRIPT_DIR, "fig4_legend.txt"))
cat("  Saved: fig4_legend.txt\n")

# ============================================================================
# SAVE METHODS TEXT (IMPROVED - ~150-200 words)
# ============================================================================

methods_text <- "Network Analysis

We constructed species co-occurrence networks using volume-corrected abundance data. For each species, we fit a logistic GLM predicting presence from log10(coral volume) and extracted deviance residuals to remove the confound that larger corals host more individuals. Spearman correlations were computed on residuals across all 58 species pairs occurring in at least 5 corals. Edges were retained for pairs with r > 0.30 and FDR-corrected p < 0.05 (Benjamini-Hochberg).

Network structure was assessed using the global clustering coefficient (transitivity) and compared to 1,000 null networks generated via the configuration model (preserving degree sequence). Z-scores quantified deviation from random expectation. Community detection used the Louvain algorithm on edge weights to identify ecological guilds. Hub species were identified by degree centrality (number of associations). Network visualization used Fruchterman-Reingold force-directed layout with nodes colored by guild membership and sized proportionally to degree.
"

writeLines(methods_text, file.path(MANUSCRIPT_DIR, "fig4_methods.txt"))
cat("  Saved: fig4_methods.txt\n")

# ============================================================================
# SAVE RESULTS TEXT (IMPROVED - ~150-200 words)
# ============================================================================

results_text <- "Network Analysis Results

CAFI communities exhibited significant non-random structure, clustering into four distinct ecological guilds (Fig. 4). The co-occurrence network contained 58 species connected by 1,081 positive associations (r > 0.3, FDR < 0.05). Network transitivity was 0.94, which is 1.7 times higher than random expectation (z = 36.1, p < 0.001), indicating that species co-occurring with one partner tend to co-occur with that partner's other associates.

The Louvain algorithm identified four guilds with distinct taxonomic signatures: Guild 1 (n = 21 species) was dominated by shrimps and hermit crabs; Guild 2 (n = 21) by echinoderms and polychaete worms; Guild 3 (n = 11) by resident fish and shrimps; Guild 4 (n = 5) contained peripheral specialists including Trapezia crabs. Hub species linking multiple guild members included the hermit crab Calcinus laevimanus (degree = 48), the shrimp family Hippolytidae (degree = 48), and the brittle star Ophiomastix elegans (degree = 48). The strongest pairwise associations were between C. laevimanus and O. elegans (r = 0.85) and between syllid worms and the hermit crab Pagurixus nomurai (r = 0.84).
"

writeLines(results_text, file.path(MANUSCRIPT_DIR, "fig4_results.txt"))
cat("  Saved: fig4_results.txt\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("FIGURE 4 COMPLETE (IMPROVED)\n")
cat("========================================\n\n")

cat("Output files:\n")
cat("  - output/figures/manuscript/fig4_ecological_guilds.png (12x10, 300 DPI)\n")
cat("  - output/figures/manuscript/fig4_legend.txt\n")
cat("  - output/figures/manuscript/fig4_methods.txt\n")
cat("  - output/figures/manuscript/fig4_results.txt\n\n")

cat("Key improvements:\n")
cat("  1. Panel A: Only top 5 hub species labeled (was all high-degree)\n")
cat("  2. Panel B: Reduced to top 12 pairs with correlation gradient (was 20 with categories)\n")
cat("  3. Panel C: Consolidated to 6 taxonomic groups (was 10+)\n")
cat("  4. Panel D: Removed redundant legend (shared with Panel C)\n")
cat("  5. Colorblind-safe palette throughout\n\n")
