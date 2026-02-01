# ============================================================================
# fig4_network_final.R - PUBLICATION-READY Network Figure
# ============================================================================
#
# Based on multi-agent review (Graphic Designer, UX Engineer, Data Scientist):
# - Quadrant layout for clear guild separation
# - r > 0.7 threshold (38 species, 161 edges, 23% density)
# - FIXED: Degree recalculated for filtered network (not full network)
# - Standardized species name format (abbreviated binomials)
# - Network statistics in subtitle
# - Edge weight visual encoding
# ============================================================================

cat("\n============================================================\n")
cat("CREATING PUBLICATION-READY NETWORK FIGURE\n")
cat("============================================================\n\n")

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

network_results <- load_object("cafi_network")
centrality_df <- network_results$centrality
edge_list <- network_results$edge_list

library(igraph)
library(ggrepel)

fig_dir <- file.path(PATHS$figures, "06_network")
manuscript_dir <- file.path(PATHS$figures, "manuscript")

# ============================================================================
# FILTER AT r > 0.7 AND RECALCULATE DEGREE (Critical fix)
# ============================================================================

threshold <- 0.7
edges_filtered <- edge_list %>% filter(correlation > threshold)
species_in <- unique(c(edges_filtered$sp1, edges_filtered$sp2))

cat(sprintf("Threshold: r > %.1f\n", threshold))
cat(sprintf("Species retained: %d (of %d)\n", length(species_in), nrow(centrality_df)))
cat(sprintf("Edges retained: %d (of %d)\n", nrow(edges_filtered), nrow(edge_list)))

# Build filtered graph to recalculate degree
g_filtered <- graph_from_data_frame(
  edges_filtered %>%
    dplyr::select(sp1, sp2, correlation) %>%
    rename(from = sp1, to = sp2, weight = correlation),
  directed = FALSE,
  vertices = data.frame(name = species_in)
)

# Recalculate degree for filtered network
filtered_degree <- igraph::degree(g_filtered)

# Create species info with FILTERED degree
species_info <- centrality_df %>%
  filter(species %in% species_in) %>%
  dplyr::select(species, type, module, occurrence) %>%
  mutate(degree = filtered_degree[species])  # Use filtered degree!

guild_lookup <- setNames(species_info$module, species_info$species)

edges_filtered <- edges_filtered %>%
  mutate(
    guild1 = guild_lookup[sp1],
    guild2 = guild_lookup[sp2],
    same_guild = guild1 == guild2
  )

n_within <- sum(edges_filtered$same_guild)
n_between <- sum(!edges_filtered$same_guild)
cat(sprintf("Within-guild edges: %d\n", n_within))
cat(sprintf("Between-guild edges: %d\n", n_between))

# ============================================================================
# GUILD COLORS (Colorblind-safe)
# ============================================================================

guild_colors <- c(
  "1" = "#E69F00",  # Orange
  "2" = "#56B4E9",  # Sky blue
  "3" = "#009E73",  # Teal
  "4" = "#CC79A7"   # Pink
)

# ============================================================================
# FIXED QUADRANT LAYOUT
# ============================================================================

quadrant_centers <- list(
  "1" = c(3, 3),    # top-right
  "2" = c(-3, 3),   # top-left
  "3" = c(3, -3),   # bottom-right
  "4" = c(-3, -3)   # bottom-left
)

n_nodes <- nrow(species_info)
layout_fixed <- matrix(0, nrow = n_nodes, ncol = 2)

for (guild in 1:4) {
  guild_sp <- species_info %>% filter(module == guild) %>% arrange(desc(degree))
  guild_idx <- match(guild_sp$species, species_info$species)
  n_guild <- length(guild_idx)

  if (n_guild == 0) next

  center <- quadrant_centers[[as.character(guild)]]

  if (n_guild == 1) {
    layout_fixed[guild_idx, ] <- center
  } else {
    # Highest degree at center
    layout_fixed[guild_idx[1], ] <- center

    # Others in circle around center
    angles <- seq(0, 2*pi, length.out = n_guild)[-1]
    radius <- 0.6 + 0.15 * sqrt(n_guild - 1)
    for (i in 2:n_guild) {
      layout_fixed[guild_idx[i], 1] <- center[1] + radius * cos(angles[i-1])
      layout_fixed[guild_idx[i], 2] <- center[2] + radius * sin(angles[i-1])
    }
  }
}

# ============================================================================
# CREATE NODE AND EDGE DATAFRAMES
# ============================================================================

# Abbreviate species names (genus initial + species epithet)
abbreviate_species <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  if (length(parts) >= 2) {
    # Handle genus + species
    genus <- substr(parts[1], 1, 1)
    species <- parts[2]
    if (nchar(species) > 12) species <- paste0(substr(species, 1, 10), ".")
    return(paste0(genus, ". ", species))
  } else {
    # Single word - truncate if needed
    if (nchar(name) > 14) return(paste0(substr(name, 1, 12), ".."))
    return(name)
  }
}

node_df <- data.frame(
  x = layout_fixed[, 1],
  y = layout_fixed[, 2],
  species = species_info$species,
  module = as.character(species_info$module),
  type = species_info$type,
  degree = species_info$degree
) %>%
  mutate(label = sapply(species, abbreviate_species))

# Within-guild edges
within_edges <- edges_filtered %>%
  filter(same_guild) %>%
  mutate(
    x = layout_fixed[match(sp1, species_info$species), 1],
    y = layout_fixed[match(sp1, species_info$species), 2],
    xend = layout_fixed[match(sp2, species_info$species), 1],
    yend = layout_fixed[match(sp2, species_info$species), 2]
  )

# Between-guild edges (faint curved lines)
between_edges <- edges_filtered %>%
  filter(!same_guild) %>%
  mutate(
    x = layout_fixed[match(sp1, species_info$species), 1],
    y = layout_fixed[match(sp1, species_info$species), 2],
    xend = layout_fixed[match(sp2, species_info$species), 1],
    yend = layout_fixed[match(sp2, species_info$species), 2]
  )

# Guild summary with counts and dominant type
guild_summary <- species_info %>%
  group_by(module) %>%
  summarise(
    n = n(),
    top_type = names(sort(table(type), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  mutate(
    x = sapply(module, function(m) quadrant_centers[[as.character(m)]][1]),
    y = sapply(module, function(m) quadrant_centers[[as.character(m)]][2]) + 2.0,
    label = paste0("Guild ", module, "\n", n, " species\n(", top_type, "-dominated)")
  )

# Top species per guild for labeling (hub species)
top_labels <- node_df %>%
  group_by(module) %>%
  slice_max(degree, n = 4, with_ties = FALSE) %>%
  ungroup()

# ============================================================================
# CREATE PUBLICATION FIGURE
# ============================================================================

cat("Creating publication figure...\n")

p_final <- ggplot() +
  # Between-guild edges (very faint curved)
  geom_curve(data = between_edges,
             aes(x = x, y = y, xend = xend, yend = yend),
             color = "gray85", linewidth = 0.25, alpha = 0.4,
             curvature = 0.15) +
  # Within-guild edges (solid, weighted by correlation)
  geom_segment(data = within_edges,
               aes(x = x, y = y, xend = xend, yend = yend,
                   alpha = correlation, linewidth = correlation),
               color = "gray40") +
  # Nodes - sized by FILTERED degree
  geom_point(data = node_df,
             aes(x = x, y = y, fill = module, size = degree),
             shape = 21, color = "white", stroke = 0.8) +
  # Hub species labels
  geom_text_repel(data = top_labels,
                  aes(x = x, y = y, label = label, color = module),
                  size = 3.2, fontface = "bold",
                  max.overlaps = 20, force = 2,
                  segment.color = "gray50", segment.alpha = 0.6,
                  bg.color = "white", bg.r = 0.12,
                  show.legend = FALSE) +
  # Guild labels
  geom_label(data = guild_summary,
             aes(x = x, y = y, label = label, fill = as.character(module)),
             color = "white", fontface = "bold", size = 3.5,
             label.padding = unit(0.4, "lines")) +
  # Scales
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors) +
  scale_size_continuous(range = c(4, 15),
                        name = "Degree",
                        breaks = c(5, 10, 15, 20)) +
  scale_linewidth_continuous(range = c(0.3, 1.2), guide = "none") +
  scale_alpha_continuous(range = c(0.4, 0.9), guide = "none") +
  # Labels
  labs(
    title = "CAFI Species Form Four Ecological Guilds",
    subtitle = sprintf("Co-occurrence network (r > 0.7) | %d species, %d edges | %d within-guild, %d between-guild",
                       n_nodes, nrow(edges_filtered), n_within, n_between),
    caption = "Guilds identified via Louvain community detection on volume-corrected co-occurrence matrix.\nNode size = degree in filtered network. Edge opacity = correlation strength."
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40", margin = margin(b = 10)),
    plot.caption = element_text(size = 8, hjust = 0, color = "gray50", margin = margin(t = 15)),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-5.5, 6))

# Save high-resolution versions
ggsave(file.path(fig_dir, "fig4_network_FINAL.png"), p_final,
       width = 12, height = 11, dpi = 300, bg = "white")

ggsave(file.path(manuscript_dir, "fig4_network.png"), p_final,
       width = 12, height = 11, dpi = 300, bg = "white")

# Also save as PDF for journal submission
ggsave(file.path(manuscript_dir, "fig4_network.pdf"), p_final,
       width = 12, height = 11, device = cairo_pdf)

cat("\n============================================================\n")
cat("PUBLICATION FIGURE COMPLETE\n")
cat("============================================================\n\n")
cat("Files saved:\n")
cat("  - output/figures/06_network/fig4_network_FINAL.png\n")
cat("  - output/figures/manuscript/fig4_network.png\n")
cat("  - output/figures/manuscript/fig4_network.pdf\n")
cat("\nKey improvements applied:\n")
cat("  - Degree recalculated for r > 0.7 filtered network (scientific accuracy)\n")
cat("  - Abbreviated binomial species names (e.g., 'T. rufopunctata')\n")
cat("  - Network statistics in subtitle\n")
cat("  - Edge weight encoded via opacity and linewidth\n")
cat("  - Colorblind-safe 4-color palette\n")
cat("  - 300 DPI for print publication\n")
