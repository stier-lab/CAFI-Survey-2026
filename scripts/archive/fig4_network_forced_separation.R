# ============================================================================
# fig4_network_forced_separation.R - Force guilds into separate quadrants
# ============================================================================

cat("\n============================================================\n")
cat("NETWORK WITH FORCED GUILD SEPARATION\n")
cat("============================================================\n\n")

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

network_results <- load_object("cafi_network")
centrality_df <- network_results$centrality
edge_list <- network_results$edge_list

library(igraph)
library(ggrepel)
library(patchwork)

fig_dir <- file.path(PATHS$figures, "06_network")

# Colors
guild_colors <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73", "4" = "#CC79A7")

# Use r > 0.7 threshold
threshold <- 0.7
edges_filtered <- edge_list %>% filter(correlation > threshold)
species_in <- unique(c(edges_filtered$sp1, edges_filtered$sp2))

species_info <- centrality_df %>%
  filter(species %in% species_in) %>%
  dplyr::select(species, type, module, degree, occurrence)

guild_lookup <- setNames(species_info$module, species_info$species)

edges_filtered <- edges_filtered %>%
  mutate(
    guild1 = guild_lookup[sp1],
    guild2 = guild_lookup[sp2],
    same_guild = guild1 == guild2
  )

# ============================================================================
# PLACE GUILDS IN FIXED QUADRANTS
# ============================================================================

cat("Placing guilds in fixed quadrants...\n")

# Guild positions (center of each quadrant)
quadrant_centers <- list(
  "1" = c(3, 3),    # top-right - hermit crabs + shrimp
  "2" = c(-3, 3),   # top-left - worms + echinoderms
  "3" = c(3, -3),   # bottom-right - fish + crabs
  "4" = c(-3, -3)   # bottom-left - peripheral specialists
)

# For each guild, arrange species in a circle around the center
n_nodes <- nrow(species_info)
layout_fixed <- matrix(0, nrow = n_nodes, ncol = 2)

for (guild in 1:4) {
  guild_sp <- species_info %>% filter(module == guild) %>% arrange(desc(degree))
  guild_idx <- match(guild_sp$species, species_info$species)
  n_guild <- length(guild_idx)

  center <- quadrant_centers[[as.character(guild)]]

  if (n_guild > 1) {
    # Arrange in circle, with highest degree species at center
    # First species (highest degree) at center
    layout_fixed[guild_idx[1], ] <- center

    if (n_guild > 1) {
      # Remaining species in circle around center
      angles <- seq(0, 2*pi, length.out = n_guild)[-1]
      radius <- 0.8 + 0.2 * sqrt(n_guild - 1)  # Scale radius with guild size
      for (i in 2:n_guild) {
        layout_fixed[guild_idx[i], 1] <- center[1] + radius * cos(angles[i-1])
        layout_fixed[guild_idx[i], 2] <- center[2] + radius * sin(angles[i-1])
      }
    }
  } else if (n_guild == 1) {
    layout_fixed[guild_idx, ] <- center
  }
}

# Build node dataframe
node_df <- data.frame(
  x = layout_fixed[, 1],
  y = layout_fixed[, 2],
  species = species_info$species,
  module = as.character(species_info$module),
  type = species_info$type,
  degree = species_info$degree
) %>%
  mutate(label = ifelse(nchar(species) > 14, paste0(substr(species, 1, 12), ".."), species))

# Build edge dataframe - within guild only for cleaner viz
within_edges <- edges_filtered %>%
  filter(same_guild) %>%
  mutate(
    x = layout_fixed[match(sp1, species_info$species), 1],
    y = layout_fixed[match(sp1, species_info$species), 2],
    xend = layout_fixed[match(sp2, species_info$species), 1],
    yend = layout_fixed[match(sp2, species_info$species), 2]
  )

between_edges <- edges_filtered %>%
  filter(!same_guild) %>%
  mutate(
    x = layout_fixed[match(sp1, species_info$species), 1],
    y = layout_fixed[match(sp1, species_info$species), 2],
    xend = layout_fixed[match(sp2, species_info$species), 1],
    yend = layout_fixed[match(sp2, species_info$species), 2]
  )

# Guild summary for labels
guild_summary <- species_info %>%
  group_by(module) %>%
  summarise(
    n = n(),
    top_type = names(sort(table(type), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  mutate(
    x = sapply(module, function(m) quadrant_centers[[as.character(m)]][1]),
    y = sapply(module, function(m) quadrant_centers[[as.character(m)]][2]) + 2.2,
    label = paste0("Guild ", module, "\n", n, " species\n(", top_type, "-dominated)")
  )

# ============================================================================
# VERSION 1: WITHIN-GUILD EDGES ONLY (cleanest)
# ============================================================================

cat("Creating version 1: Within-guild edges only...\n")

p1 <- ggplot() +
  # Within-guild edges
  geom_segment(data = within_edges,
               aes(x = x, y = y, xend = xend, yend = yend, alpha = correlation),
               color = "gray30", linewidth = 0.8) +
  # Nodes
  geom_point(data = node_df,
             aes(x = x, y = y, fill = module, size = degree),
             shape = 21, color = "white", stroke = 0.8) +
  # Species labels
  geom_text_repel(data = node_df,
                  aes(x = x, y = y, label = label, color = module),
                  size = 2.8, fontface = "bold",
                  max.overlaps = 30, force = 1,
                  segment.color = "gray60", segment.alpha = 0.5,
                  show.legend = FALSE) +
  # Guild labels
  geom_label(data = guild_summary,
             aes(x = x, y = y, label = label, fill = as.character(module)),
             color = "white", fontface = "bold", size = 3.5,
             label.padding = unit(0.5, "lines")) +
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors) +
  scale_size_continuous(range = c(3, 14), name = "Degree\n(connections)") +
  scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
  labs(
    title = "CAFI Ecological Guilds: Within-Guild Associations",
    subtitle = "Strongest co-occurrences only (r > 0.7) | Each guild in separate quadrant"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-5.5, 6))

ggsave(file.path(fig_dir, "network_quadrant_within_only.png"), p1,
       width = 14, height = 12, dpi = 300, bg = "white")

# ============================================================================
# VERSION 2: WITH BETWEEN-GUILD EDGES (faint)
# ============================================================================

cat("Creating version 2: With between-guild edges...\n")

p2 <- ggplot() +
  # Between-guild edges (very faint curved lines)
  geom_curve(data = between_edges,
             aes(x = x, y = y, xend = xend, yend = yend),
             color = "gray80", linewidth = 0.3, alpha = 0.3,
             curvature = 0.2) +
  # Within-guild edges
  geom_segment(data = within_edges,
               aes(x = x, y = y, xend = xend, yend = yend, alpha = correlation),
               color = "gray30", linewidth = 0.8) +
  # Nodes
  geom_point(data = node_df,
             aes(x = x, y = y, fill = module, size = degree),
             shape = 21, color = "white", stroke = 0.8) +
  # Species labels
  geom_text_repel(data = node_df,
                  aes(x = x, y = y, label = label, color = module),
                  size = 2.8, fontface = "bold",
                  max.overlaps = 30, force = 1,
                  segment.color = "gray60", segment.alpha = 0.5,
                  show.legend = FALSE) +
  # Guild labels
  geom_label(data = guild_summary,
             aes(x = x, y = y, label = label, fill = as.character(module)),
             color = "white", fontface = "bold", size = 3.5,
             label.padding = unit(0.5, "lines")) +
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors) +
  scale_size_continuous(range = c(3, 14), name = "Degree\n(connections)") +
  scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
  labs(
    title = "CAFI Ecological Guilds",
    subtitle = sprintf("38 species, 161 edges (r > 0.7) | %d within-guild, %d between-guild edges",
                       nrow(within_edges), nrow(between_edges))
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-5.5, 6))

ggsave(file.path(fig_dir, "network_quadrant_all_edges.png"), p2,
       width = 14, height = 12, dpi = 300, bg = "white")

# ============================================================================
# VERSION 3: SIMPLIFIED - LARGER NODES, TOP 3 LABELS PER GUILD
# ============================================================================

cat("Creating version 3: Simplified with top species only...\n")

top_labels <- node_df %>%
  group_by(module) %>%
  slice_max(degree, n = 4) %>%
  ungroup()

p3 <- ggplot() +
  # Within-guild edges
  geom_segment(data = within_edges,
               aes(x = x, y = y, xend = xend, yend = yend, alpha = correlation),
               color = "gray30", linewidth = 1) +
  # Nodes - larger
  geom_point(data = node_df,
             aes(x = x, y = y, fill = module, size = degree),
             shape = 21, color = "white", stroke = 1) +
  # Only top species labels
  geom_text_repel(data = top_labels,
                  aes(x = x, y = y, label = species, color = module),
                  size = 3.5, fontface = "bold",
                  max.overlaps = 20,
                  segment.color = "gray50",
                  bg.color = "white", bg.r = 0.15,
                  show.legend = FALSE) +
  # Guild labels
  geom_label(data = guild_summary,
             aes(x = x, y = y, label = label, fill = as.character(module)),
             color = "white", fontface = "bold", size = 4,
             label.padding = unit(0.5, "lines")) +
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors) +
  scale_size_continuous(range = c(5, 18), name = "Degree") +
  scale_alpha_continuous(range = c(0.5, 1), guide = "none") +
  labs(
    title = "CAFI Species Form Four Ecological Guilds",
    subtitle = "Hub species labeled | Node size = network connectivity"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    legend.position = "right",
    legend.title = element_text(size = 11),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-5.5, 6))

ggsave(file.path(fig_dir, "network_quadrant_simplified.png"), p3,
       width = 14, height = 12, dpi = 300, bg = "white")

# ============================================================================
# VERSION 4: COLOR BY TAXONOMIC TYPE INSTEAD OF GUILD
# ============================================================================

cat("Creating version 4: Colored by taxonomic type...\n")

type_colors <- c(
  "crab" = "#E69F00", "shrimp" = "#0072B2", "fish" = "#009E73",
  "echinoderm" = "#999999", "worm" = "#56B4E9", "snail" = "#F0E442",
  "hermit" = "#CC79A7", "squat_lobster" = "#D55E00", "amphipod" = "#666666"
)

# Simplified guild labels (just guild number)
guild_labels_simple <- guild_summary %>%
  mutate(label = paste0("Guild ", module, "\n(n=", n, ")"))

p4 <- ggplot() +
  # Within-guild edges
  geom_segment(data = within_edges,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "gray50", linewidth = 0.6, alpha = 0.6) +
  # Nodes colored by type
  geom_point(data = node_df,
             aes(x = x, y = y, fill = type, size = degree),
             shape = 21, color = "white", stroke = 0.8) +
  # Top species labels
  geom_text_repel(data = top_labels,
                  aes(x = x, y = y, label = species),
                  size = 3, fontface = "bold",
                  max.overlaps = 20,
                  segment.color = "gray50",
                  bg.color = "white", bg.r = 0.12) +
  # Guild labels
  geom_label(data = guild_labels_simple,
             aes(x = x, y = y, label = label),
             fill = "gray20", color = "white", fontface = "bold", size = 4,
             label.padding = unit(0.4, "lines")) +
  scale_fill_manual(values = type_colors, name = "Taxonomic\nGroup") +
  scale_size_continuous(range = c(4, 16), name = "Degree") +
  labs(
    title = "CAFI Ecological Guilds by Taxonomic Composition",
    subtitle = "Guilds contain different mixes of taxa | Node color = taxonomic group"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-5.5, 6))

ggsave(file.path(fig_dir, "network_quadrant_by_taxa.png"), p4,
       width = 14, height = 12, dpi = 300, bg = "white")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("FORCED SEPARATION NETWORKS COMPLETE\n")
cat("============================================================\n\n")

cat("Files saved:\n")
cat("  - network_quadrant_within_only.png (cleanest - within-guild edges only)\n")
cat("  - network_quadrant_all_edges.png (with faint between-guild edges)\n")
cat("  - network_quadrant_simplified.png (larger nodes, hub species labeled)\n")
cat("  - network_quadrant_by_taxa.png (colored by taxonomic type)\n")
