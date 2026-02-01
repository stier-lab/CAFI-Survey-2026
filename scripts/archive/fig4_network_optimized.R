# ============================================================================
# fig4_network_optimized.R - Optimize network visualization for guild structure
# ============================================================================

cat("\n============================================================\n")
cat("OPTIMIZING NETWORK VISUALIZATION FOR GUILD STRUCTURE\n")
cat("============================================================\n\n")

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

network_results <- load_object("cafi_network")
centrality_df <- network_results$centrality
edge_list <- network_results$edge_list

library(igraph)
library(ggrepel)
library(patchwork)
library(ggforce)  # For geom_mark_ellipse

fig_dir <- file.path(PATHS$figures, "06_network")

# Colors
guild_colors <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73", "4" = "#CC79A7")
guild_colors_light <- c("1" = "#E69F0030", "2" = "#56B4E930", "3" = "#009E7330", "4" = "#CC79A730")

# Use r > 0.7 threshold for cleaner network
threshold <- 0.7
edges_filtered <- edge_list %>% filter(correlation > threshold)
species_in <- unique(c(edges_filtered$sp1, edges_filtered$sp2))

cat("Using threshold r >", threshold, "\n")
cat("Species:", length(species_in), "\n")
cat("Edges:", nrow(edges_filtered), "\n\n")

# Get guild info for filtered species
species_info <- centrality_df %>%
  filter(species %in% species_in) %>%
  dplyr::select(species, type, module, degree, occurrence)

# Add guild info to edges
guild_lookup <- setNames(species_info$module, species_info$species)
edges_filtered <- edges_filtered %>%
  mutate(
    guild1 = guild_lookup[sp1],
    guild2 = guild_lookup[sp2],
    same_guild = guild1 == guild2
  )

cat("Within-guild edges:", sum(edges_filtered$same_guild), "\n")
cat("Between-guild edges:", sum(!edges_filtered$same_guild), "\n\n")

# Build graph
g <- graph_from_data_frame(
  edges_filtered %>% dplyr::select(sp1, sp2, correlation, same_guild) %>%
    rename(from = sp1, to = sp2, weight = correlation),
  directed = FALSE,
  vertices = species_info %>% rename(name = species)
)

# ============================================================================
# APPROACH 1: Force-directed with guild-based initial positions
# ============================================================================

cat("Approach 1: Force-directed with guild seeding...\n")

# Set initial positions by guild (in quadrants)
n_nodes <- vcount(g)
init_pos <- matrix(0, nrow = n_nodes, ncol = 2)

guild_centers <- list(
  "1" = c(2, -2),   # bottom-right
  "2" = c(-2, -2),  # bottom-left
  "3" = c(2, 2),    # top-right
  "4" = c(-2, 2)    # top-left
)

for (i in 1:n_nodes) {
  guild <- as.character(V(g)$module[i])
  center <- guild_centers[[guild]]
  # Add jitter around guild center
  init_pos[i, ] <- center + rnorm(2, 0, 0.5)
}

set.seed(42)
layout1 <- layout_with_fr(g, coords = init_pos, niter = 500, weights = E(g)$weight)

# ============================================================================
# APPROACH 2: Layout with stronger repulsion between guilds
# ============================================================================

cat("Approach 2: Kamada-Kawai layout...\n")

set.seed(42)
layout2 <- layout_with_kk(g, weights = E(g)$weight)

# ============================================================================
# APPROACH 3: Graphopt - good for clustered networks
# ============================================================================

cat("Approach 3: Graphopt layout...\n")

set.seed(42)
layout3 <- layout_with_graphopt(g, charge = 0.02, spring.length = 50)

# ============================================================================
# APPROACH 4: Manual guild-separated layout
# ============================================================================

cat("Approach 4: Manual guild-separated layout...\n")

# Position each guild in its own region, then use FR within each
layout4 <- matrix(0, nrow = n_nodes, ncol = 2)

for (guild in 1:4) {
  guild_idx <- which(V(g)$module == guild)
  n_guild <- length(guild_idx)

  if (n_guild > 0) {
    # Get center for this guild
    center <- guild_centers[[as.character(guild)]]

    # Create subgraph for this guild
    if (n_guild > 1) {
      # Circular layout within guild, scaled
      angles <- seq(0, 2*pi, length.out = n_guild + 1)[1:n_guild]
      radius <- sqrt(n_guild) * 0.3
      layout4[guild_idx, 1] <- center[1] + radius * cos(angles)
      layout4[guild_idx, 2] <- center[2] + radius * sin(angles)
    } else {
      layout4[guild_idx, ] <- center
    }
  }
}

# Refine with short FR iteration
set.seed(42)
layout4 <- layout_with_fr(g, coords = layout4, niter = 100, weights = E(g)$weight)

# ============================================================================
# HELPER FUNCTION TO CREATE NETWORK PLOT
# ============================================================================

make_guild_network <- function(layout, title, show_between = TRUE, show_hulls = TRUE) {

  node_df <- data.frame(
    x = layout[, 1],
    y = layout[, 2],
    species = V(g)$name,
    module = as.character(V(g)$module),
    type = V(g)$type,
    degree = igraph::degree(g)
  ) %>%
    mutate(
      # Shorten long names
      label = ifelse(nchar(species) > 15,
                     paste0(substr(species, 1, 12), ".."),
                     species)
    )

  # Edge data
  el <- as_edgelist(g, names = FALSE)
  edge_df <- data.frame(
    x = layout[el[,1], 1], y = layout[el[,1], 2],
    xend = layout[el[,2], 1], yend = layout[el[,2], 2],
    weight = E(g)$weight,
    same_guild = E(g)$same_guild
  )

  # Separate within and between guild edges
  within_edges <- edge_df %>% filter(same_guild)
  between_edges <- edge_df %>% filter(!same_guild)

  # Start plot
  p <- ggplot()

  # Add guild hulls/ellipses first (background)
  if (show_hulls) {
    p <- p +
      geom_mark_ellipse(data = node_df,
                        aes(x = x, y = y, fill = module, label = paste("Guild", module)),
                        alpha = 0.15, expand = unit(3, "mm"),
                        label.fontsize = 10, label.fill = "white",
                        con.cap = 0, show.legend = FALSE)
  }

  # Add between-guild edges (faint)
  if (show_between && nrow(between_edges) > 0) {
    p <- p +
      geom_segment(data = between_edges,
                   aes(x = x, y = y, xend = xend, yend = yend),
                   color = "gray80", linewidth = 0.3, alpha = 0.4)
  }

  # Add within-guild edges (stronger)
  p <- p +
    geom_segment(data = within_edges,
                 aes(x = x, y = y, xend = xend, yend = yend, alpha = weight),
                 color = "gray40", linewidth = 0.6)

  # Add nodes
  p <- p +
    geom_point(data = node_df,
               aes(x = x, y = y, fill = module, size = degree),
               shape = 21, color = "white", stroke = 0.8)

  # Add labels
  p <- p +
    geom_text_repel(data = node_df,
                    aes(x = x, y = y, label = label),
                    size = 2.8, max.overlaps = 25,
                    segment.color = "gray60", segment.alpha = 0.5,
                    bg.color = "white", bg.r = 0.12,
                    force = 2, force_pull = 0.5)

  # Styling
  p <- p +
    scale_fill_manual(values = guild_colors, name = "Guild") +
    scale_size_continuous(range = c(3, 12), name = "Degree") +
    scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position = "right"
    ) +
    coord_fixed()

  return(p)
}

# ============================================================================
# CREATE COMPARISON OF APPROACHES
# ============================================================================

cat("\nCreating comparison figure...\n")

p1 <- make_guild_network(layout1, "A. FR with guild seeding", show_hulls = FALSE)
p2 <- make_guild_network(layout2, "B. Kamada-Kawai", show_hulls = FALSE)
p3 <- make_guild_network(layout3, "C. Graphopt", show_hulls = FALSE)
p4 <- make_guild_network(layout4, "D. Manual separation + FR", show_hulls = FALSE)

fig_compare <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Network Layout Comparison (r > 0.7)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(fig_dir, "network_layout_comparison.png"), fig_compare,
       width = 16, height = 14, dpi = 300, bg = "white")

# ============================================================================
# BEST LAYOUT WITH GUILD ELLIPSES
# ============================================================================

cat("Creating final versions with guild ellipses...\n")

# Use the manual separation approach (layout4) - seems cleanest
p_final_hulls <- make_guild_network(layout4,
  "CAFI Co-occurrence Network: Four Ecological Guilds",
  show_between = TRUE, show_hulls = TRUE)

ggsave(file.path(fig_dir, "network_final_with_hulls.png"), p_final_hulls,
       width = 14, height = 12, dpi = 300, bg = "white")

# ============================================================================
# VERSION WITHOUT BETWEEN-GUILD EDGES (cleaner)
# ============================================================================

cat("Creating version without between-guild edges...\n")

p_within_only <- make_guild_network(layout4,
  "CAFI Co-occurrence Network: Within-Guild Associations Only",
  show_between = FALSE, show_hulls = TRUE)

ggsave(file.path(fig_dir, "network_within_guild_only.png"), p_within_only,
       width = 14, height = 12, dpi = 300, bg = "white")

# ============================================================================
# SIMPLIFIED VERSION - LARGER NODES, FEWER LABELS
# ============================================================================

cat("Creating simplified version...\n")

node_df <- data.frame(
  x = layout4[, 1],
  y = layout4[, 2],
  species = V(g)$name,
  module = as.character(V(g)$module),
  type = V(g)$type,
  degree = igraph::degree(g)
)

el <- as_edgelist(g, names = FALSE)
edge_df <- data.frame(
  x = layout4[el[,1], 1], y = layout4[el[,1], 2],
  xend = layout4[el[,2], 1], yend = layout4[el[,2], 2],
  weight = E(g)$weight,
  same_guild = E(g)$same_guild
)

# Only label top 3 species per guild
top_labels <- node_df %>%
  group_by(module) %>%
  slice_max(degree, n = 3) %>%
  ungroup() %>%
  mutate(label = ifelse(nchar(species) > 12,
                        paste0(substr(species, 1, 10), ".."),
                        species))

# Guild annotations
guild_annot <- node_df %>%
  group_by(module) %>%
  summarise(
    x = mean(x),
    y = max(y) + 0.5,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("Guild ", module, "\n(n=", n, ")"))

p_simple <- ggplot() +
  # Guild ellipses
  geom_mark_ellipse(data = node_df,
                    aes(x = x, y = y, fill = module),
                    alpha = 0.2, expand = unit(5, "mm"),
                    show.legend = FALSE) +
  # Within-guild edges only
  geom_segment(data = edge_df %>% filter(same_guild),
               aes(x = x, y = y, xend = xend, yend = yend, alpha = weight),
               color = "gray30", linewidth = 0.7) +
  # Nodes - larger
  geom_point(data = node_df,
             aes(x = x, y = y, fill = module, size = degree),
             shape = 21, color = "white", stroke = 1) +
  # Only top species labels
  geom_text_repel(data = top_labels,
                  aes(x = x, y = y, label = label),
                  size = 3.5, fontface = "bold",
                  max.overlaps = 20,
                  segment.color = "gray50",
                  bg.color = "white", bg.r = 0.15) +
  # Guild labels
  geom_label(data = guild_annot,
             aes(x = x, y = y, label = label, fill = module),
             color = "white", fontface = "bold", size = 4,
             label.padding = unit(0.4, "lines")) +
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_size_continuous(range = c(4, 15), name = "Degree\n(connectivity)") +
  scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
  labs(
    title = "CAFI Species Form Four Ecological Guilds",
    subtitle = "Network shows strongest co-occurrence associations (r > 0.7) | 38 species, 161 edges"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  coord_fixed()

ggsave(file.path(fig_dir, "network_simplified.png"), p_simple,
       width = 14, height = 12, dpi = 300, bg = "white")

# ============================================================================
# RADIAL/CIRCULAR LAYOUT BY GUILD
# ============================================================================

cat("Creating circular guild layout...\n")

# Place guilds at cardinal directions
guild_angles <- c("1" = 0, "2" = pi/2, "3" = pi, "4" = 3*pi/2)
guild_radius <- 3

layout_circular <- matrix(0, nrow = n_nodes, ncol = 2)

for (guild in 1:4) {
  guild_idx <- which(V(g)$module == guild)
  n_guild <- length(guild_idx)

  if (n_guild > 0) {
    angle <- guild_angles[[as.character(guild)]]
    center_x <- guild_radius * cos(angle)
    center_y <- guild_radius * sin(angle)

    # Arrange species in arc around guild center
    if (n_guild > 1) {
      arc_spread <- pi/3  # 60 degree arc
      arc_angles <- seq(angle - arc_spread/2, angle + arc_spread/2, length.out = n_guild)
      inner_radius <- 1.5
      layout_circular[guild_idx, 1] <- center_x + inner_radius * cos(arc_angles - angle) * sign(cos(angle) + 0.001)
      layout_circular[guild_idx, 2] <- center_y + inner_radius * sin(arc_angles - angle) * sign(sin(angle) + 0.001)
    } else {
      layout_circular[guild_idx, 1] <- center_x
      layout_circular[guild_idx, 2] <- center_y
    }
  }
}

# Refine slightly
set.seed(42)
layout_circular <- layout_with_fr(g, coords = layout_circular, niter = 50)

node_df_circ <- data.frame(
  x = layout_circular[, 1],
  y = layout_circular[, 2],
  species = V(g)$name,
  module = as.character(V(g)$module),
  type = V(g)$type,
  degree = igraph::degree(g)
)

el <- as_edgelist(g, names = FALSE)
edge_df_circ <- data.frame(
  x = layout_circular[el[,1], 1], y = layout_circular[el[,1], 2],
  xend = layout_circular[el[,2], 1], yend = layout_circular[el[,2], 2],
  weight = E(g)$weight,
  same_guild = E(g)$same_guild
)

top_labels_circ <- node_df_circ %>%
  group_by(module) %>%
  slice_max(degree, n = 3) %>%
  ungroup() %>%
  mutate(label = ifelse(nchar(species) > 12,
                        paste0(substr(species, 1, 10), ".."),
                        species))

guild_annot_circ <- node_df_circ %>%
  group_by(module) %>%
  summarise(
    x = mean(x) * 1.4,
    y = mean(y) * 1.4,
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("Guild ", module, "\n(n=", n, ")"))

p_circular <- ggplot() +
  # Between-guild edges (very faint)
  geom_segment(data = edge_df_circ %>% filter(!same_guild),
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "gray85", linewidth = 0.3, alpha = 0.5) +
  # Within-guild edges
  geom_segment(data = edge_df_circ %>% filter(same_guild),
               aes(x = x, y = y, xend = xend, yend = yend, alpha = weight),
               color = "gray30", linewidth = 0.8) +
  # Nodes
  geom_point(data = node_df_circ,
             aes(x = x, y = y, fill = module, size = degree),
             shape = 21, color = "white", stroke = 1) +
  # Labels
  geom_text_repel(data = top_labels_circ,
                  aes(x = x, y = y, label = label),
                  size = 3.2, fontface = "bold",
                  max.overlaps = 20,
                  segment.color = "gray50",
                  bg.color = "white", bg.r = 0.12) +
  # Guild labels
  geom_label(data = guild_annot_circ,
             aes(x = x, y = y, label = label, fill = module),
             color = "white", fontface = "bold", size = 4,
             label.padding = unit(0.4, "lines")) +
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_size_continuous(range = c(4, 14), name = "Degree") +
  scale_alpha_continuous(range = c(0.4, 1), guide = "none") +
  labs(
    title = "CAFI Ecological Guilds",
    subtitle = "Strongest co-occurrences (r > 0.7) | Within-guild edges emphasized"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    legend.position = "right"
  ) +
  coord_fixed()

ggsave(file.path(fig_dir, "network_circular_guilds.png"), p_circular,
       width = 14, height = 12, dpi = 300, bg = "white")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("NETWORK VISUALIZATIONS COMPLETE\n")
cat("============================================================\n\n")

cat("Files saved to", fig_dir, ":\n")
cat("  - network_layout_comparison.png (4 layout algorithms compared)\n")
cat("  - network_final_with_hulls.png (ellipses around guilds)\n")
cat("  - network_within_guild_only.png (no between-guild edges)\n")
cat("  - network_simplified.png (larger nodes, top species labeled)\n")
cat("  - network_circular_guilds.png (guilds in cardinal directions)\n")
