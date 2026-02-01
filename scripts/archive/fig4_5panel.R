# ============================================================================
# Figure 4: 5-Panel Network Figure
# Panel A: Overview circular with 4 guild groups (no species labels)
# Panels B-E: Individual guild networks with species labels (force layout)
# ============================================================================

library(tidyverse)
library(igraph)
library(ggforce)
library(ggrepel)
library(patchwork)

source("scripts/00_setup.R")

fig_dir <- file.path(PATHS$figures, "06_network")
manuscript_dir <- file.path(PATHS$figures, "manuscript")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("CREATING 5-PANEL NETWORK FIGURE\n")
cat("============================================================\n\n")

# Load data
network_results <- load_object("cafi_network")
g <- network_results$graph
communities <- network_results$communities

# Species info
sp_info <- data.frame(
  species = V(g)$name,
  guild = membership(communities),
  degree = degree(g),
  type = V(g)$type,
  stringsAsFactors = FALSE
) %>%
  mutate(guild = factor(guild))

# Guild setup
guild_names <- c(
  "1" = "Protective\nCrustaceans",
  "2" = "Core\nAssociates",
  "3" = "Fish &\nTrapezia",
  "4" = "Peripheral\nCrabs"
)

guild_names_short <- c(
  "1" = "Protective Crustaceans",
  "2" = "Core Associates",
  "3" = "Fish & Trapezia",
  "4" = "Peripheral Crabs"
)

guild_counts <- sp_info %>% count(guild)

guild_colors <- c(
  "1" = "#0072B2",
  "2" = "#D55E00",
  "3" = "#009E73",
  "4" = "#CC79A7"
)

cat("Guild sizes:\n")
for (i in 1:4) {
  cat(sprintf("  Guild %d: %d species\n", i, guild_counts$n[i]))
}

# ============================================================================
# PANEL A: Overview circular with guild aggregates
# ============================================================================

cat("\nBuilding Panel A (overview)...\n")

# Create aggregate guild network
# Count edges between guilds
edge_df <- as_data_frame(g, what = "edges")
edge_df <- edge_df %>%
  left_join(sp_info %>% dplyr::select(species, guild), by = c("from" = "species")) %>%
  rename(guild_from = guild) %>%
  left_join(sp_info %>% dplyr::select(species, guild), by = c("to" = "species")) %>%
  rename(guild_to = guild)

# Aggregate edges between guilds
guild_edges <- edge_df %>%
  mutate(
    g1 = pmin(as.numeric(guild_from), as.numeric(guild_to)),
    g2 = pmax(as.numeric(guild_from), as.numeric(guild_to))
  ) %>%
  group_by(g1, g2) %>%
  summarise(n_edges = n(), mean_weight = mean(weight), .groups = "drop")

cat(sprintf("  Guild-level edges: %d\n", nrow(guild_edges)))

# Position guilds in circle
n_guilds <- 4
guild_angles <- seq(pi/2, pi/2 - 2*pi * (1 - 1/n_guilds), length.out = n_guilds)
guild_positions <- data.frame(
  guild = factor(1:4),
  x = 3 * cos(guild_angles),
  y = 3 * sin(guild_angles),
  n = guild_counts$n,
  name = guild_names[as.character(1:4)],
  angle = guild_angles
)

# Node size by species count
guild_positions <- guild_positions %>%
  mutate(node_size = scales::rescale(n, to = c(15, 35)))

# Create bezier curves for guild-level edges
create_guild_bezier <- function(x1, y1, x2, y2, n = 40) {
  dist <- sqrt((x2-x1)^2 + (y2-y1)^2)
  pull <- 0.3
  cx <- (x1 + x2) / 2 * (1 - pull)
  cy <- (y1 + y2) / 2 * (1 - pull)
  t <- seq(0, 1, length.out = n)
  data.frame(
    x = (1-t)^2 * x1 + 2*(1-t)*t * cx + t^2 * x2,
    y = (1-t)^2 * y1 + 2*(1-t)*t * cy + t^2 * y2
  )
}

guild_bezier <- map_dfr(1:nrow(guild_edges), function(i) {
  row <- guild_edges[i,]
  p1 <- guild_positions[guild_positions$guild == row$g1,]
  p2 <- guild_positions[guild_positions$guild == row$g2,]
  pts <- create_guild_bezier(p1$x, p1$y, p2$x, p2$y)
  pts$edge_id <- i
  pts$n_edges <- row$n_edges
  pts$is_within <- row$g1 == row$g2
  pts$guild <- as.character(row$g1)
  pts
})

# Panel A plot
p_A <- ggplot() +
  # Between-guild edges (gray arcs)
  geom_path(
    data = guild_bezier %>% filter(!is_within),
    aes(x = x, y = y, group = edge_id, linewidth = n_edges),
    color = "gray60",
    alpha = 0.5
  ) +
  # Within-guild edges (colored loops - represented as arcs)
  geom_arc(
    data = guild_positions %>%
      left_join(guild_edges %>% filter(g1 == g2) %>%
                  mutate(guild = factor(g1)) %>%
                  dplyr::select(guild, n_edges),
                by = "guild"),
    aes(x0 = x * 1.3, y0 = y * 1.3, r = 0.8,
        start = angle - pi/4, end = angle + pi/4,
        color = guild, linewidth = n_edges),
    alpha = 0.6
  ) +
  # Guild nodes
  geom_point(
    data = guild_positions,
    aes(x = x, y = y, fill = guild, size = node_size),
    shape = 21, color = "white", stroke = 1.5
  ) +
  # Guild labels
  geom_text(
    data = guild_positions,
    aes(x = x * 1.55, y = y * 1.55, label = name),
    size = 3.2, fontface = "bold", lineheight = 0.85,
    color = "gray20"
  ) +
  # Species count inside nodes
  geom_text(
    data = guild_positions,
    aes(x = x, y = y, label = n),
    size = 3.5, fontface = "bold", color = "white"
  ) +
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +
  scale_linewidth_continuous(range = c(0.5, 4), guide = "none") +
  coord_fixed(xlim = c(-6, 6), ylim = c(-6, 6)) +
  labs(title = "A. Network Overview",
       subtitle = sprintf("%d species | %d co-occurrences",
                          sum(guild_counts$n), ecount(g))) +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
    plot.margin = margin(5, 5, 5, 5)
  )

# ============================================================================
# PANELS B-E: Individual guild networks with species labels
# Using force-directed layout to maximize space usage
# ============================================================================

create_guild_panel <- function(guild_id, letter) {

  guild_name <- guild_names_short[as.character(guild_id)]
  guild_color <- guild_colors[as.character(guild_id)]

  # Get species in this guild
  guild_species <- sp_info %>% filter(guild == guild_id)

  # Create subgraph for this guild
  g_sub <- induced_subgraph(g, V(g)[V(g)$name %in% guild_species$species])

  if (vcount(g_sub) == 0) {
    return(ggplot() + theme_void() +
             labs(title = paste(letter, guild_name)))
  }

  # Use Fruchterman-Reingold layout for better space usage
  set.seed(42 + guild_id)
  layout_fr <- layout_with_fr(g_sub, niter = 500)

  # Scale layout to fill space
  layout_fr[,1] <- scales::rescale(layout_fr[,1], to = c(-1, 1))
  layout_fr[,2] <- scales::rescale(layout_fr[,2], to = c(-1, 1))

  # Node data
  node_data <- data.frame(
    species = V(g_sub)$name,
    x = layout_fr[,1],
    y = layout_fr[,2],
    degree = degree(g_sub)
  ) %>%
    mutate(node_size = scales::rescale(degree, to = c(3, 12)))

  # Edge data
  if (ecount(g_sub) > 0) {
    edge_list_sub <- as_data_frame(g_sub, what = "edges")
    edge_data <- edge_list_sub %>%
      left_join(node_data %>% dplyr::select(species, x, y),
                by = c("from" = "species")) %>%
      rename(x1 = x, y1 = y) %>%
      left_join(node_data %>% dplyr::select(species, x, y),
                by = c("to" = "species")) %>%
      rename(x2 = x, y2 = y)
  } else {
    edge_data <- NULL
  }

  n_species <- nrow(node_data)
  n_edges <- ifelse(is.null(edge_data), 0, nrow(edge_data))

  # Build plot
  p <- ggplot()

  # Add edges if they exist
  if (!is.null(edge_data) && nrow(edge_data) > 0) {
    p <- p + geom_segment(
      data = edge_data,
      aes(x = x1, y = y1, xend = x2, yend = y2, alpha = weight),
      color = guild_color,
      linewidth = 0.4
    )
  }

  # Add nodes
  p <- p + geom_point(
    data = node_data,
    aes(x = x, y = y, size = node_size),
    fill = guild_color,
    shape = 21,
    color = "white",
    stroke = 0.5
  )

  # Add labels
  p <- p + geom_text_repel(
    data = node_data,
    aes(x = x, y = y, label = species),
    size = 2.4,
    fontface = "italic",
    color = "gray20",
    segment.color = "gray60",
    segment.size = 0.2,
    segment.alpha = 0.5,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = 30,
    force = 2
  )

  p <- p +
    scale_size_identity() +
    scale_alpha_continuous(range = c(0.2, 0.7), guide = "none") +
    coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +
    labs(
      title = sprintf("%s. %s", letter, guild_name),
      subtitle = sprintf("n = %d species, %d edges", n_species, n_edges)
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5,
                                color = guild_color),
      plot.subtitle = element_text(size = 8, hjust = 0.5, color = "gray50"),
      plot.margin = margin(5, 5, 5, 5),
      plot.background = element_rect(fill = "white", color = "gray85",
                                     linewidth = 0.5)
    )

  return(p)
}

cat("Building Panels B-E (individual guilds)...\n")

p_B <- create_guild_panel(1, "B")
p_C <- create_guild_panel(2, "C")
p_D <- create_guild_panel(3, "D")
p_E <- create_guild_panel(4, "E")

# ============================================================================
# COMBINE INTO 5-PANEL FIGURE
# ============================================================================

cat("Combining panels...\n")

# Layout: A on top (full width), B-E in 2x2 grid below
p_bottom <- (p_B | p_C) / (p_D | p_E)

p_final <- p_A / p_bottom +
  plot_layout(heights = c(1.2, 2)) +
  plot_annotation(
    title = "Figure 4: CAFI Co-occurrence Network Structure",
    subtitle = "Louvain community detection reveals four ecological guilds with distinct co-occurrence patterns",
    caption = paste0(
      "A: Guild-level overview showing aggregate connectivity. Numbers indicate species per guild.\n",
      "B-E: Species-level networks within each guild. Node size = degree centrality; ",
      "edge opacity = correlation strength (r > 0.7)."
    ),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 10)),
      plot.caption = element_text(size = 8, hjust = 0, color = "gray50",
                                  margin = margin(t = 10), lineheight = 1.2),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# Save
ggsave(
  file.path(fig_dir, "fig4_5panel.png"),
  p_final,
  width = 12,
  height = 14,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(manuscript_dir, "fig4_network_5panel.png"),
  p_final,
  width = 12,
  height = 14,
  dpi = 300,
  bg = "white"
)

cat("\nSaved: fig4_5panel.png\n")
cat("Saved: manuscript/fig4_network_5panel.png\n")

# ============================================================================
# ALTERNATE: Cleaner version with more space for subpanels
# ============================================================================

cat("\nCreating alternate version with compact overview...\n")

# Make overview panel more compact
p_A_compact <- ggplot() +
  geom_path(
    data = guild_bezier %>% filter(!is_within),
    aes(x = x, y = y, group = edge_id, linewidth = n_edges),
    color = "gray55", alpha = 0.4
  ) +
  geom_arc(
    data = guild_positions %>%
      left_join(guild_edges %>% filter(g1 == g2) %>%
                  mutate(guild = factor(g1)) %>%
                  dplyr::select(guild, n_edges), by = "guild"),
    aes(x0 = x * 1.25, y0 = y * 1.25, r = 0.6,
        start = angle - pi/5, end = angle + pi/5,
        color = guild, linewidth = n_edges),
    alpha = 0.5
  ) +
  geom_point(
    data = guild_positions,
    aes(x = x, y = y, fill = guild, size = node_size),
    shape = 21, color = "white", stroke = 1.2
  ) +
  geom_text(
    data = guild_positions,
    aes(x = x * 1.6, y = y * 1.6, label = paste0(gsub("\n", " ", name), "\n(", n, ")")),
    size = 2.8, fontface = "bold", lineheight = 0.9, color = "gray25"
  ) +
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +
  scale_linewidth_continuous(range = c(0.3, 3), guide = "none") +
  coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-5.5, 5.5)) +
  labs(title = "A. Guild Overview") +
  theme_void() +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2)
  )

# Larger subpanels
create_guild_panel_large <- function(guild_id, letter) {

  guild_name <- guild_names_short[as.character(guild_id)]
  guild_color <- guild_colors[as.character(guild_id)]

  guild_species <- sp_info %>% filter(guild == guild_id)
  g_sub <- induced_subgraph(g, V(g)[V(g)$name %in% guild_species$species])

  if (vcount(g_sub) == 0) return(ggplot() + theme_void())

  set.seed(42 + guild_id)
  layout_fr <- layout_with_fr(g_sub, niter = 500)
  layout_fr[,1] <- scales::rescale(layout_fr[,1], to = c(-1, 1))
  layout_fr[,2] <- scales::rescale(layout_fr[,2], to = c(-1, 1))

  node_data <- data.frame(
    species = V(g_sub)$name,
    x = layout_fr[,1],
    y = layout_fr[,2],
    degree = degree(g_sub)
  ) %>%
    mutate(node_size = scales::rescale(degree, to = c(4, 14)))

  if (ecount(g_sub) > 0) {
    edge_list_sub <- as_data_frame(g_sub, what = "edges")
    edge_data <- edge_list_sub %>%
      left_join(node_data %>% dplyr::select(species, x, y),
                by = c("from" = "species")) %>%
      rename(x1 = x, y1 = y) %>%
      left_join(node_data %>% dplyr::select(species, x, y),
                by = c("to" = "species")) %>%
      rename(x2 = x, y2 = y)
  } else {
    edge_data <- NULL
  }

  p <- ggplot()

  if (!is.null(edge_data) && nrow(edge_data) > 0) {
    p <- p + geom_segment(
      data = edge_data,
      aes(x = x1, y = y1, xend = x2, yend = y2, alpha = weight),
      color = guild_color, linewidth = 0.5
    )
  }

  p <- p +
    geom_point(
      data = node_data,
      aes(x = x, y = y, size = node_size),
      fill = guild_color, shape = 21, color = "white", stroke = 0.6
    ) +
    geom_text_repel(
      data = node_data,
      aes(x = x, y = y, label = species),
      size = 2.8,
      fontface = "italic",
      color = "gray15",
      segment.color = "gray55",
      segment.size = 0.25,
      box.padding = 0.4,
      point.padding = 0.25,
      max.overlaps = 35,
      force = 2.5
    ) +
    scale_size_identity() +
    scale_alpha_continuous(range = c(0.25, 0.75), guide = "none") +
    coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) +
    labs(title = sprintf("%s. %s (n=%d)", letter, guild_name, nrow(node_data))) +
    theme_void() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5,
                                color = guild_color),
      plot.margin = margin(3, 3, 3, 3)
    )

  return(p)
}

p_B2 <- create_guild_panel_large(1, "B")
p_C2 <- create_guild_panel_large(2, "C")
p_D2 <- create_guild_panel_large(3, "D")
p_E2 <- create_guild_panel_large(4, "E")

# Alternate layout: Overview on left, 4 panels on right in 2x2
p_right <- (p_B2 | p_C2) / (p_D2 | p_E2)

p_alt <- (p_A_compact | p_right) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title = "Figure 4: CAFI Co-occurrence Network",
    subtitle = sprintf(
      "Four ecological guilds identified via Louvain community detection (Q = %.2f) | %d species | %d edges (r > 0.7)",
      modularity(communities), vcount(g), ecount(g)
    ),
    theme = theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave(
  file.path(fig_dir, "fig4_5panel_wide.png"),
  p_alt,
  width = 16,
  height = 10,
  dpi = 300,
  bg = "white"
)

cat("Saved: fig4_5panel_wide.png\n")

cat("\n============================================================\n")
cat("COMPLETE: 5-panel figures created\n")
cat("============================================================\n")
cat("  fig4_5panel.png      - Vertical layout (12x14)\n")
cat("  fig4_5panel_wide.png - Horizontal layout (16x10)\n")
cat("============================================================\n")
