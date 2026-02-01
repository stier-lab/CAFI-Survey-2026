# ============================================================================
# Figure 4: 5-Panel Network - V3
# Spread out nodes more, reduce margins, better use of space
# ============================================================================

library(tidyverse)
library(igraph)
library(ggforce)
library(ggrepel)
library(patchwork)

source("scripts/00_setup.R")

fig_dir <- file.path(PATHS$figures, "06_network")
manuscript_dir <- file.path(PATHS$figures, "manuscript")

cat("============================================================\n")
cat("CREATING 5-PANEL FIGURE V3 - Expanded layout\n")
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
  stringsAsFactors = FALSE
) %>% mutate(guild = factor(guild))

guild_names <- c(
 "1" = "Protective Crustaceans",
 "2" = "Core Associates",
 "3" = "Fish & Trapezia",
 "4" = "Peripheral Crabs"
)

guild_colors <- c(
  "1" = "#0072B2",
  "2" = "#D55E00",
  "3" = "#009E73",
  "4" = "#CC79A7"
)

guild_counts <- sp_info %>% count(guild)

# ============================================================================
# PANEL A: Full circular network - all species, no labels
# ============================================================================

cat("Building Panel A...\n")

# Create circular layout with gaps between guilds
create_circular_layout <- function(sp_info, radius = 4) {
  sp_ordered <- sp_info %>% arrange(guild, desc(degree))

  guild_counts <- sp_ordered %>% count(guild)
  n_guilds <- nrow(guild_counts)

  # Larger gaps between guilds
  gap_angle <- 18 * pi / 180
  available_angle <- 2 * pi - n_guilds * gap_angle

  # Square root scaling for visual balance
  guild_counts <- guild_counts %>%
    mutate(
      weight = sqrt(n),
      arc_size = (weight / sum(weight)) * available_angle
    )

  current_angle <- pi/2
  guild_arcs <- list()

  for (i in 1:n_guilds) {
    g_id <- guild_counts$guild[i]
    arc_size <- guild_counts$arc_size[i]

    guild_arcs[[as.character(g_id)]] <- list(
      start = current_angle,
      end = current_angle - arc_size,
      mid = current_angle - arc_size/2
    )
    current_angle <- current_angle - arc_size - gap_angle
  }

  sp_ordered$x <- NA
  sp_ordered$y <- NA
  sp_ordered$angle <- NA

  for (g_id in unique(sp_ordered$guild)) {
    arc <- guild_arcs[[as.character(g_id)]]
    guild_sp <- which(sp_ordered$guild == g_id)
    n_sp <- length(guild_sp)

    padding <- (arc$start - arc$end) * 0.05
    angles <- seq(arc$start - padding, arc$end + padding, length.out = n_sp)

    sp_ordered$x[guild_sp] <- radius * cos(angles)
    sp_ordered$y[guild_sp] <- radius * sin(angles)
    sp_ordered$angle[guild_sp] <- angles
  }

  list(layout = sp_ordered, arcs = guild_arcs)
}

layout_A <- create_circular_layout(sp_info, radius = 4)
node_df_A <- layout_A$layout
guild_arcs <- layout_A$arcs

# Get edges
edge_df <- as_data_frame(g, what = "edges") %>%
  left_join(node_df_A %>% dplyr::select(species, guild, x, y),
            by = c("from" = "species")) %>%
  rename(x1 = x, y1 = y, guild_from = guild) %>%
  left_join(node_df_A %>% dplyr::select(species, guild, x, y),
            by = c("to" = "species")) %>%
  rename(x2 = x, y2 = y, guild_to = guild) %>%
  mutate(is_within = guild_from == guild_to)

# Create bezier curves
create_bezier <- function(x1, y1, x2, y2, n = 25) {
  pull <- 0.25
  cx <- (x1 + x2) / 2 * (1 - pull)
  cy <- (y1 + y2) / 2 * (1 - pull)
  t <- seq(0, 1, length.out = n)
  data.frame(
    x = (1-t)^2 * x1 + 2*(1-t)*t * cx + t^2 * x2,
    y = (1-t)^2 * y1 + 2*(1-t)*t * cy + t^2 * y2
  )
}

# Sample edges to reduce clutter - keep all within, sample between
within_edges <- edge_df %>% filter(is_within)
between_edges <- edge_df %>% filter(!is_within) %>%
  sample_frac(0.15)  # Only 15% of between-guild edges

sampled_edges <- bind_rows(within_edges, between_edges)

bezier_A <- map_dfr(1:nrow(sampled_edges), function(i) {
  row <- sampled_edges[i,]
  pts <- create_bezier(row$x1, row$y1, row$x2, row$y2)
  pts$edge_id <- i
  pts$is_within <- row$is_within
  pts$guild <- as.character(row$guild_from)
  pts$weight <- row$weight
  pts
})

# Node sizes
node_df_A <- node_df_A %>%
  mutate(node_size = scales::rescale(degree, to = c(2.5, 9)))

# Guild label positions
guild_labels <- data.frame(
  guild = names(guild_arcs),
  x = sapply(guild_arcs, function(a) 5.2 * cos(a$mid)),
  y = sapply(guild_arcs, function(a) 5.2 * sin(a$mid)),
  angle = sapply(guild_arcs, function(a) a$mid),
  stringsAsFactors = FALSE
) %>%
  mutate(
    name = guild_names[guild],
    n = guild_counts$n[match(guild, guild_counts$guild)],
    label = paste0(name, "\n(n=", n, ")")
  )

p_A <- ggplot() +
  # Guild arc backgrounds
  geom_arc(
    data = data.frame(
      guild = names(guild_arcs),
      start = sapply(guild_arcs, function(a) a$end),
      end = sapply(guild_arcs, function(a) a$start)
    ),
    aes(x0 = 0, y0 = 0, r = 4, start = start, end = end, color = guild),
    linewidth = 12, alpha = 0.15
  ) +
  # Between-guild edges
  geom_path(
    data = bezier_A %>% filter(!is_within),
    aes(x = x, y = y, group = edge_id),
    color = "gray70", linewidth = 0.15, alpha = 0.25
  ) +
  # Within-guild edges
  geom_path(
    data = bezier_A %>% filter(is_within),
    aes(x = x, y = y, group = edge_id, color = guild),
    linewidth = 0.3, alpha = 0.45
  ) +
  # Nodes
  geom_point(
    data = node_df_A,
    aes(x = x, y = y, fill = guild, size = node_size),
    shape = 21, color = "white", stroke = 0.4
  ) +
  # Guild labels
  geom_text(
    data = guild_labels,
    aes(x = x, y = y, label = label, color = guild),
    size = 3.2, fontface = "bold", lineheight = 0.85
  ) +
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +
  coord_fixed(xlim = c(-6.5, 6.5), ylim = c(-6.5, 6.5)) +
  labs(
    title = "A. Species Co-occurrence Network",
    subtitle = sprintf("%d species | %d co-occurrences | 4 ecological guilds",
                       nrow(node_df_A), ecount(g))
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
    plot.margin = margin(0, 0, 5, 0)  # Minimal margins
  )

# ============================================================================
# PANELS B-E: Individual guilds with EXPANDED layouts
# Key: Use more space, spread nodes out more
# ============================================================================

create_guild_panel_expanded <- function(guild_id, letter) {

  guild_name <- guild_names[as.character(guild_id)]
  guild_color <- guild_colors[as.character(guild_id)]

  guild_species <- sp_info %>% filter(guild == guild_id)
  g_sub <- induced_subgraph(g, V(g)[V(g)$name %in% guild_species$species])

  if (vcount(g_sub) == 0) return(ggplot() + theme_void())

  # Use Kamada-Kawai for better spacing, or FR with more iterations
  set.seed(42 + as.numeric(guild_id))

  # Try different layouts based on guild size
  n_species <- vcount(g_sub)

  if (n_species <= 6) {
    # Small guild - use circular layout
    layout_mat <- layout_in_circle(g_sub)
  } else {
    # Larger guild - use Kamada-Kawai for better spacing
    layout_mat <- layout_with_kk(g_sub)
  }

  # Scale to use FULL available space [-1, 1] with some padding
  layout_mat[,1] <- scales::rescale(layout_mat[,1], to = c(-0.95, 0.95))
  layout_mat[,2] <- scales::rescale(layout_mat[,2], to = c(-0.95, 0.95))

  node_data <- data.frame(
    species = V(g_sub)$name,
    x = layout_mat[,1],
    y = layout_mat[,2],
    degree = degree(g_sub)
  ) %>%
    mutate(
      node_size = scales::rescale(degree, to = c(4, 12)),
      # Abbreviated labels
      label = gsub("^([A-Z])[a-z]+ ", "\\1. ", species)
    )

  # Edges
  if (ecount(g_sub) > 0) {
    edge_data <- as_data_frame(g_sub, what = "edges") %>%
      left_join(node_data %>% dplyr::select(species, x, y),
                by = c("from" = "species")) %>%
      rename(x1 = x, y1 = y) %>%
      left_join(node_data %>% dplyr::select(species, x, y),
                by = c("to" = "species")) %>%
      rename(x2 = x, y2 = y)
  } else {
    edge_data <- NULL
  }

  n_edges <- ifelse(is.null(edge_data), 0, nrow(edge_data))

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
      aes(x = x, y = y, label = label),
      size = 2.6,
      fontface = "italic",
      color = "gray20",
      segment.color = "gray50",
      segment.size = 0.2,
      segment.alpha = 0.6,
      box.padding = 0.35,
      point.padding = 0.25,
      max.overlaps = 50,  # Allow more overlaps to be resolved
      force = 3,
      force_pull = 0.1,
      max.iter = 5000  # More iterations for better placement
    ) +
    scale_size_identity() +
    scale_alpha_continuous(range = c(0.2, 0.7), guide = "none") +
    # Use full coordinate space
    coord_cartesian(xlim = c(-1.3, 1.3), ylim = c(-1.15, 1.15), clip = "off") +
    labs(
      title = sprintf("%s. %s", letter, guild_name),
      subtitle = sprintf("%d species | %d edges", n_species, n_edges)
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5,
                                color = guild_color),
      plot.subtitle = element_text(size = 8, hjust = 0.5, color = "gray50"),
      plot.margin = margin(2, 2, 2, 2),  # Minimal margins
      plot.background = element_rect(fill = "gray98", color = "gray80",
                                     linewidth = 0.3)
    )

  return(p)
}

cat("Building Panels B-E with expanded layouts...\n")

p_B <- create_guild_panel_expanded(1, "B")
p_C <- create_guild_panel_expanded(2, "C")
p_D <- create_guild_panel_expanded(3, "D")
p_E <- create_guild_panel_expanded(4, "E")

# ============================================================================
# COMBINE - Use full width, minimal margins
# ============================================================================

cat("Combining panels...\n")

p_bottom <- (p_B | p_C) / (p_D | p_E) +
  plot_layout(widths = c(1, 1), heights = c(1, 1))

p_final <- p_A / p_bottom +
  plot_layout(heights = c(1.1, 2)) +
  plot_annotation(
    title = "Figure 4: CAFI Co-occurrence Network Structure",
    subtitle = "Louvain community detection reveals four ecological guilds with distinct association patterns",
    caption = "Node size = degree centrality. Edge opacity = correlation strength (Spearman r > 0.7, volume-corrected).",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 5)),
      plot.caption = element_text(size = 8, hjust = 0.5, color = "gray50",
                                  margin = margin(t = 8)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 5, 10, 5)  # Reduced left/right margins
    )
  )


# Save wide version to manuscript folder
ggsave(
  file.path(fig_dir, "fig4_5panel_v3_wide.png"),
  p_final,
  width = 16,
  height = 14,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(manuscript_dir, "fig4_network_5panel.png"),
  p_final,
  width = 16,
  height = 14,
  dpi = 300,
  bg = "white"
)

cat("Saved: fig4_5panel_v3_wide.png\n")
cat("Saved: manuscript/fig4_network_5panel.png\n")

cat("\n============================================================\n")
cat("COMPLETE\n")
cat("============================================================\n")
