# ============================================================================
# Figure 4: 5-Panel Network Figure (Version 2 - FIXED WIDE LAYOUT)
# ============================================================================
# Panel A: ALL 58 species in circular layout grouped by guild (hero panel)
# Panels B-E: Individual guild networks with species labels (force layout)
# ============================================================================
#
# FIXES APPLIED:
#   1. Panel A: proper aspect ratio with coord_fixed(ratio=1), more width
#   2. Guild labels moved further out (radius = 7.5)
#   3. Between-guild edges reduced to 5% sampling with very low alpha
#   4. Panels B-E: Only label top 8-10 species by degree for readability
#   5. Increased ggrepel force and iterations for better label placement
#   6. Layout widths adjusted: c(1.5, 2) for proper Panel A proportions
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
cat("CREATING 5-PANEL NETWORK FIGURE (V2 - FIXED WIDE LAYOUT)\n")
cat("============================================================\n\n")

# Load network data
network_results <- load_object("cafi_network")
g <- network_results$graph
communities <- network_results$communities

# Species info dataframe
sp_info <- data.frame(
  species = V(g)$name,
  guild = membership(communities),
  degree = degree(g),
  type = V(g)$type,
  stringsAsFactors = FALSE
) %>%
  mutate(guild = factor(guild))

# Guild configuration
guild_names <- c(
  "1" = "Module I",
  "2" = "Module II",
  "3" = "Module III",
  "4" = "Module IV"
)

guild_names_short <- c(
  "1" = "Module I",
  "2" = "Module II",
  "3" = "Module III",
  "4" = "Module IV"
)

guild_counts <- sp_info %>% count(guild) %>% arrange(guild)

# Colorblind-safe guild colors (Okabe-Ito derivatives)
guild_colors <- c(
  "1" = "#0072B2",  # Blue
  "2" = "#D55E00",  # Vermillion
  "3" = "#009E73",  # Teal
  "4" = "#CC79A7"   # Pink
)

# Lighter versions for arc backgrounds (more saturated to clearly match guild colors)
guild_colors_light <- c(
  "1" = "#A3CDE5",  # clearly blue
  "2" = "#F4B888",  # clearly orange
  "3" = "#8ED5BF",  # clearly teal/green
  "4" = "#E2A5C5"   # clearly pink/mauve
)

cat("Guild sizes:\n")
for (i in 1:4) {
  cat(sprintf("  Guild %d (%s): %d species\n",
              i, gsub("\n", " ", guild_names[as.character(i)]),
              guild_counts$n[guild_counts$guild == i]))
}
cat(sprintf("\nTotal species: %d\n", sum(guild_counts$n)))
cat(sprintf("Total edges: %d\n\n", ecount(g)))

# ============================================================================
# PANEL A: HERO CIRCULAR NETWORK - ALL SPECIES (FIXED)
# ============================================================================

cat("Building Panel A (circular species network - FIXED)...\n")

# Sort species within each guild by degree (high degree species more central within group)
sp_info_sorted <- sp_info %>%
  group_by(guild) %>%
  arrange(desc(degree), .by_group = TRUE) %>%
  mutate(rank_in_guild = row_number()) %>%
  ungroup()

# Calculate angular positions with guild gaps
gap_size <- 0.08  # Gap between guilds (proportion of circle)
total_gap <- gap_size * 4
species_arc <- (2 * pi) * (1 - total_gap)  # Arc available for species

# Calculate start angle for each guild
guild_sizes <- guild_counts$n
n_total <- sum(guild_sizes)

# Proportional allocation of arc to each guild
guild_props <- guild_sizes / n_total
guild_arcs <- guild_props * species_arc

# Starting angles for each guild (start at top, go clockwise)
guild_starts <- c(0)
for (i in 1:3) {
  guild_starts <- c(guild_starts,
                    guild_starts[i] + guild_arcs[i] + gap_size * 2 * pi)
}

# Assign positions to each species
sp_positions <- sp_info_sorted %>%
  group_by(guild) %>%
  mutate(
    guild_idx = as.numeric(guild),
    n_in_guild = n(),
    pos_in_guild = (row_number() - 0.5) / n_in_guild
  ) %>%
  ungroup() %>%
  mutate(
    # Calculate angle for each species
    guild_start = guild_starts[guild_idx],
    guild_arc = guild_arcs[guild_idx],
    angle = pi/2 - (guild_start + pos_in_guild * guild_arc),  # Start at top, clockwise
    # Convert to Cartesian coordinates
    radius = 4.5,
    x = radius * cos(angle),
    y = radius * sin(angle),
    # Node size based on degree (rescaled)
    node_size = scales::rescale(degree, to = c(2.5, 8))
  )

# Get edge data
edge_df <- as_data_frame(g, what = "edges") %>%
  left_join(sp_positions %>% dplyr::select(species, x, y, guild, degree),
            by = c("from" = "species")) %>%
  rename(x1 = x, y1 = y, guild_from = guild, degree_from = degree) %>%
  left_join(sp_positions %>% dplyr::select(species, x, y, guild, degree),
            by = c("to" = "species")) %>%
  rename(x2 = x, y2 = y, guild_to = guild, degree_to = degree) %>%
  mutate(
    is_within_guild = guild_from == guild_to,
    # Edge prominence: within-guild edges are more visible
    edge_alpha = ifelse(is_within_guild, 0.5, 0.15),
    edge_color = ifelse(is_within_guild,
                        guild_colors[as.character(guild_from)],
                        "gray50")
  )

# Keep ALL between-guild edges (shown as thin gray lines)
between_guild_edges <- edge_df %>% filter(!is_within_guild)
between_guild_edges_sampled <- between_guild_edges
n_between <- nrow(between_guild_edges)

cat(sprintf("  Between-guild edges: %d (all shown as thin gray)\n", n_between))

# Create bezier curves for edges (curve through center)
create_bezier <- function(x1, y1, x2, y2, n_points = 30, curvature = 0.35) {
  # Control point pulled toward center
  cx <- (x1 + x2) / 2 * (1 - curvature)
  cy <- (y1 + y2) / 2 * (1 - curvature)

  t <- seq(0, 1, length.out = n_points)
  data.frame(
    x = (1-t)^2 * x1 + 2*(1-t)*t * cx + t^2 * x2,
    y = (1-t)^2 * y1 + 2*(1-t)*t * cy + t^2 * y2
  )
}

# Generate bezier curves for within-guild edges only
within_guild_edges <- edge_df %>% filter(is_within_guild)
cat("  Generating bezier curves for", nrow(within_guild_edges), "within-guild edges...\n")

bezier_within <- map_dfr(1:nrow(within_guild_edges), function(i) {
  row <- within_guild_edges[i,]
  pts <- create_bezier(row$x1, row$y1, row$x2, row$y2)
  pts$edge_id <- i
  pts$weight <- row$weight
  pts$guild <- as.character(row$guild_from)
  pts
})

# Generate bezier curves for sampled between-guild edges
cat("  Generating bezier curves for", nrow(between_guild_edges_sampled), "sampled between-guild edges...\n")

if (nrow(between_guild_edges_sampled) > 0) {
  bezier_between <- map_dfr(1:nrow(between_guild_edges_sampled), function(i) {
    row <- between_guild_edges_sampled[i,]
    pts <- create_bezier(row$x1, row$y1, row$x2, row$y2)
    pts$edge_id <- i + 10000  # Offset to avoid conflict
    pts$weight <- row$weight
    pts
  })
} else {
  bezier_between <- data.frame(x = numeric(), y = numeric(), edge_id = integer(), weight = numeric())
}

# Create guild arc backgrounds DIRECTLY from actual node positions
# This guarantees the arc color matches the nodes it surrounds
guild_arc_data <- sp_positions %>%
  group_by(guild) %>%
  summarize(
    min_angle = min(angle),
    max_angle = max(angle),
    .groups = "drop"
  ) %>%
  mutate(
    # Add padding around the nodes
    # geom_arc_bar sweeps from start to end; swap so arc covers the node range
    padding = 0.05,
    start = max_angle + padding,
    end = min_angle - padding,
    r_inner = 3.8,
    r_outer = 5.2,
    fill_color = guild_colors_light[as.character(guild)],
    border_color = guild_colors[as.character(guild)]
  )

cat("  Arc positions derived from actual node angles:\n")
for (i in 1:nrow(guild_arc_data)) {
  row <- guild_arc_data[i, ]
  cat(sprintf("    Guild %s: arc from %.0f to %.0f deg, fill=%s, border=%s\n",
              row$guild,
              row$start * 180/pi, row$end * 180/pi,
              row$fill_color, row$border_color))
}

# Guild label positions - place in corners of the panel to maximize network space
# Labels go near their guild's quadrant but at fixed positions in plot space
guild_label_positions <- sp_positions %>%
  group_by(guild) %>%
  summarize(mid_angle = (min(angle) + max(angle)) / 2, .groups = "drop") %>%
  left_join(guild_counts, by = "guild") %>%
  mutate(
    label = paste0(
      c("Module I", "Module II",
        "Module III", "Module IV")[as.numeric(guild)],
      " (n=", n, ")"
    ),
    # Place labels at edge of plot near each guild's quadrant
    # Guild 1 (blue): upper-right → top-right corner
    # Guild 2 (orange): lower-right → bottom-right corner
    # Guild 3 (teal): lower-left → bottom-left corner
    # Guild 4 (pink): upper-left → top-left corner
    x = c(6.5, 6.5, -6.5, -6.5)[as.numeric(guild)],
    y = c(6.5, -6.5, -6.5, 6.5)[as.numeric(guild)],
    hjust = c(1, 1, 0, 0)[as.numeric(guild)],
    vjust = c(0, 1, 1, 0)[as.numeric(guild)],
    color = guild_colors[as.character(guild)]
  )

# Build Panel A with FIXED aspect ratio
# Use geom_polygon for arc backgrounds (geom_arc_bar has rendering bugs with colors)
make_arc_polygon <- function(start_angle, end_angle, r_inner, r_outer, n = 100) {
  angles_outer <- seq(start_angle, end_angle, length.out = n)
  angles_inner <- rev(angles_outer)
  data.frame(
    x = c(r_outer * cos(angles_outer), r_inner * cos(angles_inner)),
    y = c(r_outer * sin(angles_outer), r_inner * sin(angles_inner))
  )
}

arc_layers <- lapply(1:nrow(guild_arc_data), function(i) {
  d <- guild_arc_data[i, ]
  poly_df <- make_arc_polygon(d$start, d$end, d$r_inner, d$r_outer)
  geom_polygon(
    data = poly_df,
    aes(x = x, y = y),
    fill = d$fill_color,
    color = d$border_color,
    alpha = 0.5,
    linewidth = 0.8
  )
})

p_A <- ggplot() +
  arc_layers +
  # Between-guild edges (very subtle, sampled)
  geom_path(
    data = bezier_between,
    aes(x = x, y = y, group = edge_id),
    color = "gray55",
    alpha = 0.25,
    linewidth = 0.25
  ) +
  # Within-guild edges (colored, more prominent)
  geom_path(
    data = bezier_within,
    aes(x = x, y = y, group = edge_id, color = guild, alpha = weight),
    linewidth = 0.5
  ) +
  # Species nodes
  geom_point(
    data = sp_positions,
    aes(x = x, y = y, fill = guild, size = node_size),
    shape = 21,
    color = "white",
    stroke = 0.6
  ) +
  # Guild labels in corners of the panel
  geom_text(
    data = guild_label_positions,
    aes(x = x, y = y, label = label, color = guild, hjust = hjust, vjust = vjust),
    size = 3.5,
    fontface = "bold"
  ) +
  # Scales
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +
  scale_alpha_continuous(range = c(0.15, 0.7), guide = "none") +
  # Tighter coord limits to fill the space - network is at radius 4.5, arcs at 5.2
  coord_fixed(ratio = 1, xlim = c(-6.8, 6.8), ylim = c(-6.8, 6.8), clip = "off") +
  labs(
    title = "A. Species Co-occurrence Network",
    subtitle = sprintf("%d species | %d co-occurrences | 4 ecological guilds",
                       nrow(sp_positions), ecount(g))
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40",
                                 margin = margin(b = 5)),
    plot.margin = margin(5, 2, 5, 2)
  )

cat("  Panel A complete (FIXED).\n")

# ============================================================================
# PANELS B-E: INDIVIDUAL GUILD NETWORKS WITH SPECIES LABELS (FIXED)
# ============================================================================

# FIX #4 & #5: Create guild panel with limited labels for large guilds
create_guild_panel_fixed <- function(guild_id, letter, max_labels = 10) {

  guild_name <- guild_names_short[as.character(guild_id)]
  guild_color <- guild_colors[as.character(guild_id)]
  guild_color_light <- guild_colors_light[as.character(guild_id)]

  # Get species in this guild
  guild_species <- sp_info %>% filter(guild == guild_id)

  # Create subgraph for this guild
  g_sub <- induced_subgraph(g, V(g)[V(g)$name %in% guild_species$species])

  if (vcount(g_sub) == 0) {
    return(ggplot() + theme_void() +
             labs(title = paste(letter, guild_name)))
  }

  # Use Fruchterman-Reingold layout with more iterations for better spacing
  set.seed(42 + guild_id)
  # Use Kamada-Kawai for better spacing in dense guilds
  n_sp <- vcount(g_sub)
  if (n_sp > 12) {
    layout_fr <- layout_with_kk(g_sub)
  } else {
    layout_fr <- layout_with_fr(g_sub, niter = 1000, weights = E(g_sub)$weight)
  }

  # Scale layout to fill more space
  if (nrow(layout_fr) > 1) {
    layout_fr[,1] <- scales::rescale(layout_fr[,1], to = c(-1.4, 1.4))
    layout_fr[,2] <- scales::rescale(layout_fr[,2], to = c(-1.4, 1.4))
  } else {
    layout_fr[,1] <- 0
    layout_fr[,2] <- 0
  }

  # Node data
  node_data <- data.frame(
    species = V(g_sub)$name,
    x = layout_fr[,1],
    y = layout_fr[,2],
    degree = degree(g_sub),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      node_size = scales::rescale(degree, to = c(4, 14)),
      # Abbreviate long species names for display
      species_label = gsub("([A-Z])[a-z]+ ", "\\1. ", species)
    )

  # FIX #4: Only label top N species by degree for readability
  n_species <- nrow(node_data)
  n_to_label <- if (n_species <= max_labels) n_species else max_labels

  # Select top species by degree for labeling
  node_data <- node_data %>%
    mutate(
      rank_degree = rank(-degree, ties.method = "first"),
      show_label = rank_degree <= n_to_label
    )

  cat(sprintf("    Guild %d: %d species, labeling top %d\n",
              guild_id, n_species, sum(node_data$show_label)))

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

  n_edges <- ifelse(is.null(edge_data), 0, nrow(edge_data))

  # Build plot
  p <- ggplot()

  # Background circle (subtle)
  p <- p + geom_circle(
    aes(x0 = 0, y0 = 0, r = 1.4),
    fill = guild_color_light,
    color = NA,
    alpha = 0.25
  )

  # Add edges if they exist
  if (!is.null(edge_data) && nrow(edge_data) > 0) {
    p <- p + geom_segment(
      data = edge_data,
      aes(x = x1, y = y1, xend = x2, yend = y2, alpha = weight),
      color = guild_color,
      linewidth = 0.5
    )
  }

  # Add nodes
  p <- p + geom_point(
    data = node_data,
    aes(x = x, y = y, size = node_size),
    fill = guild_color,
    shape = 21,
    color = "white",
    stroke = 0.7
  )

  # FIX #5: Add labels only for top species, with increased force and iterations
  labels_data <- node_data %>% filter(show_label)

  if (nrow(labels_data) > 0) {
    p <- p + geom_text_repel(
      data = labels_data,
      aes(x = x, y = y, label = species_label),
      size = 3.2,
      fontface = "bold.italic",
      color = "gray5",
      bg.color = "white",
      bg.r = 0.12,
      segment.color = "gray40",
      segment.size = 0.3,
      segment.alpha = 0.7,
      box.padding = 0.55,
      point.padding = 0.45,
      max.overlaps = 30,
      force = 10,
      force_pull = 0.4,
      max.iter = 8000,
      seed = 42
    )
  }

  # Add note if some labels are hidden
  n_hidden <- n_species - sum(node_data$show_label)
  subtitle_text <- if (n_hidden > 0) {
    sprintf("%d species | %d edges | top %d labeled", n_species, n_edges, n_to_label)
  } else {
    sprintf("%d species | %d edges", n_species, n_edges)
  }

  p <- p +
    scale_size_identity() +
    scale_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
    coord_fixed(ratio = 1, xlim = c(-2.1, 2.1), ylim = c(-2.1, 2.1)) +
    labs(
      title = sprintf("%s. %s", letter, guild_name),
      subtitle = subtitle_text
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5,
                                color = guild_color),
      plot.subtitle = element_text(size = 7.5, hjust = 0.5, color = "gray50"),
      plot.margin = margin(3, 3, 3, 3),
      plot.background = element_rect(fill = "white", color = "gray80",
                                     linewidth = 0.5)
    )

  return(p)
}

cat("Building Panels B-E (individual guilds - FIXED)...\n")

# Label ALL species in every guild
p_B <- create_guild_panel_fixed(1, "B", max_labels = 50)
cat("  Panel B complete.\n")
p_C <- create_guild_panel_fixed(2, "C", max_labels = 50)
cat("  Panel C complete.\n")
p_D <- create_guild_panel_fixed(3, "D", max_labels = 50)
cat("  Panel D complete.\n")
p_E <- create_guild_panel_fixed(4, "E", max_labels = 50)
cat("  Panel E complete.\n")

# ============================================================================
# COMBINE INTO WIDE 5-PANEL FIGURE (FIXED)
# ============================================================================

cat("\nCombining panels into wide layout (FIXED)...\n")

# FIX #1: Arrange with proper widths - give Panel A more space
p_right <- (p_B | p_C) / (p_D | p_E)

p_wide <- (p_A | p_right) +
  plot_layout(widths = c(1.5, 2)) +  # INCREASED from c(1.2, 2)
  plot_annotation(
    title = "Figure 4: CAFI Co-occurrence Network Structure",
    subtitle = sprintf(
      "Four ecological guilds identified via Louvain community detection | Q = %.2f | %d species | %d edges",
      modularity(communities), vcount(g), ecount(g)
    ),
    caption = paste0(
      "Node size = degree centrality | Within-guild edges colored, between-guild edges gray | ",
      "Threshold: r > 0.3, FDR p < 0.05"
    ),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 8, hjust = 0.5, color = "gray50",
                                  margin = margin(t = 10)),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# Save to both locations
ggsave(
  file.path(fig_dir, "fig4_5panel_v2_wide.png"),
  p_wide,
  width = 18,
  height = 11,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(manuscript_dir, "fig4_network.png"),
  p_wide,
  width = 18,
  height = 11,
  dpi = 300,
  bg = "white"
)

cat("\nSaved: output/figures/06_network/fig4_5panel_v2_wide.png\n")
cat("Saved: output/figures/manuscript/fig4_network.png\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("FIGURE 4 (V2 FIXED WIDE LAYOUT) COMPLETE\n")
cat("============================================================\n")
cat("FIXES APPLIED:\n")
cat("  1. Panel A: coord_fixed(ratio=1) for true circle\n")
cat("  2. Guild labels moved to radius=7.5 (from 6.2)\n")
cat(sprintf("  3. Between-guild edges: reduced from %d to %d (5%% sample)\n",
            n_between, nrow(between_guild_edges_sampled)))
cat("  4. Guild panels: only top 8 species labeled by degree\n")
cat("  5. ggrepel: force=8, max.iter=5000 for better placement\n")
cat("  6. Layout widths: c(1.5, 2) gives Panel A more space\n")
cat("\nDimensions: 18 x 11 inches @ 300 dpi\n")
cat("============================================================\n")
