# ============================================================================
# fig4_network_alternatives.R - THREE ALTERNATIVE NETWORK VISUALIZATIONS
# ============================================================================
#
# Creates three radically different approaches to visualizing CAFI guild networks:
#   OPTION A: "Constellation" - Minimal elegance with hub species only
#   OPTION B: "Four Worlds" - Each guild as separate mini-network panel
#   OPTION C: "Circular Hierarchy" - Radial dendrogram with chord-style edges
#
# Output: output/figures/06_network/
#   - fig4_option_A_constellation.png
#   - fig4_option_B_four_worlds.png
#   - fig4_option_C_circular.png
# ============================================================================

cat("\n============================================================\n")
cat("FIGURE 4 ALTERNATIVES: Three Network Visualization Approaches\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

# Load network results
network_results <- load_object("cafi_network")
centrality_df <- network_results$centrality
edge_list <- network_results$edge_list

library(igraph)
library(ggrepel)
library(ggforce)
library(grid)

# Output directory
fig_dir <- file.path(PATHS$figures, "06_network")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# COLOR PALETTE (as specified)
# ============================================================================

guild_colors <- c(
  "1" = "#2D6A4F",  # Forest Teal
  "2" = "#1D3557",  # Navy Depth
  "3" = "#E07A5F",  # Terracotta
  "4" = "#81B29A"   # Sage
)

guild_colors_light <- c(
  "1" = "#2D6A4F30",
  "2" = "#1D355730",
  "3" = "#E07A5F30",
  "4" = "#81B29A30"
)

guild_names <- c(
  "1" = "Hermit Crabs",
  "2" = "Shrimp & Worms",
  "3" = "Fish & Trapezia",
  "4" = "Peripheral Crabs"
)

# ============================================================================
# DATA PREPARATION (shared across all three visualizations)
# ============================================================================

threshold <- 0.7
edges_filtered <- edge_list %>% filter(correlation > threshold)
species_in <- unique(c(edges_filtered$sp1, edges_filtered$sp2))

# Build filtered graph
g_filtered <- graph_from_data_frame(
  edges_filtered %>%
    dplyr::select(sp1, sp2, correlation) %>%
    rename(from = sp1, to = sp2, weight = correlation),
  directed = FALSE,
  vertices = data.frame(name = species_in)
)

filtered_degree <- igraph::degree(g_filtered)

species_info <- centrality_df %>%
  filter(species %in% species_in) %>%
  dplyr::select(species, type, module, occurrence) %>%
  mutate(
    degree = filtered_degree[species],
    guild = as.character(module)
  )

guild_lookup <- setNames(species_info$guild, species_info$species)

edges_filtered <- edges_filtered %>%
  mutate(
    guild1 = guild_lookup[sp1],
    guild2 = guild_lookup[sp2],
    same_guild = guild1 == guild2
  )

cat(sprintf("Species: %d | Edges: %d (within: %d, between: %d)\n",
            length(species_in), nrow(edges_filtered),
            sum(edges_filtered$same_guild), sum(!edges_filtered$same_guild)))

# Format species names for display
format_species_name <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  if (length(parts) >= 2) {
    genus <- substr(parts[1], 1, 1)
    epithet <- parts[2]
    return(paste0(genus, ". ", epithet))
  } else {
    return(name)
  }
}

species_info$label <- sapply(species_info$species, format_species_name)

# ============================================================================
# OPTION A: "CONSTELLATION" - Minimal Elegance
# ============================================================================
# Show ONLY top 5 hub species per guild as labeled nodes
# Tiny dots for non-hub species (no labels)
# Elegant curved edges, lots of white space
# ============================================================================

cat("\n--- Creating OPTION A: Constellation ---\n")

# Identify top 5 hubs per guild
top_hubs <- species_info %>%
  group_by(guild) %>%
  slice_max(degree, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  pull(species)

species_info_A <- species_info %>%
  mutate(
    is_hub = species %in% top_hubs,
    display_label = ifelse(is_hub, label, ""),
    node_size = ifelse(is_hub, 3 + (degree / max(degree)) * 8, 0.8)
  )

# Layout: spread guilds in four quadrants with generous spacing
layout_constellation <- function(sp_info) {
  n <- nrow(sp_info)
  layout <- matrix(0, nrow = n, ncol = 2)

  # Guild positions (quadrants with lots of white space)
  guild_centers <- list(
    "1" = c(-4, 3),    # Top-left
    "2" = c(4, 3),     # Top-right
    "3" = c(-4, -3),   # Bottom-left (focal)
    "4" = c(4, -3)     # Bottom-right
  )

  for (g in unique(sp_info$guild)) {
    center <- guild_centers[[g]]
    guild_sp <- sp_info %>% filter(guild == g) %>% arrange(desc(degree))
    n_g <- nrow(guild_sp)

    if (n_g == 0) next

    # Spiral layout within guild cluster
    angles <- seq(0, 2*pi * (1 - 1/n_g), length.out = n_g)
    radii <- seq(0.3, 1.8, length.out = n_g)  # Hubs closer to center

    idx <- match(guild_sp$species, sp_info$species)
    layout[idx, 1] <- center[1] + radii * cos(angles)
    layout[idx, 2] <- center[2] + radii * sin(angles)
  }

  return(layout)
}

layout_A <- layout_constellation(species_info_A)

node_df_A <- data.frame(
  x = layout_A[, 1],
  y = layout_A[, 2],
  species = species_info_A$species,
  guild = species_info_A$guild,
  is_hub = species_info_A$is_hub,
  label = species_info_A$display_label,
  node_size = species_info_A$node_size,
  degree = species_info_A$degree
)

# Build curved edges (only for hub-to-hub or hub-to-any connections)
edge_df_A <- edges_filtered %>%
  mutate(
    x = layout_A[match(sp1, species_info_A$species), 1],
    y = layout_A[match(sp1, species_info_A$species), 2],
    xend = layout_A[match(sp2, species_info_A$species), 1],
    yend = layout_A[match(sp2, species_info_A$species), 2],
    sp1_hub = sp1 %in% top_hubs,
    sp2_hub = sp2 %in% top_hubs,
    edge_alpha = ifelse(sp1_hub | sp2_hub, 0.4, 0.08),
    edge_width = ifelse(sp1_hub & sp2_hub, 0.6, 0.2)
  )

# Create plot
p_constellation <- ggplot() +
  # Edges - thin elegant curves
  geom_curve(
    data = edge_df_A,
    aes(x = x, y = y, xend = xend, yend = yend, alpha = edge_alpha),
    color = "gray40",
    linewidth = edge_df_A$edge_width,
    curvature = 0.2,
    show.legend = FALSE
  ) +
  # Non-hub nodes (tiny dots)
  geom_point(
    data = node_df_A %>% filter(!is_hub),
    aes(x = x, y = y, color = guild),
    size = 1.2,
    alpha = 0.4
  ) +
  # Hub nodes (larger, prominent)
  geom_point(
    data = node_df_A %>% filter(is_hub),
    aes(x = x, y = y, color = guild, size = node_size),
    alpha = 0.9
  ) +
  # Hub labels
  geom_text_repel(
    data = node_df_A %>% filter(is_hub),
    aes(x = x, y = y, label = label, color = guild),
    size = 3.2,
    fontface = "italic",
    segment.color = "gray70",
    segment.size = 0.3,
    segment.alpha = 0.6,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3,
    force = 2
  ) +
  # Guild annotations
  annotate("text", x = -4, y = 5.2, label = guild_names["1"],
           color = guild_colors["1"], size = 4, fontface = "bold") +
  annotate("text", x = 4, y = 5.2, label = guild_names["2"],
           color = guild_colors["2"], size = 4, fontface = "bold") +
  annotate("text", x = -4, y = -5.2, label = guild_names["3"],
           color = guild_colors["3"], size = 4, fontface = "bold") +
  annotate("text", x = 4, y = -5.2, label = guild_names["4"],
           color = guild_colors["4"], size = 4, fontface = "bold") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_alpha_identity() +
  scale_size_identity() +
  coord_fixed(ratio = 1, xlim = c(-7, 7), ylim = c(-6, 6)) +
  labs(
    title = "CAFI Species Network: Constellation View",
    subtitle = "Top 5 hub species per guild labeled; edges show co-occurrence (r > 0.7)"
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5, margin = margin(b = 15)),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(file.path(fig_dir, "fig4_option_A_constellation.png"),
       p_constellation, width = 12, height = 10, dpi = 300, bg = "white")
cat("  Saved: fig4_option_A_constellation.png\n")


# ============================================================================
# OPTION B: "FOUR WORLDS" - Each Guild as Separate Mini-Network
# ============================================================================
# 4 separate small network panels arranged 2x2
# Each panel shows ONLY that guild's internal structure
# Between-guild connections shown as thin lines BETWEEN panels
# Equal visual real estate regardless of guild size
# ============================================================================

cat("\n--- Creating OPTION B: Four Worlds ---\n")

# Create within-guild layouts for each guild
create_guild_layout <- function(sp_info, guild_id, center_x, center_y, radius = 1.5) {
  guild_sp <- sp_info %>%
    filter(guild == guild_id) %>%
    arrange(desc(degree))

  n <- nrow(guild_sp)
  if (n == 0) return(NULL)

  # Force-directed-like layout within the guild box
  angles <- seq(0, 2*pi * (1 - 1/max(n, 2)), length.out = n)

  # Hub species closer to center
  radii <- sqrt(seq(1, n) / n) * radius

  data.frame(
    species = guild_sp$species,
    x = center_x + radii * cos(angles),
    y = center_y + radii * sin(angles),
    guild = guild_id,
    degree = guild_sp$degree,
    label = guild_sp$label
  )
}

# Panel centers (2x2 grid with space between)
panel_centers <- list(
  "1" = c(-3.5, 3.5),   # Top-left
  "2" = c(3.5, 3.5),    # Top-right
  "3" = c(-3.5, -3.5),  # Bottom-left
  "4" = c(3.5, -3.5)    # Bottom-right
)

# Create all node positions
node_list_B <- lapply(names(panel_centers), function(g) {
  create_guild_layout(species_info, g, panel_centers[[g]][1], panel_centers[[g]][2])
})
node_df_B <- bind_rows(node_list_B)

# Match species_info order
node_df_B <- node_df_B %>%
  mutate(
    node_size = 2 + (degree / max(degree)) * 6,
    is_top5 = species %in% (species_info %>%
                              group_by(guild) %>%
                              slice_max(degree, n = 5) %>%
                              pull(species))
  )

# Create lookup for positions
pos_lookup_x <- setNames(node_df_B$x, node_df_B$species)
pos_lookup_y <- setNames(node_df_B$y, node_df_B$species)

# Within-guild edges
within_edges_B <- edges_filtered %>%
  filter(same_guild) %>%
  mutate(
    x = pos_lookup_x[sp1],
    y = pos_lookup_y[sp1],
    xend = pos_lookup_x[sp2],
    yend = pos_lookup_y[sp2],
    guild = guild1
  ) %>%
  filter(!is.na(x) & !is.na(xend))

# Between-guild edges (will be drawn with different style)
between_edges_B <- edges_filtered %>%
  filter(!same_guild) %>%
  mutate(
    x = pos_lookup_x[sp1],
    y = pos_lookup_y[sp1],
    xend = pos_lookup_x[sp2],
    yend = pos_lookup_y[sp2]
  ) %>%
  filter(!is.na(x) & !is.na(xend))

# Create panel backgrounds
panel_boxes <- data.frame(
  guild = c("1", "2", "3", "4"),
  xmin = c(-6, 1, -6, 1),
  xmax = c(-1, 6, -1, 6),
  ymin = c(1, 1, -6, -6),
  ymax = c(6, 6, -1, -1)
)

# Guild labels for panels
panel_labels <- data.frame(
  guild = c("1", "2", "3", "4"),
  x = c(-3.5, 3.5, -3.5, 3.5),
  y = c(5.7, 5.7, -0.8, -0.8),
  label = c(guild_names["1"], guild_names["2"], guild_names["3"], guild_names["4"])
)

p_four_worlds <- ggplot() +
  # Panel backgrounds
  geom_rect(
    data = panel_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = guild),
    alpha = 0.08,
    color = NA
  ) +
  # Panel borders
  geom_rect(
    data = panel_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = guild),
    fill = NA,
    linewidth = 1.2
  ) +
  # Between-guild edges (thin, crossing panels)
  geom_segment(
    data = between_edges_B,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "gray60",
    linewidth = 0.2,
    alpha = 0.3,
    linetype = "dotted"
  ) +
  # Within-guild edges
  geom_segment(
    data = within_edges_B,
    aes(x = x, y = y, xend = xend, yend = yend, color = guild),
    linewidth = 0.5,
    alpha = 0.5
  ) +
  # Nodes
  geom_point(
    data = node_df_B,
    aes(x = x, y = y, color = guild, size = node_size),
    alpha = 0.85
  ) +
  # Labels for top species only
  geom_text_repel(
    data = node_df_B %>% filter(is_top5),
    aes(x = x, y = y, label = label, color = guild),
    size = 2.8,
    fontface = "italic",
    segment.size = 0.2,
    segment.alpha = 0.5,
    max.overlaps = 15,
    box.padding = 0.25,
    point.padding = 0.2
  ) +
  # Panel titles
  geom_text(
    data = panel_labels,
    aes(x = x, y = y, label = label, color = guild),
    size = 4.5,
    fontface = "bold"
  ) +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +
  coord_fixed(ratio = 1, xlim = c(-7, 7), ylim = c(-7, 7)) +
  labs(
    title = "CAFI Species Network: Four Worlds",
    subtitle = "Each guild shown as separate community; dotted lines show between-guild co-occurrence"
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5, margin = margin(b = 15)),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(file.path(fig_dir, "fig4_option_B_four_worlds.png"),
       p_four_worlds, width = 12, height = 10, dpi = 300, bg = "white")
cat("  Saved: fig4_option_B_four_worlds.png\n")


# ============================================================================
# OPTION C: "CIRCULAR HIERARCHY" - Radial Dendrogram with Chord-Style Edges
# ============================================================================
# Species arranged in a circle, grouped by guild (4 arcs)
# Edges drawn as bezier curves through the CENTER
# Guild labels as arcs around the outside
# Like a chord diagram but for networks
# ============================================================================

cat("\n--- Creating OPTION C: Circular Hierarchy ---\n")

# Arrange species in circle, grouped by guild
create_circular_layout <- function(sp_info, radius = 4) {
  # Order by guild, then by degree within guild
  sp_ordered <- sp_info %>%
    arrange(guild, desc(degree))

  n <- nrow(sp_ordered)

  # Assign arc segments to each guild with small gaps
  guild_counts <- sp_ordered %>% count(guild) %>% arrange(guild)
  total_sp <- sum(guild_counts$n)
  gap_angle <- 8 * pi / 180  # 8 degree gaps between guilds
  available_angle <- 2 * pi - 4 * gap_angle

  # Calculate start angle for each guild
  guild_arcs <- guild_counts %>%
    mutate(
      arc_size = (n / total_sp) * available_angle,
      arc_start = cumsum(lag(arc_size, default = 0)) + (row_number() - 1) * gap_angle,
      arc_end = arc_start + arc_size,
      arc_mid = (arc_start + arc_end) / 2
    )

  # Assign positions
  layout_df <- sp_ordered
  layout_df$angle <- NA
  layout_df$x <- NA
  layout_df$y <- NA

  for (g in guild_arcs$guild) {
    arc_info <- guild_arcs %>% filter(guild == g)
    guild_sp <- sp_ordered %>% filter(guild == g)
    n_g <- nrow(guild_sp)

    if (n_g == 1) {
      angles <- arc_info$arc_mid
    } else {
      angles <- seq(arc_info$arc_start + 0.02, arc_info$arc_end - 0.02, length.out = n_g)
    }

    idx <- which(layout_df$guild == g)
    layout_df$angle[idx] <- angles
    layout_df$x[idx] <- radius * cos(angles - pi/2)  # Rotate so 0 is at top
    layout_df$y[idx] <- radius * sin(angles - pi/2)
  }

  layout_df$guild_arc <- guild_arcs$arc_mid[match(layout_df$guild, guild_arcs$guild)]

  return(list(nodes = layout_df, arcs = guild_arcs))
}

circular_data <- create_circular_layout(species_info, radius = 4)
node_df_C <- circular_data$nodes
guild_arcs <- circular_data$arcs

# Add display attributes
node_df_C <- node_df_C %>%
  mutate(
    node_size = 2 + (degree / max(degree)) * 5,
    is_top3 = species %in% (species_info %>%
                              group_by(guild) %>%
                              slice_max(degree, n = 3) %>%
                              pull(species)),
    # Label angle - text reads outward
    label_angle = (angle - pi/2) * 180 / pi,
    # Flip text on left side so it's readable
    label_angle_adj = ifelse(label_angle > 90 | label_angle < -90,
                             label_angle + 180, label_angle),
    hjust = ifelse(label_angle > 90 | label_angle < -90, 1, 0)
  )

# Position lookup
pos_lookup_x_C <- setNames(node_df_C$x, node_df_C$species)
pos_lookup_y_C <- setNames(node_df_C$y, node_df_C$species)

# Create bezier curves through center for edges
create_bezier_edge <- function(x1, y1, x2, y2, n_points = 50) {
  # Control point at center (0, 0) with some pull
  # Use quadratic bezier
  t <- seq(0, 1, length.out = n_points)

  # Pull toward center more for longer edges
  dist <- sqrt((x2-x1)^2 + (y2-y1)^2)
  pull <- 0.1 + 0.15 * (dist / 8)  # Adjust pull based on distance

  cx <- (x1 + x2) / 2 * (1 - pull)
  cy <- (y1 + y2) / 2 * (1 - pull)

  # Quadratic bezier
  x <- (1-t)^2 * x1 + 2*(1-t)*t * cx + t^2 * x2
  y <- (1-t)^2 * y1 + 2*(1-t)*t * cy + t^2 * y2

  data.frame(x = x, y = y)
}

# Generate bezier paths for all edges
edge_beziers <- edges_filtered %>%
  rowwise() %>%
  mutate(
    x1 = pos_lookup_x_C[sp1],
    y1 = pos_lookup_y_C[sp1],
    x2 = pos_lookup_x_C[sp2],
    y2 = pos_lookup_y_C[sp2]
  ) %>%
  filter(!is.na(x1) & !is.na(x2)) %>%
  ungroup()

# Create bezier curves grouped by edge
bezier_list <- vector("list", nrow(edge_beziers))
for (i in seq_len(nrow(edge_beziers))) {
  bezier_pts <- create_bezier_edge(
    edge_beziers$x1[i], edge_beziers$y1[i],
    edge_beziers$x2[i], edge_beziers$y2[i]
  )
  bezier_pts$edge_id <- i
  bezier_pts$same_guild <- edge_beziers$same_guild[i]
  bezier_pts$guild1 <- edge_beziers$guild1[i]
  bezier_pts$correlation <- edge_beziers$correlation[i]
  bezier_list[[i]] <- bezier_pts
}
bezier_df <- bind_rows(bezier_list)

# Create guild arc labels (positioned outside the circle)
label_radius <- 5.3
guild_labels_C <- guild_arcs %>%
  mutate(
    x = label_radius * cos(arc_mid - pi/2),
    y = label_radius * sin(arc_mid - pi/2),
    label = guild_names[guild],
    angle = (arc_mid - pi/2) * 180 / pi,
    angle_adj = ifelse(angle > 90 | angle < -90, angle + 180, angle)
  )

# Create arc segments for guild boundaries
arc_segments <- guild_arcs %>%
  rowwise() %>%
  mutate(
    arc_points = list({
      angles <- seq(arc_start - pi/2, arc_end - pi/2, length.out = 50)
      data.frame(
        x = 4.6 * cos(angles),
        y = 4.6 * sin(angles)
      )
    })
  ) %>%
  unnest(arc_points)

p_circular <- ggplot() +
  # Guild arcs (outer ring)
  geom_path(
    data = arc_segments,
    aes(x = x, y = y, color = guild, group = guild),
    linewidth = 8,
    alpha = 0.3,
    lineend = "round"
  ) +
  # Bezier edges through center
  geom_path(
    data = bezier_df,
    aes(x = x, y = y, group = edge_id,
        color = ifelse(same_guild, guild1, "between")),
    alpha = 0.25,
    linewidth = 0.4
  ) +
  # Nodes on the circle
  geom_point(
    data = node_df_C,
    aes(x = x, y = y, color = guild, size = node_size),
    alpha = 0.9
  ) +
  # Labels for top species (positioned outside)
  geom_text(
    data = node_df_C %>% filter(is_top3),
    aes(x = x * 1.15, y = y * 1.15, label = label,
        angle = label_angle_adj, hjust = hjust, color = guild),
    size = 2.8,
    fontface = "italic"
  ) +
  # Guild labels
  geom_label(
    data = guild_labels_C,
    aes(x = x, y = y, label = label, color = guild),
    size = 4,
    fontface = "bold",
    fill = "white",
    alpha = 0.9,
    label.padding = unit(0.25, "lines"),
    label.size = 0
  ) +
  scale_color_manual(
    values = c(guild_colors, "between" = "gray50"),
    guide = "none"
  ) +
  scale_size_identity() +
  coord_fixed(ratio = 1, xlim = c(-6.5, 6.5), ylim = c(-6.5, 6.5)) +
  labs(
    title = "CAFI Species Network: Circular Hierarchy",
    subtitle = "Species arranged by guild; edges curve through center showing co-occurrence structure"
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5, margin = margin(b = 15)),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(file.path(fig_dir, "fig4_option_C_circular.png"),
       p_circular, width = 12, height = 10, dpi = 300, bg = "white")
cat("  Saved: fig4_option_C_circular.png\n")


# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("THREE NETWORK ALTERNATIVES COMPLETE\n")
cat("============================================================\n\n")
cat("Output saved to:", fig_dir, "\n\n")
cat("OPTION A: Constellation - Minimal elegance with hub species\n")
cat("OPTION B: Four Worlds - Each guild as separate panel\n")
cat("OPTION C: Circular - Radial chord-style layout\n\n")
cat("All figures: 12x10 inches, 300 DPI\n")
