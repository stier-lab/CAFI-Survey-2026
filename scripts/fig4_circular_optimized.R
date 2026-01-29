# ============================================================================
# Figure 4: Optimized Circular Guild Network
# Fixes: label overlaps, arc balance, edge clutter, readability
# ============================================================================

library(tidyverse)
library(igraph)
library(ggforce)
library(ggrepel)

# Setup paths
source("scripts/00_setup.R")

# Output directory
fig_dir <- file.path(PATHS$figures, "06_network")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading network data...\n")

# Load network results
network_results <- tryCatch(
  load_object("cafi_network"),
  error = function(e) {
    stop("Network results not found. Run 06_network_analysis.R first.")
  }
)

g <- network_results$graph
communities <- network_results$communities

# Get species info
sp_info <- data.frame(
  species = V(g)$name,
  guild = membership(communities),
  degree = degree(g),
  type = V(g)$type,
  stringsAsFactors = FALSE
) %>%
  mutate(guild = factor(guild))

# Guild naming with counts
guild_counts <- sp_info %>% count(guild) %>% arrange(guild)
guild_names <- c(
  "1" = "Protective Crustaceans",
  "2" = "Core Associates",
  "3" = "Fish & Trapezia",
  "4" = "Peripheral Crabs"
)
guild_counts$name <- guild_names[as.character(guild_counts$guild)]
guild_counts$label <- paste0(guild_counts$name, "\n(n = ", guild_counts$n, ")")

cat("\nGuild composition:\n")
print(guild_counts)

# Guild colors - high contrast, colorblind safe
guild_colors <- c(

  "1" = "#0072B2",  # Blue - Protective Crustaceans
  "2" = "#D55E00",  # Vermillion - Core Associates
  "3" = "#009E73",  # Teal - Fish & Trapezia
  "4" = "#CC79A7"   # Pink - Peripheral Crabs
)

# ============================================================================
# OPTIMIZED CIRCULAR LAYOUT
# Key fix: Give small guilds MORE arc space for visual balance
# ============================================================================

create_balanced_circular_layout <- function(sp_info, radius = 4) {

  sp_ordered <- sp_info %>%
    arrange(guild, desc(degree))

  guild_counts <- sp_ordered %>%
    count(guild) %>%
    arrange(guild)

  n_guilds <- nrow(guild_counts)
  gap_angle <- 12 * pi / 180  # 12 degree gaps between guilds
  available_angle <- 2 * pi - n_guilds * gap_angle

  # KEY FIX: Use SQUARE ROOT of count for arc size

# This gives small guilds proportionally more space
  # Guild with 4 species gets ~50% of arc of guild with 17 (instead of ~24%)
  guild_counts <- guild_counts %>%
    mutate(
      weight = sqrt(n),  # Square root for balance
      arc_size = (weight / sum(weight)) * available_angle
    )

  # Calculate arc start/end positions
  # Start at top (90 degrees) and go clockwise
  current_angle <- pi/2
  guild_arcs <- list()

  for (i in 1:n_guilds) {
    g_id <- guild_counts$guild[i]
    arc_size <- guild_counts$arc_size[i]

    guild_arcs[[as.character(g_id)]] <- list(
      start = current_angle,
      end = current_angle - arc_size,
      mid = current_angle - arc_size/2,
      n_species = guild_counts$n[i]
    )

    current_angle <- current_angle - arc_size - gap_angle
  }

  # Position species within their guild's arc
  sp_ordered$x <- NA
  sp_ordered$y <- NA
  sp_ordered$angle <- NA

  for (g_id in unique(sp_ordered$guild)) {
    arc <- guild_arcs[[as.character(g_id)]]
    guild_sp <- which(sp_ordered$guild == g_id)
    n_sp <- length(guild_sp)

    if (n_sp == 1) {
      angles <- arc$mid
    } else {
      # Distribute species evenly within arc with padding
      padding <- (arc$start - arc$end) * 0.08
      angles <- seq(arc$start - padding, arc$end + padding, length.out = n_sp)
    }

    sp_ordered$x[guild_sp] <- radius * cos(angles)
    sp_ordered$y[guild_sp] <- radius * sin(angles)
    sp_ordered$angle[guild_sp] <- angles
  }

  list(
    layout = sp_ordered,
    arcs = guild_arcs,
    guild_counts = guild_counts
  )
}

# Create layout
layout_result <- create_balanced_circular_layout(sp_info, radius = 4)
node_df <- layout_result$layout
guild_arcs <- layout_result$arcs

# ============================================================================
# EDGE FILTERING - Reduce clutter
# Show only strongest edges OR within-guild edges
# ============================================================================

# Get edge list with weights
edge_list <- as_data_frame(g, what = "edges")

# Add source/target guild info
edge_list <- edge_list %>%
  left_join(node_df %>% dplyr::select(species, guild, x, y),
            by = c("from" = "species")) %>%
  rename(x1 = x, y1 = y, guild_from = guild) %>%
  left_join(node_df %>% dplyr::select(species, guild, x, y),
            by = c("to" = "species")) %>%
  rename(x2 = x, y2 = y, guild_to = guild)

# Classify edges
edge_list <- edge_list %>%
  mutate(
    is_within_guild = guild_from == guild_to,
    edge_type = ifelse(is_within_guild, "within", "between")
  )

cat("\nEdge breakdown:\n")
cat("  Within-guild edges:", sum(edge_list$is_within_guild), "\n")
cat("  Between-guild edges:", sum(!edge_list$is_within_guild), "\n")

# FILTER: Keep all within-guild + top 30% of between-guild edges
within_edges <- edge_list %>% filter(is_within_guild)
between_edges <- edge_list %>%
  filter(!is_within_guild) %>%
  arrange(desc(weight)) %>%
  slice_head(prop = 0.3)

filtered_edges <- bind_rows(within_edges, between_edges)

cat("  After filtering:", nrow(filtered_edges), "edges shown\n")

# ============================================================================
# CREATE CURVED EDGES (Bezier through center)
# ============================================================================

create_bezier_points <- function(x1, y1, x2, y2, n_points = 30) {
  # Calculate distance for curve pull strength
  dist <- sqrt((x2 - x1)^2 + (y2 - y1)^2)

  # Pull toward center proportional to distance
  pull <- 0.15 + 0.2 * (dist / 8)

  # Control point near center
  cx <- (x1 + x2) / 2 * (1 - pull)
  cy <- (y1 + y2) / 2 * (1 - pull)

  # Quadratic bezier
  t <- seq(0, 1, length.out = n_points)
  x <- (1-t)^2 * x1 + 2*(1-t)*t * cx + t^2 * x2
  y <- (1-t)^2 * y1 + 2*(1-t)*t * cy + t^2 * y2

  data.frame(x = x, y = y)
}

# Generate bezier curves for all filtered edges
bezier_edges <- map2_dfr(
  seq_len(nrow(filtered_edges)),
  filtered_edges$weight,
  function(i, w) {
    row <- filtered_edges[i, ]
    pts <- create_bezier_points(row$x1, row$y1, row$x2, row$y2)
    pts$edge_id <- i
    pts$weight <- w
    pts$edge_type <- row$edge_type
    pts$guild <- row$guild_from  # For within-guild coloring
    pts
  }
)

# ============================================================================
# LABEL POSITIONING - Prevent overlaps
# ============================================================================

# Calculate smart label positions
node_df <- node_df %>%
  mutate(
    # Radial label position (further out for labels)
    label_radius = 4.8,
    label_x = label_radius * cos(angle),
    label_y = label_radius * sin(angle),

    # Angle for text (readable orientation)
    angle_deg = angle * 180 / pi,

    # Flip text on left side so it reads left-to-right
    label_angle = case_when(
      angle_deg > 90 | angle_deg < -90 ~ angle_deg + 180,
      TRUE ~ angle_deg
    ),

    # Horizontal justification based on position
    hjust = case_when(
      angle_deg > 90 | angle_deg < -90 ~ 1,  # Right-align on left side
      TRUE ~ 0  # Left-align on right side
    ),

    # Size nodes by degree (normalized)
    node_size = scales::rescale(degree, to = c(3, 10))
  )

# Select which species to label (top by degree, balanced across guilds)
# Label top 5 from large guilds, all from small guilds
labels_to_show <- node_df %>%
  group_by(guild) %>%
  arrange(desc(degree)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  mutate(
    n_in_guild = ave(as.numeric(guild), guild, FUN = length),
    show_label = case_when(
      n_in_guild <= 6 ~ TRUE,  # Show all in small guilds
      rank <= 5 ~ TRUE,        # Top 5 in large guilds
      TRUE ~ FALSE
    )
  )

cat("\nLabeling", sum(labels_to_show$show_label), "species\n")

# ============================================================================
# GUILD ARC PATHS
# ============================================================================

# Create arc paths for guild backgrounds
arc_paths <- map_dfr(names(guild_arcs), function(g_id) {
  arc <- guild_arcs[[g_id]]
  n_pts <- 50
  angles <- seq(arc$start, arc$end, length.out = n_pts)

  # Inner and outer arc for ribbon
  inner_r <- 3.4
  outer_r <- 4.6

  data.frame(
    guild = g_id,
    x_inner = inner_r * cos(angles),
    y_inner = inner_r * sin(angles),
    x_outer = outer_r * cos(angles),
    y_outer = outer_r * sin(angles),
    angle = angles
  )
})

# Create guild label positions (at arc midpoint, outside)
guild_label_df <- map_dfr(names(guild_arcs), function(g_id) {
  arc <- guild_arcs[[g_id]]
  mid_angle <- arc$mid

  data.frame(
    guild = g_id,
    x = 5.4 * cos(mid_angle),
    y = 5.4 * sin(mid_angle),
    angle = mid_angle,
    angle_deg = mid_angle * 180 / pi,
    label = guild_counts$label[guild_counts$guild == g_id],
    n = arc$n_species
  )
}) %>%
  mutate(
    # Adjust label angle for readability
    label_angle = case_when(
      angle_deg > 90 | angle_deg < -90 ~ angle_deg + 180,
      TRUE ~ angle_deg
    ),
    hjust = case_when(
      angle_deg > 90 | angle_deg < -90 ~ 1,
      TRUE ~ 0
    )
  )

# ============================================================================
# BUILD THE FIGURE
# ============================================================================

cat("\nBuilding optimized circular figure...\n")

p_optimized <- ggplot() +

  # 1. Guild arc backgrounds (using geom_arc for simplicity)
  geom_arc(
    data = data.frame(
      guild = names(guild_arcs),
      start = sapply(guild_arcs, function(x) x$end),
      end = sapply(guild_arcs, function(x) x$start)
    ),
    aes(x0 = 0, y0 = 0, r = 4, start = start, end = end, color = guild),
    linewidth = 10,
    alpha = 0.2
  ) +

  # 2. Between-guild edges (gray, subtle)
  geom_path(
    data = bezier_edges %>% filter(edge_type == "between"),
    aes(x = x, y = y, group = edge_id, alpha = weight),
    color = "gray60",
    linewidth = 0.3
  ) +

  # 3. Within-guild edges (colored, more prominent)
  geom_path(
    data = bezier_edges %>% filter(edge_type == "within"),
    aes(x = x, y = y, group = edge_id, color = guild, alpha = weight),
    linewidth = 0.5
  ) +

  # 4. Nodes
  geom_point(
    data = node_df,
    aes(x = x, y = y, fill = guild, size = node_size),
    shape = 21,
    color = "white",
    stroke = 0.5
  ) +

  # 5. Species labels (only selected ones)
  geom_text(
    data = labels_to_show %>% filter(show_label),
    aes(x = label_x, y = label_y, label = species,
        angle = label_angle, hjust = hjust, color = guild),
    size = 2.5,
    fontface = "italic"
  ) +

  # 6. Guild labels (bold, outside arcs)
  geom_text(
    data = guild_label_df,
    aes(x = x, y = y, label = label, angle = label_angle, hjust = hjust),
    size = 3.5,
    fontface = "bold",
    color = "gray20",
    lineheight = 0.85
  ) +

  # Scales
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +
  scale_alpha_continuous(range = c(0.15, 0.6), guide = "none") +

  # Coordinate system
  coord_fixed(xlim = c(-7, 7), ylim = c(-7, 7)) +

  # Labels
  labs(
    title = "CAFI Co-occurrence Network: Four Ecological Guilds",
    subtitle = sprintf(
      "%d species | %d significant co-occurrences (r > 0.7) | Louvain modularity Q = %.2f",
      nrow(node_df),
      nrow(filtered_edges),
      modularity(communities)
    ),
    caption = paste0(
      "Node size = degree centrality. ",
      "Curved edges connect species with volume-corrected Spearman r > 0.7.\n",
      "Within-guild edges colored; between-guild edges gray. ",
      "Arc width reflects guild size (square-root scaled for visual balance)."
    )
  ) +

  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40",
                                 margin = margin(b = 10)),
    plot.caption = element_text(size = 8, hjust = 0.5, color = "gray50",
                                margin = margin(t = 15), lineheight = 1.2),
    plot.margin = margin(20, 20, 20, 20),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Save
ggsave(
  file.path(fig_dir, "fig4_circular_v2.png"),
  p_optimized,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)

cat("\nSaved: fig4_circular_v2.png\n")

# ============================================================================
# VERSION 2: Even cleaner with ggrepel for labels
# ============================================================================

cat("\nBuilding version with ggrepel labels...\n")

p_repel <- ggplot() +

  # Guild arc backgrounds
  geom_arc(
    data = data.frame(
      guild = names(guild_arcs),
      start = sapply(guild_arcs, function(x) x$end),
      end = sapply(guild_arcs, function(x) x$start)
    ),
    aes(x0 = 0, y0 = 0, r = 4, start = start, end = end, color = guild),
    linewidth = 8,
    alpha = 0.3
  ) +

  # Between-guild edges
  geom_path(
    data = bezier_edges %>% filter(edge_type == "between"),
    aes(x = x, y = y, group = edge_id),
    color = "gray70",
    linewidth = 0.2,
    alpha = 0.3
  ) +

  # Within-guild edges
  geom_path(
    data = bezier_edges %>% filter(edge_type == "within"),
    aes(x = x, y = y, group = edge_id, color = guild),
    linewidth = 0.4,
    alpha = 0.5
  ) +

  # Nodes
  geom_point(
    data = node_df,
    aes(x = x, y = y, fill = guild, size = node_size),
    shape = 21,
    color = "white",
    stroke = 0.6
  ) +

  # Species labels with ggrepel
  geom_text_repel(
    data = labels_to_show %>% filter(show_label),
    aes(x = x, y = y, label = species, color = guild),
    size = 2.4,
    fontface = "italic",
    segment.color = "gray70",
    segment.size = 0.3,
    segment.alpha = 0.5,
    min.segment.length = 0.3,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = 20,
    force = 2,
    force_pull = 0.5
  ) +

  # Guild labels in corners
  annotate(
    "text",
    x = c(-5.5, 5.5, -5.5, 5.5),
    y = c(5.5, 5.5, -5.5, -5.5),
    label = guild_counts$label,
    hjust = c(0, 1, 0, 1),
    vjust = c(1, 1, 0, 0),
    size = 3.2,
    fontface = "bold",
    color = as.character(guild_colors[as.character(guild_counts$guild)]),
    lineheight = 0.9
  ) +

  # Scales
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +

  coord_fixed(xlim = c(-7, 7), ylim = c(-7, 7)) +

  labs(
    title = "CAFI Co-occurrence Network",
    subtitle = sprintf(
      "Four ecological guilds | %d species | Modularity Q = %.2f (z = %.1f vs null)",
      nrow(node_df),
      modularity(communities),
      network_results$null_comparison$z_score[
        network_results$null_comparison$metric == "Modularity"]
    )
  ) +

  theme_void() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40",
                                 margin = margin(b = 5)),
    plot.margin = margin(15, 15, 15, 15),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  file.path(fig_dir, "fig4_circular_v3_repel.png"),
  p_repel,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)

cat("Saved: fig4_circular_v3_repel.png\n")

# ============================================================================
# VERSION 3: Maximum clarity - separate guild panels
# ============================================================================

cat("\nBuilding multi-panel version...\n")

# Create individual guild plots
guild_plots <- map(levels(node_df$guild), function(g_id) {

  # Filter to this guild's nodes and edges
  g_nodes <- node_df %>% filter(guild == g_id)
  g_edges_within <- filtered_edges %>%
    filter(guild_from == g_id & guild_to == g_id)

  # Get edges to other guilds
  g_edges_between <- filtered_edges %>%
    filter((guild_from == g_id | guild_to == g_id) &
             guild_from != guild_to)

  # Create local bezier curves
  local_bezier <- NULL
  if (nrow(g_edges_within) > 0) {
    local_bezier <- map2_dfr(
      seq_len(nrow(g_edges_within)),
      g_edges_within$weight,
      function(i, w) {
        row <- g_edges_within[i, ]
        pts <- create_bezier_points(row$x1, row$y1, row$x2, row$y2, n_points = 20)
        pts$edge_id <- i
        pts$weight <- w
        pts
      }
    )
  }

  guild_name <- guild_names[g_id]
  guild_n <- nrow(g_nodes)

  p <- ggplot() +
    # Edges within guild
    {if (!is.null(local_bezier) && nrow(local_bezier) > 0)
      geom_path(
        data = local_bezier,
        aes(x = x, y = y, group = edge_id, alpha = weight),
        color = guild_colors[g_id],
        linewidth = 0.6
      )
    } +
    # Nodes
    geom_point(
      data = g_nodes,
      aes(x = x, y = y, size = node_size),
      fill = guild_colors[g_id],
      shape = 21,
      color = "white",
      stroke = 0.5
    ) +
    # Labels
    geom_text_repel(
      data = g_nodes,
      aes(x = x, y = y, label = species),
      size = 2.8,
      fontface = "italic",
      color = "gray20",
      segment.color = "gray60",
      segment.size = 0.2,
      max.overlaps = 15,
      box.padding = 0.5
    ) +
    scale_size_identity() +
    scale_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
    coord_fixed() +
    labs(
      title = guild_name,
      subtitle = sprintf("n = %d species", guild_n)
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5,
                                color = guild_colors[g_id]),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray50"),
      plot.margin = margin(10, 10, 10, 10)
    )

  p
})

# Combine into 2x2 grid
library(patchwork)

p_multipanel <- (guild_plots[[1]] | guild_plots[[2]]) /
  (guild_plots[[3]] | guild_plots[[4]]) +
  plot_annotation(
    title = "CAFI Ecological Guilds: Species Composition",
    subtitle = sprintf(
      "%d total species across 4 guilds | Volume-corrected correlations r > 0.7",
      nrow(node_df)
    ),
    caption = "Node size = degree centrality (number of significant co-occurrences)",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 8, hjust = 0.5, color = "gray50")
    )
  )

ggsave(
  file.path(fig_dir, "fig4_circular_v4_panels.png"),
  p_multipanel,
  width = 12,
  height = 12,
  dpi = 300,
  bg = "white"
)

cat("Saved: fig4_circular_v4_panels.png\n")

# ============================================================================
# VERSION 4: Hybrid - Full network + panel insets
# ============================================================================

cat("\nBuilding hybrid version with insets...\n")

# Main circular plot (simplified)
p_main <- ggplot() +
  # Arc backgrounds
  geom_arc(
    data = data.frame(
      guild = names(guild_arcs),
      start = sapply(guild_arcs, function(x) x$end),
      end = sapply(guild_arcs, function(x) x$start)
    ),
    aes(x0 = 0, y0 = 0, r = 4, start = start, end = end, color = guild),
    linewidth = 12,
    alpha = 0.25
  ) +
  # All edges (subtle)
  geom_path(
    data = bezier_edges,
    aes(x = x, y = y, group = edge_id, color = guild, alpha = weight),
    linewidth = 0.3
  ) +
  # Nodes
  geom_point(
    data = node_df,
    aes(x = x, y = y, fill = guild, size = node_size),
    shape = 21,
    color = "white",
    stroke = 0.4
  ) +
  # Only label top 3 per guild
  geom_text_repel(
    data = node_df %>%
      group_by(guild) %>%
      slice_max(degree, n = 3) %>%
      ungroup(),
    aes(x = x, y = y, label = species, color = guild),
    size = 2.5,
    fontface = "italic",
    segment.size = 0.2,
    segment.alpha = 0.5,
    box.padding = 0.3,
    max.overlaps = 12
  ) +
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +
  scale_alpha_continuous(range = c(0.1, 0.5), guide = "none") +
  coord_fixed(xlim = c(-6, 6), ylim = c(-6, 6)) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

# Create legend/key panel
legend_data <- guild_counts %>%
  mutate(
    y = 4:1,
    color = guild_colors[as.character(guild)]
  )

p_legend <- ggplot(legend_data) +
  geom_point(aes(x = 0, y = y, fill = guild),
             size = 6, shape = 21, color = "white", stroke = 0.5) +
  geom_text(aes(x = 0.3, y = y,
                label = paste0(name, " (", n, ")")),
            hjust = 0, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = guild_colors, guide = "none") +
  xlim(-0.5, 4) +
  ylim(0, 5) +
  theme_void()

# Combine
p_hybrid <- p_main +
  inset_element(p_legend, left = 0.02, bottom = 0.02, right = 0.45, top = 0.25) +
  plot_annotation(
    title = "CAFI Co-occurrence Network: Four Ecological Guilds",
    subtitle = sprintf(
      "%d species | %d edges (r > 0.7) | Louvain modularity Q = %.2f",
      nrow(node_df), nrow(filtered_edges), modularity(communities)
    ),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 10))
    )
  )

ggsave(
  file.path(fig_dir, "fig4_circular_v5_hybrid.png"),
  p_hybrid,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)

cat("Saved: fig4_circular_v5_hybrid.png\n")

cat("\n============================================================\n")
cat("COMPLETE: Created 4 optimized versions\n")
cat("============================================================\n")
cat("  1. fig4_circular_v2.png       - Balanced arcs, radial labels\n")
cat("  2. fig4_circular_v3_repel.png - ggrepel labels, corner guild names\n")
cat("  3. fig4_circular_v4_panels.png - 2x2 guild panels (max clarity)\n")
cat("  4. fig4_circular_v5_hybrid.png - Full network + legend inset\n")
cat("\nKey improvements:\n")
cat("  - Arc sizes use sqrt(n) for visual balance\n")
cat("  - Only top 30% between-guild edges shown (reduces clutter)\n")
cat("  - Guild labels include species counts\n")
cat("  - Better label positioning to prevent overlaps\n")
cat("  - Improved title with network statistics\n")
cat("============================================================\n")
