# ============================================================================
# Figure 4: FINAL Publication-Ready Circular Guild Network
# Maximum clarity: minimal edges, clear labels, professional styling
# ============================================================================

library(tidyverse)
library(igraph)
library(ggforce)
library(ggrepel)
library(patchwork)

# Setup paths
source("scripts/00_setup.R")

# Output directory
fig_dir <- file.path(PATHS$figures, "06_network")
manuscript_dir <- file.path(PATHS$figures, "manuscript")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("CREATING FINAL PUBLICATION FIGURE\n")
cat("============================================================\n\n")

# Load network results
network_results <- load_object("cafi_network")
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

# Guild naming
guild_names <- c(
  "1" = "Protective Crustaceans",
  "2" = "Core Associates",
  "3" = "Fish & Trapezia",
  "4" = "Peripheral Crabs"
)

guild_counts <- sp_info %>% count(guild) %>% arrange(guild)
guild_counts$name <- guild_names[as.character(guild_counts$guild)]

cat("Guild sizes:\n")
for (i in 1:nrow(guild_counts)) {
  cat(sprintf("  %s: %d species\n", guild_counts$name[i], guild_counts$n[i]))
}

# Guild colors - colorblind safe, high contrast
guild_colors <- c(
  "1" = "#0072B2",  # Blue
  "2" = "#D55E00",  # Vermillion
  "3" = "#009E73",  # Teal
  "4" = "#CC79A7"   # Pink
)

# ============================================================================
# BALANCED CIRCULAR LAYOUT
# ============================================================================

create_balanced_layout <- function(sp_info, radius = 4) {
  sp_ordered <- sp_info %>% arrange(guild, desc(degree))

  guild_counts <- sp_ordered %>% count(guild) %>% arrange(guild)
  n_guilds <- nrow(guild_counts)

  gap_angle <- 15 * pi / 180  # 15 degree gaps
  available_angle <- 2 * pi - n_guilds * gap_angle

  # Square root scaling for visual balance
  guild_counts <- guild_counts %>%
    mutate(
      weight = sqrt(n),
      arc_size = (weight / sum(weight)) * available_angle
    )

  # Position starting at top, clockwise
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

  # Position species
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
      padding <- (arc$start - arc$end) * 0.06
      angles <- seq(arc$start - padding, arc$end + padding, length.out = n_sp)
    }

    sp_ordered$x[guild_sp] <- radius * cos(angles)
    sp_ordered$y[guild_sp] <- radius * sin(angles)
    sp_ordered$angle[guild_sp] <- angles
  }

  list(layout = sp_ordered, arcs = guild_arcs, counts = guild_counts)
}

layout_result <- create_balanced_layout(sp_info, radius = 4)
node_df <- layout_result$layout
guild_arcs <- layout_result$arcs

# ============================================================================
# EDGE FILTERING - ONLY WITHIN-GUILD + TOP BETWEEN
# ============================================================================

edge_list <- as_data_frame(g, what = "edges")

edge_list <- edge_list %>%
  left_join(node_df %>% dplyr::select(species, guild, x, y),
            by = c("from" = "species")) %>%
  rename(x1 = x, y1 = y, guild_from = guild) %>%
  left_join(node_df %>% dplyr::select(species, guild, x, y),
            by = c("to" = "species")) %>%
  rename(x2 = x, y2 = y, guild_to = guild) %>%
  mutate(is_within = guild_from == guild_to)

# Keep all within-guild + top 15% between-guild
within_edges <- edge_list %>% filter(is_within)
between_edges <- edge_list %>%
  filter(!is_within) %>%
  arrange(desc(weight)) %>%
  slice_head(prop = 0.15)

filtered_edges <- bind_rows(within_edges, between_edges)

cat(sprintf("\nEdges: %d within-guild + %d between-guild = %d total (from %d)\n",
            nrow(within_edges), nrow(between_edges), nrow(filtered_edges), nrow(edge_list)))

# ============================================================================
# CREATE BEZIER CURVES
# ============================================================================

create_bezier <- function(x1, y1, x2, y2, n = 30) {
  dist <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  pull <- 0.2 + 0.15 * (dist / 8)
  cx <- (x1 + x2) / 2 * (1 - pull)
  cy <- (y1 + y2) / 2 * (1 - pull)
  t <- seq(0, 1, length.out = n)
  data.frame(
    x = (1-t)^2 * x1 + 2*(1-t)*t * cx + t^2 * x2,
    y = (1-t)^2 * y1 + 2*(1-t)*t * cy + t^2 * y2
  )
}

bezier_edges <- map2_dfr(
  seq_len(nrow(filtered_edges)),
  filtered_edges$weight,
  function(i, w) {
    row <- filtered_edges[i, ]
    pts <- create_bezier(row$x1, row$y1, row$x2, row$y2)
    pts$edge_id <- i
    pts$weight <- w
    pts$is_within <- row$is_within
    pts$guild <- row$guild_from
    pts
  }
)

# ============================================================================
# NODE SIZING AND LABELING
# ============================================================================

node_df <- node_df %>%
  mutate(
    node_size = scales::rescale(degree, to = c(3, 12)),
    # Label positioning
    label_angle_deg = angle * 180 / pi,
    # Flip left side labels
    text_angle = case_when(
      label_angle_deg > 90 | label_angle_deg < -90 ~ label_angle_deg + 180,
      TRUE ~ label_angle_deg
    ),
    hjust_val = case_when(
      label_angle_deg > 90 | label_angle_deg < -90 ~ 1,
      TRUE ~ 0
    )
  )

# Select species to label: all from small guilds, top 6 from large guilds
labels_df <- node_df %>%
  group_by(guild) %>%
  arrange(desc(degree)) %>%
  mutate(
    rank = row_number(),
    n_guild = n(),
    show = case_when(
      n_guild <= 8 ~ TRUE,
      rank <= 6 ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(show) %>%
  ungroup()

cat(sprintf("Labeling %d of %d species\n", nrow(labels_df), nrow(node_df)))

# ============================================================================
# FINAL FIGURE
# ============================================================================

cat("\nBuilding final figure...\n")

p_final <- ggplot() +

  # 1. Arc backgrounds
  geom_arc(
    data = data.frame(
      guild = names(guild_arcs),
      start = sapply(guild_arcs, function(x) x$end),
      end = sapply(guild_arcs, function(x) x$start)
    ),
    aes(x0 = 0, y0 = 0, r = 4, start = start, end = end, color = guild),
    linewidth = 14,
    alpha = 0.18
  ) +

  # 2. Between-guild edges (very subtle gray)
  geom_path(
    data = bezier_edges %>% filter(!is_within),
    aes(x = x, y = y, group = edge_id),
    color = "gray75",
    linewidth = 0.25,
    alpha = 0.35
  ) +

  # 3. Within-guild edges (colored, stronger)
  geom_path(
    data = bezier_edges %>% filter(is_within),
    aes(x = x, y = y, group = edge_id, color = guild),
    linewidth = 0.45,
    alpha = 0.55
  ) +

  # 4. Nodes
  geom_point(
    data = node_df,
    aes(x = x, y = y, fill = guild, size = node_size),
    shape = 21,
    color = "white",
    stroke = 0.6
  ) +

  # 5. Species labels with repel
  geom_text_repel(
    data = labels_df,
    aes(x = x, y = y, label = species, color = guild),
    size = 2.6,
    fontface = "italic",
    segment.color = "gray65",
    segment.size = 0.25,
    segment.alpha = 0.5,
    min.segment.length = 0.2,
    box.padding = 0.45,
    point.padding = 0.35,
    max.overlaps = 25,
    force = 3,
    force_pull = 0.3
  ) +

  # Scales
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +

  coord_fixed(xlim = c(-6.5, 6.5), ylim = c(-6.5, 6.5)) +

  theme_void() +
  theme(
    plot.margin = margin(5, 5, 5, 5)
  )

# ============================================================================
# ADD LEGEND PANEL
# ============================================================================

legend_df <- data.frame(
  guild = as.character(1:4),
  name = guild_names[as.character(1:4)],
  n = guild_counts$n,
  y = 4:1
)

p_legend <- ggplot(legend_df) +
  geom_point(aes(x = 0, y = y, fill = guild),
             size = 8, shape = 21, color = "white", stroke = 0.8) +
  geom_text(aes(x = 0.5, y = y, label = sprintf("%s (n=%d)", name, n)),
            hjust = 0, size = 3.5, fontface = "bold", color = "gray20") +
  scale_fill_manual(values = guild_colors, guide = "none") +
  xlim(-0.5, 6) +
  ylim(0, 5) +
  theme_void()

# ============================================================================
# COMBINE WITH PATCHWORK
# ============================================================================

p_combined <- p_final +
  inset_element(p_legend, left = 0, bottom = 0, right = 0.5, top = 0.28) +
  plot_annotation(
    title = "CAFI Co-occurrence Network: Four Ecological Guilds",
    subtitle = sprintf(
      "%d species | %d edges (Spearman r > 0.7) | Louvain modularity Q = %.2f",
      nrow(node_df), nrow(filtered_edges), modularity(communities)
    ),
    caption = paste0(
      "Nodes sized by degree centrality (number of connections). ",
      "Colored edges = within-guild co-occurrences; gray = between-guild.\n",
      "Correlations volume-corrected to remove size confound."
    ),
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5,
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = 10.5, hjust = 0.5, color = "gray35",
                                   margin = margin(b = 8)),
      plot.caption = element_text(size = 8, hjust = 0.5, color = "gray50",
                                  margin = margin(t = 12), lineheight = 1.15),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# Save to network folder
ggsave(
  file.path(fig_dir, "fig4_network_final.png"),
  p_combined,
  width = 11,
  height = 11,
  dpi = 300,
  bg = "white"
)

# Also save to manuscript folder
ggsave(
  file.path(manuscript_dir, "fig4_network.png"),
  p_combined,
  width = 11,
  height = 11,
  dpi = 300,
  bg = "white"
)

cat("\nSaved: fig4_network_final.png\n")
cat("Saved: manuscript/fig4_network.png\n")

# ============================================================================
# ALTERNATE: DARK THEME VERSION
# ============================================================================

cat("\nCreating dark theme version...\n")

p_dark <- ggplot() +

  # Arc backgrounds (brighter on dark)
  geom_arc(
    data = data.frame(
      guild = names(guild_arcs),
      start = sapply(guild_arcs, function(x) x$end),
      end = sapply(guild_arcs, function(x) x$start)
    ),
    aes(x0 = 0, y0 = 0, r = 4, start = start, end = end, color = guild),
    linewidth = 14,
    alpha = 0.25
  ) +

  # Between-guild edges
  geom_path(
    data = bezier_edges %>% filter(!is_within),
    aes(x = x, y = y, group = edge_id),
    color = "gray45",
    linewidth = 0.2,
    alpha = 0.3
  ) +

  # Within-guild edges
  geom_path(
    data = bezier_edges %>% filter(is_within),
    aes(x = x, y = y, group = edge_id, color = guild),
    linewidth = 0.5,
    alpha = 0.6
  ) +

  # Nodes (white stroke visible on dark)
  geom_point(
    data = node_df,
    aes(x = x, y = y, fill = guild, size = node_size),
    shape = 21,
    color = "gray90",
    stroke = 0.5
  ) +

  # Labels
  geom_text_repel(
    data = labels_df,
    aes(x = x, y = y, label = species, color = guild),
    size = 2.6,
    fontface = "italic",
    segment.color = "gray55",
    segment.size = 0.25,
    segment.alpha = 0.6,
    box.padding = 0.45,
    point.padding = 0.35,
    max.overlaps = 25,
    force = 3
  ) +

  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +

  coord_fixed(xlim = c(-6.5, 6.5), ylim = c(-6.5, 6.5)) +

  labs(
    title = "CAFI Co-occurrence Network",
    subtitle = sprintf("%d species across 4 ecological guilds | Q = %.2f",
                       nrow(node_df), modularity(communities))
  ) +

  theme_void() +
  theme(
    plot.background = element_rect(fill = "#1a1a2e", color = NA),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                              color = "white", margin = margin(t = 10, b = 5)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray70",
                                 margin = margin(b = 10)),
    plot.margin = margin(15, 15, 15, 15)
  )

# Add dark legend
legend_dark <- ggplot(legend_df) +
  geom_point(aes(x = 0, y = y, fill = guild),
             size = 6, shape = 21, color = "gray80", stroke = 0.6) +
  geom_text(aes(x = 0.4, y = y, label = sprintf("%s (%d)", name, n)),
            hjust = 0, size = 3, fontface = "bold", color = "gray85") +
  scale_fill_manual(values = guild_colors, guide = "none") +
  xlim(-0.3, 5) +
  ylim(0, 5) +
  theme_void() +
  theme(plot.background = element_rect(fill = "#1a1a2e", color = NA))

p_dark_final <- p_dark +
  inset_element(legend_dark, left = 0, bottom = 0, right = 0.55, top = 0.28)

ggsave(
  file.path(fig_dir, "fig4_network_dark.png"),
  p_dark_final,
  width = 10,
  height = 10,
  dpi = 300
)

cat("Saved: fig4_network_dark.png\n")

# ============================================================================
# MINIMAL VERSION (for text-heavy publications)
# ============================================================================

cat("\nCreating minimal version...\n")

# Even fewer labels - only top 3 per guild
minimal_labels <- node_df %>%
  group_by(guild) %>%
  slice_max(degree, n = 3) %>%
  ungroup()

p_minimal <- ggplot() +

  # Arcs
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

  # Only within-guild edges
  geom_path(
    data = bezier_edges %>% filter(is_within),
    aes(x = x, y = y, group = edge_id, color = guild),
    linewidth = 0.35,
    alpha = 0.4
  ) +

  # Nodes
  geom_point(
    data = node_df,
    aes(x = x, y = y, fill = guild, size = node_size * 0.8),
    shape = 21,
    color = "white",
    stroke = 0.5
  ) +

  # Minimal labels
  geom_text_repel(
    data = minimal_labels,
    aes(x = x, y = y, label = species, color = guild),
    size = 2.8,
    fontface = "italic",
    segment.color = "gray60",
    segment.size = 0.2,
    box.padding = 0.5,
    max.overlaps = 15
  ) +

  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +

  coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-5.5, 5.5)) +

  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

# Compact legend
legend_compact <- data.frame(
  guild = as.character(1:4),
  label = paste0(guild_names[as.character(1:4)], " (", guild_counts$n, ")"),
  x = c(1, 2, 3, 4)
)

p_minimal_legend <- ggplot(legend_compact) +
  geom_point(aes(x = x, y = 1, fill = guild), size = 4, shape = 21,
             color = "white", stroke = 0.5) +
  geom_text(aes(x = x, y = 0.7, label = label),
            size = 2.5, fontface = "bold", hjust = 0.5) +
  scale_fill_manual(values = guild_colors, guide = "none") +
  xlim(0.5, 4.5) + ylim(0.4, 1.2) +
  theme_void()

p_minimal_final <- p_minimal / p_minimal_legend +
  plot_layout(heights = c(10, 1)) +
  plot_annotation(
    title = "CAFI Co-occurrence Network",
    subtitle = sprintf("%d species | 4 guilds | Modularity = %.2f",
                       nrow(node_df), modularity(communities)),
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40")
    )
  )

ggsave(
  file.path(fig_dir, "fig4_network_minimal.png"),
  p_minimal_final,
  width = 8,
  height = 9,
  dpi = 300,
  bg = "white"
)

cat("Saved: fig4_network_minimal.png\n")

cat("\n============================================================\n")
cat("COMPLETE: Final publication figures created\n")
cat("============================================================\n")
cat("  Main figure: fig4_network_final.png (11x11)\n")
cat("  Dark theme:  fig4_network_dark.png (10x10)\n")
cat("  Minimal:     fig4_network_minimal.png (8x9)\n")
cat("  Manuscript:  manuscript/fig4_network.png\n")
cat("============================================================\n")
