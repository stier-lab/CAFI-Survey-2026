# ============================================================================
# Figure 4: 5-Panel Network Figure (Version 2 - Optimized)
# ============================================================================
# Panel A: ALL 58 species in circular layout grouped by guild (hero panel)
# Panels B-E: Individual guild networks with species labels (force layout)
# ============================================================================
#
# Design principles applied:
#   - Clean, publication-ready aesthetic
#   - Node sizes reflect degree centrality
#   - Edge opacity/weight reflects correlation strength
#   - Guild colors are colorblind-safe (Okabe-Ito derivatives)
#   - Subtle guild arc backgrounds behind species groups
#   - Gaps between guilds for visual separation
#   - Visual hierarchy: within-guild edges more prominent
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
cat("CREATING 5-PANEL NETWORK FIGURE (V2 - OPTIMIZED)\n")
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

guild_counts <- sp_info %>% count(guild) %>% arrange(guild)

# Colorblind-safe guild colors (Okabe-Ito derivatives)
guild_colors <- c(
  "1" = "#0072B2",  # Blue
  "2" = "#D55E00",  # Vermillion
  "3" = "#009E73",  # Teal
  "4" = "#CC79A7"   # Pink
)

# Lighter versions for arc backgrounds
guild_colors_light <- c(
  "1" = "#C6DDED",
  "2" = "#F5D5C3",
  "3" = "#C6E8E0",
  "4" = "#F0DAE8"
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
# PANEL A: HERO CIRCULAR NETWORK - ALL SPECIES
# ============================================================================

cat("Building Panel A (circular species network)...\n")

# Calculate positions for each species on the circle
# Group by guild with gaps between guilds

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

# Generate bezier curves for all edges
cat("  Generating bezier curves for", nrow(edge_df), "edges...\n")

bezier_data <- map_dfr(1:nrow(edge_df), function(i) {
  row <- edge_df[i,]
  pts <- create_bezier(row$x1, row$y1, row$x2, row$y2)
  pts$edge_id <- i
  pts$weight <- row$weight
  pts$is_within_guild <- row$is_within_guild
  pts$guild <- as.character(row$guild_from)
  pts$edge_alpha <- row$edge_alpha
  pts
})

# Create guild arc backgrounds
guild_arc_data <- data.frame(
  guild = factor(1:4),
  start = pi/2 - guild_starts - guild_arcs,
  end = pi/2 - guild_starts,
  r_inner = 3.8,
  r_outer = 5.2
) %>%
  mutate(
    fill = guild_colors_light[as.character(guild)],
    color = guild_colors[as.character(guild)]
  )

# Guild label positions (outside the arc)
guild_label_positions <- data.frame(
  guild = factor(1:4),
  angle = pi/2 - (guild_starts + guild_arcs/2),
  label = c("Protective\nCrustaceans", "Core\nAssociates",
            "Fish &\nTrapezia", "Peripheral\nCrabs")
) %>%
  mutate(
    radius = 6.2,
    x = radius * cos(angle),
    y = radius * sin(angle),
    color = guild_colors[as.character(guild)]
  )

# Build Panel A
p_A <- ggplot() +
  # Guild arc backgrounds (subtle)
  geom_arc_bar(
    data = guild_arc_data,
    aes(x0 = 0, y0 = 0, r0 = r_inner, r = r_outer,
        start = start, end = end, fill = guild),
    alpha = 0.3,
    color = NA
  ) +
  # Between-guild edges (subtle gray)
  geom_path(
    data = bezier_data %>% filter(!is_within_guild),
    aes(x = x, y = y, group = edge_id, alpha = weight),
    color = "gray60",
    linewidth = 0.3
  ) +
  # Within-guild edges (colored, more prominent)
  geom_path(
    data = bezier_data %>% filter(is_within_guild),
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
  # Guild labels (outside)
  geom_text(
    data = guild_label_positions,
    aes(x = x, y = y, label = label, color = guild),
    size = 3.5,
    fontface = "bold",
    lineheight = 0.85
  ) +
  # Guild species counts
  geom_text(
    data = guild_label_positions %>%
      left_join(guild_counts, by = "guild") %>%
      mutate(
        label_n = paste0("n=", n),
        x_count = (radius - 0.7) * cos(angle),
        y_count = (radius - 0.7) * sin(angle)
      ),
    aes(x = x_count, y = y_count, label = label_n),
    size = 2.8,
    color = "gray40"
  ) +
  # Scales
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +
  scale_alpha_continuous(range = c(0.08, 0.6), guide = "none") +
  coord_fixed(xlim = c(-7.5, 7.5), ylim = c(-7.5, 7.5)) +
  labs(
    title = "A. Species Co-occurrence Network",
    subtitle = sprintf("%d species | %d co-occurrences | 4 ecological guilds",
                       nrow(sp_positions), ecount(g))
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40",
                                 margin = margin(b = 5)),
    plot.margin = margin(10, 10, 10, 10)
  )

cat("  Panel A complete.\n")

# ============================================================================
# PANELS B-E: INDIVIDUAL GUILD NETWORKS WITH SPECIES LABELS
# ============================================================================

create_guild_panel <- function(guild_id, letter) {

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
  layout_fr <- layout_with_fr(g_sub, niter = 800, weights = E(g_sub)$weight)

  # Scale layout to fill space
  if (nrow(layout_fr) > 1) {
    layout_fr[,1] <- scales::rescale(layout_fr[,1], to = c(-1, 1))
    layout_fr[,2] <- scales::rescale(layout_fr[,2], to = c(-1, 1))
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

  # Background circle (subtle)
  p <- p + geom_circle(
    aes(x0 = 0, y0 = 0, r = 1.3),
    fill = guild_color_light,
    color = NA,
    alpha = 0.3
  )

  # Add edges if they exist
  if (!is.null(edge_data) && nrow(edge_data) > 0) {
    p <- p + geom_segment(
      data = edge_data,
      aes(x = x1, y = y1, xend = x2, yend = y2, alpha = weight),
      color = guild_color,
      linewidth = 0.6
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

  # Add labels with repulsion
  p <- p + geom_text_repel(
    data = node_data,
    aes(x = x, y = y, label = species_label),
    size = 2.6,
    fontface = "italic",
    color = "gray15",
    segment.color = "gray50",
    segment.size = 0.25,
    segment.alpha = 0.6,
    box.padding = 0.35,
    point.padding = 0.3,
    max.overlaps = 40,
    force = 3,
    seed = 42
  )

  p <- p +
    scale_size_identity() +
    scale_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
    coord_fixed(xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6)) +
    labs(
      title = sprintf("%s. %s", letter, guild_name),
      subtitle = sprintf("%d species | %d edges", n_species, n_edges)
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5,
                                color = guild_color),
      plot.subtitle = element_text(size = 8, hjust = 0.5, color = "gray50"),
      plot.margin = margin(5, 5, 5, 5),
      plot.background = element_rect(fill = "white", color = "gray80",
                                     linewidth = 0.5)
    )

  return(p)
}

cat("Building Panels B-E (individual guilds)...\n")

p_B <- create_guild_panel(1, "B")
cat("  Panel B complete.\n")
p_C <- create_guild_panel(2, "C")
cat("  Panel C complete.\n")
p_D <- create_guild_panel(3, "D")
cat("  Panel D complete.\n")
p_E <- create_guild_panel(4, "E")
cat("  Panel E complete.\n")

# ============================================================================
# COMBINE INTO 5-PANEL FIGURE
# ============================================================================

cat("\nCombining panels...\n")

# Layout: A on top (hero panel), B-E in 2x2 grid below
p_bottom <- (p_B | p_C) / (p_D | p_E)

p_final <- p_A / p_bottom +
  plot_layout(heights = c(1.3, 2)) +
  plot_annotation(
    title = "Figure 4: CAFI Co-occurrence Network Structure",
    subtitle = "Louvain community detection reveals four ecological guilds with distinct co-occurrence patterns",
    caption = paste0(
      "Panel A: Full network with all species arranged by guild membership. ",
      "Node size = degree centrality; edges curve through center.\n",
      "Within-guild edges (colored) emphasized over between-guild edges (gray). ",
      "Panels B-E: Species-level detail for each guild.\n",
      "Edge threshold: Spearman r > 0.3, FDR-corrected p < 0.05. ",
      "Volume-residualized presence to control for coral size."
    ),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 15)),
      plot.caption = element_text(size = 8, hjust = 0, color = "gray50",
                                  margin = margin(t = 15), lineheight = 1.3),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# Save high-resolution versions
ggsave(
  file.path(fig_dir, "fig4_5panel_v2.png"),
  p_final,
  width = 13,
  height = 16,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(manuscript_dir, "fig4_network_5panel_v2.png"),
  p_final,
  width = 13,
  height = 16,
  dpi = 300,
  bg = "white"
)

cat("\nSaved: fig4_5panel_v2.png\n")
cat("Saved: manuscript/fig4_network_5panel_v2.png\n")

# ============================================================================
# ALTERNATE: Wide landscape version
# ============================================================================

cat("\nCreating wide landscape version...\n")

# Make Panel A slightly smaller for side-by-side layout
p_A_wide <- p_A +
  coord_fixed(xlim = c(-7, 7), ylim = c(-7, 7)) +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
    plot.margin = margin(5, 5, 5, 5)
  )

# Arrange: A on left, B-E in 2x2 on right
p_right <- (p_B | p_C) / (p_D | p_E)

p_wide <- (p_A_wide | p_right) +
  plot_layout(widths = c(1.2, 2)) +
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

ggsave(
  file.path(fig_dir, "fig4_5panel_v2_wide.png"),
  p_wide,
  width = 18,
  height = 11,
  dpi = 300,
  bg = "white"
)

cat("Saved: fig4_5panel_v2_wide.png\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("FIGURE 4 (V2) COMPLETE\n")
cat("============================================================\n")
cat("Panel A (Hero):\n")
cat(sprintf("  - %d species arranged in circular layout\n", nrow(sp_positions)))
cat(sprintf("  - %d total edges (%d within-guild, %d between-guild)\n",
            nrow(edge_df),
            sum(edge_df$is_within_guild),
            sum(!edge_df$is_within_guild)))
cat("  - Species grouped by guild with visual gaps\n")
cat("  - Node size reflects degree centrality\n")
cat("  - Bezier curves emphasize within-guild connectivity\n\n")

cat("Panels B-E (Guild Details):\n")
for (i in 1:4) {
  g_sub <- induced_subgraph(g, V(g)[V(g)$name %in% (sp_info %>% filter(guild == i) %>% pull(species))])
  cat(sprintf("  - Guild %d (%s): %d species, %d edges\n",
              i, gsub("\n", " ", guild_names[as.character(i)]),
              vcount(g_sub), ecount(g_sub)))
}

cat("\nOutput files:\n")
cat("  - output/figures/06_network/fig4_5panel_v2.png (13x16, portrait)\n")
cat("  - output/figures/06_network/fig4_5panel_v2_wide.png (18x11, landscape)\n")
cat("  - output/figures/manuscript/fig4_network_5panel_v2.png\n")
cat("============================================================\n")
