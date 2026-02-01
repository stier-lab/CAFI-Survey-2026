# ============================================================================
# fig4_network_redesign.R - PUBLICATION FIGURE (Design Agency Redesign)
# ============================================================================
#
# Based on Creative Director + Graphic Designer agent feedback:
# - Orbital layout: small guilds INNER ring, large guilds OUTER ring
# - "Deep Ocean" color palette (colorblind-safe, professional)
# - Atmospheric halos instead of box labels
# - Dashed between-guild edges at low opacity
# - Proper italicized binomial species names
# ============================================================================

cat("\n============================================================\n")
cat("FIGURE 4: REDESIGNED NETWORK (Orbital Layout)\n")
cat("============================================================\n\n")

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

network_results <- load_object("cafi_network")
centrality_df <- network_results$centrality
edge_list <- network_results$edge_list

library(igraph)
library(ggrepel)
library(ggforce)  # for geom_mark_hull

fig_dir <- file.path(PATHS$figures, "06_network")
manuscript_dir <- file.path(PATHS$figures, "manuscript")

# ============================================================================
# FILTER AT r > 0.7 AND RECALCULATE DEGREE
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
  mutate(degree = filtered_degree[species])

guild_lookup <- setNames(species_info$module, species_info$species)

edges_filtered <- edges_filtered %>%
  mutate(
    guild1 = guild_lookup[sp1],
    guild2 = guild_lookup[sp2],
    same_guild = guild1 == guild2
  )

n_within <- sum(edges_filtered$same_guild)
n_between <- sum(!edges_filtered$same_guild)

cat(sprintf("Species: %d | Edges: %d (within: %d, between: %d)\n",
            length(species_in), nrow(edges_filtered), n_within, n_between))

# Guild sizes
guild_sizes <- species_info %>% count(module)
cat("\nGuild sizes:\n")
print(guild_sizes)

# ============================================================================
# NEW COLOR PALETTE: "Deep Ocean" (Colorblind-safe, professional)
# ============================================================================

guild_colors <- c(
  "1" = "#2D6A4F",  # Forest Teal (hermit crabs)
  "2" = "#1D3557",  # Navy Depth (shrimp/worms)
  "3" = "#E07A5F",  # Terracotta (fish/Trapezia) - FOCAL
  "4" = "#81B29A"   # Sage (peripheral crabs)
)

guild_colors_light <- c(
  "1" = "#2D6A4F20",  # 12% opacity for halos
  "2" = "#1D355720",
  "3" = "#E07A5F20",
  "4" = "#81B29A20"
)

# ============================================================================
# ORBITAL LAYOUT: Small guilds INNER, Large guilds OUTER
# ============================================================================

cat("\nApplying orbital layout...\n")

# Design principle: Small guilds (3,4) get INNER ring = prime visual real estate
# Large guilds (1,2) spread on OUTER ring

orbital_layout <- function(species_info) {
  n_nodes <- nrow(species_info)
  layout <- matrix(0, nrow = n_nodes, ncol = 2)

  # Ring radii
  inner_r <- 2.0   # Small guilds (closer = more prominent)
  outer_r <- 4.5   # Large guilds (spread out)

  # Arc assignments (degrees): each guild gets ~90° of arc
  # Small guilds: right and left sides of inner ring
  # Large guilds: top and bottom of outer ring
  guild_config <- list(
    "3" = list(ring = inner_r, arc_start = -45, arc_end = 45),    # Inner right (focal)
    "4" = list(ring = inner_r, arc_start = 135, arc_end = 225),   # Inner left
    "1" = list(ring = outer_r, arc_start = 30, arc_end = 150),    # Outer top-right
    "2" = list(ring = outer_r, arc_start = 210, arc_end = 330)    # Outer bottom-left
  )

  for (g in names(guild_config)) {
    cfg <- guild_config[[g]]
    guild_sp <- species_info %>%
      filter(module == as.integer(g)) %>%
      arrange(desc(degree))

    n <- nrow(guild_sp)
    if (n == 0) next

    # Distribute evenly across arc
    angles <- seq(cfg$arc_start, cfg$arc_end, length.out = n)
    angles_rad <- angles * pi / 180

    idx <- match(guild_sp$species, species_info$species)
    layout[idx, 1] <- cfg$ring * cos(angles_rad)
    layout[idx, 2] <- cfg$ring * sin(angles_rad)
  }

  return(layout)
}

layout_matrix <- orbital_layout(species_info)

# ============================================================================
# FORMAT SPECIES NAMES: Proper italicized binomials
# ============================================================================

format_species_name <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  if (length(parts) >= 2) {
    # Genus initial + full species epithet
    genus <- substr(parts[1], 1, 1)
    epithet <- parts[2]
    return(paste0(genus, ". ", epithet))
  } else {
    # Family-level (no underscore) - return as-is, will show in roman
    return(name)
  }
}

# ============================================================================
# BUILD NODE DATAFRAME
# ============================================================================

node_df <- data.frame(
  x = layout_matrix[, 1],
  y = layout_matrix[, 2],
  species = species_info$species,
  module = as.character(species_info$module),
  type = species_info$type,
  degree = species_info$degree
) %>%
  mutate(
    label = sapply(species, format_species_name),
    # Node size scaled by degree (larger range for visual impact)
    node_size = 4 + (degree / max(degree)) * 12
  )

# ============================================================================
# BUILD EDGE DATAFRAMES
# ============================================================================

within_edges <- edges_filtered %>%
  filter(same_guild) %>%
  mutate(
    x = layout_matrix[match(sp1, species_info$species), 1],
    y = layout_matrix[match(sp1, species_info$species), 2],
    xend = layout_matrix[match(sp2, species_info$species), 1],
    yend = layout_matrix[match(sp2, species_info$species), 2],
    edge_width = 0.3 + (correlation - 0.7) * 3  # Scale 0.7-1.0 to 0.3-1.2
  )

between_edges <- edges_filtered %>%
  filter(!same_guild) %>%
  mutate(
    x = layout_matrix[match(sp1, species_info$species), 1],
    y = layout_matrix[match(sp2, species_info$species), 1],
    y = layout_matrix[match(sp1, species_info$species), 2],
    xend = layout_matrix[match(sp2, species_info$species), 1],
    yend = layout_matrix[match(sp2, species_info$species), 2]
  )

# ============================================================================
# GUILD SUMMARIES FOR LABELS
# ============================================================================

guild_descriptions <- c(
  "1" = "Hermit crabs\n& shrimp",
  "2" = "Shrimp, worms\n& echinoderms",
  "3" = "Fish & guard\ncrabs",
  "4" = "Peripheral\ncrabs"
)

guild_summary <- species_info %>%
  group_by(module) %>%
  summarise(
    n = n(),
    cx = mean(layout_matrix[match(species, species_info$species), 1]),
    cy = mean(layout_matrix[match(species, species_info$species), 2]),
    .groups = "drop"
  ) %>%
  mutate(
    description = guild_descriptions[as.character(module)],
    # Label position: offset from centroid
    label_x = cx * 1.35,
    label_y = cy * 1.35
  )

# ============================================================================
# SELECT HUB SPECIES FOR LABELING (top 3 per guild)
# ============================================================================

hub_labels <- node_df %>%
  group_by(module) %>%
  slice_max(degree, n = 3, with_ties = FALSE) %>%
  ungroup()

cat(sprintf("\nLabeling %d hub species\n", nrow(hub_labels)))

# ============================================================================
# CREATE PUBLICATION FIGURE
# ============================================================================

cat("Creating orbital network figure...\n")

p_orbital <- ggplot() +

  # LAYER 1: Atmospheric halos (guild backgrounds)
  geom_mark_hull(
    data = node_df,
    aes(x = x, y = y, fill = module),
    expand = unit(12, "mm"),
    radius = unit(8, "mm"),
    alpha = 0.08,
    color = NA,
    show.legend = FALSE
  ) +

  # LAYER 2: Between-guild edges (dashed, very faint)
  geom_curve(
    data = between_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "#CCCCCC",
    linewidth = 0.3,
    alpha = 0.25,
    linetype = "dashed",
    curvature = 0.2
  ) +

  # LAYER 3: Within-guild edges (solid, colored by guild)
  geom_segment(
    data = within_edges %>% left_join(node_df %>% dplyr::select(species, module),
                                       by = c("sp1" = "species")),
    aes(x = x, y = y, xend = xend, yend = yend,
        color = module, linewidth = edge_width, alpha = correlation),
    show.legend = FALSE
  ) +

  # LAYER 4: Nodes
  geom_point(
    data = node_df,
    aes(x = x, y = y, fill = module, size = node_size),
    shape = 21,
    color = "white",
    stroke = 1.2
  ) +

  # LAYER 5: Hub species labels
  geom_text_repel(
    data = hub_labels,
    aes(x = x, y = y, label = label, color = module),
    size = 3.5,
    fontface = "italic",
    max.overlaps = 20,
    force = 3,
    segment.color = "gray60",
    segment.alpha = 0.5,
    point.padding = 0.5,
    box.padding = 0.4,
    bg.color = "white",
    bg.r = 0.15,
    show.legend = FALSE
  ) +

  # LAYER 6: Guild labels
  geom_label(
    data = guild_summary,
    aes(x = label_x, y = label_y,
        label = paste0("Guild ", module, "\n", n, " species\n", description),
        fill = as.character(module)),
    color = "white",
    fontface = "bold",
    size = 3.2,
    label.padding = unit(0.4, "lines"),
    label.r = unit(0.3, "lines")
  ) +

  # SCALES
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_identity() +  # Use actual node_size values
  scale_linewidth_identity() +
  scale_alpha_continuous(range = c(0.4, 0.9), guide = "none") +

  # LABELS
  labs(
    title = "Coral-Associated Fauna Form Four Ecological Guilds",
    subtitle = sprintf("Co-occurrence network at r > 0.7 | %d species, %d associations",
                       nrow(node_df), nrow(edges_filtered)),
    caption = "Guilds identified via Louvain community detection on volume-corrected co-occurrence.\nNode size = network degree. Solid lines = within-guild; dashed = between-guild associations."
  ) +

  # THEME
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      color = "#666666",
      margin = margin(b = 15)
    ),
    plot.caption = element_text(
      size = 9,
      hjust = 0,
      color = "#888888",
      lineheight = 1.2,
      margin = margin(t = 15)
    ),
    plot.background = element_rect(fill = "#FAFAFA", color = NA),
    plot.margin = margin(25, 25, 25, 25)
  ) +
  coord_fixed(xlim = c(-7, 7), ylim = c(-6.5, 7))

# ============================================================================
# SAVE OUTPUTS
# ============================================================================

ggsave(file.path(fig_dir, "fig4_network_orbital.png"), p_orbital,
       width = 12, height = 11, dpi = 300, bg = "#FAFAFA")

ggsave(file.path(manuscript_dir, "fig4_network.png"), p_orbital,
       width = 12, height = 11, dpi = 300, bg = "#FAFAFA")

cat("\n============================================================\n")
cat("ORBITAL NETWORK FIGURE COMPLETE\n")
cat("============================================================\n\n")
cat("Files saved:\n")
cat("  - output/figures/06_network/fig4_network_orbital.png\n")
cat("  - output/figures/manuscript/fig4_network.png\n")
cat("\nDesign features:\n")
cat("  - Orbital layout: small guilds INNER (prominent), large guilds OUTER\n")
cat("  - Deep Ocean palette (colorblind-safe)\n")
cat("  - Atmospheric halos for guild backgrounds\n")
cat("  - Dashed between-guild edges at low opacity\n")
cat("  - Italicized species names for hub species\n")
