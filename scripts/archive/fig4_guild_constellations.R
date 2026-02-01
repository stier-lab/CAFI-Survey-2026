# ============================================================================
# fig4_guild_constellations.R - ART DIRECTOR'S TOP RECOMMENDATION
# ============================================================================
#
# "Guild Constellations" - Packed Circles with Internal Networks
# Each guild is a large circle containing its species network
# Circles sized for VISUAL BALANCE, not data proportion
# Between-guild edges cross the void between circles
#
# This design was recommended because:
# 1. Readers immediately grasp "4 groups" from 4 circles
# 2. Works at both full-page and thumbnail sizes
# 3. Can show both within-guild and between-guild structure
# 4. Unique and memorable for journal cover potential
# ============================================================================

cat("\n============================================================\n")
cat("FIGURE 4: GUILD CONSTELLATIONS (Packed Circles)\n")
cat("============================================================\n\n")

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

network_results <- load_object("cafi_network")
centrality_df <- network_results$centrality
edge_list <- network_results$edge_list

library(igraph)
library(ggrepel)
library(ggforce)

fig_dir <- file.path(PATHS$figures, "06_network")
manuscript_dir <- file.path(PATHS$figures, "manuscript")

# ============================================================================
# DATA SETUP
# ============================================================================

threshold <- 0.7
edges_filtered <- edge_list %>% filter(correlation > threshold)
species_in <- unique(c(edges_filtered$sp1, edges_filtered$sp2))

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

cat(sprintf("Species: %d | Edges: %d\n", length(species_in), nrow(edges_filtered)))

# ============================================================================
# COLOR PALETTE - Refined for elegance
# ============================================================================

guild_colors <- c(
  "1" = "#2D6A4F",  # Forest Teal
  "2" = "#1D3557",  # Navy Depth
  "3" = "#E07A5F",  # Terracotta (FOCAL - fish & Trapezia)
  "4" = "#81B29A"   # Sage
)

guild_colors_light <- c(
  "1" = "#2D6A4F15",
  "2" = "#1D355715",
  "3" = "#E07A5F18",  # Slightly more visible for focal guild
  "4" = "#81B29A15"
)

# ============================================================================
# CIRCLE DEFINITIONS - Visual balance, NOT data-proportional
# ============================================================================

# Key insight: Give small guilds LARGER circles relative to species count
# This makes each species more prominent

circles <- data.frame(
  guild = c(1, 2, 3, 4),
  x0 = c(2.8, -2.8, 2.5, -2.2),     # Circle centers
  y0 = c(2.2, 2.5, -2.5, -2.2),
  r = c(2.8, 3.0, 2.2, 2.0),        # Radii - NOT proportional to species count!
  n_species = c(13, 17, 4, 4),
  label = c("Guild 1\nHermit Crabs\n& Shrimp",
            "Guild 2\nShrimp, Worms\n& Echinoderms",
            "Guild 3\nFish & Guard\nCrabs",
            "Guild 4\nPeripheral\nCrabs")
)

# ============================================================================
# POSITION SPECIES WITHIN CIRCLES
# Use force-directed layout CONSTRAINED to each circle
# ============================================================================

cat("Computing constrained layouts within circles...\n")

position_species_in_circle <- function(guild_num, species_info, edges_filtered, circles) {
  guild_sp <- species_info %>% filter(module == guild_num)
  n <- nrow(guild_sp)

  if (n == 0) return(NULL)

  circle_info <- circles %>% filter(guild == guild_num)
  x0 <- circle_info$x0
  y0 <- circle_info$y0
  r <- circle_info$r

  # Get within-guild edges for this guild
  guild_edges <- edges_filtered %>%
    filter(guild1 == guild_num & guild2 == guild_num)

  if (nrow(guild_edges) > 0 && n > 2) {
    # Create subgraph and get force-directed layout
    g_sub <- graph_from_data_frame(
      guild_edges %>% dplyr::select(sp1, sp2, correlation),
      directed = FALSE,
      vertices = data.frame(name = guild_sp$species)
    )
    set.seed(42 + guild_num)
    layout_sub <- layout_with_fr(g_sub)

    # Scale to fit within circle (use 75% of radius for padding)
    if (max(abs(layout_sub)) > 0) {
      scale_factor <- (r * 0.75) / max(abs(layout_sub))
      layout_sub <- layout_sub * scale_factor
    }

    guild_sp$x <- x0 + layout_sub[, 1]
    guild_sp$y <- y0 + layout_sub[, 2]
  } else {
    # For small guilds or no edges, arrange in a small circle
    angles <- seq(0, 2*pi, length.out = n + 1)[1:n]
    guild_sp$x <- x0 + (r * 0.5) * cos(angles)
    guild_sp$y <- y0 + (r * 0.5) * sin(angles)
  }

  return(guild_sp)
}

# Position all species
all_positions <- bind_rows(
  position_species_in_circle(1, species_info, edges_filtered, circles),
  position_species_in_circle(2, species_info, edges_filtered, circles),
  position_species_in_circle(3, species_info, edges_filtered, circles),
  position_species_in_circle(4, species_info, edges_filtered, circles)
)

# ============================================================================
# FORMAT SPECIES NAMES
# ============================================================================

format_name <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  if (length(parts) >= 2) {
    return(paste0(substr(parts[1], 1, 1), ". ", parts[2]))
  }
  return(name)
}

all_positions <- all_positions %>%
  mutate(
    label = sapply(species, format_name),
    module = as.character(module)
  )

# ============================================================================
# BUILD EDGE DATAFRAMES WITH POSITIONS
# ============================================================================

pos_lookup_x <- setNames(all_positions$x, all_positions$species)
pos_lookup_y <- setNames(all_positions$y, all_positions$species)

within_edges <- edges_filtered %>%
  filter(same_guild) %>%
  mutate(
    x = pos_lookup_x[sp1],
    y = pos_lookup_y[sp1],
    xend = pos_lookup_x[sp2],
    yend = pos_lookup_y[sp2],
    guild = as.character(guild1)
  ) %>%
  filter(!is.na(x) & !is.na(xend))

between_edges <- edges_filtered %>%
  filter(!same_guild) %>%
  mutate(
    x = pos_lookup_x[sp1],
    y = pos_lookup_y[sp1],
    xend = pos_lookup_x[sp2],
    yend = pos_lookup_y[sp2]
  ) %>%
  filter(!is.na(x) & !is.na(xend))

cat(sprintf("Within-guild edges: %d | Between-guild edges: %d\n",
            nrow(within_edges), nrow(between_edges)))

# ============================================================================
# SELECT HUB SPECIES FOR LABELING
# ============================================================================

hub_species <- all_positions %>%
  group_by(module) %>%
  slice_max(degree, n = 3, with_ties = FALSE) %>%
  ungroup()

# ============================================================================
# CREATE THE FIGURE
# ============================================================================

cat("Creating Guild Constellations figure...\n")

p_constellations <- ggplot() +

  # LAYER 1: Guild circle backgrounds (soft fills)
  geom_circle(
    data = circles,
    aes(x0 = x0, y0 = y0, r = r, fill = as.character(guild)),
    alpha = 0.12,
    color = NA
  ) +

  # LAYER 2: Guild circle borders (elegant thin strokes)
  geom_circle(
    data = circles,
    aes(x0 = x0, y0 = y0, r = r, color = as.character(guild)),
    fill = NA,
    linewidth = 1.5,
    alpha = 0.6
  ) +

  # LAYER 3: Between-guild edges (crossing the void)
  geom_curve(
    data = between_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "#888888",
    linewidth = 0.4,
    alpha = 0.35,
    curvature = 0.15
  ) +

  # LAYER 4: Within-guild edges (inside circles)
  geom_segment(
    data = within_edges,
    aes(x = x, y = y, xend = xend, yend = yend, color = guild),
    linewidth = 0.5,
    alpha = 0.4
  ) +

  # LAYER 5: Species nodes
  geom_point(
    data = all_positions,
    aes(x = x, y = y, fill = module, size = degree),
    shape = 21,
    color = "white",
    stroke = 0.8
  ) +

  # LAYER 6: Hub species labels (elegant, minimal)
  geom_text_repel(
    data = hub_species,
    aes(x = x, y = y, label = label, color = module),
    size = 3.2,
    fontface = "italic",
    max.overlaps = 15,
    force = 2,
    segment.color = "gray50",
    segment.alpha = 0.5,
    segment.size = 0.3,
    point.padding = 0.4,
    box.padding = 0.3,
    bg.color = "white",
    bg.r = 0.1,
    show.legend = FALSE
  ) +

  # LAYER 7: Guild labels (inside circles, top area)
  geom_label(
    data = circles %>% mutate(
      label_y = y0 + r * 0.65  # Position in upper part of circle
    ),
    aes(x = x0, y = label_y, label = label, fill = as.character(guild)),
    color = "white",
    fontface = "bold",
    size = 3,
    label.padding = unit(0.35, "lines"),
    label.r = unit(0.2, "lines"),
    alpha = 0.95
  ) +

  # LAYER 8: Species counts (subtle)
  geom_text(
    data = circles %>% mutate(
      count_y = y0 - r * 0.75
    ),
    aes(x = x0, y = count_y, label = paste0("n = ", n_species)),
    color = "#666666",
    size = 3,
    fontface = "italic"
  ) +

  # SCALES

  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_size_continuous(range = c(3, 12), name = "Network\nDegree") +

  # LABELS
  labs(
    title = "Four Ecological Guilds Structure Coral-Associated Communities",
    subtitle = "Co-occurrence network (r > 0.7) | Species positioned by association strength within guilds",
    caption = "Each circle represents a guild identified via Louvain community detection.\nLines within circles = species co-occurrence; lines between circles = cross-guild associations.\nNode size reflects network connectivity (degree)."
  ) +

  # THEME
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5,
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0.5,
      color = "#555555",
      margin = margin(b = 15)
    ),
    plot.caption = element_text(
      size = 9,
      hjust = 0,
      color = "#777777",
      lineheight = 1.3,
      margin = margin(t = 20)
    ),
    legend.position = c(0.92, 0.15),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "#FAFAFA", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_fixed(xlim = c(-7, 7), ylim = c(-6.5, 6.5))

# ============================================================================
# SAVE
# ============================================================================

ggsave(file.path(fig_dir, "fig4_guild_constellations.png"), p_constellations,
       width = 13, height = 11, dpi = 300, bg = "#FAFAFA")

ggsave(file.path(manuscript_dir, "fig4_network.png"), p_constellations,
       width = 13, height = 11, dpi = 300, bg = "#FAFAFA")

cat("\n============================================================\n")
cat("GUILD CONSTELLATIONS COMPLETE\n")
cat("============================================================\n")
cat("\nSaved:\n")
cat("  - output/figures/06_network/fig4_guild_constellations.png\n")
cat("  - output/figures/manuscript/fig4_network.png\n")
cat("\nDesign features:\n")
cat("  - 4 circles with VISUAL BALANCE (not data-proportional)\n")
cat("  - Species positioned by force-directed layout within circles\n")
cat("  - Between-guild edges cross the void elegantly\n")
cat("  - Hub species labeled with italicized names\n")
cat("  - Suitable for journal cover potential\n")
