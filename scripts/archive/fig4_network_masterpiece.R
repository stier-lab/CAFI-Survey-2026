# ============================================================================
# fig4_network_masterpiece.R - NATURE/SCIENCE COVER-QUALITY VISUALIZATION
# ============================================================================
#
# DESIGN PHILOSOPHY: "Coral Reef Nebula"
#
# The challenge: 4 guilds with VERY different sizes (17 vs 4 species)
# The solution: An abstract, elegant representation where:
#   - Each guild occupies its own orbital ring around a central "coral"
#   - Species are positioned by their connectivity (degree)
#   - Ring size is VISUALLY BALANCED, not proportional to species count
#   - Small guilds get MORE space per node, making them equally prominent
#
# Visual metaphor: Like planets around a star, or reef organisms
# radiating outward from a coral host
#
# COLOR PALETTE: Deep ocean tones with bioluminescent accents
#   - Deep navy background creating depth
#   - Guilds: coral pink, tropical teal, sunset gold, deep violet
#   - Edges: gossamer silver threads
#
# Author: AI-assisted figure design
# Date: 2026-01-28
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("  CREATING MASTERPIECE: Coral Reef Nebula Network              \n")
cat("================================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

library(igraph)
library(ggforce)
library(ggrepel)
library(patchwork)

network_results <- load_object("cafi_network")
centrality_df <- network_results$centrality
edge_list <- network_results$edge_list

fig_dir <- file.path(PATHS$figures, "06_network")
manuscript_dir <- file.path(PATHS$figures, "manuscript")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)

cat("[OK] Setup complete\n\n")

# ============================================================================
# FILTER DATA AT r > 0.7 THRESHOLD
# ============================================================================

threshold <- 0.7
edges_filtered <- edge_list %>% filter(correlation > threshold)
species_in <- unique(c(edges_filtered$sp1, edges_filtered$sp2))

# Build filtered graph and recalculate degree
g_filtered <- graph_from_data_frame(
  edges_filtered %>% dplyr::select(sp1, sp2) %>% rename(from = sp1, to = sp2),
  directed = FALSE,
  vertices = data.frame(name = species_in)
)
filtered_degree <- igraph::degree(g_filtered)

# Get species info with filtered degree
species_info <- centrality_df %>%
  filter(species %in% species_in) %>%
  dplyr::select(species, type, module, occurrence) %>%
  mutate(
    degree = filtered_degree[species],
    module = as.factor(module)
  ) %>%
  arrange(module, desc(degree))

guild_lookup <- setNames(as.integer(species_info$module), species_info$species)

# Add guild info to edges
edges_filtered <- edges_filtered %>%
  mutate(
    guild1 = guild_lookup[sp1],
    guild2 = guild_lookup[sp2],
    same_guild = guild1 == guild2
  )

cat(sprintf("Network stats (r > %.1f):\n", threshold))
cat(sprintf("  Species: %d\n", length(species_in)))
cat(sprintf("  Edges: %d\n", nrow(edges_filtered)))
cat(sprintf("  Within-guild: %d\n", sum(edges_filtered$same_guild)))
cat(sprintf("  Between-guild: %d\n\n", sum(!edges_filtered$same_guild)))

# ============================================================================
# DESIGN PARAMETERS: CORAL REEF NEBULA
# ============================================================================

# Deep ocean color palette - sophisticated, nature-inspired
bg_color <- "#0A1628"          # Deep ocean night
bg_color_light <- "#FAFCFF"    # Alternative: soft paper white

# Guild colors - bioluminescent reef palette
guild_palette <- c(
  "1" = "#FF6B6B",   # Living coral / warm red
  "2" = "#4ECDC4",   # Tropical lagoon teal
  "3" = "#FFE66D",   # Bioluminescent gold (FOCAL - fish & Trapezia)
  "4" = "#C49BD8"    # Deep purple anemone
)

# Lighter versions for fills
guild_palette_soft <- c(
  "1" = "#FF6B6B30",
  "2" = "#4ECDC430",
  "3" = "#FFE66D40",
  "4" = "#C49BD830"
)

# Edge colors
edge_within_color <- "#FFFFFF20"   # Gossamer white
edge_between_color <- "#FFFFFF10"  # Even more subtle

# Type symbols
type_shapes <- c(
  "shrimp" = 21, "crab" = 22, "hermit" = 23, "fish" = 24,
  "echinoderm" = 25, "worm" = 21, "snail" = 22, "amphipod" = 23,
  "squat_lobster" = 24, "other" = 21
)

# ============================================================================
# ORBITAL RING LAYOUT - Equal visual weight per guild
# ============================================================================

# Key insight: Larger rings for larger guilds, but SAME SPACING between nodes
# This means small guilds look just as important (species spread apart)

guild_sizes <- species_info %>%
  group_by(module) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(module)

cat("Guild sizes:\n")
print(guild_sizes)

# Orbital ring radii - visual balance, not data-proportional
# Rings 1 & 2 (large guilds): larger radius, species closer together
# Rings 3 & 4 (small guilds): smaller radius but same apparent density
orbital_radii <- c(
  "1" = 3.2,   # 13 species
  "2" = 4.2,   # 17 species
  "3" = 2.2,   # 4 species (inner orbit - "close companions")
  "4" = 1.4    # 4 species (innermost - "guard crabs" concept)
)

# ============================================================================
# POSITION SPECIES ON ORBITAL RINGS
# ============================================================================

position_species <- function(species_info, orbital_radii) {
  positions <- data.frame()

  for (guild in 1:4) {
    guild_sp <- species_info %>%
      filter(module == guild) %>%
      arrange(desc(degree))  # Hub species at "prime" positions

    n <- nrow(guild_sp)
    if (n == 0) next

    r <- orbital_radii[as.character(guild)]

    # Angular positions - evenly spaced, with hub species at key angles
    # Start angle offset for visual separation between guilds
    start_angle <- switch(as.character(guild),
      "1" = pi/6,        # Top-right arc
      "2" = pi + pi/6,   # Bottom-left arc
      "3" = -pi/3,       # Right side
      "4" = 2*pi/3       # Left side
    )

    # Spread across ~120 degrees for large guilds, ~90 for small
    arc_span <- ifelse(n > 8, pi * 0.7, pi * 0.5)

    angles <- seq(start_angle, start_angle + arc_span, length.out = n)

    guild_sp$x <- r * cos(angles)
    guild_sp$y <- r * sin(angles)
    guild_sp$angle <- angles

    positions <- bind_rows(positions, guild_sp)
  }

  return(positions)
}

all_positions <- position_species(species_info, orbital_radii)

cat(sprintf("\nPositioned %d species on orbital rings\n", nrow(all_positions)))

# ============================================================================
# FORMAT SPECIES NAMES (Abbreviated binomials)
# ============================================================================

format_species_name <- function(name) {
  parts <- strsplit(name, "_")[[1]]
  if (length(parts) >= 2) {
    # Genus_species -> G. species
    genus_init <- substr(parts[1], 1, 1)
    epithet <- parts[2]
    if (nchar(epithet) > 10) epithet <- paste0(substr(epithet, 1, 9), ".")
    return(paste0(genus_init, ". ", epithet))
  }
  # Single name (family/order level)
  if (nchar(name) > 12) return(paste0(substr(name, 1, 11), "."))
  return(name)
}

all_positions <- all_positions %>%
  mutate(
    label = sapply(species, format_species_name),
    module_char = as.character(module)
  )

# ============================================================================
# BUILD EDGE DATAFRAMES
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

cat(sprintf("Edges: %d within, %d between guilds\n",
            nrow(within_edges), nrow(between_edges)))

# ============================================================================
# SELECT HUB SPECIES FOR LABELS (top 3 per guild)
# ============================================================================

hub_labels <- all_positions %>%
  group_by(module) %>%
  slice_max(degree, n = 3, with_ties = FALSE) %>%
  ungroup()

cat(sprintf("Labeling %d hub species\n\n", nrow(hub_labels)))

# ============================================================================
# ORBITAL RING ANNOTATIONS
# ============================================================================

ring_data <- data.frame(
  guild = factor(1:4),
  r = unname(orbital_radii),
  n_species = guild_sizes$n,
  label = c(
    "GUILD 1: Hermits & Shrimp",
    "GUILD 2: Worms & Echinoderms",
    "GUILD 3: Fish & Guard Crabs",
    "GUILD 4: Coral Associates"
  )
)

# Position labels at the end of each arc
ring_labels <- ring_data %>%
  mutate(
    angle = c(pi/6 + pi*0.7 + 0.15,
              pi + pi/6 + pi*0.7 + 0.15,
              -pi/3 + pi*0.5 + 0.15,
              2*pi/3 + pi*0.5 + 0.15),
    x = r * cos(angle) * 1.12,
    y = r * sin(angle) * 1.12,
    hjust = ifelse(cos(angle) > 0, 0, 1)
  )

# ============================================================================
# CREATE THE MASTERPIECE - DARK VERSION
# ============================================================================

cat("Creating dark theme masterpiece...\n")

p_dark <- ggplot() +

  # LAYER 1: Subtle orbital ring guides
  geom_circle(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r = r, color = guild),
    fill = NA,
    linewidth = 0.3,
    alpha = 0.25,
    linetype = "dotted"
  ) +

  # LAYER 2: Between-guild edges (very subtle curves)
  geom_curve(
    data = between_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "#FFFFFF",
    linewidth = 0.2,
    alpha = 0.08,
    curvature = 0.2
  ) +

  # LAYER 3: Within-guild edges (gossamer threads)
  geom_segment(
    data = within_edges,
    aes(x = x, y = y, xend = xend, yend = yend, color = guild),
    linewidth = 0.35,
    alpha = 0.25
  ) +

  # LAYER 4: Species nodes - outer glow effect
  geom_point(
    data = all_positions,
    aes(x = x, y = y, color = module_char),
    size = all_positions$degree / 2.5 + 4,
    alpha = 0.15
  ) +

  # LAYER 5: Species nodes - main
  geom_point(
    data = all_positions,
    aes(x = x, y = y, fill = module_char, size = degree),
    shape = 21,
    color = "#FFFFFF",
    stroke = 0.8,
    alpha = 0.95
  ) +

  # LAYER 6: Hub species labels
  geom_text_repel(
    data = hub_labels,
    aes(x = x, y = y, label = label, color = module_char),
    size = 3.5,
    fontface = "italic",
    max.overlaps = 20,
    force = 3,
    force_pull = 0.5,
    segment.color = "#FFFFFF50",
    segment.size = 0.3,
    segment.alpha = 0.6,
    point.padding = 0.6,
    box.padding = 0.4,
    bg.color = "#0A162880",
    bg.r = 0.15,
    show.legend = FALSE
  ) +

  # LAYER 7: Guild labels (subtle, at arc ends)
  geom_text(
    data = ring_labels,
    aes(x = x, y = y, label = paste0("Guild ", guild, "\n(n=", n_species, ")"),
        color = guild, hjust = hjust),
    size = 3.2,
    fontface = "bold",
    lineheight = 0.9,
    show.legend = FALSE
  ) +

  # LAYER 8: Central "coral" indicator
  annotate("point", x = 0, y = 0, size = 8, shape = 21,
           fill = "#FF6B6B20", color = "#FF6B6B50", stroke = 2) +
  annotate("text", x = 0, y = -0.5, label = "CORAL\nHOST",
           size = 2.5, color = "#FFFFFF60", fontface = "bold", lineheight = 0.85) +

  # SCALES
  scale_fill_manual(values = guild_palette, name = "Guild") +
  scale_color_manual(values = guild_palette, guide = "none") +
  scale_size_continuous(range = c(4, 18), name = "Degree\n(connections)") +

  # LABELS
  labs(
    title = "Coral-Associated Fauna Form Four Ecological Guilds",
    subtitle = sprintf("Co-occurrence network at r > 0.7 | %d species, %d associations | Node size = network degree",
                       nrow(all_positions), nrow(edges_filtered)),
    caption = "Guilds detected via Louvain community detection on volume-corrected correlations.\nOrbital layout emphasizes guild separation; connections show significant positive co-occurrence."
  ) +

  # THEME - Dark ocean
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = bg_color, color = NA),
    panel.background = element_rect(fill = bg_color, color = NA),
    plot.title = element_text(
      color = "#FFFFFF",
      face = "bold",
      size = 18,
      hjust = 0.5,
      margin = margin(b = 8)
    ),
    plot.subtitle = element_text(
      color = "#FFFFFF90",
      size = 11,
      hjust = 0.5,
      margin = margin(b = 20)
    ),
    plot.caption = element_text(
      color = "#FFFFFF60",
      size = 9,
      hjust = 0,
      lineheight = 1.3,
      margin = margin(t = 25)
    ),
    legend.position = c(0.92, 0.25),
    legend.background = element_rect(fill = "#FFFFFF10", color = NA),
    legend.title = element_text(color = "#FFFFFF", size = 10, face = "bold"),
    legend.text = element_text(color = "#FFFFFF90", size = 9),
    legend.key = element_rect(fill = NA, color = NA),
    plot.margin = margin(25, 25, 25, 25)
  ) +
  coord_fixed(xlim = c(-6, 6), ylim = c(-5.5, 5.5)) +
  guides(
    fill = guide_legend(override.aes = list(size = 6), order = 1),
    size = guide_legend(override.aes = list(fill = "#FFFFFF50"), order = 2)
  )

# Save dark version
ggsave(file.path(fig_dir, "fig4_masterpiece_dark.png"), p_dark,
       width = 13, height = 11, dpi = 300, bg = bg_color)

cat("  Saved: fig4_masterpiece_dark.png\n")

# ============================================================================
# CREATE LIGHT VERSION (for journal submission)
# ============================================================================

cat("Creating light theme version...\n")

# Adjusted palette for light background
guild_palette_light <- c(
  "1" = "#C73E3E",   # Deeper coral red
  "2" = "#2A8A83",   # Deeper teal
  "3" = "#D4A012",   # Rich gold
  "4" = "#8B5DA5"    # Rich purple
)

p_light <- ggplot() +

  # LAYER 1: Orbital ring guides
  geom_circle(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r = r, color = guild),
    fill = NA,
    linewidth = 0.4,
    alpha = 0.3,
    linetype = "dotted"
  ) +

  # LAYER 2: Between-guild edges
  geom_curve(
    data = between_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "#00000015",
    linewidth = 0.25,
    curvature = 0.2
  ) +

  # LAYER 3: Within-guild edges
  geom_segment(
    data = within_edges,
    aes(x = x, y = y, xend = xend, yend = yend, color = guild),
    linewidth = 0.4,
    alpha = 0.35
  ) +

  # LAYER 4: Species nodes
  geom_point(
    data = all_positions,
    aes(x = x, y = y, fill = module_char, size = degree),
    shape = 21,
    color = "#333333",
    stroke = 0.7,
    alpha = 0.9
  ) +

  # LAYER 5: Hub species labels
  geom_text_repel(
    data = hub_labels,
    aes(x = x, y = y, label = label, color = module_char),
    size = 3.3,
    fontface = "italic",
    max.overlaps = 20,
    force = 3,
    segment.color = "#55555580",
    segment.size = 0.25,
    point.padding = 0.5,
    box.padding = 0.35,
    bg.color = "#FFFFFF",
    bg.r = 0.1,
    show.legend = FALSE
  ) +

  # LAYER 6: Guild labels
  geom_text(
    data = ring_labels,
    aes(x = x, y = y, label = paste0("Guild ", guild, "\n(n=", n_species, ")"),
        color = guild, hjust = hjust),
    size = 3,
    fontface = "bold",
    lineheight = 0.9,
    show.legend = FALSE
  ) +

  # LAYER 7: Central coral indicator
  annotate("point", x = 0, y = 0, size = 6, shape = 21,
           fill = "#C73E3E20", color = "#C73E3E50", stroke = 1.5) +
  annotate("text", x = 0, y = -0.45, label = "CORAL\nHOST",
           size = 2.2, color = "#66666680", fontface = "bold", lineheight = 0.85) +

  # SCALES
  scale_fill_manual(values = guild_palette_light, name = "Guild") +
  scale_color_manual(values = guild_palette_light, guide = "none") +
  scale_size_continuous(range = c(4, 16), name = "Degree\n(connections)") +

  # LABELS
  labs(
    title = "Coral-Associated Fauna Form Four Ecological Guilds",
    subtitle = sprintf("Co-occurrence network (r > 0.7) | %d species, %d associations | Node size proportional to degree centrality",
                       nrow(all_positions), nrow(edges_filtered)),
    caption = "Guilds identified via Louvain community detection on volume-corrected Spearman correlations. Orbital layout illustrates guild structure;\nlines connect significantly co-occurring species within (colored) and between (gray) guilds."
  ) +

  # THEME - Light academic
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "#FAFCFF", color = NA),
    panel.background = element_rect(fill = "#FAFCFF", color = NA),
    plot.title = element_text(
      color = "#1a1a1a",
      face = "bold",
      size = 16,
      hjust = 0.5,
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      color = "#555555",
      size = 10,
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    plot.caption = element_text(
      color = "#777777",
      size = 8.5,
      hjust = 0,
      lineheight = 1.3,
      margin = margin(t = 20)
    ),
    legend.position = c(0.92, 0.22),
    legend.background = element_rect(fill = "#FFFFFF", color = "#DDDDDD", linewidth = 0.3),
    legend.title = element_text(color = "#333333", size = 9, face = "bold"),
    legend.text = element_text(color = "#555555", size = 8),
    legend.margin = margin(8, 8, 8, 8),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_fixed(xlim = c(-6, 6), ylim = c(-5.5, 5.5)) +
  guides(
    fill = guide_legend(override.aes = list(size = 5), order = 1),
    size = guide_legend(override.aes = list(fill = "#88888880"), order = 2)
  )

# Save light version
ggsave(file.path(fig_dir, "fig4_masterpiece_light.png"), p_light,
       width = 13, height = 11, dpi = 300, bg = "#FAFCFF")

ggsave(file.path(manuscript_dir, "fig4_network.png"), p_light,
       width = 13, height = 11, dpi = 300, bg = "#FAFCFF")

# Also save PDF for journal
ggsave(file.path(manuscript_dir, "fig4_network.pdf"), p_light,
       width = 13, height = 11, device = cairo_pdf, bg = "#FAFCFF")

cat("  Saved: fig4_masterpiece_light.png\n")
cat("  Saved: fig4_network.png (manuscript)\n")
cat("  Saved: fig4_network.pdf (manuscript)\n")

# ============================================================================
# ALTERNATIVE: COMPACT 4-QUADRANT VIEW WITH STATISTICS
# ============================================================================

cat("\nCreating compact alternative with statistics panel...\n")

# Statistics summary for annotation
stats_text <- sprintf(
  paste0("NETWORK SUMMARY\n\n",
         "Species: %d\nEdges: %d (r > 0.7)\n\n",
         "Within-guild: %d (%.0f%%)\nBetween-guild: %d (%.0f%%)\n\n",
         "Guild 1 (Hermits): %d sp.\nGuild 2 (Inverts): %d sp.\n",
         "Guild 3 (Fish+): %d sp.\nGuild 4 (Crabs): %d sp."),
  nrow(all_positions), nrow(edges_filtered),
  sum(edges_filtered$same_guild),
  100 * sum(edges_filtered$same_guild) / nrow(edges_filtered),
  sum(!edges_filtered$same_guild),
  100 * sum(!edges_filtered$same_guild) / nrow(edges_filtered),
  guild_sizes$n[1], guild_sizes$n[2], guild_sizes$n[3], guild_sizes$n[4]
)

# Main network panel (same as light version)
p_network <- p_light +
  theme(
    plot.title = element_text(size = 14),
    plot.subtitle = element_text(size = 9),
    plot.caption = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

# Legend panel
p_legend <- ggplot() +
  # Guild color legend
  geom_point(
    data = data.frame(
      guild = factor(1:4),
      y = 4:1,
      label = c("Guild 1: Hermits & Shrimp",
                "Guild 2: Worms & Echinoderms",
                "Guild 3: Fish & Guard Crabs",
                "Guild 4: Coral Associates")
    ),
    aes(x = 0.5, y = y, fill = guild),
    size = 6, shape = 21, color = "#333333", stroke = 0.5
  ) +
  geom_text(
    data = data.frame(
      guild = factor(1:4),
      y = 4:1,
      label = c("Guild 1: Hermits & Shrimp (n=13)",
                "Guild 2: Worms & Echinoderms (n=17)",
                "Guild 3: Fish & Guard Crabs (n=4)",
                "Guild 4: Coral Associates (n=4)")
    ),
    aes(x = 1, y = y, label = label),
    hjust = 0, size = 3.5, color = "#333333"
  ) +
  scale_fill_manual(values = guild_palette_light, guide = "none") +
  xlim(0, 5) + ylim(0, 5) +
  theme_void() +
  theme(plot.margin = margin(20, 10, 20, 10))

# Combine
p_combined <- p_network +
  inset_element(p_legend, 0.72, 0.02, 0.98, 0.35)

ggsave(file.path(fig_dir, "fig4_masterpiece_annotated.png"), p_combined,
       width = 13, height = 11, dpi = 300, bg = "#FAFCFF")

cat("  Saved: fig4_masterpiece_annotated.png\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("  MASTERPIECE COMPLETE                                         \n")
cat("================================================================\n\n")

cat("Files created:\n")
cat("  1. fig4_masterpiece_dark.png    - Cover-quality dark theme\n")
cat("  2. fig4_masterpiece_light.png   - Journal-ready light theme\n")
cat("  3. fig4_masterpiece_annotated.png - With embedded legend\n")
cat("  4. fig4_network.png (manuscript/)  - Final submission version\n")
cat("  5. fig4_network.pdf (manuscript/)  - Vector format for print\n\n")

cat("Design innovations:\n")
cat("  - ORBITAL LAYOUT: Species on concentric rings by guild\n")
cat("  - VISUAL BALANCE: Small guilds (4 sp.) given equal visual weight\n")
cat("  - BIOLUMINESCENT PALETTE: Deep ocean tones for elegance\n")
cat("  - HUB HIGHLIGHTING: Top 3 species per guild labeled\n")
cat("  - MINIMAL EDGE RENDERING: Focus on structure, not hairball\n")
cat("  - CENTRAL METAPHOR: Coral host as gravitational center\n\n")

cat("Why this is better than previous attempts:\n")
cat("  - Quadrant layouts made small guilds look pathetic\n")
cat("  - Force-directed layouts created illegible hairballs\n")
cat("  - This design treats all 4 guilds with equal visual respect\n")
cat("  - The orbital metaphor communicates 'community structure'\n")
cat("  - Clean enough for thumbnail; detailed enough for full-page\n\n")

cat("================================================================\n\n")
