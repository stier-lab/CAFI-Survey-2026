# ============================================================================
# fig4_topography.R - Topographic Landscape Visualization
# ============================================================================
#
# ARTISTIC VISION: Guilds as mountain ranges/islands rising from the sea
# Height = connectivity/importance
# Guild 3 (Guardians) as the highest peak
# Physical, sculptural, terrain-map aesthetic
#
# Inspired by: Relief maps, bathymetric charts, topographic art
# "The community as landscape - where ecology becomes geography"
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    TOPOGRAPHIC NETWORK VISUALIZATION\n")
cat("    'Islands of Mutualism'\n")
cat("============================================================\n\n")

# Setup
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))

# Required packages
library(ggplot2)
library(dplyr)

# Create output directory
fig_dir <- file.path(PATHS$figures, "06_network")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# THE DATA: 4 GUILDS AS GEOGRAPHIC FEATURES
# ============================================================================
#
# Guild 1: WORKERS (21 species) - Large footprint, moderate height
# Guild 2: DECOMPOSERS (21 species) - Broad base, lower elevation
# Guild 3: GUARDIANS (11 species) - THE PEAK - Highest point, dramatic
# Guild 4: SPECIALISTS (5 species) - Small island, isolated
#
# Elevation represents:
# - Mean degree centrality (connectivity)
# - Ecological importance
# - Network hub status
# ============================================================================

guilds <- data.frame(
  guild = 1:4,
  name = c("WORKERS", "DECOMPOSERS", "GUARDIANS", "SPECIALISTS"),
  subtitle = c("Hermit Crabs & Shrimp", "Worms & Echinoderms",
               "Trapezia & Hawkfish", "Rare Crabs"),
  n_species = c(21, 21, 11, 5),
  # From the actual data: mean degree by module
  mean_degree = c(43, 44, 28, 4),  # Approximate from data
  # Normalized "elevation" (0-1 scale, with artistic adjustment)
  elevation = c(0.65, 0.70, 1.0, 0.25),  # Guardians are THE PEAK
  stringsAsFactors = FALSE
)

# ============================================================================
# TOPOGRAPHIC COLOR PALETTE
# ============================================================================

# Ocean depths to mountain peaks (hypsometric tint)
topo_ocean_deep <- "#0A1628"
topo_ocean_mid <- "#1A3A5C"
topo_ocean_shallow <- "#2E6B8A"
topo_land_low <- "#4A8C6F"
topo_land_mid <- "#7DAA5E"
topo_land_high <- "#C4B264"
topo_peak <- "#F5E6D3"
topo_snow <- "#FFFFFF"

# Background: deep ocean
topo_bg <- "#0D1B2A"

# Guild-specific peak colors
peak_colors <- c(
  "#E07A5F",  # Guild 1: Workers - terracotta
  "#3D5A80",  # Guild 2: Decomposers - slate blue
  "#81B29A",  # Guild 3: Guardians - sage green (THE PEAK)
  "#F2CC8F"   # Guild 4: Specialists - sand
)

# ============================================================================
# HELPER: Generate contour lines for a "mountain"
# ============================================================================

generate_mountain_contours <- function(cx, cy, base_radius, peak_height,
                                        n_contours = 12, asymmetry = 0.15,
                                        seed = 42) {
  set.seed(seed)

  contours <- data.frame()

  for (i in 1:n_contours) {
    # Elevation fraction (0 = base, 1 = peak)
    elev_frac <- (i - 1) / (n_contours - 1)

    # Radius decreases with elevation (mountain shape)
    # Use exponential decay for more natural mountain profile
    radius <- base_radius * (1 - elev_frac)^0.7

    # Generate irregular contour
    n_points <- 100
    angles <- seq(0, 2*pi, length.out = n_points)

    # Add controlled noise for natural look
    noise <- asymmetry * radius * sin(angles * sample(2:5, 1) + runif(1, 0, 2*pi))
    noise <- noise + asymmetry * 0.5 * radius * sin(angles * sample(5:10, 1) + runif(1, 0, 2*pi))

    r_varied <- radius + noise
    r_varied[r_varied < 0.01] <- 0.01  # Prevent negative radii

    contour_df <- data.frame(
      x = cx + r_varied * cos(angles),
      y = cy + r_varied * sin(angles),
      elevation = elev_frac * peak_height,
      contour_id = i
    )
    contours <- rbind(contours, contour_df)
  }

  contours
}

# ============================================================================
# HELPER: Generate filled elevation polygons (for terrain shading)
# ============================================================================

generate_terrain_fill <- function(cx, cy, base_radius, peak_height,
                                   n_layers = 15, asymmetry = 0.12, seed = 42) {
  set.seed(seed)

  terrain <- data.frame()

  for (i in n_layers:1) {  # Draw from bottom to top
    elev_frac <- (i - 1) / (n_layers - 1)
    radius <- base_radius * (1 - elev_frac)^0.7

    n_points <- 80
    angles <- seq(0, 2*pi, length.out = n_points)

    # Consistent noise pattern across layers
    set.seed(seed + 1000)  # Same seed for smooth transitions
    noise <- asymmetry * radius * sin(angles * 3 + runif(1, 0, 2*pi))
    noise <- noise + asymmetry * 0.5 * radius * sin(angles * 7)

    r_varied <- pmax(radius + noise, 0.01)

    terrain_df <- data.frame(
      x = cx + r_varied * cos(angles),
      y = cy + r_varied * sin(angles),
      elevation = elev_frac * peak_height,
      layer = i
    )
    terrain <- rbind(terrain, terrain_df)
  }

  terrain
}

# ============================================================================
# HELPER: Generate "shadow" for 3D effect
# ============================================================================

generate_shadow <- function(cx, cy, base_radius, offset_x = 0.15, offset_y = -0.1) {
  n_points <- 80
  angles <- seq(0, 2*pi, length.out = n_points)

  shadow_df <- data.frame(
    x = cx + offset_x + base_radius * 1.05 * cos(angles),
    y = cy + offset_y + base_radius * 0.95 * sin(angles)
  )
  shadow_df
}

# ============================================================================
# POSITION THE FOUR ISLAND-MOUNTAINS
# ============================================================================

# Arrange like an archipelago
# Guardians (guild 3) in the prominent position
positions <- data.frame(
  guild = 1:4,
  cx = c(-2.5, 2.5, 0, -0.5),
  cy = c(-1.0, -1.0, 2.0, -3.5),
  base_radius = c(2.0, 2.0, 2.5, 1.0),  # Size based on n_species
  peak_height = guilds$elevation * 3,    # Scale to canvas
  stringsAsFactors = FALSE
)

# ============================================================================
# GENERATE ALL TERRAIN DATA
# ============================================================================

cat("Generating terrain data...\n")

# Shadows
all_shadows <- data.frame()
for (i in 1:4) {
  shadow <- generate_shadow(positions$cx[i], positions$cy[i],
                             positions$base_radius[i])
  shadow$guild <- i
  all_shadows <- rbind(all_shadows, shadow)
}

# Terrain fills
all_terrain <- data.frame()
for (i in 1:4) {
  terrain <- generate_terrain_fill(
    positions$cx[i], positions$cy[i],
    positions$base_radius[i],
    positions$peak_height[i],
    n_layers = 15,
    asymmetry = 0.1 + 0.05 * (i == 3),  # More variation for Guardians
    seed = 200 + i
  )
  terrain$guild <- i
  all_terrain <- rbind(all_terrain, terrain)
}

# Contour lines
all_contours <- data.frame()
for (i in 1:4) {
  contours <- generate_mountain_contours(
    positions$cx[i], positions$cy[i],
    positions$base_radius[i],
    positions$peak_height[i],
    n_contours = 10,
    asymmetry = 0.08,
    seed = 300 + i
  )
  contours$guild <- i
  all_contours <- rbind(all_contours, contours)
}

# ============================================================================
# CREATE THE HYPSOMETRIC COLOR SCALE
# ============================================================================

# Elevation-based coloring (like bathymetric/topographic maps)
# Lower = ocean colors, higher = land colors

elev_colors <- c(
  "#1A3A5C",  # Deep (0.0)
  "#2E5D7C",  # Mid-depth (0.2)
  "#4A8C6F",  # Shallow/coastal (0.4)
  "#7DAA5E",  # Low land (0.6)
  "#C4B264",  # Mid elevation (0.8)
  "#F5E6D3"   # High elevation (1.0)
)

# ============================================================================
# GENERATE OCEAN PATTERN (subtle waves)
# ============================================================================

set.seed(42)
n_wave_lines <- 30
ocean_waves <- data.frame()

for (i in 1:n_wave_lines) {
  y_base <- -5 + (i - 1) * 10 / n_wave_lines
  x_vals <- seq(-6, 6, length.out = 200)
  y_vals <- y_base + 0.1 * sin(x_vals * 2 + i * 0.5) + 0.05 * sin(x_vals * 5)

  wave_df <- data.frame(x = x_vals, y = y_vals, wave_id = i)
  ocean_waves <- rbind(ocean_waves, wave_df)
}

# ============================================================================
# CREATE THE FIGURE
# ============================================================================

cat("Creating topographic visualization...\n")

# Calculate max elevation for color scaling
max_elev <- max(all_terrain$elevation)

p_topo <- ggplot() +

  # Deep ocean background
  theme_void() +
  theme(
    plot.background = element_rect(fill = topo_bg, color = NA),
    panel.background = element_rect(fill = topo_bg, color = NA),
    plot.margin = margin(20, 20, 40, 20)
  ) +

  # Subtle ocean wave pattern
  geom_path(data = ocean_waves,
            aes(x = x, y = y, group = wave_id),
            color = topo_ocean_mid, linewidth = 0.15, alpha = 0.3) +

  # Shadows (3D effect)
  geom_polygon(data = all_shadows,
               aes(x = x, y = y, group = guild),
               fill = "#000000", alpha = 0.4) +

  # Terrain fills - each layer as a polygon
  # Color by elevation (hypsometric tinting)
  geom_polygon(data = all_terrain,
               aes(x = x, y = y, group = interaction(guild, layer),
                   fill = elevation / max_elev),
               color = NA, alpha = 0.9) +

  # Contour lines
  geom_path(data = all_contours,
            aes(x = x, y = y, group = interaction(guild, contour_id)),
            color = "#2C3E50", linewidth = 0.25, alpha = 0.5) +

  # Peak markers
  geom_point(data = positions,
             aes(x = cx, y = cy + peak_height * 0.3),
             shape = 24, size = 4, fill = topo_peak, color = "#2C3E50") +

  # Elevation scale (hypsometric)
  scale_fill_gradientn(
    colors = elev_colors,
    values = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    name = "Elevation",
    guide = "none"
  ) +

  # Guild labels (like place names on a map)
  geom_text(data = positions,
            aes(x = cx, y = cy - base_radius - 0.4,
                label = paste0(guilds$name, "\n", guilds$n_species, " spp.")),
            family = "sans", size = 3, color = "#E0E1DD",
            fontface = "italic", lineheight = 0.85) +

  # Title and cartouche
  annotate("text", x = 0, y = 4.5,
           label = "ARCHIPELAGO OF MUTUALISM",
           family = "sans", fontface = "bold",
           size = 6.5, color = "#E0E1DD") +

  annotate("text", x = 0, y = 4.0,
           label = "Coral-Associated Fauna Network Topology",
           family = "sans", fontface = "italic",
           size = 3.5, color = "#778DA9") +

  # Legend text
  annotate("text", x = 4.5, y = -4.2,
           label = "ELEVATION = NETWORK CENTRALITY",
           family = "sans", size = 2.5, color = "#778DA9") +

  annotate("text", x = 4.5, y = -4.5,
           label = "Peak: GUARDIANS (Trapezia crabs)",
           family = "sans", fontface = "italic",
           size = 2.2, color = "#81B29A") +

  # Scale bar
  annotate("segment", x = -5, xend = -3, y = -4.5, yend = -4.5,
           color = "#E0E1DD", linewidth = 0.5) +
  annotate("text", x = -4, y = -4.2,
           label = "Network distance", family = "sans",
           size = 2, color = "#778DA9") +

  # North arrow (decorative)
  annotate("text", x = 5, y = 4,
           label = "N", family = "sans", fontface = "bold",
           size = 4, color = "#E0E1DD") +
  annotate("segment", x = 5, xend = 5, y = 3.2, yend = 3.7,
           arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
           color = "#E0E1DD", linewidth = 0.5) +

  # Fixed coordinates

  coord_fixed(xlim = c(-6, 6), ylim = c(-5, 5))

# ============================================================================
# SAVE THE ARTWORK
# ============================================================================

cat("Saving topographic figure...\n")

ggsave(
  file.path(fig_dir, "fig4_topography.png"),
  p_topo,
  width = 14, height = 12, dpi = 300,
  bg = topo_bg
)

cat("\n")
cat("============================================================\n")
cat("    TOPOGRAPHIC VISUALIZATION COMPLETE\n")
cat("    Saved: output/figures/06_network/fig4_topography.png\n")
cat("============================================================\n")
cat("\n")
cat("    The four guilds as islands in an ecological sea.\n")
cat("    Height represents network centrality.\n")
cat("    The GUARDIANS rise highest - the keystone species.\n")
cat("============================================================\n\n")
