# ============================================================================
# fig4_haeckel.R - Ernst Haeckel-Inspired Scientific Illustration
# ============================================================================
#
# ARTISTIC VISION: Victorian natural history aesthetic
# Each guild becomes an intricate organic form, hand-drawn feeling
# Stippling, hatching, ornate borders, scientific elegance
#
# Inspired by Haeckel's "Kunstformen der Natur" (Art Forms in Nature, 1904)
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    HAECKEL-STYLE NETWORK VISUALIZATION\n")
cat("    'Art Forms in the Coral'\n")
cat("============================================================\n\n")

# Setup
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))

# Required packages
library(ggplot2)
library(ggforce)
library(grid)
library(gridExtra)

# Create output directory
fig_dir <- file.path(PATHS$figures, "06_network")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# THE DATA: 4 GUILDS
# ============================================================================

guilds <- data.frame(
  guild = c(1, 2, 3, 4),
  name = c("WORKERS", "DECOMPOSERS", "GUARDIANS", "SPECIALISTS"),
  full_name = c(
    "The Workers\nHermit Crabs & Shrimp",
    "The Decomposers\nWorms & Echinoderms",
    "The Guardians\nTrapezia & Hawkfish",
    "The Specialists\nRare Crabs"
  ),
  n_species = c(21, 21, 11, 5),
  hub = c("Calcinus laevimanus", "Ophiomastix elegans",
          "Liocarpilodes integerrimus", "Paguridae"),
  stringsAsFactors = FALSE
)

# Network stats
n_species_total <- 58
n_edges <- 1081  # Note: reduced display for artistic purposes

# ============================================================================
# HAECKEL COLOR PALETTE - Sepia, ochre, deep ocean tones
# ============================================================================

haeckel_bg <- "#F5F0E1"        # Aged paper
haeckel_ink <- "#2C1810"       # Deep sepia
haeckel_gold <- "#8B7355"      # Antique gold
haeckel_coral <- "#C4846C"     # Coral pink
haeckel_ocean <- "#4A6B7C"     # Deep ocean
haeckel_accent <- "#6B4423"    # Rich brown

guild_colors <- c(
  "#C4846C",   # Guild 1: Workers - coral/terracotta

  "#6B8E7C",   # Guild 2: Decomposers - sage/moss
  "#4A6B7C",   # Guild 3: Guardians - ocean blue (STARS)
  "#8B7355"    # Guild 4: Specialists - antique gold
)

# ============================================================================
# HELPER: Generate organic stippling patterns
# ============================================================================

generate_stipple_points <- function(cx, cy, radius, n_points, seed = 42) {
  set.seed(seed)
  # Generate points with density gradient (more dense toward edges for Haeckel look)
  angles <- runif(n_points, 0, 2*pi)
  # Use beta distribution to cluster points in a ring pattern
  radii <- radius * sqrt(rbeta(n_points, 1.5, 1.5))

  data.frame(
    x = cx + radii * cos(angles),
    y = cy + radii * sin(angles)
  )
}

# ============================================================================
# HELPER: Generate radial spokes/tentacles (like radiolarians)
# ============================================================================

generate_radial_form <- function(cx, cy, n_arms, inner_r, outer_r, wave_amp = 0.1) {
  # Create flowing organic tentacle-like shapes
  spokes <- data.frame()

  for (i in 1:n_arms) {
    angle <- (i - 1) * 2 * pi / n_arms

    # Create wavy line from center outward
    t <- seq(0, 1, length.out = 50)
    wave <- wave_amp * sin(t * 4 * pi) * (1 - t)  # Diminishing wave

    r <- inner_r + (outer_r - inner_r) * t
    a <- angle + wave

    spoke_df <- data.frame(
      x = cx + r * cos(a),
      y = cy + r * sin(a),
      arm = i
    )
    spokes <- rbind(spokes, spoke_df)
  }
  spokes
}

# ============================================================================
# HELPER: Generate concentric rings with organic variation
# ============================================================================

generate_organic_rings <- function(cx, cy, radii, n_points = 200, wobble = 0.02) {
  rings <- data.frame()

  for (r in radii) {
    angles <- seq(0, 2*pi, length.out = n_points)
    # Add organic wobble
    set.seed(round(r * 1000))
    r_var <- r + wobble * r * sin(angles * sample(5:12, 1))

    ring_df <- data.frame(
      x = cx + r_var * cos(angles),
      y = cy + r_var * sin(angles),
      radius = r
    )
    rings <- rbind(rings, ring_df)
  }
  rings
}

# ============================================================================
# BUILD THE COMPOSITION
# ============================================================================

cat("Creating Haeckel-style visualization...\n")

# Canvas dimensions (centered, with decorative borders)
canvas_cx <- 0
canvas_cy <- 0

# Guild positions - arranged like specimens on a scientific plate
# Central arrangement with the GUARDIANS (guild 3) prominent
positions <- data.frame(
  guild = 1:4,
  x = c(-2.8, 2.8, 0, 0),      # Workers left, Decomposers right, Guardians center-top
  y = c(-0.5, -0.5, 2.2, -3),  # Specialists at bottom
  size = c(1.8, 1.8, 2.2, 1.0) # Guardians largest
)

# ============================================================================
# GENERATE ALL GRAPHIC ELEMENTS
# ============================================================================

# 1. Stippling for each guild
all_stipples <- data.frame()
for (i in 1:4) {
  n_pts <- guilds$n_species[i] * 80  # Density proportional to species
  stipples <- generate_stipple_points(
    positions$x[i], positions$y[i],
    positions$size[i] * 0.9,
    n_pts,
    seed = 100 + i
  )
  stipples$guild <- i
  all_stipples <- rbind(all_stipples, stipples)
}

# 2. Radial spokes for each guild (like radiolarian skeletons)
all_spokes <- data.frame()
for (i in 1:4) {
  n_arms <- guilds$n_species[i]  # Arms = species count
  spokes <- generate_radial_form(
    positions$x[i], positions$y[i],
    n_arms = n_arms,
    inner_r = positions$size[i] * 0.2,
    outer_r = positions$size[i] * 1.1,
    wave_amp = 0.05 + 0.02 * (4 - i)  # More wave for larger guilds
  )
  spokes$guild <- i
  all_spokes <- rbind(all_spokes, spokes)
}

# 3. Concentric rings for each guild
all_rings <- data.frame()
for (i in 1:4) {
  n_rings <- 3 + floor(guilds$n_species[i] / 7)
  radii <- seq(positions$size[i] * 0.3, positions$size[i], length.out = n_rings)
  rings <- generate_organic_rings(positions$x[i], positions$y[i], radii, wobble = 0.03)
  rings$guild <- i
  all_rings <- rbind(all_rings, rings)
}

# 4. Connection curves between guilds (representing network connections)
# Draw elegant curved lines between guild centers
connections <- data.frame(
  from = c(1, 1, 1, 2, 2, 3),
  to = c(2, 3, 4, 3, 4, 4),
  strength = c(0.8, 0.6, 0.3, 0.7, 0.4, 0.5)  # Visual weight
)

# Generate bezier curves for connections
connection_curves <- data.frame()
for (j in 1:nrow(connections)) {
  from_g <- connections$from[j]
  to_g <- connections$to[j]

  x0 <- positions$x[from_g]
  y0 <- positions$y[from_g]
  x1 <- positions$x[to_g]
  y1 <- positions$y[to_g]

  # Control points for bezier (curved toward center)
  ctrl_x <- (x0 + x1) / 2 + (y1 - y0) * 0.2
  ctrl_y <- (y0 + y1) / 2 - (x1 - x0) * 0.2

  # Generate curve points
  t <- seq(0, 1, length.out = 100)
  curve_x <- (1-t)^2 * x0 + 2*(1-t)*t * ctrl_x + t^2 * x1
  curve_y <- (1-t)^2 * y0 + 2*(1-t)*t * ctrl_y + t^2 * y1

  curve_df <- data.frame(
    x = curve_x, y = curve_y,
    connection = j,
    strength = connections$strength[j]
  )
  connection_curves <- rbind(connection_curves, curve_df)
}

# ============================================================================
# DECORATIVE BORDER ELEMENTS
# ============================================================================

# Ornate Victorian border
border_outer <- data.frame(
  x = c(-6, 6, 6, -6, -6),
  y = c(-5.5, -5.5, 4.5, 4.5, -5.5)
)

border_inner <- data.frame(
  x = c(-5.7, 5.7, 5.7, -5.7, -5.7),
  y = c(-5.2, -5.2, 4.2, 4.2, -5.2)
)

# Corner flourishes (simplified ornamental shapes)
corner_size <- 0.4
corners <- expand.grid(
  cx = c(-5.5, 5.5),
  cy = c(-5, 4)
)

# ============================================================================
# CREATE THE FIGURE
# ============================================================================

p_haeckel <- ggplot() +

  # Background - aged paper
  theme_void() +
  theme(
    plot.background = element_rect(fill = haeckel_bg, color = NA),
    panel.background = element_rect(fill = haeckel_bg, color = NA)
  ) +

  # Double border frame
  geom_path(data = border_outer, aes(x = x, y = y),
            color = haeckel_ink, linewidth = 2) +
  geom_path(data = border_inner, aes(x = x, y = y),
            color = haeckel_ink, linewidth = 0.5) +

  # Connection curves (behind everything)
  geom_path(data = connection_curves,
            aes(x = x, y = y, group = connection, alpha = strength),
            color = haeckel_gold, linewidth = 1.5) +

  # Concentric rings for each guild
  geom_path(data = all_rings,
            aes(x = x, y = y, group = interaction(guild, radius), color = factor(guild)),
            linewidth = 0.3, alpha = 0.6) +

  # Radial spokes (like radiolarian skeletons)
  geom_path(data = all_spokes,
            aes(x = x, y = y, group = interaction(guild, arm), color = factor(guild)),
            linewidth = 0.4, alpha = 0.7) +

  # Stippling points
  geom_point(data = all_stipples,
             aes(x = x, y = y, color = factor(guild)),
             size = 0.15, alpha = 0.4) +

  # Central circles for each guild (solid core)
  geom_circle(data = positions,
              aes(x0 = x, y0 = y, r = size * 0.15, fill = factor(guild)),
              color = haeckel_ink, linewidth = 0.8, alpha = 0.9) +

  # Guild numbers in Roman numerals style
  geom_text(data = positions,
            aes(x = x, y = y, label = c("I", "II", "III", "IV")),
            family = "serif", fontface = "bold",
            size = 6, color = haeckel_bg) +

  # Guild labels below each form
  geom_text(data = positions,
            aes(x = x, y = y - size - 0.5,
                label = paste0(guilds$name, "\n", guilds$n_species, " species")),
            family = "serif", fontface = "italic",
            size = 3.2, color = haeckel_ink, lineheight = 0.85) +

  # Color scales
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_fill_manual(values = guild_colors, guide = "none") +
  scale_alpha_identity() +

  # Title (Victorian style)
  annotate("text", x = 0, y = 3.8,
           label = "CAFI COMMUNITATIS CORALLIORUM",
           family = "serif", fontface = "bold.italic",
           size = 7, color = haeckel_ink) +

  annotate("text", x = 0, y = 3.35,
           label = "The Community Forms of Coral-Associated Fauna",
           family = "serif", fontface = "italic",
           size = 4, color = haeckel_gold) +

  # Subtitle with network statistics
  annotate("text", x = 0, y = -4.8,
           label = paste0("Plate I.    N = ", n_species_total, " species    ",
                          "161 associations    4 modules"),
           family = "serif", size = 3, color = haeckel_ink) +

  # Scientific caption
  annotate("text", x = 0, y = -5.05,
           label = "Pocillopora spp., Mo'orea, French Polynesia",
           family = "serif", fontface = "italic",
           size = 2.8, color = haeckel_gold) +

  # Coordinate system

  coord_fixed(xlim = c(-6.2, 6.2), ylim = c(-5.6, 4.6))

# ============================================================================
# SAVE THE ARTWORK
# ============================================================================

cat("Saving Haeckel-style figure...\n")

ggsave(
  file.path(fig_dir, "fig4_haeckel.png"),
  p_haeckel,
  width = 14, height = 12, dpi = 300,
  bg = haeckel_bg
)

cat("\n")
cat("============================================================\n")
cat("    HAECKEL VISUALIZATION COMPLETE\n")
cat("    Saved: output/figures/06_network/fig4_haeckel.png\n")
cat("============================================================\n\n")
