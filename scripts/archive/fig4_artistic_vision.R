# ============================================================================
# fig4_artistic_vision.R - Magazine-Cover Quality Network Visualization
# ============================================================================
#
# PURPOSE: Create an artistic, abstract visualization of CAFI guild relationships
#          that could appear on the cover of Nature or Science
#
# DESIGN PHILOSOPHY:
#   - Abstract the network into elegant shapes representing guilds
#   - Use metaphor: coral reef bioluminescence, ocean depths
#   - Embrace negative space and minimal labeling
#   - Create visual hierarchy through size, color, and luminosity
#
# DATA:
#   - 4 ecological guilds (58 species total)
#   - Guild 1: 21 species (shrimp-dominated mobile commensals)
#   - Guild 2: 21 species (mobile cleaners with echinoderms & worms)
#   - Guild 3: 11 species (vertebrate associates - fish, Trapezia crabs)
#   - Guild 4: 5 species (peripheral crabs)
#   - 1081 edges connecting co-occurring species
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-28
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    ARTISTIC NETWORK VISUALIZATION\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))

# Required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)      # For geom_circle, geom_arc_bar, geom_bezier
library(grid)
library(gridExtra)
library(scales)

# Install ggforce if needed
if (!require("ggforce", quietly = TRUE)) {
  install.packages("ggforce", repos = "https://cran.rstudio.com")
  library(ggforce)
}

# Create output directory
fig_dir <- file.path(PATHS$figures, "06_network")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("[OK] Setup complete\n\n")

# ============================================================================
# COLOR PALETTE: Deep Ocean to Bioluminescent
# ============================================================================

# Primary palette - inspired by deep ocean bioluminescence
palette <- list(
  # Background gradient
  bg_deep = "#0a1628",       # Deep ocean midnight
  bg_mid = "#102840",        # Twilight zone

  # Guild colors - bioluminescent organisms
  guild1 = "#00d4aa",        # Bioluminescent teal (shrimp/mobile)
  guild2 = "#7b68ee",        # Violet (echinoderms/worms - cnidarian glow)
  guild3 = "#ff6b6b",        # Coral pink (fish/Trapezia - warm reef)
  guild4 = "#ffd93d",        # Golden (crabs - sunlit shallows)

  # Accent colors
  glow = "#00ffcc",          # Bright bioluminescent
  subtle = "#ffffff20",      # Subtle white
  edge_light = "#ffffff15",  # Edge glow
  text = "#e8f4f8"           # Soft white text
)

cat("[OK] Color palette defined\n")

# ============================================================================
# GUILD DATA
# ============================================================================

# Guild characteristics from the analysis
guilds <- data.frame(
  id = 1:4,
  name = c("Mobile Commensals", "Cryptic Cleaners",
           "Reef Guardians", "Peripheral Sentinels"),
  n_species = c(21, 21, 11, 5),
  total_abundance = c(1815, 436, 689, 65),
  hub_species = c("Calcinus laevimanus", "Ophiomastix elegans",
                  "Liocarpilodes integerrimus", "Paguridae"),
  dominant_type = c("shrimp", "shrimp", "fish", "crab"),
  mean_degree = c(42.3, 44.7, 28.6, 4.0)
)

# Between-guild connections (approximated from module structure)
# High connectivity within guilds 1-3, Guild 4 is peripheral
guild_connections <- data.frame(
  from = c(1, 1, 1, 2, 2, 3),
  to = c(2, 3, 4, 3, 4, 4),
  strength = c(0.85, 0.60, 0.20, 0.55, 0.15, 0.10)  # Relative connection strength
)

cat("[OK] Guild data prepared\n")

# ============================================================================
# DESIGN 1: "CONSTELLATION" - Organic circular arrangement
# ============================================================================

cat("\nCreating Design 1: Constellation...\n")

# Position guilds in an organic, slightly asymmetric arrangement
# Like celestial bodies with gravitational attraction
set.seed(42)

# Main positions - asymmetric for visual interest
guild_positions <- data.frame(
  id = 1:4,
  x = c(-2.5, 2.2, 0, 4.5),
  y = c(1.2, 1.8, -2.0, -0.5)
)

# Scale circle radius by sqrt of n_species for area-proportional representation
guilds <- guilds %>%
  left_join(guild_positions, by = "id") %>%
  mutate(
    radius = sqrt(n_species) / 2.5,  # Normalize to reasonable size
    color = c(palette$guild1, palette$guild2, palette$guild3, palette$guild4),
    glow_radius = radius * 1.3
  )

# Create bezier curves for connections between guilds
create_bezier_curve <- function(x1, y1, x2, y2, curvature = 0.3) {
  # Calculate control points for smooth curves
  mid_x <- (x1 + x2) / 2
  mid_y <- (y1 + y2) / 2

  # Perpendicular offset for curvature
  dx <- x2 - x1
  dy <- y2 - y1
  len <- sqrt(dx^2 + dy^2)

  # Rotate 90 degrees for perpendicular
  perp_x <- -dy / len * curvature * len
  perp_y <- dx / len * curvature * len

  ctrl_x <- mid_x + perp_x
  ctrl_y <- mid_y + perp_y

  data.frame(
    x = c(x1, ctrl_x, x2),
    y = c(y1, ctrl_y, y2)
  )
}

# Build connection data with bezier curves
connections_data <- guild_connections %>%
  left_join(guilds %>% dplyr::select(id, x, y, radius), by = c("from" = "id")) %>%
  rename(x1 = x, y1 = y, r1 = radius) %>%
  left_join(guilds %>% dplyr::select(id, x, y, radius), by = c("to" = "id")) %>%
  rename(x2 = x, y2 = y, r2 = radius)

# Generate bezier points for each connection
bezier_points <- do.call(rbind, lapply(1:nrow(connections_data), function(i) {
  row <- connections_data[i,]
  curve <- create_bezier_curve(row$x1, row$y1, row$x2, row$y2,
                                curvature = 0.4 * (1 - row$strength))
  curve$group <- i
  curve$strength <- row$strength
  curve
}))

# Create the constellation plot
p_constellation <- ggplot() +
  # Dark background
  theme_void() +
  theme(
    plot.background = element_rect(fill = palette$bg_deep, color = NA),
    panel.background = element_rect(fill = palette$bg_deep, color = NA),
    plot.margin = margin(20, 20, 20, 20)
  ) +

  # Connection curves (underneath)
  geom_path(
    data = bezier_points,
    aes(x = x, y = y, group = group, alpha = strength),
    color = palette$glow,
    linewidth = 1.2,
    lineend = "round"
  ) +

  # Outer glow circles (atmosphere effect)
  geom_circle(
    data = guilds,
    aes(x0 = x, y0 = y, r = glow_radius),
    fill = NA,
    color = "white",
    alpha = 0.05,
    linewidth = 8
  ) +
  geom_circle(
    data = guilds,
    aes(x0 = x, y0 = y, r = radius * 1.15),
    fill = NA,
    color = "white",
    alpha = 0.1,
    linewidth = 4
  ) +

  # Main guild circles with fill
  geom_circle(
    data = guilds,
    aes(x0 = x, y0 = y, r = radius, fill = color),
    color = "white",
    alpha = 0.85,
    linewidth = 0.5
  ) +

  # Inner highlight (3D effect)
  geom_circle(
    data = guilds,
    aes(x0 = x - radius * 0.2, y0 = y + radius * 0.2, r = radius * 0.3),
    fill = "white",
    alpha = 0.15,
    color = NA
  ) +

  # Species count inside circles
  geom_text(
    data = guilds,
    aes(x = x, y = y, label = n_species),
    color = "white",
    size = 8,
    fontface = "bold",
    family = "sans"
  ) +

  # Guild labels - elegant positioning
  geom_text(
    data = guilds,
    aes(x = x, y = y - radius - 0.5, label = name),
    color = palette$text,
    size = 3.5,
    fontface = "italic",
    family = "sans"
  ) +

  scale_fill_identity() +
  scale_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
  coord_fixed(xlim = c(-5, 7), ylim = c(-4, 4)) +

  # Title
  labs(
    title = "Coral Reef Guild Architecture",
    subtitle = "58 species across 4 ecological guilds connected by 1,081 co-occurrence associations"
  ) +
  theme(
    plot.title = element_text(
      color = palette$text,
      size = 18,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      color = alpha(palette$text, 0.7),
      size = 10,
      hjust = 0.5,
      margin = margin(b = 15)
    )
  )

ggsave(
  file.path(fig_dir, "fig4_artistic.png"),
  p_constellation,
  width = 12, height = 10, dpi = 300, bg = palette$bg_deep
)

cat("     Saved: fig4_artistic.png\n")

# ============================================================================
# DESIGN 2: "FOUR WORLDS" - Radial arrangement with organic connections
# ============================================================================

cat("\nCreating Design 2: Four Worlds...\n")

# Position guilds at cardinal-ish points with artistic offset
world_positions <- data.frame(
  id = 1:4,
  angle = c(135, 45, -60, -135) * pi / 180,  # Degrees to radians
  dist = c(3.2, 3.0, 3.5, 2.5)
) %>%
  mutate(
    x = dist * cos(angle),
    y = dist * sin(angle)
  )

guilds_world <- guilds %>%
  dplyr::select(-x, -y) %>%
  left_join(world_positions, by = "id")

# Create flowing connections using smoother bezier
connection_flows <- connections_data %>%
  rowwise() %>%
  mutate(
    # Create multiple points along curve for gradient effect
    points = list({
      n <- 50
      t <- seq(0, 1, length.out = n)

      # Cubic bezier with control points
      ctrl_offset <- 0.4 + 0.3 * (1 - strength)
      mid_x <- (x1 + x2) / 2
      mid_y <- (y1 + y2) / 2
      dx <- x2 - x1
      dy <- y2 - y1

      # Control point perpendicular to line
      ctrl_x <- mid_x - dy * ctrl_offset
      ctrl_y <- mid_y + dx * ctrl_offset

      # Bezier interpolation
      bx <- (1-t)^2 * x1 + 2*(1-t)*t * ctrl_x + t^2 * x2
      by <- (1-t)^2 * y1 + 2*(1-t)*t * ctrl_y + t^2 * y2

      data.frame(x = bx, y = by, t = t)
    })
  ) %>%
  ungroup()

flow_data <- do.call(rbind, lapply(1:nrow(connection_flows), function(i) {
  pts <- connection_flows$points[[i]]
  pts$group <- i
  pts$strength <- connection_flows$strength[i]
  pts$from <- connection_flows$from[i]
  pts$to <- connection_flows$to[i]
  pts
}))

# Create "nebula" effect - scattered points representing species
set.seed(123)
species_particles <- do.call(rbind, lapply(1:4, function(g) {
  guild <- guilds_world[g,]
  n <- guild$n_species * 3  # More particles for visual effect

  # Generate points in a gaussian distribution around guild center
  r <- abs(rnorm(n, 0, guild$radius * 0.6))
  theta <- runif(n, 0, 2*pi)

  data.frame(
    x = guild$x + r * cos(theta),
    y = guild$y + r * sin(theta),
    guild = g,
    size = runif(n, 0.5, 2),
    alpha = runif(n, 0.3, 0.9)
  )
}))

p_worlds <- ggplot() +
  theme_void() +
  theme(
    plot.background = element_rect(fill = palette$bg_deep, color = NA),
    panel.background = element_rect(fill = palette$bg_deep, color = NA),
    plot.margin = margin(30, 30, 30, 30)
  ) +

  # Flow connections with gradient alpha
  geom_path(
    data = flow_data,
    aes(x = x, y = y, group = group, alpha = strength * (1 - abs(t - 0.5) * 1.5)),
    color = palette$glow,
    linewidth = 2.5,
    lineend = "round"
  ) +

  # Species particles (nebula effect)
  geom_point(
    data = species_particles,
    aes(x = x, y = y,
        color = factor(guild),
        size = size,
        alpha = alpha),
    shape = 16
  ) +

  # Core circles for each guild
  geom_circle(
    data = guilds_world,
    aes(x0 = x, y0 = y, r = radius * 0.4, fill = color),
    color = "white",
    alpha = 0.9,
    linewidth = 0.3
  ) +

  # Center point (representing coral)
  annotate(
    "point",
    x = 0, y = 0,
    color = "#ff9966",
    size = 15,
    alpha = 0.8
  ) +
  annotate(
    "point",
    x = 0, y = 0,
    color = "#ffcc99",
    size = 8,
    alpha = 0.9
  ) +

  # Species counts
  geom_text(
    data = guilds_world,
    aes(x = x, y = y, label = n_species),
    color = "white",
    size = 6,
    fontface = "bold"
  ) +

  # Guild names at edge
  geom_text(
    data = guilds_world %>%
      mutate(
        label_x = x + sign(x) * (radius + 1),
        label_y = y + sign(y) * 0.3
      ),
    aes(x = label_x, y = label_y, label = name),
    color = palette$text,
    size = 3.2,
    fontface = "italic",
    hjust = ifelse(guilds_world$x > 0, 0, 1)
  ) +

  scale_color_manual(
    values = c(palette$guild1, palette$guild2, palette$guild3, palette$guild4),
    guide = "none"
  ) +
  scale_fill_identity() +
  scale_size_identity() +
  scale_alpha_identity() +
  coord_fixed(xlim = c(-6.5, 6.5), ylim = c(-6, 6)) +

  labs(
    title = "Four Worlds of the Coral Kingdom",
    subtitle = "Species communities orbiting their coral hosts"
  ) +
  theme(
    plot.title = element_text(
      color = palette$text,
      size = 20,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      color = alpha(palette$text, 0.6),
      size = 11,
      hjust = 0.5
    )
  )

ggsave(
  file.path(fig_dir, "fig4_four_worlds.png"),
  p_worlds,
  width = 12, height = 12, dpi = 300, bg = palette$bg_deep
)

cat("     Saved: fig4_four_worlds.png\n")

# ============================================================================
# DESIGN 3: "TIDAL FLOW" - Minimal, elegant, lots of negative space
# ============================================================================

cat("\nCreating Design 3: Tidal Flow...\n")

# Arrange guilds as flowing waves
wave_y <- function(x, phase) {
  0.8 * sin(x * 0.5 + phase)
}

# Position guilds along a gentle curve
wave_positions <- data.frame(
  id = 1:4,
  x = c(-4, -1.2, 1.8, 5)
) %>%
  mutate(y = wave_y(x, phase = 0.5))

guilds_wave <- guilds %>%
  dplyr::select(-x, -y) %>%
  left_join(wave_positions, by = "id")

# Create flowing connections (like water currents)
wave_connections <- connections_data %>%
  rowwise() %>%
  mutate(
    points = list({
      n <- 80
      t <- seq(0, 1, length.out = n)

      x_from <- guilds_wave$x[from]
      y_from <- guilds_wave$y[from]
      x_to <- guilds_wave$x[to]
      y_to <- guilds_wave$y[to]

      # Linear interpolation with wave modulation
      x <- x_from + t * (x_to - x_from)
      y <- y_from + t * (y_to - y_from) + 0.3 * sin(t * pi * 2) * strength

      data.frame(x = x, y = y, t = t)
    })
  ) %>%
  ungroup()

wave_flow_data <- do.call(rbind, lapply(1:nrow(wave_connections), function(i) {
  pts <- wave_connections$points[[i]]
  pts$group <- i
  pts$strength <- wave_connections$strength[i]
  pts
}))

p_tidal <- ggplot() +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#fafcfd", color = NA),
    panel.background = element_rect(fill = "#fafcfd", color = NA),
    plot.margin = margin(60, 40, 60, 40)
  ) +

  # Subtle wave lines in background
  annotate(
    "segment",
    x = rep(-6, 5), xend = rep(7, 5),
    y = seq(-1.5, 1.5, length.out = 5),
    yend = seq(-1.5, 1.5, length.out = 5),
    color = "#e8f4f8",
    linewidth = 0.3
  ) +

  # Flow connections
  geom_path(
    data = wave_flow_data,
    aes(x = x, y = y, group = group, alpha = strength),
    color = "#4a90a4",
    linewidth = 1.5,
    lineend = "round"
  ) +

  # Guild circles with soft shadows
  geom_circle(
    data = guilds_wave %>% mutate(x = x + 0.05, y = y - 0.08),
    aes(x0 = x, y0 = y, r = radius),
    fill = "#00000010",
    color = NA
  ) +

  # Main guild circles
  geom_circle(
    data = guilds_wave,
    aes(x0 = x, y0 = y, r = radius, fill = color),
    color = "white",
    linewidth = 2,
    alpha = 0.95
  ) +

  # Species count
  geom_text(
    data = guilds_wave,
    aes(x = x, y = y, label = n_species),
    color = "white",
    size = 10,
    fontface = "bold"
  ) +

  # Elegant guild labels below
  geom_text(
    data = guilds_wave,
    aes(x = x, y = y - radius - 0.7, label = name),
    color = "#2c3e50",
    size = 3.5,
    fontface = "plain",
    family = "sans"
  ) +

  scale_fill_identity() +
  scale_alpha_continuous(range = c(0.15, 0.5), guide = "none") +
  coord_fixed(xlim = c(-6.5, 7.5), ylim = c(-3.5, 3.5)) +

  # Minimal title
  labs(title = "Guild Flow") +
  theme(
    plot.title = element_text(
      color = "#2c3e50",
      size = 24,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 20)
    )
  )

ggsave(
  file.path(fig_dir, "fig4_tidal_flow.png"),
  p_tidal,
  width = 14, height = 8, dpi = 300, bg = "#fafcfd"
)

cat("     Saved: fig4_tidal_flow.png\n")

# ============================================================================
# DESIGN 4: "REEF ARCHITECTURE" - Abstract structural representation
# ============================================================================

cat("\nCreating Design 4: Reef Architecture...\n")

# Use arc segments to represent guilds as parts of a coral structure
arc_data <- guilds %>%
  mutate(
    start_angle = c(160, 70, -20, -70) * pi / 180,
    end_angle = start_angle + (n_species / 58) * 2 * pi * 0.8,
    r_inner = 2,
    r_outer = 2 + sqrt(total_abundance) / 20
  )

# Species dots along arcs
set.seed(456)
arc_species <- do.call(rbind, lapply(1:4, function(g) {
  arc <- arc_data[g,]
  angles <- seq(arc$start_angle, arc$end_angle, length.out = arc$n_species)
  r_vals <- runif(arc$n_species, arc$r_inner + 0.2, arc$r_outer - 0.2)

  data.frame(
    x = r_vals * cos(angles),
    y = r_vals * sin(angles),
    guild = g,
    angle = angles,
    species_num = 1:arc$n_species
  )
}))

# Connection arcs between guilds
create_connection_arc <- function(angle1, angle2, r, curvature = 1.2) {
  n <- 50
  angles <- seq(angle1, angle2, length.out = n)
  r_vals <- r * curvature * (1 - ((angles - mean(c(angle1, angle2))) / (angle2 - angle1 + 0.001))^2 * 0.3)

  data.frame(
    x = r_vals * cos(angles),
    y = r_vals * sin(angles)
  )
}

p_architecture <- ggplot() +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#0d1b2a", color = NA),
    panel.background = element_rect(fill = "#0d1b2a", color = NA),
    plot.margin = margin(40, 40, 40, 40)
  ) +

  # Central coral representation
  annotate("point", x = 0, y = 0, size = 40, color = "#1b263b", alpha = 0.8) +
  annotate("point", x = 0, y = 0, size = 25, color = "#415a77", alpha = 0.7) +
  annotate("point", x = 0, y = 0, size = 12, color = "#778da9", alpha = 0.8) +

  # Guild arcs
  geom_arc_bar(
    data = arc_data,
    aes(x0 = 0, y0 = 0, r0 = r_inner, r = r_outer,
        start = start_angle, end = end_angle, fill = color),
    alpha = 0.85,
    color = "white",
    linewidth = 0.3
  ) +

  # Species points along arcs
  geom_point(
    data = arc_species,
    aes(x = x, y = y, color = factor(guild)),
    size = 2.5,
    alpha = 0.9
  ) +

  # Guild labels at outer edge
  geom_text(
    data = arc_data %>%
      mutate(
        mid_angle = (start_angle + end_angle) / 2,
        label_r = r_outer + 0.8,
        label_x = label_r * cos(mid_angle),
        label_y = label_r * sin(mid_angle),
        hjust_val = ifelse(abs(mid_angle) < pi/2, 0, 1),
        angle_deg = mid_angle * 180 / pi,
        angle_text = ifelse(abs(angle_deg) > 90, angle_deg + 180, angle_deg)
      ),
    aes(x = label_x, y = label_y, label = paste0(name, " (", n_species, ")")),
    color = "#e0e1dd",
    size = 3.5,
    fontface = "italic"
  ) +

  # Central label
  annotate(
    "text", x = 0, y = 0,
    label = "CORAL\nHOST",
    color = "#e0e1dd",
    size = 4,
    fontface = "bold",
    lineheight = 0.9
  ) +

  scale_fill_identity() +
  scale_color_manual(
    values = c(palette$guild1, palette$guild2, palette$guild3, palette$guild4),
    guide = "none"
  ) +
  coord_fixed(xlim = c(-6, 6), ylim = c(-6, 6)) +

  labs(
    title = "Reef Architecture",
    subtitle = "58 coral-associated species organized into 4 ecological guilds"
  ) +
  theme(
    plot.title = element_text(
      color = "#e0e1dd",
      size = 22,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      color = alpha("#e0e1dd", 0.6),
      size = 11,
      hjust = 0.5
    )
  )

ggsave(
  file.path(fig_dir, "fig4_reef_architecture.png"),
  p_architecture,
  width = 11, height = 11, dpi = 300, bg = "#0d1b2a"
)

cat("     Saved: fig4_reef_architecture.png\n")

# ============================================================================
# DESIGN 5: "EMERGENCE" - The final magazine cover design
# ============================================================================

cat("\nCreating Design 5: Emergence (FINAL)...\n")

# This is the synthesis - combining the best elements:
# - Organic circular shapes from Constellation
# - Particle effects from Four Worlds
# - Negative space from Tidal Flow
# - Structural clarity from Reef Architecture

# Position guilds with intentional asymmetry and visual hierarchy
# Guild 3 (Fish/Trapezia = "stars") gets prominent position
emergence_positions <- data.frame(
  id = 1:4,
  x = c(-2.8, 2.5, 0.2, 4.2),
  y = c(0.8, 1.5, -1.8, -0.3)
)

guilds_final <- guilds %>%
  dplyr::select(-x, -y) %>%
  left_join(emergence_positions, by = "id") %>%
  mutate(
    # Emphasize Guild 3 (the "stars")
    radius_display = ifelse(id == 3, radius * 1.2, radius),
    glow_intensity = ifelse(id == 3, 0.25, 0.15)
  )

# Create organic edge flows
set.seed(789)
n_flow_lines <- 80

flow_lines <- do.call(rbind, lapply(1:nrow(connections_data), function(i) {
  row <- connections_data[i,]
  g1 <- guilds_final[guilds_final$id == row$from,]
  g2 <- guilds_final[guilds_final$id == row$to,]

  n_lines <- ceiling(row$strength * 15)

  do.call(rbind, lapply(1:n_lines, function(j) {
    # Randomize start/end within guild circles
    theta1 <- runif(1, 0, 2*pi)
    theta2 <- runif(1, 0, 2*pi)
    r1 <- runif(1, 0, g1$radius_display * 0.6)
    r2 <- runif(1, 0, g2$radius_display * 0.6)

    x1 <- g1$x + r1 * cos(theta1)
    y1 <- g1$y + r1 * sin(theta1)
    x2 <- g2$x + r2 * cos(theta2)
    y2 <- g2$y + r2 * sin(theta2)

    # Create curved path
    n_pts <- 30
    t <- seq(0, 1, length.out = n_pts)

    # Add curvature
    curve_amt <- runif(1, 0.2, 0.5) * sample(c(-1, 1), 1)
    mid_x <- (x1 + x2) / 2 - (y2 - y1) * curve_amt
    mid_y <- (y1 + y2) / 2 + (x2 - x1) * curve_amt

    # Quadratic bezier
    x <- (1-t)^2 * x1 + 2*(1-t)*t * mid_x + t^2 * x2
    y <- (1-t)^2 * y1 + 2*(1-t)*t * mid_y + t^2 * y2

    data.frame(
      x = x, y = y,
      group = paste(i, j, sep = "_"),
      strength = row$strength,
      t = t
    )
  }))
}))

# Create luminous particles around each guild
particles <- do.call(rbind, lapply(1:4, function(g) {
  guild <- guilds_final[g,]
  n_particles <- guild$n_species * 5

  # Particles in gaussian distribution, more dense at center
  r <- abs(rnorm(n_particles, 0, guild$radius_display * 0.7))
  theta <- runif(n_particles, 0, 2*pi)

  data.frame(
    x = guild$x + r * cos(theta),
    y = guild$y + r * sin(theta),
    guild = g,
    size = rexp(n_particles, 2) + 0.3,
    alpha = pmin(0.8, rexp(n_particles, 1.5) + 0.1)
  )
}))

# The final composition
p_emergence <- ggplot() +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "#050a14", color = NA),
    panel.background = element_rect(fill = "#050a14", color = NA),
    plot.margin = margin(50, 50, 50, 50)
  ) +

  # Subtle radial gradient effect (simulated)
  annotate(
    "point", x = 0, y = 0,
    size = 200, color = "#0a1525", alpha = 0.5
  ) +

  # Flow lines (underneath everything)
  geom_path(
    data = flow_lines,
    aes(x = x, y = y, group = group,
        alpha = strength * 0.3 * (1 - abs(t - 0.5) * 1.8)),
    color = palette$glow,
    linewidth = 0.3,
    lineend = "round"
  ) +

  # Particles (luminous species dots)
  geom_point(
    data = particles,
    aes(x = x, y = y, color = factor(guild), size = size, alpha = alpha),
    shape = 16
  ) +

  # Outer glow halos
  geom_circle(
    data = guilds_final,
    aes(x0 = x, y0 = y, r = radius_display * 1.4),
    fill = NA,
    color = "white",
    alpha = 0.03,
    linewidth = 15
  ) +
  geom_circle(
    data = guilds_final,
    aes(x0 = x, y0 = y, r = radius_display * 1.15),
    fill = NA,
    color = "white",
    alpha = 0.08,
    linewidth = 5
  ) +

  # Main guild cores (smaller, more focused)
  geom_circle(
    data = guilds_final,
    aes(x0 = x, y0 = y, r = radius_display * 0.5, fill = color),
    color = "white",
    alpha = 0.95,
    linewidth = 0.8
  ) +

  # Species counts
  geom_text(
    data = guilds_final,
    aes(x = x, y = y, label = n_species),
    color = "white",
    size = 7,
    fontface = "bold"
  ) +

  # Elegant minimal labels
  geom_text(
    data = guilds_final %>%
      mutate(
        label_y = y - radius_display - 0.6,
        short_name = c("Mobile", "Cryptic", "Guardians", "Peripheral")
      ),
    aes(x = x, y = label_y, label = short_name),
    color = alpha("white", 0.7),
    size = 3.5,
    fontface = "italic"
  ) +

  scale_color_manual(
    values = c(palette$guild1, palette$guild2, palette$guild3, palette$guild4),
    guide = "none"
  ) +
  scale_fill_identity() +
  scale_size_identity() +
  scale_alpha_identity() +
  coord_fixed(xlim = c(-6, 7), ylim = c(-5, 5)) +

  # Title treatment
  labs(
    title = "E M E R G E N C E",
    subtitle = "Four guilds of coral-associated fauna | 58 species | 1,081 co-occurrence edges"
  ) +
  theme(
    plot.title = element_text(
      color = "white",
      size = 28,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 8),
      family = "sans",
      letter_spacing = unit(0.3, "cm")
    ),
    plot.subtitle = element_text(
      color = alpha("white", 0.5),
      size = 10,
      hjust = 0.5
    )
  )

ggsave(
  file.path(fig_dir, "fig4_emergence.png"),
  p_emergence,
  width = 14, height = 12, dpi = 300, bg = "#050a14"
)

cat("     Saved: fig4_emergence.png\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    ARTISTIC VISUALIZATION COMPLETE\n")
cat("============================================================\n\n")

cat("Five artistic designs created:\n\n")

cat("  1. CONSTELLATION (fig4_artistic.png)\n")
cat("     - Organic circular arrangement with bezier connections\n")
cat("     - Bioluminescent color palette on deep ocean background\n\n")

cat("  2. FOUR WORLDS (fig4_four_worlds.png)\n")
cat("     - Nebula-style particle effects around guild cores\n")
cat("     - Central coral host with orbiting communities\n\n")

cat("  3. TIDAL FLOW (fig4_tidal_flow.png)\n")
cat("     - Minimal, elegant, lots of negative space\n")
cat("     - Light background for journal interior pages\n\n")

cat("  4. REEF ARCHITECTURE (fig4_reef_architecture.png)\n")
cat("     - Radial arc segments showing structural organization\n")
cat("     - Species dots along guild arcs\n\n")

cat("  5. EMERGENCE (fig4_emergence.png) <- RECOMMENDED FOR COVER\n")
cat("     - Synthesis of best elements from all designs\n")
cat("     - Luminous particles, organic flow lines, focused cores\n")
cat("     - Magazine-cover quality with dramatic visual hierarchy\n\n")

cat("All figures saved to:\n")
cat("  ", fig_dir, "\n\n")

cat("============================================================\n")
cat("    DESIGN PHILOSOPHY NOTES\n")
cat("============================================================\n\n")

cat("These visualizations prioritize:\n")
cat("  - ABSTRACTION over literal representation\n")
cat("  - NEGATIVE SPACE over visual clutter\n")
cat("  - METAPHOR (bioluminescence, reef structure) over raw data\n")
cat("  - VISUAL HIERARCHY that emphasizes Guild 3 (the 'stars')\n")
cat("  - COLOR STORY from deep ocean to luminous life\n\n")

cat("For a Nature/Science cover, recommend: fig4_emergence.png\n")
cat("For journal interior: fig4_tidal_flow.png\n\n")
