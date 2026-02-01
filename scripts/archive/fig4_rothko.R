# ============================================================================
# fig4_rothko.R - Mark Rothko-Inspired Abstract Color Field
# ============================================================================
#
# ARTISTIC VISION: Pure abstraction through color
# 4 luminous rectangles, soft edges, emotional resonance
# ZERO labels, ZERO nodes, ZERO text - pure visual communication
# Relationships through overlap, proximity, color intensity
#
# Inspired by Rothko's "multiform" paintings (1949-1970)
# "The people who weep before my pictures are having the same
#  religious experience I had when painting them." - Mark Rothko
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    ROTHKO-STYLE ABSTRACT VISUALIZATION\n")
cat("    'Color Field Ecology'\n")
cat("============================================================\n\n")

# Setup
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))

# Required packages
library(ggplot2)
library(grid)

# Create output directory
fig_dir <- file.path(PATHS$figures, "06_network")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# THE DATA: 4 GUILDS (expressed only through color)
# ============================================================================
#
# Guild 1: WORKERS (21 species) - Warm, industrious
# Guild 2: DECOMPOSERS (21 species) - Earth tones, grounded
# Guild 3: GUARDIANS (11 species) - THE STARS - Luminous, protective
# Guild 4: SPECIALISTS (5 species) - Subtle, mysterious
#
# The relationships become implied through:
# - Vertical position (hierarchy)
# - Width (importance)
# - Color luminosity (centrality)
# - Soft edges (connectivity/permeability)
# ============================================================================

# ============================================================================
# ROTHKO COLOR PHILOSOPHY
# ============================================================================
# Rothko used thin washes of color, layered to create luminosity
# Colors should appear to float, to breathe, to have depth
# The edges are never sharp - they bleed into each other

# Deep background - like being underwater in twilight
rothko_bg <- "#0A1628"  # Midnight blue-black

# The four color fields - each a different emotional register
# These are NOT random - they carry meaning:

# Guild 1 (Workers): Cadmium orange-red - warmth, energy, labor
color_workers <- c("#E85D04", "#DC2F02", "#9D0208")

# Guild 2 (Decomposers): Deep burgundy/maroon - earth, transformation
color_decomposers <- c("#7B2D26", "#5C1D16", "#3D0C0C")

# Guild 3 (Guardians - THE STARS): Luminous blue - protection, depth, the ocean
color_guardians <- c("#0096C7", "#00B4D8", "#48CAE4")

# Guild 4 (Specialists): Muted gold/ochre - rarity, preciousness
color_specialists <- c("#C9A227", "#A68A17", "#806B0F")

# ============================================================================
# HELPER: Create a soft-edged rectangle with layered color
# ============================================================================

create_rothko_rect <- function(xmin, xmax, ymin, ymax, colors, n_layers = 15,
                                edge_softness = 0.03, alpha_base = 0.85) {
  # Creates multiple overlapping rectangles with slight variations
  # This creates the soft, luminous Rothko edge effect

  rects <- data.frame()

  width <- xmax - xmin
  height <- ymax - ymin

  for (i in 1:n_layers) {
    # Each layer is slightly different size and position
    layer_factor <- (i - 1) / (n_layers - 1)  # 0 to 1

    # Outer layers are larger (creates soft edge)
    expansion <- edge_softness * width * (1 - layer_factor)

    # Color selection - cycle through the palette
    color_idx <- ((i - 1) %% length(colors)) + 1

    # Alpha decreases for outer layers
    alpha <- alpha_base * (0.3 + 0.7 * layer_factor)

    rect_df <- data.frame(
      xmin = xmin - expansion,
      xmax = xmax + expansion,
      ymin = ymin - expansion * (height/width),
      ymax = ymax + expansion * (height/width),
      fill = colors[color_idx],
      alpha = alpha,
      layer = i
    )
    rects <- rbind(rects, rect_df)
  }

  rects
}

# ============================================================================
# HELPER: Create a glowing gradient effect
# ============================================================================

create_glow <- function(cx, cy, radius, color, n_rings = 20) {
  # Creates concentric circles that fade out - like a glow
  glow_data <- data.frame()

  for (i in 1:n_rings) {
    r <- radius * (i / n_rings)
    alpha <- 0.3 * (1 - i/n_rings)^2  # Fades quadratically

    glow_data <- rbind(glow_data, data.frame(
      x = cx, y = cy, r = r,
      fill = color, alpha = alpha
    ))
  }
  glow_data
}

# ============================================================================
# COMPOSE THE PAINTING
# ============================================================================

cat("Creating Rothko-style visualization...\n")

# Canvas: vertical orientation like classic Rothko
# Height emphasis creates contemplative, spiritual feeling

# The four fields, arranged vertically
# GUARDIANS (guild 3) are at the top - the most luminous, aspirational position
# WORKERS below - foundational energy
# DECOMPOSERS - the dark earth
# SPECIALISTS - a sliver of gold at the bottom

# All rectangles expressed as proportions of canvas (0-1 scale)
all_rects <- rbind(
  # GUARDIANS (Guild 3) - Top, THE STARS
  # Wider, more prominent, the most luminous
  create_rothko_rect(
    xmin = 0.12, xmax = 0.88,
    ymin = 0.65, ymax = 0.92,
    colors = color_guardians,
    n_layers = 20,
    edge_softness = 0.04,
    alpha_base = 0.92
  ) %>% mutate(guild = 3),

  # WORKERS (Guild 1) - Upper middle
  # Warm energy supporting the guardians
  create_rothko_rect(
    xmin = 0.15, xmax = 0.85,
    ymin = 0.38, ymax = 0.60,
    colors = color_workers,
    n_layers = 18,
    edge_softness = 0.035,
    alpha_base = 0.88
  ) %>% mutate(guild = 1),

  # DECOMPOSERS (Guild 2) - Lower middle
  # Deep earth tones, grounding element
  create_rothko_rect(
    xmin = 0.18, xmax = 0.82,
    ymin = 0.15, ymax = 0.34,
    colors = color_decomposers,
    n_layers = 15,
    edge_softness = 0.03,
    alpha_base = 0.85
  ) %>% mutate(guild = 2),

  # SPECIALISTS (Guild 4) - Bottom sliver
  # Rare, precious, a glimpse of gold
  create_rothko_rect(
    xmin = 0.25, xmax = 0.75,
    ymin = 0.05, ymax = 0.11,
    colors = color_specialists,
    n_layers = 12,
    edge_softness = 0.025,
    alpha_base = 0.80
  ) %>% mutate(guild = 4)
)

# ============================================================================
# ADD SUBTLE LUMINOUS EFFECTS
# ============================================================================

# A gentle glow behind the Guardians field (they are the stars)
guardian_glow <- create_glow(
  cx = 0.5, cy = 0.785,
  radius = 0.4,
  color = "#48CAE4",
  n_rings = 25
)

# ============================================================================
# CREATE THE PAINTING
# ============================================================================

p_rothko <- ggplot() +

  # Deep background
  theme_void() +
  theme(
    plot.background = element_rect(fill = rothko_bg, color = NA),
    panel.background = element_rect(fill = rothko_bg, color = NA),
    plot.margin = margin(30, 30, 30, 30)
  ) +

  # Guardian glow (behind everything)
  geom_point(data = guardian_glow,
             aes(x = x, y = y, size = r, alpha = alpha),
             color = "#48CAE4", shape = 16) +
  scale_size_identity() +

  # The color fields - layered from outer (soft) to inner (solid)
  # Order by layer to get proper stacking
  geom_rect(data = all_rects %>% arrange(guild, layer),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = fill, alpha = alpha),
            color = NA) +

  # Use identity scales for pre-computed values
  scale_fill_identity() +
  scale_alpha_identity() +

  # Fixed aspect ratio (portrait like Rothko)
  coord_fixed(ratio = 1.2, xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)

# ============================================================================
# SAVE THE ARTWORK
# ============================================================================

cat("Saving Rothko-style figure...\n")

# Rothko paintings are meant to be large and immersive
ggsave(
  file.path(fig_dir, "fig4_rothko.png"),
  p_rothko,
  width = 14, height = 12, dpi = 300,
  bg = rothko_bg
)

cat("\n")
cat("============================================================\n")
cat("    ROTHKO VISUALIZATION COMPLETE\n")
cat("    Saved: output/figures/06_network/fig4_rothko.png\n")
cat("============================================================\n")
cat("\n")
cat("    VIEWING INSTRUCTIONS:\n")
cat("    View at arm's length in a dimly lit room.\n")
cat("    Let the colors wash over you.\n")
cat("    The four guilds are implied, not stated.\n")
cat("============================================================\n\n")
