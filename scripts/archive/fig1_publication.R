# ============================================================================
# fig1_publication.R - Publication-Quality Figure 1: Landscape Heterogeneity
# ============================================================================
#
# PURPOSE: Create a polished 3-panel figure emphasizing the RANGE and
#          HETEROGENEITY of landscape variables sampled across Mo'orea reefs
#
# DESIGN PHILOSOPHY:
#   - Panel A: Satellite map showing study sites (geographic context)
#   - Panel B: Colony volume distribution (log scale) - shows 3+ orders of magnitude
#   - Panel C: Neighborhood density distribution - shows 1-76 range
#
# KEY MESSAGE: This figure emphasizes VARIATION within the sampling design,
#              NOT site-to-site comparisons. Single cohesive color palette.
#
# OUTPUT: 11 x 6.5 inches, 300 dpi (compact, journal-ready)
#
# Author: CAFI Survey Analysis
# Date: 2026-01-25
# ============================================================================

cat("\n========================================\n")
cat("FIGURE 1: Publication-Quality Landscape Heterogeneity\n")
cat("========================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

# Set PROJ_LIB for spatial operations (required for sf/terra)
Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sf/proj")

# Load setup and data
source(here::here("scripts/00_setup.R"))
source(here::here("scripts/01_load_data.R"))

# Load additional packages
library(sf)
library(terra)

# Create output directory
MANUSCRIPT_DIR <- file.path(PATHS$figures, "manuscript")
dir.create(MANUSCRIPT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# COLOR PALETTE - Cohesive, publication-ready scheme
# ============================================================================

# Primary accent color (sophisticated teal-blue)
ACCENT_PRIMARY <- "#3182BD"     # Strong blue
ACCENT_SECONDARY <- "#6BAED6"   # Medium blue
ACCENT_LIGHT <- "#C6DBEF"       # Light blue fill
DENSITY_LINE <- "#08519C"       # Dark blue for density curves

# Site colors (for map markers only - Okabe-Ito)
SITE_COLORS <- c(
  "HAU" = "#E69F00",  # Orange
  "MAT" = "#0072B2",  # Blue
  "MRB" = "#009E73"   # Teal
)

cat("Data loaded successfully\n")
cat("  Corals:", nrow(coral_master), "\n")
cat("  CAFI records:", nrow(cafi_clean), "\n\n")

# ============================================================================
# CALCULATE SUMMARY STATISTICS
# ============================================================================

cat("Calculating summary statistics...\n")

# Volume statistics
vol_data <- coral_master %>% filter(!is.na(volume), volume > 0)
vol_stats <- list(
  n = nrow(vol_data),
  min = min(vol_data$volume),
  max = max(vol_data$volume),
  median = median(vol_data$volume),
  cv = sd(vol_data$volume) / mean(vol_data$volume) * 100,
  range_orders = log10(max(vol_data$volume)) - log10(min(vol_data$volume))
)

# Neighborhood statistics
neighbor_data <- coral_master %>% filter(!is.na(n_neighbors))
neigh_stats <- list(
  n = nrow(neighbor_data),
  min = min(neighbor_data$n_neighbors),
  max = max(neighbor_data$n_neighbors),
  median = median(neighbor_data$n_neighbors),
  cv = sd(neighbor_data$n_neighbors) / mean(neighbor_data$n_neighbors) * 100
)

cat("  Volume range:", round(vol_stats$min), "to",
    format(round(vol_stats$max), big.mark=","), "cm3\n")
cat("  Volume spans", round(vol_stats$range_orders, 1), "orders of magnitude\n")
cat("  Neighbors range:", neigh_stats$min, "to", neigh_stats$max,
    "(n =", neigh_stats$n, "corals)\n\n")

# ============================================================================
# ENHANCED THEME FOR FIGURE 1
# ============================================================================

# Define font family for consistency (sans-serif for publication)
FONT_FAMILY <- "Helvetica"

theme_fig1 <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      # Clean white background with subtle border
      panel.background = element_rect(fill = "white", color = "gray85", linewidth = 0.3),
      plot.background = element_rect(fill = "white", color = NA),

      # Subtle grid
      panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
      panel.grid.minor = element_blank(),

      # Clean axis styling with explicit tick length for print
      axis.line = element_line(color = "gray40", linewidth = 0.4),
      axis.ticks = element_line(color = "gray40", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text = element_text(color = "gray20", size = base_size - 1, family = FONT_FAMILY),
      axis.title = element_text(color = "gray10", size = base_size, face = "plain", family = FONT_FAMILY),

      # Panel label (B, C) - bold and prominent
      plot.title = element_text(face = "bold", size = base_size + 4,
                                color = "gray5", hjust = 0,
                                margin = margin(b = 3), family = FONT_FAMILY),
      plot.subtitle = element_text(size = base_size - 1, color = "gray50",
                                   hjust = 0, margin = margin(b = 10), family = FONT_FAMILY),

      # Optimized margins
      plot.margin = margin(10, 12, 8, 10)
    )
}

# ============================================================================
# PANEL A: SATELLITE MAP WITH STUDY SITES
# ============================================================================

cat("Creating Panel A: Satellite map...\n")

# Load satellite raster
satellite_path <- here::here("output/figures/moorea_satellite.tif")

if (file.exists(satellite_path)) {

  sat_raster <- terra::rast(satellite_path)
  sat_df <- as.data.frame(sat_raster, xy = TRUE)
  names(sat_df) <- c("x", "y", "r", "g", "b")

  # Normalize RGB
  sat_df <- sat_df %>%
    mutate(
      r = pmin(pmax(r / 255, 0), 1),
      g = pmin(pmax(g / 255, 0), 1),
      b = pmin(pmax(b / 255, 0), 1)
    )

  # Site data
  site_coords <- tibble(
    site = c("HAU", "MRB", "MAT"),
    site_name = c("Hauru", "Barrier Reef", "Maatea"),
    lat = c(-17.516, -17.475, -17.604),
    long = c(-149.922, -149.817, -149.815),
    n = c(
      sum(coral_master$site == "HAU"),
      sum(coral_master$site == "MRB"),
      sum(coral_master$site == "MAT")
    )
  )

  # Calculate bounds
  x_range <- range(sat_df$x, na.rm = TRUE)
  y_range <- range(sat_df$y, na.rm = TRUE)

  # Scale bar: 5 km at this latitude ~ 0.045 degrees longitude
  scale_bar_length <- 0.045  # ~5km
  scale_bar_x_start <- x_range[1] + 0.015
  scale_bar_y <- y_range[1] + 0.015

  panel_a <- ggplot() +
    # Satellite imagery
    geom_raster(data = sat_df, aes(x = x, y = y, fill = rgb(r, g, b))) +
    scale_fill_identity() +

    # Site markers - clean white-bordered circles
    geom_point(data = site_coords,
               aes(x = long, y = lat),
               shape = 21, size = 4.5, stroke = 1.8,
               fill = SITE_COLORS[site_coords$site], color = "white") +

    # Site labels with semi-transparent background
    # HAU (Hauru): west coast - nudge right
    # MRB (Barrier Reef): north - nudge right
    # MAT (Maatea): south tip - nudge up and left to stay on island
    geom_label(data = site_coords,
               aes(x = long, y = lat,
                   label = paste0(site_name, "\n(n=", n, ")")),
               size = 2.6, fontface = "bold", family = FONT_FAMILY,
               fill = alpha("white", 0.92),
               label.padding = unit(0.18, "lines"),
               label.r = unit(0.12, "lines"),
               linewidth = 0,
               nudge_x = c(0.045, 0.045, -0.042),
               nudge_y = c(0.012, 0.018, 0.022)) +

    # Scale bar - single clean line with end caps
    annotate("segment",
             x = scale_bar_x_start, xend = scale_bar_x_start + scale_bar_length,
             y = scale_bar_y, yend = scale_bar_y,
             color = "white", linewidth = 1.5) +
    # Left cap
    annotate("segment",
             x = scale_bar_x_start, xend = scale_bar_x_start,
             y = scale_bar_y - 0.003, yend = scale_bar_y + 0.003,
             color = "white", linewidth = 1) +
    # Right cap
    annotate("segment",
             x = scale_bar_x_start + scale_bar_length,
             xend = scale_bar_x_start + scale_bar_length,
             y = scale_bar_y - 0.003, yend = scale_bar_y + 0.003,
             color = "white", linewidth = 1) +
    # Scale label
    annotate("text",
             x = scale_bar_x_start + scale_bar_length/2,
             y = scale_bar_y + 0.012,
             label = "5 km", color = "white", size = 2.8, fontface = "bold",
             family = FONT_FAMILY) +

    # North arrow (compact)
    annotate("text", x = x_range[2] - 0.02, y = y_range[2] - 0.015,
             label = "N", color = "white", size = 3.5, fontface = "bold",
             family = FONT_FAMILY) +
    annotate("segment",
             x = x_range[2] - 0.02, xend = x_range[2] - 0.02,
             y = y_range[2] - 0.035, yend = y_range[2] - 0.02,
             arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
             color = "white", linewidth = 0.8) +

    coord_fixed(ratio = 1) +
    labs(title = "A") +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0,
                                color = "gray5", family = FONT_FAMILY,
                                margin = margin(b = 4, l = 5)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 8, 8, 8)
    )

  cat("  Satellite map created successfully\n")

} else {
  stop("Satellite image not found at: ", satellite_path)
}

# ============================================================================
# PANEL B: COLONY VOLUME DISTRIBUTION (LOG SCALE)
# ============================================================================

cat("Creating Panel B: Volume distribution...\n")

# Calculate density for y-axis scaling
vol_density <- density(log10(vol_data$volume))
y_max_b <- max(vol_density$y) * 1.15

panel_b <- ggplot(vol_data, aes(x = volume)) +
  # Histogram - use breaks aligned with log scale (linewidth optimized for print)
  geom_histogram(aes(y = after_stat(density)),
                 breaks = 10^seq(1.2, 4.7, length.out = 20),
                 fill = ACCENT_PRIMARY, color = "white",
                 alpha = 0.70, linewidth = 0.4) +
  # Density overlay - smoother curve for publication

  geom_density(color = DENSITY_LINE, linewidth = 1.1, adjust = 1.3) +
  # Rug (subtle)
  geom_rug(sides = "b", alpha = 0.25, color = ACCENT_PRIMARY, linewidth = 0.2) +

  # Range annotation - clean text box in top right
  annotate("label",
           x = 28000, y = y_max_b * 0.92,
           label = ">3 orders\nof magnitude",
           size = 2.6, color = "gray25", fontface = "italic",
           family = FONT_FAMILY,
           fill = alpha("white", 0.85), linewidth = 0,
           hjust = 1, vjust = 1, lineheight = 0.9) +

  # Log scale with clean breaks - ensure all data visible
  scale_x_log10(
    labels = scales::label_comma(),
    breaks = c(100, 1000, 10000),
    limits = c(15, 55000),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18)),
    breaks = c(0, 0.5, 1.0)
  ) +

  labs(
    title = "B",
    subtitle = paste0("n=", vol_stats$n, "  |  CV=", round(vol_stats$cv), "%"),
    x = expression("Colony Volume (cm"^3*")"),
    y = "Density"
  ) +
  theme_fig1() +
  theme(
    axis.title.x = element_text(margin = margin(t = 6)),
    panel.grid.major.x = element_blank()
  )

# ============================================================================
# PANEL C: NEIGHBORHOOD DENSITY DISTRIBUTION
# ============================================================================

cat("Creating Panel C: Neighborhood distribution...\n")

# Calculate density for y-axis scaling
neigh_density <- density(neighbor_data$n_neighbors)
y_max_c <- max(neigh_density$y) * 1.15

panel_c <- ggplot(neighbor_data, aes(x = n_neighbors)) +
  # Histogram - use boundary-aligned bins (linewidth optimized for print)
  geom_histogram(aes(y = after_stat(density)),
                 breaks = seq(0, 80, by = 5),
                 fill = ACCENT_PRIMARY, color = "white",
                 alpha = 0.70, linewidth = 0.4) +
  # Density overlay - smoother curve for publication
  geom_density(color = DENSITY_LINE, linewidth = 1.1, adjust = 1.2) +
  # Rug
  geom_rug(sides = "b", alpha = 0.30, color = ACCENT_PRIMARY, linewidth = 0.2) +

  # Median line (dashed)
  geom_vline(xintercept = neigh_stats$median,
             linetype = "dashed", color = "gray45", linewidth = 0.5) +

  # Median annotation
  annotate("text",
           x = neigh_stats$median + 2, y = y_max_c * 0.85,
           label = paste0("median=", neigh_stats$median),
           size = 2.5, color = "gray40", hjust = 0, fontface = "italic",
           family = FONT_FAMILY) +

  scale_x_continuous(
    breaks = seq(0, 80, by = 20),
    limits = c(-2, 82),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +

  labs(
    title = "C",
    subtitle = paste0("n=", neigh_stats$n, "  |  CV=", round(neigh_stats$cv), "%"),
    x = "Neighbors within 5 m",
    y = "Density"
  ) +
  theme_fig1() +
  theme(
    axis.title.x = element_text(margin = margin(t = 6)),
    panel.grid.major.x = element_blank()
  )

# ============================================================================
# COMBINE PANELS INTO FINAL FIGURE
# ============================================================================

cat("Combining panels...\n")

# Layout: A on left (1.25 width), B and C stacked on right
# Slightly reduced map width for better balance
fig1_final <- panel_a + (panel_b / panel_c) +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(
    title = "Figure 1. Landscape heterogeneity across Mo'orea reef sites",
    subtitle = expression(italic("Pocillopora")*" colonies  |  French Polynesia  |  2019"),
    caption = paste0(
      "Sites span fringing reef (Hauru, n=38), back-reef (Maatea, n=39), and barrier reef (MRB, n=35). ",
      "Colony volumes: ", round(vol_stats$min), "\u2013",
      format(round(vol_stats$max), big.mark = ","), " cm\u00B3. ",
      "Neighborhood density: ", neigh_stats$min, "\u2013", neigh_stats$max, " corals within 5 m radius."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, color = "gray5",
                                margin = margin(b = 4), family = FONT_FAMILY),
      plot.subtitle = element_text(size = 10, color = "gray40",
                                   margin = margin(b = 12), family = FONT_FAMILY),
      plot.caption = element_text(size = 8.5, color = "gray45", hjust = 0,
                                  margin = margin(t = 12), lineheight = 1.3,
                                  family = FONT_FAMILY),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(15, 18, 12, 15)
    )
  )

# ============================================================================
# SAVE FIGURE
# ============================================================================

output_path <- file.path(MANUSCRIPT_DIR, "fig1_study_design.png")

ggsave(
  output_path,
  fig1_final,
  width = 11,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

cat("\n========================================\n")
cat("FIGURE 1 COMPLETE\n")
cat("========================================\n")
cat("Saved to:", output_path, "\n\n")

cat("Design specifications:\n")
cat("  - Dimensions: 11 x 6.5 inches (Nature single-column compatible)\n")
cat("  - Resolution: 300 dpi\n")
cat("  - Layout: A (satellite map) | B (volume) / C (neighbors)\n")
cat("  - Color scheme: Cohesive blue palette (emphasizes variation)\n")
cat("  - Scale bar: Continuous with end caps (5 km)\n")
cat("  - Key data:\n")
cat("    * Total corals:", vol_stats$n, "\n")
cat("    * Volume range:", round(vol_stats$min), "-",
    format(round(vol_stats$max), big.mark=","), "cm3\n")
cat("    * Neighbor range:", neigh_stats$min, "-", neigh_stats$max, "\n\n")
