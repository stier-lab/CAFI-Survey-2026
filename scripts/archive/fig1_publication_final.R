# ============================================================================
# fig1_publication_final.R - Publication-Ready Figure 1 (Nature/Science Quality)
# ============================================================================
#
# Creates a polished 3-panel figure showing:
#   A: Satellite map of Mo'orea with study sites
#   B: Colony size distribution (histogram + density)
#   C: Neighborhood density distribution (histogram + density)
#
# IMPROVEMENTS over previous version:
#   - Real satellite imagery via maptiles package
#   - White scale bar and north arrow for visibility
#   - Refined label positioning
#   - Consistent typography
#   - Better color balance
#
# Author: CAFI Survey Analysis Pipeline
# VERSION: FINAL for publication
# ============================================================================

# Set PROJ library for sf (MUST be before loading sf/maptiles)
Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sf/proj")

# Load setup
source(here::here("scripts/00_setup.R"))
source(here::here("scripts/01_load_data.R"))

# Additional packages for mapping
library(sf)
library(maptiles)

cat("\n========================================\n")
cat("Creating Figure 1: Publication Final Version\n")
cat("========================================\n\n")

# ============================================================================
# DEFINE SITE COORDINATES (from user specification)
# ============================================================================

site_coords <- tibble(
  site = c("HAU", "MAT", "MRB"),
  site_name = c("Hauru", "Maatea", "Maharepa"),
  lat = c(-17.516, -17.604, -17.475),
  long = c(-149.922, -149.815, -149.817),
  habitat = c("North Shore", "East Shore", "Barrier Reef")
)

# Count corals per site
site_counts <- coral_master %>%
  group_by(site) %>%
  summarise(
    n_corals = n(),
    total_cafi = sum(total_cafi, na.rm = TRUE),
    .groups = "drop"
  )

site_data <- site_coords %>%
  left_join(site_counts, by = "site")

# ============================================================================
# COLOR PALETTE (Colorblind-safe - Okabe-Ito from user spec)
# ============================================================================

SITE_COLORS_FIG1 <- c(
  "HAU" = "#E69F00",  # Orange
  "MAT" = "#0072B2",  # Deep blue
  "MRB" = "#009E73"   # Teal-green
)

# Data colors for histograms
DATA_FILL <- "#4A90A4"        # Sophisticated blue-teal
DATA_FILL_ALPHA <- "#7EB3C4"  # Lighter version for histogram
DENSITY_LINE <- "#1A3A5C"     # Dark blue for density curve

# ============================================================================
# ENHANCED THEME FOR FIGURE 1
# ============================================================================

theme_fig1 <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      # Clean background
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
      panel.grid.minor = element_blank(),

      # Axis styling
      axis.line = element_line(color = "gray40", linewidth = 0.4),
      axis.ticks = element_line(color = "gray40", linewidth = 0.3),
      axis.text = element_text(color = "gray20", size = base_size - 1),
      axis.title = element_text(color = "gray10", size = base_size, face = "plain"),

      # Title styling - bold panel labels (matched to Panel A at size 15)
      plot.title = element_text(face = "bold", size = 15,
                                color = "black", hjust = 0,
                                margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size - 0.5, color = "gray40",
                                   hjust = 0, margin = margin(b = 8)),

      # Clean margins
      plot.margin = margin(10, 12, 8, 10)
    )
}

# ============================================================================
# PANEL A: SATELLITE MAP OF MO'OREA WITH STUDY SITES
# ============================================================================

cat("Creating Panel A: Satellite map with study sites...\n")

# Create bounding box for Mo'orea
bbox <- st_bbox(c(xmin = -149.98, ymin = -17.65, xmax = -149.70, ymax = -17.40),
                crs = st_crs(4326))

# Convert to sf for maptiles
bbox_sf <- st_as_sfc(bbox)

# Fetch satellite tiles
cat("  Fetching satellite imagery...\n")
sat_tiles <- tryCatch({
  get_tiles(bbox_sf, provider = "Esri.WorldImagery", zoom = 12, crop = TRUE)
}, error = function(e) {
  cat("  WARNING: Could not fetch satellite tiles:", e$message, "\n")
  NULL
})

# Convert site data to sf
sites_sf <- st_as_sf(site_data, coords = c("long", "lat"), crs = 4326)

# Label positioning - labels positioned near markers ON the island
# HAU: label to right
# MAT: label above-left (toward island interior) to stay ON the map
# MRB: label to right
label_positions <- tibble(
  site = c("HAU", "MAT", "MRB"),
  x_off = c(0.045, -0.055, 0.048),
  y_off_name = c(0.010, 0.045, 0.010),
  y_off_n = c(-0.012, 0.020, -0.012)  # not used in combined label
)

site_data <- site_data %>%
  left_join(label_positions, by = "site")

# Build the map
if (!is.null(sat_tiles)) {
  # With satellite background
  panel_a <- ggplot() +
    # Satellite basemap
    tidyterra::geom_spatraster_rgb(data = sat_tiles) +

    # Study sites - white border for visibility on satellite
    geom_sf(data = sites_sf, aes(fill = site),
            shape = 21, size = 6, color = "white", stroke = 2.5) +

    # Combined site name + sample size labels - single clean label with white border
    geom_label(data = site_data,
               aes(x = long + x_off, y = lat + y_off_name,
                   label = paste0(site_name, "\n(n=", n_corals, ")")),
               size = 3.0, fontface = "bold", color = "white",
               fill = alpha("gray15", 0.85), linewidth = 0.3, label.r = unit(0.15, "lines"),
               label.padding = unit(0.22, "lines"), lineheight = 0.85) +

    scale_fill_manual(values = SITE_COLORS_FIG1, guide = "none") +

    # Map bounds - slightly wider to prevent label cutoff
    coord_sf(xlim = c(-149.97, -149.71), ylim = c(-17.65, -17.40), expand = FALSE) +

    # White scale bar with shadow for visibility
    # Shadow
    annotate("segment", x = -149.949, xend = -149.859, y = -17.621, yend = -17.621,
             color = "black", linewidth = 2.0, alpha = 0.4) +
    # Main bar
    annotate("segment", x = -149.95, xend = -149.86, y = -17.62, yend = -17.62,
             color = "white", linewidth = 1.5) +
    annotate("segment", x = -149.95, xend = -149.95, y = -17.628, yend = -17.612,
             color = "white", linewidth = 1.2) +
    annotate("segment", x = -149.86, xend = -149.86, y = -17.628, yend = -17.612,
             color = "white", linewidth = 1.2) +
    # Scale bar text with shadow for visibility
    annotate("text", x = -149.904, y = -17.601, label = "5 km",
             size = 3.0, color = "black", fontface = "bold", alpha = 0.5) +
    annotate("text", x = -149.905, y = -17.602, label = "5 km",
             size = 3.0, color = "white", fontface = "bold") +

    # White north arrow (top right) with shadow
    annotate("text", x = -149.734, y = -17.424, label = "N",
             size = 4.5, fontface = "bold", color = "black", alpha = 0.4) +
    annotate("text", x = -149.735, y = -17.425, label = "N",
             size = 4.5, fontface = "bold", color = "white") +
    annotate("segment", x = -149.735, xend = -149.735, y = -17.46, yend = -17.435,
             arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
             color = "white", linewidth = 1.0) +

    labs(title = "A") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0,
                                margin = margin(b = 5, l = 5)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 5, 5, 5)
    )

} else {
  # Fallback: schematic map if satellite tiles unavailable
  cat("  Using schematic map (satellite unavailable)...\n")

  # Mo'orea island outline (simplified)
  moorea_coords <- matrix(c(
    -149.75, -17.46, -149.77, -17.44, -149.80, -17.44, -149.84, -17.46,
    -149.88, -17.48, -149.91, -17.51, -149.92, -17.54, -149.90, -17.58,
    -149.86, -17.60, -149.81, -17.60, -149.78, -17.57, -149.76, -17.53,
    -149.75, -17.49, -149.75, -17.46
  ), ncol = 2, byrow = TRUE)

  moorea_sf <- st_polygon(list(moorea_coords)) %>%
    st_sfc(crs = 4326) %>%
    st_sf(name = "Mo'orea")

  panel_a <- ggplot() +
    annotate("rect", xmin = -149.98, xmax = -149.70, ymin = -17.65, ymax = -17.40,
             fill = "#E8F4F8", color = NA) +
    geom_sf(data = moorea_sf, fill = "#2D7D7D", color = "#1A5050", linewidth = 0.6) +
    geom_sf(data = sites_sf, aes(fill = site),
            shape = 21, size = 6, color = "white", stroke = 2) +
    geom_label(data = site_data,
               aes(x = long + x_off, y = lat + y_off_name, label = site_name),
               size = 3.0, fontface = "bold", color = "gray10",
               fill = alpha("white", 0.85), label.size = 0) +
    geom_text(data = site_data,
              aes(x = long + x_off, y = lat + y_off_n,
                  label = paste0("(n=", n_corals, ")")),
              size = 2.5, color = "gray35") +
    scale_fill_manual(values = SITE_COLORS_FIG1, guide = "none") +
    coord_sf(xlim = c(-149.97, -149.72), ylim = c(-17.64, -17.41), expand = FALSE) +
    annotate("segment", x = -149.95, xend = -149.86, y = -17.62, yend = -17.62,
             color = "gray30", linewidth = 0.8) +
    annotate("text", x = -149.905, y = -17.605, label = "5 km",
             size = 2.5, color = "gray30") +
    labs(title = "A") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0,
                                margin = margin(b = 5, l = 5)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 5, 5, 5)
    )
}

# ============================================================================
# PANEL B: COLONY SIZE DISTRIBUTION
# ============================================================================

cat("Creating Panel B: Colony size distribution...\n")

# Calculate statistics
vol_data <- coral_master %>%
  filter(!is.na(volume), volume > 0)

vol_stats <- vol_data %>%
  summarise(
    n = n(),
    min_vol = min(volume),
    max_vol = max(volume),
    median_vol = median(volume),
    mean_vol = mean(volume),
    cv = sd(volume) / mean(volume) * 100,
    range_orders = log10(max(volume)) - log10(min(volume))
  )

cat("  Volume range:", round(vol_stats$min_vol), "to", round(vol_stats$max_vol), "cm3\n")
cat("  Range spans", round(vol_stats$range_orders, 1), "orders of magnitude\n")

# Get density max for annotation positioning
dens_obj <- density(log10(vol_data$volume))
dens_max_b <- max(dens_obj$y)

panel_b <- vol_data %>%
  ggplot(aes(x = volume)) +
  # Histogram
  geom_histogram(aes(y = after_stat(density)),
                 bins = 18, fill = DATA_FILL, color = "white",
                 alpha = 0.85, linewidth = 0.3) +
  # Density overlay
  geom_density(color = DENSITY_LINE, linewidth = 1.1, adjust = 1.1) +

  # Log scale
  scale_x_log10(
    labels = scales::label_comma(),
    breaks = c(100, 1000, 10000, 100000),
    limits = c(15, 180000)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +

  # Annotation showing range
  annotate("text", x = sqrt(vol_stats$min_vol * vol_stats$max_vol),
           y = dens_max_b * 1.15,
           label = ">3 orders\nof magnitude",
           size = 2.8, color = "gray35", fontface = "italic", lineheight = 0.9) +

  labs(
    title = "B",
    subtitle = paste0("n=", vol_stats$n, "  |  CV=", round(vol_stats$cv), "%"),
    x = expression("Colony Volume (cm"^3*")"),
    y = "Density"
  ) +
  theme_fig1() +
  theme(
    axis.title.x = element_text(margin = margin(t = 8)),
    panel.grid.major.x = element_blank()
  )

# ============================================================================
# PANEL C: NEIGHBORHOOD DENSITY DISTRIBUTION
# ============================================================================

cat("Creating Panel C: Neighborhood density distribution...\n")

# Filter to corals with neighborhood data
neighbor_data <- coral_master %>%
  filter(!is.na(n_neighbors))

neighbor_stats <- neighbor_data %>%
  summarise(
    n = n(),
    min_n = min(n_neighbors),
    max_n = max(n_neighbors),
    median_n = median(n_neighbors),
    mean_n = mean(n_neighbors),
    cv = sd(n_neighbors) / mean(n_neighbors) * 100
  )

cat("  Neighbors range:", neighbor_stats$min_n, "to", neighbor_stats$max_n, "\n")
cat("  Median:", neighbor_stats$median_n, "\n")

# Get density max
dens_max_c <- max(density(neighbor_data$n_neighbors)$y)

panel_c <- neighbor_data %>%
  ggplot(aes(x = n_neighbors)) +
  # Histogram
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 5, fill = DATA_FILL, color = "white",
                 alpha = 0.85, linewidth = 0.3) +
  # Density overlay
  geom_density(color = DENSITY_LINE, linewidth = 1.1, adjust = 1.0) +

  # Median line
  geom_vline(xintercept = neighbor_stats$median_n,
             linetype = "dashed", color = "gray45", linewidth = 0.6) +

  # Median annotation
  annotate("text", x = neighbor_stats$median_n + 3, y = dens_max_c * 0.88,
           label = paste0("median=", neighbor_stats$median_n),
           size = 2.8, color = "gray40", hjust = 0, fontface = "italic") +

  scale_x_continuous(breaks = seq(0, 80, by = 20), limits = c(-2, 85)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +

  labs(
    title = "C",
    subtitle = paste0("n=", neighbor_stats$n, "  |  CV=", round(neighbor_stats$cv), "%"),
    x = "Neighbors within 5 m",
    y = "Density"
  ) +
  theme_fig1() +
  theme(
    axis.title.x = element_text(margin = margin(t = 8)),
    panel.grid.major.x = element_blank()
  )

# ============================================================================
# COMBINE PANELS
# ============================================================================

cat("Combining panels...\n")

# Layout: A on left (map), B and C stacked on right
fig1 <- panel_a + (panel_b / panel_c) +
  plot_layout(widths = c(1.1, 1)) +
  plot_annotation(
    title = NULL,  # No overall title for journal submission
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 12, 8, 12)
    )
  )

# ============================================================================
# SAVE FIGURE
# ============================================================================

output_path <- file.path(PATHS$fig_manuscript, "fig1_study_design.png")

ggsave(
  output_path,
  fig1,
  width = 10.5,
  height = 6,
  dpi = 300,
  bg = "white"
)

cat("\n========================================\n")
cat("Figure 1 saved to:\n")
cat(output_path, "\n")
cat("========================================\n\n")

# Also save high-resolution version for print
output_path_hires <- file.path(PATHS$fig_manuscript, "fig1_study_design_600dpi.tiff")
ggsave(
  output_path_hires,
  fig1,
  width = 10.5,
  height = 6,
  dpi = 600,
  bg = "white",
  compression = "lzw"
)
cat("High-res TIFF (600dpi) saved to:\n", output_path_hires, "\n")

# Summary
cat("\nFigure specifications:\n")
cat("  - Dimensions: 10.5 x 6 inches\n")
cat("  - Resolution: 300 dpi (PNG), 600 dpi (TIFF)\n")
cat("  - Panels: A (satellite map), B (volume), C (neighbors)\n")
cat("  - Site colors: HAU=#E69F00, MAT=#0072B2, MRB=#009E73\n")
cat("  - Total corals: ", nrow(coral_master), "\n")
cat("  - Corals with neighborhood data: ", neighbor_stats$n, "\n")
