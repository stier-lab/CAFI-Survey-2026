# ============================================================================
# create_fig1_landscape.R - Publication-Quality Figure 1: Landscape Heterogeneity
# ============================================================================
#
# Creates a polished 3-panel figure showing:
#   A: Satellite map of Mo'orea with study sites
#   B: Colony size distribution (histogram + density)
#   C: Neighborhood density distribution (histogram + density)
#
# Author: CAFI Survey Analysis Pipeline
# VERSION 2 - Improved aesthetics and layout
# ============================================================================

# Set PROJ library for sf
Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sf/proj")

# Load setup
source(here::here("scripts/00_setup.R"))
source(here::here("scripts/01_load_data.R"))

# Additional packages for mapping
library(sf)
library(ggspatial)

cat("\n========================================\n")
cat("Creating Figure 1: Landscape Heterogeneity (v2)\n")
cat("========================================\n\n")

# ============================================================================
# DEFINE SITE COORDINATES
# ============================================================================

site_coords <- tibble(
  site = c("HAU", "MAT", "MRB"),
  site_name = c("Hauru", "Maatea", "Barrier Reef"),
  lat = c(-17.516, -17.604, -17.475),
  long = c(-149.922, -149.815, -149.817),
  habitat = c("North Shore\nFringing Reef", "East Shore\nBack-Reef", "Offshore\nBarrier Reef")
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
# COLOR PALETTE (Colorblind-safe - Okabe-Ito derivatives)
# ============================================================================

SITE_COLORS <- c(
  "HAU" = "#E69F00",  # Orange
  "MAT" = "#0072B2",  # Deep blue
  "MRB" = "#009E73"   # Teal-green
)

# Data colors for histograms - elegant blue-gray palette
DATA_FILL <- "#4A90A4"     # Sophisticated teal
DENSITY_LINE <- "#C44E52"  # Muted coral-red for contrast

# ============================================================================
# ENHANCED THEME FOR FIGURE 1
# ============================================================================

theme_fig1 <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      # Clean background
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray93", linewidth = 0.25),
      panel.grid.minor = element_blank(),

      # Axis styling - refined
      axis.line = element_line(color = "gray50", linewidth = 0.35),
      axis.ticks = element_line(color = "gray50", linewidth = 0.25),
      axis.text = element_text(color = "gray25", size = base_size - 1),
      axis.title = element_text(color = "gray15", size = base_size),

      # Title styling - bolder panel labels
      plot.title = element_text(face = "bold", size = base_size + 4,
                                color = "gray5", hjust = 0,
                                margin = margin(b = 3)),
      plot.subtitle = element_text(size = base_size - 0.5, color = "gray45",
                                   hjust = 0, margin = margin(b = 10)),

      # Clean margins
      plot.margin = margin(12, 15, 10, 12)
    )
}

# ============================================================================
# PANEL A: MAP OF MO'OREA WITH STUDY SITES
# ============================================================================

cat("Creating Panel A: Study site map...\n")

# Mo'orea island boundary - more accurate heart shape
moorea_coords <- matrix(c(
  -149.75, -17.46,   # NE tip
  -149.77, -17.44,   # N coast
  -149.80, -17.44,   # N coast (Temae)
  -149.84, -17.46,   # N coast
  -149.88, -17.48,   # NW coast
  -149.91, -17.51,   # W coast (Hauru bay area)
  -149.92, -17.54,   # SW coast
  -149.90, -17.58,   # S coast
  -149.86, -17.60,   # SE coast (toward Maatea)
  -149.81, -17.60,   # SE coast
  -149.78, -17.57,   # E coast
  -149.76, -17.53,   # E coast
  -149.75, -17.49,   # NE coast

  -149.75, -17.46    # Close polygon
), ncol = 2, byrow = TRUE)

# Create sf polygon
moorea_sf <- st_polygon(list(moorea_coords)) %>%
  st_sfc(crs = 4326) %>%
  st_sf(name = "Mo'orea")

# Barrier reef outline (surrounding lagoon)
reef_coords <- matrix(c(
  -149.73, -17.43,
  -149.78, -17.41,
  -149.85, -17.43,
  -149.92, -17.46,
  -149.96, -17.51,
  -149.96, -17.56,
  -149.93, -17.61,
  -149.87, -17.64,
  -149.79, -17.64,
  -149.74, -17.61,
  -149.72, -17.54,
  -149.72, -17.47,
  -149.73, -17.43
), ncol = 2, byrow = TRUE)

reef_sf <- st_polygon(list(reef_coords)) %>%
  st_sfc(crs = 4326) %>%
  st_sf(name = "Reef")

# Convert site data to sf
sites_sf <- st_as_sf(site_data, coords = c("long", "lat"), crs = 4326)

# Custom label positions for each site (adjusted to avoid edge clipping)
# Labels positioned outside island with clean white background boxes
label_offsets <- tibble(
  site = c("HAU", "MAT", "MRB"),
  x_off = c(0.045, 0.050, 0.055),  # HAU label moved RIGHT of point (was clipping left)
  y_off_name = c(0.016, 0.028, 0.018),
  y_off_n = c(-0.010, 0.006, -0.010)
)

site_data <- site_data %>%
  left_join(label_offsets, by = "site")

# Create the map
panel_a <- ggplot() +
  # Ocean background - soft blue gradient effect
  annotate("rect", xmin = -150.0, xmax = -149.68, ymin = -17.68, ymax = -17.38,
           fill = "#E8F4F8", color = NA) +

  # Lagoon area (lighter blue between reef and island)
  geom_sf(data = reef_sf, fill = "#D4E8F2", color = "#8FBFD9",
          linewidth = 0.6, linetype = "dashed") +

  # Island (land mass) - rich teal
  geom_sf(data = moorea_sf, fill = "#2D7D7D", color = "#1A5050",
          linewidth = 0.6, alpha = 1.0) +

  # Study sites - larger markers with clear borders and shadows
  geom_sf(data = sites_sf, aes(fill = site),
          shape = 21, size = 7, color = "white", stroke = 2.2) +

  # Site name labels with subtle background
  geom_label(data = site_data,
             aes(x = long + x_off, y = lat + y_off_name, label = site_name),
             size = 3.2, fontface = "bold", color = "gray10",
             fill = alpha("white", 0.85), label.size = 0,
             label.padding = unit(0.15, "lines")) +

  # Sample size labels
  geom_text(data = site_data,
            aes(x = long + x_off, y = lat + y_off_n,
                label = paste0("n = ", n_corals)),
            size = 2.8, color = "gray35", fontface = "plain") +

  scale_fill_manual(values = SITE_COLORS, guide = "none") +

  # Map bounds - slightly wider to prevent clipping
  coord_sf(xlim = c(-149.98, -149.70), ylim = c(-17.66, -17.40), expand = FALSE) +

  # Island name (centered) - white text on teal
  annotate("text", x = -149.84, y = -17.52, label = "Mo'orea",
           size = 5.5, fontface = "bold.italic", color = "white") +

  # Coordinates
  annotate("text", x = -149.84, y = -17.645,
           label = expression("17\u00B030'S, 149\u00B050'W"),
           size = 2.5, color = "gray40") +

  # Scale bar with cleaner styling
  annotate("segment", x = -149.96, xend = -149.87, y = -17.625, yend = -17.625,
           color = "gray30", linewidth = 0.8) +
  annotate("segment", x = -149.96, xend = -149.96, y = -17.63, yend = -17.62,
           color = "gray30", linewidth = 0.6) +
  annotate("segment", x = -149.87, xend = -149.87, y = -17.63, yend = -17.62,
           color = "gray30", linewidth = 0.6) +
  annotate("text", x = -149.915, y = -17.608, label = "5 km",
           size = 2.5, color = "gray30") +

  # North arrow (top right corner) - cleaner design
  annotate("text", x = -149.72, y = -17.415, label = "N",
           size = 4.2, fontface = "bold", color = "gray25") +
  annotate("segment", x = -149.72, xend = -149.72, y = -17.45, yend = -17.42,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
           color = "gray25", linewidth = 0.6) +

  labs(title = "A") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0,
                              margin = margin(b = 5, l = 8)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================================
# PANEL B: COLONY SIZE DISTRIBUTION
# ============================================================================

cat("Creating Panel B: Colony size distribution...\n")

# Calculate statistics
vol_stats <- coral_master %>%
  filter(!is.na(volume), volume > 0) %>%
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

# Calculate density max for positioning annotations
vol_data <- coral_master %>% filter(!is.na(volume), volume > 0)
dens_max_b <- max(density(log10(vol_data$volume))$y)

panel_b <- vol_data %>%
  ggplot(aes(x = volume)) +
  # Histogram - refined bins and styling
  geom_histogram(aes(y = after_stat(density)),
                 bins = 20, fill = DATA_FILL, color = "white",
                 alpha = 0.85, linewidth = 0.3) +
  # Density overlay - elegant curve
  geom_density(color = DENSITY_LINE, linewidth = 1.2, adjust = 1.2) +
  # Rug plot - subtle
  geom_rug(alpha = 0.18, color = DATA_FILL, linewidth = 0.2, sides = "b") +

  # Log scale for x-axis
  scale_x_log10(
    labels = scales::label_comma(),
    breaks = c(100, 1000, 10000, 100000),
    limits = c(18, 150000)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +

  # Elegant bracket annotation
  annotate("segment",
           x = vol_stats$min_vol * 1.1, xend = vol_stats$max_vol * 0.9,
           y = dens_max_b * 1.08, yend = dens_max_b * 1.08,
           color = "gray45", linewidth = 0.5) +
  # Left cap
  annotate("segment",
           x = vol_stats$min_vol * 1.1, xend = vol_stats$min_vol * 1.1,
           y = dens_max_b * 1.04, yend = dens_max_b * 1.08,
           color = "gray45", linewidth = 0.5) +
  # Right cap
  annotate("segment",
           x = vol_stats$max_vol * 0.9, xend = vol_stats$max_vol * 0.9,
           y = dens_max_b * 1.04, yend = dens_max_b * 1.08,
           color = "gray45", linewidth = 0.5) +
  annotate("text",
           x = sqrt(vol_stats$min_vol * vol_stats$max_vol),  # geometric mean position
           y = dens_max_b * 1.20,
           label = ">3 orders of magnitude",
           size = 3.2, color = "gray35", fontface = "italic") +

  labs(
    title = "B",
    subtitle = paste0("n = ", vol_stats$n, " colonies | CV = ", round(vol_stats$cv), "%"),
    x = expression("Colony Volume (cm"^3*")"),
    y = "Density"
  ) +
  theme_fig1() +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.grid.major.x = element_blank()
  )

# ============================================================================
# PANEL C: NEIGHBORHOOD DENSITY DISTRIBUTION
# ============================================================================

cat("Creating Panel C: Neighborhood density distribution...\n")

# Calculate statistics (only for corals with neighborhood data)
neighbor_stats <- coral_master %>%
  filter(!is.na(n_neighbors)) %>%
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

# Get density max for positioning
neighbor_data <- coral_master %>% filter(!is.na(n_neighbors))
dens_max_c <- max(density(neighbor_data$n_neighbors)$y)

panel_c <- neighbor_data %>%
  ggplot(aes(x = n_neighbors)) +
  # Histogram
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 5, fill = DATA_FILL, color = "white",
                 alpha = 0.80, linewidth = 0.25) +
  # Density overlay
  geom_density(color = DENSITY_LINE, linewidth = 1.1, adjust = 1.0) +
  # Rug plot
  geom_rug(alpha = 0.20, color = DATA_FILL, linewidth = 0.25, sides = "b") +

  # Median line - elegant dashed
  geom_vline(xintercept = neighbor_stats$median_n,
             linetype = "dashed", color = "gray45", linewidth = 0.5) +

  # Median annotation - clean placement
  annotate("text", x = neighbor_stats$median_n + 3, y = dens_max_c * 0.85,
           label = paste0("median = ", neighbor_stats$median_n),
           size = 3, color = "gray40", hjust = 0, fontface = "italic") +

  # Range annotation in corner
  annotate("text", x = 70, y = dens_max_c * 1.1,
           label = paste0("range: ", neighbor_stats$min_n, "-", neighbor_stats$max_n),
           size = 2.8, color = "gray45", hjust = 1) +

  scale_x_continuous(breaks = seq(0, 80, by = 20), limits = c(-2, 82)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +

  labs(
    title = "C",
    subtitle = paste0("n = ", neighbor_stats$n, " colonies | CV = ", round(neighbor_stats$cv), "%"),
    x = "Neighboring corals within 5m radius",
    y = "Density"
  ) +
  theme_fig1() +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.grid.major.x = element_blank()
  )

# ============================================================================
# COMBINE PANELS
# ============================================================================

cat("Combining panels...\n")

# Layout: A on left (larger), B and C stacked on right
fig1 <- panel_a + (panel_b / panel_c) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = "Figure 1. Landscape heterogeneity across Mo'orea reef sites",
    subtitle = expression(italic("Pocillopora")*" coral colonies surveyed across three reef habitats | Mo'orea, French Polynesia | Summer 2019"),
    caption = paste0(
      "Study sites span north shore fringing reef (HAU, n=", site_data$n_corals[site_data$site=="HAU"],
      "), east shore back-reef (MAT, n=", site_data$n_corals[site_data$site=="MAT"],
      "), and offshore barrier reef (MRB, n=", site_data$n_corals[site_data$site=="MRB"],
      "). Colony volumes span >3 orders of magnitude (21-42,333 cm\u00B3)."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, color = "gray5",
                                margin = margin(b = 3)),
      plot.subtitle = element_text(size = 10.5, color = "gray40",
                                   margin = margin(b = 12)),
      plot.caption = element_text(size = 8.5, color = "gray50", hjust = 0,
                                  margin = margin(t = 12)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(18, 18, 12, 18)
    )
  )

# ============================================================================
# SAVE FIGURE
# ============================================================================

output_path <- file.path(PATHS$fig_manuscript, "fig1_study_design.png")

ggsave(
  output_path,
  fig1,
  width = 12,
  height = 7.5,
  dpi = 300,
  bg = "white"
)

cat("\n========================================\n")
cat("Figure 1 saved to:\n")
cat(output_path, "\n")
cat("========================================\n\n")

# Summary
cat("Figure specifications:\n")
cat("  - Dimensions: 12 x 7.5 inches\n")
cat("  - Resolution: 300 dpi\n")
cat("  - Panels: A (map), B (size), C (neighbors)\n")
cat("  - Total corals: ", nrow(coral_master), "\n")
cat("  - Corals with neighborhood data: ", neighbor_stats$n, "\n")
