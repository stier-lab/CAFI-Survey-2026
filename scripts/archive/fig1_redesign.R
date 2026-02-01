# ============================================================================
# fig1_redesign.R - Redesigned Figure 1: Emphasize Landscape VARIATION
# ============================================================================
#
# PURPOSE: Create Figure 1 that emphasizes the RANGE and DISTRIBUTION of
#          landscape variables sampled, NOT site-to-site differences
#
# DESIGN CHANGES:
#   - Panel A: Keep satellite map (site locations)
#   - Panel B: Coral Volume distribution (pooled, single color, log scale)
#   - Panel C: Volume vs CAFI scatter (core scaling relationship)
#   - Panel D: Neighborhood density distribution (pooled)
#
# OUTPUT: 14x10 inches, 300 dpi
#
# Author: CAFI Survey Analysis
# Date: 2026-01-25
# ============================================================================

cat("\n========================================\n")
cat("FIGURE 1 REDESIGN: Landscape Variation\n")
cat("========================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

# Set PROJ_LIB for spatial operations
Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sf/proj")

# Load setup and data
source(here::here("scripts/00_setup.R"))
source(here::here("scripts/01_load_data.R"))

# Load additional packages for mapping
library(sf)
library(terra)
library(ggtext)

# Create output directory
MANUSCRIPT_DIR <- file.path(PATHS$figures, "manuscript")
dir.create(MANUSCRIPT_DIR, showWarnings = FALSE, recursive = TRUE)

# Define teal gradient color palette (cohesive single-color scheme)
TEAL_MAIN <- "#0072B2"      # Main teal-blue
TEAL_LIGHT <- "#56B4E9"     # Light sky blue
TEAL_DARK <- "#005A8C"      # Darker teal

cat("Data loaded successfully\n")
cat("  Corals:", nrow(coral_master), "\n")
cat("  CAFI records:", nrow(cafi_clean), "\n\n")

# ============================================================================
# CALCULATE SUMMARY STATISTICS FOR ANNOTATIONS
# ============================================================================

cat("Calculating summary statistics...\n")

# Volume statistics
vol_stats <- list(
  min = min(coral_master$volume, na.rm = TRUE),
  max = max(coral_master$volume, na.rm = TRUE),
  median = median(coral_master$volume, na.rm = TRUE),
  cv = sd(coral_master$volume, na.rm = TRUE) / mean(coral_master$volume, na.rm = TRUE) * 100,
  range_orders = log10(max(coral_master$volume, na.rm = TRUE)) - log10(min(coral_master$volume, na.rm = TRUE))
)

# CAFI statistics
cafi_stats <- list(
  min = min(coral_master$total_cafi, na.rm = TRUE),
  max = max(coral_master$total_cafi, na.rm = TRUE),
  median = median(coral_master$total_cafi, na.rm = TRUE),
  cv = sd(coral_master$total_cafi, na.rm = TRUE) / mean(coral_master$total_cafi, na.rm = TRUE) * 100
)

# Neighborhood statistics (subset with data)
neighbor_data <- coral_master %>% filter(!is.na(n_neighbors))
neigh_stats <- list(
  n = nrow(neighbor_data),
  min = min(neighbor_data$n_neighbors, na.rm = TRUE),
  max = max(neighbor_data$n_neighbors, na.rm = TRUE),
  median = median(neighbor_data$n_neighbors, na.rm = TRUE),
  cv = sd(neighbor_data$n_neighbors, na.rm = TRUE) / mean(neighbor_data$n_neighbors, na.rm = TRUE) * 100
)

cat("  Volume range:", round(vol_stats$min), "to", format(round(vol_stats$max), big.mark=","), "cm3\n")
cat("  Volume spans", round(vol_stats$range_orders, 1), "orders of magnitude\n")
cat("  CAFI range:", cafi_stats$min, "to", cafi_stats$max, "per coral\n")
cat("  Neighbors range:", neigh_stats$min, "to", neigh_stats$max, "(n =", neigh_stats$n, "corals)\n\n")

# ============================================================================
# PANEL A: SATELLITE MAP WITH SITE MARKERS
# ============================================================================

cat("Creating Panel A: Satellite map...\n")

# Try to load satellite image
satellite_path <- here::here("output/figures/moorea_satellite.tif")

if (file.exists(satellite_path)) {

  # Load satellite raster
  sat_raster <- terra::rast(satellite_path)

  # Convert to data frame for ggplot
  sat_df <- as.data.frame(sat_raster, xy = TRUE)
  names(sat_df) <- c("x", "y", "r", "g", "b")

  # Normalize RGB values
  sat_df <- sat_df %>%
    mutate(
      r = pmin(pmax(r / 255, 0), 1),
      g = pmin(pmax(g / 255, 0), 1),
      b = pmin(pmax(b / 255, 0), 1)
    )

  # Site coordinates
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

  # Calculate extent for scale bar
  x_range <- range(sat_df$x, na.rm = TRUE)
  y_range <- range(sat_df$y, na.rm = TRUE)

  # Create map panel
  panel_a <- ggplot() +
    geom_raster(data = sat_df, aes(x = x, y = y, fill = rgb(r, g, b))) +
    scale_fill_identity() +
    # Site markers - white circles with colored fill
    geom_point(data = site_coords,
               aes(x = long, y = lat),
               shape = 21, size = 5, stroke = 1.5,
               fill = SITE_COLORS[site_coords$site], color = "white") +
    # Site labels
    geom_label(data = site_coords,
               aes(x = long, y = lat,
                   label = paste0(site_name, "\n(n=", n, ")")),
               size = 2.8, fontface = "bold",
               fill = alpha("white", 0.85),
               label.padding = unit(0.15, "lines"),
               label.size = 0.3,
               nudge_y = c(0.025, 0.03, -0.035)) +
    # Scale bar
    annotate("segment",
             x = x_range[1] + 0.02, xend = x_range[1] + 0.02 + 0.054,  # ~6km
             y = y_range[1] + 0.02, yend = y_range[1] + 0.02,
             color = "white", linewidth = 2) +
    annotate("text",
             x = x_range[1] + 0.02 + 0.027, y = y_range[1] + 0.04,
             label = "6 km", color = "white", size = 3, fontface = "bold") +
    # North arrow
    annotate("text", x = x_range[2] - 0.03, y = y_range[2] - 0.02,
             label = "N", color = "white", size = 5, fontface = "bold") +
    annotate("segment",
             x = x_range[2] - 0.03, xend = x_range[2] - 0.03,
             y = y_range[2] - 0.04, yend = y_range[2] - 0.025,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             color = "white", linewidth = 1) +
    coord_fixed(ratio = 1) +
    labs(title = "A") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0, margin = margin(b = 5)),
      plot.margin = margin(5, 5, 5, 5)
    )

  cat("  Satellite map loaded successfully\n")

} else {
  # Fallback: schematic map if satellite not available
  cat("  Satellite image not found, creating schematic...\n")

  site_data <- tibble(
    site = c("HAU", "MAT", "MRB"),
    site_name = c("Hauru\n(North Shore)", "Maatea\n(East Shore)", "Barrier Reef\n(Offshore)"),
    x = c(1.2, 2.8, 2.8),
    y = c(1.6, 0.9, 2.2),
    n = c(
      sum(coral_master$site == "HAU"),
      sum(coral_master$site == "MRB"),
      sum(coral_master$site == "MAT")
    )
  )

  panel_a <- ggplot() +
    annotate("path",
             x = 2 + 1.6 * cos(seq(0, 2*pi, length.out = 100)),
             y = 1.5 + 1.1 * sin(seq(0, 2*pi, length.out = 100)),
             color = "#87CEEB", linewidth = 0.8, linetype = "dashed") +
    annotate("path",
             x = 2 + 1.0 * cos(seq(0, 2*pi, length.out = 100)),
             y = 1.5 + 0.7 * sin(seq(0, 2*pi, length.out = 100)),
             color = "gray50", linewidth = 1.2) +
    geom_point(data = site_data, aes(x = x, y = y),
               color = "black", fill = SITE_COLORS[site_data$site],
               shape = 21, size = 6, stroke = 1.2) +
    geom_text(data = site_data, aes(x = x, y = y - 0.32, label = site_name),
              size = 2.8, fontface = "bold", lineheight = 0.85) +
    annotate("text", x = 2, y = 1.5, label = "Mo'orea", size = 4.5, fontface = "italic") +
    coord_fixed(ratio = 1, xlim = c(0, 4), ylim = c(-0.1, 2.8)) +
    labs(title = "A") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
}

# ============================================================================
# PANEL B: CORAL VOLUME DISTRIBUTION (POOLED, LOG SCALE)
# ============================================================================

cat("Creating Panel B: Volume distribution...\n")

# Create data for density plot with rug
vol_data <- coral_master %>%
  filter(!is.na(volume), volume > 0) %>%
  mutate(log_vol = log10(volume))

panel_b <- ggplot(vol_data, aes(x = volume)) +
  # Density curve with fill
  geom_density(fill = TEAL_MAIN, color = TEAL_DARK, alpha = 0.6, linewidth = 1) +
  # Individual data points as rug
  geom_rug(sides = "b", alpha = 0.4, color = TEAL_DARK, linewidth = 0.3) +
  # Vertical lines for min, median, max
  geom_vline(xintercept = vol_stats$median, linetype = "dashed",
             color = "gray30", linewidth = 0.8) +
  # Annotations
  annotate("text", x = vol_stats$median * 1.5, y = Inf,
           label = paste0("Median: ", format(round(vol_stats$median), big.mark = ",")),
           vjust = 2, hjust = 0, size = 3, color = "gray30") +
  # Range annotation box
  annotate("label", x = 100, y = Inf,
           label = paste0("Range: ", round(vol_stats$min), " - ",
                          format(round(vol_stats$max), big.mark = ","), " cm\u00B3\n",
                          "~3 orders of magnitude"),
           vjust = 1.5, hjust = 0, size = 3, fontface = "bold",
           fill = alpha("white", 0.9), label.size = 0) +
  scale_x_log10(
    labels = scales::label_comma(),
    breaks = c(10, 100, 1000, 10000, 100000),
    limits = c(10, 100000)
  ) +
  labs(
    x = expression("Coral Volume (cm"^3*")"),
    y = "Density",
    title = "B",
    subtitle = paste0("n = ", nrow(vol_data), " corals | CV = ", round(vol_stats$cv), "%")
  ) +
  theme_publication() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    panel.grid.minor.x = element_blank()
  )

# ============================================================================
# PANEL C: VOLUME vs CAFI SCATTER (CORE SCALING RELATIONSHIP)
# ============================================================================

cat("Creating Panel C: Volume vs CAFI scatter...\n")

# Fit model for trend line
scaling_model <- MASS::glm.nb(total_cafi ~ log10(volume), data = coral_master)
beta_est <- coef(scaling_model)["log10(volume)"]
beta_ci <- confint(scaling_model, "log10(volume)")

# Create prediction data for smooth line
pred_data <- tibble(
  volume = 10^seq(log10(20), log10(50000), length.out = 100)
) %>%
  mutate(
    predicted = predict(scaling_model, newdata = ., type = "response")
  )

panel_c <- ggplot(coral_master, aes(x = volume, y = total_cafi)) +
  # Data points with transparency
  geom_point(alpha = 0.6, size = 2.5, color = TEAL_MAIN) +
  # Fitted line
  geom_line(data = pred_data, aes(x = volume, y = predicted),
            color = TEAL_DARK, linewidth = 1.2) +
  # 1:1 reference line (Field of Dreams expectation)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "gray50", linewidth = 0.6) +
  # Scaling annotation
  annotate("label", x = 50, y = 200,
           label = paste0("\u03B2 = ", round(beta_est, 2),
                          " [", round(beta_ci[1], 2), ", ", round(beta_ci[2], 2), "]"),
           size = 3.5, fontface = "bold",
           fill = alpha("white", 0.9), label.size = 0, hjust = 0) +
  annotate("text", x = 300, y = 10,
           label = "1:1 line", size = 2.8, color = "gray50", fontface = "italic") +
  scale_x_log10(
    labels = scales::label_comma(),
    breaks = c(100, 1000, 10000)
  ) +
  scale_y_log10(
    labels = scales::label_comma(),
    breaks = c(1, 10, 100)
  ) +
  labs(
    x = expression("Coral Volume (cm"^3*")"),
    y = "Total CAFI Abundance",
    title = "C",
    subtitle = "Larger corals host more fauna"
  ) +
  theme_publication() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

# ============================================================================
# PANEL D: NEIGHBORHOOD DENSITY DISTRIBUTION (POOLED)
# ============================================================================

cat("Creating Panel D: Neighborhood density distribution...\n")

# Filter to corals with neighborhood data
neigh_data <- coral_master %>%
  filter(!is.na(n_neighbors))

panel_d <- ggplot(neigh_data, aes(x = n_neighbors)) +
  # Histogram with density overlay
  geom_histogram(aes(y = after_stat(density)),
                 bins = 20, fill = TEAL_MAIN, color = "white", alpha = 0.7) +
  geom_density(color = TEAL_DARK, linewidth = 1.2) +
  # Individual data points as rug
  geom_rug(sides = "b", alpha = 0.5, color = TEAL_DARK, linewidth = 0.3) +
  # Median line
  geom_vline(xintercept = neigh_stats$median, linetype = "dashed",
             color = "gray30", linewidth = 0.8) +
  # Range annotation
  annotate("label", x = 60, y = Inf,
           label = paste0("Range: ", neigh_stats$min, " - ", neigh_stats$max, " neighbors\n",
                          "Median: ", round(neigh_stats$median)),
           vjust = 1.5, hjust = 0, size = 3, fontface = "bold",
           fill = alpha("white", 0.9), label.size = 0) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  labs(
    x = "Number of Neighbors (within 5m)",
    y = "Density",
    title = "D",
    subtitle = paste0("n = ", neigh_stats$n, " corals | CV = ", round(neigh_stats$cv), "%")
  ) +
  theme_publication() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

# ============================================================================
# COMBINE PANELS INTO FINAL FIGURE
# ============================================================================

cat("Combining panels...\n")

# Layout: A takes left half, B-C-D stack on right
fig1_redesign <- (panel_a | (panel_b / panel_c / panel_d)) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = "Figure 1: Study Design and Landscape Variation",
    subtitle = expression(italic("Pocillopora")*" coral-associated fauna survey, Mo'orea, French Polynesia (Summer 2019)"),
    caption = paste0(
      "A: Satellite imagery (Esri) showing three reef sites. ",
      "B: Coral colony sizes span 3 orders of magnitude. ",
      "C: Total CAFI abundance scales with coral volume (\u03B2 = ", round(beta_est, 2), "). ",
      "D: Neighborhood density ranges from isolated to crowded."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray30"),
      plot.caption = element_text(size = 9, color = "gray50", hjust = 0, lineheight = 1.2),
      plot.margin = margin(10, 10, 10, 10)
    )
  )

# ============================================================================
# SAVE FIGURE
# ============================================================================

output_path <- file.path(MANUSCRIPT_DIR, "fig1_study_design.png")

ggsave(output_path, fig1_redesign,
       width = 14, height = 10, dpi = 300, bg = "white")

cat("\n========================================\n")
cat("FIGURE 1 REDESIGN COMPLETE\n")
cat("========================================\n")
cat("Saved to:", output_path, "\n\n")

cat("Key messages emphasized:\n")
cat("  - Corals span 3 orders of magnitude in size (20-42,000 cm3)\n")
cat("  - CAFI abundance ranges from", cafi_stats$min, "to", cafi_stats$max, "per coral\n")
cat("  - Neighborhood density varies from isolated (0) to crowded (", neigh_stats$max, ")\n")
cat("  - Single teal color scheme emphasizes variation, not site differences\n\n")
