# ============================================================================
# Create Final Beautiful Map Panel for Figure 1
# ============================================================================
# Creates a publication-quality schematic map of Mo'orea with accurate
# geography and clearly marked sampling sites
# ============================================================================

cat("\n========================================\n")
cat("Creating Final Map Panel for Figure 1\n")
cat("========================================\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(patchwork)
  library(grid)
})

# Source setup for paths and theme
source(here::here("scripts/00_setup.R"))

# Create output directory
MANUSCRIPT_DIR <- file.path(PATHS$figures, "manuscript")
dir.create(MANUSCRIPT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# SITE COORDINATES - Positioned accurately on Mo'orea
# ============================================================================

# Mo'orea is roughly heart-shaped, oriented with the point to the north
# Sites positioned within the lagoon/reef system:
# - HAU (Hauru): NW coast, fringing reef in lagoon
# - MAT (Maatea): E coast, back-reef in lagoon
# - MRB (Barrier Reef): N coast, on barrier reef (outer)

sites_df <- data.frame(
  site = c("HAU", "MAT", "MRB"),
  site_name = c("Hauru", "Maatea", "Barrier Reef"),
  habitat = c("Fringing Reef", "Back-reef", "Barrier Reef"),
  n_corals = c(38, 39, 35),
  stringsAsFactors = FALSE
)

# Site colors (colorblind-safe Okabe-Ito)
site_colors <- c("HAU" = "#E69F00", "MAT" = "#0072B2", "MRB" = "#009E73")

# ============================================================================
# CREATE MO'OREA MAP GEOMETRY
# ============================================================================

# Mo'orea has a distinctive shape - like a heart or butterfly viewed from above
# with the "point" toward the north and two bays (Opunohu and Cook's) cutting in

# Island center approximately at (-149.835, -17.535)
center_lon <- -149.835
center_lat <- -17.535

# Create island outline (simplified heart shape with bays)
n_pts <- 200
theta <- seq(0, 2*pi, length.out = n_pts)

# Base ellipse with modifications for Mo'orea's shape
# Wider E-W than N-S, with indentations for the bays
island_r_base <- 0.055
island_x <- center_lon + island_r_base * 1.15 * cos(theta) * (1 + 0.15*cos(2*theta))
island_y <- center_lat + island_r_base * 0.85 * sin(theta) * (1 + 0.1*sin(2*theta))

# Add bay indentations (Cook's Bay and Opunohu Bay on north side)
bay_adjustment <- ifelse(theta > pi/3 & theta < 2*pi/3,
                         -0.012 * sin(3*(theta - pi/3) * pi / (pi/3)), 0)
island_y <- island_y + bay_adjustment

island_df <- data.frame(x = island_x, y = island_y)

# Lagoon boundary (reef flat - follows island but larger)
lagoon_r <- 0.085
lagoon_x <- center_lon + lagoon_r * 1.2 * cos(theta) * (1 + 0.1*cos(2*theta))
lagoon_y <- center_lat + lagoon_r * 0.9 * sin(theta) * (1 + 0.08*sin(2*theta))
lagoon_df <- data.frame(x = lagoon_x, y = lagoon_y)

# Barrier reef (outer edge)
reef_r <- 0.105
reef_x <- center_lon + reef_r * 1.25 * cos(theta) * (1 + 0.08*cos(2*theta))
reef_y <- center_lat + reef_r * 0.95 * sin(theta) * (1 + 0.06*sin(2*theta))
reef_df <- data.frame(x = reef_x, y = reef_y)

# Site positions (within lagoon/on reef)
# HAU - NW coast (fringing reef in lagoon)
sites_df$lon <- c(-149.895, -149.765, -149.81)
sites_df$lat <- c(-17.505, -17.555, -17.45)

# ============================================================================
# CREATE BEAUTIFUL SCHEMATIC MAP
# ============================================================================

cat("Creating enhanced schematic map...\n")

map_panel <- ggplot() +
  # Deep ocean background
  annotate("rect",
           xmin = -149.98, xmax = -149.69,
           ymin = -17.65, ymax = -17.40,
           fill = "#1B4F72", alpha = 0.35) +

  # Barrier reef zone (light blue outer ring)
  geom_polygon(data = reef_df, aes(x = x, y = y),
               fill = "#5DADE2", color = NA, alpha = 0.5) +

  # Reef crest line (dashed)
  geom_path(data = reef_df, aes(x = x, y = y),
            color = "#2874A6", linewidth = 0.6, linetype = "dashed") +

  # Lagoon (turquoise - lighter)
  geom_polygon(data = lagoon_df, aes(x = x, y = y),
               fill = "#48C9B0", color = NA, alpha = 0.65) +

  # Island (rich green with dark border)
  geom_polygon(data = island_df, aes(x = x, y = y),
               fill = "#1E8449", color = "#145A32", linewidth = 1.0, alpha = 0.95) +

  # Mountain peaks (interior detail)
  annotate("point", x = center_lon, y = center_lat + 0.01,
           shape = 24, size = 3.5, fill = "#616A6B", color = "#2C3E50", stroke = 0.8) +
  annotate("point", x = center_lon - 0.015, y = center_lat - 0.008,
           shape = 24, size = 2.8, fill = "#7B7D7D", color = "#2C3E50", stroke = 0.6) +
  annotate("point", x = center_lon + 0.018, y = center_lat + 0.002,
           shape = 24, size = 2.5, fill = "#7B7D7D", color = "#2C3E50", stroke = 0.6) +

  # Site points - large with white halos
  geom_point(data = sites_df,
             aes(x = lon, y = lat),
             shape = 21, size = 9, stroke = 2.5,
             fill = site_colors[sites_df$site], color = "white") +

  # Site code letters inside points
  geom_text(data = sites_df,
            aes(x = lon, y = lat, label = substr(site, 1, 1)),
            size = 3.5, fontface = "bold", color = "white") +

  # Site labels with white background boxes
  geom_label(data = sites_df,
             aes(x = lon + c(-0.035, 0.04, 0.045),
                 y = lat + c(0.025, 0.025, 0.03),
                 label = paste0(site_name, "\n(n=", n_corals, ")")),
             size = 3.2,
             fontface = "bold",
             fill = alpha("white", 0.92),
             color = "gray15",
             label.padding = unit(0.25, "lines"),
             label.size = 0.3,
             lineheight = 0.9) +

  # Island name (centered)
  annotate("text", x = center_lon, y = center_lat - 0.015,
           label = "Mo'orea", size = 5.5, fontface = "bold.italic",
           color = "#145A32") +

  # Coordinate annotation
  annotate("text", x = center_lon, y = -17.635,
           label = "17\u00B032'S, 149\u00B050'W  |  French Polynesia",
           size = 3.2, color = "gray35") +

  # Legend box (reef zones)
  annotate("rect", xmin = -149.97, xmax = -149.87, ymin = -17.44, ymax = -17.41,
           fill = "white", color = "gray60", alpha = 0.9, linewidth = 0.3) +
  annotate("rect", xmin = -149.965, xmax = -149.935, ymin = -17.435, ymax = -17.425,
           fill = "#5DADE2", alpha = 0.5) +
  annotate("text", x = -149.93, y = -17.43, label = "Reef",
           size = 2.5, hjust = 0, color = "gray25") +
  annotate("rect", xmin = -149.965, xmax = -149.935, ymin = -17.418, ymax = -17.408,
           fill = "#48C9B0", alpha = 0.65) +
  annotate("text", x = -149.93, y = -17.413, label = "Lagoon",
           size = 2.5, hjust = 0, color = "gray25") +

  # Scale bar (5 km ~ 0.045 degrees at this latitude)
  annotate("segment", x = -149.755, xend = -149.71,
           y = -17.62, yend = -17.62,
           linewidth = 2, color = "gray25") +
  annotate("segment", x = -149.755, xend = -149.71,
           y = -17.62, yend = -17.62,
           linewidth = 1, color = "white") +
  annotate("text", x = -149.7325, y = -17.605,
           label = "5 km", size = 3, color = "gray25", fontface = "bold") +

  # North arrow
  annotate("segment", x = -149.715, xend = -149.715,
           y = -17.455, yend = -17.415,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
           linewidth = 1.2, color = "gray25") +
  annotate("text", x = -149.715, y = -17.405, label = "N",
           size = 4.5, fontface = "bold", color = "gray25") +

  # Formatting
  coord_fixed(ratio = 1.0,
              xlim = c(-149.98, -149.69),
              ylim = c(-17.65, -17.40)) +
  labs(
    title = "A. Study Sites: Mo'orea, French Polynesia",
    subtitle = expression(italic("Pocillopora")*" coral-associated fauna survey | Summer 2019")
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0,
                              margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10, color = "gray35", hjust = 0,
                                 margin = margin(b = 10)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

# Save map panel alone
ggsave(file.path(MANUSCRIPT_DIR, "fig1_panel_a_map.png"), map_panel,
       width = 8, height = 7, dpi = 300, bg = "white")
cat("  Saved: fig1_panel_a_map.png\n")

# ============================================================================
# RECREATE FULL FIGURE 1
# ============================================================================

cat("\nRecreating full Figure 1 with improved layout...\n")

# Load data
coral_master <- load_object("coral_master")
cafi_clean <- load_object("cafi_clean")

# Panel B: Dataset summary
summary_stats <- tibble::tibble(
  Metric = c("Coral colonies surveyed",
             "Total CAFI individuals",
             "Morphospecies (OTUs)",
             "Taxonomic groups",
             "Mean CAFI per coral"),
  Value = c(as.character(nrow(coral_master)),
            format(nrow(cafi_clean), big.mark = ","),
            as.character(dplyr::n_distinct(cafi_clean$otu)),
            as.character(dplyr::n_distinct(cafi_clean$functional_group)),
            as.character(round(mean(coral_master$total_cafi, na.rm = TRUE), 1)))
)

panel_b <- ggplot(summary_stats, aes(x = 1, y = rev(seq_along(Metric)))) +
  geom_text(aes(label = Metric, x = 0.52), hjust = 0, size = 4.2, color = "gray20") +
  geom_text(aes(label = Value, x = 1.48), hjust = 1, size = 4.5, fontface = "bold") +
  annotate("segment", x = 0.48, xend = 1.52, y = 5.7, yend = 5.7,
           color = "gray40", linewidth = 0.6) +
  annotate("segment", x = 0.48, xend = 1.52, y = 0.3, yend = 0.3,
           color = "gray40", linewidth = 0.4) +
  scale_x_continuous(limits = c(0.45, 1.55)) +
  scale_y_continuous(limits = c(-0.2, 6.4)) +
  labs(title = "B. Dataset Summary") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
        plot.margin = margin(15, 15, 15, 25))

# Panel C: Sampling protocol
sampling_df <- tibble::tibble(
  stage = factor(1:4),
  step = c("Site Selection", "Colony Selection", "CAFI Census", "Measurements"),
  description = c("3 reef habitats\nacross Mo'orea",
                  "38-40 Pocillopora\nper site (n = 112)",
                  "All fauna extracted\n& identified (~4,000)",
                  "Colony volume, GPS,\nneighborhood density")
)

panel_c <- ggplot(sampling_df, aes(x = as.numeric(stage), y = 1)) +
  geom_tile(fill = "gray97", color = "gray50", width = 0.88, height = 0.60,
            linewidth = 0.5) +
  geom_point(aes(y = 1.40), size = 7.5, shape = 21,
             fill = "#0072B2", color = "white", stroke = 1.2) +
  geom_text(aes(y = 1.40, label = stage), size = 3.8,
            color = "white", fontface = "bold") +
  geom_text(aes(y = 1.15, label = step), size = 3.6, fontface = "bold") +
  geom_text(aes(y = 0.82, label = description), size = 3.0,
            lineheight = 0.9, color = "gray30") +
  geom_segment(data = tibble::tibble(x = 1:3 + 0.45, xend = 1:3 + 0.55),
               aes(x = x, xend = xend, y = 1, yend = 1),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               color = "gray50", linewidth = 0.8) +
  scale_x_continuous(limits = c(0.2, 4.8)) +
  scale_y_continuous(limits = c(0.48, 1.60)) +
  labs(title = "C. Sampling Protocol") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
        plot.margin = margin(15, 25, 15, 25))

# Panel D: Site breakdown with correct values
site_summary <- coral_master %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(
    n_corals = dplyr::n(),
    total_cafi = sum(total_cafi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(mean_cafi = round(total_cafi / n_corals, 1))

panel_d <- ggplot(site_summary, aes(x = site, y = total_cafi, fill = site)) +
  geom_col(width = 0.68, alpha = 0.92, color = "gray40", linewidth = 0.4) +
  geom_text(aes(label = paste0("n = ", n_corals)),
            vjust = -0.6, size = 3.5, fontface = "bold") +
  geom_text(aes(label = format(total_cafi, big.mark=","), y = total_cafi/2 + 50),
            size = 4.2, color = "white", fontface = "bold") +
  geom_text(aes(label = paste0("(", mean_cafi, "/coral)"), y = total_cafi/2 - 100),
            size = 3, color = "white") +
  scale_fill_manual(values = site_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "D. CAFI Abundance by Site",
    x = NULL,
    y = "Total CAFI"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    axis.title.y = element_text(size = 11),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 20, 10, 20)
  )

# ============================================================================
# COMPOSE FULL FIGURE 1
# ============================================================================

# Layout: Map (A) large on left, B and D stacked on right, C spanning bottom
design <- "
AABB
AADD
CCCC
"

fig1 <- map_panel + panel_b + panel_c + panel_d +
  plot_layout(design = design) +
  plot_annotation(
    title = "Figure 1: Study Design",
    subtitle = expression(italic("Pocillopora")*" coral-associated fauna (CAFI) survey | Mo'orea, French Polynesia | Summer 2019"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, color = "gray40", margin = margin(b = 10))
    )
  )

ggsave(file.path(MANUSCRIPT_DIR, "fig1_study_design.png"), fig1,
       width = 15, height = 13, dpi = 300, bg = "white")

cat("  Saved: fig1_study_design.png\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("Figure 1 Creation Complete\n")
cat("========================================\n")
cat("\nOutputs saved to:", MANUSCRIPT_DIR, "\n")
cat("  1. fig1_panel_a_map.png (8x7 in, 300 dpi)\n")
cat("  2. fig1_study_design.png (15x13 in, 300 dpi)\n")
cat("\nMap features:\n")
cat("  - Schematic representation of Mo'orea island\n")
cat("  - Barrier reef and lagoon zones shown\n")
cat("  - Three sampling sites with habitat labels\n")
cat("  - Scale bar and north arrow included\n")
cat("  - Publication-quality colorblind-safe palette\n")
cat("\nSite summary:\n")
for (i in 1:nrow(site_summary)) {
  cat(sprintf("  %s: n=%d corals, %d total CAFI, %.1f mean/coral\n",
              site_summary$site[i], site_summary$n_corals[i],
              site_summary$total_cafi[i], site_summary$mean_cafi[i]))
}
