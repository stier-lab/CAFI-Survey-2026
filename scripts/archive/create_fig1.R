# ============================================================================
# create_fig1.R - Create Publication-Quality Figure 1 for CAFI Manuscript
# ============================================================================
#
# PURPOSE: Generate Figure 1 with 4 panels:
#   A - Map of Mo'orea showing three sampling sites
#   B - Coral volume distribution by site
#   C - Neighbor density by site
#   D - Total CAFI abundance by site
#
# USAGE: source("scripts/create_fig1.R")
#
# OUTPUT: output/figures/manuscript/fig1_study_design.png
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    Creating Figure 1: Study Design\n")
cat("============================================================\n\n")

# ============================================================================
# LOAD DEPENDENCIES
# ============================================================================

source("scripts/00_setup.R")
source("scripts/01_load_data.R")

# Additional packages for mapping
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
if (!requireNamespace("rnaturalearth", quietly = TRUE)) install.packages("rnaturalearth")
if (!requireNamespace("rnaturalearthdata", quietly = TRUE)) install.packages("rnaturalearthdata")

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

cat("[OK] Dependencies loaded\n\n")

# ============================================================================
# THEME AND COLOR SETUP
# ============================================================================

# Enhanced publication theme for manuscript figures
theme_manuscript <- function(base_size = 12) {
  theme_publication(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "gray30"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      plot.margin = margin(15, 15, 15, 15)
    )
}

# Site colors - using SITE_COLORS from 00_setup.R
# HAU: #E69F00 (Orange), MAT: #56B4E9 (Blue), MRB: #009E73 (Green)
site_colors <- SITE_COLORS

# Site labels for display
site_labels <- c(
  "HAU" = "Hauru (NW)",
  "MAT" = "Maatea (E)",
  "MRB" = "Barrier Reef (N)"
)

# Ensure site order is consistent
coral_master <- coral_master %>%
  mutate(site = factor(site, levels = c("HAU", "MAT", "MRB")))

cat("[OK] Theme and colors configured\n\n")

# ============================================================================
# PANEL A: MAP OF MO'OREA WITH SAMPLING SITES
# ============================================================================

cat("Creating Panel A: Map of Mo'orea...\n")

# Define Mo'orea approximate coordinates
moorea_center <- c(lon = -149.83, lat = -17.53)

# Site coordinates (approximate - reef sites around the island)
# HAU (Hauru): NW fringing reef
# MAT (Maatea): E back-reef
# MRB (Barrier Reef): N barrier reef
sites_sf <- tibble(
  site = c("HAU", "MAT", "MRB"),
  site_name = c("Hauru", "Maatea", "Barrier Reef"),
  lon = c(-149.90, -149.75, -149.82),
  lat = c(-17.49, -17.54, -17.465),
  n_corals = c(
    sum(coral_master$site == "HAU"),
    sum(coral_master$site == "MAT"),
    sum(coral_master$site == "MRB")
  )
) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  mutate(site = factor(site, levels = c("HAU", "MAT", "MRB")))

# Create a more accurate Mo'orea outline polygon
# Mo'orea is volcanic with two prominent bays (Cook's and Opunohu) cutting into the north
moorea_coords <- matrix(c(
  # Starting from NW coast, going clockwise
  -149.91, -17.47,  # NW coast (Hauru area)
  -149.88, -17.47,  # N coast west of Opunohu
  -149.865, -17.475, # Opunohu Bay entrance west
  -149.855, -17.50,  # Inner Opunohu Bay
  -149.845, -17.475, # Opunohu Bay entrance east
  -149.83, -17.47,  # Between the bays (Pao Pao ridge)
  -149.81, -17.475, # Cook's Bay entrance west
  -149.80, -17.505, # Inner Cook's Bay
  -149.79, -17.475, # Cook's Bay entrance east
  -149.76, -17.48,  # NE coast
  -149.74, -17.50,  # E coast north (Maatea area)
  -149.74, -17.54,  # E coast
  -149.76, -17.58,  # SE coast
  -149.80, -17.60,  # S coast
  -149.85, -17.59,  # SW coast
  -149.90, -17.56,  # W coast south
  -149.91, -17.52,  # W coast (Hauru area)
  -149.91, -17.47   # Close polygon at NW
), ncol = 2, byrow = TRUE)

moorea_polygon <- st_polygon(list(moorea_coords)) %>%
  st_sfc(crs = 4326) %>%
  st_sf(name = "Mo'orea")

# Create the map
panel_a <- ggplot() +
  # Mo'orea island polygon
  geom_sf(data = moorea_polygon, fill = "gray85", color = "gray40", linewidth = 0.5) +
  # Site points
  geom_sf(data = sites_sf, aes(color = site), size = 5, shape = 19) +
  # Site labels with slight offset (manually positioned for clarity)
  geom_sf_text(data = sites_sf,
               aes(label = site_name, color = site),
               nudge_y = c(-0.025, 0.025, 0.025),
               nudge_x = c(-0.01, 0.02, 0.01),
               size = 3.5, fontface = "bold",
               show.legend = FALSE) +
  # Styling
  scale_color_manual(values = site_colors, guide = "none") +
  # Add scale bar approximation (text annotation)
  annotate("segment", x = -149.88, xend = -149.83, y = -17.61, yend = -17.61,
           linewidth = 0.8, color = "black") +
  annotate("text", x = -149.855, y = -17.625, label = "5 km", size = 3) +
  # Add north arrow
  annotate("segment", x = -149.735, xend = -149.735, y = -17.58, yend = -17.55,
           linewidth = 0.8, color = "black",
           arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  annotate("text", x = -149.735, y = -17.595, label = "N", size = 3, fontface = "bold") +
  # Coordinate labels
  coord_sf(xlim = c(-149.93, -149.72), ylim = c(-17.63, -17.44)) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_manuscript(base_size = 11) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    plot.margin = margin(10, 10, 10, 10)
  )

cat("   [OK] Panel A created\n")

# ============================================================================
# PANEL B: CORAL VOLUME BY SITE
# ============================================================================

cat("Creating Panel B: Coral volume by site...\n")

# Use volume_field if available, otherwise volume
volume_col <- if ("volume_field" %in% names(coral_master)) "volume_field" else "volume"

panel_b <- ggplot(coral_master, aes(x = site, y = .data[[volume_col]], fill = site)) +
  geom_violin(alpha = 0.6, color = NA, width = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "gray30") +
  geom_jitter(aes(color = site), width = 0.1, alpha = 0.4, size = 1.2) +
  scale_y_log10(
    labels = scales::label_comma(),
    breaks = c(100, 1000, 10000, 100000)
  ) +
  scale_fill_manual(values = site_colors, guide = "none") +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(
    x = "Site",
    y = expression(bold("Coral Volume"~(cm^3)))
  ) +
  theme_manuscript(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

cat("   [OK] Panel B created\n")

# ============================================================================
# PANEL C: NEIGHBOR DENSITY BY SITE
# ============================================================================

cat("Creating Panel C: Neighbor density by site...\n")

# Filter to corals with neighbor data
neighbor_data <- coral_master %>%
  filter(!is.na(n_neighbors))

n_with_neighbors <- nrow(neighbor_data)
cat("   Corals with neighbor data:", n_with_neighbors, "/", nrow(coral_master), "\n")

panel_c <- ggplot(neighbor_data, aes(x = site, y = n_neighbors, fill = site)) +
  geom_violin(alpha = 0.6, color = NA, width = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "gray30") +
  geom_jitter(aes(color = site), width = 0.1, alpha = 0.4, size = 1.2) +
  scale_fill_manual(values = site_colors, guide = "none") +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(
    x = "Site",
    y = expression(bold("Neighbors within 5m"))
  ) +
  annotate("text", x = 2, y = max(neighbor_data$n_neighbors, na.rm = TRUE) * 0.95,
           label = paste0("n = ", n_with_neighbors, " corals"),
           size = 3, color = "gray40", fontface = "italic") +
  theme_manuscript(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

cat("   [OK] Panel C created\n")

# ============================================================================
# PANEL D: TOTAL CAFI ABUNDANCE BY SITE
# ============================================================================

cat("Creating Panel D: Total CAFI by site...\n")

panel_d <- ggplot(coral_master, aes(x = site, y = total_cafi, fill = site)) +
  geom_violin(alpha = 0.6, color = NA, width = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", color = "gray30") +
  geom_jitter(aes(color = site), width = 0.1, alpha = 0.4, size = 1.2) +
  scale_fill_manual(values = site_colors, guide = "none") +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(
    x = "Site",
    y = expression(bold("Total CAFI Abundance"))
  ) +
  theme_manuscript(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

cat("   [OK] Panel D created\n")

# ============================================================================
# COMBINE PANELS WITH PATCHWORK
# ============================================================================

cat("\nAssembling figure...\n")

# Combine panels: (A | B) / (C | D)
fig1 <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14)
    )
  ) &
  theme(plot.tag.position = c(0.02, 0.98))

cat("   [OK] Figure assembled\n")

# ============================================================================
# SAVE FIGURE
# ============================================================================

cat("\nSaving figure...\n")

# Ensure output directory exists
dir.create(PATHS$fig_manuscript, recursive = TRUE, showWarnings = FALSE)

# Save
output_path <- file.path(PATHS$fig_manuscript, "fig1_study_design.png")
ggsave(
  output_path,
  fig1,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("   [OK] Saved to:", output_path, "\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    Figure 1 Complete\n")
cat("============================================================\n\n")

cat("Panel contents:\n")
cat("  A: Map of Mo'orea, French Polynesia with 3 sampling sites\n")
cat("  B: Coral volume distribution by site (log10 scale)\n")
cat("  C: Neighbor density by site (n =", n_with_neighbors, "corals with data)\n")
cat("  D: Total CAFI abundance by site\n\n")

cat("Site sample sizes:\n")
coral_master %>%
  group_by(site) %>%
  summarise(
    n_corals = n(),
    n_with_neighbors = sum(!is.na(n_neighbors)),
    median_volume = round(median(volume, na.rm = TRUE)),
    median_cafi = median(total_cafi),
    .groups = "drop"
  ) %>%
  print()

cat("\nOutput:", output_path, "\n")
cat("Dimensions: 10 x 8 inches at 300 dpi\n\n")
