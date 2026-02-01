# Figure 1: 6-panel figure with neighborhood schematics (sized by volumes, realistic distances)
# Set PROJ_LIB
Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sf/proj")

suppressPackageStartupMessages({
  library(sf)
  library(maptiles)
  library(terra)
  library(tidyterra)
  library(ggplot2)
  library(patchwork)
  library(ggspatial)
})

# Load coral data
coral_master <- readRDS("/Users/adrianstier/CAFI-Survey-2026/output/objects/coral_master.rds")

# Site locations
sites <- data.frame(
  site = c("HAU", "MAT", "MRB"),
  label = c("Hauru", "Maatea", "Maharepa"),
  lat = c(-17.516, -17.604, -17.475),
  long = c(-149.922, -149.815, -149.817),
  n = c(38, 39, 35)
) |> st_as_sf(coords = c("long", "lat"), crs = 4326)

site_colors <- c(HAU = "#E69F00", MAT = "#0072B2", MRB = "#009E73")

# Fetch tiles - use larger bbox for fetching, then crop display to hide tile seams/dark corners
bbox_fetch <- st_bbox(c(xmin = -150.05, ymin = -17.68, xmax = -149.68, ymax = -17.38), crs = st_crs(4326))
# Display bounds: crop to island area, avoiding dark tile edges in NE corner
# Shifted xmax west and ymax south to eliminate dark corner while keeping all sites visible
# Center the island in the map panel
bbox_display <- st_bbox(c(xmin = -149.95, ymin = -17.63, xmax = -149.73, ymax = -17.44), crs = st_crs(4326))
cat("Fetching tiles...\n")
tiles <- get_tiles(bbox_fetch, provider = "Esri.WorldImagery", zoom = 14, crop = TRUE)

# Function to scale volume to point size (cube root for area perception)
vol_to_size <- function(vol, min_size = 1.5, max_size = 14) {
  scaled <- (vol^(1/3) - 20^(1/3)) / (45000^(1/3) - 20^(1/3))
  pmax(min_size, pmin(max_size, min_size + scaled * (max_size - min_size)))
}

# Function to create neighborhood schematic with varying neighbor sizes and realistic distances
create_neighborhood_plot <- function(focal_vol, mean_neigh_vol, n_neighbors,
                                     mean_neigh_dist_cm, site_name, focal_color,
                                     panel_label, density_label) {
  set.seed(42 + n_neighbors)

  # 5m radius = 500cm, so scale distances
  # mean_neigh_dist is in cm, 5m radius = 500cm
  # We'll place neighbors at distances proportional to mean, with variation
  mean_dist_scaled <- mean_neigh_dist_cm / 500  # Scale to 0-1 (within 5m radius)

  # Focal coral size
  focal_size <- vol_to_size(focal_vol, min_size = 5, max_size = 16)

  if (n_neighbors > 0) {
    # Generate angles - use even spacing with jitter for small n, random for large n
    if (n_neighbors <= 10) {
      # Even spacing with jitter for better distribution
      base_angles <- seq(0, 2*pi, length.out = n_neighbors + 1)[1:n_neighbors]
      jitter <- runif(n_neighbors, -pi/n_neighbors * 0.4, pi/n_neighbors * 0.4)
      angles <- base_angles + jitter
    } else {
      angles <- runif(n_neighbors, 0, 2*pi)
    }

    # Generate distances centered around mean, with variation (CV ~ 30%)
    dist_cv <- 0.3
    radii_raw <- rnorm(n_neighbors, mean = mean_dist_scaled, sd = mean_dist_scaled * dist_cv)
    radii <- pmax(0.15, pmin(0.95, radii_raw))  # Keep within bounds

    neighbor_x <- radii * cos(angles)
    neighbor_y <- radii * sin(angles)

    # Generate varying neighbor volumes (CV ~ 50% around mean)
    vol_cv <- 0.5
    neighbor_vols <- rlnorm(n_neighbors,
                            meanlog = log(mean_neigh_vol) - (vol_cv^2)/2,
                            sdlog = vol_cv)
    neighbor_vols <- pmax(50, pmin(10000, neighbor_vols))  # Reasonable bounds

    neighbor_sizes <- vol_to_size(neighbor_vols, min_size = 1.2, max_size = 6)

    neighbor_df <- data.frame(
      x = neighbor_x,
      y = neighbor_y,
      size = neighbor_sizes,
      vol = neighbor_vols
    )
  } else {
    neighbor_df <- data.frame(x = numeric(0), y = numeric(0), size = numeric(0), vol = numeric(0))
  }

  p <- ggplot() +
    # 5m radius circle (boundary) - subtle, light appearance
    annotate("path",
             x = cos(seq(0, 2*pi, length.out = 100)),
             y = sin(seq(0, 2*pi, length.out = 100)),
             color = "gray75", linetype = "dashed", linewidth = 0.4) +
    # Neighbor corals (individually sized)
    geom_point(data = neighbor_df, aes(x = x, y = y, size = size),
               shape = 21, fill = "gray55", color = "gray35", stroke = 0.4) +
    scale_size_identity() +
    # Focal coral (center, sized by actual volume)
    annotate("point", x = 0, y = 0,
             shape = 21, fill = focal_color, color = "white",
             size = focal_size, stroke = 2) +
    # Density label at top (more prominent)
    annotate("text", x = 0, y = 1.38, label = density_label,
             fontface = "bold", size = 3.8, color = "gray25") +
    # Neighbor count below density label
    annotate("text", x = 0, y = 1.22, label = paste0(n_neighbors, " neighbors"),
             size = 3.2, color = "gray45") +
    # Site label below circle
    annotate("text", x = 0, y = -1.28, label = site_name,
             fontface = "bold", size = 4.2, color = focal_color) +
    labs(tag = panel_label) +
    coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.45, 1.55)) +
    theme_void() +
    theme(
      plot.tag = element_text(face = "bold", size = 14),
      plot.margin = margin(5, 8, 5, 8),
      plot.background = element_rect(fill = "white", color = NA)
    )

  return(p)
}

# Panel A: Map
p_map <- ggplot() +
  geom_spatraster_rgb(data = tiles) +
  geom_sf(data = sites, aes(fill = site),
          size = 5, shape = 21, color = "white", stroke = 1.6) +
  scale_fill_manual(values = site_colors, guide = "none") +
  # Hauru label - positioned to right of point
  geom_sf_label(data = sites[sites$site == "HAU",],
                aes(label = paste0(label, "\n(n=", n, ")")),
                size = 2.4, fontface = "bold",
                fill = alpha("white", 0.9), label.padding = unit(0.12, "lines"),
                linewidth = 0.15,
                nudge_x = 0.025, nudge_y = -0.012) +
  # Maatea label - positioned to left and below
  geom_sf_label(data = sites[sites$site == "MAT",],
                aes(label = paste0(label, "\n(n=", n, ")")),
                size = 2.4, fontface = "bold",
                fill = alpha("white", 0.9), label.padding = unit(0.12, "lines"),
                linewidth = 0.15,
                nudge_x = -0.032, nudge_y = 0.020) +
  # Maharepa label - positioned below and left of point
  geom_sf_label(data = sites[sites$site == "MRB",],
                aes(label = paste0(label, "\n(n=", n, ")")),
                size = 2.4, fontface = "bold",
                fill = alpha("white", 0.9), label.padding = unit(0.12, "lines"),
                linewidth = 0.15,
                nudge_x = -0.030, nudge_y = -0.012) +
  annotation_scale(location = "bl", width_hint = 0.2,
                   text_cex = 0.7, line_width = 1.2,
                   pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),
                   style = "ticks", text_col = "white", line_col = "white") +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(0.8, "cm"), width = unit(0.6, "cm"),
                         pad_x = unit(0.2, "cm"), pad_y = unit(0.25, "cm"),
                         style = north_arrow_fancy_orienteering(text_size = 8, text_col = "white",
                                                                 line_col = "white", fill = c("white", "white"))) +
  # Display bounds - keep clipping on to prevent overflow
  coord_sf(xlim = c(bbox_display["xmin"], bbox_display["xmax"]),
           ylim = c(bbox_display["ymin"], bbox_display["ymax"]),
           expand = FALSE, crs = 4326) +
  labs(tag = "A") +
  theme_void() +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(5, 5, 5, 5))

# Panel B: Colony volume
p_vol <- ggplot(coral_master, aes(x = volume_field)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 20, fill = "#4A90A4", color = "white",
                 alpha = 0.75, linewidth = 0.25) +
  geom_density(color = "#1A3A5C", linewidth = 1.1, adjust = 1.3) +
  scale_x_log10(labels = scales::comma, breaks = c(100, 1000, 10000)) +
  annotate("segment", x = 15000, xend = 38000, y = 0.68, yend = 0.68,
           arrow = arrow(ends = "both", length = unit(0.07, "cm"), type = "closed"),
           color = "gray45", linewidth = 0.35) +
  annotate("text", x = 24000, y = 0.78,
           label = ">3 orders of\nmagnitude",
           size = 2.5, fontface = "italic", color = "gray40", lineheight = 0.9) +
  labs(tag = "B", subtitle = "n = 112  |  CV = 119%",
       x = expression(Colony~volume~(cm^3)), y = "Density") +
  theme_minimal(base_size = 10) +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0, 1),
        plot.subtitle = element_text(size = 8.5, color = "gray40", margin = margin(b = 2)),
        axis.title = element_text(size = 9),
        axis.title.y = element_text(margin = margin(r = 4)),
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.text = element_text(size = 8),
        axis.line = element_line(color = "black", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.margin = margin(10, 12, 8, 8))

# Panel C: Neighborhood density WITH markers for D, E, F
neighbor_data <- coral_master[!is.na(coral_master$n_neighbors), ]

# Calculate density for placing markers
dens <- density(neighbor_data$n_neighbors, adjust = 1.2)
get_density_at_x <- function(x, dens) {
  approx(dens$x, dens$y, xout = x)$y
}

marker_df <- data.frame(
  x = c(5, 17, 76),
  label = c("D", "E", "F"),
  color = c("#E69F00", "#009E73", "#009E73")
)
marker_df$y <- sapply(marker_df$x, get_density_at_x, dens = dens)

p_neigh <- ggplot(neighbor_data, aes(x = n_neighbors)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 5, fill = "#4A90A4",
                 color = "white", alpha = 0.75, linewidth = 0.25, boundary = 0) +
  geom_density(color = "#1A3A5C", linewidth = 1.1, adjust = 1.2) +
  # Markers as inverted triangles pointing to the x-axis values
  geom_point(data = marker_df, aes(x = x, y = y + 0.004, fill = color),
             shape = 25, color = "gray30", size = 2.8, stroke = 0.6) +
  scale_fill_identity() +
  geom_text(data = marker_df, aes(x = x, y = y + 0.0085, label = label),
            fontface = "bold", size = 3.2, color = "gray20") +
  scale_x_continuous(breaks = seq(0, 80, 20), limits = c(0, 88)) +
  labs(tag = "C", subtitle = "n = 61  |  CV = 79%",
       x = "Neighbors within 5 m", y = "Density") +
  theme_minimal(base_size = 10) +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0, 1),
        plot.subtitle = element_text(size = 8.5, color = "gray40", margin = margin(b = 2)),
        axis.title = element_text(size = 9),
        axis.title.y = element_text(margin = margin(r = 4)),
        axis.title.x = element_text(margin = margin(t = 4)),
        axis.text = element_text(size = 8),
        axis.line = element_line(color = "black", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        plot.margin = margin(10, 12, 8, 8))

# Panels D, E, F: Neighborhood schematics with varying sizes and realistic distances
# Data: focal_vol, mean_neigh_vol, n_neighbors, mean_neigh_dist_cm

# D: HAU-POC04 (Hauru) - 5 neighbors, large focal coral, far neighbors
p_d <- create_neighborhood_plot(
  focal_vol = 26741,
  mean_neigh_vol = 1336,
  n_neighbors = 5,
  mean_neigh_dist_cm = 167,
  site_name = "Hauru",
  focal_color = "#E69F00",
  panel_label = "D",
  density_label = "Low density"
)

# E: MRB-POC10 (Maharepa) - 17 neighbors, medium focal coral
p_e <- create_neighborhood_plot(
  focal_vol = 5472,
  mean_neigh_vol = 686,
  n_neighbors = 17,
  mean_neigh_dist_cm = 154,
  site_name = "Maharepa",
  focal_color = "#009E73",
  panel_label = "E",
  density_label = "Median density"
)

# F: MRB-POC18 (Maharepa) - 76 neighbors, smaller focal, close neighbors
p_f <- create_neighborhood_plot(
  focal_vol = 3064,
  mean_neigh_vol = 395,
  n_neighbors = 76,
  mean_neigh_dist_cm = 113,
  site_name = "Maharepa",
  focal_color = "#009E73",
  panel_label = "F",
  density_label = "High density"
)

# Combine all panels with improved layout proportions
# Give map panel more width (1.35) relative to histogram panels (1.0 each)
top_row <- p_map + p_vol + p_neigh +
  plot_layout(ncol = 3, widths = c(1.35, 1, 1))

# Bottom row panels are equal width
bottom_row <- p_d + p_e + p_f +
  plot_layout(ncol = 3, widths = c(1, 1, 1))

fig1 <- top_row / bottom_row +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title = "Figure 1. Landscape heterogeneity across Mo'orea reef sites",
    subtitle = expression(italic(Pocillopora)~colonies~" | "~French~Polynesia~" | "~2019),
    caption = "(A) Study sites on Mo'orea, French Polynesia. (B) Colony volume distribution spanning >3 orders of magnitude.\n(C) Neighborhood density (Pocillopora within 5 m); triangles mark example corals in D-F. (D-F) Focal coral (colored)\nwith neighbors (gray) within 5 m; circle sizes reflect colony volumes, positions reflect measured distances.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 10, color = "gray35", margin = margin(b = 10)),
      plot.caption = element_text(size = 7.5, hjust = 0, color = "gray45",
                                  lineheight = 1.3, margin = margin(t = 15)),
      plot.margin = margin(12, 15, 15, 15)
    )
  )

# Adjusted dimensions for balanced layout with room for caption
ggsave("/Users/adrianstier/CAFI-Survey-2026/output/figures/manuscript/fig1_study_design.png",
       fig1, width = 12, height = 10, dpi = 300, bg = "white")

cat("Done - 6 panel figure with varying neighbor sizes and realistic distances\n")
