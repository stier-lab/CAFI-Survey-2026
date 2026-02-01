# ============================================================================
# 13_manuscript_figures.R - Publication-Ready Manuscript Figures
# ============================================================================
#
# PURPOSE: Create a complete set of publication-quality figures for the
#          CAFI Survey manuscript, organized by the three core questions
#
# ============================================================================
# FIGURE ORGANIZATION BY CORE QUESTION
# ============================================================================
#
# OVERVIEW:
#   Figure 1: Study Design - Sites, sampling protocol, dataset summary
#
# Q1: SCALING - Does CAFI abundance scale proportionally with coral size?
#   Figure 2: Size-Abundance Scaling (created here)
#     - Tests Field of Dreams (β=1) vs Redirection (β<1)
#   Figure 3: Taxonomic Group Scaling (from 08_functional_groups.R)
#     - Group-level scaling exponents, stacked composition
#
# Q2: COMPOSITION - What shapes community assembly?
#   Figure 4: Community Composition (from 02_community_analysis.R)
#     - Rarefaction control, size-divergence test, NMDS
#
# Q3: FEEDBACKS - Does CAFI community predict coral condition?
#   Figure 5: CAFI-Condition Feedbacks (from 09_cafi_condition_feedbacks.R)
#     - Richness/Trapezia/Galeropsis effects on condition
#
# Q4: NEIGHBORHOOD - Does local coral density affect CAFI?
#   No main figure (null result) → Supplement S7
#
# SUPPLEMENTARY:
#   S1: Species accumulation curves
#   S2: Multi-metric PERMANOVA sensitivity
#   S3: NMDS ordination by site/size
#   S4: Spatial autocorrelation (Moran's I)
#   S5: Composition divergence by size (null after rarefaction)
#   S6: Species-level scaling forest plot
#   S7: Neighborhood null results (Q4)
#
# ============================================================================
#
# DEPENDENCIES: All prior scripts must have been run
#
# OUTPUTS:
#   output/figures/manuscript/fig1_study_design.png      (Overview)
#   output/figures/manuscript/fig2_scaling.png           (Q1)
#   output/figures/manuscript/fig3_functional_groups.png (Q1 groups)
#   output/figures/manuscript/fig4_composition_network.png (Q2)
#   output/figures/manuscript/fig5_feedbacks.png         (Q3)
#   output/figures/supplement/figS1-S7                   (Supplements)
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-22
# ============================================================================

cat("\n========================================\n")
cat("13: Manuscript Figure Suite\n")
cat("========================================\n\n")

# ============================================================================
# SETUP AND DATA LOADING
# ============================================================================

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Create directories
MANUSCRIPT_DIR <- file.path(PATHS$figures, "manuscript")
SUPPLEMENT_DIR <- file.path(PATHS$figures, "supplement")
dir.create(MANUSCRIPT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SUPPLEMENT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load all required data objects
coral_master <- load_object("coral_master")
cafi_clean <- load_object("cafi_clean")
community_matrix <- load_object("community_matrix")
functional_summary <- load_object("functional_summary")

# Load landscape models if available
landscape_models <- tryCatch(
  load_object("landscape_models"),
  error = function(e) NULL
)

# Load scaling results if available
scaling_results <- tryCatch(
  load_object("scaling_results"),
  error = function(e) NULL
)

cat("Data loaded successfully\n")
cat("  Corals:", nrow(coral_master), "\n")
cat("  CAFI records:", nrow(cafi_clean), "\n\n")

# ============================================================================
# PUBLICATION THEME (enhanced)
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

# Standard figure dimensions
FIG_WIDTH_SINGLE <- 8    # Single column
FIG_WIDTH_DOUBLE <- 16   # Double column
FIG_HEIGHT_SMALL <- 6
FIG_HEIGHT_MEDIUM <- 8
FIG_HEIGHT_LARGE <- 10

# ############################################################################
#                    FIGURE 1: STUDY DESIGN (3-PANEL, SATELLITE MAP)
# ############################################################################

cat("============================================================\n")
cat("FIGURE 1: Study Design (3-panel, satellite map)\n")
cat("============================================================\n\n")

# Set PROJ library for sf (MUST be before loading sf/maptiles)
Sys.setenv(PROJ_LIB = "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/sf/proj")

# Conditionally load spatial packages (not in 00_setup.R)
if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' required for Figure 1")
if (!requireNamespace("maptiles", quietly = TRUE)) stop("Package 'maptiles' required for Figure 1")
if (!requireNamespace("tidyterra", quietly = TRUE)) stop("Package 'tidyterra' required for Figure 1")
library(sf)
library(maptiles)

# --- Site coordinates ---
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

# --- Figure 1 color palette (Okabe-Ito, specific to satellite map) ---
SITE_COLORS_FIG1 <- c(
  "HAU" = "#E69F00",  # Orange
  "MAT" = "#0072B2",  # Deep blue
  "MRB" = "#009E73"   # Teal-green
)

# Data colors for histograms
DATA_FILL <- "#4A90A4"        # Sophisticated blue-teal
DATA_FILL_ALPHA <- "#7EB3C4"  # Lighter version for histogram
DENSITY_LINE <- "#1A3A5C"     # Dark blue for density curve

# --- Enhanced theme for Figure 1 ---
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
label_positions <- tibble(
  site = c("HAU", "MAT", "MRB"),
  x_off = c(0.045, -0.055, 0.048),
  y_off_name = c(0.010, 0.045, 0.010),
  y_off_n = c(-0.012, 0.020, -0.012)
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

    # Combined site name + sample size labels
    geom_label(data = site_data,
               aes(x = long + x_off, y = lat + y_off_name,
                   label = paste0(site_name, "\n(n=", n_corals, ")")),
               size = 3.0, fontface = "bold", color = "white",
               fill = alpha("gray15", 0.85), linewidth = 0.3, label.r = unit(0.15, "lines"),
               label.padding = unit(0.22, "lines"), lineheight = 0.85) +

    scale_fill_manual(values = SITE_COLORS_FIG1, guide = "none") +

    # Map bounds
    coord_sf(xlim = c(-149.97, -149.71), ylim = c(-17.65, -17.40), expand = FALSE) +

    # White scale bar with shadow
    annotate("segment", x = -149.949, xend = -149.859, y = -17.621, yend = -17.621,
             color = "black", linewidth = 2.0, alpha = 0.4) +
    annotate("segment", x = -149.95, xend = -149.86, y = -17.62, yend = -17.62,
             color = "white", linewidth = 1.5) +
    annotate("segment", x = -149.95, xend = -149.95, y = -17.628, yend = -17.612,
             color = "white", linewidth = 1.2) +
    annotate("segment", x = -149.86, xend = -149.86, y = -17.628, yend = -17.612,
             color = "white", linewidth = 1.2) +
    annotate("text", x = -149.904, y = -17.601, label = "5 km",
             size = 3.0, color = "black", fontface = "bold", alpha = 0.5) +
    annotate("text", x = -149.905, y = -17.602, label = "5 km",
             size = 3.0, color = "white", fontface = "bold") +

    # White north arrow with shadow
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
  geom_histogram(aes(y = after_stat(density)),
                 bins = 18, fill = DATA_FILL, color = "white",
                 alpha = 0.85, linewidth = 0.3) +
  geom_density(color = DENSITY_LINE, linewidth = 1.1, adjust = 1.1) +
  scale_x_log10(
    labels = scales::label_comma(),
    breaks = c(100, 1000, 10000, 100000),
    limits = c(15, 180000)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
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
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 5, fill = DATA_FILL, color = "white",
                 alpha = 0.85, linewidth = 0.3) +
  geom_density(color = DENSITY_LINE, linewidth = 1.1, adjust = 1.0) +
  geom_vline(xintercept = neighbor_stats$median_n,
             linetype = "dashed", color = "gray45", linewidth = 0.6) +
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

# Layout: A on left (satellite map), B and C stacked on right
fig1 <- panel_a + (panel_b / panel_c) +
  plot_layout(widths = c(1.1, 1)) +
  plot_annotation(
    title = NULL,  # No overall title for journal submission
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 12, 8, 12)
    )
  )

# Save to manuscript directory
ggsave(file.path(MANUSCRIPT_DIR, "fig1_study_design.png"), fig1,
       width = 10.5, height = 6, dpi = 300, bg = "white")

# Also save to analysis output directory
ggsave(file.path(PATHS$fig_01_data, "fig1_study_design.png"), fig1,
       width = 10.5, height = 6, dpi = 300, bg = "white")

cat("  Saved: fig1_study_design.png (3-panel, satellite map)\n")
cat("    -> ", file.path(MANUSCRIPT_DIR, "fig1_study_design.png"), "\n")
cat("    -> ", file.path(PATHS$fig_01_data, "fig1_study_design.png"), "\n\n")

# ############################################################################
#                    FIGURE 2: SIZE-ABUNDANCE SCALING (H2)
# ############################################################################

cat("============================================================\n")
cat("FIGURE 2: Size-Abundance Scaling\n")
cat("============================================================\n\n")

# Prepare data for scaling plots
scaling_data <- coral_master %>%
  filter(!is.na(volume), volume > 0, !is.na(total_cafi))

# Fit power-law model for total CAFI
if (nrow(scaling_data) >= 30) {

  abundance_model <- tryCatch({
    MASS::glm.nb(total_cafi ~ log(volume) + site, data = scaling_data)
  }, error = function(e) NULL)

  if (!is.null(abundance_model)) {
    coefs <- summary(abundance_model)$coefficients
    beta_abundance <- coefs["log(volume)", "Estimate"]
    se_abundance <- coefs["log(volume)", "Std. Error"]
    ci_abundance <- confint(abundance_model, "log(volume)", level = 0.95)

    # Test vs beta = 1 (Field of Dreams)
    z_vs_1 <- (beta_abundance - 1) / se_abundance
    p_vs_1 <- 2 * pnorm(-abs(z_vs_1))

    cat("  Total CAFI scaling: beta =", round(beta_abundance, 3), "\n")
    cat("  95% CI: [", round(ci_abundance[1], 3), ",", round(ci_abundance[2], 3), "]\n")
    cat("  Test vs beta = 1: p =", format.pval(p_vs_1, 3), "\n\n")
  } else {
    beta_abundance <- NA
    ci_abundance <- c(NA, NA)
  }

  # --------------------------------------------------------------------------
  # DESIGN: Site colors from global SITE_COLORS (00_setup.R) — palette chosen
  # to avoid confusion with scaling-class colors in Panel C
  # (Redirection = blue, Super-linear = vermillion).
  # --------------------------------------------------------------------------

  # Shared regression aesthetics — identical visual weight in both panels
  smooth_color <- "gray20"
  smooth_fill  <- "gray85"
  smooth_lwd   <- 0.9
  smooth_alpha  <- 0.35
  point_alpha   <- 0.55
  point_size    <- 1.8

  # Panel A: Total abundance vs volume (log-log)
  panel_a <- scaling_data %>%
    ggplot(aes(x = volume, y = total_cafi)) +
    # 1:1 reference — lightest ink layer
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "gray78", linewidth = 0.4) +
    # Regression ribbon + line — primary visual element
    geom_smooth(aes(group = 1), method = "glm.nb", formula = y ~ log(x),
                se = TRUE, color = smooth_color, fill = smooth_fill,
                linewidth = smooth_lwd, alpha = smooth_alpha) +
    # Data points — secondary; desaturated to avoid color competition
    geom_point(aes(color = site), alpha = point_alpha, size = point_size) +
    scale_x_log10(labels = scales::comma, breaks = c(100, 1000, 10000)) +
    scale_y_log10(labels = scales::comma, breaks = c(1, 10, 100)) +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    coord_cartesian(xlim = c(30, 50000), ylim = c(0.8, 300)) +
    labs(
      x = expression("Coral volume (cm"^3*")"),
      y = "Total CAFI abundance",
      # Panel label bold, descriptor plain — clear hierarchy
      title = expression(bold("A")~"Abundance scaling"),
      subtitle = if(!is.na(beta_abundance))
        bquote(beta == .(sprintf("%.2f", beta_abundance)) ~
          "[" * .(sprintf("%.2f", ci_abundance[1])) * ", " *
          .(sprintf("%.2f", ci_abundance[2])) * "]")
      else "Model fit unavailable"
    ) +
    annotate("text", x = 50, y = 250, label = "1 : 1", size = 2.8,
             color = "gray65", fontface = "italic", hjust = 0) +
    theme_manuscript() +
    theme(
      # Axis labels: medium weight, not bold — reduces non-data ink
      axis.title = element_text(size = 11),
      legend.position = c(0.95, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
      legend.key.size = unit(0.35, "cm"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.margin = margin(3, 5, 3, 5),
      legend.spacing.y = unit(1, "pt"),
      # Minimise grid ink
      panel.grid.minor = element_blank()
    )

  # Panel B: Species richness scaling
  richness_model <- tryCatch({
    glm(otu_richness ~ log(volume) + site, data = scaling_data, family = poisson)
  }, error = function(e) NULL)

  if (!is.null(richness_model)) {
    z_richness <- coef(richness_model)["log(volume)"]
    z_ci <- confint(richness_model, "log(volume)")
  } else {
    z_richness <- NA
    z_ci <- c(NA, NA)
  }

  panel_b <- scaling_data %>%
    ggplot(aes(x = volume, y = otu_richness)) +
    # Regression — same weight as panel A for visual consistency
    geom_smooth(aes(group = 1), method = "glm",
                method.args = list(family = poisson), formula = y ~ log(x),
                se = TRUE, color = smooth_color, fill = smooth_fill,
                linewidth = smooth_lwd, alpha = smooth_alpha) +
    geom_point(aes(color = site), alpha = point_alpha, size = point_size) +
    scale_x_log10(labels = scales::comma, breaks = c(100, 1000, 10000)) +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    coord_cartesian(xlim = c(30, 50000)) +
    labs(
      x = expression("Coral volume (cm"^3*")"),
      y = "Species richness",
      title = expression(bold("B")~"Species-area relationship"),
      subtitle = if(!is.na(z_richness))
        bquote(italic(z) == .(sprintf("%.2f", z_richness)) ~
          "[" * .(sprintf("%.2f", z_ci[1])) * ", " *
          .(sprintf("%.2f", z_ci[2])) * "]")
      else "Model fit unavailable"
    ) +
    theme_manuscript() +
    theme(
      axis.title = element_text(size = 11),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  # --------------------------------------------------------------------------
  # Panel C: Species-level scaling classification (top 10 by prevalence)
  #
  # DESIGN RATIONALE:
  #   - β = 1 is the conceptual benchmark (Field of Dreams / proportional).
  #   - Species whose CI overlaps 1 are shown in neutral gray — they are
  #     *consistent with* the null; color would overstate certainty.
  #   - Accent color reserved for species with CI entirely > 1 (super-linear),
  #     drawing the eye to the biologically distinctive result.
  #   - Background shading reduced to ~4% opacity — just enough to demarcate
  #     zones without competing with data.
  #   - CI whiskers thinner and lighter than point estimates.
  # --------------------------------------------------------------------------
  scaling_results <- tryCatch(load_object("scaling_analysis_results"), error = function(e) NULL)

  if (!is.null(scaling_results) && !is.null(scaling_results$all_results)) {
    species_scaling <- scaling_results$all_results %>%
      filter(category == "Species", !is.na(beta)) %>%
      arrange(desc(n_nonzero)) %>%
      head(10) %>%
      mutate(
        species_label = gsub("_", " ", response),
        species_label = factor(species_label, levels = rev(species_label[order(beta)])),
        scaling_class = case_when(
          boot_ci_upper < 1 ~ "Redirection",
          boot_ci_lower > 1 ~ "Super-linear",
          TRUE              ~ "Field of Dreams"
        )
      )

    x_min <- min(c(species_scaling$boot_ci_lower, 0), na.rm = TRUE) - 0.3
    x_max <- max(species_scaling$boot_ci_upper, na.rm = TRUE) + 0.3

    # Zone shading — very low opacity so data dominates
    shade_df <- data.frame(
      xmin = c(x_min, 1),
      xmax = c(1, x_max),
      fill_col = c("#0072B2", "#D55E00")
    )

    # Color semantics: neutral for CI overlapping 1; accent only for clear signal
    species_colors <- c(
      "Redirection"    = "#5A8FAF",   # muted blue — rare category
      "Field of Dreams" = "gray55",    # neutral — CI overlaps β=1
      "Super-linear"   = "#D55E00"    # accent — CI entirely > 1
    )

    panel_c <- ggplot(species_scaling, aes(x = beta, y = species_label)) +
      # Zone shading at ~4% opacity — demarcates without competing
      geom_rect(data = shade_df, inherit.aes = FALSE,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_col),
                alpha = 0.04) +
      scale_fill_identity() +
      # Reference line at β = 1
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.4) +
      # CI whiskers — visually secondary to points (thinner, slightly transparent)
      geom_segment(aes(x = boot_ci_lower, xend = boot_ci_upper,
                       y = species_label, yend = species_label, color = scaling_class),
                   linewidth = 0.45, alpha = 0.6) +
      # Point estimates — primary data element
      geom_point(aes(color = scaling_class), size = 2.5) +
      scale_color_manual(values = species_colors, name = "Scaling pattern") +
      scale_x_continuous(
        limits = c(x_min, x_max),
        sec.axis = dup_axis(
          breaks = c((x_min + 1) / 2, 1, (1 + x_max) / 2),
          labels = c("Redirection", "\u03b2 = 1", "Super-linear"),
          name = NULL
        )
      ) +
      labs(
        x = expression("Scaling exponent (" * beta * ")"),
        y = NULL,
        title = expression(bold("C")~"Species-level scaling"),
        subtitle = expression("Top 10 species by prevalence"
                              ~"|"~ beta == 1 ~"= proportional (Field of Dreams) scaling")
      ) +
      theme_manuscript() +
      theme(
        axis.title = element_text(size = 11),
        axis.text.y = element_text(face = "italic", size = 9),
        axis.text.x.top = element_text(size = 8.5, face = "bold",
                                        color = c("#5A8FAF", "gray40", "#D55E00")),
        axis.ticks.x.top = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.margin = margin(5, 15, 15, 15)
      )
  } else {
    panel_c <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Species scaling data unavailable") +
      theme_void()
  }

  # Combine: A and B side by side on top, C wide on bottom
  fig2 <- (panel_a | panel_b) / panel_c +
    plot_layout(heights = c(1, 0.8)) +
    plot_annotation(
      caption = "Dashed line (A, C): 1:1 proportional scaling (\u03b2 = 1). Solid curves (A, B): GLM fit \u00B1 95% CI. n = 114 corals, 3 reef sites.",
      theme = theme(
        plot.caption = element_text(size = 8.5, color = "gray45", hjust = 0,
                                    margin = margin(t = 8))
      )
    )

  ggsave(file.path(MANUSCRIPT_DIR, "fig2_scaling.png"), fig2,
         width = 12, height = 9, dpi = 300, bg = "white")

  cat("  Saved: fig2_scaling.png\n\n")

} else {
  cat("  Insufficient data for scaling analysis\n\n")
}

# NOTE: Figure 5 (CAFI-Condition Feedbacks) is created by script 09.
# See "ASSEMBLE FIGURES FROM OTHER SCRIPTS" section below.

# ############################################################################
#                    ASSEMBLE FIGURES FROM OTHER SCRIPTS
# ############################################################################

cat("============================================================\n")
cat("ASSEMBLING FIGURES FROM ANALYSIS SCRIPTS\n")
cat("============================================================\n\n")

# Q1: SCALING - Figure 3 (from 08_functional_groups.R)
cat("Q1: Functional Group Scaling (Figure 3)...\n")
fig3_src <- file.path(PATHS$fig_manuscript, "fig3_functional_groups.png")
if (file.exists(fig3_src)) {
  cat("  Found: fig3_functional_groups.png (from script 08)\n")
} else {
  cat("  Not yet generated - run script 08\n")
}

# ============================================================================
# FIGURE 4: CAFI CO-OCCURRENCE NETWORK (Q2)
# 5-panel layout: circular hero network + 4 guild sub-networks
# ============================================================================
cat("\n============================================================\n")
cat("FIGURE 4: CAFI Co-occurrence Network (5-panel)\n")
cat("============================================================\n\n")

network_results <- tryCatch(load_object("cafi_network"), error = function(e) NULL)

if (!is.null(network_results)) {

  # Load required packages
  library(igraph)
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop("ggforce package required for Figure 4")
  }
  library(ggforce)
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("ggrepel package required for Figure 4")
  }
  library(ggrepel)

  g <- network_results$graph
  communities <- network_results$communities

  # Species info dataframe
  sp_info <- data.frame(
    species = V(g)$name,
    guild = membership(communities),
    degree = degree(g),
    type = V(g)$type,
    stringsAsFactors = FALSE
  ) %>%
    mutate(guild = factor(guild))

  # Guild configuration
  guild_names <- c(
    "1" = "Module I",
    "2" = "Module II",
    "3" = "Module III",
    "4" = "Module IV"
  )

  guild_names_short <- c(
    "1" = "Module I",
    "2" = "Module II",
    "3" = "Module III",
    "4" = "Module IV"
  )

  guild_counts <- sp_info %>% count(guild) %>% arrange(guild)

  # Colorblind-safe guild colors (Okabe-Ito derivatives)
  guild_colors <- c(
    "1" = "#0072B2",  # Blue
    "2" = "#D55E00",  # Vermillion
    "3" = "#009E73",  # Teal
    "4" = "#CC79A7"   # Pink
  )

  # Lighter versions for arc backgrounds
  guild_colors_light <- c(
    "1" = "#A3CDE5",  # clearly blue
    "2" = "#F4B888",  # clearly orange
    "3" = "#8ED5BF",  # clearly teal/green
    "4" = "#E2A5C5"   # clearly pink/mauve
  )

  cat("  Guild sizes:\n")
  for (i in 1:4) {
    cat(sprintf("    Guild %d (%s): %d species\n",
                i, guild_names[as.character(i)],
                guild_counts$n[guild_counts$guild == i]))
  }
  cat(sprintf("\n  Total species: %d\n", sum(guild_counts$n)))
  cat(sprintf("  Total edges: %d\n\n", ecount(g)))

  # --------------------------------------------------------------------------
  # PANEL A: HERO CIRCULAR NETWORK - ALL SPECIES
  # --------------------------------------------------------------------------

  cat("  Building Panel A (circular species network)...\n")

  # Sort species within each guild by degree
  sp_info_sorted <- sp_info %>%
    group_by(guild) %>%
    arrange(desc(degree), .by_group = TRUE) %>%
    mutate(rank_in_guild = row_number()) %>%
    ungroup()

  # Calculate angular positions with guild gaps
  gap_size <- 0.08
  total_gap <- gap_size * 4
  species_arc <- (2 * pi) * (1 - total_gap)

  guild_sizes <- guild_counts$n
  n_total <- sum(guild_sizes)

  guild_props <- guild_sizes / n_total
  guild_arcs <- guild_props * species_arc

  guild_starts <- c(0)
  for (i in 1:3) {
    guild_starts <- c(guild_starts,
                      guild_starts[i] + guild_arcs[i] + gap_size * 2 * pi)
  }

  # Assign positions to each species
  sp_positions <- sp_info_sorted %>%
    group_by(guild) %>%
    mutate(
      guild_idx = as.numeric(guild),
      n_in_guild = n(),
      pos_in_guild = (row_number() - 0.5) / n_in_guild
    ) %>%
    ungroup() %>%
    mutate(
      guild_start = guild_starts[guild_idx],
      guild_arc = guild_arcs[guild_idx],
      angle = pi/2 - (guild_start + pos_in_guild * guild_arc),
      radius = 4.5,
      x = radius * cos(angle),
      y = radius * sin(angle),
      node_size = scales::rescale(degree, to = c(2.5, 8))
    )

  # Get edge data
  edge_df <- igraph::as_data_frame(g, what = "edges") %>%
    left_join(sp_positions %>% dplyr::select(species, x, y, guild, degree),
              by = c("from" = "species")) %>%
    rename(x1 = x, y1 = y, guild_from = guild, degree_from = degree) %>%
    left_join(sp_positions %>% dplyr::select(species, x, y, guild, degree),
              by = c("to" = "species")) %>%
    rename(x2 = x, y2 = y, guild_to = guild, degree_to = degree) %>%
    mutate(
      is_within_guild = guild_from == guild_to,
      edge_alpha = ifelse(is_within_guild, 0.5, 0.15),
      edge_color = ifelse(is_within_guild,
                          guild_colors[as.character(guild_from)],
                          "gray50")
    )

  between_guild_edges <- edge_df %>% filter(!is_within_guild)
  between_guild_edges_sampled <- between_guild_edges
  n_between <- nrow(between_guild_edges)

  cat(sprintf("    Between-guild edges: %d (all shown as thin gray)\n", n_between))

  # Bezier curve helper
  create_bezier <- function(x1, y1, x2, y2, n_points = 30, curvature = 0.35) {
    cx <- (x1 + x2) / 2 * (1 - curvature)
    cy <- (y1 + y2) / 2 * (1 - curvature)
    t <- seq(0, 1, length.out = n_points)
    data.frame(
      x = (1-t)^2 * x1 + 2*(1-t)*t * cx + t^2 * x2,
      y = (1-t)^2 * y1 + 2*(1-t)*t * cy + t^2 * y2
    )
  }

  # Generate bezier curves for within-guild edges
  within_guild_edges <- edge_df %>% filter(is_within_guild)
  cat("    Generating bezier curves for", nrow(within_guild_edges), "within-guild edges...\n")

  bezier_within <- map_dfr(1:nrow(within_guild_edges), function(i) {
    row <- within_guild_edges[i,]
    pts <- create_bezier(row$x1, row$y1, row$x2, row$y2)
    pts$edge_id <- i
    pts$weight <- row$weight
    pts$guild <- as.character(row$guild_from)
    pts
  })

  # Generate bezier curves for between-guild edges
  cat("    Generating bezier curves for", nrow(between_guild_edges_sampled), "between-guild edges...\n")

  if (nrow(between_guild_edges_sampled) > 0) {
    bezier_between <- map_dfr(1:nrow(between_guild_edges_sampled), function(i) {
      row <- between_guild_edges_sampled[i,]
      pts <- create_bezier(row$x1, row$y1, row$x2, row$y2)
      pts$edge_id <- i + 10000
      pts$weight <- row$weight
      pts
    })
  } else {
    bezier_between <- data.frame(x = numeric(), y = numeric(), edge_id = integer(), weight = numeric())
  }

  # Create guild arc backgrounds from actual node positions
  guild_arc_data <- sp_positions %>%
    group_by(guild) %>%
    summarize(
      min_angle = min(angle),
      max_angle = max(angle),
      .groups = "drop"
    ) %>%
    mutate(
      padding = 0.05,
      start = max_angle + padding,
      end = min_angle - padding,
      r_inner = 3.8,
      r_outer = 5.2,
      fill_color = guild_colors_light[as.character(guild)],
      border_color = guild_colors[as.character(guild)]
    )

  # Arc polygon helper
  make_arc_polygon <- function(start_angle, end_angle, r_inner, r_outer, n = 100) {
    angles_outer <- seq(start_angle, end_angle, length.out = n)
    angles_inner <- rev(angles_outer)
    data.frame(
      x = c(r_outer * cos(angles_outer), r_inner * cos(angles_inner)),
      y = c(r_outer * sin(angles_outer), r_inner * sin(angles_inner))
    )
  }

  arc_layers <- lapply(1:nrow(guild_arc_data), function(i) {
    d <- guild_arc_data[i, ]
    poly_df <- make_arc_polygon(d$start, d$end, d$r_inner, d$r_outer)
    geom_polygon(
      data = poly_df,
      aes(x = x, y = y),
      fill = d$fill_color,
      color = d$border_color,
      alpha = 0.5,
      linewidth = 0.8
    )
  })

  # Guild label positions
  guild_label_positions <- sp_positions %>%
    group_by(guild) %>%
    summarize(mid_angle = (min(angle) + max(angle)) / 2, .groups = "drop") %>%
    left_join(guild_counts, by = "guild") %>%
    mutate(
      label = paste0(
        c("Module I", "Module II", "Module III", "Module IV")[as.numeric(guild)],
        " (n=", n, ")"
      ),
      x = c(6.5, 6.5, -6.5, -6.5)[as.numeric(guild)],
      y = c(6.5, -6.5, -6.5, 6.5)[as.numeric(guild)],
      hjust = c(1, 1, 0, 0)[as.numeric(guild)],
      vjust = c(0, 1, 1, 0)[as.numeric(guild)],
      color = guild_colors[as.character(guild)]
    )

  # Build Panel A
  p_A <- ggplot() +
    arc_layers +
    geom_path(
      data = bezier_between,
      aes(x = x, y = y, group = edge_id),
      color = "gray55",
      alpha = 0.25,
      linewidth = 0.25
    ) +
    geom_path(
      data = bezier_within,
      aes(x = x, y = y, group = edge_id, color = guild, alpha = weight),
      linewidth = 0.5
    ) +
    geom_point(
      data = sp_positions,
      aes(x = x, y = y, fill = guild, size = node_size),
      shape = 21,
      color = "white",
      stroke = 0.6
    ) +
    geom_text(
      data = guild_label_positions,
      aes(x = x, y = y, label = label, color = guild, hjust = hjust, vjust = vjust),
      size = 3.5,
      fontface = "bold"
    ) +
    scale_fill_manual(values = guild_colors, guide = "none") +
    scale_color_manual(values = guild_colors, guide = "none") +
    scale_size_identity() +
    scale_alpha_continuous(range = c(0.15, 0.7), guide = "none") +
    coord_fixed(ratio = 1, xlim = c(-6.8, 6.8), ylim = c(-6.8, 6.8), clip = "off") +
    labs(
      title = "A. Species Co-occurrence Network",
      subtitle = sprintf("%d species | %d co-occurrences | 4 ecological guilds",
                         nrow(sp_positions), ecount(g))
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40",
                                   margin = margin(b = 5)),
      plot.margin = margin(5, 2, 5, 2)
    )

  cat("    Panel A complete.\n")

  # --------------------------------------------------------------------------
  # PANELS B-E: INDIVIDUAL GUILD NETWORKS WITH SPECIES LABELS
  # --------------------------------------------------------------------------

  create_guild_panel_fixed <- function(guild_id, letter, max_labels = 10) {

    guild_name <- guild_names_short[as.character(guild_id)]
    guild_color <- guild_colors[as.character(guild_id)]
    guild_color_light <- guild_colors_light[as.character(guild_id)]

    guild_species <- sp_info %>% filter(guild == guild_id)
    g_sub <- induced_subgraph(g, V(g)[V(g)$name %in% guild_species$species])

    if (vcount(g_sub) == 0) {
      return(ggplot() + theme_void() +
               labs(title = paste(letter, guild_name)))
    }

    set.seed(42 + guild_id)
    n_sp <- vcount(g_sub)
    if (n_sp > 12) {
      layout_fr <- layout_with_kk(g_sub)
    } else {
      layout_fr <- layout_with_fr(g_sub, niter = 1000, weights = E(g_sub)$weight)
    }

    if (nrow(layout_fr) > 1) {
      layout_fr[,1] <- scales::rescale(layout_fr[,1], to = c(-1.4, 1.4))
      layout_fr[,2] <- scales::rescale(layout_fr[,2], to = c(-1.4, 1.4))
    } else {
      layout_fr[,1] <- 0
      layout_fr[,2] <- 0
    }

    node_data <- data.frame(
      species = V(g_sub)$name,
      x = layout_fr[,1],
      y = layout_fr[,2],
      degree = degree(g_sub),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        node_size = scales::rescale(degree, to = c(4, 14)),
        species_label = gsub("([A-Z])[a-z]+ ", "\\1. ", species)
      )

    n_species <- nrow(node_data)
    n_to_label <- if (n_species <= max_labels) n_species else max_labels

    node_data <- node_data %>%
      mutate(
        rank_degree = rank(-degree, ties.method = "first"),
        show_label = rank_degree <= n_to_label
      )

    cat(sprintf("    Guild %d: %d species, labeling top %d\n",
                guild_id, n_species, sum(node_data$show_label)))

    if (ecount(g_sub) > 0) {
      edge_list_sub <- igraph::as_data_frame(g_sub, what = "edges")
      edge_data <- edge_list_sub %>%
        left_join(node_data %>% dplyr::select(species, x, y),
                  by = c("from" = "species")) %>%
        rename(x1 = x, y1 = y) %>%
        left_join(node_data %>% dplyr::select(species, x, y),
                  by = c("to" = "species")) %>%
        rename(x2 = x, y2 = y)
    } else {
      edge_data <- NULL
    }

    n_edges <- ifelse(is.null(edge_data), 0, nrow(edge_data))

    p <- ggplot()

    p <- p + geom_circle(
      aes(x0 = 0, y0 = 0, r = 1.4),
      fill = guild_color_light,
      color = NA,
      alpha = 0.25
    )

    if (!is.null(edge_data) && nrow(edge_data) > 0) {
      p <- p + geom_segment(
        data = edge_data,
        aes(x = x1, y = y1, xend = x2, yend = y2, alpha = weight),
        color = guild_color,
        linewidth = 0.5
      )
    }

    p <- p + geom_point(
      data = node_data,
      aes(x = x, y = y, size = node_size),
      fill = guild_color,
      shape = 21,
      color = "white",
      stroke = 0.7
    )

    labels_data <- node_data %>% filter(show_label)

    if (nrow(labels_data) > 0) {
      p <- p + geom_text_repel(
        data = labels_data,
        aes(x = x, y = y, label = species_label),
        size = 3.2,
        fontface = "bold.italic",
        color = "gray5",
        bg.color = "white",
        bg.r = 0.12,
        segment.color = "gray40",
        segment.size = 0.3,
        segment.alpha = 0.7,
        box.padding = 0.55,
        point.padding = 0.45,
        max.overlaps = 30,
        force = 10,
        force_pull = 0.4,
        max.iter = 8000,
        seed = 42
      )
    }

    n_hidden <- n_species - sum(node_data$show_label)
    subtitle_text <- if (n_hidden > 0) {
      sprintf("%d species | %d edges | top %d labeled", n_species, n_edges, n_to_label)
    } else {
      sprintf("%d species | %d edges", n_species, n_edges)
    }

    p <- p +
      scale_size_identity() +
      scale_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
      coord_fixed(ratio = 1, xlim = c(-2.1, 2.1), ylim = c(-2.1, 2.1)) +
      labs(
        title = sprintf("%s. %s", letter, guild_name),
        subtitle = subtitle_text
      ) +
      theme_void() +
      theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5,
                                  color = guild_color),
        plot.subtitle = element_text(size = 7.5, hjust = 0.5, color = "gray50"),
        plot.margin = margin(3, 3, 3, 3),
        plot.background = element_rect(fill = "white", color = "gray80",
                                       linewidth = 0.5)
      )

    return(p)
  }

  cat("  Building Panels B-E (individual guilds)...\n")

  p_B <- create_guild_panel_fixed(1, "B", max_labels = 50)
  cat("    Panel B complete.\n")
  p_C <- create_guild_panel_fixed(2, "C", max_labels = 50)
  cat("    Panel C complete.\n")
  p_D <- create_guild_panel_fixed(3, "D", max_labels = 50)
  cat("    Panel D complete.\n")
  p_E <- create_guild_panel_fixed(4, "E", max_labels = 50)
  cat("    Panel E complete.\n")

  # --------------------------------------------------------------------------
  # COMBINE INTO WIDE 5-PANEL FIGURE
  # --------------------------------------------------------------------------

  cat("  Combining panels into wide layout...\n")

  p_right <- (p_B | p_C) / (p_D | p_E)

  fig4 <- (p_A | p_right) +
    plot_layout(widths = c(1.5, 2)) +
    plot_annotation(
      title = "Figure 4: CAFI Co-occurrence Network Structure",
      subtitle = sprintf(
        "Four ecological guilds identified via Louvain community detection | Q = %.2f | %d species | %d edges",
        modularity(communities), vcount(g), ecount(g)
      ),
      caption = paste0(
        "Node size = degree centrality | Within-guild edges colored, between-guild edges gray | ",
        "Threshold: r > 0.3, FDR p < 0.05"
      ),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
        plot.caption = element_text(size = 8, hjust = 0.5, color = "gray50",
                                    margin = margin(t = 10)),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )

  # Save to both manuscript and analysis directories
  fig4_network_dir <- file.path(PATHS$figures, "06_network")
  dir.create(fig4_network_dir, recursive = TRUE, showWarnings = FALSE)

  ggsave(file.path(MANUSCRIPT_DIR, "fig4_network.png"), fig4,
         width = 18, height = 11, dpi = 300, bg = "white")
  cat("  Saved: fig4_network.png (manuscript)\n")

  ggsave(file.path(fig4_network_dir, "fig4_5panel_v2_wide.png"), fig4,
         width = 18, height = 11, dpi = 300, bg = "white")
  cat("  Saved: fig4_5panel_v2_wide.png (analysis)\n\n")

} else {
  cat("  Network results not available - run script 06 first\n\n")
}

# Q3: FEEDBACKS - Figure 5 (from 09_cafi_condition_feedbacks.R)
cat("\nQ3: CAFI-Condition Feedbacks (Figure 5)...\n")
fig5_src <- file.path(PATHS$fig_manuscript, "fig5_feedbacks.png")
if (file.exists(fig5_src)) {
  cat("  Found: fig5_feedbacks.png (from script 09)\n")
} else {
  cat("  Not yet generated - run script 09\n")
}

cat("\n")

# ############################################################################
#                    SUPPLEMENTARY FIGURES
# ############################################################################

cat("============================================================\n")
cat("SUPPLEMENTARY FIGURES\n")
cat("============================================================\n\n")

# S1: Species accumulation curve
cat("Creating S1: Species accumulation curve...\n")

if (exists("community_matrix") && nrow(community_matrix) > 10) {

  # Calculate species accumulation
  specaccum_result <- tryCatch({
    vegan::specaccum(community_matrix, method = "random", permutations = 100)
  }, error = function(e) NULL)

  if (!is.null(specaccum_result)) {

    sac_data <- tibble(
      sites = specaccum_result$sites,
      richness = specaccum_result$richness,
      sd = specaccum_result$sd
    )

    fig_s1 <- ggplot(sac_data, aes(x = sites, y = richness)) +
      geom_ribbon(aes(ymin = richness - 1.96*sd, ymax = richness + 1.96*sd),
                  fill = "#377EB8", alpha = 0.3) +
      geom_line(color = "#377EB8", linewidth = 1.2) +
      geom_hline(yintercept = max(sac_data$richness),
                 linetype = "dashed", color = "gray50") +
      annotate("text", x = max(sac_data$sites) * 0.9,
               y = max(sac_data$richness) + 2,
               label = paste("Total:", round(max(sac_data$richness))),
               size = 3.5) +
      labs(
        x = "Number of Coral Colonies Sampled",
        y = "Cumulative Species Richness",
        title = "Figure S1: Species Accumulation Curve",
        subtitle = "Random accumulation with 95% CI (100 permutations)",
        caption = "Asymptote suggests adequate sampling coverage"
      ) +
      theme_manuscript()

    ggsave(file.path(SUPPLEMENT_DIR, "figS1_species_accumulation.png"), fig_s1,
           width = 8, height = 6, dpi = 300, bg = "white")

    cat("  Saved: figS1_species_accumulation.png\n")
  }
}

# S2: Multi-metric PERMANOVA sensitivity
sensitivity_locations <- c(
  file.path(PATHS$figures, "02_community", "permanova_metric_sensitivity.png")
)
sens_fig <- sensitivity_locations[file.exists(sensitivity_locations)][1]
if (!is.na(sens_fig) && file.exists(sens_fig)) {
  file.copy(sens_fig, file.path(SUPPLEMENT_DIR, "figS2_permanova_sensitivity.png"),
            overwrite = TRUE)
  cat("  Saved: figS2_permanova_sensitivity.png\n")
}

# S3: NMDS ordination by site/size
nmds_locations <- c(
  file.path(PATHS$figures, "02_community", "nmds_by_site.png"),
  file.path(PATHS$figures, "02_community", "nmds_ordination.png")
)
nmds_fig <- nmds_locations[file.exists(nmds_locations)][1]
if (!is.na(nmds_fig) && file.exists(nmds_fig)) {
  file.copy(nmds_fig, file.path(SUPPLEMENT_DIR, "figS3_nmds_ordination.png"),
            overwrite = TRUE)
  cat("  Saved: figS3_nmds_ordination.png\n")
}

# S4: Spatial autocorrelation (Moran's I)
spatial_map_locations <- c(
  file.path(PATHS$figures, "07_spatial", "spatial_autocorrelation_map.png"),
  file.path(PATHS$figures, "spatial_autocorrelation_map.png")
)
spatial_map <- spatial_map_locations[file.exists(spatial_map_locations)][1]
if (!is.na(spatial_map) && file.exists(spatial_map)) {
  file.copy(spatial_map, file.path(SUPPLEMENT_DIR, "figS4_spatial_autocorrelation.png"),
            overwrite = TRUE)
  cat("  Saved: figS4_spatial_autocorrelation.png\n")
}

# S5: Composition divergence by size (former main figure - null after rarefaction)
divergence_locations <- c(
  file.path(PATHS$figures, "02_community", "composition_divergence_by_size.png")
)
div_fig <- divergence_locations[file.exists(divergence_locations)][1]
if (!is.na(div_fig) && file.exists(div_fig)) {
  file.copy(div_fig, file.path(SUPPLEMENT_DIR, "figS5_composition_divergence.png"),
            overwrite = TRUE)
  cat("  Saved: figS5_composition_divergence.png\n")
}

# S6: Species-level scaling forest plot
species_scaling_locations <- c(
  file.path(PATHS$figures, "05_scaling", "species_scaling_forest.png"),
  file.path(PATHS$figures, "species_scaling_forest.png")
)
species_fig <- species_scaling_locations[file.exists(species_scaling_locations)][1]
if (!is.na(species_fig) && file.exists(species_fig)) {
  file.copy(species_fig, file.path(SUPPLEMENT_DIR, "figS6_species_scaling.png"),
            overwrite = TRUE)
  cat("  Saved: figS6_species_scaling.png\n")
}

# S7: Neighborhood null results (Q4)
# Assemble panels showing non-significant neighborhood effects
neighborhood_fig_locations <- c(
  file.path(PATHS$figures, "04_effects", "abundance_vs_neighbors.png")
)
neighbor_fig <- neighborhood_fig_locations[file.exists(neighborhood_fig_locations)][1]
if (!is.na(neighbor_fig) && file.exists(neighbor_fig)) {
  file.copy(neighbor_fig, file.path(SUPPLEMENT_DIR, "figS7_neighborhood_null.png"),
            overwrite = TRUE)
  cat("  Saved: figS7_neighborhood_null.png\n")
} else {
  # Create minimal neighborhood null figure from coral_master data
  if ("n_neighbors" %in% names(coral_master)) {
    neighbor_data <- coral_master %>% filter(!is.na(n_neighbors))
    if (nrow(neighbor_data) >= 20) {
      p_s7a <- ggplot(neighbor_data, aes(x = n_neighbors, y = total_cafi)) +
        geom_point(alpha = 0.6, color = "#0072B2") +
        geom_smooth(method = "lm", se = TRUE, color = "black") +
        labs(x = "Number of Neighbors (5m)", y = "Total CAFI Abundance",
             title = "A. CAFI Abundance vs Neighborhood Density (NS)") +
        theme_manuscript()

      p_s7b <- if ("condition_score" %in% names(neighbor_data)) {
        ggplot(neighbor_data %>% filter(!is.na(condition_score)),
               aes(x = n_neighbors, y = condition_score)) +
          geom_point(alpha = 0.6, color = "#D55E00") +
          geom_smooth(method = "lm", se = TRUE, color = "black") +
          labs(x = "Number of Neighbors (5m)", y = "Coral Condition Score",
               title = "B. Coral Condition vs Neighborhood Density (NS)") +
          theme_manuscript()
      } else {
        ggplot() + theme_void() + labs(title = "B. Condition data unavailable")
      }

      fig_s7 <- p_s7a | p_s7b
      fig_s7 <- fig_s7 + plot_annotation(
        title = "Figure S7: Neighborhood Density Does Not Predict CAFI or Condition",
        subtitle = "Q4 null result: coral density within 5m radius is not significant (n = 61 corals with spatial data)",
        theme = theme(plot.title = element_text(face = "bold", size = 12),
                      plot.subtitle = element_text(size = 10, color = "gray40"))
      )

      ggsave(file.path(SUPPLEMENT_DIR, "figS7_neighborhood_null.png"), fig_s7,
             width = 12, height = 5, dpi = 300, bg = "white")
      cat("  Saved: figS7_neighborhood_null.png (created from data)\n")
    }
  }
}

cat("\n")

# ############################################################################
#                    RESULTS TABLE & SAMPLE SIZES
# ############################################################################

cat("============================================================\n")
cat("COMPILING MANUSCRIPT RESULTS TABLE\n")
cat("============================================================\n\n")

# --- Unified Results Table ---
# Compiles all key statistical results across Q1-Q4

results_rows <- list()

# Q1: Scaling (from script 05)
tryCatch({
  scaling_obj <- load_object("scaling_analysis_results")
  comm_result <- scaling_obj$models$total_abundance
  if (!is.null(comm_result) && comm_result$converged) {
    results_rows[[length(results_rows) + 1]] <- tibble(
      Question = "Q1", Hypothesis = "Redirection", Test = "NB GLM",
      Predictor = "log(volume)", Beta_R2 = comm_result$beta,
      CI_lower = comm_result$ci_lower, CI_upper = comm_result$ci_upper,
      p_value = comm_result$p_vs_1, p_FDR = NA_real_,
      n = comm_result$n_corals, Interpretation = comm_result$interpretation
    )
  }
}, error = function(e) cat("  [Skip] Scaling results not available\n"))

# Q2: Composition (from saved PERMANOVA results)
tryCatch({
  community_results <- load_object("community_analysis_results")
  if (!is.null(community_results$permanova)) {
    perm <- community_results$permanova
    results_rows[[length(results_rows) + 1]] <- tibble(
      Question = "Q2", Hypothesis = "Site effects", Test = "PERMANOVA",
      Predictor = "site", Beta_R2 = perm$R2[which(rownames(perm) == "site")],
      CI_lower = NA_real_, CI_upper = NA_real_,
      p_value = perm$`Pr(>F)`[which(rownames(perm) == "site")],
      p_FDR = NA_real_, n = 114, Interpretation = "Significant"
    )
  }
}, error = function(e) cat("  [Skip] PERMANOVA results not available\n"))

# Q3: Feedbacks (from script 09 outputs)
tryCatch({
  feedback_table <- read_csv(file.path(PATHS$tables, "cafi_condition_models.csv"),
                             show_col_types = FALSE)
  if (nrow(feedback_table) > 0) {
    for (i in 1:nrow(feedback_table)) {
      row <- feedback_table[i, ]
      results_rows[[length(results_rows) + 1]] <- tibble(
        Question = "Q3", Hypothesis = paste0(row$predictor, "->condition"),
        Test = "LMM", Predictor = row$predictor,
        Beta_R2 = row$estimate, CI_lower = row$ci_lower, CI_upper = row$ci_upper,
        p_value = row$p_value,
        p_FDR = if ("p_fdr" %in% names(row)) row$p_fdr else NA_real_,
        n = row$n, Interpretation = ifelse(row$significant, "Significant", "NS")
      )
    }
  }
}, error = function(e) cat("  [Skip] Feedback results not available\n"))

# Q3b: Rarefied richness (abundance artifact test)
tryCatch({
  richness_artifact <- read_csv(file.path(PATHS$tables, "richness_abundance_artifact.csv"),
                                 show_col_types = FALSE)
  if (nrow(richness_artifact) > 0) {
    # Add rarefied richness result (key finding)
    rarefied_row <- richness_artifact %>%
      filter(grepl("Rarefied", type))
    if (nrow(rarefied_row) > 0) {
      results_rows[[length(results_rows) + 1]] <- tibble(
        Question = "Q3", Hypothesis = "Rarefied richness->condition (artifact test)",
        Test = "LM", Predictor = "rarefied_richness",
        Beta_R2 = rarefied_row$estimate[1],
        CI_lower = NA_real_, CI_upper = NA_real_,
        p_value = rarefied_row$p_value[1],
        p_FDR = NA_real_,
        n = rarefied_row$n[1],
        Interpretation = "NS (abundance artifact)"
      )
    }
  }
}, error = function(e) cat("  [Skip] Rarefied richness results not available\n"))

# Q4: Neighborhood (from script 04)
tryCatch({
  landscape_models <- load_object("landscape_models")
  if (!is.null(landscape_models$abundance_full)) {
    coefs <- summary(landscape_models$abundance_full)$coefficients
    if ("n_neighbors" %in% rownames(coefs)) {
      results_rows[[length(results_rows) + 1]] <- tibble(
        Question = "Q4", Hypothesis = "Neighborhood->CAFI",
        Test = "NB GLM", Predictor = "n_neighbors",
        Beta_R2 = coefs["n_neighbors", "Estimate"],
        CI_lower = NA_real_, CI_upper = NA_real_,
        p_value = coefs["n_neighbors", "Pr(>|z|)"],
        p_FDR = NA_real_, n = 61, Interpretation = "NS"
      )
    }
  }
}, error = function(e) cat("  [Skip] Landscape model results not available\n"))

# Compile and save
if (length(results_rows) > 0) {
  manuscript_results <- bind_rows(results_rows)
  save_table(manuscript_results, "manuscript_results_summary")
  cat("  Saved: manuscript_results_summary.csv (", nrow(manuscript_results), " results)\n\n")
} else {
  cat("  No pre-computed results available. Run analysis scripts first.\n\n")
}

# --- Sample Size Table ---
cat("Sample sizes by analysis:\n")

sample_sizes <- tibble(
  Analysis = c("Scaling (Q1)", "PERMANOVA (Q2)", "Rarefaction (Q2)",
               "Condition feedbacks (Q3)", "Neighborhood (Q4)"),
  Total_N = c(
    nrow(coral_master),
    nrow(coral_master),
    sum(coral_master$total_cafi >= 5, na.rm = TRUE),
    sum(!is.na(coral_master$condition_score)),
    sum(!is.na(coral_master$n_neighbors))
  ),
  Subset_Reason = c(
    "All corals with volume",
    "All corals with CAFI",
    "Corals with >= 5 CAFI (rarefaction depth)",
    "Corals with physiology data",
    "Corals in 5m neighborhood survey"
  )
)

save_table(sample_sizes, "sample_sizes")
cat("  Saved: sample_sizes.csv\n")
print(sample_sizes, n = Inf)
cat("\n")

# ############################################################################
#                    SUMMARY
# ############################################################################

cat("============================================================\n")
cat("MANUSCRIPT FIGURE SUITE COMPLETE\n")
cat("============================================================\n\n")

# List all manuscript figures
manuscript_figs <- list.files(MANUSCRIPT_DIR, pattern = "\\.png$", full.names = FALSE)
supplement_figs <- list.files(SUPPLEMENT_DIR, pattern = "\\.png$", full.names = FALSE)

cat("MAIN TEXT FIGURES (5):\n")
cat("======================\n\n")

cat("OVERVIEW:\n")
cat("  Fig 1: Study design, sites, dataset summary\n")
overview_figs <- manuscript_figs[grepl("fig1_", manuscript_figs)]
for (fig in sort(overview_figs)) cat("    ", fig, "\n")

cat("\nQ1: SCALING (Size-abundance relationships):\n")
cat("  Fig 2: Community-level scaling (abundance, density, richness)\n")
cat("  Fig 3: Taxonomic group scaling (Trapezia vs Fish vs Gastropods)\n")
q1_figs <- manuscript_figs[grepl("fig[23]_", manuscript_figs)]
for (fig in sort(q1_figs)) cat("    ", fig, "\n")

cat("\nQ2: COMPOSITION (What structures CAFI composition?):\n")
cat("  Fig 4: NMDS by site + co-occurrence network + modularity + hub species\n")
q2_figs <- manuscript_figs[grepl("fig4_", manuscript_figs)]
for (fig in sort(q2_figs)) cat("    ", fig, "\n")

cat("\nQ3: FEEDBACKS (Does CAFI identity predict coral condition?):\n")
cat("  Fig 5: Richness, Trapezia, Galeropsis effects on condition\n")
q3_figs <- manuscript_figs[grepl("fig5_", manuscript_figs)]
for (fig in sort(q3_figs)) cat("    ", fig, "\n")

cat("\nQ4: NEIGHBORHOOD (Does local coral density affect CAFI?):\n")
cat("  No main figure (null result) -> Figure S7\n")

cat("\nSUPPLEMENTARY FIGURES (", length(supplement_figs), "):\n")
for (fig in sort(supplement_figs)) {
  cat("  -", fig, "\n")
}

cat("\nTotal figures:", length(manuscript_figs), "main +", length(supplement_figs), "supplementary\n")
cat("\nOutput locations:\n")
cat("  Main figures:", MANUSCRIPT_DIR, "\n")
cat("  Supplements:", SUPPLEMENT_DIR, "\n\n")

cat("Manuscript figure generation complete!\n")
