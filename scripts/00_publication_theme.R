#!/usr/bin/env Rscript
# ============================================================================
# 00_publication_theme.R - Publication-Quality Figure Style Guide
#
# Purpose: Define consistent styling for all figures in the CAFI manuscript
# Target Journal: Marine Ecology Progress Series (MEPS)
#
# This file should be sourced after 00_load_libraries.R in all figure scripts
#
# Author: CAFI Analysis Pipeline
# Date: 2025-11-23
# ============================================================================

cat("Loading publication theme...\n")

# ============================================================================
# COLOR PALETTES
# ============================================================================

# Site colors - distinct, colorblind-friendly
site_colors <- c(

"HAU" = "#E69F00",
  "MAT" = "#56B4E9",
"MRB" = "#009E73"
)

# Branch architecture colors
branch_colors <- c(
  "tight" = "#D55E00",
  "wide" = "#0072B2"
)

# Taxonomic group colors
taxon_colors <- c(
  "crab" = "#CC79A7",
  "shrimp" = "#F0E442",
  "snail" = "#999999",
  "fish" = "#0072B2"
)

# Gradient palettes for continuous variables
# Use viridis for accessibility
gradient_low <- "#440154"
gradient_high <- "#FDE725"

# ============================================================================
# TYPOGRAPHY
# ============================================================================

# Base font size for figures (will scale with figure dimensions)
base_size <- 12

# Font family - use system fonts that work across platforms
# For publication, Times New Roman or Helvetica are common
font_family <- "sans"  # Use "serif" for Times-like

# ============================================================================
# PUBLICATION THEME
# ============================================================================

theme_publication <- function(base_size = 12, base_family = "sans") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Panel
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Axes
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text = element_text(size = rel(0.9), color = "black"),
      axis.title = element_text(size = rel(1), face = "plain"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),

      # Legend
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.8, "lines"),
      legend.text = element_text(size = rel(0.8)),
      legend.title = element_text(size = rel(0.9), face = "plain"),
      legend.position = "right",
      legend.margin = margin(0, 0, 0, 0),

      # Facets
      strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.3),
      strip.text = element_text(size = rel(0.9), face = "plain",
                                margin = margin(3, 3, 3, 3)),

      # Plot title and captions
      plot.title = element_text(size = rel(1.1), face = "bold",
                                hjust = 0, margin = margin(b = 8)),
      plot.subtitle = element_text(size = rel(0.9), hjust = 0,
                                   margin = margin(b = 8), color = "gray30"),
      plot.caption = element_text(size = rel(0.7), hjust = 1,
                                  color = "gray50"),
      plot.margin = margin(10, 10, 10, 10),

      # Overall
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Standard ggsave wrapper for consistent output
save_publication_figure <- function(plot, filename, width = 10, height = 8,
                                    dpi = 300, ...) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    ...
  )
}

# Scale functions for consistent sizing
scale_color_site <- function(...) {
  scale_color_manual(values = site_colors, name = "Site", ...)
}

scale_fill_site <- function(...) {
  scale_fill_manual(values = site_colors, name = "Site", ...)
}

scale_color_branch <- function(...) {
  scale_color_manual(values = branch_colors, name = "Branch\narchitecture", ...)
}

scale_fill_branch <- function(...) {
  scale_fill_manual(values = branch_colors, name = "Branch\narchitecture", ...)
}

scale_color_taxon <- function(...) {
  scale_color_manual(values = taxon_colors, name = "Taxonomic\ngroup", ...)
}

scale_fill_taxon <- function(...) {
  scale_fill_manual(values = taxon_colors, name = "Taxonomic\ngroup", ...)
}

# ============================================================================
# SET GLOBAL THEME
# ============================================================================

# Set as default theme for all ggplot2 plots
theme_set(theme_publication(base_size = base_size))

# Update default geom aesthetics
update_geom_defaults("point", list(size = 2.5, alpha = 0.7))
update_geom_defaults("line", list(linewidth = 0.8))
update_geom_defaults("smooth", list(linewidth = 1, alpha = 0.2))
update_geom_defaults("boxplot", list(outlier.size = 1.5, alpha = 0.7))
update_geom_defaults("bar", list(alpha = 0.8))

# ============================================================================
# FIGURE SIZE STANDARDS (in inches)
# ============================================================================

# MEPS column width: 80 mm (3.15 in), page width: 169 mm (6.65 in)
fig_single_col <- 3.15
fig_double_col <- 6.65
fig_full_page <- 9  # For complex multi-panel figures

# Standard figure dimensions
fig_sizes <- list(
  single = c(width = 3.15, height = 3),
  wide = c(width = 6.65, height = 4),
  square = c(width = 5, height = 5),
  panel_2x2 = c(width = 8, height = 8),
  panel_3x2 = c(width = 10, height = 7),
  full = c(width = 10, height = 8)
)

# ============================================================================
# EXPORT
# ============================================================================

cat("✓ Publication theme loaded\n")
cat("  - Site colors: HAU (orange), MAT (blue), MRB (green)\n")
cat("  - Branch colors: tight (red-orange), wide (blue)\n")
cat("  - Base font size:", base_size, "pt\n")
cat("  - Theme: theme_publication()\n\n")
