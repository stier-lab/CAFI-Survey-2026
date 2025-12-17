#!/usr/bin/env Rscript
# ============================================================================
# 00_setup.R - Setup and Configuration for CAFI Publication Analyses
# ============================================================================
#
# PURPOSE: Load all required libraries, set paths, and define publication theme
#
# This script is sourced by all analysis scripts. It provides:
#   1. Required package loading with graceful fallbacks
#   2. Standardized output paths
#   3. Publication-quality figure theme (MEPS journal specs)
#   4. Helper functions for saving outputs
#   5. Functional group definitions for CAFI analyses
#
# FUNCTIONAL GROUPS (biologically understood):
#   - Trapezia (Mutualist Defenders): protect corals from predators, sediments
#   - Resident Fishes (Nutrient Providers): enhance coral via ammonium excretion
#   - Corallivores (Drupella snails): negative effects - tissue loss, disease
#   - Other Invertebrates: included statistically, functional roles unknown
#
# Author: CAFI Analysis Pipeline (Refactored)
# Date: 2025-12-03
# ============================================================================

cat("\n========================================\n")
cat("CAFI Publication Analysis Setup\n")
cat("========================================\n\n")

# ============================================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================================

cat("Loading required packages...\n")

# IMPORTANT: Load MASS first to prevent it from masking dplyr::select
# MASS masks dplyr::select if loaded after tidyverse
suppressPackageStartupMessages({
  library(MASS)         # Negative binomial GLMs - LOAD FIRST
})

# Statistical packages (that don't mask tidyverse)
suppressPackageStartupMessages({
  library(vegan)        # Community ecology (PERMANOVA, NMDS, diversity)
  library(lme4)         # Mixed-effects models
  library(lmerTest)     # P-values for lmer
  library(emmeans)      # Estimated marginal means
  library(car)          # Type III ANOVA, VIF
  library(broom)        # Tidy model outputs
  library(broom.mixed)  # Tidy mixed model outputs
})

# Core packages - load tidyverse AFTER MASS
suppressPackageStartupMessages({
  library(tidyverse)    # Loaded after MASS so select() is from dplyr
  library(here)
  library(readxl)
  library(janitor)
})

# Visualization packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)    # Combining plots
  library(viridis)      # Colorblind-friendly palettes
  library(scales)       # Axis formatting
  library(ggrepel)      # Text label repulsion
  library(corrplot)     # Correlation matrices
})

# Network analysis
suppressPackageStartupMessages({
  library(igraph)
})

# Spatial (optional)
if (requireNamespace("sf", quietly = TRUE)) {
  suppressPackageStartupMessages(library(sf))
}

# GAMs for nonlinear relationships
if (requireNamespace("mgcv", quietly = TRUE)) {
  suppressPackageStartupMessages(library(mgcv))
}

cat("  Required packages loaded\n")

# ============================================================================
# 2. OUTPUT PATHS
# ============================================================================

cat("Setting up output directories...\n")

# Base paths
BASE_DIR <- here::here()
DATA_DIR <- here("data")
SCRIPTS_DIR <- here("scripts/publication")

# Output structure
OUTPUT_DIR <- here("output")
FIGURES_DIR <- here("output/figures/publication")
TABLES_DIR <- here("output/tables")
OBJECTS_DIR <- here("output/objects")
REPORTS_DIR <- here("output/reports")

# Create directories if needed
for (d in c(FIGURES_DIR, TABLES_DIR, OBJECTS_DIR, REPORTS_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Figure subdirectories for each publication figure
FIGURE_DIRS <- list(
  fig1 = file.path(FIGURES_DIR, "fig1_study_system"),
  fig2 = file.path(FIGURES_DIR, "fig2_landscape_scaling"),
  fig3 = file.path(FIGURES_DIR, "fig3_functional_groups"),
  fig4 = file.path(FIGURES_DIR, "fig4_composition_networks"),
  fig5 = file.path(FIGURES_DIR, "fig5_condition_landscape"),
  fig6 = file.path(FIGURES_DIR, "fig6_cafi_feedbacks"),
  supp = file.path(FIGURES_DIR, "supplementary")
)

for (d in FIGURE_DIRS) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("  Output directories configured\n")

# ============================================================================
# 3. PUBLICATION THEME (MEPS Journal Specs)
# ============================================================================

cat("Setting publication theme...\n")

# Color palettes - colorblind-friendly
SITE_COLORS <- c(
  "HAU" = "#E69F00",  # Orange
  "MAT" = "#56B4E9",  # Blue
  "MRB" = "#009E73"   # Green
)

BRANCH_COLORS <- c(
  "tight" = "#D55E00",  # Red-orange
  "wide" = "#0072B2"    # Blue
)

# Functional group colors (biologically meaningful groups)
FUNCTIONAL_COLORS <- c(
  "Trapezia" = "#CC79A7",      # Pink - Mutualist defenders
  "Resident Fish" = "#0072B2", # Blue - Nutrient providers
  "Corallivore" = "#D55E00",   # Red - Negative effects
  "Other Crab" = "#E69F00",    # Orange
  "Shrimp" = "#F0E442",        # Yellow
  "Other" = "#999999"          # Gray
)

# Taxonomic group colors (for general analyses)
TAXON_COLORS <- c(
  "crab" = "#CC79A7",

"shrimp" = "#F0E442",
  "snail" = "#999999",
  "fish" = "#0072B2"
)

# Publication theme
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

      # Facets
      strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.3),
      strip.text = element_text(size = rel(0.9), face = "plain", margin = margin(3, 3, 3, 3)),

      # Plot
      plot.title = element_text(size = rel(1.1), face = "bold", hjust = 0, margin = margin(b = 8)),
      plot.subtitle = element_text(size = rel(0.9), hjust = 0, margin = margin(b = 8), color = "gray30"),
      plot.caption = element_text(size = rel(0.7), hjust = 1, color = "gray50"),
      plot.margin = margin(10, 10, 10, 10),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# Set as default
theme_set(theme_publication(base_size = 12))

# Update geom defaults
update_geom_defaults("point", list(size = 2.5, alpha = 0.7))
update_geom_defaults("line", list(linewidth = 0.8))
update_geom_defaults("smooth", list(linewidth = 1, alpha = 0.2))

# Figure size standards (MEPS)
FIG_SIZES <- list(
  single_col = c(width = 3.15, height = 3),      # 80mm
  double_col = c(width = 6.65, height = 4),      # 169mm
  full_page = c(width = 6.65, height = 9),
  panel_2x2 = c(width = 8, height = 8),
  panel_3x2 = c(width = 10, height = 7)
)

cat("  Publication theme set\n")

# ============================================================================
# 4. HELPER FUNCTIONS
# ============================================================================

cat("Defining helper functions...\n")

# Save publication figure
save_pub_fig <- function(plot, filename, dir = FIGURES_DIR,
                         width = 10, height = 8, dpi = 300, ...) {
  full_path <- file.path(dir, filename)
  ggsave(full_path, plot, width = width, height = height, dpi = dpi,
         bg = "white", ...)
  cat("  Saved:", full_path, "\n")
  invisible(full_path)
}

# Save table as CSV
save_table <- function(df, filename, dir = TABLES_DIR) {
  full_path <- file.path(dir, filename)
  write_csv(df, full_path)
  cat("  Saved:", full_path, "\n")
  invisible(full_path)
}

# Save R object
save_object <- function(obj, filename, dir = OBJECTS_DIR) {
  full_path <- file.path(dir, filename)
  saveRDS(obj, full_path)
  cat("  Saved:", full_path, "\n")
  invisible(full_path)
}

# Load R object
load_object <- function(filename, dir = OBJECTS_DIR) {
  full_path <- file.path(dir, filename)
  if (!file.exists(full_path)) {
    stop("Object file not found: ", full_path)
  }
  readRDS(full_path)
}

# Scale functions for consistent styling
scale_color_site <- function(...) {
  scale_color_manual(values = SITE_COLORS, name = "Site", ...)
}

scale_fill_site <- function(...) {
  scale_fill_manual(values = SITE_COLORS, name = "Site", ...)
}

scale_color_functional <- function(...) {
  scale_color_manual(values = FUNCTIONAL_COLORS, name = "Functional\nGroup", ...)
}

scale_fill_functional <- function(...) {
  scale_fill_manual(values = FUNCTIONAL_COLORS, name = "Functional\nGroup", ...)
}

scale_color_taxon <- function(...) {
  scale_color_manual(values = TAXON_COLORS, name = "Type", ...)
}

scale_fill_taxon <- function(...) {
  scale_fill_manual(values = TAXON_COLORS, name = "Type", ...)
}

# ============================================================================
# 5. FUNCTIONAL GROUP DEFINITIONS
# ============================================================================

# Define functional groups based on biological understanding
# These are used for Figure 3 and Figure 6 analyses

FUNCTIONAL_GROUPS <- list(
  # Trapezia crabs - mutualist defenders
  # Protect corals from predators (Acanthaster), sediments, and fouling
  trapezia = c("Trapezia", "trapezia", "Trapezia sp", "Trapezia spp"),

  # Resident fishes - nutrient providers
  # Enhance coral condition via ammonium excretion and microflow
  resident_fish = c("fish", "Fish", "Gobiodon", "Paragobiodon"),

  # Corallivores - primarily Drupella snails
  # Negative effects: tissue loss, disease transmission
  corallivore = c("Drupella", "drupella", "Coralliophila", "snail_corallivore"),

  # Other crabs - functional role less clear
  other_crab = c("Tetralia", "Cymo", "Domecia", "Pilumnus", "crab"),

  # Shrimp - various roles
  shrimp = c("shrimp", "Shrimp", "Alpheus", "Synalpheus", "Periclimenes")
)

# Function to assign functional group
assign_functional_group <- function(species_name, type) {
  species_lower <- tolower(species_name)
  type_lower <- tolower(type)

  # Check Trapezia first (most important mutualist)
  if (grepl("trapezia", species_lower)) {
    return("Trapezia")
  }

  # Check fish
  if (type_lower == "fish") {
    return("Resident Fish")
  }

  # Check corallivores (Drupella snails)
  if (grepl("drupella|coralliophila", species_lower)) {
    return("Corallivore")
  }

  # Check other crabs
  if (type_lower == "crab" && !grepl("trapezia", species_lower)) {
    return("Other Crab")
  }

  # Check shrimp
  if (type_lower == "shrimp") {
    return("Shrimp")
  }

  # Default
  return("Other")
}

cat("  Helper functions and group definitions loaded\n")

# ============================================================================
# 6. STATISTICAL HELPER FUNCTIONS
# ============================================================================

# Function to test and report if regression line should be shown
# (only show if relationship is significant at p < 0.10)
should_show_regression <- function(model, term = NULL, alpha = 0.10) {
  tidy_model <- broom::tidy(model)
  if (!is.null(term)) {
    tidy_model <- tidy_model %>% filter(term == !!term)
  } else {
    # Exclude intercept
    tidy_model <- tidy_model %>% filter(term != "(Intercept)")
  }
  any(tidy_model$p.value < alpha, na.rm = TRUE)
}

# Function to format p-values for reporting
format_p <- function(p) {
  if (p < 0.001) return("< 0.001")
  if (p < 0.01) return(sprintf("= %.3f", p))
  if (p < 0.05) return(sprintf("= %.3f", p))
  return(sprintf("= %.2f", p))
}

# Function to extract key statistics from model
extract_model_stats <- function(model, term) {
  tidy_model <- broom::tidy(model, conf.int = TRUE)
  stats <- tidy_model %>% filter(term == !!term)
  if (nrow(stats) == 0) return(NULL)

  list(
    estimate = stats$estimate,
    se = stats$std.error,
    conf.low = stats$conf.low,
    conf.high = stats$conf.high,
    statistic = stats$statistic,
    p.value = stats$p.value,
    significant = stats$p.value < 0.05
  )
}

# ============================================================================
# SETUP COMPLETE
# ============================================================================

cat("\n========================================\n")
cat("Setup Complete\n")
cat("========================================\n\n")
cat("Paths configured:\n")
cat("  Data:", DATA_DIR, "\n")
cat("  Figures:", FIGURES_DIR, "\n")
cat("  Tables:", TABLES_DIR, "\n")
cat("  Objects:", OBJECTS_DIR, "\n\n")
cat("Ready for publication analyses!\n\n")
