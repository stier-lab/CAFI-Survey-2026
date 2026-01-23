# ============================================================================
# 00_setup.R - CAFI Survey Analysis: Packages & Configuration
# ============================================================================
#
# PURPOSE: Load required packages, set paths, and define publication theme
#
# USAGE:
#   source("scripts/00_setup.R")
#
# NOTE: Run this script before any analysis scripts
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI SURVEY ANALYSIS - Setup\n")
cat("============================================================\n\n")

# ============================================================================
# PACKAGES
# ============================================================================

# Core packages
library(here)        # Project-relative paths
library(tidyverse)   # Data manipulation & visualization
library(readxl)      # Read Excel files
library(janitor)     # Data cleaning helpers

# Statistical packages
library(MASS)        # Negative binomial GLM (glm.nb)
library(vegan)       # Community ecology (diversity, NMDS, PERMANOVA)
library(lme4)        # Mixed effects models (optional)
library(broom)       # Tidy model outputs
library(car)         # VIF and other diagnostics
library(moments)     # Skewness, kurtosis

# Visualization
library(patchwork)   # Combine ggplots
library(viridis)     # Color scales
library(scales)      # Axis formatting

cat("[OK] Packages loaded\n")

# ============================================================================
# PATHS
# ============================================================================

# Set working directory to project root
setwd(here())

# Define paths (all relative to project root via here())
PATHS <- list(
  # Data inputs
  data = here("data"),
  cafi_data = here("data", "survey_cafi_data_w_taxonomy_summer2019_v5.csv"),
  coral_data = here("data", "survey_coral_characteristics_merged_v2.csv"),
  physio_data = here("data", "survey_master_phys_data_v3.csv"),

  # Base outputs
  output = here("output"),
  objects = here("output", "objects"),
  tables = here("output", "tables"),
  reports = here("output", "reports"),

  # Script-specific figure directories
  figures = here("output", "figures"),
  fig_01_data = here("output", "figures", "01_data"),
  fig_02_community = here("output", "figures", "02_community"),
  fig_03_landscape = here("output", "figures", "03_landscape"),
  fig_04_effects = here("output", "figures", "04_effects"),
  fig_05_scaling = here("output", "figures", "05_scaling"),
  fig_manuscript = here("output", "figures", "manuscript")
)

# Create all output directories
output_dirs <- c("objects", "tables", "reports", "figures",
                 "fig_01_data", "fig_02_community", "fig_03_landscape",
                 "fig_04_effects", "fig_05_scaling", "fig_manuscript")
walk(PATHS[output_dirs], ~dir.create(.x, recursive = TRUE, showWarnings = FALSE))

cat("[OK] Paths configured\n")

# ============================================================================
# PUBLICATION THEME
# ============================================================================

# Define consistent theme for all figures
theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      # Panel
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.border = element_rect(color = "black", linewidth = 0.5),

      # Axis
      axis.text = element_text(color = "black", size = base_size - 2),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),

      # Legend
      legend.position = "right",
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),

      # Strip (for facets)
      strip.background = element_rect(fill = "gray95", color = "black"),
      strip.text = element_text(size = base_size - 1, face = "bold"),

      # Plot
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      plot.caption = element_text(size = base_size - 3, hjust = 1, color = "gray50")
    )
}

# Set as default
theme_set(theme_publication())

# Site color palette (consistent across all figures)
SITE_COLORS <- c(
  "HAU" = "#E69F00",  # Orange - Hauru (fringing)
  "MAT" = "#56B4E9",  # Blue - Maatea (back-reef)
  "MRB" = "#009E73"   # Green - Moorea Barrier Reef
)

# Functional group colors (matches functional_group column in data)
FUNC_COLORS <- c(
  "Trapezia" = "#D55E00",        # Red-orange (mutualist crabs)
  "Resident Fish" = "#0072B2",   # Blue
  "Gastropod" = "#CC79A7",  # Pink
  "Other Crab" = "#009E73",      # Teal
  "Shrimp" = "#E69F00",          # Orange
  "Other" = "#999999"            # Gray
)

cat("[OK] Publication theme set\n")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Save RDS object with standard naming
save_object <- function(obj, name) {
  path <- file.path(PATHS$objects, paste0(name, ".rds"))
  saveRDS(obj, path)
  cat("  Saved:", path, "\n")
  invisible(path)
}

# Load RDS object
load_object <- function(name) {
  path <- file.path(PATHS$objects, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Object not found: ", path)
  }
  readRDS(path)
}

# Save figure with standard settings
# script_dir: one of "01_data", "02_community", "03_landscape", "04_effects", "05_scaling", "manuscript"
save_figure <- function(plot, name, width = 8, height = 6, dpi = 300,
                       script_dir = NULL, manuscript = FALSE) {
  if (manuscript) {
    fig_dir <- PATHS$fig_manuscript
  } else if (!is.null(script_dir)) {
    dir_key <- paste0("fig_", script_dir)
    if (dir_key %in% names(PATHS)) {
      fig_dir <- PATHS[[dir_key]]
    } else {
      # Create custom directory if needed
      fig_dir <- file.path(PATHS$figures, script_dir)
      dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    }
  } else {
    fig_dir <- PATHS$figures
  }
  path <- file.path(fig_dir, paste0(name, ".png"))
  ggsave(path, plot, width = width, height = height, dpi = dpi, bg = "white")
  cat("  Saved:", path, "\n")
  invisible(path)
}

# Save table as CSV
save_table <- function(df, name) {
  path <- file.path(PATHS$tables, paste0(name, ".csv"))
  write_csv(df, path)
  cat("  Saved:", path, "\n")
  invisible(path)
}

# Calculate pseudo-R² (McFadden's)
calc_pseudo_r2 <- function(model, null_model = NULL) {
  if (is.null(null_model)) {
    # Extract from model if possible
    null_dev <- model$null.deviance
    res_dev <- model$deviance
    1 - (res_dev / null_dev)
  } else {
    1 - (logLik(model)[1] / logLik(null_model)[1])
  }
}

cat("[OK] Helper functions defined\n")

# ============================================================================
# SESSION INFO
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("Setup complete.\n")
cat("  R version:", R.version.string, "\n")
cat("  Project root:", here(), "\n")
cat("  Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("------------------------------------------------------------\n\n")

cat("Next: source('scripts/01_load_data.R')\n\n")
