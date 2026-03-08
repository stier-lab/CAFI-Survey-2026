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
# Last Updated: 2026-03-04
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI SURVEY ANALYSIS - Setup\n")
cat("============================================================\n\n")

# ============================================================================
# PACKAGES
# ============================================================================

# Pre-flight check: verify all required packages are installed
required_pkgs <- c(
  # Core
  "here", "tidyverse", "readxl", "janitor", "conflicted",
  # Statistical
  "vegan", "lme4", "MuMIn", "broom", "car", "moments",
  # Visualization
  "patchwork", "viridis", "scales",
  # Used by downstream scripts (04, 06, 07, 09)
  "igraph", "spdep", "sf", "sandwich", "lmtest", "DHARMa", "ggrepel", "cowplot",
  # Conditional but important (01, 05, 07, 09)
  "magick", "geosphere", "piecewiseSEM", "mediation"
)
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with: install.packages(c(\"", paste(missing_pkgs, collapse = "\", \""), "\"))")
}

# Core packages
library(here)        # Project-relative paths
library(tidyverse)   # Data manipulation & visualization
library(readxl)      # Read Excel files
library(janitor)     # Data cleaning helpers
library(conflicted)  # Detect & resolve namespace conflicts

# Statistical packages
# MASS loaded on-demand via MASS::glm.nb() and MASS::theta.ml()
library(vegan)       # Community ecology (diversity, NMDS, PERMANOVA)
library(lme4)        # Mixed effects models (optional)
library(MuMIn)       # R² for mixed models (r.squaredGLMM)
library(broom)       # Tidy model outputs
library(car)         # VIF and other diagnostics
library(moments)     # Skewness, kurtosis

# Visualization
library(patchwork)   # Combine ggplots
library(viridis)     # Color scales
library(scales)      # Axis formatting

# Ensure dplyr::select() always wins (safety net for any future MASS load)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("lag", "dplyr")
conflicted::conflict_prefer("recode", "dplyr")
conflicted::conflict_prefer("some", "purrr")
conflicted::conflict_prefer("discard", "purrr")
conflicted::conflict_prefer("expand", "tidyr")
conflicted::conflict_prefer("pack", "tidyr")
conflicted::conflict_prefer("unpack", "tidyr")
conflicted::conflict_prefer("chisq.test", "stats")
conflicted::conflict_prefer("fisher.test", "stats")
conflicted::conflict_prefer("col_factor", "scales")

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
  fig_06_network = here("output", "figures", "06_network"),
  fig_07_spatial = here("output", "figures", "07_spatial"),
  fig_08_functional = here("output", "figures", "functional_groups"),
  fig_supplement = here("output", "figures", "supplement"),
  fig_manuscript = here("output", "figures", "manuscript")
)

# Create all output directories (automatically walks all PATHS entries)
walk(PATHS, ~dir.create(.x, recursive = TRUE, showWarnings = FALSE))

cat("[OK] Paths configured\n")

# ============================================================================
# PUBLICATION THEME
# ============================================================================

# Define consistent theme for all figures
# Aligned with companion experimental paper (CAFI-136) for visual cohesion
# @param base_size  Base font size in pt (default 13; use 9 for embedded sub-panels)
# @param grid       If TRUE, show major gridlines (grey90); default FALSE
theme_publication <- function(base_size = 13, grid = FALSE) {
  theme_bw(base_size = base_size) +
    theme(
      # Panel
      panel.grid.major  = if(grid) element_line(color = "grey90", linewidth = 0.25) else element_blank(),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.background  = element_rect(fill = "white"),

      # Axis
      axis.title        = element_text(size = base_size, face = "bold"),
      axis.text         = element_text(size = base_size - 1, color = "black"),
      axis.ticks        = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),

      # Legend — default "none" per no-legend policy; override explicitly when needed
      legend.position    = "none",
      legend.direction   = "horizontal",
      legend.title       = element_text(size = base_size - 1, face = "bold"),
      legend.text        = element_text(size = base_size - 2),
      legend.key         = element_rect(fill = "white", color = NA),
      legend.background  = element_rect(fill = "white", color = NA),

      # Facet strips
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.8),
      strip.text       = element_text(size = base_size - 1, face = "bold"),

      # Plot-level
      plot.title    = element_text(size = base_size + 3, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      plot.caption  = element_text(size = base_size - 3, hjust = 1, color = "grey40"),
      plot.margin   = margin(0.5, 0.5, 0.5, 0.5, "cm"),
      plot.background = element_rect(fill = "white", color = NA),

      # Panel tags (for patchwork)
      plot.tag = element_text(face = "bold", size = 16)
    )
}

# Compact theme for dense multi-panel composites (e.g., Fig 5: 4-panel)
# @param base_size  Base font size in pt (default 11; smaller than theme_publication)
theme_multipanel <- function(base_size = 11) {
  theme_publication(base_size = base_size) +
    theme(
      plot.title  = element_text(size = base_size + 3),
      axis.title  = element_text(size = base_size),
      axis.text   = element_text(size = base_size - 1),
      strip.text  = element_text(size = base_size),
      legend.position = "none",
      legend.box      = "horizontal"
    )
}

# Set as default
theme_set(theme_publication())

# Site color palette (consistent across all figures)
SITE_COLORS <- c(
  "HAU" = "#9B7EB8",  # Muted purple - Hauru (fringing)
  "MAT" = "#7B9BAE",  # Cool slate - Maatea (back-reef)
  "MRB" = "#7AAC6D"   # Sage green - Moorea Barrier Reef
)
# NOTE: Palette chosen to avoid confusion with scaling-class colors in Fig 2C
# (Redirection = blue #5A8FAF, Super-linear = vermillion #D55E00).
# Previous Okabe-Ito palette (orange/blue/teal) conflicted with those semantics.

# Scaling-class palette (used in Fig 2 Panel C and species forest plots)
SCALING_COLORS <- c(
  "Redirection"     = "#5A8FAF",
  "Field of Dreams" = "gray55",
  "Super-linear"    = "#D55E00"
)

# Functional group colors (used in scripts 02 and 08 for grouped analyses)
FUNC_GROUP_COLORS <- c(
  "Trapezia" = "#D55E00",        # Red-orange (mutualist crabs)
  "Resident Fish" = "#0072B2",   # Blue
  "Gastropod" = "#CC79A7",       # Pink (all snails)
  "Galeropsis" = "#CC79A7",      # Pink (dominant coral-feeding snail)
  "Other Crab" = "#009E73",      # Teal
  "Shrimp" = "#E69F00",          # Orange
  "Other" = "#999999"            # Gray
)

# Broad taxonomic type colors (Fig 3 barchart + NMDS vectors; colorblind-safe)
# Hue families deliberately separated to avoid pink/magenta confusion
TYPE_COLORS <- c(
  "Shrimp" = "#D55E00",
  "Crabs" = "#0072B2",
  "Hermit crabs" = "#117733",
  "Fish" = "#DDCC77",
  "Gastropods" = "#CC79A7",
  "Echinoderms" = "#88CCEE",
  "Polychaetes" = "#999999",
  "Amphipods" = "#332288",
  "Squat lobsters" = "#661100",
  "Other" = "#BBBBBB"
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

# Save figure as PNG + PDF dual-save (aligned with experimental companion paper)
# path: full file path including .png extension; PDF auto-generated alongside
save_figure <- function(plot, path, width = 8, height = 6, units = "in", dpi = 600) {
  ggsave(path, plot, width = width, height = height, units = units, dpi = dpi, bg = "white")
  pdf_path <- sub("\\.[^.]+$", ".pdf", path)
  # Try cairo_pdf first (higher quality); fall back to standard pdf
  tryCatch(
    suppressWarnings(ggsave(pdf_path, plot, width = width, height = height, units = units,
           device = cairo_pdf, bg = "white")),
    error = function(e) {
      cat("  WARNING: cairo_pdf failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!file.exists(pdf_path) || file.info(pdf_path)$size == 0) {
    tryCatch(
      ggsave(pdf_path, plot, width = width, height = height, units = units,
             device = grDevices::pdf),
      error = function(e) {
        cat("  WARNING: PDF fallback also failed:", conditionMessage(e), "\n")
        NULL
      }
    )
  }
  pdf_ok <- file.exists(pdf_path)
  # PL8: Warn on suspiciously small PDF (likely corrupted partial write)
  if (pdf_ok && file.info(pdf_path)$size < 100) {
    warning(paste("PDF file suspiciously small:", pdf_path))
  }
  cat("  Saved:", path, if(pdf_ok) "(+ PDF)" else "(PNG only)", "\n")
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
# Note: pass null_model explicitly when calling from inside functions,
# because update(model, . ~ 1) can fail due to R scoping (the 'data'
# argument in the model call is a symbol that may not resolve in this scope).
calc_pseudo_r2 <- function(model, null_model = NULL) {
  if (is.null(null_model)) {
    null_model <- tryCatch(
      update(model, . ~ 1),
      error = function(e) NULL
    )
  }
  if (is.null(null_model)) {
    # Fallback: use null.deviance (available on all glm objects)
    # This gives 1 - D_model/D_null (deviance-explained, not identical to McFadden's R²)
    if (!is.null(model$null.deviance) && !is.null(model$deviance)) {
      warning("calc_pseudo_r2: Using deviance-explained (logLik failed). ",
              "Value may differ slightly from McFadden's R-squared.")
      return(1 - (model$deviance / model$null.deviance))
    }
    return(NA_real_)
  }
  ll_full <- as.numeric(logLik(model))
  ll_null <- as.numeric(logLik(null_model))
  1 - (ll_full / ll_null)
}

# Compute Cook's D and flag influential observations
flag_influential <- function(model, threshold = NULL) {
  cd <- cooks.distance(model)
  n <- length(cd)
  if (is.null(threshold)) threshold <- 4 / n
  list(
    cooks_d = cd,
    max_d = max(cd, na.rm = TRUE),
    n_influential = sum(cd > threshold, na.rm = TRUE),
    influential_ids = which(cd > threshold)
  )
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
