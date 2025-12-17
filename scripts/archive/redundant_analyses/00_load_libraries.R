#!/usr/bin/env Rscript
# ============================================================================
# 00_load_libraries.R - Load all required libraries for Survey analysis
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
# ============================================================================

cat("Loading Survey analysis libraries...\n")

# Core packages
library(tidyverse)
library(here)
library(readxl)
library(janitor)

# Statistical packages
library(vegan)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(broom)
library(broom.mixed)

# Visualization packages
library(ggplot2)
library(patchwork)
library(viridis)
library(scales)
library(ggpubr)
library(ggrepel)
library(GGally)
library(corrplot)

# Data manipulation
library(reshape2)
library(data.table)

# Spatial analysis
library(sf)
library(sp)
library(geosphere)
library(ape)  # For Moran.I function

# Network analysis
library(igraph)

# Community ecology
library(indicspecies)
# BiodiversityR is optional - load if available
if (requireNamespace("BiodiversityR", quietly = TRUE)) {
  library(BiodiversityR)
} else {
  cat("Note: BiodiversityR package not installed. Some advanced diversity analyses will be skipped.\n")
}
# iNEXT is optional - load if available
if (requireNamespace("iNEXT", quietly = TRUE)) {
  library(iNEXT)
} else {
  cat("Note: iNEXT package not installed. Some rarefaction analyses will be skipped.\n")
}
# cooccur is optional - load if available
if (requireNamespace("cooccur", quietly = TRUE)) {
  library(cooccur)
} else {
  cat("Note: cooccur package not installed. Network analyses will be skipped.\n")
}
# spatialreg is optional - load if available
if (requireNamespace("spatialreg", quietly = TRUE)) {
  library(spatialreg)
} else {
  cat("Note: spatialreg package not installed. Spatial autocorrelation analyses will be skipped.\n")
}

# Reporting
library(knitr)
library(gt)
library(DT)

# Load custom path configuration
source(here("scripts/utils/path_config.R"))

# Load publication theme for consistent figure styling
source(here("scripts/00_publication_theme.R"))

# Set up Survey-specific paths (flat structure in output/)
SURVEY_BASE    <- here("output")
SURVEY_FIGURES <- here("output", "figures")
SURVEY_TABLES  <- here("output", "tables")
SURVEY_OBJECTS <- here("output", "objects")
SURVEY_REPORTS <- here("output", "reports")

# Create directories if they don't exist
survey_dirs <- c(SURVEY_FIGURES, SURVEY_TABLES, SURVEY_OBJECTS, SURVEY_REPORTS)
for (dir in survey_dirs) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

# Note: Theme is set by 00_publication_theme.R (sourced above)
# Do not override theme_set() here

# Helper function for saving figures
save_survey_plot <- function(plot, filename, width = 10, height = 8, dpi = 300, ...) {
  full_path <- file.path(SURVEY_FIGURES, filename)
  ggsave(full_path, plot, width = width, height = height, dpi = dpi, ...)
  cat("Figure saved:", full_path, "\n")
}

cat("✓ Libraries loaded successfully\n")
cat("✓ Survey paths configured\n")
cat("  - Figures:", SURVEY_FIGURES, "\n")
cat("  - Tables:", SURVEY_TABLES, "\n")
cat("  - Objects:", SURVEY_OBJECTS, "\n")
cat("  - Reports:", SURVEY_REPORTS, "\n\n")