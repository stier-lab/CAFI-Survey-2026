#!/usr/bin/env Rscript
# ============================================================================
# 00_setup.R - Setup and Configuration
# ============================================================================
#
# PURPOSE: Load all packages, configure paths, define publication theme,
#          and provide helper functions for the CAFI Survey analysis.
#
# USAGE: This script is sourced by all other scripts in the pipeline.
#        source(here::here("scripts/00_setup.R"))
#
# CONTENTS:
#   1. Package Loading
#   2. Output Path Configuration
#   3. Color Palettes (colorblind-friendly)
#   4. Publication Theme (MEPS journal specs)
#   5. Scale Functions
#   6. Helper Functions (save, stats, functional groups)
#   7. Figure Size Standards
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-04
# ============================================================================

cat("\n========================================\n")
cat("CAFI Survey Analysis Setup\n")
cat("========================================\n\n")

# ============================================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================================

cat("Loading required packages...\n")

# IMPORTANT: Load MASS first to prevent it from masking dplyr::select
suppressPackageStartupMessages({
  library(MASS)
})

# Statistical packages
suppressPackageStartupMessages({
  library(vegan)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(car)
  library(broom)
  library(broom.mixed)
})

# Core packages - load tidyverse AFTER MASS
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(readxl)
  library(janitor)
})

# Visualization packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(scales)
  library(ggrepel)
  library(corrplot)
})

# Network analysis
suppressPackageStartupMessages({
  library(igraph)
})

# Optional packages
if (requireNamespace("indicspecies", quietly = TRUE)) {
  suppressPackageStartupMessages(library(indicspecies))
}
if (requireNamespace("ape", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ape))
}
if (requireNamespace("moments", quietly = TRUE)) {
  suppressPackageStartupMessages(library(moments))
}

cat("  Packages loaded successfully\n")

# ============================================================================
# 2. CONFIGURE OUTPUT PATHS
# ============================================================================

cat("Setting up output directories...\n")

# Base directories
DATA_DIR        <- here("data")
OUTPUT_DIR      <- here("output")
FIGURES_DIR     <- here("output/figures")
TABLES_DIR      <- here("output/tables")
OBJECTS_DIR     <- here("output/objects")
REPORTS_DIR     <- here("output/reports")
MANUSCRIPT_DIR  <- here("output/figures/manuscript")

# Manuscript figure output directory (all publication figures go here)
# Figures are named: Fig1_*.png, Fig2_*.png, etc.

# Create all directories
all_dirs <- c(DATA_DIR, OUTPUT_DIR, FIGURES_DIR, TABLES_DIR, OBJECTS_DIR,
              REPORTS_DIR, MANUSCRIPT_DIR)
for (d in all_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("  Output directories configured\n")

# ============================================================================
# 3. COLOR PALETTES (Colorblind-Friendly Okabe-Ito)
# ============================================================================

# Site colors
site_colors <- c(
  "HAU" = "#E69F00",
  "MAT" = "#56B4E9",
  "MRB" = "#009E73"
)

# Branch architecture colors
branch_colors <- c(
  "tight" = "#D55E00",
  "wide"  = "#0072B2"
)

# Taxonomic group colors (lowercase and title case for flexibility)
taxon_colors <- c(
  # Lowercase versions
  "crab"        = "#D55E00",     # vermillion
  "shrimp"      = "#0072B2",     # blue
  "fish"        = "#009E73",     # bluish green
  "snail"       = "#CC79A7",     # reddish purple
  "hermit"      = "#E69F00",     # orange
  "echinoderm"  = "#56B4E9",     # sky blue
  "worm"        = "#F0E442",     # yellow
  "amphipod"    = "#F0E442",     # yellow
  "amph"        = "#F0E442",     # yellow
  "squat_lobster" = "#CC79A7",   # reddish purple
  "other"       = "#999999",     # gray
  # Title case versions
  "Crab"        = "#D55E00",
  "Shrimp"      = "#0072B2",
  "Fish"        = "#009E73",
  "Snail"       = "#CC79A7",
  "Hermit"      = "#E69F00",
  "Echinoderm"  = "#56B4E9",
  "Worm"        = "#F0E442",
  "Amphipod"    = "#F0E442",
  "Other"       = "#999999"
)

# Functional group colors
functional_colors <- c(
  "Trapezia"      = "#CC79A7",
  "Resident Fish" = "#0072B2",
  "Corallivore"   = "#D55E00",
  "Other Crab"    = "#F0E442",
  "Shrimp"        = "#56B4E9",
  "Other"         = "#999999"
)

# Shape palettes for redundant encoding (accessibility)
site_shapes <- c("HAU" = 16, "MAT" = 17, "MRB" = 15)
taxon_shapes <- c("crab" = 16, "shrimp" = 17, "snail" = 15, "fish" = 18)
functional_shapes <- c(
  "Trapezia" = 16, "Resident Fish" = 17, "Corallivore" = 18,
  "Other Crab" = 15, "Shrimp" = 3, "Other" = 4
)

# ============================================================================
# 4. PUBLICATION THEME (MEPS Journal Specs)
# ============================================================================

cat("Setting publication theme...\n")

theme_publication <- function(base_size = 12, base_family = "sans") {
  # Font size minimums for publication: axis titles >= 12pt, axis text >= 10pt
  # base_size = 12 ensures: titles = 12pt, text = 10.8pt (0.9 * 12)
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Panel
      panel.background = element_rect(fill = "white", color = NA),
      panel.border     = element_rect(fill = NA, color = "black", linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Axes - IMPORTANT: y-axis title vertical (angle = 90)
      axis.line        = element_blank(),
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text        = element_text(size = rel(0.9), color = "black"),  # ~10.8pt
      axis.title       = element_text(size = rel(1), face = "plain"),     # 12pt
      axis.title.x     = element_text(margin = ggplot2::margin(t = 8)),
      axis.title.y     = element_text(angle = 90, vjust = 0.5, margin = ggplot2::margin(r = 8)),

      # Legend - minimum 10pt text
      legend.background = element_rect(fill = "white", color = NA),
      legend.key        = element_rect(fill = "white", color = NA),
      legend.key.size   = unit(0.9, "lines"),
      legend.text       = element_text(size = rel(0.85)),   # ~10.2pt
      legend.title      = element_text(size = rel(0.9), face = "plain"),
      legend.position   = "right",
      legend.margin     = ggplot2::margin(2, 2, 2, 2),
      legend.spacing    = unit(0.3, "cm"),

      # Facets
      strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.3),
      strip.text       = element_text(size = rel(0.9), face = "plain", margin = ggplot2::margin(4, 4, 4, 4)),

      # Plot
      plot.title      = element_text(size = rel(1.1), face = "bold", hjust = 0, margin = ggplot2::margin(b = 10)),
      plot.subtitle   = element_text(size = rel(0.9), hjust = 0, margin = ggplot2::margin(b = 8), color = "gray30"),
      plot.caption    = element_text(size = rel(0.75), hjust = 1, color = "gray50"),
      plot.margin     = ggplot2::margin(12, 12, 12, 12),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# Set as default theme
theme_set(theme_publication(base_size = 12))

# Update default geom aesthetics
update_geom_defaults("point", list(size = 2.5, alpha = 0.7))
update_geom_defaults("line", list(linewidth = 0.8))
update_geom_defaults("smooth", list(linewidth = 1, alpha = 0.2))
update_geom_defaults("boxplot", list(outlier.size = 1.5, alpha = 0.7))

cat("  Publication theme set\n")

# ============================================================================
# 5. SCALE FUNCTIONS FOR CONSISTENT STYLING
# ============================================================================

scale_color_site <- function(...) {
  scale_color_manual(values = site_colors, name = "Site", ...)
}

scale_fill_site <- function(...) {
  scale_fill_manual(values = site_colors, name = "Site", ...)
}

scale_color_branch <- function(...) {
  scale_color_manual(values = branch_colors, name = "Branch\nArchitecture", ...)
}

scale_fill_branch <- function(...) {
  scale_fill_manual(values = branch_colors, name = "Branch\nArchitecture", ...)
}

scale_color_taxon <- function(...) {
  scale_color_manual(values = taxon_colors, name = "Taxonomic\nGroup", ...)
}

scale_fill_taxon <- function(...) {
  scale_fill_manual(values = taxon_colors, name = "Taxonomic\nGroup", ...)
}

scale_color_functional <- function(...) {
  scale_color_manual(values = functional_colors, name = "Functional\nGroup", ...)
}

scale_fill_functional <- function(...) {
  scale_fill_manual(values = functional_colors, name = "Functional\nGroup", ...)
}

# Shape scales for redundant encoding
scale_shape_site <- function(...) {
  scale_shape_manual(values = site_shapes, name = "Site", ...)
}

scale_shape_taxon <- function(...) {
  scale_shape_manual(values = taxon_shapes, name = "Taxonomic\nGroup", ...)
}

scale_shape_functional <- function(...) {
  scale_shape_manual(values = functional_shapes, name = "Functional\nGroup", ...)
}

# ============================================================================
# 6. HELPER FUNCTIONS
# ============================================================================

cat("Defining helper functions...\n")

# --- File I/O Functions ---

#' Save publication-quality figure
#' @param plot ggplot object
#' @param filename filename (e.g., "Fig2_scaling.png")
#' @param dir output directory (default: MANUSCRIPT_DIR)
#' @param width figure width in inches
#' @param height figure height in inches
#' @param dpi resolution (default: 300)
save_pub_fig <- function(plot, filename, dir = MANUSCRIPT_DIR,
                         width = 10, height = 8, dpi = 300) {
  full_path <- file.path(dir, filename)
  ggsave(full_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  cat("  [MANUSCRIPT FIGURE] Saved:", full_path, "\n")
}

#' Save exploratory figure (not for publication)
save_fig <- function(plot, filename, width = 10, height = 8, dpi = 300) {
  full_path <- file.path(FIGURES_DIR, filename)
  ggsave(full_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  cat("  Saved:", full_path, "\n")
}

#' Save table to CSV
save_table <- function(df, filename) {
  full_path <- file.path(TABLES_DIR, filename)
  write_csv(df, full_path)
  cat("  Saved:", full_path, "\n")
}

#' Save R object to RDS
save_object <- function(obj, filename) {
  full_path <- file.path(OBJECTS_DIR, filename)
  saveRDS(obj, full_path)
  cat("  Saved:", full_path, "\n")
}

#' Load R object from RDS
load_object <- function(filename) {
  full_path <- file.path(OBJECTS_DIR, filename)
  if (!file.exists(full_path)) {
    stop("Object file not found: ", full_path)
  }
  readRDS(full_path)
}

# --- Statistical Helper Functions ---

#' Format p-values for display
format_p <- function(p) {
  case_when(
    is.na(p)    ~ "NA",
    p < 0.001   ~ "< 0.001",
    p < 0.01    ~ sprintf("= %.3f", p),
    p < 0.05    ~ sprintf("= %.2f", p),
    TRUE        ~ sprintf("= %.2f", p)
  )
}

#' Extract model statistics for a specific coefficient
extract_model_stats <- function(model, term) {
  tryCatch({
    tidy_model <- broom::tidy(model, conf.int = TRUE)
    coef_row <- tidy_model %>% filter(term == !!term)

    if (nrow(coef_row) == 0) return(NULL)

    list(
      estimate   = coef_row$estimate[1],
      se         = coef_row$std.error[1],
      conf.low   = coef_row$conf.low[1],
      conf.high  = coef_row$conf.high[1],
      p.value    = coef_row$p.value[1],
      significant = coef_row$p.value[1] < 0.05
    )
  }, error = function(e) NULL)
}

#' Create a standardized statistical results row
#' @param hypothesis Short hypothesis identifier (e.g., "H1", "H2a")
#' @param question Research question being addressed
#' @param test_name Name of statistical test
#' @param test_statistic Test statistic value (F, t, chi-sq, etc.)
#' @param df Degrees of freedom (can be string like "2, 108")
#' @param p_value P-value
#' @param effect_size Effect size value (R², Cohen's d, etc.)
#' @param effect_type Type of effect size (e.g., "R²", "Cohen's d", "η²")
#' @param n Sample size
#' @param notes Additional notes or interpretation
create_result_row <- function(hypothesis, question, test_name,
                              test_statistic, df, p_value,
                              effect_size, effect_type, n, notes = "") {
  # Handle NA values safely
  stat_rounded <- if (is.na(test_statistic)) NA_real_ else round(test_statistic, 3)
  effect_rounded <- if (is.na(effect_size)) NA_real_ else round(effect_size, 4)
  sig_value <- if (is.na(p_value)) NA_character_ else ifelse(p_value < 0.05, "Yes", "No")

  data.frame(
    hypothesis = hypothesis,
    question = question,
    test = test_name,
    statistic = stat_rounded,
    df = as.character(df),
    p_value = p_value,
    p_formatted = format_p(p_value),
    effect_size = effect_rounded,
    effect_type = effect_type,
    n = n,
    significant = sig_value,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

#' Initialize an empty results data frame with standard columns
init_results_df <- function() {
  data.frame(
    hypothesis = character(),
    question = character(),
    test = character(),
    statistic = numeric(),
    df = character(),
    p_value = numeric(),
    p_formatted = character(),
    effect_size = numeric(),
    effect_type = character(),
    n = integer(),
    significant = character(),
    notes = character(),
    stringsAsFactors = FALSE
  )
}

#' Save statistical results summary to CSV
#' @param results_df Data frame with statistical results
#' @param script_name Name of the script (e.g., "03_diversity")
#' @param script_title Human-readable title for the analysis
save_stats_summary <- function(results_df, script_name, script_title) {
  # Add metadata
  results_df$script <- script_name
  results_df$analysis <- script_title
  results_df$timestamp <- Sys.time()

  # Save to tables directory
  filename <- paste0(script_name, "_statistical_results.csv")
  full_path <- file.path(TABLES_DIR, filename)
  write_csv(results_df, full_path)

  cat("\n========================================\n")
  cat("STATISTICAL RESULTS SUMMARY\n")
  cat("========================================\n")
  cat("Script:", script_name, "\n")
  cat("Analysis:", script_title, "\n")
  cat("Tests performed:", nrow(results_df), "\n")
  cat("Significant results:", sum(results_df$significant == "Yes", na.rm = TRUE),
      "of", sum(!is.na(results_df$p_value)), "\n")
  cat("Saved to:", full_path, "\n")
  cat("========================================\n\n")

  invisible(results_df)
}

#' Calculate Cohen's d effect size
cohens_d <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  mx <- mean(x, na.rm = TRUE)
  my <- mean(y, na.rm = TRUE)
  sx <- sd(x, na.rm = TRUE)
  sy <- sd(y, na.rm = TRUE)

  # Pooled standard deviation
  sp <- sqrt(((nx - 1) * sx^2 + (ny - 1) * sy^2) / (nx + ny - 2))
  (mx - my) / sp
}

#' Calculate eta-squared from ANOVA
eta_squared <- function(aov_result) {
  ss <- summary(aov_result)[[1]]["Sum Sq"]
  ss_total <- sum(ss)
  ss_effect <- ss[1, 1]
  ss_effect / ss_total
}

#' Calculate partial eta-squared from ANOVA with multiple factors
partial_eta_squared <- function(aov_result, term_index = 1) {
  ss <- summary(aov_result)[[1]]["Sum Sq"]
  ss_effect <- ss[term_index, 1]
  ss_residual <- ss[nrow(ss), 1]
  ss_effect / (ss_effect + ss_residual)
}

# --- Functional Group Assignment ---

#' Assign CAFI to functional groups based on taxonomy
assign_functional_group <- function(otu, type) {
  otu_lower <- tolower(otu)
  type_lower <- tolower(type)

  case_when(
    str_detect(otu_lower, "trapezia|trapez") ~ "Trapezia",
    str_detect(otu_lower, "drupella")        ~ "Corallivore",
    type_lower == "fish"                     ~ "Resident Fish",
    type_lower == "crab"                     ~ "Other Crab",
    type_lower == "shrimp"                   ~ "Shrimp",
    TRUE                                     ~ "Other"
  )
}

# ============================================================================
# 7. FIGURE SIZE STANDARDS (MEPS Journal)
# ============================================================================

# MEPS column width: 80 mm (3.15 in), page width: 169 mm (6.65 in)
fig_sizes <- list(
  single     = c(width = 3.15, height = 3),
  wide       = c(width = 6.65, height = 4),
  square     = c(width = 5, height = 5),
  panel_2x2  = c(width = 8, height = 8),
  panel_3x2  = c(width = 10, height = 7),
  full       = c(width = 10, height = 8)
)

# ============================================================================
# SETUP COMPLETE
# ============================================================================

cat("\n========================================\n")
cat("Setup Complete\n")
cat("========================================\n\n")

cat("Paths configured:\n")
cat("  Data:", DATA_DIR, "\n")
cat("  Figures:", FIGURES_DIR, "\n")
cat("  Manuscript Figures:", MANUSCRIPT_DIR, "\n")
cat("  Tables:", TABLES_DIR, "\n")
cat("  Objects:", OBJECTS_DIR, "\n\n")

cat("Theme: theme_publication() (colorblind-friendly)\n")
cat("Sites: HAU (orange), MAT (blue), MRB (green)\n\n")

cat("Ready for analysis!\n\n")
