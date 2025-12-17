#!/usr/bin/env Rscript
# ============================================================================
# run_pipeline.R - CAFI Survey Analysis Pipeline
# ============================================================================
#
# PURPOSE: Execute the complete analysis pipeline for the CAFI manuscript
#
# USAGE:
#   source("scripts/run_pipeline.R")                    # Run all scripts
#   source("scripts/run_pipeline.R"); run_core_only()   # Core only
#
# PIPELINE STRUCTURE:
#   Core Pipeline (00-09):
#     00_setup.R               - Load packages, set paths, define theme
#     01_load_clean_data.R     - Load and clean data, create core objects
#     02_community_composition - Community composition summaries
#     03_diversity_analysis    - Alpha/beta diversity (supports Fig 4)
#     04_scaling_relationships - Size-abundance scaling (Fig 2)
#     05_coral_condition       - Condition score development (Fig 5)
#     06_network_analysis      - Co-occurrence networks (Fig 4)
#     07_neighborhood_effects  - Local neighborhood effects (Fig 6)
#     08_cafi_condition_feedbacks - CAFI -> condition (Fig 6)
#     09_functional_groups     - Functional group analysis (Fig 3)
#
#   Supplementary Pipeline (10-19):
#     10_spatial_patterns      - Geographic distribution & maps
#     11_size_biomass_scaling  - CAFI size structure & biomass
#     12_machine_learning      - Random Forest, variable importance
#     13_spatial_autocorrelation - Moran's I, spatial correlograms
#     14_composition_by_size   - How composition shifts with size
#     15_composition_by_neighborhood - Composition & neighborhood
#     16_summary_figures       - Multi-panel summary figures
#     18_manuscript_figures    - Final publication figures
#     19_figure2_abundance_richness_landscape - Landscape effects
#
# OUTPUTS:
#   output/figures/manuscript/  - Publication figures (Fig2-6)
#   output/figures/             - Exploratory figures
#   output/tables/              - Statistical results (CSV)
#   output/objects/             - R data objects (RDS)
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI SURVEY ANALYSIS PIPELINE\n")
cat("    Mo'orea, French Polynesia - Summer 2019\n")
cat("============================================================\n\n")

# Record start time
start_time <- Sys.time()

# Set working directory to project root
library(here)
setwd(here())

# Create log file
log_dir <- here("output/reports")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("pipeline_log_", format(Sys.Date(), "%Y%m%d"), ".txt"))

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Run script with error handling
run_script_safe <- function(script_path, script_name) {
  cat("\n------------------------------------------------------------\n")
  cat("Running:", script_name, "\n")
  cat("------------------------------------------------------------\n")

  tryCatch({
    source(script_path)
    cat("[OK]", script_name, "completed successfully\n")
    return(TRUE)
  }, error = function(e) {
    cat("[ERROR] in", script_name, ":\n")
    cat("  ", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# ============================================================================
# DEFINE PIPELINE SCRIPTS
# ============================================================================

# Core pipeline (required for manuscript figures)
core_scripts <- list(
  c("00_setup.R",                   "Setup & Configuration"),
  c("01_load_clean_data.R",         "Load & Clean Data"),
  c("02_community_composition.R",   "Community Composition"),
  c("03_diversity_analysis.R",      "Diversity Analysis (Fig 4)"),
  c("04_scaling_relationships.R",   "Scaling Relationships (Fig 2)"),
  c("05_coral_condition.R",         "Coral Condition (Fig 5)"),
  c("06_network_analysis.R",        "Network Analysis (Fig 4)"),
  c("07_neighborhood_effects.R",    "Neighborhood Effects (Fig 6)"),
  c("08_cafi_condition_feedbacks.R","CAFI-Condition Feedbacks (Fig 6)"),
  c("09_functional_groups.R",       "Functional Groups (Fig 3)")
)

# Supplementary analyses (optional)
supplementary_scripts <- list(
  c("10_spatial_patterns.R",        "Spatial Patterns"),
  c("11_size_biomass_scaling.R",    "Size-Biomass Scaling"),
  c("12_machine_learning.R",        "Machine Learning"),
  c("13_spatial_autocorrelation.R", "Spatial Autocorrelation"),
  c("14_composition_by_size.R",     "Composition by Size"),
  c("15_composition_by_neighborhood.R", "Composition by Neighborhood"),
  c("16_summary_figures.R",         "Summary Figures"),
  c("18_manuscript_figures.R",      "Manuscript Figures"),
  c("19_figure2_abundance_richness_landscape.R", "Landscape Effects")
)

# ============================================================================
# RUN PIPELINE
# ============================================================================

run_core_only <- function() {
  cat("\nRunning CORE pipeline only (scripts 00-09)...\n")
  run_scripts(core_scripts)
}

run_all <- function() {
  cat("\nRunning FULL pipeline (all scripts)...\n")
  all_scripts <- c(core_scripts, supplementary_scripts)
  run_scripts(all_scripts)
}

run_scripts <- function(scripts) {
  results <- list()

  for (script_info in scripts) {
    script_file <- script_info[1]
    script_name <- script_info[2]
    script_path <- here("scripts", script_file)

    if (file.exists(script_path)) {
      results[[script_file]] <- run_script_safe(script_path, script_name)
    } else {
      cat("[SKIP] Script not found:", script_file, "\n")
      results[[script_file]] <- NA
    }
  }

  # Generate summary
  generate_summary(results)

  invisible(results)
}

generate_summary <- function(results) {
  end_time <- Sys.time()
  run_time <- difftime(end_time, start_time, units = "mins")

  cat("\n")
  cat("============================================================\n")
  cat("    PIPELINE SUMMARY\n")
  cat("============================================================\n\n")

  # Count results
  n_success <- sum(unlist(results) == TRUE, na.rm = TRUE)
  n_failed <- sum(unlist(results) == FALSE, na.rm = TRUE)
  n_skipped <- sum(is.na(unlist(results)))
  n_total <- length(results)

  cat("Execution Summary:\n")
  cat("  Scripts run:", n_total, "\n")
  cat("  Successful:", n_success, "\n")
  cat("  Failed:", n_failed, "\n")
  cat("  Skipped:", n_skipped, "\n")
  cat("  Run time:", round(run_time, 1), "minutes\n\n")

  # List outputs
  cat("Generated Outputs:\n")
  cat("  Manuscript figures:", length(list.files(here("output/figures/manuscript"),
                                                  pattern = "\\.png$")), "\n")
  cat("  Exploratory figures:", length(list.files(here("output/figures"),
                                                   recursive = TRUE, pattern = "\\.png$")), "\n")
  cat("  Tables:", length(list.files(here("output/tables"), pattern = "\\.csv$")), "\n")
  cat("  Objects:", length(list.files(here("output/objects"), pattern = "\\.rds$")), "\n\n")

  # Write log file
  sink(log_file)
  cat("CAFI Survey Analysis Pipeline Log\n")
  cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
  cat("Run time:", round(run_time, 1), "minutes\n\n")

  cat("Script Results:\n")
  for (i in seq_along(results)) {
    status <- case_when(
      isTRUE(results[[i]]) ~ "[OK]",
      isFALSE(results[[i]]) ~ "[ERROR]",
      TRUE ~ "[SKIP]"
    )
    cat(sprintf("  %-35s %s\n", names(results)[i], status))
  }
  sink()

  # Final message
  if (n_failed == 0 && n_skipped == 0) {
    cat("[SUCCESS] Pipeline completed successfully!\n\n")
  } else if (n_failed > 0) {
    cat("[WARNING] Pipeline completed with", n_failed, "errors.\n")
    cat("   Check log file:", log_file, "\n\n")
  }

  cat("Output locations:\n")
  cat("  Manuscript figures: output/figures/manuscript/\n")
  cat("  All figures:        output/figures/\n")
  cat("  Tables:             output/tables/\n")
  cat("  Log file:           ", log_file, "\n\n")

  cat("============================================================\n")
  cat("    END OF PIPELINE\n")
  cat("============================================================\n\n")
}

# ============================================================================
# AUTO-RUN CORE PIPELINE
# ============================================================================

cat("Starting pipeline...\n")
cat("To run only core: run_core_only()\n")
cat("To run all:       run_all()\n\n")

# Default: run core pipeline
run_core_only()
