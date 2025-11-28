#!/usr/bin/env Rscript
# ==============================================================================
# Run Complete Survey Analysis Pipeline
# ==============================================================================

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║                   SURVEY ANALYSIS PIPELINE                     ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

start_time <- Sys.time()
errors <- list()

# Core Survey scripts in order
scripts <- c(
  # Setup
  "00_load_libraries.R",
  "01_load_clean_data.R",

  # Core analyses
  "02_community_composition.R",
  "03_spatial_patterns.R",
  "04_diversity_analysis.R",
  "05_coral_cafi_relationships.R",
  "06_network_analysis.R",
  "07_morphotype_habitat_analysis.R",

  # Advanced analyses
  "08_advanced_statistical_models.R",
  "08_size_biomass_scaling.R",
  "09_machine_learning_predictions.R",
  "10_neighborhood_arrival_comparison.R",
  "11_spatial_autocorrelation.R",

  # Comprehensive analyses
  "13_comprehensive_predictor_analysis.R",
  "14_coral_size_neighbor_effects.R",

  # Guild-specific analyses
  "17_trapezid_guild_analysis.R",

  # Visualization
  "12_visualization_suite.R",
  "15_comprehensive_visual_summary.R",
  "16_generate_all_summary_figures.R"
)

# Run each script
for (script_name in scripts) {
  script_path <- here::here("scripts", "Survey", script_name)

  if (file.exists(script_path)) {
    cat("\n▶ Running:", script_name, "\n")
    cat("─────────────────────────────────────────────────────\n")

    tryCatch({
      source(script_path, echo = FALSE, print.eval = FALSE)
      cat("✓ Completed:", script_name, "\n")
    }, error = function(e) {
      cat("✗ Error in", script_name, ":", conditionMessage(e), "\n")
      errors[[script_name]] <- conditionMessage(e)
    })
  } else {
    cat("⚠ Script not found:", script_name, "\n")
  }
}

# Summary
cat("\n════════════════════════════════════════════════════════════════\n")
cat("SURVEY ANALYSIS COMPLETE\n")
cat("════════════════════════════════════════════════════════════════\n\n")

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat("Total execution time:", round(total_time, 2), "minutes\n\n")

if (length(errors) == 0) {
  cat("✅ All scripts completed successfully!\n\n")
} else {
  cat("⚠️ Errors encountered in", length(errors), "scripts:\n")
  for (script in names(errors)) {
    cat("  -", script, "\n")
  }
  cat("\nReview the errors above.\n\n")
}

# Check output structure
cat("Output Structure:\n")
cat("─────────────────────────────────────────────────────\n")

# Count files in each subdirectory
fig_subdirs <- list.dirs("output/figures", recursive = FALSE, full.names = FALSE)
for (subdir in fig_subdirs) {
  path <- file.path("output/figures", subdir)
  if (dir.exists(path)) {
    n_files <- length(list.files(path))
    if (n_files > 0) {
      cat(sprintf("  ✓ figures/%s: %d files\n", subdir, n_files))
    }
  }
}

total_figs <- length(list.files("output/figures", recursive = TRUE, pattern = "\\.(png|pdf)$"))
total_tables <- length(list.files("output/tables", recursive = TRUE, pattern = "\\.csv$"))

cat(sprintf("\n  Total figures: %d\n", total_figs))
cat(sprintf("  Total tables: %d\n", total_tables))

cat("\n✅ Survey analysis pipeline complete!\n")