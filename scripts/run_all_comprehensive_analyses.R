#!/usr/bin/env Rscript
# ============================================================================
# run_all_comprehensive_analyses.R - Execute all Survey analyses in order
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
# ============================================================================

cat("\n========================================\n")
cat("Running Comprehensive Survey Analysis Pipeline\n")
cat("========================================\n\n")

# Record start time
start_time <- Sys.time()

# Set up error handling
options(error = function() {
  cat("\n❌ Error occurred. Check the error message above.\n")
  traceback()
})

# ============================================================================
# Run Analysis Scripts in Order
# ============================================================================

scripts_to_run <- c(
  "00_load_libraries.R",
  "01_load_clean_data.R",
  "02_community_composition.R",
  "03_spatial_patterns.R",
  "04_diversity_analysis.R",
  "05_coral_cafi_relationships.R",
  "06_network_analysis.R",
  "07_morphotype_habitat_analysis.R",
  "08_size_biomass_scaling.R",
  "08_advanced_statistical_models.R",
  "09_machine_learning_predictions.R",
  "10_neighborhood_arrival_comparison.R",
  "11_spatial_autocorrelation.R",
  "12_visualization_suite.R",
  "13_comprehensive_predictor_analysis.R",
  "14_coral_size_neighbor_effects.R"
)

# Track execution status
execution_log <- data.frame(
  script = scripts_to_run,
  status = NA,
  execution_time = NA,
  stringsAsFactors = FALSE
)

cat("Will run", length(scripts_to_run), "analysis scripts\n\n")
cat("========================================\n\n")

# Run each script
for (i in seq_along(scripts_to_run)) {
  script <- scripts_to_run[i]
  script_path <- here::here("scripts/Survey", script)

  cat(sprintf("\n[%d/%d] Running %s...\n", i, length(scripts_to_run), script))
  cat(rep("-", 50), "\n", sep = "")

  # Check if script exists
  if (!file.exists(script_path)) {
    cat("  ⚠️ Script not found:", script_path, "\n")
    execution_log$status[i] <- "Not Found"
    next
  }

  # Run the script
  script_start <- Sys.time()

  tryCatch({
    source(script_path)
    execution_log$status[i] <- "Success"
    cat("\n  ✓ Script completed successfully\n")
  }, error = function(e) {
    execution_log$status[i] <- "Failed"
    cat("\n  ❌ Script failed with error:\n")
    cat("    ", conditionMessage(e), "\n")
  })

  script_end <- Sys.time()
  execution_log$execution_time[i] <- round(difftime(script_end, script_start,
                                                    units = "secs"), 2)

  cat("  Execution time:", execution_log$execution_time[i], "seconds\n")
}

# ============================================================================
# Generate Comprehensive Report
# ============================================================================

cat("\n========================================\n")
cat("Generating Comprehensive Report\n")
cat("========================================\n\n")

report_file <- here::here("output/Survey/reports/Survey_Comprehensive_Report.Rmd")

if (file.exists(report_file)) {
  cat("Rendering HTML report...\n")

  tryCatch({
    rmarkdown::render(
      report_file,
      output_file = "Survey_Comprehensive_Report.html",
      output_dir = here::here("output/Survey/reports"),
      quiet = TRUE
    )
    cat("✓ Report generated successfully\n")
    cat("  Output: output/Survey/reports/Survey_Comprehensive_Report.html\n")
  }, error = function(e) {
    cat("❌ Report generation failed:\n")
    cat("  ", conditionMessage(e), "\n")
  })
} else {
  cat("⚠️ Report template not found\n")
}

# ============================================================================
# Summary
# ============================================================================

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

cat("\n========================================\n")
cat("Analysis Pipeline Complete\n")
cat("========================================\n\n")

# Print execution summary
cat("Execution Summary:\n")
cat("  Total scripts:", length(scripts_to_run), "\n")
cat("  Successful:", sum(execution_log$status == "Success", na.rm = TRUE), "\n")
cat("  Failed:", sum(execution_log$status == "Failed", na.rm = TRUE), "\n")
cat("  Not found:", sum(execution_log$status == "Not Found", na.rm = TRUE), "\n")
cat("  Total time:", round(total_time, 2), "minutes\n\n")

# Show failed scripts
failed_scripts <- execution_log[execution_log$status == "Failed", ]
if (nrow(failed_scripts) > 0) {
  cat("Failed Scripts:\n")
  for (i in 1:nrow(failed_scripts)) {
    cat("  -", failed_scripts$script[i], "\n")
  }
  cat("\n")
}

# Save execution log
log_file <- here::here("output/Survey/reports",
                      paste0("execution_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write.csv(execution_log, log_file, row.names = FALSE)
cat("Execution log saved to:", log_file, "\n\n")

# Print output locations
cat("Output Locations:\n")
cat("  Figures: output/Survey/figures/\n")
cat("  Tables: output/Survey/tables/\n")
cat("  Objects: output/Survey/objects/\n")
cat("  Reports: output/Survey/reports/\n\n")

cat("✅ All analyses complete!\n")