#!/usr/bin/env Rscript
# ============================================================================
# run_all_survey_analyses.R - Execute complete Survey analysis pipeline
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
# ============================================================================

cat("\n")
cat("═══════════════════════════════════════════════\n")
cat("    SURVEY DATA ANALYSIS PIPELINE             \n")
cat("    Summer 2019 Field Survey                  \n")
cat("═══════════════════════════════════════════════\n\n")

# Record start time
start_time <- Sys.time()

# Set working directory to project root
library(here)
setwd(here())

# Create log file
log_file <- file.path(here("output/reports"),
                      paste0("survey_pipeline_log_",
                            format(Sys.Date(), "%Y%m%d"), ".txt"))
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

# Function to run script with error handling
run_script_safe <- function(script_path, script_name) {
  cat("\n────────────────────────────────────\n")
  cat("Running:", script_name, "\n")
  cat("────────────────────────────────────\n")

  tryCatch({
    source(script_path)
    cat("✅", script_name, "completed successfully\n")
    return(TRUE)
  }, error = function(e) {
    cat("❌ Error in", script_name, ":\n")
    cat("  ", conditionMessage(e), "\n")
    return(FALSE)
  })
}

# Define analysis scripts in order (aligned with H1-H5 hypotheses)
scripts <- list(
  c("00_load_libraries.R", "Load Libraries"),
  c("01_load_clean_data.R", "Load and Clean Data"),
  c("02_community_composition.R", "Community Composition (H1)"),
  c("03_spatial_patterns.R", "Spatial Pattern Analysis"),
  c("04_diversity_analysis.R", "Diversity Analysis (H1, H4)"),
  c("05_coral_cafi_relationships.R", "Coral-CAFI Relationships (H2, H4)"),
  c("06_network_analysis.R", "Network Analysis (H5)"),
  c("14_local_neighborhood_effects.R", "Neighborhood Effects (H3)")
)

# Run all scripts
results <- list()
for (script_info in scripts) {
  script_file <- script_info[1]
  script_name <- script_info[2]
  script_path <- here("scripts", script_file)

  if (file.exists(script_path)) {
    results[[script_file]] <- run_script_safe(script_path, script_name)
  } else {
    cat("⚠️  Script not found:", script_path, "\n")
    results[[script_file]] <- FALSE
  }
}

# Calculate run time
end_time <- Sys.time()
run_time <- difftime(end_time, start_time, units = "mins")

# ============================================================================
# Generate Summary Report
# ============================================================================

cat("\n")
cat("═══════════════════════════════════════════════\n")
cat("    PIPELINE SUMMARY                          \n")
cat("═══════════════════════════════════════════════\n\n")

# Count successes
n_success <- sum(unlist(results))
n_total <- length(results)

cat("Execution Summary:\n")
cat("  - Scripts run:", n_total, "\n")
cat("  - Successful:", n_success, "\n")
cat("  - Failed:", n_total - n_success, "\n")
cat("  - Run time:", round(run_time, 1), "minutes\n\n")

# List outputs
cat("Generated Outputs:\n")
cat("  - Figures:", length(list.files(here("output/figures"),
                                    recursive = TRUE, pattern = "\\.(png|pdf)$")), "\n")
cat("  - Tables:", length(list.files(here("output/tables"),
                                   pattern = "\\.csv$")), "\n")
cat("  - Objects:", length(list.files(here("output/objects"),
                                    pattern = "\\.rds$")), "\n\n")

# Write summary to log file
sink(log_file)
cat("Survey Analysis Pipeline Log\n")
cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("Run time:", round(run_time, 1), "minutes\n\n")

cat("Script Results:\n")
for (i in seq_along(results)) {
  status <- ifelse(results[[i]], "✅ Success", "❌ Failed")
  cat(sprintf("  %-30s %s\n", names(results)[i], status))
}

cat("\nOutput Summary:\n")
cat("  Figures:", length(list.files(here("output/figures"),
                                    recursive = TRUE, pattern = "\\.(png|pdf)$")), "\n")
cat("  Tables:", length(list.files(here("output/tables"),
                                   pattern = "\\.csv$")), "\n")
cat("  Objects:", length(list.files(here("output/objects"),
                                    pattern = "\\.rds$")), "\n")
sink()

# ============================================================================
# Create Data Summary
# ============================================================================

if (file.exists(here("output/objects/survey_master_data.rds"))) {
  survey_data <- readRDS(here("output/objects/survey_master_data.rds"))

  cat("Data Summary:\n")
  cat("  - Total corals:", nrow(survey_data), "\n")
  cat("  - Sites:", length(unique(survey_data$site)), "\n")
  cat("  - Variables:", ncol(survey_data), "\n\n")
}

# Final message
if (n_success == n_total) {
  cat("✅ SURVEY ANALYSIS PIPELINE COMPLETED SUCCESSFULLY!\n\n")
} else {
  cat("⚠️  Pipeline completed with", n_total - n_success, "errors.\n")
  cat("   Check log file for details:", log_file, "\n\n")
}

cat("📁 All outputs saved to: output/\n")
cat("📊 View results in: output/figures/\n")
cat("📋 Statistical tables in: output/tables/\n")
cat("📝 Log file: ", log_file, "\n\n")

# ============================================================================
# Optional: Generate HTML Summary
# ============================================================================

cat("Generating HTML summary report...\n")

# Create simple HTML summary
html_content <- paste0(
  "<html><head><title>Survey Analysis Summary</title></head><body>",
  "<h1>Survey Data Analysis Summary</h1>",
  "<p>Generated: ", Sys.Date(), "</p>",
  "<h2>Pipeline Results</h2>",
  "<ul>",
  "<li>Scripts run: ", n_total, "</li>",
  "<li>Successful: ", n_success, "</li>",
  "<li>Run time: ", round(run_time, 1), " minutes</li>",
  "</ul>",
  "<h2>Outputs</h2>",
  "<ul>",
  "<li>Figures: ", length(list.files(here("output/figures"),
                                     recursive = TRUE, pattern = "\\.(png|pdf)$")), "</li>",
  "<li>Tables: ", length(list.files(here("output/tables"),
                                    pattern = "\\.csv$")), "</li>",
  "<li>Objects: ", length(list.files(here("output/objects"),
                                     pattern = "\\.rds$")), "</li>",
  "</ul>",
  "</body></html>"
)

html_file <- file.path(here("output/reports"), "survey_analysis_summary.html")
writeLines(html_content, html_file)
cat("📄 HTML summary saved:", html_file, "\n\n")

cat("═══════════════════════════════════════════════\n")
cat("    END OF SURVEY ANALYSIS PIPELINE           \n")
cat("═══════════════════════════════════════════════\n\n")