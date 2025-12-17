#!/usr/bin/env Rscript
# ============================================================================
# run_all_publication_analyses.R - Master Orchestrator for Publication Pipeline
# ============================================================================
#
# PURPOSE: Execute all publication analysis scripts in correct order to
# generate the complete 6-figure manuscript analysis.
#
# FIGURE PLAN:
#   Figure 1: Study System + CAFI Overview (conceptual - no statistics)
#   Figure 2: Landscape -> CAFI Scaling (H1, H2)
#   Figure 3: Functional/Taxonomic Group Responses
#   Figure 4: Compositional Changes + Co-occurrence Networks
#   Figure 5: Coral Condition vs Landscape (baseline)
#   Figure 6: CAFI -> Coral Condition Feedbacks (key mutualism tests)
#
# SCRIPT EXECUTION ORDER:
#   1. 00_setup.R - Load packages, set paths, define theme
#   2. 01_load_data.R - Load and clean all data, create objects
#   3. 02_fig2_landscape_scaling.R - Size/proximity -> CAFI
#   4. 03_fig3_functional_groups.R - Group-specific responses
#   5. 04_fig4_composition_networks.R - NMDS, PERMANOVA, networks
#   6. 05_fig5_condition_landscape.R - Condition vs landscape
#   7. 06_fig6_cafi_feedbacks.R - CAFI -> condition (key tests)
#
# USAGE:
#   Rscript run_all_publication_analyses.R
#   OR
#   source("scripts/publication/run_all_publication_analyses.R")
#
# Author: CAFI Analysis Pipeline
# Date: 2025-12-03
# ============================================================================

# Record start time
start_time <- Sys.time()

cat("\n")
cat("================================================================\n")
cat("  CAFI PUBLICATION ANALYSIS PIPELINE\n")
cat("  CAFIs Survey Paper - Full 6-Figure Analysis\n")
cat("================================================================\n\n")
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

# Set working directory to project root
if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
setwd(here::here())

cat("Working directory:", getwd(), "\n\n")

# ============================================================================
# DEFINE SCRIPT EXECUTION ORDER
# ============================================================================

scripts <- c(
  # Data preparation (always runs first)
  "01_load_data.R",

  # Figure analyses (in figure order)
  "02_fig2_landscape_scaling.R",
  "03_fig3_functional_groups.R",
  "04_fig4_composition_networks.R",
  "05_fig5_condition_landscape.R",
  "06_fig6_cafi_feedbacks.R"
)

script_dir <- here::here("scripts/publication")

# Track results
results <- list()
errors <- character()

# ============================================================================
# EXECUTE SCRIPTS
# ============================================================================

for (script in scripts) {

  script_path <- file.path(script_dir, script)

  cat("----------------------------------------------------------------\n")
  cat("EXECUTING:", script, "\n")
  cat("----------------------------------------------------------------\n\n")

  script_start <- Sys.time()

  # Execute script with error handling
  result <- tryCatch({
    source(script_path, local = FALSE)
    list(
      script = script,
      status = "SUCCESS",
      duration = difftime(Sys.time(), script_start, units = "secs"),
      error = NULL
    )
  }, error = function(e) {
    list(
      script = script,
      status = "FAILED",
      duration = difftime(Sys.time(), script_start, units = "secs"),
      error = conditionMessage(e)
    )
  })

  results[[script]] <- result

  if (result$status == "FAILED") {
    errors <- c(errors, script)
    cat("\n!! ERROR in", script, ":\n")
    cat("  ", result$error, "\n\n")
  } else {
    cat("\n== Completed", script, "in",
        round(as.numeric(result$duration), 1), "seconds ==\n\n")
  }
}

# ============================================================================
# SUMMARY REPORT
# ============================================================================

end_time <- Sys.time()
total_duration <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("================================================================\n")
cat("  PIPELINE EXECUTION SUMMARY\n")
cat("================================================================\n\n")

cat("Total execution time:", round(as.numeric(total_duration), 1), "minutes\n\n")

# Script-by-script results
cat("Script Results:\n")
cat("---------------\n")
for (script in scripts) {
  r <- results[[script]]
  status_icon <- if (r$status == "SUCCESS") "[OK]" else "[FAILED]"
  cat(sprintf("  %s %s (%.1f sec)\n", status_icon, r$script, as.numeric(r$duration)))
}

cat("\n")

# Output summary
if (length(errors) == 0) {
  cat("STATUS: ALL SCRIPTS COMPLETED SUCCESSFULLY\n\n")

  # Count outputs
  figures_dir <- here::here("output/figures/publication")
  tables_dir <- here::here("output/tables")
  objects_dir <- here::here("output/objects")

  n_figures <- length(list.files(figures_dir, pattern = "\\.png$", recursive = TRUE))
  n_tables <- length(list.files(tables_dir, pattern = "\\.csv$"))
  n_objects <- length(list.files(objects_dir, pattern = "\\.rds$"))

  cat("Outputs Generated:\n")
  cat("------------------\n")
  cat("  Figures:", n_figures, "PNG files\n")
  cat("  Tables:", n_tables, "CSV files\n")
  cat("  R Objects:", n_objects, "RDS files\n\n")

  cat("Output Locations:\n")
  cat("-----------------\n")
  cat("  Figures:", figures_dir, "\n")
  cat("  Tables:", tables_dir, "\n")
  cat("  Objects:", objects_dir, "\n\n")

} else {
  cat("STATUS: SOME SCRIPTS FAILED\n\n")
  cat("Failed scripts:\n")
  for (err_script in errors) {
    cat("  -", err_script, "\n")
    cat("    Error:", results[[err_script]]$error, "\n")
  }
  cat("\n")
}

# ============================================================================
# FIGURE CHECKLIST
# ============================================================================

cat("================================================================\n")
cat("  PUBLICATION FIGURE CHECKLIST\n")
cat("================================================================\n\n")

checklist <- list(
  "Figure 1" = list(
    desc = "Study System + CAFI Overview",
    status = "CONCEPTUAL (no statistics needed)",
    files = NA
  ),
  "Figure 2" = list(
    desc = "Landscape -> CAFI Scaling",
    status = if ("02_fig2_landscape_scaling.R" %in% errors) "FAILED" else "COMPLETE",
    files = c("fig2_landscape_scaling_combined.png")
  ),
  "Figure 3" = list(
    desc = "Functional Group Responses",
    status = if ("03_fig3_functional_groups.R" %in% errors) "FAILED" else "COMPLETE",
    files = c("fig3_functional_groups_combined.png")
  ),
  "Figure 4" = list(
    desc = "Composition + Networks",
    status = if ("04_fig4_composition_networks.R" %in% errors) "FAILED" else "COMPLETE",
    files = c("fig4_composition_partial.png", "fig4d_cooccurrence_network.png")
  ),
  "Figure 5" = list(
    desc = "Coral Condition vs Landscape",
    status = if ("05_fig5_condition_landscape.R" %in% errors) "FAILED" else "COMPLETE",
    files = c("fig5_condition_landscape_combined.png")
  ),
  "Figure 6" = list(
    desc = "CAFI -> Coral Condition Feedbacks",
    status = if ("06_fig6_cafi_feedbacks.R" %in% errors) "FAILED" else "COMPLETE",
    files = c("fig6_cafi_feedbacks_combined.png")
  )
)

for (fig_name in names(checklist)) {
  fig <- checklist[[fig_name]]
  status_icon <- case_when(
    fig$status == "COMPLETE" ~ "[X]",
    fig$status == "CONCEPTUAL (no statistics needed)" ~ "[~]",
    TRUE ~ "[ ]"
  )
  cat(sprintf("%s %s: %s\n", status_icon, fig_name, fig$desc))
  cat(sprintf("    Status: %s\n", fig$status))
}

cat("\n")
cat("Legend: [X] = Complete, [~] = Conceptual, [ ] = Not Complete\n\n")

# ============================================================================
# NEXT STEPS
# ============================================================================

cat("================================================================\n")
cat("  NEXT STEPS\n")
cat("================================================================\n\n")

if (length(errors) == 0) {
  cat("1. Review generated figures in output/figures/publication/\n")
  cat("2. Check statistical results in output/tables/\n")
  cat("3. Create Figure 1 (conceptual/study system diagram)\n")
  cat("4. Compile figures into manuscript layout\n")
  cat("5. Write results section based on statistical outputs\n\n")
} else {
  cat("1. Fix errors in failed scripts\n")
  cat("2. Re-run pipeline: source('scripts/publication/run_all_publication_analyses.R')\n")
  cat("3. Check data files exist in data/ directory\n\n")
}

cat("================================================================\n")
cat("  PIPELINE COMPLETE\n")
cat("================================================================\n\n")

# Return summary invisibly
invisible(list(
  start_time = start_time,
  end_time = end_time,
  duration_mins = as.numeric(total_duration),
  results = results,
  errors = errors,
  n_errors = length(errors)
))
