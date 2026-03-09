#!/usr/bin/env Rscript
# REPRODUCIBILITY NOTE: Run renv::init() and renv::snapshot() to pin package versions.
# See renv.lock for the locked environment used in analysis.
# ============================================================================
# run_full_pipeline.R - CAFI Survey Complete Analysis Pipeline
# ============================================================================
#
# PURPOSE: Execute the complete CAFI analysis pipeline in correct order with:
#   - Error handling and recovery
#   - Progress logging to output/pipeline.log
#   - Timing for each script
#   - Dependency checking
#   - Skip completed steps option
#   - Summary report at end
#
# USAGE:
#   # From R console:
#   source("scripts/run_full_pipeline.R")
#
#   # With options:
#   run_pipeline(skip_completed = TRUE, verbose = TRUE)
#
#   # From command line:
#   Rscript scripts/run_full_pipeline.R
#
# OPTIONS:
#   skip_completed - Skip scripts that have already produced output (default: FALSE)
#   verbose - Print detailed progress (default: TRUE)
#   stop_on_error - Stop pipeline on first error (default: FALSE)
#
# PIPELINE ORDER:
#   Core Pipeline:
#     00_setup.R                    - Packages, paths, theme
#     00b_data_quality_audit.R      - Data quality assessment
#     01_load_data.R                - Load and clean data
#     02_community_analysis.R       - Community composition
#     03_landscape_characterization.R - Landscape metrics
#     04_landscape_effects.R        - Landscape effects on CAFI
#     05_species_scaling_analysis.R - Species-area scaling
#
#   Extended Analyses:
#     06_cooccurrence_analysis.R    - Null-model co-occurrence
#     07_spatial_autocorrelation.R  - Spatial patterns
#     08_functional_groups.R        - Functional group analysis
#     09_cafi_condition_feedbacks.R - CAFI-condition feedbacks
#
#   Machine Learning (exploratory, not in default pipeline):
#     10_feature_engineering.R
#     11_machine_learning.R
#     12_model_evaluation.R
#
#   Note: Manuscript figures are created by their respective analysis scripts
#     (01 → Fig 1, 05 → Figs 2+3, 02 → Fig 4, 09 → Fig 5, 06 → S11, 08 → S9)
#
# OUTPUTS:
#   - output/pipeline.log          - Detailed execution log
#   - output/tables/pipeline_timing.csv - Script execution times
#   - Console output with progress bars
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-21
# ============================================================================

# ============================================================================
# CONFIGURATION
# ============================================================================

# Define pipeline scripts in execution order
PIPELINE_SCRIPTS <- list(
  # Core Pipeline
  core = list(
    list(name = "00_setup.R", desc = "Setup: Packages, paths, theme", required = TRUE),
    list(name = "00b_data_quality_audit.R", desc = "Data quality audit", required = FALSE),
    list(name = "01_load_data.R", desc = "Load and clean data", required = TRUE),
    list(name = "02_community_analysis.R", desc = "Community composition analysis", required = TRUE),
    list(name = "03_landscape_characterization.R", desc = "Landscape characterization", required = TRUE),
    list(name = "04_landscape_effects.R", desc = "Landscape effects on CAFI", required = TRUE),
    list(name = "05_species_scaling_analysis.R", desc = "Species-area scaling", required = TRUE)
  ),

  # Extended Analyses
  extended = list(
    list(name = "06_cooccurrence_analysis.R", desc = "Co-occurrence null model analysis", required = FALSE),
    list(name = "07_spatial_autocorrelation.R", desc = "Spatial autocorrelation", required = FALSE),
    list(name = "08_functional_groups.R", desc = "Functional group analysis", required = FALSE),
    list(name = "09_cafi_condition_feedbacks.R", desc = "CAFI-condition feedbacks", required = FALSE),
    list(name = "13_taxonomy_sensitivity.R", desc = "Taxonomy sensitivity analysis", required = FALSE)
  ),

  # Machine Learning (exploratory - not in default pipeline)
  ml = list(
    list(name = "archive/10_feature_engineering.R", desc = "Feature engineering", required = FALSE),
    list(name = "archive/11_machine_learning.R", desc = "Machine learning models", required = FALSE),
    list(name = "archive/12_model_evaluation.R", desc = "Model evaluation", required = FALSE)
  )
)

# Expected outputs for dependency checking
EXPECTED_OUTPUTS <- list(
  "00_setup.R" = character(0),  # Creates global variables only
  # PL6: Added 00b entry

  "00b_data_quality_audit.R" = c(
    "output/tables/data_quality_report.csv"
  ),
  "01_load_data.R" = c(
    "output/objects/coral_master.rds",
    "output/objects/cafi_clean.rds",
    "output/objects/community_matrix.rds",
    "output/objects/condition_scores.rds",
    "output/objects/cafi_pca_results.rds",
    "output/objects/cafi_by_coral.rds",
    "output/objects/otu_taxonomy.rds",
    "output/objects/taxonomy_scenario_data.rds",
    "output/objects/functional_summary.rds",
    "output/objects/overlap_genera.rds",        # PL7: was missing
    "output/figures/manuscript/fig1_study_design.png"
  ),
  "02_community_analysis.R" = c(
    "output/objects/community_analysis_results.rds",
    "output/figures/02_community/community_summary.png",
    "output/figures/manuscript/fig4_composition.png",
    "output/figures/supplement/figS1_species_accumulation.png",
    "output/figures/supplement/figS2_permanova_sensitivity.png",  # PL2
    "output/figures/supplement/figS3_nmds_ordination.png",        # PL2
    "output/figures/supplement/figS5_composition_divergence.png", # PL2
    "output/tables/permanova_subsampling_summary.csv"
  ),
  "03_landscape_characterization.R" = c(
    "output/objects/landscape_selected_predictors.rds",
    "output/tables/landscape_structure_summary.csv"
  ),
  "04_landscape_effects.R" = c(
    "output/tables/landscape_full_model_results.csv",
    "output/figures/supplement/figS7_neighborhood_null.png"       # PL2
  ),
  "05_species_scaling_analysis.R" = c(
    "output/objects/scaling_analysis_results.rds",
    "output/tables/scaling_results_all.csv",
    "output/tables/occurrence_scaling_results.csv",
    "output/figures/manuscript/fig2_scaling.png",
    "output/figures/manuscript/fig3_species_group_scaling.png",
    "output/figures/supplement/figS6_species_scaling.png",        # PL2
    "output/figures/supplement/figS15_occurrence_curves.png"
  ),
  "06_cooccurrence_analysis.R" = c(
    "output/objects/cooccurrence_results.rds",
    "output/objects/cafi_network.rds",
    "output/tables/pairwise_cooccurrence.csv",
    "output/tables/intraspecific_density.csv",
    "output/tables/size_dependent_cooccurrence.csv",
    "output/tables/network_metrics.csv",
    "output/tables/hub_species.csv",
    "output/figures/supplement/figS11_cooccurrence.png"
  ),
  "07_spatial_autocorrelation.R" = c(
    "output/tables/morans_i_results.csv",
    "output/figures/supplement/figS4_spatial_autocorrelation.png" # PL2
  ),
  "08_functional_groups.R" = c(
    "output/figures/supplement/figS9_functional_groups.png",
    "output/tables/taxonomic_group_scaling.csv",
    "output/objects/functional_analysis_results.rds"
  ),
  "09_cafi_condition_feedbacks.R" = c(
    "output/figures/manuscript/fig5_feedbacks.png",
    "output/figures/supplement/figS10_rarefaction_sensitivity.png",
    "output/figures/supplement/figS12_bef_variance_partitioning.png",
    "output/figures/supplement/figS13_condition_details.png",
    "output/figures/supplement/figS16_species_trait_heatmap.png",     # PL2
    "output/figures/supplement/figS17_species_trait_biplots.png",     # PL2
    "output/tables/cafi_condition_models.csv",
    "output/tables/key_species_effects.csv",                         # PL5
    "output/tables/species_trait_correlations.csv"                    # PL5
  ),
  "13_taxonomy_sensitivity.R" = c(
    "output/tables/taxonomy_sensitivity.csv",
    "output/tables/taxonomy_sensitivity_species_scaling.csv",
    "output/tables/network_topology_sensitivity.csv",
    "output/tables/network_edge_overlap.csv",
    "output/tables/network_hub_stability.csv",
    "output/figures/supplement/figS8_taxonomy_sensitivity.png",
    "output/figures/supplement/figS14_network_sensitivity.png"
  )
)

# ============================================================================
# LOGGING FUNCTIONS
# ============================================================================

#' Initialize pipeline log
init_log <- function(log_path) {
  # Create output directory if needed
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)

  # Write header
  header <- paste0(
    "============================================================\n",
    "CAFI SURVEY ANALYSIS PIPELINE LOG\n",
    "============================================================\n",
    "Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "R version: ", R.version.string, "\n",
    "Working directory: ", getwd(), "\n",
    "============================================================\n\n"
  )

  writeLines(header, log_path)
  return(log_path)
}

#' Write to log file
log_message <- function(log_path, msg, console = TRUE) {
  timestamp <- format(Sys.time(), "%H:%M:%S")
  full_msg <- paste0("[", timestamp, "] ", msg, "\n")
  cat(full_msg, file = log_path, append = TRUE)
  if (console) cat(full_msg)
}

#' Log section header
log_section <- function(log_path, title) {
  section <- paste0(
    "\n------------------------------------------------------------\n",
    title, "\n",
    "------------------------------------------------------------\n"
  )
  cat(section, file = log_path, append = TRUE)
  cat(section)
}

# ============================================================================
# DEPENDENCY CHECKING
# ============================================================================

#' Check if required objects exist in environment
check_objects <- function(required_objects) {
  missing <- character(0)
  for (obj in required_objects) {
    if (!exists(obj, envir = .GlobalEnv)) {
      missing <- c(missing, obj)
    }
  }
  return(missing)
}

#' Check if expected output files exist
check_outputs <- function(script_name) {
  expected <- EXPECTED_OUTPUTS[[script_name]]
  if (is.null(expected) || length(expected) == 0) return(TRUE)

  all_exist <- all(file.exists(expected))
  return(all_exist)
}

#' Verify dependencies for a script
verify_dependencies <- function(script_name) {
  deps <- list(
    "01_load_data.R" = list(objects = c("PATHS")),
    "02_community_analysis.R" = list(objects = c("coral_master", "cafi_clean", "community_matrix")),
    "03_landscape_characterization.R" = list(objects = c("coral_master")),
    "04_landscape_effects.R" = list(objects = c("coral_master")),
    "05_species_scaling_analysis.R" = list(objects = c("coral_master", "cafi_clean")),
    "06_cooccurrence_analysis.R" = list(objects = c("coral_master", "community_matrix", "cafi_clean")),
    "07_spatial_autocorrelation.R" = list(objects = c("coral_master", "community_matrix")),
    "08_functional_groups.R" = list(objects = c("coral_master", "cafi_clean")),
    "09_cafi_condition_feedbacks.R" = list(objects = c("coral_master", "cafi_clean", "community_matrix", "condition_scores")),
    "13_taxonomy_sensitivity.R" = list(objects = c("coral_master", "cafi_clean", "condition_scores", "otu_taxonomy", "taxonomy_scenario_data"))
  )

  script_deps <- deps[[script_name]]
  if (is.null(script_deps)) return(list(success = TRUE, missing = character(0)))

  missing <- check_objects(script_deps$objects)

  # PL3: Fallback — if object not in global env, check for corresponding RDS on disk
  if (length(missing) > 0) {
    still_missing <- character(0)
    for (obj in missing) {
      rds_path <- file.path("output", "objects", paste0(obj, ".rds"))
      if (file.exists(rds_path)) {
        # Load the RDS into global env so the script can use it
        assign(obj, readRDS(rds_path), envir = .GlobalEnv)
      } else {
        still_missing <- c(still_missing, obj)
      }
    }
    missing <- still_missing
  }

  return(list(
    success = length(missing) == 0,
    missing = missing
  ))
}

# ============================================================================
# SCRIPT EXECUTION
# ============================================================================

#' Execute a single script with error handling
run_script <- function(script_name, log_path, verbose = TRUE) {
  script_path <- file.path("scripts", script_name)

  # Check if script exists
  if (!file.exists(script_path)) {
    return(list(
      success = FALSE,
      error = paste("Script not found:", script_path),
      time = NA,
      skipped = TRUE
    ))
  }

  # Check dependencies
  dep_check <- verify_dependencies(script_name)
  if (!dep_check$success) {
    return(list(
      success = FALSE,
      error = paste("Missing dependencies:", paste(dep_check$missing, collapse = ", ")),
      time = NA,
      skipped = FALSE
    ))
  }

  # Execute script
  start_time <- Sys.time()

  # Use withCallingHandlers for warnings (logs but continues)
  # and tryCatch for errors (stops execution)
  warnings_collected <- character(0)

  result <- tryCatch({
    withCallingHandlers({
      source(script_path, local = FALSE, echo = FALSE)
      list(success = TRUE, error = NULL, warnings = warnings_collected)
    }, warning = function(w) {
      # Log warning but CONTINUE execution (don't return early)
      log_message(log_path, paste("WARNING:", conditionMessage(w)))
      warnings_collected <<- c(warnings_collected, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e), warnings = warnings_collected)
  })

  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "secs")

  result$time <- as.numeric(elapsed)
  result$skipped <- FALSE

  return(result)
}

#' Format elapsed time for display
format_time <- function(seconds) {
  if (is.na(seconds)) return("--")

  if (seconds < 60) {
    return(sprintf("%.1f sec", seconds))
  } else if (seconds < 3600) {
    mins <- floor(seconds / 60)
    secs <- round(seconds %% 60)
    return(sprintf("%d min %d sec", mins, secs))
  } else {
    hours <- floor(seconds / 3600)
    mins <- floor((seconds %% 3600) / 60)
    return(sprintf("%d hr %d min", hours, mins))
  }
}

# ============================================================================
# MAIN PIPELINE FUNCTION
# ============================================================================

#' Run the complete CAFI analysis pipeline
#'
#' @param skip_completed Skip scripts whose outputs already exist (default FALSE)
#' @param verbose Print detailed progress (default TRUE)
#' @param stop_on_error Stop on first error (default FALSE)
#' @param sections Which sections to run: "all", "core", "extended", "ml"
#'   Default runs core + extended (excludes ML exploration).
#'   Use sections = "all" to include ML scripts.
#' @return List with timing and status for each script
run_pipeline <- function(skip_completed = FALSE,
                         verbose = TRUE,
                         stop_on_error = FALSE,
                         sections = c("core", "extended")) {

  # Set working directory to project root
  if (requireNamespace("here", quietly = TRUE)) {
    setwd(here::here())
  }

  # Initialize logging
  log_path <- "output/pipeline.log"
  init_log(log_path)

  cat("\n")
  cat("============================================================\n")
  cat("    CAFI SURVEY ANALYSIS PIPELINE\n")
  cat("============================================================\n")
  cat("Starting:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Logging to:", log_path, "\n")
  cat("============================================================\n\n")

  # Determine which sections to run
  if ("all" %in% sections) {
    sections <- c("core", "extended", "ml")
  }

  # Collect all scripts to run
  scripts_to_run <- list()
  for (section in sections) {
    if (!is.null(PIPELINE_SCRIPTS[[section]])) {
      scripts_to_run <- c(scripts_to_run, PIPELINE_SCRIPTS[[section]])
    }
  }

  # Results tracking
  results <- list()
  total_time <- 0
  n_success <- 0
  n_failed <- 0
  n_skipped <- 0

  pipeline_start <- Sys.time()

  # Run each script
  for (i in seq_along(scripts_to_run)) {
    script_info <- scripts_to_run[[i]]
    script_name <- script_info$name
    script_desc <- script_info$desc

    log_section(log_path, sprintf("[%d/%d] %s", i, length(scripts_to_run), script_desc))
    log_message(log_path, paste("Script:", script_name))

    # Check if should skip
    if (skip_completed && check_outputs(script_name)) {
      log_message(log_path, "SKIPPED: Outputs already exist")
      n_skipped <- n_skipped + 1
      results[[script_name]] <- list(
        status = "skipped",
        time = NA,
        error = NULL
      )
      next
    }

    # Run script
    log_message(log_path, "Running...")
    result <- run_script(script_name, log_path, verbose)

    if (result$skipped) {
      log_message(log_path, paste("SKIPPED:", result$error))
      n_skipped <- n_skipped + 1
      results[[script_name]] <- list(
        status = "skipped",
        time = NA,
        error = result$error
      )
    } else if (result$success) {
      log_message(log_path, sprintf("SUCCESS (%.1f sec)", result$time))
      n_success <- n_success + 1
      total_time <- total_time + result$time
      results[[script_name]] <- list(
        status = "success",
        time = result$time,
        error = NULL
      )
    } else {
      log_message(log_path, paste("FAILED:", result$error))
      n_failed <- n_failed + 1
      results[[script_name]] <- list(
        status = "failed",
        time = result$time,
        error = result$error
      )

      # Stop on error if requested
      if (stop_on_error) {
        log_message(log_path, "Pipeline stopped due to error (stop_on_error = TRUE)")
        break
      }
    }

    cat("\n")
  }

  pipeline_end <- Sys.time()
  pipeline_duration <- difftime(pipeline_end, pipeline_start, units = "secs")

  # ============================================================================
  # SUMMARY REPORT
  # ============================================================================

  log_section(log_path, "PIPELINE SUMMARY")

  summary_msg <- sprintf(
    "Total scripts: %d\n  Successful: %d\n  Failed: %d\n  Skipped: %d\n  Total time: %s\n",
    length(scripts_to_run), n_success, n_failed, n_skipped, format_time(as.numeric(pipeline_duration))
  )
  log_message(log_path, summary_msg, console = TRUE)

  # Timing table
  timing_df <- data.frame(
    script = names(results),
    status = sapply(results, function(x) x$status),
    time_seconds = sapply(results, function(x) ifelse(is.null(x$time), NA, x$time)),
    error = sapply(results, function(x) ifelse(is.null(x$error), "", x$error)),
    stringsAsFactors = FALSE
  )

  timing_df$time_formatted <- sapply(timing_df$time_seconds, format_time)

  # Save timing table
  dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
  write.csv(timing_df, "output/tables/pipeline_timing.csv", row.names = FALSE)
  log_message(log_path, "Saved: output/tables/pipeline_timing.csv")

  # Print detailed results
  cat("\n")
  cat("Script Execution Times:\n")
  cat("-----------------------\n")
  for (i in 1:nrow(timing_df)) {
    status_icon <- switch(timing_df$status[i],
                          success = "[OK]",
                          failed = "[FAIL]",
                          skipped = "[SKIP]")
    cat(sprintf("  %s %-40s %s\n",
                status_icon,
                timing_df$script[i],
                timing_df$time_formatted[i]))
  }
  cat("\n")

  # List any failures
  if (n_failed > 0) {
    cat("FAILED SCRIPTS:\n")
    failed_scripts <- timing_df[timing_df$status == "failed", ]
    for (i in 1:nrow(failed_scripts)) {
      cat(sprintf("  - %s: %s\n", failed_scripts$script[i], failed_scripts$error[i]))
    }
    cat("\n")
  }

  # Final message
  log_message(log_path, sprintf("Pipeline completed at %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  cat("============================================================\n")
  if (n_failed == 0) {
    cat("    PIPELINE COMPLETED SUCCESSFULLY\n")
  } else {
    cat("    PIPELINE COMPLETED WITH ERRORS\n")
  }
  cat("============================================================\n\n")
  cat("Log file:", log_path, "\n")
  cat("Timing table: output/tables/pipeline_timing.csv\n\n")

  # Return results invisibly
  invisible(list(
    results = results,
    timing = timing_df,
    summary = list(
      total = length(scripts_to_run),
      success = n_success,
      failed = n_failed,
      skipped = n_skipped,
      duration = as.numeric(pipeline_duration)
    )
  ))
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

#' Run only the core pipeline
run_core_pipeline <- function(...) {
  run_pipeline(sections = "core", ...)
}

#' Run extended analyses only (assumes core is complete)
run_extended_analyses <- function(...) {
  run_pipeline(sections = "extended", ...)
}

#' Run ML exploration scripts only (not part of core hypothesis testing)
run_ml_exploration <- function(...) {
  run_pipeline(sections = "ml", ...)
}

#' Run everything including ML exploration
run_full_pipeline <- function(...) {
  run_pipeline(sections = "all", ...)
}

# ============================================================================
# FAST ITERATION FUNCTIONS
# ============================================================================

# All pipeline script names in order (core + extended)
ALL_PIPELINE_SCRIPTS <- c(
  sapply(PIPELINE_SCRIPTS$core, function(x) x$name),
  sapply(PIPELINE_SCRIPTS$extended, function(x) x$name)
)

#' Resolve a shorthand script reference to the full script name
#'
#' @param x Script reference: number ("09"), partial name ("09_cafi"), or full name
#' @return Full script filename (e.g. "09_cafi_condition_feedbacks.R")
resolve_script <- function(x) {
  # Already a full name?
  if (x %in% ALL_PIPELINE_SCRIPTS) return(x)

  # Match by number prefix (e.g. "09" -> "09_cafi_condition_feedbacks.R")
  matches <- ALL_PIPELINE_SCRIPTS[startsWith(ALL_PIPELINE_SCRIPTS, paste0(x, "_"))]
  if (length(matches) == 1) return(matches)

  # Match by partial name (e.g. "09_cafi" -> "09_cafi_condition_feedbacks.R")
  matches <- ALL_PIPELINE_SCRIPTS[startsWith(ALL_PIPELINE_SCRIPTS, x)]
  if (length(matches) == 1) return(matches)

  # Substring match
  matches <- ALL_PIPELINE_SCRIPTS[grepl(x, ALL_PIPELINE_SCRIPTS, fixed = TRUE)]
  if (length(matches) == 1) return(matches)

  if (length(matches) == 0) stop("No script matching '", x, "' found in pipeline.")
  stop("Ambiguous: '", x, "' matches multiple scripts: ", paste(matches, collapse = ", "))
}

#' Ensure setup and data are loaded (loads only what's missing)
ensure_setup_and_data <- function(script_name) {
  if (requireNamespace("here", quietly = TRUE)) setwd(here::here())

  # Always need setup

  if (!exists("PATHS", envir = .GlobalEnv)) {
    cat("Loading setup...\n")
    source("scripts/00_setup.R", local = FALSE)
  }

  # Check what this script needs
  dep_check <- verify_dependencies(script_name)
  if (!dep_check$success) {
    cat("Loading data (missing:", paste(dep_check$missing, collapse = ", "), ")...\n")
    source("scripts/01_load_data.R", local = FALSE)
  }
}

#' Run a single script with auto-loading of dependencies
#'
#' @param script Script reference: number ("09"), partial name ("09_cafi"), or full name
#' @examples
#'   run_one("09")    # CAFI-condition feedbacks (Fig 5)
#'   run_one("06")    # Co-occurrence analysis (Fig S11)
#'   run_one("13")    # Taxonomy sensitivity (Fig S8, S14)
run_one <- function(script) {
  script_name <- resolve_script(script)
  cat("\n=== Running:", script_name, "===\n\n")

  ensure_setup_and_data(script_name)

  start <- Sys.time()
  source(file.path("scripts", script_name), local = FALSE)
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  cat(sprintf("\n=== %s completed in %s ===\n", script_name, format_time(elapsed)))
  invisible(elapsed)
}

#' Run from a specific script through the end of the pipeline
#'
#' @param script Script reference: number ("06"), partial name, or full name
#' @examples
#'   run_from("06")   # Runs scripts 06, 07, 08, 09, 13 (~30 sec)
run_from <- function(script) {
  script_name <- resolve_script(script)
  idx <- match(script_name, ALL_PIPELINE_SCRIPTS)
  if (is.na(idx)) stop("Script '", script_name, "' not found in pipeline order.")

  scripts_to_run <- ALL_PIPELINE_SCRIPTS[idx:length(ALL_PIPELINE_SCRIPTS)]
  cat("\n=== Running", length(scripts_to_run), "scripts from", script_name, "===\n")
  cat("    Scripts:", paste(scripts_to_run, collapse = " -> "), "\n\n")

  ensure_setup_and_data(scripts_to_run[1])

  total_start <- Sys.time()
  for (s in scripts_to_run) {
    cat("--- Running:", s, "---\n")
    start <- Sys.time()
    source(file.path("scripts", s), local = FALSE)
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    cat(sprintf("    %s: %s\n\n", s, format_time(elapsed)))
  }

  total_elapsed <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))
  cat(sprintf("=== All %d scripts completed in %s ===\n",
              length(scripts_to_run), format_time(total_elapsed)))
  invisible(total_elapsed)
}

#' Run pipeline skipping slow scripts (02_community ~9 min, 05_scaling ~2.5 min)
#'
#' Requires that output objects from scripts 02 and 05 already exist from a
#' prior full pipeline run. Perfect for rapid figure iteration.
#'
#' @examples
#'   run_quick()  # ~1 min instead of ~12 min
run_quick <- function() {
  if (requireNamespace("here", quietly = TRUE)) setwd(here::here())

  slow_scripts <- c("02_community_analysis.R", "05_species_scaling_analysis.R")
  required_rds <- c(
    "output/objects/community_analysis_results.rds",
    "output/objects/scaling_analysis_results.rds"
  )

  missing_rds <- required_rds[!file.exists(required_rds)]
  if (length(missing_rds) > 0) {
    stop("run_quick() requires prior output from slow scripts.\n",
         "  Missing: ", paste(missing_rds, collapse = ", "), "\n",
         "  Run run_pipeline() once first to generate these files.")
  }

  fast_scripts <- setdiff(ALL_PIPELINE_SCRIPTS, slow_scripts)
  cat("\n=== QUICK PIPELINE (skipping 02 + 05) ===\n")
  cat("    Scripts:", paste(fast_scripts, collapse = " -> "), "\n\n")

  total_start <- Sys.time()
  for (s in fast_scripts) {
    cat("--- Running:", s, "---\n")
    start <- Sys.time()
    source(file.path("scripts", s), local = FALSE)
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    cat(sprintf("    %s: %s\n\n", s, format_time(elapsed)))
  }

  total_elapsed <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))
  cat(sprintf("=== Quick pipeline: %d scripts in %s (skipped 02, 05) ===\n",
              length(fast_scripts), format_time(total_elapsed)))
  invisible(total_elapsed)
}

#' Check pipeline status without running
check_pipeline_status <- function() {
  cat("\nCAFI Pipeline Status Check\n")
  cat("==========================\n\n")

  all_scripts <- c(
    unlist(lapply(PIPELINE_SCRIPTS$core, function(x) x$name)),
    unlist(lapply(PIPELINE_SCRIPTS$extended, function(x) x$name)),
    unlist(lapply(PIPELINE_SCRIPTS$ml, function(x) x$name))
  )

  for (script in all_scripts) {
    exists <- file.exists(file.path("scripts", script))
    outputs_exist <- check_outputs(script)

    status <- if (!exists) {
      "[NOT FOUND]"
    } else if (outputs_exist) {
      "[COMPLETE]"
    } else {
      "[PENDING]"
    }

    cat(sprintf("  %s %s\n", status, script))
  }
  cat("\n")
}

#' Compile all results into summary tables
compile_results <- function() {
  cat("\nCompiling pipeline results...\n")

  # This will be filled in after the pipeline runs
  # Check for scaling results
  if (file.exists("output/objects/scaling_analysis_results.rds")) {
    scaling_results <- readRDS("output/objects/scaling_analysis_results.rds")
    cat("  Loaded scaling analysis results\n")
  }

  # Check for network results
  if (file.exists("output/objects/cafi_network.rds")) {
    network_results <- readRDS("output/objects/cafi_network.rds")
    cat("  Loaded network analysis results\n")
  }

  # Check for community results
  if (file.exists("output/objects/community_analysis_results.rds")) {
    community_results <- readRDS("output/objects/community_analysis_results.rds")
    cat("  Loaded community analysis results\n")
  }

  cat("Results compilation complete.\n\n")
}

# ============================================================================
# AUTO-RUN IF CALLED AS SCRIPT
# ============================================================================

# Check if running as script (not sourced)
if (sys.nframe() == 0) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)

  skip_completed <- "--skip-completed" %in% args
  stop_on_error <- "--stop-on-error" %in% args
  quiet <- "--quiet" %in% args

  # Run pipeline
  run_pipeline(
    skip_completed = skip_completed,
    verbose = !quiet,
    stop_on_error = stop_on_error
  )
}
