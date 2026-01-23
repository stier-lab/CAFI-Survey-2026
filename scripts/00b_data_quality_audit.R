# ============================================================================
# 00b_data_quality_audit.R - Data Quality Audit for CAFI Survey
# ============================================================================
#
# PURPOSE: Comprehensive data quality assessment including:
#   1. Physiology data missingness analysis
#   2. Trait database integration
#   3. GPS/Spatial data validation
#   4. Data quality report generation
#
# USAGE:
#   source("scripts/00_setup.R")  # Must run first
#   source("scripts/00b_data_quality_audit.R")
#
# INPUT FILES:
#   - data/survey_master_phys_data_v3.csv
#   - data/traits/cafi_traits_final.csv
#   - data/survey_coral_characteristics_merged_v2.csv
#   - data/survey_cafi_data_w_taxonomy_summer2019_v5.csv
#
# OUTPUTS:
#   - output/objects/cafi_with_traits.rds
#   - output/tables/physiology_missingness_summary.csv
#   - output/tables/spatial_validation_summary.csv
#   - output/tables/data_quality_report.csv
#   - output/tables/outlier_flags.csv
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-01-21
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI SURVEY - Data Quality Audit\n")
cat("============================================================\n\n")

# Check that setup has been run
if (!exists("PATHS")) {
  source(here::here("scripts/00_setup.R"))
}

# ============================================================================
# 1. LOAD RAW DATA FILES
# ============================================================================

cat("1. Loading raw data files...\n")

physio_raw <- read_csv(PATHS$physio_data, show_col_types = FALSE) %>%
  clean_names()

coral_raw <- read_csv(PATHS$coral_data, show_col_types = FALSE) %>%
  clean_names()

cafi_raw <- read_csv(PATHS$cafi_data, show_col_types = FALSE) %>%
  clean_names()

traits_raw <- read_csv(here("data", "traits", "cafi_traits_final.csv"),
                       show_col_types = FALSE) %>%
  clean_names()

cat("   Physiology records:", nrow(physio_raw), "\n")
cat("   Coral records:", nrow(coral_raw), "\n")
cat("   CAFI records:", nrow(cafi_raw), "\n")
cat("   Trait records:", nrow(traits_raw), "\n\n")

# ============================================================================
# 2. PHYSIOLOGY DATA MISSINGNESS ASSESSMENT
# ============================================================================

cat("2. Assessing physiology data missingness...\n")

# Define key physiological variables
physio_vars <- c("protein_mg_cm2", "carb_mg_cm2", "zoox_cells_cm2", "afdw_mg_cm2")

# Calculate missingness for each variable
missingness_summary <- tibble(
  variable = physio_vars,
  n_total = nrow(physio_raw),
  n_missing = map_int(physio_vars, ~sum(is.na(physio_raw[[.x]]))),
  n_present = n_total - n_missing,
  pct_missing = round(n_missing / n_total * 100, 1),
  pct_present = round(n_present / n_total * 100, 1)
)

# Add additional physiology variables that may have missingness
additional_vars <- c("surface_area", "dry_weight", "haplotype")
additional_missingness <- tibble(
  variable = additional_vars,
  n_total = nrow(physio_raw),
  n_missing = map_int(additional_vars, ~sum(is.na(physio_raw[[.x]]))),
  n_present = n_total - n_missing,
  pct_missing = round(n_missing / n_total * 100, 1),
  pct_present = round(n_present / n_total * 100, 1)
)

missingness_summary <- bind_rows(missingness_summary, additional_missingness)

cat("   Missingness Summary:\n")
print(missingness_summary)
cat("\n")

# Identify corals missing physiology data
corals_with_physio <- physio_raw %>%
  filter(!is.na(coral_id)) %>%
  pull(coral_id) %>%
  unique()

all_corals <- coral_raw %>%
  filter(!is.na(coral_id)) %>%
  pull(coral_id) %>%
  unique()

corals_missing_physio <- setdiff(all_corals, corals_with_physio)

cat("   Corals in coral dataset:", length(all_corals), "\n")
cat("   Corals with physiology data:", length(corals_with_physio), "\n")
cat("   Corals missing physiology:", length(corals_missing_physio), "\n")

if (length(corals_missing_physio) > 0) {
  cat("   Missing coral IDs:", paste(head(corals_missing_physio, 10), collapse = ", "))
  if (length(corals_missing_physio) > 10) cat("...(", length(corals_missing_physio) - 10, " more)")
  cat("\n")
}

# Document which corals are missing each physiology variable
corals_missing_by_var <- map_dfr(physio_vars, function(var) {
  missing_ids <- physio_raw %>%
    filter(is.na(.data[[var]])) %>%
    pull(coral_id) %>%
    unique()

  tibble(
    variable = var,
    n_corals_missing = length(missing_ids),
    coral_ids = paste(missing_ids, collapse = "; ")
  )
})

# Missingness pattern analysis
# Note: Full Little's MCAR test requires the 'naniar' or 'misty' package
# Here we do a simplified assessment based on correlation with other variables

cat("\n   Missingness Pattern Analysis:\n")

# Check if missingness correlates with coral volume (site from coral_raw)
physio_with_coral <- physio_raw %>%
  left_join(coral_raw %>% dplyr::select(coral_id, volume_lab, volume_field, site, survey_type),
            by = "coral_id")

# Create missingness indicators
for (var in physio_vars) {
  physio_with_coral[[paste0(var, "_missing")]] <- is.na(physio_with_coral[[var]])
}

# Test if missingness is associated with site (potential MAR indicator)
pattern_assessment <- tibble(
  variable = physio_vars,
  pattern_type = character(length(physio_vars)),
  notes = character(length(physio_vars))
)

for (i in seq_along(physio_vars)) {
  var <- physio_vars[i]
  missing_col <- paste0(var, "_missing")

  # Check association with site
  if (sum(!is.na(physio_with_coral$site)) > 0) {
    site_table <- table(physio_with_coral$site, physio_with_coral[[missing_col]])
    if (all(dim(site_table) > 1) && min(site_table) > 0) {
      chisq_test <- tryCatch(
        chisq.test(site_table),
        error = function(e) NULL,
        warning = function(w) suppressWarnings(chisq.test(site_table))
      )

      if (!is.null(chisq_test) && chisq_test$p.value < 0.05) {
        pattern_assessment$pattern_type[i] <- "Potentially MAR"
        pattern_assessment$notes[i] <- paste("Associated with site (p =",
                                             round(chisq_test$p.value, 3), ")")
      } else {
        pattern_assessment$pattern_type[i] <- "Potentially MCAR"
        pattern_assessment$notes[i] <- "No significant association with site"
      }
    } else {
      pattern_assessment$pattern_type[i] <- "Insufficient data"
      pattern_assessment$notes[i] <- "Cannot test association"
    }
  } else {
    pattern_assessment$pattern_type[i] <- "Unknown"
    pattern_assessment$notes[i] <- "Site data not available"
  }
}

cat("   Pattern assessment:\n")
print(pattern_assessment)
cat("\n")

# Save missingness results
missingness_full <- missingness_summary %>%
  left_join(pattern_assessment, by = "variable")

save_table(missingness_full, "physiology_missingness_summary")
save_table(corals_missing_by_var, "corals_missing_physiology_by_variable")

# ============================================================================
# 3. TRAIT DATABASE INTEGRATION
# ============================================================================

cat("3. Integrating trait database...\n")

# Examine trait database structure
cat("   Trait database columns:\n")
cat("   ", paste(names(traits_raw)[1:10], collapse = ", "), "...\n\n")

# Key traits to extract
key_trait_cols <- c("taxon_input", "accepted_name", "aphia_id", "rank",
                    "pocillopora_association", "functional_role_str",
                    "max_length_final_mm", "depth_min_final", "depth_max_final",
                    "trophic_level", "trophic_se")

# Check which key columns exist
available_traits <- intersect(key_trait_cols, names(traits_raw))
missing_traits <- setdiff(key_trait_cols, names(traits_raw))

cat("   Available key traits:", length(available_traits), "\n")
if (length(missing_traits) > 0) {
  cat("   Missing traits:", paste(missing_traits, collapse = ", "), "\n")
}

# Clean trait data for merging
traits_clean <- traits_raw %>%
  dplyr::select(any_of(key_trait_cols)) %>%
  rename(
    species_name = taxon_input,
    max_length_mm = max_length_final_mm,
    depth_min = depth_min_final,
    depth_max = depth_max_final,
    functional_role = functional_role_str,
    association_type = pocillopora_association
  ) %>%
  mutate(
    depth_range = ifelse(!is.na(depth_min) & !is.na(depth_max),
                         paste0(depth_min, "-", depth_max, " m"),
                         NA_character_)
  )

cat("   Clean trait records:", nrow(traits_clean), "\n")

# Clean CAFI data for merging - use search_term as the species identifier
cafi_for_merge <- cafi_raw %>%
  filter(!is.na(coral_id), !is.na(search_term)) %>%
  mutate(
    species_name = str_trim(search_term),
    site = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(coral_id, site, type, species_name, cafi_size_mm,
         genus, species, family, lowest_level)

# Merge traits into CAFI records
cafi_with_traits <- cafi_for_merge %>%
  left_join(traits_clean, by = "species_name")

# Calculate merge success rate
n_matched <- sum(!is.na(cafi_with_traits$accepted_name))
n_total <- nrow(cafi_with_traits)
match_rate <- round(n_matched / n_total * 100, 1)

cat("   CAFI records with traits matched:", n_matched, "/", n_total,
    "(", match_rate, "%)\n")

# Summary of trait coverage by taxonomic type
trait_coverage <- cafi_with_traits %>%
  group_by(type) %>%
  summarise(
    n_records = n(),
    n_with_traits = sum(!is.na(accepted_name)),
    pct_coverage = round(n_with_traits / n_records * 100, 1),
    n_with_max_length = sum(!is.na(max_length_mm)),
    n_with_depth = sum(!is.na(depth_range)),
    n_with_trophic = sum(!is.na(trophic_level)),
    .groups = "drop"
  )

cat("\n   Trait coverage by taxonomic type:\n")
print(trait_coverage)

# Identify unmatched species
unmatched_species <- cafi_with_traits %>%
  filter(is.na(accepted_name)) %>%
  distinct(species_name) %>%
  pull(species_name)

cat("\n   Unmatched species (", length(unmatched_species), "):\n")
if (length(unmatched_species) > 0) {
  cat("   ", paste(head(unmatched_species, 15), collapse = ", "))
  if (length(unmatched_species) > 15) cat("... (", length(unmatched_species) - 15, " more)")
  cat("\n")
}

# Save enhanced dataset
save_object(cafi_with_traits, "cafi_with_traits")
save_table(trait_coverage, "trait_coverage_by_type")

# ============================================================================
# 4. GPS/SPATIAL DATA VALIDATION
# ============================================================================

cat("\n4. Validating GPS/Spatial data...\n")

# Check GPS coordinates in coral characteristics
coral_spatial <- coral_raw %>%
  mutate(
    has_lat = !is.na(lat) & lat != "NA",
    has_long = !is.na(long) & long != "NA",
    has_gps = has_lat & has_long,
    lat_num = suppressWarnings(as.numeric(lat)),
    long_num = suppressWarnings(as.numeric(long))
  )

# GPS coverage summary
gps_summary <- coral_spatial %>%
  summarise(
    n_total = n(),
    n_with_gps = sum(has_gps, na.rm = TRUE),
    n_missing_gps = sum(!has_gps, na.rm = TRUE),
    pct_with_gps = round(n_with_gps / n_total * 100, 1)
  )

cat("   GPS coverage:\n")
cat("     Total corals:", gps_summary$n_total, "\n")
cat("     With GPS:", gps_summary$n_with_gps, "(", gps_summary$pct_with_gps, "%)\n")
cat("     Missing GPS:", gps_summary$n_missing_gps, "\n")

# Document neighborhood vs size-only corals
survey_type_summary <- coral_raw %>%
  count(survey_type, name = "n_corals") %>%
  mutate(pct = round(n_corals / sum(n_corals) * 100, 1))

cat("\n   Survey type breakdown:\n")
print(survey_type_summary)

# Identify the 63 neighborhood corals and 51 size-only corals
neighborhood_corals <- coral_raw %>%
  filter(survey_type == "neighborhood") %>%
  pull(coral_id)

size_only_corals <- coral_raw %>%
  filter(survey_type == "size") %>%
  pull(coral_id)

cat("\n   Neighborhood survey corals (5m surveys):", length(neighborhood_corals), "\n")
cat("   Size-only corals:", length(size_only_corals), "\n")

# Calculate spatial extent of study
spatial_extent <- coral_spatial %>%
  filter(has_gps) %>%
  summarise(
    lat_min = min(lat_num, na.rm = TRUE),
    lat_max = max(lat_num, na.rm = TRUE),
    long_min = min(long_num, na.rm = TRUE),
    long_max = max(long_num, na.rm = TRUE),
    lat_range = lat_max - lat_min,
    long_range = long_max - long_min
  )

if (any(!is.na(spatial_extent$lat_min))) {
  cat("\n   Spatial extent of study:\n")
  cat("     Latitude range:", round(spatial_extent$lat_min, 5), "to",
      round(spatial_extent$lat_max, 5), "\n")
  cat("     Longitude range:", round(spatial_extent$long_min, 5), "to",
      round(spatial_extent$long_max, 5), "\n")

  # Approximate distance in km (rough calculation for Moorea's latitude)
  lat_km <- abs(spatial_extent$lat_range) * 111  # ~111 km per degree latitude
  long_km <- abs(spatial_extent$long_range) * 111 * cos(abs(spatial_extent$lat_min) * pi / 180)
  cat("     Approximate extent: ~", round(lat_km, 2), "km (N-S) x ~",
      round(long_km, 2), "km (E-W)\n")
}

# Spatial summary by site
spatial_by_site <- coral_spatial %>%
  filter(has_gps) %>%
  mutate(
    site_code = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ "Unknown"
    )
  ) %>%
  group_by(site_code) %>%
  summarise(
    n_corals = n(),
    lat_mean = mean(lat_num, na.rm = TRUE),
    long_mean = mean(long_num, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n   Corals with GPS by site:\n")
print(spatial_by_site)

# Save spatial validation results
spatial_validation <- bind_rows(
  tibble(metric = "Total corals", value = gps_summary$n_total),
  tibble(metric = "Corals with GPS", value = gps_summary$n_with_gps),
  tibble(metric = "Corals missing GPS", value = gps_summary$n_missing_gps),
  tibble(metric = "Neighborhood survey corals", value = length(neighborhood_corals)),
  tibble(metric = "Size-only corals", value = length(size_only_corals))
)

save_table(spatial_validation, "spatial_validation_summary")

# ============================================================================
# 5. DATA QUALITY REPORT
# ============================================================================

cat("\n5. Generating data quality report...\n")

# Summary statistics for key variables in coral data
coral_clean <- coral_raw %>%
  mutate(
    volume = coalesce(volume_lab, volume_field),
    depth_m = suppressWarnings(as.numeric(depth)),
    site_code = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(!is.na(coral_id))

# Key variable summaries
summary_stats <- coral_clean %>%
  summarise(
    n_corals = n(),

    # Volume statistics
    volume_n = sum(!is.na(volume)),
    volume_mean = round(mean(volume, na.rm = TRUE), 2),
    volume_sd = round(sd(volume, na.rm = TRUE), 2),
    volume_min = round(min(volume, na.rm = TRUE), 2),
    volume_max = round(max(volume, na.rm = TRUE), 2),
    volume_median = round(median(volume, na.rm = TRUE), 2),

    # Depth statistics
    depth_n = sum(!is.na(depth_m)),
    depth_mean = round(mean(depth_m, na.rm = TRUE), 2),
    depth_sd = round(sd(depth_m, na.rm = TRUE), 2),
    depth_min = min(depth_m, na.rm = TRUE),
    depth_max = max(depth_m, na.rm = TRUE)
  )

# Reshape for report
quality_report <- tibble(
  metric = c(
    "Total corals",
    "Corals with volume", "Volume mean (cm3)", "Volume SD",
    "Volume min", "Volume max", "Volume median",
    "Corals with depth", "Depth mean (m)", "Depth SD",
    "Depth min (m)", "Depth max (m)"
  ),
  value = c(
    summary_stats$n_corals,
    summary_stats$volume_n, summary_stats$volume_mean, summary_stats$volume_sd,
    summary_stats$volume_min, summary_stats$volume_max, summary_stats$volume_median,
    summary_stats$depth_n, summary_stats$depth_mean, summary_stats$depth_sd,
    summary_stats$depth_min, summary_stats$depth_max
  )
)

cat("   Summary statistics:\n")
print(quality_report)

# ============================================================================
# 6. OUTLIER DETECTION
# ============================================================================

cat("\n6. Detecting potential outliers...\n")

# Function to detect outliers using IQR method
detect_outliers_iqr <- function(x, multiplier = 1.5) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - multiplier * iqr
  upper <- q3 + multiplier * iqr
  x < lower | x > upper
}

# Check for outliers in key variables
outlier_flags <- coral_clean %>%
  mutate(
    volume_outlier = detect_outliers_iqr(volume),
    depth_outlier = detect_outliers_iqr(depth_m)
  ) %>%
  dplyr::select(coral_id, site_code, volume, volume_outlier, depth_m, depth_outlier) %>%
  filter(volume_outlier | depth_outlier)

cat("   Volume outliers:", sum(coral_clean$volume > 0 &
                               detect_outliers_iqr(coral_clean$volume), na.rm = TRUE), "\n")
cat("   Depth outliers:", sum(detect_outliers_iqr(coral_clean$depth_m), na.rm = TRUE), "\n")

if (nrow(outlier_flags) > 0) {
  cat("\n   Flagged records:\n")
  print(head(outlier_flags, 10))
}

# Check physiology data for outliers
physio_outliers <- physio_raw %>%
  mutate(
    protein_outlier = detect_outliers_iqr(protein_mg_cm2),
    carb_outlier = detect_outliers_iqr(carb_mg_cm2),
    zoox_outlier = detect_outliers_iqr(zoox_cells_cm2),
    afdw_outlier = detect_outliers_iqr(afdw_mg_cm2)
  ) %>%
  dplyr::select(coral_id, protein_mg_cm2, protein_outlier,
         carb_mg_cm2, carb_outlier,
         zoox_cells_cm2, zoox_outlier,
         afdw_mg_cm2, afdw_outlier) %>%
  filter(protein_outlier | carb_outlier | zoox_outlier | afdw_outlier)

cat("\n   Physiology outliers:\n")
cat("     Protein:", sum(detect_outliers_iqr(physio_raw$protein_mg_cm2), na.rm = TRUE), "\n")
cat("     Carbohydrate:", sum(detect_outliers_iqr(physio_raw$carb_mg_cm2), na.rm = TRUE), "\n")
cat("     Zooxanthellae:", sum(detect_outliers_iqr(physio_raw$zoox_cells_cm2), na.rm = TRUE), "\n")
cat("     AFDW:", sum(detect_outliers_iqr(physio_raw$afdw_mg_cm2), na.rm = TRUE), "\n")

# Save outlier reports
save_table(outlier_flags, "outlier_flags_coral")
save_table(physio_outliers, "outlier_flags_physiology")

# ============================================================================
# 7. DATA COMPLETENESS SUMMARY
# ============================================================================

cat("\n7. Data completeness summary...\n")

# Overall completeness
completeness_coral <- coral_clean %>%
  summarise(
    across(everything(), ~sum(!is.na(.)) / n() * 100)
  ) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "pct_complete") %>%
  arrange(pct_complete) %>%
  mutate(pct_complete = round(pct_complete, 1))

# Show variables with lowest completeness
low_completeness <- completeness_coral %>%
  filter(pct_complete < 90) %>%
  head(15)

if (nrow(low_completeness) > 0) {
  cat("   Variables with <90% completeness:\n")
  print(low_completeness)
}

# Save final quality report
quality_report_final <- bind_rows(
  quality_report,
  tibble(
    metric = c(
      "Physiology variables assessed",
      "Trait records available",
      "CAFI-trait match rate (%)",
      "Coral outliers flagged",
      "Physiology outliers flagged"
    ),
    value = c(
      length(physio_vars),
      nrow(traits_clean),
      match_rate,
      nrow(outlier_flags),
      nrow(physio_outliers)
    )
  )
)

save_table(quality_report_final, "data_quality_report")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    Data Quality Audit Complete\n")
cat("============================================================\n\n")

cat("Key Findings:\n")
cat("  1. Physiology data: ", missingness_summary$pct_present[1], "% complete for protein\n", sep = "")
cat("  2. Trait matching: ", match_rate, "% of CAFI records matched to traits\n", sep = "")
cat("  3. GPS coverage: ", gps_summary$pct_with_gps, "% of corals have coordinates\n", sep = "")
cat("  4. Survey types: ", length(neighborhood_corals), " neighborhood + ",
    length(size_only_corals), " size-only corals\n", sep = "")
cat("  5. Outliers: ", nrow(outlier_flags), " coral + ",
    nrow(physio_outliers), " physiology records flagged\n\n", sep = "")

cat("Outputs saved to:\n")
cat("  - ", PATHS$objects, "/cafi_with_traits.rds\n", sep = "")
cat("  - ", PATHS$tables, "/physiology_missingness_summary.csv\n", sep = "")
cat("  - ", PATHS$tables, "/spatial_validation_summary.csv\n", sep = "")
cat("  - ", PATHS$tables, "/data_quality_report.csv\n", sep = "")
cat("  - ", PATHS$tables, "/outlier_flags_coral.csv\n", sep = "")
cat("  - ", PATHS$tables, "/outlier_flags_physiology.csv\n\n", sep = "")

cat("Next: source('scripts/01_load_data.R')\n\n")
