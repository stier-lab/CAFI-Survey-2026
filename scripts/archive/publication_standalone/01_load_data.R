#!/usr/bin/env Rscript
# ============================================================================
# 01_load_data.R - Load and Prepare All Data for Publication Analyses
# ============================================================================
#
# PURPOSE: Load, clean, and prepare all datasets needed for the 6-figure
# publication plan. This script creates standardized data objects used by
# all downstream analysis scripts.
#
# OUTPUTS (saved to output/objects/):
#   - cafi_clean.rds: Clean individual CAFI records
#   - coral_master.rds: Coral characteristics with all metrics
#   - community_matrix.rds: Coral x OTU abundance matrix
#   - condition_scores.rds: Position-corrected coral condition (PC1)
#   - neighborhood_data.rds: Coral spatial neighborhood metrics
#   - functional_groups.rds: CAFI classified by functional group
#
# IMPORTANT NOTES:
#   - All corals are Pocillopora spp. (morphotypes not reliable species IDs)
#   - CAFI "species" are morphological OTUs (no genetic confirmation)
#   - Sampling position bias corrected for physiological metrics
#   - Sites: HAU (Hauru), MAT (Maatea), MRB (Moorea Barrier Reef)
#
# Author: CAFI Analysis Pipeline (Refactored)
# Date: 2025-12-03
# ============================================================================

cat("\n========================================\n")
cat("Loading and Preparing Publication Data\n")
cat("========================================\n\n")

# Load setup (paths, packages, theme)
source(here::here("scripts/publication/00_setup.R"))

# ============================================================================
# 1. LOAD RAW DATA FILES
# ============================================================================

cat("1. Loading raw data files...\n")

# CAFI individual records
cafi_raw <- read_csv(here("data/survey_cafi_data_w_taxonomy_summer2019_v5.csv"),
                     show_col_types = FALSE) %>%
  clean_names()

# Coral characteristics
coral_raw <- read_csv(here("data/survey_coral_characteristics_merged_v2.csv"),
                      show_col_types = FALSE) %>%
  clean_names()

# Physiology data
physio_raw <- read_csv(here("data/survey_master_phys_data_v3.csv"),
                       show_col_types = FALSE) %>%
  clean_names()

cat("  - CAFI records:", nrow(cafi_raw), "\n")
cat("  - Coral records:", nrow(coral_raw), "\n")
cat("  - Physiology records:", nrow(physio_raw), "\n\n")

# ============================================================================
# 2. CLEAN CAFI DATA
# ============================================================================

cat("2. Cleaning CAFI data...\n")

cafi_clean <- cafi_raw %>%
  filter(!is.na(coral_id)) %>%
  mutate(
    # Standardize type
    type = tolower(type),

    # Create OTU identifier (morphological groupings)
    # NOTE: These are NOT genetic species - interpret cautiously
    otu = case_when(
      !is.na(search_term) & search_term != "" ~ str_trim(search_term),
      !is.na(code) ~ paste(type, code, sep = "_"),
      TRUE ~ type
    ),

    # Size in mm
    size_mm = as.numeric(cafi_size_mm),

    # Extract site from coral_id prefix
    site = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ NA_character_
    ),

    # Assign functional group (biologically meaningful)
    functional_group = map2_chr(otu, type, assign_functional_group)
  ) %>%
  filter(!is.na(site))

cat("  - Clean CAFI records:", nrow(cafi_clean), "\n")
cat("  - Unique OTUs:", n_distinct(cafi_clean$otu), "\n")
cat("  - Unique corals:", n_distinct(cafi_clean$coral_id), "\n")

# Summarize by OTU and coral
cafi_by_coral <- cafi_clean %>%
  group_by(coral_id, site) %>%
  summarise(
    total_cafi = n(),
    otu_richness = n_distinct(otu),
    shannon = vegan::diversity(table(otu)),

    # Functional group counts
    n_trapezia = sum(functional_group == "Trapezia"),
    n_resident_fish = sum(functional_group == "Resident Fish"),
    n_corallivore = sum(functional_group == "Corallivore"),
    n_other_crab = sum(functional_group == "Other Crab"),
    n_shrimp = sum(functional_group == "Shrimp"),
    n_other = sum(functional_group == "Other"),

    # Taxonomic group counts
    n_crabs = sum(type == "crab"),
    n_shrimps = sum(type == "shrimp"),
    n_fish = sum(type == "fish"),
    n_snails = sum(type == "snail"),

    # Mean size
    mean_cafi_size = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  )

cat("  - CAFI summary by coral created\n\n")

# ============================================================================
# 3. CLEAN CORAL DATA
# ============================================================================

cat("3. Cleaning coral characteristics...\n")

coral_clean <- coral_raw %>%
  distinct(coral_id, .keep_all = TRUE) %>%
  mutate(
    # Extract site
    site = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ NA_character_
    ),

    # Branch width is the real measurable trait
    branch_width = if ("branch_width" %in% names(.)) tolower(branch_width) else NA_character_,

    # Convert depth
    depth_m = as.numeric(depth),

    # Calculate coral volume (best available measurement)
    volume = coalesce(volume_lab, volume_field,
                      length_lab * width_lab * height_lab),
    log_volume = log10(volume + 1),

    # Coral dimensions
    height = coalesce(height_lab, height_field),
    width = coalesce(width_lab, width_field),
    length = coalesce(length_lab, length_field),

    # Extract sampling position for bias correction
    stump_length = coalesce(nub1_stump_length, nub2_stump_length),
    nubbin_length = coalesce(nub1_nubbin_length, nub2_nubbin_length),

    # Neighborhood metrics (from field measurements)
    n_neighbors = number_of_neighbors,
    mean_neighbor_dist = mean_neighbor_distance,
    total_neighbor_volume = combined_total_volume_of_neighbors,
    mean_neighbor_volume = mean_total_volume_of_neighbors
  ) %>%
  filter(!is.na(site), !is.na(volume), volume > 0)

cat("  - Clean coral records:", nrow(coral_clean), "\n")
cat("  - Sites:", paste(unique(coral_clean$site), collapse = ", "), "\n")

# ============================================================================
# 4. CREATE MASTER CORAL DATASET
# ============================================================================

cat("\n4. Creating master coral dataset...\n")

# Merge coral characteristics with CAFI metrics
coral_master <- coral_clean %>%
  left_join(cafi_by_coral, by = c("coral_id", "site")) %>%
  mutate(
    # Replace NA CAFI counts with 0
    across(starts_with("n_") | starts_with("total_") | starts_with("otu_"),
           ~replace_na(., 0)),

    # Create size categories
    size_class = cut(volume,
                     breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE),

    # Neighborhood derived metrics
    isolation_index = mean_neighbor_dist / (volume^(1/3) + 1),
    crowding_index = total_neighbor_volume / (mean_neighbor_dist + 1),
    relative_size = volume / (mean_neighbor_volume + 1),

    # CAFI density (per unit coral volume)
    cafi_density = total_cafi / volume * 1000  # per 1000 cm3
  )

cat("  - Master dataset rows:", nrow(coral_master), "\n")
cat("  - Master dataset columns:", ncol(coral_master), "\n")

# ============================================================================
# 5. POSITION-CORRECTED CORAL CONDITION
# ============================================================================

cat("\n5. Computing position-corrected coral condition scores...\n")

# Physiology variables (surface-area normalized)
physio_vars <- c("protein_mg_cm2", "carb_mg_cm2", "zoox_cells_cm2", "afdw_mg_cm2")

# Merge physiology with coral data
physio_merged <- physio_raw %>%
  left_join(coral_clean %>% select(coral_id, site, volume, stump_length),
            by = "coral_id") %>%
  filter(!is.na(stump_length))

# Check we have enough data
if (nrow(physio_merged) > 20) {

  # Position correction: regress each trait on stump_length, use residuals
  # This removes the sampling position bias

  correction_data <- physio_merged %>%
    select(coral_id, site, volume, stump_length, any_of(physio_vars)) %>%
    drop_na()

  # Extract position-corrected residuals for each trait
  corrected_traits <- correction_data %>%
    select(coral_id, site, volume, stump_length)

  for (var in physio_vars) {
    if (var %in% names(correction_data) && sum(!is.na(correction_data[[var]])) > 20) {
      # Model: trait ~ stump_length
      model <- lm(reformulate("stump_length", response = var), data = correction_data)
      # Extract and standardize residuals
      resid_z <- scale(residuals(model))[, 1]

      # Store with shortened name
      var_name <- case_when(
        var == "protein_mg_cm2" ~ "protein_corr",
        var == "carb_mg_cm2" ~ "carb_corr",
        var == "zoox_cells_cm2" ~ "zoox_corr",
        var == "afdw_mg_cm2" ~ "afdw_corr"
      )
      corrected_traits[[var_name]] <- resid_z
    }
  }

  # PCA on corrected traits to create condition score
  corrected_vars <- c("protein_corr", "carb_corr", "zoox_corr", "afdw_corr")
  corrected_vars <- corrected_vars[corrected_vars %in% names(corrected_traits)]

  if (length(corrected_vars) >= 3) {
    pca_data <- corrected_traits %>%
      select(all_of(corrected_vars)) %>%
      drop_na()

    if (nrow(pca_data) > 10) {
      pca_result <- prcomp(pca_data, scale. = FALSE, center = TRUE)

      # PC1 loadings - flip if needed so higher = better
      loadings <- pca_result$rotation[, 1]
      if (sum(loadings < 0) > sum(loadings > 0)) {
        condition_pc1 <- -pca_result$x[, 1]
      } else {
        condition_pc1 <- pca_result$x[, 1]
      }

      # Variance explained
      var_explained <- summary(pca_result)$importance[2, 1] * 100

      # Add to corrected_traits
      corrected_traits$condition_score <- NA
      corrected_traits$condition_score[complete.cases(pca_data)] <- condition_pc1

      cat("  - Position correction applied to", length(corrected_vars), "traits\n")
      cat("  - Condition score (PC1) explains", round(var_explained, 1), "% variance\n")
    }
  }

  # Save condition scores
  condition_scores <- corrected_traits %>%
    select(coral_id, site, condition_score, any_of(corrected_vars))

  save_object(condition_scores, "condition_scores.rds")

  # Add to master dataset
  coral_master <- coral_master %>%
    left_join(condition_scores %>% select(coral_id, condition_score),
              by = "coral_id")

} else {
  cat("  - Insufficient physiology data for position correction\n")
  coral_master$condition_score <- NA
}

# ============================================================================
# 6. CREATE COMMUNITY MATRIX
# ============================================================================

cat("\n6. Creating community matrix (coral x OTU)...\n")

community_matrix <- cafi_clean %>%
  group_by(coral_id, otu) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    id_cols = coral_id,
    names_from = otu,
    values_from = count,
    values_fill = 0
  ) %>%
  column_to_rownames("coral_id") %>%
  as.matrix()

cat("  - Matrix dimensions:", nrow(community_matrix), "x", ncol(community_matrix), "\n")

# ============================================================================
# 7. CREATE FUNCTIONAL GROUP SUMMARY
# ============================================================================

cat("\n7. Creating functional group summary...\n")

functional_summary <- cafi_clean %>%
  group_by(functional_group) %>%
  summarise(
    n_individuals = n(),
    n_otus = n_distinct(otu),
    n_corals = n_distinct(coral_id),
    mean_size_mm = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    proportion = n_individuals / sum(n_individuals)
  ) %>%
  arrange(desc(n_individuals))

cat("  - Functional groups:\n")
for (i in 1:nrow(functional_summary)) {
  cat(sprintf("    %s: %d individuals (%.1f%%)\n",
              functional_summary$functional_group[i],
              functional_summary$n_individuals[i],
              functional_summary$proportion[i] * 100))
}

# ============================================================================
# 8. SAVE ALL PROCESSED DATA
# ============================================================================

cat("\n8. Saving processed data objects...\n")

save_object(cafi_clean, "cafi_clean.rds")
save_object(coral_master, "coral_master.rds")
save_object(community_matrix, "community_matrix.rds")
save_object(functional_summary, "functional_summary.rds")
save_object(cafi_by_coral, "cafi_by_coral.rds")

# Also save key tables
save_table(functional_summary, "functional_group_summary.csv")

otu_summary <- cafi_clean %>%
  group_by(otu, type, functional_group) %>%
  summarise(
    n_individuals = n(),
    n_corals = n_distinct(coral_id),
    mean_size_mm = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_individuals))

save_table(otu_summary, "otu_abundance_ranking.csv")

# ============================================================================
# DATA SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("Data Loading Summary\n")
cat("========================================\n\n")

cat("Dataset Overview:\n")
cat("  - Study: CAFI Survey Mo'orea Summer 2019\n")
cat("  - Coral species: Pocillopora spp. (all morphotypes pooled)\n")
cat("  - Sites: HAU, MAT, MRB\n\n")

cat("Sample Sizes:\n")
cat("  - Total corals:", nrow(coral_master), "\n")
cat("  - Total CAFI individuals:", nrow(cafi_clean), "\n")
cat("  - Unique OTUs:", n_distinct(cafi_clean$otu), "\n")
cat("  - Corals with condition scores:",
    sum(!is.na(coral_master$condition_score)), "\n\n")

cat("CAFI by Site:\n")
site_cafi <- coral_master %>%
  group_by(site) %>%
  summarise(
    n_corals = n(),
    total_cafi = sum(total_cafi),
    mean_cafi = round(mean(total_cafi), 1),
    mean_richness = round(mean(otu_richness), 1),
    .groups = "drop"
  )
print(site_cafi)

cat("\nData files saved to:", OBJECTS_DIR, "\n")
cat("Tables saved to:", TABLES_DIR, "\n\n")
cat("Ready for analysis scripts!\n\n")
