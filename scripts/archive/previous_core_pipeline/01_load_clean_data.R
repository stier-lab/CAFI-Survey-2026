#!/usr/bin/env Rscript
# ============================================================================
# 01_load_clean_data.R - Load and Clean Data
# ============================================================================
#
# PURPOSE: Load, clean, and prepare all datasets for the CAFI Survey analysis.
#
# USAGE: This script sources 00_setup.R and should be run before any analysis.
#        source(here::here("scripts/01_load_clean_data.R"))
#
# INPUT FILES (from data/):
#   - survey_cafi_data_w_taxonomy_summer2019_v5.csv
#   - survey_coral_characteristics_merged_v2.csv
#   - survey_master_phys_data_v3.csv
#
# OUTPUTS (saved to output/objects/):
#   - cafi_clean.rds         : Clean individual CAFI records with functional groups
#   - coral_master.rds       : Coral characteristics with all CAFI metrics
#   - community_matrix.rds   : Coral x OTU abundance matrix
#   - condition_scores.rds   : Position-corrected coral condition (PC1)
#   - functional_summary.rds : CAFI classified by functional group
#
# NOTES:
#   - All corals are Pocillopora spp. (morphotypes not reliable species IDs)
#   - CAFI "species" are morphological OTUs (no genetic confirmation)
#   - Sampling position bias corrected for physiological metrics
#   - Sites: HAU (Hauru), MAT (Maatea), MRB (Moorea Barrier Reef)
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-04
# ============================================================================

cat("\n========================================\n")
cat("Loading and Cleaning Survey Data\n")
cat("========================================\n\n")

# Load setup (paths, packages, theme)
source(here::here("scripts/00_setup.R"))

# ============================================================================
# 1. LOAD RAW DATA FILES
# ============================================================================

cat("1. Loading raw data files...\n")

cafi_raw <- read_csv(here("data/survey_cafi_data_w_taxonomy_summer2019_v5.csv"),
                     show_col_types = FALSE) %>%
  clean_names()

coral_raw <- read_csv(here("data/survey_coral_characteristics_merged_v2.csv"),
                      show_col_types = FALSE) %>%
  clean_names()

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
    type = tolower(type),

    # Create OTU identifier (morphological groupings)
    # Handle case where search_term column may not exist
    otu = case_when(
      "search_term" %in% names(.) & !is.na(search_term) & search_term != "" ~ str_trim(search_term),
      "code" %in% names(.) & !is.na(code) ~ paste(type, code, sep = "_"),
      TRUE ~ type
    ),

    species   = otu,
    size_mm   = if ("cafi_size_mm" %in% names(.)) as.numeric(cafi_size_mm) else NA_real_,
    abundance = 1,

    # Extract site from coral_id prefix
    site = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ NA_character_
    ),

    # Assign functional group
    functional_group = map2_chr(otu, type, assign_functional_group)
  ) %>%
  filter(!is.na(site))

cat("  - Clean CAFI records:", nrow(cafi_clean), "\n")
cat("  - Unique OTUs:", n_distinct(cafi_clean$otu), "\n")
cat("  - Unique corals:", n_distinct(cafi_clean$coral_id), "\n")

# Summarize CAFI by coral
cafi_by_coral <- cafi_clean %>%
  group_by(coral_id, site) %>%
  summarise(
    total_cafi     = n(),
    otu_richness   = n_distinct(otu),
    shannon        = vegan::diversity(table(otu)),

    # Functional group counts
    n_trapezia      = sum(functional_group == "Trapezia"),
    n_resident_fish = sum(functional_group == "Resident Fish"),
    n_corallivore   = sum(functional_group == "Corallivore"),
    n_other_crab    = sum(functional_group == "Other Crab"),
    n_shrimp        = sum(functional_group == "Shrimp"),
    n_other         = sum(functional_group == "Other"),

    # Taxonomic group counts
    n_crabs   = sum(type == "crab"),
    n_shrimps = sum(type == "shrimp"),
    n_fish    = sum(type == "fish"),
    n_snails  = sum(type == "snail"),

    mean_cafi_size = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  )

cat("  - Sites:", paste(unique(cafi_clean$site), collapse = ", "), "\n\n")

# ============================================================================
# 3. CLEAN CORAL DATA
# ============================================================================

cat("3. Cleaning coral characteristics...\n")

coral_clean <- coral_raw %>%
  distinct(coral_id, .keep_all = TRUE) %>%
  mutate(
    morphotype_original = tolower(morphotype),
    morphotype_display  = "Pocillopora spp.",
    morphotype          = morphotype_original,

    branch_width = if ("branch_width" %in% names(.)) tolower(branch_width) else NA_character_,
    depth_m      = as.numeric(depth),
    year         = year(mdy(date)),
    month        = month(mdy(date)),
    subsite      = site,

    # Extract site from coral_id
    site = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ NA_character_
    ),

    # Coral volume
    volume     = coalesce(volume_lab, volume_field, length_lab * width_lab * height_lab),
    log_volume = log10(volume + 1),

    # Dimensions
    height = coalesce(height_lab, height_field),
    width  = coalesce(width_lab, width_field),
    length = coalesce(length_lab, length_field),

    # Sampling position
    stump_length  = coalesce(nub1_stump_length, nub2_stump_length),
    nubbin_length = coalesce(nub1_nubbin_length, nub2_nubbin_length),

    # Neighborhood metrics
    n_neighbors          = number_of_neighbors,
    mean_neighbor_dist   = mean_neighbor_distance,
    total_neighbor_volume = combined_total_volume_of_neighbors,
    mean_neighbor_volume  = mean_total_volume_of_neighbors
  ) %>%
  filter(!is.na(site), !is.na(volume), volume > 0)

cat("  - Clean coral records:", nrow(coral_clean), "\n")
cat("  - Branch widths:", paste(na.omit(unique(coral_clean$branch_width)), collapse = ", "), "\n")
cat("  - Sites:", paste(unique(coral_clean$site), collapse = ", "), "\n\n")

# ============================================================================
# 4. CLEAN PHYSIOLOGY DATA
# ============================================================================

cat("4. Processing physiology data...\n")

physio_cols <- names(physio_raw)
metric_cols <- physio_cols[!physio_cols %in% c("coral_id", "site", "date", "notes")]

physio_clean <- physio_raw %>%
  mutate(across(all_of(metric_cols), as.numeric)) %>%
  filter(if_any(all_of(metric_cols), ~!is.na(.))) %>%
  mutate(
    site = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ NA_character_
    )
  )

cat("  - Clean records:", nrow(physio_clean), "\n")
cat("  - Physiology metrics:", length(metric_cols), "\n\n")

# ============================================================================
# 5. CREATE MASTER CORAL DATASET
# ============================================================================

cat("5. Creating master coral dataset...\n")

coral_master <- coral_clean %>%
  left_join(cafi_by_coral, by = c("coral_id", "site")) %>%
  mutate(
    across(starts_with("n_") | starts_with("total_") | starts_with("otu_"),
           ~replace_na(., 0)),
    shannon = replace_na(shannon, 0),

    size_class = cut(volume,
                     breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE),

    isolation_index = mean_neighbor_dist / (volume^(1/3) + 1),
    crowding_index  = total_neighbor_volume / (mean_neighbor_dist + 1),
    relative_size   = volume / (mean_neighbor_volume + 1),
    cafi_density    = total_cafi / volume * 1000,

    species_richness  = otu_richness,
    shannon_diversity = shannon
  )

cat("  - Master dataset:", nrow(coral_master), "rows,", ncol(coral_master), "columns\n\n")

# ============================================================================
# 6. POSITION-CORRECTED CORAL CONDITION SCORES
# ============================================================================

cat("6. Computing position-corrected condition scores...\n")

physio_vars <- c("protein_mg_cm2", "carb_mg_cm2", "zoox_cells_cm2", "afdw_mg_cm2")

physio_merged <- physio_raw %>%
  left_join(coral_clean %>% dplyr::select(coral_id, site, volume, stump_length),
            by = "coral_id") %>%
  filter(!is.na(stump_length))

if (nrow(physio_merged) > 20) {

  correction_data <- physio_merged %>%
    dplyr::select(coral_id, site, volume, stump_length, any_of(physio_vars)) %>%
    drop_na()

  corrected_traits <- correction_data %>%
    dplyr::select(coral_id, site, volume, stump_length)

  for (var in physio_vars) {
    if (var %in% names(correction_data) && sum(!is.na(correction_data[[var]])) > 20) {
      model <- lm(reformulate("stump_length", response = var), data = correction_data)
      resid_z <- scale(residuals(model))[, 1]

      var_name <- case_when(
        var == "protein_mg_cm2"  ~ "protein_corr",
        var == "carb_mg_cm2"     ~ "carb_corr",
        var == "zoox_cells_cm2"  ~ "zoox_corr",
        var == "afdw_mg_cm2"     ~ "afdw_corr"
      )
      corrected_traits[[var_name]] <- resid_z
    }
  }

  corrected_vars <- c("protein_corr", "carb_corr", "zoox_corr", "afdw_corr")
  corrected_vars <- corrected_vars[corrected_vars %in% names(corrected_traits)]

  if (length(corrected_vars) >= 3) {
    pca_data <- corrected_traits %>%
      dplyr::select(all_of(corrected_vars)) %>%
      drop_na()

    if (nrow(pca_data) > 10) {
      pca_result <- prcomp(pca_data, scale. = FALSE, center = TRUE)

      loadings <- pca_result$rotation[, 1]
      condition_pc1 <- if (sum(loadings < 0) > sum(loadings > 0)) {
        -pca_result$x[, 1]
      } else {
        pca_result$x[, 1]
      }

      var_explained <- summary(pca_result)$importance[2, 1] * 100

      corrected_traits$condition_score <- NA
      corrected_traits$condition_score[complete.cases(pca_data)] <- condition_pc1

      cat("  - Position correction applied to", length(corrected_vars), "traits\n")
      cat("  - Condition score (PC1) explains", round(var_explained, 1), "% variance\n")
    }
  }

  condition_scores <- corrected_traits %>%
    dplyr::select(coral_id, site, condition_score, any_of(corrected_vars))

  save_object(condition_scores, "condition_scores.rds")

  coral_master <- coral_master %>%
    left_join(condition_scores %>% dplyr::select(coral_id, condition_score),
              by = "coral_id")

} else {
  cat("  - Insufficient physiology data for position correction\n")
  coral_master$condition_score <- NA
  condition_scores <- tibble(coral_id = character(), condition_score = numeric())
}

cat("  - Corals with condition scores:", sum(!is.na(coral_master$condition_score)), "\n\n")

# ============================================================================
# 7. CREATE COMMUNITY MATRIX
# ============================================================================

cat("7. Creating community matrix (coral x OTU)...\n")

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

cat("  - Matrix dimensions:", nrow(community_matrix), "x", ncol(community_matrix), "\n\n")

# ============================================================================
# 8. CREATE FUNCTIONAL GROUP SUMMARY
# ============================================================================

cat("8. Creating functional group summary...\n")

functional_summary <- cafi_clean %>%
  group_by(functional_group) %>%
  summarise(
    n_individuals = n(),
    n_otus        = n_distinct(otu),
    n_corals      = n_distinct(coral_id),
    mean_size_mm  = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(proportion = n_individuals / sum(n_individuals)) %>%
  arrange(desc(n_individuals))

cat("  - Functional groups:\n")
for (i in 1:nrow(functional_summary)) {
  cat(sprintf("    %s: %d individuals (%.1f%%)\n",
              functional_summary$functional_group[i],
              functional_summary$n_individuals[i],
              functional_summary$proportion[i] * 100))
}

# ============================================================================
# 9. SAVE ALL PROCESSED DATA
# ============================================================================

cat("\n9. Saving processed data objects...\n")

save_object(cafi_clean, "cafi_clean.rds")
save_object(coral_master, "coral_master.rds")
save_object(community_matrix, "community_matrix.rds")
save_object(functional_summary, "functional_summary.rds")
save_object(cafi_by_coral, "cafi_by_coral.rds")

saveRDS(coral_clean, file.path(OBJECTS_DIR, "coral_clean.rds"))
saveRDS(physio_clean, file.path(OBJECTS_DIR, "physio_clean.rds"))

save_table(functional_summary, "functional_group_summary.csv")

otu_summary <- cafi_clean %>%
  group_by(otu, type, functional_group) %>%
  summarise(
    n_individuals = n(),
    n_corals      = n_distinct(coral_id),
    mean_size_mm  = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_individuals))

save_table(otu_summary, "otu_abundance_ranking.csv")

# ============================================================================
# DATA SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("Data Loading Complete\n")
cat("========================================\n\n")

cat("Dataset Overview:\n")
cat("  - Study: CAFI Survey Mo'orea Summer 2019\n")
cat("  - Coral species: Pocillopora spp.\n")
cat("  - Sites: HAU, MAT, MRB\n\n")

cat("Sample Sizes:\n")
cat("  - Total corals:", nrow(coral_master), "\n")
cat("  - Total CAFI individuals:", nrow(cafi_clean), "\n")
cat("  - Unique OTUs:", n_distinct(cafi_clean$otu), "\n")
cat("  - Corals with condition scores:", sum(!is.na(coral_master$condition_score)), "\n\n")

cat("CAFI by Site:\n")
site_cafi <- coral_master %>%
  group_by(site) %>%
  summarise(
    n_corals       = n(),
    total_cafi     = sum(total_cafi),
    mean_cafi      = sum(total_cafi),  # Total for display purposes
    mean_richness  = round(mean(otu_richness), 1),
    .groups = "drop"
  )
print(site_cafi)

cat("\nData files saved to:", OBJECTS_DIR, "\n")
cat("Tables saved to:", TABLES_DIR, "\n\n")
cat("Ready for analysis!\n\n")
