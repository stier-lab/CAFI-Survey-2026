# ============================================================================
# 01_load_data.R - Load and Clean Survey Data
# ============================================================================
#
# PURPOSE: Load raw data files, clean, and create analysis-ready objects
#
# USAGE:
#   source("scripts/00_setup.R")  # Must run first
#   source("scripts/01_load_data.R")
#
# INPUT FILES (data/):
#   - survey_cafi_data_w_taxonomy_summer2019_v5.csv
#   - survey_coral_characteristics_merged_v2.csv
#   - survey_master_phys_data_v3.csv
#
# OUTPUTS (output/objects/):
#   - cafi_clean.rds        : Individual CAFI records with functional groups
#   - coral_master.rds      : Coral characteristics merged with CAFI metrics
#   - community_matrix.rds  : Coral x OTU abundance matrix
#   - condition_scores.rds  : Position-corrected coral condition
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI SURVEY - Load & Clean Data\n")
cat("============================================================\n\n")

# Check that setup has been run
if (!exists("PATHS")) {
  source(here::here("scripts/00_setup.R"))
}

# ============================================================================
# FUNCTIONAL GROUP ASSIGNMENT
# ============================================================================

assign_functional_group <- function(otu, type) {
  otu_lower <- tolower(otu)
  type_lower <- tolower(type)

  case_when(
    # Mutualist defenders (Trapezia crabs)
    str_detect(otu_lower, "trapezia") ~ "Trapezia",

    # Gastropods (snails)
    str_detect(otu_lower, "coralliophila|drupella") ~ "Gastropod",
    type_lower == "snail" ~ "Gastropod",

    # Resident fish (gobies, hawkfish, damsels)
    str_detect(otu_lower, "paragobiodon|gobiodon|caracanthus|dascyllus") ~ "Resident Fish",
    type_lower == "fish" ~ "Resident Fish",

    # Shrimp (Alpheus and others)
    str_detect(otu_lower, "alpheus") ~ "Shrimp",
    type_lower == "shrimp" ~ "Shrimp",

    # Other crabs
    type_lower == "crab" ~ "Other Crab",

    # Everything else
    TRUE ~ "Other"
  )
}

# ============================================================================
# 1. LOAD RAW DATA
# ============================================================================

cat("1. Loading raw data...\n")

cafi_raw <- read_csv(PATHS$cafi_data, show_col_types = FALSE) %>%
  clean_names()

coral_raw <- read_csv(PATHS$coral_data, show_col_types = FALSE) %>%
  clean_names()

physio_raw <- read_csv(PATHS$physio_data, show_col_types = FALSE) %>%
  clean_names()

cat("   CAFI records:", nrow(cafi_raw), "\n")
cat("   Coral records:", nrow(coral_raw), "\n")
cat("   Physiology records:", nrow(physio_raw), "\n\n")

# ============================================================================
# 2. CLEAN CAFI DATA
# ============================================================================

cat("2. Cleaning CAFI data...\n")

# Pre-check column existence (can't use names(.) inside case_when/mutate)
has_search_term <- "search_term" %in% names(cafi_raw)
has_code <- "code" %in% names(cafi_raw)
has_cafi_size <- "cafi_size_mm" %in% names(cafi_raw)

cafi_clean <- cafi_raw %>%
  filter(!is.na(coral_id)) %>%
  mutate(
    type = tolower(type),

    # Create OTU identifier
    otu = case_when(
      has_search_term & !is.na(search_term) & search_term != "" ~ str_trim(search_term),
      has_code & !is.na(code) ~ paste(type, code, sep = "_"),
      TRUE ~ type
    ),

    species = otu,
    size_mm = if (has_cafi_size) as.numeric(cafi_size_mm) else NA_real_,
    abundance = 1,

    # Extract site from coral_id
    site = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ NA_character_
    ),

    # Functional group
    functional_group = map2_chr(otu, type, assign_functional_group)
  )

# Check for unmatched coral_ids before filtering
n_no_site <- sum(is.na(cafi_clean$site))
if (n_no_site > 0) {
  warning(sprintf("DATA LOSS: %d CAFI records dropped (coral_id did not match MRB/HAU/MAT pattern)",
                  n_no_site))
  cat("   WARNING:", n_no_site, "records dropped (unrecognized coral_id pattern)\n")
}
cafi_clean <- cafi_clean %>% filter(!is.na(site))

cat("   Clean records:", nrow(cafi_clean), "\n")
cat("   Unique OTUs:", n_distinct(cafi_clean$otu), "\n")
cat("   Unique corals:", n_distinct(cafi_clean$coral_id), "\n\n")

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
    n_corallivore   = sum(functional_group == "Gastropod"),
    n_galeropsis    = sum(str_detect(tolower(otu), "galeropsis"), na.rm = TRUE),
    n_other_crab    = sum(functional_group == "Other Crab"),
    n_shrimp        = sum(functional_group == "Shrimp"),
    n_other         = sum(functional_group == "Other"),

    # Taxonomic group counts
    n_crabs   = sum(type == "crab"),
    n_shrimps = sum(type == "shrimp"),
    n_fish    = sum(type == "fish"),
    n_snails  = sum(type == "snail"),

    mean_cafi_size = if (all(is.na(size_mm))) NA_real_ else mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================================
# 3. CLEAN CORAL DATA
# ============================================================================

cat("3. Cleaning coral data...\n")

coral_clean <- coral_raw %>%
  distinct(coral_id, .keep_all = TRUE) %>%
  mutate(
    morphotype_original = tolower(morphotype),
    morphotype_display  = "Pocillopora spp.",

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

    # Survey type: "neighborhood" corals (POC01-21) had full 5m surveys
    #              "size" corals (POC22+) were surveyed only for size/CAFI
    survey_type = if ("survey_type" %in% names(.)) survey_type else NA_character_,

    # Volume
    volume     = coalesce(volume_lab, volume_field, length_lab * width_lab * height_lab),
    log_volume = log10(volume),  # volume > 0 guaranteed by filter

    # Dimensions
    height = coalesce(height_lab, height_field),
    width  = coalesce(width_lab, width_field),
    length = coalesce(length_lab, length_field),

    # Sampling position
    stump_length  = coalesce(nub1_stump_length, nub2_stump_length),
    nubbin_length = coalesce(nub1_nubbin_length, nub2_nubbin_length),

    # Neighborhood metrics (only valid for survey_type == "neighborhood")
    n_neighbors           = number_of_neighbors,
    mean_neighbor_dist    = mean_neighbor_distance,
    total_neighbor_volume = combined_total_volume_of_neighbors,
    mean_neighbor_volume  = mean_total_volume_of_neighbors
  ) %>%
  filter(!is.na(site), !is.na(volume), volume > 0)

cat("   Clean records:", nrow(coral_clean), "\n")
cat("   Sites:", paste(unique(coral_clean$site), collapse = ", "), "\n")
cat("   Survey types:\n")
cat("      neighborhood (5m surveys):", sum(coral_clean$survey_type == "neighborhood", na.rm = TRUE), "\n")
cat("      size (CAFI only):", sum(coral_clean$survey_type == "size", na.rm = TRUE), "\n\n")

# ============================================================================
# 4. CREATE MASTER DATASET
# ============================================================================

cat("4. Creating master coral dataset...\n")

coral_master <- coral_clean %>%
  left_join(cafi_by_coral, by = c("coral_id", "site")) %>%
  mutate(
    # Fill NAs with 0 for CAFI counts only
    # IMPORTANT: Preserve NA for neighborhood metrics (n_neighbors, total_neighbor_volume, etc.)
    # because size-only corals (survey_type == "size") were not surveyed for neighbors
    across(c(n_trapezia, n_resident_fish, n_corallivore, n_other_crab, n_shrimp, n_other,
             n_crabs, n_shrimps, n_fish, n_snails, total_cafi, otu_richness),
           ~replace_na(., 0)),
    shannon = replace_na(shannon, 0),

    # Size class
    size_class = cut(volume,
                     breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE),

    # Derived indices
    isolation_index = mean_neighbor_dist / (volume^(1/3) + 1),
    crowding_index  = total_neighbor_volume / (mean_neighbor_dist + 1),
    relative_size   = volume / (mean_neighbor_volume + 1),
    cafi_density    = total_cafi / volume * 1000,

    # Alias names
    species_richness  = otu_richness,
    shannon_diversity = shannon
  )

# Log size class boundaries for reproducibility
size_breaks <- quantile(coral_master$volume, c(0, 0.33, 0.67, 1), na.rm = TRUE)
cat("   Size class breaks (volume cm³):", paste(round(size_breaks, 1), collapse = " | "), "\n")
cat("   Size class distribution:", paste(table(coral_master$size_class), collapse = "/"), "\n")
cat("   Master dataset:", nrow(coral_master), "corals,", ncol(coral_master), "variables\n\n")

# --- Data validation after join ---
n_expected <- nrow(coral_clean)
n_got <- nrow(coral_master)
if (n_got != n_expected) {
  warning(sprintf("JOIN VALIDATION: Expected %d rows after left_join but got %d (possible duplicate keys in cafi_by_coral)",
                  n_expected, n_got))
}
n_dup <- sum(duplicated(coral_master$coral_id))
if (n_dup > 0) {
  warning(sprintf("JOIN VALIDATION: %d duplicate coral_ids detected after join", n_dup))
}
cat("   Validation: row count", ifelse(n_got == n_expected, "OK", "MISMATCH"),
    "| duplicates:", n_dup, "\n\n")

# ============================================================================
# 5. POSITION-CORRECTED CONDITION SCORES
# ============================================================================

cat("5. Computing position-corrected condition scores...\n")

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

      cat("   Position correction applied\n")
      cat("   Condition score (PC1):", round(var_explained, 1), "% variance\n")
    }
  }

  condition_scores <- corrected_traits %>%
    dplyr::select(coral_id, site, condition_score, any_of(corrected_vars))

  save_object(condition_scores, "condition_scores")

  coral_master <- coral_master %>%
    left_join(condition_scores %>% dplyr::select(coral_id, condition_score),
              by = "coral_id")

} else {
  cat("   Insufficient physiology data\n")
  coral_master$condition_score <- NA
  condition_scores <- tibble(coral_id = character(), condition_score = numeric())
}

cat("   Corals with condition:", sum(!is.na(coral_master$condition_score)), "\n\n")

# --- Data validation after condition join ---
n_after_cond <- nrow(coral_master)
if (n_after_cond != n_expected) {
  warning(sprintf("CONDITION JOIN: Row count changed from %d to %d (duplicate coral_ids in condition_scores?)",
                  n_expected, n_after_cond))
}

# ============================================================================
# 6. CREATE COMMUNITY MATRIX
# ============================================================================

cat("6. Creating community matrix...\n")

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

# Add zero-CAFI corals (present in coral_master but absent from cafi_clean)
missing_corals <- setdiff(coral_master$coral_id, rownames(community_matrix))
if (length(missing_corals) > 0) {
  zero_rows <- matrix(0, nrow = length(missing_corals), ncol = ncol(community_matrix),
                      dimnames = list(missing_corals, colnames(community_matrix)))
  community_matrix <- rbind(community_matrix, zero_rows)
  cat("   Added", length(missing_corals), "zero-CAFI corals to community matrix\n")
}

cat("   Matrix:", nrow(community_matrix), "corals x", ncol(community_matrix), "OTUs\n\n")

# ============================================================================
# 6b. CREATE PCA-BASED CAFI COMMUNITY SCORE (PC1_CAFI)
# ============================================================================
# Following Stier et al. manuscript approach: Use PCA of community composition
# to create a composite CAFI score representing overall community structure

cat("6b. Computing PCA-based CAFI community score (PC1_CAFI)...\n")

# Use Hellinger transformation for community data (recommended for PCA)
# This down-weights rare species and handles zeros appropriately
if (nrow(community_matrix) >= 10 && ncol(community_matrix) >= 3) {

  # Hellinger transformation
  comm_hell <- vegan::decostand(community_matrix, method = "hellinger")

  # PCA on Hellinger-transformed community data
  cafi_pca <- prcomp(comm_hell, scale. = FALSE, center = TRUE)

  # Extract PC1 - orient so higher values = more/diverse CAFI
  pc1_scores <- cafi_pca$x[, 1]

  # Check if PC1 should be flipped (make it positively correlated with total abundance)
  total_cafi_vec <- rowSums(community_matrix)
  if (cor(pc1_scores, total_cafi_vec) < 0) {
    pc1_scores <- -pc1_scores
  }

  # Standardize to z-scores
  pc1_cafi <- scale(pc1_scores)[, 1]

  # Calculate variance explained
  var_explained_cafi <- summary(cafi_pca)$importance[2, 1] * 100

  # Create PC1_CAFI data frame
  pc1_cafi_df <- tibble(
    coral_id = rownames(community_matrix),
    pc1_cafi = pc1_cafi
  )

  # Also create PC2_CAFI for additional analyses
  pc2_scores <- cafi_pca$x[, 2]
  pc1_cafi_df$pc2_cafi <- scale(pc2_scores)[, 1]

  cat("   PC1_CAFI computed successfully\n")
  cat("   Variance explained: PC1 =", round(var_explained_cafi, 1), "%,",
      "PC2 =", round(summary(cafi_pca)$importance[2, 2] * 100, 1), "%\n")
  cat("   Correlation with total abundance: r =", round(cor(pc1_scores, total_cafi_vec), 3), "\n")

  # Save PCA results
  cafi_pca_results <- list(
    pca = cafi_pca,
    pc_scores = pc1_cafi_df,
    var_explained = summary(cafi_pca)$importance,
    loadings = cafi_pca$rotation
  )
  save_object(cafi_pca_results, "cafi_pca_results")

  # Merge PC1_CAFI into coral_master
  coral_master <- coral_master %>%
    left_join(pc1_cafi_df, by = "coral_id")

} else {
  cat("   Insufficient data for CAFI PCA\n")
  coral_master$pc1_cafi <- NA
  coral_master$pc2_cafi <- NA
}

cat("\n")

# ============================================================================
# 7. SAVE ALL OBJECTS
# ============================================================================

cat("7. Saving processed data...\n")

save_object(cafi_clean, "cafi_clean")
save_object(coral_master, "coral_master")
save_object(community_matrix, "community_matrix")
save_object(cafi_by_coral, "cafi_by_coral")

# Functional group summary
functional_summary <- cafi_clean %>%
  group_by(functional_group) %>%
  summarise(
    n_individuals = n(),
    n_otus        = n_distinct(otu),
    n_corals      = n_distinct(coral_id),
    proportion    = n() / nrow(cafi_clean),
    .groups = "drop"
  ) %>%
  arrange(desc(n_individuals))

save_object(functional_summary, "functional_summary")
save_table(functional_summary, "functional_group_summary")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    Data Loading Complete\n")
cat("============================================================\n\n")

cat("Study: CAFI Survey, Mo'orea, Summer 2019\n")
cat("Coral species: Pocillopora spp.\n\n")

cat("Sample sizes:\n")
cat("  Corals:", nrow(coral_master), "\n")
cat("  CAFI individuals:", nrow(cafi_clean), "\n")
cat("  Unique OTUs:", n_distinct(cafi_clean$otu), "\n")
cat("  Corals with condition (PC1_Coral):", sum(!is.na(coral_master$condition_score)), "\n")
cat("  Corals with PC1_CAFI:", sum(!is.na(coral_master$pc1_cafi)), "\n\n")

cat("By site:\n")
coral_master %>%
  group_by(site) %>%
  summarise(
    n = n(),
    total_cafi = sum(total_cafi),
    mean_richness = round(mean(otu_richness), 1),
    .groups = "drop"
  ) %>%
  print()

cat("\nObjects saved to:", PATHS$objects, "\n")
cat("\nReady for analysis!\n\n")
