#!/usr/bin/env Rscript
# ============================================================================
# 01_load_clean_data.R - Load and clean Survey field data
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
# ============================================================================

cat("\n========================================\n")
cat("Loading and Cleaning Survey Data\n")
cat("========================================\n\n")

# Load libraries
source(here::here("scripts/00_load_libraries.R"))

# ============================================================================
# Load Raw Data
# ============================================================================

cat("Loading raw data files...\n")

# Load CAFI data
cafi_data <- read_csv(here("data/1. survey_cafi_data_w_taxonomy_summer2019_v5.csv"),
                      show_col_types = FALSE) %>%
  clean_names()

# Load coral characteristics
coral_chars <- read_csv(here("data/1. survey_coral_characteristics_merged_v2.csv"),
                        show_col_types = FALSE) %>%
  clean_names()

# Load physiology data
physio_data <- read_csv(here("data/1. survey_master_phys_data_v3.csv"),
                        show_col_types = FALSE) %>%
  clean_names()

cat("✓ Data files loaded\n")
cat("  - CAFI records:", nrow(cafi_data), "\n")
cat("  - Coral records:", nrow(coral_chars), "\n")
cat("  - Physiology records:", nrow(physio_data), "\n\n")

# ============================================================================
# Clean and Process CAFI Data
# ============================================================================

cat("Processing CAFI data...\n")

# Clean CAFI data
cafi_clean <- cafi_data %>%
  # Remove empty rows
  filter(!is.na(coral_id)) %>%
  # Standardize taxonomic names
  mutate(
    type = tolower(type),
    lowest_level = str_trim(lowest_level),
    search_term = str_trim(search_term),
    # Create CAFI OTU identifier (operational taxonomic unit)
    # Note: Without genetic haplotyping, these are morphological groupings not biological species
    # Use for richness/diversity metrics but not species-level ecological inferences
    # Use search_term (the actual taxon name) preferentially, then fall back to type_code
    species = case_when(
      !is.na(search_term) & search_term != "" ~ search_term,
      !is.na(code) ~ paste(type, code, sep = "_"),
      TRUE ~ type
    ),
    # Convert size to numeric
    size_mm = as.numeric(cafi_size_mm),
    # Add abundance (each row is one individual)
    abundance = 1
  ) %>%
  # Add site information from coral_id
  # Sites: HAU (Hauru), MAT (Maatea), MRB (Moorea Barrier Reef)
  mutate(
    site = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(site != "Unknown")

# Summarize by species and coral
cafi_summary <- cafi_clean %>%
  group_by(coral_id, site, type, species) %>%
  summarise(
    count = n(),
    mean_size_mm = mean(size_mm, na.rm = TRUE),
    min_size_mm = min(size_mm, na.rm = TRUE),
    max_size_mm = max(size_mm, na.rm = TRUE),
    .groups = "drop"
  )

cat("✓ CAFI data processed\n")
cat("  - Clean records:", nrow(cafi_clean), "\n")
cat("  - Unique corals:", n_distinct(cafi_clean$coral_id), "\n")
cat("  - Unique species:", n_distinct(cafi_clean$species), "\n")
cat("  - Sites:", paste(unique(cafi_clean$site), collapse = ", "), "\n\n")

# ============================================================================
# Clean and Process Coral Data
# ============================================================================

cat("Processing coral characteristics...\n")

coral_clean <- coral_chars %>%
  # Remove duplicates
  distinct(coral_id, .keep_all = TRUE) %>%
  # Clean and prepare coral trait variables
  mutate(
    # IMPORTANT: Morphotype names (meandrina, eudoxi, verucosa) are NOT confirmed species
    # They are morphological field categories without genetic verification
    # Keep for legacy compatibility but DO NOT interpret as biological species
    morphotype_original = tolower(morphotype),
    morphotype_display = "Pocillopora spp.",  # Use this for all figures
    morphotype = morphotype_original,  # Legacy - avoid using in new analyses

    # BRANCH_WIDTH is the real measurable trait to use for analyses
    # tight = narrow branching (meandrina/verucosa morphs)
    # wide = broad branching (eudoxi morph)
    branch_width = if("branch_width" %in% names(.)) tolower(branch_width) else NA_character_,

    # Convert depth to numeric
    depth_m = as.numeric(depth),

    # Add survey year/month
    year = year(mdy(date)),
    month = month(mdy(date)),

    # Preserve sub-site information from raw data
    subsite = site,

    # Extract major site from coral_id prefix
    # HAU = Hauru, MAT = Maatea, MRB = Moorea Barrier Reef
    site = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ "Unknown"
    )
  ) %>%
  # Filter to summer 2019 survey
  filter(year == 2019, month %in% c(6, 7, 8), site != "Unknown")

cat("✓ Coral data processed\n")
cat("  - Clean corals:", nrow(coral_clean), "\n")
cat("  - Branch widths (PRIMARY TRAIT):", paste(na.omit(unique(coral_clean$branch_width)), collapse = ", "), "\n")
cat("  - Morphotypes (legacy, not species):", paste(unique(coral_clean$morphotype), collapse = ", "), "\n")
cat("\n")

# ============================================================================
# Clean and Process Physiology Data
# ============================================================================

cat("Processing physiology data...\n")

# Get column names to identify physiology metrics
physio_cols <- names(physio_data)
metric_cols <- physio_cols[!physio_cols %in% c("coral_id", "site", "date", "notes")]

physio_clean <- physio_data %>%
  # Convert metrics to numeric
  mutate(across(all_of(metric_cols), as.numeric)) %>%
  # Remove rows with all NA metrics
  filter(if_any(all_of(metric_cols), ~!is.na(.))) %>%
  # Add site if missing
  # Sites: HAU (Hauru), MAT (Maatea), MRB (Moorea Barrier Reef)
  mutate(
    site = case_when(
      str_detect(coral_id, "^MRB") ~ "MRB",
      str_detect(coral_id, "^HAU") ~ "HAU",
      str_detect(coral_id, "^MAT") ~ "MAT",
      TRUE ~ NA_character_
    )
  )

cat("✓ Physiology data processed\n")
cat("  - Clean records:", nrow(physio_clean), "\n")
cat("  - Physiology metrics:", length(metric_cols), "\n\n")

# ============================================================================
# Merge Datasets
# ============================================================================

cat("Merging datasets...\n")

# Create master dataset
survey_master <- coral_clean %>%
  # Add CAFI community data
  left_join(
    cafi_summary %>%
      pivot_wider(
        id_cols = coral_id,
        names_from = species,
        values_from = count,
        values_fn = sum,
        values_fill = list(count = 0),
        names_prefix = "sp_"
      ),
    by = "coral_id"
  ) %>%
  # Add physiology data
  left_join(
    physio_clean %>%
      select(coral_id, all_of(metric_cols)),
    by = "coral_id"
  )

# Calculate community metrics
community_metrics <- cafi_summary %>%
  group_by(coral_id) %>%
  summarise(
    total_cafi = sum(count),
    species_richness = n_distinct(species),
    shannon_diversity = -sum((count/sum(count)) * log(count/sum(count))),
    mean_cafi_size = weighted.mean(mean_size_mm, count, na.rm = TRUE),
    .groups = "drop"
  )

# Add community metrics to master
survey_master <- survey_master %>%
  left_join(community_metrics, by = "coral_id")

# Consolidate branch_width columns from multiple joins
# Joins can create .x and .y suffixes - consolidate into single column
if("branch_width.x" %in% names(survey_master)) {
  survey_master <- survey_master %>%
    mutate(
      branch_width.x = as.character(branch_width.x),
      branch_width.y = as.character(branch_width.y),
      branch_width = coalesce(branch_width.x, branch_width.y)
    ) %>%
    select(-matches("branch_width\\.[xy]$"))
}

cat("✓ Datasets merged\n")
cat("  - Master dataset rows:", nrow(survey_master), "\n")
cat("  - Master dataset columns:", ncol(survey_master), "\n\n")

# ============================================================================
# Create Analysis Datasets
# ============================================================================

cat("Creating analysis datasets...\n")

# Community matrix (corals x species)
species_cols <- names(survey_master)[str_detect(names(survey_master), "^sp_")]
community_matrix <- survey_master %>%
  select(coral_id, all_of(species_cols)) %>%
  column_to_rownames("coral_id") %>%
  as.matrix()

# Environmental data
env_cols <- c("coral_id", "site", "morphotype", "depth_m", "lat", "long", metric_cols)
if("branch_width" %in% names(survey_master)) env_cols <- c(env_cols[1:3], "branch_width", env_cols[4:length(env_cols)])
env_data <- survey_master %>%
  select(any_of(env_cols))

# Metadata - using any_of() to handle optional columns
meta_cols <- c("coral_id", "site", "subsite", "date", "survey_type", "field_obs",
               "morphotype", "branch_width", "depth_m", "lat", "long")
metadata <- survey_master %>%
  select(any_of(meta_cols))

cat("✓ Analysis datasets created\n")
cat("  - Community matrix:", dim(community_matrix)[1], "x", dim(community_matrix)[2], "\n")
cat("  - Environmental variables:", ncol(env_data) - 1, "\n\n")

# ============================================================================
# Save Processed Data
# ============================================================================

cat("Saving processed data...\n")

# Save as RDS for R analyses
saveRDS(survey_master, file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
saveRDS(cafi_clean, file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
saveRDS(coral_clean, file.path(SURVEY_OBJECTS, "coral_clean.rds"))
saveRDS(physio_clean, file.path(SURVEY_OBJECTS, "physio_clean.rds"))
saveRDS(community_matrix, file.path(SURVEY_OBJECTS, "community_matrix.rds"))
saveRDS(env_data, file.path(SURVEY_OBJECTS, "environmental_data.rds"))
saveRDS(metadata, file.path(SURVEY_OBJECTS, "metadata.rds"))

# Save summary tables as CSV
write_csv(cafi_summary, file.path(SURVEY_TABLES, "cafi_species_summary.csv"))
write_csv(community_metrics, file.path(SURVEY_TABLES, "community_metrics.csv"))

# Save data dictionary
data_dict <- tibble(
  dataset = c("survey_master", "cafi_clean", "coral_clean", "physio_clean",
              "community_matrix", "env_data", "metadata"),
  description = c(
    "Complete merged dataset with all variables",
    "Clean CAFI observations with taxonomy",
    "Clean coral characteristics",
    "Clean physiology measurements",
    "Corals x species abundance matrix",
    "Environmental and physiology variables",
    "Survey metadata"
  ),
  n_rows = c(nrow(survey_master), nrow(cafi_clean), nrow(coral_clean),
             nrow(physio_clean), nrow(community_matrix), nrow(env_data),
             nrow(metadata)),
  n_cols = c(ncol(survey_master), ncol(cafi_clean), ncol(coral_clean),
             ncol(physio_clean), ncol(community_matrix), ncol(env_data),
             ncol(metadata))
)

write_csv(data_dict, file.path(SURVEY_TABLES, "data_dictionary.csv"))

cat("✓ Data saved successfully\n\n")

# ============================================================================
# Summary Statistics
# ============================================================================

cat("========================================\n")
cat("Data Loading Summary\n")
cat("========================================\n\n")

cat("Survey Coverage:\n")
cat("  - Total corals surveyed:", n_distinct(survey_master$coral_id), "\n")
cat("  - Total CAFI observed:", sum(cafi_clean$abundance), "\n")
cat("  - Species richness:", n_distinct(cafi_clean$species), "\n")
cat("  - Sites:", n_distinct(survey_master$site), "\n")
cat("  - Date range:", min(coral_clean$date), "to", max(coral_clean$date), "\n\n")

cat("Taxonomic Breakdown:\n")
type_summary <- cafi_clean %>%
  group_by(type) %>%
  summarise(
    n_individuals = n(),
    n_species = n_distinct(species),
    .groups = "drop"
  ) %>%
  arrange(desc(n_individuals))

print(type_summary)

cat("\n✅ Data loading and cleaning complete!\n")
cat("All processed data saved to:", SURVEY_OBJECTS, "\n")