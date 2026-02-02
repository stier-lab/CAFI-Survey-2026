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
    log_volume = log(volume),  # natural log; volume > 0 guaranteed by filter

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
# Position correction removes a SAMPLING ARTIFACT caused by tissue availability:
#
# PROBLEM: Coral tissue shows a tip-to-base gradient (tips have lower tissue
# density than bases). Smaller corals had less tissue available, so we had to
# sample more of the branch (longer nubbins) to get enough material. This means
# small-coral samples include more low-density tip tissue, creating a spurious
# negative correlation between coral size and tissue metrics (especially AFDW).
#
# EVIDENCE: nubbin_length predicts AFDW even after controlling for volume
# (partial r = -0.29), confirming this is a methodological artifact, not biology.
#
# SOLUTION: Regress each physiological trait on stump_length + nubbin_length,
# then use residuals. This removes the sampling position effect while preserving
# true biological variation in coral condition.
#
# Key correlations (raw data):
#   AFDW vs volume: r = -0.33 (artifact - disappears after correction)
#   AFDW vs nubbin_length: r = -0.41 (tissue gradient effect)
#   stump vs nubbin: r = 0.25 (low collinearity, VIF = 1.06)
# ============================================================================

cat("5. Computing position-corrected condition scores...\n")

physio_vars <- c("protein_mg_cm2", "carb_mg_cm2", "zoox_cells_cm2", "afdw_mg_cm2")

physio_merged <- physio_raw %>%
  left_join(coral_clean %>% dplyr::select(coral_id, site, volume, stump_length, nubbin_length),
            by = "coral_id") %>%
  filter(!is.na(stump_length) | !is.na(nubbin_length))

if (nrow(physio_merged) > 20) {

  correction_data <- physio_merged %>%
    dplyr::select(coral_id, site, volume, stump_length, nubbin_length, any_of(physio_vars)) %>%
    drop_na()

  corrected_traits <- correction_data %>%
    dplyr::select(coral_id, site, volume, stump_length, nubbin_length)

  cat("   Correcting for: stump_length + nubbin_length (distance from tip effects)\n")

  for (var in physio_vars) {
    if (var %in% names(correction_data) && sum(!is.na(correction_data[[var]])) > 20) {
      # Correct for BOTH stump_length and nubbin_length
      model <- lm(reformulate(c("stump_length", "nubbin_length"), response = var),
                  data = correction_data)
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

# ############################################################################
# MANUSCRIPT FIGURE 1: STUDY DESIGN
# ############################################################################

# Set PROJ library for sf (MUST be before loading sf/maptiles)
# Auto-detect PROJ_LIB from sf package installation
proj_path <- tryCatch(
  system.file("proj", package = "sf"),
  error = function(e) ""
)
if (nzchar(proj_path)) Sys.setenv(PROJ_LIB = proj_path)

tryCatch({

  # --- Check for required spatial packages ---
  has_sf       <- requireNamespace("sf",       quietly = TRUE)
  has_maptiles <- requireNamespace("maptiles", quietly = TRUE)

  if (!has_sf || !has_maptiles) {
    warning("Figure 1 skipped: packages 'sf' and/or 'maptiles' not available. ",
            "Install with install.packages(c('sf','maptiles')) to enable map generation.")
  } else {

    library(sf)
    library(maptiles)

    cat("\n========================================\n")
    cat("Creating Figure 1: Study Design\n")
    cat("========================================\n\n")

    # ------------------------------------------------------------------
    # Site coordinates
    # ------------------------------------------------------------------
    site_coords <- tibble(
      site = c("HAU", "MAT", "MRB"),
      site_name = c("Hauru", "Maatea", "Maharepa"),
      lat = c(-17.516, -17.604, -17.475),
      long = c(-149.922, -149.815, -149.817),
      habitat = c("North Shore", "East Shore", "Barrier Reef")
    )

    site_counts <- coral_master %>%
      group_by(site) %>%
      summarise(
        n_corals = n(),
        total_cafi = sum(total_cafi, na.rm = TRUE),
        .groups = "drop"
      )

    site_data <- site_coords %>%
      left_join(site_counts, by = "site")

    # ------------------------------------------------------------------
    # Color palette (Okabe-Ito, colorblind-safe)
    # ------------------------------------------------------------------
    SITE_COLORS_FIG1 <- c(
      "HAU" = "#E69F00",
      "MAT" = "#0072B2",
      "MRB" = "#009E73"
    )

    DATA_FILL       <- "#4A90A4"
    DATA_FILL_ALPHA  <- "#7EB3C4"
    DENSITY_LINE     <- "#1A3A5C"

    # ------------------------------------------------------------------
    # Theme for Figure 1
    # ------------------------------------------------------------------
    theme_fig1 <- function(base_size = 10) {
      theme_minimal(base_size = base_size) +
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
          panel.grid.minor = element_blank(),
          axis.line  = element_line(color = "gray40", linewidth = 0.4),
          axis.ticks = element_line(color = "gray40", linewidth = 0.3),
          axis.text  = element_text(color = "gray20", size = base_size - 1),
          axis.title = element_text(color = "gray10", size = base_size, face = "plain"),
          plot.title = element_text(face = "bold", size = 15, color = "black",
                                    hjust = 0, margin = margin(b = 2)),
          plot.subtitle = element_text(size = base_size - 0.5, color = "gray40",
                                       hjust = 0, margin = margin(b = 8)),
          plot.margin = margin(10, 12, 8, 10)
        )
    }

    # ==================================================================
    # PANEL A: Satellite map of Mo'orea with study sites
    # ==================================================================
    cat("Creating Panel A: Satellite map with study sites...\n")

    bbox <- st_bbox(c(xmin = -149.98, ymin = -17.65, xmax = -149.70, ymax = -17.40),
                    crs = st_crs(4326))
    bbox_sf <- st_as_sfc(bbox)

    cat("  Fetching satellite imagery...\n")
    sat_tiles <- tryCatch({
      get_tiles(bbox_sf, provider = "Esri.WorldImagery", zoom = 12, crop = TRUE)
    }, error = function(e) {
      cat("  WARNING: Could not fetch satellite tiles:", e$message, "\n")
      NULL
    })

    sites_sf <- st_as_sf(site_data, coords = c("long", "lat"), crs = 4326)

    label_positions <- tibble(
      site = c("HAU", "MAT", "MRB"),
      x_off = c(0.025, 0.038, 0.035),
      y_off_name = c(0.005, 0.005, 0.005),
      y_off_n = c(-0.012, 0.020, -0.012)
    )

    site_data <- site_data %>%
      left_join(label_positions, by = "site")

    if (!is.null(sat_tiles)) {
      panel_a <- ggplot() +
        tidyterra::geom_spatraster_rgb(data = sat_tiles) +
        geom_sf(data = sites_sf, aes(fill = site),
                shape = 21, size = 6, color = "white", stroke = 2.5) +
        geom_label(data = site_data,
                   aes(x = long + x_off, y = lat + y_off_name,
                       label = paste0(site_name, "\n(n=", n_corals, ")")),
                   size = 3.0, fontface = "bold", color = "white",
                   fill = alpha("gray15", 0.85), linewidth = 0.3,
                   label.r = unit(0.15, "lines"),
                   label.padding = unit(0.22, "lines"), lineheight = 0.85) +
        scale_fill_manual(values = SITE_COLORS_FIG1, guide = "none") +
        coord_sf(xlim = c(-149.97, -149.71), ylim = c(-17.65, -17.40), expand = FALSE) +
        annotate("segment", x = -149.949, xend = -149.859, y = -17.621, yend = -17.621,
                 color = "black", linewidth = 2.0, alpha = 0.4) +
        annotate("segment", x = -149.95, xend = -149.86, y = -17.62, yend = -17.62,
                 color = "white", linewidth = 1.5) +
        annotate("segment", x = -149.95, xend = -149.95, y = -17.628, yend = -17.612,
                 color = "white", linewidth = 1.2) +
        annotate("segment", x = -149.86, xend = -149.86, y = -17.628, yend = -17.612,
                 color = "white", linewidth = 1.2) +
        annotate("text", x = -149.904, y = -17.601, label = "5 km",
                 size = 3.0, color = "black", fontface = "bold", alpha = 0.5) +
        annotate("text", x = -149.905, y = -17.602, label = "5 km",
                 size = 3.0, color = "white", fontface = "bold") +
        annotate("text", x = -149.734, y = -17.424, label = "N",
                 size = 4.5, fontface = "bold", color = "black", alpha = 0.4) +
        annotate("text", x = -149.735, y = -17.425, label = "N",
                 size = 4.5, fontface = "bold", color = "white") +
        annotate("segment", x = -149.735, xend = -149.735, y = -17.46, yend = -17.435,
                 arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
                 color = "white", linewidth = 1.0) +
        labs(title = "A") +
        theme_void() +
        theme(
          plot.title = element_text(face = "bold", size = 15, hjust = 0,
                                    margin = margin(b = 5, l = 5)),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(8, 5, 5, 5)
        )

    } else {
      cat("  Using schematic map (satellite unavailable)...\n")

      moorea_coords <- matrix(c(
        -149.75, -17.46, -149.77, -17.44, -149.80, -17.44, -149.84, -17.46,
        -149.88, -17.48, -149.91, -17.51, -149.92, -17.54, -149.90, -17.58,
        -149.86, -17.60, -149.81, -17.60, -149.78, -17.57, -149.76, -17.53,
        -149.75, -17.49, -149.75, -17.46
      ), ncol = 2, byrow = TRUE)

      moorea_sf <- st_polygon(list(moorea_coords)) %>%
        st_sfc(crs = 4326) %>%
        st_sf(name = "Mo'orea")

      panel_a <- ggplot() +
        annotate("rect", xmin = -149.98, xmax = -149.70, ymin = -17.65, ymax = -17.40,
                 fill = "#E8F4F8", color = NA) +
        geom_sf(data = moorea_sf, fill = "#2D7D7D", color = "#1A5050", linewidth = 0.6) +
        geom_sf(data = sites_sf, aes(fill = site),
                shape = 21, size = 6, color = "white", stroke = 2) +
        geom_label(data = site_data,
                   aes(x = long + x_off, y = lat + y_off_name, label = site_name),
                   size = 3.0, fontface = "bold", color = "gray10",
                   fill = alpha("white", 0.85), label.size = 0) +
        geom_text(data = site_data,
                  aes(x = long + x_off, y = lat + y_off_n,
                      label = paste0("(n=", n_corals, ")")),
                  size = 2.5, color = "gray35") +
        scale_fill_manual(values = SITE_COLORS_FIG1, guide = "none") +
        coord_sf(xlim = c(-149.97, -149.72), ylim = c(-17.64, -17.41), expand = FALSE) +
        annotate("segment", x = -149.95, xend = -149.86, y = -17.62, yend = -17.62,
                 color = "gray30", linewidth = 0.8) +
        annotate("text", x = -149.905, y = -17.605, label = "5 km",
                 size = 2.5, color = "gray30") +
        labs(title = "A") +
        theme_void() +
        theme(
          plot.title = element_text(face = "bold", size = 15, hjust = 0,
                                    margin = margin(b = 5, l = 5)),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(8, 5, 5, 5)
        )
    }

    # ==================================================================
    # PANEL B: Colony size distribution
    # ==================================================================
    cat("Creating Panel B: Colony size distribution...\n")

    vol_data <- coral_master %>%
      filter(!is.na(volume), volume > 0)

    vol_stats <- vol_data %>%
      summarise(
        n = n(),
        min_vol = min(volume),
        max_vol = max(volume),
        median_vol = median(volume),
        mean_vol = mean(volume),
        cv = sd(volume) / mean(volume) * 100,
        range_orders = log10(max(volume)) - log10(min(volume))
      )

    cat("  Volume range:", round(vol_stats$min_vol), "to", round(vol_stats$max_vol), "cm3\n")
    cat("  Range spans", round(vol_stats$range_orders, 1), "orders of magnitude\n")

    dens_obj <- density(log10(vol_data$volume))
    dens_max_b <- max(dens_obj$y)

    panel_b <- vol_data %>%
      ggplot(aes(x = volume)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 18, fill = DATA_FILL, color = "white",
                     alpha = 0.85, linewidth = 0.3) +
      geom_density(color = DENSITY_LINE, linewidth = 1.1, adjust = 1.1) +
      scale_x_log10(
        labels = scales::label_comma(),
        breaks = c(100, 1000, 10000, 100000),
        limits = c(15, 180000)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      annotate("text", x = sqrt(vol_stats$min_vol * vol_stats$max_vol),
               y = dens_max_b * 1.15,
               label = ">3 orders\nof magnitude",
               size = 2.8, color = "gray35", fontface = "italic", lineheight = 0.9) +
      labs(
        title = "B",
        subtitle = paste0("n=", vol_stats$n, "  |  CV=", round(vol_stats$cv), "%"),
        x = expression("Colony Volume (cm"^3*")"),
        y = "Density"
      ) +
      theme_fig1() +
      theme(
        axis.title.x = element_text(margin = margin(t = 8)),
        panel.grid.major.x = element_blank()
      )

    # ==================================================================
    # PANEL C: Neighborhood density distribution
    # ==================================================================
    cat("Creating Panel C: Neighborhood density distribution...\n")

    neighbor_data <- coral_master %>%
      filter(!is.na(n_neighbors))

    neighbor_stats <- neighbor_data %>%
      summarise(
        n = n(),
        min_n = min(n_neighbors),
        max_n = max(n_neighbors),
        median_n = median(n_neighbors),
        mean_n = mean(n_neighbors),
        cv = sd(n_neighbors) / mean(n_neighbors) * 100
      )

    cat("  Neighbors range:", neighbor_stats$min_n, "to", neighbor_stats$max_n, "\n")
    cat("  Median:", neighbor_stats$median_n, "\n")

    dens_max_c <- max(density(neighbor_data$n_neighbors)$y)

    # Compute density for D/E/F markers on Panel C
    dens_obj <- density(neighbor_data$n_neighbors, adjust = 1.0)
    get_density_at_x <- function(x, dens) approx(dens$x, dens$y, xout = x)$y

    marker_df <- data.frame(
      x = c(5, 17, 76),
      label = c("D", "E", "F"),
      color = c(SITE_COLORS_FIG1["HAU"], SITE_COLORS_FIG1["MRB"], SITE_COLORS_FIG1["MRB"])
    )
    marker_df$y <- sapply(marker_df$x, get_density_at_x, dens = dens_obj)

    panel_c <- neighbor_data %>%
      ggplot(aes(x = n_neighbors)) +
      geom_histogram(aes(y = after_stat(density)),
                     binwidth = 5, fill = DATA_FILL, color = "white",
                     alpha = 0.85, linewidth = 0.3) +
      geom_density(color = DENSITY_LINE, linewidth = 1.1, adjust = 1.0) +
      # D/E/F markers as inverted triangles
      geom_point(data = marker_df, aes(x = x, y = y + 0.004, fill = color),
                 shape = 25, color = "gray30", size = 2.8, stroke = 0.6) +
      scale_fill_identity() +
      geom_text(data = marker_df, aes(x = x, y = y + 0.0085, label = label),
                fontface = "bold", size = 3.2, color = "gray20") +
      scale_x_continuous(breaks = seq(0, 80, by = 20), limits = c(-2, 88)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(
        title = "C",
        subtitle = paste0("n=", neighbor_stats$n, "  |  CV=", round(neighbor_stats$cv), "%"),
        x = "Neighbors within 5 m",
        y = "Density"
      ) +
      theme_fig1() +
      theme(
        axis.title.x = element_text(margin = margin(t = 8)),
        panel.grid.major.x = element_blank()
      )

    # ==================================================================
    # PANELS D, E, F: Neighborhood schematics
    # ==================================================================
    cat("Creating Panels D-F: Neighborhood schematics...\n")

    # Scale volume to point size (cube root for area perception)
    vol_to_size <- function(vol, min_size = 1.5, max_size = 14) {
      scaled <- (vol^(1/3) - 20^(1/3)) / (45000^(1/3) - 20^(1/3))
      pmax(min_size, pmin(max_size, min_size + scaled * (max_size - min_size)))
    }

    # Neighborhood schematic function
    create_neighborhood_plot <- function(focal_vol, mean_neigh_vol, n_neighbors,
                                         mean_neigh_dist_cm, site_name, focal_color,
                                         panel_label, density_label,
                                         focal_size_override = NULL) {
      set.seed(42 + n_neighbors)
      mean_dist_scaled <- mean_neigh_dist_cm / 500
      focal_size <- if (!is.null(focal_size_override)) focal_size_override else
        vol_to_size(focal_vol, min_size = 5, max_size = 16)

      if (n_neighbors > 0) {
        if (n_neighbors <= 10) {
          base_angles <- seq(0, 2*pi, length.out = n_neighbors + 1)[1:n_neighbors]
          jitter <- runif(n_neighbors, -pi/n_neighbors * 0.4, pi/n_neighbors * 0.4)
          angles <- base_angles + jitter
        } else {
          angles <- runif(n_neighbors, 0, 2*pi)
        }
        radii_raw <- rnorm(n_neighbors, mean = mean_dist_scaled, sd = mean_dist_scaled * 0.3)
        radii <- pmax(0.15, pmin(0.95, radii_raw))
        neighbor_vols <- rlnorm(n_neighbors, meanlog = log(mean_neigh_vol) - 0.125, sdlog = 0.5)
        neighbor_vols <- pmax(50, pmin(10000, neighbor_vols))
        neighbor_sizes <- vol_to_size(neighbor_vols, min_size = 1.2, max_size = 6)
        neighbor_df <- data.frame(x = radii * cos(angles), y = radii * sin(angles),
                                  size = neighbor_sizes)
      } else {
        neighbor_df <- data.frame(x = numeric(0), y = numeric(0), size = numeric(0))
      }

      ggplot() +
        annotate("path", x = cos(seq(0, 2*pi, length.out = 100)),
                 y = sin(seq(0, 2*pi, length.out = 100)),
                 color = "gray75", linetype = "dashed", linewidth = 0.4) +
        geom_point(data = neighbor_df, aes(x = x, y = y, size = size),
                   shape = 21, fill = "gray55", color = "gray35", stroke = 0.4) +
        scale_size_identity() +
        annotate("point", x = 0, y = 0, shape = 21, fill = focal_color,
                 color = "white", size = focal_size, stroke = 2) +
        annotate("text", x = 0, y = 1.38, label = density_label,
                 fontface = "bold", size = 3.8, color = "gray25") +
        annotate("text", x = 0, y = 1.22, label = paste0(n_neighbors, " neighbors"),
                 size = 3.2, color = "gray45") +
        annotate("text", x = 0, y = -1.28, label = site_name,
                 fontface = "bold", size = 4.2, color = focal_color) +
        labs(tag = panel_label) +
        coord_fixed(xlim = c(-1.2, 1.2), ylim = c(-1.45, 1.55)) +
        theme_void() +
        theme(plot.tag = element_text(face = "bold", size = 14),
              plot.margin = margin(5, 8, 5, 8),
              plot.background = element_rect(fill = "white", color = NA))
    }

    # Fixed focal coral size across D-F (panels illustrate density, not coral size)
    FOCAL_SIZE <- 10

    # D: HAU-POC04 (Hauru) - 5 neighbors, far neighbors
    panel_d <- create_neighborhood_plot(
      focal_vol = 26741, mean_neigh_vol = 1336, n_neighbors = 5,
      mean_neigh_dist_cm = 167, site_name = "Hauru",
      focal_color = SITE_COLORS_FIG1["HAU"], panel_label = "D",
      density_label = "Low density", focal_size_override = FOCAL_SIZE)

    # E: MRB-POC10 (Maharepa) - 17 neighbors
    panel_e <- create_neighborhood_plot(
      focal_vol = 5472, mean_neigh_vol = 686, n_neighbors = 17,
      mean_neigh_dist_cm = 154, site_name = "Maharepa",
      focal_color = SITE_COLORS_FIG1["MRB"], panel_label = "E",
      density_label = "Median density", focal_size_override = FOCAL_SIZE)

    # F: MRB-POC18 (Maharepa) - 76 neighbors, close neighbors
    panel_f <- create_neighborhood_plot(
      focal_vol = 3064, mean_neigh_vol = 395, n_neighbors = 76,
      mean_neigh_dist_cm = 113, site_name = "Maharepa",
      focal_color = SITE_COLORS_FIG1["MRB"], panel_label = "F",
      density_label = "High density", focal_size_override = FOCAL_SIZE)

    # ==================================================================
    # Combine panels and save
    # ==================================================================
    cat("Combining panels...\n")

    # Top row: map + histograms; Bottom row: neighborhood schematics
    top_row <- panel_a + panel_b + panel_c +
      plot_layout(ncol = 3, widths = c(1.35, 1, 1))

    bottom_row <- panel_d + panel_e + panel_f +
      plot_layout(ncol = 3, widths = c(1, 1, 1))

    fig1 <- top_row / bottom_row +
      plot_layout(heights = c(1.3, 1)) +
      plot_annotation(
        caption = paste0(
          "(A) Study sites on Mo\u2019orea, French Polynesia. ",
          "(B) Colony volume distribution spanning >3 orders of magnitude.\n",
          "(C) Neighborhood density (Pocillopora within 5 m); triangles mark example corals in D\u2013F. ",
          "(D\u2013F) Focal coral (colored)\nwith neighbors (gray) within 5 m; ",
          "circle sizes reflect colony volumes, positions reflect measured distances."),
        theme = theme(
          plot.background = element_rect(fill = "white", color = NA),
          plot.caption = element_text(size = 7.5, hjust = 0, color = "gray45",
                                     lineheight = 1.3, margin = margin(t = 15)),
          plot.margin = margin(12, 15, 15, 15)
        )
      )

    # Save to manuscript directory
    output_path_ms <- file.path(PATHS$fig_manuscript, "fig1_study_design.png")
    ggsave(output_path_ms, fig1, width = 12, height = 10, dpi = 300, bg = "white")
    cat("  Saved:", output_path_ms, "\n")

    # Save to analysis directory
    output_path_data <- file.path(PATHS$fig_01_data, "fig1_study_design.png")
    ggsave(output_path_data, fig1, width = 12, height = 10, dpi = 300, bg = "white")
    cat("  Saved:", output_path_data, "\n")

    # ------------------------------------------------------------------
    # Write legend, methods, and results text file
    # ------------------------------------------------------------------
    legend_path <- file.path(PATHS$fig_manuscript, "fig1_legend_results.txt")

    legend_text <- paste0(
'FIGURE 1: STUDY DESIGN AND SAMPLING OVERVIEW
================================================================================

FIGURE LEGEND
-------------
Figure 1. Study design and sampling overview. (A) Satellite imagery of Mo\'orea,
French Polynesia (17\u00b030\'S, 149\u00b050\'W), showing the three reef sites: Hauru
(fringing reef, north shore; n = ', site_data$n_corals[site_data$site == "HAU"], ' corals), Maatea (lagoon/back reef, east
shore; n = ', site_data$n_corals[site_data$site == "MAT"], '), and Maharepa (barrier reef, north shore; n = ', site_data$n_corals[site_data$site == "MRB"], '). Site markers
colored by reef site (orange = Hauru, blue = Maatea, green = Maharepa).
(B) Distribution of Pocillopora colony volumes on a log10 scale, spanning >3
orders of magnitude (', round(vol_stats$min_vol), '\u2013', format(round(vol_stats$max_vol), big.mark = ","), ' cm3; CV = ', round(vol_stats$cv), '%; n = ', vol_stats$n, ' colonies with volume
data). Black curve shows kernel density estimate. (C) Distribution of
neighborhood density (number of Pocillopora colonies within 5 m; n = ', neighbor_stats$n, ' corals
with neighborhood surveys). Inverted triangles mark the three example colonies
shown in panels D\u2013F. (D\u2013F) Neighborhood schematics for three representative
colonies at low (D; 5 neighbors), median (E; 17 neighbors), and high (F; 76
neighbors) density. Colored central point = focal coral; gray points = neighbors
within 5 m. Point sizes for neighbors reflect relative colony volumes; positions
approximate measured distances from the focal coral. Dashed circle = 5 m survey
radius. Focal coral size is held constant across panels to emphasize the density
gradient.

================================================================================

METHODS
-------

Study System
------------
We surveyed 114 Pocillopora coral colonies across three reef habitats on
Mo\'orea, French Polynesia during summer 2019. ', 114 - vol_stats$n, ' colonies were excluded from
volume analyses due to missing measurements (final n = ', vol_stats$n, '). Sites were chosen
to span the range of reef environments around Mo\'orea:

  - Hauru (HAU): Fringing reef on the north shore. Characterized by moderate
    wave exposure and relatively shallow reef flat.
  - Maatea (MAT): Lagoon/back-reef habitat on the east shore. Sheltered
    environment with calmer water conditions.
  - Maharepa (MRB): Barrier reef on the north shore. Exposed to oceanic swell,
    higher wave energy, and deeper reef structure.

At each site, 35\u201339 Pocillopora colonies were haphazardly selected across a
range of sizes. Colony size was estimated as ellipsoidal volume from field
measurements of length, width, and height (V = 4/3 * pi * a * b * c, where
a, b, c are the semi-axes).

Neighborhood Surveys
--------------------
For a subset of ', neighbor_stats$n, ' corals (across all three sites), we counted and measured all
Pocillopora neighbors within a 5 m radius of the focal colony. The 5 m radius
was chosen to capture the local spatial context relevant to larval settlement
and post-settlement interactions.

Neighborhood density (n_neighbors) ranged from ', neighbor_stats$min_n, ' to ', neighbor_stats$max_n, ' colonies within 5 m
(median = ', neighbor_stats$median_n, ', mean = ', round(neighbor_stats$mean_n, 1), ', CV = ', round(neighbor_stats$cv), '%).

Example Colonies (Panels D\u2013F)
-----------------------------
Three colonies were selected to illustrate the range of neighborhood densities:

  Panel D \u2014 HAU-POC04 (Hauru):
    Volume: 22,852 cm3 | Neighbors: 5 | Mean neighbor distance: 167 cm
    Represents the low-density extreme (sparse, isolated colony).

  Panel E \u2014 MRB-POC10 (Maharepa):
    Volume: 5,472 cm3 | Neighbors: 17 | Mean neighbor distance: 154 cm
    Represents the median neighborhood density.

  Panel F \u2014 MRB-POC18 (Maharepa):
    Volume: 3,079 cm3 | Neighbors: 76 | Mean neighbor distance: 113 cm
    Represents the high-density extreme (densely packed neighborhood).

Focal coral point size is held constant across panels D\u2013F to emphasize the
density contrast rather than differences in focal colony volume. Neighbor
volumes were drawn from a lognormal distribution parameterized by measured mean
neighbor volume at each site.

================================================================================

RESULTS
-------

Colony Size Distribution (Panel B)
-----------------------------------
Pocillopora colony volumes spanned more than 3 orders of magnitude
(', round(vol_stats$min_vol), '\u2013', format(round(vol_stats$max_vol), big.mark = ","), ' cm3; n = ', vol_stats$n, '). The distribution was right-skewed on a linear scale
(mean = ', format(round(vol_stats$mean_vol), big.mark = ","), ' cm3, median = ', format(round(vol_stats$median_vol), big.mark = ","), ' cm3) but approximately normal on a log10
scale, consistent with lognormal colony-size distributions commonly observed
in coral populations. The coefficient of variation (CV = ', round(vol_stats$cv), '%) confirms
substantial size heterogeneity among surveyed colonies.

Neighborhood Density Distribution (Panel C)
--------------------------------------------
Among the ', neighbor_stats$n, ' corals with neighborhood data, density varied widely (', neighbor_stats$min_n, '\u2013', neighbor_stats$max_n, '
neighbors within 5 m; CV = ', round(neighbor_stats$cv), '%). The distribution was right-skewed, with most
corals having 5\u201325 neighbors and a long tail extending to ', neighbor_stats$max_n, '.

Site Allocation
---------------
  HAU (Hauru):     ', site_data$n_corals[site_data$site == "HAU"], ' corals
  MAT (Maatea):    ', site_data$n_corals[site_data$site == "MAT"], ' corals
  MRB (Maharepa):  ', site_data$n_corals[site_data$site == "MRB"], ' corals
  Total:          ', nrow(coral_master), ' corals

================================================================================

KEY STATISTICS
--------------

Colony Volume:
  n = ', vol_stats$n, '
  Range: ', round(vol_stats$min_vol), ' \u2013 ', format(round(vol_stats$max_vol), big.mark = ","), ' cm3
  Median: ', format(round(vol_stats$median_vol), big.mark = ","), ' cm3
  Mean: ', format(round(vol_stats$mean_vol), big.mark = ","), ' cm3
  CV: ', round(vol_stats$cv), '%
  Log10 range: ', round(vol_stats$range_orders, 1), ' orders of magnitude

Neighborhood Density:
  n = ', neighbor_stats$n, '
  Range: ', neighbor_stats$min_n, ' \u2013 ', neighbor_stats$max_n, ' neighbors
  Median: ', neighbor_stats$median_n, ' neighbors
  Mean: ', round(neighbor_stats$mean_n, 1), ' neighbors
  CV: ', round(neighbor_stats$cv), '%

Location: Mo\'orea, French Polynesia (17\u00b030\'S, 149\u00b050\'W)
Survey period: Summer 2019
Focal species: Pocillopora spp.

================================================================================

COLOR SCHEME
------------
Site palette (Panel A markers, D\u2013F focal corals):
  HAU (Hauru):    #E69F00 (Okabe-Ito orange)
  MAT (Maatea):   #0072B2 (Okabe-Ito blue)
  MRB (Maharepa): #009E73 (Okabe-Ito green)

Data distributions (Panels B\u2013C):
  Histogram fill: #4A90A4 (steel blue)
  Density line:   #1A3A5C (dark navy)

Neighborhood schematics (Panels D\u2013F):
  Focal coral:  Site color with white stroke
  Neighbors:    gray55 fill, gray35 border
  5 m radius:   gray75 dashed circle

================================================================================
Generated: ', Sys.Date(), '
Source script: scripts/01_load_data.R
')

    writeLines(legend_text, legend_path)
    cat("  Saved:", legend_path, "\n")

    cat("\nFigure 1 specifications:\n")
    cat("  - Dimensions: 12 x 10 inches\n")
    cat("  - Resolution: 300 dpi (PNG)\n")
    cat("  - Panels: A (satellite map), B (volume), C (neighbors), D-F (neighborhood schematics)\n")
    cat("  - Site colors: HAU=#E69F00, MAT=#0072B2, MRB=#009E73\n")
    cat("  - Total corals:", nrow(coral_master), "\n")
    cat("  - Corals with neighborhood data:", neighbor_stats$n, "\n")
  }

}, error = function(e) {
  warning("Figure 1 generation failed (non-fatal): ", e$message)
  cat("WARNING: Figure 1 could not be generated:", e$message, "\n")
  cat("  The data loading pipeline continues normally.\n")
})

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
