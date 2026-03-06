# ============================================================================
# 10_feature_engineering.R - Feature Engineering for ML Models
# ============================================================================
#
# PURPOSE: Create optimized feature sets for machine learning models predicting
#          CAFI abundance, richness, and community composition
#
# FEATURE CATEGORIES:
#   A. Coral Morphology Features (log_volume, branch_width, site)
#   B. Condition Features (condition_pc1 - position-corrected)
#   C. Neighborhood/Landscape Features (n_neighbors, distances, derived indices)
#   D. Derived Interaction Features (site_volume, size_class)
#   E. Trait-Based Features (community-weighted means from trait database)
#
# OUTPUTS:
#   - output/objects/feature_matrix.rds (full feature set)
#   - output/objects/feature_matrix_reduced.rds (selected features, low VIF)
#   - output/tables/feature_summary.csv
#   - output/tables/feature_correlations.csv
#   - output/tables/vif_results.csv
#   - output/figures/feature_correlations.png
#   - output/figures/feature_importance.png
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-21
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI SURVEY - Feature Engineering for ML\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

# Load setup
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))

# Additional packages for feature engineering
# Set CRAN mirror for package installation
cran_mirror <- "https://cran.rstudio.com"
if (!require("caret", quietly = TRUE)) install.packages("caret", repos = cran_mirror)
if (!require("Boruta", quietly = TRUE)) install.packages("Boruta", repos = cran_mirror)
if (!require("corrplot", quietly = TRUE)) install.packages("corrplot", repos = cran_mirror)

library(caret)
library(Boruta)
library(corrplot)

# Create output directory for feature engineering figures
FIG_DIR <- file.path(PATHS$figures, "10_features")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("1. Loading data...\n")

# Load processed data objects
coral_master <- load_object("coral_master")
cafi_clean <- load_object("cafi_clean")
community_matrix <- load_object("community_matrix")

# Check if condition_scores exists
condition_path <- file.path(PATHS$objects, "condition_scores.rds")
if (file.exists(condition_path)) {
  condition_scores <- readRDS(condition_path)
  cat("   Condition scores loaded: n =", nrow(condition_scores), "\n")
} else {
  condition_scores <- NULL
  cat("   No condition scores found\n")
}

# Load trait data
trait_path <- here::here("data/traits/cafi_traits_final.csv")
if (file.exists(trait_path)) {
  traits_db <- read_csv(trait_path, show_col_types = FALSE)
  cat("   Trait database loaded:", nrow(traits_db), "taxa\n")
} else {
  traits_db <- NULL
  cat("   No trait database found\n")
}

cat("   Coral master: n =", nrow(coral_master), "corals\n")
cat("   CAFI records: n =", nrow(cafi_clean), "individuals\n\n")

# ============================================================================
# A. CORAL MORPHOLOGY FEATURES
# ============================================================================

cat("2. Creating coral morphology features...\n")

morphology_features <- coral_master %>%
  transmute(
    coral_id = coral_id,

    # Natural log, consistent with core pipeline (scripts 01-09)
    log_volume = log(volume),

    # Raw volume for some uses
    volume = volume,

    # Branch width categorical
    branch_width = case_when(
      is.na(branch_width) ~ "unknown",
      tolower(branch_width) %in% c("tight", "narrow", "fine") ~ "tight",
      tolower(branch_width) %in% c("wide", "broad", "open") ~ "wide",
      TRUE ~ "unknown"
    ) %>% factor(levels = c("tight", "wide", "unknown")),

    # Site factor
    site = factor(site, levels = c("HAU", "MAT", "MRB"))
  )

cat("   log_volume range:", round(range(morphology_features$log_volume, na.rm = TRUE), 2), "\n")
cat("   Branch width distribution:\n")
print(table(morphology_features$branch_width, useNA = "ifany"))
cat("\n")

# ============================================================================
# B. CONDITION FEATURES
# ============================================================================

cat("3. Creating condition features...\n")

if (!is.null(condition_scores) && "condition_score" %in% names(condition_scores)) {
  condition_features <- condition_scores %>%
    dplyr::select(coral_id, condition_pc1 = condition_score) %>%
    filter(!is.na(condition_pc1))

  cat("   Corals with condition data: n =", nrow(condition_features), "\n")
  cat("   condition_pc1 range:", round(range(condition_features$condition_pc1), 2), "\n\n")
} else {
  condition_features <- tibble(coral_id = character(), condition_pc1 = numeric())
  cat("   No condition features available\n\n")
}

# ============================================================================
# C. NEIGHBORHOOD/LANDSCAPE FEATURES
# ============================================================================

cat("4. Creating neighborhood/landscape features...\n")

# Filter to corals with neighborhood data
neighborhood_features <- coral_master %>%
  filter(!is.na(n_neighbors)) %>%
  transmute(
    coral_id = coral_id,

    # Count within 5m
    n_neighbors = n_neighbors,

    # Sum of neighbor volumes
    total_neighbor_volume = total_neighbor_volume,
    log_total_neighbor_volume = log10(total_neighbor_volume + 1),

    # Average distance to neighbors (in cm)
    mean_neighbor_dist = mean_neighbor_dist,

    # Mean neighbor volume
    mean_neighbor_volume = mean_neighbor_volume,

    # Isolation index: Mean distance / focal volume^(1/3)
    # Higher values = more isolated relative to size
    isolation_index = mean_neighbor_dist / (volume^(1/3) + 1),

    # Relative size: Focal volume / mean neighbor volume
    # >1 = focal is larger than average neighbor
    relative_size = volume / (mean_neighbor_volume + 1),
    log_relative_size = log10(relative_size + 0.01),

    # Crowding index: Total neighbor volume / mean distance
    # Higher values = more crowded neighborhood
    crowding_index = total_neighbor_volume / (mean_neighbor_dist + 1),
    log_crowding_index = log10(crowding_index + 1)
  )

cat("   Corals with neighborhood data: n =", nrow(neighborhood_features), "\n")
cat("   n_neighbors range:", range(neighborhood_features$n_neighbors), "\n")
cat("   mean_neighbor_dist range (cm):", round(range(neighborhood_features$mean_neighbor_dist, na.rm = TRUE)), "\n\n")

# ============================================================================
# D. DERIVED INTERACTION FEATURES
# ============================================================================

cat("5. Creating derived interaction features...\n")

derived_features <- coral_master %>%
  transmute(
    coral_id = coral_id,

    # Size class bins (tertiles)
    size_class = cut(volume,
                     breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                     labels = c("small", "medium", "large"),
                     include.lowest = TRUE) %>% factor(),

    # Site-specific volume effect (interaction proxy)
    # Natural log, consistent with core pipeline (scripts 01-09)
    site_volume_HAU = ifelse(site == "HAU", log(volume), 0),
    site_volume_MAT = ifelse(site == "MAT", log(volume), 0),
    site_volume_MRB = ifelse(site == "MRB", log(volume), 0)
  )

cat("   Size class distribution:\n")
print(table(derived_features$size_class, useNA = "ifany"))
cat("\n")

# ============================================================================
# E. TRAIT-BASED FEATURES (Community-Weighted Means)
# ============================================================================

cat("6. Creating trait-based features...\n")

if (!is.null(traits_db)) {

  # Clean trait data
  traits_clean <- traits_db %>%
    mutate(
      # Parse functional role
      functional_role = str_split(functional_role_str, ";\\s*"),
      is_obligate = pocillopora_association == "obligate",
      is_facultative = pocillopora_association == "facultative",

      # Max length (mm)
      max_length_mm = as.numeric(max_length_final_mm),

      # Trophic level (from FishBase)
      trophic_level_num = as.numeric(trophic_level)
    ) %>%
    dplyr::select(
      taxon = taxon_input,
      accepted_name,
      max_length_mm,
      pocillopora_association,
      is_obligate,
      is_facultative,
      functional_role_str,
      trophic_level_num
    )

  # Get unique functional roles
  all_roles <- unique(unlist(str_split(traits_db$functional_role_str, ";\\s*")))
  all_roles <- all_roles[!is.na(all_roles) & all_roles != "unknown"]
  cat("   Functional roles found:", paste(head(all_roles, 5), collapse = ", "), "...\n")

  # Calculate trait-based features for each coral
  trait_features_list <- list()

  for (i in 1:nrow(coral_master)) {
    coral_id <- coral_master$coral_id[i]

    # Get CAFI for this coral
    coral_cafi <- cafi_clean %>%
      filter(coral_id == !!coral_id)

    if (nrow(coral_cafi) == 0) {
      trait_features_list[[i]] <- tibble(
        coral_id = coral_id,
        mean_body_size = NA_real_,
        functional_diversity = NA_real_,
        obligate_ratio = 0,
        trophic_diversity = NA_real_,
        n_unique_roles = 0
      )
      next
    }

    # Match CAFI to traits
    coral_traits <- coral_cafi %>%
      left_join(traits_clean, by = c("otu" = "taxon")) %>%
      filter(!is.na(accepted_name) | !is.na(otu))

    # Calculate community-weighted means
    # Mean body size (average max_length of species present)
    mean_body_size <- mean(coral_traits$max_length_mm, na.rm = TRUE)

    # Obligate ratio (proportion of obligate associates)
    n_total <- nrow(coral_traits)
    n_obligate <- sum(coral_traits$is_obligate, na.rm = TRUE)
    obligate_ratio <- n_obligate / max(n_total, 1)

    # Trophic diversity (SD of trophic levels as simple diversity measure)
    trophic_vals <- coral_traits$trophic_level_num[!is.na(coral_traits$trophic_level_num)]
    trophic_diversity <- if (length(trophic_vals) > 1) sd(trophic_vals) else 0

    # Functional diversity (number of unique functional roles)
    roles_present <- unique(unlist(str_split(coral_traits$functional_role_str, ";\\s*")))
    roles_present <- roles_present[!is.na(roles_present) & roles_present != "unknown"]
    n_unique_roles <- length(roles_present)

    # Simpson diversity of functional roles
    if (n_unique_roles > 1) {
      role_counts <- table(unlist(str_split(coral_traits$functional_role_str, ";\\s*")))
      role_counts <- role_counts[names(role_counts) != "unknown"]
      role_props <- role_counts / sum(role_counts)
      functional_diversity <- 1 - sum(role_props^2)  # Simpson diversity
    } else {
      functional_diversity <- 0
    }

    trait_features_list[[i]] <- tibble(
      coral_id = coral_id,
      mean_body_size = mean_body_size,
      functional_diversity = functional_diversity,
      obligate_ratio = obligate_ratio,
      trophic_diversity = trophic_diversity,
      n_unique_roles = n_unique_roles
    )
  }

  trait_features <- bind_rows(trait_features_list)

  cat("   Trait features calculated for n =", nrow(trait_features), "corals\n")
  cat("   mean_body_size range:", round(range(trait_features$mean_body_size, na.rm = TRUE), 1), "mm\n")
  cat("   obligate_ratio range:", round(range(trait_features$obligate_ratio, na.rm = TRUE), 2), "\n")
  cat("   functional_diversity range:", round(range(trait_features$functional_diversity, na.rm = TRUE), 2), "\n\n")

} else {
  trait_features <- coral_master %>%
    transmute(
      coral_id = coral_id,
      mean_body_size = NA_real_,
      functional_diversity = NA_real_,
      obligate_ratio = NA_real_,
      trophic_diversity = NA_real_,
      n_unique_roles = NA_integer_
    )
  cat("   No trait features available\n\n")
}

# ============================================================================
# COMBINE ALL FEATURES
# ============================================================================

cat("7. Combining all features into feature matrix...\n")

# Add response variables from coral_master
response_vars <- coral_master %>%
  dplyr::select(coral_id, total_cafi, otu_richness, shannon)

# Full feature matrix
feature_matrix <- response_vars %>%
  left_join(morphology_features, by = "coral_id") %>%
  left_join(condition_features, by = "coral_id") %>%
  left_join(neighborhood_features, by = "coral_id") %>%
  left_join(derived_features, by = "coral_id") %>%
  left_join(trait_features, by = "coral_id")

cat("   Full feature matrix: n =", nrow(feature_matrix), "observations,",
    ncol(feature_matrix), "variables\n")

# Summary of missingness
missing_summary <- feature_matrix %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(feature_matrix) * 100, 1)) %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))

cat("\n   Variables with missing values:\n")
print(missing_summary, n = 15)
cat("\n")

# ============================================================================
# HANDLE MISSING VALUES
# ============================================================================

cat("8. Handling missing values...\n")

# Document missing value strategy
missing_strategy <- tibble(
  variable = c("condition_pc1", "neighborhood features", "trait features", "branch_width"),
  n_available = c(
    sum(!is.na(feature_matrix$condition_pc1)),
    sum(!is.na(feature_matrix$n_neighbors)),
    sum(!is.na(feature_matrix$mean_body_size)),
    sum(feature_matrix$branch_width != "unknown", na.rm = TRUE)
  ),
  strategy = c(
    "Only n=84 corals sampled for physiology; exclude from models requiring condition",
    "Only n=63 corals with 5m neighborhood surveys; subset for landscape models",
    "Calculate where CAFI present and matched to traits; impute 0 for empty corals",
    "Categorical: unknown level for missing; or exclude from models using branch_width"
  )
)

cat("   Missing value handling strategy:\n")
print(missing_strategy)
cat("\n")

# For trait features on empty corals, set to sensible defaults
feature_matrix <- feature_matrix %>%
  mutate(
    # Empty corals (0 CAFI) get 0 for diversity metrics
    mean_body_size = ifelse(total_cafi == 0 & is.na(mean_body_size), 0, mean_body_size),
    functional_diversity = ifelse(total_cafi == 0 & is.na(functional_diversity), 0, functional_diversity),
    obligate_ratio = ifelse(total_cafi == 0 & is.na(obligate_ratio), 0, obligate_ratio),
    trophic_diversity = ifelse(total_cafi == 0 & is.na(trophic_diversity), 0, trophic_diversity),
    n_unique_roles = ifelse(total_cafi == 0 & is.na(n_unique_roles), 0, n_unique_roles)
  )

# ============================================================================
# SCALE CONTINUOUS FEATURES
# ============================================================================

cat("9. Scaling continuous features...\n")

# Identify numeric features to scale (exclude responses and IDs)
numeric_features <- feature_matrix %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-total_cafi, -otu_richness, -shannon, -volume) %>%
  names()

cat("   Numeric features to scale:", length(numeric_features), "\n")

# Create scaled version
scaling_params <- list()
feature_matrix_scaled <- feature_matrix

for (feat in numeric_features) {
  vals <- feature_matrix[[feat]]
  if (sum(!is.na(vals)) > 1) {
    mu <- mean(vals, na.rm = TRUE)
    sigma <- sd(vals, na.rm = TRUE)
    if (sigma > 0) {
      feature_matrix_scaled[[paste0(feat, "_scaled")]] <- (vals - mu) / sigma
      scaling_params[[feat]] <- list(mean = mu, sd = sigma)
    }
  }
}

cat("   Scaling parameters saved for", length(scaling_params), "features\n\n")

# ============================================================================
# CREATE DUMMY VARIABLES FOR FACTORS
# ============================================================================

cat("10. Creating dummy variables for factors...\n")

# Site dummies (reference: HAU)
feature_matrix_scaled <- feature_matrix_scaled %>%
  mutate(
    site_MAT = as.numeric(site == "MAT"),
    site_MRB = as.numeric(site == "MRB")
  )

# Size class dummies (reference: small)
feature_matrix_scaled <- feature_matrix_scaled %>%
  mutate(
    size_medium = as.numeric(size_class == "medium"),
    size_large = as.numeric(size_class == "large")
  )

# Branch width dummies (reference: tight)
feature_matrix_scaled <- feature_matrix_scaled %>%
  mutate(
    branch_wide = as.numeric(branch_width == "wide"),
    branch_unknown = as.numeric(branch_width == "unknown")
  )

cat("   Dummy variables created for: site, size_class, branch_width\n\n")

# ============================================================================
# CORRELATION ANALYSIS
# ============================================================================

cat("11. Checking feature correlations...\n")

# Select numeric features for correlation (excluding scaled duplicates)
cor_features <- feature_matrix %>%
  dplyr::select(
    log_volume, n_neighbors, total_neighbor_volume, mean_neighbor_dist,
    isolation_index, relative_size, crowding_index, condition_pc1,
    mean_body_size, functional_diversity, obligate_ratio
  ) %>%
  drop_na()

if (nrow(cor_features) > 10) {
  cor_matrix <- cor(cor_features, use = "pairwise.complete.obs")

  # Find high correlations
  high_cor_pairs <- which(abs(cor_matrix) > 0.8 & cor_matrix < 1, arr.ind = TRUE)

  if (nrow(high_cor_pairs) > 0) {
    high_cor_df <- data.frame(
      var1 = rownames(cor_matrix)[high_cor_pairs[,1]],
      var2 = colnames(cor_matrix)[high_cor_pairs[,2]],
      correlation = cor_matrix[high_cor_pairs]
    ) %>%
      filter(var1 < var2) %>%
      arrange(desc(abs(correlation)))

    cat("   High correlations (|r| > 0.8):\n")
    print(high_cor_df)
    cat("\n")
  } else {
    cat("   No correlations |r| > 0.8 found\n\n")
    high_cor_df <- tibble(var1 = character(), var2 = character(), correlation = numeric())
  }

  # Save correlation matrix
  cor_matrix_df <- as.data.frame(cor_matrix) %>%
    rownames_to_column("variable")

  # Create correlation plot
  png(file.path(FIG_DIR, "feature_correlations.png"), width = 10, height = 10, units = "in", res = 300)
  corrplot(cor_matrix, method = "color", type = "upper",
           order = "hclust", hclust.method = "ward.D2",
           addCoef.col = "black", number.cex = 0.7,
           tl.col = "black", tl.srt = 45,
           title = "Feature Correlation Matrix",
           mar = c(0, 0, 2, 0))
  dev.off()
  cat("   Saved: feature_correlations.png\n\n")

} else {
  cat("   Insufficient data for correlation analysis\n\n")
  cor_matrix_df <- tibble()
  high_cor_df <- tibble()
}

# ============================================================================
# VIF ANALYSIS
# ============================================================================

cat("12. Calculating Variance Inflation Factors (VIF)...\n")

# Prepare data for VIF (complete cases only for key predictors)
vif_data <- feature_matrix %>%
  filter(!is.na(n_neighbors), !is.na(log_volume)) %>%
  dplyr::select(total_cafi, log_volume, n_neighbors,
                total_neighbor_volume, mean_neighbor_dist, site)

if (nrow(vif_data) > 20) {
  # Fit model for VIF
  vif_model <- MASS::glm.nb(total_cafi ~ log_volume + n_neighbors +
                       log(total_neighbor_volume + 1) + mean_neighbor_dist + site,
                      data = vif_data)

  vif_results <- car::vif(vif_model)

  # Format VIF results
  vif_df <- data.frame(
    variable = rownames(vif_results),
    GVIF = vif_results[, "GVIF"],
    Df = vif_results[, "Df"],
    `GVIF^(1/(2*Df))` = vif_results[, "GVIF^(1/(2*Df))"]
  ) %>%
    mutate(
      flag = case_when(
        `GVIF..1..2.Df..` > 5 ~ "EXCLUDE (VIF > 5)",
        `GVIF..1..2.Df..` > 3 ~ "Caution (VIF > 3)",
        TRUE ~ "OK"
      )
    )

  cat("   VIF Results (threshold: VIF > 5 problematic):\n")
  print(vif_df)
  cat("\n")

} else {
  vif_df <- tibble(variable = character(), GVIF = numeric(), flag = character())
  cat("   Insufficient data for VIF analysis\n\n")
}

# ============================================================================
# FEATURE IMPORTANCE (Boruta Algorithm)
# ============================================================================

cat("13. Running Boruta feature importance analysis...\n")

# Prepare data for Boruta (complete cases for key features)
boruta_data <- feature_matrix %>%
  filter(!is.na(log_volume), total_cafi > 0) %>%
  dplyr::select(
    total_cafi,
    log_volume, site,
    n_neighbors, mean_neighbor_dist, total_neighbor_volume,
    mean_body_size, functional_diversity, obligate_ratio
  ) %>%
  drop_na()

if (nrow(boruta_data) > 30) {
  cat("   Running Boruta on n =", nrow(boruta_data), "observations...\n")

  set.seed(42)
  boruta_result <- tryCatch({
    Boruta(total_cafi ~ ., data = boruta_data,
           doTrace = 0, maxRuns = 100)
  }, error = function(e) {
    cat("   Boruta failed:", e$message, "\n")
    NULL
  })

  if (!is.null(boruta_result)) {
    # Get importance
    boruta_importance <- attStats(boruta_result) %>%
      rownames_to_column("feature") %>%
      arrange(desc(meanImp))

    cat("   Feature importance (Boruta):\n")
    print(boruta_importance %>% dplyr::select(feature, meanImp, decision))
    cat("\n")

    # Plot importance
    p_importance <- boruta_importance %>%
      mutate(feature = fct_reorder(feature, meanImp)) %>%
      ggplot(aes(x = meanImp, y = feature, fill = decision)) +
      geom_col(alpha = 0.8) +
      scale_fill_manual(values = c("Confirmed" = "#009E73",
                                    "Tentative" = "#E69F00",
                                    "Rejected" = "#D55E00")) +
      labs(x = "Mean Importance (Boruta Z-score)",
           y = "Feature",
           title = "Feature Importance for CAFI Abundance",
           subtitle = "Boruta algorithm with Random Forest") +
      theme_publication()

    ggsave(file.path(FIG_DIR, "feature_importance.png"), p_importance,
           width = 8, height = 6, dpi = 300, bg = "white")
    cat("   Saved: feature_importance.png\n\n")

    # Selected features
    confirmed_features <- boruta_importance %>%
      filter(decision == "Confirmed") %>%
      pull(feature)

  } else {
    boruta_importance <- tibble()
    confirmed_features <- c()
  }

} else {
  cat("   Insufficient data for Boruta analysis\n\n")
  boruta_importance <- tibble()
  confirmed_features <- c()
}

# ============================================================================
# CREATE REDUCED FEATURE SET
# ============================================================================

cat("14. Creating reduced feature set (low multicollinearity)...\n")

# Select non-redundant features based on VIF and correlation analysis
# Exclude: isolation_index (confounded with volume), crowding_index (collinear)
reduced_features <- c(
  "coral_id",
  # Response variables
  "total_cafi", "otu_richness", "shannon",
  # Core predictors
  "log_volume", "site",
  # Neighborhood (non-redundant set)
  "n_neighbors", "mean_neighbor_dist", "total_neighbor_volume",
  # Condition (when available)
  "condition_pc1",
  # Trait-based
  "mean_body_size", "functional_diversity", "obligate_ratio",
  # Derived
  "size_class"
)

# Add scaled versions of continuous features
reduced_scaled <- c(
  "log_volume_scaled", "n_neighbors_scaled", "mean_neighbor_dist_scaled",
  "total_neighbor_volume_scaled", "mean_body_size_scaled"
)

feature_matrix_reduced <- feature_matrix_scaled %>%
  dplyr::select(any_of(c(reduced_features, reduced_scaled)))

cat("   Reduced feature matrix:", ncol(feature_matrix_reduced), "variables\n")
cat("   Excluded due to collinearity: isolation_index, crowding_index, relative_size\n\n")

# ============================================================================
# FEATURE SUMMARY TABLE
# ============================================================================

cat("15. Creating feature summary table...\n")

feature_summary <- tibble(
  feature = c("log_volume", "n_neighbors", "mean_neighbor_dist",
              "total_neighbor_volume", "condition_pc1",
              "mean_body_size", "functional_diversity", "obligate_ratio"),
  category = c("Morphology", "Neighborhood", "Neighborhood",
               "Neighborhood", "Condition",
               "Trait", "Trait", "Trait"),
  description = c(
    "Natural log-transformed coral volume (cm3)",
    "Count of Pocillopora neighbors within 5m radius",
    "Mean distance to neighboring corals (cm)",
    "Sum of neighbor volumes within 5m (cm3)",
    "Position-corrected physiological condition (PC1)",
    "Community-weighted mean body size (mm)",
    "Simpson diversity of functional roles",
    "Proportion of obligate coral associates"
  ),
  n_available = c(
    sum(!is.na(feature_matrix$log_volume)),
    sum(!is.na(feature_matrix$n_neighbors)),
    sum(!is.na(feature_matrix$mean_neighbor_dist)),
    sum(!is.na(feature_matrix$total_neighbor_volume)),
    sum(!is.na(feature_matrix$condition_pc1)),
    sum(!is.na(feature_matrix$mean_body_size) & feature_matrix$mean_body_size > 0),
    sum(!is.na(feature_matrix$functional_diversity)),
    sum(!is.na(feature_matrix$obligate_ratio))
  ),
  mean_value = c(
    round(mean(feature_matrix$log_volume, na.rm = TRUE), 2),
    round(mean(feature_matrix$n_neighbors, na.rm = TRUE), 1),
    round(mean(feature_matrix$mean_neighbor_dist, na.rm = TRUE), 0),
    round(mean(feature_matrix$total_neighbor_volume, na.rm = TRUE), 0),
    round(mean(feature_matrix$condition_pc1, na.rm = TRUE), 2),
    round(mean(feature_matrix$mean_body_size, na.rm = TRUE), 1),
    round(mean(feature_matrix$functional_diversity, na.rm = TRUE), 2),
    round(mean(feature_matrix$obligate_ratio, na.rm = TRUE), 2)
  ),
  sd_value = c(
    round(sd(feature_matrix$log_volume, na.rm = TRUE), 2),
    round(sd(feature_matrix$n_neighbors, na.rm = TRUE), 1),
    round(sd(feature_matrix$mean_neighbor_dist, na.rm = TRUE), 0),
    round(sd(feature_matrix$total_neighbor_volume, na.rm = TRUE), 0),
    round(sd(feature_matrix$condition_pc1, na.rm = TRUE), 2),
    round(sd(feature_matrix$mean_body_size, na.rm = TRUE), 1),
    round(sd(feature_matrix$functional_diversity, na.rm = TRUE), 2),
    round(sd(feature_matrix$obligate_ratio, na.rm = TRUE), 2)
  )
)

cat("   Feature summary:\n")
print(feature_summary)
cat("\n")

# ============================================================================
# SAVE OUTPUTS
# ============================================================================

cat("16. Saving outputs...\n")

# Objects
save_object(feature_matrix, "feature_matrix")
save_object(feature_matrix_reduced, "feature_matrix_reduced")
save_object(feature_matrix_scaled, "feature_matrix_scaled")
save_object(scaling_params, "feature_scaling_params")

# Tables
save_table(feature_summary, "feature_summary")
if (nrow(cor_matrix_df) > 0) save_table(cor_matrix_df, "feature_correlations")
if (nrow(vif_df) > 0) save_table(vif_df, "vif_results")
if (nrow(boruta_importance) > 0) save_table(boruta_importance, "feature_importance_boruta")
save_table(missing_strategy, "missing_value_strategy")

cat("\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("============================================================\n")
cat("    FEATURE ENGINEERING COMPLETE\n")
cat("============================================================\n\n")

cat("FEATURE MATRIX SUMMARY:\n")
cat("  Total corals:", nrow(feature_matrix), "\n")
cat("  Total features:", ncol(feature_matrix), "\n")
cat("  Reduced features:", ncol(feature_matrix_reduced), "\n\n")

cat("SAMPLE SIZES BY DATA AVAILABILITY:\n")
cat("  All corals (morphology only):", nrow(feature_matrix), "\n")
cat("  With condition data:", sum(!is.na(feature_matrix$condition_pc1)), "\n")
cat("  With neighborhood data:", sum(!is.na(feature_matrix$n_neighbors)), "\n")
cat("  With trait features:", sum(!is.na(feature_matrix$mean_body_size) &
                                    feature_matrix$mean_body_size > 0), "\n\n")

cat("KEY FINDINGS:\n")
if (nrow(high_cor_df) > 0) {
  cat("  High correlations found between:\n")
  for (i in 1:min(nrow(high_cor_df), 3)) {
    cat("    -", high_cor_df$var1[i], "~", high_cor_df$var2[i],
        ": r =", round(high_cor_df$correlation[i], 2), "\n")
  }
}

if (nrow(vif_df) > 0) {
  high_vif <- vif_df %>% filter(grepl("EXCLUDE|Caution", flag))
  if (nrow(high_vif) > 0) {
    cat("\n  High VIF features (consider excluding):\n")
    for (i in 1:nrow(high_vif)) {
      cat("    -", high_vif$variable[i], ": VIF =", round(high_vif$`GVIF..1..2.Df..`[i], 2), "\n")
    }
  }
}

if (length(confirmed_features) > 0) {
  cat("\n  Boruta-confirmed important features:\n")
  cat("   ", paste(confirmed_features, collapse = ", "), "\n")
}

cat("\nOUTPUT FILES:\n")
cat("  Objects:\n")
cat("    - feature_matrix.rds (full)\n")
cat("    - feature_matrix_reduced.rds (selected)\n")
cat("    - feature_matrix_scaled.rds (with scaled versions)\n")
cat("    - feature_scaling_params.rds\n")
cat("  Tables:\n")
cat("    - feature_summary.csv\n")
cat("    - feature_correlations.csv\n")
cat("    - vif_results.csv\n")
cat("    - feature_importance_boruta.csv\n")
cat("  Figures:\n")
cat("    - feature_correlations.png\n")
cat("    - feature_importance.png\n\n")

cat("RECOMMENDED MODEL SPECIFICATIONS:\n")
cat("  Full model (n=63): log_volume + n_neighbors + mean_neighbor_dist +\n")
cat("                     log(total_neighbor_volume) + site\n")
cat("  With traits:       Add mean_body_size + functional_diversity + obligate_ratio\n")
cat("  With condition:    Add condition_pc1 (n=84 subsample)\n\n")

cat("Next: Use feature matrices for ML modeling\n\n")
