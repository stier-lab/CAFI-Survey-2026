# ============================================================================
# path_config.R - Global path configuration for CAFI Survey Analysis
# ============================================================================

# This file sets up standardized paths for all analysis scripts

library(here)

# Base output directory
OUTPUT_BASE <- here("output")

# Main subdirectories
FIGURES_DIR <- here("output", "figures")
TABLES_DIR <- here("output", "tables")
OBJECTS_DIR <- here("output", "objects")
REPORTS_DIR <- here("output", "reports")
TRAITS_DIR <- here("output", "traits")

# Data directory
DATA_DIR <- here("data")

# Scripts directory
SCRIPTS_DIR <- here("scripts")

# Create all directories if they don't exist
all_dirs <- c(
  OUTPUT_BASE,
  FIGURES_DIR,
  TABLES_DIR,
  OBJECTS_DIR,
  REPORTS_DIR,
  TRAITS_DIR
)

for (dir in all_dirs) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

cat("Path configuration loaded.\n")
