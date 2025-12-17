# Data Directory

**CAFI Survey Analysis - Mo'orea, French Polynesia**

This directory contains all raw and supporting data for the coral-associated fauna (CAFI) manuscript analysis. This README is the authoritative guide to understanding and working with the data.

---

## Quick Reference

| What you need | File to use |
|--------------|-------------|
| CAFI individual records | `survey_cafi_data_w_taxonomy_summer2019_v5.csv` |
| Coral colony attributes & GPS | `survey_coral_characteristics_merged_v2.csv` |
| Coral physiology measurements | `survey_master_phys_data_v3.csv` |
| Species-association mapping | `association_lookup.csv` |
| Obligate vs facultative classification | `facultative_obligate.csv` |

---

## Primary Data Files

### 1. survey_cafi_data_w_taxonomy_summer2019_v5.csv

**What it is**: Individual CAFI specimen records with full taxonomic classification.

| Property | Value |
|----------|-------|
| Rows | 3,989 (one per individual specimen) |
| Key join column | `coral_id` (e.g., "HAU-POC29") |
| File size | 1.6 MB |

**Key columns you'll use most**:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `coral_id` | text | Unique coral colony ID | "HAU-POC29" |
| `site` | text | Site code (may be blank - use coral_id prefix) | "HAU04" |
| `code` | text | 4-letter species code | "ACHI", "TRCY" |
| `type` | text | Taxonomic group | "crab", "shrimp", "fish" |
| `genus` | text | Genus name | "Trapezia", "Alpheus" |
| `species` | text | Full species name | "Actaeodes hirsutissimus" |
| `cafi_size_mm` | numeric | Body size in mm | 10, 15.5 |
| `measurement` | text | What was measured | "carapace_width", "rostrum_to_tail" |
| `family` | text | Taxonomic family | "Xanthidae", "Alpheidae" |

**Note on site column**: The `site` column is inconsistently populated. Extract site from `coral_id` prefix (HAU, MAT, MRB) instead.

---

### 2. survey_coral_characteristics_merged_v2.csv

**What it is**: Coral colony morphology, GPS coordinates, and pre-calculated neighborhood metrics.

| Property | Value |
|----------|-------|
| Rows | 114 (one per coral colony) |
| Key join column | `coral_id` |
| File size | 37 KB |

**Key columns you'll use most**:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `coral_id` | text | Unique coral colony ID | "MRB-POC01" |
| `site` | text | Site code | "MRB1", "HAU04" |
| `lat` | numeric | Latitude (WGS84) | -17.47428 |
| `long` | numeric | Longitude (WGS84) | -149.81492 |
| `morphotype` | text | Putative species | "meandrina", "eudoxi", "verrucosa" |
| `branch_width` | text | Branch architecture | "tight", "wide" |
| `depth` | numeric | Water depth (m) | 5, 6 |
| `volume_field` | numeric | Field-estimated volume (cm³) | 3539.5 |
| `volume_lab` | numeric | Lab-measured volume (cm³) | NA (often missing) |

**Neighborhood metrics** (pre-calculated within 5m radius):

| Column | Description |
|--------|-------------|
| `number_of_neighbors` | Count of coral colonies within 5m |
| `mean_neighbor_distance` | Mean distance to neighbors (cm) |
| `combined_total_volume_of_neighbors` | Sum of neighbor volumes (cm³) |
| `number_of_wide_branching_neighbors` | Count of wide-branch neighbors |
| `number_of_tight_branching_neighbors` | Count of tight-branch neighbors |

---

### 3. survey_master_phys_data_v3.csv

**What it is**: Coral physiological measurements from branch tip samples.

| Property | Value |
|----------|-------|
| Rows | 108 (subset with physiology data) |
| Key join column | `coral_id` |
| File size | 54 KB |

**Key columns you'll use most**:

| Column | Type | Description | Units |
|--------|------|-------------|-------|
| `coral_id` | text | Unique coral colony ID | — |
| `nub` | integer | Sample nubbin number | 1, 2 |
| `branch_width` | text | Branch architecture | "tight", "wide" |
| `protein_mg_cm2` | numeric | Tissue protein | mg/cm² |
| `carb_mg_cm2` | numeric | Tissue carbohydrate | mg/cm² |
| `zooxDensity` | numeric | Symbiont density | cells/cm² |
| `zoox_cells_cm2` | numeric | Symbiont count per area | cells/cm² |
| `afdw_mg_cm2` | numeric | Ash-free dry weight (biomass) | mg/cm² |
| `surface_area` | numeric | Sample surface area | cm² |

**Important**: Sampling position (`stump_length` - distance from base to sample) correlates with colony size. The analysis scripts apply position correction to physiology metrics.

---

## Lookup Tables

### association_lookup.csv

Maps CAFI species to their ecological association type with corals.

| Column | Description |
|--------|-------------|
| `species` | Species name or code |
| `association_type` | Ecological relationship |

### facultative_obligate.csv / .xls

Classifies species as obligate (coral-dependent) or facultative (can live elsewhere).

---

## Metadata Documentation (Excel/Word files)

These files document the original data collection protocols:

| File | Contents |
|------|----------|
| `README_survey_cafi_metadata_v1.xlsx` | Column definitions for CAFI data |
| `README_survey_coral_characteristics_metadata_v2.xlsx` | Column definitions for coral data |
| `README_survey_physio_metadata_v4.xlsx` | Column definitions for physiology data |
| `README_survey_project_overview.docx` | Project overview and field methods |
| `README_tip_stump_comparison_dec_2019.docx` | Sampling position comparison documentation |

---

## Subdirectories

### gis/

OpenStreetMap shapefiles for Mo'orea island (buildings, landuse, roads, water bodies). Used for map visualizations. ~100 MB total.

**Key files**: `gis_osm_natural_a_free_1.*` (natural features), `gis_osm_water_a_free_1.*` (water bodies)

### traits/

Trait database for CAFI species compiled from multiple sources (WoRMS, FishBase, OBIS).

**Authoritative file**: `cafi_traits_final.csv` - Use this for trait-based analyses.

Other files are intermediate extraction steps (documented in `TRAIT_EXTRACTION_GUIDE.md`).

### projects/

Internal state files for analysis pipeline. Not needed for data analysis.

---

## Working with the Data in R

```r
# Load the three core datasets
library(tidyverse)

# CAFI individual records
cafi <- read_csv("data/survey_cafi_data_w_taxonomy_summer2019_v5.csv")

# Coral colony characteristics
coral <- read_csv("data/survey_coral_characteristics_merged_v2.csv")

# Physiology data
physio <- read_csv("data/survey_master_phys_data_v3.csv")

# Join CAFI to coral (aggregate first)
cafi_summary <- cafi %>%
  group_by(coral_id) %>%
  summarise(
    cafi_abundance = n(),
    cafi_richness = n_distinct(code),
    .groups = "drop"
  )

# Merge with coral characteristics
master <- coral %>%
  left_join(cafi_summary, by = "coral_id") %>%
  left_join(physio %>% select(coral_id, protein_mg_cm2, zooxDensity, afdw_mg_cm2),
            by = "coral_id")
```

**Or use the pre-built pipeline**:
```r
source("scripts/01_load_clean_data.R")
# Creates: cafi_clean, coral_clean, physio_clean, coral_master
# Saves to: output/objects/
```

---

## Data Provenance

| Property | Value |
|----------|-------|
| Collection period | June–August 2019 |
| Location | Mo'orea, French Polynesia (17°30'S, 149°50'W) |
| Sites | HAU (Hauru), MAT (Maatea), MRB (Moorea Barrier Reef) |
| Host coral | *Pocillopora* spp. |
| Sample size | 114 colonies, 3,989 CAFI individuals, 87 OTUs |
| Collectors | Field team initials in data (JC, AP, AS, CO) |

---

## Known Data Issues

1. **Site column in CAFI data**: Often blank. Extract from `coral_id` prefix instead.
2. **Volume measurements**: `volume_lab` often missing. Use `volume_field` or calculate from dimensions.
3. **Morphotype classification**: Putative species assignments—genetic confirmation not performed.
4. **Physiology subset**: Only 108 of 114 corals have physiology data.
5. **GPS precision**: ±3m accuracy affects fine-scale neighborhood metrics.

---

## File Naming Convention

```
survey_[datatype]_[description]_v[version].csv

Examples:
- survey_cafi_data_w_taxonomy_summer2019_v5.csv
- survey_coral_characteristics_merged_v2.csv
- survey_master_phys_data_v3.csv
```

---

## Contact

**PI**: Adrian Stier (astier@ucsb.edu)
**Lab**: Stier Lab, UC Santa Barbara

---

*For detailed variable definitions, see [DATA_DICTIONARY.md](DATA_DICTIONARY.md)*
