# Data Dictionary

## CAFI Survey Analysis - Mo'orea, French Polynesia

**Version**: 1.0
**Date**: November 2025
**Contact**: Adrian Stier (astier@ucsb.edu)

---

## Overview

This data dictionary describes all data files associated with the manuscript "Landscape characteristics structure coral-associated fauna communities and their effects on coral condition."

**Study Location**: Mo'orea, French Polynesia (17°30'S, 149°50'W)
**Survey Period**: June–August 2019
**Host Coral**: *Pocillopora* spp.
**Sample Size**: 114 coral colonies, 3,989 CAFI records, 87 OTUs

---

## Primary Data Files

### 1. survey_cafi_data_w_taxonomy_summer2019_v5.csv

**Description**: Individual coral-associated fauna (CAFI) specimen records with taxonomic identification.

**Rows**: 3,989 individual records
**Location**: `data/`

| Variable | Type | Description | Units/Values | Missing |
|----------|------|-------------|--------------|---------|
| coral_id | character | Unique coral colony identifier | Format: "SITE_###" (e.g., "HAU_001") | None |
| site | character | Study site code | HAU, MAT, MRB | None |
| type | character | Broad taxonomic group | Crab, Shrimp, Fish, Gastropod, Echinoderm, Worm, Bivalve, Other | None |
| genus | character | Genus name | Latin genus name (e.g., "Trapezia") | 12 records |
| species | character | Species epithet or morphospecies code | Latin epithet or "sp." | 45 records |
| size_mm | numeric | Individual body size (carapace width for crabs, total length for shrimp/fish) | millimeters (mm) | 156 records |
| count | integer | Number of individuals of this morphotype | 1+ | None |
| notes | character | Additional identification notes | Free text | Most records |

**Taxonomic Notes**:
- Specimens identified to operational taxonomic units (OTUs) based on morphology
- Genetic confirmation not performed
- "sp." indicates identification to genus only
- Some cryptic species complexes may be present (e.g., *Trapezia* spp.)

---

### 2. survey_coral_characteristics_merged_v2.csv

**Description**: Coral colony morphological characteristics and GPS coordinates.

**Rows**: 114 colonies
**Location**: `data/`

| Variable | Type | Description | Units/Values | Missing |
|----------|------|-------------|--------------|---------|
| coral_id | character | Unique coral colony identifier | Format: "SITE-POC##" (e.g., "HAU-POC01") | None |
| site | character | Study site code | HAU, MAT, MRB | None |
| **survey_type** | character | Type of survey conducted | "neighborhood" or "size" | None |
| latitude | numeric | GPS latitude (WGS84) | Decimal degrees (negative = South) | None |
| longitude | numeric | GPS longitude (WGS84) | Decimal degrees (negative = West) | None |
| depth | numeric | Water depth at colony | meters (m) | None |
| max_diameter | numeric | Maximum colony diameter | centimeters (cm) | None |
| perp_diameter | numeric | Perpendicular diameter | centimeters (cm) | None |
| height | numeric | Colony height | centimeters (cm) | None |
| volume_field | numeric | Volume estimated in field | cubic centimeters (cm³) | 8 records |
| volume_lab | numeric | Volume measured in lab (water displacement) | cubic centimeters (cm³) | 23 records |
| surface_area_cm2 | numeric | Estimated surface area | square centimeters (cm²) | 15 records |
| branch_width | character | Branch architecture classification | "tight" or "wide" | 2 records |
| morphotype | character | Putative species morphotype | "verrucosa", "meandrina", "eydouxi" | 5 records |
| collection_date | date | Date of collection | YYYY-MM-DD | None |
| collector | character | Collector initials | 2-3 letter code | None |

**CRITICAL: Survey Type Distinction**:

The study employed two complementary sampling designs:

| survey_type | Coral IDs | Description | N |
|-------------|-----------|-------------|---|
| **neighborhood** | POC01-POC21 at each site | Full 5m radius neighborhood surveys + size + CAFI | 63 |
| **size** | POC22+ at each site | Size and CAFI measurements only (NO neighborhood census) | 51 |

⚠️ **Important**: Corals with `survey_type == "size"` have **NA values for all neighborhood metrics** (number_of_neighbors, mean_neighbor_distance, etc.). These should NOT be interpreted as "isolated" corals—they simply were not surveyed for neighbors.

**Notes**:
- Volume calculated as: V = (2/3) × π × r₁ × r₂ × h (hemi-ellipsoid)
- Branch width classified based on inter-branch spacing: tight (<10mm), wide (≥10mm)
- GPS accuracy approximately ±3 m
- 2 colonies excluded from final analysis due to incomplete data

---

### 3. survey_master_phys_data_v3.csv

**Description**: Coral physiological measurements from branch tip samples.

**Rows**: 108 samples (subset of colonies with physiology)
**Location**: `data/`

| Variable | Type | Description | Units/Values | Missing |
|----------|------|-------------|--------------|---------|
| coral_id | character | Unique coral colony identifier | Format: "SITE_###" | None |
| sample_id | character | Unique sample identifier | Format: "SITE_###_S#" | None |
| stump_length | numeric | Distance from colony base to sampling point | centimeters (cm) | None |
| protein | numeric | Tissue protein content | mg/cm² | 3 records |
| carbohydrate | numeric | Tissue carbohydrate content | mg/cm² | 3 records |
| lipid | numeric | Tissue lipid content | mg/cm² | 12 records |
| zoox_density | numeric | Symbiodiniaceae (zooxanthellae) density | cells × 10⁶/cm² | 5 records |
| chlorophyll_a | numeric | Chlorophyll a concentration | μg/cm² | 5 records |
| chlorophyll_c | numeric | Chlorophyll c concentration | μg/cm² | 8 records |
| afdw | numeric | Ash-free dry weight (tissue biomass) | mg/cm² | 2 records |

**Position Correction Note**:
Sampling position (stump_length) correlates with colony size (r = 0.565), creating a confound. Position-corrected traits are derived by:
1. Regressing each trait on stump_length
2. Extracting residuals as position-corrected values
3. Standardizing to z-scores

---

## Derived Data Files

### 4. master_analysis_data.csv

**Description**: Comprehensive merged analysis dataset with all calculated metrics.

**Rows**: 114 coral colonies
**Location**: `output/tables/`

#### Core Identifiers

| Variable | Type | Description |
|----------|------|-------------|
| coral_id | character | Unique coral colony identifier |
| site | character | Study site (HAU, MAT, MRB) |

#### Coral Characteristics

| Variable | Type | Description | Units |
|----------|------|-------------|-------|
| volume_cm3 | numeric | Colony volume | cm³ |
| log_volume | numeric | log₁₀(volume_cm3) | — |
| surface_area_cm2 | numeric | Colony surface area | cm² |
| depth | numeric | Water depth | m |
| branch_width | character | Branch architecture | tight/wide |
| latitude | numeric | GPS latitude | degrees |
| longitude | numeric | GPS longitude | degrees |

#### CAFI Community Metrics

| Variable | Type | Description | Units |
|----------|------|-------------|-------|
| cafi_abundance | integer | Total CAFI individuals per coral | count |
| cafi_richness | integer | Number of CAFI species per coral | count |
| cafi_density | numeric | Abundance / volume | individuals/cm³ |
| shannon | numeric | Shannon diversity index (H') | nats |
| simpson | numeric | Simpson diversity index (1-D) | — |
| evenness | numeric | Pielou's evenness (J') | — |

#### Taxonomic Group Abundances

| Variable | Type | Description | Units |
|----------|------|-------------|-------|
| crab_abundance | integer | Total crabs per coral | count |
| shrimp_abundance | integer | Total shrimp per coral | count |
| fish_abundance | integer | Total fish per coral | count |
| gastropod_abundance | integer | Total gastropods per coral | count |
| echinoderm_abundance | integer | Total echinoderms per coral | count |

#### Neighborhood Metrics (5m radius)

| Variable | Type | Description | Units |
|----------|------|-------------|-------|
| n_neighbors | integer | Number of coral colonies within 5m | count |
| mean_neighbor_dist | numeric | Mean distance to neighbors | m |
| total_neighbor_volume | numeric | Sum of neighbor volumes | cm³ |
| isolation_index | numeric | mean_distance / volume^(1/3) | — |
| relative_size | numeric | focal_volume / mean_neighbor_volume | ratio |
| spillover_potential | numeric | neighbor_volume / mean_distance | — |

#### Condition Scores

| Variable | Type | Description | Units |
|----------|------|-------------|-------|
| condition_score | numeric | PC1 of raw physiological traits | — |
| condition_score_corrected | numeric | PC1 of position-corrected traits | — |
| protein_corrected | numeric | Position-corrected protein | z-score |
| carbohydrate_corrected | numeric | Position-corrected carbohydrate | z-score |
| zoox_corrected | numeric | Position-corrected zooxanthellae | z-score |
| afdw_corrected | numeric | Position-corrected AFDW | z-score |

---

### 5. cafi_species_summary.csv

**Description**: Species-level summary statistics.

**Rows**: 87 species (OTUs)
**Location**: `output/tables/`

| Variable | Type | Description |
|----------|------|-------------|
| species | character | Species name (Genus species) |
| type | character | Taxonomic group |
| total_abundance | integer | Total individuals across all corals |
| n_corals | integer | Number of corals where species occurred |
| occupancy | numeric | Proportion of corals occupied (0-1) |
| mean_per_coral | numeric | Mean abundance when present |
| max_per_coral | integer | Maximum abundance on single coral |
| cv | numeric | Coefficient of variation |

---

### 6. alpha_diversity_metrics.csv

**Description**: Per-coral diversity metrics.

**Rows**: 114 corals
**Location**: `output/tables/`

| Variable | Type | Description |
|----------|------|-------------|
| coral_id | character | Coral identifier |
| site | character | Study site |
| richness | integer | Species count (S) |
| abundance | integer | Total individuals (N) |
| shannon | numeric | Shannon entropy (H') |
| simpson | numeric | Simpson index (1-D) |
| inv_simpson | numeric | Inverse Simpson (1/D) |
| evenness | numeric | Pielou's J' |
| chao1 | numeric | Chao1 richness estimator |
| ACE | numeric | ACE richness estimator |

---

### 7. cafi_network_metrics.csv

**Description**: Co-occurrence network summary metrics.

**Rows**: 1 (network-level summary)
**Location**: `output/tables/`

| Variable | Type | Description |
|----------|------|-------------|
| n_nodes | integer | Number of species in network |
| n_edges | integer | Number of significant co-occurrences |
| connectance | numeric | Edge density (realized / possible edges) |
| modularity_Q | numeric | Newman-Girvan modularity |
| n_modules | integer | Number of modules detected |
| transitivity | numeric | Global clustering coefficient |
| mean_degree | numeric | Average edges per node |
| mean_path_length | numeric | Average shortest path |
| diameter | integer | Maximum shortest path |

---

### 8. cafi_keystone_species.csv

**Description**: Network centrality metrics for each species.

**Rows**: 200 species (all in network)
**Location**: `output/tables/`

| Variable | Type | Description |
|----------|------|-------------|
| species | character | Species name |
| degree | integer | Number of co-occurrence partners |
| betweenness | numeric | Betweenness centrality (0-1) |
| closeness | numeric | Closeness centrality |
| eigenvector | numeric | Eigenvector centrality |
| module | integer | Module membership |
| within_module_degree | numeric | Z-score of within-module degree |
| participation_coef | numeric | Among-module connectivity |
| keystone_index | numeric | Composite centrality score |

---

## Data Quality Notes

### Quality Control Steps Applied

1. **Duplicate removal**: Checked for duplicate coral_ids
2. **Range validation**: Verified physiological values within expected ranges
3. **GPS validation**: Confirmed coordinates fall within Mo'orea boundaries
4. **Taxonomy standardization**: Harmonized species names across datasets
5. **Outlier flagging**: Identified statistical outliers (±3 SD) for review

### Known Limitations

- **Taxonomic resolution**: Morphological identification only; cryptic species possible
- **Physiology subset**: Only 108/112 corals have physiological data
- **GPS precision**: ±3 m accuracy may affect fine-scale neighborhood metrics
- **Temporal snapshot**: Single summer survey; seasonal variation not captured

---

## File Formats

- **CSV files**: UTF-8 encoding, comma-separated, double-quoted strings
- **RDS files**: R binary format (version 3)
- **PNG figures**: 300 DPI resolution, RGB color space

---

## Citation

If using these data, please cite:

```
Stier AC, et al. (2026). Landscape characteristics structure coral-associated
fauna communities and their effects on coral condition in Mo'orea, French Polynesia.
Marine Ecology Progress Series [submitted].
```

Data archived at: [Dryad DOI to be assigned]

---

## Contact

**Principal Investigator**: Adrian Stier
**Email**: astier@ucsb.edu
**Institution**: UC Santa Barbara, Department of Ecology, Evolution, and Marine Biology
**Lab Website**: https://stierlab.com

---

*Data Dictionary v1.1 - December 2025*
