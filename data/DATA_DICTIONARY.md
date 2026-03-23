# Data Dictionary

## CAFI Survey Analysis - Mo'orea, French Polynesia

**Version**: 1.2
**Date**: March 2026
**Contact**: Adrian Stier (astier@ucsb.edu)

---

## Overview

This data dictionary describes all data files associated with the manuscript "Sublinear scaling of coral-associated fauna reveals propagule limitation in established reef communities."

**Study Location**: Mo'orea, French Polynesia (17°30'S, 149°50'W)
**Survey Period**: June–August 2019
**Host Coral**: *Pocillopora* spp.
**Sample Size**: 114 coral colonies, 3,989 CAFI records, 243 OTUs

---

## Primary Data Files

### 1. survey_cafi_data_w_taxonomy_summer2019_v5.csv

**Description**: Individual coral-associated fauna (CAFI) specimen records with taxonomic identification.

**Rows**: 3,989 individual records
**Location**: `data/`

| Variable | Type | Description | Units/Values | Missing |
|----------|------|-------------|--------------|---------|
| coral_id | character | Unique coral colony identifier | Format: "SITE-POC##" (e.g., "HAU-POC29") | None |
| site | character | Study site code | HAU, MAT, MRB | None |
| type | character | Broad taxonomic group | Crab, Shrimp, Fish, Gastropod, Echinoderm, Worm, Bivalve, Other | None |
| genus | character | Genus name | Latin genus name (e.g., "Trapezia") | 12 records |
| species | character | Species epithet or morphospecies code | Latin epithet or "sp." | 45 records |
| cafi_size_mm | numeric | Individual body size (carapace width for crabs, total length for shrimp/fish) | millimeters (mm) | 156 records |
| notes | character | Additional identification notes | Free text | Most records |

**Taxonomic Notes**:
- Specimens identified to operational taxonomic units (OTUs) based on morphology
- Genetic confirmation not performed
- "sp." indicates identification to genus only
- Some cryptic species complexes may be present (e.g., *Trapezia* spp.)

---

### 1b. survey_coral_haplotypes_v1.csv

**Description**: Symbiodiniaceae ITS2 haplotype assignments for all 114 survey corals.

**Rows**: 114 colonies
**Location**: `data/`

| Variable | Type | Description | Units/Values | Missing |
|----------|------|-------------|--------------|---------|
| coral_id | character | Unique coral colony identifier | Format: "SITE-POC##" | None |
| site | character | Study site code | HAU, MAT, MRB | None |
| haplotype | character | Symbiodiniaceae ITS2 haplotype | 1a, 1a_Pe, 1a_Pm, 3a, 3b, 3f, 3h, 8a, 10, "Did not amplify", "No sample" | None (but 13 non-valid) |
| dna_collection_date | date | Date tissue sampled for DNA extraction | YYYY-MM-DD | 7 records |

**Haplotype Notes**:
- 101 valid haplotype assignments (88.6% amplification success)
- 11 samples: "Did not amplify" (PCR failure)
- 2 samples: "No sample" (tissue not collected)
- Dominant haplotypes: `1a_Pe` (49, *P. eydouxi*-type) and `1a_Pm` (32, *P. meandrina*-type)
- ITS2 haplotypes identify Symbiodiniaceae clades, not coral host species
- mtORF restriction fragment data available in archived metadata (CAFI_2025 repo)

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
| lat | numeric | GPS latitude (WGS84) | Decimal degrees (negative = South) | None |
| long | numeric | GPS longitude (WGS84) | Decimal degrees (negative = West) | None |
| depth | numeric | Water depth at colony | meters (m) | None |
| length_field | numeric | Maximum colony length (field) | centimeters (cm) | None |
| width_field | numeric | Colony width (field) | centimeters (cm) | None |
| height_field | numeric | Colony height (field) | centimeters (cm) | None |
| volume_field | numeric | Volume estimated in field | cubic centimeters (cm³) | 8 records |
| volume_lab | numeric | Volume measured in lab (water displacement) | cubic centimeters (cm³) | 23 records |
| branch_width | character | Branch architecture classification | "tight" or "wide" | 2 records |
| morphotype | character | Putative species morphotype | "verrucosa", "meandrina", "eydouxi" | 5 records |
| date | date | Date of collection | YYYY-MM-DD | None |
| field_obs | character | Field observer initials | 2-3 letter code | None |

**CRITICAL: Survey Type Distinction**:

The study employed two complementary sampling designs:

| survey_type | Coral IDs | Description | N |
|-------------|-----------|-------------|---|
| **neighborhood** | POC01-POC21 at each site | Full 5m radius neighborhood surveys + size + CAFI | 63 (61 after volume filter) |
| **size** | POC22+ at each site | Size and CAFI measurements only (NO neighborhood census) | 51 |

⚠️ **Important**: Corals with `survey_type == "size"` have **NA values for all neighborhood metrics** (number_of_neighbors, mean_neighbor_distance, etc.). These should NOT be interpreted as "isolated" corals—they simply were not surveyed for neighbors.

**Notes**:
- Volume calculated as: V = (4/3) × π × a × b × c (ellipsoid, where a, b, c are semi-axes)
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
| coral_id | character | Unique coral colony identifier | Format: "SITE-POC##" (e.g., "HAU-POC01") | None |
| nub | character | Nubbin identifier | Nubbin number within coral | None |
| stump_length | numeric | Distance from colony base to sampling point | centimeters (cm) | None |
| protein_mg_cm2 | numeric | Tissue protein content | mg/cm² | 3 records |
| carb_mg_cm2 | numeric | Tissue carbohydrate content | mg/cm² | 3 records |
| zooxDensity | numeric | Symbiodiniaceae (zooxanthellae) density | cells/cm² | 5 records |
| zoox_cells_cm2 | numeric | Zooxanthellae density (alternative column) | cells/cm² | 5 records |
| afdw_mg_cm2 | numeric | Ash-free dry weight (tissue biomass) | mg/cm² | 2 records |

**Position Correction Note**:
Sampling position (stump_length) correlates with colony size (r = 0.565), creating a confound. Position-corrected traits are derived by:
1. Regressing each trait on stump_length + nubbin_length
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
| log_volume | numeric | ln(volume_cm3) [natural log] | — |
| surface_area_cm2 | numeric | Colony surface area | cm² |
| depth | numeric | Water depth | m |
| branch_width | character | Branch architecture | tight/wide |
| latitude | numeric | GPS latitude | degrees |
| longitude | numeric | GPS longitude | degrees |

**Note**: All models use natural logarithms (R's `log()` function). The `log_volume` column is ln(volume), not log₁₀(volume).

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

**Rows**: 243 species (OTUs)
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

### 7. Co-occurrence outputs (from `06_cooccurrence_analysis.R`)

The current pipeline generates co-occurrence analysis outputs including `pairwise_cooccurrence.csv` (pairwise SES), `intraspecific_density.csv`, `size_dependent_cooccurrence.csv`, `network_metrics.csv`, and `hub_species.csv`. See `output/tables/` for details. The legacy files `cafi_network_metrics.csv` and `cafi_keystone_species.csv` (from the deprecated `archive/06_network_analysis_legacy.R`) are no longer generated.

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
- **Physiology subset**: Only 108/114 corals have physiological data
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
Stier AC, et al. (2026). Sublinear scaling of coral-associated fauna reveals
propagule limitation in established reef communities.
Marine Ecology Progress Series [in preparation].
```

Data archived at: [Dryad DOI to be assigned]

---

## Contact

**Principal Investigator**: Adrian Stier
**Email**: astier@ucsb.edu
**Institution**: UC Santa Barbara, Department of Ecology, Evolution, and Marine Biology
**Lab Website**: https://stierlab.com

---

*Data Dictionary v1.2 - March 2026*
