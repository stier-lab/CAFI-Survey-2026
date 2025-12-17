# CAFI Survey Analysis

**Landscape characteristics structure coral-associated fauna communities and their effects on coral condition**

Mo'orea, French Polynesia | Summer 2019 | *Pocillopora* corals

---

## What This Project Is

This repository contains a complete R analysis pipeline for a coral reef ecology manuscript. We surveyed 114 *Pocillopora* coral colonies across 3 reef sites and catalogued ~4,000 individual coral-associated fauna (CAFI) spanning 87 species. The analysis tests how coral size and neighborhood context shape CAFI community assembly, and whether CAFI communities in turn affect coral physiological condition.

**Key question**: Do landscape characteristics (coral size, local density) drive CAFI community structure, and do CAFI provide measurable benefits to their host corals?

---

## Quick Start

### Prerequisites

R 4.0+ with these packages (auto-installed by setup script):
- tidyverse, vegan, lme4, lmerTest, glmmTMB
- sf, geosphere, ggplot2, patchwork, viridis

### Run the Full Analysis

```r
# From R console in project root:
source("scripts/run_pipeline.R")
```

This runs the complete pipeline (~5-10 min) and generates:
- Manuscript figures in `output/figures/manuscript/`
- Exploratory figures in `output/figures/`
- Statistical tables in `output/tables/`
- Processed data objects in `output/objects/`

### Run Individual Scripts

```r
# 1. Load setup (required first)
source("scripts/00_setup.R")

# 2. Load and clean data (required second)
source("scripts/01_load_clean_data.R")

# 3. Run specific analysis
source("scripts/04_scaling_relationships.R")  # Fig 2: Size-abundance scaling
source("scripts/06_network_analysis.R")       # Fig 4: Networks
source("scripts/08_cafi_condition_feedbacks.R") # Fig 6: Feedbacks
```

---

## Repository Structure

```
CAFI-Survey-2026/
├── data/                          # Raw data files
│   ├── survey_cafi_data_*.csv     # CAFI specimen records (3,989 rows)
│   ├── survey_coral_*.csv         # Coral colony attributes (114 rows)
│   ├── survey_master_phys_*.csv   # Physiology measurements (108 rows)
│   ├── traits/                    # CAFI trait database
│   ├── gis/                       # Mo'orea shapefiles for maps
│   └── README.md                  # DATA DOCUMENTATION (start here)
│
├── scripts/                       # R analysis scripts (see below)
│   ├── 00_setup.R                 # Setup, packages, paths, theme
│   ├── 01_load_clean_data.R       # Data loading and cleaning
│   ├── 02-09_*.R                  # Core analyses + manuscript figures
│   ├── 10-13_*.R                  # Supplementary analyses
│   ├── run_pipeline.R             # Pipeline orchestrator
│   └── archive/                   # Deprecated scripts (reference only)
│
├── output/                        # Generated outputs
│   ├── figures/manuscript/        # Publication figures (Fig 2-6)
│   ├── figures/                   # Exploratory figures by analysis
│   ├── tables/                    # CSV statistical results
│   └── objects/                   # RDS R data objects
│
└── docs/                          # Documentation
    ├── KEY_RESULTS_SUMMARY.md     # Main findings
    ├── ANALYSIS_METHODS_SUMMARY.md # Detailed methods
    └── PRD.md                     # Product requirements
```

---

## Analysis Pipeline

### Script Overview

| Script | Purpose | Manuscript Figure |
|--------|---------|-------------------|
| `00_setup.R` | Load packages, configure paths, define theme | — |
| `01_load_clean_data.R` | Load CSVs, clean data, create RDS objects | — |
| `02_community_composition.R` | Taxonomic summaries, PERMANOVA | — |
| `03_diversity_analysis.R` | Alpha/beta diversity, NMDS | Fig 4 (ordination) |
| `04_scaling_relationships.R` | Size-abundance power-law scaling | **Fig 2** |
| `05_coral_condition.R` | Position-corrected condition scores | **Fig 5** |
| `06_network_analysis.R` | Co-occurrence networks, modularity | **Fig 4** (networks) |
| `07_neighborhood_effects.R` | Meter-scale neighborhood effects | **Fig 6** (partial) |
| `08_cafi_condition_feedbacks.R` | CAFI -> coral condition | **Fig 6** |
| `09_functional_groups.R` | Trapezia, fish, corallivore patterns | **Fig 3** |
| `10_spatial_patterns.R` | Maps, spatial extent | Supplementary |
| `11_size_biomass_scaling.R` | CAFI body size distributions | Supplementary |
| `12_machine_learning.R` | Random Forest predictions | Supplementary |
| `13_spatial_autocorrelation.R` | Moran's I, spatial regression | Supplementary |

### Dependency Chain

```
00_setup.R
    ↓
01_load_clean_data.R
    ↓
┌───┴───────────────────────────────────┐
│  02_community_composition.R           │
│  03_diversity_analysis.R     → Fig 4  │
│  04_scaling_relationships.R  → Fig 2  │
│  05_coral_condition.R        → Fig 5  │
│  06_network_analysis.R       → Fig 4  │
│  07_neighborhood_effects.R   → Fig 6  │
│  08_cafi_condition_feedbacks.R → Fig 6│
│  09_functional_groups.R      → Fig 3  │
└───────────────────────────────────────┘
```

---

## Data Overview

| Dataset | File | Records | Key Variables |
|---------|------|---------|---------------|
| CAFI specimens | `survey_cafi_data_w_taxonomy_summer2019_v5.csv` | 3,989 | coral_id, type, genus, species, size_mm |
| Coral colonies | `survey_coral_characteristics_merged_v2.csv` | 114 | coral_id, lat/long, volume, branch_width, neighbor metrics |
| Physiology | `survey_master_phys_data_v3.csv` | 108 | coral_id, protein, carbs, zoox_density, afdw |

**Join key**: All datasets link via `coral_id` (format: "SITE-POC##", e.g., "HAU-POC29")

**Full data documentation**: See [data/README.md](data/README.md)

---

## Key Research Questions

### 1. How does coral size affect CAFI abundance? (Fig 2)

The relationship between coral size (V) and CAFI abundance (N) follows a power law: **N = a × V^β**

Three competing hypotheses:

| Hypothesis | Scaling | Mechanism |
|------------|---------|-----------|
| **Field of Dreams** | β = 1 | Abundance proportional to size; passive habitat filling |
| **Propagule Redistribution** | β < 1 | Larger corals "dilute" settlers; per-capita density decreases |
| **Super-linear** | β > 1 | Larger corals disproportionately attractive; aggregation |

**Key Results**:
- **Total abundance**: β = 1.20 (super-linear) — larger corals harbor more per unit volume
- **Species richness**: β = 0.81 (redistribution) — diversity saturates in larger corals
- **Individual species**: Most show β ≈ 1 (Field of Dreams) — proportional scaling

See `output/reports/SCALING_ANALYSIS.html` for comprehensive analysis.

### 2. Do CAFI functional groups affect coral condition? (Figs 3, 6)

| Functional Group | Taxa | Expected Effect |
|-----------------|------|-----------------|
| Mutualist defenders | *Trapezia*, *Tetralia* crabs | Positive |
| Nutrient providers | *Paragobiodon*, *Caracanthus* fish | Positive |
| Corallivores | *Drupella*, *Coralliophila* snails | Negative |

**Result**: Defender abundance correlates with higher condition; corallivores with lower.

### 3. Do local neighborhoods affect CAFI? (Figs 4, 6)

**Result**: Neighbor density has weak positive effect; coral size dominates.

---

## Study Sites

| Code | Name | Environment | N |
|------|------|-------------|---|
| HAU | Hauru | Fringing reef, north shore | 38 |
| MAT | Maatea | Lagoon, back reef | 38 |
| MRB | Barrier Reef | Outer barrier, oceanic | 38 |

**Location**: Mo'orea, French Polynesia (17°30'S, 149°50'W)

---

## Working with Pre-computed Objects

After running `01_load_clean_data.R`, these RDS files are available:

```r
# Load pre-computed objects
coral_master <- readRDS("output/objects/coral_master.rds")       # Main merged dataset
cafi_clean <- readRDS("output/objects/cafi_clean.rds")           # Clean CAFI records
community_matrix <- readRDS("output/objects/community_matrix.rds") # Coral x OTU matrix
condition_scores <- readRDS("output/objects/condition_scores.rds") # Condition PC1
functional_summary <- readRDS("output/objects/functional_summary.rds") # By functional group
```

---

## Key Outputs

### Manuscript Figures (`output/figures/manuscript/`)

| Figure | Content | Script |
|--------|---------|--------|
| Fig 2 | Size-abundance-richness scaling | `04_scaling_relationships.R` |
| Fig 3 | Functional group responses | `09_functional_groups.R` |
| Fig 4 | NMDS ordination + networks | `03_diversity_analysis.R`, `06_network_analysis.R` |
| Fig 5 | Condition vs landscape | `05_coral_condition.R` |
| Fig 6 | CAFI-condition feedbacks | `07_neighborhood_effects.R`, `08_cafi_condition_feedbacks.R` |

### Statistical Tables (`output/tables/`)

Key output files:
- `scaling_model_results.csv` - Size-abundance coefficients
- `permanova_results.csv` - Community composition tests
- `network_metrics.csv` - Modularity, centrality
- `condition_model_results.csv` - Condition predictors

---

## Important Notes

### Data Quirks

1. **Site extraction**: The `site` column in CAFI data is sometimes blank. Extract from `coral_id` prefix (HAU, MAT, MRB).

2. **Volume**: Use `volume_field` (field estimate) since `volume_lab` is often missing.

3. **Position correction**: Physiology metrics require position correction—sampling position correlates with colony size. Scripts handle this automatically.

4. **Taxonomy**: CAFI are morphological OTUs, not genetically confirmed species.

### Finding Manuscript Figure Code

Each script that generates manuscript figures has a clear marker:

```r
# >>> MANUSCRIPT FIGURE X <<<
```

Search for this pattern to quickly find publication figure code.

---

## Citation

```
Stier AC, et al. (2026). Landscape characteristics structure coral-associated
fauna communities and their effects on coral condition in Mo'orea, French Polynesia.
Marine Ecology Progress Series [submitted].
```

---

## Contact

**PI**: Adrian Stier
**Email**: astier@ucsb.edu
**Lab**: Stier Lab, UC Santa Barbara
**Website**: https://stierlab.com

---

## License

MIT License

---

*Last updated: December 2025*
