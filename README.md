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
source("scripts/run_full_pipeline.R")
run_pipeline()
```

This runs the complete pipeline (~5-10 min) and generates:
- Manuscript figures in `output/figures/manuscript/`
- Exploratory figures in `output/figures/`
- Statistical tables in `output/tables/`
- Processed data objects in `output/objects/`
- Execution log in `output/pipeline.log`

**Pipeline options:**

```r
# Skip already-completed steps:
run_pipeline(skip_completed = TRUE)

# Run only core pipeline:
run_core_pipeline()

# Check pipeline status:
check_pipeline_status()
```

**From command line:**

```bash
Rscript scripts/run_full_pipeline.R
Rscript scripts/run_full_pipeline.R --skip-completed
```

### Run Individual Scripts

```r
# 1. Load setup (required first)
source("scripts/00_setup.R")

# 2. Load and clean data (required second)
source("scripts/01_load_data.R")

# 3. Run specific analysis
source("scripts/02_community_analysis.R")        # Community composition
source("scripts/03_landscape_characterization.R") # Landscape metrics
source("scripts/04_landscape_effects.R")          # Landscape effects on CAFI
source("scripts/05_species_scaling_analysis.R")   # Species-area scaling
source("scripts/06_network_analysis.R")           # Co-occurrence networks
source("scripts/07_spatial_autocorrelation.R")    # Spatial patterns
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
│   ├── run_full_pipeline.R        # MASTER PIPELINE RUNNER
│   ├── 00_setup.R                 # Setup, packages, paths, theme
│   ├── 00b_data_quality_audit.R   # Data quality assessment
│   ├── 01_load_data.R             # Data loading and cleaning
│   ├── 02_community_analysis.R    # Community composition (Q2)
│   ├── 03_landscape_characterization.R # Landscape metrics
│   ├── 04_landscape_effects.R     # Neighborhood effects (Q4)
│   ├── 05_species_scaling_analysis.R # Size-abundance scaling (Q1)
│   ├── 06_network_analysis.R      # Co-occurrence networks
│   ├── 07_spatial_autocorrelation.R # Spatial patterns (Moran's I)
│   ├── 08_functional_groups.R     # Functional group scaling
│   ├── 09_cafi_condition_feedbacks.R # CAFI-condition feedbacks (Q3)
│   ├── 10_feature_engineering.R   # Feature creation for ML
│   ├── 11_machine_learning.R      # RF + XGBoost models
│   ├── 12_model_evaluation.R      # Cross-validation, diagnostics
│   ├── 13_manuscript_figures.R    # Publication figures (Q1-Q4)
│   └── archive/                   # Deprecated scripts (reference only)
│
├── output/                        # Generated outputs
│   ├── figures/manuscript/        # Publication figures (Fig 2-6)
│   ├── figures/                   # Exploratory figures by analysis
│   ├── tables/                    # CSV statistical results
│   └── objects/                   # RDS R data objects
│
└── docs/                          # Documentation
    └── METHODS.md                 # Detailed statistical methods
```

---

## Analysis Pipeline

### Script Overview

#### Core Pipeline
| Script | Purpose | Description |
|--------|---------|-------------|
| `00_setup.R` | Setup | Load packages, configure paths, define ggplot theme |
| `00b_data_quality_audit.R` | Quality | Data quality assessment, trait integration |
| `01_load_data.R` | Data loading | Load CSVs, clean data, create RDS objects |
| `02_community_analysis.R` | Community | Taxonomic summaries, PERMANOVA, composition |
| `03_landscape_characterization.R` | Landscape | Neighborhood density, predictor selection |
| `04_landscape_effects.R` | Effects | Size and neighborhood effects on CAFI |
| `05_species_scaling_analysis.R` | Scaling | Species-area relationships, power-law |

#### Extended Analyses
| Script | Purpose | Question |
|--------|---------|----------|
| `06_network_analysis.R` | Co-occurrence networks, modularity | Q2 |
| `07_spatial_autocorrelation.R` | Moran's I, LISA, Mantel tests | Q2 |
| `08_functional_groups.R` | Trapezia, fish, corallivore patterns | Q1/Q3 |
| `09_cafi_condition_feedbacks.R` | PCA, mixed models, FDR correction | Q3/Q4 |

#### Machine Learning & Publication
| Script | Purpose |
|--------|---------|
| `10_feature_engineering.R` | Feature creation, VIF selection |
| `11_machine_learning.R` | Random Forest, XGBoost models |
| `12_model_evaluation.R` | Cross-validation, diagnostics |
| `13_manuscript_figures.R` | Publication figures organized by Q1-Q4 |

**Note**: Archived scripts are in `scripts/archive/` for reference only.

### Dependency Chain

```
run_full_pipeline.R (orchestrates all scripts)
        ↓
00_setup.R → 00b_data_quality_audit.R → 01_load_data.R
        ↓
┌─────────────────────────────────────────────┐
│  Core Analyses (Q1-Q2):                     │
│  ├─ 02_community_analysis.R      (Q2)      │
│  ├─ 03_landscape_characterization.R        │
│  ├─ 04_landscape_effects.R       (Q4)      │
│  └─ 05_species_scaling_analysis.R (Q1)     │
└─────────────────────────────────────────────┘
        ↓
┌─────────────────────────────────────────────┐
│  Extended Analyses (Q2-Q4):                 │
│  ├─ 06_network_analysis.R        (Q2)      │
│  ├─ 07_spatial_autocorrelation.R  (Q2)     │
│  ├─ 08_functional_groups.R       (Q1/Q3)   │
│  └─ 09_cafi_condition_feedbacks.R (Q3/Q4)  │
└─────────────────────────────────────────────┘
        ↓
┌─────────────────────────────────────────────┐
│  Machine Learning & Figures:                │
│  ├─ 10_feature_engineering.R                │
│  ├─ 11_machine_learning.R                   │
│  ├─ 12_model_evaluation.R                   │
│  └─ 13_manuscript_figures.R                 │
└─────────────────────────────────────────────┘
        ↓
output/pipeline.log
```

**Note**: Scripts 02-09 can be run independently after setup and data loading.

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

## Research Questions (Q1-Q4)

This study addresses four core questions linking coral habitat, CAFI communities, and coral condition:

### Q1: SCALING — How does CAFI abundance scale with coral size?

**N = a × V^β** — Three competing hypotheses:

| Hypothesis | Scaling | Mechanism |
|------------|---------|-----------|
| **Field of Dreams** | β = 1 | Abundance proportional to size; passive habitat filling |
| **Propagule Redirection** | β < 1 | Larger corals "dilute" settlers; per-capita density decreases |
| **Super-linear** | β > 1 | Larger corals disproportionately attractive |

**Result**: Total CAFI abundance β = 1.20 [1.01, 1.40] — marginally super-linear. Species richness z = 0.79 [0.69, 0.89] — sublinear. Per-capita density decreases with size (dilution slope = -0.42). Functional groups vary: Trapezia β ≈ 1.0 (linear), Fish β ≈ 1.7 (super-linear).

### Q2: COMPOSITION — Do larger corals support more distinct communities?

Tests whether community *composition* (not just richness) diverges with coral size — i.e., whether larger corals accumulate unique species assemblages rather than just more individuals.

**Methods**: PERMANOVA, betadisper (distance to centroid), rarefaction robustness check (≥5 individuals, n=105).

**Result**: Sites strongly structure composition. Raw analysis shows community distinctness decreases with size, but this is **NOT significant after rarefaction** (p=0.61). The size-composition pattern is an abundance sampling artifact, not true divergence.

### Q3: FEEDBACKS — Does CAFI community identity predict coral condition?

Tests whether the identity of CAFI (mutualists vs parasites) provides measurable physiological benefits to host corals. FDR-corrected for multiple testing.

| Predictor | p-value | p_FDR | Direction |
|-----------|---------|-------|-----------|
| Species richness | 0.008 | 0.053 (marginal) | Positive |
| Total abundance | 0.10 | 0.25 | Positive |
| PC1_CAFI | 0.55 | 0.55 | — |

**Result**: Species richness → condition is the strongest signal (p=0.008) but marginal after FDR correction (p_FDR=0.053). Suggestive that community complexity may matter for coral health, but requires larger samples to confirm.

### Q4: NEIGHBORHOOD — Does local coral density affect CAFI recruitment?

Tests whether neighborhood density (n_neighbors within 5m) acts as a source of CAFI recruitment or facilitates coral condition — the observational analog to the experimental density manipulation.

**Result**: n_neighbors NOT significant for CAFI abundance (p=0.86) or condition (p=0.93). Coral volume remains the dominant predictor. Neighborhood density does not explain CAFI or condition variation in this established community.

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

After running `01_load_data.R`, these RDS files are available:

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

| Figure | Content | Question | Script |
|--------|---------|----------|--------|
| Fig 1 | Study design, sites, dataset summary | Overview | `13_manuscript_figures.R` |
| Fig 2 | Size-abundance scaling (power law) | Q1 | `13_manuscript_figures.R`, `05` |
| Fig 3 | Functional group scaling patterns | Q1 | `08_functional_groups.R` |
| Fig 4 | Co-occurrence network analysis | Q2 | `06_network_analysis.R` |
| Fig 5 | CAFI-condition feedbacks | Q3 | `09_cafi_condition_feedbacks.R` |

### Statistical Tables (`output/tables/`)

Key output files:
- `scaling_results_all.csv` - Species-area scaling coefficients
- `landscape_effects_summary.csv` - Predictor effect sizes
- `network_metrics.csv` - Network structure (modularity, centrality)
- `morans_i_results.csv` - Spatial autocorrelation tests
- `pipeline_timing.csv` - Script execution times

### Pipeline Outputs

- `output/pipeline.log` - Detailed execution log with timing
- `output/tables/pipeline_timing.csv` - Script-by-script performance

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

*Last updated: January 2026*
