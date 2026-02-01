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
- tidyverse, vegan, MASS, car, patchwork, viridis, scales
- sf, geosphere, igraph, ggrepel
- sandwich, lmtest, DHARMa, MuMIn

### Run the Full Analysis

```r
# From R console in project root:
source("scripts/run_full_pipeline.R")
run_pipeline()
```

This runs the complete pipeline (~6-7 min) and generates:
- 5 manuscript figures in `output/figures/manuscript/`
- 7 supplementary figures in `output/figures/supplement/`
- Exploratory figures in `output/figures/` (by analysis)
- Statistical tables in `output/tables/`
- Processed data objects in `output/objects/`
- Execution log in `output/pipeline.log`

**Pipeline options:**

```r
# Skip already-completed steps:
run_pipeline(skip_completed = TRUE)

# Run only core pipeline (scripts 00-05):
run_core_pipeline()

# Run everything including ML exploration:
run_full_pipeline()

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
source("scripts/08_functional_groups.R")          # Taxonomic group scaling
source("scripts/09_cafi_condition_feedbacks.R")   # CAFI-condition feedbacks
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
│   ├── 01_load_data.R             # Data loading, cleaning, Fig 1
│   ├── 02_community_analysis.R    # Community composition (Q2)
│   ├── 03_landscape_characterization.R # Landscape metrics
│   ├── 04_landscape_effects.R     # Neighborhood effects (Q4)
│   ├── 05_species_scaling_analysis.R # Size-abundance scaling (Q1), Fig 2
│   ├── 06_network_analysis.R      # Co-occurrence networks, Fig 4
│   ├── 07_spatial_autocorrelation.R # Spatial patterns (Moran's I)
│   ├── 08_functional_groups.R     # Taxonomic group scaling, Fig 3
│   ├── 09_cafi_condition_feedbacks.R # CAFI-condition feedbacks (Q3), Fig 5
│   ├── 10-12_*.R                  # Machine learning (exploratory, not in default pipeline)
│   └── archive/                   # Deprecated scripts (reference only)
│
├── output/                        # Generated outputs (gitignored)
│   ├── figures/
│   │   ├── manuscript/            # 5 publication figures (fig1-fig5)
│   │   ├── supplement/            # 7 supplementary figures (figS1-S7)
│   │   ├── 01_data/              # Study design figures
│   │   ├── 02_community/         # Community analysis (11 figures)
│   │   ├── 03_landscape/         # Landscape characterization (3 figures)
│   │   ├── 04_effects/           # Landscape effects (14 figures)
│   │   ├── 05_scaling/           # Species-area scaling (7 figures)
│   │   ├── 06_network/           # Network analysis (1 figure)
│   │   ├── feedbacks/            # CAFI-condition feedbacks (8 figures)
│   │   └── functional_groups/    # Taxonomic group analysis (7 figures)
│   ├── tables/                    # ~46 CSV statistical results
│   └── objects/                   # 17 RDS R data objects
│
├── manuscript/                    # Manuscript drafts
├── CLAUDE.md                      # AI assistant context (detailed)
└── README.md                      # This file
```

---

## Analysis Pipeline

### Script-to-Question Mapping

Each manuscript figure is created by its source analysis script with dual saves (analysis directory + `output/figures/manuscript/`).

#### Core Pipeline (scripts 00-05)
| Script | Purpose | Output |
|--------|---------|--------|
| `00_setup.R` | Load packages, paths, ggplot theme, helpers | Global environment |
| `00b_data_quality_audit.R` | Data quality assessment, trait integration | QC tables |
| `01_load_data.R` | Load CSVs, clean data, create RDS objects | **Fig 1** (study design) |
| `02_community_analysis.R` | PERMANOVA, NMDS, betadisper, rarefaction | Figs S1-S3, S5 |
| `03_landscape_characterization.R` | Neighborhood density, predictor selection | Landscape metrics |
| `04_landscape_effects.R` | Size and neighborhood effects on CAFI | Fig S7, Q4 tables |
| `05_species_scaling_analysis.R` | NB GLM scaling, bootstrap CI | **Fig 2**, Fig S6 |

#### Extended Analyses (scripts 06-09)
| Script | Purpose | Output |
|--------|---------|--------|
| `06_network_analysis.R` | Co-occurrence networks, modularity, hubs | **Fig 4** |
| `07_spatial_autocorrelation.R` | Moran's I, LISA, Mantel tests | Fig S4 |
| `08_functional_groups.R` | Taxonomic group scaling and composition | **Fig 3** |
| `09_cafi_condition_feedbacks.R` | PCA, fixed-effect LMs, FDR correction | **Fig 5** |

#### Exploratory ML (not in default pipeline)
| Script | Purpose |
|--------|---------|
| `10_feature_engineering.R` | Feature creation, VIF selection |
| `11_machine_learning.R` | Random Forest, XGBoost models |
| `12_model_evaluation.R` | Cross-validation, diagnostics |

### Dependency Chain

```
run_full_pipeline.R (orchestrates all scripts)
        ↓
00_setup.R → 00b_data_quality_audit.R → 01_load_data.R
        ↓
┌─────────────────────────────────────────────┐
│  Core Analyses (Q1-Q2):                     │
│  ├─ 02_community_analysis.R      (Q2)      │
│  ├─ 03_landscape_characterization.R         │
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

**Result**: Total CAFI abundance β = 0.52 [0.44, 0.61] — **sublinear (Propagule Redirection)**. Species richness z = 0.35 [0.28, 0.43] — sublinear. Per-capita density decreases with size (dilution slope = -0.48). Taxonomic groups vary: Gastropods β ≈ 0.94 (Field of Dreams), Trapezia β ≈ 0.43 (Redirection), Fish β ≈ 0.74 (Redirection).

### Q2: COMPOSITION — What structures CAFI community composition?

Tests whether community *composition* (not just richness) diverges with coral size, and what structures species co-occurrence patterns.

**Methods**: PERMANOVA, betadisper (distance to centroid), rarefaction robustness check, NMDS ordination, co-occurrence network analysis (Louvain modularity, Erdos-Renyi null model).

**Result**: Sites strongly structure composition (PERMANOVA R² ~ 0.04-0.06, p < 0.01). Size-divergence is **NOT significant after rarefaction** (p = 0.61) — the size-composition pattern is an abundance artifact. Co-occurrence network reveals non-random modular assembly with identifiable hub species (modularity significantly exceeds null expectation).

### Q3: FEEDBACKS — Does CAFI community identity predict coral condition?

Tests whether CAFI identity provides measurable physiological benefits to host corals. FDR-corrected for multiple testing.

| Predictor | p-value | p_FDR | Direction |
|-----------|---------|-------|-----------|
| Species richness | 0.018 | 0.126 | Positive |
| Total abundance | 0.137 | 0.438 | Positive |
| PC1_CAFI | 0.837 | 0.837 | — |

**Result**: Species richness shows the strongest raw signal (p = 0.018) but does NOT survive FDR correction (p_FDR = 0.126). Critically, rarefied richness (controlling for sampling intensity) shows NO relationship with condition (p = 0.45), revealing the raw richness signal as an **abundance artifact**. No CAFI metric reliably predicts coral condition in this observational dataset.

### Q4: NEIGHBORHOOD — Does local coral density affect CAFI recruitment?

Tests whether neighborhood density (n_neighbors within 5m) acts as a source of CAFI recruitment or facilitates coral condition — the observational analog to the experimental density manipulation.

**Result**: n_neighbors NOT significant for CAFI abundance (p = 0.37) or condition (p = 0.78). Coral volume remains the dominant predictor. Neighborhood density does not explain CAFI or condition variation in this established community.

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
# Load using helper function (defined in 00_setup.R)
coral_master <- load_object("coral_master")       # Main merged dataset
cafi_clean <- load_object("cafi_clean")            # Clean CAFI records
community_matrix <- load_object("community_matrix") # Coral × OTU matrix
condition_scores <- load_object("condition_scores") # Condition PC1

# Or load directly
coral_master <- readRDS("output/objects/coral_master.rds")
```

---

## Key Outputs

### Manuscript Figures (`output/figures/manuscript/`)

| Figure | Content | Question | Script |
|--------|---------|----------|--------|
| Fig 1 | Study design: satellite map + volume/neighborhood histograms | Overview | `01_load_data.R` |
| Fig 2 | Size-abundance scaling (power law) + species forest plot | Q1 | `05_species_scaling_analysis.R` |
| Fig 3 | Taxonomic group scaling and composition | Q1 | `08_functional_groups.R` |
| Fig 4 | Co-occurrence network: hero layout + guild panels | Q2 | `06_network_analysis.R` |
| Fig 5 | CAFI-condition feedback models | Q3 | `09_cafi_condition_feedbacks.R` |

### Supplementary Figures (`output/figures/supplement/`)

| Figure | Content | Script |
|--------|---------|--------|
| S1 | Species accumulation curves | `02_community_analysis.R` |
| S2 | PERMANOVA metric sensitivity | `02_community_analysis.R` |
| S3 | NMDS ordination by site/size | `02_community_analysis.R` |
| S4 | Spatial autocorrelation (Moran's I) | `07_spatial_autocorrelation.R` |
| S5 | Composition divergence by size | `02_community_analysis.R` |
| S6 | Species-level scaling forest plot | `05_species_scaling_analysis.R` |
| S7 | Neighborhood null results | `04_landscape_effects.R` |

### Statistical Tables (`output/tables/`, 46 files)

Key output files:
- `scaling_results_all.csv` — Species-area scaling coefficients for all taxa
- `scaling_summary_by_category.csv` — Scaling summary grouped by taxonomic category
- `cafi_condition_models.csv` — CAFI → condition model results
- `key_species_effects.csv` — Experimental species predictions vs survey
- `network_metrics.csv` — Network structure (modularity, centrality)
- `hub_species.csv` — Hub species centrality scores
- `morans_i_results.csv` — Spatial autocorrelation tests
- `landscape_full_model_results.csv` — Full GLM results for landscape effects
- `pipeline_timing.csv` — Script execution times

### Data Objects (`output/objects/`, 17 files)

Key RDS files:
- `coral_master.rds` — Main merged dataset (coral + CAFI + physiology)
- `cafi_clean.rds` — Cleaned CAFI specimen records
- `community_matrix.rds` — Coral × species abundance matrix
- `condition_scores.rds` — PCA condition scores
- `cafi_network.rds` — Co-occurrence network + modularity results
- `scaling_analysis_results.rds` — Scaling model outputs + bootstrap CIs

---

## Analytical Quality Controls

| Issue | Fix | Script |
|-------|-----|--------|
| Multiple testing (feedbacks) | FDR correction (Benjamini-Hochberg) | `09` |
| Multiple testing (key species) | FDR across 6 species tests | `09` |
| Multiple testing (scaling) | FDR within category (species/group) | `05` |
| Multiple testing (network edges) | Pairwise cor.test p-values + FDR | `06` |
| Abundance confound (composition) | Iterated rarefaction (100 draws) | `02` |
| Abundance confound (richness-condition) | Rarefied richness artifact test | `09` |
| Volume confound (co-occurrence) | Residualized presence on log(volume) | `06` |
| Random effects (k=3 sites) | Fixed-effect site (Bolker et al. 2009) | `04`, `09` |
| Bootstrap site structure | Stratified bootstrap (`strata` argument) | `05` |
| NB convergence failure | Poisson fallback with logging | `04` |
| igraph/dplyr namespace conflict | Explicit `dplyr::filter()` throughout | `06` |
| Colorblind accessibility | Okabe-Ito palette throughout | All scripts |

---

## Important Notes

### Data Quirks

1. **Site extraction**: The `site` column in CAFI data is sometimes blank. Extract from `coral_id` prefix (HAU, MAT, MRB).

2. **Volume**: Use `volume` (field estimate). The `volume_field` alias was removed for consistency.

3. **Position correction**: Physiology metrics require position correction — sampling position correlates with colony size. Scripts handle this automatically.

4. **Taxonomy**: CAFI are morphological OTUs, not genetically confirmed species.

5. **Log base**: All models use natural log (`log()` in R). Previous versions used `log10()`, inflating coefficients by ln(10) ≈ 2.303.

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

*Last updated: February 2026*
