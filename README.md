# CAFI Survey Analysis

**Sublinear scaling of coral-associated fauna reveals propagule limitation in established reef communities**

Mo'orea, French Polynesia | Summer 2019 | *Pocillopora* corals

---

## What This Project Is

This repository contains a complete R analysis pipeline for a coral reef ecology manuscript. We surveyed 114 *Pocillopora* coral colonies across 3 reef sites and catalogued ~4,000 individual coral-associated fauna (CAFI) spanning 243 OTUs. The analysis tests how coral size structures CAFI abundance, richness, and community composition in established reef communities.

**Key finding**: CAFI abundance scales sublinearly with coral size (β = 0.52), supporting propagule dilution ("Redirection") over proportional scaling ("Field of Dreams"). Community composition is structured by site and size. Species richness predicts coral condition (Hochberg-corrected p = 0.036), supporting a biodiversity-ecosystem function (BEF) relationship.

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

This runs the complete pipeline (~12 min) and generates:
- 5 manuscript figures in `output/figures/manuscript/`
- Supplementary figures in `output/figures/supplement/`
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

# Run everything:
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
source("scripts/06_cooccurrence_analysis.R")       # Co-occurrence null models
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
│   ├── traits/                    # CAFI trait database (cafi_traits_final.csv)
│   └── README.md                  # DATA DOCUMENTATION (start here)
│
├── scripts/                       # R analysis scripts (see below)
│   ├── run_full_pipeline.R        # MASTER PIPELINE RUNNER
│   ├── 00_setup.R                 # Setup, packages, paths, theme
│   ├── 00b_data_quality_audit.R   # Data quality assessment
│   ├── 01_load_data.R             # Data loading, cleaning, Fig 1
│   ├── 02_community_analysis.R    # Community composition (Q2)
│   ├── 03_landscape_characterization.R # Landscape metrics
│   ├── 04_landscape_effects.R     # Neighborhood effects (supplement)
│   ├── 05_species_scaling_analysis.R # Size-abundance scaling (Q1), Fig 2
│   ├── 06_cooccurrence_analysis.R # Co-occurrence null models (supplement)
│   ├── 07_spatial_autocorrelation.R # Spatial patterns (Moran's I)
│   ├── 08_functional_groups.R     # Taxonomic group scaling, Fig S7
│   ├── 09_cafi_condition_feedbacks.R # CAFI-condition feedbacks (Q3), Fig 5
│   ├── 13_taxonomy_sensitivity.R  # Taxonomy sensitivity analysis, Fig S6
│   └── archive/                   # Archived scripts (not part of manuscript or pipeline)
│
├── output/                        # Generated outputs (gitignored)
│   ├── figures/
│   │   ├── manuscript/            # 5 publication figures (fig1-fig5) + legend files
│   │   ├── supplement/            # Supplementary figures (figS1-S14)
│   │   ├── 01_data/              # Study design figures
│   │   ├── 02_community/         # Community analysis
│   │   ├── 03_landscape/         # Landscape characterization
│   │   ├── 04_effects/           # Landscape effects
│   │   ├── 05_scaling/           # Species-area scaling
│   │   ├── 06_network/           # Co-occurrence analysis
│   │   ├── feedbacks/            # CAFI-condition feedbacks
│   │   └── functional_groups/    # Taxonomic group analysis
│   ├── tables/                    # CSV statistical results
│   └── objects/                   # RDS R data objects
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
| `02_community_analysis.R` | PERMANOVA, NMDS, betadisper, rarefaction | **Fig 4**, Figs S1-S3 |
| `03_landscape_characterization.R` | Neighborhood density, predictor selection | Landscape metrics |
| `04_landscape_effects.R` | Size and neighborhood effects on CAFI | Fig S5, Q4 tables |
| `05_species_scaling_analysis.R` | NB GLM scaling, bootstrap CI | **Fig 2**, **Fig 3**, Fig S4 |

#### Extended Analyses (scripts 06-09, 13)
| Script | Purpose | Output |
|--------|---------|--------|
| `06_cooccurrence_analysis.R` | Pairwise co-occurrence null models, intraspecific density | Fig S9 |
| `07_spatial_autocorrelation.R` | Moran's I, LISA, Mantel tests | Table S10 |
| `08_functional_groups.R` | Taxonomic group scaling and composition | Fig S7 |
| `09_cafi_condition_feedbacks.R` | PCA, fixed-effect LMs, Hochberg FWER correction | **Fig 5** |
| `13_taxonomy_sensitivity.R` | Taxonomy robustness across 5 OTU scenarios | Fig S6 |

#### Archived Scripts (`scripts/archive/` — NOT part of manuscript or pipeline)
| Script | Purpose |
|--------|---------|
| `archive/06_network_analysis_legacy.R` | Legacy network analysis (replaced by `06_cooccurrence_analysis.R`) |
| `archive/10_feature_engineering.R` | Exploratory ML feature creation |
| `archive/11_machine_learning.R` | Exploratory Random Forest, XGBoost |
| `archive/12_model_evaluation.R` | Exploratory cross-validation |

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
│  ├─ 06_cooccurrence_analysis.R   (Q2)      │
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

## Research Questions (Q1-Q3)

This study addresses three core questions about how coral size structures CAFI communities:

### Q1: SCALING — How does CAFI abundance scale with coral size?

**N = a × V^β** — Three competing hypotheses:

| Hypothesis | Scaling | Mechanism |
|------------|---------|-----------|
| **Field of Dreams** | β = 1 | Abundance proportional to size; passive habitat filling |
| **Propagule Redirection** | β < 1 | Larger corals "dilute" settlers; per-capita density decreases |
| **Super-linear** | β > 1 | Larger corals disproportionately attractive |

**Result**: Total CAFI abundance β = 0.52 [0.44, 0.62] — **sublinear (Propagule Redirection)**. Species richness z = 0.34 [0.27, 0.42] — sublinear. Per-capita density decreases with size (dilution slope = -0.48). Taxonomic groups vary: Gastropods β ≈ 0.94 (Field of Dreams), Trapezia β ≈ 0.43 (Redirection), Fish β ≈ 0.74 (Redirection).

### Q2: COMPOSITION — What structures CAFI community composition?

Tests whether community *composition* (not just richness) diverges with coral size, and what structures species co-occurrence patterns.

**Methods**: PERMANOVA, betadisper (distance to centroid), rarefaction robustness check, NMDS ordination, db-RDA constrained ordination, NODF nestedness test, pairwise co-occurrence null models (volume-weighted Bernoulli, 10,000 iterations).

**Result**: Sites strongly structure composition (PERMANOVA: site R² = 0.06, volume R² = 0.08, both p < 0.001). Coral volume explains 7.8% of composition (db-RDA, p = 0.001), robust to rarefaction (2.6%, p = 0.001). Categorical size-divergence is **NOT significant after rarefaction** (p = 0.61) — an abundance artifact. Communities are not nested along the size gradient (NODF p = 0.28), indicating species turnover. No pairwise co-occurrence is significant after FDR correction (0/528 pairs); intraspecific density shows mating-pair patterns in 6/15 species.

### Q3: FEEDBACKS — Does CAFI community identity predict coral condition?

Tests whether CAFI identity provides measurable physiological benefits to host corals. FDR-corrected for multiple testing.

| Predictor | p-value | p_Hochberg | Direction |
|-----------|---------|------------|-----------|
| Species richness | 0.018 | 0.036 | Positive |
| Total abundance | 0.048 | 0.048 | Marginal |

**Result**: Species richness predicts coral condition (Hochberg-corrected p = 0.036), supporting a BEF (biodiversity-ecosystem function) relationship. Variance partitioning attributes 29% uniquely to richness vs <1% to abundance. Rarefied richness (n=20) is non-significant (p = 0.50), but this test is ambiguous — rarefaction may remove the BEF mechanism itself (diversity→abundance→condition). Community identity (PC1_CAFI) does not predict condition (supplement).

### Supporting: NEIGHBORHOOD (Supplement)

Neighborhood density (n_neighbors within 5m) was tested as a predictor of CAFI abundance and coral condition but showed no evidence of effects (all p > 0.18). This analysis is underpowered (n=61 corals with neighborhood data) for detecting small effects and is presented in the supplement (Fig S5).

---

## Study Sites

| Code | Name | Environment | N |
|------|------|-------------|---|
| HAU | Hauru | Fringing reef, north shore | 38 |
| MAT | Maatea | Lagoon, back reef | 39 |
| MRB | Maharepa | Barrier reef, north shore | 35 |

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
| Fig 1 | Study design + representative CAFI: satellite map + distributions + species photos (6-panel) | Overview | `01_load_data.R` |
| Fig 2 | Scaling: abundance + richness + density dilution (3 vertical panels) | Q1 | `05_species_scaling_analysis.R` |
| Fig 3 | Species + taxonomic group scaling (2×2: curves + β forest plots) | Q1 | `05_species_scaling_analysis.R` |
| Fig 4 | Composition: NMDS ordination + taxonomic barchart by site | Q2 | `02_community_analysis.R` |
| Fig 5 | BEF diversity-condition: richness scatter + abundance scatter (2-panel) | Q3 | `09_cafi_condition_feedbacks.R` |

### Supplementary Figures (`output/figures/supplement/`)

| Figure | Content | Script |
|--------|---------|--------|
| S1 | Species accumulation curves | `02_community_analysis.R` |
| S2 | PERMANOVA metric sensitivity | `02_community_analysis.R` |
| S3 | Composition divergence by size | `02_community_analysis.R` |
| S4 | Species-level scaling forest plot | `05_species_scaling_analysis.R` |
| S5 | Neighborhood null results | `04_landscape_effects.R` |
| S6 | Taxonomy sensitivity forest plot | `13_taxonomy_sensitivity.R` |
| S7 | Taxonomic group scaling and composition | `08_functional_groups.R` |
| S8 | Rarefaction depth sensitivity | `09_cafi_condition_feedbacks.R` |
| S9 | Co-occurrence: pairwise SES heatmap + intraspecific density + size-dependent (3-panel) | `06_cooccurrence_analysis.R` |
| S10 | BEF variance partitioning + partial regression + path model | `09_cafi_condition_feedbacks.R` |
| S11 | A priori forest + rarefied richness + exploratory forest + species scatters + bidirectional | `09_cafi_condition_feedbacks.R` |
| S12 | Species occurrence probability vs. coral size | `05_species_scaling_analysis.R` |
| S13 | Species × trait heatmap: 19 prevalent species (≥5 corals) × 5 condition metrics | `09_cafi_condition_feedbacks.R` |
| S14 | Species × trait biplots: strongest associations | `09_cafi_condition_feedbacks.R` |

### Statistical Tables (`output/tables/`)

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

### Data Objects (`output/objects/`)

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
| Multiple testing (feedbacks) | Three-tier: Hochberg FWER (a priori BEF, k=2), BH-FDR (exploratory, k=4) | `09` |
| Multiple testing (key species) | Hochberg FWER across up to 10 species tests | `09` |
| Multiple testing (scaling) | FDR within category (species/group) | `05` |
| Multiple testing (co-occurrence) | FDR across pairwise SES tests | `06` |
| Abundance confound (composition) | Iterated rarefaction (100 draws) | `02` |
| Abundance confound (richness-condition) | Rarefied richness artifact test | `09` |
| Volume confound (co-occurrence) | Volume-weighted Bernoulli null model | `06` |
| Random effects (k=3 sites) | Fixed-effect site (Bolker et al. 2009) | `04`, `09` |
| Bootstrap site structure | Stratified bootstrap (`strata` argument) | `05` |
| NB convergence failure | Poisson fallback with logging | `04` |
| Colorblind accessibility | Colorblind-safe palettes (Okabe-Ito + custom site/scaling) | All scripts |

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
Stier AC, et al. (2026). Sublinear scaling of coral-associated fauna reveals
propagule limitation in established reef communities. [in preparation].
```

---

## Contact

**PI**: Adrian Stier
**Email**: astier@ucsb.edu
**Lab**: Stier Lab, UC Santa Barbara
**Website**: https://stierlab.com

---

## License

CC-BY-4.0 — see [LICENSE](LICENSE) for details.

---

*Last updated: March 2026*
