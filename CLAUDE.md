# CLAUDE.md - AI Assistant Context

This file provides essential context for AI assistants working with this codebase.

## Project Overview

**CAFI Survey Analysis** - A marine ecology research project studying how coral size and neighborhood context shape coral-associated fauna (CAFI) communities in Mo'orea, French Polynesia.

- **Survey scope**: 114 *Pocillopora* coral colonies across 3 reef sites
- **CAFI catalogued**: ~4,000 individual specimens spanning 87 species
- **Research focus**: Landscape effects on community assembly; CAFI feedbacks on coral condition

## Relationship to Experimental Study

This **observational survey** complements a parallel **experimental study** that manipulated coral density (1, 3, 6 corals per 25m² reef). The two studies test the same hypotheses using different approaches:

| Approach | Lever | Strength |
|----------|-------|----------|
| **Experiment** | Coral density (manipulated) | Causal inference, colonization dynamics |
| **Survey** | Coral size + neighborhood (natural variation) | Ecological realism, established communities |

## Four Core Research Questions (Q1-Q4)

| Question | Script | Key Result |
|----------|--------|------------|
| **Q1: SCALING** | `05_species_scaling_analysis.R` | Abundance β=1.20 (marginal super-linear); Richness z=0.79 (sublinear); density dilution |
| **Q2: COMPOSITION** | `02_community_analysis.R` | No size-divergence after rarefaction (abundance artifact) |
| **Q3: FEEDBACKS** | `09_cafi_condition_feedbacks.R` | Species richness → condition (p=0.008, p_FDR=0.053 marginal) |
| **Q4: NEIGHBORHOOD** | `04_landscape_effects.R` | n_neighbors NOT significant for CAFI or condition |

---

### Q1: Scaling — How does CAFI scale with coral size?

**Hypotheses tested:**
- **Field of Dreams (β = 1)**: Abundance scales proportionally—"if you build it, they will come"
- **Propagule Redistribution (β < 1)**: Abundance scales sublinearly—larger corals dilute per-capita colonization
- **Super-linear (β > 1)**: Larger corals attract disproportionately more settlers

**Methods:**
- Negative binomial GLM: `total_cafi ~ log10(volume) + site`
- Bootstrap 95% CI on scaling exponent (stratified by site)
- Species-level scaling for key CAFI taxa

**Key findings:**
- Total CAFI abundance: β = 1.20 [1.01, 1.40] — marginally super-linear
- Species richness (SAR): z = 0.79 [0.69, 0.89] — **sublinear, consistent with Redistribution**
- Density dilution: per-capita CAFI density decreases with size (slope = -0.42)
- Trapezia: β ≈ 0.99 (Field of Dreams); Fish: β ≈ 1.71 (super-linear)
- Interpretation: species accumulate sublinearly despite super-linear total abundance

### Q2: Composition — Do larger corals support more distinct communities?

**Hypotheses tested:**
- **Site effects**: Reef-scale differences in species pools
- **Size-dependent divergence**: Larger corals accumulate unique species (not just more individuals)

**Methods:**
- PERMANOVA: community ~ site × size × neighborhood
- Composition divergence: betadisper (distance to centroid) by size class
- Rarefaction control: re-test divergence on abundance-equalized data
- NMDS ordination with size-class trajectories

**Key findings:**
- Strong site effects on composition (PERMANOVA R² ~ 0.04-0.06, p < 0.01)
- Raw analysis: community distinctness decreases with coral size (β < 0, p < 0.05)
- **After rarefaction (n=105, depth=5): trend NOT significant (p=0.61)**
- Conclusion: size-composition pattern is an abundance artifact, not true divergence

### Q3: Feedbacks — Does CAFI identity predict coral condition?

**Hypotheses tested:**
- **Mutualism**: Beneficial CAFI (crabs, shrimp) improve coral condition
- **Community identity**: PC1_CAFI → PC1_Coral (composition matters beyond abundance)
- **Key species**: Individual species effects on condition (matching experimental results)

**Methods:**
- PCA on CAFI abundances → PC1_CAFI; PCA on physiology → PC1_Coral
- Mixed models: condition ~ CAFI predictors + (1|site)
- FDR correction (Benjamini-Hochberg) across predictor families
- Key species models: condition ~ species presence + log(volume) + (1|site)

**Key findings:**
- **Species richness → condition**: p = 0.008, p_FDR = 0.053 (marginal, strongest signal)
- Total abundance and PC1_CAFI do NOT significantly predict condition after FDR
- Key species: *Trapezia* shows positive trend but NS; *Halichoeres* positive (marginal)
- Richness-condition link suggests **community complexity** may matter for coral health
- No predictor from condition→CAFI direction survives FDR correction

### Q4: Neighborhood — Does local coral density affect CAFI recruitment?

**Hypotheses tested:**
- **Neighborhood recruitment**: More neighbors → more CAFI (spillover/source-sink)
- **Neighborhood condition**: More neighbors → better coral condition (facilitation)

**Methods:**
- NB GLM: `total_cafi ~ log_volume * n_neighbors + site`
- LMM: `condition_score ~ log_volume * n_neighbors + (1|site)`
- Available on 61/114 corals (5m survey subset)

**Key findings:**
- n_neighbors NOT significant for CAFI abundance (p = 0.86)
- n_neighbors NOT significant for coral condition (p = 0.93)
- Volume remains the dominant predictor in all models
- **Neighborhood density does not explain CAFI or condition variation**

---

## Key Data Files

| File | Description | Records |
|------|-------------|---------|
| `data/survey_cafi_data_w_taxonomy_summer2019_v5.csv` | CAFI specimen records with taxonomy | 3,989 rows |
| `data/survey_coral_characteristics_merged_v2.csv` | Coral colony attributes (size, position, neighbors) | 114 rows |
| `data/survey_master_phys_data_v3.csv` | Coral physiology measurements | 108 rows |
| `data/traits/cafi_traits_final.csv` | CAFI trait database (size, depth, trophic) | varies |

**Join key**: `coral_id` (format: "SITE-POC##", e.g., "HAU-POC29")

## Pipeline Scripts

Located in `scripts/` folder. Run in order via `run_full_pipeline.R`.

### Setup & Data (Scripts 00-01)
| Script | Purpose |
|--------|---------|
| `00_setup.R` | Packages, paths, ggplot theme, helper functions |
| `00b_data_quality_audit.R` | Data quality assessment, trait integration |
| `01_load_data.R` | Data loading, cleaning, creates core objects |

### Script-to-Question Mapping

#### Q1: SCALING
| Script | Analysis | Output |
|--------|----------|--------|
| `05_species_scaling_analysis.R` | NB GLM scaling with bootstrap CI | Field of Dreams vs Redistribution |

#### Q2: COMPOSITION
| Script | Analysis | Output |
|--------|----------|--------|
| `02_community_analysis.R` | PERMANOVA, NMDS, betadisper, rarefaction | Site/size effects on composition |
| `03_landscape_characterization.R` | Neighborhood metrics, spatial patterns | Landscape predictor variables |
| `06_network_analysis.R` | Co-occurrence networks, modularity | Non-random assembly patterns |
| `07_spatial_autocorrelation.R` | Moran's I, LISA, Mantel tests | Spatial structure diagnostics |
| `08_functional_groups.R` | Trapezia, fish, gastropod scaling patterns | Taxonomic group drivers |

#### Q3: FEEDBACKS
| Script | Analysis | Output |
|--------|----------|--------|
| `09_cafi_condition_feedbacks.R` | PCA, mixed models, FDR correction | CAFI-condition relationships |

#### Q4: NEIGHBORHOOD
| Script | Analysis | Output |
|--------|----------|--------|
| `04_landscape_effects.R` | GLMs: size + neighbors → abundance/diversity | Neighborhood effect sizes |
| `09_cafi_condition_feedbacks.R` (Part G) | Neighborhood density → CAFI & condition | Density-independent result |

### Supporting Analyses
| Script | Purpose |
|--------|---------|
| `10_feature_engineering.R` | Feature creation, VIF selection |
| `11_machine_learning.R` | Random Forest, XGBoost models |
| `12_model_evaluation.R` | Cross-validation, diagnostics |
| `13_manuscript_figures.R` | Publication figures organized by Q1-Q4 |

### Publication
| Script | Purpose |
|--------|---------|
| `13_manuscript_figures.R` | Publication figures |
| `run_full_pipeline.R` | Pipeline orchestrator with logging |

**Note**: `scripts/archive/` contains deprecated scripts for reference only.

## Key Commands

### Run Full Pipeline
```r
source("scripts/run_full_pipeline.R")
run_pipeline()  # Runs all scripts with logging
```

### Run From Command Line
```bash
Rscript scripts/run_full_pipeline.R
```

### Run Individual Script
```r
source("scripts/00_setup.R")
source("scripts/01_load_data.R")
source("scripts/04_landscape_effects.R")  # Example
```

## Analytical Quality Controls

The following quality measures are implemented:

| Issue | Fix | Script |
|-------|-----|--------|
| Multiple testing (feedbacks) | FDR correction (Benjamini-Hochberg) | `09_cafi_condition_feedbacks.R` |
| Abundance confound (composition) | Rarefaction + re-test | `02_community_analysis.R` |
| Bootstrap ignoring site | Site-stratified bootstrap | `05_species_scaling_analysis.R` |
| Log-log intercept bias | Mean of log-differences | `04_landscape_effects.R` |
| Network analysis hard failure | Graceful warning + skip | `06_network_analysis.R` |
| NB convergence failure | Poisson fallback with logging | `04_landscape_effects.R` |
| Colorblind inaccessibility | Blue/orange/teal palette | All figure scripts |

## Output Structure

```
output/
├── figures/
│   ├── manuscript/          # Publication figures (fig1-6)
│   ├── supplement/          # Supplementary figures (figS1-4)
│   ├── 02_community/        # Community analysis figures
│   ├── 03_landscape/        # Landscape characterization
│   ├── 04_effects/          # Landscape effects on CAFI
│   ├── 05_scaling/          # Species-area scaling
│   ├── 06_network/          # Network analysis
│   ├── 07_spatial/          # Spatial autocorrelation
│   ├── 10_features/         # Feature engineering
│   ├── 12_evaluation/       # Model evaluation
│   ├── feedbacks/           # CAFI-condition feedbacks
│   ├── functional_groups/   # Functional group analysis
│   └── machine_learning/    # ML model outputs
├── tables/                  # CSV outputs
├── objects/                 # RDS files (coral_master, cafi_clean, models)
└── pipeline.log             # Execution log
```

## Code Conventions

- **Join key**: `coral_id` links all datasets
- **Site codes**: HAU (Hauru), MAT (Maatea), MRB (Barrier Reef)
- **Volume**: Use `volume_field` (field estimate)
- **Packages**: Use `dplyr::select()` explicitly (MASS conflict)
- **Colors**: Colorblind-safe palette (Okabe-Ito derivatives): `#0072B2` (blue), `#D55E00` (vermillion), `#009E73` (teal), `#E69F00` (orange), `#56B4E9` (sky blue)

## Load Pre-computed Objects
```r
coral_master <- load_object("coral_master")
cafi_clean <- load_object("cafi_clean")
community_matrix <- load_object("community_matrix")
```
