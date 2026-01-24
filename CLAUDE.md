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
| **Q2: COMPOSITION** | `02_community_analysis.R` + `06_network_analysis.R` | Site pools structure composition; co-occurrence reveals non-random modular assembly |
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

### Q2: Composition — What structures CAFI community composition?

**Hypotheses tested:**
- **Site effects**: Reef-scale differences in species pools
- **Size-dependent divergence**: Larger corals accumulate unique species (not just more individuals)
- **Non-random co-occurrence**: Species co-occur in modular networks with identifiable hub species

**Methods:**
- PERMANOVA: community ~ site × size × neighborhood
- Composition divergence: betadisper (distance to centroid) by size class
- Rarefaction control: re-test divergence on abundance-equalized data
- NMDS ordination with size-class trajectories
- Co-occurrence network (r > 0.3, Louvain modularity, Erdos-Renyi null model)
- Hub species identification (degree + eigenvector centrality)

**Key findings:**
- Strong site effects on composition (PERMANOVA R² ~ 0.04-0.06, p < 0.01)
- Size-divergence NOT significant after rarefaction (p=0.61; Fig. S5) — abundance artifact
- Co-occurrence network: modular structure significantly exceeds null expectation (z >> 2)
- Hub species connect modules and may serve as keystone structuring agents
- Conclusion: site pools and non-random co-occurrence shape composition, not coral size

### Q3: Feedbacks — Does CAFI identity predict coral condition?

**Hypotheses tested:**
- **Mutualism**: Beneficial CAFI (crabs, shrimp) improve coral condition
- **Community identity**: PC1_CAFI → PC1_Coral (composition matters beyond abundance)
- **Key species**: Individual species effects on condition (matching experimental results)

**Methods:**
- PCA on CAFI abundances → PC1_CAFI; PCA on physiology → PC1_Coral
- Linear models with fixed-effect site: condition ~ CAFI predictors + log_volume + site
- FDR correction (Benjamini-Hochberg) across predictor families
- Key species models: condition ~ species presence + log(volume) + site (FDR-corrected)
- Standardized β (effect in SD units); VIF diagnostics
- Note: 3 sites insufficient for random intercepts (Bolker et al. 2009)

**Key findings:**
- **Species richness → condition**: p = 0.008, p_FDR = 0.053 (marginal, strongest signal)
- Total abundance and PC1_CAFI do NOT significantly predict condition after FDR
- Key species: *Trapezia* shows positive trend but NS; *Galeropsis* negative (tissue consumer)
- Richness-condition link suggests **community complexity** may matter for coral health
- No predictor from condition→CAFI direction survives FDR correction

**Note on *Galeropsis monodonta***: This coralliophiline snail (Muricidae) dominates gastropod assemblages (73% of all gastropods = 356/489 individuals). It feeds on coral tissue (subfamily Coralliophilinae). Used as species-level predictor `n_galeropsis` instead of all gastropods combined.

### Q4: Neighborhood — Does local coral density affect CAFI recruitment?

**Hypotheses tested:**
- **Neighborhood recruitment**: More neighbors → more CAFI (spillover/source-sink)
- **Neighborhood condition**: More neighbors → better coral condition (facilitation)

**Methods:**
- NB GLM: `total_cafi ~ log_volume * n_neighbors + site`
- LM: `condition_score ~ log_volume * n_neighbors + site` (fixed effect)
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

#### Q2: COMPOSITION (Fig 4: NMDS + Network + Modularity + Hubs)
| Script | Analysis | Output |
|--------|----------|--------|
| `02_community_analysis.R` | PERMANOVA, NMDS, betadisper, rarefaction | Site/size effects on composition |
| `06_network_analysis.R` | Co-occurrence networks, modularity, hub species | Fig 4 panels B-D |
| `03_landscape_characterization.R` | Neighborhood metrics, spatial patterns | Landscape predictor variables |
| `07_spatial_autocorrelation.R` | Moran's I, LISA, Mantel tests | Spatial structure diagnostics |
| `08_functional_groups.R` | Taxonomic group scaling (loads from script 05), composition | Fig 3 + group patterns |

#### Q3: FEEDBACKS
| Script | Analysis | Output |
|--------|----------|--------|
| `09_cafi_condition_feedbacks.R` | PCA, fixed-effect LMs, FDR correction | CAFI-condition relationships |

#### Q4: NEIGHBORHOOD
| Script | Analysis | Output |
|--------|----------|--------|
| `04_landscape_effects.R` | GLMs: size + neighbors → abundance/diversity | Neighborhood effect sizes |
| `09_cafi_condition_feedbacks.R` (Part G) | Neighborhood density → CAFI & condition | Density-independent result |

### Exploratory ML (not in default pipeline)
| Script | Purpose |
|--------|---------|
| `10_feature_engineering.R` | Feature creation, VIF selection |
| `11_machine_learning.R` | Random Forest, XGBoost models |
| `12_model_evaluation.R` | Cross-validation, diagnostics |

### Publication
| Script | Purpose |
|--------|---------|
| `13_manuscript_figures.R` | Publication figures (5 main + 7 supplementary), results table, sample sizes |
| `run_full_pipeline.R` | Pipeline orchestrator with logging |

**Note**: `scripts/archive/` contains deprecated scripts for reference only.

## Key Commands

### Run Pipeline
```r
source("scripts/run_full_pipeline.R")
run_pipeline()           # Core + extended + publication (default)
run_full_pipeline()      # Everything including ML exploration
run_ml_exploration()     # ML scripts only (exploratory)
run_core_pipeline()      # Core hypothesis-testing scripts only
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
| Multiple testing (key species) | FDR correction across 6 species tests | `09_cafi_condition_feedbacks.R` |
| Abundance confound (composition) | Iterated rarefaction (100 draws, averaged) | `02_community_analysis.R` |
| Volume confound (co-occurrence) | Residualized presence on log(volume) | `06_network_analysis.R` |
| Random effects (k=3 sites) | Fixed-effect site throughout | `04`, `09` |
| Bootstrap ignoring site | Site-stratified bootstrap | `05_species_scaling_analysis.R` |
| Log-log intercept bias | Mean of log-differences | `04_landscape_effects.R` |
| Null model (network) | Configuration model (degree-preserving) | `06_network_analysis.R` |
| NB convergence failure | Poisson fallback with logging | `04_landscape_effects.R` |
| Effect size ambiguity | Adjusted R², standardized β, VIF | `04`, `05`, `09` |
| Model diagnostics | VIF, Shapiro-Wilk, Cook's D | `09`, `04` |
| PERMANOVA Type I confound | Marginal (Type III) PERMANOVA | `02_community_analysis.R` |
| Colorblind inaccessibility | Okabe-Ito palette | All figure scripts |

## Output Structure

```
output/
├── figures/
│   ├── manuscript/          # Main text figures (fig1-fig5)
│   │   ├── fig1_study_design.png
│   │   ├── fig2_scaling.png
│   │   ├── fig3_functional_groups.png
│   │   ├── fig4_composition_network.png
│   │   └── fig5_feedbacks.png
│   ├── supplement/          # Supplementary figures (figS1-S7)
│   ├── 02_community/        # Community analysis figures
│   ├── 03_landscape/        # Landscape characterization
│   ├── 04_effects/          # Landscape effects on CAFI
│   ├── 05_scaling/          # Species-area scaling
│   ├── 06_network/          # Network analysis
│   ├── 07_spatial/          # Spatial autocorrelation
│   ├── feedbacks/           # CAFI-condition feedbacks
│   └── functional_groups/   # Taxonomic group analysis
├── tables/
│   ├── manuscript_results_summary.csv  # All Q1-Q4 results
│   ├── sample_sizes.csv               # N for each analysis
│   └── ...                            # Other CSV outputs
├── objects/                 # RDS files (coral_master, cafi_clean, models)
└── pipeline.log             # Execution log
```

## Code Conventions

- **Join key**: `coral_id` links all datasets
- **Site codes**: HAU (Hauru), MAT (Maatea), MRB (Barrier Reef)
- **Volume**: Use `volume_field` (field estimate)
- **Key columns**: `n_galeropsis` (Galeropsis count per coral), `n_corallivore` (all gastropods)
- **Packages**: Use `dplyr::select()` explicitly (MASS conflict); car::vif() for diagnostics
- **Colors**: Colorblind-safe palette (Okabe-Ito derivatives): `#0072B2` (blue), `#D55E00` (vermillion), `#009E73` (teal), `#E69F00` (orange), `#56B4E9` (sky blue)

## Load Pre-computed Objects
```r
coral_master <- load_object("coral_master")
cafi_clean <- load_object("cafi_clean")
community_matrix <- load_object("community_matrix")
```
