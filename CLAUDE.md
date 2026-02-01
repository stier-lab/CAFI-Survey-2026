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
| **Q1: SCALING** | `05_species_scaling_analysis.R` | Abundance β=0.52 (sublinear, Redirection); Richness z=0.35 (sublinear); density dilution |
| **Q2: COMPOSITION** | `02_community_analysis.R` + `06_network_analysis.R` | Site pools structure composition; co-occurrence reveals non-random modular assembly |
| **Q3: FEEDBACKS** | `09_cafi_condition_feedbacks.R` | **No CAFI metric predicts condition**; raw richness signal is abundance artifact (rarefied p=0.45) |
| **Q4: NEIGHBORHOOD** | `04_landscape_effects.R` | n_neighbors NOT significant for CAFI or condition |

---

### Q1: Scaling — How does CAFI scale with coral size?

**Hypotheses tested:**
- **Field of Dreams (β = 1)**: Abundance scales proportionally—"if you build it, they will come"
- **Propagule Redirection (β < 1)**: Abundance scales sublinearly—larger corals dilute per-capita colonization
- **Super-linear (β > 1)**: Larger corals attract disproportionately more settlers

**Methods:**
- Negative binomial GLM: `total_cafi ~ log(volume) + site` (natural log; coefficient = power-law exponent)
- Poisson GLM for species richness: `richness ~ log(volume) + site`
- Bootstrap 95% CI on scaling exponent (2,000 iterations, stratified by site)
- Species-level scaling for 10 most prevalent CAFI taxa

**Key findings:**
- Total CAFI abundance: β = 0.52 [0.44, 0.61] — **sublinear (Propagule Redirection)**
- Species richness (SAR): z = 0.35 [0.28, 0.43] — **sublinear (Redirection)**
- Density dilution: per-capita CAFI density decreases with size (slope = -0.48)
- 7/10 top species: Redirection (β < 1); 3/10: Field of Dreams (CI spans 1); 0/10: super-linear
- Obligate symbionts (Trapezia, Alpheus): consistently sublinear scaling
- Interpretation: finite propagule pools diluted across larger corals; both abundance and richness support Redirection

**NOTE**: Previous reports of β = 1.20 were based on a log-base mismatch: the model used log10(volume) as predictor with a natural-log link, inflating coefficients by ln(10) ≈ 2.303. Corrected 2026-01-30.

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
- **Position correction**: Physiological traits regressed on stump_length + nubbin_length before PCA (removes tissue gradient sampling artifact: smaller corals required sampling more of branch, including low-density tips)
- Linear models with fixed-effect site: condition ~ CAFI predictors + log_volume + site
- FDR correction (Benjamini-Hochberg) across predictor families
- Key species models: condition ~ species presence + log(volume) + site (FDR-corrected)
- **Rarefied richness test**: Expected species richness at n=20 individuals (removes abundance confound)
- Standardized β (effect in SD units); VIF diagnostics
- Note: 3 sites insufficient for random intercepts (Bolker et al. 2009)

**Key findings:**
- **Species richness → condition**: Raw richness shows nominal signal (p = 0.018), BUT this is an **ABUNDANCE ARTIFACT**
  - Richness correlates strongly with total CAFI abundance (r = 0.84)
  - **Rarefied richness** (controlling for sampling intensity) shows NO relationship (p = 0.45)
  - Correlation with abundance: raw r = 0.84, rarefied r = −0.05
  - Conclusion: apparent diversity effect disappears when abundance confound is removed
- Total abundance, PC1_CAFI, and individual CAFI metrics do NOT predict condition (all p > 0.10)
- Key species: *Trapezia* shows positive trend but NS; *Galeropsis* positive (not negative as expected)
- Functional group effect directions match expectations (Trapezia +, Fish +) but none are significant
- No predictor from condition→CAFI direction survives FDR correction
- **Landscape factors (size, neighborhood, site) also do NOT predict condition** (all p > 0.30)

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
- n_neighbors NOT significant for CAFI abundance (p = 0.37)
- n_neighbors NOT significant for coral condition (p = 0.78)
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
| `05_species_scaling_analysis.R` | NB GLM scaling with bootstrap CI | Field of Dreams vs Redirection |

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

### Publication Figures (embedded in analysis scripts — no separate figure script)

Each manuscript figure is created by its source analysis script with **dual saves**: once to the analysis output directory, once to `output/figures/manuscript/` or `output/figures/supplement/`.

| Figure | Script | Description |
|--------|--------|-------------|
| **Fig 1** | `01_load_data.R` | Study design: satellite map + volume/neighborhood histograms |
| **Fig 2** | `05_species_scaling_analysis.R` | Scaling: abundance + richness GLMs + species forest plot |
| **Fig 3** | `08_functional_groups.R` | Taxonomic group scaling and composition |
| **Fig 4** | `06_network_analysis.R` | Network: circular hero layout + 4 guild sub-networks |
| **Fig 5** | `09_cafi_condition_feedbacks.R` | CAFI-condition feedback models |
| **S1** | `02_community_analysis.R` | Species accumulation curves |
| **S2** | `02_community_analysis.R` | PERMANOVA metric sensitivity |
| **S3** | `02_community_analysis.R` | NMDS ordination by site/size |
| **S4** | `07_spatial_autocorrelation.R` | Spatial autocorrelation (Moran's I) |
| **S5** | `02_community_analysis.R` | Composition divergence by size |
| **S6** | `05_species_scaling_analysis.R` | Species-level scaling forest plot |
| **S7** | `04_landscape_effects.R` | Neighborhood null results |

| Script | Purpose |
|--------|---------|
| `run_full_pipeline.R` | Pipeline orchestrator with logging |

**Note**: `scripts/archive/` contains deprecated scripts (including former `13_manuscript_figures.R`) for reference only.

## Key Commands

### Run Pipeline
```r
source("scripts/run_full_pipeline.R")
run_pipeline()           # Core + extended (default; figures created by each script)
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
| Multiple testing (scaling) | FDR correction within category (species/group) | `05_species_scaling_analysis.R` |
| Multiple testing (network edges) | Pairwise cor.test p-values + FDR correction | `06_network_analysis.R` |
| Abundance confound (composition) | Iterated rarefaction (100 draws, averaged) | `02_community_analysis.R` |
| Volume confound (co-occurrence) | Residualized presence on log(volume) | `06_network_analysis.R` |
| Random effects (k=3 sites) | Fixed-effect site throughout | `04`, `09` |
| Bootstrap ignoring site | Site-stratified bootstrap (`strata` argument) | `05_species_scaling_analysis.R` |
| Log-log intercept bias | Mean of log-differences | `04_landscape_effects.R` |
| Null model (network) | Configuration model (degree-preserving), unweighted for comparison | `06_network_analysis.R` |
| NB convergence failure | Poisson fallback with logging | `04_landscape_effects.R` |
| Effect size ambiguity | Adjusted R², partial standardized β (z-scored), VIF | `04`, `05`, `09` |
| Model diagnostics (GLM) | DHARMa simulated residuals, proper Cook's D, VIF | `04`, `05` |
| Model diagnostics (LM) | Shapiro-Wilk, Breusch-Pagan via HC3 robust SEs | `09` |
| Poisson overdispersion | Pearson X²/df check; auto-switch to quasipoisson | `02`, `04` |
| Heteroscedasticity (count predictors) | HC3 sandwich robust standard errors | `09_cafi_condition_feedbacks.R` |
| SEM saturation | Explicit df=0 check; fit indices suppressed when uninformative | `09_cafi_condition_feedbacks.R` |
| Community matrix (zero-CAFI corals) | Added zero-abundance rows for all corals | `01_load_data.R` |
| Log volume bias | Removed +1 offset (volume > 0 guaranteed) | `01_load_data.R` |
| Tissue sampling artifact (condition) | Regress physio traits on stump_length + nubbin_length | `01_load_data.R` |
| PERMANOVA Type I confound | Marginal (Type III) PERMANOVA | `02_community_analysis.R` |
| Colorblind inaccessibility | Okabe-Ito palette; Fig 2 site colors shifted to avoid scaling-class overlap | All figure scripts |

## Output Structure

```
output/
├── figures/
│   ├── manuscript/          # Main text figures (fig1-fig5)
│   │   ├── fig1_study_design.png
│   │   ├── fig2_scaling.png
│   │   ├── fig3_functional_groups.png
│   │   ├── fig4_network.png
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
- **Volume**: Use `volume` (field estimate)
- **Key columns**: `n_galeropsis` (Galeropsis count per coral), `n_corallivore` (all gastropods)
- **Packages**: Use `dplyr::select()` explicitly (MASS conflict); car::vif(), DHARMa, sandwich/lmtest for diagnostics
- **Colors**: Two semantic palettes, chosen to avoid cross-figure confusion:

  **Site palette** (`SITE_COLORS` in `00_setup.R`; used globally across all figures):
  | Site | Hex | Description |
  |------|-----|-------------|
  | HAU | `#9B7EB8` | Muted purple — Hauru (fringing) |
  | MAT | `#7B9BAE` | Cool slate — Maatea (back-reef) |
  | MRB | `#7AAC6D` | Sage green — Barrier Reef |

  **Scaling-class palette** (Figure 2 Panel C):
  | Class | Hex | Usage |
  |-------|-----|-------|
  | Redirection (β < 1) | `#5A8FAF` | Muted blue |
  | Field of Dreams (CI overlaps 1) | `gray55` | Neutral — consistent with null |
  | Super-linear (CI > 1) | `#D55E00` | Vermillion accent — distinctive result |

  Site palette deliberately avoids blue and orange/vermillion to prevent confusion with scaling-class semantics.

## Load Pre-computed Objects
```r
coral_master <- load_object("coral_master")
cafi_clean <- load_object("cafi_clean")
community_matrix <- load_object("community_matrix")
```
