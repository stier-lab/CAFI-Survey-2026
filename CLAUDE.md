# CLAUDE.md - AI Assistant Context

This file provides essential context for AI assistants working with this codebase.

## Project Overview

**CAFI Survey Analysis** - A marine ecology research project studying how coral size and neighborhood context shape coral-associated fauna (CAFI) communities in Mo'orea, French Polynesia. **Target journal: Journal of Animal Ecology** — write for a broad ecological audience (community ecology, BEF, landscape ecology, conservation), using the coral system as a case study to illustrate general principles.

- **Survey scope**: 114 *Pocillopora* coral colonies across 3 reef sites
- **CAFI catalogued**: ~4,000 individual specimens spanning 243 OTUs (154 species-level, 89 genus/family/higher)
- **Research focus**: How coral size structures CAFI abundance, richness, and composition (scaling + community assembly)

## Relationship to Experimental Study

This **observational survey** complements a parallel **experimental study** that manipulated coral density (1, 3, 6 corals per 25m² reef). The two studies test the same hypotheses using different approaches:

| Approach | Lever | Strength |
|----------|-------|----------|
| **Experiment** | Coral density (manipulated) | Causal inference, colonization dynamics |
| **Survey** | Coral size + neighborhood (natural variation) | Ecological realism, established communities |

## Three Core Research Questions (Q1-Q3)

| Question | Script | Key Result |
|----------|--------|------------|
| **Q1: SCALING** | `05_species_scaling_analysis.R` | Abundance β=0.52 (sublinear, Redirection); Richness z=0.34 (sublinear); density dilution |
| **Q2: COMPOSITION** | `02_community_analysis.R` | Site pools + size gradient structure composition; db-RDA confirms size effect survives rarefaction |
| **Q3: FEEDBACKS** | `09_cafi_condition_feedbacks.R` | Species richness predicts condition (BEF framework, p=0.018); variance partitioning: 29% unique to richness, <1% to abundance |

### Supporting Analyses (Supplement)

| Analysis | Script | Key Result |
|----------|--------|------------|
| **Co-occurrence** | `06_cooccurrence_analysis.R` | Volume-weighted null model; pairwise associations mostly explained by volume + site (Fig S9) |
| **Neighborhood** | `04_landscape_effects.R` | n_neighbors NOT significant (underpowered, n=61); supplement only |

---

### Q1: Scaling — How does CAFI scale with coral size?

**Hypotheses tested:**
- **Field of Dreams (β = 1)**: Abundance scales proportionally—"if you build it, they will come"
- **Propagule Redirection (β < 1)**: Abundance scales sublinearly—larger corals dilute per-capita colonization
- **Super-linear (β > 1)**: Larger corals attract disproportionately more settlers

**Methods:**
- Negative binomial GLM: `total_cafi ~ log(volume) + site` (natural log; coefficient = power-law exponent)
- Poisson GLM for species richness: `richness ~ log(volume) + site`
- Bootstrap 95% CI on scaling exponent (1,000 iterations, stratified by site)
- Species-level scaling for 21 most prevalent CAFI species

**Key findings:**
- Total CAFI abundance: β = 0.52 [0.44, 0.62] — **sublinear (Propagule Redirection)**
- Species richness (SAR): z = 0.34 [0.27, 0.42] — **sublinear (Redirection)**
- Rarefied richness (n=20): slope=−0.07, p=0.50 — **not significant (passive sampling)**; raw SAR is abundance artifact
- Species occurrence curves: 14/24 prevalent species have significant size-dependent occurrence (FDR<0.05)
- Density dilution: per-capita CAFI density decreases with size (slope = -0.48)
- 11/21 species: Redirection (β < 1); 10/21: Field of Dreams (CI spans 1); 0/21: super-linear
- Obligate symbionts (Trapezia, Alpheus): consistently sublinear scaling
- Interpretation: finite propagule pools diluted across larger corals; both abundance and richness support Redirection

**NOTE**: Previous reports of β = 1.20 were based on a log-base mismatch: the model used log10(volume) as predictor with a natural-log link, inflating coefficients by ln(10) ≈ 2.303. Corrected 2026-01-30.

### Q2: Composition — What structures CAFI community composition?

**Hypotheses tested:**
- **Site effects**: Reef-scale differences in species pools
- **Size-dependent divergence**: Larger corals accumulate unique species (not just more individuals)
- **Non-random co-occurrence**: Pairwise null-model co-occurrence (Stier et al. 2012 framework)

**Methods:**
- PERMANOVA: community ~ log(volume) + site (marginal/Type III)
- Composition divergence: betadisper (distance to centroid) by size class
- Rarefaction control: re-test divergence on abundance-equalized data
- NMDS ordination with size-class trajectories
- db-RDA constrained ordination: `dbrda(community ~ log(volume) + Condition(site))`, variance partitioning, species scores via weighted averaging, rarefied control
- Pairwise co-occurrence null model: volume-weighted Bernoulli draws (10,000 iterations, SES + FDR correction)
- Intraspecific density patterns: multinomial allocation null (volume-proportional, 10,000 iterations)
- Size-dependent co-occurrence: null model run separately for 3 size classes

**Key findings:**
- Strong site effects on composition (PERMANOVA site R² ≈ 0.06, volume R² ≈ 0.08, p < 0.01; varies by permutation run)
- db-RDA: volume explains 7.8% of composition (F=9.74, p=0.001); survives rarefaction (2.4%, p=0.001); *T. punctimanus* loads most strongly toward smaller corals (score = -2.95)
- Size-divergence (categorical) NOT significant after rarefaction (p=0.61; Fig. S5) — abundance artifact
- Nestedness (NODF): 18.37, z=−1.09, p=0.277 — **not nested**; small-coral faunas are not subsets of large-coral faunas
- Pairwise co-occurrence: 0 of 528 pairs significant after FDR correction; strongest signals: H. beaupresii–P. modestus (SES = −3.43, p_FDR = 0.11), H. beaupresii–H. spinigera (SES = −3.42, p_FDR = 0.21)
- Intraspecific density patterns tested for mating-pair hypothesis (Stier et al. 2012)
- Conclusion: site pools and continuous size gradient shape composition; co-occurrence is largely explained by volume + site

### Q3: Feedbacks — Does CAFI diversity predict coral condition? (BEF framework)

**Hypotheses tested:**
- **BEF complementarity (a priori)**: More CAFI species → better condition (Tilman diversity-function)
- **Abundance pathway (a priori)**: More CAFI individuals → better condition
- **Community identity (supplement)**: PC1_CAFI → PC1_Coral (composition matters)
- **Key species**: Individual species effects on condition (matching experimental results)

**Methods:**
- PCA on CAFI abundances → PC1_CAFI; PCA on physiology → PC1_Coral
- **Position correction**: Physiological traits regressed on stump_length + nubbin_length before PCA (removes tissue gradient sampling artifact)
- Linear models with fixed-effect site: condition ~ CAFI predictors + log_volume + site
- **Multiple testing**: A priori BEF predictors (richness, abundance) tested without correction (pre-specified hypotheses); exploratory (k=4) BH-FDR; PC1_CAFI tested separately in supplement
- **OLS standard errors** (primary; Breusch-Pagan confirms homoscedasticity, BP p > 0.5)
- **HC3 robust SEs** reported as supplement sensitivity (conservative at n < 100; Long & Ervin 2000)
- **BEF analysis (Part A4)**: Partial regression, variance partitioning, path model (piecewiseSEM)
- Key species models: condition ~ species abundance + log(volume) + site (FDR-corrected, up to 10 species; those with n_present < 5 excluded at runtime)
- Note: 3 sites insufficient for random intercepts (Bolker et al. 2009)

**Key findings:**
- **Species richness → condition**: p = 0.018 — **SIGNIFICANT** (pre-specified a priori hypothesis)
  - Richness correlates strongly with total CAFI abundance (r = 0.84)
  - **Rarefied richness** (n=20) shows NO relationship (p = 0.50), but this test is **AMBIGUOUS**:
    rarefaction may remove abundance artifact OR the BEF mechanism itself (diversity→abundance→condition)
- **Total CAFI abundance → condition**: p = 0.048 — marginal (pre-specified a priori hypothesis)
- **Variance partitioning**: 29.1% unique to richness, <1% unique to abundance, 70.8% shared
- **Path model**: β(richness→condition) = 0.55, β(abundance→condition) = 0.02
- PC1_CAFI does NOT predict condition (supplement; p > 0.10)
- Exploratory functional groups: directions match expectations (Trapezia +, Fish +) but none survives BH-FDR
- Key species (10 tested): directions mostly match experiment but none survive FDR correction
- Reverse direction (condition→CAFI): non-significant
- **Landscape factors (size, neighborhood, site) also do NOT predict condition** (all p > 0.30)

**Note on *Galeropsis monodonta***: This coralliophiline snail (Muricidae) dominates gastropod assemblages (73% of all gastropods = 356/489 individuals). It feeds on coral tissue (subfamily Coralliophilinae). Used as species-level predictor `n_galeropsis` instead of all gastropods combined.

### Q4: Neighborhood — Does local coral density affect CAFI recruitment?

**Hypotheses tested:**
- **Neighborhood recruitment**: More neighbors → more CAFI (spillover/source-sink)
- **Neighborhood condition**: More neighbors → better coral condition (facilitation)

**Methods:**
- NB GLM: `total_cafi ~ log_volume * n_neighbors + site`
- LM: `condition_score ~ log_volume * n_neighbors + site` (fixed effect)
- Available on 61/114 corals (5m survey subset, after volume filter)

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
| `01_load_data.R` | Data loading, cleaning, creates core objects + taxonomy scenario data |
| `run_full_pipeline.R` | Pipeline orchestrator with logging |

### Script-to-Question Mapping

#### Q1: SCALING
| Script | Analysis | Output |
|--------|----------|--------|
| `05_species_scaling_analysis.R` | NB GLM scaling with bootstrap CI | Field of Dreams vs Redirection |

#### Q1 (cont.): SPECIES & GROUP SCALING (Fig 3)
| Script | Analysis | Output |
|--------|----------|--------|
| `05_species_scaling_analysis.R` | Species + taxonomic group NB GLMs, forest plots | Fig 3 (2×2: curves + β estimates) |

#### Q2: COMPOSITION (Fig 4: db-RDA + Composition)
| Script | Analysis | Output |
|--------|----------|--------|
| `02_community_analysis.R` | PERMANOVA, NMDS, betadisper, rarefaction, db-RDA | Fig 4 + site/size effects |
| `03_landscape_characterization.R` | Neighborhood metrics, spatial patterns | Landscape predictor variables |
| `07_spatial_autocorrelation.R` | Moran's I, LISA, Mantel tests | Spatial structure diagnostics |
| `08_functional_groups.R` | Taxonomic group scaling (loads from script 05), composition | Fig S7 + group patterns |

#### Q3: FEEDBACKS (Fig 5)
| Script | Analysis | Output |
|--------|----------|--------|
| `09_cafi_condition_feedbacks.R` | PCA, fixed-effect LMs, FDR correction | Fig 5 (BEF richness + abundance scatter) |

#### SUPPORTING (Supplement)
| Script | Analysis | Output |
|--------|----------|--------|
| `06_cooccurrence_analysis.R` | Null-model pairwise co-occurrence, intraspecific density | Fig S9 (co-occurrence) |
| `04_landscape_effects.R` | GLMs: size + neighbors → abundance/diversity | Fig S5 (neighborhood null) |

#### SENSITIVITY (cross-cuts Q1-Q4)
| Script | Analysis | Output |
|--------|----------|--------|
| `13_taxonomy_sensitivity.R` | 5 taxonomy scenarios × 7 metrics; uses pre-built data from 01_load_data.R | Fig S6 + sensitivity tables |

### Archived Scripts (`scripts/archive/` — NOT part of the manuscript or pipeline)

These scripts are retained for reference only. Do not run, edit, or cite them when building the paper.

| Script | Purpose |
|--------|---------|
| `archive/06_network_analysis_legacy.R` | Legacy network analysis (replaced by `06_cooccurrence_analysis.R`) |
| `archive/10_feature_engineering.R` | Exploratory ML feature creation, VIF selection |
| `archive/11_machine_learning.R` | Exploratory Random Forest, XGBoost models |
| `archive/12_model_evaluation.R` | Exploratory cross-validation, diagnostics |

### Publication Figures (embedded in analysis scripts — no separate figure script)

Each manuscript figure is created by its source analysis script with **dual saves**: once to the analysis output directory, once to `output/figures/manuscript/` or `output/figures/supplement/`.

| Figure | Script | Description |
|--------|--------|-------------|
| **Fig 1** | `01_load_data.R` | Study design + representative CAFI: satellite map + distributions + species photos (6-panel) |
| **Fig 2** | `05_species_scaling_analysis.R` | Scaling: (A) abundance + (B) density dilution + (C) richness (3 vertical panels) |
| **Fig 3** | `05_species_scaling_analysis.R` | Species + taxonomic group scaling (2×2: curves + β forest plots) |
| **Fig 4** | `02_community_analysis.R` | Composition: (A) db-RDA biplot + (B) taxonomic barchart by site |
| **Fig 5** | `09_cafi_condition_feedbacks.R` | BEF diversity-condition: (A) richness scatter + (B) abundance scatter (2-panel) |
| **S1** | `02_community_analysis.R` | Species accumulation curves |
| **S2** | `02_community_analysis.R` | PERMANOVA metric sensitivity |
| **S3** | `02_community_analysis.R` | Composition divergence by size |
| **S4** | `05_species_scaling_analysis.R` | Species-level scaling forest plot |
| **S5** | `04_landscape_effects.R` | Neighborhood null results |
| **S6** | `13_taxonomy_sensitivity.R` | Taxonomy sensitivity forest plot |
| **S7** | `08_functional_groups.R` | Taxonomic group scaling and composition |
| **S8** | `09_cafi_condition_feedbacks.R` | Rarefaction depth sensitivity for richness → condition |
| **S9** | `06_cooccurrence_analysis.R` | Co-occurrence: pairwise SES heatmap + intraspecific density + size-dependent (3-panel) |
| **S10** | `09_cafi_condition_feedbacks.R` | BEF variance partitioning + partial regression + path model coefficients |
| **S11** | `09_cafi_condition_feedbacks.R` | A priori forest + rarefied richness + exploratory forest + Trapezia/Galeropsis scatter + bidirectional (6-panel) |
| **S12** | `05_species_scaling_analysis.R` | Species occurrence probability vs. coral size (logistic GLM, 24 species) |
| **S13** | `09_cafi_condition_feedbacks.R` | Species × trait heatmap: 19 prevalent species (≥5 corals) × 5 condition metrics (β values + FDR p-values) |
| **S14** | `09_cafi_condition_feedbacks.R` | Species × trait biplots: 9 strongest associations (scatter + regression) |

**Note**: All archived/exploratory scripts are in `scripts/archive/` — they are NOT part of the manuscript or pipeline.

## Key Commands

### Run Pipeline
```r
source("scripts/run_full_pipeline.R")
run_pipeline()           # Full pipeline ~12 min (core + extended)
run_full_pipeline()      # Everything including archived ML scripts (if available)

# Fast iteration (use these when editing figures/code):
run_one("09")            # Single script ~5 sec (auto-loads setup + data)
run_from("06")           # Scripts 06-13 ~25 sec
run_quick()              # Skip slow scripts 02+05, ~42 sec (needs prior full run)
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
| Multiple testing (feedbacks) | A priori BEF predictors (richness, abundance) uncorrected (pre-specified hypotheses); BH-FDR (exploratory k=4); uncorrected (PC1_CAFI k=1) | `09_cafi_condition_feedbacks.R` |
| Multiple testing (scaling) | FDR correction within category (species/group) | `05_species_scaling_analysis.R` |
| Multiple testing (co-occurrence) | FDR correction across all pairwise SES tests | `06_cooccurrence_analysis.R` |
| Abundance confound (composition) | Iterated rarefaction (100 draws, averaged) | `02_community_analysis.R` |
| db-RDA rarefaction control | Rarefied db-RDA confirms size gradient is not abundance artifact | `02_community_analysis.R` |
| Volume confound (co-occurrence) | Volume-weighted Bernoulli null model (logistic GLM predicted probs) | `06_cooccurrence_analysis.R` |
| Random effects (k=3 sites) | Fixed-effect site throughout | `04`, `09` |
| Bootstrap ignoring site | Site-stratified bootstrap (`strata` argument) | `05_species_scaling_analysis.R` |
| Null model (co-occurrence) | Volume-weighted Bernoulli draws (10,000 iter); multinomial density null | `06_cooccurrence_analysis.R` |
| NB convergence failure | Poisson fallback with logging | `04_landscape_effects.R` |
| Effect size ambiguity | Adjusted R², partial standardized β (z-scored), VIF | `04`, `05`, `09` |
| Model diagnostics (GLM) | DHARMa simulated residuals, proper Cook's D, VIF | `04`, `05` |
| Model diagnostics (LM) | Shapiro-Wilk normality; Breusch-Pagan homoscedasticity (BP p > 0.5); OLS primary, HC3 supplement | `09_cafi_condition_feedbacks.R` |
| Poisson overdispersion | Pearson X²/df check; auto-switch to quasipoisson | `02`, `04` |
| Heteroscedasticity (count predictors) | OLS primary (BP confirms homoscedasticity); HC3 as supplement sensitivity | `09_cafi_condition_feedbacks.R` |
| Mediation bootstrap | 1000 bootstrap iterations via mediation::mediate() (computed but not reported in manuscript) | `09_cafi_condition_feedbacks.R` |
| Community matrix (zero-CAFI corals) | Added zero-abundance rows for all corals | `01_load_data.R` |
| Log volume bias | Removed +1 offset (volume > 0 guaranteed) | `01_load_data.R` |
| Tissue sampling artifact (condition) | Regress physio traits on stump_length + nubbin_length | `01_load_data.R` |
| PERMANOVA Type I confound | Marginal (Type III) PERMANOVA | `02_community_analysis.R` |
| Colorblind inaccessibility | Okabe-Ito palette; Fig 2 site colors shifted to avoid scaling-class overlap | All figure scripts |
| Count predictor skew | sqrt() applied to count-based CAFI predictors in condition models | `09_cafi_condition_feedbacks.R` |
| Key species FDR correction | BH-FDR for species-level tests | `09_cafi_condition_feedbacks.R` |
| Community transform sensitivity | 3 transforms + filtered PCA tested for PC1→condition | `09_cafi_condition_feedbacks.R` |
| PERMANOVA site-size robustness | Balanced site subsampling (500 iter, min_n per site) | `02_community_analysis.R` |
| Power analysis (feedbacks) | Prospective power for Q3 (n=84); medium effects detectable | `09_cafi_condition_feedbacks.R` |
| Power analysis (neighborhood) | Prospective power for Q4 (n=61); underpowered for small effects | `04_landscape_effects.R` |
| FDR family vs global | Both family-wise and global FDR reported; global as sensitivity | `09_cafi_condition_feedbacks.R` |
| BEF variance partitioning | Hierarchical R² decomposition: unique richness, unique abundance, shared | `09_cafi_condition_feedbacks.R` |
| BEF path model | piecewiseSEM DAG: volume→richness→condition, volume→abundance→condition | `09_cafi_condition_feedbacks.R` |
| BEF partial regression | Richness + abundance in same model; VIF check (VIF ≈ 6.2) | `09_cafi_condition_feedbacks.R` |
| Position correction sensitivity | Raw vs corrected condition scores compared | `09_cafi_condition_feedbacks.R` |
| Bootstrap fallback logging | BCa→percentile CI fallback logged with warning | `05_species_scaling_analysis.R` |
| Bootstrap p vs β=1 | Proportion-based p-value from bootstrap replicates (complements Wald z) | `05_species_scaling_analysis.R` |
| Cross-study sign concordance | Binomial test: species condition directions match experiment | `09_cafi_condition_feedbacks.R` |
| Cross-study power comparison | Survey powered to detect experiment's effect sizes (R²≈0.12) | `09_cafi_condition_feedbacks.R` |
| Species scaling concordance | Survey vs experiment scaling patterns for overlapping species | `05_species_scaling_analysis.R` |
| Neighborhood composition divergence | Betadisper: compositional variability by neighbor density | `04_landscape_effects.R` |

## Output Structure

All outputs are gitignored and regenerated by the pipeline.

```
output/
├── figures/
│   ├── manuscript/          # 5 main text figures + legend files (PNG + PDF dual saves via save_figure())
│   │   ├── fig1_study_design.png       (from 01_load_data.R; map + distributions + CAFI photos)
│   │   ├── fig1_legend_results.txt     (from 01_load_data.R)
│   │   ├── fig2_scaling.png            (from 05_species_scaling_analysis.R; 3 vertical panels)
│   │   ├── fig2_legend_results.txt     (from 05_species_scaling_analysis.R)
│   │   ├── fig3_species_group_scaling.png (from 05_species_scaling_analysis.R; 2×2 species + group)
│   │   ├── fig3_legend_results.txt     (from 05_species_scaling_analysis.R)
│   │   ├── fig4_composition.png        (from 02_community_analysis.R)
│   │   ├── fig4_legend_results.txt     (from 02_community_analysis.R)
│   │   ├── fig5_feedbacks.png          (from 09_cafi_condition_feedbacks.R; 2-panel BEF: richness scatter + abundance scatter)
│   │   └── fig5_legend_results.txt     (from 09_cafi_condition_feedbacks.R)
│   ├── supplement/          # Supplementary figures (figS1-S14)
│   ├── 01_data/             # Study design (1 figure)
│   ├── 02_community/        # Community analysis (~13 figures)
│   ├── 03_landscape/        # Landscape characterization (3 figures)
│   ├── 04_effects/          # Landscape effects (~12 figures + 2 HTML tables)
│   ├── 05_scaling/          # Species-area scaling (~10 figures)
│   ├── 06_network/          # Co-occurrence analysis (supplement only: figS9)
│   ├── feedbacks/           # CAFI-condition feedbacks (~13 figures)
│   └── functional_groups/   # Taxonomic group analysis (7 figures)
├── tables/                  # ~66 CSV statistical results
│   ├── scaling_results_all.csv         # All scaling coefficients
│   ├── cafi_condition_models.csv       # Feedback model results
│   ├── key_species_effects.csv         # Key species condition effects (up to 10 species)
│   ├── species_trait_correlations.csv  # Species × trait Pearson r (cf. expt Table S3)
│   ├── individual_physiology_cafi_responses.csv # Individual traits ~ CAFI predictors
│   ├── cross_study_species_comparison.csv # Survey vs experimental paper comparison
│   ├── landscape_full_model_results.csv # Q4 landscape GLMs
│   ├── network_metrics.csv             # Network structure
│   ├── hub_species.csv                 # Hub species centrality
│   ├── morans_i_results.csv            # Spatial autocorrelation
│   ├── pipeline_timing.csv             # Script execution times
│   ├── taxonomy_sensitivity.csv        # Sensitivity results (5 scenarios)
│   ├── taxonomy_sensitivity_species_scaling.csv  # Species-level sensitivity
│   └── ...                             # ~36 more analysis tables
├── objects/                 # 22 RDS files
│   ├── coral_master.rds                # Main merged dataset
│   ├── cafi_clean.rds                  # Clean CAFI records
│   ├── community_matrix.rds            # Coral × species matrix
│   ├── condition_scores.rds            # PCA condition scores
│   ├── otu_taxonomy.rds               # OTU resolution lookup (from 01_load_data.R)
│   ├── taxonomy_scenario_data.rds     # Pre-built sensitivity scenarios (from 01_load_data.R)
│   ├── cafi_network.rds                # Network + modularity results
│   ├── scaling_analysis_results.rds    # Scaling models + bootstrap
│   └── ...                             # 12 more analysis objects
└── pipeline.log             # Execution log
```

**Note**: Archived/exploratory scripts (`scripts/archive/`) are NOT part of the pipeline and do not generate any of the above outputs. Ignore them when building the paper.

## Code Conventions

- **Join key**: `coral_id` links all datasets
- **Site codes**: HAU (Hauru), MAT (Maatea), MRB (Maharepa barrier reef)
- **Volume**: Use `volume` (field estimate)
- **Key columns**: `n_galeropsis` (Galeropsis count per coral), `n_corallivore` (all gastropods)
- **Packages**: Use `dplyr::select()` explicitly (MASS conflict); car::vif(), DHARMa, sandwich/lmtest for diagnostics
- **Colors**: Two semantic palettes, chosen to avoid cross-figure confusion:

  **Site palette** (`SITE_COLORS` in `00_setup.R`; used globally across all figures):
  | Site | Hex | Description |
  |------|-----|-------------|
  | HAU | `#9B7EB8` | Muted purple — Hauru (fringing) |
  | MAT | `#7B9BAE` | Cool slate — Maatea (back-reef) |
  | MRB | `#7AAC6D` | Sage green — Maharepa (barrier reef) |

  **Scaling-class palette** (Figure 2 Panel C):
  | Class | Hex | Usage |
  |-------|-----|-------|
  | Redirection (β < 1) | `#5A8FAF` | Muted blue |
  | Field of Dreams (CI overlaps 1) | `gray55` | Neutral — consistent with null |
  | Super-linear (CI > 1) | `#D55E00` | Vermillion accent — distinctive result |

  Site palette deliberately avoids blue and orange/vermillion to prevent confusion with scaling-class semantics.

## Helper Functions & Object Loading

### Loading Objects
```r
coral_master <- load_object("coral_master")
cafi_clean <- load_object("cafi_clean")
community_matrix <- load_object("community_matrix")
```

### Utility Functions (defined in `00_setup.R`)

| Function | Purpose |
|----------|---------|
| `save_object(obj, name)` | Saves RDS to `output/objects/name.rds` |
| `load_object(name)` | Loads RDS from `output/objects/name.rds` (errors if missing) |
| `save_figure(plot, path, width, height, units, dpi)` | Saves PNG + PDF dual format |
| `save_table(df, name)` | Saves CSV to `output/tables/name.csv` |
| `calc_pseudo_r2(model, null_model)` | McFadden's pseudo-R² for GLMs. Pass null_model explicitly inside functions (scoping). |
| `flag_influential(model, threshold)` | Flags Cook's D > 4/n |
| `theme_publication(base_size)` | Publication-ready ggplot theme |
| `theme_multipanel(base_size)` | Compact theme for dense multi-panel figures |

### Global Color Palettes (defined in `00_setup.R`)

| Variable | Colors |
|----------|--------|
| `SITE_COLORS` | HAU=#9B7EB8, MAT=#7B9BAE, MRB=#7AAC6D |
| `SCALING_COLORS` | Redirection=#5A8FAF, FoD=gray55, SuperLinear=#D55E00 |
| `FUNC_GROUP_COLORS` | Trapezia, Fish, Gastropod, Shrimp, Other Crab, Other |
| `TYPE_COLORS` | Broad taxonomic groups |

## Script Dependency Matrix

Scripts 00 and 01 must always run first. After that, most scripts can run independently.

| Script | Loads Objects From | Creates Objects | Needs Prior Script? |
|--------|--------------------|-----------------|---------------------|
| `00_setup.R` | — | PATHS, themes, helpers | Always first |
| `01_load_data.R` | PATHS | coral_master, community_matrix, condition_scores, cafi_pca_results, otu_taxonomy, taxonomy_scenario_data | Always second |
| `02_community_analysis.R` | coral_master, community_matrix | community_analysis_results | No |
| `03_landscape_characterization.R` | coral_master | landscape_selected_predictors | No |
| `04_landscape_effects.R` | coral_master, community_matrix | — | No |
| `05_species_scaling_analysis.R` | coral_master | scaling_analysis_results | No |
| `06_cooccurrence_analysis.R` | coral_master, community_matrix | cooccurrence_results, cafi_network | No |
| `07_spatial_autocorrelation.R` | coral_master, landscape_selected_predictors | — | Yes: needs 03 |
| `08_functional_groups.R` | coral_master, scaling_analysis_results | — | Yes: needs 05 |
| `09_cafi_condition_feedbacks.R` | coral_master, condition_scores, cafi_pca_results | — | No |
| `13_taxonomy_sensitivity.R` | coral_master, taxonomy_scenario_data | — | No |

## Data Schema Quick Reference

### survey_cafi_data_w_taxonomy_summer2019_v5.csv (3,989 rows)
Key columns: `coral_id`, `genus`, `species`, `type`, `cafi_size_mm`, `family`
Note: `site` column unreliable — always extract site from `coral_id` prefix (first 3 chars: HAU/MAT/MRB)

### survey_coral_characteristics_merged_v2.csv (114 rows)
Key columns: `coral_id`, `site`, `volume`, `depth`, `branch_width`, `lat`, `long`, `number_of_neighbors`
Critical: Two survey types — `neighborhood` (n=61) has neighbor data, `size` (n=53) has NA for neighbor columns
Join key: `coral_id`

### survey_master_phys_data_v3.csv (108 rows)
Key columns: `coral_id`, `nub`, `protein_mg_cm2`, `carb_mg_cm2`, `zooxDensity`, `afdw_mg_cm2`, `stump_length`, `nubbin_length`
Note: Position-corrected values created in script 01 by regressing traits on stump_length + nubbin_length

### Join Logic
- All joins on `coral_id` (format: "SITE-POC##", e.g., "HAU-POC29")
- CAFI → Coral: aggregate CAFI by coral_id first, then left_join to coral characteristics
- Physiology subset: only 108/114 corals have physiology; condition PCA uses n=84 (complete cases)

## Known Issues & Troubleshooting

| Error | Cause | Solution |
|-------|-------|----------|
| `Error in dplyr::select()` / MASS conflict | MASS loaded before dplyr | Already handled in 00_setup.R via `conflicted::conflict_prefer("select", "dplyr")` |
| `object 'PATHS' not found` | Setup not sourced | Always run `source("scripts/00_setup.R")` first |
| `object 'coral_master' not found` | Data not loaded | Always run `source("scripts/01_load_data.R")` second |
| `NB convergence failure` (script 04) | Negative binomial GLM won't converge | Automatic fallback to Poisson with console warning |
| `Cairo PDF failed` | Cairo device unavailable | `save_figure()` falls back to standard PDF device |
| PERMANOVA p-values vary between runs | Uses permutation randomness | Expected; set `set.seed()` before script 02 for exact reproducibility |

### Data Quirks
- **CAFI site column unreliable**: Extract site from `coral_id` prefix, not the `site` column
- **Neighborhood data NA for 53 corals**: survey_type="size" corals have no neighborhood surveys — not "isolated," just not surveyed
- **Morphotype is putative**: Not genetically confirmed; some cryptic species complexes may be present

## Manuscript Editing Workflow

### Editing Section Text

1. Edit the section file: `manuscript/introduction.md`, `methods.md`, `results.md`, or `discussion.md`
2. Manually update `combined_manuscript.md`: replace the corresponding section between headers
3. Preserve all other sections (Abstract, Keywords, Figure Legends, References)
4. Verify inline figure paths: `![](../output/figures/manuscript/figN_*.png)`
5. Commit both files together to keep them in sync

### Regenerating Figures

1. Run the source script: `run_one("05")` (e.g., to regenerate Fig 2-3)
2. Figures auto-update in `output/figures/manuscript/` (PNG + PDF)
3. Legend files also regenerate: `figN_legend_results.txt`
4. Inline image paths in `combined_manuscript.md` don't change — they auto-reference updated PNGs

### Pending Before Submission
- 1 BibTeX placeholder: StierInReview (companion experimental paper — update when published)
- Repository URL/DOI placeholders in Methods Data Accessibility section

---

## Manuscript

### Combined Manuscript

The full assembled manuscript lives in `manuscript/combined_manuscript.md`. This is the **single authoritative document** for submission preparation. It contains all sections in JAE order:

1. Title page (with placeholder co-authors)
2. Abstract (~230 words)
3. Keywords
4. Introduction
5. Methods
6. Results
7. Discussion
8. Acknowledgements (placeholder)
9. Author Contributions (placeholder)
10. Conflict of Interest
11. Figure Legends (all 5 main figures, with inline figure images)
12. References (→ `manuscript/references.bib`, 82 BibTeX entries)
13. Supplement pointer (→ `manuscript/combined_supplement.md`)

**Figures are embedded inline** at the point of first reference in the text (markdown `![](...)` syntax with relative paths to `output/figures/manuscript/`). Figure legends also appear with images in the Legends section.

**One placeholder remains**: `[cite MCR-LTER / Edmunds / emerging work]` in the Discussion heatwave paragraph.

### Section Source Files

Individual section files remain in `manuscript/` for editing convenience. **After editing any section file, re-assemble into `combined_manuscript.md`.**

| File | Status |
|------|--------|
| `manuscript/combined_manuscript.md` | **Authoritative** — full assembly with inline figures |
| `manuscript/combined_supplement.md` | Complete supplement: 14 figures (inline) + 11 tables + methods |
| `manuscript/references.bib` | BibTeX file (81 entries, 70 with DOIs; 1 placeholder: StierInReview) |
| `manuscript/figure_index.md` | Standalone figure inventory (all 19 figures with panels + descriptions) |
| `manuscript/introduction.md` | Section source file for editing |
| `manuscript/methods.md` | Section source file for editing |
| `manuscript/results.md` | Section source file for editing |
| `manuscript/discussion.md` | Section source file for editing |
| `manuscript/nlm_literature_notes.md` | Literature scaffold (8 thematic queries from NotebookLM) |
| `manuscript/nlm_paper_index.md` | 146 source summaries (auto-generated reference) |
| `manuscript/key citations/` | 3 PDFs (Stier & Osenberg 2010, 2024a, 2024b) + experimental paper docx files |

### Figure Legend Files

Legend files in `output/figures/manuscript/` contain **publication-ready captions only** — no embedded methods, results, statistics, or color scheme sections. Format: `figN_legend_results.txt`.

### Figure Panel Order (verified from scripts)

| Figure | Panel A | Panel B | Panel C | Panel D |
|--------|---------|---------|---------|---------|
| Fig 2 | Abundance scaling | Density dilution | Species richness | — |
| Fig 3 | Species curves | Group curves | Species β forest | Group β forest |
| Fig 4 | db-RDA biplot | Taxonomic barchart | — | — |
| Fig 5 | Richness scatter | Abundance scatter | — | — |

## File Ownership (parallel work)
- `scripts/` — analysis scripts (00-13); can be edited independently
- `scripts/archive/` — archived/exploratory scripts; NOT part of manuscript or pipeline
- `manuscript/` — text files, independent from analysis code
- `data/` — READ ONLY (raw data; archived files in `data/archive/`)
- `output/` — generated files, safe to regenerate
