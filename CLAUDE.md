# CLAUDE.md - AI Assistant Context

This file provides essential context for AI assistants working with this codebase.

## Project Overview

**CAFI Survey Analysis** - A marine ecology research project studying how coral size and neighborhood context shape coral-associated fauna (CAFI) communities in Mo'orea, French Polynesia.

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
| **Q3: FEEDBACKS** | `09_cafi_condition_feedbacks.R` | **No CAFI metric predicts condition**; raw richness signal is abundance artifact (rarefied p=0.45) |

### Supporting Analyses (Supplement)

| Analysis | Script | Key Result |
|----------|--------|------------|
| **Co-occurrence** | `06_cooccurrence_analysis.R` | Volume-weighted null model; pairwise associations mostly explained by volume + site (Fig S) |
| **Neighborhood** | `04_landscape_effects.R` | n_neighbors NOT significant (underpowered, n=63); supplement only |

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
- Species-level scaling for 10 most prevalent CAFI taxa

**Key findings:**
- Total CAFI abundance: β = 0.52 [0.44, 0.61] — **sublinear (Propagule Redirection)**
- Species richness (SAR): z = 0.34 [0.28, 0.41] — **sublinear (Redirection)**
- Rarefied richness (n=20): slope=0.14, p=0.64 — **not significant (passive sampling)**; raw SAR is abundance artifact
- Species occurrence curves: 14/24 prevalent species have significant size-dependent occurrence (FDR<0.05)
- Density dilution: per-capita CAFI density decreases with size (slope = -0.48)
- 7/10 top species: Redirection (β < 1); 3/10: Field of Dreams (CI spans 1); 0/10: super-linear
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
- Strong site effects on composition (PERMANOVA R² ~ 0.04-0.06, p < 0.01)
- db-RDA: volume explains 7.8% of composition (F=9.74, p=0.001); survives rarefaction (2.6%, p=0.001); *T. punctimanus* loads most strongly toward larger corals
- Size-divergence (categorical) NOT significant after rarefaction (p=0.61; Fig. S5) — abundance artifact
- Nestedness (NODF): 18.37, z=−1.09, p=0.277 — **not nested**; small-coral faunas are not subsets of large-coral faunas
- Pairwise co-occurrence reveals significant positive and negative associations after controlling for coral volume
- Intraspecific density patterns tested for mating-pair hypothesis (Stier et al. 2012)
- Conclusion: site pools and continuous size gradient shape composition; co-occurrence is largely explained by volume + site

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
- Key species models: condition ~ species abundance + log(volume) + site (FDR-corrected, 10 species from experimental paper)
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
- Key species (10 tested, matching experimental paper): directions mostly match experiment but none survive FDR
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
| `01_load_data.R` | Data loading, cleaning, creates core objects + taxonomy scenario data |
| `run_full_pipeline.R` | Pipeline orchestrator with logging |

### Script-to-Question Mapping

#### Q1: SCALING
| Script | Analysis | Output |
|--------|----------|--------|
| `05_species_scaling_analysis.R` | NB GLM scaling with bootstrap CI | Field of Dreams vs Redirection |

#### Q2: COMPOSITION (Fig 3: NMDS + Composition)
| Script | Analysis | Output |
|--------|----------|--------|
| `02_community_analysis.R` | PERMANOVA, NMDS, betadisper, rarefaction, db-RDA | Fig 3 + site/size effects |
| `03_landscape_characterization.R` | Neighborhood metrics, spatial patterns | Landscape predictor variables |
| `07_spatial_autocorrelation.R` | Moran's I, LISA, Mantel tests | Spatial structure diagnostics |
| `08_functional_groups.R` | Taxonomic group scaling (loads from script 05), composition | Fig S9 + group patterns |

#### Q3: FEEDBACKS (Fig 4)
| Script | Analysis | Output |
|--------|----------|--------|
| `09_cafi_condition_feedbacks.R` | PCA, fixed-effect LMs, FDR correction | Fig 4 (CAFI-condition null results) |

#### SUPPORTING (Supplement only)
| Script | Analysis | Output |
|--------|----------|--------|
| `06_cooccurrence_analysis.R` | Null-model pairwise co-occurrence, intraspecific density | Fig S (co-occurrence) |
| `04_landscape_effects.R` | GLMs: size + neighbors → abundance/diversity | Fig S7 (neighborhood null) |

#### SENSITIVITY (cross-cuts Q1-Q4)
| Script | Analysis | Output |
|--------|----------|--------|
| `13_taxonomy_sensitivity.R` | 5 taxonomy scenarios × 7 metrics; uses pre-built data from 01_load_data.R | Fig S8 + sensitivity tables |

### Exploratory ML (not in default pipeline; in `scripts/exploration/`)
| Script | Purpose |
|--------|---------|
| `exploration/10_feature_engineering.R` | Feature creation, VIF selection |
| `exploration/11_machine_learning.R` | Random Forest, XGBoost models |
| `exploration/12_model_evaluation.R` | Cross-validation, diagnostics |

### Publication Figures (embedded in analysis scripts — no separate figure script)

Each manuscript figure is created by its source analysis script with **dual saves**: once to the analysis output directory, once to `output/figures/manuscript/` or `output/figures/supplement/`.

| Figure | Script | Description |
|--------|--------|-------------|
| **Fig 1** | `01_load_data.R` | Study design: satellite map + volume/neighborhood distributions + schematics (6-panel) |
| **Fig 2** | `05_species_scaling_analysis.R` | Scaling: abundance + richness GLMs + density dilution + species forest plot |
| **Fig 3** | `02_community_analysis.R` | Composition: NMDS ordination + taxonomic barchart by site |
| **Fig 4** | `09_cafi_condition_feedbacks.R` | CAFI-condition feedbacks (forest plot + key null results) |
| **S1** | `02_community_analysis.R` | Species accumulation curves |
| **S2** | `02_community_analysis.R` | PERMANOVA metric sensitivity |
| **S3** | `02_community_analysis.R` | NMDS ordination by site/size |
| **S4** | `07_spatial_autocorrelation.R` | Spatial autocorrelation (Moran's I) |
| **S5** | `02_community_analysis.R` | Composition divergence by size |
| **S6** | `05_species_scaling_analysis.R` | Species-level scaling forest plot |
| **S7** | `04_landscape_effects.R` | Neighborhood null results |
| **S8** | `13_taxonomy_sensitivity.R` | Taxonomy sensitivity forest plot |
| **S9** | `08_functional_groups.R` | Taxonomic group scaling and composition |
| **S10** | `09_cafi_condition_feedbacks.R` | Rarefaction depth sensitivity for richness → condition |
| **S-cooccurrence** | `06_cooccurrence_analysis.R` | Co-occurrence: pairwise SES heatmap + intraspecific density (demoted from main text) |
| **S11** | `06_cooccurrence_analysis.R` | Legacy co-occurrence network (circular layout) |
| **S12** | `09_cafi_condition_feedbacks.R` | Additional CAFI-condition panels |
| **S13** | `13_taxonomy_sensitivity.R` | Network sensitivity to taxonomic resolution |
| **S14** | `05_species_scaling_analysis.R` | Species occurrence probability vs. coral size (logistic GLM, 24 species) |

**Note**: Deprecated scripts in `scripts/archive/`; ML exploration scripts in `scripts/exploration/`.

## Key Commands

### Run Pipeline
```r
source("scripts/run_full_pipeline.R")
run_pipeline()           # Full pipeline ~12 min (core + extended)
run_full_pipeline()      # Everything including ML exploration

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
| Multiple testing (feedbacks) | FDR correction (Benjamini-Hochberg) | `09_cafi_condition_feedbacks.R` |
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
| Model diagnostics (LM) | Shapiro-Wilk normality test; HC3 robust SEs applied unconditionally | `09_cafi_condition_feedbacks.R` |
| Poisson overdispersion | Pearson X²/df check; auto-switch to quasipoisson | `02`, `04` |
| Heteroscedasticity (count predictors) | HC3 sandwich robust standard errors | `09_cafi_condition_feedbacks.R` |
| Mediation bootstrap | 1000 bootstrap iterations via mediation::mediate() | `09_cafi_condition_feedbacks.R` |
| Community matrix (zero-CAFI corals) | Added zero-abundance rows for all corals | `01_load_data.R` |
| Log volume bias | Removed +1 offset (volume > 0 guaranteed) | `01_load_data.R` |
| Tissue sampling artifact (condition) | Regress physio traits on stump_length + nubbin_length | `01_load_data.R` |
| PERMANOVA Type I confound | Marginal (Type III) PERMANOVA | `02_community_analysis.R` |
| Colorblind inaccessibility | Okabe-Ito palette; Fig 2 site colors shifted to avoid scaling-class overlap | All figure scripts |
| Count predictor skew | sqrt() applied to count-based CAFI predictors in condition models | `09_cafi_condition_feedbacks.R` |
| Key species FWER correction | Hochberg (FWER) for species-level tests (matches experiment) | `09_cafi_condition_feedbacks.R` |
| Community transform sensitivity | 3 transforms + filtered PCA tested for PC1→condition | `09_cafi_condition_feedbacks.R` |
| PERMANOVA site-size robustness | Balanced site subsampling (500 iter, min_n per site) | `02_community_analysis.R` |
| Power analysis (feedbacks) | Prospective power for Q3 (n=84); medium effects detectable | `09_cafi_condition_feedbacks.R` |
| Power analysis (neighborhood) | Prospective power for Q4 (n=63); underpowered for small effects | `04_landscape_effects.R` |
| FDR family vs global | Both family-wise and global FDR reported; global as sensitivity | `09_cafi_condition_feedbacks.R` |
| Position correction sensitivity | Raw vs corrected condition scores compared | `09_cafi_condition_feedbacks.R` |
| Bootstrap fallback logging | BCa→percentile CI fallback logged with warning | `05_species_scaling_analysis.R` |

## Output Structure

All outputs are gitignored and regenerated by the pipeline.

```
output/
├── figures/
│   ├── manuscript/          # 4 main text figures + legend files
│   │   ├── fig1_study_design.png       (from 01_load_data.R)
│   │   ├── fig1_legend_results.txt     (from 01_load_data.R)
│   │   ├── fig2_scaling.png            (from 05_species_scaling_analysis.R; 4-panel with density dilution)
│   │   ├── fig2_legend_results.txt     (from 05_species_scaling_analysis.R)
│   │   ├── fig3_composition.png        (from 02_community_analysis.R)
│   │   ├── fig3_legend_results.txt     (from 02_community_analysis.R)
│   │   ├── fig4_feedbacks.png          (from 09_cafi_condition_feedbacks.R)
│   │   └── fig4_legend_results.txt     (from 09_cafi_condition_feedbacks.R)
│   ├── supplement/          # Supplementary figures (figS1-S14 + co-occurrence)
│   ├── 01_data/             # Study design (1 figure)
│   ├── 02_community/        # Community analysis (11 figures)
│   ├── 03_landscape/        # Landscape characterization (3 figures)
│   ├── 04_effects/          # Landscape effects (14 figures + 2 HTML tables)
│   ├── 05_scaling/          # Species-area scaling (7 figures)
│   ├── 06_network/          # Co-occurrence analysis (fig4_cooccurrence.png)
│   ├── feedbacks/           # CAFI-condition feedbacks (8 figures)
│   └── functional_groups/   # Taxonomic group analysis (7 figures)
├── tables/                  # ~54 CSV statistical results
│   ├── scaling_results_all.csv         # All scaling coefficients
│   ├── cafi_condition_models.csv       # Feedback model results
│   ├── key_species_effects.csv         # Key species condition effects (10 species)
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
├── objects/                 # 20 RDS files
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

**Note**: ML exploratory outputs (scripts 10-12) are NOT generated by the default pipeline. They must be run separately and will create additional files in `output/`.

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

## Load Pre-computed Objects
```r
coral_master <- load_object("coral_master")
cafi_clean <- load_object("cafi_clean")
community_matrix <- load_object("community_matrix")
```

---

## File Ownership (parallel work)
- `scripts/` — individual analysis scripts can be edited independently
- `manuscript/` — text files, independent from analysis code
- `docs/` — documentation, independent
- `data/` — READ ONLY (raw data)
- `output/` — generated files, safe to regenerate
