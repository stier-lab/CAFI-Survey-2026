# CLAUDE.md - AI Assistant Context

This file provides essential context for AI assistants working with this codebase.

## Project Overview

**CAFI Survey Analysis** - A marine ecology research project studying how coral size and neighborhood context shape coral-associated fauna (CAFI) communities in Mo'orea, French Polynesia. **Target journal: Journal of Animal Ecology** — write for a broad ecological audience (community ecology, BEF, landscape ecology, conservation), using the coral system as a case study to illustrate general principles.

- **Survey scope**: 114 *Pocillopora* coral colonies across 3 reef sites
- **CAFI catalogued**: ~4,000 individual specimens spanning 243 OTUs (154 species-level, 89 genus/family/higher)
- **Research focus**: How coral size structures CAFI abundance, richness, and composition (scaling + community assembly)

## Companion Repos (CAFI Trilogy + Supporting)

This project is **Paper 3** of a three-paper series. All use the same CAFI system in Mo'orea but manipulate different habitat axes:

| Paper | Repo (local) | GitHub | Lever | Status |
|-------|-------------|--------|-------|--------|
| **1. CAFI136 / MRB Amount** (experiment) | `~/Stier-2025-CAFI136-MRB-AMOUNT` | `adrianstier/Stier-2025-CAFI136-MRB-AMOUNT` + `adrianstier/coral-cafi-density-experiment` | Coral density (1, 3, 6 per 25m²) | Publication-ready, Zenodo archived (DOI: 10.5281/zenodo.18239647) |
| **2. Maatea Size** (observational) | `~/Stier-CAFI-Size-Maatea-2026` | `stier-lab/Stier-CAFI-Size-Maatea-2026` | Coral size + branch spacing (60 colonies, 39 survivors) | Analysis pipeline built, manuscript drafting |
| **3. CAFI Survey 2026** (this repo) | `~/CAFI-Survey-2026` | `stier-lab/CAFI-Survey-2026` | Coral size + neighborhood (natural variation, 114 colonies) | Manuscript drafting |

### Supporting Repos

| Repo | Local | GitHub | Purpose | Status |
|------|-------|--------|---------|--------|
| **CAFI_2025** (monorepo) | `~/CAFI_2025` | `stier-lab/CAFI_2025` | Original monorepo with all 3 study pipelines (MRB + Maatea + Survey) | BCO-DMO data submission stage |
| **moorea-cafi-data** (data archive) | `~/moorea-cafi-data` | — | BCO-DMO compliant data package: 23 data files across all 3 studies | BCO-DMO submission-ready |

**Cross-repo notes:**
- `CAFI_2025` is the original monorepo containing all three study pipelines. This repo (`CAFI-Survey-2026`) extracted and refined the Survey component; `Stier-CAFI-Size-Maatea-2026` extracted and refined the Maatea component.
- `Stier-2025-CAFI136-MRB-AMOUNT` is the archival repo for Paper 1 (the density experiment). The manuscript references it as `StierInReview` in `references.bib`.
- `Stier-CAFI-Size-Maatea-2026` is Paper 2 — examines how coral size and interbranch distance drive CAFI community assembly at the Maatea reef site (60 *Pocillopora* colonies, with growth tracked to May 2021).
- `moorea-cafi-data` is the BCO-DMO data submission package for the NSF final report — contains cleaned, metadata-documented versions of all raw data files across all three studies.
- Shared data heritage: raw survey data in `data/` here originated from `CAFI_2025/data/Survey/`.
- Cross-study comparisons (species scaling concordance, condition sign concordance) reference Paper 1 results directly.

## Relationship to Companion Studies

The three papers form a complementary trilogy testing how habitat structure shapes CAFI communities using different approaches:

| Paper | Approach | Lever | Strength |
|-------|----------|-------|----------|
| **1. MRB Amount** | Experiment | Coral density (manipulated) | Causal inference, colonization dynamics |
| **2. Maatea Size** | Observational | Coral size + branch morphology (natural) | Morphological mechanisms, growth feedbacks |
| **3. Survey** (this repo) | Observational | Coral size + neighborhood (natural) | Ecological realism, landscape context, established communities |

## Three Core Research Questions (Q1-Q3)

| Question | Script | Key Result |
|----------|--------|------------|
| **Q1: SCALING + NEIGHBORHOOD** | `05_species_scaling_analysis.R`, `04_landscape_effects.R` | Abundance β=0.52 (sublinear, Redirection); Richness z=0.34 (sublinear); density dilution; neighbor distance predicts richness (p=0.001) and rarefied richness (p=0.005) |
| **Q2: COMPOSITION** | `02_community_analysis.R` | Site pools + size gradient structure composition; db-RDA confirms size effect survives rarefaction; coral species explains 8.3% of composition independently (p=0.001); 28/52 OTUs show genotype effects (13 genuine, 15 architecture-mediated) |
| **Q3: FEEDBACKS** | `09_cafi_condition_feedbacks.R` | Species richness predicts condition (BEF framework, p=0.018); variance partitioning: 29% unique to richness, <1% to abundance |

### Supporting Analyses (Supplement)

| Analysis | Script | Key Result |
|----------|--------|------------|
| **Co-occurrence** | `06_cooccurrence_analysis.R` | Volume-weighted null model; pairwise associations mostly explained by volume + site (Fig S9) |
| **Neighborhood (count/volume)** | `04_landscape_effects.R` | n_neighbors and total_neighbor_volume NOT significant (all p>0.37); distance effect in main text under Q1 |
| **Genotype × CAFI** | `16_genotype_cafi_analysis.R` | Coral species predicts CAFI composition (R²=0.083, p=0.001); 28/52 OTUs significant (13 genuine, 15 architecture-mediated); Galeropsis×P.verrucosa, T.serenei character displacement |
| **Community assembly** | `15_community_assembly.R` | Raup-Crick null model (deterministic vs stochastic); Mantel/partial Mantel (dispersal limitation); 3-way variation partitioning (size vs architecture vs space); beta-dispersion convergence; SAD fitting; NRI/NTI taxonomic structure; beta diversity decomposition |

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
- Host species identity (genetic): coral species explains 8.3% of composition (PERMANOVA p=0.001); 28/52 OTUs show genotype effects (13 genuine, 15 architecture-mediated); three-way variance partitioning: architecture 5.6% > volume 4.6% > site 2.9% unique
- **Main text Q2 now has 3 summary sentences** pointing to supplement for genotype/OTU/architecture details (Tables S16–S19, Figs S12–S14)
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

### Neighborhood context (reported under Q1 in main text)

**Hypotheses tested:**
- **Proximity enhancement**: Corals closer to neighbors harbour more species (shared propagule pool)
- **Density recruitment**: More neighbors → more CAFI (spillover/source-sink)

**Methods:**
- Full model: `response ~ log(volume) + n_neighbors + log(total_neighbor_volume+1) + mean_neighbor_dist + site`
- AIC backward elimination; rarefied richness vs distance; functional group distance effects
- Available on 61/114 corals (5m survey subset); rarefied richness available for 39/61

**Key findings:**
- **mean_neighbor_dist → richness: β = −0.005, z = −3.20, p = 0.001 (SIGNIFICANT)**
- **mean_neighbor_dist → Shannon: β = −0.007, t = −3.49, p = 0.001 (SIGNIFICANT)**
- **Rarefied richness → distance: β = −0.041, t = −2.97, p = 0.005 (survives rarefaction)**
- mean_neighbor_dist → abundance: p = 0.78 (NS)
- n_neighbors NOT significant for any response (all p > 0.37)
- total_neighbor_volume NOT significant for any response (all p > 0.79)
- No size × neighborhood interactions (all p_FDR > 0.22)
- No functional group shows significant distance response (all p > 0.44)
- Compositional variability: betadisper trend p = 0.53 (no continuous relationship)
- **Interpretation: proximity matters, raw count does not; genuine species accumulation (not passive sampling)**

---

## Key Data Files

| File | Description | Records |
|------|-------------|---------|
| `data/survey_cafi_data_w_taxonomy_summer2019_v5.csv` | CAFI specimen records with taxonomy | 3,989 rows |
| `data/survey_coral_characteristics_merged_v2.csv` | Coral colony attributes (size, position, neighbors) | 114 rows |
| `data/survey_master_phys_data_v3.csv` | Coral physiology measurements | 108 rows |
| `data/survey_coral_haplotypes_v1.csv` | Coral mtORF haplotype assignments (species ID) | 114 rows |
| `data/genetic_references/Johnston2022_Pocillopora_tree.nex` | Published Pocillopora phylogeny (7 tips) | — |
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
| `16_genotype_cafi_analysis.R` | Host species PERMANOVA, OTU × genotype screen, indicator species, Trapezia body size | Fig S12–S14, Tables S16–S19 |

#### Q3: FEEDBACKS (Fig 5)
| Script | Analysis | Output |
|--------|----------|--------|
| `09_cafi_condition_feedbacks.R` | PCA, fixed-effect LMs, FDR correction | Fig 5 (BEF richness + abundance scatter) |

#### SUPPORTING (Supplement)
| Script | Analysis | Output |
|--------|----------|--------|
| `06_cooccurrence_analysis.R` | Null-model pairwise co-occurrence, intraspecific density | Fig S9 (co-occurrence) |
| `04_landscape_effects.R` | GLMs: size + neighbors → abundance/diversity | Fig S5 (neighborhood null) |
| `15_community_assembly.R` | Q2 (assembly): Raup-Crick null model, dispersal limitation, variation partitioning, beta-dispersion, NRI/NTI | Fig S15–S17, Table S20 |

#### SENSITIVITY (cross-cuts Q1-Q4)
| Script | Analysis | Output |
|--------|----------|--------|
| `13_taxonomy_sensitivity.R` | 5 taxonomy scenarios × 7 metrics; uses pre-built data from 01_load_data.R | Fig S6 + sensitivity tables |
| `14_supplementary_sensitivity.R` | 11+ sensitivity analyses (4 in supplement: betapart, mediation, morphotype, missing data; Parts 8b-8g: morphotype-haplotype concordance, BEF haplotype covariate, scaling by species, composition by species, genotype-architecture mediation (branch width × species), phylogenetic distance & symbiont identity; 7 archived: envfit, IndVal, iNEXT, CWM, nonlinear BEF, C-score, MuMIn) | Tables S11–S15 |

### Archived Scripts (`scripts/archive/` — NOT part of the manuscript or pipeline)

These scripts are retained for reference only. Do not run, edit, or cite them when building the paper.

| Script | Purpose |
|--------|---------|
| `archive/06_network_analysis_legacy.R` | Legacy network analysis (replaced by `06_cooccurrence_analysis.R`) |
| `archive/10_feature_engineering.R` | Exploratory ML feature creation, VIF selection |
| `archive/11_machine_learning.R` | Exploratory Random Forest, XGBoost models |
| `archive/12_model_evaluation.R` | Exploratory cross-validation, diagnostics |
| `archive/15_reviewer_response_analyses.R` | Pre-built code for reviewer requests: mvabund, GDM, TITAN2, Bayesian SEM, LOOCV, Bayesian multilevel scaling, cross-study meta-analysis. Set `if(FALSE)` blocks to `TRUE` to run. |

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
| **S7** | `08_functional_groups.R` | *Archived* — taxonomic group scaling (redundant with Fig 3D) |
| **S8** | `09_cafi_condition_feedbacks.R` | Rarefaction depth sensitivity for richness → condition |
| **S9** | `06_cooccurrence_analysis.R` | Co-occurrence: pairwise SES heatmap + intraspecific density + size-dependent (3-panel) |
| **S10** | `09_cafi_condition_feedbacks.R` | BEF variance partitioning + partial regression + path model coefficients |
| **S11** | `09_cafi_condition_feedbacks.R` | A priori forest + rarefied richness + bidirectional (trimmed from 6 to 3 panels) |
| **S12** | `16_genotype_cafi_analysis.R` | Genotype composition: db-RDA biplot + variance partitioning |
| **S13** | `16_genotype_cafi_analysis.R` | Genotype screen forest plot (mediation classification) |
| **S14** | `16_genotype_cafi_analysis.R` | Trapezia body size: architecture filtering + character displacement |
| **S15** | `15_community_assembly.R` | Beta-dispersion convergence across coral size classes |
| **S16** | `15_community_assembly.R` | Community assembly null-model: Raup-Crick + beta-dispersion + variation partitioning (4-panel) |
| **S17** | `15_community_assembly.R` | Taxonomic structure (SES.MPD / NRI) vs coral volume |

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
| Beta diversity partitioning | betapart: 81% turnover, 19% nestedness; turnover decreases with volume | `14_supplementary_sensitivity.R` |
| Envfit species vectors | Formal permutation-tested R² for species-NMDS associations | `14_supplementary_sensitivity.R` |
| Indicator species by size class | IndVal for small/medium/large corals (41 significant) | `14_supplementary_sensitivity.R` |
| Coverage-based rarefaction (iNEXT) | Hill q=0,1,2 by size class; per-coral scaling exponents | `14_supplementary_sensitivity.R` |
| CWM obligate-mutualist fraction | Functional identity does not predict condition (p=0.54) | `14_supplementary_sensitivity.R` |
| Nonlinear BEF sensitivity | Linear best (dAIC=0 vs log dAIC=3.9, poly dAIC=1.0) | `14_supplementary_sensitivity.R` |
| Missing data characterization | Only site (MRB) predicts dropout; richness/volume NS | `14_supplementary_sensitivity.R` |
| Morphotype confound | Adding morphotype absorbs richness signal (p goes from 0.015→0.27) | `14_supplementary_sensitivity.R` |
| Mediation (richness→abundance→condition) | ACME NS (p=0.69); richness effect is direct, not abundance-mediated | `14_supplementary_sensitivity.R` |
| C-score community co-occurrence | SES=2.60, p=0.017; weak community-wide segregation under quasiswap null | `14_supplementary_sensitivity.R` |
| Morphotype-haplotype concordance | Chi-square/Fisher's exact, Cramer's V | `14_supplementary_sensitivity.R` |
| BEF haplotype covariate | Richness survives adding genetic species (p=0.010 vs 0.008) | `14_supplementary_sensitivity.R` |
| Scaling by species | Volume x species interaction NS (p=0.39); beta consistent across species | `14_supplementary_sensitivity.R` |
| Composition by species | Marginal PERMANOVA: species R²=0.098, p=0.001 | `14_supplementary_sensitivity.R` |
| Phylogenetic distance | Mantel r=0.12, p=0.013; partial Mantel r=0.13, p=0.002 | `14_supplementary_sensitivity.R` |
| Symbiont identity | PERMANOVA R²=0.038, p=0.001; symbiont->condition beta=1.28, p=0.008 | `14_supplementary_sensitivity.R` |
| Host species confound | Species-level PERMANOVA, variance partitioning, architecture mediation test | `16_genotype_cafi_analysis.R` |
| Architecture vs genotype | Mediation classification: 13 genuine vs 15 architecture-mediated OTUs | `16_genotype_cafi_analysis.R` |
| Genotype-architecture mediation | Branch width × species Fisher's test; species→richness/abundance/condition/composition controlling volume+site | `14_supplementary_sensitivity.R` |
| Body size partitioning | Trapezia species × architecture chi-square, character displacement test | `16_genotype_cafi_analysis.R` |

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
│   ├── supplement/          # Supplementary figures (figS1-S17, S7 archived)
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
│   ├── genotype_permanova_marginal.csv  # Host species PERMANOVA results
│   ├── genotype_species_screen.csv      # OTU × genotype screen with mediation
│   ├── genotype_indicator_species.csv   # Indicator species for coral hosts
│   ├── genotype_host_specificity.csv    # Host specificity index
│   ├── trapezia_body_size_genotype.csv  # Trapezia body size by host
│   ├── trapezia_pair_composition.csv    # Trapezia pair assemblages
│   ├── sensitivity_genotype_architecture.csv # Genotype→architecture mediation summary
│   ├── sensitivity_phylogenetic_symbiont.csv # Phylogenetic distance + symbiont identity
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
- **Key columns**: `n_galeropsis` (Galeropsis count per coral), `n_corallivore` (all gastropods), `poc_species` (genetic species: P. grandis, P. meandrina, P. tuahiniensis, P. verrucosa), `haplotype` (mtORF haplotype code: 1a_Pe, 1a_Pm, 10, 3a, 3b, 3f, 3h, 8a), `haplotype_valid` (TRUE for 101 successfully genotyped colonies)
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
| `14_supplementary_sensitivity.R` | coral_master, community_matrix, condition_scores, cafi_clean | — | No |
| `15_community_assembly.R` | coral_master, community_matrix | — | No |
| `16_genotype_cafi_analysis.R` | coral_master, community_matrix, cafi_clean | — | No (uses haplotype data from coral_master) |

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
Note: Position-corrected values created in script 01 by regressing traits on stump_length + nubbin_length. Also contains a `haplotype` column (mtORF), but only for 108/114 corals — use `survey_coral_haplotypes_v1.csv` for complete coverage including `poc_species` and `haplotype_valid`.

### survey_coral_haplotypes_v1.csv (114 rows)
Key columns: `coral_id`, `site`, `haplotype`, `poc_species`, `haplotype_valid`
Coral mtORF haplotype assignments for all 114 survey corals. 101 successfully genotyped (`haplotype_valid=TRUE`), 11 did not amplify, 2 no sample. Haplotype-to-species mapping: 1a_Pe→P. grandis (n=49), 1a_Pm→P. meandrina (n=34), 3a/3b/3f/3h→P. verrucosa (n=10), 10→P. tuahiniensis (n=7), 8a→P. acuta (n=1). Provenance: extracted from CAFI_2025 physiology master + original genotype file. Full metadata in `data/README_survey_haplotype_metadata_v1.md`.

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
- **Morphotype is putative but now genetically confirmed**: mtORF haplotyping resolved 101/114 colonies to 4 species (P. grandis n=49, P. meandrina n=34, P. tuahiniensis n=7, P. verrucosa n=10). Field morphotype shows 33-94% concordance with genetic species (Table S15). Haplotype data loaded in script 01.

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
12. References (→ `manuscript/references.bib`, 64 BibTeX entries)
13. Supplement pointer (→ `manuscript/combined_supplement.md`)

**Figures are embedded inline** at the point of first reference in the text (markdown `![](...)` syntax with relative paths to `output/figures/manuscript/`). Figure legends also appear with images in the Legends section.

**One placeholder remains**: `[cite MCR-LTER / Edmunds / emerging work]` in the Discussion heatwave paragraph.

### Section Source Files

Individual section files remain in `manuscript/` for editing convenience. **After editing any section file, re-assemble into `combined_manuscript.md`.**

| File | Status |
|------|--------|
| `manuscript/combined_manuscript.md` | **Authoritative** — full assembly with inline figures |
| `manuscript/combined_supplement.md` | Complete supplement: 17 figures (S1–S17, S7 archived; inline) + 20 tables (S1–S20) + methods |
| `manuscript/references.bib` | BibTeX file (64 entries; 1 placeholder: StierInReview; includes Johnston, Cunning & Burgess 2022 and Burgess et al. 2021 for haplotype/symbiont analyses) |
| `manuscript/figure_index.md` | Standalone figure inventory (all 22 figures with panels + descriptions; S12–S14 entries may be stale — combined_supplement.md is authoritative) |
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
- `scripts/` — analysis scripts (00–16); can be edited independently
- `scripts/archive/` — archived/exploratory scripts; NOT part of manuscript or pipeline
- `manuscript/` — text files, independent from analysis code
- `data/` — READ ONLY (raw data; archived files in `data/archive/`)
- `output/` — generated files, safe to regenerate
