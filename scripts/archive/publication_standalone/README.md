# CAFI Publication Analysis Scripts

This folder contains the refactored analysis pipeline for the CAFI (Coral-Associated Fauna Inhabitants) Survey Paper. The scripts are organized to match the **6-figure publication plan** with clear hypotheses and statistical approaches.

## Quick Start

```bash
# Run entire pipeline
Rscript scripts/publication/run_all_publication_analyses.R

# Or source from R console
source("scripts/publication/run_all_publication_analyses.R")
```

## Script Order & Dependencies

| Script | Purpose | Outputs |
|--------|---------|---------|
| `00_setup.R` | Load packages, set paths, define theme | Environment setup |
| `01_load_data.R` | Load and clean all data | `cafi_clean.rds`, `coral_master.rds`, `condition_scores.rds`, etc. |
| `02_fig2_landscape_scaling.R` | Size/proximity -> CAFI abundance/richness/diversity | Figure 2 panels |
| `03_fig3_functional_groups.R` | Functional group-specific responses | Figure 3 panels |
| `04_fig4_composition_networks.R` | NMDS, PERMANOVA, co-occurrence networks | Figure 4 panels |
| `05_fig5_condition_landscape.R` | Coral condition vs landscape (baseline) | Figure 5 panels |
| `06_fig6_cafi_feedbacks.R` | CAFI -> coral condition (feedbacks) | Figure 6 panels |
| `run_all_publication_analyses.R` | Master orchestrator | All figures + summary |

## Figure Plan Overview

### Figure 1: Study System + CAFI Overview (Conceptual)
- **Purpose**: Introduce study system, key taxa, conceptual framework
- **Stats**: None (conceptual figure)
- **Key elements**: Mo'orea map, Trapezia/fish/Drupella photos, scaling diagrams

### Figure 2: Landscape -> CAFI Scaling
- **Hypotheses**:
  - H1: Larger corals host more CAFI (positive size-abundance slope)
  - H2: Isolated corals host more CAFI per unit size (propagule redirection)
  - H3: Diversity increases with size but may decrease with proximity
- **Models**: NB GLM for abundance; Poisson/NB for richness; LM for Shannon
- **Key test**: Power-law exponent (expect 0.75-0.85 for sublinear scaling)

### Figure 3: Functional & Taxonomic Group Responses
- **Functional groups**:
  - Trapezia (Mutualist Defenders)
  - Resident Fishes (Nutrient Providers)
  - Corallivores (Drupella - negative effects)
  - Other Crabs, Shrimp
- **Analyses**: Group-level GLMs, species-specific slopes, heatmaps

### Figure 4: Compositional Changes + Co-occurrence Networks
- **Methods**: NMDS ordination, PERMANOVA (size, proximity, site effects)
- **Networks**: Spearman correlations, Louvain modularity detection
- **Outputs**: Ordination plots, incidence heatmaps, network graphs

### Figure 5: Coral Condition vs Landscape
- **Purpose**: Establish baseline landscape patterns BEFORE testing CAFI effects
- **Predictors**: Coral size, neighbor proximity, neighbor count
- **Note**: Position-corrected condition scores (sampling bias removed)

### Figure 6: CAFI -> Coral Condition Feedbacks
- **Key hypotheses**:
  - H1: More Trapezia -> higher condition (mutualism)
  - H2: More resident fish -> improved condition (nutrient provisioning)
  - H3: More corallivores -> lower condition
  - H4: Composition effects stronger than abundance alone
- **Models**: LM for condition ~ CAFI metrics + functional groups

## Statistical Notes

### Plotting Rule
**If relationship is nonsignificant (p > 0.10), show points only without regression line.** This prevents misleading visual interpretations.

### Model Selection
- **Abundance**: Negative binomial GLM (overdispersed counts)
- **Richness**: Poisson GLM, switching to NB if overdispersed
- **Diversity (Shannon)**: Linear model
- **Condition**: Linear model (PC1 is continuous)

### Position Correction
All physiological metrics are position-corrected:
1. Regress trait ~ stump_length (sampling position)
2. Extract residuals
3. Standardize to z-scores
4. Use PC1 of corrected traits as "condition score"

This removes the confound where larger colonies were systematically sampled higher on branches.

## Output Structure

```
output/
├── figures/
│   └── publication/
│       ├── fig1_study_system/
│       ├── fig2_landscape_scaling/
│       ├── fig3_functional_groups/
│       ├── fig4_composition_networks/
│       ├── fig5_condition_landscape/
│       ├── fig6_cafi_feedbacks/
│       └── supplementary/
├── tables/
│   ├── fig2_model_results.csv
│   ├── fig3_functional_group_effects.csv
│   ├── fig4_permanova_results.csv
│   ├── fig5_landscape_condition_results.csv
│   ├── fig6_cafi_condition_results.csv
│   └── ...
└── objects/
    ├── cafi_clean.rds
    ├── coral_master.rds
    ├── condition_scores.rds
    ├── fig2_models.rds
    └── ...
```

## Key Data Objects

| Object | Description |
|--------|-------------|
| `cafi_clean.rds` | Individual CAFI records with taxonomy and functional groups |
| `coral_master.rds` | Coral characteristics merged with CAFI metrics |
| `condition_scores.rds` | Position-corrected condition scores (PC1) |
| `community_matrix.rds` | Coral x OTU abundance matrix |
| `functional_summary.rds` | Functional group summary statistics |

## Troubleshooting

### Missing packages
```r
install.packages(c("tidyverse", "vegan", "lme4", "MASS", "igraph", "patchwork"))
```

### Data files not found
Ensure these files exist in `data/`:
- `survey_cafi_data_w_taxonomy_summer2019_v5.csv`
- `survey_coral_characteristics_merged_v2.csv`
- `survey_master_phys_data_v3.csv`

### Script errors
Check the console output for specific error messages. Common issues:
- Missing columns in data files
- Insufficient sample sizes for some analyses
- Package conflicts (restart R session)

## Citation

If using this pipeline, please cite:
- Stier Lab, UC Santa Barbara
- Mo'orea LTER (data collection)
