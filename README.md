# CAFI Survey Analysis 2026

**Nonlinear scaling and neighborhood effects structure coral-associated fauna communities**

Analysis pipeline and manuscript for investigating how coral habitat attributes drive variation in coral-associated fauna (CAFI) density and community structure in Mo'orea, French Polynesia.

## Overview

This repository contains the complete analysis pipeline for a manuscript investigating the ecological and physical factors determining CAFI community structure in *Pocillopora* corals. The work tests predictions from propagule redirection theory and marine landscape ecology.

### Theoretical Background

A critical but underappreciated pattern in marine landscape ecology is that occupant *abundance* scales nonlinearly with habitat amount, which causes occupant *density* (per unit habitat) to decrease as habitat increases. This pattern arises through **propagule redirection**—larvae settling to habitat landscapes distribute among available patches based on chemical cue strength, so isolated habitats receive disproportionately more settlers per unit area than clustered habitats.

Three key habitat attributes drive these density patterns:

1. **Habitat Amount**: Areas with more coral → higher abundance but *lower density*
2. **Habitat Size**: Larger corals → more CAFI but at *lower density*
3. **Habitat Proximity**: Isolated corals → *higher* CAFI density

These patterns have profound implications because CAFI affect coral growth and survival—beneficially through predator defense and sediment removal, or detrimentally through tissue consumption. Habitat-driven variation in CAFI density should therefore feed back to alter coral dynamics and landscape patterns.

### Study Summary

- **Location**: Mo'orea, French Polynesia (3 reef sites: HAU, MAT, MRB)
- **Host Coral**: *Pocillopora* spp.
- **Survey Date**: Summer 2019
- **Sample Size**: 114 coral colonies, 12,834 CAFI individuals, 87 species
- **Coordinates**: GPS positions enable meter-scale neighborhood analysis

### Important Analysis Note

**Morphotype Excluded from Analysis**: Initial analyses included coral morphotype (*P. verrucosa* vs *P. meandrina*) as a potential predictor. However, morphotype was found to be a weak and non-significant predictor of CAFI community patterns compared to other factors like coral size, site, and branch architecture. To focus on the most ecologically meaningful predictors and avoid overfitting, morphotype has been excluded from all final analyses and figures. This decision improves model parsimony and focuses attention on the habitat attributes that genuinely structure CAFI communities.

## Research Questions & Hypotheses

### Central Questions

1. How do coral attributes (amount, size, proximity) determine CAFI community composition and density?
2. What are the scaling relationships between coral size and CAFI abundance?
3. How do meter-scale neighborhood effects influence CAFI communities?
4. Do CAFI density patterns relate to coral physiological condition?

### Hypotheses

| Hypothesis | Prediction |
|------------|------------|
| H1 (Site effects) | CAFI community composition differs among reef sites due to variation in coral landscapes and environmental conditions |
| H2 (Size scaling) | CAFI abundance scales with coral volume following a power-law with exponent < 1, indicating larger corals have lower CAFI densities |
| H3 (Neighborhood effects) | CAFI abundance varies with local coral density and proximity, reflecting propagule redirection and spillover effects |
| H4 (Condition relationships) | Coral physiological condition positively predicts CAFI diversity, consistent with bidirectional coral-CAFI interactions |
| H5 (Network structure) | CAFI co-occurrence networks exhibit non-random modular structure with identifiable keystone species |

## Repository Structure

```
CAFI-Survey-2026/
├── agents/                 # Python research workflow agents
│   ├── all_agents.py       # 8 specialized research agents
│   ├── orchestrator_agent.py
│   └── research_prd_agent.py
│
├── data/                   # Raw input data
│   ├── survey_cafi_data_w_taxonomy_summer2019_v5.csv
│   ├── survey_coral_characteristics_merged_v2.csv
│   ├── survey_master_phys_data_v3.csv
│   └── README_*.xlsx       # Metadata documentation
│
├── scripts/                # R analysis pipeline (17 numbered scripts)
│   ├── 00_load_libraries.R
│   ├── 01_load_clean_data.R
│   ├── 02_community_composition.R
│   ├── ...
│   ├── 14_local_neighborhood_effects.R  # Meter-scale neighborhood analysis
│   └── run_all_survey_analyses.R
│
├── output/                 # Generated results
│   ├── figures/            # PNG visualizations (300 dpi)
│   │   ├── neighborhood_effects/
│   │   └── ...
│   ├── tables/             # CSV statistical results
│   ├── objects/            # RDS R data objects
│   ├── manuscript/         # MEPS-formatted manuscript
│   └── reports/            # Documentation and summaries
│
├── docs/                   # Project documentation
│   └── PRD.md              # Product Requirements Document
│
└── README.md
```

## Key Analyses

### 1. Volume-Abundance Scaling

Test whether CAFI abundance scales nonlinearly with coral size:

- **Power-law model**: log(Abundance) ~ log(Volume)
- **Expected exponent**: ~0.75 for 3D habitat (if < 1, density decreases with size)
- **GLMM**: Site random effects, branch architecture fixed effect

### 2. Local Neighborhood Effects

Quantify meter-scale spatial context effects:

- **Neighbor density**: Number of corals within 5m radius
- **Neighbor volume**: Total habitat in local neighborhood
- **Isolation index**: Mean distance to nearest neighbors
- **Relative size**: Focal coral size relative to neighbors
- **Spillover potential**: Colonist supply from nearby corals

### 3. Coral Condition Relationships

Link CAFI patterns to coral physiological metrics:

- **Condition score**: Integrated PC1 of protein, carbohydrate, zooxanthellae, chlorophyll
- **Model**: Diversity ~ Condition + Volume + Branch Width + (1|Site)

### 4. Community Structure

- **PERMANOVA**: Test site, depth, branch architecture effects on composition
- **NMDS ordination**: Visualize community patterns
- **Network analysis**: Co-occurrence patterns and modularity (Q = 0.42)

## Quick Start

### Prerequisites

**R packages** (installed automatically by `00_load_libraries.R`):

- Core: `tidyverse`, `here`, `readxl`, `janitor`
- Statistics: `vegan`, `lme4`, `lmerTest`, `glmmTMB`, `mgcv`
- Spatial: `sf`, `sp`, `geosphere`, `ape`
- Visualization: `ggplot2`, `patchwork`, `viridis`

**Python packages** (for research agents):

```bash
pip install anthropic pyyaml python-dotenv
```

### Running the Analysis

**Full pipeline:**

```r
source("scripts/run_all_survey_analyses.R")
source("scripts/run_all_comprehensive_analyses.R")
```

**Key individual analyses:**

```r
source("scripts/00_load_libraries.R")
source("scripts/01_load_clean_data.R")
source("scripts/05_coral_cafi_relationships.R")  # Scaling relationships
source("scripts/14_local_neighborhood_effects.R") # Neighborhood analysis
```

**Generate manuscript:**

```r
rmarkdown::render("output/manuscript/MANUSCRIPT.Rmd")
```

## Key Findings

### Scaling Relationships (H2 supported)

- Power-law exponent = 0.81 ± 0.12 (95% CI: 0.58–1.04)
- Not significantly different from 0.75 theoretical prediction
- Confirms that CAFI *density* decreases with coral size
- Wide-branching corals support 34% higher abundance than tight-branching

### Neighborhood Effects (H3 partially supported)

- **Neighbor density**: Positive effect on CAFI abundance (facilitation/spillover)
- **Neighbor volume**: Marginal positive effect (p = 0.068)
- **Isolation**: No significant effect at meter scales (p = 0.96)
- **Relative size**: Larger-than-neighbors corals support more CAFI

The positive neighbor density effect suggests spillover/facilitation may counteract propagule dilution at local scales.

### Spatial Structure

- Significant positive spatial autocorrelation (Moran's I = 0.23, p < 0.001)
- Strongest clustering at 10–50 m scales
- Crabs: I = 0.31 (strongest); Fish: I = 0.08 (not significant)

### Condition Relationships (H4 supported)

- Coral condition positively predicts CAFI diversity (β = 0.19, p < 0.01)
- Effect independent of coral size
- Suggests bidirectional coral-CAFI interactions

## Implications

### For Theory

- Scaling exponent supports propagule redirection predictions
- Positive neighbor effects suggest facilitation at fine scales
- Results inform models of coral-CAFI dynamics

### For Conservation

- Maintain coral populations with diverse size structures
- Coral density within patches may matter more than spacing
- Branch architecture influences CAFI more than putative species ID

### Data Notes

1. **Pocillopora taxonomy**: Species cannot be reliably distinguished morphologically. We focus on measurable architectural traits (branch width).

2. **OTU approach**: CAFI are field-identified morphological OTUs, not genetically confirmed species.

3. **Physiology subset**: Not all corals have physiological measurements.

4. **Position correction**: Sampling position on branches correlates with colony size (r = 0.565). We address this confound by regressing each physiological trait on stump length and using residuals as position-corrected values. The condition score (PC1 of corrected traits) shows minimal correlation with colony volume (|r| < 0.10).

## Sites

| Code | Name | Environment | N Corals |
|------|------|-------------|----------|
| HAU | Hauru | Fringing reef (north shore) | 68 |
| MAT | Maatea | Lagoon (back reef) | 68 |
| MRB | Moorea Barrier Reef | Barrier reef (outer) | 68 |

## Output Summary

| Type | Count | Description |
|------|-------|-------------|
| Figures | 114+ | PNG visualizations (300 dpi) |
| Tables | 77+ | CSV statistical results |
| R Objects | 11 | RDS datasets and models |
| Manuscript | 1 | MEPS-formatted publication |

## Research Workflow Agents

The `agents/` directory contains Python-based research workflow agents using Claude:

| Agent | Purpose |
|-------|---------|
| Research PRD | Structured problem definition |
| Literature | Gap analysis & synthesis |
| Data QA | Quality assessment |
| EDA | Exploratory analysis |
| Modeling | Statistical planning |
| Figure Factory | Visualization design |
| Scientific Writer | Manuscript sections |

## Citation

If you use this analysis pipeline or data, please cite:

```
Stier AC, et al. (2026). Nonlinear scaling and neighborhood effects
structure coral-associated fauna communities in Mo'orea, French Polynesia.
Marine Ecology Progress Series [submitted].
```

## Funding

This research was supported by the National Science Foundation.

## License

This project is licensed under the MIT License.

## Contact

- **Lab**: Stier Lab, Department of Ecology, Evolution, and Marine Biology, UC Santa Barbara
- **Corresponding author**: Adrian Stier (astier@ucsb.edu)
- **Issues**: [GitHub Issues](https://github.com/stier-lab/CAFI-Survey-2026/issues)

---

*Analysis pipeline for manuscript submitted to Marine Ecology Progress Series*
