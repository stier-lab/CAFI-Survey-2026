# CAFI Survey Analysis 2026

**Coral-Associated Fauna Investigation (CAFI)** - Analyzing cryptic invertebrate communities in *Pocillopora* corals from Mo'orea, French Polynesia.

## Overview

This repository contains the complete analysis pipeline for investigating the ecological relationships between *Pocillopora* coral characteristics and their associated fauna communities (CAFI = Crabs, Alpheid shrimp, Fish, and snails/Invertebrates).

### Study Summary

- **Location**: Mo'orea, French Polynesia (3 reef sites)
- **Host Coral**: *Pocillopora* spp. (morphotypes: verrucosa, meandrina)
- **Survey Date**: Summer 2019
- **Sample Size**: 114 coral colonies, 2,847 CAFI individuals, 10 morphological OTUs

### Research Questions

1. What ecological and morphological factors determine CAFI community composition?
2. How do coral size, branching morphology, depth, and physiology affect fauna abundance and diversity?
3. What are the spatial patterns of species coexistence and distribution?
4. How do neighbor effects and environmental context shape CAFI communities?

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
│   ├── 03_spatial_patterns.R
│   ├── 04_diversity_analysis.R
│   ├── 05_coral_cafi_relationships.R
│   ├── 05a_coral_characteristics.R
│   ├── 06_network_analysis.R
│   ├── 07_morphotype_habitat_analysis.R
│   ├── 08_advanced_statistical_models.R
│   ├── 08_size_biomass_scaling.R
│   ├── 09_machine_learning_predictions.R
│   ├── 10_neighborhood_arrival_comparison.R
│   ├── 11_spatial_autocorrelation.R
│   ├── 12_visualization_suite.R
│   ├── 13_comprehensive_predictor_analysis.R
│   ├── 14_coral_size_neighbor_effects.R
│   ├── 15_comprehensive_visual_summary.R
│   ├── 16_generate_all_summary_figures.R
│   ├── 17_trapezid_guild_analysis.R
│   ├── run_all_survey_analyses.R
│   └── run_all_comprehensive_analyses.R
│
├── output/                 # Generated results
│   ├── figures/            # 114 PNG visualizations
│   ├── tables/             # 77 CSV result tables
│   ├── objects/            # 11 RDS R data objects
│   └── reports/            # Documentation and summaries
│
├── docs/                   # Project documentation
│   └── PRD.md              # Product Requirements Document
│
└── README.md
```

## Quick Start

### Prerequisites

**R packages** (installed automatically by `00_load_libraries.R`):

- Core: `tidyverse`, `here`, `readxl`, `janitor`
- Statistics: `vegan`, `lme4`, `lmerTest`, `emmeans`, `car`
- Visualization: `ggplot2`, `patchwork`, `viridis`, `corrplot`
- Spatial: `sf`, `sp`, `geosphere`, `ape`
- ML: `randomForest`, `xgboost` (optional)

**Python packages** (for research agents):

```bash
pip install anthropic pyyaml
```

### Running the Analysis

**Full pipeline:**

```r
source("scripts/run_all_survey_analyses.R")      # Scripts 00-05a
source("scripts/run_all_comprehensive_analyses.R") # Scripts 06-17
```

**Individual scripts:**

```r
source("scripts/00_load_libraries.R")
source("scripts/01_load_clean_data.R")
source("scripts/02_community_composition.R")
# ... etc
```

## Analysis Pipeline

### Data Processing (Scripts 00-01)
- Load and configure 30+ R packages
- Merge 3 data sources (CAFI, coral, physiology)
- Create community matrices and environmental datasets

### Ecological Analysis (Scripts 02-05)
- **Community composition**: Species rankings, taxonomic breakdowns
- **Spatial patterns**: Depth gradients, site distributions
- **Diversity analysis**: Alpha/beta diversity, NMDS ordination, PERMANOVA
- **Coral-CAFI relationships**: Association patterns, preferences

### Advanced Analysis (Scripts 06-11)
- **Network analysis**: Co-occurrence patterns, species modules
- **Morphotype-habitat**: Branch width effects on fauna
- **Statistical models**: GLMMs, dbRDA, variance partitioning
- **Size scaling**: Allometric relationships
- **Machine learning**: Random Forest, XGBoost predictions
- **Spatial autocorrelation**: Moran's I, LISA clustering

### Synthesis (Scripts 12-17)
- **Visualization suite**: Integrated dashboards
- **Comprehensive predictors**: Full 30+ variable analysis
- **Size-neighbor effects**: Interaction modeling
- **Trapezid guild**: Crab-specific patterns

## Key Findings

### Size Effects
- Power law scaling between coral volume and CAFI abundance (~0.75 exponent)
- Surface area and volume both significant predictors

### Neighbor Effects
- Optimal intermediate distances maximize CAFI abundance
- Wide-branching neighbors increase CAFI diversity

### Spatial Patterns
- Significant spatial autocorrelation (Moran's I p < 0.001)
- Distance decay of community similarity

### Morphotype Differences
- Verrucosa morphotype hosts higher CAFI abundance
- Wide branches consistently support more diversity

## Data Notes

### Important Caveats

1. **OTU Approach**: CAFI "species" are field-identified morphological OTUs (Operational Taxonomic Units), NOT genetically confirmed species. No DNA haplotyping was performed.

2. **Coral Taxonomy**: All corals are *Pocillopora* spp. - morphotypes (verrucosa, meandrina) cannot be reliably distinguished in field conditions without genetic confirmation.

3. **Branch Width**: The primary measurable coral morphological trait used in analyses (tight vs. wide).

4. **Physiology Subset**: Not all 114 corals have physiology measurements (zooxanthellae, chlorophyll).

## Research Workflow Agents

The `agents/` directory contains Python-based research workflow agents using Claude:

| Agent | Purpose |
|-------|---------|
| Research PRD | Structured problem definition |
| Literature | Gap analysis & synthesis |
| Framework | Conceptual diagrams |
| Data QA | Quality assessment |
| EDA | Exploratory analysis |
| Modeling | Statistical planning |
| Figure Factory | Visualization design |
| Scientific Writer | Manuscript sections |

See [docs/PRD.md](docs/PRD.md) for the full Product Requirements Document.

## Output Summary

| Type | Count | Description |
|------|-------|-------------|
| Figures | 114 | PNG visualizations (300 dpi) |
| Tables | 77 | CSV statistical results |
| R Objects | 11 | RDS datasets and models |
| Reports | 13 | Markdown documentation |

## Sites

| Code | Name | Location | N Corals |
|------|------|----------|----------|
| HAU | Hauru | North shore fringing reef | 38 |
| MAT | Maatea | Interior lagoon | 39 |
| MRB | Moorea Barrier Reef | Outer barrier reef | 37 |

## Citation

If you use this analysis pipeline or data, please cite:

```
Stier Lab. (2026). CAFI Survey Analysis: Coral-associated fauna communities
in Pocillopora corals from Mo'orea, French Polynesia.
GitHub: https://github.com/stier-lab/CAFI-Survey-2026
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-analysis`)
3. Commit changes (`git commit -am 'Add new analysis'`)
4. Push to branch (`git push origin feature/new-analysis`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## Contact

- **Lab**: Stier Lab, UC Santa Barbara
- **Issues**: [GitHub Issues](https://github.com/stier-lab/CAFI-Survey-2026/issues)

---

*This analysis pipeline was developed for investigating coral-associated fauna ecology in French Polynesia.*
