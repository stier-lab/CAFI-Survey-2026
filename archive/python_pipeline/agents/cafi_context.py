"""
CAFI-specific domain context for research workflow agents.

This module provides domain knowledge about the CAFI (Coral-Associated Fauna Investigation)
survey project to inform agent prompts and outputs.
"""

CAFI_CONTEXT = """
## CAFI Project Context

### What is CAFI?
**CAFI** = Crabs, Alpheid shrimp, Fish, and snails (Invertebrates) - cryptic fauna living in Pocillopora coral colonies.

### Study System
- **Location**: Mo'orea, French Polynesia
- **Sites**:
  - HAU (Hauru) - North shore fringing reef, 38 corals
  - MAT (Maatea) - Interior lagoon, 39 corals
  - MRB (Moorea Barrier Reef) - Outer barrier reef, 37 corals
- **Host Coral**: Pocillopora spp. (morphotypes: verrucosa, meandrina, eydouxi)
- **Survey Date**: Summer 2019 (June-August)
- **Sample Size**: 204 coral colonies, 2,847 CAFI individuals, 243 OTUs

### Key Variables
**Response Variables:**
- CAFI abundance (count per coral)
- CAFI richness (number of OTUs per coral)
- Shannon diversity index
- Community composition (Bray-Curtis dissimilarity)

**Predictor Variables:**
- Coral volume (cm³) - ranges 50-5,000 cm³
- Coral surface area (cm²)
- Depth (m) - ranges 2-15 m
- Site (categorical: HAU, MAT, MRB)
- Morphotype (categorical: verrucosa, meandrina, eydouxi)
- Branch width (categorical: tight, wide)

**Physiological Variables (position-corrected):**
- Protein content (mg/cm²)
- Carbohydrate content (mg/cm²)
- Zooxanthellae density (cells/cm²)
- AFDW - ash-free dry weight (mg/cm²)
- Condition score (PC1 from PCA, explains ~60% variance)

### Key Statistical Findings

**Community Structure:**
- Top 20 OTUs comprise >80% of community
- Crabs dominate, followed by shrimp > fish > snails
- High beta diversity (Bray-Curtis ~0.55 between corals)
- Significant site effects on composition (PERMANOVA p < 0.05)

**Diversity Metrics:**
- Mean richness: 7-10 OTUs per coral
- Shannon diversity: ~2.0
- Positive spatial autocorrelation (Moran's I p < 0.05)

**Network Properties:**
- Modular structure (modularity 0.3-0.5)
- High clustering (transitivity 0.3-0.5)
- Keystone species identified via degree × betweenness

**Coral-CAFI Relationships:**
- Larger corals host more CAFI (power law scaling ~0.75)
- Branch width affects community composition
- Position-corrected condition score integrates physiology

### Important Caveats

1. **OTU Approach**: CAFI "species" are field-identified morphological OTUs, NOT genetically confirmed species. No DNA sequencing was performed.

2. **Coral Taxonomy**: All corals are Pocillopora spp. - morphotypes cannot be reliably distinguished without genetic confirmation.

3. **Branch Width**: The primary measurable coral morphological trait (tight vs. wide branching).

4. **Position Bias**: Larger corals were systematically sampled farther from base (r = 0.565). All physiological analyses use position-corrected values.

5. **Physiology Subset**: Not all 204 corals have complete physiology measurements.

### Analysis Pipeline Summary

**Data Processing (Scripts 00-01):**
- Load 3 CSV files, merge by coral_id
- Create 204 × 243 community matrix
- Calculate diversity metrics per coral

**Community Analysis (Scripts 02-04):**
- Rank-abundance curves
- NMDS ordination (stress < 0.2)
- PERMANOVA for site/morphotype effects
- Alpha/beta diversity quantification

**Coral-CAFI Relationships (Scripts 05-05a):**
- Position-corrected condition score
- CAFI abundance/richness vs coral traits
- GLMs with site controls

**Network Analysis (Script 06):**
- Co-occurrence networks (Spearman |r| > 0.3)
- Module detection (Louvain algorithm)
- Keystone species identification

**Advanced Statistics (Scripts 08-14):**
- GLMMs with site random effects
- Machine learning (Random Forest, XGBoost)
- Spatial autocorrelation (Moran's I, LISA)
- Size-neighbor interaction effects

### Output Files Available

**Key RDS Objects:**
- `survey_master_data.rds` - Complete dataset (204 corals × 300+ variables)
- `community_matrix.rds` - 204 × 243 species abundance matrix
- `coral_condition_scores.rds` - Position-corrected condition metric
- `machine_learning_models.rds` - Trained RF/XGBoost models

**Key Tables:**
- `species_abundance_ranking.csv` - OTU rankings
- `alpha_diversity_metrics.csv` - Per-coral diversity indices
- `cafi_keystone_species.csv` - Network centrality metrics
- `cafi_condition_model_coefficients.csv` - Statistical results

**Key Figures (114 total):**
- Community composition plots
- NMDS ordinations
- Network visualizations
- Model diagnostics
"""

# Specific context for different agent types
LITERATURE_CONTEXT = """
### Relevant Literature Topics

**Coral-Associated Fauna:**
- Cryptic reef biodiversity (Stella et al. 2011)
- Pocillopora-associated communities (Coles 1980)
- Guardian crabs (Trapezia spp.) mutualism (Glynn 1983)
- Coral morphology effects on fauna (Vytopil & Willis 2001)

**Community Ecology:**
- Species-area relationships on coral hosts
- Network approaches in mutualistic systems
- Beta diversity on coral reefs
- Spatial autocorrelation in marine communities

**Coral Physiology:**
- Coral health indicators (protein, carbs, zooxanthellae)
- Position effects on coral tissue
- Size-scaling of physiological traits

**Study Location:**
- Mo'orea LTER site
- French Polynesian reef ecology
- Pocillopora dominance in Pacific reefs
"""

MODELING_CONTEXT = """
### Statistical Modeling Recommendations

**Primary Models (Confirmatory):**
1. CAFI abundance ~ coral_volume + branch_width + depth + (1|site)
   - Distribution: Negative binomial (overdispersed counts)
   - Random effect: Site (3 levels)

2. CAFI richness ~ coral_volume + morphotype + (1|site)
   - Distribution: Poisson or negative binomial

3. Community composition: PERMANOVA
   - Distance: Bray-Curtis
   - Factors: site, morphotype, depth
   - Permutations: 999

**Secondary Models (Exploratory):**
4. Shannon diversity ~ condition_score + branch_width
   - Distribution: Gaussian

5. Spatial regression: SAR/CAR models
   - Address spatial autocorrelation

6. Machine learning: Random Forest
   - Feature importance ranking
   - Cross-validation for prediction

**Model Diagnostics:**
- DHARMa residuals for GLMMs
- VIF for multicollinearity
- Moran's I on residuals for spatial dependence

**Effect Reporting:**
- Effect sizes with 95% CIs
- Standardized coefficients for comparison
- R² (marginal and conditional for mixed models)
"""

FIGURE_CONTEXT = """
### Figure Design Standards

**Resolution & Format:**
- 300 dpi for publication
- PNG format
- White backgrounds
- Sans-serif fonts (10-12pt)

**Color Palettes:**
- Sites: distinct colors (e.g., #E41A1C, #377EB8, #4DAF4A)
- Morphotypes: categorical palette
- Continuous: viridis or plasma
- Functional groups: crabs (orange), shrimp (blue), fish (green), snails (purple)

**Required Main Figures:**
1. Study overview (map + design)
2. Community composition (rank-abundance + ordination)
3. Coral-CAFI relationships (scatterplots + effects)
4. Spatial patterns (Moran's I + LISA)
5. Network structure (modules + keystones)

**Panel Layout:**
- Multi-panel with (a), (b), (c) labels
- Consistent axis labels across panels
- Shared legends where possible

**Caption Elements:**
- What's plotted (variables, transformations)
- Sample sizes (N corals, N CAFI)
- Statistical results (F, p, R²)
- Main interpretation
"""

DATA_QA_CONTEXT = """
### Data Quality Checks

**Structural Checks:**
1. CAFI-coral merge: many-to-one on coral_id
2. Expected dimensions: 2,847 CAFI rows, 204 coral rows
3. Join completeness: all CAFI should match to corals

**Variable Validation:**
- coral_id: unique identifiers (format: SITE_XXX)
- depth_m: 0-20 m range
- volume_cm3: 50-5,000 cm³ (positive)
- site: categorical (HAU, MAT, MRB only)
- morphotype: categorical (verrucosa, meandrina, eydouxi)
- branch_width: categorical (tight, wide)

**Ecological Sanity:**
- Positive correlation: coral size vs CAFI abundance
- All corals should have GPS coordinates
- No corals with 0 surface area but >0 CAFI

**Missing Data:**
- Physiology: subset only (~92% complete)
- CAFI sizes: ~95% measured
- Check for systematic patterns in missingness

**Outlier Detection:**
- 3 SD threshold for continuous variables
- Flag but don't remove without ecological rationale
"""

EDA_CONTEXT = """
### Exploratory Data Analysis Plan

**Univariate:**
- Histograms: CAFI abundance, richness, coral volume
- Boxplots: by site, morphotype, branch_width
- Summary statistics: mean, SD, range, skewness

**Bivariate:**
- Scatterplots: abundance vs volume, richness vs depth
- Correlation matrix: all continuous predictors
- Chi-square: categorical associations

**Multivariate:**
- PCA: continuous coral traits
- NMDS: community composition (Bray-Curtis)
- Pairs plots: key predictor relationships

**Spatial:**
- Map of coral locations with CAFI abundance
- Variogram: spatial structure
- Moran's I: clustering test

**Patterns to Flag:**
- Non-linear relationships (need GAMs?)
- Threshold effects (breakpoints?)
- Outliers (ecological or error?)
- Multicollinearity (VIF > 5)
"""
