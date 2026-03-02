# METHODS.md - Statistical Methods for Supplementary Materials

## CAFI Survey Analysis - Detailed Statistical Methods

This document provides detailed statistical methods suitable for inclusion in supplementary materials or methods sections of the manuscript.

---

## 1. Study Design and Sampling

### 1.1 Study Sites
- **Location**: Mo'orea, French Polynesia (17°30'S, 149°50'W)
- **Sites**: 3 reef environments
  - HAU (Hauru): Fringing reef, north shore (n = 38 corals)
  - MAT (Maatea): Lagoon, back reef (n = 39 corals)
  - MRB (Maharepa): Outer barrier reef, oceanic exposure (n = 35 corals)

### 1.2 Survey Design
- **Sample size**: 114 *Pocillopora* coral colonies (112 with volume data; 2 excluded for missing measurements)
- **Survey types**:
  - Neighborhood surveys (n = 61): Full 5 m radius surveys of all neighboring Pocillopora
  - Size-only surveys (n = 53): Coral measurements without neighborhood data
- **CAFI collection**: All coral-associated fauna extracted via coral fragmentation
- **Total CAFI**: 3,989 individuals across 87 morphological OTUs

---

## 2. Data Processing

### 2.1 Coral Volume Estimation
Coral volume estimated using ellipsoid approximation:

```
V = (4/3) * pi * (L/2) * (W/2) * (H/2)
```

Where L, W, H are the maximum length, width, and height of the colony in cm. Volume ranged from 21 to 42,333 cm³ (CV = 119%).

### 2.2 Community Matrix Construction
- **Format**: Corals (rows) × Species (columns)
- **Values**: Abundance counts (number of individuals)
- **Zero-CAFI corals**: Included as zero-abundance rows for all corals (ensures community matrix dimensions match coral_master)
- **Filtering**: Species present on < 3 corals excluded from some analyses
- **Transformation**: Hellinger transformation for ordination analyses

### 2.3 Diversity Metrics
- **Species richness**: Number of unique OTUs per coral
- **Shannon diversity**: H' = -sum(p_i × ln(p_i)), where p_i is proportional abundance
- **Simpson diversity**: D = 1 - sum(p_i²)
- **Rarefied richness**: Expected species richness at n = 20 individuals (removes abundance confound from richness estimates)

### 2.4 Position Correction (Physiology)
Physiological traits (protein, carbohydrates, zooxanthellae density, AFDW) are regressed on `stump_length + nubbin_length` before PCA, removing a tissue-gradient sampling artifact: smaller corals required sampling more of the branch, including low-density tips. Residuals from these regressions are used as position-corrected condition metrics.

---

## 3. Statistical Models

### 3.1 Species-Area Scaling (Script: 05_species_scaling_analysis.R)

**Model specification**:

The relationship between habitat size and inhabitant abundance follows a power law:

```
N = a × V^β
```

In log-log space: ln(N) = ln(a) + β × ln(V)

**GLM implementation**:

```r
model <- MASS::glm.nb(abundance ~ log(volume) + site, data = data)
```

- **Distribution**: Negative binomial (handles overdispersion in count data)
- **Link function**: Log (natural log)
- **Fixed effects**: log(volume), site (fixed effect; k = 3 sites insufficient for random intercepts; Bolker et al. 2009)
- **Response variables**: Total CAFI, species richness, individual species abundance
- **Bootstrap CI**: 1,000 site-stratified bootstrap iterations for 95% CI on scaling exponent

**NOTE**: All models use natural log (`log()` in R). Previous versions used `log10()`, inflating coefficients by ln(10) ≈ 2.303. Corrected 2026-01-30.

**Hypothesis testing** (β vs 1):

```
z = (β - 1) / SE(β)
p = 2 × pnorm(-|z|)
```

- H0: β = 1 (Field of Dreams hypothesis)
- HA: β ≠ 1

**Multiple testing correction**: FDR (Benjamini-Hochberg) applied within category (species-level and group-level scaling tests separately).

**Interpretation**:
| Outcome | Interpretation |
|---------|----------------|
| β < 1 (p < 0.05) | Propagule redistribution / dilution effect |
| β ≈ 1 (p ≥ 0.05) | Field of Dreams / proportional scaling |
| β > 1 (p < 0.05) | Super-linear scaling / aggregation |

### 3.2 Landscape Effects Models (Script: 04_landscape_effects.R)

**Predictor selection**:

Non-redundant predictor subset (all VIF < 2):
1. `log(volume)` - Focal coral size (primary predictor)
2. `n_neighbors` - Number of neighbors within 5 m
3. `total_neighbor_volume` - Summed volume of all neighbors
4. `mean_neighbor_dist` - Mean distance to neighbors (cm)

**Model specification**:

```r
# CAFI abundance
model <- MASS::glm.nb(total_cafi ~ log(volume) + n_neighbors +
                       total_neighbor_volume + mean_neighbor_dist + site,
                       data = landscape_data)

# Species richness
model <- glm(otu_richness ~ log(volume) + n_neighbors +
             total_neighbor_volume + mean_neighbor_dist + site,
             family = poisson, data = landscape_data)

# Shannon diversity
model <- lm(shannon ~ log(volume) + n_neighbors +
            total_neighbor_volume + mean_neighbor_dist + site,
            data = landscape_data)
```

- **Site**: Fixed effect (k = 3 insufficient for random intercepts; Bolker et al. 2009)
- **NB convergence**: Poisson fallback with logging when NB fails to converge

**Effect size calculation**:

Standardized coefficients computed by centering and scaling continuous predictors (mean = 0, SD = 1) before model fitting. Partial standardized β reported in SD units.

### 3.3 Community Composition (Script: 02_community_analysis.R)

**PERMANOVA**:

Tests for multivariate differences in community composition among groups.

```r
model <- vegan::adonis2(community_matrix ~ site * size_class,
                        data = metadata,
                        method = "bray",
                        permutations = 999,
                        by = "margin")  # Type III (marginal) sums of squares
```

- **Distance metric**: Bray-Curtis dissimilarity
- **Permutations**: 999
- **Effect size**: R² (proportion of variance explained)
- **Type**: Marginal (Type III) to avoid order-dependence of predictors

**Composition divergence**:

```r
# betadisper: distance to group centroid by size class
bd <- vegan::betadisper(bray_dist, size_class)
permutest(bd)

# Rarefaction control: re-test on abundance-equalized data (100 draws, averaged)
```

Size-divergence is NOT significant after rarefaction (p = 0.61) — the raw pattern is an abundance artifact.

**NMDS Ordination**:

```r
nmds <- vegan::metaMDS(community_matrix,
                       distance = "bray",
                       k = 2,
                       trymax = 100,
                       autotransform = FALSE)
```

- **Dimensions**: 2
- **Convergence criterion**: stress < 0.20
- **Transformations**: Hellinger (pre-applied)

### 3.4 CAFI-Condition Feedbacks (Script: 09_cafi_condition_feedbacks.R)

**Condition metric**: PC1 from PCA on position-corrected physiology traits (protein, carbohydrates, zooxanthellae density, AFDW).

**Model specification**:

```r
# Linear model with fixed-effect site
model <- lm(condition_PC1 ~ predictor + log(volume) + site, data = data)
```

- **Response**: Continuous PC1 score (linear model appropriate; not count data)
- **Predictors tested**: Total CAFI abundance, species richness, rarefied richness, PC1_CAFI, functional group abundances, key species presence
- **Multiple testing**: FDR correction (Benjamini-Hochberg) applied within predictor families
- **Heteroscedasticity**: HC3 sandwich robust standard errors for count predictors
- **Diagnostics**: Shapiro-Wilk (normality), Breusch-Pagan (homoscedasticity), VIF (collinearity)

**Rarefied richness artifact test**:

Raw species richness shows nominal signal (p = 0.018), but richness correlates strongly with abundance (r = 0.84). Rarefied richness at n = 20 (controlling for sampling intensity) shows NO relationship with condition (p = 0.45, r with abundance = −0.05). The raw richness signal is an abundance artifact.

### 3.5 Network Analysis (Script: 06_network_analysis.R)

**Co-occurrence network construction**:

1. Convert to presence-absence matrix
2. Residualize presence on log(volume) to remove size confound
3. Filter species with ≥ 5% occurrence
4. Calculate Spearman correlations between species pairs (pairwise cor.test with FDR correction)
5. Retain edges with r > 0.3 (positive associations only; threshold confirmed by sensitivity analysis)

**Network metrics**:

| Metric | Description | Formula |
|--------|-------------|---------|
| Modularity (Q) | Strength of community structure | Q = (1/2m) × sum[(A_ij - k_i×k_j/2m) × δ(c_i, c_j)] |
| Transitivity | Clustering coefficient | C = 3 × (triangles) / (connected triplets) |
| Degree | Number of connections | k_i = sum(A_ij) |
| Betweenness | Bridge importance | B_i = sum(σ_st(i) / σ_st) |

**Community detection**:

```r
communities <- igraph::cluster_louvain(graph, weights = E(graph)$weight)
```

- **Algorithm**: Louvain (handles weighted networks, scalable, widely used in ecology)
- **Edge weights**: Correlation coefficients

**Null model comparison**:

```r
# 1000 configuration-model (degree-preserving) random networks
# Unweighted modularity used for null comparison (configuration model preserves degree, not weights)
z = (observed_modularity - mean(null)) / sd(null)
```

Hub species identified by combined degree + eigenvector centrality.

### 3.6 Spatial Autocorrelation (Script: 07_spatial_autocorrelation.R)

**Moran's I**:

Global test for spatial autocorrelation.

```r
weights <- spdep::nb2listw(neighbors, style = "W")
moran_test <- spdep::moran.mc(values, weights, nsim = 999)
```

- **Spatial weights**: Inverse distance weights (k-nearest neighbors fallback)
- **Bandwidth**: 25th percentile of pairwise distances
- **Permutations**: 999

**Interpretation**:
| Moran's I | Interpretation |
|-----------|----------------|
| I > 0 (significant) | Positive autocorrelation (clustering) |
| I ≈ 0 | Random spatial pattern |
| I < 0 (significant) | Negative autocorrelation (dispersion) |

**Local Indicators of Spatial Association (LISA)**:

```r
local_moran <- spdep::localmoran(values, weights)
```

**Cluster classification**:
- High-High: High value, high neighborhood mean
- Low-Low: Low value, low neighborhood mean
- High-Low: High value, low neighborhood mean (outlier)
- Low-High: Low value, high neighborhood mean (outlier)

**Mantel Test (Distance-Decay)**:

```r
comm_dist <- vegan::vegdist(community_matrix, method = "bray")
geo_dist <- dist(coordinates) * 111  # Convert degrees to km

mantel_result <- vegan::mantel(comm_dist, geo_dist,
                                method = "pearson",
                                permutations = 999)
```

### 3.7 Taxonomy Sensitivity Analysis (Script: 13_taxonomy_sensitivity.R)

Tests robustness of all key findings across 5 OTU resolution scenarios:

| Scenario | Description |
|----------|-------------|
| Baseline | Original 87 OTUs (species-level where possible) |
| Family-level | Aggregated to family |
| Genus-level | Aggregated to genus |
| No singletons | Remove OTUs occurring on only 1 coral |
| Presence-absence | Binary (0/1) community matrix |

**Metrics tested** (7 per scenario): abundance scaling β, richness scaling z, PERMANOVA R², betadisper F, Shannon diversity, rarefied richness, network modularity.

Pre-built scenario data created in `01_load_data.R` (stored as `taxonomy_scenario_data.rds`) and consumed by `13_taxonomy_sensitivity.R`. Results visualized as forest plot (Fig S8).

---

## 4. Model Diagnostics

### 4.1 Negative Binomial Models
- **Overdispersion**: Pearson χ²/df check; auto-switch to quasipoisson if needed
- **Residual diagnostics**: DHARMa package for simulated residuals
- **Collinearity**: VIF < 5 for all predictors
- **Influential points**: Cook's distance (proper formula for GLMs)

### 4.2 Linear Models
- **Residual normality**: Shapiro-Wilk test, Q-Q plots
- **Homoscedasticity**: Breusch-Pagan test; HC3 robust SEs when violated
- **Influential points**: Cook's distance < 1

### 4.3 PERMANOVA Assumptions
- **Multivariate dispersion**: betadisper() test
- **Sample size balance**: Site sample sizes similar (35-39)
- **Type III**: Marginal sums of squares to avoid order-dependence

---

## 5. Multiple Comparison Corrections

| Context | Method | Script |
|---------|--------|--------|
| CAFI-condition predictors | FDR (Benjamini-Hochberg) within predictor families | `09` |
| Key species effects | FDR across 6 species tests | `09` |
| Scaling (species/group) | FDR within category | `05` |
| Network edge significance | Pairwise cor.test p-values + FDR | `06` |
| Condition → CAFI direction | FDR across predictor families | `09` |

---

## 6. Software Environment

### 6.1 R Version
- R 4.0+ (tested on R 4.4.x)

### 6.2 Key Packages

| Package | Purpose |
|---------|---------|
| tidyverse | Data manipulation and visualization |
| vegan | Community ecology (PERMANOVA, NMDS, diversity) |
| MASS | Negative binomial GLMs |
| igraph | Network analysis |
| spdep | Spatial statistics (Moran's I, LISA) |
| sf | Spatial data handling |
| ggplot2 | Visualization |
| patchwork | Multi-panel figure composition |
| sandwich | HC3 robust standard errors |
| lmtest | Breusch-Pagan and coefficient tests |
| DHARMa | Simulated residual diagnostics |
| car | VIF, Type III tests |
| scales | Axis formatting |
| ggrepel | Non-overlapping text labels |

### 6.3 Reproducibility

```r
# Set random seed for reproducibility
set.seed(42)

# Session info
sessionInfo()
```

---

## 7. References

### Statistical Methods

- Anderson MJ (2001) A new method for non-parametric multivariate analysis of variance. Austral Ecology 26:32-46. [PERMANOVA]
- Anselin L (1995) Local indicators of spatial association - LISA. Geographical Analysis 27:93-115. [LISA]
- Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc B 57:289-300. [FDR correction]
- Blondel VD et al. (2008) Fast unfolding of communities in large networks. J Stat Mech P10008. [Louvain algorithm]
- Bolker BM et al. (2009) Generalized linear mixed models: a practical guide for ecology and evolution. Trends Ecol Evol 24:127-135. [RE guidance]
- Moran PAP (1950) Notes on continuous stochastic phenomena. Biometrika 37:17-23. [Moran's I]

### Ecological Theory

- Stier AC, Osenberg CW (2010) Propagule redirection: habitat availability reduces colonization and increases recruitment in reef fishes. Ecology 91:2884-2892.
- Stier AC, Osenberg CW (2024) Field of dreams or propagule redistribution? Current Biology 34:R1186-R1189.

---

*Last updated: March 2026*
*Pipeline version: 2.0*
