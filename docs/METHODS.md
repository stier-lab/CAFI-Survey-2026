# METHODS.md - Statistical Methods for Supplementary Materials

## CAFI Survey Analysis - Detailed Statistical Methods

This document provides detailed statistical methods suitable for inclusion in supplementary materials or methods sections of the manuscript.

---

## 1. Study Design and Sampling

### 1.1 Study Sites
- **Location**: Mo'orea, French Polynesia (17deg30'S, 149deg50'W)
- **Sites**: 3 reef environments
  - HAU (Hauru): Fringing reef, north shore
  - MAT (Maatea): Lagoon, back reef
  - MRB: Outer barrier reef, oceanic exposure

### 1.2 Survey Design
- **Sample size**: 114 *Pocillopora* coral colonies
- **Survey types**:
  - Neighborhood surveys (n=63): Full 5m radius surveys of all neighboring Pocillopora
  - Size-only surveys (n=51): Coral measurements without neighborhood data
- **CAFI collection**: All coral-associated fauna extracted via coral fragmentation
- **Total CAFI**: 3,989 individuals across 87 morphological OTUs

---

## 2. Data Processing

### 2.1 Coral Volume Estimation
Coral volume estimated using ellipsoid approximation:

```
V = (4/3) * pi * (L/2) * (W/2) * (H/2)
```

Where L, W, H are the maximum length, width, and height of the colony in cm.

### 2.2 Community Matrix Construction
- **Format**: Corals (rows) x Species (columns)
- **Values**: Abundance counts (number of individuals)
- **Filtering**: Species present on < 3 corals excluded from some analyses
- **Transformation**: Hellinger transformation for ordination analyses

### 2.3 Diversity Metrics
- **Species richness**: Number of unique OTUs per coral
- **Shannon diversity**: H' = -sum(p_i * ln(p_i)), where p_i is proportional abundance
- **Simpson diversity**: D = 1 - sum(p_i^2)

---

## 3. Statistical Models

### 3.1 Species-Area Scaling (Scripts: 05_species_scaling_analysis.R)

**Model specification**:

The relationship between habitat size and inhabitant abundance follows a power law:

```
N = a * V^beta
```

In log-log space: log(N) = log(a) + beta * log(V)

**GLM implementation**:

```r
model <- MASS::glm.nb(abundance ~ log10(volume) + site, data = data)
```

- **Distribution**: Negative binomial (handles overdispersion in count data)
- **Link function**: Log
- **Fixed effects**: log10(volume), site
- **Response variables**: Total CAFI, species richness, individual species abundance

**Hypothesis testing** (beta vs 1):

```
z = (beta - 1) / SE(beta)
p = 2 * pnorm(-|z|)
```

- H0: beta = 1 (Field of Dreams hypothesis)
- HA: beta != 1

**Interpretation**:
| Outcome | Interpretation |
|---------|----------------|
| beta < 1 (p < 0.05) | Propagule redistribution / dilution effect |
| beta approx 1 (p >= 0.05) | Field of Dreams / proportional scaling |
| beta > 1 (p < 0.05) | Super-linear scaling / aggregation |

### 3.2 Landscape Effects Models (Scripts: 04_landscape_effects.R)

**Predictor selection**:

Non-redundant predictor subset (all VIF < 2):
1. `log(volume)` - Focal coral size (primary predictor)
2. `n_neighbors` - Number of neighbors within 5m
3. `total_neighbor_volume` - Summed volume of all neighbors
4. `mean_neighbor_dist` - Mean distance to neighbors (cm)

**Model specification**:

```r
# CAFI abundance
model <- MASS::glm.nb(total_cafi ~ log10(volume) + n_neighbors +
                       total_neighbor_volume + mean_neighbor_dist + site,
                       data = landscape_data)

# Species richness
model <- glm(otu_richness ~ log10(volume) + n_neighbors +
             total_neighbor_volume + mean_neighbor_dist + site,
             family = poisson, data = landscape_data)

# Shannon diversity
model <- lm(shannon ~ log10(volume) + n_neighbors +
            total_neighbor_volume + mean_neighbor_dist + site,
            data = landscape_data)
```

**Effect size calculation**:

Standardized coefficients computed by centering and scaling continuous predictors (mean = 0, SD = 1) before model fitting.

### 3.3 Community Composition (Scripts: 02_community_analysis.R)

**PERMANOVA**:

Tests for multivariate differences in community composition among groups.

```r
model <- vegan::adonis2(community_matrix ~ site * size_class,
                        data = metadata,
                        method = "bray",
                        permutations = 999)
```

- **Distance metric**: Bray-Curtis dissimilarity
- **Permutations**: 999
- **Effect size**: R-squared (proportion of variance explained)

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

### 3.4 Network Analysis (Scripts: 06_network_analysis.R)

**Co-occurrence network construction**:

1. Convert to presence-absence matrix
2. Filter species with >= 5% occurrence
3. Calculate Spearman correlations between species pairs
4. Retain edges with r > 0.3 (positive associations only)

**Network metrics**:

| Metric | Description | Formula |
|--------|-------------|---------|
| Modularity (Q) | Strength of community structure | Q = (1/2m) * sum[(A_ij - k_i*k_j/2m) * delta(c_i, c_j)] |
| Transitivity | Clustering coefficient | C = 3 * (triangles) / (connected triplets) |
| Degree | Number of connections | k_i = sum(A_ij) |
| Betweenness | Bridge importance | B_i = sum(sigma_st(i) / sigma_st) |

**Community detection**:

```r
communities <- igraph::cluster_louvain(graph, weights = E(graph)$weight)
```

- **Algorithm**: Louvain (greedy modularity optimization)
- **Edge weights**: Correlation coefficients

**Null model comparison**:

```r
# Generate 1000 Erdos-Renyi random networks
for (i in 1:1000) {
  g_random <- igraph::erdos.renyi.game(n_nodes, density, type = "gnp")
  null_metrics[i] <- calculate_metrics(g_random)
}

# Z-score
z = (observed - mean(null)) / sd(null)
```

### 3.5 Spatial Autocorrelation (Scripts: 07_spatial_autocorrelation.R)

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
| I approx 0 | Random spatial pattern |
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

---

## 4. Model Diagnostics

### 4.1 Negative Binomial Models

- **Overdispersion**: theta parameter estimated from data
- **Residual diagnostics**: DHARMa package for simulated residuals
- **Collinearity**: VIF < 5 for all predictors

### 4.2 Linear Models

- **Residual normality**: Shapiro-Wilk test, Q-Q plots
- **Homoscedasticity**: Breusch-Pagan test
- **Influential points**: Cook's distance < 1

### 4.3 PERMANOVA Assumptions

- **Multivariate dispersion**: betadisper() test
- **Sample size balance**: Site sample sizes similar

---

## 5. Multiple Comparison Corrections

- **Within-analysis**: Tukey HSD for post-hoc pairwise comparisons
- **Across-hypotheses**: Results interpreted as hypothesis-generating rather than confirmatory

---

## 6. Software Environment

### 6.1 R Version
- R 4.0+ (tested on R 4.3.1)

### 6.2 Key Packages

| Package | Version | Purpose |
|---------|---------|---------|
| tidyverse | 2.0.0 | Data manipulation |
| vegan | 2.6-4 | Community ecology |
| MASS | 7.3-60 | Negative binomial GLMs |
| lme4 | 1.1-35 | Mixed effects models |
| glmmTMB | 1.1-8 | Generalized mixed models |
| igraph | 1.5.1 | Network analysis |
| spdep | 1.3-1 | Spatial statistics |
| sf | 1.0-14 | Spatial data |
| ggplot2 | 3.4.4 | Visualization |
| patchwork | 1.2.0 | Figure composition |

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
- Blondel VD et al. (2008) Fast unfolding of communities in large networks. J Stat Mech P10008. [Louvain algorithm]
- Moran PAP (1950) Notes on continuous stochastic phenomena. Biometrika 37:17-23. [Moran's I]

### Ecological Theory

- Stier AC, Osenberg CW (2010) Propagule redirection: habitat availability reduces colonization and increases recruitment in reef fishes. Ecology 91:2884-2892.
- Stier AC, Osenberg CW (2024) Field of dreams or propagule redistribution? Current Biology 34:R1186-R1189.

---

*Last updated: January 2026*
*Pipeline version: 1.0*
