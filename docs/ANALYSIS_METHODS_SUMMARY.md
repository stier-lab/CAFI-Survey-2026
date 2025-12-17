# Comprehensive Summary of Analytical Methods

## CAFI Survey Analysis Pipeline - All Analytical Approaches

This document summarizes all analytical methods applied to the coral-associated fauna (CAFI) data across 47 R scripts. Organized by analysis category.

---

## 1. DATA PREPARATION & COMMUNITY MATRICES

### Scripts: `01_load_clean_data.R`, `02_community_composition.R`

**Data Processing:**
- Load CAFI abundance data with taxonomy (3,989 individuals, 243 OTUs)
- Load coral morphology, GPS, and physiological data (204 colonies)
- Merge datasets by `coral_id`
- Create species × site community matrices
- Calculate coral volume using ellipsoid approximation

**Community Matrix Transformations:**
- Hellinger transformation (`decostand(method = "hellinger")`)
- Wisconsin double standardization
- Presence/absence conversion
- Log(x+1) transformation

---

## 2. SCALING & POWER-LAW ANALYSIS

### Scripts: `05_coral_cafi_relationships.R`, `H2_power_law_scaling_analysis.R`

**Core Hypothesis (H2):** CAFI abundance scales sublinearly with coral volume (propagule redirection)

**Methods:**
1. **Log-log linear regression:**
   ```r
   log(Abundance) ~ log(Volume)
   ```
   - Scaling exponent β with 95% CI
   - Test whether CI excludes 1 (sublinear)
   - Test whether CI includes 0.75 (theoretical prediction)

2. **Mixed-effects model with site random slopes:**
   ```r
   log(Abundance) ~ log(Volume) + (1 + log(Volume) | Site)
   ```
   - Site-specific scaling exponents
   - Likelihood ratio tests vs null model

3. **Branch architecture interaction:**
   ```r
   log(Abundance) ~ log(Volume) * Branch_Type + (1|Site)
   ```

4. **Density scaling analysis:**
   ```r
   log(Density) ~ log(Volume)
   ```
   - Tests propagule dilution hypothesis
   - Expected negative slope (β - 1)

**Key Result:** β = 0.49 (95% CI: 0.40–0.57), confirming sublinear scaling

---

## 3. NEIGHBORHOOD EFFECTS ANALYSIS

### Scripts: `14_local_neighborhood_effects.R`, `10_neighborhood_arrival_comparison.R`

**Core Hypothesis (H3):** CAFI varies with local coral density and proximity

**Neighborhood Metrics Calculated:**
1. **Neighbor count** - # corals within 5m
2. **Neighbor volume** - total habitat nearby
3. **Local density** - neighbors / (π × radius²)
4. **Crowding index** - neighbor_volume / mean_distance
5. **Isolation index** - mean_distance / focal_volume^(1/3)
6. **Relative size** - focal_volume / mean_neighbor_volume
7. **Size asymmetry** - |focal - neighbor| / (focal + neighbor)
8. **Spillover potential** - total_neighbor_volume / mean_distance

**Statistical Models:**
- GAMs with negative binomial family:
  ```r
  gam(total_cafi ~ s(n_neighbors, k=5) + s(log_volume) + site, family = nb())
  ```
- Interactive effects (tensor smooths):
  ```r
  gam(total_cafi ~ te(log_volume, crowding_index, k=c(5,5)) + site)
  ```
- Taxon-specific responses (crabs, shrimp, fish, snails)

**Key Result:** <1% additional variance explained by neighborhood after controlling for coral size

---

## 4. DIVERSITY-CONDITION ANALYSIS

### Scripts: `04_diversity_analysis.R`, `diversity_condition_airtight.R`, `richness_condition_diagnostic.R`, `18_cafi_predicts_condition.R`

**Core Hypothesis (H4):** Coral condition positively predicts CAFI diversity

**Coral Condition Assessment:**
- Position correction for sampling bias (stump length correlation)
- PCA on position-corrected traits (protein, carbohydrate, zooxanthellae, AFDW)
- PC1 explains 52% variance = condition score

**Diversity Metrics:**
1. **Raw richness** - OTU count
2. **Rarefied richness** - standardized to n=10, 15, 20 individuals
3. **Residualized richness** - richness independent of abundance
4. **Shannon H'** - entropy-based diversity
5. **Simpson's D** - probability-based dominance
6. **Pielou's evenness** - H' / log(S)

**Statistical Tests:**
- Linear regression: `Condition ~ Diversity + log(Volume)`
- Multiple regression with abundance control
- Partial correlations controlling for volume
- Functional group-specific tests (guardian crabs, snapping shrimp, fish, corallivores)

**Sampling Artifact Detection:**
- Richness ~ log(Abundance): R² = 72% (species-area artifact)
- Rarefaction eliminates artifact
- Leave-one-species-out analysis
- FDR correction for multiple testing

**Key Result:** No robust diversity-condition relationship after proper corrections

---

## 5. COMMUNITY COMPOSITION & BETA DIVERSITY ANALYSIS

### Scripts: `02_community_composition.R`, `04_diversity_analysis.R`, `99_ultra_diversity_figure.R`, `archive/07_morphotype_habitat_analysis.R`

**Core Hypothesis (H1):** Community composition differs among reef sites

### 5.1 Alpha Diversity (Within-Colony)

**Metrics Calculated:**
```r
alpha_diversity <- data.frame(
  species_richness = specnumber(community_matrix),      # OTU count
  shannon = diversity(community_matrix, "shannon"),      # H' entropy
  simpson = diversity(community_matrix, "simpson"),      # 1 - D dominance
  evenness = shannon / log(species_richness)            # Pielou's J
)
```

**Site Comparisons:**
- Boxplots by site (HAU, MAT, MRB)
- ANOVA / Kruskal-Wallis tests
- Post-hoc pairwise comparisons

### 5.2 Beta Diversity (Between-Colony Dissimilarity)

**Distance Metrics Calculated:**
```r
beta_bray <- vegdist(community_matrix, method = "bray")        # Abundance-weighted
beta_jaccard <- vegdist(community_matrix, method = "jaccard", binary = TRUE)  # Presence/absence
beta_morisita <- vegdist(community_matrix, method = "morisita") # Probabilistic, robust to sample size
```

**Beta Diversity Across All Colonies:**
- Mean pairwise Bray-Curtis dissimilarity
- Distribution of dissimilarity values
- Dissimilarity vs geographic distance (Mantel test)

**Beta Diversity Within vs Among Sites:**

1. **Within-site beta diversity:**
   ```r
   # For each site, calculate mean pairwise dissimilarity
   for (site in c("HAU", "MAT", "MRB")) {
     site_comm <- community_matrix[metadata$site == site, ]
     within_beta <- mean(vegdist(site_comm, method = "bray"))
   }
   ```

2. **Among-site beta diversity:**
   ```r
   # Calculate pairwise dissimilarity between sites
   site_centroids <- aggregate(community_matrix ~ site, FUN = mean)
   among_beta <- vegdist(site_centroids, method = "bray")
   ```

3. **PERMDISP (Multivariate Dispersion Test):**
   ```r
   betadisper_site <- betadisper(beta_bray, metadata$site)
   permutest(betadisper_site, permutations = 999)
   ```
   - Tests homogeneity of multivariate dispersions
   - Significant = sites differ in community variability (not just composition)
   - Distance to centroid per site

**Ordination Methods:**
1. **NMDS (Non-metric Multidimensional Scaling):**
   ```r
   nmds_bray <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)
   ```
   - Stress < 0.1 = excellent, < 0.2 = good
   - 95% confidence ellipses by site
   - Species scores overlay

2. **PCA on Hellinger-transformed data:**
   - Compositional PC1, PC2 scores
   - Biplot with species loadings
   - Scree plot (variance explained)

3. **db-RDA (Distance-based Redundancy Analysis):**
   ```r
   dbrda_model <- dbrda(community_matrix ~ site + depth + morphotype, distance = "bray")
   ```
   - Constrained ordination with environmental predictors
   - Permutation tests for significance
   - Variance partitioning by predictor

### 5.3 Diversity Partitioning (Alpha-Beta-Gamma)

**Multiplicative Partitioning:**
```r
# Gamma (regional) diversity
gamma_div <- diversity(colSums(community_matrix), "shannon")

# Alpha (mean local) diversity
alpha_by_site <- aggregate(shannon ~ site, data = alpha_diversity, FUN = mean)
mean_alpha <- mean(alpha_by_site$mean_alpha)

# Beta (turnover) diversity
beta_div <- gamma_div / mean_alpha  # Multiplicative: how many times more diverse is region than average site

# Proportions
proportion_within <- mean_alpha / gamma_div
proportion_between <- 1 - proportion_within
```

**Additive Partitioning (Alternative):**
```r
# Beta = Gamma - mean(Alpha)
# Proportion attributed to turnover vs local diversity
```

**Hierarchical Partitioning by Morphotype:**
```r
# Within morphotype (alpha, beta, gamma)
for (morph in unique(metadata$morphotype)) {
  morph_comm <- community_matrix[metadata$morphotype == morph, ]
  alpha_morph <- mean(diversity(morph_comm))
  beta_morph <- mean(vegdist(morph_comm, "bray"))
  gamma_morph <- diversity(colSums(morph_comm))
}
```

### 5.4 Statistical Tests for Composition

**PERMANOVA (Permutational MANOVA):**
```r
# Site effect
permanova_site <- adonis2(community_matrix ~ site, data = metadata,
                          method = "bray", permutations = 999)

# Combined model
permanova_combined <- adonis2(community_matrix ~ site + depth_m + morphotype,
                              data = metadata, method = "bray", permutations = 999)

# Pairwise comparisons
pairwise.adonis(community_matrix, metadata$site, sim.method = "bray")
```

**Key Results:**
- Site effect: R² = 0.18, p < 0.001
- Pairwise: HAU vs MAT, HAU vs MRB, MAT vs MRB all significant

**ANOSIM (Analysis of Similarities):**
```r
anosim(community_matrix, metadata$site, permutations = 999)
```

**Environmental Vector Fitting:**
```r
envfit(nmds_bray, metadata[, c("depth_m", "volume", "condition_score")], permutations = 999)
```

### 5.5 Beta Diversity vs Coral Size

**Size-Driven Compositional Turnover:**
```r
# Distance from group centroid by coral size class
size_classes <- cut(metadata$volume, breaks = quantile(metadata$volume, c(0, 0.33, 0.66, 1)))
betadisper_size <- betadisper(beta_bray, size_classes)

# Correlation: dissimilarity from centroid vs volume
dist_to_centroid <- betadisper_size$distances
cor.test(dist_to_centroid, metadata$volume)
```

**Findings:**
- Larger corals have more distinct (further from centroid) communities
- Compositional turnover: more fish/gastropods in large, more crabs in small

---

## 6. NETWORK ANALYSIS

### Scripts: `06_network_analysis.R`

**Core Hypothesis (H5):** Co-occurrence networks exhibit non-random modular structure

**Network Construction:**
- Species co-occurrence matrix (positive associations only)
- Spearman correlation threshold (ρ > 0.3)
- Minimum co-occurrence frequency filter

**Network Metrics:**
1. **Global metrics:**
   - Density (edges / possible edges)
   - Transitivity / clustering coefficient
   - Mean path length
   - Diameter
   - Modularity (Q) via Louvain algorithm

2. **Node-level metrics:**
   - Degree centrality
   - Betweenness centrality
   - Closeness centrality
   - Eigenvector centrality

3. **Module analysis:**
   - Community detection (Louvain)
   - Module membership
   - Within/between-module connectivity

**Null Model Comparisons:**
- Erdős-Rényi random graphs (1000 permutations)
- Configuration model (degree-preserving)
- Observed transitivity vs null: 3.8× higher than expected

**Key Result:** Q = 0.52, 7 modules, *Alpheus diadema* as hub species (degree = 12)

---

## 7. SPATIAL ANALYSIS

### Scripts: `03_spatial_patterns.R`, `11_spatial_autocorrelation.R`

**Spatial Autocorrelation:**
1. **Moran's I (global):**
   - k-nearest neighbors weights (k=8)
   - Tests for abundance, richness, Shannon diversity

2. **LISA (Local Indicators of Spatial Association):**
   - Local Moran's I
   - Hot-spot detection (High-High, Low-Low clusters)

3. **Spatial correlograms:**
   - Moran's I at different distance classes

**Distance-Decay Analysis:**
- **Mantel test:** community dissimilarity vs geographic distance
  ```r
  mantel(comm_dist, geo_dist, method = "pearson", permutations = 999)
  ```

**Spatial Interpolation:**
- IDW (Inverse Distance Weighting)
- Prediction surfaces for abundance

**Spatial Regression Models:**
- OLS (baseline)
- Spatial Lag Model (SAR)
- Spatial Error Model (CAR)
- Model comparison via AIC

---

## 8. MACHINE LEARNING

### Scripts: `09_machine_learning_predictions.R`

**Classification (High/Low Diversity):**
1. **Random Forest:**
   - 500 trees
   - Variable importance (Mean Decrease Gini)
   - Confusion matrix, accuracy

2. **XGBoost:**
   - Gradient boosting
   - Feature importance (Gain)

3. **Ensemble:**
   - Weighted average of RF and XGBoost predictions

**Regression (CAFI Abundance):**
- Random Forest (ranger)
- RMSE, MAE, R² metrics
- Partial dependence plots

**Species Distribution Models:**
- GLM for presence/absence by species
- AUC evaluation

**Cross-Validation:**
- 10-fold CV
- Model comparison (logistic vs RF)

---

## 9. SIZE & BIOMASS SCALING

### Scripts: `08_size_biomass_scaling.R`

**Size-to-Biomass Conversion:**
- Length-weight relationships: W = a × L^b
- Taxon-specific coefficients (crabs, shrimp, fish, snails)

**Size Spectra Analysis:**
- Normalized abundance spectra (log2 bins)
- Power-law slope fitting
- Taxonomic-specific spectra

**Biomass-Abundance (B-N) Scaling:**
- Energy equivalence test (slope = 1?)
- B-N exponent estimation

**Allometric Scaling:**
- SMA (Standardized Major Axis) regression
- Test for isometry (b = 3)

**Trophic Structure:**
- Size-based trophic level assignment
- Trophic pyramid construction

---

## 10. FUNCTIONAL GROUP ANALYSIS

### Scripts: `17_trapezid_guild_analysis.R`, `13_comprehensive_predictor_analysis.R`

**Trapezid (Guardian Crab) Analysis:**
- Trapezid abundance/richness per coral
- Size-abundance relationships by functional role
- Trapezid-physiology correlations

**Functional Group Effects on Condition:**
- Guardian crabs (*Trapezia*, *Tetralia*)
- Snapping shrimp (*Alpheus*)
- Commensal fish (Gobiidae, Blenniidae)
- Corallivores (Coralliophilidae, Muricidae)

**Group Composition Ratios:**
- Proportions by functional group
- Crab-to-shrimp ratio effects

---

## 11. CORAL CONDITION VS LANDSCAPE CHARACTERISTICS

### Scripts: `19_condition_vs_landscape.R`, `05a_coral_characteristics.R`

**Purpose:** Test whether coral physiological condition varies with landscape characteristics (size, neighborhood, isolation)

### 11.1 Response Variables

**Condition PC1 (Composite Score):**
- PCA on position-corrected physiological traits
- PC1 explains 52% variance
- Positive loadings on protein, carbohydrate, zooxanthellae, AFDW

**Individual Physiological Metrics:**
- Protein (mg/cm², position-corrected)
- Carbohydrate (mg/cm², position-corrected)
- Zooxanthellae density (cells/cm², position-corrected)
- AFDW (mg/cm², position-corrected)

### 11.2 Predictor Variables

**Coral Size:**
```r
# Does condition vary with coral size?
m_size <- lm(condition_score ~ log_volume, data = analysis_data)
```

**Neighborhood Metrics:**
```r
# All 5 neighborhood metrics tested
n_neighbors           # Count within 5m
neighbor_volume       # Total habitat nearby
isolation_index       # mean_distance / volume^(1/3)
relative_size         # focal_volume / mean_neighbor_volume
crowding_index        # neighbor_volume / mean_distance
```

**Site (Environmental Context):**
```r
# ANOVA for site differences
m_site <- lm(condition_score ~ site, data = analysis_data)
anova(m_site)
```

### 11.3 Analytical Approaches

**Simple Linear Regression:**
```r
# Each predictor separately
condition_score ~ log_volume
condition_score ~ n_neighbors
condition_score ~ neighbor_volume
condition_score ~ isolation_index
# etc.
```

**Multiple Regression:**
```r
# Combined model
condition_score ~ log_volume + site + n_neighbors
```

**ANOVA for Site Effects:**
```r
# Site comparisons
m_site <- lm(condition_score ~ site, data = analysis_data)
anova(m_site)  # F-test
```

**Individual Physiology Metrics:**
```r
# Each metric vs each predictor
protein_corrected ~ log_volume
carb_corrected ~ log_volume
zoox_corrected ~ log_volume + site
afdw_corrected ~ log_volume
```

**Correlation Matrix:**
- All condition × landscape variable correlations
- Heatmap visualization

### 11.4 Key Results (n = 84 corals with condition data)

| Response | Predictor | β | R² | p-value |
|----------|-----------|---|-----|---------|
| Condition PC1 | log(Volume) | −0.10 | 0.002 | 0.70 |
| Condition PC1 | Site | — | — | 0.51 (F=0.67) |
| Condition PC1 | n_neighbors | +0.02 | 0.033 | 0.27 |
| Condition PC1 | neighbor_volume | −0.00 | 0.010 | 0.55 |
| Condition PC1 | isolation_index | +0.05 | 0.003 | 0.73 |
| AFDW | log(Volume) | **−0.53** | **0.124** | **0.001** |
| Zooxanthellae | Site | — | — | **0.002** (F=6.98) |

**Key Finding:** Coral condition (PC1) does NOT vary with size, site, or neighborhood characteristics. Only AFDW shows significant decline with size, and zooxanthellae density differs among sites.

---

## 12. MORPHOTYPE & MICROHABITAT ANALYSIS

### Scripts: `archive/07_morphotype_habitat_analysis.R`, `05a_coral_characteristics.R`

**Morphotype-Specific Communities:**
- Distribution patterns across morphotypes
- Shannon distribution entropy
- Concentrated vs Even distribution classification

**Microhabitat Utilization:**
- Levins' B (niche breadth)
- Standardized niche breadth
- Microhabitat × morphotype interactions

**Indicator Species Analysis:**
- `multipatt()` for morphotype associations
- Indicator values with permutation tests

**Diversity Partitioning:**
- Alpha diversity (within-coral)
- Beta diversity (between-coral within morphotype)
- Gamma diversity (total morphotype)

---

## 12. ADVANCED STATISTICAL MODELS

### Scripts: `archive/08_advanced_statistical_models.R`, `13_comprehensive_predictor_analysis.R`

**Generalized Linear Mixed Models (GLMMs):**
```r
glmmTMB(abundance ~ morphotype + branch_width + depth + (1|site), family = nbinom2())
```
- Negative binomial (overdispersion)
- Zero-inflated negative binomial
- Poisson (for comparison)

**Model Diagnostics:**
- DHARMa residual simulations
- Dispersion tests
- Zero-inflation tests

**Model Selection:**
- `dredge()` for all subsets
- AICc ranking
- Model averaging (`model.avg`)

**Variance Partitioning:**
```r
varpart(community, ~ coral_characteristics, ~ spatial, ~ neighbors)
```
- Unique and shared variance components

**GAMs (Generalized Additive Models):**
```r
gam(abundance ~ s(depth, k=5) + morphotype + s(lat, long, k=10) + s(site, bs="re"), family=nb())
```
- Smooth terms for non-linear relationships
- Tensor product interactions

**Distance-Based Methods:**
- MRM (Multiple Regression on Distance Matrices)
  ```r
  MRM(comm_dist ~ geo_dist + size_dist + neighbor_dist)
  ```

---

## 13. COMPOSITIONAL TURNOVER

### Scripts: `generate_composition_size_figure.R`, `generate_composition_neighborhood_figure.R`

**Size-Driven Compositional Change:**
- Taxa proportions vs coral volume
- Fish and gastropod increase in large corals
- Crab proportion decrease in large corals

**Beta Diversity vs Size:**
- Bray-Curtis dissimilarity from centroid
- Network connectivity vs coral size

---

## 14. VISUALIZATION & FIGURE GENERATION

### Scripts: `12_visualization_suite.R`, `15_comprehensive_visual_summary.R`, `16_generate_all_summary_figures.R`, `generate_manuscript_figures.R`

**Publication-Quality Figures:**
- 300 DPI PNG output
- Consistent theme (`theme_publication()`)
- Site color palette (HAU, MAT, MRB)

**Figure Types Generated:**
- Log-log scaling plots
- NMDS ordinations with ellipses
- Network diagrams
- Coefficient forest plots
- Partial dependence plots
- Spatial maps
- Heatmaps
- Ridge density plots
- Alluvial diagrams

---

## Summary Table: All Statistical Tests Applied

| Analysis Category | Tests/Methods | Key Packages |
|-------------------|---------------|--------------|
| Scaling | Log-log regression, LMM, LRT | lme4, lmerTest |
| Neighborhood | GAM, tensor products | mgcv |
| Alpha Diversity | Richness, Shannon, Simpson, Evenness | vegan |
| Beta Diversity | Bray-Curtis, Jaccard, Morisita dissimilarity | vegan |
| Within/Among Site Beta | PERMDISP, betadisper, distance-to-centroid | vegan |
| Diversity Partitioning | Multiplicative α-β-γ, additive partitioning | vegan |
| Ordination | NMDS, PCA, db-RDA | vegan |
| Composition Tests | PERMANOVA, ANOSIM, pairwise.adonis | vegan, pairwiseAdonis |
| Diversity-Condition | Rarefaction, partial correlation, LM | vegan, car |
| Condition-Landscape | LM, ANOVA, correlation matrix | base R, corrplot |
| Network | Modularity, centrality, null models | igraph |
| Spatial | Moran's I, LISA, Mantel, SAR/CAR | spdep, spatialreg |
| Machine Learning | RF, XGBoost, CV | randomForest, xgboost, caret |
| GLMMs | NB, ZINB, model selection | glmmTMB, MuMIn |
| Variance Partition | varpart, MRM | vegan, ecodist |

---

## Data Files Generated

**Tables (CSV):** 77+ output files
**Figures (PNG):** 114+ figures at 300 DPI
**R Objects (RDS):** 11 saved model objects
**Manuscript:** MANUSCRIPT.Rmd → HTML/PDF

---

*Document generated: 2025-12-02*
*CAFI Survey Analysis Pipeline v2.0*
