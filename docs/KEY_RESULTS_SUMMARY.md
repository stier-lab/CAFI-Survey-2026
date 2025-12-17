# Key Analyses and Results Summary

## CAFI Survey Analysis - Executive Summary

**Dataset:** 114 *Pocillopora* colonies, 3 sites (HAU, MAT, MRB), 3,989 CAFI individuals, 243 OTUs

---

## HYPOTHESIS 1: Site Effects on Community Composition

**Question:** Does CAFI community composition differ among reef sites?

**Analyses:**
- PERMANOVA on Bray-Curtis dissimilarity
- NMDS ordination (k=2)
- Pairwise site comparisons
- PERMDISP (dispersion homogeneity)

**Key Results:**
| Test | Statistic | p-value |
|------|-----------|---------|
| PERMANOVA (site) | R² = 0.18 | **< 0.001** |
| HAU vs MAT | — | **< 0.05** |
| HAU vs MRB | — | **< 0.05** |
| MAT vs MRB | — | **< 0.05** |
| NMDS stress | 0.18 | (good fit) |

**Conclusion:** ✓ **SUPPORTED** - Sites harbor distinct CAFI communities; 18% of compositional variance explained by site.

---

## HYPOTHESIS 2: Size Scaling (Propagule Redirection)

**Question:** Does CAFI abundance scale sublinearly with coral volume?

**Analyses:**
- Log-log regression: log(Abundance) ~ log(Volume)
- Mixed-effects model with site random slopes
- Density scaling: log(Density) ~ log(Volume)
- Theoretical comparison to metabolic scaling (β = 0.75)

**Key Results:**
| Model | Scaling Exponent (β) | 95% CI | R² |
|-------|---------------------|--------|-----|
| Overall | **0.49** | 0.40–0.57 | 0.56 |
| By site (HAU) | 0.44 | — | — |
| By site (MAT) | 0.52 | — | — |
| By site (MRB) | 0.48 | — | — |

**Interpretation:**
- β = 0.49 means 10× larger coral → only 3.1× more CAFI (not 10×)
- CI excludes 1.0 (isometric) → **sublinear confirmed**
- CI excludes 0.75 (metabolic theory) → **steeper than predicted**
- Density declines with size (slope = −0.51)

**Conclusion:** ✓ **SUPPORTED** - Sublinear scaling consistent with propagule dilution; larger corals have lower fauna density per unit volume.

---

## HYPOTHESIS 3: Neighborhood Effects

**Question:** Does local coral neighborhood (density, isolation) affect CAFI abundance beyond focal coral size?

**Analyses:**
- 8 neighborhood metrics calculated (within 5m radius)
- GAMs with volume control
- Taxon-specific responses

**Key Results:**
| Neighborhood Metric | β | ΔR² | p-value |
|--------------------|---|-----|---------|
| Neighbor count | +0.10 | 0.006 | 0.52 |
| Relative size | +7.25 | 0.004 | 0.46 |
| Spillover potential | +9.60 | 0.008 | 0.32 |
| Isolation index | +0.04 | 0.001 | 0.76 |
| Crowding index | — | <0.01 | 0.68 |

**Interpretation:**
- All neighborhood metrics explain **<1% additional variance** after controlling for coral size
- Power analysis: could detect effects explaining ≥8% variance
- 5m radius may be too small for propagule dynamics (scales 10s–100s meters)

**Conclusion:** ✗ **NOT SUPPORTED** at measured scale - Coral size dominates; large neighborhood effects absent within 5m.

---

## HYPOTHESIS 4: CAFI Diversity Predicts Coral Condition

**Question:** Does CAFI diversity/composition predict coral physiological condition?

**Analyses:**
- Position correction for sampling bias (stump length → condition)
- PCA on corrected physiology (PC1 = condition score, 52% variance)
- Raw vs rarefied richness comparisons
- Robust correlation methods (Spearman, Kendall)
- Leave-one-species-out analysis
- Functional group-specific tests

**Key Results:**
| Diversity Metric | Controls | β | p-value | Interpretation |
|-----------------|----------|---|---------|----------------|
| Raw richness | Volume | +0.058 | 0.041 | Confounded by sampling |
| Raw richness | Volume + Abundance | +0.069 | 0.104 | Effect disappears |
| Rarefied richness (n=10) | Volume | −0.014 | **0.931** | No effect |
| Residualized richness | Volume | +0.055 | 0.170 | No pure diversity effect |
| Shannon H' | Volume | +0.087 | 0.811 | No effect |
| Evenness | Volume | −2.94 | 0.107 | No effect |
| Community PC1 (Pearson) | — | r = −0.29 | 0.007 | Driven by 3 outliers |
| Community PC1 (Spearman) | — | ρ = −0.06 | **0.61** | No robust relationship |

**Sampling Artifact:**
- Richness ~ Abundance: R² = **72%** (species-area artifact)
- Position correction verified: all corrected metrics r = 0.000 with stump length

**Conclusion:** ✗ **NOT SUPPORTED** - No robust diversity-condition relationship after proper corrections. Raw richness effect is a sampling artifact.

---

## HYPOTHESIS 5: Network Structure

**Question:** Do CAFI co-occurrence networks exhibit non-random modular structure?

**Analyses:**
- Co-occurrence network (48 species, 84 edges)
- Modularity (Louvain algorithm)
- Null model comparisons (Erdős-Rényi, 1000 permutations)
- Centrality metrics (degree, betweenness)

**Key Results:**
| Metric | Observed | Null Mean | Ratio |
|--------|----------|-----------|-------|
| Modularity (Q) | **0.52** | 0.24 | 2.2× |
| Transitivity | **0.28** | 0.074 | **3.8×** |
| n modules | 7 | — | — |
| Diameter | 8 | — | — |

**Hub Species (by degree):**
1. *Alpheus diadema* (degree = 12, betweenness = 260)
2. *Alpheus collumianus* (degree = 8)
3. *Caracanthus maculatus* (degree = 7)
4. *Trapezia serenei* (degree = 6)

**Conclusion:** ✓ **SUPPORTED** - Networks show elevated transitivity (3.8× null) and modularity; species form predictable co-occurring assemblages.

---

## ADDITIONAL FINDINGS

### Coral Condition vs Landscape Characteristics

**Question:** Does coral condition vary with size, site, or neighborhood?

| Response | Predictor | β | R² | p-value |
|----------|-----------|---|-----|---------|
| Condition PC1 | log(Volume) | −0.10 | 0.002 | **0.70** |
| Condition PC1 | Site | — | — | **0.51** |
| Condition PC1 | n_neighbors | +0.02 | 0.033 | **0.27** |
| AFDW | log(Volume) | **−0.53** | **0.124** | **0.001** |
| Zooxanthellae | Site | — | — | **0.002** |

**Conclusion:** Coral condition (PC1) does NOT vary with size, site, or neighborhood. Only AFDW shows significant size effect (larger corals have lower tissue biomass per area).

### Compositional Turnover with Size

- Larger corals: ↑ fish, ↑ gastropods
- Smaller corals: ↑ crabs (proportionally)
- Network connectivity increases with coral size
- Supports two-stage assembly model (propagule dilution → post-settlement interactions)

### Beta Diversity Patterns

- Within-site β: mean Bray-Curtis ~0.6–0.7
- Among-site β: significant (PERMANOVA R² = 0.18)
- Distance-decay: community similarity decreases with geographic distance

---

## SUMMARY: WHAT MATTERS FOR CAFI ASSEMBLY

### Factors That Matter:
1. **Coral size** - Dominant predictor; explains 56% of abundance variance
2. **Reef site** - 18% of compositional variance; distinct assemblages
3. **Network structure** - Non-random; species form predictable modules

### Factors That Don't Matter (at scales measured):
1. **Neighborhood** - <1% additional variance beyond size
2. **Coral condition** - No relationship with CAFI diversity/composition
3. **Individual diversity metrics** - Sampling artifacts when not corrected

### Practical Implications for Restoration:
- **Size distribution matters more than spacing** - many small vs few large corals
- **Mixed-size strategy** - small corals maximize density; large corals support distinct communities
- **Site selection** - different reef environments support different CAFI assemblages

---

## STATISTICAL RIGOR APPLIED

| Issue | Solution |
|-------|----------|
| Species-area sampling artifact | Rarefaction, residualization |
| Multiple testing | FDR correction |
| Outlier sensitivity | Robust correlations (Spearman, Kendall) |
| Position confound | Stump length correction for physiology |
| Non-independence | Mixed-effects models, PERMANOVA |
| Non-linear relationships | GAMs with smooth terms |
| Network null models | 1000 Erdős-Rényi permutations |

---

*Summary generated: 2025-12-02*
*CAFI Survey Analysis Pipeline v2.0*
