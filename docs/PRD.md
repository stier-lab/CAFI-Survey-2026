# Product Requirements Document (PRD)

## CAFI Survey Analysis Pipeline

**Version**: 2.0
**Date**: November 2025
**Target Journal**: Marine Ecology Progress Series (or similar)

---

## 1. Executive Summary

### Purpose

Develop a comprehensive, reproducible analysis pipeline for a manuscript testing how coral habitat attributes (size, proximity, amount) drive variation in coral-associated fauna (CAFI) density and community structure, with implications for marine landscape ecology and coral restoration.

### Theoretical Background

A critical pattern in marine landscape ecology is that occupant *abundance* scales nonlinearly with habitat amount, causing occupant *density* (per unit habitat) to decrease as habitat increases. This arises through **propagule redirection**—larvae distribute among available habitats based on chemical cue strength, so isolated habitats receive disproportionately more settlers per unit area.

Key predictions:
- **Size effect**: Larger corals → more CAFI but at *lower density*
- **Proximity effect**: Isolated corals → *higher* CAFI density
- **Amount effect**: More coral → higher abundance but *lower density*

These patterns matter because CAFI affect coral growth and survival through predator defense, sediment removal, and nutrient provision, creating bidirectional feedbacks.

### Scope

- Survey 204 coral colonies at 3 reef sites in Mo'orea
- Test scaling relationships (power-law exponents)
- Quantify meter-scale neighborhood effects
- Link CAFI patterns to coral physiological condition
- Generate publication-ready manuscript for MEPS

### Success Criteria

- Estimate scaling exponent with 95% CI
- Test neighborhood effect predictions
- Demonstrate condition-CAFI relationships
- Publication-ready figures and statistical summaries

---

## 2. Research Hypotheses

### H1: Site Effects

CAFI community composition will differ among reef sites due to variation in coral landscapes and environmental conditions.

**Test**: PERMANOVA on Bray-Curtis dissimilarity
**Metric**: R², F-statistic, pairwise comparisons

### H2: Size Scaling

CAFI abundance will scale with coral volume following a power-law relationship with exponent less than 1, indicating that larger corals have lower CAFI densities.

**Model**: log(Abundance) ~ log(Volume) + (1 + log(Volume) | Site)
**Metric**: Scaling exponent β, 95% CI, comparison to 0.75

### H3: Neighborhood Effects

CAFI abundance and diversity will vary with local coral density and proximity, reflecting propagule redirection and potential spillover effects.

**Metrics**:
- Neighbor density (count within 5m)
- Neighbor volume (total habitat nearby)
- Isolation index (distance to neighbors)
- Relative size (focal vs. neighbor volumes)

**Predictions**:
- Positive isolation effect (propagule redirection)
- Possible positive density effect (spillover/facilitation)

### H4: Condition Relationships

Coral physiological condition will positively predict CAFI diversity, consistent with bidirectional coral-CAFI interactions.

**Model**: Diversity ~ Condition + Volume + Branch Width + (1|Site)
**Metric**: β for condition, p-value, effect size

### H5: Network Structure

CAFI co-occurrence networks will exhibit non-random modular structure, with modules corresponding to functional groups or shared habitat preferences, and identifiable keystone species.

**Tests**:
- Modularity (Q) vs. null model
- Module composition by taxonomic group
- Centrality metrics (degree, betweenness)

**Predictions**:
- Significant modularity (Q > 0.3)
- Modules reflect ecological similarity
- Keystone species identifiable through centrality

---

## 3. Functional Requirements

### 3.1 Data Management

| ID | Requirement | Priority | Status |
|----|-------------|----------|--------|
| F1.1 | Load CAFI abundance data with taxonomy | Must | Complete |
| F1.2 | Load coral morphology and GPS | Must | Complete |
| F1.3 | Load coral physiology data | Must | Complete |
| F1.4 | Merge datasets by coral_id | Must | Complete |
| F1.5 | Create community matrix | Must | Complete |
| F1.6 | Calculate proximity metrics from GPS | Must | Complete |

### 3.2 Scaling Analysis (H2)

| ID | Requirement | Priority | Status |
|----|-------------|----------|--------|
| F2.1 | Power-law regression (log-log) | Must | Complete |
| F2.2 | GLMM with site random slopes | Must | Complete |
| F2.3 | Test exponent vs. theoretical 0.75 | Must | Complete |
| F2.4 | Branch architecture interaction | Should | Complete |

### 3.3 Neighborhood Analysis (H3)

| ID | Requirement | Priority | Status |
|----|-------------|----------|--------|
| F3.1 | Calculate neighbor count within 5m | Must | Complete |
| F3.2 | Calculate total neighbor volume | Must | Complete |
| F3.3 | Calculate isolation index | Must | Complete |
| F3.4 | Calculate relative size | Must | Complete |
| F3.5 | GAM models for nonlinear effects | Must | Complete |
| F3.6 | Taxon-specific responses | Should | Complete |

### 3.4 Condition Analysis (H4)

| ID | Requirement | Priority | Status |
|----|-------------|----------|--------|
| F4.1 | Diagnose sampling position bias | Must | Complete |
| F4.2 | Position correction for physiology | Must | Complete |
| F4.3 | PCA on position-corrected traits | Must | Complete |
| F4.4 | Condition score (PC1) | Must | Complete |
| F4.5 | Validate minimal size correlations | Must | Complete |
| F4.6 | Diversity ~ Condition model | Must | Complete |
| F4.7 | Control for coral size | Must | Complete |

### 3.5 Community Analysis (H1)

| ID | Requirement | Priority | Status |
|----|-------------|----------|--------|
| F5.1 | Bray-Curtis dissimilarity | Must | Complete |
| F5.2 | PERMANOVA | Must | Complete |
| F5.3 | NMDS ordination | Must | Complete |
| F5.4 | Environmental vector fitting | Should | Complete |
| F5.5 | Network modularity | Should | Complete |

### 3.6 Spatial Analysis

| ID | Requirement | Priority | Status |
|----|-------------|----------|--------|
| F6.1 | Moran's I for autocorrelation | Must | Complete |
| F6.2 | Taxon-specific autocorrelation | Should | Complete |
| F6.3 | Spatial random effects in models | Should | Complete |

### 3.7 Visualization

| ID | Requirement | Priority | Status |
|----|-------------|----------|--------|
| F7.1 | Log-log scaling plots | Must | Complete |
| F7.2 | Neighborhood effect panels | Must | Complete |
| F7.3 | NMDS ordination | Must | Complete |
| F7.4 | Network visualization | Should | Complete |
| F7.5 | Publication quality (300 dpi) | Must | Complete |

---

## 4. Key Analyses

### 4.1 Volume-Abundance Scaling

**Purpose**: Test H2—whether CAFI abundance scales sublinearly with coral volume

**Model**:
```
log(Abundance) ~ log(Volume) + Branch Width + (1 + log(Volume) | Site)
```

**Key outputs**:
- Scaling exponent β with 95% CI
- Test whether CI excludes 1 (sublinear)
- Test whether CI includes 0.75 (theoretical prediction)

**Expected result**: β ≈ 0.75–0.85

### 4.2 Local Neighborhood Effects

**Purpose**: Test H3—how meter-scale context affects CAFI

**Metrics**:
```r
local_density = n_neighbors / (π × radius²)
crowding_index = total_neighbor_volume / mean_distance
isolation_index = mean_distance / focal_volume^(1/3)
relative_size = focal_volume / mean_neighbor_volume
spillover_potential = total_neighbor_volume / mean_distance
```

**Models**: GAMs for nonlinear relationships

**Expected results**:
- Positive isolation effect (propagule redirection)
- Possible positive density effect (spillover)
- Relative size effects on CAFI abundance

### 4.3 Position Correction & Coral Condition

**Purpose**: Address sampling position bias and test H4—whether coral health predicts CAFI diversity

**Problem**: Sampling position (stump length) correlates with colony size (r = 0.565), confounding size effects with positional effects on physiology.

**Position correction methodology**:
```r
# For each physiological trait:
position_model <- lm(trait ~ stump_length)
corrected_trait <- residuals(position_model)
corrected_trait_z <- scale(corrected_trait)  # Standardize to z-scores
```

**Condition score**: PC1 from position-corrected protein, carbohydrate, zooxanthellae, AFDW (~60% variance explained)

**Validation**: Corrected traits show |r| < 0.10 with colony volume

**Model**:
```
Shannon ~ Condition + log(Volume) + Branch Width + (1|Site)
```

**Expected result**: Positive β for condition

---

## 5. Data Requirements

### Input Data

| File | Records | Variables | Description |
|------|---------|-----------|-------------|
| survey_cafi_data_w_taxonomy_summer2019_v5.csv | 2,847 | ~10 | CAFI individuals |
| survey_coral_characteristics_merged_v2.csv | 204 | ~30 | Coral morphology, GPS |
| survey_master_phys_data_v3.csv | ~100 | ~15 | Physiology subset |

### Key Variables

**Predictors**:
- `coral_volume`: cm³
- `n_neighbors`: count within 5m
- `mean_neighbor_dist`: meters
- `total_neighbor_volume`: cm³
- `branch_width`: tight/wide
- `condition_score`: PC1

**Responses**:
- `cafi_abundance`: count per coral
- `cafi_density`: abundance/volume
- `shannon_diversity`: H'
- `cafi_richness`: OTU count

### Output Data

| Type | Count | Format |
|------|-------|--------|
| Figures | 114+ | PNG (300 dpi) |
| Tables | 77+ | CSV |
| R Objects | 11 | RDS |
| Manuscript | 1 | HTML/PDF |

---

## 6. Technical Architecture

### Technology Stack

- **R 4.x**: tidyverse, vegan, lme4, glmmTMB, mgcv
- **Python 3.x**: Research agents (anthropic API)
- **Git**: Version control

### Pipeline Architecture

```
[Data Loading] → [GPS Processing] → [Community Matrices]
       ↓
[Scaling Analysis] → [Power-law Exponents] → [Fig. 2]
       ↓
[Neighborhood Analysis] → [GAM Models] → [Fig. 6]
       ↓
[Condition Analysis] → [CAFI-Physiology Links]
       ↓
[Community Analysis] → [PERMANOVA, NMDS] → [Fig. 1]
       ↓
[Manuscript] → [MEPS Format]
```

### Key Scripts

| Script | Purpose | Hypothesis |
|--------|---------|------------|
| 05_coral_cafi_relationships.R | Scaling analysis | H2 |
| 14_local_neighborhood_effects.R | Neighborhood effects | H3 |
| 04_diversity_analysis.R | Community/condition | H1, H4 |
| 06_network_analysis.R | Co-occurrence | H1 |

---

## 7. Manuscript Structure

### Target Journal

Marine Ecology Progress Series (MEPS)
- Format: Abstract, Introduction, Methods, Results, Discussion
- Style: Clean academic, Times New Roman
- Figures: Embedded at end

### Main Figures

1. **Fig. 1**: Study design and community composition (NMDS)
2. **Fig. 2**: Volume-abundance scaling (log-log)
3. **Fig. 3**: Diversity patterns across gradients
4. **Fig. 4**: Spatial autocorrelation
5. **Fig. 5**: Co-occurrence network
6. **Fig. 6**: Neighborhood effects panel

### Key Results to Report

- Scaling exponent: 0.81 ± 0.12 (95% CI: 0.58–1.04)
- Site effect: R² = 0.18, p < 0.001
- Branch architecture: 34% higher abundance in wide
- Condition effect: β = 0.19, p < 0.01
- Network modularity: Q = 0.42, p < 0.001
- Neighbor density: positive effect
- Isolation: no significant effect (p = 0.96)

---

## 8. Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| GPS coverage incomplete | Can't test proximity | Report subset (n=60) |
| Physiology subset small | Low power for H4 | Report effect sizes |
| Species ID unreliable | Can't test morphotype | Use branch architecture |
| Autocorrelation | Invalid p-values | Spatial random effects |
| Facilitation > dilution | Counter to prediction | Discuss mechanisms |

---

## 9. Timeline

### Completed

- [x] Data loading and integration
- [x] Scaling relationship analysis
- [x] Neighborhood effects analysis
- [x] Community composition (PERMANOVA, NMDS)
- [x] Network analysis
- [x] Spatial autocorrelation
- [x] Manuscript draft (MEPS format)

### In Progress

- [ ] Final figure polish
- [ ] Co-author review
- [ ] Supplement preparation

### Future

- [ ] Journal submission
- [ ] Peer review revisions
- [ ] Data archiving (Dryad/Zenodo)

---

## 10. Success Metrics

### Quantitative

| Metric | Target | Result |
|--------|--------|--------|
| Scaling exponent CI includes 0.75 | Yes | 0.58–1.04 ✓ |
| Neighborhood effects tested | 5 | 5 ✓ |
| Condition-diversity p < 0.05 | Yes | p = 0.01 ✓ |
| Figures publication-ready | Yes | Yes ✓ |
| Manuscript complete | Yes | Yes ✓ |

### Qualitative

- Results interpretable in propagule redirection framework
- Discussion connects to restoration implications
- Methods appropriate for MEPS standards

---

## 11. References

Key literature informing this work:

- Abele & Patton (1976) - Size effects on coral commensals
- Stier & Osenberg (2010) - Propagule redirection
- Hamman et al. (2018) - Spatial models of coral-CAFI
- McKeon et al. (2012) - Multiple defender effects
- Silliman et al. (2015) - Facilitation in restoration

---

*Document prepared for CAFI Survey Analysis Pipeline v2.0*
*Manuscript: Nonlinear scaling and neighborhood effects structure coral-associated fauna communities*
