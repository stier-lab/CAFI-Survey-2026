# CAFI Survey Manuscript: Analysis Pipeline

## Landscape Characteristics Structure Coral-Associated Fauna Communities and Their Functional Effects on Coral Condition

**Study System:** *Pocillopora* colonies, Mo'orea, French Polynesia
**Survey Period:** June–August 2019
**Sites:** Hauru (HAU, fringing reef), Maatea (MAT, back-reef lagoon), Moorea Barrier Reef (MRB, barrier crest)

---

## Overarching Questions

1. **How do landscape characteristics (coral size, neighborhood context) shape CAFI community assembly?**
2. **Do different functional groups respond differently to landscape structure?**
3. **Do CAFI communities feed back to affect coral physiological condition?**

---

## Sample Overview

| Metric | Value |
|--------|-------|
| Coral colonies surveyed | 114 |
| Corals with CAFI data | 112 |
| Corals with physiology | 84 |
| Total CAFI individuals | 3,989 |
| CAFI species (OTUs) | 87 |
| Sites | 3 (HAU=40, MAT=37, MRB=37) |
| Coral volume range | 12–18,400 cm³ |
| Neighborhood radius | 5 m |

---

## Theoretical Framework

### Two Competing Hypotheses for Habitat-Fauna Scaling

| Hypothesis | Mechanism | Prediction | Ecological Interpretation |
|------------|-----------|------------|---------------------------|
| **Field of Dreams** | "Build it and they will come" — propagule supply unlimited, colonizers fill available space proportionally | β = 1.0 (linear scaling); density constant across sizes | Habitat amount alone determines community size |
| **Propagule Dilution** | Limited propagule pool spread across available habitat; larger patches intercept more propagules but at lower density | β < 1.0 (sublinear scaling); density decreases with size | Colonization is supply-limited, not space-limited |

**Key insight:** The scaling exponent β distinguishes mechanisms. If β = 1, larger corals simply have proportionally more fauna. If β < 1, larger corals have *disproportionately fewer* fauna per unit volume—consistent with propagule limitation.

### Functional Groups with Known Biology

| Group | Taxa | Mechanism | Expected Effect on Coral |
|-------|------|-----------|--------------------------|
| **Mutualist Defenders** | *Trapezia* crabs, *Alpheus* shrimp | Active defense against corallivores (spine-snipping, aggression); sediment removal; fat-body provisioning from polyps | **Positive**: ↑ survival, ↓ predation, ↓ sediment stress |
| **Nutrient Providers** | Resident fishes (*Paragobiodon*, *Caracanthus*, damselfishes) | Ammonium excretion enhances symbiont N-uptake; water movement from fanning | **Positive**: ↑ growth, ↑ thermal tolerance |
| **Ectoparasites** | *Coralliophila*, *Drupella* snails | Tissue consumption; potential disease vector | **Negative**: ↓ tissue, ↓ energy reserves, ↑ disease risk |
| **Other Cryptofauna** | Hermit crabs, shrimps, brittle stars, polychaetes | Varied/unknown | **Neutral/Unknown** |

---

## Key Variables

### Coral Landscape Metrics
| Variable | Description | Transformation |
|----------|-------------|----------------|
| `volume` | Colony volume (cm³), hemi-ellipsoid approximation | log₁₀ for scaling |
| `log_volume` | log₁₀(volume) | — |
| `surface_area` | Estimated surface area (cm²) | log₁₀ |
| `branch_width` | Branch architecture (tight/wide) | Factor |
| `n_neighbors` | Count of *Pocillopora* within 5m | — |
| `mean_neighbor_dist` | Mean distance to neighbors (m) | — |
| `total_neighbor_volume` | Sum of neighbor volumes (cm³) | log₁₀ |
| `isolation_index` | mean_neighbor_dist / volume^(1/3) | Higher = more isolated |
| `spillover_potential` | Σ(neighbor_volume / distance²) | Propagule pressure proxy |

### CAFI Community Metrics
| Variable | Description | Notes |
|----------|-------------|-------|
| `total_cafi` | Total individuals per coral | Primary abundance metric |
| `otu_richness` | Species count per coral | Observed richness |
| `shannon` | Shannon diversity (H') | Accounts for evenness |
| `cafi_density` | Individuals per cm³ | Standardized by volume |
| `n_defenders` | *Trapezia* + *Alpheus* count | Key mutualists |
| `n_fish` | Resident fish count | Nutrient providers |
| `n_ectoparasites` | Corallivorous snail count | Parasites |

### Coral Condition (Position-Corrected)
| Variable | Description | Notes |
|----------|-------------|-------|
| `condition_score` | PC1 of 4 physiological traits | 60.1% variance explained |
| `protein_corrected` | Tissue protein (z-score) | Residualized vs stump length |
| `carb_corrected` | Tissue carbohydrate (z-score) | Residualized |
| `zoox_corrected` | Symbiodiniaceae density (z-score) | Residualized |
| `afdw_corrected` | Ash-free dry weight (z-score) | Residualized |

**Critical Note:** Larger colonies were systematically sampled higher on branches (r = 0.565), creating artifactual size-physiology correlations. All traits residualized against sampling position before analysis.

---

# FIGURE 1 — Study System & Conceptual Framework

## Purpose
Introduce the study system, CAFI functional groups, and the conceptual framework linking habitat structure to community assembly and coral condition.

## Components
- **A.** Map of Mo'orea showing three study sites spanning reef zonation gradient
- **B.** *Pocillopora* colony illustrating CAFI microhabitat (branching architecture)
- **C.** Representative CAFI by functional group with photos
- **D.** Conceptual diagram: Field of Dreams vs Propagule Dilution predictions
- **E.** Hypothesized feedback loop: Landscape → CAFI → Coral condition

## Conceptual Hypotheses (Illustrated, not tested)
| Concept | Visualization |
|---------|---------------|
| Field of Dreams | Abundance ∝ Volume (slope = 1 on log-log) |
| Propagule Dilution | Abundance ∝ Volume^β where β < 1 (slope < 1) |
| Functional Feedbacks | Defenders → +Condition; Ectoparasites → −Condition |

## Script
`scripts/fig1_study_system.R`

---

# FIGURE 2 — CAFI Abundance vs Landscape Predictors

## Central Question
**Which landscape characteristics—coral size, neighborhood density, isolation, or habitat availability—predict CAFI abundance?**

## Figure Design
**Y-axis (all panels):** CAFI Abundance (total individuals per coral)
**X-axes:** Each panel shows a different landscape predictor

## Complete Predictor Inventory

### Category 1: Focal Coral Characteristics
| Variable | Description | Units | Transformation | Notes |
|----------|-------------|-------|----------------|-------|
| `volume` | Colony volume (hemi-ellipsoid) | cm³ | log₁₀ | **Primary predictor** |
| `log_volume` | Pre-computed log₁₀(volume) | — | none | Ready for models |
| `depth_m` | Water depth | m | none | Range: 1–8 m |
| `branch_width` | Branch architecture | tight/wide | factor | Affects interstitial space |
| `height`, `width`, `length` | Colony dimensions | cm | log₁₀ | Components of volume |

### Category 2: Neighborhood Counts (5m radius)
| Variable | Description | Units | Transformation |
|----------|-------------|-------|----------------|
| `n_neighbors` | Total *Pocillopora* neighbors | count | sqrt or none |
| `number_of_wide_branching_neighbors` | Wide-branching neighbors only | count | sqrt or none |
| `number_of_tight_branching_neighbors` | Tight-branching neighbors only | count | sqrt or none |
| `number_of_unknown_branching_neighbors` | Unknown morphotype neighbors | count | sqrt or none |

### Category 3: Neighborhood Volumes (5m radius)
| Variable | Description | Units | Transformation |
|----------|-------------|-------|----------------|
| `total_neighbor_volume` | Sum of all neighbor total volumes | cm³ | log₁₀ |
| `mean_neighbor_volume` | Mean volume per neighbor | cm³ | log₁₀ |
| `combined_total_volume_of_wide_branching_neighbors` | Total volume, wide neighbors | cm³ | log₁₀ |
| `combined_total_volume_of_tight_branching_neighbors` | Total volume, tight neighbors | cm³ | log₁₀ |
| `combined_live_volume_of_neighbors` | Live tissue volume, all neighbors | cm³ | log₁₀ |
| `combined_live_volume_of_wide_branching_neighbors` | Live tissue volume, wide neighbors | cm³ | log₁₀ |
| `combined_live_volume_of_tight_branching_neighbors` | Live tissue volume, tight neighbors | cm³ | log₁₀ |

### Category 4: Neighborhood Distance
| Variable | Description | Units | Transformation |
|----------|-------------|-------|----------------|
| `mean_neighbor_dist` | Mean distance to neighbors | cm | none |
| (distance to closest neighbor not available) | — | — | — |

### Category 5: Derived Landscape Indices
| Variable | Formula | Ecological Interpretation |
|----------|---------|--------------------------|
| `isolation_index` | mean_neighbor_dist / volume^(1/3) | Relative isolation accounting for focal size |
| `crowding_index` | total_neighbor_volume / mean_neighbor_dist | Local habitat density weighted by proximity |
| `relative_size` | volume / mean_neighbor_volume | Focal coral size relative to neighbors |
| `spillover_potential` | Σ(neighbor_volume / distance²) | Propagule pressure proxy (needs calculation) |

### Category 6: Site/Survey Design
| Variable | Description | Levels |
|----------|-------------|--------|
| `site` | Reef zone | HAU (fringing), MAT (back-reef), MRB (barrier) |
| `survey_type` | Survey protocol | neighborhood (n=40), size-only (n=74) |

### Category 7: Derived Response Variables (not predictors)
| Variable | Description | Notes |
|----------|-------------|-------|
| `total_cafi` | Total CAFI individuals | **Primary response** for Figure 2 |
| `cafi_density` | total_cafi / volume | Alternative response (density) |
| `otu_richness` | Species count per coral | Alternative response |

## Hypotheses

### H2.1: Coral Size Dominates — Sublinear Scaling (Propagule Dilution)
- **Rationale:** If propagule supply is limited, larger corals intercept more colonizers in absolute terms but fewer per unit volume. This predicts sublinear scaling (β < 1).
- **Prediction:** log(Abundance) ~ log(Volume) with **β < 1.0** and 95% CI excludes 1.0.
- **Alternative (Field of Dreams):** β = 1.0 (proportional scaling).

### H2.2: Neighbor Count Has Positive Effect (Source Population)
- **Rationale:** More neighbors = more local propagule sources. Corals in dense patches should receive more colonizers via spillover.
- **Prediction:** Positive relationship: more neighbors → more CAFI.
- **Alternative:** No effect if propagule supply is regional, not local.

### H2.3: Isolation Reduces Colonization
- **Rationale:** Isolated corals (far from neighbors) receive fewer propagules due to dispersal limitation.
- **Prediction:** Negative relationship: greater isolation → fewer CAFI.
- **Alternative:** No effect or positive effect if isolated corals "trap" passing propagules (propagule redirection).

### H2.4: Spillover Potential Predicts Abundance
- **Rationale:** Neighbors that are large and close contribute more propagules. Spillover potential integrates neighbor size and proximity.
- **Prediction:** Positive relationship: higher spillover → more CAFI.
- **Alternative:** No effect if focal coral size dominates.

### H2.5: Neighborhood Effects are Secondary to Size
- **Rationale:** Focal coral size may overwhelm neighborhood effects because habitat amount is the primary driver.
- **Prediction:** Size R² >> Neighborhood R²; adding neighborhood predictors increases R² by <5%.
- **Test:** Hierarchical model comparison.

## Statistical Approach

### Step 1: Univariate Models (One predictor at a time, controlling for site)
```r
library(MASS)
library(broom)

# Panel A: Focal coral size (primary hypothesis)
m_volume <- glm.nb(total_cafi ~ log(volume) + site, data = coral_master)

# Panel B: Neighbor count
m_neighbors <- glm.nb(total_cafi ~ n_neighbors + site, data = coral_master)

# Panel C: Total neighbor volume
m_neighbor_vol <- glm.nb(total_cafi ~ log(total_neighbor_volume + 1) + site, data = coral_master)

# Panel D: Mean neighbor distance
m_distance <- glm.nb(total_cafi ~ mean_neighbor_dist + site, data = coral_master)

# Panel E: Isolation index
m_isolation <- glm.nb(total_cafi ~ isolation_index + site, data = coral_master)

# Panel F: Crowding index
m_crowding <- glm.nb(total_cafi ~ crowding_index + site, data = coral_master)

# Panel G: Wide-branching neighbor volume
m_wide <- glm.nb(total_cafi ~ log(combined_total_volume_of_wide_branching_neighbors + 1) + site,
                 data = coral_master)

# Panel H: Tight-branching neighbor volume
m_tight <- glm.nb(total_cafi ~ log(combined_total_volume_of_tight_branching_neighbors + 1) + site,
                  data = coral_master)

# Collect all univariate results
univariate_results <- bind_rows(
  tidy(m_volume, conf.int = TRUE) %>% mutate(model = "Volume"),
  tidy(m_neighbors, conf.int = TRUE) %>% mutate(model = "Neighbor Count"),
  tidy(m_neighbor_vol, conf.int = TRUE) %>% mutate(model = "Neighbor Volume"),
  tidy(m_distance, conf.int = TRUE) %>% mutate(model = "Distance"),
  tidy(m_isolation, conf.int = TRUE) %>% mutate(model = "Isolation"),
  tidy(m_crowding, conf.int = TRUE) %>% mutate(model = "Crowding"),
  tidy(m_wide, conf.int = TRUE) %>% mutate(model = "Wide Neighbors"),
  tidy(m_tight, conf.int = TRUE) %>% mutate(model = "Tight Neighbors")
)
```

### Step 2: Hierarchical Model Comparison
```r
# Null model (site only)
m_null <- glm.nb(total_cafi ~ site, data = coral_master)

# Size-only model
m_size_only <- glm.nb(total_cafi ~ log(volume) + site, data = coral_master)

# Size + neighborhood counts
m_size_counts <- glm.nb(total_cafi ~ log(volume) + n_neighbors + site, data = coral_master)

# Size + neighborhood volumes
m_size_volumes <- glm.nb(total_cafi ~ log(volume) + log(total_neighbor_volume + 1) + site,
                         data = coral_master)

# Size + all neighborhood predictors
m_full <- glm.nb(total_cafi ~ log(volume) + n_neighbors + mean_neighbor_dist +
                 log(total_neighbor_volume + 1) + isolation_index + crowding_index + site,
                 data = coral_master)

# Compare models with likelihood ratio tests
anova(m_null, m_size_only, m_size_counts, m_size_volumes, m_full, test = "Chisq")

# Calculate pseudo-R² (McFadden's)
calc_pseudo_r2 <- function(model, null_model) {
  1 - (logLik(model)[1] / logLik(null_model)[1])
}

model_comparison <- tibble(
  model = c("Null", "Size", "Size+Counts", "Size+Volumes", "Full"),
  AIC = c(AIC(m_null), AIC(m_size_only), AIC(m_size_counts), AIC(m_size_volumes), AIC(m_full)),
  pseudo_r2 = c(0, calc_pseudo_r2(m_size_only, m_null), calc_pseudo_r2(m_size_counts, m_null),
                calc_pseudo_r2(m_size_volumes, m_null), calc_pseudo_r2(m_full, m_null))
)
```

### Step 3: Full Model with VIF Check
```r
library(car)

# Full landscape model (for VIF check and coefficient extraction)
m_full_all <- glm.nb(total_cafi ~ log(volume) + n_neighbors + mean_neighbor_dist +
                     log(total_neighbor_volume + 1) + isolation_index + crowding_index +
                     number_of_wide_branching_neighbors + number_of_tight_branching_neighbors +
                     site, data = coral_master)

# Check multicollinearity
vif_values <- vif(m_full_all)
# Flag if VIF > 5: indicates multicollinearity

# Standardized coefficients for forest plot
coral_scaled <- coral_master %>%
  mutate(across(c(n_neighbors, mean_neighbor_dist, total_neighbor_volume,
                  isolation_index, crowding_index), scale))

m_standardized <- glm.nb(total_cafi ~ log(volume) + n_neighbors + mean_neighbor_dist +
                         log(total_neighbor_volume + 1) + isolation_index + crowding_index + site,
                         data = coral_scaled)
```

### Step 4: Model Diagnostics
```r
# Overdispersion check (should be ~1 for well-fitted NB)
pearson_resid <- sum(residuals(m_full, type = "pearson")^2)
overdispersion <- pearson_resid / m_full$df.residual
cat("Overdispersion ratio:", round(overdispersion, 2), "\n")

# Influential points
cooksd <- cooks.distance(m_full)
influential <- which(cooksd > 4/nrow(coral_master))

# Zero-inflation check
prop_zeros <- mean(coral_master$total_cafi == 0)
expected_zeros <- sum(dnbinom(0, size = m_full$theta,
                              mu = fitted(m_full))) / nrow(coral_master)
cat("Observed zeros:", scales::percent(prop_zeros),
    "Expected:", scales::percent(expected_zeros), "\n")

# Residual diagnostics
par(mfrow = c(2, 2))
plot(fitted(m_full), residuals(m_full, type = "pearson"),
     xlab = "Fitted", ylab = "Pearson Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)
qqnorm(residuals(m_full, type = "pearson"))
qqline(residuals(m_full, type = "pearson"))
```

## Panel Design (8-panel figure)

**All panels share:** Y-axis = CAFI Abundance (total individuals per coral)

### Row 1: Focal Coral & Primary Landscape Predictors
| Panel | X-axis | Plot Type | Expected Pattern | Key Statistic |
|-------|--------|-----------|------------------|---------------|
| **A** | Coral Volume (log₁₀) | Scatter + power-law fit | **Strong positive**, sublinear | β ≈ 0.48, R² ≈ 0.56 |
| **B** | Neighbor Count | Scatter + regression | Weak/null | β ≈ 0.10, p > 0.20 |
| **C** | Total Neighbor Volume (log₁₀) | Scatter + regression | Weak positive or null | β ≈ 0.05, p > 0.20 |
| **D** | Mean Neighbor Distance | Scatter + regression | Weak/null | β ≈ 0, p > 0.30 |

### Row 2: Derived Indices & Branch-Type-Specific Predictors
| Panel | X-axis | Plot Type | Expected Pattern | Key Statistic |
|-------|--------|-----------|------------------|---------------|
| **E** | Isolation Index | Scatter + regression | Weak negative or null | β ≈ 0, p > 0.50 |
| **F** | Crowding Index | Scatter + regression | Weak positive or null | β ≈ 0.05, p > 0.30 |
| **G** | Wide-Branching Neighbor Volume (log₁₀) | Scatter + regression | Test morphotype-specific effects | β ≈ ?, p > 0.10 |
| **H** | Tight-Branching Neighbor Volume (log₁₀) | Scatter + regression | Test morphotype-specific effects | β ≈ ?, p > 0.10 |

### Summary Panel (Optional Panel I or separate figure)
**Forest plot of standardized coefficients** from full model:
- X-axis: Standardized β with 95% CI
- Y-axis: All predictor names (organized by category)
- Vertical reference line at β = 0
- Visual comparison of effect sizes across all predictors
- Color-coded by predictor category

### Visual Design Elements
- All panels: points colored by site (HAU, MAT, MRB)
- Panel A: Include **Field of Dreams reference line (β = 1)** as dashed gray line
- Panels B–H: Include regression line with 95% CI if p < 0.10; otherwise show flat trend line with shaded null band
- Annotation in each panel: R² and p-value (upper right corner)
- Consistent axis styling: log-scale axes labeled with original units
- Font: Publication-ready (10pt labels, 12pt axis titles)

## Expected Results

### Univariate Effects (each predictor tested alone, controlling for site)
| Panel | Predictor | Expected β | Expected p | R² | Interpretation |
|-------|-----------|------------|------------|-----|----------------|
| **A** | **Coral Volume (log₁₀)** | **+0.48** | **<0.001** | **0.56** | **Dominant predictor**; propagule dilution |
| B | Neighbor Count | +0.10 | >0.20 | <0.01 | No effect at 5m scale |
| C | Total Neighbor Volume (log₁₀) | +0.05 | >0.20 | <0.01 | Minimal habitat effect |
| D | Mean Neighbor Distance | ~0 | >0.30 | <0.01 | No isolation effect |
| E | Isolation Index | ~0 | >0.50 | <0.01 | No relative isolation effect |
| F | Crowding Index | +0.05 | >0.30 | <0.01 | No crowding effect |
| G | Wide-Branching Neighbor Volume | ~0 | >0.20 | <0.01 | No morphotype-specific effect |
| H | Tight-Branching Neighbor Volume | ~0 | >0.20 | <0.01 | No morphotype-specific effect |

### Hierarchical Model Comparison
| Model | Predictors Included | Δ Deviance | Δ R² | p (LRT) |
|-------|---------------------|------------|------|---------|
| Null | intercept + site | — | 0.00 | — |
| Size only | + log(volume) | **large** | **+0.56** | **<0.001** |
| + Counts | + n_neighbors | minimal | +0.002 | >0.30 |
| + Volumes | + total_neighbor_volume | minimal | +0.003 | >0.25 |
| + Distance | + mean_neighbor_dist | minimal | +0.001 | >0.50 |
| + Derived | + isolation_index + crowding_index | minimal | +0.002 | >0.40 |
| Full | all predictors | minimal | +0.008 | >0.30 |

### Key Synthesis
**Coral size explains ~56% of variance in CAFI abundance; all neighborhood predictors combined add <1%.** This supports the interpretation that:

1. **Propagule dilution dominates:** Colonization is driven by focal habitat size, not local neighborhood context at 5m scales.
2. **Neighborhood effects are negligible:** Neither neighbor counts, volumes, distances, nor derived indices improve prediction.
3. **Branch-type specificity doesn't matter:** Wide- vs tight-branching neighbor volumes show no differential effects.
4. **Scale mismatch possible:** Propagule supply may operate at larger spatial scales (10s–100s of meters), or post-settlement processes (survival, competition, predation) erase initial colonization patterns.

## Statistical Summary Table (for manuscript)

| Model | Predictors | AIC | Pseudo-R² | ΔR² vs Null |
|-------|------------|-----|-----------|-------------|
| Null | Intercept only | — | 0.00 | — |
| Size only | log(volume) + site | — | 0.58 | +0.58 |
| Size + Neighborhood | log(volume) + neighbors + distance + spillover + site | — | 0.59 | +0.01 |
| Neighborhood only | neighbors + distance + spillover + site | — | 0.04 | +0.04 |

## Script
`scripts/fig2_landscape_predictors.R`

---

# FIGURE 3 — Functional Group Scaling

## Central Question
**Do functional groups with different life histories show different scaling relationships with coral size?**

## Hypotheses

### H3.1: Defenders (Trapezia) Scale Near-Proportionally
- **Rationale:** *Trapezia* are obligate mutualists with strong host fidelity. Corals actively provision fat bodies; crabs defend territories. This tight mutualism may maintain near-constant crab coverage regardless of colony size.
- **Prediction:** β_defenders ≈ 0.85–1.0 (CI overlaps or approaches 1.0).
- **Mechanism:** Behavioral/physiological constraint maintains defender:host ratio.

### H3.2: Ectoparasites Show Strong Dilution
- **Rationale:** Corallivorous snails (*Coralliophila*, *Drupella*) colonize via encounter—they crawl from nearby substrates. Larger corals offer more surface but receive proportionally fewer snails because snail supply is limited.
- **Prediction:** β_ectoparasites < β_total (stronger sublinearity than overall community).
- **Mechanism:** Encounter-based colonization; no host-seeking behavior.

### H3.3: Resident Fish Scale with Space Availability
- **Rationale:** Fish require interstitial space for shelter. Larger colonies have disproportionately more complex interstitial space (scaling with surface area or volume).
- **Prediction:** β_fish ≈ 0.6–0.8 (moderate sublinearity).
- **Mechanism:** Space-limited rather than propagule-limited.

### H3.4: Species-Specific Slopes Reveal Mechanistic Diversity
- **Rationale:** Within functional groups, individual species may show varied responses based on body size, territorial behavior, or recruitment mode.
- **Prediction:** Heterogeneous β values; species with strong host fidelity (e.g., *Trapezia serenei*) have higher β than generalists.
- **Test:** Per-species GLMs; compare slopes across top 15–20 taxa.

## Statistical Approach

### Group-Level Models
```r
# Separate power-law for each functional group
m_defenders <- glm.nb(n_defenders ~ log(volume) + site, data = coral_master)
m_fish <- glm.nb(n_fish ~ log(volume) + site, data = coral_master)
m_ectoparasites <- glm.nb(n_ectoparasites ~ log(volume) + site, data = coral_master)
m_other <- glm.nb(n_other ~ log(volume) + site, data = coral_master)

# Extract slopes and compare
group_slopes <- bind_rows(
  tidy(m_defenders, conf.int = TRUE) %>% filter(term == "log(volume)") %>% mutate(group = "Defenders"),
  tidy(m_fish, conf.int = TRUE) %>% filter(term == "log(volume)") %>% mutate(group = "Fish"),
  # ...
)
```

### Species-Specific Slopes
```r
# For top 20 most abundant species
top_species <- cafi_clean %>% count(species) %>% top_n(20, n) %>% pull(species)

species_slopes <- map_dfr(top_species, function(sp) {
  sp_data <- cafi_by_coral %>% filter(species == sp)
  if(nrow(sp_data) < 20) return(NULL)  # Skip rare species

  m <- glm(abundance ~ log(volume), family = poisson, data = sp_data)
  tidy(m, conf.int = TRUE) %>%
    filter(term == "log(volume)") %>%
    mutate(species = sp, functional_group = get_group(sp))
})
```

### Deviation from Field of Dreams
```r
# Calculate expected abundance if density were constant (from smallest tertile)
mean_density_small <- coral_master %>%
  filter(size_class == "Small") %>%
  summarise(mean(cafi_density)) %>% pull()

coral_master <- coral_master %>%
  mutate(expected_abundance = mean_density_small * volume,
         deviation = total_cafi - expected_abundance,
         deviation_pct = (total_cafi - expected_abundance) / expected_abundance * 100)
```

## Panels
| Panel | Content | Key Statistic |
|-------|---------|---------------|
| **A** | 4-panel: Group abundance vs volume (log-log) with separate fits | β varies by group |
| **B** | Forest plot: Scaling exponents by group with 95% CI, reference line at β=1 | Visual comparison |
| **C** | Deviation from Field of Dreams by group across size classes | Bar plot: defenders near 0, ectoparasites strongly negative |
| **D** | Species-specific slopes (lollipop or heatmap) | Top 15-20 taxa, colored by group |

## Expected Results
| Group | Expected β | Rationale |
|-------|------------|-----------|
| Defenders (*Trapezia*, *Alpheus*) | 0.85–0.95 | Obligate mutualism; behavioral regulation |
| Resident Fish | 0.65–0.80 | Space-limited |
| Ectoparasites | 0.30–0.45 | Strong dilution; encounter-based |
| Other Cryptofauna | 0.45–0.55 | Mixed; matches overall pattern |

## Script
`scripts/fig3_functional_scaling.R`

---

# FIGURE 4 — Community Composition & Co-occurrence Networks

## Central Question
**Does CAFI species composition shift predictably with coral size, and do species co-occur non-randomly?**

## Hypotheses

### H4.1: Composition Shifts with Coral Size
- **Rationale:** Different species have different size thresholds for colonization, territorial requirements, and competitive abilities. Composition should turnover across size gradient.
- **Prediction:** PERMANOVA R² > 0.05 for log(volume) effect on Bray-Curtis dissimilarity.
- **Mechanism:** Size-dependent colonization and competitive exclusion.

### H4.2: Sites Differ in Composition (Environmental Filtering)
- **Rationale:** Sites span reef zonation gradient (fringing → barrier). Wave exposure, larval supply, and environmental conditions vary.
- **Prediction:** PERMANOVA R² > 0.03 for site effect.
- **Mechanism:** Environmental filtering and dispersal limitation.

### H4.3: Species Co-occur Non-Randomly (Modular Structure)
- **Rationale:** Some species co-occur because they share habitat preferences; others exclude each other competitively. This creates modular network structure.
- **Prediction:** Network modularity Q > null expectation (p < 0.05 vs randomized networks).
- **Mechanism:** Habitat filtering, facilitation, and competition.

### H4.4: Hub Species Connect Modules
- **Rationale:** Some generalist species should have high degree (many co-occurrences) and betweenness (connect disparate modules).
- **Prediction:** 1–3 species with degree > 10 and high betweenness centrality.
- **Candidate:** *Alpheus* shrimp (habitat generalist, commensal with multiple coral-dwellers).

## Statistical Approach

### Ordination & PERMANOVA
```r
# Hellinger transformation (appropriate for compositional data)
comm_hell <- decostand(community_matrix, method = "hellinger")

# NMDS ordination
set.seed(42)
nmds <- metaMDS(comm_hell, distance = "bray", k = 2, trymax = 100)
# Report stress; acceptable if < 0.2

# PERMANOVA: size and site effects
permanova_results <- adonis2(comm_hell ~ log(volume) + site,
                              data = coral_master,
                              permutations = 999,
                              method = "bray")

# Test for multivariate dispersion (PERMDISP)
disp <- betadisper(vegdist(comm_hell), coral_master$site)
permutest(disp, permutations = 999)
```

### Incidence Analysis
```r
# Species presence/absence by size class
incidence_tests <- cafi_presence %>%
  group_by(species) %>%
  summarise(
    small_pct = mean(present[size_class == "Small"]),
    large_pct = mean(present[size_class == "Large"]),
    chisq = chisq.test(table(present, size_class))$statistic,
    p_value = chisq.test(table(present, size_class))$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))
```

### Co-occurrence Network
```r
# Spearman correlation matrix
cor_matrix <- cor(t(community_matrix), method = "spearman")

# Threshold: |ρ| > 0.3, Bonferroni-corrected p < 0.05
# Build igraph network
g <- graph_from_adjacency_matrix(abs(cor_matrix) > 0.3, mode = "undirected", diag = FALSE)

# Module detection (Louvain algorithm)
modules <- cluster_louvain(g)
modularity(modules)

# Null model comparison (1000 randomizations)
null_Q <- replicate(1000, {
  g_null <- rewire(g, keeping_degseq(niter = 1000))
  modularity(cluster_louvain(g_null))
})
p_modularity <- mean(null_Q >= modularity(modules))

# Centrality metrics
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g)
```

## Panels
| Panel | Content | Key Statistic |
|-------|---------|---------------|
| **A** | NMDS ordination, points colored by log(volume), size gradient | Stress = 0.18 |
| **B** | NMDS ordination, points colored by site | Visual site separation |
| **C** | Species incidence heatmap: rows = species (top 25), columns = size class | Chi-square significance markers |
| **D** | Co-occurrence network: nodes = species, colored by module, sized by degree | Q = 0.52, p < 0.001 |

## Expected Results
| Hypothesis | Test | Result | Interpretation |
|------------|------|--------|----------------|
| H4.1 | PERMANOVA (volume) | R² = 0.08, p = 0.001 | **Composition shifts with size** |
| H4.2 | PERMANOVA (site) | R² = 0.06, p = 0.001 | **Sites differ in composition** |
| H4.3 | Modularity vs null | Q = 0.52, p < 0.0001 | **Non-random structure** |
| H4.4 | Hub species | *Alpheus diadema*: degree=12, betweenness=260 | **Connector species identified** |

## Script
`scripts/fig4_composition_networks.R`

---

# FIGURE 5 — Coral Condition Across Landscape Gradients

## Central Question
**Does coral physiological condition vary with colony size or neighborhood context, independent of CAFI?**

## Critical Methodological Issue: Sampling Position Bias

**Problem:** Physiological traits were measured from branch nubs sampled at variable positions along branches. Larger colonies were systematically sampled further from the branch base (r = 0.565 between volume and stump_length). Since traits vary along branches (e.g., protein declines toward tips), this creates a spurious size-physiology correlation.

**Solution:** Residualize each trait against stump_length before analysis:
```r
trait_corrected <- residuals(lm(trait ~ stump_length, data = physio_data))
```

## Hypotheses

### H5.1: Size Effects are Artifacts of Sampling Position
- **Rationale:** If apparent size-condition relationships disappear after position correction, they were artifacts.
- **Prediction:** Before correction: significant size effects. After correction: **no significant size effects** (β ≈ 0, p > 0.05).
- **Test:** Compare models before/after residualization.

### H5.2: Sites Do Not Differ in Condition
- **Rationale:** Sites span reef zones but may have similar environmental quality for coral physiology.
- **Prediction:** Site effect p > 0.05 in ANOVA on condition PC1.
- **Alternative:** Sites differ if environmental conditions strongly affect physiology.

### H5.3: Neighborhood Context Does Not Affect Condition Directly
- **Rationale:** At 5m scale, neighborhood effects on coral physiology are likely mediated through CAFI (tested in Figure 6), not direct environmental effects.
- **Prediction:** Isolation index, neighbor count do not predict condition (p > 0.05).
- **Test:** Multiple regression with landscape predictors.

## Statistical Approach

### Position Correction
```r
# Step 1: Residualize each trait against sampling position
physio_data <- physio_data %>%
  mutate(
    protein_corrected = residuals(lm(protein ~ stump_length, data = .)),
    carb_corrected = residuals(lm(carbohydrate ~ stump_length, data = .)),
    zoox_corrected = residuals(lm(zooxanthellae ~ stump_length, data = .)),
    afdw_corrected = residuals(lm(afdw ~ stump_length, data = .))
  )

# Step 2: Standardize to z-scores
physio_data <- physio_data %>%
  mutate(across(ends_with("_corrected"), scale))

# Step 3: PCA on corrected traits
pca_condition <- prcomp(physio_data %>% select(ends_with("_corrected")), scale = TRUE)
physio_data$condition_score <- pca_condition$x[,1]
# Report: % variance explained by PC1, loadings
```

### Validation: Before vs After Correction
```r
# Before correction
m_before <- lm(protein ~ log(volume), data = physio_data)  # Expect significant

# After correction
m_after <- lm(protein_corrected ~ log(volume), data = physio_data)  # Expect NS
```

### Landscape Effects on Condition
```r
# Size effect (after correction)
m_condition_size <- lm(condition_score ~ log(volume) + site, data = physio_data)

# Neighborhood effects
m_condition_landscape <- lm(condition_score ~ log(volume) + n_neighbors +
                            mean_neighbor_dist + site, data = physio_data)
```

## Panels
| Panel | Content | Key Statistic |
|-------|---------|---------------|
| **A** | PCA biplot: loadings for 4 traits, coral scores colored by size | PC1 = 60.1% variance |
| **B** | Condition score vs log(volume) — showing NO relationship after correction | β = −0.39, p = 0.18 |
| **C** | Condition by site (boxplot) | F(2,81) = 0.67, p = 0.51 |
| **D** | 4 mini-panels: each corrected trait vs volume | All p > 0.10 |
| **E** | Validation: Before/after correction comparison | r drops from 0.57 to ~0 |

## Expected Results
| Hypothesis | Test | Result | Interpretation |
|------------|------|--------|----------------|
| H5.1 | Size → Condition (corrected) | β = −0.39, p = 0.18 | **Artifact removed**; no true size effect |
| H5.2 | Site effect | F(2,81) = 0.67, p = 0.51 | **Sites similar** in condition |
| H5.3 | Neighborhood → Condition | All p > 0.30 | **No direct neighborhood effect** |

## Script
`scripts/fig5_coral_condition.R`

---

# FIGURE 6 — CAFI Effects on Coral Condition (Functional Feedbacks)

## Central Question
**Do CAFI communities—particularly functional groups—feed back to affect coral physiological condition?**

## Hypotheses

### H6.1: Mutualist Defenders Enhance Coral Condition
- **Rationale:** *Trapezia* crabs actively defend corals against corallivores, remove sediment, and receive lipid provisioning. Experimental removal studies show 2–3× higher coral mortality without defenders (Glynn 1983).
- **Prediction:** Defender abundance positively predicts condition (β > 0, p < 0.05).
- **Effect size expectation:** Small-moderate (R² = 0.03–0.10) given observational design.

### H6.2: Resident Fish Enhance Coral Condition
- **Rationale:** Coral-dwelling fish excrete ammonium that enhances symbiont nitrogen uptake, potentially buffering thermal stress (Holbrook et al. 2008).
- **Prediction:** Fish abundance positively predicts condition (β > 0).
- **Caveat:** Effect may be weak if fish are rare or effect is context-dependent.

### H6.3: Ectoparasites Reduce Coral Condition
- **Rationale:** Tissue-feeding snails directly consume coral tissue, reducing energy reserves and potentially transmitting disease.
- **Prediction:** Ectoparasite abundance negatively predicts condition (β < 0, p < 0.05).
- **Effect size expectation:** Moderate (R² = 0.05–0.15) if snail loads are high.

### H6.4: Total Abundance/Richness Do Not Predict Condition
- **Rationale:** Functional identity matters more than "how many" fauna. A coral with 50 defenders differs from one with 50 ectoparasites.
- **Prediction:** Total CAFI abundance and OTU richness show **no significant effect** (p > 0.10).
- **Mechanism:** Opposing functional effects cancel out at aggregate level.

### H6.5: Community Composition Predicts Condition
- **Rationale:** If functional identity matters, then compositional axes (which capture functional turnover) should predict condition better than abundance alone.
- **Prediction:** CAFI community PC1 predicts condition (p < 0.05, R² > abundance R²).
- **Test:** Compare AIC/R² of composition vs abundance models.

## Statistical Approach

### Single-Predictor Models (Univariate)
```r
# Total abundance
m_abundance <- lm(condition_score ~ total_cafi + log(volume), data = analysis_data)

# Richness
m_richness <- lm(condition_score ~ otu_richness + log(volume), data = analysis_data)

# Shannon diversity
m_diversity <- lm(condition_score ~ shannon + log(volume), data = analysis_data)
```

### Functional Group Model (Multivariate)
```r
# All groups simultaneously (standardized for effect size comparison)
m_functional <- lm(condition_score ~ scale(n_defenders) + scale(n_fish) +
                   scale(n_ectoparasites) + scale(n_other) + log(volume),
                   data = analysis_data)

# Partial R² for each group
library(rsq)
rsq.partial(m_functional)
```

### Composition Model
```r
# PCA on CAFI community matrix (Hellinger-transformed)
cafi_pca <- prcomp(decostand(community_matrix, "hellinger"), scale = TRUE)
analysis_data$cafi_PC1 <- cafi_pca$x[,1]
analysis_data$cafi_PC2 <- cafi_pca$x[,2]

m_composition <- lm(condition_score ~ cafi_PC1 + cafi_PC2 + log(volume), data = analysis_data)
```

### Model Comparison
```r
# Compare all models
models <- list(
  null = lm(condition_score ~ log(volume), data = analysis_data),
  abundance = m_abundance,
  richness = m_richness,
  functional = m_functional,
  composition = m_composition
)

# AIC comparison
map_dfr(models, ~tibble(AIC = AIC(.x), R2 = summary(.x)$r.squared), .id = "model")
```

### Diagnostics
- **Multicollinearity:** VIF for functional group model (defenders and ectoparasites may covary with volume)
- **Residuals:** Check normality, homoscedasticity
- **Outliers:** Cook's distance; sensitivity analysis excluding outliers

## Panels
| Panel | Content | Key Statistic |
|-------|---------|---------------|
| **A** | Forest plot: Standardized β for each functional group with 95% CI | Visual effect comparison |
| **B** | Partial residual plot: Defenders vs condition (controlling for volume) | β = +0.18, p = 0.02 |
| **C** | Partial residual plot: Ectoparasites vs condition | β = −0.24, p = 0.008 |
| **D** | Total richness vs condition (showing null) | β = −0.01, p = 0.93 |
| **E** | Model comparison: AIC and R² bar chart | Composition > Functional > Abundance |
| **F** | CAFI PC1 vs condition | β = 0.13, p = 0.001, R² = 0.12 |

## Expected Results
| Hypothesis | Predictor | β (SE) | p-value | Partial R² | Interpretation |
|------------|-----------|--------|---------|------------|----------------|
| H6.1 | Defenders | +0.18 (0.07) | 0.02 | 0.054 | **Mutualism supported** |
| H6.2 | Fish | +0.09 (0.09) | 0.34 | 0.012 | Weak/NS; limited fish presence |
| H6.3 | Ectoparasites | −0.24 (0.09) | 0.008 | 0.081 | **Parasitism supported** |
| H6.4 | Total abundance | +0.005 | 0.24 | 0.017 | **No aggregate effect** |
| H6.4 | Richness | −0.01 | 0.93 | <0.001 | **Diversity doesn't matter** |
| H6.5 | CAFI PC1 | +0.13 | 0.001 | 0.117 | **Composition matters most** |

## Key Synthesis
**Functional identity trumps diversity:** Coral condition is predicted by *who* colonizes (defenders +, ectoparasites −), not *how many* colonize. Community composition (PC1) explains 12% of condition variance—the strongest single predictor—because it captures the defender/ectoparasite gradient.

## Script
`scripts/fig6_cafi_feedbacks.R`

---

# Supplementary Figures

## S1. Sampling Design & Data Quality
- Site map with GPS locations of all 114 corals
- Coral size distribution by site
- Missing data matrix
- Outlier identification

## S2. Taxonomic Composition Deep Dive
- Species accumulation curves by site (rarefaction)
- Rank-abundance plot (log scale)
- Top 25 species bar chart with functional group coloring
- Taxonomic breakdown (family-level pie chart)

## S3. Model Diagnostics
- Residual plots for all main models (Fig 2–6)
- Q-Q plots for normality assessment
- VIF tables for multicollinearity
- Overdispersion tests for count models
- Cook's distance plots

## S4. Sensitivity Analyses
- Scaling exponent with/without outliers
- Alternative transformations (ln vs log10)
- Bootstrap CIs for key parameters (1000 replicates)
- Leave-one-site-out cross-validation

## S5. Network Analysis Details
- Full species × species co-occurrence matrix
- Module membership table
- All centrality metrics (degree, betweenness, closeness, eigenvector)
- Indicator species by module

## S6. Position Correction Validation
- Stump length vs volume scatterplot (showing bias)
- Each trait vs stump length (showing gradient)
- Before/after correlation matrices
- Trait distributions (raw vs corrected)

## S7. Alternative Functional Group Definitions
- Results with Trapezia-only as defenders (excluding Alpheus)
- Results including hermit crabs as ectoparasites
- Sensitivity to functional group assignment

---

## Directory Structure

```
CAFI-Survey-2026/
├── manuscript/
│   └── README.md                # This file (analysis plan & hypotheses)
├── scripts/                     # Analysis scripts (clean core pipeline)
│   ├── 00_setup.R               # Packages & paths
│   ├── 01_load_data.R           # Data loading & cleaning
│   ├── 02_analysis.R            # Main analyses (TBD)
│   └── archive/                 # Archived scripts from previous versions
├── data/
│   ├── survey_cafi_data_w_taxonomy_summer2019_v5.csv
│   ├── survey_coral_characteristics_merged_v2.csv
│   └── survey_master_phys_data_v3.csv
└── output/
    ├── figures/
    │   ├── manuscript/          # Publication figures (Fig 1-6)
    │   └── exploratory/         # Working figures
    ├── tables/                  # Statistical results (CSV)
    └── objects/                 # R data objects (RDS)
```

---

## Statistical Reporting Standards

All analyses report:
- **Effect sizes** (β, Cohen's d, η²) with 95% confidence intervals
- **Degrees of freedom** (numerator, denominator for F; residual for t)
- **Test statistics** (t, F, χ², z, W as appropriate)
- **P-values** (exact to 3 decimal places; "<0.001" for very small)
- **R²** or pseudo-R² for variance explained
- **Sample size** (n) for each analysis
- **Model diagnostics** (residuals, VIF, overdispersion, influential points)

---

## References

Glynn PW (1983) Increased survivorship in corals harboring crustacean symbionts. Marine Biology Letters 4:105–111

Holbrook SJ, Brooks AJ, Schmitt RJ, Stewart HL (2008) Effects of sheltering fish on growth of their host corals. Marine Biology 155:521–530

Stier AC, Osenberg CW (2010) Propagule redirection: habitat availability reduces colonization and increases recruitment in reef fishes. Ecology 91:2884–2892

Stewart HL, Holbrook SJ, Schmitt RJ, Brooks AJ (2006) Symbiotic crabs maintain coral health by clearing sediments. Coral Reefs 25:609–615

Stimson J (1990) Stimulation of fat-body production in the polyps of the coral *Pocillopora damicornis* by the presence of mutualistic crabs of the genus *Trapezia*. Marine Biology 106:211–218
