# Supplementary Materials

## Supplementary Figures

**Figure S1.** Species accumulation curves for each reef site. Curves show rarefied species richness as a function of sampling effort (number of individuals), computed using vegan::specaccum (method = "rarefaction"). All three sites approach asymptotic richness, indicating adequate sampling coverage. Hauru supports the highest total richness, followed by Maharepa and Maatea.

**Figure S2.** PERMANOVA results are robust across five distance metrics. (A) Site and volume R² values for each metric, with significance indicated by filled versus open symbols. Site effects are significant across all five metrics; volume effects are significant for all except Raup–Crick. (B) Robustness of the composition divergence trend (distance-to-centroid vs. log volume) across metrics. Gower distance treats shared absences as similarity and should be interpreted cautiously. The site effect on composition is further validated by 500 balanced subsampling iterations across all metrics (100% of iterations significant at p < 0.05).

**Figure S3.** Composition divergence along the coral size gradient. (A) Distance-to-centroid by size class using raw Bray–Curtis dissimilarities — the apparent divergence trend is driven by larger corals having more individuals. (B) After iterated rarefaction (100 draws at equal sampling depth, averaged), the divergence trend disappears (PERMDISP: p = 0.61), confirming it is an abundance artifact.

**Figure S4.** Species-level scaling forest plot for all 21 prevalent species (≥30 individuals and ≥15% prevalence). Points show the negative binomial GLM scaling exponent (β) with 95% confidence intervals (profile CIs from GLM; bootstrap BCa CIs reported in Table S1). Species are colored by classification: blue = Redirection (β < 1, bootstrap CI excludes 1.0), gray = Field of Dreams (bootstrap CI spans 1.0), vermillion = super-linear (β > 1, CI > 1.0). Eleven of 21 species show Redirection; none shows super-linear scaling. Dashed vertical line indicates β = 1 (proportional scaling).

**Figure S5.** Neighborhood effects on CAFI abundance and diversity. Scatter plots of (A) total CAFI abundance, (B) species richness, and (C) Shannon diversity versus number of neighbors within 5 m (n = 61 corals with neighborhood data). No significant relationship is observed for abundance or richness; mean neighbor distance shows a significant negative association with richness and diversity (see main text). Power analysis indicates 65% power for medium effects (R² = 0.10) at α = 0.05.

**Figure S6.** Taxonomy sensitivity analysis. Forest plot showing seven key metrics (abundance β, richness z, PERMANOVA R² for site and volume, Shannon β, rarefied richness–condition p, network modularity Q) across five taxonomic resolution scenarios: baseline (243 OTUs), species-only (154 OTUs), merge-up, lump-down (genus-level), and rare-excluded. All scaling and composition results are robust to taxonomic resolution. Abundance β ranges from 0.515 to 0.527 across scenarios; richness z from 0.310 to 0.343; site PERMANOVA R² from 0.049 to 0.063; volume PERMANOVA R² from 0.078 to 0.098.

**Figure S7.** Taxonomic group scaling and composition. (A) Forest plot of scaling exponents for six functional groups: Trapezia crabs (β = 0.43), shrimps (β = 0.50), other crabs (β = 0.47), fish (β = 0.74), gastropods (β = 0.94, not significant), and other invertebrates (β = 0.50). Five of six groups show significant sublinear scaling; gastropods are the sole exception. (B) Group composition by site showing proportional representation.

**Figure S8.** Rarefaction depth sensitivity for the richness → condition relationship. The relationship between rarefied richness and condition was tested at rarefaction depths of n = 10, 15, 20, 25, and 30 individuals. The relationship is non-significant at all depths (all p > 0.10), confirming that the null result is not an artifact of the chosen rarefaction depth (n = 20).

**Figure S9.** Null-model co-occurrence analysis of coral-associated fauna. (A) Pairwise co-occurrence heatmap showing standardized effect sizes (SES) from a volume-weighted null model (10,000 iterations). Red = positive association (co-occurrence more frequent than expected); blue = avoidance (less frequent than expected). Species are ordered by taxonomic group. Zero of 528 pairs are significant after FDR correction, indicating that after accounting for coral volume and site, most pairwise associations are consistent with independent colonization. (B) Intraspecific density patterns for 8 focal species (of 15 tested). Red points show observed frequency of corals hosting 1, 2, 3, 4, or 5+ individuals; grey bars show null expectation (95% CI) from volume-proportional random allocation. Six species show significant mating-pair aggregation (FDR < 0.05). (C) Size-dependent co-occurrence. SES values for key species pairs across three coral size classes (terciles). Dashed line = null expectation (SES = 0); dotted lines = ±1.96 significance reference.

**Figure S10.** BEF variance partitioning and path model diagnostics. (A) Stacked bar showing decomposition of variance explained by richness and abundance (beyond volume + site): 29.1% unique to richness, <1% unique to abundance, 70.8% shared (confounded). (B) Partial regression of condition residuals on richness residuals (both after removing abundance + volume + site effects). (C) Standardized path model coefficients from a directed acyclic graph (volume → richness → condition, volume → abundance → condition); blue = significant (p < 0.05), gray = non-significant. The richness → condition path (β = 0.55) is far stronger than the abundance → condition path (β = 0.02).

**Figure S11.** Additional CAFI-condition diagnostics. (A) A priori BEF forest plot: species richness and total CAFI abundance tested as pre-specified predictors (Hochberg FWER, k = 2); richness is significant, abundance is marginal. (B) Rarefied richness (expected species at n = 20) versus coral condition — no significant relationship (p = 0.50), but this test is ambiguous (see main text). (C) Exploratory functional group forest plot (BH-FDR, k = 4): Shannon diversity, *Trapezia*, resident fish, and *Galeropsis* effects on condition; all non-significant. (D–E) Individual species scatterplots for *Trapezia* spp. and *Galeropsis monodonta* versus condition, showing no significant relationships. (F) Full bidirectional analysis across all CAFI predictors (both CAFI → condition and condition → CAFI directions).

**Figure S12.** Species occurrence probability as a function of coral size. Logistic GLM curves for 24 prevalent species showing probability of occurrence versus log(volume). Fourteen species show significant size-dependent occurrence after FDR correction (p < 0.05), with most showing increasing probability with coral size. Shaded bands show 95% confidence intervals.

**Figure S13.** Species-level condition trait heatmap. (A) Standardized regression coefficients (beta) from linear models: trait ~ sqrt(species count) + log(volume) + site, for 19 prevalent CAFI species (present on >= 5 corals) across five condition measures (condition PC1, protein, carbohydrate, zooxanthellae, AFDW; all position-corrected). Bold values with asterisks indicate raw p < 0.05. Species are organized by taxonomic group and ordered by condition PC1 beta within each group. (B) FDR-adjusted p-values (Benjamini-Hochberg across all 95 tests; grey tiles: FDR-adjusted p >= 0.10; significance threshold: FDR < 0.05). No individual species x trait test survives FDR correction, consistent with a complementarity (BEF) mechanism in which condition benefits arise from aggregate diversity rather than individual species identity. Compare to companion experiment Table S3 (Stier et al., which used Pearson correlations without volume/site controls).

**Figure S14.** Species-condition biplots for the nine strongest species x trait associations (raw p < 0.10 or |beta| > 0.40). Each panel shows the position-corrected trait value versus species abundance (sqrt-scaled x-axis) with site-colored points and linear model fit (model: trait ~ sqrt(count) + log(volume) + site). Standardized beta and raw p-value annotated. (A) *Breviturma pica* x protein (beta = +0.40, p = 0.030). (B) *Breviturma pica* x carbohydrate (+0.38, p = 0.037). (C) *Paracirrhites arcatus* x AFDW (+0.49, p = 0.043). (D) *Periclimenes* sp. x zooxanthellae (+0.36, p = 0.043). (E) *Periclimenes* sp. x condition PC1 (+0.59, p = 0.051). (F) *Harpiliopsis spinigera* x carbohydrate (-0.28, p = 0.054). (G) *Breviturma pica* x condition PC1 (+0.55, p = 0.058). (H) *Alpheus lottini* x zooxanthellae (+0.32, p = 0.060). (I) *Periclimenes* sp. x AFDW (+0.36, p = 0.062). None survives FDR correction; effect directions are consistent with mutualism hypotheses (positive for guard crabs and shrimps, variable for mobile species).

---

## Supplementary Tables

**Table S1.** Scaling analysis results for all species, taxonomic groups, and community-level metrics. For each response, the table reports the number of corals, the number with non-zero counts, the negative binomial (or Poisson) scaling exponent β, standard error, profile and bootstrap 95% CIs, the Wald z-test and bootstrap p-value against β = 1, McFadden's pseudo-R², the FDR-corrected p-value within category, and the interpretation (Redirection, Field of Dreams, or Super-linear). n = 112 corals; bootstrap = 1,000 site-stratified iterations.

**Table S2.** CAFI → condition model results. For each of seven CAFI predictors, the table reports the standardized regression coefficient (β), OLS standard error (primary), t-statistic, uncorrected p-value, corrected p-value (Hochberg for a priori BEF predictors, BH-FDR for exploratory), 95% CI, adjusted R², Breusch–Pagan p-value, and HC3-robust p-value (supplement sensitivity). Model: condition_PC1 ~ CAFI_predictor + log(volume) + site. Count predictors are square-root-transformed. n = 84 corals with physiology data. OLS is justified as primary inference by Breusch–Pagan tests confirming homoscedasticity (all BP p > 0.5).

**Table S3.** Functional group effects on condition. For Trapezia crabs, resident fish, and *Galeropsis monodonta*, the table reports the predicted direction (from ecological theory and the companion experiment), observed direction, regression coefficient, p-value, and whether the observation matches the prediction. Trapezia and fish effect directions match expectations (both positive), but *Galeropsis* shows a positive direction despite the predicted negative effect. None is significant.

**Table S4.** Cross-study species comparison: survey vs. companion experiment. For 11 species tested in the companion experimental study, the table reports experimental condition effects (β, p, Hochberg-corrected p) alongside survey condition effects (β, robust p, Hochberg-corrected p), the number of corals where the species was present in the survey, and whether effect directions match. Only *Caracanthus maculatus* was significant in the experiment after correction (p_adj = 0.017); no species is significant in the survey. Sign concordance: 2/5 testable species match direction (binomial p = 0.81).

**Table S5.** Landscape model results. Full negative binomial (abundance), Poisson (richness), and linear (Shannon diversity) model results for the 61 neighborhood-surveyed corals, including coefficients, standard errors, z/t statistics, and p-values for each predictor: log(volume), n_neighbors, total_neighbor_volume, and mean_neighbor_dist. Only log(volume) is consistently significant.

**Table S6.** Pairwise co-occurrence analysis. For all species pairs (species present on ≥10 corals), the table reports observed co-occurrence count, null model mean and SD, SES, raw p-value, FDR-corrected p-value, and direction (positive or negative). Null model: volume-weighted Bernoulli draws from logistic GLM predicted probabilities (10,000 iterations). No pair is significant after FDR correction (all p_FDR > 0.05).

**Table S7.** Intraspecific density patterns. For 15 species with sufficient abundance, the table reports total individuals, number of corals occupied, observed counts of singles/pairs/3+, null model mean and 95% CI for pair count, SES, raw and FDR-corrected p-values. Six species show significant mating-pair aggregation (FDR < 0.05), with *Synalpheus charon* (SES = 10.5) and *Alpheus lottini* (SES = 10.1) showing the strongest patterns.

**Table S8.** Cross-study scaling concordance. For 7 species with scaling data in both the survey and companion experiment, the table compares experimental scaling response (proportional or sub-proportional) with the survey scaling exponent (β), bootstrap CI, and classification (sublinear or proportional). Two species that scaled proportionally during experimental colonization show sublinear scaling in established communities (*Calcinus latens*, *Alpheus lottini*); one species was sublinear in both studies (*Trapezia serenei*); four scaled proportionally in both. The shift from proportional to sublinear is consistent with density-dependent processes emerging over time.

**Table S9.** Taxonomy sensitivity results. Seven key metrics evaluated across five taxonomic resolution scenarios. All scaling and composition metrics are robust to taxonomic resolution. Abundance β ranges from 0.515 to 0.527; richness z from 0.310 to 0.343; site PERMANOVA R² from 0.049 to 0.063.

**Table S10.** Spatial autocorrelation results. Global Moran's I for CAFI abundance, species richness, and Shannon diversity. None shows significant spatial autocorrelation (all p > 0.25).

**Table S11.** Key species scaling results (companion experiment species). Scaling exponents for species highlighted in the companion experimental study, including *Caracanthus maculatus* (β = 1.18, proportional), *Alpheus lottini* (β = 0.37, Redirection), *Trapezia serenei* (β = 0.49, Redirection), *Galeropsis monodonta* (β = 1.27, proportional), and grouped taxa (All *Alpheus* spp., harmful xanthids).

---

## Data Availability: Supplementary Table → CSV File Index

All statistical results referenced in the supplementary tables are exported as CSV files in `output/tables/`. The table below maps each supplement table to its source file(s).

| Supplement Table | CSV File(s) | Source Script |
|-----------------|-------------|---------------|
| **Table S1** (Scaling results) | `scaling_results_all.csv`, `scaling_community_level.csv`, `scaling_meta_analysis.csv` | `05_species_scaling_analysis.R` |
| **Table S2** (CAFI → condition models) | `cafi_condition_models.csv`, `breusch_pagan_diagnostics.csv`, `bef_hypothesis_corrections.csv` | `09_cafi_condition_feedbacks.R` |
| **Table S3** (Functional group effects) | `functional_effects.csv` | `09_cafi_condition_feedbacks.R` |
| **Table S4** (Cross-study species comparison) | `cross_study_species_comparison.csv`, `cross_study_sign_concordance.csv` | `09_cafi_condition_feedbacks.R` |
| **Table S5** (Landscape models) | `landscape_full_model_results.csv`, `landscape_univariate_results.csv`, `neighborhood_effects.csv` | `04_landscape_effects.R` |
| **Table S6** (Pairwise co-occurrence) | `pairwise_cooccurrence.csv` | `06_cooccurrence_analysis.R` |
| **Table S7** (Intraspecific density) | `intraspecific_density.csv` | `06_cooccurrence_analysis.R` |
| **Table S8** (Cross-study scaling concordance) | `cross_study_scaling_concordance.csv`, `key_species_scaling_experimental.csv` | `05_species_scaling_analysis.R` |
| **Table S9** (Taxonomy sensitivity) | `taxonomy_sensitivity.csv`, `taxonomy_sensitivity_species_scaling.csv` | `13_taxonomy_sensitivity.R` |
| **Table S10** (Spatial autocorrelation) | `morans_i_results.csv`, `morans_i_by_site.csv` | `07_spatial_autocorrelation.R` |
| **Table S11** (Key species scaling) | `key_species_scaling_experimental.csv`, `scaling_results_all.csv` | `05_species_scaling_analysis.R` |

### Additional CSV outputs not mapped to numbered tables

| CSV File | Description | Source Script |
|----------|-------------|---------------|
| `permanova_results.csv` | PERMANOVA Type I + Type III results | `02_community_analysis.R` |
| `permdisp_results.csv` | PERMDISP homogeneity of dispersions | `02_community_analysis.R` |
| `permanova_pairwise.csv` | Pairwise site comparisons (Bonferroni) | `02_community_analysis.R` |
| `permanova_metric_sensitivity.csv` | 5-metric PERMANOVA robustness | `02_community_analysis.R` |
| `permanova_subsampling_summary.csv` | 500-iteration balanced subsampling | `02_community_analysis.R` |
| `dbrda_variance_partitioning.csv` | db-RDA variance partitioning (volume/site) | `02_community_analysis.R` |
| `dbrda_species_scores.csv` | Top 15 species loadings on constrained axis | `02_community_analysis.R` |
| `dbrda_rarefied.csv` | Rarefied db-RDA robustness check | `02_community_analysis.R` |
| `nestedness_nodf.csv` | NODF nestedness test (z-score, p-value) | `02_community_analysis.R` |
| `composition_divergence_stats.csv` | PERMDISP by size class (raw + rarefied) | `02_community_analysis.R` |
| `density_dilution.csv` | Per-capita density slope vs. volume | `05_species_scaling_analysis.R` |
| `occurrence_scaling_results.csv` | Species occurrence logistic GLM results | `05_species_scaling_analysis.R` |
| `occurrence_summary.csv` | FDR summary: n significant/positive/negative | `05_species_scaling_analysis.R` |
| `scaling_interpretation_summary.csv` | Species scaling classification summary | `05_species_scaling_analysis.R` |
| `bef_diversity_abundance_partition.csv` | Variance partitioning: richness vs. abundance | `09_cafi_condition_feedbacks.R` |
| `bef_path_coefficients.csv` | Path model (piecewiseSEM) coefficients | `09_cafi_condition_feedbacks.R` |
| `power_analysis.csv` | Power analysis for Q3 and Q4 | `09_cafi_condition_feedbacks.R` |
| `richness_abundance_correlation.csv` | Richness–abundance Pearson r | `09_cafi_condition_feedbacks.R` |
| `rarefaction_depth_sensitivity.csv` | Richness → condition across rarefaction depths | `09_cafi_condition_feedbacks.R` |
| `key_species_effects.csv` | Individual species → condition effects | `09_cafi_condition_feedbacks.R` |
| `species_trait_correlations.csv` | Species × trait Pearson r matrix | `09_cafi_condition_feedbacks.R` |
| `species_trait_heatmap_data.csv` | Species × trait β values + FDR p-values | `09_cafi_condition_feedbacks.R` |
| `reverse_direction_models.csv` | Condition → CAFI reverse-direction tests | `09_cafi_condition_feedbacks.R` |
| `size_dependent_cooccurrence.csv` | Co-occurrence SES by size class | `06_cooccurrence_analysis.R` |
| `network_metrics.csv` | Network topology metrics | `06_cooccurrence_analysis.R` |
| `hub_species.csv` | Hub species centrality scores | `06_cooccurrence_analysis.R` |
| `taxonomic_group_scaling.csv` | Functional group scaling exponents | `08_functional_groups.R` |

---

## Supplementary Methods

### Position correction for physiological traits

Coral branch tissue exhibits gradients in protein, carbohydrate, and zooxanthellae density along the branch axis, with higher values near the base and lower values at the tips. Because smaller colonies required sampling a proportionally larger fraction of the branch (including lower-density tips), a systematic size bias exists in raw trait values. To correct for this, we regressed each physiological trait on stump length (distance from branch base to the cut) and nubbin length (total fragment length) using linear models. The residuals from these regressions represent position-corrected trait values that can be meaningfully compared across colony sizes. PCA was then performed on the four corrected traits (protein, carbohydrate, zooxanthellae density, AFDW) to derive condition PC1.

### Volume-weighted co-occurrence null model

The null model for pairwise co-occurrence follows the framework of Stier et al. (2012), extended to account for the coral size confound. For each of the species included (present on ≥10 corals), we fit a logistic GLM:

    P(species present | coral i) ~ log(volume_i) + site_i

yielding a predicted probability p_ij for species j on coral i. For each null iteration (N = 10,000), species presences were drawn independently from Bernoulli(p_ij) distributions, preserving the volume-dependent occurrence structure while randomizing pairwise associations. For each species pair, observed co-occurrence was compared to the null distribution to compute SES and a two-sided p-value. FDR correction (Benjamini–Hochberg) was applied across all pairwise tests. This approach ensures that any detected co-occurrence signal is not a trivial consequence of both species preferring large (or small) corals.

### Intraspecific density null model

To test whether conspecifics aggregate on individual corals more than expected (consistent with mating-pair formation), we used a multinomial allocation null model. For each species, N total individuals were allocated across corals proportional to coral volume (larger corals receive more individuals). For each of 10,000 null iterations, we recorded the number of corals with exactly two individuals ("pairs") and compared the observed count to the null distribution using SES. Species with SES > 1.96 (FDR-corrected) were classified as showing significant mating-pair aggregation.

### Power analysis

Post hoc power analyses were conducted using the pwr package (Champely 2020). For Q3 (CAFI–condition feedbacks; n = 84), we computed power to detect medium (R² = 0.10; power = 76%) and small (R² = 0.05; power = 43%) effects at α = 0.05 with 2 numerator df (CAFI predictor) and site + volume as covariates. The minimum detectable R² at 80% power was 0.108. To assess cross-study comparability, we computed power to detect the companion experiment's observed effect size (R² ≈ 0.12): the survey had 85% power, confirming that the null result is not an artifact of insufficient sample size. For Q4 (neighborhood effects; n = 61, 3 numerator df), power was 65% for medium effects and 35% for small effects.

### Distance metric sensitivity

To ensure that compositional results were not driven by distance metric choice, we repeated the marginal PERMANOVA (log(volume) + site) using five metrics: (1) Bray–Curtis on Hellinger-transformed abundances (primary analysis); (2) Bray–Curtis on Wisconsin-transformed abundances; (3) Jaccard on presence–absence data; (4) Gower distance on raw abundances; and (5) Raup–Crick (null-model-based). Site effects were significant across all five metrics (all p ≤ 0.001). Volume effects were significant for all metrics except Raup–Crick. Composition divergence by size class was not robust across metrics (significant for Gower only), consistent with the rarefaction result showing the divergence pattern is an abundance artifact.

### Balanced site subsampling

To verify that PERMANOVA results are robust to unequal site sample sizes (35–39 corals per site), we performed 500 balanced subsampling iterations. In each iteration, n = min(site counts) = 30 corals were randomly drawn from each site (after excluding low-abundance corals), and PERMANOVA was re-run. The site effect was significant (p < 0.05) in 100% of iterations across all metrics, confirming that the compositional difference is not an artifact of sample-size imbalance.
