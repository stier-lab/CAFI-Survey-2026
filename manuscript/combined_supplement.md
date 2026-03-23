---
title: "Supplementary Materials: Colony size drives sublinear scaling, compositional turnover, and biodiversity--condition feedbacks in coral-associated fauna"
bibliography: references.bib
csl: journal-of-animal-ecology.csl
link-citations: true
lang: en
---

# Supplementary Materials

**Colony size drives sublinear scaling, compositional turnover, and biodiversity--condition feedbacks in coral-associated fauna**

Stier, A.C. et al.

Abbreviations: CAFI, coral-associated fishes and invertebrates; BEF, biodiversity--ecosystem function; SES, standardized effect size; FDR, false discovery rate (Benjamini--Hochberg); PERMANOVA, permutational multivariate analysis of variance; PERMDISP, permutational analysis of multivariate dispersions; db-RDA, distance-based redundancy analysis; NODF, nestedness metric based on overlap and decreasing fill; OTU, operational taxonomic unit; GLM, generalized linear model; OLS, ordinary least squares; PCA, principal component analysis; AFDW, ash-free dry weight; SAR, species--area relationship.

---

## Supplementary Figures

**Figure S1.** Species accumulation curves for each reef site. Curves show rarefied species richness as a function of sampling effort (number of individuals), computed using vegan::specaccum (method = "rarefaction"). All three sites approach asymptotic richness, indicating adequate sampling coverage. Hauru supports the highest total richness, followed by Maharepa and Maatea.

![](../output/figures/supplement/figS1_species_accumulation.png)

---

**Figure S2.** PERMANOVA results are robust across five distance metrics. (A) Site and volume R² values for each metric, with significance indicated by filled versus open symbols. Site effects are significant across all five metrics; volume effects are significant for all except Gower (which treats shared absences as similarity and should be interpreted cautiously). (B) Robustness of the composition divergence trend (distance-to-centroid vs. log volume) across metrics. The site effect on composition is further validated by 500 balanced subsampling iterations across all metrics (100% of iterations significant at p < 0.05).

![](../output/figures/supplement/figS2_permanova_sensitivity.png)

---

**Figure S3.** Composition divergence along the coral size gradient. (A) Distance-to-centroid by size class using raw Bray--Curtis dissimilarities -- the apparent divergence trend is driven by larger corals having more individuals. (B) After iterated rarefaction (100 draws at equal sampling depth, averaged), the divergence trend disappears (PERMDISP: p = 0.61), confirming it is an abundance artifact. (C) Continuous convergence trend: distance-to-centroid vs. ln(colony volume), colored by site with OLS regression line (β = −0.034, R² = 0.26, p < 0.001). This strong raw trend also vanishes after rarefaction (β = 0.002, p = 0.65), confirming that the apparent compositional convergence of larger corals reflects richer samples rather than deterministic homogenization.

![](../output/figures/supplement/figS3_composition_divergence.png)

---

**Figure S4.** Species-level scaling forest plot for all 21 prevalent species (≥30 individuals and ≥15% prevalence). Points show the negative binomial GLM scaling exponent (β) with 95% confidence intervals (profile likelihood CIs from GLM; bootstrap BCa CIs reported in Table S1). Species are colored by scaling classification: blue = Redirection (β < 1, bootstrap CI excludes 1.0), gray = Field of Dreams (bootstrap CI spans 1.0), vermillion = super-linear (β > 1, CI > 1.0). Eleven of 21 species show Redirection; none shows super-linear scaling. Dashed vertical line at β = 1 indicates proportional scaling.

![](../output/figures/supplement/figS4_species_scaling.png)

---

**Figure S5.** Neighborhood effects on CAFI communities. Three response variables (rows: total CAFI abundance, species richness, Shannon diversity) plotted against two neighborhood predictors (columns: number of neighbors within 5 m, mean distance to neighbors). (A) Abundance vs. neighbor count; (B) abundance vs. neighbor distance; (C) richness vs. neighbor count; (D) richness vs. neighbor distance; (E) Shannon diversity vs. neighbor count; (F) Shannon diversity vs. neighbor distance. Trend lines shown only for significant relationships (p < 0.05): mean neighbor distance negatively predicts richness (D; β = −0.005, p = 0.001) and Shannon diversity (F; β = −0.007, p = 0.001), but not abundance (B; p = 0.78). Neighbor count does not predict any response (A, C, E; all p > 0.37). n = 61 corals with 5-m neighborhood surveys across three reef sites. Points colored by site (purple = Hauru, slate = Maatea, green = Maharepa). Power analysis indicates 65% power for medium effects (R² = 0.10) at α = 0.05.

![](../output/figures/supplement/figS5_neighborhood.png)

---

**Figure S6.** Taxonomy sensitivity analysis. Forest plot showing seven key metrics (abundance β, richness z, PERMANOVA R² for site and volume, Shannon β, rarefied richness--condition p, network modularity Q) across five taxonomic resolution scenarios: baseline (243 OTUs), species-only (154 OTUs), merge-up, lump-down (genus-level), and rare-excluded. All scaling and composition results are robust to taxonomic resolution. Abundance β ranges from 0.515 to 0.527 across scenarios; richness z from 0.310 to 0.343; site PERMANOVA R² from 0.049 to 0.063; volume PERMANOVA R² from 0.078 to 0.098.

![](../output/figures/supplement/figS6_taxonomy_sensitivity.png)

---

**Figure S7.** *Archived — available upon request. Taxonomic group scaling and composition (partially redundant with Fig. 3D main text).*

---

**Figure S8.** Rarefaction depth sensitivity for the richness → condition relationship. The relationship between rarefied richness and condition was tested at rarefaction depths of n = 10, 15, 20, 25, and 30 individuals. The relationship is non-significant at all depths (all p > 0.10), confirming that the null result is not an artifact of the chosen rarefaction depth (n = 20).

![](../output/figures/supplement/figS8_rarefaction_sensitivity.png)

---

**Figure S9.** Null-model co-occurrence analysis of coral-associated fauna. (A) Pairwise co-occurrence heatmap showing standardized effect sizes (SES) from a volume-weighted null model (10,000 iterations). Red = positive association (co-occurrence more frequent than expected); blue = avoidance (less frequent than expected). Species are ordered by taxonomic group. No pairs are significant after FDR correction (0 of 528 at FDR < 0.05), indicating that after accounting for coral volume and site, pairwise associations are consistent with independent colonization. The strongest signals are *Harpiliopsis beaupresii*--*Paragobiodon modestus* (SES = −3.43, p_FDR = 0.11) and *H. beaupresii*--*H. spinigera* (SES = −3.42, p_FDR = 0.21). (B) Intraspecific density patterns for 8 focal species (of 15 tested). Red points show observed frequency of corals hosting 1, 2, 3, 4, or 5+ individuals; gray bars show null expectation (95% CI) from volume-proportional random allocation. Six species show significant mating-pair aggregation (FDR < 0.05). (C) Size-dependent co-occurrence. SES values for key species pairs across three coral size classes (terciles). Dashed line = null expectation (SES = 0); dotted lines = ±1.96 significance reference.

![](../output/figures/supplement/figS9_cooccurrence.png)

---

**Figure S10.** BEF variance partitioning and path model diagnostics. (A) Stacked bar chart showing decomposition of variance explained by richness and abundance (beyond volume + site): 29.1% unique to richness, <1% unique to abundance, 70.8% shared (confounded). (B) Partial regression of condition residuals on richness residuals (both after removing abundance + volume + site effects). (C) Standardized path model coefficients from piecewiseSEM (volume → richness → condition, volume → abundance → condition). The richness → condition path (β = 0.55) is substantially larger in magnitude than the abundance → condition path (β = 0.02), though neither reached significance individually (see main text).

![](../output/figures/supplement/figS10_bef_variance_partitioning.png)

---

**Figure S11.** CAFI--condition diagnostics. (A) A priori BEF forest plot: species richness and total CAFI abundance tested as pre-specified predictors; richness is significant (p = 0.018), abundance is marginally significant (p = 0.048). (B) Rarefied richness (expected species at n = 20) versus coral condition -- no significant relationship (p = 0.50), but this test is ambiguous (see main text). (C) Full bidirectional analysis across all CAFI predictors (both CAFI → condition and condition → CAFI directions). *Panels showing exploratory functional group predictors and individual species scatterplots archived — available upon request.*

![](../output/figures/supplement/figS11_condition_details.png)

---

**Figures S12--S14.** *Archived — available upon request. S12: Species occurrence probability curves (24 species; logistic GLM); S13: Species × condition trait heatmap (19 species × 5 traits; 0 of 95 tests survived FDR); S14: Species × condition biplots (9 strongest associations). These analyses are fully computed in `05_species_scaling_analysis.R` and `09_cafi_condition_feedbacks.R`.*

---

## Supplementary Tables

**Table S1.** Scaling analysis results for all species, taxonomic groups, and community-level metrics. For each response, the table reports the number of corals with volume data, the number with non-zero counts, the negative binomial (or Poisson) scaling exponent β, standard error, profile likelihood and bootstrap 95% CIs, the Wald z-test and bootstrap p-value against β = 1, deviance-explained R², the FDR-corrected p-value within category, and the scaling interpretation (Redirection, Field of Dreams, or Super-linear). n = 112 corals with volume data (of 114 surveyed); bootstrap = 1,000 site-stratified iterations.

*Data: `output/tables/scaling_results_all.csv`, `scaling_community_level.csv`, `scaling_meta_analysis.csv`*

**Table S2.** CAFI--condition model results. For each of seven CAFI predictors (species richness, total CAFI abundance, Shannon diversity, *Trapezia* abundance, resident fish abundance, *Galeropsis monodonta* abundance, and community composition PC1), the table reports the standardized regression coefficient (β), OLS standard error (primary), t-statistic, uncorrected p-value, FDR-corrected p-value (Hochberg step-up for the two a priori predictors; BH-FDR for the four exploratory predictors; uncorrected for the single PC1 test), 95% CI, adjusted R², Breusch--Pagan p-value (homoscedasticity diagnostic), and HC3-robust p-value (sensitivity check). Model: condition_PC1 ~ CAFI_predictor + log(volume) + site. Count predictors are square-root-transformed. n = 84 corals with complete physiology data.

*Data: `output/tables/cafi_condition_models.csv`, `breusch_pagan_diagnostics.csv`, `bef_hypothesis_corrections.csv`*

**Table S3.** Cross-study condition concordance. (A) Functional group effects: for *Trapezia* crabs, resident fish, and *Galeropsis monodonta*, the table reports the predicted direction (from ecological theory and the companion experiment), observed direction, regression coefficient, p-value, and whether the observation matches the prediction. (B) Species-level comparison: for 11 species tested in the companion experimental study, the table reports experimental condition effects (β, p) alongside survey condition effects (β, HC3-robust p-value), the number of corals where the species was present in the survey, and whether effect directions match between studies.

*Data: `output/tables/functional_effects.csv`, `cross_study_species_comparison.csv`, `cross_study_sign_concordance.csv`*

**Table S4.** Neighborhood model results. Full negative binomial (abundance), Poisson (richness), and linear (Shannon diversity) model results for the 61 neighborhood-surveyed corals. Predictors: log(volume), number of neighbors within 5 m, log(total neighbor volume + 1), mean inter-colony distance, and site. Table reports coefficients, standard errors, z/t statistics, and p-values for each predictor.

*Data: `output/tables/landscape_full_model_results.csv`, `landscape_univariate_results.csv`, `neighborhood_effects.csv`*

**Table S5.** Pairwise co-occurrence analysis. For all species pairs (species present on ≥10 corals), the table reports observed co-occurrence count, null model mean and SD, SES, raw p-value, FDR-corrected p-value, and direction. Null model: volume-weighted Bernoulli draws from logistic GLM predicted probabilities (10,000 iterations).

*Data: `output/tables/pairwise_cooccurrence.csv`*

**Table S6.** Intraspecific density patterns. For 15 species with sufficient abundance for the multinomial null model (see Supplementary Methods), the table reports total individuals, number of corals occupied, observed counts of singles/pairs/3+, null model mean and 95% CI for pair count, SES, raw and FDR-corrected p-values.

*Data: `output/tables/intraspecific_density.csv`*

**Table S7.** Cross-study scaling concordance. For 7 species with scaling data in both the survey and companion experiment, the table reports the survey scaling exponent (β), standard error, bootstrap CI, classification, and the corresponding experimental scaling response. This table consolidates the concordance analysis and the key species scaling details into a single reference.

*Data: `output/tables/cross_study_scaling_concordance.csv`, `key_species_scaling_experimental.csv`, `scaling_results_all.csv`*

**Table S8.** Taxonomy sensitivity results. Seven key metrics (abundance β, richness z, PERMANOVA R² for site and volume, Shannon diversity slope, rarefied richness--condition p-value, and network modularity Q) evaluated across five taxonomic resolution scenarios (baseline, species-only, merge-up, lump-down, rare-excluded; see Supplementary Methods for scenario definitions). All scaling and composition metrics are robust to taxonomic resolution.

*Data: `output/tables/taxonomy_sensitivity.csv`, `taxonomy_sensitivity_species_scaling.csv`*

**Table S9.** Spatial autocorrelation results. Global Moran's I for CAFI abundance, species richness, and Shannon diversity across all corals and within each site. None shows significant spatial autocorrelation (all p > 0.28).

*Data: `output/tables/morans_i_results.csv`, `morans_i_by_site.csv`*

**Table S10.** BEF model diagnostics. (A) Variance partitioning: hierarchical R² decomposition showing variance uniquely attributable to species richness, uniquely to total abundance, and shared between them, for the condition ~ richness + abundance + log(volume) + site model (n = 84). (B) Path model coefficients: standardized path coefficients from the piecewiseSEM model: volume → richness → condition and volume → abundance → condition (n = 84).

*Data: `output/tables/bef_diversity_abundance_partition.csv`, `bef_path_coefficients.csv`*

**Table S11.** Beta diversity partitioning: turnover vs nestedness components (betapart decomposition).

*Data: `output/tables/sensitivity_betapart_decomposition.csv`*

**Table S12.** Mediation analysis: richness → abundance → condition (bootstrap ACME).

*Data: `output/tables/sensitivity_mediation_richness_abundance.csv`*

**Table S13.** Morphotype covariate sensitivity: BEF model with and without coral morphotype.

*Data: `output/tables/sensitivity_morphotype_bef.csv`*

**Table S14.** Missing data characterization: logistic regression of dropout predictors.

*Data: `output/tables/sensitivity_missing_data_predictors.csv`*

**Table S15.** Haplotype sensitivity analyses. Molecular species assignments (mtORF haplotyping per @JohnstonCunningBurgess2022) were available for 101 of 114 colonies (99 of the 112 with volume data) and used to test whether genetic species identity confounds the main results. (A) **Morphotype--haplotype concordance.** Cross-tabulation of field morphotype assignment vs. genetic species. Concordance varied among morphotypes: "eudoxi" mapped to *P. grandis* in 94% of cases (45/48), "meandrina" to *P. meandrina* in 66% (19/29), and "verucosa" to *P. verrucosa* in only 33% (7/21). Chi-square test with simulated p-value confirmed significant non-random association between morphotype and haplotype. (B) **BEF with haplotype covariate.** Adding genetic species identity as a covariate to the richness → condition model reduced the richness coefficient by only 8% (β = 0.081, p = 0.008 without haplotype; β = 0.074, p = 0.010 with haplotype; n = 74 with both physiology and haplotype data). The richness signal is robust to genetic species identity, in contrast to morphotype which absorbed the signal entirely (Table S13: p went from 0.015 to 0.27). (C) **Scaling by genetic species.** The volume × species interaction was not significant (p = 0.39), indicating that scaling exponents do not differ among species. Species-specific estimates were all sublinear: *P. grandis* β = 0.46 (n = 49), *P. meandrina* β = 0.62 (n = 32), *P. verrucosa* β = 0.73 (n = 10, not significant, p = 0.17). (D) **Composition by genetic species.** Marginal PERMANOVA (community ~ log(volume) + genetic species + site) confirmed that genetic species explains 9.8% of composition variation (F = 2.95, p = 0.001), comparable to volume (8.1%, p = 0.001) and exceeding site (6.0%, p = 0.001). (E) **Phylogenetic distance and symbiont identity.** Mantel test detected a significant correlation between phylogenetic distance and CAFI community dissimilarity (r = 0.12, p = 0.013), which strengthened after controlling for volume (partial Mantel r = 0.13, p = 0.002). Symbiont identity (*Cladocopium latusorum* vs. *C. pacificum*, assigned per @Burgess2021) explained 3.8% of composition (PERMANOVA R² = 0.038, p = 0.001) and predicted coral condition (β = 1.28, p = 0.008; n = 74).

*Data: `output/tables/sensitivity_morphotype_haplotype_concordance.csv`, `sensitivity_haplotype_bef.csv`, `sensitivity_scaling_by_species.csv`, `sensitivity_composition_by_species.csv`, `sensitivity_phylogenetic_symbiont.csv`*

---

## Supplementary Methods

### Iterated rarefaction for composition analyses

To control for abundance confounds in composition analyses, we used iterated rarefaction. For each of 100 iterations, abundances were rarefied to the minimum observed abundance per coral (drawing without replacement from each colony's species-abundance vector), and the resulting rarefied community matrix was used to recompute Bray--Curtis dissimilarities. Ordination, PERMANOVA, PERMDISP, and db-RDA were then repeated on the averaged distance matrix. This approach preserves the rank-order of species dominance while equalizing total abundance across colonies, distinguishing genuine compositional turnover from passive sampling artifacts. The rarefied db-RDA confirmed that volume explains significant composition variation (2.4% of variation, p = 0.001) even after removing the abundance confound.

### Distance metric sensitivity

To ensure that compositional results were not driven by distance metric choice, we repeated the marginal PERMANOVA (log(volume) + site) using five metrics: (1) Bray--Curtis on Hellinger-transformed abundances (primary analysis); (2) Bray--Curtis on Wisconsin-transformed abundances; (3) Jaccard on presence--absence data; (4) Gower distance on raw abundances; and (5) Raup--Crick (null-model-based). Site effects were significant across all five metrics (all p ≤ 0.001). Volume effects were significant for all metrics except Gower (Volume p = 1.0; Gower treats shared absences as similarity, inflating distances among sparse communities). Composition divergence by size class was not robust across metrics, consistent with the rarefaction result showing the divergence pattern is an abundance artifact.

### Balanced site subsampling

To verify that PERMANOVA results are robust to unequal site sample sizes (35--39 corals per site), we performed 500 balanced subsampling iterations. In each iteration, n = min(site counts) = 30 corals were randomly drawn from each site (after excluding low-abundance corals), and PERMANOVA was re-run. The site effect was significant (p < 0.05) in 100% of iterations (mean p = 0.001 for Bray--Curtis on Hellinger-transformed abundances), confirming that the site signal is robust to sample size reduction. This result was consistent across all five distance metrics (Fig. S2).

### Volume-weighted co-occurrence null model

The null model for pairwise co-occurrence follows the framework of @Stier2012, extended to account for the coral size confound. For each of the species included (present on ≥10 corals), we fit a logistic GLM:

P(species present | coral *i*) ~ log(volume~*i*~) + site~*i*~

yielding a predicted probability p~ij~ for species *j* on coral *i*. For each null iteration (N = 10,000), species presences were drawn independently from Bernoulli(p~ij~) distributions, preserving the volume-dependent occurrence structure while randomizing pairwise associations. For each species pair, observed co-occurrence was compared to the null distribution to compute SES and a two-sided p-value. FDR correction (Benjamini--Hochberg) was applied across all pairwise tests. This approach ensures that any detected co-occurrence signal is not a trivial consequence of both species preferring large (or small) corals.

### Intraspecific density null model

To test whether conspecifics aggregate on individual corals more than expected (consistent with mating-pair formation), we used a multinomial allocation null model. For each species, N total individuals were allocated across corals proportional to coral volume (larger corals receive more individuals). For each of 10,000 null iterations, we recorded the number of corals with exactly two individuals ("pairs") and compared the observed count to the null distribution using SES. Species with SES > 1.96 (FDR-corrected) were classified as showing significant mating-pair aggregation.

### Position correction for physiological traits

Coral branch tissue exhibits gradients in protein, carbohydrate, and zooxanthellae density along the branch axis, with higher values near the base and lower values at the tips. Because smaller colonies required sampling a proportionally larger fraction of the branch (including lower-density tips), a systematic size bias exists in raw trait values. To correct for this, we regressed each physiological trait on stump length (distance from branch base to the cut) and nubbin length (total fragment length) using linear models. The residuals from these regressions were standardized (z-scored) to place all traits on a common scale, yielding position-corrected trait values that can be meaningfully compared across colony sizes. PCA was then performed on the four standardized corrected traits (protein, carbohydrate, zooxanthellae density, AFDW) to derive condition PC1. A sensitivity analysis comparing corrected and raw (uncorrected) condition scores confirmed that the position correction did not qualitatively change the richness--condition relationship (both yielded the same significance outcome).

### Power analysis

Prospective power analyses were conducted using the pwr package [@Champely2020]. For the CAFI--condition feedbacks (Q3; n = 84), we computed power to detect medium (R² = 0.10; power = 76%) and small (R² = 0.05; power = 43%) effects at α = 0.05 with 2 numerator degrees of freedom (CAFI predictor + intercept) and site + volume as covariates. The minimum detectable R² at 80% power was 0.108. To assess cross-study comparability, we computed power to detect the companion experiment's observed effect size (R² ≈ 0.12): the survey had 85% power, confirming that the null result for individual predictors is not an artifact of insufficient sample size. For the neighborhood analysis (n = 61, 3 numerator degrees of freedom), power was 65% for medium effects (R² = 0.10) and 35% for small effects (R² = 0.05).

### Taxonomy sensitivity analysis

To assess whether results depended on taxonomic resolution, we repeated seven key analyses under five resolution scenarios: (1) baseline (243 OTUs at the finest available resolution), (2) species-only (154 species-level OTUs; higher-level identifications excluded), (3) merge-up (specimens identified to genus or family merged with congeneric or confamilial species-level OTUs), (4) lump-down (all OTUs lumped to genus level, reducing resolution), and (5) rare-excluded (OTUs with fewer than 3 total individuals removed). For each scenario, we recomputed: abundance scaling exponent (β), richness scaling exponent (z), PERMANOVA R² for site and volume, Shannon diversity scaling slope, rarefied richness--condition p-value, and network modularity (Q). All scenarios produced qualitatively identical results (Fig. S6; Table S8).

### Sensitivity analyses

Four supplementary sensitivity analyses addressed methodological gaps identified during pre-submission audit (Tables S11--S14; script `14_supplementary_sensitivity.R`). A fifth set of haplotype sensitivity analyses (Table S15) is described separately below. (1) **Beta diversity partitioning** (betapart; Baselga 2010) decomposed total Sørensen dissimilarity into turnover (Simpson) and nestedness-resultant components along the coral size gradient (Table S11). (2) **Mediation analysis** (mediation package; 1,000 bootstrap iterations) tested whether the richness → condition effect operates through the abundance pathway (ACME; Table S12). (3) **Morphotype covariate** sensitivity tested whether adding coral morphotype to the BEF model changed the richness → condition result (Table S13). (4) **Missing data characterization** used logistic regression to test whether the 114 → 84 dropout for condition analyses was predicted by volume, site, richness, or abundance (Table S14). Additional sensitivity analyses (envfit species vectors, indicator species by size class, iNEXT coverage-based rarefaction, community-weighted mean obligate-mutualist fraction, nonlinear BEF, C-score community co-occurrence, AIC model averaging) were computed and are available upon request.

### Haplotype sensitivity analyses

Molecular species assignments were obtained via mtORF haplotyping following @JohnstonCunningBurgess2022. PCR amplification and sequencing of the mitochondrial open reading frame successfully genotyped 101 of 114 colonies (11 PCR failures, 2 without tissue samples); 99 of the 112 colonies with volume data had valid haplotypes. Genotyped colonies resolved to four genetic species: *Pocillopora grandis* (n = 49), *P. meandrina* (n = 32), *P. tuahiniensis* (n = 7), and *P. verrucosa* (n = 10), with 1 colony unresolved. Phylogenetic distances among host species were derived from the *Pocillopora* mtORF gene tree. Symbiont identity (*Cladocopium latusorum* vs. *C. pacificum*) was assigned based on the host species--symbiont associations documented by @Burgess2021, which showed strong species-specificity in Mo'orean *Pocillopora*. Five analyses tested whether genetic species identity confounds the main results (Table S15): (A) morphotype--haplotype concordance via chi-square test with simulated p-value; (B) BEF model with genetic species as a covariate; (C) volume × species interaction test for scaling differences among genetic species; (D) marginal PERMANOVA with genetic species as a predictor of community composition; and (E) Mantel and partial Mantel tests for phylogenetic distance--community dissimilarity correlations, plus PERMANOVA and linear model tests for symbiont effects on composition and condition.

---

## Data Availability: Supplementary Table → CSV File Index

All statistical results referenced in the supplementary tables are exported as CSV files in `output/tables/`. The table below maps each supplement table to its source file(s).

| Supplement Table | CSV File(s) | Source Script |
|-----------------|-------------|---------------|
| **Table S1** (Scaling results) | `scaling_results_all.csv`, `scaling_community_level.csv`, `scaling_meta_analysis.csv` | `05_species_scaling_analysis.R` |
| **Table S2** (CAFI → condition models) | `cafi_condition_models.csv`, `breusch_pagan_diagnostics.csv`, `bef_hypothesis_corrections.csv` | `09_cafi_condition_feedbacks.R` |
| **Table S3** (Cross-study condition concordance) | `functional_effects.csv`, `cross_study_species_comparison.csv`, `cross_study_sign_concordance.csv` | `09_cafi_condition_feedbacks.R` |
| **Table S4** (Neighborhood models) | `landscape_full_model_results.csv`, `landscape_univariate_results.csv`, `neighborhood_effects.csv` | `04_landscape_effects.R` |
| **Table S5** (Pairwise co-occurrence) | `pairwise_cooccurrence.csv` | `06_cooccurrence_analysis.R` |
| **Table S6** (Intraspecific density) | `intraspecific_density.csv` | `06_cooccurrence_analysis.R` |
| **Table S7** (Cross-study scaling concordance) | `cross_study_scaling_concordance.csv`, `key_species_scaling_experimental.csv`, `scaling_results_all.csv` | `05_species_scaling_analysis.R` |
| **Table S8** (Taxonomy sensitivity) | `taxonomy_sensitivity.csv`, `taxonomy_sensitivity_species_scaling.csv` | `13_taxonomy_sensitivity.R` |
| **Table S9** (Spatial autocorrelation) | `morans_i_results.csv`, `morans_i_by_site.csv` | `07_spatial_autocorrelation.R` |
| **Table S10** (BEF model diagnostics) | `bef_diversity_abundance_partition.csv`, `bef_path_coefficients.csv` | `09_cafi_condition_feedbacks.R` |
| **Table S11** (Beta diversity partitioning) | `sensitivity_betapart_decomposition.csv` | `14_supplementary_sensitivity.R` |
| **Table S12** (Mediation: richness → abundance → condition) | `sensitivity_mediation_richness_abundance.csv` | `14_supplementary_sensitivity.R` |
| **Table S13** (Morphotype covariate sensitivity) | `sensitivity_morphotype_bef.csv` | `14_supplementary_sensitivity.R` |
| **Table S14** (Missing data predictors) | `sensitivity_missing_data_predictors.csv` | `14_supplementary_sensitivity.R` |
| **Table S15** (Haplotype sensitivity) | `sensitivity_morphotype_haplotype_concordance.csv`, `sensitivity_haplotype_bef.csv`, `sensitivity_scaling_by_species.csv`, `sensitivity_composition_by_species.csv`, `sensitivity_phylogenetic_symbiont.csv` | `14_supplementary_sensitivity.R` |

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
| `power_analysis.csv` | Power analysis for condition and neighborhood models | `09_cafi_condition_feedbacks.R` |
| `richness_abundance_correlation.csv` | Richness--abundance Pearson r | `09_cafi_condition_feedbacks.R` |
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

## References

::: {#refs}
:::
