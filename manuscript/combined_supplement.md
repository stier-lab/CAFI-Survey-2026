---
title: "Supplementary Materials: Colony size drives density dilution, species turnover, and biodiversity--condition covariation in coral-associated fauna"
bibliography: references.bib
csl: journal-of-animal-ecology.csl
link-citations: true
lang: en
---

# Supplementary Materials

**Colony size drives density dilution, species turnover, and biodiversity--condition covariation in coral-associated fauna**

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

**Figure S5.** Neighborhood effects on CAFI communities. Three response variables (rows: total CAFI abundance, species richness, Shannon diversity) plotted against two neighborhood predictors (columns: number of neighbors within 5 m, mean distance to neighbors). (A) Abundance vs. neighbor count; (B) abundance vs. neighbor distance; (C) richness vs. neighbor count; (D) richness vs. neighbor distance; (E) Shannon diversity vs. neighbor count; (F) Shannon diversity vs. neighbor distance. Trend lines shown only for significant relationships (p < 0.05): mean neighbor distance negatively predicts richness (D; β = −0.005, p = 0.001) and Shannon diversity (F; β = −0.007, p = 0.001), but not abundance (B; p = 0.78). Neighbor count does not predict any response (A, C, E; all p > 0.37). n = 61 corals with 5-m neighborhood surveys across three reef sites. Points colored by site (purple = Hauru, slate = Maatea, green = Maharepa). Power-analysis details for this reduced neighborhood subset are provided in the Supplementary Methods.

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

**Figure S12.** Host species identity and CAFI composition. (A) Distance-based redundancy analysis (db-RDA) biplot constrained by coral genetic species (*Pocillopora grandis*, *P. meandrina*, *P. tuahiniensis*, *P. verrucosa*) after partialing out colony volume and site effects; coral species explains 8.3% of compositional variation (F = 3.24, p = 0.001; n = 98 genotyped corals). Point symbols indicate coral species; vectors show top CAFI species loadings on constrained axes. (B) Variance partitioning (Hellinger-transformed community matrix) among three predictors: coral species (5.4% unique), log(volume) (4.9% unique), and site (4.1% unique), with minimal shared fractions (<0.2%), indicating largely independent axes of community variation. Total explained: 14.4%.

![](../output/figures/supplement/figS12_genotype_composition.png)

---

**Figure S13.** Systematic screen of CAFI taxa for coral species effects. Forest plot showing negative binomial GLM species-effect chi-square (Type II LR) for 52 testable OTUs (present in >=5 genotyped corals), controlling for log(volume) and site. OTUs colored by mediation classification: blue = genuine genotype effect (species effect persists after adding branch_width; n = 13); orange = architecture-mediated (species effect eliminated by branch_width; n = 15); gray = non-significant (n = 24). Dashed vertical line at FDR q = 0.10 significance threshold. Top genuine genotype effects include *Galeropsis monodonta* (corallivore, strongly elevated on *P. verrucosa*), *Harpiliopsis beaupresii* (shrimp, elevated on non-grandis hosts), and *Fennera* sp. (absent from *P. tuahiniensis*/*P. verrucosa*). Top architecture-mediated effects include *Paragobiodon modestus* (goby), *Alpheus lottini* (snapping shrimp), and *Trapezia tigrina* (guard crab).

![](../output/figures/supplement/figS13_genotype_screen.png)

---

**Figure S14.** Trapezia body size, architectural filtering, and niche partitioning. (A) Trapezia carapace width versus coral branch architecture (tight vs wide branching); *T. rufopunctata* (largest species, mean 17.0 mm) occurs almost exclusively in wide-branched corals (92%), while *T. tigrina* and *T. punctimanus* concentrate in tight-branched corals (chi-square p = 0.0002). (B) Character displacement in *T. serenei*: individuals are significantly smaller when co-occurring with other Trapezia species (6.2 vs 7.5 mm alone, p = 0.001), consistent with size-structured coexistence. (C) Body size ratios in two-species Trapezia assemblages (n = 33 corals); mean ratio = 1.66, consistent with size-ratio coexistence theory.

![](../output/figures/supplement/figS14_trapezia_body_size.png)

---

**Figure S15.** Community assembly: beta-dispersion convergence across coral size classes. Distance to group centroid (Bray–Curtis) decreases from small (0.576 ± 0.091) to medium (0.515 ± 0.073) to large (0.506 ± 0.063) corals (ANOVA F = 9.21, p < 0.001). This pattern indicates that CAFI communities are less variable on larger corals, though this trend does not survive rarefaction (Fig. S3B) and likely reflects statistical convergence of richer samples rather than deterministic filtering (see main text). n = 112 corals across three size terciles; points show individual coral distances overlaid on boxplots.

![](../output/figures/supplement/figS15_beta_dispersion.png)

---

**Figure S16.** Community assembly null-model analysis. (A) Histogram of pairwise Raup–Crick dissimilarity values (mean = 0.19, far below the null expectation of 0.5; dashed red line). Communities are significantly more similar than expected under stochastic assembly (t-test vs. 0.5: p < 2 × 10⁻¹⁶), indicating deterministic environmental filtering. (B) Raup–Crick dissimilarity by coral size class (violin plots). All three size classes show values well below 0.5 (Small: 0.23; Medium: 0.23; Large: 0.19), with no significant trend across size classes. (C) Beta-dispersion by size class (boxplots; same data as Fig. S15). (D) Variation partitioning: host architecture (branch width + haplotype; 5.6% unique), coral size (4.7%), and space (site; 2.9%) all significantly structure CAFI composition (all db-RDA p = 0.001), with 85% residual.

![](../output/figures/supplement/figS16_community_assembly.png)

---

**Figure S17.** Taxonomic community structure along the coral size gradient. Standardized effect size of mean pairwise taxonomic distance (SES.MPD, equivalent to Net Relatedness Index) plotted against log(colony volume). The global mean SES.MPD is significantly negative (−1.18, p < 10⁻¹³), indicating that co-occurring CAFI are more taxonomically diverse than expected under random assembly — consistent with limiting similarity or interspecific competition. The trend toward more negative NRI values in larger corals (β = −0.17, R² = 0.02, p = 0.11) is suggestive but not significant. Dashed lines at ±1.96 mark conventional significance thresholds. Points colored by size class (blue = Small, orange = Medium, green = Large). n = 112 corals.

![](../output/figures/supplement/figS17_taxonomic_structure.png)

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

**Table S13.** *Removed (morphotype covariate sensitivity). Field morphotype assignments showed limited concordance with genetic species (33--94%; Table S15A), so this analysis was not retained.*

**Table S14.** Missing data characterization: logistic regression of dropout predictors.

*Data: `output/tables/sensitivity_missing_data_predictors.csv`*

**Table S15.** Haplotype sensitivity analyses. Molecular species assignments (mtORF haplotyping per @JohnstonCunningBurgess2022) were available for 101 of 114 colonies (99 of the 112 with volume data) and used to test whether genetic species identity confounds the main results. (A) **Morphotype--haplotype concordance.** Cross-tabulation of field morphotype assignment vs. genetic species. Concordance varied among morphotypes: "eudoxi" mapped to *P. grandis* in 94% of cases (45/48), "meandrina" to *P. meandrina* in 66% (19/29), and "verucosa" to *P. verrucosa* in only 33% (7/21), confirming that field morphotype is an unreliable proxy for genetic species identity. (B) **BEF with haplotype covariate.** Adding genetic species identity as a covariate to the richness → condition model reduced the richness coefficient by only 8% (β = 0.081, p = 0.008 without haplotype; β = 0.074, p = 0.010 with haplotype; n = 74 with both physiology and haplotype data). The richness signal is robust to genetic species identity. (C) **Scaling by genetic species.** The volume × species interaction was not significant (p = 0.39), indicating that scaling exponents do not differ among species. Species-specific estimates were all sublinear: *P. grandis* β = 0.46 (n = 49), *P. meandrina* β = 0.62 (n = 32 of 34 with volume data), *P. verrucosa* β = 0.73 (n = 10, not significant, p = 0.17). (D) **Composition by genetic species.** Marginal PERMANOVA (community ~ log(volume) + genetic species + site) confirmed that genetic species explains 9.8% of composition variation (F = 2.95, p = 0.001), comparable to volume (8.1%, p = 0.001) and exceeding site (6.0%, p = 0.001). (E) **Phylogenetic distance and symbiont identity.** Mantel test detected a significant correlation between phylogenetic distance and CAFI community dissimilarity (r = 0.12, p = 0.013), which strengthened after controlling for volume (partial Mantel r = 0.13, p = 0.002). Putative symbiont identity (*Cladocopium latusorum* vs. *C. pacificum*, inferred from published host species--symbiont associations in @Burgess2021, not directly measured) explained 3.8% of composition (PERMANOVA R² = 0.038, p = 0.001) and predicted coral condition (β = 1.28, p = 0.008; n = 74). Because symbiont identity was imputed from host species, these results are not independent of host species effects and should be interpreted as reframing the host species signal rather than as evidence for symbiont-specific mechanisms.

*Data: `output/tables/sensitivity_morphotype_haplotype_concordance.csv`, `sensitivity_haplotype_bef.csv`, `sensitivity_scaling_by_species.csv`, `sensitivity_composition_by_species.csv`, `sensitivity_phylogenetic_symbiont.csv`*

**Table S16.** Host species identity effects on CAFI composition. (A) Marginal PERMANOVA (Bray-Curtis, 999 permutations; n = 98 genotyped corals): coral species R² = 0.083, volume R² = 0.082, site R² = 0.053, all p = 0.001. (B) Pairwise PERMANOVA between all six coral species pairs (FDR-corrected; all p < 0.008). (C) Variance partitioning (adjusted R²): coral species alone 5.4%, volume alone 4.9%, site alone 4.1%, all shared fractions <0.2%. (D) Species vs architecture: species conditioned on branch_width = 2.3% (p = 0.001); branch_width conditioned on species = -0.1% (p = 0.65). Source script: `16_genotype_cafi_analysis.R`.

*Data: `output/tables/genotype_permanova_marginal.csv`, `genotype_pairwise_permanova.csv`, `genotype_variance_partitioning.csv`*

**Table S17.** Systematic screen of 52 CAFI OTUs for coral species effects. For each OTU: Model A (species + volume + site) and Model B (+ branch_width) LR chi-square and p-values for the species term, FDR-corrected q-value, mediation classification (genuine genotype / architecture-mediated / non-significant), and direction of effect for each coral species relative to *P. grandis*. Source: `genotype_species_screen.csv`.

*Data: `output/tables/genotype_species_screen.csv`*

**Table S18.** Indicator species and host specificity. (A) Indicator species (multipatt, IndVal.g): 19 significant indicators at p <= 0.05 -- 8 for *P. verrucosa*, 5 for *P. tuahiniensis*, 2 for *P. grandis*, 4 multi-host. (B) Host specificity index (SIS): 9 specialists (SIS < -1.96, all concentrated on *P. grandis*), 1 generalist (*Breviturma pica*, SIS = +2.24), 54 neutral. (C) Co-occurrence network structure by host species. Source: `genotype_indicator_species.csv`, `genotype_host_specificity.csv`, `genotype_cooccurrence_networks.csv`.

*Data: `output/tables/genotype_indicator_species.csv`, `genotype_host_specificity.csv`, `genotype_cooccurrence_networks.csv`*

**Table S19.** Trapezia body size, architectural filtering, and pair composition. (A) Body size by Trapezia species x coral species (n = 659 individuals with carapace width). (B) Species x architecture chi-square contingency. (C) Pair composition: per-coral Trapezia assemblage characterization (n = 114 corals). (D) Character displacement: body size of *T. serenei* alone vs co-occurring. Source: `trapezia_body_size_genotype.csv`, `trapezia_pair_composition.csv`.

*Data: `output/tables/trapezia_body_size_genotype.csv`, `trapezia_pair_composition.csv`*

**Table S20.** Community assembly results summary. Test statistics for nine community assembly tests (including separate NRI and NTI entries) across 112 *Pocillopora* colonies. Results are from `scripts/15_community_assembly.R`; full details in `output/tables/assembly_*.csv`.

| Hypothesis | Test | Statistic | p-value |
|---|---|---|---|
| Deterministic vs. stochastic | Raup-Crick vs. 0.5 | mean RC = 0.19 | < 2 × 10⁻¹⁶ |
| Dispersal limitation | Mantel (Bray ~ Geo) | r = 0.10 | 0.031 |
| Dispersal (controlling for size) | Partial Mantel | r = 0.11 | 0.017 |
| Size effect (variation partitioning) | db-RDA: log(volume) | F = 9.19 | 0.001 |
| Architecture effect | db-RDA: branch + haplotype | F = 2.09 | 0.001 |
| Space effect | db-RDA: site | F = 2.87 | 0.001 |
| Community convergence | Beta-dispersion ANOVA | F = 9.21 | < 0.001 |
| Taxonomic overdispersion (NRI) | SES.MPD vs. 0 | mean = −1.18 | < 10⁻¹³ |
| Taxonomic overdispersion (NTI) | SES.MNTD vs. 0 | mean = −1.53 | < 2 × 10⁻¹⁶ |

---

## Supplementary Methods

### Multiple-testing correction and model diagnostics

We applied three levels of multiple-testing correction: (1) Hochberg step-up for the two pre-specified BEF predictors (richness and total abundance; k = 2); (2) Benjamini--Hochberg FDR for four exploratory BEF predictors (Shannon diversity, *Trapezia* abundance, resident fish abundance, and *Galeropsis monodonta* abundance; k = 4); and (3) FDR across all species in species-level analyses. General model diagnostics included simulated residuals (`DHARMa`), Cook's distance (4/n), and variance inflation factors. For OLS BEF models, Breusch--Pagan tests assessed homoscedasticity (all BP p > 0.75, confirming homoscedastic residuals). HC3-robust standard errors were computed as a sensitivity check: richness significance changed from p = 0.018 (OLS) to p = 0.14 (HC3), and abundance from p = 0.048 to p = 0.25. HC3 is known to be overly conservative at n < 100 (Long & Ervin 2000), and the Breusch--Pagan diagnostics do not indicate heteroscedasticity; we therefore retain OLS as the primary inference but report HC3 results for transparency. The PERMANOVA volume effect was also tested with site-stratified permutations (permuting within site blocks), yielding R² = 0.078, p = 0.005 — still significant but more conservative than the unstratified test (p < 0.001).

### Species- and group-level scaling analyses

The scaling analysis was repeated for 21 prevalent species (≥30 individuals and ≥15% prevalence) and six taxonomic groups (*Trapezia* crabs, shrimps, other crabs, gastropods, fish, and other invertebrates). We computed an inverse-variance-weighted mean β as a fixed-effect summary (weights = 1/SE²); because species co-occur on the same corals, this summary may underestimate uncertainty. We applied FDR correction within species-level and group-level test families and characterized size-dependent occurrence using logistic GLMs (presence ~ log(volume) + site) for 24 species with ≥15% prevalence. For community-level abundance, β = 1 was evaluated with both Wald and site-stratified bootstrap tests (1,000 iterations); bootstrap BCa confidence intervals were used when estimable, with percentile intervals as fallback.

### Additional BEF sensitivity analyses

Count-based predictors (total abundance, *Trapezia*, fish, and *Galeropsis* counts) were square-root transformed; species richness was untransformed. Beyond the main forward models, we fit joint models containing richness and √abundance, hierarchical R² variance partitioning, and a piecewise structural equation model with volume → richness → condition and volume → abundance → condition pathways (z-scored predictors; Fisher's C fit). We also tested mediation with bootstrap mediation analysis (`mediation`; 1,000 iterations; treatment = richness, mediator = √abundance, outcome = condition). As sensitivity checks, we added coral morphotype and genetic species identity (mtORF haplotype) as covariates in separate models. Species-specific contributions were tested for 19 prevalent species (≥5 corals; trait ~ √(species abundance) + log(volume) + site) across five condition measures, with FDR correction across all 95 tests. Reverse models (condition → CAFI) tested seven response variables with FDR correction, and cross-study concordance with the companion experiment was evaluated with a binomial test.

### Additional neighborhood analyses

For the 61 neighborhood-surveyed colonies, the full model included log(volume), neighbor count within 5 m, total neighbor volume, mean inter-colony distance, and site. AIC-based backward elimination identified the supported predictor subset for abundance, richness, and Shannon diversity models. We also tested size × neighborhood interactions and functional-group responses with FDR correction and assessed compositional variability across neighbor density with PERMDISP.

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

Design-sensitivity calculations were conducted using the `pwr` package (`pwr.f2.test`). For the CAFI--condition feedbacks (Q3; n = 84), power was computed to detect the incremental variance explained by a single CAFI predictor beyond the covariates (log(volume) + site), with 1 numerator degree of freedom and 79 denominator degrees of freedom (84 − 5 model parameters). Power to detect medium (ΔR² = 0.10) effects was 84%, and small (ΔR² = 0.05) effects 53%, at α = 0.05. The minimum detectable ΔR² at 80% power was 0.09. To assess cross-study comparability, we computed power to detect the companion experiment's observed effect size (ΔR² ≈ 0.12): the survey had 91% power. For the rarefied richness test (n = 47; 1 numerator df, 42 denominator df), power was substantially lower: 58% for medium effects (ΔR² = 0.10), 32% for small effects (ΔR² = 0.05), and only 26% for effects the size of that observed in the full sample (ΔR² = 0.04). The minimum detectable ΔR² at 80% power was 0.16, indicating that the rarefied subset had limited sensitivity to detect the modest effects found in the main analysis. For the neighborhood analysis (n = 61; 1 numerator df, 55 denominator df), power was 70% for medium effects (ΔR² = 0.10) and 40% for small effects (ΔR² = 0.05). Full values are exported in `power_analysis.csv`.

### Taxonomy sensitivity analysis

To assess whether results depended on taxonomic resolution, we repeated seven key analyses under five resolution scenarios: (1) baseline (243 OTUs at the finest available resolution), (2) species-only (154 species-level OTUs; higher-level identifications excluded), (3) merge-up (specimens identified to genus or family merged with congeneric or confamilial species-level OTUs), (4) lump-down (all OTUs lumped to genus level, reducing resolution), and (5) rare-excluded (OTUs with fewer than 3 total individuals removed). For each scenario, we recomputed: abundance scaling exponent (β), richness scaling exponent (z), PERMANOVA R² for site and volume, Shannon diversity scaling slope, rarefied richness--condition p-value, and network modularity (Q). All scenarios produced qualitatively identical results (Fig. S6; Table S8).

### Sensitivity analyses

Four supplementary sensitivity analyses addressed methodological gaps identified during pre-submission audit (Tables S11--S14; script `14_supplementary_sensitivity.R`). A fifth set of haplotype sensitivity analyses (Table S15) is described separately below. (1) **Beta diversity partitioning** (betapart; Baselga 2010) decomposed total Sørensen dissimilarity into turnover (Simpson) and nestedness-resultant components along the coral size gradient (Table S11). (2) **Mediation analysis** (mediation package; 1,000 bootstrap iterations) tested whether the richness → condition effect operates through the abundance pathway (ACME; Table S12). (3) **Morphotype covariate** sensitivity was planned but not retained because field morphotype assignments showed limited concordance with genetic species (33--94%; Table S15A). (4) **Missing data characterization** used logistic regression to test whether the 114 → 84 dropout for condition analyses was predicted by volume, site, richness, or abundance (Table S14). Additional sensitivity analyses (envfit species vectors, indicator species by size class, iNEXT coverage-based rarefaction, community-weighted mean obligate-mutualist fraction, nonlinear BEF, C-score community co-occurrence, AIC model averaging) were computed and are archived with the analysis code in the data repository.

### Haplotype sensitivity analyses

Molecular species assignments were obtained via mtORF haplotyping following @JohnstonCunningBurgess2022. PCR amplification and sequencing of the mitochondrial open reading frame successfully genotyped 101 of 114 colonies (11 PCR failures, 2 without tissue samples); 99 of the 112 colonies with volume data had valid haplotypes. Genotyped colonies resolved to four genetic species: *Pocillopora grandis* (n = 49), *P. meandrina* (n = 32 of 34 with volume data), *P. tuahiniensis* (n = 7), and *P. verrucosa* (n = 10), with 1 colony unresolved. Phylogenetic distances among host species were derived from the *Pocillopora* mtORF gene tree. Symbiont identity (*Cladocopium latusorum* vs. *C. pacificum*) was assigned based on the host species--symbiont associations documented by @Burgess2021, which showed strong species-specificity in Mo'orean *Pocillopora*. Five analyses tested whether genetic species identity confounds the main results (Table S15): (A) morphotype--haplotype concordance via chi-square test with simulated p-value (confirming that field morphotype is unreliable); (B) BEF model with genetic species as a covariate; (C) volume × species interaction test for scaling differences among genetic species; (D) marginal PERMANOVA with genetic species as a predictor of community composition; and (E) Mantel and partial Mantel tests for phylogenetic distance--community dissimilarity correlations, plus PERMANOVA and linear model tests for symbiont effects on composition and condition.

### Host species identity and CAFI associations

To test whether coral genetic species identity predicts CAFI community structure beyond colony volume and site, we conducted four complementary analyses on the 98 genotyped corals retained after excluding the single ambiguous *P. meandrina/grandis* haplotype and 2 colonies without volume. First, we ran marginal PERMANOVA (Bray-Curtis, 999 permutations) with coral species, log(volume), and site as predictors, followed by pairwise species comparisons (FDR-corrected). Three-way variance partitioning (vegan::varpart, Hellinger-transformed) decomposed community variation among coral species, volume, and site. We also partitioned species vs branch_width to test architecture mediation. Second, we screened all 52 CAFI OTUs present in >=5 genotyped corals using negative binomial GLMs (Poisson fallback) with Type II likelihood ratio tests for the coral species term (BH-FDR across OTUs). Each significant OTU was classified based on a covariate-sensitivity test: "species effect robust to architecture" (species term persists after adding branch_width) or "species effect absorbed by architecture" (species term eliminated by branch_width). Because species identity is upstream of branch_width in *Pocillopora* (species determines branching morphology), adding branch_width may over-control for genuinely genetic effects; these classifications should be interpreted as sensitivity tests, not causal mediation. Third, we identified indicator species using indicspecies::multipatt (func = "IndVal.g") and computed a standardized host specificity index comparing observed Simpson's diversity of host usage against a null model (999 iterations shuffling CAFI across corals while preserving coral species frequencies). Fourth, Trapezia body size (carapace width) was analyzed using LMMs (crab species + branch_width + coral species + volume + site as fixed effects, coral_id as random intercept). Trapezia species x architecture filtering was tested via chi-square with simulated p-value, and size-based habitat partitioning was assessed by comparing *T. serenei* body size when alone vs co-occurring with congeners (Wilcoxon rank-sum test).

### Co-occurrence patterns

We tested pairwise co-occurrence using a volume-weighted Bernoulli null model (10,000 iterations; restricted to species present on ≥10 coral colonies, yielding 33 species and 528 pairwise combinations; Stier et al. 2012) with FDR correction. Intraspecific density patterns and size-dependent co-occurrence were assessed with multinomial and tercile-stratified null models (see Volume-weighted co-occurrence null model and Intraspecific density null model sections above).

### Community assembly

We assessed the balance between deterministic and stochastic assembly using the Raup--Crick metric (vegan::raupcrick; 999 iterations; Chase et al. 2011), where values near 0 indicate stochastic assembly and values departing from 0 indicate deterministic processes. We tested taxonomic structure using standardized effect sizes of mean pairwise distance (SES.MPD / Net Relatedness Index; picante::ses.mpd; 999 randomizations). Taxonomic distances were computed from the Linnaean hierarchy (species, genus, family, order, class, phylum) using UPGMA clustering, providing a proxy for phylogenetic distances in the absence of a resolved molecular phylogeny for the full CAFI assemblage. Distance-decay of community similarity was tested with Mantel tests (Spearman correlation between Bray--Curtis dissimilarity and Haversine geographic distance; 9,999 permutations) on 57 colonies with GPS coordinates, with a partial Mantel controlling for coral size. We note that Mantel tests have known limitations, including inflated Type I error rates under certain spatial autocorrelation structures (Guillot & Rousset 2013) and low power for detecting nonlinear distance-decay relationships; results should be interpreted as complementary evidence alongside the variation partitioning and null-model analyses.

### Sensitivity analyses (main text)

We tested all results for robustness across five taxonomic resolution scenarios (baseline 243 taxa, species-only, merge-up, lump-down, rare-excluded; Fig. S6; Table S8) and assessed spatial autocorrelation using Moran's I (Table S9).

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
| **Table S13** (Removed -- morphotype unreliable) | — | — |
| **Table S14** (Missing data predictors) | `sensitivity_missing_data_predictors.csv` | `14_supplementary_sensitivity.R` |
| **Table S15** (Haplotype sensitivity) | `sensitivity_morphotype_haplotype_concordance.csv`, `sensitivity_haplotype_bef.csv`, `sensitivity_scaling_by_species.csv`, `sensitivity_composition_by_species.csv`, `sensitivity_phylogenetic_symbiont.csv` | `14_supplementary_sensitivity.R` |
| **Table S16** (Host species identity effects) | `genotype_permanova_marginal.csv`, `genotype_pairwise_permanova.csv`, `genotype_variance_partitioning.csv` | `16_genotype_cafi_analysis.R` |
| **Table S17** (CAFI OTU species screen) | `genotype_species_screen.csv` | `16_genotype_cafi_analysis.R` |
| **Table S18** (Indicator species and host specificity) | `genotype_indicator_species.csv`, `genotype_host_specificity.csv`, `genotype_cooccurrence_networks.csv` | `16_genotype_cafi_analysis.R` |
| **Table S19** (Trapezia body size and architecture) | `trapezia_body_size_genotype.csv`, `trapezia_pair_composition.csv` | `16_genotype_cafi_analysis.R` |
| **Table S20** (Community assembly) | `assembly_*.csv` | `15_community_assembly.R` |

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
