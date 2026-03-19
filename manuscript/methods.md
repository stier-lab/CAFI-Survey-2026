# Methods

## Study system and survey design

We surveyed 114 *Pocillopora* spp. coral colonies (morphotype identifications; not genetically confirmed) across three reef habitats on Mo'orea, French Polynesia (17°30'S, 149°50'W) during June--August 2019 (Fig. 1A). Sites were selected to span the range of reef environments: Hauru (HAU; fringing reef, north shore; n = 38), Maatea (MAT; lagoon/back reef, east shore; n = 39), and Maharepa (MRB; barrier reef, north shore; n = 37). At each site, colonies were selected non-randomly to span the full available size range, with divers targeting small, medium, and large colonies to ensure adequate coverage of the three-order-of-magnitude size gradient. Two MRB colonies lacked volume measurements and were excluded from size-dependent analyses (final n = 112; MRB n = 35).

## Coral measurements

Colony size was estimated as ellipsoidal volume: V = (4/3) × π × (L/2) × (W/2) × (H/2), where L, W, and H are the maximum length, width, and height (cm). While 3D photogrammetry provides more precise volume estimates (Curtis et al. 2023), the ellipsoidal approximation is standard for *Pocillopora* colonies, which have roughly symmetrical branching morphology. This approximation overestimates true colony volume because it includes interstitial space between branches, but the bias is approximately proportional across colony sizes and does not affect log-linear scaling estimates. Colony volumes ranged from 21 to 42,333 cm³ (Fig. 1B). <!-- CRITIQUE: METH-06 fixed: named the size-only subset and clarified exclusion -->
For a subset of 61 colonies ("neighborhood corals"), all *Pocillopora* neighbors within a 5-m radius were counted and their distances measured. The remaining 53 colonies were surveyed for CAFI and colony size only and are excluded from neighborhood analyses. Neighborhood density ranged from 1 to 76 neighbors (median = 17; Fig. 1C).

## CAFI collection and identification

All coral-associated fishes and invertebrates (CAFI) were extracted by enclosing each colony in a fine-mesh bag, carefully detaching the coral from the substrate, and preserving the entire colony and bag contents in 95% ethanol. No chemical anesthetics (e.g., clove oil) were used. Boring and deeply cryptic fauna embedded within the coral skeleton were not included. Specimens were identified to the lowest practical taxonomic level: 154 to species and 89 to genus, family, or higher, yielding 243 taxa from 3,989 total individuals. Species richness and Shannon diversity (H') were computed per colony using the full taxon set. Rarefied richness (expected species at n = 20 individuals) was computed using vegan::rarefy (Oksanen et al. 2022) to control for the well-documented correlation between abundance and raw richness.

## Coral physiological condition

Tissue samples were collected from the upper branch of each colony for physiological assays of total protein (mg cm⁻²), carbohydrate (mg cm⁻²), zooxanthellae density (cells cm⁻²), and ash-free dry weight (mg cm⁻²). Position correction removed a tissue-gradient sampling artifact: sampling nubbins from smaller colonies captured a larger fraction of the branch axis, including more of the lower-density tip tissue, creating a systematic size bias in raw trait values. We regressed each trait on stump length and nubbin length and used the residuals. Coral condition was then quantified as the first principal component (PC1) of a principal component analysis (PCA) on the four position-corrected traits (61.6% of variance; all four traits loaded positively). Physiology data were available for 108 colonies; after merging with complete CAFI and volume data, 84 colonies were retained for BEF analyses.

## Statistical analyses

<!-- CRITIQUE: METH-01 fixed: replaced placeholder DOIs with acceptance-contingent language -->
All analyses were conducted in R v4.4 (session information including package versions archived at [DOI to be assigned upon acceptance]). We used natural logarithms throughout; site was included as a fixed effect (three sites provide too few levels for reliable random intercept estimation; Bolker et al. 2009). <!-- CRITIQUE: METH-04 fixed: expanded multiple-testing correction with explicit predictor names -->
We applied three levels of multiple-testing correction: (1) Holm correction for the two pre-specified BEF predictors (richness and total abundance; k = 2); (2) Benjamini--Hochberg FDR for four exploratory CAFI predictors (Shannon diversity, *Trapezia* abundance, fish abundance, and *Galeropsis* abundance; k = 4); and (3) Benjamini--Hochberg FDR across all species tested in species-level analyses. Model fit for GLMs was assessed using simulated residuals (DHARMa package), Cook's distance (threshold: 4/n), and variance inflation factors. All tests used α = 0.05. Spatial independence was confirmed using Moran's I on inverse-distance weights (Table S10). All code is available at [URL to be provided upon acceptance].

### Q1: Scaling of CAFI with coral size

We tested how CAFI abundance scales with coral volume by fitting the power-law relationship N = a × V^β^ using a negative binomial generalized linear model (GLM; log link):

    total_cafi ~ log(volume) + site

where the coefficient on log(volume) directly estimates the scaling exponent β. We tested three hypotheses: Field of Dreams (β = 1; abundance scales proportionally with volume), Propagule Redirection (β < 1; sublinear scaling indicates density dilution; Stier & Osenberg 2010), and super-linear scaling (β > 1; larger corals disproportionately attractive). The null hypothesis (β = 1) was tested using both a Wald z-test (z = [β − 1]/SE) and a two-sided bootstrap proportion test (proportion of 1,000 site-stratified bootstrap replicates with β ≥ 1 or β ≤ 1, doubled). Bootstrap CIs are the primary inference; Wald z-tests are reported for comparability with prior studies. Confidence intervals were computed using bias-corrected and accelerated (BCa) bootstrap with site-stratified resampling (1,000 iterations; 1,000 replicates are adequate for percentile CIs, though conservative for BCa, and results were stable across replicate runs). Where BCa acceleration constants could not be estimated, percentile CIs were used as a fallback. Species richness was modeled analogously using a Poisson GLM. Rarefied richness was regressed on log(volume) using ordinary least squares. <!-- CRITIQUE: METH-03 fixed: clarified per-capita density computation and log-log slope relationship -->
Per-capita CAFI density (individuals cm⁻³; total CAFI / colony volume) was computed per colony; the log--log slope of density on volume equals β − 1 by definition and was not independently modeled.

We repeated the scaling analysis for each of the 21 most prevalent individual species (≥30 total individuals and ≥15% prevalence) and for six taxonomic groups (*Trapezia* crabs, shrimps, other crabs, gastropods, fish, other invertebrates). An inverse-variance-weighted mean β across all species was computed as a fixed-effect summary (weights = 1/SE²); because species co-occur on the same corals, this assumes independence among species-level estimates and the SE of the weighted mean may be underestimated. FDR correction was applied within species-level and group-level test families. <!-- CRITIQUE: METH-07 fixed: explained 21 vs 24 species threshold difference -->
To characterize size-dependent occurrence, we fit logistic GLMs (presence ~ log(volume) + site) for 24 species with ≥15% prevalence and applied FDR correction across all tests (Fig. S12). (These 24 species met the ≥15% prevalence threshold; the 21-species scaling set additionally required ≥30 total individuals.)

### Q2: Community composition

Community composition was analyzed on all 112 colonies with volume data. Species abundances were Hellinger-transformed (square root of relative abundances) to down-weight the influence of the most abundant species on dissimilarity calculations (Legendre & Gallagher 2001). We tested for compositional differences among sites using marginal (Type III) permutational multivariate analysis of variance (PERMANOVA; vegan::adonis2) on Bray--Curtis dissimilarities with 999 permutations, including log(volume) and site as predictors. Multivariate dispersion homogeneity was assessed using permutational analysis of multivariate dispersions (PERMDISP; vegan::betadisper). Non-metric multidimensional scaling (NMDS) ordination (2 dimensions) was computed on all 112 colonies with volume data. All statistical tests (PERMANOVA, db-RDA) used the full dataset (n = 112). Species associated with site differences were identified by fitting species vectors to the NMDS ordination (envfit; 999 permutations; R² > 0.10, p < 0.01). Pairwise site comparisons used separate PERMANOVA tests with Bonferroni correction (k = 3). PERMANOVA results were validated across five distance metrics (Bray--Curtis on Hellinger and Wisconsin transformations, Jaccard on presence--absence, Gower, and Raup--Crick) and via 500 balanced site-subsampling iterations (Fig. S2).

To test whether composition shifts continuously along the coral size gradient, we used distance-based redundancy analysis (db-RDA; vegan::dbrda) with log(volume) as the constrained variable and site partialed out (Condition). The db-RDA F-statistic is algebraically identical to the marginal PERMANOVA F when both condition on the same covariates. Variance partitioning (vegan::varpart) decomposed community variation into volume-only, site-only, shared, and residual fractions. Species scores on the constrained axis were computed as weighted averages of site scores. Robustness was verified by repeating the db-RDA on iterated-rarefied distances: for each of 100 iterations, communities were rarefied to the minimum observed abundance among colonies with ≥5 individuals, a Bray--Curtis distance matrix was computed, and the resulting matrices were averaged element-wise before ordination and hypothesis testing. Nestedness along the size gradient was tested using NODF (nestedness metric based on overlap and decreasing fill; Almeida-Neto et al. 2008) against 999 quasiswap null matrices (an algorithm preserving row and column totals; Miklós & Podani 2004).

### Co-occurrence patterns (Supplement)

Pairwise co-occurrence was tested using a volume-weighted Bernoulli null model (Stier et al. 2012). For each species, a logistic GLM predicted per-coral occurrence probability from log(volume) + site. For each of 10,000 null iterations, species presence was drawn independently from Bernoulli distributions with predicted probabilities. The standardized effect size (SES) of observed co-occurrence was computed as (observed − null mean) / null SD, with false discovery rate (FDR) correction across all pairs. Species with ≥10 occurrences were included. Intraspecific density patterns were tested with a multinomial allocation null model to assess whether conspecifics aggregate on fewer corals than expected (consistent with mating-pair formation; Stier et al. 2012). To assess whether co-occurrence patterns vary with coral size, we repeated the volume-weighted null model separately for three coral size classes (volume terciles), fitting logistic regressions within each class to obtain class-specific occurrence probabilities (Fig. S9C).

### Q3: CAFI diversity and coral condition (BEF framework)

We tested whether CAFI diversity predicts coral condition using a biodiversity--ecosystem function (BEF) framework (Tilman et al. 2014). The complementarity hypothesis predicts that more CAFI species provide non-redundant functional benefits (e.g., defense, cleaning, nutrient cycling), generating a positive diversity--condition relationship beyond what abundance alone predicts.

Forward models (CAFI → condition) took the form:

    condition_PC1 ~ CAFI_predictor + log(volume) + site

with ordinary least squares (OLS) standard errors as primary inference. Breusch--Pagan tests confirmed homoscedasticity for all models (BP p > 0.5); HC3 heteroscedasticity-consistent standard errors (sandwich package) are reported as a supplement sensitivity check (conservative at n < 100; Long & Ervin 2000). Count-based CAFI predictors (total abundance, *Trapezia* count, fish count, *Galeropsis* count) were square-root transformed prior to fitting to reduce right skew. Species richness was not transformed because it has a more linear relationship with condition and lower skewness than raw abundance counts.

Community composition (PC1_CAFI) was tested separately in the supplement.

To disentangle diversity from abundance effects, we conducted three complementary analyses:
1. **Partial regression**: condition ~ richness + √(abundance) + log(volume) + site, with variance inflation factor (VIF) diagnostics to assess collinearity.
2. **Variance partitioning**: Hierarchical R² decomposition quantifying variance uniquely attributable to richness, uniquely to abundance, and shared (confounded) between them.
3. **Path model**: A piecewise structural equation model (piecewiseSEM; Lefcheck 2016) with three component OLS models: abundance ~ log(volume) + site, richness ~ log(volume) + site, and condition ~ richness + abundance + log(volume) + site. Predictors were z-scored for standardized path coefficients. Model fit was evaluated using Fisher's C statistic.

<!-- CRITIQUE: METH-02 fixed: explained n=47 drop for rarefied richness comparison -->
To test whether the richness → condition relationship was an abundance artifact, we compared models using raw versus rarefied richness (expected species at n = 20; vegan::rarefy). Rarefied richness was undefined for colonies with fewer than 20 CAFI individuals, reducing the available sample from 84 to 47 for this comparison. We address the interpretive ambiguity of this test in the Discussion.

To test species-specific contributions to individual condition traits, we fitted linear models for each of the 19 most abundant species (present on ≥5 corals in the physiological dataset, n = 84): trait ~ √(species abundance) + log(volume) + site, where trait was one of five measures (protein, carbohydrate, zooxanthellae density, AFDW, or condition PC1). Standardized β coefficients were compiled in a heatmap (Fig. S13), and the nine strongest associations were visualized as scatter plots (Fig. S14). FDR correction was applied across all 95 species × trait tests.

Reverse models (condition → CAFI) tested seven response variables: total CAFI abundance, species richness, Shannon diversity, *Trapezia* abundance, fish abundance, *Galeropsis monodonta* abundance, and community composition (PC1_CAFI). Count responses used negative binomial or Poisson GLMs; continuous responses used OLS. BH-FDR correction was applied across all seven tests. We tested the 8 key species from the companion experiment that had ≥5 occurrences in our physiological dataset (of 10 originally defined) individually for condition effects, with FDR correction. Cross-study concordance of effect directions was evaluated with a one-sided binomial test (H₀: 50% concordance by chance). A sensitivity power analysis (minimum detectable R² increment at 80% power, α = 0.05, n = 84) assessed whether the survey sample was sufficient to detect the companion experiment's observed effect sizes.

### Neighborhood effects (Supplement)

We tested whether neighborhood density affects CAFI abundance and coral condition using the 61 neighborhood-surveyed corals. Full models included log(volume), n_neighbors, total_neighbor_volume, mean_neighbor_dist, and site as predictors. Power analysis indicated 65% power for medium effects (R² = 0.10) at α = 0.05. Compositional variability across neighborhood density categories was assessed using PERMDISP (betadisper) on Bray--Curtis dissimilarities.

### Sensitivity analyses

Results were tested for robustness across five taxonomic resolution scenarios: baseline (243 taxa), species-only (154 taxa), merge-up (higher-level IDs merged to nearest named taxon), lump-down (all IDs resolved to genus), and rare-excluded (singletons removed). Six metrics were evaluated per scenario: abundance scaling exponent (β), richness scaling exponent (z), Shannon scaling slope, PERMANOVA R² (site and volume, counted separately), and rarefied richness--condition slope (Fig. S6; Table S9). Spatial autocorrelation was assessed using Moran's I on CAFI abundance, richness, and Shannon diversity (Table S10).

### Ethical statement

Research was conducted under permits issued by the Haut-commissariat de la République en Polynésie française and in collaboration with the Université de la Polynésie française (CRIOBE). All protocols followed institutional and local regulations for scientific collection on the reefs of Mo'orea. Coral collections were conducted as part of a broader research programme at the UC Berkeley Richard B. Gump South Pacific Research Station, <!-- CRITIQUE: METH-09 fixed: clarified IACUC applicability for invertebrate vs coral protocols -->
Invertebrate collections did not require IACUC review; coral tissue sampling protocols were approved under broader programme oversight at UCSB.

### Data accessibility

Data and analysis code are available at [DOI to be assigned upon acceptance]. Data and code will be archived in the Dryad Digital Repository upon acceptance.
