# Methods

## Study system and survey design

We surveyed 114 *Pocillopora* coral colonies across three reef habitats on Mo'orea, French Polynesia (17°30'S, 149°50'W) during June–August 2019 (Fig. 1A). Sites were selected to span the range of reef environments: Hauru (HAU; fringing reef, north shore; n = 38), Maatea (MAT; lagoon/back reef, east shore; n = 39), and Maharepa (MRB; barrier reef, north shore; n = 35). At each site, colonies were haphazardly selected across a range of sizes. Two colonies lacked volume measurements and were excluded from size-dependent analyses (final n = 112). This survey complements a companion experimental study (Stier et al. in review) that manipulated coral density (1, 3, or 6 *Pocillopora* colonies per 25 m² reef) over 17 months; the survey tests whether scaling and feedback patterns observed under experimental colonization conditions hold in established, naturally assembled communities.

## Coral measurements

Colony size was estimated as ellipsoidal volume: V = (4/3) × π × (L/2) × (W/2) × (H/2), where L, W, and H are the maximum length, width, and height (cm). Colony volumes spanned more than three orders of magnitude (21–42,333 cm³; CV = 119%; Fig. 1B). For a subset of 61 colonies ("neighborhood corals"), all *Pocillopora* neighbors within a 5-m radius were counted and their distances measured. Neighborhood density ranged from 1 to 76 neighbors (median = 17; Fig. 1C).

## CAFI collection and identification

All coral-associated fauna and invertebrates (CAFI) were extracted by enclosing each colony in a fine-mesh bag, detaching the coral from the substrate, and preserving the contents in 95% ethanol. Specimens were identified to the lowest practical taxonomic level: 154 to species and 89 to genus, family, or higher, yielding 243 operational taxonomic units (OTUs) from 3,989 total individuals. Species richness and Shannon diversity (H') were computed per colony using the full OTU set. Rarefied richness (expected species at n = 20 individuals) was computed using vegan::rarefy (Oksanen et al. 2022) to control for the well-documented correlation between abundance and raw richness.

## Coral physiological condition

Tissue samples were collected from the upper branch of each colony for physiological assays of total protein (mg cm⁻²), carbohydrate (mg cm⁻²), zooxanthellae density (cells cm⁻²), and ash-free dry weight (mg cm⁻²). Coral condition was quantified as the first principal component (PC1) of the four position-corrected physiological traits, where position correction involved regressing each trait on stump length and nubbin length to remove a tissue-gradient sampling artifact (smaller corals required sampling more of the branch, including lower-density tips). Physiology data were available for 84 colonies.

## Statistical analyses

All analyses were conducted in R v4.4.x. We used natural logarithms throughout; site was included as a fixed effect (three sites are insufficient for random intercepts; Bolker et al. 2009). Multiple testing was corrected using a tiered approach: a priori BEF predictors (richness and abundance) were corrected with Hochberg's step-up procedure (FWER, k = 2); exploratory CAFI predictors were corrected with Benjamini–Hochberg FDR (k = 4); key species tests used Hochberg FWER to match the companion experiment. A global FDR across all test families is reported as a sensitivity check. All code is available at [repository URL].

### Q1: Scaling of CAFI with coral size

We tested how CAFI abundance scales with coral volume by fitting the power-law relationship N = a × V^β using a negative binomial generalized linear model (GLM):

    total_cafi ~ log(volume) + site

where the coefficient on log(volume) directly estimates the scaling exponent β. We tested three hypotheses: Field of Dreams (β = 1; abundance scales proportionally with volume), Propagule Redirection (β < 1; sublinear scaling indicates density dilution; Stier & Osenberg 2010), and super-linear scaling (β > 1; larger corals disproportionately attractive). The null hypothesis (β = 1) was tested using both a Wald z-test (z = [β − 1]/SE) and a bootstrap proportion test (two-sided p from 1,000 site-stratified bootstrap replicates). Confidence intervals were computed using bias-corrected and accelerated (BCa) bootstrap with site-stratified resampling (1,000 iterations). Species richness was modeled analogously using a Poisson GLM. Rarefied richness was regressed on log(volume) using ordinary least squares. Per-capita CAFI density (individuals cm⁻³) was computed to visualize the density dilution prediction. We repeated the scaling analysis for each of the 21 most prevalent individual species (≥30 total individuals and ≥15% prevalence) and for six taxonomic groups (Trapezia crabs, shrimps, other crabs, gastropods, fish, other invertebrates). FDR correction was applied within species-level and group-level test families. To characterize size-dependent occurrence, we fit logistic GLMs (presence ~ log(volume) + site) for 24 species with ≥15% prevalence and applied FDR correction across all tests (Fig. S15).

### Q2: Community composition

Community composition was analyzed on all 112 colonies with volume data. Species abundances were Hellinger-transformed to down-weight rare species (Legendre & Gallagher 2001). We tested for compositional differences among sites using marginal (Type III) PERMANOVA (vegan::adonis2) on Bray–Curtis dissimilarities with 999 permutations, including log(volume) and site as predictors. Multivariate dispersion homogeneity was assessed using PERMDISP (vegan::betadisper). Community structure was visualized using non-metric multidimensional scaling (NMDS; 2 dimensions) on a subset of 97 colonies with ≥10 CAFI individuals to ensure stable ordination positions; statistical tests (PERMANOVA, db-RDA) used the full n = 112. Species associated with ordination axes were identified using envfit permutation tests (p < 0.01, R² > 0.10). PERMANOVA results were validated across five distance metrics (Bray–Curtis on Hellinger and Wisconsin transformations, Jaccard on presence–absence, Gower, and Raup–Crick) and via 500 balanced site-subsampling iterations (Fig. S2).

To test whether composition shifts continuously along the coral size gradient, we used distance-based redundancy analysis (db-RDA; vegan::dbrda) with log(volume) as the constrained variable and site partialed out (Condition). Variance partitioning (vegan::varpart) decomposed community variation into volume-only, site-only, shared, and residual fractions. Robustness was verified by repeating the db-RDA on iterated-rarefied distances (100 draws, averaged). Nestedness along the size gradient was tested using NODF (Almeida-Neto et al. 2008) against 999 quasiswap null matrices.

### Co-occurrence patterns (Supplement)

Pairwise co-occurrence was tested using a volume-weighted Bernoulli null model (Stier et al. 2012). For each species, a logistic GLM predicted per-coral occurrence probability from log(volume) + site. For each of 10,000 null iterations, species presence was drawn independently from Bernoulli distributions with predicted probabilities. The standardized effect size (SES) of observed co-occurrence was computed as (observed − null mean) / null SD, with FDR correction across all pairs. Species with ≥10 occurrences were included. Intraspecific density patterns were tested with a multinomial allocation null model to assess whether conspecifics aggregate on fewer corals than expected (consistent with mating-pair formation; Stier et al. 2012). To assess whether co-occurrence patterns vary with coral size, we repeated the volume-weighted null model separately for three coral size classes (volume terciles), fitting logistic regressions within each class to obtain class-specific occurrence probabilities (Fig. S11C).

### Q3: CAFI diversity and coral condition (BEF framework)

We tested whether CAFI diversity predicts coral condition using a biodiversity–ecosystem function (BEF) framework (Tilman et al. 2014). The complementarity hypothesis predicts that more CAFI species provide non-redundant functional benefits (e.g., defense, cleaning, nutrient cycling), generating a positive diversity–condition relationship beyond what abundance alone predicts.

Forward models (CAFI → condition) took the form:

    condition_PC1 ~ CAFI_predictor + log(volume) + site

with ordinary least squares (OLS) standard errors as primary inference. Breusch–Pagan tests confirmed homoscedasticity for all models (BP p > 0.5); HC3 heteroscedasticity-consistent standard errors (sandwich package) are reported as a supplement sensitivity check (conservative at n < 100; Long & Ervin 2000). Count predictors were square-root transformed prior to fitting.

Predictors were organized into three testing tiers with distinct correction procedures:
- **A priori BEF** (k = 2; Hochberg FWER): Species richness and total CAFI abundance, predicted by BEF theory as the diversity and abundance pathways to ecosystem function.
- **Exploratory** (k = 4; Benjamini–Hochberg FDR): Shannon diversity, *Trapezia* abundance, resident fish abundance, and *Galeropsis monodonta* abundance.
- **Supplement composition** (k = 1; uncorrected): PC1_CAFI (community identity, a distinct question from BEF).

To disentangle diversity from abundance effects, we conducted three complementary analyses:
1. **Partial regression**: condition ~ richness + √(abundance) + log(volume) + site, with VIF diagnostics to assess collinearity.
2. **Variance partitioning**: Hierarchical R² decomposition quantifying variance uniquely attributable to richness, uniquely to abundance, and shared (confounded) between them.
3. **Path model**: A directed acyclic graph (piecewiseSEM) with volume → richness → condition and volume → abundance → condition pathways, using z-scored predictors for standardized coefficients.

To test whether the richness → condition relationship was an abundance artifact, we compared models using raw versus rarefied richness (expected species at n = 20; vegan::rarefy). However, rarefaction is an ambiguous test: it may remove either an abundance artifact or the BEF mechanism itself (diversity → abundance → condition), where the diversity-mediated abundance pathway is the core prediction of complementarity theory.

To test species-specific contributions to individual condition traits, we fitted linear models for each of the 19 most abundant species (present on ≥5 corals in the physiological dataset, n = 84): trait ~ √(species abundance) + log(volume) + site, where trait was one of five measures (protein, carbohydrate, zooxanthellae density, AFDW, or condition PC1). Standardized β coefficients were compiled in a heatmap (Fig. S16), and the nine strongest associations were visualized as scatter plots (Fig. S17). FDR correction was applied across all 95 species × trait tests.

Reverse models (condition → CAFI) used negative binomial or Poisson GLMs with condition_PC1 as predictor. Key species from the companion experimental study (up to 10 defined; 8 with sufficient survey occurrences after filtering n_present >= 5) were tested individually for condition effects, with Hochberg FWER correction. Cross-study concordance of effect directions was evaluated with a one-sided binomial test (H₀: 50% concordance by chance). A post hoc power analysis confirmed that with n = 84 corals, the survey had 85% power to detect the experiment's effect size (R² ≈ 0.12), ensuring that null results for exploratory predictors are credible rather than underpowered.

### Neighborhood effects (Supplement)

We tested whether neighborhood density affects CAFI abundance and coral condition using the 61 neighborhood-surveyed corals. Full models included log(volume), n_neighbors, total_neighbor_volume, mean_neighbor_dist, and site as predictors. Power analysis indicated 65% power for medium effects (R² = 0.10) at α = 0.05. Compositional variability across neighborhood density categories was assessed using PERMDISP (betadisper) on Bray–Curtis dissimilarities.

### Sensitivity analyses

Results were tested for robustness across five taxonomic resolution scenarios: baseline (243 OTUs), species-only (154 OTUs), merge-up (higher-level IDs merged to nearest named taxon), lump-down (all IDs resolved to genus), and rare-excluded (singletons removed). Seven metrics were evaluated per scenario: abundance β, richness z, PERMANOVA R², betadisper F, Shannon diversity, rarefied richness–condition relationship, and co-occurrence network modularity (Fig. S8). Spatial autocorrelation was assessed using Moran's I on CAFI abundance, richness, and Shannon diversity (Fig. S4).
