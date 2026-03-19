---
title: "Colony size drives sublinear scaling, compositional turnover, and biodiversity--condition feedbacks in coral-associated fauna"
bibliography: references.bib
csl: journal-of-animal-ecology.csl
link-citations: true
lang: en
---

# Colony size drives sublinear scaling, compositional turnover, and biodiversity--condition feedbacks in coral-associated fauna

Adrian C. Stier^1,2^*, Alexander Primo^3^, Joseph S. Curtis^1,4^, Craig W. Osenberg^3^

^1^Department of Ecology, Evolution, and Marine Biology, University of California Santa Barbara, Santa Barbara, CA 93106, USA

^2^Marine Science Institute, University of California Santa Barbara, Santa Barbara, CA 93106, USA

^3^Odum School of Ecology, University of Georgia, 140 E Green St, Athens, GA 30602, USA

^4^Department of Marine Science, University of Otago, PO Box 56, Dunedin, Otago 9054, New Zealand

*Corresponding author: Adrian C. Stier (astier@ucsb.edu)

---

## Abstract

The abundance of organisms in biogenic habitats scales predictably with habitat size, yet whether that scaling is proportional or diluted across larger patches carries consequences for the habitats themselves. Whether sublinear scaling generates compositional turnover, and whether that turnover feeds back on habitat condition, has not been tested. We surveyed 114 *Pocillopora* coral colonies spanning three orders of magnitude in size across three reef habitats in Mo'orea, French Polynesia, cataloguing 3,989 coral-associated fishes and invertebrates (CAFI) from 243 taxa. Total CAFI abundance scaled sublinearly with colony volume (β = 0.52 [0.44, 0.62]): larger corals harbour more occupants but at progressively lower per-capita densities, consistent with Propagule Redirection (density dilution as habitat competes for a finite colonist pool), contrasting with proportional scaling during initial experimental colonization that implicates post-settlement processes. Colony volume explained 7.8% of community composition after accounting for reef site (db-RDA, p = 0.001), with species turnover -- not nestedness -- along the size gradient. Species richness was positively associated with coral physiological condition (p = 0.018), consistent with biodiversity--ecosystem function theory, though the strong richness--abundance correlation (r = 0.84) and modest effect size prevent definitive causal attribution. These results trace a feedback from habitat quantity through community assembly to habitat condition -- a chain relevant wherever biogenic habitats are restructured by disturbance or restoration -- and suggest that density dilution may compound the vulnerability of large habitat-forming organisms to climate stress.

**Keywords:** biodiversity--ecosystem function, community assembly, coral-associated fauna, density dilution, *Pocillopora*, Propagule Redirection, species--area relationship, species turnover

---

## Introduction

The relationship between habitat quantity and occupant abundance is among the most general patterns in ecology [@ConnorMcCoy1979; @Preston1962; @Rosenzweig1995], yet its functional form -- whether occupant abundance increases proportionally with habitat or less than proportionally -- remains poorly resolved in most systems. In ecosystems structured by habitat-forming foundation species -- trees, corals, kelps, oysters -- this distinction carries consequences well beyond the occupants themselves [@Dayton1972; @Hanski1998; @MacArthurWilson1967]. Proportional scaling maintains constant per-unit density across the habitat landscape; sublinear scaling dilutes occupants in larger patches and concentrates them in smaller ones.

The functional form matters most when the habitat is biogenic and responsive to occupant density. Organisms residing in biogenic habitats alter the growth, survival, and condition of their hosts through feeding [@Silliman2005], nutrient processing [@MeyerSchultz1985], predator defense [@Glynn1983], and interactions with other residents [@Angelini2011]. Per-capita density therefore determines the intensity of these occupant--host interactions, linking habitat configuration to habitat condition. This link has been recognized theoretically [@Hamman2018; @Moeller2023] but tested in only a handful of systems.

The response of colonists to variation in biogenic habitat ranges from proportional tracking to active dilution. The Field of Dreams hypothesis proposes that increased habitat leads to proportionate increases in colonists, maintaining constant per-unit density [@Palmer1997; @ResetaritsBinckley2009]. The Propagule Redirection hypothesis posits that habitat patches compete for a shared colonist pool, so that increasing habitat dilutes occupants and reduces per-unit density [@CarrHixon1997; @Osenberg2002; @StierOsenberg2010]. Both regimes can operate simultaneously -- juveniles may redistribute among patches while adult abundances track habitat proportionally -- and the outcome depends on colonist longevity, density-dependent mortality, and migration frequency [@Keller2017; @Hamman2018].

These scaling regimes carry distinct consequences for habitat-forming organisms. Under Propagule Redirection, small or isolated patches accumulate disproportionately high occupant densities, enhancing condition if occupants provide benefits but suppressing it if they impose costs [@Bissett2025; @StierOsenberg2024a]. Under Field of Dreams, density is invariant across the landscape, decoupling habitat configuration from occupant-mediated effects. The distinction matters for predicting ecosystem trajectories, particularly as ongoing habitat loss restructures the occupant communities that influence recovery [@Beck2011; @Byrnes2011; @Kimbro2017; @Rennick2022].

If the scaling regime determines which species occupy a given patch, then community variation along the habitat size gradient could itself influence habitat condition. Biodiversity-ecosystem function (BEF) theory predicts that functionally diverse communities enhance ecosystem performance through complementarity -- species contribute non-redundantly to key processes [@Hooper2005; @LoreauHector2001]. Whether BEF principles extend to habitat-occupant systems, where diverse residents feed back on the biogenic habitat that supports them, has not been formally considered.

Branching corals in the genus *Pocillopora* offer a tractable model for linking habitat quantity, community assembly, and habitat condition. *Pocillopora* colonies create complex three-dimensional structures that support dense assemblages of coral-associated fishes and invertebrates (CAFI), sometimes exceeding 65 individuals from over 20 species on a single colony [@Curtis2023]. Guard crabs (*Trapezia* spp.) and snapping shrimp (*Alpheus* spp.) defend corals from predatory sea stars and remove sediment [@Glynn1983; @Stewart2006; @Stier2012], while coral-dwelling fishes subsidize coral nutrition through ammonia-rich excretion [@Chase2020; @Holbrook2008; @Shantz2023]. Corallivorous gastropods, by contrast, consume tissue and facilitate disease [@Boucher1986; @Nicolet2013]. Because many of these species are obligate symbionts that occupy corals as territorial mating pairs [@Castro1978; @HuberColes1986], intraspecific density on a single colony is effectively capped. As total occupancy increases on larger corals, growth comes predominantly from adding new species rather than more individuals of the same species -- a process termed "species stacking" [@Stier2012]. Species stacking generates a compositional prediction distinct from passive sampling: if intraspecific ceilings force diversification along the size gradient, small-coral faunas should not be nested subsets of large-coral faunas. Instead, colony size should structure the identity of occupants, not just their quantity.

CAFI abundance increases with coral size but less than proportionally [@Gotelli1985; @Curtis2023], a pattern consistent across coral genera and reef systems [@Harborne2011; @Brush2026] that implies density dilution under Propagule Redirection. Scaling slopes vary widely among taxa -- from 0.06 to 0.64 across decapod species [@Gotelli1985] -- and colony-scale variation can explain up to 84% of cryptofaunal richness variance [@Counsell2018]. A companion experiment that manipulated coral density on reef patches (1, 3, or 6 colonies per 25 m^2^) found proportional scaling during active colonization [@StierInReview], raising the question of whether post-settlement processes shift the scaling regime over time.

Three links in the habitat--community--condition chain remain untested in established communities: whether abundance and richness scale sublinearly with colony size, whether composition turns over or merely nests along the size gradient, and whether the resulting community variation covaries with coral condition independently of abundance. Here, we survey 114 *Pocillopora* colonies spanning a natural size gradient across three reef habitats on Mo'orea, French Polynesia. We catalogued 3,989 individuals from 243 taxa and address three questions:

**Q1 (Scaling):** How do CAFI abundance, richness, and per-capita density scale with coral colony volume? We predicted sublinear scaling consistent with Propagule Redirection and density dilution in larger colonies, and that the raw species--area relationship would be primarily a passive sampling artifact (the "more individuals" hypothesis) [@SrivastavaLawton1998; @Wright1983].

**Q2 (Composition):** Does CAFI community composition change along the coral size gradient through species turnover or nestedness? If species stacking drives assembly, small-coral faunas should not be nested subsets of large-coral faunas; instead, colony size should structure the identity -- not just the quantity -- of occupants.

**Q3 (Condition):** Does CAFI species richness predict coral physiological condition independently of abundance? We predicted that complementarity among functionally distinct residents -- guard crabs, nutrient-subsidizing fishes, corallivorous gastropods -- would generate a diversity--condition relationship beyond the richness--abundance correlation.

Together, these three questions address successive links in a hypothesized causal chain from habitat quantity through community assembly to habitat condition, extending the scaling patterns observed during experimental colonization [@StierInReview] to established communities on a natural size gradient. If occupant density and composition depend predictably on habitat size, then habitat loss and restoration restructure the communities that feed back on habitat condition -- a prediction as relevant to oyster reefs and forest canopies as to coral branches.

---

## Methods

### Study system and survey design

We surveyed 114 *Pocillopora* spp. coral colonies (morphotype identifications; not genetically confirmed) across three reef habitats on Mo'orea, French Polynesia (17°30'S, 149°50'W) during June--August 2019 (Fig. 1A). Sites were selected to span the range of reef environments: Hauru (HAU; fringing reef, north shore; n = 38), Maatea (MAT; lagoon/back reef, east shore; n = 39), and Maharepa (MRB; barrier reef, north shore; n = 37). At each site, colonies were selected non-randomly to span the full available size range, with divers targeting small, medium, and large colonies to ensure adequate coverage of the three-order-of-magnitude size gradient. Two MRB colonies lacked volume measurements and were excluded from size-dependent analyses (final n = 112; MRB n = 35).

### Coral measurements

Colony size was estimated as ellipsoidal volume: V = (4/3) × π × (L/2) × (W/2) × (H/2), where L, W, and H are the maximum length, width, and height (cm). While 3D photogrammetry provides more precise volume estimates [@Curtis2023], the ellipsoidal approximation is standard for *Pocillopora* colonies, which have roughly symmetrical branching morphology. This approximation overestimates true colony volume because it includes interstitial space between branches, but the bias is approximately proportional across colony sizes and does not affect log-linear scaling estimates. Colony volumes ranged from 21 to 42,333 cm³ (Fig. 1B). For a subset of 61 colonies ("neighborhood corals"), all *Pocillopora* neighbors within a 5-m radius were counted and their distances measured. The remaining 53 colonies were surveyed for CAFI and colony size only and are excluded from neighborhood analyses. Neighborhood density ranged from 1 to 76 neighbors (median = 17; Fig. 1C).

![](../output/figures/manuscript/fig1_study_design.png)

### CAFI collection and identification

All coral-associated fishes and invertebrates (CAFI) were extracted by enclosing each colony in a fine-mesh bag, carefully detaching the coral from the substrate, and preserving the entire colony and bag contents in 95% ethanol. No chemical anesthetics (e.g., clove oil) were used. Boring and deeply cryptic fauna embedded within the coral skeleton were not included. Specimens were identified to the lowest practical taxonomic level: 154 to species and 89 to genus, family, or higher, yielding 243 taxa from 3,989 total individuals. Species richness and Shannon diversity (H') were computed per colony using the full taxon set. Rarefied richness (expected species at n = 20 individuals) was computed using vegan::rarefy [@Oksanen2022] to control for the well-documented correlation between abundance and raw richness.

### Coral physiological condition

Tissue samples were collected from the upper branch of each colony for physiological assays of total protein (mg cm⁻²), carbohydrate (mg cm⁻²), zooxanthellae density (cells cm⁻²), and ash-free dry weight (mg cm⁻²). Position correction removed a tissue-gradient sampling artifact: sampling nubbins from smaller colonies captured a larger fraction of the branch axis, including more of the lower-density tip tissue, creating a systematic size bias in raw trait values. We regressed each trait on stump length and nubbin length and used the residuals. Coral condition was then quantified as the first principal component (PC1) of a principal component analysis (PCA) on the four position-corrected traits (61.6% of variance; all four traits loaded positively). Physiology data were available for 108 colonies; after merging with complete CAFI and volume data, 84 colonies were retained for BEF analyses.

### Statistical analyses

All analyses were conducted in R v4.4 (session information including package versions archived at [DOI to be assigned upon acceptance]). We used natural logarithms throughout; site was included as a fixed effect (three sites provide too few levels for reliable random intercept estimation; Bolker et al. 2009). We applied three levels of multiple-testing correction: (1) Holm correction for the two pre-specified BEF predictors (richness and total abundance; k = 2); (2) Benjamini--Hochberg FDR for four exploratory CAFI predictors (Shannon diversity, *Trapezia* abundance, fish abundance, and *Galeropsis* abundance; k = 4); and (3) Benjamini--Hochberg FDR across all species tested in species-level analyses. Model fit for GLMs was assessed using simulated residuals (DHARMa package), Cook's distance (threshold: 4/n), and variance inflation factors. All tests used α = 0.05. Spatial independence was confirmed using Moran's I on inverse-distance weights (Table S10). All code is available at [URL to be provided upon acceptance].

#### Q1: Scaling of CAFI with coral size

We tested how CAFI abundance scales with coral volume by fitting the power-law relationship N = a × V^β^ using a negative binomial generalized linear model (GLM; log link):

    total_cafi ~ log(volume) + site

where the coefficient on log(volume) directly estimates the scaling exponent β. We tested three hypotheses: Field of Dreams (β = 1; abundance scales proportionally with volume), Propagule Redirection (β < 1; sublinear scaling indicates density dilution; Stier & Osenberg 2010), and super-linear scaling (β > 1; larger corals disproportionately attractive). The null hypothesis (β = 1) was tested using both a Wald z-test (z = [β − 1]/SE) and a two-sided bootstrap proportion test (proportion of 1,000 site-stratified bootstrap replicates with β ≥ 1 or β ≤ 1, doubled). Bootstrap CIs are the primary inference; Wald z-tests are reported for comparability with prior studies. Confidence intervals were computed using bias-corrected and accelerated (BCa) bootstrap with site-stratified resampling (1,000 iterations; 1,000 replicates are adequate for percentile CIs, though conservative for BCa, and results were stable across replicate runs). Where BCa acceleration constants could not be estimated, percentile CIs were used as a fallback. Species richness was modeled analogously using a Poisson GLM. Rarefied richness was regressed on log(volume) using ordinary least squares. Per-capita CAFI density (individuals cm⁻³; total CAFI / colony volume) was computed per colony; the log--log slope of density on volume equals β − 1 by definition and was not independently modeled.

We repeated the scaling analysis for each of the 21 most prevalent individual species (≥30 total individuals and ≥15% prevalence) and for six taxonomic groups (*Trapezia* crabs, shrimps, other crabs, gastropods, fish, other invertebrates). An inverse-variance-weighted mean β across all species was computed as a fixed-effect summary (weights = 1/SE²); because species co-occur on the same corals, this assumes independence among species-level estimates and the SE of the weighted mean may be underestimated. FDR correction was applied within species-level and group-level test families. To characterize size-dependent occurrence, we fit logistic GLMs (presence ~ log(volume) + site) for 24 species with ≥15% prevalence and applied FDR correction across all tests (Fig. S12). (These 24 species met the ≥15% prevalence threshold; the 21-species scaling set additionally required ≥30 total individuals.)

#### Q2: Community composition

Community composition was analyzed on all 112 colonies with volume data. Species abundances were Hellinger-transformed (square root of relative abundances) to down-weight the influence of the most abundant species on dissimilarity calculations [@LegendreGallagher2001]. We tested for compositional differences among sites using marginal (Type III) permutational multivariate analysis of variance (PERMANOVA; vegan::adonis2) on Bray--Curtis dissimilarities with 999 permutations, including log(volume) and site as predictors. Multivariate dispersion homogeneity was assessed using permutational analysis of multivariate dispersions (PERMDISP; vegan::betadisper). Non-metric multidimensional scaling (NMDS) ordination (2 dimensions) was computed on all 112 colonies with volume data. All statistical tests (PERMANOVA, db-RDA) used the full dataset (n = 112). Species associated with site differences were identified by fitting species vectors to the NMDS ordination (envfit; 999 permutations; R² > 0.10, p < 0.01). Pairwise site comparisons used separate PERMANOVA tests with Bonferroni correction (k = 3). PERMANOVA results were validated across five distance metrics (Bray--Curtis on Hellinger and Wisconsin transformations, Jaccard on presence--absence, Gower, and Raup--Crick) and via 500 balanced site-subsampling iterations (Fig. S2).

To test whether composition shifts continuously along the coral size gradient, we used distance-based redundancy analysis (db-RDA; vegan::dbrda) with log(volume) as the constrained variable and site partialed out (Condition). The db-RDA F-statistic is algebraically identical to the marginal PERMANOVA F when both condition on the same covariates. Variance partitioning (vegan::varpart) decomposed community variation into volume-only, site-only, shared, and residual fractions. Species scores on the constrained axis were computed as weighted averages of site scores. Robustness was verified by repeating the db-RDA on iterated-rarefied distances: for each of 100 iterations, communities were rarefied to the minimum observed abundance among colonies with ≥5 individuals, a Bray--Curtis distance matrix was computed, and the resulting matrices were averaged element-wise before ordination and hypothesis testing. Nestedness along the size gradient was tested using NODF (nestedness metric based on overlap and decreasing fill; Almeida-Neto et al. 2008) against 999 quasiswap null matrices (an algorithm preserving row and column totals; Miklós & Podani 2004).

#### Co-occurrence patterns (Supplement)

Pairwise co-occurrence was tested using a volume-weighted Bernoulli null model [@Stier2012]. For each species, a logistic GLM predicted per-coral occurrence probability from log(volume) + site. For each of 10,000 null iterations, species presence was drawn independently from Bernoulli distributions with predicted probabilities. The standardized effect size (SES) of observed co-occurrence was computed as (observed − null mean) / null SD, with false discovery rate (FDR) correction across all pairs. Species with ≥10 occurrences were included. Intraspecific density patterns were tested with a multinomial allocation null model to assess whether conspecifics aggregate on fewer corals than expected (consistent with mating-pair formation; Stier et al. 2012). To assess whether co-occurrence patterns vary with coral size, we repeated the volume-weighted null model separately for three coral size classes (volume terciles), fitting logistic regressions within each class to obtain class-specific occurrence probabilities (Fig. S9C).

#### Q3: CAFI diversity and coral condition (BEF framework)

We tested whether CAFI diversity predicts coral condition using a biodiversity--ecosystem function (BEF) framework [@Tilman2014]. The complementarity hypothesis predicts that more CAFI species provide non-redundant functional benefits (e.g., defense, cleaning, nutrient cycling), generating a positive diversity--condition relationship beyond what abundance alone predicts.

Forward models (CAFI → condition) took the form:

    condition_PC1 ~ CAFI_predictor + log(volume) + site

with ordinary least squares (OLS) standard errors as primary inference. Breusch--Pagan tests confirmed homoscedasticity for all models (BP p > 0.5); HC3 heteroscedasticity-consistent standard errors (sandwich package) are reported as a supplement sensitivity check (conservative at n < 100; Long & Ervin 2000). Count-based CAFI predictors (total abundance, *Trapezia* count, fish count, *Galeropsis* count) were square-root transformed prior to fitting to reduce right skew. Species richness was not transformed because it has a more linear relationship with condition and lower skewness than raw abundance counts.

Community composition (PC1_CAFI) was tested separately in the supplement.

To disentangle diversity from abundance effects, we conducted three complementary analyses:
1. **Partial regression**: condition ~ richness + √(abundance) + log(volume) + site, with variance inflation factor (VIF) diagnostics to assess collinearity.
2. **Variance partitioning**: Hierarchical R² decomposition quantifying variance uniquely attributable to richness, uniquely to abundance, and shared (confounded) between them.
3. **Path model**: A piecewise structural equation model (piecewiseSEM; Lefcheck 2016) with three component OLS models: abundance ~ log(volume) + site, richness ~ log(volume) + site, and condition ~ richness + abundance + log(volume) + site. Predictors were z-scored for standardized path coefficients. Model fit was evaluated using Fisher's C statistic.

To test whether the richness → condition relationship was an abundance artifact, we compared models using raw versus rarefied richness (expected species at n = 20; vegan::rarefy). Rarefied richness was undefined for colonies with fewer than 20 CAFI individuals, reducing the available sample from 84 to 47 for this comparison. We address the interpretive ambiguity of this test in the Discussion.

To test species-specific contributions to individual condition traits, we fitted linear models for each of the 19 most abundant species (present on ≥5 corals in the physiological dataset, n = 84): trait ~ √(species abundance) + log(volume) + site, where trait was one of five measures (protein, carbohydrate, zooxanthellae density, AFDW, or condition PC1). Standardized β coefficients were compiled in a heatmap (Fig. S13), and the nine strongest associations were visualized as scatter plots (Fig. S14). FDR correction was applied across all 95 species × trait tests.

Reverse models (condition → CAFI) tested seven response variables: total CAFI abundance, species richness, Shannon diversity, *Trapezia* abundance, fish abundance, *Galeropsis monodonta* abundance, and community composition (PC1_CAFI). Count responses used negative binomial or Poisson GLMs; continuous responses used OLS. BH-FDR correction was applied across all seven tests. We tested the 8 key species from the companion experiment that had ≥5 occurrences in our physiological dataset (of 10 originally defined) individually for condition effects, with FDR correction. Cross-study concordance of effect directions was evaluated with a one-sided binomial test (H₀: 50% concordance by chance). A sensitivity power analysis (minimum detectable R² increment at 80% power, α = 0.05, n = 84) assessed whether the survey sample was sufficient to detect the companion experiment's observed effect sizes.

#### Neighborhood effects (Supplement)

We tested whether neighborhood density affects CAFI abundance and coral condition using the 61 neighborhood-surveyed corals. Full models included log(volume), n_neighbors, total_neighbor_volume, mean_neighbor_dist, and site as predictors. Power analysis indicated 65% power for medium effects (R² = 0.10) at α = 0.05. Compositional variability across neighborhood density categories was assessed using PERMDISP (betadisper) on Bray--Curtis dissimilarities.

#### Sensitivity analyses

Results were tested for robustness across five taxonomic resolution scenarios: baseline (243 taxa), species-only (154 taxa), merge-up (higher-level IDs merged to nearest named taxon), lump-down (all IDs resolved to genus), and rare-excluded (singletons removed). Six metrics were evaluated per scenario: abundance scaling exponent (β), richness scaling exponent (z), Shannon scaling slope, PERMANOVA R² (site and volume, counted separately), and rarefied richness--condition slope (Fig. S6; Table S9). Spatial autocorrelation was assessed using Moran's I on CAFI abundance, richness, and Shannon diversity (Table S10).

#### Ethical statement

Research was conducted under permits issued by the Haut-commissariat de la République en Polynésie française and in collaboration with the Université de la Polynésie française (CRIOBE). All protocols followed institutional and local regulations for scientific collection on the reefs of Mo'orea. Coral collections were conducted as part of a broader research programme at the UC Berkeley Richard B. Gump South Pacific Research Station, Invertebrate collections did not require IACUC review; coral tissue sampling protocols were approved under broader programme oversight at UCSB.

#### Data accessibility

Data and analysis code are available at [DOI to be assigned upon acceptance]. Data and code will be archived in the Dryad Digital Repository upon acceptance.

---

## Results

### Q1: CAFI abundance and richness scale sublinearly with coral size

We found that total CAFI abundance scaled sublinearly with colony volume, consistent with Propagule Redirection (Fig. 2A). The negative binomial GLM estimated a scaling exponent of β = 0.52 (95% bootstrap CI [0.44, 0.62]; Wald z vs. β = 1: z = −11.45, p < 0.001; bootstrap p < 0.001), decisively rejecting the Field of Dreams hypothesis (β = 1). Species richness also scaled sublinearly (Poisson GLM: z = 0.34 [0.27, 0.42]; z vs. 1 = −18.76, p < 0.001; Fig. 2C). However, rarefied richness (expected species at n = 20 individuals) showed no relationship with volume (OLS: slope = 0.14, SE = 0.30, t₆₄ = 0.47, p = 0.64, n = 68; colonies with <20 CAFI excluded). We return to the interpretation of this pattern in the Discussion.

Per-capita CAFI density declined with colony size (log--log slope = −0.48; Fig. 2B). Because density equals abundance divided by volume, this slope is the algebraic complement of the abundance exponent (0.52 − 1 = −0.48), illustrating that sublinear scaling necessarily produces density dilution -- larger corals harbor fewer CAFI per unit volume.

![Figure 2. CAFI abundance and richness scale sublinearly with coral volume.](../output/figures/manuscript/fig2_scaling.png)

We detected sublinear scaling at the species level: 11 of 21 prevalent species (52%) showed Redirection (bootstrap CI excluding 1.0), 10 showed scaling consistent with Field of Dreams (CI spanning 1.0), and none showed super-linear scaling (Fig. 3A,C; Table S1). The inverse-variance-weighted mean β across 21 species was 0.51 [0.45, 0.56], significantly below 1.0 (z = −18.8, p < 0.001; Fig. S4). This average treats species as independent estimates, although species co-occurring on the same corals are not strictly independent. Obligate symbionts -- *Trapezia* crabs (β = 0.43) and snapping shrimp *Alpheus lottini* (β = 0.37) -- consistently showed sublinear scaling, while facultative associates such as *Caracanthus maculatus* (β = 1.18, CI overlapping 1.0) and the coralliophiline snail *Galeropsis monodonta* (β = 1.27, CI overlapping 1.0) scaled proportionally.

Among the six taxonomic groups, five showed significant sublinear scaling: *Trapezia* crabs (β = 0.43 [0.35, 0.52], p < 0.001), shrimps (β = 0.50 [0.40, 0.60], p < 0.001), other crabs (β = 0.47 [0.24, 0.71], p < 0.001), fish (β = 0.74 [0.58, 0.92], p = 0.004), and other invertebrates (β = 0.50 [0.33, 0.67], p < 0.001). Gastropods were the sole exception (β = 0.94 [0.73, 1.17], p vs. 1 = 0.60), scaling proportionally (Fig. 3B,D; Fig. S7). All results were robust to taxonomic resolution: abundance β varied by <0.01 across five sensitivity scenarios, and richness z remained significantly sublinear in all cases (Fig. S6).

![Figure 3. Scaling of individual species and taxonomic groups with coral volume.](../output/figures/manuscript/fig3_species_group_scaling.png)

Species occurrence probability was size-dependent for 14 of 24 prevalent species (logistic GLM, FDR < 0.05; Fig. S12), with most showing increasing occurrence probability with coral size. Size-dependent scaling is thus common across the community, though 10 species showed no significant size dependence. (Spatial autocorrelation was negligible for all metrics: Moran's I = −0.004–0.024, all p > 0.28; Table S10.)

### Q2: Site pools and coral size structure community composition

Community composition differed significantly among the three reef sites (marginal PERMANOVA on Bray--Curtis dissimilarity of Hellinger-transformed abundances: volume F₁,₁₀₈ = 9.74, R² = 0.08, p = 0.001; site F₂,₁₀₈ = 3.74, R² = 0.06, p = 0.001; n = 112; dispersion homogeneity confirmed; PERMDISP F = 0.89, p = 0.42; sampling adequacy confirmed, Fig. S1). Together, site and volume explained ~14% of compositional variation, with 86% unexplained by the measured predictors. All pairwise site comparisons were significant after Bonferroni correction (all p_adj ≤ 0.006). The site effect was robust across all five distance metrics tested (Fig. S2). Balanced site subsampling (n = 35 per site, 500 iterations) confirmed the site effect was robust, with 100% of iterations significant at p < 0.05 (mean p = 0.001).

The three sites exhibited distinct taxonomic signatures (Fig. 4B). Maatea was characterized by high hermit crab abundance (33% of all CAFI, primarily *Calcinus latens*). Maharepa was dominated by obligate coral symbionts, with shrimps and crabs together comprising 71% of the assemblage. Hauru supported the most taxonomically balanced community and the highest proportion of coral-dwelling fishes (12%). Seven species drove the largest compositional differences among sites (envfit R² > 0.10, p < 0.01), including *Paragobiodon modestus* (R² = 0.46), *Fennera* spp. (R² = 0.39), *Calcinus latens* (R² = 0.39), and *Alpheus lottini* (R² = 0.27).

We also found that community composition shifted continuously along the coral size gradient. Distance-based redundancy analysis (db-RDA; Fig. 4A), after partialing out site effects, confirmed that coral volume explained 7.8% of compositional variation (p = 0.001). This size--composition gradient was robust to rarefaction (2.4%, F₁,₁₀₈ = 2.64, p = 0.001), confirming genuine compositional turnover rather than an abundance artifact.

![Figure 4. CAFI community composition differs among reef sites and along the coral size gradient.](../output/figures/manuscript/fig4_composition.png)

Species with the strongest loadings on the db-RDA size axis included *Trapezia punctimanus* (score = −2.95) and *Luniella pugil* (−1.71), both associated with smaller corals, while *Euplica varians* and *Trapezia flavopunctata* (scores = 1.06) were most strongly associated with larger corals. Variance partitioning (adjusted R²) attributed 4.8% to volume alone, 3.9% shared between volume and site, with negligible unique site variation (−0.2%) and 91.5% residual. Compositional divergence among size classes was not significant after controlling for abundance through rarefaction (betadisper: p = 0.61; Fig. S3).

Community nestedness along the size gradient was not significant (NODF = 18.4, z = −1.09, p = 0.28). NODF ranges from 0 (no nestedness) to 100 (perfect nestedness); the observed value of 18.4 indicates low nestedness, confirming that small-coral faunas are not subsets of large-coral faunas. Combined with the significant db-RDA, this pattern is consistent with species turnover -- not passive accumulation -- along the size gradient: some species are gained while others are lost as coral size increases. This non-nested pattern is consistent with the species-stacking constraint: intraspecific density ceilings drive diversification along the size gradient rather than passive accumulation.


### Q3: CAFI species richness positively predicts coral condition

We found that corals harboring more CAFI species were in better physiological condition (n = 84 colonies with physiological data). Species richness positively predicted coral condition (standardized β = 0.27, t₇₉ = 2.42, p = 0.018, Hochberg-corrected for two a priori tests; Fig. 5A; Fig. S11A), and total CAFI abundance showed a similar but weaker association (standardized β = 0.32, t₇₉ = 2.01, p = 0.048; Fig. 5B). Variance partitioning indicated that the richness signal was largely independent of abundance: of the incremental R² = 0.07 explained by CAFI predictors beyond volume and site, 29.1% was uniquely attributable to richness (~2% of total variance), less than 1% uniquely to abundance (<0.1% of total variance), and 70.8% was shared between them (~5% of total variance; Fig. S10A). The richness-unique fraction, while small in absolute terms, was roughly 30x larger than the abundance-unique fraction.

The overall model (richness + volume + site) explained modest variance (adjusted R² = 0.04). Because richness and abundance were strongly correlated (r = 0.84), partial regression including both predictors produced high collinearity (VIF = 6.2 for richness, 6.8 for abundance), complicating attribution of variance to diversity per se versus abundance.

In the path model, the standardized coefficient from richness to condition (β = 0.55) was larger in magnitude than from abundance to condition (β = 0.02), though neither path reached significance individually (p = 0.20 and p > 0.50, respectively; Fig. S10C). The path model did not fit significantly better than the null (Fisher's C p = 0.20), so these coefficients should be interpreted cautiously.

Rarefied richness (expected species at n = 20 individuals) showed no relationship with condition (β = −0.07, p = 0.50, n = 47; colonies with fewer than 20 CAFI excluded from rarefaction; Fig. S8; Fig. S11B). We address the interpretive ambiguity of this test in the Discussion.

![Figure 5. CAFI diversity and abundance as predictors of coral physiological condition (BEF framework).](../output/figures/manuscript/fig5_feedbacks.png)

We found that no exploratory predictor -- Shannon diversity, *Trapezia* abundance, resident fish abundance, or *Galeropsis monodonta* abundance -- survived BH-FDR correction (all p_FDR > 0.80; Fig. S11C), though effect directions mostly matched expectations (Trapezia positive, fish positive). The exception was *Galeropsis monodonta*, which showed an unexpected positive association with condition despite being a known corallivore (Fig. S11D--E). Community composition (PC1_CAFI) did not predict condition (p > 0.10; tested separately in the supplement).

The reverse direction (condition → CAFI) showed one nominally significant result (condition → total CAFI abundance, raw p = 0.037), but this did not survive BH-FDR correction across the seven reverse-direction tests (Fig. S11F). Key species from the companion experimental study (8 of 10 with sufficient survey occurrences) were tested individually; no species was individually significant after FDR correction (Table S4). The survey had 80% power to detect a small-to-medium effect size (R² ≈ 0.09), indicating that null results for individual species are credible rather than underpowered.

Species-level analysis of individual condition traits (protein, carbohydrate, zooxanthellae, AFDW) revealed no significant associations after FDR correction across 95 species × trait tests (Fig. S13), suggesting that no single species drives the richness--condition relationship. The nine strongest individual associations (Fig. S14) showed effect directions consistent with mutualism hypotheses (positive for guard crabs and shrimps). In the neighborhood/landscape models, coral volume, neighborhood density, and site did not predict condition (all p > 0.30; Fig. S5).

### Supporting results (Supplement)

#### Neighborhood density does not predict CAFI abundance

Neighborhood analyses were restricted to the 61 colonies with 5-m survey data (Fig. S5; Table S5). Neighborhood density (n_neighbors) did not significantly predict CAFI abundance (β = −0.005, z = −0.89, p = 0.37) or species richness (β = −0.003, z = −0.88, p = 0.38) after controlling for coral volume and site. Mean neighbor distance showed a significant negative association with richness (β = −0.005, z = −3.20, p = 0.001) and Shannon diversity (β = −0.007, t = −3.49, p = 0.001), suggesting that corals closer to neighbors support slightly higher diversity, though this pattern was not observed for total abundance. No size × neighborhood interactions were significant (all p > 0.10). The power analysis indicated 65% power for medium effects (R² = 0.10), so null results for density may reflect insufficient power rather than a true absence of neighborhood effects.

#### Pairwise co-occurrences are explained by volume and site

We detected no significant pairwise co-occurrences using volume-weighted null models after FDR correction (0 of 528 pairs at FDR < 0.05; Fig. S9; Table S6). The two strongest signals were *Harpiliopsis beaupresii*--*Paragobiodon modestus* (SES = −3.43, p_FDR = 0.11) and *H. beaupresii*--*H. spinigera* (SES = −3.42, p_FDR = 0.21), both negative but non-significant after correction, indicating that pairwise associations are largely explained by volume and site. Intraspecific density analysis identified six species with significant mating-pair aggregation (FDR < 0.05): *Synalpheus charon* (SES = 10.3), *Alpheus lottini* (SES = 10.3), *Harpiliopsis beaupresii* (SES = 5.4), *Trapezia serenei* (SES = 4.4), *Alpheus pachychirus* (SES = 3.9), and *Trapezia bidentata* (SES = 3.8; Table S7).

#### Survey and experimental scaling exponents are broadly concordant

Seven species had scaling data in both the survey and companion experiment. Of these, two showed sublinear scaling in the survey but proportional scaling during experimental colonization (*Calcinus latens*, *Alpheus lottini*), one was sublinear in both studies (*Trapezia serenei*), and four scaled proportionally in both (*Caracanthus maculatus*, *Harpiliopsis spinigera*, *Galeropsis monodonta*, *Periclimenes watamuae*; Table S8).

---

## Discussion

CAFI communities on *Pocillopora* trace a partial feedback loop: colony size shapes community diversity through sublinear occupant scaling, and that diversity covaries -- modestly but consistently -- with coral physiological condition. Abundance scales at half-power with colony volume, density dilutes as corals grow, composition turns over (without nesting) along the size gradient, and the resulting community variation covaries with coral condition (Figs. 2, 4, 5). This loop is not empirically closed -- our cross-sectional design demonstrates associations at each link but cannot establish causality for any of them. We interpret our findings alongside the companion density-manipulation experiment [@StierInReview] and the broader experimental literature.

### Sublinear scaling and species stacking

Sublinear scaling of occupant abundance with habitat size recurs across biogenic habitat systems -- arthropods on trees [@Lawton1983; @Schlinkert2015], invertebrates in kelp holdfasts [@Anderson2005], parasites on hosts [@Kuris2008] -- pointing to a common constraint: the supply of colonists is finite relative to the habitat competing for them. Our community-level exponent (β = 0.52; Fig. 2A) is comparable to Abele and Patton's (1976) foundational estimate from *Pocillopora damicornis* (β = 0.62, area-based). Converting their exponent to a volume basis under roughly isometric growth yields ~0.41 -- directionally consistent with our estimate, though the conversion is approximate because *Pocillopora* colonies become proportionally flatter as they grow. The pattern was robust across taxonomic resolution scenarios (Fig. S6).

In *Pocillopora*, at least three mechanisms reinforce this constraint. Chemical settlement cues emanate from the colony surface and thus scale with surface area (proportional to radius squared) rather than volume (proportional to radius cubed). This geometric mismatch predicts sublinear volume-based scaling even if per-unit-surface-area colonization rates are constant, because signal halos whose interception area grows more slowly than habitat capacity [@Hamman2018; @StierOsenberg2010]. Obligate symbionts impose per-colony ceilings through territorial pair-bonding [@Castro1978; @HuberColes1986], and these protection mutualists reshape the broader assemblage -- *Trapezia* and *Alpheus* presence explains 20% of composition variation during colonization [@CounsellDonahue2021]. Increasing total CAFI on a colony therefore favours adding new species over additional conspecifics -- the "species stacking" constraint [@Stier2012] that links density dilution to community-level diversity.

The taxon-specificity of scaling exponents reinforces this mechanistic picture (Fig. 3). Obligate pair-forming symbionts show the most strongly sublinear scaling because pair-bonding imposes a hard density ceiling regardless of colony size. Gastropods, by contrast, scale near-proportionally, suggesting that *Galeropsis monodonta* colonises primarily by crawling between adjacent corals instead of pelagic larval delivery -- though whether *Galeropsis* recruitment is primarily crawl-mediated or larval remains untested. The one functional group whose colonisation biology plausibly bypasses the mechanism behaves accordingly. Analogous taxon-specificity appears in plant--herbivore and host--parasite systems, where specialist or directly-transmitted species track habitat size more tightly than generalists or trophically-transmitted ones [@Otway2005; @Poulin2007].

The companion experiment found proportional scaling during active colonization, while our survey of established communities shows sublinear scaling. This temporal shift suggests that density-dependent territorial exclusion progressively reduces per-capita survival as communities mature [@MacArthurWilson1967]. Several post-settlement processes could drive the transition: priority effects by established residents, density-dependent mortality at high occupancy, and ontogenetic changes in coral architecture that alter niche availability as colonies grow. These filters progressively reduce per-capita colonization success on larger, more densely occupied colonies. The two species that shifted from proportional to sublinear -- *Alpheus lottini* and *Calcinus latens* -- are both territorial, which is the expected pattern if pair-bonding drives the transition. However, only 2 of 7 overlapping species showed the predicted shift, so this remains a hypothesis awaiting broader confirmation. If the pattern holds, it has practical implications: initial CAFI recruitment to restoration outplants will attenuate as communities mature, and planning that assumes proportional returns will overestimate long-term gains.

The relative importance of Field of Dreams versus Propagule Redirection depends on ambient settlement intensity [@StierOsenberg2010]. At low settlement, the propagule pool is exhausted quickly and Redirection dominates; at high settlement, the system saturates and Field of Dreams takes over. The strong sublinear scaling we observe places Mo'orea in the Redirection-dominated regime, yielding a testable prediction: reefs with higher larval supply should show exponents closer to 1.0.

### Composition, co-occurrence, and the BEF question

Small-coral communities are not impoverished subsets of large-coral communities: the rejection of nestedness confirms that species turn over along the size gradient, generating qualitatively different assemblages on different-sized colonies. This parallels patterns in terrestrial fragments [@Haddad2015] and ant--myrmecophyte systems [@FonsecaGanade1996] and implies that colony size determines not just the quantity of occupant services but which services are available. A small *Pocillopora* with high per-capita densities of *Trapezia* pairs receives active predator defence and sediment removal; a large colony with lower densities of those obligate mutualists but higher representation of gastropods and xanthid crabs receives a different -- and potentially less favourable -- balance of effects. Size-dependent composition likely reflects both biotic and abiotic gradients: smaller colonies offer fewer microhabitat niches and simpler branching architecture, while larger colonies provide greater interstitial space, reduced flow in the colony interior, and lower predation risk in the canopy. Depth partitioning among *Trapezia* species along the colony size axis has been documented in *Pocillopora* [@Counsell2018]. This turnover survived rarefaction control (Fig. 4A), confirming genuine species replacement independent of passive sampling.

Yet site and volume together explain only ~14% of composition. The remaining 86% points to stochastic colonisation of discrete habitat patches [@Laurance2002; @MetaxasScheibling1993]. Volume-weighted null models detected no significant pairwise associations after FDR correction (Fig. S9), while six species showed significant mating-pair aggregation, concentrating conspecific pairs on a subset of corals and leaving others unoccupied. This non-random mosaic of well-defended and poorly-defended colonies -- arising from pair-aggregation within a broadly stochastic backdrop -- creates exactly the kind of among-colony variation in functional composition that a diffuse BEF mechanism would act on.

Whether CAFI diversity enhances coral condition or merely correlates with it remains ambiguous: richness predicted condition (Fig. 5A), but the BEF signal is modest and entangled with abundance. Most explanatory power was shared between richness and abundance (70.8%), with little uniquely attributable to either, and the overall effect was small (adjusted R² = 0.04), leaving open the possibility that the richness signal is an abundance artifact. The rarefied richness null result (Fig. S8) either confirms this interpretation or reflects that rarefaction strips out the mechanism through which diversity operates -- complementary effects require the additional individuals that diverse assemblages generate. Several lines of evidence favour complementarity as the dominant pathway: variance partitioning attributes 30x more unique explanatory power to richness than to abundance (Fig. S10A), the path model estimates a larger richness-to-condition coefficient (Fig. S10C), and no individual species survived FDR correction, implicating many species contributing small, additive effects [@Hector1999] instead of a single dominant mutualist. An equally parsimonious interpretation, however, is that richness is simply a proxy for total occupancy and the condition signal reflects the abundance pathway alone.

Our observational design cannot resolve this ambiguity. A further alternative is that richness and condition both respond to unmeasured drivers -- colony age, site-specific larval supply, or microhabitat complexity -- without any direct richness-to-condition effect. The unique richness contribution (~2% of total variance) is small enough that even a well-powered study would struggle to detect it reliably in partial regression -- our sample of 84 colonies had 80% power to detect R² ≈ 0.09 but not effects this small. The experimental literature provides mechanistic plausibility for complementarity -- defence [@Glynn1983; @Stier2012], nutrient subsidies [@Holbrook2008], bleaching buffering [@Chase2018], synergistic defence [@McKeonMoore2014; @McKeon2012], and a growing catalogue of facultative mutualists [@Clements2024; @Honeycutt2023; @SpadaroButler2021] -- but plausibility is not evidence. A manipulative experiment that holds abundance constant while varying richness would be required to establish the causal link.

The null result for neighborhood density (n = 61, 65% power for medium effects) likely reflects insufficient power. The significant negative association between mean neighbor distance and richness (p = 0.001) suggests spatial proximity to neighbors may matter more than raw count, pointing toward a source-sink dynamic that warrants replication at larger sample sizes. This result parallels the experimental finding that coral density affects CAFI colonization [@StierInReview]: the survey's smaller sample and coarser metric may lack the resolution to detect the same signal.

One unexpected result was *Galeropsis monodonta's* positive association with condition despite being a tissue-feeding corallivore. This likely reflects reverse causation -- healthier corals sustaining larger corallivore populations -- not a mutualistic effect. More broadly, mutualisms are context-dependent [@Bronstein1994; @Palmer2008; @StierOsenberg2024a], and warming can shift coral--CAFI relationships from mutualistic toward antagonistic as obligate symbionts decline [@Stella2014] and vermetids impose greater costs [@Shima2010; @Zill2017].

### Colony size, heatwaves, and the feedback loop

Density dilution has a physiological consequence: the mutualist services that coral colonies receive per unit volume diminish as they grow. The three questions we address form a hypothesized feedback in which coral size structures CAFI communities (Q1), whose composition turns over along the size gradient (Q2), and whose diversity may feed back on coral condition (Q3). This feedback weakens as corals grow -- the sublinear exponent (β = 0.52) means that a doubling of colony volume increases CAFI abundance by only ~43% -- and species turnover means the functional composition changes along the size gradient: small and large corals receive different services from partially non-overlapping species pools.

Density dilution may compound the well-documented size-dependent mortality of corals during marine heatwaves [@Speare2022]. Size-selective mortality is typically attributed to thermal physiology -- larger colonies have lower surface-area-to-volume ratios, reduced internal circulation, and greater accumulated bleaching susceptibility. Bleaching can itself disrupt CAFI mutualisms: @Stella2011 documented declines in crab density on bleached corals, and prolonged bleaching triggers emigration of obligate symbionts [@Stella2014], potentially flipping the interaction from mutualistic to neutral or antagonistic. We hypothesise that density dilution represents an additional, previously unrecognised contributing factor to size-dependent mortality. Large colonies receive mutualist services -- nutrient subsidies from damselfishes [@Chase2020; @Shantz2023], sediment removal by *Trapezia* [@Stewart2006; @Stier2012] -- diminishing to ~72% of previous intensity per doubling, a direct arithmetic consequence of the density-dilution exponent. The lower per-capita densities of beneficial CAFI on large corals could reduce biotic insurance under thermal stress. This hypothesis remains untested: it requires demonstrating that large colonies retaining high CAFI densities for their size class survive heatwaves at higher rates than those with depleted communities.

Isolated or small outplants should attract disproportionately high per-capita CAFI densities -- a prediction that restoration programmes rarely incorporate. Outplanting programmes optimise colony density for growth and survival but seldom consider how spatial design shapes the CAFI communities that colonise outplants. Whether elevated per-capita densities enhance or degrade outplant performance depends on which species dominate the local colonist pool and requires the species-specific "service catalogue" that Stier and Osenberg (2024b) advocate.

### Future directions

The central unresolved question is whether CAFI diversity causally improves coral condition or merely tracks it. Three lines of work would address this. First, a manipulative BEF experiment assembling CAFI communities of controlled richness (e.g. 1, 3, 6, 10 species) on standardised colonies, holding total abundance constant, would directly test whether diversity enhances condition independent of abundance. The species-stacking constraint provides a natural design: each added species occupies a distinct niche, making richness manipulations ecologically realistic. Second, longitudinal tracking of marked colonies with annual CAFI censuses before and after thermal stress events would test whether the diversity--condition association persists under warming or reverses as obligate symbionts decline (Stella et al. 2011, 2014). Colonies that naturally vary in CAFI density for their size class would provide a within-size-class test of the biotic insurance hypothesis. Third, comparative scaling analyses across coral morphologies (branching vs. massive), biogeographic regions, and settlement intensities would test whether β ≈ 0.5 is a conserved community property or varies predictably with ecological context -- with three-dimensional photogrammetry [@Curtis2023] testing whether structural complexity is a stronger predictor of scaling than bulk volume.

Whether the occupants are ants in acacia thorns, invertebrates in kelp holdfasts, or crustaceans in coral branches, the architecture is the same: habitat quantity shapes community assembly through predictable scaling rules, and the resulting communities have the potential to feed back on habitat condition. As biogenic habitats are restructured by disturbance and restoration, quantifying these feedbacks -- not just documenting the associations -- becomes essential for predicting what comes next.

---

## Acknowledgements

We thank the many people who contributed to establishing and maintaining the field surveys, as well as M. Brzezinski and D. Cryan for assistance with data management and laboratory work and C. Wall, C. Bove, and J. Baumann for guidance on coral physiological protocols. We are grateful to the staff of the University of California Gump Research Station for their logistical and technical support, and to the Moorea Coral Reef LTER, UCSB Ocean Recoveries Lab, and UGA Osenberg Lab for their insights. Research was conducted under permits issued by the Territorial Government of French Polynesia (Délégation à la Recherche) and the Haut-commissariat de la République en Polynésie Française (DTRT) (Protocole d'Accueil), whose continued support we greatly appreciate. This work was supported by the National Science Foundation (OCE-1851510 and OCE-1851032).

---

## Author Contributions

ACS and CWO conceived the project. AP and JSC conducted the majority of the fieldwork, laboratory, and image analyses. ACS performed the analyses and drafted the figures with input from CWO. ACS wrote the first draft, and ACS and CWO contributed substantially to revisions.

---

## Conflict of Interest

The authors declare no conflicts of interest.

---

## Figure Legends

**Figure 1.** Study design and representative coral-associated fauna (CAFI). (A) Satellite imagery of Mo'orea, French Polynesia (17°30'S, 149°50'W), showing the three reef sites: Hauru (fringing reef, north shore; n = 38 corals), Maatea (lagoon/back reef, east shore; n = 39), and Maharepa (barrier reef, north shore; n = 35; 2 additional colonies lacked volume data, total site n = 37). Site markers colored by reef site (purple = Hauru, slate = Maatea, sage green = Maharepa). (B) Distribution of *Pocillopora* colony volumes on a log₁₀ scale, spanning more than three orders of magnitude (21--42,333 cm³; n = 112 colonies). Black curve shows kernel density estimate. (C) Distribution of neighborhood density (number of *Pocillopora* colonies within 5 m; n = 61 corals with neighborhood surveys). (D--F) Representative CAFI species from three major functional groups: (D) *Trapezia* sp. guard crab (Trapeziidae), an obligate symbiont that defends corals from predators; (E) *Harpiliopsis spinigera* (Palaemonidae), an obligate *Pocillopora* symbiont; (F) *Neocirrhites armatus* (flame hawkfish), a coral-dwelling reef fish. All photographs from Mo'orea field collections (2019).

![Figure 1. Study design and representative CAFI.](../output/figures/manuscript/fig1_study_design.png)

**Figure 2.** CAFI abundance and richness scale sublinearly with coral volume. (A) Total CAFI abundance versus colony volume (log--log scale). Solid curve: negative binomial GLM fit (β = 0.52 [0.44, 0.62]); shaded band: 95% CI. Dashed line: proportional scaling (β = 1, Field of Dreams hypothesis). (B) Per-capita CAFI density (individuals per cm³) versus colony volume (log--log scale). Dashed horizontal line: mean density. Log--log slope = −0.48, confirming density dilution in larger corals. (C) Species richness versus colony volume. Solid curve: Poisson GLM (z = 0.34 [0.27, 0.42]). Points colored by site (purple = Hauru, slate = Maatea, green = Maharepa). All panels share the same x-axis. n = 112 *Pocillopora* corals across three reef sites.

![Figure 2. CAFI scaling with coral volume.](../output/figures/manuscript/fig2_scaling.png)

**Figure 3.** Scaling of individual species and taxonomic groups with coral volume. (A) Abundance versus colony volume for the 10 most prevalent CAFI species, fitted with negative binomial GLMs (log--log scale). (B) Abundance versus colony volume for four taxonomic groups (crabs, shrimps, fishes, snails), with direct-labeled NB GLM fits. (C) Species-level scaling exponents (β) with 95% bootstrap CI (1,000 site-stratified iterations). Blue: Redirection (β < 1); grey: Field of Dreams (CI spans 1.0); vermillion: super-linear (β > 1). (D) Taxonomic group scaling exponents with 95% bootstrap CI. Dashed vertical line: proportional scaling (β = 1). FDR correction (Benjamini--Hochberg) applied within species and group categories. n = 112 *Pocillopora* corals across three reef sites.

![Figure 3. Species and group scaling.](../output/figures/manuscript/fig3_species_group_scaling.png)

**Figure 4.** CAFI community composition differs among reef sites and along the coral size gradient. (A) Distance-based redundancy analysis (db-RDA) biplot of Hellinger-transformed community data (community ~ log(volume) + site). Each point represents one coral colony (n = 112); point size is proportional to coral volume. Ellipses show 80% confidence intervals for each site. Vectors indicate the top five species driving compositional variation (weighted-average scores on constrained axes), colored by taxonomic group. (B) Relative abundance of major taxonomic groups at each site. Percentages shown for groups comprising >10% of the site assemblage. Site abbreviations: Hauru (n = 38), Maatea (n = 39), Maharepa (n = 35).

![Figure 4. Community composition.](../output/figures/manuscript/fig4_composition.png)

**Figure 5.** CAFI diversity and abundance as predictors of coral physiological condition (BEF framework). (A) Species richness versus coral condition score (PC1). (B) Total CAFI abundance (square-root-scaled x-axis) versus condition. Points colored by site; solid lines show marginal fits from covariate-adjusted models (condition ~ predictor + log(volume) + site). Shading shows 95% CI. n = 84 colonies with physiological data.

![Figure 5. BEF diversity-condition.](../output/figures/manuscript/fig5_feedbacks.png)

---

## References

::: {#refs}
:::

---

## Supplementary Materials

See separate document: `combined_supplement.md` (14 supplementary figures, 11 supplementary tables, supplementary methods).
