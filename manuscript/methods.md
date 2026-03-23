# Methods

## Study system and survey design

We surveyed 114 *Pocillopora* coral colonies across three reef habitats on Mo'orea, French Polynesia (17°30'S, 149°50'W) during June--August 2019 (Fig. 1A). Colonies were identified to species using mitochondrial open reading frame (mtORF) haplotyping via PCR and restriction fragment analysis (101 of 114 successfully genotyped; 11 PCR failures, 2 missing samples). Haplotypes were assigned to species following @JohnstonCunningBurgess2022: *P. grandis* (n = 49), *P. meandrina* (n = 34), *P. tuahiniensis* (n = 7), and *P. verrucosa* (n = 10). Field morphotype identifications showed limited concordance with genetic species assignments (Table S15), consistent with the well-documented cryptic species problem in *Pocillopora* [@Burgess2021]. We selected sites to span the range of reef environments: Hauru (HAU; fringing reef, north shore; n = 38), Maatea (MAT; lagoon/back reef, east shore; n = 39), and Maharepa (MRB; barrier reef, north shore; n = 37). At each site, divers targeted small, medium, and large colonies to ensure coverage of the three-order-of-magnitude size gradient. Two MRB colonies lacked volume measurements and were excluded from size-dependent analyses (n = 112; MRB n = 35).

## Coral measurements

We estimated colony size as ellipsoidal volume: V = (4/3) × π × (L/2) × (W/2) × (H/2), where L, W, and H are maximum length, width, and height (cm). The ellipsoidal approximation is standard for *Pocillopora*, which have roughly symmetrical branching morphology; it includes interstitial space, but this proportional bias does not affect log-linear scaling estimates. Colony volumes ranged from 21 to 42,333 cm³ (Fig. 1B).

For 61 of the 114 colonies ("neighborhood corals"), we counted all *Pocillopora* neighbors within a 5-m radius and measured inter-colony distances. We surveyed the remaining 53 colonies for CAFI and colony volume only. Neighbor count ranged from 1 to 76 (median = 17; Fig. 1C).

## CAFI collection and identification

We extracted CAFI by enclosing each colony in a fine-mesh bag, detaching it from the substrate, and preserving all contents in 95% ethanol; macroboring and deeply cryptic fauna embedded within the skeleton were excluded. Specimens were identified to the lowest practical taxonomic level (154 to species, 89 to genus or higher), yielding 243 operational taxonomic units (OTUs; two unresolved fish morphotypes were pooled as a single OTU) from 3,989 individuals. Species richness and Shannon diversity (H') were computed per colony. Rarefied richness (expected species at n = 20 individuals; vegan::rarefy; Oksanen et al. 2022) controlled for the correlation between abundance and raw richness.

## Coral physiological condition

We collected tissue samples from the upper branch of each colony and assayed total protein (mg cm⁻²), carbohydrate (mg cm⁻²), zooxanthellae density (cells cm⁻²), and ash-free dry weight (AFDW; mg cm⁻²). Because trait values vary with the position of the tissue sample along the branch (see Supplementary Methods), we regressed each trait on stump length and nubbin length and used the residuals as position-corrected values. We summarized condition as PCA axis 1 of the four corrected traits (61.6% of variance; all loadings positive). Physiology data were available for 108 colonies; after removing incomplete cases, 84 were retained for BEF analyses. Position measurements were unavailable for 10 Maharepa colonies, reducing site balance in the BEF subsample.

## Statistical analyses

All analyses were conducted in R v4.4 (session information archived at [DOI to be assigned upon acceptance]). Natural logarithms were used throughout; site was included as a fixed effect (three sites provide too few levels for random intercepts; Bolker et al. 2009).
We applied three levels of multiple-testing correction: (1) Hochberg step-up for two pre-specified BEF predictors (richness, total abundance; k = 2); (2) Benjamini--Hochberg FDR for four exploratory predictors (Shannon diversity, *Trapezia*, fish, and *Galeropsis* Blainville, 1832 abundance; k = 4); and (3) FDR across all species in species-level analyses. GLM diagnostics included simulated residuals (DHARMa), Cook's distance (4/n), and variance inflation factors. All tests used α = 0.05. Code is available at [URL to be provided upon acceptance].

### Q1: Scaling of CAFI with coral size

We tested how CAFI abundance scales with coral volume by fitting the power-law relationship N = aV^β^ using a negative binomial GLM (log link; negative binomial to accommodate overdispersion in count data):

    total_cafi ~ log(volume) + site

where the coefficient on log(volume) directly estimates β. We tested three hypotheses: Field of Dreams (β = 1; proportional scaling), Propagule Redirection (β < 1; sublinear, density dilution; Stier & Osenberg 2010), and super-linear scaling (β > 1). The null (β = 1) was tested using a Wald test (*z*~Wald~ = [β − 1]/SE) and a two-sided bootstrap proportion test (1,000 site-stratified replicates). Bootstrap BCa confidence intervals served as the primary inference; percentile CIs were substituted where BCa acceleration constants could not be estimated. Species richness was modeled analogously using a Poisson GLM, yielding the species--area relationship exponent *z* (S = c × V^*z*^; Preston 1962). Throughout, *z* refers to this SAR exponent; Wald z-statistics are denoted *z*~Wald~. Rarefied richness (expected species at n = 20 individuals; available for 68 of 112 colonies with ≥20 CAFI) was regressed on log(volume) using OLS.
Per-capita CAFI density (total CAFI / colony volume) equals β − 1 on the log--log scale by definition and was not independently modeled.

The scaling analysis was repeated for 21 prevalent species (≥30 individuals and ≥15% prevalence) and six taxonomic groups (*Trapezia* crabs, shrimps, other crabs, gastropods, fish, other invertebrates). We computed an inverse-variance-weighted mean β as a fixed-effect summary (weights = 1/SE²); because species co-occur on the same corals, the SE may be underestimated. We applied FDR correction within species-level and group-level test families and characterized size-dependent occurrence using logistic GLMs (presence ~ log(volume) + site) for 24 species with ≥15% prevalence.

### Q2: Community composition

We analyzed community composition for the 112 colonies with volume data. Species abundances were Hellinger-transformed (Legendre & Gallagher 2001). We tested for compositional differences using marginal (Type III) PERMANOVA (vegan::adonis2) on Bray--Curtis dissimilarities with 999 permutations, including log(volume) and site. We assessed dispersion homogeneity with permutational analysis of multivariate dispersions (PERMDISP; vegan::betadisper). We used two-dimensional NMDS ordination for visualization and identified species associated with site differences using envfit (R² > 0.10, p < 0.01). We compared sites pairwise using Bonferroni-corrected PERMANOVA (k = 3) and validated robustness across five distance metrics and 500 balanced site-subsampling iterations (Fig. S2; Supplementary Methods).

To test whether composition shifts continuously along the size gradient, we used distance-based redundancy analysis (db-RDA, a constrained ordination partitioning community variation by environmental predictors; vegan::dbrda) with log(volume) as the constrained variable and site partialed out. Variance partitioning (vegan::varpart) decomposed community variation into volume-only, site-only, shared, and residual fractions. We computed species scores on the constrained axis as weighted averages of site scores and verified robustness by repeating the db-RDA on iterated-rarefied distances (100 iterations at minimum observed abundance; Supplementary Methods). We tested nestedness using the nestedness metric based on overlap and decreasing fill (NODF; Almeida-Neto et al. 2008) against 999 quasiswap null matrices (Miklós & Podani 2004). We decomposed total beta diversity into turnover (Simpson dissimilarity) and nestedness-resultant components using the betapart package (Baselga 2010).

### Q3: CAFI species richness and coral condition (BEF framework)

We tested whether CAFI diversity predicts coral condition using a biodiversity--ecosystem function (BEF) framework (Tilman et al. 2014). Forward models testing CAFI effects on condition took the form:

    condition_PC1 ~ CAFI_predictor + log(volume) + site

with OLS standard errors as primary inference. Breusch--Pagan tests confirmed homoscedasticity (BP p > 0.5); HC3 standard errors are reported in the supplement as a sensitivity check (Long & Ervin 2000). Count-based predictors (total abundance, *Trapezia*, fish, *Galeropsis* counts) were square-root transformed; species richness was untransformed.

To disentangle diversity from abundance effects, we used partial regression (richness + √abundance in the same model), hierarchical R² variance partitioning, and a piecewise structural equation model (piecewiseSEM; Lefcheck 2016) with volume → richness → condition and volume → abundance → condition pathways (z-scored predictors; Fisher's C fit). To test whether the richness--condition relationship was an abundance artifact, we compared models using raw versus rarefied richness (n = 20; reducing available sample from 84 to 47). The interpretive ambiguity of this comparison is addressed in the Discussion.

We tested whether the richness effect operates through abundance using bootstrap mediation analysis (mediation package; 1,000 iterations; treatment = richness, mediator = √abundance, outcome = condition). As sensitivity checks, we added coral morphotype and genetic species identity (mtORF haplotype) as covariates in separate models to assess whether host identity confounds the richness--condition association (Tables S13, S15). Species-specific contributions were tested for 19 prevalent species (≥5 corals; trait ~ √(species abundance) + log(volume) + site) across five condition measures, with FDR correction across all 95 tests. Reverse models (condition → CAFI) tested seven response variables with FDR correction. Cross-study concordance with the companion experiment was evaluated with a binomial test, and power analysis confirmed sufficient sample size to detect the experiment's effect sizes (Supplementary Methods).

### Neighborhood context

Using the 61 neighborhood-surveyed colonies, we fitted full models with four predictors: log(volume), neighbor count (*Pocillopora* within 5 m), total neighbor volume, and mean inter-colony distance, plus site. Abundance was modeled with negative binomial GLMs, richness with Poisson GLMs, and Shannon diversity with OLS; AIC-based backward elimination identified the best predictor subset. The distance--richness relationship was retested using rarefied richness (n = 20; available for 39 colonies with ≥20 individuals). Size × neighborhood interactions and functional group responses were tested with FDR correction. Compositional variability across neighbor density was assessed with PERMDISP. Power analysis indicated 65% power for medium effects (R² = 0.10) at n = 61 (Supplementary Methods).

### Co-occurrence patterns (Supplement)

We tested pairwise co-occurrence using a volume-weighted Bernoulli null model (10,000 iterations; species with ≥10 occurrences; Stier et al. 2012) with FDR correction across all pairs. We tested intraspecific density patterns with a multinomial allocation null model and assessed size-dependent co-occurrence by repeating the null model within volume terciles. Full null-model algorithms are described in Supplementary Methods.

### Sensitivity analyses

We tested all results for robustness across five taxonomic resolution scenarios (baseline 243 taxa, species-only, merge-up, lump-down, rare-excluded; Fig. S6; Table S8) and assessed spatial autocorrelation using Moran's I (Table S9).

### Ethical statement

Research was conducted under permits issued by the Haut-commissariat de la République en Polynésie française, in collaboration with the Université de la Polynésie française (CRIOBE). All protocols followed institutional and local regulations for scientific collection on the reefs of Mo'orea. Invertebrate collections did not require IACUC review; coral tissue sampling protocols were approved under broader program oversight at the authors' institution.

### Data accessibility

Data and analysis code are available at [DOI to be assigned upon acceptance]. Data and code will be archived in the Dryad Digital Repository upon acceptance.
