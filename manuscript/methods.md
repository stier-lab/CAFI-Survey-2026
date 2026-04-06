## Methods

### Study system and survey design

We surveyed 114 *Pocillopora* coral colonies across three reef habitats on Mo'orea, French Polynesia (17°30'S, 149°50'W) during June--August 2019 (Fig. 1A). Colonies were identified to species using mitochondrial open reading frame (mtORF) haplotyping via PCR and restriction fragment analysis (101 of 114 successfully genotyped; 11 PCR failures, 2 missing samples). We assigned haplotypes to species following @JohnstonCunningBurgess2022: *P. grandis* (n = 49), *P. meandrina* (n = 34), *P. tuahiniensis* (n = 7), and *P. verrucosa* (n = 10). Field morphotype concordance with genetic species varied among morphotypes (94% for "eydouxi"/*P. grandis*, 66% for "meandrina"/*P. meandrina*, 33% for "verrucosa"/*P. verrucosa*; Table S15A), consistent with the well-documented cryptic species problem in *Pocillopora* [@Burgess2021]. We selected sites to span the range of reef environments: Hauru (HAU; fringing reef, north shore; n = 38), Maatea (MAT; lagoon/back reef, east shore; n = 39), and Maharepa (MRB; barrier reef, north shore; n = 37). At each site, divers targeted small, medium, and large colonies to ensure coverage of the three-order-of-magnitude size gradient. Two MRB colonies lacked volume measurements and were excluded from size-dependent analyses (n = 112; MRB n = 35).

### Coral measurements

We estimated colony size as ellipsoidal volume: V = (4/3) × π × (L/2) × (W/2) × (H/2), where L, W, and H are maximum length, width, and height (cm). The ellipsoidal approximation is standard for *Pocillopora*, which have roughly symmetrical branching morphology; it includes interstitial space, but this proportional bias does not affect log-linear scaling estimates. Colony volumes ranged from 21 to 42,333 cm³ (Fig. 1B).

For 61 of the 114 colonies ("neighborhood corals"), we counted all *Pocillopora* neighbors within a 5-m radius and measured inter-colony distances. We surveyed the remaining 53 colonies for CAFI and colony volume only. Neighbor count ranged from 1 to 76 (median = 17; Fig. 1C).

### CAFI collection and identification

We extracted coral-associated fishes and invertebrates (CAFI) by enclosing each colony in a fine-mesh bag, detaching it from the substrate, and preserving all contents in 95% ethanol; macroboring and deeply cryptic fauna embedded within the skeleton were excluded (this may underrepresent boring taxa, which increase with colony size). We identified specimens to the lowest practical taxonomic level (154 to species, 89 to genus or higher), yielding 243 operational taxonomic units (OTUs; two unresolved fish morphotypes were pooled as a single OTU) from 3,989 individuals. We computed species richness and Shannon diversity (H') per colony. We computed rarefied richness (expected species at n = 20 individuals; vegan::rarefy; Oksanen et al. 2022) to control for the abundance--richness correlation.

### Coral physiological condition

We collected tissue samples from the upper branch of each colony and assayed total protein (mg cm⁻²), carbohydrate (mg cm⁻²), zooxanthellae density (cells cm⁻²), and ash-free dry weight (AFDW; mg cm⁻²). Because trait values vary with the position of the tissue sample along the branch (see Supplementary Methods), we regressed each trait on stump length and nubbin length and used the residuals as position-corrected values. We summarized condition as PCA axis 1 of the four corrected traits (61.6% of variance; all loadings positive). Physiology data were available for 108 colonies; after removing incomplete cases, 84 were retained for BEF analyses (missingness was predicted by site but not by volume, richness, or abundance; Table S14). Position measurements were unavailable for 10 Maharepa colonies, reducing site balance in the BEF subsample.

### Statistical analyses

All analyses were conducted in R v4.4 (package versions archived with code). Natural logarithms were used throughout; site was included as a fixed effect because three sites provide too few levels for reliable random-intercept estimation [@Bolker2009]. PERMANOVA results depend on permutation randomness; we set seeds before all permutation-based tests for reproducibility. Sample sizes vary across analyses: n = 114 surveyed, n = 112 with volume data (2 MRB colonies excluded), n = 84 with complete physiology for BEF analyses (Table S14), n = 61 with neighborhood surveys, and n = 68 (Q1) or n = 47 (Q3, restricted to the physiology subset) with ≥20 CAFI for rarefied richness. Additional details on multiple-testing correction, diagnostics, bootstrap implementation, and secondary analyses are provided in the Supplementary Methods.

#### Q1: Scaling of CAFI with coral size

We estimated the scaling exponent β from a negative binomial GLM (log link; chosen over Poisson to accommodate overdispersion): `total_cafi ~ log(volume) + site`, where the coefficient on log(volume) directly estimates β. We contrasted proportional (β = 1), sublinear (β < 1), and super-linear (β > 1) scaling. Species richness was modeled with a Poisson GLM, yielding the species--area exponent *z*. Rarefied richness (expected species at n = 20 individuals, chosen to retain the majority of colonies [68 of 112] while adequately equalizing abundance) was regressed on log(volume) using OLS. Species- and group-level analyses are in the supplement.

#### Q2: Community composition

We analyzed community composition for the 112 colonies with volume data. Species abundances were Hellinger-transformed [@LegendreGallagher2001]. We tested for compositional differences using marginal (Type III) PERMANOVA on Bray--Curtis dissimilarities with log(volume) and site and assessed dispersion homogeneity with PERMDISP. Two-dimensional NMDS provided visualization.

To test whether composition shifts continuously along the size gradient, we used distance-based redundancy analysis (db-RDA) with log(volume) as the constrained variable and site partialed out. We also tested nestedness using NODF against quasiswap null matrices [@AlmeidaNeto2008; @MiklosPodani2004] and decomposed total beta diversity into turnover and nestedness-resultant components using `betapart` [@Baselga2010].

Additional analyses of site robustness, iterated rarefaction, host architecture, and community assembly are described in the Supplementary Methods.

#### Q3: CAFI species richness and coral condition (BEF framework)

We tested whether CAFI diversity predicts coral condition using a biodiversity--ecosystem function (BEF) framework (Tilman et al. 2014). Forward models testing CAFI effects on condition took the form:

    condition_PC1 ~ CAFI_predictor + log(volume) + site

with OLS error structure (Breusch--Pagan confirmed homoscedasticity; HC3-robust SEs reported as sensitivity). Count-based predictors were square-root transformed; richness was untransformed. We corrected the two a priori predictors (richness, abundance) with Hochberg step-up (k = 2), four exploratory predictors with BH-FDR (k = 4), and PC1 independently. To disentangle diversity from abundance, we used variance partitioning and a piecewise structural equation model [@Lefcheck2016], and compared raw versus rarefied richness (n = 20; reducing n from 84 to 47). Full correction details, diagnostics, and sensitivity analyses are in the Supplementary Methods.

#### Neighborhood context

Using the 61 neighborhood-surveyed colonies, we modeled abundance with negative binomial GLMs, richness with Poisson GLMs, and Shannon diversity with OLS as functions of log(volume), neighbor count (*Pocillopora* within 5 m), total neighbor volume, mean inter-colony distance, and site. The distance--richness relationship was retested using rarefied richness (n = 20; available for 39 colonies with ≥20 individuals). Additional neighborhood, co-occurrence, host-identity, community-assembly, and robustness analyses are described in the Supplementary Methods.

---
