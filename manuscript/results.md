# Results

## Q1: CAFI abundance and richness scale sublinearly with coral size

<!-- CRITIQUE: PROSE-01 fixed: added investigator voice to paragraph openings -->
We found that total CAFI abundance scaled sublinearly with colony volume, consistent with Propagule Redirection (Fig. 2A). The negative binomial GLM estimated a scaling exponent of β = 0.52 (95% bootstrap CI [0.44, 0.62]; Wald z vs. β = 1: z = −11.45, p < 0.001; bootstrap p < 0.001), decisively rejecting the Field of Dreams hypothesis (β = 1). Species richness also scaled sublinearly (Poisson GLM: z = 0.34 [0.27, 0.42]; z vs. 1 = −18.76, p < 0.001; Fig. 2C). However, rarefied richness (expected species at n = 20 individuals) showed no relationship with volume (OLS: slope = 0.14, SE = 0.30, t₆₄ = 0.47, p = 0.64, n = 68; colonies with <20 CAFI excluded). <!-- CRITIQUE: RES-01a fixed: removed interpretation creep from Q1 rarefied richness -->
We return to the interpretation of this pattern in the Discussion.

Per-capita CAFI density declined with colony size (log--log slope = −0.48; Fig. 2B). Because density equals abundance divided by volume, this slope is the algebraic complement of the abundance exponent (0.52 − 1 = −0.48), illustrating that sublinear scaling necessarily produces density dilution -- larger corals harbor fewer CAFI per unit volume.

We detected sublinear scaling at the species level: 11 of 21 prevalent species (52%) showed Redirection (bootstrap CI excluding 1.0), 10 showed scaling consistent with Field of Dreams (CI spanning 1.0), and none showed super-linear scaling (Fig. 3A,C; Table S1). The inverse-variance-weighted mean β across 21 species was 0.51 [0.45, 0.56], significantly below 1.0 (z = −18.8, p < 0.001; Fig. S4). This average treats species as independent estimates, although species co-occurring on the same corals are not strictly independent. Obligate symbionts -- *Trapezia* crabs (β = 0.43) and snapping shrimp *Alpheus lottini* (β = 0.37) -- consistently showed sublinear scaling, while facultative associates such as *Caracanthus maculatus* (β = 1.18, CI overlapping 1.0) and the coralliophiline snail *Galeropsis monodonta* (β = 1.27, CI overlapping 1.0) scaled proportionally.

Among the six taxonomic groups, five showed significant sublinear scaling: *Trapezia* crabs (β = 0.43 [0.35, 0.52], p < 0.001), shrimps (β = 0.50 [0.40, 0.60], p < 0.001), other crabs (β = 0.47 [0.24, 0.71], p < 0.001), fish (β = 0.74 [0.58, 0.92], p = 0.004), and other invertebrates (β = 0.50 [0.33, 0.67], p < 0.001). Gastropods were the sole exception (β = 0.94 [0.73, 1.17], p vs. 1 = 0.60), scaling proportionally (Fig. 3B,D; Fig. S7). All results were robust to taxonomic resolution: abundance β varied by <0.01 across five sensitivity scenarios, and richness z remained significantly sublinear in all cases (Fig. S6).

Species occurrence probability was size-dependent for 14 of 24 prevalent species (logistic GLM, FDR < 0.05; Fig. S12), with most showing increasing occurrence probability with coral size. Size-dependent scaling is thus common across the community, though 10 species showed no significant size dependence. (Spatial autocorrelation was negligible for all metrics: Moran's I = −0.004–0.024, all p > 0.28; Table S10.)

## Q2: Site pools and coral size structure community composition

<!-- CRITIQUE: RES-02a fixed: removed methods recapitulation; Fig. S1 moved to parenthetical -->
<!-- CRITIQUE: RES-09 fixed: biology-first opening (NMDS stress removed; leads with site structure finding) -->
Community composition differed significantly among the three reef sites (marginal PERMANOVA on Bray--Curtis dissimilarity of Hellinger-transformed abundances: volume F₁,₁₀₈ = 9.74, R² = 0.08, p = 0.001; site F₂,₁₀₈ = 3.74, R² = 0.06, p = 0.001; n = 112; <!-- CRITIQUE: RES-02b fixed: PERMDISP compressed to parenthetical -->
dispersion homogeneity confirmed; PERMDISP F = 0.89, p = 0.42; sampling adequacy confirmed, Fig. S1). <!-- CRITIQUE: RES-01b fixed: replaced interpretation with empirical statement -->
Together, site and volume explained ~14% of compositional variation, with 86% unexplained by the measured predictors. All pairwise site comparisons were significant after Bonferroni correction (all p_adj ≤ 0.006). The site effect was robust across all five distance metrics tested (Fig. S2). Balanced site subsampling (n = 35 per site, 500 iterations) confirmed the site effect was robust, with 100% of iterations significant at p < 0.05 (mean p = 0.001).

The three sites exhibited distinct taxonomic signatures (Fig. 4B). Maatea was characterized by high hermit crab abundance (33% of all CAFI, primarily *Calcinus latens*). Maharepa was dominated by obligate coral symbionts, with shrimps and crabs together comprising 71% of the assemblage. Hauru supported the most taxonomically balanced community and the highest proportion of coral-dwelling fishes (12%). Seven species drove the largest compositional differences among sites (envfit R² > 0.10, p < 0.01), including *Paragobiodon modestus* (R² = 0.46), *Fennera* spp. (R² = 0.39), *Calcinus latens* (R² = 0.39), and *Alpheus lottini* (R² = 0.27).

We also found that community composition shifted continuously along the coral size gradient. Distance-based redundancy analysis (db-RDA; Fig. 4A), after partialing out site effects, confirmed that <!-- CRITIQUE: RES-04 fixed: removed duplicate F-statistic (already reported in PERMANOVA) -->
coral volume explained 7.8% of compositional variation (p = 0.001). This size--composition gradient was robust to rarefaction (2.4%, F₁,₁₀₈ = 2.64, p = 0.001), confirming genuine compositional turnover rather than an abundance artifact.

Species with the strongest loadings on the db-RDA size axis included *Trapezia punctimanus* (score = −2.95) and *Luniella pugil* (−1.71), both associated with smaller corals, while *Euplica varians* and *Trapezia flavopunctata* (scores = 1.06) were most strongly associated with larger corals. Variance partitioning (adjusted R²) attributed 4.8% to volume alone, 3.9% shared between volume and site, with negligible unique site variation (−0.2%) and 91.5% residual. Compositional divergence among size classes was not significant after controlling for abundance through rarefaction (betadisper: p = 0.61; Fig. S3).

Community nestedness along the size gradient was not significant (NODF = 18.4, z = −1.09, p = 0.28). NODF ranges from 0 (no nestedness) to 100 (perfect nestedness); the observed value of 18.4 indicates low nestedness, confirming that small-coral faunas are not subsets of large-coral faunas. Combined with the significant db-RDA, this pattern is consistent with species turnover -- not passive accumulation -- along the size gradient: some species are gained while others are lost as coral size increases. <!-- CRITIQUE: CROSS-01 fixed: linked nestedness to Introduction's species-stacking prediction -->
This non-nested pattern is consistent with the species-stacking constraint: intraspecific density ceilings drive diversification along the size gradient rather than passive accumulation.

<!-- CRITIQUE: RES-02c fixed: spatial autocorrelation moved to end of Q1 as global diagnostic -->

<!-- CRITIQUE: PROSE-06 fixed: stronger directional header -->
## Q3: CAFI species richness positively predicts coral condition

<!-- CRITIQUE: RES-03 fixed: restructured Q3 to lead with independence-from-abundance -->
We found that corals harboring more CAFI species were in better physiological condition (n = 84 colonies with physiological data). <!-- CRITIQUE: CROSS-08 fixed: added Hochberg correction visibility -->
Species richness positively predicted coral condition (standardized β = 0.27, t₇₉ = 2.42, p = 0.018, Hochberg-corrected for two a priori tests; Fig. 5A; Fig. S11A), and total CAFI abundance showed a similar but weaker association (standardized β = 0.32, t₇₉ = 2.01, p = 0.048; Fig. 5B). Variance partitioning indicated that the richness signal was largely independent of abundance: of the incremental R² = 0.07 explained by CAFI predictors beyond volume and site, 29.1% was uniquely attributable to richness (~2% of total variance), less than 1% uniquely to abundance (<0.1% of total variance), and 70.8% was shared between them (~5% of total variance; Fig. S10A). The richness-unique fraction, while small in absolute terms, was roughly 30x larger than the abundance-unique fraction.

The overall model (richness + volume + site) explained modest variance (adjusted R² = 0.04). Because richness and abundance were strongly correlated (r = 0.84), partial regression including both predictors produced high collinearity (VIF = 6.2 for richness, 6.8 for abundance), complicating attribution of variance to diversity per se versus abundance.

In the path model, the standardized coefficient from richness to condition (β = 0.55) was larger in magnitude than from abundance to condition (β = 0.02), though neither path reached significance individually (p = 0.20 and p > 0.50, respectively; Fig. S10C). The path model did not fit significantly better than the null (Fisher's C p = 0.20), so these coefficients should be interpreted cautiously.
<!-- CRITIQUE: RES-01c fixed: removed interpretation creep from path model paragraph (Discussion handles this) -->

<!-- CRITIQUE: RES-06 fixed: added sample size and exclusion note for rarefied richness -->
Rarefied richness (expected species at n = 20 individuals) showed no relationship with condition (β = −0.07, p = 0.50, n = 47; colonies with fewer than 20 CAFI excluded from rarefaction; Fig. S8; Fig. S11B). We address the interpretive ambiguity of this test in the Discussion.

We found that no exploratory predictor -- Shannon diversity, *Trapezia* abundance, resident fish abundance, or *Galeropsis monodonta* abundance -- survived BH-FDR correction (all p_FDR > 0.80; Fig. S11C), though effect directions mostly matched expectations (Trapezia positive, fish positive). The exception was *Galeropsis monodonta*, which showed an unexpected positive association with condition despite being a known corallivore (Fig. S11D--E). Community composition (PC1_CAFI) did not predict condition (p > 0.10; tested separately in the supplement).

The reverse direction (condition → CAFI) showed one nominally significant result (condition → total CAFI abundance, raw p = 0.037), but this did not survive BH-FDR correction across the seven reverse-direction tests (Fig. S11F). Key species from the companion experimental study (8 of 10 with sufficient survey occurrences) were tested individually; no species was individually significant after FDR correction (Table S4). The survey had 80% power to detect a small-to-medium effect size (R² ≈ 0.09), indicating that null results for individual species are credible rather than underpowered.

Species-level analysis of individual condition traits (protein, carbohydrate, zooxanthellae, AFDW) revealed no significant associations after FDR correction across 95 species × trait tests (Fig. S13), suggesting that no single species drives the richness--condition relationship. The nine strongest individual associations (Fig. S14) showed effect directions consistent with mutualism hypotheses (positive for guard crabs and shrimps). In the neighborhood/landscape models, coral volume, neighborhood density, and site did not predict condition (all p > 0.30; Fig. S5).

## Supporting results (Supplement)

<!-- CRITIQUE: RES-05 fixed: renamed supplement subheadings to informative result statements -->
### Neighborhood density does not predict CAFI abundance

Neighborhood analyses were restricted to the 61 colonies with 5-m survey data (Fig. S5; Table S5). Neighborhood density (n_neighbors) did not significantly predict CAFI abundance (β = −0.005, z = −0.89, p = 0.37) or species richness (β = −0.003, z = −0.88, p = 0.38) after controlling for coral volume and site. Mean neighbor distance showed a significant negative association with richness (β = −0.005, z = −3.20, p = 0.001) and Shannon diversity (β = −0.007, t = −3.49, p = 0.001), suggesting that corals closer to neighbors support slightly higher diversity, though this pattern was not observed for total abundance. No size × neighborhood interactions were significant (all p > 0.10). The power analysis indicated 65% power for medium effects (R² = 0.10), so null results for density may reflect insufficient power rather than a true absence of neighborhood effects.

### Pairwise co-occurrences are explained by volume and site

We detected no significant pairwise co-occurrences using volume-weighted null models after FDR correction (0 of 528 pairs at FDR < 0.05; Fig. S9; Table S6). The two strongest signals were *Harpiliopsis beaupresii*--*Paragobiodon modestus* (SES = −3.43, p_FDR = 0.11) and *H. beaupresii*--*H. spinigera* (SES = −3.42, p_FDR = 0.21), both negative but non-significant after correction, indicating that pairwise associations are largely explained by volume and site. Intraspecific density analysis identified six species with significant mating-pair aggregation (FDR < 0.05): *Synalpheus charon* (SES = 10.3), *Alpheus lottini* (SES = 10.3), *Harpiliopsis beaupresii* (SES = 5.4), *Trapezia serenei* (SES = 4.4), *Alpheus pachychirus* (SES = 3.9), and *Trapezia bidentata* (SES = 3.8; Table S7).

### Survey and experimental scaling exponents are broadly concordant

Seven species had scaling data in both the survey and companion experiment. Of these, two showed sublinear scaling in the survey but proportional scaling during experimental colonization (*Calcinus latens*, *Alpheus lottini*), one was sublinear in both studies (*Trapezia serenei*), and four scaled proportionally in both (*Caracanthus maculatus*, *Harpiliopsis spinigera*, *Galeropsis monodonta*, *Periclimenes watamuae*; Table S8).
