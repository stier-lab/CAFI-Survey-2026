# Results

## Q1: CAFI abundance and richness scale sublinearly with coral size

Total CAFI abundance scaled sublinearly with colony volume (n = 112; β = 0.52, 95% bootstrap CI [0.44, 0.62]; p < 0.001 vs β = 1; Fig. 2A), consistent with Propagule Redirection: corals doubling in volume gained only ~43% more CAFI. Per-capita density declined accordingly (slope = −0.48; Fig. 2B). Species richness also scaled sublinearly (z = 0.34 [0.27, 0.42], p < 0.001; Fig. 2C), but rarefied richness showed no relationship with volume (p = 0.64, n = 68), indicating that the raw richness signal is tightly coupled to the abundance gradient.

Sublinear scaling extended to individual species: 11 of 21 prevalent species showed Redirection, 10 were consistent with Field of Dreams, and none showed super-linear scaling (Fig. 3A,C; Fig. S4; Table S1). Obligate symbionts---*Trapezia* (β = 0.43) and *Alpheus lottini* (β = 0.37)---consistently scaled sublinearly, while facultative associates such as *Galeropsis monodonta* (β = 1.27) scaled proportionally.

Five of six taxonomic groups scaled sublinearly (β range 0.43–0.74; Fig. 3B,D); gastropods were the sole exception (β = 0.94, CI spanning 1.0). Size-dependent occurrence reinforced these patterns: 14 of 24 prevalent species showed significant occurrence–volume relationships (logistic GLM, FDR < 0.05). All results were robust to taxonomic resolution (Fig. S6).

### Neighborhood context

Among the 61 colonies with neighborhood surveys (5-m radius), closer neighbors predicted higher richness (β = −0.005, z = −3.20, p = 0.001) and Shannon diversity (p = 0.001), but not total abundance (p = 0.78; Fig. S5; Table S4). Neither neighbor count nor total neighbor volume predicted any response (all p > 0.37). The distance–richness relationship survived rarefaction (n = 39, p = 0.005), confirming genuine species accumulation rather than a passive sampling artifact.

Beyond shaping how many species a colony hosts, colony size also structures which species are present.

## Q2: Site pools and coral size structure community composition

Marginal PERMANOVA confirmed that both site and colony size structured CAFI composition (volume R² = 0.08, site R² = 0.06, both p = 0.001; n = 112; PERMDISP p = 0.42), together explaining ~14% of variation. The site effect was robust across five distance metrics and balanced subsampling (Fig. S2). The three sites exhibited distinct taxonomic signatures (Fig. 4B): Maatea was characterized by hermit crabs (33% of CAFI), Maharepa by obligate symbionts (71%), and Hauru by fishes (12%). Species vector fitting (envfit) identified *Trapezia punctimanus* (R² = 0.53), *Harpiliopsis beaupresii* (R² = 0.30), and *Paragobiodon modestus* (R² = 0.25) as the strongest drivers.

Community composition also shifted continuously along the coral size gradient. Distance-based redundancy analysis (db-RDA; n = 112), after partialing out site effects, confirmed that coral volume explained 7.8% of compositional variation (p = 0.001; Fig. 4A). This size–composition gradient was robust to rarefaction (mean R² = 2.4% across 100 rarefaction draws, F₁,₁₀₈ = 2.64, p = 0.001), confirming genuine compositional turnover rather than an abundance artifact.

*Trapezia punctimanus* (db-RDA score = −2.95) and *Luniella pugil* (−1.71) loaded toward smaller corals, while *Euplica varians* and *Trapezia flavopunctata* (1.06) loaded toward larger corals. Variance partitioning attributed 4.8% to volume alone, 3.9% shared with site, and 91.5% residual. Compositional divergence among size classes was not significant after rarefaction (PERMDISP: p = 0.61; Fig. S3), confirming that apparent convergence of large-coral assemblages reflects passive sampling rather than deterministic assembly.

Nestedness was not significant (NODF = 18.4, p = 0.28); beta diversity partitioning confirmed that species replacement accounted for 81% of total dissimilarity (Table S11). Combined with the significant db-RDA, this indicates species turnover---not passive accumulation---along the size gradient.

Community variability declined from small to large corals (distance-to-centroid ANOVA F = 9.21, p < 0.001; Fig. S15), consistent with compositional convergence as communities mature. Because *Pocillopora* comprises multiple cryptic species with distinct branching architectures, we tested whether host species identity contributes to compositional variation beyond colony size (Supplementary Methods). Coral genetic species predicted CAFI richness (p = 0.007) and composition (PERMANOVA R² = 0.08, p = 0.001) but not total abundance (p = 0.24; Table S16). Three-way variance partitioning attributed comparable unique variation to host architecture (5.6%) and colony volume (4.6%), both exceeding site (2.9%; Table S16). The corallivore *Galeropsis monodonta* showed strong host specificity to *P. verrucosa* (7-fold enrichment), while obligate *Trapezia* did not differ among coral species (p = 0.51), suggesting that guard crabs colonize based on colony geometry rather than host identity (Tables S17--S18; Figs. S12--S13).

### Pairwise co-occurrences are explained by volume and site

After accounting for coral volume and site, no species pair co-occurred more or less often than expected (0 of 528 pairs significant after FDR correction; Fig. S9; Table S5), indicating that pairwise associations are largely explained by habitat size and location. However, six species showed significant mating-pair aggregation -- concentrating conspecific pairs on a subset of corals (Table S6) -- consistent with the territorial pair-bonding that underpins species stacking.

Because composition turns over along the size gradient rather than merely accumulating species, variation in CAFI richness reflects partially non-overlapping functional assemblages -- the prerequisite for testing whether diversity predicts coral condition.

## Q3: CAFI species richness positively predicts coral condition

Species richness positively predicted coral condition (n = 84; standardized β = 0.36, 95% CI [0.06, 0.66], p = 0.018, Hochberg-adjusted p = 0.036; Fig. 5A), and total abundance showed a weaker association (β = 0.32, p = 0.048; Fig. 5B; Table S2). Richness and abundance were pre-specified a priori hypotheses; we therefore report these tests at nominal significance levels (global BH-FDR across all six predictors shown for transparency in Table S2). Variance partitioning attributed 29.1% of the incremental R² to richness uniquely, less than 1% to abundance uniquely, and 70.8% shared (Table S10; Fig. S10A).

The overall model explained limited variance (adjusted R² = 0.04). Mediation analysis confirmed that the richness effect operates directly rather than through abundance (ACME p = 0.69; Table S12). The path model suggested a stronger richness-to-condition path (β = 0.55) than abundance-to-condition (β = 0.02), though neither reached significance individually (Table S10; Fig. S10C). Adding coral morphotype as a covariate absorbed the richness signal (p = 0.27 vs 0.015; Table S13), suggesting morphotype-associated variation may partly drive the association (see Discussion). In contrast, adding genetic species identity (mtORF haplotype) did not absorb the signal (p = 0.010; Table S15), indicating that morphotype captures architectural variation beyond species identity alone.

Rarefied richness showed no relationship with condition (p = 0.50, n = 47; minimum detectable R² = 0.15 at 80% power; Fig. S8). No exploratory predictor survived FDR correction (Table S2), and reverse-direction tests were non-significant (Fig. S11C; Table S2).
