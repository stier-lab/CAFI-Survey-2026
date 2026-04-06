---
title: "Sublinear scaling drives density dilution and compositional turnover in coral-associated fauna"
bibliography: references.bib
csl: journal-of-animal-ecology.csl
link-citations: true
lang: en
---

<!-- LEGACY ASSEMBLED EXPORT. Canonical source: manuscript/comprehensive_manuscript.md; see manuscript/README.md. -->

# Sublinear scaling drives density dilution and compositional turnover in coral-associated fauna

Adrian C. Stier^1,2^*, Alexander Primo^3^, Joseph S. Curtis^1,4^, Craig W. Osenberg^3^

^1^Department of Ecology, Evolution, and Marine Biology, University of California Santa Barbara, Santa Barbara, CA 93106, USA
^2^Marine Science Institute, University of California Santa Barbara, Santa Barbara, CA 93106, USA
^3^Odum School of Ecology, University of Georgia, 140 E Green St, Athens, GA 30602, USA
^4^Department of Marine Science, University of Otago, PO Box 56, Dunedin, Otago 9054, New Zealand

*Corresponding author: Adrian C. Stier (astier@ucsb.edu)

---

## Abstract

1. Animals that occupy biogenic habitats -- corals, trees, kelps, oysters -- can alter the condition of their living hosts through defence, nutrient subsidies, or tissue consumption. Whether occupant abundance scales proportionally or sublinearly with habitat size determines per-capita density and thus the strength of these interactions. Yet whether sublinear scaling also reshapes occupant community composition, and whether that compositional change covaries with host condition, has not been tested in a single system.

2. We surveyed 114 *Pocillopora* coral colonies spanning three orders of magnitude in volume across three reef habitats in Mo'orea, French Polynesia, cataloguing 3,989 coral-associated fishes and invertebrates (CAFI) from 243 taxa to test a three-link chain: scaling, compositional turnover, and diversity--condition covariation.

3. Abundance scaled sublinearly with colony volume (β = 0.52 [95% CI: 0.44, 0.62]): corals doubling in volume gained only 43% more occupants. Territorial pair-bonding capped intraspecific density, so larger colonies accumulated new species rather than more conspecifics. Composition shifted through species replacement rather than nestedness, with colony volume explaining 7.8% of compositional variation after controlling for site. Proximity to neighbouring colonies independently predicted higher richness, an effect that survived rarefaction (p = 0.005). Species richness positively covaried with coral physiological condition (p = 0.018, adjusted R² = 0.036), though rarefied richness did not, leaving the causal pathway unresolved.

4. The first two links -- sublinear scaling and compositional turnover -- are robust across sensitivity analyses and generalise to any biogenic habitat where occupants draw from finite colonist pools. The diversity--condition link requires experimental confirmation. This chain applies wherever habitat loss or climate stress restructures the animal communities on which foundation species depend.

**Keywords:** biodiversity--ecosystem function, community assembly, coral reef, density dilution, foundation species, habitat structure, species--area relationship, sublinear scaling

---

## Introduction

Animals that inhabit biogenic structures -- crabs sheltering among coral branches, beetles tunnelling beneath bark, mites nestled in moss cushions -- depend on habitat-forming organisms for shelter, food, and mating sites (Dayton 1972; Ellison et al. 2005). A central question for these systems is how occupant abundance scales with habitat size: proportional scaling maintains constant per-unit density, whereas sublinear scaling means larger habitats support fewer occupants per unit volume (Connor & McCoy 1979; MacArthur & Wilson 1967). Where the habitat is itself a living organism, the scaling regime determines per-unit occupant density and thus the cumulative magnitude of occupant--host interactions (Bronstein 1994). Occupants alter host condition -- here defined as the host's integrated physiological state, reflected in tissue reserves and symbiont density -- sometimes positively through predator defence and nutrient subsidies (Glynn 1983; Holbrook et al. 2008), sometimes negatively through tissue consumption or growth suppression (Silliman et al. 2005; Shima et al. 2010). In ecosystems structured by foundation species, changes in habitat-former size can therefore alter both the richness and composition of associated animal communities and the physiological condition of the habitat former (Angelini et al. 2011; Ellison et al. 2005), yet empirical tests tracing the full chain -- from habitat scaling through community assembly to host condition -- remain rare.

Several non-exclusive frameworks predict how occupants scale with habitat size. Under proportional scaling (the "Field of Dreams" scenario), colonists arrive in proportion to available habitat, maintaining constant per-unit density. Under Propagule Redirection, habitat patches compete for a shared colonist pool and larger patches dilute per-unit density (Stier & Osenberg 2010) -- a concept rooted in the attraction-versus-production framework for reef fish habitat (Carr & Hixon 1997). The More Individuals Hypothesis offers a complementary prediction for species richness: larger habitats support more individuals, and richness increases as a passive sampling artifact of greater abundance (Srivastava & Lawton 1998). Habitat heterogeneity provides a fourth mechanism, in which larger habitats contain more microhabitat types and thereby accommodate more species independently of abundance (Lawton 1983). Which regime predominates depends on colonist-supply dynamics, density-dependent mortality, and community maturity (Hamman et al. 2018). The prevailing regime determines the consequences of habitat change: under Propagule Redirection, small habitats accumulate disproportionately high per-unit occupant densities, intensifying occupant--host interactions whose net sign depends on the balance of mutualists and antagonists (McKeon et al. 2012; Stier & Osenberg 2024).

Density alone, however, cannot predict host condition if the occupant assemblage changes in composition -- not just abundance -- along the habitat-size gradient. In systems where occupants defend territories and exclude conspecifics, intraspecific density ceilings cause larger habitats to accumulate new species rather than more conspecifics (Castro 1978; Huber & Coles 1986). Because occupants range from net mutualists (branch defenders, nutrient subsidisers) to net antagonists (corallivores, tissue consumers), the composition of the assemblage should influence host condition alongside its abundance. Experimental work in coral systems supports this expectation: guard crabs defend corals against predators (Glynn 1983; McKeon et al. 2012), resident fishes enhance coral growth through nutrient subsidies (Holbrook et al. 2008), and vermetid gastropods suppress coral growth (Shima et al. 2010). Biodiversity--ecosystem function (BEF) theory predicts that more diverse assemblages produce greater net ecosystem function through niche complementarity (Tilman et al. 2014; Loreau & Hector 2001), suggesting that corals hosting richer CAFI assemblages should receive stronger net mutualist benefits. However, under sublinear scaling, richness and abundance are tightly correlated, making it difficult to separate diversity effects from abundance effects in observational data. To our knowledge, no prior study has simultaneously evaluated habitat quantity, community compositional turnover, and host condition within a single biogenic-habitat system.

Branching corals in the genus *Pocillopora* provide a tractable field system for testing these linked predictions. *Pocillopora* colonies create densely branching three-dimensional structures that support assemblages of coral-associated fishes and invertebrates (CAFI), with the most densely occupied colonies exceeding 65 individuals from over 20 species (Curtis et al. 2023; Counsell et al. 2018). Many CAFI are obligate symbionts that occupy colonies as territorial mating pairs, defending branch junctions against conspecifics while tolerating heterospecifics (Castro 1978; Huber & Coles 1986). This behavioural rule caps intraspecific density while permitting interspecific coexistence, producing "species stacking": as colony size increases, increases in total occupancy come predominantly from adding new species rather than more conspecifics (Stier et al. 2012; Counsell & Donahue 2021). Species stacking predicts that small-coral assemblages should not be nested subsets of large-coral assemblages, and that colony volume should predict occupant composition even after controlling for abundance. Mo'orean *Pocillopora* comprises at least four cryptic genetic species (Johnston et al. 2022) that differ in bleaching susceptibility and colony morphology (Burgess et al. 2021), a potential confound we control for statistically. Across coral genera and reef systems, CAFI abundance increases with colony volume but less than proportionally (Abele & Patton 1976; Gotelli et al. 1985; Curtis et al. 2023; Brush et al. 2026), and a companion experiment found proportional scaling during initial colonisation but sublinear scaling in established communities (Stier et al. in review), suggesting that density-dependent processes -- including territorial exclusion -- shift the scaling regime as communities mature.

Together, these observations motivate a linked set of predictions connecting habitat quantity, through density dilution and community assembly, to host condition. We surveyed 114 *Pocillopora* colonies spanning three orders of magnitude in volume across three reef habitats on Mo'orea, French Polynesia, cataloguing 3,989 individuals from 243 taxa. First, we predicted that CAFI abundance scales sublinearly with colony volume, consistent with Propagule Redirection, and that the raw species--area relationship reflects a passive sampling artifact -- more individuals on larger corals yield more species by chance -- rather than active species accumulation beyond what abundance predicts. Second, if sublinear scaling promotes species stacking, composition should turn over along the size gradient through species replacement rather than nestedness, and colony volume should predict occupant composition even after controlling for abundance differences among corals. Third, we predicted that CAFI species richness covaries with coral physiological condition after accounting for colony size and site -- a prediction grounded in experimental evidence that CAFI mutualists enhance coral condition (Glynn 1983; McKeon et al. 2012; Holbrook et al. 2008) and in BEF theory, which holds that diverse assemblages deliver complementary services (Tilman et al. 2014). Because richness and abundance are strongly correlated under sublinear scaling, any observed association could reflect abundance variation, diversity per se, or both; we used variance partitioning and path analysis to evaluate these pathways. We also tested whether spatial proximity to neighbouring *Pocillopora* colonies predicts CAFI richness, as expected if colonists disperse among nearby patches (Hanski 1998).

If this chain holds -- even partially -- it predicts that climate-driven reductions in habitat-former size will disproportionately reduce associated biodiversity relative to habitat area lost. Large coral colonies are disproportionately vulnerable to marine heatwaves (Speare et al. 2022), and the scaling exponent determines how colony-size shifts translate into projected changes in associated biodiversity. Whether analogous scaling--assembly--condition dynamics operate in other biogenic-habitat systems -- oyster beds, kelp forests, forest canopies -- will depend on colonist-supply structure and occupant functional diversity in each system, parameters that targeted surveys can quantify.


## Methods

### Survey design

We surveyed 114 *Pocillopora* colonies across three reef habitats on Mo'orea, French Polynesia (17°30'S, 149°50'W) during June--August 2019 (Fig. 1A). Sites spanned the range of fore-reef environments at 3--8 m depth (mean 5.4 m): Hauru (fringing reef, 17.516°S 149.922°W; n = 38), Maatea (lagoon/back reef; n = 39), and Maharepa (barrier reef, 17.475°S 149.817°W; n = 37). At each site, divers haphazardly selected colonies across the full size range, targeting approximately equal numbers of small, medium, and large colonies to ensure coverage of the three-order-of-magnitude volume gradient. This size-stratified convenience sample overrepresents the tails of the natural size distribution; our scaling and composition results therefore reflect size effects conditional on this sampling design rather than population-level size-frequency relationships.

Colonies were identified to species via mtORF haplotyping (101 of 114 successfully genotyped; Johnston et al. 2022): *P. grandis* (n = 49), *P. meandrina* (n = 34), *P. tuahiniensis* (n = 7), and *P. verrucosa* (n = 10); 13 colonies failed amplification and were excluded from haplotype-dependent analyses. Colony size was estimated as ellipsoidal volume (V = 4/3 × π × L/2 × W/2 × H/2, where L, W, and H are maximum length, width, and height of the colony). Dimensions were measured in situ with a flexible ruler to the nearest centimetre. Where available, laboratory-measured volumes from water displacement (n = 93) replaced field estimates (n = 19); two colonies lacked all dimension data (n = 112 for size-dependent analyses). Volume ranged from 21 to 42,333 cm³ (Fig. 1B).

For the first 21 colonies at each site (63 total; 61 with volume data), we additionally counted all *Pocillopora* neighbours within a 5-m radius of the focal colony's centre and measured inter-colony distances (Fig. 1C). The remaining 51 colonies were surveyed for size and CAFI only (no neighbourhood census), reflecting logistical constraints on the time-intensive neighbour surveys. Colonies lacking neighbourhood data should not be interpreted as isolated; they simply were not surveyed for neighbours.

### CAFI collection

Each colony was enclosed in a fine-mesh bag (1 mm mesh) underwater. Clove oil (eugenol) was applied inside the bag to anaesthetise mobile fauna, the colony was detached from the substrate with a chisel, and the sealed bag was transported to the laboratory. In the laboratory, 68 of 114 colonies were additionally broken apart to extract fauna from deep within the branch matrix ("clove and smash" protocol); the remaining 46 were processed with clove oil extraction only. Macroboring and deeply cryptic fauna embedded within the coral skeleton were excluded. All specimens were preserved in 95% ethanol.

Specimens were identified to the lowest practical taxonomic level using regional taxonomic references (fishes: Randall 2005; decapods: Poupin 2010; molluscs: Kay 1979; general invertebrates: WoRMS), yielding 243 operational taxonomic units (OTUs) from 3,989 individuals (154 to species, 89 to genus or higher). Representative specimens were photographed and retained as vouchers. Species accumulation curves confirmed adequate sampling coverage (Fig. S1). We computed species richness, Shannon diversity, and rarefied richness (expected species at n = 20 individuals; vegan::rarefy; available for 68 of 112 colonies with ≥20 individuals).

### Coral physiological condition

We assayed four physiological traits from upper-branch tissue samples (one branch per colony): total protein (Bradford assay, BSA standard, 595 nm), total carbohydrate (anthrone method, glucose standard, 630 nm), zooxanthellae density (haemocytometer cell counts), and tissue biomass (ash-free dry weight: dried at 60°C for 48 h, combusted at 450°C for 4 h). All traits were normalised to coral surface area (cm²), determined by the paraffin wax-dipping method. Each trait was corrected for sampling position along the coral branch by extracting residuals from regressions on the lengths of the remaining branch stump and the detached tissue fragment (nubbin); position measures were uncorrelated with colony volume (both |r| < 0.15). Condition was summarised as PCA axis 1 of the four standardised (z-scored) residuals, computed with centering but without rescaling (residuals were pre-standardised). PCA axis 1 captured 61.6% of variance with all loadings positive (protein 0.55, carbohydrate 0.52, zooxanthellae density 0.48, AFDW 0.44). Of 114 colonies, 84 retained complete data for BEF analyses; missingness was predicted by site but not by volume, richness, or abundance (Table S14).

### Statistical analyses

All analyses used R v4.5.2 (R Core Team 2024) with key packages vegan 2.7.2, MASS 7.3.65, piecewiseSEM 2.3.1, DHARMa 0.4.7, boot 1.3.32, sandwich 3.1.1, betapart 1.6.1, and mediation 4.5.1 (full session information archived with data). Significance was assessed at α = 0.05 throughout. Natural logarithms were used for all log transformations; because the negative binomial and Poisson GLMs use a natural-log link, the log(volume) coefficient directly estimates the power-law scaling exponent without unit conversion. Site was included as a fixed effect because three levels provide too few for reliable random-effect variance estimation (Bolker et al. 2009); inferences therefore apply to these three sites specifically.

Multiple-testing corrections: Hochberg step-up (k = 2) for the two pre-specified BEF predictors (species richness and total CAFI abundance, chosen a priori based on the experimental evidence reviewed in the Introduction); Benjamini--Hochberg FDR (k = 4) for four exploratory predictors (Shannon diversity, *Trapezia* abundance, resident fish abundance, *Galeropsis monodonta* abundance); FDR across species in species-level tests. GLM diagnostics included DHARMa simulated residuals, Cook's distance (threshold 4/n), and variance inflation factors (VIF > 10 indicating serious multicollinearity; Dormann et al. 2013). Overdispersion in Poisson models was assessed via the Pearson χ²/df ratio, with automatic fallback to quasipoisson if the ratio exceeded 1.5. All 112 colonies with volume data were included in scaling and composition analyses, including colonies with zero CAFI (n = 0 colonies had zero CAFI, but the analytical pipeline handles this case). Sample sizes for specific analyses: 112 with volume, 84 with complete physiology, 68 with rarefied richness (≥20 individuals), 47 with both physiology and rarefied richness, 61 with neighbourhood surveys, and 39 with rarefied richness plus neighbourhood data.

### Q1: Scaling

CAFI abundance was modelled as a power law (N = aV^β^) using negative binomial GLMs with log link (total CAFI ~ log(volume) + site); the log(volume) coefficient directly estimates β. We tested β = 1 via Wald test and a bootstrap proportion test (1,000 site-stratified replicates with BCa confidence intervals; percentile CIs used as fallback when BCa acceleration failed; Efron 1987). Species richness was modelled with Poisson GLMs (Pearson χ²/df < 1.5; no overdispersion detected); the log(volume) coefficient yields the species--area relationship (SAR) exponent *z*. Rarefied richness (n = 20) was regressed on log(volume) via ordinary least-squares (OLS) regression. Scaling was repeated for 21 prevalent species (≥30 total individuals and ≥15% prevalence across colonies) and six taxonomic groups, with FDR correction within each category (Supplementary Methods).

### Q2: Composition

Compositional differences were tested using marginal PERMANOVA (community ~ log(volume) + site; Bray--Curtis distances on Hellinger-transformed abundances; vegan::adonis2, by = "margin"; 9,999 permutations). Hellinger transformation was applied prior to Bray--Curtis distance computation to down-weight dominant species while preserving abundance information (Legendre & Gallagher 2001). Dispersion homogeneity was assessed with PERMDISP (vegan::betadisper). Continuous size effects were tested with partial distance-based redundancy analysis (db-RDA; vegan::dbrda; same Bray--Curtis on Hellinger-transformed distances; site conditioned out), with variance partitioning (vegan::varpart) and species scores computed as weighted averages (WA scaling). Robustness was verified on iterated-rarefied distances (100 draws, rarefied to the minimum observed abundance per colony). Nestedness was tested with NODF (Almeida-Neto et al. 2008) against quasiswap null matrices (1,000 iterations; Miklós & Podani 2004). Beta diversity was decomposed into Sørensen-based turnover and nestedness-resultant components (betapart; Baselga 2010). Community variability was assessed with PERMDISP within equal-count volume terciles.

### Q3: Diversity and condition

We tested whether CAFI diversity predicts coral condition using OLS models (condition PC1 ~ predictor + log(volume) + site). Count predictors (abundance, functional group counts) were square-root-transformed to reduce leverage of high-abundance colonies while preserving zeros; richness was untransformed. To disentangle diversity from abundance effects, we used three complementary approaches: (1) partial regression (condition residuals after removing abundance + volume + site, regressed on richness residuals after the same removal, and vice versa); (2) hierarchical R² decomposition comparing full (richness + abundance) and reduced (each predictor alone) models; and (3) piecewiseSEM path analysis (Lefcheck 2016) with the following directed acyclic graph: volume → richness, volume → abundance, and richness + abundance + volume + site → condition. All predictors were z-scored for standardised path coefficients. Rarefied richness was compared to raw richness as a control for abundance confounding. Mediation analysis (1,000 bootstrap iterations; Tingley et al. 2014) tested whether the richness → condition effect operates through the abundance pathway (average causal mediation effect, ACME). Genetic species identity (mtORF haplotype) was added as a sensitivity covariate (Table S15). Reverse models tested condition as a predictor of CAFI metrics. Species-specific condition effects were tested for up to 10 key species (present on ≥5 colonies) with FDR correction across species.

### Neighbourhood context

For the 61 neighbourhood-surveyed colonies, we fitted full models with log(volume), neighbour count, log(total neighbour volume + 1), mean neighbour distance (colony centre to colony centre), and site. Abundance: negative binomial GLM; richness: Poisson GLM; Shannon diversity: OLS. Terms were dropped sequentially via AIC backward elimination (site protected from removal; terms dropped if removal reduced AIC), and the selected model was compared to the full model. The distance--richness relationship was retested with rarefied richness (n = 39). Supplementary analyses (co-occurrence null models, community assembly diagnostics, host species identity effects, and sensitivity analyses) are described in the Supplementary Methods.

### Ethics and permits

Research was conducted under permits issued by the Haut-commissariat de la République en Polynésie française, in collaboration with the Université de la Polynésie française (CRIOBE). Collections of invertebrates and small reef fishes did not require institutional animal care committee (IACUC) approval under University of California policy, which exempts invertebrates and cold-blooded vertebrates collected in the field.

### Data accessibility

Data and analysis code will be archived in the Dryad Digital Repository upon acceptance [DOI to be assigned]. Code is also available at [GitHub URL to be provided]. Raw specimen-level CAFI data (3,989 records), coral characteristics, and physiology measurements are included in the archived dataset.


## Results

### Q1: CAFI abundance and richness scale sublinearly with coral size

Total CAFI abundance scaled sublinearly with colony volume (n = 112; β = 0.52, 95% bootstrap CI [0.44, 0.62]; z_Wald = −11.45, p < 0.001 vs β = 1; Fig. 2A), consistent with Propagule Redirection: corals doubling in volume gained only ~43% more CAFI. Per-capita density declined accordingly: because density equals abundance divided by volume, a scaling exponent of 0.52 yields a density slope of −0.48 on the log--log scale (Fig. 2B). Species richness also scaled sublinearly (z = 0.34, 95% bootstrap CI [0.27, 0.42], z_Wald = −18.76, p < 0.001 vs z = 1; Fig. 2C), but rarefied richness showed no relationship with volume (p = 0.64, n = 68), indicating that the raw richness signal is tightly coupled to the abundance gradient.

Sublinear scaling extended to individual species: 11 of 21 prevalent species showed Propagule Redirection, 10 matched Field of Dreams, and none showed super-linear scaling (Fig. 3A,C; Fig. S4; Table S1). Territorial pair-bonders -- *Trapezia* (β = 0.43) and *Alpheus lottini* (a snapping shrimp; β = 0.37) -- consistently scaled sublinearly, while non-territorial species such as *Galeropsis monodonta* (β = 1.27, CI spanning 1.0) scaled proportionally.

Five of six taxonomic groups scaled sublinearly (β range 0.43–0.74; Fig. 3B,D); gastropods were the sole exception (β = 0.94, CI spanning 1.0). Size-dependent occurrence reinforced these patterns: 14 of 24 prevalent species showed significant occurrence–volume relationships (logistic GLM, FDR-adjusted p < 0.05). All results were robust to taxonomic resolution (Fig. S6). Cross-study comparison with the companion experiment showed concordant scaling classification for overlapping species; two species (*Alpheus lottini*, a snapping shrimp, and *Calcinus latens*, a hermit crab) shifted from proportional to sublinear scaling in established communities (Table S7).

### Neighbourhood context

Among the 61 colonies with neighbourhood surveys (5-m radius), closer neighbours predicted higher richness (β = −0.005, z = −3.20, p = 0.001; each additional metre of mean inter-colony distance was associated with 0.5% fewer species) and Shannon diversity (β = −0.007, t = −3.49, p = 0.001), but not total abundance (p = 0.78; Fig. S5; Table S4). Neither neighbour count nor total neighbour volume predicted any response (all p > 0.37). The distance–richness relationship survived rarefaction (β = −0.041, t = −2.97, p = 0.005, n = 39), confirming genuine species accumulation rather than a passive sampling artifact.

### Q2: Site pools and coral size structure community composition

Colony size structures not just how many species a coral hosts but which ones. Marginal PERMANOVA confirmed that both site and colony size structured CAFI composition (volume R² = 0.078, F₁,₁₀₈ = 9.74, p = 0.001; site R² = 0.060, F₂,₁₀₈ = 3.74, p = 0.001; n = 112; PERMDISP F₂,₁₀₉ = 0.89, p = 0.42), together explaining ~14% of variation. The site effect was robust across five distance metrics and balanced subsampling (Fig. S2). The three sites exhibited distinct taxonomic signatures (Fig. 4B): Maatea was characterised by hermit crabs (33% of CAFI), Maharepa by obligate symbionts (71%), and Hauru by fishes (12%).

Community composition also shifted continuously along the coral size gradient. Distance-based redundancy analysis (db-RDA; n = 112), after partialing out site effects, confirmed that coral volume explained 7.8% of compositional variation (F₁,₁₀₈ = 9.74, p = 0.001; Fig. 4A; this F-statistic equals the marginal PERMANOVA value because both tests partition the same Bray--Curtis distance matrix with respect to volume). The db-RDA size–composition gradient was robust to rarefaction (mean R² = 2.4% across 100 rarefaction draws, F₁,₁₀₈ = 2.64, p = 0.001), confirming genuine compositional turnover rather than an abundance artifact. This turnover reflects species stacking: territorial pair-bonders that exhibit the strongest sublinear scaling (e.g., *Trapezia punctimanus*, db-RDA score = −2.95) are associated with smaller corals, while other taxa replace them on larger colonies. Species vector fitting (envfit) confirmed *Trapezia punctimanus* (R² = 0.53), *Harpiliopsis beaupresii* (a palaemonid shrimp; R² = 0.30), and *Paragobiodon modestus* (a coral goby; R² = 0.25) as the strongest compositional drivers.

*Luniella pugil* (a xanthid crab; db-RDA score = −1.71) also loaded toward smaller corals, while *Euplica varians* (a columbellid snail; 1.06) and *Trapezia flavopunctata* (1.06) loaded toward larger corals. Two-way variance partitioning (volume and site) attributed 4.8% to volume alone, 3.9% shared, a negligible fraction to site alone (negative adjusted R²), and 91.5% residual. Compositional divergence among size classes was not significant after rarefaction (PERMDISP: p = 0.61; Fig. S3), suggesting that apparent convergence of large-coral assemblages reflects passive sampling rather than deterministic assembly.

Small-coral assemblages were not nested subsets of large-coral assemblages (NODF = 18.4, z = −1.09, p = 0.28; n = 112); most between-colony differences reflected species replacement, which accounted for 81% of total dissimilarity (Table S11). Combined with the significant db-RDA, this confirms species turnover -- not passive accumulation -- across colony sizes.

Before rarefaction, community variability declined continuously from small to large corals (distance-to-centroid ANOVA F₂,₁₀₉ = 9.21, p < 0.001; Fig. S15). Because the variability decline did not survive rarefaction (above), the decline likely reflects richer samples converging statistically rather than deterministic assembly, though size-mediated environmental filtering operating through abundance cannot be ruled out.

Supplementary community assembly analyses supported non-random structure: communities were more similar than expected under stochastic assembly (Raup-Crick; Fig. S16; Table S20), co-occurring CAFI were more taxonomically diverse than expected (net relatedness index; Fig. S17), and community similarity declined with geographic distance (Mantel r = 0.10, p = 0.031; partial Mantel controlling for size r = 0.11, p = 0.017; n = 57).

Because *Pocillopora* comprises multiple cryptic genetic species that differ in colony morphology and branching architecture (Johnston et al. 2022; Burgess et al. 2021), we tested whether host species identity contributes to compositional variation beyond colony size (Supplementary Methods). Coral genetic species predicted CAFI richness (likelihood ratio test, p = 0.007) and composition (marginal PERMANOVA R² = 0.083, p = 0.001; n = 98) but not total abundance (p = 0.24; Table S16). Three-way variance partitioning attributed comparable unique variation to host architecture (5.6%) and colony volume (4.7%), both exceeding site (2.9%; Fig. S16D; Table S20). The corallivore *Galeropsis monodonta* showed strong host specificity to *P. verrucosa* (7-fold enrichment), while obligate *Trapezia* as a genus did not differ among coral species (p = 0.51). Within *Trapezia*, however, individual species sorted by body size: across six species (n ≥ 10), host usage differed significantly (χ² = 53.98, df = 15, p < 0.001), with larger-bodied species concentrating on wide-branched *P. grandis* (Spearman r = 0.89, p = 0.019; Fig. 4C; Table S19).

### Pairwise co-occurrences are explained by volume and site

After accounting for coral volume and site, no species pair co-occurred more or less often than expected (0 of 528 pairs significant after FDR correction; Fig. S9; Table S5), indicating that pairwise associations are largely explained by habitat size and location. However, six species showed significant mating-pair aggregation -- concentrating conspecific pairs on a subset of corals (Table S6) -- as predicted by the territorial pair-bonding that underpins species stacking.

Because species replace rather than merely accumulate with increasing volume, corals differing in richness host partially non-overlapping assemblages -- a prerequisite for testing whether diversity predicts coral condition.

### Q3: CAFI species richness and coral condition

Species richness positively predicted coral condition (n = 84; standardised β = 0.36, 95% CI [0.06, 0.66], p = 0.018, Hochberg-adjusted p = 0.036; Fig. 5A), but the association was small (adjusted R² = 0.036, with the unique richness fraction explaining approximately 1% of condition variation) and did not survive abundance control: rarefied richness showed no association (β = −0.07, 95% CI [−0.29, 0.15], p = 0.50, n = 47; Fig. S8). Applying the same interpretive standard used in Q1 -- where the rarefied richness null was taken as evidence of passive sampling -- the most parsimonious interpretation is that the richness--condition association reflects abundance variation rather than species identity per se. Total abundance showed a weaker association (β = 0.32, p = 0.048, Hochberg p = 0.048; Fig. 5B; Table S2).

Variance partitioning attributed 29.1% of the model's incremental R² to richness uniquely (absolute ΔR² ≈ 0.010) and less than 1% to abundance uniquely, but neither fraction reached significance (richness F₁,₇₉ = 1.68, p = 0.20), reflecting the tight richness--abundance correlation (r ≈ 0.84; VIF = 6.2). Path model and mediation results are detailed in the supplement (Fig. S10; Tables S10, S12) but do not decisively separate complementarity from an abundance confound.

Adding genetic species identity (mtORF haplotype) as a covariate did not absorb the richness signal (p = 0.010, n = 74; Table S15), indicating that the richness--condition association is robust to host species identity.

The rarefied richness null result (reported above) should therefore be interpreted alongside the design-sensitivity calculations reported in the Supplementary Methods. Species-specific contributions to condition were tested for 10 key species (≥5 corals present); no species survived FDR correction (Table S3). Effect directions for functional groups and individual species broadly matched the companion experiment (binomial test; Table S3). No exploratory predictor survived FDR correction (Table S2), and reverse-direction tests were non-significant (Fig. S11C; Table S2).


## Discussion

Sublinear scaling in biogenic habitats does not merely reduce per-capita occupant density -- it restructures community composition through species replacement, generating qualitatively different assemblages on different-sized habitats. Our survey of 114 *Pocillopora* colonies documents this chain: abundance scales sublinearly with colony volume and per-unit density dilutes as corals grow (Q1), composition turns over through replacement rather than nestedness (Q2), and the resulting community variation covaries with coral condition -- though the most informative Q3 result is that this association is fragile (Figs. 2, 4, 5). We interpret our findings alongside the companion manipulation experiment (Stier et al. in review) and the broader literature on habitat-size effects. Because our cross-sectional design documents associations at each link but cannot establish causality, we flag interpretive ambiguities at each step below.

### Sublinear scaling and species stacking

Positive effects of host or habitat size on associated assemblages are widespread -- arthropods on trees (Lawton 1983; Schlinkert et al. 2015), invertebrates in tree holes (Srivastava & Lawton 1998), and fauna in kelp holdfasts (Anderson et al. 2005). Sublinear abundance scaling has been documented specifically in branching-coral systems (Abele & Patton 1976; Gotelli et al. 1985), reflecting finite colonist supply relative to available habitat. Our community-level exponent (β = 0.52; Fig. 2A) is broadly consistent with Abele & Patton's (1976) foundational estimate from *Pocillopora damicornis* (β ≈ 0.41 on a volume basis). Arthropod scaling exponents on trees are typically 0.2--0.5 on a crown-volume basis (Lawton 1983; Schlinkert et al. 2015), within this range. Srivastava & Lawton (1998) reported similarly sublinear richness scaling (z ≈ 0.20--0.30) for invertebrates in water-filled tree holes. The convergence of scaling exponents across these structurally disparate systems -- branching corals, tree canopies, phytotelmata -- suggests that propagule limitation relative to habitat volume may be a general constraint on occupancy in biogenic habitats, regardless of occupant identity or habitat architecture. Two reinforcing mechanisms support this pattern: colonization limited by propagule supply and landscape configuration rather than available habitat volume alone (Stier & Osenberg 2010; Hamman et al. 2018), and obligate symbionts imposing per-colony ceilings through territorial pair-bonding (Castro 1978; Huber & Coles 1986).

Increasing total CAFI on a colony therefore favours adding new species over additional conspecifics, reflecting species stacking, wherein intraspecific density ceilings imposed by territorial exclusion force diversification (Stier et al. 2012). Our cross-sectional design cannot separate size from age; larger colonies may be older, giving colonization and competitive exclusion more time to shape the scaling regime.

Taxon-specificity of scaling exponents reinforces this mechanistic picture. *Trapezia* scaled most strongly sublinearly, while gastropods were near proportional.

The companion experiment found proportional scaling during active colonization, while our survey shows the sublinear pattern, suggesting that territorial exclusion progressively reduces per-capita survival as communities mature.

At the landscape scale, spatial proximity to neighboring corals predicted higher richness and Shannon diversity but not total abundance. Nearby colonies likely serve as stepping-stone sources for rare taxa, and neighbor count and total neighbor volume had no detectable effect. Design-sensitivity calculations for this subset are reported in the Supplementary Methods.

### Composition, turnover, and co-occurrence

The rejection of nestedness indicates that species replace rather than merely accumulate across colony sizes, generating qualitatively different assemblages on different-sized colonies.

Compositional turnover survived rarefaction control, confirming genuine species replacement independent of passive sampling. Site and volume together explain only about 14% of composition, leaving substantial variation attributable to host identity and unmeasured environmental factors. Supplementary analyses point to host architecture as an additional driver, while pairwise co-occurrence is mostly size- and site-driven.

### Occupant diversity and coral condition

The richness--condition association survives genetic species control, but the signal remains entangled with abundance and host architecture. Supplementary analyses of variance partitioning, path models, mediation, and sensitivity checks are reported in Fig. S10, Fig. S11, and Tables S10--S15.

### Colony size, heatwaves, and the hypothesised chain

The three questions form a hypothesised chain: coral size structures occupant communities, composition turns over with increasing volume, and diversity shows a fragile association with condition. Whether these patterns persist across recruitment cycles and thermal stress events remains unknown. Per-capita density dilution could interact with size-dependent coral mortality during marine heatwaves, reducing biotic insurance when thermal stress demands it most. Longitudinal tracking of marked colonies would provide a direct test.


---

We thank the many people who contributed to establishing and maintaining the field surveys, as well as M. Brzezinski and D. Cryan for assistance with data management and laboratory work and C. Wall, C. Bove, and J. Baumann for guidance on coral physiological protocols. We are grateful to the staff of the University of California Gump Research Station, the Moorea Coral Reef LTER, and the UGA Osenberg Lab. Research was conducted under permits issued by the Territorial Government of French Polynesia and the Haut-commissariat de la Republique en Polynesie Francaise. Supported by NSF OCE-1851510 and OCE-1851032.

---

## Author Contributions

ACS and CWO conceived the project. AP and JSC conducted fieldwork and laboratory analyses. ACS performed analyses and drafted figures. ACS wrote the first draft; ACS and CWO revised.

---

## Conflict of Interest

The authors declare no conflicts of interest.

---

## Figure Legends

**Fig. 1.** Study design and representative CAFI. (A) Mo'orea reef sites: Hauru (n = 38), Maatea (n = 39), Maharepa (n = 37); 114 colonies total. (B) Colony volume distribution (log10 scale; 21--42,333 cm³; n = 112 with volume data). (C) Neighbourhood density (conspecifics within 5 m; n = 61). (D--F) Representative CAFI: *Trapezia* sp. guard crab, *Harpiliopsis spinigera* shrimp, *Neocirrhites armatus* hawkfish.

![](../output/figures/manuscript/fig1_study_design.png)

**Fig. 2.** CAFI scaling with coral volume (n = 112; log--log scale). (A) Total abundance scales sublinearly (β = 0.52 [0.44, 0.62]; dashed: proportional). (B) Per-capita density declines (slope = −0.48). (C) Species richness scales sublinearly (z = 0.34 [0.27, 0.42]). Points coloured by site.

![](../output/figures/manuscript/fig2_scaling.png)

**Fig. 3.** Species and group scaling. (A) Abundance vs volume for 10 prevalent species. (B) Four broad taxonomic groups. (C) Species-level β with 95% CI (blue: Redirection; grey: Field of Dreams). (D) Group-level exponents. Dashed: β = 1.

![](../output/figures/manuscript/fig3_species_group_scaling.png)

**Fig. 4.** Community composition and architectural filtering. (A) db-RDA biplot (site partialed); volume explains 7.8% (p = 0.001). Point size proportional to volume; ellipses: 80% concentration by site. (B) Taxonomic group relative abundance by site. (C) *Trapezia* guard-crab species composition on tight- vs wide-branched corals, showing body-size-mediated architectural sorting (χ² = 53.98, p < 0.001).

![](../output/figures/manuscript/fig4_composition.png)

**Fig. 5.** Occupant diversity and coral condition (n = 84). Condition = PC1 of four physiological traits. (A) Richness predicted condition (Hochberg p = 0.036). (B) Total abundance: weaker association (p = 0.048). Points by site; lines: marginal fits; shading: 95% CI.

![](../output/figures/manuscript/fig5_feedbacks.png)

---

## References

Abele, Lawrence G., and William K. Patton. 1976. "The Size of Coral Heads and the Community Biology of Associated Decapod Crustaceans." *Journal of Biogeography* 3: 35--47.

Almeida-Neto, M., et al. 2008. "A Consistent Metric for Nestedness Analysis in Ecological Systems." *Oikos* 117: 1227--39.

Anderson, M.J., et al. 2005. "Consistency and Variation in Kelp Holdfast Assemblages." *Journal of Experimental Marine Biology and Ecology* 320: 35--56.

Angelini, C., et al. 2011. "Interactions Among Foundation Species and Their Consequences for Community Organization, Biodiversity, and Conservation." *BioScience* 61: 782--89.

Baselga, A. 2010. "Partitioning the Turnover and Nestedness Components of Beta Diversity." *Global Ecology and Biogeography* 19: 134--43.

Bolker, B.M., et al. 2009. "Generalized Linear Mixed Models: A Practical Guide for Ecology and Evolution." *Trends in Ecology & Evolution* 24: 127--35.

Bronstein, J.L. 1994. "Conditional Outcomes in Mutualistic Interactions." *Trends in Ecology & Evolution* 9: 214--17.

Brush, E.G., et al. 2026. "Habitat Characteristics and Priority Effects Shape Fish and Invertebrate Assemblages Inhabiting the Coral *Pocillopora grandis* in Hawai'i." *Oikos*, e11247.

Burgess, S.C., et al. 2021. "Response Diversity in Corals: Hidden Differences in Bleaching Mortality Among Cryptic *Pocillopora* Species." *Ecology* 102: e03324.

Burns, K.C. and J. Dawson. 2005. "Patterns in the Diversity and Distribution of Epiphytes and Vines in a New Zealand Forest." *Austral Ecology* 30: 883--91.

Carr, M.H. and M.A. Hixon. 1997. "Artificial Reefs: The Importance of Comparisons with Natural Reefs." *Fisheries* 22: 28--33.

Castro, P. 1978. "Movements Between Coral Colonies in *Trapezia ferruginea*." *Marine Biology* 46: 237--45.

Chase, J.M., et al. 2011. "Using Null Models to Disentangle Variation in Community Dissimilarity from Variation in alpha-Diversity." *Ecosphere* 2: art24.

Chase, T.J., et al. 2018. "Coral-Dwelling Fish Moderate Bleaching Susceptibility of Coral Hosts." *PLoS ONE* 13: e0208545.

Connor, E.F. and E.D. McCoy. 1979. "The Statistics and Biology of the Species-Area Relationship." *The American Naturalist* 113: 791--833.

Counsell, C.W.W. and M.J. Donahue. 2021. "Protection Mutualists Affect Colonization and Establishment of Host-Associated Species in a Coral Reef Cryptofauna Community." *Oikos* 130: 1823--33.

Counsell, C.W.W., et al. 2018. "Variation in Coral-Associated Cryptofaunal Communities Across Spatial Scales and Environmental Gradients." *Coral Reefs* 37: 827--40.

Curtis, J.S., et al. 2023. "3D Photogrammetry Improves Measurement of Growth and Biodiversity Patterns in Branching Corals." *Coral Reefs* 42: 623--27.

Dayton, P.K. 1972. "Toward an Understanding of Community Resilience and the Potential Effects of Enrichments to the Benthos at McMurdo Sound, Antarctica." In *Proceedings of the Colloquium on Conservation Problems in Antarctica*, edited by B.C. Parker. Allen Press, pp. 81--95.

Dormann, C.F., et al. 2013. "Collinearity: A Review of Methods to Deal with It and a Simulation Study Evaluating Their Performance." *Ecography* 36: 27--46.

Efron, B. 1987. "Better Bootstrap Confidence Intervals." *Journal of the American Statistical Association* 82: 171--85.

Ellison, A.M., et al. 2005. "Loss of Foundation Species: Consequences for the Structure and Dynamics of Forested Ecosystems." *Frontiers in Ecology and the Environment* 3: 479--86.

Fonseca, C.R. and G. Ganade. 1996. "Asymmetries, Compartments and Null Interactions in an Amazonian Ant-Plant Community." *Journal of Animal Ecology* 65: 339--47.

Glynn, P.W. 1983. "Increased Survivorship in Corals Harboring Crustacean Symbionts." *Marine Biology Letters* 4: 105--11.

Gotelli, N.J., S.L. Gilchrist, and L.G. Abele. 1985. "Population Biology of *Trapezia* spp. and Other Coral-Associated Decapods." *Marine Ecology Progress Series* 21: 89--98.

Haddad, N.M., et al. 2015. "Habitat Fragmentation and Its Lasting Impact on Earth's Ecosystems." *Science Advances* 1: e1500052.

Hamman, E.A., et al. 2018. "Landscape Configuration Drives Persistent Spatial Patterns of Occupant Distributions." *Theoretical Ecology* 11: 111--27.

Hanski, I. 1998. "Metapopulation Dynamics." *Nature* 396: 41--49.

Hartig, F. 2022. *DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models*. R package version 0.4.6.

Hector, A., et al. 1999. "Plant Diversity and Productivity Experiments in European Grasslands." *Science* 286: 1123--27.

Hochberg, Y. 1988. "A Sharper Bonferroni Procedure for Multiple Tests of Significance." *Biometrika* 75: 800--802.

Holbrook, S.J., et al. 2008. "Effects of Sheltering Fish on Growth of Their Host Corals." *Marine Biology* 155: 521--30.

Huber, M.E. and S.L. Coles. 1986. "Resource Utilization and Competition Among the Five Hawaiian Species of *Trapezia*." *Marine Ecology Progress Series* 30: 21--31.

Johnston, E.C., R. Cunning, and S.C. Burgess. 2022. "Cophylogeny and Specificity Between Cryptic Coral Species (*Pocillopora* spp.) at Mo'orea and Their Symbionts (Symbiodiniaceae)." *Molecular Ecology* 31: 5368--85.

Lawton, J.H. 1983. "Plant Architecture and the Diversity of Phytophagous Insects." *Annual Review of Entomology* 28: 23--39.

Lefcheck, J.S. 2016. "piecewiseSEM: Piecewise Structural Equation Modelling in R." *Methods in Ecology and Evolution* 7: 573--79.

Loreau, M. and A. Hector. 2001. "Partitioning Selection and Complementarity in Biodiversity Experiments." *Nature* 412: 72--76.

MacArthur, R.H. and E.O. Wilson. 1967. *The Theory of Island Biogeography*. Princeton University Press.

McKeon, C.S., et al. 2012. "Multiple Defender Effects: Synergistic Coral Defense by Mutualist Crustaceans." *Oecologia* 169: 1095--103.

Miklós, I. and J. Podani. 2004. "Randomization of Presence--Absence Matrices: Comments and New Algorithms." *Ecology* 85: 86--92.

Oksanen, J., et al. 2022. *vegan: Community Ecology Package*. R package version 2.6-2.

Otway, S.J., A. Hector, and J.H. Lawton. 2005. "Resource Dilution Effects on Specialist Insect Herbivores in a Grassland Biodiversity Experiment." *Journal of Animal Ecology* 74: 234--40.

Poulin, R. 1997. "Species Richness of Parasite Assemblages: Evolution and Patterns." *Annual Review of Ecology and Systematics* 28: 341--58.

R Core Team. 2024. *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna, Austria.

Schlinkert, H., et al. 2015. "Plant Size as Determinant of Species Richness of Herbivores, Natural Enemies and Pollinators Across 21 Brassicaceae Species." *PLoS ONE* 10: e0135928.

Shima, J.S., C.W. Osenberg, and A.C. Stier. 2010. "The Vermetid Gastropod *Dendropoma maximum* Reduces Coral Growth and Survival." *Biology Letters* 6: 815--18.

Silliman, B.R., et al. 2005. "Drought, Snails, and Large-Scale Die-Off of Southern U.S. Salt Marshes." *Science* 310: 1803--6.

Speare, K.E., et al. 2022. "Size-Dependent Mortality of Corals During Marine Heatwave Erodes Recovery Capacity of a Coral Reef." *Global Change Biology* 28: 1342--58.

Srivastava, D.S. and J.H. Lawton. 1998. "Why More Productive Sites Have More Species: An Experimental Test of Theory Using Tree-Hole Communities." *The American Naturalist* 152: 510--29.

Stella, J.S., et al. 2014. "From Cooperation to Combat: Adverse Effect of Thermal Stress in a Symbiotic Coral-Crustacean Community." *Oecologia* 174: 1187--95.

Stier, A.C., et al. 2012. "Housekeeping Mutualisms: Do More Symbionts Facilitate Host Performance?" *PLoS ONE* 7: e32079.

Stier, A.C. and C.W. Osenberg. 2010. "Propagule Redirection: Habitat Availability Reduces Colonization and Increases Recruitment in Reef Fishes." *Ecology* 91: 2826--32.

Stier, A.C. and C.W. Osenberg. 2024. "Coral Guard Crabs." *Current Biology* 34: R5--7.

Stier, A.C., C. Primo, J.S. Curtis, and C.W. Osenberg. In review. "Habitat Quantity Drives Community Assembly and Biodiversity-Ecosystem Function Relationships in Coral Reef Fauna."

Tilman, D., F. Isbell, and J.M. Cowles. 2014. "Biodiversity and Ecosystem Functioning." *Annual Review of Ecology, Evolution, and Systematics* 45: 471--93.

Tingley, D., et al. 2014. "mediation: R Package for Causal Mediation Analysis." *Journal of Statistical Software* 59(5): 1--38.

Yachi, S. and M. Loreau. 1999. "Biodiversity and Ecosystem Productivity in a Fluctuating Environment: The Insurance Hypothesis." *Proceedings of the National Academy of Sciences* 96: 1463--68.

---

## Supplementary Materials

See separate document: `combined_supplement.md` (17 supplementary figures, 20 supplementary tables, supplementary methods).
