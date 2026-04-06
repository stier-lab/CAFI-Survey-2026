# Figure Index

**Colony size drives sublinear scaling, compositional turnover, and biodiversity--condition covariation in coral-associated fauna**

Stier, A.C. et al. — *Journal of Animal Ecology*

---

## Main Text Figures

| Fig | File | Panels | Description | Source Script |
|-----|------|--------|-------------|---------------|
| 1 | `fig1_study_design.png` | A--F | **(A)** Satellite map of Mo'orea showing three reef sites (Hauru, Maatea, Maharepa). **(B)** Distribution of *Pocillopora* colony volumes (log10 scale, 21--42,333 cm^3, n = 112). **(C)** Distribution of neighborhood density (neighbors within 5 m, n = 61). **(D--F)** Representative CAFI species: *Trapezia* sp. guard crab, *Harpiliopsis spinigera* shrimp, *Neocirrhites armatus* flame hawkfish. | `01_load_data.R` |
| 2 | `fig2_scaling.png` | A--C | **(A)** Total CAFI abundance vs. colony volume (NB GLM, beta = 0.52 [0.44, 0.62]; sublinear = Propagule Redirection). **(B)** Per-capita CAFI density vs. volume (slope = -0.48; density dilution). **(C)** Species richness vs. volume (Poisson GLM, z = 0.34 [0.27, 0.42]). Three vertical panels, points colored by site. | `05_species_scaling_analysis.R` |
| 3 | `fig3_species_group_scaling.png` | A--D | **(A)** Abundance vs. volume for 10 most prevalent species (NB GLM curves). **(B)** Abundance vs. volume for 4 taxonomic groups (crabs, shrimps, fish, snails). **(C)** Species-level scaling exponents (beta +/- 95% bootstrap CI, n = 21 species). **(D)** Group-level scaling exponents. Colors: blue = Redirection, grey = Field of Dreams, vermillion = super-linear. | `05_species_scaling_analysis.R` |
| 4 | `fig4_composition.png` | A--B | **(A)** db-RDA biplot (Hellinger-transformed; volume explains 7.8%, p = 0.001). Points sized by coral volume, 80% site ellipses, top-5 species vectors. **(B)** Relative abundance of taxonomic groups at each site. | `02_community_analysis.R` |
| 5 | `fig5_feedbacks.png` | A--B | **(A)** Species richness vs. coral condition PC1 (p = 0.018). **(B)** Total CAFI abundance (sqrt-scaled) vs. condition (p = 0.048). Points colored by site, linear fits with 95% CI. n = 84 colonies with physiology data. | `09_cafi_condition_feedbacks.R` |

All main text figures are saved in `output/figures/manuscript/` as both PNG (300 dpi) and PDF (cairo_pdf). Legend text files: `output/figures/manuscript/figN_legend_results.txt`.

---

## Supplementary Figures

| Fig | File | Panels | Description | Source Script |
|-----|------|--------|-------------|---------------|
| S1 | `figS1_species_accumulation.png` | 1 | Species accumulation (rarefaction) curves by site. All sites approach asymptotic richness; Hauru has highest total richness. | `02_community_analysis.R` |
| S2 | `figS2_permanova_sensitivity.png` | A--B | **(A)** PERMANOVA R^2 for site and volume across 5 distance metrics. **(B)** Composition divergence trend robustness across metrics. Site effect significant for all 5 metrics; volume for all except Gower. | `02_community_analysis.R` |
| S3 | `figS3_composition_divergence.png` | A--B | **(A)** Distance-to-centroid by size class (raw Bray--Curtis). **(B)** After rarefaction (100 draws): divergence trend disappears (PERMDISP p = 0.61) -- abundance artifact. | `02_community_analysis.R` |
| S4 | `figS4_species_scaling.png` | 1 | Species-level scaling forest plot (21 species). beta +/- 95% CI. 11 Redirection, 10 Field of Dreams, 0 super-linear. FDR-corrected within species category. | `05_species_scaling_analysis.R` |
| S5 | `figS5_neighborhood.png` | A--F | **(A,C,E)** Total CAFI abundance, species richness, and Shannon diversity vs. number of neighbors; **(B,D,F)** the same responses vs. mean neighbor distance (n = 61). Neighbor count effects were null, but shorter mean neighbor distance predicted higher richness and Shannon diversity. | `04_landscape_effects.R` |
| S6 | `figS6_taxonomy_sensitivity.png` | 1 | Forest plot of 7 metrics across 5 taxonomy scenarios (baseline, species-only, merge-up, lump-down, rare-excluded). All results robust: abundance beta = 0.515--0.527; richness z = 0.310--0.343. | `13_taxonomy_sensitivity.R` |
| S7 | `figS7_functional_groups.png` | A--B | **(A)** Taxonomic group scaling exponents (6 groups). 5/6 sublinear; gastropods the sole exception (beta = 0.94, n.s.). **(B)** Group composition by site. | `08_functional_groups.R` |
| S8 | `figS8_rarefaction_sensitivity.png` | 1 | Rarefied richness vs. condition at 5 rarefaction depths (n = 10, 15, 20, 25, 30). Non-significant at all depths (all p > 0.10). | `09_cafi_condition_feedbacks.R` |
| S9 | `figS9_cooccurrence.png` | A--C | **(A)** Pairwise SES heatmap (volume-weighted null, 10,000 iter). 1/528 pairs significant after FDR. **(B)** Intraspecific density: 6 species show mating-pair aggregation (FDR < 0.05). **(C)** Size-dependent co-occurrence across 3 size classes. | `06_cooccurrence_analysis.R` |
| S10 | `figS10_bef_variance_partitioning.png` | A--C | **(A)** Variance decomposition: 29.1% unique to richness, <1% to abundance, 70.8% shared. **(B)** Partial regression of condition on richness residuals. **(C)** Path model coefficients (richness beta = 0.55, abundance beta = 0.02). | `09_cafi_condition_feedbacks.R` |
| S11 | `figS11_condition_details.png` | A--F | **(A)** A priori BEF forest plot. **(B)** Rarefied richness vs. condition (p = 0.50). **(C)** Exploratory group forest plot. **(D--E)** *Trapezia* and *Galeropsis* scatterplots. **(F)** Bidirectional analysis (CAFI <-> condition). | `09_cafi_condition_feedbacks.R` |
| S12 | `figS12_occurrence_curves.png` | 24 | Logistic GLM occurrence probability vs. log(volume) for 24 species. 14/24 significant after FDR; all positive (larger corals = higher occurrence). | `05_species_scaling_analysis.R` |
| S13 | `figS13_species_trait_heatmap.png` | A--B | **(A)** Standardized beta for 19 species x 5 condition traits. **(B)** FDR-adjusted p-values. No species x trait test survives FDR -- consistent with complementarity (BEF). | `09_cafi_condition_feedbacks.R` |
| S14 | `figS14_species_trait_biplots.png` | A--I | 9 strongest species x trait associations (raw p < 0.10). Position-corrected traits vs. sqrt(abundance). Includes *Breviturma pica* x protein, *Paracirrhites arcatus* x AFDW, *Periclimenes* sp. x zooxanthellae, others. | `09_cafi_condition_feedbacks.R` |
| S15 | `figS15_beta_dispersion.png` | 1 | Distance-to-centroid by size class. Small > Medium > Large (ANOVA F = 9.21, p < 0.001). Communities converge with coral size. | `15_community_assembly.R` |
| S16 | `figS16_community_assembly.png` | A--D | **(A)** Raup-Crick histogram (mean = 0.19, null = 0.5). **(B)** Raup-Crick by size class (violins). **(C)** Beta-dispersion by size class. **(D)** Variation partitioning bar chart (architecture 5.6%, size 4.6%, space 2.9%). | `15_community_assembly.R` |
| S17 | `figS17_taxonomic_structure.png` | 1 | SES.MPD (NRI) vs. log(volume). Taxonomic overdispersion: mean NRI = −1.18 (p < 10⁻¹³). Trend with volume suggestive but NS (β = −0.17, p = 0.11). | `15_community_assembly.R` |

All supplement figures are saved in `output/figures/supplement/` as both PNG (300 dpi) and PDF (cairo_pdf).

---

## Summary

| Category | Count |
|----------|-------|
| Main text figures | 5 (14 panels total) |
| Supplement figures | 17 (46+ panels total) |
| Supplement tables | 12 |
| CSV output files | ~66 |
| Legend text files | 5 main + 3 supplement |
