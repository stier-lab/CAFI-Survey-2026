## Methods Diagnosis

### Score: 33/40

### Rubric

| Dimension | Score (1-5) | Notes |
|-----------|-------------|-------|
| 1. Study system description | 5 | Site coordinates, dates, habitat context, morphotype caveat, non-random selection rationale — all present. Colony volume range stated. |
| 2. Design clarity | 4 | Two-survey-type design (neighborhood n=61 vs. size-only n=53) is explained for the neighborhood subset but the size-only n=53 is never stated explicitly; only "final n=112" and "61 colonies" are given. Sampling unit (colony) is unambiguous. Exclusion of 2 MRB colonies is documented. |
| 3. Collection protocols | 5 | Bag extraction, ethanol preservation, no anesthetic, exclusion of boring fauna, taxonomic resolution levels — all stated. Volume formula spelled out. |
| 4. Statistical descriptions | 3 | Core models are well-described. Three problems: (a) the multiple-testing scheme is scattered across two paragraphs and partially self-referential; (b) the per-capita density "slope" appears in Results as a regression result but Methods says it was only "computed to visualize"; (c) the F-statistic identity between PERMANOVA volume and db-RDA (both = 9.74) is noted in Methods but the explanation is buried, risking reader confusion when Results presents them as separate results. |
| 5. Sample sizes | 3 | Most N values stated correctly. Two gaps: (a) rarefied richness–condition analysis uses n=47 in Results but Methods never explains the drop from 84→47; (b) the "size-only" subset (n=53) is never named as a distinct analysis group even though it is implicitly excluded from all neighborhood analyses. |
| 6. Software citation | 3 | R v4.4 cited. vegan (Oksanen et al. 2022) cited with author-year but no version number. piecewiseSEM (Lefcheck 2016) cited. sandwich, lmtest, DHARMa, car packages used in analyses (confirmed from CLAUDE.md and script descriptions) but not cited in Methods at all. |
| 7. Analysis–results order match | 4 | Q1 → Q2 → Q3 order consistent with Introduction and Results. Co-occurrence and neighborhood are correctly labeled as supplement analyses. One ordering ambiguity: Moran's I spatial autocorrelation appears at the end of Q2 Results but the Methods section for spatial autocorrelation is placed under "Sensitivity analyses," not under Q2. |
| 8. Data accessibility | 1 | Both the repository DOI and code URL are placeholders ("[repository DOI]", "[repository URL]"). Package versions are incomplete. This is a CRITICAL gap for submission to JAE, which requires data archiving per ARRIVE/data policy. |

---

### PASS/FAIL: FAIL

Score = 33/40, but Dimension 8 (Data Accessibility) = 1/5, which is below the minimum-3 threshold. Additionally, Dimensions 4 and 5 are borderline and carry ripple risk into Results.

---

### Issues Found

---

**[METH-01]** severity: **critical**
- Problem: Both data and code placeholders remain unfilled. Methods states "All analyses were conducted in R v4.4 (session information including package versions archived at [repository DOI])" and closes with "Data and analysis code are available at [repository URL]." Neither placeholder has been replaced with an actual DOI or URL.
- Location: Statistical analyses preamble paragraph; Data accessibility section (final line).
- Proposed fix: Replace both placeholders before submission. If archiving is pending, use the standard JAE language: "Data and code will be deposited in the Dryad Digital Repository (DOI to be assigned upon acceptance) and are available from the corresponding author on reasonable request." Do not submit with bare brackets.
- Ripple risk: Abstract and combined_manuscript.md will also need these updated.

---

**[METH-02]** severity: **major**
- Problem: The rarefied richness–condition analysis uses n=47 in Results ("β = −0.07, p = 0.50, n = 47") but Methods never explains the reduction from 84 to 47. The 84-colony BEF dataset is defined ("84 colonies were retained for BEF analyses") but the additional exclusions for the rarefied richness analysis are undocumented. Likely cause: colonies with <20 CAFI are excluded (analogous to the rarefied richness–volume analysis which states "colonies with <20 CAFI excluded"). This must be stated in Methods.
- Location: Q3 statistical analysis subsection, rarefied richness test paragraph. Also check: "To test whether the richness → condition relationship was an abundance artifact, we compared models using raw versus rarefied richness..."
- Proposed fix: Add one sentence after the rarefied richness comparison sentence: "Rarefied richness was undefined for colonies with <20 CAFI individuals, reducing the available sample from 84 to 47 for this comparison."
- Ripple risk: Results Q3 rarefied richness paragraph is already correct (states n=47); no Results change needed.

---

**[METH-03]** severity: **major**
- Problem: Per-capita CAFI density is described in Methods as computed "to visualize the density dilution prediction" — implying a visualization aid only. But Results presents it as a regression result: "Per-capita CAFI density declined with colony size (log–log slope = −0.48)." This implies a regression was fitted; the Methods should either (a) state that the log–log slope was estimated via OLS or (b) clarify that the −0.48 value is the algebraic complement of β (which is what the Results text actually does say, parenthetically). If no separate regression was run, the Methods language "computed to visualize" is accurate — but then Results should not present the slope as if it were an independently fitted parameter.
- Location: Q1 statistical analysis subsection, last sentence of first paragraph: "Per-capita CAFI density (individuals cm⁻³) was computed to visualize the density dilution prediction."
- Proposed fix (option A — preferred, since Results already explains the algebra): Replace the Methods sentence with: "Per-capita CAFI density (individuals cm⁻³; total CAFI / colony volume) was computed per colony; the log–log slope of density on volume equals β − 1 by definition and was not independently modeled." This aligns Methods with the explanation already in Results.
- Proposed fix (option B): If an OLS regression was actually fitted, state it: "We regressed log(per-capita density) on log(volume) + site using OLS to quantify density dilution."
- Ripple risk: Results Q1 second paragraph is consistent with option A; no change needed if option A is adopted.

---

**[METH-04]** severity: **major**
- Problem: The multiple-testing correction scheme is described in two places (preamble paragraph and Q3 subsection) but the Q3 subsection says "Correction procedures for the a priori and exploratory predictors follow the three-tier scheme described above" — a self-referential cross-reference that forces the reader to scroll back. The preamble also uses "Hochberg's step-up procedure (k=2)" for the a priori BEF predictors. This is non-standard terminology: the standard correction for k=2 a priori tests is simply Holm (which is equivalent to Hochberg when k=2 and more familiar to JAE reviewers). If Hochberg is correct, it should be justified; if Holm is intended, use the standard name. More importantly, the number of tests per family is inconsistently stated: "k=2" for a priori BEF, "k=4" for exploratory, but then a third tier for species-level tests is introduced without a stated k.
- Location: Statistical analyses preamble: "We used natural logarithms throughout; site was included as a fixed effect...The two a priori BEF predictors (richness and total abundance) were corrected with Hochberg's step-up procedure (k = 2); exploratory CAFI predictors with Benjamini–Hochberg FDR (k = 4); key species tests with FDR correction across all species tested."
- Proposed fix: Rewrite the preamble multiple-testing sentence for clarity and JAE audience: "We applied three levels of multiple-testing correction: (1) Holm correction for the two pre-specified BEF predictors (richness and total abundance; k = 2); (2) Benjamini–Hochberg FDR for four exploratory CAFI predictors (Shannon diversity, Trapezia, fish, and Galeropsis abundance; k = 4); and (3) Benjamini–Hochberg FDR across all species tested in species-level analyses. Significance threshold α = 0.05 throughout." Then remove the redundant cross-reference sentence in Q3.
- Ripple risk: Results Q3 refers to "BH-FDR correction" and "pre-specified a priori hypotheses" — verify these labels match after any revision.

---

**[METH-05]** severity: **major**
- Problem: Key packages used in analyses are not cited: DHARMa (simulated residuals), sandwich/lmtest (HC3 robust SEs and Breusch–Pagan tests), car (VIF), and mediation (bootstrap mediation — computed but not reported per CLAUDE.md). These are all non-base R packages whose authors deserve citation and whose version numbers are needed for reproducibility. The Methods preamble states "session information including package versions archived at [repository DOI]" — but this is currently a placeholder (METH-01) and even if filled, JAE convention expects at-use citations for major analytical packages.
- Location: Statistical analyses preamble; Q3 subsection ("Breusch–Pagan tests confirmed homoscedasticity...HC3 heteroscedasticity-consistent standard errors (sandwich package)").
- Proposed fix: Add citations for DHARMa (Hartig 2022), sandwich (Zeileis et al. 2020), lmtest (Zeileis & Hothorn 2002), and car (Fox & Weisberg 2019) at first use. These should be added to references.bib. The current text says "(sandwich package)" without author/year — replace with "(sandwich; Zeileis et al. 2020)".
- Ripple risk: references.bib will need 3–4 new entries.

---

**[METH-06]** severity: **minor**
- Problem: The size-only survey subset (n=53 colonies that had volume data but no neighborhood surveys) is never explicitly named or sized in Methods. Methods only introduces "a subset of 61 colonies ('neighborhood corals')" but does not state that the remaining 53 colonies are "size-only" corals with no neighborhood data — a distinction that matters because it explains why neighborhood analyses use n=61 while all other analyses use n=112. A reader unfamiliar with the design would not know the complement of 112−61=51 (or 112−61=51 with 2 excluded, implying size-only n=53 before exclusions) represents a distinct survey type.
- Location: Coral measurements section: "For a subset of 61 colonies ('neighborhood corals'), all Pocillopora neighbors within a 5-m radius were counted..."
- Proposed fix: Add one sentence: "The remaining 53 colonies were surveyed for CAFI and size only ('size-only corals') and are not included in neighborhood analyses." This makes the two-group structure explicit.
- Ripple risk: Neighborhood effects supplement section in Methods already correctly states n=61; no change needed there.

---

**[METH-07]** severity: **minor**
- Problem: The 21-species scaling analysis and the 24-species occurrence probability analysis use different prevalence cutoffs. Scaling uses "≥30 total individuals and ≥15% prevalence" (21 species) while logistic occurrence models use "≥15% prevalence" only (24 species). The difference is explained by the "≥30 individuals" requirement for scaling (which would exclude very rare but commonly present species), but this discrepancy is not explained in Methods and will draw reviewer attention. Results also notes "21 prevalent species" vs. "24 prevalent species" without explaining why the sets differ.
- Location: Q1 statistical analysis subsection: "We repeated the scaling analysis for each of the 21 most prevalent individual species (≥30 total individuals and ≥15% prevalence)" and "we fit logistic GLMs...for 24 species with ≥15% prevalence."
- Proposed fix: Add a parenthetical after the occurrence analysis sentence: "(These 24 species satisfy the ≥15% prevalence threshold; the 21-species scaling set requires the additional criterion of ≥30 total individuals, excluding three species that are prevalent but sparse in counts.)" Alternatively, define both sets upfront: "For species-level analyses, we defined two prevalence thresholds: 24 species with ≥15% occurrence prevalence (used for logistic occurrence models) and 21 species that additionally met a count minimum of ≥30 total individuals (used for scaling GLMs)."
- Ripple risk: Results Q1 references both 21 and 24 — once Methods is clarified, Results language can remain as is.

---

**[METH-08]** severity: **minor**
- Problem: Spatial autocorrelation (Moran's I) is described under the "Sensitivity analyses" subsection in Methods but appears in Results under Q2 (community composition), presented as confirming spatial independence for the entire dataset ("CAFI abundance showed no significant spatial autocorrelation...Table S10"). Placing Moran's I under sensitivity analyses implies it is a robustness check for one analysis; in practice it is a data-level quality control that precedes all spatial analyses. This ordering mismatch causes a minor structure break between Methods and Results.
- Location: Sensitivity analyses subsection, last sentence: "Spatial autocorrelation was assessed using Moran's I on CAFI abundance, richness, and Shannon diversity (Table S10)."
- Proposed fix: Move the Moran's I sentence from "Sensitivity analyses" to the end of the Statistical analyses preamble paragraph, where it would read alongside the other global quality controls (DHARMa, Cook's D, VIF). Or add it as a final sentence to the Q2 subsection since that is where Results reports it. Either placement is defensible; the key fix is that it should not be buried under "Sensitivity analyses."
- Ripple risk: The supplementary Table S10 reference would carry over unchanged.

---

**[METH-09]** severity: **minor**
- Problem: The Ethical statement says "IACUC...where applicable." This hedge is vague and reviewers at JAE, which requires explicit ethical compliance statements, may flag it. Coral collection that does not involve vertebrates may not require IACUC review — if true, state this directly: "Collection of invertebrate specimens did not require IACUC oversight; coral tissue sampling was conducted under the collection permits described above." If IACUC approval was obtained, state the protocol number.
- Location: Ethical statement subsection: "with oversight from the University of California, Santa Barbara Institutional Animal Care and Use Committee (IACUC) where applicable."
- Proposed fix: Either remove the IACUC clause if invertebrates are not covered by UCSB IACUC, or replace with: "Coral tissue collection for physiological assays was conducted under UCSB IACUC protocol [NUMBER]."
- Ripple risk: None.

---

**[METH-10]** severity: **minor**
- Problem: The position correction for physiological traits is described clearly but the rationale for using residuals (rather than, e.g., ANCOVA or standardization) is not given. A sentence of justification is warranted because this is a non-obvious analytical choice that directly affects all Q3 inferences: "We regressed each trait on stump length and nubbin length and used the residuals." Why residuals and not the raw values or ANCOVA? The Methods states the goal (remove tissue-gradient sampling artifact) but not the statistical rationale for this specific approach.
- Location: Coral physiological condition section: "We regressed each trait on stump length and nubbin length and used the residuals."
- Proposed fix: Add a justification clause: "...and used the residuals as position-corrected trait values, an approach equivalent to partial regression on the remaining predictors (Cohen et al. 2003)." Or cite a published protocol that uses this method. If the method is from a prior Stier lab paper, cite it.
- Ripple risk: None; Q3 analyses are unaffected.

---

### Diagnostic Checklist

**Design**
- [x] System, site, coordinates, dates
- [x] Sampling/experimental design (replication, controls described)
- [ ] Sampling unit defined — colony is unambiguous, but the two survey types (size-only vs. neighborhood) not explicitly distinguished as separate groups with distinct N
- [x] Sample sizes (total n=114, per site n=38/39/37, neighborhood n=61 stated)
- [ ] n=53 size-only corals never named

**Collection**
- [x] All variables defined with units (volume formula, CAFI count, physiology traits)
- [x] Protocols described (bag extraction, ethanol, no anesthetic, boring fauna excluded)
- [x] Deviations: 2 MRB colonies excluded — documented
- [x] Temporal extent: June–August 2019

**Analysis**
- [x] Model type for each analysis (NB GLM, Poisson, OLS, logistic, PERMANOVA, db-RDA, piecewiseSEM)
- [x] Response and predictor variables named for each model
- [x] Transformations: Hellinger for community, log for volume, sqrt for count predictors — documented with justifications
- [ ] Multiple testing correction fully clear — Hochberg vs. Holm terminology ambiguous; self-referential cross-reference in Q3 (METH-04)
- [ ] Software versions incomplete — R v4.4 stated; vegan cited without version; DHARMa, sandwich, lmtest, car not cited (METH-05)
- [x] Significance threshold: α = 0.05 stated
- [x] Model diagnostics described (DHARMa, Cook's D, VIF, Shapiro–Wilk, Breusch–Pagan)
- [ ] Missing data handling: n=84 for BEF defined, but n=47 for rarefied richness–condition not explained (METH-02)
- [ ] Order matches Results — Moran's I placed under Sensitivity in Methods but reported under Q2 in Results (METH-08)

**Data Accessibility**
- [ ] Repository and accession number — PLACEHOLDER only (METH-01)
- [ ] Code availability — PLACEHOLDER only (METH-01)
