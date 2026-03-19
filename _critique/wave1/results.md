## Results Diagnosis

### Score: 29/40

### Rubric

| Dimension | Score (1–5) | Notes |
|-----------|-------------|-------|
| 1. Biology-first reporting | 4 | Most findings lead with biology; a few open with statistical apparatus (e.g., "Volume-weighted null models detected..."; "Species accumulation curves confirmed..."). |
| 2. Effect sizes + CIs | 4 | CIs reported consistently for key GLM exponents. Gaps: Q3 adjusted R² = 0.04 lacks CI; path-model path coefficients lack SEs; rarefied-richness slope lacks CI. |
| 3. Non-significant results | 4 | Null results for co-occurrence, nestedness, neighborhood, and rarefied richness are all reported — this is strong. Minor gap: reverse-direction FDR context unclear. |
| 4. Figure/table refs | 3 | Most figures are cited. Critical gaps: Fig. S10 (variance partitioning table/figure) referenced but the exact panel labels are inconsistent; Fig. S11C cited for exploratory null results but the reference comes after the result is stated, not inline. Several supplement tables (S1, S2, S4–S8) cited without the reader being guided to what they show. |
| 5. Parallel structure with Intro | 3 | Q1→Q2→Q3 order is maintained, which is good. However, the Intro also posits a density dilution PREDICTION (per-capita density declining) as a consequence of sublinear scaling — Results reports it but buries it mid-paragraph 1 of Q1 rather than as a separate beat. The Intro's "species stacking" prediction (small faunas not nested) is answered at the end of Q2, not at its opening. The Intro's Q3 explicitly invokes "complementarity among functionally distinct residents" — Results mentions this in Q3 but the complementarity framing is not the opening claim; instead the opening foregrounds correlation. |
| 6. Subheading quality | 3 | The three main subheadings are finding-based ("scale sublinearly," "structure community composition," "is associated with coral condition") — correct. But "Supporting results (Supplement)" and its two sub-subheadings ("Neighborhood effects," "Co-occurrence patterns") are method-named, not finding-based. "Cross-study scaling concordance" is also method-framed. |
| 7. No interpretation creep | 2 | Several passages belong in the Discussion: (a) "The most parsimonious interpretation is that the raw species–area relationship reflects passive sampling" — this is interpretation; (b) "illustrating that sublinear scaling necessarily produces density dilution" — mechanistic explanation; (c) "consistent with high stochasticity in community assembly" — interpretation of residual variance; (d) the path-model paragraph ends with "provide only weak directional support given the non-significant model fit" — this is inferential conclusion appropriate to Discussion; (e) "suggesting that pairwise associations are largely explained by volume and site" — explanatory claim. Some of this is in the Stier voice style ("integrate results with interpretation"), but the degree goes beyond orienting statements into Discussion-level explanation. |
| 8. No methods recapitulation | 2 | Multiple passages recapitulate Methods: (a) Q2 opens with "Species accumulation curves confirmed adequate sampling coverage" — this tests a methods assumption, not a biological finding; (b) "Multivariate dispersion did not differ among sites (PERMDISP: F = 0.89, p = 0.42), confirming that the PERMANOVA result reflects true compositional differences rather than dispersion heterogeneity" — this is a diagnostic check, not a result; (c) The spatial autocorrelation paragraph at the end of Q2 ("CAFI abundance showed no significant spatial autocorrelation...") is a methods validation, not a result about coral communities; (d) "This average treats species as independent estimates, although species co-occurring on the same corals are not strictly independent" — methods caveat embedded in Results. These should be one sentence each, if they appear at all, or moved to Supplement. |

### PASS/FAIL: FAIL (29/40; two dimensions below 3)

---

### Issues Found

1. **[RES-01]** severity: **critical**
   - Problem: Interpretation creep throughout Q1, Q2, and Q3 — Results explains mechanisms and makes inferential conclusions that belong in the Discussion. The clearest violations: "The most parsimonious interpretation is that the raw species–area relationship reflects passive sampling" (Q1 ¶1), "consistent with high stochasticity in community assembly" (Q2 ¶1), and the path-model paragraph's closing sentence: "The variance partitioning and the directional asymmetry in the path model are consistent with the interpretation that the richness signal is not merely a by-product of abundance, but provide only weak directional support given the non-significant model fit."
   - Location: Q1 ¶1 (sentence 4–5), Q2 ¶1 (last sentence), Q3 ¶3 (last sentence)
   - Proposed fix: In Q1, replace "The most parsimonious interpretation is that the raw species–area relationship reflects passive sampling: larger corals accumulate more individuals and therefore more species, but do not support higher diversity per capita. However, this test is interpretively ambiguous..." with: "Rarefied richness (expected species at n = 20 individuals) showed no relationship with volume (OLS: slope = 0.14, SE = 0.30, t₆₄ = 0.47, p = 0.64; colonies with <20 CAFI excluded; Fig. S8). We address the interpretation of this pattern in the Discussion." In Q2, replace "consistent with high stochasticity in community assembly" with "with 86% of compositional variation unexplained by the predictors measured." In Q3 ¶3, cut the final sentence ("The variance partitioning and the directional asymmetry...") — it belongs in Discussion and is already there.
   - Ripple risk: Discussion (some text would move there; avoid duplication — Discussion already contains the interpretive language so simply cutting from Results is sufficient).

2. **[RES-02]** severity: **critical**
   - Problem: Methods recapitulation in Results — three passages report methodological diagnostics rather than biological findings: (a) the PERMDISP sentence in Q2 ¶2 ("Multivariate dispersion did not differ..."), (b) the spatial autocorrelation paragraph at the end of Q2, and (c) the opening of Q2 ("Species accumulation curves confirmed adequate sampling coverage").
   - Location: Q2 ¶1 sentence 1; Q2 ¶2 sentence 3; Q2 ¶5 (entire final paragraph of Q2)
   - Proposed fix: (a) "Species accumulation curves confirmed adequate sampling coverage" — delete this sentence; it tests a methods assumption and the figure (Fig. S1) already conveys it. (b) The PERMDISP sentence should be compressed to a parenthetical: "Site effects were robust to dispersion heterogeneity (PERMDISP: F = 0.89, p = 0.42)." (c) The spatial autocorrelation paragraph should be moved to a one-sentence footnote or a supplementary methods validation note. If it must remain, reframe as a single parenthetical in the Q1 opening: "Spatial independence was confirmed (Moran's I = −0.004–0.024, all p > 0.28; Table S10)." Do not devote an entire paragraph to it in Q2.
   - Ripple risk: none (content goes to Supplement or is deleted; Fig. S1 and Table S10 citations are maintained).

3. **[RES-03]** severity: **major**
   - Problem: The Intro's Q3 predicts "complementarity among functionally distinct residents" generating a "diversity–condition relationship beyond the richness–abundance correlation." The Results Q3 opening sentence foregrounds the association but buries the complementarity framing. The opening claim is: "Corals harboring more CAFI species were in better physiological condition." This is biologically correct but does not directly answer the Intro's question, which asked whether richness predicts condition INDEPENDENTLY of abundance. The independent-of-abundance test (variance partitioning) is reported in ¶2, not ¶1. The reader reaches ¶1 wondering "but does this survive abundance control?" and the answer comes a paragraph later.
   - Location: Q3 ¶1 (opening sentence and paragraph structure)
   - Proposed fix: Restructure Q3 opening to lead with the independent-of-abundance finding: "CAFI species richness positively predicted coral physiological condition independently of abundance (partial regression: standardized β = 0.27, t₇₉ = 2.42, p = 0.018; n = 84; Fig. 5A), with variance partitioning attributing 29.1% of the CAFI-explained variance uniquely to richness and less than 1% uniquely to abundance (Fig. S10A). Total CAFI abundance was marginally associated with condition (standardized β = 0.32, t₇₉ = 2.01, p = 0.048; Fig. 5B). However, because richness and abundance were strongly correlated (r = 0.84; VIF = 6.2 and 6.8 respectively), variance partitioning explains only 7% additional variance beyond volume and site — most explanatory power was shared (70.8%)."
   - Ripple risk: Discussion ¶3 (Composition/BEF section) — consistent, no changes needed there.

4. **[RES-04]** severity: **major**
   - Problem: The PERMANOVA F-statistic is duplicated with an error — Q2 ¶1 reports "volume F₁,₁₀₈ = 9.74, R² = 0.08" and the db-RDA ¶ also reports "F₁,₁₀₈ = 9.74, p = 0.001." These are stated to be "algebraically identical" in Methods (and the CLAUDE.md confirms this), but the reader sees identical F-values for two apparently separate tests with different reported R² values (0.08 for PERMANOVA, 7.8% for db-RDA). Without explanation, this looks like a reporting error. The text should clarify this identity or use only one test at point of first report.
   - Location: Q2 ¶1 (PERMANOVA) and Q2 ¶3 (db-RDA)
   - Proposed fix: In the db-RDA paragraph, explicitly note: "The db-RDA F-statistic (F₁,₁₀₈ = 9.74, p = 0.001) is algebraically identical to the PERMANOVA marginal F reported above; the R² differs because PERMANOVA R² is unadjusted while db-RDA variance partitioning reports adjusted fractions." Or, since this is already stated in Methods, simply drop the redundant F from the db-RDA sentence: "db-RDA confirmed that coral volume explained 7.8% of compositional variation (p = 0.001; Fig. 4A)."
   - Ripple risk: Methods already contains the justification; no changes needed there.

5. **[RES-05]** severity: **major**
   - Problem: Subheadings for supplementary results are method-based, not finding-based. "Neighborhood effects," "Co-occurrence patterns," and "Cross-study scaling concordance" tell the reader nothing about what was found.
   - Location: "Supporting results (Supplement)" subsection headers
   - Proposed fix:
     - "Neighborhood effects" → "Neighborhood density does not predict CAFI abundance"
     - "Co-occurrence patterns" → "Pairwise CAFI co-occurrences are explained by volume and site"
     - "Cross-study scaling concordance" → "Survey and experimental scaling exponents are broadly concordant for seven species"
   - Ripple risk: none (headers only).

6. **[RES-06]** severity: **major**
   - Problem: The rarefied richness result in Q3 has an inconsistency between Q1 and Q3: Q1 reports the rarefied richness–volume test as n = 68 (t₆₄), slope = 0.14, p = 0.64. Q3 reports the rarefied richness–condition test as β = −0.07, p = 0.50, n = 47. These are different samples (n = 68 vs. n = 47) and different responses (volume vs. condition). The reader is not told why n differs; no CI is given for the Q3 result. Additionally, the Q3 figure citation for this result is "Fig. S8; Fig. S11B" — two separate figure references for one result, which is unusual. If both figures are genuinely needed, explain what each shows.
   - Location: Q1 ¶1 (rarefied richness vs. volume), Q3 ¶4 (rarefied richness vs. condition)
   - Proposed fix: Add a parenthetical explaining the n difference: "(n = 47; colonies with <20 CAFI and missing physiological data excluded; see Fig. S8 for rarefaction depth sensitivity)." Add CI: β = −0.07 [−0.28, 0.13], p = 0.50. Remove the duplicate Fig. S11B citation (Fig. S8 covers it, or consolidate).
   - Ripple risk: none.

7. **[RES-07]** severity: **major**
   - Problem: Q1 ¶1 mentions interpretive ambiguity of the rarefaction test with "see Discussion" but the same caveat is repeated verbatim in Q3 ¶4: "We address the interpretive ambiguity of this test in the Discussion." Identical hedge language in two separate Results paragraphs is redundant and weakens both sentences.
   - Location: Q1 ¶1 last sentence; Q3 ¶4 last sentence
   - Proposed fix: In Q1, cut the hedge entirely (the data speak; Discussion is the natural place for interpretation). In Q3, keep one brief forward pointer: "We return to this ambiguity in the Discussion." Alternatively, consolidate both observations into a single sentence at the start of the Q3 rarefied richness result, since that is where the causal question lands.
   - Ripple risk: none.

8. **[RES-08]** severity: **major**
   - Problem: The adjusted R² = 0.04 in Q3 lacks a CI and lacks a benchmark. For a BEF claim in JAE, reviewers will want to know the 95% CI on this R² and how it compares to effect sizes in the BEF literature. Also, adjusted R² for a model with 4 predictors (richness + volume + site₂ + site₃) over n = 84 is reported as 0.04, but the "incremental R² = 0.07 beyond volume and site" reported in ¶2 implies the full-model R² is higher. The reader is left with two R² values that appear contradictory (0.04 for the richness model vs. 0.07 incremental for CAFI predictors).
   - Location: Q3 ¶1 last sentence; Q3 ¶2 first sentence
   - Proposed fix: Clarify which model has R² = 0.04 (the richness-only forward model: condition ~ richness + log(volume) + site) versus which model's increment is 0.07 (the joint CAFI predictor model). State clearly: "The richness model (richness + log(volume) + site) explained adjusted R² = 0.04 of condition variance; adding abundance to this model explained an additional incremental R² = 0.03, for a total CAFI-attributable increment of R² = 0.07 beyond volume and site alone." If the numbers don't work out that way, verify against the script output.
   - Ripple risk: none (Discussion already acknowledges the modest effect; no changes needed there).

9. **[RES-09]** severity: **minor**
   - Problem: The opening sentence of Q2 ¶1 is a statistical frame, not a biological claim: "NMDS ordination (2 dimensions, stress = 0.157) showed site-level clustering along the first axis." NMDS stress is a methods diagnostic; it does not belong in the first sentence of a biology-first paragraph.
   - Location: Q2 ¶1 sentence 1
   - Proposed fix: Lead with the biological finding: "Reef site was the dominant structuring axis of CAFI community composition (Fig. 4A,B; PERMANOVA: site R² = 0.06, p = 0.001; volume R² = 0.08, p = 0.001; n = 112). NMDS stress = 0.157 indicates adequate ordination fit." Or bury the stress value in a parenthetical.
   - Ripple risk: none.

10. **[RES-10]** severity: **minor**
    - Problem: The Discussion opens with "Colony size is the dominant axis structuring CAFI communities." This exact framing is also implicit in the Results Q2 subheading ("Site pools and coral size structure community composition"), and the Discussion's synthesis repeats several Results sentences nearly verbatim. This redundancy is partly a Stier voice feature ("integrate results with interpretation") but reviewers may flag it as filler.
    - Location: Discussion ¶1 (overlap with Results Q1–Q3 summary)
    - Proposed fix: This is borderline — the issue is in the Discussion, not the Results. Flag for Discussion critic. In Results, ensure no sentence ends with an interpretive summary that pre-empts the Discussion's opening synthesis.

11. **[RES-11]** severity: **minor**
    - Problem: Missing CI for the weighted-mean β across 21 species in Q1 ¶3. The text reports β = 0.51 [0.45, 0.56] — this is present. Good. However, the Wald z-test is also reported (z = −18.8, p < 0.001) for this weighted mean, while the independence caveat follows immediately. The structure is: result → caveat → next sentence. The caveat is important but placed mid-paragraph disrupts flow.
    - Location: Q1 ¶3
    - Proposed fix: Move the independence caveat to a parenthetical: "...significantly below 1.0 (z = −18.8, p < 0.001; note species co-occur on the same corals so SE of weighted mean may be underestimated; Fig. S4)." This preserves the caveat without breaking paragraph flow.
    - Ripple risk: none.

12. **[RES-12]** severity: **minor**
    - Problem: The Intro's Q2 prediction states "small-coral faunas should not be nested subsets of large-coral faunas." The Results answers this at the end of Q2 ¶4 (nestedness result), not at the opening of Q2. The nestedness answer should come earlier — it is the direct answer to the Intro's compositional prediction — and be linked to the db-RDA result as its complement.
    - Location: Q2 ¶4 (nestedness result)
    - Proposed fix: Move the nestedness finding to immediately follow the db-RDA paragraph, and frame it as the answer to the Q2 prediction: "Consistent with the species-stacking prediction, community nestedness along the size gradient was not significant (NODF = 18.4, z = −1.09, p = 0.28), confirming that small-coral faunas are not subsets of large-coral faunas. Combined with significant compositional turnover in db-RDA, this pattern indicates species replacement rather than passive accumulation as colony size increases."
    - Ripple risk: none.

---

### Diagnostic Checklist

**From write-results skill — "When critiquing" section:**

- [x] Every Introduction question answered? — Q1, Q2, Q3 all answered. The density-dilution sub-prediction (per-capita density declines) is answered but not foregrounded. The "species stacking → no nesting" prediction is answered at the end of Q2 rather than its opening.
- [ ] Findings stated as biology first? — MOSTLY yes; violations at Q2 ¶1 opener (NMDS stress), Q2 ¶2 opener (species accumulation curves), Q2 ¶2 sentence 3 (PERMDISP as finding). See RES-02 and RES-09.
- [ ] Effect sizes AND CIs reported? — CIs present for β exponents. Missing: (a) CI for rarefied richness–condition β; (b) CI for adjusted R² = 0.04; (c) SE or CI for path model path coefficients β = 0.55 and β = 0.02. See RES-06, RES-08.
- [x] Non-significant results reported? — Yes: rarefied richness (both tests), co-occurrence, nestedness, neighborhood, reverse BEF, exploratory predictors, individual species. Strong.
- [ ] Interpretation creeping in? — Yes, in multiple locations. See RES-01 (critical). The passages "most parsimonious interpretation," "consistent with high stochasticity," and the path-model closing sentence are the three clearest violations.
- [x] Figure/table references present at point of relevance? — Mostly yes; a few late references (Fig. S11C mentioned after result; supplement table citations could be more informative). Adequate but not optimal.

**Additional checks (scientific-writing-principles):**

- [ ] One idea per paragraph? — Q3 ¶1 covers two distinct findings (richness association AND abundance association AND their correlation AND the overall R²). Overloaded. See RES-03.
- [x] Active voice dominant? — Yes, predominantly active throughout.
- [ ] Consistent terminology? — "CAFI abundance" and "total CAFI abundance" used interchangeably. "Compositional variation" used for both raw % explained and adjusted R² fractions without distinguishing.
- [x] Parallel structure Q1→Q2→Q3 maintained? — Yes at the section level; broken within sections (species stacking prediction answered at end of Q2, not opening; see RES-12).
- [ ] No methods recapitulation? — Violated (RES-02): PERMDISP, species accumulation curves, spatial autocorrelation paragraph.
- [x] Subheadings biology-based for main sections? — Yes for Q1–Q3; no for supplement subsections (RES-05).
- [ ] No duplicate identical hedges? — Violated (RES-07): "interpretive ambiguity" caveat appears twice.

**Cross-section parallelism checks:**

- [x] Intro Q1 (scaling) → Results Q1 (scaling): aligned. β values answer the FoD vs. Redirection framing explicitly.
- [ ] Intro Q2 (composition: turnover vs. nesting) → Results Q2: the nesting answer is buried at end rather than leading (RES-12).
- [ ] Intro Q3 (richness predicts condition INDEPENDENTLY of abundance) → Results Q3: independence is not the opening claim (RES-03).
- [x] Abstract numerics match Results: β = 0.52 [0.44, 0.62] ✓; p = 0.018 ✓; r = 0.84 ✓; 243 taxa ✓.
- [x] Methods → Results order: Q1 (scaling) → Q2 (composition) → Q3 (condition) matches Methods section order.
- [x] Discussion opens with correct synthesis of Results: all three Q findings referenced in Discussion ¶1. ✓
