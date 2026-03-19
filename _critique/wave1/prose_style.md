## Prose Style Diagnosis

### Score: 31/40

### Rubric

| Dimension | Score (1-5) | Notes |
|-----------|-------------|-------|
| 1. Active voice ratio (>70% non-Methods) | 4 | Introduction (3 passives/1000 words), Results (~5), Discussion (~2.5) all within acceptable range; Methods correctly passive-heavy |
| 2. Hedging control (<3/1000 words, no double-hedges) | 3 | Discussion hits 3.5 hedges/1000 words, just above threshold; no double-hedges found; "consistent with" used 14× across paper as a soft hedge substitute |
| 3. Filler word density | 3 | "rather than" used 12× in Discussion alone (overuse creates a hedged, contrastive refrain); no "it is important to note," no "in order to"; isolated "very" and "quite" instances are minor |
| 4. Sentence variety (median ~25 words, range 8–45) | 2 | All four sections fall below the 25-word median target: Introduction median 17, Methods 19, Results 21, Discussion 18; writing runs consistently short and declarative without the long analytical sentences that balance punches |
| 5. Paragraph structure (claims, not announcements; no one-sentence paragraphs) | 3 | Five weak Discussion topic sentences; two genuine one-sentence paragraphs (Methods code display Para, Results spatial autocorrelation Para); Q3 Results section header uses passive ("is associated with") |
| 6. Transitions quality | 3 | No mechanical "Additionally/Furthermore" openers (good); Discussion Para 12 opens with structural announcement ("The three questions we address form...") rather than a conceptual bridge; Discussion Para 13 opens with a known fact ("Size-dependent coral mortality during marine heatwaves is well-documented") instead of the new claim |
| 7. Conciseness | 3 | "rather than" repeated 12× in Discussion creates a rhythmic tick; "consistent with" used 14× across paper (9 in Results + Discussion) as a hedged claim substitute; "These results are consistent with a hypothesized feedback loop" (Discussion P1) delays the claim the reader wants |
| 8. Voice consistency with Stier patterns | 4 | "we" absent from Results (0.5/1000 words — effectively absent); Discussion uses "we" appropriately; named hypotheses present; cross-system analogies strong; honest epistemic signaling present; "We posit" and "Here, we" are underused in Discussion; no banned words |

### PASS/FAIL: FAIL (score 31/40, Dimension 4 below 3)

---

### Issues Found

1. **[PROSE-01]** severity: **critical**
   - Problem: The Results section uses "we" only once across 1,934 words (0.5/1000 words). Most paragraphs are written in the third person or use passive constructions, stripping agency and voice: "Total CAFI abundance scaled sublinearly," "Multivariate dispersion did not differ," "Community nestedness along the size gradient was not significant," "Variance partitioning clarified this ambiguity." The Stier voice rule is "we" throughout.
   - Location: Results, every paragraph
   - Proposed fix: Re-lead each paragraph with the investigator's interpretive frame: "We found that total CAFI abundance scaled sublinearly..." / "Variance partitioning clarified this ambiguity: we attribute 29.1% of the explainable increment uniquely to richness..." / "We detected no significant pairwise co-occurrences after FDR correction..."
   - Ripple risk: none — stylistic change only, no content shift

2. **[PROSE-02]** severity: **critical**
   - Problem: Sentence lengths run persistently short. Median = 17 words (Introduction), 19 (Methods), 21 (Results), 18 (Discussion). The Stier target is ~25 words. The writing has developed a staccato reflex: short declarative sentences pile up without analytical sentences that weave mechanism into finding. Example from Discussion: "Species richness predicted coral condition (Fig. 5A), but the BEF signal is modest and entangled with abundance. Most explanatory power was shared between richness and abundance (70.8%) rather than uniquely attributable to either, and the overall effect was small (adjusted R² = 0.04), leaving open the possibility that the richness signal is an abundance artifact." These could be synthesized into a single analytical sentence that carries the interpretive weight.
   - Location: All sections
   - Proposed fix: Identify sequences of 3+ short declarative sentences and fuse them into one long analytical sentence. Target at least 4–6 sentences per section with 40+ words to expand the range. Example rewrite: "Species richness predicted coral condition (Fig. 5A), but the signal is modest — adjusted R² = 0.04 — and most explanatory power is shared between richness and abundance (70.8%), leaving unresolved whether the richness effect is a genuine diversity signal or an artifact of the richness–abundance correlation (r = 0.84)."
   - Ripple risk: none

3. **[PROSE-03]** severity: **major**
   - Problem: "consistent with" appears 14 times across Results (6) and Discussion (7). This construction is the paper's primary hedge — it asserts pattern but disclaims mechanism. Used once or twice, it signals appropriate caution; used 14 times, it becomes a verbal tic that undermines confidence and produces tautological sentences: "these results are consistent with a hypothesized feedback loop" (Discussion P1) and "this temporal shift is consistent with density-dependent processes emerging over time" (Discussion P5).
   - Location: Results throughout; Discussion Para 1, 2, 4, 5, 8, 9, 11
   - Proposed fix: Replace with direct claims where the data support them; reserve "consistent with" for genuinely ambiguous cases (2–3 uses maximum).
     - "These results are consistent with a hypothesized feedback loop" → "Together, these results outline a feedback loop in which coral size shapes occupant communities that, in turn, may affect coral performance."
     - "This temporal shift is consistent with density-dependent processes emerging over time" → "This temporal shift suggests that density-dependent territorial exclusion progressively reduces per-capita survival on larger, more occupied colonies."
   - Ripple risk: none

4. **[PROSE-04]** severity: **major**
   - Problem: "rather than" appears 12 times in the Discussion alone, creating a repetitive contrastive refrain: "composition turns over rather than nests," "directional rather than quantitative agreement," "proportional to radius squared rather than volume," "adding new species rather than more individuals," "genuine species replacement rather than a passive sampling artifact," "species contributing small, additive effects rather than a single dominant mutualist," etc. The construction is legitimate in each instance but its frequency creates a monotonous rhetorical pattern.
   - Location: Discussion throughout
   - Proposed fix: Cut or rephrase at least half. Many of the "X rather than Y" constructions can become positive claims: "composition turns over along the size gradient — small-coral faunas are not nested subsets of large-coral faunas" (already stated explicitly later); "the community-level finding reflects the summed contributions of many species, not a single dominant mutualist."
   - Ripple risk: none

5. **[PROSE-05]** severity: **major**
   - Problem: Five Discussion topic sentences are announcements or fact-statements rather than interpretive claims:
     - P9: "Species richness predicted coral condition (Fig. 5A), but the BEF signal is modest and entangled with abundance." — This is a Results re-statement, not a Discussion claim.
     - P12: "The three questions we address form a hypothesized feedback in which coral size structures CAFI communities (Q1)..." — structural meta-commentary, not a new claim.
     - P13: "Size-dependent coral mortality during marine heatwaves is well-documented (Speare et al. 2022)..." — opens with known background, buries the new hypothesis.
     - P14: "These scaling results also have implications for coral restoration." — "have implications for" is the classic filler setup.
     - P15: "Three lines of work would resolve the key ambiguities." — structural announcement.
   - Location: Discussion paragraphs 9, 12, 13, 14, 15
   - Proposed fixes:
     - P9: "Richness and abundance both predict coral condition, but the richness signal is modest, partially shared with abundance, and mechanistically ambiguous — it cannot be attributed to complementarity without a manipulative experiment."
     - P12: "Density dilution has a direct physiological implication: per-capita mutualist services decline as colonies grow, potentially amplifying size-dependent vulnerability to thermal stress."
     - P13: "We hypothesise that density dilution of beneficial CAFI represents an underappreciated contributing factor to the disproportionate mortality of large corals during marine heatwaves."
     - P14: "Restoration programmes that optimise colony density for growth while ignoring CAFI communities will underestimate the service benefits accruing to small, isolated outplants."
     - P15: "Three manipulative experiments would resolve what this observational survey cannot."
   - Ripple risk: none

6. **[PROSE-06]** severity: **major**
   - Problem: Results section header "## Q3: CAFI species richness is associated with coral condition" uses passive construction and the weak verb "is associated with." Section headers are the highest-profile sentences in the paper — reviewers read them first.
   - Location: Results, Q3 header
   - Proposed fix: "## Q3: CAFI species richness positively predicts coral physiological condition"
   - Ripple risk: none

7. **[PROSE-07]** severity: **minor**
   - Problem: Two genuine one-sentence paragraphs in Results: (1) "CAFI abundance showed no significant spatial autocorrelation (Moran's I = −0.004, z = 0.24, p = 0.51), nor did richness (I = 0.024, p = 0.28) or Shannon diversity (I = −0.039, p = 0.77; Table S10), indicating that spatial non-independence does not confound our analyses." — This single sentence sits as a standalone paragraph. (2) In Methods, the model formula display "condition_PC1 ~ CAFI_predictor + log(volume) + site" is isolated as a paragraph.
   - Location: Results Q2 section (last paragraph); Methods Q3 subsection
   - Proposed fix: For Results: append to preceding paragraph about PERMANOVA robustness, or fuse with the following Q3 paragraph using a bridge ("With spatial structure ruled out, we turn to the diversity–condition question..."). For Methods: integrate the formula into the surrounding sentence rather than isolating it.
   - Ripple risk: none

8. **[PROSE-08]** severity: **minor**
   - Problem: Discussion Para 8 opens with "Yet site and volume together explain only ~14% of composition, with 86% residual..." The "Yet" opener signals a pivot from the previous paragraph, but the previous paragraph (Para 7) makes a positive claim about nestedness and species turnover — the pivot to residual stochasticity is abrupt. The reader has not been set up to expect a qualification.
   - Location: Discussion Para 8
   - Proposed fix: The last sentence of Para 7 should bridge to Para 8: add a sentence at the end of Para 7 such as "This turnover is real, but colony identity accounts for only a fraction of the variance, with most composition driven by processes beyond size." Para 8 can then open with the stochasticity claim as a positive assertion: "Stochastic colonisation dominates: site and volume together explain only ~14% of composition..."
   - Ripple risk: none

9. **[PROSE-09]** severity: **minor**
   - Problem: Discussion Para 10 opens "We cannot resolve this ambiguity with our observational design." This is a legitimate honest-epistemic statement, but leading a paragraph with a negative capability claim buries the actual content of the paragraph (which provides evidence favouring complementarity). The paragraph should lead with what we can say, not what we cannot.
   - Location: Discussion Para 10
   - Proposed fix: "Several lines of evidence favour complementarity over a pure abundance pathway, though the observational design cannot resolve the ambiguity." Then develop the evidence. Close with: "A manipulative experiment holding abundance constant while varying richness is required to establish the causal link."
   - Ripple risk: none

10. **[PROSE-10]** severity: **minor**
    - Problem: Discussion Para 6 ("The relative importance of Field of Dreams versus Propagule Redirection depends on ambient settlement intensity...") is only 3 sentences and 69 words. It introduces a strong predictive claim — that reef systems with higher larval supply should show exponents closer to 1.0 — but abandons it immediately without developing the mechanism or naming a testable prediction. For a 3-sentence paragraph in the Discussion, this is close to under-development.
    - Location: Discussion Para 6
    - Proposed fix: Develop this paragraph to 5–6 sentences: explain why Mo'orea is in the Redirection regime (lower larval supply relative to habitat capacity), cite any settlement-rate data for Mo'orea if available, and sharpen the comparative prediction.
    - Ripple risk: none

---

### Per-Section Audit

#### Introduction
- Passive count: 3 (~2.7/1000 words; acceptable)
  - "This link has been recognized theoretically..." — could be active but conventional usage is fine
  - Two instances of "passive sampling" are content terms, not passive constructions
- Hedges: 2 (1.8/1000 words; within target)
  - "juveniles may redistribute among patches"
  - "community variation along the habitat size gradient could itself influence habitat condition"
- Fillers: 2 minor ("rather than" ×1; "very" appears in a citation parenthetical, not prose)
- Mechanical transitions: none
- Notes: Introduction is the strongest section for voice. Topic sentences are genuine claims. "Here, we" pivot is present. Named hypotheses used throughout. Cross-system breadth (trees, kelps, oysters) is a Stier signature move.

#### Methods
- Passive count: 57 (25.4/1000 words) — expected and appropriate for Methods
- Hedges: 2 (0.9/1000 words) — well-controlled
  - "BCa acceleration constants could not be estimated" — appropriate conditional
  - "the SE of the weighted mean may be underestimated" — appropriate caveat
- Fillers: 1 minor ("very" appears in "FDR correction" running text, not as an intensifier)
- Mechanical transitions: none
- Notes: Methods are correctly passive in construction ("were extracted," "were computed"). The use of "We tested," "We repeated," and "We used" for analytical steps is good. Three one-sentence fragments exist: the formula display for the Q3 model (isolated paragraph), "All code is available at [repository URL]" (standalone), and "Species with ≥10 occurrences were included" (standalone). The last two could be folded into adjacent sentences.

#### Results
- Passive count: 8 (~4.1/1000 words; slightly elevated)
  - "Maatea was characterized by high hermit crab abundance" — active alternative: "Maatea supported high hermit crab abundance"
  - "Maharepa was dominated by obligate coral symbionts" — active: "Obligate coral symbionts dominated at Maharepa"
  - "Key species...were tested individually" — active: "We tested key species individually"
  - "should be interpreted cautiously" — acceptable passive
  - "CAFI species richness is associated with coral condition" (header) — must be fixed (see PROSE-06)
- Hedges: 2 (1.0/1000 words)
  - "rarefaction could also strip out a genuine size–diversity pathway" — appropriate
  - "null results for density may reflect insufficient power" — appropriate
- Fillers: "rather than" ×4 — each is locally justified but the pattern is building toward overuse
- Mechanical transitions: none
- Notes: Critical issue is the near-absence of "we" (1 instance in 1,934 words). Every paragraph reads as a self-narrating data summary rather than an investigator making and explaining findings. This is the single most important fix in the manuscript.

#### Discussion
- Passive count: 5 (~2.5/1000 words; acceptable)
  - "At low settlement, the propagule pool is exhausted quickly" — acceptable (process description)
  - "Depth partitioning...has been documented" — active option: "Counsell et al. (2018) documented depth partitioning..."
  - "Most explanatory power was shared between richness and abundance" — active: "Richness and abundance shared most explanatory power (70.8%)"
  - "A manipulative experiment...would be required to establish the causal link" — acceptable modal construction
  - "As biogenic habitats are restructured by disturbance and restoration" — acceptable gerund construction
- Hedges: 7 (3.5/1000 words — **above 3/1000 threshold**)
  - "those communities may feed back to affect coral performance" (Para 1)
  - "Several post-settlement processes could drive the transition" (Para 5)
  - "potentially less favourable balance of effects" (Para 7)
  - "whose diversity may feed back on coral condition (Q3)" (Para 12)
  - "potentially flipping the interaction from mutualistic" (Para 12)
  - "could reduce biotic insurance under thermal stress" (Para 13)
  - "the potential to feed back on habitat condition" (Para 16)
  - Note: Several of these are appropriate (the feedback loop is genuinely untested); but two or three can be converted to direct claims.
- Fillers: "rather than" ×12 — see PROSE-04
- "consistent with" ×7 — see PROSE-03
- Mechanical transitions: none
- Banned words: none found

---

### Stier Voice Checklist

- [x] "we" used throughout — **FAIL in Results** (0.5/1000 words); adequate in Introduction, Methods, Discussion
- [x] Named hypotheses present — Field of Dreams and Propagule Redirection used consistently
- [x] Cross-system analogies — arthropods on trees, invertebrates in kelp holdfasts, parasites on hosts, ants in acacia thorns (Discussion closer is excellent)
- [x] Honest epistemic signaling — present in Discussion ("this remains a hypothesis rather than an established finding," "We cannot resolve this ambiguity," "plausibility is not evidence")
- [x] Quantitative claims in Discussion — present and strong (β = 0.52, 2^0.52^ ≈ 1.43, 2^−0.48^ ≈ 72%, R² ≈ 0.09)
- [x] No banned words — confirmed ("novel," "groundbreaking," "unprecedented," "paradigm shift" all absent)
- [ ] "We posit" underused — Discussion proposes multiple mechanisms without using this phrase; would strengthen the speculation-signaling
- [ ] Future directions developed — three priorities are described but the paragraph structure is a numbered list in narrative disguise; each priority could be a full developed paragraph rather than 2–3 sentences

---

### Priority Fixes (in order of impact)

1. **PROSE-01** — Add "we" back into Results. This is a 30-minute pass that transforms the voice of the entire section.
2. **PROSE-05** — Rewrite 5 Discussion topic sentences as claims (see rewrites above).
3. **PROSE-02** — Extend sentence length variation: target 4–6 sentences >40 words in Results and Discussion by fusing adjacent short declaratives into analytical sentences.
4. **PROSE-03 + PROSE-04** — Cut "consistent with" from 14 to ≤4; cut "rather than" in Discussion from 12 to ≤5.
5. **PROSE-06** — Fix Q3 Results section header.
6. **PROSE-07 through PROSE-10** — Minor structural fixes.
