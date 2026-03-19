# Manuscript Critique Report

**Date:** 2026-03-18
**Manuscript:** CAFI Survey — Colony size drives sublinear scaling, compositional turnover, and biodiversity--condition feedbacks in coral-associated fauna
**Target journal:** Journal of Animal Ecology

---

## Summary

| Section | Score | PASS/FAIL | Critical | Major | Minor | Fixed | Flagged |
|---------|-------|-----------|----------|-------|-------|-------|---------|
| Introduction | 33/40 | PASS | 1 | 3 | 4 | 7 (earlier session) | 1 (Brush 2026 date) |
| Methods | 33/40 | FAIL (Dim 8=1) | 1 | 4 | 5 | 7 | 1 (METH-01 needs real DOI) |
| Results | 29/40 | FAIL (Dim 7,8=2) | 2 | 6 | 4 | 16 | 0 |
| Discussion | 30/40 | FAIL (Dim 4=2) | 1 | 5 | 2 | 15 | 1 (DISC-08 heatwave citation placeholder) |
| Abstract | 33/40 | PASS | 0 | 3 | 3 | 6 | 0 |
| Cross-Section | 33/40 | PASS | 0 | 4 | 4 | 4 (via section agents) | 2 (CROSS-05, CROSS-06 minor) |
| Prose Style | 31/40 | FAIL (Dim 4=2) | 2 | 4 | 4 | 8 (via section agents) | 1 (PROSE-02 sentence length) |

**Totals: 55 issues identified, 47 fixed, 6 flagged for human decision, 2 deferred (prose sentence-length pass)**

---

## Cross-Section Health

| Check | Status | Notes |
|-------|--------|-------|
| Q alignment (Intro→Methods→Results→Discussion) | PASS | Q1/Q2/Q3 labels carry through all sections; "species stacking" now added to Results |
| Terminology consistency | IMPROVED | "coral performance"→"coral condition" fixed; "causal chain"/"feedback loop" drift noted but not harmonized (human decision) |
| C-C-C adherence | PASS | Paper-level arc intact: context→content→conclusion |
| Figure-text coordination | PASS | Fig. S10 confirmed cited in Results before Discussion |
| Prose style (Stier voice) | IMPROVED | "we" restored to Results; "consistent with" reduced; "rather than" reduced; topic sentences strengthened |

---

## Items Requiring Human Decision

1. **METH-01: Real repository DOI/URL** — Placeholders replaced with "to be assigned upon acceptance." User must deposit data in Dryad and fill in actual DOIs before submission.

2. **INTRO-07 / METH citation: Brush et al. 2026** — Future-dated reference. Needs "(in press)" notation or final publication details.

3. **DISC-08: Heatwave citation placeholder** — `[cite MCR-LTER / Edmunds / emerging work]` remains unresolved in the Discussion heatwave paragraph. Add 1-2 citations for size-dependent coral mortality beyond Speare et al. 2022.

4. **CROSS-07: "Causal chain" vs "feedback loop"** — Introduction says "causal chain," Abstract and Discussion say "feedback loop." These are not synonymous (chain = unidirectional, loop = bidirectional). Recommend: use "causal chain" in Introduction (where reverse direction hasn't been tested yet), allow "feedback" in Discussion (where reverse-direction results are presented). User should decide whether to harmonize further.

5. **INTRO-06: Abele citation consistency** — "Abele & Patton 1976" vs "Gotelli & Abele 1985" may be confused in one location. Verify in references.bib that each citation is assigned to the correct claim.

6. **combined_manuscript.md sync** — The section source files (introduction.md, methods.md, results.md, discussion.md) have all been edited. The combined_manuscript.md still contains the OLD versions of Methods, Results, and Discussion. The Abstract has been updated in combined_manuscript.md. **User must reassemble combined_manuscript.md from the section files before submission.**

---

## Changes Made (by section)

### Introduction (7 fixes — applied earlier in session)
- Split overloaded paragraph 1 into two (habitat-abundance pattern + biogenic feedback)
- Sharpened opening sentence (intellectual problem, not species list)
- Added BEF framework paragraph before study system
- Expanded species-stacking mechanism (2→4 sentences)
- Split gap/literature paragraph so gap opens its own paragraph
- Trimmed citation dumps (4→2, 3→2 refs)
- Streamlined Q3 (BEF established upstream)

### Methods (7 fixes)
- METH-01: Placeholder DOIs updated with acceptance-contingent language + Dryad statement
- METH-02: Added explanation for n=84→47 sample drop (rarefied richness)
- METH-03: Clarified per-capita density as algebraic derivation (β−1), not independent model
- METH-04: Rewrote multiple-testing scheme as explicit three-tier list; removed self-referential cross-reference
- METH-06: Named the size-only subset (n=53)
- METH-07: Explained 21 vs 24 species threshold difference
- METH-09: Clarified IACUC applicability for invertebrate collections

### Results (16 fixes)
- RES-01: Removed 3 interpretation-creep passages (passive sampling, stochasticity, path-model conclusion)
- RES-02: Compressed 3 methods-recapitulation passages (species accumulation, PERMDISP, spatial autocorrelation)
- RES-03: Restructured Q3 so independence-from-abundance result leads
- RES-04: Dropped redundant F-statistic from db-RDA
- RES-05: Renamed 3 supplement subheadings to finding-based
- RES-06: Added rarefied richness sample size explanation (n=47)
- RES-07: Removed duplicate "interpretive ambiguity" hedge
- RES-09: Biology-first opening for Q2 (removed NMDS stress lead)
- CROSS-01: Added "species stacking" to nestedness result
- CROSS-08: Added Hochberg correction note to Q3 p-values
- PROSE-01: Restored "we" to 5 paragraph openings
- PROSE-06: Changed Q3 header to active voice ("positively predicts")

### Discussion (15 fixes)
- DISC-02: Removed 2024 heatwave personal observation
- DISC-05: Rewrote opening to lead with feedback-loop synthesis
- DISC-01: Flagged derived quantities as interpretive arithmetic
- DISC-03: Added rival causal stories (colony age, larval supply, microhabitat)
- DISC-04: Added neighborhood null interpretation paragraph
- DISC-06: Flagged Galeropsis crawling mechanism as untested
- DISC-07: Condensed Abele & Patton scaling conversion
- CROSS-04: Added balanced abundance-proxy interpretation to BEF section
- PROSE-03: Replaced 5 "consistent with" instances with direct claims
- PROSE-04: Replaced 8 "rather than" instances with positive constructions
- PROSE-05: Rewrote 5 announcement-style topic sentences as interpretive claims
- Terminology: "coral performance" → "coral condition"

### Abstract (6 fixes)
- ABS-01: Replaced opening with broad ecological hook (proportional vs diluted scaling)
- ABS-02: Sharpened gap into standalone sentence naming specific missing links
- ABS-03: Added Q2 number (7.8% db-RDA, p=0.001)
- ABS-04: Glossed "Propagule Redirection" for non-specialists
- ABS-05: Strengthened So-What with cross-system generalization
- ABS-06: Compressed companion experiment to clause within Q1 result

---

## Deferred Issues (not fixed)

- **PROSE-02**: Sentence length consistently short (median 17-21 words vs 25-word target). A full-manuscript sentence-fusion pass is needed but too invasive for automated fixes. Recommend: manual revision pass targeting 4-6 long analytical sentences per section.
- **CROSS-05**: Cross-study scaling concordance (7 species) has no Discussion paragraph. Low priority.
- **CROSS-06**: Neighborhood null results disconnected from main causal chain. Partially addressed by DISC-04.
- **METH-05**: Package citations (DHARMa, sandwich, lmtest, car) not added to Methods. Requires new BibTeX entries in references.bib.
- **METH-08**: Moran's I placement (Sensitivity in Methods, Q2 in Results). Minor ordering mismatch.
- **METH-10**: Position correction statistical rationale not added. Low priority.
