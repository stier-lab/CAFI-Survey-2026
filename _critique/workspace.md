# Manuscript Critique Workspace
Generated: 2026-03-18
Manuscript: manuscript/ (split section files)

---

## Wave 1 Summary

| Section | Score | Status | Critical/Major Issues |
|---------|-------|--------|----------------------|
| Introduction | 33/40 | PASS | combined_manuscript.md out of sync with introduction.md |
| Methods | 33/40 | FAIL (Dim 8 = 1) | Placeholder DOIs; n=47 drop unexplained; package citations missing |
| Results | 29/40 | FAIL (Dim 7,8 = 2) | Interpretation creep (3 passages); methods recapitulation (3 passages); Q3 buries independence-from-abundance |
| Discussion | 30/40 | FAIL (Dim 4 = 2) | Personal observation as evidence; derived quantities not in Results; limitations lack rival causal stories; neighborhood null missing |
| Abstract | 33/40 | PASS | Opening leads with mechanism; gap buried; Q2 lacks numbers |
| Cross-Section | 33/40 | PASS | "Species stacking" absent from Methods/Results; Q2+Q3 merged in Discussion; terminology drift |
| Prose Style | 31/40 | FAIL (Dim 4 = 2) | "we" absent from Results; sentence length short; "consistent with" ×14; "rather than" ×12 |

Full diagnoses: `_critique/wave1/*.md`

---

## Wave 2: Reconciled Edit Plan

### Conflict Resolution

1. **CONFLICT-01: Introduction version mismatch**
   - The introduction.md (standalone) was edited earlier in this session and is tighter than the combined_manuscript.md version. All section fixes should be applied to the standalone section files. combined_manuscript.md will be synced at the end.

2. **CONFLICT-02: Interpretation in Results vs. Stier voice**
   - Results agent flags "interpretation creep" (RES-01). Prose agent notes Stier voice "integrates results with interpretation." Resolution: Brief orienting statements ("consistent with Propagule Redirection") are OK. Remove full mechanistic explanations and inferential conclusions. The three flagged passages (passive sampling interpretation, stochasticity claim, path-model closing) all cross the line into Discussion territory.

3. **CONFLICT-03: Derived quantities in Discussion (DISC-01)**
   - Discussion agent flags 2^0.52^ ≈ 1.43 and 2^-0.48^ ≈ 72% as "new results." These are algebraic derivations from reported β, not new analyses. Resolution: Flag as derivations in Discussion ("The sublinear exponent means that...") rather than adding to Results. This is interpretive arithmetic, not new data.

### Ripple Map

| Source Change | Affected Sections | Required Update |
|--------------|-------------------|-----------------|
| RES-01: Remove interpretation from Results | Discussion | Already has content — just cut from Results |
| RES-03: Restructure Q3 opening | Discussion | No change needed |
| CROSS-01: Add "species stacking" to Results Q2 | Methods Q2, Results Q2 | Add term at nestedness result |
| CROSS-02: BEF label consistency | Results Q3 header | Add "(BEF framework)" |
| CROSS-07: "causal chain" vs "feedback loop" | Intro, Abstract, Discussion | Use "causal chain" in Intro; "feedback" in Discussion after reverse-direction tested |
| CROSS-08: Hochberg correction visibility | Results Q3 | Add correction note to p-values |
| PROSE-01: Add "we" to Results | Results | Systematic pass |
| DISC-02: Remove personal observation | Discussion heatwave para | Soften to published literature only |
| METH-01: Placeholder DOIs | Abstract, combined_manuscript | User action before submission |
| ABS-03: Add Q2 number to Abstract | Abstract | After Results finalized |

### Edit Order (dependency-aware)

1. **Results fixes** (RES-01 through RES-12, CROSS-01, CROSS-08, PROSE-01, PROSE-06) — highest-impact section, most issues
2. **Discussion fixes** (DISC-01 through DISC-08, CROSS-04, PROSE-03/04/05) — depends on knowing what was removed from Results
3. **Methods fixes** (METH-02 through METH-10) — independent of Results/Discussion changes
4. **Abstract rewrite** (ABS-01 through ABS-06) — depends on all above
5. **Prose style final pass** — all sections, after structural fixes
6. **Sync combined_manuscript.md** — final step

### Cross-Section Issues (new from combining diagnoses)

1. **CROSS-NEW-01**: combined_manuscript.md must be synced with all section files after Wave 3. The introduction.md has already been edited; combined_manuscript.md has the old version.
2. **CROSS-NEW-02**: Fig. S10A/S10C ARE cited in Results Q3 ¶2-¶3 before Discussion cites them — confirmed not orphaned. The cross-section agent's concern was precautionary.
3. **CROSS-NEW-03**: The "species stacking" gap (CROSS-01) and the "BEF label" gap (CROSS-02) are the same underlying issue: the conceptual vocabulary established in the Introduction doesn't carry through Methods/Results. Fixing both in one pass is most efficient.
