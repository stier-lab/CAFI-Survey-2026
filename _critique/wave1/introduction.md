## Introduction Diagnosis

### Score: 33/40

### Rubric

| Dimension | Score (1-5) | Notes |
|-----------|-------------|-------|
| Opening hook | 4 | Strong intellectual opening on functional form of habitat-occupant scaling. First sentence is substantive. Minor issue: the combined_manuscript.md version has a weaker, more generic opener ("In many ecosystems, habitat-forming foundation species...") — the introduction.md standalone version is sharper. The two versions diverge; the combined version should use the standalone opener. |
| Literature curation | 4 | Well-curated argument structure. Each citation does work. One over-citation cluster in ¶4 (Beck et al. 2011; Byrnes et al. 2011; Kimbro et al. 2017; Rennick et al. 2022) lists four references for a single motivational claim about habitat loss; any two would suffice. The Abele/Gotelli/Counsell/Curtis citation cluster in ¶6 also runs long. |
| Gap clarity | 4 | The gap is stated explicitly in ¶7: "Three links in the habitat--community--condition chain remain untested in established communities." This is clear and sharp. Minor problem: the gap paragraph jumps immediately from gap to methods ("Here, we survey...") without a sentence framing WHY this gap matters now, making the transition feel abrupt. |
| Hypothesis specificity | 4 | Q1-Q3 have clear directional predictions. Q1 predicts sublinear scaling and passive-sampling as the richness mechanism — named hypotheses used well (Field of Dreams, Propagule Redirection). Q2 predicts turnover not nestedness. Q3 predicts a diversity-condition relationship beyond abundance. The named hypothesis apparatus is strong. Minor weakness: Q3 prediction is framed in the combined version with a hedge ("alternatively, any apparent diversity signal may reflect...") that softens the prediction in the opening before the reader even enters the Results. |
| Scope signaling | 5 | The closing paragraph ("Together, these three questions...") explicitly signals scope, situates within companion experiment, and lands the generalization to oyster reefs and forest canopies. Well done. |
| Paragraph logic | 3 | The funnel is largely sound but ¶4 (consequences of scaling regimes) and ¶5 (*Pocillopora* system) are doing double work: ¶4 ends with habitat loss as motivation, then ¶5 also opens with habitat loss before pivoting to the system description. The transition is repetitive. More critically, ¶6 (prior CAFI scaling literature) and ¶7 (gap + objectives) are partially fused in the combined version, creating a run-on paragraph that encompasses the entire evidence-to-gap-to-objectives transition — a structural load-bearing section that deserves cleaner separation. Within ¶5, the "species stacking" mechanism is introduced and explained at length (4 sentences), then repeated conceptually in the Q2 prediction. |
| Citation quality | 4 | Primary literature throughout. Named citations for key empirical examples (Resetarits & Binckley 2009, Stier & Osenberg 2010, Keller et al. 2017). One placeholder: Brush et al. 2026 — needs updating before submission. One flagged citation: "Abele 1976" vs. "Gotelli & Abele 1985" appears in both combined_manuscript.md and results — check that the correct reference is used throughout (Gotelli & Abele 1985 is the species-level exponent paper; Abele & Patton 1976 is the community-level estimate cited in Discussion). The Discussion cites "Abele and Patton (1976)" while the Introduction cites "Gotelli & Abele 1985" — these are different papers; ensure the correct one is cited for each specific claim. |
| Length proportion | 4 | The introduction runs 7 paragraphs (standalone version) or the equivalent in the combined file. Word count is appropriate for a JAE manuscript (~650 words in standalone; the combined version is longer at ~900 words due to expanded ¶2-¶4). The combined version is at the outer edge of the ≤20% rule. Standalone version is well-proportioned. |

### PASS/FAIL: PASS (33/40; no dimension below 3)

---

### Issues Found

**1. [INTRO-01]** severity: critical
- Problem: The combined_manuscript.md Introduction (the authoritative submission version) opens with a generic foundation-species survey sentence: "In many ecosystems, habitat-forming foundation species -- trees, corals, kelps, oysters -- define the structure and composition of associated animal communities through their size, density, and spatial arrangement." This violates the anti-pattern rule ("Never open with 'Coral reefs are important ecosystems...'") and the Stier voice rule ("Start with the intellectual problem"). The standalone introduction.md opens with the stronger, intellectually-leading sentence: "The relationship between habitat quantity and occupant abundance is among the most general patterns in ecology... yet its functional form... remains poorly resolved in most systems." The combined version should use the standalone opener.
- Location: combined_manuscript.md line 35, first sentence of Introduction section
- Proposed fix: Replace the combined_manuscript.md opening sentence with the standalone version's opener: "The relationship between habitat quantity and occupant abundance is among the most general patterns in ecology (Connor & McCoy 1979; Preston 1962; Rosenzweig 1995), yet its functional form -- whether occupant abundance increases proportionally with habitat or less than proportionally -- remains poorly resolved in most systems. In ecosystems structured by habitat-forming foundation species -- trees, corals, kelps, oysters -- this distinction carries consequences well beyond the occupants themselves (Dayton 1972; Hanski 1998; MacArthur & Wilson 1967)."
- Ripple risk: Abstract (which references the "chain from habitat size through community assembly to habitat condition" framing already set up correctly) — no ripple. The combined_manuscript.md Introduction section needs to be resynchronized with the standalone introduction.md.

**2. [INTRO-02]** severity: major
- Problem: Paragraphs 4 and 5 (combined version) both open with habitat-loss motivation before pivoting. ¶4 ends: "particularly as ongoing habitat loss restructures the occupant communities that influence recovery." ¶5 opens: "The ongoing decline of biogenic habitats in marine systems makes this scaling question urgent, because losses of foundation species simultaneously restructure the occupant communities that influence habitat recovery." This is near-verbatim repetition of the same motivational claim across consecutive paragraphs, diluting both. The habitat-loss motivation belongs in one place — either as ¶4's closing sentence or as ¶5's opener, not both.
- Location: Introduction ¶4 closing sentence and ¶5 opening sentence (combined_manuscript.md lines 39-41)
- Proposed fix: Cut the habitat-loss sentence from ¶4's close ("The distinction matters for predicting ecosystem trajectories... remains limited to a few systems."). Consolidate the motivation into ¶5's opener. Rewrite ¶4's closing as a conceptual bridge to ¶5 by ending on the self-reinforcing-spatial-patterns claim ("These feedbacks can generate self-reinforcing spatial patterns...") and removing the duplicated conservation motivation.
- Ripple risk: None — purely within Introduction.

**3. [INTRO-03]** severity: major
- Problem: The species-stacking mechanism is introduced and fully explained in ¶5 ("This process -- termed 'species stacking' -- generates per-colony density ceilings..."), then restated in the Q2 prediction: "If species stacking drives assembly, small-coral faunas should not be nested subsets of large-coral faunas." The mechanism explanation in ¶5 already includes the example ("because obligate symbionts such as *Trapezia* crabs and *Alpheus* shrimp occupy corals as territorial pairs, adding abundance beyond one pair requires adding a new species"). This creates redundancy between the system paragraph and the objectives paragraph. The objectives paragraph should name the prediction without re-explaining the mechanism.
- Location: ¶5 system description (combined_manuscript.md line 41, last 4 sentences) and Q2 prediction (combined_manuscript.md line 47)
- Proposed fix: Shorten ¶5's species-stacking explanation to 2 sentences (the named concept + one example). The Q2 prediction already handles the consequence ("If species stacking drives assembly, small-coral faunas should not be nested subsets...") — the reader has enough context from ¶5 and does not need the mechanism restated. Trim: "For example, because obligate symbionts such as *Trapezia* crabs and *Alpheus* shrimp occupy corals as territorial pairs, adding abundance beyond one pair requires adding a new species -- so that larger corals accumulate diversity through species replacement rather than density increases within species." can be condensed to one clause within the prior sentence.
- Ripple risk: None.

**4. [INTRO-04]** severity: major
- Problem: ¶6 (combined version) is structurally overloaded. It opens with prior CAFI scaling evidence, embeds two specific quantitative examples (Gotelli & Abele 1985 exponent range; Counsell 2018 R² claim), pivots to the companion experiment, then transitions directly into the gap statement and the study's methods — all in one paragraph. This compresses three logical beats (prior evidence → gap → what we did) into a single paragraph, obscuring the gap statement that should land with impact.
- Location: Combined_manuscript.md lines 43-44 (the long paragraph from "CAFI abundance increases..." through "...Here, we present a survey")
- Proposed fix: Split into two paragraphs. First paragraph: prior evidence (Gotelli, Counsell, Curtis) + companion experiment. Second paragraph: gap statement opening with "Yet three links...remain untested" → "Here, we survey..." This separation gives the gap statement its own paragraph and the impact it deserves. The standalone introduction.md already implements this split — the combined version should match.
- Ripple risk: None — structural only.

**5. [INTRO-05]** severity: minor
- Problem: Q3 prediction in the combined version ends with an alternative framing: "alternatively, any apparent diversity signal may reflect the strong richness--abundance correlation." This is technically accurate but is an anti-prediction — it presents the null result as a co-equal hypothesis alongside the directional prediction, which softens the hypothesis before the reader reaches Results. Introductions should state the directional prediction; the null alternative belongs in the Discussion.
- Location: Combined_manuscript.md line 49, Q3 prediction last clause
- Proposed fix: Remove "alternatively, any apparent diversity signal may reflect the strong richness--abundance correlation." from the Q3 prediction. The Q3 prediction should end at "BEF theory predicts that complementarity among species enhances ecosystem performance (Hooper et al. 2005; Loreau & Hector 2001)." The rarefied richness null and the interpretation ambiguity are already handled in the Results and Discussion.
- Ripple risk: Discussion (¶3 of Composition section already addresses the ambiguity correctly) — no change needed there.

**6. [INTRO-06]** severity: minor
- Problem: Citation inconsistency between Introduction and Discussion for the foundational Abele paper. The Introduction (combined version, ¶6) cites "Abele 1976" as first in the citation cluster. The Discussion cites "Abele and Patton (1976)" with the specific β = 0.62 estimate. These appear to be the same paper, but "Abele 1976" and "Abele and Patton 1976" are different citation formats — it may be two different papers (Abele 1976 alone vs. Abele & Patton 1976). The Results cites "Gotelli & Abele 1985" for the 0.06–0.64 exponent range. Verify that Abele (1976), Abele & Patton (1976), and Gotelli & Abele (1985) are all distinct references correctly assigned to distinct claims, and that the BibTeX entries are consistent.
- Location: Introduction ¶6 (combined version); Discussion ¶1 of "Sublinear scaling" section
- Proposed fix: Check references.bib. If "Abele 1976" is the same paper as "Abele & Patton 1976," consolidate to one key. Verify that each citation in the Introduction cluster (@Abele1976; @Gotelli1985; @Counsell2018; @Curtis2023) is assigned to the correct specific claim.
- Ripple risk: References section — BibTeX cleanup only.

**7. [INTRO-07]** severity: minor
- Problem: "Brush et al. 2026" appears in ¶6 as a citation supporting cross-system generality of sublinear CAFI scaling. This is a future-dated reference, presumably a paper in press or in review. Before submission to JAE, this either needs a final publication date or the citation should note "(in press)" or "(in review)." JAE does not allow unreferenced future citations.
- Location: Introduction ¶6 / combined_manuscript.md line 43
- Proposed fix: Update to "Brush et al. (in press)" or supply the 2026 volume/page numbers once available. If not yet accepted, note "(in review)" and verify JAE's policy on citing in-review manuscripts.
- Ripple risk: References.bib — update BibTeX entry.

**8. [INTRO-08]** severity: minor
- Problem: The transition from ¶3 to ¶4 (both versions) uses "These scaling regimes carry distinct consequences for habitat-forming organisms" (standalone) or "These scaling regimes carry distinct consequences when the habitat itself is biogenic and responsive to occupant density" (combined). The combined version's transition is more precise and should be retained, but it lacks a conceptual bridge from ¶3's key mechanistic point (that both regimes can operate simultaneously) to ¶4's consequences. The standalone version transitions more cleanly; the combined version opens ¶4 by restating what was just established rather than advancing the argument.
- Location: Combined_manuscript.md line 39, ¶4 opening sentence
- Proposed fix: Rewrite the bridge: "The regime governing a given system carries distinct consequences for the habitat-forming organisms themselves." This presupposes the regime distinction established in ¶3 rather than re-introducing it.
- Ripple risk: None.

---

### Diagnostic Checklist

- [x] Gap statement is findable and explicit ("Three links in the habitat--community--condition chain remain untested in established communities")
- [x] Every middle paragraph advances toward the gap (broadly yes, with noted redundancy in ¶4-¶5 transition)
- [ ] No paragraph is a literature tangent (¶4 habitat-loss close partially repeats ¶5 motivation opener — see INTRO-02)
- [x] Hypotheses stated with directional predictions (Q1: sublinear; Q2: turnover not nestedness; Q3: diversity → condition)
- [ ] No alternative hypothesis pre-stated in the Introduction (Q3 in combined version includes the null alternative — see INTRO-05)
- [x] A neighboring-field reader (community ecologist) can follow the logic without Googling
- [x] Opening sentence is substantive, not definitional (standalone version; combined version fails this — see INTRO-01)
- [x] Named hypotheses used as cognitive hooks (Field of Dreams, Propagule Redirection, species stacking)
- [x] Cross-system analogies present (oyster reefs, forest canopies, kelp holdfasts)
- [x] Introduction ends with clear objectives and predictions
- [ ] No redundancy between system paragraph and objectives paragraph (species-stacking mechanism restated — see INTRO-03)
- [x] Gap paragraph does not rush methods (the study sentence follows gap cleanly)
- [ ] ¶6 (prior evidence → gap → methods) is split correctly (combined version fuses three beats; see INTRO-04)
- [x] Length ≤20% of manuscript (borderline in combined version; within range)
- [x] No results previewed
- [x] Closing sentence generalizes beyond the focal system
- [ ] Citation consistency verified (Abele 1976 vs Abele & Patton 1976 — see INTRO-06)
- [ ] Placeholder citations resolved (Brush et al. 2026 — see INTRO-07)
- [x] Active voice dominant throughout
- [x] No mechanical connectors ("Additionally," "Furthermore") as paragraph openers
- [x] Paragraph topic sentences make claims, not announcements

---

### Summary for Wave 2 Reconciliation

The Introduction is structurally sound and passes the rubric (33/40). The primary issue is a **version mismatch**: the standalone `introduction.md` is tighter and better structured than the corresponding section in `combined_manuscript.md`. The combined version weakens the opening hook (INTRO-01), fuses the prior-evidence/gap/methods paragraph (INTRO-04), and repeats the habitat-loss motivation across consecutive paragraphs (INTRO-02). These are all structural edits, not content changes. The named hypothesis apparatus (Field of Dreams, Propagule Redirection, species stacking) is a genuine strength and should be preserved exactly. One content fix is needed across versions: remove the alternative hypothesis from the Q3 prediction (INTRO-05). Two citation issues need resolution before submission (INTRO-06, INTRO-07).

**Cross-section signals for Wave 2:**
- Results ¶1 (Q1) correctly mirrors Introduction Q1 prediction — parallel structure intact.
- Results ¶3 (Q2) correctly references nestedness prediction from Introduction — parallel structure intact.
- Discussion ¶3 (BEF) correctly handles the rarefied richness ambiguity that should NOT appear in Introduction Q3 — the fix in INTRO-05 requires no Discussion change.
- The "Abele" citation discrepancy (INTRO-06) may also affect the Discussion's quantitative comparison and the Methods' regression to prior estimates.
