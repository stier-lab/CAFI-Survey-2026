# CAFI Survey 2026 — Comprehensive Repo Audit
## Sunday March 8, 2026

Six parallel review agents covering: code quality, statistical methods, figures, manuscript text, data flow/pipeline, and documentation accuracy.

**De-duplicated totals across all agents:**

| Domain | Critical | High | Medium | Low | Total |
|--------|----------|------|--------|-----|-------|
| Code quality + bugs | 2 | 6 | 7 | 6 | 21 |
| Statistical methods | 0 | 1 | 7 | 6 | 14 |
| Figures | 0 | 1 | 6 | 4 | 11 |
| Manuscript text | 0 | 1 | 3 | 5 | 9 |
| Data flow + pipeline | 0 | 0 | 4 | 6 | 10 |
| Documentation (CLAUDE.md) | 0 | 2 | 4 | 3 | 9 |
| **Total** | **2** | **11** | **31** | **30** | **74** |

---

## 1. Code Quality + Bugs

### CRITICAL

| # | Script | Lines | Issue | Fix |
|---|--------|-------|-------|-----|
| CQ1 | `09_feedbacks.R` | 3021-3033 | `if_else()`/`scale()` mismatch in Part H: scalar `if()` inside dplyr `if_else()` with `scale()` returning full-length vector. Works by accident but fragile — will break if data shape changes. | Compute scaled values outside `mutate()` or use `as.numeric(scale(x))` pattern. |
| CQ2 | `05_scaling.R` | 1337 | `comm_pa[coral_id, sp]` matrix indexing uses coral_id as row names without validating they exist in `rownames(comm_pa)`. Silent NAs if mismatch. | Add `stopifnot(all(coral_data$coral_id %in% rownames(comm_pa)))`. |

### HIGH

| # | Script | Lines | Issue | Fix |
|---|--------|-------|-------|-----|
| CQ3 | `07_spatial.R` | 558-609 | Missing `set.seed()` before site-specific `moran.mc()` loop. P-values non-reproducible. | Add `set.seed(42)` before or inside loop. |
| CQ4 | `05_scaling.R` | 944-949 | `species_functional` may have duplicate rows per OTU if functional_group/type inconsistent. Join at L949 would inflate counts. | Add `stopifnot(!any(duplicated(species_functional$otu)))`. |
| CQ5 | `08_groups.R` | 214, 667 | Prevalence denominator uses `n_distinct(cafi_clean$coral_id)` — excludes zero-CAFI corals. Systematically inflates prevalence. | Replace with `nrow(coral_master)` or 114. |
| CQ6 | `08_groups.R` | 686 | `coral_data` variable overwritten in Part C (gastropods). Downstream Parts D+ see gastropod-specific data. | Use distinct name like `gastropod_data`. |
| CQ7 | `09_feedbacks.R` | 438-439 | Positional `pred_row` index from OLS applied to standardized model. Fragile if formula changes. | Use named coefficient extraction. |
| CQ8 | `09_feedbacks.R` | 2607-2610 | Same `if()`/`scale()` pattern as CQ1, different location. Works but fragile. | Same fix as CQ1. |

### MEDIUM

| # | Script | Lines | Issue | Fix |
|---|--------|-------|-------|-----|
| CQ9 | `05_scaling.R` | 1951 | `make_species_panel()` references `species_scaling_data` from closure, not as argument. | Pass as explicit parameter. |
| CQ10 | `05_scaling.R` | 1793, 1965 | `geom_smooth(method = MASS::glm.nb)` can crash on sparse species without `tryCatch`. | Wrap in error handler. |
| CQ11 | `09_feedbacks.R` | 3237-3300 | AIC comparison table mixes M1 (n~112) with M2-M5 (n~61). Table misleads even though code correctly restricts selection. | Add header note about non-comparable sample sizes. |
| CQ12 | `05_scaling.R` | ~178 | Bootstrap failed iteration count not logged. Silent CI degradation if >10% fail. | Log and warn if failure rate >10%. |
| CQ13 | `09_feedbacks.R` | 2480-2493 | Hardcoded experimental paper values in `tribble()`. Will be silently stale if companion paper revises. | Add source citation comment. |
| CQ14 | `01_load.R` | 332-397 | Cascading `drop_na()` + `complete.cases()` could silently drop many corals without warning. | Add warning if dropped count is large. |
| CQ15 | `01_load.R` | 446-451 | `community_matrix` row order not guaranteed after `rbind` of zero-CAFI rows. | Sort by coral_id after rbind. |

### LOW

| # | Script | Issue |
|---|--------|-------|
| CQ16 | `09_feedbacks.R` L3031 | `depth_scaled` computed but never used — dead code. |
| CQ17 | `05_scaling.R` L1514 | Variable named `s14_legend` but figure is S15 (renumbering artifact). |
| CQ18 | `06_cooccurrence.R` L744 | `pretty() %>% .[. == floor(.)]` may return empty vector for small ranges. |
| CQ19 | `08_groups.R` L462 | `functional_group` column assumed without existence check. |
| CQ20 | `05_scaling.R` various | Inconsistent `theme_publication()` vs `theme_minimal()` across supplement figures. |
| CQ21 | `08_groups.R` L323+ | In-plot legends on diagnostic figures (acceptable for non-manuscript). |

---

## 2. Statistical Methods

### HIGH

| # | Script | Issue | Fix |
|---|--------|-------|-----|
| ST1 | `02_community.R` | db-RDA uses Bray-Curtis on Hellinger-transformed data (double transformation). Non-standard; tested via sensitivity but underdocumented. | Document explicitly as "Bray-Curtis on Hellinger-transformed abundances." |

### MEDIUM

| # | Script | Lines | Issue | Fix |
|---|--------|-------|-------|-----|
| ST2 | `05_scaling.R` | 395-405 | Bootstrap for richness always fits Poisson even when main model used quasipoisson. CIs may be anti-conservative. | Match bootstrap family to main model. |
| ST3 | `09_feedbacks.R` | 327-329 | Power analysis denominator df: `v = n - 4` but model has 4 non-intercept params → correct is `v = n - 5`. | Use `v = n_physio - 5`. |
| ST4 | `09_feedbacks.R` | 748 | Abundance subsample model omits `log_volume` — inconsistent with all other models. | Add `log_volume` to formula. |
| ST5 | `09_feedbacks.R` | 2328-2335 | Species-trait correlations use Pearson on raw counts (zero-inflated, skewed). No FDR correction. | Use Spearman or sqrt-transform; add FDR. |
| ST6 | `06_cooccurrence.R` | 464-500 | Size-dependent co-occurrence uses different species pools per size class. SES not directly comparable. | Use common species pool or present as descriptive. |
| ST7 | `04_landscape.R` | 126-127 | Power analysis denominator df off: correct is `v = n - 7`, not `n - 5`. | Fix to `v = n_landscape - 7`. |

### LOW

| # | Script | Issue |
|---|--------|-------|
| ST8 | `05_scaling.R` L453 | Wald z-test for quasipoisson (should use t-test with residual df). |
| ST9 | `05_scaling.R` L1082 | NB GLM with only 5 non-zero obs — flag as exploratory. |
| ST10 | `02_community.R` L570 | Interaction PERMANOVA uses sequential Type I (order-dependent). |
| ST11 | `04_landscape.R` L660 | Functional group landscape models uncorrected for multiple testing. |
| ST12 | `06_cooccurrence.R` L360 | Intraspecific density null model ignores site structure. |
| ST13 | `13_sensitivity.R` L117 | Overdispersion threshold (2.0) differs from main analysis (1.5). |

---

## 3. Figures

### HIGH

| # | Figure | Issue | Fix |
|---|--------|-------|-----|
| FG1 | Fig 3 legend text | Factual error: "Of the 10 species shown: 11 Redirection, 10 Field of Dreams" — counts are from all 21 species, not the 10 shown. Self-contradictory. | Filter to shown species before counting, or change to "Across all N species tested." |

### MEDIUM

| # | Figure | Issue | Fix |
|---|--------|-------|-----|
| FG2 | Fig 2 Panel A | In-plot site legend violates no-legend policy. | Set `legend.position = "none"`. |
| FG3 | Fig 3 Panels A+B | In-plot legends (species names, group names). Panel A has 10 species — hard to direct-label. | Move to legend text file or bottom legend. |
| FG4 | Fig S4 (spatial) | HAU panel has overlapping longitude values on x-axis (~149.922-149.921). | Fewer decimal places or rotated labels. |
| FG5 | Fig S11 Panel A | Species names overlapping on x-axis of heatmap. | Increase figure width or reduce font size. |
| FG6 | Fig S16 Panel B | Appears nearly empty — only cells with FDR p < 0.10 shown, most are grey. | Add subtitle: "Only cells with FDR p < 0.10 shown." |

### LOW

| # | Figure | Issue |
|---|--------|-------|
| FG7 | Fig 5 | Panel tags at bottom-left instead of conventional top-left. |
| FG8 | Fig 4 Panel B | In-plot legend for stacked bar chart (arguably necessary). |
| FG9 | Supplement figs | 14/17 supplement figures lack companion legend text files. |
| FG10 | Fig S14 Panel C | Species names very small, may be clipped on y-axis. |

---

## 4. Manuscript Text

### HIGH

| # | File | Issue | Fix |
|---|------|-------|-----|
| MS1 | `results.md` L29 | Total CAFI abundance β for condition = 0.24 (wrong). Actual: **0.194** from `cafi_condition_models.csv`. | Change to β = 0.19. |

### MEDIUM

| # | File | Issue | Fix |
|---|------|-------|-----|
| MS2 | `results.md` L33 | Rarefied richness p = 0.45 (wrong). Actual: **p = 0.50** from output tables. | Change to p = 0.50. |
| MS3 | `results.md` L37, `methods.md` L65 | "10 key species tested" but `key_species_effects.csv` has **8 species/groups**. | Change to 8. |
| MS4 | `methods.md` L63, `supplement.md` L35 | Species-trait test count: methods says 100 (20 spp), supplement says 125 (25 spp). Actual: **95 tests (19 spp × 5 traits)**. | Reconcile all to 19 species, 95 tests. |

### LOW

| # | File | Issue |
|---|------|-------|
| MS5 | `results.md` L17 | Site PERMDISP (F=0.89, p=0.42) not in any output table — unverifiable. |
| MS6 | `results.md` L31 | Path model richness→condition β=0.55 but **p=0.20** (not significant). Not reported. |
| MS7 | Various | Neighborhood n=61 (manuscript) vs n=63 (CLAUDE.md) — minor but inconsistent. |
| MS8 | `supplement.md` L35 | FigS16 caption: grey threshold (0.10) vs significance threshold (0.05) — confusing. |
| MS9 | CLAUDE.md vs results.md | Richness-abundance r = 0.77 (CLAUDE.md) vs 0.84 (results.md, correct). CLAUDE.md needs update. |

---

## 5. Data Flow + Pipeline

### MEDIUM

| # | File | Issue | Fix |
|---|------|-------|-----|
| PL1 | `run_full_pipeline.R` L247 | `verify_dependencies` for script 08 requires `scaling_analysis_results` in global env, but `ensure_setup_and_data` only sources 01 (not 05). `run_one("08")` fails if 05 hasn't run. | Remove from verify_dependencies (script 08 has its own RDS check). |
| PL2 | `run_full_pipeline.R` L99-168 | EXPECTED_OUTPUTS missing 9/17 supplement figures. Undermines `skip_completed`. | Add missing figure paths. |
| PL3 | `run_full_pipeline.R` L238-261 | `verify_dependencies` only checks global env, not RDS files. Forces unnecessary re-sourcing. | Add `file.exists()` fallback for RDS files. |
| PL4 | `01_load_data.R` L853 | `tidyterra` package used but not checked with `requireNamespace`. Cryptic error if missing. | Add to package guard. |

### LOW

| # | File | Issue |
|---|------|-------|
| PL5 | `run_full_pipeline.R` L152-158 | Script 09 EXPECTED_OUTPUTS lists only 1 of 14 tables. |
| PL6 | `run_full_pipeline.R` L99 | `00b_data_quality_audit.R` has no EXPECTED_OUTPUTS entry. |
| PL7 | `run_full_pipeline.R` L101-111 | `overlap_genera.rds` saved but not in EXPECTED_OUTPUTS. |
| PL8 | `00_setup.R` L263 | `save_figure` PDF fallback could report success on corrupted partial writes. |
| PL9 | `01_load_data.R` L160 | Galeropsis count relies on string match — silently returns 0 if OTU name changes. |
| PL10 | `run_full_pipeline.R` L624 | `run_one("13")` docstring says "S13" but script 13 produces S14. |

---

## 6. Documentation (CLAUDE.md + README)

### HIGH

| # | File | Issue | Fix |
|---|------|-------|-----|
| DC1 | `run_full_pipeline.R` L92-94 | ML script paths reference `exploration/` but files are in `archive/`. `run_full_pipeline()` silently skips. | Update paths to `archive/`. |
| DC2 | `DATA_DICTIONARY.md` L148 | Says `log_volume = log₁₀(volume)` but code uses natural log. Same class of bug that caused wrong scaling exponents. | Fix to `log(volume)` (natural log). |

### MEDIUM

| # | File | Issue | Fix |
|---|------|-------|-----|
| DC3 | CLAUDE.md L346 | Output table count: says "~54 CSV". Actual: **66 CSV files**. | Update count. |
| DC4 | CLAUDE.md L361 | Output objects count: says "20 RDS". Actual: **22 RDS files**. | Update count. |
| DC5 | CLAUDE.md fig dirs | Figure subdirectory counts wrong: community 11→13, effects 14→12, scaling 7→9, feedbacks 8→13, network has ~9 legacy files. | Update counts or remove specifics. |
| DC6 | `DATA_DICTIONARY.md` L33-74 | Multiple wrong column names: coral_id format (`SITE_###` vs `SITE-POC##`), physiology columns, CAFI `count` column doesn't exist. | Comprehensive dictionary update. |

### LOW

| # | File | Issue |
|---|------|-------|
| DC7 | CLAUDE.md L82 | "PERMANOVA R² ~ 0.06-0.08" — runtime-dependent, may not match last run. |
| DC8 | CLAUDE.md L327-336 | Manuscript figures listed as `.png` only; PDFs also exist. |
| DC9 | `DATA_DICTIONARY.md` L239-277 | Describes legacy network outputs that no longer exist. |

---

## Top 10 Priorities for Fixing

1. **CQ1** (CRITICAL): `if_else`/`scale()` mismatch in script 09 — fragile, will break
2. **CQ2** (CRITICAL): Unvalidated matrix indexing in script 05 — silent wrong results
3. **MS1** (HIGH): Wrong β value in results.md (0.24 → 0.19)
4. **FG1** (HIGH): Fig 3 legend text self-contradictory species counts
5. **CQ5** (HIGH): Prevalence denominator excludes zero-CAFI corals in script 08
6. **DC2** (HIGH): DATA_DICTIONARY says log10 but code uses ln
7. **ST3+ST7** (MEDIUM): Power analysis denominator df off in scripts 04 and 09
8. **MS4** (MEDIUM): Species-trait test counts inconsistent across manuscript files
9. **ST4** (MEDIUM): Abundance subsample model missing log_volume
10. **PL1** (MEDIUM): `run_one("08")` broken due to dependency check design
