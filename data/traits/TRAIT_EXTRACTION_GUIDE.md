# CAFI Traits Database - Complete Extraction Guide

**Generated**: 2025-11-02
**Status**: ✅ OBIS extraction complete | 📦 R script ready for execution

---

## Overview

This document explains the complete trait extraction pipeline for the 233 CAFI taxa, including the logic behind each step, data sources, and how to use the provided scripts.

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    CAFI TAXA LIST (233)                     │
│         survey_cafi_data_w_taxonomy_summer2019_v5.csv       │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│                   STEP 1: TAXONOMY                          │
│                   build_cafi_traits_database.py              │
│                   Source: WoRMS REST API                     │
│                   Output: cafi_traits.jsonl                  │
│                   Coverage: 99.6% (232/233)                  │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│              STEP 2: ASSOCIATION TYPES                       │
│              update_traits_with_associations.py              │
│              Source: Lab reference sheet (2012)              │
│              Output: cafi_traits_updated.jsonl               │
│              Coverage: 34.8% (81/233)                        │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│               STEP 3: DEPTH RANGES (OBIS)                    │
│               add_obis_depth_traits.py                       │
│               Source: OBIS occurrence API                    │
│               Output: cafi_traits_with_obis.jsonl            │
│               Coverage: 74.2% (173/233)                      │
│               ✅ COMPLETE                                    │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│          STEP 4: LIFE HISTORY TRAITS (FishBase/SLB)         │
│          extract_traits_rfishbase.R                          │
│          Source: FishBase + SeaLifeBase                      │
│          Output: cafi_traits_fishbase.csv                    │
│          Coverage: TBD (requires R execution)                │
│          📦 READY TO RUN                                     │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│                FINAL INTEGRATED DATABASE                     │
│                Used in network analysis (Script 05, 06)      │
└─────────────────────────────────────────────────────────────┘
```

---

## Data Sources and Logic

### 1. WoRMS (World Register of Marine Species)
**What it provides**: Taxonomic backbone (accepted names, AphiaIDs, classification)

**Logic**:
- Query WoRMS REST API for each taxon name
- Match to accepted species name (resolves synonyms)
- Obtain AphiaID (unique taxonomic identifier)
- Extract full classification (Phylum → Species)

**Coverage**: 99.6% (232/233 taxa resolved)

**Confidence**: HIGH - Authoritative taxonomic database, expert-curated

---

### 2. OBIS (Ocean Biodiversity Information System)
**What it provides**: Depth ranges from occurrence records

**Logic**:
1. OBIS aggregates species occurrence records from global marine surveys
2. Each occurrence may include depth measurements (min/max depth in meters)
3. For each species, we query all occurrence records
4. Extract:
   - `depth_min_m`: Shallowest observation (min of all minimumDepthInMeters)
   - `depth_max_m`: Deepest observation (max of all maximumDepthInMeters)
   - `n_obis_records`: Number of records with depth data (confidence indicator)

**Interpretation**:
- Depth ranges represent **observed depths**, not physiological limits
- Species with n_records < 5 should be considered low confidence
- Some outliers possible (measurement errors in source data)
- Absence of data means no OBIS records with depth, NOT that species has no depth preference

**Coverage**: 74.2% (173/233 taxa with depth data)

**Confidence**:
- HIGH: ≥20 records
- MEDIUM: 5-19 records
- LOW: 1-4 records

**API details**:
- Base URL: `https://api.obis.org/occurrence`
- No API key required
- Rate limit: 0.5s delay between requests (conservative)

**Example depth data quality**:
```
Caracanthus maculatus: 0-20.7 m (n=6 records) → MEDIUM confidence
Alpheus diadema: 3-5 m (n=1 records) → LOW confidence
Trapezia serenei: 12-None m (n=1 records) → LOW, incomplete range
Paragobiodon modestus: No data → Species exists but OBIS has no depth records
```

---

### 3. FishBase / SeaLifeBase (via rfishbase)
**What it provides**: Comprehensive life history traits

**Logic**:
1. **FishBase** covers fish species; **SeaLifeBase** covers invertebrates
2. Query multiple trait tables:
   - `species()`: Max length, common length, weight
   - `ecology()`: Depth range (min/max), trophic level, diet
   - `reproduction()`: Spawning mode, fertilization, fecundity
3. Merge tables by SpecCode (unique identifier)
4. **Reconcile depth data** between FishBase/SeaLifeBase and OBIS:
   - If FishBase/SeaLifeBase has depth: USE IT (more authoritative, based on expert review)
   - If missing: Use OBIS depth (occurrence-based)
   - If both available: Prefer FishBase/SeaLifeBase but flag discrepancies

**Why FishBase/SeaLifeBase depth is better than OBIS**:
- Expert-reviewed and curated (not raw occurrence data)
- Represents ecological depth range, not just observation bias
- Includes common depth range (typical habitat depth)
- References primary literature

**Coverage**: TBD (requires R execution)

**Expected coverage** (based on database size):
- Fish species: ~30-40% (FishBase has good coverage for Indo-Pacific reef fish)
- Invertebrates: ~20-30% (SeaLifeBase less comprehensive than FishBase)
- Cryptic CAFI: Lower coverage (many rare/undescribed species)

**Confidence**: HIGH - Peer-reviewed, species-specific data

**Requires**: R + rfishbase package

---

## Script Usage

### Script 1: OBIS Depth Extraction (Python) ✅ COMPLETE

**File**: `scripts/add_obis_depth_traits.py`

**Status**: ✅ Executed successfully (2025-11-02)

**Results**:
- Total processed: 233/233 taxa
- With depth data: 173/233 (74.2%)
- Without depth data: 60/233 (25.8%)
- Processing time: 312.6 seconds (~5 minutes)
- Output: `output/Survey/traits/cafi_traits_with_obis.jsonl`
- Output: `output/Survey/traits/cafi_traits_with_obis.csv`

**To re-run** (not necessary unless updating data):
```bash
python3 scripts/add_obis_depth_traits.py
```

**Test mode** (5 species sample):
```python
# Edit line 48 in add_obis_depth_traits.py:
TEST_MODE = True
```

---

### Script 2: FishBase/SeaLifeBase Extraction (R) 📦 READY TO RUN

**File**: `scripts/extract_traits_rfishbase.R`

**Status**: 📦 Script created, ready for execution

**Prerequisites**:
1. Install R (if not already installed): https://www.r-project.org/
2. Install required packages:
```r
install.packages("rfishbase")
install.packages("tidyverse")
```

**To run**:
```bash
Rscript scripts/extract_traits_rfishbase.R
```

**Expected runtime**: 5-10 minutes (depends on internet speed)

**Outputs**:
- `output/Survey/traits/cafi_traits_fishbase.csv` - Final trait database with FishBase/SeaLifeBase data
- Console log with coverage statistics

**What the script does**:
1. Loads CAFI taxa from `cafi_traits_with_obis.csv`
2. Queries FishBase for fish species
3. Queries SeaLifeBase for invertebrate species
4. Merges trait tables (species, ecology, reproduction)
5. Reconciles depth data (FishBase/SeaLifeBase > OBIS)
6. Exports final integrated trait database
7. Prints coverage summary

**Troubleshooting**:
- If `rfishbase` installation fails: Update R to version ≥4.0
- If species not found: Normal - many CAFI are cryptic/rare species not in database
- If API timeout: Increase timeout in rfishbase (see package documentation)

---

## Depth Data Reconciliation Logic

When both OBIS and FishBase/SeaLifeBase have depth data, we prioritize FishBase/SeaLifeBase. Here's why:

| Aspect | OBIS | FishBase/SeaLifeBase |
|--------|------|---------------------|
| **Data source** | Raw occurrence records | Expert-reviewed literature |
| **Depth type** | Observation depth (where sampled) | Ecological depth range (where species lives) |
| **Bias** | Sampling bias (e.g., shallow sites overrepresented) | Ecological synthesis across studies |
| **Outliers** | Possible measurement errors | Curated, outliers removed |
| **Confidence** | Depends on n_records | HIGH (peer-reviewed) |
| **Example** | "Found at 3m, 5m, 200m" → range 3-200m (outlier?) | "Lives 0-30m, common 5-15m" (ecological) |

**Reconciliation algorithm** (in `extract_traits_rfishbase.R`):
```r
depth_min_final = case_when(
  !is.na(DepthRangeShallow) ~ DepthRangeShallow,  # 1. Prefer FishBase/SeaLifeBase
  !is.na(depth_min_m) ~ depth_min_m,              # 2. Fall back to OBIS
  TRUE ~ NA_real_                                 # 3. Missing
)
```

**Flag large discrepancies** (optional analysis):
```r
# After running R script, check for discrepancies:
cafi_traits_fishbase %>%
  filter(!is.na(DepthRangeDeep) & !is.na(depth_max_m)) %>%
  mutate(depth_diff = abs(DepthRangeDeep - depth_max_m)) %>%
  filter(depth_diff > 50) %>%  # More than 50m difference
  select(accepted_name, DepthRangeDeep, depth_max_m, depth_diff)
```

---

## Current Coverage Summary

| Trait | Source | Coverage | Status |
|-------|--------|----------|--------|
| **Taxonomy** | WoRMS | 99.6% (232/233) | ✅ Complete |
| **Pocillopora association** | Lab reference | 34.8% (81/233) | ✅ Complete |
| **Functional role** | Literature | 24.0% (56/233) | ✅ Complete |
| **Depth range** | OBIS | 74.2% (173/233) | ✅ Complete |
| **Max length** | FishBase/SeaLifeBase | TBD | 📦 Ready to run |
| **Trophic level** | FishBase/SeaLifeBase | TBD | 📦 Ready to run |
| **Reproduction** | FishBase/SeaLifeBase | TBD | 📦 Ready to run |

---

## Expected Final Coverage (after R script)

Based on database sizes and CAFI taxonomic composition:

| Trait | Expected Coverage | Confidence |
|-------|------------------|------------|
| Taxonomy | 99.6% | HIGH |
| Pocillopora association | 34.8% | HIGH (evidence-based) |
| Functional role | 40-50% | MEDIUM (after rfishbase) |
| Depth range | 80-85% | HIGH (OBIS + FishBase/SeaLifeBase combined) |
| Max length | 35-45% | HIGH (for matched species) |
| Trophic level | 25-35% | HIGH (for matched species) |
| Reproduction | 15-25% | HIGH (for matched species) |

**Species with complete trait coverage** (all traits): ~15-20% (35-45 taxa)

**Species with partial trait coverage** (≥3 traits): ~60-70% (140-165 taxa)

**Species with minimal trait coverage** (taxonomy only): ~15-20% (35-45 taxa)

---

## Missing Data Strategy

### High Priority (Network Keystone Species)
Target manual extraction for these 5 keystone species:
1. **Alpheus diadema** (keystone, degree=12)
2. **Caracanthus maculatus** (keystone predator, degree=7)
3. **Alpheus collumianus** (keystone, degree=8)
4. **Breviturma pica** (keystone, degree=6)
5. **Periclimenes** sp. (keystone, degree=7)

**Method**: Targeted literature review + FishBase/SeaLifeBase manual lookup

---

### Medium Priority (Abundant Obligates)
These species are ecologically important (obligate, abundant):
1. **Trapezia serenei** (n=452 individuals)
2. **Paragobiodon modestus** (n=69)
3. All other Trapezia spp. (obligate mutualists)

**Method**: Already have good functional role data; add size/depth if missing

---

### Low Priority (Rare Taxa)
~70 taxa with <5 individuals observed

**Method**:
- Accept trait gaps for now
- Mark as "low priority" in database
- Revisit if species becomes abundant in future surveys

---

## Data Quality Flags

Each trait should be flagged with confidence level:

| Confidence | Criteria | Example |
|-----------|----------|---------|
| **HIGH** | Species-specific, peer-reviewed data | FishBase depth for *Caracanthus maculatus* |
| **MEDIUM** | Genus-level inference, or OBIS with ≥5 records | Alpheidae functional role, OBIS depth n=8 |
| **LOW** | Family-level inference, or OBIS with <5 records | Unknown functional role, OBIS depth n=1 |
| **UNKNOWN** | No data available | No depth data for *Paragobiodon modestus* |

**In trait database**:
```json
"confidence": {
  "per_trait": {
    "taxonomy": "high",
    "functional_role": "medium",
    "depth_range": "high",
    "max_length": "high",
    "trophic_level": "unknown"
  }
}
```

---

## Integration with Network Analysis

Once traits are extracted, use in network analyses:

### Test 1: Trait Effects on Network Position
```r
# Do larger species have higher degree?
lm(degree ~ max_length_final_mm, data = network_traits)

# Do deeper species have different module membership?
glm(module ~ depth_max_final, family = binomial, data = network_traits)
```

### Test 2: Association Type and Co-occurrence
```r
# Do obligate species co-occur more than facultative?
lm(cooccurrence ~ pocillopora_association, data = network_traits)
```

### Test 3: Functional Role and Network Centrality
```r
# Are mutualists more central than commensals?
aov(betweenness ~ functional_role_str, data = network_traits)
```

---

## Next Steps

### Immediate (this session):
1. ✅ OBIS depth extraction (COMPLETE)
2. ✅ Create R script for rfishbase (COMPLETE)
3. ⏳ Document pipeline (IN PROGRESS)

### Follow-up (requires user action):
1. **Run R script**: `Rscript scripts/extract_traits_rfishbase.R`
2. **Review coverage**: Check `cafi_traits_fishbase.csv` for completeness
3. **Integrate with network analysis**: Merge traits with Script 05, 06 outputs
4. **Manual extraction**: Target keystone species with missing traits

### Long-term:
1. Add trait data to network visualizations (color by functional role, size by max length)
2. Test trait-based hypotheses about network structure
3. Update trait database annually (as new literature published)

---

## File Reference

| File | Purpose | Status |
|------|---------|--------|
| `build_cafi_traits_database.py` | Initial taxonomy + functional roles | ✅ Used |
| `update_traits_with_associations.py` | Add reference sheet associations | ✅ Used |
| `add_obis_depth_traits.py` | Extract OBIS depth ranges | ✅ Complete |
| `extract_traits_rfishbase.R` | Extract FishBase/SeaLifeBase traits | 📦 Ready |
| `cafi_traits.jsonl` | Taxonomy + functional roles | ✅ Generated |
| `cafi_traits_updated.jsonl` | + Reference sheet associations | ✅ Generated |
| `cafi_traits_with_obis.jsonl` | + OBIS depth data | ✅ Generated |
| `cafi_traits_with_obis.csv` | CSV version for analysis | ✅ Generated |
| `cafi_traits_fishbase.csv` | Final database with all traits | 📦 Pending R execution |
| `association_lookup.csv` | Reference sheet (parsed) | ✅ Generated |
| `trait_database_summary.md` | Summary documentation | ✅ Generated |
| `api_databases_reference.md` | API documentation | ✅ Generated |
| `TRAIT_EXTRACTION_GUIDE.md` | This guide | ✅ Generated |

---

**Last updated**: 2025-11-02
**Author**: Claude Code (Anthropic)
**Project**: CAFI Network Analysis, Mo'orea 2019 Survey
**PI**: Adrian Stier Lab
