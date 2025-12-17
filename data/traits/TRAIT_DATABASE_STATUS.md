# CAFI Traits Database - Current Status & Next Steps

**Date**: 2025-11-01
**Total Taxa**: 233
**Processed**: 100% taxonomy, 34.8% associations, 24% functional roles

---

## ✅ COMPLETED

### 1. Efficient Batch Processing (vs. Manual Web Search)
- Built Python pipeline processing 233 taxa in ~4 minutes
- 99.6% WoRMS taxonomy matches (232/233)
- Avoided days of manual lookup

### 2. Pocillopora Association Assignments
- **Source**: Your lab's authoritative reference sheet (`facultative_obligate.xls`, 2012)
- **Obligate**: 17 taxa (Trapezia crabs, obligate gobies, select shrimp)
- **Facultative**: 64 taxa (most Alpheidae, Xanthidae, fish)
- **Evidence-based**: High confidence for common/keystone species

### 3. Functional Role Classification
- 56/233 taxa with assigned roles
- Key groups: mutualist/defenders, commensals, bioeroders, predators, filter feeders
- Genus-level inference where species data unavailable

### 4. Database Integration
- JSON Lines format (full structured data)
- CSV format (analysis-ready)
- Complete provenance (WoRMS AphiaIDs, citations, confidence levels)
- Ready for network analysis integration

---

## 📊 API DATABASES IDENTIFIED

### Immediately Available (Python-friendly)
1. **WoRMS** - ✅ Already using
2. **OBIS** - ⏳ Can add depth ranges
3. **GBIF** - ⏳ Can verify Moorea occurrences

### High-Value (Requires R)
4. **FishBase/SeaLifeBase (rfishbase)** - Best source for size, depth, trophic, diet, reproduction
5. **EOL TraitBank** - 11M trait records across 330+ attributes

### Specialized
6. **Coral Trait Database** - Pocillopora host traits (not CAFI)
7. **Pelagic Species DB** - Not relevant (CAFI are benthic)

**See**: `api_databases_reference.md` for complete documentation

---

## ❌ STILL MISSING

### Life History Traits (0% coverage)
- Max size/length
- Depth range (min/max)
- Trophic level
- Diet composition
- Reproduction mode
- Diel activity

**Why**: No direct Python API for FishBase/SeaLifeBase (R package only)

### Functional Roles (76% unknown)
- Many cryptic/rare species lack literature
- Genus-level inference may be inaccurate
- Need species-specific ecology studies

---

## 🎯 RECOMMENDED NEXT STEPS

### Option A: Quick Python Enhancement (30 minutes)
**Add OBIS depth extraction to pipeline**
- Use OBIS API to get depth ranges from occurrence records
- Validate with GBIF Moorea records
- Improves depth coverage for common species

**Pros**: Fast, Python-only, immediate improvement
**Cons**: Depth data from occurrences (less accurate than FishBase ranges)

**Code ready**: See `api_databases_reference.md` for implementation

---

### Option B: R-Based Trait Extraction (Best Quality, 1-2 hours)
**Install R + rfishbase package and batch query**

1. Install R and rfishbase:
```r
install.packages("rfishbase")
library(rfishbase)
```

2. Run batch queries for all 233 taxa:
```r
# Species info
species_data <- species(our_species_list)

# Ecology (depth, trophic)
ecology_data <- ecology(our_species_list,
                       fields=c("DepthRangeShallow", "DepthRangeDeep",
                               "FoodTroph", "DietTroph"))

# Size data
# Already in species table

# For invertebrates
inverts <- fb_tbl("species", "sealifebase")
```

3. Export to CSV, merge with our database

**Pros**: Best quality, comprehensive traits, authoritative source
**Cons**: Requires R installation and learning curve

**I can create the R script** if you want to run it

---

### Option C: Targeted Manual Extraction (2-3 hours)
**Focus on ~20 keystone/abundant species**

Priority species (from network analysis):
1. *Alpheus diadema* (keystone, degree=12)
2. *Caracanthus maculatus* (keystone)
3. *Alpheus collumianus* (keystone)
4. *Breviturma pica* (keystone)
5. *Trapezia serenei* (most abundant, n=452)
6. *Paragobiodon modestus* (abundant goby, n=69)
7. *Harpiliopsis spinigera* (obligate associate)
8. *Synalpheus charon* (obligate shrimp)
9-20. Other abundant/keystone taxa from network

**Method**: Manual lookup on FishBase/SeaLifeBase web interface
**Pros**: Targets most important species, guaranteed results
**Cons**: Tedious, not scalable, no automation

---

## 💡 MY RECOMMENDATION

**Start with Option A (OBIS depth extraction) NOW**
- Quick win, improves database immediately
- No new tools required
- Can run in this session

**Then Option B (R rfishbase) as follow-up**
- Best long-term solution
- I can create the R script for you to run
- Comprehensive trait coverage

**Skip Option C unless R fails**
- Only as last resort for critical species

---

## WHAT WOULD YOU LIKE TO DO?

**Option 1**: Add OBIS/GBIF depth extraction to Python pipeline (I can do now)
**Option 2**: Create R script for rfishbase batch queries (I create, you run)
**Option 3**: Both (Option 1 now, Option 2 as deliverable for later)
**Option 4**: Something else?

**Your call!**

