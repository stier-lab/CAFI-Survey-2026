# Marine Species Trait Databases with API Access

**Updated**: 2025-11-01
**Purpose**: Identify databases with programmatic access for CAFI trait extraction

---

## Summary of Available APIs

| Database | Taxa Coverage | Traits Available | API Type | Access | Cost |
|----------|--------------|------------------|----------|---------|------|
| **WoRMS** | All marine | Taxonomy only | REST API | ✅ Open | Free |
| **OBIS** | All marine | Occurrence, depth, lat/lon | REST API | ✅ Open | Free |
| **GBIF** | All life | Occurrence, traits (limited) | REST API | ✅ Open | Free |
| **FishBase/SeaLifeBase** | Fish + inverts | Size, depth, trophic, diet, reproduction | R package (rfishbase) | ⚠️ R only | Free |
| **EOL TraitBank** | All life | 330+ attributes, 11M records | R package (traits) | ⚠️ R only | Free |
| **Coral Trait Database** | Scleractinian corals | 150+ traits, life history | R package (traits) / Web | ✅ Open | Free |
| **Pelagic Species DB** | Pelagic fish/inverts | 529 species, functional traits | Download | ⚠️ Static | Free |

---

## 1. WoRMS (World Register of Marine Species)

**URL**: https://www.marinespecies.org/rest
**Coverage**: All marine organisms (~240,000 accepted names)
**Traits**: Taxonomy, synonyms, classification, expert curators

### API Details
- **Type**: REST API (JSON)
- **No authentication required**
- **Rate limit**: Please respect (0.5-1s between requests)
- **Python-friendly**: Direct HTTP requests

### Endpoints (Already Using ✅)
```
GET /AphiaRecordsByName/{scientificname}
GET /AphiaRecordByAphiaID/{aphiaid}
GET /AphiaClassificationByAphiaID/{aphiaid}
GET /AphiaSynonymsByAphiaID/{aphiaid}
```

### Example
```python
import requests
response = requests.get("https://www.marinespecies.org/rest/AphiaRecordsByName/Trapezia%20serenei?like=false&marine_only=true")
data = response.json()
# Returns: AphiaID, taxonomic classification, status, etc.
```

**Status**: ✅ **Already integrated in our database**

---

## 2. OBIS (Ocean Biodiversity Information System)

**URL**: https://api.obis.org/
**Coverage**: Marine species occurrences (global)
**Traits**: Depth range, geographic distribution, habitat associations

### API Details
- **Type**: REST API (JSON)
- **No API key required**
- **Python-friendly**: Direct HTTP requests

### Useful Endpoints
```
GET /occurrence?scientificname={name}
GET /taxon/{aphiaid}
```

### Fields Available
- `minimumDepthInMeters`, `maximumDepthInMeters`
- `decimalLatitude`, `decimalLongitude`
- `habitat` (reef, soft bottom, etc.)
- `establishmentMeans` (native, invasive)

### Example
```python
import requests
url = "https://api.obis.org/occurrence"
params = {
    'scientificname': 'Alpheus diadema',
    'fields': 'minimumDepthInMeters,maximumDepthInMeters'
}
response = requests.get(url, params=params)
data = response.json()
# Returns: depth observations for the species
```

**Usefulness**: ⭐⭐⭐⭐
- Can extract depth ranges for species
- Verify Moorea occurrence records
- Get habitat associations

**Action**: Add OBIS depth extraction to trait pipeline

---

## 3. GBIF (Global Biodiversity Information Facility)

**URL**: https://api.gbif.org/v1/
**Coverage**: All taxa (140M+ occurrence records)
**Traits**: Occurrences, some trait data, links to other databases

### API Details
- **Type**: REST API (JSON)
- **No API key required** (but rate limits apply)
- **Python-friendly**

### Endpoints
```
GET /species/match?name={scientificname}
GET /species/{key}
GET /occurrence/search?taxonKey={key}
```

### Example
```python
import requests
# Match species
match = requests.get("https://api.gbif.org/v1/species/match",
                    params={'name': 'Trapezia serenei'}).json()
taxon_key = match['usageKey']

# Get occurrences
occs = requests.get("https://api.gbif.org/v1/occurrence/search",
                   params={'taxonKey': taxon_key, 'country': 'PF'}).json()
# PF = French Polynesia
```

**Usefulness**: ⭐⭐⭐
- Good for verifying species presence in Moorea/French Polynesia
- Links to FishBase, WoRMS, EOL via taxonKey
- Limited trait data directly

---

## 4. FishBase / SeaLifeBase (via rfishbase)

**URL**: https://fishbase.ropensci.org/ (data hosting)
**Coverage**:
- FishBase: 35,000+ fish species
- SeaLifeBase: Invertebrates, aquatic non-fish

**Traits**: Size (TL, SL), depth, trophic level, diet, reproduction, habitat

### API Details
- **Type**: R package (`rfishbase`) accessing HuggingFace parquet files
- **No API key required**
- **⚠️ R-only** (no direct Python HTTP API documented)

### Available Functions (R)
```r
library(rfishbase)

# Species info
species(c("Alpheus diadema"))

# Ecology (depth, habitat)
ecology("Alpheus diadema",
        fields=c("SpecCode", "DietTroph", "FoodTroph", "DepthRangeShallow", "DepthRangeDeep"))

# Length-weight parameters
popgrowth("Alpheus diadema")

# Reproduction
reproduction("Alpheus diadema")

# For invertebrates, use SeaLifeBase
fb_tbl("species", "sealifebase")
```

### Fields Available
- **Size**: `Length` (max total length), `CommonLength`, length-weight `a`, `b`
- **Depth**: `DepthRangeShallow`, `DepthRangeDeep`
- **Trophic**: `FoodTroph`, `DietTroph`, trophic SE
- **Diet**: Diet composition, food items
- **Reproduction**: Spawning, fecundity, sexual system
- **Habitat**: Reef-associated, substrate preference

**Usefulness**: ⭐⭐⭐⭐⭐
- **BEST SOURCE** for fish/invertebrate life history traits
- Comprehensive, curated, peer-reviewed

**Challenge**: Requires R environment or manual web scraping

**Solutions**:
1. **Install R + rfishbase** and run batch queries
2. **Web scraping** (less preferred, against terms of service)
3. **Manual lookup** for ~20 keystone species (feasible)

**Action**: Install R, run rfishbase queries for our 233 taxa

---

## 5. Encyclopedia of Life (EOL) TraitBank

**URL**: https://eol.org/traitbank
**Coverage**: 1.7M taxa, 11M trait records
**Traits**: 330+ attributes (size, diet, habitat, behavior, morphology)

### API Details
- **Type**: SPARQL endpoint + R package (`traits`)
- **No API key required**
- **R package available**: `traits::traitbank()`

### Access Methods

**Option 1: R package `traits`**
```r
library(traits)

# Search for species traits
result <- traitbank(query = "Trapezia serenei")

# Filter by trait type
size_traits <- traitbank_trait_search("body size")
```

**Option 2: SPARQL query (advanced)**
```
# EOL provides SPARQL endpoint for custom queries
# More complex but allows precise trait filtering
```

### Available Traits (Examples)
- Body size, body mass, length
- Habitat (reef, rocky, sand)
- Diet (carnivore, omnivore, herbivore)
- Diel activity (diurnal, nocturnal)
- Reproduction mode
- Lifespan
- Trophic level

**Usefulness**: ⭐⭐⭐⭐
- Aggregates data from multiple sources
- Good coverage for common species
- Some gaps for cryptic/rare CAFI

**Challenge**: R-only access via package

**Action**: Use R `traits` package to query EOL for CAFI taxa

---

## 6. Coral Trait Database

**URL**: https://coraltraits.org
**Coverage**: 1,500+ coral species (Scleractinia)
**Traits**: 150+ traits including growth, reproduction, symbiosis, morphology

### API Details
- **Type**: Web interface + R package (`traits`)
- **Download**: Full database available as CSV
- **R access**: `coraltraits::coral_traits()`

### Example
```r
library(coraltraits)

# Get all traits for Pocillopora
poc_traits <- coral_traits(trait_name = NULL,
                           coral_species = "Pocillopora")
```

**Usefulness for CAFI**: ⭐⭐
- **Corals only**, not CAFI species
- Useful for understanding **host traits** (Pocillopora characteristics)
- Could correlate CAFI networks with Pocillopora traits

**Action**: Extract Pocillopora traits to test host trait effects on CAFI

---

## 7. Pelagic Species Trait Database (2024)

**URL**: https://www.nature.com/articles/s41597-023-02689-9
**Coverage**: 529 pelagic fish + invertebrate species
**Traits**: Functional traits for pelagic taxa

### Access
- **Static dataset** (Zenodo download)
- No API

**Usefulness**: ⭐
- Limited to pelagic species
- CAFI are benthic/reef-associated (not pelagic)
- Unlikely overlap with our taxa

---

## Recommended Implementation Strategy

### Phase 1: Immediate (Python-friendly APIs)
1. ✅ **WoRMS** - Already done
2. ⏳ **OBIS** - Add depth range extraction
3. ⏳ **GBIF** - Cross-reference Moorea occurrences

### Phase 2: R-based Trait Extraction (High Value)
4. ⏳ **rfishbase** - Install R, query FishBase/SeaLifeBase for:
   - Size (max length)
   - Depth range (min/max)
   - Trophic level
   - Diet summary
   - Reproduction
5. ⏳ **EOL TraitBank** (via `traits` package) - Additional trait coverage

### Phase 3: Manual Augmentation
6. ⏳ **Web scraping** for keystone species (~20 taxa) if R queries incomplete
7. ⏳ **Literature review** for species-specific ecology

---

## Code Example: OBIS Depth Extraction

```python
import requests
import pandas as pd
import time

def get_obis_depth(species_name):
    """Get depth range from OBIS for a species."""
    url = "https://api.obis.org/occurrence"
    params = {
        'scientificname': species_name,
        'fields': 'minimumDepthInMeters,maximumDepthInMeters'
    }

    try:
        time.sleep(0.5)  # Rate limiting
        response = requests.get(url, params=params, timeout=10)
        data = response.json()

        if 'results' in data and len(data['results']) > 0:
            depths = [r for r in data['results'] if
                     'minimumDepthInMeters' in r or 'maximumDepthInMeters' in r]

            if depths:
                min_depths = [d.get('minimumDepthInMeters') for d in depths if d.get('minimumDepthInMeters')]
                max_depths = [d.get('maximumDepthInMeters') for d in depths if d.get('maximumDepthInMeters')]

                return {
                    'depth_min_m': min(min_depths) if min_depths else None,
                    'depth_max_m': max(max_depths) if max_depths else None,
                    'n_records': len(depths)
                }
    except Exception as e:
        print(f"  Error for {species_name}: {e}")

    return {'depth_min_m': None, 'depth_max_m': None, 'n_records': 0}

# Example usage
depth_data = get_obis_depth("Alpheus diadema")
print(depth_data)
```

---

## Next Steps

**Immediate Action** (this session):
1. Add OBIS depth extraction to trait pipeline
2. Cross-reference GBIF for Moorea occurrence validation

**Follow-up** (requires R installation):
1. Install R + rfishbase package
2. Create R script to batch query all 233 CAFI taxa
3. Export results to CSV
4. Merge with existing trait database

**Would you like me to:**
- A) Add OBIS/GBIF queries to the Python pipeline now (quick)
- B) Create an R script for rfishbase batch queries (requires R install)
- C) Focus on manual extraction for ~20 keystone species via web scraping

