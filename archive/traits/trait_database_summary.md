# CAFI Traits Database - Final Summary

**Generated**: 2025-11-01
**Database version**: Updated with Stier Lab reference sheet (2012)
**Total taxa**: 233

---

## Coverage Summary

### Taxonomy (WoRMS)
- **Resolved**: 232/233 (99.6%)
- **Unresolved**: 1/233 (0.4%)
- **AphiaIDs**: 232 assigned
- **Taxonomic status**: All accepted or valid

### Pocillopora Association
- **Obligate**: 17/233 (7.3%)
- **Facultative**: 64/233 (27.5%)
- **Unknown**: 152/233 (65.2%)

**Evidence source**: Stier Lab reference sheet (facultative_obligate.xls, 2012) + literature

### Functional Roles
- **Assigned**: 56/233 (24.0%)
- **Unknown**: 177/233 (76.0%)

**Top functional roles**:
1. Commensal; bioeroder (24 taxa) - Alpheidae shrimp
2. Mutualist/defender; grazer/herbivore (10 taxa) - Trapeziidae crabs
3. Commensal (7 taxa) - Various palaemonid shrimp
4. Commensal; filter_feeder (4 taxa) - Galatheidae squat lobsters
5. Detritivore; scavenger (4 taxa) - Ophiuroidea brittle stars

---

## Key Taxonomic Groups

### OBLIGATE POCILLOPORA ASSOCIATES (17 taxa)

**Trapeziidae (Guard Crabs)** - 10 species
- *Trapezia serenei* (n=452 individuals) - **Most abundant**
- *T. rufopunctata*, *T. bidentata*, *T. tigrina*, *T. bella*
- *T. guttata*, *T. speciosa*, *T. areolata*, *T. flavopunctata*, *T. punctimanus*
- **Role**: Mutualist/defender + grazer/herbivore
- **Evidence**: Defend corals against Drupella, Culcita; graze coral lipids (Stewart et al. 2014)

**Palaemonidae (Coral Shrimp)** - 4 species
- *Harpiliopsis spinigera*, *H. depressa*
- *Fennera chacei*
- **Role**: Commensal
- **Evidence**: Free-living in coral branches

**Alpheidae (Snapping Shrimp)** - 2 species
- *Alpheus lottini*, *Synalpheus charon*
- **Role**: Commensal; bioeroder
- **Note**: Other Alpheus spp. coded as **facultative**

**Gobiidae (Coral Gobies)** - 1 species
- *Paragobiodon modestus* (n=69)
- **Role**: Commensal; mucus feeder
- **Evidence**: Obligate coral dwellers (Munday 2004)

---

### FACULTATIVE POCILLOPORA ASSOCIATES (64 taxa)

**Alpheidae (Snapping Shrimp)** - 22 species
- **Keystone**: *Alpheus diadema* (n=16, keystone index=9.61)
- Also: *A. collumianus* (keystone), *A. pachychirus*, *A. parvirostris*, etc.
- **Role**: Commensal; bioeroder
- **Association**: Facultative (found in many coral types)

**Xanthidae (Crabs)** - 12 species
- Including *Cymo*, *Domecia*, *Pilodius*, *Chlorodiella*
- **Role**: Mostly commensal
- **Association**: Facultative coral associates

**Palaemonidae (Shrimp)** - 6 species
- Including *Periclimenes*, *Cuapetes*, *Chlorocurtis*
- **Role**: Commensal
- **Association**: Facultative (found in various coral/crevice habitats)

**Galatheidae (Squat Lobsters)** - 3 species
- *Galathea mauritiana*, *Coralliogalathea humilis*, *Phylladiorhynchus*
- **Role**: Filter feeder; some corallivores
- **Association**: Facultative

**Fish** - Various families
- **Cirrhitidae**: *Paracirrhites arcatus*, *Caracanthus maculatus* (keystone)
- **Gobiidae**: *Paragobiodon*, *Eviota*
- **Role**: Predators (hawkfish), mucus feeders (gobies)

---

## Network Keystone Species (from Script 06)

| Rank | Species | Functional Role | Association | Network Metrics |
|------|---------|----------------|-------------|-----------------|
| 1 | *Alpheus diadema* | Commensal; bioeroder | **Facultative** | Degree=12, Betweenness=260 |
| 2 | *Caracanthus maculatus* | Predator | **Facultative** | Degree=7, Betweenness=231 |
| 3 | *Alpheus collumianus* | Commensal; bioeroder | **Facultative** | Degree=8, Betweenness=170 |
| 4 | *Breviturma pica* | Unknown | **Unknown** | Degree=6, Betweenness=242 |
| 5 | *Periclimenes* sp. | Commensal | **Unknown** | Degree=7, Betweenness=190 |

**Key insight**: Most keystone species are **facultative** associates, not obligate. This suggests network structure is driven by **generalist** species that occur across many coral/habitat types.

---

## Data Quality & Confidence

### HIGH Confidence (99 taxa)
- **Taxonomy**: All WoRMS-matched taxa
- **Association** (exact match): Trapeziidae, obligate Palaemonidae, select Alpheidae
- **Functional role**: Trapezia (mutualist/defender), Paragobiodon (mucus feeder)

### MEDIUM Confidence (64 taxa)
- **Association** (genus-level): Most Alpheidae, some Xanthidae
- **Functional role** (inferred from family): Galatheidae, Ophiuroidea, general Alpheidae

### LOW Confidence (70 taxa)
- **Unknown associations**: Taxa not in reference sheet
- **Unknown functional roles**: Taxa without genus/family-level data

---

## Missing Data & Limitations

### Life History Traits (0% coverage)
- ❌ **Size** (max_length_mm): No data
- ❌ **Depth range**: No data
- ❌ **Trophic level**: No data
- ❌ **Reproduction**: No data
- ❌ **Diel activity**: No data

**Reason**: FishBase/SeaLifeBase require manual web scraping; not accessible via public API

**Recommendation**: Prioritize keystone species (~20 taxa) for manual trait extraction from FishBase/SeaLifeBase

### Functional Roles (76% unknown)
- Many species lack specific ecological literature
- Inferred roles based on family/genus may be inaccurate
- Need species-specific studies for cryptic/rare taxa

**Recommendation**:
1. Cross-reference with Moorea Biocode for habitat associations
2. Consult with taxonomic experts (Adrian Stier, Chris Gotschalk)
3. Flag high-priority unknowns for targeted literature review

---

## Data Files

### Primary Outputs
- **`cafi_traits_updated.jsonl`** - Full trait records (JSON Lines format)
- **`cafi_traits_updated.csv`** - Flattened table for analysis
- **`association_lookup.csv`** - Reference sheet (obligate/facultative assignments)

### Supporting Files
- **`trait_database_summary.md`** - This summary
- **`pocillopora_association_evidence.md`** - Evidence review & coding scheme
- **`facultative_obligate.csv`** - Converted reference sheet

---

## Usage in Network Analysis

The traits database is now integrated with network analysis (Script 06):

1. **Association type** (obligate vs. facultative) can be tested as predictor of:
   - Network degree
   - Module membership
   - Co-occurrence patterns

2. **Functional roles** can be used to:
   - Test if mutualists have different network positions than commensals
   - Examine functional redundancy within modules
   - Predict effects of species loss

3. **Taxonomic groups** (Alpheidae, Trapeziidae, etc.) provide:
   - Family-level network patterns
   - Phylogenetic signal in co-occurrence
   - Functional group diversity metrics

---

## Next Steps

### Immediate (High Priority)
1. ✅ Integrate reference sheet associations → **DONE**
2. ⏳ Extract size/depth for keystone species (manual FishBase/SeaLifeBase)
3. ⏳ Add functional role evidence for abundant/keystone taxa
4. ⏳ Validate obligate/facultative assignments with domain expert

### Future (Medium Priority)
1. Cross-reference Moorea Biocode for additional trait data
2. Systematic literature review for species-specific ecology
3. Add length-weight parameters for biomass calculations
4. Compile diet/trophic data for food web analysis

### Long-term (Low Priority)
1. Complete trait coverage for all 233 taxa
2. Add reproductive mode, larval type, dispersal potential
3. Integrate with phylogenetic data (if available)
4. Create interactive trait visualization dashboard

---

## Citation

**Data sources**:
- World Register of Marine Species (WoRMS) - Taxonomy
- Stier Lab reference sheet (2012) - Pocillopora associations
- Stewart et al. 2014 (PeerJ) - Trapezia defensive behavior
- Munday et al. 2004 - Coral goby ecology
- Hernandez et al. 2021 (MDPI Diversity) - Trapezia-Pocillopora symbiosis

**Database compiled by**: Claude Code (Anthropic)
**Project**: CAFI Network Analysis, Mo'orea 2019 Survey
**PI**: Adrian Stier Lab

---

**Database status**: ✅ Ready for network analysis integration
**Last updated**: 2025-11-01
