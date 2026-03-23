# Metadata: survey_coral_haplotypes_v1.csv

## Overview

Symbiodiniaceae ITS2 haplotype assignments for 114 *Pocillopora* coral colonies surveyed across 3 reef sites in Mo'orea, French Polynesia (Summer 2019). Haplotypes were determined via PCR amplification of the ITS2 region and mitochondrial open reading frame (mtORF) restriction fragment analysis (Sac1, PocHis_XhoI digests).

## Provenance

- **Source**: Extracted from `CAFI_2025/data/bco-dmo-ready/moorea_survey_physiology_master_2019_v3.csv` (haplotype, dna_collection_date columns) with cross-validation against the original genotype file (`genotype_9_2_2022.csv`)
- **Original genotyping**: PCR amplification of ITS2 barcode region, completed September 2022
- **DNA collection period**: July–August 2019
- **Lab processing**: DNA extraction January 2021, PCR plate assignment documented in `metadata_for_genetic_samples_v3.xlsx` (archived in CAFI_2025)

## Column Descriptions

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `coral_id` | character | Unique coral colony identifier (matches all other survey data files) | "HAU-POC01" |
| `site` | character | Reef site code: HAU (Hauru), MAT (Maatea), MRB (Maharepa barrier reef) | "HAU" |
| `haplotype` | character | Symbiodiniaceae ITS2 haplotype assignment (see Haplotype Key below) | "1a_Pm" |
| `dna_collection_date` | date | Date tissue sample was collected for DNA extraction (YYYY-MM-DD) | "2019-07-04" |

## Haplotype Key

| Haplotype | Description | N (this dataset) |
|-----------|-------------|-----------------|
| `1a_Pe` | Clade C, *Pocillopora eydouxi*-type symbiont | 49 |
| `1a_Pm` | Clade C, *Pocillopora meandrina*-type symbiont | 32 |
| `10` | Clade C, haplotype 10 | 7 |
| `3b` | Clade C, haplotype 3b | 6 |
| `8a` | Clade C, haplotype 8a | 2 |
| `3f` | Clade C, haplotype 3f | 2 |
| `1a` | Clade C, haplotype 1a (unresolved Pe/Pm) | 1 |
| `3a` | Clade C, haplotype 3a | 1 |
| `3h` | Clade C, haplotype 3h | 1 |
| `Did not amplify` | PCR amplification failed | 11 |
| `No sample` | No tissue sample collected | 2 |

## Data Quality

- **114 corals total** (all corals in the survey)
- **101 with valid haplotype** (88.6% success rate)
- **11 did not amplify** (PCR failure)
- **2 no sample** (tissue not collected)
- All coral_ids match `survey_coral_characteristics_merged_v2.csv` and can be joined via `coral_id`

## Site Distribution

| Site | N corals | Valid haplotypes |
|------|----------|-----------------|
| HAU (Hauru) | 38 | 27 |
| MAT (Maatea) | 39 | 33 |
| MRB (Maharepa) | 37 | 41* |

*Counts reflect successful amplification per site.

## Usage in R

```r
haplotypes <- read_csv("data/survey_coral_haplotypes_v1.csv")

# Join to coral characteristics
coral <- read_csv("data/survey_coral_characteristics_merged_v2.csv") %>%
  left_join(haplotypes, by = c("coral_id", "site"))

# Filter to valid haplotypes only
coral_with_hap <- coral %>%
  filter(!haplotype %in% c("Did not amplify", "No sample"))
```

## Related Files

- `survey_coral_characteristics_merged_v2.csv` — coral morphology and GPS (join by `coral_id`)
- `survey_master_phys_data_v3.csv` — coral physiology (also contains haplotype column, but only 108 of 114 corals)
- `survey_cafi_data_w_taxonomy_summer2019_v5.csv` — CAFI specimen records (join by `coral_id`)

## Citation

Stier AC, Primo A, Curtis JS, Osenberg CW (2026). Colony size drives sublinear scaling, compositional turnover, and biodiversity-condition feedbacks in coral-associated fauna. [in preparation]

## Contact

Adrian Stier | astier@ucsb.edu | Stier Lab, UC Santa Barbara
