# Product Requirements Document (PRD)

## CAFI Survey Analysis Pipeline

**Version**: 1.0
**Date**: November 2025
**Project**: Coral-Associated Fauna Investigation (CAFI) Survey Analysis

---

## 1. Executive Summary

### Purpose
Develop a comprehensive, reproducible analysis pipeline to investigate ecological relationships between *Pocillopora* coral characteristics and their cryptic fauna communities (CAFI) across reef sites in Mo'orea, French Polynesia.

### Scope
- Analyze 114 coral colonies with 2,847 associated fauna individuals
- Implement 17+ statistical and machine learning analyses
- Generate publication-ready figures and tables
- Support research workflow with AI-assisted agents

### Success Criteria
- Complete analytical pipeline with documented outputs
- Reproducible results from raw data to final figures
- Clear identification of ecological patterns and predictors
- Publication-ready visualizations and statistical summaries

---

## 2. Problem Statement

### Background
Coral reefs harbor diverse cryptic fauna communities that are poorly understood due to sampling challenges and taxonomic complexity. Understanding the factors structuring these communities is essential for reef conservation and ecosystem function research.

### Current Challenges
1. **Data complexity**: Multiple data sources (CAFI counts, coral morphology, physiology) require careful integration
2. **Spatial structure**: Fauna distributions show spatial autocorrelation that standard models ignore
3. **Predictor uncertainty**: 30+ potential predictors make it unclear which factors matter most
4. **Taxonomic resolution**: Field identifications lack genetic confirmation

### Target Users
- Coral reef ecologists
- Graduate students analyzing CAFI data
- Collaborators needing reproducible methods

---

## 3. Objectives and Goals

### Primary Objectives

| ID | Objective | Metric |
|----|-----------|--------|
| O1 | Characterize CAFI community structure | Species richness, diversity indices, composition |
| O2 | Identify key coral predictors | Variable importance from multiple models |
| O3 | Quantify spatial patterns | Moran's I, distance decay coefficients |
| O4 | Build predictive models | R², AIC, cross-validation accuracy |
| O5 | Generate publication figures | 100+ figures across 15 analytical themes |

### Secondary Objectives

- Document all methods for reproducibility
- Create reusable analysis templates
- Train ML models for fauna prediction
- Compare statistical approaches (GLMM vs GAM vs ML)

---

## 4. User Stories

### Researcher Stories

1. **As an ecologist**, I want to load and merge all data sources automatically so I can begin analysis without data wrangling.

2. **As a statistician**, I want to compare multiple modeling approaches (GLMM, GAM, RF) so I can identify the most appropriate method.

3. **As a graduate student**, I want documented code with clear comments so I can learn and modify analyses.

4. **As a reviewer**, I want to see reproducible outputs so I can verify the methodology.

### Agent Stories

1. **As the Research PRD Agent**, I need structured problem definitions to guide analysis design.

2. **As the Modeling Agent**, I need data quality reports to recommend appropriate statistical approaches.

3. **As the Figure Factory Agent**, I need analysis results to create publication visualizations.

---

## 5. Functional Requirements

### 5.1 Data Management

| ID | Requirement | Priority |
|----|-------------|----------|
| F1.1 | Load CAFI abundance data with taxonomy | Must |
| F1.2 | Load coral morphology measurements | Must |
| F1.3 | Load coral physiology data | Must |
| F1.4 | Merge datasets by coral_id | Must |
| F1.5 | Create species × site community matrix | Must |
| F1.6 | Handle missing values appropriately | Must |
| F1.7 | Generate data dictionary | Should |

### 5.2 Community Analysis

| ID | Requirement | Priority |
|----|-------------|----------|
| F2.1 | Calculate alpha diversity (richness, Shannon, Simpson) | Must |
| F2.2 | Calculate beta diversity (Bray-Curtis, Jaccard) | Must |
| F2.3 | Perform NMDS ordination | Must |
| F2.4 | Test community differences (PERMANOVA) | Must |
| F2.5 | Identify indicator species | Should |
| F2.6 | Rarefaction curves | Should |

### 5.3 Statistical Modeling

| ID | Requirement | Priority |
|----|-------------|----------|
| F3.1 | Fit GLMMs for abundance with site random effects | Must |
| F3.2 | Fit GLMMs for richness | Must |
| F3.3 | Distance-based redundancy analysis (dbRDA) | Must |
| F3.4 | Spatial regression models | Should |
| F3.5 | GAMs for non-linear relationships | Should |
| F3.6 | Model comparison (AIC, R²) | Must |

### 5.4 Machine Learning

| ID | Requirement | Priority |
|----|-------------|----------|
| F4.1 | Random Forest for variable importance | Must |
| F4.2 | XGBoost predictions | Should |
| F4.3 | Cross-validation | Must |
| F4.4 | Feature importance plots | Must |

### 5.5 Spatial Analysis

| ID | Requirement | Priority |
|----|-------------|----------|
| F5.1 | Moran's I for spatial autocorrelation | Must |
| F5.2 | LISA clusters | Should |
| F5.3 | Distance decay analysis | Must |
| F5.4 | Variograms | Could |
| F5.5 | Spatial weights matrices | Must |

### 5.6 Visualization

| ID | Requirement | Priority |
|----|-------------|----------|
| F6.1 | Publication-quality figures (300 dpi) | Must |
| F6.2 | Consistent color schemes | Must |
| F6.3 | White backgrounds | Must |
| F6.4 | Faceted plots by site/morphotype | Should |
| F6.5 | Summary dashboards | Should |

### 5.7 Network Analysis

| ID | Requirement | Priority |
|----|-------------|----------|
| F7.1 | Co-occurrence networks | Should |
| F7.2 | Module detection | Could |
| F7.3 | Keystone species identification | Could |

---

## 6. Non-Functional Requirements

### 6.1 Performance
- Pipeline should complete in < 30 minutes on standard laptop
- Individual scripts should complete in < 5 minutes

### 6.2 Reliability
- Scripts should handle missing data gracefully
- Failed analyses should not block subsequent scripts
- Error messages should be informative

### 6.3 Usability
- Clear script numbering and naming
- Comprehensive comments in code
- README with quick start instructions

### 6.4 Maintainability
- Modular script design
- Centralized path configuration
- Version control with git

### 6.5 Reproducibility
- Fixed random seeds where applicable
- Session info logging
- Complete dependency documentation

---

## 7. Data Requirements

### Input Data

| File | Records | Variables | Description |
|------|---------|-----------|-------------|
| survey_cafi_data_w_taxonomy_summer2019_v5.csv | 2,847 | ~10 | CAFI individuals with taxonomy |
| survey_coral_characteristics_merged_v2.csv | 114 | ~30 | Coral morphology and metadata |
| survey_master_phys_data_v3.csv | ~50 | ~15 | Coral physiology subset |

### Key Variables

**Response Variables:**
- CAFI abundance (count)
- CAFI richness (count)
- Shannon diversity (continuous)

**Predictor Variables:**
- Coral volume (cm³)
- Coral surface area (cm²)
- Depth (m)
- Site (categorical: HAU, MAT, MRB)
- Morphotype (categorical: verrucosa, meandrina)
- Branch width (categorical: tight, wide)
- Zooxanthellae density (cells/cm²)
- Chlorophyll concentration (µg/cm²)

### Output Data

| Type | Count | Format |
|------|-------|--------|
| Figures | 114 | PNG (300 dpi) |
| Tables | 77 | CSV |
| R Objects | 11 | RDS |
| Reports | 13 | Markdown/HTML |

---

## 8. Technical Architecture

### Technology Stack

**Analysis:**
- R 4.x
- tidyverse ecosystem
- vegan (community ecology)
- lme4 (mixed models)
- randomForest, xgboost (ML)

**Agents:**
- Python 3.x
- Anthropic Claude API
- YAML configuration

**Version Control:**
- Git
- GitHub

### Directory Structure

```
CAFI-Survey-2026/
├── agents/          # Python research agents
├── data/            # Raw input data
├── scripts/         # R analysis pipeline
│   └── utils/       # Helper functions
├── output/          # Generated results
│   ├── figures/
│   ├── tables/
│   ├── objects/
│   └── reports/
└── docs/            # Documentation
```

### Pipeline Architecture

```
[Data Loading] → [Data Cleaning] → [Feature Engineering]
       ↓
[Community Analysis] → [Diversity Metrics]
       ↓
[Statistical Models] → [Model Comparison]
       ↓
[Machine Learning] → [Variable Importance]
       ↓
[Spatial Analysis] → [Autocorrelation Tests]
       ↓
[Visualization] → [Publication Figures]
```

---

## 9. Agent Integration

### Research Workflow Agents

| Agent | Input | Output | Integration Point |
|-------|-------|--------|-------------------|
| Research PRD | Research question | Structured PRD | Project initiation |
| Literature | Topic keywords | Gap analysis | Background research |
| Framework | Hypotheses | Conceptual diagram | Theory development |
| Data QA | Raw data | Quality report | Before analysis |
| EDA | Clean data | Exploratory visualizations | Initial analysis |
| Modeling | EDA results | Model recommendations | Statistical planning |
| Figure Factory | Model results | Publication figures | Visualization |
| Scientific Writer | All outputs | Manuscript sections | Writing |

### Agent Modifications Needed

1. **Data paths**: Update to use project-specific paths
2. **Output formats**: Align with R pipeline outputs
3. **Context**: Add CAFI-specific domain knowledge
4. **Prompts**: Customize for coral ecology terminology

---

## 10. Risks and Mitigations

| Risk | Impact | Probability | Mitigation |
|------|--------|-------------|------------|
| Missing R packages | Analysis fails | Medium | Graceful package checks |
| Spatial autocorrelation | Invalid p-values | High | Use spatial models |
| Morphotype collapse | No comparisons | High | Restore original categories |
| Memory limits | Pipeline crashes | Low | Process in chunks |
| API rate limits | Agent delays | Medium | Implement retries |

---

## 11. Timeline and Milestones

### Phase 1: Foundation (Complete)
- [x] Data loading and cleaning
- [x] Community composition analysis
- [x] Basic diversity metrics

### Phase 2: Core Analysis (In Progress)
- [x] GLMM modeling
- [x] Spatial pattern analysis
- [ ] Complete network analysis
- [ ] Fix comprehensive predictor analysis

### Phase 3: Advanced (Pending)
- [ ] Full ML pipeline
- [ ] Agent integration
- [ ] Manuscript preparation

### Phase 4: Publication (Future)
- [ ] Peer review
- [ ] Revisions
- [ ] Data archiving

---

## 12. Success Metrics

### Quantitative Metrics

| Metric | Target | Current |
|--------|--------|---------|
| Scripts completed | 17/17 | 13/17 |
| Figures generated | 120 | 114 |
| Tables generated | 80 | 77 |
| Test coverage | 80% | 0% |
| Documentation pages | 15 | 13 |

### Qualitative Metrics

- Code is readable and well-commented
- Results are reproducible from raw data
- Visualizations are publication-ready
- Methods are statistically appropriate

---

## 13. Appendices

### A. Package Dependencies

See `scripts/00_load_libraries.R` for complete list.

### B. Variable Definitions

See `output/tables/data_dictionary.csv`.

### C. Analysis Log

See `output/reports/PIPELINE_STATUS.md`.

---

*Document prepared for CAFI Survey Analysis Pipeline v1.0*
