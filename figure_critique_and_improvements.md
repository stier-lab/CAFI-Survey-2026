# CAFI Survey Figure Critique and Improvement Report

## Executive Summary
This document provides a comprehensive critique of all figures in the CAFI Survey repository against professional graphic design principles, followed by specific improvements to enhance clarity, impact, and publication readiness.

## Figure Design Principles Applied

### 1. Visual Hierarchy & Focus
- Clear primary message/finding immediately apparent
- Strategic use of size, color intensity, and position

### 2. Color Theory & Accessibility
- Colorblind-friendly palettes (viridis, Cividis)
- Consistent color schemes across figures
- Sufficient contrast for print and screen

### 3. Typography & Labeling
- Minimum 10-12pt font for readability
- Sans-serif fonts for clarity
- Clear, descriptive titles and axis labels

### 4. Data-Ink Ratio (Tufte's Principle)
- Maximize data-ink, minimize non-data-ink
- Remove chartjunk and unnecessary elements

### 5. Clarity & Simplicity
- One main message per figure
- Direct labeling when possible
- Strategic use of white space

### 6. Statistical Integrity
- Show confidence intervals/error bars
- Appropriate scales
- Clear sample sizes

### 7. Consistency & Professional Polish
- Uniform style across all figures
- Aligned elements
- Publication-ready resolution (300 DPI)

### 8. Context & Interpretation
- Reference lines and benchmarks
- Annotations for key findings
- Clear units of measurement

## Detailed Figure Critiques by Category

### 1. Community Composition Figures

#### Current Issues:
- **Taxonomic pie chart**: Outdated visualization method, difficult to compare segments
- **Stacked bar charts**: Hard to compare non-baseline categories
- **Species rank abundance**: Missing confidence intervals
- **Heatmaps**: Default color schemes not colorblind-friendly

#### Improvements Needed:
- Replace pie charts with horizontal bar charts
- Use grouped bars instead of stacked where appropriate
- Add confidence intervals to all statistical plots
- Implement viridis color palette consistently

### 2. Diversity Analysis Figures

#### Current Issues:
- **NMDS ordinations**: Overlapping text labels, unclear groupings
- **PCA biplots**: Arrows too small, loadings hard to read
- **Alpha diversity boxplots**: No statistical significance indicators
- **Rarefaction curves**: Line colors too similar

#### Improvements Needed:
- Use ggrepel for non-overlapping labels
- Increase arrow size and add loading labels
- Add significance brackets between groups
- Use distinct, colorblind-friendly line styles

### 3. Coral-CAFI Relationship Figures

#### Current Issues:
- **Scatter plots**: Missing trend lines and confidence bands
- **Correlation matrices**: Hard to distinguish strength levels
- **Multi-panel figures**: Inconsistent axis scales
- **Regression plots**: No R² or p-values displayed

#### Improvements Needed:
- Add loess smoothers with 95% CI
- Use diverging color palette for correlations
- Standardize scales or clearly indicate when different
- Add statistical annotations directly on plots

### 4. Spatial Pattern Figures

#### Current Issues:
- **Maps**: Insufficient contrast between sites
- **Depth distributions**: Bins too wide
- **Site boundaries**: Unclear legend
- **Interpolation maps**: Color scale not intuitive

#### Improvements Needed:
- Increase point size and use distinct shapes
- Optimize bin width using Freedman-Diaconis rule
- Improve legend clarity with better labels
- Use sequential color palette for continuous data

### 5. Network Analysis Figures

#### Current Issues:
- **Network plots**: Node labels overlap
- **Degree distributions**: Missing expected distribution overlay
- **Module visualization**: Colors too similar
- **Keystone ranking**: Bar chart instead of lollipop

#### Improvements Needed:
- Implement force-directed layout with repulsion
- Add power-law fit overlay
- Use qualitative palette with maximum distinction
- Convert to lollipop chart for cleaner appearance

### 6. Machine Learning Figures

#### Current Issues:
- **Variable importance**: Not sorted by importance
- **Prediction plots**: Missing 1:1 reference line
- **Model comparison**: Box plots overlap
- **Feature importance**: No error bars

#### Improvements Needed:
- Sort features by importance descending
- Add diagonal reference line and R²
- Use violin plots with jittered points
- Add bootstrap confidence intervals

### 7. Comprehensive Summary Figures

#### Current Issues:
- **Dashboard layout**: Too crowded, panels too small
- **Text size**: Too small for publication
- **Color consistency**: Different palettes across panels
- **Annotations**: Missing key findings highlights

#### Improvements Needed:
- Reduce number of panels per figure
- Increase base font size to 12pt minimum
- Unify color scheme across all panels
- Add callout boxes for key results

## Implementation Plan

### Priority 1: Critical Issues (Immediate)
1. Fix all text readability issues (font size < 10pt)
2. Replace non-colorblind-friendly palettes
3. Add missing statistical indicators
4. Fix overlapping labels

### Priority 2: Enhancement (Short-term)
1. Standardize color schemes across all figures
2. Add confidence intervals and error bars
3. Improve legend clarity
4. Optimize plot layouts

### Priority 3: Polish (Medium-term)
1. Add annotations for key findings
2. Implement consistent themes
3. Create multi-panel compositions
4. Generate high-resolution exports

## Color Palette Recommendations

### Primary Palette (Colorblind-Safe)
- Site HAU: #E69F00 (Orange)
- Site MAT: #56B4E9 (Sky Blue)
- Site MRB: #009E73 (Teal)

### Taxonomic Groups
- Crab: #CC79A7 (Pink)
- Shrimp: #F0E442 (Yellow)
- Snail: #999999 (Gray)
- Fish: #0072B2 (Blue)

### Sequential Data
- Low to High: Viridis palette
- Diverging: RdBu (reversed)

### Background
- Plot background: White (#FFFFFF)
- Grid lines: Light gray (#F0F0F0)
- Text: Dark gray (#333333)

## Technical Specifications

### Resolution
- Screen: 100 DPI minimum
- Print: 300 DPI minimum
- Vector formats when possible (PDF, SVG)

### Dimensions
- Single column: 80mm (3.15 inches)
- Double column: 169mm (6.65 inches)
- Full page: 229mm (9 inches)

### File Formats
- Primary: PNG for raster, PDF for vector
- Supplementary: SVG for web, TIFF for submission

## Conclusion

The current figures demonstrate strong analytical content but require systematic improvements in visual design for publication readiness. Key priorities are:

1. **Accessibility**: Implement colorblind-friendly palettes throughout
2. **Readability**: Increase font sizes and reduce overlapping elements
3. **Consistency**: Unify themes and color schemes
4. **Clarity**: Add statistical annotations and simplify complex visualizations
5. **Impact**: Highlight key findings with strategic annotations

With these improvements, the figures will meet publication standards for top-tier journals while maximizing clarity and impact for readers.