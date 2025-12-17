# **CAFIs Survey Paper – Full Figure Plan with Hypotheses & Statistical Approaches**

Revised to reflect the functional groups **we confidently understand biologically**:
- **Trapezia (Mutualist Defenders)** – protect corals from predators, sediments, fouling.
- **Resident Fishes (Nutrient Providers)** – enhance coral condition via ammonium excretion and microflow.
- **Corallivores (Snail predators: e.g., *Drupella*)** – negative effects: tissue loss, disease transmission.
- **Other Invertebrates / Cryptofauna** – included statistically, but functional roles unknown.

Below is the *publication-ready figure plan*, now including:
- Clear hypotheses for each figure
- Full statistical approaches (model structures, transformations, and decisions)
- Alignment with functional groups we know

---

# **FIGURE 1 — Study System + CAFI Overview + Conceptual Scaling**

### **Purpose**
Introduce the study system, key CAFI taxa, and the conceptual framework.

### **Hypotheses**
1. Coral size and local coral density influence CAFI colonization.
2. Two competing conceptual expectations:
   - **H1: Field‑of‑Dreams** — CAFI abundance rises *proportionally* with coral size.
   - **H2: Propagule Redirection** — isolated corals or small patches receive *disproportionately more* CAFI.

### **Stats / Methods**
- Conceptual only; no analytics.
- Choose functional exemplars: **Trapezia**, resident fishes, corallivores.

---

# **FIGURE 2 — Landscape Characteristics → Scaling (Abundance, Richness, Diversity)**

### **Purpose**
Test how coral size and proximity shape total CAFI abundance, species richness, and Shannon diversity.

### **Hypotheses**
- **H1:** Larger corals host more CAFI (positive size–abundance slope).
- **H2:** Isolated corals host more CAFI (propagule redirection prediction).
- **H3:** Diversity increases with coral size but may decrease with proximity if redirection strongly limits colonization.

### **Statistical Approaches**
- **Response variables:** total abundance, richness, Shannon.
- **Predictors:** coral size (continuous), nearest‑neighbor distance, site (random or fixed factor).
- **Model:** GLM or negative binomial GLM for abundance; Gaussian/GLM for richness/diversity.
- **Formula:**
  - Abundance: `NB glm: CAFI_abund ~ log(size) + proximity + site`
  - Richness: `Poisson/NB glm: richness ~ log(size) + proximity`
  - Diversity: `lm: Shannon ~ size + proximity`
- **Plotting rule:** *If the relationship is nonsignificant, remove regression line and show raw trend only.*

---

# **FIGURE 3 — Functional & Taxonomic Group Responses + Taxon-Specific Slopes**

### **Purpose**
Identify which taxa and functional groups drive landscape patterns.

### **Functional Groups Used in Analysis:**
- **Trapezia (Defenders)** – expected positive effect on coral health; may show weak redirection.
- **Resident Fishes** – nutrient providers; may increase with size (more space) rather than isolation.
- **Corallivores (Drupella)** – expected to increase with coral size (more food) but decrease with isolation.
- **Other taxa** – included without functional assumptions.

### **Hypotheses**
- **H1:** Trapezia abundance increases with coral size; may not increase with isolation.
- **H2:** Resident fishes show positive size scaling but weak proximity effects.
- **H3:** Corallivores may prefer larger colonies but not isolated ones.
- **H4:** Taxon‑specific slopes reveal mechanistic scaling differences.

### **Statistical Approaches**
- **Group-level models:**
  - `NB glm: group_abund ~ log(size) + proximity`

- **Species-specific slopes:**
  - Fit GLMs per species (top 12–20 taxa):
    `glm or NB glm: species_abundance ~ log(size) + proximity`
  - Extract slope coefficients.

- **Heatmap:** species × slope estimates.

---

# **FIGURE 4 — Species Compositional Changes + Incidence + Co-occurrence Networks**

### **Purpose**
Understand how CAFI species composition changes across landscape gradients.

### **Hypotheses**
- **H1:** Species composition shifts strongly with coral size and proximity.
- **H2:** Sites may differ in composition but not strongly (shared habitat context).
- **H3:** Co-occurrence networks reveal nonrandom species associations (e.g., Trapezia clusters with certain cryptofauna; corallivores occur without defenders).

### **Statistical Approaches**
- **Ordination:** NMDS or PCoA using **Hellinger-transformed abundance**.
- **Tests:** PERMANOVA (`adonis2`) to test effects of size, proximity, site.
- **Incidence plots:** frequency of species presence across size bins.
- **Co-occurrence network:**
  - Use probabilistic co-occurrence or correlation-based edges.
  - Identify clusters via modularity.

---

# **FIGURE 5 — Coral Condition Across Landscape Characteristics**

### **Purpose**
Test whether coral physiological condition varies with size or landscape position, independent of CAFI.

### **Hypotheses**
- **H1:** Larger corals may show higher biomass but lower symbiont density (common density-dependent pattern).
- **H2:** Isolated corals may experience better/worse condition depending on CAFI-mediated processes.
- **H3:** PCA summarizing condition metrics reveals overall quality gradients.

### **Statistical Approaches**
- **PCA:** protein, carbohydrate, AFDM, zooxanthellae.
- **Models:**
  - Condition PC1 ~ size + proximity
  - Individual metrics ~ size + proximity
- **Model type:** Gaussian linear models.
- **Visualization:** scatterplots with regression when significant.

---

# **FIGURE 6 — Linking CAFI Abundance, Diversity & Composition to Coral Condition (Feedbacks)**

### **Purpose**
Directly test whether CAFI communities influence coral condition.

### **Hypotheses**
- **H1:** More defenders (Trapezia) → higher condition.
- **H2:** More resident fishes → improved condition via nutrient provisioning.
- **H3:** More corallivores → lower condition.
- **H4:** CAFI community composition (PC1) predicts coral condition more strongly than abundance alone.
- **H5:** Landscape may modulate these effects (e.g., mutualist benefits weaken on very large colonies).

### **Statistical Approaches**
- **Models:**
  - Condition PC1 ~ total CAFI abundance + CAFI diversity
  - Condition PC1 ~ defender abundance + fish abundance + corallivore abundance
  - Condition PC1 ~ CAFI composition (PC1, PC2)

- **SEM (optional):**
  - Size → CAFI → Condition
  - Proximity → CAFI → Condition

- **Visualization:**
  - Lollipop plots of functional group effect sizes.
  - Partial residual plots.

---

# **Summary of Analytical Logic**
1. **Landscape structure** → affects CAFI abundance & diversity.
2. **Landscape structure** → affects which CAFI taxa arrive and assemble.
3. **CAFIs** → affect coral condition based on their functional roles.
4. **Feedbacks** emerge where coral state influences future CAFI assembly.

This figure plan now reflects:
- the biology we actually know,
- clean mechanistic hypotheses,
- matched statistical workflows,
- publishable, logically sequenced analyses.

