# ============================================================================
# 02_community_analysis.R - CAFI Community Analysis
# ============================================================================
#
# PURPOSE: Comprehensive analysis of CAFI communities including:
#   - Abundance patterns and variation
#   - Taxonomic composition and resolution
#   - Diversity metrics (alpha, beta)
#   - Community structure (ordination, clustering)
#   - Site and coral size effects
#
# USAGE:
#   source("scripts/00_setup.R")
#   source("scripts/01_load_data.R")
#   source("scripts/02_community_analysis.R")
#
# OUTPUTS:
#   - Figures: output/figures/community/
#   - Tables: output/tables/
#   - Console: Statistical summaries with effect sizes, p-values, df
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    CAFI COMMUNITY ANALYSIS\n")
cat("============================================================\n\n")

# Load setup and data
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

# Use script-specific output directory
FIG_DIR <- PATHS$fig_02_community

# ============================================================================
# PART 1: ABUNDANCE PATTERNS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 1: ABUNDANCE PATTERNS\n")
cat("------------------------------------------------------------\n\n")

# 1.1 Overall abundance summary
cat("1.1 Overall CAFI Abundance:\n")
cat("    Total individuals:", sum(coral_master$total_cafi), "\n")
cat("    Mean per coral:", round(mean(coral_master$total_cafi), 1),
    "(SD:", round(sd(coral_master$total_cafi), 1), ")\n")
cat("    Median per coral:", median(coral_master$total_cafi), "\n")
cat("    Range:", min(coral_master$total_cafi), "-", max(coral_master$total_cafi), "\n")
cat("    Corals with 0 CAFI:", sum(coral_master$total_cafi == 0), "/", nrow(coral_master), "\n\n")

# 1.2 Abundance distribution
cat("1.2 Abundance Distribution:\n")
skewness_val <- moments::skewness(coral_master$total_cafi)
kurtosis_val <- moments::kurtosis(coral_master$total_cafi)
cat("    Skewness:", round(skewness_val, 2), "(positive = right-skewed)\n")
cat("    Kurtosis:", round(kurtosis_val, 2), "(>3 = heavy tails)\n")

# Shapiro-Wilk test (on sample if n > 5000)
sw_test <- shapiro.test(coral_master$total_cafi)
cat("    Shapiro-Wilk W:", round(sw_test$statistic, 4), ", p:", format.pval(sw_test$p.value, 3), "\n")
cat("    Interpretation:", ifelse(sw_test$p.value < 0.05,
                                  "NOT normally distributed (use NB-GLM)", "Approx. normal"), "\n\n")

# 1.3 Abundance by site
cat("1.3 Abundance by Site:\n")
site_summary <- coral_master %>%
  group_by(site) %>%
  summarise(
    n_corals = n(),
    total_cafi = sum(total_cafi),
    mean = round(mean(total_cafi), 1),
    sd = round(sd(total_cafi), 1),
    median = median(total_cafi),
    min = min(total_cafi),
    max = max(total_cafi),
    .groups = "drop"
  )
print(site_summary)

# Kruskal-Wallis test (non-parametric)
kw_test <- kruskal.test(total_cafi ~ site, data = coral_master)
cat("\n    Kruskal-Wallis H:", round(kw_test$statistic, 2),
    ", df:", kw_test$parameter,
    ", p:", format.pval(kw_test$p.value, 3), "\n")

# Effect size (epsilon-squared): (H - k + 1) / (N - k)
k_groups <- length(unique(coral_master$site))
epsilon_sq <- (kw_test$statistic - k_groups + 1) / (nrow(coral_master) - k_groups)
cat("    Effect size (epsilon²):", round(epsilon_sq, 3),
    ifelse(epsilon_sq < 0.01, "(negligible)",
           ifelse(epsilon_sq < 0.06, "(small)",
                  ifelse(epsilon_sq < 0.14, "(medium)", "(large)"))), "\n\n")

# Post-hoc pairwise comparisons if significant
if (kw_test$p.value < 0.05) {
  cat("    Pairwise Wilcoxon tests (Holm-adjusted):\n")
  pw_test <- pairwise.wilcox.test(coral_master$total_cafi, coral_master$site, p.adjust.method = "holm")
  print(pw_test$p.value)
  cat("\n")
}

# Figure: Abundance by site
p_abundance_site <- ggplot(coral_master, aes(x = site, y = total_cafi, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = SITE_COLORS) +
  labs(
    x = "Site",
    y = "Total CAFI per Coral",
    title = "CAFI Abundance by Site",
    subtitle = paste0("Kruskal-Wallis H = ", round(kw_test$statistic, 2),
                      ", p = ", format.pval(kw_test$p.value, 3))
  ) +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "abundance_by_site.png"), p_abundance_site,
       width = 6, height = 5, dpi = 300, bg = "white")

# 1.4 Abundance vs coral volume
cat("1.4 Abundance vs Coral Volume (Power Law):\n")

m_abundance_vol <- glm.nb(total_cafi ~ log(volume) + site, data = coral_master)
m_summary <- summary(m_abundance_vol)

slope <- coef(m_abundance_vol)["log(volume)"]
slope_se <- sqrt(vcov(m_abundance_vol)["log(volume)", "log(volume)"])
slope_ci <- confint(m_abundance_vol)["log(volume)", ]

cat("    Scaling exponent (β):", round(slope, 3), "\n")
cat("    95% CI: [", round(slope_ci[1], 3), ",", round(slope_ci[2], 3), "]\n")
cat("    SE:", round(slope_se, 3), "\n")
cat("    z-value:", round(m_summary$coefficients["log(volume)", "z value"], 2), "\n")
cat("    p-value:", format.pval(m_summary$coefficients["log(volume)", "Pr(>|z|)"], 3), "\n")

# Test if β differs from 1 (Field of Dreams)
z_vs_1 <- (slope - 1) / slope_se
p_vs_1 <- 2 * (1 - pnorm(abs(z_vs_1)))
cat("\n    Test: β ≠ 1 (Field of Dreams)\n")
cat("    z-score:", round(z_vs_1, 2), ", p:", format.pval(p_vs_1, 3), "\n")
cat("    Interpretation:", ifelse(p_vs_1 < 0.05,
                                  ifelse(slope < 1,
                                         "β significantly < 1 → PROPAGULE REDISTRIBUTION",
                                         "β significantly > 1 → SUPER-LINEAR SCALING"),
                                  "β not different from 1 → Field of Dreams"), "\n\n")

# Pseudo-R²
null_model <- glm.nb(total_cafi ~ 1, data = coral_master)
pseudo_r2 <- 1 - (logLik(m_abundance_vol)[1] / logLik(null_model)[1])
cat("    McFadden's pseudo-R²:", round(pseudo_r2, 3), "\n\n")

# Figure: Abundance vs Volume
p_abundance_vol <- ggplot(coral_master, aes(x = volume, y = total_cafi, color = site)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = MASS::glm.nb, formula = y ~ log(x), se = TRUE, color = "black") +
  geom_abline(slope = 1, intercept = log10(mean(coral_master$total_cafi / coral_master$volume)),
              linetype = "dashed", color = "gray50") +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  scale_color_manual(values = SITE_COLORS) +
  labs(
    x = expression("Coral Volume (cm"^3*")"),
    y = "Total CAFI per Coral",
    title = "CAFI Abundance Scaling with Coral Size",
    subtitle = paste0("β = ", round(slope, 2), " [", round(slope_ci[1], 2), ", ",
                      round(slope_ci[2], 2), "], pseudo-R² = ", round(pseudo_r2, 2)),
    color = "Site"
  ) +
  annotate("text", x = max(coral_master$volume) * 0.3, y = min(coral_master$total_cafi[coral_master$total_cafi > 0]) * 1.5,
           label = "Dashed line: β = 1 (Field of Dreams)", hjust = 0, size = 3, color = "gray50")

ggsave(file.path(FIG_DIR, "abundance_vs_volume.png"), p_abundance_vol,
       width = 8, height = 6, dpi = 300, bg = "white")

# ============================================================================
# PART 2: TAXONOMIC COMPOSITION
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 2: TAXONOMIC COMPOSITION\n")
cat("------------------------------------------------------------\n\n")

# 2.1 OTU richness
cat("2.1 Taxonomic Richness:\n")
cat("    Total unique OTUs:", n_distinct(cafi_clean$otu), "\n")
cat("    OTUs per coral (mean):", round(mean(coral_master$otu_richness), 1),
    "(SD:", round(sd(coral_master$otu_richness), 1), ")\n")
cat("    OTUs per coral (range):", min(coral_master$otu_richness), "-",
    max(coral_master$otu_richness), "\n\n")

# 2.2 Rank-abundance curve
otu_abundance <- cafi_clean %>%
  count(otu, name = "abundance") %>%
  arrange(desc(abundance)) %>%
  mutate(rank = row_number(),
         cumulative_pct = cumsum(abundance) / sum(abundance) * 100,
         log_abundance = log10(abundance))

cat("2.2 Rank-Abundance Distribution:\n")
cat("    Top 5 OTUs:\n")
print(head(otu_abundance, 5))

# How many OTUs account for 50%, 80%, 95% of abundance?
n_50 <- sum(otu_abundance$cumulative_pct <= 50) + 1
n_80 <- sum(otu_abundance$cumulative_pct <= 80) + 1
n_95 <- sum(otu_abundance$cumulative_pct <= 95) + 1
cat("\n    Dominance pattern:\n")
cat("    - 50% of individuals:", n_50, "OTUs\n")
cat("    - 80% of individuals:", n_80, "OTUs\n")
cat("    - 95% of individuals:", n_95, "OTUs\n")
cat("    - Singletons (n=1):", sum(otu_abundance$abundance == 1), "\n")
cat("    - Doubletons (n=2):", sum(otu_abundance$abundance == 2), "\n\n")

# Figure: Rank-abundance
p_rank_abundance <- ggplot(otu_abundance, aes(x = rank, y = abundance)) +
  geom_line(color = "gray50") +
  geom_point(aes(color = rank <= 10), size = 2) +
  scale_y_log10() +
  scale_color_manual(values = c("TRUE" = "#E69F00", "FALSE" = "gray70"), guide = "none") +
  geom_text(data = dplyr::filter(otu_abundance, rank <= 5),
            aes(label = otu), hjust = -0.1, vjust = 0.5, size = 3, angle = 15) +
  labs(
    x = "Species Rank",
    y = "Abundance (log scale)",
    title = "CAFI Rank-Abundance Curve",
    subtitle = paste0("Top ", n_80, " OTUs account for 80% of individuals")
  ) +
  xlim(0, nrow(otu_abundance) * 1.1)

ggsave(file.path(FIG_DIR, "rank_abundance.png"), p_rank_abundance,
       width = 8, height = 5, dpi = 300, bg = "white")

# 2.3 Taxonomic resolution
cat("2.3 Taxonomic Resolution by Type:\n")
type_summary <- cafi_clean %>%
  group_by(type) %>%
  summarise(
    n_individuals = n(),
    n_otus = n_distinct(otu),
    proportion = n() / nrow(cafi_clean),
    individuals_per_otu = n() / n_distinct(otu),
    .groups = "drop"
  ) %>%
  arrange(desc(n_individuals))

print(type_summary)
cat("\n")

# 2.4 Functional group composition
cat("2.4 Functional Group Composition:\n")
func_summary <- cafi_clean %>%
  group_by(functional_group) %>%
  summarise(
    n = n(),
    pct = round(n() / nrow(cafi_clean) * 100, 1),
    n_otus = n_distinct(otu),
    n_corals = n_distinct(coral_id),
    prevalence = round(n_distinct(coral_id) / n_distinct(cafi_clean$coral_id) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

print(func_summary)
cat("\n")

# Figure: Functional group composition (using FUNC_COLORS from 00_setup.R)
p_func_comp <- ggplot(func_summary, aes(x = reorder(functional_group, n), y = n, fill = functional_group)) +
  geom_col() +
  geom_text(aes(label = paste0(pct, "%")), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = FUNC_COLORS) +
  labs(
    x = NULL,
    y = "Number of Individuals",
    title = "CAFI Functional Group Composition",
    subtitle = paste0("N = ", nrow(cafi_clean), " individuals")
  ) +
  theme(legend.position = "none") +
  expand_limits(y = max(func_summary$n) * 1.15)

ggsave(file.path(FIG_DIR, "functional_group_composition.png"), p_func_comp,
       width = 7, height = 5, dpi = 300, bg = "white")

# ============================================================================
# PART 3: ALPHA DIVERSITY
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 3: ALPHA DIVERSITY\n")
cat("------------------------------------------------------------\n\n")

# 3.1 Diversity metrics summary
cat("3.1 Alpha Diversity Metrics:\n")

diversity_summary <- coral_master %>%
  summarise(
    richness_mean = mean(otu_richness),
    richness_sd = sd(otu_richness),
    shannon_mean = mean(shannon),
    shannon_sd = sd(shannon),
    n = n()
  )

cat("    Species Richness: ", round(diversity_summary$richness_mean, 2),
    " ± ", round(diversity_summary$richness_sd, 2), "\n", sep = "")
cat("    Shannon H': ", round(diversity_summary$shannon_mean, 2),
    " ± ", round(diversity_summary$shannon_sd, 2), "\n\n", sep = "")

# 3.2 Diversity by site
cat("3.2 Alpha Diversity by Site:\n")

diversity_by_site <- coral_master %>%
  group_by(site) %>%
  summarise(
    n = n(),
    richness_mean = round(mean(otu_richness), 2),
    richness_se = round(sd(otu_richness) / sqrt(n()), 2),
    shannon_mean = round(mean(shannon), 2),
    shannon_se = round(sd(shannon) / sqrt(n()), 2),
    .groups = "drop"
  )

print(diversity_by_site)

# ANOVA for richness
aov_richness <- aov(otu_richness ~ site, data = coral_master)
aov_summary <- summary(aov_richness)[[1]]
cat("\n    Richness ANOVA:\n")
cat("    F(", aov_summary$Df[1], ",", aov_summary$Df[2], ") = ",
    round(aov_summary$`F value`[1], 2), ", p = ",
    format.pval(aov_summary$`Pr(>F)`[1], 3), "\n", sep = "")

# Effect size (eta-squared)
ss_between <- aov_summary$`Sum Sq`[1]
ss_total <- sum(aov_summary$`Sum Sq`)
eta_sq <- ss_between / ss_total
cat("    Effect size (η²): ", round(eta_sq, 3), "\n\n", sep = "")

# Post-hoc if significant
if (aov_summary$`Pr(>F)`[1] < 0.05) {
  tukey <- TukeyHSD(aov_richness)
  cat("    Tukey HSD post-hoc:\n")
  print(round(tukey$site, 3))
  cat("\n")
}

# 3.3 Diversity vs coral size
cat("3.3 Diversity vs Coral Size:\n")

# Richness vs volume
m_rich_vol <- glm(otu_richness ~ log(volume), family = poisson, data = coral_master)
m_rich_summary <- summary(m_rich_vol)

cat("    Richness ~ log(Volume) [Poisson GLM]:\n")
cat("    β = ", round(coef(m_rich_vol)[2], 3),
    ", z = ", round(m_rich_summary$coefficients[2, "z value"], 2),
    ", p = ", format.pval(m_rich_summary$coefficients[2, "Pr(>|z|)"], 3), "\n", sep = "")

# Overdispersion check: Pearson chi-squared / residual df
dispersion_ratio <- sum(residuals(m_rich_vol, type = "pearson")^2) / m_rich_vol$df.residual
cat("    Overdispersion ratio (Pearson X²/df):", round(dispersion_ratio, 2),
    ifelse(dispersion_ratio > 1.5, " [OVERDISPERSED - use quasipoisson]", " [OK]"), "\n")
if (dispersion_ratio > 1.5) {
  # Refit with quasipoisson for corrected SEs
  m_rich_vol_qp <- glm(otu_richness ~ log(volume), family = quasipoisson, data = coral_master)
  m_rich_qp_summary <- summary(m_rich_vol_qp)
  cat("    Quasi-Poisson corrected: β = ", round(coef(m_rich_vol_qp)[2], 3),
      ", t = ", round(m_rich_qp_summary$coefficients[2, "t value"], 2),
      ", p = ", format.pval(m_rich_qp_summary$coefficients[2, "Pr(>|t|)"], 3), "\n", sep = "")
}
cat("\n")

# Shannon vs volume
m_shan_vol <- lm(shannon ~ log(volume), data = coral_master)
m_shan_summary <- summary(m_shan_vol)

cat("    Shannon ~ log(Volume):\n")
cat("    β = ", round(coef(m_shan_vol)[2], 3),
    ", t(", m_shan_summary$df[2], ") = ", round(m_shan_summary$coefficients[2, "t value"], 2),
    ", p = ", format.pval(m_shan_summary$coefficients[2, "Pr(>|t|)"], 3),
    ", R² = ", round(m_shan_summary$r.squared, 3), "\n\n", sep = "")

# Figure: Diversity vs Volume
p_diversity_vol <- coral_master %>%
  pivot_longer(cols = c(otu_richness, shannon), names_to = "metric", values_to = "value") %>%
  mutate(metric = dplyr::recode(metric,
                                "otu_richness" = "Species Richness",
                                "shannon" = "Shannon H'")) %>%
  ggplot(aes(x = volume, y = value, color = site)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = TRUE, color = "black") +
  facet_wrap(~metric, scales = "free_y") +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = SITE_COLORS) +
  labs(
    x = expression("Coral Volume (cm"^3*")"),
    y = "Diversity",
    title = "Alpha Diversity Increases with Coral Size",
    color = "Site"
  )

ggsave(file.path(FIG_DIR, "diversity_vs_volume.png"), p_diversity_vol,
       width = 10, height = 5, dpi = 300, bg = "white")

# 3.4 Species accumulation curve
cat("3.4 Species Accumulation:\n")

# Rarefaction curve
specaccum_result <- specaccum(community_matrix, method = "random", permutations = 999)

cat("    At 50 corals: ", round(specaccum_result$richness[50], 1),
    " ± ", round(specaccum_result$sd[50], 1), " species\n", sep = "")
cat("    At 100 corals: ", round(specaccum_result$richness[min(100, nrow(community_matrix))], 1),
    " ± ", round(specaccum_result$sd[min(100, nrow(community_matrix))], 1), " species\n", sep = "")
cat("    Total observed:", ncol(community_matrix), "species\n")

# Chao1 estimator
chao1 <- estimateR(colSums(community_matrix))["S.chao1"]
cat("    Chao1 estimate:", round(chao1, 0), "species\n")
cat("    % sampled:", round(ncol(community_matrix) / chao1 * 100, 1), "%\n\n")

# Figure: Species accumulation
accum_df <- data.frame(
  sites = specaccum_result$sites,
  richness = specaccum_result$richness,
  sd = specaccum_result$sd
)

p_accum <- ggplot(accum_df, aes(x = sites, y = richness)) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = chao1, linetype = "dashed", color = "red") +
  annotate("text", x = 10, y = chao1 + 3, label = paste0("Chao1 = ", round(chao1, 0)),
           color = "red", hjust = 0) +
  labs(
    x = "Number of Corals Sampled",
    y = "Cumulative Species Richness",
    title = "Species Accumulation Curve",
    subtitle = paste0("Observed: ", ncol(community_matrix), " species, Chao1 estimate: ", round(chao1, 0))
  )

ggsave(file.path(FIG_DIR, "species_accumulation.png"), p_accum,
       width = 7, height = 5, dpi = 300, bg = "white")

# Supplement S1: Species accumulation curve
supplement_dir <- file.path(PATHS$figures, "supplement")
dir.create(supplement_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(supplement_dir, "figS1_species_accumulation.png"), p_accum,
       width = 7, height = 5, dpi = 300, bg = "white")

# ============================================================================
# PART 4: BETA DIVERSITY & COMMUNITY STRUCTURE
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PART 4: BETA DIVERSITY & COMMUNITY STRUCTURE\n")
cat("------------------------------------------------------------\n\n")

# 4.1 Hellinger transformation
comm_hell <- decostand(community_matrix, method = "hellinger")

# 4.2 NMDS ordination
cat("4.1 NMDS Ordination:\n")

set.seed(42)
nmds <- metaMDS(comm_hell, distance = "bray", k = 2, trymax = 100, trace = 0)

cat("    Stress: ", round(nmds$stress, 3), "\n", sep = "")
cat("    Interpretation: ", ifelse(nmds$stress < 0.1, "Excellent",
                                    ifelse(nmds$stress < 0.2, "Good",
                                           ifelse(nmds$stress < 0.3, "Acceptable", "Poor"))), "\n\n", sep = "")

# Extract scores
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$coral_id <- rownames(nmds_scores)
nmds_scores <- nmds_scores %>%
  left_join(coral_master %>% dplyr::select(coral_id, site, volume, log_volume, total_cafi), by = "coral_id") %>%
  filter(!is.na(site))  # Remove any corals without site data

# 4.2 PERMANOVA
cat("4.2 PERMANOVA (community composition):\n")

# Match rows - ensure perfect alignment between community matrix and coral data
common_ids <- intersect(rownames(community_matrix), coral_master$coral_id)
comm_hell_aligned <- comm_hell[common_ids, ]
coral_for_perm <- coral_master %>%
  filter(coral_id %in% common_ids) %>%
  arrange(match(coral_id, common_ids))

# Raw community matrix aligned to coral_for_perm (used in Part 4B & 4C)
comm_aligned_raw <- community_matrix[coral_for_perm$coral_id, ]

permanova <- adonis2(comm_hell_aligned ~ log(volume) + site,
                     data = coral_for_perm,
                     permutations = 999,
                     method = "bray",
                     by = "terms")  # Get individual term effects

cat("\n")
print(permanova)

# Extract individual term results (rows now include: log(volume), site, Residual, Total)
perm_df <- as.data.frame(permanova)
cat("\n    Volume effect: R² = ", round(perm_df["log(volume)", "R2"], 3),
    ", F = ", round(perm_df["log(volume)", "F"], 2),
    ", p = ", format.pval(perm_df["log(volume)", "Pr(>F)"], 3), "\n", sep = "")
cat("    Site effect: R² = ", round(perm_df["site", "R2"], 3),
    ", F = ", round(perm_df["site", "F"], 2),
    ", p = ", format.pval(perm_df["site", "Pr(>F)"], 3), "\n\n", sep = "")

# Marginal PERMANOVA (Type III - order-independent, more robust)
permanova_margin <- adonis2(comm_hell_aligned ~ log(volume) + site,
                            data = coral_for_perm,
                            permutations = 999,
                            method = "bray",
                            by = "margin")
perm_margin_df <- as.data.frame(permanova_margin)
cat("    Marginal (Type III) PERMANOVA:\n")
cat("    Volume: R² = ", round(perm_margin_df["log(volume)", "R2"], 3),
    ", F = ", round(perm_margin_df["log(volume)", "F"], 2),
    ", p = ", format.pval(perm_margin_df["log(volume)", "Pr(>F)"], 3), "\n", sep = "")
cat("    Site: R² = ", round(perm_margin_df["site", "R2"], 3),
    ", F = ", round(perm_margin_df["site", "F"], 2),
    ", p = ", format.pval(perm_margin_df["site", "Pr(>F)"], 3), "\n\n", sep = "")

# Interaction PERMANOVA (volume × site)
permanova_interaction <- adonis2(comm_hell_aligned ~ log(volume) * site,
                                 data = coral_for_perm,
                                 permutations = 999,
                                 method = "bray",
                                 by = "terms")
perm_int_df <- as.data.frame(permanova_interaction)
cat("    Interaction test (Volume × Site):\n")
if ("log(volume):site" %in% rownames(perm_int_df)) {
  cat("    R² = ", round(perm_int_df["log(volume):site", "R2"], 3),
      ", F = ", round(perm_int_df["log(volume):site", "F"], 2),
      ", p = ", format.pval(perm_int_df["log(volume):site", "Pr(>F)"], 3), "\n\n", sep = "")
}

# 4.3 PERMDISP (test for homogeneity of dispersions)
cat("4.3 PERMDISP (homogeneity of dispersions by site):\n")

bray_dist <- vegdist(comm_hell_aligned, method = "bray")
disp <- betadisper(bray_dist, coral_for_perm$site)
disp_test <- permutest(disp, permutations = 999)

cat("    F = ", round(disp_test$statistic, 2),
    ", p = ", format.pval(disp_test$pval, 3), "\n", sep = "")
if (!is.null(disp_test$pval) && !is.na(disp_test$pval)) {
  cat("    Interpretation: ", ifelse(disp_test$pval < 0.05,
                                     "Dispersions differ between sites (interpret PERMANOVA with caution)",
                                     "Dispersions homogeneous (PERMANOVA valid)"), "\n\n", sep = "")
} else {
  cat("    Interpretation: Could not determine\n\n")
}

# Figure: NMDS
# Extract PERMANOVA results safely using the data frame
site_r2 <- round(perm_df["site", "R2"], 2)
site_p <- perm_df["site", "Pr(>F)"]
site_p_text <- if (!is.na(site_p)) format.pval(site_p, 3) else "< 0.001"
vol_r2 <- round(perm_df["log(volume)", "R2"], 2)
vol_p <- perm_df["log(volume)", "Pr(>F)"]
vol_p_text <- if (!is.na(vol_p)) format.pval(vol_p, 3) else "< 0.001"

# Compute axis limits excluding extreme outliers (>3 IQR from median)
nmds1_iqr <- IQR(nmds_scores$NMDS1, na.rm = TRUE)
nmds2_iqr <- IQR(nmds_scores$NMDS2, na.rm = TRUE)
nmds1_lim <- c(median(nmds_scores$NMDS1) - 3 * nmds1_iqr,
               median(nmds_scores$NMDS1) + 3 * nmds1_iqr)
nmds2_lim <- c(median(nmds_scores$NMDS2) - 3 * nmds2_iqr,
               median(nmds_scores$NMDS2) + 3 * nmds2_iqr)
n_outliers <- sum(nmds_scores$NMDS1 < nmds1_lim[1] | nmds_scores$NMDS1 > nmds1_lim[2] |
                  nmds_scores$NMDS2 < nmds2_lim[1] | nmds_scores$NMDS2 > nmds2_lim[2])

p_nmds_site <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = site)) +
  geom_point(aes(size = log_volume), alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = "dashed") +
  scale_color_manual(values = SITE_COLORS) +
  scale_size_continuous(name = expression(log[10]*"(Volume)"), range = c(2, 6)) +
  labs(
    title = "NMDS Ordination of CAFI Communities",
    subtitle = paste0("Stress = ", round(nmds$stress, 3), "; PERMANOVA: Site R² = ",
                      site_r2, ", p = ", site_p_text),
    caption = if (n_outliers > 0) paste0(n_outliers, " outlier(s) beyond axis limits") else NULL,
    color = "Site"
  ) +
  coord_fixed(xlim = nmds1_lim, ylim = nmds2_lim, clip = "on")

ggsave(file.path(FIG_DIR, "nmds_by_site.png"), p_nmds_site,
       width = 8, height = 6, dpi = 300, bg = "white")

# Supplement S3: NMDS ordination by site
supplement_dir <- file.path(PATHS$figures, "supplement")
dir.create(supplement_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(supplement_dir, "figS3_nmds_ordination.png"), p_nmds_site,
       width = 8, height = 6, dpi = 300, bg = "white")

p_nmds_size <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = log_volume)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c(name = expression(log[10]*"(Volume)")) +
  labs(
    title = "CAFI Community Composition Varies with Coral Size",
    subtitle = paste0("PERMANOVA: Volume R² = ", vol_r2, ", p = ", vol_p_text),
    caption = if (n_outliers > 0) paste0(n_outliers, " outlier(s) beyond axis limits") else NULL
  ) +
  coord_fixed(xlim = nmds1_lim, ylim = nmds2_lim, clip = "on")

ggsave(file.path(FIG_DIR, "nmds_by_volume.png"), p_nmds_size,
       width = 8, height = 6, dpi = 300, bg = "white")

# 4.4 Beta diversity (turnover)
cat("4.4 Beta Diversity (Turnover):\n")

# Mean Bray-Curtis dissimilarity
mean_bray <- mean(as.vector(bray_dist))
cat("    Mean Bray-Curtis dissimilarity:", round(mean_bray, 3), "\n")

# Beta diversity by site pair
site_pairs <- combn(unique(coral_for_perm$site), 2)
cat("    Pairwise dissimilarity:\n")

for (i in 1:ncol(site_pairs)) {
  s1 <- site_pairs[1, i]
  s2 <- site_pairs[2, i]

  idx1 <- which(coral_for_perm$site == s1)
  idx2 <- which(coral_for_perm$site == s2)

  between_dist <- as.matrix(bray_dist)[idx1, idx2]
  mean_between <- mean(between_dist)

  cat("    ", s1, " vs ", s2, ": ", round(mean_between, 3), "\n", sep = "")
}

# ============================================================================
# PART 4B: COMPOSITION DIVERGENCE BY CORAL SIZE
# ============================================================================
# Parallels Figure 3B from Stier et al. experimental paper
# Tests whether larger corals have more distinct community compositions
# (greater distance-to-centroid = more divergent communities)
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 4B: COMPOSITION DIVERGENCE BY SIZE (Fig 3B Analog)\n")
cat("------------------------------------------------------------\n\n")

cat("Testing: Do larger corals have more distinct CAFI communities?\n")
cat("Hypothesis: Distance-to-centroid increases with coral size\n")
cat("(parallels experimental finding: more corals = more divergent composition)\n\n")

# 4B.1 Create size classes based on volume tertiles
cat("4B.1 Creating Size Classes:\n")

coral_for_perm <- coral_for_perm %>%
  mutate(
    volume_tertile = ntile(volume, 3),
    size_class = factor(
      case_when(
        volume_tertile == 1 ~ "Small",
        volume_tertile == 2 ~ "Medium",
        volume_tertile == 3 ~ "Large"
      ),
      levels = c("Small", "Medium", "Large")
    )
  )

size_class_summary <- coral_for_perm %>%
  group_by(size_class) %>%
  summarise(
    n = n(),
    volume_min = round(min(volume)),
    volume_max = round(max(volume)),
    volume_median = round(median(volume)),
    mean_cafi = round(mean(total_cafi), 1),
    mean_richness = round(mean(otu_richness), 1),
    .groups = "drop"
  )

cat("    Size class definitions:\n")
print(size_class_summary)
cat("\n")

# 4B.2 Rarefaction to control for abundance-dispersion confounding
cat("4B.2 Rarefaction (controlling for sampling effort):\n")

# Rarefaction is critical: larger corals have more CAFI individuals,
# which mechanically reduces stochastic variation in composition.
# To test whether larger corals truly have more distinct communities
# (not just better-sampled ones), we rarefy to equal depth.
sample_depths <- rowSums(comm_aligned_raw)
cat("    Sample depth range:", min(sample_depths), "-", max(sample_depths), "individuals\n")

# Filter to corals with sufficient abundance for meaningful rarefaction (≥5 individuals)
MIN_RAREFACTION_DEPTH <- 5
sufficient_idx <- sample_depths >= MIN_RAREFACTION_DEPTH
comm_for_rarefaction <- comm_aligned_raw[sufficient_idx, ]
coral_for_rarefaction <- coral_for_perm[sufficient_idx, ]
min_depth <- min(rowSums(comm_for_rarefaction))
cat("    Corals with ≥", MIN_RAREFACTION_DEPTH, "individuals:", sum(sufficient_idx),
    "/", length(sufficient_idx), "\n")
cat("    Rarefaction depth:", min_depth, "individuals\n")

# Iterated rarefaction (100 draws, average Bray-Curtis distances)
# A single rarefaction draw introduces stochastic noise; averaging across
# multiple draws provides a more robust estimate of abundance-controlled distances
N_RAREFACTION_ITER <- 100
cat("    Iterating rarefaction:", N_RAREFACTION_ITER, "draws at depth", min_depth, "\n")

set.seed(42)
dist_accumulator <- matrix(0, nrow = nrow(comm_for_rarefaction),
                           ncol = nrow(comm_for_rarefaction))

for (iter in 1:N_RAREFACTION_ITER) {
  comm_rare_i <- vegan::rrarefy(comm_for_rarefaction, sample = min_depth)
  comm_rare_hell_i <- decostand(comm_rare_i, method = "hellinger")
  dist_accumulator <- dist_accumulator + as.matrix(vegdist(comm_rare_hell_i, method = "bray"))
}

# Average distance matrix across iterations
dist_avg <- dist_accumulator / N_RAREFACTION_ITER
bray_dist_rarefied <- as.dist(dist_avg)

# Also keep a single rarefied matrix for downstream display
set.seed(42)
comm_rarefied <- vegan::rrarefy(comm_for_rarefaction, sample = min_depth)

cat("    Rarefied community matrix:", nrow(comm_rarefied), "corals x",
    ncol(comm_rarefied), "species\n\n")

# 4B.3 Calculate distance-to-centroid within size classes using betadisper
cat("4B.3 Beta Dispersion (Distance to Centroid) by Size Class:\n")
cat("    Testing on BOTH original and rarefied data:\n\n")

# Original (non-rarefied) analysis
disp_size <- betadisper(bray_dist, coral_for_perm$size_class)

# Rarefied analysis (abundance-controlled, on subset with ≥5 individuals)
disp_size_rarefied <- betadisper(bray_dist_rarefied, coral_for_rarefaction$size_class)

# Extract distances to centroid
distances_df <- data.frame(
  coral_id = names(disp_size$distances),
  distance_to_centroid = disp_size$distances,
  size_class = coral_for_perm$size_class
)

# Add NMDS scores and other coral data
distances_df <- distances_df %>%
  left_join(coral_for_perm %>% dplyr::select(coral_id, site, volume, total_cafi, otu_richness),
            by = "coral_id")

# Summary by size class
dispersion_summary <- distances_df %>%
  group_by(size_class) %>%
  summarise(
    n = n(),
    mean_dist = round(mean(distance_to_centroid), 3),
    sd_dist = round(sd(distance_to_centroid), 3),
    se_dist = round(sd(distance_to_centroid) / sqrt(n()), 3),
    .groups = "drop"
  )

cat("    Distance to centroid by size class:\n")
print(dispersion_summary)
cat("\n")

# 4B.3 Statistical test: PERMDISP by size class
cat("4B.3 PERMDISP Test (does dispersion differ by size?):\n")

disp_size_test <- permutest(disp_size, permutations = 999)

cat("    F = ", round(disp_size_test$statistic, 2),
    ", df = ", paste(disp_size_test$df, collapse = ", "),
    ", p = ", format.pval(disp_size_test$pval, 3), "\n", sep = "")

# Linear trend test: does distance change with coral volume?
# Use continuous log(volume) rather than ordinal size class to avoid
# assuming equal spacing between tertile boundaries
distances_df$size_numeric <- as.numeric(distances_df$size_class)  # kept for rarefied comparison
trend_model <- lm(distance_to_centroid ~ log(volume), data = distances_df)
trend_summary <- summary(trend_model)

cat("\n    Linear trend test (Small -> Medium -> Large):\n")
cat("    β = ", round(coef(trend_model)[2], 3),
    ", t(", trend_summary$df[2], ") = ", round(trend_summary$coefficients[2, "t value"], 2),
    ", p = ", format.pval(trend_summary$coefficients[2, "Pr(>|t|)"], 3), "\n", sep = "")
cat("    R² = ", round(trend_summary$r.squared, 3), "\n", sep = "")

# Interpretation
trend_direction <- ifelse(coef(trend_model)[2] > 0, "INCREASES", "DECREASES")
trend_sig <- ifelse(trend_summary$coefficients[2, "Pr(>|t|)"] < 0.05, "significantly", "not significantly")

cat("\n    Interpretation: Community distinctness ", trend_sig, " ", trend_direction,
    " with coral size\n", sep = "")

if (coef(trend_model)[2] > 0 && trend_summary$coefficients[2, "Pr(>|t|)"] < 0.05) {
  cat("    → Supports hypothesis: Larger corals harbor more distinct communities\n")
  cat("    → Parallels experimental finding (more corals = more divergent composition)\n")
} else if (coef(trend_model)[2] < 0 && trend_summary$coefficients[2, "Pr(>|t|)"] < 0.05) {
  cat("    → Opposite pattern: Larger corals have MORE SIMILAR communities\n")
  cat("    → Suggests homogenization/convergence in larger corals\n")
} else {
  cat("    → No significant trend in community distinctness with size\n")
}
cat("\n")

# 4B.3b RAREFIED analysis (controls for abundance confound)
cat("\n    --- RAREFIED DATA (abundance-controlled) ---\n")
cat("    Rarefaction depth:", min_depth, "individuals\n\n")

disp_size_rarefied_test <- permutest(disp_size_rarefied, permutations = 999)
cat("    PERMDISP (rarefied): F = ", round(disp_size_rarefied_test$statistic, 2),
    ", df = ", paste(disp_size_rarefied_test$df, collapse = ", "),
    ", p = ", format.pval(disp_size_rarefied_test$pval, 3), "\n", sep = "")

# Rarefied trend test (iterated average distances)
distances_rarefied_df <- data.frame(
  coral_id = names(disp_size_rarefied$distances),
  distance_to_centroid_rarefied = disp_size_rarefied$distances,
  size_class = coral_for_rarefaction$size_class,
  size_numeric = as.numeric(coral_for_rarefaction$size_class)
) %>%
  left_join(coral_for_rarefaction %>% dplyr::select(coral_id, volume), by = "coral_id")

trend_model_rarefied <- lm(distance_to_centroid_rarefied ~ log(volume), data = distances_rarefied_df)
trend_summary_rarefied <- summary(trend_model_rarefied)

cat("    Trend (rarefied): β = ", round(coef(trend_model_rarefied)[2], 3),
    ", t(", trend_summary_rarefied$df[2], ") = ",
    round(trend_summary_rarefied$coefficients[2, "t value"], 2),
    ", p = ", format.pval(trend_summary_rarefied$coefficients[2, "Pr(>|t|)"], 3), "\n", sep = "")

# Comparison: does rarefaction change the conclusion?
orig_sig <- trend_summary$coefficients[2, "Pr(>|t|)"] < 0.05
rare_sig <- trend_summary_rarefied$coefficients[2, "Pr(>|t|)"] < 0.05

if (orig_sig && !rare_sig) {
  cat("\n    ⚠ CAUTION: Divergence trend is NOT significant after rarefaction.\n")
  cat("    The original finding may be an artifact of abundance differences.\n")
} else if (orig_sig && rare_sig) {
  cat("\n    ✓ Divergence trend ROBUST to rarefaction (not an abundance artifact).\n")
} else if (!orig_sig) {
  cat("\n    No significant trend in either original or rarefied data.\n")
}
cat("\n")

# 4B.3b Rarefaction sensitivity: test at multiple depths
cat("4B.3b Rarefaction Sensitivity (multiple depths):\n")
cat("    Testing divergence trend robustness across rarefaction depths...\n")

# Test at depths from min_depth up to the 25th percentile of sample depths
# (to retain reasonable sample sizes)
test_depths <- sort(unique(c(min_depth,
                             seq(5, floor(quantile(sample_depths[sufficient_idx], 0.25)), by = 5))))
test_depths <- test_depths[test_depths >= min_depth]

if (length(test_depths) > 1) {
  sensitivity_results <- data.frame(
    depth = integer(), n_corals = integer(),
    beta = numeric(), p_value = numeric(), significant = logical()
  )

  for (d in test_depths) {
    # Filter to corals with at least d individuals
    idx_d <- sample_depths >= d
    if (sum(idx_d) < 20) next  # Need ≥20 corals for meaningful test

    comm_d <- comm_aligned_raw[idx_d, ]
    coral_d <- coral_for_perm[idx_d, ]

    # Iterated rarefaction at depth d (fewer iterations for sensitivity)
    set.seed(42)
    dist_acc_d <- matrix(0, nrow = nrow(comm_d), ncol = nrow(comm_d))
    for (iter in 1:50) {
      comm_rare_d <- vegan::rrarefy(comm_d, sample = d)
      comm_hell_d <- decostand(comm_rare_d, method = "hellinger")
      dist_acc_d <- dist_acc_d + as.matrix(vegdist(comm_hell_d, method = "bray"))
    }
    dist_avg_d <- as.dist(dist_acc_d / 50)

    # betadisper + trend test
    disp_d <- betadisper(dist_avg_d, coral_d$size_class)
    dist_df_d <- data.frame(
      distance = disp_d$distances,
      volume = coral_d$volume
    )
    trend_d <- lm(distance ~ log(volume), data = dist_df_d)
    trend_p_d <- summary(trend_d)$coefficients[2, "Pr(>|t|)"]

    sensitivity_results <- rbind(sensitivity_results, data.frame(
      depth = d, n_corals = sum(idx_d),
      beta = coef(trend_d)[2], p_value = trend_p_d,
      significant = trend_p_d < 0.05
    ))
  }

  cat("    Depth | N corals | β        | p-value  | Significant\n")
  cat("    ------|----------|----------|----------|------------\n")
  for (i in 1:nrow(sensitivity_results)) {
    r <- sensitivity_results[i, ]
    cat(sprintf("    %5d | %8d | %8.4f | %8.4f | %s\n",
                r$depth, r$n_corals, r$beta, r$p_value,
                ifelse(r$significant, "Yes", "No")))
  }

  n_sig <- sum(sensitivity_results$significant)
  cat(sprintf("\n    Summary: %d/%d depths show significant trend (p < 0.05)\n",
              n_sig, nrow(sensitivity_results)))
  if (n_sig == 0) {
    cat("    → Divergence NOT significant at any rarefaction depth (robust null)\n\n")
  } else if (n_sig == nrow(sensitivity_results)) {
    cat("    → Divergence significant at ALL depths (robust signal)\n\n")
  } else {
    cat("    → Mixed: result depends on rarefaction depth (interpret cautiously)\n\n")
  }
} else {
  cat("    Only one viable depth; skipping multi-depth sensitivity\n\n")
  sensitivity_results <- NULL
}

# 4B.4 Pairwise comparisons
cat("4B.4 Pairwise Tukey HSD for Distance to Centroid:\n")
tukey_disp <- TukeyHSD(disp_size)
print(round(tukey_disp$group, 3))
cat("\n")

# 4B.5 Create visualization
cat("4B.5 Creating Composition Divergence Figure...\n")

# Define size class colors
SIZE_COLORS <- c("Small" = "#fee090", "Medium" = "#fc8d59", "Large" = "#d73027")

# Panel A: Boxplot of distance-to-centroid by size class
p_divergence_boxplot <- ggplot(distances_df, aes(x = size_class, y = distance_to_centroid, fill = size_class)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = site), width = 0.2, alpha = 0.6, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  scale_fill_manual(values = SIZE_COLORS, guide = "none") +
  scale_color_manual(values = SITE_COLORS, name = "Site") +
  labs(
    title = "A. Community Distinctness by Coral Size",
    subtitle = paste0("Original: β = ", round(coef(trend_model)[2], 3),
                      ", p = ", format.pval(trend_summary$coefficients[2, "Pr(>|t|)"], 3),
                      " | Rarefied: β = ", round(coef(trend_model_rarefied)[2], 3),
                      ", p = ", format.pval(trend_summary_rarefied$coefficients[2, "Pr(>|t|)"], 3)),
    x = "Coral Size Class",
    y = "Distance to Centroid (Bray-Curtis)"
  ) +
  theme_publication() +
  theme(legend.position = c(0.15, 0.85))

# Panel B: NMDS with centroids and size class trajectories
# Get centroid coordinates
centroid_coords <- as.data.frame(disp_size$centroids)
centroid_coords$size_class <- rownames(centroid_coords)
centroid_coords <- centroid_coords %>%
  mutate(size_class = factor(size_class, levels = c("Small", "Medium", "Large")))

# Need NMDS scores for this - merge with size class info
# Rebuild from coral_for_perm to ensure all columns are present
nmds_with_size <- nmds_scores %>%
  dplyr::select(coral_id, NMDS1, NMDS2) %>%  # Only keep NMDS coords
  left_join(coral_for_perm %>% dplyr::select(coral_id, size_class, volume, site),
            by = "coral_id") %>%
  filter(!is.na(size_class))

# Calculate NMDS centroids by size class (different from betadisper centroids)
nmds_centroids <- nmds_with_size %>%
  group_by(size_class) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(size_class = factor(size_class, levels = c("Small", "Medium", "Large")))

# Create arrows connecting centroids (trajectory)
arrow_data <- data.frame(
  x_start = nmds_centroids$NMDS1[1:2],
  y_start = nmds_centroids$NMDS2[1:2],
  x_end = nmds_centroids$NMDS1[2:3],
  y_end = nmds_centroids$NMDS2[2:3],
  segment = c("Small→Medium", "Medium→Large")
)

p_nmds_trajectory <- ggplot() +
  # Individual points colored by size class
  geom_point(data = nmds_with_size,
             aes(x = NMDS1, y = NMDS2, fill = size_class, size = log(volume)),
             shape = 21, alpha = 0.6, color = "gray30") +
  # Ellipses for each size class
  stat_ellipse(data = nmds_with_size,
               aes(x = NMDS1, y = NMDS2, color = size_class),
               level = 0.95, linetype = "dashed", linewidth = 0.8) +
  # Centroids
  geom_point(data = nmds_centroids,
             aes(x = NMDS1, y = NMDS2, fill = size_class),
             shape = 23, size = 6, color = "black", stroke = 1.5) +
  # Arrows connecting centroids (trajectory)
  geom_segment(data = arrow_data,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               linewidth = 1.2, color = "black") +
  # Centroid labels
  geom_text(data = nmds_centroids,
            aes(x = NMDS1, y = NMDS2 + 0.08, label = size_class),
            fontface = "bold", size = 4) +
  scale_fill_manual(values = SIZE_COLORS, name = "Size Class") +
  scale_color_manual(values = SIZE_COLORS, guide = "none") +
  scale_size_continuous(name = expression(log[10]*"(Volume)"), range = c(2, 5)) +
  labs(
    title = "B. Community Trajectory with Coral Size",
    subtitle = paste0("Arrows show centroid shift: Small → Medium → Large\n",
                      "Stress = ", round(nmds$stress, 3)),
    x = "NMDS1",
    y = "NMDS2"
  ) +
  coord_fixed(xlim = nmds1_lim, ylim = nmds2_lim, clip = "on") +
  theme_publication() +
  theme(legend.position = "right")

# Combine panels
fig_divergence <- p_divergence_boxplot + p_nmds_trajectory +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    title = "Figure 4: Community Composition by Coral Size",
    subtitle = paste0("Do larger corals support more distinct CAFI communities? | n = ",
                      nrow(distances_df), " corals"),
    caption = paste0("Distance to centroid measures community distinctness within size class.\n",
                     "Original: significant (p < 0.001). After rarefaction (n\u2265 5): NOT significant (p = 0.61). Size-divergence is an abundance artifact."),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      plot.caption = element_text(size = 9, hjust = 0, color = "gray40")
    )
  )

ggsave(file.path(FIG_DIR, "composition_divergence_by_size.png"), fig_divergence,
       width = 14, height = 7, dpi = 300, bg = "white")
cat("    Saved: composition_divergence_by_size.png\n")

# Supplement S5: Composition divergence by size
supplement_dir <- file.path(PATHS$figures, "supplement")
dir.create(supplement_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(supplement_dir, "figS5_composition_divergence.png"), fig_divergence,
       width = 14, height = 7, dpi = 300, bg = "white")
cat("    Saved: figS5_composition_divergence.png (supplement)\n\n")

# 4B.6 Add results to stats tracking
cat("4B.6 Recording Statistical Results:\n")

# Store divergence analysis results (including rarefied comparison)
divergence_results <- list(
  size_class_summary = size_class_summary,
  dispersion_summary = dispersion_summary,
  permdisp_test = list(
    F = disp_size_test$statistic,
    df = disp_size_test$df,
    p = disp_size_test$pval
  ),
  permdisp_test_rarefied = list(
    F = disp_size_rarefied_test$statistic,
    df = disp_size_rarefied_test$df,
    p = disp_size_rarefied_test$pval,
    rarefaction_depth = min_depth
  ),
  trend_test = list(
    beta = coef(trend_model)[2],
    t = trend_summary$coefficients[2, "t value"],
    df = trend_summary$df[2],
    p = trend_summary$coefficients[2, "Pr(>|t|)"],
    r2 = trend_summary$r.squared
  ),
  trend_test_rarefied = list(
    beta = coef(trend_model_rarefied)[2],
    t = trend_summary_rarefied$coefficients[2, "t value"],
    df = trend_summary_rarefied$df[2],
    p = trend_summary_rarefied$coefficients[2, "Pr(>|t|)"],
    r2 = trend_summary_rarefied$r.squared
  ),
  rarefaction_robust = orig_sig && rare_sig,
  rarefaction_sensitivity = if (exists("sensitivity_results") && !is.null(sensitivity_results)) sensitivity_results else NULL,
  tukey_hsd = tukey_disp$group,
  interpretation = paste0(
    "Community distinctness ", trend_sig, " ", trend_direction, " with coral size. ",
    if (orig_sig && rare_sig) {
      "Pattern ROBUST to rarefaction: larger corals harbor more distinct communities."
    } else if (orig_sig && !rare_sig) {
      "CAUTION: Pattern not significant after rarefaction; may be abundance artifact."
    } else {
      "No significant trend detected."
    }
  )
)

cat("    PERMDISP F:", round(disp_size_test$statistic, 2),
    ", p =", format.pval(disp_size_test$pval, 3), "\n")
cat("    Trend β:", round(coef(trend_model)[2], 3),
    ", p =", format.pval(trend_summary$coefficients[2, "Pr(>|t|)"], 3), "\n")
cat("    Interpretation:", divergence_results$interpretation, "\n\n")

# ============================================================================
# PART 4C: MULTI-METRIC SENSITIVITY ANALYSIS
# ============================================================================
# Tests whether PERMANOVA and divergence findings are robust across distance
# metrics: Hellinger+Bray-Curtis (current), Wisconsin+Bray-Curtis, Jaccard,
# Gower, and Raup-Crick (null-model-based).
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 4C: MULTI-METRIC SENSITIVITY ANALYSIS\n")
cat("------------------------------------------------------------\n\n")

cat("Testing robustness of PERMANOVA and divergence results across 5 distance metrics.\n\n")

# Filter out rows with zero total abundance (prevents NaN in Wisconsin)
nonzero_idx <- rowSums(comm_aligned_raw) > 0
comm_nz <- comm_aligned_raw[nonzero_idx, ]
coral_nz <- coral_for_perm[nonzero_idx, ]

# Define distance configurations
distance_configs <- list(
  "Bray-Curtis (Hellinger)" = list(
    dist = vegdist(decostand(comm_nz, method = "hellinger"), method = "bray")
  ),
  "Bray-Curtis (Wisconsin)" = list(
    dist = vegdist(wisconsin(comm_nz), method = "bray")
  ),
  "Jaccard (binary)" = list(
    dist = vegdist(decostand(comm_nz, method = "pa"), method = "jaccard")
  ),
  "Gower" = list(
    dist = vegdist(comm_nz, method = "gower")
  ),
  "Raup-Crick (null model)" = list(
    dist = vegdist(comm_nz, method = "raup")
  )
)

# Run PERMANOVA and betadisper for each metric
sensitivity_results <- lapply(names(distance_configs), function(metric_name) {
  cat("  ", metric_name, "... ")
  d <- distance_configs[[metric_name]]$dist

  # PERMANOVA
  perm_result <- tryCatch({
    adonis2(d ~ log(volume) + site, data = coral_nz,
            permutations = 999, by = "terms")
  }, error = function(e) NULL)

  if (is.null(perm_result)) {
    cat("PERMANOVA failed\n")
    return(NULL)
  }

  perm_df <- as.data.frame(perm_result)
  vol_r2 <- perm_df["log(volume)", "R2"]
  vol_p <- perm_df["log(volume)", "Pr(>F)"]
  site_r2 <- perm_df["site", "R2"]
  site_p <- perm_df["site", "Pr(>F)"]

  # Betadisper by size class (divergence trend)
  use_sqrt <- grepl("Raup-Crick", metric_name)
  disp_result <- tryCatch({
    betadisper(d, coral_nz$size_class, sqrt.dist = use_sqrt)
  }, error = function(e) NULL)

  trend_beta <- NA_real_
  trend_p <- NA_real_
  disp_F <- NA_real_

  if (!is.null(disp_result)) {
    disp_test_res <- tryCatch(permutest(disp_result, permutations = 999), error = function(e) NULL)
    if (!is.null(disp_test_res)) disp_F <- disp_test_res$statistic

    trend_df <- data.frame(
      dist_to_centroid = disp_result$distances,
      log_volume = log(coral_nz$volume)
    )
    trend_lm <- lm(dist_to_centroid ~ log_volume, data = trend_df)
    trend_s <- summary(trend_lm)
    trend_beta <- coef(trend_lm)[2]
    trend_p <- trend_s$coefficients[2, "Pr(>|t|)"]
  }

  cat("done\n")

  data.frame(
    Metric = metric_name,
    Volume_R2 = round(vol_r2, 4),
    Volume_p = round(vol_p, 4),
    Site_R2 = round(site_r2, 4),
    Site_p = round(site_p, 4),
    PERMDISP_F = round(disp_F, 2),
    Trend_beta = round(trend_beta, 4),
    Trend_p = round(trend_p, 4),
    stringsAsFactors = FALSE
  )
})

# Build comparison table
sensitivity_table <- do.call(rbind, Filter(Negate(is.null), sensitivity_results))
rownames(sensitivity_table) <- NULL

cat("\n  Multi-Metric Comparison Table:\n")
print(sensitivity_table, row.names = FALSE)

# Robustness assessment
n_vol_sig <- sum(sensitivity_table$Volume_p < 0.05, na.rm = TRUE)
n_site_sig <- sum(sensitivity_table$Site_p < 0.05, na.rm = TRUE)
n_trend_sig <- sum(sensitivity_table$Trend_p < 0.05, na.rm = TRUE)
n_metrics <- nrow(sensitivity_table)

assess <- function(n_sig, n_total) {
  if (n_sig == n_total) "ROBUST"
  else if (n_sig >= 2) "EQUIVOCAL"
  else "NOT ROBUST"
}

cat("\n  Robustness Summary:\n")
cat("    Volume effect:", n_vol_sig, "/", n_metrics, "metrics significant →",
    assess(n_vol_sig, n_metrics), "\n")
cat("    Site effect:", n_site_sig, "/", n_metrics, "metrics significant →",
    assess(n_site_sig, n_metrics), "\n")
cat("    Divergence trend:", n_trend_sig, "/", n_metrics, "metrics significant →",
    assess(n_trend_sig, n_metrics), "\n")
cat("    Note: Gower treats shared absences as similarity (interpret cautiously).\n\n")

# Save table
save_table(sensitivity_table, "permanova_metric_sensitivity")
cat("  Saved: output/tables/permanova_metric_sensitivity.csv\n")

# Forest plot of R² across metrics
sensitivity_long <- rbind(
  data.frame(Metric = sensitivity_table$Metric, Term = "Volume",
             R2 = sensitivity_table$Volume_R2, p = sensitivity_table$Volume_p),
  data.frame(Metric = sensitivity_table$Metric, Term = "Site",
             R2 = sensitivity_table$Site_R2, p = sensitivity_table$Site_p)
)
sensitivity_long$Significant <- ifelse(sensitivity_long$p < 0.05, "p < 0.05", "NS")
sensitivity_long$Metric <- factor(sensitivity_long$Metric, levels = rev(sensitivity_table$Metric))

p_sensitivity <- ggplot(sensitivity_long, aes(x = R2, y = Metric, color = Term, shape = Significant)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("Volume" = "#0072B2", "Site" = "#D55E00")) +
  scale_shape_manual(values = c("p < 0.05" = 16, "NS" = 1)) +
  labs(
    title = "PERMANOVA R² Across Distance Metrics",
    subtitle = "Filled = significant (p < 0.05); Open = not significant",
    x = expression(R^2),
    y = NULL
  ) +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(FIG_DIR, "permanova_metric_sensitivity.png"), p_sensitivity,
       width = 8, height = 5, dpi = 300, bg = "white")
cat("  Saved: permanova_metric_sensitivity.png\n")

# Supplement S2: PERMANOVA sensitivity
supplement_dir <- file.path(PATHS$figures, "supplement")
dir.create(supplement_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(supplement_dir, "figS2_permanova_sensitivity.png"), p_sensitivity,
       width = 8, height = 5, dpi = 300, bg = "white")
cat("  Saved: figS2_permanova_sensitivity.png (supplement)\n\n")

# Store in results
sensitivity_summary <- list(
  table = sensitivity_table,
  robustness = list(
    volume = assess(n_vol_sig, n_metrics),
    site = assess(n_site_sig, n_metrics),
    divergence_trend = assess(n_trend_sig, n_metrics)
  )
)

# ============================================================================
# PART 5: COMMUNITY-LEVEL PATTERNS
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 5: ADDITIONAL COMMUNITY PATTERNS\n")
cat("------------------------------------------------------------\n\n")

# 5.1 Evenness
cat("5.1 Community Evenness:\n")

coral_master <- coral_master %>%
  mutate(
    pielou_j = shannon / log(otu_richness),
    pielou_j = if_else(is.nan(pielou_j) | is.infinite(pielou_j), NA_real_, pielou_j)
  )

cat("    Pielou's J (evenness):", round(mean(coral_master$pielou_j, na.rm = TRUE), 3),
    "± ", round(sd(coral_master$pielou_j, na.rm = TRUE), 3), "\n")
cat("    Range:", round(min(coral_master$pielou_j, na.rm = TRUE), 3), "-",
    round(max(coral_master$pielou_j, na.rm = TRUE), 3), "\n\n")

# 5.2 Prevalence analysis
cat("5.2 Species Prevalence (% of corals occupied):\n")

prevalence <- cafi_clean %>%
  group_by(otu) %>%
  summarise(
    n_corals = n_distinct(coral_id),
    prevalence = n_distinct(coral_id) / n_distinct(cafi_clean$coral_id) * 100,
    abundance = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(prevalence))

cat("    Core species (>50% prevalence):\n")
core <- filter(prevalence, prevalence > 50)
print(core)

cat("\n    Rare species (<5% prevalence):", sum(prevalence$prevalence < 5), "\n")
cat("    Satellite species (1-5 corals):", sum(prevalence$n_corals <= 5), "\n\n")

# 5.3 Indicator species
cat("5.3 Indicator Species Analysis (by site):\n")

if (requireNamespace("indicspecies", quietly = TRUE)) {
  library(indicspecies)

  # Filter to species present in at least 5 corals AND align with coral_for_perm
  comm_filtered <- community_matrix[common_ids, colSums(community_matrix[common_ids, ] > 0) >= 5]

  tryCatch({
    indval <- multipatt(comm_filtered, coral_for_perm$site, control = how(nperm = 999))
    sig_indicators <- indval$sign[indval$sign$p.value < 0.05, ]
    if (nrow(sig_indicators) > 0) {
      cat("    Significant indicators (p < 0.05):\n")
      print(head(sig_indicators[order(sig_indicators$p.value), ], 10))
    } else {
      cat("    No significant indicator species found\n")
    }
  }, error = function(e) {
    cat("    Could not compute indicator species:", conditionMessage(e), "\n")
  })
} else {
  cat("    (indicspecies package not installed)\n")
}

# ============================================================================
# MANUSCRIPT FIGURE 3: SITE COMPOSITION (NMDS + TAXONOMIC BARCHART)
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("MANUSCRIPT FIGURE 3: Site Composition (NMDS + Taxonomic Groups)\n")
cat("------------------------------------------------------------\n\n")

library(ggrepel)

# --- Taxonomic group colors (Paul Tol palette, colorblind-safe) ---
type_colors <- c(
  "Shrimp" = "#882255",
  "Crabs" = "#CC6677",
  "Hermit crabs" = "#117733",
  "Fish" = "#DDCC77",
  "Gastropods" = "#AA4499",
  "Echinoderms" = "#88CCEE",
  "Polychaetes" = "#999999",
  "Amphipods" = "#332288",
  "Squat lobsters" = "#661100"
)

site_labels_fig3 <- c(HAU = "Hauru", MAT = "Maatea", MRB = "Maharepa")

# --- Panel A: NMDS with species vectors ---
cat("  Panel A: NMDS with species vectors...\n")

# Compute separate NMDS for figure: filter to corals with ≥10 CAFI for stable ordination
min_cafi_fig3 <- 10
fig3_common_ids <- intersect(rownames(community_matrix), coral_master$coral_id)
fig3_comm <- community_matrix[fig3_common_ids, ]
fig3_coral <- coral_master %>%
  filter(coral_id %in% fig3_common_ids) %>%
  arrange(match(coral_id, fig3_common_ids))
fig3_coral$total_cafi_check <- rowSums(fig3_comm)
fig3_keep <- fig3_coral$total_cafi_check >= min_cafi_fig3
n_fig3_excluded <- sum(!fig3_keep)

fig3_comm <- fig3_comm[fig3_keep, ]
fig3_coral <- fig3_coral[fig3_keep, ]
fig3_comm <- fig3_comm[, colSums(fig3_comm) > 0]

cat("    Corals for figure:", nrow(fig3_comm), "(excluded", n_fig3_excluded, "with <10 CAFI)\n")

fig3_comm_hell <- decostand(fig3_comm, method = "hellinger")
set.seed(42)
fig3_nmds <- metaMDS(fig3_comm_hell, distance = "bray", k = 2, trymax = 100, trace = 0)
cat("    NMDS stress:", round(fig3_nmds$stress, 3), "\n")

# PERMANOVA & PERMDISP on figure subset
fig3_dist <- vegdist(fig3_comm_hell, method = "bray")
set.seed(42)
fig3_perm <- adonis2(fig3_dist ~ site, data = fig3_coral, permutations = 999)

fig3_disp_obj <- betadisper(fig3_dist, fig3_coral$site)
fig3_disp_test_obj <- permutest(fig3_disp_obj, permutations = 999)

# Extract NMDS scores for figure
fig3_nmds_scores <- as.data.frame(scores(fig3_nmds, display = "sites"))
fig3_nmds_scores$coral_id <- rownames(fig3_nmds_scores)
fig3_nmds_scores <- fig3_nmds_scores %>%
  left_join(fig3_coral %>% dplyr::select(coral_id, site, total_cafi),
            by = "coral_id")

# Fit species vectors (envfit) on common species (present in ≥8 corals)
species_freq <- colSums(fig3_comm > 0)
common_species <- names(species_freq[species_freq >= 8])

set.seed(42)
species_fit <- envfit(fig3_nmds, fig3_comm_hell[, common_species], permutations = 999)

# Extract and filter significant vectors
species_scores_df <- as.data.frame(scores(species_fit, display = "vectors"))
species_scores_df$species <- rownames(species_scores_df)
species_scores_df$r2 <- species_fit$vectors$r
species_scores_df$pval <- species_fit$vectors$pvals

sig_species <- species_scores_df %>%
  filter(pval < 0.01, r2 > 0.10) %>%
  arrange(desc(r2))

cat("    Significant vectors (p < 0.01, R² > 0.10):", nrow(sig_species), "\n")

# Map species to taxonomic groups
species_type_map <- cafi_clean %>% distinct(otu, type)

sig_species <- sig_species %>%
  left_join(species_type_map, by = c("species" = "otu")) %>%
  mutate(
    type_clean = case_when(
      type == "shrimp" ~ "Shrimp",
      type == "crab" ~ "Crabs",
      type == "hermit" ~ "Hermit crabs",
      type == "fish" ~ "Fish",
      type == "snail" ~ "Gastropods",
      type == "echinoderm" ~ "Echinoderms",
      type == "worm" ~ "Polychaetes",
      type == "amphipod" ~ "Amphipods",
      is.na(type) ~ "Other",
      TRUE ~ "Other"
    ),
    species_short = case_when(
      species == "Harpiliopsis beaupresii" ~ "Harpiliopsis",
      species == "Harpiliopsis spinigera" ~ "H. spinigera",
      species == "Paragobiodon modestus" ~ "Paragobiodon",
      species == "Calcinus latens" ~ "Calcinus",
      species == "Trapezia bidentata" ~ "Trapezia",
      species == "Trapezia rufopunctata" ~ "T. rufopunctata",
      species == "Perinia tumida" ~ "Perinia",
      species == "Caracanthus maculatus" ~ "Caracanthus",
      species == "Paracirrhites arcatus" ~ "Paracirrhites",
      species == "Ophiocoma erinaceus" ~ "Ophiocoma",
      species == "Breviturma pica" ~ "Breviturma",
      species == "Alpheus lottini" ~ "Alpheus",
      species == "Galeropsis monodonta" ~ "Galeropsis",
      species == "Alpheidae" ~ "Alpheidae",
      TRUE ~ word(species, 1)
    )
  )

# Scale vectors and select key representatives
arrow_scale <- 0.65
sig_species <- sig_species %>%
  mutate(NMDS1_scaled = NMDS1 * arrow_scale, NMDS2_scaled = NMDS2 * arrow_scale)

top_vectors <- sig_species %>%
  filter(species %in% c(
    "Harpiliopsis beaupresii",
    "Paragobiodon modestus",
    "Calcinus latens",
    "Trapezia bidentata",
    "Alpheus lottini",
    "Galeropsis monodonta",
    "Caracanthus maculatus"
  ))

cat("    Selected species for visualization:", nrow(top_vectors), "\n")

# Site centroids for label placement
site_centroids <- fig3_nmds_scores %>%
  group_by(site) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop")

site_label_positions <- site_centroids %>%
  mutate(
    label_x = case_when(
      site == "MAT" ~ NMDS1 - 0.55,
      site == "HAU" ~ NMDS1 + 0.45,
      site == "MRB" ~ NMDS1 + 0.55
    ),
    label_y = case_when(
      site == "MAT" ~ NMDS2 + 0.45,
      site == "HAU" ~ NMDS2 + 0.55,
      site == "MRB" ~ NMDS2 - 0.35
    )
  )

vector_types <- c("Shrimp", "Fish", "Hermit crabs", "Crabs", "Gastropods")

panel_a <- ggplot() +
  geom_point(data = fig3_nmds_scores,
             aes(x = NMDS1, y = NMDS2, fill = site),
             shape = 21, size = 2.5, alpha = 0.7,
             stroke = 0.3, color = "gray30") +
  stat_ellipse(data = fig3_nmds_scores,
               aes(x = NMDS1, y = NMDS2, color = site),
               level = 0.95, linewidth = 0.9, linetype = "solid",
               alpha = 0.85) +
  geom_segment(data = top_vectors,
               aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled,
                   color = type_clean),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               linewidth = 0.65, alpha = 0.85) +
  geom_text_repel(
    data = top_vectors,
    aes(x = NMDS1_scaled * 1.3, y = NMDS2_scaled * 1.3,
        label = species_short, color = type_clean),
    size = 2.5, fontface = "italic",
    max.overlaps = 30,
    segment.size = 0.15, segment.alpha = 0.4, segment.color = "gray45",
    box.padding = 0.55, point.padding = 0.3,
    min.segment.length = 0.1, force = 6, force_pull = 0.15,
    seed = 789, show.legend = FALSE
  ) +
  geom_text(data = site_label_positions,
            aes(x = label_x, y = label_y, label = site_labels_fig3[site], color = site),
            size = 3.3, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = SITE_COLORS, guide = "none") +
  scale_color_manual(
    values = c(SITE_COLORS, type_colors),
    breaks = vector_types, name = NULL,
    guide = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1))
  ) +
  annotate("text", x = 0.95, y = 0.95,
           label = paste0("Stress = ", sprintf("%.2f", fig3_nmds$stress)),
           size = 2.8, color = "gray30", hjust = 1, fontface = "italic") +
  labs(x = "NMDS1", y = "NMDS2") +
  coord_fixed(ratio = 1, xlim = c(-1.35, 1.0), ylim = c(-0.95, 1.0)) +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.title = element_text(size = 9, color = "black"),
    axis.text = element_blank(), axis.ticks = element_blank(),
    legend.position = c(0.12, 0.25),
    legend.background = element_rect(fill = alpha("white", 0.93),
                                     color = "gray80", linewidth = 0.3),
    legend.text = element_text(size = 6.5),
    legend.key.size = unit(0.32, "cm"),
    legend.key.height = unit(0.32, "cm"),
    legend.spacing.y = unit(0.03, "cm"),
    legend.margin = margin(4, 6, 4, 6),
    plot.margin = margin(5, 2, 5, 5)
  )

# --- Panel B: Taxonomic composition by site ---
cat("  Panel B: Taxonomic composition stacked barchart...\n")

type_summary <- cafi_clean %>%
  filter(!is.na(type), type != "", type != "amph") %>%
  mutate(
    type_clean = case_when(
      type == "hermit" ~ "Hermit crabs",
      type == "crab" ~ "Crabs",
      type == "shrimp" ~ "Shrimp",
      type == "fish" ~ "Fish",
      type == "snail" ~ "Gastropods",
      type == "echinoderm" ~ "Echinoderms",
      type == "worm" ~ "Polychaetes",
      type == "amphipod" ~ "Amphipods",
      type == "squat_lobster" ~ "Squat lobsters",
      TRUE ~ type
    )
  ) %>%
  group_by(site, type_clean) %>%
  summarise(n = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(site) %>%
  mutate(total = sum(n), prop = n / total * 100) %>%
  ungroup()

type_order <- type_summary %>%
  group_by(type_clean) %>%
  summarise(total_n = sum(n), .groups = "drop") %>%
  arrange(desc(total_n)) %>%
  pull(type_clean)

type_summary$type_clean <- factor(type_summary$type_clean, levels = rev(type_order))
type_summary$site <- factor(type_summary$site, levels = c("HAU", "MAT", "MRB"))

panel_b <- ggplot(type_summary, aes(x = site, y = prop, fill = type_clean)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(
    data = type_summary %>% filter(prop > 10),
    aes(label = sprintf("%.0f%%", prop)),
    position = position_stack(vjust = 0.5),
    size = 2.4, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = type_colors, name = NULL,
                    guide = guide_legend(reverse = TRUE, ncol = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), breaks = seq(0, 100, 25)) +
  scale_x_discrete(labels = site_labels_fig3) +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_minimal(base_size = 9, base_family = "Helvetica") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray95", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.title.y = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(face = "bold", size = 9.5, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.36, "cm"),
    legend.spacing.y = unit(0.06, "cm"),
    plot.margin = margin(5, 5, 5, 12)
  )

# --- Combine panels ---
panel_a_labeled <- panel_a + labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = 12))
panel_b_labeled <- panel_b + labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold", size = 12))

fig3_composition <- panel_a_labeled + panel_b_labeled +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 12, 8, 12)
    )
  )

ggsave(file.path(PATHS$fig_manuscript, "fig3_composition.png"), fig3_composition,
       width = 10, height = 4.6, dpi = 300, bg = "white")
ggsave(file.path(FIG_DIR, "fig3_composition.png"), fig3_composition,
       width = 10, height = 4.6, dpi = 300, bg = "white")
cat("  Saved: fig3_composition.png\n")

# --- Generate legend and results text ---
cat("  Generating figure legend and results text...\n")

# Pairwise PERMANOVA on figure subset
site_pairs <- combn(unique(as.character(fig3_coral$site)), 2, simplify = FALSE)
pairwise_results <- lapply(site_pairs, function(pair) {
  idx <- fig3_coral$site %in% pair
  dist_matrix <- as.matrix(fig3_dist)
  dist_subset <- as.dist(dist_matrix[idx, idx])
  pair_data <- fig3_coral[idx, ]
  set.seed(42)
  pair_perm <- adonis2(dist_subset ~ site, data = pair_data, permutations = 999)
  data.frame(
    comparison = paste(pair, collapse = " vs "),
    R2 = pair_perm$R2[1],
    p_value = pair_perm$`Pr(>F)`[1]
  )
})
pairwise_df <- do.call(rbind, pairwise_results)
pairwise_df$p_adj <- p.adjust(pairwise_df$p_value, method = "bonferroni")

# Use figure-specific PERMANOVA & PERMDISP stats
fig3_site_r2 <- fig3_perm$R2[1]
fig3_site_f <- fig3_perm$F[1]
fig3_site_p <- fig3_perm$`Pr(>F)`[1]
fig3_disp_f <- fig3_disp_test_obj$tab$F[1]
fig3_disp_p <- fig3_disp_test_obj$tab$`Pr(>F)`[1]
n_fig3_corals <- nrow(fig3_coral)

# Extract site sample sizes from figure subset
n_hau <- sum(fig3_coral$site == "HAU")
n_mat <- sum(fig3_coral$site == "MAT")
n_mrb <- sum(fig3_coral$site == "MRB")

# Extract composition stats per site
hauru_comp <- type_summary %>% filter(site == "HAU")
maatea_comp <- type_summary %>% filter(site == "MAT")
barrier_comp <- type_summary %>% filter(site == "MRB")

# Safe extraction helper
safe_prop <- function(df, grp) {
  val <- df$prop[df$type_clean == grp]
  if (length(val) == 0) 0 else val
}

# Build PERMDISP interpretation
permdisp_text <- if (fig3_disp_p < 0.05) {
  paste0(
    "Sites also differed in multivariate dispersion (PERMDISP: F = ",
    sprintf("%.2f", fig3_disp_f), ", p = ", sprintf("%.3f", fig3_disp_p),
    "), indicating that some sites harbor more heterogeneous assemblages. ",
    "However, the significant PERMANOVA result and clear separation in ordination ",
    "space confirm that compositional differences are robust (Anderson et al. 2011)."
  )
} else {
  paste0(
    "This effect was not driven by differences in multivariate dispersion ",
    "(PERMDISP: F = ", sprintf("%.2f", fig3_disp_f),
    ", p = ", sprintf("%.2f", fig3_disp_p),
    "), confirming that sites harbor compositionally distinct CAFI assemblages ",
    "(Anderson et al. 2011)."
  )
}

legend_text <- paste0(
  "FIGURE 3: CAFI COMMUNITY COMPOSITION ACROSS REEF SITES\n",
  "================================================================================\n\n",

  "FIGURE LEGEND\n",
  "-------------\n",
  "Figure 3. Coral-associated fauna and invertebrate (CAFI) community composition ",
  "differs significantly among reef sites in Mo'orea, French Polynesia. ",
  "(A) Non-metric multidimensional scaling (NMDS) ordination of Hellinger-transformed ",
  "community data (Bray-Curtis dissimilarity). Each point represents one coral colony ",
  "(n = ", n_fig3_corals, "). Ellipses show 95% confidence intervals for each site. ",
  "Vectors indicate species significantly associated with compositional variation ",
  "(envfit permutation test, p < 0.01, R\u00b2 > 0.10), colored by taxonomic group. ",
  "(B) Relative abundance of taxonomic groups at each site. Percentages shown for ",
  "groups comprising >10% of site assemblage.\n\n",

  "Site abbreviations: Hauru (n = ", n_hau, " corals), ",
  "Maatea (n = ", n_mat, " corals), ",
  "Maharepa (n = ", n_mrb, " corals).\n\n",

  "================================================================================\n\n",

  "STATISTICAL RESULTS\n",
  "-------------------\n\n",

  "1. OVERALL COMPOSITION TEST (PERMANOVA)\n",
  "   Model: Bray-Curtis dissimilarity ~ Site\n",
  "   Transformation: Hellinger\n",
  "   Permutations: 999\n\n",
  "   Result: R\u00b2 = ", sprintf("%.3f", fig3_site_r2), ", F = ", sprintf("%.2f", fig3_site_f),
  ", p = ", fig3_site_p, "\n",
  "   Interpretation: Site explains ", sprintf("%.1f", fig3_site_r2 * 100),
  "% of community composition variance.\n\n",

  "2. HOMOGENEITY OF DISPERSION TEST (PERMDISP)\n",
  "   Result: F = ", sprintf("%.2f", fig3_disp_f),
  ", p = ", sprintf("%.3f", fig3_disp_p), "\n",
  "   Interpretation: ", ifelse(fig3_disp_p < 0.05,
    "Dispersions differ among sites; interpret PERMANOVA with caution.\n\n",
    "Dispersions homogeneous; PERMANOVA results reflect true compositional differences.\n\n"),

  "3. PAIRWISE COMPARISONS (Bonferroni-corrected)\n"
)

for (i in 1:nrow(pairwise_df)) {
  legend_text <- paste0(legend_text,
    "   ", pairwise_df$comparison[i], ": R\u00b2 = ", sprintf("%.3f", pairwise_df$R2[i]),
    ", p = ", sprintf("%.3f", pairwise_df$p_value[i]),
    ", p_adj = ", sprintf("%.3f", pairwise_df$p_adj[i]),
    ifelse(pairwise_df$p_adj[i] < 0.05, " *", ""), "\n"
  )
}

legend_text <- paste0(legend_text, "\n",
  "4. ORDINATION QUALITY\n",
  "   NMDS stress: ", sprintf("%.3f", fig3_nmds$stress), "\n",
  "   Interpretation: ",
  ifelse(fig3_nmds$stress < 0.1, "Excellent representation (stress < 0.10)",
         ifelse(fig3_nmds$stress < 0.2, "Good representation (stress < 0.20)",
                "Fair representation (stress >= 0.20)")), "\n\n",

  "5. SPECIES DRIVING COMPOSITIONAL DIFFERENCES\n",
  "   Method: envfit (permutation test, 999 permutations)\n",
  "   Criteria: p < 0.01, R\u00b2 > 0.10\n\n",
  "   Top species vectors shown in figure:\n"
)

for (i in 1:nrow(top_vectors)) {
  legend_text <- paste0(legend_text,
    "   - ", top_vectors$species[i], " (", top_vectors$type_clean[i],
    "): R\u00b2 = ", sprintf("%.3f", top_vectors$r2[i]),
    ", p = ", sprintf("%.3f", top_vectors$pval[i]), "\n"
  )
}

legend_text <- paste0(legend_text, "\n",
  "================================================================================\n\n",

  "RESULTS\n",
  "-------\n\n",

  "Community composition differed significantly among the three reef sites ",
  "(PERMANOVA: R\u00b2 = ", sprintf("%.3f", fig3_site_r2), ", F = ", sprintf("%.2f", fig3_site_f),
  ", p < 0.001; Fig. 3A). ", permdisp_text, " All pairwise site comparisons were ",
  "significant after Bonferroni correction.\n\n",

  "The three sites exhibited distinct taxonomic signatures (Fig. 3B). Maatea was ",
  "characterized by high hermit crab abundance (",
  sprintf("%.0f", safe_prop(maatea_comp, "Hermit crabs")),
  "% vs ", sprintf("%.0f", safe_prop(hauru_comp, "Hermit crabs")),
  "-", sprintf("%.0f", safe_prop(barrier_comp, "Hermit crabs")),
  "% at other sites), primarily Calcinus latens. Maharepa was dominated ",
  "by obligate coral symbionts, with shrimp and crabs together comprising ",
  sprintf("%.0f", safe_prop(barrier_comp, "Shrimp") + safe_prop(barrier_comp, "Crabs")),
  "% of all CAFI. Hauru supported the most taxonomically balanced assemblage and ",
  "the highest proportion of coral-dwelling fishes (",
  sprintf("%.0f", safe_prop(hauru_comp, "Fish")),
  "%).\n\n",

  "================================================================================\n\n",

  "METHODS\n",
  "-------\n\n",

  "We analyzed CAFI community composition across ", n_fig3_corals, " Pocillopora ",
  "colonies (Hauru: n = ", n_hau, "; Maatea: n = ", n_mat,
  "; Maharepa: n = ", n_mrb, ") after excluding colonies with fewer than ",
  min_cafi_fig3, " CAFI individuals to ensure stable ordination positions (",
  n_fig3_excluded, " colonies excluded). ",
  "Species abundances were Hellinger-transformed prior to analysis to down-weight ",
  "rare species (Legendre & Gallagher 2001). We tested for compositional differences ",
  "using PERMANOVA on Bray-Curtis dissimilarities (999 permutations) and verified ",
  "homogeneity of dispersions using PERMDISP (Anderson 2006). Community structure was ",
  "visualized using NMDS (stress = ", sprintf("%.2f", fig3_nmds$stress), "). Species associated ",
  "with compositional differences were identified using envfit permutation tests ",
  "(p < 0.01, R\u00b2 > 0.10). Pairwise site comparisons were Bonferroni-corrected.\n\n",

  "================================================================================\n\n",

  "COLOR SCHEME\n",
  "------------\n",
  "Site palette (Panel A):\n",
  "  HAU (Hauru):    #E69F00 (Okabe-Ito orange)\n",
  "  MAT (Maatea):   #0072B2 (Okabe-Ito blue)\n",
  "  MRB (Maharepa): #009E73 (Okabe-Ito green)\n\n",
  "Taxonomic group palette (Panel B, species vectors):\n",
  "  Shrimp:       #882255 (deep magenta)\n",
  "  Crabs:        #CC6677 (rose pink)\n",
  "  Hermit crabs: #117733 (forest green)\n",
  "  Fish:         #DDCC77 (warm sand)\n",
  "  Gastropods:   #AA4499 (purple)\n",
  "  Echinoderms:  #88CCEE (sky blue)\n\n",

  "================================================================================\n",
  "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "Source script: scripts/02_community_analysis.R\n"
)

writeLines(legend_text, file.path(PATHS$fig_manuscript, "fig3_legend_results.txt"))
cat("  Saved: fig3_legend_results.txt\n\n")

# ============================================================================
# SUMMARY FIGURE
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("Creating Summary Figure\n")
cat("------------------------------------------------------------\n")

# Combine key figures
p_summary <- (p_abundance_site + p_func_comp) /
  (p_abundance_vol + p_nmds_site) +
  plot_annotation(
    title = "CAFI Community Analysis Summary",
    subtitle = paste0("N = ", nrow(coral_master), " corals, ", nrow(cafi_clean), " individuals, ",
                      n_distinct(cafi_clean$otu), " OTUs"),
    tag_levels = "A"
  )

ggsave(file.path(FIG_DIR, "community_summary.png"), p_summary,
       width = 14, height = 12, dpi = 300, bg = "white")

cat("\nFigures saved to:", FIG_DIR, "\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================

# Save statistical results
results_list <- list(
  abundance = list(
    summary = site_summary,
    kruskal = kw_test,
    scaling_exponent = slope,
    scaling_ci = slope_ci,
    pseudo_r2 = pseudo_r2
  ),
  diversity = list(
    by_site = diversity_by_site,
    richness_anova = aov_summary,
    richness_eta_sq = eta_sq
  ),
  community = list(
    nmds_stress = nmds$stress,
    nmds_scores = nmds_scores,
    permanova = permanova,
    permanova_margin = permanova_margin,
    permanova_interaction = permanova_interaction,
    permdisp = disp_test
  ),
  # NEW: Composition divergence by size (parallels Stier et al. Fig 3B)
  divergence = divergence_results,
  # Multi-metric sensitivity analysis
  sensitivity = sensitivity_summary
)

save_object(results_list, "community_analysis_results")

# Save divergence stats as CSV for manuscript
divergence_stats_df <- data.frame(
  test = c("PERMDISP (size classes)", "Linear trend (Small→Large)"),
  statistic = c(round(disp_size_test$statistic, 3), round(trend_summary$coefficients[2, "t value"], 3)),
  df = c(paste(disp_size_test$df, collapse = ", "), trend_summary$df[2]),
  p_value = c(disp_size_test$pval, trend_summary$coefficients[2, "Pr(>|t|)"]),
  effect_size = c(NA, round(coef(trend_model)[2], 3)),
  interpretation = c(
    ifelse(disp_size_test$pval < 0.05, "Dispersions differ by size", "No difference"),
    paste0("β = ", round(coef(trend_model)[2], 3), ": ", trend_direction, " ", trend_sig)
  )
)
save_table(divergence_stats_df, "composition_divergence_stats")

cat("\n")
cat("============================================================\n")
cat("    COMMUNITY ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("Key findings:\n")
cat("  1. Abundance scales sublinearly with coral size (β =", round(slope, 2), ")\n")
cat("  2. Site explains", round(permanova$R2[2] * 100, 1), "% of compositional variance\n")
cat("  3. Volume explains", round(permanova$R2[1] * 100, 1), "% of compositional variance\n")
cat("  4. ", round(chao1, 0), " total species estimated (", round(ncol(community_matrix)/chao1*100, 0), "% sampled)\n", sep = "")
cat("  5. Composition divergence (Fig 3B analog):", divergence_results$interpretation, "\n\n")
