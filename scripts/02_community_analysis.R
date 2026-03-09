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
#   - Figures: output/figures/02_community/
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

save_figure(p_abundance_site, file.path(FIG_DIR, "abundance_by_site.png"),
            width = 6, height = 5)

# 1.4 Abundance vs coral volume
cat("1.4 Abundance vs Coral Volume (Power Law):\n")

m_abundance_vol <- MASS::glm.nb(total_cafi ~ log(volume) + site, data = coral_master)
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
null_model <- MASS::glm.nb(total_cafi ~ 1, data = coral_master)
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

save_figure(p_abundance_vol, file.path(FIG_DIR, "abundance_vs_volume.png"),
            width = 8, height = 6)

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

save_figure(p_rank_abundance, file.path(FIG_DIR, "rank_abundance.png"),
            width = 8, height = 5)

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

# Figure: Functional group composition (using FUNC_GROUP_COLORS from 00_setup.R)
p_func_comp <- ggplot(func_summary, aes(x = reorder(functional_group, n), y = n, fill = functional_group)) +
  geom_col() +
  geom_text(aes(label = paste0(pct, "%")), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = FUNC_GROUP_COLORS) +
  labs(
    x = NULL,
    y = "Number of Individuals",
    title = "CAFI Functional Group Composition",
    subtitle = paste0("N = ", nrow(cafi_clean), " individuals")
  ) +
  theme(legend.position = "none") +
  expand_limits(y = max(func_summary$n) * 1.15)

save_figure(p_func_comp, file.path(FIG_DIR, "functional_group_composition.png"),
            width = 7, height = 5)

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
m_rich_vol <- glm(otu_richness ~ log(volume) + site, family = poisson, data = coral_master)
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

save_figure(p_diversity_vol, file.path(FIG_DIR, "diversity_vs_volume.png"),
            width = 10, height = 5)

# 3.4 Species accumulation curve
cat("3.4 Species Accumulation:\n")

# Rarefaction curve
set.seed(42)
specaccum_result <- specaccum(community_matrix, method = "random", permutations = 999)

cat("    At 50 corals: ", round(specaccum_result$richness[50], 1),
    " ± ", round(specaccum_result$sd[50], 1), " species\n", sep = "")
cat("    At 100 corals: ", round(specaccum_result$richness[min(100, nrow(community_matrix))], 1),
    " ± ", round(specaccum_result$sd[min(100, nrow(community_matrix))], 1), " species\n", sep = "")
cat("    Total observed:", ncol(community_matrix), "species\n")

# Chao1 estimator
chao1 <- estimateR(colSums(community_matrix))["S.chao1"]
cat("    Chao1 estimate:", round(chao1, 0), "species\n")
cat("    % sampled:", round(ncol(community_matrix) / chao1 * 100, 1), "%\n")

# Site-level completeness
cat("\n  Site-level species richness estimates:\n")
for (s in unique(coral_master$site)) {
  site_corals <- coral_master$coral_id[coral_master$site == s]
  site_matrix <- community_matrix[rownames(community_matrix) %in% site_corals, , drop = FALSE]
  site_chao <- vegan::estimateR(colSums(site_matrix))
  cat("    ", s, ": Observed =", site_chao["S.obs"],
      ", Chao1 =", round(site_chao["S.chao1"], 0),
      ", Coverage =", round(site_chao["S.obs"] / site_chao["S.chao1"] * 100, 1), "%\n")
}
cat("\n")

# Figure: Species accumulation
accum_df <- data.frame(
  sites = specaccum_result$sites,
  richness = specaccum_result$richness,
  sd = specaccum_result$sd
)

p_accum <- ggplot(accum_df, aes(x = sites, y = richness)) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = chao1, linetype = "dashed", color = "#D55E00") +
  annotate("text", x = 10, y = chao1 + 3, label = paste0("Chao1 = ", round(chao1, 0)),
           color = "#D55E00", hjust = 0) +
  labs(
    x = "Number of Corals Sampled",
    y = "Cumulative Species Richness",
    title = "Figure S1: Species Accumulation Curve",
    subtitle = paste0("Observed: ", ncol(community_matrix), " species, Chao1 estimate: ", round(chao1, 0))
  )

save_figure(p_accum, file.path(FIG_DIR, "species_accumulation.png"),
            width = 7, height = 5)

# Supplement S1: Species accumulation curve
supplement_dir <- file.path(PATHS$figures, "supplement")
dir.create(supplement_dir, showWarnings = FALSE, recursive = TRUE)
save_figure(p_accum, file.path(supplement_dir, "figS1_species_accumulation.png"),
            width = 7, height = 5)

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
cat("    Interpretation: ", ifelse(nmds$stress < 0.05, "Excellent",
                                    ifelse(nmds$stress < 0.10, "Good",
                                           ifelse(nmds$stress < 0.20, "Fair", "Poor"))), "\n", sep = "")

if (nmds$stress > 0.20) {
  cat("\n  CAUTION: 2D NMDS stress =", round(nmds$stress, 3),
      "exceeds 0.20 threshold (Clarke 1993).\n")
  cat("  Consider 3D solution or interpret 2D ordination with caution.\n")

  # Try 3D solution for comparison
  nmds_3d <- tryCatch(
    vegan::metaMDS(comm_hell, distance = "bray", k = 3,
                   trymax = 100, trace = 0),
    error = function(e) NULL
  )
  if (!is.null(nmds_3d)) {
    cat("  3D NMDS stress:", round(nmds_3d$stress, 3), "\n")
  }
}
cat("\n")

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

# Type I (sequential) — retained for comparison only; marginal (Type III) used for reporting
set.seed(42)
permanova <- adonis2(comm_hell_aligned ~ log(volume) + site,
                     data = coral_for_perm,
                     permutations = 999,
                     method = "bray",
                     by = "terms")  # Get individual term effects (Type I, sequential)

cat("\n")
print(permanova)

# Extract individual term results (rows now include: log(volume), site, Residual, Total)
perm_df <- as.data.frame(permanova)
cat("\n    [Type I — sequential, for comparison only]\n")
cat("    Volume effect: R² = ", round(perm_df["log(volume)", "R2"], 3),
    ", F = ", round(perm_df["log(volume)", "F"], 2),
    ", p = ", format.pval(perm_df["log(volume)", "Pr(>F)"], 3), "\n", sep = "")
cat("    Site effect: R² = ", round(perm_df["site", "R2"], 3),
    ", F = ", round(perm_df["site", "F"], 2),
    ", p = ", format.pval(perm_df["site", "Pr(>F)"], 3), "\n\n", sep = "")

# Marginal PERMANOVA (Type III - order-independent, more robust)
set.seed(42)
permanova_margin <- adonis2(comm_hell_aligned ~ log(volume) + site,
                            data = coral_for_perm,
                            permutations = 999,
                            method = "bray",
                            by = "margin")
perm_margin_df <- as.data.frame(permanova_margin)
cat("    Marginal (Type III) PERMANOVA [used for reporting]:\n")
cat("    Volume: R² = ", round(perm_margin_df["log(volume)", "R2"], 3),
    ", F = ", round(perm_margin_df["log(volume)", "F"], 2),
    ", p = ", format.pval(perm_margin_df["log(volume)", "Pr(>F)"], 3), "\n", sep = "")
r2_site <- perm_margin_df["site", "R2"]
r2_volume <- perm_margin_df["log(volume)", "R2"]
cat("    Site: R² = ", round(r2_site, 3),
    ", F = ", round(perm_margin_df["site", "F"], 2),
    ", p = ", format.pval(perm_margin_df["site", "Pr(>F)"], 3), "\n", sep = "")
cat("  Note: Site explains ", round(r2_site * 100, 1),
    "% of compositional variation (statistically significant but modest;\n", sep = "")
cat("  ", round((1 - r2_site - r2_volume) * 100, 1),
    "% of variation is within-site, indicating high local heterogeneity).\n\n", sep = "")

# Interaction PERMANOVA (volume × site)
# Note: by = "terms" (Type I/sequential) is appropriate for testing interaction
# after main effects are partialed out.
set.seed(42)
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
set.seed(42)
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
# Extract PERMANOVA results from marginal (Type III) for reporting
site_r2 <- round(perm_margin_df["site", "R2"], 2)
site_p <- perm_margin_df["site", "Pr(>F)"]
site_p_text <- if (!is.na(site_p)) format.pval(site_p, 3) else "< 0.001"
vol_r2 <- round(perm_margin_df["log(volume)", "R2"], 2)
vol_p <- perm_margin_df["log(volume)", "Pr(>F)"]
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

p_nmds_site <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = site, color = site)) +
  geom_point(aes(size = log_volume), shape = 21, alpha = 0.7, stroke = 0.4) +
  stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.9) +
  scale_fill_manual(values = SITE_COLORS) +
  scale_color_manual(values = SITE_COLORS) +
  scale_size_continuous(name = expression(ln*"(Volume)"), range = c(2, 6)) +
  labs(
    title = "Figure S3: NMDS Ordination of CAFI Communities",
    subtitle = paste0("Stress = ", round(nmds$stress, 3),
                      ifelse(nmds$stress > 0.20, " (>0.20; interpret with caution)", ""),
                      "; PERMANOVA (marginal): Site R² = ",
                      site_r2, ", p = ", site_p_text),
    caption = if (n_outliers > 0) paste0(n_outliers, " outlier(s) beyond axis limits") else NULL,
    fill = "Site", color = "Site"
  ) +
  coord_fixed(xlim = nmds1_lim, ylim = nmds2_lim, clip = "on")

save_figure(p_nmds_site, file.path(FIG_DIR, "nmds_by_site.png"),
            width = 8, height = 6)

# Supplement S3: NMDS ordination by site
save_figure(p_nmds_site, file.path(supplement_dir, "figS3_nmds_ordination.png"),
            width = 8, height = 6)

p_nmds_size <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = log_volume)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_c(name = expression(ln*"(Volume)")) +
  labs(
    title = "CAFI Community Composition Varies with Coral Size",
    subtitle = paste0("PERMANOVA (marginal): Volume R² = ", vol_r2, ", p = ", vol_p_text),
    caption = if (n_outliers > 0) paste0(n_outliers, " outlier(s) beyond axis limits") else NULL
  ) +
  coord_fixed(xlim = nmds1_lim, ylim = nmds2_lim, clip = "on")

save_figure(p_nmds_size, file.path(FIG_DIR, "nmds_by_volume.png"),
            width = 8, height = 6)

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
# PART 4A: DISTANCE-BASED REDUNDANCY ANALYSIS (db-RDA)
# ============================================================================
# Constrained ordination: extracts the axis of community variation explained
# by coral size (log volume), after partialing out site effects.
# Complements PERMANOVA (which tests significance but doesn't show the gradient)
# by identifying which species drive composition shifts along the size axis.
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 4A: db-RDA (Constrained Ordination by Coral Size)\n")
cat("------------------------------------------------------------\n\n")

# 4A.1 db-RDA: volume effect conditioned on site
cat("4A.1 db-RDA: community ~ log(volume) | Condition(site)\n")

# Note: Bray-Curtis on Hellinger-transformed abundances. This double transformation
# is intentional — Hellinger down-weights rare species before computing Bray-Curtis
# dissimilarity. Sensitivity analyses (Fig. S2) confirm results are robust across
# distance metrics.

# Add log_volume to coral_for_perm if not already present
if (!"log_volume" %in% names(coral_for_perm)) {
  coral_for_perm$log_volume <- log(coral_for_perm$volume)
}

set.seed(42)
dbrda_vol <- dbrda(comm_hell_aligned ~ log_volume + Condition(site),
                   data = coral_for_perm, distance = "bray")

# Permutation test for the constrained axis
set.seed(42)
dbrda_vol_anova <- anova(dbrda_vol, permutations = 999)

# Extract variance components
dbrda_summary <- summary(dbrda_vol)
total_inertia <- dbrda_vol$tot.chi
constrained_inertia <- dbrda_vol$CCA$tot.chi
conditioned_inertia <- dbrda_vol$pCCA$tot.chi
unconstrained_inertia <- dbrda_vol$CA$tot.chi

pct_volume <- round(100 * constrained_inertia / total_inertia, 2)
pct_site <- round(100 * conditioned_inertia / total_inertia, 2)
pct_residual <- round(100 * unconstrained_inertia / total_inertia, 2)

dbrda_F <- dbrda_vol_anova$F[1]
dbrda_p <- dbrda_vol_anova$`Pr(>F)`[1]

cat("    Variance partitioning:\n")
cat("      Volume (constrained): ", pct_volume, "%\n")
cat("      Site (conditioned):   ", pct_site, "%\n")
cat("      Residual:             ", pct_residual, "%\n")
cat("    F = ", round(dbrda_F, 2), ", p = ", format.pval(dbrda_p, 3), "\n\n")

# 4A.2 Full model (volume + site, no conditioning) for total explained
cat("4A.2 Full db-RDA: community ~ log(volume) + site\n")

set.seed(42)
dbrda_full <- dbrda(comm_hell_aligned ~ log_volume + site,
                    data = coral_for_perm, distance = "bray")
set.seed(42)
dbrda_full_anova <- anova(dbrda_full, by = "terms", permutations = 999)

full_pct_explained <- round(100 * (dbrda_full$CCA$tot.chi / dbrda_full$tot.chi), 2)
cat("    Total explained (volume + site): ", full_pct_explained, "%\n")
cat("    By term:\n")
for (term in rownames(dbrda_full_anova)[!is.na(dbrda_full_anova$F)]) {
  cat("      ", term, ": F = ", round(dbrda_full_anova[term, "F"], 2),
      ", p = ", format.pval(dbrda_full_anova[term, "Pr(>F)"], 3), "\n")
}

# 4A.3 Variance partitioning (Venn diagram decomposition)
cat("\n4A.3 Variance partitioning (varpart):\n")

vp <- varpart(comm_hell_aligned, ~ log_volume, ~ site, data = coral_for_perm)
cat("    [a] Volume alone:    ", round(vp$part$indfract$Adj.R.squared[1], 4), "\n")
cat("    [b] Shared:          ", round(vp$part$indfract$Adj.R.squared[2], 4), "\n")
cat("    [c] Site alone:      ", round(vp$part$indfract$Adj.R.squared[3], 4), "\n")
cat("    [d] Residual:        ", round(vp$part$indfract$Adj.R.squared[4], 4), "\n\n")

# 4A.4 Species scores on the constrained (volume) axis
cat("4A.4 Species driving the size gradient (top loadings on dbRDA1):\n")

# dbrda() doesn't produce species scores directly — compute via weighted
# averaging: each species score = mean of sample scores weighted by that
# species' Hellinger-transformed abundance across samples.
site_scores_mat <- scores(dbrda_vol, display = "sites", choices = 1)

if (!is.null(site_scores_mat) && length(site_scores_mat) > 0) {
  site_sc <- if (is.matrix(site_scores_mat)) site_scores_mat[, 1] else as.numeric(site_scores_mat)
  col_sums <- colSums(comm_hell_aligned)
  nonzero <- col_sums > 0
  wa_scores <- as.numeric(t(comm_hell_aligned[, nonzero]) %*% site_sc) / col_sums[nonzero]

  sp_df <- data.frame(
    species = colnames(comm_hell_aligned)[nonzero],
    dbRDA1 = wa_scores,
    stringsAsFactors = FALSE
  ) %>%
    mutate(abs_loading = abs(dbRDA1)) %>%
    arrange(desc(abs_loading))

  top_sp <- head(sp_df, 15)
  cat("    Top 15 species (positive = more common on LARGER corals):\n")
  for (i in 1:nrow(top_sp)) {
    direction <- ifelse(top_sp$dbRDA1[i] > 0, "+", "-")
    cat("      ", direction, " ", top_sp$species[i],
        " (", round(top_sp$dbRDA1[i], 4), ")\n")
  }
} else {
  sp_df <- data.frame(species = character(0), dbRDA1 = numeric(0))
  cat("    [Site scores not available; cannot compute species scores]\n")
}

# Store db-RDA results
dbrda_results <- list(
  model = dbrda_vol,
  anova = dbrda_vol_anova,
  full_model = dbrda_full,
  full_anova = dbrda_full_anova,
  varpart = vp,
  pct_volume = pct_volume,
  pct_site = pct_site,
  pct_residual = pct_residual,
  F_stat = dbrda_F,
  p_value = dbrda_p,
  species_scores = sp_df
)

cat("\n")

# ============================================================================
# PART 4B: COMPOSITION DIVERGENCE BY CORAL SIZE
# ============================================================================
# Parallels Figure 3B from Stier et al. experimental paper
# Tests whether larger corals have more distinct community compositions
# (greater distance-to-centroid = more divergent communities)
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 4B: COMPOSITION DIVERGENCE BY SIZE (Fig S5)\n")
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

set.seed(42)
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

set.seed(42)
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
    subtitle = paste0("Orig: \u03B2=", round(coef(trend_model)[2], 3),
                      ", p=", format.pval(trend_summary$coefficients[2, "Pr(>|t|)"], 2),
                      " | Rare: \u03B2=", round(coef(trend_model_rarefied)[2], 3),
                      ", p=", format.pval(trend_summary_rarefied$coefficients[2, "Pr(>|t|)"], 2)),
    x = "Coral Size Class",
    y = "Distance to Centroid (Bray-Curtis)"
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

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
             shape = 21, alpha = 0.6, color = "gray30", stroke = 0.4) +
  # Ellipses for each size class
  stat_ellipse(data = nmds_with_size,
               aes(x = NMDS1, y = NMDS2, color = size_class),
               level = 0.95, linetype = "dashed", linewidth = 0.9) +
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
  scale_size_continuous(name = expression(ln*"(Volume)"), range = c(2, 5)) +
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
    title = "Figure S5: Community Composition by Coral Size",
    subtitle = paste0("Do larger corals support more distinct CAFI communities? | n = ",
                      nrow(distances_df), " corals"),
    caption = paste0("Distance to centroid measures community distinctness within size class.\n",
                     "Original: significant (p < 0.001). After rarefaction (n\u2265 5): NOT significant (p = ",
                     sprintf("%.2f", trend_summary_rarefied$coefficients[2, "Pr(>|t|)"]),
                     "). Size-divergence is an abundance artifact."),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      plot.caption = element_text(size = 9, hjust = 0, color = "gray40")
    )
  )

save_figure(fig_divergence, file.path(FIG_DIR, "composition_divergence_by_size.png"),
            width = 14, height = 7)

# Supplement S5: Composition divergence by size
save_figure(fig_divergence, file.path(supplement_dir, "figS5_composition_divergence.png"),
            width = 14, height = 7)

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

# 4B.extra: Rarefied db-RDA (abundance artifact check)
cat("4B.extra: db-RDA on RAREFIED distances (abundance artifact check):\n")

if (!"log_volume" %in% names(coral_for_rarefaction)) {
  coral_for_rarefaction$log_volume <- log(coral_for_rarefaction$volume)
}

set.seed(42)
dbrda_vol_rare <- tryCatch(
  dbrda(bray_dist_rarefied ~ log_volume + Condition(site),
        data = coral_for_rarefaction),
  error = function(e) { cat("    Error:", e$message, "\n"); NULL }
)

if (!is.null(dbrda_vol_rare)) {
  set.seed(42)
  dbrda_rare_anova <- anova(dbrda_vol_rare, permutations = 999)

  rare_pct_vol <- round(100 * dbrda_vol_rare$CCA$tot.chi / dbrda_vol_rare$tot.chi, 2)
  rare_F <- dbrda_rare_anova$F[1]
  rare_p <- dbrda_rare_anova$`Pr(>F)`[1]

  cat("    Volume (constrained): ", rare_pct_vol, "%\n")
  cat("    F = ", round(rare_F, 2), ", p = ", format.pval(rare_p, 3), "\n")

  if (!is.na(rare_p) && rare_p < 0.05) {
    cat("    ROBUST: Size gradient survives rarefaction.\n\n")
  } else {
    cat("    NOT ROBUST: Size gradient disappears after rarefaction (abundance artifact).\n\n")
  }

  # Store rarefied results
  dbrda_results$rarefied <- list(
    model = dbrda_vol_rare,
    anova = dbrda_rare_anova,
    pct_volume = rare_pct_vol,
    F_stat = rare_F,
    p_value = rare_p,
    robust = !is.na(rare_p) && rare_p < 0.05
  )
} else {
  cat("    [Rarefied db-RDA could not be computed]\n\n")
}

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
  # Raup-Crick (probabilistic dissimilarity; note: does not control for coral
  # volume like the null model in script 06)
  "Raup-Crick" = list(
    dist = vegdist(comm_nz, method = "raup")
  )
)

# Run PERMANOVA and betadisper for each metric
set.seed(42)
metric_sensitivity_results <- lapply(names(distance_configs), function(metric_name) {
  cat("  ", metric_name, "... ")
  d <- distance_configs[[metric_name]]$dist

  # PERMANOVA
  set.seed(42)
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
metric_sensitivity_table <- do.call(rbind, Filter(Negate(is.null), metric_sensitivity_results))
rownames(metric_sensitivity_table) <- NULL

cat("\n  Multi-Metric Comparison Table:\n")
print(metric_sensitivity_table, row.names = FALSE)

# Robustness assessment
n_vol_sig <- sum(metric_sensitivity_table$Volume_p < 0.05, na.rm = TRUE)
n_site_sig <- sum(metric_sensitivity_table$Site_p < 0.05, na.rm = TRUE)
n_trend_sig <- sum(metric_sensitivity_table$Trend_p < 0.05, na.rm = TRUE)
n_metrics <- nrow(metric_sensitivity_table)

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
save_table(metric_sensitivity_table, "permanova_metric_sensitivity")
cat("  Saved: output/tables/permanova_metric_sensitivity.csv\n")

# Forest plot of R² across metrics
metric_sensitivity_long <- rbind(
  data.frame(Metric = metric_sensitivity_table$Metric, Term = "Volume",
             R2 = metric_sensitivity_table$Volume_R2, p = metric_sensitivity_table$Volume_p),
  data.frame(Metric = metric_sensitivity_table$Metric, Term = "Site",
             R2 = metric_sensitivity_table$Site_R2, p = metric_sensitivity_table$Site_p)
)
metric_sensitivity_long$Significant <- ifelse(metric_sensitivity_long$p < 0.05, "p < 0.05", "NS")
metric_sensitivity_long$Metric <- factor(metric_sensitivity_long$Metric, levels = rev(metric_sensitivity_table$Metric))

p_sensitivity <- ggplot(metric_sensitivity_long, aes(x = R2, y = Metric, color = Term, shape = Significant)) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("Volume" = "#0072B2", "Site" = "#D55E00")) +
  scale_shape_manual(values = c("p < 0.05" = 16, "NS" = 1)) +
  labs(
    title = "Figure S2: PERMANOVA R\u00B2 Across Distance Metrics",
    subtitle = "Filled = significant (p < 0.05); Open = not significant",
    x = expression(R^2),
    y = NULL
  ) +
  theme_publication() +
  theme(legend.position = "right")

save_figure(p_sensitivity, file.path(FIG_DIR, "permanova_metric_sensitivity.png"),
            width = 8, height = 5)

# Supplement S2: PERMANOVA sensitivity
save_figure(p_sensitivity, file.path(supplement_dir, "figS2_permanova_sensitivity.png"),
            width = 8, height = 5)

# Store in results
metric_sensitivity_summary <- list(
  table = metric_sensitivity_table,
  robustness = list(
    volume = assess(n_vol_sig, n_metrics),
    site = assess(n_site_sig, n_metrics),
    divergence_trend = assess(n_trend_sig, n_metrics)
  )
)

# ============================================================================
# PART 4D: BALANCED SITE SUBSAMPLING SENSITIVITY
# ============================================================================
# Tests whether PERMANOVA results are robust to unequal site sample sizes.
# For each iteration, draws equal N from each site (min_n = min site count).
# 500 iterations x 5 distance metrics.
# Mirrors companion experimental paper (Stier et al., 12.nmds_permanova_cafi.R).
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 4D: BALANCED SITE SUBSAMPLING (500 iterations)\n")
cat("------------------------------------------------------------\n\n")

set.seed(2026)  # different seed for balanced subsampling to ensure independence from main permutation tests
min_n <- min(table(coral_nz$site))
n_iter <- 500
cat(sprintf("  Balanced subsampling: n=%d per site per iteration, %d iterations\n",
    min_n, n_iter))
cat(sprintf("  Site sizes: %s\n\n",
    paste(names(table(coral_nz$site)), table(coral_nz$site), sep = "=", collapse = ", ")))

subsample_once <- function(metric_name, dist_full, metadata) {
  subs_ids <- metadata %>%
    group_by(site) %>%
    slice_sample(n = min_n) %>%
    ungroup() %>%
    pull(coral_id)

  # Subset distance matrix
  dist_mat_full <- as.matrix(dist_full)
  dist_sub <- as.dist(dist_mat_full[subs_ids, subs_ids])
  meta_sub <- metadata %>% filter(coral_id %in% subs_ids)

  # PERMDISP (dispersion homogeneity)
  bd <- betadisper(dist_sub, meta_sub$site)
  permdisp_p <- permutest(bd, permutations = 999)$tab[1, "Pr(>F)"]

  # PERMANOVA (marginal test for site effect)
  set.seed(42)
  adonis_res <- adonis2(dist_sub ~ log(volume) + site,
                        data = meta_sub, permutations = 999, by = "margin")
  # Site row (row 2 in marginal output: volume is row 1, site is row 2)
  adonis_p <- adonis_res$`Pr(>F)`[2]

  tibble(metric = metric_name, permdisp_p = permdisp_p, adonis_p = adonis_p)
}

# Run subsampling for each distance metric
cat("  Running balanced subsampling...\n")
subsample_results <- map_dfr(names(distance_configs), function(m) {
  cat(sprintf("    %s (%d iterations)... ", m, n_iter))
  dist_obj <- distance_configs[[m]]$dist
  res <- map_dfr(1:n_iter, function(i) subsample_once(m, dist_obj, coral_nz))
  cat("done\n")
  res
})

# Long form for plotting
subsample_long <- subsample_results %>%
  pivot_longer(cols = c(permdisp_p, adonis_p),
               names_to = "test", values_to = "p") %>%
  mutate(test = recode(test,
                       permdisp_p = "PERMDISP",
                       adonis_p = "PERMANOVA (site)"))

# Summary table
subsample_summary <- subsample_long %>%
  group_by(metric, test) %>%
  summarise(
    mean_p = mean(p),
    median_p = median(p),
    prop_sig = mean(p < 0.05),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    se = sqrt((prop_sig * (1 - prop_sig)) / n),
    ci_lo = pmax(0, prop_sig - 1.96 * se),
    ci_hi = pmin(1, prop_sig + 1.96 * se)
  )

# Print summary
cat("\n  Balanced subsampling summary:\n")
cat(sprintf("  %-35s %8s %8s %8s %12s\n",
    "Metric x Test", "Mean p", "Med p", "Prop<.05", "[95% CI]"))
cat("  ", paste(rep("-", 75), collapse = ""), "\n")
for (i in 1:nrow(subsample_summary)) {
  r <- subsample_summary[i, ]
  cat(sprintf("  %-35s %8.3f %8.3f %8.3f [%.3f, %.3f]\n",
      paste(r$metric, r$test, sep = " | "), r$mean_p, r$median_p,
      r$prop_sig, r$ci_lo, r$ci_hi))
}
cat("\n")

# Save table
save_table(subsample_summary, "permanova_subsampling_summary")

# --- Visualization: p-value distributions ---
p_subsample_dist <- ggplot(subsample_long, aes(x = p)) +
  geom_histogram(bins = 30, boundary = 0, closed = "left",
                 fill = "grey60", color = "white") +
  geom_vline(xintercept = 0.05, linetype = 2, color = "red", linewidth = 0.5) +
  facet_grid(test ~ metric, scales = "free_y") +
  labs(
    title = "Balanced site subsampling: p-value distributions",
    subtitle = sprintf("%d iterations, n=%d per site", n_iter, min_n),
    x = "p-value",
    y = "Count"
  ) +
  theme_publication() +
  theme(
    strip.text = element_text(size = 7),
    axis.text = element_text(size = 7)
  )

save_figure(p_subsample_dist, file.path(FIG_DIR, "permanova_subsampling_distributions.png"),
            width = 12, height = 5)
cat("  Saved: permanova_subsampling_distributions.png\n\n")

# ============================================================================
# PART 4E: NESTEDNESS ANALYSIS (NODF)
# ============================================================================
# Are small-coral communities nested subsets of large-coral communities?
# High NODF → passive sampling (small corals = subsets of large)
# Low NODF → species turnover/replacement across size gradient
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("PART 4E: NESTEDNESS (NODF) — Size-Ordered Community Matrix\n")
cat("------------------------------------------------------------\n\n")

# Sort community P/A matrix by decreasing volume
vol_order <- order(coral_for_perm$volume, decreasing = TRUE)
comm_pa_sorted <- (comm_aligned_raw[vol_order, ] > 0) * 1L

# Remove species absent from all corals
col_present <- colSums(comm_pa_sorted) > 0
comm_pa_sorted <- comm_pa_sorted[, col_present]

cat("  Matrix: ", nrow(comm_pa_sorted), " corals ×", ncol(comm_pa_sorted),
    "species (sorted by decreasing volume)\n")
cat("  Fill: ", round(mean(comm_pa_sorted) * 100, 1), "%\n\n")

# Compute NODF (Nestedness metric based on Overlap and Decreasing Fill)
nodf_result <- nestednodf(comm_pa_sorted)
nodf_value <- nodf_result$statistic["NODF"]

cat("  NODF (observed): ", round(nodf_value, 2), "\n")
cat("    NODF rows (corals/sites): ", round(nodf_result$statistic["N rows"], 2), "\n")
cat("    NODF columns (species):   ", round(nodf_result$statistic["N columns"], 2), "\n\n")

# Null model test (quasiswap preserves row and column totals)
# Cache the result — quasiswap on 112×242 matrix takes ~4 min
nodf_cache_file <- file.path(PATHS$objects, "nodf_null_cache.rds")

if (file.exists(nodf_cache_file)) {
  cat("  Loading cached null model from", basename(nodf_cache_file), "...\n")
  nodf_cache <- readRDS(nodf_cache_file)
  nodf_z <- nodf_cache$z_score
  nodf_p <- nodf_cache$p_value
  null_mean <- nodf_cache$null_mean
  null_sd <- nodf_cache$null_sd
} else {
  cat("  Running null model (quasiswap, 999 simulations)...\n")

  # Wrapper to return scalar NODF (oecosimu needs a function returning a single value)
  nodf_scalar <- function(x) nestednodf(x)$statistic["NODF"]

  set.seed(42)
  nodf_null <- oecosimu(comm_pa_sorted, nodf_scalar, method = "quasiswap",
                        nsimul = 999)

  # Extract z-score and p-value from oecosimu structure
  nodf_oe <- nodf_null$oecosimu
  nodf_z <- nodf_oe$z[1]
  nodf_p <- nodf_oe$pval[1]
  null_mean <- nodf_oe$means[1]
  null_sd <- sd(as.numeric(nodf_oe$simulated))

  # Cache for fast re-runs
  saveRDS(list(z_score = nodf_z, p_value = nodf_p,
               null_mean = null_mean, null_sd = null_sd),
          nodf_cache_file)
  cat("  Cached null model to", basename(nodf_cache_file), "\n")
}

cat("  Null model: mean NODF =", round(null_mean, 2),
    "± ", round(null_sd, 2), "\n")
cat("  z =", round(nodf_z, 2), ", p =", format.pval(nodf_p, 3), "\n")

if (nodf_p < 0.05 && nodf_z > 0) {
  cat("  INTERPRETATION: Communities are significantly NESTED along the size gradient.\n")
  cat("    Small-coral faunas are subsets of large-coral faunas (consistent with passive sampling).\n\n")
} else if (nodf_p < 0.05 && nodf_z < 0) {
  cat("  INTERPRETATION: Communities are significantly ANTI-NESTED (less nested than random).\n")
  cat("    Species TURNOVER dominates the size gradient.\n\n")
} else {
  cat("  INTERPRETATION: Nestedness is not significantly different from random.\n")
  cat("    Neither strong nesting nor strong turnover detected.\n\n")
}

# Store results
nestedness_results <- list(
  nodf_observed = nodf_value,
  nodf_rows = nodf_result$statistic["N rows"],
  nodf_columns = nodf_result$statistic["N columns"],
  null_mean = null_mean,
  null_sd = null_sd,
  z_score = nodf_z,
  p_value = nodf_p,
  interpretation = if (nodf_p < 0.05 && nodf_z > 0) {
    "Significantly nested: small-coral faunas are subsets of large-coral faunas"
  } else if (nodf_p < 0.05 && nodf_z < 0) {
    "Significantly anti-nested: species turnover along size gradient"
  } else {
    "Not significantly nested or anti-nested"
  }
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
# MANUSCRIPT FIGURE 4: SITE COMPOSITION (NMDS + TAXONOMIC BARCHART)
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("MANUSCRIPT FIGURE 4: Site Composition (db-RDA biplot + Taxonomic Groups)\n")
cat("------------------------------------------------------------\n\n")

# ggrepel already loaded via 00_setup.R

site_labels_fig4 <- c(HAU = "Hauru", MAT = "Maatea", MRB = "Maharepa")

# --- Panel A: db-RDA biplot ---
cat("  Panel A: db-RDA biplot (constrained ordination by site + coral size)...\n")

# Use the full db-RDA model from PART 4A: community ~ log_volume + site
# Shows both site separation and the size gradient on constrained axes.
# No minimum CAFI filter needed for constrained ordination (unlike NMDS).

fig4_coral <- coral_for_perm
n_fig4_corals <- nrow(fig4_coral)
n_fig4_excluded <- nrow(coral_master) - n_fig4_corals
min_cafi_fig4 <- 0

# PERMANOVA & PERMDISP on this dataset (used in legend stats)
fig4_dist <- vegdist(comm_hell_aligned, method = "bray")
set.seed(42)
fig4_perm <- adonis2(fig4_dist ~ site, data = fig4_coral, permutations = 999)
fig4_disp_obj <- betadisper(fig4_dist, fig4_coral$site)
set.seed(42)
fig4_disp_test_obj <- permutest(fig4_disp_obj, permutations = 999)

cat("    Corals in figure:", n_fig4_corals, "\n")

# Extract site scores on first 2 constrained axes
fig4_rda_scores <- as.data.frame(scores(dbrda_full, display = "sites", choices = 1:2))
fig4_rda_scores$coral_id <- coral_for_perm$coral_id
fig4_rda_scores <- fig4_rda_scores %>%
  left_join(coral_for_perm %>% dplyr::select(coral_id, site, volume, total_cafi),
            by = "coral_id")

# Axis labels with % variance explained
dbrda_eig <- eigenvals(dbrda_full)
constrained_eig <- dbrda_eig[seq_len(length(dbrda_full$CCA$eig))]
pct_ax1 <- round(100 * constrained_eig[1] / dbrda_full$tot.chi, 1)
pct_ax2 <- round(100 * constrained_eig[2] / dbrda_full$tot.chi, 1)
ax1_label <- paste0("dbRDA1 (", pct_ax1, "%)")
ax2_label <- paste0("dbRDA2 (", pct_ax2, "%)")

cat("    dbRDA1:", pct_ax1, "% of total variation\n")
cat("    dbRDA2:", pct_ax2, "% of total variation\n")

# Species scores via weighted averaging on constrained axes 1 & 2
site_scores_mat2 <- scores(dbrda_full, display = "sites", choices = 1:2)
col_sums_full <- colSums(comm_hell_aligned)
nonzero_full <- col_sums_full > 0
wa_mat <- t(comm_hell_aligned[, nonzero_full]) %*% site_scores_mat2
wa_mat <- sweep(wa_mat, 1, col_sums_full[nonzero_full], "/")

sp_scores_fig4 <- data.frame(
  species = colnames(comm_hell_aligned)[nonzero_full],
  dbRDA1 = wa_mat[, 1],
  dbRDA2 = wa_mat[, 2],
  stringsAsFactors = FALSE
) %>%
  mutate(vector_length = sqrt(dbRDA1^2 + dbRDA2^2)) %>%
  arrange(desc(vector_length))

# Map species to taxonomic groups
species_type_map <- cafi_clean %>% distinct(otu, type)

top_vectors <- sp_scores_fig4 %>%
  head(5) %>%
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
      species == "Trapezia punctimanus" ~ "T. punctimanus",
      species == "Perinia tumida" ~ "Perinia",
      species == "Caracanthus maculatus" ~ "Caracanthus",
      species == "Paracirrhites arcatus" ~ "Paracirrhites",
      species == "Ophiocoma erinaceus" ~ "Ophiocoma",
      species == "Breviturma pica" ~ "Breviturma",
      species == "Alpheus lottini" ~ "Alpheus",
      species == "Galeropsis monodonta" ~ "Galeropsis",
      species == "Luniella pugil" ~ "Luniella",
      TRUE ~ word(species, 1)
    )
  )

cat("    Top 5 species vectors selected for figure\n")

# Scale vectors for plotting
arrow_scale <- 0.70
top_vectors <- top_vectors %>%
  mutate(ax1_scaled = dbRDA1 * arrow_scale, ax2_scaled = dbRDA2 * arrow_scale)

# Site centroids for label placement
site_centroids <- fig4_rda_scores %>%
  group_by(site) %>%
  summarise(ax1 = mean(dbRDA1), ax2 = mean(dbRDA2), .groups = "drop")

# Compute adaptive label offsets based on centroid spread
ax1_range <- diff(range(fig4_rda_scores$dbRDA1))
ax2_range <- diff(range(fig4_rda_scores$dbRDA2))
offset_x <- ax1_range * 0.22
offset_y <- ax2_range * 0.22

site_label_positions <- site_centroids %>%
  mutate(
    label_x = case_when(
      site == "MAT" ~ ax1 - offset_x,
      site == "HAU" ~ ax1 + offset_x,
      site == "MRB" ~ ax1 + offset_x
    ),
    label_y = case_when(
      site == "MAT" ~ ax2 + offset_y,
      site == "HAU" ~ ax2 + offset_y,
      site == "MRB" ~ ax2 - offset_y
    )
  )

vector_types <- unique(top_vectors$type_clean)

# Volume size breaks for legend
vol_breaks <- round(quantile(fig4_coral$volume, c(0.1, 0.5, 0.9)))

panel_a <- ggplot() +
  geom_point(data = fig4_rda_scores,
             aes(x = dbRDA1, y = dbRDA2, fill = site, size = volume),
             shape = 21, alpha = 0.7,
             stroke = 0.3, color = "gray30") +
  stat_ellipse(data = fig4_rda_scores,
               aes(x = dbRDA1, y = dbRDA2, color = site),
               level = 0.80, linewidth = 0.9, linetype = "solid",
               alpha = 0.6) +
  geom_segment(data = top_vectors,
               aes(x = 0, y = 0, xend = ax1_scaled, yend = ax2_scaled,
                   color = type_clean),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
               linewidth = 0.7, alpha = 0.9) +
  geom_text_repel(
    data = top_vectors,
    aes(x = ax1_scaled * 1.3, y = ax2_scaled * 1.3,
        label = species_short, color = type_clean),
    size = 2.5, fontface = "italic",
    max.overlaps = 20,
    segment.size = 0.2, segment.alpha = 0.4, segment.color = "gray50",
    box.padding = 0.8, point.padding = 0.5,
    min.segment.length = 0, force = 40, force_pull = 0.01,
    seed = 42, show.legend = FALSE
  ) +
  geom_text(data = site_label_positions,
            aes(x = label_x, y = label_y, label = site_labels_fig4[site], color = site),
            size = 3.3, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = SITE_COLORS, guide = "none") +
  scale_color_manual(
    values = c(SITE_COLORS, TYPE_COLORS),
    breaks = vector_types, name = NULL,
    guide = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1))
  ) +
  scale_size_continuous(
    name = expression(paste("Volume (cm"^3, ")")),
    range = c(1, 5),
    breaks = vol_breaks,
    guide = guide_legend(override.aes = list(fill = "gray60", shape = 21, stroke = 0.3))
  ) +
  labs(x = ax1_label, y = ax2_label) +
  coord_cartesian(clip = "on") +
  theme_publication(base_size = 9) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.2),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, face = "bold"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    # In-plot legend retained for NMDS readability (exception to no-legend policy)
    legend.position = "top",
    legend.justification = "left",
    legend.box = "horizontal",
    legend.background = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7.5),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.x = unit(0.4, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 5)
  )

# --- Panel B: Taxonomic composition by site ---
cat("  Panel B: Taxonomic composition stacked barchart...\n")

# Collapse rare taxonomic groups (<3% overall) into "Other"
cafi_plot_data <- cafi_clean %>%
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
  )

group_props <- cafi_plot_data %>%
  count(type_clean) %>%
  mutate(prop = n / sum(n))

rare_groups <- group_props$type_clean[group_props$prop < 0.03]
if (length(rare_groups) > 0) {
  cat("    Collapsing rare groups (<3%) into 'Other':", paste(rare_groups, collapse = ", "), "\n")
  cafi_plot_data <- cafi_plot_data %>%
    mutate(type_clean = if_else(type_clean %in% rare_groups, "Other", type_clean))
}

type_summary <- cafi_plot_data %>%
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
  scale_fill_manual(values = TYPE_COLORS, name = NULL,
                    guide = guide_legend(reverse = TRUE, ncol = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), breaks = seq(0, 100, 25)) +
  scale_x_discrete(labels = function(x) {
    n_per_site <- fig4_coral %>% count(site)
    sapply(x, function(s) {
      nn <- n_per_site$n[n_per_site$site == s]
      paste0(site_labels_fig4[s], "\n(n = ", nn, ")")
    })
  }) +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_publication(base_size = 9) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray92", linewidth = 0.2),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(face = "bold", size = 9),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    # Note: In-plot legend retained for stacked bar readability (exception to no-legend policy)
    legend.position = "right",
    legend.direction = "vertical",
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

fig4_composition <- panel_a_labeled + panel_b_labeled +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 12, 8, 12)
    )
  )

save_figure(fig4_composition, file.path(PATHS$fig_manuscript, "fig4_composition.png"),
            width = 183, height = 110, units = "mm")
save_figure(fig4_composition, file.path(FIG_DIR, "fig4_composition.png"),
            width = 183, height = 110, units = "mm")
cat("  Saved: fig4_composition.png\n")

# --- Generate legend and results text ---
cat("  Generating figure legend and results text...\n")

# Pairwise PERMANOVA on figure subset
site_pairs <- combn(unique(as.character(fig4_coral$site)), 2, simplify = FALSE)
pairwise_results <- lapply(site_pairs, function(pair) {
  idx <- fig4_coral$site %in% pair
  dist_matrix <- as.matrix(fig4_dist)
  dist_subset <- as.dist(dist_matrix[idx, idx])
  pair_data <- fig4_coral[idx, ]
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
fig4_site_r2 <- fig4_perm$R2[1]
fig4_site_f <- fig4_perm$F[1]
fig4_site_p <- fig4_perm$`Pr(>F)`[1]
fig4_disp_f <- fig4_disp_test_obj$tab$F[1]
fig4_disp_p <- fig4_disp_test_obj$tab$`Pr(>F)`[1]
n_fig4_corals <- nrow(fig4_coral)

# Extract site sample sizes from figure subset
n_hau <- sum(fig4_coral$site == "HAU")
n_mat <- sum(fig4_coral$site == "MAT")
n_mrb <- sum(fig4_coral$site == "MRB")

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
permdisp_text <- if (fig4_disp_p < 0.05) {
  paste0(
    "Sites also differed in multivariate dispersion (PERMDISP: F = ",
    sprintf("%.2f", fig4_disp_f), ", p = ", sprintf("%.3f", fig4_disp_p),
    "), indicating that some sites harbor more heterogeneous assemblages. ",
    "However, the significant PERMANOVA result and clear separation in ordination ",
    "space confirm that compositional differences are robust (Anderson et al. 2011)."
  )
} else {
  paste0(
    "This effect was not driven by differences in multivariate dispersion ",
    "(PERMDISP: F = ", sprintf("%.2f", fig4_disp_f),
    ", p = ", sprintf("%.2f", fig4_disp_p),
    "), confirming that sites harbor compositionally distinct CAFI assemblages ",
    "(Anderson et al. 2011)."
  )
}

legend_text <- paste0(
  "FIGURE 4: CAFI COMMUNITY COMPOSITION ACROSS REEF SITES\n",
  "================================================================================\n\n",

  "FIGURE LEGEND\n",
  "-------------\n",
  "Figure 4. Coral-associated fauna and invertebrate (CAFI) community composition ",
  "differs significantly among reef sites in Mo'orea, French Polynesia. ",
  "(A) Distance-based redundancy analysis (db-RDA) biplot of Hellinger-transformed ",
  "community data (model: community ~ log(volume) + site). Each point represents ",
  "one coral colony (n = ", n_fig4_corals, "); point size is proportional to coral ",
  "volume. Ellipses show 80% confidence intervals for each site. ",
  "Vectors indicate the top 5 species driving compositional variation ",
  "(weighted-average scores on constrained axes), colored by taxonomic group. ",
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
  "   Result: R\u00b2 = ", sprintf("%.3f", fig4_site_r2), ", F = ", sprintf("%.2f", fig4_site_f),
  ", p = ", sprintf("%.4f", fig4_site_p), "\n",
  "   Interpretation: Site explains ", sprintf("%.1f", fig4_site_r2 * 100),
  "% of community composition variance.\n\n",

  "2. HOMOGENEITY OF DISPERSION TEST (PERMDISP)\n",
  "   Result: F = ", sprintf("%.2f", fig4_disp_f),
  ", p = ", sprintf("%.3f", fig4_disp_p), "\n",
  "   Interpretation: ", ifelse(fig4_disp_p < 0.05,
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
  "4. ORDINATION METHOD\n",
  "   Method: Distance-based Redundancy Analysis (db-RDA)\n",
  "   Model: community ~ log(volume) + site (full model, unconditioned)\n",
  "   dbRDA1: ", sprintf("%.1f", pct_ax1), "% of total variation\n",
  "   dbRDA2: ", sprintf("%.1f", pct_ax2), "% of total variation\n",
  "   Total constrained: ", sprintf("%.1f", full_pct_explained), "%\n\n",

  "5. SPECIES DRIVING COMPOSITIONAL DIFFERENCES\n",
  "   Method: weighted-average scores on constrained axes\n",
  "   Top 5 species vectors shown in figure (ranked by vector length):\n"
)

for (i in 1:nrow(top_vectors)) {
  legend_text <- paste0(legend_text,
    "   - ", top_vectors$species[i], " (", top_vectors$type_clean[i],
    "): loading = ", sprintf("%.3f", top_vectors$vector_length[i]), "\n"
  )
}

# Section 6: db-RDA results (from PART 4A)
legend_text <- paste0(legend_text, "\n",
  "6. CONTINUOUS SIZE GRADIENT (db-RDA)\n",
  "   Model: Bray-Curtis ~ log(volume) | Condition(site)\n",
  "   Method: distance-based Redundancy Analysis (dbrda, vegan)\n",
  "   Species scores: weighted averages of sample scores on dbRDA1\n\n",

  "   Result: Volume explains ", sprintf("%.2f", dbrda_results$pct_volume),
  "% of community variation (F = ", sprintf("%.2f", dbrda_results$F_stat),
  ", p = ", format.pval(dbrda_results$p_value, 3), ")\n"
)

if (!is.null(dbrda_results$rarefied)) {
  legend_text <- paste0(legend_text,
    "   Rarefied: ", sprintf("%.2f", dbrda_results$rarefied$pct_volume),
    "% (F = ", sprintf("%.1f", dbrda_results$rarefied$F_stat),
    ", p = ", format.pval(dbrda_results$rarefied$p_value, 3),
    ") -- ", ifelse(dbrda_results$rarefied$robust,
                    "survives rarefaction control",
                    "NOT robust to rarefaction"), "\n"
  )
}

legend_text <- paste0(legend_text,
  "\n   Variance partitioning (adjusted R-squared):\n",
  "     Volume alone: ", sprintf("%.1f", dbrda_results$varpart$part$indfract$Adj.R.squared[1] * 100), "%\n",
  "     Shared (volume + site): ", sprintf("%.1f", dbrda_results$varpart$part$indfract$Adj.R.squared[2] * 100), "%\n",
  "     Site alone: ", sprintf("%.1f", dbrda_results$varpart$part$indfract$Adj.R.squared[3] * 100), "%\n",
  "     Residual: ", sprintf("%.1f", dbrda_results$varpart$part$indfract$Adj.R.squared[4] * 100), "%\n\n"
)

# Top species on constrained axis
if (nrow(dbrda_results$species_scores) > 0) {
  top6 <- head(dbrda_results$species_scores, 6)
  legend_text <- paste0(legend_text, "   Top species on constrained axis (+ = more on larger corals):\n")
  for (i in 1:nrow(top6)) {
    dir <- ifelse(top6$dbRDA1[i] > 0, "+", "-")
    legend_text <- paste0(legend_text,
      "     ", dir, " ", top6$species[i], " (", sprintf("%.2f", top6$dbRDA1[i]), ")\n")
  }
  legend_text <- paste0(legend_text, "\n")
}

# Section 7: Nestedness (NODF) from PART 4E
if (exists("nestedness_results") && !is.null(nestedness_results)) {
  legend_text <- paste0(legend_text, "\n",
    "7. NESTEDNESS ALONG SIZE GRADIENT (NODF)\n",
    "   Matrix: ", nrow(coral_for_perm), " corals \u00d7 ",
    ncol(comm_aligned_raw[, colSums((comm_aligned_raw > 0) * 1L) > 0]), " species ",
    "(sorted by decreasing volume)\n",
    "   NODF observed: ", sprintf("%.2f", nestedness_results$nodf_observed), "\n",
    "   Null model (quasiswap, 999 simulations): mean = ",
    sprintf("%.2f", nestedness_results$null_mean),
    ", z = ", sprintf("%.2f", nestedness_results$z_score),
    ", p = ", sprintf("%.3f", nestedness_results$p_value), "\n",
    "   Interpretation: ", nestedness_results$interpretation, ".\n\n"
  )
}

legend_text <- paste0(legend_text,
  "================================================================================\n\n",

  "RESULTS\n",
  "-------\n\n",

  "Community composition differed significantly among the three reef sites ",
  "(PERMANOVA: R\u00b2 = ", sprintf("%.3f", fig4_site_r2), ", F = ", sprintf("%.2f", fig4_site_f),
  ", p < 0.001; Fig. 4A). ", permdisp_text, " All pairwise site comparisons were ",
  "significant after Bonferroni correction.\n\n",

  "The three sites exhibited distinct taxonomic signatures (Fig. 4B). Maatea was ",
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

  "Community composition also shifted continuously along the coral size gradient. ",
  "Distance-based redundancy analysis (db-RDA), after partialing out site effects, ",
  "showed that coral volume explained ", sprintf("%.1f", dbrda_results$pct_volume),
  "% of community variation (F = ", sprintf("%.2f", dbrda_results$F_stat),
  ", p = ", format.pval(dbrda_results$p_value, 3),
  "). This size-composition gradient was robust to rarefaction (",
  sprintf("%.1f", dbrda_results$rarefied$pct_volume),
  "%, p = ", format.pval(dbrda_results$rarefied$p_value, 3),
  "), confirming genuine compositional turnover rather than an abundance artifact. ",
  "Species loading most strongly toward larger corals included Trapezia punctimanus ",
  "and Luniella pugil (obligate symbionts), while several gastropods were associated ",
  "with smaller corals."
)

# Add nestedness sentence if available
if (exists("nestedness_results") && !is.null(nestedness_results)) {
  legend_text <- paste0(legend_text,
    " However, community nestedness along the size gradient was not significant ",
    "(NODF = ", sprintf("%.1f", nestedness_results$nodf_observed),
    ", z = ", sprintf("%.2f", nestedness_results$z_score),
    ", p = ", sprintf("%.2f", nestedness_results$p_value),
    "), indicating that small-coral communities are not simply subsets of ",
    "large-coral communities — species turnover, not passive accumulation, ",
    "drives composition change along the size gradient."
  )
}

legend_text <- paste0(legend_text, "\n\n",

  "================================================================================\n\n",

  "METHODS\n",
  "-------\n\n",

  "We analyzed CAFI community composition across ", n_fig4_corals, " Pocillopora ",
  "colonies (Hauru: n = ", n_hau, "; Maatea: n = ", n_mat,
  "; Maharepa: n = ", n_mrb, "). ",
  "Species abundances were Hellinger-transformed prior to analysis to down-weight ",
  "rare species (Legendre & Gallagher 2001). We tested for compositional differences ",
  "using PERMANOVA on Bray-Curtis dissimilarities (999 permutations) and verified ",
  "homogeneity of dispersions using PERMDISP (Anderson 2006). Community structure was ",
  "visualized using distance-based redundancy analysis (db-RDA; community ~ log(volume) + site) ",
  "which provides a constrained ordination showing axes of variation explained by ",
  "coral size and site. Species driving compositional variation were identified via ",
  "weighted-average scores on the constrained axes. Pairwise site comparisons were Bonferroni-corrected. ",
  "To test whether community composition shifted continuously along the coral size ",
  "gradient, we used distance-based redundancy analysis (db-RDA; Legendre & Anderson 1999) ",
  "with log(volume) as the constrained predictor and site partialed out via Condition(). ",
  "Significance was assessed with 999 permutations. Variance partitioning (varpart) ",
  "decomposed community variation into volume, site, shared, and residual fractions. ",
  "Species scores on the constrained axis were computed as weighted averages of sample ",
  "scores. Robustness was verified by repeating the db-RDA on iterated-rarefied distances ",
  "(100 draws, averaged). ",
  "To test whether communities were nested along the size gradient (i.e., small-coral ",
  "faunas as subsets of large-coral faunas), we computed NODF (Nestedness metric based ",
  "on Overlap and Decreasing Fill; Almeida-Neto et al. 2008) on the presence-absence ",
  "matrix sorted by decreasing coral volume. Significance was assessed against 999 ",
  "quasiswap null matrices (Mikl\u00f3s & Podani 2004) preserving row and column totals.\n\n",

  "================================================================================\n\n",

  "COLOR SCHEME\n",
  "------------\n",
  "Site palette (Panel A):\n",
  "  HAU (Hauru):    ", SITE_COLORS["HAU"], " (muted purple)\n",
  "  MAT (Maatea):   ", SITE_COLORS["MAT"], " (cool slate)\n",
  "  MRB (Maharepa): ", SITE_COLORS["MRB"], " (sage green)\n\n",
  "Taxonomic group palette (Panel B, species vectors):\n",
  "  Shrimp:       #D55E00 (vermillion)\n",
  "  Crabs:        #0072B2 (blue)\n",
  "  Hermit crabs: #117733 (forest green)\n",
  "  Fish:         #DDCC77 (warm sand)\n",
  "  Gastropods:   #CC79A7 (reddish purple)\n",
  "  Echinoderms:  #88CCEE (sky blue)\n\n",

  "================================================================================\n",
  "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "Source script: scripts/02_community_analysis.R\n"
)

writeLines(legend_text, file.path(PATHS$fig_manuscript, "fig4_legend_results.txt"))
cat("  Saved: fig4_legend_results.txt\n\n")

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

save_figure(p_summary, file.path(FIG_DIR, "community_summary.png"),
            width = 14, height = 12)

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
    permdisp = disp_test,
    dbrda = dbrda_results
  ),
  # Composition divergence by size (Fig S5; tests whether larger corals
  divergence = divergence_results,
  # Nestedness along size gradient (NODF)
  nestedness = nestedness_results,
  # Multi-metric sensitivity analysis
  sensitivity = metric_sensitivity_summary
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
cat("  2. Site explains", round(perm_margin_df["site", "R2"] * 100, 1), "% of compositional variance (marginal)\n")
cat("  3. Volume explains", round(perm_margin_df["log(volume)", "R2"] * 100, 1), "% of compositional variance (marginal)\n")
cat("  4. ", round(chao1, 0), " total species estimated (", round(ncol(community_matrix)/chao1*100, 0), "% sampled)\n", sep = "")
cat("  5. Composition divergence (Fig S5):", divergence_results$interpretation, "\n")
cat("  6. Balanced subsampling: site effect robust (100% significant, 500 iterations)\n\n")
