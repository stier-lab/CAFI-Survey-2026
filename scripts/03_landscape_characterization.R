# ============================================================================
# 03_landscape_characterization.R - Pocillopora Landscape Structure
# ============================================================================
#
# PURPOSE: Comprehensive characterization of Pocillopora landscape structure
#          including base metrics, extended analyses, and predictor selection
#
# IMPORTANT: Survey Design Distinction
#   - "neighborhood" corals (n=63): POC01-POC21 at each site - full 5m radius
#     surveys of all neighboring Pocillopora within 5m
#   - "size" corals (n=51): POC22+ at each site - surveyed only for size/CAFI,
#     NO neighborhood data collected (n_neighbors = NA)
#
#   This script uses ONLY corals with actual neighborhood data
#
# SECTIONS:
#   PART 1: Base Landscape Characterization
#     - Focal coral size distribution
#     - Neighborhood density
#     - Neighbor distance
#     - Neighborhood habitat volume
#   PART 2: Correlation & Predictor Selection
#     - Correlation structure
#     - VIF analysis
#     - Predictor selection
#   PART 3: Extended Analyses
#     - Size-isolation relationship
#     - Site-specific patterns
#     - Neighborhood size structure
#     - Landscape typology (clustering)
#
# OUTPUTS:
#   - Figures: output/figures/landscape/
#   - Tables: output/tables/
#   - Objects: selected predictors for 04_landscape_effects.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-10
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    POCILLOPORA LANDSCAPE CHARACTERIZATION\n")
cat("============================================================\n\n")

# Load setup and data
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

# Use script-specific output directory
FIG_DIR <- PATHS$fig_03_landscape

# ============================================================================
# PREPARE DATA
# ============================================================================

cat("Preparing landscape data...\n\n")

# Survey design reminder
cat("Survey Design:\n")
cat("   Total corals in dataset:", nrow(coral_master), "\n")
cat("   Neighborhood surveys (5m radius):", sum(!is.na(coral_master$n_neighbors)), "\n")
cat("   Size-only surveys (no neighbors):", sum(is.na(coral_master$n_neighbors)), "\n\n")

# Filter to neighborhood-surveyed corals ONLY
landscape_data <- coral_master %>%
  filter(
    !is.na(volume), volume > 0,
    !is.na(total_cafi),
    !is.na(n_neighbors)  # Filters to ONLY neighborhood-surveyed corals
  ) %>%
  mutate(
    log_volume = log10(volume),
    has_neighbors = n_neighbors > 0,

    # Recalculate derived indices
    isolation_index = mean_neighbor_dist / (volume^(1/3) + 1),
    crowding_index = total_neighbor_volume / (mean_neighbor_dist + 1),
    relative_size = volume / (mean_neighbor_volume + 1),

    # Size categories
    size_tertile = cut(volume,
                       breaks = quantile(volume, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                       labels = c("Small", "Medium", "Large"),
                       include.lowest = TRUE),

    # Isolation categories
    isolation_cat = case_when(
      n_neighbors == 0 ~ "Isolated (no neighbors within 5m)",
      mean_neighbor_dist > quantile(mean_neighbor_dist[n_neighbors > 0], 0.67, na.rm = TRUE) ~ "Distant",
      mean_neighbor_dist > quantile(mean_neighbor_dist[n_neighbors > 0], 0.33, na.rm = TRUE) ~ "Moderate",
      TRUE ~ "Close"
    ),

    # Neighborhood density categories
    density_cat = case_when(
      n_neighbors == 0 ~ "Isolated",
      n_neighbors <= 3 ~ "Sparse",
      n_neighbors <= 10 ~ "Moderate",
      TRUE ~ "Dense"
    ) %>% factor(levels = c("Isolated", "Sparse", "Moderate", "Dense")),

    # Relative position in neighborhood
    relative_position = case_when(
      n_neighbors == 0 ~ "Solitary",
      volume > mean_neighbor_volume * 1.5 ~ "Dominant",
      volume < mean_neighbor_volume * 0.67 ~ "Subordinate",
      TRUE ~ "Similar"
    ) %>% factor(levels = c("Solitary", "Subordinate", "Similar", "Dominant"))
  )

cat("ANALYSIS DATASET (neighborhood-surveyed corals only):\n")
cat("   Sample size:", nrow(landscape_data), "corals\n")
cat("   Sites:", paste(unique(landscape_data$site), collapse = ", "), "\n")
cat("   Truly isolated (0 neighbors in 5m):", sum(landscape_data$n_neighbors == 0), "\n")
cat("   Corals with neighbors:", sum(landscape_data$has_neighbors), "\n\n")

# ############################################################################
#                    PART 1: BASE LANDSCAPE CHARACTERIZATION
# ############################################################################

cat("############################################################\n")
cat("PART 1: BASE LANDSCAPE CHARACTERIZATION\n")
cat("############################################################\n\n")

# ============================================================================
# 1.1 FOCAL CORAL SIZE DISTRIBUTION
# ============================================================================

cat("------------------------------------------------------------\n")
cat("1.1 FOCAL CORAL SIZE DISTRIBUTION\n")
cat("------------------------------------------------------------\n\n")

vol_summary <- landscape_data %>%
  summarise(
    n = n(),
    min = min(volume),
    q10 = quantile(volume, 0.10),
    q25 = quantile(volume, 0.25),
    median = median(volume),
    mean = mean(volume),
    q75 = quantile(volume, 0.75),
    q90 = quantile(volume, 0.90),
    max = max(volume),
    sd = sd(volume),
    cv = sd(volume) / mean(volume) * 100,
    log_range = log10(max(volume)) - log10(min(volume))
  )

cat("Volume Distribution (cm³):\n")
cat("    N:", vol_summary$n, "corals\n")
cat("    Min:", round(vol_summary$min), "\n")
cat("    10th %ile:", round(vol_summary$q10), "\n")
cat("    25th %ile:", round(vol_summary$q25), "\n")
cat("    Median:", round(vol_summary$median), "\n")
cat("    Mean:", round(vol_summary$mean), "\n")
cat("    75th %ile:", round(vol_summary$q75), "\n")
cat("    90th %ile:", round(vol_summary$q90), "\n")
cat("    Max:", round(vol_summary$max), "\n")
cat("    SD:", round(vol_summary$sd), "\n")
cat("    CV:", round(vol_summary$cv, 1), "%\n")
cat("    Log range:", round(vol_summary$log_range, 2), "orders of magnitude\n\n")

# By site
vol_by_site <- landscape_data %>%
  group_by(site) %>%
  summarise(
    n = n(),
    median_vol = median(volume),
    mean_vol = mean(volume),
    cv = sd(volume) / mean(volume) * 100,
    .groups = "drop"
  )

cat("By Site:\n")
print(vol_by_site)
cat("\n")

# Size distribution figures
p_vol_hist <- ggplot(landscape_data, aes(x = volume, fill = site)) +
  geom_histogram(bins = 25, alpha = 0.7, position = "identity") +
  scale_x_log10(labels = scales::comma,
                breaks = c(100, 1000, 10000)) +
  scale_fill_manual(values = SITE_COLORS) +
  labs(x = expression("Coral Volume (cm"^3*", log scale)"),
       y = "Count",
       title = "A. Coral Size Distribution",
       subtitle = paste0("CV = ", round(vol_summary$cv, 0), "%; ",
                        round(vol_summary$log_range, 1), " orders of magnitude")) +
  theme(legend.position = c(0.85, 0.75))

p_vol_box <- ggplot(landscape_data, aes(x = site, y = volume, fill = site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_manual(values = SITE_COLORS) +
  labs(x = "Site",
       y = expression("Coral Volume (cm"^3*", log scale)"),
       title = "B. Size by Site") +
  theme(legend.position = "none")

# ============================================================================
# 1.2 NEIGHBORHOOD DENSITY (NUMBER OF NEIGHBORS)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("1.2 NEIGHBORHOOD DENSITY (5m radius)\n")
cat("------------------------------------------------------------\n\n")

neigh_count_summary <- landscape_data %>%
  summarise(
    n_isolated = sum(n_neighbors == 0),
    pct_isolated = mean(n_neighbors == 0) * 100,
    n_sparse = sum(n_neighbors >= 1 & n_neighbors <= 2),
    pct_sparse = mean(n_neighbors >= 1 & n_neighbors <= 2) * 100,
    n_moderate = sum(n_neighbors >= 3 & n_neighbors <= 10),
    pct_moderate = mean(n_neighbors >= 3 & n_neighbors <= 10) * 100,
    n_dense = sum(n_neighbors > 10),
    pct_dense = mean(n_neighbors > 10) * 100,
    mean_neighbors = mean(n_neighbors),
    median_neighbors = median(n_neighbors),
    max_neighbors = max(n_neighbors)
  )

cat("Neighbor Count Distribution:\n")
cat("    Isolated (0 neighbors):", neigh_count_summary$n_isolated,
    "(", round(neigh_count_summary$pct_isolated, 1), "%)\n")
cat("    Sparse (1-2 neighbors):", neigh_count_summary$n_sparse,
    "(", round(neigh_count_summary$pct_sparse, 1), "%)\n")
cat("    Moderate (3-10 neighbors):", neigh_count_summary$n_moderate,
    "(", round(neigh_count_summary$pct_moderate, 1), "%)\n")
cat("    Dense (>10 neighbors):", neigh_count_summary$n_dense,
    "(", round(neigh_count_summary$pct_dense, 1), "%)\n")
cat("    Mean:", round(neigh_count_summary$mean_neighbors, 1), "\n")
cat("    Median:", neigh_count_summary$median_neighbors, "\n")
cat("    Max:", neigh_count_summary$max_neighbors, "\n\n")

# Neighbor count figure
p_neigh_count <- ggplot(landscape_data, aes(x = n_neighbors, fill = site)) +
  geom_histogram(binwidth = 2, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = SITE_COLORS) +
  labs(x = "Number of Neighbors (5m radius)",
       y = "Count",
       title = "C. Neighborhood Density",
       subtitle = paste0(round(neigh_count_summary$pct_isolated, 0), "% isolated, ",
                        round(neigh_count_summary$pct_dense, 0), "% in dense patches")) +
  theme(legend.position = "none")

# ============================================================================
# 1.3 NEIGHBOR DISTANCE (ISOLATION)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("1.3 MEAN DISTANCE TO NEIGHBORS\n")
cat("------------------------------------------------------------\n\n")

# Only for corals with neighbors
dist_data <- landscape_data %>%
  filter(n_neighbors > 0, !is.na(mean_neighbor_dist))

dist_summary <- dist_data %>%
  summarise(
    n = n(),
    min_cm = min(mean_neighbor_dist),
    median_cm = median(mean_neighbor_dist),
    mean_cm = mean(mean_neighbor_dist),
    max_cm = max(mean_neighbor_dist),
    min_m = min(mean_neighbor_dist) / 100,
    median_m = median(mean_neighbor_dist) / 100,
    mean_m = mean(mean_neighbor_dist) / 100,
    max_m = max(mean_neighbor_dist) / 100
  )

cat("Distance Distribution (among corals with neighbors, n =", dist_summary$n, "):\n")
cat("    Min:", round(dist_summary$min_m, 2), "m (", round(dist_summary$min_cm), "cm)\n")
cat("    Median:", round(dist_summary$median_m, 2), "m (", round(dist_summary$median_cm), "cm)\n")
cat("    Mean:", round(dist_summary$mean_m, 2), "m (", round(dist_summary$mean_cm), "cm)\n")
cat("    Max:", round(dist_summary$max_m, 2), "m (", round(dist_summary$max_cm), "cm)\n\n")

# Distance figure
p_neigh_dist <- dist_data %>%
  ggplot(aes(x = mean_neighbor_dist / 100, fill = site)) +
  geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = SITE_COLORS) +
  labs(x = "Mean Distance to Neighbors (m)",
       y = "Count",
       title = "D. Neighbor Distance",
       subtitle = paste0("Mean = ", round(dist_summary$mean_m, 2), " m")) +
  theme(legend.position = "none")

# ============================================================================
# 1.4 NEIGHBORHOOD HABITAT VOLUME
# ============================================================================

cat("------------------------------------------------------------\n")
cat("1.4 NEIGHBORHOOD HABITAT AVAILABILITY\n")
cat("------------------------------------------------------------\n\n")

# Only for corals with neighbors
vol_neigh_data <- landscape_data %>%
  filter(total_neighbor_volume > 0)

neigh_vol_summary <- vol_neigh_data %>%
  summarise(
    n = n(),
    min = min(total_neighbor_volume),
    median = median(total_neighbor_volume),
    mean = mean(total_neighbor_volume),
    max = max(total_neighbor_volume),
    cv = sd(total_neighbor_volume) / mean(total_neighbor_volume) * 100
  )

cat("Total Neighbor Volume (among corals with neighbors, n =", neigh_vol_summary$n, "):\n")
cat("    Min:", round(neigh_vol_summary$min), "cm³\n")
cat("    Median:", round(neigh_vol_summary$median), "cm³\n")
cat("    Mean:", round(neigh_vol_summary$mean), "cm³\n")
cat("    Max:", round(neigh_vol_summary$max), "cm³\n")
cat("    CV:", round(neigh_vol_summary$cv, 1), "%\n\n")

# Neighbor volume figure
p_neigh_vol <- vol_neigh_data %>%
  ggplot(aes(x = total_neighbor_volume, fill = site)) +
  geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = SITE_COLORS) +
  labs(x = expression("Total Neighbor Volume (cm"^3*", log scale)"),
       y = "Count",
       title = "E. Neighborhood Habitat",
       subtitle = paste0("CV = ", round(neigh_vol_summary$cv, 0), "%")) +
  theme(legend.position = "none")

# ############################################################################
#                    PART 2: CORRELATION & PREDICTOR SELECTION
# ############################################################################

cat("############################################################\n")
cat("PART 2: CORRELATION & PREDICTOR SELECTION\n")
cat("############################################################\n\n")

# ============================================================================
# 2.1 CORRELATION STRUCTURE AMONG PREDICTORS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("2.1 CORRELATION STRUCTURE AMONG LANDSCAPE PREDICTORS\n")
cat("------------------------------------------------------------\n\n")

# Select candidate predictors
cor_vars <- landscape_data %>%
  dplyr::select(
    volume, n_neighbors, total_neighbor_volume,
    mean_neighbor_volume, mean_neighbor_dist, isolation_index,
    crowding_index, relative_size
  ) %>%
  drop_na()

cat("Correlation analysis based on n =", nrow(cor_vars), "corals with complete data\n\n")

cor_mat <- cor(cor_vars, use = "pairwise.complete.obs")

# High correlations
cat("High Correlations (|r| > 0.70):\n")
high_cor <- which(abs(cor_mat) > 0.70 & cor_mat < 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  high_cor_df <- data.frame(
    var1 = rownames(cor_mat)[high_cor[,1]],
    var2 = colnames(cor_mat)[high_cor[,2]],
    r = cor_mat[high_cor]
  ) %>%
    filter(var1 < var2) %>%
    arrange(desc(abs(r)))

  for (i in 1:nrow(high_cor_df)) {
    cat(sprintf("    %s ~ %s: r = %.2f\n",
                high_cor_df$var1[i], high_cor_df$var2[i], high_cor_df$r[i]))
  }
} else {
  cat("    None found\n")
}

# Correlations with focal volume
cat("\nCorrelation with Focal Coral Volume:\n")
vol_cors <- cor_mat["volume", ]
vol_cors <- vol_cors[names(vol_cors) != "volume"]
vol_cors <- sort(vol_cors, decreasing = TRUE)
for (nm in names(vol_cors)) {
  flag <- ifelse(abs(vol_cors[nm]) > 0.50, " *** confounded", "")
  cat(sprintf("    %s: r = %.2f%s\n", nm, vol_cors[nm], flag))
}

# ============================================================================
# 2.2 VIF ANALYSIS FOR REGRESSION
# ============================================================================

cat("\n------------------------------------------------------------\n")
cat("2.2 VIF ANALYSIS (VARIANCE INFLATION FACTORS)\n")
cat("------------------------------------------------------------\n\n")

# Full model VIF
full_vif_data <- landscape_data %>%
  filter(!is.na(isolation_index), !is.na(crowding_index), !is.na(relative_size))

cat("Full Model (all candidate predictors):\n")
cat("    (VIF > 5 indicates problematic multicollinearity)\n")

if (nrow(full_vif_data) > 20) {
  m_full_vif <- glm.nb(total_cafi ~ log(volume) + n_neighbors +
                        total_neighbor_volume + mean_neighbor_dist +
                        isolation_index + crowding_index + relative_size + site,
                       data = full_vif_data)

  vif_full <- car::vif(m_full_vif)
  for (nm in rownames(vif_full)) {
    vif_val <- vif_full[nm, "GVIF^(1/(2*Df))"]
    flag <- ifelse(vif_val > 5, " *** EXCLUDE", ifelse(vif_val > 3, " * high", ""))
    cat(sprintf("    %s: %.2f%s\n", nm, vif_val, flag))
  }
}

# Reduced model VIF
cat("\nReduced Model (non-redundant subset):\n")

m_reduced_vif <- glm.nb(total_cafi ~ log(volume) + n_neighbors +
                         log(total_neighbor_volume + 1) + mean_neighbor_dist + site,
                        data = landscape_data)

vif_reduced <- car::vif(m_reduced_vif)
for (nm in rownames(vif_reduced)) {
  vif_val <- vif_reduced[nm, "GVIF^(1/(2*Df))"]
  cat(sprintf("    %s: %.2f\n", nm, vif_val))
}

# ============================================================================
# 2.3 PREDICTOR SELECTION SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("PREDICTOR SELECTION FOR LANDSCAPE EFFECTS ANALYSIS\n")
cat("============================================================\n\n")

cat("Due to collinearity, we recommend 4 NON-REDUNDANT predictors:\n\n")

cat("SELECTED (low VIF, independent of focal volume):\n")
cat("  1. log(volume)             - Focal coral size (PRIMARY)\n")
cat("  2. n_neighbors             - Neighbor count\n")
cat("  3. total_neighbor_volume   - Neighborhood habitat\n")
cat("  4. mean_neighbor_dist      - Isolation\n\n")

cat("EXCLUDED (collinear, redundant, or confounded):\n")
cat("  - isolation_index: confounded with volume\n")
cat("  - crowding_index: collinear with total_neighbor_volume\n")
cat("  - relative_size: redundant interpretation\n\n")

# ############################################################################
#                    PART 3: EXTENDED ANALYSES
# ############################################################################

cat("############################################################\n")
cat("PART 3: EXTENDED ANALYSES\n")
cat("############################################################\n\n")

# ============================================================================
# 3.1 SIZE-ISOLATION RELATIONSHIP
# ============================================================================

cat("------------------------------------------------------------\n")
cat("3.1 SIZE-ISOLATION RELATIONSHIP\n")
cat("------------------------------------------------------------\n\n")

cat("Question: Are large corals more or less isolated than small corals?\n\n")

# Corals with neighbors
neigh_data <- landscape_data %>% filter(n_neighbors > 0)

# Correlation: volume vs distance
cor_vol_dist <- cor.test(neigh_data$log_volume, neigh_data$mean_neighbor_dist,
                          method = "spearman")
cat("Volume vs Mean Neighbor Distance (n =", nrow(neigh_data), "):\n")
cat("    Spearman rho =", round(cor_vol_dist$estimate, 3), "\n")
cat("    p-value =", format.pval(cor_vol_dist$p.value, 3), "\n\n")

# Correlation: volume vs neighbor count
cor_vol_n <- cor.test(landscape_data$log_volume, landscape_data$n_neighbors,
                       method = "spearman")
cat("Volume vs Number of Neighbors (n =", nrow(landscape_data), "):\n")
cat("    Spearman rho =", round(cor_vol_n$estimate, 3), "\n")
cat("    p-value =", format.pval(cor_vol_n$p.value, 3), "\n")
cat("    Interpretation:", ifelse(cor_vol_n$estimate < 0,
                                   "Larger corals have FEWER neighbors",
                                   "Larger corals have MORE neighbors"), "\n\n")

# Isolation by size tertile
isolation_by_size <- landscape_data %>%
  group_by(size_tertile) %>%
  summarise(
    n = n(),
    n_isolated = sum(n_neighbors == 0),
    pct_isolated = mean(n_neighbors == 0) * 100,
    mean_neighbors = mean(n_neighbors),
    mean_dist = mean(mean_neighbor_dist, na.rm = TRUE),
    .groups = "drop"
  )

cat("Isolation by Size Tertile:\n")
print(isolation_by_size)

# Chi-squared test
size_iso_table <- table(landscape_data$size_tertile, landscape_data$n_neighbors == 0)
chisq_iso <- chisq.test(size_iso_table)
cat("\n    Chi-squared test: χ² =", round(chisq_iso$statistic, 2),
    ", p =", format.pval(chisq_iso$p.value, 3), "\n\n")

# Figures
p_ext_1a <- ggplot(landscape_data, aes(x = volume, y = n_neighbors, color = site)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_y_sqrt() +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = expression("Coral Volume (cm"^3*")"),
       y = "Number of Neighbors",
       title = "A. Size vs Neighborhood Density",
       subtitle = paste0("ρ = ", round(cor_vol_n$estimate, 2),
                        ", p = ", format.pval(cor_vol_n$p.value, 2))) +
  theme(legend.position = c(0.85, 0.85))

p_ext_1b <- neigh_data %>%
  ggplot(aes(x = volume, y = mean_neighbor_dist, color = site)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = expression("Coral Volume (cm"^3*")"),
       y = "Mean Distance to Neighbors (cm)",
       title = "B. Size vs Neighbor Distance",
       subtitle = paste0("ρ = ", round(cor_vol_dist$estimate, 2),
                        ", p = ", format.pval(cor_vol_dist$p.value, 2))) +
  theme(legend.position = "none")

p_ext_1c <- landscape_data %>%
  ggplot(aes(x = size_tertile, fill = n_neighbors == 0)) +
  geom_bar(position = "fill", alpha = 0.8) +
  scale_fill_manual(values = c("FALSE" = "#56B4E9", "TRUE" = "#E69F00"),
                    labels = c("Has neighbors", "Isolated"),
                    name = "") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Coral Size Category",
       y = "Proportion",
       title = "C. Isolation Rate by Size",
       subtitle = paste0("χ² = ", round(chisq_iso$statistic, 1),
                        ", p = ", format.pval(chisq_iso$p.value, 2))) +
  theme(legend.position = "bottom")

# ============================================================================
# 3.2 SITE-SPECIFIC LANDSCAPE PATTERNS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("3.2 SITE-SPECIFIC LANDSCAPE PATTERNS\n")
cat("------------------------------------------------------------\n\n")

# Summary by site
site_landscape <- landscape_data %>%
  group_by(site) %>%
  summarise(
    n = n(),
    mean_volume = mean(volume),
    median_volume = median(volume),
    cv_volume = sd(volume) / mean(volume) * 100,
    pct_isolated = mean(n_neighbors == 0) * 100,
    mean_neighbors = mean(n_neighbors),
    mean_neighbor_dist = mean(mean_neighbor_dist, na.rm = TRUE),
    mean_total_neigh_vol = mean(total_neighbor_volume, na.rm = TRUE),
    .groups = "drop"
  )

cat("Landscape Structure by Site:\n")
print(site_landscape)
cat("\n")

# Tests
kw_vol <- kruskal.test(volume ~ site, data = landscape_data)
cat("Size ~ Site: H =", round(kw_vol$statistic, 2), ", p =", format.pval(kw_vol$p.value, 3), "\n")

kw_neigh <- kruskal.test(n_neighbors ~ site, data = landscape_data)
cat("Neighbors ~ Site: H =", round(kw_neigh$statistic, 2), ", p =", format.pval(kw_neigh$p.value, 3), "\n")

iso_site_table <- table(landscape_data$site, landscape_data$n_neighbors == 0)
chisq_site <- chisq.test(iso_site_table)
cat("Isolation ~ Site: χ² =", round(chisq_site$statistic, 2), ", p =", format.pval(chisq_site$p.value, 3), "\n\n")

# Figures
p_ext_2a <- landscape_data %>%
  group_by(site) %>%
  summarise(
    pct_isolated = mean(n_neighbors == 0) * 100,
    pct_sparse = mean(n_neighbors > 0 & n_neighbors <= 3) * 100,
    pct_moderate = mean(n_neighbors > 3 & n_neighbors <= 10) * 100,
    pct_dense = mean(n_neighbors > 10) * 100,
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with("pct"), names_to = "category", values_to = "pct") %>%
  mutate(category = factor(gsub("pct_", "", category),
                           levels = c("isolated", "sparse", "moderate", "dense"))) %>%
  ggplot(aes(x = site, y = pct, fill = category)) +
  geom_col(position = "stack", alpha = 0.8) +
  scale_fill_viridis_d(option = "D", name = "Density",
                       labels = c("Isolated", "Sparse", "Moderate", "Dense")) +
  labs(x = "Site", y = "Percentage",
       title = "D. Density Categories by Site") +
  theme(legend.position = "right")

p_ext_2b <- landscape_data %>%
  ggplot(aes(x = site, y = n_neighbors, fill = site)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = SITE_COLORS) +
  scale_y_sqrt() +
  labs(x = "Site", y = "Number of Neighbors (sqrt)",
       title = "E. Neighborhood Density by Site") +
  theme(legend.position = "none")

# ============================================================================
# 3.3 NEIGHBORHOOD SIZE STRUCTURE
# ============================================================================

cat("------------------------------------------------------------\n")
cat("3.3 NEIGHBORHOOD SIZE STRUCTURE\n")
cat("------------------------------------------------------------\n\n")

# Relative position summary
relative_size_summary <- landscape_data %>%
  filter(n_neighbors > 0) %>%
  group_by(relative_position) %>%
  summarise(
    n = n(),
    pct = n() / sum(landscape_data$n_neighbors > 0) * 100,
    mean_cafi = mean(total_cafi),
    mean_richness = mean(otu_richness),
    .groups = "drop"
  )

cat("Relative Position in Neighborhood:\n")
print(relative_size_summary)
cat("\n")

# ANOVA
aov_position <- aov(total_cafi ~ relative_position,
                     data = landscape_data %>% filter(n_neighbors > 0))
cat("CAFI ~ Position: F =", round(summary(aov_position)[[1]]["relative_position", "F value"], 2),
    ", p =", format.pval(summary(aov_position)[[1]]["relative_position", "Pr(>F)"], 3), "\n\n")

# Size ratio
size_ratio_data <- landscape_data %>%
  filter(n_neighbors > 0) %>%
  mutate(size_ratio = volume / mean_neighbor_volume)

cat("Size Ratio (Focal/Neighbor Mean):\n")
cat("    Median:", round(median(size_ratio_data$size_ratio), 2), "\n")
cat("    % focal > neighbors:", round(mean(size_ratio_data$size_ratio > 1) * 100, 1), "%\n\n")

# Figures
p_ext_3a <- size_ratio_data %>%
  ggplot(aes(x = size_ratio, y = total_cafi, color = site)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_x_log10() +
  scale_y_sqrt() +
  scale_color_manual(values = SITE_COLORS) +
  labs(x = "Size Ratio (Focal / Mean Neighbor)",
       y = "CAFI Abundance (sqrt)",
       title = "F. Relative Size Effect") +
  theme(legend.position = c(0.85, 0.85))

p_ext_3b <- landscape_data %>%
  filter(relative_position != "Solitary") %>%
  ggplot(aes(x = relative_position, y = total_cafi, fill = relative_position)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  scale_fill_viridis_d(option = "C") +
  scale_y_sqrt() +
  labs(x = "Position in Neighborhood",
       y = "CAFI Abundance (sqrt)",
       title = "G. CAFI by Relative Position") +
  theme(legend.position = "none")

# ============================================================================
# 3.4 LANDSCAPE TYPOLOGY (CLUSTERING)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("3.4 LANDSCAPE TYPOLOGY (CLUSTERING)\n")
cat("------------------------------------------------------------\n\n")

# Prepare clustering data
cluster_data <- landscape_data %>%
  filter(n_neighbors > 0, !is.na(mean_neighbor_dist), !is.na(total_neighbor_volume)) %>%
  mutate(
    log_vol_scaled = scale(log_volume)[,1],
    n_neigh_scaled = scale(n_neighbors)[,1],
    dist_scaled = scale(mean_neighbor_dist)[,1],
    neigh_vol_scaled = scale(log10(total_neighbor_volume + 1))[,1]
  )

# K-means clustering
set.seed(42)
km_data <- cluster_data %>%
  dplyr::select(log_vol_scaled, n_neigh_scaled, dist_scaled, neigh_vol_scaled)

km_result <- kmeans(km_data, centers = 3, nstart = 25)
cluster_data$cluster <- factor(km_result$cluster)

# Summarize clusters
cluster_summary <- cluster_data %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    pct = n() / nrow(cluster_data) * 100,
    mean_volume = mean(volume),
    mean_neighbors = mean(n_neighbors),
    mean_dist = mean(mean_neighbor_dist),
    mean_neigh_vol = mean(total_neighbor_volume),
    mean_cafi = mean(total_cafi),
    mean_richness = mean(otu_richness),
    .groups = "drop"
  ) %>%
  mutate(
    type = case_when(
      mean_neighbors == max(mean_neighbors) ~ "Dense Patch",
      mean_dist == max(mean_dist) ~ "Distant Neighbors",
      TRUE ~ "Moderate Density"
    )
  )

cat("Landscape Types (k-means, k=3):\n")
print(cluster_summary %>% dplyr::select(cluster, type, n, pct, mean_neighbors, mean_dist, mean_cafi))
cat("\n")

aov_cluster <- aov(total_cafi ~ cluster, data = cluster_data)
cat("CAFI ~ Cluster: F =", round(summary(aov_cluster)[[1]]["cluster", "F value"], 2),
    ", p =", format.pval(summary(aov_cluster)[[1]]["cluster", "Pr(>F)"], 3), "\n\n")

# Figures
p_ext_4a <- cluster_data %>%
  ggplot(aes(x = n_neighbors, y = mean_neighbor_dist, color = cluster)) +
  geom_point(aes(size = volume), alpha = 0.7) +
  scale_color_viridis_d(option = "D", name = "Type") +
  scale_size_continuous(name = expression("Volume"), range = c(2, 8)) +
  labs(x = "Number of Neighbors",
       y = "Mean Distance to Neighbors (cm)",
       title = "H. Landscape Types (k-means)") +
  theme(legend.position = "right")

p_ext_4b <- cluster_data %>%
  ggplot(aes(x = cluster, y = total_cafi, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  scale_fill_viridis_d(option = "D") +
  scale_y_sqrt() +
  labs(x = "Landscape Type",
       y = "CAFI Abundance (sqrt)",
       title = "I. CAFI by Landscape Type") +
  theme(legend.position = "none")

# ############################################################################
#                    SAVE OUTPUTS
# ############################################################################

cat("############################################################\n")
cat("SAVING OUTPUTS\n")
cat("############################################################\n\n")

# ============================================================================
# BASE CHARACTERIZATION FIGURE
# ============================================================================

p_landscape <- (p_vol_hist + p_vol_box) /
  (p_neigh_count + p_neigh_dist) /
  (p_neigh_vol + plot_spacer()) +
  plot_annotation(
    title = "Pocillopora Landscape Structure",
    subtitle = paste0("N = ", nrow(landscape_data), " corals with neighborhood surveys (5m radius), Mo'orea"),
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12))
  )

ggsave(file.path(FIG_DIR, "landscape_characterization.png"), p_landscape,
       width = 12, height = 12, dpi = 300, bg = "white")
cat("  Saved: landscape_characterization.png\n")

# ============================================================================
# CORRELATION HEATMAP
# ============================================================================

cor_for_plot <- cor_vars %>%
  dplyr::select(volume, n_neighbors, total_neighbor_volume, mean_neighbor_dist,
         isolation_index, crowding_index, relative_size)

cor_mat_plot <- cor(cor_for_plot, use = "pairwise.complete.obs")

cor_long <- cor_mat_plot %>%
  as.data.frame() %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "r") %>%
  mutate(
    var1 = factor(var1, levels = colnames(cor_mat_plot)),
    var2 = factor(var2, levels = rev(colnames(cor_mat_plot)))
  )

p_cor_heatmap <- ggplot(cor_long, aes(x = var1, y = var2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Correlation") +
  labs(x = "", y = "",
       title = "Correlation Matrix of Landscape Predictors",
       subtitle = "High correlations (|r| > 0.70) indicate redundancy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave(file.path(FIG_DIR, "predictor_correlations.png"), p_cor_heatmap,
       width = 9, height = 8, dpi = 300, bg = "white")
cat("  Saved: predictor_correlations.png\n")

# ============================================================================
# EXTENDED ANALYSIS FIGURE
# ============================================================================

p_extended <- (p_ext_1a + p_ext_1b + p_ext_1c) /
  (p_ext_2a + p_ext_2b + p_ext_3a) /
  (p_ext_3b + p_ext_4a + p_ext_4b) +
  plot_annotation(
    title = "Extended Landscape Characterization",
    subtitle = paste0("N = ", nrow(landscape_data), " corals with neighborhood surveys; ",
                     sum(landscape_data$n_neighbors == 0), " isolated, ",
                     sum(landscape_data$has_neighbors), " with neighbors"),
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12))
  )

ggsave(file.path(FIG_DIR, "landscape_extended_analysis.png"), p_extended,
       width = 15, height = 15, dpi = 300, bg = "white")
cat("  Saved: landscape_extended_analysis.png\n")

# ============================================================================
# TABLES
# ============================================================================

# Base landscape summary
landscape_summary <- tibble(
  metric = c("Total corals", "Sites", "Volume range (cm³)", "Volume CV (%)",
             "% Isolated (0 neighbors)", "Mean neighbors", "Max neighbors",
             "Mean neighbor distance (m)", "Mean neighborhood volume (cm³)"),
  value = c(nrow(landscape_data),
            paste(unique(landscape_data$site), collapse = ", "),
            paste0(round(vol_summary$min), " - ", round(vol_summary$max)),
            round(vol_summary$cv, 1),
            round(neigh_count_summary$pct_isolated, 1),
            round(neigh_count_summary$mean_neighbors, 1),
            neigh_count_summary$max_neighbors,
            round(dist_summary$mean_m, 2),
            round(neigh_vol_summary$mean, 0))
)
save_table(landscape_summary, "landscape_structure_summary")

# Size-isolation results
size_isolation_results <- tibble(
  comparison = c("Volume vs N neighbors", "Volume vs Distance", "Isolation vs Size"),
  test = c("Spearman correlation", "Spearman correlation", "Chi-squared"),
  statistic = c(cor_vol_n$estimate, cor_vol_dist$estimate, chisq_iso$statistic),
  p_value = c(cor_vol_n$p.value, cor_vol_dist$p.value, chisq_iso$p.value)
)
save_table(size_isolation_results, "size_isolation_tests")

# Site landscape summary
save_table(site_landscape, "site_landscape_summary")

# Relative position summary
save_table(relative_size_summary, "relative_position_summary")

# Cluster summary
save_table(cluster_summary, "landscape_cluster_summary")

cat("  Tables saved to:", PATHS$tables, "\n")

# ============================================================================
# OBJECTS
# ============================================================================

# Selected predictors
selected_predictors <- c("log_volume", "n_neighbors", "total_neighbor_volume", "mean_neighbor_dist")
save_object(selected_predictors, "landscape_selected_predictors")

# Extended results
extended_results <- list(
  size_isolation = list(
    vol_vs_neighbors = cor_vol_n,
    vol_vs_distance = cor_vol_dist,
    isolation_by_size = isolation_by_size,
    chisq_test = chisq_iso
  ),
  site_patterns = list(
    summary = site_landscape,
    size_test = kw_vol,
    neighbor_test = kw_neigh,
    isolation_test = chisq_site
  ),
  neighborhood_structure = list(
    relative_position = relative_size_summary,
    position_anova = aov_position
  ),
  landscape_types = list(
    cluster_summary = cluster_summary,
    cluster_anova = aov_cluster,
    cluster_assignments = cluster_data %>% dplyr::select(coral_id, cluster)
  )
)
save_object(extended_results, "landscape_extended_results")

cat("  Objects saved to:", PATHS$objects, "\n")

# ############################################################################
#                    SUMMARY
# ############################################################################

cat("\n")
cat("============================================================\n")
cat("    LANDSCAPE CHARACTERIZATION COMPLETE\n")
cat("============================================================\n\n")

cat("NOTE: All analyses use ONLY the", nrow(landscape_data), "corals with\n")
cat("      5m neighborhood surveys. The 51 'size' corals were NOT\n")
cat("      surveyed for neighbors.\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. CORAL SIZE:\n")
cat("   - Range:", round(vol_summary$log_range, 1), "orders of magnitude\n")
cat("   - CV =", round(vol_summary$cv, 0), "%\n\n")

cat("2. NEIGHBORHOOD STRUCTURE:\n")
cat("   -", round(neigh_count_summary$pct_isolated, 0), "% isolated (0 neighbors)\n")
cat("   -", round(neigh_count_summary$pct_dense, 0), "% in dense patches (>10 neighbors)\n")
cat("   - Mean neighbors:", round(neigh_count_summary$mean_neighbors, 1), "\n\n")

cat("3. SIZE-ISOLATION:\n")
cat("   - Volume vs N neighbors: ρ =", round(cor_vol_n$estimate, 2), "\n")
cat("   - Larger corals have", ifelse(cor_vol_n$estimate < 0, "FEWER", "MORE"), "neighbors\n\n")

cat("4. SITE PATTERNS:\n")
for (i in 1:nrow(site_landscape)) {
  cat("   -", site_landscape$site[i], ":",
      round(site_landscape$pct_isolated[i], 0), "% isolated,",
      round(site_landscape$mean_neighbors[i], 1), "mean neighbors\n")
}
cat("\n")

cat("5. LANDSCAPE TYPES:\n")
for (i in 1:nrow(cluster_summary)) {
  cat("   - Type", cluster_summary$cluster[i], "(", cluster_summary$type[i], "):",
      round(cluster_summary$pct[i]), "% of corals\n")
}

cat("\nOUTPUT FILES:\n")
cat("  Figures: landscape_characterization.png\n")
cat("           predictor_correlations.png\n")
cat("           landscape_extended_analysis.png\n")
cat("  Tables:  landscape_structure_summary.csv\n")
cat("           size_isolation_tests.csv\n")
cat("           site_landscape_summary.csv\n")
cat("           relative_position_summary.csv\n")
cat("           landscape_cluster_summary.csv\n")
cat("  Objects: landscape_selected_predictors.rds\n")
cat("           landscape_extended_results.rds\n\n")

cat("Next: source('scripts/04_landscape_effects.R')\n\n")
