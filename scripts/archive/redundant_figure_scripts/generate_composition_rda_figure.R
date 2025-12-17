#!/usr/bin/env Rscript
# =============================================================================
# Figure: RDA analysis of community composition vs. ALL coral & neighborhood traits
# Shows how size, proximity, neighbor count, and neighbor volume predict composition
# =============================================================================

library(tidyverse)
library(vegan)
library(patchwork)

# Publication theme
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(face = "bold", size = base_size - 1),
      plot.title = element_text(face = "bold"),
      legend.background = element_rect(fill = "white", color = NA)
    )
}

cat("Loading and preparing data...\n")

# Load data
cafi_data <- read_csv("data/survey_cafi_data_w_taxonomy_summer2019_v5.csv",
                      show_col_types = FALSE)
coral_data <- read_csv("data/survey_coral_characteristics_merged_v2.csv",
                       show_col_types = FALSE)

# Standardize type names
cafi_data <- cafi_data %>%
  mutate(type_clean = case_when(
    type == "crab" ~ "Crab",
    type == "shrimp" ~ "Shrimp",
    type == "fish" ~ "Fish",
    type == "hermit" ~ "Hermit",
    type == "snail" ~ "Snail",
    type == "echinoderm" ~ "Echinoderm",
    TRUE ~ "Other"
  ))

# Create community matrix (corals x functional groups) - counts
comm_matrix <- cafi_data %>%
  group_by(coral_id, type_clean) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = type_clean, values_from = count, values_fill = 0)

coral_ids <- comm_matrix$coral_id
comm_mat <- as.matrix(comm_matrix[, -1])
rownames(comm_mat) <- coral_ids

# Get environmental data - ALL neighborhood metrics
env_data <- coral_data %>%
  select(coral_id, volume_field, number_of_neighbors, mean_neighbor_distance,
         combined_total_volume_of_neighbors, site) %>%
  filter(coral_id %in% coral_ids) %>%
  mutate(
    log_vol = log10(volume_field),
    log_neighbor_vol = log10(combined_total_volume_of_neighbors + 1),
    # Scale mean_neighbor_distance for better interpretation
    neighbor_proximity = 1 / (mean_neighbor_distance + 1)  # Higher = closer neighbors
  )

# Filter to corals with complete data
complete_ids <- env_data %>%
  filter(!is.na(volume_field),
         !is.na(number_of_neighbors),
         !is.na(mean_neighbor_distance),
         !is.na(combined_total_volume_of_neighbors)) %>%
  pull(coral_id)

comm_mat_filtered <- comm_mat[rownames(comm_mat) %in% complete_ids, ]
env_filtered <- env_data %>% filter(coral_id %in% complete_ids)

# Ensure matching order
env_filtered <- env_filtered[match(rownames(comm_mat_filtered),
                                    env_filtered$coral_id), ]

cat("Corals with complete data:", nrow(comm_mat_filtered), "\n")

# Print variable correlations
cat("\n=== Predictor Correlations ===\n")
pred_cors <- env_filtered %>%
  select(log_vol, number_of_neighbors, mean_neighbor_distance, log_neighbor_vol) %>%
  cor(use = "complete.obs")
print(round(pred_cors, 3))

# =============================================================================
# Run db-RDA (distance-based RDA) - constrained ordination with ALL predictors
# Using Jaccard distance on presence/absence to emphasize composition
# =============================================================================

cat("\nRunning db-RDA analysis with ALL predictors...\n")
cat("Using Jaccard dissimilarity on presence/absence data\n")

# Convert to presence/absence
comm_pa <- (comm_mat_filtered > 0) * 1

# Run db-RDA using capscale with Jaccard distance
# capscale is the proper function for distance-based RDA
rda_full <- capscale(comm_pa ~ log_vol + number_of_neighbors +
                       mean_neighbor_distance + log_neighbor_vol,
                     data = env_filtered,
                     distance = "jaccard")

# Test significance
cat("\n=== FULL MODEL (all 4 predictors) ===\n")
rda_anova_full <- anova(rda_full, permutations = 999)
cat("Overall model significance:\n")
print(rda_anova_full)

rda_terms <- anova(rda_full, by = "terms", permutations = 999)
cat("\nSequential term significance (Type I SS):\n")
print(rda_terms)

rda_margin <- anova(rda_full, by = "margin", permutations = 999)
cat("\nMarginal term significance (Type III SS - unique contribution):\n")
print(rda_margin)

# Model selection using forward selection
cat("\n=== FORWARD SELECTION ===\n")
rda_null <- capscale(comm_pa ~ 1, data = env_filtered, distance = "jaccard")
rda_forward <- ordistep(rda_null,
                         scope = formula(rda_full),
                         direction = "forward",
                         permutations = 999,
                         trace = TRUE)
cat("\nBest model from forward selection:\n")
print(rda_forward$call)

# Also try backward elimination
cat("\n=== BACKWARD ELIMINATION ===\n")
rda_backward <- ordistep(rda_full,
                          direction = "backward",
                          permutations = 999,
                          trace = TRUE)
cat("\nBest model from backward elimination:\n")
print(rda_backward$call)

# Use the FULL model for visualization to show all predictors
rda_result <- rda_full

# Summary
cat("\n=== FULL MODEL SUMMARY ===\n")
print(summary(rda_full))

# =============================================================================
# Variance partitioning - separate size vs neighborhood effects
# Using Jaccard distance for consistency
# =============================================================================

cat("\n=== VARIANCE PARTITIONING ===\n")

# Compute Jaccard distance matrix
jacc_dist <- vegdist(comm_pa, method = "jaccard", binary = TRUE)

# Two-way partition: Size vs All Neighborhood
var_part_2way <- varpart(jacc_dist,
                          ~ log_vol,
                          ~ number_of_neighbors + mean_neighbor_distance + log_neighbor_vol,
                          data = env_filtered)
cat("\n--- Size vs Neighborhood (combined) ---\n")
print(var_part_2way)

# Four-way partition: each predictor separately
var_part_4way <- varpart(jacc_dist,
                          ~ log_vol,
                          ~ number_of_neighbors,
                          ~ mean_neighbor_distance,
                          ~ log_neighbor_vol,
                          data = env_filtered)
cat("\n--- All 4 predictors separately ---\n")
print(var_part_4way)

# =============================================================================
# Extract scores for plotting
# =============================================================================

# Site scores - capscale uses CAP1, CAP2 instead of RDA1, RDA2
site_scores <- as.data.frame(scores(rda_result, display = "sites", scaling = 2))
# Rename columns to RDA1, RDA2 for consistency
if("CAP1" %in% colnames(site_scores)) {
  colnames(site_scores) <- gsub("CAP", "RDA", colnames(site_scores))
}
site_scores$coral_id <- rownames(site_scores)
site_scores <- site_scores %>% left_join(env_filtered, by = "coral_id")

# Species scores
species_scores <- as.data.frame(scores(rda_result, display = "species", scaling = 2))
if("CAP1" %in% colnames(species_scores)) {
  colnames(species_scores) <- gsub("CAP", "RDA", colnames(species_scores))
}
species_scores$taxon <- rownames(species_scores)

# Environmental vectors (biplot scores)
env_scores <- as.data.frame(scores(rda_result, display = "bp", scaling = 2))
if("CAP1" %in% colnames(env_scores)) {
  colnames(env_scores) <- gsub("CAP", "RDA", colnames(env_scores))
}
env_scores$variable <- rownames(env_scores)
env_scores$label <- c("log(Volume)", "# Neighbors", "Mean Distance", "log(Neighbor Vol)")

# Variance explained
eig <- rda_result$CCA$eig
var_rda1 <- eig[1] / sum(eig) * 100
var_rda2 <- eig[2] / sum(eig) * 100
var_rda3 <- if(length(eig) >= 3) eig[3] / sum(eig) * 100 else 0
var_rda4 <- if(length(eig) >= 4) eig[4] / sum(eig) * 100 else 0
total_constrained <- sum(eig) / rda_result$tot.chi * 100

cat(sprintf("\nTotal constrained variance: %.1f%%\n", total_constrained))
cat(sprintf("RDA1: %.1f%%, RDA2: %.1f%%, RDA3: %.1f%%, RDA4: %.1f%%\n",
            var_rda1, var_rda2, var_rda3, var_rda4))

# Colors
taxon_colors <- c(
  "Crab" = "#E74C3C",
  "Shrimp" = "#3498DB",
  "Fish" = "#2ECC71",
  "Hermit" = "#F39C12",
  "Snail" = "#9B59B6",
  "Echinoderm" = "#1ABC9C",
  "Other" = "#95A5A6"
)

# =============================================================================
# Panel A: RDA biplot with all environmental vectors
# =============================================================================

cat("\nCreating Panel A: RDA biplot with all predictors...\n")

# Scale for arrows
sp_scale <- 1.5
env_scale <- 0.35

p_a <- ggplot() +
  # Site points colored by size
  geom_point(data = site_scores,
             aes(x = RDA1, y = RDA2, color = log_vol, size = volume_field),
             alpha = 0.7) +
  # Species arrows
  geom_segment(data = species_scores,
               aes(x = 0, y = 0, xend = RDA1 * sp_scale, yend = RDA2 * sp_scale),
               arrow = arrow(length = unit(0.15, "cm")),
               color = "gray30", linewidth = 0.7) +
  geom_text(data = species_scores,
            aes(x = RDA1 * sp_scale * 1.12, y = RDA2 * sp_scale * 1.12,
                label = taxon),
            size = 2.8, fontface = "bold", color = "gray20") +
  # Environmental vectors (all 4)
  geom_segment(data = env_scores,
               aes(x = 0, y = 0, xend = RDA1 * env_scale, yend = RDA2 * env_scale),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "#C0392B", linewidth = 1.1) +
  geom_text(data = env_scores,
            aes(x = RDA1 * env_scale * 1.15, y = RDA2 * env_scale * 1.15, label = label),
            size = 3, fontface = "bold", color = "#C0392B") +
  scale_color_viridis_c(name = expression("log"[10]*"(Vol)"),
                        option = "plasma", end = 0.9) +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  labs(title = "A. RDA biplot: composition vs. coral & neighborhood traits",
       subtitle = sprintf("Model explains %.1f%% of variance (p = %.3f)",
                         total_constrained, rda_anova_full$`Pr(>F)`[1]),
       x = sprintf("RDA1 (%.1f%%)", var_rda1),
       y = sprintf("RDA2 (%.1f%%)", var_rda2)) +
  theme_pub(base_size = 11) +
  theme(plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        legend.position = c(0.12, 0.82),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# =============================================================================
# Panel B: Predictor importance (marginal R²)
# =============================================================================

cat("Creating Panel B: Predictor importance...\n")

# Get marginal (unique) contributions from variance partitioning
# Use the 2-way partition for cleaner visualization
adj_r2_size <- var_part_2way$part$indfract$Adj.R.squared[1]
adj_r2_shared <- var_part_2way$part$indfract$Adj.R.squared[3]
adj_r2_neighborhood <- var_part_2way$part$indfract$Adj.R.squared[2]

# Also get individual predictor marginal significance from rda_margin
# Note: capscale uses SumOfSqs instead of Variance
pred_importance <- data.frame(
  predictor = c("Coral Size", "# Neighbors", "Mean Distance", "Neighbor Volume"),
  F_value = rda_margin$F[1:4],
  p_value = rda_margin$`Pr(>F)`[1:4],
  variance = rda_margin$SumOfSqs[1:4]
)
pred_importance$prop_variance <- pred_importance$variance / sum(pred_importance$variance) * 100
pred_importance$sig <- ifelse(pred_importance$p_value < 0.001, "***",
                               ifelse(pred_importance$p_value < 0.01, "**",
                                      ifelse(pred_importance$p_value < 0.05, "*", "")))
pred_importance <- pred_importance %>%
  arrange(desc(F_value)) %>%
  mutate(predictor = factor(predictor, levels = predictor))

p_b <- ggplot(pred_importance, aes(x = predictor, y = F_value)) +
  geom_col(aes(fill = p_value < 0.05), width = 0.7, alpha = 0.85) +
  geom_text(aes(label = sig), vjust = -0.3, size = 5, fontface = "bold") +
  geom_text(aes(label = sprintf("p=%.3f", p_value)), vjust = 1.5, size = 2.8,
            color = "white", fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#27AE60", "FALSE" = "#BDC3C7"),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "B. Predictor importance (marginal F-test)",
       subtitle = "Unique contribution controlling for other predictors",
       x = NULL,
       y = "F-statistic") +
  theme_pub(base_size = 11) +
  theme(plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        axis.text.x = element_text(angle = 25, hjust = 1))

# =============================================================================
# Panel C: Variance partitioning - Size vs Neighborhood
# =============================================================================

cat("Creating Panel C: Variance partitioning...\n")

var_df <- data.frame(
  component = c("Size\n(unique)", "Shared", "Neighborhood\n(unique)", "Residual"),
  value = c(max(0, adj_r2_size),
            max(0, adj_r2_shared),
            max(0, adj_r2_neighborhood),
            1 - max(0, adj_r2_size) - max(0, adj_r2_shared) - max(0, adj_r2_neighborhood)) * 100,
  fill_col = c("#3498DB", "#9B59B6", "#E74C3C", "#BDC3C7")
)
var_df$component <- factor(var_df$component, levels = var_df$component)

p_c <- ggplot(var_df, aes(x = component, y = value, fill = fill_col)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%", value)), vjust = -0.3,
            size = 3.5, fontface = "bold") +
  scale_fill_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     limits = c(0, 100)) +
  labs(title = "C. Variance partitioning",
       subtitle = "Size vs. all neighborhood metrics combined",
       x = NULL,
       y = "Variance explained (%)") +
  theme_pub(base_size = 11) +
  theme(plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"))

# =============================================================================
# Panel D: Species loadings on RDA1
# =============================================================================

cat("Creating Panel D: Species loadings...\n")

species_scores_plot <- species_scores %>%
  arrange(RDA1) %>%
  mutate(taxon = factor(taxon, levels = taxon))

p_d <- ggplot(species_scores_plot, aes(x = taxon, y = RDA1, fill = taxon)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = taxon_colors, guide = "none") +
  coord_flip() +
  labs(title = "D. Taxon associations with primary gradient",
       subtitle = "RDA1 dominated by coral size effect",
       x = NULL,
       y = "RDA1 score") +
  theme_pub(base_size = 11) +
  theme(plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        axis.text.y = element_text(face = "bold"))

# =============================================================================
# Panel E: RDA1 vs each predictor
# =============================================================================

cat("Creating Panel E: RDA1 vs predictors...\n")

# Correlation of RDA1 with each predictor
cors_rda1 <- data.frame(
  predictor = c("log_vol", "number_of_neighbors", "mean_neighbor_distance", "log_neighbor_vol"),
  label = c("log(Volume)", "# Neighbors", "Mean Distance", "log(Neighbor Vol)")
)
cors_rda1$r <- sapply(cors_rda1$predictor, function(x) {
  cor(site_scores$RDA1, site_scores[[x]])
})
cors_rda1$p <- sapply(cors_rda1$predictor, function(x) {
  cor.test(site_scores$RDA1, site_scores[[x]])$p.value
})

# Long format for faceted plot
site_long <- site_scores %>%
  select(coral_id, RDA1, log_vol, number_of_neighbors,
         mean_neighbor_distance, log_neighbor_vol) %>%
  pivot_longer(cols = c(log_vol, number_of_neighbors,
                        mean_neighbor_distance, log_neighbor_vol),
               names_to = "predictor",
               values_to = "value") %>%
  left_join(cors_rda1, by = "predictor") %>%
  mutate(label = factor(label, levels = c("log(Volume)", "# Neighbors",
                                          "Mean Distance", "log(Neighbor Vol)")))

p_e <- ggplot(site_long, aes(x = value, y = RDA1)) +
  geom_point(alpha = 0.5, size = 1.5, color = "#3498DB") +
  geom_smooth(method = "lm", se = TRUE, color = "#2C3E50", linewidth = 0.8) +
  geom_text(data = cors_rda1 %>%
              mutate(label = factor(label, levels = c("log(Volume)", "# Neighbors",
                                                      "Mean Distance", "log(Neighbor Vol)"))),
            aes(x = -Inf, y = Inf,
                label = sprintf("r=%.2f%s", r,
                               ifelse(p < 0.001, "***",
                                      ifelse(p < 0.01, "**",
                                             ifelse(p < 0.05, "*", ""))))),
            hjust = -0.1, vjust = 1.5, size = 3, fontface = "italic",
            inherit.aes = FALSE) +
  facet_wrap(~ label, scales = "free_x", nrow = 1) +
  labs(title = "E. RDA1 correlation with each predictor",
       subtitle = "Coral size shows strongest association with compositional gradient",
       x = "Predictor value",
       y = "RDA1 score") +
  theme_pub(base_size = 10) +
  theme(plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        strip.text = element_text(size = 9))

# =============================================================================
# Panel F: Environmental vector angles
# =============================================================================

cat("Creating Panel F: Environmental vector directions...\n")

# Calculate vector lengths and angles
env_scores$length <- sqrt(env_scores$RDA1^2 + env_scores$RDA2^2)
env_scores$angle <- atan2(env_scores$RDA2, env_scores$RDA1) * 180 / pi

p_f <- ggplot(env_scores, aes(x = reorder(label, -length), y = length)) +
  geom_col(fill = "#C0392B", alpha = 0.8, width = 0.6) +
  geom_text(aes(label = sprintf("%.0f°", angle)), vjust = -0.3, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "F. Environmental vector strength",
       subtitle = "Length = effect magnitude; angle shown in degrees",
       x = NULL,
       y = "Vector length (RDA space)") +
  theme_pub(base_size = 11) +
  theme(plot.title = element_text(size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        axis.text.x = element_text(angle = 25, hjust = 1))

# =============================================================================
# Combine panels
# =============================================================================

cat("Combining panels...\n")

fig_rda <- (p_a | p_b | p_c) / (p_d | p_f) / p_e +
  plot_layout(heights = c(1.2, 1, 0.8)) +
  plot_annotation(
    title = "Figure X. Multivariate analysis: Coral size dominates, neighborhood effects minimal",
    subtitle = sprintf("db-RDA on Jaccard dissimilarity (presence/absence) | Variance explained: %.1f%% (p = %.3f)",
                       total_constrained, rda_anova_full$`Pr(>F)`[1]),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray30")
    )
  )

# Save
dir.create("output/figures/composition", showWarnings = FALSE, recursive = TRUE)
ggsave("output/figures/composition/rda_full_composition_analysis.png", fig_rda,
       width = 15, height = 13, dpi = 300, bg = "white")

cat("\n✓ Figure saved to: output/figures/composition/rda_full_composition_analysis.png\n")

# Simplified version for manuscript
fig_simple <- (p_a + theme(legend.position = "right")) / (p_b | p_c) / p_e +
  plot_layout(heights = c(1.3, 1, 0.8)) +
  plot_annotation(
    title = "Figure X. Coral size is the primary driver of CAFI composition",
    subtitle = sprintf("db-RDA on Jaccard dissimilarity | %.1f%% variance explained",
                       total_constrained),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray30")
    )
  )

ggsave("output/figures/manuscript/Figure8_rda_composition.png", fig_simple,
       width = 13, height = 11, dpi = 300, bg = "white")

cat("✓ Main text figure saved to: output/figures/manuscript/Figure8_rda_composition.png\n")

# =============================================================================
# Print comprehensive summary
# =============================================================================

cat("\n")
cat("================================================================================\n")
cat("                         COMPREHENSIVE SUMMARY                                  \n")
cat("================================================================================\n")

cat(sprintf("\nTotal variance explained by full model: %.1f%%\n", total_constrained))
cat(sprintf("Model significance: F = %.2f, p = %.4f\n",
            rda_anova_full$F[1], rda_anova_full$`Pr(>F)`[1]))

cat("\n--- Predictor Significance (Marginal/Type III) ---\n")
print(pred_importance %>% select(predictor, F_value, p_value, sig))

cat("\n--- Variance Partitioning (Size vs Neighborhood) ---\n")
cat(sprintf("  Unique to coral size: %.1f%%\n", max(0, adj_r2_size) * 100))
cat(sprintf("  Shared: %.1f%%\n", max(0, adj_r2_shared) * 100))
cat(sprintf("  Unique to neighborhood: %.1f%%\n", max(0, adj_r2_neighborhood) * 100))

cat("\n--- RDA1 Correlations with Predictors ---\n")
print(cors_rda1 %>% select(label, r, p) %>% arrange(desc(abs(r))))

cat("\n--- Species Associations with RDA1 ---\n")
print(species_scores_plot %>% select(taxon, RDA1) %>% arrange(desc(RDA1)))

cat("\n--- Environmental Vector Lengths (effect size in RDA space) ---\n")
print(env_scores %>% select(label, length, angle) %>% arrange(desc(length)))

cat("\n================================================================================\n")
cat("CONCLUSION: Coral size is the dominant predictor of CAFI composition.\n")
cat("Neighborhood metrics (# neighbors, proximity, neighbor volume) show\n")
cat("minimal unique contributions after accounting for coral size.\n")
cat("================================================================================\n")
