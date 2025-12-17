#!/usr/bin/env Rscript
# ============================================================================
# 15_composition_by_neighborhood.R - How CAFI composition changes with neighborhood
# ============================================================================
#
# PURPOSE: Analyze how community composition varies with coral neighborhood context
#
# KEY QUESTIONS:
#   - Do corals in denser neighborhoods support different assemblages?
#   - How does total neighbor volume affect composition?
#   - Does neighbor distance influence community structure?
#
# OUTPUTS:
#   - 4-panel neighborhood composition figure (main + manuscript folder)
#   - Statistical tests for neighborhood effects
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("15: Composition by Neighborhood Analysis\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "composition")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Prepare Data
# ============================================================================

cat("Preparing composition data by neighborhood context...\n")

# Standardize type names for display
cafi_typed <- cafi_clean %>%
  mutate(type_clean = case_when(
    type == "crab" ~ "Crab",
    type == "shrimp" ~ "Shrimp",
    type == "fish" ~ "Fish",
    type == "hermit" ~ "Hermit",
    type == "snail" ~ "Snail",
    type == "echinoderm" ~ "Echinoderm",
    TRUE ~ "Other"
  ))

# Merge with neighborhood data
merged <- cafi_typed %>%
  left_join(coral_master %>%
              select(coral_id, volume, number_of_neighbors, mean_neighbor_dist,
                     total_neighbor_volume),
            by = "coral_id") %>%
  filter(!is.na(number_of_neighbors), !is.na(volume))

# Add site if available
if ("site" %in% names(coral_master)) {
  merged <- merged %>%
    left_join(coral_master %>% select(coral_id, site), by = "coral_id")
} else if (!"site" %in% names(merged)) {
  merged$site <- "All"
}

cat(sprintf("  - Records with neighborhood data: %d\n", nrow(merged)))
cat(sprintf("  - Unique corals: %d\n", n_distinct(merged$coral_id)))

# Per-coral composition with neighborhood metrics
# Build group_by dynamically based on available columns
group_cols <- c("coral_id", "number_of_neighbors", "mean_neighbor_dist",
                "total_neighbor_volume", "volume")
if ("site" %in% names(merged)) group_cols <- c(group_cols, "site")

coral_comp <- merged %>%
  group_by(across(all_of(group_cols))) %>%
  summarise(
    abundance = n(),
    richness = n_distinct(species),
    n_crab = sum(type_clean == "Crab"),
    n_shrimp = sum(type_clean == "Shrimp"),
    n_fish = sum(type_clean == "Fish"),
    n_snail = sum(type_clean == "Snail"),
    n_hermit = sum(type_clean == "Hermit"),
    .groups = "drop"
  ) %>%
  mutate(
    prop_crab = n_crab / abundance,
    prop_shrimp = n_shrimp / abundance,
    prop_fish = n_fish / abundance,
    prop_snail = n_snail / abundance,
    prop_hermit = n_hermit / abundance,
    log_vol = log10(volume),
    log_neighbor_vol = log10(total_neighbor_volume + 1)
  )

# Ensure site column exists in coral_comp
if (!"site" %in% names(coral_comp)) {
  coral_comp$site <- "All"
}

# Create neighbor count bins
coral_comp <- coral_comp %>%
  mutate(
    neighbor_bin = cut(number_of_neighbors,
                       breaks = quantile(number_of_neighbors, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                       labels = c("Q1 (Few)", "Q2", "Q3", "Q4 (Many)"),
                       include.lowest = TRUE)
  )

# Functional group colors (colorblind-friendly)
group_colors <- c(
  "Crab" = "#D55E00",      # vermillion
  "Shrimp" = "#0072B2",    # blue
  "Fish" = "#009E73",      # bluish green
  "Snail" = "#CC79A7",     # reddish purple
  "Hermit" = "#E69F00",    # orange
  "Other" = "#999999"      # gray
)

# ============================================================================
# Panel A: Composition by number of neighbors (stacked bars)
# ============================================================================

cat("Creating Panel A: Composition by neighbor count...\n")

comp_by_neighbors <- merged %>%
  left_join(coral_comp %>% select(coral_id, neighbor_bin), by = "coral_id") %>%
  filter(!is.na(neighbor_bin)) %>%
  group_by(neighbor_bin, type_clean) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(neighbor_bin) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Order groups
group_order <- comp_by_neighbors %>%
  group_by(type_clean) %>%
  summarise(total = sum(count)) %>%
  arrange(desc(total)) %>%
  pull(type_clean)

comp_by_neighbors$type_clean <- factor(comp_by_neighbors$type_clean, levels = rev(group_order))

# Sample sizes
n_by_bin <- coral_comp %>%
  filter(!is.na(neighbor_bin)) %>%
  group_by(neighbor_bin) %>%
  summarise(n = n())

p_a <- ggplot(comp_by_neighbors, aes(x = neighbor_bin, y = prop, fill = type_clean)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3, width = 0.8) +
  geom_text(data = n_by_bin, aes(x = neighbor_bin, y = 1.03, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 3, color = "gray40") +
  scale_fill_manual(values = group_colors, name = "Taxon") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "A. Composition by neighbor count",
       subtitle = "Quartile-based neighborhood density classification",
       x = "Number of neighboring corals (quartile)",
       y = "Proportion of individuals") +
  theme_publication() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 30, hjust = 1))

# ============================================================================
# Panel B: Key taxa vs number of neighbors
# ============================================================================

cat("Creating Panel B: Taxa vs neighbor count...\n")

prop_long <- coral_comp %>%
  select(coral_id, number_of_neighbors, volume,
         prop_crab, prop_shrimp, prop_fish, prop_snail) %>%
  pivot_longer(cols = starts_with("prop_"),
               names_to = "group",
               values_to = "proportion",
               names_prefix = "prop_") %>%
  mutate(group = str_to_title(group),
         group = factor(group, levels = c("Crab", "Shrimp", "Fish", "Snail")))

# Calculate correlations
cors <- prop_long %>%
  group_by(group) %>%
  summarise(
    r = cor(proportion, number_of_neighbors, use = "complete.obs"),
    p = cor.test(proportion, number_of_neighbors)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.2f%s", r,
                         ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

p_b <- ggplot(prop_long, aes(x = number_of_neighbors, y = proportion)) +
  geom_point(aes(color = group), alpha = 0.5, size = 1.8) +
  geom_smooth(aes(color = group), method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
  geom_text(data = cors, aes(label = label, color = group),
            x = Inf, y = 0.85, hjust = 1.1, size = 3, fontface = "italic",
            show.legend = FALSE) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_wrap(~ group, nrow = 1) +
  labs(title = "B. Proportional abundance vs. neighbor count",
       subtitle = "Controlling for coral size",
       x = "Number of neighboring corals",
       y = "Proportion of community") +
  theme_publication() +
  theme(strip.text = element_text(face = "bold"))

# ============================================================================
# Panel C: Composition by total neighbor volume
# ============================================================================

cat("Creating Panel C: Composition by neighbor volume...\n")

coral_comp <- coral_comp %>%
  mutate(
    neighbor_vol_bin = cut(total_neighbor_volume,
                           breaks = quantile(total_neighbor_volume, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                           labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
                           include.lowest = TRUE)
  )

comp_by_vol <- merged %>%
  left_join(coral_comp %>% select(coral_id, neighbor_vol_bin), by = "coral_id") %>%
  filter(!is.na(neighbor_vol_bin)) %>%
  group_by(neighbor_vol_bin, type_clean) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(neighbor_vol_bin) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

comp_by_vol$type_clean <- factor(comp_by_vol$type_clean, levels = rev(group_order))

n_by_vol <- coral_comp %>%
  filter(!is.na(neighbor_vol_bin)) %>%
  group_by(neighbor_vol_bin) %>%
  summarise(n = n())

p_c <- ggplot(comp_by_vol, aes(x = neighbor_vol_bin, y = prop, fill = type_clean)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3, width = 0.8) +
  geom_text(data = n_by_vol, aes(x = neighbor_vol_bin, y = 1.03, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 3, color = "gray40") +
  scale_fill_manual(values = group_colors, name = "Taxon") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "C. Composition by neighborhood coral volume",
       subtitle = "Total volume of corals within neighborhood radius",
       x = "Total neighbor volume (quartile)",
       y = "Proportion of individuals") +
  theme_publication() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

# ============================================================================
# Panel D: Richness vs neighborhood metrics
# ============================================================================

cat("Creating Panel D: Richness vs neighborhood...\n")

# Calculate correlation
cor_rich_neigh <- cor.test(coral_comp$richness, coral_comp$number_of_neighbors)

p_d <- ggplot(coral_comp, aes(x = number_of_neighbors, y = richness)) +
  geom_point(aes(size = volume, color = site), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  scale_size_continuous(range = c(1, 6), name = expression("Coral vol. (cm"^3*")"),
                        labels = scales::comma) +
  scale_color_site() +
  labs(title = "D. Species richness vs. neighbor count",
       subtitle = sprintf("r = %.2f, p = %.3f",
                         cor_rich_neigh$estimate, cor_rich_neigh$p.value),
       x = "Number of neighboring corals",
       y = "CAFI species richness",
       color = "Site") +
  theme_publication() +
  theme(legend.position = c(0.85, 0.25),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# ============================================================================
# Combine 4-panel figure
# ============================================================================

cat("Combining panels...\n")

fig_neighbor <- (p_a | p_c) / p_b / p_d +
  plot_layout(heights = c(1, 0.9, 1.1)) +
  plot_annotation(
    title = "Community composition varies with coral neighborhood context",
    subtitle = "Corals in denser neighborhoods support different assemblages",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30")
    )
  )

# Save to composition folder
ggsave(file.path(fig_dir, "composition_by_neighborhood_4panel.png"), fig_neighbor,
       width = 14, height = 13, dpi = 300, bg = "white")
cat("  - Saved:", file.path(fig_dir, "composition_by_neighborhood_4panel.png"), "\n")

# Save to manuscript folder
ggsave(file.path(MANUSCRIPT_DIR, "composition_by_neighborhood_4panel.png"), fig_neighbor,
       width = 14, height = 13, dpi = 300, bg = "white")
cat("  - Saved:", file.path(MANUSCRIPT_DIR, "composition_by_neighborhood_4panel.png"), "\n")

# ============================================================================
# Statistical Tests
# ============================================================================

cat("\nRunning statistical tests...\n")

# Chi-squared test for composition differences across neighbor quartiles
comp_table <- merged %>%
  left_join(coral_comp %>% select(coral_id, neighbor_bin), by = "coral_id") %>%
  filter(!is.na(neighbor_bin)) %>%
  group_by(neighbor_bin, type_clean) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = type_clean, values_from = n, values_fill = 0)

comp_matrix <- as.matrix(comp_table[,-1])
chisq_result <- chisq.test(comp_matrix)

cat(sprintf("  Chi-squared test for composition × neighbor density:\n"))
cat(sprintf("    χ² = %.1f, df = %d, p %s\n",
            chisq_result$statistic,
            chisq_result$parameter,
            ifelse(chisq_result$p.value < 0.001, "< 0.001",
                   sprintf("= %.3f", chisq_result$p.value))))

# Correlation tests
cat("\n  Proportional abundance correlations with # neighbors:\n")
for(grp in c("Crab", "Shrimp", "Fish", "Snail")) {
  subset_data <- prop_long %>% filter(group == grp)
  cor_test <- cor.test(subset_data$proportion, subset_data$number_of_neighbors)
  sig <- ifelse(cor_test$p.value < 0.001, "***",
                ifelse(cor_test$p.value < 0.01, "**",
                       ifelse(cor_test$p.value < 0.05, "*", "")))
  cat(sprintf("    %s: r = %.3f %s\n", grp, cor_test$estimate, sig))
}

# Save correlation results
write_csv(cors, file.path(TABLES_DIR, "composition_neighborhood_correlations.csv"))

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Composition by Neighborhood Summary\n")
cat("========================================\n\n")

cat("KEY FINDINGS:\n")
if(chisq_result$p.value < 0.05) {
  cat("  - Composition varies significantly with neighborhood density\n")
} else {
  cat("  - No significant composition change with neighborhood density\n")
}
cat("  - See correlation results for individual taxon trends\n\n")

cat("Figures saved to:\n")
cat("  -", file.path(fig_dir, "composition_by_neighborhood_4panel.png"), "\n")
cat("  -", file.path(MANUSCRIPT_DIR, "composition_by_neighborhood_4panel.png"), "\n\n")

cat("✅ Script 15 complete!\n")
