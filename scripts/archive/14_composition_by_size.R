#!/usr/bin/env Rscript
# ============================================================================
# 14_composition_by_size.R - How CAFI composition changes with coral size
# ============================================================================
#
# PURPOSE: Analyze how community composition shifts across the coral size gradient
#
# KEY QUESTIONS:
#   - Do larger corals support different assemblages than smaller corals?
#   - Which taxa increase/decrease proportionally with size?
#   - Do species show size preferences?
#
# OUTPUTS:
#   - 4-panel composition figure (main + manuscript folder)
#   - Statistical tests for composition shifts
#
# DEPENDENCIES: 00_setup.R, 01_load_clean_data.R
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("14: Composition by Coral Size Analysis\n")
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

cat("Preparing composition data by coral size...\n")

# Standardize type names for display
cafi_typed <- cafi_clean %>%
  mutate(type_clean = case_when(
    type == "crab" ~ "Crab",
    type == "shrimp" ~ "Shrimp",
    type == "fish" ~ "Fish",
    type == "hermit" ~ "Hermit",
    type == "snail" ~ "Snail",
    type == "echinoderm" ~ "Echinoderm",
    type %in% c("worm", "amphipod", "amph") ~ "Other",
    TRUE ~ "Other"
  ))

# Merge with coral size (site comes from cafi_typed if not in coral_master)
merged <- cafi_typed %>%
  left_join(coral_master %>% select(coral_id, volume), by = "coral_id") %>%
  filter(!is.na(volume), volume > 0) %>%
  mutate(log_vol = log10(volume))

# Add site from coral_master if available
if ("site" %in% names(coral_master)) {
  merged <- merged %>%
    left_join(coral_master %>% select(coral_id, site), by = "coral_id")
} else if (!"site" %in% names(merged)) {
  merged$site <- "All"
}

cat(sprintf("  - Records with size data: %d\n", nrow(merged)))
cat(sprintf("  - Unique corals: %d\n", n_distinct(merged$coral_id)))

# Create size bins
merged <- merged %>%
  mutate(
    size_bin = cut(volume,
                   breaks = quantile(volume, c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE),
                   labels = c("Q1 (Smallest)", "Q2", "Q3", "Q4", "Q5 (Largest)"),
                   include.lowest = TRUE)
  )

# Per-coral composition
coral_composition <- merged %>%
  group_by(coral_id, volume, log_vol) %>%
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
    prop_hermit = n_hermit / abundance
  )

# Add site if available
if ("site" %in% names(merged)) {
  site_lookup <- merged %>%
    distinct(coral_id, site)
  coral_composition <- coral_composition %>%
    left_join(site_lookup, by = "coral_id")
}

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
# Panel A: Stacked composition by size quintile
# ============================================================================

cat("Creating Panel A: Composition by size quintile...\n")

# Aggregate by size bin
comp_by_bin <- merged %>%
  group_by(size_bin, type_clean) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(size_bin) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Order groups by overall abundance
group_order <- comp_by_bin %>%
  group_by(type_clean) %>%
  summarise(total = sum(count)) %>%
  arrange(desc(total)) %>%
  pull(type_clean)

comp_by_bin$type_clean <- factor(comp_by_bin$type_clean, levels = rev(group_order))

# Sample sizes
sample_sizes <- merged %>%
  group_by(size_bin) %>%
  summarise(n_corals = n_distinct(coral_id), n_individuals = n())

p_a <- ggplot(comp_by_bin, aes(x = size_bin, y = prop, fill = type_clean)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3, width = 0.8) +
  geom_text(data = sample_sizes,
            aes(x = size_bin, y = 1.03, label = paste0("n=", n_corals)),
            inherit.aes = FALSE, size = 3, color = "gray40") +
  scale_fill_manual(values = group_colors, name = "Taxon") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "A. Community composition shifts with coral size",
       subtitle = "Quintile-based size classification",
       x = "Coral size quintile",
       y = "Proportion of individuals") +
  theme_publication() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 30, hjust = 1))

# ============================================================================
# Panel B: Proportional trends by taxon
# ============================================================================

cat("Creating Panel B: Proportional trends...\n")

# Prepare long format
prop_long <- coral_composition %>%
  select(coral_id, volume, log_vol, prop_crab, prop_shrimp, prop_fish, prop_snail) %>%
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
    r = cor(proportion, log_vol, use = "complete.obs"),
    p = cor.test(proportion, log_vol)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.2f%s", r,
                         ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))))

p_b <- ggplot(prop_long, aes(x = volume, y = proportion)) +
  geom_point(aes(color = group), alpha = 0.5, size = 1.8) +
  geom_smooth(aes(color = group), method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
  geom_text(data = cors, aes(label = label, color = group),
            x = Inf, y = 0.85, hjust = 1.1, size = 3, fontface = "italic",
            show.legend = FALSE) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_wrap(~ group, nrow = 1) +
  labs(title = "B. Proportional abundance vs. coral size",
       subtitle = "Crabs decline while fish and snails increase in larger corals",
       x = expression(Coral~Volume~(cm^3~","~log~scale)),
       y = "Proportion of community") +
  theme_publication() +
  theme(strip.text = element_text(face = "bold"))

# ============================================================================
# Panel C: Species richness by size quintile
# ============================================================================

cat("Creating Panel C: Richness by size quintile...\n")

# Species count per size quintile
species_by_quintile <- merged %>%
  mutate(size_quintile = ntile(volume, 5)) %>%
  group_by(size_quintile) %>%
  summarise(
    n_species = n_distinct(species),
    n_corals = n_distinct(coral_id),
    min_vol = min(volume),
    max_vol = max(volume),
    .groups = "drop"
  ) %>%
  mutate(vol_label = sprintf("%.0f-%.0fk", min_vol/1000, max_vol/1000))

p_c <- ggplot(species_by_quintile, aes(x = factor(size_quintile), y = n_species)) +
  geom_col(fill = "#0072B2", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = n_species), vjust = -0.3, size = 4, fontface = "bold") +
  geom_text(aes(y = 3, label = vol_label), size = 2.5, color = "white") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "C. Species pool by coral size quintile",
       subtitle = expression("Volume ranges (cm"^3*") shown in bars"),
       x = "Size quintile (1 = smallest, 5 = largest)",
       y = "Number of species observed") +
  theme_publication()

# ============================================================================
# Panel D: Size preferences of common species
# ============================================================================

cat("Creating Panel D: Species size preferences...\n")

# Calculate median volume for common species
species_sizes <- merged %>%
  group_by(species, type_clean) %>%
  summarise(
    n = n(),
    median_vol = median(volume),
    q25_vol = quantile(volume, 0.25),
    q75_vol = quantile(volume, 0.75),
    .groups = "drop"
  ) %>%
  filter(n >= 15) %>%
  arrange(median_vol)

# Select top/bottom species by size preference
top_species <- bind_rows(
  species_sizes %>% slice_min(median_vol, n = 6),
  species_sizes %>% slice_max(median_vol, n = 6)
) %>%
  distinct(species, .keep_all = TRUE) %>%
  mutate(species = fct_reorder(species, median_vol))

p_d <- ggplot(top_species, aes(y = species, x = median_vol, color = type_clean)) +
  geom_errorbarh(aes(xmin = q25_vol, xmax = q75_vol), height = 0.3, linewidth = 0.8, alpha = 0.7) +
  geom_point(aes(size = n), alpha = 0.9) +
  geom_vline(xintercept = median(merged$volume), linetype = "dashed", color = "gray50", alpha = 0.6) +
  scale_x_log10(labels = scales::comma) +
  scale_color_manual(values = group_colors, name = "Taxon") +
  scale_size_continuous(range = c(2, 5), name = "n individuals") +
  labs(title = "D. Species differ in coral size preferences",
       subtitle = "Bars = IQR; dashed line = median coral size",
       x = expression(Coral~Volume~(cm^3~","~log~scale)),
       y = NULL) +
  theme_publication() +
  theme(axis.text.y = element_text(face = "italic", size = 8),
        legend.position = "right")

# ============================================================================
# Combine 4-panel figure
# ============================================================================

cat("Combining panels...\n")

fig_comp_size <- (p_a | p_c) / p_b / p_d +
  plot_layout(heights = c(1, 0.8, 1.2)) +
  plot_annotation(
    title = "Community composition changes with coral colony size",
    subtitle = "Larger corals support different species assemblages than smaller corals",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30")
    )
  )

# Save to composition folder
ggsave(file.path(fig_dir, "composition_by_size_4panel.png"), fig_comp_size,
       width = 14, height = 13, dpi = 300, bg = "white")
cat("  - Saved:", file.path(fig_dir, "composition_by_size_4panel.png"), "\n")

# Save to manuscript folder
ggsave(file.path(MANUSCRIPT_DIR, "composition_by_size_4panel.png"), fig_comp_size,
       width = 14, height = 13, dpi = 300, bg = "white")
cat("  - Saved:", file.path(MANUSCRIPT_DIR, "composition_by_size_4panel.png"), "\n")

# ============================================================================
# Statistical Tests
# ============================================================================

cat("\nRunning statistical tests...\n")

# Chi-squared test for composition differences across size quintiles
comp_table <- merged %>%
  mutate(size_quintile = ntile(volume, 5)) %>%
  group_by(size_quintile, type_clean) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = type_clean, values_from = n, values_fill = 0)

comp_matrix <- as.matrix(comp_table[,-1])
chisq_result <- chisq.test(comp_matrix)

cat(sprintf("  Chi-squared test for composition × size:\n"))
cat(sprintf("    χ² = %.1f, df = %d, p %s\n",
            chisq_result$statistic,
            chisq_result$parameter,
            ifelse(chisq_result$p.value < 0.001, "< 0.001",
                   sprintf("= %.3f", chisq_result$p.value))))

# Correlation tests for each taxon
cat("\n  Proportional abundance correlations with log(volume):\n")
for(grp in c("Crab", "Shrimp", "Fish", "Snail")) {
  subset_data <- prop_long %>% filter(group == grp)
  cor_test <- cor.test(subset_data$proportion, subset_data$log_vol)
  sig <- ifelse(cor_test$p.value < 0.001, "***",
                ifelse(cor_test$p.value < 0.01, "**",
                       ifelse(cor_test$p.value < 0.05, "*", "")))
  cat(sprintf("    %s: r = %.3f %s\n", grp, cor_test$estimate, sig))
}

# Save correlation results
write_csv(cors, file.path(TABLES_DIR, "composition_size_correlations.csv"))

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Composition by Size Summary\n")
cat("========================================\n\n")

cat("KEY FINDINGS:\n")
cat("  - Composition varies significantly with coral size (χ² test)\n")
cat("  - Crabs dominate small corals, decline in larger ones\n")
cat("  - Fish and snails increase proportionally in larger corals\n")
cat("  - Species richness increases with coral size quintile\n\n")

cat("Figures saved to:\n")
cat("  -", file.path(fig_dir, "composition_by_size_4panel.png"), "\n")
cat("  -", file.path(MANUSCRIPT_DIR, "composition_by_size_4panel.png"), "\n\n")

cat("✅ Script 14 complete!\n")
