#!/usr/bin/env Rscript
# ============================================================================
# 16_summary_figures.R - Generate Summary Figures Across All Analyses
# ============================================================================
#
# PURPOSE: Create comprehensive summary figures combining key results from
#          all previous analysis scripts into publication-ready visualizations
#
# OUTPUTS:
#   - Multi-panel summary figures
#   - Key results table
#   - Manuscript-ready composite figures
#
# DEPENDENCIES: 00_setup.R, all previous scripts must have run
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("16: Summary Figures Generation\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")

# Create figure subdirectory
fig_dir <- file.path(FIGURES_DIR, "summary")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Panel 1: Study Overview
# ============================================================================

cat("Creating study overview panel...\n")

# Site summary
site_summary <- coral_master %>%
  group_by(site) %>%
  summarise(
    n_corals = n(),
    total_cafi = sum(total_cafi, na.rm = TRUE),
    mean_richness = mean(otu_richness, na.rm = TRUE),
    mean_volume = mean(volume, na.rm = TRUE),
    .groups = "drop"
  )

p_sites <- ggplot(site_summary, aes(x = site, y = n_corals, fill = site)) +
  geom_col(width = 0.7, alpha = 0.8) +
  geom_text(aes(label = n_corals), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_site() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "A. Corals surveyed by site",
       x = NULL, y = "Number of corals") +
  theme_publication() +
  theme(legend.position = "none")

# Volume distribution
p_volume <- ggplot(coral_master, aes(x = volume)) +
  geom_histogram(bins = 30, fill = "#0072B2", color = "white", alpha = 0.8) +
  scale_x_log10(labels = scales::comma) +
  labs(title = "B. Coral size distribution",
       subtitle = sprintf("Median = %.0f cm³, Range = %.0f–%.0f cm³",
                         median(coral_master$volume, na.rm = TRUE),
                         min(coral_master$volume, na.rm = TRUE),
                         max(coral_master$volume, na.rm = TRUE)),
       x = expression(Coral~Volume~(cm^3~","~log~scale)),
       y = "Count") +
  theme_publication()

# CAFI abundance distribution
p_cafi <- ggplot(coral_master, aes(x = total_cafi)) +
  geom_histogram(bins = 30, fill = "#D55E00", color = "white", alpha = 0.8) +
  scale_x_sqrt(breaks = c(0, 25, 100, 225, 400)) +
  labs(title = "C. CAFI abundance distribution",
       subtitle = sprintf("Median = %.0f, Mean = %.1f ± %.1f",
                         median(coral_master$total_cafi, na.rm = TRUE),
                         mean(coral_master$total_cafi, na.rm = TRUE),
                         sd(coral_master$total_cafi, na.rm = TRUE)),
       x = "CAFI abundance (sqrt scale)",
       y = "Count") +
  theme_publication()

# OTU richness distribution
p_richness <- ggplot(coral_master, aes(x = otu_richness)) +
  geom_histogram(bins = 20, fill = "#009E73", color = "white", alpha = 0.8) +
  labs(title = "D. OTU richness distribution",
       subtitle = sprintf("Median = %.0f, Mean = %.1f ± %.1f",
                         median(coral_master$otu_richness, na.rm = TRUE),
                         mean(coral_master$otu_richness, na.rm = TRUE),
                         sd(coral_master$otu_richness, na.rm = TRUE)),
       x = "OTU richness",
       y = "Count") +
  theme_publication()

# Combine overview
fig_overview <- (p_sites | p_volume) / (p_cafi | p_richness) +
  plot_annotation(
    title = "Study Overview: CAFI on Pocillopora spp. at Mo'orea",
    subtitle = sprintf("n = %d corals surveyed across 3 sites | %d total CAFI individuals | %d OTUs",
                      nrow(coral_master),
                      sum(coral_master$total_cafi, na.rm = TRUE),
                      n_distinct(cafi_clean$species)),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30")
    )
  )

ggsave(file.path(fig_dir, "study_overview.png"), fig_overview,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("  - Saved: study_overview.png\n")

# ============================================================================
# Panel 2: Key Relationships Summary
# ============================================================================

cat("Creating key relationships summary...\n")

# Size-abundance relationship
p_size_abund <- ggplot(coral_master, aes(x = volume, y = total_cafi)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_y_sqrt() +
  scale_color_site() +
  labs(title = "A. Size-abundance scaling",
       x = expression(Volume~(cm^3)), y = "CAFI abundance") +
  theme_publication() +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.8)))

# Size-richness relationship
p_size_rich <- ggplot(coral_master, aes(x = volume, y = otu_richness)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
  scale_x_log10(labels = scales::comma) +
  scale_color_site() +
  labs(title = "B. Size-richness relationship",
       x = expression(Volume~(cm^3)), y = "OTU richness") +
  theme_publication() +
  theme(legend.position = "none")

# Abundance-richness relationship
p_abund_rich <- ggplot(coral_master, aes(x = total_cafi, y = otu_richness)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
  scale_x_sqrt() +
  scale_color_site() +
  labs(title = "C. Abundance-richness relationship",
       x = "CAFI abundance (sqrt)", y = "OTU richness") +
  theme_publication() +
  theme(legend.position = "none")

# Taxonomic composition
tax_comp <- cafi_clean %>%
  count(type) %>%
  mutate(prop = n / sum(n),
         type = str_to_title(type)) %>%
  arrange(desc(prop)) %>%
  mutate(type = fct_reorder(type, prop))

p_tax <- ggplot(tax_comp, aes(x = type, y = prop, fill = type)) +
  geom_col(width = 0.7, alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", prop * 100)), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_taxon() +
  labs(title = "D. Taxonomic composition",
       x = NULL, y = "Proportion") +
  theme_publication() +
  theme(legend.position = "none")

fig_relationships <- (p_size_abund | p_size_rich) / (p_abund_rich | p_tax) +
  plot_annotation(
    title = "Key Ecological Relationships",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(fig_dir, "key_relationships.png"), fig_relationships,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("  - Saved: key_relationships.png\n")

# ============================================================================
# Panel 3: Condition Effects Summary
# ============================================================================

cat("Creating condition effects summary...\n")

# coral_master already has condition_score - filter for corals with condition data
coral_condition <- coral_master %>%
  filter(!is.na(condition_score))

if (nrow(coral_condition) > 10) {

  # CAFI vs condition
  p_cafi_cond <- ggplot(coral_condition, aes(x = total_cafi, y = condition_score)) +
    geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_x_sqrt() +
    scale_color_site() +
    labs(title = "A. CAFI abundance vs. coral condition",
         x = "CAFI abundance (sqrt scale)",
         y = "Condition score (position-corrected)") +
    theme_publication() +
    theme(legend.position = c(0.85, 0.85))

  # Richness vs condition
  p_rich_cond <- ggplot(coral_condition, aes(x = otu_richness, y = condition_score)) +
    geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_site() +
    labs(title = "B. OTU richness vs. coral condition",
         x = "OTU richness",
         y = "Condition score") +
    theme_publication() +
    theme(legend.position = "none")

  fig_condition <- (p_cafi_cond | p_rich_cond) +
    plot_annotation(
      title = "CAFI-Coral Condition Relationships",
      subtitle = "Condition = position-corrected PC1 of protein, carbs, zoox, AFDW",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray30")
      )
    )

  ggsave(file.path(fig_dir, "condition_effects.png"), fig_condition,
         width = 12, height = 6, dpi = 300, bg = "white")
  cat("  - Saved: condition_effects.png\n")
}

# ============================================================================
# Panel 4: Neighborhood Effects Summary
# ============================================================================

cat("Creating neighborhood effects summary...\n")

if ("number_of_neighbors" %in% names(coral_master) &&
    sum(!is.na(coral_master$number_of_neighbors)) > 10) {

  coral_neigh <- coral_master %>%
    filter(!is.na(number_of_neighbors))

  # Abundance vs neighbors
  p_neigh_abund <- ggplot(coral_neigh, aes(x = number_of_neighbors, y = total_cafi)) +
    geom_point(aes(color = site, size = volume), alpha = 0.6) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
    scale_y_sqrt() +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    scale_color_site() +
    labs(title = "A. CAFI abundance vs. neighbor count",
         x = "Number of neighboring corals",
         y = "CAFI abundance (sqrt)") +
    theme_publication() +
    theme(legend.position = c(0.85, 0.85))

  # Richness vs neighbors
  p_neigh_rich <- ggplot(coral_neigh, aes(x = number_of_neighbors, y = otu_richness)) +
    geom_point(aes(color = site, size = volume), alpha = 0.6) +
    geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 1) +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    scale_color_site() +
    labs(title = "B. OTU richness vs. neighbor count",
         x = "Number of neighboring corals",
         y = "OTU richness") +
    theme_publication() +
    theme(legend.position = "none")

  fig_neighborhood <- (p_neigh_abund | p_neigh_rich) +
    plot_annotation(
      title = "Neighborhood Context Effects",
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )

  ggsave(file.path(fig_dir, "neighborhood_effects.png"), fig_neighborhood,
         width = 12, height = 6, dpi = 300, bg = "white")
  cat("  - Saved: neighborhood_effects.png\n")
}

# ============================================================================
# Key Results Table
# ============================================================================

cat("\nCompiling key results table...\n")

key_results <- tibble(
  analysis = c(
    "Sample size (corals)",
    "Sample size (CAFI individuals)",
    "Number of OTUs",
    "Sites surveyed",
    "Mean CAFI per coral",
    "Mean OTU richness per coral",
    "Size-abundance correlation (r)",
    "Size-richness correlation (r)",
    "Abundance-richness correlation (r)"
  ),
  value = c(
    nrow(coral_master),
    sum(coral_master$total_cafi, na.rm = TRUE),
    n_distinct(cafi_clean$species),
    n_distinct(coral_master$site),
    round(mean(coral_master$total_cafi, na.rm = TRUE), 1),
    round(mean(coral_master$otu_richness, na.rm = TRUE), 1),
    round(cor(log(coral_master$volume), coral_master$total_cafi, use = "complete.obs"), 3),
    round(cor(log(coral_master$volume), coral_master$otu_richness, use = "complete.obs"), 3),
    round(cor(coral_master$total_cafi, coral_master$otu_richness, use = "complete.obs"), 3)
  )
)

write_csv(key_results, file.path(TABLES_DIR, "key_results_summary.csv"))

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Summary Figures Complete\n")
cat("========================================\n\n")

cat("Figures saved to:", fig_dir, "\n")
cat("  - study_overview.png\n")
cat("  - key_relationships.png\n")
cat("  - condition_effects.png\n")
cat("  - neighborhood_effects.png\n\n")

cat("Key results table saved to:", file.path(TABLES_DIR, "key_results_summary.csv"), "\n\n")

cat("✅ Script 16 complete!\n")
