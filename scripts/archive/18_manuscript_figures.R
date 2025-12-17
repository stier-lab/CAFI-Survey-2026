#!/usr/bin/env Rscript
# ============================================================================
# 18_manuscript_figures.R - Compile All Manuscript-Ready Figures
# ============================================================================
#
# PURPOSE: Generate final publication-quality figures for the manuscript,
#          combining key results from all analyses into cohesive multi-panel
#          figures suitable for peer-reviewed publication.
#
# OUTPUTS:
#   - Figure 1: Study overview and methods
#   - Figure 2: Size-abundance scaling (H2 test)
#   - Figure 3: Community composition patterns
#   - Figure 4: Neighborhood context effects
#   - Figure 5: CAFI-coral condition relationships
#   - Supplementary figures
#
# DEPENDENCIES: 00_setup.R, all previous scripts must have run
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2025-12-05
# ============================================================================

cat("\n========================================\n")
cat("18: Manuscript Figure Compilation\n")
cat("========================================\n\n")

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Load processed data objects
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")

# Ensure manuscript directory exists
dir.create(MANUSCRIPT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# FIGURE 1: Study Overview
# ============================================================================

cat("Creating Figure 1: Study Overview...\n")

# Panel A: Site map placeholder (sampling locations)
site_summary <- coral_master %>%
  group_by(site) %>%
  summarise(
    n_corals = n(),
    total_cafi = sum(total_cafi, na.rm = TRUE),
    mean_volume = mean(volume, na.rm = TRUE),
    .groups = "drop"
  )

fig1_a <- ggplot(site_summary, aes(x = site, y = n_corals, fill = site)) +
  geom_col(width = 0.7, alpha = 0.85, color = "black", linewidth = 0.3) +
  geom_text(aes(label = n_corals), vjust = -0.3, fontface = "bold", size = 4) +
  scale_fill_site() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "A. Sampling effort",
       x = NULL, y = "Corals surveyed") +
  theme_publication() +
  theme(legend.position = "none")

# Panel B: Coral volume distribution
fig1_b <- ggplot(coral_master, aes(x = volume, fill = site)) +
  geom_histogram(bins = 25, alpha = 0.7, color = "white", position = "identity") +
  scale_x_log10(labels = scales::comma) +
  scale_fill_site() +
  labs(title = "B. Coral size distribution",
       x = expression(Volume~(cm^3)), y = "Count") +
  theme_publication() +
  theme(legend.position = c(0.85, 0.75),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# Panel C: CAFI abundance distribution
fig1_c <- ggplot(coral_master, aes(x = total_cafi)) +
  geom_histogram(bins = 25, fill = "#D55E00", alpha = 0.8, color = "white") +
  scale_x_sqrt(breaks = c(0, 25, 100, 225, 400)) +
  labs(title = "C. CAFI abundance distribution",
       x = "Total CAFI (sqrt scale)", y = "Count") +
  theme_publication()

# Panel D: Taxonomic composition
tax_comp <- cafi_clean %>%
  count(type) %>%
  mutate(prop = n / sum(n), type = str_to_title(type)) %>%
  filter(prop > 0.01) %>%
  arrange(desc(prop)) %>%
  mutate(type = fct_reorder(type, prop))

fig1_d <- ggplot(tax_comp, aes(x = type, y = prop, fill = type)) +
  geom_col(width = 0.7, alpha = 0.85, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", prop * 100)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_taxon() +
  labs(title = "D. Community composition",
       x = NULL, y = "Proportion") +
  theme_publication() +
  theme(legend.position = "none")

# Combine Figure 1
fig1 <- (fig1_a | fig1_b) / (fig1_c | fig1_d) +
  plot_annotation(
    title = "Figure 1. Study overview: CAFI communities on Pocillopora spp. at Mo'orea",
    subtitle = sprintf("n = %d corals | %d CAFI individuals | %d taxa",
                      nrow(coral_master),
                      sum(coral_master$total_cafi),
                      n_distinct(cafi_clean$species)),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray30")
    )
  )

ggsave(file.path(MANUSCRIPT_DIR, "Figure1_study_overview.png"), fig1,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("  - Saved: Figure1_study_overview.png\n")

# ============================================================================
# FIGURE 2: Size-Abundance Scaling
# ============================================================================

cat("Creating Figure 2: Size-Abundance Scaling...\n")

scaling_data <- coral_master %>%
  filter(!is.na(volume), volume > 0) %>%
  mutate(
    log_vol = log(volume),
    log_abund = log(total_cafi + 1),
    cafi_density = total_cafi / volume
  )

# Fit model for statistics
m_scaling <- lm(log_abund ~ log_vol, data = scaling_data)
beta <- coef(m_scaling)[2]
beta_ci <- confint(m_scaling)[2,]
r2 <- summary(m_scaling)$r.squared

# Panel A: Main scaling plot
fig2_a <- ggplot(scaling_data, aes(x = volume, y = total_cafi)) +
  geom_point(aes(color = site), alpha = 0.65, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1.2, se = TRUE) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  scale_color_site() +
  annotation_logticks() +
  labs(title = "A. Power-law scaling",
       subtitle = sprintf("β = %.2f [%.2f, %.2f] | R² = %.2f",
                         beta, beta_ci[1], beta_ci[2], r2),
       x = expression(Coral~Volume~(cm^3)),
       y = "CAFI Abundance",
       color = "Site") +
  theme_publication() +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# Panel B: Density scaling
fig2_b <- scaling_data %>%
  filter(cafi_density > 0) %>%
  ggplot(aes(x = volume, y = cafi_density)) +
  geom_point(aes(color = site), alpha = 0.65, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1.2, se = TRUE) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  scale_color_site() +
  annotation_logticks() +
  labs(title = "B. Density scaling (propagule dilution)",
       subtitle = "Larger corals have lower CAFI per unit volume",
       x = expression(Coral~Volume~(cm^3)),
       y = expression(CAFI~Density~(ind./cm^3))) +
  theme_publication() +
  theme(legend.position = "none")

# Panel C: Richness scaling
fig2_c <- ggplot(scaling_data, aes(x = volume, y = otu_richness)) +
  geom_point(aes(color = site), alpha = 0.65, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1.2, se = TRUE) +
  scale_x_log10(labels = scales::comma) +
  scale_color_site() +
  labs(title = "C. Species-area relationship",
       x = expression(Coral~Volume~(cm^3)),
       y = "OTU Richness") +
  theme_publication() +
  theme(legend.position = "none")

# Panel D: Site comparison
site_slopes <- scaling_data %>%
  group_by(site) %>%
  summarise(
    beta = coef(lm(log_abund ~ log_vol))[2],
    se = summary(lm(log_abund ~ log_vol))$coefficients[2, 2],
    .groups = "drop"
  ) %>%
  mutate(ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se)

fig2_d <- ggplot(site_slopes, aes(x = site, y = beta, color = site)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = beta, linetype = "solid", alpha = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 1) +
  geom_point(size = 4) +
  scale_color_site() +
  labs(title = "D. Site-specific scaling",
       subtitle = "Red = isometric; Black = overall",
       x = NULL, y = "Scaling exponent (β)") +
  theme_publication() +
  theme(legend.position = "none")

# Combine Figure 2
fig2 <- (fig2_a | fig2_b) / (fig2_c | fig2_d) +
  plot_annotation(
    title = "Figure 2. CAFI abundance scales sublinearly with coral size",
    subtitle = "Hypothesis H2: Propagule redirection predicts β < 1",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray30")
    )
  )

ggsave(file.path(MANUSCRIPT_DIR, "Figure2_size_scaling.png"), fig2,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("  - Saved: Figure2_size_scaling.png\n")

# ============================================================================
# FIGURE 3: Community Composition
# ============================================================================

cat("Creating Figure 3: Community Composition...\n")

# Prepare data with type labels
# cafi_clean already has site, only need volume from coral_master
cafi_typed <- cafi_clean %>%
  mutate(type_label = str_to_title(type)) %>%
  left_join(coral_master %>% select(coral_id, volume), by = "coral_id") %>%
  filter(!is.na(volume))

# Size quintiles
cafi_typed <- cafi_typed %>%
  mutate(size_quintile = ntile(volume, 5))

# Composition by size
comp_by_size <- cafi_typed %>%
  group_by(size_quintile, type_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(size_quintile) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

group_order <- comp_by_size %>%
  group_by(type_label) %>%
  summarise(total = sum(n)) %>%
  arrange(desc(total)) %>%
  pull(type_label)

comp_by_size$type_label <- factor(comp_by_size$type_label, levels = rev(group_order))

# Group colors
group_cols <- c("Crab" = "#D55E00", "Shrimp" = "#0072B2", "Fish" = "#009E73",
                "Snail" = "#CC79A7", "Hermit" = "#E69F00", "Other" = "#999999")

fig3_a <- ggplot(comp_by_size, aes(x = factor(size_quintile), y = prop, fill = type_label)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3, width = 0.8) +
  scale_fill_manual(values = group_cols, name = "Taxon") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "A. Composition by coral size quintile",
       x = "Size quintile (1 = smallest)", y = "Proportion") +
  theme_publication()

# By site
comp_by_site <- cafi_typed %>%
  group_by(site, type_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(site) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

comp_by_site$type_label <- factor(comp_by_site$type_label, levels = rev(group_order))

fig3_b <- ggplot(comp_by_site, aes(x = site, y = prop, fill = type_label)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = group_cols, name = "Taxon") +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "B. Composition by site",
       x = NULL, y = "Proportion") +
  theme_publication() +
  theme(legend.position = "none")

# Combine Figure 3
fig3 <- (fig3_a | fig3_b) +
  plot_layout(widths = c(1.5, 1)) +
  plot_annotation(
    title = "Figure 3. Community composition varies with coral size and site",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(MANUSCRIPT_DIR, "Figure3_community_composition.png"), fig3,
       width = 14, height = 6, dpi = 300, bg = "white")
cat("  - Saved: Figure3_community_composition.png\n")

# ============================================================================
# FIGURE 4: CAFI-Condition Relationships
# ============================================================================

cat("Creating Figure 4: CAFI-Condition Relationships...\n")

# coral_master already has condition_score - filter for corals with data
cafi_condition <- coral_master %>%
  filter(!is.na(condition_score))

if (nrow(cafi_condition) > 10) {
  fig4_a <- ggplot(cafi_condition, aes(x = total_cafi, y = condition_score)) +
    geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", color = "black", linewidth = 1.2, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_x_sqrt() +
    scale_color_site() +
    labs(title = "A. CAFI abundance vs. coral condition",
         x = "CAFI abundance (sqrt)", y = "Condition score") +
    theme_publication() +
    theme(legend.position = c(0.85, 0.15))

  fig4_b <- ggplot(cafi_condition, aes(x = otu_richness, y = condition_score)) +
    geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
    geom_smooth(method = "lm", color = "black", linewidth = 1.2, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_site() +
    labs(title = "B. OTU richness vs. coral condition",
         x = "OTU richness", y = "Condition score") +
    theme_publication() +
    theme(legend.position = "none")

  fig4 <- (fig4_a | fig4_b) +
    plot_annotation(
      title = "Figure 4. CAFI-coral condition relationships",
      subtitle = "Condition = position-corrected PC1 (protein, carbs, zooxanthellae, AFDW)",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "gray30")
      )
    )

  ggsave(file.path(MANUSCRIPT_DIR, "Figure4_condition_relationships.png"), fig4,
         width = 12, height = 6, dpi = 300, bg = "white")
  cat("  - Saved: Figure4_condition_relationships.png\n")
}

# ============================================================================
# List All Manuscript Figures
# ============================================================================

cat("\n========================================\n")
cat("Manuscript Figures Summary\n")
cat("========================================\n\n")

ms_figs <- list.files(MANUSCRIPT_DIR, pattern = "\\.png$", full.names = FALSE)
cat(sprintf("Generated %d manuscript figures:\n", length(ms_figs)))
for (f in sort(ms_figs)) {
  cat(sprintf("  - %s\n", f))
}

cat("\nAll figures saved to:", MANUSCRIPT_DIR, "\n\n")

cat("✅ Script 18 complete!\n")
