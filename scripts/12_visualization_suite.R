#!/usr/bin/env Rscript
# ============================================================================
# 12_visualization_suite.R - Comprehensive visualization suite for Survey data
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
# ============================================================================

cat("\n========================================\n")
cat("Comprehensive Visualization Suite\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/Survey/00_load_libraries.R"))
library(patchwork)
library(ggridges)
library(ggrepel)
library(viridis)
library(RColorBrewer)

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "visualization_suite")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Set consistent theme
theme_survey <- theme_minimal() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "gray50"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

theme_set(theme_survey)

# ============================================================================
# Color Palettes
# ============================================================================

# Define consistent color palettes
morphotype_colors <- c("verrucosa" = "#E41A1C", "meandrina" = "#377EB8",
                       "other" = "#4DAF4A")
type_colors <- c("crab" = "#FF7F00", "shrimp" = "#6A3D9A",
                "fish" = "#1F78B4", "snail" = "#B2DF8A")
site_colors <- viridis(length(unique(metadata$site)))
names(site_colors) <- unique(metadata$site)

# ============================================================================
# 1. Overview Dashboard
# ============================================================================

cat("Creating overview dashboard...\n")

# Prepare summary data
summary_data <- metadata %>%
  left_join(
    cafi_clean %>%
      group_by(coral_id) %>%
      summarise(
        total_cafi = n(),
        species_richness = n_distinct(species),
        shannon = vegan::diversity(table(species)),
        .groups = "drop"
      ),
    by = "coral_id"
  ) %>%
  mutate(across(c(total_cafi:shannon), ~replace_na(., 0)))

# A. Distribution plots
p1_abundance <- ggplot(summary_data, aes(x = total_cafi)) +
  geom_histogram(fill = "steelblue", alpha = 0.7, bins = 30) +
  geom_vline(xintercept = median(summary_data$total_cafi),
             linetype = "dashed", color = "red") +
  labs(title = "A. CAFI Abundance Distribution",
       x = "Total CAFI per Coral",
       y = "Number of Corals")

p1_richness <- ggplot(summary_data, aes(x = species_richness)) +
  geom_histogram(fill = "darkgreen", alpha = 0.7, bins = 20) +
  geom_vline(xintercept = median(summary_data$species_richness),
             linetype = "dashed", color = "red") +
  labs(title = "B. CAFI OTU Richness Distribution",
       x = "Number of OTUs per Coral",
       y = "Number of Corals")

p1_diversity <- ggplot(summary_data, aes(x = shannon)) +
  geom_histogram(fill = "purple", alpha = 0.7, bins = 30) +
  geom_vline(xintercept = median(summary_data$shannon),
             linetype = "dashed", color = "red") +
  labs(title = "C. Shannon Diversity Distribution",
       x = "Shannon Diversity Index",
       y = "Number of Corals")

# B. Morphotype comparison
p1_morph <- summary_data %>%
  pivot_longer(cols = c(total_cafi, species_richness, shannon),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric,
                        levels = c("total_cafi", "species_richness", "shannon"),
                        labels = c("Abundance", "Richness", "Shannon"))) %>%
  ggplot(aes(x = morphotype, y = value, fill = morphotype)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = morphotype_colors) +
  labs(title = "D. Community Metrics by Morphotype",
       x = "Morphotype",
       y = "Value") +
  theme(legend.position = "none")

# Combine dashboard
dashboard <- (p1_abundance | p1_richness | p1_diversity) /
             p1_morph +
  plot_annotation(
    title = "Survey Data Overview Dashboard",
    subtitle = paste("N =", nrow(summary_data), "corals,",
                    n_distinct(cafi_clean$species), "OTUs (morphotype-based groupings),",
                    sum(summary_data$total_cafi), "total CAFI")
  )

ggsave(file.path(fig_dir, "overview_dashboard.png"),
       dashboard, width = 16, height = 12, dpi = 300)

cat("✓ Overview dashboard created\n\n")

# ============================================================================
# 2. Community Composition Visualization
# ============================================================================

cat("Creating community composition visualizations...\n")

# Top OTUs by abundance
top_species <- cafi_clean %>%
  count(species) %>%
  top_n(20, n) %>%
  pull(species)

# A. OTU rank abundance
p2_rank <- cafi_clean %>%
  count(species) %>%
  arrange(desc(n)) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x = rank, y = n)) +
  geom_point(alpha = 0.6) +
  geom_line(alpha = 0.3) +
  scale_y_log10() +
  scale_x_log10() +
  labs(title = "CAFI OTU Rank-Abundance Curve",
       subtitle = "Morphotype-based groupings (not genetic species)",
       x = "OTU Rank (log scale)",
       y = "Abundance (log scale)")

# B. Top OTU heatmap by morphotype
# Note: These are abundance patterns, not ecological specialization
species_morph_matrix <- cafi_clean %>%
  filter(species %in% top_species) %>%
  left_join(metadata %>% select(coral_id, morphotype), by = "coral_id") %>%
  count(species, morphotype) %>%
  pivot_wider(names_from = morphotype, values_from = n, values_fill = 0) %>%
  column_to_rownames("species") %>%
  as.matrix()

# Scale by row (OTU)
species_morph_scaled <- t(scale(t(species_morph_matrix)))

# Create heatmap data
heatmap_data <- expand.grid(
  species = rownames(species_morph_scaled),
  morphotype = colnames(species_morph_scaled)
) %>%
  mutate(
    value = as.vector(species_morph_scaled),
    species = factor(species, levels = rev(rownames(species_morph_scaled)))
  )

p2_heatmap <- ggplot(heatmap_data, aes(x = morphotype, y = species, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                      midpoint = 0, name = "Z-score") +
  labs(title = "CAFI OTU Distribution Across Morphotypes",
       subtitle = "Abundance patterns (community-level)",
       x = "Morphotype",
       y = "CAFI OTU") +
  theme(axis.text.y = element_text(size = 8))

# C. Community composition stacked bar
p2_stacked <- cafi_clean %>%
  filter(species %in% top_species) %>%
  left_join(metadata %>% select(coral_id, site), by = "coral_id") %>%
  count(site, species) %>%
  group_by(site) %>%
  mutate(proportion = n / sum(n)) %>%
  ggplot(aes(x = site, y = proportion, fill = species)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d(guide = guide_legend(ncol = 2)) +
  labs(title = "Community Composition by Site",
       x = "Site",
       y = "Proportion",
       fill = "CAFI OTU") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 8))

# Combine community plots
species_viz <- p2_rank / (p2_heatmap | p2_stacked) +
  plot_annotation(title = "Community Composition Analysis")

ggsave(file.path(fig_dir, "community_composition_analysis.png"),
       species_viz, width = 16, height = 14, dpi = 300)

cat("✓ Community composition visualizations created\n\n")

# ============================================================================
# 3. Morphotype-Specific Patterns
# ============================================================================

cat("Creating morphotype-specific visualizations...\n")

# A. Ridge plots by morphotype
p3_ridges <- cafi_clean %>%
  filter(!is.na(size_mm), size_mm > 0) %>%
  left_join(metadata %>% select(coral_id, morphotype), by = "coral_id") %>%
  filter(type %in% c("crab", "shrimp")) %>%
  ggplot(aes(x = size_mm, y = morphotype, fill = morphotype)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +
  facet_wrap(~type, scales = "free_x") +
  scale_fill_manual(values = morphotype_colors) +
  labs(title = "Size Distributions by Morphotype and Type",
       x = "Size (mm)",
       y = "Morphotype") +
  theme(legend.position = "none")

# B. Morphotype × Type interaction
interaction_data <- cafi_clean %>%
  left_join(metadata %>% select(coral_id, morphotype), by = "coral_id") %>%
  count(morphotype, type) %>%
  group_by(morphotype) %>%
  mutate(proportion = n / sum(n))

p3_interaction <- ggplot(interaction_data,
                         aes(x = morphotype, y = proportion, fill = type)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = type_colors) +
  labs(title = "CAFI Type Composition by Morphotype",
       x = "Morphotype",
       y = "Proportion",
       fill = "Type")

# C. Species accumulation by morphotype
if (nrow(community_matrix) > 10) {
  accum_results <- list()

  for (morph in unique(metadata$morphotype)) {
    if (!is.na(morph)) {
      morph_corals <- metadata %>%
        filter(morphotype == morph) %>%
        pull(coral_id)

      morph_comm <- community_matrix[rownames(community_matrix) %in% morph_corals, ]

      if (nrow(morph_comm) > 5) {
        spec_accum <- specaccum(morph_comm, method = "random", permutations = 100)
        accum_results[[morph]] <- data.frame(
          morphotype = morph,
          sites = spec_accum$sites,
          richness = spec_accum$richness,
          sd = spec_accum$sd
        )
      }
    }
  }

  if (length(accum_results) > 0) {
    accum_df <- bind_rows(accum_results)

    p3_accum <- ggplot(accum_df, aes(x = sites, y = richness, color = morphotype)) +
      geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd,
                     fill = morphotype), alpha = 0.2, color = NA) +
      geom_line(size = 1.5) +
      scale_color_manual(values = morphotype_colors) +
      scale_fill_manual(values = morphotype_colors) +
      labs(title = "OTU Accumulation by Morphotype",
           x = "Number of Corals Sampled",
           y = "OTU Richness",
           color = "Morphotype",
           fill = "Morphotype")
  }
}

# Combine morphotype plots
if (exists("p3_accum")) {
  morphotype_viz <- p3_ridges / (p3_interaction | p3_accum) +
    plot_annotation(title = "Morphotype-Specific Patterns")
} else {
  morphotype_viz <- p3_ridges / p3_interaction +
    plot_annotation(title = "Morphotype-Specific Patterns")
}

ggsave(file.path(fig_dir, "morphotype_patterns.png"),
       morphotype_viz, width = 16, height = 12, dpi = 300)

cat("✓ Morphotype visualizations created\n\n")

# ============================================================================
# 4. Spatial Patterns
# ============================================================================

cat("Creating spatial visualizations...\n")

# Check for coordinates
if (all(c("lat", "long") %in% names(metadata))) {
  # A. Spatial distribution of corals
  p4_spatial <- metadata %>%
    left_join(
      cafi_clean %>%
        group_by(coral_id) %>%
        summarise(total_cafi = n(), .groups = "drop"),
      by = "coral_id"
    ) %>%
    mutate(total_cafi = replace_na(total_cafi, 0)) %>%
    ggplot(aes(x = long, y = lat)) +
    geom_point(aes(size = total_cafi, color = morphotype), alpha = 0.6) +
    scale_size_continuous(range = c(2, 10), name = "Total CAFI") +
    scale_color_manual(values = morphotype_colors) +
    labs(title = "Spatial Distribution of Survey Corals",
         x = "Longitude",
         y = "Latitude") +
    coord_quickmap()

  # B. Site-level summaries
  site_summary <- metadata %>%
    left_join(
      cafi_clean %>%
        group_by(coral_id) %>%
        summarise(
          total_cafi = n(),
          species_richness = n_distinct(species),
          .groups = "drop"
        ),
      by = "coral_id"
    ) %>%
    mutate(across(c(total_cafi, species_richness), ~replace_na(., 0))) %>%
    group_by(site) %>%
    summarise(
      mean_lat = mean(lat, na.rm = TRUE),
      mean_long = mean(long, na.rm = TRUE),
      mean_abundance = mean(total_cafi),
      mean_richness = mean(species_richness),
      n_corals = n(),
      .groups = "drop"
    )

  p4_sites <- ggplot(site_summary, aes(x = mean_long, y = mean_lat)) +
    geom_point(aes(size = mean_abundance, color = mean_richness), alpha = 0.8) +
    geom_text_repel(aes(label = site), size = 3) +
    scale_size_continuous(range = c(5, 15), name = "Mean\nAbundance") +
    scale_color_viridis_c(name = "Mean\nRichness") +
    labs(title = "Site-Level Community Metrics",
         x = "Longitude",
         y = "Latitude") +
    coord_quickmap()

  spatial_viz <- p4_spatial | p4_sites
} else {
  # Non-spatial site comparison if no coordinates
  p4_sites <- summary_data %>%
    ggplot(aes(x = site, y = total_cafi)) +
    geom_boxplot(aes(fill = site), alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    scale_fill_manual(values = site_colors) +
    labs(title = "CAFI Abundance by Site",
         x = "Site",
         y = "Total CAFI") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  spatial_viz <- p4_sites
}

ggsave(file.path(fig_dir, "spatial_patterns.png"),
       spatial_viz, width = 14, height = 8, dpi = 300)

cat("✓ Spatial visualizations created\n\n")

# ============================================================================
# 5. Publication-Ready Figure
# ============================================================================

cat("Creating publication-ready figure...\n")

# Prepare data for publication figure
pub_data <- summary_data %>%
  filter(!is.na(morphotype))

# A. Main abundance comparison
p_pub_abundance <- ggplot(pub_data, aes(x = morphotype, y = sqrt(total_cafi))) +
  geom_violin(aes(fill = morphotype), alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  scale_fill_manual(values = morphotype_colors) +
  labs(x = "Coral Morphotype",
       y = expression(sqrt("Total CAFI"))) +
  theme_classic() +
  theme(legend.position = "none")

# B. Richness comparison
p_pub_richness <- ggplot(pub_data, aes(x = morphotype, y = species_richness)) +
  geom_violin(aes(fill = morphotype), alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  scale_fill_manual(values = morphotype_colors) +
  labs(x = "Coral Morphotype",
       y = "OTU Richness") +
  theme_classic() +
  theme(legend.position = "none")

# C. Top OTU distribution (community-level patterns only)
top10_species <- cafi_clean %>%
  count(species) %>%
  top_n(10, n) %>%
  pull(species)

species_morph_summary <- cafi_clean %>%
  filter(species %in% top10_species) %>%
  left_join(metadata %>% select(coral_id, morphotype), by = "coral_id") %>%
  filter(!is.na(morphotype)) %>%
  count(species, morphotype) %>%
  group_by(species) %>%
  mutate(total = sum(n),
         proportion = n / total) %>%
  arrange(desc(total))

p_pub_species <- ggplot(species_morph_summary,
                        aes(x = reorder(species, total), y = proportion,
                            fill = morphotype)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = morphotype_colors, name = "Morphotype") +
  coord_flip() +
  labs(x = "CAFI OTU",
       y = "Proportion of Occurrences") +
  theme_classic() +
  theme(legend.position = "bottom")

# D. Ordination
if (nrow(community_matrix) > 10) {
  set.seed(123)
  nmds <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)

  nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
    rownames_to_column("coral_id") %>%
    left_join(metadata %>% select(coral_id, morphotype), by = "coral_id") %>%
    filter(!is.na(morphotype))

  p_pub_ord <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(color = morphotype), size = 2, alpha = 0.6) +
    stat_ellipse(aes(color = morphotype), level = 0.95) +
    scale_color_manual(values = morphotype_colors, name = "Morphotype") +
    labs(subtitle = paste("Stress =", round(nmds$stress, 3))) +
    theme_classic() +
    theme(legend.position = "bottom")
}

# Combine for publication
if (exists("p_pub_ord")) {
  pub_figure <- (p_pub_abundance + p_pub_richness) /
                (p_pub_species + p_pub_ord) +
    plot_annotation(tag_levels = "A")
} else {
  pub_figure <- (p_pub_abundance + p_pub_richness) /
                p_pub_species +
    plot_annotation(tag_levels = "A")
}

ggsave(file.path(fig_dir, "publication_figure.png"),
       pub_figure, width = 12, height = 12, dpi = 300)

ggsave(file.path(fig_dir, "publication_figure.pdf"),
       pub_figure, width = 12, height = 12)

cat("✓ Publication figure created\n\n")

# ============================================================================
# 6. Interactive Summary Table
# ============================================================================

cat("Creating summary tables...\n")

# Create comprehensive summary table
summary_table <- summary_data %>%
  group_by(morphotype, site) %>%
  summarise(
    n_corals = n(),
    mean_abundance = round(mean(total_cafi), 1),
    sd_abundance = round(sd(total_cafi), 1),
    mean_richness = round(mean(species_richness), 1),
    sd_richness = round(sd(species_richness), 1),
    mean_shannon = round(mean(shannon), 2),
    sd_shannon = round(sd(shannon), 2),
    .groups = "drop"
  ) %>%
  mutate(
    abundance = paste0(mean_abundance, " ± ", sd_abundance),
    richness = paste0(mean_richness, " ± ", sd_richness),
    shannon = paste0(mean_shannon, " ± ", sd_shannon)
  ) %>%
  select(morphotype, site, n_corals, abundance, richness, shannon)

write_csv(summary_table,
          file.path(SURVEY_TABLES, "comprehensive_summary_table.csv"))

# OTU occurrence table
species_table <- cafi_clean %>%
  count(species, type) %>%
  arrange(desc(n)) %>%
  mutate(
    rank = row_number(),
    percentage = round(100 * n / sum(n), 2)
  ) %>%
  select(rank, species, type, n, percentage)

write_csv(species_table,
          file.path(SURVEY_TABLES, "otu_occurrence_table.csv"))

cat("✓ Summary tables created\n\n")

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Visualization Suite Summary\n")
cat("========================================\n\n")

cat("Visualizations Created:\n")
cat("  1. Overview dashboard (4 panels)\n")
cat("  2. Community composition analysis\n")
cat("  3. Morphotype-specific patterns\n")
cat("  4. Spatial distribution maps\n")
cat("  5. Publication-ready figure\n")
cat("  6. Summary tables\n\n")

cat("Key Findings:\n")
cat("  - Total corals analyzed:", nrow(metadata), "\n")
cat("  - Total CAFI individuals:", nrow(cafi_clean), "\n")
cat("  - Number of OTUs (morphotype-based):", n_distinct(cafi_clean$species), "\n")
cat("  - Mean abundance per coral:", round(mean(summary_data$total_cafi), 1), "\n")
cat("  - Mean OTU richness:", round(mean(summary_data$species_richness), 1), "\n\n")

cat("✅ Visualization suite complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", SURVEY_TABLES, "\n")