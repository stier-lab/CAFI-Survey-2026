#!/usr/bin/env Rscript
# ============================================================================
# 02_community_composition.R - Community Composition Analysis (H1)
#
# Hypothesis H1: CAFI community composition differs among reef sites due to
# variation in coral landscapes and environmental conditions.
#
# Theoretical Background:
#   Propagule redirection occurs at multiple scales. At the landscape scale,
#   sites with different coral abundance, size distributions, and spatial
#   configurations should support different CAFI communities due to variation
#   in settler supply and post-settlement processes.
#
# Tests:
#   - PERMANOVA on Bray-Curtis dissimilarity (site effect)
#   - NMDS ordination visualization
#   - Pairwise site comparisons
#
# IMPORTANT TAXONOMIC NOTES:
# - CAFI "species" are morphological OTUs (Operational Taxonomic Units)
# - NO genetic haplotyping was performed - these are field identifications
# - Use for richness/diversity metrics but NOT species-level inferences
# - Coral "morphotypes" (meandrina/eydoxi/verucosa) are NOT confirmed species
# - branch_width (tight/wide) is the actual measurable coral trait
#
# Sites: HAU = Hauru, MAT = Maatea, MRB = Moorea Barrier Reef
#
# Author: CAFI Analysis Pipeline
# Date: 2025-10-31
# ============================================================================

cat("\n========================================\n")
cat("Survey Community Composition Analysis\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/00_load_libraries.R"))

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "community_composition")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Overall Community Composition
# ============================================================================

cat("Analyzing overall community composition...\n")

# CAFI OTU abundance ranking
# NOTE: "species" here refers to morphological OTUs (Operational Taxonomic Units)
# These are field-identified morphotypes, NOT genetically confirmed species
# Use for community metrics (richness, diversity) but interpret cautiously
species_abundance <- cafi_clean %>%
  group_by(species, type) %>%
  summarise(
    total_abundance = n(),
    n_corals = n_distinct(coral_id),
    mean_size = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_abundance)) %>%
  mutate(
    rank = row_number(),
    proportion = total_abundance / sum(total_abundance),
    cumulative_prop = cumsum(proportion)
  )

# Save species abundance table
write_csv(species_abundance,
          file.path(SURVEY_TABLES, "species_abundance_ranking.csv"))

# Top 20 most abundant species
top20_species <- species_abundance %>%
  slice_head(n = 20)

# Plot CAFI OTU rank abundance
p_rank_abundance <- ggplot(species_abundance, aes(x = rank, y = total_abundance)) +
  geom_point(aes(color = type), size = 3, alpha = 0.7) +
  geom_line(alpha = 0.3) +
  scale_y_log10(labels = scales::comma) +
  scale_color_taxon() +
  labs(
    title = "CAFI OTU Rank Abundance Curve",
    subtitle = "Field Survey Summer 2019 (morphotype-based groupings)",
    x = "OTU Rank",
    y = "Total Abundance (log scale)",
    color = "Functional\nGroup"
  ) +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "species_rank_abundance.png"),
       p_rank_abundance, width = 10, height = 6, dpi = 300)

# ============================================================================
# Community Composition by Site
# ============================================================================

cat("Analyzing site-level patterns...\n")

# Site-level community summary
# Comparing 3 major sites:
# - HAU (Hauru): 38 corals, north shore fringing reef
# - MAT (Maatea): 39 corals, interior lagoon
# - MRB (Moorea Barrier Reef): 37 corals, outer barrier reef
site_summary <- cafi_clean %>%
  group_by(site, species, type) %>%
  summarise(
    abundance = n(),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  mutate(
    total_site_abundance = sum(abundance),
    proportion = abundance / total_site_abundance
  )

# Species richness by site
site_richness <- cafi_clean %>%
  group_by(site) %>%
  summarise(
    n_corals = n_distinct(coral_id),
    species_richness = n_distinct(species),
    total_abundance = n(),
    shannon = vegan::diversity(table(species)),
    simpson = vegan::diversity(table(species), index = "simpson"),
    .groups = "drop"
  )

write_csv(site_richness,
          file.path(SURVEY_TABLES, "site_diversity_metrics.csv"))

# Plot site diversity
p_site_diversity <- site_richness %>%
  pivot_longer(cols = c(species_richness, shannon, simpson),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("species_richness", "shannon", "simpson"),
                         labels = c("Species Richness", "Shannon Diversity", "Simpson Diversity"))) %>%
  ggplot(aes(x = site, y = value, fill = site)) +
  geom_col() +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_site() +
  labs(
    title = "Diversity Metrics by Site",
    x = "Site",
    y = "Value"
  ) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "site_diversity_metrics.png"),
       p_site_diversity, width = 12, height = 6, dpi = 300)

# ============================================================================
# Taxonomic Composition
# ============================================================================

cat("Analyzing taxonomic composition...\n")

# Taxonomic breakdown
taxonomic_summary <- cafi_clean %>%
  group_by(type) %>%
  summarise(
    total_abundance = n(),
    n_species = n_distinct(species),
    n_corals = n_distinct(coral_id),
    mean_size = mean(size_mm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(proportion = total_abundance / sum(total_abundance))

write_csv(taxonomic_summary,
          file.path(SURVEY_TABLES, "taxonomic_composition.csv"))

# Plot taxonomic composition pie chart
p_taxonomic_pie <- ggplot(taxonomic_summary,
                          aes(x = "", y = proportion, fill = type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_taxon() +
  labs(
    title = "Taxonomic Composition of CAFI Community",
    fill = "Taxonomic Group"
  ) +
  theme_void() +
  theme(legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA)) +
  geom_text(aes(label = paste0(round(proportion * 100, 1), "%")),
            position = position_stack(vjust = 0.5))

ggsave(file.path(fig_dir, "taxonomic_composition_pie.png"),
       p_taxonomic_pie, width = 8, height = 6, dpi = 300)

# ============================================================================
# Top Species Analysis
# ============================================================================

cat("Analyzing dominant species...\n")

# Top 10 species by site
top10_by_site <- site_summary %>%
  group_by(site) %>%
  slice_max(order_by = abundance, n = 10) %>%
  ungroup()

# Heatmap of top species across sites
top_species_matrix <- top20_species$species
site_species_matrix <- site_summary %>%
  filter(species %in% top_species_matrix) %>%
  select(site, species, proportion) %>%
  pivot_wider(names_from = site,
              values_from = proportion,
              values_fn = sum,
              values_fill = 0)

# Create heatmap
heatmap_data <- site_species_matrix %>%
  column_to_rownames("species") %>%
  as.matrix()

png(file.path(fig_dir, "top_species_heatmap.png"),
    width = 10, height = 8, units = "in", res = 300, bg = "white")
heatmap(heatmap_data,
        scale = "row",
        main = "Top 20 CAFI OTU Distribution Across Sites",
        sub = "(Morphotype-based groupings, not genetic species)",
        xlab = "Site",
        ylab = "CAFI OTU",
        col = viridis(100))
dev.off()

# ============================================================================
# Species Accumulation Curves
# ============================================================================

cat("Calculating species accumulation curves...\n")

# Prepare data for accumulation curves
comm_by_site <- list()
sites <- unique(metadata$site)

for (s in sites) {
  site_corals <- metadata %>%
    filter(site == s) %>%
    pull(coral_id)

  if (length(site_corals) > 0) {
    comm_by_site[[s]] <- community_matrix[rownames(community_matrix) %in% site_corals, ]
    # Remove empty species columns
    comm_by_site[[s]] <- comm_by_site[[s]][, colSums(comm_by_site[[s]]) > 0]
  }
}

# Calculate species accumulation
spec_accum_list <- list()
for (s in names(comm_by_site)) {
  if (nrow(comm_by_site[[s]]) > 1) {
    spec_accum_list[[s]] <- specaccum(comm_by_site[[s]], method = "random", permutations = 100)
  }
}

# Plot species accumulation curves
png(file.path(fig_dir, "species_accumulation_curves.png"),
    width = 10, height = 6, units = "in", res = 300, bg = "white")
par(mfrow = c(1, 1), bg = "white")
plot(1, type = "n", xlim = c(0, max(sapply(spec_accum_list, function(x) max(x$sites)))),
     ylim = c(0, max(sapply(spec_accum_list, function(x) max(x$richness)))),
     xlab = "Number of Corals Sampled",
     ylab = "Species Richness",
     main = "Species Accumulation Curves by Site")

colors <- viridis(length(spec_accum_list))
for (i in seq_along(spec_accum_list)) {
  if (!is.null(spec_accum_list[[i]])) {
    lines(spec_accum_list[[i]]$sites, spec_accum_list[[i]]$richness,
          col = colors[i], lwd = 2)
  }
}
legend("bottomright", names(spec_accum_list), col = colors, lwd = 2)
dev.off()

# ============================================================================
# Community Structure Visualization
# ============================================================================

cat("Visualizing community structure...\n")

# Stacked bar chart of community composition by site
p_community_stack <- site_summary %>%
  group_by(site) %>%
  slice_max(order_by = proportion, n = 15) %>%
  ggplot(aes(x = site, y = proportion, fill = species)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d() +
  labs(
    title = "Community Composition by Site",
    subtitle = "Top 15 species per site",
    x = "Site",
    y = "Proportion of Community",
    fill = "Species"
  ) +
  theme_publication() +
  theme(legend.position = "right",
        legend.text = element_text(size = 8))

ggsave(file.path(fig_dir, "community_composition_stacked.png"),
       p_community_stack, width = 12, height = 8, dpi = 300)

# ============================================================================
# Size Distribution Analysis
# ============================================================================

cat("Analyzing size distributions...\n")

# Size distribution by taxonomic group
p_size_distribution <- cafi_clean %>%
  filter(!is.na(size_mm)) %>%
  ggplot(aes(x = size_mm, fill = type)) +
  geom_histogram(binwidth = 2, alpha = 0.7, position = "identity") +
  facet_wrap(~type, scales = "free_y") +
  scale_fill_taxon() +
  labs(
    title = "Size Distribution of CAFI by Taxonomic Group",
    x = "Size (mm)",
    y = "Count",
    fill = "Type"
  ) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "size_distribution_by_type.png"),
       p_size_distribution, width = 12, height = 8, dpi = 300)

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Community Composition Summary\n")
cat("========================================\n\n")

cat("Overall Community:\n")
cat("  - Total species:", n_distinct(cafi_clean$species), "\n")
cat("  - Total individuals:", nrow(cafi_clean), "\n")
cat("  - Most abundant species:", top20_species$species[1],
    "(", top20_species$total_abundance[1], "individuals)\n")
cat("  - Dominant taxonomic group:", taxonomic_summary$type[1],
    "(", round(taxonomic_summary$proportion[1] * 100, 1), "%)\n\n")

cat("Site Comparison:\n")
cat("  - Highest richness:", site_richness$site[which.max(site_richness$species_richness)],
    "(", max(site_richness$species_richness), "species)\n")
cat("  - Highest abundance:", site_richness$site[which.max(site_richness$total_abundance)],
    "(", max(site_richness$total_abundance), "individuals)\n")
cat("  - Highest diversity:", site_richness$site[which.max(site_richness$shannon)],
    "(H' =", round(max(site_richness$shannon), 2), ")\n\n")

cat("✅ Community composition analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", SURVEY_TABLES, "\n")