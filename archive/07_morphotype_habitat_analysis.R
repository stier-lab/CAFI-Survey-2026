#!/usr/bin/env Rscript
# ============================================================================
# 07_morphotype_habitat_analysis.R - Detailed morphotype and microhabitat analysis
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
# ============================================================================

cat("\n========================================\n")
cat("Morphotype and Microhabitat Analysis\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/00_load_libraries.R"))
library(ggridges)
library(ggalluvial)

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))

# Add branch_width as NA if it doesn't exist (for compatibility)
if(!"branch_width" %in% names(metadata)) {
  metadata$branch_width <- NA_character_
  cat("Note: branch_width column not in data - using NA values\n")
}

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "morphotype_habitat")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Morphotype-Specific Community Structure
# ============================================================================

cat("Analyzing morphotype-specific communities...\n")

# Calculate community metrics by morphotype
# Note: Focus on functional groups (crabs, shrimp, fish, snails) rather than individual OTUs
morphotype_communities <- cafi_clean %>%
  left_join(metadata %>% select(coral_id, morphotype, any_of("branch_width")), by = "coral_id") %>%
  filter(!is.na(morphotype)) %>%
  group_by(morphotype, species, type) %>%
  summarise(
    abundance = n(),
    mean_size = mean(size_mm, na.rm = TRUE),
    occurrence_freq = n_distinct(coral_id),
    .groups = "drop"
  )

# Calculate distribution patterns (community-level, not ecological specialization)
# Without genetic IDs, we cannot infer true ecological preferences
total_by_species <- morphotype_communities %>%
  group_by(species) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop")

morphotype_distribution_patterns <- morphotype_communities %>%
  left_join(total_by_species, by = "species") %>%
  mutate(
    proportion = abundance / total_abundance,
    distribution_entropy = -proportion * log(proportion)
  ) %>%
  group_by(species) %>%
  summarise(
    shannon_distribution = sum(distribution_entropy),
    n_morphotypes = n(),
    most_common_morphotype = morphotype[which.max(proportion)],
    max_proportion = max(proportion),
    .groups = "drop"
  ) %>%
  mutate(
    distribution_score = 1 - (shannon_distribution / log(n_morphotypes)),
    distribution_pattern = case_when(
      max_proportion > 0.75 ~ "Concentrated",
      max_proportion > 0.5 ~ "Skewed",
      TRUE ~ "Even"
    )
  )

write_csv(morphotype_distribution_patterns,
          file.path(SURVEY_TABLES, "morphotype_distribution_patterns.csv"))

# ============================================================================
# Microhabitat Utilization Patterns
# ============================================================================

cat("Analyzing microhabitat utilization...\n")

# Define microhabitat based on morphotype + branch width
microhabitat_data <- metadata %>%
  mutate(
    microhabitat = paste(morphotype, branch_width, sep = "_"),
    habitat_complexity = case_when(
      morphotype %in% c("verrucosa", "cespitosa") & branch_width == "tight" ~ "High",
      morphotype == "meandrina" ~ "Medium",
      morphotype == "plating" ~ "Low",
      TRUE ~ "Medium"
    )
  )

# CAFI distribution across microhabitats
microhabitat_communities <- cafi_clean %>%
  left_join(microhabitat_data %>% select(coral_id, microhabitat, habitat_complexity),
            by = "coral_id") %>%
  filter(!is.na(microhabitat))

# Calculate microhabitat distribution breadth (abundance patterns, not true niche breadth)
# Without genetic species IDs, these are community distribution metrics
microhabitat_distribution <- microhabitat_communities %>%
  group_by(species, microhabitat) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(
    prop = count / sum(count),
    total_count = sum(count)
  ) %>%
  summarise(
    levins_B = 1 / sum(prop^2),
    n_microhabitats = n(),
    standardized_B = (levins_B - 1) / (n_microhabitats - 1),
    total_abundance = first(total_count),
    .groups = "drop"
  )

write_csv(microhabitat_distribution,
          file.path(SURVEY_TABLES, "otu_microhabitat_distribution.csv"))

# ============================================================================
# Morphotype × Branch Width Interaction
# ============================================================================

cat("Analyzing morphotype × branch width interactions...\n")

# Two-way analysis
morphotype_branch_matrix <- microhabitat_data %>%
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
  mutate(
    total_cafi = replace_na(total_cafi, 0),
    species_richness = replace_na(species_richness, 0)
  )

# Interaction plot
p_interaction <- morphotype_branch_matrix %>%
  filter(!is.na(morphotype), !is.na(branch_width)) %>%
  group_by(morphotype, branch_width) %>%
  summarise(
    mean_cafi = mean(total_cafi),
    se_cafi = sd(total_cafi) / sqrt(n()),
    mean_richness = mean(species_richness),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = branch_width, y = mean_cafi, color = morphotype, group = morphotype)) +
  geom_errorbar(aes(ymin = mean_cafi - se_cafi, ymax = mean_cafi + se_cafi),
                width = 0.1, position = position_dodge(0.1)) +
  geom_line(size = 1.5, position = position_dodge(0.1)) +
  geom_point(size = 3, position = position_dodge(0.1)) +
  scale_color_viridis_d() +
  labs(
    title = "CAFI Abundance: Morphotype × Branch Width Interaction",
    x = "Branch Width",
    y = "Mean CAFI per Coral (±SE)",
    color = "Morphotype"
  )

ggsave(file.path(fig_dir, "morphotype_branch_interaction.png"),
       p_interaction, width = 10, height = 6, dpi = 300)

# ============================================================================
# Habitat Complexity Gradient
# ============================================================================

cat("Analyzing habitat complexity gradient...\n")

# Complexity effects on community
complexity_effects <- microhabitat_data %>%
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
  mutate(
    habitat_complexity = factor(habitat_complexity, levels = c("Low", "Medium", "High")),
    total_cafi = replace_na(total_cafi, 0),
    species_richness = replace_na(species_richness, 0)
  )

# Ridge plot of abundance distribution
p_complexity_ridge <- ggplot(complexity_effects,
                             aes(x = total_cafi, y = habitat_complexity, fill = habitat_complexity)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5, rel_min_height = 0.01) +
  scale_fill_viridis_d() +
  scale_x_sqrt() +
  labs(
    title = "CAFI Abundance Distribution Across Habitat Complexity",
    x = "Total CAFI per Coral (sqrt scale)",
    y = "Habitat Complexity"
  ) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "habitat_complexity_distributions.png"),
       p_complexity_ridge, width = 10, height = 6, dpi = 300)

# ============================================================================
# Species-Morphotype Association Network
# ============================================================================

cat("Building species-morphotype association network...\n")

# Calculate association strength
morphotype_associations <- morphotype_communities %>%
  group_by(morphotype) %>%
  mutate(morphotype_total = sum(abundance)) %>%
  ungroup() %>%
  mutate(
    expected = (sum(abundance) * morphotype_total) / sum(abundance)^2,
    association_strength = (abundance - expected) / sqrt(expected),
    significant = abs(association_strength) > 2
  ) %>%
  filter(significant)

# Create bipartite network visualization
if (nrow(morphotype_associations) > 0) {
  # Prepare data for alluvial plot
  alluvial_data <- morphotype_associations %>%
    select(morphotype, species, abundance, type) %>%
    group_by(morphotype, type) %>%
    slice_max(order_by = abundance, n = 5) %>%
    ungroup()

  p_alluvial <- ggplot(alluvial_data,
                       aes(axis1 = morphotype, axis2 = species, y = abundance)) +
    geom_alluvium(aes(fill = type), width = 1/12) +
    geom_stratum(width = 1/8, fill = "grey90", color = "black") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              size = 3, angle = 0) +
    scale_fill_viridis_d() +
    labs(
      title = "Species-Morphotype Associations",
      subtitle = "Top 5 species per taxonomic group",
      y = "Abundance",
      fill = "Taxonomic Group"
    ) +
    theme_void() +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, "species_morphotype_alluvial.png"),
         p_alluvial, width = 12, height = 8, dpi = 300)
}

# ============================================================================
# Indicator OTUs by Morphotype (Community Patterns)
# ============================================================================

cat("Identifying morphotype-associated OTUs (community-level patterns)...\n")

# Indicator OTU analysis for each morphotype
# NOTE: These identify abundance associations, not ecological specialization
morphotype_indicators <- list()

for (morph in unique(metadata$morphotype[!is.na(metadata$morphotype)])) {
  # Get corals of this morphotype
  morph_corals <- metadata %>%
    filter(morphotype == morph) %>%
    pull(coral_id)

  # Create indicator variable
  indicator_var <- ifelse(rownames(community_matrix) %in% morph_corals, 1, 0)

  # Run indicator analysis if enough samples
  if (sum(indicator_var == 1) >= 5 && sum(indicator_var == 0) >= 5) {
    indval_result <- multipatt(community_matrix, indicator_var,
                               control = how(nperm = 999))

    # Extract significant indicators
    sig_indicators <- indval_result$sign %>%
      filter(p.value < 0.05) %>%
      mutate(
        species = rownames(.),
        morphotype = morph
      )

    if (nrow(sig_indicators) > 0) {
      morphotype_indicators[[morph]] <- sig_indicators
    }
  }
}

# Combine all indicators
if (length(morphotype_indicators) > 0) {
  all_indicators <- bind_rows(morphotype_indicators)
  write_csv(all_indicators,
            file.path(SURVEY_TABLES, "morphotype_associated_otus.csv"))

  # Plot indicator values
  p_indicators <- all_indicators %>%
    ggplot(aes(x = reorder(species, stat), y = stat, fill = morphotype)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~morphotype, scales = "free_y") +
    scale_fill_viridis_d() +
    labs(
      title = "Morphotype-Associated OTUs",
      subtitle = "Abundance associations (community patterns)",
      x = "CAFI OTU",
      y = "Indicator Value",
      fill = "Morphotype"
    )

  ggsave(file.path(fig_dir, "morphotype_associated_otus.png"),
         p_indicators, width = 12, height = 10, dpi = 300)
}

# ============================================================================
# Morphotype Distribution Profiles (Community-Level)
# ============================================================================

cat("Creating morphotype distribution profiles...\n")

# Calculate distribution patterns for common OTUs (abundance only, not true preferences)
top_species_distribution <- morphotype_communities %>%
  group_by(species) %>%
  filter(sum(abundance) >= 20) %>%  # Focus on common OTUs
  group_by(species, type) %>%
  mutate(
    total_sp_abundance = sum(abundance),
    distribution_score = abundance / total_sp_abundance
  ) %>%
  ungroup()

# Heatmap of distribution patterns
# Aggregate by species and morphotype first to avoid duplicate rows
distribution_matrix <- top_species_distribution %>%
  group_by(species, morphotype) %>%
  summarise(distribution_score = sum(abundance) / sum(total_sp_abundance), .groups = "drop") %>%
  pivot_wider(names_from = morphotype, values_from = distribution_score, values_fill = 0) %>%
  column_to_rownames("species") %>%
  as.matrix()

if (nrow(distribution_matrix) > 5) {
  png(file.path(fig_dir, "morphotype_distribution_heatmap.png"),
      width = 12, height = 10, units = "in", res = 300)
  heatmap(distribution_matrix,
          col = viridis(100),
          main = "OTU Distribution Across Morphotypes",
          sub = "Abundance patterns (community-level)",
          xlab = "Morphotype",
          ylab = "CAFI OTU",
          cexRow = 0.7,
          cexCol = 1)
  dev.off()
}

# ============================================================================
# Microhabitat Diversity Patterns
# ============================================================================

cat("Analyzing microhabitat diversity patterns...\n")

# Calculate alpha and beta diversity within morphotypes
morphotype_diversity <- metadata %>%
  filter(!is.na(morphotype)) %>%
  group_by(morphotype) %>%
  summarise(
    n_corals = n(),
    .groups = "drop"
  ) %>%
  filter(n_corals >= 5)

diversity_by_morphotype <- list()

for (morph in morphotype_diversity$morphotype) {
  morph_corals <- metadata %>%
    filter(morphotype == morph) %>%
    pull(coral_id)

  if (length(morph_corals) > 5) {
    morph_comm <- community_matrix[rownames(community_matrix) %in% morph_corals, ]
    morph_comm <- morph_comm[, colSums(morph_comm) > 0]

    if (nrow(morph_comm) > 2) {
      # Calculate diversity metrics
      alpha_div <- mean(vegan::diversity(morph_comm))
      beta_div <- mean(vegdist(morph_comm, method = "bray"))
      gamma_div <- vegan::diversity(colSums(morph_comm))

      diversity_by_morphotype[[morph]] <- data.frame(
        morphotype = morph,
        n_corals = length(morph_corals),
        mean_alpha = alpha_div,
        mean_beta = beta_div,
        gamma = gamma_div,
        beta_ratio = beta_div / alpha_div
      )
    }
  }
}

if (length(diversity_by_morphotype) > 0) {
  morphotype_div_summary <- bind_rows(diversity_by_morphotype)
  write_csv(morphotype_div_summary,
            file.path(SURVEY_TABLES, "morphotype_diversity_partitioning.csv"))

  # Plot diversity partitioning
  p_div_partition <- morphotype_div_summary %>%
    pivot_longer(cols = c(mean_alpha, mean_beta, gamma),
                 names_to = "diversity_component",
                 values_to = "value") %>%
    ggplot(aes(x = morphotype, y = value, fill = diversity_component)) +
    geom_col(position = "dodge") +
    scale_fill_viridis_d(labels = c("Alpha", "Beta", "Gamma")) +
    labs(
      title = "Diversity Partitioning by Morphotype",
      x = "Morphotype",
      y = "Diversity Value",
      fill = "Component"
    )

  ggsave(file.path(fig_dir, "morphotype_diversity_partitioning.png"),
         p_div_partition, width = 10, height = 6, dpi = 300)
}

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Morphotype and Microhabitat Analysis Summary\n")
cat("========================================\n\n")

cat("Morphotype Analysis:\n")
cat("  - Morphotypes analyzed:", length(unique(metadata$morphotype[!is.na(metadata$morphotype)])), "\n")
cat("  - Total microhabitats:", length(unique(microhabitat_data$microhabitat)), "\n\n")

if (exists("morphotype_distribution_patterns")) {
  cat("Distribution Patterns:\n")
  dist_summary <- table(morphotype_distribution_patterns$distribution_pattern)
  cat("  - Concentrated:", dist_summary["Concentrated"], "OTUs\n")
  cat("  - Skewed:", dist_summary["Skewed"], "OTUs\n")
  cat("  - Even:", dist_summary["Even"], "OTUs\n\n")
}

if (exists("microhabitat_distribution")) {
  cat("Microhabitat Distribution:\n")
  cat("  - Mean Levins' B:", round(mean(microhabitat_distribution$levins_B), 2), "\n")
  cat("  - Most concentrated:", microhabitat_distribution$species[which.min(microhabitat_distribution$levins_B)],
      "(B =", round(min(microhabitat_distribution$levins_B), 2), ")\n")
  cat("  - Most widespread:", microhabitat_distribution$species[which.max(microhabitat_distribution$levins_B)],
      "(B =", round(max(microhabitat_distribution$levins_B), 2), ")\n\n")
}

if (exists("all_indicators")) {
  cat("Morphotype-Associated OTUs:\n")
  cat("  - Total OTUs identified:", nrow(all_indicators), "\n")
  indicator_counts <- table(all_indicators$morphotype)
  for (morph in names(indicator_counts)) {
    cat("  -", morph, ":", indicator_counts[morph], "OTUs\n")
  }
}

cat("\n✅ Morphotype and microhabitat analysis complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", SURVEY_TABLES, "\n")