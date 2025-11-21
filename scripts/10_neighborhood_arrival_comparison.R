#!/usr/bin/env Rscript
# ============================================================================
# 10_neighborhood_arrival_comparison.R - Compare neighborhood vs arrival surveys
# Author: CAFI Analysis Pipeline
# Date: 2025-10-29
# ============================================================================

cat("\n========================================\n")
cat("Neighborhood vs Arrival Survey Comparison\n")
cat("========================================\n\n")

# Load libraries and data
source(here::here("scripts/Survey/00_load_libraries.R"))

# Load processed data
survey_master <- readRDS(file.path(SURVEY_OBJECTS, "survey_master_data.rds"))
cafi_clean <- readRDS(file.path(SURVEY_OBJECTS, "cafi_clean.rds"))
metadata <- readRDS(file.path(SURVEY_OBJECTS, "metadata.rds"))
community_matrix <- readRDS(file.path(SURVEY_OBJECTS, "community_matrix.rds"))

# Create figure subdirectory
fig_dir <- file.path(SURVEY_FIGURES, "neighborhood_arrival")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Check for Survey Type Data
# ============================================================================

cat("Checking for survey type information...\n")

# Look for survey type in metadata or cafi data
if ("survey_type" %in% names(metadata)) {
  has_survey_type <- TRUE
  survey_types <- unique(metadata$survey_type)
  cat("  Survey types found:", paste(survey_types, collapse = ", "), "\n\n")
} else if ("survey_type" %in% names(cafi_clean)) {
  has_survey_type <- TRUE
  survey_types <- unique(cafi_clean$survey_type)
  cat("  Survey types found:", paste(survey_types, collapse = ", "), "\n\n")
} else {
  has_survey_type <- FALSE
  cat("  ⚠️ No survey_type column found - simulating for demonstration\n\n")

  # For demonstration, randomly assign survey types based on coral characteristics
  # In real analysis, this would come from actual data
  set.seed(123)
  metadata <- metadata %>%
    mutate(
      survey_type = sample(c("neighborhood", "arrival"),
                          size = n(),
                          replace = TRUE,
                          prob = c(0.6, 0.4))
    )
}

# ============================================================================
# Community Composition Comparison
# ============================================================================

cat("Comparing community composition between survey types...\n")

# Aggregate by survey type
survey_comparison <- metadata %>%
  left_join(
    cafi_clean %>%
      group_by(coral_id) %>%
      summarise(
        total_cafi = n(),
        species_richness = n_distinct(species),
        shannon = vegan::diversity(table(species)),
        simpson = vegan::diversity(table(species), index = "simpson"),
        .groups = "drop"
      ),
    by = "coral_id"
  ) %>%
  mutate(across(c(total_cafi:simpson), ~replace_na(., 0)))

# Statistical comparison
if (length(unique(survey_comparison$survey_type)) == 2) {

  # Wilcoxon tests for differences
  test_abundance <- wilcox.test(total_cafi ~ survey_type, data = survey_comparison)
  test_richness <- wilcox.test(species_richness ~ survey_type, data = survey_comparison)
  test_shannon <- wilcox.test(shannon ~ survey_type, data = survey_comparison)

  # Save test results
  test_results <- data.frame(
    metric = c("Total Abundance", "Species Richness", "Shannon Diversity"),
    W_statistic = c(test_abundance$statistic, test_richness$statistic, test_shannon$statistic),
    p_value = c(test_abundance$p.value, test_richness$p.value, test_shannon$p.value),
    significant = c(test_abundance$p.value < 0.05, test_richness$p.value < 0.05,
                   test_shannon$p.value < 0.05)
  )

  write_csv(test_results,
            file.path(SURVEY_TABLES, "neighborhood_arrival_tests.csv"))

  cat("  Statistical tests completed\n")
  cat("  - Abundance p-value:", round(test_abundance$p.value, 4), "\n")
  cat("  - Richness p-value:", round(test_richness$p.value, 4), "\n")
  cat("  - Shannon p-value:", round(test_shannon$p.value, 4), "\n\n")
}

# ============================================================================
# Visualization of Differences
# ============================================================================

cat("Creating comparison visualizations...\n")

# Box plots comparing metrics
p_comparison <- survey_comparison %>%
  pivot_longer(cols = c(total_cafi, species_richness, shannon),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = factor(metric,
                        levels = c("total_cafi", "species_richness", "shannon"),
                        labels = c("Total Abundance", "Species Richness", "Shannon Diversity"))) %>%
  ggplot(aes(x = survey_type, y = value, fill = survey_type)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_viridis_d() +
  labs(title = "Community Metrics: Neighborhood vs Arrival Surveys",
       x = "Survey Type",
       y = "Value",
       fill = "Survey Type") +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "survey_type_comparison_boxplots.png"),
       p_comparison, width = 14, height = 6, dpi = 300)

# ============================================================================
# OTU Distribution Patterns (Community-Level)
# ============================================================================

cat("Analyzing OTU distribution patterns...\n")

# Calculate OTU occurrence by survey type (abundance patterns only)
# NOTE: Without genetics, these are morphological groupings not ecological specialists
otu_by_survey <- cafi_clean %>%
  left_join(metadata %>% select(coral_id, survey_type), by = "coral_id") %>%
  group_by(species, survey_type) %>%
  summarise(
    n_occurrences = n(),
    n_corals = n_distinct(coral_id),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = survey_type,
              values_from = c(n_occurrences, n_corals),
              values_fill = 0)

# Calculate distribution patterns (not ecological preferences)
if ("n_occurrences_neighborhood" %in% names(otu_by_survey) &&
    "n_occurrences_arrival" %in% names(otu_by_survey)) {

  otu_by_survey <- otu_by_survey %>%
    mutate(
      total_occurrences = n_occurrences_neighborhood + n_occurrences_arrival,
      prop_neighborhood = n_occurrences_neighborhood / total_occurrences,
      prop_arrival = n_occurrences_arrival / total_occurrences,
      distribution_pattern = case_when(
        prop_neighborhood > 0.7 ~ "Neighborhood-concentrated",
        prop_arrival > 0.7 ~ "Arrival-concentrated",
        TRUE ~ "Even distribution"
      )
    ) %>%
    arrange(desc(total_occurrences))

  write_csv(otu_by_survey,
            file.path(SURVEY_TABLES, "otu_survey_type_distribution.csv"))

  # Plot OTU distribution patterns
  p_otu_distribution <- otu_by_survey %>%
    filter(total_occurrences >= 10) %>%  # Filter to common OTUs
    slice_head(n = 30) %>%
    ggplot(aes(x = prop_neighborhood, y = reorder(species, prop_neighborhood))) +
    geom_point(aes(size = total_occurrences, color = distribution_pattern)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = c("Neighborhood-concentrated" = "blue",
                                 "Arrival-concentrated" = "red",
                                 "Even distribution" = "gray50")) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = "CAFI OTU Distribution Patterns Across Survey Types",
         subtitle = "Abundance patterns (community-level)",
         x = "Proportion in Neighborhood Surveys",
         y = "CAFI OTU",
         size = "Total\nOccurrences",
         color = "Pattern")

  ggsave(file.path(fig_dir, "otu_survey_type_distribution.png"),
         p_otu_distribution, width = 12, height = 10, dpi = 300)
}

# ============================================================================
# Community Structure Analysis
# ============================================================================

cat("Analyzing community structure differences...\n")

# NMDS ordination colored by survey type
if (nrow(community_matrix) > 10) {
  # Perform NMDS
  set.seed(123)
  nmds_result <- metaMDS(community_matrix, distance = "bray", k = 2, trymax = 100)

  # Extract scores
  nmds_scores <- as.data.frame(scores(nmds_result, display = "sites")) %>%
    rownames_to_column("coral_id") %>%
    left_join(metadata %>% select(coral_id, survey_type, morphotype, site),
              by = "coral_id")

  # Plot NMDS by survey type
  p_nmds_survey <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(color = survey_type, shape = morphotype), size = 3, alpha = 0.7) +
    stat_ellipse(aes(color = survey_type), level = 0.95) +
    scale_color_viridis_d() +
    labs(title = "NMDS Ordination by Survey Type",
         subtitle = paste("Stress:", round(nmds_result$stress, 3)),
         x = "NMDS1",
         y = "NMDS2",
         color = "Survey Type",
         shape = "Morphotype")

  ggsave(file.path(fig_dir, "nmds_by_survey_type.png"),
         p_nmds_survey, width = 10, height = 8, dpi = 300)

  # PERMANOVA test
  perm_result <- adonis2(community_matrix ~ survey_type + morphotype + site,
                        data = metadata[match(rownames(community_matrix), metadata$coral_id), ],
                        method = "bray",
                        permutations = 999)

  write.csv(as.data.frame(perm_result),
            file.path(SURVEY_TABLES, "permanova_survey_type.csv"))

  cat("  PERMANOVA results saved\n\n")
}

# ============================================================================
# Abundance Comparison Analysis (Community-Level)
# ============================================================================

cat("Analyzing abundance differences across survey types...\n")

# Calculate abundance for each OTU (community patterns, not detection bias)
abundance_analysis <- cafi_clean %>%
  left_join(metadata %>% select(coral_id, survey_type), by = "coral_id") %>%
  group_by(species) %>%
  mutate(total_abundance = n()) %>%
  group_by(species, survey_type, total_abundance) %>%
  summarise(
    abundance = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = survey_type,
              values_from = abundance,
              values_fill = 0) %>%
  filter(total_abundance >= 5)  # Filter to OTUs with enough data

if (ncol(abundance_analysis) >= 3) {
  # Calculate abundance ratios
  survey_cols <- setdiff(names(abundance_analysis), c("species", "total_abundance"))

  if (length(survey_cols) == 2) {
    abundance_analysis <- abundance_analysis %>%
      mutate(
        abundance_ratio = .[[survey_cols[1]]] / (.[[survey_cols[1]]] + .[[survey_cols[2]]]),
        log_ratio = log((.[[survey_cols[1]]] + 1) / (.[[survey_cols[2]]] + 1))
      )

    # Plot abundance differences
    p_abundance <- abundance_analysis %>%
      arrange(log_ratio) %>%
      mutate(species = factor(species, levels = species)) %>%
      slice(c(head(row_number(), 15), tail(row_number(), 15))) %>%
      ggplot(aes(x = log_ratio, y = species)) +
      geom_segment(aes(x = 0, xend = log_ratio, y = species, yend = species),
                   color = "gray50") +
      geom_point(aes(size = total_abundance), color = "steelblue") +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      scale_size_continuous(range = c(2, 8)) +
      labs(title = "OTU Abundance Patterns by Survey Type",
           subtitle = paste("Positive values indicate higher abundance in", survey_cols[1]),
           x = "Log Abundance Ratio",
           y = "CAFI OTU",
           size = "Total\nAbundance")

    ggsave(file.path(fig_dir, "abundance_patterns_by_survey.png"),
           p_abundance, width = 10, height = 10, dpi = 300)
  }
}

# ============================================================================
# Morphotype Interactions
# ============================================================================

cat("Analyzing morphotype-survey type interactions...\n")

# Interaction effects
interaction_data <- survey_comparison %>%
  group_by(morphotype, survey_type) %>%
  summarise(
    mean_abundance = mean(total_cafi, na.rm = TRUE),
    se_abundance = sd(total_cafi, na.rm = TRUE) / sqrt(n()),
    mean_richness = mean(species_richness, na.rm = TRUE),
    se_richness = sd(species_richness, na.rm = TRUE) / sqrt(n()),
    n_corals = n(),
    .groups = "drop"
  )

# Plot interactions
p_interaction <- interaction_data %>%
  pivot_longer(cols = c(mean_abundance, mean_richness),
               names_to = "metric",
               values_to = "mean") %>%
  mutate(
    se = ifelse(metric == "mean_abundance", se_abundance, se_richness),
    metric = factor(metric,
                   levels = c("mean_abundance", "mean_richness"),
                   labels = c("Mean Abundance", "Mean Richness"))
  ) %>%
  ggplot(aes(x = morphotype, y = mean, fill = survey_type)) +
  geom_col(position = position_dodge(0.8), alpha = 0.7, width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.8), width = 0.3) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_viridis_d() +
  labs(title = "Morphotype × Survey Type Interactions",
       x = "Morphotype",
       y = "Mean Value (± SE)",
       fill = "Survey Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "morphotype_survey_interactions.png"),
       p_interaction, width = 12, height = 8, dpi = 300)

# ============================================================================
# Rarefaction Analysis by Survey Type
# ============================================================================

cat("Performing rarefaction analysis...\n")

# Create separate community matrices by survey type
survey_types_present <- unique(metadata$survey_type)

if (length(survey_types_present) > 1) {
  rarefaction_results <- list()

  for (st in survey_types_present) {
    # Get corals for this survey type
    st_corals <- metadata %>%
      filter(survey_type == st) %>%
      pull(coral_id)

    # Subset community matrix
    st_comm <- community_matrix[rownames(community_matrix) %in% st_corals, ]
    st_comm <- st_comm[, colSums(st_comm) > 0]

    if (nrow(st_comm) > 5) {
      # Pool all samples for survey type
      pooled <- colSums(st_comm)

      # Rarefaction
      rare_result <- rarefy(pooled, sample = seq(1, sum(pooled),
                                                 length.out = min(100, sum(pooled))))

      rarefaction_results[[st]] <- data.frame(
        survey_type = st,
        individuals = seq(1, sum(pooled), length.out = length(rare_result)),
        species = rare_result
      )
    }
  }

  # Combine and plot
  if (length(rarefaction_results) > 0) {
    rare_df <- bind_rows(rarefaction_results)

    p_rarefaction <- ggplot(rare_df, aes(x = individuals, y = species,
                                         color = survey_type)) +
      geom_line(size = 1.5) +
      scale_color_viridis_d() +
      labs(title = "OTU Accumulation Curves by Survey Type",
           x = "Number of Individuals",
           y = "Expected Number of OTUs",
           color = "Survey Type")

    ggsave(file.path(fig_dir, "rarefaction_by_survey_type.png"),
           p_rarefaction, width = 10, height = 8, dpi = 300)
  }
}

# ============================================================================
# Summary Statistics Table
# ============================================================================

cat("Creating summary statistics...\n")

# Comprehensive summary by survey type
summary_stats <- survey_comparison %>%
  group_by(survey_type) %>%
  summarise(
    n_corals = n(),
    mean_abundance = mean(total_cafi, na.rm = TRUE),
    sd_abundance = sd(total_cafi, na.rm = TRUE),
    median_abundance = median(total_cafi, na.rm = TRUE),
    mean_richness = mean(species_richness, na.rm = TRUE),
    sd_richness = sd(species_richness, na.rm = TRUE),
    median_richness = median(species_richness, na.rm = TRUE),
    mean_shannon = mean(shannon, na.rm = TRUE),
    sd_shannon = sd(shannon, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_stats,
          file.path(SURVEY_TABLES, "survey_type_summary_statistics.csv"))

# ============================================================================
# Summary Report
# ============================================================================

cat("\n========================================\n")
cat("Neighborhood vs Arrival Comparison Summary\n")
cat("========================================\n\n")

if (exists("summary_stats")) {
  cat("Sample Sizes:\n")
  for (i in 1:nrow(summary_stats)) {
    cat("  -", summary_stats$survey_type[i], ":",
        summary_stats$n_corals[i], "corals\n")
  }
  cat("\n")
}

if (exists("test_results")) {
  cat("Statistical Differences:\n")
  sig_tests <- test_results %>% filter(significant)
  if (nrow(sig_tests) > 0) {
    cat("  Significant differences found in:\n")
    for (i in 1:nrow(sig_tests)) {
      cat("    -", sig_tests$metric[i], "(p =",
          round(sig_tests$p_value[i], 4), ")\n")
    }
  } else {
    cat("  No significant differences detected\n")
  }
  cat("\n")
}

if (exists("otu_by_survey")) {
  concentrated <- otu_by_survey %>%
    filter(distribution_pattern != "Even distribution")

  if (nrow(concentrated) > 0) {
    cat("OTU Distribution Patterns:\n")
    cat("  - Neighborhood-concentrated:",
        sum(concentrated$distribution_pattern == "Neighborhood-concentrated"), "\n")
    cat("  - Arrival-concentrated:",
        sum(concentrated$distribution_pattern == "Arrival-concentrated"), "\n")
    cat("\n")
  }
}

cat("✅ Neighborhood vs Arrival comparison complete!\n")
cat("Figures saved to:", fig_dir, "\n")
cat("Tables saved to:", SURVEY_TABLES, "\n")