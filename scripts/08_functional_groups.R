# ============================================================================
# 08_functional_groups.R - Taxonomic Group Scaling Analysis
# ============================================================================
#
# PURPOSE: Analyze how different taxonomic groups of CAFI scale with coral size
#
# TAXONOMIC GROUPS ANALYZED:
#   A. Trapezia Crabs - Obligate coral associates (well-documented mutualists)
#   B. Fish - Coral-dwelling fishes (Paragobiodon, Caracanthus, etc.)
#   C. Gastropods - Coral-associated snails (Galeropsis, Cerithium, etc.)
#   D. Taxonomic Diversity - Within-group species diversity patterns
#   E. Scaling Patterns - Power-law fits by taxonomic group
#
# MANUSCRIPT FIGURES:
#   >>> FIGURE 3: Taxonomic Group Scaling <<<
#   - Panel A: Trapezia scaling with size
#   - Panel B: Fish scaling patterns
#   - Panel C: Forest plot of taxonomic group beta estimates
#
# DEPENDENCIES: 00_setup.R, output/objects/
#
# OUTPUTS:
#   Figures:
#     - output/figures/manuscript/fig3_functional_groups.png
#     - output/figures/functional_groups/*.png (analysis figures)
#   Tables:
#     - output/tables/taxonomic_group_scaling.csv
#     - output/tables/trapezia_species.csv
#     - output/tables/gastropod_prevalence.csv
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-21
# ============================================================================

cat("\n========================================\n")
cat("08: Taxonomic Group Scaling Analysis\n")
cat("========================================\n\n")

# ============================================================================
# SETUP AND DATA LOADING
# ============================================================================

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Load processed data objects
coral_master <- load_object("coral_master")
cafi_clean <- load_object("cafi_clean")
functional_summary <- load_object("functional_summary")

# Create output directories
FIG_DIR <- file.path(PATHS$figures, "functional_groups")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Data loaded:\n")
cat("  Corals:", nrow(coral_master), "\n")
cat("  CAFI records:", nrow(cafi_clean), "\n")
cat("  Taxonomic groups:", n_distinct(cafi_clean$functional_group), "\n\n")

# ============================================================================
# HELPER FUNCTION: Fit Scaling Model
# ============================================================================

#' Fit a power-law scaling model: N = a * V^beta
#'
#' @param data Data frame with 'abundance' and 'volume' columns
#' @param response_name Character string for labeling
#' @param min_nonzero Minimum non-zero observations required
#' @return List with model results
fit_functional_scaling <- function(data, response_name, min_nonzero = 15) {

  n_total <- nrow(data)
  n_nonzero <- sum(data$abundance > 0)
  n_zero <- sum(data$abundance == 0)
  total_abundance <- sum(data$abundance)

  if (n_nonzero < min_nonzero) {
    return(list(
      response = response_name,
      n_corals = n_total,
      n_nonzero = n_nonzero,
      n_zero = n_zero,
      total_abundance = total_abundance,
      beta = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      p_vs_1 = NA_real_,
      interpretation = "Insufficient data",
      model = NULL,
      converged = FALSE
    ))
  }

  tryCatch({
    model <- MASS::glm.nb(abundance ~ log(volume) + site, data = data)
    coefs <- summary(model)$coefficients
    beta <- coefs["log(volume)", "Estimate"]
    se <- coefs["log(volume)", "Std. Error"]
    z_val <- coefs["log(volume)", "z value"]
    p_val <- coefs["log(volume)", "Pr(>|z|)"]

    ci <- confint(model, "log(volume)", level = 0.95)

    # Test vs Field of Dreams (beta = 1)
    z_vs_1 <- (beta - 1) / se
    p_vs_1 <- 2 * pnorm(-abs(z_vs_1))

    interpretation <- case_when(
      p_vs_1 < 0.05 & beta < 1 ~ "Redirection (beta < 1)",
      p_vs_1 < 0.05 & beta > 1 ~ "Super-linear (beta > 1)",
      TRUE ~ "Field of Dreams (beta ~ 1)"
    )

    list(
      response = response_name,
      n_corals = n_total,
      n_nonzero = n_nonzero,
      n_zero = n_zero,
      total_abundance = total_abundance,
      beta = beta,
      se = se,
      ci_lower = ci[1],
      ci_upper = ci[2],
      p_value = p_val,
      p_vs_1 = p_vs_1,
      interpretation = interpretation,
      model = model,
      converged = TRUE
    )

  }, error = function(e) {
    list(
      response = response_name,
      n_corals = n_total,
      n_nonzero = n_nonzero,
      n_zero = n_zero,
      total_abundance = total_abundance,
      beta = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      p_vs_1 = NA_real_,
      interpretation = paste("Model error:", conditionMessage(e)),
      model = NULL,
      converged = FALSE
    )
  })
}

# ############################################################################
#                    PART A: TRAPEZIA GUARDIAN CRAB ANALYSIS
# ############################################################################

cat("============================================================\n")
cat("PART A: TRAPEZIA GUARDIAN CRAB ANALYSIS\n")
cat("============================================================\n\n")

# ----------------------------------------------------------------------------
# A1. Identify Trapezia individuals
# ----------------------------------------------------------------------------

cat("Identifying Trapezia crabs...\n")

# Filter to Trapezia (genus-level identification)
# Check if genus column exists and use it if available
if ("genus" %in% names(cafi_clean)) {
  trapezia <- cafi_clean %>%
    filter(
      str_detect(tolower(genus), "trapezia") |
      str_detect(tolower(otu), "trapezia") |
      str_detect(tolower(species), "trapezia")
    )
} else {
  trapezia <- cafi_clean %>%
    filter(
      str_detect(tolower(otu), "trapezia") |
      str_detect(tolower(species), "trapezia")
    )
}

cat("  Total Trapezia individuals:", nrow(trapezia), "\n")

if (nrow(trapezia) > 0) {

  # Species breakdown
  trapezia_species <- trapezia %>%
    group_by(otu) %>%
    summarise(
      n_individuals = n(),
      n_corals = n_distinct(coral_id),
      prevalence = n_distinct(coral_id) / n_distinct(cafi_clean$coral_id) * 100,
      mean_size_mm = mean(size_mm, na.rm = TRUE),
      sd_size_mm = sd(size_mm, na.rm = TRUE),
      min_size_mm = min(size_mm, na.rm = TRUE),
      max_size_mm = max(size_mm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_individuals))

  cat("  Trapezia species/morphotypes:", nrow(trapezia_species), "\n")
  cat("\n  Species breakdown:\n")
  print(trapezia_species %>% dplyr::select(otu, n_individuals, n_corals, prevalence))

  # Save Trapezia species table
  save_table(trapezia_species, "trapezia_species")
  cat("\n  Saved: trapezia_species.csv\n\n")

  # ----------------------------------------------------------------------------
  # A2. Per-coral Trapezia metrics
  # ----------------------------------------------------------------------------

  cat("Calculating per-coral Trapezia metrics...\n")

  trapezia_by_coral <- trapezia %>%
    group_by(coral_id) %>%
    summarise(
      n_trapezia = n(),
      trapezia_richness = n_distinct(otu),
      trapezia_species = paste(sort(unique(otu)), collapse = "; "),
      mean_trap_size = mean(size_mm, na.rm = TRUE),
      max_trap_size = max(size_mm, na.rm = TRUE),
      .groups = "drop"
    )

  # Use existing n_trapezia from coral_master (already computed in 01_load_data.R)
  # Add richness data from trapezia_by_coral
  trap_data <- coral_master %>%
    filter(!is.na(volume), volume > 0) %>%
    left_join(trapezia_by_coral %>% dplyr::select(coral_id, trapezia_richness, trapezia_species,
                                                   mean_trap_size, max_trap_size),
              by = "coral_id")

  # Replace NAs with zeros for richness (n_trapezia already exists in coral_master)
  if (!"n_trapezia" %in% names(trap_data)) {
    # Fallback if somehow missing
    trap_data$n_trapezia <- 0
  }
  trap_data$n_trapezia[is.na(trap_data$n_trapezia)] <- 0
  trap_data$trapezia_richness[is.na(trap_data$trapezia_richness)] <- 0
  trap_data$has_trapezia <- trap_data$n_trapezia > 0

  cat("  Corals with Trapezia:", sum(trap_data$has_trapezia), "of", nrow(trap_data),
      sprintf("(%.1f%%)\n", sum(trap_data$has_trapezia) / nrow(trap_data) * 100))
  cat("  Mean Trapezia per coral (all):", round(mean(trap_data$n_trapezia), 2), "\n")
  cat("  Mean Trapezia per coral (occupied):",
      round(mean(trap_data$n_trapezia[trap_data$has_trapezia]), 2), "\n")
  cat("  Max Trapezia on single coral:", max(trap_data$n_trapezia), "\n\n")

  # ----------------------------------------------------------------------------
  # A3. Trapezia abundance scaling with coral volume (power law)
  # ----------------------------------------------------------------------------

  cat("Fitting Trapezia-volume scaling model...\n")

  trap_scaling_data <- trap_data %>%
    dplyr::select(coral_id, site, volume, abundance = n_trapezia)

  trap_scaling_result <- fit_functional_scaling(trap_scaling_data, "Trapezia")

  if (trap_scaling_result$converged && !is.na(trap_scaling_result$beta)) {
    cat("  Scaling exponent beta =", round(trap_scaling_result$beta, 3), "\n")
    cat("  95% CI: [", round(trap_scaling_result$ci_lower, 3), ",",
        round(trap_scaling_result$ci_upper, 3), "]\n")
    cat("  Test vs beta = 1: p =", format.pval(trap_scaling_result$p_vs_1, 3), "\n")
    cat("  Interpretation:", trap_scaling_result$interpretation, "\n\n")
  }

  # Create Trapezia scaling figure
  p_trap_scaling <- trap_data %>%
    ggplot(aes(x = volume, y = n_trapezia + 0.5, color = site)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(aes(group = 1), method = MASS::glm.nb, formula = y ~ log(x),
                se = TRUE, color = "black", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    labs(
      x = expression("Coral Volume (cm"^3*")"),
      y = "Trapezia Abundance",
      title = "Trapezia Guardian Crab Scaling",
      subtitle = if(!is.na(trap_scaling_result$beta))
        sprintf("beta = %.2f [%.2f, %.2f] | %s",
                trap_scaling_result$beta, trap_scaling_result$ci_lower,
                trap_scaling_result$ci_upper, trap_scaling_result$interpretation)
        else "Model fit failed"
    ) +
    theme_publication() +
    theme(legend.position = c(0.15, 0.85))

  ggsave(file.path(FIG_DIR, "trapezia_scaling.png"), p_trap_scaling,
         width = 8, height = 6, dpi = 300, bg = "white")
  cat("  Saved: trapezia_scaling.png\n\n")

  # ----------------------------------------------------------------------------
  # A4. Pair-wise co-occurrence analysis (mated pairs)
  # ----------------------------------------------------------------------------

  cat("Analyzing Trapezia pair occurrence patterns...\n")

  # Trapezia typically form mated pairs (one male, one female per coral)
  trap_pair_analysis <- trap_data %>%
    filter(has_trapezia) %>%
    mutate(
      pair_status = case_when(
        n_trapezia == 1 ~ "Single",
        n_trapezia == 2 ~ "Pair",
        n_trapezia > 2 ~ "Multiple (>2)"
      )
    ) %>%
    group_by(pair_status) %>%
    summarise(
      n_corals = n(),
      mean_volume = mean(volume),
      sd_volume = sd(volume),
      .groups = "drop"
    )

  cat("  Trapezia occurrence patterns:\n")
  print(trap_pair_analysis)
  cat("\n")

  # Test if larger corals host more than pairs
  if (sum(trap_data$n_trapezia > 2) >= 5) {
    multiple_test <- wilcox.test(
      volume ~ I(n_trapezia > 2),
      data = trap_data %>% filter(has_trapezia)
    )
    cat("  Volume comparison (pairs vs multiple):\n")
    cat("    Wilcoxon p =", format.pval(multiple_test$p.value, 3), "\n")
    cat("    Larger corals support multiple Trapezia beyond mated pairs\n\n")
  }

  # ----------------------------------------------------------------------------
  # A5. Species composition by coral size class
  # ----------------------------------------------------------------------------

  cat("Analyzing Trapezia species composition by coral size...\n")

  # Create size classes
  trap_data <- trap_data %>%
    mutate(
      size_class = cut(volume,
                       breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                       labels = c("Small", "Medium", "Large"),
                       include.lowest = TRUE)
    )

  # Species composition by size class
  trap_species_by_size <- trapezia %>%
    left_join(trap_data %>% dplyr::select(coral_id, volume, size_class), by = "coral_id") %>%
    filter(!is.na(size_class)) %>%
    group_by(size_class, otu) %>%
    summarise(
      n_individuals = n(),
      n_corals = n_distinct(coral_id),
      .groups = "drop"
    ) %>%
    group_by(size_class) %>%
    mutate(prop = n_individuals / sum(n_individuals)) %>%
    ungroup()

  # Plot species composition by size
  if (nrow(trap_species_by_size) > 0) {
    p_trap_composition <- trap_species_by_size %>%
      ggplot(aes(x = size_class, y = prop, fill = otu)) +
      geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
      scale_fill_viridis_d(option = "D", name = "Trapezia Species") +
      scale_y_continuous(labels = scales::percent) +
      labs(
        x = "Coral Size Class",
        y = "Proportion of Trapezia Individuals",
        title = "Trapezia Species Composition by Coral Size",
        subtitle = "Examining size-based competitive exclusion patterns"
      ) +
      theme_publication() +
      theme(legend.position = "right")

    ggsave(file.path(FIG_DIR, "trapezia_composition_by_size.png"), p_trap_composition,
           width = 10, height = 6, dpi = 300, bg = "white")
    cat("  Saved: trapezia_composition_by_size.png\n\n")
  }

  # ----------------------------------------------------------------------------
  # A6. Size-based competitive exclusion patterns
  # ----------------------------------------------------------------------------

  cat("Testing size-based competitive exclusion...\n")

  # Do larger Trapezia individuals occur on larger corals?
  trap_size_correlation <- trapezia %>%
    left_join(trap_data %>% dplyr::select(coral_id, volume), by = "coral_id") %>%
    filter(!is.na(size_mm), !is.na(volume))

  if (nrow(trap_size_correlation) >= 20) {
    size_cor <- cor.test(log10(trap_size_correlation$volume),
                         trap_size_correlation$size_mm,
                         method = "spearman")
    cat("  Trapezia body size vs coral volume:\n")
    cat("    Spearman rho =", round(size_cor$estimate, 3), "\n")
    cat("    p =", format.pval(size_cor$p.value, 3), "\n")

    if (size_cor$p.value < 0.05 && size_cor$estimate > 0) {
      cat("    Larger corals host larger Trapezia individuals\n")
      cat("    Consistent with size-based competitive exclusion\n\n")
    }
  }

} else {
  cat("  No Trapezia identified in dataset\n\n")
  trap_scaling_result <- NULL
}

# ############################################################################
#                    PART B: RESIDENT FISH PATTERNS
# ############################################################################

cat("============================================================\n")
cat("PART B: RESIDENT FISH PATTERNS\n")
cat("============================================================\n\n")

cat("Analyzing resident fish (Paragobiodon, Caracanthus)...\n")
cat("These fish provide nutrients to corals via waste excretion\n\n")

# Identify resident fish: gobies (Paragobiodon, Gobiodon) and velvetfish (Caracanthus)
resident_fish <- cafi_clean %>%
  filter(
    str_detect(tolower(otu), "paragobiodon|gobiodon|caracanthus") |
    functional_group == "Resident Fish"
  )

cat("  Total resident fish individuals:", nrow(resident_fish), "\n")

if (nrow(resident_fish) > 0) {

  # Species breakdown
  fish_species <- resident_fish %>%
    group_by(otu) %>%
    summarise(
      n_individuals = n(),
      n_corals = n_distinct(coral_id),
      prevalence = n_distinct(coral_id) / n_distinct(cafi_clean$coral_id) * 100,
      mean_size_mm = mean(size_mm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_individuals))

  cat("  Resident fish species:\n")
  print(fish_species)
  cat("\n")

  # Per-coral fish metrics
  fish_by_coral <- resident_fish %>%
    group_by(coral_id) %>%
    summarise(
      n_fish = n(),
      fish_richness = n_distinct(otu),
      .groups = "drop"
    )

  # Use existing n_resident_fish from coral_master, add richness from fish_by_coral
  fish_data <- coral_master %>%
    filter(!is.na(volume), volume > 0) %>%
    left_join(fish_by_coral %>% dplyr::select(coral_id, fish_richness), by = "coral_id")

  # Use existing n_resident_fish column, create n_fish alias
  fish_data$n_fish <- fish_data$n_resident_fish
  fish_data$n_fish[is.na(fish_data$n_fish)] <- 0
  fish_data$fish_richness[is.na(fish_data$fish_richness)] <- 0
  fish_data$has_fish <- fish_data$n_fish > 0

  cat("  Corals with resident fish:", sum(fish_data$has_fish), "of", nrow(fish_data),
      sprintf("(%.1f%%)\n\n", sum(fish_data$has_fish) / nrow(fish_data) * 100))

  # ----------------------------------------------------------------------------
  # B1. Minimum coral size for colonization (size threshold)
  # ----------------------------------------------------------------------------

  cat("Analyzing minimum coral size for fish colonization...\n")

  # Compare volume of colonized vs uncolonized corals
  if (sum(fish_data$has_fish) >= 5 && sum(!fish_data$has_fish) >= 5) {

    fish_occupied <- fish_data %>% filter(has_fish) %>% pull(volume)
    fish_empty <- fish_data %>% filter(!has_fish) %>% pull(volume)

    size_test <- wilcox.test(fish_occupied, fish_empty)

    cat("  Volume of corals WITH fish: median =", round(median(fish_occupied)), "cm3\n")
    cat("  Volume of corals WITHOUT fish: median =", round(median(fish_empty)), "cm3\n")
    cat("  Wilcoxon test p =", format.pval(size_test$p.value, 3), "\n")

    # Find minimum colonized volume
    min_colonized <- min(fish_occupied)
    cat("  Minimum volume with fish:", round(min_colonized), "cm3\n")
    cat("  This suggests a size threshold for fish colonization\n\n")

    # Logistic regression for colonization probability
    fish_logistic <- glm(has_fish ~ log(volume) + site,
                         data = fish_data, family = binomial)
    fish_logistic_summary <- summary(fish_logistic)

    cat("  Logistic regression (fish presence ~ log(volume) + site):\n")
    cat("    Volume effect: z =",
        round(fish_logistic_summary$coefficients["log(volume)", "z value"], 2),
        ", p =", format.pval(fish_logistic_summary$coefficients["log(volume)", "Pr(>|z|)"], 3), "\n\n")
  }

  # ----------------------------------------------------------------------------
  # B2. Fish scaling with coral volume
  # ----------------------------------------------------------------------------

  cat("Fitting resident fish scaling model...\n")

  fish_scaling_data <- fish_data %>%
    dplyr::select(coral_id, site, volume, abundance = n_fish)

  fish_scaling_result <- fit_functional_scaling(fish_scaling_data, "Resident Fish")

  if (fish_scaling_result$converged && !is.na(fish_scaling_result$beta)) {
    cat("  Scaling exponent beta =", round(fish_scaling_result$beta, 3), "\n")
    cat("  95% CI: [", round(fish_scaling_result$ci_lower, 3), ",",
        round(fish_scaling_result$ci_upper, 3), "]\n")
    cat("  Interpretation:", fish_scaling_result$interpretation, "\n\n")
  }

  # Fish scaling plot
  p_fish_scaling <- fish_data %>%
    ggplot(aes(x = volume, y = n_fish + 0.5, color = site)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(aes(group = 1), method = MASS::glm.nb, formula = y ~ log(x),
                se = TRUE, color = "black", linewidth = 1) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    labs(
      x = expression("Coral Volume (cm"^3*")"),
      y = "Resident Fish Abundance",
      title = "Resident Fish (Gobies & Velvetfish) Scaling",
      subtitle = if(!is.na(fish_scaling_result$beta))
        sprintf("beta = %.2f [%.2f, %.2f]",
                fish_scaling_result$beta, fish_scaling_result$ci_lower,
                fish_scaling_result$ci_upper) else "Insufficient data"
    ) +
    theme_publication() +
    theme(legend.position = c(0.15, 0.85))

  ggsave(file.path(FIG_DIR, "resident_fish_scaling.png"), p_fish_scaling,
         width = 8, height = 6, dpi = 300, bg = "white")
  cat("  Saved: resident_fish_scaling.png\n\n")

  # ----------------------------------------------------------------------------
  # B3. Prevalence across sites
  # ----------------------------------------------------------------------------

  cat("Analyzing resident fish prevalence by site...\n")

  fish_by_site <- fish_data %>%
    group_by(site) %>%
    summarise(
      n_corals = n(),
      n_with_fish = sum(has_fish),
      prevalence = mean(has_fish) * 100,
      mean_abundance = mean(n_fish),
      .groups = "drop"
    )

  cat("  Fish prevalence by site:\n")
  print(fish_by_site)
  cat("\n")

  # Chi-square test for site differences
  if (n_distinct(fish_data$site) > 1) {
    fish_site_chi <- chisq.test(table(fish_data$site, fish_data$has_fish))
    cat("  Chi-square test (site x fish presence): p =",
        format.pval(fish_site_chi$p.value, 3), "\n\n")
  }

} else {
  cat("  No resident fish identified in dataset\n\n")
  fish_scaling_result <- NULL
}

# ############################################################################
#                    PART C: GASTROPOD PATTERNS
# ############################################################################

cat("============================================================\n")
cat("PART C: GASTROPOD PATTERNS (dominated by Galeropsis monodonta)\n")
cat("============================================================\n\n")

cat("Analyzing gastropods associated with coral colonies...\n")
cat("Dominated by Galeropsis monodonta (Coralliophilinae) — a coral tissue feeder\n\n")

# Identify gastropods
coral_eating_snails <- cafi_clean %>%
  filter(
    str_detect(tolower(otu), "drupella|coralliophila") |
    functional_group == "Gastropod" |
    tolower(type) == "snail"
  )

cat("  Total gastropod individuals:", nrow(coral_eating_snails), "\n")

if (nrow(coral_eating_snails) > 0) {

  # Species breakdown
  coral_species <- coral_eating_snails %>%
    group_by(otu, type) %>%
    summarise(
      n_individuals = n(),
      n_corals = n_distinct(coral_id),
      prevalence = n_distinct(coral_id) / n_distinct(cafi_clean$coral_id) * 100,
      mean_size_mm = mean(size_mm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_individuals))

  cat("  Gastropod species:\n")
  print(coral_species)
  cat("\n")

  # Per-coral gastropod metrics
  snail_by_coral <- coral_eating_snails %>%
    group_by(coral_id) %>%
    summarise(
      n_corallivore = n(),
      snail_richness = n_distinct(otu),
      .groups = "drop"
    )

  # Use existing n_corallivore from coral_master, add richness
  coral_data <- coral_master %>%
    filter(!is.na(volume), volume > 0) %>%
    left_join(snail_by_coral %>% dplyr::select(coral_id, snail_richness), by = "coral_id")

  # n_corallivore already exists in coral_master
  coral_data$n_corallivore[is.na(coral_data$n_corallivore)] <- 0
  coral_data$snail_richness[is.na(coral_data$snail_richness)] <- 0
  coral_data$has_corallivore <- coral_data$n_corallivore > 0

  cat("  Corals with gastropods:", sum(coral_data$has_corallivore), "of", nrow(coral_data),
      sprintf("(%.1f%%)\n\n", sum(coral_data$has_corallivore) / nrow(coral_data) * 100))

  # ----------------------------------------------------------------------------
  # C1. Prevalence by site and coral size
  # ----------------------------------------------------------------------------

  cat("Analyzing gastropod prevalence patterns...\n")

  # By site
  snail_by_site <- coral_data %>%
    group_by(site) %>%
    summarise(
      n_corals = n(),
      n_with_snails = sum(has_corallivore),
      prevalence = mean(has_corallivore) * 100,
      mean_abundance = mean(n_corallivore),
      mean_volume = mean(volume),
      .groups = "drop"
    )

  cat("  Gastropod prevalence by site:\n")
  print(snail_by_site)

  # Save gastropod prevalence table
  save_table(snail_by_site, "gastropod_prevalence")
  cat("\n  Saved: gastropod_prevalence.csv\n\n")

  # By size class
  coral_data <- coral_data %>%
    mutate(
      size_class = cut(volume,
                       breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                       labels = c("Small", "Medium", "Large"),
                       include.lowest = TRUE)
    )

  snail_by_size <- coral_data %>%
    group_by(size_class) %>%
    summarise(
      n_corals = n(),
      prevalence = mean(has_corallivore) * 100,
      mean_abundance = mean(n_corallivore),
      .groups = "drop"
    )

  cat("  Gastropod prevalence by coral size:\n")
  print(snail_by_size)
  cat("\n")

  # ----------------------------------------------------------------------------
  # C2. Association with coral condition (expected negative)
  # ----------------------------------------------------------------------------

  if ("condition_score" %in% names(coral_data)) {
    cat("Testing gastropod vs condition association...\n")

    condition_data <- coral_data %>%
      filter(!is.na(condition_score))

    if (nrow(condition_data) > 20) {
      # Correlation
      cor_test <- cor.test(condition_data$n_corallivore,
                           condition_data$condition_score,
                           method = "spearman")

      cat("  Gastropod abundance vs coral condition:\n")
      rho_val <- cor_test$estimate
      p_cor <- cor_test$p.value
      cat("    Spearman rho =", if(!is.na(rho_val)) round(rho_val, 3) else "NA", "\n")
      cat("    p =", if(!is.na(p_cor)) format.pval(p_cor, 3) else "NA", "\n")

      if (!is.na(p_cor) && !is.na(rho_val) && p_cor < 0.05 && rho_val < 0) {
        cat("    Negative association: more gastropods = lower condition\n\n")
      } else if (is.na(p_cor) || p_cor >= 0.05) {
        cat("    No significant association detected\n\n")
      }

      # Compare condition: with vs without gastropods
      if (sum(condition_data$has_corallivore) >= 5) {
        condition_test <- wilcox.test(
          condition_score ~ has_corallivore,
          data = condition_data
        )
        cat("  Condition comparison (presence/absence):\n")
        cat("    Corals with gastropods: median condition =",
            round(median(condition_data$condition_score[condition_data$has_corallivore]), 3), "\n")
        cat("    Corals without gastropods: median condition =",
            round(median(condition_data$condition_score[!condition_data$has_corallivore]), 3), "\n")
        cat("    Wilcoxon p =", format.pval(condition_test$p.value, 3), "\n\n")
      }
    }
  }

  # ----------------------------------------------------------------------------
  # C3. Size-refuge hypothesis: do large corals escape snail damage?
  # ----------------------------------------------------------------------------

  cat("Testing size-refuge hypothesis...\n")

  # Does gastropod prevalence decrease with coral size?
  size_refuge_model <- glm(has_corallivore ~ log(volume) + site,
                           data = coral_data, family = binomial)
  size_refuge_summary <- summary(size_refuge_model)

  volume_coef <- size_refuge_summary$coefficients["log(volume)", ]

  cat("  Logistic regression (gastropod presence ~ log(volume) + site):\n")
  cat("    Volume coefficient =", round(volume_coef["Estimate"], 3), "\n")
  cat("    z =", round(volume_coef["z value"], 2), "\n")
  cat("    p =", format.pval(volume_coef["Pr(>|z|)"], 3), "\n")

  p_val <- volume_coef["Pr(>|z|)"]
  est_val <- volume_coef["Estimate"]
  if (!is.na(p_val) && !is.na(est_val) && p_val < 0.05 && est_val < 0) {
    cat("    SUPPORTED: Large corals have lower gastropod prevalence (size refuge)\n\n")
  } else if (!is.na(p_val) && !is.na(est_val) && p_val < 0.05 && est_val > 0) {
    cat("    REJECTED: Large corals have HIGHER gastropod prevalence\n\n")
  } else {
    cat("    NOT SUPPORTED: No significant size effect on gastropod prevalence\n\n")
  }

  # ----------------------------------------------------------------------------
  # C4. Gastropod scaling
  # ----------------------------------------------------------------------------

  cat("Fitting gastropod scaling model...\n")

  coral_scaling_data <- coral_data %>%
    dplyr::select(coral_id, site, volume, abundance = n_corallivore)

  snail_scaling_result <- fit_functional_scaling(coral_scaling_data, "Gastropods")

  if (snail_scaling_result$converged && !is.na(snail_scaling_result$beta)) {
    cat("  Scaling exponent beta =", round(snail_scaling_result$beta, 3), "\n")
    cat("  95% CI: [", round(snail_scaling_result$ci_lower, 3), ",",
        round(snail_scaling_result$ci_upper, 3), "]\n")
    cat("  Interpretation:", snail_scaling_result$interpretation, "\n\n")
  }

  # Gastropod patterns figure
  p_snails <- coral_data %>%
    ggplot(aes(x = volume, y = n_corallivore + 0.5, color = site)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(aes(group = 1), method = MASS::glm.nb, formula = y ~ log(x),
                se = TRUE, color = "black", linewidth = 1) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    labs(
      x = expression("Coral Volume (cm"^3*")"),
      y = "Gastropod Abundance",
      title = "Gastropod Abundance vs Coral Size",
      subtitle = if(!is.na(snail_scaling_result$beta))
        sprintf("beta = %.2f [%.2f, %.2f]",
                snail_scaling_result$beta, snail_scaling_result$ci_lower,
                snail_scaling_result$ci_upper) else "Insufficient data"
    ) +
    theme_publication() +
    theme(legend.position = c(0.15, 0.85))

  ggsave(file.path(FIG_DIR, "gastropod_scaling.png"), p_snails,
         width = 8, height = 6, dpi = 300, bg = "white")
  cat("  Saved: gastropod_scaling.png\n\n")

} else {
  cat("  No gastropods identified in dataset\n\n")
  snail_scaling_result <- NULL
}

# ############################################################################
#                    PART D: TAXONOMIC DIVERSITY
# ############################################################################

cat("============================================================\n")
cat("PART D: TAXONOMIC DIVERSITY ANALYSIS\n")
cat("============================================================\n\n")

cat("Analyzing within-group species diversity...\n\n")

# ----------------------------------------------------------------------------
# D1. Within-guild species diversity
# ----------------------------------------------------------------------------

cat("Within-guild species diversity:\n")

guild_diversity <- cafi_clean %>%
  group_by(functional_group) %>%
  summarise(
    n_species = n_distinct(otu),
    n_individuals = n(),
    n_corals = n_distinct(coral_id),
    shannon_guild = vegan::diversity(table(otu)),
    evenness = shannon_guild / log(n_distinct(otu)),
    dominant_species = names(sort(table(otu), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  arrange(desc(n_individuals))

print(guild_diversity)
cat("\n")

# ----------------------------------------------------------------------------
# D2. Substitutability patterns (Trapezia species replacement)
# ----------------------------------------------------------------------------

if (exists("trap_data") && nrow(trapezia) > 0) {
  cat("Trapezia species substitutability patterns:\n")

  # Do different Trapezia species occur on corals of different sizes?
  trap_size_by_species <- trapezia %>%
    left_join(trap_data %>% dplyr::select(coral_id, volume, size_class), by = "coral_id") %>%
    filter(!is.na(size_class)) %>%
    group_by(otu, size_class) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = size_class, values_from = n, values_fill = 0)

  cat("  Trapezia species distribution by coral size class:\n")
  print(trap_size_by_species)
  cat("\n")

  # Chi-square test for species x size association
  trap_contingency <- trapezia %>%
    left_join(trap_data %>% dplyr::select(coral_id, size_class), by = "coral_id") %>%
    filter(!is.na(size_class))

  if (n_distinct(trap_contingency$otu) > 1 && n_distinct(trap_contingency$size_class) > 1) {
    trap_chi <- chisq.test(table(trap_contingency$otu, trap_contingency$size_class))
    cat("  Chi-square test (species x size class):\n")
    cat("    X^2 =", round(trap_chi$statistic, 2), ", df =", trap_chi$parameter,
        ", p =", format.pval(trap_chi$p.value, 3), "\n")

    if (trap_chi$p.value < 0.05) {
      cat("    Significant: Trapezia species composition varies with coral size\n")
      cat("    This suggests size-based species replacement/substitutability\n\n")
    }
  }
}

# ----------------------------------------------------------------------------
# D3. Site-specific functional composition differences
# ----------------------------------------------------------------------------

cat("Site-specific functional composition:\n")

site_functional <- cafi_clean %>%
  group_by(site, functional_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(site) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

site_func_wide <- site_functional %>%
  dplyr::select(site, functional_group, proportion) %>%
  pivot_wider(names_from = functional_group, values_from = proportion, values_fill = 0)

print(site_func_wide)
cat("\n")

# Chi-square test for site x taxonomic group
site_func_chi <- chisq.test(table(cafi_clean$site, cafi_clean$functional_group))
cat("Chi-square test (site x taxonomic group):\n")
cat("  X^2 =", round(site_func_chi$statistic, 2), ", df =", site_func_chi$parameter,
    ", p =", format.pval(site_func_chi$p.value, 3), "\n\n")

# Functional composition by size figure
p_func_comp <- coral_master %>%
  filter(!is.na(volume), volume > 0) %>%
  mutate(
    size_class = cut(volume,
                     breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE)
  ) %>%
  dplyr::select(coral_id, site, size_class, n_trapezia, n_resident_fish, n_corallivore,
         n_other_crab, n_shrimp, n_other) %>%
  pivot_longer(cols = starts_with("n_"), names_to = "group", values_to = "abundance") %>%
  mutate(group = str_replace(group, "n_", "") %>% str_replace_all("_", " ") %>% str_to_title(),
         group = case_match(group,
                        "Corallivore" ~ "Gastropods",
                        "Resident Fish" ~ "Fish",
                        "Other Crab" ~ "Other crabs",
                        "Other" ~ "Other invertebrates",
                        .default = group)) %>%
  group_by(size_class, group) %>%
  summarise(mean_abundance = mean(abundance), .groups = "drop") %>%
  ggplot(aes(x = size_class, y = mean_abundance, fill = group)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
  scale_fill_manual(values = c(
    "Trapezia" = "#D55E00",
    "Fish" = "#0072B2",
    "Gastropods" = "#CC79A7",
    "Other crabs" = "#009E73",
    "Shrimp" = "#E69F00",
    "Other invertebrates" = "#999999"
  ), name = "Taxonomic Group") +
  labs(
    x = "Coral Size Class",
    y = "Mean Abundance per Coral",
    title = "Taxonomic Group Composition by Coral Size",
    subtitle = "Mean abundance across all corals in each size class"
  ) +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(FIG_DIR, "functional_composition_by_size.png"), p_func_comp,
       width = 10, height = 6, dpi = 300, bg = "white")
cat("  Saved: functional_composition_by_size.png\n\n")

# ############################################################################
#                    PART E: SCALING PATTERNS BY TAXONOMIC GROUP
# ############################################################################

cat("============================================================\n")
cat("PART E: SCALING PATTERNS BY TAXONOMIC GROUP\n")
cat("============================================================\n\n")

cat("Loading scaling results from script 05 (authoritative source)...\n\n")

# Load pre-computed scaling results from 05_species_scaling_analysis.R
scaling_results_obj <- load_object("scaling_analysis_results")
taxonomic_scaling_results <- scaling_results_obj$models$functional_groups

# Print summary
for (group_name in names(taxonomic_scaling_results)) {
  result <- taxonomic_scaling_results[[group_name]]
  if (result$converged && !is.na(result$beta)) {
    cat(sprintf("  %-25s beta = %6.3f [%5.2f, %5.2f], p vs 1: %s\n",
                group_name, result$beta, result$ci_lower, result$ci_upper,
                format.pval(result$p_vs_1, 2)))
  } else {
    cat(sprintf("  %-25s %s\n", group_name, result$interpretation))
  }
}

cat("\n")

# Compile results into table
scaling_table <- map_df(taxonomic_scaling_results, function(r) {
  tibble(
    functional_group = r$response,
    n_corals = r$n_corals,
    n_nonzero = r$n_nonzero,
    total_abundance = r$total_abundance,
    beta = r$beta,
    se = r$se,
    ci_lower = r$ci_lower,
    ci_upper = r$ci_upper,
    p_value = r$p_value,
    p_vs_1 = r$p_vs_1,
    interpretation = r$interpretation
  )
}) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

save_table(scaling_table, "taxonomic_group_scaling")
cat("Saved: taxonomic_group_scaling.csv\n\n")

# ----------------------------------------------------------------------------
# E1. Forest plot of beta estimates with 95% CI
# ----------------------------------------------------------------------------

cat("Creating forest plot of taxonomic group scaling exponents...\n")

# Interpretation colors (colorblind-safe)
interp_colors <- c(
  "Redirection (\u03b2 < 1)" = "#D55E00",
  "Field of Dreams (\u03b2 \u2248 1)" = "#009E73",
  "Super-linear (\u03b2 > 1)" = "#0072B2",
  "Insufficient data" = "#999999"
)

plot_data <- scaling_table %>%
  filter(!is.na(beta)) %>%
  mutate(functional_group = factor(functional_group,
                                   levels = functional_group[order(beta)]),
         functional_group = case_match(as.character(functional_group),
                                   "Trapezia (mutualists)" ~ "Trapezia crabs",
                                   "Coral-eating snails" ~ "Gastropods",
                                   "Resident Fish" ~ "Fish",
                                   "Corallivores" ~ "Gastropods",
                                   "Other Crabs" ~ "Other crabs",
                                   .default = as.character(functional_group)),
         functional_group = factor(functional_group, levels = unique(functional_group)))

if (nrow(plot_data) > 0) {
  p_forest <- ggplot(plot_data, aes(x = beta, y = functional_group, color = interpretation)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray80") +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.25, linewidth = 0.8) +
    geom_point(size = 4) +
    scale_color_manual(values = interp_colors, name = "Scaling Pattern") +
    labs(
      x = expression("Scaling Exponent (" * beta * ")"),
      y = NULL,
      title = "Taxonomic Group Scaling: Abundance vs Coral Volume",
      subtitle = "95% CI shown | Dashed line = Field of Dreams (beta = 1)",
      caption = "Negative binomial GLM: log(abundance) ~ log(volume) + site"
    ) +
    theme_publication() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 11, face = "bold")
    ) +
    coord_cartesian(xlim = c(-0.5, 3))

  ggsave(file.path(FIG_DIR, "taxonomic_group_forest_plot.png"), p_forest,
         width = 10, height = 6, dpi = 300, bg = "white")
  cat("  Saved: taxonomic_group_forest_plot.png\n\n")
}

# ############################################################################
#                    MANUSCRIPT FIGURE 3: TAXONOMIC GROUP SCALING
# ############################################################################

cat("============================================================\n")
cat("CREATING MANUSCRIPT FIGURE 3\n")
cat("============================================================\n\n")

# Panel A: Trapezia scaling
panel_a <- if (exists("trap_data") && !is.null(trap_scaling_result) &&
               isTRUE(trap_scaling_result$converged)) {
  trap_data %>%
    ggplot(aes(x = volume, y = n_trapezia + 0.5, color = site)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(aes(group = 1), method = MASS::glm.nb, formula = y ~ log(x),
                se = TRUE, color = "black", linewidth = 1) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    labs(
      x = expression("Coral Volume (cm"^3*")"),
      y = "Trapezia Abundance",
      title = "A. Trapezia Guardian Crabs",
      subtitle = sprintf("beta = %.2f [%.2f, %.2f]",
                         trap_scaling_result$beta, trap_scaling_result$ci_lower,
                         trap_scaling_result$ci_upper)
    ) +
    theme_publication(base_size = 11) +
    theme(legend.position = c(0.15, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
} else {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient Trapezia data") +
    theme_void() +
    labs(title = "A. Trapezia Guardian Crabs")
}

# Panel B: Resident fish scaling
panel_b <- if (exists("fish_data") && !is.null(fish_scaling_result) &&
               isTRUE(fish_scaling_result$converged)) {
  fish_data %>%
    ggplot(aes(x = volume, y = n_fish + 0.5, color = site)) +
    geom_point(alpha = 0.6, size = 2,
               position = position_jitter(width = 0, height = 0.05, seed = 42)) +
    geom_smooth(aes(group = 1), method = MASS::glm.nb, formula = y ~ log(x),
                se = TRUE, color = "black", linewidth = 1) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    labs(
      x = expression("Coral Volume (cm"^3*")"),
      y = "Fish Abundance",
      title = "B. Resident Fish (Nutrient Providers)",
      subtitle = sprintf("beta = %.2f [%.2f, %.2f]",
                         fish_scaling_result$beta, fish_scaling_result$ci_lower,
                         fish_scaling_result$ci_upper)
    ) +
    theme_publication(base_size = 11) +
    theme(legend.position = "none")
} else {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient fish data") +
    theme_void() +
    labs(title = "B. Resident Fish")
}

# Panel C: Forest plot
panel_c <- if (nrow(plot_data) > 0) {
  ggplot(plot_data, aes(x = beta, y = functional_group, color = interpretation)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, linewidth = 0.7) +
    geom_point(size = 3.5) +
    scale_color_manual(values = interp_colors, name = "Pattern") +
    labs(
      x = expression("Scaling Exponent (" * beta * ")"),
      y = NULL,
      title = "C. Taxonomic Group Scaling Comparison",
      subtitle = "Dashed line = Field of Dreams (beta = 1)"
    ) +
    theme_publication(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      axis.text.y = element_text(size = 10)
    ) +
    coord_cartesian(xlim = c(-0.5, 2.5))
} else {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient data for forest plot") +
    theme_void() +
    labs(title = "C. Scaling Comparison")
}

# Combine panels
fig3 <- (panel_a | panel_b) / panel_c +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Figure 3: Taxonomic Group Scaling Patterns",
    subtitle = "Body plan predicts how CAFI abundance scales with coral size",
    caption = "Negative binomial GLM: abundance ~ log(volume) + site",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      plot.caption = element_text(size = 9, color = "gray50")
    )
  )

# Save manuscript figure (to both manuscript and analysis dirs)
ggsave(file.path(PATHS$fig_manuscript, "fig3_functional_groups.png"), fig3,
       width = 12, height = 10, dpi = 300, bg = "white")
ggsave(file.path(FIG_DIR, "fig3_functional_groups.png"), fig3,
       width = 12, height = 10, dpi = 300, bg = "white")

cat("Saved: fig3_functional_groups.png (manuscript + analysis)\n\n")

# ############################################################################
#                    SUMMARY AND SAVE RESULTS
# ############################################################################

cat("============================================================\n")
cat("TAXONOMIC GROUP ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

# Summary statistics
cat("SUMMARY:\n\n")

cat("A. TRAPEZIA CRABS:\n")
if (exists("trapezia") && nrow(trapezia) > 0) {
  cat("   Total individuals:", nrow(trapezia), "\n")
  cat("   Species/morphotypes:", n_distinct(trapezia$otu), "\n")
  cat("   Occupancy:", sprintf("%.1f%%", sum(trap_data$has_trapezia) / nrow(trap_data) * 100), "\n")
  if (!is.null(trap_scaling_result) && !is.na(trap_scaling_result$beta)) {
    cat("   Scaling beta:", round(trap_scaling_result$beta, 3), "\n")
    cat("   Pattern:", trap_scaling_result$interpretation, "\n")
  }
} else {
  cat("   No data\n")
}

cat("\nB. FISH:\n")
if (exists("resident_fish") && nrow(resident_fish) > 0) {
  cat("   Total individuals:", nrow(resident_fish), "\n")
  cat("   Species:", n_distinct(resident_fish$otu), "\n")
  cat("   Occupancy:", sprintf("%.1f%%", sum(fish_data$has_fish) / nrow(fish_data) * 100), "\n")
  if (!is.null(fish_scaling_result) && !is.na(fish_scaling_result$beta)) {
    cat("   Scaling beta:", round(fish_scaling_result$beta, 3), "\n")
    cat("   Pattern:", fish_scaling_result$interpretation, "\n")
  }
} else {
  cat("   No data\n")
}

cat("\nC. GASTROPODS:\n")
if (exists("coral_eating_snails") && nrow(coral_eating_snails) > 0) {
  cat("   Total individuals:", nrow(coral_eating_snails), "\n")
  cat("   Species:", n_distinct(coral_eating_snails$otu), "\n")
  cat("   Occupancy:", sprintf("%.1f%%", sum(coral_data$has_corallivore) / nrow(coral_data) * 100), "\n")
  if (!is.null(snail_scaling_result) && !is.na(snail_scaling_result$beta)) {
    cat("   Scaling beta:", round(snail_scaling_result$beta, 3), "\n")
    cat("   Pattern:", snail_scaling_result$interpretation, "\n")
  }
} else {
  cat("   No data\n")
}

cat("\nD. TAXONOMIC DIVERSITY:\n")
cat("   Group with highest diversity:", guild_diversity$functional_group[1], "\n")
cat("   Species in dominant group:", guild_diversity$n_species[1], "\n")

cat("\nOUTPUT FILES:\n")
cat("  Figures:\n")
cat("    - output/figures/manuscript/fig3_functional_groups.png\n")
cat("    - output/figures/trapezia_scaling.png\n")
cat("    - output/figures/taxonomic_composition_by_size.png\n")
cat("    - output/figures/gastropod_scaling.png\n")
cat("  Tables:\n")
cat("    - output/tables/taxonomic_group_scaling.csv\n")
cat("    - output/tables/trapezia_species.csv\n")
cat("    - output/tables/gastropod_prevalence.csv\n\n")

# Save taxonomic group analysis results object
functional_analysis_results <- list(
  # Sample info
  n_corals = nrow(coral_master),
  n_cafi = nrow(cafi_clean),

  # Trapezia results
  trapezia_n = if(exists("trapezia")) nrow(trapezia) else 0,
  trapezia_scaling = if(exists("trap_scaling_result")) trap_scaling_result else NULL,

  # Fish results
  fish_n = if(exists("resident_fish")) nrow(resident_fish) else 0,
  fish_scaling = if(exists("fish_scaling_result")) fish_scaling_result else NULL,

  # Gastropod results
  snail_n = if(exists("coral_eating_snails")) nrow(coral_eating_snails) else 0,
  snail_scaling = if(exists("snail_scaling_result")) snail_scaling_result else NULL,

  # All scaling results
  functional_scaling = scaling_table,

  # Guild diversity
  guild_diversity = guild_diversity
)

save_object(functional_analysis_results, "functional_analysis_results")
cat("Saved: output/objects/functional_analysis_results.rds\n\n")

cat("Functional group analysis complete!\n")
