#!/usr/bin/env Rscript
# ============================================================================
# 03_fig3_functional_groups.R - FIGURE 3: Functional & Taxonomic Group Responses
# ============================================================================
#
# PURPOSE: Identify which taxa and functional groups drive landscape patterns.
# Test taxon-specific scaling relationships with coral size and proximity.
#
# FUNCTIONAL GROUPS (biologically understood):
#   - Trapezia (Mutualist Defenders): protect corals from predators, sediments
#   - Resident Fishes (Nutrient Providers): enhance coral via ammonium excretion
#   - Corallivores (Drupella snails): negative effects - tissue loss, disease
#   - Other Crabs: included statistically, functional role less clear
#   - Shrimp: various roles, included statistically
#
# HYPOTHESES:
#   H1: Trapezia abundance increases with coral size; may not increase with isolation
#   H2: Resident fishes show positive size scaling but weak proximity effects
#   H3: Corallivores prefer larger colonies but not isolated ones
#   H4: Taxon-specific slopes reveal mechanistic scaling differences
#
# STATISTICAL APPROACHES:
#   - Group-level models: NB GLM for each functional group
#   - Species-specific slopes: GLMs for top 15-20 taxa
#   - Heatmap: Species x slope coefficients for visualization
#
# Author: CAFI Analysis Pipeline
# Date: 2025-12-03
# ============================================================================

cat("\n========================================\n")
cat("FIGURE 3: Functional Group Responses\n")
cat("========================================\n\n")

# Load setup and data
source(here::here("scripts/publication/00_setup.R"))

# Load pre-processed data
coral_master <- load_object("coral_master.rds")
cafi_clean <- load_object("cafi_clean.rds")

# Output directory
FIG_DIR <- FIGURE_DIRS$fig3
cat("Outputs will be saved to:", FIG_DIR, "\n\n")

# ============================================================================
# PREPARE ANALYSIS DATA
# ============================================================================

cat("Preparing analysis data...\n")

# Prepare coral-level data with functional group counts
analysis_data <- coral_master %>%
  filter(
    !is.na(volume), volume > 0,
    !is.na(mean_neighbor_dist)
  ) %>%
  mutate(
    log_volume = log10(volume),
    proximity_m = mean_neighbor_dist / 100,

    # Standardize for model stability
    volume_z = scale(log_volume)[, 1],
    proximity_z = scale(proximity_m)[, 1]
  )

cat("  - Sample size:", nrow(analysis_data), "corals\n\n")

# ============================================================================
# ANALYSIS 1: FUNCTIONAL GROUP SCALING WITH SIZE
# ============================================================================

cat("1. Testing functional group scaling with coral size...\n\n")

# List of functional groups to analyze
functional_groups <- c("n_trapezia", "n_resident_fish", "n_corallivore",
                       "n_other_crab", "n_shrimp")

group_labels <- c(
  "n_trapezia" = "Trapezia (Defenders)",
  "n_resident_fish" = "Resident Fish",
  "n_corallivore" = "Corallivores (Drupella)",
  "n_other_crab" = "Other Crabs",
  "n_shrimp" = "Shrimp"
)

# Fit models for each group
group_size_models <- list()
group_size_results <- list()

for (grp in functional_groups) {
  grp_label <- group_labels[grp]

  # Check if group has sufficient data
  n_present <- sum(analysis_data[[grp]] > 0, na.rm = TRUE)

  if (n_present >= 10) {
    # Fit negative binomial GLM with fallback to Poisson
    model <- tryCatch({
      suppressWarnings(
        glm.nb(reformulate(c("log_volume", "site"), response = grp),
               data = analysis_data, maxit = 100)
      )
    }, error = function(e) {
      # Fall back to Poisson if NB fails
      tryCatch({
        glm(reformulate(c("log_volume", "site"), response = grp),
            data = analysis_data, family = poisson())
      }, error = function(e2) NULL)
    })

    if (!is.null(model)) {
      group_size_models[[grp]] <- model

      # Extract size coefficient
      coef_info <- extract_model_stats(model, "log_volume")

      if (!is.null(coef_info)) {
        group_size_results[[grp]] <- tibble(
          group = grp_label,
          n_present = n_present,
          size_effect = coef_info$estimate,
          se = coef_info$se,
          conf_low = coef_info$conf.low,
          conf_high = coef_info$conf.high,
          p_value = coef_info$p.value,
          significant = coef_info$significant
        )

        cat(sprintf("  %s: size effect = %.3f, p %s (n = %d corals)\n",
                    grp_label, coef_info$estimate, format_p(coef_info$p.value), n_present))
      } else {
        cat(sprintf("  %s: model fitted but coefficient extraction failed (n = %d)\n",
                    grp_label, n_present))
      }
    } else {
      cat(sprintf("  %s: model fitting failed (n = %d)\n", grp_label, n_present))
    }
  } else {
    cat(sprintf("  %s: insufficient data (n = %d)\n", grp_label, n_present))
  }
}

group_size_df <- bind_rows(group_size_results)

# ============================================================================
# ANALYSIS 2: FUNCTIONAL GROUP SCALING WITH PROXIMITY
# ============================================================================

cat("\n2. Testing functional group scaling with proximity...\n\n")

group_prox_models <- list()
group_prox_results <- list()

for (grp in functional_groups) {
  grp_label <- group_labels[grp]
  n_present <- sum(analysis_data[[grp]] > 0, na.rm = TRUE)

  if (n_present >= 10) {
    model <- tryCatch({
      suppressWarnings(
        glm.nb(reformulate(c("proximity_m", "log_volume", "site"), response = grp),
               data = analysis_data, maxit = 100)
      )
    }, error = function(e) {
      tryCatch({
        glm(reformulate(c("proximity_m", "log_volume", "site"), response = grp),
            data = analysis_data, family = poisson())
      }, error = function(e2) NULL)
    })

    if (!is.null(model)) {
      group_prox_models[[grp]] <- model
      coef_info <- extract_model_stats(model, "proximity_m")

      if (!is.null(coef_info)) {
        group_prox_results[[grp]] <- tibble(
          group = grp_label,
          n_present = n_present,
          prox_effect = coef_info$estimate,
          se = coef_info$se,
          conf_low = coef_info$conf.low,
          conf_high = coef_info$conf.high,
          p_value = coef_info$p.value,
          significant = coef_info$significant
        )

        cat(sprintf("  %s: proximity effect = %.3f, p %s\n",
                    grp_label, coef_info$estimate, format_p(coef_info$p.value)))
      }
    }
  }
}

group_prox_df <- bind_rows(group_prox_results)

# ============================================================================
# FIGURE 3A: FUNCTIONAL GROUP SIZE EFFECTS (FOREST PLOT)
# ============================================================================

cat("\n3. Creating Figure 3A: Functional group size effects...\n")

p_3a <- group_size_df %>%
  mutate(group = factor(group, levels = rev(group))) %>%
  ggplot(aes(x = size_effect, y = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high),
                 height = 0.2, linewidth = 1) +
  geom_point(aes(color = significant, size = n_present)) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "#0072B2"),
                     labels = c("Not significant", "p < 0.05"),
                     name = "") +
  scale_size_continuous(range = c(3, 8), name = "N corals\nwith group") +
  labs(
    title = "A. Size Scaling by Functional Group",
    subtitle = "Effect of log(coral volume) on group abundance | Controlling for site",
    x = "Size Effect (log-scale coefficient)",
    y = ""
  ) +
  theme_publication() +
  theme(legend.position = "right")

save_pub_fig(p_3a, "fig3a_functional_size_effects.png", dir = FIG_DIR,
             width = 10, height = 6)

# ============================================================================
# FIGURE 3B: FUNCTIONAL GROUP PROXIMITY EFFECTS (FOREST PLOT)
# ============================================================================

cat("4. Creating Figure 3B: Functional group proximity effects...\n")

p_3b <- group_prox_df %>%
  mutate(group = factor(group, levels = rev(group))) %>%
  ggplot(aes(x = prox_effect, y = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high),
                 height = 0.2, linewidth = 1) +
  geom_point(aes(color = significant, size = n_present)) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "#D55E00"),
                     labels = c("Not significant", "p < 0.05"),
                     name = "") +
  scale_size_continuous(range = c(3, 8), name = "N corals\nwith group") +
  labs(
    title = "B. Proximity Effects by Functional Group",
    subtitle = "Effect of neighbor distance on group abundance | Controlling for size + site",
    x = "Proximity Effect (coefficient per meter)",
    y = ""
  ) +
  theme_publication() +
  theme(legend.position = "right")

save_pub_fig(p_3b, "fig3b_functional_proximity_effects.png", dir = FIG_DIR,
             width = 10, height = 6)

# ============================================================================
# ANALYSIS 3: SPECIES-SPECIFIC SCALING SLOPES
# ============================================================================

cat("\n5. Calculating species-specific scaling slopes...\n")

# Get top 20 most abundant OTUs
top_otus <- cafi_clean %>%
  count(otu, functional_group, sort = TRUE) %>%
  slice_head(n = 20)

cat("  Analyzing top", nrow(top_otus), "OTUs...\n")

# Create OTU abundance matrix per coral
otu_by_coral <- cafi_clean %>%
  filter(otu %in% top_otus$otu) %>%
  group_by(coral_id, otu) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(id_cols = coral_id, names_from = otu,
              values_from = count, values_fill = 0)

# Merge with coral characteristics
otu_analysis <- analysis_data %>%
  select(coral_id, log_volume, proximity_m, site) %>%
  left_join(otu_by_coral, by = "coral_id") %>%
  mutate(across(any_of(top_otus$otu), ~replace_na(., 0)))

# Fit models for each OTU
otu_size_results <- list()
otu_prox_results <- list()

for (i in 1:nrow(top_otus)) {
  otu_name <- top_otus$otu[i]
  func_grp <- top_otus$functional_group[i]

  if (!otu_name %in% names(otu_analysis)) next

  n_present <- sum(otu_analysis[[otu_name]] > 0, na.rm = TRUE)

  if (n_present >= 5) {
    # Size model - use quasipossion for overdispersion robustness
    size_model <- tryCatch({
      suppressWarnings(
        glm(reformulate(c("log_volume", "site"), response = paste0("`", otu_name, "`")),
            data = otu_analysis, family = quasipoisson())
      )
    }, error = function(e) NULL)

    if (!is.null(size_model)) {
      size_coef <- tryCatch({
        extract_model_stats(size_model, "log_volume")
      }, error = function(e) NULL)

      if (!is.null(size_coef)) {
        otu_size_results[[otu_name]] <- tibble(
          otu = otu_name,
          functional_group = func_grp,
          n_present = n_present,
          total_abundance = top_otus$n[i],
          size_slope = size_coef$estimate,
          size_se = size_coef$se,
          size_p = size_coef$p.value
        )
      }
    }

    # Proximity model
    prox_model <- tryCatch({
      suppressWarnings(
        glm(reformulate(c("proximity_m", "log_volume", "site"),
                        response = paste0("`", otu_name, "`")),
            data = otu_analysis, family = quasipoisson())
      )
    }, error = function(e) NULL)

    if (!is.null(prox_model)) {
      prox_coef <- tryCatch({
        extract_model_stats(prox_model, "proximity_m")
      }, error = function(e) NULL)

      if (!is.null(prox_coef)) {
        otu_prox_results[[otu_name]] <- tibble(
          otu = otu_name,
          functional_group = func_grp,
          prox_slope = prox_coef$estimate,
          prox_se = prox_coef$se,
          prox_p = prox_coef$p.value
        )
      }
    }
  }
}

otu_size_df <- bind_rows(otu_size_results)
otu_prox_df <- bind_rows(otu_prox_results)

# Merge size and proximity results
otu_slopes <- otu_size_df %>%
  left_join(otu_prox_df %>% select(otu, prox_slope, prox_se, prox_p),
            by = "otu")

cat("  Calculated slopes for", nrow(otu_slopes), "OTUs\n")

# ============================================================================
# FIGURE 3C: SPECIES-SPECIFIC SLOPES HEATMAP
# ============================================================================

cat("\n6. Creating Figure 3C: Species-specific slopes heatmap...\n")

heatmap_data <- otu_slopes %>%
  mutate(
    # Significance markers
    size_sig = ifelse(size_p < 0.05, "*", ""),
    prox_sig = ifelse(prox_p < 0.05, "*", "")
  ) %>%
  select(otu, functional_group, size_slope, prox_slope) %>%
  pivot_longer(cols = c(size_slope, prox_slope),
               names_to = "predictor",
               values_to = "slope") %>%
  mutate(
    predictor = case_when(
      predictor == "size_slope" ~ "Coral Size",
      predictor == "prox_slope" ~ "Neighbor Proximity"
    )
  )

# Order OTUs by functional group and size slope
otu_order <- otu_slopes %>%
  arrange(functional_group, desc(size_slope)) %>%
  pull(otu)

p_3c <- heatmap_data %>%
  mutate(otu = factor(otu, levels = rev(otu_order))) %>%
  ggplot(aes(x = predictor, y = otu, fill = slope)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", slope)),
            color = ifelse(abs(heatmap_data$slope) > 1, "white", "black"),
            size = 3) +
  scale_fill_gradient2(
    low = "#D55E00", mid = "white", high = "#0072B2",
    midpoint = 0, name = "Slope",
    limits = c(-max(abs(heatmap_data$slope), na.rm = TRUE),
               max(abs(heatmap_data$slope), na.rm = TRUE))
  ) +
  facet_grid(functional_group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "C. Species-Specific Scaling Slopes",
    subtitle = "Effect of coral size and neighbor proximity on individual OTU abundance",
    x = "Predictor",
    y = "OTU (morphotype)"
  ) +
  theme_publication() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 9),
    legend.position = "right"
  )

save_pub_fig(p_3c, "fig3c_species_slopes_heatmap.png", dir = FIG_DIR,
             width = 10, height = 12)

# ============================================================================
# FIGURE 3D: SIZE SLOPES VS PROXIMITY SLOPES
# ============================================================================

cat("\n7. Creating Figure 3D: Size vs proximity slopes scatter...\n")

p_3d <- otu_slopes %>%
  ggplot(aes(x = size_slope, y = prox_slope)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = functional_group, size = total_abundance),
             alpha = 0.7) +
  geom_text_repel(aes(label = otu), size = 2.5, max.overlaps = 15) +
  scale_color_functional() +
  scale_size_continuous(range = c(2, 8), name = "Total\nabundance") +
  labs(
    title = "D. Taxon-Specific Scaling Patterns",
    subtitle = "Each point is one OTU | Size vs proximity scaling relationship",
    x = "Size Slope (log-volume coefficient)",
    y = "Proximity Slope (distance coefficient)",
    color = "Functional\nGroup"
  ) +
  theme_publication()

save_pub_fig(p_3d, "fig3d_size_vs_proximity_slopes.png", dir = FIG_DIR,
             width = 10, height = 8)

# ============================================================================
# COMBINED FIGURE 3 PANEL
# ============================================================================

cat("\n8. Creating combined Figure 3 panel...\n")

fig3_combined <- (p_3a + p_3b) /
  (p_3c + p_3d) +
  plot_annotation(
    title = "Figure 3. Functional and Taxonomic Group Responses to Landscape",
    subtitle = "Group-level and species-specific scaling with coral size and neighbor proximity",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray30")
    )
  ) +
  plot_layout(heights = c(1, 1.5))

save_pub_fig(fig3_combined, "fig3_functional_groups_combined.png",
             dir = FIG_DIR, width = 16, height = 18)

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n9. Saving results...\n")

# Save tables
save_table(group_size_df, "fig3_functional_group_size_effects.csv")
save_table(group_prox_df, "fig3_functional_group_proximity_effects.csv")
save_table(otu_slopes, "fig3_otu_specific_slopes.csv")

# Save models
save_object(list(
  group_size_models = group_size_models,
  group_prox_models = group_prox_models
), "fig3_models.rds")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("Figure 3 Analysis Summary\n")
cat("========================================\n\n")

cat("FUNCTIONAL GROUP SIZE SCALING:\n")
for (i in 1:nrow(group_size_df)) {
  sig <- if (group_size_df$significant[i]) " *" else ""
  cat(sprintf("  %s: %.2f%s\n",
              group_size_df$group[i],
              group_size_df$size_effect[i], sig))
}

cat("\nFUNCTIONAL GROUP PROXIMITY EFFECTS:\n")
for (i in 1:nrow(group_prox_df)) {
  sig <- if (group_prox_df$significant[i]) " *" else ""
  cat(sprintf("  %s: %.2f%s\n",
              group_prox_df$group[i],
              group_prox_df$prox_effect[i], sig))
}

cat("\nKEY FINDINGS:\n")

# Find groups with strongest size effects
strongest_size <- group_size_df %>%
  filter(significant) %>%
  arrange(desc(abs(size_effect)))

if (nrow(strongest_size) > 0) {
  cat(sprintf("  - Strongest size effect: %s (%.2f)\n",
              strongest_size$group[1], strongest_size$size_effect[1]))
}

# Find groups with strongest proximity effects
strongest_prox <- group_prox_df %>%
  filter(significant) %>%
  arrange(desc(abs(prox_effect)))

if (nrow(strongest_prox) > 0) {
  cat(sprintf("  - Strongest proximity effect: %s (%.2f)\n",
              strongest_prox$group[1], strongest_prox$prox_effect[1]))
}

cat("\nFigures saved to:", FIG_DIR, "\n")
cat("Tables saved to:", TABLES_DIR, "\n\n")
cat("Figure 3 analysis complete!\n\n")
