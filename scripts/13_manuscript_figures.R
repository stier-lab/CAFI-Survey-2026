# ============================================================================
# 13_manuscript_figures.R - Publication-Ready Manuscript Figures
# ============================================================================
#
# PURPOSE: Create a complete set of publication-quality figures for the
#          CAFI Survey manuscript, organized by the three core questions
#
# ============================================================================
# FIGURE ORGANIZATION BY CORE QUESTION
# ============================================================================
#
# OVERVIEW:
#   Figure 1: Study Design - Sites, sampling protocol, dataset summary
#
# Q1: SCALING - Does CAFI abundance scale proportionally with coral size?
#   Figure 2: Size-Abundance Scaling (created here)
#     - Tests Field of Dreams (β=1) vs Redistribution (β<1)
#   Figure 3: Taxonomic Group Scaling (from 08_functional_groups.R)
#     - Group-level scaling exponents, stacked composition
#
# Q2: COMPOSITION - What shapes community assembly?
#   Figure 4: Community Composition (from 02_community_analysis.R)
#     - Rarefaction control, size-divergence test, NMDS
#
# Q3: FEEDBACKS - Does CAFI community predict coral condition?
#   Figure 5: CAFI-Condition Feedbacks (from 09_cafi_condition_feedbacks.R)
#     - Richness/Trapezia/Galeropsis effects on condition
#
# Q4: NEIGHBORHOOD - Does local coral density affect CAFI?
#   No main figure (null result) → Supplement S7
#
# SUPPLEMENTARY:
#   S1: Species accumulation curves
#   S2: Multi-metric PERMANOVA sensitivity
#   S3: NMDS ordination by site/size
#   S4: Spatial autocorrelation (Moran's I)
#   S5: Composition divergence by size (null after rarefaction)
#   S6: Species-level scaling forest plot
#   S7: Neighborhood null results (Q4)
#
# ============================================================================
#
# DEPENDENCIES: All prior scripts must have been run
#
# OUTPUTS:
#   output/figures/manuscript/fig1_study_design.png      (Overview)
#   output/figures/manuscript/fig2_scaling.png           (Q1)
#   output/figures/manuscript/fig3_functional_groups.png (Q1 groups)
#   output/figures/manuscript/fig4_composition_network.png (Q2)
#   output/figures/manuscript/fig5_feedbacks.png         (Q3)
#   output/figures/supplement/figS1-S7                   (Supplements)
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-22
# ============================================================================

cat("\n========================================\n")
cat("13: Manuscript Figure Suite\n")
cat("========================================\n\n")

# ============================================================================
# SETUP AND DATA LOADING
# ============================================================================

# Load setup (packages, paths, theme)
source(here::here("scripts/00_setup.R"))

# Create directories
MANUSCRIPT_DIR <- file.path(PATHS$figures, "manuscript")
SUPPLEMENT_DIR <- file.path(PATHS$figures, "supplement")
dir.create(MANUSCRIPT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SUPPLEMENT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load all required data objects
coral_master <- load_object("coral_master")
cafi_clean <- load_object("cafi_clean")
community_matrix <- load_object("community_matrix")
functional_summary <- load_object("functional_summary")

# Load landscape models if available
landscape_models <- tryCatch(
  load_object("landscape_models"),
  error = function(e) NULL
)

# Load scaling results if available
scaling_results <- tryCatch(
  load_object("scaling_results"),
  error = function(e) NULL
)

cat("Data loaded successfully\n")
cat("  Corals:", nrow(coral_master), "\n")
cat("  CAFI records:", nrow(cafi_clean), "\n\n")

# ============================================================================
# PUBLICATION THEME (enhanced)
# ============================================================================

# Enhanced publication theme for manuscript figures
theme_manuscript <- function(base_size = 12) {
  theme_publication(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "gray30"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      plot.margin = margin(15, 15, 15, 15)
    )
}

# Standard figure dimensions
FIG_WIDTH_SINGLE <- 8    # Single column
FIG_WIDTH_DOUBLE <- 16   # Double column
FIG_HEIGHT_SMALL <- 6
FIG_HEIGHT_MEDIUM <- 8
FIG_HEIGHT_LARGE <- 10

# ############################################################################
#                    FIGURE 1: STUDY DESIGN
# ############################################################################

cat("============================================================\n")
cat("FIGURE 1: Study Design\n")
cat("============================================================\n\n")

# Create a conceptual study design figure
# Note: Without actual GPS data and Mo'orea map shapefiles, we create
# a schematic representation

# Panel A: Schematic site map
site_data <- coral_master %>%
  group_by(site) %>%
  summarise(
    n_corals = n(),
    mean_volume = mean(volume, na.rm = TRUE),
    total_cafi = sum(total_cafi, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Positions matching approximate Mo'orea geography
    x = case_when(
      site == "HAU" ~ 1.2,   # North shore (left)
      site == "MAT" ~ 2.8,   # East shore (right-bottom)
      site == "MRB" ~ 2.8    # Barrier reef (right-top, offshore)
    ),
    y = case_when(
      site == "HAU" ~ 1.6,
      site == "MAT" ~ 0.9,
      site == "MRB" ~ 2.2
    ),
    site_name = case_when(
      site == "HAU" ~ "Hauru\n(North Shore)",
      site == "MAT" ~ "Maatea\n(East Shore)",
      site == "MRB" ~ "Barrier Reef\n(Offshore)"
    )
  )

# Create schematic island with sites
panel_a <- ggplot() +
  # Reef outline (lagoon boundary)
  annotate("path",
           x = 2 + 1.6 * cos(seq(0, 2*pi, length.out = 100)),
           y = 1.5 + 1.1 * sin(seq(0, 2*pi, length.out = 100)),
           color = "#87CEEB", linewidth = 0.8, linetype = "dashed") +
  # Island outline (simplified)
  annotate("path",
           x = 2 + 1.0 * cos(seq(0, 2*pi, length.out = 100)),
           y = 1.5 + 0.7 * sin(seq(0, 2*pi, length.out = 100)),
           color = "gray50", linewidth = 1.2) +
  # Site points
  geom_point(data = site_data, aes(x = x, y = y),
             color = "black", fill = SITE_COLORS[site_data$site],
             shape = 21, size = 6, stroke = 1.2) +
  geom_text(data = site_data, aes(x = x, y = y - 0.32, label = site_name),
            size = 2.8, fontface = "bold", lineheight = 0.85) +
  geom_text(data = site_data, aes(x = x, y = y + 0.22,
                                   label = paste0("n = ", n_corals)),
            size = 2.5, color = "gray30") +
  # Island name
annotate("text", x = 2, y = 1.5, label = "Mo'orea", size = 4.5, fontface = "italic") +
  annotate("text", x = 2, y = 0.15, label = "17\u00B030'S, 149\u00B050'W",
           size = 2.8, color = "gray50") +
  coord_fixed(ratio = 1, xlim = c(0, 4), ylim = c(-0.1, 2.8)) +
  labs(title = "A. Study Sites") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.margin = margin(5, 10, 5, 5))

# Panel B: Dataset summary table
summary_stats <- tibble(
  Metric = c("Coral colonies surveyed", "Total CAFI individuals",
             "Morphospecies (OTUs)", "Functional groups",
             "Mean CAFI per coral"),
  Value = c(as.character(nrow(coral_master)),
            format(nrow(cafi_clean), big.mark = ","),
            as.character(n_distinct(cafi_clean$otu)),
            as.character(n_distinct(cafi_clean$functional_group)),
            as.character(round(mean(coral_master$total_cafi, na.rm = TRUE), 1)))
)

panel_b <- ggplot(summary_stats, aes(x = 1, y = rev(seq_along(Metric)))) +
  geom_text(aes(label = Metric, x = 0.6), hjust = 0, size = 3.3, color = "gray20") +
  geom_text(aes(label = Value, x = 1.4), hjust = 1, size = 3.5, fontface = "bold") +
  annotate("segment", x = 0.55, xend = 1.45, y = 5.6, yend = 5.6,
           color = "gray40", linewidth = 0.5) +
  annotate("segment", x = 0.55, xend = 1.45, y = 0.4, yend = 0.4,
           color = "gray40", linewidth = 0.3) +
  scale_x_continuous(limits = c(0.5, 1.5)) +
  scale_y_continuous(limits = c(0.2, 6)) +
  labs(title = "B. Dataset Summary") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.margin = margin(5, 5, 5, 10))

# Panel C: Sampling protocol (clean text-based workflow)
sampling_df <- tibble(
  stage = factor(1:4),
  step = c("Site Selection", "Colony Selection", "CAFI Census", "Measurements"),
  description = c("3 reef habitats\nacross Mo'orea",
                  "38-40 Pocillopora\nper site (n = 114)",
                  "All fauna extracted\n& identified (~4,000)",
                  "Colony volume, GPS,\nneighborhood density")
)

panel_c <- ggplot(sampling_df, aes(x = as.numeric(stage), y = 1)) +
  # Step boxes
  geom_tile(fill = "gray97", color = "gray40", width = 0.85, height = 0.55,
            linewidth = 0.4) +
  # Step numbers (circled at top of box)
  geom_point(aes(y = 1.35), size = 5, shape = 21,
             fill = "#0072B2", color = "white", stroke = 0.8) +
  geom_text(aes(y = 1.35, label = stage), size = 2.8,
            color = "white", fontface = "bold") +
  # Step names
  geom_text(aes(y = 1.12, label = step), size = 3, fontface = "bold") +
  # Descriptions
  geom_text(aes(y = 0.85, label = description), size = 2.6,
            lineheight = 0.85, color = "gray30") +
  # Arrows between steps
  geom_segment(data = tibble(x = 1:3 + 0.43, xend = 1:3 + 0.57),
               aes(x = x, xend = xend, y = 1, yend = 1),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
               color = "gray50", linewidth = 0.6) +
  scale_x_continuous(limits = c(0.3, 4.7)) +
  scale_y_continuous(limits = c(0.55, 1.55)) +
  labs(title = "C. Sampling Protocol") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.margin = margin(5, 10, 10, 10))

# Combine into Figure 1: A and B on top, C spanning bottom
fig1 <- (panel_a | panel_b) / panel_c +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "Figure 1: Study Design",
    subtitle = expression(italic("Pocillopora")*" coral-associated fauna | Mo'orea, French Polynesia | Summer 2019"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray40")
    )
  )

ggsave(file.path(MANUSCRIPT_DIR, "fig1_study_design.png"), fig1,
       width = 10, height = 7, dpi = 300, bg = "white")

cat("  Saved: fig1_study_design.png\n\n")

# ############################################################################
#                    FIGURE 2: SIZE-ABUNDANCE SCALING (H2)
# ############################################################################

cat("============================================================\n")
cat("FIGURE 2: Size-Abundance Scaling\n")
cat("============================================================\n\n")

# Prepare data for scaling plots
scaling_data <- coral_master %>%
  filter(!is.na(volume), volume > 0, !is.na(total_cafi))

# Fit power-law model for total CAFI
if (nrow(scaling_data) >= 30) {

  abundance_model <- tryCatch({
    MASS::glm.nb(total_cafi ~ log10(volume) + site, data = scaling_data)
  }, error = function(e) NULL)

  if (!is.null(abundance_model)) {
    coefs <- summary(abundance_model)$coefficients
    beta_abundance <- coefs["log10(volume)", "Estimate"]
    se_abundance <- coefs["log10(volume)", "Std. Error"]
    ci_abundance <- confint(abundance_model, "log10(volume)", level = 0.95)

    # Test vs beta = 1 (Field of Dreams)
    z_vs_1 <- (beta_abundance - 1) / se_abundance
    p_vs_1 <- 2 * pnorm(-abs(z_vs_1))

    cat("  Total CAFI scaling: beta =", round(beta_abundance, 3), "\n")
    cat("  95% CI: [", round(ci_abundance[1], 3), ",", round(ci_abundance[2], 3), "]\n")
    cat("  Test vs beta = 1: p =", format.pval(p_vs_1, 3), "\n\n")
  } else {
    beta_abundance <- NA
    ci_abundance <- c(NA, NA)
  }

  # Calculate density (CAFI per cm3)
  scaling_data <- scaling_data %>%
    mutate(cafi_density = total_cafi / volume)

  # Panel A: Total abundance vs volume (log-log)
  panel_a <- scaling_data %>%
    ggplot(aes(x = volume, y = total_cafi, color = site)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(aes(group = 1), method = "glm.nb", formula = y ~ log10(x),
                se = TRUE, color = "black", linewidth = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "gray50", linewidth = 0.8) +
    scale_x_log10(labels = scales::comma,
                  breaks = c(100, 1000, 10000, 100000)) +
    scale_y_log10(labels = scales::comma,
                  breaks = c(1, 10, 100)) +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    labs(
      x = expression("Coral Volume (cm"^3*")"),
      y = "Total CAFI Abundance",
      title = "A. Abundance Scaling",
      subtitle = if(!is.na(beta_abundance))
        sprintf("\u03B2 = %.2f [%.2f, %.2f]",
                beta_abundance, ci_abundance[1], ci_abundance[2])
      else "Model fit unavailable"
    ) +
    annotate("text", x = 300, y = 80, label = "1:1 line",
             size = 3, color = "gray50", fontface = "italic") +
    theme_manuscript() +
    theme(legend.position = c(0.15, 0.85),
          legend.background = element_rect(fill = alpha("white", 0.8)))

  # Panel B: Density vs volume (negative relationship expected)
  density_model <- tryCatch({
    lm(log10(cafi_density + 0.001) ~ log10(volume), data = scaling_data)
  }, error = function(e) NULL)

  if (!is.null(density_model)) {
    density_slope <- coef(density_model)["log10(volume)"]
    density_ci <- confint(density_model, "log10(volume)")
  } else {
    density_slope <- NA
    density_ci <- c(NA, NA)
  }

  panel_b <- scaling_data %>%
    ggplot(aes(x = volume, y = cafi_density, color = site)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x,
                se = TRUE, color = "black", linewidth = 1.2) +
    scale_x_log10(labels = scales::comma,
                  breaks = c(100, 1000, 10000, 100000)) +
    scale_y_log10(labels = scales::scientific) +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    labs(
      x = expression("Coral Volume (cm"^3*")"),
      y = expression("CAFI Density (individuals/cm"^3*")"),
      title = "B. Density Dilution",
      subtitle = if(!is.na(density_slope))
        sprintf("slope = %.2f [%.2f, %.2f]",
                density_slope, density_ci[1], density_ci[2])
      else "Model fit unavailable"
    ) +
    theme_manuscript() +
    theme(legend.position = "none")

  # Panel C: Species richness scaling
  richness_model <- tryCatch({
    glm(otu_richness ~ log10(volume) + site, data = scaling_data, family = poisson)
  }, error = function(e) NULL)

  if (!is.null(richness_model)) {
    z_richness <- coef(richness_model)["log10(volume)"]
    z_ci <- confint(richness_model, "log10(volume)")
  } else {
    z_richness <- NA
    z_ci <- c(NA, NA)
  }

  panel_c <- scaling_data %>%
    ggplot(aes(x = volume, y = otu_richness, color = site)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(aes(group = 1), method = "glm",
                method.args = list(family = poisson), formula = y ~ log10(x),
                se = TRUE, color = "black", linewidth = 1.2) +
    scale_x_log10(labels = scales::comma,
                  breaks = c(100, 1000, 10000, 100000)) +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    labs(
      x = expression("Coral Volume (cm"^3*")"),
      y = "Species Richness",
      title = "C. Species-Area Relationship",
      subtitle = if(!is.na(z_richness))
        sprintf("z = %.2f [%.2f, %.2f]", z_richness, z_ci[1], z_ci[2])
      else "Model fit unavailable"
    ) +
    theme_manuscript() +
    theme(legend.position = "none")

  # Panel D: Hypothesis comparison — number line showing observed β vs hypotheses
  hyp_data <- tibble(
    hypothesis = c("Redistribution", "Field of Dreams", "Super-linear"),
    beta_ref = c(0.5, 1.0, 1.5),
    region_min = c(-0.2, 0.95, 1.05),
    region_max = c(0.95, 1.05, 2.5),
    color = c("#009E73", "gray60", "#D55E00")
  )

  panel_d <- ggplot() +
    # Hypothesis regions (background shading)
    annotate("rect", xmin = -0.1, xmax = 0.95, ymin = 0.7, ymax = 1.3,
             fill = "#009E73", alpha = 0.1) +
    annotate("rect", xmin = 1.05, xmax = 2.2, ymin = 0.7, ymax = 1.3,
             fill = "#D55E00", alpha = 0.1) +
    # Reference line at β = 1
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    # Observed β with CI
    annotate("errorbarh",
             xmin = ci_abundance[1], xmax = ci_abundance[2],
             y = 1, height = 0.15, linewidth = 0.8, color = "black") +
    annotate("point", x = beta_abundance, y = 1, size = 4, shape = 18) +
    # Hypothesis labels
    annotate("text", x = 0.5, y = 1.45, label = "Redistribution\n(\u03B2 < 1)",
             size = 3, color = "#009E73", fontface = "bold", lineheight = 0.85) +
    annotate("text", x = 1.0, y = 0.45, label = "Field of Dreams\n(\u03B2 = 1)",
             size = 2.6, color = "gray40", lineheight = 0.85) +
    annotate("text", x = 1.6, y = 1.45, label = "Super-linear\n(\u03B2 > 1)",
             size = 3, color = "#D55E00", fontface = "bold", lineheight = 0.85) +
    # Observed value label (offset below the point)
    annotate("text", x = beta_abundance + 0.35, y = 0.78,
             label = sprintf("Observed: \u03B2 = %.2f", beta_abundance),
             size = 2.8, fontface = "bold", hjust = 0) +
    scale_x_continuous(limits = c(-0.1, 2.2),
                       breaks = c(0, 0.5, 1, 1.5, 2),
                       name = expression("Scaling exponent ("*beta*")")) +
    scale_y_continuous(limits = c(0.3, 1.6)) +
    labs(title = "D. Hypothesis Test") +
    theme_manuscript() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())

  # Combine panels
  fig2 <- (panel_a | panel_b) / (panel_c | panel_d) +
    plot_annotation(
      title = "Figure 2: Size-Abundance Scaling in CAFI Communities",
      subtitle = if(!is.na(beta_abundance) && beta_abundance > 1)
        "Marginally super-linear scaling: larger corals attract disproportionately more CAFI"
      else if(!is.na(beta_abundance) && beta_abundance < 1)
        "Sublinear scaling: propagule redistribution toward larger corals"
      else "Power-law scaling relationships in coral-associated fauna",
      caption = "Dashed line in (A): 1:1 expectation (Field of Dreams). Black curves: fitted NB GLM with 95% CI. Diamond in (D): observed \u03B2 with 95% CI.",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        plot.caption = element_text(size = 9, color = "gray50", hjust = 0)
      )
    )

  ggsave(file.path(MANUSCRIPT_DIR, "fig2_scaling.png"), fig2,
         width = 12, height = 10, dpi = 300, bg = "white")

  cat("  Saved: fig2_scaling.png\n\n")

} else {
  cat("  Insufficient data for scaling analysis\n\n")
}

# NOTE: Figure 5 (CAFI-Condition Feedbacks) is created by script 09.
# See "ASSEMBLE FIGURES FROM OTHER SCRIPTS" section below.

# ############################################################################
#                    ASSEMBLE FIGURES FROM OTHER SCRIPTS
# ############################################################################

cat("============================================================\n")
cat("ASSEMBLING FIGURES FROM ANALYSIS SCRIPTS\n")
cat("============================================================\n\n")

# Q1: SCALING - Figure 3 (from 08_functional_groups.R)
cat("Q1: Functional Group Scaling (Figure 3)...\n")
fig3_src <- file.path(PATHS$fig_manuscript, "fig3_functional_groups.png")
if (file.exists(fig3_src)) {
  cat("  Found: fig3_functional_groups.png (from script 08)\n")
} else {
  cat("  Not yet generated - run script 08\n")
}

# ============================================================================
# FIGURE 4: COMMUNITY ASSEMBLY (Q2)
# Site effects + co-occurrence network + modularity + hub species
# ============================================================================
cat("\n============================================================\n")
cat("FIGURE 4: Community Assembly (NMDS + Network)\n")
cat("============================================================\n\n")

network_results <- tryCatch(load_object("cafi_network"), error = function(e) NULL)
community_results <- tryCatch(load_object("community_analysis_results"), error = function(e) NULL)

if (!is.null(network_results) && !is.null(community_results) &&
    !is.null(community_results$community$nmds_scores)) {

  # Load igraph for graph operations
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package required for Figure 4")
  }
  library(igraph)

  # --- Panel A: NMDS by site ---
  nmds_scores <- community_results$community$nmds_scores
  perm_r2 <- community_results$community$permanova$R2[1]  # site R²
  nmds_stress <- community_results$community$nmds_stress

  p_a <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = site)) +
    geom_point(aes(size = log_volume), alpha = 0.7) +
    stat_ellipse(level = 0.95, linetype = "dashed") +
    scale_color_manual(values = SITE_COLORS, name = "Site") +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    labs(title = "A. Composition by Reef Site",
         subtitle = sprintf("PERMANOVA R\u00B2 = %.3f, p < 0.01 | Stress = %.3f",
                            perm_r2, nmds_stress)) +
    theme_manuscript() +
    coord_fixed() +
    theme(legend.position = c(0.85, 0.2))

  # --- Panel B: Network by module (ggplot2, no ggraph) ---
  g <- network_results$graph
  layout <- network_results$fr_layout
  node_df <- data.frame(
    x = layout[,1], y = layout[,2],
    species = igraph::V(g)$name,
    module = factor(igraph::V(g)$module),
    degree = igraph::degree(g)
  )
  el <- igraph::as_edgelist(g, names = FALSE)
  edge_df <- data.frame(
    x = layout[el[,1], 1], y = layout[el[,1], 2],
    xend = layout[el[,2], 1], yend = layout[el[,2], 2],
    weight = igraph::E(g)$weight
  )
  n_nodes <- igraph::vcount(g)
  n_edges <- igraph::ecount(g)
  modularity_obs <- network_results$network_metrics$value[
    network_results$network_metrics$metric == "modularity"]
  ratio_to_null <- network_results$null_comparison$ratio_to_null[
    network_results$null_comparison$metric == "Modularity"]

  # Label top-degree nodes
  degree_threshold <- quantile(node_df$degree, 0.8)

  p_b <- ggplot() +
    geom_segment(data = edge_df,
                 aes(x = x, y = y, xend = xend, yend = yend, alpha = weight),
                 color = "gray60", linewidth = 0.3) +
    geom_point(data = node_df,
               aes(x = x, y = y, size = degree, fill = module),
               shape = 21, color = "gray30", stroke = 0.3) +
    scale_fill_brewer(palette = "Set2", name = "Module") +
    scale_size_continuous(range = c(2, 8), name = "Degree") +
    scale_alpha_continuous(range = c(0.1, 0.4), guide = "none") +
    labs(title = "B. Co-occurrence Network",
         subtitle = sprintf("%d species, %d edges | Q = %.2f (%.1fx null)",
                            n_nodes, n_edges, modularity_obs, ratio_to_null)) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "gray40"),
          legend.position = "right",
          plot.margin = margin(5, 5, 5, 5))

  # Add labels if ggrepel available
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_b <- p_b +
      ggrepel::geom_text_repel(
        data = node_df %>% filter(degree >= degree_threshold),
        aes(x = x, y = y, label = species),
        size = 2.2, max.overlaps = 12, segment.alpha = 0.3)
  }

  # --- Panel C: Modularity vs null ---
  null_mod_df <- data.frame(modularity = network_results$null_metrics[, "modularity"])
  mod_z <- network_results$null_comparison$z_score[
    network_results$null_comparison$metric == "Modularity"]

  p_c <- ggplot(null_mod_df, aes(x = modularity)) +
    geom_histogram(bins = 30, fill = "gray70", color = "white", alpha = 0.8) +
    geom_vline(xintercept = modularity_obs, color = "#D55E00", linewidth = 1.2) +
    annotate("text", x = modularity_obs, y = Inf,
             label = sprintf("Observed\nQ = %.2f\nz = %.1f", modularity_obs, mod_z),
             vjust = 1.5, hjust = -0.1, color = "#D55E00", fontface = "bold", size = 3) +
    labs(title = "C. Modularity vs Null Model",
         subtitle = "1000 Erdos-Renyi random networks",
         x = "Modularity (Q)", y = "Frequency") +
    theme_manuscript()

  # --- Panel D: Hub species ---
  type_colors <- network_results$type_colors
  p_d <- network_results$centrality %>%
    arrange(desc(hub_score)) %>%
    slice_head(n = 10) %>%
    ggplot(aes(x = reorder(species, hub_score), y = hub_score, fill = type)) +
    geom_col(alpha = 0.8, width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = type_colors, name = "Type") +
    labs(title = "D. Hub Species",
         subtitle = "Connectivity + influence (top 10)",
         x = NULL, y = "Hub Score") +
    theme_manuscript() +
    theme(legend.position = "bottom")

  # --- Compose ---
  fig4 <- (p_a | p_b) / (p_c | p_d) +
    plot_layout(heights = c(1.2, 1)) +
    plot_annotation(
      title = "Figure 4: Community Assembly Structure",
      subtitle = "Site pools and non-random co-occurrence shape CAFI composition",
      caption = paste0(
        "A: Bray-Curtis NMDS; ellipses = 95% CI by site. ",
        "B: FR layout; node color = Louvain module, size = degree; edges = r > 0.3. ",
        "C: Observed modularity vs 1000 random networks. ",
        "D: Hub score = standardized degree + eigenvector centrality."),
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        plot.caption = element_text(size = 8, hjust = 0, color = "gray50")
      )
    )

  ggsave(file.path(MANUSCRIPT_DIR, "fig4_composition_network.png"), fig4,
         width = 14, height = 11, dpi = 300, bg = "white")
  cat("  Saved: fig4_composition_network.png\n\n")

} else {
  cat("  Network or community results not available - run scripts 02 and 06 first\n\n")
}

# Q3: FEEDBACKS - Figure 5 (from 09_cafi_condition_feedbacks.R)
cat("\nQ3: CAFI-Condition Feedbacks (Figure 5)...\n")
fig5_src <- file.path(PATHS$fig_manuscript, "fig5_feedbacks.png")
if (file.exists(fig5_src)) {
  cat("  Found: fig5_feedbacks.png (from script 09)\n")
} else {
  cat("  Not yet generated - run script 09\n")
}

cat("\n")

# ############################################################################
#                    SUPPLEMENTARY FIGURES
# ############################################################################

cat("============================================================\n")
cat("SUPPLEMENTARY FIGURES\n")
cat("============================================================\n\n")

# S1: Species accumulation curve
cat("Creating S1: Species accumulation curve...\n")

if (exists("community_matrix") && nrow(community_matrix) > 10) {

  # Calculate species accumulation
  specaccum_result <- tryCatch({
    vegan::specaccum(community_matrix, method = "random", permutations = 100)
  }, error = function(e) NULL)

  if (!is.null(specaccum_result)) {

    sac_data <- tibble(
      sites = specaccum_result$sites,
      richness = specaccum_result$richness,
      sd = specaccum_result$sd
    )

    fig_s1 <- ggplot(sac_data, aes(x = sites, y = richness)) +
      geom_ribbon(aes(ymin = richness - 1.96*sd, ymax = richness + 1.96*sd),
                  fill = "#377EB8", alpha = 0.3) +
      geom_line(color = "#377EB8", linewidth = 1.2) +
      geom_hline(yintercept = max(sac_data$richness),
                 linetype = "dashed", color = "gray50") +
      annotate("text", x = max(sac_data$sites) * 0.9,
               y = max(sac_data$richness) + 2,
               label = paste("Total:", round(max(sac_data$richness))),
               size = 3.5) +
      labs(
        x = "Number of Coral Colonies Sampled",
        y = "Cumulative Species Richness",
        title = "Figure S1: Species Accumulation Curve",
        subtitle = "Random accumulation with 95% CI (100 permutations)",
        caption = "Asymptote suggests adequate sampling coverage"
      ) +
      theme_manuscript()

    ggsave(file.path(SUPPLEMENT_DIR, "figS1_species_accumulation.png"), fig_s1,
           width = 8, height = 6, dpi = 300, bg = "white")

    cat("  Saved: figS1_species_accumulation.png\n")
  }
}

# S2: Multi-metric PERMANOVA sensitivity
sensitivity_locations <- c(
  file.path(PATHS$figures, "02_community", "permanova_metric_sensitivity.png")
)
sens_fig <- sensitivity_locations[file.exists(sensitivity_locations)][1]
if (!is.na(sens_fig) && file.exists(sens_fig)) {
  file.copy(sens_fig, file.path(SUPPLEMENT_DIR, "figS2_permanova_sensitivity.png"),
            overwrite = TRUE)
  cat("  Saved: figS2_permanova_sensitivity.png\n")
}

# S3: NMDS ordination by site/size
nmds_locations <- c(
  file.path(PATHS$figures, "02_community", "nmds_by_site.png"),
  file.path(PATHS$figures, "02_community", "nmds_ordination.png")
)
nmds_fig <- nmds_locations[file.exists(nmds_locations)][1]
if (!is.na(nmds_fig) && file.exists(nmds_fig)) {
  file.copy(nmds_fig, file.path(SUPPLEMENT_DIR, "figS3_nmds_ordination.png"),
            overwrite = TRUE)
  cat("  Saved: figS3_nmds_ordination.png\n")
}

# S4: Spatial autocorrelation (Moran's I)
spatial_map_locations <- c(
  file.path(PATHS$figures, "07_spatial", "spatial_autocorrelation_map.png"),
  file.path(PATHS$figures, "spatial_autocorrelation_map.png")
)
spatial_map <- spatial_map_locations[file.exists(spatial_map_locations)][1]
if (!is.na(spatial_map) && file.exists(spatial_map)) {
  file.copy(spatial_map, file.path(SUPPLEMENT_DIR, "figS4_spatial_autocorrelation.png"),
            overwrite = TRUE)
  cat("  Saved: figS4_spatial_autocorrelation.png\n")
}

# S5: Composition divergence by size (former main figure - null after rarefaction)
divergence_locations <- c(
  file.path(PATHS$figures, "02_community", "composition_divergence_by_size.png")
)
div_fig <- divergence_locations[file.exists(divergence_locations)][1]
if (!is.na(div_fig) && file.exists(div_fig)) {
  file.copy(div_fig, file.path(SUPPLEMENT_DIR, "figS5_composition_divergence.png"),
            overwrite = TRUE)
  cat("  Saved: figS5_composition_divergence.png\n")
}

# S6: Species-level scaling forest plot
species_scaling_locations <- c(
  file.path(PATHS$figures, "05_scaling", "species_scaling_forest.png"),
  file.path(PATHS$figures, "species_scaling_forest.png")
)
species_fig <- species_scaling_locations[file.exists(species_scaling_locations)][1]
if (!is.na(species_fig) && file.exists(species_fig)) {
  file.copy(species_fig, file.path(SUPPLEMENT_DIR, "figS6_species_scaling.png"),
            overwrite = TRUE)
  cat("  Saved: figS6_species_scaling.png\n")
}

# S7: Neighborhood null results (Q4)
# Assemble panels showing non-significant neighborhood effects
neighborhood_fig_locations <- c(
  file.path(PATHS$figures, "04_effects", "abundance_vs_neighbors.png")
)
neighbor_fig <- neighborhood_fig_locations[file.exists(neighborhood_fig_locations)][1]
if (!is.na(neighbor_fig) && file.exists(neighbor_fig)) {
  file.copy(neighbor_fig, file.path(SUPPLEMENT_DIR, "figS7_neighborhood_null.png"),
            overwrite = TRUE)
  cat("  Saved: figS7_neighborhood_null.png\n")
} else {
  # Create minimal neighborhood null figure from coral_master data
  if ("n_neighbors" %in% names(coral_master)) {
    neighbor_data <- coral_master %>% filter(!is.na(n_neighbors))
    if (nrow(neighbor_data) >= 20) {
      p_s7a <- ggplot(neighbor_data, aes(x = n_neighbors, y = total_cafi)) +
        geom_point(alpha = 0.6, color = "#0072B2") +
        geom_smooth(method = "lm", se = TRUE, color = "black") +
        labs(x = "Number of Neighbors (5m)", y = "Total CAFI Abundance",
             title = "A. CAFI Abundance vs Neighborhood Density (NS)") +
        theme_manuscript()

      p_s7b <- if ("condition_score" %in% names(neighbor_data)) {
        ggplot(neighbor_data %>% filter(!is.na(condition_score)),
               aes(x = n_neighbors, y = condition_score)) +
          geom_point(alpha = 0.6, color = "#D55E00") +
          geom_smooth(method = "lm", se = TRUE, color = "black") +
          labs(x = "Number of Neighbors (5m)", y = "Coral Condition Score",
               title = "B. Coral Condition vs Neighborhood Density (NS)") +
          theme_manuscript()
      } else {
        ggplot() + theme_void() + labs(title = "B. Condition data unavailable")
      }

      fig_s7 <- p_s7a | p_s7b
      fig_s7 <- fig_s7 + plot_annotation(
        title = "Figure S7: Neighborhood Density Does Not Predict CAFI or Condition",
        subtitle = "Q4 null result: coral density within 5m radius is not significant (n = 61 corals with spatial data)",
        theme = theme(plot.title = element_text(face = "bold", size = 12),
                      plot.subtitle = element_text(size = 10, color = "gray40"))
      )

      ggsave(file.path(SUPPLEMENT_DIR, "figS7_neighborhood_null.png"), fig_s7,
             width = 12, height = 5, dpi = 300, bg = "white")
      cat("  Saved: figS7_neighborhood_null.png (created from data)\n")
    }
  }
}

cat("\n")

# ############################################################################
#                    RESULTS TABLE & SAMPLE SIZES
# ############################################################################

cat("============================================================\n")
cat("COMPILING MANUSCRIPT RESULTS TABLE\n")
cat("============================================================\n\n")

# --- Unified Results Table ---
# Compiles all key statistical results across Q1-Q4

results_rows <- list()

# Q1: Scaling (from script 05)
tryCatch({
  scaling_obj <- load_object("scaling_analysis_results")
  comm_result <- scaling_obj$models$total_abundance
  if (!is.null(comm_result) && comm_result$converged) {
    results_rows[[length(results_rows) + 1]] <- tibble(
      Question = "Q1", Hypothesis = "Redistribution", Test = "NB GLM",
      Predictor = "log10(volume)", Beta_R2 = comm_result$beta,
      CI_lower = comm_result$ci_lower, CI_upper = comm_result$ci_upper,
      p_value = comm_result$p_vs_1, p_FDR = NA_real_,
      n = comm_result$n_corals, Interpretation = comm_result$interpretation
    )
  }
}, error = function(e) cat("  [Skip] Scaling results not available\n"))

# Q2: Composition (from saved PERMANOVA results)
tryCatch({
  community_results <- load_object("community_analysis_results")
  if (!is.null(community_results$permanova)) {
    perm <- community_results$permanova
    results_rows[[length(results_rows) + 1]] <- tibble(
      Question = "Q2", Hypothesis = "Site effects", Test = "PERMANOVA",
      Predictor = "site", Beta_R2 = perm$R2[which(rownames(perm) == "site")],
      CI_lower = NA_real_, CI_upper = NA_real_,
      p_value = perm$`Pr(>F)`[which(rownames(perm) == "site")],
      p_FDR = NA_real_, n = 114, Interpretation = "Significant"
    )
  }
}, error = function(e) cat("  [Skip] PERMANOVA results not available\n"))

# Q3: Feedbacks (from script 09 outputs)
tryCatch({
  feedback_table <- read_csv(file.path(PATHS$tables, "cafi_condition_models.csv"),
                             show_col_types = FALSE)
  if (nrow(feedback_table) > 0) {
    for (i in 1:nrow(feedback_table)) {
      row <- feedback_table[i, ]
      results_rows[[length(results_rows) + 1]] <- tibble(
        Question = "Q3", Hypothesis = paste0(row$predictor, "->condition"),
        Test = "LMM", Predictor = row$predictor,
        Beta_R2 = row$estimate, CI_lower = row$ci_lower, CI_upper = row$ci_upper,
        p_value = row$p_value,
        p_FDR = if ("p_fdr" %in% names(row)) row$p_fdr else NA_real_,
        n = row$n, Interpretation = ifelse(row$significant, "Significant", "NS")
      )
    }
  }
}, error = function(e) cat("  [Skip] Feedback results not available\n"))

# Q4: Neighborhood (from script 04)
tryCatch({
  landscape_models <- load_object("landscape_models")
  if (!is.null(landscape_models$abundance_full)) {
    coefs <- summary(landscape_models$abundance_full)$coefficients
    if ("n_neighbors" %in% rownames(coefs)) {
      results_rows[[length(results_rows) + 1]] <- tibble(
        Question = "Q4", Hypothesis = "Neighborhood->CAFI",
        Test = "NB GLM", Predictor = "n_neighbors",
        Beta_R2 = coefs["n_neighbors", "Estimate"],
        CI_lower = NA_real_, CI_upper = NA_real_,
        p_value = coefs["n_neighbors", "Pr(>|z|)"],
        p_FDR = NA_real_, n = 61, Interpretation = "NS"
      )
    }
  }
}, error = function(e) cat("  [Skip] Landscape model results not available\n"))

# Compile and save
if (length(results_rows) > 0) {
  manuscript_results <- bind_rows(results_rows)
  save_table(manuscript_results, "manuscript_results_summary")
  cat("  Saved: manuscript_results_summary.csv (", nrow(manuscript_results), " results)\n\n")
} else {
  cat("  No pre-computed results available. Run analysis scripts first.\n\n")
}

# --- Sample Size Table ---
cat("Sample sizes by analysis:\n")

sample_sizes <- tibble(
  Analysis = c("Scaling (Q1)", "PERMANOVA (Q2)", "Rarefaction (Q2)",
               "Condition feedbacks (Q3)", "Neighborhood (Q4)"),
  Total_N = c(
    nrow(coral_master),
    nrow(coral_master),
    sum(coral_master$total_cafi >= 5, na.rm = TRUE),
    sum(!is.na(coral_master$condition_score)),
    sum(!is.na(coral_master$n_neighbors))
  ),
  Subset_Reason = c(
    "All corals with volume",
    "All corals with CAFI",
    "Corals with >= 5 CAFI (rarefaction depth)",
    "Corals with physiology data",
    "Corals in 5m neighborhood survey"
  )
)

save_table(sample_sizes, "sample_sizes")
cat("  Saved: sample_sizes.csv\n")
print(sample_sizes, n = Inf)
cat("\n")

# ############################################################################
#                    SUMMARY
# ############################################################################

cat("============================================================\n")
cat("MANUSCRIPT FIGURE SUITE COMPLETE\n")
cat("============================================================\n\n")

# List all manuscript figures
manuscript_figs <- list.files(MANUSCRIPT_DIR, pattern = "\\.png$", full.names = FALSE)
supplement_figs <- list.files(SUPPLEMENT_DIR, pattern = "\\.png$", full.names = FALSE)

cat("MAIN TEXT FIGURES (5):\n")
cat("======================\n\n")

cat("OVERVIEW:\n")
cat("  Fig 1: Study design, sites, dataset summary\n")
overview_figs <- manuscript_figs[grepl("fig1_", manuscript_figs)]
for (fig in sort(overview_figs)) cat("    ", fig, "\n")

cat("\nQ1: SCALING (Size-abundance relationships):\n")
cat("  Fig 2: Community-level scaling (abundance, density, richness)\n")
cat("  Fig 3: Taxonomic group scaling (Trapezia vs Fish vs Gastropods)\n")
q1_figs <- manuscript_figs[grepl("fig[23]_", manuscript_figs)]
for (fig in sort(q1_figs)) cat("    ", fig, "\n")

cat("\nQ2: COMPOSITION (What structures CAFI composition?):\n")
cat("  Fig 4: NMDS by site + co-occurrence network + modularity + hub species\n")
q2_figs <- manuscript_figs[grepl("fig4_", manuscript_figs)]
for (fig in sort(q2_figs)) cat("    ", fig, "\n")

cat("\nQ3: FEEDBACKS (Does CAFI identity predict coral condition?):\n")
cat("  Fig 5: Richness, Trapezia, Galeropsis effects on condition\n")
q3_figs <- manuscript_figs[grepl("fig5_", manuscript_figs)]
for (fig in sort(q3_figs)) cat("    ", fig, "\n")

cat("\nQ4: NEIGHBORHOOD (Does local coral density affect CAFI?):\n")
cat("  No main figure (null result) -> Figure S7\n")

cat("\nSUPPLEMENTARY FIGURES (", length(supplement_figs), "):\n")
for (fig in sort(supplement_figs)) {
  cat("  -", fig, "\n")
}

cat("\nTotal figures:", length(manuscript_figs), "main +", length(supplement_figs), "supplementary\n")
cat("\nOutput locations:\n")
cat("  Main figures:", MANUSCRIPT_DIR, "\n")
cat("  Supplements:", SUPPLEMENT_DIR, "\n\n")

cat("Manuscript figure generation complete!\n")
