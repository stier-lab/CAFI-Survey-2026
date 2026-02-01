# ============================================================================
# Figure 2: Landscape Effects on CAFI - Size-Abundance Scaling
# ============================================================================
#
# 3-panel figure showing how coral volume affects CAFI abundance and diversity
#
# A: Abundance ~ Volume (NB GLM fit) - tests Field of Dreams hypothesis
# B: Richness ~ Volume (GLM fit) - species-area relationship
# C: Top 10 species scaling (forest plot) - heterogeneous responses
#
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
  library(patchwork)
  library(scales)
  library(ggtext)   # For italicized species names with richtext
  library(cowplot)  # For plot composition and legend extraction
})

# Load data
coral_master <- readRDS("output/objects/coral_master.rds")
scaling_results <- readRDS("output/objects/scaling_analysis_results.rds")

# Site colors (colorblind-safe Okabe-Ito palette)
site_colors <- c(HAU = "#E69F00", MAT = "#0072B2", MRB = "#009E73")

# Theme for publication
theme_fig2 <- function() {

  theme_minimal(base_size = 11) +
    theme(
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(-0.02, 1.02),
      plot.title = element_text(face = "bold", size = 11, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 9, color = "gray40", margin = margin(b = 5)),
      axis.title = element_text(size = 10),
      axis.title.y = element_text(margin = margin(r = 5)),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.text = element_text(size = 9, color = "gray20"),
      axis.line = element_line(color = "gray30", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
      legend.position = "bottom",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.margin = margin(t = 0),
      plot.margin = margin(5, 8, 5, 5)
    )
}

# ============================================================================
# PANEL A: Abundance ~ Volume
# ============================================================================
cat("Building Panel A: Abundance scaling...\n")

# Fit NB GLM
mod_abundance <- glm.nb(total_cafi ~ log10(volume_field) + site, data = coral_master)

# Create prediction data
pred_data_a <- expand.grid(

  volume_field = 10^seq(log10(min(coral_master$volume_field, na.rm = TRUE)),
                        log10(max(coral_master$volume_field, na.rm = TRUE)),
                        length.out = 100),
  site = "HAU"  # Reference site for main line
)

# Predict with SE
pred_a <- predict(mod_abundance, newdata = pred_data_a, type = "link", se.fit = TRUE)
pred_data_a$fit <- exp(pred_a$fit)
pred_data_a$lwr <- exp(pred_a$fit - 1.96 * pred_a$se.fit)
pred_data_a$upr <- exp(pred_a$fit + 1.96 * pred_a$se.fit)

# Get scaling exponent from results
abund_row <- scaling_results$all_results[scaling_results$all_results$response == "Total CAFI Abundance", ]
beta_abund_val <- sprintf("%.2f", abund_row$beta)
beta_abund_ci <- sprintf("[%.2f, %.2f]", abund_row$ci_lower, abund_row$ci_upper)

# 1:1 reference line data (Field of Dreams expectation)
ref_line_a <- data.frame(
  volume_field = 10^seq(1.5, 4.5, length.out = 50)
)
# Scale so it passes through median point
median_vol <- median(coral_master$volume_field, na.rm = TRUE)
median_cafi <- median(coral_master$total_cafi, na.rm = TRUE)
ref_line_a$cafi <- median_cafi * (ref_line_a$volume_field / median_vol)^1

p_a <- ggplot(coral_master, aes(x = volume_field, y = total_cafi)) +
  # 1:1 reference line (Field of Dreams expectation)
  geom_line(data = ref_line_a, aes(x = volume_field, y = cafi),
            linetype = "dashed", color = "gray55", linewidth = 0.7) +
  # Confidence ribbon
  geom_ribbon(data = pred_data_a, aes(x = volume_field, y = fit, ymin = lwr, ymax = upr),
              fill = "gray60", alpha = 0.3, inherit.aes = FALSE) +
  # Fitted line
  geom_line(data = pred_data_a, aes(x = volume_field, y = fit),
            color = "gray15", linewidth = 1.1, inherit.aes = FALSE) +
  # Data points with better visibility
  geom_point(aes(fill = site), shape = 21, size = 2.8, alpha = 0.85, stroke = 0.4, color = "white") +
  scale_fill_manual(values = site_colors, name = "Site") +
  scale_x_log10(
    labels = comma,
    breaks = c(100, 1000, 10000),
    limits = c(40, 50000)
  ) +
  scale_y_log10(
    labels = comma,
    limits = c(0.8, 400),
    breaks = c(1, 10, 100)
  ) +
  # Beta annotation - use two separate annotations for clean formatting
  annotate("text", x = 55, y = 280,
           label = paste0("beta == ", beta_abund_val), parse = TRUE,
           size = 4, hjust = 0, fontface = "bold") +
  annotate("text", x = 140, y = 280,
           label = beta_abund_ci,
           size = 3.5, hjust = 0, fontface = "plain") +
  # 1:1 reference label - clearer positioning
  annotate("text", x = 25000, y = 3.5, label = "1:1", color = "gray45",
           size = 3.2, fontface = "italic") +
  labs(
    tag = "A",
    title = "Abundance Scaling",
    subtitle = "Marginally super-linear (beta > 1)",
    x = expression(Coral~Volume~(cm^3)),
    y = "Total CAFI Abundance"
  ) +
  theme_fig2() +
  theme(
    legend.position = c(0.85, 0.18),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7)
  )

# ============================================================================
# PANEL B: Richness ~ Volume
# ============================================================================
cat("Building Panel B: Richness scaling...\n")

# Fit Poisson GLM for richness
mod_richness <- glm(otu_richness ~ log10(volume_field) + site,
                    data = coral_master, family = poisson)

# Create prediction data
pred_data_b <- expand.grid(
  volume_field = 10^seq(log10(min(coral_master$volume_field, na.rm = TRUE)),
                        log10(max(coral_master$volume_field, na.rm = TRUE)),
                        length.out = 100),
  site = "HAU"
)

pred_b <- predict(mod_richness, newdata = pred_data_b, type = "link", se.fit = TRUE)
pred_data_b$fit <- exp(pred_b$fit)
pred_data_b$lwr <- exp(pred_b$fit - 1.96 * pred_b$se.fit)
pred_data_b$upr <- exp(pred_b$fit + 1.96 * pred_b$se.fit)

# Get scaling exponent
rich_row <- scaling_results$all_results[scaling_results$all_results$response == "Species Richness", ]
z_rich_val <- sprintf("%.2f", rich_row$beta)
z_rich_ci <- sprintf("[%.2f, %.2f]", rich_row$ci_lower, rich_row$ci_upper)

# 1:1 reference line
ref_line_b <- data.frame(
  volume_field = 10^seq(1.5, 4.5, length.out = 50)
)
median_rich <- median(coral_master$otu_richness, na.rm = TRUE)
ref_line_b$richness <- median_rich * (ref_line_b$volume_field / median_vol)^1

p_b <- ggplot(coral_master, aes(x = volume_field, y = otu_richness)) +
  # 1:1 reference line
  geom_line(data = ref_line_b, aes(x = volume_field, y = richness),
            linetype = "dashed", color = "gray55", linewidth = 0.7) +
  # Confidence ribbon
  geom_ribbon(data = pred_data_b, aes(x = volume_field, y = fit, ymin = lwr, ymax = upr),
              fill = "gray60", alpha = 0.3, inherit.aes = FALSE) +
  # Fitted line
  geom_line(data = pred_data_b, aes(x = volume_field, y = fit),
            color = "gray15", linewidth = 1.1, inherit.aes = FALSE) +
  # Data points
  geom_point(aes(fill = site), shape = 21, size = 2.8, alpha = 0.85, stroke = 0.4, color = "white") +
  scale_fill_manual(values = site_colors, name = "Site") +
  scale_x_log10(
    labels = comma,
    breaks = c(100, 1000, 10000),
    limits = c(40, 50000)
  ) +
  scale_y_log10(
    limits = c(0.8, 100),
    breaks = c(1, 10, 100)
  ) +
  # z annotation (species-area exponent) - split for clean formatting
  annotate("text", x = 55, y = 70,
           label = paste0("italic(z) == ", z_rich_val), parse = TRUE,
           size = 4, hjust = 0, fontface = "bold") +
  annotate("text", x = 165, y = 70,
           label = z_rich_ci,
           size = 3.5, hjust = 0, fontface = "plain") +
  # 1:1 reference label
  annotate("text", x = 25000, y = 2.8, label = "1:1", color = "gray45",
           size = 3.2, fontface = "italic") +
  labs(
    tag = "B",
    title = "Species-Area Relationship",
    subtitle = "Sublinear (z < 1, Redirection)",
    x = expression(Coral~Volume~(cm^3)),
    y = "Species Richness"
  ) +
  theme_fig2() +
  theme(legend.position = "none")

# Colors for scaling interpretation (used in Panel C)
# Using same palette but with distinct associations for clarity
# ASCII-safe legend labels to avoid encoding issues
interp_colors <- c(
  "Super-linear (beta > 1)" = "#0072B2",      # Blue - aggregation/attraction
  "Field of Dreams (beta ~ 1)" = "#009E73",   # Teal - neutral/proportional
  "Redirection (beta < 1)" = "#D55E00"     # Vermillion - dilution
)

# ============================================================================
# PANEL C: Top 10 Species Scaling (Forest Plot)
# ============================================================================
cat("Building Panel C: Species scaling forest plot...\n")

# Get species-level results and select top 10 by abundance
species_results <- scaling_results$all_results %>%
  filter(category == "Species") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 10) %>%
  mutate(
    # Clean species labels - add "sp." for genus-only names
    species_label = case_when(
      !grepl(" ", response) ~ paste0(response, " sp."),
      TRUE ~ response
    ),
    species_label = factor(species_label, levels = rev(species_label)),
    # Map interpretation to ASCII-safe labels (handle all possible values)
    interpretation_safe = case_when(
      grepl("Super-linear", interpretation) ~ "Super-linear (beta > 1)",
      grepl("Field of Dreams", interpretation) ~ "Field of Dreams (beta ~ 1)",
      grepl("Redirection", interpretation) ~ "Redirection (beta < 1)",
      TRUE ~ "Field of Dreams (beta ~ 1)"  # Default fallback
    )
  )


p_c <- ggplot(species_results, aes(x = beta, y = species_label)) +
  # Shaded hypothesis regions - softer, more professional colors
  annotate("rect", xmin = -0.3, xmax = 0.85, ymin = -Inf, ymax = Inf,
           fill = "#D55E00", alpha = 0.08) +
  annotate("rect", xmin = 0.85, xmax = 1.15, ymin = -Inf, ymax = Inf,
           fill = "#009E73", alpha = 0.10) +
  annotate("rect", xmin = 1.15, xmax = 4.2, ymin = -Inf, ymax = Inf,
           fill = "#0072B2", alpha = 0.08) +
  # Region labels at bottom (use -Inf for discrete y-axis)
  annotate("text", x = 0.27, y = -Inf, label = "Redirection",
           color = "#D55E00", size = 3, fontface = "bold", hjust = 0.5, vjust = 1.5) +
  annotate("text", x = 1, y = -Inf, label = "Field of\nDreams",
           color = "#009E73", size = 2.8, fontface = "bold", hjust = 0.5, vjust = 1.3, lineheight = 0.85) +
  annotate("text", x = 2.7, y = -Inf, label = "Super-linear",
           color = "#0072B2", size = 3, fontface = "bold", hjust = 0.5, vjust = 1.5) +
  # Reference line at beta = 1 (Field of Dreams)
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray35", linewidth = 0.7) +
  # Error bars with slightly thicker lines - use orientation instead of deprecated geom_errorbarh
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper, color = interpretation_safe),
                width = 0.35, linewidth = 0.8, orientation = "y") +
  # Points with better contrast
  geom_point(aes(fill = interpretation_safe), size = 3.5, shape = 21, color = "white", stroke = 1) +
  scale_fill_manual(
    values = interp_colors,
    name = "Scaling Pattern",
    drop = TRUE  # Only show patterns present in data
  ) +
  scale_color_manual(values = interp_colors, guide = "none") +
  scale_x_continuous(
    limits = c(-0.3, 4.2),
    breaks = seq(0, 4, 1),
    expand = c(0, 0)
  ) +
  labs(
    tag = "C",
    title = "Species-Level Scaling (top 10 by abundance)",
    x = expression(paste("Scaling Exponent (", beta, ")")),
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_fig2() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "italic", size = 9.5),
    plot.margin = margin(t = 10, r = 10, b = 35, l = 10)  # Extra bottom margin for region labels
  )

# ============================================================================
# COMBINE PANELS
# ============================================================================
cat("Combining panels...\n")

# Clean layout using patchwork - legend is now inside Panel A
p_a_clean <- p_a + theme(plot.tag.position = "topleft")
p_b_clean <- p_b + theme(plot.tag.position = "topleft")
p_c_clean <- p_c + theme(plot.tag.position = "topleft")

# Layout: A and B on top row, C spanning bottom
fig2_combined <- (p_a_clean | p_b_clean) / p_c_clean +
  plot_layout(heights = c(1, 0.9))

# Build caption with line breaks for proper wrapping
caption_text <- paste0(
  "(A) Total CAFI abundance vs coral volume; beta = ", beta_abund_val, " ", beta_abund_ci,
  " indicates marginally super-linear scaling.\n",
  "(B) Species richness vs volume; z = ", z_rich_val, " ", z_rich_ci,
  " indicates sublinear accumulation (Redirection).\n",
  "(C) Top 10 species by abundance show heterogeneous scaling patterns. ",
  "Dashed line = beta = 1 (Field of Dreams hypothesis). n = 114 corals."
)

# Add overall title, subtitle, and caption
fig2 <- cowplot::ggdraw() +
  cowplot::draw_plot(fig2_combined, x = 0, y = 0.07, width = 1, height = 0.87) +
  cowplot::draw_label(
    "Figure 2: Size-Abundance Scaling in CAFI Communities",
    x = 0.02, y = 0.98, hjust = 0, vjust = 1,
    fontface = "bold", size = 13
  ) +
  cowplot::draw_label(
    "Larger corals host more CAFI, but species accumulate sublinearly; individual species show heterogeneous responses",
    x = 0.02, y = 0.95, hjust = 0, vjust = 1,
    color = "gray30", size = 9
  ) +
  cowplot::draw_label(
    caption_text,
    x = 0.02, y = 0.01, hjust = 0, vjust = 0,
    color = "gray45", size = 7.5, lineheight = 1.2
  )

# Save at publication quality - wider aspect ratio reduces whitespace
ggsave("output/figures/manuscript/fig2_scaling.png", fig2,
       width = 12, height = 10, dpi = 300, bg = "white")

# Also save PDF for journal submission (skip cairo if not available)
tryCatch({
  ggsave("output/figures/manuscript/fig2_scaling.pdf", fig2,
         width = 12, height = 10, device = cairo_pdf, bg = "white")
  cat("  - output/figures/manuscript/fig2_scaling.pdf\n")
}, error = function(e) {
  # Fallback to standard PDF if cairo not available
  ggsave("output/figures/manuscript/fig2_scaling.pdf", fig2,
         width = 12, height = 10, bg = "white")
  cat("  - output/figures/manuscript/fig2_scaling.pdf (standard PDF)\n")
})

cat("\nDone! Figure saved to:\n")
cat("  - output/figures/manuscript/fig2_scaling.png\n")
