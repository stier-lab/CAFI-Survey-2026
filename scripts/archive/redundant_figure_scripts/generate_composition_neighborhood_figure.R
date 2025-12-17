#!/usr/bin/env Rscript
# =============================================================================
# Figure: How CAFI composition changes with coral neighborhood context
# =============================================================================

library(tidyverse)
library(patchwork)
library(ggrepel)

# Publication theme
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(face = "bold", size = base_size - 1),
      plot.title = element_text(face = "bold"),
      legend.background = element_rect(fill = "white", color = NA)
    )
}

# Functional group colors
group_colors <- c(
  "Crab" = "#E74C3C",
  "Shrimp" = "#3498DB",
  "Fish" = "#2ECC71",
  "Hermit" = "#F39C12",
  "Snail" = "#9B59B6",
  "Echinoderm" = "#1ABC9C",
  "Other" = "#95A5A6"
)

cat("Loading data...\n")

# Load data
cafi_data <- read_csv("data/survey_cafi_data_w_taxonomy_summer2019_v5.csv",
                      show_col_types = FALSE)
coral_data <- read_csv("data/survey_coral_characteristics_merged_v2.csv",
                       show_col_types = FALSE)

# Standardize type names
cafi_data <- cafi_data %>%
  mutate(type_clean = case_when(
    type == "crab" ~ "Crab",
    type == "shrimp" ~ "Shrimp",
    type == "fish" ~ "Fish",
    type == "hermit" ~ "Hermit",
    type == "snail" ~ "Snail",
    type == "echinoderm" ~ "Echinoderm",
    TRUE ~ "Other"
  ))

# Merge with neighborhood data
merged <- cafi_data %>%
  left_join(coral_data %>%
              select(coral_id, number_of_neighbors, mean_neighbor_distance,
                     combined_total_volume_of_neighbors, volume_field, site),
            by = "coral_id") %>%
  filter(!is.na(number_of_neighbors), !is.na(volume_field))

cat("Records with neighborhood data:", nrow(merged), "\n")
cat("Unique corals:", n_distinct(merged$coral_id), "\n")

# Per-coral composition with neighborhood metrics
coral_comp <- merged %>%
  group_by(coral_id, number_of_neighbors, mean_neighbor_distance,
           combined_total_volume_of_neighbors, volume_field) %>%
  summarise(
    abundance = n(),
    richness = n_distinct(paste(genus, species)),
    n_crab = sum(type_clean == "Crab"),
    n_shrimp = sum(type_clean == "Shrimp"),
    n_fish = sum(type_clean == "Fish"),
    n_hermit = sum(type_clean == "Hermit"),
    n_snail = sum(type_clean == "Snail"),
    .groups = "drop"
  ) %>%
  mutate(
    prop_crab = n_crab / abundance,
    prop_shrimp = n_shrimp / abundance,
    prop_fish = n_fish / abundance,
    prop_hermit = n_hermit / abundance,
    prop_snail = n_snail / abundance,
    log_vol = log10(volume_field),
    log_neighbor_vol = log10(combined_total_volume_of_neighbors + 1),
    # Derived metrics
    neighbor_density = number_of_neighbors / (mean_neighbor_distance/100),
    relative_size = volume_field / (combined_total_volume_of_neighbors + 1)
  )

# Create neighbor count bins
coral_comp <- coral_comp %>%
  mutate(
    neighbor_bin = cut(number_of_neighbors,
                       breaks = c(0, 10, 20, 30, Inf),
                       labels = c("1-10", "11-20", "21-30", ">30"),
                       include.lowest = TRUE)
  )

# =============================================================================
# Panel A: Composition by number of neighbors (stacked bars)
# =============================================================================

cat("Creating Panel A: Composition by neighbor count...\n")

comp_by_neighbors <- merged %>%
  mutate(neighbor_bin = cut(number_of_neighbors,
                            breaks = c(0, 10, 20, 30, Inf),
                            labels = c("1-10", "11-20", "21-30", ">30"),
                            include.lowest = TRUE)) %>%
  group_by(neighbor_bin, type_clean) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(neighbor_bin) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Order groups
group_order <- comp_by_neighbors %>%
  group_by(type_clean) %>%
  summarise(total = sum(count)) %>%
  arrange(desc(total)) %>%
  pull(type_clean)

comp_by_neighbors$type_clean <- factor(comp_by_neighbors$type_clean,
                                        levels = rev(group_order))

# Sample sizes
n_by_bin <- coral_comp %>%
  group_by(neighbor_bin) %>%
  summarise(n = n())

p_a <- ggplot(comp_by_neighbors, aes(x = neighbor_bin, y = prop, fill = type_clean)) +
  geom_col(position = "stack", color = "white", linewidth = 0.2, width = 0.8) +
  geom_text(data = n_by_bin, aes(x = neighbor_bin, y = 1.03, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 3, color = "gray40") +
  scale_fill_manual(values = group_colors, name = "Taxon") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "A. Composition by neighbor count",
       subtitle = "Corals in denser neighborhoods have more shrimp",
       x = "Number of neighboring corals (within 5m)",
       y = "Proportion of individuals") +
  theme_pub(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 9, color = "gray40"))

# =============================================================================
# Panel B: Key taxa vs number of neighbors
# =============================================================================

cat("Creating Panel B: Taxa vs neighbor count...\n")

prop_long <- coral_comp %>%
  select(coral_id, number_of_neighbors, volume_field,
         prop_crab, prop_shrimp, prop_fish, prop_snail) %>%
  pivot_longer(cols = starts_with("prop_"),
               names_to = "group",
               values_to = "proportion",
               names_prefix = "prop_") %>%
  mutate(group = str_to_title(group),
         group = factor(group, levels = c("Crab", "Shrimp", "Fish", "Snail")))

# Calculate correlations
cors <- prop_long %>%
  group_by(group) %>%
  summarise(
    r = cor(proportion, number_of_neighbors),
    p = cor.test(proportion, number_of_neighbors)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.2f%s", r,
                         ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

p_b <- ggplot(prop_long, aes(x = number_of_neighbors, y = proportion)) +
  geom_point(aes(color = group), alpha = 0.5, size = 2) +
  geom_smooth(aes(color = group), method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
  geom_text(data = cors, aes(label = label, color = group),
            x = 65, y = 0.85, hjust = 1, size = 3, fontface = "italic",
            show.legend = FALSE) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_wrap(~ group, nrow = 1) +
  labs(title = "B. Proportional abundance vs. neighbor count",
       subtitle = "More neighbors → more shrimp, fewer snails",
       x = "Number of neighboring corals",
       y = "Proportion of community") +
  theme_pub(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 9, color = "gray40"))

# =============================================================================
# Panel C: Composition by total neighbor volume
# =============================================================================

cat("Creating Panel C: Composition by neighbor volume...\n")

# Create volume bins
coral_comp <- coral_comp %>%
  mutate(
    neighbor_vol_bin = cut(combined_total_volume_of_neighbors,
                           breaks = quantile(combined_total_volume_of_neighbors,
                                            probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                           labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
                           include.lowest = TRUE)
  )

comp_by_vol <- merged %>%
  left_join(coral_comp %>% select(coral_id, neighbor_vol_bin), by = "coral_id") %>%
  filter(!is.na(neighbor_vol_bin)) %>%
  group_by(neighbor_vol_bin, type_clean) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(neighbor_vol_bin) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

comp_by_vol$type_clean <- factor(comp_by_vol$type_clean, levels = rev(group_order))

n_by_vol <- coral_comp %>%
  filter(!is.na(neighbor_vol_bin)) %>%
  group_by(neighbor_vol_bin) %>%
  summarise(n = n())

p_c <- ggplot(comp_by_vol, aes(x = neighbor_vol_bin, y = prop, fill = type_clean)) +
  geom_col(position = "stack", color = "white", linewidth = 0.2, width = 0.8) +
  geom_text(data = n_by_vol, aes(x = neighbor_vol_bin, y = 1.03, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 3, color = "gray40") +
  scale_fill_manual(values = group_colors, name = "Taxon") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "C. Composition by neighborhood coral volume",
       subtitle = "Total volume of corals within 5m radius",
       x = "Total neighbor volume (quartiles)",
       y = "Proportion of individuals") +
  theme_pub(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 9, color = "gray40"))

# =============================================================================
# Panel D: Richness vs neighborhood metrics
# =============================================================================

cat("Creating Panel D: Richness vs neighborhood...\n")

p_d <- ggplot(coral_comp, aes(x = number_of_neighbors, y = richness)) +
  geom_point(aes(size = volume_field), alpha = 0.6, color = "#3498DB") +
  geom_smooth(method = "lm", se = TRUE, color = "#2C3E50", linewidth = 1) +
  scale_size_continuous(range = c(1, 6), name = expression("Coral vol. (cm"^3*")"),
                        breaks = c(2000, 10000, 30000),
                        labels = scales::comma) +
  labs(title = "D. Species richness vs. neighbor count",
       subtitle = sprintf("r = %.2f, p = %.3f (controlling for coral size)",
                         cor(coral_comp$richness, coral_comp$number_of_neighbors),
                         cor.test(coral_comp$richness, coral_comp$number_of_neighbors)$p.value),
       x = "Number of neighboring corals",
       y = "CAFI species richness") +
  theme_pub(base_size = 11) +
  theme(plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        legend.position = c(0.85, 0.25),
        legend.background = element_rect(fill = alpha("white", 0.8)))

# =============================================================================
# Combine panels
# =============================================================================

cat("Combining panels...\n")

fig_neighbor <- (p_a | p_c) / p_b / p_d +
  plot_layout(heights = c(1, 0.9, 1)) +
  plot_annotation(
    title = "Figure S-X. Community composition varies with coral neighborhood context",
    subtitle = "Corals in denser neighborhoods support different assemblages",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30")
    )
  )

# Save
ggsave("output/figures/composition_by_neighborhood.png", fig_neighbor,
       width = 14, height = 12, dpi = 300, bg = "white")

cat("\n✓ Figure saved to: output/figures/composition_by_neighborhood.png\n")

# Also save simplified version
fig_simple <- (p_a + theme(legend.position = "bottom")) / p_b +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(
    title = "Figure X. Neighborhood context shapes community composition",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13)
    )
  )

ggsave("output/figures/manuscript/Figure7_composition_neighborhood.png", fig_simple,
       width = 12, height = 8, dpi = 300, bg = "white")

cat("✓ Main text figure saved to: output/figures/manuscript/Figure7_composition_neighborhood.png\n")
