#!/usr/bin/env Rscript
# =============================================================================
# Figure: How CAFI composition changes with coral size
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
  "Worm" = "#E67E22",
  "Amphipod" = "#95A5A6",
  "Other" = "#BDC3C7"
)

cat("Loading data...\n")

# Load data
cafi_data <- read_csv("data/survey_cafi_data_w_taxonomy_summer2019_v5.csv", show_col_types = FALSE)
coral_data <- read_csv("data/survey_coral_characteristics_merged_v2.csv", show_col_types = FALSE)

# Standardize type names
cafi_data <- cafi_data %>%
  mutate(type_clean = case_when(
    type == "crab" ~ "Crab",
    type == "shrimp" ~ "Shrimp",
    type == "fish" ~ "Fish",
    type == "hermit" ~ "Hermit",
    type == "snail" ~ "Snail",
    type == "echinoderm" ~ "Echinoderm",
    type == "worm" ~ "Worm",
    type %in% c("amphipod", "amph") ~ "Amphipod",
    type == "squat_lobster" ~ "Other",
    TRUE ~ "Other"
  ))

# Merge with coral size
merged <- cafi_data %>%
  left_join(coral_data %>% select(coral_id, volume_field), by = "coral_id") %>%
  filter(!is.na(volume_field)) %>%
  mutate(log_vol = log10(volume_field))

# Create size bins for stacked bar chart
merged <- merged %>%
  mutate(
    size_bin = cut(volume_field,
                   breaks = c(0, 1000, 3000, 7000, 15000, Inf),
                   labels = c("<1000", "1000-3000", "3000-7000", "7000-15000", ">15000")),
    size_bin = factor(size_bin, levels = c("<1000", "1000-3000", "3000-7000", "7000-15000", ">15000"))
  )

# Per-coral composition
coral_composition <- merged %>%
  group_by(coral_id, volume_field, log_vol) %>%
  summarise(
    abundance = n(),
    richness = n_distinct(paste(genus, species)),
    n_crab = sum(type_clean == "Crab"),
    n_shrimp = sum(type_clean == "Shrimp"),
    n_fish = sum(type_clean == "Fish"),
    n_hermit = sum(type_clean == "Hermit"),
    n_snail = sum(type_clean == "Snail"),
    n_echinoderm = sum(type_clean == "Echinoderm"),
    .groups = "drop"
  ) %>%
  mutate(
    prop_crab = n_crab / abundance,
    prop_shrimp = n_shrimp / abundance,
    prop_fish = n_fish / abundance,
    prop_hermit = n_hermit / abundance,
    prop_snail = n_snail / abundance,
    prop_echinoderm = n_echinoderm / abundance
  )

# =============================================================================
# Panel A: Stacked area chart showing composition across size gradient
# =============================================================================

cat("Creating Panel A: Stacked composition by size...\n")

# Aggregate by size bin
comp_by_bin <- merged %>%
  group_by(size_bin, type_clean) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(size_bin) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Order groups by overall abundance
group_order <- comp_by_bin %>%
  group_by(type_clean) %>%
  summarise(total = sum(count)) %>%
  arrange(desc(total)) %>%
  pull(type_clean)

comp_by_bin$type_clean <- factor(comp_by_bin$type_clean, levels = rev(group_order))

# Get sample sizes for labels
sample_sizes <- merged %>%
  group_by(size_bin) %>%
  summarise(n_corals = n_distinct(coral_id), n_individuals = n())

p_a <- ggplot(comp_by_bin, aes(x = size_bin, y = prop, fill = type_clean)) +
  geom_col(position = "stack", color = "white", linewidth = 0.2, width = 0.85) +
  geom_text(data = sample_sizes,
            aes(x = size_bin, y = 1.02, label = paste0("n=", n_corals)),
            inherit.aes = FALSE, size = 2.8, color = "gray40") +
  scale_fill_manual(values = group_colors, name = "Taxon") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(title = "A. Community composition shifts with coral size",
       x = expression("Coral volume (cm"^3*")"),
       y = "Proportion of individuals") +
  theme_pub(base_size = 11) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(size = 12))

# =============================================================================
# Panel B: Scatter plots for key groups with trend lines
# =============================================================================

cat("Creating Panel B: Proportional trends...\n")

# Prepare long format for faceted plot
prop_long <- coral_composition %>%
  select(coral_id, volume_field, log_vol, prop_crab, prop_shrimp, prop_fish, prop_snail) %>%
  pivot_longer(cols = starts_with("prop_"),
               names_to = "group",
               values_to = "proportion",
               names_prefix = "prop_") %>%
  mutate(group = str_to_title(group),
         group = factor(group, levels = c("Crab", "Shrimp", "Fish", "Snail")))

# Calculate correlations for labels
cors <- prop_long %>%
  group_by(group) %>%
  summarise(
    r = cor(proportion, log_vol),
    p = cor.test(proportion, log_vol)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.2f%s", r, ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))))

p_b <- ggplot(prop_long, aes(x = volume_field, y = proportion)) +
  geom_point(aes(color = group), alpha = 0.5, size = 1.5) +
  geom_smooth(aes(color = group), method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
  geom_text(data = cors, aes(label = label, color = group),
            x = 35000, y = 0.85, hjust = 1, size = 3, fontface = "italic",
            show.legend = FALSE) +
  scale_x_log10(labels = scales::comma,
                breaks = c(100, 500, 2000, 10000, 40000)) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_wrap(~ group, nrow = 1) +
  labs(title = "B. Proportional abundance vs. coral size",
       subtitle = "Crabs decline while fish and snails increase in larger corals",
       x = expression("Coral volume (cm"^3*", log scale)"),
       y = "Proportion of community") +
  theme_pub(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 9, color = "gray40"))

# =============================================================================
# Panel C: Species-specific size preferences
# =============================================================================

cat("Creating Panel C: Species size preferences...\n")

# Calculate median volume for common species
species_sizes <- merged %>%
  mutate(species_name = paste(genus, species)) %>%
  group_by(species_name, type_clean) %>%
  summarise(
    n = n(),
    median_vol = median(volume_field),
    mean_vol = mean(volume_field),
    q25_vol = quantile(volume_field, 0.25),
    q75_vol = quantile(volume_field, 0.75),
    .groups = "drop"
  ) %>%
  filter(n >= 20) %>%
  arrange(median_vol) %>%
  # Clean up species names
  mutate(
    species_label = str_replace(species_name, "^([A-Z][a-z]+) \\1 ", "\\1 "),
    species_label = str_replace(species_label, " NA$", " sp."),
    species_label = str_replace(species_label, "NA NA", "Unknown"),
    species_label = str_replace(species_label, "#N/A #N/A", "Unknown fish")
  )

# Top and bottom species by median size
top_species <- bind_rows(
  species_sizes %>% slice_min(median_vol, n = 8),
  species_sizes %>% slice_max(median_vol, n = 8)
) %>%
  distinct(species_name, .keep_all = TRUE) %>%
  mutate(species_label = fct_reorder(species_label, median_vol))

p_c <- ggplot(top_species, aes(y = species_label, x = median_vol, color = type_clean)) +
  geom_errorbarh(aes(xmin = q25_vol, xmax = q75_vol), height = 0.3, linewidth = 0.8, alpha = 0.6) +
  geom_point(aes(size = n), alpha = 0.9) +
  geom_vline(xintercept = median(merged$volume_field), linetype = "dashed", color = "gray50") +
  annotate("text", x = median(merged$volume_field), y = 0.5,
           label = "Median coral size", hjust = -0.1, vjust = -0.5,
           size = 2.8, color = "gray40", fontface = "italic") +
  scale_x_log10(labels = scales::comma, breaks = c(500, 2000, 5000, 15000, 40000)) +
  scale_color_manual(values = group_colors, name = "Taxon") +
  scale_size_continuous(range = c(2, 6), name = "n individuals", breaks = c(25, 50, 100, 200)) +
  labs(title = "C. Species differ in coral size preferences",
       subtitle = "Bars show interquartile range; vertical line = median coral size",
       x = expression("Coral volume (cm"^3*", log scale)"),
       y = NULL) +
  theme_pub(base_size = 11) +
  theme(axis.text.y = element_text(face = "italic", size = 9),
        plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        legend.position = "right",
        legend.box = "vertical")

# =============================================================================
# Panel D: Beta diversity / turnover with size
# =============================================================================

cat("Creating Panel D: Community turnover...\n")

# Create size classes for beta diversity
coral_for_beta <- coral_composition %>%
  mutate(size_class = ntile(volume_field, 5)) %>%
  arrange(size_class)

# For each size class, calculate mean pairwise Bray-Curtis with other classes
# Simplified: show proportion of unique species in each size class

species_by_class <- merged %>%
  mutate(size_class = ntile(volume_field, 5),
         species_name = paste(genus, species)) %>%
  group_by(size_class) %>%
  summarise(
    n_species = n_distinct(species_name),
    species_list = list(unique(species_name)),
    .groups = "drop"
  )

# Find unique species per class
all_species <- unique(merged$species_name <- paste(merged$genus, merged$species))

unique_to_class <- sapply(1:5, function(i) {
  class_species <- species_by_class$species_list[[i]]
  other_species <- unique(unlist(species_by_class$species_list[-i]))
  sum(!class_species %in% other_species)
})

species_by_class$unique_species <- unique_to_class
species_by_class$prop_unique <- unique_to_class / species_by_class$n_species

# Size labels
size_labels <- coral_composition %>%
  mutate(size_class = ntile(volume_field, 5)) %>%
  group_by(size_class) %>%
  summarise(
    min_vol = min(volume_field),
    max_vol = max(volume_field),
    mean_vol = mean(volume_field),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%d–%dk", round(min_vol/1000), round(max_vol/1000)))

species_by_class <- species_by_class %>%
  left_join(size_labels, by = "size_class")

p_d <- ggplot(species_by_class, aes(x = factor(size_class), y = n_species)) +
  geom_col(fill = "#3498DB", alpha = 0.7, width = 0.7) +
  geom_text(aes(label = n_species), vjust = -0.3, size = 3.5, fontface = "bold") +
  geom_text(aes(y = 5, label = label), size = 2.8, color = "white") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "D. Species richness by coral size quintile",
       subtitle = expression("Volume ranges (cm"^3*") shown in bars"),
       x = "Size quintile (1 = smallest, 5 = largest)",
       y = "Number of species") +
  theme_pub(base_size = 11) +
  theme(plot.title = element_text(size = 12),
        plot.subtitle = element_text(size = 9, color = "gray40"))

# =============================================================================
# Combine panels
# =============================================================================

cat("Combining panels...\n")

fig_comp_size <- (p_a | p_d) / p_b / p_c +
  plot_layout(heights = c(1, 0.8, 1.2)) +
  plot_annotation(
    title = "Figure S-X. Community composition changes with coral colony size",
    subtitle = "Larger corals support different species assemblages than smaller corals",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray30")
    )
  )

# Save
ggsave("output/figures/composition_by_size.png", fig_comp_size,
       width = 14, height = 12, dpi = 300, bg = "white")

cat("\n✓ Figure saved to: output/figures/composition_by_size.png\n")

# Also save a simplified 2-panel version for main text
fig_simple <- (p_a + theme(legend.position = "bottom", legend.box = "horizontal")) /
  (p_b + theme(strip.text = element_text(size = 10))) +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(
    title = "Figure X. Community composition shifts with coral size",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13)
    )
  )

ggsave("output/figures/manuscript/Figure6_composition_size.png", fig_simple,
       width = 12, height = 8, dpi = 300, bg = "white")

cat("✓ Main text figure saved to: output/figures/manuscript/Figure6_composition_size.png\n")
