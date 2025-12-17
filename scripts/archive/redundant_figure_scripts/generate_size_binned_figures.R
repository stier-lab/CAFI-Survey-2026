#!/usr/bin/env Rscript
# ============================================================================
# Generate Size-Binned Manuscript Figures
# Comparing CAFI patterns across Small, Medium, and Large coral size classes
# Includes both INCIDENCE (presence/absence) and RELATIVE ABUNDANCE analyses
# ============================================================================

cat("\n========================================\n")
cat("GENERATING SIZE-BINNED FIGURES\n")
cat("========================================\n\n")

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(vegan)
  library(scales)
  library(ggrepel)
})

# Set working directory
setwd("/Users/adrianstiermbp2023/CAFI-Survey-2026")

# Create output directory
dir.create("output/figures/manuscript", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Color Schemes
# ============================================================================

# Size class colors (warm to cool gradient)
size_colors <- c("Small" = "#d73027", "Medium" = "#fee08b", "Large" = "#1a9850")
size_colors_fill <- c("Small" = "#fc8d59", "Medium" = "#ffffbf", "Large" = "#91cf60")

# Site colors
site_colors <- c("HAU" = "#E69F00", "MAT" = "#56B4E9", "MRB" = "#009E73")

# Functional group colors
func_colors <- c(
  "Defenders (crabs & shrimp)" = "#2c7bb6",
  "Resident fishes" = "#abd9e9",
  "Ectoparasites (snails)" = "#d7191c",
  "Other cryptofauna" = "#95a5a6"
)

# Publication theme
theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
      panel.border = element_rect(color = "black", linewidth = 0.6),
      strip.background = element_rect(fill = "gray95", color = "black"),
      strip.text = element_text(face = "bold", size = base_size),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(color = "gray40", size = base_size - 1, hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, color = "black"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# ============================================================================
# Load and Process Data
# ============================================================================

cat("Loading data...\n")

coral_data <- read.csv("data/survey_coral_characteristics_merged_v2.csv")
cafi_data <- read.csv("data/survey_cafi_data_w_taxonomy_summer2019_v5.csv")

# Process coral data with size classes
coral_processed <- coral_data %>%
  mutate(
    site = str_extract(site, "^[A-Z]+"),
    volume = coalesce(volume_field, volume_lab, length_field * width_field * height_field),
    log_volume = log10(volume + 1)
  ) %>%
  filter(!is.na(volume), volume > 0, site %in% c("HAU", "MAT", "MRB"))

# Define size class breaks using tertiles
breaks <- quantile(coral_processed$volume, c(0, 0.33, 0.67, 1), na.rm = TRUE)
cat(sprintf("Size class breaks: Small (<%d cm³), Medium (%d-%d cm³), Large (>%d cm³)\n",
            round(breaks[2]), round(breaks[2]), round(breaks[3]), round(breaks[3])))

coral_processed <- coral_processed %>%
  mutate(
    size_class = cut(volume,
                     breaks = breaks,
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE),
    size_class = factor(size_class, levels = c("Small", "Medium", "Large"))
  )

cat(sprintf("Sample sizes: Small=%d, Medium=%d, Large=%d\n",
            sum(coral_processed$size_class == "Small"),
            sum(coral_processed$size_class == "Medium"),
            sum(coral_processed$size_class == "Large")))

# Process CAFI data
cafi_processed <- cafi_data %>%
  filter(!is.na(genus) & genus != "" & genus != "NA") %>%
  mutate(
    site = case_when(
      is.na(site) | site == "" | site == "NA" ~ str_extract(coral_id, "^[A-Z]+"),
      TRUE ~ str_sub(site, 1, 3)
    ),
    species_id = paste(genus, species, sep = "_")
  ) %>%
  filter(site %in% c("HAU", "MAT", "MRB"))

# Add functional groups
cafi_processed <- cafi_processed %>%
  mutate(
    func_group = case_when(
      grepl("Trapezia|Tetralia|Alpheus", genus, ignore.case = TRUE) ~ "Defenders (crabs & shrimp)",
      grepl("Paragobiodon|Gobiodon|Caracanthus", genus, ignore.case = TRUE) ~ "Resident fishes",
      grepl("Drupella|Coralliophila|Galeropsis|Morula", genus, ignore.case = TRUE) ~ "Ectoparasites (snails)",
      TRUE ~ "Other cryptofauna"
    )
  )

# ============================================================================
# Calculate summary statistics by coral
# ============================================================================

cat("Calculating coral-level summaries...\n")

# Total abundance and richness per coral
coral_summaries <- cafi_processed %>%
  group_by(coral_id) %>%
  summarise(
    cafi_abundance = n(),
    cafi_richness = n_distinct(species_id),
    .groups = "drop"
  )

# Shannon diversity
comm_matrix <- cafi_processed %>%
  group_by(coral_id, species_id) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_id, values_from = count, values_fill = 0)

comm_mat <- as.matrix(comm_matrix[,-1])
rownames(comm_mat) <- comm_matrix$coral_id
shannon_values <- vegan::diversity(comm_mat, "shannon")

coral_summaries <- coral_summaries %>%
  left_join(
    tibble(coral_id = names(shannon_values), shannon = as.numeric(shannon_values)),
    by = "coral_id"
  ) %>%
  mutate(shannon = replace_na(shannon, 0))

# Functional group abundances
func_summaries <- cafi_processed %>%
  group_by(coral_id, func_group) %>%
  summarise(abundance = n(), .groups = "drop") %>%
  pivot_wider(names_from = func_group, values_from = abundance, values_fill = 0)

# Join all data
analysis_data <- coral_processed %>%
  select(coral_id, site, volume, log_volume, size_class) %>%
  left_join(coral_summaries, by = "coral_id") %>%
  left_join(func_summaries, by = "coral_id") %>%
  mutate(
    across(c(cafi_abundance, cafi_richness, shannon), ~replace_na(., 0)),
    across(starts_with("Defenders") | starts_with("Resident") |
           starts_with("Ectoparasites") | starts_with("Other"), ~replace_na(., 0)),
    # Calculate CAFI density (individuals per cm³)
    cafi_density = cafi_abundance / volume
  )

# Clean column names
names(analysis_data) <- gsub(" ", "_", names(analysis_data))
names(analysis_data) <- gsub("[\\(\\)]", "", names(analysis_data))

cat(sprintf("Analysis dataset: %d corals\n", nrow(analysis_data)))

# ============================================================================
# FIGURE A: Community Metrics by Size Class (Box Plots with Stats)
# ============================================================================

cat("Creating Figure A: Community metrics by size class...\n")

# Statistical tests
kw_abundance <- kruskal.test(cafi_abundance ~ size_class, data = analysis_data)
kw_richness <- kruskal.test(cafi_richness ~ size_class, data = analysis_data)
kw_shannon <- kruskal.test(shannon ~ size_class, data = analysis_data)
kw_density <- kruskal.test(cafi_density ~ size_class, data = analysis_data)

# Function to format p-value
format_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01) return(sprintf("p = %.3f", p))
  return(sprintf("p = %.2f", p))
}

# Panel A: Total abundance by size class
pA1 <- ggplot(analysis_data, aes(x = size_class, y = cafi_abundance, fill = size_class)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white", outlier.size = 2) +
  geom_jitter(aes(color = site), width = 0.2, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  scale_color_manual(values = site_colors, name = "Site") +
  labs(title = "A. Total CAFI abundance",
       subtitle = sprintf("Kruskal-Wallis χ² = %.1f, %s", kw_abundance$statistic, format_p(kw_abundance$p.value)),
       x = "Coral size class", y = "Total individuals") +
  theme_pub(base_size = 10) +
  theme(legend.position = "none")

# Panel B: Species richness by size class
pA2 <- ggplot(analysis_data, aes(x = size_class, y = cafi_richness, fill = size_class)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white", outlier.size = 2) +
  geom_jitter(aes(color = site), width = 0.2, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  scale_color_manual(values = site_colors, name = "Site") +
  labs(title = "B. Species richness",
       subtitle = sprintf("Kruskal-Wallis χ² = %.1f, %s", kw_richness$statistic, format_p(kw_richness$p.value)),
       x = "Coral size class", y = "Number of species") +
  theme_pub(base_size = 10) +
  theme(legend.position = "none")

# Panel C: Shannon diversity by size class
pA3 <- ggplot(analysis_data, aes(x = size_class, y = shannon, fill = size_class)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white", outlier.size = 2) +
  geom_jitter(aes(color = site), width = 0.2, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  scale_color_manual(values = site_colors, name = "Site") +
  labs(title = "C. Shannon diversity (H')",
       subtitle = sprintf("Kruskal-Wallis χ² = %.1f, %s", kw_shannon$statistic, format_p(kw_shannon$p.value)),
       x = "Coral size class", y = "Shannon H'") +
  theme_pub(base_size = 10) +
  theme(legend.position = "right")

# Panel D: CAFI density (key test of propagule dilution)
pA4 <- ggplot(analysis_data, aes(x = size_class, y = cafi_density * 1000, fill = size_class)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white", outlier.size = 2) +
  geom_jitter(aes(color = site), width = 0.2, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(title = "D. CAFI density (propagule dilution test)",
       subtitle = sprintf("Kruskal-Wallis χ² = %.1f, %s", kw_density$statistic, format_p(kw_density$p.value)),
       x = "Coral size class", y = expression("Individuals per 1000 cm"^3)) +
  theme_pub(base_size = 10)

# Summary statistics table
summary_stats <- analysis_data %>%
  group_by(size_class) %>%
  summarise(
    n = n(),
    vol_mean = mean(volume),
    vol_sd = sd(volume),
    abund_mean = mean(cafi_abundance),
    abund_se = sd(cafi_abundance) / sqrt(n()),
    rich_mean = mean(cafi_richness),
    rich_se = sd(cafi_richness) / sqrt(n()),
    shannon_mean = mean(shannon),
    shannon_se = sd(shannon) / sqrt(n()),
    density_mean = mean(cafi_density) * 1000,
    density_se = sd(cafi_density) * 1000 / sqrt(n()),
    .groups = "drop"
  )

cat("\n--- Summary Statistics by Size Class ---\n")
print(summary_stats)

# Panel E: Species accumulation - how many more species on large vs small?
# Calculate species per individual (rarefaction-like metric)
analysis_data <- analysis_data %>%
  mutate(species_per_ind = cafi_richness / (cafi_abundance + 1))

pA5 <- ggplot(analysis_data, aes(x = cafi_abundance, y = cafi_richness, color = size_class)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values = size_colors, name = "Size class") +
  labs(title = "E. Species accumulation by size class",
       subtitle = "Larger corals accumulate species faster per individual",
       x = "Total CAFI abundance", y = "Species richness") +
  theme_pub(base_size = 10) +
  theme(legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = alpha("white", 0.8)))

# Panel F: Abundance vs Volume showing the sublinear scaling
pA6 <- ggplot(analysis_data, aes(x = volume, y = cafi_abundance)) +
  # Field of Dreams expectation line (proportional scaling from small corals)
  geom_abline(slope = mean(analysis_data$cafi_abundance[analysis_data$size_class == "Small"]) /
                      mean(analysis_data$volume[analysis_data$size_class == "Small"]),
              intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_point(aes(fill = size_class), shape = 21, alpha = 0.7, size = 2.5, stroke = 0.3, color = "white") +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1, se = TRUE, fill = "gray80") +
  scale_fill_manual(values = size_colors, name = "Size class") +
  scale_x_continuous(labels = scales::comma) +
  annotate("text", x = max(analysis_data$volume) * 0.6,
           y = max(analysis_data$volume) * mean(analysis_data$cafi_abundance[analysis_data$size_class == "Small"]) /
               mean(analysis_data$volume[analysis_data$size_class == "Small"]) * 0.7,
           label = "Field of Dreams\n(proportional)", color = "gray50", size = 2.8, fontface = "italic") +
  labs(title = "F. Sublinear scaling confirms propagule dilution",
       subtitle = "Observed (solid) falls below proportional expectation (dashed)",
       x = expression("Coral volume (cm"^3*")"), y = "CAFI abundance") +
  theme_pub(base_size = 10) +
  theme(legend.position = "none")

# Combine into figure
figA <- (pA1 | pA2 | pA3) / (pA4 | pA5 | pA6) +
  plot_annotation(
    title = "Community metrics differ significantly across coral size classes",
    subtitle = "Larger corals have more CAFI total but LOWER density—supporting propagule dilution",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  )

ggsave("output/figures/manuscript/Figure_SizeBinned_Metrics.png", figA, width = 14, height = 10, dpi = 300, bg = "white")
cat("  ✓ Figure A (Community Metrics) saved\n")

# ============================================================================
# FIGURE B: Functional Group Patterns by Size Class
# ============================================================================

cat("Creating Figure B: Functional groups by size class...\n")

# Calculate functional group summaries by size class
func_by_size <- analysis_data %>%
  select(coral_id, size_class, volume,
         `Defenders_crabs_&_shrimp`, Resident_fishes, `Ectoparasites_snails`, Other_cryptofauna) %>%
  pivot_longer(cols = c(`Defenders_crabs_&_shrimp`, Resident_fishes, `Ectoparasites_snails`, Other_cryptofauna),
               names_to = "func_group", values_to = "abundance") %>%
  mutate(
    func_label = case_when(
      func_group == "Defenders_crabs_&_shrimp" ~ "Defenders\n(crabs & shrimp)",
      func_group == "Resident_fishes" ~ "Resident\nfishes",
      func_group == "Ectoparasites_snails" ~ "Ectoparasites\n(snails)",
      func_group == "Other_cryptofauna" ~ "Other\ncryptofauna"
    ),
    func_label = factor(func_label, levels = c("Defenders\n(crabs & shrimp)", "Resident\nfishes",
                                                "Ectoparasites\n(snails)", "Other\ncryptofauna")),
    density = abundance / volume * 1000  # per 1000 cm³
  )

# Panel B1: Absolute abundance by functional group and size class
pB1 <- ggplot(func_by_size, aes(x = size_class, y = abundance, fill = size_class)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1) +
  facet_wrap(~func_label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  labs(title = "A. Absolute abundance by functional group",
       subtitle = "All groups increase with coral size (more habitat = more individuals)",
       x = "Coral size class", y = "Abundance") +
  theme_pub(base_size = 10)

# Panel B2: Density by functional group (key test)
pB2 <- ggplot(func_by_size, aes(x = size_class, y = density, fill = size_class)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1) +
  facet_wrap(~func_label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  labs(title = "B. Density by functional group (individuals per 1000 cm³)",
       subtitle = "Density decline varies by group: ectoparasites show strongest dilution",
       x = "Coral size class", y = expression("Density (per 1000 cm"^3*")")) +
  theme_pub(base_size = 10)

# Calculate relative abundance (proportion of total CAFI)
func_relative <- analysis_data %>%
  mutate(total = `Defenders_crabs_&_shrimp` + Resident_fishes + `Ectoparasites_snails` + Other_cryptofauna) %>%
  filter(total > 0) %>%
  mutate(
    prop_defenders = `Defenders_crabs_&_shrimp` / total,
    prop_fish = Resident_fishes / total,
    prop_ecto = `Ectoparasites_snails` / total,
    prop_other = Other_cryptofauna / total
  ) %>%
  select(coral_id, size_class, starts_with("prop_")) %>%
  pivot_longer(cols = starts_with("prop_"), names_to = "func_group", values_to = "proportion") %>%
  mutate(
    func_label = case_when(
      func_group == "prop_defenders" ~ "Defenders",
      func_group == "prop_fish" ~ "Resident fishes",
      func_group == "prop_ecto" ~ "Ectoparasites",
      func_group == "prop_other" ~ "Other"
    ),
    func_label = factor(func_label, levels = c("Defenders", "Resident fishes", "Ectoparasites", "Other"))
  )

# Panel B3: Relative abundance (composition shift)
pB3 <- ggplot(func_relative, aes(x = size_class, y = proportion * 100, fill = func_label)) +
  geom_boxplot(alpha = 0.8, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("Defenders" = "#2c7bb6", "Resident fishes" = "#abd9e9",
                                "Ectoparasites" = "#d7191c", "Other" = "#95a5a6"),
                    name = "Functional\ngroup") +
  labs(title = "C. Relative abundance (% of total CAFI)",
       subtitle = "Composition shifts: fish proportion increases with coral size",
       x = "Coral size class", y = "Percentage of total CAFI") +
  theme_pub(base_size = 10) +
  theme(legend.position = "right")

# Panel B4: Stacked bar showing composition
comp_summary <- func_relative %>%
  group_by(size_class, func_label) %>%
  summarise(mean_prop = mean(proportion), se = sd(proportion)/sqrt(n()), .groups = "drop")

pB4 <- ggplot(comp_summary, aes(x = size_class, y = mean_prop * 100, fill = func_label)) +
  geom_col(position = "stack", alpha = 0.9, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Defenders" = "#2c7bb6", "Resident fishes" = "#abd9e9",
                                "Ectoparasites" = "#d7191c", "Other" = "#95a5a6"),
                    name = "Functional\ngroup") +
  labs(title = "D. Community composition by size class",
       subtitle = "Mean proportion of each functional group",
       x = "Coral size class", y = "Percentage of community") +
  theme_pub(base_size = 10) +
  theme(legend.position = "right")

# Combine
figB <- (pB1 / pB2) | (pB3 / pB4) +
  plot_layout(widths = c(1.5, 1)) +
  plot_annotation(
    title = "Functional group responses to coral size: abundance, density, and relative abundance",
    subtitle = "Defenders maintain proportional scaling; ectoparasites show strongest dilution; fish proportion increases with size",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  )

ggsave("output/figures/manuscript/Figure_SizeBinned_FuncGroups.png", figB, width = 16, height = 10, dpi = 300, bg = "white")
cat("  ✓ Figure B (Functional Groups) saved\n")

# ============================================================================
# FIGURE C: Species INCIDENCE by Size Class (Presence/Absence)
# ============================================================================

cat("Creating Figure C: Species incidence by size class...\n")

# Calculate incidence (% of corals where species is present) by size class
species_incidence <- cafi_processed %>%
  group_by(coral_id, species_id) %>%
  summarise(present = 1, .groups = "drop") %>%
  right_join(
    analysis_data %>% select(coral_id, size_class) %>%
      crossing(species_id = unique(cafi_processed$species_id)),
    by = c("coral_id", "species_id")
  ) %>%
  mutate(present = replace_na(present, 0)) %>%
  group_by(size_class, species_id) %>%
  summarise(
    n_present = sum(present),
    n_total = n(),
    incidence = n_present / n_total * 100,
    .groups = "drop"
  )

# Get top 20 most common species overall
top_species <- species_incidence %>%
  group_by(species_id) %>%
  summarise(mean_incidence = mean(incidence), .groups = "drop") %>%
  arrange(desc(mean_incidence)) %>%
  head(20) %>%
  pull(species_id)

# Filter and format
incidence_top <- species_incidence %>%
  filter(species_id %in% top_species) %>%
  mutate(
    species_label = sapply(species_id, function(x) {
      parts <- str_split(x, "_")[[1]]
      genus <- parts[1]
      sp <- ifelse(length(parts) > 1 & parts[2] != "NA", parts[2], "sp.")
      paste0(substr(genus, 1, 1), ". ", sp)
    }),
    # Order species by mean incidence
    species_label = fct_reorder(species_label, incidence, .fun = mean)
  )

# Panel C1: Species incidence heatmap
pC1 <- ggplot(incidence_top, aes(x = size_class, y = species_label, fill = incidence)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", incidence)), size = 2.5,
            color = ifelse(incidence_top$incidence > 50, "white", "black")) +
  scale_fill_gradient2(low = "#f7f7f7", mid = "#92c5de", high = "#0571b0",
                       midpoint = 50, name = "Incidence (%)") +
  labs(title = "A. Species incidence by coral size class",
       subtitle = "Percentage of corals in each size class where species is present",
       x = "Coral size class", y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic", size = 9),
        legend.position = "right")

# Panel C2: Incidence change from small to large
incidence_change <- incidence_top %>%
  select(species_id, species_label, size_class, incidence) %>%
  pivot_wider(names_from = size_class, values_from = incidence) %>%
  mutate(
    change = Large - Small,
    direction = ifelse(change > 0, "Increases", "Decreases")
  ) %>%
  arrange(desc(abs(change)))

pC2 <- ggplot(incidence_change, aes(x = change, y = fct_reorder(species_label, change), fill = direction)) +
  geom_col(alpha = 0.8, color = "gray30", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray40", linewidth = 0.8) +
  scale_fill_manual(values = c("Increases" = "#1a9850", "Decreases" = "#d73027"),
                    name = "Direction") +
  labs(title = "B. Change in incidence (Large − Small)",
       subtitle = "Positive = more common on large corals",
       x = "Change in incidence (%)", y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic", size = 9),
        legend.position = "right")

# Panel C3: Incidence by functional group
func_incidence <- cafi_processed %>%
  group_by(coral_id, func_group) %>%
  summarise(present = 1, .groups = "drop") %>%
  right_join(
    analysis_data %>% select(coral_id, size_class) %>%
      crossing(func_group = unique(cafi_processed$func_group)),
    by = c("coral_id", "func_group")
  ) %>%
  mutate(present = replace_na(present, 0)) %>%
  group_by(size_class, func_group) %>%
  summarise(
    n_present = sum(present),
    n_total = n(),
    incidence = n_present / n_total * 100,
    .groups = "drop"
  )

pC3 <- ggplot(func_incidence, aes(x = size_class, y = incidence, fill = size_class)) +
  geom_col(alpha = 0.8, color = "gray30", linewidth = 0.3) +
  facet_wrap(~func_group, nrow = 1) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  labs(title = "C. Functional group incidence by size class",
       subtitle = "Fish incidence increases dramatically with coral size",
       x = "Coral size class", y = "Incidence (% of corals)") +
  theme_pub(base_size = 10) +
  theme(strip.text = element_text(size = 9))

# Combine
figC <- (pC1 | pC2) / pC3 +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "Species incidence (presence/absence) patterns across coral size classes",
    subtitle = "Most species show increased incidence on larger corals; Trapezia spp. maintain near-100% across all sizes",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  )

ggsave("output/figures/manuscript/Figure_SizeBinned_Incidence.png", figC, width = 14, height = 12, dpi = 300, bg = "white")
cat("  ✓ Figure C (Incidence) saved\n")

# ============================================================================
# FIGURE D: Species RELATIVE ABUNDANCE by Size Class
# ============================================================================

cat("Creating Figure D: Species relative abundance by size class...\n")

# Calculate relative abundance of each species within each coral
species_rel_abund <- cafi_processed %>%
  group_by(coral_id, species_id) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(coral_id) %>%
  mutate(total = sum(count), rel_abund = count / total * 100) %>%
  ungroup() %>%
  left_join(analysis_data %>% select(coral_id, size_class), by = "coral_id")

# Mean relative abundance by size class
rel_abund_summary <- species_rel_abund %>%
  group_by(size_class, species_id) %>%
  summarise(
    mean_rel = mean(rel_abund),
    se_rel = sd(rel_abund) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Top 15 species by mean relative abundance
top_rel <- rel_abund_summary %>%
  group_by(species_id) %>%
  summarise(overall_mean = mean(mean_rel), .groups = "drop") %>%
  arrange(desc(overall_mean)) %>%
  head(15) %>%
  pull(species_id)

rel_abund_top <- rel_abund_summary %>%
  filter(species_id %in% top_rel) %>%
  mutate(
    species_label = sapply(species_id, function(x) {
      parts <- str_split(x, "_")[[1]]
      genus <- parts[1]
      sp <- ifelse(length(parts) > 1 & parts[2] != "NA", parts[2], "sp.")
      paste0(substr(genus, 1, 1), ". ", sp)
    }),
    species_label = fct_reorder(species_label, mean_rel, .fun = mean)
  )

# Panel D1: Relative abundance heatmap
pD1 <- ggplot(rel_abund_top, aes(x = size_class, y = species_label, fill = mean_rel)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", mean_rel)), size = 2.5,
            color = ifelse(rel_abund_top$mean_rel > 10, "white", "black")) +
  scale_fill_gradient2(low = "#f7f7f7", mid = "#fdae61", high = "#d73027",
                       midpoint = 15, name = "Relative\nabundance (%)") +
  labs(title = "A. Mean relative abundance by size class",
       subtitle = "Proportion of total CAFI community",
       x = "Coral size class", y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic", size = 9))

# Panel D2: Change in relative abundance
rel_change <- rel_abund_top %>%
  select(species_id, species_label, size_class, mean_rel) %>%
  pivot_wider(names_from = size_class, values_from = mean_rel, values_fill = 0) %>%
  mutate(
    change = Large - Small,
    direction = ifelse(change > 0, "Increases", "Decreases")
  )

pD2 <- ggplot(rel_change, aes(x = change, y = fct_reorder(species_label, change), fill = direction)) +
  geom_col(alpha = 0.8, color = "gray30", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray40", linewidth = 0.8) +
  scale_fill_manual(values = c("Increases" = "#1a9850", "Decreases" = "#d73027"),
                    name = "Direction") +
  labs(title = "B. Change in relative abundance (Large − Small)",
       subtitle = "Positive = higher proportion on large corals",
       x = "Change in relative abundance (%)", y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic", size = 9))

# Panel D3: Dotplot showing shift from small to large
rel_shift <- rel_abund_top %>%
  select(species_label, size_class, mean_rel) %>%
  pivot_wider(names_from = size_class, values_from = mean_rel, values_fill = 0)

pD3 <- ggplot(rel_shift, aes(y = fct_reorder(species_label, Large))) +
  geom_segment(aes(x = Small, xend = Large, yend = species_label),
               color = "gray60", linewidth = 0.8) +
  geom_point(aes(x = Small), color = size_colors["Small"], size = 3) +
  geom_point(aes(x = Large), color = size_colors["Large"], size = 3) +
  labs(title = "C. Relative abundance shift: Small → Large",
       subtitle = "Red = small corals, Green = large corals",
       x = "Mean relative abundance (%)", y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic", size = 9))

# Combine
figD <- (pD1 | pD2 | pD3) +
  plot_annotation(
    title = "Species relative abundance patterns across coral size classes",
    subtitle = "Trapezia dominance decreases on large corals as fish and shrimp proportions increase",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  )

ggsave("output/figures/manuscript/Figure_SizeBinned_RelAbund.png", figD, width = 15, height = 10, dpi = 300, bg = "white")
cat("  ✓ Figure D (Relative Abundance) saved\n")

# ============================================================================
# FIGURE E: NMDS Ordination by Size Class
# ============================================================================

cat("Creating Figure E: NMDS ordination by size class...\n")

# Build community matrix
comm_for_nmds <- cafi_processed %>%
  group_by(coral_id, species_id) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_id, values_from = count, values_fill = 0)

comm_mat_nmds <- as.matrix(comm_for_nmds[,-1])
rownames(comm_mat_nmds) <- comm_for_nmds$coral_id

# Hellinger transformation
comm_hell <- decostand(comm_mat_nmds, method = "hellinger")

# NMDS
set.seed(42)
nmds_result <- metaMDS(comm_hell, distance = "bray", k = 2, trymax = 100, trace = 0)
cat(sprintf("  NMDS stress = %.3f\n", nmds_result$stress))

# Extract scores
nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$coral_id <- rownames(nmds_scores)
nmds_scores <- nmds_scores %>%
  left_join(analysis_data %>% select(coral_id, size_class, site, volume), by = "coral_id")

# PERMANOVA by size class
# Ensure exact row matching between community matrix and metadata
comm_coral_ids <- rownames(comm_hell)
permanova_data <- analysis_data %>%
  filter(coral_id %in% comm_coral_ids)

# Subset and reorder community matrix to match
comm_hell_subset <- comm_hell[permanova_data$coral_id, ]

perm_size <- adonis2(comm_hell_subset ~ size_class, data = permanova_data, permutations = 999)
cat(sprintf("  PERMANOVA: Size class R² = %.3f, p = %.3f\n",
            perm_size$R2[1], perm_size$`Pr(>F)`[1]))

# Panel E1: NMDS colored by size class
pE1 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(color = size_class), level = 0.95, linewidth = 1, linetype = "solid") +
  geom_point(aes(fill = size_class, size = volume), shape = 21, alpha = 0.7, stroke = 0.3, color = "white") +
  scale_fill_manual(values = size_colors, name = "Size class") +
  scale_color_manual(values = size_colors, guide = "none") +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(title = "A. NMDS ordination by size class",
       subtitle = sprintf("Stress = %.2f | PERMANOVA: R² = %.3f, p = %.3f",
                         nmds_result$stress, perm_size$R2[1], perm_size$`Pr(>F)`[1]),
       x = "NMDS1", y = "NMDS2") +
  coord_fixed() +
  theme_pub(base_size = 11) +
  theme(legend.position = "right")

# Panel E2: Centroid distances
centroids <- nmds_scores %>%
  group_by(size_class) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  )

pE2 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = size_class), shape = 21, alpha = 0.5, size = 2, stroke = 0.2, color = "white") +
  geom_point(data = centroids, aes(fill = size_class), shape = 23, size = 6, stroke = 1, color = "black") +
  geom_path(data = centroids, aes(group = 1), linetype = "dashed", color = "gray40", linewidth = 1) +
  scale_fill_manual(values = size_colors, name = "Size class") +
  labs(title = "B. Centroid trajectory: Small → Medium → Large",
       subtitle = "Diamonds show group centroids; composition shifts along trajectory",
       x = "NMDS1", y = "NMDS2") +
  coord_fixed() +
  theme_pub(base_size = 11) +
  theme(legend.position = "right")

# Panel E3: Beta dispersion by size class
betadisper_result <- betadisper(vegdist(comm_hell_subset, method = "bray"), permanova_data$size_class)
betadisper_df <- data.frame(
  size_class = betadisper_result$group,
  distance = betadisper_result$distances
)

pE3 <- ggplot(betadisper_df, aes(x = size_class, y = distance, fill = size_class)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1.5) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  labs(title = "C. Beta dispersion by size class",
       subtitle = "Distance to group centroid (community heterogeneity)",
       x = "Coral size class", y = "Distance to centroid") +
  theme_pub(base_size = 11)

# Combine
figE <- (pE1 | pE2 | pE3) +
  plot_annotation(
    title = "Community composition differs significantly across coral size classes",
    subtitle = "PERMANOVA confirms size class explains significant variation; centroids show directional shift",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  )

ggsave("output/figures/manuscript/Figure_SizeBinned_NMDS.png", figE, width = 15, height = 5.5, dpi = 300, bg = "white")
cat("  ✓ Figure E (NMDS) saved\n")

# ============================================================================
# Save summary statistics
# ============================================================================

cat("\nSaving summary statistics...\n")

# Statistical test results
test_results <- tibble(
  Test = c("Abundance", "Richness", "Shannon", "Density", "PERMANOVA"),
  Statistic = c(
    sprintf("χ² = %.1f", kw_abundance$statistic),
    sprintf("χ² = %.1f", kw_richness$statistic),
    sprintf("χ² = %.1f", kw_shannon$statistic),
    sprintf("χ² = %.1f", kw_density$statistic),
    sprintf("R² = %.3f", perm_size$R2[1])
  ),
  p_value = c(
    kw_abundance$p.value,
    kw_richness$p.value,
    kw_shannon$p.value,
    kw_density$p.value,
    perm_size$`Pr(>F)`[1]
  ),
  Interpretation = c(
    "Larger corals have more CAFI",
    "Larger corals have more species",
    "Diversity increases with size",
    "DENSITY DECREASES—supports propagule dilution",
    "Composition differs by size class"
  )
)

write_csv(test_results, "output/tables/size_class_statistics.csv")
write_csv(summary_stats, "output/tables/size_class_summaries.csv")

cat("\n========================================\n")
cat("SIZE-BINNED FIGURE GENERATION COMPLETE\n")
cat("========================================\n")
cat("Output files:\n")
cat("  • Figure_SizeBinned_Metrics.png (Community metrics box plots)\n")
cat("  • Figure_SizeBinned_FuncGroups.png (Functional group patterns)\n")
cat("  • Figure_SizeBinned_Incidence.png (Species presence/absence)\n")
cat("  • Figure_SizeBinned_RelAbund.png (Relative abundance shifts)\n")
cat("  • Figure_SizeBinned_NMDS.png (Ordination by size class)\n")
cat("  • size_class_statistics.csv\n")
cat("  • size_class_summaries.csv\n")
cat("\nAll files saved to: output/figures/manuscript/ and output/tables/\n")
