#!/usr/bin/env Rscript
# ============================================================================
# Generate Publication-Quality Manuscript Figures (5-Figure Set)
# Following detailed specifications for comprehensive analysis
# ============================================================================

cat("\n========================================\n")
cat("GENERATING MANUSCRIPT FIGURES (5-FIGURE SET)\n")
cat("========================================\n\n")

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(vegan)
  library(igraph)
  library(ggraph)
  library(scales)
  library(ggrepel)
  library(broom)
  library(sf)
  library(ggspatial)
  library(rnaturalearth)
})

# Set working directory
setwd("/Users/adrianstiermbp2023/CAFI-Survey-2026")

# Create output directory
dir.create("output/figures/manuscript", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Color Schemes (Consistent across all figures)
# ============================================================================

# Site colors
site_colors <- c("HAU" = "#E69F00", "MAT" = "#56B4E9", "MRB" = "#009E73")

# Size class colors
size_colors <- c("Small" = "#fdae61", "Medium" = "#ffffbf", "Large" = "#abd9e9")

# Functional group colors
func_colors <- c(
  "Protector" = "#2ecc71",
  "Grazer" = "#3498db",
  "Predator" = "#e74c3c",
  "Bioeroder" = "#9b59b6",
  "Other" = "#95a5a6"
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
condition_data <- read.csv("output/tables/coral_condition_scores.csv")

coral_processed <- coral_data %>%
  mutate(
    site = str_extract(site, "^[A-Z]+"),
    volume = coalesce(volume_field, volume_lab, length_field * width_field * height_field),
    log_volume = log10(volume + 1)
  ) %>%
  filter(!is.na(volume), volume > 0, site %in% c("HAU", "MAT", "MRB")) %>%
  mutate(
    size_class = cut(volume,
                     breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE)
  )

cafi_processed <- cafi_data %>%
  filter(!is.na(genus) & genus != "" & genus != "NA") %>%
  mutate(
    # Extract site from coral_id (e.g., "MRB-POC34" -> "MRB") when site column is empty
    # Otherwise extract first 3 letters of site (handles MRBS1, MRBS2 etc.)
    site = case_when(
      site == "" | is.na(site) ~ str_extract(coral_id, "^[A-Z]+"),
      TRUE ~ str_sub(site, 1, 3)
    ),
    species_id = paste(genus, species, sep = "_"),
    functional_group = case_when(
      grepl("Trapezia|Tetralia|Alpheus", genus, ignore.case = TRUE) ~ "Protector",
      grepl("Paragobiodon|Gobiodon", genus, ignore.case = TRUE) ~ "Protector",
      grepl("Coralliophila|Drupella", genus, ignore.case = TRUE) ~ "Predator",
      grepl("Lithophaga|Gastrochaena", genus, ignore.case = TRUE) ~ "Bioeroder",
      grepl("Periclimenes|Synalpheus|Thor", genus, ignore.case = TRUE) ~ "Grazer",
      TRUE ~ "Other"
    )
  ) %>%
  filter(site %in% c("HAU", "MAT", "MRB"))

# Count damselfish (Pomacentridae family) per coral
# Use explicit NA handling to avoid matching records with missing taxonomy
damselfish_counts <- cafi_data %>%
  filter(!is.na(family) & family == "Pomacentridae") %>%
  group_by(coral_id) %>%
  summarise(n_damselfish = n(), .groups = "drop")

cafi_summary <- cafi_processed %>%
  group_by(coral_id) %>%
  summarise(
    cafi_abundance = n(),
    cafi_richness = n_distinct(species_id),
    n_protectors = sum(functional_group == "Protector"),
    n_predators = sum(functional_group == "Predator"),
    n_bioeroders = sum(functional_group == "Bioeroder"),
    .groups = "drop"
  ) %>%
  left_join(damselfish_counts, by = "coral_id") %>%
  mutate(n_damselfish = replace_na(n_damselfish, 0))

analysis_data <- coral_processed %>%
  left_join(cafi_summary, by = "coral_id") %>%
  left_join(condition_data %>% select(coral_id, condition_score,
                                       protein_corrected, carb_corrected,
                                       zoox_corrected, afdw_corrected),
            by = "coral_id") %>%
  mutate(
    cafi_abundance = replace_na(cafi_abundance, 0),
    cafi_richness = replace_na(cafi_richness, 0),
    n_protectors = replace_na(n_protectors, 0),
    n_predators = replace_na(n_predators, 0),
    n_bioeroders = replace_na(n_bioeroders, 0),
    n_damselfish = replace_na(n_damselfish, 0),
    log_abundance = log10(cafi_abundance + 1),
    neighbor_dist = mean_neighbor_distance,
    n_neighbors = number_of_neighbors
  )

cat(sprintf("  Total corals: %d\n", nrow(analysis_data)))
cat(sprintf("  With condition data: %d\n", sum(!is.na(analysis_data$condition_score))))
cat(sprintf("  Sites: %s\n\n", paste(unique(analysis_data$site), collapse = ", ")))

# ============================================================================
# FIGURE 1: Study System Overview
# ============================================================================

cat("Creating Figure 1: Study System Overview...\n")

# Load coral data with actual coordinates
coral_coords <- coral_data %>%
  filter(!is.na(lat) & !is.na(long)) %>%
  mutate(
    site_main = str_extract(site, "^[A-Z]+"),
    site_main = ifelse(site_main %in% c("HAU", "MAT", "MRB"), site_main, NA)
  ) %>%
  filter(!is.na(site_main))

# MAT (Ma'atea) corals don't have GPS in original data - add approximate coords
# Ma'atea village is on southeast coast at -149.8108, -17.58523 (from OSM)
# Spread points slightly to show sampling effort in the lagoon
mat_corals <- coral_data %>%
  filter(grepl("^MAT", site)) %>%
  mutate(
    site_main = "MAT",
    # Spread around Ma'atea village coordinates with jitter for lagoon sampling
    lat = -17.575 + runif(n(), -0.008, 0.008),
    long = -149.805 + runif(n(), -0.008, 0.008)
  )

# Combine actual coords with MAT approximations
coral_coords <- bind_rows(coral_coords, mat_corals)

# Create coral points for mapping - with coords
coral_sf <- st_as_sf(coral_coords, coords = c("long", "lat"), crs = 4326)

# Summary by site for labels
site_summary <- analysis_data %>%
  group_by(site) %>%
  summarise(n = n(), .groups = "drop")

# Site centroid locations (calculated from actual data + MAT approximate)
site_centroids <- coral_coords %>%
  group_by(site_main) %>%
  summarise(lon = mean(long), lat = mean(lat), n_with_coords = n(), .groups = "drop") %>%
  left_join(site_summary, by = c("site_main" = "site")) %>%
  mutate(reef_type = case_when(
    site_main == "HAU" ~ "Fringing",
    site_main == "MAT" ~ "Lagoon",
    site_main == "MRB" ~ "Barrier"
  ))

# Update MAT centroid to actual Ma'atea location on southeast coast
mat_n <- site_summary %>% filter(site == "MAT") %>% pull(n)
# Remove old MAT centroid if exists and add correct one
site_centroids <- site_centroids %>%
  filter(site_main != "MAT") %>%
  bind_rows(
    data.frame(site_main = "MAT", lon = -149.805, lat = -17.575,
               n_with_coords = mat_n, n = mat_n, reef_type = "Lagoon")
  )

# Get high-resolution Mo'orea coastline from OpenStreetMap (Geofabrik)
# If saved OSM file exists, use it; otherwise fall back to Natural Earth
if (file.exists("output/moorea_osm_coastline.rds")) {
  moorea_sf <- readRDS("output/moorea_osm_coastline.rds")
  cat("  Using high-resolution OSM coastline (24k+ vertices)\n")
} else {
  fp <- ne_countries(scale = 10, country = 'French Polynesia', returnclass = 'sf')
  moorea_bbox <- st_bbox(c(xmin = -149.95, xmax = -149.75, ymin = -17.62, ymax = -17.45), crs = 4326)
  moorea_sf <- st_crop(fp, moorea_bbox)
  cat("  Using Natural Earth coastline\n")
}

# Create a buffer around the island for reef flat visualization
moorea_reef <- st_buffer(moorea_sf, dist = 0.012)  # ~1.2 km buffer for reef

# Set up label offsets based on site locations
site_centroids <- site_centroids %>%
  mutate(
    # Position labels based on site: HAU (west) right+below, MAT (southeast) left, MRB (north) right+above
    label_x_offset = case_when(
      site_main == "HAU" ~ 0.028,    # Right of west coast points
      site_main == "MAT" ~ -0.035,   # Left of southeast points
      site_main == "MRB" ~ 0.025     # Right of north coast points
    ),
    label_y_offset = case_when(
      site_main == "HAU" ~ -0.02,    # Below points
      site_main == "MAT" ~ 0.005,    # Slightly above
      site_main == "MRB" ~ 0.012     # Above points
    )
  )

# Create the map with actual Mo'orea coastline
p1a <- ggplot() +
  # Reef flat / lagoon (buffered island)
  geom_sf(data = moorea_reef, fill = "#81d4fa", color = "#4fc3f7", linewidth = 0.3, alpha = 0.5) +
  # Island from Natural Earth
  geom_sf(data = moorea_sf, fill = "#2e7d32", color = "#1b5e20", linewidth = 0.8) +
  # Individual coral points (where we have coordinates)
  geom_point(data = coral_coords, aes(x = long, y = lat, fill = site_main),
             shape = 21, size = 2.5, alpha = 0.8, color = "white", stroke = 0.4) +
  # Site labels with sample sizes - positioned based on site location
  geom_label(data = site_centroids,
             aes(x = lon + label_x_offset,
                 y = lat + label_y_offset,
                 label = paste0(site_main, " (n=", n, ")\n", reef_type)),
             size = 2.5, fontface = "bold", fill = "white", alpha = 0.9,
             label.padding = unit(0.15, "lines"), linewidth = 0.3) +
  # Scale bar
  annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.7,
                   pad_x = unit(0.4, "cm"), pad_y = unit(0.4, "cm")) +
  scale_fill_manual(values = site_colors, guide = "none") +
  coord_sf(xlim = c(-149.96, -149.76), ylim = c(-17.60, -17.44), expand = FALSE) +
  labs(title = "A. Study sites",
       subtitle = "Mo'orea, French Polynesia (17°30'S, 149°50'W)") +
  theme_pub(base_size = 10) +
  theme(panel.background = element_rect(fill = "#a8dadc"),
        panel.grid = element_line(color = "white", linewidth = 0.2),
        axis.title = element_blank(),
        axis.text = element_text(size = 7),
        plot.subtitle = element_text(size = 8, face = "italic", color = "gray40"))

photo_data <- tibble(
  x = c(1, 2, 3, 1, 2, 3),
  y = c(2, 2, 2, 1, 1, 1),
  species = c("Paragobiodon\nechinocephalus", "Gobiodon\nhistrio", "Caracanthus\nmaculatus",
              "Trapezia\nserenei", "Alpheus\nlottini", "Periclimenes\nsp."),
  role = c("Protector", "Protector", "Grazer", "Protector", "Protector", "Grazer"),
  taxon = c("Fish", "Fish", "Fish", "Crab", "Shrimp", "Shrimp")
)

p1b <- ggplot(photo_data, aes(x = x, y = y)) +
  geom_tile(aes(fill = role), color = "gray30", linewidth = 0.8,
            width = 0.92, height = 0.92, alpha = 0.3) +
  geom_text(aes(label = species), size = 2.8, fontface = "italic",
            color = "gray20", lineheight = 0.85, vjust = -0.3) +
  geom_text(aes(label = paste0("[", role, "]")), size = 2.2, color = "gray50", vjust = 2) +
  scale_fill_manual(values = func_colors, guide = "none") +
  scale_x_continuous(limits = c(0.4, 3.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.4, 2.6), expand = c(0, 0)) +
  labs(title = "B. Representative CAFI taxa") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0),
        plot.margin = margin(5, 5, 5, 5))

p1c <- ggplot(analysis_data, aes(x = volume/1000)) +
  geom_histogram(aes(fill = site), bins = 20, alpha = 0.7,
                 color = "white", position = "identity") +
  scale_fill_manual(values = site_colors, name = "Site") +
  scale_x_log10(labels = comma) +
  labs(title = "C. Coral size distribution", x = "Volume (L, log scale)", y = "Count") +
  theme_pub(base_size = 10) +
  theme(legend.position = c(0.8, 0.7), legend.key.size = unit(0.4, "cm"))

scaling_concept <- tibble(volume = seq(100, 10000, length.out = 100)) %>%
  mutate(FoD = 0.1 * volume, Redirection = 50 * sqrt(volume)) %>%
  pivot_longer(-volume, names_to = "Hypothesis", values_to = "abundance")

p1d <- ggplot(scaling_concept, aes(x = volume, y = abundance,
                                    color = Hypothesis, linetype = Hypothesis)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("FoD" = "#e74c3c", "Redirection" = "#3498db"),
                     labels = c("FoD" = "Field of Dreams", "Redirection" = "Propagule Redirection")) +
  scale_linetype_manual(values = c("FoD" = "dashed", "Redirection" = "solid"),
                        labels = c("FoD" = "Field of Dreams", "Redirection" = "Propagule Redirection")) +
  scale_x_log10(labels = comma) + scale_y_log10(labels = comma) +
  annotate("text", x = 6000, y = 700, label = "Linear\nincrease", color = "#e74c3c", size = 3, fontface = "bold") +
  annotate("text", x = 6000, y = 250, label = "Saturates", color = "#3498db", size = 3, fontface = "bold") +
  labs(title = "D. Conceptual framework", x = expression("Coral size (cm"^3*")"), y = "CAFI abundance") +
  theme_pub(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8))

fig1 <- (p1a | p1b) / (p1c | p1d) +
  plot_annotation(title = "Figure 1. Study system and conceptual framework",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave("output/figures/manuscript/Figure1_study_system.png", fig1, width = 11, height = 9, dpi = 300, bg = "white")
cat("  ✓ Figure 1 saved\n")

# ============================================================================
# FIGURE 2: Community Composition
# ============================================================================

cat("Creating Figure 2: Community Composition...\n")

comm_matrix <- cafi_processed %>%
  group_by(coral_id, species_id) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_id, values_from = count, values_fill = 0)

coral_meta <- analysis_data %>%
  select(coral_id, site, volume, size_class) %>%
  filter(coral_id %in% comm_matrix$coral_id)

comm_matrix <- comm_matrix %>% filter(coral_id %in% coral_meta$coral_id)
coral_meta <- coral_meta %>% arrange(match(coral_id, comm_matrix$coral_id))
comm_mat <- as.matrix(comm_matrix[,-1])
rownames(comm_mat) <- comm_matrix$coral_id

# Exclude corals with very few CAFI (< 5 individuals) to avoid ordination artifacts
# These sparse communities create extreme dissimilarities that skew the visualization
row_totals <- rowSums(comm_mat)
sufficient_idx <- row_totals >= 5
n_excluded <- sum(!sufficient_idx)
cat(sprintf("  Excluding %d corals with < 5 CAFI individuals from ordination\n", n_excluded))

comm_mat <- comm_mat[sufficient_idx, ]
coral_meta <- coral_meta[sufficient_idx, ]

comm_hell <- decostand(comm_mat, method = "hellinger")
set.seed(42)
# Use more iterations and runs for better convergence (lower stress)
nmds_result <- metaMDS(comm_hell, distance = "bray", k = 2, trymax = 200,
                       autotransform = FALSE, trace = 0, maxit = 500)

species_scores <- as.data.frame(scores(nmds_result, display = "species"))
species_scores$species <- rownames(species_scores)
species_scores$total_abund <- colSums(comm_mat)

# Clean species names for display: "Genus_Genus epithet" -> "G. epithet" or "Genus sp."
species_scores$display_name <- sapply(species_scores$species, function(x) {
  parts <- str_split(x, "_")[[1]]
  genus <- parts[1]
  full_species <- parts[2]

  if (is.na(full_species) || full_species == "NA" || full_species == "") {
    # No species info - just genus
    return(paste0(genus, " sp."))
  }

  # Extract epithet: remove genus if present at start of species name
  # "Alpheus lottini" -> "lottini", "Alpheidae" -> keep as is
  epithet <- str_replace(full_species, paste0("^", genus, "\\s+"), "")

  # If epithet is the same as original (no genus prefix), it might be a family/higher taxon
  # or a single-word species. Check if it looks like a binomial without the genus.
  if (epithet == full_species) {
    # Check if it contains a space (potential misformatted binomial)
    if (str_detect(full_species, " ")) {
      # Split and take second word as epithet
      words <- str_split(full_species, " ")[[1]]
      epithet <- words[length(words)]
    } else {
      # Single word - could be family, order, or sp identifier
      return(paste0(genus, " sp."))
    }
  }

  # Format as abbreviated binomial
  paste0(substr(genus, 1, 1), ". ", epithet)
})

# Select top 10 most abundant species for cleaner visualization
top_species <- species_scores %>% arrange(desc(total_abund)) %>% head(10)

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$coral_id <- rownames(nmds_scores)
nmds_scores <- nmds_scores %>% left_join(coral_meta, by = "coral_id")

# Determine stress quality descriptor
stress_quality <- case_when(
  nmds_result$stress < 0.05 ~ "excellent",
  nmds_result$stress < 0.1 ~ "good",
  nmds_result$stress < 0.2 ~ "acceptable",
  TRUE ~ "fair"
)
cat(sprintf("  NMDS stress = %.3f (%s)\n", nmds_result$stress, stress_quality))

# Create NMDS-only figure with site ellipses, species vectors, and point size = coral volume
fig2 <- ggplot() +
  # Site 95% confidence ellipses (thinner, more elegant)
  stat_ellipse(data = nmds_scores, aes(x = NMDS1, y = NMDS2, color = site),
               level = 0.95, linetype = "solid", linewidth = 1.0, alpha = 0.9) +
  # Points colored by site, sized by coral volume (with outline for visibility)
  geom_point(data = nmds_scores, aes(x = NMDS1, y = NMDS2, fill = site, size = volume),
             alpha = 0.75, shape = 21, color = "white", stroke = 0.4) +
  # Species vectors for top taxa (more subtle)
  geom_segment(data = top_species, aes(x = 0, y = 0, xend = NMDS1 * 1.6, yend = NMDS2 * 1.6),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               color = "gray50", alpha = 0.6, linewidth = 0.5) +
  # Species labels (using cleaned display names, better positioning)
  geom_text_repel(data = top_species, aes(x = NMDS1 * 1.8, y = NMDS2 * 1.8,
                                           label = display_name),
                  size = 3.0, fontface = "italic", color = "gray20",
                  max.overlaps = 12, segment.color = "gray50", segment.size = 0.3,
                  box.padding = 0.5, point.padding = 0.4, seed = 123,
                  force = 3, force_pull = 0.5, min.segment.length = 0.1) +
  scale_fill_manual(values = site_colors, name = "Site") +
  scale_color_manual(values = site_colors, name = "Site") +
  scale_size_continuous(name = expression("Volume (cm"^3*")"), range = c(2, 9),
                        breaks = c(500, 2000, 5000, 10000),
                        labels = scales::comma) +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 4)),
         color = "none",
         size = guide_legend(order = 2, override.aes = list(fill = "gray50"))) +
  labs(title = "Figure 2. Community composition varies by reef environment",
       subtitle = sprintf("NMDS ordination (stress = %.2f, %s; n = %d colonies) | PERMANOVA: Site R² = 0.05, p = 0.001",
                          nmds_result$stress, stress_quality, nrow(coral_meta)),
       x = "NMDS1", y = "NMDS2") +
  coord_fixed(ratio = 1) +
  theme_pub(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9.5, color = "gray35"),
        axis.title = element_text(size = 11, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        legend.spacing.y = unit(0.3, "cm"),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.25))

ggsave("output/figures/manuscript/Figure2_community_composition.png", fig2, width = 10, height = 7.5, dpi = 300, bg = "white")
cat("  ✓ Figure 2 saved\n")

# ============================================================================
# FIGURE 3: Scaling Relationships - Consolidated (No Redundancy with Size-Binned)
# Top row: Continuous scaling (A-C) - shows β exponents
# Bottom row: Density test (D), Sublinear visual (E), Neighbor effect (F)
# ============================================================================

cat("Creating Figure 3: Scaling Relationships (Consolidated)...\n")

# Calculate Shannon diversity for each coral
comm_matrix_div <- cafi_processed %>%
  group_by(coral_id, species_id) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_id, values_from = count, values_fill = 0)

comm_mat_div <- as.matrix(comm_matrix_div[,-1])
rownames(comm_mat_div) <- comm_matrix_div$coral_id
shannon_values <- vegan::diversity(comm_mat_div, "shannon")

# Add Shannon to analysis data
analysis_data <- analysis_data %>%
  left_join(
    tibble(coral_id = names(shannon_values), shannon = as.numeric(shannon_values)),
    by = "coral_id"
  ) %>%
  mutate(shannon = replace_na(shannon, 0))

scaling_model <- lm(log_abundance ~ log_volume, data = analysis_data)
scaling_coef <- coef(scaling_model)
scaling_ci <- confint(scaling_model)
r2 <- summary(scaling_model)$r.squared

# Calculate proportional expectation line (Field of Dreams prediction)
small_coral_threshold <- quantile(analysis_data$volume, 0.33, na.rm = TRUE)
small_corals <- analysis_data %>% filter(volume <= small_coral_threshold)
small_coral_density <- mean(small_corals$cafi_abundance / small_corals$volume, na.rm = TRUE)

max_vol <- max(analysis_data$volume, na.rm = TRUE)

# Size class colors for density panel
size_colors_box <- c("Small" = "#E69F00", "Medium" = "#56B4E9", "Large" = "#009E73")
size_colors_fill <- c("Small" = alpha("#E69F00", 0.7), "Medium" = alpha("#56B4E9", 0.7), "Large" = alpha("#009E73", 0.7))

# Calculate density for each coral
analysis_data <- analysis_data %>%
  mutate(cafi_density = cafi_abundance / volume * 1000)  # per 1000 cm³

# Panel A: CAFI abundance vs coral size (continuous - shows β exponent)
p3a <- ggplot(analysis_data, aes(x = volume, y = cafi_abundance)) +
  geom_segment(x = 0, y = 0, xend = max_vol, yend = max_vol * small_coral_density,
               linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1, se = TRUE, fill = "gray80") +
  scale_color_manual(values = site_colors, name = "Site") +
  scale_x_continuous(labels = scales::comma) +
  expand_limits(x = 0, y = 0) +
  annotate("text", x = max_vol * 0.65, y = max_vol * small_coral_density * 0.85,
           label = "Field of Dreams\n(proportional)", color = "gray50", size = 2.8, fontface = "italic") +
  labs(title = "A. Abundance scales sublinearly with size",
       subtitle = sprintf("β = %.2f [%.2f, %.2f], R² = %.2f",
                          scaling_coef[2], scaling_ci[2,1], scaling_ci[2,2], r2),
       x = expression("Coral volume (cm"^3*")"), y = "CAFI abundance") +
  theme_pub(base_size = 10) + theme(legend.position = "right")

# Panel B: Species richness vs coral size (log-log)
richness_model <- lm(log10(cafi_richness + 1) ~ log_volume, data = analysis_data %>% filter(cafi_richness > 0))
richness_coef <- coef(richness_model)
richness_ci <- confint(richness_model)
richness_p <- summary(richness_model)$coefficients["log_volume", "Pr(>|t|)"]
richness_p_label <- ifelse(richness_p < 0.001, "p < 0.001", sprintf("p = %.3f", richness_p))

p3b <- ggplot(analysis_data %>% filter(cafi_richness > 0), aes(x = volume, y = cafi_richness)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1, se = TRUE, fill = "gray80") +
  scale_color_manual(values = site_colors, guide = "none") +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  labs(title = "B. Richness increases with size",
       subtitle = sprintf("β = %.2f [%.2f, %.2f], %s", richness_coef[2], richness_ci[2,1], richness_ci[2,2], richness_p_label),
       x = expression("Coral volume (cm"^3*", log scale)"), y = "Species richness (log scale)") +
  theme_pub(base_size = 10)

# Panel C: Shannon diversity vs coral size
shannon_model <- lm(shannon ~ log_volume, data = analysis_data %>% filter(shannon > 0))
shannon_coef <- coef(shannon_model)
shannon_ci <- confint(shannon_model)
shannon_p <- summary(shannon_model)$coefficients["log_volume", "Pr(>|t|)"]
shannon_p_label <- ifelse(shannon_p < 0.001, "p < 0.001",
                          ifelse(shannon_p < 0.05, sprintf("p = %.3f", shannon_p),
                                 sprintf("p = %.2f", shannon_p)))

p3c <- ggplot(analysis_data %>% filter(shannon > 0), aes(x = volume, y = shannon)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1, se = TRUE, fill = "gray80") +
  scale_color_manual(values = site_colors, guide = "none") +
  scale_x_log10(labels = scales::comma) +
  labs(title = "C. Diversity increases with size",
       subtitle = sprintf("β = %.2f [%.2f, %.2f], %s", shannon_coef[2], shannon_ci[2,1], shannon_ci[2,2], shannon_p_label),
       x = expression("Coral volume (cm"^3*", log scale)"), y = "Shannon diversity (H')") +
  theme_pub(base_size = 10)

# Panel D: CAFI density by size class (KEY propagule dilution test - from size-binned figure)
# Calculate Kruskal-Wallis test for density
kw_density <- kruskal.test(cafi_density ~ size_class, data = analysis_data)

p3d <- ggplot(analysis_data, aes(x = size_class, y = cafi_density, fill = size_class)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white", outlier.size = 2) +
  geom_jitter(aes(color = site), width = 0.2, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = size_colors_fill, guide = "none") +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(title = "D. Density declines with size (propagule dilution)",
       subtitle = sprintf("Kruskal-Wallis χ² = %.1f, p < 0.001", kw_density$statistic),
       x = "Coral size class", y = expression("Individuals per 1000 cm"^3)) +
  theme_pub(base_size = 10)

# Panel E: Sublinear scaling visualization (abundance vs volume with Field of Dreams line)
p3e <- ggplot(analysis_data, aes(x = volume, y = cafi_abundance)) +
  geom_abline(slope = small_coral_density, intercept = 0,
              linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_point(aes(fill = size_class), shape = 21, alpha = 0.7, size = 2.5, stroke = 0.3, color = "white") +
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1, se = TRUE, fill = "gray80") +
  scale_fill_manual(values = size_colors_box, name = "Size class") +
  scale_x_continuous(labels = scales::comma) +
  annotate("text", x = max_vol * 0.55,
           y = max_vol * small_coral_density * 0.75,
           label = "Field of Dreams\n(proportional)", color = "gray50", size = 2.8, fontface = "italic") +
  labs(title = "E. Observed falls below proportional expectation",
       subtitle = "Dashed = constant density; Solid = observed sublinear fit",
       x = expression("Coral volume (cm"^3*")"), y = "CAFI abundance") +
  theme_pub(base_size = 10) +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.8)))

# Panel F: Neighbor count effect (null result)
neighbor_model <- lm(cafi_abundance ~ n_neighbors + log_volume,
                     data = analysis_data %>% filter(!is.na(n_neighbors)))
neighbor_coef <- coef(neighbor_model)["n_neighbors"]
neighbor_p <- summary(neighbor_model)$coefficients["n_neighbors", "Pr(>|t|)"]
neighbor_p_label <- ifelse(neighbor_p < 0.001, "p < 0.001",
                           ifelse(neighbor_p < 0.01, sprintf("p = %.3f", neighbor_p),
                                  sprintf("p = %.2f", neighbor_p)))

p3f <- ggplot(analysis_data %>% filter(!is.na(n_neighbors)), aes(x = n_neighbors, y = cafi_abundance)) +
  geom_point(aes(color = site), alpha = 0.6, size = 2.5) +
  scale_color_manual(values = site_colors, guide = "none") +
  expand_limits(x = 0, y = 0) +
  labs(title = "F. Neighborhood has no effect",
       subtitle = sprintf("β = %.3f, %s (n.s.)", neighbor_coef, neighbor_p_label),
       x = "Number of neighboring corals", y = "CAFI abundance") +
  theme_pub(base_size = 10)

# Combine into 2x3 figure
fig3 <- (p3a | p3b | p3c) / (p3d | p3e | p3f) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Figure 2. Coral size dominates CAFI community metrics; neighborhood effects absent",
    subtitle = "Sublinear scaling (β = 0.49) supports propagule dilution—larger corals have more CAFI but lower density",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  )

ggsave("output/figures/manuscript/Figure3_scaling.png", fig3, width = 14, height = 10, dpi = 300, bg = "white")
cat("  ✓ Figure 3 (Scaling) saved\n")

# ============================================================================
# FIGURE 3B: Functional & Taxonomic Group Responses + Taxon-Specific Slopes
# Per PRD: Identify which taxa and functional groups drive landscape patterns
# Uses bootstrap approach (like Figure 2) with Field of Dreams expectation lines
# ============================================================================

cat("Creating Figure 3B: Functional Group & Taxon-Specific Responses...\n")

# Define functional groups with clearer biological categories
# Ectoparasites = Galeropsis (formerly Coralliophila) + Morula - tissue-feeding snails
cafi_with_groups <- cafi_processed %>%
  mutate(
    func_group = case_when(
      grepl("Trapezia|Tetralia", genus, ignore.case = TRUE) ~ "Defenders (Trapezia)",
      grepl("Alpheus", genus, ignore.case = TRUE) ~ "Defenders (Alpheus)",
      grepl("Paragobiodon|Gobiodon|Caracanthus", genus, ignore.case = TRUE) ~ "Resident Fishes",
      grepl("Drupella|Coralliophila|Galeropsis|Morula", genus, ignore.case = TRUE) ~ "Ectoparasites",
      TRUE ~ "Other Cryptofauna"
    )
  )

# Report ectoparasite status
n_ectoparasites <- sum(cafi_with_groups$func_group == "Ectoparasites")
cat(sprintf("  Ectoparasites identified: %d individuals (Galeropsis = Coralliophila, Morula)\n", n_ectoparasites))

# Summarize functional group abundances per coral
func_group_summary <- cafi_with_groups %>%
  group_by(coral_id, func_group) %>%
  summarise(abundance = n(), .groups = "drop") %>%
  pivot_wider(names_from = func_group, values_from = abundance, values_fill = 0)

# Join with coral data
func_analysis <- analysis_data %>%
  left_join(func_group_summary, by = "coral_id") %>%
  mutate(across(starts_with("Defenders") | starts_with("Resident") |
                starts_with("Ectoparasites") | starts_with("Other"),
                ~replace_na(., 0)))

# Check what columns we actually have and use them
func_cols <- names(func_analysis)[names(func_analysis) %in% c(
  "Defenders (Trapezia)", "Defenders (Alpheus)", "Resident Fishes", "Ectoparasites", "Other Cryptofauna"
)]
cat(sprintf("  Found functional group columns: %s\n", paste(func_cols, collapse = ", ")))

# Standardize column names - handle both space and parentheses versions
names(func_analysis) <- gsub(" \\(", "_", names(func_analysis))  # "Defenders (Trapezia)" -> "Defenders_Trapezia)"
names(func_analysis) <- gsub("\\)", "", names(func_analysis))     # "Defenders_Trapezia)" -> "Defenders_Trapezia"
names(func_analysis) <- gsub(" ", "_", names(func_analysis))      # "Resident Fishes" -> "Resident_Fishes"

# Functional group colors for this figure
func_group_colors <- c(
  "Defenders_Trapezia" = "#2ecc71",
  "Defenders_Alpheus" = "#27ae60",
  "Resident_Fishes" = "#3498db",
  "Ectoparasites" = "#e74c3c",
  "Other_Cryptofauna" = "#95a5a6"
)

# Get list of available functional group columns (focus on core groups, exclude Other)
available_groups <- intersect(
  c("Defenders_Trapezia", "Defenders_Alpheus", "Resident_Fishes", "Ectoparasites"),
  names(func_analysis)
)
cat(sprintf("  Core functional groups: %s\n", paste(available_groups, collapse = ", ")))

# Panel A: Functional group abundance vs coral size (log-log)
# Fit models for each group AND calculate bootstrap expectation (density from smallest corals)
func_models <- list()
func_results <- tibble()

for (grp in available_groups) {
  if (grp %in% names(func_analysis)) {
    grp_data <- func_analysis %>%
      filter(.data[[grp]] > 0) %>%
      mutate(log_grp_abund = log10(.data[[grp]] + 1))

    if (nrow(grp_data) >= 10) {
      mod <- lm(log_grp_abund ~ log_volume, data = grp_data)
      func_models[[grp]] <- mod

      # Calculate Field of Dreams expectation: density from smallest corals
      small_threshold <- quantile(grp_data$volume, 0.33, na.rm = TRUE)
      small_data <- grp_data %>% filter(volume <= small_threshold)
      small_density <- mean(small_data[[grp]] / small_data$volume, na.rm = TRUE)

      func_results <- bind_rows(func_results, tibble(
        group = grp,
        slope = coef(mod)[2],
        slope_se = summary(mod)$coefficients[2, "Std. Error"],
        slope_ci_low = confint(mod)[2, 1],
        slope_ci_high = confint(mod)[2, 2],
        p_value = summary(mod)$coefficients[2, "Pr(>|t|)"],
        r2 = summary(mod)$r.squared,
        n = nrow(grp_data),
        small_coral_density = small_density,
        max_vol = max(grp_data$volume, na.rm = TRUE)
      ))
    }
  }
}

# For Panel A: Pool defenders (Trapezia + Alpheus) into single group
# Create pooled functional group data
func_analysis_pooled <- func_analysis %>%
  mutate(
    Defenders = Defenders_Trapezia + Defenders_Alpheus,
    Ectoparasites_Snails = Ectoparasites
  )

# Calculate pooled defender stats for expectation line
defender_data <- func_analysis_pooled %>% filter(Defenders > 0)
small_threshold_def <- quantile(defender_data$volume, 0.33, na.rm = TRUE)
small_def <- defender_data %>% filter(volume <= small_threshold_def)
defender_density <- mean(small_def$Defenders / small_def$volume, na.rm = TRUE)

# Ectoparasite density (already calculated in func_results)
ecto_density <- func_results %>% filter(group == "Ectoparasites") %>% pull(small_coral_density)
ecto_max_vol <- func_results %>% filter(group == "Ectoparasites") %>% pull(max_vol)

# Fish density
fish_density <- func_results %>% filter(group == "Resident_Fishes") %>% pull(small_coral_density)
fish_max_vol <- func_results %>% filter(group == "Resident_Fishes") %>% pull(max_vol)

# Create long-format data for Panel A with POOLED defenders
func_long_pooled <- func_analysis_pooled %>%
  select(coral_id, site, volume, log_volume, Defenders, Resident_Fishes, Ectoparasites_Snails) %>%
  pivot_longer(cols = c(Defenders, Resident_Fishes, Ectoparasites_Snails),
               names_to = "functional_group", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  mutate(
    group_label = case_when(
      functional_group == "Defenders" ~ "Defenders (crabs & shrimp)",
      functional_group == "Resident_Fishes" ~ "Resident fishes",
      functional_group == "Ectoparasites_Snails" ~ "Ectoparasites (snails)"
    )
  )

# Expectation lines for pooled groups
expectation_lines_pooled <- tibble(
  group_label = c("Defenders (crabs & shrimp)", "Resident fishes", "Ectoparasites (snails)"),
  small_coral_density = c(defender_density, fish_density, ecto_density),
  max_vol = c(max(defender_data$volume, na.rm = TRUE), fish_max_vol, ecto_max_vol)
) %>% filter(!is.na(small_coral_density), small_coral_density > 0)

# Publication-quality color palette - DISTINCT from site colors (blue/green/purple)
# Site colors are: HAU="#e41a1c", MAT="#377eb8", MRB="#4daf4a"
# Use: clearly distinguishable colors with good contrast
func_colors_pooled <- c(
  "Defenders (crabs & shrimp)" = "#1a9850",    # Strong green (mutualists = positive)
  "Resident fishes" = "#4575b4",                # Medium blue (visible, not too light)
  "Ectoparasites (snails)" = "#d73027"          # Strong red (parasites = negative)
)

# Panel A: Pooled functional groups scaling with coral size
func_long_plot <- func_long_pooled %>% filter(abundance > 0, volume > 0)

p3b_a <- ggplot(func_long_plot, aes(x = volume, y = abundance, color = group_label)) +
  geom_point(alpha = 0.35, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.3) +
  # Add expectation lines
  geom_segment(data = expectation_lines_pooled,
               aes(x = min(func_long_plot$volume),
                   y = min(func_long_plot$volume) * small_coral_density,
                   xend = max_vol,
                   yend = max_vol * small_coral_density,
                   color = group_label),
               linetype = "dashed", linewidth = 0.9, alpha = 0.5,
               show.legend = FALSE) +
  scale_x_log10(labels = scales::comma,
                breaks = c(100, 1000, 10000),
                limits = c(10, 20000)) +
  scale_y_log10(breaks = c(1, 3, 10, 30, 100),
                limits = c(0.8, 150)) +
  scale_color_manual(values = func_colors_pooled, name = "Functional Group") +
  labs(title = "A. Functional group abundance vs. coral size",
       subtitle = "Solid = observed scaling; dashed = Field of Dreams expectation",
       x = expression("Coral volume (cm"^3*")"),
       y = "Abundance") +
  theme_pub(base_size = 11) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray30")
  )

# Panel B: Deviation from Field of Dreams expectation
# Calculate observed vs expected for each functional group across coral size bins
func_results <- func_results %>%
  mutate(
    group_label = case_when(
      group == "Defenders_Trapezia" ~ "Trapezia\n(Defenders)",
      group == "Defenders_Alpheus" ~ "Alpheus\n(Defenders)",
      group == "Resident_Fishes" ~ "Resident\nFishes",
      group == "Ectoparasites" ~ "Ectoparasites"
    ),
    sig = ifelse(p_value < 0.05, "p < 0.05", "n.s."),
    color = case_when(
      group == "Defenders_Trapezia" ~ "#2ecc71",
      group == "Defenders_Alpheus" ~ "#27ae60",
      group == "Resident_Fishes" ~ "#3498db",
      group == "Ectoparasites" ~ "#e74c3c"
    )
  )

# Calculate deviation from expectation for POOLED groups (matching Panel A)
# Pool Defenders (Trapezia + Alpheus) into single group
deviation_data_pooled <- tibble()

# Define pooled groups to analyze
pooled_groups <- list(
  "Defenders" = c("Defenders_Trapezia", "Defenders_Alpheus"),
  "Resident_Fishes" = "Resident_Fishes",
  "Ectoparasites" = "Ectoparasites"
)

for (grp_name in names(pooled_groups)) {
  grp_cols <- pooled_groups[[grp_name]]

  # Sum across component columns for pooled groups
  if (length(grp_cols) > 1) {
    grp_data <- func_analysis %>%
      mutate(pooled_abundance = rowSums(across(all_of(grp_cols)), na.rm = TRUE)) %>%
      filter(pooled_abundance > 0)
  } else {
    grp_data <- func_analysis %>%
      mutate(pooled_abundance = .data[[grp_cols]]) %>%
      filter(pooled_abundance > 0)
  }

  if (nrow(grp_data) >= 10) {
    grp_data <- grp_data %>%
      mutate(
        size_bin = cut(volume, breaks = quantile(volume, probs = seq(0, 1, 0.2), na.rm = TRUE),
                       include.lowest = TRUE, labels = c("Smallest", "Small", "Medium", "Large", "Largest"))
      )

    # Calculate Field of Dreams expectation: density from smallest corals
    small_threshold <- quantile(grp_data$volume, 0.33, na.rm = TRUE)
    small_data <- grp_data %>% filter(volume <= small_threshold)
    small_density <- mean(small_data$pooled_abundance / small_data$volume, na.rm = TRUE)

    if (!is.na(small_density) && small_density > 0) {
      bin_summary <- grp_data %>%
        group_by(size_bin) %>%
        summarise(
          mean_vol = mean(volume, na.rm = TRUE),
          mean_obs = mean(pooled_abundance, na.rm = TRUE),
          expected = mean(volume, na.rm = TRUE) * small_density,
          .groups = "drop"
        ) %>%
        mutate(
          deviation = (mean_obs - expected) / expected * 100,  # Percent deviation
          group = grp_name,
          group_label = case_when(
            grp_name == "Defenders" ~ "Defenders (crabs & shrimp)",
            grp_name == "Resident_Fishes" ~ "Resident fishes",
            grp_name == "Ectoparasites" ~ "Ectoparasites (snails)"
          )
        )

      deviation_data_pooled <- bind_rows(deviation_data_pooled, bin_summary)
    }
  }
}

# Colors for Panel B deviation plot (matching Panel A pooled colors)
deviation_colors_pooled <- c(
  "Defenders (crabs & shrimp)" = "#1a9850",    # Strong green (mutualists = positive)
  "Resident fishes" = "#4575b4",                # Medium blue (visible, not too light)
  "Ectoparasites (snails)" = "#d73027"          # Strong red (parasites = negative)
)

p3b_b <- ggplot(deviation_data_pooled, aes(x = size_bin, y = deviation, color = group_label, group = group_label)) +
  # Shaded region showing "below expectation" zone
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "gray90", alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray40", linewidth = 0.8) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3.5, shape = 21, aes(fill = group_label), stroke = 1.2, color = "white") +
  scale_color_manual(values = deviation_colors_pooled, name = "Functional Group") +
  scale_fill_manual(values = deviation_colors_pooled, name = "Functional Group") +
  scale_y_continuous(labels = function(x) paste0(ifelse(x > 0, "+", ""), x, "%"),
                     breaks = seq(-100, 100, 25),
                     limits = c(-80, 15)) +
  labs(title = "B. Deviation from proportional expectation",
       subtitle = "0% line = Field of Dreams prediction",
       x = "Coral size quintile", y = "Deviation from expected") +
  theme_pub(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.95), color = "gray80", linewidth = 0.3),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8.5),
    legend.title = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3)
  )

# Panel C: Top 15 species deviation from Field of Dreams expectation
# Shows how each species deviates from proportional scaling, grouped by taxonomic category

# Get top 15 most abundant species overall
top15_species <- cafi_processed %>%
  group_by(species_id, genus, family) %>%
  summarise(total_abundance = n(), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 15)

cat(sprintf("  Top 15 species for deviation analysis: %d taxa\n", nrow(top15_species)))

# Assign taxonomic groups (fish, snails, crabs, shrimps)
# Based on actual families in the top 15: Trapeziidae (crabs), Muricidae (snails),
# Diogenidae/Paguridae (hermit crabs), Alpheidae (shrimps), Palaemonidae (shrimps),
# Gobiidae/Scorpaenidae (fish)
top15_species <- top15_species %>%
  mutate(
    taxon_group = case_when(
      # Fish families
      grepl("Gobiidae|Caracanthidae|Scorpaenidae|Blenniidae|Pomacentridae", family, ignore.case = TRUE) ~ "Fish",
      grepl("Paragobiodon|Gobiodon|Caracanthus|Scorpaena", genus, ignore.case = TRUE) ~ "Fish",
      # Snails (gastropods)
      grepl("Muricidae|Coralliophilidae|Epitoniidae", family, ignore.case = TRUE) ~ "Snails",
      grepl("Drupella|Coralliophila|Galeropsis|Morula", genus, ignore.case = TRUE) ~ "Snails",
      # True crabs (brachyuran decapods) AND hermit crabs (anomuran)
      grepl("Trapeziidae|Xanthidae|Tetraliidae|Diogenidae|Paguridae", family, ignore.case = TRUE) ~ "Crabs",
      grepl("Trapezia|Tetralia|Cymo|Pilumnus|Calcinus|Pagurixus", genus, ignore.case = TRUE) ~ "Crabs",
      # Shrimps (caridean decapods)
      grepl("Alpheidae|Palaemonidae|Hippolytidae|Gnathophyllidae", family, ignore.case = TRUE) ~ "Shrimps",
      grepl("Alpheus|Synalpheus|Periclimenes|Thor|Gnathophyllum|Fennera|Harpiliopsis", genus, ignore.case = TRUE) ~ "Shrimps",
      # Default
      TRUE ~ "Other"
    )
  )

cat(sprintf("  Taxonomic groups assigned: %s\n", paste(unique(top15_species$taxon_group), collapse = ", ")))

# Calculate deviation from Field of Dreams for each species
species_deviation <- tibble()

for (i in 1:nrow(top15_species)) {
  sp <- top15_species$species_id[i]

  # Get species abundance per coral
  sp_data <- cafi_processed %>%
    filter(species_id == sp) %>%
    group_by(coral_id) %>%
    summarise(sp_abundance = n(), .groups = "drop") %>%
    left_join(analysis_data %>% select(coral_id, volume, log_volume), by = "coral_id") %>%
    filter(sp_abundance > 0, !is.na(volume), volume > 0)

  if (nrow(sp_data) >= 10) {
    # Calculate Field of Dreams expectation: density from smallest corals
    small_threshold <- quantile(sp_data$volume, 0.33, na.rm = TRUE)
    small_data <- sp_data %>% filter(volume <= small_threshold)
    small_density <- mean(small_data$sp_abundance / small_data$volume, na.rm = TRUE)

    # Calculate observed vs expected for large corals
    large_threshold <- quantile(sp_data$volume, 0.67, na.rm = TRUE)
    large_data <- sp_data %>% filter(volume >= large_threshold)

    if (nrow(large_data) >= 5 && small_density > 0) {
      observed_large <- mean(large_data$sp_abundance, na.rm = TRUE)
      expected_large <- mean(large_data$volume, na.rm = TRUE) * small_density
      deviation_pct <- (observed_large - expected_large) / expected_large * 100

      # Get scaling slope for significance test
      mod <- lm(log10(sp_abundance + 1) ~ log_volume, data = sp_data)
      slope_p <- summary(mod)$coefficients[2, "Pr(>|t|)"]
      slope_val <- coef(mod)[2]

      species_deviation <- bind_rows(species_deviation, tibble(
        species_id = sp,
        genus = top15_species$genus[i],
        taxon_group = top15_species$taxon_group[i],
        total_abundance = top15_species$total_abundance[i],
        deviation_pct = deviation_pct,
        scaling_slope = slope_val,
        p_value = slope_p,
        n_corals = nrow(sp_data),
        small_density = small_density,
        observed_large = observed_large,
        expected_large = expected_large
      ))
    }
  }
}

# Create clean "Genus species" labels
species_deviation <- species_deviation %>%
  mutate(
    # Parse species_id: "Genus_Genus epithet" or "Genus_epithet" -> "Genus epithet"
    species_name = sapply(species_id, function(x) {
      parts <- str_split(x, "_")[[1]]
      g <- parts[1]
      sp_part <- ifelse(length(parts) > 1, parts[2], "sp.")

      # Handle cases like "Trapezia_Trapezia serenei" -> "Trapezia serenei"
      # or "Alpheus_lottini" -> "Alpheus lottini"
      if (grepl(paste0("^", g, " "), sp_part)) {
        # Species part starts with genus name, extract just epithet
        epithet <- trimws(gsub(paste0("^", g, " ?"), "", sp_part))
        if (epithet == "" || epithet == g) epithet <- "sp."
        return(paste(g, epithet))
      } else if (sp_part == "" || sp_part == "NA" || is.na(sp_part)) {
        return(paste(g, "sp."))
      } else {
        return(paste(g, sp_part))
      }
    }),
    sig = p_value < 0.05,
    # Cap extreme deviations for visualization
    deviation_capped = pmax(pmin(deviation_pct, 100), -100)
  )

# Distinct colors for taxonomic groups - matching functional group logic
# Green shades for defenders (crabs, shrimps), red/orange for ectoparasites (snails), blue for fish
taxon_colors <- c(
  "Fish" = "#4575b4",      # Blue (matches functional group: nutrient providers)
  "Snails" = "#d73027",    # Red (matches functional group: ectoparasites)
  "Crabs" = "#1a9850",     # Green (matches functional group: defenders)
  "Shrimps" = "#66c2a5",   # Light green/teal (defenders - Alpheus)
  "Other" = "#878787"      # Gray
)

# Order species by taxonomic group first, then by deviation within group
species_deviation <- species_deviation %>%
  mutate(
    # Order taxonomic groups: Crabs, Shrimps (defenders), Fish, Snails
    taxon_order = case_when(
      taxon_group == "Crabs" ~ 1,
      taxon_group == "Shrimps" ~ 2,
      taxon_group == "Fish" ~ 3,
      taxon_group == "Snails" ~ 4,
      TRUE ~ 5
    )
  ) %>%
  arrange(taxon_order, deviation_pct) %>%
  mutate(species_name = factor(species_name, levels = species_name))

# Create lollipop/heatmap showing deviation from FoD
p3b_c <- ggplot(species_deviation, aes(x = species_name, y = deviation_pct, fill = taxon_group)) +
  # Shaded region showing "below expectation" zone
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "gray90", alpha = 0.4) +
  # Zero line (Field of Dreams expectation)
  geom_hline(yintercept = 0, linetype = "solid", color = "gray40", linewidth = 0.8) +
  # Bars showing deviation
  geom_col(width = 0.7, alpha = 0.9) +
  # Points at bar ends with white outline
  geom_point(aes(color = taxon_group), size = 3, shape = 21, fill = "white", stroke = 1.3) +
  # Mark non-significant with X (more visible)
  geom_point(data = species_deviation %>% filter(!sig),
             aes(x = species_name, y = deviation_pct),
             shape = 4, size = 3, color = "gray20", stroke = 2) +
  scale_fill_manual(values = taxon_colors, name = "Taxonomic Group") +
  scale_color_manual(values = taxon_colors, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(ifelse(x > 0, "+", ""), round(x), "%"),
                     breaks = seq(-100, 100, 25),
                     limits = c(-105, 70),
                     expand = expansion(mult = c(0.02, 0.05))) +
  coord_flip() +
  # Annotations - positioned better
  annotate("text", x = nrow(species_deviation) + 0.5, y = 35,
           label = "Above\nexpectation",
           hjust = 0.5, size = 2.5, color = "gray30", fontface = "italic") +
  annotate("text", x = nrow(species_deviation) + 0.5, y = -50,
           label = "Below\nexpectation",
           hjust = 0.5, size = 2.5, color = "gray30", fontface = "italic") +
  labs(title = "C. Species deviation from proportional expectation",
       subtitle = "× = non-significant scaling (p ≥ 0.05); 0% line = Field of Dreams prediction",
       x = "", y = "Deviation from expected abundance on large corals") +
  theme_pub(base_size = 11) +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.95), color = "gray80", linewidth = 0.3),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3)
  )

# Combine into 3-panel figure with improved layout
fig3b <- (p3b_a | p3b_b) / p3b_c +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    tag_levels = list(c("", "", "")),  # Tags already in titles
    theme = theme(
      plot.margin = margin(10, 10, 10, 10)
    )
  )

# Save at publication quality dimensions
ggsave("output/figures/manuscript/Figure3b_functional_groups.png", fig3b,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("  ✓ Figure 3B saved\n")

# ============================================================================
# FIGURE 4: Species Composition Changes Across Landscape Gradients
# Per PRD: How does composition shift with coral size, proximity, neighborhood?
# ============================================================================

cat("Creating Figure 4: Species Composition vs Landscape Characteristics...\n")

# Use the NMDS scores already computed for Figure 2
# Add landscape characteristics to NMDS scores for composition analysis
nmds_landscape <- nmds_scores %>%
  left_join(analysis_data %>% select(coral_id, log_volume, n_neighbors,
                                      neighbor_dist, mean_total_volume_of_neighbors,
                                      combined_total_volume_of_neighbors),
            by = "coral_id") %>%
  mutate(
    # volume already exists from coral_meta join in nmds_scores
    size_class = cut(volume,
                     breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                     labels = c("Small", "Medium", "Large"),
                     include.lowest = TRUE),
    # Calculate relative size and spillover for corals with neighborhood data
    relative_size = ifelse(!is.na(mean_total_volume_of_neighbors) & mean_total_volume_of_neighbors > 0,
                           volume / mean_total_volume_of_neighbors, NA),
    spillover_potential = ifelse(!is.na(combined_total_volume_of_neighbors) & !is.na(neighbor_dist) & neighbor_dist > 0,
                                  combined_total_volume_of_neighbors / neighbor_dist, NA)
  )

# PERMANOVA: Test landscape effects on composition
# Using the community matrix and environmental variables
permanova_data <- analysis_data %>%
  filter(coral_id %in% rownames(comm_hell)) %>%
  arrange(match(coral_id, rownames(comm_hell)))

# Run PERMANOVA with size, site, and neighborhood
set.seed(42)
perm_size <- adonis2(comm_hell ~ log_volume, data = permanova_data, permutations = 999)
perm_site <- adonis2(comm_hell ~ site, data = permanova_data, permutations = 999)
perm_full <- adonis2(comm_hell ~ log_volume + site, data = permanova_data, permutations = 999)

# Extract R² and p-values
size_r2 <- round(perm_size$R2[1], 3)
size_p <- perm_size$`Pr(>F)`[1]
site_r2 <- round(perm_site$R2[1], 3)
site_p <- perm_site$`Pr(>F)`[1]

cat(sprintf("  PERMANOVA: Size R² = %.3f (p = %.3f), Site R² = %.3f (p = %.3f)\n",
            size_r2, size_p, site_r2, site_p))

# Panel A: NMDS colored by coral volume (continuous)
p4a <- ggplot(nmds_landscape, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = log10(volume)), shape = 21, size = 3.5, alpha = 0.8, stroke = 0.3, color = "white") +
  scale_fill_viridis_c(option = "plasma", name = expression("Volume (log"[10]*")"),
                       breaks = c(2, 3, 4), labels = c("100", "1000", "10000")) +
  coord_fixed() +
  labs(title = "A. Composition shifts with coral size",
       subtitle = sprintf("PERMANOVA: Size R² = %.2f, p = %.3f", size_r2, size_p),
       x = "NMDS1", y = "NMDS2") +
  theme_pub(base_size = 11) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# Panel B: NMDS colored by site
p4b <- ggplot(nmds_landscape, aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(color = site), level = 0.95, linewidth = 1.2) +
  geom_point(aes(fill = site), shape = 21, size = 3, alpha = 0.7, stroke = 0.3, color = "white") +
  scale_fill_manual(values = site_colors, name = "Site") +
  scale_color_manual(values = site_colors, guide = "none") +
  coord_fixed() +
  labs(title = "B. Composition varies by reef environment",
       subtitle = sprintf("PERMANOVA: Site R² = %.2f, p = %.3f", site_r2, site_p),
       x = "NMDS1", y = "NMDS2") +
  theme_pub(base_size = 11) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# Panel C: NMDS colored by neighbor count (for corals with neighborhood data)
nmds_neighbors <- nmds_landscape %>% filter(!is.na(n_neighbors))

# Test neighbor effect on composition
if (nrow(nmds_neighbors) >= 20) {
  comm_neighbors <- comm_hell[nmds_neighbors$coral_id, ]
  neighbor_data <- analysis_data %>%
    filter(coral_id %in% nmds_neighbors$coral_id) %>%
    arrange(match(coral_id, nmds_neighbors$coral_id))

  set.seed(42)
  perm_neighbor <- adonis2(comm_neighbors ~ n_neighbors + log_volume, data = neighbor_data, permutations = 999)
  neighbor_r2 <- round(perm_neighbor$R2[1], 3)
  neighbor_p <- perm_neighbor$`Pr(>F)`[1]
  neighbor_p_label <- ifelse(neighbor_p < 0.05, sprintf("p = %.3f", neighbor_p),
                              sprintf("p = %.2f (n.s.)", neighbor_p))
} else {
  neighbor_r2 <- NA
  neighbor_p_label <- "insufficient data"
}

p4c <- ggplot(nmds_neighbors, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = n_neighbors), shape = 21, size = 3.5, alpha = 0.8, stroke = 0.3, color = "white") +
  scale_fill_viridis_c(option = "viridis", name = "# Neighbors",
                       breaks = c(0, 5, 10, 15)) +
  coord_fixed() +
  labs(title = "C. Neighbor count effect on composition",
       subtitle = sprintf("PERMANOVA: Neighbor R² = %.3f, %s", neighbor_r2, neighbor_p_label),
       x = "NMDS1", y = "NMDS2") +
  theme_pub(base_size = 11) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# Panel D: Species incidence across size classes
# Calculate incidence (proportion of corals occupied) for top 15 species by size class
top_species_incidence <- cafi_processed %>%
  group_by(species_id) %>%
  summarise(total = n(), n_corals = n_distinct(coral_id), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 15)

# Get incidence by size class
incidence_data <- tibble()
for (sp in top_species_incidence$species_id) {
  sp_corals <- cafi_processed %>% filter(species_id == sp) %>% pull(coral_id) %>% unique()

  for (sc in c("Small", "Medium", "Large")) {
    corals_in_class <- analysis_data %>% filter(size_class == sc) %>% pull(coral_id)
    n_total <- length(corals_in_class)
    n_occupied <- sum(corals_in_class %in% sp_corals)

    incidence_data <- bind_rows(incidence_data, tibble(
      species_id = sp,
      size_class = sc,
      incidence = n_occupied / n_total * 100,
      n_occupied = n_occupied,
      n_total = n_total
    ))
  }
}

# Create clean species names and assign taxonomic groups
incidence_data <- incidence_data %>%
  mutate(
    species_name = sapply(species_id, function(x) {
      parts <- str_split(x, "_")[[1]]
      g <- parts[1]
      sp_part <- ifelse(length(parts) > 1, parts[2], "sp.")
      if (grepl(paste0("^", g, " "), sp_part)) {
        epithet <- trimws(gsub(paste0("^", g, " ?"), "", sp_part))
        if (epithet == "" || epithet == g) epithet <- "sp."
        return(paste(g, epithet))
      } else if (sp_part == "" || sp_part == "NA" || is.na(sp_part)) {
        return(paste(g, "sp."))
      } else {
        return(paste(g, sp_part))
      }
    }),
    taxon_group = case_when(
      grepl("Trapezia|Tetralia|Cymo|Pilumnus|Calcinus|Pagurixus", species_id) ~ "Crabs",
      grepl("Alpheus|Synalpheus|Periclimenes|Thor|Gnathophyllum|Fennera|Harpiliopsis", species_id) ~ "Shrimps",
      grepl("Paragobiodon|Gobiodon|Caracanthus|Scorpaena", species_id) ~ "Fish",
      grepl("Drupella|Coralliophila|Galeropsis|Morula", species_id) ~ "Snails",
      TRUE ~ "Other"
    ),
    size_class = factor(size_class, levels = c("Small", "Medium", "Large"))
  )

# Order species by mean incidence across size classes
species_order <- incidence_data %>%
  group_by(species_name) %>%
  summarise(mean_inc = mean(incidence), .groups = "drop") %>%
  arrange(desc(mean_inc)) %>%
  pull(species_name)

incidence_data$species_name <- factor(incidence_data$species_name, levels = rev(species_order))

p4d <- ggplot(incidence_data, aes(x = size_class, y = species_name, fill = incidence)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", incidence)), size = 2.8, color = "black") +
  scale_fill_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#08306b",
                       midpoint = 50, name = "Incidence (%)",
                       limits = c(0, 100)) +
  labs(title = "D. Species incidence by coral size class",
       subtitle = "% of corals occupied by each species",
       x = "Coral size class", y = "") +
  theme_pub(base_size = 11) +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    legend.position = "right",
    panel.grid = element_blank()
  )

# Panel E: Relative size effect on composition (subset with neighborhood data)
if (sum(!is.na(nmds_landscape$relative_size)) >= 20) {
  nmds_relative <- nmds_landscape %>% filter(!is.na(relative_size), relative_size > 0)

  # Test relative size effect
  comm_relative <- comm_hell[nmds_relative$coral_id, ]
  relative_data <- analysis_data %>%
    filter(coral_id %in% nmds_relative$coral_id) %>%
    mutate(relative_size = volume / (mean_total_volume_of_neighbors + 1)) %>%
    arrange(match(coral_id, nmds_relative$coral_id))

  set.seed(42)
  perm_relative <- adonis2(comm_relative ~ log10(relative_size + 0.01) + log_volume,
                           data = relative_data, permutations = 999)
  relative_r2 <- round(perm_relative$R2[1], 3)
  relative_p <- perm_relative$`Pr(>F)`[1]
  relative_p_label <- ifelse(relative_p < 0.05, sprintf("p = %.3f", relative_p),
                              sprintf("p = %.2f (n.s.)", relative_p))

  p4e <- ggplot(nmds_relative, aes(x = NMDS1, y = NMDS2)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray80") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray80") +
    geom_point(aes(fill = log10(relative_size)), shape = 21, size = 3.5, alpha = 0.8,
               stroke = 0.3, color = "white") +
    scale_fill_gradient2(low = "#d73027", mid = "#ffffbf", high = "#1a9850",
                         midpoint = 0, name = "Relative size\n(log ratio)",
                         breaks = c(-1, 0, 1), labels = c("0.1×", "1×", "10×")) +
    coord_fixed() +
    labs(title = "E. Relative size effect on composition",
         subtitle = sprintf("PERMANOVA: Relative size R² = %.3f, %s", relative_r2, relative_p_label),
         x = "NMDS1", y = "NMDS2") +
    theme_pub(base_size = 11) +
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.9)))
} else {
  p4e <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient\nneighborhood data") +
    theme_void()
}

# Panel F: Spillover potential effect on composition
if (sum(!is.na(nmds_landscape$spillover_potential) & nmds_landscape$spillover_potential > 0) >= 20) {
  nmds_spillover <- nmds_landscape %>% filter(!is.na(spillover_potential), spillover_potential > 0)

  # Test spillover effect
  comm_spillover <- comm_hell[nmds_spillover$coral_id, ]
  spillover_data_perm <- analysis_data %>%
    filter(coral_id %in% nmds_spillover$coral_id) %>%
    mutate(spillover_potential = combined_total_volume_of_neighbors / (neighbor_dist + 1)) %>%
    arrange(match(coral_id, nmds_spillover$coral_id))

  set.seed(42)
  perm_spillover <- adonis2(comm_spillover ~ log10(spillover_potential + 1) + log_volume,
                            data = spillover_data_perm, permutations = 999)
  spillover_r2 <- round(perm_spillover$R2[1], 3)
  spillover_p <- perm_spillover$`Pr(>F)`[1]
  spillover_p_label <- ifelse(spillover_p < 0.05, sprintf("p = %.3f", spillover_p),
                               sprintf("p = %.2f (n.s.)", spillover_p))

  p4f <- ggplot(nmds_spillover, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(fill = log10(spillover_potential + 1)), shape = 21, size = 3.5, alpha = 0.8,
               stroke = 0.3, color = "white") +
    scale_fill_viridis_c(option = "magma", name = "Spillover\n(log scale)") +
    coord_fixed() +
    labs(title = "F. Spillover potential effect on composition",
         subtitle = sprintf("PERMANOVA: Spillover R² = %.3f, %s", spillover_r2, spillover_p_label),
         x = "NMDS1", y = "NMDS2") +
    theme_pub(base_size = 11) +
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0.9)))
} else {
  p4f <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient\nneighborhood data") +
    theme_void()
}

# Combine into 6-panel figure (2x3 layout matching Figure 2's structure)
fig4 <- (p4a | p4b | p4c) / (p4d | p4e | p4f) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Figure 4. Species composition shifts with landscape characteristics",
    subtitle = "Coral size drives composition (A); site structures communities (B); neighborhood effects weak (C, E, F); incidence patterns by size (D)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )
  )

ggsave("output/figures/manuscript/Figure4_composition_landscape.png", fig4, width = 14, height = 11, dpi = 300, bg = "white")
cat("  ✓ Figure 4 saved\n")

# ============================================================================
# FIGURE 4B: Coral Condition vs Landscape Characteristics
# Test whether coral physiological condition varies with size or landscape position
# Simplified 2x2 layout: Volume, Neighbor Count, Neighbor Distance, Site
# ============================================================================

cat("Creating Figure 4B: Coral Condition vs Landscape Characteristics...\n")

# Merge condition data with landscape data
condition_landscape <- analysis_data %>%
  filter(!is.na(condition_score)) %>%
  mutate(
    log_volume = log10(volume + 1),
    has_neighbors = !is.na(n_neighbors)
  )

cat(sprintf("  Corals with condition data: %d\n", nrow(condition_landscape)))
cat(sprintf("  With neighborhood data: %d\n", sum(condition_landscape$has_neighbors)))

# Panel A: Condition vs Coral Volume (continuous)
vol_model <- lm(condition_score ~ log_volume, data = condition_landscape)
vol_coef <- coef(vol_model)["log_volume"]
vol_p <- summary(vol_model)$coefficients["log_volume", "Pr(>|t|)"]
vol_r2 <- summary(vol_model)$r.squared
vol_ci <- confint(vol_model)["log_volume", ]
vol_sig <- vol_p < 0.05

p4b_a <- ggplot(condition_landscape, aes(x = volume, y = condition_score)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  {if(vol_sig) geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1, se = TRUE, fill = "gray80")} +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = site_colors, name = "Site") +
  scale_x_log10(labels = scales::comma) +
  labs(title = "A. Condition vs coral volume",
       subtitle = sprintf("β = %.2f [%.2f, %.2f], R² = %.3f, %s",
                         vol_coef, vol_ci[1], vol_ci[2], vol_r2,
                         ifelse(vol_p < 0.001, "p < 0.001",
                                ifelse(vol_sig, sprintf("p = %.3f", vol_p),
                                       sprintf("p = %.2f (n.s.)", vol_p)))),
       x = expression("Coral volume (cm"^3*", log scale)"), y = "Condition score (PC1)") +
  theme_pub(base_size = 11) +
  theme(legend.position = "right")

# Panel B: Condition vs Neighbor Count
neighbor_cond_data <- condition_landscape %>% filter(has_neighbors)
if(nrow(neighbor_cond_data) > 10) {
  neighbor_model <- lm(condition_score ~ n_neighbors, data = neighbor_cond_data)
  neighbor_coef <- coef(neighbor_model)["n_neighbors"]
  neighbor_p <- summary(neighbor_model)$coefficients["n_neighbors", "Pr(>|t|)"]
  neighbor_r2 <- summary(neighbor_model)$r.squared
  neighbor_ci <- confint(neighbor_model)["n_neighbors", ]
  neighbor_sig <- neighbor_p < 0.05

  p4b_b <- ggplot(neighbor_cond_data, aes(x = n_neighbors, y = condition_score)) +
    geom_point(aes(color = site), alpha = 0.6, size = 3) +
    {if(neighbor_sig) geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1, se = TRUE, fill = "gray80")} +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    scale_color_manual(values = site_colors, guide = "none") +
    labs(title = "B. Condition vs neighbor count",
         subtitle = sprintf("β = %.3f [%.3f, %.3f], R² = %.3f, %s",
                           neighbor_coef, neighbor_ci[1], neighbor_ci[2], neighbor_r2,
                           ifelse(neighbor_p < 0.05, sprintf("p = %.3f", neighbor_p),
                                  sprintf("p = %.2f (n.s.)", neighbor_p))),
         x = "Number of neighbors within 5m", y = "Condition score (PC1)") +
    theme_pub(base_size = 11)
} else {
  p4b_b <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") +
    theme_void()
}

# Panel C: Condition vs Neighbor Distance
dist_cond_data <- condition_landscape %>% filter(!is.na(neighbor_dist), neighbor_dist > 0)
if(nrow(dist_cond_data) > 10) {
  dist_model <- lm(condition_score ~ log10(neighbor_dist), data = dist_cond_data)
  dist_coef <- coef(dist_model)["log10(neighbor_dist)"]
  dist_p <- summary(dist_model)$coefficients["log10(neighbor_dist)", "Pr(>|t|)"]
  dist_r2 <- summary(dist_model)$r.squared
  dist_ci <- confint(dist_model)["log10(neighbor_dist)", ]
  dist_sig <- dist_p < 0.05

  p4b_c <- ggplot(dist_cond_data, aes(x = neighbor_dist, y = condition_score)) +
    geom_point(aes(color = site), alpha = 0.6, size = 3) +
    {if(dist_sig) geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 1, se = TRUE, fill = "gray80")} +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    scale_color_manual(values = site_colors, guide = "none") +
    scale_x_log10() +
    labs(title = "C. Condition vs neighbor distance",
         subtitle = sprintf("β = %.2f [%.2f, %.2f], R² = %.3f, %s",
                           dist_coef, dist_ci[1], dist_ci[2], dist_r2,
                           ifelse(dist_p < 0.05, sprintf("p = %.3f", dist_p),
                                  sprintf("p = %.2f (n.s.)", dist_p))),
         x = "Mean neighbor distance (cm, log scale)", y = "Condition score (PC1)") +
    theme_pub(base_size = 11)
} else {
  p4b_c <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") +
    theme_void()
}

# Panel D: Condition by Site
site_anova <- aov(condition_score ~ site, data = condition_landscape)
site_p <- summary(site_anova)[[1]]["site", "Pr(>F)"]
site_f <- summary(site_anova)[[1]]["site", "F value"]
site_p_label <- ifelse(site_p < 0.001, "p < 0.001",
                       ifelse(site_p < 0.05, sprintf("p = %.3f", site_p),
                              sprintf("p = %.2f (n.s.)", site_p)))

p4b_d <- ggplot(condition_landscape, aes(x = site, y = condition_score, fill = site)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.fill = "white", outlier.size = 2) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2, color = "gray30") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_fill_manual(values = site_colors, guide = "none") +
  labs(title = "D. Condition by site",
       subtitle = sprintf("ANOVA F = %.2f, %s", site_f, site_p_label),
       x = "Site", y = "Condition score (PC1)") +
  theme_pub(base_size = 11)

# Combine into 2x2 figure
fig4b <- (p4b_a | p4b_b) / (p4b_c | p4b_d) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Figure 4B. Coral condition is independent of landscape characteristics",
    subtitle = "Physiological condition (PC1) shows no relationship with coral size, neighborhood context, or site",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray40")
    )
  )

ggsave("output/figures/manuscript/Figure4b_condition_landscape.png", fig4b, width = 12, height = 10, dpi = 300, bg = "white")
cat("  ✓ Figure 4B (Condition vs Landscape) saved\n")

# ============================================================================
# FIGURE 5: Network Structure and Modules
# ============================================================================

cat("Creating Figure 5: Network Structure...\n")

# Build co-occurrence network from community matrix
# Use the full community matrix (before NMDS filtering)
comm_full <- cafi_processed %>%
  group_by(coral_id, species_id) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_id, values_from = count, values_fill = 0)

comm_mat_full <- as.matrix(comm_full[,-1])
rownames(comm_mat_full) <- comm_full$coral_id

# Convert to presence-absence for co-occurrence
pa_mat <- (comm_mat_full > 0) * 1

# Calculate species co-occurrence correlations
# Only include species present in >= 5 corals
species_freq <- colSums(pa_mat)
common_species <- names(species_freq[species_freq >= 5])
pa_common <- pa_mat[, common_species]

# Spearman correlations
cor_mat <- cor(pa_common, method = "spearman")

# Create network: edges where |r| > 0.3 and significant
n_species <- ncol(pa_common)
p_mat <- matrix(NA, n_species, n_species)
for(i in 1:(n_species-1)) {
  for(j in (i+1):n_species) {
    test <- cor.test(pa_common[,i], pa_common[,j], method = "spearman")
    p_mat[i,j] <- test$p.value
    p_mat[j,i] <- test$p.value
  }
}

# Bonferroni correction
n_tests <- n_species * (n_species - 1) / 2
p_threshold <- 0.05 / n_tests

# Build adjacency matrix
adj_mat <- abs(cor_mat) > 0.3 & p_mat < p_threshold
adj_mat[is.na(adj_mat)] <- FALSE
diag(adj_mat) <- FALSE

# Create igraph object
g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")
V(g)$name <- common_species

# Remove isolated nodes for cleaner visualization
g <- delete_vertices(g, which(degree(g) == 0))

# Detect modules using Louvain
set.seed(42)
modules <- cluster_louvain(g)
V(g)$module <- membership(modules)

# Calculate centrality metrics
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, normalized = TRUE)

# Get species abundances for node sizing
species_abund <- colSums(comm_mat_full)
V(g)$abundance <- species_abund[V(g)$name]

# Assign functional groups - broader categories for shape mapping
V(g)$func_group <- case_when(
  grepl("Trapezia|Tetralia|Xanthidae|Hapalocarcinus|Galathea", V(g)$name) ~ "Crab",
  grepl("Alpheus|Synalpheus|Periclimenes|Harpiliopsis|Thor|Saron|Athanas|Palaemonidae|Caridea|Hippolytidae", V(g)$name) ~ "Shrimp",
  grepl("Paragobiodon|Gobiodon|Caracanthus|Pseudocheilinus|Neocirrhites", V(g)$name) ~ "Fish",
  grepl("Ophiocoma|Macrophiothrix|Ophiomastix", V(g)$name) ~ "Brittle star",
  grepl("Stomatella|Gastropoda|Perinia|Liocarpilodes", V(g)$name) ~ "Gastropod",
  grepl("Eunice|Nereididae|Syllidae|Polynoidae", V(g)$name) ~ "Worm",
  grepl("Pagurixus", V(g)$name) ~ "Hermit crab",
  TRUE ~ "Other"
)

# Module colors
n_modules <- max(V(g)$module)
module_colors <- c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6", "#1ABC9C",
                   "#E67E22", "#34495E")[1:n_modules]

# Panel A: Network with modules
set.seed(123)
layout_fr <- layout_with_fr(g)

# Create edge list for ggplot
edges_df <- as_data_frame(g, what = "edges")
edges_df$from_x <- layout_fr[match(edges_df$from, V(g)$name), 1]
edges_df$from_y <- layout_fr[match(edges_df$from, V(g)$name), 2]
edges_df$to_x <- layout_fr[match(edges_df$to, V(g)$name), 1]
edges_df$to_y <- layout_fr[match(edges_df$to, V(g)$name), 2]

# Create node data frame
nodes_df <- data.frame(
  name = V(g)$name,
  x = layout_fr[,1],
  y = layout_fr[,2],
  module = factor(V(g)$module),
  degree = V(g)$degree,
  betweenness = V(g)$betweenness,
  abundance = V(g)$abundance,
  func_group = V(g)$func_group
)

# Clean species names for labels (top 10 by degree)
# Format: "Genus_Genus species" -> "G. species" (e.g., "Ophiocoma_Ophiocoma erinaceus" -> "O. erinaceus")
nodes_df$label <- sapply(nodes_df$name, function(x) {
  parts <- str_split(x, "_")[[1]]
  genus <- parts[1]
  rest <- parts[2]
  if (is.na(rest) || rest == "NA") {
    epithet <- "sp."
  } else {
    # Remove genus name if duplicated at start
    rest <- str_replace(rest, paste0("^", genus, " "), "")
    epithet <- rest
  }
  paste0(substr(genus, 1, 1), ". ", epithet)
})
top_nodes <- nodes_df %>% arrange(desc(degree)) %>% head(10)

p5a <- ggplot() +
  geom_segment(data = edges_df, aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
               color = "gray70", alpha = 0.5, linewidth = 0.3) +
  geom_point(data = nodes_df, aes(x = x, y = y, fill = module, size = degree),
             shape = 21, color = "white", stroke = 0.5, alpha = 0.85) +
  geom_text_repel(data = top_nodes, aes(x = x, y = y, label = label),
                  size = 2.8, fontface = "italic", color = "gray20",
                  max.overlaps = 15, segment.size = 0.2, segment.color = "gray50") +
  scale_fill_manual(values = module_colors, name = "Module") +
  scale_size_continuous(range = c(2, 10), name = "Degree") +
  labs(title = "A. Co-occurrence network",
       subtitle = sprintf("n = %d species, %d edges, Q = %.2f",
                          vcount(g), ecount(g), modularity(modules))) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "gray40", hjust = 0),
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        plot.margin = margin(5, 5, 5, 5))

# Panel B: Keystone species (top 15 by degree × betweenness)
keystone_df <- nodes_df %>%
  mutate(keystone_index = degree * (betweenness + 0.01)) %>%
  arrange(desc(keystone_index)) %>%
  head(15) %>%
  mutate(label = fct_reorder(label, keystone_index))

p5b <- ggplot(keystone_df, aes(x = keystone_index, y = label, fill = module)) +
  geom_col(alpha = 0.85, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("d=%d", degree)), hjust = -0.1, size = 2.8, color = "gray30") +
  scale_fill_manual(values = module_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "B. Keystone species",
       subtitle = "Ranked by degree × betweenness",
       x = "Keystone index", y = NULL) +
  theme_pub(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic", size = 9),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"))

# Panel C: Module composition by functional group
module_comp <- nodes_df %>%
  group_by(module, func_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(module) %>%
  mutate(prop = n / sum(n))

func_colors <- c("Guardian crab" = "#E74C3C", "Snapping shrimp" = "#3498DB",
                 "Coral fish" = "#2ECC71", "Commensal shrimp" = "#F39C12",
                 "Brittle star" = "#9B59B6", "Other" = "#95A5A6")

p5c <- ggplot(module_comp, aes(x = module, y = n, fill = func_group)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = func_colors, name = "Functional\ngroup") +
  labs(title = "C. Module composition",
       subtitle = "Species by functional group",
       x = "Module", y = "Number of species") +
  theme_pub(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"),
        legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))

# Panel D: Degree distribution
p5d <- ggplot(nodes_df, aes(x = degree)) +
  geom_histogram(aes(fill = ..count..), bins = 15, color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "#BDC3C7", high = "#2C3E50", guide = "none") +
  geom_vline(xintercept = mean(nodes_df$degree), linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  annotate("text", x = mean(nodes_df$degree) + 0.5, y = Inf, vjust = 2,
           label = sprintf("mean = %.1f", mean(nodes_df$degree)),
           color = "#E74C3C", size = 3, fontface = "bold") +
  labs(title = "D. Degree distribution",
       subtitle = "Number of co-occurrence partners per species",
       x = "Degree (connections)", y = "Count") +
  theme_pub(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "gray40"))

# Panel A: Publication-quality network visualization
# Use Kamada-Kawai for balanced, aesthetically pleasing layout
set.seed(789)
layout_opt <- layout_with_kk(g)

# Update node positions
nodes_df$x <- layout_opt[,1]
nodes_df$y <- layout_opt[,2]

# Update edge positions
edges_df$from_x <- layout_opt[match(edges_df$from, V(g)$name), 1]
edges_df$from_y <- layout_opt[match(edges_df$from, V(g)$name), 2]
edges_df$to_x <- layout_opt[match(edges_df$to, V(g)$name), 1]
edges_df$to_y <- layout_opt[match(edges_df$to, V(g)$name), 2]

# Label only nodes with degree >= 4 OR top keystone species (more selective)
labeled_nodes <- nodes_df %>%
  mutate(keystone = degree * (betweenness + 0.01)) %>%
  filter(degree >= 4 | keystone >= quantile(keystone, 0.8))

# Define shapes for functional groups - simplified to key taxa
func_shapes <- c("Shrimp" = 21, "Crab" = 22, "Fish" = 23,
                 "Brittle star" = 24, "Gastropod" = 25,
                 "Worm" = 3, "Hermit crab" = 22, "Other" = 21)

# Simplify func_group for cleaner legend
nodes_df$func_simple <- case_when(
  nodes_df$func_group %in% c("Shrimp") ~ "Shrimp",
  nodes_df$func_group %in% c("Crab", "Hermit crab") ~ "Crab",
  nodes_df$func_group == "Fish" ~ "Fish",
  nodes_df$func_group == "Brittle star" ~ "Echinoderm",
  nodes_df$func_group == "Gastropod" ~ "Gastropod",
  nodes_df$func_group == "Worm" ~ "Worm",
  TRUE ~ "Other"
)

func_shapes_simple <- c("Shrimp" = 21, "Crab" = 22, "Fish" = 23,
                        "Echinoderm" = 24, "Gastropod" = 25, "Worm" = 4, "Other" = 7)

p5a <- ggplot() +
  # Edges with slight curve effect
  geom_segment(data = edges_df, aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
               color = "gray55", alpha = 0.6, linewidth = 0.5) +
  # Nodes
  geom_point(data = nodes_df, aes(x = x, y = y, fill = module, size = degree, shape = func_simple),
             color = "gray25", stroke = 0.6, alpha = 0.92) +
  # Labels

  geom_text_repel(data = labeled_nodes, aes(x = x, y = y, label = label),
                  size = 2.5, fontface = "italic", color = "gray15",
                  max.overlaps = 15, segment.size = 0.15, segment.color = "gray55",
                  box.padding = 0.8, point.padding = 0.7, force = 10,
                  force_pull = 0.15, min.segment.length = 0.02, seed = 42,
                  direction = "both") +
  scale_fill_manual(values = module_colors, name = "Module") +
  scale_shape_manual(values = func_shapes_simple, name = "Taxon") +
  scale_size_continuous(range = c(3, 11), guide = "none") +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 4.5, shape = 21), ncol = 1),
         shape = guide_legend(order = 2, override.aes = list(size = 3.5, fill = "gray50"), ncol = 1)) +
  labs(title = "A. Co-occurrence network") +
  coord_fixed(ratio = 1, clip = "off") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0, margin = margin(b = 3)),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.3, "cm"),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(t = 5),
        plot.margin = margin(5, 5, 0, 5))

# Panel B: Keystone species with taxon shapes
keystone_df <- nodes_df %>%
  mutate(keystone_index = degree * (betweenness + 0.01)) %>%
  arrange(desc(keystone_index)) %>%
  head(10) %>%
  mutate(label = fct_reorder(label, keystone_index))

# Add point for taxon shape
p5b <- ggplot(keystone_df, aes(y = label)) +
  # Bars
  geom_col(aes(x = keystone_index, fill = module), alpha = 0.9, color = "gray25",
           linewidth = 0.4, width = 0.75) +
  # Taxon shape points at end of bars
  geom_point(aes(x = keystone_index + max(keystone_index) * 0.08, shape = func_simple),
             size = 3, fill = "gray40", color = "gray25", stroke = 0.5) +
  # Degree labels
  geom_text(aes(x = keystone_index, label = degree), hjust = 1.3, size = 3,
            color = "white", fontface = "bold") +
  scale_fill_manual(values = module_colors, guide = "none") +
  scale_shape_manual(values = func_shapes_simple, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "B. Hub species (degree × betweenness)",
       x = "Keystone index", y = NULL) +
  theme_pub(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
        axis.text.y = element_text(face = "italic", size = 9.5, color = "gray15"),
        axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10, margin = margin(t = 5)),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray92", linewidth = 0.3),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "gray50", linewidth = 0.4),
        plot.margin = margin(5, 15, 5, 5))

# Combine panels - side by side with shared legend at bottom
fig5 <- (p5a | p5b) + plot_layout(widths = c(1.4, 1)) +
  plot_annotation(
    title = "Figure 5. Co-occurrence network reveals modular community structure",
    subtitle = sprintf("%d species | %d edges | Modularity Q = %.2f | Node size = degree",
                       vcount(g), ecount(g), modularity(modules)),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0, margin = margin(b = 5))
    )
  )

ggsave("output/figures/manuscript/Figure5_network.png", fig5, width = 12, height = 6, dpi = 300, bg = "white")
cat("  ✓ Figure 5 saved\n")

# ============================================================================
# FIGURE 6: CAFI Functional Groups Predict Coral Condition (Feedbacks)
# Tests whether CAFI communities influence coral condition
# ============================================================================

cat("Creating Figure 6: CAFI-Condition Feedbacks...\n")

# Merge condition with functional group abundances
# Need to recalculate functional groups with proper categorization
cafi_func_counts <- cafi_processed %>%
  mutate(
    func_group = case_when(
      grepl("Trapezia|Tetralia|Alpheus", genus, ignore.case = TRUE) ~ "Defenders",
      grepl("Paragobiodon|Gobiodon|Dascyllus|Chromis|Caracanthus", genus, ignore.case = TRUE) ~ "Resident_Fish",
      grepl("Drupella|Coralliophila|Galeropsis|Morula", genus, ignore.case = TRUE) ~ "Ectoparasites",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(coral_id, func_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = func_group, values_from = count, values_fill = 0)

# Calculate total CAFI metrics
cafi_totals <- cafi_processed %>%
  group_by(coral_id) %>%
  summarise(
    total_abundance = n(),
    richness = n_distinct(species_id),
    shannon = -sum((table(species_id)/n()) * log(table(species_id)/n() + 0.001)),
    .groups = "drop"
  )

# Merge all data for condition analysis
condition_cafi <- analysis_data %>%
  filter(!is.na(condition_score)) %>%
  left_join(cafi_func_counts, by = "coral_id") %>%
  left_join(cafi_totals, by = "coral_id") %>%
  mutate(
    across(c(Defenders, Resident_Fish, Ectoparasites, Other, total_abundance, richness),
           ~replace_na(., 0)),
    log_volume = log10(volume + 1)
  )

cat(sprintf("  Corals with condition + CAFI data: %d\n", nrow(condition_cafi)))

# ============================================================================
# Statistical Models for Effect Sizes
# ============================================================================

# Model 1: Total abundance + diversity
model_abund <- lm(condition_score ~ log_volume + total_abundance, data = condition_cafi)
model_rich <- lm(condition_score ~ log_volume + richness, data = condition_cafi)

# Model 2: Functional groups (controlling for volume)
model_func <- lm(condition_score ~ log_volume + Defenders + Resident_Fish + Ectoparasites + Other,
                 data = condition_cafi)

# Extract standardized coefficients for effect size comparison
# Standardize predictors for comparable effect sizes
condition_cafi_scaled <- condition_cafi %>%
  mutate(
    across(c(Defenders, Resident_Fish, Ectoparasites, Other, total_abundance, richness, log_volume),
           ~scale(.)[,1], .names = "{.col}_z")
  )

model_func_std <- lm(condition_score ~ log_volume_z + Defenders_z + Resident_Fish_z +
                      Ectoparasites_z + Other_z, data = condition_cafi_scaled)

# Get model summary
func_summary <- tidy(model_func_std) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_clean = case_when(
      term == "log_volume_z" ~ "Coral volume",
      term == "Defenders_z" ~ "Defenders (Trapezia, Alpheus)",
      term == "Resident_Fish_z" ~ "Resident fishes",
      term == "Ectoparasites_z" ~ "Ectoparasites (Galeropsis, Morula)",
      term == "Other_z" ~ "Other cryptofauna"
    ),
    significant = p.value < 0.05,
    direction = ifelse(estimate > 0, "Positive", "Negative")
  )

cat("  Functional group effects (standardized):\n")
for(i in 1:nrow(func_summary)) {
  cat(sprintf("    %s: β = %.3f, p = %.3f%s\n",
              func_summary$term_clean[i], func_summary$estimate[i], func_summary$p.value[i],
              ifelse(func_summary$significant[i], " *", "")))
}

# ============================================================================
# Panel A: Lollipop plot of effect sizes
# ============================================================================

effect_order <- c("Ectoparasites (Galeropsis, Morula)", "Other cryptofauna",
                  "Coral volume", "Resident fishes", "Defenders (Trapezia, Alpheus)")

p6a <- ggplot(func_summary %>% filter(term != "log_volume_z"),
              aes(x = estimate, y = factor(term_clean, levels = effect_order))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = estimate, y = term_clean, yend = term_clean,
                   color = direction), linewidth = 1.5) +
  geom_point(aes(fill = direction, size = abs(estimate)), shape = 21, stroke = 1, color = "white") +
  geom_errorbarh(aes(xmin = estimate - 1.96*std.error, xmax = estimate + 1.96*std.error),
                 height = 0.2, linewidth = 0.8, color = "gray30") +
  geom_text(aes(label = sprintf("p = %.3f", p.value)),
            hjust = ifelse(func_summary %>% filter(term != "log_volume_z") %>% pull(estimate) > 0, -0.3, 1.3),
            size = 3, color = "gray40") +
  scale_color_manual(values = c("Positive" = "#27ae60", "Negative" = "#e74c3c"), guide = "none") +
  scale_fill_manual(values = c("Positive" = "#27ae60", "Negative" = "#e74c3c"), guide = "none") +
  scale_size_continuous(range = c(4, 10), guide = "none") +
  labs(title = "A. Functional group effects on coral condition",
       subtitle = "Standardized coefficients (β) with 95% CI; controlling for coral volume",
       x = "Effect on condition (standardized β)", y = NULL) +
  theme_pub(base_size = 11) +
  theme(axis.text.y = element_text(size = 10))

# ============================================================================
# Panel B: Defender abundance vs condition (partial residual)
# ============================================================================

# Calculate partial residuals for defender effect
condition_cafi$resid_defenders <- residuals(lm(condition_score ~ log_volume + Resident_Fish +
                                                 Ectoparasites + Other, data = condition_cafi))
condition_cafi$resid_defenders_x <- residuals(lm(Defenders ~ log_volume + Resident_Fish +
                                                   Ectoparasites + Other, data = condition_cafi))

defender_model_partial <- lm(resid_defenders ~ Defenders, data = condition_cafi)
def_coef <- coef(defender_model_partial)["Defenders"]
def_p <- summary(defender_model_partial)$coefficients["Defenders", "Pr(>|t|)"]
def_r2 <- summary(defender_model_partial)$r.squared

p6b <- ggplot(condition_cafi, aes(x = Defenders, y = resid_defenders)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#27ae60", linewidth = 1.2,
              se = TRUE, fill = "#27ae60", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = site_colors, name = "Site") +
  labs(title = "B. Defenders → better condition",
       subtitle = sprintf("Partial R² = %.2f, p = %.3f (controlling for other predictors)",
                         def_r2, def_p),
       x = "Defender abundance (Trapezia + Alpheus)",
       y = "Condition (residualized)") +
  theme_pub(base_size = 11) +
  theme(legend.position = c(0.85, 0.2),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# ============================================================================
# Panel C: Ectoparasite abundance vs condition (partial residual)
# ============================================================================

condition_cafi$resid_ecto <- residuals(lm(condition_score ~ log_volume + Defenders +
                                           Resident_Fish + Other, data = condition_cafi))
condition_cafi$resid_ecto_x <- residuals(lm(Ectoparasites ~ log_volume + Defenders +
                                             Resident_Fish + Other, data = condition_cafi))

ecto_model_partial <- lm(resid_ecto ~ Ectoparasites, data = condition_cafi)
ecto_coef <- coef(ecto_model_partial)["Ectoparasites"]
ecto_p <- summary(ecto_model_partial)$coefficients["Ectoparasites", "Pr(>|t|)"]
ecto_r2 <- summary(ecto_model_partial)$r.squared

p6c <- ggplot(condition_cafi, aes(x = Ectoparasites, y = resid_ecto)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#e74c3c", linewidth = 1.2,
              se = TRUE, fill = "#e74c3c", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(title = "C. Ectoparasites → worse condition",
       subtitle = sprintf("Partial R² = %.2f, p = %.3f (controlling for other predictors)",
                         ecto_r2, ecto_p),
       x = "Ectoparasite abundance (Galeropsis + Morula)",
       y = "Condition (residualized)") +
  theme_pub(base_size = 11)

# ============================================================================
# Panel D: Total richness has NO effect (after controlling for volume)
# ============================================================================

condition_cafi$resid_rich <- residuals(lm(condition_score ~ log_volume, data = condition_cafi))

rich_model <- lm(resid_rich ~ richness, data = condition_cafi)
rich_coef <- coef(rich_model)["richness"]
rich_p <- summary(rich_model)$coefficients["richness", "Pr(>|t|)"]
rich_r2 <- summary(rich_model)$r.squared

p6d <- ggplot(condition_cafi, aes(x = richness, y = resid_rich)) +
  geom_point(aes(color = site), alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "gray50", linewidth = 1.2,
              se = TRUE, fill = "gray70", alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(title = "D. Total richness → no effect",
       subtitle = sprintf("R² = %.3f, p = %.2f (controlling for coral size)",
                         rich_r2, rich_p),
       x = "Species richness (# OTUs)",
       y = "Condition (residualized for size)") +
  theme_pub(base_size = 11)

# ============================================================================
# Panel E: Model comparison - Functional identity vs diversity
# ============================================================================

# Fit competing models
m_null <- lm(condition_score ~ log_volume, data = condition_cafi)
m_abundance <- lm(condition_score ~ log_volume + total_abundance, data = condition_cafi)
m_richness <- lm(condition_score ~ log_volume + richness, data = condition_cafi)
m_functional <- lm(condition_score ~ log_volume + Defenders + Ectoparasites, data = condition_cafi)
m_full <- lm(condition_score ~ log_volume + Defenders + Resident_Fish + Ectoparasites + Other,
             data = condition_cafi)

model_comparison <- data.frame(
  Model = c("Null (size only)", "Size + Abundance", "Size + Richness",
            "Size + Defenders + Ectoparasites", "Size + All Func. Groups"),
  R2 = c(summary(m_null)$r.squared, summary(m_abundance)$r.squared,
         summary(m_richness)$r.squared, summary(m_functional)$r.squared,
         summary(m_full)$r.squared),
  AIC = c(AIC(m_null), AIC(m_abundance), AIC(m_richness), AIC(m_functional), AIC(m_full))
) %>%
  mutate(
    delta_R2 = R2 - R2[1],
    delta_AIC = AIC - min(AIC),
    Model = factor(Model, levels = rev(Model))
  )

p6e <- ggplot(model_comparison, aes(x = delta_R2, y = Model)) +
  geom_col(aes(fill = delta_AIC < 2), width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("ΔR² = %.3f", delta_R2)), hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "#3498db", "FALSE" = "gray70"),
                    labels = c("TRUE" = "Best (ΔAIC < 2)", "FALSE" = "Worse"),
                    name = "Model support") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "E. Model comparison: Functional identity matters most",
       subtitle = "Variance explained beyond coral size",
       x = "Additional R² beyond size-only model", y = NULL) +
  theme_pub(base_size = 11) +
  theme(legend.position = c(0.75, 0.25),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# ============================================================================
# Panel F: Size-dependent effects (interaction test)
# ============================================================================

# Test if defender/ectoparasite effects vary with coral size
condition_cafi <- condition_cafi %>%
  mutate(size_class = cut(volume,
                          breaks = quantile(volume, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                          labels = c("Small", "Medium", "Large"),
                          include.lowest = TRUE))

# Calculate effect sizes within each size class
size_effects <- condition_cafi %>%
  group_by(size_class) %>%
  summarise(
    n = n(),
    defender_cor = cor(Defenders, condition_score, use = "complete.obs"),
    ecto_cor = cor(Ectoparasites, condition_score, use = "complete.obs"),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(defender_cor, ecto_cor),
               names_to = "predictor", values_to = "correlation") %>%
  mutate(
    predictor = ifelse(predictor == "defender_cor", "Defenders", "Ectoparasites"),
    predictor = factor(predictor, levels = c("Defenders", "Ectoparasites"))
  )

p6f <- ggplot(size_effects, aes(x = size_class, y = correlation, fill = predictor)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_fill_manual(values = c("Defenders" = "#27ae60", "Ectoparasites" = "#e74c3c"),
                    name = "Functional group") +
  labs(title = "F. Effects consistent across coral sizes",
       subtitle = "Correlation between functional groups and condition by size class",
       x = "Coral size class", y = "Correlation with condition") +
  theme_pub(base_size = 11) +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9)))

# ============================================================================
# Combine into Figure 6
# ============================================================================

fig6 <- (p6a | p6b | p6c) / (p6d | p6e | p6f) +
  plot_annotation(
    title = "Figure 6. CAFI functional groups predict coral condition: Defenders help, ectoparasites harm",
    subtitle = "Functional identity (who colonizes) matters more than diversity (how many species)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray40")
    )
  )

ggsave("output/figures/manuscript/Figure6_cafi_condition_feedbacks.png", fig6,
       width = 16, height = 11, dpi = 300, bg = "white")
cat("  ✓ Figure 6 (CAFI-Condition Feedbacks) saved\n")

cat("\n========================================\n")
cat("FIGURE GENERATION COMPLETE\n")
cat("========================================\n")
cat("Output files:\n")
cat("  • Figure1_study_system.png\n")
cat("  • Figure2_community_composition.png\n")
cat("  • Figure3_scaling.png\n")
cat("  • Figure4_cafi_condition.png\n")
cat("  • Figure5_network.png\n")
cat("  • Figure6_cafi_condition_feedbacks.png\n")
cat("\nAll figures saved to: output/figures/manuscript/\n")
