# ============================================================================
# fig4_guild_heatmap.R - Ecological Guild Structure Visualization (Figure 4)
# ============================================================================
#
# PURPOSE: Create a clear visualization of CAFI species guild structure using
#          heatmaps and clean summary plots instead of dense network diagrams.
#
# RATIONALE: The co-occurrence network has 58 species with 1,081 edges (65%
#            density), making traditional network visualization uninformative.
#            This approach uses:
#            - Panel A: Co-occurrence heatmap ordered by guild (shows modular structure)
#            - Panel B: Guild composition treemap/bar (what's in each guild)
#            - Panel C: Strongest species associations (top ecological signals)
#
# OUTPUTS:
#   - output/figures/manuscript/fig4_ecological_guilds.png (main figure)
#   - output/figures/supplement/figS_guild_details.png (extended version)
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-01-28
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    FIGURE 4: ECOLOGICAL GUILD STRUCTURE\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

# Load setup and data
if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

# Load network results
network_results <- tryCatch(
  load_object("cafi_network"),
  error = function(e) {
    cat("  [WARNING] Network results not found. Running network analysis first.\n")
    source(here::here("scripts/06_network_analysis.R"))
    load_object("cafi_network")
  }
)

# Additional packages
if (!require("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = "https://cran.rstudio.com")
}
library(pheatmap)

if (!require("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer", repos = "https://cran.rstudio.com")
}
library(RColorBrewer)

if (!require("treemapify", quietly = TRUE)) {
  install.packages("treemapify", repos = "https://cran.rstudio.com")
}
library(treemapify)

# Create output directories
fig_dir_manuscript <- file.path(PATHS$figures, "manuscript")
fig_dir_supplement <- file.path(PATHS$figures, "supplement")
dir.create(fig_dir_manuscript, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir_supplement, recursive = TRUE, showWarnings = FALSE)

cat("[OK] Setup complete\n\n")

# ============================================================================
# EXTRACT DATA FROM NETWORK RESULTS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("EXTRACTING DATA\n")
cat("------------------------------------------------------------\n\n")

# Extract key objects
g <- network_results$graph
centrality_df <- network_results$centrality
module_summary <- network_results$module_summary
edge_list <- network_results$edge_list
type_colors <- network_results$type_colors

# Get species-level info with modules
species_info <- centrality_df %>%
  dplyr::select(species, type, functional_group, module, degree, abundance, occurrence,
                eigenvector, hub_score) %>%
  arrange(module, desc(degree))

n_species <- nrow(species_info)
n_modules <- length(unique(species_info$module))

cat("  Species in network:", n_species, "\n")
cat("  Number of modules:", n_modules, "\n")
cat("  Total edges:", nrow(edge_list), "\n\n")

# ============================================================================
# PANEL A: CO-OCCURRENCE HEATMAP BY GUILD
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PANEL A: CO-OCCURRENCE HEATMAP\n")
cat("------------------------------------------------------------\n\n")

# Create correlation matrix from edge list
all_species <- unique(c(edge_list$sp1, edge_list$sp2))
n_sp <- length(all_species)

# Initialize correlation matrix
cor_mat <- matrix(0, nrow = n_sp, ncol = n_sp)
rownames(cor_mat) <- all_species
colnames(cor_mat) <- all_species

# Fill in correlations from edge list
for (i in 1:nrow(edge_list)) {
  sp1 <- edge_list$sp1[i]
  sp2 <- edge_list$sp2[i]
  cor_val <- edge_list$correlation[i]
  cor_mat[sp1, sp2] <- cor_val
  cor_mat[sp2, sp1] <- cor_val
}

# Set diagonal to 1
diag(cor_mat) <- 1

# Order species by module (to reveal block structure)
species_order <- species_info %>%
  arrange(module, desc(degree)) %>%
  pull(species)

# Filter to species in our matrix
species_order <- species_order[species_order %in% rownames(cor_mat)]

# Reorder correlation matrix
cor_mat_ordered <- cor_mat[species_order, species_order]

# Create annotation data for modules
module_annotation <- species_info %>%
  filter(species %in% species_order) %>%
  dplyr::select(species, module, type) %>%
  distinct() %>%
  arrange(match(species, species_order))

rownames(module_annotation) <- module_annotation$species
module_annotation$module <- factor(module_annotation$module)
module_annotation$Guild <- paste0("Guild ", module_annotation$module)
annotation_df <- module_annotation %>%
  dplyr::select(Guild, type) %>%
  rename(Type = type)

# Define annotation colors
guild_colors <- c(
  "Guild 1" = "#E69F00",  # Orange
  "Guild 2" = "#56B4E9",  # Sky blue
  "Guild 3" = "#009E73",  # Teal
  "Guild 4" = "#D55E00"   # Vermillion
)

type_colors_expanded <- c(
  "crab" = "#D55E00",
  "shrimp" = "#0072B2",
  "fish" = "#009E73",
  "snail" = "#E69F00",
  "hermit" = "#CC79A7",
  "echinoderm" = "#999999",
  "worm" = "#56B4E9",
  "squat_lobster" = "#F0E442",
  "amphipod" = "#999999",
  "unknown" = "#999999"
)

annotation_colors <- list(
  Guild = guild_colors[1:n_modules],
  Type = type_colors_expanded
)

cat("  Creating heatmap with", length(species_order), "species\n")

# Create custom color palette (blue-white-red for correlations)
heatmap_colors <- colorRampPalette(c("#0072B2", "white", "#D55E00"))(100)

# Generate heatmap (save as grob for later combination)
png(file.path(fig_dir_manuscript, "fig4_panel_a_heatmap.png"),
    width = 10, height = 9, units = "in", res = 300, bg = "white")

pheatmap(
  cor_mat_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = heatmap_colors,
  breaks = seq(0, 1, length.out = 101),
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "A. Species Co-occurrence by Guild",
  fontsize = 10,
  fontsize_row = 6,
  fontsize_col = 6,
  border_color = NA,
  legend_breaks = c(0, 0.3, 0.5, 0.7, 1),
  legend_labels = c("0", "0.3", "0.5", "0.7", "1"),
  annotation_legend = TRUE,
  annotation_names_row = TRUE,
  annotation_names_col = FALSE
)

dev.off()

cat("  Saved: fig4_panel_a_heatmap.png\n\n")

# ============================================================================
# PANEL B: GUILD COMPOSITION
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PANEL B: GUILD COMPOSITION\n")
cat("------------------------------------------------------------\n\n")

# Calculate guild composition
guild_composition <- species_info %>%
  group_by(module) %>%
  summarise(
    n_species = n(),
    total_abundance = sum(abundance),
    types = paste(unique(type), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    guild_label = paste0("Guild ", module),
    pct_species = round(100 * n_species / sum(n_species), 1)
  )

cat("  Guild composition:\n")
print(guild_composition)
cat("\n")

# Guild by type breakdown
guild_by_type <- species_info %>%
  group_by(module, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(guild_label = paste0("Guild ", module))

# Create bar chart of guild composition by type
p_guild_comp <- ggplot(guild_by_type, aes(x = guild_label, y = n, fill = type)) +
  geom_col(width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = type_colors_expanded, name = "Taxa") +
  labs(
    x = NULL,
    y = "Number of Species",
    title = "B. Guild Composition by Taxa"
  ) +
  theme_publication(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold")
  ) +
  coord_flip()

ggsave(file.path(fig_dir_manuscript, "fig4_panel_b_guild_comp.png"),
       p_guild_comp, width = 6, height = 4, dpi = 300, bg = "white")
cat("  Saved: fig4_panel_b_guild_comp.png\n\n")

# ============================================================================
# PANEL C: STRONGEST SPECIES ASSOCIATIONS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PANEL C: STRONGEST ASSOCIATIONS\n")
cat("------------------------------------------------------------\n\n")

# Get top associations
top_associations <- edge_list %>%
  arrange(desc(correlation)) %>%
  head(20) %>%
  mutate(
    # Get guild for each species
    guild1 = species_info$module[match(sp1, species_info$species)],
    guild2 = species_info$module[match(sp2, species_info$species)],
    type1 = species_info$type[match(sp1, species_info$species)],
    type2 = species_info$type[match(sp2, species_info$species)],
    # Create pair label
    pair_label = paste0(sp1, " - ", sp2),
    # Is it within or between guild?
    association_type = ifelse(guild1 == guild2, "Within guild", "Between guild")
  )

cat("  Top 10 strongest associations:\n")
print(top_associations %>% head(10) %>% dplyr::select(sp1, sp2, correlation, association_type))
cat("\n")

# Create lollipop plot
p_top_assoc <- top_associations %>%
  head(15) %>%
  ggplot(aes(x = reorder(pair_label, correlation), y = correlation, color = association_type)) +
  geom_segment(aes(xend = pair_label, y = 0, yend = correlation), linewidth = 1) +
  geom_point(size = 4) +
  coord_flip() +
  scale_color_manual(
    values = c("Within guild" = "#0072B2", "Between guild" = "#D55E00"),
    name = "Association"
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  labs(
    x = NULL,
    y = "Co-occurrence Correlation",
    title = "C. Strongest Species Associations"
  ) +
  theme_publication(base_size = 10) +
  theme(
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(fill = "white", color = "gray80"),
    axis.text.y = element_text(size = 8)
  )

ggsave(file.path(fig_dir_manuscript, "fig4_panel_c_associations.png"),
       p_top_assoc, width = 7, height = 5, dpi = 300, bg = "white")
cat("  Saved: fig4_panel_c_associations.png\n\n")

# ============================================================================
# PANEL D: HUB SPECIES BY GUILD (simplified network or bar)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("PANEL D: HUB SPECIES\n")
cat("------------------------------------------------------------\n\n")

# Get top hub species per guild
hub_by_guild <- species_info %>%
  group_by(module) %>%
  slice_max(order_by = hub_score, n = 3) %>%
  ungroup() %>%
  mutate(
    guild_label = paste0("Guild ", module),
    # Shorten species names if too long
    species_short = ifelse(nchar(species) > 20,
                          paste0(substr(species, 1, 18), "..."),
                          species)
  )

cat("  Top hub species by guild:\n")
print(hub_by_guild %>% dplyr::select(guild_label, species, type, degree, hub_score))
cat("\n")

# Create hub species plot
p_hubs <- hub_by_guild %>%
  ggplot(aes(x = reorder(species_short, hub_score), y = hub_score, fill = guild_label)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_text(aes(label = paste0("d=", degree)),
            hjust = -0.1, size = 2.8, color = "gray30") +
  coord_flip() +
  facet_wrap(~guild_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = guild_colors, guide = "none") +
  labs(
    x = NULL,
    y = "Hub Score (connectivity + influence)",
    title = "D. Key Species by Guild"
  ) +
  theme_publication(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8, face = "italic")
  )

ggsave(file.path(fig_dir_manuscript, "fig4_panel_d_hubs.png"),
       p_hubs, width = 8, height = 5, dpi = 300, bg = "white")
cat("  Saved: fig4_panel_d_hubs.png\n\n")

# ============================================================================
# COMBINED FIGURE (3-PANEL LAYOUT)
# ============================================================================

cat("------------------------------------------------------------\n")
cat("CREATING COMBINED FIGURE\n")
cat("------------------------------------------------------------\n\n")

# For the combined figure, we'll use a modified approach
# Panel A (heatmap) needs special handling since pheatmap doesn't return a ggplot

# Create a ggplot version of the heatmap for patchwork compatibility
cor_mat_long <- as.data.frame(cor_mat_ordered) %>%
  mutate(species1 = rownames(cor_mat_ordered)) %>%
  pivot_longer(-species1, names_to = "species2", values_to = "correlation") %>%
  mutate(
    # Add guild info
    guild1 = species_info$module[match(species1, species_info$species)],
    guild2 = species_info$module[match(species2, species_info$species)],
    # Order factors by guild
    species1 = factor(species1, levels = species_order),
    species2 = factor(species2, levels = species_order)
  )

# Calculate guild information for annotations
guild_counts <- table(species_info$module[match(species_order, species_info$species)])
guild_breaks <- cumsum(guild_counts)
guild_centers <- c(0, guild_breaks[-length(guild_breaks)]) + diff(c(0, guild_breaks))/2

# Create guild labels with species counts
guild_labels <- paste0("Guild ", 1:n_modules, "\n(n=", as.numeric(guild_counts), ")")

# Create ggplot heatmap with clear axes
p_heatmap <- ggplot(cor_mat_long, aes(x = species2, y = species1, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#0072B2",
    mid = "white",
    high = "#D55E00",
    midpoint = 0.5,
    limits = c(0, 1),
    name = "Co-occurrence\nCorrelation"
  ) +
  # Add guild separator lines (more prominent)
  geom_vline(xintercept = guild_breaks[-n_modules] + 0.5,
             color = "black", linewidth = 1.2) +
  geom_hline(yintercept = guild_breaks[-n_modules] + 0.5,
             color = "black", linewidth = 1.2) +
  # Add guild labels at appropriate positions
  scale_x_discrete(
    breaks = species_order[round(guild_centers)],
    labels = guild_labels,
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    breaks = species_order[round(guild_centers)],
    labels = guild_labels,
    expand = c(0, 0)
  ) +
  labs(
    x = "Species (ordered by guild membership)",
    y = "Species (ordered by guild membership)",
    title = "A. Species Co-occurrence Matrix by Guild"
  ) +
  theme_publication(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 1, face = "bold"),
    axis.text.y = element_text(size = 8, hjust = 1, vjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 8)),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold")
  ) +
  coord_fixed()

# ============================================================================
# FINAL COMBINED FIGURE - 3 Panel Version
# ============================================================================

# Create cleaner version of guild composition (stacked percentage)
guild_by_type_pct <- guild_by_type %>%
  group_by(guild_label) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup()

p_guild_comp_clean <- ggplot(guild_by_type_pct,
                              aes(x = guild_label, y = pct, fill = type)) +
  geom_col(position = "fill", width = 0.75, alpha = 0.9) +
  scale_fill_manual(values = type_colors_expanded, name = "Taxa") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Proportion of Species",
    title = "B. Guild Taxonomic Composition"
  ) +
  theme_publication(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(face = "bold", size = 10)
  ) +
  guides(fill = guide_legend(nrow = 2))

# Cleaner top associations panel (only top 12)
p_top_assoc_clean <- top_associations %>%
  head(12) %>%
  mutate(
    # Create shorter labels
    sp1_short = gsub(" ", "\n", sp1),
    sp2_short = gsub(" ", "\n", sp2),
    pair_short = paste0(
      ifelse(nchar(sp1) > 12, paste0(substr(sp1, 1, 10), ".."), sp1),
      " - ",
      ifelse(nchar(sp2) > 12, paste0(substr(sp2, 1, 10), ".."), sp2)
    )
  ) %>%
  ggplot(aes(x = reorder(pair_short, correlation), y = correlation,
             fill = association_type)) +
  geom_col(width = 0.7, alpha = 0.9) +
  coord_flip() +
  scale_fill_manual(
    values = c("Within guild" = "#0072B2", "Between guild" = "#D55E00"),
    name = NULL
  ) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
  labs(
    x = NULL,
    y = "Correlation Strength",
    title = "C. Strongest Species Pairs"
  ) +
  theme_publication(base_size = 10) +
  theme(
    legend.position = c(0.75, 0.15),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
    legend.text = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  )

# ============================================================================
# ASSEMBLE FINAL FIGURE - 3 Panel Layout
# ============================================================================

# Layout: Large heatmap on left, two smaller panels stacked on right
layout <- "
AABB
AACC
"

p_combined <- p_heatmap +
  p_guild_comp_clean +
  p_top_assoc_clean +
  plot_layout(design = layout, heights = c(1, 1)) +
  plot_annotation(
    title = "Figure 4. CAFI Species Guild Structure",
    subtitle = sprintf("%d species form %d distinct ecological guilds based on co-occurrence patterns (n = %d coral hosts)",
                      n_species, n_modules, nrow(community_matrix)),
    caption = "Co-occurrence measured as volume-residualized Spearman correlation (r > 0.3, FDR p < 0.05).\nGuilds identified via Louvain community detection.",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 10, hjust = 0),
      plot.caption = element_text(size = 8, hjust = 0, color = "gray40")
    )
  )

ggsave(file.path(fig_dir_manuscript, "fig4_ecological_guilds.png"),
       p_combined, width = 12, height = 8, dpi = 300, bg = "white")
cat("  Saved: fig4_ecological_guilds.png\n\n")

# ============================================================================
# SUPPLEMENTARY FIGURE - EXTENDED GUILD DETAILS
# ============================================================================

cat("------------------------------------------------------------\n")
cat("SUPPLEMENTARY FIGURE: EXTENDED GUILD DETAILS\n")
cat("------------------------------------------------------------\n\n")

# S1: All species by guild
all_species_by_guild <- species_info %>%
  mutate(guild_label = paste0("Guild ", module)) %>%
  group_by(guild_label) %>%
  slice_max(order_by = degree, n = 15) %>%
  ungroup()

p_all_species <- ggplot(all_species_by_guild,
                        aes(x = reorder(species, degree), y = degree, fill = type)) +
  geom_col(width = 0.7, alpha = 0.9) +
  coord_flip() +
  facet_wrap(~guild_label, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = type_colors_expanded, name = "Taxa") +
  labs(
    x = NULL,
    y = "Degree (number of associations)",
    title = "Top Species by Guild (ranked by connectivity)"
  ) +
  theme_publication(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 7, face = "italic"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 2))

# S2: Guild connectivity matrix
guild_edges <- edge_list %>%
  mutate(
    guild1 = species_info$module[match(sp1, species_info$species)],
    guild2 = species_info$module[match(sp2, species_info$species)]
  ) %>%
  count(guild1, guild2) %>%
  mutate(
    guild1 = paste0("Guild ", guild1),
    guild2 = paste0("Guild ", guild2)
  )

# Make symmetric
guild_mat <- guild_edges %>%
  pivot_wider(names_from = guild2, values_from = n, values_fill = 0)

p_guild_matrix <- guild_edges %>%
  ggplot(aes(x = guild1, y = guild2, fill = n)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = n), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#0072B2", name = "Edges") +
  labs(
    x = NULL,
    y = NULL,
    title = "Inter-Guild Connectivity",
    subtitle = "Number of significant co-occurrence edges"
  ) +
  theme_publication(base_size = 11) +
  theme(
    axis.text = element_text(face = "bold"),
    panel.grid = element_blank()
  ) +
  coord_fixed()

# S3: Degree distribution by guild
p_degree_dist <- ggplot(species_info,
                        aes(x = degree, fill = factor(module))) +
  geom_histogram(binwidth = 2, alpha = 0.7, position = "identity") +
  facet_wrap(~paste0("Guild ", module), ncol = 2) +
  scale_fill_manual(values = unname(guild_colors[1:n_modules]), guide = "none") +
  labs(
    x = "Degree (number of associations)",
    y = "Number of species",
    title = "Degree Distribution by Guild"
  ) +
  theme_publication(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

# S4: Guild summary statistics
guild_stats <- species_info %>%
  group_by(module) %>%
  summarise(
    n_species = n(),
    mean_degree = round(mean(degree), 1),
    max_degree = max(degree),
    total_abundance = sum(abundance),
    dominant_type = names(which.max(table(type))),
    hub_species = species[which.max(hub_score)],
    .groups = "drop"
  ) %>%
  mutate(Guild = paste0("Guild ", module))

# Convert to long format for table display (convert all to character first)
guild_table_plot <- guild_stats %>%
  dplyr::select(Guild, n_species, mean_degree, dominant_type, hub_species) %>%
  mutate(
    n_species = as.character(n_species),
    mean_degree = as.character(mean_degree)
  ) %>%
  rename(
    `Species` = n_species,
    `Mean Degree` = mean_degree,
    `Dominant Taxa` = dominant_type,
    `Hub Species` = hub_species
  ) %>%
  pivot_longer(-Guild, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Guild, y = Metric)) +
  geom_tile(fill = "gray95", color = "white") +
  geom_text(aes(label = Value), size = 3) +
  labs(
    x = NULL,
    y = NULL,
    title = "Guild Summary Statistics"
  ) +
  theme_publication(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.text.y = element_text(face = "bold", size = 9)
  )

# Combine supplementary panels
p_supplement <- (p_all_species | (p_guild_matrix / p_degree_dist)) +
  plot_annotation(
    title = "Figure S4. Extended Guild Analysis",
    subtitle = "Detailed breakdown of CAFI ecological guild structure",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
  )

ggsave(file.path(fig_dir_supplement, "figS4_guild_details.png"),
       p_supplement, width = 14, height = 10, dpi = 300, bg = "white")
cat("  Saved: figS4_guild_details.png\n\n")

# ============================================================================
# SAVE SUMMARY TABLE
# ============================================================================

cat("------------------------------------------------------------\n")
cat("SAVING SUMMARY TABLES\n")
cat("------------------------------------------------------------\n\n")

# Guild summary with ecological interpretation
guild_interpretation <- guild_stats %>%
  mutate(
    interpretation = case_when(
      dominant_type == "shrimp" ~ "Shrimp-dominated: mobile cleaners/commensals",
      dominant_type == "crab" ~ "Crab-dominated: territorial residents",
      dominant_type == "worm" ~ "Cryptofauna: sessile/tube-dwelling",
      dominant_type == "fish" ~ "Vertebrate associates",
      TRUE ~ "Mixed assemblage"
    )
  )

save_table(guild_interpretation, "guild_ecological_summary")
cat("  Saved: guild_ecological_summary.csv\n")

# Top associations table
top_assoc_table <- top_associations %>%
  dplyr::select(sp1, sp2, correlation, type1, type2, guild1, guild2, association_type) %>%
  rename(
    species_1 = sp1,
    species_2 = sp2,
    type_1 = type1,
    type_2 = type2,
    guild_1 = guild1,
    guild_2 = guild2
  )

save_table(top_assoc_table, "top_species_associations")
cat("  Saved: top_species_associations.csv\n\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    FIGURE 4 GENERATION COMPLETE\n")
cat("============================================================\n\n")

cat("Key findings visualized:\n")
cat("  1. Species cluster into", n_modules, "distinct ecological guilds\n")
cat("  2. Guild structure visible as block-diagonal pattern in heatmap\n")
cat("  3. Most associations are within-guild (taxonomic clustering)\n")
cat("  4. Hub species identified for each guild\n\n")

cat("Files generated:\n")
cat("  Main figure:\n")
cat("    - output/figures/manuscript/fig4_ecological_guilds.png\n")
cat("  Supplementary:\n")
cat("    - output/figures/supplement/figS4_guild_details.png\n")
cat("  Individual panels:\n")
cat("    - output/figures/manuscript/fig4_panel_a_heatmap.png\n")
cat("    - output/figures/manuscript/fig4_panel_b_guild_comp.png\n")
cat("    - output/figures/manuscript/fig4_panel_c_associations.png\n")
cat("    - output/figures/manuscript/fig4_panel_d_hubs.png\n")
cat("  Tables:\n")
cat("    - output/tables/guild_ecological_summary.csv\n")
cat("    - output/tables/top_species_associations.csv\n\n")

cat("============================================================\n\n")
