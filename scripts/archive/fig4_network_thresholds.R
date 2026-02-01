# ============================================================================
# fig4_network_thresholds.R - Compare network at different correlation thresholds
# ============================================================================

cat("\n============================================================\n")
cat("NETWORK VISUALIZATIONS AT DIFFERENT THRESHOLDS\n")
cat("============================================================\n\n")

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

network_results <- load_object("cafi_network")
centrality_df <- network_results$centrality
edge_list <- network_results$edge_list

library(igraph)
library(ggrepel)
library(patchwork)

fig_dir <- file.path(PATHS$figures, "06_network")

# Color palette
guild_colors <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73", "4" = "#CC79A7")
type_colors <- c(
  "crab" = "#E69F00", "shrimp" = "#0072B2", "fish" = "#009E73",
  "echinoderm" = "#999999", "worm" = "#56B4E9", "snail" = "#F0E442",
  "hermit" = "#CC79A7", "squat_lobster" = "#D55E00", "amphipod" = "#666666"
)

# Function to create network plot at given threshold
make_network_plot <- function(edge_list, centrality_df, threshold, title) {

  # Filter edges
  edges_filtered <- edge_list %>% filter(correlation > threshold)

  if (nrow(edges_filtered) == 0) {
    return(ggplot() + labs(title = paste(title, "- No edges")))
  }

  # Get species in filtered network
  species_in <- unique(c(edges_filtered$sp1, edges_filtered$sp2))

  # Build graph
  g <- graph_from_data_frame(
    edges_filtered %>% dplyr::select(sp1, sp2, correlation) %>%
      rename(from = sp1, to = sp2, weight = correlation),
    directed = FALSE,
    vertices = centrality_df %>%
      filter(species %in% species_in) %>%
      dplyr::select(species, type, module, degree) %>%
      rename(name = species)
  )

  # Layout
  set.seed(42)
  layout <- layout_with_fr(g, weights = E(g)$weight)

  # Node data
  node_df <- data.frame(
    x = layout[, 1],
    y = layout[, 2],
    species = V(g)$name,
    module = as.character(V(g)$module),
    type = V(g)$type,
    degree = igraph::degree(g)
  )

  # Edge data
  el <- as_edgelist(g, names = FALSE)
  edge_df <- data.frame(
    x = layout[el[,1], 1], y = layout[el[,1], 2],
    xend = layout[el[,2], 1], yend = layout[el[,2], 2],
    weight = E(g)$weight
  )

  # Labels for top species
  top_sp <- node_df %>% slice_max(degree, n = min(10, nrow(node_df)))

  n_sp <- nrow(node_df)
  n_edges <- nrow(edge_df)
  density <- round(n_edges / (n_sp * (n_sp - 1) / 2), 2)

  p <- ggplot() +
    geom_segment(data = edge_df,
                 aes(x = x, y = y, xend = xend, yend = yend, alpha = weight),
                 color = "gray50", linewidth = 0.4) +
    geom_point(data = node_df,
               aes(x = x, y = y, fill = module, size = degree),
               shape = 21, color = "white", stroke = 0.5) +
    geom_text_repel(data = top_sp,
                    aes(x = x, y = y, label = species),
                    size = 2.5, max.overlaps = 15, segment.alpha = 0.5,
                    bg.color = "white", bg.r = 0.1) +
    scale_fill_manual(values = guild_colors, name = "Guild") +
    scale_size_continuous(range = c(2, 10), name = "Degree") +
    scale_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
    labs(title = title,
         subtitle = sprintf("%d species, %d edges, density = %.0f%%", n_sp, n_edges, density * 100)) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9, color = "gray40"),
          legend.position = "none")

  return(p)
}

# Create plots at different thresholds
cat("Creating network at r > 0.3 (current)...\n")
p_03 <- make_network_plot(edge_list, centrality_df, 0.3, "A. r > 0.3 (current)")

cat("Creating network at r > 0.5...\n")
p_05 <- make_network_plot(edge_list, centrality_df, 0.5, "B. r > 0.5")

cat("Creating network at r > 0.6...\n")
p_06 <- make_network_plot(edge_list, centrality_df, 0.6, "C. r > 0.6")

cat("Creating network at r > 0.7...\n")
p_07 <- make_network_plot(edge_list, centrality_df, 0.7, "D. r > 0.7")

# Combine
fig_thresholds <- (p_03 | p_05) / (p_06 | p_07) +
  plot_annotation(
    title = "Network Visualization at Different Correlation Thresholds",
    subtitle = "Higher thresholds show only strongest associations, making guild structure clearer",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(fig_dir, "network_threshold_comparison.png"), fig_thresholds,
       width = 14, height = 12, dpi = 300, bg = "white")
cat("\nSaved: network_threshold_comparison.png\n")

# Also create standalone versions of the cleaner networks
cat("\nCreating standalone r > 0.6 network...\n")
p_06_big <- make_network_plot(edge_list, centrality_df, 0.6, "Co-occurrence Network (r > 0.6)") +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "network_r06.png"), p_06_big,
       width = 12, height = 10, dpi = 300, bg = "white")

cat("Creating standalone r > 0.7 network...\n")
p_07_big <- make_network_plot(edge_list, centrality_df, 0.7, "Co-occurrence Network (r > 0.7)") +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "network_r07.png"), p_07_big,
       width = 12, height = 10, dpi = 300, bg = "white")

# ============================================================================
# NETWORK COLORED BY TAXONOMIC TYPE (at r > 0.6)
# ============================================================================

cat("\nCreating network colored by taxonomic type (r > 0.6)...\n")

edges_06 <- edge_list %>% filter(correlation > 0.6)
species_in <- unique(c(edges_06$sp1, edges_06$sp2))

g <- graph_from_data_frame(
  edges_06 %>% dplyr::select(sp1, sp2, correlation) %>%
    rename(from = sp1, to = sp2, weight = correlation),
  directed = FALSE,
  vertices = centrality_df %>%
    filter(species %in% species_in) %>%
    dplyr::select(species, type, module, degree) %>%
    rename(name = species)
)

set.seed(42)
layout <- layout_with_fr(g, weights = E(g)$weight)

node_df <- data.frame(
  x = layout[, 1], y = layout[, 2],
  species = V(g)$name,
  module = as.character(V(g)$module),
  type = V(g)$type,
  degree = igraph::degree(g)
)

el <- as_edgelist(g, names = FALSE)
edge_df <- data.frame(
  x = layout[el[,1], 1], y = layout[el[,1], 2],
  xend = layout[el[,2], 1], yend = layout[el[,2], 2],
  weight = E(g)$weight
)

top_sp <- node_df %>% slice_max(degree, n = 12)

p_type <- ggplot() +
  geom_segment(data = edge_df,
               aes(x = x, y = y, xend = xend, yend = yend, alpha = weight),
               color = "gray50", linewidth = 0.5) +
  geom_point(data = node_df,
             aes(x = x, y = y, fill = type, size = degree),
             shape = 21, color = "white", stroke = 0.6) +
  geom_text_repel(data = top_sp,
                  aes(x = x, y = y, label = species),
                  size = 3, max.overlaps = 20, segment.alpha = 0.5,
                  bg.color = "white", bg.r = 0.15) +
  scale_fill_manual(values = type_colors, name = "Taxa") +
  scale_size_continuous(range = c(3, 12), name = "Degree") +
  scale_alpha_continuous(range = c(0.3, 0.9), guide = "none") +
  labs(title = "Co-occurrence Network by Taxonomic Group",
       subtitle = sprintf("r > 0.6 | %d species, %d edges | Node size = connectivity",
                          nrow(node_df), nrow(edge_df))) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "right")

ggsave(file.path(fig_dir, "network_by_taxa_r06.png"), p_type,
       width = 14, height = 10, dpi = 300, bg = "white")

# ============================================================================
# COMBINED: NETWORK + GUILD SUMMARY
# ============================================================================

cat("\nCreating combined network + guild summary...\n")

# Guild summary for filtered network
guild_summary <- node_df %>%
  group_by(module) %>%
  summarise(
    n = n(),
    taxa = paste(names(sort(table(type), decreasing = TRUE))[1:min(2, n_distinct(type))], collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(guild = paste0("Guild ", module, "\n(n=", n, ")\n", taxa))

p_network_clean <- ggplot() +
  geom_segment(data = edge_df,
               aes(x = x, y = y, xend = xend, yend = yend, alpha = weight),
               color = "gray50", linewidth = 0.5) +
  geom_point(data = node_df,
             aes(x = x, y = y, fill = module, size = degree),
             shape = 21, color = "white", stroke = 0.6) +
  geom_text_repel(data = top_sp,
                  aes(x = x, y = y, label = species),
                  size = 3, max.overlaps = 20, segment.alpha = 0.5,
                  bg.color = "white", bg.r = 0.15) +
  scale_fill_manual(values = guild_colors,
                    labels = guild_summary$guild,
                    name = "Guild") +
  scale_size_continuous(range = c(3, 12), name = "Degree") +
  scale_alpha_continuous(range = c(0.3, 0.9), guide = "none") +
  labs(title = "CAFI Co-occurrence Network",
       subtitle = sprintf("Strongest associations only (r > 0.6) | %d species, %d edges",
                          nrow(node_df), nrow(edge_df))) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        legend.position = "right",
        legend.text = element_text(size = 9))

ggsave(file.path(fig_dir, "network_clean_r06.png"), p_network_clean,
       width = 14, height = 10, dpi = 300, bg = "white")

cat("\n============================================================\n")
cat("DONE - Networks saved to output/figures/06_network/\n")
cat("============================================================\n")
cat("\nFiles created:\n")
cat("  - network_threshold_comparison.png (4-panel comparison)\n")
cat("  - network_r06.png (r > 0.6, by guild)\n")
cat("  - network_r07.png (r > 0.7, by guild)\n")
cat("  - network_by_taxa_r06.png (r > 0.6, by taxonomic type)\n")
cat("  - network_clean_r06.png (r > 0.6, with guild descriptions)\n")
