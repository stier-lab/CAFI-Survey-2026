# ============================================================================
# fig4_all_options.R - Generate multiple figure options for guild analysis
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    GENERATING ALL FIGURE OPTIONS FOR GUILD ANALYSIS\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP
# ============================================================================

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))
if (!exists("coral_master")) source(here::here("scripts/01_load_data.R"))

network_results <- tryCatch(
  load_object("cafi_network"),
  error = function(e) {
    source(here::here("scripts/06_network_analysis.R"))
    load_object("cafi_network")
  }
)

library(patchwork)
library(ggrepel)
library(igraph)

# Output directory for all options
fig_dir <- file.path(PATHS$figures, "guild_options")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("[OK] Setup complete\n")
cat("Output directory:", fig_dir, "\n\n")

# ============================================================================
# EXTRACT AND PREPARE DATA
# ============================================================================

g <- network_results$graph
centrality_df <- network_results$centrality
edge_list <- network_results$edge_list

# Species info
species_info <- centrality_df %>%
  dplyr::select(species, type, module, degree, hub_score, occurrence) %>%
  arrange(module, desc(degree))

n_species <- nrow(species_info)
n_modules <- length(unique(species_info$module))

# Build correlation matrix
all_species <- unique(c(edge_list$sp1, edge_list$sp2))
cor_mat <- matrix(0, nrow = length(all_species), ncol = length(all_species))
rownames(cor_mat) <- colnames(cor_mat) <- all_species
for (i in 1:nrow(edge_list)) {
  cor_mat[edge_list$sp1[i], edge_list$sp2[i]] <- edge_list$correlation[i]
  cor_mat[edge_list$sp2[i], edge_list$sp1[i]] <- edge_list$correlation[i]
}
diag(cor_mat) <- 1

# Guild-sorted order
species_guild <- species_info %>%
  filter(species %in% rownames(cor_mat)) %>%
  arrange(module, desc(degree)) %>%
  pull(species)
cor_mat_guild <- cor_mat[species_guild, species_guild]

# Guild lookup
guild_lookup <- setNames(species_info$module, species_info$species)

# All pairwise correlations
all_pairs <- expand.grid(sp1 = species_guild, sp2 = species_guild, stringsAsFactors = FALSE) %>%
  filter(sp1 < sp2) %>%
  mutate(
    correlation = mapply(function(s1, s2) cor_mat_guild[s1, s2], sp1, sp2),
    guild1 = guild_lookup[sp1],
    guild2 = guild_lookup[sp2],
    type = ifelse(guild1 == guild2, "Within guild", "Between guilds")
  )

within_mean <- mean(all_pairs$correlation[all_pairs$type == "Within guild"])
between_mean <- mean(all_pairs$correlation[all_pairs$type == "Between guilds"])

# Color palettes
type_colors <- c(
  "crab" = "#E69F00", "shrimp" = "#0072B2", "fish" = "#009E73",
  "echinoderm" = "#999999", "worm" = "#56B4E9", "snail" = "#F0E442",
  "hermit" = "#CC79A7", "squat_lobster" = "#D55E00", "amphipod" = "#666666",
  "other" = "#BBBBBB"
)

guild_colors <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73", "4" = "#CC79A7")

cat("Data prepared:", n_species, "species,", n_modules, "guilds\n\n")

# ============================================================================
# OPTION 1: HEATMAP + DISTRIBUTIONS + COMPOSITION + TOP PAIRS (4-panel)
# ============================================================================

cat("Creating Option 1: 4-panel with heatmap + distributions...\n")

# Heatmap
cor_long <- as.data.frame(cor_mat_guild) %>%
  mutate(species1 = rownames(cor_mat_guild)) %>%
  tidyr::pivot_longer(-species1, names_to = "species2", values_to = "correlation") %>%
  mutate(
    species1 = factor(species1, levels = species_guild),
    species2 = factor(species2, levels = species_guild)
  )

guild_counts <- species_info %>%
  filter(species %in% species_guild) %>%
  group_by(module) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(module) %>%
  mutate(cumsum = cumsum(n))

guild_breaks <- guild_counts$cumsum
guild_centers <- c(0, guild_breaks[-n_modules]) + diff(c(0, guild_breaks)) / 2
guild_labels_vec <- paste0("Guild ", guild_counts$module, "\n(n=", guild_counts$n, ")")

p1_heatmap <- ggplot(cor_long, aes(x = species2, y = species1, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0.5, limits = c(0, 1), name = "Correlation") +
  geom_vline(xintercept = guild_breaks[-n_modules] + 0.5, color = "black", linewidth = 1) +
  geom_hline(yintercept = guild_breaks[-n_modules] + 0.5, color = "black", linewidth = 1) +
  scale_x_discrete(breaks = species_guild[round(guild_centers)], labels = guild_labels_vec, expand = c(0,0)) +
  scale_y_discrete(breaks = species_guild[round(guild_centers)], labels = guild_labels_vec, expand = c(0,0)) +
  labs(title = "A. Co-occurrence Matrix", x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_text(size = 7, face = "bold"), panel.grid = element_blank(),
        legend.position = "right", plot.title = element_text(face = "bold")) +
  coord_fixed()

p1_dist <- ggplot(all_pairs, aes(x = correlation, fill = type)) +
  geom_density(alpha = 0.7, color = "white") +
  geom_vline(xintercept = within_mean, color = "#D55E00", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = between_mean, color = "#0072B2", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("Within guild" = "#D55E00", "Between guilds" = "#0072B2"), name = NULL) +
  labs(title = "B. Within vs Between Guild", x = "Correlation", y = "Density",
       subtitle = sprintf("Within mean=%.2f, Between mean=%.2f", within_mean, between_mean)) +
  theme_minimal(base_size = 9) +
  theme(legend.position = c(0.8, 0.8), plot.title = element_text(face = "bold"))

guild_by_type <- species_info %>%
  filter(species %in% species_guild) %>%
  group_by(module, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(guild_label = paste0("Guild ", module))

p1_comp <- ggplot(guild_by_type, aes(x = guild_label, y = n, fill = type)) +
  geom_col(position = "fill", width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = type_colors, name = "Taxa") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(title = "C. Guild Composition", x = NULL, y = "Proportion") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold")) +
  guides(fill = guide_legend(nrow = 2))

top_pairs <- edge_list %>%
  arrange(desc(correlation)) %>%
  head(10) %>%
  mutate(
    guild1 = guild_lookup[sp1], guild2 = guild_lookup[sp2],
    same_guild = guild1 == guild2,
    pair = paste0(substr(sp1, 1, 10), "..-", substr(sp2, 1, 10), "..")
  )

p1_pairs <- ggplot(top_pairs, aes(x = reorder(pair, correlation), y = correlation, fill = same_guild)) +
  geom_col(alpha = 0.85, width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"),
                    labels = c("TRUE" = "Within", "FALSE" = "Between"), name = NULL) +
  labs(title = "D. Strongest Pairs", x = NULL, y = "Correlation") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

fig_opt1 <- (p1_heatmap | p1_dist) / (p1_comp | p1_pairs) +
  plot_annotation(title = "Option 1: Heatmap + Distributions",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(fig_dir, "option1_heatmap_distributions.png"), fig_opt1,
       width = 12, height = 10, dpi = 300, bg = "white")

# ============================================================================
# OPTION 2: GUILD MEMBERSHIP PANELS (who's in each guild)
# ============================================================================

cat("Creating Option 2: Guild membership panels...\n")

guild_panels <- list()
for (g in 1:4) {
  gd <- species_info %>% filter(module == g, species %in% species_guild) %>%
    mutate(sp_label = ifelse(nchar(species) > 18, paste0(substr(species, 1, 16), ".."), species))

  dominant <- names(sort(table(gd$type), decreasing = TRUE))[1]

  guild_panels[[g]] <- ggplot(gd, aes(x = reorder(sp_label, degree), y = degree, fill = type)) +
    geom_col(alpha = 0.85, width = 0.8) +
    coord_flip() +
    scale_fill_manual(values = type_colors, guide = "none") +
    labs(title = paste0("Guild ", g, " (n=", nrow(gd), ") - ", dominant),
         x = NULL, y = "Degree") +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 7))
}

fig_opt2 <- (guild_panels[[1]] | guild_panels[[2]]) / (guild_panels[[3]] | guild_panels[[4]]) +
  plot_annotation(title = "Option 2: Guild Membership (Species in Each Guild)",
                  subtitle = "Bar color = taxonomic group, ordered by network degree",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(fig_dir, "option2_guild_membership.png"), fig_opt2,
       width = 14, height = 12, dpi = 300, bg = "white")

# ============================================================================
# OPTION 3: NETWORK DIAGRAM (nodes colored by guild)
# ============================================================================

cat("Creating Option 3: Network diagram...\n")

# Rebuild graph from edge list to ensure it's a proper igraph object
g_new <- graph_from_data_frame(
  edge_list %>% dplyr::select(sp1, sp2, correlation) %>% rename(from = sp1, to = sp2, weight = correlation),
  directed = FALSE,
  vertices = species_info %>% filter(species %in% c(edge_list$sp1, edge_list$sp2)) %>%
    dplyr::select(species, type, module, degree) %>% rename(name = species)
)

set.seed(42)
layout_fr <- layout_with_fr(g_new)

node_df <- data.frame(
  x = layout_fr[, 1],
  y = layout_fr[, 2],
  species = V(g_new)$name,
  module = as.character(V(g_new)$module),
  degree = igraph::degree(g_new),
  type = centrality_df$type[match(V(g_new)$name, centrality_df$species)]
)

el <- igraph::as_edgelist(g_new, names = FALSE)
edge_df <- data.frame(
  x = layout_fr[el[,1], 1], y = layout_fr[el[,1], 2],
  xend = layout_fr[el[,2], 1], yend = layout_fr[el[,2], 2]
)

# Label top species per guild
top_species <- node_df %>%
  group_by(module) %>%
  slice_max(degree, n = 3) %>%
  ungroup()

p_network <- ggplot() +
  geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend),
               color = "gray70", alpha = 0.15, linewidth = 0.3) +
  geom_point(data = node_df, aes(x = x, y = y, fill = module, size = degree),
             shape = 21, color = "white", stroke = 0.5) +
  geom_text_repel(data = top_species, aes(x = x, y = y, label = species),
                  size = 2.5, max.overlaps = 20, segment.alpha = 0.5) +
  scale_fill_manual(values = guild_colors, name = "Guild") +
  scale_size_continuous(range = c(2, 10), name = "Degree") +
  labs(title = "Option 3: Network Diagram",
       subtitle = "Node color = guild, size = degree (connectivity)") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "right")

ggsave(file.path(fig_dir, "option3_network_diagram.png"), p_network,
       width = 12, height = 10, dpi = 300, bg = "white")

# ============================================================================
# OPTION 4: HEATMAP + GUILD MEMBERSHIP (2 panels, horizontal)
# ============================================================================

cat("Creating Option 4: Heatmap + Membership side by side...\n")

# Simplified membership - one panel with all guilds
all_members <- species_info %>%
  filter(species %in% species_guild) %>%
  mutate(
    guild = paste0("Guild ", module),
    sp_label = ifelse(nchar(species) > 15, paste0(substr(species, 1, 13), ".."), species)
  ) %>%
  arrange(module, desc(degree))

p4_members <- ggplot(all_members, aes(x = reorder(sp_label, -as.numeric(factor(guild)) * 100 + degree),
                                       y = degree, fill = type)) +
  geom_col(width = 0.8, alpha = 0.85) +
  facet_wrap(~guild, scales = "free", ncol = 2) +
  coord_flip() +
  scale_fill_manual(values = type_colors, name = "Taxa") +
  labs(title = "B. Guild Members", x = NULL, y = "Network Degree") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 6),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))

fig_opt4 <- p1_heatmap + p4_members +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(title = "Option 4: Heatmap + Guild Members",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(fig_dir, "option4_heatmap_plus_members.png"), fig_opt4,
       width = 16, height = 10, dpi = 300, bg = "white")

# ============================================================================
# OPTION 5: SIMPLE 3-PANEL (Heatmap, Composition, Top Pairs)
# ============================================================================

cat("Creating Option 5: Simple 3-panel...\n")

fig_opt5 <- p1_heatmap + p1_comp + p1_pairs +
  plot_layout(widths = c(1.2, 0.8, 0.8)) +
  plot_annotation(title = "Option 5: Simple 3-Panel",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(fig_dir, "option5_simple_3panel.png"), fig_opt5,
       width = 14, height = 6, dpi = 300, bg = "white")

# ============================================================================
# OPTION 6: NETWORK + COMPOSITION + TOP PAIRS
# ============================================================================

cat("Creating Option 6: Network + Composition + Pairs...\n")

p6_network <- ggplot() +
geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend),
               color = "gray70", alpha = 0.1, linewidth = 0.2) +
  geom_point(data = node_df, aes(x = x, y = y, fill = module, size = degree),
             shape = 21, color = "white", stroke = 0.3) +
  scale_fill_manual(values = guild_colors, name = "Guild") +
  scale_size_continuous(range = c(1.5, 8), guide = "none") +
  labs(title = "A. Co-occurrence Network") +
  theme_void() +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

fig_opt6 <- p6_network + p1_comp + p1_pairs +
  plot_layout(widths = c(1.2, 0.8, 0.8)) +
  plot_annotation(title = "Option 6: Network + Composition + Pairs",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(fig_dir, "option6_network_comp_pairs.png"), fig_opt6,
       width = 14, height = 6, dpi = 300, bg = "white")

# ============================================================================
# OPTION 7: GUILD PROFILES (radar/polar style for each guild)
# ============================================================================

cat("Creating Option 7: Guild taxonomic profiles...\n")

guild_type_matrix <- species_info %>%
  filter(species %in% species_guild) %>%
  group_by(module, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(module) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  mutate(guild = paste0("Guild ", module))

p7_profiles <- ggplot(guild_type_matrix, aes(x = type, y = pct, fill = type)) +
  geom_col(alpha = 0.85, width = 0.7) +
  facet_wrap(~guild, nrow = 1) +
  scale_fill_manual(values = type_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  coord_flip() +
  labs(title = "Option 7: Guild Taxonomic Profiles",
       subtitle = "What types of organisms dominate each guild?",
       x = NULL, y = "Proportion of Species") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 11))

ggsave(file.path(fig_dir, "option7_guild_profiles.png"), p7_profiles,
       width = 14, height = 6, dpi = 300, bg = "white")

# ============================================================================
# OPTION 8: COMPREHENSIVE 6-PANEL
# ============================================================================

cat("Creating Option 8: Comprehensive 6-panel...\n")

fig_opt8 <- (p1_heatmap | p1_dist) / (p6_network | p1_comp) / (guild_panels[[1]] | guild_panels[[2]]) +
  plot_annotation(title = "Option 8: Comprehensive View",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(fig_dir, "option8_comprehensive.png"), fig_opt8,
       width = 14, height = 16, dpi = 300, bg = "white")

# ============================================================================
# OPTION 9: ANNOTATED HEATMAP (with species names visible)
# ============================================================================

cat("Creating Option 9: Annotated heatmap with species names...\n")

# Smaller subset - top 30 species by degree
top30 <- species_info %>%
  filter(species %in% species_guild) %>%
  slice_max(degree, n = 30) %>%
  arrange(module, desc(degree)) %>%
  pull(species)

cor_mat_top30 <- cor_mat_guild[top30, top30]

cor_long_top30 <- as.data.frame(cor_mat_top30) %>%
  mutate(species1 = rownames(cor_mat_top30)) %>%
  tidyr::pivot_longer(-species1, names_to = "species2", values_to = "correlation") %>%
  mutate(
    species1 = factor(species1, levels = top30),
    species2 = factor(species2, levels = top30),
    guild1 = guild_lookup[as.character(species1)],
    guild2 = guild_lookup[as.character(species2)]
  )

# Get guild breaks for top 30
guild_counts_30 <- data.frame(species = top30) %>%
  mutate(module = guild_lookup[species]) %>%
  group_by(module) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(module) %>%
  mutate(cumsum = cumsum(n))

breaks_30 <- guild_counts_30$cumsum

p9_heatmap <- ggplot(cor_long_top30, aes(x = species2, y = species1, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                       midpoint = 0.5, limits = c(0, 1), name = "Correlation") +
  geom_vline(xintercept = breaks_30[-length(breaks_30)] + 0.5, color = "black", linewidth = 0.8) +
  geom_hline(yintercept = breaks_30[-length(breaks_30)] + 0.5, color = "black", linewidth = 0.8) +
  labs(title = "Option 9: Top 30 Species Heatmap (readable names)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 12)) +
  coord_fixed()

ggsave(file.path(fig_dir, "option9_annotated_heatmap.png"), p9_heatmap,
       width = 12, height = 10, dpi = 300, bg = "white")

# ============================================================================
# OPTION 10: MINIMAL 2-PANEL (Heatmap + Key Insight)
# ============================================================================

cat("Creating Option 10: Minimal 2-panel...\n")

p10_insight <- ggplot(all_pairs, aes(x = type, y = correlation, fill = type)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
  scale_fill_manual(values = c("Within guild" = "#D55E00", "Between guilds" = "#0072B2"), guide = "none") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  labs(title = "B. Within-Guild Species\nCo-occur More Strongly",
       subtitle = sprintf("Mean: %.2f vs %.2f (%.1fx difference)", within_mean, between_mean, within_mean/between_mean),
       x = NULL, y = "Co-occurrence Correlation") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

fig_opt10 <- p1_heatmap + p10_insight +
  plot_layout(widths = c(1.2, 0.8)) +
  plot_annotation(title = "Option 10: Minimal 2-Panel",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(file.path(fig_dir, "option10_minimal_2panel.png"), fig_opt10,
       width = 12, height = 6, dpi = 300, bg = "white")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    ALL OPTIONS GENERATED\n")
cat("============================================================\n\n")

cat("Saved to:", fig_dir, "\n\n")

cat("Options:\n")
cat("  1. option1_heatmap_distributions.png - 4-panel: heatmap + distributions + composition + pairs\n")
cat("  2. option2_guild_membership.png - Who's in each guild (species names)\n")
cat("  3. option3_network_diagram.png - Force-directed network colored by guild\n")
cat("  4. option4_heatmap_plus_members.png - Heatmap + guild members side by side\n")
cat("  5. option5_simple_3panel.png - Heatmap + composition + top pairs\n")
cat("  6. option6_network_comp_pairs.png - Network + composition + pairs\n")
cat("  7. option7_guild_profiles.png - Taxonomic profiles for each guild\n")
cat("  8. option8_comprehensive.png - 6-panel comprehensive view\n")
cat("  9. option9_annotated_heatmap.png - Top 30 species with readable names\n")
cat(" 10. option10_minimal_2panel.png - Heatmap + boxplot comparison\n")

cat("\n============================================================\n")
