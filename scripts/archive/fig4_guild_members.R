# ============================================================================
# fig4_guild_members.R - Show WHO is in each guild
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    FIGURE 4: GUILD MEMBERSHIP\n")
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

fig_dir <- file.path(PATHS$figures, "manuscript")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# EXTRACT GUILD MEMBERSHIP DATA
# ============================================================================

centrality_df <- network_results$centrality
edge_list <- network_results$edge_list

# Species info by guild
species_by_guild <- centrality_df %>%
  dplyr::select(species, type, module, degree, hub_score, occurrence) %>%
  arrange(module, desc(degree)) %>%
  mutate(
    guild = paste0("Guild ", module),
    # Clean up species names for display
    species_label = gsub("_", " ", species),
    species_label = ifelse(nchar(species_label) > 20,
                           paste0(substr(species_label, 1, 18), ".."),
                           species_label)
  )

# Guild summaries
guild_summary <- species_by_guild %>%
  group_by(guild, module) %>%
  summarise(
    n_species = n(),
    dominant_type = names(sort(table(type), decreasing = TRUE))[1],
    types = paste(unique(type), collapse = ", "),
    .groups = "drop"
  )

print(guild_summary)

# ============================================================================
# COLOR PALETTE FOR TAXA
# ============================================================================

type_colors <- c(
  "crab" = "#E69F00",
  "shrimp" = "#0072B2",
  "fish" = "#009E73",
  "echinoderm" = "#999999",
  "worm" = "#56B4E9",
  "snail" = "#F0E442",
  "hermit" = "#CC79A7",
  "squat_lobster" = "#D55E00",
  "amphipod" = "#666666",
  "other" = "#BBBBBB"
)

# ============================================================================
# PANEL: GUILD MEMBERSHIP VISUALIZATION
# ============================================================================

cat("\nCreating guild membership panels...\n")

# Create a panel for each guild showing member species
guild_panels <- list()

for (g in 1:4) {
  guild_data <- species_by_guild %>%
    filter(module == g) %>%
    arrange(desc(degree)) %>%
    mutate(rank = row_number())

  n_sp <- nrow(guild_data)
  dominant <- names(sort(table(guild_data$type), decreasing = TRUE))[1]

  # Title with description
  guild_title <- paste0("Guild ", g, " (n = ", n_sp, ")")
  guild_subtitle <- paste0("Dominated by: ", dominant)

  p <- ggplot(guild_data, aes(x = reorder(species_label, degree), y = degree, fill = type)) +
    geom_col(alpha = 0.85, width = 0.8) +
    geom_text(aes(label = type), hjust = -0.1, size = 2.5, color = "gray30") +
    coord_flip() +
    scale_fill_manual(values = type_colors, guide = "none") +
    labs(
      title = guild_title,
      subtitle = guild_subtitle,
      x = NULL,
      y = "Network Degree"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      axis.text.y = element_text(size = 8),
      panel.grid.major.y = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3)))

  guild_panels[[g]] <- p
}

# ============================================================================
# COMBINE INTO FIGURE
# ============================================================================

cat("Combining panels...\n")

# 2x2 layout of guild panels
fig4 <- (guild_panels[[1]] | guild_panels[[2]]) /
        (guild_panels[[3]] | guild_panels[[4]]) +
  plot_annotation(
    title = "Figure 4. CAFI Ecological Guilds: Species Membership",
    subtitle = "58 species cluster into 4 guilds based on co-occurrence patterns | Species ordered by network degree (connectivity)",
    caption = "Guilds identified via Louvain community detection on volume-corrected co-occurrence network (r > 0.3, FDR < 0.05).\nNetwork transitivity = 0.94 (z = 36.1 vs null, p < 0.001). Bar color = taxonomic group.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      plot.caption = element_text(size = 8, hjust = 0, color = "gray50")
    )
  )

ggsave(
  file.path(fig_dir, "fig4_guild_members.png"),
  fig4, width = 14, height = 14, dpi = 300, bg = "white"
)

cat("\n[OK] Saved: fig4_guild_members.png\n")

# ============================================================================
# ALSO CREATE A SUMMARY TABLE
# ============================================================================

guild_table <- species_by_guild %>%
  group_by(guild) %>%
  summarise(
    n_species = n(),
    species_list = paste(species, collapse = ", "),
    taxa = paste(names(sort(table(type), decreasing = TRUE)), collapse = " > "),
    hub_species = species[which.max(degree)],
    .groups = "drop"
  )

cat("\n")
cat("============================================================\n")
cat("GUILD SUMMARY\n")
cat("============================================================\n")

for (g in 1:4) {
  guild_data <- species_by_guild %>% filter(module == g)
  cat("\n--- Guild", g, "(n =", nrow(guild_data), "species) ---\n")
  cat("Dominant taxa:", names(sort(table(guild_data$type), decreasing = TRUE))[1], "\n")
  cat("Hub species:", guild_data$species[1], "(degree =", guild_data$degree[1], ")\n")
  cat("Members:\n")
  for (i in 1:nrow(guild_data)) {
    cat("  ", guild_data$species[i], "(", guild_data$type[i], ")\n")
  }
}

cat("\n============================================================\n")
cat("    COMPLETE\n")
cat("============================================================\n")
