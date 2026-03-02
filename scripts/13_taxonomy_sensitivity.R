# ============================================================================
# 13_taxonomy_sensitivity.R - Taxonomic Resolution Sensitivity Analysis
# ============================================================================
#
# PURPOSE: Test how species inclusion/exclusion criteria and taxonomic
#          resolution affect the main results (Q1-Q4)
#
# SCENARIOS (defined in 01_load_data.R, pre-built as taxonomy_scenario_data):
#   1. Baseline    — All OTUs as-is (current pipeline)
#   2. Species-only — Exclude genus/family/higher OTUs (keep species-level only)
#   3. Merge-up    — Merge genus-only records into the most common species
#                     within that genus; keep genus bin if no species exist
#   4. Lump-down   — Collapse all species within a genus to genus-level
#   5. Rare-excluded — Remove OTUs present on fewer than 3 corals
#
# ANALYSES RE-RUN PER SCENARIO:
#   - Total CAFI abundance scaling (NB GLM beta)
#   - Species richness scaling (Poisson GLM z)
#   - Shannon diversity scaling
#   - Top-10 species scaling (where applicable)
#   - PERMANOVA R-squared for site effect
#   - Network modularity Q and number of modules
#   - Rarefied richness -> condition relationship
#
# DEPENDS ON: 01_load_data.R (otu_taxonomy, taxonomy_scenario_data,
#             coral_master, condition_scores)
#
# OUTPUTS:
#   - output/tables/taxonomy_sensitivity.csv
#   - output/tables/taxonomy_sensitivity_species_scaling.csv
#   - output/figures/supplement/figS8_taxonomy_sensitivity.png
#
# Author: CAFI Survey Analysis Pipeline
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("    TAXONOMIC RESOLUTION SENSITIVITY ANALYSIS\n")
cat("============================================================\n\n")

# ============================================================================
# SETUP — Load pre-built scenario data from 01_load_data.R
# ============================================================================

if (!exists("PATHS")) source(here::here("scripts/00_setup.R"))

suppressPackageStartupMessages(library(igraph))

# Load pre-computed objects (scenario data built in 01_load_data.R section 6a)
coral_master     <- load_object("coral_master")
cafi_clean       <- load_object("cafi_clean")
condition_scores <- load_object("condition_scores")
otu_taxonomy     <- load_object("otu_taxonomy")
scenario_data    <- load_object("taxonomy_scenario_data")

fig_dir <- file.path(PATHS$figures, "supplement")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Scenario ordering (for figures)
scenario_names <- names(scenario_data)

cat("[OK] Loaded pre-built scenario data from 01_load_data.R\n")
cat("     Scenarios:", paste(scenario_names, collapse = ", "), "\n")
cat("     Corals:", nrow(coral_master), "\n")
cat("     CAFI records:", nrow(cafi_clean), "\n")
cat("     OTU taxonomy:", nrow(otu_taxonomy), "OTUs\n\n")

# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================
# These re-fit simplified versions of models from scripts 05, 02, 06, 09.
# They are intentionally self-contained here (not shared with the full scripts)
# because the sensitivity versions use stripped-down specifications.
# ============================================================================

# --- A. Total CAFI abundance scaling (NB GLM) ---
fit_abundance_scaling <- function(metrics_df) {
  tryCatch({
    model <- MASS::glm.nb(total_cafi ~ log(volume) + site, data = metrics_df)
    coefs <- summary(model)$coefficients
    beta <- coefs["log(volume)", "Estimate"]
    se <- coefs["log(volume)", "Std. Error"]
    p <- coefs["log(volume)", "Pr(>|z|)"]
    ci <- tryCatch(confint(model, "log(volume)"),
                   error = function(e) beta + c(-1.96, 1.96) * se)
    list(beta = beta, se = se, ci_lower = ci[1], ci_upper = ci[2], p = p)
  }, error = function(e) {
    list(beta = NA, se = NA, ci_lower = NA, ci_upper = NA, p = NA)
  })
}

# --- B. Species richness scaling (Poisson GLM) ---
fit_richness_scaling <- function(metrics_df) {
  tryCatch({
    model <- glm(otu_richness ~ log(volume) + site, family = poisson, data = metrics_df)

    # Check overdispersion
    overdispersion <- sum(residuals(model, type = "pearson")^2) / model$df.residual
    if (overdispersion > 2) {
      model <- glm(otu_richness ~ log(volume) + site, family = quasipoisson, data = metrics_df)
    }

    coefs <- summary(model)$coefficients
    z <- coefs["log(volume)", "Estimate"]
    se <- coefs["log(volume)", "Std. Error"]
    p_col <- intersect(colnames(coefs), c("Pr(>|z|)", "Pr(>|t|)"))
    p <- if (length(p_col) > 0) coefs["log(volume)", p_col[1]] else NA_real_
    ci <- z + c(-1.96, 1.96) * se
    list(z = z, se = se, ci_lower = ci[1], ci_upper = ci[2], p = p,
         overdispersion = overdispersion)
  }, error = function(e) {
    cat("    Richness scaling error:", e$message, "\n")
    list(z = NA, se = NA, ci_lower = NA, ci_upper = NA, p = NA,
         overdispersion = NA)
  })
}

# --- C. Shannon diversity scaling (Gaussian GLM) ---
fit_shannon_scaling <- function(metrics_df) {
  tryCatch({
    model <- lm(shannon ~ log(volume) + site, data = metrics_df)
    coefs <- summary(model)$coefficients
    beta <- coefs["log(volume)", "Estimate"]
    se <- coefs["log(volume)", "Std. Error"]
    p <- coefs["log(volume)", "Pr(>|t|)"]
    ci <- confint(model, "log(volume)")
    r2 <- summary(model)$adj.r.squared
    list(beta = beta, se = se, ci_lower = ci[1], ci_upper = ci[2], p = p, r2 = r2)
  }, error = function(e) {
    list(beta = NA, se = NA, ci_lower = NA, ci_upper = NA, p = NA, r2 = NA)
  })
}

# --- D. PERMANOVA for site effect ---
fit_permanova <- function(comm_mat, coral_master_df) {
  tryCatch({
    common_ids <- intersect(rownames(comm_mat), coral_master_df$coral_id)
    comm_sub <- comm_mat[common_ids, , drop = FALSE]
    coral_sub <- coral_master_df %>% filter(coral_id %in% common_ids)

    nonzero <- rowSums(comm_sub) > 0
    comm_sub <- comm_sub[nonzero, , drop = FALSE]
    coral_sub <- coral_sub %>% filter(coral_id %in% rownames(comm_sub))
    comm_sub <- comm_sub[, colSums(comm_sub) > 0, drop = FALSE]

    if (nrow(comm_sub) < 10 || ncol(comm_sub) < 3) {
      return(list(r2_site = NA, f_site = NA, p_site = NA,
                  r2_volume = NA, f_volume = NA, p_volume = NA))
    }

    comm_hell <- vegan::decostand(comm_sub, method = "hellinger")

    set.seed(42)
    perm_result <- vegan::adonis2(
      comm_hell ~ log(volume) + site,
      data = coral_sub,
      permutations = 999,
      method = "bray",
      by = "margin"
    )
    perm_df <- as.data.frame(perm_result)

    list(
      r2_site = perm_df["site", "R2"],
      f_site = perm_df["site", "F"],
      p_site = perm_df["site", "Pr(>F)"],
      r2_volume = perm_df["log(volume)", "R2"],
      f_volume = perm_df["log(volume)", "F"],
      p_volume = perm_df["log(volume)", "Pr(>F)"]
    )
  }, error = function(e) {
    cat("    PERMANOVA error:", e$message, "\n")
    list(r2_site = NA, f_site = NA, p_site = NA,
         r2_volume = NA, f_volume = NA, p_volume = NA)
  })
}

# --- E. Network modularity ---
fit_network <- function(comm_mat, coral_master_df, cafi_df) {
  tryCatch({
    common_ids <- intersect(rownames(comm_mat), coral_master_df$coral_id)
    comm_sub <- comm_mat[common_ids, , drop = FALSE]

    comm_binary <- comm_sub
    comm_binary[comm_binary > 0] <- 1

    min_occ <- ceiling(nrow(comm_binary) * 0.05)
    species_keep <- colSums(comm_binary) >= min_occ
    comm_filt <- comm_binary[, species_keep, drop = FALSE]

    n_species <- ncol(comm_filt)
    if (n_species < 5) {
      return(list(modularity_Q = NA, n_modules = NA, n_species_network = NA,
                  n_edges = NA, mean_degree = NA))
    }

    volume_vec <- coral_master_df$volume[match(rownames(comm_filt), coral_master_df$coral_id)]
    log_vol <- log(volume_vec)

    resid_mat <- matrix(NA, nrow = nrow(comm_filt), ncol = ncol(comm_filt))
    colnames(resid_mat) <- colnames(comm_filt)
    rownames(resid_mat) <- rownames(comm_filt)

    for (sp in seq_len(ncol(comm_filt))) {
      y <- comm_filt[, sp]
      fit <- tryCatch(
        glm(y ~ log_vol, family = binomial(link = "logit")),
        warning = function(w) glm(y ~ log_vol, family = binomial(link = "logit")),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        resid_mat[, sp] <- residuals(fit, type = "deviance")
      } else {
        resid_mat[, sp] <- y - mean(y)
      }
    }

    cor_mat <- cor(resid_mat, method = "spearman", use = "pairwise.complete.obs")
    cor_mat[is.na(cor_mat)] <- 0

    threshold <- 0.3
    n_sp <- ncol(resid_mat)
    p_mat <- matrix(1, nrow = n_sp, ncol = n_sp)

    for (i in 1:(n_sp - 1)) {
      for (j in (i + 1):n_sp) {
        if (abs(cor_mat[i, j]) > 0.2) {
          ct <- tryCatch(
            cor.test(resid_mat[, i], resid_mat[, j], method = "spearman"),
            error = function(e) list(p.value = 1)
          )
          p_mat[i, j] <- ct$p.value
          p_mat[j, i] <- ct$p.value
        }
      }
    }

    edge_idx_raw <- which(upper.tri(cor_mat) & cor_mat > threshold, arr.ind = TRUE)
    if (nrow(edge_idx_raw) == 0) {
      return(list(modularity_Q = NA, n_modules = NA, n_species_network = 0,
                  n_edges = 0, mean_degree = NA))
    }

    raw_pvals <- sapply(seq_len(nrow(edge_idx_raw)), function(k)
      p_mat[edge_idx_raw[k, 1], edge_idx_raw[k, 2]])
    fdr_pvals <- p.adjust(raw_pvals, method = "BH")
    sig_mask <- fdr_pvals < 0.05
    edge_idx <- edge_idx_raw[sig_mask, , drop = FALSE]

    if (nrow(edge_idx) == 0) {
      return(list(modularity_Q = NA, n_modules = NA, n_species_network = 0,
                  n_edges = 0, mean_degree = NA))
    }

    edge_list <- data.frame(
      sp1 = colnames(cor_mat)[edge_idx[, 1]],
      sp2 = colnames(cor_mat)[edge_idx[, 2]],
      weight = cor_mat[edge_idx],
      stringsAsFactors = FALSE
    )

    g <- igraph::graph_from_data_frame(edge_list[, c("sp1", "sp2")], directed = FALSE)
    igraph::E(g)$weight <- edge_list$weight

    set.seed(42)
    comm_det <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
    mod_Q <- igraph::modularity(comm_det)
    n_mod <- length(unique(igraph::membership(comm_det)))

    list(
      modularity_Q = mod_Q,
      n_modules = n_mod,
      n_species_network = igraph::vcount(g),
      n_edges = igraph::ecount(g),
      mean_degree = mean(igraph::degree(g))
    )
  }, error = function(e) {
    cat("    Network error:", e$message, "\n")
    list(modularity_Q = NA, n_modules = NA, n_species_network = NA,
         n_edges = NA, mean_degree = NA)
  })
}

# --- F. Rarefied richness -> condition ---
fit_rarefied_condition <- function(comm_mat, condition_df) {
  tryCatch({
    common_ids <- intersect(rownames(comm_mat), condition_df$coral_id)
    comm_sub <- comm_mat[common_ids, , drop = FALSE]
    cond_sub <- condition_df %>% filter(coral_id %in% common_ids)
    comm_sub <- comm_sub[, colSums(comm_sub) > 0, drop = FALSE]

    abund <- rowSums(comm_sub)
    rarefy_depth <- 20
    has_enough <- abund >= rarefy_depth

    if (sum(has_enough) < 20) {
      return(list(beta_rare = NA, se_rare = NA, p_rare = NA,
                  beta_raw = NA, p_raw = NA, n = sum(has_enough)))
    }

    comm_eligible <- comm_sub[has_enough, , drop = FALSE]
    rare_rich <- vegan::rarefy(comm_eligible, sample = rarefy_depth)

    analysis_df <- cond_sub %>%
      filter(coral_id %in% names(rare_rich)) %>%
      mutate(
        rarefied_richness = rare_rich[coral_id],
        otu_richness = rowSums(comm_sub[coral_id, , drop = FALSE] > 0),
        log_volume = log(volume)
      ) %>%
      filter(!is.na(condition_score), !is.na(volume))

    if (nrow(analysis_df) < 15) {
      return(list(beta_rare = NA, se_rare = NA, p_rare = NA,
                  beta_raw = NA, p_raw = NA, n = nrow(analysis_df)))
    }

    m_raw <- lm(condition_score ~ otu_richness + log_volume + site, data = analysis_df)
    s_raw <- summary(m_raw)
    p_raw <- s_raw$coefficients["otu_richness", "Pr(>|t|)"]
    b_raw <- coef(m_raw)["otu_richness"]

    m_rare <- lm(condition_score ~ rarefied_richness + log_volume + site, data = analysis_df)
    s_rare <- summary(m_rare)
    p_rare <- s_rare$coefficients["rarefied_richness", "Pr(>|t|)"]
    b_rare <- coef(m_rare)["rarefied_richness"]
    se_rare <- s_rare$coefficients["rarefied_richness", "Std. Error"]

    list(beta_rare = b_rare, se_rare = se_rare, p_rare = p_rare,
         beta_raw = b_raw, p_raw = p_raw, n = nrow(analysis_df))
  }, error = function(e) {
    cat("    Condition model error:", e$message, "\n")
    list(beta_rare = NA, se_rare = NA, p_rare = NA,
         beta_raw = NA, p_raw = NA, n = NA)
  })
}

# --- G. Top-10 species scaling ---
fit_species_scaling <- function(comm_mat, coral_master_df) {
  sp_abund <- colSums(comm_mat)
  sp_abund <- sort(sp_abund, decreasing = TRUE)
  top_sp <- names(sp_abund)[1:min(10, length(sp_abund))]

  results <- list()
  for (sp in top_sp) {
    sp_data <- coral_master_df %>%
      filter(coral_id %in% rownames(comm_mat)) %>%
      mutate(abundance = comm_mat[coral_id, sp])

    if (sum(sp_data$abundance > 0) < 15) {
      results[[sp]] <- tibble(
        species = sp, beta = NA, se = NA, ci_lower = NA, ci_upper = NA,
        p = NA, total_n = sum(sp_data$abundance), n_occupied = sum(sp_data$abundance > 0)
      )
      next
    }

    tryCatch({
      model <- MASS::glm.nb(abundance ~ log(volume) + site, data = sp_data)
      coefs <- summary(model)$coefficients
      beta <- coefs["log(volume)", "Estimate"]
      se <- coefs["log(volume)", "Std. Error"]
      p <- coefs["log(volume)", "Pr(>|z|)"]
      ci <- tryCatch(confint(model, "log(volume)"),
                     error = function(e) beta + c(-1.96, 1.96) * se)

      results[[sp]] <- tibble(
        species = sp, beta = beta, se = se, ci_lower = ci[1], ci_upper = ci[2],
        p = p, total_n = sum(sp_data$abundance), n_occupied = sum(sp_data$abundance > 0)
      )
    }, error = function(e) {
      results[[sp]] <<- tibble(
        species = sp, beta = NA, se = NA, ci_lower = NA, ci_upper = NA,
        p = NA, total_n = sum(sp_data$abundance), n_occupied = sum(sp_data$abundance > 0)
      )
    })
  }

  bind_rows(results)
}

# ============================================================================
# PREPARE CONDITION DATA
# ============================================================================

condition_data <- condition_scores %>%
  dplyr::select(coral_id, site, condition_score) %>%
  left_join(coral_master %>% dplyr::select(coral_id, volume), by = "coral_id") %>%
  filter(!is.na(condition_score), !is.na(volume)) %>%
  mutate(site = factor(site))

cat("  Condition data: n =", nrow(condition_data), "corals\n\n")

# ============================================================================
# RUN ALL SCENARIOS (using pre-built data from 01_load_data.R)
# ============================================================================

cat("============================================================\n")
cat("RUNNING SENSITIVITY SCENARIOS\n")
cat("============================================================\n\n")

all_results <- list()
all_species_scaling <- list()

for (scenario_name in scenario_names) {
  sd <- scenario_data[[scenario_name]]

  cat("------------------------------------------------------------\n")
  cat("SCENARIO:", scenario_name,
      "| Records:", sd$n_records,
      "| OTUs:", sd$n_otus,
      "| Matrix:", nrow(sd$community_matrix), "x", ncol(sd$community_matrix), "\n")
  cat("------------------------------------------------------------\n\n")

  # Use pre-built objects
  comm_mat     <- sd$community_matrix
  metrics      <- sd$metrics
  cafi_modified <- sd$cafi_modified

  cat("  Running analyses...\n")

  # A. Abundance scaling
  cat("    [A] Abundance scaling (NB GLM)... ")
  abund_result <- fit_abundance_scaling(metrics)
  cat("beta =", round(abund_result$beta, 3), "\n")

  # B. Richness scaling
  cat("    [B] Richness scaling (Poisson GLM)... ")
  rich_result <- fit_richness_scaling(metrics)
  cat("z =", round(rich_result$z, 3), "\n")

  # C. Shannon scaling
  cat("    [C] Shannon scaling... ")
  shan_result <- fit_shannon_scaling(metrics)
  cat("beta =", round(shan_result$beta, 3), "\n")

  # D. PERMANOVA
  cat("    [D] PERMANOVA (site effect)... ")
  perm_result <- fit_permanova(comm_mat, coral_master)
  cat("R2_site =", round(perm_result$r2_site, 3), "\n")

  # E. Network modularity
  cat("    [E] Network modularity... ")
  net_result <- fit_network(comm_mat, coral_master, cafi_modified)
  cat("Q =", round(net_result$modularity_Q, 3),
      ", modules =", net_result$n_modules, "\n")

  # F. Rarefied richness -> condition
  cat("    [F] Rarefied richness -> condition... ")
  cond_result <- fit_rarefied_condition(comm_mat, condition_data)
  cat("p_rare =", round(cond_result$p_rare, 3), "\n")

  # G. Species-level scaling (top 10)
  cat("    [G] Top-10 species scaling...\n")
  sp_scaling <- fit_species_scaling(comm_mat, coral_master)
  sp_scaling$scenario <- scenario_name
  all_species_scaling[[scenario_name]] <- sp_scaling

  # Compile results
  all_results[[scenario_name]] <- tibble(
    scenario = scenario_name,
    n_records = sd$n_records,
    n_otus = sd$n_otus,
    n_corals_with_cafi = n_distinct(cafi_modified$coral_id),
    n_comm_species = ncol(comm_mat),

    abundance_beta = abund_result$beta,
    abundance_se = abund_result$se,
    abundance_ci_lower = abund_result$ci_lower,
    abundance_ci_upper = abund_result$ci_upper,
    abundance_p = abund_result$p,

    richness_z = rich_result$z,
    richness_se = rich_result$se,
    richness_ci_lower = rich_result$ci_lower,
    richness_ci_upper = rich_result$ci_upper,
    richness_p = rich_result$p,

    shannon_beta = shan_result$beta,
    shannon_se = shan_result$se,
    shannon_ci_lower = shan_result$ci_lower,
    shannon_ci_upper = shan_result$ci_upper,
    shannon_p = shan_result$p,
    shannon_r2 = shan_result$r2,

    permanova_r2_site = perm_result$r2_site,
    permanova_f_site = perm_result$f_site,
    permanova_p_site = perm_result$p_site,
    permanova_r2_volume = perm_result$r2_volume,
    permanova_f_volume = perm_result$f_volume,
    permanova_p_volume = perm_result$p_volume,

    network_modularity_Q = net_result$modularity_Q,
    network_n_modules = net_result$n_modules,
    network_n_species = net_result$n_species_network,
    network_n_edges = net_result$n_edges,
    network_mean_degree = net_result$mean_degree,

    condition_beta_rarefied = cond_result$beta_rare,
    condition_se_rarefied = cond_result$se_rare,
    condition_p_rarefied = cond_result$p_rare,
    condition_beta_raw = cond_result$beta_raw,
    condition_p_raw = cond_result$p_raw,
    condition_n = cond_result$n
  )

  cat("\n")
}

# ============================================================================
# COMPILE AND SAVE RESULTS
# ============================================================================

cat("============================================================\n")
cat("COMPILING RESULTS\n")
cat("============================================================\n\n")

results_df <- bind_rows(all_results)
species_scaling_df <- bind_rows(all_species_scaling)

save_table(results_df, "taxonomy_sensitivity")
save_table(species_scaling_df, "taxonomy_sensitivity_species_scaling")

# --- Formatted results table ---
cat("=== TAXONOMY SENSITIVITY RESULTS ===\n\n")

cat(sprintf("%-15s %5s %5s | %6s %6s | %6s %6s | %5s %6s | %5s %4s | %6s %6s\n",
            "Scenario", "Recs", "OTUs",
            "AbundB", "   SE",
            "RichZ", "   SE",
            "R2_st", "p_site",
            "ModQ", "nMod",
            "RarB", "p_rare"))
cat(paste(rep("-", 110), collapse = ""), "\n")

for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  cat(sprintf("%-15s %5d %5d | %6.3f %6.3f | %6.3f %6.3f | %5.3f %6.3f | %5.3f %4d | %6.4f %6.3f\n",
              r$scenario, r$n_records, r$n_otus,
              r$abundance_beta, r$abundance_se,
              r$richness_z, r$richness_se,
              r$permanova_r2_site, r$permanova_p_site,
              r$network_modularity_Q, r$network_n_modules,
              r$condition_beta_rarefied, r$condition_p_rarefied))
}

cat("\n")

# --- Qualitative stability assessment ---
cat("=== QUALITATIVE STABILITY ASSESSMENT ===\n\n")

baseline <- results_df %>% filter(scenario == "Baseline")

for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  if (r$scenario == "Baseline") next

  cat("  ", r$scenario, "vs Baseline:\n")

  abund_sublinear <- !is.na(r$abundance_ci_upper) && r$abundance_ci_upper < 1
  cat("    Q1 Abundance scaling: beta =", round(r$abundance_beta, 3),
      "[", round(r$abundance_ci_lower, 3), ",", round(r$abundance_ci_upper, 3), "]",
      ifelse(abund_sublinear, "-> SUBLINEAR (consistent)", "-> CI overlaps 1 (CHANGED)"), "\n")

  rich_sublinear <- !is.na(r$richness_ci_upper) && r$richness_ci_upper < 1
  cat("    Q1 Richness scaling:  z =", round(r$richness_z, 3),
      "[", round(r$richness_ci_lower, 3), ",", round(r$richness_ci_upper, 3), "]",
      ifelse(rich_sublinear, "-> SUBLINEAR (consistent)", "-> CI overlaps 1 (CHANGED)"), "\n")

  site_sig <- !is.na(r$permanova_p_site) && r$permanova_p_site < 0.05
  cat("    Q2 Site effect:       R2 =", round(r$permanova_r2_site, 3),
      ", p =", round(r$permanova_p_site, 3),
      ifelse(site_sig, "-> SIGNIFICANT (consistent)", "-> NS (CHANGED)"), "\n")

  cond_ns <- !is.na(r$condition_p_rarefied) && r$condition_p_rarefied > 0.05
  cat("    Q3 Rarefied -> cond:  p =", round(r$condition_p_rarefied, 3),
      ifelse(cond_ns, "-> NS (consistent)", "-> SIGNIFICANT (CHANGED)"), "\n")

  cat("\n")
}

# ============================================================================
# FIGURE: Forest-Plot Style Comparison (Fig S8)
# ============================================================================

cat("============================================================\n")
cat("GENERATING SENSITIVITY FIGURE\n")
cat("============================================================\n\n")

# Panel data
panel_a_data <- results_df %>%
  dplyr::select(scenario, beta = abundance_beta, ci_lower = abundance_ci_lower,
         ci_upper = abundance_ci_upper) %>%
  mutate(scenario = factor(scenario, levels = rev(scenario_names)))

panel_b_data <- results_df %>%
  dplyr::select(scenario, beta = richness_z, ci_lower = richness_ci_lower,
         ci_upper = richness_ci_upper) %>%
  mutate(scenario = factor(scenario, levels = rev(scenario_names)))

panel_c_data <- results_df %>%
  dplyr::select(scenario, beta = permanova_r2_site) %>%
  mutate(ci_lower = NA_real_, ci_upper = NA_real_,
         scenario = factor(scenario, levels = rev(scenario_names)))

panel_d_data <- results_df %>%
  dplyr::select(scenario, beta = network_modularity_Q) %>%
  mutate(ci_lower = NA_real_, ci_upper = NA_real_,
         scenario = factor(scenario, levels = rev(scenario_names)))

scenario_colors <- c(
  "Baseline" = "gray30",
  "Species-only" = "#E69F00",
  "Merge-up" = "#56B4E9",
  "Lump-down" = "#009E73",
  "Rare-excluded" = "#CC79A7"
)

p_a <- ggplot(panel_a_data, aes(x = beta, y = scenario, color = scenario)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, linewidth = 0.6,
                orientation = "y") +
  scale_color_manual(values = scenario_colors, guide = "none") +
  scale_x_continuous(breaks = seq(0.3, 1.1, 0.1)) +
  coord_cartesian(xlim = c(0.35, 1.1)) +
  labs(x = expression(beta~"(dashed = Field of Dreams)"), y = NULL,
       title = "A. Abundance scaling") +
  theme_publication(base_size = 10) +
  theme(plot.margin = margin(5, 10, 5, 5, "mm"))

p_b <- ggplot(panel_b_data, aes(x = beta, y = scenario, color = scenario)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, linewidth = 0.6,
                orientation = "y") +
  scale_color_manual(values = scenario_colors, guide = "none") +
  scale_x_continuous(breaks = seq(0.2, 1.1, 0.1)) +
  coord_cartesian(xlim = c(0.2, 1.1)) +
  labs(x = expression(italic(z)~"(dashed = proportional)"), y = NULL,
       title = "B. Richness scaling") +
  theme_publication(base_size = 10) +
  theme(plot.margin = margin(5, 10, 5, 5, "mm"))

r2_range <- range(panel_c_data$beta, na.rm = TRUE)
r2_pad <- max(diff(r2_range) * 0.5, 0.01)

p_c <- ggplot(panel_c_data %>% filter(!is.na(beta)),
              aes(x = beta, y = scenario, color = scenario)) +
  geom_point(size = 3) +
  scale_color_manual(values = scenario_colors, guide = "none") +
  coord_cartesian(xlim = c(max(0, r2_range[1] - r2_pad), r2_range[2] + r2_pad)) +
  labs(x = expression(R^2), y = NULL,
       title = expression("C. PERMANOVA"~R^2~"(site)")) +
  theme_publication(base_size = 10) +
  theme(plot.margin = margin(5, 10, 5, 5, "mm"))

q_range <- range(panel_d_data$beta, na.rm = TRUE)
q_pad <- max(diff(q_range) * 0.5, 0.005)

p_d <- ggplot(panel_d_data %>% filter(!is.na(beta)),
              aes(x = beta, y = scenario, color = scenario)) +
  geom_point(size = 3) +
  scale_color_manual(values = scenario_colors, guide = "none") +
  coord_cartesian(xlim = c(max(0, q_range[1] - q_pad), q_range[2] + q_pad)) +
  labs(x = "Modularity Q", y = NULL,
       title = "D. Network modularity") +
  theme_publication(base_size = 10) +
  theme(plot.margin = margin(5, 10, 5, 5, "mm"))

p_combined <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = "Sensitivity of Key Results to Taxonomic Resolution",
    subtitle = paste0("5 scenarios tested across ", nrow(coral_master), " corals, ",
                      nrow(cafi_clean), " CAFI records"),
    theme = theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  )

ggsave(file.path(fig_dir, "figS8_taxonomy_sensitivity.png"),
       p_combined, width = 10, height = 7, dpi = 300, bg = "white")
cat("  Saved: figS8_taxonomy_sensitivity.png\n\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("============================================================\n")
cat("TAXONOMY SENSITIVITY ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

n_scenarios <- nrow(results_df) - 1
n_stable_abundance <- sum(results_df$abundance_ci_upper[results_df$scenario != "Baseline"] < 1, na.rm = TRUE)
n_stable_richness <- sum(results_df$richness_ci_upper[results_df$scenario != "Baseline"] < 1, na.rm = TRUE)
n_stable_site <- sum(results_df$permanova_p_site[results_df$scenario != "Baseline"] < 0.05, na.rm = TRUE)
n_stable_condition <- sum(results_df$condition_p_rarefied[results_df$scenario != "Baseline"] > 0.05, na.rm = TRUE)

cat("STABILITY SUMMARY:\n")
cat("  Q1 Abundance sublinear:  ", n_stable_abundance, "/", n_scenarios, "scenarios agree\n")
cat("  Q1 Richness sublinear:   ", n_stable_richness, "/", n_scenarios, "scenarios agree\n")
cat("  Q2 Site effect signif:   ", n_stable_site, "/", n_scenarios, "scenarios agree\n")
cat("  Q3 Rarefied richness NS: ", n_stable_condition, "/", n_scenarios, "scenarios agree\n")

cat("\nOutputs:\n")
cat("  Table: output/tables/taxonomy_sensitivity.csv\n")
cat("  Table: output/tables/taxonomy_sensitivity_species_scaling.csv\n")
cat("  Figure: output/figures/supplement/figS8_taxonomy_sensitivity.png\n\n")

cat("Script 13 complete.\n\n")
