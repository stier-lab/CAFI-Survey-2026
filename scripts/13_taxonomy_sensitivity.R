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
#   - output/tables/network_topology_sensitivity.csv
#   - output/tables/network_edge_overlap.csv
#   - output/tables/network_hub_stability.csv
#   - output/figures/supplement/figS6_taxonomy_sensitivity.png
#
# Author: CAFI Survey Analysis Pipeline
# Last Updated: 2026-03-04
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

# Resolve igraph namespace conflicts
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflict_prefer("union", "base")
  conflicted::conflict_prefer("intersect", "base")
  conflicted::conflict_prefer("setdiff", "base")
  conflicted::conflict_prefer("groups", "dplyr")  # needed: igraph also exports groups()
  conflicted::conflict_prefer("as_data_frame", "igraph")
}

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
                   error = function(e) {
                     cat("    NOTE: confint() failed for abundance scaling, using Wald CI\n")
                     beta + c(-1.96, 1.96) * se
                   })
    list(beta = beta, se = se, ci_lower = ci[1], ci_upper = ci[2], p = p)
  }, error = function(e) {
    cat("    Abundance scaling error:", e$message, "\n")
    list(beta = NA, se = NA, ci_lower = NA, ci_upper = NA, p = NA)
  })
}

# --- B. Species richness scaling (Poisson GLM) ---
fit_richness_scaling <- function(metrics_df) {
  tryCatch({
    model <- glm(otu_richness ~ log(volume) + site, family = poisson, data = metrics_df)

    # Check overdispersion
    # Overdispersion threshold matches main analysis (script 05)
    overdispersion <- sum(residuals(model, type = "pearson")^2) / model$df.residual
    if (overdispersion > 1.5) {
      model <- glm(otu_richness ~ log(volume) + site, family = quasipoisson, data = metrics_df)
    }

    coefs <- summary(model)$coefficients
    beta <- coefs["log(volume)", "Estimate"]
    se <- coefs["log(volume)", "Std. Error"]
    p_col <- intersect(colnames(coefs), c("Pr(>|z|)", "Pr(>|t|)"))
    p <- if (length(p_col) > 0) coefs["log(volume)", p_col[1]] else NA_real_
    ci <- beta + c(-1.96, 1.96) * se
    list(beta = beta, se = se, ci_lower = ci[1], ci_upper = ci[2], p = p,
         overdispersion = overdispersion)
  }, error = function(e) {
    cat("    Richness scaling error:", e$message, "\n")
    list(beta = NA, se = NA, ci_lower = NA, ci_upper = NA, p = NA,
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
    cat("    Shannon scaling error:", e$message, "\n")
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
fit_network <- function(comm_mat, coral_master_df, return_full = FALSE) {
  empty_full <- list(modularity_Q = NA, n_modules = NA, n_species_network = NA,
                     n_edges = NA, mean_degree = NA, density = NA, transitivity = NA,
                     graph = NULL, edge_list = NULL, centrality = NULL,
                     membership = NULL, species_in_network = character(0))
  empty_summary <- list(modularity_Q = NA, n_modules = NA, n_species_network = NA,
                        n_edges = NA, mean_degree = NA)
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
      return(if (return_full) empty_full else empty_summary)
    }

    volume_vec <- coral_master_df$volume[match(rownames(comm_filt), coral_master_df$coral_id)]
    log_vol <- log(volume_vec)
    # Include site as covariate (consistent with main co-occurrence analysis in script 06)
    site_fac <- factor(coral_master_df$site[match(rownames(comm_filt), coral_master_df$coral_id)])

    resid_mat <- matrix(NA, nrow = nrow(comm_filt), ncol = ncol(comm_filt))
    colnames(resid_mat) <- colnames(comm_filt)
    rownames(resid_mat) <- rownames(comm_filt)

    for (sp in seq_len(ncol(comm_filt))) {
      y <- comm_filt[, sp]
      # Try volume + site first; fall back to volume-only if site causes separation
      fit <- tryCatch(
        glm(y ~ log_vol + site_fac, family = binomial(link = "logit")),
        warning = function(w) suppressWarnings(glm(y ~ log_vol + site_fac, family = binomial(link = "logit"))),
        error = function(e) NULL
      )
      if (!is.null(fit) && any(!is.finite(coef(fit)))) fit <- NULL
      if (is.null(fit)) {
        fit <- tryCatch(
          glm(y ~ log_vol, family = binomial(link = "logit")),
          warning = function(w) glm(y ~ log_vol, family = binomial(link = "logit")),
          error = function(e) NULL
        )
      }
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
      result <- list(modularity_Q = NA, n_modules = NA, n_species_network = 0,
                     n_edges = 0, mean_degree = NA)
      if (return_full) {
        result$density <- NA; result$transitivity <- NA
        result$graph <- NULL; result$edge_list <- NULL; result$centrality <- NULL
        result$membership <- NULL; result$species_in_network <- character(0)
      }
      return(result)
    }

    raw_pvals <- sapply(seq_len(nrow(edge_idx_raw)), function(k)
      p_mat[edge_idx_raw[k, 1], edge_idx_raw[k, 2]])
    fdr_pvals <- p.adjust(raw_pvals, method = "BH")
    sig_mask <- fdr_pvals < 0.05
    edge_idx <- edge_idx_raw[sig_mask, , drop = FALSE]

    if (nrow(edge_idx) == 0) {
      result <- list(modularity_Q = NA, n_modules = NA, n_species_network = 0,
                     n_edges = 0, mean_degree = NA)
      if (return_full) {
        result$density <- NA; result$transitivity <- NA
        result$graph <- NULL; result$edge_list <- NULL; result$centrality <- NULL
        result$membership <- NULL; result$species_in_network <- character(0)
      }
      return(result)
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

    result <- list(
      modularity_Q = mod_Q,
      n_modules = n_mod,
      n_species_network = igraph::vcount(g),
      n_edges = igraph::ecount(g),
      mean_degree = mean(igraph::degree(g))
    )

    if (return_full) {
      result$density <- igraph::edge_density(g)
      result$transitivity <- igraph::transitivity(g, type = "global")
      result$graph <- g
      result$edge_list <- edge_list
      # Centrality metrics
      deg <- igraph::degree(g)
      eig <- tryCatch(
        igraph::eigen_centrality(g, weights = igraph::E(g)$weight)$vector,
        error = function(e) setNames(rep(NA, igraph::vcount(g)), igraph::V(g)$name)
      )
      result$centrality <- data.frame(
        species = igraph::V(g)$name,
        degree = deg,
        eigenvector_centrality = eig[igraph::V(g)$name],
        stringsAsFactors = FALSE
      )
      result$membership <- igraph::membership(comm_det)
      result$species_in_network <- igraph::V(g)$name
    }

    result
  }, error = function(e) {
    cat("    Network error:", e$message, "\n")
    if (return_full) empty_full else empty_summary
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
                     error = function(e) {
                       cat("      NOTE: confint() failed for", sp, ", using Wald CI\n")
                       beta + c(-1.96, 1.96) * se
                     })

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
network_objects <- list()

for (scenario_name in scenario_names) {
  scen_data <- scenario_data[[scenario_name]]

  cat("------------------------------------------------------------\n")
  cat("SCENARIO:", scenario_name,
      "| Records:", scen_data$n_records,
      "| OTUs:", scen_data$n_otus,
      "| Matrix:", nrow(scen_data$community_matrix), "x", ncol(scen_data$community_matrix), "\n")
  cat("------------------------------------------------------------\n\n")

  # Use pre-built objects
  comm_mat     <- scen_data$community_matrix
  metrics      <- scen_data$metrics
  cafi_modified <- scen_data$cafi_modified

  cat("  Running analyses...\n")

  # A. Abundance scaling
  cat("    [A] Abundance scaling (NB GLM)... ")
  abund_result <- fit_abundance_scaling(metrics)
  cat("beta =", round(abund_result$beta, 3), "\n")

  # B. Richness scaling
  cat("    [B] Richness scaling (Poisson GLM)... ")
  rich_result <- fit_richness_scaling(metrics)
  cat("z =", round(rich_result$beta, 3), "\n")

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
  net_result <- fit_network(comm_mat, coral_master, return_full = TRUE)
  network_objects[[scenario_name]] <- net_result
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
    n_records = scen_data$n_records,
    n_otus = scen_data$n_otus,
    n_corals_with_cafi = n_distinct(cafi_modified$coral_id),
    n_comm_species = ncol(comm_mat),

    abundance_beta = abund_result$beta,
    abundance_se = abund_result$se,
    abundance_ci_lower = abund_result$ci_lower,
    abundance_ci_upper = abund_result$ci_upper,
    abundance_p = abund_result$p,

    richness_z = rich_result$beta,
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
# NETWORK SENSITIVITY DEEP ANALYSIS
# ============================================================================
# Detailed comparison of co-occurrence networks across taxonomy scenarios:
#   A. Topology comparison (density, transitivity, etc.)
#   B. Edge overlap (Jaccard similarity between scenario networks)
#   C. Hub stability (do the same species remain central?)
#   D. Degree distribution comparison
# ============================================================================

cat("============================================================\n")
cat("NETWORK SENSITIVITY DEEP ANALYSIS\n")
cat("============================================================\n\n")

# --- A. Topology comparison table ---
cat("  [A] Building topology comparison table...\n")

topology_df <- bind_rows(lapply(scenario_names, function(sc) {
  nr <- network_objects[[sc]]
  tibble(
    scenario = sc,
    n_species = nr$n_species_network %||% NA_integer_,
    n_edges = nr$n_edges %||% NA_integer_,
    density = nr$density %||% NA_real_,
    modularity_Q = nr$modularity_Q %||% NA_real_,
    n_modules = nr$n_modules %||% NA_integer_,
    mean_degree = nr$mean_degree %||% NA_real_,
    transitivity = nr$transitivity %||% NA_real_
  )
}))

save_table(topology_df, "network_topology_sensitivity")
cat("    Saved: network_topology_sensitivity.csv\n")

# Print topology comparison
cat("\n    Network Topology Across Scenarios:\n")
cat(sprintf("    %-15s %5s %5s %6s %6s %4s %6s %6s\n",
            "Scenario", "Nodes", "Edges", "Dens.", "ModQ", "nMod", "MnDeg", "Trans."))
cat("    ", paste(rep("-", 60), collapse = ""), "\n")
for (i in seq_len(nrow(topology_df))) {
  r <- topology_df[i, ]
  cat(sprintf("    %-15s %5s %5s %6s %6s %4s %6s %6s\n",
              r$scenario,
              ifelse(is.na(r$n_species), "NA", as.character(r$n_species)),
              ifelse(is.na(r$n_edges), "NA", as.character(r$n_edges)),
              ifelse(is.na(r$density), "NA", sprintf("%.3f", r$density)),
              ifelse(is.na(r$modularity_Q), "NA", sprintf("%.3f", r$modularity_Q)),
              ifelse(is.na(r$n_modules), "NA", as.character(r$n_modules)),
              ifelse(is.na(r$mean_degree), "NA", sprintf("%.1f", r$mean_degree)),
              ifelse(is.na(r$transitivity), "NA", sprintf("%.3f", r$transitivity))))
}
cat("\n")

# --- B. Edge overlap (Jaccard similarity) ---
cat("  [B] Computing edge overlap (Jaccard similarity)...\n")

# Create edge set strings for each scenario
edge_sets <- lapply(scenario_names, function(sc) {
  el <- network_objects[[sc]]$edge_list
  if (is.null(el) || nrow(el) == 0) return(character(0))
  # Sort species names within each edge for consistent comparison
  apply(el[, c("sp1", "sp2")], 1, function(x) paste(sort(x), collapse = " -- "))
})
names(edge_sets) <- scenario_names

# Compute pairwise Jaccard
n_sc <- length(scenario_names)
jaccard_mat <- matrix(NA_real_, nrow = n_sc, ncol = n_sc,
                      dimnames = list(scenario_names, scenario_names))

for (i in seq_len(n_sc)) {
  for (j in seq_len(n_sc)) {
    s_i <- edge_sets[[i]]
    s_j <- edge_sets[[j]]
    if (length(s_i) == 0 && length(s_j) == 0) {
      jaccard_mat[i, j] <- NA_real_
    } else if (length(s_i) == 0 || length(s_j) == 0) {
      jaccard_mat[i, j] <- 0
    } else {
      inter <- length(base::intersect(s_i, s_j))
      union_n <- length(base::union(s_i, s_j))
      jaccard_mat[i, j] <- inter / union_n
    }
  }
}

jaccard_df <- as.data.frame(jaccard_mat) %>%
  mutate(scenario = rownames(jaccard_mat)) %>%
  dplyr::select(scenario, everything())
save_table(jaccard_df, "network_edge_overlap")
cat("    Saved: network_edge_overlap.csv\n")

cat("    Jaccard edge overlap:\n")
for (i in seq_len(n_sc)) {
  for (j in i:n_sc) {
    if (i != j) {
      cat(sprintf("      %s vs %s: J = %.3f\n",
                  scenario_names[i], scenario_names[j],
                  jaccard_mat[i, j]))
    }
  }
}
cat("\n")

# --- C. Hub stability ---
cat("  [C] Analyzing hub species stability...\n")

# Collect centrality from all scenarios
hub_list <- list()
for (sc in scenario_names) {
  cent <- network_objects[[sc]]$centrality
  if (!is.null(cent) && nrow(cent) > 0) {
    cent$scenario <- sc
    hub_list[[sc]] <- cent
  }
}

if (length(hub_list) > 0) {
  hub_all <- bind_rows(hub_list)

  # Get top 10 species from baseline by eigenvector centrality
  baseline_hubs <- hub_all %>%
    filter(scenario == "Baseline") %>%
    arrange(desc(eigenvector_centrality)) %>%
    head(10)

  cat("    Top 10 hub species (Baseline, by eigenvector centrality):\n")
  for (i in seq_len(nrow(baseline_hubs))) {
    sp <- baseline_hubs$species[i]
    # Count how many scenarios this species appears in as a network node
    n_appear <- sum(sapply(scenario_names, function(sc) {
      sp %in% (network_objects[[sc]]$species_in_network %||% character(0))
    }))
    cat(sprintf("      %2d. %-30s eig=%.3f  deg=%d  in %d/%d scenarios\n",
                i, sp, baseline_hubs$eigenvector_centrality[i],
                baseline_hubs$degree[i], n_appear, length(scenario_names)))
  }

  save_table(hub_all, "network_hub_stability")
  cat("    Saved: network_hub_stability.csv\n\n")
} else {
  hub_all <- tibble()
  cat("    No networks produced centrality data.\n\n")
}

# --- D. Degree distribution data ---
cat("  [D] Extracting degree distributions...\n")

degree_list <- list()
for (sc in scenario_names) {
  g <- network_objects[[sc]]$graph
  if (!is.null(g)) {
    deg <- igraph::degree(g)
    degree_list[[sc]] <- tibble(
      scenario = sc,
      species = names(deg),
      degree = as.integer(deg)
    )
  }
}

if (length(degree_list) > 0) {
  degree_all <- bind_rows(degree_list)
  cat("    Degree distributions extracted for", length(degree_list), "scenarios\n\n")
} else {
  degree_all <- tibble()
  cat("    No degree data available.\n\n")
}

# ============================================================================
# FIGURE S14: Network Sensitivity (4-panel)
# ============================================================================

cat("============================================================\n")
cat("GENERATING NETWORK SENSITIVITY FIGURE (S14)\n")
cat("============================================================\n\n")

scenario_colors <- c(
  "Baseline" = "gray30",
  "Species-only" = "#E69F00",
  "Merge-up" = "#56B4E9",
  "Lump-down" = "#009E73",
  "Rare-excluded" = "#CC79A7"
)

# --- Panel A: Topology comparison (faceted bars) ---
topo_long <- topology_df %>%
  dplyr::select(scenario, n_species, n_edges, mean_degree) %>%
  tidyr::pivot_longer(-scenario, names_to = "metric", values_to = "value") %>%
  mutate(
    scenario = factor(scenario, levels = scenario_names),
    metric = factor(metric,
                    levels = c("n_species", "n_edges", "mean_degree"),
                    labels = c("Species", "Edges", "Mean degree"))
  )

p_topo <- ggplot(topo_long, aes(x = scenario, y = value, fill = scenario)) +
  geom_col(color = "gray30", linewidth = 0.3) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = scenario_colors, guide = "none") +
  labs(x = NULL, y = NULL, title = "A. Network topology") +
  theme_publication(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 9),
    plot.margin = margin(5, 5, 5, 5, "mm")
  )

# --- Panel B: Edge Jaccard heatmap ---
jaccard_long <- jaccard_mat %>%
  as.data.frame() %>%
  mutate(scenario1 = factor(rownames(.), levels = scenario_names)) %>%
  tidyr::pivot_longer(-scenario1, names_to = "scenario2", values_to = "jaccard") %>%
  mutate(scenario2 = factor(scenario2, levels = scenario_names))

p_jaccard <- ggplot(jaccard_long, aes(x = scenario1, y = scenario2, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(jaccard), "", sprintf("%.2f", jaccard))),
            size = 3, color = "black") +
  scale_fill_viridis_c(option = "D", na.value = "grey90", limits = c(0, 1),
                       name = "Jaccard", breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = NULL, y = NULL, title = "B. Edge overlap (Jaccard)") +
  theme_publication(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.35, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(5, 5, 5, 5, "mm")
  )

# --- Panel C: Hub stability dot plot ---
if (nrow(hub_all) > 0 && nrow(baseline_hubs) > 0) {
  # Show top 10 baseline hubs across all scenarios
  target_species <- baseline_hubs$species

  hub_plot_data <- hub_all %>%
    filter(species %in% target_species) %>%
    mutate(
      species = factor(species, levels = rev(target_species)),
      scenario = factor(scenario, levels = scenario_names)
    )

  p_hubs <- ggplot(hub_plot_data, aes(x = eigenvector_centrality, y = species,
                                       color = scenario)) +
    geom_point(shape = 21, aes(fill = scenario), size = 2.5, stroke = 0.4,
               position = position_dodge(width = 0.6)) +
    scale_color_manual(values = scenario_colors, name = "Scenario") +
    scale_fill_manual(values = scenario_colors, name = "Scenario") +
    labs(x = "Eigenvector centrality", y = NULL,
         title = "C. Hub species stability") +
    theme_publication(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 8, face = "italic"),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
} else {
  p_hubs <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient network data") +
    labs(title = "C. Hub species stability") +
    theme_void()
}

# --- Panel D: Degree distribution ---
if (nrow(degree_all) > 0) {
  degree_plot_data <- degree_all %>%
    mutate(scenario = factor(scenario, levels = scenario_names))

  p_degree <- ggplot(degree_plot_data, aes(x = degree, fill = scenario, color = scenario)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    scale_fill_manual(values = scenario_colors, name = "Scenario") +
    scale_color_manual(values = scenario_colors, name = "Scenario") +
    labs(x = "Node degree", y = "Density",
         title = "D. Degree distributions") +
    theme_publication(base_size = 10) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
} else {
  p_degree <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Insufficient network data") +
    labs(title = "D. Degree distributions") +
    theme_void()
}

# (S14 figure removed — network diagnostics retained for internal use only)

# ============================================================================
# FIGURE: Forest-Plot Style Comparison (Fig S6)
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

# scenario_colors already defined above (line 814); reuse it here

p_a <- ggplot(panel_a_data, aes(x = beta, y = scenario, color = scenario)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.4) +
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
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray45", linewidth = 0.4) +
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
    title = "Figure S6: Sensitivity of Key Results to Taxonomic Resolution",
    subtitle = paste0("5 scenarios tested across ", nrow(coral_master), " corals, ",
                      nrow(cafi_clean), " CAFI records"),
    theme = theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  )

save_figure(p_combined, file.path(fig_dir, "figS6_taxonomy_sensitivity.png"),
            width = 10, height = 7)
cat("  Saved: figS6_taxonomy_sensitivity.png\n\n")

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
cat("  Table: output/tables/network_topology_sensitivity.csv\n")
cat("  Table: output/tables/network_edge_overlap.csv\n")
cat("  Table: output/tables/network_hub_stability.csv\n")
cat("  Figure: output/figures/supplement/figS6_taxonomy_sensitivity.png\n\n")

cat("Script 13 complete.\n\n")
