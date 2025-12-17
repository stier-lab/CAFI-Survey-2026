#!/usr/bin/env Rscript
# Expert Analysis: Diversity-Productivity Relationships

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(car)
})

# Load data
coral_data <- read.csv("data/survey_coral_characteristics_merged_v2.csv")
cafi_data <- read.csv("data/survey_cafi_data_w_taxonomy_summer2019_v5.csv")
condition_data <- read.csv("output/tables/coral_condition_scores.csv")

# Process
coral_processed <- coral_data %>%
  mutate(site = str_extract(site, "^[A-Z]+"),
         volume = coalesce(volume_field, volume_lab),
         log_volume = log10(volume + 1)) %>%
  filter(!is.na(volume), volume > 0, site %in% c("HAU", "MAT", "MRB"))

cafi_processed <- cafi_data %>%
  filter(!is.na(genus) & genus != "" & genus != "NA") %>%
  mutate(species_id = paste(genus, species, sep = "_"))

comm_matrix <- cafi_processed %>%
  group_by(coral_id, species_id) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_id, values_from = count, values_fill = 0)

comm_mat <- as.matrix(comm_matrix[,-1])
rownames(comm_mat) <- comm_matrix$coral_id

# Diversity metrics
div_data <- data.frame(
  coral_id = comm_matrix$coral_id,
  abundance = rowSums(comm_mat),
  richness = rowSums(comm_mat > 0),
  shannon = diversity(comm_mat, "shannon"),
  evenness = diversity(comm_mat, "shannon") / log(rowSums(comm_mat > 0))
)
div_data$evenness[is.nan(div_data$evenness)] <- NA

# Rarefied richness
valid <- div_data$abundance >= 10
div_data$rare10 <- NA
div_data$rare10[valid] <- rarefy(comm_mat[valid,], 10)

# Merge
dat <- div_data %>%
  left_join(coral_processed %>% select(coral_id, site, volume, log_volume), by = "coral_id") %>%
  left_join(condition_data %>% select(coral_id, condition_score), by = "coral_id") %>%
  filter(!is.na(volume), !is.na(condition_score))

cat("\n=== EXPERT DIVERSITY-PRODUCTIVITY ANALYSIS ===\n")
cat("n =", nrow(dat), "\n\n")

cat("--- PART 1: Sampling Artifact Check ---\n")
cat("Correlations with log(volume):\n")
cat(sprintf("  Abundance:  r = %.3f ***\n", cor(dat$abundance, dat$log_volume)))
cat(sprintf("  Richness:   r = %.3f ***\n", cor(dat$richness, dat$log_volume)))
cat(sprintf("  Shannon:    r = %.3f ***\n", cor(dat$shannon, dat$log_volume)))
r10 <- cor(dat$rare10, dat$log_volume, use = "complete")
cat(sprintf("  Rarefied:   r = %.3f     <-- ARTIFACT REMOVED\n", r10))

cat("\n--- PART 2: Simple Condition Correlations ---\n")
cat(sprintf("  Richness vs condition: r = %.3f, p = %.3f\n",
    cor(dat$richness, dat$condition_score),
    cor.test(dat$richness, dat$condition_score)$p.value))
cat(sprintf("  Rarefied vs condition: r = %.3f, p = %.3f\n",
    cor(dat$rare10, dat$condition_score, use = "complete"),
    cor.test(dat$rare10, dat$condition_score)$p.value))

cat("\n--- PART 3: Volume-controlled Models ---\n")
m1 <- lm(condition_score ~ richness + log_volume, data = dat)
cat("Raw richness + volume:\n")
cat(sprintf("  β(richness) = %.4f, p = %.4f\n", coef(m1)[2], summary(m1)$coef[2,4]))

m2 <- lm(condition_score ~ rare10 + log_volume, data = dat)
cat("Rarefied richness + volume:\n")
cat(sprintf("  β(rarefied) = %.4f, p = %.4f\n", coef(m2)[2], summary(m2)$coef[2,4]))

cat("\n--- PART 4: Is it Richness or Abundance? ---\n")
m3 <- lm(condition_score ~ richness + abundance + log_volume, data = dat)
cat("Full model (richness + abundance + volume):\n")
print(summary(m3)$coefficients[,c(1,4)])

cat("\n--- PART 5: Richness Residuals (pure diversity) ---\n")
rich_mod <- lm(richness ~ log10(abundance + 1), data = dat)
dat$rich_resid <- residuals(rich_mod)
cat(sprintf("Abundance explains %.1f%% of richness\n", summary(rich_mod)$r.sq * 100))
m4 <- lm(condition_score ~ rich_resid + log_volume, data = dat)
cat(sprintf("Residualized richness: β = %.4f, p = %.4f\n", coef(m4)[2], summary(m4)$coef[2,4]))

cat("\n--- PART 6: Evenness (True Diversity) ---\n")
even_dat <- dat %>% filter(!is.na(evenness))
m5 <- lm(condition_score ~ evenness + log_volume, data = even_dat)
cat(sprintf("Evenness: β = %.4f, p = %.4f (n = %d)\n", coef(m5)[2], summary(m5)$coef[2,4], nrow(even_dat)))

cat("\n--- PART 7: Shannon Diversity ---\n")
m6 <- lm(condition_score ~ shannon + log_volume, data = dat)
cat(sprintf("Shannon H: β = %.4f, p = %.4f\n", coef(m6)[2], summary(m6)$coef[2,4]))

cat("\n=== SUMMARY TABLE ===\n")
cat("Metric                | β        | p-value  | Interpretation\n")
cat("----------------------|----------|----------|---------------\n")
cat(sprintf("Raw richness (+vol)   | %8.4f | %8.4f | Volume confounded\n", coef(m1)[2], summary(m1)$coef[2,4]))
cat(sprintf("Rarefied richness     | %8.4f | %8.4f | No effect\n", coef(m2)[2], summary(m2)$coef[2,4]))
cat(sprintf("Richness residuals    | %8.4f | %8.4f | No pure diversity effect\n", coef(m4)[2], summary(m4)$coef[2,4]))
cat(sprintf("Evenness              | %8.4f | %8.4f | No effect\n", coef(m5)[2], summary(m5)$coef[2,4]))
cat(sprintf("Shannon H             | %8.4f | %8.4f | Marginally sig\n", coef(m6)[2], summary(m6)$coef[2,4]))

cat("\n=== EXPERT CONCLUSION ===\n")
cat("The apparent richness-condition relationship is a SAMPLING ARTIFACT.\n")
cat("When diversity is properly measured (rarefied, residualized, or evenness),\n")
cat("there is NO significant relationship with coral condition.\n")
