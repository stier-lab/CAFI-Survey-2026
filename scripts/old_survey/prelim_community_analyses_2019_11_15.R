rm(list=ls())

### CAFI Community Analysis - Stier Lab Moorea Summer 2019
### Creation Date: November 15, 2019
### Author: Joseph Curtis, UCSB

# Load libraries
library(here)
library(tidyverse)
library(vegan)
library(reshape)
library(psych)

# Load data

# Load cleaned and treated data, remove unnecessary 

#source("survey/code/coral_condition_data_prep_2019_10_21.R")

cafi_combined <- read_csv("~/Google Drive/Stier Lab/People/Kai Kopecky/Projects/CAFI/survey/old/prelim_cafi_data_surveys_summer2019.csv")

######################
# Community analyses

# Turn long form data into community matrix, create clean metadata key
cafi_matrix <- cafi_combined %>% 
  count(code, coral_id = coral_id) %>% 
  spread(code,n) %>% 
  mutate_all(list(~replace_na(.,0))) %>% 
  column_to_rownames(var = "coral_id")

# Remove rare species

cafi_matrix_trim <- cafi_matrix %>% bind_rows(summarise_each(.,'sum'))
cafi_matrix_trim <- bind_rows(prop.table(cafi_matrix_trim))

cafi_coral_metadat <- cafi_w_physio %>% drop_na(cafi_present) %>% 
  select(-date, -morphotype, -depth, -h_substrate, -cafi, -tissue_vol_tot, -sa_cm2, -obs)

# Combine population matrix and metadata

cafi_matrix_w_metadat <- cafi_matrix <- cafi_combined %>% 
  count(code, coral_id = coral_id) %>% 
  spread(code,n) %>% 
  mutate_all(list(~replace_na(.,0))) %>% 
  full_join(cafi_coral_metadat, by = "coral_id") %>% 
  select(-nub, -col_w_nub, -cafi_present) %>% 
  unite(notes, notes, notes.x, notes.y, sep = ";")
  

# Convert absolute abundance to relative abundance and create dissimilarity matrix using Bray-Curtis

cafi_matrix_dist <- cafi_matrix %>% 
  decostand(method = "total") %>% 
  vegdist(method = "bray") %>% 
  as.matrix(labels = T)

# Run NDMS

cafi_nmds <- metaMDS(cafi_matrix, distance = "bray", k = 2, maxit = 999, trymax = 500, wascores = TRUE, weakties = FALSE)
stressplot(cafi_nmds)
plot(cafi_nmds, "sites")
orditorp(cafi_nmds, "sites")
ordihull(cafi_nmds, groups = cafi_coral_metadat$reef, draw = "polygon", label = T)

# Create clean metadata for all CAFI corals

cafi_coral_metadat <- fe_poc_size %>% 
  filter(cafi == "empty") %>% 
  bind_rows(survey_poc_size) %>% 
  select(-lat, -long, -field_obs, -ends_with("lab"), -position, -size_class, -site, -h_substrate, -cafi, -morphotype, -date, -depth) %>%
  separate(coral_id, c("reef", NA), remove = FALSE)


# Ordination analyses -----------------------------------------------------

## Principal components analysis on species composition 

# Subset data to include only variables to be used in PCA
cafi_species_only <- cafi_combined %>% 
  select(coral_id, reef, code) %>% 
  group_by(coral_id, reef) %>% 
  count(code) %>% 
  pivot_wider(names_from = code, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  rownames(cafi_species_only) <- cafi_species_only$coral_id

# Vectors of only Coral ID's or Reefs
coral_ID <- cafi_species_only$coral_id
reef <- cafi_species_only$reef

# Run PCA
cafi_PCA <- princomp(cafi_matrix, cor = FALSE)

# Plot the principal coordinates
plot(cafi_PCA$scores, pch = 16)
#legend(0, 0.4, c("Field experiment", "Haumi", "Ma'atea", "MRB"), pch = 16, col = c("black", "red", "blue", "grey"))

# Interpretaion
summary(cafi_PCA) # Shows numerically the proportion of variance explained by each principal component
screeplot(cafi_PCA, type = 'lines') # Shows visually the proportion of variance explained by each principal component
loadings(cafi_PCA) # Shows numerically how strongly correlated (negatively or postively) each variable is to each principal component
biplot(cafi_PCA) # Shows visually how strongly correlated each variable is to the first two principal components; longer arrows indicate tighter correlations

################

## Canonical analysis of principal coordinates
# Hypothesis: variation in amount, size, and isolation of habitat will drive variation in species composition

cafi_CAP <- CAPdiscrim(formula, 
                       data = cafi_PCA, 
                       dist = "bray", 
                       axes = 2, 
                       m = 0, 
                       add = FALSE)

###################

### How do different habitat features predict variation in abundance, diversity, and species composition? 
# Create data frame containing only variables associated with coral neighborhood characteristics and CAFI species richness/abundance

cafi_neighborhoods <- cafi_coral_metadat %>% 
  select(-sw, -col_w_nub, -survey_type, -notes, -nub, -notes.x, -notes.y, -thickness_g_cm2) %>% 
  drop_na(num_neighbors, avg_dist, closest, volume_tot, volume_alive)

## Pairwise comparisons

# Number/richness of CAFI as a function of total neighborhood coral volume
ggplot(cafi_neighborhoods, aes(x = volume_tot, y = num_cafi)) +
  geom_point()

ggplot(cafi_neighborhoods, aes(x = volume_tot, y = cafi_richness)) +
  geom_point()

# Number/richness of CAFI as a function of number of neighbors
ggplot(cafi_neighborhoods, aes(x = num_neighbors, y = num_cafi)) +
  geom_point()

ggplot(cafi_neighborhoods, aes(x = num_neighbors, y = cafi_richness)) +
  geom_point()

# Number/richness of CAFI as a function of focal coral volume
ggplot(cafi_neighborhoods, aes(x = coral_volume, y = num_cafi)) +
  geom_point()

ggplot(cafi_neighborhoods, aes(x = coral_volume, y = cafi_richness)) +
  geom_point()


## Linear regression

## Multiple linear regression:

# CAFI richness as a function of focal coral and neighborhood characteristics
cafi_richness_lm <- lm(cafi_richness ~ branch_width + coral_volume + num_neighbors + avg_dist + volume_tot, data = cafi_neighborhoods)

# Summary outputs and table
cafi_richness_lm
summary(cafi_richness_lm)
stargazer(cafi_richness_lm, type = "html")


# CAFI abundance as a function of focal coral and neighborhood characteristics
cafi_abundnace_lm <- lm(num_cafi ~ branch_width + coral_volume + num_neighbors + avg_dist + volume_tot, data = cafi_neighborhoods)

# Summary outputs and table
cafi_abundnace_lm
summary(cafi_abundnace_lm)
stargazer(cafi_abundnace_lm, type = "html")


##########################

### Checking for collinearity of predictor variables (neighborhood characteristics)

# Histograms of predictor variables

hist(cafi_neighborhoods$coral_volume)
hist(cafi_neighborhoods$num_neighbors)
hist(cafi_neighborhoods$avg_dist)
hist(cafi_neighborhoods$closest)
hist(cafi_neighborhoods$volume_tot)
hist(cafi_neighborhoods$volume_alive)

# Pairwise comparisons

neighborhood_metrics <- cafi_neighborhoods %>% 
  select(coral_volume, num_neighbors, avg_dist, closest, volume_tot, volume_alive) %>% 
  drop_na()

pairs.panels(neighborhood_metrics)

# PCA of neighborhood characteristics

neighborhood_PCA <- princomp(neighborhood_metrics, cor = TRUE)
summary(neighborhood_PCA)
str(neighborhood_PCA)

# Plot the principal coordinates
plot(neighborhood_PCA$scores, pch = 16)

# Interpretaion
screeplot(neighborhood_PCA, type = 'lines') # Shows visually the amount of variance explained by each principal component
loadings(neighborhood_PCA) # Shows numerically how strongly correlated (negatively or postively) each variable is to each principal component
biplot(neighborhood_PCA, pc.biplot = TRUE) # Shows visually how strongly correlated each variable is to the first two principal components; longer arrows indicate tighter correlations



