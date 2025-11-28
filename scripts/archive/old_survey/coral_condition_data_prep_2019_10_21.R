rm(list=ls())

### Coral Condition Data Preliminary Analysis - Stier Lab Moorea Summer 2019
### Creation Date: October 21, 2019
### Author: Joseph Curtis, UCSB

# Load libraries
library(here)
library(tidyverse)
library(vegan)
library(reshape)

# Load data

afdm <- read.csv(here("laboratory/tissue_removal/data/ash_free_dry_mass_moorea_summer2019.csv"))
slurry_vol <- read.csv(here("laboratory/tissue_removal/data/coral_slurry_data_moorea_summer2019.csv"))
fe_poc_size <- read.csv(here("experiment/size_and_removal_maatea/data/field_experiment_colony_measurements_moorea_summer2019.csv"))
survey_poc_size <- read.csv(here("survey/data/survey_coral_characteristics_moorea_summer2019.csv"))
wax_surface_area <- read.csv(here("laboratory/tissue_removal/data/wax_dipping_data_moorea_summer2019.csv"))
neighborhood <- read.csv(here("survey/data/neighborhood_survey_data_summer2019.csv"))
cafi_surveys <- read.csv(here("survey/data/revised_cafi_data_moorea_summer2019_11_27.csv"))
cafi_field <- read.csv(here("experiment/size_and_removal_maatea/data/prelim_cafi_counts_moorea_summer2019.csv"))
nub_size <- read.csv(here("survey/data/coral_nub_lengths_summer2019.csv"))
## Data cleaning and tissue thickness calculations 

# Exclude "backup" or unusable samples from AFDM, surface area, and slurry volume datasets

afdm <- subset(afdm, use == "y")
slurry_vol <- subset(slurry_vol, use == "y")
wax_surface_area <- subset(wax_surface_area, use == "y")
nub_size <- subset(nub_size, use == "y")

# Average replicate AFDW measurements

afdm_avg <- add_count(afdm, coral_id) %>%
  group_by(coral_id) %>%
  summarize(afdm_g_ml = mean(afdw_g_ml), rep_sd = sd(afdw_g_ml), rep_num = max(n)) %>% 
  replace(is.na(.),0)

# Calculate total mass of tissue in slurry
afdm_w_vol <- full_join(afdm_avg, slurry_vol, by = "coral_id") %>% select(coral_id, afdm_g_ml, volume, notes)
afdm_w_vol$tissue_vol_tot <- afdm_w_vol$afdm_g_ml*afdm_w_vol$volume

# Convert surface area of nubs to cm2 from mm2

wax_surface_area$sa_cm2 <- wax_surface_area$surface_area*0.01
thickness_calc <- full_join(afdm_w_vol, wax_surface_area, by = "coral_id") %>% select(coral_id, tissue_vol_tot, sa_cm2, nub, notes.x, notes.y)
thickness_calc$thickness_g_cm2 <- thickness_calc$tissue_vol_tot/thickness_calc$sa_cm2

# Condense neighborhood data

neighborhood_sum <- group_by(neighborhood, coral_id) %>%
  summarise(obs = first(obs), num_neighbors = n(),
            avg_dist = mean(distance_cm), closest = min(distance_cm), 
            volume_tot = sum(volume_tot), volume_alive = sum(volume_alive)) %>%
  mutate(num_neighbors = replace(num_neighbors, volume_tot == "0", 0)) %>%
  replace(is.na(.),250)

# Combine size data from surveys and experiment with physio and neighborhood data

fe_poc_size$survey_type <- "experiment"
comb_survey_data <- bind_rows(survey_poc_size,fe_poc_size) %>%
  full_join(thickness_calc, by = "coral_id") %>%
  select(-lat, -long, -field_obs, -ends_with("lab"), -position, -size_class, -site) %>%
  separate(coral_id, c("reef", NA), remove = FALSE) %>%
  full_join(neighborhood_sum, by = "coral_id")

comb_survey_data <- comb_survey_data %>%  mutate(max_dim = pmax(comb_survey_data$length_field, comb_survey_data$width_field, comb_survey_data$height_field)) %>%
  dplyr::rename(coral_volume = volume_field) %>% 
  select(-ends_with("field"))

# Add "corrected" nub size based on relative length to total branch

nub_size <- nub_size %>% unite(col_w_nub, coral_id, sample_number, sep = "_nub_")
comb_survey_data <- comb_survey_data %>% unite(col_w_nub, coral_id, nub, sep = "_nub_", remove = "FALSE")

nub_size$total_branch <- nub_size$nubbin_length + nub_size$stump_length
nub_size$prop_nub <- nub_size$nubbin_length/nub_size$total_branch
comb_data_w_nub <- inner_join(comb_survey_data, nub_size, by = "col_w_nub")

### Incorporate CAFI data

# Uncount CAFI data from field experiment

cafi_field_long <- uncount(cafi_field, count) %>% select(order, coral_id, type, code, cafi_size_mm, known_unknown) %>%
  separate(coral_id, c("reef", NA), remove = FALSE)

cafi_surveys <- cafi_surveys %>% select(sort, coral_id, type, code, cafi_size_mm, known_unknown) %>%
  separate(coral_id, c("reef", NA), remove = FALSE)

# While identification and data cleaning is still underway, we will only analyze known organisms

cafi_combined <- bind_rows(cafi_surveys, cafi_field_long) %>% 
  filter(known_unknown == "known")

# Calculate CAFI richness, abundance, and diversity (shannon weiner) for each coral
cafi_summarized <- group_by(cafi_combined, coral_id) %>%
  summarise(num_cafi = n(), cafi_richness = length(unique(code)), cafi_present = paste(sort(unique(code)), collapse = ";"))

cafi_summarized$sw <- cafi_combined %>% 
  count(code, coral_id = coral_id) %>% 
  spread(code,n) %>% 
  mutate_all(list(~replace_na(.,0))) %>% 
  select(-coral_id) %>% 
  diversity(index = "shannon")

# Calculate Trapezid abundance and richness for each coral

its_a_trap <- cafi_combined %>% filter(str_detect(code, "^TR"))
trap_summarized <- group_by(its_a_trap, coral_id) %>%
  summarise(num_cafi = n(), cafi_richness = length(unique(code)), cafi_present = paste(sort(unique(code)), collapse = ";"))

# Combine with physio data, add 0s for corals with no crabs
trap_w_physio <- trap_summarized %>% full_join(comb_survey_data, by = "coral_id") %>% 
  filter(cafi != "CAFI" | survey_type != "experiment")

trap_w_physio$num_trap <- trap_w_physio$num_cafi %>% replace_na(0) %>% rename("num_trap")
trap_w_physio$trap_richness <- trap_w_physio$cafi_richness %>% replace_na(0) %>% rename("trap_richness")
trap_w_physio <- trap_w_physio %>% select(-cafi_richness, -num_cafi)
  
# Calculate CAFI richness and abundance for each coral by CAFI type

cafi_sum_type <- cafi_combined %>% group_by(coral_id, type) %>%
  summarise(num_cafi = n(), cafi_richness = length(unique(code)), cafi_present = paste(sort(unique(code)), collapse = ";")) %>%
  full_join(comb_survey_data, by = "coral_id")


# Merge with physio and morphological data

cafi_w_physio <- cafi_summarized %>% full_join(comb_survey_data, by = "coral_id")


# Separate survey types for future analyses
neighborhood_only <- filter(cafi_w_physio, survey_type == "neighborhood")
size_only <- filter(cafi_w_physio, survey_type == "size")
experiment_only <- filter(comb_survey_data, survey_type == "experiment")

# Extract corals from neighborhood surveys to include in size surveys where neighborhood was "medium" (between 1st and 3rd quantiles)

size_complete <- neighborhood_only %>% filter(volume_alive >=(summary(volume_alive)[[2]]) & volume_alive <= (summary(volume_alive)[[5]])) %>%
  bind_rows(size_only)

# Clean up environment

rm(afdm, afdm_avg, afdm_w_vol, cafi_field, cafi_field_long, cafi_summarized, cafi_surveys, fe_poc_size, neighborhood_sum, slurry_vol, survey_poc_size, thickness_calc, wax_surface_area)

######################################################################
###Visualizations

# Plot coral tissue biomass with colony volume

ggplot(comb_survey_data, aes(coral_volume,thickness_g_cm2)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_classic() +
  ggtitle("Coral Tissue Biomass vs Coral Volume") + 
  xlab(bquote("Coral Volume" ~cm^3)) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

ggplot(comb_survey_data, aes(coral_volume,thickness_g_cm2)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  facet_wrap(~survey_type) +
  theme_classic() +
  ggtitle("Coral Tissue Biomass vs Coral Volume ~survey") + 
  xlab(bquote("Coral Volume" ~cm^3)) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

ggplot(comb_survey_data, aes(coral_volume,thickness_g_cm2)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  facet_wrap(~branch_width) +
  theme_classic() +
  ggtitle("Coral Tissue Biomass vs Coral Volume ~Branch Width") + 
  xlab(bquote("Coral Volume" ~cm^3)) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

# Plot coral tissue biomass with coral type

ggplot(comb_survey_data, aes(branch_width, thickness_g_cm2)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("Coral Tissue Biomass vs Coral Volume ~Branch Width") + 
  xlab(bquote("Branch Width")) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))


# Plot coral biomass vs nub size

ggplot(comb_survey_data, aes(sa_cm2,thickness_g_cm2)) +
  geom_point() +
  geom_smooth() +
  theme_classic() +
  ggtitle("Coral Nub Surface Area vs. Coral Tissue Biomass") + 
  xlab(bquote("Nub Surface Area"~(cm^2))) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

ggplot(comb_survey_data, aes(sa_cm2,thickness_g_cm2)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~branch_width) +
  theme_classic() +
  ggtitle("Coral Nub Surface Area vs. Coral Tissue Biomass ~Branch Width") + 
  xlab(bquote("Nub Surface Area"~(cm^2))) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

ggplot(comb_survey_data, aes(sa_cm2,thickness_g_cm2)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~survey_type) +
  theme_classic() +
  ggtitle("Coral Nub Surface Area vs. Coral Tissue Biomass ~Survey Type") + 
  xlab(bquote("Nub Surface Area"~(cm^2))) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

# Plot volume vs nub size
ggplot(comb_survey_data, aes(coral_volume, sa_cm2)) +
  geom_point() +
  geom_smooth() +
  theme_classic() +
  ggtitle("Coral Volume vs. Nub Surface Area") + 
  xlab(bquote("Coral Volume" ~cm^3)) +
  ylab(bquote("Nub Surface Area"~(cm^2)))

# Plot volume vs nub size, faceted by branch width
ggplot(comb_survey_data, aes(coral_volume, sa_cm2)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  facet_wrap(~branch_width) +
  theme_classic() +
  ggtitle("Coral Volume vs. Nub Surface Area ~Branch Width") + 
  xlab(bquote("Coral Volume" ~cm^3)) +
  ylab(bquote("Nub Surface Area"~(cm^2)))

# Plot volume vs nub size, faceted by survey type
ggplot(comb_survey_data, aes(coral_volume, sa_cm2)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  facet_wrap(~survey_type) +
  theme_classic() +
  ggtitle("Coral Volume vs. Nub Surface Area ~Survey Type") + 
  xlab(bquote("Coral Volume" ~cm^3)) +
  ylab(bquote("Nub Surface Area"~(cm^2)))

### CAFI Plots

# Coral size vs richness, abundance, diversity

# ggplot(cafi_w_physio, aes(coral_volume, num_cafi)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   theme_classic() +
#   ggtitle("Coral Volume vs CAFI Abundance") + 
#   xlab(bquote("Coral Volume" ~cm^3)) +
#   ylab("CAFI Abundance")

ggplot(cafi_w_physio, aes(coral_volume, num_cafi)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  facet_wrap(~branch_width) +
  ggtitle("Coral Volume vs CAFI Abundance~ Branch Width") + 
  xlab(bquote("Coral Volume" ~cm^3)) +
  ylab("CAFI Abundance")

# CAFI Abundance ~ Branch Width

# ggplot(cafi_w_physio, aes(branch_width, num_cafi)) +
#  geom_boxplot() +
#  theme_classic() +
#  ggtitle("CAFI Abundance~Branch Width") + 
#  xlab(bquote("Branch Width")) +
#  ylab("CAFI Abundance")

# ggplot(cafi_w_physio, aes(coral_volume, cafi_richness)) +
#   geom_point() +
#   geom_smooth(method = "lm")  +
#   theme_classic() +
#   ggtitle("Coral Volume vs CAFI Richness") + 
#   xlab(bquote("Coral Volume" ~cm^3)) +
#   ylab("CAFI Richness")

ggplot(cafi_w_physio, aes(coral_volume, cafi_richness)) +
  geom_point() +
  geom_smooth(method = "lm")  +
  theme_classic() +
  facet_wrap(~branch_width) +
  ggtitle("Coral Volume vs CAFI Richness ~Branch Width") + 
  xlab(bquote("Coral Volume" ~cm^3)) +
  ylab("CAFI Richness")

# ggplot(cafi_w_physio, aes(coral_volume, sw)) +
#   geom_point() +
#   geom_smooth(method = "lm")  +
#   theme_classic() +
#   ggtitle("Coral Volume vs CAFI Diversity (SW)") + 
#   xlab(bquote("Coral Volume" ~cm^3)) +
#   ylab("CAFI Diversity (SW Index)")

ggplot(cafi_w_physio, aes(coral_volume, sw)) +
  geom_point() +
  geom_smooth(method = "lm")  +
  theme_classic() +
  facet_wrap(~branch_width) +
  ggtitle("Coral Volume vs CAFI Diversity (SW)") + 
  xlab(bquote("Coral Volume" ~cm^3)) +
  ylab("CAFI Diversity (SW Index)")

ggplot(cafi_w_physio, aes(reef, num_cafi)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("CAFI Abundance~Reef") + 
  xlab(bquote("Reef")) +
  ylab("CAFI Abundance")

ggplot(cafi_w_physio, aes(reef, cafi_richness)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("CAFI Richness~Reef") + 
  xlab(bquote("Reef")) +
  ylab("CAFI Richness")

ggplot(cafi_w_physio, aes(reef, sw)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("CAFI Diversity~Reef") + 
  xlab(bquote("Reef")) +
  ylab("CAFI Diversity (SW Index)")

#Plot CAFI abundance vs thickness
ggplot(cafi_w_physio, aes(num_cafi, thickness_g_cm2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggtitle("CAFI Abundance vs Coral Tissue Biomass") + 
  xlab(bquote("CAFI Abundance")) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

## Plot abundance, richness, and diversity of CAFI vs coral volume in consistent, medium neighborhoods
# ggplot(size_complete, aes(coral_volume, num_cafi)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~branch_width) 
# 
# ggplot(size_complete, aes(coral_volume, cafi_richness)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~branch_width)
# 
# ggplot(size_complete, aes(coral_volume, sw)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~branch_width)

## Same plots as above, but with trapezids only

ggplot(trap_w_physio, aes(coral_volume, num_trap)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  facet_wrap(~branch_width) +
  ggtitle("Coral Volume vs trap Abundance~ Branch Width") + 
  xlab(bquote("Coral Volume" ~cm^3)) +
  ylab("trap Abundance")

# trap Abundance ~ Branch Width

# ggplot(trap_w_physio, aes(branch_width, num_trap)) +
#  geom_boxplot() +
#  theme_classic() +
#  ggtitle("trap Abundance~Branch Width") + 
#  xlab(bquote("Branch Width")) +
#  ylab("trap Abundance")

# ggplot(trap_w_physio, aes(coral_volume, trap_richness)) +
#   geom_point() +
#   geom_smooth(method = "lm")  +
#   theme_classic() +
#   ggtitle("Coral Volume vs trap Richness") + 
#   xlab(bquote("Coral Volume" ~cm^3)) +
#   ylab("trap Richness")

ggplot(trap_w_physio, aes(coral_volume, trap_richness)) +
  geom_point() +
  geom_smooth(method = "lm")  +
  theme_classic() +
  facet_wrap(~branch_width) +
  ggtitle("Coral Volume vs trap Richness ~Branch Width") + 
  xlab(bquote("Coral Volume" ~cm^3)) +
  ylab("trap Richness")

ggplot(trap_w_physio, aes(reef, num_trap)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("trap Abundance~Reef") + 
  xlab(bquote("Reef")) +
  ylab("trap Abundance")

ggplot(trap_w_physio, aes(reef, trap_richness)) +
  geom_boxplot() +
  theme_classic() +
  ggtitle("trap Richness~Reef") + 
  xlab(bquote("Reef")) +
  ylab("trap Richness")


#Plot trap abundance vs thickness
ggplot(trap_w_physio, aes(num_trap, thickness_g_cm2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggtitle("trap Abundance vs Coral Tissue Biomass") + 
  xlab(bquote("trap Abundance")) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

## Plot abundance, richness, and diversity of CAFI vs neighborhood characteristics

# vs Volume of live coral in 2m radius

ggplot(neighborhood_only, aes(volume_alive, num_cafi)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggtitle("CAFI Abundance vs Volume of Live Coral in 2m Radius") + 
  xlab(bquote("Volume of Live Coral in 2m Radius")) +
  ylab(bquote("CAFI Abundance"))

ggplot(neighborhood_only, aes(volume_alive, cafi_richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggtitle("CAFI Richness vs Volume of Live Coral in 2m Radius") + 
  xlab(bquote("Volume of Live Coral in 2m Radius")) +
  ylab(bquote("CAFI Richness"))

ggplot(neighborhood_only, aes(volume_alive, sw)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggtitle("CAFI Diversity vs Volume of Live Coral in 2m Radius") + 
  xlab(bquote("Volume of Live Coral in 2m Radius")) +
  ylab(bquote("CAFI Diversity (SW)"))

# vs Volume of live coral in 2m radius faceted by type

# ggplot(cafi_sum_type, aes(volume_alive, num_cafi)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~type)
# 
# ggplot(cafi_sum_type, aes(volume_alive, cafi_richness)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~type)

# vs number of corals in 2m radius

ggplot(neighborhood_only, aes(num_neighbors, num_cafi)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggtitle("CAFI Abundance vs Number of Corals in 2m Radius") + 
  xlab(bquote("Number of Corals in 2m Radius")) +
  ylab(bquote("CAFI Abundance"))

ggplot(neighborhood_only, aes(num_neighbors, cafi_richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggtitle("CAFI Richness vs Number of Corals in 2m Radius") + 
  xlab(bquote("Number of Corals in 2m Radius")) +
  ylab(bquote("CAFI Richness"))

ggplot(neighborhood_only, aes(num_neighbors, sw)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggtitle("CAFI Diversity vs Number of Corals in 2m Radius") + 
  xlab(bquote("Number of Corals in 2m Radius")) +
  ylab(bquote("CAFI Diversity (SW)"))

# vs number of corals in 2m radius - faceted by CAFI type

# ggplot(cafi_sum_type, aes(num_neighbors, num_cafi)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~type)
# 
# ggplot(cafi_sum_type, aes(num_neighbors, cafi_richness)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~type)

# vs minimum distance to nearest coral  

# ggplot(neighborhood_only, aes(closest, num_cafi)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(neighborhood_only, aes(closest, cafi_richness)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(neighborhood_only, aes(closest, sw)) +
#   geom_point() +
#   geom_smooth(method = "lm")

# vs minimum distance to nearest coral faceted by type 

# ggplot(cafi_sum_type, aes(closest, num_cafi)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~type)
# 
# ggplot(cafi_sum_type, aes(closest, cafi_richness)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~type)

# vs average distance to nearest coral  
# 
# ggplot(neighborhood_only, aes(avg_dist, num_cafi)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(neighborhood_only, aes(avg_dist, cafi_richness)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# 
# ggplot(neighborhood_only, aes(avg_dist, sw)) +
#   geom_point() +
#   geom_smooth(method = "lm")

# vs average distance to nearest coral faceted by type 

# ggplot(cafi_sum_type, aes(avg_dist, num_cafi)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~type)
# 
# ggplot(cafi_sum_type, aes(avg_dist, cafi_richness)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~type)
# 
# 
# ggplot(neighborhood_only, aes(avg_dist, sw)) +
#   geom_point() +
#   geom_smooth(method = "lm")

# CAFI Abundance per reef faceted by type

#ggplot(cafi_sum_type, aes(reef, num_cafi)) +
 # geom_boxplot() +
  #theme_classic() +
  #ggtitle("CAFI Abundance~Reef") + 
  #facet_wrap(~type) +
  #xlab(bquote("Reef")) +
  #ylab("CAFI Abundance")

#### Physio Plots with nub size corrected

ggplot(comb_data_w_nub, aes(total_branch,thickness_g_cm2)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_classic() +
  ggtitle("Coral Tissue Biomass vs Total Branch Length") + 
  xlab(bquote("Branch Length (mm)")) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

ggplot(comb_data_w_nub, aes(nubbin_length,thickness_g_cm2)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_classic() +
  ggtitle("Coral Tissue Biomass vs Nub Length") + 
  xlab(bquote("Nub Length (mm)")) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

ggplot(comb_data_w_nub, aes(stump.length,thickness_g_cm2)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_classic() +
  ggtitle("Coral Tissue Biomass vs Stump Length") + 
  xlab(bquote("Nub Length (mm)")) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

ggplot(comb_data_w_nub, aes(nubbin_length,coral_volume)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_classic() +
  ggtitle("Nub Length vs Coral Volume") + 
  xlab(bquote("Nub Length (mm)")) +
  ylab(bquote("Coral Volume" ~(cm^3)))

ggplot(comb_data_w_nub, aes(nubbin_length,max_dim)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_classic() +
  ggtitle("Nub Length vs Max Linear Dimension") + 
  xlab(bquote("Nub Length (mm)")) +
  ylab(bquote("Largest Linear Dimension (cm)"))

ggplot(comb_data_w_nub, aes(max_dim, thickness_g_cm2)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_classic() +
  ggtitle("Maximum Linear Dimension vs Tissue Biomass") + 
  xlab(bquote("Largest Linear Dimension (cm)")) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

ggplot(comb_data_w_nub, aes(prop_nub, coral_volume)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_classic() +
  ggtitle("Maximum Linear Dimension vs Tissue Biomass") + 
  xlab(bquote("Largest Linear Dimension (cm)")) +
  ylab(bquote("Coral Tissue Biomass" ~g/cm^2))

