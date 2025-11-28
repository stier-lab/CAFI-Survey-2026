# Unique CAFI List Creator

# Author: Joseph Curtis
# Date: November 27, 2019

# Clear environment and load libraries and data

rm(list=ls())
library(here)
library(tidyverse)

survey_cafi <- read.csv(here("survey/data/revised_cafi_data_moorea_summer2019_11_27.csv"))
experiment_cafi <- read.csv(here("experiment/size_and_removal_maatea/data/prelim_cafi_counts_moorea_summer2019.csv"))

survey_cafi_trim <- survey_cafi %>% select(code, known_unknown, type)
experiment_cafi_trim <- experiment_cafi %>% select(code, known_unknown, type)

complete_cafi_list <- bind_rows(survey_cafi_trim, experiment_cafi_trim)

final_list <- complete_cafi_list %>% distinct() %>% arrange(code) %>% arrange (known_unknown) %>% arrange(type)
qa_check <- final_list %>% group_by(code) %>% summarise(count = n()) %>% filter(count>1) #should equal 0 if there are no errors

write.csv(final_list, here("resources/original_cafi_species_list_summer2019.csv"))
