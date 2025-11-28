rm(list=ls())

# Updating CAFI Data from Summer Surveys (MRB, HAU, MAT-POC colonies) with second opinion IDs by JC
# Date: May 7, 2020

# Load libraries, import data

library(here)
library(tidyverse)

cafi_data <- read.csv(here("survey/data/intermediate_cafi_data_before_updated_ids_2020_5_7.csv"))
jc_updates <- read.csv(here("survey/data/ap_identifications_w_JC_second_opinions_2020_5_7.csv"))

# Look for changes to AP IDs by JC

jc_updates$AP_Code <- as.character(jc_updates$AP_Code)
jc_updates$JC_Opinion <- as.character(jc_updates$JC_Opinion)
jc_updates <- add_column(jc_updates, updated_id = if_else(jc_updates$AP_Code == jc_updates$JC_Opinion, 0,1), .after = "JC_Opinion")
jc_updates <- jc_updates %>% rename(AP_list_num = sort)

# Filter out corrected data only from where change was made based on JC second opinion, but notes from all entries

updated_only <- filter(jc_updates, updated_id == 1) %>% select(-JC.notes)
id_notes_only <- select(jc_updates, JC.notes, AP_list_num)
cafi_data_w_updates <- left_join(cafi_data,updated_only, by = "AP_list_num") %>% left_join(id_notes_only, by = "AP_list_num")

# duplicate_check <- group_by(cafi_data_w_updates, AP_list_num) %>% summarize(count = n()) %>% filter(count ==) - One duplicate found but corrected


# Need to: Replace old CAFI code with updated second opinion - Done
#          Add date to revision_date - Done
#          Swap identified by to JC - Done
#          Add JC to revised by - Done
#          Swap ZSP with AP - Done
#          Add notes to ID Notes - Done
#          Add ID to previous_ID - Done

cafi_data_w_updates <- cafi_data_w_updates %>% rename(identification_notes = JC.notes)
cafi_data_w_updates$updated_id <- replace_na(cafi_data_w_updates$updated_id, 0)
cafi_data_w_updates$code <- as.character(cafi_data_w_updates$code)
cafi_data_w_updates$revised_by <- str_replace_all(cafi_data_w_updates$revised_by, "ZSP", "AP")

cafi_data_w_updates$code <- if_else(cafi_data_w_updates$updated_id == 1, 
                                    cafi_data_w_updates$JC_Opinion, 
                                    cafi_data_w_updates$code)

cafi_data_w_updates$revised_by <- if_else(cafi_data_w_updates$updated_id == 1, 
                                    paste(cafi_data_w_updates$revised_by, "JC", sep = "; "), 
                                    cafi_data_w_updates$revised_by)

cafi_data_w_updates$revised_by <- if_else(cafi_data_w_updates$updated_id == 1, 
                                          paste(cafi_data_w_updates$revised_by, "JC", sep = "; "), 
                                          cafi_data_w_updates$revised_by)

cafi_data$code <- as.character(cafi_data$code)
cafi_data_w_updates$previous_name <- as.character(cafi_data_w_updates$previous_name)

cafi_data_w_updates$previous_name <- if_else(cafi_data_w_updates$updated_id == 1, 
                                          paste(cafi_data_w_updates$previous_name, cafi_data$code, sep = "; "), 
                                          cafi_data_w_updates$previous_name)

cafi_data_w_updates$identified_by <- as.character(cafi_data_w_updates$identified_by)

cafi_data_w_updates$identified_by <- if_else(cafi_data_w_updates$updated_id == 1, 
                                             "JC", cafi_data_w_updates$identified_by)

cafi_data_w_updates$revision_date <- as.character(cafi_data_w_updates$revision_date)

cafi_data_w_updates$revision_date <- if_else(cafi_data_w_updates$updated_id == 1, 
                                             paste(cafi_data_w_updates$revision_date, print(Sys.Date()), sep = "; "), 
                                             cafi_data_w_updates$revision_date)

cafi_data_w_updates <- cafi_data_w_updates %>% select(c(names(cafi_data)), identification_notes)

write.csv(cafi_data_w_updates, here("survey/data/cafi_data_summer2019_intermediate_ID_update_2020_5_7.csv"))

