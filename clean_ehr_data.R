library(tidyverse)
library(janitor)
library(readxl)
library(skimr)
library(readr)

demo <- read_excel("./raw_data/CombinedDataSet-LM.xlsx", sheet = "Demographics")
fragility <- read_excel("./raw_data/CombinedDataSet-LM.xlsx", sheet = "FRAP")
tb <- read_excel("./raw_data/CombinedDataSet-LM.xlsx", sheet = "TB")
tb2 <- read_excel("./raw_data/TinaBachas_Data-LM Sept 22 FINAL send.xlsx", sheet = "Forms A&B-copyv3 FINAL")



# clean demo data
demo <-demo %>% clean_names() %>% 
  separate(household_income,
           into = c("income_min","income_max"),
           sep = (" to | than | or ")) %>% 
  mutate(income_min = ifelse(income_min == "Less",
                             0,
                             income_min),
         income_max = ifelse(income_max == "more",
                             999999999999999,
                             income_max)) %>% 
  mutate_at(7:8, parse_number)

# clean fragility dat 

fragility <- fragility %>% 
  clean_names() %>% 
  separate(identifier,
           into = c("study_id", "timepoint")) %>% 
  mutate(study_id = as.numeric(study_id))


# clean tb dat
tb <- tb[,1:12] %>% 
  clean_names() %>% 
  filter(!is.na(as.numeric(record_id)))

# Tb2 has a column about second transplant that is absent from tb,
# cleaning to join that column appropriately
tb2 <- tb2[,1:11] %>% 
  clean_names() %>% 
  filter(!is.na(as.numeric(record_id)))

tb <- tb %>% 
  left_join(tb2) %>% 
  rename("study_id" = "record_id",
         "months_after_tx" = "event_name") %>% 
  # create numeric event data
  mutate(months_after_tx = ifelse(months_after_tx == "6 Months",
                                  6,
                                  ifelse(months_after_tx == "1 Year",
                                         12,
                                         NA)),
         # Add a timepoint to separate data from other timepoints
         timepoint= ifelse(months_after_tx == 6,
                           "C",
                           ifelse(months_after_tx == 12,
                                  "D",
                                  NA)),
         #study_id should be numeric
         study_id = as.numeric(study_id),
         # recode "yes" and "no" as 1s and 0s, then factor
         relapse = ifelse(relapse == "Yes",
                          1,
                          ifelse(relapse == "No",
                                 0,
                                 NA)),
         chemotherapy_since_tx = ifelse(chemotherapy_since_tx == "Yes",
                          1,
                          ifelse(chemotherapy_since_tx == "No",
                                 0,
                                 NA)),
         non_chemotherapy_since_tx = ifelse(non_chemotherapy_since_tx == "Yes",
                                        1,
                                        ifelse(non_chemotherapy_since_tx == "No",
                                               0,
                                               NA))) %>% 
  mutate_at(9:11, as.factor)
  
  
# Join statement
first <- demo %>% 
  full_join(fragility) 

second <- demo %>% 
  full_join(tb) 

ehr_joined <- first %>% 
  full_join(second) %>%
  relocate(timepoint, .after = study_id)


write_csv(ehr_joined, "./clean_data/ehr_full.csv")

  


