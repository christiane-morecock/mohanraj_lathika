## NormFinder Analysis
## 4/26/22
## Amy Olex

## Following instructions from Bei and Jacob Graham
library(tidyverse)
library(janitor)

setwd("./mcc_bnfo/miRNA/mohanraj_lathika/normfinder")

source('./normfinder/Normfinder_R.R')

### Clean data 

dat0 <- read_tsv("../normfinder/all_samples_bg-mean1SD_0922.txt")
dat0 %>% head(10)
dat0 %>% tail(10)

dat<- dat0 %>% 
  clean_names() %>% 
  select(-annotation,
         -accession_number,
         -ns_probe_id,
         -class_name,
         -analyte_type,
         -percent_samples_above_threshold,
         -percent_cv) %>% 
  mutate(probe_name = ifelse(x20220818_210034211023_1_031a_10_rcc == "A",
                             "timepoint",
                             probe_name)) %>% 
  arrange(probe_name == "timepoint") 





### filter data
# According to this tech note, you want to done, sort
# background subtracted data by average counts across all
# samples, and delete all miRNAs expressed below 50 mean
# counts when averaged across all samples. (This is definitely too high for this dataset. Unsure which threshold should be used)
# TN_MK1342_Plasma-Serum_r6. Google Docs. Accessed September 22, 2022. https://docs.google.com/document/d/1O-bWRuzHbG00-tw6oRCpt5hWQDujtQfV0IBWaBY-NWA/edit?usp=drive_web&ouid=108696562158233134703

# save the grouping variable 
char_row <- dat %>% 
  subset(probe_name == "timepoint")

# convert datafram to numeric to calculate the rowmeans 
num_dat <- dat %>% 
  subset(probe_name != "timepoint") %>%
  mutate_at(2:47, as.numeric)

dat1 <- num_dat %>% 
  # calculate row means
  rowwise() %>% 
  mutate(total = mean(c_across(2:47))) %>% 
  # filter rows where mean is less than 50
  filter(total >= 50 ) %>% 
  # remove the total column
  select(-total) %>% 
  # convert columns back to character
  mutate_all(as.character)



datnorm <- dat1 %>% 
  rows_append(char_row) 



write_tsv(datnorm, "housekeeping_all_samples_bg-mean1SD_0922_nf_format.txt",)

## since this is linear count data and not ct-values from qtPCR I set the ctVal to FALSE
results <- Normfinder("housekeeping_all_samples_bg-mean1SD_0922_nf_format.txt", ctVal=FALSE)



file_ext = str_c("bg-mean1SD_",Sys.Date(), ".csv")

write.table(results$Ordered, file= str_c("NormFinder_OrderedResults_", file_ext), row.names = TRUE, quote = FALSE, sep=",")
write.table(results$UnOrdered, file=str_c("NormFinder_UnOrderedResults_", file_ext), row.names = TRUE, quote = FALSE, sep=",")
write.table(results$PairOfGenes, file= str_c("NormFinder_PairOfGenesResults", file_ext), row.names = TRUE, quote = FALSE, sep=",")

#####CV review

candidates <- read_csv("./normfinder/NormFinder_OrderedResults_bg-mean1SD_2022-09-27.csv") %>% 
  filter(Stability < 0.25)

raw_data <-dat0 %>% 
  clean_names() %>% 
  select(-annotation,
         -accession_number,
         -ns_probe_id,
         -class_name,
         -analyte_type,
         -percent_samples_above_threshold) %>% 
  mutate(probe_name = ifelse(x20220818_210034211023_1_031a_10_rcc == "A",
                             "timepoint",
                             probe_name)) %>% 
  arrange(probe_name == "timepoint") %>% 
  column_to_rownames("probe_name")


candidates <- candidates %>% 
  pull("...1")

cv_data <- raw_data[candidates, "percent_cv", drop=FALSE]

firstpass <-cv_data %>% 
  arrange(percent_cv) %>% 
  filter(percent_cv < min(percent_cv) + 20) 

write.csv(firstpass, "housekeeping_genes.csv")


results$Ordered[firstpass,] %>% rownames_to_column("gene") %>% 
  arrange(gene)
