## NormFinder Analysis
## 4/26/22
## Amy Olex

## Following instructions from Bei and Jacob Graham
library(tidyverse)
library(janitor)

source('./normfinder/Normfinder_R.R')

file_ext = str_c("threshold1_22",Sys.Date(), ".csv")
directory = ("./threshold1_22/")


### load in data 

dat0 <- read_csv("./nSolver_outputs/background_threshold/threshold_22_1.csv")
dat0 %>% head(10)
dat0 %>% tail(10)

### Unique data cleaning 
dat<- dat0[-2,] %>% 
  clean_names() %>% 
  select(-percent_cv) %>% 
  mutate(probe_name = ifelse(x20220818_210034211023_1_031a_10_rcc == "A",
                             "timepoint",
                             probe_name)) %>% 
  arrange(probe_name == "timepoint") 

dat %>% tail()





### filter data
# According to this tech note, you want to done, sort
# background subtracted data by average counts across all
# samples, and delete all miRNAs expressed below 50 mean
# counts when averaged across all samples.
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


### save the 
write_tsv(datnorm, "normfinder_format.txt",)

## since this is linear count data and not ct-values from qtPCR I set the ctVal to FALSE
results <- Normfinder("normfinder_format.txt", ctVal=FALSE)


results$Ordered

write.csv(results$Ordered, file= str_c(directory, "NormFinder_OrderedResults_", file_ext), row.names = TRUE, sep=",")
write.csv(results$UnOrdered, file=str_c(directory, "NormFinder_UnOrderedResults_", file_ext), row.names = TRUE, quote = FALSE, sep=",")
write.csv(results$PairOfGenes, file= str_c(directory, "NormFinder_PairOfGenesResults_", file_ext), row.names = TRUE, quote = FALSE, sep=",")

##### CV review
# filtered for stability under 0.25 
candidates <- read_csv("./threshold1_22/NormFinder_OrderedResults_threshold1_222022-10-10.csv") %>% 
  filter(Stability < 0.25)

candidates
# 83 candidates 

raw_data <-dat0[-2,] %>% 
  clean_names() %>% 
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
