## NormFinder Analysis
## 4/26/22
## Amy Olex
### Edited by Christiane Morecock for project

## Following instructions from Bei and Jacob Graham
library(tidyverse)
library(janitor)

source('./normfinder/Normfinder_R.R')

file_ext = str_c("threshold2_28",Sys.Date(), ".csv")
directory = ("./threshold2_28/")


### load in data 

dat0 <- read_csv("./nSolver_outputs/background_threshold/threshold_28_2.csv")
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
write_tsv(datnorm, str_c(directory,"normfinder_format.txt"))

## since this is linear count data and not ct-values from qtPCR I set the ctVal to FALSE
results <- Normfinder(str_c(directory, "normfinder_format.txt"), ctVal=FALSE)


results$Ordered

write.csv(results$Ordered, file= str_c(directory, "NormFinder_OrderedResults_", file_ext), row.names = TRUE)
write.csv(results$UnOrdered, file=str_c(directory, "NormFinder_UnOrderedResults_", file_ext), row.names = TRUE, quote = FALSE)
write.csv(results$PairOfGenes, file= str_c(directory, "NormFinder_PairOfGenesResults_", file_ext), row.names = TRUE, quote = FALSE)

##### CV review
# filtered for stability under 0.25 
candidates <- read_csv(str_c(directory, "NormFinder_OrderedResults_threshold2_282023-02-24.csv")) %>% 
  filter(Stability < 0.25)

candidates
# 86 candidates 

# Include CV in new datatable for probes
raw_data <-dat0[-2,] %>% 
  clean_names() %>% 
  mutate(probe_name = ifelse(x20220818_210034211023_1_031a_10_rcc == "A",
                             "timepoint",
                             probe_name)) %>% 
  arrange(probe_name == "timepoint") %>% 
  column_to_rownames("probe_name")

# Pull just the candidate genes from the datatable
candidates <- candidates %>% 
  pull("...1")
cv_data <- raw_data[candidates, "percent_cv", drop=FALSE]


# Filter from minimum CV% to +20 CV%
firstpass <-cv_data %>% 
  arrange(percent_cv) %>% 
  filter(percent_cv < min(percent_cv) + 20) 

# Join with all information 
housekeepingGenesInfo <- results$Ordered[rownames(firstpass),] %>% rownames_to_column("gene") %>% cbind(firstpass) %>% as_tibble()


write_csv(housekeepingGenesInfo, str_c(directory, "housekeeping_genes.csv"))

### 19 housekeeping genes 
n <- dim(dat)[[2]]

# Check how many probes have under 100 for average counts
dat %>% 
  filter(probe_name %in% housekeepingGenesInfo$gene) %>% 
  mutate_at(2:n, as.numeric) %>% 
  rowwise(probe_name) %>% 
  summarise(average = mean(c_across(everything()))) %>% 
  arrange(average)


housekeepingGenesInfo %>% arrange(gene)
