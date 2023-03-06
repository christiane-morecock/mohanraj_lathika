# original script: wgcna/correlation_tests.R
# Author: Christiane Morecock
# Purpose: clean clinical data into data subsets
#
#
#
#
#
#
#
#
#
#
########################################

# Load libraries
library(tidyverse)
library(janitor)
library(skimr)
library(gghighlight)



## Define your directory
directory <- "./threshold1_18_removeNegC/"
## Define your file name
filename <- "threshold1_18_removeNegC" 
## Define your threshold
threshold <- 18

### load in data 
data <- read_csv(str_c(directory, filename, "_normalized.csv"))[-1:-2,] %>% 
  as_tibble()
### Explore the data

head(data)
skim(data)


### Unique cleaning of sample names #########
# Run only one time 
colnames(data)[-1] <- colnames(data)[-1] %>% 
str_extract("\\d\\d\\d(A|B)")
# Look at your data
head(data)

# Function to turn threshold into NA
threshold_na <- function(x){
  na_if(x,threshold)
}

#  log2(data+1)
plus1_log2 <- function(x){
  log2(x+1)
}

#log2 of all the data+1 to eliminate all 0 
log2_norm <- data %>%
  mutate_at(2:47, as.numeric) %>% 
  # convert threshold to NA
  mutate_if(is.numeric, threshold_na) %>% 
  mutate_if(is.numeric,plus1_log2) %>% 
  column_to_rownames("Probe Name") %>% 
  as.data.frame()

# Remove rows that contain ALL NAs
log2_norm <- (log2_norm[rowSums(is.na(log2_norm)) != ncol(log2_norm), ])



# A timepoint and B timepoints 

rename_fn <- function(x){
  str_remove(x,"(A|B)")
}
log2_normA <-  
  log2_norm %>% 
  select(ends_with("A")) %>% 
  rename_with(.fn = rename_fn, .cols = everything()) %>%
  t() %>%  
  as.data.frame() %>% 
  rownames_to_column("study_id") %>% 
  mutate(study_id = as.numeric(study_id)) %>% 
  column_to_rownames("study_id") 



log2_normB <- log2_norm %>% 
  select(ends_with("B")) %>% 
  rename_with(.fn = rename_fn, .cols = everything()) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("study_id") %>% 
  mutate(study_id = as.numeric(study_id)) %>% 
  column_to_rownames("study_id")

# Cleaned expression data:
log2_normA
log2_normB
log2_norm