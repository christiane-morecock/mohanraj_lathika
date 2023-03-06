
# Load libraries
library(tidyverse)
library(janitor)
library(skimr)
library(gghighlight)

## Define your directory
directory <- "./raw_data/"
## Define your file name
filename <- "all_samples" 
## Define your threshold
# No threshold for Raw data

### load in data 
dat <- read_csv(str_c(directory, filename,"_raw.csv")

### Explore the data
head(dat)
data <- dat[-1,] %>% clean_names()
head(data)
skim(data)
raw_data <- data


### Unique cleaning of sample names #########
# Run only one time 
colnames(data)[-1:-3] <- colnames(data)[-1:-3] %>%
   str_extract("\\d\\d\\d(a|b)")

 
colnames(data)[4:49] <- colnames(data)[4:49] %>% 
  str_extract("(a|b)") %>% 
  str_c(colnames(data)[4:49])
# Look at your data
head(data)

# Create Metadata #
df <- data.frame( colnames(dat)[-1:-3])

colnames(df) <- c("file")

str(df)

df <- df %>% 
  mutate(sample = str_extract(file, 
                              "\\d\\d\\d(A|B)"),
         timepoint = str_extract(file, "(A|B)"),
         cart = str_extract(file, "21003(\\d{3})\\d+"))


############# transpose the data  ############
dim(data)
# transpose retaining the probenames 
intermediate <- data %>% 
  select(-annotation,
         -accession_number) %>% 
  column_to_rownames("probe_name") %>% 
  t() 

# recapture the sample names
names <- rownames(intermediate)


# create a transposed tibble  that can use tidyverse functions 
trans_dat <- intermediate %>% 
  as_tibble() %>% 
  clean_names()

trans_dat$names <- names

trans_dat

############ Plot out the spike-ins ##############


# Graph the spike-ins over each sample 
trans_dat %>% 
  mutate(max_spike = osa_mi_r414 == max(osa_mi_r414),
         min_spike = osa_mi_r414 == min(osa_mi_r414)) %>% 
  ggplot(aes(y= osa_mi_r414, x = fct_reorder(names, osa_mi_r414)))+
  geom_point()+
  gghighlight::gghighlight(max_spike ==TRUE |  min_spike == TRUE,
                           label_key = osa_mi_r414)+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Plotting the Spike-in Control across samples")



#### Calculating thresholds ##############

dim(trans_dat)
class(trans_dat)
str(trans_dat)


# gather means of all negative controls 
means_vec <- trans_dat %>% 
  select(starts_with("neg", ignore.case = TRUE)) %>% 
  colMeans() 

means_vec_rmNegC <- means_vec[-3]


# mean of mean + 1sd 
threshold1 <-  round(mean(means_vec) + sd(means_vec), digits = 1)
threshold1
# 22
threshold1_nc <-  round(mean(means_vec_rmNegC) + sd(means_vec_rmNegC), digits = 1)
threshold1_nc
# 18
# mean of mean + 2sd 
threshold2 <-  round(mean(means_vec) + 2 * sd(means_vec), digits = 1)
threshold2
# 28
threshold2_nc <-  round(mean(means_vec_rmNegC) + 2 * sd(means_vec_rmNegC), digits = 1)
threshold2_nc
# 22

# mean of mean + 3sd 
threshold3 <-  round(mean(means_vec) + 3 * sd(means_vec), digits = 1)
threshold3
# 33
threshold3_nc <-  round(mean(means_vec_rmNegC) + 3 * sd(means_vec_rmNegC), digits = 1)
threshold3_nc
# 26

# mean of mean + 4sd 
threshold4 <-  round(mean(means_vec) + 4 * sd(means_vec), digits = 1)
threshold4
# 39
threshold4_nc <-  round(mean(means_vec_rmNegC) + 4 * sd(means_vec_rmNegC), digits = 1)
threshold4_nc
# 30




### log2 normalize data

log2_raw <- data %>% 
  mutate_at(4:49, log2)

############# transpose the log2 data  ############
dim(log2_raw)
# transpose retaining the probenames 
intermediate <- log2_raw %>% 
  select(-annotation,
         -accession_number) %>% 
  column_to_rownames("probe_name") %>% 
  t() 

# recapture the sample names
names <- rownames(intermediate)


# create a transposed tibble  that can use tidyverse functions 
log2_trans_dat <- intermediate %>% 
  as_tibble() %>% 
  clean_names()

log2_trans_dat$names <- names

log2_trans_dat <- log2_trans_dat %>% 
  relocate(names)
log2_trans_dat

###### Plotting all the data #######

raw_dat_long <- data %>% 
  pivot_longer(cols =  4:49,
               names_to = "sample",
               values_to = "raw_counts") %>% 
  mutate(group = str_extract(sample, "(a|b)"))

log_dat_long <- log2_raw %>% 
  pivot_longer(cols =  4:49,
               names_to = "sample",
               values_to = "raw_counts") %>% 
  mutate(group = str_extract(sample, "(a|b)"))


raw_dat_long %>% 
  # filter(raw_counts < 200) %>% 
  ggplot(aes(x= probe_name, y = raw_counts,
             color = group))+
  geom_jitter(size = 1,
              alpha = 0.5) +
  scale_colour_grey()+
  geom_hline(yintercept = log2(threshold1),
             color ="#FF5733",
             size = 0.5) +
  geom_hline(yintercept = log2(threshold2),
             color ="#FF7833",
             size = 0.5)+
  geom_hline(yintercept = log2(threshold3),
             color ="#FF9F33",
             size = 0.5) +
  geom_hline(yintercept = log2(threshold4),
             color ="#FFCF33",
             size = 0.5)+
  theme(axis.text.x = element_blank())

#### Write datasets for use in the report ### 
save(log2_raw, raw_data, means_vec, means_vec_rmNegC, df, trans_dat, raw_dat_long, log_dat_long, file = "raw_data_transfigurations.RData")



