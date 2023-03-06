library(tidyverse)
library(dendextend)
library(heatmaply)
library(tibble)
library(dynamicTreeCut)
library(gower)
library(WGCNA)
suppressPackageStartupMessages(library(gplots))
# dendextend documentation: https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html#introduction

directory <- "./threshold1_22/wgcna/"
load(str_c(directory, "clean_expr_clin_dat.RData"))


#### Corr matrix
# generate correlation matrix of data 
data_cor <- cor(log2_normA, use = 'pairwise.complete.obs')


# if names were numeric
# data_cor[ order(as.numeric(row.names(data_cor))), ]
hc <- data_cor %>%
  # calculate a distance matrix using the euclidean method (standard),
  dist(method = "euclidean") %>%
  # compute hierarchical clustering
  hclust(method = "average")

dend <- hc %>%
  # turn that object into a dendrogram.
  as.dendrogram()

# Let AI pick the clusters
clusters <- cutreeDynamic(hc, deepSplit = TRUE)
# Save number of clusters
n_branch <- length(clusters %>% unique())

str(clusters)
data_clusters<- data_cor %>% as.data.frame()
data_clusters$clusters <- clusters
#######~ ~ ~ ~ ~ ~ ~ ~ ~   W G C N A  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ##############

#### 1. Engraftment

# Define numbers of genes and samples
nGenes = ncol(log2_normA);
nSamples = nrow(log2_normA);
# Recalculate MEs with color labels
MEs0A = moduleEigengenes(log2_normA, clusters, impute = FALSE)$eigengenes
MEsA = orderMEs(MEs0A)
moduleTraitCor = cor(MEsA, 
                     engraftment_vars, 
                     use = "p", 
                     method = "spearman");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
colnames(engraftment_vars) <- c("eng_day_anc", "eng_day_plt")

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(engraftment_vars),
               yLabels = names(MEsA),
               ySymbols = names(MEsA),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Cluster 6 has some significance with ANC engraftment day... Let's remove 
# study_id 8 
# made no difference

engraftment_anc_mirna <- data_clusters %>% 
  filter(clusters ==6 ) %>% 
  rownames()

engraftment_anc_mirna <- log2_normA %>% 
  dplyr::select(all_of(engraftment_anc_mirna))

########################################################################  
# dend <- dend %>%
#   set("labels_colors", k = n_branch) %>%
#   set("branches_k_color", k = n_branch) %>%
#   set("branches_lwd", c(0.5)) %>%
#   ladderize

#################################################################
# library(RColorBrewer)
# 
# col <- colorRampPalette(c("white", "red")) (n=20)
# 
# heatmap.2(data_cor,
#           distfun =function(x) dist(x,method = 'euclidean'),
#           # hclustfun=function(x) hclust(x,method = 'average'),
#           Rowv = dend,
#           Colv = dend,
#           labRow = NA,
#           labCol = NA,
#           density.info="none",
#           dendrogram="both",
#           trace="none",
#           scale="none",
#           # RowSideColors = side_col,
#           # ColSideColors = side_col,
#           col = col)

### 2. Blood Cells Variables ####
# a_bVarsB <- a_bVars[rows_8,]


moduleTraitCor = cor(MEsA, 
                     a_bVars, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(a_bVars),
               yLabels = names(MEsA),
               ySymbols = names(MEsA),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


####
sig_fun <- function(x){
  x < (0.05/14)
}

moduleTraitPvalue %>%
  as.data.frame() %>% 
  mutate_all(sig_fun)

### non_chemotherapy_since_tx_C is significantly correlated with cluster 1,
# non_chemotherapy_since_tx_D is significantly correlated with cluster 1


blood_cellsA <- data_clusters %>% 
  filter(clusters ==1 ) %>% 
  rownames()


blood_cellsA <- log2_normA %>% 
  dplyr::select(all_of(blood_cellsA))

write.csv(blood_cellsA, str_c(directory, "non_chemo_tpA.csv"),row.names = TRUE)

##### ~ ~ ~ ~ ~ Identifying genes in ME 1 ~ ~ ~ ~ ~ ~ ~ ######

# Define variable weight containing the non_chemotherapy_since_tx_Cis  column of datTrait
weight = as.data.frame(a_bVars$non_chemotherapy_since_tx_C )
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEsA), 3)

geneModuleMembership = as.data.frame(cor(log2_normA, MEsA, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(log2_normA, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

### 1 MOD
module = "1"
column = match(module, modNames);
moduleGenes = clusters==module;

# par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for monocyte_count_C ",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


df <- data.frame(x = abs(geneModuleMembership[moduleGenes, column]) ) 
df$y <- abs(geneTraitSignificance[moduleGenes, 1])

df %>% 
  ggplot(aes(x = x, y=y))+
  geom_point() +
  xlab(paste("Module Membership in", module, "module"))+ 
  ylab("Gene significance for monocyte_count_C")+
  ggtitle("Module membership vs. gene significance\n")


### 3. Clinical Scoring  ####

moduleTraitCor = cor(MEsA, 
                     Clin_scoreVars, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(Clin_scoreVars),
               yLabels = names(MEsA),
               ySymbols = names(MEsA),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

####
sig_fun <- function(x){
  x < (0.05/8)
}

moduleTraitPvalue %>%
  as.data.frame() %>% 
  mutate_all(sig_fun)

### frailty_score_A with cluster 7


clinical_scores <- data_clusters %>% 
  filter(clusters ==7 ) %>% 
  rownames()


clinical_scoresA <- log2_normB %>% 
  dplyr::select(all_of(clinical_scores))



#### TIME POINT B ##########################
# data has sd = 0 problems for the miRNAs thresholded to 22, eliminating these counts from the matrix. 
test <- skimr::skim(log2_normB)
dim(log2_normB)


# Filtering out miRNAs that are all the same number (threshold log2(22))
filter_vec <- test %>% 
  filter(numeric.sd != 0) %>% pull(skim_variable)
length(filter_vec)

rows_8 <- rownames(log2_normB)


log2_normB_drop_sd0 <- log2_normB %>% 
  select(all_of(filter_vec)) %>% 
  as_tibble()

# generate correlation matrix of data 
data_cor <- cor(log2_normB_drop_sd0, use = 'pairwise.complete.obs')


# if names were numeric
# data_cor[ order(as.numeric(row.names(data_cor))), ]
hc <- data_cor %>%
  # calculate a distance matrix using the euclidean method (standard),
  dist(method = "euclidean") %>%
  # compute hierarchical clustering
  hclust(method = "average")

dend <- hc %>%
  # turn that object into a dendrogram.
  as.dendrogram()


clusters <- cutreeDynamic(hc, deepSplit = TRUE)
n_branch <- length(clusters %>% unique())
n_branch
data_clusters<- data_cor %>% as.data.frame()
data_clusters$clusters <- clusters
#######~ ~ ~ ~ ~ ~ ~ ~ ~   W G C N A  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ##############

### 1. Engraftment Variables 
engraftment_varsB <- engraftment_vars[rows_8,]

# Define numbers of genes and samples
nGenes = ncol(log2_normB_drop_sd0);
nSamples = nrow(log2_normB_drop_sd0);
# Recalculate MEs with color labels
MEs0B = moduleEigengenes(log2_normB_drop_sd0, clusters, impute = FALSE)$eigengenes
MEsB = orderMEs(MEs0B)
moduleTraitCor = cor(MEsB, 
                     engraftment_varsB, use = "p", method = "spearman");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
colnames(engraftment_varsB) <- c("eng_day_anc", "eng_day_plt")

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(engraftment_varsB),
               yLabels = names(MEsB),
               ySymbols = names(MEsB),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Cluster 8 has some significance with Plt engraftment day
# made no difference

engraftment_plt_mirnaB <- data_clusters %>% 
  filter(clusters ==8 ) %>% 
  rownames()


engraftment_plt_mirnaB <- log2_normB %>% 
  dplyr::select(all_of(engraftment_plt_mirnaB))

### 2. Blood Cells Variables ####
a_bVarsB <- a_bVars[rows_8,]


moduleTraitCor = cor(MEsB, 
                     a_bVarsB, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(a_bVarsB),
               yLabels = names(MEsB),
               ySymbols = names(MEsB),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

####
sig_fun <- function(x){
  x < (0.05/14)
}

moduleTraitPvalue %>%
  as.data.frame() %>% 
  mutate_all(sig_fun)

### Monocyte count C is significantly correlated with cluster 0,
# which is the un=clusterable genes... so that doesn't look good
  

blood_cellsB <- data_clusters %>% 
  filter(clusters ==0 ) %>% 
  rownames()


blood_cellsB <- log2_normB %>% 
  dplyr::select(all_of(blood_cellsB))

##### ~ ~ ~ ~ ~ Identifying genes in ME Turquoise ~ ~ ~ ~ ~ ~ ~ ######

# Define variable weight containing the monocyte_count_C  column of datTrait
weight = as.data.frame(a_bVarsB$monocyte_count_C )
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEsB), 3)

geneModuleMembership = as.data.frame(cor(log2_normB_drop_sd0, MEsB, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(log2_normB_drop_sd0, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

### 0 MOD
module = "0"
column = match(module, modNames);
moduleGenes = clusters==module;

# par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for monocyte_count_C ",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(log2_normA)[moduleColors==module]
df <- data.frame(x = abs(geneModuleMembership[moduleGenes, column]) ) 
df$y <- abs(geneTraitSignificance[moduleGenes, 1])

df %>% 
  ggplot(aes(x = x, y=y))+
  geom_point() +
  xlab(paste("Module Membership in", module, "module"))+ 
  ylab("Gene significance for monocyte_count_C")+
  ggtitle("Module membership vs. gene significance\n")


### 3. Clinical Scoring  ####
Clin_scoreVarsB <- Clin_scoreVars[rows_8,]


moduleTraitCor = cor(MEsB, 
                     Clin_scoreVarsB, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(Clin_scoreVarsB),
               yLabels = names(MEsB),
               ySymbols = names(MEsB),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

####
sig_fun <- function(x){
  x < (0.05/8)
}

moduleTraitPvalue %>%
  as.data.frame() %>% 
  mutate_all(sig_fun)

### frailty_score_A with cluster 8 


clinical_scores <- data_clusters %>% 
  filter(clusters ==8 ) %>% 
  rownames()


clinical_scoresB <- log2_normB %>% 
  dplyr::select(all_of(clinical_scores))

##### ~ ~ ~ ~ ~ Identifying genes in ME 8 ~ ~ ~ ~ ~ ~ ~ ######

# Define variable weight containing the monocyte_count_C  column of datTrait
weight = as.data.frame(Clin_scoreVarsB$frailty_score_A )
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEsB), 3)

geneModuleMembership = as.data.frame(cor(log2_normB_drop_sd0, MEsB, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(log2_normB_drop_sd0, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

### 8 MOD
module = "8"
column = match(module, modNames);
moduleGenes = clusters==module;

# par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for monocyte_count_C ",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


df <- data.frame(x = abs(geneModuleMembership[moduleGenes, column]) ) 
df$y <- abs(geneTraitSignificance[moduleGenes, 1])

df %>% 
  ggplot(aes(x = x, y=y))+
  geom_point() +
  xlab(paste("Module Membership in", module, "module"))+ 
  ylab("Gene significance for monocyte_count_C")+
  ggtitle("Module membership vs. gene significance\n")


##### NEW CODE 02-2023
# Spearman’s rank correlation coefficient, or Spearman’s correlation coefficient, as the name suggests, is a nonparametric approach to measuring correlation using rank values.

log2_norm %>% 
  t() %>% as.data.frame() %>% 
  rownames()


orderMEs(log2_normA)
moduleTraitCor = stats::cor(log2_normA, 
                     engraftment_vars, method = "spearman");

## remove all the NAs
TraitCor_dropNA<- moduleTraitCor[rowSums(is.na(moduleTraitCor)) !=2,]


# Calculate P-values based on sample size
moduleTraitPvalue = corPvalueStudent(TraitCor_dropNA, nSamples)

pvalue_threshold <- 0.2

notSignificant <- moduleTraitPvalue %>%
  as.data.frame() %>% 
  filter(anc_engraftment_day >= pvalue_threshold & plt_engraftment_day>= pvalue_threshold) %>% 
  row.names()

TraitCor_dropNA <- TraitCor_dropNA[!( rownames(TraitCor_dropNA) %in% notSignificant),]
moduleTraitPvalue <- moduleTraitPvalue[!( rownames(moduleTraitPvalue) %in% notSignificant),]
# Will display correlations and their p-values
textMatrix =  paste(signif(TraitCor_dropNA, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(TraitCor_dropNA)
# Display the correlation values within a heatmap plot

library(textshape)

labeledHeatmap(Matrix = cluster_matrix(TraitCor_dropNA),
               xLabels = names(engraftment_vars),
               yLabels = rownames(TraitCor_dropNA),
               ySymbols = rownames(TraitCor_dropNA),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

####
sig_fun <- function(x){
  x <= (0.05/2)
}

moduleTraitPvalue %>%
  as.data.frame() %>% 
  mutate_all(sig_fun) %>% 
  filter(anc_engraftment_day ==TRUE | plt_engraftment_day == TRUE)

engraftment_vars

data("cars")

cars %>% 
ggplot(aes(x=speed, y=dist)) + 
  geom_point(color='#2980B9', size = 4) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color='#2C3E50')
corr<- cor.test(cars$speed, cars$dist, method = "spearman")

cor(x= cars, method = "spearman")
corr

rank(cars$speed)
