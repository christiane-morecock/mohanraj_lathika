### WGCNA 
# Weighted gene correlation network analysis

# BiocManager::install("WGCNA")
library(WGCNA)
library(tidyverse)
library(dendextend)
library(heatmaply)
library(tibble)
library(dynamicTreeCut)
# install.packages("gower")
library(gower)


# set working directory
directory <- "./threshold1_22/wgcna/"
load("./raw_data_transfigurations.RData")
load("./threshold1_22/normalized_data_configurations.RData")
ehr_dat <- read_csv("./clean_data/ehr_full.csv")

head(ehr_dat)
colnames(ehr_dat)
threshold = 22


### Clean clinical traits and expression data to proper format
# Samples as rows
# columns as traits
# as.data.frame 

# ALL eHR dataPoints



### Cleaning 
needed_vars  <- ehr_dat %>% 
  filter(timepoint=="C" | timepoint == "D") %>% 
  select(study_id, neutrophil_count, monocyte_count, platelet_count, hemoglobin_level, timepoint, chemotherapy_since_tx, non_chemotherapy_since_tx, relapse) %>% 
  pivot_wider(names_from = "timepoint",
              values_from = c("neutrophil_count",
                              "monocyte_count",
                              "platelet_count",
                              "hemoglobin_level",
                              "chemotherapy_since_tx",
                              "non_chemotherapy_since_tx",
                              "relapse")) 

  clinical_scoring <- ehr_dat %>% 
  filter(timepoint=="A" | timepoint == "B") %>% 
  select(study_id, timepoint, total_moca_score, promis_score, frailty_score, factbmt_score) %>% 
  pivot_wider(names_from = "timepoint",
              values_from = c("total_moca_score",
                              "promis_score",
                              "frailty_score",
                              "factbmt_score"
                              )) 
  
  engraftment_vars <- ehr_dat %>% 
    select(study_id, anc_engraftment_day, plt_engraftment_day) %>%
    distinct() %>% 
    ### Filter STep, removing outlier
    # filter(study_id != 8) %>% 
    #
    column_to_rownames("study_id") %>% 
    as.data.frame()
  
  a_bVars <- needed_vars %>% 
    column_to_rownames("study_id") %>%
    distinct() %>% 
    as.data.frame()
  
  Clin_scoreVars <- clinical_scoring %>% 
    column_to_rownames("study_id") %>% 
    as.data.frame()
  
# Cleaned eHR datasets
  engraftment_vars
  a_bVars
  Clin_scoreVars
 
  
####################################################################   
## Cleaning the expression data
  data <- read_csv("./threshold1_22/threshold1_22_normalized.csv")[-1:-2,] %>% 
    as_tibble()
  
  colnames(data)[-1] <- colnames(data)[-1] %>% str_extract("\\d\\d\\d(A|B)")
  
  
  threshold_na <- function(x){
    na_if(x,threshold)
  }
  
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
  log2_norm <- (log2_norm[rowSums(is.na(log2_norm)) < (ncol(log2_norm)-2), ])


  dim(log2_norm)
  head(log2_norm)
  skimr::skim(log2_norm)
  
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
dim(log2_normA)

# Remove rows that contain  ALL NAs
log2_normA <- (log2_normA[colSums(is.na(log2_normA)) < (nrow(log2_normA)-2)])
dim(log2_normA) 
  

### Time point B #####
log2_normB <- log2_norm %>% 
  select(ends_with("B")) %>% 
  rename_with(.fn = rename_fn, .cols = everything()) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("study_id") %>% 
  mutate(study_id = as.numeric(study_id)) %>% 
  column_to_rownames("study_id")
dim(log2_normB)

# Remove rows that contain ALL NAs
log2_normB <- (log2_normB[colSums(is.na(log2_normB)) < (nrow(log2_normB)-2)])
dim(log2_normB) 

  
# Cleaned expression data:
log2_normA
log2_normB

### save WGCNA data

save(  engraftment_vars,
       a_bVars,
       Clin_scoreVars,
       log2_normA,
       log2_norm,
       log2_normB,
       file = str_c(directory, "clean_expr_clin_dat.RData"))

### Tutorial for the following code:
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.R

# START by gathering the "module eigengenes" (first principal component of a gene cluster) 
# for both A and B timepoints

load(str_c(directory, "clean_expr_clin_dat.RData"))


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(log2_normA, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.50,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



# By the looks of those graphs we pick  for power
power = 9

net = blockwiseModules(log2_normA, power = power,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       verbose = 3)


# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEsA = net$MEs;
geneTree = net$dendrograms[[1]]

#~~~~~~~~~~~~~~~~~~~ TIMEPOINT B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(log2_normB, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.60,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



# By the looks of those graphs we pick  for power
power = 9

net = blockwiseModules(log2_normB, power = power,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       verbose = 3)


# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEsB = net$MEs;
geneTree = net$dendrograms[[1]]




# 1. Engraftment data

# Define numbers of genes and samples
nGenes = ncol(log2_normA);
nSamples = nrow(log2_normA);
# Recalculate MEs with color labels
MEs0A = moduleEigengenes(log2_normA, moduleColors, impute = FALSE)$eigengenes
MEsA = orderMEs(MEs0A)
moduleTraitCor = cor(MEsA, 
                     engraftment_vars, use = "p", method = "spearman");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
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

##### nothing Significant for engraftment data with timepoint A

#~~~~~~~~~~~~~~~~~~~ TIMEPOINT B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# B has 2 less observations
rowsB <- rownames(log2_normB)
engraftment_varsB <- engraftment_vars[rowsB,]


# Define numbers of genes and samples
nGenes = ncol(log2_normB);
nSamples = nrow(log2_normB);
# Recalculate MEs with color labels
MEs0B = moduleEigengenes(log2_normB, moduleColors)$eigengenes
MEsB = orderMEs(MEs0B)

moduleTraitCor = cor(MEsB, 
                     engraftment_varsB, use = "p", method = "spearman");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
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

# nothing significant for B 

# 2. Blood Cell Counts 

# Define numbers of genes and samples
nGenes = ncol(log2_normA);
nSamples = nrow(log2_normA);
# Recalculate MEs with color labels
MEs0A = moduleEigengenes(log2_normA, moduleColors)$eigengenes
MEsA = orderMEs(MEs0A)
#### Changing the correlation? 
moduleTraitCor = cor(MEsA, 
                     a_bVars, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
 par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
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

##### ME Turquoise is significantly correlated to:
 # relapse at timepoint D
 # Therapy since tx at timepoint C
 # Therapy Since tx at timepoint D
 
##### ~ ~ ~ ~ ~ Identifying genes in ME Turquoise ~ ~ ~ ~ ~ ~ ~ ######
 
 # Define variable weight containing the chemotherapy_since_tx_C column of datTrait
 weight = as.data.frame(a_bVars$chemotherapy_since_tx_C)
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

 ### TURQUOISE MOD
 module = "turquoise"
 column = match(module, modNames);
 moduleGenes = moduleColors==module;

 # par(mfrow = c(1,1))
 verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                    abs(geneTraitSignificance[moduleGenes, 1]),
                    xlab = paste("Module Membership in", module, "module"),
                    ylab = "Gene significance for chemotherapy_since_tx_C",
                    main = paste("Module membership vs. gene significance\n"),
                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
 
 names(log2_normA)[moduleColors==module]
 

#~~~~~~~~~~~~~~~~~~~ TIMEPOINT B ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# B has 2 less observations
rowsB <- rownames(log2_normB)
a_bVarsB <- a_bVars[rowsB,]


# Define numbers of genes and samples
nGenes = ncol(log2_normB);
nSamples = nrow(log2_normB);
# Recalculate MEs with color labels
MEs0B = moduleEigengenes(log2_normB, moduleColors)$eigengenes
MEsB = orderMEs(MEs0B)

moduleTraitCor = cor(MEsB, 
                     a_bVarsB, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
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

# chemotherapy_since_tx_C is the only significant finding

##### ~ ~ ~ ~ ~ Identifying genes in ME Turquoise ~ ~ ~ ~ ~ ~ ~ ######

# Define variable weight containing the chemotherapy_since_tx_C column of datTrait
weight = as.data.frame(a_bVarsB$chemotherapy_since_tx_C)
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEsB), 3)

geneModuleMembership = as.data.frame(cor(log2_normB, MEsB, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(log2_normB, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

### TURQUOISE MOD
module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;

# par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for chemotherapy_since_tx_C",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(log2_normB)[moduleColors==module]

