library(tidyverse)
# library(dendextend)
library(heatmaply)
library(textshape)
library(gplots)
library(WGCNA)
library(MASS)
library(knitr)

directory <- "./threshold1_22/wgcna/"
load(str_c(directory, "clean_expr_clin_dat.RData"))


##### NEW CODE 02-2023
# Define numbers of genes and samples
nGenes = ncol(log2_normA);
nSamples = nrow(log2_normA);

# Spearman’s rank correlation coefficient, or Spearman’s correlation coefficient, as the name suggests, is a nonparametric approach to measuring correlation using rank values.
moduleTraitCor = stats::cor(log2_normA, 
                            engraftment_vars, 
                            method = "spearman")# spearman: rank correlation 

## remove all the NAs
TraitCor_dropNA<- moduleTraitCor[rowSums(is.na(moduleTraitCor)) !=2,]
TraitCor_dropNA_reorder <- cluster_matrix(TraitCor_dropNA)

# Calculate P-values based on sample size (package: WGCNA)
TraitPvalue = corPvalueStudent(TraitCor_dropNA_reorder, nSamples)

#### trim the insignificant genes to make the figure more readable 
pvalue_threshold <- 0.15
# identify non significant genes according to threshold
notSignificant <- TraitPvalue %>%
  as.data.frame() %>% 
  filter(anc_engraftment_day >= pvalue_threshold & plt_engraftment_day >= pvalue_threshold) %>% 
  row.names()
# remove genes from data sets
TraitCorFig <- TraitCor_dropNA_reorder[!( rownames(TraitCor_dropNA_reorder) %in% notSignificant),]
TraitPvalueFig <- TraitPvalue[!( rownames(TraitPvalue) %in% notSignificant),]

# TextMatrix will display correlations and their p-values
textMatrix =  paste(signif(TraitCorFig, 2), "\n(",
                    signif(TraitPvalueFig, 1), ")", sep = "");
dim(textMatrix) = dim(TraitCorFig)


# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = TraitCorFig,
               xLabels = names(engraftment_vars),
               yLabels = rownames(TraitCorFig),
               ySymbols = rownames(TraitCorFig),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#### Useless function.. Just go with it. 
## Logical, identifying significance at valuewritten in function.
sig_fun <- function(x){
  x <= (0.05/2)
}

TraitPvalueFig %>%
  as.data.frame() %>% 
  mutate_all(sig_fun) %>% 
  filter(anc_engraftment_day ==TRUE | plt_engraftment_day == TRUE) %>% 
  kable()


##### Step AIC #### 
# THIS IS A DEAD END. This is not working. 

engraftmentVars_tib <-engraftment_vars %>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  as_tibble()


Clin_scoreVars_tib <-Clin_scoreVars %>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  as_tibble()


bloodCellVars_tib <-a_bVars %>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  as_tibble()

datEngraftmentA <- log2_normA %>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  as_tibble() %>% 
  left_join(engraftmentVars_tib) %>% 
  mutate_if(is.numeric, rank)%>% 
  column_to_rownames("ID") %>% 
  dplyr::select(-starts_with("POS"), -starts_with("NEG"), -starts_with("LIG"))

datEngraftmentB <- log2_normB %>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  as_tibble() %>% 
  left_join(engraftmentVars_tib) %>% 
  mutate_if(is.numeric, rank)%>% 
  column_to_rownames("ID") %>% 
  dplyr::select(-starts_with("POS"), -starts_with("NEG"), -starts_with("LIG"))

### full equation
full_model <- lm(anc_engraftment_day ~ . - plt_engraftment_day,
                  data = datEngraftmentA, na.action = na.omit) # Run regression with Poisson distribution
# Empty equation
null_model <- lm(anc_engraftment_day  ~ 1, data = datEngraftmentA)

model_sel <- stepAIC(null_model,
                     scope=list(lower = null_model, upper = full_model),
                     direction = "forward",
                     trace = 0)
model_sel$anova

lm.fit <- lm(formula = anc_engraftment_day ~ `hsa-miR-301b-5p` + `hsa-miR-3168` + 
               `hsa-miR-660-5p` + `hsa-miR-185-5p` + `hsa-miR-6720-3p` + 
               `hsa-miR-2110` + `hsa-miR-4443` + `hsa-miR-1203` + `hsa-miR-1286` + 
               `hsa-miR-497-5p` + `hsa-miR-378f` + `hsa-miR-196a-5p` + `hsa-miR-188-3p` + 
               `hsa-miR-29c-3p` + `hsa-miR-5001-5p` + `hsa-miR-3144-3p` + 
               `hsa-miR-197-3p` + `hsa-miR-656-3p` + `hsa-miR-126-3p` + 
               `hsa-miR-203a-5p` + `hsa-miR-1537-3p` + `hsa-miR-186-5p`, data = datEngraftmentA)

kable(summary(lm.fit)$coefficients, digits = 6)
summary(lm.fit)$adj.r.squared
library(car)
vif(lm.fit)
## Variables removed for multicolinearity:
# `hsa-miR-660-5p`
# `hsa-miR-6720-3p`
# `hsa-miR-2110`
# `hsa-miR-1203`
# `hsa-miR-1537-3p`

lm.fit <- lm(formula = anc_engraftment_day ~ `hsa-miR-301b-5p` + `hsa-miR-3168` +
               `hsa-miR-185-5p` + `hsa-miR-4443` + `hsa-miR-1286` + 
               `hsa-miR-497-5p` + `hsa-miR-378f` + `hsa-miR-196a-5p` + `hsa-miR-188-3p` + 
               `hsa-miR-29c-3p` + `hsa-miR-5001-5p` + `hsa-miR-3144-3p` + 
               `hsa-miR-197-3p` + `hsa-miR-656-3p` + `hsa-miR-126-3p` + 
               `hsa-miR-203a-5p` + `hsa-miR-186-5p`, data = datEngraftmentA)

kable(summary(lm.fit)$coefficients, digits = 6)

### miRNA with pvalue under 0.05
# `hsa-miR-301b-5p`
# `hsa-miR-185-5p`



### full equation
full_model <- lm(plt_engraftment_day ~ . - anc_engraftment_day,
                  data = datEngraftmentA, na.action = na.omit) # Run regression with Poisson distribution
# Empty equation
null_model <- lm(anc_engraftment_day  ~ 1, data = datEngraftmentA)

model_sel <- stepAIC(null_model,
                     scope=list(lower = null_model, upper = full_model),
                     direction = "forward",
                     trace = 0)
model_sel$anova

lm.fit <- lm(formula = anc_engraftment_day ~ `hsa-miR-301b-5p` + `hsa-miR-3168` + 
                `hsa-miR-660-5p` + `hsa-miR-185-5p` + `hsa-miR-6720-3p` + 
                `hsa-miR-2110` + `hsa-miR-4443` + `hsa-miR-1203` + `hsa-miR-1286` + 
                `hsa-miR-497-5p` + `hsa-miR-378f` + `hsa-miR-196a-5p` + `hsa-miR-188-3p` + 
                `hsa-miR-29c-3p` + `hsa-miR-5001-5p` + `hsa-miR-3144-3p` + 
                `hsa-miR-197-3p` + `hsa-miR-656-3p` + `hsa-miR-126-3p` + 
                `hsa-miR-203a-5p` + `hsa-miR-1537-3p` + `hsa-miR-186-5p` , data = datEngraftmentA)


vif(full_model)
## Variables removed for multicolinearity:
# `hsa-miR-1287-3p`
# `hsa-miR-1203`
# `hsa-miR-6720-3p`
# `hsa-miR-126-3p`
# `hsa-miR-1537-3p`
# `hsa-miR-5001-5p`

lm.fit <- lm(formula = anc_engraftment_day ~ `hsa-miR-301b-5p` + `hsa-miR-3168` + 
                `hsa-miR-660-5p` + `hsa-miR-185-5p` + `hsa-miR-2110` +`hsa-miR-4443` +
                `hsa-miR-1286` + `hsa-miR-497-5p` + `hsa-miR-378f` +
                `hsa-miR-196a-5p` + `hsa-miR-188-3p` + `hsa-miR-29c-3p` + 
                `hsa-miR-3144-3p` + `hsa-miR-197-3p` + `hsa-miR-656-3p` + 
                `hsa-miR-203a-5p` + `hsa-miR-186-5p` , data = datEngraftmentA)

kable(summary(lm.fit)$coefficients, digits = 6)
summary(lm.fit)$adj.r.squared

# significant indicators:
# `hsa-miR-301b-5p` 
# `hsa-miR-660-5p` 
# `hsa-miR-185-5p`  
# `hsa-miR-4443`
# `hsa-miR-2110` 

full_model <- lm(formula = anc_engraftment_day ~ `hsa-miR-301b-5p` + `hsa-miR-3168` + 
     `hsa-miR-660-5p` + `hsa-miR-185-5p` + `hsa-miR-2110` +`hsa-miR-4443` +
     `hsa-miR-1286` + `hsa-miR-497-5p` + `hsa-miR-378f` +
     `hsa-miR-196a-5p` + `hsa-miR-188-3p` + `hsa-miR-29c-3p` + 
     `hsa-miR-3144-3p` + `hsa-miR-197-3p` + `hsa-miR-656-3p` + 
     `hsa-miR-203a-5p` + `hsa-miR-186-5p` , data = datEngraftmentA)

null_model <- lm(anc_engraftment_day  ~ 1, data = datEngraftmentA)

model_sel <- stepAIC(full_model,
                     scope=list(lower = null_model, upper = full_model),
                     direction = "backward",
                     trace = 0)
model_sel$anova


lm.fit <- lm(formula = anc_engraftment_day ~ `hsa-miR-301b-5p` + `hsa-miR-3168` + 
     `hsa-miR-660-5p` + `hsa-miR-185-5p` + `hsa-miR-2110` + `hsa-miR-4443` + 
     `hsa-miR-378f` + `hsa-miR-196a-5p` + `hsa-miR-188-3p` + `hsa-miR-29c-3p` + 
     `hsa-miR-3144-3p` + `hsa-miR-656-3p` + `hsa-miR-203a-5p` + 
     `hsa-miR-186-5p`, data = datEngraftmentA)

summary(lm.fit)
