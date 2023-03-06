# Figure Generation


# Read in Libraries
library(dplyr)
library(janitor)
library("DESeq2")
library("ComplexHeatmap")
library("NMF")
library(RNASeqBits)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(viridis)
library(plotly)
library(patchwork)

load("./raw_data_transfigurations.RData")
load("./threshold1_22/normalized_data_configurations.RData")

log2_raw <- log2_raw %>% 
  column_to_rownames("probe_name") %>% 
  select(-accession_number, -annotation)

## PCA Analysis
data.scaled <- na.omit( as.data.frame(t(scale(t(as.matrix(log2_raw))))) )  ##also removes all rows that were assigned a value of 20 for all samples as these will have a NaN value.

pca <-  prcomp(t(data.scaled))

## The percent variance Explained:
# data.frame(summary(pca)$importance)[, 1:min(5, ncol(summary(pca)$importance))]

colorby <- "condition"
a <- ggplot(data = data.frame(pca$x, samples, samples = samples$condition, stringsAsFactors = F), 
       aes(x = as.numeric(PC1), y = as.numeric(PC2), label = name)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca)$importance[2,2] * 100, digits = 2), "% variability" ))+
  scale_color_viridis(discrete = TRUE)

data.scaled <- na.omit( as.data.frame(t(scale(t(as.matrix(log2_norm))))) )  ##also removes all rows that were assigned a value of 20 for all samples as these will have a NaN value.

pca_scaled <-  prcomp(t(data.scaled))

## The percent variance Explained:
# data.frame(summary(pca)$importance)[, 1:min(5, ncol(summary(pca)$importance))]

colorby <- "condition"
b <- ggplot(data = data.frame(pca_scaled$x, samples, samples = samples$condition, stringsAsFactors = F), 
       aes(x = as.numeric(PC1), y = as.numeric(PC2), label = name)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", size=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", size=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca_scaled)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca_scaled)$importance[2,2] * 100, digits = 2), "% variability" )) +
  scale_color_viridis(discrete = TRUE)


c <- a+ b


ggsave(c, 
       filename = "./threshold1_22/pca_comp.pdf",
       device = "pdf",
       dpi = 600,width = 16,height = 8,units = "in")
