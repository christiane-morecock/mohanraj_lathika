## DEG Analysis on normalized data
## Amy Olex
## 4/26/22
## This script was written by Amy Olex, but modified by Christiane Morecock for Lathika Mohanraj data. 
### this script uses normalized data by the threshold and generates: 

# install.packages("remotes")
# remotes::install_github("AmyOlex/RNASeqBits")

library("DESeq2")
library("ComplexHeatmap")
library("NMF")
library(RNASeqBits)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(corrplot)
library(gplots)
library(dendextend)
library(heatmaply)
suppressPackageStartupMessages(library(gplots))
# dendextend documentation: https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html#introduction


directory <- "threshold2_28/"
filename <- "threshold2_28"
threshold <- 28
#data with removed timepoint
data <- read_csv(str_c(directory, filename, "_normalized.csv"))[-1:-2,] %>% 
  as_tibble()
#same data, retained timepoint
data_norm <- read_csv(str_c(directory, filename, "_normalized.csv"))[-2,]
# Capturing sample name
colnames(data_norm)[-1] <- colnames(data_norm)[-1] %>% str_extract("\\d\\d\\d(A|B)")


############# transpose the data  ############
# count the number of rows 
n <- dim(data_norm)[[1]]
# transpose retaining the probenames 
intermediate2 <- data_norm %>% 
  mutate( `Probe Name` = ifelse(is.na(`Probe Name`),
                                "Timepoint",
                                `Probe Name`)) %>% 
  column_to_rownames("Probe Name") %>% 
  t() 

# recapture the sample names
names2 <- rownames(intermediate2)

# create a transposed tibble  that can use tidyverse functions 
trans_norm <- intermediate2 %>% 
  as_tibble() %>% 
  mutate_at(2:n, as.numeric)

# names as the first column
trans_norm$names<-names2

trans_norm<- trans_norm %>% 
  relocate(names)
##################################


n2 <- dim(data)[[2]]
## Form matrix with all numeric values
data <- data %>% 
  column_to_rownames("Probe Name")%>%
  mutate_at(1:(n2-1), as.numeric)



# create metadata for each sample 
metadata <- data %>% 
  colnames() %>% 
  as.data.frame() 

colnames(metadata) <- "sample"
samples <- metadata %>% 
  # extract time point
  mutate(condition = str_extract(sample, "(A|B)"),
         # extract sample name
         name = str_extract(sample, "[0-9]{3}(A|B)"),
         # extract cartridge ID
         cart = str_extract(sample, "21003(\\d{3})\\d+"))

# write_csv(samples, "./raw_sample_metadata.csv")

########################################################################
# Continue Script

# Capture sample name from file name, rename columns 
colnames(data) <- colnames(data) %>% str_extract("[0-9]{3}(A|B)") 

# I'm not sure if this is necessary yet 
# row names of metadata transformed to sample names. 
all(row.names(samples) == names(data))
row.names(samples) <- colnames(data) 

#log2 of all the data+1 to eliminate all 0 
log2_norm <- log2(data+1)

###############################################################################

###### Correlation plot ######

###############################################################################

# generate correlation matrix of data 
data_cor <- cor(log2_norm)


# if names were numeric use the following line of code:
# data_cor[ order(as.numeric(row.names(data_cor))), ]

##### Create an Hclust object ####### 
hc <- data_cor %>%
# calculate a distance matrix using the euclidean method (standard),
   dist(method = "euclidean") %>%
# compute hierarchical clustering
   hclust(method = "average")

###### Create a dendrogram object from the hclust object ######
  dend <- hc %>%
    # turn that object into a dendrogram.
     as.dendrogram()

##### Color the branches, and change the width ####
dend <- dend %>%
  # set("labels_colors", cutree(dend, k = 6)) %>%
   set("branches_k_color", k = 6) %>%
   set("branches_lwd", c(0.5)) %>%
   ladderize


##### Heatmap correlation plot ######
col <- colorRampPalette(c("white", "red")) (n=20)
groups <- data_cor %>% colnames %>%  str_detect("A")

side_col <- ifelse(groups == TRUE, "yellow", "purple")
  
heatmap.2(data_cor,
          distfun =function(x) dist(x,method = 'euclidean'),
         # hclustfun=function(x) hclust(x,method = 'average'),
         Rowv = dend,
         Colv = dend,
         density.info="none",
         dendrogram="both",
         trace="none",
         scale="none",
         RowSideColors = side_col,
         ColSideColors = side_col,
         col = col)




# FDR cutoff.
MIN_FDR = 0.05

# Set the margins
MARGINS = c(10, 10)

# Relative heights of the rows in the plot.
LHEI = c(1, 5)

# log2 norm as a matrix
values <- log2_norm %>% 
  as.matrix() 

# shouldn't need jitter, it eliminates 0s without 
values = jitter(values, factor = 1, amount = 0.00001)

zscores = NULL
for (i in 1 : nrow(values)) {
  row = values[i,]
  zrow = (row - mean(row, na.rm = TRUE)) / sd(row, na.rm = TRUE)
  zscores = rbind(zscores, zrow)
}
 

# Set the row names on the zscores.
row.names(zscores) = row.names(values)

# Turn the data into a matrix for heatmap2.
zscores = as.matrix(zscores) 


# Set the color palette.
col = greenred

# create color string
test <- zscores %>% colnames %>% str_detect("A")

# Generate heatmap
heatmap_save <- heatmap.2(zscores,
                          na.rm=TRUE,
                          col=col,
                          density.info="none",
                          dendrogram="both",
                          trace="none",
                          ColSideColors=ifelse(test == TRUE, "yellow", "purple"),
                          labRow=FALSE,
                          margins=MARGINS,
                          lhei=LHEI
)




## PCA Analysis
data.scaled <- na.omit( as.data.frame(t(scale(t(as.matrix(log2_norm))))) )  ##also removes all rows that were assigned a value of 20 for all samples as these will have a NaN value.

pca <-  prcomp(t(data.scaled))

## The percent variance Explained:
# data.frame(summary(pca)$importance)[, 1:min(5, ncol(summary(pca)$importance))]

colorby <- "condition"
ggplot(data = data.frame(pca$x, samples, samples = samples$condition, stringsAsFactors = F), 
       aes(x = as.numeric(PC1), y = as.numeric(PC2), label = name)) +
  theme(plot.title = element_text(lineheight = 0.8, face="bold")) +
  ggtitle(paste("PCA Analysis, coloring by ", colorby)) +
  geom_point(aes(color = eval(parse(text = colorby))), size = 3) +
  geom_text_repel(colour = "black", size = 3) +
  geom_hline(yintercept = 0, colour = "lightgrey", linewidth=.1) +
  geom_vline(xintercept = 0, colour = "lightgrey", linewidth=.1) +
  labs(color = colorby) +
  theme_bw() +
  scale_x_continuous(name = paste0("PC1, ", round(summary(pca)$importance[2,1] * 100, digits = 2), "% variability" )) +
  scale_y_continuous(name = paste0("PC2, ", round(summary(pca)$importance[2,2] * 100, digits = 2), "% variability" ))


 ### TESTING HEATMAPS?
## Expression heatmap
centered_colors <- center.palette(data.scaled, palette_length = 100)

#gene_vars = apply(data.scaled, 1, var)

#top2000 = data.scaled[names(sort(gene_vars, decreasing = TRUE))[1:2000],,drop=FALSE]

#gene_sym_top2000 <- human_TPM_genes[row.names(top2000),,drop=FALSE]
#all(row.names(gene_sym_top2000)==row.names(top2000))
#row.names(top2000) <- gene_sym_top2000$symbol

#centered_colors <- center.palette(top2000, palette_length = 100)
#my_annot <- HeatmapAnnotation(df = annot_data_50m, col = list(Tissue=c("Brain"="red", "Liver" = "blue", "Lung" = "black") ) )

# png(filename=paste("Heatmap_log2_NormalizedCount_prelim.png", sep=""), width=2000, height=3000, res=150)
# 
# Heatmap(as.matrix(data.scaled), col = centered_colors$colors,  
#         column_dend_height = unit(3, "cm"), show_row_dend = TRUE, row_dend_width = unit(3, "cm"),
#         clustering_method_columns="ward.D2", clustering_method_rows="ward.D2",
#         show_row_names=FALSE, show_column_names=TRUE, column_names_gp = gpar(fontsize = 30),  heatmap_legend_param = list(title = "Scaled Log2(NormCount)", 
#                                                                                                                           title_gp = gpar(fontsize=20), 
#                                                                                                                          labels_gp = gpar(fontsize = 20),
                                                                                                                          # legend_height = unit(4, "in")) ) 

# dev.off() 




### Box Plot of Log2 Normalized Data
library(reshape2)

mtx_to_plot <- melt(log2_norm)
colnames(mtx_to_plot) <- c("Sample", "Exp")
#mtx_to_plot <- left_join(mtx_to_plot, annot)

# png(filename=paste("BoxPlot_log2_NormalizedCount_prelim.png", sep=""), width=800, height=800, res=300)

mtx_to_plot %>% 
  ggplot(aes(x = Sample, y = Exp))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("log2-Expression")

# dev.off()

####################################################################################################

### DESeq2 Analysis for timepoint A vs. B
## extract selected samples
A_B <- samples[samples$condition %in% c("A", "B"),]
A_B$condition <- factor(A_B$condition)

## Extract the data columns that match these sample names
selected_data <- data[,row.names(A_B)]

## Sanity Check, should be TRUE
stopifnot(all(names(selected_data) == row.names(A_B)))

## import into DESeq2
dds <- DESeqDataSetFromMatrix(countData=round(selected_data, digits=0), colData = A_B, design = ~ condition)

############# Contrast timepoint A (control) vs timepoint B

contrast <- "Avs.B"
control_cond="A"
experiment_cond="B"
date = "final"

### Only need to do the following as long as MGT is used as the control ###
## re-level dds so that the control is the reference
dds$condition <- relevel(dds$condition, ref=control_cond)

#Perform differential expression analysis
dds1 <- DESeq(dds)

res <- results(dds1, contrast = c("condition", experiment_cond, control_cond))


# jpeg(filename = paste("MAPlot", contrast, date, ".jpg", sep = "_"), width = 800, height = 600, units = "px", res = 100)
  plotMA(res, main = paste("DESeq2", contrast, date, sep = " "))
# dev.off()

degs <- na.omit(as.data.frame(res))
write.table(degs, file = paste(directory, "DESeq2_DEGs_", contrast, "_", date, ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.csv(degs, file = paste(directory, "DESeq2_DEGs_", contrast, "_", date, ".csv", sep = ""), row.names=TRUE)

head(degs)

sig_degs_AB <- degs[degs$pvalue <= 0.05,,drop=FALSE]

degs[degs$padj <= 0.05,,drop=FALSE]

DE_AB_dat <- data %>% 
  rownames_to_column("rownames") %>% 
  right_join(sig_degs_AB%>% 
              rownames_to_column("rownames")) %>% 
  column_to_rownames("rownames") 

values <- DE_AB_dat %>% 
  select(-47:-52)%>% 
  as.matrix() 

values = jitter(values, factor = 1, amount = 0.00001)

zscores = NULL
for (i in 1 : nrow(values)) {
  row = values[i,]
  zrow = (row - mean(row, na.rm = TRUE)) / sd(row, na.rm = TRUE)
  zscores = rbind(zscores, zrow)
}


# Set the row names on the zscores.
row.names(zscores) = row.names(values)

# Turn the data into a matrix for heatmap2.
zscores = as.matrix(zscores) 


# Set the color palette.
col = greenred


test <- zscores %>% colnames %>% str_detect("B")
ifelse(test == TRUE, "yellow", "purple")


# Generate heatmap
heatmap_save <- heatmap.2(zscores,
                          na.rm=TRUE,
                          col=col,
                          density.info="none",
                          dendrogram="both",
                          trace="none",
                          ColSideColors=ifelse(test == TRUE, "yellow", "purple"),
                          # labRow=FALSE,
                          margins=MARGINS,
                          lhei=LHEI
)




save(trans_norm, samples, log2_norm, threshold, DE_AB_dat, data_cor, file = str_c(directory, "normalized_data_configurations.RData"))
