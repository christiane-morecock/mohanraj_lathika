## Load Libraries

library(tidyverse)
library(VennDiagram)

## Read in Data
threshold18 <- read.csv("./threshold1_18_removeNegC/DESeq2_DEGs_Avs.B_final.csv")
threshold22 <- read.csv("./threshold1_22/DESeq2_DEGs__Avs.B_prelim.csv")
threshold28 <- read.csv("./threshold2_28/DESeq2_DEGs_Avs.B_final.csv")


sigThreshold18 <- threshold18[threshold18$pvalue <= 0.05,,drop=FALSE] %>% as_tibble()
sigThreshold22 <- threshold22[threshold22$pvalue <= 0.05,,drop=FALSE] %>% as_tibble()
sigThreshold28 <- threshold28[threshold28$pvalue <= 0.05,,drop=FALSE] %>% as_tibble()


sigThreshold18 <- sigThreshold18 %>% 
  filter(padj < 0.05)
# 4 genes 

sigThreshold22 <- sigThreshold22 %>% 
  filter(padj < 0.05)
# 4 genes

sigThreshold28 <- sigThreshold28 %>% 
  filter(padj < 0.05)
# 5 genes


# Load library
library(VennDiagram)

# Generate 3 sets of 200 words
set1 <- sigThreshold18$X
set2 <- sigThreshold22$X
set3 <- sigThreshold28$X

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Threshold: 18" , "Threshold: 22" , "Threshold: 28"),
  filename = 'compareThresholdsVenn_padj.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 600,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.2,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


intersect(sigThreshold18$X, sigThreshold22$X) %>% 
  intersect(sigThreshold28$X)

