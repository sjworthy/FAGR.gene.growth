#### Co-Expresion Analysis ####

# load libraries
library(tidyverse)
library(WGCNA)
library(ggdendro)
library(gplots)
library(grid)
library(goseq)
library(matrixStats)
library(pvclust)
library(edgeR)

# Using WCGNA to relate differences in expression between years to differences in growth between years
# Doing this: 
# 1) by subtracting log gene expression from TMM file to get log2FC values, 
# 2) same thing as 1, but only significantly differentially expressed genes are considered, not all genes in WGCNA 

#### Subtract TMM values among years in data sets #####



# Subset TMM by significantly differently expressed genes or maybe all genes and take difference between samples,
# so differences in log2FC change between samples between years.
