# Code for building coexpression networks

# installing WGCNA
install.packages("BiocManager")
BiocManager::install("WGCNA", force = T)

# Load Libraries
library(tidyverse)
library(pvclust)
library(WGCNA)
library(ggdendro)
library(gplots)
library(grid)
library(goseq)
library(matrixStats)

# Calculating Variance function
calc.cv <- function(x, na.rm=TRUE) { 
  if(na.rm==TRUE) x <- na.omit(x)
  result <- sd(x) / mean(x) 
  result <- abs(result) 
  return(result)
}

# Read in our TPM normalized RNA-seq expression data to R
# This is just the counts data that has been normalized previously
TPM_ExprData <- read.csv("./Raw.Data/HTseq-master-counts.csv", row.names = 1)

# Calculate the Variance of each Gene across all of our samples

TPM_ExprData <-rownames_to_column(TPM_ExprData, var = "GeneID")
TPM_ExprData <- TPM_ExprData %>% 
  rowwise() %>% 
  mutate(cv = calc.cv(c_across(-GeneID)))%>% 
  ungroup() %>% select(GeneID, cv, everything())

# filter genes with NA for cv value

TPM_ExprData_noNA = subset(TPM_ExprData, !is.na(TPM_ExprData$cv))

# not sure the standard number of most variable genes to select
# going to try John's top 40% and Julin's top 1000 genes

# Filter the data frame to contain the top 40% or so most variable genes
n <- 40
TPM_ExprData_40 = TPM_ExprData_noNA[TPM_ExprData_noNA$cv > quantile(TPM_ExprData_noNA$cv,prob=1-n/100),]

# Filter the data frame for top 1000 genes
TPM_ExprData_1000 = TPM_ExprData_noNA %>%
  arrange(desc(cv)) %>%
  slice(1:1000)

# scale and center the data so each gene has mean of 0 and sd of 1
# prevents genes with high expression from having undue influence
# convert to matrix for downstream functions

TPM_ExprData_40_matrix <- TPM_ExprData_40 %>% 
  select(-GeneID, -cv) %>% 
  as.matrix() %>% 
  scale()

TPM_ExprData_1000_matrix <- TPM_ExprData_1000 %>% 
  select(-GeneID, -cv) %>% 
  as.matrix() %>% 
  scale()
# looks like 1000 doesn't work b/c many have the same cv value for NA for mean and sd

# calculate distance, larger distance = more dissimilar RNA expression value
# euclidean distance

# cluster by genes
gene_hclust_row_40 <- TPM_ExprData_40_matrix %>% dist() %>% hclust()
ggdendrogram(gene_hclust_row_40) # looks weird

gene_hclust_row_1000 <- TPM_ExprData_1000_matrix %>% dist() %>% hclust()
ggdendrogram(gene_hclust_row_1000) # looks weird also

# cluster by sample
gene_hclust_col_40 <- TPM_ExprData_40_matrix %>% t() %>% dist() %>% hclust()
ggdendrogram(gene_hclust_col_40) # looks like 4 clusters

ggdendrogram(gene_hclust_col_40) + 
  theme(text = element_text(size = 5))

ddata <- dendro_data(as.dendrogram(gene_hclust_col_40))

gene_hclust_col_1000 <- TPM_ExprData_1000_matrix %>% t() %>% dist() %>% hclust()
ggdendrogram(gene_hclust_col_1000) # not working b/c of NA

# define subclusters within the tree
plot(gene_hclust_col_40, cex=.6) #redraw the tree everytime before adding the rectangles
rect.hclust(gene_hclust_col_40, k = 6, border = "red")
# looks like mainly spring, summer, fall clusters

plot(gene_hclust_col_40, cex=.6) #redraw the tree everytime before adding the rectangles
rect.hclust(gene_hclust_col_40, k = 5, border = "red") # looks the best

plot(gene_hclust_col_40, cex=.6) #redraw the tree everytime before adding the rectangles
rect.hclust(gene_hclust_col_40, k = 4, border = "red")

# use pvclust to get p-values for each cluster

fit_40 <- pvclust(TPM_ExprData_40_matrix, method.hclust = "ward.D", method.dist = "euclidean", nboot = 50)
plot(fit_40, print.num=FALSE)

fit_40_1000 <- pvclust(TPM_ExprData_40_matrix, method.hclust = "ward.D", method.dist = "euclidean", nboot = 1000)
plot(fit_40_1000, print.num=FALSE)

# Deselect the cv column and flip our data frame to contain sample along the 
# rows and genes along the columns

TPM_ExprData_40 <- select(TPM_ExprData_40, -cv)
TPM_ExprData_40 <- column_to_rownames(TPM_ExprData_40,var = "GeneID")
TPM_ExprData_40<- as.data.frame(t(TPM_ExprData_40))

TPM_ExprData_1000 <- select(TPM_ExprData_1000, -cv)
TPM_ExprData_1000 <- column_to_rownames(TPM_ExprData_1000,var = "GeneID")
TPM_ExprData_1000<- as.data.frame(t(TPM_ExprData_1000))




#Load descriptions

Hydrothermal.description <-
  read_csv("Hydrothermal_Description.csv")
datExpr <- TMM_ExprData
datExpr


