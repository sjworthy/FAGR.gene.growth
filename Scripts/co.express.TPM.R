#### Co-Expresion using TPMs ####

# load libraries
library(tidyverse)
library(WGCNA)
library(ggdendro)
library(gplots)
library(grid)
library(goseq)
library(matrixStats)
library(pvclust)

# Perform WGCNA using the TPM normalized data

# Calculating Coefficient of variation function
calc.cv <- function(x, na.rm=TRUE) { 
  if(na.rm==TRUE) x <- na.omit(x)
  result <- sd(x) / mean(x) 
  result <- abs(result) 
  return(result)
}

#### Read in our TPM normalized RNA-seq expression data ####

TPM_ExprData <- read_csv("./Raw.Data/HTseq-master-counts.csv")
# 295 samples

# split samples into time points
TPM_2017_Spring = TPM_ExprData[,c(1:21,108:127)] # 40 samples
TPM_2017_Summer = TPM_ExprData[,c(1,22:39,128:146)] # 37 samples
TPM_2017_Fall = TPM_ExprData[,c(1,40:59,147:166)] # 40 samples

TPM_2018_Spring = TPM_ExprData[,c(1,60:72,167:186)] # 33 samples
TPM_2018_Summer = TPM_ExprData[,c(1,73:89,187:206)] # 37 samples
TPM_2018_Fall = TPM_ExprData[,c(1,90:107,207:223)] # 35 samples

TPM_2019_Spring = TPM_ExprData[,c(1,224:240)] # 17 samples
TPM_2019_Summer = TPM_ExprData[,c(1,241:259)] # 19 samples
TPM_2019_Fall = TPM_ExprData[,c(1,260:276)] # 17 samples

TPM_2020_Spring = TPM_ExprData[,c(1,277:296)] # 20 samples

TPM_2017 = TPM_ExprData[,c(1:59,108:166)] # 117 samples
TPM_2018 = TPM_ExprData[,c(1,60:107,167:223)] # 105 samples
TPM_2019 = TPM_ExprData[,c(1,224:276)] # 53 samples


#### Calculate the Variance of each Gene across all of our samples ####
TPM_ExprData_cv <- TPM_ExprData %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Fall_cv <- TPM_2017_Fall %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Spring_cv <- TPM_2017_Spring %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Summer_cv <- TPM_2017_Summer %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Fall_cv <- TPM_2018_Fall %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Spring_cv <- TPM_2018_Spring %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Summer_cv <- TPM_2018_Summer %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2019_Fall_cv <- TPM_2019_Fall %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2019_Spring_cv <- TPM_2019_Spring %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2019_Summer_cv <- TPM_2019_Summer %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2020_Spring_cv <- TPM_2020_Spring %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

# Filter the data frame to contain the top 30% most variable genes)
# ties are kept together
TPM_ExprData_30  <- TPM_ExprData_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Fall_30  <- TPM_2017_Fall_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Spring_30  <- TPM_2017_Spring_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Summer_30  <- TPM_2017_Summer_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Fall_30  <- TPM_2018_Fall_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Spring_30  <- TPM_2018_Spring_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Summer_30  <- TPM_2018_Summer_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2019_Fall_30  <- TPM_2019_Fall_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2019_Spring_30  <- TPM_2019_Spring_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2019_Summer_30  <- TPM_2019_Summer_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2020_Spring_30  <- TPM_2020_Spring_cv  %>% slice_max(order_by = cv, prop = .30)

# Deselect the cv column and flip our data frame to contain sample along the rows and genes along the columns
TPM_ExprData_30 <- select(TPM_ExprData_30, -cv)
TPM_ExprData_30 <- column_to_rownames(TPM_ExprData_30,var = "Gene_ID")
TPM_ExprData_30 <- as.data.frame(t(TPM_ExprData_30))

TPM_2017_Fall_30 <- select(TPM_2017_Fall_30, -cv)
TPM_2017_Fall_30 <- column_to_rownames(TPM_2017_Fall_30,var = "Gene_ID")
TPM_2017_Fall_30 <- as.data.frame(t(TPM_2017_Fall_30))

TPM_2017_Spring_30 <- select(TPM_2017_Spring_30, -cv)
TPM_2017_Spring_30 <- column_to_rownames(TPM_2017_Spring_30,var = "Gene_ID")
TPM_2017_Spring_30 <- as.data.frame(t(TPM_2017_Spring_30))

TPM_2017_Summer_30 <- select(TPM_2017_Summer_30, -cv)
TPM_2017_Summer_30 <- column_to_rownames(TPM_2017_Summer_30,var = "Gene_ID")
TPM_2017_Summer_30 <- as.data.frame(t(TPM_2017_Summer_30))

TPM_2018_Fall_30 <- select(TPM_2018_Fall_30, -cv)
TPM_2018_Fall_30 <- column_to_rownames(TPM_2018_Fall_30,var = "Gene_ID")
TPM_2018_Fall_30 <- as.data.frame(t(TPM_2018_Fall_30))

TPM_2018_Spring_30 <- select(TPM_2018_Spring_30, -cv)
TPM_2018_Spring_30 <- column_to_rownames(TPM_2018_Spring_30,var = "Gene_ID")
TPM_2018_Spring_30 <- as.data.frame(t(TPM_2018_Spring_30))

TPM_2018_Summer_30 <- select(TPM_2018_Summer_30, -cv)
TPM_2018_Summer_30 <- column_to_rownames(TPM_2018_Summer_30,var = "Gene_ID")
TPM_2018_Summer_30 <- as.data.frame(t(TPM_2018_Summer_30))

TPM_2019_Fall_30 <- select(TPM_2019_Fall_30, -cv)
TPM_2019_Fall_30 <- column_to_rownames(TPM_2019_Fall_30,var = "Gene_ID")
TPM_2019_Fall_30 <- as.data.frame(t(TPM_2019_Fall_30))

TPM_2019_Spring_30 <- select(TPM_2019_Spring_30, -cv)
TPM_2019_Spring_30 <- column_to_rownames(TPM_2019_Spring_30,var = "Gene_ID")
TPM_2019_Spring_30 <- as.data.frame(t(TPM_2019_Spring_30))

TPM_2019_Summer_30 <- select(TPM_2019_Summer_30, -cv)
TPM_2019_Summer_30 <- column_to_rownames(TPM_2019_Summer_30,var = "Gene_ID")
TPM_2019_Summer_30 <- as.data.frame(t(TPM_2019_Summer_30))

TPM_2020_Spring_30 <- select(TPM_2020_Spring_30, -cv)
TPM_2020_Spring_30 <- column_to_rownames(TPM_2020_Spring_30,var = "Gene_ID")
TPM_2020_Spring_30 <- as.data.frame(t(TPM_2020_Spring_30))

# Quickly save gene set to use on DGE analysis
# not done right now
#Top_30_Most_Variable_Genes <- colnames(TPM_ExprData_30)
#saveRDS(Top_30_Most_Variable_Genes, file = "./Data/Top_30_Most_Variable_Genes.R")

#Top_30_Most_Variable_Genes_2017_Fall <- rownames(TPM_2017_Fall_30)
#saveRDS(Top_30_Most_Variable_Genes_2017_Fall, file = "./Data/Top_30_Most_Variable_Genes_2017_Fall.R")
#Top_30_Most_Variable_Genes_2017_Spring <- rownames(TPM_2017_Spring_30)
#saveRDS(Top_30_Most_Variable_Genes_2017_Spring, file = "./Data/Top_30_Most_Variable_Genes_2017_Spring.R")
#Top_30_Most_Variable_Genes_2017_Summer <- rownames(TPM_2017_Summer_30)
#saveRDS(Top_30_Most_Variable_Genes_2017_Summer, file = "./Data/Top_30_Most_Variable_Genes_2017_Summer.R")
#Top_30_Most_Variable_Genes_2018_Fall <- rownames(TPM_2018_Fall_30)
#saveRDS(Top_30_Most_Variable_Genes_2018_Fall, file = "./Data/Top_30_Most_Variable_Genes_2018_Fall.R")
#Top_30_Most_Variable_Genes_2018_Spring <- rownames(TPM_2018_Spring_30)
#saveRDS(Top_30_Most_Variable_Genes_2018_Spring, file = "./Data/Top_30_Most_Variable_Genes_2018_Spring.R")
#Top_30_Most_Variable_Genes_2018_Summer <- rownames(TPM_2018_Summer_30)
#saveRDS(Top_30_Most_Variable_Genes_2018_Summer, file = "./Data/Top_30_Most_Variable_Genes_2018_Summer.R")
#Top_30_Most_Variable_Genes_2019_Fall <- rownames(TPM_2019_Fall_30)
#saveRDS(Top_30_Most_Variable_Genes_2019_Fall, file = "./Data/Top_30_Most_Variable_Genes_2019_Fall.R")
#Top_30_Most_Variable_Genes_2019_Spring <- rownames(TPM_2019_Spring_30)
#saveRDS(Top_30_Most_Variable_Genes_2019_Spring, file = "./Data/Top_30_Most_Variable_Genes_2019_Spring.R")
#Top_30_Most_Variable_Genes_2019_Summer <- rownames(TPM_2019_Summer_30)
#saveRDS(Top_30_Most_Variable_Genes_2019_Summer, file = "./Data/Top_30_Most_Variable_Genes_2019_Summer.R")
#Top_30_Most_Variable_Genes_2020_Spring <- rownames(TPM_2020_Spring_30)
#saveRDS(Top_30_Most_Variable_Genes_2020_Spring, file = "./Data/Top_30_Most_Variable_Genes_2020_Spring.R")

# Continue with formatting
Sample_Description <-
  read_csv("./Formatted.Data/FAGR.description.csv")

# rename dataframe
datExpr <- TPM_ExprData_30
datExpr
datExpr_2017_Fall <- TPM_2017_Fall_30
datExpr_2017_Spring <- TPM_2017_Spring_30
datExpr_2017_Summer <- TPM_2017_Summer_30
datExpr_2018_Fall <- TPM_2018_Fall_30
datExpr_2018_Spring <- TPM_2018_Spring_30
datExpr_2018_Summer <- TPM_2018_Summer_30
datExpr_2019_Fall <- TPM_2019_Fall_30
datExpr_2019_Spring <- TPM_2019_Spring_30
datExpr_2019_Summer <- TPM_2019_Summer_30
datExpr_2020_Spring <- TPM_2020_Spring_30

#### QC ####
# iterative filtering of samples and genes with too many missing entries
gsg = goodSamplesGenes(TPM_ExprData_30, verbose = 3);
gsg$allOK

gsg_2017_Fall = goodSamplesGenes(TPM_2017_Fall_30, verbose = 3);
gsg_2017_Fall$allOK
gsg_2017_Spring = goodSamplesGenes(TPM_2017_Spring_30, verbose = 3);
gsg_2017_Spring$allOK
gsg_2017_Summer = goodSamplesGenes(TPM_2017_Summer_30, verbose = 3);
gsg_2017_Summer$allOK

gsg_2018_Fall = goodSamplesGenes(TPM_2018_Fall_30, verbose = 3);
gsg_2018_Fall$allOK
gsg_2018_Spring = goodSamplesGenes(TPM_2018_Spring_30, verbose = 3);
gsg_2018_Spring$allOK
gsg_2018_Summer = goodSamplesGenes(TPM_2018_Summer_30, verbose = 3);
gsg_2018_Summer$allOK

gsg_2019_Fall = goodSamplesGenes(TPM_2019_Fall_30, verbose = 3);
gsg_2019_Fall$allOK
gsg_2019_Spring = goodSamplesGenes(TPM_2019_Spring_30, verbose = 3);
gsg_2019_Spring$allOK
gsg_2019_Summer = goodSamplesGenes(TPM_2019_Summer_30, verbose = 3);
gsg_2019_Summer$allOK

gsg_2020_Spring = goodSamplesGenes(TPM_2020_Spring_30, verbose = 3);
gsg_2020_Spring$allOK

#### First Tree ####
datExpr2 <- datExpr
sampleTree = hclust(dist(datExpr2), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12, 9)
par(cex = 0.6,
    mar = c(0, 4, 2, 0))
jpeg("./Plots/Full_Basic_Tree.jpeg",height = 1080, width = 1920)
plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)
#Plot a line to show the cut
abline(h = 180, col = "red") # what does this mean

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 180, minSize = 10)
table(clust)
dev.off()

#### Second Tree ####
datExpr2 <- datExpr
sampleTree <- hclust(dist(datExpr2), method = "average")
jpeg("./Plots/Full_Boxed_Basic_Tree.jpeg",height = 1080, width = 1920)
ggdendrogram(sampleTree)
plot(sampleTree, cex=.6) #redraw the tree everytime before adding the rectangles
rect.hclust(sampleTree, k = 5, border = "red")
dev.off()

#### HeatMap ####

sampleTree <- hclust(dist(datExpr2), method = "average")
heatmap.2(as.matrix(datExpr2), Rowv = as.dendrogram(sampleTree),  
          density.info="none", trace="none", margins = c(10,5))

sampleTree_2017_Fall <- hclust(dist(datExpr_2017_Fall), method = "average")
heatmap.2(as.matrix(datExpr_2017_Fall), Rowv = as.dendrogram(sampleTree_2017_Fall),  
          density.info="none", trace="none", margins = c(10,5))
sampleTree_2017_Spring <- hclust(dist(datExpr_2017_Spring), method = "average")
heatmap.2(as.matrix(datExpr_2017_Spring), Rowv = as.dendrogram(sampleTree_2017_Spring),  
          density.info="none", trace="none", margins = c(10,5))
sampleTree_2017_Summer <- hclust(dist(datExpr_2017_Summer), method = "average")
heatmap.2(as.matrix(datExpr_2017_Summer), Rowv = as.dendrogram(sampleTree_2017_Summer),  
          density.info="none", trace="none", margins = c(10,5))
sampleTree_2018_Fall <- hclust(dist(datExpr_2018_Fall), method = "average")
heatmap.2(as.matrix(datExpr_2018_Fall), Rowv = as.dendrogram(sampleTree_2018_Fall),  
          density.info="none", trace="none", margins = c(10,5))
sampleTree_2018_Spring <- hclust(dist(datExpr_2018_Spring), method = "average")
heatmap.2(as.matrix(datExpr_2018_Spring), Rowv = as.dendrogram(sampleTree_2018_Spring),  
          density.info="none", trace="none", margins = c(10,5))
sampleTree_2018_Summer <- hclust(dist(datExpr_2018_Summer), method = "average")
heatmap.2(as.matrix(datExpr_2018_Summer), Rowv = as.dendrogram(sampleTree_2018_Summer),  
          density.info="none", trace="none", margins = c(10,5))
sampleTree_2019_Fall <- hclust(dist(datExpr_2019_Fall), method = "average")
heatmap.2(as.matrix(datExpr_2019_Fall), Rowv = as.dendrogram(sampleTree_2019_Fall),  
          density.info="none", trace="none", margins = c(10,5))
sampleTree_2019_Spring <- hclust(dist(datExpr_2019_Spring), method = "average")
heatmap.2(as.matrix(datExpr_2019_Spring), Rowv = as.dendrogram(sampleTree_2019_Spring),  
          density.info="none", trace="none", margins = c(10,5))
sampleTree_2019_Summer <- hclust(dist(datExpr_2019_Summer), method = "average")
heatmap.2(as.matrix(datExpr_2019_Summer), Rowv = as.dendrogram(sampleTree_2019_Summer),  
          density.info="none", trace="none", margins = c(10,5))
sampleTree_2020_Spring <- hclust(dist(datExpr_2020_Spring), method = "average")
heatmap.2(as.matrix(datExpr_2020_Spring), Rowv = as.dendrogram(sampleTree_2020_Spring),  
          density.info="none", trace="none", margins = c(10,5))

#### Soft Threshold ####
# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft_2017_Fall <- pickSoftThreshold(datExpr_2017_Fall, powerVector = powers, verbose = 5)
sft_2017_Spring <- pickSoftThreshold(datExpr_2017_Spring, powerVector = powers, verbose = 5)
sft_2017_Summer <- pickSoftThreshold(datExpr_2017_Summer, powerVector = powers, verbose = 5)
sft_2018_Fall <- pickSoftThreshold(datExpr_2018_Fall, powerVector = powers, verbose = 5)
sft_2018_Spring <- pickSoftThreshold(datExpr_2018_Spring, powerVector = powers, verbose = 5)
sft_2018_Summer <- pickSoftThreshold(datExpr_2018_Summer, powerVector = powers, verbose = 5)
sft_2019_Fall <- pickSoftThreshold(datExpr_2019_Fall, powerVector = powers, verbose = 5)
sft_2019_Spring <- pickSoftThreshold(datExpr_2019_Spring, powerVector = powers, verbose = 5)
sft_2019_Summer <- pickSoftThreshold(datExpr_2019_Summer, powerVector = powers, verbose = 5)
sft_2020_Spring <- pickSoftThreshold(datExpr_2020_Spring, powerVector = powers, verbose = 5)

#### Soft Threshold Plotting ####
# Plot the results:
jpeg("./Plots/soft-thresholding.jpeg",height = 1080, width = 1920)
par(mfrow = c(1, 2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit,signed R^2",
  type = "n",
  main = paste("Scale independence")
)

text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)

# this line corresponds to using an R^2 cut-off of h
abline(h = 0.90, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)
dev.off()

#### Adjacency ####
# signed hybrid sets all negatively correlated genes to 0
# choice of power comes from soft threshold analysis above

# only all data, 2018_Spring, 2020_Spring had values above 0.90

adjacency = adjacency(datExpr, power = 3,type = "signed hybrid")
adjacency_2017_Fall = adjacency(datExpr_2017_Fall, power = 5,type = "signed hybrid")
adjacency_2017_Spring = adjacency(datExpr_2017_Spring, power = 7,type = "signed hybrid")
adjacency_2017_Summer = adjacency(datExpr_2017_Summer, power = 5,type = "signed hybrid")
adjacency_2018_Fall = adjacency(datExpr_2018_Fall, power = 5,type = "signed hybrid")
adjacency_2018_Spring = adjacency(datExpr_2018_Spring, power = 7,type = "signed hybrid")
adjacency_2018_Summer = adjacency(datExpr_2018_Summer, power = 5,type = "signed hybrid")
adjacency_2019_Fall = adjacency(datExpr_2019_Fall, power = 5,type = "signed hybrid")
adjacency_2019_Spring = adjacency(datExpr_2019_Spring, power = 4,type = "signed hybrid")
adjacency_2019_Summer = adjacency(datExpr_2019_Summer, power = 6,type = "signed hybrid")
adjacency_2020_Spring = adjacency(datExpr_2020_Spring, power = 3,type = "signed hybrid")

#### Topological overlap matrix (TOM) ####
# need to rm adjacency matrices and TOM after use to make space

TOM = TOMsimilarity(adjacency, TOMType = "signed Nowick")
dissTOM = 1 - TOM
#saveRDS(TOM, file = "./Data/TOM.R")
TOM_2017_Fall = TOMsimilarity(adjacency_2017_Fall, TOMType = "signed Nowick")
dissTOM_2017_Fall = 1 - TOM_2017_Fall
#saveRDS(TOM_2017_Fall, file = "./Data/TOM_2017_Fall.R")
TOM_2017_Spring = TOMsimilarity(adjacency_2017_Spring, TOMType = "signed Nowick")
dissTOM_2017_Spring = 1 - TOM_2017_Spring
#saveRDS(TOM_2017_Spring, file = "./Data/TOM_2017_Spring.R")
TOM_2017_Summer = TOMsimilarity(adjacency_2017_Summer, TOMType = "signed Nowick")
dissTOM_2017_Summer = 1 - TOM_2017_Summer
#saveRDS(TOM_2017_Summer, file = "./Data/TOM_2017_Summer.R")
TOM_2018_Fall = TOMsimilarity(adjacency_2018_Fall, TOMType = "signed Nowick")
dissTOM_2018_Fall = 1 - TOM_2018_Fall
#saveRDS(TOM_2018_Fall, file = "./Data/TOM_2018_Fall.R")
TOM_2018_Spring = TOMsimilarity(adjacency_2018_Spring, TOMType = "signed Nowick")
dissTOM_2018_Spring = 1 - TOM_2018_Spring
#saveRDS(TOM_2018_Spring, file = "./Data/TOM_2018_Spring.R")
TOM_2018_Summer = TOMsimilarity(adjacency_2018_Summer, TOMType = "signed Nowick")
dissTOM_2018_Summer = 1 - TOM_2018_Summer
#saveRDS(TOM_2018_Summer, file = "./Data/TOM_2018_Summer.R")
TOM_2019_Fall = TOMsimilarity(adjacency_2019_Fall, TOMType = "signed Nowick")
dissTOM_2019_Fall = 1 - TOM_2019_Fall
#saveRDS(TOM_2019_Fall, file = "./Data/TOM_2019_Fall.R")
TOM_2019_Spring = TOMsimilarity(adjacency_2019_Spring, TOMType = "signed Nowick")
dissTOM_2019_Spring = 1 - TOM_2019_Spring
#saveRDS(TOM_2019_Spring, file = "./Data/TOM_2019_Spring.R")
TOM_2019_Summer = TOMsimilarity(adjacency_2019_Summer, TOMType = "signed Nowick")
dissTOM_2019_Summer = 1 - TOM_2019_Summer
#saveRDS(TOM_2019_Summer, file = "./Data/TOM_2019_Summer.R")
TOM_2020_Spring = TOMsimilarity(adjacency_2020_Spring, TOMType = "signed Nowick")
dissTOM_2020_Spring = 1 - TOM_2020_Spring
#saveRDS(TOM_2020_Spring, file = "./Data/TOM_2020_Spring.R")

#### Gene Clustering Plot ####
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12, 9)
jpeg("./Plots/dissTOM-tree.jpeg",height = 1080, width = 1920)
plot(
  geneTree,
  xlab = "",
  sub = "",
  main = "Gene clustering on TOM-based dissimilarity",
  labels = FALSE,
  hang = 0.04
)

#### Build Moduldes ####
# change module size depending on results?
# We like large modules, so we set the minimum module size relatively high:
minModuleSize <- 30

# Module identification using dynamic tree cut:
# adaptive branch pruning of hierarchical clustering dendrograms
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  method = "hybrid",
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
)

table(dynamicMods)

#### Plot Tree with Modules ####
# Convert numeric labels into colors
set.seed(1232)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath

sizeGrWindow(8, 6)
jpeg("./Plots/dissTOM-tree-blocks.jpeg",height = 1080, width = 1920)
plotDendroAndColors(
  geneTree,
  dynamicColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)

#### Eigen Clustering ####
# Calculate eigengenes
# calculates module eigengenes (1st PC) of modules in a given dataset
# softPower needs to be changed here to match adjacency analysis above
MEList <- moduleEigengenes(datExpr,
                           colors = dynamicColors,
                           softPower = 3)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(7, 6)
jpeg("./Plots/ME-tree.jpeg",height = 1080, width = 1920)

plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")

MEDissThres <-  0.25
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")

# Call an automatic merging function
merge <- mergeCloseModules(datExpr,
                           dynamicColors,
                           cutHeight = MEDissThres,
                           verbose = 3)
# The merged module colors
mergedColors <- merge$colors

# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs
dev.off()

sizeGrWindow(12, 9)
pdf(file = "./Plots/geneDendro.pdf", wi = 9, he = 6)
plotDendroAndColors(
  geneTree,
  cbind(dynamicColors, mergedColors),
  c("Dynamic Tree Cut", "Merged dynamic"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()

#### Color Renaming ####
# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))

moduleLabels <- match(moduleColors, colorOrder) - 1

MEs <- mergedMEs

#### TOM plotting ####
#plotTOM <- TOM

# Set diagonal to NA for a nicer plot
#diag(plotTOM) <- NA

# Call the plot function
#sizeGrWindow(9, 9)
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

#### Creating a Module Heatmap of gene Expression ####
Module_Heatmap <- datExpr2 %>%
  rownames_to_column(var = "sample") %>%
  #mutate(Sample_Group = str_replace(string = sample, "(.*)_(.*)_(.*)", "\\1_\\2")) %>%
  select(-sample) %>%
  #group_by(Sample_Group) %>%
  summarize(across(.cols = everything(),
                   .fns = mean)) %>%
  column_to_rownames(var = "Sample_Group")


ModuleGenes_Df <- data.frame(gene = NULL, color = NULL)
for(color in unique(moduleColors)) {
  temp <- data.frame(names(datExpr2)[moduleColors == color])
  names(temp)[1] <- "gene"
  temp$color <- color
  ModuleGenes_Df <- rbind(ModuleGenes_Df, temp)
}

#Module_Heatmap <- Module_Heatmap %>%
  #mutate(Treatment = str_replace(rownames(.), "(.*)_(.*)", "\\2")) %>% 
  #arrange(Treatment) %>%
  #select(-Treatment)

Module_Matrix <- as.matrix(Module_Heatmap)

reorder.idx <- match(colnames(Module_Matrix), ModuleGenes_Df$gene)
ModuleGenes_Df <- ModuleGenes_Df[reorder.idx, ]

my.pallete <-
  scales::div_gradient_pal(low = "red", mid = "white", high = "blue")
my.colors <- my.pallete(seq(0, 1, length.out = 9))
#my.breaks <- c(-1, -.5, .5, 1)

# play with margins here
jpeg("./Plots/ModuleGene_Heatmap.jpeg",
     height = 1080,
     width = 1920)
heatmap.2(Module_Matrix,
          scale = "col",
          density.info="none",
          trace="none",
          dendrogram = "col",
          col=my.colors,
          #breaks = my.breaks,
          Colv = as.dendrogram(geneTree),
          margins = c(40,40),
          Rowv = NULL,
          cexRow = 2.3,
          cexCol,
          ColSideColors = moduleColors,
          labCol = F)
dev.off()

#### formatting eigen clusters with sample descriptions ####
growth_MEs <- MEs %>%
  rownames_to_column(var = "sample") %>%
  as_tibble()

exp <-  Sample_Description %>%
  select(sample.description, Site, Year, Season, Tree_ID)
exp <- exp[match(growth_MEs$sample, exp$sample.description),]

growth_MEs
exp

#### Bring in growth data ####
growth.data <- read_csv("./Formatted.Data/Dendro_FAGR.csv") %>%
  #filter(sample %in% exp$sample) %>%
  select(SITE, YEAR, TREE_ID, RGR, GR, Max.growth.day)

exp <- exp %>%
  left_join(growth.data, by = "sample.description")
exp

#### Exporting Eigengene values ####
growth_MEs <- growth_MEs %>% 
  rename_all(~ str_replace(., "ME","ME_"))
  #left_join(exp[, c("sample.description")], by = "sample") %>%
  #select(sample, sample.description, everything())
write_csv(growth_MEs,"./Data/all.samples.MEs.csv")






