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

TPM_ExprData <- read.csv("./Raw.Data/HTseq-master-counts.csv", row.names = 1)
# 297 samples

# split samples into time points

TPM_2017_Fall
TPM_2017_Spring
TPM_2017_Summer

TPM_2018_Fall
TPM_2018_Spring
TPM_2018_Summer

TPM_2019_Fall
TPM_2019_Spring
TPM_2019_Summer

TPM_2020_Spring

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
Top_30_Most_Variable_Genes <- rownames(TPM_ExprData_30)
saveRDS(Top_30_Most_Variable_Genes, file = "./Data/Top_30_Most_Variable_Genes.R")

Top_30_Most_Variable_Genes_2017_Fall <- rownames(TPM_2017_Fall_30)
saveRDS(Top_30_Most_Variable_Genes_2017_Fall, file = "./Data/Top_30_Most_Variable_Genes_2017_Fall.R")
Top_30_Most_Variable_Genes_2017_Spring <- rownames(TPM_2017_Spring_30)
saveRDS(Top_30_Most_Variable_Genes_2017_Spring, file = "./Data/Top_30_Most_Variable_Genes_2017_Spring.R")
Top_30_Most_Variable_Genes_2017_Summer <- rownames(TPM_2017_Summer_30)
saveRDS(Top_30_Most_Variable_Genes_2017_Summer, file = "./Data/Top_30_Most_Variable_Genes_2017_Summer.R")
Top_30_Most_Variable_Genes_2018_Fall <- rownames(TPM_2018_Fall_30)
saveRDS(Top_30_Most_Variable_Genes_2018_Fall, file = "./Data/Top_30_Most_Variable_Genes_2018_Fall.R")
Top_30_Most_Variable_Genes_2018_Spring <- rownames(TPM_2018_Spring_30)
saveRDS(Top_30_Most_Variable_Genes_2018_Spring, file = "./Data/Top_30_Most_Variable_Genes_2018_Spring.R")
Top_30_Most_Variable_Genes_2018_Summer <- rownames(TPM_2018_Summer_30)
saveRDS(Top_30_Most_Variable_Genes_2018_Summer, file = "./Data/Top_30_Most_Variable_Genes_2018_Summer.R")
Top_30_Most_Variable_Genes_2019_Fall <- rownames(TPM_2019_Fall_30)
saveRDS(Top_30_Most_Variable_Genes_2019_Fall, file = "./Data/Top_30_Most_Variable_Genes_2019_Fall.R")
Top_30_Most_Variable_Genes_2019_Spring <- rownames(TPM_2019_Spring_30)
saveRDS(Top_30_Most_Variable_Genes_2019_Spring, file = "./Data/Top_30_Most_Variable_Genes_2019_Spring.R")
Top_30_Most_Variable_Genes_2019_Summer <- rownames(TPM_2019_Summer_30)
saveRDS(Top_30_Most_Variable_Genes_2019_Summer, file = "./Data/Top_30_Most_Variable_Genes_2019_Summer.R")
Top_30_Most_Variable_Genes_2020_Spring <- rownames(TPM_2020_Spring_30)
saveRDS(Top_30_Most_Variable_Genes_2020_Spring, file = "./Data/Top_30_Most_Variable_Genes_2020_Spring.R")

# Continue with formatting
Sample_Description <-
  read_csv("./FAGR.description.csv")

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
rownames(datExpr2) <- Sample_Description %>%
  filter(sample %in% rownames(datExpr2)) %>%
  pull(sample.description)
sampleTree = hclust(dist(datExpr2), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12, 9)
par(cex = 0.6,
    mar = c(0, 4, 2, 0))
jpeg("../Plots/Full_Basic_Tree.jpeg",height = 1080, width = 1920)
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
abline(h = 180, col = "red")

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 180, minSize = 10)
table(clust)
dev.off()

#### Second Tree ####
datExpr2 <- datExpr %>%
  filter(grepl("^a|b|c|d", rownames(.)))
rownames(datExpr2) <- Sample_Description %>%
  filter(sample %in% rownames(datExpr2)) %>%
  pull(sample.description)
sampleTree <- hclust(dist(datExpr2), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12, 9)
par(cex = 0.6,
    mar = c(0, 4, 2, 0))
jpeg("../Plots/Caulanthus_Basic_Tree.jpeg",height = 1080, width = 1920)
plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)
dev.off()
