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
TMM_2017.18 = read.csv("./Data/DE.data/TMM_NormData_LogCPM_2017_2018.csv")
TMM_2017.18.sub = TMM_2017.18[,c(46:89,141:191)] - TMM_2017.18[c(2:45,90:140)]
TMM_2017.18.sub = cbind(TMM_2017.18$Gene_ID,TMM_2017.18.sub)
colnames(TMM_2017.18.sub)[1] = "Gene_ID"

TMM_2018.19 = read.csv("./Data/DE.data/TMM_NormData_LogCPM_2018_2019.csv")
TMM_2018.19.sub = TMM_2018.19[,c(50:97)] - TMM_2018.19[,c(2:49)]
TMM_2018.19.sub = cbind(TMM_2018.19$Gene_ID,TMM_2018.19.sub)
colnames(TMM_2018.19.sub)[1] = "Gene_ID"

TMM_2017.19 = read.csv("./Data/DE.data/TMM_NormData_LogCPM_2017_2019.csv")
TMM_2017.19.sub = TMM_2017.19[,c(49:95)] - TMM_2017.19[,c(2:48)]
TMM_2017.19.sub = cbind(TMM_2017.19$Gene_ID,TMM_2017.19.sub)
colnames(TMM_2017.19.sub)[1] = "Gene_ID"

TMM_2019.20 = read.csv("./Data/DE.data/TMM_NormData_LogCPM_2019_2020.csv")
TMM_2019.20.sub = TMM_2019.20[,c(14:25)] - TMM_2019.20[,c(2:13)]
TMM_2019.20.sub = cbind(TMM_2019.20$Gene_ID,TMM_2019.20.sub)
colnames(TMM_2019.20.sub)[1] = "Gene_ID"

TMM_2017.18.HF = read.csv("./Data/DE.data/TMM_NormData_LogCPM_2017_2018.HF.csv")
TMM_2017.18.HF.sub = TMM_2017.18.HF[,c(46:89)] - TMM_2017.18.HF[,c(2:45)]
TMM_2017.18.HF.sub = cbind(TMM_2017.18.HF$Gene_ID,TMM_2017.18.HF.sub)
colnames(TMM_2017.18.HF.sub)[1] = "Gene_ID"

TMM_2017.18.SERC = read.csv("./Data/DE.data/TMM_NormData_LogCPM_2017_2018.SERC.csv")
TMM_2017.18.SERC.sub = TMM_2017.18.SERC[,c(53:103)] - TMM_2017.18.SERC[,c(2:52)]
TMM_2017.18.SERC.sub = cbind(TMM_2017.18.SERC$Gene_ID,TMM_2017.18.SERC.sub)
colnames(TMM_2017.18.SERC.sub)[1] = "Gene_ID"

TMM_2017.20 = read.csv("./Data/DE.data/TMM_NormData_LogCPM_2017_2020.csv")
TMM_2017.20.sub = TMM_2017.20[,c(15:27)] - TMM_2017.20[,c(2:14)]
TMM_2017.20.sub = cbind(TMM_2017.20$Gene_ID,TMM_2017.20.sub)
colnames(TMM_2017.20.sub)[1] = "Gene_ID"

TMM_2018.20 = read.csv("./Data/DE.data/TMM_NormData_LogCPM_2018_2020.csv")
TMM_2018.20.sub = TMM_2018.20[,c(15:27)] - TMM_2018.20[,c(2:14)]
TMM_2018.20.sub = cbind(TMM_2018.20$Gene_ID,TMM_2018.20.sub)
colnames(TMM_2018.20.sub)[1] = "Gene_ID"

TMM_2017.18.Spring = read.csv("./Data/DE.data/TMM_NormData_LogCPM_2017_2018.Spring.csv")
TMM_2017.18.Spring.sub = TMM_2017.18.Spring[,c(14:25,44:61)] - TMM_2017.18.Spring[,c(2:13,26:43)]
TMM_2017.18.Spring.sub = cbind(TMM_2017.18.Spring$Gene_ID,TMM_2017.18.Spring.sub)
colnames(TMM_2017.18.Spring.sub)[1] = "Gene_ID"

TMM_2017.18.Summer = read.csv("./Data/DE.data/TMM_NormData_LogCPM_2017_2018.Summer.csv")
TMM_2017.18.Summer.sub = TMM_2017.18.Summer[,c(18:33,51:67)] - TMM_2017.18.Summer[,c(2:17,34:50)]
TMM_2017.18.Summer.sub = cbind(TMM_2017.18.Summer$Gene_ID,TMM_2017.18.Summer.sub)
colnames(TMM_2017.18.Summer.sub)[1] = "Gene_ID"

TMM_2017.18.Fall= read.csv("./Data/DE.data/TMM_NormData_LogCPM_2017_2018.Fall.csv")
TMM_2017.18.Fall.sub = TMM_2017.18.Fall[,c(18:33,50:65)] - TMM_2017.18.Fall[,c(2:17,34:49)]
TMM_2017.18.Fall.sub = cbind(TMM_2017.18.Fall$Gene_ID,TMM_2017.18.Fall.sub)
colnames(TMM_2017.18.Fall.sub)[1] = "Gene_ID"

#### Calculate the Variance of each Gene across all of our samples ####

# Calculating Coefficient of variation function
calc.cv <- function(x, na.rm=TRUE) { 
  if(na.rm==TRUE) x <- na.omit(x)
  result <- sd(x) / mean(x) 
  result <- abs(result) 
  return(result)
}

TMM_2017.18_cv <- TMM_2017.18.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2018.19_cv <- TMM_2018.19.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2017.19_cv <- TMM_2017.19.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2019.20_cv <- TMM_2019.20.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2017.18.HF_cv <- TMM_2017.18.HF.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2017.18.SERC_cv <- TMM_2017.18.SERC.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2017.20_cv <- TMM_2017.20.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2018.20_cv <- TMM_2018.20.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2017.18.Spring_cv <- TMM_2017.18.Spring.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2017.18.Summer_cv <- TMM_2017.18.Summer.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2017.18.Fall_cv <- TMM_2017.18.Fall.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

# Filter the data frame to contain the top 30% most variable genes)
# ties are kept together
TMM_2017.18_30  <- TMM_2017.18_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2018.19_30  <- TMM_2018.19_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2017.19_30  <- TMM_2017.19_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2019.20_30  <- TMM_2019.20_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2017.18.HF_30  <- TMM_2017.18.HF_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2017.18.SERC_30  <- TMM_2017.18.SERC_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2017.20_30  <- TMM_2017.20_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2018.20_30  <- TMM_2018.20_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2017.18.Spring_30  <- TMM_2017.18.Spring_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2017.18.Summer_30  <- TMM_2017.18.Summer_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2017.18.Fall_30  <- TMM_2017.18.Fall_cv  %>% slice_max(order_by = cv, prop = .30)

# Deselect the cv column and flip our data frame to contain sample along the rows and genes along the columns
TMM_2017.18_30 <- select(TMM_2017.18_30, -cv)
TMM_2017.18_30 <- column_to_rownames(TMM_2017.18_30,var = "Gene_ID")
TMM_2017.18_30 <- as.data.frame(t(TMM_2017.18_30))

TMM_2018.19_30 <- select(TMM_2018.19_30, -cv)
TMM_2018.19_30 <- column_to_rownames(TMM_2018.19_30,var = "Gene_ID")
TMM_2018.19_30 <- as.data.frame(t(TMM_2018.19_30))

TMM_2017.19_30 <- select(TMM_2017.19_30, -cv)
TMM_2017.19_30 <- column_to_rownames(TMM_2017.19_30,var = "Gene_ID")
TMM_2017.19_30 <- as.data.frame(t(TMM_2017.19_30))

TMM_2019.20_30 <- select(TMM_2019.20_30, -cv)
TMM_2019.20_30 <- column_to_rownames(TMM_2019.20_30,var = "Gene_ID")
TMM_2019.20_30 <- as.data.frame(t(TMM_2019.20_30))

TMM_2017.18.HF_30 <- select(TMM_2017.18.HF_30, -cv)
TMM_2017.18.HF_30 <- column_to_rownames(TMM_2017.18.HF_30,var = "Gene_ID")
TMM_2017.18.HF_30 <- as.data.frame(t(TMM_2017.18.HF_30))

TMM_2017.18.SERC_30 <- select(TMM_2017.18.SERC_30, -cv)
TMM_2017.18.SERC_30 <- column_to_rownames(TMM_2017.18.SERC_30,var = "Gene_ID")
TMM_2017.18.SERC_30 <- as.data.frame(t(TMM_2017.18.SERC_30))

TMM_2017.20_30 <- select(TMM_2017.20_30, -cv)
TMM_2017.20_30 <- column_to_rownames(TMM_2017.20_30,var = "Gene_ID")
TMM_2017.20_30 <- as.data.frame(t(TMM_2017.20_30))

TMM_2018.20_30 <- select(TMM_2018.20_30, -cv)
TMM_2018.20_30 <- column_to_rownames(TMM_2018.20_30,var = "Gene_ID")
TMM_2018.20_30 <- as.data.frame(t(TMM_2018.20_30))

TMM_2017.18.Spring_30 <- select(TMM_2017.18.Spring_30, -cv)
TMM_2017.18.Spring_30 <- column_to_rownames(TMM_2017.18.Spring_30,var = "Gene_ID")
TMM_2017.18.Spring_30 <- as.data.frame(t(TMM_2017.18.Spring_30))

TMM_2017.18.Summer_30 <- select(TMM_2017.18.Summer_30, -cv)
TMM_2017.18.Summer_30 <- column_to_rownames(TMM_2017.18.Summer_30,var = "Gene_ID")
TMM_2017.18.Summer_30 <- as.data.frame(t(TMM_2017.18.Summer_30))

TMM_2017.18.Fall_30 <- select(TMM_2017.18.Fall_30, -cv)
TMM_2017.18.Fall_30 <- column_to_rownames(TMM_2017.18.Fall_30,var = "Gene_ID")
TMM_2017.18.Fall_30 <- as.data.frame(t(TMM_2017.18.Fall_30))

# Continue with formatting
Sample_Description = read_csv("./Formatted.Data/FAGR.description.csv")

#### QC ####
# iterative filtering of samples and genes with too many missing entries
gsg_2017.18 = goodSamplesGenes(TMM_2017.18_30, verbose = 3);
gsg_2017.18$allOK
gsg_2018.19 = goodSamplesGenes(TMM_2018.19_30, verbose = 3);
gsg_2018.19$allOK # FALSE, some samples with 0 log2FC, 43 genes exlcluded for missing samples or zero variance
# filter bad genes
TMM_2018.19_30_2 = TMM_2018.19_30[gsg_2018.19$goodSamples, gsg_2018.19$goodGenes]
# rerun QC
gsg_2018.19 = goodSamplesGenes(TMM_2018.19_30_2, verbose = 3);
gsg_2018.19$allOK

gsg_2017.19 = goodSamplesGenes(TMM_2017.19_30, verbose = 3);
gsg_2017.19$allOK # FALSE, some samples with 0 log2FC, 92 genes excluded for missing samples or zero variance
# filter bad genes
TMM_2017.19_30_2 = TMM_2017.19_30[gsg_2017.19$goodSamples, gsg_2017.19$goodGenes]
# rerun QC
gsg_2017.19 = goodSamplesGenes(TMM_2017.19_30_2, verbose = 3);
gsg_2017.19$allOK

gsg_2019.20 = goodSamplesGenes(TMM_2019.20_30, verbose = 3);
gsg_2019.20$allOK # FALSE, some samples with 0 log2FC, 222 genes excluded for missing samples or zero variance
# filter bad genes
TMM_2019.20_30_2 = TMM_2019.20_30[gsg_2019.20$goodSamples, gsg_2019.20$goodGenes]
# rerun QC
gsg_2019.20 = goodSamplesGenes(TMM_2019.20_30_2, verbose = 3);
gsg_2019.20$allOK

gsg_2017.18.HF = goodSamplesGenes(TMM_2017.18.HF_30, verbose = 3);
gsg_2017.18.HF$allOK
gsg_2017.18.SERC = goodSamplesGenes(TMM_2017.18.SERC_30, verbose = 3);
gsg_2017.18.SERC$allOK
gsg_2017.20 = goodSamplesGenes(TMM_2017.20_30, verbose = 3);
gsg_2017.20$allOK # FALSE, some samples with 0 log2FC, 267 genes excluded for missing samples or zero variance
# filter bad genes
TMM_2017.20_30_2 = TMM_2017.20_30[gsg_2017.20$goodSamples, gsg_2017.20$goodGenes]
# rerun QC
gsg_2017.20 = goodSamplesGenes(TMM_2017.20_30_2, verbose = 3);
gsg_2017.20$allOK

gsg_2018.20 = goodSamplesGenes(TMM_2018.20_30, verbose = 3);
gsg_2018.20$allOK # FALSE, some samples with 0 log2FC, 241 genes excluded for missing samples or zero variance
# filter bad genes
TMM_2018.20_30_2 = TMM_2018.20_30[gsg_2018.20$goodSamples, gsg_2018.20$goodGenes]
# rerun QC
gsg_2018.20 = goodSamplesGenes(TMM_2018.20_30_2, verbose = 3);
gsg_2018.20$allOK

gsg_2017.18.Spring = goodSamplesGenes(TMM_2017.18.Spring_30, verbose = 3);
gsg_2017.18.Spring$allOK # FALSE, some samples with 0 log2FC, 52 genes excluded for missing samples or zero variance
# filter bad genes
TMM_2017.18.Spring_30_2 = TMM_2017.18.Spring_30[gsg_2017.18.Spring$goodSamples, gsg_2017.18.Spring$goodGenes]
# rerun QC
gsg_2017.18.Spring = goodSamplesGenes(TMM_2017.18.Spring_30_2, verbose = 3);
gsg_2017.18.Spring$allOK

gsg_2017.18.Summer = goodSamplesGenes(TMM_2017.18.Summer_30, verbose = 3);
gsg_2017.18.Summer$allOK
gsg_2017.18.Fall = goodSamplesGenes(TMM_2017.18.Fall_30, verbose = 3);
gsg_2017.18.Fall$allOK # FALSE, some samples with 0 log2FC, 100 genes excluded for missing samples or zero variance
# filter bad genes
TMM_2017.18.Fall_30_2 = TMM_2017.18.Fall_30[gsg_2017.18.Fall$goodSamples, gsg_2017.18.Fall$goodGenes]
# rerun QC
gsg_2017.18.Fall = goodSamplesGenes(TMM_2017.18.Fall_30_2, verbose = 3);
gsg_2017.18.Fall$allOK

#### Soft Threshold ####
# Choose a set of soft-thresholding powers
# value is the power to which co-expression similarity is raised to calculate adjacency
# chose the value closest to 0.90
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft_2017.18 <- pickSoftThreshold(TMM_2017.18_30, powerVector = powers, verbose = 5)
sft_2018.19 <- pickSoftThreshold(TMM_2018.19_30_2, powerVector = powers, verbose = 5)
sft_2017.19 <- pickSoftThreshold(TMM_2017.19_30_2, powerVector = powers, verbose = 5)
sft_2019.20 <- pickSoftThreshold(TMM_2019.20_30_2, powerVector = powers, verbose = 5)
sft_2017.18.HF <- pickSoftThreshold(TMM_2017.18.HF_30, powerVector = powers, verbose = 5)
sft_2017.18.SERC <- pickSoftThreshold(TMM_2017.18.SERC_30, powerVector = powers, verbose = 5)
sft_2017.20 <- pickSoftThreshold(TMM_2017.20_30_2, powerVector = powers, verbose = 5)
sft_2018.20 <- pickSoftThreshold(TMM_2018.20_30_2, powerVector = powers, verbose = 5)
sft_2017.18.Spring <- pickSoftThreshold(TMM_2017.18.Spring_30_2, powerVector = powers, verbose = 5)
sft_2017.18.Summer <- pickSoftThreshold(TMM_2017.18.Summer_30, powerVector = powers, verbose = 5)
sft_2017.18.Fall <- pickSoftThreshold(TMM_2017.18.Fall_30_2, powerVector = powers, verbose = 5)

#### Soft Threshold Plotting ####
# plotting code repeated for each dataset, changing sft to appropriate name
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
# power value was chosen as value closest to 90%

adjacency_2017.18 = adjacency(TMM_2017.18_30, power = 7,type = "signed hybrid")
adjacency_2018.19 = adjacency(TMM_2018.19_30_2, power = 5,type = "signed hybrid")
adjacency_2017.19 = adjacency(TMM_2017.19_30_2, power = 6,type = "signed hybrid")
adjacency_2019.20 = adjacency(TMM_2019.20_30_2, power = 10,type = "signed hybrid")
adjacency_2017.18.HF = adjacency(TMM_2017.18.HF_30, power = 12,type = "signed hybrid")
adjacency_2017.18.SERC = adjacency(TMM_2017.18.SERC_30, power = 7,type = "signed hybrid")
adjacency_2017.20 = adjacency(TMM_2017.20_30_2, power = 7,type = "signed hybrid")
adjacency_2018.20 = adjacency(TMM_2018.20_30_2, power = 8,type = "signed hybrid")
adjacency_2017.18.Spring = adjacency(TMM_2017.18.Spring_30_2, power = 6,type = "signed hybrid")
adjacency_2017.18.Summer = adjacency(TMM_2017.18.Summer_30, power = 6,type = "signed hybrid")
adjacency_2017.18.Fall = adjacency(TMM_2017.18.Fall_30_2, power = 7,type = "signed hybrid")

#### Topological overlap matrix (TOM) ####
# to minimize effects of noise and spurious associations, tranform adjacency into TOM, calculate dissimilarity
# need to rm adjacency matrices and TOM after use to make space
TOM_2017.18 = TOMsimilarity(adjacency_2017.18, TOMType = "signed Nowick")
dissTOM_2017.18 = 1 - TOM_2017.18
rm(TOM_2017.18)
rm(adjacency_2017.18)

TOM_2018.19 = TOMsimilarity(adjacency_2018.19, TOMType = "signed Nowick")
dissTOM_2018.19 = 1 - TOM_2018.19
rm(TOM_2018.19)
rm(adjacency_2018.19)

TOM_2017.19 = TOMsimilarity(adjacency_2017.19, TOMType = "signed Nowick")
dissTOM_2017.19 = 1 - TOM_2017.19
rm(TOM_2017.19)
rm(adjacency_2017.19)

TOM_2019.20 = TOMsimilarity(adjacency_2019.20, TOMType = "signed Nowick")
dissTOM_2019.20 = 1 - TOM_2019.20
rm(TOM_2019.20)
rm(adjacency_2019.20)

TOM_2017.18.HF = TOMsimilarity(adjacency_2017.18.HF, TOMType = "signed Nowick")
dissTOM_2017.18.HF = 1 - TOM_2017.18.HF
rm(TOM_2017.18.HF)
rm(adjacency_2017.18.HF)

TOM_2017.18.SERC = TOMsimilarity(adjacency_2017.18.SERC, TOMType = "signed Nowick")
dissTOM_2017.18.SERC = 1 - TOM_2017.18.SERC
rm(TOM_2017.18.SERC)
rm(adjacency_2017.18.SERC)

TOM_2017.20 = TOMsimilarity(adjacency_2017.20, TOMType = "signed Nowick")
dissTOM_2017.20 = 1 - TOM_2017.20
rm(TOM_2017.20)
rm(adjacency_2017.20)

TOM_2018.20 = TOMsimilarity(adjacency_2018.20, TOMType = "signed Nowick")
dissTOM_2018.20 = 1 - TOM_2018.20
rm(TOM_2018.20)
rm(adjacency_2018.20)

TOM_2017.18.Spring = TOMsimilarity(adjacency_2017.18.Spring, TOMType = "signed Nowick")
dissTOM_2017.18.Spring = 1 - TOM_2017.18.Spring
rm(TOM_2017.18.Spring)
rm(adjacency_2017.18.Spring)

TOM_2017.18.Summer = TOMsimilarity(adjacency_2017.18.Summer, TOMType = "signed Nowick")
dissTOM_2017.18.Summer = 1 - TOM_2017.18.Summer
rm(TOM_2017.18.Summer)
rm(adjacency_2017.18.Summer)

TOM_2017.18.Fall = TOMsimilarity(adjacency_2017.18.Fall, TOMType = "signed Nowick")
dissTOM_2017.18.Fall = 1 - TOM_2017.18.Fall
rm(TOM_2017.18.Fall)
rm(adjacency_2017.18.Fall)

# Subset TMM by significantly differently expressed genes or maybe all genes and take difference between samples,
# so differences in log2FC change between samples between years.

#### Gene Clustering Plot ####
# Code from this section until the end was repeated for each dissTom element.
# For each repeat of dataset change: name of dissTOM, soft power number, TPM_2017_Fall_30, ect.
# Call the hierarchical clustering function
# hclust dendrogram of genes
geneTree = hclust(as.dist(dissTOM_2017.18), method = "average")

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12, 9)
jpeg("./Plots/DE/dissTOM_2017.18-tree.jpeg",height = 1080, width = 1920)
plot(
  geneTree,
  xlab = "",
  sub = "",
  main = "Gene clustering on TOM-based dissimilarity",
  labels = FALSE,
  hang = 0.04
)
# in plot, each vertical line is a gene, branches group highly co-expressed genes.
# module identification amounts to identification of individual branches

#### Build Moduldes ####
# change module size depending on results?
# We like large modules, so we set the minimum module size relatively high:
minModuleSize <- 30

# Module identification using dynamic tree cut:
# adaptive branch pruning of hierarchical clustering dendrograms
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM_2017.18,
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
jpeg("./Plots/DE/dissTOM_2017.18-tree-blocks.jpeg",height = 1080, width = 1920)
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
# quantify co-expression similarity of entire modules to see which to merge
MEList <- moduleEigengenes(TMM_2017.18_30,
                           colors = dynamicColors,
                           softPower = 7)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(7, 6)
jpeg("./Plots/ME-tree_2017.18.jpeg",height = 1080, width = 1920)

plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")

MEDissThres <-  0.25 # threshold for merging modules, corresponds to module correlation of 0.75

# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")

# Call an automatic merging function
merge <- mergeCloseModules(TMM_2017.18_30,
                           dynamicColors,
                           cutHeight = MEDissThres,
                           verbose = 3)
# The merged module colors
mergedColors <- merge$colors

# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs
dev.off()

sizeGrWindow(12, 9)
pdf(file = "./Plots/DE/geneDendro_2017.18.pdf", wi = 9, he = 6)
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
# plot shows what merging of modules did

#### Color Renaming ####
# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))

moduleLabels <- match(moduleColors, colorOrder) - 1

MEs <- mergedMEs

#### Creating a Module Heatmap of gene Expression ####
Module_Heatmap <- TMM_2017.18_30 %>%
  rownames_to_column(var = "sample") %>%
  column_to_rownames(var = "sample")


ModuleGenes_Df <- data.frame(gene = NULL, color = NULL)
for(color in unique(moduleColors)) {
  temp <- data.frame(names(TMM_2017.18_30)[moduleColors == color])
  names(temp)[1] <- "gene"
  temp$color <- color
  ModuleGenes_Df <- rbind(ModuleGenes_Df, temp)
}

Module_Matrix <- as.matrix(Module_Heatmap)

reorder.idx <- match(colnames(Module_Matrix), ModuleGenes_Df$gene)
ModuleGenes_Df <- ModuleGenes_Df[reorder.idx, ]

my.pallete <-
  scales::div_gradient_pal(low = "red", mid = "white", high = "blue")
my.colors <- my.pallete(seq(0, 1, length.out = 9))
#my.breaks <- c(-1, -.5, .5, 1)

# play with margins here
jpeg("./Plots/DE/ModuleGene_Heatmap_2017.18.jpeg",
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

#write_csv(growth_MEs,"./Data/DE.MEs/time.2017.18.MEs.csv")

#### Models #####
# had to manual take the difference in RGR and growth between years
# added this to ME dataframes
# Can only do single season comparisons because growth would be repeated between seasons

##### Spring 2017-2018 #####
Spring.2017.18.MEs = read.csv("./Data/DE.MEs/time.2017.18.Spring.MEs.csv") %>%
  mutate(growth.change.log = log(growth.change+1),
         RGR.change.log = log(RGR.change+1)) %>%
  rename_all(~ str_replace(., "ME","ME_"))

Spring.2017.18.MEs.2 = Spring.2017.18.MEs %>%
  pivot_longer(starts_with("ME"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,eigen_value,RGR.change,growth.change,RGR.change.log,growth.change.log))

Spring.2017.18.MEs.RGR.change <- Spring.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR.change ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.18.MEs.RGR.change <- Spring.2017.18.MEs.RGR.change %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR.change)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Spring.2017.18.MEs.RGR.change %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Spring.2017.18.MEs.growth.change <- Spring.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth.change ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.18.MEs.growth.change <- Spring.2017.18.MEs.growth.change %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth.change)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Spring.2017.18.MEs.growth.change %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Spring.2017.18.MEs.RGR.change.log <- Spring.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR.change.log ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.18.MEs.RGR.change.log <- Spring.2017.18.MEs.RGR.change.log %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR.change.log)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Spring.2017.18.MEs.RGR.change.log %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Spring.2017.18.MEs.growth.change.log <- Spring.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth.change.log ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.18.MEs.growth.change.log <- Spring.2017.18.MEs.growth.change.log %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth.change.log)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Spring.2017.18.MEs.growth.change.log %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

#### Summer 2017-2018 ####
