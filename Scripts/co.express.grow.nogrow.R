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

#### Read in our TMM normalized RNA-seq expression data ####

#### TMM Data Saving ####
TMM.spring.2017 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2017.csv")
TMM.summer.2017 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2017.csv") 
TMM.fall.2017 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2017.csv") 

TMM.spring.2018 = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.spring.2018.csv")
TMM.summer.2018 = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.summer.2018.csv") 
TMM.fall.2018 = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.fall.2018.csv") 

TMM.spring.2019 = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.spring.2019.csv")
TMM.summer.2019 = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.summer.2019.csv") 
TMM.fall.2019 = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.fall.2019.csv") 

TMM.spring.2020 = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.spring.2020.csv")

TMM.spring.2017.HF = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.spring.2017.HF.csv")
TMM.summer.2017.HF = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.summer.2017.HF.csv") 
TMM.fall.2017.HF = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.fall.2017.HF.csv") 
TMM.spring.2017.SERC = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.spring.2017.SERC.csv")
TMM.summer.2017.SERC = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.summer.2017.SERC.csv") 
TMM.fall.2017.SERC = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.fall.2017.SERC.csv") 

TMM.spring.2018.HF = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.spring.2018.HF.csv")
TMM.summer.2018.HF = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.summer.2018.HF.csv") 
TMM.fall.2018.HF = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.fall.2018.HF.csv") 
TMM.spring.2018.SERC = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.spring.2018.SERC.csv")
TMM.summer.2018.SERC = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.summer.2018.SERC.csv") 
TMM.fall.2018.SERC = read.csv("./Data/grow.nogrow/TMM.NormLogCPM.fall.2018.SERC.csv")

#### Calculate the Variance of each Gene across all of our samples ####

# Calculating Coefficient of variation function
calc.cv <- function(x, na.rm=TRUE) { 
  if(na.rm==TRUE) x <- na.omit(x)
  result <- sd(x) / mean(x) 
  result <- abs(result) 
  return(result)
}

TMM.spring.2017.cv <- TMM.spring.2017 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.summer.2017.cv <- TMM.summer.2017 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.fall.2017.cv <- TMM.fall.2017 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

TMM.spring.2018.cv <- TMM.spring.2018 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.summer.2018.cv <- TMM.summer.2018 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.fall.2018.cv <- TMM.fall.2018 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

TMM.spring.2019.cv <- TMM.spring.2019 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.summer.2019.cv <- TMM.summer.2019 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.fall.2019.cv <- TMM.fall.2019 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

TMM.spring.2020.cv <- TMM.spring.2020 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

TMM.spring.2017.HF.cv <- TMM.spring.2017.HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.summer.2017.HF.cv <- TMM.summer.2017.HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.fall.2017.HF.cv <- TMM.fall.2017.HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.spring.2017.SERC.cv <- TMM.spring.2017.SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.summer.2017.SERC.cv <- TMM.summer.2017.SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.fall.2017.SERC.cv <- TMM.fall.2017.SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

TMM.spring.2018.HF.cv <- TMM.spring.2018.HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.summer.2018.HF.cv <- TMM.summer.2018.HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.fall.2018.HF.cv <- TMM.fall.2018.HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.spring.2018.SERC.cv <- TMM.spring.2018.SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.summer.2018.SERC.cv <- TMM.summer.2018.SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.fall.2018.SERC.cv <- TMM.fall.2018.SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

# Filter the data frame to contain the top 30% most variable genes)
# ties are kept together

TMM.spring.2017.30  <- TMM.spring.2017.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2017.30  <- TMM.summer.2017.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2017.30  <- TMM.fall.2017.cv  %>% slice_max(order_by = cv, prop = .30)

TMM.spring.2018.30  <- TMM.spring.2018.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2018.30  <- TMM.summer.2018.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2018.30  <- TMM.fall.2018.cv  %>% slice_max(order_by = cv, prop = .30)

TMM.spring.2019.30  <- TMM.spring.2019.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2019.30  <- TMM.summer.2019.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2019.30  <- TMM.fall.2019.cv  %>% slice_max(order_by = cv, prop = .30)

TMM.spring.2020.30  <- TMM.spring.2020.cv  %>% slice_max(order_by = cv, prop = .30)

TMM.spring.2017.HF.30  <- TMM.spring.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2017.HF.30  <- TMM.summer.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2017.HF.30  <- TMM.fall.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.spring.2017.HF.30  <- TMM.spring.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2017.HF.30  <- TMM.summer.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2017.HF.30  <- TMM.fall.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)

TMM.spring.2018.HF.30  <- TMM.spring.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2018.HF.30  <- TMM.summer.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2018.HF.30  <- TMM.fall.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.spring.2018.HF.30  <- TMM.spring.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2018.HF.30  <- TMM.summer.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2018.HF.30  <- TMM.fall.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)

# Deselect the cv column and flip our data frame to contain sample along the rows and genes along the columns
TMM.spring.2017.30 <- select(TMM.spring.2017.30, -cv)
TMM.spring.2017.30 <- column_to_rownames(TMM.spring.2017.30,var = "Gene_ID")
TMM.spring.2017.30 <- as.data.frame(t(TMM.spring.2017.30))

TMM.summer.2017.30 <- select(TMM.summer.2017.30, -cv)
TMM.summer.2017.30 <- column_to_rownames(TMM.summer.2017.30,var = "Gene_ID")
TMM.summer.2017.30 <- as.data.frame(t(TMM.summer.2017.30))

TMM.fall.2017.30 <- select(TMM.fall.2017.30, -cv)
TMM.fall.2017.30 <- column_to_rownames(TMM.fall.2017.30,var = "Gene_ID")
TMM.fall.2017.30 <- as.data.frame(t(TMM.fall.2017.30))

TMM.spring.2018.30 <- select(TMM.spring.2018.30, -cv)
TMM.spring.2018.30 <- column_to_rownames(TMM.spring.2018.30,var = "Gene_ID")
TMM.spring.2018.30 <- as.data.frame(t(TMM.spring.2018.30))

TMM.summer.2018.30 <- select(TMM.summer.2018.30, -cv)
TMM.summer.2018.30 <- column_to_rownames(TMM.summer.2018.30,var = "Gene_ID")
TMM.summer.2018.30 <- as.data.frame(t(TMM.summer.2018.30))

TMM.fall.2018.30 <- select(TMM.fall.2018.30, -cv)
TMM.fall.2018.30 <- column_to_rownames(TMM.fall.2018.30,var = "Gene_ID")
TMM.fall.2018.30 <- as.data.frame(t(TMM.fall.2018.30))

TMM.spring.2019.30 <- select(TMM.spring.2019.30, -cv)
TMM.spring.2019.30 <- column_to_rownames(TMM.spring.2019.30,var = "Gene_ID")
TMM.spring.2019.30 <- as.data.frame(t(TMM.spring.2019.30))

TMM.summer.2019.30 <- select(TMM.summer.2019.30, -cv)
TMM.summer.2019.30 <- column_to_rownames(TMM.summer.2019.30,var = "Gene_ID")
TMM.summer.2019.30 <- as.data.frame(t(TMM.summer.2019.30))

TMM.fall.2019.30 <- select(TMM.fall.2019.30, -cv)
TMM.fall.2019.30 <- column_to_rownames(TMM.fall.2019.30,var = "Gene_ID")
TMM.fall.2019.30 <- as.data.frame(t(TMM.fall.2019.30))

TMM.spring.2020.30 <- select(TMM.spring.2020.30, -cv)
TMM.spring.2020.30 <- column_to_rownames(TMM.spring.2020.30,var = "Gene_ID")
TMM.spring.2020.30 <- as.data.frame(t(TMM.spring.2020.30))

TMM.spring.2017.HF.30 <- select(TMM.spring.2017.HF.30, -cv)
TMM.spring.2017.HF.30 <- column_to_rownames(TMM.spring.2017.HF.30,var = "Gene_ID")
TMM.spring.2017.HF.30 <- as.data.frame(t(TMM.spring.2017.HF.30))

TMM.summer.2017.HF.30 <- select(TMM.summer.2017.HF.30, -cv)
TMM.summer.2017.HF.30 <- column_to_rownames(TMM.summer.2017.HF.30,var = "Gene_ID")
TMM.summer.2017.HF.30 <- as.data.frame(t(TMM.summer.2017.HF.30))

TMM.fall.2017.HF.30 <- select(TMM.fall.2017.HF.30, -cv)
TMM.fall.2017.HF.30 <- column_to_rownames(TMM.fall.2017.HF.30,var = "Gene_ID")
TMM.fall.2017.HF.30 <- as.data.frame(t(TMM.fall.2017.HF.30))

TMM.spring.2017.SERC.30 <- select(TMM.spring.2017.SERC.30, -cv)
TMM.spring.2017.SERC.30 <- column_to_rownames(TMM.spring.2017.SERC.30,var = "Gene_ID")
TMM.spring.2017.SERC.30 <- as.data.frame(t(TMM.spring.2017.SERC.30))

TMM.summer.2017.SERC.30 <- select(TMM.summer.2017.SERC.30, -cv)
TMM.summer.2017.SERC.30 <- column_to_rownames(TMM.summer.2017.SERC.30,var = "Gene_ID")
TMM.summer.2017.SERC.30 <- as.data.frame(t(TMM.summer.2017.SERC.30))

TMM.fall.2017.SERC.30 <- select(TMM.fall.2017.SERC.30, -cv)
TMM.fall.2017.SERC.30 <- column_to_rownames(TMM.fall.2017.SERC.30,var = "Gene_ID")
TMM.fall.2017.SERC.30 <- as.data.frame(t(TMM.fall.2017.SERC.30))

TMM.spring.2018.HF.30 <- select(TMM.spring.2018.HF.30, -cv)
TMM.spring.2018.HF.30 <- column_to_rownames(TMM.spring.2018.HF.30,var = "Gene_ID")
TMM.spring.2018.HF.30 <- as.data.frame(t(TMM.spring.2018.HF.30))

TMM.summer.2018.HF.30 <- select(TMM.summer.2018.HF.30, -cv)
TMM.summer.2018.HF.30 <- column_to_rownames(TMM.summer.2018.HF.30,var = "Gene_ID")
TMM.summer.2018.HF.30 <- as.data.frame(t(TMM.summer.2018.HF.30))

TMM.fall.2018.HF.30 <- select(TMM.fall.2018.HF.30, -cv)
TMM.fall.2018.HF.30 <- column_to_rownames(TMM.fall.2018.HF.30,var = "Gene_ID")
TMM.fall.2018.HF.30 <- as.data.frame(t(TMM.fall.2018.HF.30))

TMM.spring.2018.SERC.30 <- select(TMM.spring.2018.SERC.30, -cv)
TMM.spring.2018.SERC.30 <- column_to_rownames(TMM.spring.2018.SERC.30,var = "Gene_ID")
TMM.spring.2018.SERC.30 <- as.data.frame(t(TMM.spring.2018.SERC.30))

TMM.summer.2018.SERC.30 <- select(TMM.summer.2018.SERC.30, -cv)
TMM.summer.2018.SERC.30 <- column_to_rownames(TMM.summer.2018.SERC.30,var = "Gene_ID")
TMM.summer.2018.SERC.30 <- as.data.frame(t(TMM.summer.2018.SERC.30))

TMM.fall.2018.SERC.30 <- select(TMM.fall.2018.SERC.30, -cv)
TMM.fall.2018.SERC.30 <- column_to_rownames(TMM.fall.2018.SERC.30,var = "Gene_ID")
TMM.fall.2018.SERC.30 <- as.data.frame(t(TMM.fall.2018.SERC.30))

#### QC ####
# iterative filtering of samples and genes with too many missing entries
gsg.spring.2017 = goodSamplesGenes(TMM.spring.2017.30, verbose = 3);
gsg.spring.2017$allOK
gsg.summer.2017 = goodSamplesGenes(TMM.summer.2017.30, verbose = 3);
gsg.summer.2017$allOK
gsg.fall.2017 = goodSamplesGenes(TMM.fall.2017.30, verbose = 3);
gsg.fall.2017$allOK
gsg.spring.2018 = goodSamplesGenes(TMM.spring.2018.30, verbose = 3);
gsg.spring.2018$allOK
gsg.summer.2018 = goodSamplesGenes(TMM.summer.2018.30, verbose = 3);
gsg.summer.2018$allOK
gsg.fall.2018 = goodSamplesGenes(TMM.fall.2018.30, verbose = 3);
gsg.fall.2018$allOK
gsg.spring.2019 = goodSamplesGenes(TMM.spring.2019.30, verbose = 3);
gsg.spring.2019$allOK
gsg.summer.2019 = goodSamplesGenes(TMM.summer.2019.30, verbose = 3);
gsg.summer.2019$allOK
gsg.fall.2019 = goodSamplesGenes(TMM.fall.2019.30, verbose = 3);
gsg.fall.2019$allOK
gsg.spring.2020 = goodSamplesGenes(TMM.spring.2020.30, verbose = 3);
gsg.spring.2020$allOK
gsg.spring.2017.HF = goodSamplesGenes(TMM.spring.2017.HF.30, verbose = 3);
gsg.spring.2017.HF$allOK
gsg.summer.2017.HF = goodSamplesGenes(TMM.summer.2017.HF.30, verbose = 3);
gsg.summer.2017.HF$allOK
gsg.fall.2017.HF = goodSamplesGenes(TMM.fall.2017.HF.30, verbose = 3);
gsg.fall.2017.HF$allOK
gsg.spring.2017.SERC = goodSamplesGenes(TMM.spring.2017.SERC.30, verbose = 3);
gsg.spring.2017.SERC$allOK
gsg.summer.2017.SERC = goodSamplesGenes(TMM.summer.2017.SERC.30, verbose = 3);
gsg.summer.2017.SERC$allOK
gsg.fall.2017.SERC = goodSamplesGenes(TMM.fall.2017.SERC.30, verbose = 3);
gsg.fall.2017.SERC$allOK
gsg.spring.2018.HF = goodSamplesGenes(TMM.spring.2018.HF.30, verbose = 3);
gsg.spring.2018.HF$allOK
gsg.summer.2018.HF = goodSamplesGenes(TMM.summer.2018.HF.30, verbose = 3);
gsg.summer.2018.HF$allOK
gsg.fall.2018.HF = goodSamplesGenes(TMM.fall.2018.HF.30, verbose = 3);
gsg.fall.2018.HF$allOK
gsg.spring.2018.SERC = goodSamplesGenes(TMM.spring.2018.SERC.30, verbose = 3);
gsg.spring.2018.SERC$allOK
gsg.summer.2018.SERC = goodSamplesGenes(TMM.summer.2018.SERC.30, verbose = 3);
gsg.summer.2018.SERC$allOK
gsg.fall.2018.SERC = goodSamplesGenes(TMM.fall.2018.SERC.30, verbose = 3);
gsg.fall.2018.SERC$allOK

#### Soft Threshold ####
# Choose a set of soft-thresholding powers
# value is the power to which co-expression similarity is raised to calculate adjacency
# chose the value closest to 0.90
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft.spring.2017 <- pickSoftThreshold(TMM.spring.2017.30, powerVector = powers, verbose = 5)
sft.summer.2017 <- pickSoftThreshold(TMM.summer.2017.30, powerVector = powers, verbose = 5)
sft.fall.2017 <- pickSoftThreshold(TMM.fall.2017.30, powerVector = powers, verbose = 5)
sft.spring.2018 <- pickSoftThreshold(TMM.spring.2018.30, powerVector = powers, verbose = 5)
sft.summer.2018 <- pickSoftThreshold(TMM.summer.2018.30, powerVector = powers, verbose = 5)
sft.fall.2018 <- pickSoftThreshold(TMM.fall.2018.30, powerVector = powers, verbose = 5)
sft.spring.2019 <- pickSoftThreshold(TMM.spring.2019.30, powerVector = powers, verbose = 5)
sft.summer.2019 <- pickSoftThreshold(TMM.summer.2019.30, powerVector = powers, verbose = 5)
sft.fall.2019 <- pickSoftThreshold(TMM.fall.2019.30, powerVector = powers, verbose = 5)
sft.spring.2020 <- pickSoftThreshold(TMM.spring.2020.30, powerVector = powers, verbose = 5)
sft.spring.2017.HF <- pickSoftThreshold(TMM.spring.2017.HF.30, powerVector = powers, verbose = 5)
sft.summer.2017.HF <- pickSoftThreshold(TMM.summer.2017.HF.30, powerVector = powers, verbose = 5)
sft.fall.2017.HF <- pickSoftThreshold(TMM.fall.2017.HF.30, powerVector = powers, verbose = 5)
sft.spring.2017.SERC <- pickSoftThreshold(TMM.spring.2017.SERC.30, powerVector = powers, verbose = 5)
sft.summer.2017.SERC <- pickSoftThreshold(TMM.summer.2017.SERC.30, powerVector = powers, verbose = 5)
sft.fall.2017.SERC <- pickSoftThreshold(TMM.fall.2017.SERC.30, powerVector = powers, verbose = 5)
sft.spring.2018.HF <- pickSoftThreshold(TMM.spring.2018.HF.30, powerVector = powers, verbose = 5)
sft.summer.2018.HF <- pickSoftThreshold(TMM.summer.2018.HF.30, powerVector = powers, verbose = 5)
sft.fall.2018.HF <- pickSoftThreshold(TMM.fall.2018.HF.30, powerVector = powers, verbose = 5)
sft.spring.2018.SERC <- pickSoftThreshold(TMM.spring.2018.SERC.30, powerVector = powers, verbose = 5)
sft.summer.2018.SERC <- pickSoftThreshold(TMM.summer.2018.SERC.30, powerVector = powers, verbose = 5)
sft.fall.2018.SERC <- pickSoftThreshold(TMM.fall.2018.SERC.30, powerVector = powers, verbose = 5)

#### Soft Threshold Plotting ####
# plotting code repeated for each dataset, changing sft to appropriate name

# Plot the results:
jpeg("./Plots/grow.nogrow/soft.thresholding.spring.2017.jpeg",height = 1080, width = 1920)
par(mfrow = c(1, 2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power

plot(
  sft.spring.2017$fitIndices[, 1],
  -sign(sft.spring.2017$fitIndices[, 3]) * sft.spring.2017$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit,signed R^2",
  type = "n",
  main = paste("Scale independence")
)

text(
  sft.spring.2017$fitIndices[, 1],
  -sign(sft.spring.2017$fitIndices[, 3]) * sft.spring.2017$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)

# this line corresponds to using an R^2 cut-off of h
abline(h = 0.90, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(
  sft.spring.2017$fitIndices[, 1],
  sft.spring.2017$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  sft.spring.2017$fitIndices[, 1],
  sft.spring.2017$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)
dev.off()

#### Adjacency ####
# signed hybrid sets all negatively correlated genes to 0
# choice of power comes from soft threshold analysis above
# power value was chosen as value closest to 90%

adjacency.spring.2017 = adjacency(TMM.spring.2017.30, power = 5,type = "signed hybrid")
adjacency.summer.2017 = adjacency(TMM.summer.2017.30, power = 5,type = "signed hybrid")
adjacency.fall.2017 = adjacency(TMM.fall.2017.30, power = 5,type = "signed hybrid")
adjacency.spring.2018 = adjacency(TMM.spring.2018.30, power = 5,type = "signed hybrid")
adjacency.summer.2018 = adjacency(TMM.summer.2018.30, power = 5,type = "signed hybrid")
adjacency.fall.2018 = adjacency(TMM.fall.2018.30, power = 5,type = "signed hybrid")
adjacency.spring.2019 = adjacency(TMM.spring.2019.30, power = 5,type = "signed hybrid")
adjacency.summer.2019 = adjacency(TMM.summer.2019.30, power = 5,type = "signed hybrid")
adjacency.fall.2019 = adjacency(TMM.fall.2019.30, power = 5,type = "signed hybrid")
adjacency.spring.2020 = adjacency(TMM.spring.2020.30, power = 5,type = "signed hybrid")
adjacency.spring.2017.HF = adjacency(TMM.spring.2017.HF.30, power = 5,type = "signed hybrid")
adjacency.summer.2017.HF = adjacency(TMM.summer.2017.HF.30, power = 5,type = "signed hybrid")
adjacency.fall.2017.HF = adjacency(TMM.fall.2017.HF.30, power = 5,type = "signed hybrid")
adjacency.spring.2017.SERC = adjacency(TMM.spring.2017.SERC.30, power = 5,type = "signed hybrid")
adjacency.summer.2017.SERC = adjacency(TMM.summer.2017.SERC.30, power = 5,type = "signed hybrid")
adjacency.fall.2017.SERC = adjacency(TMM.fall.2017.SERC.30, power = 5,type = "signed hybrid")
adjacency.spring.2018.HF = adjacency(TMM.spring.2018.HF.30, power = 5,type = "signed hybrid")
adjacency.summer.2018.HF = adjacency(TMM.summer.2018.HF.30, power = 5,type = "signed hybrid")
adjacency.fall.2018.HF = adjacency(TMM.fall.2018.HF.30, power = 5,type = "signed hybrid")
adjacency.spring.2018.SERC = adjacency(TMM.spring.2018.SERC.30, power = 5,type = "signed hybrid")
adjacency.summer.2018.SERC = adjacency(TMM.summer.2018.SERC.30, power = 5,type = "signed hybrid")
adjacency.fall.2018.SERC = adjacency(TMM.fall.2018.SERC.30, power = 5,type = "signed hybrid")

#### Topological overlap matrix (TOM) ####
# to minimize effects of noise and spurious associations, transform adjacency into TOM, calculate dissimilarity
# need to rm adjacency matrices and TOM after use to make space

TOM.spring.2017 = TOMsimilarity(adjacency.spring.2017, TOMType = "signed Nowick")
dissTOM.spring.2017 = 1 - TOM.spring.2017
rm(TOM.spring.2017)
rm(adjacency.spring.2017)
TOM.summer.2017 = TOMsimilarity(adjacency.summer.2017, TOMType = "signed Nowick")
dissTOM.summer.2017 = 1 - TOM.summer.2017
rm(TOM.summer.2017)
rm(adjacency.summer.2017)
TOM.fall.2017 = TOMsimilarity(adjacency.fall.2017, TOMType = "signed Nowick")
dissTOM.fall.2017 = 1 - TOM.fall.2017
rm(TOM.fall.2017)
rm(adjacency.fall.2017)

TOM.spring.2018 = TOMsimilarity(adjacency.spring.2018, TOMType = "signed Nowick")
dissTOM.spring.2018 = 1 - TOM.spring.2018
rm(TOM.spring.2018)
rm(adjacency.spring.2018)
TOM.summer.2018 = TOMsimilarity(adjacency.summer.2018, TOMType = "signed Nowick")
dissTOM.summer.2018 = 1 - TOM.summer.2018
rm(TOM.summer.2018)
rm(adjacency.summer.2018)
TOM.fall.2018 = TOMsimilarity(adjacency.fall.2018, TOMType = "signed Nowick")
dissTOM.fall.2018 = 1 - TOM.fall.2018
rm(TOM.fall.2018)
rm(adjacency.fall.2018)

TOM.spring.2019 = TOMsimilarity(adjacency.spring.2019, TOMType = "signed Nowick")
dissTOM.spring.2019 = 1 - TOM.spring.2019
rm(TOM.spring.2019)
rm(adjacency.spring.2019)
TOM.summer.2019 = TOMsimilarity(adjacency.summer.2019, TOMType = "signed Nowick")
dissTOM.summer.2019 = 1 - TOM.summer.2019
rm(TOM.summer.2019)
rm(adjacency.summer.2019)
TOM.fall.2019 = TOMsimilarity(adjacency.fall.2019, TOMType = "signed Nowick")
dissTOM.fall.2019 = 1 - TOM.fall.2019
rm(TOM.fall.2019)
rm(adjacency.fall.2019)

TOM.spring.2020 = TOMsimilarity(adjacency.spring.2020, TOMType = "signed Nowick")
dissTOM.spring.2020 = 1 - TOM.spring.2020
rm(TOM.spring.2020)
rm(adjacency.spring.2020)

TOM.spring.2017.HF = TOMsimilarity(adjacency.spring.2017.HF, TOMType = "signed Nowick")
dissTOM.spring.2017.HF = 1 - TOM.spring.2017.HF
rm(TOM.spring.2017.HF)
rm(adjacency.spring.2017.HF)
TOM.summer.2017.HF = TOMsimilarity(adjacency.summer.2017.HF, TOMType = "signed Nowick")
dissTOM.summer.2017.HF = 1 - TOM.summer.2017.HF
rm(TOM.summer.2017.HF)
rm(adjacency.summer.2017.HF)
TOM.fall.2017.HF = TOMsimilarity(adjacency.fall.2017.HF, TOMType = "signed Nowick")
dissTOM.fall.2017.HF = 1 - TOM.fall.2017.HF
rm(TOM.fall.2017.HF)
rm(adjacency.fall.2017.HF)
TOM.spring.2017.SERC = TOMsimilarity(adjacency.spring.2017.SERC, TOMType = "signed Nowick")
dissTOM.spring.2017.SERC = 1 - TOM.spring.2017.SERC
rm(TOM.spring.2017.SERC)
rm(adjacency.spring.2017.SERC)
TOM.summer.2017.SERC = TOMsimilarity(adjacency.summer.2017.SERC, TOMType = "signed Nowick")
dissTOM.summer.2017.SERC = 1 - TOM.summer.2017.SERC
rm(TOM.summer.2017.SERC)
rm(adjacency.summer.2017.SERC)
TOM.fall.2017.SERC = TOMsimilarity(adjacency.fall.2017.SERC, TOMType = "signed Nowick")
dissTOM.fall.2017.SERC = 1 - TOM.fall.2017.SERC
rm(TOM.fall.2017.SERC)
rm(adjacency.fall.2017.SERC)

TOM.spring.2018.HF = TOMsimilarity(adjacency.spring.2018.HF, TOMType = "signed Nowick")
dissTOM.spring.2018.HF = 1 - TOM.spring.2018.HF
rm(TOM.spring.2018.HF)
rm(adjacency.spring.2018.HF)
TOM.summer.2018.HF = TOMsimilarity(adjacency.summer.2018.HF, TOMType = "signed Nowick")
dissTOM.summer.2018.HF = 1 - TOM.summer.2018.HF
rm(TOM.summer.2018.HF)
rm(adjacency.summer.2018.HF)
TOM.fall.2018.HF = TOMsimilarity(adjacency.fall.2018.HF, TOMType = "signed Nowick")
dissTOM.fall.2018.HF = 1 - TOM.fall.2018.HF
rm(TOM.fall.2018.HF)
rm(adjacency.fall.2018.HF)
TOM.spring.2018.SERC = TOMsimilarity(adjacency.spring.2018.SERC, TOMType = "signed Nowick")
dissTOM.spring.2018.SERC = 1 - TOM.spring.2018.SERC
rm(TOM.spring.2018.SERC)
rm(adjacency.spring.2018.SERC)
TOM.summer.2018.SERC = TOMsimilarity(adjacency.summer.2018.SERC, TOMType = "signed Nowick")
dissTOM.summer.2018.SERC = 1 - TOM.summer.2018.SERC
rm(TOM.summer.2018.SERC)
rm(adjacency.summer.2018.SERC)
TOM.fall.2018.SERC = TOMsimilarity(adjacency.fall.2018.SERC, TOMType = "signed Nowick")
dissTOM.fall.2018.SERC = 1 - TOM.fall.2018.SERC
rm(TOM.fall.2018.SERC)
rm(adjacency.fall.2018.SERC)


#### Gene Clustering Plot ####
# Code from this section until the end was repeated for each dissTom element.
# For each repeat of dataset change: name of dissTOM, soft power number, TPM_2017_Fall_30, ect.

# Call the hierarchical clustering function
# hclust dendrogram of genes
geneTree = hclust(as.dist(dissTOM.spring.2017), method = "average")

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12, 9)
jpeg("./Plots/grow.nogrow/dissTOM.spring.2017-tree.jpeg",height = 1080, width = 1920)
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

# We like large modules, so we set the minimum module size relatively high:
minModuleSize <- 30

# Module identification using dynamic tree cut:
# adaptive branch pruning of hierarchical clustering dendrograms
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM.spring.2017,
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
jpeg("./Plots/grow.nogrow/dissTOM.spring.2017-tree-blocks.jpeg",height = 1080, width = 1920)
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

MEList <- moduleEigengenes(TMM.spring.2017.30,
                           colors = dynamicColors,
                           softPower = 5)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(7, 6)
jpeg("./Plots/grow.nogrow/ME-tree.spring.2017.jpeg",height = 1080, width = 1920)

plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")

MEDissThres <-  0.25 # threshold for merging modules, corresponds to module correlation of 0.75

# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")

# Call an automatic merging function
merge <- mergeCloseModules(TMM.spring.2017.30,
                           dynamicColors,
                           cutHeight = MEDissThres,
                           verbose = 3)
# The merged module colors
mergedColors <- merge$colors

# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs
dev.off()

sizeGrWindow(12, 9)
pdf(file = "./Plots/grow.nogrow/geneDendro.spring.2017.pdf", wi = 9, he = 6)
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
Module_Heatmap <- TMM.spring.2017.30 %>%
  rownames_to_column(var = "sample") %>%
  column_to_rownames(var = "sample")


ModuleGenes_Df <- data.frame(gene = NULL, color = NULL)
for(color in unique(moduleColors)) {
  temp <- data.frame(names(TMM.spring.2017.30)[moduleColors == color])
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
jpeg("./Plots/grow.nogrow/ModuleGene.Heatmap.spring.2017.jpeg",
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
growth.MEs <- MEs %>%
  rownames_to_column(var = "sample") %>%
  as_tibble()

# Bring in growth data

growth.data = read.csv("./Formatted.Data/all.growth.data.csv") %>%
  filter(YEAR == 2017) %>%
  select(SITE, YEAR, TREE_ID, RGR, growth)
#growth.data.2 = subset(growth.data, growth.data$RGR != "NA")

# bring in Sample.Description

Sample.Description = read_csv("./Formatted.Data/FAGR.description.csv")

sample.sub = subset(Sample.Description, Sample.Description$sample.description %in% growth.MEs$sample)
growth.data.3 = left_join(growth.data.2,sample_sub, by = c("TREE_ID" = "Tree_ID"))

#### Exporting Eigengene values ####
growth.MEs = left_join(growth.MEs, growth.data.3, by = c("sample" = "sample.description"))

growth.MEs <- growth.MEs %>% 
  rename_all(~ str_replace(., "ME","ME_"))
# write_csv(growth.MEs,"./Data/grow.nogrow/MEs/spring.2017.MEs.csv")

#### Models ####

spring.2017.MEs = read.csv("./Data/grow.nogrow/MEs/spring.2017.MEs.csv") %>%
  mutate(growth.state.binary = ifelse(grep("High", growth_state), 1, 0)) # binarize categorical variable

spring.2017.MEs.2 = spring.2017.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,growth.state.binary))

Spring.2017.MEs.3 <- Spring.2017.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ growth.state.binary, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.MEs.3 <- Spring.2017.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = growth.state.binary, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2017.MEs.3 <- Spring.2017.MEs.3 %>%
  filter(term == "growth.state.binary0") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# Correlation

nSamples = nrow(spring.2017.MEs) # define number of samples
nGenes = ncol(spring.2017.MEs) # define number of genes

# correlation between module eigenegens and growth states
module.trait.corr <- cor(spring.2017.MEs[,c()], spring.2017.MEs$growth.state.binary, use = 'p') # pearson correlation
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples) # p-values for correlations

# heat map of correlation
heatmap.data = merge(spring.2017.MEs[,c()], spring.2017.MEs$growth.state.binary, by  = row.names) # combining data into one dataframe
# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[19:23], y= names(heatmap.data)[1:18],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *
# extract genes from modules with significance with high growth











