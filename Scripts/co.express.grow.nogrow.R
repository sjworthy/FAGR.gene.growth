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
devtools::install_github("kevinblighe/CorLevelPlot")
library(CorLevelPlot)

#### Read in our TMM normalized RNA-seq expression data ####

TMM.spring.2017 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2017.csv")
TMM.summer.2017 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2017.csv") 
TMM.fall.2017 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2017.csv") 

TMM.spring.2018 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2018.csv")
TMM.summer.2018 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2018.csv") 
TMM.fall.2018 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2018.csv") 

TMM.spring.2019 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2019.csv")
TMM.summer.2019 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2019.csv") 
TMM.fall.2019 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2019.csv") 

#TMM.spring.2020 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2020.csv")

TMM.spring.2017.HF = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2017.HF.csv")
TMM.summer.2017.HF = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2017.HF.csv") 
TMM.fall.2017.HF = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2017.HF.csv") 
TMM.spring.2017.SERC = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2017.SERC.csv")
TMM.summer.2017.SERC = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2017.SERC.csv") 
TMM.fall.2017.SERC = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2017.SERC.csv") 

TMM.spring.2018.HF = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2018.HF.csv")
TMM.summer.2018.HF = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2018.HF.csv") 
TMM.fall.2018.HF = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2018.HF.csv") 
TMM.spring.2018.SERC = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2018.SERC.csv")
TMM.summer.2018.SERC = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2018.SERC.csv") 
TMM.fall.2018.SERC = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2018.SERC.csv")

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

# TMM.spring.2020.cv <- TMM.spring.2020 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

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

#TMM.spring.2020.30  <- TMM.spring.2020.cv  %>% slice_max(order_by = cv, prop = .30)

TMM.spring.2017.HF.30  <- TMM.spring.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2017.HF.30  <- TMM.summer.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2017.HF.30  <- TMM.fall.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.spring.2017.HF.30  <- TMM.spring.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2017.HF.30  <- TMM.summer.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2017.HF.30  <- TMM.fall.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.spring.2017.SERC.30  <- TMM.spring.2017.SERC.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2017.SERC.30  <- TMM.summer.2017.SERC.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2017.SERC.30  <- TMM.fall.2017.SERC.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.spring.2017.SERC.30  <- TMM.spring.2017.SERC.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2017.SERC.30  <- TMM.summer.2017.SERC.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2017.SERC.30  <- TMM.fall.2017.SERC.cv  %>% slice_max(order_by = cv, prop = .30)

TMM.spring.2018.HF.30  <- TMM.spring.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2018.HF.30  <- TMM.summer.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2018.HF.30  <- TMM.fall.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.spring.2018.HF.30  <- TMM.spring.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2018.HF.30  <- TMM.summer.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2018.HF.30  <- TMM.fall.2018.HF.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.spring.2018.SERC.30  <- TMM.spring.2018.SERC.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2018.SERC.30  <- TMM.summer.2018.SERC.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2018.SERC.30  <- TMM.fall.2018.SERC.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.spring.2018.SERC.30  <- TMM.spring.2018.SERC.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2018.SERC.30  <- TMM.summer.2018.SERC.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.fall.2018.SERC.30  <- TMM.fall.2018.SERC.cv  %>% slice_max(order_by = cv, prop = .30)

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

#TMM.spring.2020.30 <- select(TMM.spring.2020.30, -cv)
#TMM.spring.2020.30 <- column_to_rownames(TMM.spring.2020.30,var = "Gene_ID")
#TMM.spring.2020.30 <- as.data.frame(t(TMM.spring.2020.30))

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
#sft.spring.2020 <- pickSoftThreshold(TMM.spring.2020.30, powerVector = powers, verbose = 5)
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

adjacency.spring.2017 = adjacency(TMM.spring.2017.30, power = 6,type = "signed hybrid")
adjacency.summer.2017 = adjacency(TMM.summer.2017.30, power = 5,type = "signed hybrid")
adjacency.fall.2017 = adjacency(TMM.fall.2017.30, power = 6,type = "signed hybrid")
adjacency.spring.2018 = adjacency(TMM.spring.2018.30, power = 4,type = "signed hybrid")
adjacency.summer.2018 = adjacency(TMM.summer.2018.30, power = 8,type = "signed hybrid")
adjacency.fall.2018 = adjacency(TMM.fall.2018.30, power = 8,type = "signed hybrid")
adjacency.spring.2019 = adjacency(TMM.spring.2019.30, power = 10,type = "signed hybrid")
adjacency.summer.2019 = adjacency(TMM.summer.2019.30, power = 8,type = "signed hybrid")
adjacency.fall.2019 = adjacency(TMM.fall.2019.30, power = 9,type = "signed hybrid")
#adjacency.spring.2020 = adjacency(TMM.spring.2020.30, power = 5,type = "signed hybrid")
adjacency.spring.2017.HF = adjacency(TMM.spring.2017.HF.30, power = 8,type = "signed hybrid")
adjacency.summer.2017.HF = adjacency(TMM.summer.2017.HF.30, power = 7,type = "signed hybrid")
adjacency.fall.2017.HF = adjacency(TMM.fall.2017.HF.30, power = 10,type = "signed hybrid")
adjacency.spring.2017.SERC = adjacency(TMM.spring.2017.SERC.30, power = 9,type = "signed hybrid")
adjacency.summer.2017.SERC = adjacency(TMM.summer.2017.SERC.30, power = 9,type = "signed hybrid")
adjacency.fall.2017.SERC = adjacency(TMM.fall.2017.SERC.30, power = 7,type = "signed hybrid")
adjacency.spring.2018.HF = adjacency(TMM.spring.2018.HF.30, power = 10,type = "signed hybrid")
adjacency.summer.2018.HF = adjacency(TMM.summer.2018.HF.30, power = 6,type = "signed hybrid")
adjacency.fall.2018.HF = adjacency(TMM.fall.2018.HF.30, power = 5,type = "signed hybrid")
adjacency.spring.2018.SERC = adjacency(TMM.spring.2018.SERC.30, power = 6,type = "signed hybrid")
adjacency.summer.2018.SERC = adjacency(TMM.summer.2018.SERC.30, power = 9,type = "signed hybrid")
adjacency.fall.2018.SERC = adjacency(TMM.fall.2018.SERC.30, power = 6,type = "signed hybrid")

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

#TOM.spring.2020 = TOMsimilarity(adjacency.spring.2020, TOMType = "signed Nowick")
#dissTOM.spring.2020 = 1 - TOM.spring.2020
#rm(TOM.spring.2020)
#rm(adjacency.spring.2020)

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

#write.csv(ModuleGenes_Df, file = "./Data/module.genes/summer.2017.module.genes.csv")

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

# Bring in growth parameters

growth.data = read.csv("./Formatted.Data/Parameter_table_FAGUS-2.csv") %>%
  filter(YEAR == 2017,
         BAND_NUM == 1)

# bring in Sample.Description

Sample.Description = read_csv("./Formatted.Data/FAGR.description.csv")

# subset growth data to match samples we have

growth.data.2 = subset(growth.data, growth.data$TREE_ID %in% Sample.Description$Tree_ID)

# subset to samples in growth.MEs
sample.sub = subset(Sample.Description, Sample.Description$sample.description %in% growth.MEs$sample)

growth.data.3 = left_join(growth.data.2,sample.sub, by = c("TREE_ID" = "Tree_ID"))

#### Exporting Eigengene values ####
growth.MEs = left_join(growth.MEs, growth.data.3, by = c("sample" = "sample.description"))

growth.MEs <- growth.MEs %>% 
  rename_all(~ str_replace(., "ME","ME_"))
# write_csv(growth.MEs,"./Data/grow.nogrow.MEs/spring.2017.MEs.csv")

#### Models ####
##### spring.2017 #####
spring.2017.MEs = read.csv("./Data/grow.nogrow.MEs/spring.2017.MEs.csv")
spring.2017.MEs.2 = spring.2017.MEs[,c(1:29,45,60,66,68:70)]

spring.2017.MEs.3 = spring.2017.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.SPAD,Average.LA))

###### growth signal ######
spring.2017.MEs.4 <- spring.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site as variable
spring.2017.MEs.4 <- spring.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR ######

spring.2017.MEs.4 <- spring.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

spring.2017.MEs.4 <- spring.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######

spring.2017.MEs.4 <- spring.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

spring.2017.MEs.4 <- spring.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


###### LA ######

spring.2017.MEs.4 <- spring.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

spring.2017.MEs.4 <- spring.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


###### SPAD ######

spring.2017.MEs.4 <- spring.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

spring.2017.MEs.4 <- spring.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.MEs.4 <- spring.2017.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


###### turn growth into quantiles ######

traits = spring.2017.MEs[,c(45,60,68:70)]
rownames(traits) = spring.2017.MEs$sample
traits[is.na(traits)] <- 0
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)


###### Correlation ######

nSamples = nrow(spring.2017.MEs) # define number of samples
nGenes = ncol(spring.2017.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = spring.2017.MEs[,c(2:26)]
row.names(module_eigengenes) = spring.2017.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[26:32], y= names(heatmap.data)[1:25],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *
# extract genes from modules with significance with high growth

##### summer.2017 #####
summer.2017.MEs = read.csv("./Data/grow.nogrow.MEs/summer.2017.MEs.csv")
summer.2017.MEs.2 = summer.2017.MEs[,c(1:35,51,66,72,74:76)]
summer.2017.MEs.2[is.na(summer.2017.MEs.2)] <- 0

summer.2017.MEs.3 = summer.2017.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### Growth Signal #####
summer.2017.MEs.4 <- summer.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable
summer.2017.MEs.4 <- summer.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR ######
summer.2017.MEs.4 <- summer.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value, color = SITE)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# two modules, turquoise and lightyellow

sig.mods = summer.2017.MEs.4[c(19,20,47,48),]
sig.mods$Plot[[2]]
sig.mods$Plot[[4]]

# with site variable

summer.2017.MEs.4 <- summer.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######
summer.2017.MEs.4 <- summer.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value, color = SITE)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# three modules: turquoise, darkgrey, green

sig.mods = summer.2017.MEs.4[c(19,20,21,22,51,52),]
sig.mods$Plot[[2]]
sig.mods$Plot[[4]]
sig.mods$Plot[[6]]

# with site variable

summer.2017.MEs.4 <- summer.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


###### LA ######
summer.2017.MEs.4 <- summer.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable
summer.2017.MEs.4 <- summer.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD ######
summer.2017.MEs.4 <- summer.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable
summer.2017.MEs.4 <- summer.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.MEs.4 <- summer.2017.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


###### turn growth into quantiles ######

traits = summer.2017.MEs[,c(51,66,74:76)]
rownames(traits) = summer.2017.MEs$sample
traits[is.na(traits)] <- 0
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(summer.2017.MEs) # define number of samples
nGenes = ncol(summer.2017.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = summer.2017.MEs[,c(2:32)]
row.names(module_eigengenes) = summer.2017.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

module.trait.corr.pvals.3 = module.trait.corr.pvals.2 %>%
  filter(fdr.GR < 0.01 | fdr.SLA < 0.01)

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[32:38], y= names(heatmap.data)[1:31],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

# extract genes from modules with significance with high growth
# turquoise is low growth
# lightyellow is high growth

##### fall.2017 #####
fall.2017.MEs = read.csv("./Data/grow.nogrow.MEs/fall.2017.MEs.csv")
fall.2017.MEs.2 = fall.2017.MEs[,c(1:36,37,38,39,55,70,76,78:80)]
fall.2017.MEs.2[is.na(fall.2017.MEs.2)] <- 0

fall.2017.MEs.3 = fall.2017.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### Growth Signal #####
fall.2017.MEs.4 <- fall.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site predictor
fall.2017.MEs.4 <- fall.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


####### GR ######
fall.2017.MEs.4 <- fall.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site predictor

fall.2017.MEs.4 <- fall.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA #####
fall.2017.MEs.4 <- fall.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value, color = SITE)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# two modules: cyan, midnightblue

sig.mods = fall.2017.MEs.4[c(1,2,25,26),]
sig.mods$Plot[[2]]
sig.mods$Plot[[4]]

# with site predictor

fall.2017.MEs.4 <- fall.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value, color = SITE)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA #####
fall.2017.MEs.4 <- fall.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site predictor
fall.2017.MEs.4 <- fall.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


###### SPAD #####
fall.2017.MEs.4 <- fall.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site predictor
fall.2017.MEs.4 <- fall.2017.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.MEs.4 <- fall.2017.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


###### turn growth into quantiles #####

traits = fall.2017.MEs[,c(55,70,78:80)]
rownames(traits) = fall.2017.MEs$sample
traits[is.na(traits)] <- 0
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(fall.2017.MEs) # define number of samples
nGenes = ncol(fall.2017.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = fall.2017.MEs[,c(2:36)]
row.names(module_eigengenes) = fall.2017.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[36:42], y= names(heatmap.data)[1:35],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### spring.2018 #####
spring.2018.MEs = read.csv("./Data/grow.nogrow.MEs/spring.2018.MEs.csv")
# remove second S181Fg04 associated with band 2
spring.2018.MEs = spring.2018.MEs[c(1:16,18:34),]

spring.2018.MEs.2 = spring.2018.MEs[,c(1:21,22,23,24,40,55,61,63:65)]

spring.2018.MEs.3 = spring.2018.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### Growth signal #####
spring.2018.MEs.4 <- spring.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

spring.2018.MEs.4 <- spring.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR #####
spring.2018.MEs.4 <- spring.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

spring.2018.MEs.4 <- spring.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


###### SLA #####
spring.2018.MEs.4 <- spring.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site predictor

spring.2018.MEs.4 <- spring.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA #####
spring.2018.MEs.4 <- spring.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

spring.2018.MEs.4 <- spring.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD #####
spring.2018.MEs.4 <- spring.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

spring.2018.MEs.4 <- spring.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.MEs.4 <- spring.2018.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### turn growth into quantiles #####

traits = spring.2018.MEs[,c(40,55,63:65)]
rownames(traits) = spring.2018.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation #####

nSamples = nrow(spring.2018.MEs) # define number of samples
nGenes = ncol(spring.2018.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = spring.2018.MEs[,c(2:21)]
row.names(module_eigengenes) = spring.2018.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA= p.adjust(Average.SLA, method = "fdr"),
         fdr.LA= p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD= p.adjust(Average.SPAD, method = "fdr"))

# fdr > 0.05 for high.grow quantile purple module

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[21:27], y= names(heatmap.data)[1:20],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *
# extract genes from modules with significance with high growth

##### summer.2018 #####
summer.2018.MEs = read.csv("./Data/grow.nogrow.MEs/summer.2018.MEs.csv")
# remove second S181Fg04 associated with band 2
summer.2018.MEs = summer.2018.MEs[c(1:20,22:38),]

summer.2018.MEs.2 = summer.2018.MEs[,c(1:20,36,51,57,59:61)]

summer.2018.MEs.3 = summer.2018.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### Growth Signal #####
summer.2018.MEs.4 <- summer.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

summer.2018.MEs.4 <- summer.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


###### GR #####
summer.2018.MEs.4 <- summer.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable
summer.2018.MEs.4 <- summer.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .05)

###### SLA #####
summer.2018.MEs.4 <- summer.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

summer.2018.MEs.4 <- summer.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


###### LA ######
summer.2018.MEs.4 <- summer.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

summer.2018.MEs.4 <- summer.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .05)
#
###### SPAD ######
summer.2018.MEs.4 <- summer.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

summer.2018.MEs.4 <- summer.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.MEs.4 <- summer.2018.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


###### turn growth into quantiles ######

traits = summer.2018.MEs[,c(36,51,59:61)]
rownames(traits) = summer.2018.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

##### Correlation #####

nSamples = nrow(summer.2018.MEs) # define number of samples
nGenes = ncol(summer.2018.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = summer.2018.MEs[,c(2:17)]
row.names(module_eigengenes) = summer.2018.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[17:23], y= names(heatmap.data)[1:16],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### fall.2018 #####
fall.2018.MEs = read.csv("./Data/grow.nogrow.MEs/fall.2018.MEs.csv")
fall.2018.MEs.2 = fall.2018.MEs[,c(1:21,37,52,58,60:62)]

fall.2018.MEs.3 = fall.2018.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### growth signal #####
fall.2018.MEs.4 <- fall.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable
fall.2018.MEs.4 <- fall.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR #####
fall.2018.MEs.4 <- fall.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

fall.2018.MEs.4 <- fall.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######
fall.2018.MEs.4 <- fall.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable
fall.2018.MEs.4 <- fall.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA ######
fall.2018.MEs.4 <- fall.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

fall.2018.MEs.4 <- fall.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD ######
fall.2018.MEs.4 <- fall.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# with site variable

fall.2018.MEs.4 <- fall.2018.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD + SITE, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.MEs.4 <- fall.2018.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles ######

traits = fall.2018.MEs[,c(37,52,60:62)]
rownames(traits) = fall.2018.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(fall.2018.MEs) # define number of samples
nGenes = ncol(fall.2018.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = fall.2018.MEs[,c(2:18)]
row.names(module_eigengenes) = fall.2018.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[18:24], y= names(heatmap.data)[1:17],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *


##### spring.2019 #####
spring.2019.MEs = read.csv("./Data/grow.nogrow.MEs/spring.2019.MEs.csv")
# remove second S191Fg18 associated with band 2
spring.2019.MEs = spring.2019.MEs[c(1:8,10:18),]

spring.2019.MEs.2 = spring.2019.MEs[,c(1:26,42,57,63,65:67)]

spring.2019.MEs.3 = spring.2019.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### Growth signal #####
spring.2019.MEs.4 <- spring.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2019.MEs.4 <- spring.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2019.MEs.4 <- spring.2019.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# 2 modeuls grey60, darkred
sig.modules = spring.2019.MEs.4[c(23,24,31,32),]
sig.modules$Plot[[2]]
sig.modules$Plot[[4]]

###### GR ######
spring.2019.MEs.4 <- spring.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2019.MEs.4 <- spring.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2019.MEs.4 <- spring.2019.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######
spring.2019.MEs.4 <- spring.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2019.MEs.4 <- spring.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2019.MEs.4 <- spring.2019.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA ######
spring.2019.MEs.4 <- spring.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2019.MEs.4 <- spring.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2019.MEs.4 <- spring.2019.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
###### SPAD ######
spring.2019.MEs.4 <- spring.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2019.MEs.4 <- spring.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2019.MEs.4 <- spring.2019.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles ######

traits = spring.2019.MEs[,c(42,57,65:67)]
rownames(traits) = spring.2019.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))

traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(spring.2019.MEs) # define number of samples
nGenes = ncol(spring.2019.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = spring.2019.MEs[,c(2:23)]
row.names(module_eigengenes) = spring.2019.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

module.trait.corr.pvals.3 = module.trait.corr.pvals.2 %>%
  filter(fdr.GS < 0.01)

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[23:29], y= names(heatmap.data)[1:22],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### summer.2019 #####
summer.2019.MEs = read.csv("./Data/grow.nogrow.MEs/summer.2019.MEs.csv")
# remove second S191Fg18 associated with band 2
summer.2019.MEs = summer.2019.MEs[c(1:9,11:18),]

summer.2019.MEs.2 = summer.2019.MEs[,c(1:43,59,74,80,82:84)]

summer.2019.MEs.3 = summer.2019.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### growth signal #####
summer.2019.MEs.4 <- summer.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2019.MEs.4 <- summer.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2019.MEs.4 <- summer.2019.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR #####
summer.2019.MEs.4 <- summer.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2019.MEs.4 <- summer.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2019.MEs.4 <- summer.2019.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA #####
summer.2019.MEs.4 <- summer.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2019.MEs.4 <- summer.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2019.MEs.4 <- summer.2019.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA #####
summer.2019.MEs.4 <- summer.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2019.MEs.4 <- summer.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2019.MEs.4 <- summer.2019.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD #####
summer.2019.MEs.4 <- summer.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2019.MEs.4 <- summer.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2019.MEs.4 <- summer.2019.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles #####

traits = summer.2019.MEs[,c(59,74,82:84)]
rownames(traits) = summer.2019.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)


###### Correlation ######

nSamples = nrow(summer.2019.MEs) # define number of samples
nGenes = ncol(summer.2019.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = summer.2019.MEs[,c(2:40)]
row.names(module_eigengenes) = summer.2019.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[40:46], y= names(heatmap.data)[1:39],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### fall.2019 #####
fall.2019.MEs = read.csv("./Data/grow.nogrow.MEs/fall.2019.MEs.csv")
fall.2019.MEs.2 = fall.2019.MEs[,c(1:30,46,61,67,69:71)]

fall.2019.MEs.3 = fall.2019.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### Growth signal #####
fall.2019.MEs.4 <- fall.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2019.MEs.4 <- fall.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2019.MEs.4 <- fall.2019.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# 1 module: royalblue
sig.modules = fall.2019.MEs.4[c(45,46),]
sig.modules$Plot[[2]]

###### GR #####
fall.2019.MEs.4 <- fall.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2019.MEs.4 <- fall.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2019.MEs.4 <- fall.2019.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA #####
fall.2019.MEs.4 <- fall.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2019.MEs.4 <- fall.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2019.MEs.4 <- fall.2019.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA #####
fall.2019.MEs.4 <- fall.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2019.MEs.4 <- fall.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2019.MEs.4 <- fall.2019.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
###### SPAD #####
fall.2019.MEs.4 <- fall.2019.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2019.MEs.4 <- fall.2019.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2019.MEs.4 <- fall.2019.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles #####

traits = fall.2019.MEs[,c(46,61,69:71)]
rownames(traits) = fall.2019.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(fall.2019.MEs) # define number of samples
nGenes = ncol(fall.2019.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = fall.2019.MEs[,c(2:27)]
row.names(module_eigengenes) = fall.2019.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

module.trait.corr.pvals.3 = module.trait.corr.pvals.2 %>%
  filter(fdr.GS < 0.01)

# royalblue

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[27:33], y= names(heatmap.data)[1:26],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### spring.2017.HF #####
spring.2017.HF.MEs = read.csv("./Data/grow.nogrow.MEs/spring.2017.HF.MEs.csv")
spring.2017.HF.MEs.2 = spring.2017.HF.MEs[,c(1:40,56,71,77,79:81)]
spring.2017.HF.MEs.2[is.na(spring.2017.HF.MEs.2)] <- 0

spring.2017.HF.MEs.3 = spring.2017.HF.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### Growth signal ######
spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR ######

spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# get 1 significant module: tan

sig.modules = spring.2017.HF.MEs.4[c(39,40),]
sig.modules$Plot[[2]]

###### SLA ######
spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA ######

spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD ######

spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.HF.MEs.4 <- spring.2017.HF.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### turn growth into quantiles ######

traits = spring.2017.HF.MEs[,c(56,71,79:81)]
rownames(traits) = spring.2017.HF.MEs$sample
traits[is.na(traits)] <- 0
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)


###### Correlation ######

nSamples = nrow(spring.2017.HF.MEs) # define number of samples
nGenes = ncol(spring.2017.HF.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = spring.2017.HF.MEs[,c(2:37)]
row.names(module_eigengenes) = spring.2017.HF.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

module.trait.corr.pvals.3 = module.trait.corr.pvals.2 %>%
  filter(fdr.GR < 0.01)

# tan

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[37:43], y= names(heatmap.data)[1:36],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### summer.2017.HF #####
summer.2017.HF.MEs = read.csv("./Data/grow.nogrow.MEs/summer.2017.HF.MEs.csv")
summer.2017.HF.MEs.2 = summer.2017.HF.MEs[,c(1:50,66,81,87,89:91)]

summer.2017.HF.MEs.3 = summer.2017.HF.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### Growth signal #####
summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR #####
summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA #####
summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA #####
summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD #####
summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.HF.MEs.4 <- summer.2017.HF.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles #####

traits = summer.2017.HF.MEs[,c(66,81,89:91)]
rownames(traits) = summer.2017.HF.MEs$sample
traits[is.na(traits)] <- 0
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation #####

nSamples = nrow(summer.2017.HF.MEs) # define number of samples
nGenes = ncol(summer.2017.HF.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = summer.2017.HF.MEs[,c(2:47)]
row.names(module_eigengenes) = summer.2017.HF.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[47:53], y= names(heatmap.data)[1:46],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### fall.2017.HF #####
fall.2017.HF.MEs = read.csv("./Data/grow.nogrow.MEs/fall.2017.HF.MEs.csv")
fall.2017.HF.MEs.2 = fall.2017.HF.MEs[,c(1:28,44,59,65,67:69)]

fall.2017.HF.MEs.3 = fall.2017.HF.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### growth signal #####
fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .05)

###### GR ######
fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .05)

###### SLA ######
fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .05)

# 5 significant modules: black, tan, green, purple, salmon

sig.mods = fall.2017.HF.MEs.4[c(9,10,11,12,15,16,17,18,19,20),]
sig.mods$Plot[[2]] # green
sig.mods$Plot[[4]] # salmon
sig.mods$Plot[[6]] # purple
sig.mods$Plot[[8]] # black
sig.mods$Plot[[10]] # tan

###### LA ######
fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .05)

###### SPAD ######
fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.HF.MEs.4 <- fall.2017.HF.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .05)

###### turn growth into quantiles ######

traits = fall.2017.HF.MEs[,c(44,59,67:69)]
rownames(traits) = fall.2017.HF.MEs$sample
traits[is.na(traits)] <- 0
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(fall.2017.HF.MEs) # define number of samples
nGenes = ncol(fall.2017.HF.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = fall.2017.HF.MEs[,c(2:25)]
row.names(module_eigengenes) = fall.2017.HF.MEs$sample

module.trait.corr <- cor(module_eigengenes,traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[25:31], y= names(heatmap.data)[1:24],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### spring.2017.SERC #####
spring.2017.SERC.MEs = read.csv("./Data/grow.nogrow.MEs/spring.2017.SERC.MEs.csv")
spring.2017.SERC.MEs.2 = spring.2017.SERC.MEs[,c(1:57,73,88,94,96:98)]

spring.2017.SERC.MEs.3 = spring.2017.SERC.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### growth signal ######
spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR ######
spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######
spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA ######
spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD ######
spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2017.SERC.MEs.4 <- spring.2017.SERC.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles ######

traits = spring.2017.SERC.MEs[,c(73,88,96:98)]
rownames(traits) = spring.2017.SERC.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(spring.2017.SERC.MEs) # define number of samples
nGenes = ncol(spring.2017.SERC.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = spring.2017.SERC.MEs[,c(2:54)]
row.names(module_eigengenes) = spring.2017.SERC.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

module.trait.corr.pvals.3 = module.trait.corr.pvals.2 %>%
  filter(fdr.high.grow < 0.01)

# lightsteelblue1

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[54:60], y= names(heatmap.data)[1:53],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *
# extract genes from modules with significance with high growth

##### summer.2017.SERC #####
summer.2017.SERC.MEs = read.csv("./Data/grow.nogrow.MEs/summer.2017.SERC.MEs.csv")
summer.2017.SERC.MEs.2 = summer.2017.SERC.MEs[,c(1:39,55,70,76,78:80)]

summer.2017.SERC.MEs.3 = summer.2017.SERC.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### growth signal #######
summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR #######
summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA #######
summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA #######
summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD #######
summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2017.SERC.MEs.4 <- summer.2017.SERC.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### turn growth into quantiles #####

traits = summer.2017.SERC.MEs[,c(55,70,78:80)]
rownames(traits) = summer.2017.SERC.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

##### Correlation #####

nSamples = nrow(summer.2017.SERC.MEs) # define number of samples
nGenes = ncol(summer.2017.SERC.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = summer.2017.SERC.MEs[,c(2:36)]
row.names(module_eigengenes) = summer.2017.SERC.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[36:42], y= names(heatmap.data)[1:35],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### fall.2017.SERC #####
fall.2017.SERC.MEs = read.csv("./Data/grow.nogrow.MEs/fall.2017.SERC.MEs.csv")
fall.2017.SERC.MEs.2 = fall.2017.SERC.MEs[,c(1:58,74,89,95,97:99)]

fall.2017.SERC.MEs.3 = fall.2017.SERC.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### growth signal ######
fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR ######
fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######
fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA ######
fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD ######
fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2017.SERC.MEs.4 <- fall.2017.SERC.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles ######

traits = fall.2017.SERC.MEs[,c(74,89,97:99)]
rownames(traits) = fall.2017.SERC.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(fall.2017.SERC.MEs) # define number of samples
nGenes = ncol(fall.2017.SERC.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = fall.2017.SERC.MEs[,c(2:55)]
row.names(module_eigengenes) = fall.2017.SERC.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[55:61], y= names(heatmap.data)[1:54],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### spring.2018.HF #####
spring.2018.HF.MEs = read.csv("./Data/grow.nogrow.MEs/spring.2018.HF.MEs.csv")
spring.2018.HF.MEs.2 = spring.2018.HF.MEs[,c(1:44,60,75,81,83:85)]

spring.2018.HF.MEs.3 = spring.2018.HF.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

# All growth signal values are 1 so can't do growth signal

###### GR ######
spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######
spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA ######
spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD ######
spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.HF.MEs.4 <- spring.2018.HF.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles ######

traits = as.data.frame(spring.2018.HF.MEs[,c(60,83:85)])
rownames(traits) = spring.2018.HF.MEs$sample
colnames(traits)[1] = "GR"
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4)],traits.2)

###### Correlation ######

nSamples = nrow(spring.2018.HF.MEs) # define number of samples
nGenes = ncol(spring.2018.HF.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = spring.2018.HF.MEs[,c(2:41)]
row.names(module_eigengenes) = spring.2018.HF.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[41:46], y= names(heatmap.data)[1:40],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### summer.2018.HF #####
summer.2018.HF.MEs = read.csv("./Data/grow.nogrow.MEs/summer.2018.HF.MEs.csv")
summer.2018.HF.MEs.2 = summer.2018.HF.MEs[,c(1:27,43,58,64,66:68)]

summer.2018.HF.MEs.3 = summer.2018.HF.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

# All values are 1 so can't do growth signal

###### GR ######
summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######
summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA ######
summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD ######
summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.HF.MEs.4 <- summer.2018.HF.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### turn growth into quantiles #####

traits = as.data.frame(summer.2018.HF.MEs[,c(43,66:68)])
rownames(traits) = summer.2018.HF.MEs$sample
colnames(traits)[1] = "GR"
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4)],traits.2)

###### Correlation ######

nSamples = nrow(summer.2018.HF.MEs) # define number of samples
nGenes = ncol(summer.2018.HF.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = summer.2018.HF.MEs[,c(2:24)]
row.names(module_eigengenes) = summer.2018.HF.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[24:29], y= names(heatmap.data)[1:23],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *
# extract genes from modules with significance with high growth

##### fall.2018.HF #####
fall.2018.HF.MEs = read.csv("./Data/grow.nogrow.MEs/fall.2018.HF.MEs.csv")
fall.2018.HF.MEs.2 = fall.2018.HF.MEs[,c(1:52,68,83,89,91:93)]

fall.2018.HF.MEs.3 = fall.2018.HF.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### growth signal ######
fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR ######
fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######
fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA ######
fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD ######
fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.HF.MEs.4 <- fall.2018.HF.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles ######

traits = fall.2018.HF.MEs[,c(68,83,91:93)]
rownames(traits) = fall.2018.HF.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(fall.2018.HF.MEs) # define number of samples
nGenes = ncol(fall.2018.HF.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = fall.2018.HF.MEs[,c(2:49)]
row.names(module_eigengenes) = fall.2018.HF.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[49:55], y= names(heatmap.data)[1:48],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### spring.2018.SERC #####
spring.2018.SERC.MEs = read.csv("./Data/grow.nogrow.MEs/spring.2018.SERC.MEs.csv")
# remove Band 2 sample
spring.2018.SERC.MEs = spring.2018.SERC.MEs[c(1:3,5:21),]

spring.2018.SERC.MEs.2 = spring.2018.SERC.MEs[,c(1:30,46,61,67,69:71)]

spring.2018.SERC.MEs.3 = spring.2018.SERC.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### growth signal ######
spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR ######
spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######
spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA ######
spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD ######
spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_spring.2018.SERC.MEs.4 <- spring.2018.SERC.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles ######

traits = spring.2018.SERC.MEs[,c(46,61,69:71)]
rownames(traits) = spring.2018.SERC.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(spring.2018.SERC.MEs) # define number of samples
nGenes = ncol(spring.2018.SERC.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = spring.2018.SERC.MEs[,c(2:27)]
row.names(module_eigengenes) = spring.2018.SERC.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[27:33], y= names(heatmap.data)[1:26],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### summer.2018.SERC #####
summer.2018.SERC.MEs = read.csv("./Data/grow.nogrow.MEs/summer.2018.SERC.MEs.csv")
# remove Band 2 sample
summer.2018.SERC.MEs = summer.2018.SERC.MEs[c(1:3,5:21),]

summer.2018.SERC.MEs.2 = summer.2018.SERC.MEs[,c(1:21,37,52,58,60:62)]

summer.2018.SERC.MEs.3 = summer.2018.SERC.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### Growth signal ######
summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### GR######
summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######
summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA ######
summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD ######
summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_summer.2018.SERC.MEs.4 <- summer.2018.SERC.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles ######

traits = summer.2018.SERC.MEs[,c(37,52,60:62)]
rownames(traits) = summer.2018.SERC.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(summer.2018.SERC.MEs) # define number of samples
nGenes = ncol(summer.2018.SERC.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = summer.2018.SERC.MEs[,c(2:18)]
row.names(module_eigengenes) = summer.2018.SERC.MEs$sample

module.trait.corr <- cor(module_eigengenes, traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[18:24], y= names(heatmap.data)[1:17],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *

##### fall.2018.SERC #####
fall.2018.SERC.MEs = read.csv("./Data/grow.nogrow.MEs/fall.2018.SERC.MEs.csv")
fall.2018.SERC.MEs.2 = fall.2018.SERC.MEs[,c(1:60,76,91,97,99:101)]

fall.2018.SERC.MEs.3 = fall.2018.SERC.MEs.2 %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,eigen_value,GROWTH_SIGNAL,GR,Average.SLA,Average.LA,Average.SPAD))

###### growth signal ######
fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GROWTH_SIGNAL, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GROWTH_SIGNAL, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.4 %>%
  filter(term == "GROWTH_SIGNAL") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
###### GR ######
fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ GR, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = GR, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.4 %>%
  filter(term == "GR") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SLA ######
fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SLA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SLA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.4 %>%
  filter(term == "Average.SLA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### LA ######
fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.LA, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.LA, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.4 %>%
  filter(term == "Average.LA") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### SPAD ######
fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.3 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ Average.SPAD, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.4 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = Average.SPAD, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_fall.2018.SERC.MEs.4 <- fall.2018.SERC.MEs.4 %>%
  filter(term == "Average.SPAD") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

###### turn growth into quantiles ######

traits = fall.2018.SERC.MEs[,c(76,91,99:101)]
rownames(traits) = fall.2018.SERC.MEs$sample
traits = traits %>% 
  mutate(quantilegroup = ntile(GR, 3)) 
traits$quantilegroup = factor(traits$quantilegroup, levels = c("1","2","3"))
traits.2 = binarizeCategoricalColumns(traits$quantilegroup,
                                      includePairwise = FALSE,
                                      includeLevelVsAll = TRUE,
                                      minCount = 1)

traits.3 = cbind(traits[,c(1,2,3,4,5)],traits.2)

###### Correlation ######

nSamples = nrow(fall.2018.SERC.MEs) # define number of samples
nGenes = ncol(fall.2018.SERC.MEs) # define number of genes

# correlation between module eigenegens and growth states
module_eigengenes = fall.2018.SERC.MEs[,c(2:57)]
row.names(module_eigengenes) = fall.2018.SERC.MEs$sample

module.trait.corr <- cor(module_eigengenes,traits.3, use = 'p') # pearson correlation
module.trait.corr.pvals <- as.data.frame(corPvalueStudent(module.trait.corr, nSamples)) # p-values for correlations

module.trait.corr.pvals.2 = module.trait.corr.pvals %>%
  mutate(fdr.GR = p.adjust(GR, method = "fdr"),
         fdr.GS = p.adjust(GROWTH_SIGNAL, method = "fdr"),
         fdr.mid.grow = p.adjust(data.2.vs.all, method = "fdr"),
         fdr.high.grow = p.adjust(data.3.vs.all, method = "fdr"),
         fdr.SLA = p.adjust(Average.SLA, method = "fdr"),
         fdr.LA = p.adjust(Average.LA, method = "fdr"),
         fdr.SPAD = p.adjust(Average.SPAD, method = "fdr"))

# heat map of correlation
heatmap.data = cbind(module_eigengenes, traits.3) # combining data into one dataframe

# specify columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[57:63], y= names(heatmap.data)[1:56],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *
