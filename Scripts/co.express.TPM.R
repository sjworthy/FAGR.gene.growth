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
library(edgeR)

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

#### Expression Filtering ####
# Remove low expression genes
# Need at least 10 reads in at least 3 samples

TPM_ExprData <- TPM_ExprData[rowSums(TPM_ExprData[,-1] > 10) >= 3,]

#### split samples into time points ####
# also remove samples where we don't have growth data

# split TPMS for each year by each season
TPM_2017_Spring = TPM_ExprData[,c(1:16,18:21,108:127)] # 39 samples
TPM_2017_Summer = TPM_ExprData[,c(1,22:34,36:39,128:146)] # 36 samples
TPM_2017_Fall = TPM_ExprData[,c(1,40:54,56:59,147:166)] # 39 samples
TPM_2018_Spring = TPM_ExprData[,c(1,60:72,167:168,170:186)] # 32 samples
TPM_2018_Summer = TPM_ExprData[,c(1,73:89,187:188,190:206)] # 36 samples
TPM_2018_Fall = TPM_ExprData[,c(1,90:101,103:107,207:223)] # 34 samples
TPM_2019_Spring = TPM_ExprData[,c(1,224:230,232:240)] # 16 samples
TPM_2019_Summer = TPM_ExprData[,c(1,241:248,250:259)] # 18 samples
TPM_2019_Fall = TPM_ExprData[,c(1,260:276)] # 17 samples
TPM_2020_Spring = TPM_ExprData[,c(1,277,279:284,286:291,294:295)] # 15 samples

# split TPMS for each year by each season by each site
# only for 2017 and 2018 where we have both sites

TPM_2017_Spring_HF = TPM_2017_Spring[,c(1:20)]
TPM_2017_Spring_SERC = TPM_2017_Spring[,c(1,21:40)]
TPM_2018_Spring_HF = TPM_2018_Spring[,c(1:14)]
TPM_2018_Spring_SERC = TPM_2018_Spring[,c(1,15:33)]
TPM_2017_Summer_HF = TPM_2017_Summer[,c(1:18)]
TPM_2017_Summer_SERC = TPM_2017_Summer[,c(1,19:37)]
TPM_2018_Summer_HF = TPM_2018_Summer[,c(1:18)]
TPM_2018_Summer_SERC = TPM_2018_Summer[,c(1,19:37)]
TPM_2017_Fall_HF = TPM_2017_Fall[,c(1:20)]
TPM_2017_Fall_SERC = TPM_2017_Fall[,c(1,21:40)]
TPM_2018_Fall_HF = TPM_2018_Fall[,c(1:18)]
TPM_2018_Fall_SERC = TPM_2018_Fall[,c(1,19:35)]

#### log transform read counts ####
# convert counts to matrix 
TPM_2017_Spring_mat <- TPM_2017_Spring %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2017_Spring_mat) <- TPM_2017_Spring$Gene_ID

TPM_2017_Summer_mat <- TPM_2017_Summer %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2017_Summer_mat) <- TPM_2017_Summer$Gene_ID

TPM_2017_Fall_mat <- TPM_2017_Fall %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2017_Fall_mat) <- TPM_2017_Fall$Gene_ID

TPM_2018_Spring_mat <- TPM_2018_Spring %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2018_Spring_mat) <- TPM_2018_Spring$Gene_ID

TPM_2018_Summer_mat <- TPM_2018_Summer %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2018_Summer_mat) <- TPM_2018_Summer$Gene_ID

TPM_2018_Fall_mat <- TPM_2018_Fall %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2018_Fall_mat) <- TPM_2018_Fall$Gene_ID

TPM_2019_Spring_mat <- TPM_2019_Spring %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2019_Spring_mat) <- TPM_2019_Spring$Gene_ID

TPM_2019_Summer_mat <- TPM_2019_Summer %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2019_Summer_mat) <- TPM_2019_Summer$Gene_ID

TPM_2019_Fall_mat <- TPM_2019_Fall %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2019_Fall_mat) <- TPM_2019_Fall$Gene_ID

TPM_2020_Spring_mat <- TPM_2020_Spring %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2020_Spring_mat) <- TPM_2020_Spring$Gene_ID

TPM_2017_Spring_HF_mat <- TPM_2017_Spring_HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2017_Spring_HF_mat) <- TPM_2017_Spring_HF$Gene_ID

TPM_2017_Summer_HF_mat <- TPM_2017_Summer_HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2017_Summer_HF_mat) <- TPM_2017_Summer_HF$Gene_ID

TPM_2017_Fall_HF_mat <- TPM_2017_Fall_HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2017_Fall_HF_mat) <- TPM_2017_Fall_HF$Gene_ID

TPM_2018_Spring_HF_mat <- TPM_2018_Spring_HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2018_Spring_HF_mat) <- TPM_2018_Spring_HF$Gene_ID

TPM_2018_Summer_HF_mat <- TPM_2018_Summer_HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2018_Summer_HF_mat) <- TPM_2018_Summer_HF$Gene_ID

TPM_2018_Fall_HF_mat <- TPM_2018_Fall_HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2018_Fall_HF_mat) <- TPM_2018_Fall_HF$Gene_ID

TPM_2017_Spring_SERC_mat <- TPM_2017_Spring_SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2017_Spring_SERC_mat) <- TPM_2017_Spring_SERC$Gene_ID

TPM_2017_Summer_SERC_mat <- TPM_2017_Summer_SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2017_Summer_SERC_mat) <- TPM_2017_Summer_SERC$Gene_ID

TPM_2017_Fall_SERC_mat <- TPM_2017_Fall_SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2017_Fall_SERC_mat) <- TPM_2017_Fall_SERC$Gene_ID

TPM_2018_Spring_SERC_mat <- TPM_2018_Spring_SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2018_Spring_SERC_mat) <- TPM_2018_Spring_SERC$Gene_ID

TPM_2018_Summer_SERC_mat <- TPM_2018_Summer_SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2018_Summer_SERC_mat) <- TPM_2018_Summer_SERC$Gene_ID

TPM_2018_Fall_SERC_mat <- TPM_2018_Fall_SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(TPM_2018_Fall_SERC_mat) <- TPM_2018_Fall_SERC$Gene_ID

# log2 transform using cpm, convert to df, and add Gene_ID back as a column
TPM_2017_Spring_log = as.data.frame(cpm(TPM_2017_Spring_mat, log = TRUE))
TPM_2017_Spring_log$Gene_ID = rownames(TPM_2017_Spring_log)
TPM_2017_Summer_log = as.data.frame(cpm(TPM_2017_Summer_mat, log = TRUE))
TPM_2017_Summer_log$Gene_ID = rownames(TPM_2017_Summer_log)
TPM_2017_Fall_log = as.data.frame(cpm(TPM_2017_Fall_mat, log = TRUE))
TPM_2017_Fall_log$Gene_ID = rownames(TPM_2017_Fall_log)
TPM_2018_Spring_log = as.data.frame(cpm(TPM_2018_Spring_mat, log = TRUE))
TPM_2018_Spring_log$Gene_ID = rownames(TPM_2018_Spring_log)
TPM_2018_Summer_log = as.data.frame(cpm(TPM_2018_Summer_mat, log = TRUE))
TPM_2018_Summer_log$Gene_ID = rownames(TPM_2018_Summer_log)
TPM_2018_Fall_log = as.data.frame(cpm(TPM_2018_Fall_mat, log = TRUE))
TPM_2018_Fall_log$Gene_ID = rownames(TPM_2018_Fall_log)
TPM_2019_Spring_log = as.data.frame(cpm(TPM_2019_Spring_mat, log = TRUE))
TPM_2019_Spring_log$Gene_ID = rownames(TPM_2019_Spring_log)
TPM_2019_Summer_log = as.data.frame(cpm(TPM_2019_Summer_mat, log = TRUE))
TPM_2019_Summer_log$Gene_ID = rownames(TPM_2019_Summer_log)
TPM_2019_Fall_log = as.data.frame(cpm(TPM_2019_Fall_mat, log = TRUE))
TPM_2019_Fall_log$Gene_ID = rownames(TPM_2019_Fall_log)
TPM_2020_Spring_log = as.data.frame(cpm(TPM_2020_Spring_mat, log = TRUE))
TPM_2020_Spring_log$Gene_ID = rownames(TPM_2020_Spring_log)

TPM_2017_Spring_HF_log = as.data.frame(cpm(TPM_2017_Spring_HF_mat, log = TRUE))
TPM_2017_Spring_HF_log$Gene_ID = rownames(TPM_2017_Spring_HF_log)
TPM_2017_Summer_HF_log = as.data.frame(cpm(TPM_2017_Summer_HF_mat, log = TRUE))
TPM_2017_Summer_HF_log$Gene_ID = rownames(TPM_2017_Summer_HF_log)
TPM_2017_Fall_HF_log = as.data.frame(cpm(TPM_2017_Fall_HF_mat, log = TRUE))
TPM_2017_Fall_HF_log$Gene_ID = rownames(TPM_2017_Fall_HF_log)
TPM_2018_Spring_HF_log = as.data.frame(cpm(TPM_2018_Spring_HF_mat, log = TRUE))
TPM_2018_Spring_HF_log$Gene_ID = rownames(TPM_2018_Spring_HF_log)
TPM_2018_Summer_HF_log = as.data.frame(cpm(TPM_2018_Summer_HF_mat, log = TRUE))
TPM_2018_Summer_HF_log$Gene_ID = rownames(TPM_2018_Summer_HF_log)
TPM_2018_Fall_HF_log = as.data.frame(cpm(TPM_2018_Fall_HF_mat, log = TRUE))
TPM_2018_Fall_HF_log$Gene_ID = rownames(TPM_2018_Fall_HF_log)
TPM_2017_Spring_SERC_log = as.data.frame(cpm(TPM_2017_Spring_SERC_mat, log = TRUE))
TPM_2017_Spring_SERC_log$Gene_ID = rownames(TPM_2017_Spring_SERC_log)
TPM_2017_Summer_SERC_log = as.data.frame(cpm(TPM_2017_Summer_SERC_mat, log = TRUE))
TPM_2017_Summer_SERC_log$Gene_ID = rownames(TPM_2017_Summer_SERC_log)
TPM_2017_Fall_SERC_log = as.data.frame(cpm(TPM_2017_Fall_SERC_mat, log = TRUE))
TPM_2017_Fall_SERC_log$Gene_ID = rownames(TPM_2017_Fall_SERC_log)
TPM_2018_Spring_SERC_log = as.data.frame(cpm(TPM_2018_Spring_SERC_mat, log = TRUE))
TPM_2018_Spring_SERC_log$Gene_ID = rownames(TPM_2018_Spring_SERC_log)
TPM_2018_Summer_SERC_log = as.data.frame(cpm(TPM_2018_Summer_SERC_mat, log = TRUE))
TPM_2018_Summer_SERC_log$Gene_ID = rownames(TPM_2018_Summer_SERC_log)
TPM_2018_Fall_SERC_log = as.data.frame(cpm(TPM_2018_Fall_SERC_mat, log = TRUE))
TPM_2018_Fall_SERC_log$Gene_ID = rownames(TPM_2018_Fall_SERC_log)

#### Calculate the Variance of each Gene across all of our samples ####
TPM_2017_Fall_cv <- TPM_2017_Fall_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Spring_cv <- TPM_2017_Spring_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Summer_cv <- TPM_2017_Summer_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Fall_cv <- TPM_2018_Fall_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Spring_cv <- TPM_2018_Spring_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Summer_cv <- TPM_2018_Summer_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2019_Fall_cv <- TPM_2019_Fall_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2019_Spring_cv <- TPM_2019_Spring_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2019_Summer_cv <- TPM_2019_Summer_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2020_Spring_cv <- TPM_2020_Spring_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

TPM_2017_Fall_HF_cv <- TPM_2017_Fall_HF_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Spring_HF_cv <- TPM_2017_Spring_HF_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Summer_HF_cv <- TPM_2017_Summer_HF_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Fall_SERC_cv <- TPM_2017_Fall_SERC_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Spring_SERC_cv <- TPM_2017_Spring_SERC_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Summer_SERC_cv <- TPM_2017_Summer_SERC_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Fall_HF_cv <- TPM_2018_Fall_HF_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Spring_HF_cv <- TPM_2018_Spring_HF_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Summer_HF_cv <- TPM_2018_Summer_HF_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Fall_SERC_cv <- TPM_2018_Fall_SERC_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Spring_SERC_cv <- TPM_2018_Spring_SERC_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Summer_SERC_cv <- TPM_2018_Summer_SERC_log %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

# Filter the data frame to contain the top 30% most variable genes)
# ties are kept together
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

TPM_2017_Fall_HF_30  <- TPM_2017_Fall_HF_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Spring_HF_30  <- TPM_2017_Spring_HF_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Summer_HF_30  <- TPM_2017_Summer_HF_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Fall_HF_30  <- TPM_2018_Fall_HF_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Spring_HF_30  <- TPM_2018_Spring_HF_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Summer_HF_30  <- TPM_2018_Summer_HF_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Fall_SERC_30  <- TPM_2017_Fall_SERC_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Spring_SERC_30  <- TPM_2017_Spring_SERC_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Summer_SERC_30  <- TPM_2017_Summer_SERC_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Fall_SERC_30  <- TPM_2018_Fall_SERC_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Spring_SERC_30  <- TPM_2018_Spring_SERC_cv  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Summer_SERC_30  <- TPM_2018_Summer_SERC_cv  %>% slice_max(order_by = cv, prop = .30)

# Deselect the cv column and flip our data frame to contain sample along the rows and genes along the columns
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

TPM_2017_Fall_HF_30 <- select(TPM_2017_Fall_HF_30, -cv)
TPM_2017_Fall_HF_30 <- column_to_rownames(TPM_2017_Fall_HF_30,var = "Gene_ID")
TPM_2017_Fall_HF_30 <- as.data.frame(t(TPM_2017_Fall_HF_30))

TPM_2017_Spring_HF_30 <- select(TPM_2017_Spring_HF_30, -cv)
TPM_2017_Spring_HF_30 <- column_to_rownames(TPM_2017_Spring_HF_30,var = "Gene_ID")
TPM_2017_Spring_HF_30 <- as.data.frame(t(TPM_2017_Spring_HF_30))

TPM_2017_Summer_HF_30 <- select(TPM_2017_Summer_HF_30, -cv)
TPM_2017_Summer_HF_30 <- column_to_rownames(TPM_2017_Summer_HF_30,var = "Gene_ID")
TPM_2017_Summer_HF_30 <- as.data.frame(t(TPM_2017_Summer_HF_30))

TPM_2018_Fall_HF_30 <- select(TPM_2018_Fall_HF_30, -cv)
TPM_2018_Fall_HF_30 <- column_to_rownames(TPM_2018_Fall_HF_30,var = "Gene_ID")
TPM_2018_Fall_HF_30 <- as.data.frame(t(TPM_2018_Fall_HF_30))

TPM_2018_Spring_HF_30 <- select(TPM_2018_Spring_HF_30, -cv)
TPM_2018_Spring_HF_30 <- column_to_rownames(TPM_2018_Spring_HF_30,var = "Gene_ID")
TPM_2018_Spring_HF_30 <- as.data.frame(t(TPM_2018_Spring_HF_30))

TPM_2018_Summer_HF_30 <- select(TPM_2018_Summer_HF_30, -cv)
TPM_2018_Summer_HF_30 <- column_to_rownames(TPM_2018_Summer_HF_30,var = "Gene_ID")
TPM_2018_Summer_HF_30 <- as.data.frame(t(TPM_2018_Summer_HF_30))

TPM_2017_Fall_SERC_30 <- select(TPM_2017_Fall_SERC_30, -cv)
TPM_2017_Fall_SERC_30 <- column_to_rownames(TPM_2017_Fall_SERC_30,var = "Gene_ID")
TPM_2017_Fall_SERC_30 <- as.data.frame(t(TPM_2017_Fall_SERC_30))

TPM_2017_Spring_SERC_30 <- select(TPM_2017_Spring_SERC_30, -cv)
TPM_2017_Spring_SERC_30 <- column_to_rownames(TPM_2017_Spring_SERC_30,var = "Gene_ID")
TPM_2017_Spring_SERC_30 <- as.data.frame(t(TPM_2017_Spring_SERC_30))

TPM_2017_Summer_SERC_30 <- select(TPM_2017_Summer_SERC_30, -cv)
TPM_2017_Summer_SERC_30 <- column_to_rownames(TPM_2017_Summer_SERC_30,var = "Gene_ID")
TPM_2017_Summer_SERC_30 <- as.data.frame(t(TPM_2017_Summer_SERC_30))

TPM_2018_Fall_SERC_30 <- select(TPM_2018_Fall_SERC_30, -cv)
TPM_2018_Fall_SERC_30 <- column_to_rownames(TPM_2018_Fall_SERC_30,var = "Gene_ID")
TPM_2018_Fall_SERC_30 <- as.data.frame(t(TPM_2018_Fall_SERC_30))

TPM_2018_Spring_SERC_30 <- select(TPM_2018_Spring_SERC_30, -cv)
TPM_2018_Spring_SERC_30 <- column_to_rownames(TPM_2018_Spring_SERC_30,var = "Gene_ID")
TPM_2018_Spring_SERC_30 <- as.data.frame(t(TPM_2018_Spring_SERC_30))

TPM_2018_Summer_SERC_30 <- select(TPM_2018_Summer_SERC_30, -cv)
TPM_2018_Summer_SERC_30 <- column_to_rownames(TPM_2018_Summer_SERC_30,var = "Gene_ID")
TPM_2018_Summer_SERC_30 <- as.data.frame(t(TPM_2018_Summer_SERC_30))

# Continue with formatting
Sample_Description = read_csv("./Formatted.Data/FAGR.description.csv")

#### QC ####
# iterative filtering of samples and genes with too many missing entries
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

gsg_2017_Fall_HF = goodSamplesGenes(TPM_2017_Fall_HF_30, verbose = 3);
gsg_2017_Fall_HF$allOK
gsg_2017_Spring_HF = goodSamplesGenes(TPM_2017_Spring_HF_30, verbose = 3);
gsg_2017_Spring_HF$allOK
gsg_2017_Summer_HF = goodSamplesGenes(TPM_2017_Summer_HF_30, verbose = 3);
gsg_2017_Summer_HF$allOK
gsg_2018_Fall_HF = goodSamplesGenes(TPM_2018_Fall_HF_30, verbose = 3);
gsg_2018_Fall_HF$allOK
gsg_2018_Spring_HF = goodSamplesGenes(TPM_2018_Spring_HF_30, verbose = 3);
gsg_2018_Spring_HF$allOK
gsg_2018_Summer_HF = goodSamplesGenes(TPM_2018_Summer_HF_30, verbose = 3);
gsg_2018_Summer_HF$allOK
gsg_2017_Fall_SERC = goodSamplesGenes(TPM_2017_Fall_SERC_30, verbose = 3);
gsg_2017_Fall_SERC$allOK
gsg_2017_Spring_SERC = goodSamplesGenes(TPM_2017_Spring_SERC_30, verbose = 3);
gsg_2017_Spring_SERC$allOK
gsg_2017_Summer_SERC = goodSamplesGenes(TPM_2017_Summer_SERC_30, verbose = 3);
gsg_2017_Summer_SERC$allOK
gsg_2018_Fall_SERC = goodSamplesGenes(TPM_2018_Fall_SERC_30, verbose = 3);
gsg_2018_Fall_SERC$allOK
gsg_2018_Spring_SERC = goodSamplesGenes(TPM_2018_Spring_SERC_30, verbose = 3);
gsg_2018_Spring_SERC$allOK
gsg_2018_Summer_SERC = goodSamplesGenes(TPM_2018_Summer_SERC_30, verbose = 3);
gsg_2018_Summer_SERC$allOK

#### Soft Threshold ####
# Choose a set of soft-thresholding powers
# value is the power to which co-expression similarity is raised to calculate adjacency
# chose the value closest to 0.90
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft_2017_Fall <- pickSoftThreshold(TPM_2017_Fall_30, powerVector = powers, verbose = 5)
sft_2017_Spring <- pickSoftThreshold(TPM_2017_Spring_30, powerVector = powers, verbose = 5)
sft_2017_Summer <- pickSoftThreshold(TPM_2017_Summer_30, powerVector = powers, verbose = 5)
sft_2018_Fall <- pickSoftThreshold(TPM_2018_Fall_30, powerVector = powers, verbose = 5)
sft_2018_Spring <- pickSoftThreshold(TPM_2018_Spring_30, powerVector = powers, verbose = 5)
sft_2018_Summer <- pickSoftThreshold(TPM_2018_Summer_30, powerVector = powers, verbose = 5)
sft_2019_Fall <- pickSoftThreshold(TPM_2019_Fall_30, powerVector = powers, verbose = 5)
sft_2019_Spring <- pickSoftThreshold(TPM_2019_Spring_30, powerVector = powers, verbose = 5)
sft_2019_Summer <- pickSoftThreshold(TPM_2019_Summer_30, powerVector = powers, verbose = 5)
sft_2020_Spring <- pickSoftThreshold(TPM_2020_Spring_30, powerVector = powers, verbose = 5)

sft_2017_Fall_HF <- pickSoftThreshold(TPM_2017_Fall_HF_30, powerVector = powers, verbose = 5)
sft_2017_Spring_HF <- pickSoftThreshold(TPM_2017_Spring_HF_30, powerVector = powers, verbose = 5)
sft_2017_Summer_HF <- pickSoftThreshold(TPM_2017_Summer_HF_30, powerVector = powers, verbose = 5)
sft_2018_Fall_HF <- pickSoftThreshold(TPM_2018_Fall_HF_30, powerVector = powers, verbose = 5)
sft_2018_Spring_HF <- pickSoftThreshold(TPM_2018_Spring_HF_30, powerVector = powers, verbose = 5)
sft_2018_Summer_HF <- pickSoftThreshold(TPM_2018_Summer_HF_30, powerVector = powers, verbose = 5)
sft_2017_Fall_SERC <- pickSoftThreshold(TPM_2017_Fall_SERC_30, powerVector = powers, verbose = 5)
sft_2017_Spring_SERC <- pickSoftThreshold(TPM_2017_Spring_SERC_30, powerVector = powers, verbose = 5)
sft_2017_Summer_SERC <- pickSoftThreshold(TPM_2017_Summer_SERC_30, powerVector = powers, verbose = 5)
sft_2018_Fall_SERC <- pickSoftThreshold(TPM_2018_Fall_SERC_30, powerVector = powers, verbose = 5)
sft_2018_Spring_SERC <- pickSoftThreshold(TPM_2018_Spring_SERC_30, powerVector = powers, verbose = 5)
sft_2018_Summer_SERC <- pickSoftThreshold(TPM_2018_Summer_SERC_30, powerVector = powers, verbose = 5)

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

adjacency_2017_Fall = adjacency(TPM_2017_Fall_30, power = 5,type = "signed hybrid")
adjacency_2017_Spring = adjacency(TPM_2017_Spring_30, power = 7,type = "signed hybrid")
adjacency_2017_Summer = adjacency(TPM_2017_Summer_30, power = 8,type = "signed hybrid")
adjacency_2018_Fall = adjacency(TPM_2018_Fall_30, power = 8,type = "signed hybrid")
adjacency_2018_Spring = adjacency(TPM_2018_Spring_30, power = 4,type = "signed hybrid")
adjacency_2018_Summer = adjacency(TPM_2018_Summer_30, power = 6,type = "signed hybrid")
adjacency_2019_Fall = adjacency(TPM_2019_Fall_30, power = 8,type = "signed hybrid")
adjacency_2019_Spring = adjacency(TPM_2019_Spring_30, power = 12,type = "signed hybrid")
adjacency_2019_Summer = adjacency(TPM_2019_Summer_30, power = 10,type = "signed hybrid")
adjacency_2020_Spring = adjacency(TPM_2020_Spring_30, power = 5,type = "signed hybrid")
# make TOMS below with the above adjacney matrices first and remove to save space.

adjacency_2017_Fall_HF = adjacency(TPM_2017_Fall_HF_30, power = 12,type = "signed hybrid")
adjacency_2017_Spring_HF = adjacency(TPM_2017_Spring_HF_30, power = 6,type = "signed hybrid")
adjacency_2017_Summer_HF = adjacency(TPM_2017_Summer_HF_30, power = 7,type = "signed hybrid")
adjacency_2018_Fall_HF = adjacency(TPM_2018_Fall_HF_30, power = 5,type = "signed hybrid")
adjacency_2018_Spring_HF = adjacency(TPM_2018_Spring_HF_30, power = 12,type = "signed hybrid")
adjacency_2018_Summer_HF = adjacency(TPM_2018_Summer_HF_30, power = 7,type = "signed hybrid")
adjacency_2017_Fall_SERC = adjacency(TPM_2017_Fall_SERC_30, power = 6,type = "signed hybrid")
adjacency_2017_Spring_SERC = adjacency(TPM_2017_Spring_SERC_30, power = 10,type = "signed hybrid")
adjacency_2017_Summer_SERC = adjacency(TPM_2017_Summer_SERC_30, power = 12,type = "signed hybrid")
adjacency_2018_Fall_SERC = adjacency(TPM_2018_Fall_SERC_30, power = 6,type = "signed hybrid")
adjacency_2018_Spring_SERC = adjacency(TPM_2018_Spring_SERC_30, power = 7,type = "signed hybrid")
adjacency_2018_Summer_SERC = adjacency(TPM_2018_Summer_SERC_30, power = 5,type = "signed hybrid")

#### Topological overlap matrix (TOM) ####
# to minimize effects of noise and spurious associations, tranform adjacency into TOM, calculate dissimilarity
# need to rm adjacency matrices and TOM after use to make space

TOM_2017_Fall = TOMsimilarity(adjacency_2017_Fall, TOMType = "signed Nowick")
dissTOM_2017_Fall = 1 - TOM_2017_Fall
#saveRDS(TOM_2017_Fall, file = "./Data/TOM_2017_Fall.R")
rm(TOM_2017_Fall)
rm(adjacency_2017_Fall)
TOM_2017_Spring = TOMsimilarity(adjacency_2017_Spring, TOMType = "signed Nowick")
dissTOM_2017_Spring = 1 - TOM_2017_Spring
#saveRDS(TOM_2017_Spring, file = "./Data/TOM_2017_Spring.R")
rm(TOM_2017_Spring)
rm(adjacency_2017_Spring)
TOM_2017_Summer = TOMsimilarity(adjacency_2017_Summer, TOMType = "signed Nowick")
dissTOM_2017_Summer = 1 - TOM_2017_Summer
#saveRDS(TOM_2017_Summer, file = "./Data/TOM_2017_Summer.R")
rm(TOM_2017_Summer)
rm(adjacency_2017_Summer)
TOM_2018_Fall = TOMsimilarity(adjacency_2018_Fall, TOMType = "signed Nowick")
dissTOM_2018_Fall = 1 - TOM_2018_Fall
#saveRDS(TOM_2018_Fall, file = "./Data/TOM_2018_Fall.R")
rm(TOM_2018_Fall)
rm(adjacency_2018_Fall)
TOM_2018_Spring = TOMsimilarity(adjacency_2018_Spring, TOMType = "signed Nowick")
dissTOM_2018_Spring = 1 - TOM_2018_Spring
#saveRDS(TOM_2018_Spring, file = "./Data/TOM_2018_Spring.R")
rm(TOM_2018_Spring)
rm(adjacency_2018_Spring)
TOM_2018_Summer = TOMsimilarity(adjacency_2018_Summer, TOMType = "signed Nowick")
dissTOM_2018_Summer = 1 - TOM_2018_Summer
#saveRDS(TOM_2018_Summer, file = "./Data/TOM_2018_Summer.R")
rm(TOM_2018_Summer)
rm(adjacency_2018_Summer)
TOM_2019_Fall = TOMsimilarity(adjacency_2019_Fall, TOMType = "signed Nowick")
dissTOM_2019_Fall = 1 - TOM_2019_Fall
#saveRDS(TOM_2019_Fall, file = "./Data/TOM_2019_Fall.R")
rm(TOM_2019_Fall)
rm(adjacency_2019_Fall)
TOM_2019_Spring = TOMsimilarity(adjacency_2019_Spring, TOMType = "signed Nowick")
dissTOM_2019_Spring = 1 - TOM_2019_Spring
#saveRDS(TOM_2019_Spring, file = "./Data/TOM_2019_Spring.R")
rm(TOM_2019_Spring)
rm(adjacency_2019_Spring)
TOM_2019_Summer = TOMsimilarity(adjacency_2019_Summer, TOMType = "signed Nowick")
dissTOM_2019_Summer = 1 - TOM_2019_Summer
#saveRDS(TOM_2019_Summer, file = "./Data/TOM_2019_Summer.R")
rm(TOM_2019_Summer)
rm(adjacency_2019_Summer)
TOM_2020_Spring = TOMsimilarity(adjacency_2020_Spring, TOMType = "signed Nowick")
dissTOM_2020_Spring = 1 - TOM_2020_Spring
#saveRDS(TOM_2020_Spring, file = "./Data/TOM_2020_Spring.R")
rm(TOM_2020_Spring)
rm(adjacency_2020_Spring)

TOM_2017_Fall_HF = TOMsimilarity(adjacency_2017_Fall_HF, TOMType = "signed Nowick")
dissTOM_2017_Fall_HF = 1 - TOM_2017_Fall_HF
#saveRDS(TOM_2017_Fall_HF, file = "./Data/TOM_2017_Fall_HF.R")
rm(TOM_2017_Fall_HF)
rm(adjacency_2017_Fall_HF)
TOM_2017_Spring_HF = TOMsimilarity(adjacency_2017_Spring_HF, TOMType = "signed Nowick")
dissTOM_2017_Spring_HF = 1 - TOM_2017_Spring_HF
#saveRDS(TOM_2017_Spring_HF, file = "./Data/TOM_2017_Spring_HF.R")
rm(TOM_2017_Spring_HF)
rm(adjacency_2017_Spring_HF)
TOM_2017_Summer_HF = TOMsimilarity(adjacency_2017_Summer_HF, TOMType = "signed Nowick")
dissTOM_2017_Summer_HF = 1 - TOM_2017_Summer_HF
#saveRDS(TOM_2017_Summer_HF, file = "./Data/TOM_2017_Summer_HF.R")
rm(TOM_2017_Summer_HF)
rm(adjacency_2017_Summer_HF)
TOM_2018_Fall_HF = TOMsimilarity(adjacency_2018_Fall_HF, TOMType = "signed Nowick")
dissTOM_2018_Fall_HF = 1 - TOM_2018_Fall_HF
#saveRDS(TOM_2018_Fall_HF, file = "./Data/TOM_2018_Fall_HF.R")
rm(TOM_2018_Fall_HF)
rm(adjacency_2018_Fall_HF)
TOM_2018_Spring_HF = TOMsimilarity(adjacency_2018_Spring_HF, TOMType = "signed Nowick")
dissTOM_2018_Spring_HF = 1 - TOM_2018_Spring_HF
#saveRDS(TOM_2018_Spring_HF, file = "./Data/TOM_2018_Spring_HF.R")
rm(TOM_2018_Spring_HF)
rm(adjacency_2018_Spring_HF)
TOM_2018_Summer_HF = TOMsimilarity(adjacency_2018_Summer_HF, TOMType = "signed Nowick")
dissTOM_2018_Summer_HF = 1 - TOM_2018_Summer_HF
#saveRDS(TOM_2018_Summer_HF, file = "./Data/TOM_2018_Summer_HF.R")
rm(TOM_2018_Summer_HF)
rm(adjacency_2018_Summer_HF)
TOM_2017_Fall_SERC = TOMsimilarity(adjacency_2017_Fall_SERC, TOMType = "signed Nowick")
dissTOM_2017_Fall_SERC = 1 - TOM_2017_Fall_SERC
#saveRDS(TOM_2017_Fall_SERC, file = "./Data/TOM_2017_Fall_SERC.R")
rm(TOM_2017_Fall_SERC)
rm(adjacency_2017_Fall_SERC)
TOM_2017_Spring_SERC = TOMsimilarity(adjacency_2017_Spring_SERC, TOMType = "signed Nowick")
dissTOM_2017_Spring_SERC = 1 - TOM_2017_Spring_SERC
#saveRDS(TOM_2017_Spring_SERC, file = "./Data/TOM_2017_Spring_SERC.R")
rm(TOM_2017_Spring_SERC)
rm(adjacency_2017_Spring_SERC)
TOM_2017_Summer_SERC = TOMsimilarity(adjacency_2017_Summer_SERC, TOMType = "signed Nowick")
dissTOM_2017_Summer_SERC = 1 - TOM_2017_Summer_SERC
#saveRDS(TOM_2017_Summer_SERC, file = "./Data/TOM_2017_Summer_SERC.R")
rm(TOM_2017_Summer_SERC)
rm(adjacency_2017_Summer_SERC)
TOM_2018_Fall_SERC = TOMsimilarity(adjacency_2018_Fall_SERC, TOMType = "signed Nowick")
dissTOM_2018_Fall_SERC = 1 - TOM_2018_Fall_SERC
#saveRDS(TOM_2018_Fall_SERC, file = "./Data/TOM_2018_Fall_SERC.R")
rm(TOM_2018_Fall_SERC)
rm(adjacency_2018_Fall_SERC)
TOM_2018_Spring_SERC = TOMsimilarity(adjacency_2018_Spring_SERC, TOMType = "signed Nowick")
dissTOM_2018_Spring_SERC = 1 - TOM_2018_Spring_SERC
#saveRDS(TOM_2018_Spring_SERC, file = "./Data/TOM_2018_Spring_SERC.R")
rm(TOM_2018_Spring_SERC)
rm(adjacency_2018_Spring_SERC)
TOM_2018_Summer_SERC = TOMsimilarity(adjacency_2018_Summer_SERC, TOMType = "signed Nowick")
dissTOM_2018_Summer_SERC = 1 - TOM_2018_Summer_SERC
#saveRDS(TOM_2018_Summer_SERC, file = "./Data/TOM_2018_Summer_SERC.R")
rm(TOM_2018_Summer_SERC)
rm(adjacency_2018_Summer_SERC)

#### Gene Clustering Plot ####
# Code from this section until the end was repeated for each dissTom element.
# For each repeat of dataset change: name of dissTOM, soft power number, TPM_2017_Fall_30, ect.
# Call the hierarchical clustering function
# hclust dendrogram of genes
geneTree = hclust(as.dist(dissTOM_2017_Fall), method = "average")

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12, 9)
jpeg("./Plots/dissTOM_2017_Fall-tree.jpeg",height = 1080, width = 1920)
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
  distM = dissTOM_2017_Fall,
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
jpeg("./Plots/dissTOM_2017_Fall-tree-blocks.jpeg",height = 1080, width = 1920)
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
MEList <- moduleEigengenes(TPM_2017_Fall_30,
                           colors = dynamicColors,
                           softPower = 5)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(7, 6)
jpeg("./Plots/ME-tree_2017_Fall.jpeg",height = 1080, width = 1920)

plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")

MEDissThres <-  0.25 # threshold for merging modules, corresponds to module correlation of 0.75

# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")

# Call an automatic merging function
merge <- mergeCloseModules(TPM_2017_Fall_30,
                           dynamicColors,
                           cutHeight = MEDissThres,
                           verbose = 3)
# The merged module colors
mergedColors <- merge$colors

# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs
dev.off()

sizeGrWindow(12, 9)
pdf(file = "./Plots/geneDendro_2017_Fall.pdf", wi = 9, he = 6)
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

#### TOM plotting ####
#plotTOM <- TOM

# Set diagonal to NA for a nicer plot
#diag(plotTOM) <- NA

# Call the plot function
#sizeGrWindow(9, 9)
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

#### Creating a Module Heatmap of gene Expression ####
Module_Heatmap <- TPM_2017_Fall_30 %>%
  rownames_to_column(var = "sample") %>%
  column_to_rownames(var = "sample")


ModuleGenes_Df <- data.frame(gene = NULL, color = NULL)
for(color in unique(moduleColors)) {
  temp <- data.frame(names(TPM_2017_Fall_30)[moduleColors == color])
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
jpeg("./Plots/ModuleGene_Heatmap_2017_Fall.jpeg",
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

# Bring in growth data

growth.data = read.csv("./Formatted.Data/all.growth.data.csv") %>%
  filter(YEAR == 2017) %>%
  select(SITE, YEAR, TREE_ID, RGR, growth)
growth.data.2 = subset(growth.data, growth.data$RGR != "NA")

sample_sub = subset(Sample_Description, Sample_Description$sample.description %in% growth_MEs$sample)
growth.data.3 = left_join(growth.data.2,sample_sub, by = c("TREE_ID" = "Tree_ID"))

#### Exporting Eigengene values ####
growth_MEs = left_join(growth_MEs, growth.data.3, by = c("sample" = "sample.description"))

growth_MEs <- growth_MEs %>% 
  rename_all(~ str_replace(., "ME","ME_"))
#write_csv(growth_MEs,"./Data/new.MEs/Spring.2017.MEs.csv")


#### Growth Models ####
##### Spring 2017 ######
Spring.2017.MEs = read.csv("./Data/new.MEs/Spring.2017.MEs.csv")

Spring.2017.MEs.2 = Spring.2017.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2017.MEs.3 <- Spring.2017.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.MEs.3 <- Spring.2017.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2017.MEs.3 <- Spring.2017.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### Summer 2017 ######
Summer.2017.MEs = read.csv("./Data/new.MEs/Summer.2017.MEs.csv")

Summer.2017.MEs.2 = Summer.2017.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2017.MEs.3 <- Summer.2017.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2017.MEs.3 <- Summer.2017.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2017.MEs.3 <- Summer.2017.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### Fall 2017 #####
Fall.2017.MEs = read.csv("./Data/new.MEs/Fall.2017.MEs.csv")

Fall.2017.MEs.2 = Fall.2017.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2017.MEs.3 <- Fall.2017.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2017.MEs.3 <- Fall.2017.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2017.MEs.3 <- Fall.2017.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# 11 significant modules:
sig.mod = subset(Fall.2017.MEs.3, Fall.2017.MEs.3$module %in% filtered_Fall.2017.MEs.3$module)
sig.mod$Plot[[20]]
# outlier negative growth value
 
##### Spring 2018 #####
Spring.2018.MEs = read.csv("./Data/new.MEs/Spring.2018.MEs.csv")

Spring.2018.MEs.2 = Spring.2018.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2018.MEs.3 <- Spring.2018.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2018.MEs.3 <- Spring.2018.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2018.MEs.3 <- Spring.2018.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### Summer 2018 #####
Summer.2018.MEs = read.csv("./Data/new.MEs/Summer.2018.MEs.csv")

Summer.2018.MEs.2 = Summer.2018.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2018.MEs.3 <- Summer.2018.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2018.MEs.3 <- Summer.2018.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2018.MEs.3 <- Summer.2018.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### Fall 2018 #####
Fall.2018.MEs = read.csv("./Data/new.MEs/Fall.2018.MEs.csv")

Fall.2018.MEs.2 = Fall.2018.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2018.MEs.3 <- Fall.2018.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2018.MEs.3 <- Fall.2018.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2018.MEs.3 <- Fall.2018.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Spring 2019 #####
Spring.2019.MEs = read.csv("./Data/new.MEs/Spring.2019.MEs.csv")

Spring.2019.MEs.2 = Spring.2019.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2019.MEs.3 <- Spring.2019.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2019.MEs.3 <- Spring.2019.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2019.MEs.3 <- Spring.2019.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Summer 2019 #####
Summer.2019.MEs = read.csv("./Data/new.MEs/Summer.2019.MEs.csv")

Summer.2019.MEs.2 = Summer.2019.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2019.MEs.3 <- Summer.2019.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2019.MEs.3 <- Summer.2019.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2019.MEs.3 <- Summer.2019.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Fall 2019 #####
Fall.2019.MEs = read.csv("./Data/new.MEs/Fall.2019.MEs.csv")

Fall.2019.MEs.2 = Fall.2019.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2019.MEs.3 <- Fall.2019.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2019.MEs.3 <- Fall.2019.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2019.MEs.3 <- Fall.2019.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Spring 2020 #####
Spring.2020.MEs = read.csv("./Data/new.MEs/Spring.2020.MEs.csv")

Spring.2020.MEs.2 = Spring.2020.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2020.MEs.3 <- Spring.2020.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2020.MEs.3 <- Spring.2020.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2020.MEs.3 <- Spring.2020.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Spring HF 2017 #####
Spring.2017.HF.MEs = read.csv("./Data/new.MEs/Spring_HF.2017.MEs.csv")

Spring.2017.HF.MEs.2 = Spring.2017.HF.MEs %>%
  dplyr::select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2017.HF.MEs.3 <- Spring.2017.HF.MEs.2 %>%
  dplyr::mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.HF.MEs.3 <- Spring.2017.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2017.HF.MEs.3 <- Spring.2017.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# 1 significant modules:
sig.mod = subset(Spring.2017.HF.MEs.3, Spring.2017.HF.MEs.3$module %in% filtered_Spring.2017.HF.MEs.3$module)
sig.mod$Plot[[2]]

# outlier growth value
##### Summer HF 2017 #####
Summer.2017.HF.MEs = read.csv("./Data/new.MEs/Summer_HF.2017.MEs.csv")

Summer.2017.HF.MEs.2 = Summer.2017.HF.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2017.HF.MEs.3 <- Summer.2017.HF.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2017.HF.MEs.3 <- Summer.2017.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2017.HF.MEs.3 <- Summer.2017.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Fall HF 2017 #####
Fall.2017.HF.MEs = read.csv("./Data/new.MEs/Fall_HF.2017.MEs.csv")

Fall.2017.HF.MEs.2 = Fall.2017.HF.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2017.HF.MEs.3 <- Fall.2017.HF.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2017.HF.MEs.3 <- Fall.2017.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2017.HF.MEs.3 <- Fall.2017.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Spring SERC 2017 #####
Spring.2017.SERC.MEs = read.csv("./Data/new.MEs/Spring_SERC.2017.MEs.csv")

Spring.2017.SERC.MEs.2 = Spring.2017.SERC.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2017.SERC.MEs.3 <- Spring.2017.SERC.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.SERC.MEs.3 <- Spring.2017.SERC.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2017.SERC.MEs.3 <- Spring.2017.SERC.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Summer SERC 2017 #####
Summer.2017.SERC.MEs = read.csv("./Data/new.MEs/Summer_SERC.2017.MEs.csv")

Summer.2017.SERC.MEs.2 = Summer.2017.SERC.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2017.SERC.MEs.3 <- Summer.2017.SERC.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2017.SERC.MEs.3 <- Summer.2017.SERC.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2017.SERC.MEs.3 <- Summer.2017.SERC.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Fall SERC 2017 #####
Fall.2017.SERC.MEs = read.csv("./Data/new.MEs/Fall_SERC.2017.MEs.csv")

Fall.2017.SERC.MEs.2 = Fall.2017.SERC.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2017.SERC.MEs.3 <- Fall.2017.SERC.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2017.SERC.MEs.3 <- Fall.2017.SERC.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2017.SERC.MEs.3 <- Fall.2017.SERC.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# 1 significant modules:
sig.mod = subset(Fall.2017.SERC.MEs.3, Fall.2017.SERC.MEs.3$module %in% filtered_Fall.2017.SERC.MEs.3$module)
sig.mod$Plot[[2]]
# outlier growth value

##### Spring HF 2018 #####
Spring.2018.HF.MEs = read.csv("./Data/new.MEs/Spring_HF.2018.MEs.csv")

Spring.2018.HF.MEs.2 = Spring.2018.HF.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2018.HF.MEs.3 <- Spring.2018.HF.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2018.HF.MEs.3 <- Spring.2018.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2018.HF.MEs.3 <- Spring.2018.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### Summer HF 2018 #####
Summer.2018.HF.MEs = read.csv("./Data/new.MEs/Summer_HF.2018.MEs.csv")

Summer.2018.HF.MEs.2 = Summer.2018.HF.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2018.HF.MEs.3 <- Summer.2018.HF.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2018.HF.MEs.3 <- Summer.2018.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2018.HF.MEs.3 <- Summer.2018.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Fall HF 2018 #####
Fall.2018.HF.MEs = read.csv("./Data/new.MEs/Fall_HF.2018.MEs.csv")

Fall.2018.HF.MEs.2 = Fall.2018.HF.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2018.HF.MEs.3 <- Fall.2018.HF.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2018.HF.MEs.3 <- Fall.2018.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2018.HF.MEs.3 <- Fall.2018.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Spring SERC 2018 #####
Spring.2018.SERC.MEs = read.csv("./Data/new.MEs/Spring_SERC.2018.MEs.csv")

Spring.2018.SERC.MEs.2 = Spring.2018.SERC.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2018.SERC.MEs.3 <- Spring.2018.SERC.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2018.SERC.MEs.3 <- Spring.2018.SERC.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2018.SERC.MEs.3 <- Spring.2018.SERC.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### Summer SERC 2018 #####
Summer.2018.SERC.MEs = read.csv("./Data/new.MEs/Summer_SERC.2018.MEs.csv")

Summer.2018.SERC.MEs.2 = Summer.2018.SERC.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2018.SERC.MEs.3 <- Summer.2018.SERC.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2018.SERC.MEs.3 <- Summer.2018.SERC.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2018.SERC.MEs.3 <- Summer.2018.SERC.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Fall SERC 2018 #####
Fall.2018.SERC.MEs = read.csv("./Data/new.MEs/Fall_SERC.2018.MEs.csv")

Fall.2018.SERC.MEs.2 = Fall.2018.SERC.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2018.SERC.MEs.3 <- Fall.2018.SERC.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2018.SERC.MEs.3 <- Fall.2018.SERC.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2018.SERC.MEs.3 <- Fall.2018.SERC.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
#### RGR Models ####
##### Spring 2017 ######
Spring.2017.MEs = read.csv("./Data/new.MEs/Spring.2017.MEs.csv")

Spring.2017.MEs.2 = Spring.2017.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2017.MEs.3 <- Spring.2017.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.MEs.3 <- Spring.2017.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2017.MEs.3 <- Spring.2017.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Summer 2017 ######
Summer.2017.MEs = read.csv("./Data/new.MEs/Summer.2017.MEs.csv")

Summer.2017.MEs.2 = Summer.2017.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2017.MEs.3 <- Summer.2017.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2017.MEs.3 <- Summer.2017.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2017.MEs.3 <- Summer.2017.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)





##### Fall 2017 ######
Fall.2017.MEs = read.csv("./Data/new.MEs/Fall.2017.MEs.csv")

Fall.2017.MEs.2 = Fall.2017.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2017.MEs.3 <- Fall.2017.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2017.MEs.3 <- Fall.2017.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2017.MEs.3 <- Fall.2017.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

# 20 significant modules:
sig.mod = subset(Fall.2017.MEs.3, Fall.2017.MEs.3$module %in% filtered_Fall.2017.MEs.3$module)
sig.mod$lm_glance
sig.mod$Plot[[10]]

# significant but outlier negative growth value

##### Spring 2018 ######
Spring.2018.MEs = read.csv("./Data/new.MEs/Spring.2018.MEs.csv")

Spring.2018.MEs.2 = Spring.2018.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2018.MEs.3 <- Spring.2018.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2018.MEs.3 <- Spring.2018.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2018.MEs.3 <- Spring.2018.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Summer 2018 ######
Summer.2018.MEs = read.csv("./Data/new.MEs/Summer.2018.MEs.csv")

Summer.2018.MEs.2 = Summer.2018.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2018.MEs.3 <- Summer.2018.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2018.MEs.3 <- Summer.2018.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2018.MEs.3 <- Summer.2018.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Fall 2018 ######
Fall.2018.MEs = read.csv("./Data/new.MEs/Fall.2018.MEs.csv")

Fall.2018.MEs.2 = Fall.2018.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2018.MEs.3 <- Fall.2018.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2018.MEs.3 <- Fall.2018.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2018.MEs.3 <- Fall.2018.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Spring 2019 ######
Spring.2019.MEs = read.csv("./Data/new.MEs/Spring.2019.MEs.csv")

Spring.2019.MEs.2 = Spring.2019.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2019.MEs.3 <- Spring.2019.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2019.MEs.3 <- Spring.2019.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2019.MEs.3 <- Spring.2019.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Summer 2019 ######
Summer.2019.MEs = read.csv("./Data/new.MEs/Summer.2019.MEs.csv")

Summer.2019.MEs.2 = Summer.2019.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2019.MEs.3 <- Summer.2019.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2019.MEs.3 <- Summer.2019.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2019.MEs.3 <- Summer.2019.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Fall 2019 ######
Fall.2019.MEs = read.csv("./Data/new.MEs/Fall.2019.MEs.csv")

Fall.2019.MEs.2 = Fall.2019.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2019.MEs.3 <- Fall.2019.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2019.MEs.3 <- Fall.2019.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2019.MEs.3 <- Fall.2019.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Spring 2020 ######
Spring.2020.MEs = read.csv("./Data/new.MEs/Spring.2020.MEs.csv")

Spring.2020.MEs.2 = Spring.2020.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2020.MEs.3 <- Spring.2020.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2020.MEs.3 <- Spring.2020.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2020.MEs.3 <- Spring.2020.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Spring HF 2017 #####
Spring.2017.HF.MEs = read.csv("./Data/new.MEs/Spring_HF.2017.MEs.csv")

Spring.2017.HF.MEs.2 = Spring.2017.HF.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2017.HF.MEs.3 <- Spring.2017.HF.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.HF.MEs.3 <- Spring.2017.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2017.HF.MEs.3 <- Spring.2017.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### Summer HF 2017 #####
Summer.2017.HF.MEs = read.csv("./Data/new.MEs/Summer_HF.2017.MEs.csv")

Summer.2017.HF.MEs.2 = Summer.2017.HF.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2017.HF.MEs.3 <- Summer.2017.HF.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2017.HF.MEs.3 <- Summer.2017.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2017.HF.MEs.3 <- Summer.2017.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Fall HF 2017 #####
Fall.2017.HF.MEs = read.csv("./Data/new.MEs/Fall_HF.2017.MEs.csv")

Fall.2017.HF.MEs.2 = Fall.2017.HF.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2017.HF.MEs.3 <- Fall.2017.HF.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2017.HF.MEs.3 <- Fall.2017.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2017.HF.MEs.3 <- Fall.2017.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### Spring HF 2018 #####
Spring.2018.HF.MEs = read.csv("./Data/new.MEs/Spring_HF.2018.MEs.csv")

Spring.2018.HF.MEs.2 = Spring.2018.HF.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Spring.2018.HF.MEs.3 <- Spring.2018.HF.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2018.HF.MEs.3 <- Spring.2018.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Spring.2018.HF.MEs.3 <- Spring.2018.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

##### Summer HF 2018 #####
Summer.2018.HF.MEs = read.csv("./Data/new.MEs/Summer_HF.2018.MEs.csv")

Summer.2018.HF.MEs.2 = Summer.2018.HF.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Summer.2018.HF.MEs.3 <- Summer.2018.HF.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2018.HF.MEs.3 <- Summer.2018.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Summer.2018.HF.MEs.3 <- Summer.2018.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
##### Fall HF 2018 #####
Fall.2018.HF.MEs = read.csv("./Data/new.MEs/Fall_HF.2018.MEs.csv")

Fall.2018.HF.MEs.2 = Fall.2018.HF.MEs %>%
  select(-Site,-Year) %>%
  pivot_longer(starts_with("ME_"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,SITE,YEAR,Season,TREE_ID,RGR,growth,eigen_value))

Fall.2018.HF.MEs.3 <- Fall.2018.HF.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2018.HF.MEs.3 <- Fall.2018.HF.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

filtered_Fall.2018.HF.MEs.3 <- Fall.2018.HF.MEs.3 %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)


