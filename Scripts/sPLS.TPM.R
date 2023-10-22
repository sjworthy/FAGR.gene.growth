#### Sparse Partial Least Squares using TPMs ####

library(tidyverse)
library(pls)
library(caTools)

# https://cran.r-project.org/web/packages/spls/vignettes/spls-example.pdf
# The main principle of this methodology is to impose sparsity within the context 
# of partial least squares and thereby carry out dimension reduction and variable 
# selection simultaneously.
# As part of pre-processing, the predictors are centered and scaled and the 
# responses are centered automatically as default by the package ‘spls’.

#### Read in Data ####
TPM_ExprData = read.csv("./Raw.Data/HTseq-master-counts.csv")
Growth_data = read.csv("./Formatted.Data/all.growth.data.csv")
Sample_Description = read.csv("./Formatted.Data/FAGR.description.csv")

#### Subset Growth Data ####

# Remove samples with NA for growth or RGR

Growth_data_2 = subset(Growth_data, Growth_data$RGR != "NA")

growth_2017 = Growth_data_2 %>%
  filter(YEAR == 2017) %>%
  droplevels(.)

growth_2018 = Growth_data_2 %>%
  filter(YEAR == 2018) %>%
  droplevels(.)
  
growth_2019 = Growth_data_2 %>%
  filter(YEAR == 2019) %>%
  droplevels(.)

growth_2020 = Growth_data_2 %>%
  filter(YEAR == 2020) %>%
  droplevels(.)

# get sample for each tree ID, to be used to subset expression data to match growth data

growth_sub_2017 = subset(Sample_Description, Sample_Description$Tree_ID %in% growth_2017$TREE_ID)
growth_sub_2018 = subset(Sample_Description, Sample_Description$Tree_ID %in% growth_2018$TREE_ID)
growth_sub_2019 = subset(Sample_Description, Sample_Description$Tree_ID %in% growth_2019$TREE_ID)
growth_sub_2020 = subset(Sample_Description, Sample_Description$Tree_ID %in% growth_2020$TREE_ID)

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
TPM_2017_Fall_SREC = TPM_2017_Fall[,c(1,21:40)]
TPM_2018_Fall_HF = TPM_2018_Fall[,c(1:18)]
TPM_2018_Fall_SREC = TPM_2018_Fall[,c(1,19:35)]

# split TPMS for each year including all seasons
TPM_2017 = TPM_ExprData[,c(1:16,18:34,36:54,56:59,108:166)] # 114 samples
TPM_2018 = TPM_ExprData[,c(1,60:101,103:107,167:168,170:188,190:223)] # 102 samples
TPM_2019 = TPM_ExprData[,c(1,224:230,232:248,250:276)] # 51 samples

# split TPMS for each year including all seasons by site

TPM_2017_HF = TPM_2017[,c(1:56)]
TPM_2017_SERC = TPM_2017[,c(1,57:115)]
TPM_2018_HF = TPM_2018[,c(1:48)]
TPM_2018_SERC = TPM_2018[,c(1,49:103)]

#### Calculate the Variance of each Gene across all of our samples ####
# same process done in WGCNA analysis of TPM

# Calculating Coefficient of variation function
calc.cv <- function(x, na.rm=TRUE) { 
  if(na.rm==TRUE) x <- na.omit(x)
  result <- sd(x) / mean(x) 
  result <- abs(result) 
  return(result)
}

TPM_2017_Spring <- TPM_2017_Spring %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Summer <- TPM_2017_Summer %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Fall <- TPM_2017_Fall %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Spring <- TPM_2018_Spring %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Summer <- TPM_2018_Summer %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Fall <- TPM_2018_Fall %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2019_Spring <- TPM_2019_Spring %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2019_Summer <- TPM_2019_Summer %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2019_Fall <- TPM_2019_Fall %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2020_Spring <- TPM_2020_Spring %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

TPM_2017_Spring_HF <- TPM_2017_Spring_HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Spring_SERC <- TPM_2017_Spring_SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Spring_HF <- TPM_2018_Spring_HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Spring_SERC <- TPM_2018_Spring_SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Summer_HF <- TPM_2017_Spring_HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Summer_SERC <- TPM_2017_Spring_SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Summer_HF <- TPM_2018_Spring_HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Summer_SERC <- TPM_2018_Spring_SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Fall_HF <- TPM_2017_Spring_HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_Fall_SERC <- TPM_2017_Spring_SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Fall_HF <- TPM_2018_Spring_HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_Fall_SERC <- TPM_2018_Spring_SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

TPM_2017 <- TPM_2017 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018 <- TPM_2018 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2019 <- TPM_2019 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

TPM_2017_HF <- TPM_2017_HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2017_SERC <- TPM_2017_SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_HF <- TPM_2018_HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TPM_2018_SERC <- TPM_2018_SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

# Filter the data frame to contain the top 30% most variable genes
# ties are kept together
TPM_2017_Spring_30  <- TPM_2017_Spring  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Summer_30  <- TPM_2017_Summer  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Fall_30  <- TPM_2017_Fall  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Spring_30  <- TPM_2018_Spring  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Summer_30  <- TPM_2018_Summer  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Fall_30  <- TPM_2018_Fall  %>% slice_max(order_by = cv, prop = .30)
TPM_2019_Spring_30  <- TPM_2019_Spring  %>% slice_max(order_by = cv, prop = .30)
TPM_2019_Summer_30  <- TPM_2019_Summer  %>% slice_max(order_by = cv, prop = .30)
TPM_2019_Fall_30  <- TPM_2019_Fall  %>% slice_max(order_by = cv, prop = .30)
TPM_2020_Spring_30  <- TPM_2020_Spring  %>% slice_max(order_by = cv, prop = .30)

TPM_2017_Spring_HF_30  <- TPM_2017_Spring_HF  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Spring_SERC_30  <- TPM_2017_Spring_SERC  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Spring_HF_30  <- TPM_2018_Spring_HF  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Spring_SERC_30  <- TPM_2018_Spring_SERC  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Summer_HF_30  <- TPM_2017_Summer_HF  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Summer_SERC_30  <- TPM_2017_Summer_SERC  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Summer_HF_30  <- TPM_2018_Summer_HF  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Summer_SERC_30  <- TPM_2018_Summer_SERC  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Fall_HF_30  <- TPM_2017_Fall_HF  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_Fall_SERC_30  <- TPM_2017_Fall_SERC  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Fall_HF_30  <- TPM_2018_Fall_HF  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_Fall_SERC_30  <- TPM_2018_Fall_SERC  %>% slice_max(order_by = cv, prop = .30)

TPM_2017_30 <- TPM_2017  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_30 <- TPM_2018  %>% slice_max(order_by = cv, prop = .30)
TPM_2019_30 <- TPM_2019  %>% slice_max(order_by = cv, prop = .30)

TPM_2017_HF_30 <- TPM_2017_HF  %>% slice_max(order_by = cv, prop = .30)
TPM_2017_SERC_30 <- TPM_2017_SERC  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_HF_30 <- TPM_2018_HF  %>% slice_max(order_by = cv, prop = .30)
TPM_2018_SERC_30 <- TPM_2018_SERC  %>% slice_max(order_by = cv, prop = .30)

# Deselect the cv column and flip our data frame to contain sample along the rows and genes along the columns
TPM_2017_Spring_30 <- select(TPM_2017_Spring_30, -cv)
TPM_2017_Spring_30 <- column_to_rownames(TPM_2017_Spring_30,var = "Gene_ID")
TPM_2017_Spring_30 <- as.matrix(t(TPM_2017_Spring_30))

TPM_2017_Summer_30 <- select(TPM_2017_Summer_30, -cv)
TPM_2017_Summer_30 <- column_to_rownames(TPM_2017_Summer_30,var = "Gene_ID")
TPM_2017_Summer_30 <- as.matrix(t(TPM_2017_Summer_30))

TPM_2017_Fall_30 <- select(TPM_2017_Fall_30, -cv)
TPM_2017_Fall_30 <- column_to_rownames(TPM_2017_Fall_30,var = "Gene_ID")
TPM_2017_Fall_30 <- as.matrix(t(TPM_2017_Fall_30))

TPM_2018_Spring_30 <- select(TPM_2018_Spring_30, -cv)
TPM_2018_Spring_30 <- column_to_rownames(TPM_2018_Spring_30,var = "Gene_ID")
TPM_2018_Spring_30 <- as.matrix(t(TPM_2018_Spring_30))

TPM_2018_Summer_30 <- select(TPM_2018_Summer_30, -cv)
TPM_2018_Summer_30 <- column_to_rownames(TPM_2018_Summer_30,var = "Gene_ID")
TPM_2018_Summer_30 <- as.matrix(t(TPM_2018_Summer_30))

TPM_2018_Fall_30 <- select(TPM_2018_Fall_30, -cv)
TPM_2018_Fall_30 <- column_to_rownames(TPM_2018_Fall_30,var = "Gene_ID")
TPM_2018_Fall_30 <- as.matrix(t(TPM_2018_Fall_30))

TPM_2019_Spring_30 <- select(TPM_2019_Spring_30, -cv)
TPM_2019_Spring_30 <- column_to_rownames(TPM_2019_Spring_30,var = "Gene_ID")
TPM_2019_Spring_30 <- as.matrix(t(TPM_2019_Spring_30))

TPM_2019_Summer_30 <- select(TPM_2019_Summer_30, -cv)
TPM_2019_Summer_30 <- column_to_rownames(TPM_2019_Summer_30,var = "Gene_ID")
TPM_2019_Summer_30 <- as.matrix(t(TPM_2019_Summer_30))

TPM_2019_Fall_30 <- select(TPM_2019_Fall_30, -cv)
TPM_2019_Fall_30 <- column_to_rownames(TPM_2019_Fall_30,var = "Gene_ID")
TPM_2019_Fall_30 <- as.matrix(t(TPM_2019_Fall_30))

TPM_2020_Spring_30 <- select(TPM_2020_Spring_30, -cv)
TPM_2020_Spring_30 <- column_to_rownames(TPM_2020_Spring_30,var = "Gene_ID")
TPM_2020_Spring_30 <- as.matrix(t(TPM_2020_Spring_30))

TPM_2017_Spring_HF_30 <- select(TPM_2017_Spring_HF_30, -cv)
TPM_2017_Spring_HF_30 <- column_to_rownames(TPM_2017_Spring_HF_30,var = "Gene_ID")
TPM_2017_Spring_HF_30 <- as.matrix(t(TPM_2017_Spring_HF_30))

TPM_2017_Spring_SERC_30 <- select(TPM_2017_Spring_SERC_30, -cv)
TPM_2017_Spring_SERC_30 <- column_to_rownames(TPM_2017_Spring_SERC_30,var = "Gene_ID")
TPM_2017_Spring_SERC_30 <- as.matrix(t(TPM_2017_Spring_SERC_30))

TPM_2018_Spring_HF_30 <- select(TPM_2018_Spring_HF_30, -cv)
TPM_2018_Spring_HF_30 <- column_to_rownames(TPM_2018_Spring_HF_30,var = "Gene_ID")
TPM_2018_Spring_HF_30 <- as.matrix(t(TPM_2018_Spring_HF_30))

TPM_2018_Spring_SERC_30 <- select(TPM_2018_Spring_SERC_30, -cv)
TPM_2018_Spring_SERC_30 <- column_to_rownames(TPM_2018_Spring_SERC_30,var = "Gene_ID")
TPM_2018_Spring_SERC_30 <- as.matrix(t(TPM_2018_Spring_SERC_30))

TPM_2017_Summer_HF_30 <- select(TPM_2017_Summer_HF_30, -cv)
TPM_2017_Summer_HF_30 <- column_to_rownames(TPM_2017_Summer_HF_30,var = "Gene_ID")
TPM_2017_Summer_HF_30 <- as.matrix(t(TPM_2017_Summer_HF_30))

TPM_2017_Summer_SERC_30 <- select(TPM_2017_Summer_SERC_30, -cv)
TPM_2017_Summer_SERC_30 <- column_to_rownames(TPM_2017_Summer_SERC_30,var = "Gene_ID")
TPM_2017_Summer_SERC_30 <- as.matrix(t(TPM_2017_Summer_SERC_30))

TPM_2018_Summer_HF_30 <- select(TPM_2018_Summer_HF_30, -cv)
TPM_2018_Summer_HF_30 <- column_to_rownames(TPM_2018_Summer_HF_30,var = "Gene_ID")
TPM_2018_Summer_HF_30 <- as.matrix(t(TPM_2018_Summer_HF_30))

TPM_2018_Summer_SERC_30 <- select(TPM_2018_Summer_SERC_30, -cv)
TPM_2018_Summer_SERC_30 <- column_to_rownames(TPM_2018_Summer_SERC_30,var = "Gene_ID")
TPM_2018_Summer_SERC_30 <- as.matrix(t(TPM_2018_Summer_SERC_30))

TPM_2017_Fall_HF_30 <- select(TPM_2017_Fall_HF_30, -cv)
TPM_2017_Fall_HF_30 <- column_to_rownames(TPM_2017_Fall_HF_30,var = "Gene_ID")
TPM_2017_Fall_HF_30 <- as.matrix(t(TPM_2017_Fall_HF_30))

TPM_2017_Fall_SERC_30 <- select(TPM_2017_Fall_SERC_30, -cv)
TPM_2017_Fall_SERC_30 <- column_to_rownames(TPM_2017_Fall_SERC_30,var = "Gene_ID")
TPM_2017_Fall_SERC_30 <- as.matrix(t(TPM_2017_Fall_SERC_30))

TPM_2018_Fall_HF_30 <- select(TPM_2018_Fall_HF_30, -cv)
TPM_2018_Fall_HF_30 <- column_to_rownames(TPM_2018_Fall_HF_30,var = "Gene_ID")
TPM_2018_Fall_HF_30 <- as.matrix(t(TPM_2018_Fall_HF_30))

TPM_2018_Fall_SERC_30 <- select(TPM_2018_Fall_SERC_30, -cv)
TPM_2018_Fall_SERC_30 <- column_to_rownames(TPM_2018_Fall_SERC_30,var = "Gene_ID")
TPM_2018_Fall_SERC_30 <- as.matrix(t(TPM_2018_Fall_SERC_30))

TPM_2017_30 <- select(TPM_2017_30, -cv)
TPM_2017_30 <- column_to_rownames(TPM_2017_30,var = "Gene_ID")
TPM_2017_30 <- as.matrix(t(TPM_2017_30))

TPM_2018_30 <- select(TPM_2018_30, -cv)
TPM_2018_30 <- column_to_rownames(TPM_2018_30,var = "Gene_ID")
TPM_2018_30 <- as.matrix(t(TPM_2018_30))

TPM_2019_30 <- select(TPM_2019_30, -cv)
TPM_2019_30 <- column_to_rownames(TPM_2019_30,var = "Gene_ID")
TPM_2019_30 <- as.matrix(t(TPM_2019_30))

TPM_2017_HF_30 <- select(TPM_2017_HF_30, -cv)
TPM_2017_HF_30 <- column_to_rownames(TPM_2017_HF_30,var = "Gene_ID")
TPM_2017_HF_30 <- as.matrix(t(TPM_2017_HF_30))

TPM_2017_SERC_30 <- select(TPM_2017_SERC_30, -cv)
TPM_2017_SERC_30 <- column_to_rownames(TPM_2017_SERC_30,var = "Gene_ID")
TPM_2017_SERC_30 <- as.matrix(t(TPM_2017_SERC_30))

TPM_2018_HF_30 <- select(TPM_2018_HF_30, -cv)
TPM_2018_HF_30 <- column_to_rownames(TPM_2018_HF_30,var = "Gene_ID")
TPM_2018_HF_30 <- as.matrix(t(TPM_2018_HF_30))

TPM_2018_SERC_30 <- select(TPM_2018_SERC_30, -cv)
TPM_2018_SERC_30 <- column_to_rownames(TPM_2018_SERC_30,var = "Gene_ID")
TPM_2018_SERC_30 <- as.matrix(t(TPM_2018_SERC_30))

#### Analyses ####

# How well does 2017 spring gene expression predict 2017 growth and RGR
# TPM_2017_Spring_30, growth_2017 (39)

# order growth data to match TPM order
# subset Sample_Description by names in TPM

sample_sub = subset(Sample_Description, Sample_Description$sample.description %in% rownames(TPM_2017_Spring_30))
growth_2017_2 = growth_2017[order(match(growth_2017$TREE_ID, sample_sub$Tree_ID)),]

# merge expression and trait data
all.dat.2017.Spring = cbind(TPM_2017_Spring_30,growth_2017_2)

# split data into train and test
split = sample.split(Y = all.dat.2017.Spring$growth, SplitRatio = 0.7)
train = all.dat.2017.Spring[split,]
test = all.dat.2017.Spring[!split,]

# subset expression data and make a matrix
train_exp = as.matrix(train[,1:9913])

# make a dataframe for model fitting
model.df = train[,c(9933,9949)]
model.df$Z = as.numeric(train_exp)
model.df$growth.scale = scale(log(model.df$growth))

growth_Spring_2017 <- plsr(growth.scale ~ Z, data = model.df, validation = "LOO", scale = FALSE)
RGR_Spring_2017 <- plsr(RGR ~ Z, data = model.df, validation = "LOO", scale = FALSE)

summary(growth_Spring_2017)
plot(RMSEP(growth_Spring_2017), legendpos = "topright") # only going up 
validationplot(growth_Spring_2017) # only going up 
validationplot(growth_Spring_2017, val.type="MSEP") # only going up 
validationplot(growth_Spring_2017, val.type="R2") # only going down 

output = RMSEP(growth_Spring_2017, estimate = "CV")
plot(output)
ncomp.permut <- selectNcomp(growth_Spring_2017, method = "randomization", plot = TRUE)
ncomp.permut.rgr <- selectNcomp(RGR_Spring_2017, method = "randomization", plot = TRUE)







#### Unused ####
# determine K
min(ncol(TPM_2017_Spring_30),0.9*nrow(growth_2017_2)) #35.1

cv_2017_Spring <- cv.spls(x = TPM_2017_Spring_30, y = growth_2017_2$growth, eta = seq(0.1,0.9,0.1), K = c(1:34),
                          scale.x = FALSE)
test = spls::spls(x = TPM_2017_Spring_30, y = growth_2017_2$growth, K = cv_2017_Spring$K.opt, eta = cv_2017_Spring$eta.opt)
test.2=coef(test)
plot.spls(test, yvar=1 )


pls.result <- pls(TPM_2017_Spring_30, growth_2017_2$growth, mode = "regression", multilevel = NULL, all.outputs = TRUE,
                  scale = TRUE) # multivariate
pls.result <- pls(X, Y.2) #univariate
plotIndiv(pls.result)
plotVar(pls.result)

spls.result <- spls(TPM_2017_Spring_30, growth_2017_2$growth, mode = "regression", scale = TRUE, multilevel = NULL, all.outputs = TRUE)

test <- plsr(growth_2017_2$growth ~ TPM_2017_Spring_30, validation = "LOO")
summary(test)
plot(RMSEP(test), legendpos = "topright")
ncomp.permut <- selectNcomp(test, method = "randomization", plot = TRUE)

## install mixOmics 
BiocManager::install('mixOmics')
library(mixOmics)

# Read in our TPM normalized RNA-seq expression data to R
# This is just the counts data that has been normalized previously
TPM_ExprData <- read.csv("./Raw.Data/HTseq-master-counts.csv", row.names = 1)


# Read in growth data
# samples as rownames, growth variables as column names
growth.data = read.csv("./Formatted.Data/Dendro_FAGR.csv")


# needs to be gene_ID as columns and samples as rows
# should be the same time 30% most variable to match WGCNA?

#### Example Code ####
#http://mixomics.org/methods/spls/
  
# PLS
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
Y.2 = Y$Cholesterol.mg.dL.
# may need to add scale so can compare among data sets
pls.result <- pls(X, Y.2, mode = "regression", multilevel = NULL, all.outputs = TRUE,
                  scale = TRUE) # multivariate
pls.result <- pls(X, Y.2) #univariate
plotIndiv(pls.result)
plotVar(pls.result)





# sPLS
spls.result <- spls(X, Y, keepX = c(10, 20), keepY = c(3, 2))  # run the method
plotIndiv(spls.result) # plot the samples
plotVar(spls.result)   # plot the variables
selectVar(spls.result, comp = 1)$X$name # extract the variables used to construct the first latent component
plotLoadings(spls.result, method = 'mean', contrib = 'max') # depict weight assigned to each of these variables


#### sPLS ####

install.packages("spls")
library(spls)

R> set.seed(1)
cv <- cv.spls( X, Y.2, eta = seq(0.1,0.9,0.1), K = c(5:10) )
test = spls::spls(X,Y.2, K = cv$K.opt, eta = cv$eta.opt)
test.2=coef(test)
plot.spls(test, yvar=1 )

### pls ####

install.packages("pls")
library(pls)
data(gasoline)

gasTrain <- gasoline[1:50,]
gasTest <- gasoline[51:60,]

gas1 <- plsr(octane ~ NIR, ncomp = 10, data = gasTrain, validation = "LOO")



