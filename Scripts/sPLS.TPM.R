#### Sparse Partial Least Squares using TPMs ####

library(spls)
# https://cran.r-project.org/web/packages/spls/vignettes/spls-example.pdf
# The main principle of this methodology is to impose sparsity within the context 
# of partial least squares and thereby carry out dimension reduction and variable 
# selection simultaneously.
# As part of pre-processing, the predictors are centered and scaled and the 
# responses are centered automatically as default by the package ‘spls’.

#### Read in Data ####
TPM_ExprData = read_csv("./Raw.Data/HTseq-master-counts.csv")
Growth_data = read_csv("./Formatted.Data/all.growth.data.csv")

#### Subset Growth Data ####

# Remove samples with NA for growth or RGR

Growth_data_2 = subset(Growth_data, Growth_data$RGR != "NA")

growth_2017 = Growth_data %>%
  filter(YEAR == 2017) %>%
  filter(growth !(NA))
  droplevels(.)


#### split samples into time points ####

# split TPMS for each year by each season
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

# split TPMS for each year by each season by each site
# only for 2017 and 2018 where we have both sites

TPM_2017_Spring_HF = TPM_2017_Spring[,c(1:21)]
TPM_2017_Spring_SERC = TPM_2017_Spring[,c(1,22:41)]
TPM_2018_Spring_HF = TPM_2018_Spring[,c(1:14)]
TPM_2018_Spring_SERC = TPM_2018_Spring[,c(1,15:34)]
TPM_2017_Summer_HF = TPM_2017_Summer[,c(1:19)]
TPM_2017_Summer_SERC = TPM_2017_Summer[,c(1,20:38)]
TPM_2018_Summer_HF = TPM_2018_Summer[,c(1:18)]
TPM_2018_Summer_SERC = TPM_2018_Summer[,c(1,19:38)]
TPM_2017_Fall_HF = TPM_2017_Fall[,c(1:21)]
TPM_2017_Fall_SREC = TPM_2017_Fall[,c(1,22:41)]
TPM_2018_Fall_HF = TPM_2018_Fall[,c(1:19)]
TPM_2018_Fall_SREC = TPM_2018_Fall[,c(1,20:36)]

# split TPMS for each year including all seasons
TPM_2017 = TPM_ExprData[,c(1:59,108:166)] # 117 samples
TPM_2018 = TPM_ExprData[,c(1,60:107,167:223)] # 105 samples
TPM_2019 = TPM_ExprData[,c(1,224:276)] # 53 samples

# split TPMS for each year including all seasons by site

TPM_2017_HF = TPM_2017[,c(1:59)]
TPM_2017_SERC = TPM_2017[,c(1,60:118)]
TPM_2018_HF = TPM_2018[,c(1:49)]
TPM_2018_SERC = TPM_2018[,c(1,50:106)]

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
TPM_2017_30 <- as.data.frame(t(TPM_2017_30))

TPM_2018_30 <- select(TPM_2018_30, -cv)
TPM_2018_30 <- column_to_rownames(TPM_2018_30,var = "Gene_ID")
TPM_2018_30 <- as.data.frame(t(TPM_2018_30))

TPM_2019_30 <- select(TPM_2019_30, -cv)
TPM_2019_30 <- column_to_rownames(TPM_2019_30,var = "Gene_ID")
TPM_2019_30 <- as.data.frame(t(TPM_2019_30))

TPM_2017_HF_30 <- select(TPM_2017_HF_30, -cv)
TPM_2017_HF_30 <- column_to_rownames(TPM_2017_HF_30,var = "Gene_ID")
TPM_2017_HF_30 <- as.data.frame(t(TPM_2017_HF_30))

TPM_2017_SERC_30 <- select(TPM_2017_SERC_30, -cv)
TPM_2017_SERC_30 <- column_to_rownames(TPM_2017_SERC_30,var = "Gene_ID")
TPM_2017_SERC_30 <- as.data.frame(t(TPM_2017_SERC_30))

TPM_2018_HF_30 <- select(TPM_2018_HF_30, -cv)
TPM_2018_HF_30 <- column_to_rownames(TPM_2018_HF_30,var = "Gene_ID")
TPM_2018_HF_30 <- as.data.frame(t(TPM_2018_HF_30))

TPM_2018_SERC_30 <- select(TPM_2018_SERC_30, -cv)
TPM_2018_SERC_30 <- column_to_rownames(TPM_2018_SERC_30,var = "Gene_ID")
TPM_2018_SERC_30 <- as.data.frame(t(TPM_2018_SERC_30))

#### Analyses ####
# first subset the growth data to match each TPM
# second tune the parameters using cv.spls(), 10-fold cross-validation
# eta = sparcity tuning parameter, value between 0 and 1
# K = number of hidden (latent) variables, range from 1 to min {p,(v − 1)n/v}, 
# where p is the number of predictors and n is the sample size
# for 10-fold, min{p,0.9n}

# TPM_2017_Spring_30
growth_2017_Spring = Growth_data %>%
  filter(YEAR == 2017) %>%
  droplevels(.)

# determine K
min(9912,0.9*40) #36

TPM_2017_Spring
x=genes, y=grwoth

cv <- cv.spls(x = TPM_2017_Spring_30, y = growth_2017_Spring$growth, eta = seq(0.1,0.9,0.1), K = c(5:10) )
test = spls::spls(X,Y.2, K = cv$K.opt, eta = cv$eta.opt)
test.2=coef(test)
plot.spls(test, yvar=1 )




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



