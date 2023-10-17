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
TPM_2017_Fall_HF = TPM_2017_Fall[,c()]
TPM_2017_Fall_SREC = TPM_2017_Fall[,c()]

# split TPMS for each year including all seasons
TPM_2017 = TPM_ExprData[,c(1:59,108:166)] # 117 samples
TPM_2018 = TPM_ExprData[,c(1,60:107,167:223)] # 105 samples
TPM_2019 = TPM_ExprData[,c(1,224:276)] # 53 samples

# split TPMS for each year including all seasons by site









cv <- cv.spls( X, Y.2, eta = seq(0.1,0.9,0.1), K = c(5:10) )
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



