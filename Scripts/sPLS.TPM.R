#### sPLS using TPMs ####

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



