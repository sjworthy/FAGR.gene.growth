#### sPLS using TPMs ####

install.packages("mixOmics")

# Read in our TPM normalized RNA-seq expression data to R
# This is just the counts data that has been normalized previously
TPM_ExprData <- read.csv("./Raw.Data/HTseq-master-counts.csv", row.names = 1)

# Read in growth data

growth.data = read.csv("./Formatted.Data/Dendro_FAGR.csv")

