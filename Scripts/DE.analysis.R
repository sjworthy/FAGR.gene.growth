#### Differential Expression Analysis ####

# comparisons to make
# 2017-2018, 2018-2019, 2019-2020, 2017-2019
# need to make sure have matching trees in each set for comparison

#Libraries
library(tidyverse)
library(magrittr)
library(randomcoloR)
library(edgeR)

#### Read in our TPM normalized RNA-seq expression data ####

TPM_ExprData <- read_csv("./Raw.Data/HTseq-master-counts.csv")
# 295 samples

#### Expression Filtering ####
# Remove low expression genes
# Need at least 10 reads in at least 3 samples

Gene_Counts <- TPM_ExprData[rowSums(TPM_ExprData[,-1] > 10) >= 3,]

#### Filtering Samples ####
#read in sample descriptions

# need to make sure we have two time points for each sample with matching season
# all samples from both years with matching season time points
time.2017.18 = Gene_Counts[,c(1:223)] # 222 samples
time.2017.18 = time.2017.18[,c(1:7,10,12,13,15,17,18,20:22,24:40,42:45,47:107,
                              108:148,150,152:154,156:198,200:223)] # 209 samples, both pops
time.2018.19 = Gene_Counts[,c(1,60:107,167:223,224:276)] # 158 samples
time.2018.19 = time.2018.19[,c(1,50:52,54:62,64:67,69:80,82:95,97:144,146:159)] # 104 samples, only SERC
time.2019.20 = Gene_Counts[,c(1,224:276,277:296)] # 73 samples
time.2019.20 = time.2019.20[,c(1:18,55:57,59:67,69:72,74)] # 34 samples, only SERC, only Spring
time.2017.19 = Gene_Counts[,c(1:59,108:166,224:276)] # 170 samples
time.2017.19 = time.2017.19[,c(1,60:62,64:72,74:77,79:90,92:102,104:106,109:146,148:171)] # 104 samples, only SERC

# may be interesting to compare DE among season to seasonal growth
# growth file only has yearly comparisons
# maybe split two populations?

#### Sample Description ####

Sample_Description = read_csv("./Formatted.Data/FAGR.description.csv")
# convert year and site to factor as it will be used to group samples
Sample_Description$Year = as.factor(Sample_Description$Year)
Sample_Description$Site = as.factor(Sample_Description$Site)
Sample_Description = Sample_Description %>%
  unite(group,c(Site,Year),sep = "_", remove = FALSE)

Sample_Description.2017.2018 = Sample_Description %>%
  filter(sample.description %in% colnames(time.2017.18)) %>%
  droplevels(.)
Sample_Description.2018.2019 = Sample_Description %>%
  filter(sample.description %in% colnames(time.2018.19)) %>%
  droplevels(.)
Sample_Description.2017.2019 = Sample_Description %>%
  filter(sample.description %in% colnames(time.2017.19)) %>%
  droplevels(.)
Sample_Description.2019.2020 = Sample_Description %>%
  filter(sample.description %in% colnames(time.2019.20)) %>%
  droplevels(.)

#### Create DGE data ####
# Create a numeric matrix of our count data to serve as an input to DGElist.
Data.matrix.2017.18 <- time.2017.18 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.2017.18) <- time.2017.18$Gene_ID

Data.matrix.2018.19 <- time.2018.19 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.2018.19) <- time.2018.19$Gene_ID

Data.matrix.2017.19 <- time.2017.19 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.2017.19) <- time.2017.19$Gene_ID

Data.matrix.2019.20 <- time.2019.20 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.2019.20) <- time.2019.20$Gene_ID

#The DGEList function creates a DGElist object for differential expression analysis. 
# we specify counts to be the data matrix we created above that holds our count data. 
# Group is specified from the Year in Sample Description we created above which contains information about the species and treatment of each sample.
DGE.data.2017.18 = DGEList(counts = Data.matrix.2017.18, group = Sample_Description.2017.2018$Year)
DGE.data.2018.19 = DGEList(counts = Data.matrix.2018.19, group = Sample_Description.2018.2019$Year)
DGE.data.2017.19 = DGEList(counts = Data.matrix.2017.19, group = Sample_Description.2017.2019$Year)
DGE.data.2019.20 = DGEList(counts = Data.matrix.2019.20, group = Sample_Description.2019.2020$Year)

# group as site and year
DGE.data.2017.18.site = DGEList(counts = Data.matrix.2017.18, group = Sample_Description.2017.2018$group)

#We that normalize the data in our DGElist object using the TMM method
DGE.data.2017.18 <- calcNormFactors(DGE.data.2017.18, method = "TMM")
DGE.data.2018.19 <- calcNormFactors(DGE.data.2018.19, method = "TMM")
DGE.data.2017.19 <- calcNormFactors(DGE.data.2017.19, method = "TMM")
DGE.data.2019.20 <- calcNormFactors(DGE.data.2019.20, method = "TMM")

DGE.data.2017.18.site <- calcNormFactors(DGE.data.2017.18.site, method = "TMM")

#### Plotting ####
# bcv is the square root of the dispersion of the negative binomial distribution
plotMDS(DGE.data.2017.18, method = "bcv") 
plotMDS(DGE.data.2018.19, method = "bcv") 
plotMDS(DGE.data.2017.18, method = "bcv") 
plotMDS(DGE.data.2019.20, method = "bcv") 
plotMDS(DGE.data.2017.18.site, method = "bcv") 


#### TMM Data Saving ####
TMM_DATA_2017.18 <- cpm(DGE.data.2017.18, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2017.18,"./Data/DE.data/TMM_NormData_LogCPM_2017_2018.csv")
TMM_DATA_2018.19 <- cpm(DGE.data.2018.19, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2018.19,"./Data/DE.data/TMM_NormData_LogCPM_2018_2019.csv")
TMM_DATA_2017.19 <- cpm(DGE.data.2017.19, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2017.19,"./Data/DE.data/TMM_NormData_LogCPM_2017_2019.csv")
TMM_DATA_2019.20 <- cpm(DGE.data.2019.20, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2019.20,"./Data/DE.data/TMM_NormData_LogCPM_2019_2020.csv")
TMM_DATA_2017.18.site <- cpm(DGE.data.2017.18.site, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2017.18.site,"./Data/DE.data/TMM_NormData_LogCPM_2017_2018.site.csv")

#Create our design matrix which considers every interaction between species and treatment.
Design.2017.18 <-
  model.matrix(~ Year, data = Sample_Description.2017.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2017.18) <- Sample_Description.2017.2018$sample.description

Design.2018.19 <-
  model.matrix(~ Year, data = Sample_Description.2018.2019) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2018.19) <- Sample_Description.2018.2019$sample.description

Design.2017.19 <-
  model.matrix(~ Year, data = Sample_Description.2017.2019) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2017.19) <- Sample_Description.2017.2019$sample.description

Design.2019.20 <-
  model.matrix(~ Year, data = Sample_Description.2019.2020) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2019.20) <- Sample_Description.2019.2020$sample.description

Design.2017.18.site <-
  model.matrix(~ Year*Site, data = Sample_Description.2017.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2017.18.site) <- Sample_Description.2017.2018$sample.description

#### Calculate Dispersion Factors ####
# To estimate common dispersion:
DGE.data.2017.18 <- estimateGLMCommonDisp(DGE.data.2017.18, Design.2017.18, verbose = TRUE)
# Disp = 0.89034, BCV = 0.9436
#To estimate trended dispersions:
DGE.data.2017.18 <- estimateGLMTrendedDisp(DGE.data.2017.18, Design.2017.18)
#To estimate tagwise dispersions:
DGE.data.2017.18 <- estimateGLMTagwiseDisp(DGE.data.2017.18, Design.2017.18)
plotBCV(DGE.data.2017.18)

DGE.data.2018.19 <- estimateGLMCommonDisp(DGE.data.2018.19, Design.2018.19, verbose = TRUE)
# Disp = 0.89034, BCV = 0.812
#To estimate trended dispersions:
DGE.data.2018.19 <- estimateGLMTrendedDisp(DGE.data.2018.19, Design.2018.19)
#To estimate tagwise dispersions:
DGE.data.2018.19 <- estimateGLMTagwiseDisp(DGE.data.2018.19, Design.2018.19)
plotBCV(DGE.data.2018.19)

DGE.data.2017.19 <- estimateGLMCommonDisp(DGE.data.2017.19, Design.2017.19, verbose = TRUE)
# Disp = 0.48847, BCV = 0.6989
#To estimate trended dispersions:
DGE.data.2017.19 <- estimateGLMTrendedDisp(DGE.data.2017.19, Design.2017.19)
#To estimate tagwise dispersions:
DGE.data.2017.19 <- estimateGLMTagwiseDisp(DGE.data.2017.19, Design.2017.19)
plotBCV(DGE.data.2017.19)

DGE.data.2019.20 <- estimateGLMCommonDisp(DGE.data.2019.20, Design.2019.20, verbose = TRUE)
# Disp = 0.2233, BCV = 0.4725
#To estimate trended dispersions:
DGE.data.2019.20 <- estimateGLMTrendedDisp(DGE.data.2019.20, Design.2019.20)
#To estimate tagwise dispersions:
DGE.data.2019.20 <- estimateGLMTagwiseDisp(DGE.data.2019.20, Design.2019.20)
plotBCV(DGE.data.2019.20)

DGE.data.2017.18.site <- estimateGLMCommonDisp(DGE.data.2017.18.site, Design.2017.18.site, verbose = TRUE)
# Disp = 0.83162, BCV = 0.9119
#To estimate trended dispersions:
DGE.data.2017.18.site <- estimateGLMTrendedDisp(DGE.data.2017.18.site, Design.2017.18.site)
#To estimate tagwise dispersions:
DGE.data.2017.18.site <- estimateGLMTagwiseDisp(DGE.data.2017.18.site, Design.2017.18.site)
plotBCV(DGE.data.2017.18.site)

### Plot log transformed TMM normalized count data ####
#Print box plots of the log transformed TMM normalized count data in the DGElist object termed DGE.data
#We can use this to exclude sample outliers with abnormally low expression levels

jpeg(
  "./Plots/Time_2017_18_Site_Normalized_Gene_Counts.jpeg",
  width = 1920,
  height = 1080
)
Norml.Count.Data <- cpm(DGE.data.2017.18.site, log = TRUE)

boxplot(Norml.Count.Data, main = "Count Data after transformation", ylab = "log2(cpm)")
dev.off()

### Fit the model ####
# fit a negative binomial generalized log-linear model into our normalized count data 
# using the full design matrix for gene wise statistical tests.

fit.2017.18 <- glmFit(DGE.data.2017.18,Design.2017.18)
fit.2017.18.lrt <- glmLRT(fit.2017.18,coef = "Year2018")

fit.2017.18.site <- glmFit(DGE.data.2017.18.site,Design.2017.18.site)
fit.2017.18.site.lrt <- glmLRT(fit.2017.18.site,coef = "Year2018:SiteSERC")


#### Get top expressed genes ####
# logFC is the log2 fold change between years.
# logFC of 2 would indicate that the gene is expressed 4 times higher in 2018 than 2017. 
# logCPM is the average expression across all samples
# LR is the likelihood ratio L(Full Model)/L(small model)
# PValue is unadjusted p-value
# FDR is the false discovery rate (p-value adjusted for multiple testing) # use this one!

topTags(fit.2017.18.lrt)
topTags(fit.2017.18.site.lrt)

#### Summmary of DGEs ####
#This uses the FDR.  0.05 would be OK also.

# number of down and up regulated genes in 2018 compared to 2017
summary(decideTestsDGE(fit.2017.18.lrt,p.value=0.01)) 
# p 0.01 = 1945 Down, 4001 Up 15696 NotSig
summary(decideTestsDGE(fit.2017.18.lrt,p.value=0.05)) 
# p 0.05 = 2862 Down, 5189 Up 13591 NotSig

summary(decideTestsDGE(fit.2017.18.site.lrt,p.value=0.01)) 
# p 0.01 = 137 Down, 796 Up 20709 NotSig
summary(decideTestsDGE(fit.2017.18.site.lrt,p.value=0.05)) 
# p 0.05 = 534 Down, 1701 Up 19407 NotSig

#Extract genes with a FDR < 0.01 (could also use 0.05)
DEgene.2017.18 <- topTags(fit.2017.18.lrt,n = Inf,p.value = 0.01)$table
DEgene.2017.18.site <- topTags(fit.2017.18.site.lrt,n = Inf,p.value = 0.01)$table


#save to a file
#write.csv(DEgene.2017.18,"./Data/DE.data/DEgenes.2017.2018.csv")


plotDE(rownames(DEgene.2017.18)[1:9],DGE.data.2017.18,Sample_Description.2017.2018)

