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
time.2017.18 = time.2017.18[,c(1:7,10,12,15,18,20:22,24:34,36:40,42:45,47:53,56:68,70:84,86:101,104:107,
                              108:109,111:113,115:129,131:133,135:148,150,152,154,156:168,170:172,
                              174:188,190:192,194:198,200:210,212:223)] # 190 samples, both pops
time.2017.18.HF = time.2017.18[c(1:89)] # 88 samples
time.2017.18.SERC = time.2017.18[c(1,90:191)] # 102 samples

time.2018.19 = Gene_Counts[,c(1,60:107,167:223,224:276)] # 158 samples
time.2018.19 = time.2018.19[,c(1,50,51,54:57,59:62,64:67,69:71,73:77,79:80,82:95,97:108,110:113,
                               115:125,127:131,133:144,146:159)] # 96 samples, only SERC

time.2019.20 = Gene_Counts[,c(1,224:276,277:296)] # 73 samples
time.2019.20 = time.2019.20[,c(1,2,4:8,10:14,17,55,57,59:62,64:67,69,72)] # 24 samples, only SERC, only Spring

time.2017.19 = Gene_Counts[,c(1:59,108:166,224:276)] # 170 samples
time.2017.19 = time.2017.19[,c(1,60:62,64,65,67,69:72,74:77,79:85,87,89:90,92:102,104,
                               106,109:123,125,127:141,143,145:146,148:159,161:171)] # 94 samples, only SERC

time.2017.20 = Gene_Counts[,c(1,108:127,277:296)] # 40 samples
time.2017.20 = time.2017.20[,c(1:2,4:7,9,11:14,16,19:20,22,24:27,29,31:34,36,39:40)] # 26 samples, only SERC

time.2018.20 = Gene_Counts[,c(1,167:186,277:296)] # 40 samples
time.2018.20 = time.2018.20[,c(1:2,5:9,11:14,16,19:20,22,25:29,31:34,36,39:40)] # 26 samples, only SERC

# divide into seasons: do differences in spring expression between years explain growth differences between years

time.2017.18.Spring = time.2017.18[c(1:13,46:57,90:107,141:158)] # 60 samples
time.2017.18.Summer= time.2017.18[c(1,14:29,58:73,108:124,159:175)] # 66 samples
time.2017.18.Fall = time.2017.18[c(1,30:45,74:89,125:140,176:191)] # 64 samples

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
Sample_Description.2017.2018.HF = Sample_Description %>%
  filter(sample.description %in% colnames(time.2017.18.HF)) %>%
  droplevels(.)
Sample_Description.2017.2018.SERC = Sample_Description %>%
  filter(sample.description %in% colnames(time.2017.18.SERC)) %>%
  droplevels(.)
Sample_Description.2017.2020 = Sample_Description %>%
  filter(sample.description %in% colnames(time.2017.20)) %>%
  droplevels(.)
Sample_Description.2018.2020 = Sample_Description %>%
  filter(sample.description %in% colnames(time.2018.20)) %>%
  droplevels(.)
Sample_Description.2017.2018.Spring = Sample_Description %>%
  filter(sample.description %in% colnames(time.2017.18.Spring)) %>%
  droplevels(.)
Sample_Description.2017.2018.Summer = Sample_Description %>%
  filter(sample.description %in% colnames(time.2017.18.Summer)) %>%
  droplevels(.)
Sample_Description.2017.2018.Fall = Sample_Description %>%
  filter(sample.description %in% colnames(time.2017.18.Fall)) %>%
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

Data.matrix.2017.18.HF <- time.2017.18.HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.2017.18.HF) <- time.2017.18.HF$Gene_ID

Data.matrix.2017.18.SERC <- time.2017.18.SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.2017.18.SERC) <- time.2017.18.SERC$Gene_ID

Data.matrix.2017.20 <- time.2017.20 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.2017.20) <- time.2017.20$Gene_ID

Data.matrix.2018.20 <- time.2018.20 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.2018.20) <- time.2018.20$Gene_ID

Data.matrix.2017.18.Spring <- time.2017.18.Spring %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.2017.18.Spring) <- time.2017.18.Spring$Gene_ID

Data.matrix.2017.18.Summer <- time.2017.18.Summer %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.2017.18.Summer) <- time.2017.18.Summer$Gene_ID

Data.matrix.2017.18.Fall <- time.2017.18.Fall %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.2017.18.Fall) <- time.2017.18.Fall$Gene_ID


#The DGEList function creates a DGElist object for differential expression analysis. 
# we specify counts to be the data matrix we created above that holds our count data. 
# Group is specified from the Year in Sample Description we created above which contains information about the species and treatment of each sample.
DGE.data.2017.18 = DGEList(counts = Data.matrix.2017.18, group = Sample_Description.2017.2018$Year)
DGE.data.2018.19 = DGEList(counts = Data.matrix.2018.19, group = Sample_Description.2018.2019$Year)
DGE.data.2017.19 = DGEList(counts = Data.matrix.2017.19, group = Sample_Description.2017.2019$Year)
DGE.data.2019.20 = DGEList(counts = Data.matrix.2019.20, group = Sample_Description.2019.2020$Year)
DGE.data.2017.18.HF = DGEList(counts = Data.matrix.2017.18.HF, group = Sample_Description.2017.2018.HF$Year)
DGE.data.2017.18.SERC = DGEList(counts = Data.matrix.2017.18.SERC, group = Sample_Description.2017.2018.SERC$Year)
DGE.data.2017.20 = DGEList(counts = Data.matrix.2017.20, group = Sample_Description.2017.2020$Year)
DGE.data.2018.20 = DGEList(counts = Data.matrix.2018.20, group = Sample_Description.2018.2020$Year)
DGE.data.2017.18.Spring = DGEList(counts = Data.matrix.2017.18.Spring, group = Sample_Description.2017.2018.Spring$Year)
DGE.data.2017.18.Summer = DGEList(counts = Data.matrix.2017.18.Summer, group = Sample_Description.2017.2018.Summer$Year)
DGE.data.2017.18.Fall = DGEList(counts = Data.matrix.2017.18.Fall, group = Sample_Description.2017.2018.Fall$Year)

# for 2017-2018 use group as site and year
DGE.data.2017.18.site = DGEList(counts = Data.matrix.2017.18, group = Sample_Description.2017.2018$group)

#### Normalize the data #####
# using the TMM method, trimmed mean of M-values
DGE.data.2017.18 <- calcNormFactors(DGE.data.2017.18, method = "TMM")
DGE.data.2018.19 <- calcNormFactors(DGE.data.2018.19, method = "TMM")
DGE.data.2017.19 <- calcNormFactors(DGE.data.2017.19, method = "TMM")
DGE.data.2019.20 <- calcNormFactors(DGE.data.2019.20, method = "TMM")
DGE.data.2017.18.HF <- calcNormFactors(DGE.data.2017.18.HF, method = "TMM")
DGE.data.2017.18.SERC <- calcNormFactors(DGE.data.2017.18.SERC, method = "TMM")
DGE.data.2017.20 <- calcNormFactors(DGE.data.2017.20, method = "TMM")
DGE.data.2018.20 <- calcNormFactors(DGE.data.2018.20, method = "TMM")
DGE.data.2017.18.Spring <- calcNormFactors(DGE.data.2017.18.Spring, method = "TMM")
DGE.data.2017.18.Summer <- calcNormFactors(DGE.data.2017.18.Summer, method = "TMM")
DGE.data.2017.18.Fall <- calcNormFactors(DGE.data.2017.18.Fall, method = "TMM")

DGE.data.2017.18.site <- calcNormFactors(DGE.data.2017.18.site, method = "TMM")

#### Plotting ####
# bcv is the square root of the dispersion of the negative binomial distribution
plotMDS(DGE.data.2017.18, method = "bcv") 
plotMDS(DGE.data.2018.19, method = "bcv") 
plotMDS(DGE.data.2017.19, method = "bcv") 
plotMDS(DGE.data.2019.20, method = "bcv") 
plotMDS(DGE.data.2017.18.HF, method = "bcv") 
plotMDS(DGE.data.2017.18.SERC, method = "bcv") 
plotMDS(DGE.data.2017.20, method = "bcv") 
plotMDS(DGE.data.2018.20, method = "bcv") 
plotMDS(DGE.data.2017.18.Spring, method = "bcv")
plotMDS(DGE.data.2017.18.Summer, method = "bcv") 
plotMDS(DGE.data.2017.18.Fall, method = "bcv") 
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
TMM_DATA_2017.18.HF <- cpm(DGE.data.2017.18.HF, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2017.18.HF,"./Data/DE.data/TMM_NormData_LogCPM_2017_2018.HF.csv")
TMM_DATA_2017.18.SERC <- cpm(DGE.data.2017.18.SERC, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2017.18.SERC,"./Data/DE.data/TMM_NormData_LogCPM_2017_2018.SERC.csv")
TMM_DATA_2017.20 <- cpm(DGE.data.2017.20, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2017.20,"./Data/DE.data/TMM_NormData_LogCPM_2017_2020.csv")
TMM_DATA_2018.20 <- cpm(DGE.data.2018.20, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2018.20,"./Data/DE.data/TMM_NormData_LogCPM_2018_2020.csv")
TMM_DATA_2017.18.Spring <- cpm(DGE.data.2017.18.Spring, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2017.18.Spring,"./Data/DE.data/TMM_NormData_LogCPM_2017_2018.Spring.csv")
TMM_DATA_2017.18.Summer <- cpm(DGE.data.2017.18.Summer, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2017.18.Summer,"./Data/DE.data/TMM_NormData_LogCPM_2017_2018.Summer.csv")
TMM_DATA_2017.18.Fall <- cpm(DGE.data.2017.18.Fall, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2017.18.Fall,"./Data/DE.data/TMM_NormData_LogCPM_2017_2018.Fall.csv")
TMM_DATA_2017.18.site <- cpm(DGE.data.2017.18.site, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM_DATA_2017.18.site,"./Data/DE.data/TMM_NormData_LogCPM_2017_2018.site.csv")

#Create our design matrix which considers every interaction between species and treatment, in this case Year.
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

Design.2017.18.HF <-
  model.matrix(~ Year, data = Sample_Description.2017.2018.HF) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2017.18.HF) <- Sample_Description.2017.2018.HF$sample.description

Design.2017.18.SERC <-
  model.matrix(~ Year, data = Sample_Description.2017.2018.SERC) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2017.18.SERC) <- Sample_Description.2017.2018.SERC$sample.description

Design.2017.20 <-
  model.matrix(~ Year, data = Sample_Description.2017.2020) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2017.20) <- Sample_Description.2017.2020$sample.description

Design.2018.20 <-
  model.matrix(~ Year, data = Sample_Description.2018.2020) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2018.20) <- Sample_Description.2018.2020$sample.description

Design.2017.18.Spring <-
  model.matrix(~ Year, data = Sample_Description.2017.2018.Spring) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2017.18.Spring) <- Sample_Description.2017.2018.Spring$sample.description

Design.2017.18.Summer <-
  model.matrix(~ Year, data = Sample_Description.2017.2018.Summer) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2017.18.Summer) <- Sample_Description.2017.2018.Summer$sample.description

Design.2017.18.Fall <-
  model.matrix(~ Year, data = Sample_Description.2017.2018.Fall) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2017.18.Fall) <- Sample_Description.2017.2018.Fall$sample.description

# testing additive and multiplicative year/site models
Design.2017.18.site.additive <-
  model.matrix(~ Year + Site, data = Sample_Description.2017.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2017.18.site.additive) <- Sample_Description.2017.2018$sample.description

Design.2017.18.site <-
  model.matrix(~ Year*Site, data = Sample_Description.2017.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.2017.18.site) <- Sample_Description.2017.2018$sample.description

#### Calculate Dispersion Factors ####
# To estimate common dispersion:
DGE.data.2017.18 <- estimateGLMCommonDisp(DGE.data.2017.18, Design.2017.18, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2017.18 <- estimateGLMTrendedDisp(DGE.data.2017.18, Design.2017.18)
#To estimate tagwise dispersions:
DGE.data.2017.18 <- estimateGLMTagwiseDisp(DGE.data.2017.18, Design.2017.18)
plotBCV(DGE.data.2017.18)

# To estimate common dispersion:
DGE.data.2018.19 <- estimateGLMCommonDisp(DGE.data.2018.19, Design.2018.19, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2018.19 <- estimateGLMTrendedDisp(DGE.data.2018.19, Design.2018.19)
#To estimate tagwise dispersions:
DGE.data.2018.19 <- estimateGLMTagwiseDisp(DGE.data.2018.19, Design.2018.19)
plotBCV(DGE.data.2018.19)

# To estimate common dispersion:
DGE.data.2017.19 <- estimateGLMCommonDisp(DGE.data.2017.19, Design.2017.19, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2017.19 <- estimateGLMTrendedDisp(DGE.data.2017.19, Design.2017.19)
#To estimate tagwise dispersions:
DGE.data.2017.19 <- estimateGLMTagwiseDisp(DGE.data.2017.19, Design.2017.19)
plotBCV(DGE.data.2017.19)

# To estimate common dispersion:
DGE.data.2019.20 <- estimateGLMCommonDisp(DGE.data.2019.20, Design.2019.20, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2019.20 <- estimateGLMTrendedDisp(DGE.data.2019.20, Design.2019.20)
#To estimate tagwise dispersions:
DGE.data.2019.20 <- estimateGLMTagwiseDisp(DGE.data.2019.20, Design.2019.20)
plotBCV(DGE.data.2019.20)

# To estimate common dispersion:
DGE.data.2017.18.HF <- estimateGLMCommonDisp(DGE.data.2017.18.HF, Design.2017.18.HF, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2017.18.HF <- estimateGLMTrendedDisp(DGE.data.2017.18.HF, Design.2017.18.HF)
#To estimate tagwise dispersions:
DGE.data.2017.18.HF <- estimateGLMTagwiseDisp(DGE.data.2017.18.HF, Design.2017.18.HF)
plotBCV(DGE.data.2017.18.HF)

# To estimate common dispersion:
DGE.data.2017.18.SERC <- estimateGLMCommonDisp(DGE.data.2017.18.SERC, Design.2017.18.SERC, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2017.18.SERC <- estimateGLMTrendedDisp(DGE.data.2017.18.SERC, Design.2017.18.SERC)
#To estimate tagwise dispersions:
DGE.data.2017.18.SERC <- estimateGLMTagwiseDisp(DGE.data.2017.18.SERC, Design.2017.18.SERC)
plotBCV(DGE.data.2017.18.SERC)

# To estimate common dispersion:
DGE.data.2017.20 <- estimateGLMCommonDisp(DGE.data.2017.20, Design.2017.20, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2017.20 <- estimateGLMTrendedDisp(DGE.data.2017.20, Design.2017.20)
#To estimate tagwise dispersions:
DGE.data.2017.20 <- estimateGLMTagwiseDisp(DGE.data.2017.20, Design.2017.20)
plotBCV(DGE.data.2017.20)

# To estimate common dispersion:
DGE.data.2018.20 <- estimateGLMCommonDisp(DGE.data.2018.20, Design.2018.20, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2018.20 <- estimateGLMTrendedDisp(DGE.data.2018.20, Design.2018.20)
#To estimate tagwise dispersions:
DGE.data.2018.20 <- estimateGLMTagwiseDisp(DGE.data.2018.20, Design.2018.20)
plotBCV(DGE.data.2018.20)

# To estimate common dispersion:
DGE.data.2017.18.Spring <- estimateGLMCommonDisp(DGE.data.2017.18.Spring, Design.2017.18.Spring, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2017.18.Spring <- estimateGLMTrendedDisp(DGE.data.2017.18.Spring, Design.2017.18.Spring)
#To estimate tagwise dispersions:
DGE.data.2017.18.Spring <- estimateGLMTagwiseDisp(DGE.data.2017.18.Spring, Design.2017.18.Spring)
plotBCV(DGE.data.2017.18.Spring)

# To estimate common dispersion:
DGE.data.2017.18.Summer <- estimateGLMCommonDisp(DGE.data.2017.18.Summer, Design.2017.18.Summer, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2017.18.Summer <- estimateGLMTrendedDisp(DGE.data.2017.18.Summer, Design.2017.18.Summer)
#To estimate tagwise dispersions:
DGE.data.2017.18.Summer <- estimateGLMTagwiseDisp(DGE.data.2017.18.Summer, Design.2017.18.Summer)
plotBCV(DGE.data.2017.18.Summer)

# To estimate common dispersion:
DGE.data.2017.18.Fall <- estimateGLMCommonDisp(DGE.data.2017.18.Fall, Design.2017.18.Fall, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2017.18.Fall <- estimateGLMTrendedDisp(DGE.data.2017.18.Fall, Design.2017.18.Fall)
#To estimate tagwise dispersions:
DGE.data.2017.18.Fall <- estimateGLMTagwiseDisp(DGE.data.2017.18.Fall, Design.2017.18.Fall)
plotBCV(DGE.data.2017.18.Fall)

# additive design
DGE.data.2017.18.site <- estimateGLMCommonDisp(DGE.data.2017.18.site, Design.2017.18.site.additive, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2017.18.site <- estimateGLMTrendedDisp(DGE.data.2017.18.site, Design.2017.18.site.additive)
#To estimate tagwise dispersions:
DGE.data.2017.18.site <- estimateGLMTagwiseDisp(DGE.data.2017.18.site, Design.2017.18.site.additive)
plotBCV(DGE.data.2017.18.site)

# multiplicative design
DGE.data.2017.18.site <- estimateGLMCommonDisp(DGE.data.2017.18.site, Design.2017.18.site, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.2017.18.site <- estimateGLMTrendedDisp(DGE.data.2017.18.site, Design.2017.18.site)
#To estimate tagwise dispersions:
DGE.data.2017.18.site <- estimateGLMTagwiseDisp(DGE.data.2017.18.site, Design.2017.18.site)
plotBCV(DGE.data.2017.18.site)


### Plot log transformed TMM normalized count data ####
#Print box plots of the log transformed TMM normalized count data in the DGElist object termed DGE.data
#We can use this to exclude sample outliers with abnormally low expression levels
# repeat code for each sample grouping
jpeg(
  "./Plots/DE/Time_2017_18_Spring_Normalized_Gene_Counts.jpeg",
  width = 1920,
  height = 1080
)
Norml.Count.Data <- cpm(DGE.data.2017.18.Spring, log = TRUE)

boxplot(Norml.Count.Data, main = "Count Data after transformation", ylab = "log2(cpm)")
dev.off()

### Fit the model ####
# fit a negative binomial generalized log-linear model to our normalized count data 
# using the full design matrix for gene wise statistical tests.

fit.2017.18 <- glmFit(DGE.data.2017.18,Design.2017.18)
fit.2017.18.lrt <- glmLRT(fit.2017.18,coef = "Year2018")

fit.2018.19 <- glmFit(DGE.data.2018.19,Design.2018.19)
fit.2018.19.lrt <- glmLRT(fit.2018.19,coef = "Year2019")

fit.2017.19 <- glmFit(DGE.data.2017.19,Design.2017.19)
fit.2017.19.lrt <- glmLRT(fit.2017.19,coef = "Year2019")

fit.2019.20 <- glmFit(DGE.data.2019.20,Design.2019.20)
fit.2019.20.lrt <- glmLRT(fit.2019.20,coef = "Year2020")

fit.2017.18.HF <- glmFit(DGE.data.2017.18.HF,Design.2017.18.HF)
fit.2017.18.HF.lrt <- glmLRT(fit.2017.18.HF,coef = "Year2018")

fit.2017.18.SERC <- glmFit(DGE.data.2017.18.SERC,Design.2017.18.SERC)
fit.2017.18.SERC.lrt <- glmLRT(fit.2017.18.SERC,coef = "Year2018")

fit.2017.20 <- glmFit(DGE.data.2017.20,Design.2017.20)
fit.2017.20.lrt <- glmLRT(fit.2017.20,coef = "Year2020")

fit.2018.20 <- glmFit(DGE.data.2018.20,Design.2018.20)
fit.2018.20.lrt <- glmLRT(fit.2018.20,coef = "Year2020")

fit.2017.18.Spring <- glmFit(DGE.data.2017.18.Spring,Design.2017.18.Spring)
fit.2017.18.Spring.lrt <- glmLRT(fit.2017.18.Spring,coef = "Year2018")

fit.2017.18.Summer <- glmFit(DGE.data.2017.18.Summer,Design.2017.18.Summer)
fit.2017.18.Summer.lrt <- glmLRT(fit.2017.18.Summer,coef = "Year2018")

fit.2017.18.Fall <- glmFit(DGE.data.2017.18.Fall,Design.2017.18.Fall)
fit.2017.18.Fall.lrt <- glmLRT(fit.2017.18.Fall,coef = "Year2018")

# additive model
fit.2017.18.site.additive <- glmFit(DGE.data.2017.18.site,Design.2017.18.site.additive)
fit.2017.18.additive.YEAR.lrt <- glmLRT(fit.2017.18.site.additive,coef = "Year2018")
fit.2017.18.additive.SITE.lrt <- glmLRT(fit.2017.18.site.additive,coef = "SiteSERC")

# multiplicative model
fit.2017.18.site <- glmFit(DGE.data.2017.18.site,Design.2017.18.site)
fit.2017.18.site.multiplicative.lrt <- glmLRT(fit.2017.18.site,coef = "Year2018:SiteSERC")

#### Get top expressed genes ####
# logFC is the log2 fold change between years.
# logFC of 2 would indicate that the gene is expressed 4 times higher in 2018 than 2017. 
# logCPM is the average expression across all samples
# LR is the likelihood ratio L(Full Model)/L(small model)
# PValue is unadjusted p-value
# FDR is the false discovery rate (p-value adjusted for multiple testing) # use this one!

topTags(fit.2017.18.lrt)
topTags(fit.2017.18.site.lrt)
topTags(fit.2017.18.site.YEAR.lrt)

#### Summmary of DGEs ####
#This uses the FDR of 0.01

# number of down and up regulated genes in 2018 compared to 2017
summary(decideTestsDGE(fit.2017.18.lrt,p.value=0.01)) 
# p 0.01 = 1991 Down, 4090 Up 15561 NotSig
summary(decideTestsDGE(fit.2018.19.lrt,p.value=0.01)) 
# p 0.01 = 4905 Down, 2545 Up 14192 NotSig
summary(decideTestsDGE(fit.2019.20.lrt,p.value=0.01)) 
# p 0.01 = 2558 Down, 2995 Up 15979 NotSig
summary(decideTestsDGE(fit.2017.18.HF.lrt,p.value=0.01)) 
# p 0.01 = 977 Down, 1064 Up 19601 NotSig
summary(decideTestsDGE(fit.2017.18.SERC.lrt,p.value=0.01)) 
# p 0.01 = 1567 Down, 3715 Up 16360 NotSig
summary(decideTestsDGE(fit.2017.20.lrt,p.value=0.01)) 
# p 0.01 = 1567 Down, 1805 Up 18270 NotSig
summary(decideTestsDGE(fit.2018.20.lrt,p.value=0.01)) 
# p 0.01 = 4334 Down, 3968 Up 13340 NotSig
summary(decideTestsDGE(fit.2017.18.Spring.lrt,p.value=0.01)) 
# p 0.01 = 3693 Down, 4488 Up 13461 NotSig
summary(decideTestsDGE(fit.2017.18.Summer.lrt,p.value=0.01)) 
# p 0.01 = 665 Down, 369 Up 20608 NotSig
summary(decideTestsDGE(fit.2017.18.Fall.lrt,p.value=0.01)) 
# p 0.01 = 1056 Down, 1369 Up 19217 NotSig

# additive model
summary(decideTestsDGE(fit.2017.18.additive.YEAR.lrt,p.value=0.01)) 
# p 0.01 = 2004 Down, 3979 Up 15659 NotSig
summary(decideTestsDGE(fit.2017.18.additive.SITE.lrt,p.value=0.01)) 
# p 0.01 = 1989 Down, 1730 Up 17923 NotSig

# multiplicative
summary(decideTestsDGE(fit.2017.18.site.multiplicative.lrt,p.value=0.01)) 
# p 0.01 = 75 Down, 431 Up 21136 NotSig


#### Extract genes with a FDR < 0.01 ####
DEgene.2017.18 <- topTags(fit.2017.18.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2017.18,"./Data/DE.data/DEgenes.2017.2018.csv")
DEgene.2018.19 <- topTags(fit.2018.19.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2018.19,"./Data/DE.data/DEgenes.2018.2019.csv")
DEgene.2017.19 <- topTags(fit.2017.19.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2017.19,"./Data/DE.data/DEgenes.2017.2019.csv")
DEgene.2019.20 <- topTags(fit.2019.20.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2019.20,"./Data/DE.data/DEgenes.2019.2020.csv")
DEgene.2017.18.HF <- topTags(fit.2017.18.HF.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2017.18.HF,"./Data/DE.data/DEgenes.2017.2018.HF.csv")
DEgene.2017.18.SERC <- topTags(fit.2017.18.SERC.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2017.18.SERC,"./Data/DE.data/DEgenes.2017.2018.SERC.csv")
DEgene.2017.20 <- topTags(fit.2017.20.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2017.20,"./Data/DE.data/DEgenes.2017.2020.csv")
DEgene.2018.20 <- topTags(fit.2018.20.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2018.20,"./Data/DE.data/DEgenes.2018.2020.csv")
DEgene.2017.18.Spring <- topTags(fit.2017.18.Spring.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2017.18.Spring,"./Data/DE.data/DEgenes.2017.2018.Spring.csv")
DEgene.2017.18.Summer <- topTags(fit.2017.18.Summer.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2017.18.Summer,"./Data/DE.data/DEgenes.2017.2018.Summer.csv")
DEgene.2017.18.Fall <- topTags(fit.2017.18.Fall.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2017.18.Fall,"./Data/DE.data/DEgenes.2017.2018.Fall.csv")
# additive model
DEgene.2017.18.additive.Year <- topTags(fit.2017.18.additive.YEAR.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2017.18.additive.Year,"./Data/DE.data/DEgenes.2017.2018.additive.Year.csv")
DEgene.2017.18.additive.Site <- topTags(fit.2017.18.additive.SITE.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2017.18.additive.Site,"./Data/DE.data/DEgenes.2017.2018.additive.Site.csv")
# multiplicative
DEgene.2017.18.multiply.site.year <- topTags(fit.2017.18.site.multiplicative.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.2017.18.multiply.site.year,"./Data/DE.data/DEgenes.2017.2018.multiply.Site.Year.csv")

#### Plotting top nine more differentially expressed genes in each set #####
plotDE <- function(genes, dge, sample.description) {
  require(ggplot2)
  tmp.data <- t(log2(cpm(dge[genes,])+1))
  tmp.data <- tmp.data %>%
    as.data.frame() %>%
    rownames_to_column("sample.description") %>%
    left_join(sample.description,by="sample.description")
  tmp.data <- tmp.data %>%
    pivot_longer(cols=starts_with("FAGR"), values_to = "log2_cpm", names_to = "gene")
  pl <- ggplot(tmp.data,aes(x=Year,y=log2_cpm))
  pl <- pl + facet_wrap( ~ gene)
  pl <- pl + ylab("log2(cpm)") + xlab("genotype")
  pl <- pl + geom_boxplot()
  pl + theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1))
}
plotDE_fill <- function(genes, dge, sample.description) {
  require(ggplot2)
  tmp.data <- t(log2(cpm(dge[genes,])+1))
  tmp.data <- tmp.data %>%
    as.data.frame() %>%
    rownames_to_column("sample.description") %>%
    left_join(sample.description,by="sample.description")
  tmp.data <- tmp.data %>%
    pivot_longer(cols=starts_with("FAGR"), values_to = "log2_cpm", names_to = "gene")
  pl <- ggplot(tmp.data,aes(x=Year,y=log2_cpm,fill=Site))
  pl <- pl + facet_wrap( ~ gene)
  pl <- pl + ylab("log2(cpm)") + xlab("genotype")
  pl <- pl + geom_boxplot()
  pl + theme(axis.text.x  = element_text(angle=45, vjust=1,hjust=1))
}

plotDE(rownames(DEgene.2017.18)[1:9],DGE.data.2017.18,Sample_Description.2017.2018)
plotDE(rownames(DEgene.2018.19)[1:9],DGE.data.2018.19,Sample_Description.2018.2019)
plotDE(rownames(DEgene.2017.19)[1:9],DGE.data.2017.19,Sample_Description.2017.2019)
plotDE(rownames(DEgene.2019.20)[1:9],DGE.data.2019.20,Sample_Description.2019.2020)
plotDE(rownames(DEgene.2017.18.HF)[1:9],DGE.data.2017.18.HF,Sample_Description.2017.2018.HF)
plotDE(rownames(DEgene.2017.18.SERC)[1:9],DGE.data.2017.18.SERC,Sample_Description.2017.2018.SERC)
plotDE(rownames(DEgene.2017.20)[1:9],DGE.data.2017.20,Sample_Description.2017.2020)
plotDE(rownames(DEgene.2018.20)[1:9],DGE.data.2018.20,Sample_Description.2018.2020)
plotDE(rownames(DEgene.2017.18.Spring)[1:9],DGE.data.2017.18.Spring,Sample_Description.2017.2018.Spring)
plotDE(rownames(DEgene.2017.18.Summer)[1:9],DGE.data.2017.18.Summer,Sample_Description.2017.2018.Summer)
plotDE(rownames(DEgene.2017.18.Fall)[1:9],DGE.data.2017.18.Fall,Sample_Description.2017.2018.Fall)
# additive models
plotDE_fill(rownames(DEgene.2017.18.additive.Year)[1:9],DGE.data.2017.18.site,Sample_Description.2017.2018)
# change Year to site for X axis in plotting function
plotDE_fill(rownames(DEgene.2017.18.additive.Site)[1:9],DGE.data.2017.18.site,Sample_Description.2017.2018)
# multiplicative
plotDE_fill(rownames(DEgene.2017.18.multiply.site.year)[1:9],DGE.data.2017.18.site,Sample_Description.2017.2018)

#### Subset TMM by DEGs for use in WGCNA ####

TMM_DATA_2017.18_DEG = subset(TMM_DATA_2017.18, TMM_DATA_2017.18$Gene_ID %in% rownames(DEgene.2017.18))
#write.csv(TMM_DATA_2017.18_DEG,"./Data/DE.data/TMM_2017.18_DEG.csv")
TMM_DATA_2018.19_DEG = subset(TMM_DATA_2018.19, TMM_DATA_2018.19$Gene_ID %in% rownames(DEgene.2018.19))
#write.csv(TMM_DATA_2018.19_DEG,"./Data/DE.data/TMM_2018.19_DEG.csv")
TMM_DATA_2017.19_DEG = subset(TMM_DATA_2017.19, TMM_DATA_2017.19$Gene_ID %in% rownames(DEgene.2017.19))
#write.csv(TMM_DATA_2017.19_DEG,"./Data/DE.data/TMM_2017.19_DEG.csv")
TMM_DATA_2019.20_DEG = subset(TMM_DATA_2019.20, TMM_DATA_2019.20$Gene_ID %in% rownames(DEgene.2019.20))
#write.csv(TMM_DATA_2019.20_DEG,"./Data/DE.data/TMM_2019.20_DEG.csv")
TMM_DATA_2017.18.HF_DEG = subset(TMM_DATA_2017.18.HF, TMM_DATA_2017.18.HF$Gene_ID %in% rownames(DEgene.2017.18.HF))
#write.csv(TMM_DATA_2017.18.HF_DEG,"./Data/DE.data/TMM_2017.18.HF_DEG.csv")
TMM_DATA_2017.18.SERC_DEG = subset(TMM_DATA_2017.18.SERC, TMM_DATA_2017.18.SERC$Gene_ID %in% rownames(DEgene.2017.18.SERC))
#write.csv(TMM_DATA_2017.18.SERC_DEG,"./Data/DE.data/TMM_2017.18.SERC_DEG.csv")
TMM_DATA_2017.20_DEG = subset(TMM_DATA_2017.20, TMM_DATA_2017.20$Gene_ID %in% rownames(DEgene.2017.20))
#write.csv(TMM_DATA_2017.20_DEG,"./Data/DE.data/TMM_2017.20_DEG.csv")
TMM_DATA_2018.20_DEG = subset(TMM_DATA_2018.20, TMM_DATA_2018.20$Gene_ID %in% rownames(DEgene.2018.20))
#write.csv(TMM_DATA_2018.20_DEG,"./Data/DE.data/TMM_2018.20_DEG.csv")
TMM_DATA_2017.18.Spring_DEG = subset(TMM_DATA_2017.18.Spring, TMM_DATA_2017.18.Spring$Gene_ID %in% rownames(DEgene.2017.18.Spring))
#write.csv(TMM_DATA_2017.18.Spring_DEG,"./Data/DE.data/TMM_2017.18.Spring_DEG.csv")
TMM_DATA_2017.18.Summer_DEG = subset(TMM_DATA_2017.18.Summer, TMM_DATA_2017.18.Summer$Gene_ID %in% rownames(DEgene.2017.18.Summer))
#write.csv(TMM_DATA_2017.18.Summer_DEG,"./Data/DE.data/TMM_2017.18.Summer_DEG.csv")
TMM_DATA_2017.18.Fall_DEG = subset(TMM_DATA_2017.18.Fall, TMM_DATA_2017.18.Fall$Gene_ID %in% rownames(DEgene.2017.18.Fall))
#write.csv(TMM_DATA_2017.18.Fall_DEG,"./Data/DE.data/TMM_2017.18.Fall_DEG.csv")
# additive
TMM_DATA_2017.18.additive.Site_DEG = subset(TMM_DATA_2017.18.site, TMM_DATA_2017.18.site$Gene_ID %in% rownames(DEgene.2017.18.additive.Site))
#write.csv(TMM_DATA_2017.18.additive.Site_DEG,"./Data/DE.data/TMM_2017.18.additive.Site_DEG.csv")
TMM_DATA_2017.18.additive.Year_DEG = subset(TMM_DATA_2017.18.site, TMM_DATA_2017.18.site$Gene_ID %in% rownames(DEgene.2017.18.additive.Year))
#write.csv(TMM_DATA_2017.18.additive.Year_DEG,"./Data/DE.data/TMM_2017.18.additive.Year_DEG.csv")
#multiplicative
TMM_DATA_2017.18.multiply.Site.Year_DEG = subset(TMM_DATA_2017.18.site, TMM_DATA_2017.18.site$Gene_ID %in% rownames(DEgene.2017.18.multiply.site.year))
#write.csv(TMM_DATA_2017.18.multiply.Site.Year_DEG,"./Data/DE.data/TMM_2017.18.multiply.Site.Year_DEG.csv")

