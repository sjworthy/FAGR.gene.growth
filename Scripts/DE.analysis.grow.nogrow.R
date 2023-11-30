#### Differential Expression Analysis ####

# comparisons to make
# grow versus no grow
# 2017, 2018 at HF
# 2017, 2018, 2019, 2020 at SERC
# need to make sure have matching trees in each set for comparison

#Libraries
library(tidyverse)
library(magrittr)
library(randomcoloR)
library(edgeR)

### Subset parameter table with expression samples ####

params = read.csv("./Formatted.Data/Parameter_table_FAGUS-1.csv")
Sample.Description = read_csv("./Formatted.Data/FAGR.description.csv")

params.2 = subset(params, params$TREE_ID %in% Sample.Description$Tree_ID)
test = params.2 %>%
  select(SITE,YEAR,TREE_ID,GROWTH_SIGNAL)
write.csv(test, file = "param.expression.match.csv")

summary(test)

table(Sample.Description$growth.signal, Sample.Description$Year, Sample.Description$Season)

#### Read in our TPM normalized RNA-seq expression data ####

TPM.ExprData <- read_csv("./Raw.Data/HTseq-master-counts.csv")
# 295 samples

#### Expression Filtering ####
# Remove low expression genes
# Need at least 10 reads in at least 3 samples

Gene.Counts <- TPM.ExprData[rowSums(TPM.ExprData[,-1] > 10) >= 3,]

#### Filtering Samples ####
#read in sample descriptions

# need to make sure we have two time points for each sample with matching season
# all samples from both years with matching season time points
spring.2017 = Gene.Counts[,c(1:21,108:127)]
summer.2017 = Gene.Counts[,c(1,22:39,128:146)]
fall.2017 = Gene.Counts[,c(1,40:59,147:166)]

spring.2018 = Gene.Counts[,c(1,60:72,167:186)]
summer.2018 = Gene.Counts[,c(1,73:89,187:206)]
fall.2018 = Gene.Counts[,c(1,90:107,207:223)]

spring.2019 = Gene.Counts[,c(1,224:240)]
summer.2019 = Gene.Counts[,c(1,241:259)]
fall.2019 = Gene.Counts[,c(1,260:276)]

spring.2020 = Gene.Counts[,c(1,277:296)]

spring.2017.HF = Gene.Counts[,c(1:21)]
summer.2017.HF = Gene.Counts[,c(1,22:39)]
fall.2017.HF = Gene.Counts[,c(1,40:59)]

spring.2017.SERC = Gene.Counts[,c(1,108:127)]
summer.2017.SERC = Gene.Counts[,c(1,128:146)]
fall.2017.SERC = Gene.Counts[,c(1, 147:166)]

spring.2018.HF = Gene.Counts[,c(1,60:72)]
summer.2018.HF = Gene.Counts[,c(1,73:89)]
fall.2018.HF = Gene.Counts[,c(1,90:107)]

spring.2018.SERC = Gene.Counts[,c(1,167:186)]
summer.2018.SERC = Gene.Counts[,c(1,187:206)]
fall.2018.SERC = Gene.Counts[,c(1,207:223)]

#### Sample Description ####

Sample.Description = read_csv("./Formatted.Data/FAGR.description.csv")
# convert grow.nogrow to factor to group samples
Sample.Description$Year = as.factor(Sample.Description$growth.signal)

Sample.Description.spring.2017 = Sample.Description %>%
  filter(sample.description %in% colnames(spring.2017)) %>%
  droplevels(.)
Sample.Description.summer.2017 = Sample.Description %>%
  filter(sample.description %in% colnames(summer.2017)) %>%
  droplevels(.)
Sample.Description.fall.2017 = Sample.Description %>%
  filter(sample.description %in% colnames(fall.2017)) %>%
  droplevels(.)
Sample.Description.spring.2018 = Sample.Description %>%
  filter(sample.description %in% colnames(spring.2018)) %>%
  droplevels(.)
Sample.Description.summer.2018 = Sample.Description %>%
  filter(sample.description %in% colnames(summer.2018)) %>%
  droplevels(.)
Sample.Description.fall.2018 = Sample.Description %>%
  filter(sample.description %in% colnames(fall.2018)) %>%
  droplevels(.)
Sample.Description.spring.2019 = Sample.Description %>%
  filter(sample.description %in% colnames(spring.2019)) %>%
  droplevels(.)
Sample.Description.summer.2019 = Sample.Description %>%
  filter(sample.description %in% colnames(summer.2019)) %>%
  droplevels(.)
Sample.Description.fall.2019 = Sample.Description %>%
  filter(sample.description %in% colnames(fall.2019)) %>%
  droplevels(.)
Sample.Description.spring.2020 = Sample.Description %>%
  filter(sample.description %in% colnames(spring.2020)) %>%
  droplevels(.)

Sample.Description.spring.2017.HF = Sample.Description %>%
  filter(sample.description %in% colnames(spring.2017.HF)) %>%
  droplevels(.)
Sample.Description.summer.2017.HF = Sample.Description %>%
  filter(sample.description %in% colnames(summer.2017.HF)) %>%
  droplevels(.)
Sample.Description.fall.2017.HF = Sample.Description %>%
  filter(sample.description %in% colnames(fall.2017.HF)) %>%
  droplevels(.)
Sample.Description.spring.2017.SERC = Sample.Description %>%
  filter(sample.description %in% colnames(spring.2017.SERC)) %>%
  droplevels(.)
Sample.Description.summer.2017.SERC = Sample.Description %>%
  filter(sample.description %in% colnames(summer.2017.SERC)) %>%
  droplevels(.)
Sample.Description.fall.2017.SERC = Sample.Description %>%
  filter(sample.description %in% colnames(fall.2017.SERC)) %>%
  droplevels(.)

Sample.Description.spring.2018.HF = Sample.Description %>%
  filter(sample.description %in% colnames(spring.2018.HF)) %>%
  droplevels(.)
Sample.Description.summer.2018.HF = Sample.Description %>%
  filter(sample.description %in% colnames(summer.2018.HF)) %>%
  droplevels(.)
Sample.Description.fall.2018.HF = Sample.Description %>%
  filter(sample.description %in% colnames(fall.2018.HF)) %>%
  droplevels(.)
Sample.Description.spring.2018.SERC = Sample.Description %>%
  filter(sample.description %in% colnames(spring.2018.SERC)) %>%
  droplevels(.)
Sample.Description.summer.2018.SERC = Sample.Description %>%
  filter(sample.description %in% colnames(summer.2018.SERC)) %>%
  droplevels(.)
Sample.Description.fall.2018.SERC = Sample.Description %>%
  filter(sample.description %in% colnames(fall.2018.SERC)) %>%
  droplevels(.)

#### Create DGE data ####
# Create a numeric matrix of our count data to serve as an input to DGElist.
Data.matrix.spring.2017 <- spring.2017 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.spring.2017) <- spring.2017$Gene_ID
Data.matrix.summer.2017 <- summer.2017 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.summer.2017) <- summer.2017$Gene_ID
Data.matrix.fall.2017 <- fall.2017 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.fall.2017) <- fall.2017$Gene_ID

Data.matrix.spring.2018 <- spring.2018 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.spring.2018) <- spring.2018$Gene_ID
Data.matrix.summer.2018 <- summer.2018 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.summer.2018) <- summer.2018$Gene_ID
Data.matrix.fall.2018 <- fall.2018 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.fall.2018) <- fall.2018$Gene_ID

Data.matrix.spring.2019 <- spring.2019 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.spring.2019) <- spring.2019$Gene_ID
Data.matrix.summer.2019 <- summer.2019 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.summer.2019) <- summer.2019$Gene_ID
Data.matrix.fall.2019 <- fall.2019 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.fall.2019) <- fall.2019$Gene_ID

Data.matrix.spring.2020 <- spring.2020 %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.spring.2020) <- spring.2020$Gene_ID

Data.matrix.spring.2017.HF <- spring.2017.HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.spring.2017.HF) <- spring.2017.HF$Gene_ID
Data.matrix.summer.2017.HF <- summer.2017.HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.summer.2017.HF) <- summer.2017.HF$Gene_ID
Data.matrix.fall.2017.HF <- fall.2017.HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.fall.2017.HF) <- fall.2017.HF$Gene_ID
Data.matrix.spring.2017.SERC <- spring.2017.SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.spring.2017.SERC) <- spring.2017.SERC$Gene_ID
Data.matrix.summer.2017.SERC <- summer.2017.SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.summer.2017.SERC) <- summer.2017.SERC$Gene_ID
Data.matrix.fall.2017.SERC <- fall.2017.SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.fall.2017.SERC) <- fall.2017.SERC$Gene_ID

Data.matrix.spring.2018.HF <- spring.2018.HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.spring.2018.HF) <- spring.2018.HF$Gene_ID
Data.matrix.summer.2018.HF <- summer.2018.HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.summer.2018.HF) <- summer.2018.HF$Gene_ID
Data.matrix.fall.2018.HF <- fall.2018.HF %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.fall.2018.HF) <- fall.2018.HF$Gene_ID

Data.matrix.spring.2018.SERC <- spring.2018.SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.spring.2018.SERC) <- spring.2018.SERC$Gene_ID
Data.matrix.summer.2018.SERC <- summer.2018.SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.summer.2018.SERC) <- summer.2018.SERC$Gene_ID
Data.matrix.fall.2018.SERC <- fall.2018.SERC %>%
  select(-Gene_ID) %>%
  as.matrix()
rownames(Data.matrix.fall.2018.SERC) <- fall.2018.SERC$Gene_ID

#### Make DGEList ####
# The DGEList function creates a DGElist object for differential expression analysis. 
# we specify counts to be the data matrix we created above that holds our count data. 
# Group is specified from the grow.nogrow variable

DGE.data.spring.2017 = DGEList(counts = Data.matrix.spring.2017, group = Sample.Description.spring.2017$growth.signal)
DGE.data.summer.2017 = DGEList(counts = Data.matrix.summer.2017, group = Sample.Description.summer.2017$growth.signal)
DGE.data.fall.2017 = DGEList(counts = Data.matrix.fall.2017, group = Sample.Description.fall.2017$growth.signal)

DGE.data.spring.2018 = DGEList(counts = Data.matrix.spring.2018, group = Sample.Description.spring.2018$growth.signal)
DGE.data.summer.2018 = DGEList(counts = Data.matrix.summer.2018, group = Sample.Description.summer.2018$growth.signal)
DGE.data.fall.2018 = DGEList(counts = Data.matrix.fall.2018, group = Sample.Description.fall.2018$growth.signal)

DGE.data.spring.2019 = DGEList(counts = Data.matrix.spring.2019, group = Sample.Description.spring.2019$growth.signal)
DGE.data.summer.2019 = DGEList(counts = Data.matrix.summer.2019, group = Sample.Description.summer.2019$growth.signal)
DGE.data.fall.2019 = DGEList(counts = Data.matrix.fall.2019, group = Sample.Description.fall.2019$growth.signal)

DGE.data.spring.2020 = DGEList(counts = Data.matrix.spring.2020, group = Sample.Description.spring.2020$growth.signal)

DGE.data.spring.2017.HF = DGEList(counts = Data.matrix.spring.2017.HF, group = Sample.Description.spring.2017.HF$growth.signal)
DGE.data.summer.2017.HF = DGEList(counts = Data.matrix.summer.2017.HF, group = Sample.Description.summer.2017.HF$growth.signal)
DGE.data.fall.2017.HF = DGEList(counts = Data.matrix.fall.2017.HF, group = Sample.Description.fall.2017.HF$growth.signal)

DGE.data.spring.2017.SERC = DGEList(counts = Data.matrix.spring.2017.SERC, group = Sample.Description.spring.2017.SERC$growth.signal)
DGE.data.summer.2017.SERC = DGEList(counts = Data.matrix.summer.2017.SERC, group = Sample.Description.summer.2017.SERC$growth.signal)
DGE.data.fall.2017.SERC = DGEList(counts = Data.matrix.fall.2017.SERC, group = Sample.Description.fall.2017.SERC$growth.signal)

DGE.data.spring.2018.HF = DGEList(counts = Data.matrix.spring.2018.HF, group = Sample.Description.spring.2018.HF$growth.signal)
DGE.data.summer.2018.HF = DGEList(counts = Data.matrix.summer.2018.HF, group = Sample.Description.summer.2018.HF$growth.signal)
DGE.data.fall.2018.HF = DGEList(counts = Data.matrix.fall.2018.HF, group = Sample.Description.fall.2018.HF$growth.signal)

DGE.data.spring.2018.SERC = DGEList(counts = Data.matrix.spring.2018.SERC, group = Sample.Description.spring.2018.SERC$growth.signal)
DGE.data.summer.2018.SERC = DGEList(counts = Data.matrix.summer.2018.SERC, group = Sample.Description.summer.2018.SERC$growth.signal)
DGE.data.fall.2018.SERC = DGEList(counts = Data.matrix.fall.2018.SERC, group = Sample.Description.fall.2018.SERC$growth.signal)

#### Normalize the data #####
# using the TMM method, trimmed mean of M-values
DGE.data.spring.2017 <- edgeR::calcNormFactors(DGE.data.spring.2017, method = "TMM")
DGE.data.summer.2017 <- calcNormFactors(DGE.data.summer.2017, method = "TMM")
DGE.data.fall.2017 <- calcNormFactors(DGE.data.fall.2017, method = "TMM")

DGE.data.spring.2018 <- calcNormFactors(DGE.data.spring.2018, method = "TMM")
DGE.data.summer.2018 <- calcNormFactors(DGE.data.summer.2018, method = "TMM")
DGE.data.fall.2018 <- calcNormFactors(DGE.data.fall.2018, method = "TMM")

DGE.data.spring.2019 <- calcNormFactors(DGE.data.spring.2019, method = "TMM")
DGE.data.summer.2019 <- calcNormFactors(DGE.data.summer.2019, method = "TMM")
DGE.data.fall.2019 <- calcNormFactors(DGE.data.fall.2019, method = "TMM")

DGE.data.spring.2020 <- calcNormFactors(DGE.data.spring.2020, method = "TMM")

DGE.data.spring.2017.HF <- calcNormFactors(DGE.data.spring.2017.HF, method = "TMM")
DGE.data.summer.2017.HF <- calcNormFactors(DGE.data.summer.2017.HF, method = "TMM")
DGE.data.fall.2017.HF <- calcNormFactors(DGE.data.fall.2017.HF, method = "TMM")

DGE.data.spring.2017.SERC <- calcNormFactors(DGE.data.spring.2017.SERC, method = "TMM")
DGE.data.summer.2017.SERC <- calcNormFactors(DGE.data.summer.2017.SERC, method = "TMM")
DGE.data.fall.2017.SERC <- calcNormFactors(DGE.data.fall.2017.SERC, method = "TMM")

DGE.data.spring.2018.HF <- calcNormFactors(DGE.data.spring.2018.HF, method = "TMM")
DGE.data.summer.2018.HF <- calcNormFactors(DGE.data.summer.2018.HF, method = "TMM")
DGE.data.fall.2018.HF <- calcNormFactors(DGE.data.fall.2018.HF, method = "TMM")

DGE.data.spring.2018.SERC <- calcNormFactors(DGE.data.spring.2018.SERC, method = "TMM")
DGE.data.summer.2018.SERC <- calcNormFactors(DGE.data.summer.2018.SERC, method = "TMM")
DGE.data.fall.2018.SERC <- calcNormFactors(DGE.data.fall.2018.SERC, method = "TMM")

#### Plotting ####
# bcv is the square root of the dispersion of the negative binomial distribution
plotMDS(DGE.data.spring.2017, method = "bcv") 
plotMDS(DGE.data.summer.2017, method = "bcv") 
plotMDS(DGE.data.fall.2017, method = "bcv") 
plotMDS(DGE.data.spring.2018, method = "bcv") 
plotMDS(DGE.data.summer.2018, method = "bcv") 
plotMDS(DGE.data.fall.2018, method = "bcv") 
plotMDS(DGE.data.spring.2019, method = "bcv") 
plotMDS(DGE.data.summer.2019, method = "bcv") 
plotMDS(DGE.data.fall.2019, method = "bcv") 
plotMDS(DGE.data.spring.2020, method = "bcv")

plotMDS(DGE.data.spring.2017.HF, method = "bcv") 
plotMDS(DGE.data.summer.2017.HF, method = "bcv") 
plotMDS(DGE.data.fall.2017.HF, method = "bcv") 
plotMDS(DGE.data.spring.2017.SERC, method = "bcv") 
plotMDS(DGE.data.summer.2017.SERC, method = "bcv") 
plotMDS(DGE.data.fall.2017.SERC, method = "bcv")

plotMDS(DGE.data.spring.2018.HF, method = "bcv") 
plotMDS(DGE.data.summer.2018.HF, method = "bcv") 
plotMDS(DGE.data.fall.2018.HF, method = "bcv") 
plotMDS(DGE.data.spring.2018.SERC, method = "bcv") 
plotMDS(DGE.data.summer.2018.SERC, method = "bcv") 
plotMDS(DGE.data.fall.2018.SERC, method = "bcv")

#### TMM Data Saving ####
TMM.DATA.spring.2017 <- cpm(DGE.data.spring.2017, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.spring.2017,"./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2017.csv")
TMM.DATA.summer.2017 <- cpm(DGE.data.summer.2017, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.summer.2017,"./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2017.csv") 
TMM.DATA.fall.2017 <- cpm(DGE.data.fall.2017, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.fall.2017,"./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2017.csv") 

TMM.DATA.spring.2018 <- cpm(DGE.data.spring.2018, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.spring.2018,"./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2018.csv")
TMM.DATA.summer.2018 <- cpm(DGE.data.summer.2018, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.summer.2018,"./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2018.csv") 
TMM.DATA.fall.2018 <- cpm(DGE.data.fall.2018, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.fall.2018,"./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2018.csv") 

TMM.DATA.spring.2019 <- cpm(DGE.data.spring.2019, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.spring.2019,"./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2019.csv")
TMM.DATA.summer.2019 <- cpm(DGE.data.summer.2019, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.summer.2019,"./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2019.csv") 
TMM.DATA.fall.2019 <- cpm(DGE.data.fall.2019, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.fall.2019,"./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2019.csv") 

TMM.DATA.spring.2020 <- cpm(DGE.data.spring.2020, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.spring.2020,"./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2020.csv")

TMM.DATA.spring.2017.HF <- cpm(DGE.data.spring.2017.HF, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.spring.2017.HF,"./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2017.HF.csv")
TMM.DATA.summer.2017.HF <- cpm(DGE.data.summer.2017.HF, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.summer.2017.HF,"./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2017.HF.csv") 
TMM.DATA.fall.2017.HF <- cpm(DGE.data.fall.2017.HF, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.fall.2017.HF,"./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2017.HF.csv") 
TMM.DATA.spring.2017.SERC <- cpm(DGE.data.spring.2017.SERC, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.spring.2017.SERC,"./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2017.SERC.csv")
TMM.DATA.summer.2017.SERC <- cpm(DGE.data.summer.2017.SERC, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.summer.2017.SERC,"./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2017.SERC.csv") 
TMM.DATA.fall.2017.SERC <- cpm(DGE.data.fall.2017.SERC, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.fall.2017.SERC,"./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2017.SERC.csv") 

TMM.DATA.spring.2018.HF <- cpm(DGE.data.spring.2018.HF, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.spring.2018.HF,"./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2018.HF.csv")
TMM.DATA.summer.2018.HF <- cpm(DGE.data.summer.2018.HF, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.summer.2018.HF,"./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2018.HF.csv") 
TMM.DATA.fall.2018.HF <- cpm(DGE.data.fall.2018.HF, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.fall.2018.HF,"./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2018.HF.csv") 
TMM.DATA.spring.2018.SERC <- cpm(DGE.data.spring.2018.SERC, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.spring.2018.SERC,"./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2018.SERC.csv")
TMM.DATA.summer.2018.SERC <- cpm(DGE.data.summer.2018.SERC, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.summer.2018.SERC,"./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2018.SERC.csv") 
TMM.DATA.fall.2018.SERC <- cpm(DGE.data.fall.2018.SERC, log = T) %>%
  data.frame() %>%
  rownames_to_column(var = "Gene_ID")
#write_csv(TMM.DATA.fall.2018.SERC,"./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2018.SERC.csv") 

#### Create design matrix ####
#Create our design matrix which considers every interaction between species and treatment, in this case Grow/NoGrow
Design.spring.2017 <-
  model.matrix(~ grow.nogrow, data = Sample.Description.spring.2017) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2017) <- Sample.Description.spring.2017$sample.description
Design.summer.2017 <-
  model.matrix(~ grow.nogrow, data = Sample.Description.summer.2017) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.summer.2017) <- Sample.Description.summer.2017$sample.description
Design.fall.2017 <-
  model.matrix(~ grow.nogrow, data = Sample.Description.fall.2017) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.fall.2017) <- Sample.Description.fall.2017$sample.description

Design.spring.2018 <-
  model.matrix(~ grow.nogrow, data = Sample.Description.spring.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2018) <- Sample.Description.spring.2018$sample.description
Design.summer.2018 <-
  model.matrix(~ grow.nogrow, data = Sample.Description.summer.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.summer.2018) <- Sample.Description.summer.2018$sample.description
Design.fall.2018 <-
  model.matrix(~ grow.nogrow, data = Sample.Description.fall.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.fall.2018) <- Sample.Description.fall.2018$sample.description

Design.spring.2019 <-
  model.matrix(~ grow.nogrow, data = Sample.Description.spring.2019) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2019) <- Sample.Description.spring.2019$sample.description
Design.summer.2019 <-
  model.matrix(~ grow.nogrow, data = Sample.Description.summer.2019) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.summer.2019) <- Sample.Description.summer.2019$sample.description
Design.fall.2019 <-
  model.matrix(~ grow.nogrow, data = Sample.Description.fall.2019) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.fall.2019) <- Sample.Description.fall.2019$sample.description

Design.spring.2020 <-
  model.matrix(~ grow.nogrow, data = Sample.Description.spring.2020) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2020) <- Sample.Description.spring.2020$sample.description

Design.spring.2017.HF <-
  model.matrix(~ grow.nogrow, data = Sample.Description.spring.2017.HF) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2017.HF) <- Sample.Description.spring.2017.HF$sample.description
Design.summer.2017.HF <-
  model.matrix(~ grow.nogrow, data = Sample.Description.summer.2017.HF) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.summer.2017.HF) <- Sample.Description.summer.2017.HF$sample.description
Design.fall.2017.HF <-
  model.matrix(~ grow.nogrow, data = Sample.Description.fall.2017.HF) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.fall.2017.HF) <- Sample.Description.fall.2017.HF$sample.description
Design.spring.2017.SERC <-
  model.matrix(~ grow.nogrow, data = Sample.Description.spring.2017.SERC) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2017.SERC) <- Sample.Description.spring.2017.SERC$sample.description
Design.summer.2017.SERC <-
  model.matrix(~ grow.nogrow, data = Sample.Description.summer.2017.SERC) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.summer.2017.SERC) <- Sample.Description.summer.2017.SERC$sample.description
Design.fall.2017.SERC <-
  model.matrix(~ grow.nogrow, data = Sample.Description.fall.2017.SERC) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.fall.2017.SERC) <- Sample.Description.fall.2017.SERC$sample.description

Design.spring.2018.HF <-
  model.matrix(~ grow.nogrow, data = Sample.Description.spring.2018.HF) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2018.HF) <- Sample.Description.spring.2018.HF$sample.description
Design.summer.2018.HF <-
  model.matrix(~ grow.nogrow, data = Sample.Description.summer.2018.HF) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.summer.2018.HF) <- Sample.Description.summer.2018.HF$sample.description
Design.fall.2018.HF <-
  model.matrix(~ grow.nogrow, data = Sample.Description.fall.2018.HF) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.fall.2018.HF) <- Sample.Description.fall.2018.HF$sample.description
Design.spring.2018.SERC <-
  model.matrix(~ grow.nogrow, data = Sample.Description.spring.2018.SERC) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2018.SERC) <- Sample.Description.spring.2018.SERC$sample.description
Design.summer.2018.SERC <-
  model.matrix(~ grow.nogrow, data = Sample.Description.summer.2018.SERC) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.summer.2018.SERC) <- Sample.Description.summer.2018.SERC$sample.description
Design.fall.2018.SERC <-
  model.matrix(~ grow.nogrow, data = Sample.Description.fall.2018.SERC) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.fall.2018.SERC) <- Sample.Description.fall.2018.SERC$sample.description

#### Additive design matrix: grow.nogrow + site ####
# testing additive and multiplicative year/site models
Design.spring.2017.additive <-
  model.matrix(~ grow.nogrow + Site, data = Sample.Description.spring.2017) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2017.additive) <- Sample.Description.spring.2017$sample.description
Design.summer.2017.additive <-
  model.matrix(~ grow.nogrow + Site, data = Sample.Description.summer.2017) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.summer.2017.additive) <- Sample.Description.summer.2017$sample.description
Design.fall.2017.additive <-
  model.matrix(~ grow.nogrow + Site, data = Sample.Description.fall.2017) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.fall.2017.additive) <- Sample.Description.fall.2017$sample.description

Design.spring.2018.additive <-
  model.matrix(~ grow.nogrow + Site, data = Sample.Description.spring.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2018.additive) <- Sample.Description.spring.2018$sample.description
Design.summer.2018.additive <-
  model.matrix(~ grow.nogrow + Site, data = Sample.Description.summer.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.summer.2018.additive) <- Sample.Description.summer.2018$sample.description
Design.fall.2018.additive <-
  model.matrix(~ grow.nogrow + Site, data = Sample.Description.fall.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.fall.2018.additive) <- Sample.Description.fall.2018$sample.description

#### Multiplicative design matrix: grow.nogrow*site ####
# sample size probably too small
Design.spring.2017.mulitply <-
  model.matrix(~ grow.nogrow*Site, data = Sample.Description.spring.2017) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2017.mulitply) <- Sample.Description.spring.2017$sample.description
Design.summer.2017.mulitply <-
  model.matrix(~ grow.nogrow*Site, data = Sample.Description.summer.2017) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.summer.2017.mulitply) <- Sample.Description.summer.2017$sample.description
Design.fall.2017.mulitply <-
  model.matrix(~ grow.nogrow*Site, data = Sample.Description.fall.2017) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.fall.2017.mulitply) <- Sample.Description.fall.2017$sample.description

Design.spring.2018.mulitply <-
  model.matrix(~ grow.nogrow*Site, data = Sample.Description.spring.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.spring.2018.mulitply) <- Sample.Description.spring.2018$sample.description
Design.summer.2018.mulitply <-
  model.matrix(~ grow.nogrow*Site, data = Sample.Description.summer.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.summer.2018.mulitply) <- Sample.Description.summer.2018$sample.description
Design.fall.2018.mulitply <-
  model.matrix(~ grow.nogrow*Site, data = Sample.Description.fall.2018) %>%
  as.data.frame() %>%
  as.matrix()
rownames(Design.fall.2018.mulitply) <- Sample.Description.fall.2018$sample.description

#### Calculate Dispersion Factors ####
# To estimate common dispersion:
DGE.data.spring.2017 <- estimateGLMCommonDisp(DGE.data.spring.2017, Design.spring.2017, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2017 <- estimateGLMTrendedDisp(GE.data.spring.2017, Design.spring.2017)
#To estimate tagwise dispersions:
DGE.data.spring.2017 <- estimateGLMTagwiseDisp(GE.data.spring.2017, Design.spring.2017)
plotBCV(DGE.data.spring.2017)
# To estimate common dispersion:
DGE.data.summer.2017 <- estimateGLMCommonDisp(DGE.data.summer.2017, Design.summer.2017, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.summer.2017 <- estimateGLMTrendedDisp(GE.data.summer.2017, Design.summer.2017)
#To estimate tagwise dispersions:
DGE.data.summer.2017 <- estimateGLMTagwiseDisp(GE.data.summer.2017, Design.summer.2017)
plotBCV(DGE.data.summer.2017)
# To estimate common dispersion:
DGE.data.fall.2017 <- estimateGLMCommonDisp(DGE.data.fall.2017, Design.fall.2017, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.fall.2017 <- estimateGLMTrendedDisp(GE.data.fall.2017, Design.fall.2017)
#To estimate tagwise dispersions:
DGE.data.fall.2017 <- estimateGLMTagwiseDisp(GE.data.fall.2017, Design.fall.2017)
plotBCV(DGE.data.fall.2017)

# To estimate common dispersion:
DGE.data.spring.2018 <- estimateGLMCommonDisp(DGE.data.spring.2018, Design.spring.2018, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2018 <- estimateGLMTrendedDisp(GE.data.spring.2018, Design.spring.2018)
#To estimate tagwise dispersions:
DGE.data.spring.2018 <- estimateGLMTagwiseDisp(GE.data.spring.2018, Design.spring.2018)
plotBCV(DGE.data.spring.2018)
# To estimate common dispersion:
DGE.data.summer.2018 <- estimateGLMCommonDisp(DGE.data.summer.2018, Design.summer.2018, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.summer.2018 <- estimateGLMTrendedDisp(GE.data.summer.2018, Design.summer.2018)
#To estimate tagwise dispersions:
DGE.data.summer.2018 <- estimateGLMTagwiseDisp(GE.data.summer.2018, Design.summer.2018)
plotBCV(DGE.data.summer.2018)
# To estimate common dispersion:
DGE.data.fall.2018 <- estimateGLMCommonDisp(DGE.data.fall.2018, Design.fall.2018, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.fall.2018 <- estimateGLMTrendedDisp(GE.data.fall.2018, Design.fall.2018)
#To estimate tagwise dispersions:
DGE.data.fall.2018 <- estimateGLMTagwiseDisp(GE.data.fall.2018, Design.fall.2018)
plotBCV(DGE.data.fall.2018)

DGE.data.spring.2019 <- estimateGLMCommonDisp(DGE.data.spring.2019, Design.spring.2019, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2019 <- estimateGLMTrendedDisp(GE.data.spring.2019, Design.spring.2019)
#To estimate tagwise dispersions:
DGE.data.spring.2019 <- estimateGLMTagwiseDisp(GE.data.spring.2019, Design.spring.2019)
plotBCV(DGE.data.spring.2019)
# To estimate common dispersion:
DGE.data.summer.2019 <- estimateGLMCommonDisp(DGE.data.summer.2019, Design.summer.2019, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.summer.2019 <- estimateGLMTrendedDisp(GE.data.summer.2019, Design.summer.2019)
#To estimate tagwise dispersions:
DGE.data.summer.2019 <- estimateGLMTagwiseDisp(GE.data.summer.2019, Design.summer.2019)
plotBCV(DGE.data.summer.2019)
# To estimate common dispersion:
DGE.data.fall.2019 <- estimateGLMCommonDisp(DGE.data.fall.2019, Design.fall.2019, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.fall.2019 <- estimateGLMTrendedDisp(GE.data.fall.2019, Design.fall.2019)
#To estimate tagwise dispersions:
DGE.data.fall.2019 <- estimateGLMTagwiseDisp(GE.data.fall.2019, Design.fall.2019)
plotBCV(DGE.data.fall.2019)

# To estimate common dispersion:
DGE.data.spring.2020 <- estimateGLMCommonDisp(DGE.data.spring.2020, Design.spring.2020, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2020 <- estimateGLMTrendedDisp(GE.data.spring.2020, Design.spring.2020)
#To estimate tagwise dispersions:
DGE.data.spring.2020 <- estimateGLMTagwiseDisp(GE.data.spring.2020, Design.spring.2020)
plotBCV(DGE.data.spring.2020)

DGE.data.spring.2017.HF <- estimateGLMCommonDisp(DGE.data.spring.2017.HF, Design.spring.2017.HF, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2017.HF <- estimateGLMTrendedDisp(GE.data.spring.2017.HF, Design.spring.2017.HF)
#To estimate tagwise dispersions:
DGE.data.spring.2017.HF <- estimateGLMTagwiseDisp(GE.data.spring.2017.HF, Design.spring.2017.HF)
plotBCV(DGE.data.spring.2017.HF)
# To estimate common dispersion:
DGE.data.summer.2017.HF <- estimateGLMCommonDisp(DGE.data.summer.2017.HF, Design.summer.2017.HF, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.summer.2017.HF <- estimateGLMTrendedDisp(GE.data.summer.2017.HF, Design.summer.2017.HF)
#To estimate tagwise dispersions:
DGE.data.summer.2017.HF <- estimateGLMTagwiseDisp(GE.data.summer.2017.HF, Design.summer.2017.HF)
plotBCV(DGE.data.summer.2017.HF)
# To estimate common dispersion:
DGE.data.fall.2017.HF <- estimateGLMCommonDisp(DGE.data.fall.2017.HF, Design.fall.2017.HF, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.fall.2017.HF <- estimateGLMTrendedDisp(GE.data.fall.2017.HF, Design.fall.2017.HF)
#To estimate tagwise dispersions:
DGE.data.fall.2017.HF <- estimateGLMTagwiseDisp(GE.data.fall.2017.HF, Design.fall.2017.HF)
plotBCV(DGE.data.fall.2017.HF)

# To estimate common dispersion:
DGE.data.spring.2018.HF <- estimateGLMCommonDisp(DGE.data.spring.2018.HF, Design.spring.2018.HF, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2018.HF <- estimateGLMTrendedDisp(GE.data.spring.2018.HF, Design.spring.2018.HF)
#To estimate tagwise dispersions:
DGE.data.spring.2018.HF <- estimateGLMTagwiseDisp(GE.data.spring.2018.HF, Design.spring.2018.HF)
plotBCV(DGE.data.spring.2018.HF)
# To estimate common dispersion:
DGE.data.summer.2018.HF <- estimateGLMCommonDisp(DGE.data.summer.2018.HF, Design.summer.2018.HF, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.summer.2018.HF <- estimateGLMTrendedDisp(GE.data.summer.2018.HF, Design.summer.2018.HF)
#To estimate tagwise dispersions:
DGE.data.summer.2018.HF <- estimateGLMTagwiseDisp(GE.data.summer.2018.HF, Design.summer.2018.HF)
plotBCV(DGE.data.summer.2018.HF)
# To estimate common dispersion:
DGE.data.fall.2018.HF <- estimateGLMCommonDisp(DGE.data.fall.2018.HF, Design.fall.2018.HF, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.fall.2018.HF <- estimateGLMTrendedDisp(GE.data.fall.2018.HF, Design.fall.2018.HF)
#To estimate tagwise dispersions:
DGE.data.fall.2018.HF <- estimateGLMTagwiseDisp(GE.data.fall.2018.HF, Design.fall.2018.HF)
plotBCV(DGE.data.fall.2018.HF)

DGE.data.spring.2017.SERC <- estimateGLMCommonDisp(DGE.data.spring.2017.SERC, Design.spring.2017.SERC, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2017.SERC <- estimateGLMTrendedDisp(GE.data.spring.2017.SERC, Design.spring.2017.SERC)
#To estimate tagwise dispersions:
DGE.data.spring.2017.SERC <- estimateGLMTagwiseDisp(GE.data.spring.2017.SERC, Design.spring.2017.SERC)
plotBCV(DGE.data.spring.2017.SERC)
# To estimate common dispersion:
DGE.data.summer.2017.SERC <- estimateGLMCommonDisp(DGE.data.summer.2017.SERC, Design.summer.2017.SERC, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.summer.2017.SERC <- estimateGLMTrendedDisp(GE.data.summer.2017.SERC, Design.summer.2017.SERC)
#To estimate tagwise dispersions:
DGE.data.summer.2017.SERC <- estimateGLMTagwiseDisp(GE.data.summer.2017.SERC, Design.summer.2017.SERC)
plotBCV(DGE.data.summer.2017.SERC)
# To estimate common dispersion:
DGE.data.fall.2017.SERC <- estimateGLMCommonDisp(DGE.data.fall.2017.SERC, Design.fall.2017.SERC, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.fall.2017.SERC <- estimateGLMTrendedDisp(GE.data.fall.2017.SERC, Design.fall.2017.SERC)
#To estimate tagwise dispersions:
DGE.data.fall.2017.SERC <- estimateGLMTagwiseDisp(GE.data.fall.2017.SERC, Design.fall.2017.SERC)
plotBCV(DGE.data.fall.2017.SERC)

# To estimate common dispersion:
DGE.data.spring.2018.SERC <- estimateGLMCommonDisp(DGE.data.spring.2018.SERC, Design.spring.2018.SERC, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2018.SERC <- estimateGLMTrendedDisp(GE.data.spring.2018.SERC, Design.spring.2018.SERC)
#To estimate tagwise dispersions:
DGE.data.spring.2018.SERC <- estimateGLMTagwiseDisp(GE.data.spring.2018.SERC, Design.spring.2018.SERC)
plotBCV(DGE.data.spring.2018.SERC)
# To estimate common dispersion:
DGE.data.summer.2018.SERC <- estimateGLMCommonDisp(DGE.data.summer.2018.SERC, Design.summer.2018.SERC, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.summer.2018.SERC <- estimateGLMTrendedDisp(GE.data.summer.2018.SERC, Design.summer.2018.SERC)
#To estimate tagwise dispersions:
DGE.data.summer.2018.SERC <- estimateGLMTagwiseDisp(GE.data.summer.2018.SERC, Design.summer.2018.SERC)
plotBCV(DGE.data.summer.2018.SERC)
# To estimate common dispersion:
DGE.data.fall.2018.SERC <- estimateGLMCommonDisp(DGE.data.fall.2018.SERC, Design.fall.2018.SERC, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.fall.2018.SERC <- estimateGLMTrendedDisp(GE.data.fall.2018.SERC, Design.fall.2018.SERC)
#To estimate tagwise dispersions:
DGE.data.fall.2018.SERC <- estimateGLMTagwiseDisp(GE.data.fall.2018.SERC, Design.fall.2018.SERC)
plotBCV(DGE.data.fall.2018.SERC)

#### additive dispersion factors ####
# additive design
DGE.data.spring.2017.additive <- estimateGLMCommonDisp(DGE.data.spring.2017, Design.spring.2017.additive, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2017.additive <- estimateGLMTrendedDisp(DGE.data.spring.2017, Design.spring.2017.additive)
#To estimate tagwise dispersions:
DGE.data.spring.2017.additive <- estimateGLMTagwiseDisp(DGE.data.spring.2017, Design.spring.2017.additive)
plotBCV(DGE.data.spring.2017.additive)
DGE.data.summer.2017.additive <- estimateGLMCommonDisp(DGE.data.summer.2017, Design.summer.2017.additive, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.summer.2017.additive <- estimateGLMTrendedDisp(DGE.data.summer.2017, Design.summer.2017.additive)
#To estimate tagwise dispersions:
DGE.data.summer.2017.additive <- estimateGLMTagwiseDisp(DGE.data.summer.2017, Design.summer.2017.additive)
plotBCV(DGE.data.summer.2017.additive)
DGE.data.fall.2017.additive <- estimateGLMCommonDisp(DGE.data.fall.2017, Design.fall.2017.additive, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.fall.2017.additive <- estimateGLMTrendedDisp(DGE.data.fall.2017, Design.fall.2017.additive)
#To estimate tagwise dispersions:
DGE.data.fall.2017.additive <- estimateGLMTagwiseDisp(DGE.data.fall.2017, Design.fall.2017.additive)
plotBCV(DGE.data.fall.2017.additive)

DGE.data.spring.2018.additive <- estimateGLMCommonDisp(DGE.data.spring.2018, Design.spring.2018.additive, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2018.additive <- estimateGLMTrendedDisp(DGE.data.spring.2018, Design.spring.2018.additive)
#To estimate tagwise dispersions:
DGE.data.spring.2018.additive <- estimateGLMTagwiseDisp(DGE.data.spring.2018, Design.spring.2018.additive)
plotBCV(DGE.data.spring.2018.additive)
DGE.data.summer.2018.additive <- estimateGLMCommonDisp(DGE.data.summer.2018, Design.summer.2018.additive, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.summer.2018.additive <- estimateGLMTrendedDisp(DGE.data.summer.2018, Design.summer.2018.additive)
#To estimate tagwise dispersions:
DGE.data.summer.2018.additive <- estimateGLMTagwiseDisp(DGE.data.summer.2018, Design.summer.2018.additive)
plotBCV(DGE.data.summer.2018.additive)
DGE.data.fall.2018.additive <- estimateGLMCommonDisp(DGE.data.fall.2018, Design.fall.2018.additive, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.fall.2018.additive <- estimateGLMTrendedDisp(DGE.data.fall.2018, Design.fall.2018.additive)
#To estimate tagwise dispersions:
DGE.data.fall.2018.additive <- estimateGLMTagwiseDisp(DGE.data.fall.2018, Design.fall.2018.additive)
plotBCV(DGE.data.fall.2018.additive)

#### multiplicative dispersion factors ####
DGE.data.spring.2017.multiply <- estimateGLMCommonDisp(DGE.data.spring.2017, Design.spring.2017.multiply, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2017.multiply <- estimateGLMTrendedDisp(DGE.data.spring.2017, Design.spring.2017.multiply)
#To estimate tagwise dispersions:
DGE.data.spring.2017.multiply <- estimateGLMTagwiseDisp(DGE.data.spring.2017, Design.spring.2017.multiply)
plotBCV(DGE.data.spring.2017.multiply)
DGE.data.summer.2017.multiply <- estimateGLMCommonDisp(DGE.data.summer.2017, Design.summer.2017.multiply, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.summer.2017.multiply <- estimateGLMTrendedDisp(DGE.data.summer.2017, Design.summer.2017.multiply)
#To estimate tagwise dispersions:
DGE.data.summer.2017.multiply <- estimateGLMTagwiseDisp(DGE.data.summer.2017, Design.summer.2017.multiply)
plotBCV(DGE.data.summer.2017.multiply)
DGE.data.fall.2017.multiply <- estimateGLMCommonDisp(DGE.data.fall.2017, Design.fall.2017.multiply, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.fall.2017.multiply <- estimateGLMTrendedDisp(DGE.data.fall.2017, Design.fall.2017.multiply)
#To estimate tagwise dispersions:
DGE.data.fall.2017.multiply <- estimateGLMTagwiseDisp(DGE.data.fall.2017, Design.fall.2017.multiply)
plotBCV(DGE.data.fall.2017.multiply)

DGE.data.spring.2018.multiply <- estimateGLMCommonDisp(DGE.data.spring.2018, Design.spring.2018.multiply, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.spring.2018.multiply <- estimateGLMTrendedDisp(DGE.data.spring.2018, Design.spring.2018.multiply)
#To estimate tagwise dispersions:
DGE.data.spring.2018.multiply <- estimateGLMTagwiseDisp(DGE.data.spring.2018, Design.spring.2018.multiply)
plotBCV(DGE.data.spring.2018.multiply)
DGE.data.summer.2018.multiply <- estimateGLMCommonDisp(DGE.data.summer.2018, Design.summer.2018.multiply, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.summer.2018.multiply <- estimateGLMTrendedDisp(DGE.data.summer.2018, Design.summer.2018.multiply)
#To estimate tagwise dispersions:
DGE.data.summer.2018.multiply <- estimateGLMTagwiseDisp(DGE.data.summer.2018, Design.summer.2018.multiply)
plotBCV(DGE.data.summer.2018.multiply)
DGE.data.fall.2018.multiply <- estimateGLMCommonDisp(DGE.data.fall.2018, Design.fall.2018.multiply, verbose = TRUE)
#To estimate trended dispersions:
DGE.data.fall.2018.multiply <- estimateGLMTrendedDisp(DGE.data.fall.2018, Design.fall.2018.multiply)
#To estimate tagwise dispersions:
DGE.data.fall.2018.multiply <- estimateGLMTagwiseDisp(DGE.data.fall.2018, Design.fall.2018.multiply)
plotBCV(DGE.data.fall.2018.multiply)

### Plot log transformed TMM normalized count data ####
#Print box plots of the log transformed TMM normalized count data in the DGElist object termed DGE.data
#We can use this to exclude sample outliers with abnormally low expression levels
# repeat code for each sample grouping
jpeg(
  "./Plots/grow.nogrow/DE/spring.2017.normalized.gene.counts.jpeg",
  width = 1920,
  height = 1080
)
Norml.Count.Data <- cpm(DGE.data.spring.2017, log = TRUE)

boxplot(Norml.Count.Data, main = "Count Data after transformation", ylab = "log2(cpm)")
dev.off()

### Fit the model ####
# fit a negative binomial generalized log-linear model to our normalized count data 
# using the full design matrix for gene wise statistical tests.

# CHANGE COEF NAME TO MATCH

fit.spring.2017 <- glmFit(DGE.data.spring.2017,Design.spring.2017)
fit.spring.2017.lrt <- glmLRT(fit.spring.2017,coef = "grow.nogrow")
fit.summer.2017 <- glmFit(DGE.data.summer.2017,Design.summer.2017)
fit.summer.2017.lrt <- glmLRT(fit.summer.2017,coef = "grow.nogrow")
fit.fall.2017 <- glmFit(DGE.data.fall.2017,Design.fall.2017)
fit.fall.2017.lrt <- glmLRT(fit.fall.2017,coef = "grow.nogrow")

fit.spring.2018 <- glmFit(DGE.data.spring.2018,Design.spring.2018)
fit.spring.2018.lrt <- glmLRT(fit.spring.2018,coef = "grow.nogrow")
fit.summer.2018 <- glmFit(DGE.data.summer.2018,Design.summer.2018)
fit.summer.2018.lrt <- glmLRT(fit.summer.2018,coef = "grow.nogrow")
fit.fall.2018 <- glmFit(DGE.data.fall.2018,Design.fall.2018)
fit.fall.2018.lrt <- glmLRT(fit.fall.2018,coef = "grow.nogrow")

fit.spring.2019 <- glmFit(DGE.data.spring.2019,Design.spring.2019)
fit.spring.2019.lrt <- glmLRT(fit.spring.2019,coef = "grow.nogrow")
fit.summer.2019 <- glmFit(DGE.data.summer.2019,Design.summer.2019)
fit.summer.2019.lrt <- glmLRT(fit.summer.2019,coef = "grow.nogrow")
fit.fall.2019 <- glmFit(DGE.data.fall.2019,Design.fall.2019)
fit.fall.2019.lrt <- glmLRT(fit.fall.2019,coef = "grow.nogrow")

fit.spring.2020 <- glmFit(DGE.data.spring.2020,Design.spring.2020)
fit.spring.2020.lrt <- glmLRT(fit.spring.2020,coef = "grow.nogrow")

fit.spring.2017.HF <- glmFit(DGE.data.spring.2017.HF,Design.spring.2017.HF)
fit.spring.2017.HF.lrt <- glmLRT(fit.spring.2017.HF,coef = "grow.nogrow")
fit.summer.2017.HF <- glmFit(DGE.data.summer.2017.HF,Design.summer.2017.HF)
fit.summer.2017.HF.lrt <- glmLRT(fit.summer.2017.HF,coef = "grow.nogrow")
fit.fall.2017.HF <- glmFit(DGE.data.fall.2017.HF,Design.fall.2017.HF)
fit.fall.2017.HF.lrt <- glmLRT(fit.fall.2017.HF,coef = "grow.nogrow")

fit.spring.2018.HF <- glmFit(DGE.data.spring.2018.HF,Design.spring.2018.HF)
fit.spring.2018.HF.lrt <- glmLRT(fit.spring.2018.HF,coef = "grow.nogrow")
fit.summer.2018.HF <- glmFit(DGE.data.summer.2018.HF,Design.summer.2018.HF)
fit.summer.2018.HF.lrt <- glmLRT(fit.summer.2018.HF,coef = "grow.nogrow")
fit.fall.2018.HF <- glmFit(DGE.data.fall.2018.HF,Design.fall.2018.HF)
fit.fall.2018.HF.lrt <- glmLRT(fit.fall.2018.HF,coef = "grow.nogrow")

fit.spring.2017.SERC <- glmFit(DGE.data.spring.2017.SERC,Design.spring.2017.SERC)
fit.spring.2017.SERC.lrt <- glmLRT(fit.spring.2017.SERC,coef = "grow.nogrow")
fit.summer.2017.SERC <- glmFit(DGE.data.summer.2017.SERC,Design.summer.2017.SERC)
fit.summer.2017.SERC.lrt <- glmLRT(fit.summer.2017.SERC,coef = "grow.nogrow")
fit.fall.2017.SERC <- glmFit(DGE.data.fall.2017.SERC,Design.fall.2017.SERC)
fit.fall.2017.SERC.lrt <- glmLRT(fit.fall.2017.SERC,coef = "grow.nogrow")

fit.spring.2018.SERC <- glmFit(DGE.data.spring.2018.SERC,Design.spring.2018.SERC)
fit.spring.2018.SERC.lrt <- glmLRT(fit.spring.2018.SERC,coef = "grow.nogrow")
fit.summer.2018.SERC <- glmFit(DGE.data.summer.2018.SERC,Design.summer.2018.SERC)
fit.summer.2018.SERC.lrt <- glmLRT(fit.summer.2018.SERC,coef = "grow.nogrow")
fit.fall.2018.SERC <- glmFit(DGE.data.fall.2018.SERC,Design.fall.2018.SERC)
fit.fall.2018.SERC.lrt <- glmLRT(fit.fall.2018.SERC,coef = "grow.nogrow")

#### Fit additive model ####
# additive model
fit.spring.2017.additive <- glmFit(DGE.data.spring.2017.additive,Design.spring.2017.additive)
fit.spring.2017.additive.grow.nogrow.lrt <- glmLRT(fit.spring.2017.additive,coef = "grow.nogrow")
fit.spring.2017.additive.site.lrt <- glmLRT(fit.spring.2017.additive,coef = "SiteSERC")
fit.summer.2017.additive <- glmFit(DGE.data.summer.2017.additive,Design.summer.2017.additive)
fit.summer.2017.additive.grow.nogrow.lrt <- glmLRT(fit.summer.2017.additive,coef = "grow.nogrow")
fit.summer.2017.additive.site.lrt <- glmLRT(fit.summer.2017.additive,coef = "SiteSERC")
fit.fall.2017.additive <- glmFit(DGE.data.fall.2017.additive,Design.fall.2017.additive)
fit.fall.2017.additive.grow.nogrow.lrt <- glmLRT(fit.fall.2017.additive,coef = "grow.nogrow")
fit.fall.2017.additive.site.lrt <- glmLRT(fit.fall.2017.additive,coef = "SiteSERC")

fit.spring.2018.additive <- glmFit(DGE.data.spring.2018.additive,Design.spring.2018.additive)
fit.spring.2018.additive.grow.nogrow.lrt <- glmLRT(fit.spring.2018.additive,coef = "grow.nogrow")
fit.spring.2018.additive.site.lrt <- glmLRT(fit.spring.2018.additive,coef = "SiteSERC")
fit.summer.2018.additive <- glmFit(DGE.data.summer.2018.additive,Design.summer.2018.additive)
fit.summer.2018.additive.grow.nogrow.lrt <- glmLRT(fit.summer.2018.additive,coef = "grow.nogrow")
fit.summer.2018.additive.site.lrt <- glmLRT(fit.summer.2018.additive,coef = "SiteSERC")
fit.fall.2018.additive <- glmFit(DGE.data.fall.2018.additive,Design.fall.2018.additive)
fit.fall.2018.additive.grow.nogrow.lrt <- glmLRT(fit.fall.2018.additive,coef = "grow.nogrow")
fit.fall.2018.additive.site.lrt <- glmLRT(fit.fall.2018.additive,coef = "SiteSERC")

#### Fit multiplicative model ####
# multiplicative model
fit.spring.2017.multiply <- glmFit(DGE.data.spring.2017.multiply,Design.spring.2017.mulitply)
fit.spring.2017.multiply.lrt <- glmLRT(fit.spring.2017.multiply,coef = "grow.nogrow:SiteSERC")
fit.summer.2017.multiply <- glmFit(DGE.data.summer.2017.multiply,Design.summer.2017.mulitply)
fit.summer.2017.multiply.lrt <- glmLRT(fit.summer.2017.multiply,coef = "grow.nogrow:SiteSERC")
fit.fall.2017.multiply <- glmFit(DGE.data.fall.2017.multiply,Design.fall.2017.mulitply)
fit.fall.2017.multiply.lrt <- glmLRT(fit.fall.2017.multiply,coef = "grow.nogrow:SiteSERC")

fit.spring.2018.multiply <- glmFit(DGE.data.spring.2018.multiply,Design.spring.2018.mulitply)
fit.spring.2018.multiply.lrt <- glmLRT(fit.spring.2018.multiply,coef = "grow.nogrow:SiteSERC")
fit.summer.2018.multiply <- glmFit(DGE.data.summer.2018.multiply,Design.summer.2018.mulitply)
fit.summer.2018.multiply.lrt <- glmLRT(fit.summer.2018.multiply,coef = "grow.nogrow:SiteSERC")
fit.fall.2018.multiply <- glmFit(DGE.data.fall.2018.multiply,Design.fall.2018.mulitply)
fit.fall.2018.multiply.lrt <- glmLRT(fit.fall.2018.multiply,coef = "grow.nogrow:SiteSERC")




#### Get top expressed genes ####
# logFC is the log2 fold change between grow.nogrow individuals
# logFC of 2 would indicate that the gene is expressed 4 times higher in trees that grew than trees that did not grow 
# logCPM is the average expression across all samples
# LR is the likelihood ratio L(Full Model)/L(small model)
# PValue is unadjusted p-value
# FDR is the false discovery rate (p-value adjusted for multiple testing) # use this one!

topTags(fit.spring.2017.lrt)
topTags(fit.summer.2017.lrt)
topTags(fit.fall.2017.lrt)
topTags(fit.spring.2018.lrt)
topTags(fit.summer.2018.lrt)
topTags(fit.fall.2018.lrt)
topTags(fit.spring.2019.lrt)
topTags(fit.summer.2019.lrt)
topTags(fit.fall.2019.lrt)
topTags(fit.spring.2020.lrt)

topTags(fit.spring.2017.HF.lrt)
topTags(fit.summer.2017.HF.lrt)
topTags(fit.fall.2017.HF.lrt)
topTags(fit.spring.2018.HF.lrt)
topTags(fit.summer.2018.HF.lrt)
topTags(fit.fall.2018.HF.lrt)

topTags(fit.spring.2017.SERC.lrt)
topTags(fit.summer.2017.SERCv)
topTags(fit.fall.2017.SERC.lrt)
topTags(fit.spring.2018.SERC.lrt)
topTags(fit.summer.2018.SERC.lrt)
topTags(fit.fall.2018.SERC.lrt)

# additive models

topTags(fit.spring.2017.additive.grow.nogrow.lrt)
topTags(fit.summer.2017.additive.grow.nogrow.lrt)
topTags(fit.fall.2017.additive.grow.nogrow.lrt)
topTags(fit.spring.2017.additive.site.lrt)
topTags(fit.summer.2017.additive.site.lrt)
topTags(fit.fall.2017.additive.site.lrt)

topTags(fit.spring.2018.additive.grow.nogrow.lrt)
topTags(fit.summer.2018.additive.grow.nogrow.lrt)
topTags(fit.fall.2018.additive.grow.nogrow.lrt)
topTags(fit.spring.2018.additive.site.lrt)
topTags(fit.summer.2018.additive.site.lrt)
topTags(fit.fall.2018.additive.site.lrt)

# multiplicative models

topTags(fit.spring.2017.multiply.lrt)
topTags(fit.summer.2017.multiply.lrt)
topTags(fit.fall.2017.multiply.lrt)

topTags(fit.spring.2018.multiply.lrt)
topTags(fit.summer.2018.multiply.lrt)
topTags(fit.fall.2018.multiply.lrt)

#### Summmary of DGEs ####
#This uses the FDR of 0.01

# number of down and up regulated genes in grow compared to no.grow
summary(decideTestsDGE(fit.spring.2017.lrt,p.value=0.01)) 
# p 0.01 = 1991 Down, 4090 Up 15561 NotSig
summary(decideTestsDGE(fit.summer.2017.lrt,p.value=0.01)) 
# p 0.01 = 4905 Down, 2545 Up 14192 NotSig
summary(decideTestsDGE(fit.fall.2017.lrt,p.value=0.01)) 
# p 0.01 = 2558 Down, 2995 Up 15979 NotSig
# number of down and up regulated genes in grow compared to no.grow
summary(decideTestsDGE(fit.spring.2018.lrt,p.value=0.01)) 
# p 0.01 = 1991 Down, 4090 Up 15561 NotSig
summary(decideTestsDGE(fit.summer.2018.lrt,p.value=0.01)) 
# p 0.01 = 4905 Down, 2545 Up 14192 NotSig
summary(decideTestsDGE(fit.fall.2018.lrt,p.value=0.01)) 
# p 0.01 = 2558 Down, 2995 Up 15979 NotSig
# number of down and up regulated genes in grow compared to no.grow
summary(decideTestsDGE(fit.spring.2019.lrt,p.value=0.01)) 
# p 0.01 = 1991 Down, 4090 Up 15561 NotSig
summary(decideTestsDGE(fit.summer.2019.lrt,p.value=0.01)) 
# p 0.01 = 4905 Down, 2545 Up 14192 NotSig
summary(decideTestsDGE(fit.fall.2019.lrt,p.value=0.01)) 
# p 0.01 = 2558 Down, 2995 Up 15979 NotSig
summary(decideTestsDGE(fit.spring.2020.lrt,p.value=0.01)) 
# p 0.01 = 1991 Down, 4090 Up 15561 NotSig

summary(decideTestsDGE(fit.spring.2017.HF.lrt,p.value=0.01)) 
# p 0.01 = 1991 Down, 4090 Up 15561 NotSig
summary(decideTestsDGE(fit.summer.2017.HF.lrt,p.value=0.01)) 
# p 0.01 = 4905 Down, 2545 Up 14192 NotSig
summary(decideTestsDGE(fit.fall.2017.HF.lrt,p.value=0.01)) 
# p 0.01 = 2558 Down, 2995 Up 15979 NotSig
# number of down and up regulated genes in grow compared to no.grow
summary(decideTestsDGE(fit.spring.2018.HF.lrt,p.value=0.01)) 
# p 0.01 = 1991 Down, 4090 Up 15561 NotSig
summary(decideTestsDGE(fit.summer.2018.HF.lrt,p.value=0.01)) 
# p 0.01 = 4905 Down, 2545 Up 14192 NotSig
summary(decideTestsDGE(fit.fall.2018.HF.lrt,p.value=0.01)) 
# p 0.01 = 2558 Down, 2995 Up 15979 NotSig

summary(decideTestsDGE(fit.spring.2017.SERC.lrt,p.value=0.01)) 
# p 0.01 = 1991 Down, 4090 Up 15561 NotSig
summary(decideTestsDGE(fit.summer.2017.SERC.lrt,p.value=0.01)) 
# p 0.01 = 4905 Down, 2545 Up 14192 NotSig
summary(decideTestsDGE(fit.fall.2017.SERC.lrt,p.value=0.01)) 
# p 0.01 = 2558 Down, 2995 Up 15979 NotSig
# number of down and up regulated genes in grow compared to no.grow
summary(decideTestsDGE(fit.spring.2018.SERC.lrt,p.value=0.01)) 
# p 0.01 = 1991 Down, 4090 Up 15561 NotSig
summary(decideTestsDGE(fit.summer.2018.SERC.lrt,p.value=0.01)) 
# p 0.01 = 4905 Down, 2545 Up 14192 NotSig
summary(decideTestsDGE(fit.fall.2018.SERC.lrt,p.value=0.01)) 
# p 0.01 = 2558 Down, 2995 Up 15979 NotSig

# additive model
summary(decideTestsDGE(fit.spring.2017.additive.site.lrt,p.value=0.01)) 
# p 0.01 = 2004 Down, 3979 Up 15659 NotSig
summary(decideTestsDGE(fit.spring.2017.additive.grow.nogrow.lrt,p.value=0.01)) 
# p 0.01 = 1989 Down, 1730 Up 17923 NotSig
summary(decideTestsDGE(fit.summer.2017.additive.site.lrt,p.value=0.01)) 
# p 0.01 = 2004 Down, 3979 Up 15659 NotSig
summary(decideTestsDGE(fit.summer.2017.additive.grow.nogrow.lrt,p.value=0.01)) 
# p 0.01 = 1989 Down, 1730 Up 17923 NotSig
summary(decideTestsDGE(fit.fall.2017.additive.site.lrt,p.value=0.01)) 
# p 0.01 = 2004 Down, 3979 Up 15659 NotSig
summary(decideTestsDGE(fit.fall.2017.additive.grow.nogrow.lrt,p.value=0.01)) 
# p 0.01 = 1989 Down, 1730 Up 17923 NotSig

summary(decideTestsDGE(fit.spring.2018.additive.site.lrt,p.value=0.01)) 
# p 0.01 = 2004 Down, 3979 Up 15659 NotSig
summary(decideTestsDGE(fit.spring.2018.additive.grow.nogrow.lrt,p.value=0.01)) 
# p 0.01 = 1989 Down, 1730 Up 17923 NotSig
summary(decideTestsDGE(fit.summer.2018.additive.site.lrt,p.value=0.01)) 
# p 0.01 = 2004 Down, 3979 Up 15659 NotSig
summary(decideTestsDGE(fit.summer.2018.additive.grow.nogrow.lrt,p.value=0.01)) 
# p 0.01 = 1989 Down, 1730 Up 17923 NotSig
summary(decideTestsDGE(fit.fall.2018.additive.site.lrt,p.value=0.01)) 
# p 0.01 = 2004 Down, 3979 Up 15659 NotSig
summary(decideTestsDGE(fit.fall.2018.additive.grow.nogrow.lrt,p.value=0.01)) 
# p 0.01 = 1989 Down, 1730 Up 17923 NotSig

# multiplicative
summary(decideTestsDGE(fit.spring.2017.multiply.lrt,p.value=0.01)) 
# p 0.01 = 75 Down, 431 Up 21136 NotSig
summary(decideTestsDGE(fit.summer.2017.multiply.lrt,p.value=0.01)) 
# p 0.01 = 75 Down, 431 Up 21136 NotSig
summary(decideTestsDGE(fit.fall.2017.multiply.lrt,p.value=0.01)) 
# p 0.01 = 75 Down, 431 Up 21136 NotSig

summary(decideTestsDGE(fit.spring.2018.multiply.lrt,p.value=0.01)) 
# p 0.01 = 75 Down, 431 Up 21136 NotSig
summary(decideTestsDGE(fit.summer.2018.multiply.lrt,p.value=0.01)) 
# p 0.01 = 75 Down, 431 Up 21136 NotSig
summary(decideTestsDGE(fit.fall.2018.multiply.lrt,p.value=0.01)) 
# p 0.01 = 75 Down, 431 Up 21136 NotSig


#### Extract genes with a FDR < 0.01 ####
DEgene.spring.2017 <- topTags(fit.spring.2017.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2017,"./Data/DE.data/grow.nogrow/DEgene.spring.2017.csv")
DEgene.summer.2017 <- topTags(fit.summer.2017.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2017,"./Data/DE.data/grow.nogrow/DEgene.summer.2017.csv")
DEgene.fall.2017 <- topTags(fit.fall.2017.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2017,"./Data/DE.data/grow.nogrow/DEgene.fall.2017.csv")

DEgene.spring.2018 <- topTags(fit.spring.2018.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2018,"./Data/DE.data/grow.nogrow/DEgene.spring.2018.csv")
DEgene.summer.2018 <- topTags(fit.summer.2018.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2018,"./Data/DE.data/grow.nogrow/DEgene.summer.2018.csv")
DEgene.fall.2018 <- topTags(fit.fall.2018.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2018,"./Data/DE.data/grow.nogrow/DEgene.fall.2018.csv")

DEgene.spring.2019 <- topTags(fit.spring.2019.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2019,"./Data/DE.data/grow.nogrow/DEgene.spring.2019.csv")
DEgene.summer.2019 <- topTags(fit.summer.2019.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2019,"./Data/DE.data/grow.nogrow/DEgene.summer.2019.csv")
DEgene.fall.2019 <- topTags(fit.fall.2019.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2019,"./Data/DE.data/grow.nogrow/DEgene.fall.2019.csv")

DEgene.spring.2020 <- topTags(fit.spring.2020.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2020,"./Data/DE.data/grow.nogrow/DEgene.spring.2020.csv")

DEgene.spring.2017.HF <- topTags(fit.spring.2017.HF.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2017.HF,"./Data/DE.data/grow.nogrow/DEgene.spring.2017.HF.csv")
DEgene.summer.2017.HF <- topTags(fit.summer.2017.HF.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2017.HF,"./Data/DE.data/grow.nogrow/DEgene.summer.2017.HF.csv")
DEgene.fall.2017.HF <- topTags(fit.fall.2017.HF.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2017.HF,"./Data/DE.data/grow.nogrow/DEgene.fall.2017.HF.csv")

DEgene.spring.2018.HF <- topTags(fit.spring.2018.HF.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2018.HF,"./Data/DE.data/grow.nogrow/DEgene.spring.2018.HF.csv")
DEgene.summer.2018.HF <- topTags(fit.summer.2018.HF.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2018.HF,"./Data/DE.data/grow.nogrow/DEgene.summer.2018.HF.csv")
DEgene.fall.2018.HF <- topTags(fit.fall.2018.HF.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2018.HF,"./Data/DE.data/grow.nogrow/DEgene.fall.2018.HF.csv")

DEgene.spring.2017.SERC <- topTags(fit.spring.2017.SERC.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2017.SERC,"./Data/DE.data/grow.nogrow/DEgene.spring.2017.SERC.csv")
DEgene.summer.2017.SERC <- topTags(fit.summer.2017.SERC.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2017.SERC,"./Data/DE.data/grow.nogrow/DEgene.summer.2017.SERC.csv")
DEgene.fall.2017.SERC <- topTags(fit.fall.2017.SERC.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2017.SERC,"./Data/DE.data/grow.nogrow/DEgene.fall.2017.SERC.csv")

DEgene.spring.2018.SERC <- topTags(fit.spring.2018.SERC.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2018.SERC,"./Data/DE.data/grow.nogrow/DEgene.spring.2018.SERC.csv")
DEgene.summer.2018.SERC <- topTags(fit.summer.2018.SERC.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2018.SERC,"./Data/DE.data/grow.nogrow/DEgene.summer.2018.SERC.csv")
DEgene.fall.2018.SERC <- topTags(fit.fall.2018.SERC.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2018.SERC,"./Data/DE.data/grow.nogrow/DEgene.fall.2018.SERC.csv")

# additive model
DEgene.spring.2017.additive.site <- topTags(fit.spring.2017.additive.site.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2017.additive.site,"./Data/DE.data/grow.nogrow/DEgene.spring.2017.additive.site.csv")
DEgene.spring.2017.additive.grow.nogrow <- topTags(fit.spring.2017.additive.grow.nogrow.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2017.additive.grow.nogrow,"./Data/DE.data/grow.nogrow/DEgene.spring.2017.additive.grow.nogrow")
DEgene.summer.2017.additive.site <- topTags(fit.summer.2017.additive.site.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2017.additive.site,"./Data/DE.data/grow.nogrow/DEgene.summer.2017.additive.site.csv")
DEgene.summer.2017.additive.grow.nogrow <- topTags(fit.summer.2017.additive.grow.nogrow.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2017.additive.grow.nogrow,"./Data/DE.data/grow.nogrow/DEgene.summer.2017.additive.grow.nogrow")
DEgene.fall.2017.additive.site <- topTags(fit.fall.2017.additive.site.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2017.additive.site,"./Data/DE.data/grow.nogrow/DEgene.fall.2017.additive.site.csv")
DEgene.fall.2017.additive.grow.nogrow <- topTags(fit.fall.2017.additive.grow.nogrow.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2017.additive.grow.nogrow,"./Data/DE.data/grow.nogrow/DEgene.fall.2017.additive.grow.nogrow")

DEgene.spring.2018.additive.site <- topTags(fit.spring.2018.additive.site.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2018.additive.site,"./Data/DE.data/grow.nogrow/DEgene.spring.2018.additive.site.csv")
DEgene.spring.2018.additive.grow.nogrow <- topTags(fit.spring.2018.additive.grow.nogrow.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2018.additive.grow.nogrow,"./Data/DE.data/grow.nogrow/DEgene.spring.2018.additive.grow.nogrow.csv")
DEgene.summer.2018.additive.site <- topTags(fit.summer.2018.additive.site.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2018.additive.site,"./Data/DE.data/grow.nogrow/DEgene.summer.2018.additive.site.csv")
DEgene.summer.2018.additive.grow.nogrow <- topTags(fit.summer.2018.additive.grow.nogrow.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2018.additive.grow.nogrow,"./Data/DE.data/grow.nogrow/DEgene.summer.2018.additive.grow.nogrow.csv")
DEgene.fall.2018.additive.site <- topTags(fit.fall.2018.additive.site.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2018.additive.site,"./Data/DE.data/grow.nogrow/DEgene.fall.2018.additive.site.csv")
DEgene.fall.2018.additive.grow.nogrow <- topTags(fit.fall.2018.additive.grow.nogrow.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2018.additive.grow.nogrow,"./Data/DE.data/grow.nogrow/DEgene.fall.2018.additive.grow.nogrow.csv")

# multiplicative
DEgene.spring.2017.multiply <- topTags(fit.spring.2017.multiply.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2017.multiply,"./Data/DE.data/grow.nogrow/DEgene.spring.2017.multiply.csv")
DEgene.summer.2017.multiply <- topTags(fit.summer.2017.multiply.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2017.multiply,"./Data/DE.data/grow.nogrow/DEgene.summer.2017.multiply.csv")
DEgene.fall.2017.multiply <- topTags(fit.fall.2017.multiply.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2017.multiply,"./Data/DE.data/grow.nogrow/DEgene.fall.2017.multiply.csv")

DEgene.spring.2018.multiply <- topTags(fit.spring.2018.multiply.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.spring.2018.multiply,"./Data/DE.data/grow.nogrow/DEgene.spring.2018.multiply.csv")
DEgene.summer.2018.multiply <- topTags(fit.summer.2018.multiply.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.summer.2018.multiply,"./Data/DE.data/grow.nogrow/DEgene.summer.2018.multiply.csv")
DEgene.fall.2018.multiply <- topTags(fit.fall.2018.multiply.lrt,n = Inf,p.value = 0.01)$table
#write.csv(DEgene.fall.2018.multiply,"./Data/DE.data/grow.nogrow/DEgene.fall.2018.multiply.csv")

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

plotDE(rownames(DEgene.spring.2017)[1:9],DGE.data.spring.2017,Sample.Description.spring.2017)
plotDE(rownames(DEgene.summer.2017)[1:9],DGE.data.summer.2017,Sample.Description.summer.2017)
plotDE(rownames(DEgene.fall.2017)[1:9],DGE.data.fall.2017,Sample.Description.fall.2017)
plotDE(rownames(DEgene.spring.2018)[1:9],DGE.data.spring.2018,Sample.Description.spring.2018)
plotDE(rownames(DEgene.summer.2018)[1:9],DGE.data.summer.2018,Sample.Description.summer.2018)
plotDE(rownames(DEgene.fall.2018)[1:9],DGE.data.fall.2018,Sample.Description.fall.2018)
plotDE(rownames(DEgene.spring.2019)[1:9],DGE.data.spring.2019,Sample.Description.spring.2019)
plotDE(rownames(DEgene.summer.2019)[1:9],DGE.data.summer.2019,Sample.Description.summer.2019)
plotDE(rownames(DEgene.fall.2019)[1:9],DGE.data.fall.2019,Sample.Description.fall.2019)
plotDE(rownames(DEgene.spring.2020)[1:9],DGE.data.spring.2020,Sample.Description.spring.2020)

plotDE(rownames(DEgene.spring.2017.HF)[1:9],DGE.data.spring.2017.HF,Sample.Description.spring.2017.HF)
plotDE(rownames(DEgene.summer.2017.HF)[1:9],DGE.data.summer.2017.HF,Sample.Description.summer.2017.HF)
plotDE(rownames(DEgene.fall.2017.HF)[1:9],DGE.data.fall.2017.HF,Sample.Description.fall.2017.HF)
plotDE(rownames(DEgene.spring.2018.HF)[1:9],DGE.data.spring.2018.HF,Sample.Description.spring.2018.HF)
plotDE(rownames(DEgene.summer.2018.HF)[1:9],DGE.data.summer.2018.HF,Sample.Description.summer.2018.HF)
plotDE(rownames(DEgene.fall.2018.HF)[1:9],DGE.data.fall.2018.HF,Sample.Description.fall.2018.HF)

plotDE(rownames(DEgene.spring.2017.SERC)[1:9],DGE.data.spring.2017.SERC,Sample.Description.spring.2017.SERC)
plotDE(rownames(DEgene.summer.2017.SERC)[1:9],DGE.data.summer.2017.SERC,Sample.Description.summer.2017.SERC)
plotDE(rownames(DEgene.fall.2017.SERC)[1:9],DGE.data.fall.2017.SERC,Sample.Description.fall.2017.SERC)
plotDE(rownames(DEgene.spring.2018.SERC)[1:9],DGE.data.spring.2018.SERC,Sample.Description.spring.2018.SERC)
plotDE(rownames(DEgene.summer.2018.SERC)[1:9],DGE.data.summer.2018.SERC,Sample.Description.summer.2018.SERC)
plotDE(rownames(DEgene.fall.2018.SERC)[1:9],DGE.data.fall.2018.SERC,Sample.Description.fall.2018.SERC)

# additive models
plotDE_fill(rownames(DEgene.spring.2018.additive.site)[1:9],DGE.data.spring.2017.additive,Sample.Description.spring.2017)
# change Year to site for X axis in plotting function
plotDE_fill(rownames(DEgene.2017.18.additive.Site)[1:9],DGE.data.2017.18.site,Sample_Description.2017.2018)
# multiplicative
plotDE_fill(rownames(DEgene.2017.18.multiply.site.year)[1:9],DGE.data.2017.18.site,Sample_Description.2017.2018)
