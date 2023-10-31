#### Co-Expresion Analysis using just Significantly Differentiallly Expressed Genes ####

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

# only doing time steps where we have growth in this analysis

#### Subtract TMM values among years in data sets #####

TMM_2019.20 = read.csv("./Data/DE.data/TMM_2019.20_DEG.csv")
TMM_2019.20.sub = TMM_2019.20[,c(14:25)] - TMM_2019.20[,c(2:13)]
TMM_2019.20.sub = cbind(TMM_2019.20$Gene_ID,TMM_2019.20.sub)
colnames(TMM_2019.20.sub)[1] = "Gene_ID"

TMM_2017.18.Spring = read.csv("./Data/DE.data/TMM_2017.18.Spring_DEG.csv")
TMM_2017.18.Spring.sub = TMM_2017.18.Spring[,c(14:25,44:61)] - TMM_2017.18.Spring[,c(2:13,26:43)]
TMM_2017.18.Spring.sub = cbind(TMM_2017.18.Spring$Gene_ID,TMM_2017.18.Spring.sub)
colnames(TMM_2017.18.Spring.sub)[1] = "Gene_ID"

TMM_2017.18.Summer = read.csv("./Data/DE.data/TMM_2017.18.Summer_DEG.csv")
TMM_2017.18.Summer.sub = TMM_2017.18.Summer[,c(18:33,51:67)] - TMM_2017.18.Summer[,c(2:17,34:50)]
TMM_2017.18.Summer.sub = cbind(TMM_2017.18.Summer$Gene_ID,TMM_2017.18.Summer.sub)
colnames(TMM_2017.18.Summer.sub)[1] = "Gene_ID"

TMM_2017.18.Fall= read.csv("./Data/DE.data/TMM_2017.18.Fall_DEG.csv")
TMM_2017.18.Fall.sub = TMM_2017.18.Fall[,c(18:33,50:65)] - TMM_2017.18.Fall[,c(2:17,34:49)]
TMM_2017.18.Fall.sub = cbind(TMM_2017.18.Fall$Gene_ID,TMM_2017.18.Fall.sub)
colnames(TMM_2017.18.Fall.sub)[1] = "Gene_ID"

#### Calculate the Variance of each Gene across all of our samples ####

# Calculating Coefficient of variation function
calc.cv <- function(x, na.rm=TRUE) { 
  if(na.rm==TRUE) x <- na.omit(x)
  result <- sd(x) / mean(x) 
  result <- abs(result) 
  return(result)
}

TMM_2019.20_cv <- TMM_2019.20.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2017.18.Spring_cv <- TMM_2017.18.Spring.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2017.18.Summer_cv <- TMM_2017.18.Summer.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM_2017.18.Fall_cv <- TMM_2017.18.Fall.sub %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())

# Filter the data frame to contain the top 30% most variable genes)
# ties are kept together
TMM_2019.20_30  <- TMM_2019.20_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2017.18.Spring_30  <- TMM_2017.18.Spring_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2017.18.Summer_30  <- TMM_2017.18.Summer_cv  %>% slice_max(order_by = cv, prop = .30)
TMM_2017.18.Fall_30  <- TMM_2017.18.Fall_cv  %>% slice_max(order_by = cv, prop = .30)

# Deselect the cv column and flip our data frame to contain sample along the rows and genes along the columns
TMM_2019.20_30 <- select(TMM_2019.20_30, -cv)
TMM_2019.20_30 <- column_to_rownames(TMM_2019.20_30,var = "Gene_ID")
TMM_2019.20_30 <- as.data.frame(t(TMM_2019.20_30))

TMM_2017.18.Spring_30 <- select(TMM_2017.18.Spring_30, -cv)
TMM_2017.18.Spring_30 <- column_to_rownames(TMM_2017.18.Spring_30,var = "Gene_ID")
TMM_2017.18.Spring_30 <- as.data.frame(t(TMM_2017.18.Spring_30))

TMM_2017.18.Summer_30 <- select(TMM_2017.18.Summer_30, -cv)
TMM_2017.18.Summer_30 <- column_to_rownames(TMM_2017.18.Summer_30,var = "Gene_ID")
TMM_2017.18.Summer_30 <- as.data.frame(t(TMM_2017.18.Summer_30))

TMM_2017.18.Fall_30 <- select(TMM_2017.18.Fall_30, -cv)
TMM_2017.18.Fall_30 <- column_to_rownames(TMM_2017.18.Fall_30,var = "Gene_ID")
TMM_2017.18.Fall_30 <- as.data.frame(t(TMM_2017.18.Fall_30))

# Continue with formatting
Sample_Description = read_csv("./Formatted.Data/FAGR.description.csv")

#### QC ####
gsg_2019.20 = goodSamplesGenes(TMM_2019.20_30, verbose = 3);
gsg_2019.20$allOK
gsg_2017.18.Spring = goodSamplesGenes(TMM_2017.18.Spring_30, verbose = 3);
gsg_2017.18.Spring$allOK
gsg_2017.18.Summer = goodSamplesGenes(TMM_2017.18.Summer_30, verbose = 3);
gsg_2017.18.Summer$allOK
gsg_2017.18.Fall = goodSamplesGenes(TMM_2017.18.Fall_30, verbose = 3);
gsg_2017.18.Fall$allOK

#### Soft Threshold ####
# Choose a set of soft-thresholding powers
# value is the power to which co-expression similarity is raised to calculate adjacency
# chose the value closest to 0.80
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft_2019.20 <- pickSoftThreshold(TMM_2019.20_30, powerVector = powers, verbose = 5) # 10 at 0.80
sft_2017.18.Spring <- pickSoftThreshold(TMM_2017.18.Spring_30, powerVector = powers, verbose = 5) # 5 at 0.90
sft_2017.18.Summer <- pickSoftThreshold(TMM_2017.18.Summer_30, powerVector = powers, verbose = 5) # 5 at 0.892
sft_2017.18.Fall <- pickSoftThreshold(TMM_2017.18.Fall_30, powerVector = powers, verbose = 5) # 8 at 0.90

#### Adjacency ####
# signed hybrid sets all negatively correlated genes to 0
# choice of power comes from soft threshold analysis above
# power value was chosen as value closest to 90%
adjacency_2019.20 = adjacency(TMM_2019.20_30, power = 10,type = "signed hybrid")
adjacency_2017.18.Spring = adjacency(TMM_2017.18.Spring_30, power = 6,type = "signed hybrid")
adjacency_2017.18.Summer = adjacency(TMM_2017.18.Summer_30, power = 6,type = "signed hybrid")
adjacency_2017.18.Fall = adjacency(TMM_2017.18.Fall_30, power = 7,type = "signed hybrid")

#### Topological overlap matrix (TOM) ####
# to minimize effects of noise and spurious associations, tranform adjacency into TOM, calculate dissimilarity
# need to rm adjacency matrices and TOM after use to make space
TOM_2019.20 = TOMsimilarity(adjacency_2019.20, TOMType = "signed Nowick")
dissTOM_2019.20 = 1 - TOM_2019.20
rm(TOM_2019.20)
rm(adjacency_2019.20)

TOM_2017.18.Spring = TOMsimilarity(adjacency_2017.18.Spring, TOMType = "signed Nowick")
dissTOM_2017.18.Spring = 1 - TOM_2017.18.Spring
rm(TOM_2017.18.Spring)
rm(adjacency_2017.18.Spring)

TOM_2017.18.Summer = TOMsimilarity(adjacency_2017.18.Summer, TOMType = "signed Nowick")
dissTOM_2017.18.Summer = 1 - TOM_2017.18.Summer
rm(TOM_2017.18.Summer)
rm(adjacency_2017.18.Summer)

TOM_2017.18.Fall = TOMsimilarity(adjacency_2017.18.Fall, TOMType = "signed Nowick")
dissTOM_2017.18.Fall = 1 - TOM_2017.18.Fall
rm(TOM_2017.18.Fall)
rm(adjacency_2017.18.Fall)

#### Gene Clustering Plot ####
# Code from this section until the end was repeated for each dissTom element.
# For each repeat of dataset change: name of dissTOM, soft power number, TPM_2017_Fall_30, ect.
# Call the hierarchical clustering function
# hclust dendrogram of genes
geneTree = hclust(as.dist(dissTOM_2019.20), method = "average")

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12, 9)
jpeg("./Plots/DE/dissTOM_2019.20.sigDEG-tree.jpeg",height = 1080, width = 1920)
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
  distM = dissTOM_2019.20,
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
jpeg("./Plots/DE/dissTOM_2019.20.sigDEG-tree-blocks.jpeg",height = 1080, width = 1920)
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
MEList <- moduleEigengenes(TMM_2019.20_30,
                           colors = dynamicColors,
                           softPower = 10)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(7, 6)
jpeg("./Plots/DE/ME-tree_2019.20.sigDEG.jpeg",height = 1080, width = 1920)

plot(METree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "")

MEDissThres <-  0.25 # threshold for merging modules, corresponds to module correlation of 0.75

# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")

# Call an automatic merging function
merge <- mergeCloseModules(TMM_2019.20_30,
                           dynamicColors,
                           cutHeight = MEDissThres,
                           verbose = 3)
# The merged module colors
mergedColors <- merge$colors

# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs
dev.off()

sizeGrWindow(12, 9)
pdf(file = "./Plots/DE/geneDendro_2019.20.sigDEG.pdf", wi = 9, he = 6)
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
Module_Heatmap <- TMM_2019.20_30 %>%
  rownames_to_column(var = "sample") %>%
  column_to_rownames(var = "sample")


ModuleGenes_Df <- data.frame(gene = NULL, color = NULL)
for(color in unique(moduleColors)) {
  temp <- data.frame(names(TMM_2019.20_30)[moduleColors == color])
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
jpeg("./Plots/DE/ModuleGene_Heatmap_2019.20.sigDEG.jpeg",
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

#write_csv(growth_MEs,"./Data/DE.MEs/time.2019.20.sigDEG.MEs.csv")









#### Models #####
# had to manual take the difference in RGR and growth between years
# added this to ME dataframes
# Can only do single season comparisons because growth would be repeated between seasons

##### Spring 2017-2018 #####
Spring.2017.18.MEs = read.csv("./Data/DE.MEs/time.2017.18.Spring.sigDEG.MEs.csv") %>%
  mutate(growth.change.log = log(growth.change+1),
         RGR.change.log = log(RGR.change+1)) %>%
  rename_all(~ str_replace(., "ME","ME_"))

Spring.2017.18.MEs.2 = Spring.2017.18.MEs %>%
  pivot_longer(starts_with("ME"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,eigen_value,RGR.change,growth.change,RGR.change.log,growth.change.log))

Spring.2017.18.MEs.RGR.change <- Spring.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR.change ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.18.MEs.RGR.change <- Spring.2017.18.MEs.RGR.change %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR.change)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Spring.2017.18.MEs.RGR.change %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Spring.2017.18.MEs.growth.change <- Spring.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth.change ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.18.MEs.growth.change <- Spring.2017.18.MEs.growth.change %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth.change)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Spring.2017.18.MEs.growth.change %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Spring.2017.18.MEs.RGR.change.log <- Spring.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR.change.log ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.18.MEs.RGR.change.log <- Spring.2017.18.MEs.RGR.change.log %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR.change.log)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Spring.2017.18.MEs.RGR.change.log %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Spring.2017.18.MEs.growth.change.log <- Spring.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth.change.log ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.18.MEs.growth.change.log <- Spring.2017.18.MEs.growth.change.log %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth.change.log)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Spring.2017.18.MEs.growth.change.log %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

#### Summer 2017-2018 ####
Summer.2017.18.MEs = read.csv("./Data/DE.MEs/time.2017.18.Summer.sigDEG.MEs.csv") %>%
  mutate(growth.change.log = log(growth.change+1),
         RGR.change.log = log(RGR.change+1)) %>%
  rename_all(~ str_replace(., "ME","ME_"))

Summer.2017.18.MEs.2 = Summer.2017.18.MEs %>%
  pivot_longer(starts_with("ME"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,eigen_value,RGR.change,growth.change,RGR.change.log,growth.change.log))

Summer.2017.18.MEs.RGR.change <- Summer.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR.change ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2017.18.MEs.RGR.change <- Summer.2017.18.MEs.RGR.change %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR.change)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Summer.2017.18.MEs.RGR.change %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Summer.2017.18.MEs.growth.change <- Summer.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth.change ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2017.18.MEs.growth.change <- Summer.2017.18.MEs.growth.change %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth.change)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Summer.2017.18.MEs.growth.change %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Summer.2017.18.MEs.RGR.change.log <- Summer.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR.change.log ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2017.18.MEs.RGR.change.log <- Summer.2017.18.MEs.RGR.change.log %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR.change.log)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Summer.2017.18.MEs.RGR.change.log %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Summer.2017.18.MEs.growth.change.log <- Summer.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth.change.log ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Summer.2017.18.MEs.growth.change.log <- Summer.2017.18.MEs.growth.change.log %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth.change.log)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Summer.2017.18.MEs.growth.change.log %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

#### Fall 2017-2018 ####
Fall.2017.18.MEs = read.csv("./Data/DE.MEs/time.2017.18.Fall.sigDEG.MEs.csv") %>%
  mutate(growth.change.log = log(growth.change+1),
         RGR.change.log = log(RGR.change+1)) %>%
  rename_all(~ str_replace(., "ME","ME_"))

Fall.2017.18.MEs.2 = Fall.2017.18.MEs %>%
  pivot_longer(starts_with("ME"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,eigen_value,RGR.change,growth.change,RGR.change.log,growth.change.log))

Fall.2017.18.MEs.RGR.change <- Fall.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR.change ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2017.18.MEs.RGR.change <- Fall.2017.18.MEs.RGR.change %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR.change)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Fall.2017.18.MEs.RGR.change %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Fall.2017.18.MEs.growth.change <- Fall.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth.change ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2017.18.MEs.growth.change <- Fall.2017.18.MEs.growth.change %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth.change)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Fall.2017.18.MEs.growth.change %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Fall.2017.18.MEs.RGR.change.log <- Fall.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR.change.log ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2017.18.MEs.RGR.change.log <- Fall.2017.18.MEs.RGR.change.log %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR.change.log)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Fall.2017.18.MEs.RGR.change.log %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

Fall.2017.18.MEs.growth.change.log <- Fall.2017.18.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth.change.log ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Fall.2017.18.MEs.growth.change.log <- Fall.2017.18.MEs.growth.change.log %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth.change.log)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

Fall.2017.18.MEs.growth.change.log %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

#### 2019-2020 ####
time.2019.20.MEs = read.csv("./Data/DE.MEs/time.2019.20.sigDEG.MEs.csv") %>%
  mutate(growth.change.log = log(growth.change+1),
         RGR.change.log = log(RGR.change+1)) %>%
  rename_all(~ str_replace(., "ME","ME_"))

time.2019.20.MEs.2 = time.2019.20.MEs %>%
  pivot_longer(starts_with("ME"),names_to = "module", values_to = "eigen_value") %>%
  nest(Data = c(sample,eigen_value,RGR.change,growth.change,RGR.change.log,growth.change.log))

time.2019.20.MEs.RGR.change <- time.2019.20.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR.change ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

time.2019.20.MEs.RGR.change <- time.2019.20.MEs.RGR.change %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR.change)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

time.2019.20.MEs.RGR.change %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

time.2019.20.MEs.growth.change <- time.2019.20.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth.change ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

time.2019.20.MEs.growth.change <- time.2019.20.MEs.growth.change %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth.change)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

time.2019.20.MEs.growth.change %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

time.2019.20.MEs.RGR.change.log <- time.2019.20.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(RGR.change.log ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

time.2019.20.MEs.RGR.change.log <- time.2019.20.MEs.RGR.change.log %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = RGR.change.log)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

time.2019.20.MEs.RGR.change.log %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)

time.2019.20.MEs.growth.change.log <- time.2019.20.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(growth.change.log ~ eigen_value, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

time.2019.20.MEs.growth.change.log <- time.2019.20.MEs.growth.change.log %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = eigen_value, y = growth.change.log)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)

time.2019.20.MEs.growth.change.log %>%
  filter(term == "eigen_value") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)
