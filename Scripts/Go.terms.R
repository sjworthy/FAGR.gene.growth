# matching genes in significant modules to GO terms

library(goseq)

#### read in gene/color modules ####

summer.2017 = read.csv("./Data/module.genes/summer.2017.module.genes.csv", row.names = 1)
spring.2019 = read.csv("./Data/module.genes/summer.2019.module.genes.csv", row.names = 1)
fall.2019 = read.csv("./Data/module.genes/fall.2019.module.genes.csv", row.names = 1)
spring.2017.HF = read.csv("./Data/module.genes/spring.2017.HF.module.genes.csv", row.names = 1)
spring.2017.SERC = read.csv("./Data/module.genes/spring.2017.SERC.module.genes.csv", row.names = 1)

#### filter genes that belong to module of interest ####

summer.2017.turq = summer.2017 %>%
  filter(color == "turquoise")
# 616 genes
summer.2017.ltyell = summer.2017 %>%
  filter(color == "lightyellow")
# 132 genes
spring.2019.grey = spring.2019 %>%
  filter(color == "grey60")
#169 genes
spring.2019.dkred = spring.2019 %>%
  filter(color == "darkred")
# 149 genes
fall.2019.ryblue = fall.2019 %>%
  filter(color == "royalblue") 
# 140 genes
spring.2017.HF.tan = spring.2017.HF %>%
  filter(color == "tan")
# 181 genes
spring.2017.SERC.ltblue = spring.2017.SERC %>%
  filter(color == "lightsteelblue1")
# 76 genes


#### Import Go terms ####
# had to remove the B:,C:,M: before GO terms
Goterms = read_tsv("./Raw.Data/Fagr_gfacs_GOs_cats.txt")
colnames(Goterms)[1]="gene"

#### Format data for GOseq ####

# need to reduce the gene.length data to only contain entries for those genes in our expressed.genes set. 
# also need this as a vector

# Calculating Coefficient of variation function
calc.cv <- function(x, na.rm=TRUE) { 
  if(na.rm==TRUE) x <- na.omit(x)
  result <- sd(x) / mean(x) 
  result <- abs(result) 
  return(result)
}

TMM.summer.2017 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2017.csv")
TMM.summer.2017.cv <- TMM.summer.2017 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.summer.2017.30  <- TMM.summer.2017.cv  %>% slice_max(order_by = cv, prop = .30)

summer.2017.gene.lengths.vector <- Goterms$Length[Goterms$gene %in% TMM.summer.2017.30$Gene_ID]
names(summer.2017.gene.lengths.vector) <- Goterms$gene[Goterms$gene %in% TMM.summer.2017.30$Gene_ID]
head(summer.2017.gene.lengths.vector)

# format GO terms for GOseq want them in list format, and need to separate the terms into separate elements.

# split GO terms into separate columns

go.list <- strsplit(Goterms$`GO IDs`,split=";")
names(go.list) <- Goterms$gene
head(go.list)

# for each gene in expressed gene, return FALSE if it is not in our module and TRUE if it is.
modgenes <- TMM.summer.2017.30$Gene_ID %in% summer.2017.turq$gene
table(modgenes)
names(modgenes) <- TMM.summer.2017.30$Gene_ID
head(modgenes)

modgenes <- as.numeric(modgenes) #convert to 0s and 1s
head(modgenes)
sum(modgenes) # number of genes in module

#### Calculate over-representation ####
#determines if there is bias due to gene length.  The plot shows the relationship.
nullp.result <- nullp(DEgenes = modgenes, bias.data = summer.2017.gene.lengths.vector)

#calculate p-values for each GO term
rownames(nullp.result) <- names(summer.2017.gene.lengths.vector) #because of a bug in nullp()
GO.out <- goseq(pwf = nullp.result, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

#list over-represented GO terms (p < 0.05)
GO.out[GO.out$over_represented_pvalue < 0.05,]

# fdr correction

GO.out.2 = GO.out %>%
  mutate(fdr.over.pvalue = p.adjust(over_represented_pvalue, method = "fdr"))
# all 1?

#### Join module genes and Go Terms ####

summer.2017.turq.GO = left_join(summer.2017.turq, Goterms, by = "gene")
summer.2017.ltyell.GO = left_join(summer.2017.ltyell, Goterms, by = "gene")
spring.2019.grey.GO = left_join(spring.2019.grey, Goterms, by = "gene")
spring.2019.dkred.GO = left_join(spring.2019.dkred, Goterms, by = "gene")
fall.2019.ryblue.GO = left_join(fall.2019.ryblue, Goterms, by = "gene")
spring.2017.HF.tan.GO = left_join(spring.2017.HF.tan, Goterms, by = "gene")
spring.2017.SERC.ltblue.GO = left_join(spring.2017.SERC.ltblue, Goterms, by = "gene")




#### Identify hub genes ####
# identify highly connected intramodular hub genes
# calculate correlation of module eigenenes with gene expression profile

# Calculating Coefficient of variation function
calc.cv <- function(x, na.rm=TRUE) { 
  if(na.rm==TRUE) x <- na.omit(x)
  result <- sd(x) / mean(x) 
  result <- abs(result) 
  return(result)
}

##### summer 2017 #####
summer.2017.MEs = read.csv("./Data/grow.nogrow.MEs/summer.2017.MEs.csv")
module_eigengenes = summer.2017.MEs[,c(2:32)]
row.names(module_eigengenes) = summer.2017.MEs$sample

TMM.summer.2017 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.summer.2017.csv")
TMM.summer.2017.cv <- TMM.summer.2017 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.summer.2017.30  <- TMM.summer.2017.cv  %>% slice_max(order_by = cv, prop = .30)
TMM.summer.2017.30 <- select(TMM.summer.2017.30, -cv)
TMM.summer.2017.30 <- column_to_rownames(TMM.summer.2017.30,var = "Gene_ID")
TMM.summer.2017.30 <- as.data.frame(t(TMM.summer.2017.30))

nSamples = nrow(TMM.summer.2017.30)

module.membership.measure = cor(module_eigengenes, TMM.summer.2017.30, use = "p")
module.membership.measure.pvals = as.data.frame(corPvalueStudent(module.membership.measure, nSamples))

# subset for significant modules

turq.mod = t(subset(module.membership.measure.pvals,rownames(module.membership.measure.pvals)=="ME_turquoise"))
turq.mod.2 = as.data.frame(turq.mod) %>%
  mutate(fdr.p = p.adjust(ME_turquoise)) %>%
  filter(fdr.p < 0.01)

# 129 of 616 genes considered significant hub genes
# most correlated gene is FAGR_UCONN_27225







expressed_genes_match <-
  datExpr[,colnames(datExpr) %in% names(Gene.length.vector)]
#Do the reverse to make sure everything matches up(reduce the  expressed.genes set  to only contain entries for those genes in our gene.length.vector)


GO.list <- strsplit(Goterms$GOs, split = ",")
names(GO.list) <- Goterms$Gene_ID
#We need our Goterms in a list format and in separate elements







