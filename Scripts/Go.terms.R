# matching genes in significant modules to GO terms

library(tidyverse)
library(goseq)
# https://bioconductor.org/packages/devel/bioc/vignettes/goseq/inst/doc/goseq.pdf

#### read in gene/color modules ####

summer.2017.module.genes = read.csv("./Data/module.genes/summer.2017.module.genes.csv", row.names = 1)
fall.2017.module.genes = read.csv("./Data/module.genes/fall.2017.module.genes.csv", row.names = 1)
spring.2019.module.genes = read.csv("./Data/module.genes/spring.2019.module.genes.csv", row.names = 1)
fall.2019.module.genes = read.csv("./Data/module.genes/fall.2019.module.genes.csv", row.names = 1)
spring.2017.HF.module.genes = read.csv("./Data/module.genes/spring.2017.HF.module.genes.csv", row.names = 1)
fall.2017.HF.module.genes = read.csv("./Data/module.genes/fall.2017.HF.module.genes.csv", row.names = 1)
spring.2017.SERC.module.genes = read.csv("./Data/module.genes/spring.2017.SERC.module.genes.csv", row.names = 1)

#### filter genes that belong to module of interest ####

summer.2017.turq = summer.2017.module.genes %>%
  filter(color == "turquoise")
# 616 genes
summer.2017.ltyell = summer.2017.module.genes %>%
  filter(color == "lightyellow")
# 132 genes
summer.2017.dkgrey = summer.2017.module.genes %>%
  filter(color == "darkgrey")
# 83 genes
summer.2017.green = summer.2017.module.genes %>%
  filter(color == "green")
# 429 genes
fall.2017.cyan = fall.2017.module.genes %>%
  filter(color == "cyan")
# 129 genes
fall.2017.mdntblue = fall.2017.module.genes %>%
  filter(color == "midnightblue")
# 128 genes
spring.2019.grey = spring.2019.module.genes %>%
  filter(color == "grey60")
#139 genes
spring.2019.dkred = spring.2019.module.genes %>%
  filter(color == "darkred")
# 105 genes
fall.2019.ryblue = fall.2019.module.genes %>%
  filter(color == "royalblue") 
# 140 genes
spring.2017.HF.tan = spring.2017.HF.module.genes %>%
  filter(color == "tan")
# 181 genes
fall.2017.HF.green = fall.2017.HF.module.genes %>%
  filter(color == "green")
# 239 genes
fall.2017.HF.salmon = fall.2017.HF.module.genes %>%
  filter(color == "salmon")
# 180 genes
fall.2017.HF.purple = fall.2017.HF.module.genes %>%
  filter(color == "purple")
# 206 genes
fall.2017.HF.black = fall.2017.HF.module.genes %>%
  filter(color == "black")
# 233 genes
fall.2017.HF.tan = fall.2017.HF.module.genes %>%
  filter(color == "tan")
# 195 genes
spring.2017.SERC.ltblue = spring.2017.SERC.module.genes %>%
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

TMM.fall.2017 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2017.csv")
TMM.fall.2017.cv <- TMM.fall.2017 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.fall.2017.30  <- TMM.fall.2017.cv  %>% slice_max(order_by = cv, prop = .30)

fall.2017.gene.lengths.vector <- Goterms$Length[Goterms$gene %in% TMM.fall.2017.30$Gene_ID]
names(fall.2017.gene.lengths.vector) <- Goterms$gene[Goterms$gene %in% TMM.fall.2017.30$Gene_ID]
head(fall.2017.gene.lengths.vector)

TMM.spring.2019 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2019.csv")
TMM.spring.2019.cv <- TMM.spring.2019 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.spring.2019.30  <- TMM.spring.2019.cv  %>% slice_max(order_by = cv, prop = .30)

spring.2019.gene.lengths.vector <- Goterms$Length[Goterms$gene %in% TMM.spring.2019.30$Gene_ID]
names(spring.2019.gene.lengths.vector) <- Goterms$gene[Goterms$gene %in% TMM.spring.2019.30$Gene_ID]
head(spring.2019.gene.lengths.vector)

TMM.fall.2019 = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2019.csv")
TMM.fall.2019.cv <- TMM.fall.2019 %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.fall.2019.30  <- TMM.fall.2019.cv  %>% slice_max(order_by = cv, prop = .30)

fall.2019.gene.lengths.vector <- Goterms$Length[Goterms$gene %in% TMM.fall.2019.30$Gene_ID]
names(fall.2019.gene.lengths.vector) <- Goterms$gene[Goterms$gene %in% TMM.fall.2019.30$Gene_ID]
head(fall.2019.gene.lengths.vector)

TMM.spring.2017.HF = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2017.HF.csv")
TMM.spring.2017.HF.cv <- TMM.spring.2017.HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.spring.2017.HF.30  <- TMM.spring.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)

spring.2017.HF.gene.lengths.vector <- Goterms$Length[Goterms$gene %in% TMM.spring.2017.HF.30$Gene_ID]
names(spring.2017.HF.gene.lengths.vector) <- Goterms$gene[Goterms$gene %in% TMM.spring.2017.HF.30$Gene_ID]
head(spring.2017.HF.gene.lengths.vector)

TMM.fall.2017.HF = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.fall.2017.HF.csv")
TMM.fall.2017.HF.cv <- TMM.fall.2017.HF %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.fall.2017.HF.30  <- TMM.fall.2017.HF.cv  %>% slice_max(order_by = cv, prop = .30)

fall.2017.HF.gene.lengths.vector <- Goterms$Length[Goterms$gene %in% TMM.fall.2017.HF.30$Gene_ID]
names(fall.2017.HF.gene.lengths.vector) <- Goterms$gene[Goterms$gene %in% TMM.fall.2017.HF.30$Gene_ID]
head(fall.2017.HF.gene.lengths.vector)

TMM.spring.2017.SERC = read.csv("./Data/grow.nogrow/TMM.NormData.LogCPM.spring.2017.SERC.csv")
TMM.spring.2017.SERC.cv <- TMM.spring.2017.SERC %>% rowwise() %>% mutate(cv = calc.cv(c_across(-Gene_ID))) %>% ungroup() %>% select(Gene_ID, cv, everything())
TMM.spring.2017.SERC.30  <- TMM.spring.2017.SERC.cv  %>% slice_max(order_by = cv, prop = .30)

spring.2017.SERC.gene.lengths.vector <- Goterms$Length[Goterms$gene %in% TMM.spring.2017.SERC.30$Gene_ID]
names(spring.2017.SERC.gene.lengths.vector) <- Goterms$gene[Goterms$gene %in% TMM.spring.2017.SERC.30$Gene_ID]
head(spring.2017.SERC.gene.lengths.vector)

# format GO terms for GOseq want them in list format, and need to separate the terms into separate elements.

# split GO terms into separate columns

go.list <- strsplit(Goterms$`GO IDs`,split=";")
names(go.list) <- Goterms$gene
head(go.list)

# for each gene in expressed gene, return FALSE if it is not in our module and TRUE if it is.
modgenes.summer.2017.turq <- TMM.summer.2017.30$Gene_ID %in% summer.2017.turq$gene
table(modgenes.summer.2017.turq)
names(modgenes.summer.2017.turq) <- TMM.summer.2017.30$Gene_ID
head(modgenes.summer.2017.turq)

modgenes.summer.2017.turq <- as.numeric(modgenes.summer.2017.turq) #convert to 0s and 1s
head(modgenes.summer.2017.turq)
sum(modgenes.summer.2017.turq) # number of genes in module

modgenes.summer.2017.ltyell <- TMM.summer.2017.30$Gene_ID %in% summer.2017.ltyell$gene
table(modgenes.summer.2017.ltyell)
names(modgenes.summer.2017.ltyell) <- TMM.summer.2017.30$Gene_ID
head(modgenes.summer.2017.ltyell)

modgenes.summer.2017.ltyell <- as.numeric(modgenes.summer.2017.ltyell) #convert to 0s and 1s
head(modgenes.summer.2017.ltyell)
sum(modgenes.summer.2017.ltyell) # number of genes in module

modgenes.summer.2017.dkgrey <- TMM.summer.2017.30$Gene_ID %in% summer.2017.dkgrey$gene
table(modgenes.summer.2017.dkgrey)
names(modgenes.summer.2017.dkgrey) <- TMM.summer.2017.30$Gene_ID
head(modgenes.summer.2017.dkgrey)

modgenes.summer.2017.dkgrey <- as.numeric(modgenes.summer.2017.dkgrey) #convert to 0s and 1s
head(modgenes.summer.2017.dkgrey)
sum(modgenes.summer.2017.dkgrey) # number of genes in module

modgenes.summer.2017.green <- TMM.summer.2017.30$Gene_ID %in% summer.2017.green$gene
table(modgenes.summer.2017.green)
names(modgenes.summer.2017.green) <- TMM.summer.2017.30$Gene_ID
head(modgenes.summer.2017.green)

modgenes.summer.2017.green <- as.numeric(modgenes.summer.2017.green) #convert to 0s and 1s
head(modgenes.summer.2017.green)
sum(modgenes.summer.2017.green) # number of genes in module

modgenes.fall.2017.cyan <- TMM.fall.2017.30$Gene_ID %in% fall.2017.cyan$gene
table(modgenes.fall.2017.cyan)
names(modgenes.fall.2017.cyan) <- TMM.fall.2017.30$Gene_ID
head(modgenes.fall.2017.cyan)

modgenes.fall.2017.cyan <- as.numeric(modgenes.fall.2017.cyan) #convert to 0s and 1s
head(modgenes.fall.2017.cyan)
sum(modgenes.fall.2017.cyan) # number of genes in module

modgenes.fall.2017.mdntblue <- TMM.fall.2017.30$Gene_ID %in% fall.2017.mdntblue$gene
table(modgenes.fall.2017.mdntblue)
names(modgenes.fall.2017.mdntblue) <- TMM.fall.2017.30$Gene_ID
head(modgenes.fall.2017.mdntblue)

modgenes.fall.2017.mdntblue <- as.numeric(modgenes.fall.2017.mdntblue) #convert to 0s and 1s
head(modgenes.fall.2017.mdntblue)
sum(modgenes.fall.2017.mdntblue) # number of genes in module

modgenes.spring.2019.grey <- TMM.spring.2019.30$Gene_ID %in% spring.2019.grey$gene
table(modgenes.spring.2019.grey)
names(modgenes.spring.2019.grey) <- TMM.spring.2019.30$Gene_ID
head(modgenes.spring.2019.grey)

modgenes.spring.2019.grey <- as.numeric(modgenes.spring.2019.grey) #convert to 0s and 1s
head(modgenes.spring.2019.grey)
sum(modgenes.spring.2019.grey) # number of genes in module

modgenes.spring.2019.dkred <- TMM.spring.2019.30$Gene_ID %in% spring.2019.dkred$gene
table(modgenes.spring.2019.dkred)
names(modgenes.spring.2019.dkred) <- TMM.spring.2019.30$Gene_ID
head(modgenes.spring.2019.dkred)

modgenes.spring.2019.dkred <- as.numeric(modgenes.spring.2019.dkred) #convert to 0s and 1s
head(modgenes.spring.2019.dkred)
sum(modgenes.spring.2019.dkred) # number of genes in module

modgenes.fall.2019.ryblue <- TMM.fall.2019.30$Gene_ID %in% fall.2019.ryblue$gene
table(modgenes.fall.2019.ryblue)
names(modgenes.fall.2019.ryblue) <- TMM.fall.2019.30$Gene_ID
head(modgenes.fall.2019.ryblue)

modgenes.fall.2019.ryblue <- as.numeric(modgenes.fall.2019.ryblue) #convert to 0s and 1s
head(modgenes.fall.2019.ryblue)
sum(modgenes.fall.2019.ryblue) # number of genes in module

modgenes.spring.2017.HF.tan <- TMM.spring.2017.HF.30$Gene_ID %in% spring.2017.HF.tan$gene
table(modgenes.spring.2017.HF.tan)
names(modgenes.spring.2017.HF.tan) <- TMM.spring.2017.HF.30$Gene_ID
head(modgenes.spring.2017.HF.tan)

modgenes.spring.2017.HF.tan <- as.numeric(modgenes.spring.2017.HF.tan) #convert to 0s and 1s
head(modgenes.spring.2017.HF.tan)
sum(modgenes.spring.2017.HF.tan) # number of genes in module

modgenes.fall.2017.HF.green <- TMM.fall.2017.HF.30$Gene_ID %in% fall.2017.HF.green$gene
table(modgenes.fall.2017.HF.green)
names(modgenes.fall.2017.HF.green) <- TMM.fall.2017.HF.30$Gene_ID
head(modgenes.fall.2017.HF.green)

modgenes.fall.2017.HF.green <- as.numeric(modgenes.fall.2017.HF.green) #convert to 0s and 1s
head(modgenes.fall.2017.HF.green)
sum(modgenes.fall.2017.HF.green) # number of genes in module

modgenes.fall.2017.HF.salmon <- TMM.fall.2017.HF.30$Gene_ID %in% fall.2017.HF.salmon$gene
table(modgenes.fall.2017.HF.salmon)
names(modgenes.fall.2017.HF.salmon) <- TMM.fall.2017.HF.30$Gene_ID
head(modgenes.fall.2017.HF.salmon)

modgenes.fall.2017.HF.salmon <- as.numeric(modgenes.fall.2017.HF.salmon) #convert to 0s and 1s
head(modgenes.fall.2017.HF.salmon)
sum(modgenes.fall.2017.HF.salmon) # number of genes in module

modgenes.fall.2017.HF.purple <- TMM.fall.2017.HF.30$Gene_ID %in% fall.2017.HF.purple$gene
table(modgenes.fall.2017.HF.purple)
names(modgenes.fall.2017.HF.purple) <- TMM.fall.2017.HF.30$Gene_ID
head(modgenes.fall.2017.HF.purple)

modgenes.fall.2017.HF.purple <- as.numeric(modgenes.fall.2017.HF.purple) #convert to 0s and 1s
head(modgenes.fall.2017.HF.purple)
sum(modgenes.fall.2017.HF.purple) # number of genes in module

modgenes.fall.2017.HF.black <- TMM.fall.2017.HF.30$Gene_ID %in% fall.2017.HF.black$gene
table(modgenes.fall.2017.HF.black)
names(modgenes.fall.2017.HF.black) <- TMM.fall.2017.HF.30$Gene_ID
head(modgenes.fall.2017.HF.black)

modgenes.fall.2017.HF.black <- as.numeric(modgenes.fall.2017.HF.black) #convert to 0s and 1s
head(modgenes.fall.2017.HF.black)
sum(modgenes.fall.2017.HF.black) # number of genes in module

modgenes.fall.2017.HF.tan <- TMM.fall.2017.HF.30$Gene_ID %in% fall.2017.HF.tan$gene
table(modgenes.fall.2017.HF.tan)
names(modgenes.fall.2017.HF.tan) <- TMM.fall.2017.HF.30$Gene_ID
head(modgenes.fall.2017.HF.tan)

modgenes.fall.2017.HF.tan <- as.numeric(modgenes.fall.2017.HF.tan) #convert to 0s and 1s
head(modgenes.fall.2017.HF.tan)
sum(modgenes.fall.2017.HF.tan) # number of genes in module

modgenes.spring.2017.SERC.ltblue <- TMM.spring.2017.SERC.30$Gene_ID %in% spring.2017.SERC.ltblue$gene
table(modgenes.spring.2017.SERC.ltblue)
names(modgenes.spring.2017.SERC.ltblue) <- TMM.spring.2017.SERC.30$Gene_ID
head(modgenes.spring.2017.SERC.ltblue)

modgenes.spring.2017.SERC.ltblue <- as.numeric(modgenes.spring.2017.SERC.ltblue) #convert to 0s and 1s
head(modgenes.spring.2017.SERC.ltblue)
sum(modgenes.spring.2017.SERC.ltblue) # number of genes in module

#### Calculate over-representation ####
#determines if there is bias due to gene length.  The plot shows the relationship.
nullp.result.summer.2017.turq <- nullp(DEgenes = modgenes.summer.2017.turq, bias.data = summer.2017.gene.lengths.vector)
nullp.result.summer.2017.ltyell <- nullp(DEgenes = modgenes.summer.2017.ltyell, bias.data = summer.2017.gene.lengths.vector)
nullp.result.summer.2017.dkgrey <- nullp(DEgenes = modgenes.summer.2017.dkgrey, bias.data = summer.2017.gene.lengths.vector)
nullp.result.summer.2017.green <- nullp(DEgenes = modgenes.summer.2017.green, bias.data = summer.2017.gene.lengths.vector)
nullp.result.fall.2017.cyan <- nullp(DEgenes = modgenes.fall.2017.cyan, bias.data = fall.2017.gene.lengths.vector)
nullp.result.fall.2017.mdntblue <- nullp(DEgenes = modgenes.fall.2017.mdntblue, bias.data = fall.2017.gene.lengths.vector)
nullp.result.spring.2019.grey <- nullp(DEgenes = modgenes.spring.2019.grey, bias.data = spring.2019.gene.lengths.vector)
# In pcls(G) : initial point very close to some inequality constraints
nullp.result.spring.2019.dkred <- nullp(DEgenes = modgenes.spring.2019.dkred, bias.data = spring.2019.gene.lengths.vector)
nullp.result.fall.2019.ryblue <- nullp(DEgenes = modgenes.fall.2019.ryblue, bias.data = fall.2019.gene.lengths.vector)
nullp.result.spring.2017.HF.tan <- nullp(DEgenes = modgenes.spring.2017.HF.tan, bias.data = spring.2017.HF.gene.lengths.vector)
nullp.result.fall.2017.HF.green <- nullp(DEgenes = modgenes.fall.2017.HF.green, bias.data = fall.2017.HF.gene.lengths.vector)
# In pcls(G) : initial point very close to some inequality constraints
nullp.result.fall.2017.HF.salmon <- nullp(DEgenes = modgenes.fall.2017.HF.salmon, bias.data = fall.2017.HF.gene.lengths.vector)
nullp.result.fall.2017.HF.purple <- nullp(DEgenes = modgenes.fall.2017.HF.purple, bias.data = fall.2017.HF.gene.lengths.vector)
# In pcls(G) : initial point very close to some inequality constraints
nullp.result.fall.2017.HF.black <- nullp(DEgenes = modgenes.fall.2017.HF.black, bias.data = fall.2017.HF.gene.lengths.vector)
nullp.result.fall.2017.HF.tan <- nullp(DEgenes = modgenes.fall.2017.HF.tan, bias.data = fall.2017.HF.gene.lengths.vector)
# In pcls(G) : initial point very close to some inequality constraints
nullp.result.spring.2017.SERC.ltblue <- nullp(DEgenes = modgenes.spring.2017.SERC.ltblue, bias.data = spring.2017.SERC.gene.lengths.vector)
# In pcls(G) : initial point very close to some inequality constraints

#calculate p-values for each GO term
rownames(nullp.result.summer.2017.turq) <- names(summer.2017.gene.lengths.vector) #because of a bug in nullp()
GO.out.summer.2017.turq <- goseq(pwf = nullp.result.summer.2017.turq, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.summer.2017.ltyell) <- names(summer.2017.gene.lengths.vector) #because of a bug in nullp()
GO.out.summer.2017.ltyell <- goseq(pwf = nullp.result.summer.2017.ltyell, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.summer.2017.dkgrey) <- names(summer.2017.gene.lengths.vector) #because of a bug in nullp()
GO.out.summer.2017.dkgrey <- goseq(pwf = nullp.result.summer.2017.dkgrey, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.summer.2017.green) <- names(summer.2017.gene.lengths.vector) #because of a bug in nullp()
GO.out.summer.2017.green <- goseq(pwf = nullp.result.summer.2017.green, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.fall.2017.cyan) <- names(fall.2017.gene.lengths.vector) #because of a bug in nullp()
GO.out.fall.2017.cyan <- goseq(pwf = nullp.result.fall.2017.cyan, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.fall.2017.mdntblue) <- names(fall.2017.gene.lengths.vector) #because of a bug in nullp()
GO.out.fall.2017.mdntblue <- goseq(pwf = nullp.result.fall.2017.mdntblue, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.spring.2019.grey) <- names(spring.2019.gene.lengths.vector) #because of a bug in nullp()
GO.out.spring.2019.grey <- goseq(pwf = nullp.result.spring.2019.grey, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.spring.2019.dkred) <- names(spring.2019.gene.lengths.vector) #because of a bug in nullp()
GO.out.spring.2019.dkred <- goseq(pwf = nullp.result.spring.2019.dkred, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.fall.2019.ryblue) <- names(fall.2019.gene.lengths.vector) #because of a bug in nullp()
GO.out.fall.2019.ryblue <- goseq(pwf = nullp.result.fall.2019.ryblue, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.spring.2017.HF.tan) <- names(spring.2017.HF.gene.lengths.vector) #because of a bug in nullp()
GO.out.spring.2017.HF.tan <- goseq(pwf = nullp.result.spring.2017.HF.tan, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.fall.2017.HF.green) <- names(fall.2017.HF.gene.lengths.vector) #because of a bug in nullp()
GO.out.fall.2017.HF.green <- goseq(pwf = nullp.result.fall.2017.HF.green, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.fall.2017.HF.salmon) <- names(fall.2017.HF.gene.lengths.vector) #because of a bug in nullp()
GO.out.fall.2017.HF.salmon <- goseq(pwf = nullp.result.fall.2017.HF.salmon, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.fall.2017.HF.purple) <- names(fall.2017.HF.gene.lengths.vector) #because of a bug in nullp()
GO.out.fall.2017.HF.purple <- goseq(pwf = nullp.result.fall.2017.HF.purple, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.fall.2017.HF.black) <- names(fall.2017.HF.gene.lengths.vector) #because of a bug in nullp()
GO.out.fall.2017.HF.black <- goseq(pwf = nullp.result.fall.2017.HF.black, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.fall.2017.HF.tan) <- names(fall.2017.HF.gene.lengths.vector) #because of a bug in nullp()
GO.out.fall.2017.HF.tan <- goseq(pwf = nullp.result.fall.2017.HF.tan, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

rownames(nullp.result.spring.2017.SERC.ltblue) <- names(spring.2017.SERC.gene.lengths.vector) #because of a bug in nullp()
GO.out.spring.2017.SERC.ltblue <- goseq(pwf = nullp.result.spring.2017.SERC.ltblue, gene2cat = go.list, test.cats=c("GO:CC", "GO:BP", "GO:MF"))

#list over-represented GO terms (p < 0.05)
GO.out.summer.2017.turq[GO.out.summer.2017.turq$over_represented_pvalue < 0.05,]
GO.out.summer.2017.ltyell[GO.out.summer.2017.ltyell$over_represented_pvalue < 0.05,]
GO.out.summer.2017.dkgrey[GO.out.summer.2017.dkgrey$over_represented_pvalue < 0.05,]
GO.out.summer.2017.green[GO.out.summer.2017.green$over_represented_pvalue < 0.05,]

GO.out.fall.2017.cyan[GO.out.fall.2017.cyan$over_represented_pvalue < 0.05,]
GO.out.fall.2017.mdntblue[GO.out.fall.2017.mdntblue$over_represented_pvalue < 0.05,]

GO.out.spring.2019.grey[GO.out.spring.2019.grey$over_represented_pvalue < 0.05,]
GO.out.spring.2019.dkred[GO.out.spring.2019.dkred$over_represented_pvalue < 0.05,]

GO.out.fall.2019.ryblue[GO.out.fall.2019.ryblue$over_represented_pvalue < 0.05,]
GO.out.spring.2017.HF.tan[GO.out.spring.2017.HF.tan$over_represented_pvalue < 0.05,]

GO.out.fall.2017.HF.green[GO.out.fall.2017.HF.green$over_represented_pvalue < 0.05,]
GO.out.fall.2017.HF.salmon[GO.out.fall.2017.HF.salmon$over_represented_pvalue < 0.05,]
GO.out.fall.2017.HF.purple[GO.out.fall.2017.HF.purple$over_represented_pvalue < 0.05,]
GO.out.fall.2017.HF.black[GO.out.fall.2017.HF.black$over_represented_pvalue < 0.05,]
GO.out.fall.2017.HF.tan[GO.out.fall.2017.HF.tan$over_represented_pvalue < 0.05,]

GO.out.spring.2017.SERC.ltblue[GO.out.spring.2017.SERC.ltblue$over_represented_pvalue < 0.05,]

write.table(GO.out.summer.2017.turq[GO.out.summer.2017.turq$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/summer.2017.turq.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.summer.2017.ltyell[GO.out.summer.2017.ltyell$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/summer.2017.ltyell.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.summer.2017.dkgrey[GO.out.summer.2017.dkgrey$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/summer.2017.dkgrey.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.summer.2017.green[GO.out.summer.2017.green$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/summer.2017.green.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.fall.2017.cyan[GO.out.fall.2017.cyan$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/fall.2017.cyan.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.fall.2017.mdntblue[GO.out.fall.2017.mdntblue$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/fall.2017.mdntblue.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.spring.2019.grey[GO.out.spring.2019.grey$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/spring.2019.grey.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.spring.2019.dkred[GO.out.spring.2019.dkred$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/spring.2019.dkred.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.fall.2019.ryblue[GO.out.fall.2019.ryblue$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/fall.2019.ryblue.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.spring.2017.HF.tan[GO.out.spring.2017.HF.tan$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/spring.2017.HF.tan.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.fall.2017.HF.green[GO.out.fall.2017.HF.green$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/fall.2017.HF.green.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.fall.2017.HF.salmon[GO.out.fall.2017.HF.salmon$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/fall.2017.HF.salmon.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.fall.2017.HF.purple[GO.out.fall.2017.HF.purple$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/fall.2017.HF.purple.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.fall.2017.HF.black[GO.out.fall.2017.HF.black$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/fall.2017.HF.black.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.fall.2017.HF.tan[GO.out.fall.2017.HF.tan$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/fall.2017.HF.tan.GO.terms.txt", quote = FALSE,col.names = TRUE)
write.table(GO.out.spring.2017.SERC.ltblue[GO.out.spring.2017.SERC.ltblue$over_represented_pvalue < 0.05,],row.names=FALSE,
            file="./Data/GO.terms/spring.2017.SERC.ltblue.GO.terms.txt", quote = FALSE,col.names = TRUE)

write.table(GO.out.summer.2017.turq[GO.out.summer.2017.turq$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.summer.2017.turq.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.summer.2017.ltyell[GO.out.summer.2017.ltyell$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.summer.2017.ltyell.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.summer.2017.dkgrey[GO.out.summer.2017.dkgrey$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.summer.2017.dkgrey.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.summer.2017.green[GO.out.summer.2017.green$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.summer.2017.green.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.fall.2017.cyan[GO.out.fall.2017.cyan$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.fall.2017.cyan.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.fall.2017.mdntblue[GO.out.fall.2017.mdntblue$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.fall.2017.mdntblue.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.spring.2019.grey[GO.out.spring.2019.grey$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.spring.2019.grey.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.spring.2019.dkred[GO.out.spring.2019.dkred$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.spring.2019.dkred.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.fall.2019.ryblue[GO.out.fall.2019.ryblue$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.fall.2019.ryblue.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.spring.2017.HF.tan[GO.out.spring.2017.HF.tan$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.spring.2017.HF.tan.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.fall.2017.HF.green[GO.out.fall.2017.HF.green$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.fall.2017.HF.green.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.fall.2017.HF.salmon[GO.out.fall.2017.HF.salmon$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.fall.2017.HF.salmon.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.fall.2017.HF.purple[GO.out.fall.2017.HF.purple$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.fall.2017.HF.purple.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.fall.2017.HF.black[GO.out.fall.2017.HF.black$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.fall.2017.HF.black.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.fall.2017.HF.tan[GO.out.fall.2017.HF.tan$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.fall.2017.HF.tan.GO.terms.txt", quote = FALSE,col.names = FALSE)
write.table(GO.out.spring.2017.SERC.ltblue[GO.out.spring.2017.SERC.ltblue$over_represented_pvalue < 0.05,1:2],row.names=FALSE,
            file="./Data/GO.terms/small.spring.2017.SERC.ltblue.GO.terms.txt", quote = FALSE,col.names = FALSE)


# fdr correction

GO.out.summer.2017.turq.2 = GO.out.summer.2017.turq %>%
  mutate(fdr.over.pvalue = p.adjust(over_represented_pvalue, method = "fdr"))
# all 1?
GO.out.summer.2017.ltyell.2 = GO.out.summer.2017.ltyell %>%
  mutate(fdr.over.pvalue = p.adjust(over_represented_pvalue, method = "fdr"))
# all 1?
GO.out.spring.2019.grey.2 = GO.out.spring.2019.grey %>%
  mutate(fdr.over.pvalue = p.adjust(over_represented_pvalue, method = "fdr"))
# all 1?
GO.out.spring.2019.dkred.2 = GO.out.spring.2019.dkred %>%
  mutate(fdr.over.pvalue = p.adjust(over_represented_pvalue, method = "fdr"))
# all 1?
GO.out.fall.2019.ryblue.2 = GO.out.fall.2019.ryblue %>%
  mutate(fdr.over.pvalue = p.adjust(over_represented_pvalue, method = "fdr"))
# all 1?
GO.out.spring.2017.HF.tan.2 = GO.out.spring.2017.HF.tan %>%
  mutate(fdr.over.pvalue = p.adjust(over_represented_pvalue, method = "fdr"))
# all 1?
GO.out.spring.2017.SERC.ltblue.2 = GO.out.spring.2017.SERC.ltblue %>%
  mutate(fdr.over.pvalue = p.adjust(over_represented_pvalue, method = "fdr"))
# all 1?

# Open text file that is created and paste into Revigo with default settings. Look at the TreeMap tab or table. Cherry pick terms and blocks
# http://revigo.irb.hr/

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







