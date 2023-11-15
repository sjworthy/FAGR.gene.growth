# With low and high growth categories
# think of them as treatments
# follow same WGCNA pathway
# with output modules, just reverse lm

# https://www.youtube.com/watch?v=mzXIxjPr_Mc
# https://deneflab.github.io/HNA_LNA_productivity/WGCNA_analysis.html#31_quantifying_module-trait_associations

# Do eigen genes significantly different among growth levels?
# Do genes or clusters of genes have significant association with high growth trees?

# binarize categorical variables
trait.Data %>%
  mutate(growth.state.binary = ifelse(grep("High", growth_state), 1, 0))
# make new column where high = 1 and low = 0

Spring.2017.MEs.3 <- Spring.2017.MEs.2 %>%
  mutate(
    lm = map(Data, ~ lm(eigen_value ~ growth.state.binary, data = .)),
    lm_glance = map(lm, broom::glance),
    lm_tidy = map(lm, broom::tidy))

Spring.2017.MEs.3 <- Spring.2017.MEs.3 %>%
  mutate(Plot = map(Data, function(.x) {
    ggplot(.x, aes(x = growth.state.binary, y = eigen_value)) +
      geom_point() +
      stat_smooth(method = "lm", col = "blue")
  })) %>%
  unnest(lm_tidy)
filtered_Spring.2017.MEs.3 <- Spring.2017.MEs.3 %>%
  filter(term == "growth.state.binary0") %>%
  arrange(p.value) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  select(module, term, estimate, p.value, fdr) %>%
  filter(fdr < .01)




#### Correlation

# binarize categorical variables
trait.Data %>%
  mutate(growth.state.binary = ifelse(grep("High", growth_state), 1, 0))
# make new column where high = 1 and low = 0

nSamples = nrow(norm.counts) # define number of samples
nGenes = ncol(norm.counts) # define number of genes

# correlation between module eigenegens and growth states
module.trait.corr <- cor(module_eigengenes, trait.Data$growth.state.binary, use = 'p') # pearson correlation
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples) # p-values for correlations

# heat map of correlation
heatmap.data = merge(modules.eigengenes, trait.Data$growth.state.binary, by  = row.names) # combining data into one dataframe
# specificy columns we need, x needs all the trait data columns, y needs all eigengene names 
CorLevelPlot(heatmap.data, x = names(heatmap.data)[19:23], y= names(heatmap.data)[1:18],
             col = c("blue","skyblue","white","pink","red"))
# level of significance indicated by *
# extract genes from modules with significance with high growth

# or use a t.test ME vs high/low
aov(eigen_gene ~ growth.state.binary)

# https://www.polarmicrobes.org/weighted-gene-correlation-network-analysis-wgcna-applied-to-microbial-communities/

# https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-022-03677-8
# Figure 1: correlation matrix
# Tables of genes correlated with high growth
  # Venn diagram of genes associated with each grwoth group
# Tables of genes differentially expressed in high growth individuals
  # plot of these genes using Julin's code
  # Revigo tree map: file:///Users/sjworthy/Documents/GitHub/BIS180L/assignment10-sjworthy/scripts/Assignment_10.html

# methods:
#https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-06809-2
#https://onlinelibrary.wiley.com/doi/10.1111/mec.15289

