#### Evaluating relationships between genotypes and growth ####

library(vcfR)
library(adegenet)
library(dartR)
## http://georges.biomatix.org/storage/app/media/uploaded-files/dartrguidetopreparatoryanalysis12.pdf
library(vegan)

# building my own distance matrix from the vcf file
# also going to try Uzay's IBS distance matrix from Tassle

#### vcftools code ####

## Filters
# remove indels
# minimum allele frequency of 0.05 (5%)
# maximum missing data 50%
# minimum read depth 5

#cd /Users/samanthaworthy/Desktop/Fagus/Data
#vcftools --gzvcf Fg-sort-indel-rm.vcf.gz --remove-indels --minDP 5  
#--maf 0.1 --max-missing 0.5  --stdout --recode --recode-INFO-all | gzip -c > Gf-sort-indel-rm-maf01.vcf.gz

### Read in vcf ####
# read vcf into R with filtered snps
vcf.fagus=read.vcfR("./Formatted.Data/Gf-sort-indel-rm-maf01.vcf.gz") 
# 55,989 snps

# convert vcf into genlight objects

fagus.genlight=vcfR2genlight(vcf.fagus) #567 loci with more than 2 alleles removed
# 40 genotypes, 55,422 binary SNPs, 29.18% missing data

# Assign population to samples
pop(fagus.genlight)=as.factor(c("HF","HF","HF","HF","HF","HF","HF","HF",
                                "HF","HF","HF","HF","HF","HF","HF","HF","HF",
                                "HF","HF","HF","SERC","SERC","SERC","SERC","SERC","SERC",
                                "SERC","SERC","SERC","SERC","SERC","SERC","SERC",
                                "SERC","SERC","SERC","SERC","SERC","SERC","SERC"))
popNames(fagus.genlight) # check population names
ploidy(fagus.genlight)=2 # set ploidy

HF.genlight=fagus.genlight[1:20,] # HF only genlight
SERC.genlight=fagus.genlight[21:40,] # SERC only genlight

# checks genlight object for compliance with dartR
fagus.genlight <- gl.compliance.check(fagus.genlight)
# data okay
HF.genlight <- gl.compliance.check(HF.genlight) 
# contains monomorphic loci
HF.genlight.2 = gl.filter.monomorphs(HF.genlight) # remove monomorphic loci and loci with all missing data
HF.genlight.2 <- gl.compliance.check(HF.genlight.2) 
# data okay
SERC.genlight <- gl.compliance.check(SERC.genlight)
# contains monomorphic loci, contains loci with all missing data
SERC.genlight.2 = gl.filter.monomorphs(SERC.genlight) # data include loci that are score NA across all individuals
SERC.genlight.2 = gl.filter.allna(SERC.genlight.2)
SERC.genlight.2 = gl.compliance.check(SERC.genlight.2)
# data okay

# smear plot
gl.smearplot(fagus.genlight)
gl.smearplot(HF.genlight.2)
gl.smearplot(SERC.genlight.2)

#### Genetic Distance #####
all.gen.dist = gl.dist.ind(fagus.genlight, method = "euclidean")
HF.gen.dist=gl.dist.ind(HF.genlight.2, method = "euclidean")
SERC.gen.dist=gl.dist.ind(SERC.genlight.2, method = "euclidean")

#### Heatmaps ####
gl.plot.heatmap(all.gen.dist) 
gl.plot.heatmap(HF.gen.dist) 
gl.plot.heatmap(SERC.gen.dist) 

#### PCOA #####
pcoa.all <- gl.pcoa(fagus.genlight)
gl.pcoa.plot(pcoa.all, fagus.genlight) 

pcoa.HF <- gl.pcoa(HF.genlight.2) 
gl.pcoa.plot(pcoa.HF, HF.genlight.2) 

pcoa.SERC <- gl.pcoa(SERC.genlight.2)
gl.pcoa.plot(pcoa.SERC, SERC.genlight.2) 


# remove samples where we don't have growth

#### Growth Distance ####


# Physical Distance
trait.subset=read.csv("Fagus.subset.csv", header = T)
HF.coords=trait.subset[1:20, c(5:7, 11:14)]
SERC.coords=trait.subset[21:40, c(5:7, 11:14)]
row.names(HF.coords)=HF.coords[,2]
row.names(SERC.coords)=SERC.coords[,2]

# physical distance for UTMs

HF.phy.dist=vegdist(HF.coords[,4:5], method = "euclidean")
SERC.phy.dist=vegdist(SERC.coords[,4:5], method="euclidean")

write.csv(as.matrix(HF.phy.dist), file="HF.phy.dist.csv")
write.csv(as.matrix(SERC.phy.dist), file="SERC.phy.dist.csv")


# alternative methods to get physical distance not used
library(raster)
HF.coord.mat1=as.matrix(HF.coords[4:5]) # make UTMs into matrix
HF.dst.phy=pointDistance(HF.coord.mat1,lonlat=FALSE, type="Euclidean")
row.names(HF.dst.phy)=row.names(HF.coord.mat1)
colnames(HF.dst.phy)=row.names(HF.coord.mat1)

SERC.coord.mat1=as.matrix(SERC.coords[4:5])
SERC.dst.phy=pointDistance(SERC.coord.mat1,lonlat=FALSE, type="Euclidean")
row.names(SERC.dst.phy)=row.names(SERC.coord.mat1)
colnames(SERC.dst.phy)=row.names(SERC.coord.mat1)

install.packages("gdata")
library(gdata)

upperTriangle(HF.dst.phy)=NA
upperTriangle(SERC.dst.phy)=NA

# Another way to calculate distance for UTMs
coordinates(HF.coords)=~Easting+Northing # make into spatialdataframe
proj4string(HF.coords) <- tmaptools::get_proj4("utm18")$proj4string #change utm to zone 18
distmatrix.2=gDistance(HF.coords, byid=TRUE) # create distance matrix

# Calculate great circle distance for long/lat
library(fields)
HF.coord.mat1=as.matrix(HF.coords[7:6])
rownames(HF.coord.mat1)=row.names(HF.coords)
SERC.coord.mat1=as.matrix(SERC.coords[7:6])
rownames(SERC.coord.mat1)=row.names(SERC.coords)
HF.great.circle=rdist.earth(HF.coord.mat1)
SERC.great.circle=rdist.earth(SERC.coord.mat1)

upperTriangle(HF.great.circle)=NA
upperTriangle(SERC.great.circle)=NA

# Calculate trait Distances
setwd("~/Documents/Fagus/Data")
trait.subset=read.csv("Fagus.Subset.csv", header=T)
HF.traits=trait.subset[1:20,]
SERC.traits=trait.subset[21:40,]

# Check for normality of traits
shapiro.test(HF.traits$Average.SLA)
shapiro.test(HF.traits$Average.SPAD)
shapiro.test(SERC.traits$Average.SLA)
shapiro.test(SERC.traits$Average.SPAD)
shapiro.test(HF.traits$DBH)
shapiro.test(SERC.traits$DBH) # not normal

install.packages("usedist")
library(usedist)

# Make trait distance matrices

HF.trait.dist=vegdist(HF.traits[,32:33], method="euclidean")
HF.trait.dist=dist_setNames(HF.trait.dist, HF.traits$Sample.ID)

HF.mat.trait.all=as.matrix(HF.dist.trait.all)
row.names(HF.mat.trait.all)=HF.traits$Sample.ID
colnames(HF.mat.trait.all)=HF.traits$Sample.ID

SERC.trait.dist=vegdist(SERC.traits[,32:33], method="euclidean")
SERC.trait.dist=dist_setNames(SERC.trait.dist, SERC.traits$Sample.ID)

SERC.mat.trait.all=as.matrix(SERC.dist.trait.all)
row.names(SERC.mat.trait.all)=SERC.traits$Sample.ID
colnames(SERC.mat.trait.all)=SERC.traits$Sample.ID

HF.dist.SLA=vegdist(HF.traits[,32], method="euclidean")
HF.dist.SLA=dist_setNames(HF.dist.SLA, HF.traits$Sample.ID)

HF.mat.SLA=as.matrix(HF.dist.SLA)
row.names(HF.mat.SLA)=HF.traits$Sample.ID
colnames(HF.mat.SLA)=HF.traits$Sample.ID

HF.dist.DBH=vegdist(HF.traits$DBH, method="euclidean")
HF.dist.DBH=dist_setNames(HF.dist.DBH, HF.traits$Sample.ID)

HF.dist.RGR=vegdist(HF.traits$RGR, method = "euclidean")
HF.dist.RGR=dist_setNames(HF.dist.RGR, HF.traits$Sample.ID)

SERC.dist.SLA=vegdist(SERC.traits[,32], method="euclidean")
SERC.dist.SLA=dist_setNames(SERC.dist.SLA, SERC.traits$Sample.ID)

SERC.mat.SLA=as.matrix(SERC.dist.SLA)
row.names(SERC.mat.SLA)=SERC.traits$Sample.ID
colnames(SERC.mat.SLA)=SERC.traits$Sample.ID

SERC.dist.DBH=vegdist(SERC.traits$DBH, method="euclidean")
SERC.dist.DBH=dist_setNames(SERC.dist.DBH, SERC.traits$Sample.ID)

SERC.dist.RGR=vegdist(SERC.traits$RGR, method="euclidean")
SERC.dist.RGR=dist_setNames(SERC.dist.RGR, SERC.traits$Sample.ID)

write.csv(as.matrix(HF.trait.dist), file="HF.trait.dist.csv")
write.csv(as.matrix(SERC.trait.dist), file="SERC.trait.dist.csv")
write.csv(as.matrix(HF.dist.SLA), file="HF.dist.SLA.csv")
write.csv(as.matrix(SERC.dist.SLA), file="SERC.dist.SLA.csv")
write.csv(as.matrix(HF.dist.DBH), file="HF.dist.DBH.csv")
write.csv(as.matrix(SERC.dist.DBH), file="SERC.dist.DBH.csv")
write.csv(as.matrix(HF.dist.RGR), file="HF.dist.RGR.csv")
write.csv(as.matrix(SERC.dist.RGR), file="SERC.dist.RGR.csv")

# Prepare data.frame
all.dist.mat=matrix(data=NA, nrow=190, ncol=11)
colnames(all.dist.mat)=c("HF.gen", "HF.traits", "HF.SLA", "HF.DBH",
                         "HF.spat", "SERC.gen", "SERC.traits", "SERC.SLA",
                         "SERC.DBH", "SERC.spat", "new.SERC.gen")
HF.dist.gen.table=as.table(HF.gen.dist)
all.dist.mat[,1]=HF.dist.gen.table
row.names(all.dist.mat)=row.names(HF.dist.gen.table)
SERC.dist.gen.table=as.table(SERC.gen.dist)
all.dist.mat[,6]=SERC.dist.gen.table
HF.dist.trait.all.table=as.table(HF.trait.dist)
all.dist.mat[,2]=HF.dist.trait.all.table
HF.dist.SLA.table=as.table(HF.dist.SLA)
all.dist.mat[,3]=HF.dist.SLA.table
SERC.dist.trait.all.table=as.table(SERC.trait.dist)
all.dist.mat[,7]=SERC.dist.trait.all.table
SERC.dist.SLA.table=as.table(SERC.dist.SLA)
all.dist.mat[,8]=SERC.dist.SLA.table
HF.dist.DBH.table=as.table(HF.dist.DBH)
all.dist.mat[,4]=HF.dist.DBH
SERC.dist.DBH.table=as.table(SERC.dist.DBH)
all.dist.mat[,9]=SERC.dist.DBH
HF.phy.dist.table=as.table(HF.phy.dist)
all.dist.mat[,5]=HF.phy.dist.table
SERC.phy.dist.table=as.table(SERC.phy.dist)
all.dist.mat[,10]=SERC.phy.dist.table
new.SERC.dist.gen.table=as.table(new.SERC.gen.dist)
all.dist.mat[,11]=new.SERC.dist.gen.table
colnames(all.dist.mat)[11]="new.SERC.gen"
write.csv(all.dist.mat, file="all.dist.mat.2.csv")

HF.dist.DBH.table=as.table(HF.dist.DBH)
write.csv(HF.dist.DBH.table, file="new.HF.dist.DBH.csv")
SERC.dist.DBH.table=as.table(SERC.dist.DBH)
write.csv(SERC.dist.DBH.table, file="new.SERC.dist.DBH.csv")
HF.dist.RGR.table=as.table(HF.dist.RGR)
write.csv(HF.dist.RGR.table, file="new.HF.dist.RGR.csv")
SERC.dist.RGR.table=as.table(SERC.dist.RGR)
write.csv(SERC.dist.RGR.table, file="new.SERC.dist.RGR.csv")

all.dist.mat.3=read.csv("all.dist.mat.2.csv", header=T, row.names = 1)

# Correlations that are significant
all.dist.df=as.data.frame(all.dist.mat.3)
cor.test(all.dist.df$SERC.gen,all.dist.df$SERC.spat) # p = 0.01197 r = 0.18
cor.test(all.dist.df$SERC.gen,all.dist.df$SERC.traits) # p = 0.1018 r = 0.12
cor.test(all.dist.df$SERC.gen,all.dist.df$SERC.SLA) # p = 0.1032 r = 0.12
cor.test(all.dist.df$SERC.gen,all.dist.df$SERC.DBH) # p < 0.0001 r = 0.32
cor.test(all.dist.df$SERC.spat,all.dist.df$SERC.traits) # p = 0.00018, r = 0.27
cor.test(all.dist.df$SERC.spat,all.dist.df$SERC.SLA) # p = 0.00018, r = 0.27
cor.test(all.dist.df$SERC.spat, all.dist.df$SERC.DBH) # p = 0.1992, r = 0.09
cor.test(all.dist.df$SERC.DBH, all.dist.df$SERC.traits) # p = 0.04059, r = 0.15
cor.test(all.dist.df$SERC.RGR, all.dist.df$SERC.gen) # p = 0.0846, r = -0.13
cor.test(all.dist.df$SERC.RGR, all.dist.df$SERC.spat) # p = 0.01859, r = 0.17
cor.test(all.dist.df$SERC.RGR, all.dist.df$SERC.traits) # p = 0.3604, r = 0.07
cor.test(all.dist.df$SERC.RGR, all.dist.df$SERC.SLA) # p = 0.3646, r = 0.07


cor.test(all.dist.df$HF.gen,all.dist.df$HF.spat) # p = 0.1168 r = -0.11
cor.test(all.dist.df$HF.gen,all.dist.df$HF.traits) # p = 0.3207 r = -0.07
cor.test(all.dist.df$HF.gen,all.dist.df$HF.SLA) # p = 0.3195 r = -0.07
cor.test(all.dist.df$HF.gen,all.dist.df$HF.DBH) # p = 0.0783 r = 0.13
cor.test(all.dist.df$HF.spat,all.dist.df$HF.traits) # p = 0.4437, r = -0.06
cor.test(all.dist.df$HF.spat,all.dist.df$HF.SLA) # p = 0.4511, r = -0.05
cor.test(all.dist.df$HF.spat, all.dist.df$HF.DBH) # p = 0.4606, r = 0.05
cor.test(all.dist.df$HF.DBH, all.dist.df$HF.traits) # p = 0.014, r = 0.18
cor.test(all.dist.df$HF.RGR, all.dist.df$HF.gen) # p = 0.0288, r = 0.17
cor.test(all.dist.df$HF.RGR, all.dist.df$HF.spat) # p = 0.02822, r = 0.16
cor.test(all.dist.df$HF.RGR, all.dist.df$HF.traits) # p = 0.966, r = -0.003
cor.test(all.dist.df$HF.RGR, all.dist.df$HF.SLA) # p = 0.9555, r = -0.004

glPlot(fagus.genlight)
glPlot(HF.genlight)
glPlot(new.SERC.genlight)


### Generate dendrograms

dendro=hclust(all.gen.dist, method="average")
dendro.ward=hclust(all.gen.dist, method = "complete")

plot(dendro)
plot(dendro.ward)
plot(dendro, hang=-1)
plot.hclust(dendro)
plot(all.gen.dist)

all.gen.mat=as.matrix(all.gen.dist)
library(gplots)
heatmap.2(all.gen.mat)

library(fields)

dim=ncol(all.gen.mat)
image.plot(1:dim, 1:dim, all.gen.mat,axes = FALSE, xlab="", ylab="", col = hcl.colors(12, "Mint",rev = TRUE), legend.lab="Distance")
axis(1, 1:dim, colnames(all.gen.mat), cex.axis = 0.5, las=3)
axis(2, 1:dim, colnames(all.gen.mat), cex.axis = 0.5, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", all.gen.mat), cex=0.2)


#### Unused ####
# check later
rm.loci.SERC = c("S120_524834","S120_528576","S120_571139","S120_577009","S120_577098",
                 "S120_577137","S120_578521","S120_578530","S120_580853","S120_580859",
                 "S120_581394","S276_318559","S386_743734","S738_2388749","S1450_404182")

new.SERC.genlight=gl.drop.loc(SERC.genlight, loc.list = rm.loci.SERC)
