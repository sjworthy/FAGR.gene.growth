library(tidyverse)

#### Read in Growth Data ####
SERC.growth = read_csv("./Formatted.Data/Dendro_FAGR.csv")
HF.growth = read_csv("./Formatted.Data/Parameters_HF.csv")
Sample_Description = read_csv("./Formatted.Data/FAGR.description.csv")

#### subset to just the years and tree IDs we need ####

SERC_samples = subset(Sample_Description, Sample_Description$Site == "SERC")

SERC.growth.sub = SERC.growth %>%
  filter(YEAR %in% c(2017,2018,2019,2020))%>%
  filter(TREE_ID %in% SERC_samples$Tree_ID)

HF_samples = subset(Sample_Description, Sample_Description$Site == "HF")

HF.growth.sub = HF.growth %>%
  filter(YEAR %in% c(2017,2018))%>%
  filter(TREE_ID %in% HF_samples$Tree_ID)

#### Calculate growth ####
# b is the upper diameter of the year
# a is the estimate of the year’s starting diameter 

SERC.growth.sub = SERC.growth.sub %>%
  mutate(growth = b - a)

# need to calculate RGR for HF also
HF.growth.sub = HF.growth.sub %>%
  mutate(growth = b - a,
         RGR = log(b) - log(a))

#write.csv(SERC.growth.sub, file = "./Formatted.Data/SERC.growth.sub.csv")
#write.csv(HF.growth.sub, file = "./Formatted.Data/HF.growth.sub.csv")

#### Merge growth data ####
all.growth.data = full_join(SERC.growth.sub,HF.growth.sub)

# write.csv(all.growth.data, file = "./Formatted.Data/all.growth.data.csv")
