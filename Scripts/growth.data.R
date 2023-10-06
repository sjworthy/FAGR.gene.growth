#### Growth Data ####

SERC.growth = read_csv("./Formatted.Data/Dendro_FAGR.csv")
Sample_Description = read_csv("./Formatted.Data/FAGR.description.csv")

# subset to just the years and tree IDs we need

SERC_samples = subset(Sample_Description, Sample_Description$Site == "SERC")

SERC.growth.sub = SERC.growth %>%
  filter(YEAR %in% c(2017,2018,2019,2020))%>%
  filter(TREE_ID %in% SERC_samples$Tree_ID)

write.csv(SERC.growth.sub, file = "./Formatted.Data/SERC.growth.sub.csv")
