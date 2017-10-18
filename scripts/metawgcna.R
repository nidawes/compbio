meta <- read.csv("/home/arcturus/Documents/fall2017/compbio/origcsv/Pvirgatum_photoperiodMetadata.csv")

#Key- 
#Tissue- base = 0 tip = 1
#Timepoint- Day = 0 Night = 1
#Genotype- Upland = 0 Lowland = 1
#Photoperiod- Long = 0 Short = 1

#Assign new Levels prior
levels(meta$tissue) <- c(levels(meta$tissue), 0, 1)
levels(meta$timepoint) <- c(levels(meta$timepoint), 0, 1)
levels(meta$genotype) <- c(levels(meta$genotype), 0, 1)
levels(meta$photoperiod) <- c(levels(meta$photoperiod), 0, 1)
	
#Tissue
meta$tissue[meta$tissue == "Base"] <- 0
meta$tissue[meta$tissue == "Tip"] <- 1

#Timepoint
meta$timepoint[meta$timepoint == "Day"] <- 0
meta$timepoint[meta$timepoint == "Night"] <- 1

#Genotype
meta$genotype[meta$genotype == "VS16" | meta$genotype == "DAC"] <- 0
meta$genotype[meta$genotype == "AP13" | meta$genotype == "WBC"] <- 1

#Photoperiod
meta$photoperiod[meta$photoperiod == "LD"] <- 0
meta$photoperiod[meta$photoperiod == "SD"] <- 1

#Output modified Metadata file
write.csv(meta, file = "/home/arcturus/Documents/fall2017/compbio/origcsv/metadata_wgcna.csv")