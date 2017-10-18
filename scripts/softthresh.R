# Load WGCNA and flashClust libraries every time you open R
library(WGCNA)
library(flashClust)
enableWGCNAThreads()
#perform soft threshold check
#perform WGCNA module creation
#Correlate results to traits

# Uploading data into R and formatting it for WGCNA
# This creates an object called "datExpr" that contains the normalized counts file output from DESeq2
datExpr = read.csv("/home/arcturus/Documents/fall2017/compbio/dds_csv/deseq2_dds_allsamp.csv")
# "head" the file to preview it
head(datExpr) # You see that genes are listed in a column named "X" and samples are in columns

 
# Manipulate file so it matches the format WGCNA needs
row.names(datExpr) = datExpr$X
datExpr$X = NULL
datExpr = as.data.frame(t(datExpr)) # now samples are rows and genes are columns
dim(datExpr) # 48 samples and 1000 genes (you will have many more genes in reality)
 
 
# Run this to check if there are gene outliers
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
 
 
#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
if (!gsg$allOK)
   {if (sum(!gsg$goodGenes)>0)
       printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
       if (sum(!gsg$goodSamples)>0)
           printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
       datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
       }
 
 
#Create an object called "datTraits" that contains your trait data
datTraits = read.csv("/home/arcturus/Documents/fall2017/compbio/origcsv/Pvirgatum_photoperiodMetadata.csv")
head(datTraits)
#form a data frame analogous to expression data that will hold the clinical traits.
rownames(datTraits) = datTraits$Sample
datTraits$Sample = NULL
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)
 
 
# You have finished uploading and formatting expression and trait data
# Expression data is in datExpr, corresponding traits are datTraits
 
 
#save(datExpr, datTraits, file="SamplesAndTraits.RData")

#load("SamplesAndTraits.RData")
powers <- c(12:25)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed" )
sizeGrWindow(9, 5)
#20 returns .896