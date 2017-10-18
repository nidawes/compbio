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
datTraits = read.csv("/home/arcturus/Documents/fall2017/compbio/origcsv/metadata_wgcna.csv")
head(datTraits)
#form a data frame analogous to expression data that will hold the clinical traits.
rownames(datTraits) = datTraits$sampleGroup
datTraits$sampleGroup = NULL
datTraits$location = NULL
datTraits$library = NULL
datTraits$grp = NULL
datTraits$sampleReplicate = NULL
datTraits$sampleID = NULL
datTraits$sampleTissue = NULL
datTraits$grpreplicate = NULL
datTraits$sampleName = NULL
datTraits$replicate = NULL
datTraits$lib. = NULL
datTraits$X = NULL
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)
 
 
# You have finished uploading and formatting expression and trait data
# Expression data is in datExpr, corresponding traits are datTraits
 
 
save(datExpr, datTraits, file="/home/arcturus/Documents/fall2017/compbio/rdata/SamplesAndTraits.RData")

#load("SamplesAndTraits.RData")
#powers <- c(18:25)
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed" )
#sizeGrWindow(9, 5)
#20 returns .896 without cutoff, 17 returns 0.894 with cutoff

#build a adjacency "correlation" matrix
softPower = 17
adjacency = adjacency(datExpr, power = softPower, type = "signed") #specify network type
head(adjacency)
 
# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM
 
# Generate Modules --------------------------------------------------------
 
 
# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
#This sets the minimum number of genes to cluster into a module
minModuleSize = 75
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = 17)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
save(dynamicMods, MEList, MEs, MEDiss, METree, file= "/home/arcturus/Documents/fall2017/compbio/rdata/Network_allSamples_signed_RLDfiltered.RData")
 
 
#plots tree showing how the eigengenes cluster together
png(file="/home/arcturus/Documents/fall2017/compbio/figures/clusterwithoutmodulecolors.png")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
MEDissThres = 0
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
dev.off()
 
#plot dendrogram with module colors below it
png(file="/home/arcturus/Documents/fall2017/compbio/figures/cluster.png", width = 2000, height = 1000, units = "px")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
dev.off()
 
save(MEs, moduleLabels, moduleColors, geneTree, file= "/home/arcturus/Documents/fall2017/compbio/rdata/Network_allSamples_signed_nomerge_RLDfiltered.RData")

# Correlate traits --------------------------------------------------------
 
 
#Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
 
 
#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))
 
 
#display the corelation values with a heatmap plot
png(file="/home/arcturus/Documents/fall2017/compbio/figures/heatmap.png", width = 1200, height = 1200, units = "px")
labeledHeatmap(Matrix= moduleTraitCor,
            xLabels= names(datTraits),
            yLabels= names(MEs),
            ySymbols= names(MEs),
            colorLabels= FALSE,
            colors= blueWhiteRed(50),
            textMatrix= textMatrix,
            setStdMargins= FALSE,
            cex.text= 0.8,
            zlim= c(-1,1),
            main= paste("Module-trait relationships"))
dev.off()