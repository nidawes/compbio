datExpr = read.csv("/home/arcturus/Documents/fall2017/compbio/dds_csv/deseq2_dds_allsamp.csv")
load("/home/arcturus/Documents/fall2017/compbio/rdata/Network_allSamples_signed_nomerge_RLDfiltered.RData")
library(WGCNA)

row.names(datExpr) = datExpr$X
datExpr$X = NULL
datExpr = as.data.frame(t(datExpr)) # now samples are rows and genes are columns
dim(datExpr)

gsg = goodSamplesGenes(datExpr, verbose = 3)
genes <- colnames(datExpr)
genenum <- gsg$goodGenes
genedat <- c()
genedat$genesaved <- genenum
genedat$genenumber <- (1:length(genenum))
genedat$genename <- genes

genedat <- genedat