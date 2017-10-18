library("DESeq2")
library("ggplot2")
#library("pheatmap")
#library("gridExtra")


#Use "deseqRawXX" for read_counts
read_count <- read.csv("/home/arcturus/Documents/fall2017/compbio/origcsv/deseqRawBD.csv", header=TRUE, row.names=1) #Replace with raw data permutations
#Use "deseqMetaXX" for read_metadata
read_meta <- read.csv("/home/arcturus/Documents/fall2017/compbio/origcsv/deseqMetaBD.csv", header=TRUE, row.names=1) #Replace with metadata permutations

resultsfunc <- function(read_counts, read_metadata){
  #definte DESeq dataset, specify design, perform DESeq function
  dds <- DESeqDataSetFromMatrix(read_counts, read_metadata, ~genotype + photoperiod + genotype:photoperiod)         
  #dds <- DESeq(dds) #Redundant probably?
  dds$group <- factor(paste0(dds$genotype, dds$photoperiod))
  design(dds) <- ~group #design and dds overwritten with permutations
  dds <- DESeq(dds)

  #Produce results of a contrast of SD and LD photoperiod within genotype
  geno = 'DAC'
  res <- results(dds, contrast = c('group', paste(geno, 'SD', sep=""), paste(geno, 'LD', sep=""))) 
  resOrdered <- res[order(res$padj),]
  #Produce a new 'dds' that only has genes deemed significant in resOrdered
  sig_genes <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.05,]  #Look at P value usage?
  dds2 <- dds[row.names(sig_genes)] #dds2 = dds redefined using only genes found with p value
  res <- results(dds2, contrast = c('group', paste(geno, 'SD', sep=""), paste(geno, 'LD', sep=""))) 
  resOrdered <- res[order(res$padj),]
  write.csv(resOrdered, "/home/arcturus/Documents/fall2017/compbio/deseqresults/DACBD_sig.csv") #Adjust Timepoint/Tissue combination and file path
  #write.csv(resOrdered, "deseq2_AP13SD_v_LD_TN.csv") #write csv from ordered results
}

sig_all_func <- function(){
  ap13_sig_genes <- read.csv("/home/arcturus/Documents/fall2017/compbio/deseqresults/AP13BD_sig.csv", header=TRUE, row.names=1) #Replace with raw data permutations
  wbc_sig_genes <- read.csv("/home/arcturus/Documents/fall2017/compbio/deseqresults/WBCBD_sig.csv", header=TRUE, row.names=1) #Replace with raw data permutations
  vs16_sig_genes <- read.csv("/home/arcturus/Documents/fall2017/compbio/deseqresults/VS16BD_sig.csv", header=TRUE, row.names=1) #Replace with raw data permutations
  dac_sig_genes <- read.csv("/home/arcturus/Documents/fall2017/compbio/deseqresults/DACBD_sig.csv", header=TRUE, row.names=1) #Replace with raw data permutations

  sig_all_vect <- c()
  for (gene in rownames(dac_sig_genes)) { # For every gene in dac sig list
    if (gene %in% rownames(ap13_sig_genes) & gene %in% rownames(vs16_sig_genes) &
        gene %in% rownames(wbc_sig_genes)) {
        sig_all_vect = c(sig_all_vect, gene)
    }
  }
  sig_all_vect
  write.csv(sig_all_vect, "/home/arcturus/Documents/fall2017/compbio/deseqresults/sig_all_p05_nolog.csv")
}

resrun <- "n" #if y, run Results CSV writing
sig_all_run <- "y"
if (resrun  == "y"){
  resultsfunc(read_count, read_meta)
}
if (sig_all_run == "y"){
  sig_all_func()
}
