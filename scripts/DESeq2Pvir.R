library("DESeq2")
library("ggplot2")
library("pheatmap")
library("gridExtra")

#Use "deseqRawXX" for read_counts
read_counts <- read.csv("/home/arcturus/Documents/fall2017/compbio/origcsv/deseqRawBD.csv", header=TRUE, row.names=1) #Replace with raw data permutations
#Use "deseqMetaXX" for read_metadata
read_metadata <- read.csv("/home/arcturus/Documents/fall2017/compbio/origcsv/deseqMetaBD.csv", header=TRUE, row.names=1) #Replace with metadata permutations

#definte DESeq dataset, specify design, perform DESeq function
dds <- DESeqDataSetFromMatrix(read_counts, read_metadata, ~genotype + photoperiod + genotype:photoperiod)         
#dds <- DESeq(dds) #Redundant probably?
dds$group <- factor(paste0(dds$genotype, dds$photoperiod))
design(dds) <- ~group #design and dds overwritten with permutations
dds <- DESeq(dds)

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
  sig_genes <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.05,]  
  dds2 <- dds[row.names(sig_genes)] #dds2 = dds redefined using only genes found with p value
  res <- results(dds2, contrast = c('group', paste(geno, 'SD', sep=""), paste(geno, 'LD', sep=""))) 
  resOrdered <- res[order(res$padj),]
  write.csv(resOrdered, "/home/arcturus/Documents/fall2017/compbio/deseqresults/DACBD_sig.csv") #Adjust Timepoint/Tissue combination and file path
  #write.csv(resOrdered, "deseq2_AP13SD_v_LD_TN.csv") #write csv from ordered results
  
}
sig_all_func <- function(){
  ap13_sig_genes <- read.csv("/home/arcturus/Documents/fall2017/compbio/deseqresults/AP13BD_sig.csv", header=TRUE, row.names=1)
  wbc_sig_genes <- read.csv("/home/arcturus/Documents/fall2017/compbio/deseqresults/WBCBD_sig.csv", header=TRUE, row.names=1) 
  vs16_sig_genes <- read.csv("/home/arcturus/Documents/fall2017/compbio/deseqresults/VS16BD_sig.csv", header=TRUE, row.names=1)
  dac_sig_genes <- read.csv("/home/arcturus/Documents/fall2017/compbio/deseqresults/DACBD_sig.csv", header=TRUE, row.names=1)
  
  sig_all_vect <- c()
  for (gene in rownames(dac_sig_genes)) { # For every gene in dac sig list
    if (gene %in% rownames(ap13_sig_genes) & gene %in% rownames(vs16_sig_genes) &
        gene %in% rownames(wbc_sig_genes)) {
      sig_all_vect = c(sig_all_vect, gene)
    }
  }
  sig_all_vect
  write.csv(sig_all_vect, "/home/arcturus/Documents/fall2017/compbio/deseqresults/sig_all_p05_nolog.csv")
} #Uses lists of sig p<0.05 genes for each genotype to find similar genes across all
exp_hist_func <- function(dds2){
  eco_mean_diff <- numeric(0)

    samplenames <- colnames(dds2)
    tissue = 'B'
    timepoint = 'D'
    geno <- c('A', 'W', 'V', 'D')
    grep_samp <- paste("[SL]", timepoint, geno[1], tissue, "*", sep="")
    grep_a <- grep(grep_samp, samplenames, perl=TRUE, value=TRUE)

    grep_samp <- paste("[SL]", timepoint, geno[2], tissue, "*", sep="")
    grep_w <- grep(grep_samp, samplenames, perl=TRUE, value=TRUE)

    grep_samp <- paste("[SL]", timepoint, geno[3], tissue, "*", sep="")
    grep_v <- grep(grep_samp, samplenames, perl=TRUE, value=TRUE)

    grep_samp <- paste("[SL]", timepoint, geno[4], tissue, "*", sep="")
    grep_d <- grep(grep_samp, samplenames, perl=TRUE, value=TRUE)

  for (curr_gene in deg_list) {
    d_AP13 <- plotCounts(dds2[,grep_a], gene= curr_gene, intgroup = 'photoperiod', returnData=TRUE)
    d_WBC <- plotCounts(dds2[,grep_w], gene=curr_gene, intgroup = 'photoperiod', returnData=TRUE)
    d_VS16 <- plotCounts(dds2[,grep_v], gene=curr_gene, intgroup = 'photoperiod', returnData=TRUE)
    d_DAC <- plotCounts(dds2[,grep_d], gene=curr_gene, intgroup = 'photoperiod', returnData=TRUE)

    lm_a <- lm(d_AP13$count ~ d_AP13$photoperiod)
    lm_w <- lm(d_WBC$count ~ d_WBC$photoperiod)
    lm_v <- lm(d_VS16$count ~ d_VS16$photoperiod)
    lm_d <- lm(d_DAC$count ~ d_DAC$photoperiod)

    lm_a_slope <- lm_a$coefficients[2]
    lm_w_slope <- lm_w$coefficients[2]
    lm_v_slope <- lm_v$coefficients[2]
    lm_d_slope <- lm_d$coefficients[2]
    
    up_mean <- ((lm_v_slope + lm_d_slope) / 2)
    low_mean <- ((lm_a_slope + lm_w_slope) / 2)
    eco_mean_diff <- c(eco_mean_diff, log10(abs(up_mean - low_mean)))

#pick out genes >10 diff between slope means
  }
  png(filename="~/Documents/fall2017/compbio/figures/eco_mean_diff.png")
  hist(eco_mean_diff, xlab = "Ecotype Expression Difference", main="Dist. of Interaction Difference Between Ecotype", right=F, breaks= 80)
  dev.off()
} #Generates histogram of difference in graph slope mean across ecotype
deg_slope_func <- function(sig_all_dds){ #Generates new DEG list, requires dataframe input, DDS or other gene list
  deg_list <- c()
  for (gene in rownames(sig_all_dds)) { #START GENE ITERATION FOR LOOP
    gene_name <- gene
    samplenames <- colnames(sig_all_dds)
    tissue = 'B'
    timepoint = 'D'
    geno <- c('A', 'W', 'V', 'D')
    grep_samp <- paste("[SL]", timepoint, geno[1], tissue, "*", sep="")
    grep_a <- grep(grep_samp, samplenames, perl=TRUE, value=TRUE)

    grep_samp <- paste("[SL]", timepoint, geno[2], tissue, "*", sep="")
    grep_w <- grep(grep_samp, samplenames, perl=TRUE, value=TRUE)

    grep_samp <- paste("[SL]", timepoint, geno[3], tissue, "*", sep="")
    grep_v <- grep(grep_samp, samplenames, perl=TRUE, value=TRUE)

    grep_samp <- paste("[SL]", timepoint, geno[4], tissue, "*", sep="")
    grep_d <- grep(grep_samp, samplenames, perl=TRUE, value=TRUE)

    d_AP13 <- plotCounts(sig_all_dds[,grep_a], gene= gene_name, intgroup = 'photoperiod', returnData=TRUE)
    d_WBC <- plotCounts(sig_all_dds[,grep_w], gene=gene_name, intgroup = 'photoperiod', returnData=TRUE)
    d_VS16 <- plotCounts(sig_all_dds[,grep_v], gene=gene_name, intgroup = 'photoperiod', returnData=TRUE)
    d_DAC <- plotCounts(sig_all_dds[,grep_d], gene=gene_name, intgroup = 'photoperiod', returnData=TRUE)
    lm_a <- lm(d_AP13$count ~ d_AP13$photoperiod)
    lm_w <- lm(d_WBC$count ~ d_WBC$photoperiod)
    lm_v <- lm(d_VS16$count ~ d_VS16$photoperiod)
    lm_d <- lm(d_DAC$count ~ d_DAC$photoperiod)

    lm_a_slope <- lm_a$coefficients[2]
    lm_w_slope <- lm_w$coefficients[2]
    lm_v_slope <- lm_v$coefficients[2]
    lm_d_slope <- lm_d$coefficients[2]

    #Add cutoff of 10 between ecotype
    #Key- A = Ap13, B = WBC, C = VS16, D = DAC
    if ((lm_a_slope < 0 & lm_w_slope > 0 & lm_v_slope > 0 & lm_d_slope > 0) | (lm_a_slope > 0 & lm_w_slope < 0 & lm_v_slope < 0 & lm_d_slope < 0) | #A- BCD, A+ BCD
    (lm_a_slope > 0 & lm_w_slope < 0 & lm_v_slope > 0 & lm_d_slope > 0) | (lm_a_slope < 0 & lm_w_slope > 0 & lm_v_slope < 0 & lm_d_slope < 0) | #B- ACD, B+ ACD
    (lm_a_slope > 0 & lm_w_slope > 0 & lm_v_slope < 0 & lm_d_slope > 0) | (lm_a_slope < 0 & lm_w_slope < 0 & lm_v_slope > 0 & lm_d_slope < 0) | #C- ABD C+ ABD
    (lm_a_slope > 0 & lm_w_slope > 0 & lm_v_slope > 0 & lm_d_slope < 0) | (lm_a_slope < 0 & lm_w_slope < 0 & lm_v_slope < 0 & lm_d_slope > 0) | #D- ABC D+ ABC
    (lm_a_slope < 0 & lm_w_slope < 0 & lm_v_slope > 0 & lm_d_slope > 0) | (lm_a_slope > 0 & lm_w_slope > 0 & lm_v_slope < 0 & lm_d_slope < 0) | #AB- CD, AB+ CD
    (lm_a_slope < 0 & lm_w_slope > 0 & lm_v_slope < 0 & lm_d_slope > 0) | (lm_a_slope > 0 & lm_w_slope < 0 & lm_v_slope > 0 & lm_d_slope < 0) | #AC- BD, AC+ BD
    (lm_a_slope < 0 & lm_w_slope > 0 & lm_v_slope > 0 & lm_d_slope < 0) | (lm_a_slope > 0 & lm_w_slope < 0 & lm_v_slope < 0 & lm_d_slope > 0)) { #AD- BC, AD+ BC
        deg_list <- c(deg_list,gene_name)
        print(gene_name)
        }
  } #END GENE ITERATION FOR LOOP

} #Generates list of genes where slope differs in sign in at least one genotype

#resultsfunc(read_counts, read_metadata) #Run Results function to generate deseq results per genotype
#sig_all_func() #run sig funtion to generate list of genes significant in all genotypes

sig_all_vect <- read.csv("/home/arcturus/Documents/fall2017/compbio/deseqresults/sig_all_p05_nolog.csv", header=TRUE, row.names=1) #Read csv of DEGs in all genotypes
sig_all_dds <- dds[row.names(sig_all_vect)] #Subset dds by only genes found in sig_all_vect



deg_slope_func(sig_all_dds)
sig_slope_deg_list
#Use single var hist instead, difference between upland and lowland genotypes
