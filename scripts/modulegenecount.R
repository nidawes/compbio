load("Network_allSamples_signed_nomerge_RLDfiltered.RData")
load("Network_allSamples_signed_RLDfiltered.RData")
modcolors <- c()
for(i in colnames(MEs)) {
	modcolors <- paste(modcolors, substr(colnames(MEs)[i], 3, nchar(colnames(MEs)[i])), sep = "")
}
modcolors


#for(i in MEList$validColors) {
	