library("edgeR")

args <- commandArgs(trailingOnly = TRUE)
countFile <- args[1]
group <- args[2]
outputFile <- args[3]
# nfactor <- args[3]
# outputFile <- args[4]

count <- read.table(countFile, sep = "\t", header = T, row.names = 1)
count <- as.matrix(count)
group <- factor(unlist(strsplit(group, ",")))
# nfactor <- as.double(unlist(strsplit(nfactor, ",")))
design <- model.matrix(~group)
d <- DGEList(counts = count, group = group)
# for (i in 1:length(nfactor)){
# 	d$samples[i,3] <- nfactor[[i]]
# }
d <- calcNormFactors(d)
if (length(group) == 2) 
{
	d$common.dispersion = 0.125
} else {
	d <- estimateDisp(d, design)
}

et <- exactTest(d)
write.table(topTags(et, n=Inf, adjust.method="BH", sort.by="PValue"), file=outputFile)



