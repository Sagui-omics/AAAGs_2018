#### Rscript for creating QC plots before and after filtering
#### Written by Genevieve Housman and M. Nieves Colon - July 2018
#### For AGAR 2018 workshop

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE) 

#Install and load library
library(vcfR)

#Load data
#dir <- "G:/My Drive/Documents/AAAG/plotting data/"
#vcf <- read.vcfR(file = paste0(dir,"filtered_vcf/human_v37_MT.gatk.called.filt_snps.vcf"), verbose = TRUE)
vcf <- read.vcfR(args[1], verbose=T)
state <- args[2]

# Save all plots as pdf
pdf(paste0("filtered_VCF/QCplots/VCF.", state, ".pdf"))

#create chromR object
chrom <- create.chromR(vcf, name = "CHROM", seq = NULL, ann = NULL,  verbose = TRUE)
chrom <- proc.chromR(chrom)

#plot data summaries as histograms or as a function of chromosomal position
plot(chrom)
chromoqc(chrom)

#boxplot of read depth per sample 
dp <- extract.gt(chrom, element="DP", as.numeric=TRUE) 
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Depth", las=2)

#plot heatmap of read depth of each snp in all samples
heatmap.bp(dp)

dev.off()

