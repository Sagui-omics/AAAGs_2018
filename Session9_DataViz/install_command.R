install.packages(c("ggdendro", "ggplot2", "gplots", "grid", "gridExtra", "gtable", "RColorBrewer", "reshape2","RMariaDB", "scales", "tidyverse", "vcfR","viridisLite"))

source("https://bioconductor.org/biocLite.R")
biocLite(c("BSgenome.Hsapiens.UCSC.hg19","GenomicFeatures","GenomicRanges", "Gviz", "rtracklayer","TxDb.Hsapiens.UCSC.hg19.knownGene"))
