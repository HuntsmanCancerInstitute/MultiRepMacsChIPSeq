#!/usr/bin/env Rscript --no-save --no-restore  

install.packages(c("optparse","RColorBrewer","reshape2","ggplot2","pheatmap","UpSetR"), repos="http://cran.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "DESeq2", "normR"))
