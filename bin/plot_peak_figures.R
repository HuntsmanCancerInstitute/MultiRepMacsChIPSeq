#!/usr/bin/env Rscript

# script for plotting collected peak information

# Timothy J. Parnell, PhD
# Huntsman Cancer Institute
# University of Utah
# Salt Lake City, UT 84112
#  
# This package is free software; you can redistribute it and/or modify
# it under the terms of the Artistic License 2.0.  
# 
# Updated versions of this file may be found in the repository
# https://github.com/HuntsmanCancerInstitute/MultiRepMacsChIPSeq


suppressPackageStartupMessages(library("optparse"))

opts <-  list(
    make_option(c("-i", "--input"), default="NA",
         help="Path and basename to the multirep_macs2_pipeline combined output"),
    make_option(c("-n","--min"), default=-4,
         help="Minimum log2 Fold Enrichment value to plot, default -4"),
    make_option(c("-x","--max"), default=4,
         help="Maximum log2 Fold Enrichment value to plot, default 4"),
    make_option(c("-q","--qmax"), default=30,
         help="Maximum q-value to plot, default 30"),
    make_option(c("-r","--fmax"), default = 4,
         help="Maximum fragment value to plot, default 4"),
    make_option(c("-p","--palette"), default = "Set1",
         help="RColorBrewer palette for samples, default Set1"),
    make_option(c("-f","--format"), default="png",
         help="Format of output file: png, pdf, default png")
)

parser <- OptionParser(option_list=opts, description = "
This script generates a number of heat maps and plots for identified peak calls, 
including the following:
  - Scatter plot of before and after duplication counts versus non-duplicate
  - PCA plot between all sample replicates based on fragment counts in peak
  - Multiple pairwise correlation (Pearson, Spearman, and Euclidean distance) 
    heat maps and clusters of the counts between sample replicates
  - Heat map and cluster of the number of peak intersections
  - Heat map and cluster of jaccard (spatial overlap) statistic between peaks
  - Pie chart of spatial overlap fraction of total merged peak coverage for 
    each sample
  - Bar chart of the fraction of fragment counts in corresponding peaks 
    for each replicate, a measure of ChIP efficiency
  - Heat map of the mean q-value scores for each ChIP over merged peaks
  - Heat map of the mean log2 fold enrichment for each ChIP over merged peaks
    with k-means clustering (6, 8, and 10 clusters)
  - Profile heat map of the fragment density over the midpoint of merged
    peaks for all samples, with and without k-means (4) clustering
  - Profile heat map of the log2 Fold Enrichment over the midpoint of merged
    peaks for all samples, with and without k-means (4) clustering
  - Mean profile line plot of fragment density over merged peaks
  - Mean profile line plot of log2 Fold Enrichment over merged peaks

Samples are identified by color palette: Try Set1, Set2, Set3, Spectral, Dark2, 
or any other named palette in RColorBrewer. Note that excessive sample numbers 
may exceed the number the colors in a given palette (8 to 12). 
")

opt <- parse_args(parser)
if( opt$input == "NA" ){
 print_help(parser)
 quit(status=1)
}


library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

# annotation colors for different numbers of clusters
# http://colorbrewer2.org/#type=diverging&scheme=BrBG&n=4
color4 = c("1" = "#a6611a", "2" = "#dfc27d", "3" = "#80cdc1", "4" = "#018571")
color6 = c("1" = "#8c510a", "2" = "#d8b365", "3" = "#f6e8c3",
           "4" = "#c7eae5", "5" = "#5ab4ac", "6" = "#01665e")
color8 = c("1" = "#8c510a", "2" = "#bf812d", "3" = "#dfc27d",
           "4" = "#f6e8c3", "5" = "#c7eae5", "6" = "#80cdc1",
           "7" = "#35978f", "8" = "#01665e")
color10 = c("1" = "#543005", "2" = "#8c510a", "3" = "#bf812d",
            "4" = "#dfc27d", "5" = "#f6e8c3", "6" = "#c7eae5",
            "7" = "#80cdc1", "8" = "#35978f", "9" = "#01665e",
            "10" = "#003c30")



######## Functions #############

makekmeans <-function(rdata, k) {
  set.seed(100)
  n <- ncol(rdata)
  o <- apply(rdata, 1, mean)
  rdata <- rdata[order(o, decreasing = T),]
  kresult <- kmeans(rdata, k, nstart = 100, iter.max = 100)
  kdata <- cbind(rdata, kresult$cluster)
  o <- order(kdata[, n+1])
  kdata <- kdata[o, ]
  rowAnnot <- data.frame(row.names = rownames(kdata),
                         "Cluster"=as.character(kdata$`kresult$cluster`))
  klist <- list("data" = kdata, "anno" = rowAnnot)
  return(klist)
}

plot_mean_k_hm <-function(rdata, hm_min, hm_max, k, hm_color, rowColor, figmain, outbase) {
  n <- ncol(rdata)
  kdata <- makekmeans(rdata, k)
  
  pheatmap(kdata$data[,1:n],  
           border_color = NA, breaks = seq(hm_min,hm_max,length=255), 
           cluster_rows = F, cluster_cols = T, 
           show_rownames = F, show_colnames = T, 
           color = hm_color,
           annotation_row = kdata$anno, 
           annotation_colors = list("Cluster" = rowColor), 
           main = figmain, 
           filename = paste0(outbase, ".", opt$format), width = 8, height = 10)
  write.table(kdata$anno, file = paste0(outbase, "_klist.txt"), 
              sep = "\t", quote = F)
}

plot_profile_k_hm <-function(rdata, hm_min, hm_max, k, hm_color, rowColor, colAnno, figmain, outbase) {
  n <- ncol(rdata)
  kdata <- makekmeans(rdata, k)
  nms <- as.vector(unique(colAnno[,1]))
  clrs <- brewer.pal(length(nms), opt$palette)
  names(clrs) <- nms
  pheatmap(kdata$data[,1:n], 
           border_color = NA, breaks = seq(hm_min,hm_max,length=255), 
           cluster_rows = F, cluster_cols = F, 
           show_rownames = F, show_colnames = F, 
           color = hm_color, 
           annotation_row = kdata$anno, 
           annotation_col = colAnno, 
           annotation_colors = list("Cluster" = rowColor, "Dataset" = clrs), 
           main = figmain,filename = paste0(outbase, ".png"), width = 8, height = 10)
  write.table(kdata$anno, file = paste0(outbase, "_klist.txt"), 
              sep = "\t", quote = F)
}

plot_profile_hm <-function(rdata, hm_min, hm_max, hm_color, colAnno, figmain, outbase) {
  n <- ncol(rdata)
  o <- apply(rdata, 1, mean)
  rdata <- rdata[order(o, decreasing = T),]
  nms <- as.vector(unique(colAnno[,1]))
  clrs <- brewer.pal(length(nms), opt$palette)
  names(clrs) <- nms
  pheatmap(rdata, 
           border_color = NA, breaks = seq(hm_min,hm_max,length=255), 
           cluster_rows = F, cluster_cols = F, 
           show_rownames = F, show_colnames = F, 
           color = hm_color, 
           annotation_col = colAnno, 
           annotation_colors = list("Dataset" = clrs), 
           main = figmain,filename = paste0(outbase, ".png"), width = 8, height = 10)
}

plot_mean_hm <-function(rdata, hm_min, hm_max, hm_color, figmain, outbase) {
  o <- apply(rdata, 1, mean)
  rdata <- rdata[order(o, decreasing = T),]
  pheatmap(rdata, cluster_cols = T, color = hm_color, 
           border_color = NA, breaks = seq(hm_min,hm_max,length=256), 
           cluster_rows = F, show_rownames = F, show_colnames = T, main = figmain,
           filename = paste0(outbase, ".", opt$format), width = 8, height = 10)
}

makeDist <- function(rdata, cdata, outfile) {
  sampleDists = dist(t(rdata))
  sampleDistMatrix <- as.matrix( sampleDists )
  hmclrs = colorRampPalette( rev(brewer.pal(9, 'Blues')) )(255)
  nms <- as.vector(unique(cdata[,1]))
  anclrs <- brewer.pal(length(nms), opt$palette)
  names(anclrs) <- nms
  pheatmap(sampleDistMatrix, color=hmclrs, annotation_row = cdata,
           annotation_colors = list("Dataset" = anclrs),
           cluster_cols = T, cluster_rows = T,
           main = "Distance Correlation of Sample Counts",
           filename = outfile, width = 10, height = 8)
}

makePearCorr <- function(rdata, cdata, outfile) {
  sampleCorr = cor(rdata, method = "pearson")
  sampleMatrix <- as.matrix( sampleCorr )
  hmclrs = colorRampPalette(brewer.pal(9, 'Blues') )(255)
  nms <- as.vector(unique(cdata[,1]))
  anclrs <- brewer.pal(length(nms), opt$palette)
  names(anclrs) <- nms
  pheatmap(sampleMatrix, color = hmclrs, annotation_row = cdata, 
           annotation_colors = list("Dataset" = anclrs),
           breaks = seq(0,1,length=255),
           cluster_cols = T, cluster_rows = T,
           main = "Pearson Correlation of Sample Counts",
           filename = outfile, width = 10, height = 8)
}

makeSpearCorr <- function(rdata, cdata, outfile) {
  sampleCorr = cor(rdata, method = "spearman")
  sampleMatrix <- as.matrix( sampleCorr )
  hmclrs = colorRampPalette(brewer.pal(9, 'Blues') )(255)
  nms <- as.vector(unique(cdata[,1]))
  anclrs <- brewer.pal(length(nms), opt$palette)
  names(anclrs) <- nms
  pheatmap(sampleMatrix, color = hmclrs, annotation_row = cdata,
           annotation_colors = list("Dataset" = anclrs),
           breaks = seq(0,1,length=255),
           cluster_cols = T, cluster_rows = T,
           main = "Spearman Correlation of Sample Counts",
           filename = outfile, width = 10, height = 8)
}

makePCA <- function(rdata, cdata, outfile) {
  if (ncol(rdata) != nrow(cdata)) {
    cdata <- subset(cdata, row.names(cdata) %in% colnames(rdata))
  }
  n  <- apply(rdata, 1, stats::var)
  x <- head(rdata[order(n, decreasing = T),], 1000)
  p <- prcomp(t(x))
  percentVar <- round(p$sdev^2/sum(p$sdev^2) * 100, 1)
  ggdat <- cbind(cdata, as.data.frame(p$x))
  p1 <- ggplot(ggdat, aes(PC1, PC2, color=Dataset, label=rownames(ggdat))) + 
    geom_point(size=3) + 
    scale_fill_brewer(palette=opt$palette) +
    ggrepel::geom_text_repel() +
    xlab(paste0("PC1: ", percentVar[1],"% variance")) + 
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    ggtitle("Principle Component Analysis of Sample Replicates")
  ggsave(filename = outfile, plot = p1, width = 6, height = 6)
} 

plotMid <- function(gdata, ylabel, figmain, outbase) {
  sampleList <- colnames(gdata)[2:ncol(gdata)]
  gdata_long <- melt(gdata, id="Midpoint", measure.vars = sampleList, 
                     variable.name = "Dataset")
  p <- ggplot(gdata_long, aes(Midpoint, value)) + 
    geom_line(aes(color = Dataset)) + 
    scale_fill_brewer(palette=opt$palette) +
    scale_x_continuous(name = "Relative distance") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    ylab(ylabel) +
    ggtitle(figmain)
  ggsave(plot = p, filename = paste0(outbase, '.', opt$format), width = 6, height = 4)
}





######## Script #############

# count data
countfile <- paste0(opt$input,"_counts.txt")
samplefile <- paste0(opt$input,"_samples.txt")
if(file.exists(countfile)) {
  # read input data
  countdata <- read.table(countfile, header=TRUE, sep="\t", row.names = 1, 
                          check.names = FALSE, na.strings = ".")
  sampledata <- read.table(samplefile, header=TRUE, sep="\t", row.names = 1, 
                           check.names = FALSE)
  
  if ( colnames(countdata)[1] == "Name" ) {
    # exclude this column name
    countdata <- countdata[,2:ncol(countdata)]
  }
  countdata[is.na(countdata)] <- 0
  
  # plot
  makeDist(countdata,sampledata,paste0(opt$input,"_distance.", opt$format))
  makePearCorr(countdata,sampledata,paste0(opt$input,"_pearson.", opt$format))
  makeSpearCorr(countdata,sampledata,paste0(opt$input,"_spearman.", opt$format))
  makePCA(countdata,sampledata,paste0(opt$input,"_PCA.", opt$format))
}



# jaccard table
jaccardfile <- paste0(opt$input,'.jaccard.txt')
if(file.exists(jaccardfile)) {
    jdata <- read.table(jaccardfile,header=TRUE,sep="\t", row.names = 1, check.names = F)
    pheatmap(jdata, col=colorRampPalette(brewer.pal(9, 'Blues'))(255), 
             main = "Spatial Overlap Between Peaks", 
             filename = paste0(opt$input, '.jaccard.', opt$format), 
             width = 8, height = 8)
}



# number of intersections
intersectfile <- paste0(opt$input,'.n_intersection.txt')
if(file.exists(intersectfile)) {
    ndata <- read.table(intersectfile,header=TRUE,sep="\t", 
                        row.names = 1, check.names = F)
    pheatmap(ndata, col=colorRampPalette(brewer.pal(9, 'Greens'))(255), 
             main = "Number of Peak Intersections", 
             filename = paste0(opt$input, '.n_intersection.', opt$format), 
             width = 8, height = 8)
}



# spatial venn pie chart
svennfile <- paste0(opt$input, '.spatialVenn.txt')
if(file.exists(svennfile)) {
  svenndata <- read.table(svennfile, header=TRUE, sep="\t", 
                          check.names = F, na.strings = '.')
  colnames(svenndata)[1] <- "Dataset"
  svenndata <- head(svenndata[order(svenndata$Fraction, decreasing = T),], 11)
  svenndata <- subset(svenndata, svenndata$Fraction >= 0.02)
  p <- ggplot(svenndata, aes(x="", y=Fraction, fill=Dataset)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette=opt$palette) +
    theme_minimal() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid=element_blank(),
          axis.text.x=element_blank()) +
    ggtitle("Spatial occupancy fraction for each sample")
  ggsave(plot = p, filename = paste0(opt$input, '.spatialVenn.pie.', opt$format), 
         width = 6, height = 4)
  
}



# chip efficiency
efffile <- paste0(opt$input, '.chip_efficiency.txt')
if(file.exists(efffile)) {
  effdata <- read.table(efffile, header=TRUE, sep = "\t",
                        check.names = F, na.strings = '.')
  effdata$Replicate <- make.unique(as.character(effdata$Replicate))
  p <- ggplot(effdata, aes(x = Replicate, y = Efficiency)) +
    geom_col(aes(fill = Dataset)) +
    scale_fill_brewer(palette=opt$palette) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Fraction of total fragments in respective called peaks")
  ggsave(plot = p, filename = paste0(opt$input, '.chip_efficiency.', opt$format), 
         width = 6, height = 4)
}



# deduplication plot
ddupfile <- paste0(opt$input, '.dedup-stats.txt')
if(file.exists(ddupfile)) {
  ddupdata <- read.table(ddupfile, header=TRUE, sep="\t", 
                        check.names = F, na.strings = '.')
  colnames(ddupdata)[4] <- "Before"
  colnames(ddupdata)[7] <- "After"
  ddupdata_long <- melt(ddupdata, id="NonDuplicateCount", 
                       measure.vars = c("Before", "After"), 
                       variable.name = "DeDuplication")
  p <- ggplot(ddupdata_long, aes(x = NonDuplicateCount, y = value, color = DeDuplication)) +
    geom_point(size=3) +
    xlab("NonDuplicate Count") +
    ylab("Duplicate Count") +
    ggtitle("Number of Duplicate Fragments After De-duplication")
  ggsave(plot = p, filename = paste0(opt$input, '.dedup-stats.', opt$format),
         width = 6, height = 4)
}



# qvalue heat map table table
qfile <- paste0(opt$input, "_meanQvalue.txt")
if(file.exists(qfile)) {
    qdata = read.table(qfile, header=TRUE, sep="\t", 
                       row.names = 1, check.names = F, na.strings = '.')
    if ( colnames(qdata)[1] == "Name" ) {
        # exclude this column name
        qdata <- qdata[,2:ncol(qdata)]
    }
    if (ncol(qdata > 1)) {
      # only make a heat map file if there is more than one sample
      qdata[is.na(qdata)] <- 0
      plot_mean_hm(qdata, 0, opt$qmax, colorRampPalette(brewer.pal(9, 'Reds'))(255), 
             "Mean Q-Value over merged peaks",
             paste0(opt$input,"_qvalue_hm"))
    }
 }



# log2 Fold Enrichment table
lfefile <- paste0(opt$input,"_meanLog2FE.txt")
if(file.exists(lfefile)) {
    lfedata = read.table(lfefile,header=TRUE,sep="\t", 
                         row.names = 1, check.names = F, na.strings = '.')
    if ( colnames(lfedata)[1] == "Name" ) {
        # exclude this column name
        lfedata <- lfedata[,2:ncol(lfedata)]
    }
    if (ncol(lfedata) > 1) {
      # only make heat map files if there is more than one sample
      lfedata[is.na(lfedata)] <- 0
      clrs <- colorRampPalette(rev(brewer.pal(9, 'RdBu')))(255)
      plot_mean_hm(lfedata, opt$min, opt$max, clrs, 
             "Mean Log2 Fold Enrichment over merged peaks",
             paste0(opt$input,"_log2FE_hm"))
      plot_mean_k_hm(lfedata, opt$min, opt$max, 6, clrs, color6, 
                 "Mean Log2 Fold Enrichment over merged peaks, 6 Clusters",
                 paste0(opt$input,"_log2FE_hm_K6"))
      plot_mean_k_hm(lfedata, opt$min, opt$max, 8, clrs, color8, 
                 "Mean Log2 Fold Enrichment over merged peaks, 8 Clusters",
                 paste0(opt$input,"_log2FE_hm_K8"))
      plot_mean_k_hm(lfedata, opt$min, opt$max, 10, clrs, color10, 
                 "Mean Log2 Fold Enrichment over merged peaks, 10 Clusters",
                 paste0(opt$input,"_log2FE_hm_K10"))
    }
}



# fragment profile heat map
fproffile <- paste0(opt$input, "_profile_fragment.txt")
if(file.exists(fproffile)) {
  fprofdata = read.table(fproffile,header=TRUE,sep="\t", 
                         row.names = 1, check.names = F, na.strings = '.')
  if ( colnames(fprofdata)[1] == "Name" ) {
    # exclude this column name
    fprofdata <- fprofdata[,2:ncol(fprofdata)]
  }
  grpdata = read.table(paste0(opt$input, "_profile_fragment.groups.txt"), header=TRUE, sep="\t", 
                       row.names = 1, check.names = F)
  clrs <- colorRampPalette(brewer.pal(9, 'Reds'))(255)
  plot_profile_hm(fprofdata, 0, opt$fmax, clrs, grpdata,
         "ChIP fragment density profile around merged peak midpoints",
         paste0(opt$input,"_profile_fragment_hm"))
  plot_profile_k_hm(fprofdata, 0, opt$fmax, 4, clrs, color4, grpdata, 
                    "ChIP fragment density profile around merged peak midpoints, 4 Clusters",
                    paste0(opt$input, "_profile_fragment_hm_K4"))
}



# log2FE profile heat map
lfeproffile <- paste0(opt$input, "_profile_log2FE.txt")
if(file.exists(lfeproffile)) {
  lfeprofdata = read.table(lfeproffile,header=TRUE,sep="\t", 
                         row.names = 1, check.names = F, na.strings = '.')
  if ( colnames(lfeprofdata)[1] == "Name" ) {
    # exclude this column name
    lfeprofdata <- lfeprofdata[,2:ncol(lfeprofdata)]
  }
  grpdata = read.table(paste0(opt$input, "_profile_log2FE.groups.txt"), header=TRUE, sep="\t", 
                       row.names = 1, check.names = F)
  clrs <- colorRampPalette(rev(brewer.pal(9, 'RdBu')))(255)
  plot_profile_hm(lfeprofdata, opt$min, opt$max, clrs, grpdata,
         "ChIP Log2 Fold Enrichment profile around merged peak midpoints",
         paste0(opt$input,"_profile_log2FE_hm"))
  plot_profile_k_hm(lfeprofdata, opt$min, opt$max, 4, clrs, color4, grpdata, 
                    "ChIP Log2 Fold Enrichment profile around merged peak midpoints, 4 Clusters",
                    paste0(opt$input, "_profile_log2FE_hm_K4"))
}



# plot fragment mean profile line plot
fprofsumfile <- paste0(opt$input, "_profile_fragment_summary.txt")
if(file.exists(fprofsumfile)) {
  fprofsumdata <- read.table(fprofsumfile, header=TRUE, sep="\t", 
                            row.names = 1, check.names = F, na.strings = '.')
  fprofsumdata[is.na(fprofsumdata)] <- 0
  plotMid(fprofsumdata, "Fragment Density", "Mean fragment density profile around merged peak midpoints", 
          paste0(opt$input, "_profile_fragment_summary"))
  
}



# plot log2FE mean profile line plot
lfeprofsumfile <- paste0(opt$input, "_profile_log2FE_summary.txt")
if(file.exists(fprofsumfile)) {
  lfeprofsumdata <- read.table(lfeprofsumfile, header=TRUE, sep="\t", 
                             row.names = 1, check.names = F, na.strings = '.')
  lfeprofsumdata[is.na(lfeprofsumdata)] <- 0
  plotMid(lfeprofsumdata, "Log2 Fold Enrichment", "Mean Log2 Fold Enrichment profile around merged peak midpoints", 
          paste0(opt$input, "_profile_log2FE_summary"))
  
}


