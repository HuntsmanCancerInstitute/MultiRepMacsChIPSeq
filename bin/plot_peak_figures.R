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
    make_option(c("--mincluster"), default = 500,
         help="Minimum number of peaks to generate k-means plots, default 500"),
    make_option(c("--spatial"), type="logical", action="store_true", default = FALSE,
         help="Generate spatial UpSet plots (compute intensive), default F"),
    make_option(c("-p","--palette"), default = "Set1",
         help="RColorBrewer palette for samples, default Set1"),
    make_option(c("-f","--format"), default="png",
         help="Format of output file: png, pdf, default png")
)

parser <- OptionParser(option_list=opts, description = "
This script generates a number of heat maps and plots for identified peak calls, 
including the following:
  - Stacked bar chart of unique and duplicate fragment counts
  - Scatter plot of before and after duplication counts versus non-duplicate
  - PCA plot between all sample replicates based on fragment counts in peak
  - Multiple pairwise correlation (Pearson, Spearman, and Euclidean distance) 
    heat maps and clusters of the counts between sample replicates
  - Bar chart plots for peak number and genomic space for each sample peaks
  - Box and whisker plots for peak length distributions
  - Line plot for ranked peak signals for each sample replicate
  - Heat map and cluster of the number of peak intersections
  - Heat map and cluster of jaccard (spatial overlap) statistic between peaks
  - Pie chart of spatial overlap fraction of total merged peak coverage for 
    each sample
  - UpSet plot of peak interactions
  - UpSet plot of peak spatial overlaps (optional, compute intensive)
  - Bar chart of the fraction of fragment counts in corresponding peaks 
    for each replicate, a measure of ChIP efficiency
  - Heat map of the mean q-value scores for each ChIP over merged peaks
  - Heat map of the mean log2 fold enrichment for each ChIP over merged peaks
    with k-means clustering (4, 6, and 8 clusters)
  - Profile heat map of the fragment density over the midpoint of merged
    peaks for all samples, with and without k-means (4) clustering, or
    clustered by sample calls
  - Profile heat map of the log2 Fold Enrichment over the midpoint of merged
    peaks for all samples, with and without k-means (4) clustering, or
    clustered by sample calls
  - Mean profile line plot of fragment density over merged peaks, including
    for each group sample calls
  - Mean profile line plot of log2 Fold Enrichment over merged peaks,
    including for each group of sample calls

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
library(UpSetR)
library(grid)
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
           cluster_rows = FALSE, cluster_cols = TRUE, 
           show_rownames = FALSE, show_colnames = TRUE, 
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
  if (length(nms) == 2) {
    clrs <- c("#984EA3", "#4DAF4A")
  }else{
    clrs <- brewer.pal(length(nms), opt$palette)
  }
  names(clrs) <- nms
  pheatmap(kdata$data[,1:n], 
           border_color = NA, breaks = seq(hm_min,hm_max,length=255), 
           cluster_rows = FALSE, cluster_cols = FALSE, 
           show_rownames = FALSE, show_colnames = FALSE, 
           color = hm_color, 
           annotation_row = kdata$anno, 
           annotation_col = colAnno, 
           annotation_colors = list("Cluster" = rowColor, "Dataset" = clrs), 
           main = figmain,filename = paste0(outbase, ".png"), width = 8, height = 10)
  write.table(kdata$anno, file = paste0(outbase, "_klist.txt"), 
              sep = "\t", quote = F)
}

plot_profile_hm <-function(rdata, hm_min, hm_max, hm_color, colAnno, figmain, outbase) {
  o <- apply(rdata, 1, mean)
  rdata <- rdata[order(o, decreasing = T),]
  nms <- as.vector(unique(colAnno[,1]))
  if (length(nms) == 2) {
    clrs <- c("#984EA3", "#4DAF4A")
  }else{
    clrs <- brewer.pal(length(nms), opt$palette)
  }
  names(clrs) <- nms
  pheatmap(rdata, 
           border_color = NA, breaks = seq(hm_min,hm_max,length=255), 
           cluster_rows = FALSE, cluster_cols = FALSE, 
           show_rownames = FALSE, show_colnames = FALSE, 
           color = hm_color, 
           annotation_col = colAnno, 
           annotation_colors = list("Dataset" = clrs), 
           main = figmain, filename = paste0(outbase, ".png"), width = 8, height = 10)
}

plot_sorted_profile_hm <-function(rdata, hm_min, hm_max, hm_color, rowAnno, colAnno, figmain, outbase) {
  nms <- as.vector(unique(colAnno[,1]))
  if (length(nms) == 2) {
    clrs <- c("#984EA3", "#4DAF4A")
  }else{
    clrs <- brewer.pal(length(nms), opt$palette)
  }
  names(clrs) <- nms
  annoColors <- list("Dataset" = clrs)
  for (i in 1:ncol(rowAnno)) {
    n <- colnames(rowAnno)[i]
    annoColors[[n]] <- c(Y = "red", N = "blue")
  }
  pheatmap(rdata, 
           border_color = NA, breaks = seq(hm_min,hm_max,length=255), 
           cluster_rows = FALSE, cluster_cols = FALSE, 
           show_rownames = FALSE, show_colnames = FALSE, 
           color = hm_color, 
           annotation_col = colAnno, annotation_row = rowAnno,
           annotation_colors = annoColors, 
           main = figmain, filename = paste0(outbase, ".png"), width = 8, height = 10)
}



plot_mean_hm <-function(rdata, hm_min, hm_max, hm_color, figmain, outbase) {
  o <- apply(rdata, 1, mean)
  rdata <- rdata[order(o, decreasing = TRUE),]
  pheatmap(rdata, cluster_cols = TRUE, color = hm_color, 
           border_color = NA, breaks = seq(hm_min,hm_max,length=256), 
           cluster_rows = FALSE, show_rownames = FALSE, show_colnames = TRUE, main = figmain,
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
           cluster_cols = TRUE, cluster_rows = TRUE,
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
           cluster_cols = TRUE, cluster_rows = TRUE,
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
           cluster_cols = TRUE, cluster_rows = TRUE,
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

# max allowed number of columns for palette
allowed <- brewer.pal.info[opt$palette,]$maxcolors

# count data
countfile <- paste0(opt$input,"_counts.txt.gz")
samplefile <- paste0(opt$input,"_samples.txt")
if(file.exists(countfile)) {
  # read input data
  countdata <- read.table(countfile, header=TRUE, sep="\t", row.names = 1, 
                          check.names = FALSE, na.strings = ".")
  sampledata <- read.table(samplefile, header=TRUE, sep="\t", row.names = 1, 
                           check.names = FALSE)

  # discard name
  if ( colnames(countdata)[1] == "Name" ) {
    # exclude this column name
    countdata <- countdata[,2:ncol(countdata)]
  }
  
  # only plot correlations if we have at least 3 samples
  if (ncol(countdata) > 2) {
    countdata[is.na(countdata)] <- 0
    
    # plot
    makeDist(countdata,sampledata,paste0(opt$input,"_distance.", opt$format))
    makePearCorr(countdata,sampledata,paste0(opt$input,"_pearson.", opt$format))
    makeSpearCorr(countdata,sampledata,paste0(opt$input,"_spearman.", opt$format))
    makePCA(countdata,sampledata,paste0(opt$input,"_PCA.", opt$format))
  }
  
  # peak score distributions
  count_long <- melt(countdata, measure.vars = colnames(countdata), 
                     variable.name = "Samples")
  p <- ggplot(count_long, aes(Samples, value)) +
    geom_boxplot() +
    geom_boxplot(outlier.alpha = 0.1) +
    scale_y_log10() +
    coord_flip() +
    ylab("Normalized Signal Count") +
    ggtitle("Peak Signal Count distribution")
  ggsave(plot = p, filename = paste0(opt$input, '_signalDistribution.', opt$format),
         width = 8, height = 4)
  
  
  # ranked peak scores
  d <- data.frame("Rank" = 1:nrow(countdata))
  for (i in 1:ncol(countdata)){
    n <-names(countdata)[i]
    d[n] <- sort(countdata[,i], decreasing = FALSE)
  }
  dlong <- melt(d, id="Rank", 
                measure.vars = names(d)[2:ncol(d)], 
                variable.name = "Samples")
  p <- ggplot(data = dlong, aes(Rank, value)) +
    geom_line(aes(color=Samples), linewidth = 0.25) + 
    scale_fill_brewer(palette=opt$palette) +
    scale_x_continuous(name = "Individual Peak Rank") + 
    scale_y_log10() +
    ylab("Normalized Signal Count") +
    ggtitle("Ranked Peak Signal For Each Sample")
  ggsave(plot = p, filename = paste0(opt$input, '_rankedSignal.', opt$format), 
         width = 8, height = 4)

}



# jaccard table
jaccardfile <- paste0(opt$input,'.jaccard.txt')
if(file.exists(jaccardfile)) {
    jdata <- read.table(jaccardfile,header=TRUE, sep="\t", row.names = 1, 
                        check.names = FALSE, na.strings = ".")
    jdata[is.na(jdata)] <- 0
    pheatmap(jdata, col=colorRampPalette(brewer.pal(9, 'Blues'))(255), 
             main = "Spatial Overlap Between Peaks", 
             breaks = seq(0,1,length=255),
             display_numbers = TRUE, number_format = "%.2f",
             number_color = 'red',
             filename = paste0(opt$input, '.jaccard.', opt$format), 
             width = 8, height = 8)
}



# number of intersections
intersectfile <- paste0(opt$input,'.n_intersection.txt')
if(file.exists(intersectfile)) {
    ndata <- read.table(intersectfile,header=TRUE,sep="\t", 
                        row.names = 1, check.names = FALSE)
    maxval <- max(ndata)
    pheatmap(ndata, col=colorRampPalette(brewer.pal(9, 'Blues'))(255), 
             main = "Number of Peak Intersections", 
             breaks = seq(0,maxval,length=255),
             display_numbers = TRUE, number_format = "%.1e",
             number_color = 'red',
             filename = paste0(opt$input, '.n_intersection.', opt$format), 
             width = 8, height = 8)
}



# Intersection UpSet plots
interfile <- paste0(opt$input, '.intersection.txt')
if (file.exists(interfile)) {
  interdata <- read.table(interfile, header=TRUE, sep="\t", 
                          check.names = FALSE, na.strings = '.')
  
  # change names to be &-delimited strings for UpSetR
  n <- interdata$Peaks
  for (i in 1:length(n)) {
    n[i] <- paste(unlist(strsplit(interdata[i,1], ",")), collapse = "&")
  }
  
  # collect data for upset plots
  intercounts <- as.vector(interdata$Count)
  names(intercounts) <- n
  interspace <- as.vector(interdata$BasePairs)
  names(interspace) <- n
  
  # UpSet plots upset based on output
  if (opt$format == "png") {
    
    png(filename = paste0(opt$input, '.intersection_upset.png'), width = 6, height = 5, units = "in", res = 300)
    print(
      upset(fromExpression(intercounts), order.by = "freq", mainbar.y.label = "Peak Intersection Counts", 
            sets.x.label = "Peak Counts", number.angles = 30, nsets = length(n))
    )
    grid.text("Counts of Peak Intersections",x = 0.65, y=0.97, gp=gpar(fontsize=10))
    dev.off()
    
    if (opt$spatial == TRUE) {
      png(filename = paste0(opt$input, '.spatial_upset.png'), width = 6, height = 5, units = "in", res = 300)
      print(
        upset(data = fromExpression(interspace), order.by = "freq", mainbar.y.label = "Spatial Overlap (bp)", 
              sets.x.label = "Peak Total Size (bp)", number.angles = 30, nsets = length(n))
      )
      grid.text("Spatial Overlap of Peak Intersections",x = 0.65, y=0.97, gp=gpar(fontsize=10))
      dev.off()
    }
  } else if (opt$format == "pdf") {
    
    pdf(file = paste0(opt$input, '.intersection_upset.pdf'), width = 6, height = 5, onefile = F)
    print(
      upset(fromExpression(intercounts), order.by = "freq", mainbar.y.label = "Peak Intersection Counts", 
            sets.x.label = "Peak Counts", number.angles = 30, nsets = length(n))
    )
    grid.text("Counts of Peak Intersections",x = 0.65, y=0.97, gp=gpar(fontsize=10))
    dev.off()
    
    if (opt$spatial == TRUE) {
      pdf(file = paste0(opt$input, '.spatial_upset.pdf'), width = 6, height = 5, onefile = F)
      print(
        upset(data = fromExpression(interspace), order.by = "freq", mainbar.y.label = "Spatial Overlap", 
              sets.x.label = "Peak Total Size", number.angles = 30, nsets = length(n))
      )
      grid.text("Spatial Overlap of Peak Intersections",x = 0.65, y=0.97, gp=gpar(fontsize=10))
      dev.off()
    }
  }
  
  # Spatial intersection pie chart
  svenndata <- subset(interdata, interdata$BasePairFraction >= 0.02)
  svenndata <- head(svenndata[order(svenndata$BasePairFraction, decreasing = T),], 9)
  p <- ggplot(svenndata, aes(x="", y=BasePairFraction, fill=Peaks)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette=opt$palette) +
    theme_classic() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid=element_blank(),
          axis.line = element_blank(),
          axis.text.x=element_blank()) +
    ggtitle("Spatial Fraction for Each Peak Intersection")
  ggsave(plot = p, filename = paste0(opt$input, '.spatial_pie.', opt$format), 
         width = 6, height = 4)
  
  # Count intersection pie chart
  p <- ggplot(svenndata, aes(x="", y=CountFraction, fill=Peaks)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette=opt$palette) +
    theme_classic() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid=element_blank(),
          axis.line = element_blank(),
          axis.text.x=element_blank()) +
    ggtitle("Count Fraction for Each Peak Intersection")
  ggsave(plot = p, filename = paste0(opt$input, '.intersection_pie.', opt$format), 
         width = 6, height = 4)
  
}



# spatial venn pie chart - OLD
svennfile <- paste0(opt$input, '.spatialVenn.txt')
if(file.exists(svennfile)) {
  svenndata <- read.table(svennfile, header=TRUE, sep="\t", 
                          check.names = FALSE, na.strings = '.')
  colnames(svenndata)[1] <- "Dataset"
  svenndata <- subset(svenndata, svenndata$Fraction >= 0.02)
  svenndata <- head(svenndata[order(svenndata$Fraction, decreasing = T),], 9)
  p <- ggplot(svenndata, aes(x="", y=Fraction, fill=Dataset)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_brewer(palette=opt$palette) +
    theme_classic() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid=element_blank(),
          axis.line = element_blank(),
          axis.text.x=element_blank()) +
    ggtitle("Spatial occupancy fraction for each sample")
  ggsave(plot = p, filename = paste0(opt$input, '.spatialVenn.pie.', opt$format), 
         width = 6, height = 4)
  
}



# chip efficiency
efffile <- paste0(opt$input, '.chip_efficiency.txt')
if(file.exists(efffile)) {
  effdata <- read.table(efffile, header=TRUE, sep = "\t",
                        check.names = FALSE, na.strings = '.')
  effdata$Replicate <- make.unique(as.character(effdata$Replicate))
  p <- ggplot(effdata, aes(x = Replicate, y = Efficiency)) +
    geom_col(aes(fill = Dataset)) +
    scale_fill_brewer(palette=opt$palette) +
    coord_flip() +
    ggtitle("Fraction of total fragments in respective called peaks")
  ggsave(plot = p, filename = paste0(opt$input, '.chip_efficiency.', opt$format), 
         width = 6, height = 4)
}



# peak length statistics
lengthfile <- paste0(opt$input, '.lengthStats.txt')
if(file.exists(lengthfile)) {
  lengthdata <- read.table(lengthfile, header=TRUE, sep = "\t",
                           check.names = FALSE, na.strings = '.')
  
  # peak number
  p <- ggplot(lengthdata, aes(x = File, y = Count)) + 
    geom_col() +
    scale_fill_brewer(palette=opt$palette) +
    coord_flip() +
    ggtitle("Number of peaks")
  ggsave(plot = p, filename = paste0(opt$input, '.peak_number.', opt$format), 
         width = 6, height = 4)
  
  # total peak space
  p <- ggplot(lengthdata, aes(x = File, y = Sum)) + 
    geom_col() +
    scale_fill_brewer(palette=opt$palette) +
    coord_flip() +
    ggtitle("Total genomic space of peaks (peak length sum)")
  ggsave(plot = p, filename = paste0(opt$input, '.peak_space.', opt$format), 
         width = 6, height = 4)
  
  # peak length distribution
  p <- ggplot(lengthdata, aes(x = File)) + 
    geom_boxplot(aes(ymin = Percentile5, lower = FirstQuartile, middle = Median, 
                     upper = ThirdQuartile, ymax = Percentile95),
                 stat = "identity") + 
    scale_colour_brewer(palette=opt$palette) +
    coord_flip() +
    ggtitle("Peak length distribution (5-95% range)")
  ggsave(plot = p, filename = paste0(opt$input, '.peak_lengths.', opt$format), 
         width = 6, height = 4)
}



# deduplication plot
ddupfile <- paste0(opt$input, '.dedup-stats.txt')
if(file.exists(ddupfile)) {
  ddupdata <- read.table(ddupfile, header=TRUE, sep="\t", 
                        check.names = FALSE, na.strings = '.')
  
  # XY scatter comparison
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
  ggsave(plot = p, filename = paste0(opt$input, '.dedup-comparison.', opt$format),
         width = 6, height = 4)
  
  # stacked bar chart of counts
  colnames(ddupdata)[3] <- "Optical"
  colnames(ddupdata)[5] <- "Unique"
  colnames(ddupdata)[7] <- "Retained"
  ddupdata$Discarded <- ddupdata$Before - ddupdata$Retained
  ddupdata_long <- melt(ddupdata, id="File", 
                        measure.vars = c("Optical", "Discarded", "Retained", "Unique"),
                        variable.name = "Duplication", value.name = "Count")
  p <- ggplot(ddupdata_long, aes(File, y = Count)) +
    geom_col(aes(fill=Duplication), position = "stack") +
    ggtitle("Fragment Duplication Counts") + 
    coord_flip()
  ggsave(plot = p, filename = paste0(opt$input, '.duplicate-counts.', opt$format),
         width = 6, height = 4)
}



# qvalue heat map table table
qfile <- paste0(opt$input, "_maxQvalue.txt.gz")
if(file.exists(qfile)) {
    qdata = read.table(qfile, header=TRUE, sep="\t", 
                       row.names = 1, check.names = FALSE, na.strings = '.')
    if ( colnames(qdata)[1] == "Name" ) {
      # exclude this column name
      qdata <- qdata[,2:ncol(qdata)]
    }
    if (ncol(qdata) >= 2 && ncol(qdata) <= allowed) {
      # only make a heat map file if there is more than one sample
      qdata[is.na(qdata)] <- 0
      plot_mean_k_hm(qdata, 0, opt$qmax, 4, colorRampPalette(brewer.pal(9, 'Reds'))(255), 
             color4, "Mean Q-Value over merged peaks, 4 Clusters",
             paste0(opt$input,"_qvalue_hm"))
    }
 }



# log2 Fold Enrichment table
lfefile <- paste0(opt$input,"_meanLog2FE.txt.gz")
if(file.exists(lfefile)) {
    lfedata = read.table(lfefile,header=TRUE,sep="\t", 
                         row.names = 1, check.names = FALSE, na.strings = '.')
    if ( colnames(lfedata)[1] == "Name" ) {
      # exclude this column name
      lfedata <- lfedata[,2:ncol(lfedata)]
    }
    if (ncol(lfedata) >= 2 && ncol(lfedata) <= allowed) {
      # only make heat map files if there is more than one sample
      lfedata[is.na(lfedata)] <- 0
      clrs <- colorRampPalette(rev(brewer.pal(9, 'RdBu')))(255)
      plot_mean_hm(lfedata, opt$min, opt$max, clrs, 
             "Mean Log2 Fold Enrichment over merged peaks",
             paste0(opt$input,"_log2FE_hm"))
      if ( nrow(lfedata) >= opt$mincluster) {
        plot_mean_k_hm(lfedata, opt$min, opt$max, 4, clrs, color4, 
                  "Mean Log2 Fold Enrichment over merged peaks, 4 Clusters",
                  paste0(opt$input,"_log2FE_hm_K4"))
        plot_mean_k_hm(lfedata, opt$min, opt$max, 6, clrs, color6, 
                 "Mean Log2 Fold Enrichment over merged peaks, 6 Clusters",
                 paste0(opt$input,"_log2FE_hm_K6"))
        plot_mean_k_hm(lfedata, opt$min, opt$max, 8, clrs, color8, 
                 "Mean Log2 Fold Enrichment over merged peaks, 8 Clusters",
                 paste0(opt$input,"_log2FE_hm_K8"))
      }
    }
}



# mean fragment profile heat map
mfproffile <- paste0(opt$input, "_profile_mean_fragment.txt.gz")
if(file.exists(mfproffile)) {
  fprofdata = read.table(mfproffile,header=TRUE,sep="\t", 
                         row.names = 1, check.names = FALSE, na.strings = '.')
  if ( colnames(fprofdata)[1] == "Name" ) {
    # exclude this column name
    fprofdata <- fprofdata[,2:ncol(fprofdata)]
  }
  fprofdata[is.na(fprofdata)] <- 0
  grpdata = read.table(paste0(opt$input, "_profile_mean_fragment.col_groups.txt"), header=TRUE, sep="\t", 
                       row.names = 1, check.names = FALSE)
  if (length(as.vector(unique(grpdata[,1]))) <= allowed) {
    # only plot if number of groups is tolerated by by color palette
    clrs <- colorRampPalette(brewer.pal(9, 'YlOrRd'))(255)
    plot_profile_hm(fprofdata, 0, opt$fmax, clrs, grpdata,
                    "ChIP fragment density profile around merged peak midpoints",
                    paste0(opt$input,"_profile_mean_fragment_hm"))
    if (nrow(fprofdata) >= opt$mincluster) {
      plot_profile_k_hm(fprofdata, 0, opt$fmax, 4, clrs, color4, grpdata, 
                        "ChIP fragment density profile around merged peak midpoints, 4 Clusters",
                        paste0(opt$input, "_profile_mean_fragment_hm_K4"))
    }
  }
}

# pre-sorted mean fragment profile heat map
mfsproffile <- paste0(opt$input, "_profile_mean_fragment.sorted.txt.gz")
if(file.exists(mfsproffile)) {
  fprofdata = read.table(mfsproffile,header=TRUE,sep="\t", 
                         row.names = 1, check.names = FALSE, na.strings = '.')
  if ( colnames(fprofdata)[1] == "Name" ) {
    # exclude this column name
    fprofdata <- fprofdata[,2:ncol(fprofdata)]
  }
  fprofdata[is.na(fprofdata)] <- 0
  cgrpdata = read.table(paste0(opt$input, "_profile_mean_fragment.col_groups.txt"), header=TRUE, sep="\t", 
                       row.names = 1, check.names = FALSE)
  rgrpdata = read.table(paste0(opt$input, "_profile_mean_fragment.sorted.row_groups.txt"), header=TRUE,
                        sep="\t", row.names = 1, check.names = FALSE)
  if (length(as.vector(unique(grpdata[,1]))) <= allowed) {
    # only plot if number of groups is tolerated by by color palette
    clrs <- colorRampPalette(brewer.pal(9, 'YlOrRd'))(255)
    plot_sorted_profile_hm(fprofdata, 0, opt$fmax, clrs, rgrpdata, cgrpdata,
                    "ChIP fragment density profile around merged peak midpoints",
                    paste0(opt$input,"_profile_mean_fragment_sorted_hm"))
  }
}



# replicate fragment profile heat map
rfproffile <- paste0(opt$input, "_profile_replicate_fragment.txt.gz")
if(file.exists(rfproffile)) {
  fprofdata = read.table(rfproffile,header=TRUE,sep="\t", 
                         row.names = 1, check.names = F, na.strings = '.')
  if ( colnames(fprofdata)[1] == "Name" ) {
    # exclude this column name
    fprofdata <- fprofdata[,2:ncol(fprofdata)]
  }
  fprofdata[is.na(fprofdata)] <- 0
  grpdata = read.table(paste0(opt$input, "_profile_replicate_fragment.col_groups.txt"), header=TRUE, sep="\t", 
                       row.names = 1, check.names = FALSE)
  if (length(as.vector(unique(grpdata[,1]))) <= allowed) {
    # only plot if number of groups is tolerated by color palette
    clrs <- colorRampPalette(brewer.pal(9, 'YlOrRd'))(255)
    plot_profile_hm(fprofdata, 0, opt$fmax, clrs, grpdata,
                    "ChIP fragment density profile around merged peak midpoints",
                    paste0(opt$input,"_profile_replicate_fragment_hm"))
    if (nrow(fprofdata) >= opt$mincluster) {
      plot_profile_k_hm(fprofdata, 0, opt$fmax, 4, clrs, color4, grpdata, 
                        "ChIP fragment density profile around merged peak midpoints, 4 Clusters",
                        paste0(opt$input, "_profile_replicate_fragment_hm_K4"))
    }
  }
}

# pre-sorted replicate fragment profile heat map
rfsproffile <- paste0(opt$input, "_profile_replicate_fragment.sorted.txt.gz")
if(file.exists(rfsproffile)) {
  fprofdata = read.table(rfsproffile,header=TRUE,sep="\t", 
                         row.names = 1, check.names = F, na.strings = '.')
  if ( colnames(fprofdata)[1] == "Name" ) {
    # exclude this column name
    fprofdata <- fprofdata[,2:ncol(fprofdata)]
  }
  fprofdata[is.na(fprofdata)] <- 0
  cgrpdata = read.table(paste0(opt$input, "_profile_replicate_fragment.col_groups.txt"), header=TRUE, sep="\t", 
                       row.names = 1, check.names = FALSE)
  rgrpdata = read.table(paste0(opt$input, "_profile_replicate_fragment.sorted.row_groups.txt"), header=TRUE,
                        sep="\t", row.names = 1, check.names = FALSE)
  if (length(as.vector(unique(grpdata[,1]))) <= allowed) {
    # only plot if number of groups is tolerated by color palette
    clrs <- colorRampPalette(brewer.pal(9, 'YlOrRd'))(255)
    plot_sorted_profile_hm(fprofdata, 0, opt$fmax, clrs, rgrpdata, cgrpdata,
                    "ChIP fragment density profile around merged peak midpoints",
                    paste0(opt$input,"_profile_replicate_fragment_sorted_hm"))
  }
}


# log2FE profile heat map
lfeproffile <- paste0(opt$input, "_profile_log2FE.txt.gz")
if(file.exists(lfeproffile)) {
  lfeprofdata = read.table(lfeproffile,header=TRUE,sep="\t", 
                         row.names = 1, check.names = F, na.strings = '.')
  if ( colnames(lfeprofdata)[1] == "Name" ) {
    # exclude this column name
    lfeprofdata <- lfeprofdata[,2:ncol(lfeprofdata)]
  }
  lfeprofdata[is.na(lfeprofdata)] <- 0
  grpdata = read.table(paste0(opt$input, "_profile_log2FE.col_groups.txt"), header=TRUE, sep="\t", 
                       row.names = 1, check.names = F)
  if (length(as.vector(unique(grpdata[,1]))) <= allowed) {
    # only plot if number of groups tolerated by palette
    clrs <- colorRampPalette(rev(brewer.pal(9, 'RdBu')))(255)
    plot_profile_hm(lfeprofdata, opt$min, opt$max, clrs, grpdata,
                    "ChIP Log2 Fold Enrichment profile around merged peak midpoints",
                    paste0(opt$input,"_profile_log2FE_hm"))
    if (nrow(lfeprofdata) >= opt$mincluster) {
      plot_profile_k_hm(lfeprofdata, opt$min, opt$max, 4, clrs, color4, grpdata, 
                        "ChIP Log2 Fold Enrichment profile around merged peak midpoints, 4 Clusters",
                        paste0(opt$input, "_profile_log2FE_hm_K4"))
    }
  }
}

# pre-sorted log2FE profile heat map
lfesproffile <- paste0(opt$input, "_profile_log2FE.sorted.txt.gz")
if(file.exists(lfesproffile)) {
  lfeprofdata = read.table(lfesproffile,header=TRUE,sep="\t", 
                           row.names = 1, check.names = F, na.strings = '.')
  if ( colnames(lfeprofdata)[1] == "Name" ) {
    # exclude this column name
    lfeprofdata <- lfeprofdata[,2:ncol(lfeprofdata)]
  }
  lfeprofdata[is.na(lfeprofdata)] <- 0
  cgrpdata = read.table(paste0(opt$input, "_profile_log2FE.col_groups.txt"), header=TRUE, sep="\t", 
                       row.names = 1, check.names = F)
  rgrpdata = read.table(paste0(opt$input, "_profile_log2FE.sorted.row_groups.txt"), header=TRUE,
                        sep="\t", row.names = 1, check.names = F)
  if (length(as.vector(unique(grpdata[,1]))) <= allowed) {
    # only plot if number of groups tolerated by palette
    clrs <- colorRampPalette(rev(brewer.pal(9, 'RdBu')))(255)
    plot_sorted_profile_hm(lfeprofdata, opt$min, opt$max, clrs, rgrpdata, cgrpdata,
                    "ChIP Log2 Fold Enrichment profile around merged peak midpoints",
                    paste0(opt$input,"_profile_log2FE_hm"))
  }
}



# plot fragment mean profile line plot
inpath <- dirname(opt$input)
inbase <- basename(opt$input)
fprofsumfiles <- list.files(path = inpath, full.names = TRUE,
                            pattern = paste0(inbase, '_profile_mean_fragment.*_summary\\.txt'))
for (fprofsumfile in fprofsumfiles) {
  fprofsumdata <- read.table(fprofsumfile, header=TRUE, sep="\t", 
                            row.names = 1, check.names = F, na.strings = '.')
  fprofsumdata[is.na(fprofsumdata)] <- 0
  plotMid(fprofsumdata, "Fragment Density", "Mean fragment density profile around merged peak midpoints", 
          gsub("\\.txt", "", fprofsumfile) )
}


# plot replicate fragment mean profile line plot
repfprofsumfiles <- list.files(path = inpath, full.names = TRUE,
                               pattern = paste0(inbase, '_profile_replicate_fragment.*_summary\\.txt'))
for (fprofsumfile in repfprofsumfiles) {
  fprofsumdata <- read.table(fprofsumfile, header=TRUE, sep="\t", 
                             row.names = 1, check.names = F, na.strings = '.')
  fprofsumdata[is.na(fprofsumdata)] <- 0
  plotMid(fprofsumdata, "Fragment Density", "Mean fragment density profile around merged peak midpoints", 
          gsub("\\.txt", "", fprofsumfile) )
}


# plot log2FE mean profile line plot
lfeprofsumfiles <- list.files(path = inpath, full.names = TRUE,
                              pattern = paste0(inbase, "_profile_log2FE.*_summary\\.txt"))
for (lfeprofsumfile in lfeprofsumfiles) {
  lfeprofsumdata <- read.table(lfeprofsumfile, header=TRUE, sep="\t", 
                             row.names = 1, check.names = F, na.strings = '.')
  lfeprofsumdata[is.na(lfeprofsumdata)] <- 0
  plotMid(lfeprofsumdata, "Log2 Fold Enrichment", "Mean Log2 Fold Enrichment profile around merged peak midpoints", 
          gsub("\\.txt", "", lfeprofsumfile) )
  
}


